#!/usr/bin/env python3
"""Compare overlapping accepted rows from direct and restart H9 production runs.

This is a physical-state regression, not a byte-hash test: OpenMP reductions
may differ by ULP while a valid restart must preserve accepted-step sequence
and all global macroscopic ledgers within the requested tolerance.
"""

from __future__ import print_function

import argparse
import hashlib
import math
import os
import sys


DEFAULT_FIELDS = (
    "time_s", "gauss_charge_residual", "N_e_after", "N_b_after", "U_E", "K_e", "K_b",
    "P_bkg", "P_tail", "P_combined", "K_tail", "K_combined",
    "N_tail_after", "N_combined_after", "tail_number_balance_error",
    "tail_particle_count", "collision_reservoir", "conversion_N", "conversion_Px",
    "conversion_K", "collision_flux_export_N", "collision_flux_export_K",
)

# ``gauss_linf`` is the dimensional, unnormalised residual of the discrete
# Poisson operator.  Its raw magnitude depends on cancellation order and is
# not a restart-equivalence norm.  Keep it in the result as a diagnostic, but
# gate restart consistency on the dimensionless charge compatibility instead.
DIAGNOSTIC_FIELDS = ("gauss_linf",)
ABSOLUTE_ERROR_FIELDS = ("gauss_charge_residual",)


def read_table(root):
    path = os.path.join(root, "vpfp_step_diagnostics.dat")
    if not os.path.isfile(path):
        return [], []
    with open(path, "r") as stream:
        raw = [line.split() for line in stream
               if line.strip() and not line.lstrip().startswith("#")]
    if len(raw) < 2:
        return raw[0] if raw else [], []
    header = raw[0]
    rows = []
    for fields in raw[1:]:
        if len(fields) != len(header):
            continue
        row = {}
        for key, value in zip(header, fields):
            try:
                row[key] = float(value)
            except ValueError:
                row[key] = value
        if row.get("accepted") == 1.0 and "step" in row:
            rows.append(row)
    return header, rows


def nonempty_failure(root):
    path = os.path.join(root, "vpfp_failure.dat")
    return os.path.isfile(path) and os.path.getsize(path) > 0


def checkpoint_hashes(root):
    if not os.path.isdir(root):
        raise RuntimeError("checkpoint directory does not exist: %s" % root)
    result = {}
    for name in sorted(os.listdir(root)):
        if not name.startswith("rank_") or not name.endswith(".bin"):
            continue
        path = os.path.join(root, name)
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        result[name] = digest.hexdigest()
    if not result:
        raise RuntimeError("no rank_*.bin files found in %s" % root)
    return result


def relative_error(left, right, atol):
    return abs(left - right) / max(atol, abs(left), abs(right))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct", required=True)
    parser.add_argument("--restart", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--atol", type=float, default=1.0e-300)
    parser.add_argument("--gauss-charge-atol", type=float, default=1.0e-12)
    parser.add_argument("--min-common-steps", type=int, default=1)
    parser.add_argument("--require-tail-nonempty", action="store_true")
    parser.add_argument("--direct-checkpoint")
    parser.add_argument("--restart-checkpoint")
    parser.add_argument("--require-checkpoint-hash", action="store_true")
    parser.add_argument("--checkpoint-hash-informational", action="store_true")
    args = parser.parse_args()

    direct_header, direct_rows = read_table(args.direct)
    restart_header, restart_rows = read_table(args.restart)
    direct = {int(row["step"]): row for row in direct_rows}
    restart = {int(row["step"]): row for row in restart_rows}
    common_steps = sorted(set(direct) & set(restart))
    if args.require_checkpoint_hash and args.checkpoint_hash_informational:
        parser.error("checkpoint hash cannot be both required and informational")
    common_fields = [name for name in DEFAULT_FIELDS
                     if name in direct_header and name in restart_header]
    diagnostic_fields = [name for name in DIAGNOSTIC_FIELDS
                         if name in direct_header and name in restart_header]
    errors = {}
    finite = True
    for name in common_fields:
        worst = 0.0
        for step in common_steps:
            left = direct[step].get(name)
            right = restart[step].get(name)
            if not isinstance(left, (int, float)) or not isinstance(right, (int, float)) or \
                    not math.isfinite(left) or not math.isfinite(right):
                finite = False
                worst = math.inf
                break
            if name in ABSOLUTE_ERROR_FIELDS:
                worst = max(worst, abs(left - right))
            else:
                worst = max(worst, relative_error(left, right, args.atol))
        errors[name] = worst
    diagnostic_errors = {}
    for name in diagnostic_fields:
        worst = 0.0
        for step in common_steps:
            left = direct[step].get(name)
            right = restart[step].get(name)
            if not isinstance(left, (int, float)) or not isinstance(right, (int, float)) or \
                    not math.isfinite(left) or not math.isfinite(right):
                finite = False
                worst = math.inf
                break
            worst = max(worst, relative_error(left, right, args.atol))
        diagnostic_errors[name] = worst
    gated_errors = [value for name, value in errors.items()
                    if name not in ABSOLUTE_ERROR_FIELDS]
    max_error = max(gated_errors) if gated_errors else 0.0
    gauss_charge_error = errors.get("gauss_charge_residual", 0.0)
    direct_tail = max((row.get("tail_particle_count", 0.0) for row in direct_rows), default=0.0)
    restart_tail = max((row.get("tail_particle_count", 0.0) for row in restart_rows), default=0.0)
    checkpoint_hash_equal = None
    checkpoint_hash_error = "none"
    hash_checked = args.require_checkpoint_hash or args.checkpoint_hash_informational
    if hash_checked:
        try:
            if not args.direct_checkpoint or not args.restart_checkpoint:
                raise RuntimeError("--require-checkpoint-hash requires both checkpoint paths")
            checkpoint_hash_equal = (checkpoint_hashes(args.direct_checkpoint) ==
                                     checkpoint_hashes(args.restart_checkpoint))
        except RuntimeError as error:
            checkpoint_hash_equal = False
            checkpoint_hash_error = str(error)
    gates = {
        "direct_rows": bool(direct_rows),
        "restart_rows": bool(restart_rows),
        "common_steps": len(common_steps) >= args.min_common_steps,
        "same_last_step": bool(direct_rows) and bool(restart_rows) and
                          int(direct_rows[-1]["step"]) == int(restart_rows[-1]["step"]),
        "no_direct_failure": not nonempty_failure(args.direct),
        "no_restart_failure": not nonempty_failure(args.restart),
        "fields_present": bool(common_fields),
        "finite": finite,
        "tolerance": math.isfinite(max_error) and max_error <= args.rtol,
        "gauss_charge": math.isfinite(gauss_charge_error) and
                        gauss_charge_error <= args.gauss_charge_atol,
        "tail_nonempty": (not args.require_tail_nonempty) or
                         (direct_tail > 0.0 and restart_tail > 0.0),
        "checkpoint_hash": (not args.require_checkpoint_hash) or
                           checkpoint_hash_equal,
    }
    passed = all(gates.values())
    with open(args.result, "w") as out:
        out.write("direct_accepted_steps=%d\n" % len(direct_rows))
        out.write("restart_accepted_steps=%d\n" % len(restart_rows))
        out.write("common_step_count=%d\n" % len(common_steps))
        out.write("common_first_step=%d\n" % (common_steps[0] if common_steps else -1))
        out.write("common_last_step=%d\n" % (common_steps[-1] if common_steps else -1))
        out.write("direct_tail_particles_max=%.17g\n" % direct_tail)
        out.write("restart_tail_particles_max=%.17g\n" % restart_tail)
        out.write("max_relative_error=%.17g\n" % max_error)
        out.write("absolute_error_gauss_charge_residual=%.17g\n" % gauss_charge_error)
        out.write("checkpoint_hash_checked=%d\n" % (1 if hash_checked else 0))
        out.write("checkpoint_hash_is_informational=%d\n" %
                  (1 if args.checkpoint_hash_informational else 0))
        out.write("checkpoint_hash_equal=%s\n" %
                  ("not_checked" if checkpoint_hash_equal is None else
                   ("1" if checkpoint_hash_equal else "0")))
        out.write("checkpoint_hash_error=%s\n" % checkpoint_hash_error)
        for name in common_fields:
            label = "absolute_error" if name in ABSOLUTE_ERROR_FIELDS else "relative_error"
            out.write("%s_%s=%.17g\n" % (label, name, errors[name]))
        for name in diagnostic_fields:
            out.write("diagnostic_relative_error_%s=%.17g\n" %
                      (name, diagnostic_errors[name]))
        for name in sorted(gates):
            out.write("gate_%s=%d\n" % (name, 1 if gates[name] else 0))
        out.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
