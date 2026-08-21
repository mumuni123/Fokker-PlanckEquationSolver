#!/usr/bin/env python3
"""Validate one accepted H9 flux-interface production run.

This tool deliberately reads only files written after accepted VPFP steps.  It
does not replay an operator, inspect trial states, or change the solver.  The
output is a single key=value result file suitable for a batch-script gate.
"""

from __future__ import print_function

import argparse
import glob
import math
import os
import sys


def read_table(path):
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
        rows.append(row)
    return header, rows


def finite_values(rows, names):
    values = []
    for row in rows:
        for name in names:
            if name in row and isinstance(row[name], (int, float)):
                values.append(row[name])
    return bool(values) and all(math.isfinite(value) for value in values)


def maximum(rows, name, default=math.nan):
    values = [row[name] for row in rows
              if name in row and isinstance(row[name], (int, float)) and
              math.isfinite(row[name])]
    return max(values) if values else default


def maximum_absolute(rows, name, default=math.nan):
    values = [abs(row[name]) for row in rows
              if name in row and isinstance(row[name], (int, float)) and
              math.isfinite(row[name])]
    return max(values) if values else default


def latest(rows, name, default=math.nan):
    values = [row[name] for row in rows
              if name in row and isinstance(row[name], (int, float)) and
              math.isfinite(row[name])]
    return values[-1] if values else default


def total(rows, name, default=math.nan):
    values = [row[name] for row in rows
              if name in row and isinstance(row[name], (int, float)) and
              math.isfinite(row[name])]
    return sum(values) if values else default


def parse_manifest(path):
    values = {}
    if not os.path.isfile(path):
        return values
    with open(path, "r") as stream:
        for line in stream:
            fields = line.split(None, 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1].strip()
    return values


def has_nonempty_failure(root):
    path = os.path.join(root, "vpfp_failure.dat")
    return os.path.isfile(path) and os.path.getsize(path) > 0


def checkpoint_manifests(root):
    # Runtime checkpoints use the space-delimited manifest.txt contract.
    # Snapshot diagnostics have a separate manifest.dat and must not be used
    # for restart/configuration validation.  Accept manifest.dat only as a
    # backward-compatible fallback for an older checkpoint writer.
    current = glob.glob(os.path.join(root, "checkpoint_*", "manifest.txt"))
    if current:
        return current
    return glob.glob(os.path.join(root, "checkpoint_*", "manifest.dat"))


def snapshot_has_threshold_files(root):
    pattern = os.path.join(root, "snapshot_*", "tail_threshold_interface_rank*.dat")
    return bool(glob.glob(pattern))


def bool_text(value):
    return "1" if value else "0"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, help="accepted production output directory")
    parser.add_argument("--result", required=True)
    parser.add_argument("--mode", choices=("beam3", "beam12", "nobeam40", "beam40"),
                        required=True)
    parser.add_argument("--min-accepted-steps", type=int, default=1)
    parser.add_argument("--max-gauss-linf", type=float, default=1.0e-8)
    parser.add_argument("--max-tail-balance", type=float, default=1.0e-9)
    parser.add_argument("--max-conversion-residual", type=float, default=1.0e-9)
    parser.add_argument("--max-tail-particles", type=float, default=-1.0)
    parser.add_argument("--max-local-tail-particles", type=float, default=-1.0)
    parser.add_argument("--max-last-wall-s", type=float, default=-1.0)
    parser.add_argument("--zero-conversion-tol", type=float, default=1.0e-300)
    parser.add_argument("--tail-growth-rtol", type=float, default=1.0e-12)
    parser.add_argument("--require-tail", choices=("auto", "yes", "no"), default="auto")
    parser.add_argument("--require-threshold-snapshot", action="store_true")
    args = parser.parse_args()

    step_header, steps = read_table(os.path.join(args.run, "vpfp_step_diagnostics.dat"))
    flux_header, flux = read_table(os.path.join(args.run, "bulk_tail_flux_accepted_steps.dat"))
    accepted = [row for row in steps if row.get("accepted") == 1.0]
    rejected = [row for row in steps if row.get("accepted") != 1.0]
    final_tail = latest(accepted, "tail_particle_count")
    max_tail = maximum(accepted, "tail_particle_count")
    max_local_tail = maximum(accepted, "tail_particles_local_max")
    # gauss_linf is the dimensional residual of the finite-difference
    # Poisson operator.  Its magnitude follows the charge normalization and
    # is not comparable to a dimensionless tolerance.  The solver's
    # boundary_charge_residual is written as gauss_charge_residual and is the
    # normalized compatibility quantity used by the production Gauss gate.
    max_gauss = maximum(accepted, "gauss_linf")
    max_gauss_charge_residual = maximum_absolute(
        accepted, "gauss_charge_residual")
    max_balance = maximum(accepted, "tail_number_balance_error")
    tail_balance_abs = max((abs(row.get("tail_number_balance_error", math.nan))
                            for row in accepted
                            if isinstance(row.get("tail_number_balance_error"), (int, float)) and
                            math.isfinite(row.get("tail_number_balance_error"))),
                           default=math.nan)
    last100 = accepted[-100:]
    max_last_wall = maximum(last100, "wall_s")
    max_conversion = max(
        maximum(flux, "conversion_number_residual"),
        maximum(flux, "conversion_px_residual"),
        maximum(flux, "conversion_energy_residual"),
        maximum(flux, "conversion_jx_residual"),
        maximum(flux, "conversion_pixx_residual"),
        maximum(flux, "conversion_piperp_residual"))
    duplicate_count = total(flux, "duplicate_count")
    duplicate_id_count = total(flux, "duplicate_id_count")
    face_mismatch = total(flux, "face_ledger_mismatch_count")
    static_extractor = total(flux, "static_extractor_call_count")
    tail_source = total(flux, "particles_created")
    tail_present = math.isfinite(max_tail) and max_tail > 0.0
    zero_conversion_tail_growth = 0.0
    for row in accepted:
        removed = row.get("conversion_N")
        before = row.get("N_tail_before")
        after = row.get("N_tail_after")
        if not all(isinstance(value, (int, float)) and math.isfinite(value)
                   for value in (removed, before, after)):
            continue
        if abs(removed) <= args.zero_conversion_tol:
            zero_conversion_tail_growth = max(zero_conversion_tail_growth,
                                              max(0.0, after - before))
    zero_conversion_growth_limit = args.tail_growth_rtol * max(
        1.0, maximum(accepted, "N_tail_before", 0.0))
    if args.require_tail == "auto":
        require_tail = args.mode in ("beam12", "beam40")
    else:
        require_tail = args.require_tail == "yes"

    manifests = [parse_manifest(path) for path in checkpoint_manifests(args.run)]
    manifest_ok = bool(manifests) and all(
        item.get("tail_conversion_mode") == "flux-interface" and
        item.get("collision_induced_conversion") == "1" and
        item.get("tail_collision_weight_algorithm") == "sentoku-kemp-bounded-v1" and
        item.get("population_control_interval") == "0"
        for item in manifests)
    required_step_columns = ("accepted", "gauss_linf", "tail_particle_count",
                             "tail_number_balance_error", "wall_s")
    required_flux_columns = ("duplicate_count", "duplicate_id_count",
                             "face_ledger_mismatch_count", "static_extractor_call_count",
                             "conversion_number_residual", "conversion_px_residual",
                             "conversion_energy_residual", "conversion_jx_residual",
                             "conversion_pixx_residual", "conversion_piperp_residual")
    metrics_finite = finite_values(accepted, required_step_columns) and \
        finite_values(flux, required_flux_columns)
    no_failure = not has_nonempty_failure(args.run)
    gates = {
        "step_file": bool(step_header) and bool(steps),
        "flux_file": bool(flux_header) and bool(flux),
        "accepted_steps": len(accepted) >= args.min_accepted_steps and not rejected,
        "no_failure": no_failure,
        "finite": metrics_finite,
        "gauss": math.isfinite(max_gauss_charge_residual) and
                 max_gauss_charge_residual <= args.max_gauss_linf,
        "tail_balance": math.isfinite(tail_balance_abs) and
                        tail_balance_abs <= args.max_tail_balance,
        "conversion_residual": math.isfinite(max_conversion) and
                               max_conversion <= args.max_conversion_residual,
        "unique_interfaces": duplicate_count == 0.0 and duplicate_id_count == 0.0 and
                             face_mismatch == 0.0,
        "no_static_extractor": static_extractor == 0.0,
        "no_collision_induced_tail_growth":
            zero_conversion_tail_growth <= zero_conversion_growth_limit,
        "tail_requirement": tail_present if require_tail else True,
        "tail_global_budget": args.max_tail_particles < 0.0 or
                              (math.isfinite(max_tail) and max_tail <= args.max_tail_particles),
        "tail_local_budget": args.max_local_tail_particles < 0.0 or
                             (math.isfinite(max_local_tail) and
                              max_local_tail <= args.max_local_tail_particles),
        "wall_budget": args.max_last_wall_s < 0.0 or
                       (math.isfinite(max_last_wall) and max_last_wall <= args.max_last_wall_s),
        "manifest": manifest_ok,
        "threshold_snapshot": (not args.require_threshold_snapshot) or
                              snapshot_has_threshold_files(args.run),
    }
    passed = all(gates.values())
    with open(args.result, "w") as out:
        out.write("mode=%s\n" % args.mode)
        out.write("accepted_step_count=%d\n" % len(accepted))
        out.write("rejected_step_count=%d\n" % len(rejected))
        out.write("first_accepted_step=%s\n" % latest(accepted[::-1], "step"))
        out.write("last_accepted_step=%s\n" % latest(accepted, "step"))
        out.write("max_gauss_linf=%.17g\n" % max_gauss)
        out.write("max_gauss_charge_residual=%.17g\n" %
                  max_gauss_charge_residual)
        out.write("max_tail_number_balance_error=%.17g\n" % tail_balance_abs)
        out.write("max_conversion_residual=%.17g\n" % max_conversion)
        out.write("tail_particles_final=%.17g\n" % final_tail)
        out.write("tail_particles_max=%.17g\n" % max_tail)
        out.write("tail_particles_local_max=%.17g\n" % max_local_tail)
        out.write("tail_particles_created_total=%.17g\n" % tail_source)
        out.write("zero_conversion_tail_growth=%.17g\n" % zero_conversion_tail_growth)
        out.write("zero_conversion_tail_growth_limit=%.17g\n" %
                  zero_conversion_growth_limit)
        out.write("last100_wall_s_max=%.17g\n" % max_last_wall)
        out.write("checkpoint_manifest_count=%d\n" % len(manifests))
        out.write("threshold_snapshot_present=%s\n" %
                  bool_text(snapshot_has_threshold_files(args.run)))
        for name in sorted(gates):
            out.write("gate_%s=%s\n" % (name, bool_text(gates[name])))
        out.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
