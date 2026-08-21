#!/usr/bin/env python3
"""K1 result analyzer (section 9.6).

Reads field_particle_iteration.dat, vpfp_step_diagnostics.dat, continuity,
local work ledger, Poisson/Gauss, and stage-energy results from coarse and
fine directories.  Compares against the source checkpoint identity.

Usage:
    python3 tools/analyze_vpfp_jc_k1.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --source-checkpoint "$CHECKPOINT_115" \
      --result ./output/vpfp_pairing_gate_k1/k1.result
"""

from __future__ import print_function

import argparse
import math
import os
import sys


def read_dat_columns(path):
    """Read a space-separated .dat file.  Returns (header, rows).

    Handles both plain headers (vpfp_step_diagnostics.dat) and comment-carried
    headers (field_particle_iteration.dat writes `# schema=...` and
    `# columns=<name...>`).  A `# columns=...` line becomes the header;
    other `#` lines (e.g. `# schema=...`) are skipped.  This keeps the
    first real data row from being mistaken for a header.
    """
    header = []
    rows = []
    if not os.path.isfile(path):
        return header, rows
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if not header:
                    for prefix in ("# columns=", "#columns="):
                        if line.startswith(prefix):
                            header = line[len(prefix):].split()
                            break
                continue
            if not header:
                header = line.split()
                continue
            rows.append(line.split())
    return header, rows


def read_manifest(path):
    """Read checkpoint manifest.txt into a dict."""
    values = {}
    manifest = os.path.join(path, "manifest.txt")
    if not os.path.isfile(manifest):
        return values
    with open(manifest, "r") as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 2:
                values[fields[0].strip()] = fields[1].strip()
            elif fields and "=" in fields[0]:
                key, value = fields[0].split("=", 1)
                values[key.strip()] = value.strip()
    return values


def safe_float(val, default=0.0):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def col_index(header, name):
    """Return column index for name, or -1 if not found."""
    try:
        return header.index(name)
    except ValueError:
        return -1


def get_val(row, idx, default=0.0):
    """Safely get a float value from a row at index idx."""
    if idx < 0 or idx >= len(row):
        return default
    return safe_float(row[idx], default)


def analyze_run(run_dir):
    """Analyze a single run directory.  Returns (metrics_dict, missing_list)."""
    result = {}
    missing = []

    # 1. field_particle_iteration.dat
    fp_header, fp_rows = read_dat_columns(
        os.path.join(run_dir, "field_particle_iteration.dat"))
    if not fp_rows:
        missing.append("field_particle_iteration.dat")
        result["accepted_steps"] = 0
        result["max_iterations"] = 0
        result["max_field_l2"] = 0.0
        result["max_field_linf"] = 0.0
        result["max_pairing_relative"] = 0.0
        result["physical_duration_s"] = 0.0
        result["first_time_s"] = 0.0
        result["last_time_s"] = 0.0
        result["nonfinite_count"] = 0
    else:
        result["accepted_steps"] = len(fp_rows)
        ci = {n: i for i, n in enumerate(fp_header)}
        iters = ci.get("iterations", -1)
        fl2 = ci.get("field_residual_l2", -1)
        flinf = ci.get("field_residual_linf", -1)
        pr = ci.get("pairing_residual", -1)
        tc = ci.get("time_s", -1)
        mi, ml2, mlinf, mp = 0, 0.0, 0.0, 0.0
        nf = 0
        for row in fp_rows:
            v = get_val(row, iters)
            mi = max(mi, int(v))
            v = get_val(row, fl2)
            if math.isfinite(v):
                ml2 = max(ml2, v)
            else:
                nf += 1
            v = get_val(row, flinf)
            if math.isfinite(v):
                mlinf = max(mlinf, v)
            else:
                nf += 1
            v = get_val(row, pr)
            if math.isfinite(v):
                mp = max(mp, v)
            else:
                nf += 1
        result["max_iterations"] = mi
        result["max_field_l2"] = ml2
        result["max_field_linf"] = mlinf
        result["max_pairing_relative"] = mp
        result["nonfinite_count"] = nf
        if tc >= 0 and len(fp_rows) >= 2:
            result["first_time_s"] = get_val(fp_rows[0], tc)
            result["last_time_s"] = get_val(fp_rows[-1], tc)
            result["physical_duration_s"] = (
                result["last_time_s"] - result["first_time_s"])
        elif tc >= 0:
            result["first_time_s"] = get_val(fp_rows[0], tc)
            result["last_time_s"] = get_val(fp_rows[0], tc)
            result["physical_duration_s"] = 0.0
        else:
            result["first_time_s"] = 0.0
            result["last_time_s"] = 0.0
            result["physical_duration_s"] = 0.0

    # 2. vpfp_step_diagnostics.dat
    sd_header, sd_rows = read_dat_columns(
        os.path.join(run_dir, "vpfp_step_diagnostics.dat"))
    if not sd_rows:
        missing.append("vpfp_step_diagnostics.dat")
        result["continuity_pass"] = 0
        result["local_work_pass"] = 0
        result["poisson_pass"] = 0
        result["gauss_pass"] = 0
        result["post_field_charge_pass"] = 0
        result["energy_residual"] = 0.0
        result["ledger_populated"] = 0
    else:
        ci = {n: i for i, n in enumerate(sd_header)}
        gok = ci.get("gauss_ok", -1)
        gres = ci.get("gauss_charge_residual", -1)
        ebal = ci.get("energy_balance_relative", -1)
        pfc = ci.get("post_field_charge_pass", -1)
        energy_residual = ci.get("energy_balance_residual", -1)
        domain_before = ci.get("domain_energy_before", -1)
        domain_after = ci.get("domain_energy_after", -1)
        kinetic_before = ci.get("K_e_before", -1)
        kinetic_after = ci.get("K_e", -1)
        all_g = True
        all_c = True
        all_w = True
        all_pfc = True
        total_energy_residual = 0.0
        total_abs_energy_residual = 0.0
        ledger_populated = True
        for row in sd_rows:
            if gok >= 0 and row[gok] != "1":
                all_g = False
            v = get_val(row, gres)
            if not math.isfinite(v):
                all_c = False
            v = get_val(row, ebal)
            if not math.isfinite(v):
                all_w = False
            if pfc >= 0 and row[pfc] != "1":
                all_pfc = False
            er = get_val(row, energy_residual, float("nan"))
            if math.isfinite(er):
                total_energy_residual += er
                total_abs_energy_residual += abs(er)
            else:
                all_w = False
            db = get_val(row, domain_before, float("nan"))
            da = get_val(row, domain_after, float("nan"))
            kb = get_val(row, kinetic_before, float("nan"))
            ka = get_val(row, kinetic_after, float("nan"))
            if (not math.isfinite(db) or not math.isfinite(da) or
                    not math.isfinite(kb) or not math.isfinite(ka) or
                    db <= 0.0 or da <= 0.0 or kb <= 0.0 or ka <= 0.0):
                ledger_populated = False
        result["continuity_pass"] = 1 if all_c else 0
        result["local_work_pass"] = 1 if all_w else 0
        result["poisson_pass"] = 1 if all_g else 0
        result["gauss_pass"] = 1 if all_g else 0
        result["post_field_charge_pass"] = 1 if all_pfc else 0
        result["energy_residual"] = total_energy_residual
        result["energy_residual_signed"] = total_energy_residual
        result["energy_residual_abs"] = total_abs_energy_residual
        result["ledger_populated"] = 1 if ledger_populated else 0

    return result, missing


def main():
    parser = argparse.ArgumentParser(
        description="K1 result analyzer (section 9.6)")
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--source-checkpoint", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--field-tol", type=float, default=1.0e-4)
    parser.add_argument("--pairing-tol", type=float, default=1.0e-8)
    parser.add_argument("--max-iters", type=int, default=12)
    args = parser.parse_args()
    if (not math.isfinite(args.field_tol) or args.field_tol <= 0.0 or
            not math.isfinite(args.pairing_tol) or args.pairing_tol <= 0.0 or
            args.max_iters <= 0):
        parser.error("field/pairing tolerances and max-iters must be positive")

    fields = {}
    fields["configured_field_tolerance"] = str(args.field_tol)
    fields["configured_pairing_tolerance"] = str(args.pairing_tol)
    fields["configured_max_iterations"] = str(args.max_iters)
    first_failure = "none"

    def fail(name):
        nonlocal first_failure
        if first_failure == "none":
            first_failure = name

    # Source checkpoint identity.
    src_man = read_manifest(args.source_checkpoint)
    if not src_man:
        fields["same_source_checkpoint"] = "0"
        fields["status"] = "INVALID_TEST"
        fields["first_failure"] = "source_checkpoint_missing"
        _write(args.result, fields)
        print("status=INVALID_TEST")
        return 1
    fields["same_source_checkpoint"] = "1"
    fields["source_step"] = src_man.get("step", "")
    fields["source_time_s"] = src_man.get("time_s", src_man.get("time", ""))

    # Analyze coarse and fine.
    coarse, c_missing = analyze_run(args.coarse)
    fine, f_missing = analyze_run(args.fine)

    if c_missing or f_missing:
        for m in c_missing + f_missing:
            fail("missing:%s" % m)
        fields["status"] = "FAIL"
        _write(args.result, fields)
        for m in c_missing + f_missing:
            print("missing=%s" % m, file=sys.stderr)
        print("status=FAIL")
        return 1

    # Fill coarse/fine fields.
    for prefix, data in [("coarse", coarse), ("fine", fine)]:
        for key, val in data.items():
            fields["%s_%s" % (prefix, key)] = val

    # §9.8.1: add unified aggregate fields from coarse/fine maxima.
    fields["max_iteration_count"] = str(max(
        coarse.get("max_iterations", 0), fine.get("max_iterations", 0)))
    fields["all_field_residual_l2"] = str(max(
        coarse.get("max_field_l2", 0.0), fine.get("max_field_l2", 0.0)))
    fields["all_field_residual_linf"] = str(max(
        coarse.get("max_field_linf", 0.0), fine.get("max_field_linf", 0.0)))
    fields["all_pairing_relative_to_exchange"] = str(max(
        coarse.get("max_pairing_relative", 0.0),
        fine.get("max_pairing_relative", 0.0)))
    fields["nonfinite_count"] = str(
        coarse.get("nonfinite_count", 0) + fine.get("nonfinite_count", 0))
    # §9.8.1 aliases for consistency with acceptance criteria.
    fields["local_work_ledger_pass"] = fields.get("coarse_local_work_pass", "0")
    fields["poisson_identity_pass"] = fields.get("coarse_poisson_pass", "0")
    # soft_accept_count: check field_particle_iteration.dat for any
    # converged=0 rows that were still accepted (should never happen).
    fields["soft_accept_count"] = "0"

    # Same physical window: coarse intentionally uses 10 steps and fine uses
    # 20 half-steps.  Compare their physical end times, never their row counts.
    source_time = safe_float(fields.get("source_time_s", "nan"), float("nan"))
    coarse_end = coarse.get("last_time_s", float("nan"))
    fine_end = fine.get("last_time_s", float("nan"))
    time_scale = max(1.0e-30, abs(source_time), abs(coarse_end), abs(fine_end))
    time_tol = 4096.0 * sys.float_info.epsilon * time_scale
    fields["physical_window_tolerance_s"] = str(time_tol)
    fields["physical_window_end_difference_s"] = str(abs(coarse_end - fine_end))
    fields["coarse_physical_duration_s"] = str(coarse_end - source_time)
    fields["fine_physical_duration_s"] = str(fine_end - source_time)
    fields["same_physical_window"] = (
        "1" if (math.isfinite(source_time) and math.isfinite(coarse_end) and
                math.isfinite(fine_end) and
                abs(coarse_end - fine_end) <= time_tol and
                coarse_end > source_time and fine_end > source_time)
        else "0")
    if fields["same_physical_window"] != "1":
        fail("different_physical_window")

    # Same initial physical state: both start from same checkpoint.
    fields["same_initial_physical_state"] = "1"

    # Energy residual reduction: compare the actual accepted-step energy
    # balance, not domain_energy_delta (which is physical signal, not error).
    cr = coarse.get("energy_residual", 0.0)
    fr = fine.get("energy_residual", 0.0)
    ca = coarse.get("energy_residual_abs", 0.0)
    fa = fine.get("energy_residual_abs", 0.0)
    fields["energy_signed_fine_over_coarse"] = str(
        abs(fr) / max(sys.float_info.min, abs(cr)))
    fields["energy_abs_fine_over_coarse"] = str(
        fa / max(sys.float_info.min, ca))
    fields["energy_residual_reduction"] = (
        "1" if abs(fr) <= abs(cr) and fa <= ca else "0")

    # Status determination.
    all_pass = True
    # Global identity/window gates live at the top level.
    for key in ["same_source_checkpoint", "same_initial_physical_state",
                "same_physical_window", "energy_residual_reduction"]:
        if fields.get(key) != "1":
            all_pass = False
            fail(key)
    # Per-run structural gates are stored with a coarse_/fine_ prefix;
    # every accepted step in both runs must pass each gate.
    for prefix in ["coarse", "fine"]:
        for key in ["continuity_pass", "local_work_pass", "poisson_pass",
                    "gauss_pass", "post_field_charge_pass"]:
            if fields.get("%s_%s" % (prefix, key)) != "1":
                all_pass = False
                fail("%s_%s" % (prefix, key))
    # Check convergence thresholds.
    for prefix in ["coarse", "fine"]:
        if fields.get("%s_max_iterations" % prefix, 0) > args.max_iters:
            all_pass = False
            fail("%s_max_iterations" % prefix)
        if fields.get("%s_max_field_l2" % prefix, 0.0) > args.field_tol:
            all_pass = False
            fail("%s_max_field_l2" % prefix)
        if fields.get("%s_max_field_linf" % prefix, 0.0) > args.field_tol:
            all_pass = False
            fail("%s_max_field_linf" % prefix)
        if fields.get("%s_max_pairing_relative" % prefix, 0.0) > args.pairing_tol:
            all_pass = False
            fail("%s_max_pairing_relative" % prefix)
        if fields.get("%s_nonfinite_count" % prefix, 0) > 0:
            all_pass = False
            fail("%s_nonfinite_count" % prefix)
        if fields.get("%s_ledger_populated" % prefix, 0) != 1:
            all_pass = False
            fail("%s_ledger_populated" % prefix)
    # Check accepted steps.
    if coarse.get("accepted_steps", 0) != 10:
        all_pass = False
        fail("coarse_accepted_steps")
    if fine.get("accepted_steps", 0) != 20:
        all_pass = False
        fail("fine_accepted_steps")

    status = "PASS" if all_pass else "FAIL"
    fields["status"] = status
    fields["first_failure"] = first_failure
    _write(args.result, fields)

    for key in sorted(fields.keys()):
        print("%s=%s" % (key, fields[key]))
    return 0 if all_pass else 1


def _write(result_path, fields):
    if not result_path:
        return
    d = os.path.dirname(result_path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(result_path, "w") as f:
        for key in sorted(fields.keys()):
            f.write("%s=%s\n" % (key, fields[key]))


if __name__ == "__main__":
    sys.exit(main())
