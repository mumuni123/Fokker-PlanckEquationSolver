#!/usr/bin/env python3
"""K1 §15.14.9 step 2 analyzer: direct face dual audit.

Verifies the two required scalar identities:
    R_face == dual_in_domain_work
    R_face + boundary_force_work == current_pair_residual
where R_face is evaluated only on owner faces using the recorded quadrature
weight.  It also reports owner/duplicate/nonfinite counts and region splits.

Usage:
    python3 tools/analyze_vpfp_k1_dual_face.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --result ./output/vpfp_pairing_gate_k1/dual_face.result
"""

from __future__ import print_function

import argparse
import glob
import math
import os
import re
import sys

PAIRING_FILE = "field_particle_power_pairing.dat"
DUAL_FACE_GLOB = "field_particle_power_dual_face_*_rank*.dat"


def safe_float(val, default=0.0):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def col_index(header, name):
    try:
        return header.index(name)
    except ValueError:
        return -1


def read_dat_columns(path):
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


def analyze_run(run_dir):
    """Analyze one run directory.  Returns (metrics, errors)."""
    metrics = {}
    errors = []

    header, rows = read_dat_columns(os.path.join(run_dir, PAIRING_FILE))
    if not rows:
        errors.append("%s missing/empty" % PAIRING_FILE)
        return metrics, errors
    ci = {n: i for i, n in enumerate(header)}
    step_i = ci.get("step", -1)
    cpr_i = ci.get("current_pair_residual", -1)
    dt_i = ci.get("dt_s", -1)
    l5_i = ci.get("dual_left5_integral", -1)
    c90_i = ci.get("dual_core90_integral", -1)
    r5_i = ci.get("dual_right5_integral", -1)
    dual_in_i = ci.get("dual_in_domain_work", -1)
    boundary_i = ci.get("boundary_force_work", -1)

    # Require the corrected v3 scalar contract.  Older pairing files are
    # intentionally rejected instead of being interpreted with guessed column
    # offsets, because that would silently invalidate the dual audit.
    if cpr_i < 0 or l5_i < 0 or dt_i < 0 or dual_in_i < 0 or boundary_i < 0:
        errors.append("pairing columns missing (needs dt_s, current_pair_residual, "
                      "dual region/in-domain and boundary work columns)")
        return metrics, errors

    accepted = []
    for r in rows:
        if len(r) < len(header):
            continue
        if len(r) > step_i and r[step_i] == "#":
            continue
        accepted.append(r)

    metrics["accepted_steps"] = len(accepted)

    # Scalar reconstruction from the globalized region integrals.
    recon_by_step = {}
    region_by_step = {}
    for r in accepted:
        step = int(get_val_r(r, step_i, 0))
        dt = safe_float(get_val_r(r, dt_i, float("nan")), float("nan"))
        cpr = safe_float(get_val_r(r, cpr_i, float("nan")), float("nan"))
        l5 = safe_float(get_val_r(r, l5_i, float("nan")), float("nan"))
        c90 = safe_float(get_val_r(r, c90_i, float("nan")), float("nan"))
        r5 = safe_float(get_val_r(r, r5_i, float("nan")), float("nan"))
        dual_in = safe_float(get_val_r(r, dual_in_i, float("nan")), float("nan"))
        boundary = safe_float(get_val_r(r, boundary_i, float("nan")), float("nan"))
        recon = l5 + c90 + r5
        recon_by_step[step] = (dt, cpr, recon, dual_in, boundary)
        region_by_step[step] = (l5, c90, r5)

    # Direct face integration across the per-rank face files.  Each rank owns
    # its local faces; the global sum over all ranks is the global integral.
    # Shared-face ownership matches the scalar integral: a rank's right face
    # carries weight 0 unless it is the global physical right endpoint, so a
    # face file's LAST row is skipped unless its global_face == nx_global.
    face_recon_by_step = {}
    face_count_by_step = {}
    face_region_by_step = {}
    face_owner_count_by_step = {}
    face_duplicate_count_by_step = {}
    face_nonfinite_count_by_step = {}
    face_seen_by_step = {}
    face_owner_seen_by_step = {}
    face_owner_duplicate_count_by_step = {}
    face_owner_rank_error_by_step = {}
    # Determine nx_global as the max global_face seen across all face files.
    nx_global = -1
    face_paths = glob.glob(os.path.join(run_dir, DUAL_FACE_GLOB))
    for path in face_paths:
        fh, frows = read_dat_columns(path)
        if not frows or not fh:
            continue
        fci = {n: i for i, n in enumerate(fh)}
        g_i = fci.get("global_face", -1)
        if g_i < 0:
            continue
        for fr in frows:
            nx_global = max(nx_global,
                            int(safe_float(get_val_r(fr, g_i, 0.0), 0.0)))
    for path in face_paths:
        fh, frows = read_dat_columns(path)
        if not frows or not fh:
            continue
        m = re.search(r"_(\d+)_rank", os.path.basename(path))
        if not m:
            continue
        step = int(m.group(1))
        rank_match = re.search(r"_rank(\d+)\.dat$", os.path.basename(path))
        file_rank = int(rank_match.group(1)) if rank_match else -1
        fci = {n: i for i, n in enumerate(fh)}
        e_i = fci.get("E_pair", -1)
        d_i = fci.get("J_charge_face_minus_Gstar_J_force_face", -1)
        r_i = fci.get("region_id", -1)
        g_i = fci.get("global_face", -1)
        owner_i = fci.get("face_is_owner", -1)
        owner_rank_i = fci.get("face_owner_rank", -1)
        weight_i = fci.get("quadrature_weight", -1)
        if e_i < 0 or d_i < 0 or owner_i < 0 or owner_rank_i < 0 or weight_i < 0:
            errors.append("step %d face schema missing owner/weight columns" % step)
            continue
        s = face_recon_by_step.setdefault(step, 0.0)
        count = face_count_by_step.setdefault(step, 0)
        regs = face_region_by_step.setdefault(
            step, {0: 0.0, 1: 0.0, 2: 0.0})
        seen = face_seen_by_step.setdefault(step, {})
        owner_count = face_owner_count_by_step.setdefault(step, 0)
        duplicate_count = face_duplicate_count_by_step.setdefault(step, 0)
        nonfinite_count = face_nonfinite_count_by_step.setdefault(step, 0)
        owner_seen = face_owner_seen_by_step.setdefault(step, set())
        owner_duplicate_count = face_owner_duplicate_count_by_step.setdefault(step, 0)
        owner_rank_error = face_owner_rank_error_by_step.setdefault(step, 0)
        dt = recon_by_step.get(step, (float("nan"),))[0]
        if not math.isfinite(dt):
            errors.append("step %d missing/nonfinite dt_s" % step)
            continue
        for fr in frows:
            if len(fr) <= d_i:
                continue
            gf = int(safe_float(get_val_r(fr, g_i, 0.0), 0.0))
            seen[gf] = seen.get(gf, 0) + 1
            if seen[gf] > 1:
                duplicate_count += 1
            e_val = safe_float(get_val_r(fr, e_i, float("nan")), float("nan"))
            d_val = safe_float(get_val_r(fr, d_i, float("nan")), float("nan"))
            w_val = safe_float(get_val_r(fr, weight_i, float("nan")), float("nan"))
            owner = int(safe_float(get_val_r(fr, owner_i, -1.0), -1.0))
            owner_rank = int(safe_float(get_val_r(fr, owner_rank_i, -1.0), -1.0))
            if not (math.isfinite(e_val) and math.isfinite(d_val) and
                    math.isfinite(w_val)):
                nonfinite_count += 1
                continue
            if owner != 1:
                if owner_rank == file_rank:
                    owner_rank_error += 1
                continue
            if owner_rank != file_rank:
                owner_rank_error += 1
            owner_count += 1
            if gf in owner_seen:
                owner_duplicate_count += 1
            owner_seen.add(gf)
            reg = int(safe_float(get_val_r(fr, r_i, 1.0), 1.0))
            contribution = -dt * e_val * d_val * w_val
            s += contribution
            regs[reg] += contribution
            count += 1
        face_recon_by_step[step] = s
        face_count_by_step[step] = count
        face_owner_count_by_step[step] = owner_count
        face_duplicate_count_by_step[step] = duplicate_count
        face_nonfinite_count_by_step[step] = nonfinite_count
        face_owner_duplicate_count_by_step[step] = owner_duplicate_count
        face_owner_rank_error_by_step[step] = owner_rank_error

    # Compare scalar reconstruction and direct face reconstruction.
    max_rel_cpr = 0.0
    max_rel_face = 0.0
    max_face_abs = 0.0
    steps_checked = 0
    region_total = {0: 0.0, 1: 0.0, 2: 0.0}
    region_abs = {0: 0.0, 1: 0.0, 2: 0.0}
    max_rel_dual = 0.0
    max_rel_current = 0.0
    max_owner_count = 0
    max_duplicate_count = 0
    max_nonfinite_count = 0
    max_missing_owner_count = 0
    max_owner_duplicate_count = 0
    max_owner_rank_error = 0
    direct_region_total = {0: 0.0, 1: 0.0, 2: 0.0}
    direct_region_abs = {0: 0.0, 1: 0.0, 2: 0.0}
    for step in sorted(recon_by_step.keys()):
        dt, cpr, recon, dual_in, boundary = recon_by_step[step]
        if not all(math.isfinite(v) for v in (dt, cpr, recon, dual_in, boundary)):
            errors.append("step %d non-finite scalar recon" % step)
            continue
        scale = max(1.0, abs(cpr), abs(recon), abs(dual_in), abs(boundary))
        # The region sum is the in-domain dual, not the full current-pair
        # residual; the latter additionally contains boundary force work.
        rel_scalar = abs(recon - dual_in) / scale
        max_rel_cpr = max(max_rel_cpr, rel_scalar)
        fr = face_recon_by_step.get(step, float("nan"))
        if math.isfinite(fr):
            rel_dual = abs(fr - dual_in) / scale
            rel_current = abs(fr + boundary - cpr) / scale
            max_rel_dual = max(max_rel_dual, rel_dual)
            max_rel_current = max(max_rel_current, rel_current)
            max_rel_face = max(max_rel_face, rel_current)
            max_face_abs = max(max_face_abs, abs(fr + boundary - cpr))
            if rel_dual > 1.0e-10:
                errors.append("step %d direct_face_vs_in_domain" % step)
            if rel_current > 1.0e-10:
                errors.append("step %d direct_face_plus_boundary_vs_current" % step)
        else:
            errors.append("step %d missing face reconstruction" % step)
        expected_faces = nx_global + 1 if nx_global >= 0 else 0
        missing = max(0, expected_faces -
                      len(face_owner_seen_by_step.get(step, set())))
        max_owner_count = max(max_owner_count,
                              face_owner_count_by_step.get(step, 0))
        max_duplicate_count = max(max_duplicate_count,
                                  face_duplicate_count_by_step.get(step, 0))
        max_nonfinite_count = max(max_nonfinite_count,
                                  face_nonfinite_count_by_step.get(step, 0))
        max_missing_owner_count = max(max_missing_owner_count, missing)
        max_owner_duplicate_count = max(
            max_owner_duplicate_count,
            face_owner_duplicate_count_by_step.get(step, 0))
        max_owner_rank_error = max(
            max_owner_rank_error, face_owner_rank_error_by_step.get(step, 0))
        if missing > 0:
            errors.append("step %d missing owner faces=%d" % (step, missing))
        if face_owner_duplicate_count_by_step.get(step, 0) > 0:
            errors.append("step %d duplicate owner faces=%d" % (
                step, face_owner_duplicate_count_by_step[step]))
        regs = region_by_step.get(step, (0.0, 0.0, 0.0))
        for reg, val in enumerate(regs):
            region_total[reg] += val
            region_abs[reg] += abs(val)
        direct_regs = face_region_by_step.get(step, {0: 0.0, 1: 0.0, 2: 0.0})
        for reg in range(3):
            direct_region_total[reg] += direct_regs.get(reg, 0.0)
            direct_region_abs[reg] += abs(direct_regs.get(reg, 0.0))
        steps_checked += 1

    metrics["steps_checked"] = steps_checked
    metrics["max_scalar_recon_relative"] = max_rel_cpr
    metrics["max_face_recon_relative"] = max_rel_face
    metrics["max_dual_in_domain_relative"] = max_rel_dual
    metrics["max_current_pair_reconstruction_relative"] = max_rel_current
    metrics["max_face_recon_abs"] = max_face_abs
    metrics["recon_left5_integral"] = region_total[0]
    metrics["recon_core90_integral"] = region_total[1]
    metrics["recon_right5_integral"] = region_total[2]
    metrics["recon_left5_abs"] = region_abs[0]
    metrics["recon_core90_abs"] = region_abs[1]
    metrics["recon_right5_abs"] = region_abs[2]
    metrics["direct_face_left5_integral"] = direct_region_total[0]
    metrics["direct_face_core90_integral"] = direct_region_total[1]
    metrics["direct_face_right5_integral"] = direct_region_total[2]
    metrics["direct_face_left5_abs"] = direct_region_abs[0]
    metrics["direct_face_core90_abs"] = direct_region_abs[1]
    metrics["direct_face_right5_abs"] = direct_region_abs[2]
    metrics["owner_face_count_max"] = max_owner_count
    metrics["duplicate_face_count_max"] = max_duplicate_count
    metrics["nonfinite_face_count_max"] = max_nonfinite_count
    metrics["missing_owner_face_count_max"] = max_missing_owner_count
    metrics["owner_duplicate_face_count_max"] = max_owner_duplicate_count
    metrics["owner_rank_error_count_max"] = max_owner_rank_error
    for step in sorted(recon_by_step.keys()):
        metrics["step_%d_owner_face_count" % step] = \
            face_owner_count_by_step.get(step, 0)
        metrics["step_%d_duplicate_face_count" % step] = \
            face_duplicate_count_by_step.get(step, 0)
        metrics["step_%d_owner_duplicate_count" % step] = \
            face_owner_duplicate_count_by_step.get(step, 0)
        metrics["step_%d_owner_rank_error" % step] = \
            face_owner_rank_error_by_step.get(step, 0)
        metrics["step_%d_nonfinite_face_count" % step] = \
            face_nonfinite_count_by_step.get(step, 0)
        metrics["step_%d_direct_left5" % step] = \
            face_region_by_step.get(step, {}).get(0, 0.0)
        metrics["step_%d_direct_core90" % step] = \
            face_region_by_step.get(step, {}).get(1, 0.0)
        metrics["step_%d_direct_right5" % step] = \
            face_region_by_step.get(step, {}).get(2, 0.0)
    metrics["dual_recon_pass"] = (
        1 if (steps_checked == metrics["accepted_steps"] and
              steps_checked > 0 and
              max_rel_cpr <= 1.0e-10 and
              max_rel_dual <= 1.0e-10 and
              max_rel_current <= 1.0e-10 and
              max_owner_duplicate_count == 0 and
              max_owner_rank_error == 0 and
              max_nonfinite_count == 0 and
              max_missing_owner_count == 0) else 0)

    # Left5/core90/right5 relative dominance (signed, like the pairing file).
    core90_scale = max(1.0, abs(region_total[1]),
                       max(abs(region_total[0]), abs(region_total[2])))
    metrics["core90_dominance_ratio"] = (
        abs(region_total[1]) / core90_scale)
    metrics["left5_dominance_ratio"] = (
        abs(region_total[0]) / core90_scale)
    metrics["right5_dominance_ratio"] = (
        abs(region_total[2]) / core90_scale)
    return metrics, errors


def get_val_r(row, idx, default=0.0):
    if idx < 0 or idx >= len(row):
        return default
    return safe_float(row[idx], default)


def write_result(path, fields):
    with open(path, "w") as f:
        for key in sorted(fields.keys()):
            f.write("%s=%s\n" % (key, fields[key]))


def main():
    parser = argparse.ArgumentParser(
        description="K1 dual-face direct audit analyzer (section 15.13.4)")
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    fields = {}
    first_failure = "none"

    def fail(name):
        nonlocal first_failure
        if first_failure == "none":
            first_failure = name

    for prefix, run_dir in [("coarse", args.coarse), ("fine", args.fine)]:
        metrics, errors = analyze_run(run_dir)
        for key, val in metrics.items():
            fields["%s_%s" % (prefix, key)] = val
        if errors:
            for err in errors:
                fail("%s_%s" % (prefix, err.replace(" ", "_")))
        if metrics.get("accepted_steps", 0) <= 0:
            fail("%s_accepted_steps" % prefix)
        if metrics.get("dual_recon_pass", 0) != 1:
            fail("%s_dual_recon_pass" % prefix)

    fields["first_failure"] = first_failure
    fields["status"] = "PASS" if first_failure == "none" else "FAIL"
    write_result(args.result, fields)
    for key in sorted(fields.keys()):
        print("%s=%s" % (key, fields[key]))
    return 0 if first_failure == "none" else 1


if __name__ == "__main__":
    sys.exit(main())
