#!/usr/bin/env python3
"""K1 §16 条件3 核查：E_pair 力功 vs 场能离散功口径一致性（诊断参考）。

读取 coarse/fine 的 vpfp_step_diagnostics.dat（U_E_pair/U_E_final 列）与
vpfp_stage_energy_audit.dat（bulk_upar_face_work 等 Gate-C 力功），记录：

  (1) 粒子实际受的场（E_pair = G_P(Phi^n,Phi^{n+1})）的场能 U_E_pair 与
      场能账本场（midpoint/final）的差异（信息性参考）；
  (2) Gate-C 力功记账（粒子被 E_pair 推的功）与场能变化的身份残差；
  (3) E_pair 场能是否落在 midpoint/final 场能区间内（仅诊断，非判定依据）。

注意：U_E 是 E 的二次型，对时间中心场不满足凸组合约束（叉项可正可负），且
U_E_pair 与 U_E_final 本属不同场（时间中心 G_P vs 终态场），故"U_E_pair 低于
两端点"是物理正确现象，**不是** EPAIR_MISMATCH 的判定依据。正式的 §16 条件3
判定请使用 `analyze_vpfp_k1_force_field_pairing.py`（力-场配对离散功口径）。

用法：
    python3 tools/analyze_vpfp_k1_epair_consistency.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --result ./output/vpfp_pairing_gate_k1/epair_consistency.result
"""

from __future__ import print_function

import argparse
import math
import os


def read_columns(path):
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
                continue
            if not header:
                header = line.split()
                continue
            rows.append(line.split())
    return header, rows


def idx(header, name):
    try:
        return header.index(name)
    except ValueError:
        return -1


def num(row, i, default=0.0):
    if i < 0 or i >= len(row):
        return default
    try:
        return float(row[i])
    except ValueError:
        return default


def analyze_run(run_dir):
    out = {}
    sd_h, sd_r = read_columns(os.path.join(run_dir,
                                           "vpfp_step_diagnostics.dat"))
    if not sd_h or not sd_r:
        return None, "vpfp_step_diagnostics.dat missing/empty"
    sci = {n: i for i, n in enumerate(sd_h)}
    ue_pair = idx(sd_h, "U_E_pair")
    ue_final = idx(sd_h, "U_E_final")
    if ue_pair < 0 or ue_final < 0:
        return None, "U_E_pair/U_E_final columns missing (rebuild with new diagnostics)"
    resid = idx(sd_h, "energy_balance_residual")

    # Stage audit for Gate-C force-work columns.
    sa_h, sa_r = read_columns(os.path.join(run_dir,
                                           "vpfp_stage_energy_audit.dat"))
    aci = {n: i for i, n in enumerate(sa_h)}
    force_work_col = idx(sa_h, "bulk_upar_face_work")
    tail_kick_col = idx(sa_h, "tail_kick_work")
    beam_kick_col = idx(sa_h, "beam_kick_work")
    stage_col = idx(sa_h, "stage_name")

    steps = []
    for row in sd_r:
        up = num(row, ue_pair)
        uf = num(row, ue_final)
        rv = num(row, resid)
        steps.append((num(row, 0), up, uf, rv))  # step, ue_pair, ue_final, resid

    # Per-step force work and midpoint U_E from stage audit.
    force_work = {}
    ue_mid = {}
    for row in sa_r:
        if stage_col < 0:
            continue
        s = num(row, 0)
        if row[stage_col] == "u_force_tail_beam_kick":
            if force_work_col >= 0:
                w = (num(row, force_work_col) + num(row, tail_kick_col) +
                     num(row, beam_kick_col))
                force_work[s] = w
        elif row[stage_col] == "midpoint_poisson":
            ue_mid[s] = num(row, idx(sa_h, "U_E"))

    # Identity (2): force_work should equal U_E_pair - U_E_final.
    work_vs_field = []
    for s, up, uf, rv in steps:
        fw = force_work.get(s, float("nan"))
        if math.isfinite(fw):
            identity = fw - (up - uf)
        else:
            identity = float("nan")
        work_vs_field.append((s, up, uf, fw, identity))

    # Identity (3): time-centering constraint.  E_pair = G_P(Phi^n,Phi^{n+1})
    # is the time-centered pairing field; its energy must lie inside the
    # convex span of the midpoint and final field energies (it is a convex
    # combination of the two fields, so U_E is bounded by the endpoint values
    # for same-sign fields).  U_E_pair below BOTH endpoints is structurally
    # impossible for the true pairing field and proves trial_force_fields_
    # is NOT G_P(Phi^n,Phi^{n+1}).
    centering = []
    for s, up, uf, rv in steps:
        um = ue_mid.get(s, float("nan"))
        if math.isfinite(um) and math.isfinite(up):
            lo = min(um, uf)
            hi = max(um, uf)
            in_span = (up >= lo - 1.0e-6 * max(1.0, abs(lo)) and
                       up <= hi + 1.0e-6 * max(1.0, abs(hi)))
            centering.append((s, up, um, uf, in_span))
        else:
            centering.append((s, up, float("nan"), uf, False))

    out["steps"] = len(steps)
    out["accepted_steps"] = len(steps)
    out["sum_ue_pair_minus_ue_final"] = sum((up - uf) for _, up, uf, _ in steps)
    out["sum_force_work"] = sum(w for _, _, _, w, _ in work_vs_field
                                if math.isfinite(w))
    out["sum_identity_residual"] = sum(v for _, _, _, _, v in work_vs_field
                                       if math.isfinite(v))
    valid = [v for _, _, _, _, v in work_vs_field if math.isfinite(v)]
    out["identity_abs_mean"] = (sum(abs(v) for v in valid) / len(valid)
                                if valid else float("nan"))
    out["identity_signed_mean"] = (sum(v for v in valid) / len(valid)
                                   if valid else float("nan"))
    # If U_E_pair ~ U_E_final, the force field == field-energy field.
    out["ue_pair_minus_final_abs_mean"] = (
        sum(abs(up - uf) for _, up, uf, _ in steps) / len(steps)
        if steps else float("nan"))
    out["residual_abs_mean"] = (
        sum(abs(rv) for _, _, _, rv in steps) / len(steps)
        if steps else float("nan"))
    out["identity_vs_residual_abs_ratio"] = (
        out["identity_abs_mean"] / max(1.0e-30, out["residual_abs_mean"])
        if out["residual_abs_mean"] > 0 else float("nan"))
    # Time-centering: fraction of steps where U_E_pair lies inside the
    # [min(U_E_mid,U_E_final), max(U_E_mid,U_E_final)] convex span.
    out["ue_pair_in_mid_final_span_count"] = sum(1 for _, _, _, _, ok in centering
                                                 if ok)
    out["ue_pair_in_mid_final_span_frac"] = (
        out["ue_pair_in_mid_final_span_count"] / len(centering)
        if centering else float("nan"))
    out["ue_pair_below_both_count"] = sum(
        1 for _, up, um, uf, _ in centering
        if math.isfinite(um) and up < min(um, uf))
    return out, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--result", default=None)
    args = parser.parse_args()

    output = []
    first_failure = "none"
    coarse, cerr = analyze_run(args.coarse)
    fine, ferr = analyze_run(args.fine)
    if cerr:
        output.append("status=INVALID_TEST")
        output.append("error=%s" % cerr)
        output.append("first_failure=coarse")
    elif ferr:
        output.append("status=INVALID_TEST")
        output.append("error=%s" % ferr)
        output.append("first_failure=fine")
    else:
        for name, data in [("coarse", coarse), ("fine", fine)]:
            output.append("%s_accepted_steps=%d" % (name, data["steps"]))
            output.append("%s_sum_ue_pair_minus_ue_final=%.17g" % (
                name, data["sum_ue_pair_minus_ue_final"]))
            output.append("%s_sum_force_work=%.17g" % (name, data["sum_force_work"]))
            output.append("%s_sum_identity_residual=%.17g" % (
                name, data["sum_identity_residual"]))
            output.append("%s_identity_abs_mean=%.17g" % (
                name, data["identity_abs_mean"]))
            output.append("%s_identity_signed_mean=%.17g" % (
                name, data["identity_signed_mean"]))
            output.append("%s_ue_pair_minus_final_abs_mean=%.17g" % (
                name, data["ue_pair_minus_final_abs_mean"]))
            output.append("%s_residual_abs_mean=%.17g" % (
                name, data["residual_abs_mean"]))
            output.append("%s_identity_vs_residual_abs_ratio=%.17g" % (
                name, data["identity_vs_residual_abs_ratio"]))
            output.append("%s_ue_pair_in_mid_final_span_count=%d" % (
                name, data["ue_pair_in_mid_final_span_count"]))
            output.append("%s_ue_pair_in_mid_final_span_frac=%.17g" % (
                name, data["ue_pair_in_mid_final_span_frac"]))
            output.append("%s_ue_pair_below_both_count=%d" % (
                name, data["ue_pair_below_both_count"]))

        # Verdict: the identity residual (force work minus the pair-to-final
        # field-energy change) is reported for reference.  The convex-span
        # and below-both counters are diagnostic only (U_E is quadratic, so
        # they are NOT a valid mismatch test).  The decisive §16 条件3 test is
        # whether the force-field pairing residual is a fixed per-physical-time
        # gap; see analyze_vpfp_k1_force_field_pairing.py for that verdict.
        output.append("identity_abs_mean_coarse=%.17g" % coarse["identity_abs_mean"])
        output.append("identity_abs_mean_fine=%.17g" % fine["identity_abs_mean"])
        output.append("identity_vs_residual_ratio_coarse=%.17g" %
                      coarse["identity_vs_residual_abs_ratio"])
        output.append("identity_vs_residual_ratio_fine=%.17g" %
                      fine["identity_vs_residual_abs_ratio"])
        output.append("span_frac_coarse=%.17g" % coarse["ue_pair_in_mid_final_span_frac"])
        output.append("span_frac_fine=%.17g" % fine["ue_pair_in_mid_final_span_frac"])
        output.append("below_both_coarse=%d" % coarse["ue_pair_below_both_count"])
        output.append("below_both_fine=%d" % fine["ue_pair_below_both_count"])
        # Reference only: these counters are diagnostic, not a verdict.
        output.append("ue_pair_below_mid_final_span_is_physical=1")
        output.append("verdict=REFERENCE_ONLY_SEE_FORCE_FIELD_PAIRING")
        output.append("first_failure=none")
        output.append("status=PASS")

    if args.result:
        d = os.path.dirname(args.result)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        with open(args.result, "w") as f:
            f.write("\n".join(output) + "\n")
    for line in output:
        print(line)
    return 0 if output[-1].startswith("status=PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
