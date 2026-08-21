#!/usr/bin/env python3
"""K1 §16 条件3 核查：力功与场能变化的离散功口径一致性（修正版）。

先前的 epair_consistency 用"U_E 凸组合约束"判定，但 U_E 是 E 的二次型，
不满足凸组合约束（叉项可正可负），且 U_E_pair 与 U_E_final 本属不同场
（时间中心 G_P vs 终态场），故该判定产生误报。

本工具聚焦 §16 条件3 真正可核查的恒等式：

  residual_step = force_pair + poisson_pair + other

其中（来自 stage audit，每步）：
  force_pair   = u_force 阶段 dK_bulk+dK_tail+dK_beam + conversion dK
                 （粒子被 G_P 场推的动能变化）
  poisson_pair = midpoint_poisson dU_E + final_poisson dU_E
                 （场能变化）
  other        = x remap / collision / H10 各 stage_balance 之和

物理守恒要求 force_pair = -poisson_pair（粒子功来自场能释放），
residual = force_pair + poisson_pair 应≈0。非零 residual 即能量缺口。

判定标准（与 stage_scale 一致的单位物理时间尺度）：
  - 若 |residual| 的每单位物理时间速率 coarse/fine 接近 1，则缺口是
    "每单位物理时间恒定的口径差异"，指向力功场 vs 场能场不同源；
  - 若速率随 dt 下降，则缺口是时间离散截断误差。

用法：
    python3 tools/analyze_vpfp_k1_force_field_pairing.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --result ./output/vpfp_pairing_gate_k1/force_field_pairing.result
"""

from __future__ import print_function

import argparse
import math
import os

STAGES_FORCE = ("u_force_tail_beam_kick", "conversion_after_force")
STAGES_POISSON = ("midpoint_poisson", "final_poisson")
STAGES_OTHER = ("x_half1", "x_half2", "collision_half1", "collision_half2",
                "conversion_after_collision", "tail_bulk_return")


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
                if not header and line.startswith("# columns="):
                    header = line[len("# columns="):].split()
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
    sa_h, sa_r = read_columns(os.path.join(run_dir,
                                           "vpfp_stage_energy_audit.dat"))
    if not sa_h or not sa_r:
        return None, "vpfp_stage_energy_audit.dat missing/empty"
    sci = {n: i for i, n in enumerate(sa_h)}
    step_idx = idx(sa_h, "step")
    stage_idx = idx(sa_h, "stage_name")
    dkb = idx(sa_h, "dK_bulk")
    dkt = idx(sa_h, "dK_tail")
    dkb_beam = idx(sa_h, "dK_beam")
    due = idx(sa_h, "dU_E")
    sb = idx(sa_h, "stage_balance")
    resid_col = idx(sa_h, "energy_balance_residual")
    time_idx = idx(sa_h, "time_s")
    if any(i < 0 for i in [step_idx, stage_idx, dkb, dkt, dkb_beam, due, sb,
                           resid_col, time_idx]):
        return None, "required stage-audit columns missing"

    by_step = {}
    for row in sa_r:
        by_step.setdefault(num(row, step_idx), []).append(row)

    steps = []
    for step in sorted(by_step):
        rows = by_step[step]
        srows = {r[stage_idx]: r for r in rows}
        force_pair = 0.0
        for s in STAGES_FORCE:
            if s in srows:
                r = srows[s]
                force_pair += (num(r, dkb) + num(r, dkt) + num(r, dkb_beam))
        poisson_pair = 0.0
        for s in STAGES_POISSON:
            if s in srows:
                poisson_pair += num(srows[s], due)
        other = 0.0
        for s in STAGES_OTHER:
            if s in srows:
                other += num(srows[s], sb)
        resid = num(rows[0], resid_col)
        t = num(rows[0], time_idx)
        steps.append((step, t, force_pair, poisson_pair, other,
                      force_pair + poisson_pair + other, resid))

    out = {}
    out["accepted_steps"] = len(steps)
    out["sum_force_pair"] = sum(s[2] for s in steps)
    out["sum_poisson_pair"] = sum(s[3] for s in steps)
    out["sum_other"] = sum(s[4] for s in steps)
    out["sum_reconstructed"] = sum(s[5] for s in steps)
    out["sum_residual"] = sum(s[6] for s in steps)
    out["abs_mean_reconstructed"] = sum(abs(s[5]) for s in steps) / len(steps)
    out["abs_mean_residual"] = sum(abs(s[6]) for s in steps) / len(steps)
    out["reconstruction_vs_residual_abs_ratio"] = (
        out["abs_mean_reconstructed"] / max(1.0e-30, out["abs_mean_residual"])
        if out["abs_mean_residual"] > 0 else float("nan"))
    # Per-unit-physical-time rates.
    if len(steps) >= 2:
        t0 = steps[0][1]
        t1 = steps[-1][1]
        dt = t1 - t0
        out["window_s"] = dt
        out["rate_signed"] = out["sum_reconstructed"] / max(1e-300, dt)
        out["rate_abs"] = sum(abs(s[5]) for s in steps) / max(1e-300, dt)
        out["rate_residual_signed"] = out["sum_residual"] / max(1e-300, dt)
    else:
        out["window_s"] = 0.0
        out["rate_signed"] = float("nan")
        out["rate_abs"] = float("nan")
        out["rate_residual_signed"] = float("nan")
    # Dominant component of reconstructed residual.
    comps = {"force_pair": out["sum_force_pair"],
             "poisson_pair": out["sum_poisson_pair"],
             "other": out["sum_other"]}
    out["dominant_component"] = max(comps, key=lambda k: abs(comps[k]))
        # §16 条件3 时间中心缺口量化（方案 D，纯诊断）：residual 与力功的相关性及比例。
        # 若 residual ≈ c * force_pair（c 固定，corr→±1），缺口是力-场离散功的
        # 固定比例不闭合；corr 接近 0 则缺口与当场力功无关。
        # residual_vs_force_pair_ratio = 最小二乘斜率 c（residual = c * force_pair），
        # 即文档 §15.12 中的 -2.7~2.8% 固定比例缺口。
    fp_list = [s[2] for s in steps]
    res_list = [s[6] for s in steps]
    n = len(steps)
    if n >= 3:
        mean_fp = sum(fp_list) / n
        mean_res = sum(res_list) / n
        num_ = sum((fp_list[i] - mean_fp) * (res_list[i] - mean_res)
                   for i in range(n))
        den_fp = math.sqrt(sum((fp_list[i] - mean_fp) ** 2 for i in range(n)))
        den_res = math.sqrt(sum((res_list[i] - mean_res) ** 2 for i in range(n)))
        corr = num_ / max(1e-300, den_fp * den_res)
        out["corr_force_pair_residual"] = corr
        # residual/|force_pair| mean ratio (signed relative to force_pair).
        ratio_sum = 0.0
        ratio_n = 0
        for fp, res in zip(fp_list, res_list):
            if abs(fp) > 1.0e-30:
                ratio_sum += res / fp
                ratio_n += 1
        out["residual_over_force_pair_mean"] = (ratio_sum / ratio_n
                                                if ratio_n else float("nan"))
        # Slope of residual vs force_pair (least squares through origin):
        # residual ≈ residual_vs_force_pair_ratio * force_pair.
        slope = sum(fp_list[i] * res_list[i] for i in range(n)) / max(
            1e-300, sum(fp_list[i] ** 2 for i in range(n)))
        out["residual_vs_force_pair_ratio"] = slope
    else:
        out["corr_force_pair_residual"] = float("nan")
        out["residual_over_force_pair_mean"] = float("nan")
        out["residual_vs_force_pair_ratio"] = float("nan")
    return out, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--result", default=None)
    args = parser.parse_args()

    output = []
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
            output.append("%s_accepted_steps=%d" % (name, data["accepted_steps"]))
            output.append("%s_sum_force_pair=%.17g" % (name, data["sum_force_pair"]))
            output.append("%s_sum_poisson_pair=%.17g" % (name, data["sum_poisson_pair"]))
            output.append("%s_sum_other=%.17g" % (name, data["sum_other"]))
            output.append("%s_sum_reconstructed=%.17g" % (name, data["sum_reconstructed"]))
            output.append("%s_sum_residual=%.17g" % (name, data["sum_residual"]))
            output.append("%s_reconstruction_vs_residual_abs_ratio=%.17g" % (
                name, data["reconstruction_vs_residual_abs_ratio"]))
            output.append("%s_window_s=%.17g" % (name, data["window_s"]))
            output.append("%s_rate_signed=%.17g" % (name, data["rate_signed"]))
            output.append("%s_rate_abs=%.17g" % (name, data["rate_abs"]))
            output.append("%s_rate_residual_signed=%.17g" % (
                name, data["rate_residual_signed"]))
            output.append("%s_dominant_component=%s" % (
                name, data["dominant_component"]))
            output.append("%s_corr_force_pair_residual=%.17g" % (
                name, data["corr_force_pair_residual"]))
            output.append("%s_residual_over_force_pair_mean=%.17g" % (
                name, data["residual_over_force_pair_mean"]))
            output.append("%s_residual_vs_force_pair_ratio=%.17g" % (
                name, data["residual_vs_force_pair_ratio"]))

        # Verdict: reconstruction must exactly reproduce the ledger residual
        # (internal self-consistency of the stage decomposition).  If it does
        # and the signed rate is unchanged by dt, the gap is a fixed-magnitude
        # force-field pairing mismatch, not time truncation.
        c_ratio = coarse["reconstruction_vs_residual_abs_ratio"]
        f_ratio = fine["reconstruction_vs_residual_abs_ratio"]
        c_rate = coarse["rate_signed"]
        f_rate = fine["rate_signed"]
        reconstruct_ok = (0.99 < c_ratio < 1.01 and 0.99 < f_ratio < 1.01)
        if not reconstruct_ok:
            output.append("decomposition_self_consistent=0")
            output.append("verdict=DECOMPOSITION_BROKEN")
            first_failure = "decomposition_broken"
        else:
            output.append("decomposition_self_consistent=1")
            # Rate ratio near 1 => fixed per-physical-time gap.
            rate_ratio = f_rate / max(1.0e-30, c_rate) if abs(c_rate) > 1.0e-30 else float("nan")
            output.append("signed_rate_fine_over_coarse=%.17g" % rate_ratio)
            if 0.5 <= rate_ratio <= 2.0:
                output.append("verdict=FIXED_FORCE_FIELD_GAP")
                first_failure = "force_field_gap"
            else:
                output.append("verdict=TIME_TRUNCATION")
                first_failure = "none"
        output.append("first_failure=%s" % first_failure)
        output.append("status=PASS" if first_failure == "none" else "status=FAIL")

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
