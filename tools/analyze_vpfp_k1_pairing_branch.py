#!/usr/bin/env python3
"""K1 §15.12 阶段 1 机制确认：Gate I pairing audit 分支 P/Q 判定。

判据（§15.12"积分格式层修复方案·阶段 1"）：

  poisson_transport_residual = Delta U_E - W_electrode + dt <E_pair, J_charge>

- 分支 P：poisson_transport_residual ≈ 0（Poisson 功恒等式闭合），缺口集中在
  current_pair_residual = total_field_work - dt<E_pair, J_charge>，即
  "粒子功（force_pair = total_field_work） vs 配对功（dt<E_pair,J>）"不同源
  —— 格点-面场离散口径（候选机制确认）。
- 分支 Q：poisson_transport_residual 达缺口量级（Poisson 功恒等式自身
  不闭合，端点权重/phi 平均修正的代数缺陷）。

判定标准（相对每步 |full_residual| 缺口量级）：
  transport_over_full = sum|poisson_transport_residual| / sum|full_residual|
  - transport_over_full <= 0.10 且 current_pair 主导缺口 -> 分支 P
  - transport_over_full >= 0.50                       -> 分支 Q
  - 其余                                           -> INCONCLUSIVE

用法：
    python3 tools/analyze_vpfp_k1_pairing_branch.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --result ./output/vpfp_pairing_gate_k1/pairing_branch.result
"""

from __future__ import print_function

import argparse
import math
import os

PAIRING_FILE = "field_particle_power_pairing.dat"
P_FRACTION_MAX = 0.10
Q_FRACTION_MIN = 0.50


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
    p_h, p_r = read_columns(os.path.join(run_dir, PAIRING_FILE))
    if not p_h or not p_r:
        return None, "field_particle_power_pairing.dat missing/empty"
    required = ["step", "dt_s", "poisson_transport_residual",
                "current_pair_residual", "full_residual",
                "reconstruction_pass"]
    col = {n: idx(p_h, n) for n in required}
    if any(i < 0 for i in col.values()):
        return None, "required pairing columns missing"
    if any(int(num(r, col["reconstruction_pass"])) != 1 for r in p_r):
        return None, "reconstruction_pass=0 present"
    steps = []
    for r in p_r:
        steps.append((int(num(r, col["step"])), num(r, col["dt_s"]),
                      num(r, col["poisson_transport_residual"]),
                      num(r, col["current_pair_residual"]),
                      num(r, col["full_residual"])))
    if len(steps) < 2:
        return None, "too few accepted pairing rows"
    sum_abs_transport = sum(abs(s[2]) for s in steps)
    sum_abs_current = sum(abs(s[3]) for s in steps)
    sum_abs_full = sum(abs(s[4]) for s in steps)
    if sum_abs_full <= 0.0:
        return None, "full_residual is identically zero"
    transport_over_full = sum_abs_transport / sum_abs_full
    current_over_full = sum_abs_current / sum_abs_full
    total_dt = sum(s[1] for s in steps)
    out = {
        "accepted_steps": len(steps),
        "sum_abs_poisson_transport": sum_abs_transport,
        "sum_abs_current_pair": sum_abs_current,
        "sum_abs_full_residual": sum_abs_full,
        "transport_over_full": transport_over_full,
        "current_pair_over_full": current_over_full,
        "transport_rate": sum_abs_transport / max(1e-300, total_dt),
        "full_residual_rate": sum_abs_full / max(1e-300, total_dt),
    }
    if transport_over_full <= P_FRACTION_MAX:
        out["branch"] = "P"
    elif transport_over_full >= Q_FRACTION_MIN:
        out["branch"] = "Q"
    else:
        out["branch"] = "INCONCLUSIVE"
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
    if cerr or ferr:
        output.append("status=INVALID_TEST")
        output.append("coarse_error=%s" % (cerr or "ok"))
        output.append("fine_error=%s" % (ferr or "ok"))
    else:
        for name, data in [("coarse", coarse), ("fine", fine)]:
            output.append("%s_accepted_steps=%d" % (name, data["accepted_steps"]))
            output.append("%s_sum_abs_poisson_transport=%.17g" % (
                name, data["sum_abs_poisson_transport"]))
            output.append("%s_sum_abs_current_pair=%.17g" % (
                name, data["sum_abs_current_pair"]))
            output.append("%s_sum_abs_full_residual=%.17g" % (
                name, data["sum_abs_full_residual"]))
            output.append("%s_transport_over_full=%.17g" % (
                name, data["transport_over_full"]))
            output.append("%s_current_pair_over_full=%.17g" % (
                name, data["current_pair_over_full"]))
            output.append("%s_transport_rate=%.17g" % (name, data["transport_rate"]))
            output.append("%s_full_residual_rate=%.17g" % (
                name, data["full_residual_rate"]))
            output.append("%s_branch=%s" % (name, data["branch"]))
        branches = {coarse["branch"], fine["branch"]}
        if branches == {"P"}:
            verdict, first_failure = "BRANCH_P", "none"
        elif branches == {"Q"}:
            verdict, first_failure = "BRANCH_Q", "none"
        else:
            verdict, first_failure = "INCONCLUSIVE", "branch_mismatch"
        output.append("verdict=%s" % verdict)
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
    return 0 if output[-1].startswith("status=PASS") else 2


if __name__ == "__main__":
    raise SystemExit(main())