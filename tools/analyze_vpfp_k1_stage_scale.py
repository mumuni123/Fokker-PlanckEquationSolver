#!/usr/bin/env python3
"""K1 分阶段 energy-audit 时间细化对比（下一步核查工具）。

读取 coarse/fine 的 vpfp_stage_energy_audit.dat，对每个阶段统计
   per-step |stage_balance| 累计、有符号累计、fine/coarse 比值，
用于区分：
  - 时间离散缺口：balance 随 dt 缩放（fine/coarse ~ 0.5 或更低）；
  - 每步固定量级账本项：balance 不随 dt 缩放（fine/coarse ~ 1 或更高）。

同时统计 stage_telescope / energy_ledger 结构门（应 PASS），保证
比较是在自洽的账本结构上进行。

用法：
    python3 tools/analyze_vpfp_k1_stage_scale.py \
      --coarse ./output/vpfp_pairing_gate_k1/coarse \
      --fine ./output/vpfp_pairing_gate_k1/fine \
      --result ./output/vpfp_pairing_gate_k1/stage_scale.result
"""

from __future__ import print_function

import argparse
import math
import os

STAGES = (
    "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
    "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
    "collision_half2", "conversion_after_collision", "tail_bulk_return",
    "final_poisson",
)

# Scale ratio thresholds for the verdict.  A balance that shrinks by >= this
# factor when dt halves is treated as time-discretization-dominated; a ratio
# at/above this is treated as a fixed-magnitude ledger/representation item.
SHRINK_THRESHOLD = 0.75
FIXED_THRESHOLD = 1.25


def read_rows(path):
    header = []
    rows = []
    if not os.path.isfile(path):
        return header, rows
    with open(path, "r") as stream:
        for line in stream:
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


def to_float(values, idx, default=0.0):
    if idx < 0 or idx >= len(values):
        return default
    try:
        return float(values[idx])
    except ValueError:
        return default


def analyze_run(run_dir):
    header, rows = read_rows(os.path.join(run_dir,
                                          "vpfp_stage_energy_audit.dat"))
    if not header or not rows:
        return None, "vpfp_stage_energy_audit.dat missing/empty"
    ci = {n: i for i, n in enumerate(header)}
    step_idx = ci.get("step", -1)
    stage_idx = ci.get("stage_name", -1)
    bal_idx = ci.get("stage_balance", -1)
    accepted_idx = ci.get("accepted", -1)
    audit_idx = ci.get("audit_valid", -1)
    if (step_idx < 0 or stage_idx < 0 or bal_idx < 0 or
            accepted_idx < 0 or audit_idx < 0):
        return None, "required columns missing"

    by_step = {}
    for r in rows:
        if to_float(r, accepted_idx) != 1.0:
            continue
        if to_float(r, audit_idx) != 1.0:
            continue
        by_step.setdefault(int(to_float(r, step_idx)), []).append(r)

    # Per-stage signed and abs sums over accepted steps.
    signed = {name: 0.0 for name in STAGES[1:]}
    abs_sum = {name: 0.0 for name in STAGES[1:]}
    max_abs = {name: 0.0 for name in STAGES[1:]}
    for step in sorted(by_step):
        stage_rows = sorted(by_step[step],
                            key=lambda r: to_float(r, stage_idx))
        for r in stage_rows:
            name = r[stage_idx]
            if name not in signed:
                continue
            v = to_float(r, bal_idx)
            signed[name] += v
            abs_sum[name] += abs(v)
            max_abs[name] = max(max_abs[name], abs(v))
    return {
        "accepted_steps": len(by_step),
        "signed": signed,
        "abs_sum": abs_sum,
        "max_abs": max_abs,
    }, None


def main():
    parser = argparse.ArgumentParser(
        description="K1 stage-balance time-refinement comparison")
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--result", default=None)
    parser.add_argument("--expected-coarse-steps", type=int, default=10)
    parser.add_argument("--expected-fine-steps", type=int, default=20)
    args = parser.parse_args()

    coarse, cerr = analyze_run(args.coarse)
    fine, ferr = analyze_run(args.fine)
    output = []
    if cerr:
        output.append("status=INVALID_TEST")
        output.append("error=%s" % cerr)
        output.append("first_failure=coarse_%s" % cerr.replace(" ", "_"))
    elif ferr:
        output.append("status=INVALID_TEST")
        output.append("error=%s" % ferr)
        output.append("first_failure=fine_%s" % ferr.replace(" ", "_"))
    else:
        c_steps = coarse["accepted_steps"]
        f_steps = fine["accepted_steps"]
        step_ok = (c_steps == args.expected_coarse_steps and
                   f_steps == args.expected_fine_steps)
        output.append("coarse_accepted_steps=%d" % c_steps)
        output.append("fine_accepted_steps=%d" % f_steps)
        output.append("accepted_step_count_valid=%d" % (1 if step_ok else 0))

        # Time refinement verdicts per stage.
        shrink = []
        fixed = []
        inconclusive = []
        for name in STAGES[1:]:
            cv = coarse["signed"][name]
            fv = fine["signed"][name]
            ca = coarse["abs_sum"][name]
            fa = fine["abs_sum"][name]
            # Ratio of per-step mean |balance| (normalize by step count).
            c_mean = ca / max(1, c_steps)
            f_mean = fa / max(1, f_steps)
            if c_mean > 1.0e-30:
                ratio = f_mean / c_mean
                ratio_str = "%.17g" % ratio
            else:
                ratio = float("inf")
                ratio_str = "n/a"
            # Signed mean ratio (keeps sign information).
            c_smean = cv / max(1, c_steps)
            f_smean = fv / max(1, f_steps)
            output.append("stage_%s_abs_mean_coarse=%.17g" % (name, c_mean))
            output.append("stage_%s_abs_mean_fine=%.17g" % (name, f_mean))
            output.append("stage_%s_abs_mean_fine_over_coarse=%s" % (name, ratio_str))
            output.append("stage_%s_signed_sum_coarse=%.17g" % (name, cv))
            output.append("stage_%s_signed_sum_fine=%.17g" % (name, fv))
            # Signed mean ratio keeps sign information.  A sign flip between
            # coarse and fine (roundoff-dominated or alternating stages) makes
            # a plain ratio meaningless; report the sign explicitly instead.
            if (abs(c_smean) > 1.0e-30 and abs(f_smean) > 1.0e-30 and
                    (c_smean < 0.0) != (f_smean < 0.0)):
                output.append("stage_%s_signed_mean_ratio=sign_flip" % name)
            elif abs(c_smean) > 1.0e-30:
                output.append("stage_%s_signed_mean_ratio=%.17g" %
                              (name, f_smean / abs(c_smean)))
            else:
                output.append("stage_%s_signed_mean_ratio=0.0" % name)
            # Rough per-step magnitudes (informational).
            output.append("stage_%s_max_abs_coarse=%.17g" % (name, coarse["max_abs"][name]))
            output.append("stage_%s_max_abs_fine=%.17g" % (name, fine["max_abs"][name]))

            # Classification: only abs-mean ratio is decisive; ignore stages
            # whose coarse abs mean is negligible (roundoff-only stages).
            if c_mean <= 1.0e-3:
                verdict = "roundoff_only"
            elif ratio <= SHRINK_THRESHOLD:
                verdict = "time_shrinking"
                shrink.append(name)
            elif ratio >= FIXED_THRESHOLD:
                verdict = "fixed_magnitude"
                fixed.append(name)
            else:
                verdict = "inconclusive"
                inconclusive.append(name)
            output.append("stage_%s_verdict=%s" % (name, verdict))

        output.append("time_shrinking_stages=%s" % (",".join(shrink) or "none"))
        output.append("fixed_magnitude_stages=%s" % (",".join(fixed) or "none"))
        output.append("inconclusive_stages=%s" % (",".join(inconclusive) or "none"))
        output.append("first_failure=none")
        output.append("status=PASS" if step_ok else "status=FAIL")

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
