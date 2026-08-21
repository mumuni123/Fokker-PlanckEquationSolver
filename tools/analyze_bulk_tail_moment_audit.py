#!/usr/bin/env python3
# Section 7.11.16B / 16A: analyse the real-moment audit main table and the
# velocity-cell histogram.  Cumulative ratios are computed from the RAW
# center_l1_X / delta_l1_X columns (never the arithmetic mean of per-step
# ratios).  Outputs audit_summary.result with the real event/cell counts,
# cumulative six-moment ratios, maximum cell ratio, below-threshold number
# fraction, feasibility counts and decision=GREEN|GRAY|RED|INVALID.
#
# Usage:
#   python tools/analyze_bulk_tail_moment_audit.py \
#     --audit <accepted_steps.dat> \
#     --velocity-histogram <velocity_histogram.dat> \
#     --result <audit_summary.result>

import argparse
import math
import sys

COMPONENTS = ["N", "Px", "Jx", "K", "Pixx", "Piperp"]


def rel(a, b):
    """r(a,b): a/b if b>0; 0 if a==0,b==0; +inf if a>0,b==0."""
    if b > 0.0:
        return a / b, 1
    if a == 0.0:
        return 0.0, 1
    return float("inf"), 0


def _is_number(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def parse_table(path, errors):
    header = None
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("INVALID_EVENT_LAYOUT"):
                errors.append("event layout inconsistency in " + path)
                continue
            parts = line.split()
            if header is None:
                header = parts
                continue
            rows.append([float(x) if _is_number(x) else x for x in parts])
    return header, rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--audit", required=True)
    p.add_argument("--velocity-histogram", required=True)
    p.add_argument("--result", required=True)
    args = p.parse_args()

    errors = []
    header, rows = parse_table(args.audit, errors)
    required = ["accepted_step", "time_fs", "conversion_location",
                "positive_request_cell_count",
                "below_threshold_number_fraction",
                "volume_target_feasible_count",
                "volume_target_failed_count",
                "eligible_target_feasible_count",
                "eligible_target_failed_count"]
    for m in range(6):
        required += ["center_%d" % m, "volume_%d" % m,
                     "eligible_raw_%d" % m, "eligible_normalized_%d" % m,
                     "delta_signed_%d" % m, "center_l1_%d" % m,
                     "delta_l1_%d" % m, "rel_l1_%d" % m,
                     "rel_signed_%d" % m, "max_cell_rel_%d" % m,
                     "relative_defined_%d" % m]
    if header is None:
        errors.append("missing audit header")
    else:
        for col in required:
            if col not in header:
                errors.append("missing column " + col)

    sum_center_l1 = [0.0] * 6
    sum_delta_l1 = [0.0] * 6
    max_cell_rel = [0.0] * 6
    undefined_relative = [0] * 6
    sum_center_n = 0.0
    sum_eligible_raw_n = 0.0
    volume_feasible = 0
    volume_failed = 0
    eligible_feasible = 0
    eligible_failed = 0
    seen_events = set()
    event_count = 0
    total_cells = 0
    finite_ok = True
    if header is not None:
        col = {name: i for i, name in enumerate(header)}
        for row in rows:
            if len(row) != len(header):
                errors.append("malformed audit row (column count)")
                continue
            step = row[col["accepted_step"]]
            loc = row[col["conversion_location"]]
            key = (step, loc)
            if key in seen_events:
                errors.append("duplicate event (step %s, location %s)" %
                              (step, loc))
            seen_events.add(key)
            event_count += 1
            total_cells += int(row[col["positive_request_cell_count"]])
            sum_center_n += float(row[col["center_0"]])
            sum_eligible_raw_n += float(row[col["eligible_raw_0"]])
            volume_feasible += int(row[col["volume_target_feasible_count"]])
            volume_failed += int(row[col["volume_target_failed_count"]])
            eligible_feasible += int(row[col["eligible_target_feasible_count"]])
            eligible_failed += int(row[col["eligible_target_failed_count"]])
            for m in range(6):
                c = float(row[col["center_l1_%d" % m]])
                d = float(row[col["delta_l1_%d" % m]])
                mr = float(row[col["max_cell_rel_%d" % m]])
                if not (math.isfinite(c) and math.isfinite(d)):
                    finite_ok = False
                sum_center_l1[m] += c
                sum_delta_l1[m] += d
                max_cell_rel[m] = max(max_cell_rel[m], mr)
                undefined_relative[m] += int(
                    row[col["relative_defined_%d" % m]] == 0)
    if not finite_ok:
        errors.append("non-finite raw center_l1/delta_l1")

    hist_cells = 0
    seen_bins = set()
    with open(args.velocity_histogram, "r") as f:
        first = True
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if first:
                first = False
                if parts != ["accepted_step", "time_fs",
                             "conversion_location", "iv", "imu",
                             "request_cell_count", "request_number"]:
                    errors.append("invalid velocity histogram header")
                continue
            try:
                step = parts[0]
                loc = parts[2]
                iv = parts[3]
                imu = parts[4]
                number = float(parts[6])
            except (IndexError, ValueError):
                errors.append("invalid velocity histogram row")
                continue
            if number <= 0.0:
                errors.append("zero-mass velocity histogram record")
            key = (step, loc, iv, imu)
            if key in seen_bins:
                errors.append("duplicate velocity bin (%s,%s,%s,%s)" % key)
            seen_bins.add(key)
            hist_cells += 1

    ratios = []
    for m in range(6):
        r, defined = rel(sum_delta_l1[m], sum_center_l1[m])
        ratios.append((r, defined))
    below_threshold_fraction = (
        1.0 - sum_eligible_raw_n / sum_center_n
        if sum_center_n > 0.0 else 0.0)

    invalid = bool(errors)
    if invalid:
        decision = "INVALID"
    else:
        max_ratio = max(r for r, _ in ratios)
        if max_ratio <= 1.0e-6:
            decision = "GREEN"
        elif max_ratio >= 1.0e-2:
            decision = "RED"
        else:
            decision = "GRAY"

    with open(args.result, "w") as out:
        out.write("input_valid=%d\n" % (0 if invalid else 1))
        out.write("real_events=%d\n" % event_count)
        # real_cells: the sum of positive-mass conversion cells over the real
        # accepted events (the section 7.11.16A "743 cells" quantity);
        # velocity_bin_rows: the number of (event, iv, imu) histogram rows.
        out.write("real_cells=%d\n" % total_cells)
        out.write("velocity_bin_rows=%d\n" % hist_cells)
        out.write("total_request_cells=%d\n" % total_cells)
        out.write("below_threshold_number_fraction=%.17g\n" %
                  below_threshold_fraction)
        for m in range(6):
            out.write("cumulative_R_L1_%s=%.17g\n" %
                      (COMPONENTS[m], ratios[m][0]))
            out.write("max_cell_rel_%s=%.17g\n" %
                      (COMPONENTS[m], max_cell_rel[m]))
            out.write("undefined_relative_count_%s=%d\n" %
                      (COMPONENTS[m], undefined_relative[m]))
        out.write("volume_target_feasible_count=%d\n" % volume_feasible)
        out.write("volume_target_failed_count=%d\n" % volume_failed)
        out.write("eligible_target_feasible_count=%d\n" % eligible_feasible)
        out.write("eligible_target_failed_count=%d\n" % eligible_failed)
        out.write("decision=%s\n" % decision)
        if errors:
            for e in errors:
                out.write("error=%s\n" % e)

    print("input_valid=%d real_events=%d real_cells=%d decision=%s" %
          (0 if invalid else 1, event_count, total_cells, decision))
    return 0 if not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
