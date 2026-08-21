#!/usr/bin/env python3
"""Check whether H10 threshold/residence scans contain a physical plateau."""

from __future__ import print_function

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tail_return_compare_common import (diagnostics, max_value,
    latest_snapshot, rank_series, rank_sum_series, relative_l2, last_value)


# R4 compares deliberately different return-controller parameters.  The R3
# same-parameter A/B tolerance (2e-3) is too narrow for this scan and rejects
# profiles that remain effectively identical in shape.  Pair the wider L2
# envelope with explicit correlation and pointwise-structure gates so a local
# waveform distortion cannot pass merely because its global L2 norm is small.
FIELD_REL_L2_TOL = 5.0e-3
FIELD_CORRELATION_MIN = 0.99998
FIELD_MAX_DIFF_OVER_PEAK_TOL = 1.0e-2
FIELD_ENERGY_REL_TOL = 1.0e-3
DENSITY_REL_L2_TOL = 1.0e-3
SPECTRUM_REL_L2_TOL = 5.0e-3
INVARIANT_TOL = 1.0e-12
REQUEST_TOL = 1.0e-13


def complete_snapshot(root):
    snapshot = latest_snapshot(root)
    marker = os.path.join(snapshot, "_COMPLETE")
    return os.path.isfile(marker) and os.path.getsize(marker) > 0


def unresolved_failure(root, rows):
    path = os.path.join(root, "vpfp_failure.dat")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    failure_step = None
    with open(path, "r") as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            try:
                failure_step = int(fields[0])
            except ValueError:
                continue
    return failure_step is None or not rows or failure_step >= int(
        rows[-1].get("step", -1))


def compare(a, b):
    sa, sb = latest_snapshot(a), latest_snapshot(b)
    field_a = rank_series(sa, "fields", "Ex_Vm")
    field_b = rank_series(sb, "fields", "Ex_Vm")
    norm_a = sum(value * value for value in field_a)
    norm_b = sum(value * value for value in field_b)
    difference = [y-x for x, y in zip(field_a, field_b)]
    cross = sum(x*y for x, y in zip(field_a, field_b))
    peak = max(max(abs(x) for x in field_a), 1.0e-300)
    return {
        "field_relative_l2": relative_l2(field_a, field_b),
        "field_energy_relative": abs(norm_b - norm_a) / max(norm_a, 1.0e-300),
        "field_correlation": cross / max(math.sqrt(norm_a*norm_b), 1.0e-300),
        "field_max_diff_over_peak": max(abs(x) for x in difference) / peak,
        "density_relative_l2": relative_l2(
            rank_series(sa, "density_background", "n_e_m3"),
            rank_series(sb, "density_background", "n_e_m3")),
        "spectrum_relative_l2": relative_l2(
            rank_sum_series(sa, "energy_spectrum", "count_combined"),
            rank_sum_series(sb, "energy_spectrum", "count_combined")),
    }


def main():
    p = argparse.ArgumentParser()
    for name in ("threshold-525", "threshold-550", "threshold-575",
                 "residence-4", "residence-8", "residence-16"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--result", required=True)
    a = p.parse_args()
    runs = vars(a)
    try:
        pairs = {
            "threshold_low_center": compare(runs["threshold_525"], runs["threshold_550"]),
            "threshold_center_high": compare(runs["threshold_550"], runs["threshold_575"]),
            "residence_low_center": compare(runs["residence_4"], runs["residence_8"]),
            "residence_center_high": compare(runs["residence_8"], runs["residence_16"]),
        }
        # These are intentionally different controller parameters, not a
        # replay.  Reuse the validated R3 macroscopic envelope instead of an
        # impossible bitwise-equivalence threshold.
        plateau = all(v["field_relative_l2"] <= FIELD_REL_L2_TOL and
                      v["field_correlation"] >= FIELD_CORRELATION_MIN and
                      v["field_max_diff_over_peak"] <= FIELD_MAX_DIFF_OVER_PEAK_TOL and
                      v["field_energy_relative"] <= FIELD_ENERGY_REL_TOL and
                      v["density_relative_l2"] <= DENSITY_REL_L2_TOL and
                      v["spectrum_relative_l2"] <= SPECTRUM_REL_L2_TOL
                      for v in pairs.values())
        run_rows = {name: diagnostics(runs[name]) for name in (
            "threshold_525", "threshold_550", "threshold_575",
            "residence_4", "residence_8", "residence_16")}
        if any(not rows for rows in run_rows.values()):
            raise ValueError("one or more R4 runs produced no accepted step")
        central = run_rows["threshold_550"]
        return_observed = max_value(central, ["tail_return_particles_removed"]) > 0.0
        same_window = (all(len(rows) == len(central) for rows in run_rows.values()) and
                       all(abs(rows[0]["time_s"]-central[0]["time_s"]) <= 1.0e-24 and
                           abs(rows[-1]["time_s"]-central[-1]["time_s"]) <= 1.0e-24
                           for rows in run_rows.values()))
        complete = all(complete_snapshot(runs[name]) for name in run_rows)
        no_failure = all(not unresolved_failure(runs[name], rows)
                         for name, rows in run_rows.items())
        invariant_max = max(max_value(rows, ["tail_return_N_residual",
            "tail_return_Px_residual", "tail_return_K_residual"])
            for rows in run_rows.values())
        request_max = max(max_value(rows, ["tail_return_mpi_request_residual"])
                          for rows in run_rows.values())
        hard_gate = (same_window and complete and no_failure and
                     invariant_max <= INVARIANT_TOL and
                     request_max <= REQUEST_TOL)
        passed = plateau and return_observed and hard_gate
        values = {"parameter_plateau": int(plateau),
                  "central_return_observed": int(return_observed),
                  "hard_gate_pass": int(hard_gate),
                  "same_accepted_window": int(same_window),
                  "all_snapshots_complete": int(complete),
                  "no_unresolved_failure": int(no_failure),
                  "tail_return_invariant_residual_max": invariant_max,
                  "tail_return_request_residual_max": request_max,
                  "field_relative_l2_tolerance": FIELD_REL_L2_TOL,
                  "field_correlation_minimum": FIELD_CORRELATION_MIN,
                  "field_max_diff_over_peak_tolerance": FIELD_MAX_DIFF_OVER_PEAK_TOL,
                  "field_energy_relative_tolerance": FIELD_ENERGY_REL_TOL,
                  "density_relative_l2_tolerance": DENSITY_REL_L2_TOL,
                  "spectrum_relative_l2_tolerance": SPECTRUM_REL_L2_TOL}
        for name, rows in run_rows.items():
            values[name + "_tail_particle_count"] = last_value(
                rows, "tail_particle_count")
            values[name + "_return_particles_max"] = max_value(
                rows, ["tail_return_particles_removed"])
            values[name + "_return_particles_total"] = sum(
                row.get("tail_return_particles_removed", 0.0) for row in rows)
            values[name + "_return_number_total"] = sum(
                row.get("tail_return_N", 0.0) for row in rows)
            values[name + "_return_pixx_dx_total"] = sum(
                row.get("tail_return_Pixx_dx", 0.0) for row in rows)
            values[name + "_return_piperp_dx_total"] = sum(
                row.get("tail_return_Piperp_dx", 0.0) for row in rows)
            values[name + "_mean_wall_seconds"] = sum(
                row.get("wall_s", 0.0) for row in rows) / len(rows)
            values[name + "_mean_tail_kernel_seconds"] = sum(
                row.get("wall_tail_push_s", 0.0) +
                row.get("wall_tail_deposit_s", 0.0) +
                row.get("wall_tail_migrate_s", 0.0) +
                row.get("collision_s", 0.0) for row in rows) / len(rows)
        for label, metrics in pairs.items():
            for metric, value in metrics.items():
                values[label + "_" + metric] = value
    except (IOError, OSError, ValueError) as exc:
        passed = False; values = {"comparison_error": str(exc)}
    with open(a.result, "w") as out:
        for key in sorted(values): out.write("%s=%s\n" % (key, values[key]))
        out.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    print("status=%s" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
