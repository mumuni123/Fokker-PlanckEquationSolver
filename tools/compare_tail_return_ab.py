#!/usr/bin/env python3
"""Compare H10 disabled/enabled runs using physical and performance gates."""

from __future__ import print_function

import argparse
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tail_return_compare_common import (diagnostics, last_value,
    latest_snapshot, max_value, rank_series, rank_sum_series, relative_l2)


# R3a compares two evolved representations, not a bitwise replay.  N/Px/K
# remain strict invariants; these bounds limit the macroscopic perturbation
# over the short 0.3 fs window.
DENSITY_REL_L2_TOL = 1.0e-3
FIELD_REL_L2_TOL = 2.0e-3
FIELD_ENERGY_REL_TOL = 1.0e-3
SPECTRUM_REL_L2_TOL = 5.0e-3
WAVE_ENVELOPE_REL_L2_TOL = 1.0e-3
ENVELOPE_BLOCK_CELLS = 100


def parse_window(text):
    if not text:
        return None
    values = [float(x) for x in text.split(",")]
    if len(values) != 2 or values[1] <= values[0]:
        raise argparse.ArgumentTypeError("window must be start,end in fs")
    return values


def mean_kernel(rows, window):
    names = ("wall_tail_push_s", "wall_tail_deposit_s",
             "wall_tail_migrate_s", "collision_s")
    selected = []
    for row in rows:
        time_fs = row.get("time_s", 0.0) * 1.0e15
        if window and not (window[0] <= time_fs <= window[1]):
            continue
        selected.append(sum(row.get(name, 0.0) for name in names))
    return sum(selected) / len(selected) if selected else math.nan


def block_rms(values, block_size):
    return [math.sqrt(sum(x*x for x in values[i:i+block_size]) /
                      len(values[i:i+block_size]))
            for i in range(0, len(values), block_size)]


def squared_norm_relative_change(baseline, candidate):
    base = sum(x*x for x in baseline)
    cand = sum(x*x for x in candidate)
    return abs(cand-base) / max(base, 1.0e-300)


def unresolved_failure(root, accepted_rows):
    """Return true only when a failure was not superseded by later steps."""
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
    if failure_step is None or not accepted_rows:
        return True
    return failure_step >= int(accepted_rows[-1].get("step", -1))


def require_complete_snapshot(snapshot):
    marker = os.path.join(snapshot, "_COMPLETE")
    if not os.path.isfile(marker) or os.path.getsize(marker) == 0:
        raise ValueError("snapshot is incomplete or from a stale pre-marker run: " +
                         snapshot)
    marker_values = {}
    with open(marker, "r") as stream:
        for line in stream:
            fields = line.strip().split("=", 1)
            if len(fields) == 2:
                marker_values[fields[0]] = fields[1]
    try:
        ranks = int(marker_values["ranks"])
    except (KeyError, ValueError):
        raise ValueError("snapshot completion marker has no valid rank count: " +
                         snapshot)
    required = ("manifest_rank", "fields_rank", "density_background_rank",
                "energy_spectrum_rank")
    for stem in required:
        paths = glob.glob(os.path.join(snapshot, stem + "*.dat"))
        if len(paths) != ranks or any(os.path.getsize(path) == 0 for path in paths):
            raise ValueError("snapshot completion marker/file mismatch for %s: %s" %
                             (stem, snapshot))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--performance-baseline",
                        help="controller-off run used only for timing/count gates")
    parser.add_argument("--performance-window-fs", type=parse_window)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()
    values = {}
    try:
        base = diagnostics(args.baseline)
        cand = diagnostics(args.candidate)
        perf_base = (diagnostics(args.performance_baseline)
                     if args.performance_baseline else base)
        values["baseline_accepted_step_count"] = len(base)
        values["candidate_accepted_step_count"] = len(cand)
        base_snapshots = [p for p in glob.glob(
            os.path.join(args.baseline, "snapshot_*")) if os.path.isdir(p)]
        candidate_snapshots = [p for p in glob.glob(
            os.path.join(args.candidate, "snapshot_*")) if os.path.isdir(p)]
        values["baseline_snapshot_count"] = len(base_snapshots)
        values["candidate_snapshot_count"] = len(candidate_snapshots)
        if not base:
            raise ValueError("baseline produced no accepted step: " +
                             args.baseline)
        if not cand:
            raise ValueError("candidate produced no accepted step: " +
                             args.candidate)
        if not base_snapshots:
            raise ValueError("baseline produced no snapshot: " +
                             args.baseline)
        if not candidate_snapshots:
            raise ValueError("candidate produced no snapshot: " +
                             args.candidate)
        bs, cs = latest_snapshot(args.baseline), latest_snapshot(args.candidate)
        require_complete_snapshot(bs)
        require_complete_snapshot(cs)
        if (len(base) != len(cand) or
                abs(base[0]["time_s"]-cand[0]["time_s"]) > 1.0e-24 or
                abs(base[-1]["time_s"]-cand[-1]["time_s"]) > 1.0e-24):
            raise ValueError("baseline/candidate accepted-step windows differ")
        density_tol = DENSITY_REL_L2_TOL
        field_tol = FIELD_REL_L2_TOL
        field_energy_tol = FIELD_ENERGY_REL_TOL
        spectrum_tol = SPECTRUM_REL_L2_TOL
        envelope_tol = WAVE_ENVELOPE_REL_L2_TOL
        base_density = rank_series(bs, "density_background", "n_e_m3")
        cand_density = rank_series(cs, "density_background", "n_e_m3")
        base_field = rank_series(bs, "fields", "Ex_Vm")
        cand_field = rank_series(cs, "fields", "Ex_Vm")
        density_l2 = relative_l2(base_density, cand_density)
        field_l2 = relative_l2(base_field, cand_field)
        field_energy_relative = squared_norm_relative_change(
            base_field, cand_field)
        field_envelope_l2 = relative_l2(
            block_rms(base_field, ENVELOPE_BLOCK_CELLS),
            block_rms(cand_field, ENVELOPE_BLOCK_CELLS))
        density_envelope_l2 = relative_l2(
            block_rms(base_density, ENVELOPE_BLOCK_CELLS),
            block_rms(cand_density, ENVELOPE_BLOCK_CELLS))
        spectrum_l2 = relative_l2(rank_sum_series(bs, "energy_spectrum", "count_combined"),
                                  rank_sum_series(cs, "energy_spectrum", "count_combined"))
        invariant_residual_names = ["tail_return_N_residual",
            "tail_return_Px_residual", "tail_return_K_residual"]
        representation_residual_names = ["tail_return_Jx_residual",
            "tail_return_Pixx_residual", "tail_return_Piperp_residual"]
        invariant_max = max_value(cand, invariant_residual_names)
        representation_max = max_value(cand, representation_residual_names)
        request_max = max_value(cand, ["tail_return_mpi_request_residual"])
        removed = max_value(cand, ["tail_return_particles_removed"])
        conservation_pass = (not unresolved_failure(args.candidate, cand) and
                             invariant_max <= 1.0e-12
                             and request_max <= 1.0e-13)
        spectrum_pass = spectrum_l2 <= spectrum_tol
        field_pass = (field_l2 <= field_tol and
                      field_energy_relative <= field_energy_tol)
        density_pass = density_l2 <= density_tol
        wave_pass = (field_envelope_l2 <= envelope_tol and
                     density_envelope_l2 <= envelope_tol)
        physical_base_tail = last_value(base, "tail_particle_count")
        base_tail = last_value(perf_base, "tail_particle_count")
        cand_tail = last_value(cand, "tail_particle_count")
        particle_reduction = (removed > 0.0 and math.isfinite(base_tail) and
                              math.isfinite(cand_tail) and cand_tail <= base_tail)
        base_wall = mean_kernel(perf_base, args.performance_window_fs)
        cand_wall = mean_kernel(cand, args.performance_window_fs)
        reduction = ((base_wall-cand_wall)/base_wall
                     if math.isfinite(base_wall) and base_wall > 0.0 else math.nan)
        performance_required = args.performance_window_fs is not None
        performance_pass = (not performance_required or
                            (particle_reduction and reduction >= 0.20))
        passed = (conservation_pass and spectrum_pass and field_pass and
                  density_pass and
                  wave_pass and particle_reduction and performance_pass)
        values.update({
            "comparison_valid": 1,
            "conservation_pass": int(conservation_pass),
            "spectrum_pass": int(spectrum_pass), "field_pass": int(field_pass),
            "density_pass": int(density_pass),
            "wave_envelope_pass": int(wave_pass),
            "particle_reduction_observed": int(particle_reduction),
            "performance_pass": int(performance_pass),
            "tail_return_invariant_residual_max": invariant_max,
            "tail_return_representation_residual_max": representation_max,
            "tail_return_request_residual_max": request_max,
            "density_relative_l2": density_l2, "field_relative_l2": field_l2,
            "field_energy_relative_difference": field_energy_relative,
            "field_envelope_relative_l2": field_envelope_l2,
            "density_envelope_relative_l2": density_envelope_l2,
            "spectrum_relative_l2": spectrum_l2,
            "density_relative_l2_tolerance": density_tol,
            "field_relative_l2_tolerance": field_tol,
            "field_energy_relative_tolerance": field_energy_tol,
            "spectrum_relative_l2_tolerance": spectrum_tol,
            "wave_envelope_relative_l2_tolerance": envelope_tol,
            "tail_particle_count_baseline": base_tail,
            "tail_particle_count_physical_baseline": physical_base_tail,
            "tail_particle_count_candidate": cand_tail,
            "performance_baseline_is_separate": int(
                args.performance_baseline is not None),
            "tail_kernel_wall_reduction": reduction,
        })
    except (IOError, OSError, ValueError) as exc:
        passed = False
        values["comparison_valid"] = 0
        values["comparison_error"] = str(exc)
    with open(args.result, "w") as out:
        for key in sorted(values):
            out.write("%s=%s\n" % (key, values[key]))
        out.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    print("status=%s" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
