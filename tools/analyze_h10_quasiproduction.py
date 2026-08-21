#!/usr/bin/env python3
"""Evaluate an H10 quasi-production run without changing solver state.

The report separates hard numerical/conservation gates from performance,
accepted-step energy accounting, and long-window macroscopic equivalence.
Older controller-off baselines may lack the new energy-ledger columns; only
the H10 run is required to contain them.
"""

from __future__ import print_function

import argparse
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tail_return_compare_common import (diagnostics, failed, latest_snapshot,
    last_value, max_value, rank_series, rank_sum_series, read_table,
    relative_l2)


INVARIANT_TOL = 1.0e-12
REQUEST_TOL = 1.0e-13
ENERGY_REL_TOL = 1.0e-3
FIELD_REL_L2_TOL = 2.0e-2
DENSITY_REL_L2_TOL = 1.0e-2
SPECTRUM_REL_L2_TOL = 2.0e-2
MOMENTUM_REL_L2_TOL = 2.0e-2
TAIL_REDUCTION_MIN = 0.20
WALL_REDUCTION_MIN = 0.20
PERFORMANCE_MIN_WINDOW_FS = 10.0


ENERGY_COLUMNS = (
    "domain_energy_before", "domain_energy_after", "domain_energy_delta",
    "accounted_energy_source", "energy_balance_residual",
    "energy_balance_relative", "electrostatic_boundary_work",
    "background_boundary_energy_net", "beam_boundary_energy_net")


def complete_snapshot(root):
    snapshot = latest_snapshot(root)
    marker = os.path.join(snapshot, "_COMPLETE")
    return os.path.isfile(marker) and os.path.getsize(marker) > 0


def mean(rows, name):
    return sum(row.get(name, 0.0) for row in rows) / max(1, len(rows))


def sum_column(rows, name):
    return sum(row.get(name, 0.0) for row in rows)


def energy_summary(rows, header):
    complete = all(column in header for column in ENERGY_COLUMNS)
    if not complete:
        return {"complete": False}
    delta = sum_column(rows, "domain_energy_delta")
    source = sum_column(rows, "accounted_energy_source")
    residual = sum_column(rows, "energy_balance_residual")
    scale = max(1.0, abs(delta) + abs(source))
    return {
        "complete": True,
        "delta": delta,
        "source": source,
        "residual": residual,
        "scale": scale,
        "relative": abs(residual) / scale,
    }


def rank_sum_columns(snapshot, stem, columns):
    totals = {}
    for path in glob.glob(os.path.join(snapshot, stem + "_rank*.dat")):
        header, rows = read_table(path)
        if any(column not in header for column in columns):
            raise ValueError("missing momentum column in " + path)
        x_name = header[0]
        for row in rows:
            x = row[x_name]
            totals[x] = totals.get(x, 0.0) + sum(row[c] for c in columns)
    if not totals:
        raise ValueError("no %s rank files under %s" % (stem, snapshot))
    return [totals[x] for x in sorted(totals)]


def macroscopic_comparison(reference_root, return_root):
    reference = latest_snapshot(reference_root)
    returned = latest_snapshot(return_root)
    reference_field = rank_series(reference, "fields", "Ex_Vm")
    returned_field = rank_series(returned, "fields", "Ex_Vm")
    reference_upar = rank_sum_columns(reference, "momentum_distribution_upar",
        ("dN_du_par_bg", "dN_du_par_tail", "dN_du_par_beam"))
    returned_upar = rank_sum_columns(returned, "momentum_distribution_upar",
        ("dN_du_par_bg", "dN_du_par_tail", "dN_du_par_beam"))
    # The current snapshot schema stores only the Eulerian perpendicular
    # marginal. Keep its name explicit; it is not a combined-tail metric.
    return {
        "field_relative_l2": relative_l2(reference_field, returned_field),
        "density_relative_l2": relative_l2(
            rank_series(reference, "density_background", "n_e_m3"),
            rank_series(returned, "density_background", "n_e_m3")),
        "spectrum_relative_l2": relative_l2(
            rank_sum_series(reference, "energy_spectrum", "count_combined"),
            rank_sum_series(returned, "energy_spectrum", "count_combined")),
        "upar_combined_relative_l2": relative_l2(reference_upar, returned_upar),
        "uperp_bulk_relative_l2": relative_l2(
            rank_sum_series(reference, "momentum_distribution_uperp",
                            "dN_du_perp_bg"),
            rank_sum_series(returned, "momentum_distribution_uperp",
                            "dN_du_perp_bg")),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--return-run", required=True)
    parser.add_argument("--baseline-run")
    parser.add_argument("--result", required=True)
    args = parser.parse_args()
    values = {}
    try:
        rows = diagnostics(args.return_run)
        if not rows:
            raise ValueError("return run has no accepted steps")
        header, _ = read_table(os.path.join(
            args.return_run, "vpfp_step_diagnostics.dat"))
        return_energy = energy_summary(rows, header)
        energy_complete = return_energy["complete"]
        values["accepted_step_count"] = len(rows)
        values["first_step"] = int(rows[0]["step"])
        values["last_step"] = int(rows[-1]["step"])
        values["first_time_s"] = rows[0]["time_s"]
        values["last_time_s"] = rows[-1]["time_s"]
        values["elapsed_time_fs"] = (
            values["last_time_s"] - values["first_time_s"]) * 1.0e15
        values["snapshot_complete"] = int(complete_snapshot(args.return_run))
        values["unresolved_failure"] = int(failed(args.return_run))
        values["split_count"] = int(sum_column(rows, "split"))
        values["collision_flux_rollback_count"] = int(sum_column(
            rows, "collision_flux_rollback_count"))
        values["return_invariant_residual_max"] = max_value(rows, (
            "tail_return_N_residual", "tail_return_Px_residual",
            "tail_return_K_residual"))
        values["return_request_residual_max"] = max_value(
            rows, ("tail_return_mpi_request_residual",))
        values["tail_particle_count_final"] = last_value(
            rows, "tail_particle_count")
        values["tail_particles_local_max_final"] = last_value(
            rows, "tail_particles_local_max")
        values["mean_wall_seconds"] = mean(rows, "wall_s")
        values["last100_mean_wall_seconds"] = mean(rows[-100:], "wall_s")
        values["return_particles_removed_total"] = sum_column(
            rows, "tail_return_particles_removed")
        values["return_number_total"] = sum_column(rows, "tail_return_N")
        values["return_representation_incompatible_cells_total"] = (
            int(sum_column(rows,
                "tail_return_projection_representation_incompatible_cells"))
            if "tail_return_projection_representation_incompatible_cells"
            in header else 0)
        for short, column in (
                ("N", "tail_return_N_difference"),
                ("Px", "tail_return_Px_difference"),
                ("Jx_dx", "tail_return_Jx_dx_difference"),
                ("K", "tail_return_K_difference"),
                ("Pixx_dx", "tail_return_Pixx_dx_difference"),
                ("Piperp_dx", "tail_return_Piperp_dx_difference")):
            values["return_%s_signed_difference_total" % short] = (
                sum_column(rows, column) if column in header else math.nan)
        piperp_added = sum_column(rows, "tail_return_Piperp_dx")
        piperp_difference = values[
            "return_Piperp_dx_signed_difference_total"]
        piperp_removed = piperp_added - piperp_difference
        values["return_Piperp_dx_signed_difference_relative"] = (
            abs(piperp_difference) /
            max(abs(piperp_removed), sys.float_info.min))

        values["energy_ledger_complete"] = int(energy_complete)
        if energy_complete:
            values["domain_energy_delta_total"] = return_energy["delta"]
            values["accounted_energy_source_total"] = return_energy["source"]
            values["electrostatic_boundary_work_total"] = sum_column(
                rows, "electrostatic_boundary_work")
            values["energy_balance_residual_total"] = return_energy["residual"]
            values["energy_balance_relative_cumulative"] = return_energy[
                "relative"]
            values["energy_balance_relative_step_max"] = max_value(
                rows, ("energy_balance_relative",))
        else:
            values["energy_balance_relative_cumulative"] = math.nan

        numerical_gate = (not failed(args.return_run) and
            values["snapshot_complete"] == 1 and
            values["split_count"] == 0 and
            values["collision_flux_rollback_count"] == 0)
        conservation_gate = (
            values["return_invariant_residual_max"] <= INVARIANT_TOL and
            values["return_request_residual_max"] <= REQUEST_TOL)
        absolute_energy_gate = (energy_complete and
            values["energy_balance_relative_cumulative"] <= ENERGY_REL_TOL)
        energy_gate = absolute_energy_gate
        values["numerical_gate_pass"] = int(numerical_gate)
        values["conservation_gate_pass"] = int(conservation_gate)
        values["absolute_energy_gate_pass"] = int(absolute_energy_gate)

        performance_gate = True
        physical_gate = True
        if args.baseline_run:
            baseline = diagnostics(args.baseline_run)
            if not baseline:
                raise ValueError("baseline run has no accepted steps")
            baseline_header, _ = read_table(os.path.join(
                args.baseline_run, "vpfp_step_diagnostics.dat"))
            baseline_energy = energy_summary(baseline, baseline_header)
            values["baseline_energy_ledger_complete"] = int(
                baseline_energy["complete"])
            if energy_complete and baseline_energy["complete"]:
                incremental_residual = (return_energy["residual"] -
                                        baseline_energy["residual"])
                incremental_scale = max(return_energy["scale"],
                                        baseline_energy["scale"], 1.0)
                incremental_relative = (abs(incremental_residual) /
                                        incremental_scale)
                values["baseline_energy_balance_residual_total"] = (
                    baseline_energy["residual"])
                values["baseline_energy_balance_relative_cumulative"] = (
                    baseline_energy["relative"])
                values["return_incremental_energy_residual_total"] = (
                    incremental_residual)
                values["return_incremental_energy_relative"] = (
                    incremental_relative)
                # H10 is responsible only for the A/B increment. The common
                # absolute energy defect remains an explicit production
                # blocker and is never hidden by this operator-level gate.
                energy_gate = incremental_relative <= ENERGY_REL_TOL
            baseline_tail = last_value(baseline, "tail_particle_count")
            baseline_wall = mean(baseline, "wall_s")
            tail_reduction = 1.0 - values["tail_particle_count_final"] / max(
                baseline_tail, 1.0)
            wall_reduction = 1.0 - values["mean_wall_seconds"] / max(
                baseline_wall, 1.0e-300)
            values["baseline_tail_particle_count_final"] = baseline_tail
            values["baseline_mean_wall_seconds"] = baseline_wall
            values["tail_particle_reduction"] = tail_reduction
            values["wall_time_reduction"] = wall_reduction
            performance_evaluated = (
                values["elapsed_time_fs"] >= PERFORMANCE_MIN_WINDOW_FS)
            values["performance_gate_evaluated"] = int(
                performance_evaluated)
            # A few-femtosecond A/B establishes numerical and physical
            # equivalence, but cannot measure the asymptotic Tail population
            # or amortized wall-time benefit. Keep those values informational
            # until the run spans a representative production window.
            performance_gate = (not performance_evaluated or
                (tail_reduction >= TAIL_REDUCTION_MIN and
                 wall_reduction >= WALL_REDUCTION_MIN))
            macro = macroscopic_comparison(args.baseline_run, args.return_run)
            values.update(macro)
            physical_gate = (
                macro["field_relative_l2"] <= FIELD_REL_L2_TOL and
                macro["density_relative_l2"] <= DENSITY_REL_L2_TOL and
                macro["spectrum_relative_l2"] <= SPECTRUM_REL_L2_TOL and
                macro["upar_combined_relative_l2"] <= MOMENTUM_REL_L2_TOL and
                macro["uperp_bulk_relative_l2"] <= MOMENTUM_REL_L2_TOL)
        values["performance_gate_pass"] = int(performance_gate)
        values["physical_gate_pass"] = int(physical_gate)
        values["energy_gate_pass"] = int(energy_gate)
        passed = (numerical_gate and conservation_gate and energy_gate and
                  performance_gate and physical_gate)
        values["production_ready"] = int(passed and absolute_energy_gate)
    except (IOError, OSError, ValueError, KeyError) as exc:
        passed = False
        values["analysis_error"] = str(exc)

    with open(args.result, "w") as stream:
        for key in sorted(values):
            stream.write("%s=%s\n" % (key, values[key]))
        stream.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    print("status=%s" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
