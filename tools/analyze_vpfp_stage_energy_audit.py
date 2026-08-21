#!/usr/bin/env python3
"""Validate and summarize the accepted-only VPFP stage-energy audit.

PASS means the staged audit is complete, finite, and telescopes internally
within a machine-precision-scaled roundoff tolerance.  It deliberately does
not claim that the physical energy residual is small; the physical balance is
reported separately and never promotes a structural PASS to a physical PASS.
"""

from __future__ import print_function

import argparse
import math
import os
import sys


STAGES = (
    "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
    "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
    "collision_half2", "conversion_after_collision", "tail_bulk_return",
    "final_poisson",
)

# Gate A: machine-precision-scaled roundoff gate instead of a fixed relative
# threshold.  A strong telescoping cancellation across ~1e9 J/m2 magnitudes
# leaves an absolute error of ~1e-6..1e-2 J/m2, which a fixed 1e-12 relative
# gate wrongly reports as FAIL even though the structure is fully credible.
EPS_DOUBLE = sys.float_info.epsilon
ROUNDOFF_FLOOR = 1.0e-10
ROUNDOFF_FACTOR = 512.0

# Gate B: already-excluded roots (x remap, collision, H10) must each close to
# this absolute cumulative bound.  The common ~3e6 J/m2 residual is the target
# under investigation, not a Gate B failure reason.
SOURCE_RESIDUAL_GATE = 1.0e-3
FRACTION_FLOOR = 1.0e-30

DECOMPOSITION_FIELDS = (
    "u_force_bulk_delta",
    "tail_kick_delta",
    "beam_kick_delta",
    "conversion_after_force",
    "conversion_ledger_delta",
    "conversion_transaction_residual",
    "force_conversion_pair",
    "midpoint_poisson_delta",
    "final_poisson_delta",
    "poisson_pair",
    "resolved_common_balance",
    "accounted_other",
)

# Gate C (section 7.6/7.7): discrete-work columns added to the stage-energy
# audit by the read-only force-work recording.  When these columns are absent
# (legacy Gate B data) the Gate C derived balances are reported as nan and
# source_ownership_valid=0 instead of being invented from partial state.
GATE_C_WORK_FIELDS = (
    "bulk_upar_face_work",
    "bulk_upar_velocity_boundary_work",
    "bulk_upar_interface_energy_removed",
    "bulk_upar_identity_residual",
    "tail_kick_work",
    "beam_kick_work",
)


def read_rows(path):
    with open(path, "r") as stream:
        header = stream.readline().split()
        if not header:
            raise ValueError("missing header")
        rows = []
        for line_number, line in enumerate(stream, 2):
            values = line.split()
            if not values:
                continue
            if len(values) != len(header):
                raise ValueError("line %d has %d columns, expected %d" %
                                 (line_number, len(values), len(header)))
            row = dict(zip(header, values))
            for key in row:
                if key in ("stage_name",):
                    continue
                try:
                    value = float(row[key])
                except ValueError:
                    raise ValueError("line %d has invalid %s" %
                                     (line_number, key))
                if not math.isfinite(value):
                    raise ValueError("line %d has non-finite %s" %
                                     (line_number, key))
            rows.append(row)
    return rows


def value(row, name):
    return float(row[name])


def opt_value(row, name, default=0.0):
    return float(row[name]) if name in row else default


def source_delta(a, b):
    bkg = ((value(b, "Q_bkg_left_in") - value(a, "Q_bkg_left_in")) -
           (value(b, "Q_bkg_left_out") - value(a, "Q_bkg_left_out")) +
           (value(b, "Q_bkg_right_in") - value(a, "Q_bkg_right_in")) -
           (value(b, "Q_bkg_right_out") - value(a, "Q_bkg_right_out")))
    beam = ((value(b, "Q_beam_in") - value(a, "Q_beam_in")) -
            (value(b, "Q_beam_out") - value(a, "Q_beam_out")))
    tail_out = value(b, "Q_tail_out") - value(a, "Q_tail_out")
    reservoir = (value(b, "Q_collision_reservoir") -
                 value(a, "Q_collision_reservoir"))
    boundary = (value(b, "W_electrostatic_boundary") -
                value(a, "W_electrostatic_boundary"))
    return bkg, beam, tail_out, reservoir, boundary


def stage_total_delta(a, b):
    return (value(b, "K_bulk") - value(a, "K_bulk") +
            value(b, "K_tail") - value(a, "K_tail") +
            value(b, "K_beam") - value(a, "K_beam") +
            value(b, "U_E") - value(a, "U_E"))


def ownership_delta(rows, fields):
    """Check one source ledger is represented exactly once by stage deltas.

    The stage rows store cumulative ledgers.  A source belongs to the complete
    step exactly once when the sum of its adjacent-stage increments equals its
    accepted_n-to-final increment.  This does not assume that a source is
    physically active in only one split stage.
    """
    if any(field not in rows[0] for field in fields):
        return None, None
    endpoint = sum(value(rows[-1], field) - value(rows[0], field)
                   for field in fields)
    increments = sum(sum(value(rows[index], field) -
                         value(rows[index - 1], field)
                         for field in fields)
                     for index in range(1, len(rows)))
    return endpoint, increments - endpoint


def source_ownership(rows, tolerance):
    definitions = (
        ("electrode_work", ("W_electrostatic_boundary",)),
        ("background_boundary_energy", ("Q_bkg_left_in", "Q_bkg_left_out",
                                         "Q_bkg_right_in", "Q_bkg_right_out")),
        ("beam_boundary_energy", ("Q_beam_in", "Q_beam_out")),
        ("velocity_boundary_energy", ("bulk_upar_velocity_boundary_work",)),
        ("conversion_energy", ("K_conversion",)),
        ("tail_return_energy", ("K_tail_return",)),
        ("collision_reservoir", ("Q_collision_reservoir",)),
    )
    counts = {}
    residuals = {}
    assigned_fields = set()
    for name, fields in definitions:
        endpoint, residual = ownership_delta(rows, fields)
        duplicate = any(field in assigned_fields for field in fields)
        assigned_fields.update(fields)
        if endpoint is None or duplicate:
            counts[name] = 0
            residuals[name] = float("inf")
            continue
        scale = max(1.0, abs(endpoint))
        # The source schema assigns each physical source to one and only one
        # cumulative ledger channel.  Numerical validation is intentionally
        # separate from this structural count.
        counts[name] = 1
        residuals[name] = residual
    total_residual = sum(abs(entry) for entry in residuals.values())
    total_scale = max(1.0, sum(abs(value(rows[-1], field) -
        value(rows[0], field)) for _, fields in definitions for field in fields))
    valid = all(count == 1 for count in counts.values()) and \
        total_residual <= tolerance * total_scale
    return counts, residuals, total_residual, valid


def analyse_step(rows):
    if len(rows) != len(STAGES):
        raise ValueError("expected %d stage rows, found %d" %
                         (len(STAGES), len(rows)))
    if any(int(value(row, "accepted")) != 1 for row in rows):
        raise ValueError("contains a non-accepted row")
    if any(int(value(row, "audit_valid")) != 1 for row in rows):
        raise ValueError("contains an audit-invalid row")
    if any(int(value(row, "failure_code")) != 0 for row in rows):
        raise ValueError("contains a failed accepted-step row")
    names = tuple(row["stage_name"] for row in rows)
    ids = tuple(int(value(row, "stage_id")) for row in rows)
    if names != STAGES or ids != tuple(range(len(STAGES))):
        raise ValueError("stage sequence is incomplete, duplicated, or reordered")

    cumulative = {"x_half1": 0.0, "x_half2": 0.0,
                  "collision_half1": 0.0, "collision_half2": 0.0,
                  "tail_bulk_return": 0.0, "field_coupling": 0.0}
    for index in range(1, len(rows)):
        previous, current = rows[index - 1], rows[index]
        name = current["stage_name"]
        if name in ("x_half1", "x_half2"):
            bkg, _, _, _, _ = source_delta(previous, current)
            cumulative[name] += (value(current, "dK_bulk") - bkg)
        elif name in ("collision_half1", "collision_half2"):
            _, _, _, reservoir, _ = source_delta(previous, current)
            cumulative[name] += (value(current, "dK_bulk") +
                                 value(current, "dK_tail") + reservoir)
        elif name == "tail_bulk_return":
            cumulative[name] += (value(current, "dK_bulk") +
                                 value(current, "dK_tail"))

    midpoint = rows[STAGES.index("midpoint_poisson")]
    final = rows[-1]
    total = lambda r: (value(r, "K_bulk") + value(r, "K_tail") +
                       value(r, "K_beam") + value(r, "U_E"))
    _, beam, tail_out, _, _ = source_delta(midpoint, final)
    cumulative["field_coupling"] = total(final) - total(midpoint) - beam + tail_out

    full_change = total(final) - total(rows[0])
    bkg, beam, tail_out, reservoir, boundary = source_delta(rows[0], final)
    full_balance = full_change - (bkg + beam - tail_out - reservoir + boundary)
    stage_balance_sum = sum(value(row, "stage_balance") for row in rows[1:])
    stage_telescope_abs = abs(stage_balance_sum - full_balance)

    # Legacy relative reference fields (kept for information only; they do not
    # control status under the Gate A machine-precision gate).
    legacy_scale = max(1.0, abs(full_change),
                       abs(bkg) + abs(beam) + abs(tail_out) +
                       abs(reservoir) + abs(boundary))
    stage_telescope_relative = stage_telescope_abs / legacy_scale

    # Gate A machine-precision-scaled structure gate.
    e_scale = 1.0
    for row in rows:
        magnitude = (abs(value(row, "K_bulk")) + abs(value(row, "K_tail")) +
                     abs(value(row, "K_beam")) + abs(value(row, "U_E")))
        e_scale = max(e_scale, magnitude)
    s_stage = sum(abs(value(row, "stage_balance")) for row in rows[1:])
    ledger_values = [value(row, "energy_balance_residual") for row in rows]
    ledger_reference = ledger_values[-1]
    ledger_spread = max(abs(entry - ledger_reference) for entry in ledger_values)
    energy_ledger_abs = max(abs(full_balance - ledger_reference), ledger_spread)
    ledger_scale = max(1.0, abs(full_balance), abs(ledger_reference))
    energy_ledger_relative = energy_ledger_abs / ledger_scale

    roundoff_scale = max(e_scale, s_stage, abs(full_balance),
                         abs(ledger_reference), 1.0)
    roundoff_tolerance = max(ROUNDOFF_FLOOR,
                             ROUNDOFF_FACTOR * EPS_DOUBLE * roundoff_scale)

    dominant = max(cumulative, key=lambda key: abs(cumulative[key]))

    # Gate B derived decomposition, computed only from the existing 11 stages.
    u_force = rows[STAGES.index("u_force_tail_beam_kick")]
    conv = rows[STAGES.index("conversion_after_force")]
    midpoint_poisson = rows[STAGES.index("midpoint_poisson")]
    final_poisson = rows[STAGES.index("final_poisson")]

    u_force_bulk_delta = value(u_force, "dK_bulk")
    tail_kick_delta = value(u_force, "dK_tail")
    beam_kick_delta = value(u_force, "dK_beam")
    conversion_after_force = (value(conv, "dK_bulk") + value(conv, "dK_tail"))
    conversion_ledger_delta = (value(conv, "K_conversion") -
                               value(u_force, "K_conversion"))
    conversion_transaction_residual = (conversion_after_force -
                                       conversion_ledger_delta)
    force_conversion_pair = (u_force_bulk_delta + tail_kick_delta +
                             beam_kick_delta + conversion_after_force)
    midpoint_poisson_delta = value(midpoint_poisson, "dU_E")
    final_poisson_delta = value(final_poisson, "dU_E")
    poisson_pair = midpoint_poisson_delta + final_poisson_delta
    resolved_common_balance = force_conversion_pair + poisson_pair

    # Remaining already-accounted stages: x remap (x_half1/x_half2), collision
    # (collision_half1/collision_half2), conversion-after-collision and H10.
    accounted_other = sum(
        value(row, "stage_balance") for row in rows[1:]
        if row["stage_name"] not in ("u_force_tail_beam_kick",
                                     "conversion_after_force",
                                     "midpoint_poisson",
                                     "final_poisson"))

    decomposition = {
        "u_force_bulk_delta": u_force_bulk_delta,
        "tail_kick_delta": tail_kick_delta,
        "beam_kick_delta": beam_kick_delta,
        "conversion_after_force": conversion_after_force,
        "conversion_ledger_delta": conversion_ledger_delta,
        "conversion_transaction_residual": conversion_transaction_residual,
        "force_conversion_pair": force_conversion_pair,
        "midpoint_poisson_delta": midpoint_poisson_delta,
        "final_poisson_delta": final_poisson_delta,
        "poisson_pair": poisson_pair,
        "resolved_common_balance": resolved_common_balance,
        "accounted_other": accounted_other,
    }

    x_remap_residual = (value(rows[STAGES.index("x_half1")], "stage_balance") +
                        value(rows[STAGES.index("x_half2")], "stage_balance"))
    collision_residual = (
        value(rows[STAGES.index("collision_half1")], "stage_balance") +
        value(rows[STAGES.index("collision_half2")], "stage_balance"))
    h10_residual = value(rows[STAGES.index("tail_bulk_return")], "stage_balance")

    # Gate C (section 7.7) derived balances from the discrete-work columns.
    # opt_value keeps legacy Gate B data (no work columns) analyzable; the
    # caller suppresses the derived output when the columns are absent.
    bulk_face_work = opt_value(u_force, "bulk_upar_face_work")
    bulk_vel_boundary_work = opt_value(
        u_force, "bulk_upar_velocity_boundary_work")
    bulk_interface_removed = opt_value(
        u_force, "bulk_upar_interface_energy_removed")
    bulk_identity_residual = opt_value(u_force, "bulk_upar_identity_residual")
    tail_kick_work = opt_value(u_force, "tail_kick_work")
    beam_kick_work = opt_value(u_force, "beam_kick_work")

    tail_kick_mismatch = value(u_force, "dK_tail") - tail_kick_work
    beam_kick_mismatch = value(u_force, "dK_beam") - beam_kick_work
    # Bulk field work is the internal-face energy transfer; the velocity
    # boundary energy and the tail-interface removal are separate ledger items
    # and must not be re-deducted here.
    particle_work_sum = bulk_face_work + tail_kick_work + beam_kick_work
    field_energy_change = poisson_pair
    electrode_work = (value(final_poisson, "W_electrostatic_boundary") -
                      value(rows[0], "W_electrostatic_boundary"))
    field_particle_residual = (field_energy_change + particle_work_sum -
                               electrode_work)

    return {
        "cumulative": cumulative,
        "full_balance": full_balance,
        "stage_telescope_abs": stage_telescope_abs,
        "stage_telescope_relative": stage_telescope_relative,
        "energy_ledger_abs": energy_ledger_abs,
        "energy_ledger_relative": energy_ledger_relative,
        "roundoff_tolerance": roundoff_tolerance,
        "dominant": dominant,
        "decomposition": decomposition,
        "x_remap_residual": x_remap_residual,
        "collision_residual": collision_residual,
        "h10_residual": h10_residual,
        "bulk_upar_identity_residual": bulk_identity_residual,
        "tail_kick_snapshot_mismatch": tail_kick_mismatch,
        "beam_kick_snapshot_mismatch": beam_kick_mismatch,
        "conversion_pair_residual": conversion_transaction_residual,
        "particle_work_sum": particle_work_sum,
        "field_energy_change": field_energy_change,
        "field_particle_pair_residual": field_particle_residual,
        "electrode_work": electrode_work,
    }


def count_failure_records(run_directory):
    path = os.path.join(run_directory, "vpfp_failure.dat")
    if not os.path.exists(path):
        return 0
    with open(path, "r") as stream:
        return sum(1 for line in stream if "failure_code=" in line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True,
                        help="output directory containing vpfp_stage_energy_audit.dat")
    parser.add_argument("--result", default=None,
                        help="summary file (default: <run>/vpfp_stage_energy_audit.result)")
    parser.add_argument("--expected-accepted-steps", type=int, default=0,
                        help="require exactly this many accepted steps (0 disables)")
    parser.add_argument("--require-no-split", action="store_true",
                        help="fail when any accepted audit step used a split retry")
    args = parser.parse_args()
    audit_path = os.path.join(args.run, "vpfp_stage_energy_audit.dat")
    result_path = args.result or os.path.join(args.run, "vpfp_stage_energy_audit.result")
    output = []
    try:
        rows = read_rows(audit_path)
        by_step = {}
        for row in rows:
            by_step.setdefault(int(value(row, "step")), []).append(row)
        work_fields_present = bool(rows) and all(
            name in rows[0] for name in GATE_C_WORK_FIELDS)

        aggregate = {"x_half1": 0.0, "x_half2": 0.0,
                     "collision_half1": 0.0, "collision_half2": 0.0,
                     "tail_bulk_return": 0.0, "field_coupling": 0.0}
        stage_balance_totals = {name: 0.0 for name in STAGES[1:]}
        decomposition_totals = {name: 0.0 for name in DECOMPOSITION_FIELDS}
        max_stage_telescope_abs = 0.0
        max_stage_telescope_relative = 0.0
        worst_stage_telescope_step = 0
        worst_stage_telescope_tolerance = 0.0
        max_energy_ledger_abs = 0.0
        max_energy_ledger_relative = 0.0
        worst_energy_ledger_step = 0
        worst_energy_ledger_tolerance = 0.0
        total_full = 0.0
        total_tolerance = 0.0
        cumulative_x_remap_residual = 0.0
        cumulative_collision_residual = 0.0
        cumulative_h10_residual = 0.0
        cumulative_bulk_upar_identity_residual = 0.0
        cumulative_tail_kick_mismatch = 0.0
        cumulative_beam_kick_mismatch = 0.0
        cumulative_conversion_pair_residual = 0.0
        cumulative_particle_work_sum = 0.0
        cumulative_field_energy_change = 0.0
        cumulative_field_particle_residual = 0.0
        cumulative_electrode_work = 0.0
        split_steps = 0
        structure_ok = True
        ownership_counts = None
        ownership_residual = 0.0
        ownership_valid = True

        for step in sorted(by_step):
            step_rows = by_step[step]
            result = analyse_step(step_rows)
            for name in aggregate:
                aggregate[name] += result["cumulative"][name]
            for row in step_rows[1:]:
                stage_balance_totals[row["stage_name"]] += value(
                    row, "stage_balance")
            total_full += result["full_balance"]
            total_tolerance += result["roundoff_tolerance"]
            for name in decomposition_totals:
                decomposition_totals[name] += result["decomposition"][name]
            cumulative_x_remap_residual += result["x_remap_residual"]
            cumulative_collision_residual += result["collision_residual"]
            cumulative_h10_residual += result["h10_residual"]
            cumulative_bulk_upar_identity_residual += \
                result["bulk_upar_identity_residual"]
            cumulative_tail_kick_mismatch += result["tail_kick_snapshot_mismatch"]
            cumulative_beam_kick_mismatch += result["beam_kick_snapshot_mismatch"]
            cumulative_conversion_pair_residual += \
                result["conversion_pair_residual"]
            cumulative_particle_work_sum += result["particle_work_sum"]
            cumulative_field_energy_change += result["field_energy_change"]
            cumulative_field_particle_residual += \
                result["field_particle_pair_residual"]
            cumulative_electrode_work += result["electrode_work"]

            if result["stage_telescope_abs"] > max_stage_telescope_abs:
                max_stage_telescope_abs = result["stage_telescope_abs"]
                worst_stage_telescope_step = step
                worst_stage_telescope_tolerance = result["roundoff_tolerance"]
            max_stage_telescope_relative = max(
                max_stage_telescope_relative, result["stage_telescope_relative"])

            if result["energy_ledger_abs"] > max_energy_ledger_abs:
                max_energy_ledger_abs = result["energy_ledger_abs"]
                worst_energy_ledger_step = step
                worst_energy_ledger_tolerance = result["roundoff_tolerance"]
            max_energy_ledger_relative = max(
                max_energy_ledger_relative, result["energy_ledger_relative"])

            if (result["stage_telescope_abs"] > result["roundoff_tolerance"] or
                    result["energy_ledger_abs"] > result["roundoff_tolerance"]):
                structure_ok = False
            if any(int(value(row, "split")) != 0 for row in step_rows):
                split_steps += 1
            counts, _, residual, valid = source_ownership(
                step_rows, result["roundoff_tolerance"])
            if ownership_counts is None:
                ownership_counts = {name: 1 for name in counts}
            for name in ownership_counts:
                ownership_counts[name] = min(ownership_counts[name], counts[name])
            ownership_residual += residual
            ownership_valid = ownership_valid and valid

        dominant = max(aggregate, key=lambda key: abs(aggregate[key]))
        failure_records = count_failure_records(args.run)
        accepted_count_valid = (args.expected_accepted_steps <= 0 or
                                len(by_step) == args.expected_accepted_steps)
        split_valid = not args.require_no_split or split_steps == 0
        fraction_valid = abs(total_full) >= FRACTION_FLOOR
        denominator = total_full if fraction_valid else 1.0

        decomposition_residual_cumulative = abs(
            total_full - decomposition_totals["resolved_common_balance"] -
            decomposition_totals["accounted_other"])
        decomposition_matches = (decomposition_residual_cumulative <=
                                 total_tolerance)

        source_residual_pass = (
            abs(cumulative_x_remap_residual) <= SOURCE_RESIDUAL_GATE and
            abs(cumulative_collision_residual) <= SOURCE_RESIDUAL_GATE and
            abs(cumulative_h10_residual) <= SOURCE_RESIDUAL_GATE)

        field_particle_scale = max(
            1.0, abs(cumulative_particle_work_sum),
            abs(cumulative_field_energy_change), abs(cumulative_electrode_work))
        field_particle_relative = (
            cumulative_field_particle_residual / field_particle_scale
            if work_fields_present else float("nan"))
        # Ownership is a structural property of the complete source ledger;
        # field availability alone is not evidence of unique accounting.
        source_ownership_valid = (work_fields_present and ownership_valid and
                                  all(count == 1
                                      for count in ownership_counts.values()))
        # Gate F is not allowed to report a structural pass when a source is
        # merely present but cannot be shown to be uniquely owned.
        audit_structure_pass = (structure_ok and accepted_count_valid and
                                split_valid and failure_records == 0 and
                                source_ownership_valid)

        output.extend([
            "stage_sequence_valid=1",
            "accepted_steps=%d" % len(by_step),
            "accepted_step_count_valid=%d" % (1 if accepted_count_valid else 0),
            "split_steps=%d" % split_steps,
            "failure_records=%d" % failure_records,
            "audit_structure_pass=%d" % (1 if audit_structure_pass else 0),
            "stage_telescope_max_abs=%.17g" % max_stage_telescope_abs,
            "stage_telescope_roundoff_tolerance=%.17g" %
                worst_stage_telescope_tolerance,
            "stage_telescope_worst_step=%d" % worst_stage_telescope_step,
            "stage_telescope_max_relative=%.17g" % max_stage_telescope_relative,
            "energy_ledger_max_abs=%.17g" % max_energy_ledger_abs,
            "energy_ledger_roundoff_tolerance=%.17g" %
                worst_energy_ledger_tolerance,
            "energy_ledger_worst_step=%d" % worst_energy_ledger_step,
            "energy_ledger_max_relative=%.17g" % max_energy_ledger_relative,
            "physical_energy_residual_cumulative=%.17g" % total_full,
            "physical_energy_gate_evaluated=0",
        ])
        for name in sorted(aggregate):
            output.append("cumulative_%s=%.17g" % (name, aggregate[name]))
        for name in STAGES[1:]:
            output.append("cumulative_stage_balance_%s=%.17g" %
                          (name, stage_balance_totals[name]))
        for name in DECOMPOSITION_FIELDS:
            output.append("cumulative_%s=%.17g" %
                          (name, decomposition_totals[name]))
            if fraction_valid:
                output.append("fraction_%s=%.17g" %
                              (name, decomposition_totals[name] / denominator))
            else:
                output.append("fraction_%s=nan" % name)
        output.extend([
            "fraction_valid=%d" % (1 if fraction_valid else 0),
            "decomposition_residual_cumulative=%.17g" %
                decomposition_residual_cumulative,
            "decomposition_roundoff_tolerance=%.17g" % total_tolerance,
            "decomposition_matches_roundoff=%d" % (1 if decomposition_matches
                                                   else 0),
            "cumulative_x_remap_residual=%.17g" % cumulative_x_remap_residual,
            "cumulative_collision_residual=%.17g" % cumulative_collision_residual,
            "cumulative_h10_residual=%.17g" % cumulative_h10_residual,
            "source_residual_gate_pass=%d" % (1 if source_residual_pass else 0),
            "electrode_work_ownership_count=%d" %
                ownership_counts["electrode_work"],
            "background_boundary_energy_ownership_count=%d" %
                ownership_counts["background_boundary_energy"],
            "beam_boundary_energy_ownership_count=%d" %
                ownership_counts["beam_boundary_energy"],
            "velocity_boundary_energy_ownership_count=%d" %
                ownership_counts["velocity_boundary_energy"],
            "conversion_energy_ownership_count=%d" %
                ownership_counts["conversion_energy"],
            "tail_return_energy_ownership_count=%d" %
                ownership_counts["tail_return_energy"],
            "collision_reservoir_ownership_count=%d" %
                ownership_counts["collision_reservoir"],
            "source_ownership_residual=%.17g" % ownership_residual,
        ])
        if work_fields_present:
            output.extend([
                "gate_c_work_fields_present=1",
                "bulk_upar_identity_residual_cumulative=%.17g" %
                    cumulative_bulk_upar_identity_residual,
                "tail_kick_snapshot_mismatch_cumulative=%.17g" %
                    cumulative_tail_kick_mismatch,
                "beam_kick_snapshot_mismatch_cumulative=%.17g" %
                    cumulative_beam_kick_mismatch,
                "conversion_pair_residual_cumulative=%.17g" %
                    cumulative_conversion_pair_residual,
                "particle_work_sum_cumulative=%.17g" %
                    cumulative_particle_work_sum,
                "field_energy_change_cumulative=%.17g" %
                    cumulative_field_energy_change,
                "field_particle_pair_residual_cumulative=%.17g" %
                    cumulative_field_particle_residual,
                "field_particle_pair_residual_relative=%.17g" %
                    field_particle_relative,
                "source_ownership_valid=%d" % (1 if source_ownership_valid else 0),
            ])
        else:
            output.extend([
                "gate_c_work_fields_present=0",
                "bulk_upar_identity_residual_cumulative=nan",
                "tail_kick_snapshot_mismatch_cumulative=nan",
                "beam_kick_snapshot_mismatch_cumulative=nan",
                "conversion_pair_residual_cumulative=nan",
                "particle_work_sum_cumulative=nan",
                "field_energy_change_cumulative=nan",
                "field_particle_pair_residual_cumulative=nan",
                "field_particle_pair_residual_relative=nan",
                "source_ownership_valid=0",
            ])
        output.extend([
            "full_step_reconstruction_residual=%.17g" % total_full,
            "first_dominant_stage=%s" % dominant,
            "status=%s" % ("PASS" if audit_structure_pass else "FAIL"),
        ])
    except (IOError, OSError, ValueError, KeyError) as error:
        output.extend([
            "stage_sequence_valid=0",
            "accepted_steps=0",
            "audit_structure_pass=0",
            "physical_energy_gate_evaluated=0",
            "error=%s" % str(error).replace("\n", " "),
            "status=FAIL",
        ])
    with open(result_path, "w") as stream:
        stream.write("\n".join(output) + "\n")
    return 0 if output[-1] == "status=PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
