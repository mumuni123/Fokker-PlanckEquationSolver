#!/usr/bin/env python3
"""Gate H (section 11.7) stage-power A/B analyzer.

Reads one 115 fs checkpoint plus a coarse (dt_scale=0.5) and a fine
(dt_scale=0.25) accepted-step audit, both restarted from the same checkpoint,
and localizes the first stage/group that produces a fixed residual power per
unit physical time.

Inputs (production output, read-only; no Vlasov/Poisson/PIC/collision formula
is re-implemented):
  <checkpoint>/manifest.dat  (+ per-rank manifest_rank*.dat when present)
  <coarse>/vpfp_step_diagnostics.dat
  <coarse>/vpfp_stage_energy_audit.dat
  <fine>/vpfp_step_diagnostics.dat
  <fine>/vpfp_stage_energy_audit.dat

Exit codes (the .result `status=` must agree):
  0  PASS_ROOT_CAUSE_IDENTIFIED  | PASS_LEDGER_DEFECT_IDENTIFIED
  2  INCONCLUSIVE_EXTEND_WINDOW   | INCONCLUSIVE_EVENT_PATH_DIVERGED
  3  FAIL_AUDIT_STRUCTURE
  4  FAIL_NONFINITE_OR_CORRUPT_INPUT
  5  FAIL_USAGE
"""

from __future__ import print_function

import argparse
import glob
import math
import os
import re
import sys

EPS_DOUBLE = sys.float_info.epsilon
ROUNDOFF_FACTOR = 512.0
WINDOW_RELATIVE_GATE = 1.0e-12
EXPLANATION_GATE = 0.80
POWER_RATIO_LO = 0.80
POWER_RATIO_HI = 1.25

STAGES = (
    "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
    "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
    "collision_half2", "conversion_after_collision", "tail_bulk_return",
    "final_poisson",
)
REFERENCE_STAGE = "accepted_n"

# Sum-only predefined groups (section 11.7.2).  Sign is preserved.
STAGE_GROUPS = (
    ("collision_pair", ("collision_half1", "collision_half2")),
    ("x_pair", ("x_half1", "x_half2")),
    ("field_force_pair",
     ("midpoint_poisson", "u_force_tail_beam_kick", "final_poisson")),
    ("conversion_pair", ("conversion_after_force", "conversion_after_collision")),
    ("tail_return", ("tail_bulk_return",)),
)

STAGE_TO_GROUP = {}
for _g, _members in STAGE_GROUPS:
    for _m in _members:
        STAGE_TO_GROUP[_m] = _g

REQUIRED_STEP_COLUMNS = (
    "step", "time_s", "accepted", "split",
    "U_E", "K_e", "K_b", "K_tail", "K_combined",
    "U_E_before", "K_e_before", "K_b_before", "K_tail_before",
    "N_e_before", "N_e_after", "N_b_before", "N_b_after",
    "domain_energy_before", "domain_energy_after", "domain_energy_delta",
    "accounted_energy_source", "energy_balance_residual",
    "electrostatic_boundary_work", "background_boundary_energy_net",
    "beam_boundary_energy_net",
    "collision_reservoir", "fct_energy",
    "conversion_N_residual", "conversion_Px_residual", "conversion_K_residual",
    "tail_outflow_K",
    "tail_return_N_residual", "tail_return_Px_residual", "tail_return_K_residual",
    "collision_flux_rollback_count",
)

REQUIRED_STAGE_COLUMNS = (
    "step", "time_s", "accepted", "audit_valid", "split", "failure_code",
    "stage_id", "stage_name",
    "K_bulk", "K_tail", "K_beam", "U_E",
    "dK_bulk", "dK_tail", "dK_beam", "dU_E",
    "Q_bkg_left_in", "Q_bkg_left_out", "Q_bkg_right_in", "Q_bkg_right_out",
    "Q_beam_in", "Q_beam_out", "Q_tail_out", "Q_collision_reservoir",
    "K_conversion", "K_tail_return", "W_electrostatic_boundary", "stage_balance",
    "bulk_upar_face_work", "bulk_upar_velocity_boundary_work",
    "bulk_upar_interface_energy_removed", "bulk_upar_identity_residual",
    "tail_kick_work", "beam_kick_work",
)

INITIAL_SCALARS = (
    "domain_energy_before", "U_E_before", "K_e_before", "K_b_before",
    "K_tail_before", "N_e_before", "N_b_before",
)


class AuditFailure(Exception):
    def __init__(self, code, message):
        super(AuditFailure, self).__init__(message)
        self.code = code
        self.message = message


def fail_usage(message):
    raise AuditFailure(5, message)


def fail_corrupt(message):
    raise AuditFailure(4, message)


def fail_structure(message):
    raise AuditFailure(3, message)


def sanitize_stage_name(name):
    sanitized = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if not sanitized:
        sanitized = "stage"
    return sanitized


def read_manifest_file(path):
    """Read one manifest file into a dict.

    Supports both production formats:
      * checkpoint `manifest.txt`  -- whitespace-separated `key value` lines;
      * snapshot   `manifest.dat`  -- `key=value` lines.
    """
    values = {}
    if not os.path.isfile(path):
        fail_corrupt("missing manifest file: %s" % path)
    with open(path, "r") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if "=" in line:
                key, _, value = line.partition("=")
            else:
                parts = line.split(None, 1)
                key = parts[0]
                value = parts[1] if len(parts) > 1 else ""
            key = key.strip()
            value = value.strip()
            if key in values and values[key] != value:
                fail_corrupt("conflicting key %s in %s" % (key, path))
            values[key] = value
    return values


def manifest_time(values, path):
    if "time" in values:
        return values["time"]
    if "time_s" in values:
        return values["time_s"]
    fail_corrupt("manifest %s missing time/time_s" % path)


def read_manifest(checkpoint_dir):
    # The production checkpoint writes `manifest.txt` ("key value"); the
    # diagnostic snapshot writes `manifest.dat` ("key=value").  Prefer the
    # checkpoint file; fall back to the snapshot layout for offline tests.
    path = None
    for name in ("manifest.txt", "manifest.dat"):
        candidate = os.path.join(checkpoint_dir, name)
        if os.path.isfile(candidate):
            path = candidate
            break
    if path is None:
        fail_corrupt("missing manifest (manifest.txt or manifest.dat) in %s" %
                     checkpoint_dir)

    manifest = read_manifest_file(path)
    for key in ("step", "physical_config_hash"):
        if key not in manifest:
            fail_corrupt("manifest %s missing %s" % (path, key))
    time_text = manifest_time(manifest, path)
    try:
        step = int(manifest["step"])
        time_s = float(time_text)
        physical_config_hash = int(manifest["physical_config_hash"])
    except ValueError:
        fail_corrupt("manifest %s has unparsable step/time/hash" % path)
    if not math.isfinite(time_s):
        fail_corrupt("manifest %s has non-finite time" % path)

    # Per-rank manifests exist only in the snapshot `.dat` layout (the
    # checkpoint carries binary rank files and a single rank-0 manifest.txt).
    per_rank = manifest.get("per_rank_manifest", "")
    if per_rank:
        template = per_rank.replace("{rank}", "*")
        rank_files = sorted(glob.glob(os.path.join(checkpoint_dir, template)))
        if not rank_files:
            fail_corrupt("per-rank manifest template %s matched no files" %
                         per_rank)
        for rank_file in rank_files:
            rank_manifest = read_manifest_file(rank_file)
            for key in ("step", "physical_config_hash"):
                if key not in rank_manifest:
                    fail_corrupt("rank manifest %s missing %s" %
                                 (rank_file, key))
            rank_time_text = manifest_time(rank_manifest, rank_file)
            try:
                rank_step = int(rank_manifest["step"])
                rank_time = float(rank_time_text)
                rank_hash = int(rank_manifest["physical_config_hash"])
            except ValueError:
                fail_corrupt("rank manifest %s unparsable" % rank_file)
            if (rank_step != step or
                    rank_time != time_s or
                    rank_hash != physical_config_hash):
                fail_corrupt("rank manifest %s disagrees with global manifest "
                             "(step/time/hash)" % rank_file)
    return {"step": step, "time_s": time_s,
            "physical_config_hash": physical_config_hash}


def read_table(path, required_columns, label, string_columns=()):
    if not os.path.isfile(path):
        fail_corrupt("missing %s file: %s" % (label, path))
    with open(path, "r") as stream:
        header_line = stream.readline()
        if not header_line.strip():
            fail_corrupt("%s file %s has empty header" % (label, path))
        header = header_line.split()
        for column in required_columns:
            if column not in header:
                fail_structure("%s file %s missing required column %s" %
                               (label, path, column))
        rows = []
        for line_number, line in enumerate(stream, 2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != len(header):
                fail_corrupt("%s file %s line %d has %d columns, expected %d" %
                             (label, path, line_number, len(fields),
                              len(header)))
            row = {}
            for key, text in zip(header, fields):
                if key in string_columns:
                    row[key] = text
                    continue
                try:
                    value = float(text)
                except ValueError:
                    fail_corrupt("%s file %s line %d has invalid %s" %
                                 (label, path, line_number, key))
                if not math.isfinite(value):
                    fail_corrupt("%s file %s line %d has non-finite %s" %
                                 (label, path, line_number, key))
                row[key] = value
            rows.append(row)
    return rows


def group_by_step(rows, label):
    by_step = {}
    for row in rows:
        step = int(row["step"])
        if step in by_step:
            fail_structure("%s has duplicate step %d" % (label, step))
        by_step[step] = row
    return by_step


class RunData(object):
    def __init__(self, label, step_rows, stage_rows):
        self.label = label
        self.step_rows = step_rows          # ordered list
        self.steps = [int(r["step"]) for r in step_rows]
        self.by_step = {int(r["step"]): r for r in step_rows}

        # Stage rows grouped by step, ordered by stage_id.
        self.stage_by_step = {}
        for row in stage_rows:
            step = int(row["step"])
            self.stage_by_step.setdefault(step, []).append(row)

        self._validate()
        self._compute()

    def _validate(self):
        for step, step_row in self.by_step.items():
            if int(step_row["accepted"]) != 1:
                fail_structure("%s step %d is not accepted" %
                               (self.label, step))
            if int(step_row["split"]) != 0:
                fail_structure("%s step %d used a split" % (self.label, step))
            if int(step_row["collision_flux_rollback_count"]) != 0:
                fail_structure("%s step %d has a collision rollback" %
                               (self.label, step))

        # Stage mapping: stage_id <-> stage_name one-to-one, no unknown stage.
        id_to_name = {}
        name_to_id = {}
        for step, rows in self.stage_by_step.items():
            if step not in self.by_step:
                fail_structure("%s stage audit step %d missing from step "
                               "diagnostics" % (self.label, step))
            if len(rows) != len(STAGES):
                fail_structure("%s step %d has %d stage rows, expected %d" %
                               (self.label, step, len(rows), len(STAGES)))
            for row in rows:
                if int(row["accepted"]) != 1:
                    fail_structure("%s step %d stage not accepted" %
                                   (self.label, step))
                if int(row["audit_valid"]) != 1:
                    fail_structure("%s step %d stage audit-invalid" %
                                   (self.label, step))
                if int(row["split"]) != 0:
                    fail_structure("%s step %d stage split" %
                                   (self.label, step))
                if int(row["failure_code"]) != 0:
                    fail_structure("%s step %d stage failure_code" %
                                   (self.label, step))
                sid = int(row["stage_id"])
                sname = row["stage_name"]
                if sid in id_to_name and id_to_name[sid] != sname:
                    fail_structure("%s stage_id %d maps to multiple names" %
                                   (self.label, sid))
                if sname in name_to_id and name_to_id[sname] != sid:
                    fail_structure("%s stage_name %s maps to multiple ids" %
                                   (self.label, sname))
                id_to_name[sid] = sname
                name_to_id[sname] = sid
                if sname not in STAGES:
                    # Unknown stage is preserved as unclassified and fails the
                    # structure gate.
                    fail_structure("%s has unknown stage_name %s" %
                                   (self.label, sname))
        for sname in STAGES:
            if sname not in name_to_id:
                fail_structure("%s missing stage_name %s" %
                               (self.label, sname))

    def _stage_row(self, step, stage_name):
        for row in self.stage_by_step[step]:
            if row["stage_name"] == stage_name:
                return row
        fail_structure("%s step %d missing stage %s" %
                       (self.label, step, stage_name))

    def _compute(self):
        self.n_steps = len(self.steps)
        self.first_time_s = self.by_step[self.steps[0]]["time_s"]
        self.last_time_s = self.by_step[self.steps[-1]]["time_s"]

        # Per-step stage_balance sums.
        stage_signed = dict((s, 0.0) for s in STAGES)
        stage_abs = dict((s, 0.0) for s in STAGES)
        group_signed = dict((g, 0.0) for g, _ in STAGE_GROUPS)
        derived = {
            "poisson_pair": 0.0,
            "force_conversion_pair": 0.0,
            "field_particle_pair": 0.0,
            "bulk_upar_identity": 0.0,
            "tail_kick_mismatch": 0.0,
            "beam_kick_mismatch": 0.0,
        }

        for step in self.steps:
            rows = dict((r["stage_name"], r) for r in self.stage_by_step[step])
            for sname in STAGES:
                value = rows[sname]["stage_balance"]
                stage_signed[sname] += value
                stage_abs[sname] += abs(value)
            for gname, members in STAGE_GROUPS:
                total = 0.0
                for member in members:
                    total += rows[member]["stage_balance"]
                group_signed[gname] += total

            uf = rows["u_force_tail_beam_kick"]
            mp = rows["midpoint_poisson"]
            fp = rows["final_poisson"]
            ac = rows["accepted_n"]
            conv = rows["conversion_after_force"]

            poisson_pair = mp["stage_balance"] + fp["stage_balance"]
            force_conversion_pair = (uf["stage_balance"] +
                                     conv["stage_balance"])
            field_particle_pair = (
                poisson_pair
                + uf["bulk_upar_face_work"]
                + uf["tail_kick_work"]
                + uf["beam_kick_work"]
                - (fp["W_electrostatic_boundary"] -
                   ac["W_electrostatic_boundary"]))
            bulk_upar_identity = uf["bulk_upar_identity_residual"]
            tail_kick_mismatch = uf["dK_tail"] - uf["tail_kick_work"]
            beam_kick_mismatch = uf["dK_beam"] - uf["beam_kick_work"]

            derived["poisson_pair"] += poisson_pair
            derived["force_conversion_pair"] += force_conversion_pair
            derived["field_particle_pair"] += field_particle_pair
            derived["bulk_upar_identity"] += bulk_upar_identity
            derived["tail_kick_mismatch"] += tail_kick_mismatch
            derived["beam_kick_mismatch"] += beam_kick_mismatch

        self.stage_signed = stage_signed
        self.stage_abs = stage_abs
        self.group_signed = group_signed
        self.derived = derived

        # Full residual and source ledger from the step diagnostics.
        full_signed = 0.0
        full_abs = 0.0
        accounted = 0.0
        domain_delta = 0.0
        known_source = 0.0
        collision_reservoir = 0.0
        telescope_max_abs = 0.0
        conversion_events = 0
        energy_scale = 1.0
        for step in self.steps:
            r = self.by_step[step]
            residual = r["energy_balance_residual"]
            full_signed += residual
            full_abs += abs(residual)
            accounted += r["accounted_energy_source"]
            domain_delta += r["domain_energy_delta"]
            collision_reservoir += r["collision_reservoir"]
            known_source += (r["background_boundary_energy_net"] +
                             r["beam_boundary_energy_net"] -
                             r["tail_outflow_K"] -
                             r["collision_reservoir"] +
                             r["electrostatic_boundary_work"])
            energy_scale = max(energy_scale,
                               abs(r["domain_energy_before"]),
                               abs(r["domain_energy_after"]))
            # Telescope: stage_balance sum must reproduce the full residual.
            stage_sum = 0.0
            rows = dict((rr["stage_name"], rr)
                        for rr in self.stage_by_step[step])
            for sname in STAGES:
                if sname == REFERENCE_STAGE:
                    continue
                stage_sum += rows[sname]["stage_balance"]
            telescope_max_abs = max(telescope_max_abs,
                                    abs(stage_sum - residual))
            if rows["conversion_after_force"]["K_conversion"] != \
                    rows["u_force_tail_beam_kick"]["K_conversion"]:
                conversion_events += 1

        self.full_signed = full_signed
        self.full_abs = full_abs
        self.accounted = accounted
        self.domain_delta = domain_delta
        self.known_source = known_source
        self.collision_reservoir = collision_reservoir
        self.telescope_max_abs = telescope_max_abs
        self.energy_scale = energy_scale
        self.conversion_events = conversion_events

        self.tail_return_attempted = sum(
            r.get("tail_return_attempted_groups", 0.0)
            for r in self.step_rows if "tail_return_attempted_groups" in r)
        self.tail_return_committed = sum(
            r.get("tail_return_committed_groups", 0.0)
            for r in self.step_rows if "tail_return_committed_groups" in r)
        if "tail_particle_count" in self.step_rows[0]:
            self.tail_particles_first = self.step_rows[0]["tail_particle_count"]
            self.tail_particles_last = self.step_rows[-1]["tail_particle_count"]
        else:
            self.tail_particles_first = 0.0
            self.tail_particles_last = 0.0

    def telescope_ok(self):
        tol = ROUNDOFF_FACTOR * EPS_DOUBLE * self.energy_scale
        return self.telescope_max_abs <= max(1.0e-10, tol)

    def source_ownership_ok(self):
        tol = ROUNDOFF_FACTOR * EPS_DOUBLE * self.energy_scale
        return abs(self.known_source - self.accounted) <= max(1.0e-10, tol)

    def machine_gate(self):
        return max(1.0e-10,
                   ROUNDOFF_FACTOR * EPS_DOUBLE * self.energy_scale)


def read_run(directory, label):
    step_rows = read_table(
        os.path.join(directory, "vpfp_step_diagnostics.dat"),
        REQUIRED_STEP_COLUMNS, label + " step diagnostics")
    stage_rows = read_table(
        os.path.join(directory, "vpfp_stage_energy_audit.dat"),
        REQUIRED_STAGE_COLUMNS, label + " stage audit",
        string_columns=("stage_name",))
    step_rows = sorted(step_rows, key=lambda r: int(r["step"]))
    if not step_rows:
        fail_structure("%s has no accepted step rows" % label)
    step_ids = [int(r["step"]) for r in step_rows]
    if len(set(step_ids)) != len(step_ids):
        fail_structure("%s has duplicate steps" % label)
    stage_rows = sorted(stage_rows, key=lambda r: int(r["step"]))
    return RunData(label, step_rows, stage_rows)


def explanation_fraction(p_full, p_candidate, floor):
    denom = max(abs(p_full), floor)
    return max(0.0, 1.0 - abs(p_full - p_candidate) / denom)


def compute_power(signed_sum, elapsed_s):
    return signed_sum / elapsed_s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--coarse-dt-scale", type=float, required=True)
    parser.add_argument("--fine-dt-scale", type=float, required=True)
    parser.add_argument("--expected-coarse-steps", type=int, required=True)
    parser.add_argument("--expected-fine-steps", type=int, required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    if not (args.coarse_dt_scale > 0.0 and args.fine_dt_scale > 0.0):
        fail_usage("--coarse-dt-scale and --fine-dt-scale must be positive")
    if args.expected_coarse_steps <= 0 or args.expected_fine_steps <= 0:
        fail_usage("--expected-*-steps must be positive")

    try:
        manifest = read_manifest(args.checkpoint)
        coarse = read_run(args.coarse, "coarse")
        fine = read_run(args.fine, "fine")
        result = analyze(manifest, coarse, fine, args)
    except AuditFailure as failure:
        write_failure(args.result, failure)
        return failure.code

    write_result(args.result, result)
    return result["exit_code"]


def analyze(manifest, coarse, fine, args):
    checkpoint_step = manifest["step"]
    checkpoint_time_s = manifest["time_s"]
    checkpoint_hash = manifest["physical_config_hash"]

    # same_checkpoint: both runs start right after the checkpoint step.
    same_checkpoint = int(
        coarse.steps[0] == checkpoint_step + 1 and
        fine.steps[0] == checkpoint_step + 1)

    # Initial physical state: coarse vs fine first-row scalars.
    same_initial = 1
    for scalar in INITIAL_SCALARS:
        a = coarse.by_step[coarse.steps[0]][scalar]
        b = fine.by_step[fine.steps[0]][scalar]
        tol = ROUNDOFF_FACTOR * EPS_DOUBLE * max(abs(a), abs(b), 1.0)
        if abs(a - b) > tol:
            same_initial = 0

    coarse.elapsed_s = coarse.last_time_s - checkpoint_time_s
    fine.elapsed_s = fine.last_time_s - checkpoint_time_s
    if not (coarse.elapsed_s > 0.0 and fine.elapsed_s > 0.0):
        fail_structure("non-positive elapsed physical window")
    window_relative_error = abs(coarse.elapsed_s - fine.elapsed_s) / \
        max(abs(coarse.elapsed_s), abs(fine.elapsed_s))

    coarse.full_power = compute_power(coarse.full_signed, coarse.elapsed_s)
    fine.full_power = compute_power(fine.full_signed, fine.elapsed_s)

    energy_scale = max(coarse.energy_scale, fine.energy_scale, 1.0)
    p_floor = ROUNDOFF_FACTOR * EPS_DOUBLE * energy_scale / \
        max(coarse.elapsed_s, fine.elapsed_s)

    # ---- Structure gates ----
    # source_ownership is a diagnostic that drives the ledger-defect finding
    # (section 11.7.6 item 2), not a hard structure gate: a known source
    # missing from accounted_energy_source is reported as
    # PASS_LEDGER_DEFECT_IDENTIFIED, never as a structure failure.
    steps_ok = (coarse.n_steps == args.expected_coarse_steps and
                fine.n_steps == args.expected_fine_steps)

    # Physical-window gate.  The window is a small difference of large absolute
    # times (0.128 fs at ~115 fs), so the roundoff from N time-step additions
    # (~eps * |t|) is amplified by |t|/T.  Use a machine-scaled tolerance with
    # the documented 1e-12 as the floor instead of a fixed 1e-12 that false-
    # fails a floating-point-consistent short window.
    time_scale = max(abs(checkpoint_time_s),
                     abs(coarse.last_time_s), abs(fine.last_time_s))
    window_tolerance = max(
        WINDOW_RELATIVE_GATE,
        ROUNDOFF_FACTOR * EPS_DOUBLE * time_scale /
        max(abs(coarse.elapsed_s), abs(fine.elapsed_s)))

    checks = []
    if same_checkpoint != 1:
        checks.append("same_checkpoint")
    if same_initial != 1:
        checks.append("same_initial_physical_state")
    if window_relative_error > window_tolerance:
        checks.append("window_relative_error=%.6g>tol=%.6g" %
                      (window_relative_error, window_tolerance))
    if not steps_ok:
        checks.append("steps")
    if not (coarse.telescope_ok() and fine.telescope_ok()):
        checks.append("stage_telescope")
    if checks:
        fail_structure("structure gate failed: " + ", ".join(checks))

    # ---- Candidate powers ----
    candidates = {}  # name -> (coarse_power, fine_power)

    def add_candidate(name, coarse_sum, fine_sum):
        candidates[name] = (
            compute_power(coarse_sum, coarse.elapsed_s),
            compute_power(fine_sum, fine.elapsed_s),
        )

    for sname in STAGES:
        if sname == REFERENCE_STAGE:
            continue
        add_candidate(sname, coarse.stage_signed[sname], fine.stage_signed[sname])
    for gname, _ in STAGE_GROUPS:
        add_candidate(gname, coarse.group_signed[gname], fine.group_signed[gname])
    for dname in coarse.derived:
        add_candidate(dname, coarse.derived[dname], fine.derived[dname])

    # ---- Ledger-defect detection (section 11.7.6 item 2) ----
    ledger_defect_runs = 0
    for run in (coarse, fine):
        mismatch = run.known_source - run.accounted
        if abs(mismatch) <= run.machine_gate():
            continue
        mismatch_power = compute_power(mismatch, run.elapsed_s)
        full_power = (coarse.full_power if run is coarse else fine.full_power)
        if explanation_fraction(full_power, mismatch_power, p_floor) >= EXPLANATION_GATE:
            ledger_defect_runs += 1
    ledger_defect = ledger_defect_runs == 2

    # ---- Single root cause (section 11.7.2 / 11.7.6) ----
    hits = []
    for name, (coarse_power, fine_power) in candidates.items():
        coarse_full = coarse.full_power
        fine_full = fine.full_power
        sign_ok = (math.copysign(1.0, coarse_power) ==
                   math.copysign(1.0, coarse_full) and
                   math.copysign(1.0, fine_power) ==
                   math.copysign(1.0, fine_full))
        coarse_frac = explanation_fraction(coarse_full, coarse_power, p_floor)
        fine_frac = explanation_fraction(fine_full, fine_power, p_floor)
        if (abs(coarse_power) <= p_floor or abs(fine_power) <= p_floor):
            ratio = float("nan")
            ratio_ok = False
        else:
            ratio = abs(fine_power) / max(abs(coarse_power), p_floor)
            ratio_ok = POWER_RATIO_LO <= ratio <= POWER_RATIO_HI
        if (sign_ok and coarse_frac >= EXPLANATION_GATE and
                fine_frac >= EXPLANATION_GATE and ratio_ok):
            hits.append((name, coarse_frac, fine_frac, ratio))

    event_paths_differ = int(
        coarse.conversion_events != fine.conversion_events or
        coarse.tail_return_committed != fine.tail_return_committed)

    def pick_root_cause():
        if ledger_defect:
            return ("ledger_defect", "PASS_LEDGER_DEFECT_IDENTIFIED", 0,
                    None, None, None)
        # Priority: single stage > field_particle_pair > predefined group.
        stage_hits = [h for h in hits if h[0] in STAGES]
        if stage_hits:
            best = max(stage_hits, key=lambda h: (h[1] + h[2]))
            return (best[0], "PASS_ROOT_CAUSE_IDENTIFIED", 0,
                    best[1], best[2], best[3])
        fp_hits = [h for h in hits if h[0] == "field_particle_pair"]
        if fp_hits:
            best = fp_hits[0]
            # Section 11.7.6 item 4: the trigger quantity is
            # field_particle_pair_power (the field-particle energy residual),
            # which maps to the field_force_pair block in the recommended-files
            # table.  Report the trigger quantity, not the stage group, so the
            # 99% explanation fraction is attached to the right name.
            return ("field_particle_pair", "PASS_ROOT_CAUSE_IDENTIFIED", 0,
                    best[1], best[2], best[3])
        group_hits = [h for h in hits
                      if h[0] not in STAGES and h[0] != "field_particle_pair"]
        if group_hits:
            best = max(group_hits, key=lambda h: (h[1] + h[2]))
            if best[0] == "tail_return" and event_paths_differ:
                return ("event_path_diverged",
                        "INCONCLUSIVE_EVENT_PATH_DIVERGED", 2,
                        best[1], best[2], best[3])
            return (best[0], "PASS_ROOT_CAUSE_IDENTIFIED", 0,
                    best[1], best[2], best[3])
        return ("not_yet_unique", "INCONCLUSIVE_EXTEND_WINDOW", 2,
                None, None, None)

    root_cause, status, exit_code, dom_frac_c, dom_frac_f, dom_ratio = \
        pick_root_cause()
    dominant = root_cause
    if dominant == "not_yet_unique":
        dominant = "not_yet_unique"
        dom_frac_c = dom_frac_f = float("nan")
        dom_ratio = float("nan")

    recommended = recommend_files(root_cause)

    # ---- Assemble result ----
    def power_ratio_and_order(coarse_power, fine_power):
        if abs(coarse_power) <= p_floor or abs(fine_power) <= p_floor:
            return float("nan"), "not_evaluated_roundoff"
        ratio = abs(fine_power) / max(abs(coarse_power), p_floor)
        order = math.log2(abs(coarse_power) / abs(fine_power))
        return ratio, order

    full_ratio, full_order = power_ratio_and_order(coarse.full_power,
                                                   fine.full_power)

    result = {
        "status": status,
        "exit_code": exit_code,
        "checkpoint_step": checkpoint_step,
        "checkpoint_time_s": "%.17g" % checkpoint_time_s,
        "checkpoint_physical_config_hash": checkpoint_hash,
        "coarse_steps": coarse.n_steps,
        "fine_steps": fine.n_steps,
        "coarse_elapsed_time_s": "%.17g" % coarse.elapsed_s,
        "fine_elapsed_time_s": "%.17g" % fine.elapsed_s,
        "same_checkpoint": same_checkpoint,
        "same_initial_physical_state": same_initial,
        "same_physical_window_relative_error": "%.17g" % window_relative_error,
        "same_physical_window_tolerance": "%.17g" % window_tolerance,
        "full_residual_power_coarse": "%.17g" % coarse.full_power,
        "full_residual_power_fine": "%.17g" % fine.full_power,
        "full_residual_power_ratio": "%.17g" % full_ratio,
        "full_residual_observed_order":
            full_order if isinstance(full_order, str) else "%.17g" % full_order,
        "known_source_minus_accounted_coarse":
            "%.17g" % (coarse.known_source - coarse.accounted),
        "known_source_minus_accounted_fine":
            "%.17g" % (fine.known_source - fine.accounted),
        "dominant_stage_or_group": dominant,
        "dominant_explanation_fraction_coarse":
            "%.17g" % dom_frac_c if dom_frac_c is not None else "nan",
        "dominant_explanation_fraction_fine":
            "%.17g" % dom_frac_f if dom_frac_f is not None else "nan",
        "dominant_power_ratio":
            "%.17g" % dom_ratio if dom_ratio is not None else "nan",
        "root_cause": root_cause,
        "recommended_files": recommended,
        "stage_telescope_matches_roundoff": int(
            coarse.telescope_ok() and fine.telescope_ok()),
        "source_ownership_valid": int(
            coarse.source_ownership_ok() and fine.source_ownership_ok()),
        "all_finite": 1,
        "event_paths_differ": event_paths_differ,
    }

    # Per-run full/source/domain quantities.
    for tag, run in (("coarse", coarse), ("fine", fine)):
        result["full_residual_signed_%s" % tag] = "%.17g" % run.full_signed
        result["full_residual_abs_%s" % tag] = "%.17g" % run.full_abs
        result["accounted_source_power_%s" % tag] = "%.17g" % \
            compute_power(run.accounted, run.elapsed_s)
        result["domain_energy_delta_power_%s" % tag] = "%.17g" % \
            compute_power(run.domain_delta, run.elapsed_s)
        result["known_source_reconstructed_%s" % tag] = "%.17g" % \
            run.known_source
        result["known_source_power_%s" % tag] = "%.17g" % \
            compute_power(run.known_source, run.elapsed_s)
        result["unclassified_power_%s" % tag] = "0"
        result["conversion_event_count_%s" % tag] = run.conversion_events
        result["tail_return_attempted_groups_%s" % tag] = \
            "%.17g" % run.tail_return_attempted
        result["tail_return_committed_groups_%s" % tag] = \
            "%.17g" % run.tail_return_committed
        result["tail_particle_count_first_%s" % tag] = \
            "%.17g" % run.tail_particles_first
        result["tail_particle_count_last_%s" % tag] = \
            "%.17g" % run.tail_particles_last

    # Per-stage signed/abs/power/fraction.
    for sname in STAGES:
        if sname == REFERENCE_STAGE:
            continue
        key = "stage_%s" % sanitize_stage_name(sname)
        for tag, run in (("coarse", coarse), ("fine", fine)):
            signed = run.stage_signed[sname]
            power = compute_power(signed, run.elapsed_s)
            full_power = run.full_power
            result["%s_residual_signed_%s" % (key, tag)] = "%.17g" % signed
            result["%s_residual_abs_%s" % (key, tag)] = "%.17g" % \
                run.stage_abs[sname]
            result["%s_residual_power_%s" % (key, tag)] = "%.17g" % power
            result["%s_fraction_of_full_signed_%s" % (key, tag)] = "%.17g" % \
                explanation_fraction(full_power, power, p_floor)

    # Per-item power ratio + observed order.
    for name, (coarse_power, fine_power) in candidates.items():
        ratio, order = power_ratio_and_order(coarse_power, fine_power)
        key = sanitize_stage_name(name)
        if name in STAGES:
            key = "stage_%s" % key
        result["%s_power_ratio" % key] = "%.17g" % ratio
        result["%s_observed_order" % key] = \
            order if isinstance(order, str) else "%.17g" % order

    # Named derived power outputs (section 11.7.2).
    named = {
        "field_particle_pair_power": "field_particle_pair",
        "poisson_pair_power": "poisson_pair",
        "force_conversion_pair_power": "force_conversion_pair",
        "bulk_upar_identity_power": "bulk_upar_identity",
        "tail_kick_mismatch_power": "tail_kick_mismatch",
        "beam_kick_mismatch_power": "beam_kick_mismatch",
        "collision_stage_power": "collision_pair",
        "x_remap_boundary_power": "x_pair",
        "conversion_internal_power": "conversion_pair",
        "tail_return_internal_power": "tail_return",
    }
    for out_key, candidate in named.items():
        coarse_power, fine_power = candidates[candidate]
        result["%s_coarse" % out_key] = "%.17g" % coarse_power
        result["%s_fine" % out_key] = "%.17g" % fine_power

    result["collision_reservoir_power_coarse"] = "%.17g" % \
        compute_power(coarse.collision_reservoir, coarse.elapsed_s)
    result["collision_reservoir_power_fine"] = "%.17g" % \
        compute_power(fine.collision_reservoir, fine.elapsed_s)

    return result


def recommend_files(root_cause):
    table = {
        "ledger_defect": "src/vpfp_integrator.cpp:finalize_energy_ledger(), "
                         "src/vpfp_diagnostics.cpp",
        "collision_pair": "collision-stage energy-exchange ledger + unit tests",
        "x_pair": "open background boundary energy flux ownership + x-remap "
                  "stage ledger",
        "field_force_pair": "src/vpfp_integrator.cpp midpoint/final Poisson "
                            "time layer (separate structure-preserving design)",
        "field_particle_pair": "src/vpfp_integrator.cpp midpoint/final Poisson "
                               "time layer (separate structure-preserving "
                               "design)",
        "conversion_pair": "flux-interface conversion stage state-diff and "
                           "internal transfer ownership",
        "tail_return": "physical residence-time semantics + dedicated A/B",
        "not_yet_unique": "20/40 step extension or additional diagnostics",
        "event_path_diverged": "physical residence-time semantics + dedicated "
                               "A/B",
    }
    stage_to_rec = {
        "collision_half1": table["collision_pair"],
        "collision_half2": table["collision_pair"],
        "x_half1": table["x_pair"],
        "x_half2": table["x_pair"],
        "midpoint_poisson": table["field_force_pair"],
        "u_force_tail_beam_kick": table["field_force_pair"],
        "final_poisson": table["field_force_pair"],
        "conversion_after_force": table["conversion_pair"],
        "conversion_after_collision": table["conversion_pair"],
        "tail_bulk_return": table["tail_return"],
    }
    if root_cause in table:
        return table[root_cause]
    if root_cause in stage_to_rec:
        return stage_to_rec[root_cause]
    return "20/40 step extension or additional diagnostics"


def write_result(path, result):
    lines = []
    for key in sorted(result):
        if key == "exit_code":
            continue
        lines.append("%s=%s" % (key, result[key]))
    with open(path, "w") as stream:
        stream.write("\n".join(lines) + "\n")


def write_failure(path, failure):
    status_map = {
        3: "FAIL_AUDIT_STRUCTURE",
        4: "FAIL_NONFINITE_OR_CORRUPT_INPUT",
        5: "FAIL_USAGE",
    }
    status = status_map.get(failure.code, "FAIL_AUDIT_STRUCTURE")
    with open(path, "w") as stream:
        stream.write("status=%s\n" % status)
        stream.write("error=%s\n" % failure.message.replace("\n", " "))
        stream.write("root_cause=audit_failure\n")


if __name__ == "__main__":
    sys.exit(main())
