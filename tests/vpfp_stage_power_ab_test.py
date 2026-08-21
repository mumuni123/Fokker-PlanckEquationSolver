#!/usr/bin/env python3
"""Gate H (section 11.7.1) unit test for analyze_vpfp_stage_power_ab.py.

Builds minimal valid step/stage tables in temporary directories and drives the
production analyzer as a subprocess.  No production data file is required.

Usage:
  python3 tests/vpfp_stage_power_ab_test.py \
      --analyzer ./tools/analyze_vpfp_stage_power_ab.py \
      --result ./output/vpfp_stage_power_ab_unit.result
The last stdout line is always "status=PASS" or "status=FAIL".
"""

from __future__ import print_function

import argparse
import os
import shutil
import subprocess
import sys
import uuid


_TMP_BASE = os.path.join(os.getcwd(), "output", ".gate_h_test_tmp")


def mkdtemp(prefix):
    """Workspace-local temp dir (system temp may be unwritable in sandbox)."""
    if not os.path.isdir(_TMP_BASE):
        os.makedirs(_TMP_BASE, exist_ok=True)
    path = os.path.join(_TMP_BASE, prefix + uuid.uuid4().hex[:10])
    os.makedirs(path)
    return path

STAGES = (
    "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
    "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
    "collision_half2", "conversion_after_collision", "tail_bulk_return",
    "final_poisson",
)

STEP_COLS = [
    "step", "time_s", "accepted", "split",
    "U_E", "K_e", "K_b", "K_tail", "K_combined",
    "U_E_before", "K_e_before", "K_b_before", "K_tail_before",
    "N_e_before", "N_e_after", "N_b_before", "N_b_after",
    "domain_energy_before", "domain_energy_after", "domain_energy_delta",
    "accounted_energy_source", "energy_balance_residual",
    "electrostatic_boundary_work", "background_boundary_energy_net",
    "beam_boundary_energy_net", "collision_reservoir", "fct_energy",
    "conversion_N_residual", "conversion_Px_residual", "conversion_K_residual",
    "tail_outflow_K",
    "tail_return_N_residual", "tail_return_Px_residual", "tail_return_K_residual",
    "collision_flux_rollback_count",
    "tail_return_attempted_groups", "tail_return_committed_groups",
    "tail_particle_count",
]

STAGE_COLS = [
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
]


def write_manifest(directory, step, time_s, physical_hash):
    # Production checkpoint layout: manifest.txt, whitespace-separated
    # "key value" lines with `time` (not `time_s`).
    with open(os.path.join(directory, "manifest.txt"), "w") as out:
        out.write("solver_kind vpfp-open-v4\n")
        out.write("version 5\n")
        out.write("step %d\n" % step)
        out.write("time %.17g\n" % time_s)
        out.write("physical_config_hash %d\n" % physical_hash)


def write_step_file(directory, n_steps, first_step, times, residual,
                    accounted, sources, extra=None, drop_col=None):
    """sources: dict of column -> per-step value (list or scalar)."""
    path = os.path.join(directory, "vpfp_step_diagnostics.dat")
    header = STEP_COLS if drop_col is None else [c for c in STEP_COLS
                                                 if c != drop_col]
    with open(path, "w") as out:
        out.write(" ".join(header) + "\n")
        for i in range(n_steps):
            step = first_step + i
            t = times[i]
            r = residual[i] if isinstance(residual, list) else residual
            a = accounted[i] if isinstance(accounted, list) else accounted
            row = dict((c, 0.0) for c in header)
            row["step"] = step
            row["time_s"] = t
            row["accepted"] = 1
            row["split"] = 0
            row["domain_energy_before"] = 5.0e9
            row["domain_energy_after"] = 5.0e9
            row["domain_energy_delta"] = 0.0
            row["accounted_energy_source"] = a
            row["energy_balance_residual"] = r
            row["collision_flux_rollback_count"] = 0
            row["tail_particle_count"] = 1.0e6 + i
            for key, value in sources.items():
                if key not in row:
                    continue
                row[key] = value[i] if isinstance(value, list) else value
            if extra:
                for key, value in extra.items():
                    if key in row:
                        row[key] = value[i] if isinstance(value, list) else value
            out.write(" ".join(str(row[c]) for c in header) + "\n")


def write_stage_file(directory, n_steps, first_step, times, stage_balances,
                     residual, drop_col=None):
    """stage_balances: dict stage_name -> per-step value (list or scalar)."""
    path = os.path.join(directory, "vpfp_stage_energy_audit.dat")
    header = STAGE_COLS if drop_col is None else [c for c in STAGE_COLS
                                                  if c != drop_col]
    with open(path, "w") as out:
        out.write(" ".join(header) + "\n")
        for i in range(n_steps):
            step = first_step + i
            t = times[i]
            for sid, sname in enumerate(STAGES):
                r = residual[i] if isinstance(residual, list) else residual
                value = 0.0
                if sname in stage_balances:
                    value = stage_balances[sname]
                    if isinstance(value, list):
                        value = value[i]
                row = dict((c, 0.0) for c in header)
                row["step"] = step
                row["time_s"] = t
                row["accepted"] = 1
                row["audit_valid"] = 1
                row["split"] = 0
                row["failure_code"] = 0
                row["stage_id"] = sid
                row["stage_name"] = sname
                row["stage_balance"] = value
                out.write(" ".join(str(row[c]) for c in header) + "\n")


def run_analyzer(analyzer, checkpoint, coarse, fine, expected_c, expected_f,
                 result_path, coarse_dt=0.5, fine_dt=0.25):
    command = [
        sys.executable, analyzer,
        "--checkpoint", checkpoint,
        "--coarse", coarse,
        "--fine", fine,
        "--coarse-dt-scale", str(coarse_dt),
        "--fine-dt-scale", str(fine_dt),
        "--expected-coarse-steps", str(expected_c),
        "--expected-fine-steps", str(expected_f),
        "--result", result_path,
    ]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    _, _ = proc.communicate()
    return proc.returncode


def read_result(path):
    result = {}
    with open(path, "r") as stream:
        for line in stream:
            line = line.strip()
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            result[key.strip()] = value.strip()
    return result


class Case:
    def __init__(self, root):
        self.root = root
        self.checkpoint = os.path.join(root, "checkpoint")
        self.coarse = os.path.join(root, "coarse")
        self.fine = os.path.join(root, "fine")
        for directory in (self.checkpoint, self.coarse, self.fine):
            os.makedirs(directory)
        # checkpoint step = 100; runs start at 101.
        write_manifest(self.checkpoint, 100, 1.0e-13, 123456789)


def make_times(n_steps, last_time, start=1.0e-13):
    # Monotonic times ending exactly at last_time; only first/last are used.
    times = []
    for i in range(n_steps):
        frac = (i + 1.0) / float(n_steps)
        times.append(start + frac * (last_time - start))
    return times


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analyzer", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    analyzer = os.path.abspath(args.analyzer)
    if not os.path.isfile(analyzer):
        print("status=FAIL")
        return 2

    results = []

    # Case 1: single stage explains 90%, same power -> PASS_ROOT_CAUSE.
    results.append(("case1_single_stage", case1(analyzer)))

    # Case 2: known source missing from accounted_energy_source -> ledger.
    results.append(("case2_ledger_defect", case2(analyzer)))

    # Case 3: max explanation below 80% -> INCONCLUSIVE_EXTEND_WINDOW.
    results.append(("case3_not_unique", case3(analyzer)))

    # Case 4: inconsistent physical window -> FAIL_AUDIT_STRUCTURE.
    results.append(("case4_window", case4(analyzer)))

    # Case 5: hard failures (missing col, duplicate step, nonfinite, split,
    # rollback).
    results.append(("case5a_missing_col", case5a(analyzer)))
    results.append(("case5b_duplicate_step", case5b(analyzer)))
    results.append(("case5c_nonfinite", case5c(analyzer)))
    results.append(("case5d_split", case5d(analyzer)))
    results.append(("case5e_rollback", case5e(analyzer)))

    # Case 6: roundoff-level stage -> order=not_evaluated_roundoff.
    results.append(("case6_roundoff_order", case6(analyzer)))

    # Case 7: 2r/r per-step same power -> observed order near 0.
    results.append(("case7_zero_order", case7(analyzer)))

    all_pass = all(ok for _, ok in results)
    with open(args.result, "w") as out:
        for name, ok in results:
            out.write("%s=%d\n" % (name, 1 if ok else 0))
        out.write("status=%s\n" % ("PASS" if all_pass else "FAIL"))
    for name, ok in results:
        print("%s=%s" % (name, "PASS" if ok else "FAIL"))
    print("status=%s" % ("PASS" if all_pass else "FAIL"))
    return 0 if all_pass else 1


def case1(analyzer):
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        stage_balances = {
            "final_poisson": 0.9,   # will be scaled by residual below
            "x_half1": 0.1,
        }
        # build_simple_case scales stage_balances?  No: pass absolute values.
        # Use per-step absolute values via list form.
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {
            "final_poisson": [0.9 * 2.0 * r] * 10,
            "x_half1": [0.1 * 2.0 * r] * 10,
        }
        # Fine uses half the per-step residual.
        sb_fine = {
            "final_poisson": [0.9 * r] * 20,
            "x_half1": [0.1 * r] * 20,
        }
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_fine, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        ok = (code == 0 and res.get("status") == "PASS_ROOT_CAUSE_IDENTIFIED"
              and res.get("root_cause") == "final_poisson"
              and res.get("dominant_stage_or_group") == "final_poisson")
        return ok
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case2(analyzer):
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e5
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        # Residual entirely from a boundary flux omitted from the ledger.
        coarse_residual = [r] * 10
        fine_residual = [r * 0.5] * 20
        sb_c = {"final_poisson": [r] * 10}
        sb_f = {"final_poisson": [r * 0.5] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10,
                        {"background_boundary_energy_net": [r] * 10})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20,
                        {"background_boundary_energy_net": [r * 0.5] * 20})
        write_stage_file(case.coarse, 10, 101, times_c, sb_c, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        ok = (code == 0 and
              res.get("status") == "PASS_LEDGER_DEFECT_IDENTIFIED" and
              res.get("root_cause") == "ledger_defect" and
              res.get("source_ownership_valid") == "0")
        return ok
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case3(analyzer):
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        # Spread across four stages of four different groups: no single stage
        # or predefined group reaches 80%.
        sb_c = {
            "collision_half1": [0.25 * 2.0 * r] * 10,
            "x_half1": [0.25 * 2.0 * r] * 10,
            "midpoint_poisson": [0.25 * 2.0 * r] * 10,
            "conversion_after_force": [0.25 * 2.0 * r] * 10,
        }
        sb_f = {
            "collision_half1": [0.25 * r] * 20,
            "x_half1": [0.25 * r] * 20,
            "midpoint_poisson": [0.25 * r] * 20,
            "conversion_after_force": [0.25 * r] * 20,
        }
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb_c, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        ok = (code == 2 and
              res.get("status") == "INCONCLUSIVE_EXTEND_WINDOW" and
              res.get("root_cause") == "not_yet_unique")
        return ok
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case4(analyzer):
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        # Coarse and fine end at different physical times.
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.02e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb_c = {"final_poisson": [0.9 * 2.0 * r] * 10,
                "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb_c, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        ok = (code == 3 and res.get("status") == "FAIL_AUDIT_STRUCTURE")
        return ok
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case5a(analyzer):  # missing column
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {"final_poisson": [0.9 * 2.0 * r] * 10,
              "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {}, drop_col="accounted_energy_source")
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        return code == 3 and res.get("status") == "FAIL_AUDIT_STRUCTURE"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case5b(analyzer):  # duplicate step
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {"final_poisson": [0.9 * 2.0 * r] * 10,
              "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        # Write a duplicate step row by hand.
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        with open(os.path.join(case.coarse, "vpfp_step_diagnostics.dat"),
                  "a") as out:
            cols = STEP_COLS
            row = dict((c, 0.0) for c in cols)
            row["step"] = 101
            row["time_s"] = 1.001e-13
            row["accepted"] = 1
            row["split"] = 0
            row["domain_energy_before"] = 5.0e9
            row["domain_energy_after"] = 5.0e9
            row["energy_balance_residual"] = 2.0 * r
            out.write(" ".join(str(row[c]) for c in cols) + "\n")
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        return code == 3 and res.get("status") == "FAIL_AUDIT_STRUCTURE"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case5c(analyzer):  # non-finite
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {"final_poisson": [0.9 * 2.0 * r] * 10,
              "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {"energy_balance_residual": [float("nan")] * 10})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        return code == 4 and res.get("status") == "FAIL_NONFINITE_OR_CORRUPT_INPUT"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case5d(analyzer):  # split
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {"final_poisson": [0.9 * 2.0 * r] * 10,
              "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {"split": [0] * 10})
        # Force split on the first coarse row.
        _rewrite_single_col(case.coarse, "vpfp_step_diagnostics.dat", 0,
                            "split", "1")
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        return code == 3 and res.get("status") == "FAIL_AUDIT_STRUCTURE"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case5e(analyzer):  # rollback
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb = {"final_poisson": [0.9 * 2.0 * r] * 10,
              "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10,
                        {"collision_flux_rollback_count": [0] * 10})
        _rewrite_single_col(case.coarse, "vpfp_step_diagnostics.dat", 0,
                            "collision_flux_rollback_count", "1")
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        return code == 3 and res.get("status") == "FAIL_AUDIT_STRUCTURE"
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _rewrite_single_col(directory, filename, row_index, column, value):
    path = os.path.join(directory, filename)
    with open(path, "r") as stream:
        lines = stream.read().rstrip("\n").split("\n")
    header = lines[0].split()
    col = header.index(column)
    fields = lines[1 + row_index].split()
    fields[col] = value
    lines[1 + row_index] = " ".join(fields)
    with open(path, "w") as stream:
        stream.write("\n".join(lines) + "\n")


def case6(analyzer):  # roundoff-level stage -> order=not_evaluated_roundoff
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        # Dominant final_poisson 90%, x_half1 10%, and a tiny roundoff stage
        # collision_half2 ~ 1e-8 (power below the machine floor).
        sb_c = {
            "final_poisson": [0.9 * 2.0 * r] * 10,
            "x_half1": [0.1 * 2.0 * r] * 10,
            "collision_half2": [1.0e-8] * 10,
        }
        sb_f = {
            "final_poisson": [0.9 * r] * 20,
            "x_half1": [0.1 * r] * 20,
            "collision_half2": [1.0e-8] * 20,
        }
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb_c, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        ok = (code == 0 and
              res.get("stage_collision_half2_observed_order") ==
              "not_evaluated_roundoff")
        return ok
    finally:
        shutil.rmtree(root, ignore_errors=True)


def case7(analyzer):  # 2r/r per-step, same power -> observed order near 0
    root = mkdtemp("gate_h_")
    try:
        r = 1.0e4
        case = Case(root)
        times_c = make_times(10, 1.01e-13)
        times_f = make_times(20, 1.01e-13)
        coarse_residual = [2.0 * r] * 10
        fine_residual = [r] * 20
        sb_c = {"final_poisson": [0.9 * 2.0 * r] * 10,
                "x_half1": [0.1 * 2.0 * r] * 10}
        sb_f = {"final_poisson": [0.9 * r] * 20,
                "x_half1": [0.1 * r] * 20}
        write_step_file(case.coarse, 10, 101, times_c, coarse_residual,
                        [0.0] * 10, {})
        write_step_file(case.fine, 20, 101, times_f, fine_residual,
                        [0.0] * 20, {})
        write_stage_file(case.coarse, 10, 101, times_c, sb_c, coarse_residual)
        write_stage_file(case.fine, 20, 101, times_f, sb_f, fine_residual)
        result_path = os.path.join(root, "out.result")
        code = run_analyzer(analyzer, case.checkpoint, case.coarse, case.fine,
                            10, 20, result_path)
        res = read_result(result_path)
        order_text = res.get("full_residual_observed_order", "")
        try:
            order = float(order_text)
        except ValueError:
            return False
        return code == 0 and abs(order) < 0.05
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
