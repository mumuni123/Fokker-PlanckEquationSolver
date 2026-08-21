#!/usr/bin/env python3
"""Regression fixture for the offline H9 production validation tools."""

from __future__ import print_function

import os
import shutil
import subprocess
import sys
import tempfile


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VALIDATE = os.path.join(ROOT, "tools", "validate_h9_production.py")
COMPARE = os.path.join(ROOT, "tools", "compare_h9_restart_production.py")


STEP_HEADER = (
    "step time_s accepted gauss_linf gauss_charge_residual N_e_after N_b_after U_E K_e K_b "
    "wall_s tail_particle_count tail_particles_local_max P_bkg P_tail "
    "P_combined K_tail K_combined N_tail_before N_tail_after "
    "N_combined_after tail_number_balance_error collision_reservoir "
    "conversion_N conversion_Px conversion_K collision_flux_export_N "
    "collision_flux_export_K\n")
FLUX_HEADER = (
    "accepted_step time_s duplicate_count duplicate_id_count "
    "face_ledger_mismatch_count static_extractor_call_count "
    "conversion_number_residual conversion_px_residual "
    "conversion_energy_residual conversion_jx_residual "
    "conversion_pixx_residual conversion_piperp_residual\n")


def write_run(root):
    os.makedirs(root)
    with open(os.path.join(root, "vpfp_step_diagnostics.dat"), "w") as out:
        out.write(STEP_HEADER)
        # Deliberately retain a large dimensional gauss_linf: the validator
        # must gate the normalized charge residual instead.
        out.write("1 1e-17 1 491520 1e-15 1 0 2 3 4 0.1 8 5 0 0 0 0 0 0 1 1 0 0 1 0 1 0 0\n")
        out.write("2 2e-17 1 983040 -2e-15 1 0 2 3 4 0.2 8 5 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0\n")
    with open(os.path.join(root, "bulk_tail_flux_accepted_steps.dat"), "w") as out:
        out.write(FLUX_HEADER)
        out.write("1 1e-17 0 0 0 0 0 0 0 0 0 0\n")
        out.write("2 2e-17 0 0 0 0 0 0 0 0 0 0\n")
    checkpoint = os.path.join(root, "checkpoint_target2fs_t2fs_step2")
    os.makedirs(checkpoint)
    with open(os.path.join(checkpoint, "manifest.txt"), "w") as out:
        out.write("tail_conversion_mode flux-interface\n")
        out.write("collision_induced_conversion 1\n")
        out.write("tail_collision_weight_algorithm sentoku-kemp-bounded-v1\n")
        out.write("population_control_interval 0\n")
    for rank in range(2):
        with open(os.path.join(checkpoint, "rank_%06d.bin" % rank), "wb") as out:
            out.write(b"fixture-rank-%d" % rank)
    snapshot = os.path.join(root, "snapshot_t2fs_step2")
    os.makedirs(snapshot)
    with open(os.path.join(snapshot, "tail_threshold_interface_rank0.dat"), "w") as out:
        out.write("fixture\n")
    return checkpoint


def require_pass(path):
    with open(path, "r") as stream:
        content = stream.read()
    if "status=PASS" not in content:
        raise RuntimeError("expected PASS:\n%s" % content)


def main():
    work = tempfile.mkdtemp(prefix="h9_production_validator_")
    try:
        direct = os.path.join(work, "direct")
        restart = os.path.join(work, "restart")
        direct_checkpoint = write_run(direct)
        restart_checkpoint = write_run(restart)
        direct_result = os.path.join(work, "direct.result")
        restart_result = os.path.join(work, "restart.result")
        compare_result = os.path.join(work, "compare.result")
        common = [sys.executable, VALIDATE, "--mode", "beam12",
                  "--min-accepted-steps", "2", "--require-tail", "yes",
                  "--require-threshold-snapshot"]
        subprocess.check_call(common + ["--run", direct, "--result", direct_result])
        subprocess.check_call(common + ["--run", restart, "--result", restart_result])
        subprocess.check_call([sys.executable, COMPARE, "--direct", direct,
                               "--restart", restart, "--result", compare_result,
                               "--min-common-steps", "2", "--require-tail-nonempty",
                               "--direct-checkpoint", direct_checkpoint,
                               "--restart-checkpoint", restart_checkpoint,
                               "--require-checkpoint-hash"])
        require_pass(direct_result)
        require_pass(restart_result)
        require_pass(compare_result)
        print("status=PASS")
        return 0
    finally:
        shutil.rmtree(work)


if __name__ == "__main__":
    sys.exit(main())
