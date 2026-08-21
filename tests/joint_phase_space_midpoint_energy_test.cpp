#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "joint_phase_space_midpoint.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "species.h"
#include "vpfp_integrator.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result_path;
    std::string test_case = "smooth-background";
    bool parsed = true;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) result_path = argv[++i];
        else parsed = false;
    }
    if (!parsed || (test_case != "smooth-background" &&
                    test_case != "smooth-perturbed-background")) {
        if (rank == 0)
            std::cerr << "usage: joint_phase_space_midpoint_energy_test "
                         "--case smooth-background|smooth-perturbed-background "
                         "[--result path]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(rank, size, 4, 4.0e-6);
    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = BackgroundXBoundaryType::PERIODIC;
    boundary_config.right_type = BackgroundXBoundaryType::PERIODIC;
    boundary_config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    boundary_config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    OpenBackgroundBoundary background_boundary(boundary_config);
    ElectrostaticBoundary field_boundary = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid, field_boundary);
    ZeroCollisionCoefficients provider;
    CylindricalFokkerPlanckCollision collision(
        provider, CollisionIntegratorType::BACKWARD_EULER);
    VpfpIntegrator integrator(background_boundary, field_solver, collision);
    integrator.set_beam_enabled(false);
    integrator.set_background_phase_space_mode(
        BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY);
    integrator.init(grid);

    Species electrons;
    electrons.init("bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                   Const::me, Param::dens, Param::temperature_e, false, grid);
    std::vector<double> electron_profile(
        static_cast<size_t>(grid.nx_local), Param::dens);
    if (test_case == "smooth-perturbed-background") {
        for (int ix = 0; ix < grid.nx_local; ++ix) {
            const int ig = grid.ix_start + ix;
            electron_profile[static_cast<size_t>(ix)] =
                Param::dens * (1.0 + 1.0e-4 * std::cos(
                    2.0 * Const::pi *
                    (static_cast<double>(ig) + 0.5) /
                    static_cast<double>(grid.nx_global)));
        }
        electrons.initialize_maxwellian_profile(electron_profile, 0.0);
    } else {
        electrons.initialize_maxwellian(0.0);
    }
    electrons.compute_moments();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    // Match the fixed-ion density to the actual discrete Maxwellian moment;
    // otherwise the test starts with a roundoff-scale net charge but a finite
    // Poisson field, and the first-step energy gate measures that preparation
    // defect instead of the J1 joint flux.
    std::vector<double> ion_density =
        test_case == "smooth-background"
            ? electrons.number_density
            : std::vector<double>(static_cast<size_t>(grid.nx_local),
                                  Param::dens);
    std::vector<double> empty_tail(static_cast<size_t>(grid.nx_local), 0.0);
    std::vector<double> empty_beam(static_cast<size_t>(grid.nx_local), 0.0);
    fields.set_charge_density(electrons, empty_tail, empty_beam, ion_density);
    double electron_density_min_local = 0.0;
    double electron_density_max_local = 0.0;
    double ion_density_min_local = 0.0;
    double ion_density_max_local = 0.0;
    if (!electrons.number_density.empty()) {
        electron_density_min_local = electrons.number_density[0];
        electron_density_max_local = electrons.number_density[0];
    }
    if (!ion_density.empty()) {
        ion_density_min_local = ion_density[0];
        ion_density_max_local = ion_density[0];
    }
    for (size_t i = 1; i < electrons.number_density.size(); ++i) {
        electron_density_min_local = std::min(
            electron_density_min_local, electrons.number_density[i]);
        electron_density_max_local = std::max(
            electron_density_max_local, electrons.number_density[i]);
    }
    for (size_t i = 1; i < ion_density.size(); ++i) {
        ion_density_min_local = std::min(ion_density_min_local, ion_density[i]);
        ion_density_max_local = std::max(ion_density_max_local, ion_density[i]);
    }
    double electron_density_min = 0.0;
    double electron_density_max = 0.0;
    double ion_density_min = 0.0;
    double ion_density_max = 0.0;
    MPI_Allreduce(&electron_density_min_local, &electron_density_min, 1,
                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&electron_density_max_local, &electron_density_max, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ion_density_min_local, &ion_density_min, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&ion_density_max_local, &ion_density_max, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    double initial_rho_linf_local = 0.0;
    for (size_t i = 0; i < fields.rho.size(); ++i)
        initial_rho_linf_local = std::max(
            initial_rho_linf_local, std::fabs(fields.rho[i]));
    double initial_rho_linf = 0.0;
    MPI_Allreduce(&initial_rho_linf_local, &initial_rho_linf, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    field_solver.solve(fields, rank, size);
    double initial_E_linf_local = 0.0;
    for (size_t i = 0; i < fields.Ex_face.size(); ++i)
        initial_E_linf_local = std::max(
            initial_E_linf_local, std::fabs(fields.Ex_face[i]));
    double initial_E_linf = 0.0;
    MPI_Allreduce(&initial_E_linf_local, &initial_E_linf, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    const double initial_gauss_linf = field_solver.diagnostics().residual_linf;

    const double dt = 1.0e-18;
    const VpfpStepResult step = integrator.advance(
        electrons, beam, fields, ion_density, 0.0, dt, rank, size);
    int local_pass = step.accepted && step.finite && step.gauss_ok &&
        step.joint_midpoint_enabled && step.joint_midpoint_converged &&
        step.joint_midpoint_pairing_field_built &&
        !step.split_used && step.failure_code == 0;
    int global_pass = 0;
    MPI_Allreduce(&local_pass, &global_pass, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::ostream* out = &std::cout;
        std::ofstream file;
        if (!result_path.empty()) {
            file.open(result_path.c_str(), std::ios::trunc);
            if (!file) global_pass = 0;
            else out = &file;
        }
        *out << std::setprecision(17)
             << "status=" << (global_pass ? "PASS" : "FAIL") << "\n"
             << "case=" << test_case << "\n"
             << "background_phase_space_mode=joint-midpoint-energy\n"
             << "accepted=" << (step.accepted ? 1 : 0) << "\n"
             << "finite=" << (step.finite ? 1 : 0) << "\n"
             << "gauss_ok=" << (step.gauss_ok ? 1 : 0) << "\n"
             << "converged=" << (step.joint_midpoint_converged ? 1 : 0) << "\n"
             << "production_advance_called=1\n"
             << "joint_midpoint_pairing_field_built="
             << (step.joint_midpoint_pairing_field_built ? 1 : 0) << "\n"
             << "joint_midpoint_delta_k_x="
             << step.joint_midpoint_delta_k_x << "\n"
             << "joint_midpoint_delta_k_u="
             << step.joint_midpoint_delta_k_u << "\n"
             << "joint_midpoint_u_face_work="
             << step.joint_midpoint_u_face_work << "\n"
             << "joint_midpoint_force_current_work="
             << step.joint_midpoint_force_current_work << "\n"
             << "joint_midpoint_charge_current_work="
             << step.joint_midpoint_charge_current_work << "\n"
             << "joint_midpoint_charge_current_work_interior="
             << step.joint_midpoint_charge_current_work_interior << "\n"
             << "joint_midpoint_charge_current_work_endpoint="
             << step.joint_midpoint_charge_current_work_endpoint << "\n"
             << "joint_midpoint_poisson_potential_charge_work="
             << step.joint_midpoint_poisson_potential_charge_work << "\n"
             << "joint_midpoint_poisson_transport_residual="
             << step.joint_midpoint_poisson_transport_residual << "\n"
             << "joint_midpoint_current_pair_residual="
             << step.joint_midpoint_current_pair_residual << "\n"
             << "joint_midpoint_force_charge_residual="
             << step.joint_midpoint_force_charge_residual << "\n"
             << "joint_midpoint_domain_energy_change="
             << step.joint_midpoint_domain_energy_change << "\n"
             << "delta_k_x=" << step.joint_midpoint_delta_k_x << "\n"
             << "delta_k_u=" << step.joint_midpoint_delta_k_u << "\n"
             << "u_face_work=" << step.joint_midpoint_u_face_work << "\n"
             << "charge_current_work="
             << step.joint_midpoint_charge_current_work << "\n"
             << "field_energy_change="
             << step.joint_midpoint_field_energy_change << "\n"
             << "electrode_work=" << step.joint_midpoint_electrode_work << "\n"
             << "poisson_potential_charge_work="
             << step.joint_midpoint_poisson_potential_charge_work << "\n"
             << "poisson_transport_residual="
             << step.joint_midpoint_poisson_transport_residual << "\n"
             << "current_pair_residual="
             << step.joint_midpoint_current_pair_residual << "\n"
             << "domain_energy_change="
             << step.joint_midpoint_domain_energy_change << "\n"
             << "min_mass=" << step.joint_midpoint_min_mass << "\n"
             << "max_mass=" << step.joint_midpoint_max_mass << "\n"
             << "iterations=" << step.joint_midpoint_iterations << "\n"
             << "residual_linf=" << step.joint_midpoint_residual_linf << "\n"
             << "poisson_residual_linf="
             << step.joint_midpoint_poisson_residual_linf << "\n"
             << "energy_residual=" << step.joint_midpoint_energy_residual << "\n"
             << "failure_stage=" << step.failure_stage << "\n"
             << "electron_density_max=" << electron_density_max << "\n"
             << "electron_density_min=" << electron_density_min << "\n"
             << "ion_density_max=" << ion_density_max << "\n"
             << "ion_density_min=" << ion_density_min << "\n"
             << "initial_rho_linf=" << initial_rho_linf << "\n"
             << "initial_E_linf=" << initial_E_linf << "\n"
             << "initial_gauss_linf=" << initial_gauss_linf << "\n"
             << "failure_code=" << step.failure_code << "\n"
             << "iteration_log_count="
             << step.joint_midpoint_iterations_log.size() << "\n";
        for (size_t i = 0; i < step.joint_midpoint_iterations_log.size(); ++i) {
            const JointPhaseSpaceIterationRecord& record =
                step.joint_midpoint_iterations_log[i];
            *out << "iteration_" << i << "_number=" << record.iteration << "\n"
                 << "iteration_" << i << "_gmres_dimension="
                 << record.gmres_dimension << "\n"
                 << "iteration_" << i << "_residual_linf="
                 << record.residual_linf << "\n"
                 << "iteration_" << i << "_poisson_linf="
                 << record.phi_residual_linf << "\n"
                 << "iteration_" << i << "_line_search_alpha="
                 << record.line_search_alpha << "\n"
                 << "iteration_" << i << "_trial_min_mass="
                 << record.trial_min_mass << "\n"
                 << "iteration_" << i << "_accepted=" << record.accepted << "\n"
                 << "iteration_" << i << "_failure_code="
                 << record.failure_code << "\n";
        }
    }
    MPI_Finalize();
    return global_pass ? 0 : 1;
}
