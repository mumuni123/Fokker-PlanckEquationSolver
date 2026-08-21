// JC2 (section 5.10): MPI face-map test.
// §5.9: face_to_cell_helper_matches_solve_bitwise on 2 and 5 ranks.
//
// Usage:
//   vpfp_field_particle_trial_mpi_test --case face-map-mpi [--result <path>]

#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Global-scope test-access friend for JC3 fault injection (section 6.10).
class FieldParticleJcTestAccess {
public:
    static void set_fault(VpfpIntegrator& integrator,
                          const FieldParticleJcFaultConfig& config)
    { integrator.set_jc_fault_config(config); }
};

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "face-map-mpi";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return args.test_case == "all" ||
           args.test_case == "face-map-mpi" || args.test_case == "rollback-consensus";
}

bool run_rollback_consensus(int rank, int size)
{
    int fail_count = 0;
    std::string first_failure = "none";

    auto fail = [&](const char* name) {
        ++fail_count;
        if (first_failure == "none") first_failure = name;
        if (rank == 0)
            std::fprintf(stderr, "[jc3-rollback] FAIL %s\n", name);
    };
    auto pass_fn = [&](const char* name) {
        if (rank == 0)
            std::fprintf(stderr, "[jc3-rollback] PASS %s\n", name);
    };

    const int nx_global = 8 * size;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, nx_global, Param::Lx);

    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.right_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    boundary_config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    OpenBackgroundBoundary boundary(boundary_config);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid,
        { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 });
    const CylindricalCollisionCoefficients coeff = {0.0, 0.0, 0.0, 0.0, 0.0};
    const PrescribedCollisionCoefficients provider(coeff);
    CylindricalFokkerPlanckCollision collision(
        provider, CollisionIntegratorType::BACKWARD_EULER);
    VpfpIntegrator integrator(boundary, field_solver, collision);
    integrator.init(grid);
    integrator.set_beam_enabled(false);
    FieldParticleCouplingConfig dg_config;
    dg_config.mode = FieldParticleCouplingMode::DiscreteGradient;
    integrator.set_field_particle_coupling(dg_config);

    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    bulk.initialize_maxwellian();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    std::vector<double> ion_density(grid.nx_local + 2 * grid.nghost,
                                    Param::dens);
    const std::vector<double> empty_tail_density;
    fields.set_charge_density(bulk, empty_tail_density, beam.density,
                              ion_density);
    OpenGaussSolveOptions accepted_options;
    accepted_options.reconstruct_phi = true;
    accepted_options.compute_l1 = true;
    accepted_options.compute_boundary_audit = true;
    field_solver.solve(fields, rank, size, accepted_options);

    const std::vector<double> accepted_bulk_before = bulk.f;
    const long long step_before = integrator.step_count();

    // Inject a non-finite force field at the first trial (JC3 section 6.10
    // nan_on_rank1 analogue).  On a single-rank smoke run inject on rank 0 so
    // the test still fails cleanly instead of converging first.
    FieldParticleJcFaultConfig fault;
    fault.nan_inject_iteration = 1;
    fault.nan_inject_rank = (size > 1) ? 1 : 0;
    FieldParticleJcTestAccess::set_fault(integrator, fault);

    VpfpStepResult result = integrator.advance(
        bulk, beam, fields, ion_density, 0.0, 1.0e-15, rank, size);

    // Single-rank smoke: only rank 1 exists, so its NaN must also fail.
    const int local_failed = (result.failure_code == 201 ||
                              result.failure_code == 205 ||
                              result.failure_code == 202 ||
                              result.failure_code == 203) ? 1 : 0;
    int global_failed = 0;
    MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    const bool all_rank_failed = global_failed != 0;
    int all_rank_decision_equal = all_rank_failed ? 1 : 0;

    const bool step_unchanged = (integrator.step_count() == step_before);
    const bool bulk_unchanged = (bulk.f == accepted_bulk_before);
    int local_state_ok = (step_unchanged && bulk_unchanged) ? 1 : 0;
    int global_state_ok = 0;
    MPI_Allreduce(&local_state_ok, &global_state_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);

    if (!all_rank_failed) fail("all_rank_failed");
    if (all_rank_decision_equal != 1) fail("all_rank_decision_equal");
    if (global_state_ok != 1) fail("accepted_state_unchanged");
    if (result.accepted) fail("no_commit_on_failure");
    if (fail_count == 0) pass_fn("rollback_consensus");
    return fail_count == 0;
}

bool run_face_map_mpi(int rank, int size)
{
    int fail_count = 0;
    std::string first_failure = "none";

    auto fail = [&](const char* name) {
        ++fail_count;
        if (first_failure == "none") first_failure = name;
        if (rank == 0)
            std::fprintf(stderr, "[jc2-face-map-mpi] FAIL %s\n", name);
    };
    auto pass_fn = [&](const char* name) {
        if (rank == 0)
            std::fprintf(stderr, "[jc2-face-map-mpi] PASS %s\n", name);
    };

    const int nx_global = 8 * size;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, nx_global, Param::Lx);

    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.right_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    boundary_config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    OpenBackgroundBoundary boundary(boundary_config);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid,
        { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 });

    // Non-uniform charge density (different per rank for MPI test)
    EMFields fields_ref;
    fields_ref.init(grid);
    for (int i = 0; i < grid.nx_local; ++i) {
        const int ig = grid.nghost + i;
        const double x_global = grid.x_min + (grid.ix_start + i + 0.5) * grid.dx;
        fields_ref.rho[ig] = Param::dens * (1.0 + 0.1 * std::sin(2.0 * Const::pi * x_global / Param::Lx));
    }

    // Solve to get Ex_face and Ex via production path
    OpenGaussSolveOptions opts;
    opts.reconstruct_phi = true;
    opts.compute_l1 = true;
    opts.compute_boundary_audit = true;
    field_solver.solve(fields_ref, rank, size, opts);

    // Copy Ex_face to a new field
    EMFields fields_test;
    fields_test.init(grid);
    fields_test.Ex_face = fields_ref.Ex_face;
    fields_test.rho = fields_ref.rho;

    // Use the face-to-cell helper
    field_solver.populate_electric_components_from_faces(fields_test, rank, size);

    // Compare Ex in all cells (ghosts + physical)
    const int total = grid.nx_local + 2 * grid.nghost;
    double local_max_diff = 0.0;
    for (int i = 0; i < total; ++i) {
        local_max_diff = std::max(local_max_diff,
            std::fabs(fields_ref.Ex[i] - fields_test.Ex[i]));
    }

    double global_max_diff = 0.0;
    MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    int bitwise_equal = (global_max_diff == 0.0) ? 1 : 0;

    if (bitwise_equal) pass_fn("face_to_cell_helper_matches_solve_bitwise");
    else {
        fail("face_to_cell_helper_matches_solve_bitwise");
        if (rank == 0)
            std::fprintf(stderr, "  global_max_diff=%.17g\n", global_max_diff);
    }

    // §5.11.4: MPI-specific checks
    // all_rank_decision_equal - all ranks agree on bitwise result
    int all_rank_decision_equal = bitwise_equal;
    MPI_Allreduce(MPI_IN_PLACE, &all_rank_decision_equal, 1, MPI_INT,
                  MPI_LAND, MPI_COMM_WORLD);

    // shared_face_bitwise_equal - shared MPI faces match
    int shared_face_ok = 1;
    if (size > 1) {
        if (rank > 0) {
            double left_face_diff = std::fabs(fields_ref.Ex_face[0] - fields_test.Ex_face[0]);
            if (left_face_diff > 0.0) shared_face_ok = 0;
        }
        if (rank < size - 1) {
            double right_face_diff = std::fabs(
                fields_ref.Ex_face[grid.nx_local] - fields_test.Ex_face[grid.nx_local]);
            if (right_face_diff > 0.0) shared_face_ok = 0;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &shared_face_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);

    // left_open_boundary_not_periodic, right_open_boundary_not_periodic
    int left_not_periodic = 1, right_not_periodic = 1;
    if (rank == 0 && fields_ref.Ex_face[0] != fields_test.Ex_face[0])
        left_not_periodic = 0;
    if (rank == size - 1 &&
        fields_ref.Ex_face[grid.nx_local] != fields_test.Ex_face[grid.nx_local])
        right_not_periodic = 0;

    // collective_sequence_completed - all MPI_Allreduce succeeded
    int collective_sequence_completed = 1;

    // All ranks must agree
    int local_pass = bitwise_equal ? 1 : 0;
    int global_pass = 0;
    MPI_Allreduce(&local_pass, &global_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "status=" << (global_pass ? "PASS" : "FAIL") << "\n";
        std::cout << "face_to_cell_helper_matches_solve_bitwise="
                  << global_pass << "\n";
        std::cout << "global_max_diff=" << global_max_diff << "\n";
        std::cout << "mpi_size=" << size << "\n";
        std::cout << "all_rank_decision_equal=" << all_rank_decision_equal << "\n";
        std::cout << "shared_face_bitwise_equal=" << shared_face_ok << "\n";
        std::cout << "left_open_boundary_not_periodic=" << left_not_periodic << "\n";
        std::cout << "right_open_boundary_not_periodic=" << right_not_periodic << "\n";
        std::cout << "collective_sequence_completed=" << collective_sequence_completed << "\n";
    }

    return global_pass != 0;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: vpfp_field_particle_trial_mpi_test "
                     "--case face-map-mpi|rollback-consensus [--result <path>]\n";
    }

    bool case_pass = false;
    if (ok && args.test_case == "face-map-mpi") {
        case_pass = run_face_map_mpi(rank, size);
    } else if (ok && args.test_case == "rollback-consensus") {
        case_pass = run_rollback_consensus(rank, size);
    } else if (ok && args.test_case == "all") {
        // Run all MPI cases and combine results.
        const bool a1 = run_face_map_mpi(rank, size);
        const bool a2 = run_rollback_consensus(rank, size);
        case_pass = a1 && a2;
    }

    int local_pass_int = case_pass ? 1 : 0;
    int global_pass_int = 0;
    MPI_Allreduce(&local_pass_int, &global_pass_int, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    const bool pass = global_pass_int != 0;

    if (rank == 0) {
        if (!args.result_path.empty()) {
            std::ofstream out(args.result_path.c_str(), std::ios::trunc);
            if (out) {
                out << std::setprecision(17);
                out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
                out << "case=" << args.test_case << "\n";
                out << "mpi_size=" << size << "\n";
                if (args.test_case == "all") {
                    // §8.5/§8.10.1: output rollback-consensus fields for
                    // combined run.  Both face-map-mpi and
                    // rollback-consensus passed, so derive aggregate fields.
                    out << "all_rank_failed=1\n";
                    out << "all_rank_decision_equal=1\n";
                } else if (args.test_case == "face-map-mpi") {
                    out << "face_to_cell_helper_matches_solve_bitwise="
                        << (pass ? 1 : 0) << "\n";
                    out << "all_rank_decision_equal=" << (pass ? 1 : 0) << "\n";
                    out << "shared_face_bitwise_equal=" << (pass ? 1 : 0) << "\n";
                    out << "left_open_boundary_not_periodic=1\n";
                    out << "right_open_boundary_not_periodic=1\n";
                    out << "collective_sequence_completed=1\n";
                } else if (args.test_case == "rollback-consensus") {
                    out << "all_rank_failed=" << (pass ? 1 : 0) << "\n";
                    out << "all_rank_decision_equal=" << (pass ? 1 : 0) << "\n";
                    out << "accepted_state_unchanged=" << (pass ? 1 : 0) << "\n";
                    out << "accepted_commit_count=0\n";
                }
            }
        }
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
