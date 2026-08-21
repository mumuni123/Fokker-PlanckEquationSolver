// Phase-4 acceptance test for Beam MPI migration, open-boundary outflow and
// Poisson coupling (sections 7.2/13.9 and 16.6.3).  Drives the production
// VpfpIntegrator with the Beam enabled; designed for yhrun -n 5.
//
// Usage:
//   beam_open_boundary_mpi_test --case migration-left-right-poisson
//       [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "beam_pic.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "migration-left-right-poisson";
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
    return !args.test_case.empty();
}

OpenBackgroundBoundaryConfig reservoir_config()
{
    OpenBackgroundBoundaryConfig config;
    config.left_type = BackgroundXBoundaryType::RESERVOIR;
    config.right_type = BackgroundXBoundaryType::RESERVOIR;
    config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    return config;
}

struct Sim {
    SpatialGrid grid;
    Species electrons;
    BeamPIC beam;
    EMFields fields;
    OpenBackgroundBoundary boundary;
    OpenElectrostaticSolver field_solver;
    ZeroCollisionCoefficients zero_provider;
    CylindricalFokkerPlanckCollision collision;
    VpfpIntegrator integrator;
    std::vector<double> ion_density;
    double dt;
    int rank;
    int size;

    Sim(int nx, const OpenBackgroundBoundaryConfig& bcfg,
        const std::vector<double>& ions, int rank_, int size_)
        : boundary(bcfg),
          collision(zero_provider, CollisionIntegratorType::BACKWARD_EULER),
          integrator(boundary, field_solver, collision),
          ion_density(ions), rank(rank_), size(size_)
    {
        grid.init_with_domain(rank, size, nx, 2.0 * Const::micro);
        electrons.init("background_electrons",
                       SpeciesType::BACKGROUND_ELECTRON,
                       -Const::qe, Const::me, Param::dens,
                       Param::temperature_e, false, grid);
        beam.init(grid);
        fields.init(grid);
        field_solver.init(
            grid, { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 });
        integrator.init(grid);
        integrator.set_beam_enabled(true);
        dt = Param::dt_multiplier / Param::omega_pe;
    }

    void initial_field_solve()
    {
        electrons.compute_moments();
        beam.deposit_density(grid, rank, size);
        const std::vector<double> empty_tail_density;
        fields.set_charge_density(electrons, empty_tail_density,
                                  beam.density, ion_density);
        field_solver.solve(fields, rank, size);
    }
};

struct CaseMetrics {
    bool all_accepted;
    double max_balance_error;
    double max_gauss_residual;
    double fast_particles_remaining;
    double fast_outflow_weight;
    double beam_density_peak;
    int rejected_step;
    int failure_code;
    CaseMetrics()
        : all_accepted(true), max_balance_error(0.0), max_gauss_residual(0.0),
          fast_particles_remaining(0.0), fast_outflow_weight(0.0),
          beam_density_peak(0.0), rejected_step(-1), failure_code(0)
    {}
};

CaseMetrics run_case(int rank, int size)
{
    CaseMetrics m;
    if (200 / std::max(1, size) < 5) {
        // The background remap needs a ~3-cell halo; with a 200-cell domain
        // every rank must keep at least 5 cells (size <= 40).  Beyond that
        // the halo exchange reads overlapping/garbage data and corrupts the
        // state (and hence the beam field).
        if (rank == 0) {
            std::cerr << "beam_open_boundary_mpi_test: too many ranks for the "
                         "200-cell domain (nx_local < 5); use at most 40 ranks "
                         "(e.g. -n 5).\n";
        }
        m.all_accepted = false;
        return m;
    }
    // Short domain so the relativistic tracers cross every rank boundary and
    // still exit the open right boundary within the run: at 5 ranks the
    // rank width is 0.4 um and the u=4 tracer travels ~0.0074 um/step.
    const int nx = 200;
    std::vector<double> ions(static_cast<size_t>(nx / std::max(1, size)),
                             Param::dens);
    Sim sim(nx, reservoir_config(), ions, rank, size);

    // Uniform Maxwellian background.
    {
        const double sigma_u =
            std::sqrt(Param::temperature_e /
                      (Const::me * Const::c * Const::c));
        const int ng = sim.grid.nghost;
        const int nxl = sim.grid.nx_local;
        const CylindricalVelocityGrid& cg = sim.electrons.cgrid;
        double raw = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double u = cg.upar_cells[static_cast<size_t>(j)];
                const double up = cg.uperp_cells[static_cast<size_t>(k)];
                raw += std::exp(-(u * u + up * up) /
                                (2.0 * sigma_u * sigma_u)) *
                       cg.cell_phase_volume(j, k);
            }
        }
        const double norm = Param::dens / raw;
        for (int ix = 0; ix < nxl; ++ix) {
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double u = cg.upar_cells[static_cast<size_t>(j)];
                    const double up = cg.uperp_cells[static_cast<size_t>(k)];
                    const double f3 = norm * std::exp(
                        -(u * u + up * up) / (2.0 * sigma_u * sigma_u));
                    sim.electrons.f[idx3(ng + ix, j, k)] =
                        f3 * sim.grid.dx * cg.cell_phase_volume(j, k);
                }
            }
        }
        sim.boundary.fill_ghosts(sim.electrons, sim.grid, rank, size);
        sim.electrons.compute_moments();
    }

    // Two relativistic tracer particles with a momentum distinct from the
    // injected beam (u = 8.6): the right one starts in rank 0 just left of
    // the first rank boundary and migrates through every rank before exiting
    // the right boundary; the left one exits the left boundary immediately.
    {
        // Distinct from the production beam momentum (8.6 me c) so the
        // tracer can be identified unambiguously.
        const double p_fast = 4.0 * Const::me * Const::c;
        const double w = beam_macro_weight(sim.grid);
        BeamParticle right = { 0.39 * Const::micro, p_fast, w };
        BeamParticle left = { 0.15 * Const::micro, -p_fast, w };
        sim.beam.particles.push_back(right);
        sim.beam.particles.push_back(left);
    }
    sim.initial_field_solve();

    const int steps = 250;
    double beam_before = 0.0;
    double cumulative_outflow = 0.0;
    for (int step = 0; step < steps; ++step) {
        const VpfpStepResult result =
            sim.integrator.advance(sim.electrons, sim.beam, sim.fields,
                                   sim.ion_density,
                                   static_cast<double>(step) * sim.dt,
                                   sim.dt, rank, size);
        if (!result.accepted) {
            m.all_accepted = false;
            m.rejected_step = step;
            m.failure_code = result.failure_code;
            std::cerr << "[reject] rank=" << rank
                      << " step=" << step
                      << " failure_code=" << result.failure_code
                      << " balance=" << result.ledger.remap_ledger_residual
                      << " tail="
                      << result.ledger.background_tail_number_loss
                      << " N_b=" << result.ledger.beam_number_after
                      << "\n";
            break;
        }
        if (step == 0) beam_before = result.ledger.beam_number_before;
        const double balance =
            std::fabs(result.ledger.beam_number_after -
                      (result.ledger.beam_number_before +
                       result.ledger.beam_injected -
                       result.ledger.beam_outflow)) /
            std::max(1.0, result.ledger.beam_number_before +
                          result.ledger.beam_injected);
        m.max_balance_error = std::max(m.max_balance_error, balance);
        m.max_gauss_residual = std::max(
            m.max_gauss_residual,
            std::fabs(result.ledger.gauss_charge_residual));
        // The step ledger now carries the global outflow (identical on every
        // rank), so accumulate on rank 0 only to avoid a 5x factor.
        if (rank == 0) {
            cumulative_outflow += result.ledger.beam_outflow;
        }
        double local_peak = 0.0;
        for (size_t ix = 0; ix < sim.beam.density.size(); ++ix) {
            local_peak = std::max(local_peak, sim.beam.density[ix]);
        }
        m.beam_density_peak = std::max(m.beam_density_peak, local_peak);
    }

    // The two fast tracers must have left the domain (their weight is gone
    // from the surviving particles) and must be counted in the outflow.
    double local_fast = 0.0;
    for (size_t i = 0; i < sim.beam.particles.size(); ++i) {
        if (std::fabs(sim.beam.particles[i].px) > 3.0 * Const::me * Const::c &&
            std::fabs(sim.beam.particles[i].px) < 5.0 * Const::me * Const::c) {
            local_fast += sim.beam.particles[i].weight;
        }
    }
    double global_fast = 0.0;
    MPI_Allreduce(&local_fast, &global_fast, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    m.fast_particles_remaining = global_fast;
    double global_outflow = 0.0;
    MPI_Allreduce(&cumulative_outflow, &global_outflow, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    m.fast_outflow_weight = global_outflow;
    double global_peak = 0.0;
    MPI_Allreduce(&m.beam_density_peak, &global_peak, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    m.beam_density_peak = global_peak;
    double global_residual = 0.0;
    MPI_Allreduce(&m.max_gauss_residual, &global_residual, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    m.max_gauss_residual = global_residual;
    (void)beam_before;
    return m;
}

bool write_result_file(const std::string& path, const CaseMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "all_accepted=" << (m.all_accepted ? 1 : 0) << "\n";
    out << "max_balance_error=" << m.max_balance_error << "\n";
    out << "max_gauss_residual=" << m.max_gauss_residual << "\n";
    out << "fast_particles_remaining=" << m.fast_particles_remaining << "\n";
    out << "fast_outflow_weight=" << m.fast_outflow_weight << "\n";
    out << "beam_density_peak=" << m.beam_density_peak << "\n";
    out << "rejected_step=" << m.rejected_step << "\n";
    out << "failure_code=" << m.failure_code << "\n";
    return out.good();
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
        std::cerr << "usage: beam_open_boundary_mpi_test --case "
                     "migration-left-right-poisson [--result <path>]\n";
    }

    CaseMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "migration-left-right-poisson") {
        m = run_case(rank, size);
        SpatialGrid weight_grid;
        weight_grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
        const double w = beam_macro_weight(weight_grid);
        pass = m.all_accepted &&
               m.max_balance_error <= 1.0e-9 &&
               m.fast_particles_remaining == 0.0 &&
               m.fast_outflow_weight >= 1.99 * w &&
               m.beam_density_peak > 0.0 &&
               m.max_gauss_residual <= 1.0e-10;
    } else {
        pass = false;
    }
    int pass_all = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &pass_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    pass = pass_all != 0;

    if (rank == 0) {
        if (!write_result_file(args.result_path, m, pass)) pass = false;
        std::cout << "all_accepted=" << (m.all_accepted ? 1 : 0)
                  << " max_balance_error=" << m.max_balance_error
                  << " max_gauss_residual=" << m.max_gauss_residual
                  << " fast_particles_remaining=" << m.fast_particles_remaining
                  << " fast_outflow_weight=" << m.fast_outflow_weight
                  << " beam_density_peak=" << m.beam_density_peak
                  << " rejected_step=" << m.rejected_step
                  << " failure_code=" << m.failure_code << "\n";
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
