// Phase-2A acceptance test for the open-boundary remap sampling interface
// (sections 13.5, 14 phase 2A and 16.3.2).  Exercises
// OpenBackgroundBoundary::incoming_cell_average/is_incoming plus a light
// integration check through the production ConservativePpmRemap.
//
// Usage:
//   open_boundary_remap_test --case reservoir-and-absorbing [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "conservative_ppm_remap.h"
#include "grid.h"
#include "open_boundary.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "reservoir-and-absorbing";
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

Species make_species(const SpatialGrid& grid)
{
    Species sp;
    sp.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
            Const::me, Param::dens, Param::temperature_e, false, grid);
    return sp;
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

OpenBackgroundBoundaryConfig absorbing_config()
{
    OpenBackgroundBoundaryConfig config;
    config.left_type = BackgroundXBoundaryType::ABSORBING;
    config.right_type = BackgroundXBoundaryType::ABSORBING;
    config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    return config;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        // This test drives the boundary/remap as a single full-domain process
        // (rank=0, size=1); MPI_COMM_WORLD reductions would otherwise sum
        // identical per-rank states.
        std::cerr << "open_boundary_remap_test must run with exactly 1 rank; "
                     "use plain ./build/open_boundary_remap_test (no "
                     "yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }
    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: open_boundary_remap_test --case "
                     "reservoir-and-absorbing [--result <path>]\n";
    }

    double reservoir_sum_error = 0.0;
    double absorbing_max_lambda = 0.0;
    double ghost_consistency_error = 0.0;
    double steady_state_diff_rel = 0.0;
    double absorbing_ledger_error = 0.0;
    double absorbing_final_mass = 0.0;

    if (ok) {
        SpatialGrid grid;
        grid.init_with_domain(0, 1, 100, 40.0 * Const::micro);
        Species sp = make_species(grid);
        const int ng = grid.nghost;
        const int nxl = grid.nx_local;

        // 1. Reservoir line densities integrate to the reservoir density.
        OpenBackgroundBoundary reservoir(reservoir_config());
        double sum_lambda = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                sum_lambda += reservoir.incoming_cell_average(
                    PhysicalSide::LEFT, iv, imu, 0.0, sp);
            }
        }
        reservoir_sum_error =
            std::fabs(sum_lambda - Param::dens) / Param::dens;

        // 2. is_incoming direction logic.
        const double vx_pos = 0.1 * Const::c;
        const double vx_neg = -0.1 * Const::c;
        ok = ok && reservoir.is_incoming(PhysicalSide::LEFT, vx_pos) &&
             !reservoir.is_incoming(PhysicalSide::LEFT, vx_neg) &&
             !reservoir.is_incoming(PhysicalSide::LEFT, 0.0) &&
             reservoir.is_incoming(PhysicalSide::RIGHT, vx_neg) &&
             !reservoir.is_incoming(PhysicalSide::RIGHT, vx_pos) &&
             !reservoir.is_incoming(PhysicalSide::RIGHT, 0.0);

        // 3. Absorbing inflow is zero.
        OpenBackgroundBoundary absorbing(absorbing_config());
        double max_lambda = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                max_lambda = std::max(
                    max_lambda,
                    std::fabs(absorbing.incoming_cell_average(
                        PhysicalSide::LEFT, iv, imu, 0.0, sp)));
                max_lambda = std::max(
                    max_lambda,
                    std::fabs(absorbing.incoming_cell_average(
                        PhysicalSide::RIGHT, iv, imu, 0.0, sp)));
            }
        }
        absorbing_max_lambda = max_lambda;
        ok = ok && max_lambda == 0.0;

        // 4. fill_ghosts reservoir values agree with incoming_cell_average.
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double lambda = reservoir.incoming_cell_average(
                    PhysicalSide::LEFT, iv, imu, 0.0, sp);
                for (int i = 0; i < nxl; ++i) {
                    sp.f[idx3(ng + i, iv, imu)] = lambda * grid.dx;
                }
            }
        }
        reservoir.fill_ghosts(sp, grid, 0, 1);
        for (int g = 0; g < grid.nghost; ++g) {
            const int left_ghost = ng - 1 - g;
            const int right_ghost = ng + nxl + g;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double lambda = reservoir.incoming_cell_average(
                        PhysicalSide::LEFT, iv, imu, 0.0, sp);
                    if (sp.cgrid.vx[idx2(iv, imu)] > 0.0) {
                        ghost_consistency_error = std::max(
                            ghost_consistency_error,
                            std::fabs(sp.f[idx3(left_ghost, iv, imu)] -
                                      lambda * grid.dx));
                    }
                    if (sp.cgrid.vx[idx2(iv, imu)] < 0.0) {
                        ghost_consistency_error = std::max(
                            ghost_consistency_error,
                            std::fabs(sp.f[idx3(right_ghost, iv, imu)] -
                                      lambda * grid.dx));
                    }
                }
            }
        }
        ghost_consistency_error /=
            std::max(1.0, Param::dens * grid.dx);

        // 5. Light integration: reservoir steady state through the remap.
        Species out = make_species(grid);
        ConservativePpmRemap remap;
        remap.init(grid, sp.cgrid);
        const double dt = 4.0 * grid.dx / (0.3 * Const::c);
        remap.advect_x(sp, out, dt, 0.0, reservoir, 0, 1);
        double max_diff = 0.0;
        for (int i = 0; i < nxl; ++i) {
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    max_diff = std::max(
                        max_diff,
                        std::fabs(out.f[idx3(ng + i, iv, imu)] -
                                  sp.f[idx3(ng + i, iv, imu)]));
                }
            }
        }
        steady_state_diff_rel =
            max_diff / std::max(1.0, Param::dens * grid.dx);

        // 6. Light integration: absorbing outflow, no wrap.
        std::fill(sp.f.begin(), sp.f.end(), 0.0);
        // Use the fastest available slice for a quick exit.
        double vmax = 0.0;
        int best_j = 0, best_k = 0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                if (sp.cgrid.vx[idx2(iv, imu)] > vmax) {
                    vmax = sp.cgrid.vx[idx2(iv, imu)];
                    best_j = iv;
                    best_k = imu;
                }
            }
        }
        const double lambda0 = 1.0e20;
        for (int i = 0; i < nxl; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * grid.dx;
            if (x >= 2.0 * Const::micro && x <= 6.0 * Const::micro) {
                sp.f[idx3(ng + i, best_j, best_k)] = lambda0 * grid.dx;
            }
        }
        const double dt_absorb = 8.0 * grid.dx / vmax;
        double before = 0.0;
        for (int i = 0; i < nxl; ++i) before += sp.f[idx3(ng + i, best_j, best_k)];
        double max_ledger = 0.0;
        for (int step = 0; step < 30; ++step) {
            const RemapDiagnostics d =
                remap.advect_x(sp, sp, dt_absorb, 0.0, absorbing, 0, 1);
            double after = 0.0;
            for (int i = 0; i < nxl; ++i) after += sp.f[idx3(ng + i, best_j, best_k)];
            max_ledger = std::max(
                max_ledger,
                std::fabs((after - before) -
                          (d.inflow_number - d.outflow_number)) /
                    std::max(1.0, before));
            before = after;
        }
        absorbing_ledger_error = max_ledger;
        absorbing_final_mass = 0.0;
        for (int i = 0; i < nxl; ++i) {
            absorbing_final_mass += sp.f[idx3(ng + i, best_j, best_k)];
        }
        absorbing_final_mass /= std::max(1.0, lambda0 * 4.0 * Const::micro);
    }

    bool pass = ok &&
        reservoir_sum_error <= 1.0e-10 &&
        ghost_consistency_error <= 1.0e-10 &&
        steady_state_diff_rel <= 1.0e-10 &&
        absorbing_ledger_error <= 1.0e-10 &&
        absorbing_final_mass <= 1.0e-10;

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::trunc);
        if (out) {
            out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
            out << "reservoir_sum_error=" << reservoir_sum_error << "\n";
            out << "absorbing_max_lambda=" << absorbing_max_lambda << "\n";
            out << "ghost_consistency_error=" << ghost_consistency_error << "\n";
            out << "steady_state_diff_rel=" << steady_state_diff_rel << "\n";
            out << "absorbing_ledger_error=" << absorbing_ledger_error << "\n";
            out << "absorbing_final_mass_rel=" << absorbing_final_mass << "\n";
        } else {
            pass = false;
        }
    }
    std::cout << "reservoir_sum_error=" << reservoir_sum_error
              << " absorbing_max_lambda=" << absorbing_max_lambda
              << " ghost_consistency_error=" << ghost_consistency_error
              << " steady_state_diff_rel=" << steady_state_diff_rel
              << " absorbing_ledger_error=" << absorbing_ledger_error
              << " absorbing_final_mass_rel=" << absorbing_final_mass << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
