// Phase-2A MPI consistency test for the conservative x remap (sections 13.6,
// 14 phase 2A and 16.3.3).  A profile crosses rank interfaces and departs over
// multiple cells; the multi-rank result must match an independent single-rank
// reference solve of the same problem.
//
// Usage (1 or 5 ranks):
//   conservative_x_remap_mpi_test --case crossing-and-multicell [--result <path>]
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
    args.test_case = "crossing-and-multicell";
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

size_t find_slice_near(const Species& sp, double target_vx)
{
    size_t best = 0;
    double best_error = std::numeric_limits<double>::infinity();
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const size_t q = idx2(iv, imu);
            const double error = std::fabs(sp.cgrid.vx[q] - target_vx);
            if (error < best_error) {
                best_error = error;
                best = q;
            }
        }
    }
    return best;
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

void fill_gaussian(Species& sp, const SpatialGrid& grid, size_t q,
                   double center, double sigma, double amplitude)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int j = static_cast<int>(q) / Param::Nmu;
    const int k = static_cast<int>(q) % Param::Nmu;
    std::fill(sp.f.begin(), sp.f.end(), 0.0);
    for (int i = 0; i < nxl; ++i) {
        const double x =
            (static_cast<double>(grid.ix_start + i) + 0.5) * grid.dx;
        sp.f[idx3(ng + i, j, k)] =
            amplitude * std::exp(-(x - center) * (x - center) /
                                 (2.0 * sigma * sigma)) * grid.dx;
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: conservative_x_remap_mpi_test --case "
                     "crossing-and-multicell [--result <path>]\n";
    }

    double max_field_diff_rel = 0.0;
    double mass_rel_error = 0.0;
    double ledger_rel_error = 0.0;
    double min_cell_mass = 0.0;
    double mass_scale = 1.0;

    if (ok) {
        const int nx = 500;
        const double length = 40.0 * Const::micro;
        SpatialGrid grid;
        grid.init_with_domain(rank, size, nx, length);
        Species sp = make_species(grid);
        Species out = make_species(grid);
        OpenBackgroundBoundary boundary(absorbing_config());
        const size_t q = find_slice_near(sp, 0.4 * Const::c);
        const int j = static_cast<int>(q) / Param::Nmu;
        const int k = static_cast<int>(q) % Param::Nmu;
        const double v = sp.cgrid.vx[q];
        const double sigma = 2.0 * Const::micro;
        const double center = 14.0 * Const::micro;
        const double amplitude = 1.0e21;
        const double distance = 0.8 * Const::micro;  // 10-cell departure
        fill_gaussian(sp, grid, q, center, sigma, amplitude);

        ConservativePpmRemap remap;
        remap.init(grid, sp.cgrid);
        const double dt = distance / v;
        const RemapDiagnostics diag =
            remap.advect_x(sp, out, dt, 0.0, boundary, rank, size);

        // Independent single-rank reference solve of the same problem.
        SpatialGrid ref_grid;
        ref_grid.init_with_domain(0, 1, nx, length);
        Species ref_sp = make_species(ref_grid);
        Species ref_out = make_species(ref_grid);
        fill_gaussian(ref_sp, ref_grid, q, center, sigma, amplitude);
        ConservativePpmRemap ref_remap;
        ref_remap.init(ref_grid, ref_sp.cgrid);
        ref_remap.advect_x(ref_sp, ref_out, dt, 0.0, boundary, 0, 1);

        // Compare the local slice against the reference global cells.
        const int ng = grid.nghost;
        const int nxl = grid.nx_local;
        double max_diff = 0.0;
        double scale = 0.0;
        double worst_diff = 0.0;
        int worst_global_cell = -1;
        double worst_value = 0.0;
        double worst_ref = 0.0;
        for (int i = 0; i < nxl; ++i) {
            const int ref_ix = ref_grid.nghost + grid.ix_start + i;
            const double a = out.f[idx3(ng + i, j, k)];
            const double b = ref_out.f[idx3(ref_ix, j, k)];
            const double d = std::fabs(a - b);
            if (d > worst_diff) {
                worst_diff = d;
                worst_global_cell = grid.ix_start + i;
                worst_value = a;
                worst_ref = b;
            }
            max_diff = std::max(max_diff, d);
            scale = std::max(scale, std::fabs(b));
        }
        // Relative field error is normalized by the global peak of the
        // reference (the physically meaningful scale).  The old
        // max(1.0, local_scale) normalization degenerates to an absolute
        // threshold on ranks whose slice is dominated by the profile tail.
        max_field_diff_rel = max_diff;

        // Failure localization: every rank reports its own worst cell.
        std::cout << "[mpi-debug r=" << rank << "] global_cell="
                  << worst_global_cell
                  << " local=" << worst_value
                  << " ref=" << worst_ref
                  << " diff=" << worst_diff << "\n";
        std::cout.flush();
        MPI_Barrier(MPI_COMM_WORLD);

        // Global mass and ledger checks (diag is already global).
        mass_rel_error =
            std::fabs(diag.number_after - diag.number_before) /
            std::max(1.0, diag.number_before);
        ledger_rel_error =
            std::fabs((diag.number_after - diag.number_before) -
                      (diag.inflow_number - diag.outflow_number)) /
            std::max(1.0, diag.number_before);
        min_cell_mass = diag.minimum_cell_mass;
        mass_scale = std::max(1.0, scale);
    }

    double global_diff = 0.0, global_mass = 0.0, global_ledger = 0.0;
    double global_min_mass = 0.0, global_mass_scale = 1.0;
    MPI_Allreduce(&max_field_diff_rel, &global_diff, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&mass_rel_error, &global_mass, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&ledger_rel_error, &global_ledger, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&min_cell_mass, &global_min_mass, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&mass_scale, &global_mass_scale, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    global_diff /= std::max(1.0, global_mass_scale);

    const double roundoff = 1.0e-10;
    bool pass = ok &&
        global_diff <= roundoff &&
        global_mass <= 1.0e-12 &&
        global_ledger <= roundoff &&
        global_min_mass >= -128.0 * std::numeric_limits<double>::epsilon() *
                               global_mass_scale;

    if (rank == 0) {
        if (!args.result_path.empty()) {
            std::ofstream out(args.result_path.c_str(), std::ios::trunc);
            if (out) {
                out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
                out << "max_field_diff_rel=" << global_diff << "\n";
                out << "mass_rel_error=" << global_mass << "\n";
                out << "ledger_rel_error=" << global_ledger << "\n";
                out << "min_cell_mass=" << global_min_mass << "\n";
            } else {
                pass = false;
            }
        }
        std::cout << "max_field_diff_rel=" << global_diff
                  << " mass_rel_error=" << global_mass
                  << " ledger_rel_error=" << global_ledger
                  << " min_cell_mass=" << global_min_mass << "\n";
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
