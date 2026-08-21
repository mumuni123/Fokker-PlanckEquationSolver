// Phase-1 acceptance test for the runtime Beam macro-particle weight
// (sections 13.1, 8.2 and 16.2.2).  The injected total weight over a fixed
// physical time interval must be determined by the beam flux alone; changing
// dx must not change the physical injected charge.
//
// Usage:
//   beam_weight_grid_test --dx-um 0.005 --particles-per-cell 1000
//                         [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "beam_pic.h"
#include "grid.h"
#include "parameters.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct TestArgs {
    double dx_um;
    int particles_per_cell;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.dx_um = 0.005;
    args.particles_per_cell = 1000;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--dx-um") {
            if (i + 1 >= argc) return false;
            args.dx_um = std::strtod(argv[++i], NULL);
        } else if (arg == "--particles-per-cell") {
            if (i + 1 >= argc) return false;
            args.particles_per_cell = std::atoi(argv[++i]);
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return args.dx_um > 0.0 && args.particles_per_cell > 0;
}

double injected_number_over_interval(const SpatialGrid& grid, double interval_dt)
{
    BeamPIC beam;
    const BeamInjectionSchedule schedule =
        beam.generate_injection_schedule(grid, 0.0, interval_dt, 0);
    double total = 0.0;
    for (size_t i = 0; i < schedule.events.size(); ++i) {
        total += schedule.events[i].weight;
    }
    return total;
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
        std::cerr << "usage: beam_weight_grid_test --dx-um 0.005 "
                     "--particles-per-cell 1000 [--result <path>]\n";
    }
    if (ok && args.particles_per_cell != Param::beam_macro_particles_per_cell) {
        ok = false;
        if (rank == 0) {
            std::cerr << "--particles-per-cell must match the production "
                         "Param::beam_macro_particles_per_cell ("
                      << Param::beam_macro_particles_per_cell << ")\n";
        }
    }

    double new_weight = 0.0, old_weight = 0.0;
    double new_injected = 0.0, old_injected = 0.0;
    double physical_injected = 0.0;
    double new_dx = 0.0, old_dx = 0.0;

    if (ok) {
        const double length = 40.0 * Const::micro;
        const double new_dx_expected = args.dx_um * Const::micro;
        const double old_dx_expected = 0.002 * Const::micro;

        // Production grid (nx from Param) and the legacy 0.002 um grid, both
        // over the same 40 um open domain.
        SpatialGrid grid_new;
        grid_new.init_with_domain(0, 1, Param::nx, length);
        SpatialGrid grid_old;
        grid_old.init_with_domain(0, 1, 20000, length);
        new_dx = grid_new.dx;
        old_dx = grid_old.dx;

        const double tolerance = 64.0 *
            std::numeric_limits<double>::epsilon();
        if (std::fabs(new_dx - new_dx_expected) >
                tolerance * std::max(1.0, new_dx_expected) ||
            std::fabs(old_dx - old_dx_expected) >
                tolerance * std::max(1.0, old_dx_expected)) {
            ok = false;
        }

        // Runtime macro weight is densb * grid.dx / particles_per_cell.
        new_weight = beam_macro_weight(grid_new);
        old_weight = beam_macro_weight(grid_old);
        const double new_weight_expected =
            Param::densb * new_dx /
            static_cast<double>(Param::beam_macro_particles_per_cell);
        const double old_weight_expected =
            Param::densb * old_dx /
            static_cast<double>(Param::beam_macro_particles_per_cell);
        if (std::fabs(new_weight - new_weight_expected) >
                tolerance * std::max(1.0, new_weight_expected) ||
            std::fabs(old_weight - old_weight_expected) >
                tolerance * std::max(1.0, old_weight_expected)) {
            ok = false;
        }

        // Injected total weight over one physical interval: determined by the
        // beam flux densb * beam_v0 * T, independent of dx.  A 3 fs interval
        // keeps the macro-particle quantization error far below 1e-4.
        const double interval = 3.0 * Const::femto;
        new_injected = injected_number_over_interval(grid_new, interval);
        old_injected = injected_number_over_interval(grid_old, interval);
        physical_injected = Param::densb * Param::beam_v0 * interval;
        const double scale = std::max(1.0, physical_injected);
        const double new_rel = std::fabs(new_injected - physical_injected) / scale;
        const double old_rel = std::fabs(old_injected - physical_injected) / scale;
        const double grid_rel = std::fabs(new_injected - old_injected) / scale;
        if (new_rel > 1.0e-3 || old_rel > 1.0e-3 || grid_rel > 1.0e-4) {
            ok = false;
        }
    }

    int ok_all = ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &ok_all, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    ok = ok_all != 0;

    if (rank == 0) {
        if (!args.result_path.empty()) {
            std::ofstream out(args.result_path.c_str(), std::ios::trunc);
            if (out) {
                out << "status=" << (ok ? "PASS" : "FAIL") << "\n";
                out << "new_dx_um=" << (new_dx / Const::micro) << "\n";
                out << "old_dx_um=" << (old_dx / Const::micro) << "\n";
                out << "new_macro_weight=" << new_weight << "\n";
                out << "old_macro_weight=" << old_weight << "\n";
                out << "new_injected_number=" << new_injected << "\n";
                out << "old_injected_number=" << old_injected << "\n";
                out << "physical_injected_number=" << physical_injected << "\n";
                out.close();
            } else {
                ok = false;
            }
        }
        std::cout << "new_dx_um=" << (new_dx / Const::micro)
                  << " old_dx_um=" << (old_dx / Const::micro)
                  << " new_injected=" << new_injected
                  << " old_injected=" << old_injected
                  << " physical=" << physical_injected << "\n";
        std::cout << "status=" << (ok ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return ok ? 0 : 1;
}
