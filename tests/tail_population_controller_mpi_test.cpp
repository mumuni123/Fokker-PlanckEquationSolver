// Stage H5 MPI consistency for TailPopulationController (sections 7.10 and
// 15 H5): the controller runs per rank on locally owned tail particles after
// the final drift (no MPI-in-flight particle is merged), so the acceptance
// checks that
//   * no MPI deadlock / failure at 1, 2 and 5 ranks;
//   * the global seven-moment residual stays within the compression
//     tolerance;
//   * the global macro-particle count decreases when bins exceed the cap;
//   * local IDs stay unique and carry a valid creation rank.
//
// Usage:
//   tail_population_controller_mpi_test [--result <path>]
// Run with yhrun -n 1 / -n 2 / -n 5 (different --result paths).
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include "tail_moment_constraint.h"
#include "tail_population_controller.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

void add_bin_particle(BackgroundTailPIC& tail, const SpatialGrid& grid,
                      int cell, int variant, double weight,
                      std::uint64_t id)
{
    BackgroundTailParticle p;
    const double xc = (static_cast<double>(cell) + 0.5) * grid.dx;
    // Offsets in [0.15, 0.35] dx keep every particle inside `cell`.
    p.x = xc + (0.15 + 0.05 * static_cast<double>(variant % 5)) * grid.dx;
    p.ux = 7.0 + 0.01 * static_cast<double>(variant % 3);
    p.uy = 0.3 + 0.1 * static_cast<double>(variant % 4);
    p.uz = 0.2 + 0.05 * static_cast<double>(variant % 5);
    p.weight = weight;
    p.id = id;
    tail.particles.push_back(p);
}

void tail_moments_sum(const BackgroundTailPIC& tail, TailMoment7& m)
{
    m = TailMoment7();
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        TailMoment7 pm;
        tail_particle_moments(p.weight, p.x, p.ux, p.uy, p.uz, pm);
        m.n += pm.n;
        m.px += pm.px;
        m.jx += pm.jx;
        m.ke += pm.ke;
        m.pixx += pm.pixx;
        m.piperp += pm.piperp;
        m.xw += pm.xw;
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    TestArgs args;
    const bool parsed = parse_args(argc, argv, args);
    int parsed_all = parsed ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &parsed_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!parsed_all) {
        if (rank == 0) {
            std::cerr << "usage: tail_population_controller_mpi_test "
                         "[--result <path>]\n";
        }
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(rank, mpi_size, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);

    BackgroundTailPIC tail;
    tail.init(grid);
    // Three bins of 20 particles in interior cells owned by this rank.
    std::uint64_t counter = 0;
    const int base_cell = grid.ix_start + 2;
    for (int g = 0; g < 3; ++g) {
        for (int i = 0; i < 20; ++i) {
            add_bin_particle(
                tail, grid, base_cell + g, i,
                1.0e14 + 5.0e12 * static_cast<double>(i % 7),
                (static_cast<std::uint64_t>(rank) << 32) | counter++);
        }
    }

    TailMoment7 ref;
    tail_moments_sum(tail, ref);
    const std::uint64_t count_before = tail.particles.size();

    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 8;
    config.max_particles_per_phase_bin = 16;
    config.max_weight_ratio = 8.0;
    config.max_support = 7;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, rank);

    TailMoment7 got;
    tail_moments_sum(tail, got);
    const std::uint64_t count_after = tail.particles.size();

    // Local ID uniqueness and creation-rank validity.
    bool ids_unique = true;
    std::set<std::uint64_t> ids;
    bool creation_rank_ok = true;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const std::uint64_t id = tail.particles[i].id;
        if (!ids.insert(id).second) ids_unique = false;
        if ((id >> 32) >= static_cast<std::uint64_t>(mpi_size)) {
            creation_rank_ok = false;
        }
    }

    // Global sums of the seven moments.
    const double local_moments[7] = {
        got.n - ref.n, got.px - ref.px, got.jx - ref.jx,
        got.ke - ref.ke, got.pixx - ref.pixx, got.piperp - ref.piperp,
        got.xw - ref.xw
    };
    double global_moments[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_moments, global_moments, 7, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double local_scales[7] = {
        std::max(1.0, std::fabs(ref.n)),
        std::max(1.0, std::fabs(ref.px)),
        std::max(1.0, std::fabs(ref.jx)),
        std::max(1.0, std::fabs(ref.ke)),
        std::max(1.0, std::fabs(ref.pixx)),
        std::max(1.0, std::fabs(ref.piperp)),
        std::max(1.0, std::fabs(ref.xw))
    };
    double global_scales[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_scales, global_scales, 7, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double residuals[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    for (int i = 0; i < 7; ++i) {
        residuals[i] = std::fabs(global_moments[i]) /
                       std::max(1.0, global_scales[i]);
    }

    std::uint64_t global_before = 0;
    std::uint64_t global_after = 0;
    MPI_Allreduce(&count_before, &global_before, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&count_after, &global_after, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);

    int local_pass =
        d.applied && d.groups_compressed >= 1 &&
        d.compression_fallback_count == 0 && ids_unique && creation_rank_ok &&
        residuals[0] <= 1.0e-9 && residuals[1] <= 1.0e-9 &&
        residuals[2] <= 1.0e-9 && residuals[3] <= 1.0e-9 &&
        residuals[4] <= 1.0e-9 && residuals[5] <= 1.0e-9 &&
        residuals[6] <= 1.0e-9 && global_after < global_before &&
        d.particles_after_local <= d.particles_before_local;
    int pass = local_pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &pass, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "global_before=" << global_before
                  << " global_after=" << global_after
                  << " groups_compressed_per_rank=" << d.groups_compressed
                  << " fallbacks=" << d.compression_fallback_count
                  << " residual_n=" << std::setprecision(17) << residuals[0]
                  << " residual_px=" << residuals[1]
                  << " residual_jx=" << residuals[2]
                  << " residual_ke=" << residuals[3]
                  << " residual_pixx=" << residuals[4]
                  << " residual_piperp=" << residuals[5]
                  << " residual_xw=" << residuals[6] << "\n";
        if (!args.result_path.empty()) {
            std::ofstream out(args.result_path.c_str(), std::ios::app);
            if (out) {
                out << "ranks=" << mpi_size
                    << " pass=" << (pass ? 1 : 0)
                    << " global_before=" << global_before
                    << " global_after=" << global_after
                    << " residual_n=" << std::setprecision(17) << residuals[0]
                    << " residual_px=" << residuals[1]
                    << " residual_jx=" << residuals[2]
                    << " residual_ke=" << residuals[3]
                    << " residual_pixx=" << residuals[4]
                    << " residual_piperp=" << residuals[5]
                    << " residual_xw=" << residuals[6] << "\n";
            }
        }
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
