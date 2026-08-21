#include "tail_bulk_return_test_common.h"

#include <mpi.h>
#include <iostream>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result;
    bool pass = size == 1 &&
        tail_return_test::parse_result_arg(argc, argv, result);
    SpatialGrid grid; Species bulk; HybridVelocityPartition partition;
    BackgroundTailPIC tail;
    tail_return_test::init_state(rank, size, grid, bulk, partition, tail);
    const std::pair<int, int> slot =
        tail_return_test::safe_velocity_slot(bulk, partition);
    if (slot.first < 0) pass = false;
    if (pass) {
        tail_return_test::add_representable_cloud(
            tail, bulk, partition, slot.first, slot.second,
            10.5 * grid.dx, 2.0e12, 100, 1);
        tail_return_test::add_representable_cloud(
            tail, bulk, partition, slot.first, slot.second,
            12.0 * grid.dx, 3.0e12, 200, 1);
    }
    TailBulkReturnConfig config;
    config.enabled = true; config.return_energy_mev = 5.5;
    config.residence_steps = 1; config.max_stencil_radius = 3;
    TailBulkReturnDiagnostics d;
    const double ideal_representation_tolerance =
        tail_bulk_return_representation_tolerance(0.0);
    const double audited_grid_representation_tolerance =
        tail_bulk_return_representation_tolerance(3.33354740672649e-2);
    const double capped_representation_tolerance =
        tail_bulk_return_representation_tolerance(2.0e-1);
    pass = pass &&
        std::fabs(ideal_representation_tolerance - 5.0e-3) <= 1.0e-15 &&
        audited_grid_representation_tolerance > 3.33354740672649e-2 &&
        audited_grid_representation_tolerance < 3.6e-2 &&
        std::fabs(capped_representation_tolerance - 5.0e-2) <= 1.0e-15;
    pass = pass && TailBulkReturn(config).apply(
        bulk, tail, grid, partition, 1, 0, 1, d);
    const double worst = tail_return_test::worst_residual(d);
    const double signed_invariant_relative = std::max(
        std::fabs(d.number_difference) / std::max(1.0, std::fabs(d.number)),
        std::max(
            std::fabs(d.px_difference) / std::max(1.0, std::fabs(d.px)),
            std::fabs(d.energy_difference) /
                std::max(1.0, std::fabs(d.energy))));
    const bool nonnegative = std::find_if(bulk.f.begin(), bulk.f.end(),
        [](double x) { return !std::isfinite(x) || x < 0.0; }) == bulk.f.end();
    pass = pass && tail.particles.empty() && d.particles_removed >= 12 &&
        d.committed_groups == 2 && d.deferred_infeasible_groups == 0 &&
        worst <= 1.0e-12 && signed_invariant_relative <= 1.0e-12 &&
        nonnegative;

    // Radius zero is intentionally infeasible and must defer without changing
    // either representation (apart from the already-satisfied residence age).
    Species bulk2; BackgroundTailPIC tail2; SpatialGrid grid2;
    HybridVelocityPartition partition2;
    tail_return_test::init_state(0, 1, grid2, bulk2, partition2, tail2);
    tail2.particles.push_back(tail_return_test::make_particle(
        bulk2, slot.first, slot.second, 14.5 * grid2.dx, 4.0e12, 21,
        std::numeric_limits<std::uint32_t>::max()));
    const std::vector<double> f2 = bulk2.f;
    const std::vector<BackgroundTailParticle> p2 = tail2.particles;
    config.max_stencil_radius = 0;
    TailBulkReturnDiagnostics deferred;
    const bool defer_ok = TailBulkReturn(config).apply(
        bulk2, tail2, grid2, partition2, 1, 0, 1, deferred);
    const bool deferred_unchanged = tail_return_test::equal_doubles(f2, bulk2.f) &&
        tail_return_test::equal_particles(p2, tail2.particles);
    pass = pass && defer_ok && deferred.deferred_infeasible_groups == 1 &&
        deferred.particles_removed == 0 && deferred_unchanged;

    // A real checkpoint group contains a broad velocity cloud.  Its support
    // cannot be represented by one small stencil around the mean velocity.
    // Use continuous off-grid particle velocities spread well beyond radius
    // three. This is the case for which exact six-moment matching on fixed
    // Eulerian cell centres is generally impossible, while N/Px/K remain the
    // exact representation invariants.
    Species bulk3; BackgroundTailPIC tail3; SpatialGrid grid3;
    HybridVelocityPartition partition3;
    tail_return_test::init_state(0, 1, grid3, bulk3, partition3, tail3);
    const int j0 = Param::Nv / 2;
    const int k0 = Param::Nmu / 5;
    const int dj_values[9] = { -20, -15, -10, -5, 0, 5, 10, 15, 20 };
    const int dk_values[5] = { -8, -4, 0, 4, 8 };
    std::uint64_t broad_id = 1000;
    for (int a = 0; a < 9; ++a) {
        for (int b = 0; b < 5; ++b) {
            const int j = j0 + dj_values[a];
            const int k = k0 + dk_values[b];
            if (j < 0 || j >= Param::Nv || k < 0 || k >= Param::Nmu)
                continue;
            const size_t id = idx2(j, k);
            if (partition3.bulk_owned_cell[id] == 0 ||
                partition3.kinetic_energy[id] >= 5.0e6 * Const::eV)
                continue;
            BackgroundTailParticle particle = tail_return_test::make_particle(
                bulk3, j, k, 18.25 * grid3.dx,
                1.0e12 * static_cast<double>(1 + a + b), broad_id++,
                std::numeric_limits<std::uint32_t>::max());
            const int jn = std::min(Param::Nv - 1, j + 1);
            const int kn = std::min(Param::Nmu - 1, k + 1);
            particle.ux = 0.63 * particle.ux +
                          0.37 * bulk3.cgrid.upar_cells[jn];
            particle.uy = 0.71 * particle.uy +
                          0.29 * bulk3.cgrid.uperp_cells[kn];
            tail3.particles.push_back(particle);
        }
    }
    const size_t broad_count = tail3.particles.size();
    std::vector<double> broad_histogram_reference(Param::Nvmu, 0.0);
    for (size_t p = 0; p < tail3.particles.size(); ++p) {
        int best_j = 0, best_k = 0;
        for (int j = 1; j < Param::Nv; ++j)
            if (std::fabs(bulk3.cgrid.upar_cells[j] - tail3.particles[p].ux) <
                std::fabs(bulk3.cgrid.upar_cells[best_j] - tail3.particles[p].ux))
                best_j = j;
        const double uperp = std::sqrt(tail3.particles[p].uy *
                                       tail3.particles[p].uy +
                                       tail3.particles[p].uz *
                                       tail3.particles[p].uz);
        for (int k = 1; k < Param::Nmu; ++k)
            if (std::fabs(bulk3.cgrid.uperp_cells[k] - uperp) <
                std::fabs(bulk3.cgrid.uperp_cells[best_k] - uperp))
                best_k = k;
        broad_histogram_reference[idx2(best_j, best_k)] +=
            tail3.particles[p].weight;
    }
    config.max_stencil_radius = 3;
    TailBulkReturnDiagnostics broad;
    const bool broad_ok = broad_count > 24 && TailBulkReturn(config).apply(
        bulk3, tail3, grid3, partition3, 1, 0, 1, broad);
    const bool broad_nonnegative = std::find_if(
        bulk3.f.begin(), bulk3.f.end(),
        [](double x) { return !std::isfinite(x) || x < 0.0; }) == bulk3.f.end();
    std::vector<double> broad_histogram_result(Param::Nvmu, 0.0);
    for (int ix = 0; ix < grid3.nx_local; ++ix) {
        const int sx = grid3.nghost + ix;
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                broad_histogram_result[idx2(j, k)] +=
                    bulk3.f[idx3(sx, j, k)];
    }
    double broad_histogram_l1 = 0.0, broad_histogram_scale = 0.0;
    for (size_t q = 0; q < broad_histogram_result.size(); ++q) {
        broad_histogram_l1 += std::fabs(
            broad_histogram_result[q] - broad_histogram_reference[q]);
        broad_histogram_scale += std::fabs(broad_histogram_reference[q]);
    }
    const double broad_histogram_relative_l1 = broad_histogram_l1 /
        std::max(1.0, broad_histogram_scale);
    const double broad_piperp_removed =
        broad.piperp_dx - broad.piperp_difference;
    const double broad_piperp_relative_expected =
        std::fabs(broad.piperp_difference) /
        std::max(std::fabs(broad_piperp_removed),
                 std::numeric_limits<double>::min());
    const bool broad_piperp_normalization_ok =
        std::fabs(broad.piperp_residual - broad_piperp_relative_expected) <=
        128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, broad_piperp_relative_expected);
    pass = pass && broad_ok && tail3.particles.empty() &&
        broad.particles_removed == broad_count &&
        broad.committed_groups == 1 &&
        tail_return_test::invariant_residual(broad) <= 1.0e-12 &&
        tail_return_test::representation_residual(broad) <= 5.0e-3 &&
        broad_histogram_relative_l1 <= 1.0e-1 &&
        broad_piperp_normalization_ok &&
        broad_nonnegative;

    // Real R2 groups may contain only one continuous-velocity particle. Its
    // moment point is not a velocity-cell centre, so this exercises the
    // deterministic geometric N/Px/K fallback rather than a grid-aligned
    // representable cloud.
    Species bulk4; BackgroundTailPIC tail4; SpatialGrid grid4;
    HybridVelocityPartition partition4;
    tail_return_test::init_state(0, 1, grid4, bulk4, partition4, tail4);
    BackgroundTailParticle continuous;
    continuous.x = 21.25 * grid4.dx;
    continuous.ux = 8.0137;
    continuous.uy = 8.2019;
    continuous.uz = 0.7311;
    continuous.weight = 7.0e12;
    continuous.id = 9001;
    continuous.return_residence_steps =
        std::numeric_limits<std::uint32_t>::max();
    tail4.particles.push_back(continuous);
    TailBulkReturnDiagnostics single;
    const bool single_ok = TailBulkReturn(config).apply(
        bulk4, tail4, grid4, partition4, 1, 0, 1, single);
    pass = pass && single_ok && tail4.particles.empty() &&
        single.particles_removed == 1 && single.committed_groups == 1 &&
        tail_return_test::invariant_residual(single) <= 1.0e-12 &&
        tail_return_test::representation_residual(single) <= 5.0e-3;
    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"moment_residual_max", worst},
            {"signed_invariant_difference_relative_max",
             signed_invariant_relative},
            {"signed_jx_difference", d.jx_difference},
            {"signed_pixx_difference", d.pixx_difference},
            {"signed_piperp_difference", d.piperp_difference},
            {"nonnegative_bulk", nonnegative ? 1.0 : 0.0},
            {"infeasible_deferred", deferred_unchanged ? 1.0 : 0.0},
            {"broad_cloud_particles", static_cast<double>(broad_count)},
            {"broad_cloud_removed", static_cast<double>(broad.particles_removed)},
            {"broad_cloud_invariant_residual_max",
             tail_return_test::invariant_residual(broad)},
            {"broad_cloud_representation_residual_max",
             tail_return_test::representation_residual(broad)},
            {"broad_cloud_velocity_histogram_relative_l1",
             broad_histogram_relative_l1},
            {"broad_cloud_piperp_relative_expected",
             broad_piperp_relative_expected},
            {"broad_cloud_piperp_normalization_ok",
             broad_piperp_normalization_ok ? 1.0 : 0.0},
            {"continuous_single_removed",
             static_cast<double>(single.particles_removed)},
            {"continuous_single_invariant_residual_max",
             tail_return_test::invariant_residual(single)},
            {"continuous_single_representation_residual_max",
             tail_return_test::representation_residual(single)},
            {"ideal_representation_tolerance",
             ideal_representation_tolerance},
            {"audited_grid_representation_tolerance",
             audited_grid_representation_tolerance},
            {"capped_representation_tolerance",
             capped_representation_tolerance}}, pass);
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
