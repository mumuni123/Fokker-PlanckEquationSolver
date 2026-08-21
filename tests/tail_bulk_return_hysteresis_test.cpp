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
    if (pass) tail_return_test::add_representable_cloud(
        tail, bulk, partition, slot.first, slot.second, 10.5 * grid.dx,
        3.0e12, 100, 0);

    TailBulkReturnConfig config;
    config.enabled = true;
    config.return_energy_mev = 5.5;
    config.residence_steps = 3;
    config.max_stencil_radius = 3;
    TailBulkReturn op(config);
    TailBulkReturnDiagnostics d;
    pass = pass && op.apply(bulk, tail, grid, partition, 1, 0, 1, d) &&
        !tail.particles.empty() &&
        tail.particles[0].return_residence_steps == 1;
    pass = pass && op.apply(bulk, tail, grid, partition, 2, 0, 1, d) &&
        !tail.particles.empty() &&
        tail.particles[0].return_residence_steps == 2;
    std::vector<double> original_ux(tail.particles.size(), 0.0);
    if (pass) for (size_t i = 0; i < tail.particles.size(); ++i) {
        original_ux[i] = tail.particles[i].ux;
        tail.particles[i].ux = 20.0;
    }
    pass = pass && op.apply(bulk, tail, grid, partition, 3, 0, 1, d) &&
        !tail.particles.empty() &&
        tail.particles[0].return_residence_steps == 0;
    if (pass) for (size_t i = 0; i < tail.particles.size(); ++i)
        tail.particles[i].ux = original_ux[i];
    pass = pass && op.apply(bulk, tail, grid, partition, 4, 0, 1, d) &&
        !tail.particles.empty();
    pass = pass && op.apply(bulk, tail, grid, partition, 5, 0, 1, d) &&
        !tail.particles.empty();
    pass = pass && op.apply(bulk, tail, grid, partition, 6, 0, 1, d) &&
        tail.particles.empty() && d.particles_removed >= 6 && d.committed;
    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"reset_above_threshold", 1.0},
            {"removed_after_residence", tail.particles.empty() ? 1.0 : 0.0},
            {"resident_particles", static_cast<double>(d.resident_particles)},
            {"attempted_groups", static_cast<double>(d.attempted_groups)},
            {"deferred_groups", static_cast<double>(d.deferred_infeasible_groups)},
            {"moment_residual_max", tail_return_test::worst_residual(d)}}, pass);
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
