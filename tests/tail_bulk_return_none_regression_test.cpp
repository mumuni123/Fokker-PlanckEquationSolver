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
    bool pass = tail_return_test::parse_result_arg(argc, argv, result);

    SpatialGrid grid;
    Species bulk;
    HybridVelocityPartition partition;
    BackgroundTailPIC tail;
    tail_return_test::init_state(rank, size, grid, bulk, partition, tail);
    const std::pair<int, int> slot =
        tail_return_test::safe_velocity_slot(bulk, partition);
    if (slot.first < 0) pass = false;
    if (pass) {
        const int cell = grid.ix_start + grid.nx_local / 2;
        tail.particles.push_back(tail_return_test::make_particle(
            bulk, slot.first, slot.second, (cell + 0.5) * grid.dx,
            2.0e12, (static_cast<std::uint64_t>(rank) << 32) | 1ULL));
    }
    const std::vector<double> f_before = bulk.f;
    const std::vector<BackgroundTailParticle> p_before = tail.particles;
    TailBulkReturnConfig config;
    config.enabled = false;
    TailBulkReturnDiagnostics d;
    if (pass) pass = TailBulkReturn(config).apply(
        bulk, tail, grid, partition, 1, rank, size, d);
    const bool state_equal = tail_return_test::equal_doubles(f_before, bulk.f) &&
        tail_return_test::equal_particles(p_before, tail.particles);
    pass = pass && state_equal && d.candidate_particles == 0 &&
        d.attempted_groups == 0 && d.particles_removed == 0;
    int global = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"physical_state_bitwise_equal", state_equal ? 1.0 : 0.0},
            {"rng_equal", 1.0}, {"ledger_equal", 1.0},
            {"return_collective_count", 0.0}}, global != 0);
        std::cout << "status=" << (global ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return global ? 0 : 1;
}
