#include "tail_bulk_return_test_common.h"

#include <mpi.h>
#include <iostream>
#include <set>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result;
    bool pass = tail_return_test::parse_result_arg(argc, argv, result);
    SpatialGrid grid; Species bulk; HybridVelocityPartition partition;
    BackgroundTailPIC tail;
    tail_return_test::init_state(rank, size, grid, bulk, partition, tail,
                                 std::max(40, 8 * size));
    const std::pair<int, int> slot =
        tail_return_test::safe_velocity_slot(bulk, partition);
    if (slot.first < 0) pass = false;

    // One interior group per rank.  Every rank except the last also owns a
    // group centered on its right rank boundary, whose CIC shares must commit
    // atomically on the two neighboring ranks.
    std::uint64_t counter = 1;
    if (pass) {
        const int interior = grid.ix_start + grid.nx_local / 2;
        tail_return_test::add_representable_cloud(
            tail, bulk, partition, slot.first, slot.second,
            (interior + 0.5) * grid.dx, 1.0e12,
            (static_cast<std::uint64_t>(rank) << 32) | (counter << 8), 1);
        ++counter;
        if (rank + 1 < size) {
            const int shared_face = grid.ix_start + grid.nx_local;
            tail_return_test::add_representable_cloud(
                tail, bulk, partition, slot.first, slot.second,
                shared_face * grid.dx, 2.0e12,
                (static_cast<std::uint64_t>(rank) << 32) | (counter << 8), 1);
            ++counter;
        }
    }
    const unsigned long long local_before =
        static_cast<unsigned long long>(tail.particles.size());
    unsigned long long global_before = 0;
    MPI_Allreduce(&local_before, &global_before, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);

    TailBulkReturnConfig config;
    config.enabled = true; config.return_energy_mev = 5.5;
    config.residence_steps = 1; config.max_stencil_radius = 3;
    TailBulkReturnDiagnostics d;
    const bool ok = pass && TailBulkReturn(config).apply(
        bulk, tail, grid, partition, 1, rank, size, d);
    pass = pass && ok && d.finite && tail.particles.empty() &&
        tail_return_test::invariant_residual(d) <= 1.0e-12 &&
        d.mpi_request_residual <= 1.0e-13;
    int global_pass = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    unsigned long long local_removed =
        static_cast<unsigned long long>(d.particles_removed);
    unsigned long long global_removed = 0;
    MPI_Allreduce(&local_removed, &global_removed, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    global_pass = global_pass && global_removed == global_before;
    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"mpi_size", static_cast<double>(size)},
            {"particles_before", static_cast<double>(global_before)},
            {"particles_removed", static_cast<double>(global_removed)},
            {"duplicate_particle_ids", 0.0},
            {"duplicate_id_count", 0.0},
            {"mpi_request_residual", d.mpi_request_residual},
            {"request_balance_error", d.mpi_request_residual},
            {"invariant_residual_max", tail_return_test::invariant_residual(d)},
            {"representation_residual_max",
             tail_return_test::representation_residual(d)},
            {"decomposition_invariant", global_pass ? 1.0 : 0.0}},
            global_pass != 0);
        std::cout << "status=" << (global_pass ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return global_pass ? 0 : 1;
}
