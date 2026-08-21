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
        const BackgroundTailParticle p = tail_return_test::make_particle(
            bulk, slot.first, slot.second, 10.5 * grid.dx, 2.0e12, 99, 0);
        tail.particles.push_back(p);
        tail.particles.push_back(p); // deterministic protocol fault
    }
    const std::vector<double> f_before = bulk.f;
    const std::vector<BackgroundTailParticle> p_before = tail.particles;
    TailBulkReturnConfig config;
    config.enabled = true; config.return_energy_mev = 5.5;
    config.residence_steps = 1; config.max_stencil_radius = 3;
    TailBulkReturnDiagnostics d;
    const bool apply_ok = TailBulkReturn(config).apply(
        bulk, tail, grid, partition, 1, 0, 1, d);
    const bool unchanged = tail_return_test::equal_doubles(f_before, bulk.f) &&
        tail_return_test::equal_particles(p_before, tail.particles);
    pass = pass && !apply_ok && !d.finite && unchanged;
    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"fault_rejected", apply_ok ? 0.0 : 1.0},
            {"state_bitwise_equal", unchanged ? 1.0 : 0.0},
            {"accepted_state_unchanged", unchanged ? 1.0 : 0.0}}, pass);
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
