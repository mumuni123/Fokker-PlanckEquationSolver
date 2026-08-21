#include "tail_bulk_return_test_common.h"

#include <mpi.h>

#include <algorithm>
#include <iostream>

namespace {

std::pair<int, int> slot_in_energy_band(
    const Species& bulk, const HybridVelocityPartition& partition,
    double low_mev, double high_mev)
{
    const double low = low_mev * 1.0e6 * Const::eV;
    const double high = high_mev * 1.0e6 * Const::eV;
    for (int j = 1; j < Param::Nv - 1; ++j) {
        for (int k = 1; k < Param::Nmu - 1; ++k) {
            const size_t slot = idx2(j, k);
            if (partition.bulk_owned_cell[slot] != 0 &&
                partition.kinetic_energy[slot] >= low &&
                partition.kinetic_energy[slot] < high) {
                return std::make_pair(j, k);
            }
        }
    }
    return std::make_pair(-1, -1);
}

TailBulkReturnConfig config(double energy_mev, std::uint32_t residence)
{
    TailBulkReturnConfig result;
    result.enabled = true;
    result.return_energy_mev = energy_mev;
    result.residence_steps = residence;
    result.max_stencil_radius = 3;
    result.moment_tolerance = 1.0e-12;
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result;
    bool pass = size == 1 &&
        tail_return_test::parse_result_arg(argc, argv, result);

    SpatialGrid grid;
    Species bulk;
    HybridVelocityPartition partition;
    BackgroundTailPIC tail;
    tail_return_test::init_state(rank, size, grid, bulk, partition, tail);
    const std::pair<int, int> band = slot_in_energy_band(
        bulk, partition, 5.25, 5.50);
    pass = pass && band.first >= 0;

    // A particle that inherited age under 5.5 MeV must immediately lose that
    // age when a branch lowers K_in below the particle energy.
    bool lower_threshold_resets = false;
    if (pass) {
        tail.particles.push_back(tail_return_test::make_particle(
            bulk, band.first, band.second, 10.5 * grid.dx,
            2.0e12, 100, 7));
        TailBulkReturn op(config(5.25, 100));
        TailBulkReturnDiagnostics d;
        pass = pass && op.apply(bulk, tail, grid, partition, 1, 0, 1, d);
        lower_threshold_resets = tail.particles.size() == 1 &&
            tail.particles[0].return_residence_steps == 0 &&
            d.candidate_particles == 0 && d.particles_removed == 0;
        pass = pass && lower_threshold_resets;
    }

    // Raising K_in makes the same velocity band newly eligible.  Its age must
    // begin at one accepted step, never inherit a synthetic warm history.
    bool raised_threshold_starts_at_one = false;
    if (pass) {
        tail.particles[0].return_residence_steps = 0;
        TailBulkReturn op(config(5.75, 100));
        TailBulkReturnDiagnostics d;
        pass = pass && op.apply(bulk, tail, grid, partition, 2, 0, 1, d);
        raised_threshold_starts_at_one = tail.particles.size() == 1 &&
            tail.particles[0].return_residence_steps == 1 &&
            d.candidate_particles == 1 && d.resident_particles == 0;
        pass = pass && raised_threshold_starts_at_one;
    }

    // N_res changes use the checkpointed count itself.  A larger value keeps
    // the cloud; switching to a smaller value may return it on the next
    // accepted step without altering the stored age beforehand.
    bool residence_switch_uses_saved_age = false;
    if (pass) {
        tail.particles.clear();
        std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
        const std::pair<int, int> safe =
            tail_return_test::safe_velocity_slot(bulk, partition);
        pass = pass && safe.first >= 0;
        if (pass) tail_return_test::add_representable_cloud(
            tail, bulk, partition, safe.first, safe.second,
            10.5 * grid.dx, 3.0e12, 200, 7);
        TailBulkReturn op(config(5.5, 16));
        TailBulkReturnDiagnostics d16;
        pass = pass && op.apply(bulk, tail, grid, partition, 3, 0, 1, d16) &&
            !tail.particles.empty() &&
            tail.particles[0].return_residence_steps == 8;
        op.set_config(config(5.5, 4));
        TailBulkReturnDiagnostics d4;
        pass = pass && op.apply(bulk, tail, grid, partition, 4, 0, 1, d4);
        residence_switch_uses_saved_age = tail.particles.empty() &&
            d4.particles_removed > 0 && d4.committed;
        pass = pass && residence_switch_uses_saved_age;
    }

    if (rank == 0) {
        tail_return_test::write_result(result, {
            {"lower_threshold_resets_ineligible_age",
             lower_threshold_resets ? 1.0 : 0.0},
            {"raised_threshold_starts_age_at_one",
             raised_threshold_starts_at_one ? 1.0 : 0.0},
            {"residence_switch_uses_checkpointed_age",
             residence_switch_uses_saved_age ? 1.0 : 0.0}}, pass);
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << '\n';
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
