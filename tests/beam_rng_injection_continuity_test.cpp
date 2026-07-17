#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <string>
#include <vector>

namespace {

void advance_beam_one_step(BeamPIC& beam, const SpatialGrid& grid,
                           const EMFields& fields, double dt, double time,
                           int rank, int size)
{
    beam.begin_step(grid, dt);
    beam.inject(grid, fields, dt, time, rank, size);
    beam.push(grid, fields, dt, rank, size);
    beam.deposit_density(grid, rank, size);
    beam.finalize_charge_conserving_current(grid, dt, rank, size);
}

double max_vector_difference(const std::vector<double>& a,
                             const std::vector<double>& b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result = std::max(result, std::fabs(a[i] - b[i]));
    }
    return result;
}

double max_particle_difference(const std::vector<BeamParticle>& a,
                               const std::vector<BeamParticle>& b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result = std::max(result, std::fabs(a[i].x - b[i].x));
        result = std::max(result, std::fabs(a[i].px - b[i].px));
        result = std::max(result, std::fabs(a[i].weight - b[i].weight));
    }
    return result;
}

}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SpatialGrid grid;
    grid.init(rank, size);
    Species background;
    background.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                    Const::me, Param::dens, Param::temperature_e, false, grid);
    background.initialize_maxwellian();
    EMFields fields;
    fields.init(grid);
    fields.sync_cell_ex_from_faces(rank, size);

    // Produce a completed accepted-like state first.  This moves the RNG and
    // leaves a nontrivial injection remainder before the checkpoint is made.
    const double dt = 0.05 * Param::dx / Const::c;
    BeamPIC seed;
    seed.init(grid);
    advance_beam_one_step(seed, grid, fields, dt, dt, rank, size);
    const size_t checkpoint_particle_count = seed.particles.size();

    const CheckpointControlState control = {1, dt, dt, 0.0, 0, 0.0};
    std::string error;
    bool ok = write_checkpoint("beam_rng_injection_continuity_tmp", control,
                               background, seed, fields, grid, rank, size,
                               error);

    Species restored_background;
    restored_background.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
                            -Const::qe, Const::me, Param::dens,
                            Param::temperature_e, false, grid);
    BeamPIC restored;
    restored.init(grid);
    EMFields restored_fields;
    restored_fields.init(grid);
    CheckpointControlState restored_control;
    ok = ok && read_checkpoint("beam_rng_injection_continuity_tmp",
                               restored_control, restored_background, restored,
                               restored_fields, grid, rank, size, error);

    // Beam density is derived from particles after restart, exactly as in the
    // production restart path.  It must reproduce the completed checkpoint.
    restored.deposit_density(grid, rank, size);
    double local_difference = max_vector_difference(seed.density, restored.density);
    const BeamPersistentState seed_state = seed.export_persistent_state();
    const BeamPersistentState restored_checkpoint_state =
        restored.export_persistent_state();
    if (std::memcmp(&seed_state, &restored_checkpoint_state,
                    sizeof(BeamPersistentState)) != 0) {
        local_difference = std::numeric_limits<double>::infinity();
    }

    BeamPIC direct = seed;
    advance_beam_one_step(direct, grid, fields, dt, 2.0 * dt, rank, size);
    advance_beam_one_step(restored, grid, restored_fields, restored_control.dt_s,
                          restored_control.time_s + restored_control.dt_s,
                          rank, size);

    local_difference = std::max(local_difference,
                                max_particle_difference(direct.particles,
                                                        restored.particles));
    local_difference = std::max(local_difference,
                                max_vector_difference(direct.density,
                                                      restored.density));
    local_difference = std::max(local_difference,
                                max_vector_difference(direct.current_face_x,
                                                      restored.current_face_x));
    local_difference = std::max(local_difference,
                                max_vector_difference(direct.current_x,
                                                      restored.current_x));
    const BeamPersistentState direct_state = direct.export_persistent_state();
    const BeamPersistentState restored_state = restored.export_persistent_state();
    if (std::memcmp(&direct_state, &restored_state,
                    sizeof(BeamPersistentState)) != 0) {
        local_difference = std::numeric_limits<double>::infinity();
    }

    const long long direct_new_local =
        static_cast<long long>(direct.particles.size()) -
        static_cast<long long>(checkpoint_particle_count);
    const long long restored_new_local =
        static_cast<long long>(restored.particles.size()) -
        static_cast<long long>(checkpoint_particle_count);
    long long direct_new_global = 0;
    long long restored_new_global = 0;
    MPI_Allreduce(&direct_new_local, &direct_new_global, 1, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&restored_new_local, &restored_new_global, 1, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &local_difference, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    const int local_ok = ok && direct_new_local == restored_new_local &&
                         direct.last_injected_number() == restored.last_injected_number() &&
                         direct.last_injected_current() == restored.last_injected_current() &&
                         local_difference == 0.0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf(
            "beam_rng_injection_continuity injected_macro_particles=%lld "
            "restart_injected_macro_particles=%lld N_in=%.17e "
            "J_beam=%.17e max_difference=%.17e status=%s\n",
            direct_new_global, restored_new_global, direct.last_injected_number(),
            direct.last_injected_current(), local_difference,
            global_ok ? "PASS" : "FAIL");
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
