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
    double difference = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        difference = std::max(difference, std::fabs(a[i] - b[i]));
    return difference;
}

}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
    SpatialGrid grid; grid.init(rank, size);
    Species before; before.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                                Const::me, Param::dens, Param::temperature_e, false, grid);
    before.initialize_maxwellian();
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const size_t base = static_cast<size_t>(grid.nghost + ix) * Param::Nvmu;
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            before.f[base + q] *= 1.0 + 1.0e-6 *
                static_cast<double>((grid.ix_start + ix + static_cast<int>(q % 7)) % 11);
        }
    }
    BeamPIC beam_before; beam_before.init(grid);
    EMFields fields_before; fields_before.init(grid);
    for (int i=0; i<grid.nx_local; ++i) fields_before.Ex_face[static_cast<size_t>(i)] = 0.125 * (rank + 1) * (i + 1);
    fields_before.sync_cell_ex_from_faces(rank, size);
    const double beam_dt = 0.05 * Param::dx / Const::c;
    advance_beam_one_step(beam_before, grid, fields_before, beam_dt, beam_dt,
                          rank, size);
    CheckpointControlState c = {7, 2.5, 0.125, 3.0, 6, 9.0};
    std::string error;
    const std::string dir = "checkpoint_roundtrip_tmp";
    bool ok = write_checkpoint(dir, c, before, beam_before, fields_before, grid, rank, size, error);
    Species after; after.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                              Const::me, Param::dens, Param::temperature_e, false, grid);
    BeamPIC beam_after; beam_after.init(grid); EMFields fields_after; fields_after.init(grid);
    CheckpointControlState restored;
    ok = ok && read_checkpoint(dir, restored, after, beam_after, fields_after, grid, rank, size, error);
    double f_difference = 0.0;
    for (int ix=0; ix<grid.nx_local; ++ix) for (size_t q=0; q<Param::Nvmu; ++q)
        f_difference = std::max(f_difference, std::fabs(before.f[static_cast<size_t>(grid.nghost+ix)*Param::Nvmu+q] - after.f[static_cast<size_t>(grid.nghost+ix)*Param::Nvmu+q]));
    double ex_difference = 0.0;
    for (int i=0; i<grid.nx_local; ++i) ex_difference = std::max(ex_difference, std::fabs(fields_before.Ex_face[i]-fields_after.Ex_face[i]));
    const BeamPersistentState before_persistent = beam_before.export_persistent_state();
    const BeamPersistentState after_persistent = beam_after.export_persistent_state();
    const bool control_match = restored.step == c.step && restored.time_s == c.time_s &&
        restored.dt_s == c.dt_s && restored.next_snapshot_s == c.next_snapshot_s &&
        restored.last_snapshot_step == c.last_snapshot_step &&
        restored.cumulative_collision_energy_delta == c.cumulative_collision_energy_delta;
    if (beam_after.particles.size() != beam_before.particles.size() || !control_match ||
        std::memcmp(&before_persistent, &after_persistent,
                    sizeof(BeamPersistentState)) != 0) ok = false;
    double particle_difference = 0.0;
    for (size_t i = 0; i < beam_before.particles.size(); ++i) {
        particle_difference = std::max(particle_difference,
            std::fabs(beam_before.particles[i].x - beam_after.particles[i].x));
        particle_difference = std::max(particle_difference,
            std::fabs(beam_before.particles[i].px - beam_after.particles[i].px));
        particle_difference = std::max(particle_difference,
            std::fabs(beam_before.particles[i].weight - beam_after.particles[i].weight));
    }
    beam_after.deposit_density(grid, rank, size);
    const double density_difference = max_vector_difference(beam_before.density,
                                                            beam_after.density);
    double differences[4] = {f_difference, ex_difference, particle_difference,
                             density_difference};
    MPI_Allreduce(MPI_IN_PLACE, differences, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const double local_max = std::max(std::max(differences[0], differences[1]),
                                      std::max(differences[2], differences[3]));
    if (particle_difference != 0.0) ok = false;
    int int_ok = ok ? 1 : 0, all_ok = 0; MPI_Allreduce(&int_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) std::printf(
        "checkpoint_roundtrip f_difference=%.17e Ex_face_difference=%.17e "
        "beam_particle_difference=%.17e beam_density_reconstruction_difference=%.17e "
        "persistent_state_match=%d control_state_match=%d status=%s\n",
        differences[0], differences[1], differences[2], differences[3],
        std::memcmp(&before_persistent, &after_persistent,
                    sizeof(BeamPersistentState)) == 0 ? 1 : 0,
        control_match ? 1 : 0, (all_ok && local_max == 0.0) ? "PASS" : "FAIL");
    MPI_Finalize(); return (all_ok && local_max == 0.0) ? 0 : 1;
}
