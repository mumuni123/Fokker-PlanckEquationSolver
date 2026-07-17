#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <vector>

namespace {
double max_owned_difference(const Species& a, const Species& b,
                            const SpatialGrid& grid)
{
    double value = 0.0;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const size_t begin = static_cast<size_t>(grid.nghost + ix) * Param::Nvmu;
        for (size_t q = 0; q < Param::Nvmu; ++q)
            value = std::max(value, std::fabs(a.f[begin + q] - b.f[begin + q]));
    }
    return value;
}

double max_vector_difference(const std::vector<double>& a,
                             const std::vector<double>& b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double value = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        value = std::max(value, std::fabs(a[i] - b[i]));
    return value;
}

double max_vector_magnitude(const std::vector<double>& values)
{
    double value = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
        value = std::max(value, std::fabs(values[i]));
    return value;
}

double max_owned_magnitude(const Species& species, const SpatialGrid& grid)
{
    double value = 0.0;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const size_t begin = static_cast<size_t>(grid.nghost + ix) * Param::Nvmu;
        for (size_t q = 0; q < Param::Nvmu; ++q)
            value = std::max(value, std::fabs(species.f[begin + q]));
    }
    return value;
}

double max_particle_difference(const BeamPIC& a, const BeamPIC& b)
{
    if (a.particles.size() != b.particles.size())
        return std::numeric_limits<double>::infinity();
    double value = 0.0;
    for (size_t i = 0; i < a.particles.size(); ++i) {
        value = std::max(value, std::fabs(a.particles[i].x - b.particles[i].x));
        value = std::max(value, std::fabs(a.particles[i].px - b.particles[i].px));
        value = std::max(value, std::fabs(a.particles[i].weight - b.particles[i].weight));
    }
    return value;
}

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
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    SpatialGrid grid; grid.init(rank, size);
    Species seed;
    seed.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
              Param::dens, Param::temperature_e, false, grid);
    seed.initialize_maxwellian();
    EMFields field_seed; field_seed.init(grid);
    for (int i = 0; i < grid.nx_local; ++i) {
        const double phase = 2.0 * Const::pi * (grid.ix_start + i + 0.5) / Param::nx;
        field_seed.Ex_face[static_cast<size_t>(i)] = 1.0e6 * std::sin(phase);
    }
    field_seed.sync_cell_ex_from_faces(rank, size);
    BeamPIC beam_seed; beam_seed.init(grid);
    VlasovAmpereMidpointSolver solver;
    solver.set_nonuniform_high_order_enabled(true);
    solver.set_fct_enabled(true);
    solver.set_max_midpoint_iterations(40);
    solver.synchronize_background_ghosts(seed, grid, rank, size);
    const double dt = 0.05 * Param::dx / Const::c;
    // Exercise a nonempty Beam state and a non-initial RNG/remainder before
    // checkpointing; the next production step also performs a fresh injection.
    advance_beam_one_step(beam_seed, grid, field_seed, dt, dt, rank, size);

    Species direct_state = seed;
    EMFields direct_field = field_seed;
    const VlasovAmpereMidpointSolver::Result direct =
        solver.advance_background_and_fields(direct_state, beam_seed, direct_field,
                                             grid, dt, 2.0 * dt, rank, size);

    const CheckpointControlState control = {1, dt, dt, 0.0, 0, 0.0};
    std::string error;
    bool ok = write_checkpoint("restart_one_step_equivalence_tmp", control, seed,
                               beam_seed, field_seed, grid, rank, size, error);
    Species restarted;
    restarted.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
                   Param::dens, Param::temperature_e, false, grid);
    BeamPIC restarted_beam; restarted_beam.init(grid);
    EMFields restarted_field; restarted_field.init(grid);
    CheckpointControlState restored;
    ok = ok && read_checkpoint("restart_one_step_equivalence_tmp", restored,
                               restarted, restarted_beam, restarted_field,
                               grid, rank, size, error);
    solver.synchronize_background_ghosts(restarted, grid, rank, size);
    restarted_beam.deposit_density(grid, rank, size);
    const VlasovAmpereMidpointSolver::Result resumed =
        solver.advance_background_and_fields(restarted, restarted_beam,
                                             restarted_field, grid, restored.dt_s,
                                             restored.time_s + restored.dt_s,
                                             rank, size);
    double differences[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0};
    differences[0] = max_owned_difference(direct.species_np1,
                                           resumed.species_np1, grid);
    for (int i = 0; i < grid.nx_local; ++i) {
        differences[1] = std::max(differences[1],
            std::fabs(direct.fields_np1.Ex_face[static_cast<size_t>(i)] -
                      resumed.fields_np1.Ex_face[static_cast<size_t>(i)]));
    }
    differences[2] = max_particle_difference(direct.beam_np1, resumed.beam_np1);
    differences[3] = max_vector_difference(direct.j_bkg_face_mid,
                                            resumed.j_bkg_face_mid);
    differences[4] = max_vector_difference(direct.j_beam_face_mid,
                                            resumed.j_beam_face_mid);
    differences[5] = max_vector_difference(direct.j_total_face_mid,
                                            resumed.j_total_face_mid);
    differences[6] = std::fabs(direct.species_np1.total_particle_number() -
                               resumed.species_np1.total_particle_number());
    differences[7] = std::fabs(direct.beam_np1.total_particle_number(grid) -
                               resumed.beam_np1.total_particle_number(grid));
    differences[8] = std::fabs(
        direct.species_np1.total_kinetic_energy() + direct.beam_np1.total_kinetic_energy() +
        direct.fields_np1.total_energy() - resumed.species_np1.total_kinetic_energy() -
        resumed.beam_np1.total_kinetic_energy() - resumed.fields_np1.total_energy());
    differences[9] = std::max(std::fabs(direct.delta_ke_bkg - resumed.delta_ke_bkg),
                              std::fabs(direct.delta_ke_beam - resumed.delta_ke_beam));
    double scales[10] = {
        max_owned_magnitude(direct.species_np1, grid),
        max_vector_magnitude(direct.fields_np1.Ex_face),
        1.0,
        max_vector_magnitude(direct.j_bkg_face_mid),
        max_vector_magnitude(direct.j_beam_face_mid),
        max_vector_magnitude(direct.j_total_face_mid),
        std::fabs(direct.species_np1.total_particle_number()),
        std::fabs(direct.beam_np1.total_particle_number(grid)),
        std::fabs(direct.species_np1.total_kinetic_energy()) +
            std::fabs(direct.beam_np1.total_kinetic_energy()) +
            std::fabs(direct.fields_np1.total_energy()),
        std::max(std::fabs(direct.delta_ke_bkg),
                 std::fabs(direct.delta_ke_beam))
    };
    MPI_Allreduce(MPI_IN_PLACE, differences, 10, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, scales, 10, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    double relative_difference = 0.0;
    for (int i = 0; i < 10; ++i) {
        relative_difference = std::max(
            relative_difference, differences[i] / std::max(1.0, scales[i]));
    }
    const int local_ok = ok && direct.state_advanced && resumed.state_advanced &&
                         !direct.failed && !resumed.failed &&
                         direct.nonlinear_iterations == resumed.nonlinear_iterations &&
                         direct.converged == resumed.converged &&
                         direct.soft_accepted == resumed.soft_accepted &&
                         direct.soft_unconverged == resumed.soft_unconverged &&
                         differences[2] == 0.0 &&
                         relative_difference <= 1.0e-13;
    double global_totals[6] = {
        direct.species_np1.total_particle_number(),
        resumed.species_np1.total_particle_number(),
        direct.beam_np1.total_particle_number(grid),
        resumed.beam_np1.total_particle_number(grid),
        direct.species_np1.total_kinetic_energy() +
            direct.beam_np1.total_kinetic_energy() + direct.fields_np1.total_energy(),
        resumed.species_np1.total_kinetic_energy() +
            resumed.beam_np1.total_kinetic_energy() + resumed.fields_np1.total_energy()
    };
    MPI_Allreduce(MPI_IN_PLACE, global_totals, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf(
            "restart_one_step_equivalence f_difference=%.17e Ex_face_difference=%.17e "
            "beam_particle_difference=%.17e JN_difference=%.17e Jbeam_difference=%.17e "
            "Jtotal_difference=%.17e N_bkg_difference=%.17e N_beam_difference=%.17e "
            "total_energy_difference=%.17e delta_K_difference=%.17e "
            "N_bkg_direct=%.17e N_bkg_restart=%.17e N_beam_direct=%.17e "
            "N_beam_restart=%.17e energy_direct=%.17e energy_restart=%.17e "
            "nonlinear_iterations=%d acceptance_match=%d relative_difference=%.17e status=%s\n",
            differences[0], differences[1], differences[2], differences[3],
            differences[4], differences[5], differences[6], differences[7],
            differences[8], differences[9], global_totals[0], global_totals[1],
            global_totals[2], global_totals[3], global_totals[4], global_totals[5],
            direct.nonlinear_iterations,
            (direct.converged == resumed.converged &&
             direct.soft_accepted == resumed.soft_accepted &&
             direct.soft_unconverged == resumed.soft_unconverged) ? 1 : 0,
            relative_difference, global_ok ? "PASS" : "FAIL");
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
