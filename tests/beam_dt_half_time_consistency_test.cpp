#include "beam_pic.h"
#include "maxwell.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>

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
    EMFields fields;
    fields.init(grid);
    fields.sync_cell_ex_from_faces(rank, size);
    const double dt = 0.05 * Param::dx / Const::c;
    const double half_dt = 0.5 * dt;

    BeamPIC full;
    full.init(grid);
    advance_beam_one_step(full, grid, fields, dt, dt, rank, size);
    const double full_n_in = full.last_injected_number();
    const double full_j_beam = full.last_injected_current();
    const double full_time = dt;

    BeamPIC split;
    split.init(grid);
    advance_beam_one_step(split, grid, fields, half_dt, half_dt, rank, size);
    const double split_first_n_in = split.last_injected_number();
    const double split_first_j_beam = split.last_injected_current();
    advance_beam_one_step(split, grid, fields, half_dt, dt, rank, size);
    const double split_second_n_in = split.last_injected_number();
    const double split_second_j_beam = split.last_injected_current();
    const double split_time = half_dt + half_dt;
    const double split_total_n_in = split_first_n_in + split_second_n_in;
    const double split_cumulative_increment =
        split.cumulative_injected_energy() /
        (Const::me * Const::c * Const::c *
         (std::sqrt(1.0 + Param::gambetab * Param::gambetab) - 1.0));

    const double time_difference = std::fabs(full_time - split_time);
    const double full_split_injection_difference =
        std::fabs(full_n_in - split_total_n_in);
    const double injection_ledger_difference =
        std::fabs(split_total_n_in - split_cumulative_increment);
    const double current_impulse_difference = std::fabs(
        full_j_beam * dt -
        (split_first_j_beam + split_second_j_beam) * half_dt);
    const double injection_scale = std::max(
        1.0, std::max(std::fabs(full_n_in), std::max(std::fabs(split_total_n_in),
                                                      std::fabs(split_cumulative_increment))));
    const double current_impulse_scale = std::max(
        1.0, std::max(std::fabs(full_j_beam * dt),
                       std::fabs((split_first_j_beam + split_second_j_beam) * half_dt)));
    const int local_ok = std::isfinite(full_n_in) && std::isfinite(full_j_beam) &&
        std::isfinite(split_first_n_in) && std::isfinite(split_second_n_in) &&
        std::isfinite(split_first_j_beam) && std::isfinite(split_second_j_beam) &&
        time_difference <= 16.0 * std::numeric_limits<double>::epsilon() * dt &&
        full_split_injection_difference <=
            4096.0 * std::numeric_limits<double>::epsilon() * injection_scale &&
        injection_ledger_difference <=
            4096.0 * std::numeric_limits<double>::epsilon() * injection_scale &&
        current_impulse_difference <=
            4096.0 * std::numeric_limits<double>::epsilon() * current_impulse_scale;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream result("beam_dt_half_time_consistency.result");
        result << std::setprecision(17)
            << "time_A " << full_time << "\n"
            << "time_D " << split_time << "\n"
            << "time_difference " << time_difference << "\n"
            << "A_N_in " << full_n_in << "\n"
            << "D_N_in_half1 " << split_first_n_in << "\n"
            << "D_N_in_half2 " << split_second_n_in << "\n"
            << "D_N_in_total " << split_total_n_in << "\n"
            << "A_D_N_in_difference " << full_split_injection_difference << "\n"
            << "D_ledger_difference " << injection_ledger_difference << "\n"
            << "A_Jbeam " << full_j_beam << "\n"
            << "D_Jbeam_half1 " << split_first_j_beam << "\n"
            << "D_Jbeam_half2 " << split_second_j_beam << "\n"
            << "current_impulse_difference " << current_impulse_difference << "\n"
            << "status " << (global_ok ? "PASS" : "FAIL") << "\n";
        std::printf(
            "beam_dt_half_time_consistency time_A=%.17e time_D=%.17e "
            "time_difference=%.17e A_N_in=%.17e D_N_in_half1=%.17e "
            "D_N_in_half2=%.17e D_N_in_total=%.17e A_D_N_in_difference=%.17e "
            "D_ledger_difference=%.17e A_Jbeam=%.17e D_Jbeam_half1=%.17e "
            "D_Jbeam_half2=%.17e current_impulse_difference=%.17e status=%s\n",
            full_time, split_time, time_difference, full_n_in, split_first_n_in,
            split_second_n_in, split_total_n_in, full_split_injection_difference,
            injection_ledger_difference, full_j_beam, split_first_j_beam,
            split_second_j_beam, current_impulse_difference,
            global_ok ? "PASS" : "FAIL");
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
