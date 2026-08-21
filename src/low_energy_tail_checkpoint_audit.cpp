#include "beam_pic.h"
#include "grid.h"
#include "low_energy_tail_audit.h"
#include "maxwell.h"
#include "species.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string checkpoint;
    std::string result;
    std::vector<double> thresholds;
    Options()
    {
        thresholds.push_back(5.5);
        thresholds.push_back(5.75);
        thresholds.push_back(6.0);
    }
};

bool parse_thresholds(const char* text, std::vector<double>& values)
{
    std::stringstream input(text);
    std::string token;
    values.clear();
    while (std::getline(input, token, ',')) {
        char* end = NULL;
        const double value = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || !std::isfinite(value) ||
            !(value > 0.0)) return false;
        values.push_back(value);
    }
    std::sort(values.begin(), values.end());
    return !values.empty() &&
           std::adjacent_find(values.begin(), values.end()) == values.end();
}

bool parse_options(int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if ((arg == "--checkpoint" || arg == "--result" ||
             arg == "--thresholds-mev") && i + 1 >= argc) return false;
        if (arg == "--checkpoint") options.checkpoint = argv[++i];
        else if (arg == "--result") options.result = argv[++i];
        else if (arg == "--thresholds-mev") {
            if (!parse_thresholds(argv[++i], options.thresholds)) return false;
        } else return false;
    }
    return !options.checkpoint.empty() && !options.result.empty();
}

std::uint64_t particle_hash(
    const std::vector<BackgroundTailParticle>& particles)
{
    const unsigned char* bytes = particles.empty() ? NULL :
        reinterpret_cast<const unsigned char*>(&particles[0]);
    const size_t count = particles.size() * sizeof(BackgroundTailParticle);
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t i = 0; i < count; ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

void reduce_result(LowEnergyTailAuditResult& result)
{
    double totals[2] = {result.total_number, result.total_kinetic_energy};
    MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.total_number = totals[0];
    result.total_kinetic_energy = totals[1];
    unsigned long long macro =
        static_cast<unsigned long long>(result.total_macro_particles);
    MPI_Allreduce(MPI_IN_PLACE, &macro, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    result.total_macro_particles = static_cast<std::uint64_t>(macro);
    int finite = result.finite ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &finite, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    result.finite = finite != 0;

    for (size_t k = 0; k < result.thresholds.size(); ++k) {
        LowEnergyTailThresholdStats& stats = result.thresholds[k];
        double values[4] = {stats.number, stats.kinetic_energy,
                            stats.outside_shape_number,
                            stats.outside_shape_energy};
        MPI_Allreduce(MPI_IN_PLACE, values, 4, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        stats.number = values[0];
        stats.kinetic_energy = values[1];
        stats.outside_shape_number = values[2];
        stats.outside_shape_energy = values[3];
        unsigned long long threshold_macro =
            static_cast<unsigned long long>(stats.macro_particles);
        MPI_Allreduce(MPI_IN_PLACE, &threshold_macro, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        stats.macro_particles =
            static_cast<std::uint64_t>(threshold_macro);
        MPI_Allreduce(MPI_IN_PLACE, stats.cell_number.data(),
                      static_cast<int>(stats.cell_number.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, stats.cell_kinetic_energy.data(),
                      static_cast<int>(stats.cell_kinetic_energy.size()),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, stats.cell_macro_supports.data(),
                      static_cast<int>(stats.cell_macro_supports.size()),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Options options;
    int parsed = parse_options(argc, argv, options) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &parsed, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!parsed) {
        if (rank == 0) std::cerr << "invalid low-energy Tail audit arguments\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(rank, size, Param::nx, Param::Lx);
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    VpfpCheckpointControl control;
    VpfpCheckpointTailState tail_state;
    std::string error;
    int read_ok = read_vpfp_checkpoint(options.checkpoint, control, electrons,
                                       beam, fields, grid, &tail_state,
                                       rank, size, error) && tail_state.present;
    MPI_Allreduce(MPI_IN_PLACE, &read_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!read_ok) {
        if (rank == 0) std::cerr << "low-energy Tail checkpoint read failed: "
                                 << error << "\n";
        MPI_Finalize();
        return 3;
    }

    const std::uint64_t hash_before = particle_hash(tail_state.tail.particles);
    LowEnergyTailAuditResult result = audit_low_energy_tail_local(
        tail_state.tail.particles, grid, options.thresholds);
    const std::uint64_t hash_after = particle_hash(tail_state.tail.particles);
    int unchanged = hash_before == hash_after ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &unchanged, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    reduce_result(result);

    if (rank == 0) {
        std::ofstream summary(options.result.c_str());
        std::ofstream cells((options.result + ".cells.dat").c_str());
        if (!summary || !cells) {
            std::cerr << "cannot write low-energy Tail audit output\n";
            MPI_Abort(MPI_COMM_WORLD, 4);
        }
        summary << std::setprecision(17)
                << "schema=low_energy_tail_checkpoint_audit_v1\n"
                << "checkpoint_step=" << control.step << "\n"
                << "checkpoint_time_s=" << control.time << "\n"
                << "mpi_size=" << size << "\n"
                << "audit_read_only_state_unchanged=" << unchanged << "\n"
                << "finite=" << (result.finite ? 1 : 0) << "\n"
                << "total_tail_number=" << result.total_number << "\n"
                << "total_tail_kinetic_energy="
                << result.total_kinetic_energy << "\n"
                << "total_tail_macro_particles="
                << result.total_macro_particles << "\n";
        for (size_t k = 0; k < result.thresholds.size(); ++k) {
            const LowEnergyTailThresholdStats& stats = result.thresholds[k];
            double deposited_number = stats.outside_shape_number;
            double deposited_energy = stats.outside_shape_energy;
            double core_number = 0.0;
            double max_cell_density = 0.0;
            std::uint64_t nonzero_cells = 0;
            const int core_begin = grid.nx_global / 10;
            const int core_end = grid.nx_global - core_begin;
            for (int ix = 0; ix < grid.nx_global; ++ix) {
                const double cell_number =
                    stats.cell_number[static_cast<size_t>(ix)];
                deposited_number += cell_number;
                deposited_energy +=
                    stats.cell_kinetic_energy[static_cast<size_t>(ix)];
                if (ix >= core_begin && ix < core_end)
                    core_number += cell_number;
                if (cell_number > 0.0) ++nonzero_cells;
                max_cell_density = std::max(
                    max_cell_density, cell_number / grid.dx);
            }
            summary << "threshold_" << k << "_mev=" << stats.threshold_mev
                    << "\nthreshold_" << k << "_number=" << stats.number
                    << "\nthreshold_" << k << "_number_fraction="
                    << stats.number / std::max(result.total_number, 1.0e-300)
                    << "\nthreshold_" << k << "_kinetic_energy="
                    << stats.kinetic_energy
                    << "\nthreshold_" << k << "_energy_fraction="
                    << stats.kinetic_energy /
                       std::max(result.total_kinetic_energy, 1.0e-300)
                    << "\nthreshold_" << k << "_macro_particles="
                    << stats.macro_particles
                    << "\nthreshold_" << k << "_outside_shape_number="
                    << stats.outside_shape_number
                    << "\nthreshold_" << k << "_outside_shape_energy="
                    << stats.outside_shape_energy
                    << "\nthreshold_" << k << "_shape_number_residual="
                    << (deposited_number - stats.number) /
                       std::max(stats.number, 1.0e-300)
                    << "\nthreshold_" << k << "_shape_energy_residual="
                    << (deposited_energy - stats.kinetic_energy) /
                       std::max(stats.kinetic_energy, 1.0e-300)
                    << "\nthreshold_" << k << "_core_number="
                    << core_number
                    << "\nthreshold_" << k << "_core_number_fraction="
                    << core_number / std::max(stats.number, 1.0e-300)
                    << "\nthreshold_" << k << "_nonzero_cell_count="
                    << nonzero_cells
                    << "\nthreshold_" << k << "_max_cell_density_m3="
                    << max_cell_density << "\n";
        }
        summary << "status=" <<
            (result.finite && unchanged ? "PASS" : "FAIL") << "\n";

        cells << std::setprecision(17)
              << "global_ix x_um threshold_mev number_m2 "
                 "density_m3 kinetic_energy_J_m2 energy_density_J_m3 "
                 "macro_supports\n";
        for (size_t k = 0; k < result.thresholds.size(); ++k) {
            const LowEnergyTailThresholdStats& stats = result.thresholds[k];
            for (int ix = 0; ix < grid.nx_global; ++ix) {
                cells << ix << " " << (ix + 0.5) * grid.dx * 1.0e6 << " "
                      << stats.threshold_mev << " "
                      << stats.cell_number[static_cast<size_t>(ix)] << " "
                      << stats.cell_number[static_cast<size_t>(ix)] /
                         grid.dx << " "
                      << stats.cell_kinetic_energy[static_cast<size_t>(ix)]
                      << " "
                      << stats.cell_kinetic_energy[static_cast<size_t>(ix)] /
                         grid.dx << " "
                      << stats.cell_macro_supports[static_cast<size_t>(ix)]
                      << "\n";
            }
        }
    }
    MPI_Finalize();
    return result.finite && unchanged ? 0 : 5;
}
