// Stage H2 acceptance: multi-cell aggregation and deterministic sparse
// compression (sections 7.3, 7.5, 7.6 and 15 H2).  A deterministic
// high-energy bulk state spanning many conversion cells (both u_parallel
// signs and non-zero u_perp, several x cells) is converted; the same global
// state must produce the same global created moments for 1, 2 and 5 ranks,
// and the six constraint moments must be independent of the aggregation-bin
// resolution while the number of particles follows the compression.
//
// Usage:
//   bulk_tail_multibin_test [--result <path>]
// Run plainly (1 rank) or with yhrun -n 2 / -n 5 (different --result paths).
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "grid.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

// Deterministic global bulk state: mass in every conversion cell of the
// local physical cells, depending only on the global cell index and the
// velocity cell, so 1/2/5-rank decompositions produce identical global data.
void fill_conversion_packet(Species& bulk, const SpatialGrid& grid,
                            const HybridVelocityPartition& partition)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
    for (int il = 0; il < nxl; ++il) {
        const int ixg = grid.ix_start + il;
        const double x_factor = 1.0 + static_cast<double>(ixg % 3) * 0.5;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                if (!partition.is_conversion(j, k)) continue;
                const double pattern =
                    1.0 + 0.1 * static_cast<double>((j * 7 + k * 13) % 5);
                bulk.f[idx3(ng + il, j, k)] =
                    1.0e19 * x_factor * pattern;
            }
        }
    }
}

struct RunResult {
    bool all_ok;
    bool all_complete;
    bool all_conservative;
    bool all_fidelity;
    bool all_finite;
    std::uint64_t global_particles_created;
    std::uint64_t global_fallback_count;
    double global_number_created;
    double global_px_created;
    double global_energy_created;
    double global_jx_created;
    double global_pixx_created;
    double global_piperp_created;
    double global_number_removed;
    double global_px_removed;
    double global_energy_removed;
    double global_jx_removed;
    double global_pixx_removed;
    double global_piperp_removed;
    double global_number_scale;
    double global_px_scale;
    double global_energy_scale;
    double global_jx_scale;
    double global_pixx_scale;
    double global_piperp_scale;
    double max_residual_rel;
    RunResult()
        : all_ok(false), all_complete(false), all_conservative(false),
          all_fidelity(false), all_finite(false), global_particles_created(0),
          global_fallback_count(0), global_number_created(0.0),
          global_px_created(0.0), global_energy_created(0.0),
          global_jx_created(0.0), global_pixx_created(0.0),
          global_piperp_created(0.0), global_number_removed(0.0),
          global_px_removed(0.0), global_energy_removed(0.0),
          global_jx_removed(0.0), global_pixx_removed(0.0),
          global_piperp_removed(0.0), global_number_scale(0.0),
          global_px_scale(0.0), global_energy_scale(0.0),
          global_jx_scale(0.0), global_pixx_scale(0.0),
          global_piperp_scale(0.0), max_residual_rel(0.0)
    {}
};

RunResult run_conversion(int rank, const HybridVelocityPartition& p,
                         Species& bulk)
{
    RunResult r;
    const SpatialGrid& grid = *bulk.sgrid;
    Species trial_bulk = bulk;
    BackgroundTailPIC trial_tail;
    trial_tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d =
        converter.extract_after_substep(trial_bulk, trial_tail, grid, p, 1,
                                        ConversionLocation::AFTER_U_SUBSTEP,
                                        rank);
    const int flags[5] = {
        d.complete ? 1 : 0, d.conservative ? 1 : 0,
        d.fidelity_ok ? 1 : 0, d.finite ? 1 : 0, 1
    };
    int global_flags[5] = { 0, 0, 0, 0, 0 };
    MPI_Allreduce(flags, global_flags, 5, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    r.all_complete = global_flags[0] != 0;
    r.all_conservative = global_flags[1] != 0;
    r.all_fidelity = global_flags[2] != 0;
    r.all_finite = global_flags[3] != 0;
    r.all_ok = r.all_complete && r.all_conservative && r.all_fidelity &&
               r.all_finite;

    const double local_created[6] = {
        d.number_created, d.px_created, d.energy_created,
        d.jx_dx_created, d.pixx_dx_created, d.piperp_dx_created
    };
    const double local_removed[6] = {
        d.number_removed, d.px_removed, d.energy_removed,
        d.jx_dx_removed, d.pixx_dx_removed, d.piperp_dx_removed
    };
    const double local_scales[6] = {
        d.number_scale, d.px_scale, d.energy_scale,
        d.jx_scale, d.pixx_scale, d.piperp_scale
    };
    double global_created[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double global_removed[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_created, global_created, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_removed, global_removed, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double global_scales[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_scales, global_scales, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    r.global_number_created = global_created[0];
    r.global_px_created = global_created[1];
    r.global_energy_created = global_created[2];
    r.global_jx_created = global_created[3];
    r.global_pixx_created = global_created[4];
    r.global_piperp_created = global_created[5];
    r.global_number_removed = global_removed[0];
    r.global_px_removed = global_removed[1];
    r.global_energy_removed = global_removed[2];
    r.global_jx_removed = global_removed[3];
    r.global_pixx_removed = global_removed[4];
    r.global_piperp_removed = global_removed[5];
    r.global_number_scale = global_scales[0];
    r.global_px_scale = global_scales[1];
    r.global_energy_scale = global_scales[2];
    r.global_jx_scale = global_scales[3];
    r.global_pixx_scale = global_scales[4];
    r.global_piperp_scale = global_scales[5];
    const double global_created_arr[6] = {
        r.global_number_created, r.global_px_created,
        r.global_energy_created, r.global_jx_created,
        r.global_pixx_created, r.global_piperp_created
    };
    const double global_scales_arr[6] = {
        r.global_number_scale, r.global_px_scale,
        r.global_energy_scale, r.global_jx_scale,
        r.global_pixx_scale, r.global_piperp_scale
    };
    r.max_residual_rel = 0.0;
    for (int i = 0; i < 6; ++i) {
        r.max_residual_rel = std::max(
            r.max_residual_rel,
            std::fabs(global_created_arr[i] - global_removed[i]) /
                std::max(1.0, global_scales_arr[i]));
    }

    std::uint64_t local_particles = d.particles_created;
    std::uint64_t local_fallback = d.compression_fallback_count;
    MPI_Allreduce(&local_particles, &r.global_particles_created, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_fallback, &r.global_fallback_count, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return r;
}

struct Metrics {
    bool all_ok;
    double max_residual_rel;
    double moment_bin_convergence_max;
    std::uint64_t particles_fine_bins;
    std::uint64_t particles_coarse_bins;
    std::uint64_t fallback_fine_bins;
    std::uint64_t fallback_coarse_bins;
    Metrics()
        : all_ok(false), max_residual_rel(0.0),
          moment_bin_convergence_max(0.0), particles_fine_bins(0),
          particles_coarse_bins(0), fallback_fine_bins(0),
          fallback_coarse_bins(0)
    {}
};

Metrics run_case(int rank, int size)
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, 12, 1.2 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);

    // Fine aggregation bins.
    HybridVelocityPartition p_fine;
    p_fine.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    fill_conversion_packet(bulk, grid, p_fine);
    const RunResult fine = run_conversion(rank, p_fine, bulk);
    m.all_ok = fine.all_ok;
    m.particles_fine_bins = fine.global_particles_created;
    m.fallback_fine_bins = fine.global_fallback_count;

    m.max_residual_rel = fine.max_residual_rel;

    // Coarse bins: the six constraint moments must be bin-resolution
    // independent (section 7.3: the spectrum converges; the moments are
    // exactly preserved), while the particle count is bounded by the
    // compression.
    HybridVelocityPartition p_coarse;
    p_coarse.init(bulk.cgrid, 6.0, 1.0, 1, 1);
    fill_conversion_packet(bulk, grid, p_coarse);
    const RunResult coarse = run_conversion(rank, p_coarse, bulk);
    m.particles_coarse_bins = coarse.global_particles_created;
    m.fallback_coarse_bins = coarse.global_fallback_count;
    m.all_ok = m.all_ok && coarse.all_ok;

    const double fine_moments[6] = {
        fine.global_number_created, fine.global_px_created,
        fine.global_energy_created, fine.global_jx_created,
        fine.global_pixx_created, fine.global_piperp_created
    };
    const double coarse_moments[6] = {
        coarse.global_number_created, coarse.global_px_created,
        coarse.global_energy_created, coarse.global_jx_created,
        coarse.global_pixx_created, coarse.global_piperp_created
    };
    double max_diff = 0.0;
    for (int i = 0; i < 6; ++i) {
        max_diff = std::max(
            max_diff,
            std::fabs(fine_moments[i] - coarse_moments[i]) /
                std::max(1.0, std::fabs(fine_moments[i])));
    }
    m.moment_bin_convergence_max = max_diff;

    // The coarse groups exercise the <=7-quartet compression: per spatial
    // cell, sign and threshold-aware energy edge bin there are at most 7
    // quartet supports (section 7.11.4 branch B: the production energy
    // grouping uses the explicit conversion_energy_edges, so the "coarse"
    // partition with energy_bins=1 still resolves the threshold fine bins).
    const std::uint64_t coarse_energy_bins = static_cast<std::uint64_t>(
        std::max(1, static_cast<int>(
                        p_coarse.conversion_energy_edges.size()) -
                        1));
    const std::uint64_t max_expected_particles =
        static_cast<std::uint64_t>(12) * 2 * coarse_energy_bins * 7 * 4;
    if (coarse.global_particles_created > max_expected_particles) {
        m.all_ok = false;
    }
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "all_ok=" << (m.all_ok ? 1 : 0) << "\n";
    out << "max_residual_rel=" << m.max_residual_rel << "\n";
    out << "moment_bin_convergence_max=" << m.moment_bin_convergence_max
        << "\n";
    out << "particles_fine_bins=" << m.particles_fine_bins << "\n";
    out << "particles_coarse_bins=" << m.particles_coarse_bins << "\n";
    out << "fallback_fine_bins=" << m.fallback_fine_bins << "\n";
    out << "fallback_coarse_bins=" << m.fallback_coarse_bins << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 5) {
        if (rank == 0) {
            std::cerr << "bulk_tail_multibin_test supports 1, 2 or 5 ranks.\n";
        }
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: bulk_tail_multibin_test [--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case(rank, size);
    bool pass = ok && m.all_ok &&
                m.max_residual_rel <= 1.0e-10 &&
                m.moment_bin_convergence_max <= 1.0e-9;
    int pass_all = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &pass_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    pass = pass_all != 0;

    if (rank == 0) {
        if (!write_result_file(args.result_path, m, pass)) pass = false;
        std::cout << std::setprecision(17);
        std::cout << "all_ok=" << (m.all_ok ? 1 : 0)
                  << " max_residual_rel=" << m.max_residual_rel
                  << " moment_bin_convergence_max="
                  << m.moment_bin_convergence_max
                  << " particles_fine_bins=" << m.particles_fine_bins
                  << " particles_coarse_bins=" << m.particles_coarse_bins
                  << " fallback_fine_bins=" << m.fallback_fine_bins
                  << " fallback_coarse_bins=" << m.fallback_coarse_bins
                  << "\n";
    }
    // Propagate the rank-0 result-file write outcome so every rank returns
    // the same exit code (section 7.11.3: no partial-rank exit on any
    // failure, including I/O).
    int final_pass = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &final_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    pass = final_pass != 0;
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
