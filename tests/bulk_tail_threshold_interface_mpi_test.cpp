// Section 7.11.3 (second step): MPI rank-consistency test for the
// bulk-to-tail threshold interface loading policies.  The same global
// analytic state is filled deterministically on every rank, converted with
// the selected policy, and the global six moments, per-bin spectrum and
// particle count are reduced and compared against the analytic reference.
// No random numbers participate.  Any rank failure is reduced collectively
// so every rank returns nonzero; only rank 0 writes the .result.
//
// Usage:
//   bulk_tail_threshold_interface_mpi_test
//     --case all --policy <golden|threshold-aware>
//     --bin-width-mev <value> [--result <path>]
// The last stdout line on every rank is "status=PASS" or "status=FAIL".

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
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string test_case;
    std::string policy;
    double bin_width_mev;
    std::string result_path;
    Args() : test_case("all"), policy("threshold-aware"), bin_width_mev(0.1) {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--policy") {
            if (i + 1 >= argc) return false;
            args.policy = argv[++i];
        } else if (arg == "--bin-width-mev") {
            if (i + 1 >= argc) return false;
            args.bin_width_mev = std::atof(argv[++i]);
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return (args.test_case == "all" || args.test_case == "smooth-exp" ||
            args.test_case == "broad-gaussian" ||
            args.test_case == "anisotropic-drift" ||
            args.test_case == "near-axis-narrow") &&
           (args.policy == "golden" ||
            args.policy == "threshold-aware") &&
           args.bin_width_mev > 0.0;
}

BulkTailLoadingPolicy policy_from_name(const std::string& name)
{
    return (name == "golden")
               ? BulkTailLoadingPolicy::GOLDEN_QUARTETS_NO_COMPRESSION
               : BulkTailLoadingPolicy::THRESHOLD_AWARE_COMPRESSION;
}

std::vector<double> make_edges(double ke_min, double ke_max, double bw_mev)
{
    std::vector<double> edges;
    const double bw = bw_mev * 1.0e6 * Const::eV;
    double e = ke_min;
    while (e < ke_max) {
        edges.push_back(e);
        e += bw;
    }
    edges.push_back(ke_max);
    return edges;
}

int bin_index(const std::vector<double>& edges, double ke)
{
    const size_t b = static_cast<size_t>(std::upper_bound(
        edges.begin(), edges.end(), ke) - edges.begin());
    return static_cast<int>(b) - 1;
}

void fill_case(Species& bulk, const SpatialGrid& grid,
               const HybridVelocityPartition& partition, int case_index)
{
    std::vector<double> shape(static_cast<size_t>(Param::Nvmu), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double ke = partition.kinetic_energy[idx2(j, k)];
            const double vol = bulk.cgrid.cell_phase_volume(j, k);
            const double upar = bulk.cgrid.upar_cells[j];
            const double uperp = bulk.cgrid.uperp_cells[k];
            double g = 0.0;
            if (case_index == 0) {
                g = std::exp(-ke / (1.0e6 * Const::eV));
            } else if (case_index == 1) {
                const double d =
                    (ke - 6.0e6 * Const::eV) / (1.2e6 * Const::eV);
                g = std::exp(-0.5 * d * d);
            } else if (case_index == 2) {
                const double r2 = (upar - 12.0) * (upar - 12.0) +
                                  uperp * uperp;
                g = std::exp(-r2 / 8.0);
            } else {
                const double ek = (ke - 6.0e6 * Const::eV) /
                                  (0.08e6 * Const::eV);
                const double uk = uperp / 0.18;
                g = std::exp(-0.5 * ek * ek - 0.5 * uk * uk) *
                    (1.0 + 0.02 * upar);
            }
            shape[idx2(j, k)] = g * vol;
        }
    }
    double conv_sum = 0.0;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            if (partition.is_conversion(j, k)) {
                conv_sum += shape[idx2(j, k)];
            }
        }
    }
    const double scale =
        (conv_sum > 0.0) ? 1.0e20 / conv_sum : 0.0;
    const int ng = grid.nghost;
    for (int il = 0; il < grid.nx_local; ++il) {
        const size_t xbase = static_cast<size_t>(ng + il) * Param::Nvmu;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                bulk.f[xbase + idx2(j, k)] =
                    scale * shape[idx2(j, k)];
            }
        }
    }
}

struct GlobalResult {
    bool ok;
    double created[6];
    double removed[6];
    double scale[6];
    std::uint64_t particles_created;
    std::uint64_t fallback_count;
    double spectrum_l1_rel;
    double spectrum_linf_rel;
    std::vector<double> global_tail_spectrum;
    std::vector<double> global_ref_spectrum;
};

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    Args args;
    if (!parse_args(argc, argv, args)) {
        if (mpi_rank == 0) {
            std::cerr << "usage: bulk_tail_threshold_interface_mpi_test "
                         "--case all --policy <golden|threshold-aware> "
                         "--bin-width-mev <value> [--result <path>]\n";
        }
        MPI_Finalize();
        return 2;
    }

    // The global analytic state is the same for any rank count: 100 spatial
    // cells, evenly divisible by 1/2/5.
    SpatialGrid grid;
    grid.init_with_domain(mpi_rank, mpi_size, 100, 100 * Param::dx);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    const double ke_min = partition.min_conversion_energy;
    const double ke_max = partition.max_conversion_energy;
    const std::vector<double> edges =
        make_edges(ke_min, ke_max, args.bin_width_mev);
    const BulkTailLoadingPolicy policy = policy_from_name(args.policy);

    bool global_pass = true;
    std::ostringstream report;
    report << std::setprecision(10);
    report << "threshold_interface_mpi policy=" << args.policy
           << " bin_width_mev=" << args.bin_width_mev
           << " ranks=" << mpi_size << "\n";

    const char* case_names[4] = {
        "smooth-exp", "broad-gaussian", "anisotropic-drift", "near-axis-narrow"
    };
    int executed_cases = 0;
    for (int ci = 0; ci < 4; ++ci) {
        if (args.test_case != "all" && args.test_case != case_names[ci]) continue;
        ++executed_cases;
        Species ref_bulk = bulk;
        fill_case(ref_bulk, grid, partition, ci);

        // Local reference spectrum and removed moments (the converter's
        // extraction is the only source of truth; the reference spectrum is
        // the conversion-cell mass).
        std::vector<double> ref_spectrum(edges.size() - 1, 0.0);
        const std::vector<double>& edge_edges =
            partition.conversion_energy_edges;
        std::vector<double> ref_edge_spectrum(edge_edges.size() - 1, 0.0);
        for (int il = 0; il < grid.nx_local; ++il) {
            const size_t xbase =
                static_cast<size_t>(grid.nghost + il) * Param::Nvmu;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    if (!partition.is_conversion(j, k)) continue;
                    const double mass = ref_bulk.f[xbase + idx2(j, k)];
                    if (!(mass > 0.0)) continue;
                    const int b = bin_index(
                        edges, partition.kinetic_energy[idx2(j, k)]);
                    if (b >= 0 &&
                        b < static_cast<int>(ref_spectrum.size())) {
                        ref_spectrum[static_cast<size_t>(b)] += mass;
                    }
                    const int eb = bin_index(
                        edge_edges, partition.kinetic_energy[idx2(j, k)]);
                    if (eb >= 0 &&
                        eb < static_cast<int>(ref_edge_spectrum.size())) {
                        ref_edge_spectrum[static_cast<size_t>(eb)] += mass;
                    }
                }
            }
        }

        Species trial = ref_bulk;
        BackgroundTailPIC tail;
        tail.init(grid);
        BulkTailConverter converter;
        converter.set_loading_policy(policy);
        const BulkTailConversionDiagnostics d =
            converter.extract_after_substep(
                trial, tail, grid, partition, 1,
                ConversionLocation::AFTER_U_SUBSTEP, mpi_rank);

        std::vector<double> tail_spectrum(edges.size() - 1, 0.0);
        std::vector<double> tail_edge_spectrum(edge_edges.size() - 1, 0.0);
        bool negative_weight = false;
        for (size_t p = 0; p < tail.particles.size(); ++p) {
            const BackgroundTailParticle& q = tail.particles[p];
            if (!(q.weight > 0.0)) negative_weight = true;
            const double gamma = std::sqrt(
                1.0 + q.ux * q.ux + q.uy * q.uy + q.uz * q.uz);
            const double ke = Const::me * Const::c * Const::c * (gamma - 1.0);
            const int b = bin_index(edges, ke);
            if (b >= 0 && b < static_cast<int>(tail_spectrum.size())) {
                tail_spectrum[static_cast<size_t>(b)] += q.weight;
            }
            const int eb = bin_index(edge_edges, ke);
            if (eb >= 0 &&
                eb < static_cast<int>(tail_edge_spectrum.size())) {
                tail_edge_spectrum[static_cast<size_t>(eb)] += q.weight;
            }
        }
        const bool local_ok =
            d.finite && d.conservative && d.complete && !negative_weight;

        // Collective reductions (identical order on every rank).
        const double created[6] = {
            d.number_created, d.px_created, d.jx_dx_created,
            d.energy_created, d.pixx_dx_created, d.piperp_dx_created
        };
        const double removed[6] = {
            d.number_removed, d.px_removed, d.jx_dx_removed,
            d.energy_removed, d.pixx_dx_removed, d.piperp_dx_removed
        };
        const double scale[6] = {
            d.number_scale, d.px_scale, d.energy_scale,
            d.jx_scale, d.pixx_scale, d.piperp_scale
        };
        double g_created[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double g_removed[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double g_scale[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        MPI_Allreduce(created, g_created, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(removed, g_removed, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(scale, g_scale, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        const size_t nbins = edges.size() - 1;
        std::vector<double> g_ref(nbins, 0.0);
        std::vector<double> g_tail(nbins, 0.0);
        MPI_Allreduce(ref_spectrum.data(), g_ref.data(),
                      static_cast<int>(nbins), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(tail_spectrum.data(), g_tail.data(),
                      static_cast<int>(nbins), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        const size_t nedge_bins = edge_edges.size() - 1;
        std::vector<double> g_ref_edge(nedge_bins, 0.0);
        std::vector<double> g_tail_edge(nedge_bins, 0.0);
        MPI_Allreduce(ref_edge_spectrum.data(), g_ref_edge.data(),
                      static_cast<int>(nedge_bins), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(tail_edge_spectrum.data(), g_tail_edge.data(),
                      static_cast<int>(nedge_bins), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        unsigned long long local_particles = d.particles_created;
        unsigned long long g_particles = 0;
        unsigned long long local_fallback = d.compression_fallback_count;
        unsigned long long g_fallback = 0;
        MPI_Allreduce(&local_particles, &g_particles, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_fallback, &g_fallback, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        int local_fail = local_ok ? 0 : 1;
        int global_fail = 0;
        MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_LOR,
                      MPI_COMM_WORLD);

        // Global residuals (created vs removed, normalised by the global
        // L1 scale) and spectrum fidelity vs the analytic reference.
        double residual_max = 0.0;
        for (int m = 0; m < 6; ++m) {
            residual_max = std::max(
                residual_max,
                std::fabs(g_created[m] - g_removed[m]) /
                    std::max(1.0, g_scale[m]));
        }
        double ref_denom = 0.0;
        double l1_diff = 0.0;
        double ref_max = 0.0;
        double linf_diff = 0.0;
        for (size_t i = 0; i < nbins; ++i) {
            ref_denom += std::fabs(g_ref[i]);
            l1_diff += std::fabs(g_tail[i] - g_ref[i]);
            ref_max = std::max(ref_max, std::fabs(g_ref[i]));
            linf_diff = std::max(linf_diff,
                                 std::fabs(g_tail[i] - g_ref[i]));
        }
        const double l1_rel = l1_diff / std::max(1.0e-300, ref_denom);
        const double linf_rel = linf_diff / std::max(1.0e-300, ref_max);
        double edge_ref_denom = 0.0;
        double edge_l1_diff = 0.0;
        for (size_t i = 0; i < nedge_bins; ++i) {
            edge_ref_denom += std::fabs(g_ref_edge[i]);
            edge_l1_diff += std::fabs(g_tail_edge[i] - g_ref_edge[i]);
        }
        const double edge_l1_rel =
            edge_l1_diff / std::max(1.0e-300, edge_ref_denom);
        const bool fine_spectrum_required =
            args.policy == "golden" || ci == 3;
        const bool fine_spectrum_pass =
            !fine_spectrum_required || l1_rel <= 1.0e-10;
        const bool edge_spectrum_pass =
            args.policy != "threshold-aware" || edge_l1_rel <= 1.0e-10;
        const bool fidelity_warning =
            args.policy == "threshold-aware" && !fine_spectrum_required &&
            l1_rel > 1.0e-10;
        const bool case_pass =
            global_fail == 0 &&
            residual_max <= 1.0e-10 &&
            fine_spectrum_pass && edge_spectrum_pass &&
            ref_denom > 0.0;
        if (!case_pass) global_pass = false;

        report << "case=" << case_names[ci]
               << " global_residual_max=" << residual_max
               << " spectrum_L1_rel=" << l1_rel
               << " spectrum_Linf_rel=" << linf_rel
               << " edge_spectrum_L1_rel=" << edge_l1_rel
               << " fine_spectrum_required="
               << (fine_spectrum_required ? 1 : 0)
               << " fine_spectrum_pass=" << (fine_spectrum_pass ? 1 : 0)
               << " edge_spectrum_pass=" << (edge_spectrum_pass ? 1 : 0)
               << " fidelity_warning=" << (fidelity_warning ? 1 : 0)
               << " particles_created=" << g_particles
               << " fallback_count=" << g_fallback
               << " case_pass=" << (case_pass ? 1 : 0) << "\n";
    }
    if (executed_cases == 0) global_pass = false;

    // Collective status so all ranks agree on the exit code (section
    // 7.11.3: no partial-rank early exit).
    int local_status = global_pass ? 0 : 1;
    int global_status = 0;
    MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);
    global_pass = global_status == 0;

    if (mpi_rank == 0) {
        std::cout << report.str();
        if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::trunc);
            if (out) {
                out << report.str();
                out << "status=" << (global_pass ? "PASS" : "FAIL") << "\n";
            }
        }
        std::cout << "status=" << (global_pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return global_pass ? 0 : 1;
}
