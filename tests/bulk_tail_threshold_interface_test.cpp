// Section 7.11.3 (second step of the bulk-to-tail threshold interface
// audit): single-rank A/B/C loading-policy interface test.  Three analytic
// smooth, strictly non-negative cell-integrated distributions are built over
// multiple spatial cells and each is converted with the three loading
// policies (golden quartets without compression, current production
// compression, threshold-aware compression).  The pre-conversion discrete
// reference spectrum, the created-tail spectrum and the six moments are
// compared on one set of explicit energy edges.  No random numbers are used.
//
// Usage:
//   bulk_tail_threshold_interface_test
//     --case <smooth-exp|broad-gaussian|anisotropic-drift|near-axis-narrow|all>
//     --policy <golden|current|threshold-aware|all>
//     --bin-width-mev <value> [--result <path>]
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
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string test_case;
    std::string policy;
    double bin_width_mev;
    std::string result_path;
    Args() : test_case("all"), policy("all"), bin_width_mev(0.1) {}
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
    const bool valid_case =
        args.test_case == "all" || args.test_case == "smooth-exp" ||
        args.test_case == "broad-gaussian" ||
        args.test_case == "anisotropic-drift" ||
        args.test_case == "near-axis-narrow";
    const bool valid_policy =
        args.policy == "all" || args.policy == "golden" ||
        args.policy == "current" || args.policy == "threshold-aware";
    return valid_case && valid_policy && args.bin_width_mev > 0.0;
}

bool case_selected(const Args& args, const char* name)
{
    return args.test_case == "all" || args.test_case == name;
}

bool policy_selected(const Args& args, const char* name)
{
    return args.policy == "all" || args.policy == name;
}

// Explicit spectrum edges over [ke_min, ke_max] with a constant width
// (Joules); the last edge is clamped to ke_max so every conversion cell and
// every created particle falls into a bin (section 7.11.3: one shared set of
// explicit boundaries).
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

double spectrum_l1_rel(const std::vector<double>& a,
                       const std::vector<double>& b)
{
    double denom = 0.0;
    double diff = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        denom += std::fabs(b[i]);
        diff += std::fabs(a[i] - b[i]);
    }
    return diff / std::max(1.0e-300, denom);
}

double spectrum_linf_rel(const std::vector<double>& a,
                         const std::vector<double>& b)
{
    double ref_max = 0.0;
    double diff_max = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        ref_max = std::max(ref_max, std::fabs(b[i]));
        diff_max = std::max(diff_max, std::fabs(a[i] - b[i]));
    }
    return diff_max / std::max(1.0e-300, ref_max);
}

// Max |d^2 ln(dN/dK)/dK^2| over the bins with positive density in the
// threshold window [K_out, K_out + 0.8 MeV]: a large negative second
// difference flags an isolated artificial valley.
double log_curvature_max(const std::vector<double>& dndk, int first_bin,
                         int last_bin)
{
    double worst = 0.0;
    for (int i = first_bin + 1; i + 1 <= last_bin; ++i) {
        const double y0 = dndk[i - 1];
        const double y1 = dndk[i];
        const double y2 = dndk[i + 1];
        if (!(y0 > 0.0) || !(y1 > 0.0) || !(y2 > 0.0)) continue;
        const double d2 = std::log(y0) - 2.0 * std::log(y1) +
                          std::log(y2);
        worst = std::max(worst, std::fabs(d2));
    }
    return worst;
}

// The section 7.11.1 valley signature: the middle of the three threshold
// bins [K_out, K_out+0.2], [K_out+0.2, K_out+0.4], [K_out+0.4, K_out+0.6]
// relative to the maximum of its two neighbours.
double threshold_valley_ratio(const std::vector<double>& dndk, int first_bin)
{
    const double middle = dndk[first_bin + 1];
    const double neigh = std::max(dndk[first_bin], dndk[first_bin + 2]);
    return middle / std::max(1.0e-300, neigh);
}

// ---- analytic distributions (cell-integrated masses, m^-2 per cell) ----

// Common shape normalisation: scale so the conversion-region mass per cell
// equals the production-scale target, keeping every cell mass above the
// converter roundoff floor.
double normalize_conversion(Species& bulk, const SpatialGrid& grid,
                            const HybridVelocityPartition& partition,
                            const std::vector<double>& shape_per_slot)
{
    const int ng = grid.nghost;
    double conv_sum = 0.0;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            if (partition.is_conversion(j, k)) {
                conv_sum += shape_per_slot[static_cast<size_t>(j) *
                                           Param::Nmu + k];
            }
        }
    }
    if (!(conv_sum > 0.0)) return 0.0;
    const double target = 1.0e20;
    const double scale = target / conv_sum;
    for (int il = 0; il < grid.nx_local; ++il) {
        const size_t xbase = static_cast<size_t>(ng + il) * Param::Nvmu;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                bulk.f[xbase + idx2(j, k)] =
                    scale * shape_per_slot[static_cast<size_t>(j) *
                                           Param::Nmu + k];
            }
        }
    }
    return scale;
}

void fill_smooth_exp(Species& bulk, const SpatialGrid& grid,
                     const HybridVelocityPartition& partition)
{
    const double t_j = 1.0e6 * Const::eV;
    std::vector<double> shape(static_cast<size_t>(Param::Nvmu), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double ke = partition.kinetic_energy[idx2(j, k)];
            const double vol = bulk.cgrid.cell_phase_volume(j, k);
            shape[idx2(j, k)] = std::exp(-ke / t_j) * vol;
        }
    }
    normalize_conversion(bulk, grid, partition, shape);
}

void fill_broad_gaussian(Species& bulk, const SpatialGrid& grid,
                         const HybridVelocityPartition& partition)
{
    const double center = 6.0e6 * Const::eV;
    const double sigma = 1.2e6 * Const::eV;
    std::vector<double> shape(static_cast<size_t>(Param::Nvmu), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double ke = partition.kinetic_energy[idx2(j, k)];
            const double vol = bulk.cgrid.cell_phase_volume(j, k);
            const double d = (ke - center) / sigma;
            shape[idx2(j, k)] = std::exp(-0.5 * d * d) * vol;
        }
    }
    normalize_conversion(bulk, grid, partition, shape);
}

void fill_anisotropic_drift(Species& bulk, const SpatialGrid& grid,
                            const HybridVelocityPartition& partition)
{
    const double u_drift = 12.0;
    const double sigma_u = 2.0;
    std::vector<double> shape(static_cast<size_t>(Param::Nvmu), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double upar = bulk.cgrid.upar_cells[j];
            const double uperp = bulk.cgrid.uperp_cells[k];
            const double vol = bulk.cgrid.cell_phase_volume(j, k);
            const double r2 = (upar - u_drift) * (upar - u_drift) +
                              uperp * uperp;
            shape[idx2(j, k)] =
                std::exp(-r2 / (2.0 * sigma_u * sigma_u)) * vol;
        }
    }
    normalize_conversion(bulk, grid, partition, shape);
}

void fill_near_axis_narrow(Species& bulk, const SpatialGrid& grid,
                           const HybridVelocityPartition& partition)
{
    std::vector<double> shape(static_cast<size_t>(Param::Nvmu), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double upar = bulk.cgrid.upar_cells[j];
            const double uperp = bulk.cgrid.uperp_cells[k];
            const double ke = partition.kinetic_energy[idx2(j, k)];
            const double threshold = 6.0e6 * Const::eV;
            const double ek = (ke - threshold) / (0.08e6 * Const::eV);
            const double uk = uperp / 0.18;
            shape[idx2(j, k)] = std::exp(-0.5 * ek * ek - 0.5 * uk * uk) *
                                (1.0 + 0.02 * upar) *
                                bulk.cgrid.cell_phase_volume(j, k);
        }
    }
    normalize_conversion(bulk, grid, partition, shape);
}

struct RunMetrics {
    bool ok;
    double n_res, px_res, jx_res, k_res, pixx_res, piperp_res;
    double rho_l2, rho_linf;
    double spectrum_l1_rel, spectrum_linf_rel;
    std::vector<double> near_threshold_bin_relative;
    double log_curvature_max;
    double threshold_valley_ratio;
    std::uint64_t particles_created;
    std::uint64_t compression_fallback_count;
    // Edge-resolution spectra (section 7.11.4 branch B): mass and kinetic
    // energy per explicit conversion_energy_edges bin, used to verify that
    // each threshold fine bin's N and K do not cross groups (section 7.11.6
    // item 2).
    std::vector<double> edge_spectrum_n;
    std::vector<double> edge_spectrum_k;
    RunMetrics()
        : ok(false), n_res(0.0), px_res(0.0), jx_res(0.0), k_res(0.0),
          pixx_res(0.0), piperp_res(0.0), rho_l2(0.0), rho_linf(0.0),
          spectrum_l1_rel(0.0), spectrum_linf_rel(0.0),
          log_curvature_max(0.0), threshold_valley_ratio(0.0),
          particles_created(0), compression_fallback_count(0)
    {}
};

RunMetrics run_case_policy(Species bulk, const SpatialGrid& grid,
                           const HybridVelocityPartition& partition,
                           BulkTailLoadingPolicy policy,
                           const std::vector<double>& edges,
                           const std::vector<double>& ref_spectrum,
                           const std::vector<double>& edge_edges)
{
    RunMetrics m;
    BackgroundTailPIC tail;
    tail.init(grid);
    BulkTailConverter converter;
    converter.set_loading_policy(policy);
    const BulkTailConversionDiagnostics d =
        converter.extract_after_substep(
            bulk, tail, grid, partition, 1,
            ConversionLocation::AFTER_U_SUBSTEP, 0);
    m.n_res = d.number_residual_rel;
    m.px_res = d.px_residual_rel;
    m.jx_res = d.jx_residual_rel;
    m.k_res = d.energy_residual_rel;
    m.pixx_res = d.pixx_residual_rel;
    m.piperp_res = d.piperp_residual_rel;
    m.rho_l2 = d.rho_l2_before_after;
    m.rho_linf = d.rho_linf_before_after;
    m.particles_created = d.particles_created;
    m.compression_fallback_count = d.compression_fallback_count;

    // Created-tail spectrum on the same explicit edges (raw count and mass,
    // then dN/dK = mass per bin / bin width).
    std::vector<double> tail_spectrum(edges.size() - 1, 0.0);
    std::vector<long> tail_raw(edges.size() - 1, 0);
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
            tail_raw[static_cast<size_t>(b)] += 1;
        }
    }
    // dN/dK arrays (mass per MeV).
    std::vector<double> ref_dndk(edges.size() - 1, 0.0);
    std::vector<double> tail_dndk(edges.size() - 1, 0.0);
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        const double width_mev =
            (edges[i + 1] - edges[i]) / (1.0e6 * Const::eV);
        ref_dndk[i] = ref_spectrum[i] / std::max(width_mev, 1.0e-300);
        tail_dndk[i] = tail_spectrum[i] / std::max(width_mev, 1.0e-300);
    }
    m.spectrum_l1_rel = spectrum_l1_rel(tail_spectrum, ref_spectrum);
    m.spectrum_linf_rel = spectrum_linf_rel(tail_spectrum, ref_spectrum);

    // Near-threshold per-bin relative errors over [K_out, K_out+0.8 MeV].
    const double k_out = partition.min_conversion_energy;
    const int first_bin = bin_index(edges, k_out);
    const int last_bin = bin_index(edges, k_out + 0.8e6 * Const::eV);
    for (int b = first_bin; b <= last_bin; ++b) {
        if (b < 0 || b >= static_cast<int>(ref_spectrum.size())) continue;
        const double ref_max =
            std::max(ref_spectrum[static_cast<size_t>(b)], 1.0e-300);
        m.near_threshold_bin_relative.push_back(
            std::fabs(tail_spectrum[static_cast<size_t>(b)] -
                      ref_spectrum[static_cast<size_t>(b)]) / ref_max);
    }
    m.log_curvature_max = log_curvature_max(tail_dndk, first_bin, last_bin);
    m.threshold_valley_ratio = threshold_valley_ratio(tail_dndk, first_bin);
    m.ok = d.finite && d.conservative && d.complete && !negative_weight;

    // Edge-resolution spectra on partition.conversion_energy_edges.
    m.edge_spectrum_n.assign(edge_edges.size() - 1, 0.0);
    m.edge_spectrum_k.assign(edge_edges.size() - 1, 0.0);
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        const BackgroundTailParticle& q = tail.particles[p];
        const double gamma = std::sqrt(
            1.0 + q.ux * q.ux + q.uy * q.uy + q.uz * q.uz);
        const double ke = Const::me * Const::c * Const::c * (gamma - 1.0);
        const int b = bin_index(edge_edges, ke);
        if (b >= 0 &&
            b < static_cast<int>(m.edge_spectrum_n.size())) {
            m.edge_spectrum_n[static_cast<size_t>(b)] += q.weight;
            m.edge_spectrum_k[static_cast<size_t>(b)] += q.weight * ke;
        }
    }
    return m;
}

// Relative N (mass) and K (mass x KE) error of `b` against the golden
// reference `a`, maximised over edge bins with nonzero reference mass.
double edge_bin_max_rel(const std::vector<double>& a,
                        const std::vector<double>& b)
{
    double worst = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!(std::fabs(a[i]) > 0.0)) continue;
        worst = std::max(worst,
                         std::fabs(b[i] - a[i]) / std::fabs(a[i]));
    }
    return worst;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "bulk_tail_threshold_interface_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: bulk_tail_threshold_interface_test "
                     "--case <smooth-exp|broad-gaussian|anisotropic-drift|all> "
                     "--policy <golden|current|threshold-aware|all> "
                     "--bin-width-mev <value> [--result <path>]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 20, 20 * Param::dx);
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
    const std::vector<double>& edge_edges = partition.conversion_energy_edges;

    bool pass = true;
    int executed_cases = 0;
    std::ostringstream report;
    report << std::setprecision(10);
    report << "threshold_interface bin_width_mev=" << args.bin_width_mev
           << " convert_energy_mev=6\n";

    const char* case_names[4] = {
        "smooth-exp", "broad-gaussian", "anisotropic-drift", "near-axis-narrow"
    };
    const char* policy_names[3] = { "golden", "current", "threshold-aware" };
    const BulkTailLoadingPolicy policies[3] = {
        BulkTailLoadingPolicy::GOLDEN_QUARTETS_NO_COMPRESSION,
        BulkTailLoadingPolicy::CURRENT_PRODUCTION_COMPRESSION,
        BulkTailLoadingPolicy::THRESHOLD_AWARE_COMPRESSION
    };

    for (int ci = 0; ci < 4; ++ci) {
        if (!case_selected(args, case_names[ci])) continue;
        ++executed_cases;
        // Fresh analytic reference bulk.
        Species ref_bulk = bulk;
        if (ci == 0) fill_smooth_exp(ref_bulk, grid, partition);
        else if (ci == 1) fill_broad_gaussian(ref_bulk, grid, partition);
        else if (ci == 2) fill_anisotropic_drift(ref_bulk, grid, partition);
        else fill_near_axis_narrow(ref_bulk, grid, partition);

        // Pre-conversion discrete reference spectrum (conversion-cell mass).
        std::vector<double> ref_spectrum(edges.size() - 1, 0.0);
        std::vector<long> ref_raw(edges.size() - 1, 0);
        std::vector<double> ref_edge_n(edge_edges.size() - 1, 0.0);
        std::vector<double> ref_edge_k(edge_edges.size() - 1, 0.0);
        double ref_threshold_mass = 0.0;
        for (int il = 0; il < grid.nx_local; ++il) {
            const size_t xbase =
                static_cast<size_t>(grid.nghost + il) * Param::Nvmu;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    if (!partition.is_conversion(j, k)) continue;
                    const double mass =
                        ref_bulk.f[xbase + idx2(j, k)];
                    if (!(mass > 0.0)) continue;
                    const double ke = partition.kinetic_energy[idx2(j, k)];
                    const int b = bin_index(edges, ke);
                    if (b >= 0 &&
                        b < static_cast<int>(ref_spectrum.size())) {
                        ref_spectrum[static_cast<size_t>(b)] += mass;
                        ref_raw[static_cast<size_t>(b)] += 1;
                        if (ke <= ke_min + 0.8e6 * Const::eV) {
                            ref_threshold_mass += mass;
                        }
                    }
                    const int eb = bin_index(edge_edges, ke);
                    if (eb >= 0 &&
                        eb < static_cast<int>(ref_edge_n.size())) {
                        ref_edge_n[static_cast<size_t>(eb)] += mass;
                        ref_edge_k[static_cast<size_t>(eb)] += mass * ke;
                    }
                }
            }
        }
        report << "case=" << case_names[ci]
               << " reference_threshold_mass=" << ref_threshold_mass << "\n";

        std::vector<double> golden_edge_n;
        std::vector<double> golden_edge_k;
        for (int pi = 0; pi < 3; ++pi) {
            if (!policy_selected(args, policy_names[pi])) continue;
            const RunMetrics m = run_case_policy(
                ref_bulk, grid, partition, policies[pi], edges,
                ref_spectrum, edge_edges);
            if (pi == 0) {
                golden_edge_n = m.edge_spectrum_n;
                golden_edge_k = m.edge_spectrum_k;
            }
            const double edge_l1 = spectrum_l1_rel(
                m.edge_spectrum_n, ref_edge_n);
            const double edge_bin_n_rel =
                golden_edge_n.empty()
                    ? edge_bin_max_rel(ref_edge_n, m.edge_spectrum_n)
                    : edge_bin_max_rel(golden_edge_n, m.edge_spectrum_n);
            const double edge_bin_k_rel =
                golden_edge_k.empty()
                    ? edge_bin_max_rel(ref_edge_k, m.edge_spectrum_k)
                    : edge_bin_max_rel(golden_edge_k, m.edge_spectrum_k);
            report << "case=" << case_names[ci]
                   << " policy=" << policy_names[pi]
                   << " ok=" << (m.ok ? 1 : 0)
                   << " N=" << m.n_res << " Px=" << m.px_res
                   << " Jx=" << m.jx_res << " K=" << m.k_res
                   << " Pixx=" << m.pixx_res << " Piperp=" << m.piperp_res
                   << " rho_l2=" << m.rho_l2
                   << " rho_linf=" << m.rho_linf
                   << " spectrum_L1_rel=" << m.spectrum_l1_rel
                   << " spectrum_Linf_rel=" << m.spectrum_linf_rel
                   << " log_curvature_max=" << m.log_curvature_max
                   << " threshold_valley_ratio=" << m.threshold_valley_ratio
                   << " particles_created=" << m.particles_created
                   << " compression_fallback_count="
                   << m.compression_fallback_count
                   << " edge_spectrum_L1_rel=" << edge_l1
                   << " edge_bin_N_max_rel=" << edge_bin_n_rel
                   << " edge_bin_K_max_rel=" << edge_bin_k_rel << "\n";
            report << "case=" << case_names[ci]
                   << " policy=" << policy_names[pi]
                   << " near_threshold_bin_relative:";
            for (size_t i = 0; i < m.near_threshold_bin_relative.size();
                 ++i) {
                report << " " << m.near_threshold_bin_relative[i];
            }
            report << "\n";
            if (!m.ok) pass = false;
            // Golden reference must be exactly conservative (six moments at
            // or below the section 7.11.6 item-1 gate).
            if (pi == 0) {
                const double golden_max = std::max(
                    std::max(m.n_res, m.px_res),
                    std::max(m.jx_res, std::max(m.k_res,
                                                 std::max(m.pixx_res,
                                                          m.piperp_res))));
                if (!(golden_max <= 1.0e-10)) pass = false;
            }
            // Section 7.11.6 item 2: the threshold-aware policy must match
            // the golden reference (fallback: the pre-conversion reference,
            // which the golden reproduces to roundoff) per threshold fine
            // (edge) bin in N and K and in the full edge-resolution
            // spectrum.
            if (pi == 2) {
                const std::vector<double>& n_ref =
                    golden_edge_n.empty() ? ref_edge_n : golden_edge_n;
                const std::vector<double>& k_ref =
                    golden_edge_k.empty() ? ref_edge_k : golden_edge_k;
                const double n_rel =
                    edge_bin_max_rel(n_ref, m.edge_spectrum_n);
                const double k_rel =
                    edge_bin_max_rel(k_ref, m.edge_spectrum_k);
                if (!(edge_l1 <= 1.0e-10 && n_rel <= 1.0e-10 &&
                      k_rel <= 1.0e-10)) {
                    pass = false;
                }
            }
        }
        if (!(ref_threshold_mass > 0.0)) {
            report << "case=" << case_names[ci]
                   << " ERROR: no reference mass near the threshold\n";
            pass = false;
        }
    }
    if (executed_cases == 0) pass = false;

    std::cout << report.str();
    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::trunc);
        if (out) {
            out << report.str();
            out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
