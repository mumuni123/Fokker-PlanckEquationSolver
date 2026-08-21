#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "grid.h"
#include "tail_subcell_quadrature.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
struct Args {
    std::string test_case;
    std::string result;
    Args() : test_case("near-axis-narrow") {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) args.test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) args.result = argv[++i];
        else return false;
    }
    return args.test_case == "near-axis-narrow";
}

int bin(const std::vector<double>& edges, double e)
{
    if (e < edges.front() || e > edges.back()) return -1;
    if (e == edges.back()) return static_cast<int>(edges.size()) - 2;
    for (size_t b = 0; b + 1 < edges.size(); ++b)
        if (e >= edges[b] && e < edges[b + 1]) return static_cast<int>(b);
    return -1;
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks = 1; MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Args args;
    if (ranks != 1 || !parse_args(argc, argv, args)) {
        std::cerr << "usage: bulk_tail_subcell_loading_test --case near-axis-narrow "
                     "[--result path]\n";
        MPI_Finalize(); return 2;
    }
    SpatialGrid grid; grid.init_with_domain(0, 1, 4, 4.0 * Param::dx);
    Species bulk;
    bulk.init("subcell_loading", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
              Const::me, Param::dens, Param::temperature_e, false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    int iv = -1, imu = -1;
    const double threshold_window_hi =
        partition.min_conversion_energy + 0.2e6 * Const::eV;
    for (int k = 0; k < std::min(8, Param::Nmu) && iv < 0; ++k) {
        for (int j = 0; j < Param::Nv; ++j) {
            if (!partition.is_conversion(j, k)) continue;
            const std::vector<TailSubcellNode> nodes =
                TailSubcellQuadrature::nodes(bulk.cgrid, j, k);
            bool all_supported = !nodes.empty();
            bool in_threshold_window = false;
            for (size_t q = 0; q < nodes.size(); ++q) {
                all_supported = all_supported &&
                    nodes[q].kinetic_energy >= partition.min_conversion_energy;
                in_threshold_window = in_threshold_window ||
                    (nodes[q].kinetic_energy <= threshold_window_hi);
            }
            if (all_supported && in_threshold_window) {
                iv = j; imu = k; break;
            }
        }
    }
    bool pass = iv >= 0;
    BulkTailConversionDiagnostics d;
    BackgroundTailPIC tail;
    std::vector<double> reference;
    std::vector<double> loaded;
    if (pass) {
        const double m0 = 1.0e20;
        bulk.f[idx3(grid.nghost, iv, imu)] = m0;
        bulk.compute_moments();
        tail.init(grid);
        BulkTailConverter converter;
        converter.set_subcell_loading_enabled(true);
        d = converter.extract_after_substep(
            bulk, tail, grid, partition, 1,
            ConversionLocation::AFTER_U_SUBSTEP, 0);
        const std::vector<TailSubcellNode> nodes =
            TailSubcellQuadrature::nodes(bulk.cgrid, iv, imu);
        const double widths_mev[] = { 0.05, 0.1, 0.2 };
        std::vector<double> spectrum_l1_rel(3, 0.0);
        std::ostringstream spectrum_rows;
        for (int width_index = 0; width_index < 3; ++width_index) {
            const double lo = 5.6e6 * Const::eV;
            const double hi = 6.4e6 * Const::eV;
            const double width = widths_mev[width_index] * 1.0e6 * Const::eV;
            std::vector<double> edges;
            for (double e = lo; e < hi; e += width) edges.push_back(e);
            edges.push_back(hi);
            reference.assign(edges.size() - 1, 0.0);
            loaded.assign(edges.size() - 1, 0.0);
            for (size_t q = 0; q < nodes.size(); ++q) {
                const int b = bin(edges, nodes[q].kinetic_energy);
                if (b >= 0) reference[static_cast<size_t>(b)] +=
                    m0 * nodes[q].mass_fraction;
            }
            for (size_t p = 0; p < tail.particles.size(); ++p) {
                const BackgroundTailParticle& particle = tail.particles[p];
                const double gamma = std::sqrt(1.0 + particle.ux * particle.ux +
                    particle.uy * particle.uy + particle.uz * particle.uz);
                const int b = bin(edges, Const::me * Const::c * Const::c *
                                   (gamma - 1.0));
                if (b >= 0) loaded[static_cast<size_t>(b)] += particle.weight;
            }
            double l1 = 0.0, scale = 0.0;
            for (size_t b = 0; b < reference.size(); ++b) {
                l1 += std::fabs(loaded[b] - reference[b]);
                scale += std::fabs(reference[b]);
            }
            spectrum_l1_rel[static_cast<size_t>(width_index)] =
                l1 / std::max(1.0, scale);
            spectrum_rows << "bin_width_mev=" << widths_mev[width_index]
                          << " spectrum_L1_rel="
                          << spectrum_l1_rel[static_cast<size_t>(width_index)]
                          << "\n";
        }
        const double spectrum_l1_rel_max = std::max(spectrum_l1_rel[0],
            std::max(spectrum_l1_rel[1], spectrum_l1_rel[2]));
        pass = d.finite && d.complete && d.conservative && d.fidelity_ok &&
               d.subcell_fallback_count == 0 && d.subcell_cells_loaded > 0 &&
               spectrum_l1_rel_max <= 0.05;
        std::ostringstream out;
        out << std::setprecision(17)
            << "case=near-axis-narrow iv=" << iv << " imu=" << imu
            << " subcell_cells_loaded=" << d.subcell_cells_loaded
            << " subcell_support_count=" << d.subcell_support_count
            << " subcell_fallback_count=" << d.subcell_fallback_count
            << " N_residual=" << d.number_residual_rel
            << " Px_residual=" << d.px_residual_rel
            << " Jx_residual=" << d.jx_residual_rel
            << " K_residual=" << d.energy_residual_rel
            << " Pixx_residual=" << d.pixx_residual_rel
            << " Piperp_residual=" << d.piperp_residual_rel
            << " spectrum_L1_rel_max=" << spectrum_l1_rel_max << "\n"
            << spectrum_rows.str();
        for (size_t q = 0; q < d.subcell_fallbacks.size(); ++q) {
            const BulkTailConversionDiagnostics::SubcellFallbackRecord& f =
                d.subcell_fallbacks[q];
            out << "fallback ix=" << f.ix_global << " iv=" << f.iv
                << " imu=" << f.imu << " reason=" << f.reason
                << " particles=" << f.fallback_particles
                << " target_N=" << f.number_target
                << " target_Px=" << f.px_target
                << " target_Jx=" << f.jx_target
                << " target_K=" << f.energy_target
                << " target_Pixx=" << f.pixx_target
                << " target_Piperp=" << f.piperp_target << "\n";
        }
        out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
        std::cout << out.str();
        if (!args.result.empty()) { std::ofstream file(args.result.c_str(), std::ios::trunc); file << out.str(); }
    }
    if (!pass && iv < 0) std::cout << "status=FAIL no_near_axis_straddling_cell\n";
    MPI_Finalize(); return pass ? 0 : 1;
}
