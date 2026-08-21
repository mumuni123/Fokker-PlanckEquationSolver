#include "grid.h"
#include "species.h"
#include "tail_subcell_quadrature.h"

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
    double convert_energy_mev;
    std::vector<int> rings;
    std::vector<double> widths;
    std::string result;
    Args() : convert_energy_mev(6.0) {}
};

bool parse_list(const char* text, std::vector<double>& values)
{
    std::stringstream in(text); std::string value;
    while (std::getline(in, value, ',')) {
        if (value.empty()) return false;
        values.push_back(std::atof(value.c_str()));
    }
    return !values.empty();
}

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--convert-energy-mev" && i + 1 < argc) {
            args.convert_energy_mev = std::atof(argv[++i]);
        } else if (arg == "--uperp-rings" && i + 1 < argc) {
            std::vector<double> raw; if (!parse_list(argv[++i], raw)) return false;
            for (size_t r = 0; r < raw.size(); ++r) args.rings.push_back(static_cast<int>(raw[r]));
        } else if (arg == "--bin-widths-mev" && i + 1 < argc) {
            if (!parse_list(argv[++i], args.widths)) return false;
        } else if (arg == "--result" && i + 1 < argc) {
            args.result = argv[++i];
        } else return false;
    }
    if (args.rings.empty()) for (int k = 0; k < 8 && k < Param::Nmu; ++k) args.rings.push_back(k);
    if (args.widths.empty()) { args.widths.push_back(0.05); args.widths.push_back(0.1); args.widths.push_back(0.2); }
    return args.convert_energy_mev > 0.0;
}

std::vector<double> edges(double lo, double hi, double width_mev)
{
    std::vector<double> result; const double width = width_mev * 1.0e6 * Const::eV;
    for (double e = lo; e < hi; e += width) result.push_back(e);
    result.push_back(hi); return result;
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks = 1; MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Args args;
    if (ranks != 1 || !parse_args(argc, argv, args)) {
        std::cerr << "usage: bulk_tail_cell_volume_spectrum_test --convert-energy-mev 6 "
                     "--uperp-rings 0,1,2,3,4,5,6,7 --bin-widths-mev 0.05,0.1,0.2 [--result path]\n";
        MPI_Finalize(); return 2;
    }
    SpatialGrid spatial; spatial.init_with_domain(0, 1, 1, Param::dx);
    Species bulk;
    bulk.init("subcell_audit", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
              Const::me, Param::dens, Param::temperature_e, false, spatial);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, args.convert_energy_mev, 1.0, 4, 4);
    const double threshold = partition.min_conversion_energy;
    const double lo = std::max(0.0, threshold - 1.0e6 * Const::eV);
    const double hi = std::min(partition.max_conversion_energy,
                               threshold + 1.0e6 * Const::eV);
    bool pass = true;
    std::ostringstream report; report << std::setprecision(17);
    report << "cell_volume_spectrum convert_energy_mev=" << args.convert_energy_mev << "\n";
    for (size_t wi = 0; wi < args.widths.size(); ++wi) {
        const std::vector<double> bin_edges = edges(lo, hi, args.widths[wi]);
        TailSubcellSpectrum spectrum(bin_edges.size() - 1);
        std::vector<double> center(bin_edges.size() - 1, 0.0);
        double threshold_window_center_n = 0.0;
        double threshold_window_volume_n = 0.0;
        double selected_threshold_window_volume_n = 0.0;
        double selected_threshold_window_volume_px = 0.0;
        double selected_threshold_window_volume_k = 0.0;
        int selected_cell_count = 0;
        int selected_threshold_support_cell_count = 0;
        const double window_low_mev = args.convert_energy_mev;
        const double window_high_mev = args.convert_energy_mev + 0.2;
        for (size_t ri = 0; ri < args.rings.size(); ++ri) {
            const int k = args.rings[ri];
            if (k < 0 || k >= Param::Nmu) { pass = false; continue; }
            for (int j = 0; j < Param::Nv; ++j) {
                // Unit cell-constant distribution: mass equals the actual
                // cylindrical cell phase volume, exposing geometric support.
                const double mass = bulk.cgrid.cell_phase_volume(j, k);
                TailSubcellQuadrature::accumulate_cell(bulk.cgrid, j, k, mass,
                                                        bin_edges, threshold, spectrum);
                const int b = TailSubcellQuadrature::energy_bin(
                    bin_edges, bulk.cgrid.kinetic_energy[idx2(j, k)]);
                if (b >= 0) center[static_cast<size_t>(b)] += mass;
                const double centre_mev = bulk.cgrid.kinetic_energy[idx2(j, k)] /
                                          (1.0e6 * Const::eV);
                const bool selected =
                    partition.is_conversion_cell[idx2(j, k)] != 0;
                if (selected) ++selected_cell_count;
                if (centre_mev >= window_low_mev && centre_mev < window_high_mev)
                    threshold_window_center_n += mass;
                const std::vector<TailSubcellNode> nodes =
                    TailSubcellQuadrature::nodes(bulk.cgrid, j, k);
                double below_fraction = 0.0;
                double above_fraction = 0.0;
                double n = 0.0, px = 0.0, ke = 0.0;
                double selected_window_n_cell = 0.0;
                double selected_window_px_cell = 0.0;
                double selected_window_k_cell = 0.0;
                for (size_t q = 0; q < nodes.size(); ++q) {
                    const double qm = mass * nodes[q].mass_fraction;
                    if (nodes[q].kinetic_energy < threshold)
                        below_fraction += nodes[q].mass_fraction;
                    else
                        above_fraction += nodes[q].mass_fraction;
                    n += qm;
                    px += Const::me * Const::c * qm * nodes[q].upar;
                    ke += qm * nodes[q].kinetic_energy;
                    const double node_mev = nodes[q].kinetic_energy /
                                            (1.0e6 * Const::eV);
                    if (node_mev >= window_low_mev &&
                        node_mev < window_high_mev) {
                        threshold_window_volume_n += qm;
                        if (selected) {
                            selected_window_n_cell += qm;
                            selected_window_px_cell +=
                                Const::me * Const::c * qm * nodes[q].upar;
                            selected_window_k_cell += qm * nodes[q].kinetic_energy;
                        }
                    }
                }
                if (selected_window_n_cell > 0.0) {
                    ++selected_threshold_support_cell_count;
                    selected_threshold_window_volume_n += selected_window_n_cell;
                    selected_threshold_window_volume_px += selected_window_px_cell;
                    selected_threshold_window_volume_k += selected_window_k_cell;
                    if (wi == 0) {
                        report << "threshold_support_cell uperp_index=" << k
                               << " upar_index=" << j
                               << " upar_low=" << bulk.cgrid.upar_faces[j]
                               << " upar_high=" << bulk.cgrid.upar_faces[j + 1]
                               << " uperp_low=" << bulk.cgrid.uperp_faces[k]
                               << " uperp_high=" << bulk.cgrid.uperp_faces[k + 1]
                               << " N_window=" << selected_window_n_cell
                               << " Px_window=" << selected_window_px_cell
                               << " K_window=" << selected_window_k_cell << "\n";
                    }
                }
                if (wi == 0 && below_fraction > 0.0 && above_fraction > 0.0) {
                    report << "straddling_cell uperp_index=" << k
                           << " upar_index=" << j
                           << " upar_low=" << bulk.cgrid.upar_faces[j]
                           << " upar_high=" << bulk.cgrid.upar_faces[j + 1]
                           << " uperp_low=" << bulk.cgrid.uperp_faces[k]
                           << " uperp_high=" << bulk.cgrid.uperp_faces[k + 1]
                           << " below_volume_fraction=" << below_fraction
                           << " above_volume_fraction=" << above_fraction
                           << " N=" << n << " Px=" << px << " K=" << ke << "\n";
                }
            }
        }
        double histogram_n = 0.0;
        for (size_t b = 0; b < center.size(); ++b) {
            histogram_n += spectrum.number[b];
            report << "bin_width_mev=" << args.widths[wi]
                   << " bin_low_mev=" << bin_edges[b] / (1.0e6 * Const::eV)
                   << " bin_high_mev=" << bin_edges[b + 1] / (1.0e6 * Const::eV)
                   << " cell_center_histogram_N=" << center[b]
                   << " cell_volume_histogram_N=" << spectrum.number[b]
                   << " cell_volume_histogram_K=" << spectrum.energy[b] << "\n";
        }
        const double nrel = std::fabs(spectrum.number_residual) /
                            std::max(1.0, spectrum.input_number);
        report << "bin_width_mev=" << args.widths[wi]
               << " cell_volume_N_residual=" << spectrum.number_residual
                << " cell_center_to_volume_K_delta=" << spectrum.energy_residual
               << " straddling_cell_count=" << spectrum.straddling_cell_count
                << " straddling_cell_mass=" << spectrum.straddling_cell_mass
                << " threshold_window_6p0_6p2_center_N=" << threshold_window_center_n
                << " threshold_window_6p0_6p2_volume_N=" << threshold_window_volume_n
                << " selected_cell_count=" << selected_cell_count
                << " selected_threshold_support_cell_count="
                << selected_threshold_support_cell_count
                << " selected_threshold_window_6p0_6p2_volume_N="
                << selected_threshold_window_volume_n
                << " selected_threshold_window_6p0_6p2_volume_Px="
                << selected_threshold_window_volume_px
                << " selected_threshold_window_6p0_6p2_volume_K="
                << selected_threshold_window_volume_k
                << " represented_histogram_N=" << histogram_n << "\n";
        // A threshold may coincide with a velocity-cell face, in which case
        // no positive-volume cell straddles it.  The relevant feasibility
        // gate is whether cells already selected by the production
        // conversion mask contain subcell support in the threshold window.
        if (!(std::isfinite(nrel) && nrel <= 1.0e-12 &&
              selected_cell_count > 0 &&
              selected_threshold_support_cell_count > 0 &&
              std::isfinite(selected_threshold_window_volume_n) &&
              selected_threshold_window_volume_n > 0.0)) {
            pass = false;
        }
    }
    std::cout << report.str();
    if (!args.result.empty()) {
        std::ofstream out(args.result.c_str(), std::ios::trunc);
        out << report.str() << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize(); return pass ? 0 : 1;
}
