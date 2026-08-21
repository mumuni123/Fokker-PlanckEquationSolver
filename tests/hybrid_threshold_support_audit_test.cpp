// Section 7.11.2 (first step of the bulk-to-tail threshold interface
// audit): a read-only audit of the production velocity-grid support near
// the conversion threshold K_out.  The test directly constructs the
// production CylindricalVelocityGrid and HybridVelocityPartition (it never
// re-derives the grid formulas), reports, for K_out +/- 1 MeV and three
// diagnostic bin widths (0.05/0.1/0.2 MeV), the number of bulk-resolved and
// conversion cell centres per bin, the conversion centre energy range, the
// largest uncovered energy gap and the conversion phase-volume sum, and
// additionally dumps the sorted unique conversion cell-centre energies with
// their adjacent gaps.  The test advances no equation, creates no particle
// and changes no production parameter.
//
// Usage:
//   hybrid_threshold_support_audit_test
//     --convert-energy-mev <value> --uperp-rings 0,1,...
//     --bin-widths-mev 0.05,0.1,0.2 [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

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
    double convert_energy_mev;
    std::vector<int> uperp_rings;
    std::vector<double> bin_widths_mev;
    std::string result_path;
    Args() : convert_energy_mev(6.0) {}
};

template <typename T>
bool parse_csv(const char* text, std::vector<T>& values)
{
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) return false;
        values.push_back(static_cast<T>(std::atof(item.c_str())));
    }
    return !values.empty();
}

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--convert-energy-mev") {
            if (i + 1 >= argc) return false;
            args.convert_energy_mev = std::atof(argv[++i]);
        } else if (arg == "--uperp-rings") {
            if (i + 1 >= argc || !parse_csv(argv[++i], args.uperp_rings)) return false;
        } else if (arg == "--bin-widths-mev") {
            if (i + 1 >= argc || !parse_csv(argv[++i], args.bin_widths_mev)) return false;
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    if (args.uperp_rings.empty()) {
        for (int k = 0; k < 8 && k < Param::Nmu; ++k) args.uperp_rings.push_back(k);
    }
    if (args.bin_widths_mev.empty()) {
        args.bin_widths_mev.push_back(0.05);
        args.bin_widths_mev.push_back(0.1);
        args.bin_widths_mev.push_back(0.2);
    }
    return args.convert_energy_mev > 0.0;
}

struct BinAudit {
    int bulk_center_count;
    int conversion_center_count;
    double min_center_energy_mev;
    double max_center_energy_mev;
    double max_uncovered_energy_gap_mev;
    double phase_volume_sum;
    BinAudit()
        : bulk_center_count(0), conversion_center_count(0),
          min_center_energy_mev(-1.0), max_center_energy_mev(-1.0),
          max_uncovered_energy_gap_mev(-1.0), phase_volume_sum(0.0)
    {}
};

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "hybrid_threshold_support_audit_test: single-rank "
                     "read-only audit\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: hybrid_threshold_support_audit_test "
                     "--convert-energy-mev <value> [--uperp-rings 0,1] "
                     "[--bin-widths-mev 0.05,0.1,0.2] [--result <path>]\n";
        MPI_Finalize();
        return 2;
    }

    // Production velocity grid and partition: constructed directly, never
    // re-derived (section 7.11.2).  The spatial grid is only needed to host
    // the Species that builds the production cylindrical velocity grid.
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 8 * Param::dx);
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    HybridVelocityPartition partition;
    partition.init(electrons.cgrid, args.convert_energy_mev, 1.0, 4, 4);

    const double k_out = partition.min_conversion_energy;
    const double k_out_mev = k_out / (1.0e6 * Const::eV);
    const double lo_mev = k_out_mev - 1.0;
    const double hi_mev = k_out_mev + 1.0;

    bool pass = true;
    std::ostringstream report;
    report << std::setprecision(10);
    report << "audit convert_energy_mev=" << args.convert_energy_mev
           << " window_mev=" << lo_mev << ".." << hi_mev << "\n";

    // Full sorted unique conversion cell-centre energies and gaps.
    std::vector<double> conv_centres_mev;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            if (!partition.is_conversion(j, k)) continue;
            const double ke_mev =
                partition.kinetic_energy[idx2(j, k)] /
                (1.0e6 * Const::eV);
            if (ke_mev < lo_mev || ke_mev > hi_mev) continue;
            conv_centres_mev.push_back(ke_mev);
        }
    }
    std::sort(conv_centres_mev.begin(), conv_centres_mev.end());
    conv_centres_mev.erase(
        std::unique(conv_centres_mev.begin(), conv_centres_mev.end()),
        conv_centres_mev.end());
    report << "conversion_center_count_window=" << conv_centres_mev.size()
           << "\n";
    report << "conversion_centres_mev:";
    for (size_t i = 0; i < conv_centres_mev.size(); ++i) {
        report << " " << conv_centres_mev[i];
    }
    report << "\n";
    report << "adjacent_gaps_mev:";
    for (size_t i = 1; i < conv_centres_mev.size(); ++i) {
        report << " " << (conv_centres_mev[i] - conv_centres_mev[i - 1]);
    }
    report << "\n";

    // Three read-only diagnostic bin widths (section 7.11.2).
    for (size_t w = 0; w < args.bin_widths_mev.size(); ++w) {
        const double bw = args.bin_widths_mev[w];
        const int nbins = static_cast<int>(
            std::ceil((hi_mev - lo_mev) / bw));
        report << "diagnostic_bin_width_mev=" << bw << "\n";
        for (int b = 0; b < nbins; ++b) {
            const double bin_lo = lo_mev + b * bw;
            const double bin_hi = std::min(hi_mev, bin_lo + bw);
            BinAudit a;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double ke_mev =
                        partition.kinetic_energy[idx2(j, k)] /
                        (1.0e6 * Const::eV);
                    if (ke_mev < bin_lo || ke_mev >= bin_hi) continue;
                    if (partition.is_conversion(j, k)) {
                        ++a.conversion_center_count;
                        if (a.min_center_energy_mev < 0.0) {
                            a.min_center_energy_mev = ke_mev;
                        }
                        a.min_center_energy_mev =
                            std::min(a.min_center_energy_mev, ke_mev);
                        a.max_center_energy_mev =
                            std::max(a.max_center_energy_mev, ke_mev);
                        a.phase_volume_sum +=
                            electrons.cgrid.cell_phase_volume(j, k);
                    } else {
                        ++a.bulk_center_count;
                    }
                }
            }
            if (a.conversion_center_count > 0) {
                // Largest uncovered gap between adjacent conversion centres
                // inside this diagnostic bin.
                std::vector<double> centres;
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        if (!partition.is_conversion(j, k)) continue;
                        const double ke_mev =
                            partition.kinetic_energy[idx2(j, k)] /
                            (1.0e6 * Const::eV);
                        if (ke_mev < bin_lo || ke_mev >= bin_hi) continue;
                        centres.push_back(ke_mev);
                    }
                }
                std::sort(centres.begin(), centres.end());
                for (size_t i = 1; i < centres.size(); ++i) {
                    a.max_uncovered_energy_gap_mev =
                        std::max(a.max_uncovered_energy_gap_mev,
                                 centres[i] - centres[i - 1]);
                }
            }
            report << "bin_low_mev=" << bin_lo
                   << " bin_high_mev=" << bin_hi
                   << " bulk_center_count=" << a.bulk_center_count
                   << " conversion_center_count="
                   << a.conversion_center_count
                   << " min_center_energy_mev="
                   << a.min_center_energy_mev
                   << " max_center_energy_mev="
                   << a.max_center_energy_mev
                   << " max_uncovered_energy_gap_mev="
                   << a.max_uncovered_energy_gap_mev
                   << " phase_volume_sum=" << a.phase_volume_sum << "\n";
        }
    }

    // Per-ring threshold crossings are the required near-axis support
    // audit.  Occupancy is represented by the production phase volume; no
    // synthetic distribution or particle loading participates here.
    for (size_t r = 0; r < args.uperp_rings.size(); ++r) {
        const int k = args.uperp_rings[r];
        if (k < 0 || k >= Param::Nmu) { pass = false; continue; }
        double below = -1.0;
        double above = -1.0;
        double occupancy = 0.0;
        double pre_conversion_number = 0.0;
        double pre_conversion_energy = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            const double ke_mev = partition.kinetic_energy[idx2(j, k)] /
                                  (1.0e6 * Const::eV);
            if (ke_mev <= k_out_mev && (below < 0.0 || ke_mev > below)) below = ke_mev;
            if (ke_mev >= k_out_mev && (above < 0.0 || ke_mev < above)) above = ke_mev;
            if (partition.is_conversion(j, k)) {
                occupancy += electrons.cgrid.cell_phase_volume(j, k);
                const double mass = electrons.f[static_cast<size_t>(grid.nghost) *
                                                Param::Nvmu + idx2(j, k)];
                pre_conversion_number += mass;
                pre_conversion_energy += mass * partition.kinetic_energy[idx2(j, k)];
            }
        }
        const double gap = (below >= 0.0 && above >= 0.0) ? above - below : -1.0;
        report << "uperp_index=" << k
               << " uperp_center=" << electrons.cgrid.uperp_cells[k]
               << " nearest_below_Kout_mev=" << below
               << " nearest_above_Kout_mev=" << above
               << " crossing_gap_mev=" << gap
               << " occupancy_phase_volume=" << occupancy
               << " pre_conversion_number=" << pre_conversion_number
               << " pre_conversion_energy_J=" << pre_conversion_energy
               << "\n";
        if (!(below >= 0.0 && above >= 0.0 && std::isfinite(gap))) pass = false;
    }

    std::cout << report.str();
    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str());
        if (out) {
            out << report.str();
            out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
