#include "grid.h"
#include "tail_subcell_quadrature.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<double> energy_edges(double width_mev, double max_energy)
{
    const double width = width_mev * 1.0e6 * Const::eV;
    const size_t bins = static_cast<size_t>(std::ceil(max_energy / width)) + 1;
    std::vector<double> edges(bins + 1, 0.0);
    for (size_t i = 0; i < edges.size(); ++i)
        edges[i] = static_cast<double>(i) * width;
    return edges;
}

std::vector<double> dense(
    const std::vector<TailEnergyBinFraction>& fractions, size_t bins)
{
    std::vector<double> result(bins, 0.0);
    for (size_t i = 0; i < fractions.size(); ++i) {
        if (fractions[i].bin >= 0 &&
            static_cast<size_t>(fractions[i].bin) < bins)
            result[static_cast<size_t>(fractions[i].bin)] +=
                fractions[i].mass_fraction;
    }
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    std::string result_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result" && i + 1 < argc) result_path = argv[++i];
        else {
            std::cerr << "usage: tail_energy_shell_overlap_test "
                         "[--result path]\n";
            return 2;
        }
    }

    CylindricalVelocityGrid grid;
    grid.init(Param::momentum_umax);
    const double max_upar = std::max(std::fabs(grid.upar_faces.front()),
                                     std::fabs(grid.upar_faces.back()));
    const double max_uperp = grid.uperp_faces.back();
    const double max_energy =
        (std::sqrt(1.0 + max_upar * max_upar + max_uperp * max_uperp) - 1.0) *
        Const::me * Const::c * Const::c * 1.001;
    const std::vector<double> fine_edges = energy_edges(0.05, max_energy);
    const std::vector<double> coarse_edges = energy_edges(0.1, max_energy);

    double max_sum_error = 0.0;
    double min_fraction = 1.0;
    double max_symmetry_error = 0.0;
    double max_rebin_error = 0.0;
    int empty_cell_count = 0;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const std::vector<TailEnergyBinFraction> fine =
                TailSubcellQuadrature::energy_bin_fractions(
                    grid, iv, imu, fine_edges);
            const std::vector<TailEnergyBinFraction> coarse =
                TailSubcellQuadrature::energy_bin_fractions(
                    grid, iv, imu, coarse_edges);
            if (fine.empty() || coarse.empty()) ++empty_cell_count;
            double sum = 0.0;
            for (size_t q = 0; q < fine.size(); ++q) {
                sum += fine[q].mass_fraction;
                min_fraction = std::min(min_fraction, fine[q].mass_fraction);
            }
            max_sum_error = std::max(max_sum_error, std::fabs(sum - 1.0));

            const std::vector<double> fine_dense =
                dense(fine, fine_edges.size() - 1);
            const std::vector<double> coarse_dense =
                dense(coarse, coarse_edges.size() - 1);
            for (size_t b = 0; b < coarse_dense.size(); ++b) {
                const double aggregated =
                    (2 * b < fine_dense.size() ? fine_dense[2 * b] : 0.0) +
                    (2 * b + 1 < fine_dense.size() ?
                        fine_dense[2 * b + 1] : 0.0);
                max_rebin_error = std::max(
                    max_rebin_error, std::fabs(aggregated - coarse_dense[b]));
            }

            const int mirror_iv = Param::Nv - 1 - iv;
            const std::vector<double> mirror_dense = dense(
                TailSubcellQuadrature::energy_bin_fractions(
                    grid, mirror_iv, imu, fine_edges),
                fine_edges.size() - 1);
            for (size_t b = 0; b < fine_dense.size(); ++b)
                max_symmetry_error = std::max(max_symmetry_error,
                    std::fabs(fine_dense[b] - mirror_dense[b]));
        }
    }

    const bool pass = empty_cell_count == 0 && min_fraction >= -1.0e-13 &&
        max_sum_error <= 2.0e-12 && max_symmetry_error <= 1.0e-10 &&
        max_rebin_error <= 2.0e-12;
    std::ostringstream report;
    report << std::setprecision(17)
        << "empty_cell_count=" << empty_cell_count << "\n"
        << "min_fraction=" << min_fraction << "\n"
        << "max_sum_error=" << max_sum_error << "\n"
        << "max_symmetry_error=" << max_symmetry_error << "\n"
        << "max_rebin_error=" << max_rebin_error << "\n"
        << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    std::cout << report.str();
    if (!result_path.empty()) {
        std::ofstream out(result_path.c_str());
        if (!out) return 2;
        out << report.str();
    }
    return pass ? 0 : 1;
}
