// Section 7.11.16C / 16B: control cases for the fixed-(192,64) offline
// prototype.  Every case calls the production functions
// (build_tail_interface_grid_candidate / replay_tail_interface_histogram /
// tail_interface_remap_masses); the test never re-implements the monitor
// integral, the face inversion or the cylindrical overlap rule.  Hand
// computed references are limited to simple exact fractions (half a cell,
// unit partition) and to the shared production moment formulas.
//
// Usage:
//   tail_interface_grid_design_test --case all [--result path]

#include "grid.h"
#include "tail_interface_grid_design.h"
#include "tail_moment_constraint.h"
#include "tail_subcell_quadrature.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string test_case;
    std::string result;
    Args() : test_case("all") {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) args.test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) args.result = argv[++i];
        else return false;
    }
    return args.test_case == "all" ||
           args.test_case == "g0-identity" ||
           args.test_case == "single-cell-overlap" ||
           args.test_case == "constant-cells" ||
           args.test_case == "symmetric-pair" ||
           args.test_case == "maxwellian-roundtrip" ||
           args.test_case == "drift-maxwellian-roundtrip";
}

std::array<double, BULK_TAIL_MOMENT_COUNT> unit_moments(double upar,
                                                        double uperp)
{
    std::array<double, BULK_TAIL_MOMENT_COUNT> r = {};
    mass_cell_moments(1.0, upar, uperp, r[0], r[1], r[3], r[2], r[4], r[5]);
    return r;
}

double conversion_energy_j()
{
    return 6.0e6 * Const::eV;
}

// Build a small request histogram around the real production conversion
// cell (iv=185, u_parallel ~= 12.7, low u_perp).
std::vector<BulkTailVelocityBinAudit> conversion_histogram(
    const CylindricalVelocityGrid& cgrid)
{
    std::vector<BulkTailVelocityBinAudit> bins;
    int iv = 185;
    int imu = 2;
    if (iv >= static_cast<int>(cgrid.upar_cells.size()) ||
        imu >= static_cast<int>(cgrid.uperp_cells.size())) {
        iv = static_cast<int>(cgrid.upar_cells.size()) - 2;
        imu = 0;
    }
    const double masses[] = {1.0e7, 1.5e7, 2.0e7, 2.5e7};
    const int imus[] = {1, 2, 3, 4};
    for (int q = 0; q < 4; ++q) {
        BulkTailVelocityBinAudit bin;
        bin.iv = iv;
        bin.imu = (imus[q] < static_cast<int>(cgrid.uperp_cells.size()))
                      ? imus[q] : imu;
        bin.request_cell_count = 1;
        bin.request_number = masses[q];
        bins.push_back(bin);
    }
    return bins;
}

bool run_g0_identity(const CylindricalVelocityGrid& cgrid,
                     std::ostringstream& report)
{
    const std::vector<BulkTailVelocityBinAudit> histogram =
        conversion_histogram(cgrid);
    TailInterfaceGridDesignConfig config;
    const TailInterfaceGridCandidate candidate =
        build_tail_interface_grid_candidate(cgrid, histogram, config, "G0");
    const TailInterfaceReplayResult result = replay_tail_interface_histogram(
        cgrid, candidate, histogram, conversion_energy_j());
    // Direct reference from the input histogram using the shared production
    // moment formulas (no overlap formula is re-implemented here).
    std::array<long double, BULK_TAIL_MOMENT_COUNT> center_ref = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> volume_ref = {};
    for (size_t b = 0; b < histogram.size(); ++b) {
        const int iv = histogram[b].iv;
        const int imu = histogram[b].imu;
        const double m = histogram[b].request_number;
        const std::array<double, BULK_TAIL_MOMENT_COUNT> center =
            unit_moments(cgrid.upar_cells[iv], cgrid.uperp_cells[imu]);
        for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x)
            center_ref[static_cast<size_t>(x)] += center[x] * m;
        const std::vector<TailSubcellNode> nodes =
            TailSubcellQuadrature::nodes(cgrid, iv, imu);
        for (size_t q = 0; q < nodes.size(); ++q) {
            const std::array<double, BULK_TAIL_MOMENT_COUNT> col =
                unit_moments(nodes[q].upar, nodes[q].uperp);
            for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x)
                volume_ref[static_cast<size_t>(x)] +=
                    col[x] * m * nodes[q].mass_fraction;
        }
    }
    bool pass = candidate.valid &&
                candidate.upar_faces == cgrid.upar_faces &&
                candidate.uperp_faces == cgrid.uperp_faces &&
                result.g0_identity_ok &&
                result.number_residual <= 1.0e-14 &&
                result.negative_mass_cells == 0;
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
        const double scale =
            std::max(1.0, std::max(std::fabs(result.center_total[x]),
                                   std::fabs(static_cast<double>(center_ref[x]))));
        pass = pass &&
               std::fabs(result.center_total[x] -
                         static_cast<double>(center_ref[x])) <=
                   1.0e-12 * scale &&
               std::fabs(result.volume_total[x] -
                         static_cast<double>(volume_ref[x])) <=
                   1.0e-12 * scale;
    }
    report << "case=g0-identity valid=" << (candidate.valid ? 1 : 0)
           << " g0_identity_ok=" << (result.g0_identity_ok ? 1 : 0)
           << " number_residual=" << std::setprecision(17)
           << result.number_residual
           << " R_L1_Piperp=" << result.r_l1[BULK_TAIL_MOMENT_PIPERP]
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_single_cell_overlap(const CylindricalVelocityGrid& cgrid,
                             std::ostringstream& report)
{
    const int nv = static_cast<int>(cgrid.upar_faces.size()) - 1;
    const int nmu = static_cast<int>(cgrid.uperp_faces.size()) - 1;
    const int j0 = 185;
    const int k0 = 2;
    const double mass = 2.5;
    // Hand-crafted candidate: move the right face of cell j0 to the exact
    // midpoint of that cell.  The old cell then overlaps two new cells with
    // the exact fractions 1/2 and 1/2 (analytic, no overlap formula here).
    TailInterfaceGridCandidate candidate;
    candidate.grid_name = "analytic-split";
    candidate.upar_faces = cgrid.upar_faces;
    candidate.uperp_faces = cgrid.uperp_faces;
    const double lo = candidate.upar_faces[static_cast<size_t>(j0)];
    const double hi = candidate.upar_faces[static_cast<size_t>(j0) + 1];
    const double mid = 0.5 * (lo + hi);
    candidate.upar_faces[static_cast<size_t>(j0) + 1] = mid;
    candidate.valid = true;
    std::vector<double> from_mass(static_cast<size_t>(nv) * nmu, 0.0);
    from_mass[static_cast<size_t>(j0) * nmu + k0] = mass;
    double partition_error = 0.0;
    const std::vector<double> to_mass = tail_interface_remap_masses(
        cgrid.upar_faces, cgrid.uperp_faces, candidate.upar_faces,
        candidate.uperp_faces, from_mass, &partition_error);
    double left_mass = 0.0, right_mass = 0.0, total = 0.0;
    for (int k = 0; k < nmu; ++k) {
        left_mass += to_mass[static_cast<size_t>(j0) * nmu + k];
        right_mass += to_mass[static_cast<size_t>(j0 + 1) * nmu + k];
        total += to_mass[static_cast<size_t>(j0) * nmu + k];
        total += to_mass[static_cast<size_t>(j0 + 1) * nmu + k];
    }
    const bool pass = partition_error <= 1.0e-14 &&
                      std::fabs(left_mass - 0.5 * mass) <=
                          1.0e-14 * mass &&
                      std::fabs(right_mass - 0.5 * mass) <=
                          1.0e-14 * mass &&
                      std::fabs(total - mass) <= 1.0e-14 * mass;
    report << "case=single-cell-overlap left=" << std::setprecision(17)
           << left_mass << " right=" << right_mass << " total=" << total
           << " partition_error=" << partition_error
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_constant_cells(const CylindricalVelocityGrid& cgrid,
                        std::ostringstream& report)
{
    const int nv = static_cast<int>(cgrid.upar_cells.size());
    const int nmu = static_cast<int>(cgrid.uperp_cells.size());
    std::vector<BulkTailVelocityBinAudit> histogram;
    const double m0 = 3.0;
    for (int j = 64; j < 129 && j < nv; ++j) {
        for (int k = 0; k < 9 && k < nmu; ++k) {
            BulkTailVelocityBinAudit bin;
            bin.iv = j;
            bin.imu = k;
            bin.request_cell_count = 1;
            bin.request_number = m0;
            histogram.push_back(bin);
        }
    }
    TailInterfaceGridDesignConfig config;
    config.ax = 1.0;
    config.aperp = 1.0;
    config.sigma_x_cells = 2.0;
    config.sigma_perp_cells = 2.0;
    const TailInterfaceGridCandidate candidate =
        build_tail_interface_grid_candidate(cgrid, histogram, config, "G2");
    const std::vector<double> agg =
        tail_interface_aggregate_histogram(histogram, nv, nmu);
    double input_total = 0.0;
    for (size_t s = 0; s < agg.size(); ++s) input_total += agg[s];
    double partition_error = 0.0;
    const std::vector<double> to_mass = tail_interface_remap_masses(
        cgrid.upar_faces, cgrid.uperp_faces, candidate.upar_faces,
        candidate.uperp_faces, agg, &partition_error);
    double output_total = 0.0;
    double max_negative = 0.0;
    for (size_t s = 0; s < to_mass.size(); ++s) {
        output_total += to_mass[s];
        if (to_mass[s] < 0.0) max_negative = std::min(max_negative, to_mass[s]);
    }
    // Control case 3 (section 7.11.16C): only the remap conservation is
    // gated; the monitor candidate may legitimately violate the width
    // constraints for this synthetic histogram (reported for information).
    const bool pass = partition_error <= 1.0e-14 &&
                      std::fabs(output_total - input_total) <=
                          1.0e-14 * input_total &&
                      max_negative == 0.0;
    report << "case=constant-cells valid=" << (candidate.valid ? 1 : 0)
           << " input_total=" << std::setprecision(17) << input_total
           << " output_total=" << output_total
           << " partition_error=" << partition_error
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_symmetric_pair(const CylindricalVelocityGrid& cgrid,
                        std::ostringstream& report)
{
    const int nv = static_cast<int>(cgrid.upar_cells.size());
    const int nmu = static_cast<int>(cgrid.uperp_cells.size());
    const int j0 = 185;
    const int mirror = nv - 1 - j0;
    const int k0 = 2;
    const double mass = 1.75;
    std::vector<BulkTailVelocityBinAudit> histogram(2);
    histogram[0].iv = j0; histogram[0].imu = k0;
    histogram[0].request_cell_count = 1; histogram[0].request_number = mass;
    histogram[1].iv = mirror; histogram[1].imu = k0;
    histogram[1].request_cell_count = 1; histogram[1].request_number = mass;
    TailInterfaceGridDesignConfig config;
    config.ax = 1.0;
    config.sigma_x_cells = 2.0;
    const TailInterfaceGridCandidate candidate =
        build_tail_interface_grid_candidate(cgrid, histogram, config, "Gx");
    const std::vector<double> agg =
        tail_interface_aggregate_histogram(histogram, nv, nmu);
    double partition_error = 0.0;
    const std::vector<double> to_mass = tail_interface_remap_masses(
        cgrid.upar_faces, cgrid.uperp_faces, candidate.upar_faces,
        candidate.uperp_faces, agg, &partition_error);
    double px_signed = 0.0, px_l1 = 0.0;
    double jx_signed = 0.0, jx_l1 = 0.0;
    double mirror_error = 0.0;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = to_mass[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            const double upar = cgrid.upar_cells[static_cast<size_t>(j)];
            const double uperp = cgrid.uperp_cells[static_cast<size_t>(k)];
            const double gamma = std::sqrt(1.0 + upar * upar + uperp * uperp);
            px_signed += Const::me * Const::c * m * upar;
            px_l1 += std::fabs(Const::me * Const::c * m * upar);
            jx_signed += -Const::qe * Const::c * m * upar / gamma;
            jx_l1 += std::fabs(Const::qe * Const::c * m * upar / gamma);
            const int jm = nv - 1 - j;
            mirror_error = std::max(
                mirror_error,
                std::fabs(m - to_mass[static_cast<size_t>(jm) * nmu + k]));
        }
    }
    const bool pass = partition_error <= 1.0e-14 &&
                      std::fabs(px_signed) <= 2.0e-13 * px_l1 &&
                      std::fabs(jx_signed) <= 2.0e-13 * jx_l1 &&
                      px_l1 > 0.0 && jx_l1 > 0.0 &&
                      mirror_error <= 1.0e-14 * mass;
    report << "case=symmetric-pair valid=" << (candidate.valid ? 1 : 0)
           << " signed_Px=" << std::setprecision(17) << px_signed
           << " Px_l1=" << px_l1 << " signed_Jx=" << jx_signed
           << " Jx_l1=" << jx_l1 << " mirror_error=" << mirror_error
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

std::vector<BulkTailVelocityBinAudit> maxwellian_histogram(
    const CylindricalVelocityGrid& cgrid, double drift_u)
{
    const int nv = static_cast<int>(cgrid.upar_cells.size());
    const int nmu = static_cast<int>(cgrid.uperp_cells.size());
    const double u_th2 =
        Param::temperature_e / (Const::me * Const::c * Const::c);
    double normalization = 0.0;
    std::vector<double> mass(static_cast<size_t>(nv) * nmu, 0.0);
    for (int j = 0; j < nv; ++j) {
        const double du = cgrid.upar_cells[static_cast<size_t>(j)] - drift_u;
        for (int k = 0; k < nmu; ++k) {
            const double dp = cgrid.uperp_cells[static_cast<size_t>(k)];
            const double u2 = du * du + dp * dp;
            mass[static_cast<size_t>(j) * nmu + k] =
                cgrid.cell_phase_volume(j, k) *
                std::exp(-u2 / (2.0 * u_th2));
            normalization += mass[static_cast<size_t>(j) * nmu + k];
        }
    }
    std::vector<BulkTailVelocityBinAudit> histogram;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = mass[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            BulkTailVelocityBinAudit bin;
            bin.iv = j;
            bin.imu = k;
            bin.request_cell_count = 1;
            bin.request_number = m / normalization;
            histogram.push_back(bin);
        }
    }
    return histogram;
}

bool run_maxwellian_roundtrip(const CylindricalVelocityGrid& cgrid,
                              bool drift, std::ostringstream& report)
{
    const int nv = static_cast<int>(cgrid.upar_cells.size());
    const int nmu = static_cast<int>(cgrid.uperp_cells.size());
    const double drift_u = drift ? 3.0 : 0.0;
    const std::vector<BulkTailVelocityBinAudit> histogram =
        maxwellian_histogram(cgrid, drift_u);
    const std::vector<double> agg =
        tail_interface_aggregate_histogram(histogram, nv, nmu);
    TailInterfaceGridDesignConfig config;
    config.ax = 1.0;
    config.aperp = 1.0;
    config.sigma_x_cells = 2.0;
    config.sigma_perp_cells = 2.0;
    const TailInterfaceGridCandidate candidate =
        build_tail_interface_grid_candidate(cgrid, histogram, config, "G2");
    double partition_error = 0.0;
    const std::vector<double> forward = tail_interface_remap_masses(
        cgrid.upar_faces, cgrid.uperp_faces, candidate.upar_faces,
        candidate.uperp_faces, agg, &partition_error);
    const std::vector<double> back = tail_interface_remap_masses(
        candidate.upar_faces, candidate.uperp_faces, cgrid.upar_faces,
        cgrid.uperp_faces, forward, &partition_error);
    double input_total = 0.0, back_total = 0.0, max_negative = 0.0;
    for (size_t s = 0; s < agg.size(); ++s) input_total += agg[s];
    for (size_t s = 0; s < back.size(); ++s) {
        back_total += back[s];
        if (back[s] < 0.0) max_negative = std::min(max_negative, back[s]);
    }
    // Six-moment round-trip errors are reported, not enforced.
    std::array<long double, BULK_TAIL_MOMENT_COUNT> ref = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> got = {};
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const size_t slot = static_cast<size_t>(j) * nmu + k;
            const std::array<double, BULK_TAIL_MOMENT_COUNT> col =
                unit_moments(cgrid.upar_cells[static_cast<size_t>(j)],
                             cgrid.uperp_cells[static_cast<size_t>(k)]);
            for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
                ref[static_cast<size_t>(x)] += col[x] * agg[slot];
                got[static_cast<size_t>(x)] += col[x] * back[slot];
            }
        }
    }
    std::array<double, BULK_TAIL_MOMENT_COUNT> rel = {};
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
        rel[static_cast<size_t>(x)] =
            std::fabs(static_cast<double>(got[static_cast<size_t>(x)] -
                                          ref[static_cast<size_t>(x)])) /
            std::max(1.0, std::fabs(static_cast<double>(
                              ref[static_cast<size_t>(x)])));
    }
    const bool pass = std::fabs(back_total - input_total) <=
                          1.0e-12 * input_total &&
                      max_negative >= -1.0e-15 * input_total;
    const char* comp[] = {"N", "Px", "Jx", "K", "Pixx", "Piperp"};
    report << "case=" << (drift ? "drift-maxwellian-roundtrip"
                                : "maxwellian-roundtrip")
           << " valid=" << (candidate.valid ? 1 : 0)
           << " input_total=" << std::setprecision(17) << input_total
           << " back_total=" << back_total
           << " number_rel_error="
           << std::fabs(back_total - input_total) /
                  std::max(input_total, 1.0e-300);
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x)
        report << " " << comp[x] << "_rel_error=" << std::setprecision(17)
               << rel[static_cast<size_t>(x)];
    report << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

} // namespace

int main(int argc, char** argv)
{
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: tail_interface_grid_design_test --case all "
                     "|g0-identity|single-cell-overlap|constant-cells|"
                     "symmetric-pair|maxwellian-roundtrip|"
                     "drift-maxwellian-roundtrip [--result path]\n";
        return 2;
    }
    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    const bool all = args.test_case == "all";
    std::ostringstream report;
    bool pass = true;
    if (all || args.test_case == "g0-identity")
        pass = run_g0_identity(cgrid, report) && pass;
    if (all || args.test_case == "single-cell-overlap")
        pass = run_single_cell_overlap(cgrid, report) && pass;
    if (all || args.test_case == "constant-cells")
        pass = run_constant_cells(cgrid, report) && pass;
    if (all || args.test_case == "symmetric-pair")
        pass = run_symmetric_pair(cgrid, report) && pass;
    if (all || args.test_case == "maxwellian-roundtrip")
        pass = run_maxwellian_roundtrip(cgrid, false, report) && pass;
    if (all || args.test_case == "drift-maxwellian-roundtrip")
        pass = run_maxwellian_roundtrip(cgrid, true, report) && pass;
    report << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    std::cout << report.str();
    if (!args.result.empty()) {
        std::ofstream out(args.result.c_str(), std::ios::trunc);
        if (!out) {
            std::cerr << "cannot write result file: " << args.result << "\n";
            return 1;
        }
        out << report.str();
    }
    return pass ? 0 : 1;
}
