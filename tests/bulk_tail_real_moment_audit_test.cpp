#include "bulk_tail_moment_audit.h"
#include "bulk_tail_moment_audit_io.h"
#include "species.h"
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
#include <mpi.h>

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
    return args.test_case == "all" || args.test_case == "single-cell-known" ||
           args.test_case == "symmetric-pair" ||
           args.test_case == "threshold-face" ||
           args.test_case == "clipped-cell" ||
           args.test_case == "accepted-only" ||
           args.test_case == "zero-denominator";
}

std::array<double, BULK_TAIL_MOMENT_COUNT> moments(double mass,
                                                    double upar,
                                                    double uperp)
{
    std::array<double, BULK_TAIL_MOMENT_COUNT> r = {};
    mass_cell_moments(mass, upar, uperp, r[0], r[1], r[3], r[2], r[4], r[5]);
    return r;
}

bool close(double lhs, double rhs, double tolerance = 2.0e-13)
{
    return std::fabs(lhs - rhs) <=
           tolerance * std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
}

bool close_array(const std::array<double, BULK_TAIL_MOMENT_COUNT>& lhs,
                 const std::array<double, BULK_TAIL_MOMENT_COUNT>& rhs,
                 double tolerance = 2.0e-13)
{
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
        if (!close(lhs[m], rhs[m], tolerance)) return false;
    return true;
}

void add(std::array<double, BULK_TAIL_MOMENT_COUNT>& dst,
         const std::array<double, BULK_TAIL_MOMENT_COUNT>& src)
{
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) dst[m] += src[m];
}

bool find_conversion_pair(const Species& bulk,
                          const HybridVelocityPartition& partition,
                          int& iv, int& mirror, int& imu)
{
    for (int k = 0; k < Param::Nmu; ++k) {
        for (int j = 0; j < Param::Nv; ++j) {
            if (!(bulk.cgrid.upar_cells[j] > 0.0) ||
                !partition.is_conversion(j, k)) continue;
            int best = -1;
            double error = 1.0e300;
            for (int q = 0; q < Param::Nv; ++q) {
                const double candidate =
                    std::fabs(bulk.cgrid.upar_cells[q] + bulk.cgrid.upar_cells[j]);
                if (candidate < error) { error = candidate; best = q; }
            }
            if (best >= 0 && partition.is_conversion(best, k) &&
                error <= 1.0e-12 * std::max(1.0, std::fabs(bulk.cgrid.upar_cells[j]))) {
                iv = j; mirror = best; imu = k; return true;
            }
        }
    }
    return false;
}

bool run_single_cell(const Species& bulk,
                     const HybridVelocityPartition& partition,
                     int iv, int imu, std::ostringstream& report)
{
    const double mass = 2.75;
    BulkTailMomentAuditRequest request = {3, iv, imu, mass};
    const MomentRepresentationAudit audit = bulk_tail_audit_conversion_requests(
        bulk.cgrid, partition, std::vector<BulkTailMomentAuditRequest>(1, request),
        0, 8);
    const std::array<double, BULK_TAIL_MOMENT_COUNT> center =
        moments(mass, bulk.cgrid.upar_cells[iv], bulk.cgrid.uperp_cells[imu]);
    std::array<double, BULK_TAIL_MOMENT_COUNT> volume = {};
    std::array<double, BULK_TAIL_MOMENT_COUNT> eligible = {};
    double eligible_mass = 0.0;
    const std::vector<TailSubcellNode> nodes =
        TailSubcellQuadrature::nodes(bulk.cgrid, iv, imu);
    for (size_t q = 0; q < nodes.size(); ++q) {
        const double weight = mass * nodes[q].mass_fraction;
        add(volume, moments(weight, nodes[q].upar, nodes[q].uperp));
        if (nodes[q].kinetic_energy >= partition.min_conversion_energy) {
            add(eligible, moments(weight, nodes[q].upar, nodes[q].uperp));
            eligible_mass += weight;
        }
    }
    std::array<double, BULK_TAIL_MOMENT_COUNT> normalized = eligible;
    if (eligible_mass > 0.0)
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
            normalized[m] *= mass / eligible_mass;
    const bool pass = audit.finite && audit.request_cell_count == 1 &&
        audit.positive_request_cell_count == 1 && close_array(audit.center, center) &&
        close_array(audit.volume, volume) && close_array(audit.eligible_raw, eligible) &&
        close_array(audit.eligible_normalized, normalized) &&
        !audit.top_cells.empty() && close_array(audit.top_cells[0].center, center) &&
        close_array(audit.top_cells[0].volume, volume);
    report << "case=single-cell-known nonunit_mass=" << mass
           << " request_cells=" << audit.request_cell_count
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_symmetric_pair(const Species& bulk,
                        const HybridVelocityPartition& partition,
                        int iv, int mirror, int imu,
                        std::ostringstream& report)
{
    const double mass = 1.75;
    std::vector<BulkTailMomentAuditRequest> requests(2);
    requests[0] = BulkTailMomentAuditRequest{2, iv, imu, mass};
    requests[1] = BulkTailMomentAuditRequest{7, mirror, imu, mass};
    const MomentRepresentationAudit audit = bulk_tail_audit_conversion_requests(
        bulk.cgrid, partition, requests, 0, 8);
    const double px_scale = std::max(1.0, audit.center_l1[BULK_TAIL_MOMENT_PX]);
    const double jx_scale = std::max(1.0, audit.center_l1[BULK_TAIL_MOMENT_JX]);
    const bool pass = audit.finite && audit.positive_request_cell_count == 2 &&
        std::fabs(audit.center[BULK_TAIL_MOMENT_PX]) <= 2.0e-13 * px_scale &&
        std::fabs(audit.center[BULK_TAIL_MOMENT_JX]) <= 2.0e-13 * jx_scale &&
        audit.center_l1[BULK_TAIL_MOMENT_PX] > 0.0 &&
        audit.center_l1[BULK_TAIL_MOMENT_JX] > 0.0;
    report << "case=symmetric-pair signed_Px="
           << audit.center[BULK_TAIL_MOMENT_PX]
           << " Px_l1=" << audit.center_l1[BULK_TAIL_MOMENT_PX]
           << " signed_Jx=" << audit.center[BULK_TAIL_MOMENT_JX]
           << " Jx_l1=" << audit.center_l1[BULK_TAIL_MOMENT_JX]
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_threshold_case(const Species& bulk,
                        const HybridVelocityPartition& base_partition,
                        int iv, int imu, bool clipped,
                        std::ostringstream& report)
{
    HybridVelocityPartition partition = base_partition;
    const std::vector<TailSubcellNode> nodes =
        TailSubcellQuadrature::nodes(bulk.cgrid, iv, imu);
    if (nodes.empty()) return false;
    double min_energy = nodes[0].kinetic_energy;
    double max_energy = nodes[0].kinetic_energy;
    for (size_t q = 1; q < nodes.size(); ++q) {
        min_energy = std::min(min_energy, nodes[q].kinetic_energy);
        max_energy = std::max(max_energy, nodes[q].kinetic_energy);
    }
    partition.min_conversion_energy = clipped
        ? 0.5 * (min_energy + max_energy) : min_energy;
    const double mass = clipped ? 3.25 : 0.375;
    BulkTailMomentAuditRequest request = {5, iv, imu, mass};
    const MomentRepresentationAudit audit = bulk_tail_audit_conversion_requests(
        bulk.cgrid, partition, std::vector<BulkTailMomentAuditRequest>(1, request),
        0, 8);
    double eligible = 0.0;
    for (size_t q = 0; q < nodes.size(); ++q)
        if (nodes[q].kinetic_energy >= partition.min_conversion_energy)
            eligible += mass * nodes[q].mass_fraction;
    bool pass = audit.finite && close(audit.eligible_raw[BULK_TAIL_MOMENT_N], eligible) &&
                close(audit.eligible_number_fraction, eligible / mass);
    if (clipped) {
        pass = pass && eligible > 0.0 && eligible < mass &&
               close(audit.eligible_normalized[BULK_TAIL_MOMENT_N], mass) &&
               audit.below_threshold_number_fraction > 0.0;
    } else {
        pass = pass && close(eligible, mass) &&
               close(audit.below_threshold_number_fraction, 0.0);
    }
    report << "case=" << (clipped ? "clipped-cell" : "threshold-face")
           << " input_N=" << mass << " eligible_raw_N=" << eligible
           << " eligible_normalized_N="
           << audit.eligible_normalized[BULK_TAIL_MOMENT_N]
           << " below_fraction=" << audit.below_threshold_number_fraction
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_accepted_only(std::ostringstream& report)
{
    // A zero-request accepted event must reduce to has_positive_requests=false
    // (no row written); a non-empty accepted event must reduce to true.
    MomentRepresentationAudit empty;
    empty.enabled = true;
    BulkTailMomentAuditGlobal empty_global;
    const BulkTailMomentAuditIoStatus empty_status =
        bulk_tail_moment_audit_reduce_event(empty, 0, 8, 0, 1, empty_global);
    MomentRepresentationAudit nonempty;
    nonempty.enabled = true;
    nonempty.positive_request_cell_count = 1;
    nonempty.request_cell_count = 1;
    nonempty.center[BULK_TAIL_MOMENT_N] = 1.0;
    nonempty.center_l1[BULK_TAIL_MOMENT_N] = 1.0;
    BulkTailMomentAuditGlobal nonempty_global;
    const BulkTailMomentAuditIoStatus nonempty_status =
        bulk_tail_moment_audit_reduce_event(nonempty, 0, 8, 0, 1,
                                            nonempty_global);
    const bool pass =
        empty_status == BulkTailMomentAuditIoStatus::OK &&
        !empty_global.has_positive_requests &&
        nonempty_status == BulkTailMomentAuditIoStatus::OK &&
        nonempty_global.has_positive_requests &&
        nonempty_global.positive_request_cell_count == 1;
    report << "case=accepted-only empty_has_requests="
           << (empty_global.has_positive_requests ? 1 : 0)
           << " nonempty_has_requests="
           << (nonempty_global.has_positive_requests ? 1 : 0)
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool run_zero_denominator(std::ostringstream& report)
{
    // Section 7.11.16B: a component whose center L1 is exactly zero with a
    // nonzero delta L1 must produce relative_defined=0 (the ratio is +inf),
    // never a dimensional-constant denominator.
    MomentRepresentationAudit audit;
    audit.enabled = true;
    audit.positive_request_cell_count = 1;
    audit.request_cell_count = 1;
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        audit.center_l1[m] = 1.0;
        audit.delta_l1[m] = 0.0;
    }
    audit.center_l1[BULK_TAIL_MOMENT_PX] = 0.0;
    audit.delta_l1[BULK_TAIL_MOMENT_PX] = 5.0;
    BulkTailMomentAuditGlobal g;
    const BulkTailMomentAuditIoStatus status =
        bulk_tail_moment_audit_reduce_event(audit, 0, 8, 0, 1, g);
    bool pass = status == BulkTailMomentAuditIoStatus::OK &&
                g.has_positive_requests &&
                g.relative_defined[BULK_TAIL_MOMENT_PX] == 0;
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        if (m == BULK_TAIL_MOMENT_PX) continue;
        pass = pass && g.relative_defined[m] == 1;
    }
    report << "case=zero-denominator px_defined="
           << static_cast<int>(g.relative_defined[BULK_TAIL_MOMENT_PX])
           << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    if (size != 1 || !parse_args(argc, argv, args)) {
        if (rank == 0)
            std::cerr << "usage: bulk_tail_real_moment_audit_test "
                         "--case all|single-cell-known|symmetric-pair|"
                         "threshold-face|clipped-cell|accepted-only "
                         "[--result path]\n";
        MPI_Finalize(); return 2;
    }
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 20, 20.0 * Param::dx);
    Species bulk;
    bulk.init("moment_audit", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
              Const::me, Param::dens, Param::temperature_e, false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    int iv = -1, mirror = -1, imu = -1;
    const bool setup_ok = find_conversion_pair(bulk, partition, iv, mirror, imu);
    bool pass = setup_ok;
    std::ostringstream report;
    report << std::setprecision(17);
    if (setup_ok && (args.test_case == "all" || args.test_case == "single-cell-known"))
        pass = run_single_cell(bulk, partition, iv, imu, report) && pass;
    if (setup_ok && (args.test_case == "all" || args.test_case == "symmetric-pair"))
        pass = run_symmetric_pair(bulk, partition, iv, mirror, imu, report) && pass;
    if (setup_ok && (args.test_case == "all" || args.test_case == "threshold-face"))
        pass = run_threshold_case(bulk, partition, iv, imu, false, report) && pass;
    if (setup_ok && (args.test_case == "all" || args.test_case == "clipped-cell"))
        pass = run_threshold_case(bulk, partition, iv, imu, true, report) && pass;
    if (args.test_case == "all" || args.test_case == "accepted-only")
        pass = run_accepted_only(report) && pass;
    if (args.test_case == "all" || args.test_case == "zero-denominator")
        pass = run_zero_denominator(report) && pass;
    report << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    std::cout << report.str();
    if (!args.result.empty()) {
        std::ofstream out(args.result.c_str(), std::ios::trunc);
        if (!out) {
            std::cerr << "cannot write result file: " << args.result << "\n";
            MPI_Finalize(); return 3;
        }
        out << report.str();
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
