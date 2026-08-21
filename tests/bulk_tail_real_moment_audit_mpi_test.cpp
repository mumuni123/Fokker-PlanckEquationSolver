// Section 7.11.16B / 16A: MPI rank-consistency test for the real-moment
// audit.  A fixed set of physical conversion requests is partitioned across
// ranks by spatial ownership (each rank never duplicates a request), the
// production reduction (bulk_tail_moment_audit_io) is used directly, and
// the global moments, counts, max-cell relatives, velocity bins and the
// global top-N set are compared against the single-rank reference computed
// with the same production functions.  Deliberately inconsistent event
// layouts (event count and conversion location) must return
// INVALID_EVENT_LAYOUT on every rank without deadlock.
//
// Usage:
//   bulk_tail_real_moment_audit_mpi_test
//     --case all|decomposition [--result path]

#include "bulk_tail_moment_audit.h"
#include "bulk_tail_moment_audit_io.h"
#include "species.h"

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
    return args.test_case == "all" || args.test_case == "decomposition";
}

bool close(double lhs, double rhs, double tolerance = 5.0e-13)
{
    return std::fabs(lhs - rhs) <= tolerance *
        std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
}

bool owns(const SpatialGrid& grid, int ix_global)
{
    return ix_global >= grid.ix_start &&
           ix_global < grid.ix_start + grid.nx_local;
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
            if (best >= 0 && partition.is_conversion(best, k)) {
                iv = j; mirror = best; imu = k; return true;
            }
        }
    }
    return false;
}

bool same_bin(const std::vector<BulkTailVelocityBinAudit>& a,
              const std::vector<BulkTailVelocityBinAudit>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].iv != b[i].iv || a[i].imu != b[i].imu ||
            a[i].request_cell_count != b[i].request_cell_count ||
            !close(a[i].request_number, b[i].request_number)) {
            return false;
        }
    }
    return true;
}

bool same_top(const std::vector<BulkTailMomentAuditTopCell>& a,
              const std::vector<BulkTailMomentAuditTopCell>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].ix_global != b[i].ix_global || a[i].iv != b[i].iv ||
            a[i].imu != b[i].imu || !close(a[i].score, b[i].score)) {
            return false;
        }
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            if (!close(a[i].center[m], b[i].center[m]) ||
                !close(a[i].volume[m], b[i].volume[m])) {
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    if (!parse_args(argc, argv, args)) {
        if (rank == 0)
            std::cerr << "usage: bulk_tail_real_moment_audit_mpi_test "
                         "--case all|decomposition [--result path]\n";
        MPI_Finalize(); return 2;
    }

    const int nx_global = 20;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, nx_global, nx_global * Param::dx);
    Species bulk;
    bulk.init("moment_audit_mpi", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    int iv = -1, mirror = -1, imu = -1;
    bool setup_ok = find_conversion_pair(bulk, partition, iv, mirror, imu);

    // Fixed physical requests partitioned across ranks by spatial ownership.
    const int global_ix[] = {1, 6, 12, 18};
    const double masses[] = {0.25, 2.0, 3.5, 0.75};
    std::vector<BulkTailMomentAuditRequest> all_requests;
    std::vector<BulkTailMomentAuditRequest> local_requests;
    if (setup_ok) {
        for (int q = 0; q < 4; ++q) {
            BulkTailMomentAuditRequest request;
            request.ix_global = global_ix[q];
            request.iv = (q % 2 == 0) ? iv : mirror;
            request.imu = imu;
            request.mass = masses[q];
            all_requests.push_back(request);
            if (owns(grid, request.ix_global)) local_requests.push_back(request);
        }
    }
    const MomentRepresentationAudit local = bulk_tail_audit_conversion_requests(
        bulk.cgrid, partition, local_requests, rank, 64);
    const MomentRepresentationAudit reference = bulk_tail_audit_conversion_requests(
        bulk.cgrid, partition, all_requests, 0, 64);

    BulkTailMomentAuditGlobal global;
    const BulkTailMomentAuditIoStatus reduce_status =
        bulk_tail_moment_audit_reduce_event(local, 0, 64, rank, size, global);
    bool pass = setup_ok && local.finite &&
                reduce_status == BulkTailMomentAuditIoStatus::OK &&
                global.has_positive_requests;
    if (pass) {
        for (int a = 0; a < 7; ++a) {
            const double* g = NULL;
            const double* r = NULL;
            if (a == 0) { g = global.center.data(); r = reference.center.data(); }
            else if (a == 1) { g = global.volume.data(); r = reference.volume.data(); }
            else if (a == 2) { g = global.eligible_raw.data(); r = reference.eligible_raw.data(); }
            else if (a == 3) { g = global.eligible_normalized.data(); r = reference.eligible_normalized.data(); }
            else if (a == 4) { g = global.delta_signed.data(); r = reference.delta_signed.data(); }
            else if (a == 5) { g = global.delta_l1.data(); r = reference.delta_l1.data(); }
            else { g = global.center_l1.data(); r = reference.center_l1.data(); }
            for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
                pass = close(g[m], r[m]) && pass;
        }
        pass = (global.request_cell_count == reference.request_cell_count) &&
               (global.positive_request_cell_count ==
                reference.positive_request_cell_count) &&
               (global.volume_target_feasible_count ==
                reference.volume_target_feasible_count) &&
               (global.volume_target_failed_count ==
                reference.volume_target_failed_count) &&
               (global.eligible_target_feasible_count ==
                reference.eligible_target_feasible_count) &&
               (global.eligible_target_failed_count ==
                reference.eligible_target_failed_count) && pass;
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
            pass = close(global.max_cell_relative[m],
                         reference.max_cell_relative[m]) && pass;
        // The velocity bins and the global top-N set are assembled on rank 0
        // only; compare them there (the reduced moments/counts above are
        // available on every rank).
        if (rank == 0) {
            pass = same_bin(global.velocity_bins, reference.velocity_bins) &&
                   pass;
            pass = same_top(global.top_cells, reference.top_cells) && pass;
        }
    }

    // Deliberately inconsistent event layouts: all ranks must take the same
    // INVALID_EVENT_LAYOUT path without deadlock.  For a single rank the
    // checks are vacuous (OK).
    const int fabricated_count = 1 + (rank == 0 ? 1 : 0);
    const BulkTailMomentAuditIoStatus count_status =
        bulk_tail_moment_audit_check_event_layout(
            static_cast<size_t>(fabricated_count), rank, size);
    const bool count_invalid_ok =
        (size > 1)
            ? count_status == BulkTailMomentAuditIoStatus::INVALID_EVENT_LAYOUT
            : count_status == BulkTailMomentAuditIoStatus::OK;
    MomentRepresentationAudit location_audit;
    location_audit.enabled = true;
    location_audit.positive_request_cell_count = 1;
    location_audit.request_cell_count = 1;
    BulkTailMomentAuditGlobal location_global;
    const BulkTailMomentAuditIoStatus location_status =
        bulk_tail_moment_audit_reduce_event(
            location_audit, (rank == 0) ? 1 : 0, 64, rank, size,
            location_global);
    const bool location_invalid_ok =
        (size > 1)
            ? location_status == BulkTailMomentAuditIoStatus::INVALID_EVENT_LAYOUT
            : location_status == BulkTailMomentAuditIoStatus::OK;
    pass = pass && count_invalid_ok && location_invalid_ok;

    int local_pass = pass ? 1 : 0, all_pass = 0;
    MPI_Allreduce(&local_pass, &all_pass, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    std::ostringstream report;
    if (rank == 0) {
        report << std::setprecision(17)
               << "case=decomposition ranks=" << size
               << " nx_global=" << nx_global
               << " request_cells=" << global.request_cell_count
               << " global_N=" << global.center[BULK_TAIL_MOMENT_N]
               << " reference_N=" << reference.center[BULK_TAIL_MOMENT_N]
               << " velocity_bins=" << global.velocity_bins.size()
               << " top_cells=" << global.top_cells.size()
               << " count_invalid_ok=" << (count_invalid_ok ? 1 : 0)
               << " location_invalid_ok=" << (location_invalid_ok ? 1 : 0)
               << " status=" << (all_pass ? "PASS" : "FAIL") << "\n"
               << "status=" << (all_pass ? "PASS" : "FAIL") << "\n";
        std::cout << report.str();
        if (!args.result.empty()) {
            std::ofstream out(args.result.c_str(), std::ios::trunc);
            if (out) out << report.str();
            else all_pass = 0;
        }
    }
    MPI_Bcast(&all_pass, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!all_pass && rank == 0 && !args.result.empty())
        std::cerr << "MPI audit failed or cannot write result file: "
                  << args.result << "\n";
    MPI_Finalize();
    return all_pass ? 0 : 1;
}
