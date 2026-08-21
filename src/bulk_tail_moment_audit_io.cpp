#include "bulk_tail_moment_audit_io.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mpi.h>

namespace {

// Section 7.11.16B item 5 r-rule used for every reported relative quantity
// (recomputed from the global raw L1s; the per-cell flags live in the audit
// top cells).
double relative_ratio(double a, double b, unsigned char& defined)
{
    if (b > 0.0) {
        defined = 1;
        return a / b;
    }
    if (a == 0.0) {
        defined = 1;
        return 0.0;
    }
    defined = 0;
    return std::numeric_limits<double>::infinity();
}

} // namespace

BulkTailMomentAuditGlobal::BulkTailMomentAuditGlobal()
    : has_positive_requests(false), request_cell_count(0),
      positive_request_cell_count(0), volume_target_feasible_count(0),
      volume_target_failed_count(0), eligible_target_feasible_count(0),
      eligible_target_failed_count(0), center(), volume(), eligible_raw(),
      eligible_normalized(), delta_signed(), delta_l1(), center_l1(),
      max_cell_relative(), relative_defined(), eligible_number_fraction(0.0),
      threshold_window_number(0.0), velocity_bins(), top_cells()
{}

BulkTailMomentAuditIoStatus bulk_tail_moment_audit_check_event_layout(
    size_t local_event_count, int mpi_rank, int mpi_size)
{
    (void)mpi_rank;
    if (mpi_size <= 1) return BulkTailMomentAuditIoStatus::OK;
    const int local = static_cast<int>(local_event_count);
    int global_min = 0;
    int global_max = 0;
    MPI_Allreduce(&local, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return (global_min == global_max)
               ? BulkTailMomentAuditIoStatus::OK
               : BulkTailMomentAuditIoStatus::INVALID_EVENT_LAYOUT;
}

BulkTailMomentAuditIoStatus bulk_tail_moment_audit_reduce_event(
    const MomentRepresentationAudit& audit, int conversion_location,
    size_t top_cell_limit, int mpi_rank, int mpi_size,
    BulkTailMomentAuditGlobal& global)
{
    global = BulkTailMomentAuditGlobal();
    // 1) conversion-location MIN/MAX consistency.
    if (mpi_size > 1) {
        int local_location = conversion_location;
        int loc_min = 0;
        int loc_max = 0;
        MPI_Allreduce(&local_location, &loc_min, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&local_location, &loc_max, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
        if (loc_min != loc_max) {
            return BulkTailMomentAuditIoStatus::INVALID_EVENT_LAYOUT;
        }
    }
    // 2) global positive-request count; zero -> no further collective.
    unsigned long long local_positive =
        static_cast<unsigned long long>(audit.positive_request_cell_count);
    unsigned long long global_positive = 0;
    MPI_Allreduce(&local_positive, &global_positive, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    global.positive_request_cell_count = global_positive;
    if (global_positive == 0) return BulkTailMomentAuditIoStatus::OK;
    global.has_positive_requests = true;

    // 3) fixed-length SUM: raw moments and L1s.
    const int kSumLen = 7 * BULK_TAIL_MOMENT_COUNT + 2;
    double local_sum[7 * BULK_TAIL_MOMENT_COUNT + 2] = {};
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        local_sum[m] = audit.center[m];
        local_sum[6 + m] = audit.volume[m];
        local_sum[12 + m] = audit.eligible_raw[m];
        local_sum[18 + m] = audit.eligible_normalized[m];
        local_sum[24 + m] = audit.delta_signed[m];
        local_sum[30 + m] = audit.delta_l1[m];
        local_sum[36 + m] = audit.center_l1[m];
    }
    local_sum[42] = audit.eligible_number_fraction *
                    audit.center[BULK_TAIL_MOMENT_N];
    local_sum[43] = audit.threshold_window_number;
    double global_sum[7 * BULK_TAIL_MOMENT_COUNT + 2] = {};
    MPI_Allreduce(local_sum, global_sum, kSumLen, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        global.center[m] = global_sum[m];
        global.volume[m] = global_sum[6 + m];
        global.eligible_raw[m] = global_sum[12 + m];
        global.eligible_normalized[m] = global_sum[18 + m];
        global.delta_signed[m] = global_sum[24 + m];
        global.delta_l1[m] = global_sum[30 + m];
        global.center_l1[m] = global_sum[36 + m];
    }
    global.threshold_window_number = global_sum[43];
    const double center_n = global.center[BULK_TAIL_MOMENT_N];
    global.eligible_number_fraction =
        center_n > 0.0
            ? global.eligible_raw[BULK_TAIL_MOMENT_N] / center_n
            : 0.0;

    // 4) per-component cell-max.
    double local_max[BULK_TAIL_MOMENT_COUNT] = {};
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        local_max[m] = audit.max_cell_relative[m];
    }
    MPI_Allreduce(local_max, global.max_cell_relative.data(),
                  BULK_TAIL_MOMENT_COUNT, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    // Event-level relative-defined flags from the global raw L1s.
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        unsigned char defined = 1;
        relative_ratio(global.delta_l1[m], global.center_l1[m], defined);
        global.relative_defined[m] = defined;
    }

    // 5) global counts (feasibility etc.).
    unsigned long long local_counts[6] = {
        static_cast<unsigned long long>(audit.request_cell_count),
        static_cast<unsigned long long>(audit.positive_request_cell_count),
        static_cast<unsigned long long>(audit.volume_target_feasible_count),
        static_cast<unsigned long long>(audit.volume_target_failed_count),
        static_cast<unsigned long long>(audit.eligible_target_feasible_count),
        static_cast<unsigned long long>(audit.eligible_target_failed_count)
    };
    unsigned long long global_counts[6] = {};
    MPI_Allreduce(local_counts, global_counts, 6, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    global.request_cell_count = global_counts[0];
    global.positive_request_cell_count = global_counts[1];
    global.volume_target_feasible_count = global_counts[2];
    global.volume_target_failed_count = global_counts[3];
    global.eligible_target_feasible_count = global_counts[4];
    global.eligible_target_failed_count = global_counts[5];

    // 6) velocity bins: two length-Nvmu SUM reductions (allocated only when
    // the audit is enabled), then rank 0 emits the non-zero bins sorted by
    // (iv, imu).
    if (audit.enabled) {
        std::vector<unsigned long long> local_count(
            static_cast<size_t>(Param::Nvmu), 0ULL);
        std::vector<double> local_number(static_cast<size_t>(Param::Nvmu),
                                         0.0);
        for (size_t b = 0; b < audit.velocity_bins.size(); ++b) {
            const BulkTailVelocityBinAudit& bin = audit.velocity_bins[b];
            const size_t key = static_cast<size_t>(bin.iv) * Param::Nmu +
                               static_cast<size_t>(bin.imu);
            local_count[key] = bin.request_cell_count;
            local_number[key] = bin.request_number;
        }
        std::vector<unsigned long long> global_count(
            static_cast<size_t>(Param::Nvmu), 0ULL);
        std::vector<double> global_number(static_cast<size_t>(Param::Nvmu),
                                          0.0);
        MPI_Allreduce(local_count.data(), global_count.data(), Param::Nvmu,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local_number.data(), global_number.data(), Param::Nvmu,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t key = static_cast<size_t>(iv) * Param::Nmu +
                                   static_cast<size_t>(imu);
                if (global_count[key] == 0) continue;
                BulkTailVelocityBinAudit bin;
                bin.iv = iv;
                bin.imu = imu;
                bin.request_cell_count = global_count[key];
                bin.request_number = global_number[key];
                global.velocity_bins.push_back(bin);
            }
        }
    }

    // 7) top cells: every rank keeps at most top_cell_limit local
    // candidates; rank 0 gathers, sorts deterministically by
    // (score desc, ix_global, iv, imu) and truncates again to
    // top_cell_limit (section 7.11.16B item 4: the second global truncation
    // is mandatory).
    const int local_top_count =
        static_cast<int>(std::min(audit.top_cells.size(), top_cell_limit));
    int top_slots = 0;
    MPI_Allreduce(&local_top_count, &top_slots, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    const int kTopI = 4;
    const int kTopD = 1 + 2 * BULK_TAIL_MOMENT_COUNT;
    std::vector<int> local_top_i(static_cast<size_t>(kTopI * top_slots), -1);
    std::vector<double> local_top_d(static_cast<size_t>(kTopD * top_slots),
                                    0.0);
    for (int q = 0; q < local_top_count; ++q) {
        const BulkTailMomentAuditTopCell& cell = audit.top_cells[static_cast<size_t>(q)];
        local_top_i[static_cast<size_t>(kTopI * q)] = cell.ix_global;
        local_top_i[static_cast<size_t>(kTopI * q + 1)] = cell.iv;
        local_top_i[static_cast<size_t>(kTopI * q + 2)] = cell.imu;
        local_top_i[static_cast<size_t>(kTopI * q + 3)] = cell.rank;
        local_top_d[static_cast<size_t>(kTopD * q)] = cell.score;
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            local_top_d[static_cast<size_t>(kTopD * q + 1 + m)] =
                cell.center[m];
            local_top_d[static_cast<size_t>(kTopD * q + 1 +
                                           BULK_TAIL_MOMENT_COUNT + m)] =
                cell.volume[m];
        }
    }
    std::vector<int> gathered_top_i;
    std::vector<double> gathered_top_d;
    if (mpi_rank == 0) {
        gathered_top_i.resize(
            static_cast<size_t>(kTopI * top_slots * mpi_size));
        gathered_top_d.resize(
            static_cast<size_t>(kTopD * top_slots * mpi_size));
    }
    MPI_Gather(local_top_i.data(), kTopI * top_slots, MPI_INT,
               mpi_rank == 0 ? gathered_top_i.data() : NULL, kTopI * top_slots,
               MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(local_top_d.data(), kTopD * top_slots, MPI_DOUBLE,
               mpi_rank == 0 ? gathered_top_d.data() : NULL, kTopD * top_slots,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        for (int q = 0; q < top_slots * mpi_size; ++q) {
            const size_t io = static_cast<size_t>(kTopI * q);
            const size_t od = static_cast<size_t>(kTopD * q);
            if (gathered_top_i[io] < 0) continue;
            BulkTailMomentAuditTopCell cell;
            cell.ix_global = gathered_top_i[io];
            cell.iv = gathered_top_i[io + 1];
            cell.imu = gathered_top_i[io + 2];
            cell.rank = gathered_top_i[io + 3];
            cell.score = gathered_top_d[od];
            for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
                cell.center[m] = gathered_top_d[od + 1 + m];
                cell.volume[m] =
                    gathered_top_d[od + 1 + BULK_TAIL_MOMENT_COUNT + m];
            }
            global.top_cells.push_back(cell);
        }
        std::sort(global.top_cells.begin(), global.top_cells.end(),
                  [](const BulkTailMomentAuditTopCell& a,
                     const BulkTailMomentAuditTopCell& b) {
                      if (a.score != b.score) return a.score > b.score;
                      if (a.ix_global != b.ix_global)
                          return a.ix_global < b.ix_global;
                      if (a.iv != b.iv) return a.iv < b.iv;
                      return a.imu < b.imu;
                  });
        if (global.top_cells.size() > top_cell_limit) {
            global.top_cells.resize(top_cell_limit);
        }
    }
    return BulkTailMomentAuditIoStatus::OK;
}
