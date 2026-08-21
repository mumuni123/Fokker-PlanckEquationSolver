#ifndef BULK_TAIL_MOMENT_AUDIT_IO_H
#define BULK_TAIL_MOMENT_AUDIT_IO_H

#include "bulk_tail_moment_audit.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Section 7.11.16B item 3: the accepted-step global reduction is a
// separately testable production function, called both by
// VpfpDiagnostics::write_bulk_tail_moment_audit_accepted_step and by the
// real-moment-audit regression tests.  Tests never re-implement the MPI
// formulas.

enum class BulkTailMomentAuditIoStatus {
    OK = 0,
    INVALID_EVENT_LAYOUT = 1
};

struct BulkTailMomentAuditGlobal {
    bool has_positive_requests;
    std::uint64_t request_cell_count;
    std::uint64_t positive_request_cell_count;
    std::uint64_t volume_target_feasible_count;
    std::uint64_t volume_target_failed_count;
    std::uint64_t eligible_target_feasible_count;
    std::uint64_t eligible_target_failed_count;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center;
    std::array<double, BULK_TAIL_MOMENT_COUNT> volume;
    std::array<double, BULK_TAIL_MOMENT_COUNT> eligible_raw;
    std::array<double, BULK_TAIL_MOMENT_COUNT> eligible_normalized;
    std::array<double, BULK_TAIL_MOMENT_COUNT> delta_signed;
    std::array<double, BULK_TAIL_MOMENT_COUNT> delta_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> max_cell_relative;
    std::array<unsigned char, BULK_TAIL_MOMENT_COUNT> relative_defined;
    double eligible_number_fraction;
    double threshold_window_number;
    std::vector<BulkTailVelocityBinAudit> velocity_bins;
    std::vector<BulkTailMomentAuditTopCell> top_cells;
    BulkTailMomentAuditGlobal();
};

// Section 7.11.16B item 4: before the per-event loop, MIN/MAX reduce the
// local event count; a mismatch returns INVALID_EVENT_LAYOUT on every rank
// (no rank may proceed into a per-event collective).
BulkTailMomentAuditIoStatus bulk_tail_moment_audit_check_event_layout(
    size_t local_event_count, int mpi_rank, int mpi_size);

// Section 7.11.16B item 4: reduce one event's local audit across ranks in
// the fixed collective order (location MIN/MAX check, positive-count SUM,
// raw moments/L1 SUM, per-component cell-max MAX, velocity bins, top cells).
// When the global positive-request count is zero the function performs no
// further collective and returns OK with has_positive_requests=false.
BulkTailMomentAuditIoStatus bulk_tail_moment_audit_reduce_event(
    const MomentRepresentationAudit& audit, int conversion_location,
    size_t top_cell_limit, int mpi_rank, int mpi_size,
    BulkTailMomentAuditGlobal& global);

#endif
