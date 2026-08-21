#ifndef BULK_TAIL_MOMENT_AUDIT_H
#define BULK_TAIL_MOMENT_AUDIT_H

#include "grid.h"

#include <array>
#include <cstdint>
#include <vector>

enum BulkTailMomentComponent {
    BULK_TAIL_MOMENT_N = 0,
    BULK_TAIL_MOMENT_PX = 1,
    BULK_TAIL_MOMENT_JX = 2,
    BULK_TAIL_MOMENT_K = 3,
    BULK_TAIL_MOMENT_PIXX = 4,
    BULK_TAIL_MOMENT_PIPERP = 5,
    BULK_TAIL_MOMENT_COUNT = 6
};

struct BulkTailMomentAuditRequest {
    int ix_global;
    int iv;
    int imu;
    double mass;
};

struct BulkTailMomentAuditTopCell {
    int ix_global;
    int iv;
    int imu;
    int rank;
    double score;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center;
    std::array<double, BULK_TAIL_MOMENT_COUNT> volume;
    // Section 7.11.16B: per-component relative ratio defined flag (1 when
    // the ratio |delta|/|center| is defined, 0 when it is +inf).
    std::array<unsigned char, BULK_TAIL_MOMENT_COUNT> relative_defined;
    BulkTailMomentAuditTopCell();
};

// Section 7.11.16B: one velocity cell (iv, imu) of the real-event velocity
// histogram.  request_number is the sum of request masses (m^-2), not the
// PIC macro-particle count.
struct BulkTailVelocityBinAudit {
    int iv;
    int imu;
    std::uint64_t request_cell_count;
    double request_number;
    BulkTailVelocityBinAudit()
        : iv(-1), imu(-1), request_cell_count(0), request_number(0.0)
    {}
};

// Read-only event diagnostic.  Every sum is local until the diagnostics
// layer performs its accepted-step MPI reduction.
struct MomentRepresentationAudit {
    bool enabled;
    bool finite;
    std::uint64_t request_cell_count;
    std::uint64_t positive_request_cell_count;
    std::uint64_t volume_target_feasible_count;
    std::uint64_t volume_target_failed_count;
    std::uint64_t eligible_target_feasible_count;
    std::uint64_t eligible_target_failed_count;
    double eligible_number_fraction;
    double below_threshold_number_fraction;
    double threshold_window_number;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center;
    std::array<double, BULK_TAIL_MOMENT_COUNT> volume;
    std::array<double, BULK_TAIL_MOMENT_COUNT> eligible_raw;
    std::array<double, BULK_TAIL_MOMENT_COUNT> eligible_normalized;
    std::array<double, BULK_TAIL_MOMENT_COUNT> delta_signed;
    std::array<double, BULK_TAIL_MOMENT_COUNT> delta_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> max_cell_relative;
    // Section 7.11.16B: event-level relative-ratio defined flag per
    // component (r(a,b) rule; 0 only when the denominator is zero with a
    // nonzero numerator, i.e. the ratio is +inf).
    std::array<unsigned char, BULK_TAIL_MOMENT_COUNT> relative_defined;
    std::vector<BulkTailVelocityBinAudit> velocity_bins;
    std::vector<BulkTailMomentAuditTopCell> top_cells;
    MomentRepresentationAudit();
};

MomentRepresentationAudit bulk_tail_audit_conversion_requests(
    const CylindricalVelocityGrid& cgrid,
    const HybridVelocityPartition& partition,
    const std::vector<BulkTailMomentAuditRequest>& requests, int rank,
    size_t top_cell_count);

// Shared accepted-state gate used by the production diagnostics and its
// independent regression test.  Trial/rejected states must never reach the
// accepted-event output path.
bool bulk_tail_moment_audit_event_is_writable(
    bool step_accepted, const MomentRepresentationAudit& audit);

#endif
