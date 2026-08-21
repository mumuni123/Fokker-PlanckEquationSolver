#include "bulk_tail_moment_audit.h"

#include "tail_moment_constraint.h"
#include "tail_subcell_quadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
void add_moments(std::array<double, BULK_TAIL_MOMENT_COUNT>& dst,
                 const std::array<double, BULK_TAIL_MOMENT_COUNT>& src,
                 double scale)
{
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) dst[m] += scale * src[m];
}

std::array<double, BULK_TAIL_MOMENT_COUNT> unit_moments(double upar,
                                                         double uperp)
{
    std::array<double, BULK_TAIL_MOMENT_COUNT> r = {};
    mass_cell_moments(1.0, upar, uperp, r[0], r[1], r[3], r[2], r[4], r[5]);
    return r;
}

bool feasible(const std::vector<std::array<double, BULK_TAIL_MOMENT_COUNT> >& cols,
              const std::array<double, BULK_TAIL_MOMENT_COUNT>& target,
              const std::vector<double>& prior)
{
    std::vector<std::vector<double> > matrix(cols.size(),
                                               std::vector<double>(BULK_TAIL_MOMENT_COUNT));
    std::vector<double> ref(BULK_TAIL_MOMENT_COUNT, 0.0), result;
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) ref[m] = target[m];
    for (size_t q = 0; q < cols.size(); ++q)
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) matrix[q][m] = cols[q][m];
    return tail_solve_nonnegative_moment_weights(matrix, ref, prior, result, 1.0e-10);
}

// Section 7.11.16B item 5: relative ratio r(a,b):
//   a/b            if b > 0,
//   0              if a == 0 and b == 0,
//   +inf           if a > 0 and b == 0.
// `defined` is 0 exactly when the ratio is +inf.
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
}

BulkTailMomentAuditTopCell::BulkTailMomentAuditTopCell()
    : ix_global(-1), iv(-1), imu(-1), rank(-1), score(0.0), center(), volume(),
      relative_defined()
{}

MomentRepresentationAudit::MomentRepresentationAudit()
    : enabled(true), finite(true), request_cell_count(0), positive_request_cell_count(0),
      volume_target_feasible_count(0), volume_target_failed_count(0),
      eligible_target_feasible_count(0), eligible_target_failed_count(0),
      eligible_number_fraction(0.0), below_threshold_number_fraction(0.0),
      threshold_window_number(0.0), center(), volume(), eligible_raw(),
      eligible_normalized(), delta_signed(), delta_l1(), center_l1(),
      max_cell_relative(), relative_defined(), velocity_bins(), top_cells()
{}

MomentRepresentationAudit bulk_tail_audit_conversion_requests(
    const CylindricalVelocityGrid& cgrid,
    const HybridVelocityPartition& partition,
    const std::vector<BulkTailMomentAuditRequest>& requests, int rank,
    size_t top_cell_count)
{
    MomentRepresentationAudit audit;
    audit.request_cell_count = requests.size();
    std::array<long double, BULK_TAIL_MOMENT_COUNT> center_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> volume_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> eligible_raw_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> eligible_norm_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> delta_signed_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> delta_l1_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> center_l1_sum = {};
    // Section 7.11.16B item 7: aggregate the real requests by velocity cell
    // (iv, imu) with key = iv*Param::Nmu + imu; only positive-mass cells.
    std::vector<size_t> bin_map(static_cast<size_t>(Param::Nvmu),
                                static_cast<size_t>(-1));
    std::vector<BulkTailVelocityBinAudit> bins;
    bins.reserve(requests.size());
    for (size_t c = 0; c < requests.size(); ++c) {
        const BulkTailMomentAuditRequest& request = requests[c];
        if (!std::isfinite(request.mass)) {
            audit.finite = false;
            continue;
        }
        if (!(request.mass > 0.0)) continue;
        if (request.iv < 0 || request.iv >= Param::Nv ||
            request.imu < 0 || request.imu >= Param::Nmu) {
            audit.finite = false;
            continue;
        }
        ++audit.positive_request_cell_count;
        const size_t key = static_cast<size_t>(request.iv) * Param::Nmu +
                           static_cast<size_t>(request.imu);
        if (bin_map[key] == static_cast<size_t>(-1)) {
            bin_map[key] = bins.size();
            BulkTailVelocityBinAudit bin;
            bin.iv = request.iv;
            bin.imu = request.imu;
            bins.push_back(bin);
        }
        BulkTailVelocityBinAudit& bin = bins[bin_map[key]];
        ++bin.request_cell_count;
        bin.request_number += request.mass;
        const std::array<double, BULK_TAIL_MOMENT_COUNT> center =
            unit_moments(cgrid.upar_cells[request.iv], cgrid.uperp_cells[request.imu]);
        std::vector<TailSubcellNode> nodes = TailSubcellQuadrature::nodes(
            cgrid, request.iv, request.imu);
        std::array<double, BULK_TAIL_MOMENT_COUNT> volume = {};
        std::array<double, BULK_TAIL_MOMENT_COUNT> eligible = {};
        std::vector<std::array<double, BULK_TAIL_MOMENT_COUNT> > all_columns;
        std::vector<double> all_prior;
        std::vector<std::array<double, BULK_TAIL_MOMENT_COUNT> > eligible_columns;
        std::vector<double> eligible_prior;
        double eligible_mass = 0.0;
        for (size_t q = 0; q < nodes.size(); ++q) {
            const std::array<double, BULK_TAIL_MOMENT_COUNT> column =
                unit_moments(nodes[q].upar, nodes[q].uperp);
            const double weight = request.mass * nodes[q].mass_fraction;
            add_moments(volume, column, weight);
            all_columns.push_back(column); all_prior.push_back(weight);
            if (nodes[q].kinetic_energy >= partition.min_conversion_energy) {
                add_moments(eligible, column, weight);
                eligible_columns.push_back(column); eligible_prior.push_back(weight);
                eligible_mass += weight;
                if (nodes[q].kinetic_energy <= partition.min_conversion_energy +
                    0.2e6 * Const::eV) audit.threshold_window_number += weight;
            }
        }
        std::array<double, BULK_TAIL_MOMENT_COUNT> integrated_center = {};
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
            integrated_center[m] = request.mass * center[m];
        add_moments(audit.volume, volume, 1.0);
        add_moments(audit.eligible_raw, eligible, 1.0);
        if (eligible_mass > 0.0) {
            add_moments(audit.eligible_normalized, eligible,
                        request.mass / eligible_mass);
        }
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            const long double center_value = integrated_center[m];
            const long double delta =
                static_cast<long double>(volume[m]) - center_value;
            center_sum[m] += center_value;
            volume_sum[m] += volume[m];
            eligible_raw_sum[m] += eligible[m];
            eligible_norm_sum[m] += (eligible_mass > 0.0)
                ? static_cast<long double>(eligible[m]) * request.mass / eligible_mass
                : 0.0L;
            delta_signed_sum[m] += delta;
            delta_l1_sum[m] += std::fabs(delta);
            center_l1_sum[m] += std::fabs(center_value);
        }
        if (feasible(all_columns, integrated_center, all_prior))
            ++audit.volume_target_feasible_count;
        else ++audit.volume_target_failed_count;
        if (feasible(eligible_columns, integrated_center, eligible_prior))
            ++audit.eligible_target_feasible_count;
        else ++audit.eligible_target_failed_count;
        BulkTailMomentAuditTopCell top;
        top.ix_global = request.ix_global; top.iv = request.iv; top.imu = request.imu;
        top.rank = rank; top.center = integrated_center; top.volume = volume;
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            const double delta = volume[m] - integrated_center[m];
            unsigned char defined = 1;
            const double r = relative_ratio(
                std::fabs(delta), std::fabs(integrated_center[m]), defined);
            top.relative_defined[m] = defined;
            top.score += r;
            audit.max_cell_relative[m] =
                std::max(audit.max_cell_relative[m], r);
        }
        audit.top_cells.push_back(top);
        // Section 7.11.16B item 5: an undefined (inf) relative ratio is a
        // diagnostic, not a physics failure; finite only checks the raw
        // moments.
        bool cell_finite = true;
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            cell_finite = cell_finite &&
                          std::isfinite(integrated_center[m]) &&
                          std::isfinite(volume[m]);
        }
        audit.finite = audit.finite && cell_finite;
    }
    std::sort(bins.begin(), bins.end(),
              [](const BulkTailVelocityBinAudit& a,
                 const BulkTailVelocityBinAudit& b) {
                  const size_t ka = static_cast<size_t>(a.iv) * Param::Nmu +
                                    static_cast<size_t>(a.imu);
                  const size_t kb = static_cast<size_t>(b.iv) * Param::Nmu +
                                    static_cast<size_t>(b.imu);
                  return ka < kb;
              });
    audit.velocity_bins = bins;
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        audit.center[m] = static_cast<double>(center_sum[m]);
        audit.volume[m] = static_cast<double>(volume_sum[m]);
        audit.eligible_raw[m] = static_cast<double>(eligible_raw_sum[m]);
        audit.eligible_normalized[m] = static_cast<double>(eligible_norm_sum[m]);
        audit.delta_signed[m] = static_cast<double>(delta_signed_sum[m]);
        audit.delta_l1[m] = static_cast<double>(delta_l1_sum[m]);
        audit.center_l1[m] = static_cast<double>(center_l1_sum[m]);
    }
    const double total = audit.center[BULK_TAIL_MOMENT_N];
    audit.eligible_number_fraction = total > 0.0
        ? audit.eligible_raw[BULK_TAIL_MOMENT_N] / total : 0.0;
    audit.below_threshold_number_fraction = total > 0.0
        ? std::max(0.0, 1.0 - audit.eligible_number_fraction) : 0.0;
    // Event-level relative-defined flags (r rule on the cumulative L1s).
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        unsigned char defined = 1;
        relative_ratio(audit.delta_l1[m], audit.center_l1[m], defined);
        audit.relative_defined[m] = defined;
    }
    std::sort(audit.top_cells.begin(), audit.top_cells.end(),
              [](const BulkTailMomentAuditTopCell& a, const BulkTailMomentAuditTopCell& b) {
                  if (a.score != b.score) return a.score > b.score;
                  if (a.ix_global != b.ix_global) return a.ix_global < b.ix_global;
                  if (a.iv != b.iv) return a.iv < b.iv;
                  return a.imu < b.imu;
              });
    if (audit.top_cells.size() > top_cell_count) audit.top_cells.resize(top_cell_count);
    return audit;
}

bool bulk_tail_moment_audit_event_is_writable(
    bool step_accepted, const MomentRepresentationAudit& audit)
{
    return step_accepted && audit.enabled;
}
