#include "bulk_tail_flux_parcel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {
void zero_moments(double a[6])
{
    for (int i = 0; i < 6; ++i) a[i] = 0.0;
}
}

BulkTailFluxParcel::BulkTailFluxParcel()
    : ix_local(-1), ix_global(-1),
      direction(VelocityFaceDirection::U_PARALLEL), face_index(-1),
      transverse_index(-1), operator_stage(0), face_number(0.0),
      audit_node_failure_reason(0),
      audit_reconstructed_target(std::numeric_limits<double>::quiet_NaN()),
      audit_node_sum(std::numeric_limits<double>::quiet_NaN()),
      number(0.0), px(0.0),
      jx_dx(0.0), kinetic_energy(0.0), pixx_dx(0.0), piperp_dx(0.0)
{}

bool BulkTailFluxParcel::finite_nonnegative() const
{
    const double values[7] = { face_number, number, px, jx_dx,
                               kinetic_energy, pixx_dx, piperp_dx };
    for (int i = 0; i < 7; ++i) {
        if (!std::isfinite(values[i])) return false;
    }
    if (face_number < 0.0) return false;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (!std::isfinite(nodes[i].upar) || !std::isfinite(nodes[i].uperp) ||
            !std::isfinite(nodes[i].mass) || nodes[i].mass < 0.0) return false;
    }
    return true;
}

void BulkTailFluxParcel::recompute_moments()
{
    number = px = jx_dx = kinetic_energy = pixx_dx = piperp_dx = 0.0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double n = 0.0, p = 0.0, e = 0.0, j = 0.0, xx = 0.0, pp = 0.0;
        mass_cell_moments(nodes[i].mass, nodes[i].upar, nodes[i].uperp,
                          n, p, e, j, xx, pp);
        number += n;
        px += p;
        kinetic_energy += e;
        jx_dx += j;
        pixx_dx += xx;
        piperp_dx += pp;
    }
}

BulkTailFluxBatch::BulkTailFluxBatch()
    : apply_interface_sink(false), finite(true), nonnegative(true),
      quadrature_error_max(0.0),
      duplicate_count(0), below_threshold_number(0.0), face_audit_count(0),
      face_audit_face_abs_sum(0.0), face_audit_parcel_abs_sum(0.0),
      face_audit_abs_error_sum(0.0), face_audit_max_relative(0.0),
      face_audit_abs_at_max_relative(0.0), face_audit_max_valid(false)
{}

void BulkTailFluxBatch::clear()
{
    // The sink policy is part of the caller's step configuration, not batch
    // data.  Preserve it while discarding the parcels from a previous trial;
    // otherwise flux-interface silently degrades to flux-audit.
    const bool keep_interface_sink = apply_interface_sink;
    parcels.clear();
    apply_interface_sink = keep_interface_sink;
    finite = true;
    nonnegative = true;
    quadrature_error_max = 0.0;
    duplicate_count = 0;
    below_threshold_number = 0.0;
    face_audit_count = 0;
    face_audit_face_abs_sum = 0.0;
    face_audit_parcel_abs_sum = 0.0;
    face_audit_abs_error_sum = 0.0;
    face_audit_max_relative = 0.0;
    face_audit_abs_at_max_relative = 0.0;
    face_audit_max_valid = false;
}

void BulkTailFluxBatch::recompute(double threshold_energy)
{
    finite = true;
    nonnegative = true;
    quadrature_error_max = 0.0;
    duplicate_count = 0;
    below_threshold_number = 0.0;
    face_audit_count = 0;
    face_audit_face_abs_sum = 0.0;
    face_audit_parcel_abs_sum = 0.0;
    face_audit_abs_error_sum = 0.0;
    face_audit_max_relative = 0.0;
    face_audit_abs_at_max_relative = 0.0;
    face_audit_max_valid = false;
    std::vector<std::uint64_t> keys;
    keys.reserve(parcels.size());
    for (size_t i = 0; i < parcels.size(); ++i) {
        parcels[i].recompute_moments();
        const bool parcel_ok = parcels[i].finite_nonnegative();
        finite = finite && parcel_ok;
        if (!parcel_ok) nonnegative = false;
        BulkTailFluxFaceAudit audit;
        audit.ix_global = parcels[i].ix_global;
        audit.direction = parcels[i].direction ==
            VelocityFaceDirection::U_PARALLEL ? 0 : 1;
        audit.face_index = parcels[i].face_index;
        audit.transverse_index = parcels[i].transverse_index;
        audit.operator_stage = parcels[i].operator_stage;
        // Keep the originating final face transfer distinct from the
        // quadrature sum; otherwise this audit would become tautological.
        audit.face_number = parcels[i].face_number;
        audit.parcel_number = parcels[i].number;
        audit.node_failure_reason = parcels[i].audit_node_failure_reason;
        audit.reconstructed_target = parcels[i].audit_reconstructed_target;
        audit.node_sum = parcels[i].audit_node_sum;
        const double absolute_error = std::fabs(
            audit.face_number - audit.parcel_number);
        const double relative_error = absolute_error / std::max(
            std::max(std::fabs(audit.face_number),
                     std::fabs(audit.parcel_number)), 1.0e-300);
        ++face_audit_count;
        face_audit_face_abs_sum += std::fabs(audit.face_number);
        face_audit_parcel_abs_sum += std::fabs(audit.parcel_number);
        face_audit_abs_error_sum += absolute_error;
        if (!face_audit_max_valid || relative_error > face_audit_max_relative) {
            face_audit_max_valid = true;
            face_audit_max_relative = relative_error;
            face_audit_abs_at_max_relative = absolute_error;
            face_audit_max = audit;
        }
        const std::uint64_t direction =
            parcels[i].direction == VelocityFaceDirection::U_PARALLEL ? 0ULL : 1ULL;
        // The key is the physical interface event, not the support node.  A
        // duplicate face export would double-remove bulk mass.
        std::uint64_t key = direction;
        key = key * 1000003ULL + static_cast<std::uint64_t>(parcels[i].ix_global + 1);
        key = key * 1000003ULL + static_cast<std::uint64_t>(parcels[i].face_index + 1);
        key = key * 1000003ULL + static_cast<std::uint64_t>(parcels[i].transverse_index + 1);
        key = key * 1000003ULL + static_cast<std::uint64_t>(parcels[i].operator_stage + 1);
        keys.push_back(key);
        for (size_t q = 0; q < parcels[i].nodes.size(); ++q) {
            const double gamma = std::sqrt(
                1.0 + parcels[i].nodes[q].upar * parcels[i].nodes[q].upar +
                parcels[i].nodes[q].uperp * parcels[i].nodes[q].uperp);
            const double ke = Const::me * Const::c * Const::c * (gamma - 1.0);
            if (ke < threshold_energy) below_threshold_number +=
                parcels[i].nodes[q].mass;
        }
    }
    std::sort(keys.begin(), keys.end());
    for (size_t i = 1; i < keys.size(); ++i)
        if (keys[i] == keys[i - 1]) ++duplicate_count;
}

bool bulk_tail_parcel_add_node(BulkTailFluxParcel& parcel,
                               double upar, double uperp, double mass)
{
    if (!std::isfinite(upar) || !std::isfinite(uperp) ||
        !std::isfinite(mass) || mass < 0.0) return false;
    FluxParcelNode node;
    node.upar = upar;
    node.uperp = uperp;
    node.mass = mass;
    parcel.nodes.push_back(node);
    parcel.recompute_moments();
    return true;
}

void bulk_tail_batch_moments(const BulkTailFluxBatch& batch,
                             double out[6])
{
    zero_moments(out);
    for (size_t p = 0; p < batch.parcels.size(); ++p) {
        out[0] += batch.parcels[p].number;
        out[1] += batch.parcels[p].px;
        out[2] += batch.parcels[p].jx_dx;
        out[3] += batch.parcels[p].kinetic_energy;
        out[4] += batch.parcels[p].pixx_dx;
        out[5] += batch.parcels[p].piperp_dx;
    }
}
