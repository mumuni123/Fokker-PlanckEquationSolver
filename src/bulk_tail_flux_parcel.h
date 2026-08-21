#ifndef BULK_TAIL_FLUX_PARCEL_H
#define BULK_TAIL_FLUX_PARCEL_H

#include "grid.h"

#include <cstddef>
#include <vector>

// A positive phase-space mass parcel exported by one final bulk face flux.
// The mass has the same units as Species::f: electrons per transverse area;
// it does not include dx.  Nodes retain the swept interval information that
// is lost by a scalar face flux.
struct FluxParcelNode {
    double upar;
    double uperp;
    double mass;
};

struct BulkTailFluxParcel {
    int ix_local;
    int ix_global;
    VelocityFaceDirection direction;
    int face_index;
    int transverse_index;
    int operator_stage;
    // Positive amount exported by the originating final interface face.
    double face_number;
    std::vector<FluxParcelNode> nodes;
    // Read-only u_parallel audit provenance.  These values never enter the
    // physical tail conversion, but identify a face/parcel mismatch without
    // retaining every quadrature node in diagnostic output.
    int audit_node_failure_reason;
    double audit_reconstructed_target;
    double audit_node_sum;
    double number;
    double px;
    double jx_dx;
    double kinetic_energy;
    double pixx_dx;
    double piperp_dx;

    BulkTailFluxParcel();
    bool finite_nonnegative() const;
    void recompute_moments();
};

struct BulkTailFluxFaceAudit {
    int ix_global;
    int direction;
    int face_index;
    int transverse_index;
    int operator_stage;
    double face_number;
    double parcel_number;
    int node_failure_reason;
    double reconstructed_target;
    double node_sum;
};

struct BulkTailFluxBatch {
    std::vector<BulkTailFluxParcel> parcels;
    // Exporting a batch and applying the one-way interface sink are separate
    // operations.  flux-audit leaves this false so the production state is
    // unchanged; flux-interface sets it true before the remap call.
    bool apply_interface_sink;
    bool finite;
    bool nonnegative;
    double quadrature_error_max;
    std::size_t duplicate_count;
    double below_threshold_number;
    std::size_t face_audit_count;
    double face_audit_face_abs_sum;
    double face_audit_parcel_abs_sum;
    double face_audit_abs_error_sum;
    double face_audit_max_relative;
    double face_audit_abs_at_max_relative;
    bool face_audit_max_valid;
    BulkTailFluxFaceAudit face_audit_max;

    BulkTailFluxBatch();
    void clear();
    void recompute(double threshold_energy);
};

// Adds one node and updates all six moments through the shared production
// mass_cell_moments() formula.  The function rejects non-finite or negative
// mass instead of silently taking an absolute value.
bool bulk_tail_parcel_add_node(BulkTailFluxParcel& parcel,
                               double upar, double uperp, double mass);

// Six-moment difference between a batch and a supplied reference.  The order
// is N, Px, Jx, K, Pixx and Piperp.
void bulk_tail_batch_moments(const BulkTailFluxBatch& batch,
                             double out[6]);

#endif
