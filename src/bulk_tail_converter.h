#ifndef BULK_TAIL_CONVERTER_H
#define BULK_TAIL_CONVERTER_H

#include "background_tail_pic.h"
#include "bulk_tail_moment_audit.h"
#include "bulk_tail_flux_parcel.h"
#include "grid.h"
#include "species.h"
#include "tail_subcell_quadrature.h"

#include <cstdint>
#include <map>
#include <vector>

// Bulk-to-tail conservative converter (sections 7 and 17).  The production
// flux-interface path consumes the final remap/collision face parcels; the
// static extractor remains isolated to explicit A/B modes.
//
// Loading rules (sections 7.3-7.5 and 7.7):
//  * conversion-cell mass is aggregated per local spatial cell by
//    (u_parallel sign, |u_parallel| bin, energy bin);
//  * the golden reference per velocity cell is one azimuthal quartet
//    (section 7.4): phi in {0, pi/2, pi, 3pi/2}, weight M/4;
//  * production loading compresses each group to at most 7 quartet supports
//    with deterministic Caratheodory/null-space elimination (section 7.5);
//  * all created particles are placed at their spatial cell centre so the
//    CIC deposit returns exactly the removed density (section 7.7);
//  * bulk mass is only removed through Species::extract_conversion_masses
//    (section 14.3), never by direct f mutation.

enum class ConversionLocation {
    AFTER_U_SUBSTEP,
    AFTER_COLLISION_HALF
};

enum class TailConversionMode {
    STATIC_CELL,
    FLUX_AUDIT,
    FLUX_INTERFACE
};

inline const char* tail_conversion_mode_name(TailConversionMode mode)
{
    switch (mode) {
        case TailConversionMode::STATIC_CELL: return "static-cell";
        case TailConversionMode::FLUX_AUDIT: return "flux-audit";
        case TailConversionMode::FLUX_INTERFACE: return "flux-interface";
    }
    return "unknown";
}

// Section 7.11.3: test-selectable loading policies for the A/B/C threshold
// interface audit.  The policy only changes how candidate supports are
// compressed within each group; it never changes bulk mass extraction, the
// conversion time layer, spatial position, particle IDs, deposition or
// transaction semantics.  Threshold-aware grouping is the production
// default after the section 7.11.9 cell-volume support gate has passed.
enum class BulkTailLoadingPolicy {
    GOLDEN_QUARTETS_NO_COMPRESSION,
    CURRENT_PRODUCTION_COMPRESSION,
    THRESHOLD_AWARE_COMPRESSION
};

inline const char* bulk_tail_loading_policy_name(
    BulkTailLoadingPolicy policy)
{
    switch (policy) {
        case BulkTailLoadingPolicy::GOLDEN_QUARTETS_NO_COMPRESSION:
            return "golden";
        case BulkTailLoadingPolicy::CURRENT_PRODUCTION_COMPRESSION:
            return "current";
        case BulkTailLoadingPolicy::THRESHOLD_AWARE_COMPRESSION:
            return "threshold-aware";
    }
    return "unknown";
}

struct BulkTailConversionDiagnostics {
    struct SubcellFallbackRecord {
        int ix_global;
        int iv;
        int imu;
        int reason;
        std::uint64_t fallback_particles;
        double number_target;
        double px_target;
        double jx_target;
        double energy_target;
        double pixx_target;
        double piperp_target;

        SubcellFallbackRecord()
            : ix_global(-1), iv(-1), imu(-1), reason(0),
              fallback_particles(0), number_target(0.0), px_target(0.0),
              jx_target(0.0), energy_target(0.0), pixx_target(0.0),
              piperp_target(0.0)
        {}
    };
    int accepted_step;
    ConversionLocation location;
    double number_removed;
    double number_created;
    double px_removed;
    double px_created;
    double energy_removed;
    double energy_created;
    double jx_dx_removed;
    double jx_dx_created;
    double pixx_dx_removed;
    double pixx_dx_created;
    double piperp_dx_removed;
    double piperp_dx_created;
    // L1 scales of the removed contributions (sum of |per-request moments|)
    // used to normalise the residuals: the signed px/jx totals can cancel
    // between the two u_parallel signs, so gates are relative to the removed
    // contribution scale (section 19.1).
    double number_scale;
    double px_scale;
    double energy_scale;
    double jx_scale;
    double pixx_scale;
    double piperp_scale;
    // Relative residuals (created - removed) normalised by the removed scale.
    double number_residual_rel;
    double px_residual_rel;
    double energy_residual_rel;
    double jx_residual_rel;
    double pixx_residual_rel;
    double piperp_residual_rel;
    // Conversion-induced density change: ||n_before - n_after|| /
    // ||n_before|| over the converted local cells (section 7.7).  Cell-centre
    // loading makes this vanish to deposition roundoff.
    double rho_l2_before_after;
    double rho_linf_before_after;
    std::uint64_t particles_created;
    std::uint64_t compression_fallback_count;
    // Hard-contract diagnostics for the flux loader.  A fallback that leaves
    // more than max_supports active supports, a reused particle id, or a
    // parcel whose node moments do not match its originating face transfer
    // is a failed transaction, never a silently accepted approximation.
    std::uint64_t support_limit_violation_count;
    std::uint64_t duplicate_id_count;
    std::uint64_t face_ledger_mismatch_count;
    // Positive face-flux parcels below the grid-scaled floating-point floor
    // are not physical tail transfers.  They are discarded before support
    // compression so subnormal collision residues cannot create a failed
    // conversion transaction.
    double roundoff_discarded_number;
    std::uint64_t subcell_cells_loaded;
    std::uint64_t subcell_support_count;
    std::uint64_t subcell_fallback_count;
    std::vector<SubcellFallbackRecord> subcell_fallbacks;
    MomentRepresentationAudit moment_audit;
    // Threshold-interface source ledger.  All four spectra share these
    // explicit energy edges.  This is an event record, not a snapshot: an
    // accepted step may contain more than one conversion location.
    std::vector<double> conversion_source_energy_edges;
    std::vector<double> pre_extraction_bulk_number_spectrum;
    std::vector<double> pre_extraction_bulk_energy_spectrum;
    std::vector<double> removed_bulk_number_spectrum;
    std::vector<double> removed_bulk_energy_spectrum;
    std::vector<double> created_tail_number_spectrum;
    std::vector<double> created_tail_energy_spectrum;
    // Gate I spatial transfer ledger [m^-2 per local x cell].  These are the
    // exact masses removed/created by this accepted conversion transaction.
    std::vector<double> removed_bulk_number_by_x;
    std::vector<double> created_tail_number_by_x;
    bool finite;
    bool conservative;   // N/Px/K residuals within the hard gate
    bool fidelity_ok;    // Jx/Pixx/Piperp residuals within the fidelity gate
    bool complete;       // extraction applied and post-conversion scan clean

    BulkTailConversionDiagnostics()
        : accepted_step(-1), location(ConversionLocation::AFTER_U_SUBSTEP),
          number_removed(0.0), number_created(0.0),
          px_removed(0.0), px_created(0.0),
          energy_removed(0.0), energy_created(0.0),
          jx_dx_removed(0.0), jx_dx_created(0.0),
          pixx_dx_removed(0.0), pixx_dx_created(0.0),
          piperp_dx_removed(0.0), piperp_dx_created(0.0),
          number_scale(0.0), px_scale(0.0), energy_scale(0.0),
          jx_scale(0.0), pixx_scale(0.0), piperp_scale(0.0),
          number_residual_rel(0.0), px_residual_rel(0.0),
          energy_residual_rel(0.0), jx_residual_rel(0.0),
          pixx_residual_rel(0.0), piperp_residual_rel(0.0),
          rho_l2_before_after(0.0), rho_linf_before_after(0.0),
          particles_created(0), compression_fallback_count(0),
          support_limit_violation_count(0), duplicate_id_count(0),
          face_ledger_mismatch_count(0), roundoff_discarded_number(0.0),
          subcell_cells_loaded(0), subcell_support_count(0),
          subcell_fallback_count(0),
          finite(false), conservative(false), fidelity_ok(false),
          complete(false)
    {
        conversion_source_energy_edges.clear();
        pre_extraction_bulk_number_spectrum.clear();
        pre_extraction_bulk_energy_spectrum.clear();
        removed_bulk_number_spectrum.clear();
        removed_bulk_energy_spectrum.clear();
        created_tail_number_spectrum.clear();
        created_tail_energy_spectrum.clear();
    }
};

class BulkTailConverter {
public:
    BulkTailConverter();

    // Section 7.11.3: test-only loading-policy selector.  Not a compile
    // macro and not a global: the policy is per-converter state.  The
    // production default is CURRENT_PRODUCTION_COMPRESSION.
    void set_loading_policy(BulkTailLoadingPolicy policy)
    {
        loading_policy_ = policy;
    }
    BulkTailLoadingPolicy loading_policy() const { return loading_policy_; }
    // Section 7.11.10 production switch.  It is enabled after the read-only
    // cell-volume audit demonstrates threshold-window support; a nonzero
    // fallback count still keeps the subcell-loading acceptance gate open.
    void set_subcell_loading_enabled(bool enabled)
    { subcell_loading_enabled_ = enabled; }
    bool subcell_loading_enabled() const { return subcell_loading_enabled_; }
    void set_moment_audit_enabled(bool enabled, size_t top_cell_count = 64)
    {
        moment_audit_enabled_ = enabled;
        moment_audit_top_cell_count_ = top_cell_count;
    }
    bool moment_audit_enabled() const { return moment_audit_enabled_; }

    // Section 7.1 suggested interface plus the creating MPI rank (needed by
    // the counter-based tail ID scheme, section 6.4).  Operates on the trial
    // bulk/tail only; on any failure neither object is modified and the
    // caller must discard the trial.
    BulkTailConversionDiagnostics extract_after_substep(
        Species& bulk_trial, BackgroundTailPIC& tail_trial,
        const SpatialGrid& grid, const HybridVelocityPartition& partition,
        int accepted_step, ConversionLocation location, int mpi_rank);

    // Section 7.11.17.5: deterministic loading of parcels emitted by the
    // final u_parallel finite-volume interface flux.  This is transactional:
    // a failed conversion leaves tail_trial, its ID counter and its ledgers
    // unchanged.
    BulkTailConversionDiagnostics convert_flux_batch(
        const BulkTailFluxBatch& batch,
        BackgroundTailPIC& tail_trial,
        const SpatialGrid& grid,
        const HybridVelocityPartition& partition,
        int accepted_step,
        ConversionLocation location,
        int mpi_rank,
        size_t max_supports = 7);

    // One-time initialization for a fresh flux-interface run.  This is not a
    // stepwise static extractor: it packages initial tail-owned bulk masses
    // into the same positive parcel batch and commits bulk/tail only after the
    // normal loader and Species mass-removal transaction both succeed.
    BulkTailConversionDiagnostics convert_initial_tail_cells(
        Species& bulk_trial, BackgroundTailPIC& tail_trial,
        const SpatialGrid& grid, const HybridVelocityPartition& partition,
        int accepted_step, int mpi_rank, size_t max_supports = 7);

private:
    BulkTailLoadingPolicy loading_policy_;
    bool subcell_loading_enabled_;
    bool moment_audit_enabled_;
    size_t moment_audit_top_cell_count_;

    struct Moment6 {
        double n;
        double px;
        double jx;
        double ke;
        double pixx;
        double piperp;
        Moment6()
            : n(0.0), px(0.0), jx(0.0), ke(0.0), pixx(0.0), piperp(0.0)
        {}
    };
    struct CellMassRef {
        int iv;
        int imu;
        double mass;
    };
    struct SubcellGeometry {
        double upar;
        double uperp;
        double mass_fraction;
        Moment6 column;
    };
    struct GroupKey {
        int ix_global;
        int sign;
        int upar_bin;
        int energy_bin;
        bool operator<(const GroupKey& other) const
        {
            if (ix_global != other.ix_global)
                return ix_global < other.ix_global;
            if (sign != other.sign) return sign < other.sign;
            if (upar_bin != other.upar_bin)
                return upar_bin < other.upar_bin;
            return energy_bin < other.energy_bin;
        }
    };

    // Geometry and unit-mass moment columns depend only on the velocity grid
    // and conversion threshold.  Cache them once per (iv, imu), so the hot
    // per-x conversion path only scales columns and solves the local system.
    void prepare_subcell_geometry(const CylindricalVelocityGrid& cgrid,
                                  double min_conversion_energy);
    std::vector<std::vector<SubcellGeometry> > subcell_geometry_;
    double subcell_geometry_threshold_;
    bool subcell_geometry_valid_;

    // Section 7.5: deterministic Caratheodory/null-space elimination that
    // reduces the nonnegative quartet weights to at most max_support active
    // supports while preserving cols^T w = ref.  On failure returns false
    // and leaves w unchanged (the caller keeps the reference quartets).
    static bool compress_quartets(const std::vector<Moment6>& cols,
                                  std::vector<double>& w, const Moment6& ref,
                                  size_t max_support, double tolerance);
};

#endif
