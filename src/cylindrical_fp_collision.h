#ifndef CYLINDRICAL_FP_COLLISION_H
#define CYLINDRICAL_FP_COLLISION_H

#include "collision_coefficients.h"
#include "species.h"

#include <cstddef>
#include <vector>

enum class CollisionIntegratorType {
    BACKWARD_EULER,
    TR_BDF2
};

enum class BulkCollisionIntegrator {
    BGK_VALIDATION,
    CHANG_COOPER_FLUX
};

// Integrated conservative transfers for one collision half-step.  A positive
// face value is directed toward increasing velocity.  Values already include
// the collision substep dt and the cylindrical face area, so they can be used
// directly for FV residuals and for bulk-to-tail parcel export.
struct CollisionFaceFluxes {
    int nx_local;
    int nv;
    int nmu;
    std::vector<double> upar_flux;
    std::vector<double> uperp_flux;
    std::vector<double> cross_upar_flux;
    std::vector<double> cross_uperp_flux;

    CollisionFaceFluxes()
        : nx_local(0), nv(0), nmu(0) {}

    void resize(int nx_local_value, int nv_value, int nmu_value)
    {
        nx_local = nx_local_value;
        nv = nv_value;
        nmu = nmu_value;
        const std::size_t nx = static_cast<std::size_t>(nx_local);
        const std::size_t nvsz = static_cast<std::size_t>(nv);
        const std::size_t nmusz = static_cast<std::size_t>(nmu);
        upar_flux.assign(nx * (nvsz + 1) * nmusz, 0.0);
        uperp_flux.assign(nx * nvsz * (nmusz + 1), 0.0);
        cross_upar_flux.assign(nx * (nvsz + 1) * nmusz, 0.0);
        cross_uperp_flux.assign(nx * nvsz * (nmusz + 1), 0.0);
    }

    std::size_t upar_index(int ix, int face, int imu) const
    {
        return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(nv + 1) +
                static_cast<std::size_t>(face)) *
               static_cast<std::size_t>(nmu) + static_cast<std::size_t>(imu);
    }

    std::size_t uperp_index(int ix, int iv, int face) const
    {
        return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(nv) +
                static_cast<std::size_t>(iv)) *
               static_cast<std::size_t>(nmu + 1) + static_cast<std::size_t>(face);
    }
};

struct CollisionDiagnostics {
    bool success;
    bool order_reduced;
    bool unsupported_cross_diffusion;
    double mass_change;
    double reservoir_energy_change;
    double reservoir_momentum_change;
    double tail_mass;
    double interface_export_number;
    double interface_export_energy;
    double implicit_flux_residual_linf;
    double cross_flux_pair_residual_linf;
    double flux_divergence_sum;
    double mass_flux_balance_residual;
    double interface_inward_clipped_number;
    std::size_t interface_parcel_count;
    int transaction_rollback_count;
    CollisionDiagnostics()
        : success(true), order_reduced(false), unsupported_cross_diffusion(false),
          mass_change(0.0), reservoir_energy_change(0.0),
          reservoir_momentum_change(0.0), tail_mass(0.0),
          interface_export_number(0.0), interface_export_energy(0.0),
          implicit_flux_residual_linf(0.0),
          cross_flux_pair_residual_linf(0.0),
          flux_divergence_sum(0.0),
          mass_flux_balance_residual(0.0),
          interface_inward_clipped_number(0.0), interface_parcel_count(0),
          transaction_rollback_count(0) {}
};

// A deliberately separate cylindrical collision implementation.  It never
// calls legacy collision.cpp, whose spherical-grid coefficients have different
// units and conservation semantics.
class CylindricalFokkerPlanckCollision {
public:
    CylindricalFokkerPlanckCollision(const CollisionCoefficientProvider& provider,
                                     CollisionIntegratorType type);
    // bulk_mask (optional): per velocity-slot flag (1 = the bulk collision
    // owns the cell, 0 = excluded).  Excluded cells (the bulk->tail
    // conversion region, section 19.3) are left unchanged and act as
    // zero-flux walls: the Chang-Cooper sweeps must not diffuse across the
    // sharp conversion boundary, which otherwise produces negative
    // overshoot and fails the positivity check (H9 beam-40fs failure).  A
    // NULL mask means every cell is owned by the bulk.
    CollisionDiagnostics apply(Species& electrons, const SpatialGrid& grid,
                               double time, double dt,
                               const std::vector<unsigned char>* bulk_mask = NULL) const;
    CollisionDiagnostics apply_with_flux(
        Species& electrons, const SpatialGrid& grid, double time, double dt,
        const std::vector<unsigned char>* bulk_mask,
        const HybridVelocityPartition* partition,
        CollisionFaceFluxes& fluxes) const;
    void set_bulk_integrator(BulkCollisionIntegrator integrator)
    { bulk_integrator_ = integrator; }
    BulkCollisionIntegrator bulk_integrator() const
    { return bulk_integrator_; }
    // Stage-H8: the same provider feeds the HybridCollisionStep (identical
    // time-layer coefficients for bulk and tail, section 10.4).
    const CollisionCoefficientProvider& provider() const { return provider_; }
    // True when the provider is the collisionless path (collision=none), so
    // the integrator can skip the trial copy for the first Strang half
    // (section 14.7 item 8: with a real collision the first half must act on
    // a trial buffer, never on the accepted species).
    bool is_trivial() const;

private:
    const CollisionCoefficientProvider& provider_;
    CollisionIntegratorType type_;
    BulkCollisionIntegrator bulk_integrator_;
};

#endif
