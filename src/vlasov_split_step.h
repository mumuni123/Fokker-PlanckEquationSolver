#ifndef VLASOV_SPLIT_STEP_H
#define VLASOV_SPLIT_STEP_H

#include "conservative_ppm_remap.h"
#include "grid.h"
#include "open_boundary.h"
#include "species.h"

struct XFaceTransportAudit;

// Per-step ledgers of the three collisionless substeps (section 13.7).
struct VlasovStepDiagnostics {
    RemapDiagnostics x_first;
    RemapDiagnostics u_full;
    RemapDiagnostics x_second;
};

// Collisionless background-electron Strang split (sections 6-7, 13.7):
//   T_x(dt/2) -> P[rho^{n+1/2}] -> T_u(E^{n+1/2}, dt) -> T_x(dt/2)
// The caller owns the persistent state buffers (VpfpIntegrator::state_*);
// this class owns only the conservative remap workspace.  Every substep
// writes a fresh output buffer, so a failed step never mutates the accepted
// state; the caller swaps the final buffer in only after validation.
class VlasovSplitStep {
public:
    VlasovSplitStep();

    void init(const SpatialGrid& grid, const Species& prototype,
              const OpenBackgroundBoundary& boundary);

    void set_x_transport_velocity_mode(XTransportVelocityMode mode)
    { remap_.set_x_transport_velocity_mode(mode); }
    XTransportVelocityMode x_transport_velocity_mode() const
    { return remap_.x_transport_velocity_mode(); }

    // T_x(dt/2): state_n -> state_x_half.  Refreshes ghosts and moments.
    // The optional Gate I audit pointer is passed through verbatim to
    // advect_x; the two half-steps must use two different objects (section
    // 4.2).
    bool first_x_half(const Species& state_n, Species& state_x_half,
                      double time, double half_dt,
                      VlasovStepDiagnostics& diag,
                      XFaceTransportAudit* audit = NULL);

    // T_u(E_mid, dt): state_x_half -> state_u_full on the non-uniform
    // u_parallel grid.  Both velocity ends default to zero inflow; outflow
    // mass/energy enter the tail ledger (checked by the caller against the
    // production threshold).  Ghost x-columns are transformed by the same
    // column operator so MPI halos stay consistent; the next T_x regenerates
    // the physical-boundary halos.  local_delta_ke_by_x is the optional
    // Gate I per-cell bulk work (section 4.4), passed through verbatim.
    bool u_full(const Species& state_x_half, Species& state_u_full,
                const EMFields& field_mid, double time_mid, double dt,
                VlasovStepDiagnostics& diag,
                const HybridVelocityPartition* partition = NULL,
                BulkTailFluxBatch* exported_flux = NULL,
                int quadrature_order = 4,
                std::vector<double>* local_delta_ke_by_x = NULL);

    // T_x(dt/2): state_u_full -> state_np1.  Refreshes ghosts and moments.
    bool second_x_half(const Species& state_u_full, Species& state_np1,
                       double time_mid, double half_dt,
                       VlasovStepDiagnostics& diag,
                       XFaceTransportAudit* audit = NULL);

private:
    SpatialGrid grid_;
    const OpenBackgroundBoundary* boundary_;
    ConservativePpmRemap remap_;
    int mpi_rank_;
    int mpi_size_;
    bool initialized_;
};

#endif
