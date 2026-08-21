#include "background_tail_collision.h"
#include "background_tail_collision_sde.h"
#include "background_tail_nanbu_perez.h"

bool advance_tail_collision(
    BackgroundTailPIC& tail, const SpatialGrid& grid,
    const TailCollisionRequest& request,
    const CollisionCoefficientProvider* provider,
    const std::vector<LocalCollisionMoments>* cell_moments,
    TailCollisionDiagnostics& diagnostics)
{
    switch (request.kernel) {
        case TailCollisionKernel::None:
            diagnostics = TailCollisionDiagnostics();
            diagnostics.success = true;
            return true;
        case TailCollisionKernel::CoulombLandauNanbuPerez:
            return nanbu_perez_collide(tail, grid, request, diagnostics);
        case TailCollisionKernel::KramersMoyalSDE:
            if (provider == NULL || cell_moments == NULL) {
                diagnostics = TailCollisionDiagnostics();
                diagnostics.success = false;
                diagnostics.failure_reason =
                    TailCollisionFailureReason::InvalidRequest;
                return false;
            }
            return sde_collide(tail, grid, request, *provider, *cell_moments,
                               diagnostics);
        case TailCollisionKernel::TraceStationaryBackground:
            // Optional EPOCH-comparison backend (section 10.3.4): marked
            // collision_approximation=trace_stationary_background with no
            // background reaction; not a production path.  Implemented as a
            // dedicated test backend, not reachable from production.
            diagnostics = TailCollisionDiagnostics();
            diagnostics.success = false;
            diagnostics.failure_reason =
                TailCollisionFailureReason::UnsupportedKernel;
            return false;
    }
    diagnostics = TailCollisionDiagnostics();
    diagnostics.success = false;
    diagnostics.failure_reason =
        TailCollisionFailureReason::UnsupportedKernel;
    return false;
}
