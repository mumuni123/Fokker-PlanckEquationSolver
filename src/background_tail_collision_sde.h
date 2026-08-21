#ifndef BACKGROUND_TAIL_COLLISION_SDE_H
#define BACKGROUND_TAIL_COLLISION_SDE_H

#include "background_tail_collision.h"
#include "collision_coefficients.h"

// Stage-H8 Kramers-Moyal SDE backend (section 10.3.2): advances each tail
// particle with
//   du = A_u dt + B_u dW,  B_u B_u^T = D_u,
// where A_u and D_u come from the same-species coefficient provider used by
// the bulk operator.  The drift includes the relativistic Ito/Jacobian
// correction term (1/2) dD/d u.
bool sde_collide(BackgroundTailPIC& tail, const SpatialGrid& grid,
                 const TailCollisionRequest& request,
                 const CollisionCoefficientProvider& provider,
                 const std::vector<LocalCollisionMoments>& cell_moments,
                 TailCollisionDiagnostics& diagnostics);

#endif
