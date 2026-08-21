#ifndef BACKGROUND_TAIL_NANBU_PEREZ_H
#define BACKGROUND_TAIL_NANBU_PEREZ_H

#include "background_tail_collision.h"

// Stage-H8 Coulomb/Landau tail--tail backend (section 10.3.1): pairwise
// relativistic Monte-Carlo scattering.  Equal-weight pairs use the strict
// two-body COM update; unequal weights use virtual_weight_split (section
// 10.3.3) with the residual weight materialised as a new macro particle.
bool nanbu_perez_collide(BackgroundTailPIC& tail, const SpatialGrid& grid,
                         const TailCollisionRequest& request,
                         TailCollisionDiagnostics& diagnostics);

#endif
