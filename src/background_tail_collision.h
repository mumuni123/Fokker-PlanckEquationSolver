#ifndef BACKGROUND_TAIL_COLLISION_H
#define BACKGROUND_TAIL_COLLISION_H

#include "background_tail_pic.h"
#include "collision_coefficients.h"
#include "grid.h"

#include <cstdint>
#include <string>
#include <vector>

// Stage-H8 unified tail collision interface (section 10.3).  The kernel
// selects the backend: a pairwise relativistic Monte-Carlo for the
// Coulomb/Landau limit, a Kramers-Moyal SDE for general kernels, or the
// EPOCH-style trace-stationary-background approximation (reference only,
// never a complete-production claim).
enum class TailCollisionKernel {
    None,
    CoulombLandauNanbuPerez,
    KramersMoyalSDE,
    TraceStationaryBackground
};

enum class TailCollisionWeightMode {
    EqualStrata,
    // Legacy CLI/checkpoint name: "virtual-split". The implementation uses
    // a bounded Sentoku--Kemp unequal-weight correction and does not create
    // persistent residual macro particles.
    VirtualSplit
};

enum class TailCollisionFailureReason {
    None = 0,
    InvalidRequest = 1,
    SubstepLimit = 2,
    ParticleGrowthBudget = 3,
    DiffusionFactorization = 4,
    UnsupportedKernel = 5
};

inline const char* tail_collision_failure_reason_name(
    TailCollisionFailureReason reason)
{
    switch (reason) {
        case TailCollisionFailureReason::None: return "none";
        case TailCollisionFailureReason::InvalidRequest:
            return "invalid-request";
        case TailCollisionFailureReason::SubstepLimit:
            return "substep-limit";
        case TailCollisionFailureReason::ParticleGrowthBudget:
            return "particle-growth-budget";
        case TailCollisionFailureReason::DiffusionFactorization:
            return "diffusion-factorization";
        case TailCollisionFailureReason::UnsupportedKernel:
            return "unsupported-kernel";
    }
    return "unknown";
}

inline const char* tail_collision_weight_mode_name(
    TailCollisionWeightMode mode)
{
    switch (mode) {
        case TailCollisionWeightMode::EqualStrata: return "equal-strata";
        case TailCollisionWeightMode::VirtualSplit: return "virtual-split";
    }
    return "unknown";
}

inline const char* tail_collision_kernel_name(TailCollisionKernel kernel)
{
    switch (kernel) {
        case TailCollisionKernel::None: return "none";
        case TailCollisionKernel::CoulombLandauNanbuPerez:
            return "coulomb-nanbu-perez";
        case TailCollisionKernel::KramersMoyalSDE:
            return "kramers-moyal-sde";
        case TailCollisionKernel::TraceStationaryBackground:
            return "trace-stationary-background";
    }
    return "unknown";
}

// Pair mask of the same collision half-step (section 10.2.2/10.2.3).  Beam
// stays collisionless by default.
struct CollisionPairMask {
    bool bulk_bulk;
    bool bulk_tail;
    bool tail_bulk;
    bool tail_tail;
    bool electron_ion;
    CollisionPairMask()
        : bulk_bulk(true), bulk_tail(true), tail_bulk(true), tail_tail(true),
          electron_ion(false)
    {}
};

struct TailCollisionRequest {
    TailCollisionKernel kernel;
    double dt;                 // collision substep duration
    int accepted_step;
    int collision_half;        // 0 = first Strang half, 1 = second
    double coulomb_log;
    double cold_velocity_cutoff;   // regularises the 1/u^3 Coulomb rate
    std::uint64_t rng_seed_base;
    int mpi_rank;
    CollisionPairMask pairs;
    TailCollisionWeightMode weight_mode;
    int max_substeps;
    // Defensive upper bound for backends that can create particles. The
    // bounded Nanbu--Perez weight correction requires no growth allowance.
    double max_particle_growth;
    TailCollisionRequest()
        : kernel(TailCollisionKernel::None), dt(0.0), accepted_step(0),
          collision_half(0), coulomb_log(20.0),
          cold_velocity_cutoff(0.05), rng_seed_base(0), mpi_rank(0),
          weight_mode(TailCollisionWeightMode::EqualStrata),
          max_substeps(1024), max_particle_growth(0.0)
    {}
};

struct TailCollisionDiagnostics {
    bool success;
    TailCollisionFailureReason failure_reason;
    double number_before;
    double number_after;
    double px_before;
    double px_after;
    double py_before;
    double py_after;
    double pz_before;
    double pz_after;
    double ke_before;
    double ke_after;
    double max_s12;
    double large_angle_fraction;
    std::uint64_t weight_split_count;
    double max_du;
    int collision_substeps;
    std::uint64_t particle_count_before;
    std::uint64_t particle_count_attempted;
    std::uint64_t particle_count_limit;
    double wall_seconds;
    TailCollisionDiagnostics()
        : success(false), failure_reason(TailCollisionFailureReason::None),
          number_before(0.0), number_after(0.0),
          px_before(0.0), px_after(0.0), py_before(0.0), py_after(0.0),
          pz_before(0.0), pz_after(0.0), ke_before(0.0), ke_after(0.0),
          max_s12(0.0), large_angle_fraction(0.0), weight_split_count(0),
          max_du(0.0), collision_substeps(1), particle_count_before(0),
          particle_count_attempted(0), particle_count_limit(0),
          wall_seconds(0.0)
    {}
};

// Unified production entry: dispatches to the backend selected by the
// request kernel (section 10.3).  Operates on the tail trial only; a failed
// backend leaves the particle list unchanged.
bool advance_tail_collision(BackgroundTailPIC& tail, const SpatialGrid& grid,
                            const TailCollisionRequest& request,
                            const CollisionCoefficientProvider* provider,
                            const std::vector<LocalCollisionMoments>* cell_moments,
                            TailCollisionDiagnostics& diagnostics);

#endif
