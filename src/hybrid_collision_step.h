#ifndef HYBRID_COLLISION_STEP_H
#define HYBRID_COLLISION_STEP_H

#include "background_tail_collision.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "species.h"

// Stage-H8 unified collision half-step (sections 10.4, 10.6 and 10.7):
// `advance` receives the bulk trial and the tail trial and performs, inside
// one collision half-step, the pair-masked bulk FP update, the tail--tail
// backend and the tail--bulk SDE with the explicit field-particle reaction.
// Only a combined ledger (tail-bulk pair balance + conservation) passes.
enum class HybridCollisionFailureReason {
    None = 0,
    InvalidConfig = 1,
    BulkBulk = 2,
    TailTail = 3,
    TailBulk = 4,
    BulkReaction = 5
};

inline const char* hybrid_collision_failure_reason_name(
    HybridCollisionFailureReason reason)
{
    switch (reason) {
        case HybridCollisionFailureReason::None: return "none";
        case HybridCollisionFailureReason::InvalidConfig:
            return "invalid-config";
        case HybridCollisionFailureReason::BulkBulk: return "bulk-bulk";
        case HybridCollisionFailureReason::TailTail: return "tail-tail";
        case HybridCollisionFailureReason::TailBulk: return "tail-bulk";
        case HybridCollisionFailureReason::BulkReaction:
            return "bulk-reaction";
    }
    return "unknown";
}

struct HybridCollisionConfig {
    CollisionCoefficientProvider* bulk_provider;
    TailCollisionKernel requested_kernel;
    TailCollisionKernel tail_tail_kernel;
    TailCollisionKernel tail_bulk_kernel;
    CollisionPairMask pairs;
    TailCollisionWeightMode weight_mode;
    int max_substeps;
    double time;
    double dt;
    int accepted_step;
    int collision_half;
    double coulomb_log;
    std::uint64_t rng_seed_base;
    int mpi_rank;
    double max_particle_growth;
    // Optional per-velocity-slot mask for the bulk FP (1 = bulk-owned,
    // 0 = tail-owned conversion region treated as a zero-flux wall,
    // section 19.3).  NULL = no mask.
    const std::vector<unsigned char>* bulk_velocity_mask;
    CollisionFaceFluxes* bulk_face_fluxes;
    const HybridVelocityPartition* bulk_partition;
    BulkCollisionIntegrator bulk_integrator;
    HybridCollisionConfig()
        : bulk_provider(NULL), requested_kernel(TailCollisionKernel::None),
          tail_tail_kernel(TailCollisionKernel::None),
          tail_bulk_kernel(TailCollisionKernel::None),
          weight_mode(TailCollisionWeightMode::EqualStrata),
          max_substeps(1024), time(0.0), dt(0.0), accepted_step(0),
          collision_half(0), coulomb_log(20.0),
          rng_seed_base(0), mpi_rank(0), max_particle_growth(1.0),
          bulk_velocity_mask(NULL), bulk_face_fluxes(NULL), bulk_partition(NULL),
          bulk_integrator(BulkCollisionIntegrator::BGK_VALIDATION)
    {}
};

struct HybridCollisionDiagnostics {
    bool success;
    HybridCollisionFailureReason failure_reason;
    int failure_cell;
    double bulk_mass_change;
    double bulk_energy_change;
    double bulk_bulk_energy_change;
    double tail_px_change;
    double tail_energy_change;
    double reaction_px_balance;
    double reaction_energy_balance;
    double tail_balance_error;
    double combined_energy_change;
    double tail_tail_px_change;
    double tail_tail_energy_change;
    double bulk_reaction_px_change;
    double bulk_reaction_energy_change;
    bool bulk_bulk_applied;
    bool tail_tail_applied;
    bool tail_bulk_applied;
    bool bulk_reaction_applied;
    TailCollisionDiagnostics tail_tail_diag;
    TailCollisionDiagnostics tail_bulk_diag;
    CollisionDiagnostics bulk_diag;
    int reaction_cells;
    int reaction_fallback_cells;
    HybridCollisionDiagnostics()
        : success(false), failure_reason(HybridCollisionFailureReason::None),
          failure_cell(-1), bulk_mass_change(0.0), bulk_energy_change(0.0),
          bulk_bulk_energy_change(0.0),
          tail_px_change(0.0), tail_energy_change(0.0),
          reaction_px_balance(0.0), reaction_energy_balance(0.0),
          tail_balance_error(0.0), combined_energy_change(0.0),
          tail_tail_px_change(0.0), tail_tail_energy_change(0.0),
          bulk_reaction_px_change(0.0), bulk_reaction_energy_change(0.0),
          bulk_bulk_applied(false), tail_tail_applied(false),
          tail_bulk_applied(false), bulk_reaction_applied(false),
          reaction_cells(0), reaction_fallback_cells(0)
    {}
};

class HybridCollisionStep {
public:
    // Operates on the trials only; a failed half-step leaves both the bulk
    // trial and the tail trial unchanged.
    HybridCollisionDiagnostics advance(Species& bulk_trial,
                                       BackgroundTailPIC& tail_trial,
                                       const SpatialGrid& grid,
                                       const HybridCollisionConfig& config);
};

#endif
