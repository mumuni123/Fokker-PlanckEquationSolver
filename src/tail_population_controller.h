#ifndef TAIL_POPULATION_CONTROLLER_H
#define TAIL_POPULATION_CONTROLLER_H

#include "background_tail_pic.h"
#include "grid.h"

#include <cstdint>
#include <vector>

// Tail macro-particle population control (section 7.10, stage H5).  Only
// compressing "this step's newly converted particles" is not enough to keep
// the 120 fs particle budget bounded; this controller works on the whole
// trial tail representation inside one physical cell, the same CIC stencil
// and adjacent tail phase-space bins.
//
// Rules (section 7.10):
//  * below target_particles_per_phase_bin: no operation;
//  * above max_particles_per_phase_bin: deterministic conservative
//    compression (at most `max_support` supports, shared 7-moment module);
//  * weight ratio above max_weight_ratio: equal-weight splitting that keeps
//    the total weight, position moment and velocity moments unchanged;
//  * boundary cells, MPI-in-flight particles and particles that will exit
//    next step do not participate (the caller invokes the controller after
//    the final spatial half-drift, so every particle is locally owned and
//    in-domain; boundary cells are skipped inside);
//  * each compression reports the N/Px/Jx/K/Pixx/Piperp/Xw residuals;
//  * on failure the original particles are kept; a low-order moment match is
//    never accepted.
//
// The controller is representation compression, not a collision: it must
// not change the spectrum or produce entropy production.
class TailPopulationController {
public:
    struct Config {
        bool enabled;
        int control_interval;                 // apply every N accepted steps
        int target_particles_per_phase_bin;   // below: no operation
        int max_particles_per_phase_bin;      // above: compress
        double max_weight_ratio;              // above: equal-weight split
        int max_support;                      // compression support bound
        Config()
            : enabled(false), control_interval(1),
              target_particles_per_phase_bin(64),
              max_particles_per_phase_bin(1024),
              max_weight_ratio(8.0), max_support(7)
        {}
    };

    struct Diagnostics {
        bool applied;
        std::uint64_t particles_before_local;
        std::uint64_t particles_after_local;
        int groups_considered;
        int groups_below_target;
        int groups_compressed;
        int groups_split;
        int boundary_groups_skipped;
        int compression_fallback_count;
        int split_capped_count;
        double max_residual[7];   // N Px Jx K Pixx Piperp Xw (relative)
        Diagnostics()
            : applied(false), particles_before_local(0),
              particles_after_local(0), groups_considered(0),
              groups_below_target(0), groups_compressed(0),
              groups_split(0), boundary_groups_skipped(0),
              compression_fallback_count(0), split_capped_count(0)
        {
            for (int i = 0; i < 7; ++i) max_residual[i] = 0.0;
        }
    };

    TailPopulationController();
    void configure(const Config& config);
    bool enabled() const { return config_.enabled; }
    const Config& config() const { return config_; }
    bool active_step(int step) const;

    // Section 7.10: representation control on the tail trial particles
    // only.  Never touches Poisson, Beam, bulk remap or accepted state
    // (section 14.5.1).  Deterministic per rank; no MPI collectives are
    // called inside.
    Diagnostics apply(BackgroundTailPIC& tail_trial,
                      const SpatialGrid& grid,
                      const HybridVelocityPartition& partition,
                      int step, int mpi_rank);

private:
    Config config_;
};

#endif
