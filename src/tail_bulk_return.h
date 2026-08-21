#ifndef TAIL_BULK_RETURN_H
#define TAIL_BULK_RETURN_H

#include "background_tail_pic.h"
#include "grid.h"
#include "species.h"

#include <cstdint>
#include <vector>

// Written by checkpoint audits so an output file proves which mathematical
// projection contract was linked into the executable.
const char* tail_bulk_return_projection_schema();

// Current and pressure moments of continuous PIC velocities cannot generally
// be represented to a fixed relative tolerance by nonnegative masses at the
// Eulerian velocity-cell centres.  This helper turns the error of the native
// nearest-cell representation into a bounded, grid-aware acceptance limit.
// It is public only so the policy can be regression-tested independently of a
// full MPI return transaction.
double tail_bulk_return_representation_tolerance(
    double nearest_cell_relative_residual);

// H10 tail-to-bulk representation conversion.  This is deliberately a
// separate transactional operator: it owns no accepted state and is only
// handed the end-of-step trial bulk/tail representations by VpfpIntegrator.
struct TailBulkReturnConfig {
    bool enabled;
    double return_energy_mev;
    std::uint32_t residence_steps;
    int max_stencil_radius;
    double moment_tolerance;
    TailBulkReturnConfig()
        : enabled(false), return_energy_mev(0.0), residence_steps(0),
          max_stencil_radius(3), moment_tolerance(1.0e-12)
    {}
};

struct TailBulkReturnDiagnostics {
    std::uint64_t candidate_particles;
    std::uint64_t resident_particles;
    std::uint64_t attempted_groups;
    std::uint64_t committed_groups;
    std::uint64_t deferred_infeasible_groups;
    std::uint64_t deferred_rank_boundary_groups;
    std::uint64_t projection_invalid_input_cells;
    std::uint64_t projection_insufficient_support_cells;
    std::uint64_t projection_infeasible_invariant_cells;
    std::uint64_t projection_representation_incompatible_cells;
    std::uint64_t particles_removed;
    double number;
    double px;
    double jx_dx;
    double energy;
    double pixx_dx;
    double piperp_dx;
    double number_residual;
    double px_residual;
    double jx_residual;
    double energy_residual;
    double pixx_residual;
    double piperp_residual;
    // Signed Eulerian-added minus PIC-removed representation errors.  N, Px
    // and K are transaction invariants; Jx and pressure moments are retained
    // as diagnostics so a persistent projection bias can be accumulated
    // offline without changing the accepted state.
    double number_difference;
    double px_difference;
    double jx_difference;
    double energy_difference;
    double pixx_difference;
    double piperp_difference;
    // Gate I local representation-transfer source [m^-2 per local x cell].
    // The same value is added to bulk and removed from Tail, so the combined
    // density is unchanged by a committed H10 representation transaction.
    std::vector<double> returned_number_by_x;
    double mpi_request_residual;
    double wall_seconds;
    bool finite;
    bool committed;
    TailBulkReturnDiagnostics()
        : candidate_particles(0), resident_particles(0), attempted_groups(0),
          committed_groups(0), deferred_infeasible_groups(0),
          deferred_rank_boundary_groups(0),
          projection_invalid_input_cells(0),
          projection_insufficient_support_cells(0),
          projection_infeasible_invariant_cells(0),
          projection_representation_incompatible_cells(0),
          particles_removed(0), number(0.0), px(0.0), jx_dx(0.0),
          energy(0.0), pixx_dx(0.0), piperp_dx(0.0),
          number_residual(0.0), px_residual(0.0), jx_residual(0.0),
          energy_residual(0.0), pixx_residual(0.0), piperp_residual(0.0),
          number_difference(0.0), px_difference(0.0), jx_difference(0.0),
          energy_difference(0.0), pixx_difference(0.0),
          piperp_difference(0.0),
          mpi_request_residual(0.0), wall_seconds(0.0), finite(true),
          committed(false)
    {}
};

class TailBulkReturn {
public:
    explicit TailBulkReturn(const TailBulkReturnConfig& config =
                                TailBulkReturnConfig())
        : config_(config)
    {}

    void set_config(const TailBulkReturnConfig& config) { config_ = config; }
    const TailBulkReturnConfig& config() const { return config_; }

    // Returns false only for a transaction integrity failure.  A physically
    // infeasible local projection is represented by a deferred group and is
    // not a step failure.
    bool apply(Species& bulk_trial, BackgroundTailPIC& tail_trial,
               const SpatialGrid& grid,
               const HybridVelocityPartition& partition,
               long long accepted_step, int mpi_rank, int mpi_size,
               TailBulkReturnDiagnostics& diagnostics) const;

private:
    TailBulkReturnConfig config_;
};

#endif
