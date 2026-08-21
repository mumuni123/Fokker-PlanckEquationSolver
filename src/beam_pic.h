#ifndef BEAM_PIC_H
#define BEAM_PIC_H

#include "grid.h"
#include "maxwell.h"
#include <cstddef>
#include <vector>

struct BeamParticle {
    double x;
    double px;
    double weight;
};

struct BeamInjectionEvent {
    double crossing_time;
    double px;
    double weight;
};

// Generated once per physical step.  The midpoint and final passes consume
// the same immutable event list; only an accepted final pass commits its
// RNG/remainder state.
struct BeamInjectionSchedule {
    std::vector<BeamInjectionEvent> events;
    unsigned long long rng_state_after_generation;
    double remainder_after_generation;
};

struct BeamPersistentState {
    double injection_remainder;
    double cumulative_injected_energy;
    double cumulative_outflow_energy;
    double last_injected_energy;
    double last_outflow_energy;
    double last_injected_number;
    double last_outflow_number;
    double last_injected_current;
    double last_outflow_current;
    double last_field_work;
    double step_dt;
    double step_signed_outflow_number;
    double interval_injected_number;
    double interval_left_outflow_signed_number;
    double interval_right_outflow_number;
    double interval_left_guard_path_number;
    double interval_right_guard_path_number;
    double last_continuity_abs_l1_residual;
    double last_continuity_abs_linf_residual;
    double last_continuity_l1_error;
    double last_continuity_linf_error;
    double last_boundary_flux_error;
    double last_trajectory_reconstruction_error;
    unsigned long long rng_state;
};

// Per-step working state of the symmetric drift-kick-drift split
// (section 13.9).  predict_to_midpoint fills it from the accepted state and
// hands the midpoint particle list to the working BeamPIC; the final pass
// consumes the per-particle durations and the first-half ledgers.
struct BeamMidpointState {
    // Midpoint particle list (x^{n+1/2}): kept in the working BeamPIC for
    // the midpoint density deposit; this field documents the state and holds
    // the list until it is moved to the working beam.
    std::vector<BeamParticle> particles;
    // Second-half drift duration per particle (dt/2 for every particle that
    // exists at the midpoint).
    std::vector<double> active_remaining_dt;
    // Field-kick duration per particle: dt for particles present since the
    // step start; (t^{n+1} - t_cross) for particles injected before the
    // midpoint (section 13.9).
    std::vector<double> kick_dt;
    // First-half (pure-drift) ledger accumulations.
    double left_outflow_number;
    double right_outflow_number;
    double left_outflow_energy;
    double right_outflow_energy;
    double injected_number;
    double injected_energy;
    std::vector<double> trajectory_delta;
    double trajectory_send_left;
    double trajectory_send_right;
    std::vector<double> boundary_source_delta;
};

class BeamPIC {
public:
    std::vector<BeamParticle> particles;
    std::vector<double> density;
    std::vector<double> current_x;
    std::vector<double> current_face_x;
    std::vector<double> density_step_start;

    BeamPIC();

    void init(const SpatialGrid& sg);
    // Reset per-step ledgers and snapshot the step-start density.  The caller
    // copies the accepted density into `density` before calling so the
    // continuity diagnostics compare against the previous accepted state.
    void begin_step(const SpatialGrid& sg, double dt);
    void begin_current_interval(const SpatialGrid& sg);
    BeamInjectionSchedule generate_injection_schedule(const SpatialGrid& sg,
                                                       double step_start,
                                                       double dt,
                                                       int mpi_rank) const;
    void commit_injection_schedule(const BeamInjectionSchedule& schedule,
                                   int mpi_rank);

    // Symmetric leapfrog split (sections 7.2/13.9).  The first half is a pure
    // drift x^{n+1/2} = x^n + (dt/2) v(p^n) with exact boundary-crossing
    // removal; particles injected at t_cross <= t_mid enter the midpoint
    // density at their real drift position v(p0)*(t_mid - t_cross).  The
    // accepted beam is const; the working midpoint state is written into
    // `midpoint` (no full particle-array copies).
    void predict_to_midpoint(const BeamInjectionSchedule& schedule,
                             const SpatialGrid& sg, const EMFields& field_n,
                             double time, double dt, int mpi_rank,
                             int mpi_size, BeamPIC& midpoint) const;
    // Consumes the stored midpoint state: kicks every midpoint particle with
    // E^{n+1/2} at x^{n+1/2} (kick duration dt, or tau for injected
    // particles), drifts the second half, creates particles injected after
    // the midpoint (pushed for tau = t^{n+1} - t_cross), removes open-boundary
    // outflow and migrates across MPI ranks.  After this call the working
    // beam holds the accepted candidate state and the step ledgers.
    // Gate I (section 4.4): when local_delta_ke_by_x is non-NULL, the kick
    // kinetic-energy change is additionally distributed to the CIC cells at
    // x^{n+1/2} (nx_local entries, caller-owned); last_field_work() semantics
    // are unchanged and the beam is never re-pushed.  Out-of-domain CIC shares
    // are never renormalized into the domain; when local_delta_ke_boundary is
    // non-NULL they are accumulated there so sum(cell) + boundary equals
    // last_field_work() (section I3).
    void finish_from_midpoint(const BeamInjectionSchedule& schedule,
                              const SpatialGrid& sg, const EMFields& field_mid,
                              double time, double dt, int mpi_rank,
                              int mpi_size,
                              std::vector<double>* local_delta_ke_by_x = NULL,
                              double* local_delta_ke_boundary = NULL);
    void deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size);
    void finalize_charge_conserving_current(const SpatialGrid& sg,
                                            double elapsed_dt,
                                            int mpi_rank, int mpi_size);
    // Gate I read-only first-drift snapshot.  It reconstructs the same
    // charge-conserving face current from the stored midpoint trajectory
    // density difference; no particle, RNG or accepted ledger is modified.
    void snapshot_midpoint_trajectory_current(
        const SpatialGrid& sg, double elapsed_dt, int mpi_rank, int mpi_size,
        std::vector<double>& current_face) const;
    // Transactional accept: exchange the full beam state without copying the
    // particle arrays.
    void swap_state(BeamPIC& other);

    double total_particle_number(const SpatialGrid& sg) const;
    double total_kinetic_energy() const;
    double last_injected_energy() const { return last_injected_energy_; }
    double cumulative_injected_energy() const { return cumulative_injected_energy_; }
    double last_outflow_energy() const { return last_outflow_energy_; }
    double cumulative_outflow_energy() const { return cumulative_outflow_energy_; }
    double last_injected_number() const { return last_injected_number_; }
    double last_outflow_number() const { return last_outflow_number_; }
    double last_injected_current() const { return last_injected_current_; }
    double last_outflow_current() const { return last_outflow_current_; }
    double last_field_work() const { return last_field_work_; }
    double last_continuity_abs_l1_residual() const {
        return last_continuity_abs_l1_residual_;
    }
    double last_continuity_abs_linf_residual() const {
        return last_continuity_abs_linf_residual_;
    }
    double last_continuity_l1_error() const { return last_continuity_l1_error_; }
    double last_continuity_linf_error() const { return last_continuity_linf_error_; }
    double last_boundary_flux_error() const { return last_boundary_flux_error_; }
    double last_boundary_source_error() const {
        return last_boundary_source_error_;
    }
    double last_open_face_error() const { return last_open_face_error_; }
    double last_trajectory_reconstruction_error() const {
        return last_trajectory_reconstruction_error_;
    }
    // Gate I read-only source ledger.  Entries are number-density changes
    // [m^-3] created by injection/removal bookkeeping and excluded from the
    // trajectory-current divergence reconstruction.
    const std::vector<double>& boundary_source_density_delta() const {
        return boundary_source_density_delta_;
    }
    BeamPersistentState export_persistent_state() const;
    void import_persistent_state(const BeamPersistentState& state,
                                 const SpatialGrid& sg);

private:
    double injection_remainder_;
    double last_injected_energy_;
    double cumulative_injected_energy_;
    double last_outflow_energy_;
    double cumulative_outflow_energy_;
    double last_injected_number_;
    double last_outflow_number_;
    double last_injected_current_;
    double last_outflow_current_;
    double last_field_work_;
    double step_dt_;
    double step_signed_outflow_number_;
    double interval_injected_number_;
    double interval_left_outflow_signed_number_;
    double interval_right_outflow_number_;
    double interval_left_guard_path_number_;
    double interval_right_guard_path_number_;
    double last_continuity_abs_l1_residual_;
    double last_continuity_abs_linf_residual_;
    double last_continuity_l1_error_;
    double last_continuity_linf_error_;
    double last_boundary_flux_error_;
    double last_boundary_source_error_;
    double last_open_face_error_;
    double last_trajectory_reconstruction_error_;
    unsigned long long rng_state_;
    BeamMidpointState midpoint_state_;
    std::vector<BeamParticle> send_left_;
    std::vector<BeamParticle> send_right_;
    std::vector<BeamParticle> keep_;
    std::vector<BeamParticle> recv_left_;
    std::vector<BeamParticle> recv_right_;
    std::vector<std::vector<BeamParticle> > thread_keep_;
    std::vector<std::vector<BeamParticle> > thread_send_left_;
    std::vector<std::vector<BeamParticle> > thread_send_right_;
    std::vector<std::vector<double> > thread_density_;
    std::vector<std::vector<double> > thread_trajectory_density_delta_;
    std::vector<double> thread_send_left_density_;
    std::vector<double> thread_send_right_density_;
    std::vector<double> thread_send_left_trajectory_;
    std::vector<double> thread_send_right_trajectory_;
    std::vector<double> trajectory_density_delta_;
    std::vector<double> boundary_source_density_delta_;
    std::vector<double> reconstructed_current_face_x_;
    std::vector<size_t> keep_offsets_;
    std::vector<size_t> left_offsets_;
    std::vector<size_t> right_offsets_;

    void exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size);
    void add_shape_density(const SpatialGrid& sg,
                           std::vector<double>& local,
                           double& send_left,
                           double& send_right,
                           double x, double weight) const;
    void exchange_trajectory_density(const SpatialGrid& sg,
                                     int mpi_rank, int mpi_size,
                                     double send_left, double send_right);
};

#endif
