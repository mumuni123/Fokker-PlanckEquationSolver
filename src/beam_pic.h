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

class BeamPIC {
public:
    std::vector<BeamParticle> particles;
    std::vector<double> density;
    std::vector<double> current_x;
    std::vector<double> current_face_x;
    std::vector<double> density_step_start;

    BeamPIC();

    void init(const SpatialGrid& sg);
    void begin_step(const SpatialGrid& sg, double dt);
    void begin_current_interval(const SpatialGrid& sg);
    void inject(const SpatialGrid& sg, const EMFields& fields,
                double dt, double time, int mpi_rank, int mpi_size);
    void push(const SpatialGrid& sg, const EMFields& fields, double dt,
              int mpi_rank, int mpi_size);
    void deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size);
    void finalize_charge_conserving_current(const SpatialGrid& sg,
                                            double elapsed_dt,
                                            int mpi_rank, int mpi_size);

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
    // Diagnostic components are transient and intentionally excluded from
    // BeamPersistentState; they are recomputed by every accepted Beam step.
    double last_boundary_source_error_;
    double last_open_face_error_;
    double last_trajectory_reconstruction_error_;
    unsigned long long rng_state_;
    size_t injected_begin_;
    std::vector<double> injected_push_dt_;
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
