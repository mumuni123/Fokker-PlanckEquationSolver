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

class BeamPIC {
public:
    std::vector<BeamParticle> particles;
    std::vector<double> density;
    std::vector<double> current_x;
    std::vector<double> current_face_x;
    std::vector<double> density_step_start;
    std::vector<double> path_density_delta;
    std::vector<double> source_density_delta;
    std::vector<double> source_current_x;
    std::vector<double>& source_current_delta;

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
    double last_continuity_l1_error() const { return last_continuity_l1_error_; }
    double last_continuity_linf_error() const { return last_continuity_linf_error_; }

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
    double left_boundary_number_flux_;
    double last_continuity_l1_error_;
    double last_continuity_linf_error_;
    unsigned long long rng_state_;
    std::vector<BeamParticle> send_left_;
    std::vector<BeamParticle> send_right_;
    std::vector<BeamParticle> keep_;
    std::vector<BeamParticle> recv_left_;
    std::vector<BeamParticle> recv_right_;
    std::vector<std::vector<BeamParticle> > thread_keep_;
    std::vector<std::vector<BeamParticle> > thread_send_left_;
    std::vector<std::vector<BeamParticle> > thread_send_right_;
    std::vector<std::vector<double> > thread_path_density_delta_;
    std::vector<std::vector<double> > thread_path_send_left_density_;
    std::vector<std::vector<double> > thread_path_send_right_density_;
    std::vector<std::vector<double> > thread_density_;
    std::vector<std::vector<double> > thread_current_;
    std::vector<double> thread_send_left_density_;
    std::vector<double> thread_send_right_density_;
    std::vector<double> thread_send_left_current_;
    std::vector<double> thread_send_right_current_;
    std::vector<double> source_send_left_density_;
    std::vector<double> source_send_right_density_;
    std::vector<double> source_recv_left_density_;
    std::vector<double> source_recv_right_density_;
    std::vector<double> source_send_left_current_;
    std::vector<double> source_send_right_current_;
    std::vector<double> source_recv_left_current_;
    std::vector<double> source_recv_right_current_;
    std::vector<double> path_send_left_density_;
    std::vector<double> path_send_right_density_;
    std::vector<double> path_recv_left_density_;
    std::vector<double> path_recv_right_density_;
    std::vector<double> all_path_numbers_;
    std::vector<size_t> keep_offsets_;
    std::vector<size_t> left_offsets_;
    std::vector<size_t> right_offsets_;

    void exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size);
    void exchange_continuity_contributions(const SpatialGrid& sg,
                                           int mpi_rank, int mpi_size);
    void reset_continuity_exchange_buffers(const SpatialGrid& sg);
    void add_path_density_delta(const SpatialGrid& sg,
                                double x0, double x1, double weight);
    void add_path_density_delta_to(const SpatialGrid& sg,
                                   std::vector<double>& local,
                                   std::vector<double>& send_left,
                                   std::vector<double>& send_right,
                                   double x0, double x1, double weight);
    void add_density_to(const SpatialGrid& sg,
                        std::vector<double>& local,
                        std::vector<double>& send_left,
                        std::vector<double>& send_right,
                        double x, double weight);
    void add_number_to_cell(const SpatialGrid& sg,
                            std::vector<double>& local,
                            std::vector<double>& send_left,
                            std::vector<double>& send_right,
                            int target_ig,
                            double value);
};

#endif
