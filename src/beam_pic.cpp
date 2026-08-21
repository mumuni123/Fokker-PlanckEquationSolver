#include "beam_pic.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <mpi.h>
#include <omp.h>

BeamPIC::BeamPIC()
    : injection_remainder_(0.0),
      last_injected_energy_(0.0),
      cumulative_injected_energy_(0.0),
      last_outflow_energy_(0.0),
      cumulative_outflow_energy_(0.0),
      last_injected_number_(0.0),
      last_outflow_number_(0.0),
      last_injected_current_(0.0),
      last_outflow_current_(0.0),
      last_field_work_(0.0),
      step_dt_(0.0),
      step_signed_outflow_number_(0.0),
      interval_injected_number_(0.0),
      interval_left_outflow_signed_number_(0.0),
      interval_right_outflow_number_(0.0),
      interval_left_guard_path_number_(0.0),
      interval_right_guard_path_number_(0.0),
      last_continuity_abs_l1_residual_(0.0),
      last_continuity_abs_linf_residual_(0.0),
      last_continuity_l1_error_(0.0),
      last_continuity_linf_error_(0.0),
      last_boundary_flux_error_(0.0),
      last_boundary_source_error_(0.0),
      last_open_face_error_(0.0),
      last_trajectory_reconstruction_error_(0.0),
      rng_state_(0x9e3779b97f4a7c15ULL)
{}

BeamPersistentState BeamPIC::export_persistent_state() const
{
    // Checkpoint serialization and restart hashes include this POD object.
    // Zero its padding so equivalent Beam states have byte-identical images.
    BeamPersistentState s = {};
    s.injection_remainder = injection_remainder_;
    s.cumulative_injected_energy = cumulative_injected_energy_;
    s.cumulative_outflow_energy = cumulative_outflow_energy_;
    s.last_injected_energy = last_injected_energy_;
    s.last_outflow_energy = last_outflow_energy_;
    s.last_injected_number = last_injected_number_;
    s.last_outflow_number = last_outflow_number_;
    s.last_injected_current = last_injected_current_;
    s.last_outflow_current = last_outflow_current_;
    s.last_field_work = last_field_work_;
    s.step_dt = step_dt_;
    s.step_signed_outflow_number = step_signed_outflow_number_;
    s.interval_injected_number = interval_injected_number_;
    s.interval_left_outflow_signed_number = interval_left_outflow_signed_number_;
    s.interval_right_outflow_number = interval_right_outflow_number_;
    s.interval_left_guard_path_number = interval_left_guard_path_number_;
    s.interval_right_guard_path_number = interval_right_guard_path_number_;
    s.last_continuity_abs_l1_residual = last_continuity_abs_l1_residual_;
    s.last_continuity_abs_linf_residual = last_continuity_abs_linf_residual_;
    s.last_continuity_l1_error = last_continuity_l1_error_;
    s.last_continuity_linf_error = last_continuity_linf_error_;
    s.last_boundary_flux_error = last_boundary_flux_error_;
    s.last_trajectory_reconstruction_error = last_trajectory_reconstruction_error_;
    s.rng_state = rng_state_;
    return s;
}

void BeamPIC::import_persistent_state(const BeamPersistentState& s,
                                      const SpatialGrid& sg)
{
    injection_remainder_ = s.injection_remainder;
    cumulative_injected_energy_ = s.cumulative_injected_energy;
    cumulative_outflow_energy_ = s.cumulative_outflow_energy;
    last_injected_energy_ = s.last_injected_energy;
    last_outflow_energy_ = s.last_outflow_energy;
    last_injected_number_ = s.last_injected_number;
    last_outflow_number_ = s.last_outflow_number;
    last_injected_current_ = s.last_injected_current;
    last_outflow_current_ = s.last_outflow_current;
    last_field_work_ = s.last_field_work;
    step_dt_ = s.step_dt;
    step_signed_outflow_number_ = s.step_signed_outflow_number;
    interval_injected_number_ = s.interval_injected_number;
    interval_left_outflow_signed_number_ = s.interval_left_outflow_signed_number;
    interval_right_outflow_number_ = s.interval_right_outflow_number;
    interval_left_guard_path_number_ = s.interval_left_guard_path_number;
    interval_right_guard_path_number_ = s.interval_right_guard_path_number;
    last_continuity_abs_l1_residual_ = s.last_continuity_abs_l1_residual;
    last_continuity_abs_linf_residual_ = s.last_continuity_abs_linf_residual;
    last_continuity_l1_error_ = s.last_continuity_l1_error;
    last_continuity_linf_error_ = s.last_continuity_linf_error;
    last_boundary_flux_error_ = s.last_boundary_flux_error;
    last_boundary_source_error_ = s.last_boundary_flux_error;
    last_open_face_error_ = 0.0;
    last_trajectory_reconstruction_error_ = s.last_trajectory_reconstruction_error;
    rng_state_ = s.rng_state;
    midpoint_state_ = BeamMidpointState();
    send_left_.clear(); send_right_.clear(); keep_.clear();
    recv_left_.clear(); recv_right_.clear(); thread_keep_.clear();
    thread_send_left_.clear(); thread_send_right_.clear();
    trajectory_density_delta_.assign(static_cast<size_t>(sg.nx_local), 0.0);
    boundary_source_density_delta_.assign(static_cast<size_t>(sg.nx_local), 0.0);
}

namespace {
size_t initial_particle_capacity(const SpatialGrid& sg)
{
    const double active_time = std::max(
        0.0, std::min(Param::t_end, Param::t_inject_end) - Param::t_inject_start);
    const double downstream_length = sg.dx * static_cast<double>(sg.nx_global);
    const double active_length = std::min(downstream_length, Param::beam_v0 * active_time);
    const int active_cells = std::max(1, static_cast<int>(std::ceil(active_length / sg.dx)) + 4);
    const int local_capacity_cells = std::max(1, std::min(sg.nx_local + 2 * sg.nghost,
                                                          active_cells + 2 * sg.nghost));
    return static_cast<size_t>(local_capacity_cells)
         * static_cast<size_t>(Param::beam_macro_particles_per_cell);
}

int particle_byte_count(int particle_count)
{
    return static_cast<int>(static_cast<size_t>(particle_count) * sizeof(BeamParticle));
}

inline double gather_staggered_ex(double x,
                                  const SpatialGrid& sg,
                                  const EMFields& fields)
{
    const double length = sg.dx * static_cast<double>(sg.nx_global);
    if (x < 0.0 || x >= length) return 0.0;

    // Ex is x-face centered. Relative to the nearest face, the shifted CIC
    // weights are h_x(-1), h_x(0), h_x(1); only two are nonzero.
    const double face_coordinate = x * (1.0 / sg.dx);
    const int center_global =
        static_cast<int>(std::floor(face_coordinate + 0.5));
    const double shift = face_coordinate - center_global;
    const double h_center = 1.0 - std::fabs(shift);
    const int center_local = center_global - sg.ix_start;
    const size_t face_count = fields.Ex_face.size();

    double ex = 0.0;
    if (center_local >= 0 &&
        static_cast<size_t>(center_local) < face_count) {
        ex = h_center * fields.Ex_face[static_cast<size_t>(center_local)];
    }
    if (shift < 0.0) {
        const int left_local = center_local - 1;
        if (left_local >= 0 && static_cast<size_t>(left_local) < face_count) {
            ex += (-shift) * fields.Ex_face[static_cast<size_t>(left_local)];
        }
    } else if (shift > 0.0) {
        const int right_local = center_local + 1;
        if (right_local >= 0 &&
            static_cast<size_t>(right_local) < face_count) {
            ex += shift * fields.Ex_face[static_cast<size_t>(right_local)];
        }
    }
    return ex;
}

inline double beam_kinetic_energy_per_particle(double px)
{
    const double pnorm = px / (Const::me * Const::c);
    const double gamma = std::sqrt(1.0 + pnorm * pnorm);
    return (gamma - 1.0) * Const::me * Const::c * Const::c;
}

inline double beam_kinetic_energy(const BeamParticle& p)
{
    return p.weight * beam_kinetic_energy_per_particle(p.px);
}

inline double beam_vx(double px)
{
    const double u = px / (Const::me * Const::c);
    return Const::c * u / std::sqrt(1.0 + u * u);
}

// Pure drift at constant vx with an exact boundary-crossing time.  Returns
// true if the particle stays inside; on a crossing sets side (-1 left,
// +1 right) and pushed_dt (time to the boundary, within [0, dt]).
inline bool drift_step(double x0, double vx, double dt, double length,
                       int& side, double& pushed_dt)
{
    const double x1 = x0 + vx * dt;
    side = 0;
    pushed_dt = dt;
    if (x1 < 0.0 || x1 >= length) {
        side = (x1 < 0.0) ? -1 : 1;
        const double boundary = (side < 0) ? 0.0 : length;
        pushed_dt = (vx != 0.0)
            ? std::max(0.0, std::min(dt, (boundary - x0) / vx)) : 0.0;
        return false;
    }
    return true;
}

inline double random_unit(unsigned long long& state)
{
    state = state * 2862933555777941757ULL + 3037000493ULL;
    return static_cast<double>(state >> 11) *
           (1.0 / 9007199254740992.0);
}

void resize_or_zero(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    } else {
        std::fill(values.begin(), values.end(), 0.0);
    }
}

void resize_or_zero(std::vector<size_t>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0);
    } else {
        std::fill(values.begin(), values.end(), 0);
    }
}

} // namespace

void BeamPIC::init(const SpatialGrid& sg)
{
    particles.clear();
    density.assign(sg.nx_local, 0.0);
    current_x.assign(sg.nx_local, 0.0);
    current_face_x.assign(sg.nx_local + 1, 0.0);
    density_step_start.assign(sg.nx_local, 0.0);
    boundary_source_density_delta_.assign(sg.nx_local, 0.0);
    injection_remainder_ = 0.0;
    last_injected_energy_ = 0.0;
    cumulative_injected_energy_ = 0.0;
    last_outflow_energy_ = 0.0;
    cumulative_outflow_energy_ = 0.0;
    last_injected_number_ = 0.0;
    last_outflow_number_ = 0.0;
    last_injected_current_ = 0.0;
    last_outflow_current_ = 0.0;
    last_field_work_ = 0.0;
    step_dt_ = 0.0;
    step_signed_outflow_number_ = 0.0;
    interval_injected_number_ = 0.0;
    interval_left_outflow_signed_number_ = 0.0;
    interval_right_outflow_number_ = 0.0;
    interval_left_guard_path_number_ = 0.0;
    interval_right_guard_path_number_ = 0.0;
    last_continuity_abs_l1_residual_ = 0.0;
    last_continuity_abs_linf_residual_ = 0.0;
    last_continuity_l1_error_ = 0.0;
    last_continuity_linf_error_ = 0.0;
    last_boundary_flux_error_ = 0.0;
    last_boundary_source_error_ = 0.0;
    last_open_face_error_ = 0.0;
    last_trajectory_reconstruction_error_ = 0.0;
    rng_state_ = 0x9e3779b97f4a7c15ULL;
    midpoint_state_ = BeamMidpointState();

    const size_t capacity = initial_particle_capacity(sg);
    const size_t boundary_capacity = std::max(
        static_cast<size_t>(8 * Param::beam_macro_particles_per_cell),
        static_cast<size_t>((2 * sg.nghost + 4) * Param::beam_macro_particles_per_cell));
    particles.reserve(capacity);
    keep_.reserve(capacity);
    send_left_.reserve(std::min(capacity, boundary_capacity));
    send_right_.reserve(std::min(capacity, boundary_capacity));
    recv_left_.reserve(std::min(capacity, boundary_capacity));
    recv_right_.reserve(std::min(capacity, boundary_capacity));
}

void BeamPIC::begin_step(const SpatialGrid& sg, double dt)
{
    step_dt_ = dt;
    last_injected_energy_ = 0.0;
    last_outflow_energy_ = 0.0;
    last_injected_number_ = 0.0;
    last_outflow_number_ = 0.0;
    last_injected_current_ = 0.0;
    last_outflow_current_ = 0.0;
    last_field_work_ = 0.0;
    step_signed_outflow_number_ = 0.0;
    interval_injected_number_ = 0.0;
    interval_left_outflow_signed_number_ = 0.0;
    interval_right_outflow_number_ = 0.0;
    interval_left_guard_path_number_ = 0.0;
    interval_right_guard_path_number_ = 0.0;
    last_continuity_abs_l1_residual_ = 0.0;
    last_continuity_abs_linf_residual_ = 0.0;
    last_continuity_l1_error_ = 0.0;
    last_continuity_linf_error_ = 0.0;
    last_boundary_flux_error_ = 0.0;
    last_boundary_source_error_ = 0.0;
    last_open_face_error_ = 0.0;
    last_trajectory_reconstruction_error_ = 0.0;
    midpoint_state_ = BeamMidpointState();
    begin_current_interval(sg);
}

void BeamPIC::begin_current_interval(const SpatialGrid& sg)
{
    density_step_start = density;
    if (current_face_x.size() != static_cast<size_t>(sg.nx_local + 1)) {
        current_face_x.assign(sg.nx_local + 1, 0.0);
        current_x.assign(sg.nx_local, 0.0);
    }
    resize_or_zero(trajectory_density_delta_,
                   static_cast<size_t>(sg.nx_local));
    resize_or_zero(boundary_source_density_delta_,
                   static_cast<size_t>(sg.nx_local));
    if (reconstructed_current_face_x_.size() !=
        static_cast<size_t>(sg.nx_local + 1)) {
        reconstructed_current_face_x_.resize(
            static_cast<size_t>(sg.nx_local + 1));
    }
    interval_injected_number_ = 0.0;
    interval_left_outflow_signed_number_ = 0.0;
    interval_right_outflow_number_ = 0.0;
    interval_left_guard_path_number_ = 0.0;
    interval_right_guard_path_number_ = 0.0;
}

BeamInjectionSchedule BeamPIC::generate_injection_schedule(
    const SpatialGrid& sg, double step_start, double dt, int mpi_rank) const
{
    BeamInjectionSchedule schedule;
    schedule.rng_state_after_generation = rng_state_;
    schedule.remainder_after_generation = injection_remainder_;
    if (mpi_rank != 0 || dt <= 0.0) return schedule;

    const double active_start = std::max(step_start, Param::t_inject_start);
    const double active_end = std::min(step_start + dt, Param::t_inject_end);
    const double active_dt = active_end - active_start;
    if (!(active_dt > 0.0)) return schedule;

    const double depth_cells = schedule.remainder_after_generation +
        Param::beam_v0 * active_dt / sg.dx;
    const double macro_depth = depth_cells * Param::beam_macro_particles_per_cell;
    const int count = static_cast<int>(macro_depth);
    schedule.remainder_after_generation = depth_cells -
        static_cast<double>(count) / Param::beam_macro_particles_per_cell;
    unsigned long long rng = schedule.rng_state_after_generation;
    schedule.events.reserve(static_cast<size_t>(std::max(0, count)));
    for (int i = 0; i < count; ++i) {
        BeamInjectionEvent event;
        event.crossing_time = active_start + random_unit(rng) * active_dt;
        event.px = Param::beam_p0;
        event.weight = beam_macro_weight(sg);
        schedule.events.push_back(event);
    }
    schedule.rng_state_after_generation = rng;
    return schedule;
}

void BeamPIC::commit_injection_schedule(const BeamInjectionSchedule& schedule,
                                        int mpi_rank)
{
    if (mpi_rank != 0) return;
    rng_state_ = schedule.rng_state_after_generation;
    injection_remainder_ = schedule.remainder_after_generation;
}

void BeamPIC::predict_to_midpoint(const BeamInjectionSchedule& schedule,
                                  const SpatialGrid& sg, const EMFields& field_n,
                                  double time, double dt, int mpi_rank,
                                  int mpi_size, BeamPIC& midpoint) const
{
    (void)field_n;  // The first half is a pure drift; no kick before E^{n+1/2}.
    (void)mpi_size;
    BeamMidpointState& ms = midpoint.midpoint_state_;
    ms.particles.clear();
    ms.active_remaining_dt.clear();
    ms.kick_dt.clear();
    ms.left_outflow_number = 0.0;
    ms.right_outflow_number = 0.0;
    ms.left_outflow_energy = 0.0;
    ms.right_outflow_energy = 0.0;
    ms.injected_number = 0.0;
    ms.injected_energy = 0.0;
    ms.trajectory_send_left = 0.0;
    ms.trajectory_send_right = 0.0;
    resize_or_zero(ms.trajectory_delta,
                   static_cast<size_t>(sg.nx_local));
    resize_or_zero(ms.boundary_source_delta,
                   static_cast<size_t>(sg.nx_local));

    const double half = 0.5 * dt;
    const double t_mid = time + half;
    const double length = sg.dx * static_cast<double>(sg.nx_global);

    // Existing particles: pure drift x^{n+1/2} = x^n + (dt/2) v(p^n), with
    // exact boundary-crossing removal and trajectory shape-change accounting.
    for (size_t i = 0; i < particles.size(); ++i) {
        const BeamParticle& p = particles[i];
        const double vx = beam_vx(p.px);
        const double x_mid = p.x + vx * half;
        int side = 0;
        double pushed_dt = half;
            if (x_mid < 0.0 || x_mid >= length) {
                drift_step(p.x, vx, half, length, side, pushed_dt);
                const double x_cross = p.x + vx * pushed_dt;
                add_shape_density(sg, ms.trajectory_delta,
                              ms.trajectory_send_left,
                              ms.trajectory_send_right,
                              p.x, -p.weight);
            add_shape_density(sg, ms.trajectory_delta,
                              ms.trajectory_send_left,
                              ms.trajectory_send_right,
                              x_cross, p.weight);
            const double ke = beam_kinetic_energy(p);
            if (side < 0) {
                ms.left_outflow_number += p.weight;
                ms.left_outflow_energy += ke;
            } else {
                ms.right_outflow_number += p.weight;
                ms.right_outflow_energy += ke;
            }
            continue;
        }
        BeamParticle mid = p;
        mid.x = x_mid;
        ms.particles.push_back(mid);
        ms.active_remaining_dt.push_back(half);
        ms.kick_dt.push_back(dt);
        add_shape_density(sg, ms.trajectory_delta,
                          ms.trajectory_send_left,
                          ms.trajectory_send_right,
                          p.x, -p.weight);
        add_shape_density(sg, ms.trajectory_delta,
                          ms.trajectory_send_left,
                          ms.trajectory_send_right,
                          x_mid, p.weight);
    }

    // Particles injected at t_cross <= t_mid: the midpoint deposit position
    // is the real drift from x = 0 to t_mid (section 13.9).
    if (mpi_rank == 0) {
        for (size_t i = 0; i < schedule.events.size(); ++i) {
            const BeamInjectionEvent& event = schedule.events[i];
            if (event.crossing_time > t_mid || !(event.weight > 0.0)) continue;
            const double v0 = beam_vx(event.px);
            const double x_mid = v0 * (t_mid - event.crossing_time);
            int side = 0;
            double pushed_dt = t_mid - event.crossing_time;
            if (x_mid < 0.0 || x_mid >= length) {
                drift_step(0.0, v0, t_mid - event.crossing_time, length,
                           side, pushed_dt);
                const double x_cross = v0 * pushed_dt;
                add_shape_density(sg, ms.trajectory_delta,
                                  ms.trajectory_send_left,
                                  ms.trajectory_send_right,
                                  0.0, -event.weight);
                add_shape_density(sg, ms.trajectory_delta,
                                  ms.trajectory_send_left,
                                  ms.trajectory_send_right,
                                  x_cross, event.weight);
                const double ke = event.weight *
                    beam_kinetic_energy_per_particle(event.px);
                if (side < 0) {
                    ms.left_outflow_number += event.weight;
                    ms.left_outflow_energy += ke;
                } else {
                    ms.right_outflow_number += event.weight;
                    ms.right_outflow_energy += ke;
                }
            } else {
                BeamParticle mid;
                mid.x = x_mid;
                mid.px = event.px;
                mid.weight = event.weight;
                ms.particles.push_back(mid);
                ms.active_remaining_dt.push_back(half);
                ms.kick_dt.push_back((time + dt) - event.crossing_time);
                add_shape_density(sg, ms.trajectory_delta,
                                  ms.trajectory_send_left,
                                  ms.trajectory_send_right,
                                  0.0, -event.weight);
                add_shape_density(sg, ms.trajectory_delta,
                                  ms.trajectory_send_left,
                                  ms.trajectory_send_right,
                                  x_mid, event.weight);
            }
            ms.injected_number += event.weight;
            ms.injected_energy += event.weight *
                beam_kinetic_energy_per_particle(event.px);
            add_shape_density(sg, ms.boundary_source_delta,
                              ms.trajectory_send_left,
                              ms.trajectory_send_right,
                              0.0, event.weight);
        }
    }
    // Move the midpoint particle list into the working beam (zero copy) so
    // the midpoint density deposit sees it; the per-particle metadata stays
    // in the midpoint state for the final pass.
    midpoint.particles.swap(ms.particles);
}

void BeamPIC::finish_from_midpoint(const BeamInjectionSchedule& schedule,
                                   const SpatialGrid& sg,
                                   const EMFields& field_mid, double time,
                                   double dt, int mpi_rank, int mpi_size,
                                   std::vector<double>* local_delta_ke_by_x,
                                   double* local_delta_ke_boundary)
{
    if (local_delta_ke_by_x != NULL) {
        local_delta_ke_by_x->assign(static_cast<size_t>(sg.nx_local), 0.0);
    }
    BeamMidpointState& ms = midpoint_state_;
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    // Self-contained sizing: the production flow calls begin_step first, but
    // the direct predict/finish unit tests do not.
    resize_or_zero(trajectory_density_delta_, nxl);
    resize_or_zero(boundary_source_density_delta_, nxl);
    // Fold the first-half (pure-drift) accumulations into the step ledgers.
    if (trajectory_density_delta_.size() == nxl &&
        ms.trajectory_delta.size() == nxl) {
        for (size_t ix = 0; ix < nxl; ++ix) {
            trajectory_density_delta_[ix] += ms.trajectory_delta[ix];
        }
    }
    if (boundary_source_density_delta_.size() == nxl &&
        ms.boundary_source_delta.size() == nxl) {
        for (size_t ix = 0; ix < nxl; ++ix) {
            boundary_source_density_delta_[ix] += ms.boundary_source_delta[ix];
        }
    }
    interval_injected_number_ += ms.injected_number;
    last_injected_number_ += ms.injected_number;
    last_injected_energy_ += ms.injected_energy;
    cumulative_injected_energy_ += ms.injected_energy;
    interval_left_outflow_signed_number_ -= ms.left_outflow_number;
    interval_right_outflow_number_ += ms.right_outflow_number;
    last_outflow_number_ += ms.left_outflow_number + ms.right_outflow_number;
    step_signed_outflow_number_ +=
        ms.right_outflow_number - ms.left_outflow_number;
    last_outflow_energy_ += ms.left_outflow_energy + ms.right_outflow_energy;
    cumulative_outflow_energy_ +=
        ms.left_outflow_energy + ms.right_outflow_energy;
    if (mpi_rank == 0) {
        interval_left_guard_path_number_ += ms.trajectory_send_left * sg.dx;
    }
    if (mpi_rank == mpi_size - 1) {
        interval_right_guard_path_number_ += ms.trajectory_send_right * sg.dx;
    }
    double trajectory_send_left = ms.trajectory_send_left;
    double trajectory_send_right = ms.trajectory_send_right;

    const double length = sg.dx * static_cast<double>(sg.nx_global);
    const double t_mid = time + 0.5 * dt;
    const double x_left = sg.ix_start * sg.dx;
    const double x_right = (sg.ix_start + sg.nx_local) * sg.dx;
    const long long np = static_cast<long long>(particles.size());
    const int nthreads = std::max(1, omp_get_max_threads());

    if (thread_keep_.size() != static_cast<size_t>(nthreads)) {
        thread_keep_.resize(nthreads);
        thread_send_left_.resize(nthreads);
        thread_send_right_.resize(nthreads);
        thread_trajectory_density_delta_.resize(nthreads);
        thread_send_left_trajectory_.resize(nthreads);
        thread_send_right_trajectory_.resize(nthreads);
    }
    const size_t reserve_each = static_cast<size_t>(np / nthreads + 1);
    for (int t = 0; t < nthreads; ++t) {
        thread_keep_[t].clear();
        thread_send_left_[t].clear();
        thread_send_right_[t].clear();
        thread_keep_[t].reserve(reserve_each);
        thread_send_left_[t].reserve(std::max<size_t>(1, reserve_each / 16));
        thread_send_right_[t].reserve(std::max<size_t>(1, reserve_each / 16));
        resize_or_zero(thread_trajectory_density_delta_[t], nxl);
        thread_send_left_trajectory_[t] = 0.0;
        thread_send_right_trajectory_[t] = 0.0;
    }

    double field_work = 0.0;
    double boundary_work = 0.0;
    double outflow_energy = 0.0;
    double left_outflow_number = 0.0;
    double right_outflow_number = 0.0;
    double second_half_injected_number = 0.0;
    double second_half_injected_energy = 0.0;

    #pragma omp parallel reduction(+:field_work,boundary_work,outflow_energy, \
                                   left_outflow_number,right_outflow_number, \
                                   second_half_injected_number, \
                                   second_half_injected_energy)
    {
        const int tid = omp_get_thread_num();
        std::vector<BeamParticle>& local_keep = thread_keep_[tid];
        std::vector<BeamParticle>& local_left = thread_send_left_[tid];
        std::vector<BeamParticle>& local_right = thread_send_right_[tid];
        std::vector<double>& local_trajectory =
            thread_trajectory_density_delta_[tid];
        double local_trajectory_send_left = 0.0;
        double local_trajectory_send_right = 0.0;

        // Kick with E^{n+1/2} at x^{n+1/2}, then the second half-drift.
        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            BeamParticle p = particles[static_cast<size_t>(i)];
            const double x_mid = p.x;
            const double kick_dt =
                (static_cast<size_t>(i) < ms.kick_dt.size())
                    ? ms.kick_dt[static_cast<size_t>(i)] : dt;
            const double rem_dt =
                (static_cast<size_t>(i) < ms.active_remaining_dt.size())
                    ? ms.active_remaining_dt[static_cast<size_t>(i)] : 0.5 * dt;
            const double ex = gather_staggered_ex(x_mid, sg, field_mid);
            const double p_before = p.px;
            p.px += (-Const::qe) * ex * kick_dt;
            const double delta_ke =
                (beam_kinetic_energy_per_particle(p.px) -
                 beam_kinetic_energy_per_particle(p_before)) * p.weight;
            field_work += delta_ke;
            if (local_delta_ke_by_x != NULL) {
                // Gate I (section 4.4): CIC cell attribution at x^{n+1/2}.
                // Out-of-domain shares go to the boundary work audit (section
                // I3), never renormalized into the domain.
                const double s = x_mid * (1.0 / sg.dx) - 0.5;
                const int i0 = static_cast<int>(std::floor(s));
                const double frac = s - static_cast<double>(i0);
                const int c0 = i0 - sg.ix_start;
                const int c1 = i0 + 1 - sg.ix_start;
                if (c0 >= 0 && c0 < sg.nx_local) {
                    #pragma omp atomic
                    (*local_delta_ke_by_x)[static_cast<size_t>(c0)] +=
                        delta_ke * (1.0 - frac);
                } else {
                    boundary_work += delta_ke * (1.0 - frac);
                }
                if (c1 >= 0 && c1 < sg.nx_local) {
                    #pragma omp atomic
                    (*local_delta_ke_by_x)[static_cast<size_t>(c1)] +=
                        delta_ke * frac;
                } else {
                    boundary_work += delta_ke * frac;
                }
            }
            const double vx = beam_vx(p.px);
            const double x_new = x_mid + vx * rem_dt;
            int side = 0;
            double pushed_dt = rem_dt;
            if (x_new < 0.0 || x_new >= length) {
                drift_step(x_mid, vx, rem_dt, length, side, pushed_dt);
                const double x_cross = x_mid + vx * pushed_dt;
                add_shape_density(sg, local_trajectory,
                                  local_trajectory_send_left,
                                  local_trajectory_send_right,
                                  x_mid, -p.weight);
                add_shape_density(sg, local_trajectory,
                                  local_trajectory_send_left,
                                  local_trajectory_send_right,
                                  x_cross, p.weight);
                outflow_energy += beam_kinetic_energy(p);
                if (side < 0) left_outflow_number += p.weight;
                else right_outflow_number += p.weight;
                continue;
            }
            p.x = x_new;
            add_shape_density(sg, local_trajectory,
                              local_trajectory_send_left,
                              local_trajectory_send_right,
                              x_mid, -p.weight);
            add_shape_density(sg, local_trajectory,
                              local_trajectory_send_left,
                              local_trajectory_send_right,
                              x_new, p.weight);

            if (p.weight <= 0.0) continue;
            if (p.x < x_left && mpi_rank > 0) {
                local_left.push_back(p);
            } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
                local_right.push_back(p);
            } else if (p.x >= x_left && p.x < x_right) {
                local_keep.push_back(p);
            }
        }

        // Particles injected after the midpoint: created in the second half
        // and pushed for tau = t^{n+1} - t_cross (local first-order source
        // discretization, documented in section 13.9).
        #pragma omp single
        if (mpi_rank == 0) {
            for (size_t i = 0; i < schedule.events.size(); ++i) {
                const BeamInjectionEvent& event = schedule.events[i];
                if (event.crossing_time <= t_mid ||
                    !(event.weight > 0.0)) continue;
                const double tau = (time + dt) - event.crossing_time;
                BeamParticle p;
                p.x = 0.0;
                p.px = event.px;
                p.weight = event.weight;
                const double ex = gather_staggered_ex(0.0, sg, field_mid);
                const double p_before = p.px;
                p.px += (-Const::qe) * ex * tau;
                const double delta_ke =
                    (beam_kinetic_energy_per_particle(p.px) -
                     beam_kinetic_energy_per_particle(p_before)) * p.weight;
                field_work += delta_ke;
                if (local_delta_ke_by_x != NULL) {
                    // Late-injection kick at the left boundary x=0.  The
                    // out-of-domain CIC share goes to the boundary work audit.
                    const double s = -0.5;
                    const int i0 = static_cast<int>(std::floor(s));
                    const double frac = s - static_cast<double>(i0);
                    const int c0 = i0 - sg.ix_start;
                    const int c1 = i0 + 1 - sg.ix_start;
                    if (c0 >= 0 && c0 < sg.nx_local) {
                        (*local_delta_ke_by_x)[static_cast<size_t>(c0)] +=
                            delta_ke * (1.0 - frac);
                    } else {
                        boundary_work += delta_ke * (1.0 - frac);
                    }
                    if (c1 >= 0 && c1 < sg.nx_local) {
                        (*local_delta_ke_by_x)[static_cast<size_t>(c1)] +=
                            delta_ke * frac;
                    } else {
                        boundary_work += delta_ke * frac;
                    }
                }
                const double vx = beam_vx(p.px);
                const double x_new = vx * tau;
                int side = 0;
                double pushed_dt = tau;
                if (x_new < 0.0 || x_new >= length) {
                    drift_step(0.0, vx, tau, length, side, pushed_dt);
                    const double x_cross = vx * pushed_dt;
                    add_shape_density(sg, local_trajectory,
                                      local_trajectory_send_left,
                                      local_trajectory_send_right,
                                      0.0, -p.weight);
                    add_shape_density(sg, local_trajectory,
                                      local_trajectory_send_left,
                                      local_trajectory_send_right,
                                      x_cross, p.weight);
                    outflow_energy += beam_kinetic_energy(p);
                    if (side < 0) left_outflow_number += p.weight;
                    else right_outflow_number += p.weight;
                } else {
                    p.x = x_new;
                    add_shape_density(sg, local_trajectory,
                                      local_trajectory_send_left,
                                      local_trajectory_send_right,
                                      0.0, -p.weight);
                    add_shape_density(sg, local_trajectory,
                                      local_trajectory_send_left,
                                      local_trajectory_send_right,
                                      x_new, p.weight);
                    if (p.x < x_left && mpi_rank > 0) {
                        local_left.push_back(p);
                    } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
                        local_right.push_back(p);
                    } else if (p.x >= x_left && p.x < x_right) {
                        local_keep.push_back(p);
                    }
                }
                second_half_injected_number += event.weight;
                second_half_injected_energy +=
                    event.weight * beam_kinetic_energy_per_particle(event.px);
                add_shape_density(sg, boundary_source_density_delta_,
                                  local_trajectory_send_left,
                                  local_trajectory_send_right,
                                  0.0, event.weight);
            }
        }
        thread_send_left_trajectory_[tid] = local_trajectory_send_left;
        thread_send_right_trajectory_[tid] = local_trajectory_send_right;
    }

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double trajectory_sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            trajectory_sum +=
                thread_trajectory_density_delta_[t][static_cast<size_t>(ix)];
        }
        trajectory_density_delta_[static_cast<size_t>(ix)] += trajectory_sum;
    }

    for (int t = 0; t < nthreads; ++t) {
        trajectory_send_left += thread_send_left_trajectory_[t];
        trajectory_send_right += thread_send_right_trajectory_[t];
    }
    if (mpi_rank == 0) {
        interval_left_guard_path_number_ += trajectory_send_left * sg.dx;
    }
    if (mpi_rank == mpi_size - 1) {
        interval_right_guard_path_number_ += trajectory_send_right * sg.dx;
    }
    exchange_trajectory_density(sg, mpi_rank, mpi_size,
                                trajectory_send_left,
                                trajectory_send_right);

    size_t keep_total = 0;
    size_t left_total = 0;
    size_t right_total = 0;
    resize_or_zero(keep_offsets_, static_cast<size_t>(nthreads + 1));
    resize_or_zero(left_offsets_, static_cast<size_t>(nthreads + 1));
    resize_or_zero(right_offsets_, static_cast<size_t>(nthreads + 1));
    for (int t = 0; t < nthreads; ++t) {
        keep_offsets_[t + 1] = keep_offsets_[t] + thread_keep_[t].size();
        left_offsets_[t + 1] = left_offsets_[t] + thread_send_left_[t].size();
        right_offsets_[t + 1] = right_offsets_[t] + thread_send_right_[t].size();
    }
    keep_total = keep_offsets_[nthreads];
    left_total = left_offsets_[nthreads];
    right_total = right_offsets_[nthreads];

    keep_.clear();
    send_left_.clear();
    send_right_.clear();
    keep_.resize(keep_total);
    send_left_.resize(left_total);
    send_right_.resize(right_total);

    #pragma omp parallel for schedule(static)
    for (int t = 0; t < nthreads; ++t) {
        std::copy(thread_keep_[t].begin(), thread_keep_[t].end(),
                  keep_.begin() + static_cast<std::ptrdiff_t>(keep_offsets_[t]));
        std::copy(thread_send_left_[t].begin(), thread_send_left_[t].end(),
                  send_left_.begin() + static_cast<std::ptrdiff_t>(left_offsets_[t]));
        std::copy(thread_send_right_[t].begin(), thread_send_right_[t].end(),
                  send_right_.begin() + static_cast<std::ptrdiff_t>(right_offsets_[t]));
    }

    particles.swap(keep_);
    exchange_particles(sg, mpi_rank, mpi_size);

    last_injected_number_ += second_half_injected_number;
    last_injected_energy_ += second_half_injected_energy;
    cumulative_injected_energy_ += second_half_injected_energy;
    interval_injected_number_ += second_half_injected_number;
    last_injected_current_ = (step_dt_ > 0.0)
        ? -Const::qe * last_injected_number_ / step_dt_ : 0.0;

    const double total_left_outflow =
        ms.left_outflow_number + left_outflow_number;
    const double total_right_outflow =
        ms.right_outflow_number + right_outflow_number;
    last_outflow_number_ += left_outflow_number + right_outflow_number;
    step_signed_outflow_number_ += right_outflow_number - left_outflow_number;
    last_outflow_energy_ += outflow_energy;
    cumulative_outflow_energy_ += outflow_energy;
    last_outflow_current_ = (step_dt_ > 0.0)
        ? -Const::qe * step_signed_outflow_number_ / step_dt_ : 0.0;
    last_field_work_ += field_work;
    if (local_delta_ke_boundary != NULL) {
        *local_delta_ke_boundary = boundary_work;
    }
    interval_left_outflow_signed_number_ -= left_outflow_number;
    interval_right_outflow_number_ += right_outflow_number;
    if (!boundary_source_density_delta_.empty()) {
        if (mpi_rank == 0) {
            double source_guard_left = 0.0;
            double source_guard_right = 0.0;
            add_shape_density(sg, boundary_source_density_delta_,
                              source_guard_left, source_guard_right,
                              0.0, -total_left_outflow);
        }
        if (mpi_rank == mpi_size - 1) {
            double source_guard_left = 0.0;
            double source_guard_right = 0.0;
            add_shape_density(sg, boundary_source_density_delta_,
                              source_guard_left, source_guard_right,
                              sg.dx * static_cast<double>(sg.nx_global),
                              -total_right_outflow);
        }
    }
    midpoint_state_ = BeamMidpointState();
}

void BeamPIC::swap_state(BeamPIC& other)
{
    particles.swap(other.particles);
    density.swap(other.density);
    current_x.swap(other.current_x);
    current_face_x.swap(other.current_face_x);
    density_step_start.swap(other.density_step_start);
    std::swap(injection_remainder_, other.injection_remainder_);
    std::swap(last_injected_energy_, other.last_injected_energy_);
    std::swap(cumulative_injected_energy_, other.cumulative_injected_energy_);
    std::swap(last_outflow_energy_, other.last_outflow_energy_);
    std::swap(cumulative_outflow_energy_, other.cumulative_outflow_energy_);
    std::swap(last_injected_number_, other.last_injected_number_);
    std::swap(last_outflow_number_, other.last_outflow_number_);
    std::swap(last_injected_current_, other.last_injected_current_);
    std::swap(last_outflow_current_, other.last_outflow_current_);
    std::swap(last_field_work_, other.last_field_work_);
    std::swap(step_dt_, other.step_dt_);
    std::swap(step_signed_outflow_number_, other.step_signed_outflow_number_);
    std::swap(interval_injected_number_, other.interval_injected_number_);
    std::swap(interval_left_outflow_signed_number_,
              other.interval_left_outflow_signed_number_);
    std::swap(interval_right_outflow_number_, other.interval_right_outflow_number_);
    std::swap(interval_left_guard_path_number_,
              other.interval_left_guard_path_number_);
    std::swap(interval_right_guard_path_number_,
              other.interval_right_guard_path_number_);
    std::swap(last_continuity_abs_l1_residual_,
              other.last_continuity_abs_l1_residual_);
    std::swap(last_continuity_abs_linf_residual_,
              other.last_continuity_abs_linf_residual_);
    std::swap(last_continuity_l1_error_, other.last_continuity_l1_error_);
    std::swap(last_continuity_linf_error_, other.last_continuity_linf_error_);
    std::swap(last_boundary_flux_error_, other.last_boundary_flux_error_);
    std::swap(last_boundary_source_error_, other.last_boundary_source_error_);
    std::swap(last_open_face_error_, other.last_open_face_error_);
    std::swap(last_trajectory_reconstruction_error_,
              other.last_trajectory_reconstruction_error_);
    std::swap(rng_state_, other.rng_state_);
    std::swap(midpoint_state_, other.midpoint_state_);
    trajectory_density_delta_.swap(other.trajectory_density_delta_);
    boundary_source_density_delta_.swap(other.boundary_source_density_delta_);
    reconstructed_current_face_x_.swap(other.reconstructed_current_face_x_);
}

void BeamPIC::exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    (void)sg;

    if (mpi_size == 1) {
        send_left_.clear();
        send_right_.clear();
        return;
    }

    int left = mpi_rank - 1;
    int right = mpi_rank + 1;
    int send_left_count = static_cast<int>(send_left_.size());
    int send_right_count = static_cast<int>(send_right_.size());
    int recv_left_count = 0;
    int recv_right_count = 0;

    MPI_Request count_reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(&send_left_count, 1, MPI_INT, left, 301, MPI_COMM_WORLD, &count_reqs[nreq++]);
        MPI_Irecv(&recv_left_count, 1, MPI_INT, left, 302, MPI_COMM_WORLD, &count_reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right_count, 1, MPI_INT, right, 302, MPI_COMM_WORLD, &count_reqs[nreq++]);
        MPI_Irecv(&recv_right_count, 1, MPI_INT, right, 301, MPI_COMM_WORLD, &count_reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, count_reqs, MPI_STATUSES_IGNORE);

    recv_left_.resize(recv_left_count);
    recv_right_.resize(recv_right_count);

    MPI_Request data_reqs[4];
    nreq = 0;
    if (left >= 0 && send_left_count > 0) {
        MPI_Isend(send_left_.data(), particle_byte_count(send_left_count), MPI_BYTE,
                  left, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (left >= 0 && recv_left_count > 0) {
        MPI_Irecv(recv_left_.data(), particle_byte_count(recv_left_count), MPI_BYTE,
                  left, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (right < mpi_size && send_right_count > 0) {
        MPI_Isend(send_right_.data(), particle_byte_count(send_right_count), MPI_BYTE,
                  right, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (right < mpi_size && recv_right_count > 0) {
        MPI_Irecv(recv_right_.data(), particle_byte_count(recv_right_count), MPI_BYTE,
                  right, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, data_reqs, MPI_STATUSES_IGNORE);

    particles.reserve(particles.size() + recv_left_count + recv_right_count);
    particles.insert(particles.end(), recv_left_.begin(), recv_left_.end());
    particles.insert(particles.end(), recv_right_.begin(), recv_right_.end());
}

void BeamPIC::add_shape_density(const SpatialGrid& sg,
                                std::vector<double>& local,
                                double& send_left,
                                double& send_right,
                                double x, double weight) const
{
    const double inv_dx = 1.0 / sg.dx;
    const double s = x * inv_dx - 0.5;
    const int i0 = static_cast<int>(std::floor(s));
    const double frac = s - i0;
    const int targets[2] = { i0, i0 + 1 };
    const double values[2] = {
        weight * (1.0 - frac) * inv_dx,
        weight * frac * inv_dx
    };

    for (int k = 0; k < 2; ++k) {
        const int il = targets[k] - sg.ix_start;
        if (il >= 0 && il < sg.nx_local) {
            local[static_cast<size_t>(il)] += values[k];
        } else if (il < 0) {
            send_left += values[k];
        } else {
            send_right += values[k];
        }
    }
}

void BeamPIC::exchange_trajectory_density(const SpatialGrid& sg,
                                          int mpi_rank, int mpi_size,
                                          double send_left,
                                          double send_right)
{
    if (mpi_size == 1) return;

    double recv_left = 0.0;
    double recv_right = 0.0;
    MPI_Request reqs[4];
    int nreq = 0;
    if (mpi_rank > 0) {
        MPI_Isend(&send_left, 1, MPI_DOUBLE, mpi_rank - 1, 451,
                  MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, mpi_rank - 1, 452,
                  MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (mpi_rank + 1 < mpi_size) {
        MPI_Isend(&send_right, 1, MPI_DOUBLE, mpi_rank + 1, 452,
                  MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, mpi_rank + 1, 451,
                  MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (sg.nx_local > 0) {
        trajectory_density_delta_[0] += recv_left;
        trajectory_density_delta_[static_cast<size_t>(sg.nx_local - 1)] +=
            recv_right;
    }
}

void BeamPIC::deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    if (particles.empty() && mpi_size == 1) {
        resize_or_zero(density, static_cast<size_t>(sg.nx_local));
        return;
    }

    const int nthreads = std::max(1, omp_get_max_threads());
    if (thread_density_.size() != static_cast<size_t>(nthreads)) {
        thread_density_.resize(nthreads);
        thread_send_left_density_.resize(nthreads);
        thread_send_right_density_.resize(nthreads);
    }
    for (int t = 0; t < nthreads; ++t) {
        if (thread_density_[t].size() != static_cast<size_t>(sg.nx_local)) {
            thread_density_[t].assign(sg.nx_local, 0.0);
        } else {
            std::fill(thread_density_[t].begin(), thread_density_[t].end(), 0.0);
        }
        thread_send_left_density_[t] = 0.0;
        thread_send_right_density_[t] = 0.0;
    }

    double send_left_density = 0.0;
    double send_right_density = 0.0;
    const long long np = static_cast<long long>(particles.size());
    const double inv_dx = 1.0 / sg.dx;

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::vector<double>& local_density = thread_density_[tid];
        double local_send_left_density = 0.0;
        double local_send_right_density = 0.0;

        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            const BeamParticle& p = particles[static_cast<size_t>(i)];
            const double s = p.x * inv_dx - 0.5;
            const int i0 = static_cast<int>(std::floor(s));
            const double frac = s - i0;

            double contribution = p.weight * (1.0 - frac) * inv_dx;
            int target_ig = i0;
            int il = target_ig - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                local_density[il] += contribution;
            } else if (target_ig < sg.ix_start && mpi_rank > 0) {
                local_send_left_density += contribution;
            } else if (target_ig >= sg.ix_start + sg.nx_local &&
                       mpi_rank + 1 < mpi_size) {
                local_send_right_density += contribution;
            }

            contribution = p.weight * frac * inv_dx;
            target_ig = i0 + 1;
            il = target_ig - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                local_density[il] += contribution;
            } else if (target_ig < sg.ix_start && mpi_rank > 0) {
                local_send_left_density += contribution;
            } else if (target_ig >= sg.ix_start + sg.nx_local &&
                       mpi_rank + 1 < mpi_size) {
                local_send_right_density += contribution;
            }
        }

        thread_send_left_density_[tid] = local_send_left_density;
        thread_send_right_density_[tid] = local_send_right_density;
    }

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double density_sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            density_sum += thread_density_[t][ix];
        }
        density[ix] = density_sum;
    }

    for (int t = 0; t < nthreads; ++t) {
        send_left_density += thread_send_left_density_[t];
        send_right_density += thread_send_right_density_[t];
    }

    double recv_left_density = 0.0;
    double recv_right_density = 0.0;
    int left = mpi_rank - 1;
    int right = mpi_rank + 1;

    MPI_Request reqs[4];
    int nreq = 0;
    if (mpi_size == 1) {
        return;
    }
    if (left >= 0) {
        MPI_Isend(&send_left_density, 1, MPI_DOUBLE, left, 401, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_left_density, 1, MPI_DOUBLE, left, 402, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right_density, 1, MPI_DOUBLE, right, 402, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_right_density, 1, MPI_DOUBLE, right, 401, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (left >= 0 && sg.nx_local > 0) {
        density[0] += recv_left_density;
    }
    if (right < mpi_size && sg.nx_local > 0) {
        density[sg.nx_local - 1] += recv_right_density;
    }
}

void BeamPIC::finalize_charge_conserving_current(const SpatialGrid& sg,
                                                 double elapsed_dt,
                                                 int mpi_rank,
                                                 int mpi_size)
{
    if (elapsed_dt <= 0.0) return;
    if (current_face_x.size() != static_cast<size_t>(sg.nx_local + 1)) {
        current_face_x.assign(sg.nx_local + 1, 0.0);
    }

    double local_changes[2] = { 0.0, 0.0 };
    double local_corrections[2] = { 0.0, 0.0 };
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const size_t slot = static_cast<size_t>(ix);
        const double density_change =
            (density[slot] - density_step_start[slot]
             - boundary_source_density_delta_[slot]) * sg.dx;
        const double density_corrected =
            density_change - local_corrections[0];
        const double density_updated = local_changes[0] + density_corrected;
        local_corrections[0] =
            (density_updated - local_changes[0]) - density_corrected;
        local_changes[0] = density_updated;

        const double trajectory_change =
            trajectory_density_delta_[slot] * sg.dx;
        const double trajectory_corrected =
            trajectory_change - local_corrections[1];
        const double trajectory_updated =
            local_changes[1] + trajectory_corrected;
        local_corrections[1] =
            (trajectory_updated - local_changes[1]) - trajectory_corrected;
        local_changes[1] = trajectory_updated;
    }

    double preceding_changes[2] = { 0.0, 0.0 };
    if (mpi_size > 1) {
        MPI_Exscan(local_changes, preceding_changes,
                   2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            preceding_changes[0] = 0.0;
            preceding_changes[1] = 0.0;
        }
    }

    double local_boundary_source_number = 0.0;
    if (!boundary_source_density_delta_.empty()) {
        local_boundary_source_number =
            boundary_source_density_delta_.front() * sg.dx;
        if (boundary_source_density_delta_.size() > 1) {
            local_boundary_source_number +=
                boundary_source_density_delta_.back() * sg.dx;
        }
    }
    const double local_boundary_numbers[6] = {
        interval_injected_number_,
        interval_left_outflow_signed_number_,
        interval_right_outflow_number_,
        local_boundary_source_number,
        interval_left_guard_path_number_,
        interval_right_guard_path_number_
    };
    double global_boundary_numbers[6] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    MPI_Allreduce(local_boundary_numbers, global_boundary_numbers,
                  6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double inv_dt = 1.0 / elapsed_dt;
    const double dx_over_dt = sg.dx * inv_dt;
    const double inv_dx = 1.0 / sg.dx;

    double trajectory_flux =
        (-global_boundary_numbers[4] - preceding_changes[1]) * inv_dt;
    double trajectory_flux_correction = 0.0;
    current_face_x[0] = -Const::qe * trajectory_flux;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const size_t slot = static_cast<size_t>(ix);
        const double flux_change =
            -trajectory_density_delta_[slot] * dx_over_dt;
        const double corrected =
            flux_change - trajectory_flux_correction;
        const double updated = trajectory_flux + corrected;
        trajectory_flux_correction =
            (updated - trajectory_flux) - corrected;
        trajectory_flux = updated;
        current_face_x[slot + 1] = -Const::qe * trajectory_flux;
    }

    double reconstructed_flux =
        (-global_boundary_numbers[4] - preceding_changes[0]) * inv_dt;
    double reconstructed_flux_correction = 0.0;
    reconstructed_current_face_x_[0] = -Const::qe * reconstructed_flux;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const size_t slot = static_cast<size_t>(ix);
        const double flux_change =
            -(density[slot] - density_step_start[slot]
              - boundary_source_density_delta_[slot]) * dx_over_dt;
        const double corrected =
            flux_change - reconstructed_flux_correction;
        const double updated = reconstructed_flux + corrected;
        reconstructed_flux_correction =
            (updated - reconstructed_flux) - corrected;
        reconstructed_flux = updated;
        reconstructed_current_face_x_[slot + 1] =
            -Const::qe * reconstructed_flux;
    }

    if (mpi_size > 1) {
        const double send_left_faces[2] = {
            current_face_x[0], reconstructed_current_face_x_[0]
        };
        double recv_right_faces[2] = { 0.0, 0.0 };
        MPI_Request reqs[2];
        int nreq = 0;
        if (mpi_rank > 0) {
            MPI_Isend(send_left_faces, 2, MPI_DOUBLE, mpi_rank - 1, 681,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (mpi_rank + 1 < mpi_size) {
            MPI_Irecv(recv_right_faces, 2, MPI_DOUBLE, mpi_rank + 1, 681,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
        if (mpi_rank + 1 < mpi_size) {
            const size_t right_face = static_cast<size_t>(sg.nx_local);
            current_face_x[right_face] = recv_right_faces[0];
            reconstructed_current_face_x_[right_face] = recv_right_faces[1];
        }
    }
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        current_x[static_cast<size_t>(ix)] =
            0.5 * (current_face_x[static_cast<size_t>(ix)]
                 + current_face_x[static_cast<size_t>(ix + 1)]);
    }

    last_boundary_flux_error_ = 0.0;
    last_trajectory_reconstruction_error_ = 0.0;
    const double reference_current_scale = std::max(1.0, std::fabs(Param::jb));
    const double reference_number_flux_scale =
        reference_current_scale / Const::qe;
    const double current_scale = std::max(
        reference_current_scale,
        Const::qe * std::fabs(global_boundary_numbers[0]) * inv_dt);
    for (size_t iface = 0; iface < current_face_x.size(); ++iface) {
        const double scale = std::max(
            current_scale,
            std::max(std::fabs(current_face_x[iface]),
                     std::fabs(reconstructed_current_face_x_[iface])));
        last_trajectory_reconstruction_error_ = std::max(
            last_trajectory_reconstruction_error_,
            std::fabs(current_face_x[iface] -
                      reconstructed_current_face_x_[iface]) / scale);
    }
    // In 1D the endpoint-density continuity recurrence is the exact
    // charge-conserving trajectory current (up to the independently supplied
    // open-boundary face constant).  The direct shape-change accumulator can
    // lose the destination-cell split when a particle endpoint crosses an MPI
    // subdomain by more than one CIC share.  Keep its discrepancy above as a
    // trajectory diagnostic, but make the conservative recurrence the
    // authoritative current used by downstream production/audit consumers.
    current_face_x = reconstructed_current_face_x_;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        current_x[static_cast<size_t>(ix)] =
            0.5 * (current_face_x[static_cast<size_t>(ix)] +
                   current_face_x[static_cast<size_t>(ix + 1)]);
    }
    const size_t n = std::min(density.size(), density_step_start.size());
    double interval_abs_l1 = 0.0;
    double interval_abs_linf = 0.0;
    for (size_t ix = 0; ix < n; ++ix) {
        const double delta_rho =
            -Const::qe * (density[ix] - density_step_start[ix]) * inv_dt;
        const double div_j =
            (current_face_x[ix + 1] - current_face_x[ix]) * inv_dx;
        const double source_rho =
            -Const::qe * boundary_source_density_delta_[ix] * inv_dt;
        const double err = delta_rho + div_j - source_rho;
        const double abs_err = std::fabs(err);
        interval_abs_l1 += abs_err * sg.dx;
        interval_abs_linf = std::max(interval_abs_linf, abs_err);
    }

    const double interval_injection_current =
        Const::qe * global_boundary_numbers[0] * inv_dt;
    const double injection_current_scale =
        std::max(std::fabs(Param::jb), std::fabs(interval_injection_current));
    const double interval_weight =
        (step_dt_ > 0.0) ? elapsed_dt / step_dt_ : 1.0;
    last_continuity_abs_l1_residual_ += interval_abs_l1 * interval_weight;
    last_continuity_abs_linf_residual_ =
        std::max(last_continuity_abs_linf_residual_, interval_abs_linf);
    last_continuity_l1_error_ +=
        interval_abs_l1 / injection_current_scale * interval_weight;
    last_continuity_linf_error_ = std::max(
        last_continuity_linf_error_,
        interval_abs_linf * sg.dx / injection_current_scale);

    const double expected_source_number = 0.5 * (
        global_boundary_numbers[0] + global_boundary_numbers[1]
        - global_boundary_numbers[2]);
    const double source_number_scale = std::max(
        reference_number_flux_scale * elapsed_dt,
        std::max(std::fabs(expected_source_number),
                 std::fabs(global_boundary_numbers[3])));
    last_boundary_source_error_ =
        std::fabs(global_boundary_numbers[3] - expected_source_number) /
        source_number_scale;
    last_boundary_flux_error_ = last_boundary_source_error_;
    last_open_face_error_ = 0.0;

    if (mpi_rank == mpi_size - 1) {
        const double expected_right_flux = global_boundary_numbers[5] * inv_dt;
        const double flux_scale = std::max(
            reference_number_flux_scale,
            std::max(std::fabs(global_boundary_numbers[0]) * inv_dt,
                          std::max(std::fabs(trajectory_flux),
                                   std::max(std::fabs(reconstructed_flux),
                                            std::fabs(expected_right_flux)))));
        last_open_face_error_ = std::max(
            std::fabs(trajectory_flux - expected_right_flux),
            std::fabs(reconstructed_flux - expected_right_flux)) / flux_scale;
        last_boundary_flux_error_ =
            std::max(last_boundary_flux_error_, last_open_face_error_);
    }
}

void BeamPIC::snapshot_midpoint_trajectory_current(
    const SpatialGrid& sg, double elapsed_dt, int mpi_rank, int mpi_size,
    std::vector<double>& current_face) const
{
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    current_face.assign(nxl + 1, 0.0);
    if (!(elapsed_dt > 0.0) ||
        midpoint_state_.trajectory_delta.size() != nxl) return;

    double local_change = 0.0;
    double correction = 0.0;
    for (size_t i = 0; i < nxl; ++i) {
        const double term = midpoint_state_.trajectory_delta[i] * sg.dx;
        const double corrected = term - correction;
        const double updated = local_change + corrected;
        correction = (updated - local_change) - corrected;
        local_change = updated;
    }
    double preceding_change = 0.0;
    if (mpi_size > 1) {
        MPI_Exscan(&local_change, &preceding_change, 1, MPI_DOUBLE, MPI_SUM,
                   MPI_COMM_WORLD);
        if (mpi_rank == 0) preceding_change = 0.0;
    }
    double local_left_guard = mpi_rank == 0
        ? midpoint_state_.trajectory_send_left * sg.dx : 0.0;
    double global_left_guard = 0.0;
    MPI_Allreduce(&local_left_guard, &global_left_guard, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    const double inv_dt = 1.0 / elapsed_dt;
    double number_flux = (-global_left_guard - preceding_change) * inv_dt;
    current_face[0] = -Const::qe * number_flux;
    double flux_correction = 0.0;
    for (size_t i = 0; i < nxl; ++i) {
        const double delta =
            -midpoint_state_.trajectory_delta[i] * sg.dx * inv_dt;
        const double corrected = delta - flux_correction;
        const double updated = number_flux + corrected;
        flux_correction = (updated - number_flux) - corrected;
        number_flux = updated;
        current_face[i + 1] = -Const::qe * number_flux;
    }

    if (mpi_size > 1) {
        const double send_left = current_face[0];
        double recv_right = 0.0;
        MPI_Request reqs[2];
        int nreq = 0;
        if (mpi_rank > 0) {
            MPI_Isend(&send_left, 1, MPI_DOUBLE, mpi_rank - 1, 682,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (mpi_rank + 1 < mpi_size) {
            MPI_Irecv(&recv_right, 1, MPI_DOUBLE, mpi_rank + 1, 682,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
        if (mpi_rank + 1 < mpi_size) current_face[nxl] = recv_right;
    }
}

double BeamPIC::total_particle_number(const SpatialGrid& sg) const
{
    (void)sg;
    double total = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        total += particles[i].weight;
    }
    return total;
}

double BeamPIC::total_kinetic_energy() const
{
    double total = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        total += beam_kinetic_energy(particles[i]);
    }
    return total;
}
