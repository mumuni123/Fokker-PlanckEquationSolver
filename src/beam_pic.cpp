#include "beam_pic.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <mpi.h>
#include <omp.h>

BeamPIC::BeamPIC()
    : source_current_delta(source_current_x),
      injection_remainder_(0.0),
      last_injected_energy_(0.0),
      cumulative_injected_energy_(0.0),
      last_outflow_energy_(0.0),
      cumulative_outflow_energy_(0.0),
      last_injected_number_(0.0),
      last_outflow_number_(0.0),
      last_injected_current_(0.0),
      last_outflow_current_(0.0),
      last_field_work_(0.0),
      left_boundary_number_flux_(0.0),
      last_continuity_l1_error_(0.0),
      last_continuity_linf_error_(0.0)
{}

namespace {
size_t initial_particle_capacity(const SpatialGrid& sg)
{
    const double active_time = std::max(
        0.0, std::min(Param::t_end, Param::t_inject_end) - Param::t_inject_start);
    const double downstream_length = std::max(0.0, Param::Lx - Param::beam_source_x_start);
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

inline double gather_cic_ex(double x, const SpatialGrid& sg, const EMFields& fields)
{
    if (x < 0.0 || x >= Param::Lx) return 0.0;

    double s = x / sg.dx - 0.5;
    int i0 = static_cast<int>(std::floor(s));
    double frac = s - i0;
    double w0 = 1.0 - frac;
    double w1 = frac;
    double ex = 0.0;

    int target_ig = i0;
    int il = target_ig - sg.ix_start;
    int ix = il + sg.nghost;
    if (ix >= 0 && ix < sg.nx_total) {
        ex += w0 * fields.Ex[ix];
    }

    target_ig = i0 + 1;
    il = target_ig - sg.ix_start;
    ix = il + sg.nghost;
    if (ix >= 0 && ix < sg.nx_total) {
        ex += w1 * fields.Ex[ix];
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

inline double beam_velocity_from_px(double px)
{
    const double pnorm = px / (Const::me * Const::c);
    const double gamma = std::sqrt(1.0 + pnorm * pnorm);
    return px / (gamma * Const::me);
}

void resize_or_zero(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    } else {
        std::fill(values.begin(), values.end(), 0.0);
    }
}

}

void BeamPIC::init(const SpatialGrid& sg)
{
    particles.clear();
    density.assign(sg.nx_local, 0.0);
    current_x.assign(sg.nx_local, 0.0);
    current_face_x.assign(sg.nx_local + 1, 0.0);
    density_step_start.assign(sg.nx_local, 0.0);
    path_density_delta.assign(sg.nx_local, 0.0);
    source_density_delta.assign(sg.nx_local, 0.0);
    source_current_x.assign(sg.nx_local, 0.0);
    injection_remainder_ = 0.0;
    last_injected_energy_ = 0.0;
    cumulative_injected_energy_ = 0.0;
    last_outflow_energy_ = 0.0;
    cumulative_outflow_energy_ = 0.0;
    last_outflow_number_ = 0.0;
    last_outflow_current_ = 0.0;
    last_field_work_ = 0.0;
    left_boundary_number_flux_ = 0.0;
    last_continuity_l1_error_ = 0.0;
    last_continuity_linf_error_ = 0.0;

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
    (void)dt;
    last_injected_energy_ = 0.0;
    last_outflow_energy_ = 0.0;
    last_injected_number_ = 0.0;
    last_outflow_number_ = 0.0;
    last_injected_current_ = 0.0;
    last_outflow_current_ = 0.0;
    last_field_work_ = 0.0;
    left_boundary_number_flux_ = 0.0;
    last_continuity_l1_error_ = 0.0;
    last_continuity_linf_error_ = 0.0;
    density_step_start = density;
    if (source_density_delta.size() != static_cast<size_t>(sg.nx_local)) {
        density_step_start.assign(sg.nx_local, 0.0);
        path_density_delta.assign(sg.nx_local, 0.0);
        source_density_delta.assign(sg.nx_local, 0.0);
        source_current_x.assign(sg.nx_local, 0.0);
        current_face_x.assign(sg.nx_local + 1, 0.0);
        current_x.assign(sg.nx_local, 0.0);
    } else {
        std::fill(path_density_delta.begin(), path_density_delta.end(), 0.0);
        std::fill(source_density_delta.begin(), source_density_delta.end(), 0.0);
        std::fill(source_current_x.begin(), source_current_x.end(), 0.0);
        std::fill(current_face_x.begin(), current_face_x.end(), 0.0);
    }
}

void BeamPIC::inject(const SpatialGrid& sg, const EMFields& fields,
                     double dt, double time, int mpi_rank, int mpi_size)
{
    const double step_start = time - dt;
    const double active_start = std::max(step_start, Param::t_inject_start);
    const double active_end = std::min(time, Param::t_inject_end);
    const double active_dt = active_end - active_start;
    if (active_dt <= 0.0) return;
    reset_continuity_exchange_buffers(sg);
    send_left_.clear();
    send_right_.clear();

    const double x_left = sg.ix_start * sg.dx;
    const double x_right = (sg.ix_start + sg.nx_local) * sg.dx;
    const double x_src = Param::enable_beam_boundary_injection
        ? 0.0
        : std::max(0.0, std::min(Param::Lx, Param::beam_source_x_start));
    const bool owns_source = Param::enable_beam_boundary_injection
        ? (mpi_rank == 0)
        : ((x_src >= x_left && x_src < x_right)
           || (mpi_rank == mpi_size - 1 && x_src == Param::Lx));

    int n_new = 0;
    if (owns_source) {
        const double physical_per_area =
            Param::densb * Param::beam_v0 * active_dt + injection_remainder_;
        n_new = static_cast<int>(physical_per_area / Param::beam_macro_weight);
        injection_remainder_ = physical_per_area
            - n_new * Param::beam_macro_weight;
    }

    particles.reserve(particles.size() + static_cast<size_t>(std::max(0, n_new)));
    const double source_ke_per_particle =
        beam_kinetic_energy_per_particle(Param::beam_p0);
    double injected_energy_call = 0.0;
    double outflow_energy_call = 0.0;
    double field_work_call = 0.0;
    for (int i = 0; i < n_new; ++i) {
        const double crossing_time =
            active_start + (i + 0.5) * active_dt / static_cast<double>(n_new);
        const double x_birth = x_src;
        const double tau = std::max(0.0, time - crossing_time);
        const double ex = gather_cic_ex(x_birth, sg, fields);
        const double px0 = Param::beam_p0;
        const double px1 = px0 + (-Const::qe) * ex * tau;
        const double vx_mid = beam_velocity_from_px(0.5 * (px0 + px1));
        const double x_end_raw = x_birth + vx_mid * tau;
        const bool leaves_domain = (x_end_raw < 0.0 || x_end_raw >= Param::Lx);

        if (Param::enable_beam_boundary_injection) {
            left_boundary_number_flux_ += Param::beam_macro_weight;
        } else {
            add_source_density_and_current(sg, x_birth, Param::beam_macro_weight,
                                           vx_mid);
        }
        if (leaves_domain) {
            if (!Param::enable_beam_boundary_injection) {
                add_density_to(sg, path_density_delta,
                               path_send_left_density_,
                               path_send_right_density_,
                               x_birth, -Param::beam_macro_weight);
            }
        } else if (Param::enable_beam_boundary_injection) {
            add_density_to(sg, path_density_delta,
                           path_send_left_density_,
                           path_send_right_density_,
                           x_end_raw, Param::beam_macro_weight);
        } else {
            add_path_density_delta(sg, x_birth, x_end_raw,
                                   Param::beam_macro_weight);
        }

        BeamParticle p;
        p.x = x_end_raw;
        p.px = px1;
        p.weight = Param::beam_macro_weight;

        last_injected_number_ += p.weight;
        last_injected_current_ += -Const::qe * vx_mid * p.weight;
        injected_energy_call += p.weight * source_ke_per_particle;
        field_work_call += p.weight
            * (beam_kinetic_energy_per_particle(px1)
               - beam_kinetic_energy_per_particle(px0));

        if (leaves_domain) {
            outflow_energy_call += beam_kinetic_energy(p);
            last_outflow_number_ += p.weight;
            last_outflow_current_ += -Const::qe * vx_mid * p.weight;
        } else if (p.weight <= 0.0) {
            continue;
        } else if (p.x < x_left && mpi_rank > 0) {
            send_left_.push_back(p);
        } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
            send_right_.push_back(p);
        } else if (p.x >= x_left && p.x < x_right) {
            particles.push_back(p);
        }
    }
    last_injected_energy_ += injected_energy_call;
    cumulative_injected_energy_ += injected_energy_call;
    last_outflow_energy_ += outflow_energy_call;
    cumulative_outflow_energy_ += outflow_energy_call;
    last_field_work_ += field_work_call;

    exchange_continuity_contributions(sg, mpi_rank, mpi_size);
    exchange_particles(sg, mpi_rank, mpi_size);
}

void BeamPIC::push(const SpatialGrid& sg, const EMFields& fields, double dt,
                   int mpi_rank, int mpi_size)
{
    reset_continuity_exchange_buffers(sg);
    const double x_left = sg.ix_start * sg.dx;
    const double x_right = (sg.ix_start + sg.nx_local) * sg.dx;
    const long long np = static_cast<long long>(particles.size());
    const int nthreads = std::max(1, omp_get_max_threads());

    if (thread_keep_.size() != static_cast<size_t>(nthreads) ||
        thread_path_density_delta_.size() != static_cast<size_t>(nthreads)) {
        thread_keep_.resize(nthreads);
        thread_send_left_.resize(nthreads);
        thread_send_right_.resize(nthreads);
        thread_path_density_delta_.resize(nthreads);
        thread_path_send_left_density_.resize(nthreads);
        thread_path_send_right_density_.resize(nthreads);
    }
    const size_t reserve_each = static_cast<size_t>(np / nthreads + 1);
    for (int t = 0; t < nthreads; ++t) {
        thread_keep_[t].clear();
        thread_send_left_[t].clear();
        thread_send_right_[t].clear();
        thread_keep_[t].reserve(reserve_each);
        thread_send_left_[t].reserve(std::max<size_t>(1, reserve_each / 16));
        thread_send_right_[t].reserve(std::max<size_t>(1, reserve_each / 16));
        resize_or_zero(thread_path_density_delta_[t],
                       static_cast<size_t>(sg.nx_local));
        resize_or_zero(thread_path_send_left_density_[t],
                       static_cast<size_t>(sg.nghost));
        resize_or_zero(thread_path_send_right_density_[t],
                       static_cast<size_t>(sg.nghost));
    }

    double field_work = 0.0;
    double outflow_energy = 0.0;
    double outflow_number = 0.0;
    double outflow_current = 0.0;

    #pragma omp parallel reduction(+:field_work,outflow_energy,outflow_number,outflow_current)
    {
        const int tid = omp_get_thread_num();
        std::vector<BeamParticle>& local_keep = thread_keep_[tid];
        std::vector<BeamParticle>& local_left = thread_send_left_[tid];
        std::vector<BeamParticle>& local_right = thread_send_right_[tid];
        std::vector<double>& local_path = thread_path_density_delta_[tid];
        std::vector<double>& local_path_left = thread_path_send_left_density_[tid];
        std::vector<double>& local_path_right = thread_path_send_right_density_[tid];

        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            BeamParticle p = particles[static_cast<size_t>(i)];
            const double x_old = p.x;
            const double ex = gather_cic_ex(p.x, sg, fields);
            const double ke_before = beam_kinetic_energy(p);

            p.px += (-Const::qe) * ex * dt;
            const double vx = beam_velocity_from_px(p.px);
            p.x += vx * dt;
            const bool leaves_domain = (p.x < 0.0 || p.x >= Param::Lx);
            const double ke_after = beam_kinetic_energy(p);
            field_work += ke_after - ke_before;

            if (leaves_domain) {
                add_density_to(sg, local_path, local_path_left,
                               local_path_right, x_old, -p.weight);
                outflow_energy += ke_after;
                outflow_number += p.weight;
                outflow_current += -Const::qe * vx * p.weight;
                continue;
            } else {
                add_path_density_delta_to(sg, local_path, local_path_left,
                                          local_path_right, x_old, p.x,
                                          p.weight);
            }

            if (p.weight <= 0.0) continue;

            if (p.x < x_left && mpi_rank > 0) {
                local_left.push_back(p);
            } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
                local_right.push_back(p);
            } else if (p.x >= x_left && p.x < x_right) {
                local_keep.push_back(p);
            }
        }
    }

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double path_sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            path_sum += thread_path_density_delta_[t][static_cast<size_t>(ix)];
        }
        path_density_delta[static_cast<size_t>(ix)] += path_sum;
    }

    for (int g = 0; g < sg.nghost; ++g) {
        double path_left_sum = 0.0;
        double path_right_sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            const size_t slot = static_cast<size_t>(g);
            path_left_sum += thread_path_send_left_density_[t][slot];
            path_right_sum += thread_path_send_right_density_[t][slot];
        }
        const size_t slot = static_cast<size_t>(g);
        path_send_left_density_[slot] += path_left_sum;
        path_send_right_density_[slot] += path_right_sum;
    }

    size_t keep_total = 0;
    size_t left_total = 0;
    size_t right_total = 0;
    std::vector<size_t> keep_offsets(nthreads + 1, 0);
    std::vector<size_t> left_offsets(nthreads + 1, 0);
    std::vector<size_t> right_offsets(nthreads + 1, 0);
    for (int t = 0; t < nthreads; ++t) {
        keep_offsets[t + 1] = keep_offsets[t] + thread_keep_[t].size();
        left_offsets[t + 1] = left_offsets[t] + thread_send_left_[t].size();
        right_offsets[t + 1] = right_offsets[t] + thread_send_right_[t].size();
    }
    keep_total = keep_offsets[nthreads];
    left_total = left_offsets[nthreads];
    right_total = right_offsets[nthreads];

    keep_.clear();
    send_left_.clear();
    send_right_.clear();
    keep_.resize(keep_total);
    send_left_.resize(left_total);
    send_right_.resize(right_total);

    #pragma omp parallel for schedule(static)
    for (int t = 0; t < nthreads; ++t) {
        std::copy(thread_keep_[t].begin(), thread_keep_[t].end(),
                  keep_.begin() + static_cast<std::ptrdiff_t>(keep_offsets[t]));
        std::copy(thread_send_left_[t].begin(), thread_send_left_[t].end(),
                  send_left_.begin() + static_cast<std::ptrdiff_t>(left_offsets[t]));
        std::copy(thread_send_right_[t].begin(), thread_send_right_[t].end(),
                  send_right_.begin() + static_cast<std::ptrdiff_t>(right_offsets[t]));
    }

    particles.swap(keep_);

    exchange_particles(sg, mpi_rank, mpi_size);
    exchange_continuity_contributions(sg, mpi_rank, mpi_size);
    last_outflow_energy_ += outflow_energy;
    cumulative_outflow_energy_ += outflow_energy;
    last_outflow_number_ += outflow_number;
    last_outflow_current_ += outflow_current;
    last_field_work_ += field_work;
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

void BeamPIC::add_source_to_cell(const SpatialGrid& sg, int target_ig,
                                 double density_delta, double current_delta)
{
    add_number_to_cell(sg, source_density_delta, source_send_left_density_,
                       source_send_right_density_, target_ig, density_delta);
    add_number_to_cell(sg, source_current_x, source_send_left_current_,
                       source_send_right_current_, target_ig, current_delta);
}

void BeamPIC::add_number_to_cell(const SpatialGrid& sg,
                                 std::vector<double>& local,
                                 std::vector<double>& send_left,
                                 std::vector<double>& send_right,
                                 int target_ig,
                                 double value)
{
    if (value == 0.0) return;

    const int ng = sg.nghost;
    const int il = target_ig - sg.ix_start;
    if (il >= 0 && il < sg.nx_local) {
        local[static_cast<size_t>(il)] += value;
    } else if (il < 0 && il >= -ng) {
        send_left[static_cast<size_t>(il + ng)] += value;
    } else if (il >= sg.nx_local && il < sg.nx_local + ng) {
        send_right[static_cast<size_t>(il - sg.nx_local)] += value;
    }
}

void BeamPIC::add_source_density_and_current(const SpatialGrid& sg,
                                             double x,
                                             double weight,
                                             double vx)
{
    if (x < 0.0 || x >= Param::Lx) return;

    const double s = x / sg.dx - 0.5;
    const int i0 = static_cast<int>(std::floor(s));
    const double frac = s - i0;
    const double density0 = weight * (1.0 - frac) / sg.dx;
    const double density1 = weight * frac / sg.dx;
    add_source_to_cell(sg, i0, density0, -Const::qe * vx * density0);
    add_source_to_cell(sg, i0 + 1, density1, -Const::qe * vx * density1);
}

void BeamPIC::add_density_to(const SpatialGrid& sg,
                             std::vector<double>& local,
                             std::vector<double>& send_left,
                             std::vector<double>& send_right,
                             double x, double weight)
{
    if (x < 0.0 || x > Param::Lx) return;

    const double s = x / sg.dx - 0.5;
    const int i0 = static_cast<int>(std::floor(s));
    const double frac = s - i0;
    add_number_to_cell(sg, local, send_left, send_right,
                       i0, weight * (1.0 - frac) / sg.dx);
    add_number_to_cell(sg, local, send_left, send_right,
                       i0 + 1, weight * frac / sg.dx);
}

void BeamPIC::add_path_density_delta(const SpatialGrid& sg,
                                     double x0, double x1, double weight)
{
    add_path_density_delta_to(sg, path_density_delta,
                              path_send_left_density_,
                              path_send_right_density_,
                              x0, x1, weight);
}

void BeamPIC::add_path_density_delta_to(const SpatialGrid& sg,
                                        std::vector<double>& local,
                                        std::vector<double>& send_left,
                                        std::vector<double>& send_right,
                                        double x0, double x1, double weight)
{
    if (x0 < 0.0 || x0 > Param::Lx) return;
    if (x1 < 0.0 || x1 > Param::Lx) return;

    const double s0 = x0 / sg.dx - 0.5;
    const int i00 = static_cast<int>(std::floor(s0));
    const double f0 = s0 - i00;
    add_number_to_cell(sg, local, send_left, send_right,
                       i00, -weight * (1.0 - f0) / sg.dx);
    add_number_to_cell(sg, local, send_left, send_right,
                       i00 + 1, -weight * f0 / sg.dx);

    const double s1 = x1 / sg.dx - 0.5;
    const int i10 = static_cast<int>(std::floor(s1));
    const double f1 = s1 - i10;
    add_number_to_cell(sg, local, send_left, send_right,
                       i10, weight * (1.0 - f1) / sg.dx);
    add_number_to_cell(sg, local, send_left, send_right,
                       i10 + 1, weight * f1 / sg.dx);
}

void BeamPIC::reset_continuity_exchange_buffers(const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    if (source_send_left_density_.size() != static_cast<size_t>(ng)) {
        source_send_left_density_.assign(ng, 0.0);
        source_send_right_density_.assign(ng, 0.0);
        source_recv_left_density_.assign(ng, 0.0);
        source_recv_right_density_.assign(ng, 0.0);
        source_send_left_current_.assign(ng, 0.0);
        source_send_right_current_.assign(ng, 0.0);
        source_recv_left_current_.assign(ng, 0.0);
        source_recv_right_current_.assign(ng, 0.0);
        path_send_left_density_.assign(ng, 0.0);
        path_send_right_density_.assign(ng, 0.0);
        path_recv_left_density_.assign(ng, 0.0);
        path_recv_right_density_.assign(ng, 0.0);
    } else {
        std::fill(source_send_left_density_.begin(), source_send_left_density_.end(), 0.0);
        std::fill(source_send_right_density_.begin(), source_send_right_density_.end(), 0.0);
        std::fill(source_recv_left_density_.begin(), source_recv_left_density_.end(), 0.0);
        std::fill(source_recv_right_density_.begin(), source_recv_right_density_.end(), 0.0);
        std::fill(source_send_left_current_.begin(), source_send_left_current_.end(), 0.0);
        std::fill(source_send_right_current_.begin(), source_send_right_current_.end(), 0.0);
        std::fill(source_recv_left_current_.begin(), source_recv_left_current_.end(), 0.0);
        std::fill(source_recv_right_current_.begin(), source_recv_right_current_.end(), 0.0);
        std::fill(path_send_left_density_.begin(), path_send_left_density_.end(), 0.0);
        std::fill(path_send_right_density_.begin(), path_send_right_density_.end(), 0.0);
        std::fill(path_recv_left_density_.begin(), path_recv_left_density_.end(), 0.0);
        std::fill(path_recv_right_density_.begin(), path_recv_right_density_.end(), 0.0);
    }
}

void BeamPIC::exchange_continuity_contributions(const SpatialGrid& sg,
                                                int mpi_rank, int mpi_size)
{
    const int ng = sg.nghost;
    if (ng <= 0) return;
    if (mpi_size == 1) {
        return;
    }

    const int left = mpi_rank - 1;
    const int right = mpi_rank + 1;
    MPI_Request reqs[12];
    int nreq = 0;

    if (left >= 0) {
        MPI_Isend(source_send_left_density_.data(), ng, MPI_DOUBLE,
                  left, 501, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(source_recv_left_density_.data(), ng, MPI_DOUBLE,
                  left, 502, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(source_send_left_current_.data(), ng, MPI_DOUBLE,
                  left, 503, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(source_recv_left_current_.data(), ng, MPI_DOUBLE,
                  left, 504, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(path_send_left_density_.data(), ng, MPI_DOUBLE,
                  left, 505, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(path_recv_left_density_.data(), ng, MPI_DOUBLE,
                  left, 506, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(source_send_right_density_.data(), ng, MPI_DOUBLE,
                  right, 502, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(source_recv_right_density_.data(), ng, MPI_DOUBLE,
                  right, 501, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(source_send_right_current_.data(), ng, MPI_DOUBLE,
                  right, 504, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(source_recv_right_current_.data(), ng, MPI_DOUBLE,
                  right, 503, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(path_send_right_density_.data(), ng, MPI_DOUBLE,
                  right, 506, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(path_recv_right_density_.data(), ng, MPI_DOUBLE,
                  right, 505, MPI_COMM_WORLD, &reqs[nreq++]);
    }

    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (left >= 0) {
        for (int g = 0; g < ng && g < sg.nx_local; ++g) {
            const size_t slot = static_cast<size_t>(g);
            source_density_delta[slot] += source_recv_left_density_[slot];
            source_current_x[slot] += source_recv_left_current_[slot];
            path_density_delta[slot] += path_recv_left_density_[slot];
        }
    }
    if (right < mpi_size) {
        for (int g = 0; g < ng; ++g) {
            const int il = sg.nx_local - ng + g;
            if (il < 0 || il >= sg.nx_local) continue;
            const size_t dst = static_cast<size_t>(il);
            const size_t src = static_cast<size_t>(g);
            source_density_delta[dst] += source_recv_right_density_[src];
            source_current_x[dst] += source_recv_right_current_[src];
            path_density_delta[dst] += path_recv_right_density_[src];
        }
    }
}

void BeamPIC::deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    const int nthreads = std::max(1, omp_get_max_threads());
    if (thread_density_.size() != static_cast<size_t>(nthreads)) {
        thread_density_.resize(nthreads);
        thread_current_.resize(nthreads);
        thread_send_left_density_.resize(nthreads);
        thread_send_right_density_.resize(nthreads);
        thread_send_left_current_.resize(nthreads);
        thread_send_right_current_.resize(nthreads);
    }
    if (thread_current_.size() != static_cast<size_t>(nthreads)) {
        thread_current_.resize(nthreads);
        thread_send_left_current_.resize(nthreads);
        thread_send_right_current_.resize(nthreads);
    }
    for (int t = 0; t < nthreads; ++t) {
        if (thread_density_[t].size() != static_cast<size_t>(sg.nx_local) ||
            thread_current_[t].size() != static_cast<size_t>(sg.nx_local)) {
            thread_density_[t].assign(sg.nx_local, 0.0);
            thread_current_[t].assign(sg.nx_local, 0.0);
        } else {
            std::fill(thread_density_[t].begin(), thread_density_[t].end(), 0.0);
            std::fill(thread_current_[t].begin(), thread_current_[t].end(), 0.0);
        }
        thread_send_left_density_[t] = 0.0;
        thread_send_right_density_[t] = 0.0;
        thread_send_left_current_[t] = 0.0;
        thread_send_right_current_[t] = 0.0;
    }

    double send_left_density = 0.0;
    double send_right_density = 0.0;
    double send_left_current = 0.0;
    double send_right_current = 0.0;
    const long long np = static_cast<long long>(particles.size());

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::vector<double>& local_density = thread_density_[tid];
        std::vector<double>& local_current = thread_current_[tid];
        double local_send_left_density = 0.0;
        double local_send_right_density = 0.0;
        double local_send_left_current = 0.0;
        double local_send_right_current = 0.0;

        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            const BeamParticle& p = particles[static_cast<size_t>(i)];
            double s = p.x / sg.dx - 0.5;
            int i0 = static_cast<int>(std::floor(s));
            double frac = s - i0;
            const double pnorm = p.px / (Const::me * Const::c);
            const double gamma = std::sqrt(1.0 + pnorm * pnorm);
            const double vx = p.px / (gamma * Const::me);

            double contribution = p.weight * (1.0 - frac) / sg.dx;
            double current_contribution = -Const::qe * vx * contribution;
            int target_ig = i0;
            int il = target_ig - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                local_density[il] += contribution;
                local_current[il] += current_contribution;
            } else if (target_ig < sg.ix_start && mpi_rank > 0) {
                local_send_left_density += contribution;
                local_send_left_current += current_contribution;
            } else if (target_ig >= sg.ix_start + sg.nx_local &&
                       mpi_rank + 1 < mpi_size) {
                local_send_right_density += contribution;
                local_send_right_current += current_contribution;
            }

            contribution = p.weight * frac / sg.dx;
            current_contribution = -Const::qe * vx * contribution;
            target_ig = i0 + 1;
            il = target_ig - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                local_density[il] += contribution;
                local_current[il] += current_contribution;
            } else if (target_ig < sg.ix_start && mpi_rank > 0) {
                local_send_left_density += contribution;
                local_send_left_current += current_contribution;
            } else if (target_ig >= sg.ix_start + sg.nx_local &&
                       mpi_rank + 1 < mpi_size) {
                local_send_right_density += contribution;
                local_send_right_current += current_contribution;
            }
        }

        thread_send_left_density_[tid] = local_send_left_density;
        thread_send_right_density_[tid] = local_send_right_density;
        thread_send_left_current_[tid] = local_send_left_current;
        thread_send_right_current_[tid] = local_send_right_current;
    }

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double density_sum = 0.0;
        double current_sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            density_sum += thread_density_[t][ix];
            current_sum += thread_current_[t][ix];
        }
        density[ix] = density_sum;
        current_x[ix] = current_sum;
    }

    for (int t = 0; t < nthreads; ++t) {
        send_left_density += thread_send_left_density_[t];
        send_right_density += thread_send_right_density_[t];
        send_left_current += thread_send_left_current_[t];
        send_right_current += thread_send_right_current_[t];
    }

    double send_left_payload[2] = { send_left_density, send_left_current };
    double send_right_payload[2] = { send_right_density, send_right_current };
    double recv_left_payload[2] = { 0.0, 0.0 };
    double recv_right_payload[2] = { 0.0, 0.0 };
    int left = mpi_rank - 1;
    int right = mpi_rank + 1;

    MPI_Request reqs[4];
    int nreq = 0;
    if (mpi_size == 1) {
        return;
    }
    if (left >= 0) {
        MPI_Isend(send_left_payload, 2, MPI_DOUBLE, left, 401, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_left_payload, 2, MPI_DOUBLE, left, 402, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(send_right_payload, 2, MPI_DOUBLE, right, 402, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_right_payload, 2, MPI_DOUBLE, right, 401, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (left >= 0 && sg.nx_local > 0) {
        density[0] += recv_left_payload[0];
        current_x[0] += recv_left_payload[1];
    }
    if (right < mpi_size && sg.nx_local > 0) {
        density[sg.nx_local - 1] += recv_right_payload[0];
        current_x[sg.nx_local - 1] += recv_right_payload[1];
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

    double local_path_number = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_path_number += path_density_delta[static_cast<size_t>(ix)] * sg.dx;
    }

    if (all_path_numbers_.size() != static_cast<size_t>(mpi_size)) {
        all_path_numbers_.assign(static_cast<size_t>(mpi_size), 0.0);
    }
    MPI_Allgather(&local_path_number, 1, MPI_DOUBLE,
                  all_path_numbers_.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    double global_left_boundary_number = 0.0;
    MPI_Allreduce(&left_boundary_number_flux_, &global_left_boundary_number,
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double number_flux = global_left_boundary_number / elapsed_dt;
    for (int r = 0; r < mpi_rank; ++r) {
        number_flux -= all_path_numbers_[static_cast<size_t>(r)] / elapsed_dt;
    }

    current_face_x[0] = -Const::qe * number_flux;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        number_flux -= path_density_delta[static_cast<size_t>(ix)]
                     * sg.dx / elapsed_dt;
        current_face_x[static_cast<size_t>(ix + 1)] = -Const::qe * number_flux;
    }

    for (int ix = 0; ix < sg.nx_local; ++ix) {
        current_x[static_cast<size_t>(ix)] =
            0.5 * (current_face_x[static_cast<size_t>(ix)]
                 + current_face_x[static_cast<size_t>(ix + 1)]);
    }

    last_continuity_l1_error_ = 0.0;
    last_continuity_linf_error_ = 0.0;
    const size_t n = std::min(density.size(), density_step_start.size());
    for (size_t ix = 0; ix < n; ++ix) {
        const double expected =
            density_step_start[ix] + path_density_delta[ix]
            + source_density_delta[ix];
        const double err = density[ix] - expected;
        const double abs_err = std::fabs(err);
        last_continuity_l1_error_ += abs_err * sg.dx;
        last_continuity_linf_error_ =
            std::max(last_continuity_linf_error_, abs_err);
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
