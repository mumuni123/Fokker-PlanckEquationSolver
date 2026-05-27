#include "beam_pic.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <mpi.h>
#include <omp.h>

BeamPIC::BeamPIC()
    : injection_remainder_(0.0)
{}

namespace {
size_t initial_particle_capacity(const SpatialGrid& sg)
{
    const double active_time = std::max(
        0.0, std::min(Param::t_end, Param::t_inject_end) - Param::t_inject_start);
    const double active_length = std::min(Param::Lx, Param::beam_v0 * active_time);
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
    double s = x / sg.dx - 0.5;
    int i0 = static_cast<int>(std::floor(s));
    double frac = s - i0;
    double w0 = 1.0 - frac;
    double w1 = frac;
    double ex = 0.0;

    int target_ig = std::max(0, std::min(i0, sg.nx_global - 1));
    int il = target_ig - sg.ix_start;
    int ix = il + sg.nghost;
    if (ix >= 0 && ix < sg.nx_total) {
        ex += w0 * fields.Ex[ix];
    }

    target_ig = std::max(0, std::min(i0 + 1, sg.nx_global - 1));
    il = target_ig - sg.ix_start;
    ix = il + sg.nghost;
    if (ix >= 0 && ix < sg.nx_total) {
        ex += w1 * fields.Ex[ix];
    }

    return ex;
}
}

void BeamPIC::init(const SpatialGrid& sg)
{
    particles.clear();
    density.assign(sg.nx_local, 0.0);
    current_x.assign(sg.nx_local, 0.0);
    injection_remainder_ = 0.0;

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

void BeamPIC::inject(const SpatialGrid& sg, double dt, double time, int mpi_rank)
{
    const double step_start = time - dt;
    const double active_start = std::max(step_start, Param::t_inject_start);
    const double active_end = std::min(time, Param::t_inject_end);
    const double active_dt = active_end - active_start;
    if (active_dt <= 0.0) return;
    if (mpi_rank != 0 || sg.ix_start != 0) return;

    double physical_per_area = Param::densb * Param::beam_v0 * active_dt
                             + injection_remainder_;
    int n_new = static_cast<int>(physical_per_area / Param::beam_macro_weight);
    injection_remainder_ = physical_per_area - n_new * Param::beam_macro_weight;
    if (n_new <= 0) return;

    particles.reserve(particles.size() + static_cast<size_t>(n_new));
    for (int i = 0; i < n_new; ++i) {
        // Boundary-source particles start outside and cross x_min during this step.
        const double crossing_time =
            active_start + (i + 0.5) * active_dt / static_cast<double>(n_new);
        BeamParticle p;
        p.x = -Param::beam_v0 * (crossing_time - step_start);
        p.px = Param::beam_p0;
        p.weight = Param::beam_macro_weight;
        particles.push_back(p);
    }
}

void BeamPIC::push(const SpatialGrid& sg, const EMFields& fields, double dt,
                   int mpi_rank, int mpi_size)
{
    const double x_left = sg.ix_start * sg.dx;
    const double x_right = (sg.ix_start + sg.nx_local) * sg.dx;
    const double domain_left = 0.0;
    const double domain_right = Param::Lx;
    const long long np = static_cast<long long>(particles.size());
    const int nthreads = std::max(1, omp_get_max_threads());

    if (thread_keep_.size() != static_cast<size_t>(nthreads)) {
        thread_keep_.resize(nthreads);
        thread_send_left_.resize(nthreads);
        thread_send_right_.resize(nthreads);
    }
    const size_t reserve_each = static_cast<size_t>(np / nthreads + 1);
    for (int t = 0; t < nthreads; ++t) {
        thread_keep_[t].clear();
        thread_send_left_[t].clear();
        thread_send_right_[t].clear();
        thread_keep_[t].reserve(reserve_each);
        thread_send_left_[t].reserve(std::max<size_t>(1, reserve_each / 16));
        thread_send_right_[t].reserve(std::max<size_t>(1, reserve_each / 16));
    }

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::vector<BeamParticle>& local_keep = thread_keep_[tid];
        std::vector<BeamParticle>& local_left = thread_send_left_[tid];
        std::vector<BeamParticle>& local_right = thread_send_right_[tid];

        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            BeamParticle p = particles[static_cast<size_t>(i)];
            const double ex = gather_cic_ex(p.x, sg, fields);

            p.px += (-Const::qe) * ex * dt;
            const double pnorm = p.px / (Const::me * Const::c);
            const double gamma = std::sqrt(1.0 + pnorm * pnorm);
            const double vx = p.px / (gamma * Const::me);
            p.x += vx * dt;

            if (p.x < domain_left || p.x >= domain_right) {
                continue;
            }
            if (p.x < x_left && mpi_rank > 0) {
                local_left.push_back(p);
            } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
                local_right.push_back(p);
            } else if (p.x >= x_left && p.x < x_right) {
                local_keep.push_back(p);
            }
        }
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
}

void BeamPIC::exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    (void)sg;

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
    if (left >= 0) {
        if (send_left_count > 0) {
            MPI_Isend(send_left_.data(), particle_byte_count(send_left_count), MPI_BYTE,
                      left, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
        }
        if (recv_left_count > 0) {
            MPI_Irecv(recv_left_.data(), particle_byte_count(recv_left_count), MPI_BYTE,
                      left, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
        }
    }
    if (right < mpi_size) {
        if (send_right_count > 0) {
            MPI_Isend(send_right_.data(), particle_byte_count(send_right_count), MPI_BYTE,
                      right, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
        }
        if (recv_right_count > 0) {
            MPI_Irecv(recv_right_.data(), particle_byte_count(recv_right_count), MPI_BYTE,
                      right, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
        }
    }
    if (nreq > 0) MPI_Waitall(nreq, data_reqs, MPI_STATUSES_IGNORE);

    particles.reserve(particles.size() + recv_left_count + recv_right_count);
    particles.insert(particles.end(), recv_left_.begin(), recv_left_.end());
    particles.insert(particles.end(), recv_right_.begin(), recv_right_.end());
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
            if (target_ig >= 0 && target_ig < sg.nx_global) {
                if (il >= 0 && il < sg.nx_local) {
                    local_density[il] += contribution;
                    local_current[il] += current_contribution;
                } else if (target_ig < sg.ix_start) {
                    local_send_left_density += contribution;
                    local_send_left_current += current_contribution;
                } else {
                    local_send_right_density += contribution;
                    local_send_right_current += current_contribution;
                }
            }

            contribution = p.weight * frac / sg.dx;
            current_contribution = -Const::qe * vx * contribution;
            target_ig = i0 + 1;
            il = target_ig - sg.ix_start;
            if (target_ig >= 0 && target_ig < sg.nx_global) {
                if (il >= 0 && il < sg.nx_local) {
                    local_density[il] += contribution;
                    local_current[il] += current_contribution;
                } else if (target_ig < sg.ix_start) {
                    local_send_left_density += contribution;
                    local_send_left_current += current_contribution;
                } else {
                    local_send_right_density += contribution;
                    local_send_right_current += current_contribution;
                }
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
        double pnorm = particles[i].px / (Const::me * Const::c);
        double gamma = std::sqrt(1.0 + pnorm * pnorm);
        total += particles[i].weight * (gamma - 1.0) * Const::me * Const::c * Const::c;
    }
    return total;
}
