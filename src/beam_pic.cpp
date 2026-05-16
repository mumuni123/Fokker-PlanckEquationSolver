#include "beam_pic.h"
#include <algorithm>
#include <cmath>
#include <mpi.h>

BeamPIC::BeamPIC()
    : injection_remainder_(0.0)
{}

namespace {
double cic_weight_to_cell(double x, int ig, double dx)
{
    double x_cell = (ig + 0.5) * dx;
    double r = std::fabs(x - x_cell) / dx;
    return (r < 1.0) ? (1.0 - r) : 0.0;
}

double gather_cic_ex(double x, const SpatialGrid& sg, const EMFields& fields)
{
    double s = x / sg.dx - 0.5;
    int i0 = static_cast<int>(std::floor(s));
    int i1 = i0 + 1;
    double ex = 0.0;

    for (int n = 0; n < 2; ++n) {
        int ig = (n == 0) ? i0 : i1;
        int target_ig = std::max(0, std::min(ig, sg.nx_global - 1));
        double w = cic_weight_to_cell(x, ig, sg.dx);
        int il = target_ig - sg.ix_start;
        int ix = il + sg.nghost;
        if (ix >= 0 && ix < sg.nx_total) {
            ex += w * fields.Ex[ix];
        }
    }

    return ex;
}
}

void BeamPIC::init(const SpatialGrid& sg)
{
    particles.clear();
    density.assign(sg.nx_local, 0.0);
    injection_remainder_ = 0.0;
}

void BeamPIC::inject(const SpatialGrid& sg, double dt, double time, int mpi_rank)
{
    if (time > Param::t_inject_end) return;
    if (mpi_rank != 0 || sg.ix_start != 0) return;

    double physical_per_area = Param::densb * Param::beam_v0 * dt + injection_remainder_;
    int n_new = static_cast<int>(physical_per_area / Param::beam_macro_weight);
    injection_remainder_ = physical_per_area - n_new * Param::beam_macro_weight;

    double injection_length = std::min(Param::beam_v0 * dt, sg.dx);
    for (int i = 0; i < n_new; ++i) {
        BeamParticle p;
        p.x = (i + 0.5) * injection_length / std::max(n_new, 1);
        p.px = Param::beam_p0;
        p.weight = Param::beam_macro_weight;
        particles.push_back(p);
    }
}

void BeamPIC::push(const SpatialGrid& sg, const EMFields& fields, double dt,
                   int mpi_rank, int mpi_size)
{
    for (size_t i = 0; i < particles.size(); ++i) {
        BeamParticle& p = particles[i];
        double ex = gather_cic_ex(p.x, sg, fields);

        p.px += (-Const::qe) * ex * dt;
        double gamma = std::sqrt(1.0 + (p.px / (Const::me * Const::c)) *
                                       (p.px / (Const::me * Const::c)));
        double vx = p.px / (gamma * Const::me);
        p.x += vx * dt;
    }

    // Remove particles that left the global physical domain.
    double domain_left = 0.0;
    double domain_right = Param::Lx;
    particles.erase(std::remove_if(particles.begin(), particles.end(),
        [domain_left, domain_right](const BeamParticle& p) {
            return p.x < domain_left || p.x >= domain_right;
        }), particles.end());

    exchange_particles(sg, mpi_rank, mpi_size);
}

void BeamPIC::exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    double x_left = sg.ix_start * sg.dx;
    double x_right = (sg.ix_start + sg.nx_local) * sg.dx;

    std::vector<BeamParticle> send_left;
    std::vector<BeamParticle> send_right;
    std::vector<BeamParticle> keep;
    keep.reserve(particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        const BeamParticle& p = particles[i];
        if (p.x < x_left && mpi_rank > 0) {
            send_left.push_back(p);
        } else if (p.x >= x_right && mpi_rank + 1 < mpi_size) {
            send_right.push_back(p);
        } else if (p.x >= x_left && p.x < x_right) {
            keep.push_back(p);
        }
    }
    particles.swap(keep);

    int left = mpi_rank - 1;
    int right = mpi_rank + 1;
    int send_left_count = static_cast<int>(send_left.size());
    int send_right_count = static_cast<int>(send_right.size());
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

    std::vector<double> send_left_buf(3 * send_left.size());
    std::vector<double> send_right_buf(3 * send_right.size());
    std::vector<double> recv_left_buf(3 * recv_left_count);
    std::vector<double> recv_right_buf(3 * recv_right_count);

    for (size_t i = 0; i < send_left.size(); ++i) {
        send_left_buf[3*i + 0] = send_left[i].x;
        send_left_buf[3*i + 1] = send_left[i].px;
        send_left_buf[3*i + 2] = send_left[i].weight;
    }
    for (size_t i = 0; i < send_right.size(); ++i) {
        send_right_buf[3*i + 0] = send_right[i].x;
        send_right_buf[3*i + 1] = send_right[i].px;
        send_right_buf[3*i + 2] = send_right[i].weight;
    }

    MPI_Request data_reqs[4];
    nreq = 0;
    if (left >= 0) {
        MPI_Isend(send_left_buf.data(), 3 * send_left_count, MPI_DOUBLE, left, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
        MPI_Irecv(recv_left_buf.data(), 3 * recv_left_count, MPI_DOUBLE, left, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(send_right_buf.data(), 3 * send_right_count, MPI_DOUBLE, right, 304, MPI_COMM_WORLD, &data_reqs[nreq++]);
        MPI_Irecv(recv_right_buf.data(), 3 * recv_right_count, MPI_DOUBLE, right, 303, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, data_reqs, MPI_STATUSES_IGNORE);

    for (int i = 0; i < recv_left_count; ++i) {
        BeamParticle p;
        p.x = recv_left_buf[3*i + 0];
        p.px = recv_left_buf[3*i + 1];
        p.weight = recv_left_buf[3*i + 2];
        particles.push_back(p);
    }
    for (int i = 0; i < recv_right_count; ++i) {
        BeamParticle p;
        p.x = recv_right_buf[3*i + 0];
        p.px = recv_right_buf[3*i + 1];
        p.weight = recv_right_buf[3*i + 2];
        particles.push_back(p);
    }
}

void BeamPIC::deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size)
{
    std::fill(density.begin(), density.end(), 0.0);
    double send_left = 0.0;
    double send_right = 0.0;

    for (size_t i = 0; i < particles.size(); ++i) {
        const BeamParticle& p = particles[i];
        double s = p.x / sg.dx - 0.5;
        int i0 = static_cast<int>(std::floor(s));

        for (int n = 0; n < 2; ++n) {
            int ig = i0 + n;
            int target_ig = std::max(0, std::min(ig, sg.nx_global - 1));
            double contribution = p.weight * cic_weight_to_cell(p.x, ig, sg.dx) / sg.dx;
            if (contribution == 0.0) continue;

            int il = target_ig - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                density[il] += contribution;
            } else if (target_ig < sg.ix_start) {
                send_left += contribution;
            } else {
                send_right += contribution;
            }
        }
    }

    double recv_left = 0.0;
    double recv_right = 0.0;
    int left = mpi_rank - 1;
    int right = mpi_rank + 1;

    MPI_Request reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(&send_left, 1, MPI_DOUBLE, left, 401, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, left, 402, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right, 1, MPI_DOUBLE, right, 402, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, right, 401, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (left >= 0 && sg.nx_local > 0) density[0] += recv_left;
    if (right < mpi_size && sg.nx_local > 0) density[sg.nx_local - 1] += recv_right;
}

double BeamPIC::total_particle_number(const SpatialGrid& sg) const
{
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
