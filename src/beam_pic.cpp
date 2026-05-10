#include "beam_pic.h"
#include <algorithm>
#include <cmath>
#include <mpi.h>

BeamPIC::BeamPIC()
    : injection_remainder_(0.0)
{}

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

    double x0 = 0.5 * sg.dx;
    for (int i = 0; i < n_new; ++i) {
        BeamParticle p;
        p.x = x0;
        p.px = Param::beam_p0;
        p.weight = Param::beam_macro_weight;
        particles.push_back(p);
    }
}

void BeamPIC::push(const SpatialGrid& sg, const EMFields& fields, double dt,
                   int mpi_rank, int mpi_size)
{
    int ng = sg.nghost;
    for (size_t i = 0; i < particles.size(); ++i) {
        BeamParticle& p = particles[i];
        int ig = static_cast<int>(p.x / sg.dx);
        int il = ig - sg.ix_start;
        if (il < 0) il = 0;
        if (il >= sg.nx_local) il = sg.nx_local - 1;
        double ex = fields.Ex[ng + il];

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

void BeamPIC::deposit_density(const SpatialGrid& sg)
{
    std::fill(density.begin(), density.end(), 0.0);
    for (size_t i = 0; i < particles.size(); ++i) {
        int ig = static_cast<int>(particles[i].x / sg.dx);
        int il = ig - sg.ix_start;
        if (il >= 0 && il < sg.nx_local) {
            density[il] += particles[i].weight / sg.dx;
        }
    }
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
