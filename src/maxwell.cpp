#include "maxwell.h"
#include "species.h"
#include <algorithm>
#include <mpi.h>
#include <vector>

void EMFields::init(const SpatialGrid& sg)
{
    nx_total = sg.nx_total;
    dx = sg.dx;
    Ex.assign(nx_total, 0.0);
    rho.assign(nx_total, 0.0);
}

void EMFields::zero_currents()
{
    std::fill(rho.begin(), rho.end(), 0.0);
}

void EMFields::accumulate_moments(const Species& sp)
{
    int ng = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;
    for (int ix = 0; ix < nxl; ++ix) {
        rho[ix + ng] += sp.charge_density[ix];
    }
}

void EMFields::set_charge_density(const Species& electrons,
                                  const std::vector<double>& beam_density)
{
    int ng = electrons.sgrid->nghost;
    int nxl = electrons.sgrid->nx_local;
    std::fill(rho.begin(), rho.end(), 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        double ion_charge_density = Const::qe * Param::dens;
        double electron_charge_density = -Const::qe * electrons.number_density[ix];
        double beam_charge_density = -Const::qe * beam_density[ix];
        rho[ix + ng] = ion_charge_density + electron_charge_density + beam_charge_density;
    }
}

void EMFields::solve_poisson(int mpi_rank, int mpi_size)
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;

    // Non-periodic electrostatic solve. Use Ex(left boundary) = 0 as the
    // integration reference and enforce dEx/dx = rho/eps0 directly.
    double local_prefix_end = 0.0;
    for (int i = ng; i < ng + nxl; ++i) {
        local_prefix_end += rho[i] * dx / Const::eps0;
    }

    std::vector<double> rank_offsets(mpi_size, 0.0);
    MPI_Allgather(&local_prefix_end, 1, MPI_DOUBLE,
                  rank_offsets.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    double offset = 0.0;
    for (int r = 0; r < mpi_rank; ++r) offset += rank_offsets[r];

    double e = offset;
    for (int i = ng; i < ng + nxl; ++i) {
        e += 0.5 * rho[i] * dx / Const::eps0;
        Ex[i] = e;
        e += 0.5 * rho[i] * dx / Const::eps0;
    }

    exchange_ghosts(mpi_rank, mpi_size);
}

void EMFields::exchange_ghosts(int mpi_rank, int mpi_size)
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;
    int left = mpi_rank - 1;
    int right = mpi_rank + 1;

    std::vector<double> send_l(ng), send_r(ng), recv_l(ng), recv_r(ng);
    for (int g = 0; g < ng; ++g) {
        send_l[g] = Ex[ng + g];
        send_r[g] = Ex[ng + nxl - ng + g];
    }

    MPI_Request reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(send_l.data(), ng, MPI_DOUBLE, left, 201, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_l.data(), ng, MPI_DOUBLE, left, 202, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(send_r.data(), ng, MPI_DOUBLE, right, 202, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_r.data(), ng, MPI_DOUBLE, right, 201, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) {
        MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
    }

    for (int g = 0; g < ng; ++g) {
        Ex[g] = (left >= 0) ? recv_l[g] : Ex[ng];
        Ex[ng + nxl + g] = (right < mpi_size) ? recv_r[g] : Ex[ng + nxl - 1];
    }
}

double EMFields::total_energy() const
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;
    double energy = 0.0;
    for (int i = ng; i < ng + nxl; ++i) {
        energy += 0.5 * Const::eps0 * Ex[i] * Ex[i];
    }
    return energy * dx;
}
