#include "maxwell.h"
#include "species.h"
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <vector>

namespace {
void prepare_mpi_layout(EMFields& fields, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;

    if (fields.counts_mpi_size == mpi_size &&
        fields.counts_nx_local == nxl &&
        fields.global_rhs.size() == static_cast<size_t>(Param::nx)) {
        return;
    }

    fields.counts_mpi_size = mpi_size;
    fields.counts_nx_local = nxl;
    fields.counts.assign(static_cast<size_t>(mpi_size), 0);
    fields.displs.assign(static_cast<size_t>(mpi_size), 0);
    MPI_Allgather(&nxl, 1, MPI_INT, fields.counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int r = 1; r < mpi_size; ++r) {
        fields.displs[static_cast<size_t>(r)] =
            fields.displs[static_cast<size_t>(r - 1)] +
            fields.counts[static_cast<size_t>(r - 1)];
    }

    fields.local_rhs.assign(static_cast<size_t>(nxl), 0.0);
    fields.global_rhs.assign(static_cast<size_t>(Param::nx), 0.0);
    fields.global_ex.assign(static_cast<size_t>(Param::nx), 0.0);
    fields.global_phi.assign(static_cast<size_t>(Param::nx), 0.0);
}

void gather_rho(EMFields& fields, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;
    prepare_mpi_layout(fields, mpi_size);

    for (int ix = 0; ix < nxl; ++ix) {
        fields.local_rhs[static_cast<size_t>(ix)] = fields.rho[ng + ix];
    }

    MPI_Gatherv(fields.local_rhs.data(), nxl, MPI_DOUBLE,
                fields.global_rhs.data(), fields.counts.data(),
                fields.displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void remove_global_mean_charge_source(EMFields& fields)
{
    if (!Param::poisson_remove_global_mean_charge ||
        fields.global_rhs.empty()) {
        return;
    }

    double sum = 0.0;
    for (size_t i = 0; i < fields.global_rhs.size(); ++i) {
        sum += fields.global_rhs[i];
    }
    const double mean = sum / static_cast<double>(fields.global_rhs.size());
    for (size_t i = 0; i < fields.global_rhs.size(); ++i) {
        fields.global_rhs[i] -= mean;
    }
}

void compute_dirichlet_poisson(EMFields& fields,
                               bool compute_ex,
                               bool compute_phi)
{
    const int n = Param::nx;
    if (n <= 0) return;

    fields.tri_work_a.assign(static_cast<size_t>(n), -1.0);
    fields.tri_work_b.assign(static_cast<size_t>(n), 2.0);
    fields.tri_work_c.assign(static_cast<size_t>(n), -1.0);
    fields.tri_work_d.assign(static_cast<size_t>(n), 0.0);
    fields.tri_work_a[0] = 0.0;
    fields.tri_work_c[static_cast<size_t>(n - 1)] = 0.0;

    const double rhs_scale = fields.dx * fields.dx / Const::eps0;
    for (int i = 0; i < n; ++i) {
        fields.tri_work_d[static_cast<size_t>(i)] =
            fields.global_rhs[static_cast<size_t>(i)] * rhs_scale;
    }

    for (int i = 1; i < n; ++i) {
        const size_t im = static_cast<size_t>(i - 1);
        const size_t ii = static_cast<size_t>(i);
        const double m = fields.tri_work_a[ii] / fields.tri_work_b[im];
        fields.tri_work_b[ii] -= m * fields.tri_work_c[im];
        fields.tri_work_d[ii] -= m * fields.tri_work_d[im];
    }

    fields.global_phi.assign(static_cast<size_t>(n), 0.0);
    fields.global_phi[static_cast<size_t>(n - 1)] =
        fields.tri_work_d[static_cast<size_t>(n - 1)] /
        fields.tri_work_b[static_cast<size_t>(n - 1)];
    for (int i = n - 2; i >= 0; --i) {
        const size_t ii = static_cast<size_t>(i);
        fields.global_phi[ii] =
            (fields.tri_work_d[ii] -
             fields.tri_work_c[ii] * fields.global_phi[ii + 1]) /
            fields.tri_work_b[ii];
    }

    if (compute_ex) {
        fields.global_ex.assign(static_cast<size_t>(n), 0.0);
        if (n > 1) {
            fields.global_ex[0] = -fields.global_phi[1] / fields.dx;
            for (int i = 1; i < n - 1; ++i) {
                fields.global_ex[static_cast<size_t>(i)] =
                    -(fields.global_phi[static_cast<size_t>(i + 1)] -
                      fields.global_phi[static_cast<size_t>(i - 1)]) /
                    (2.0 * fields.dx);
            }
            fields.global_ex[static_cast<size_t>(n - 1)] =
                fields.global_phi[static_cast<size_t>(n - 2)] / fields.dx;
        }
    }

    if (!compute_phi) {
        std::fill(fields.global_phi.begin(), fields.global_phi.end(), 0.0);
    }
}

void exchange_scalar_ghosts(EMFields& fields,
                            std::vector<double>& a,
                            int tag_base,
                            int mpi_rank,
                            int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;

    if (fields.send_left.size() != static_cast<size_t>(ng)) {
        fields.send_left.resize(ng);
        fields.send_right.resize(ng);
        fields.recv_left.resize(ng);
        fields.recv_right.resize(ng);
    }

    const int left = mpi_rank - 1;
    const int right = mpi_rank + 1;

    if (mpi_size == 1) {
        for (int g = 0; g < ng; ++g) {
            a[g] = a[ng];
            a[ng + nxl + g] = a[ng + nxl - 1];
        }
        return;
    }

    for (int g = 0; g < ng; ++g) {
        fields.send_left[g] = a[ng + g];
        fields.send_right[g] = a[ng + nxl - ng + g];
    }

    MPI_Request reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(fields.send_left.data(), ng, MPI_DOUBLE, left, tag_base,
                  MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(fields.recv_left.data(), ng, MPI_DOUBLE, left, tag_base + 1,
                  MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(fields.send_right.data(), ng, MPI_DOUBLE, right, tag_base + 1,
                  MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(fields.recv_right.data(), ng, MPI_DOUBLE, right, tag_base,
                  MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    for (int g = 0; g < ng; ++g) {
        a[g] = (left >= 0) ? fields.recv_left[g] : a[ng];
        a[ng + nxl + g] =
            (right < mpi_size) ? fields.recv_right[g] : a[ng + nxl - 1];
    }
}
}

void EMFields::init(const SpatialGrid& sg)
{
    nx_total = sg.nx_total;
    counts_mpi_size = 0;
    counts_nx_local = -1;
    dx = sg.dx;
    Ex.assign(nx_total, 0.0);
    phi.assign(nx_total, 0.0);
    rho.assign(nx_total, 0.0);
    send_left.assign(sg.nghost, 0.0);
    send_right.assign(sg.nghost, 0.0);
    recv_left.assign(sg.nghost, 0.0);
    recv_right.assign(sg.nghost, 0.0);
}

void EMFields::zero_currents()
{
    std::fill(rho.begin(), rho.end(), 0.0);
}

void EMFields::accumulate_moments(const Species& sp)
{
    const int ng = sp.sgrid->nghost;
    const int nxl = sp.sgrid->nx_local;
    for (int ix = 0; ix < nxl; ++ix) {
        rho[ix + ng] += sp.charge_density[ix];
    }
}

void EMFields::set_charge_density(const Species& electrons,
                                  const std::vector<double>& beam_density,
                                  const std::vector<double>& ion_density_profile)
{
    const int ng = electrons.sgrid->nghost;
    const int nxl = electrons.sgrid->nx_local;
    std::fill(rho.begin(), rho.end(), 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        const double zni = (ix < static_cast<int>(ion_density_profile.size()))
                         ? ion_density_profile[static_cast<size_t>(ix)]
                         : 0.0;
        const double nb = (ix < static_cast<int>(beam_density.size()))
                        ? beam_density[static_cast<size_t>(ix)]
                        : 0.0;
        rho[ix + ng] = Const::qe * (zni - electrons.number_density[ix] - nb);
    }
}

void EMFields::solve_poisson(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    gather_rho(*this, mpi_size);

    if (mpi_rank == 0) {
        remove_global_mean_charge_source(*this);
        compute_dirichlet_poisson(*this, true, false);
    }

    MPI_Scatterv(global_ex.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 local_rhs.data(), nxl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int ix = 0; ix < nxl; ++ix) Ex[ng + ix] = local_rhs[static_cast<size_t>(ix)];
    exchange_ex_ghosts(mpi_rank, mpi_size);
}

void EMFields::compute_potential(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    gather_rho(*this, mpi_size);

    if (mpi_rank == 0) {
        remove_global_mean_charge_source(*this);
        compute_dirichlet_poisson(*this, false, true);
    }

    MPI_Scatterv(global_phi.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 local_rhs.data(), nxl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int ix = 0; ix < nxl; ++ix) phi[ng + ix] = local_rhs[static_cast<size_t>(ix)];
    exchange_phi_ghosts(mpi_rank, mpi_size);
}

void EMFields::exchange_ex_ghosts(int mpi_rank, int mpi_size)
{
    exchange_scalar_ghosts(*this, Ex, 201, mpi_rank, mpi_size);
}

void EMFields::exchange_phi_ghosts(int mpi_rank, int mpi_size)
{
    exchange_scalar_ghosts(*this, phi, 211, mpi_rank, mpi_size);
}

double EMFields::total_energy() const
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    double energy = 0.0;
    for (int i = ng; i < ng + nxl; ++i) {
        energy += 0.5 * Const::eps0 * Ex[i] * Ex[i];
    }
    return energy * dx;
}
