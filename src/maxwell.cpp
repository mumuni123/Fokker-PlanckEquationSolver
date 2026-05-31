#include "maxwell.h"
#include "fft_utils.h"
#include "species.h"
#include <algorithm>
#include <cmath>
#include <complex>
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

void enforce_discrete_neutrality(EMFields& fields)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;

    double local_sum = 0.0;
    for (int ix = 0; ix < nxl; ++ix) local_sum += fields.rho[ng + ix];

    double global_sum = 0.0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    const double rho_mean = global_sum / static_cast<double>(Param::nx);

    for (int ix = 0; ix < nxl; ++ix) fields.rho[ng + ix] -= rho_mean;
}

void gather_neutral_rho(EMFields& fields, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;
    prepare_mpi_layout(fields, mpi_size);
    enforce_discrete_neutrality(fields);

    for (int ix = 0; ix < nxl; ++ix) {
        fields.local_rhs[static_cast<size_t>(ix)] = fields.rho[ng + ix];
    }

    MPI_Gatherv(fields.local_rhs.data(), nxl, MPI_DOUBLE,
                fields.global_rhs.data(), fields.counts.data(),
                fields.displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void compute_fft_periodic_poisson(EMFields& fields,
                                  bool compute_ex,
                                  bool compute_phi)
{
    const int n = Param::nx;
    const double length = static_cast<double>(n) * fields.dx;
    std::vector<std::complex<double> > rho_hat(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        rho_hat[static_cast<size_t>(i)] =
            std::complex<double>(fields.global_rhs[static_cast<size_t>(i)], 0.0);
    }

    fft_any(rho_hat, false);
    rho_hat[0] = std::complex<double>(0.0, 0.0);

    std::vector<std::complex<double> > ex_hat;
    std::vector<std::complex<double> > phi_hat;
    if (compute_ex) {
        ex_hat.assign(static_cast<size_t>(n), std::complex<double>(0.0, 0.0));
    }
    if (compute_phi) {
        phi_hat.assign(static_cast<size_t>(n), std::complex<double>(0.0, 0.0));
    }

    for (int k = 1; k < n; ++k) {
        const bool nyquist = ((n % 2) == 0 && k == n / 2);
        const int mode = (k <= n / 2) ? k : k - n;
        const double wave_number =
            2.0 * Const::pi * static_cast<double>(mode) / length;
        const double theta = wave_number * fields.dx;
        const double laplace_symbol =
            4.0 * std::sin(0.5 * theta) * std::sin(0.5 * theta)
            / (fields.dx * fields.dx);
        if (!(laplace_symbol > 0.0)) continue;
        const std::complex<double> phi_k =
            rho_hat[static_cast<size_t>(k)] /
            (Const::eps0 * laplace_symbol);
        if (compute_phi) phi_hat[static_cast<size_t>(k)] = phi_k;
        if (compute_ex && !nyquist) {
            const double grad_symbol = std::sin(theta) / fields.dx;
            ex_hat[static_cast<size_t>(k)] =
                std::complex<double>(0.0, -grad_symbol) * phi_k;
        }
    }

    if (compute_ex) {
        fft_any(ex_hat, true);
        double mean_ex = 0.0;
        for (int i = 0; i < n; ++i) {
            const double value = ex_hat[static_cast<size_t>(i)].real();
            fields.global_ex[static_cast<size_t>(i)] = value;
            mean_ex += value;
        }
        mean_ex /= static_cast<double>(n);
        for (int i = 0; i < n; ++i) fields.global_ex[static_cast<size_t>(i)] -= mean_ex;
    }

    if (compute_phi) {
        fft_any(phi_hat, true);
        double mean_phi = 0.0;
        for (int i = 0; i < n; ++i) {
            const double value = phi_hat[static_cast<size_t>(i)].real();
            fields.global_phi[static_cast<size_t>(i)] = value;
            mean_phi += value;
        }
        mean_phi /= static_cast<double>(n);
        for (int i = 0; i < n; ++i) fields.global_phi[static_cast<size_t>(i)] -= mean_phi;
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

    const int left = (mpi_rank - 1 + mpi_size) % mpi_size;
    const int right = (mpi_rank + 1) % mpi_size;

    if (mpi_size == 1) {
        for (int g = 0; g < ng; ++g) {
            a[g] = a[ng + nxl - ng + g];
            a[ng + nxl + g] = a[ng + g];
        }
        return;
    }

    for (int g = 0; g < ng; ++g) {
        fields.send_left[g] = a[ng + g];
        fields.send_right[g] = a[ng + nxl - ng + g];
    }

    MPI_Request reqs[4];
    int nreq = 0;
    MPI_Isend(fields.send_left.data(), ng, MPI_DOUBLE, left, tag_base,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(fields.recv_left.data(), ng, MPI_DOUBLE, left, tag_base + 1,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Isend(fields.send_right.data(), ng, MPI_DOUBLE, right, tag_base + 1,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(fields.recv_right.data(), ng, MPI_DOUBLE, right, tag_base,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    for (int g = 0; g < ng; ++g) {
        a[g] = fields.recv_left[g];
        a[ng + nxl + g] = fields.recv_right[g];
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
    gather_neutral_rho(*this, mpi_size);

    if (mpi_rank == 0) compute_fft_periodic_poisson(*this, true, false);

    MPI_Scatterv(global_ex.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 local_rhs.data(), nxl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int ix = 0; ix < nxl; ++ix) Ex[ng + ix] = local_rhs[static_cast<size_t>(ix)];
    exchange_ex_ghosts(mpi_rank, mpi_size);
}

void EMFields::compute_potential(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    gather_neutral_rho(*this, mpi_size);

    if (mpi_rank == 0) compute_fft_periodic_poisson(*this, false, true);

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
