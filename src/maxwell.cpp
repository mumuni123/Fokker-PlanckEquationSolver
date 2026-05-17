#include "maxwell.h"
#include "species.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <vector>

namespace {
void dirichlet_coefficients(int i, int n, double& a, double& b, double& c)
{
    a = (i == 0) ? 0.0 : -1.0;
    b = (i == 0 || i == n - 1) ? 3.0 : 2.0;
    c = (i == n - 1) ? 0.0 : -1.0;
}

void solve_tridiagonal(std::vector<double> a,
                       std::vector<double> b,
                       std::vector<double> c,
                       std::vector<double> d,
                       std::vector<double>& x)
{
    const int n = static_cast<int>(d.size());
    x.assign(n, 0.0);
    if (n == 0) return;

    for (int i = 1; i < n; ++i) {
        const double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

void solve_tridiagonal_cached(const std::vector<double>& a,
                              const std::vector<double>& b,
                              const std::vector<double>& c,
                              const std::vector<double>& d,
                              std::vector<double>& x,
                              std::vector<double>& aw,
                              std::vector<double>& bw,
                              std::vector<double>& cw,
                              std::vector<double>& dw)
{
    aw = a;
    bw = b;
    cw = c;
    dw = d;
    const int n = static_cast<int>(dw.size());
    x.assign(n, 0.0);
    if (n == 0) return;

    for (int i = 1; i < n; ++i) {
        const double m = aw[i] / bw[i - 1];
        bw[i] -= m * cw[i - 1];
        dw[i] -= m * dw[i - 1];
    }

    x[n - 1] = dw[n - 1] / bw[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (dw[i] - cw[i] * x[i + 1]) / bw[i];
    }
}

void solve_dense(std::vector<double> a,
                 std::vector<double> b,
                 int n,
                 std::vector<double>& x)
{
    x.assign(n, 0.0);
    if (n == 0) return;

    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double pivot_abs = std::fabs(a[static_cast<size_t>(k) * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double v = std::fabs(a[static_cast<size_t>(i) * n + k]);
            if (v > pivot_abs) {
                pivot_abs = v;
                pivot = i;
            }
        }
        if (pivot != k) {
            for (int j = k; j < n; ++j) {
                std::swap(a[static_cast<size_t>(k) * n + j],
                          a[static_cast<size_t>(pivot) * n + j]);
            }
            std::swap(b[k], b[pivot]);
        }

        const double diag = a[static_cast<size_t>(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double m = a[static_cast<size_t>(i) * n + k] / diag;
            if (m == 0.0) continue;
            a[static_cast<size_t>(i) * n + k] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                a[static_cast<size_t>(i) * n + j] -=
                    m * a[static_cast<size_t>(k) * n + j];
            }
            b[i] -= m * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double rhs = b[i];
        for (int j = i + 1; j < n; ++j) {
            rhs -= a[static_cast<size_t>(i) * n + j] * x[j];
        }
        x[i] = rhs / a[static_cast<size_t>(i) * n + i];
    }
}

void build_global_dirichlet_system(const std::vector<double>& rhs,
                                   std::vector<double>& a,
                                   std::vector<double>& b,
                                   std::vector<double>& c,
                                   std::vector<double>& d)
{
    const int n = static_cast<int>(rhs.size());
    a.resize(n);
    b.resize(n);
    c.resize(n);
    d = rhs;
    for (int i = 0; i < n; ++i) {
        dirichlet_coefficients(i, n, a[i], b[i], c[i]);
    }
}

void gather_global_rhs(EMFields& fields,
                       const std::vector<int>& counts,
                       const std::vector<int>& displs,
                       int nxl,
                       std::vector<double>& global_rhs)
{
    const int ng = Param::Nghost;
    fields.local_rhs.resize(nxl);
    const double scale = fields.dx * fields.dx / Const::eps0;
    for (int i = 0; i < nxl; ++i) {
        fields.local_rhs[i] = scale * fields.rho[ng + i];
    }

    global_rhs.assign(Param::nx, 0.0);
    MPI_Allgatherv(fields.local_rhs.data(), nxl, MPI_DOUBLE,
                   global_rhs.data(), counts.data(), displs.data(),
                   MPI_DOUBLE, MPI_COMM_WORLD);
}

void copy_global_phi_to_local(EMFields& fields,
                              const std::vector<int>& displs,
                              int mpi_rank,
                              int nxl,
                              const std::vector<double>& global_phi)
{
    const int ng = Param::Nghost;
    const int start = displs[mpi_rank];
    for (int i = 0; i < nxl; ++i) {
        fields.phi[ng + i] = global_phi[start + i];
    }
}

void solve_parallel_cyclic_reduction_system(const std::vector<double>& rhs,
                                            std::vector<double>& phi)
{
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
    build_global_dirichlet_system(rhs, a, b, c, d);

    const int n = static_cast<int>(rhs.size());
    std::vector<double> a_old(n);
    std::vector<double> b_old(n);
    std::vector<double> c_old(n);
    std::vector<double> d_old(n);

    for (int stride = 1; stride < n; stride <<= 1) {
        a_old = a;
        b_old = b;
        c_old = c;
        d_old = d;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            const int il = i - stride;
            const int ir = i + stride;
            const double alpha = (il >= 0) ? a_old[i] / b_old[il] : 0.0;
            const double gamma = (ir < n) ? c_old[i] / b_old[ir] : 0.0;

            b[i] = b_old[i];
            d[i] = d_old[i];
            a[i] = 0.0;
            c[i] = 0.0;

            if (il >= 0) {
                b[i] -= alpha * c_old[il];
                d[i] -= alpha * d_old[il];
                a[i] = -alpha * a_old[il];
            }
            if (ir < n) {
                b[i] -= gamma * a_old[ir];
                d[i] -= gamma * d_old[ir];
                c[i] = -gamma * c_old[ir];
            }
        }
    }

    phi.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        phi[i] = d[i] / b[i];
    }
}

void apply_dirichlet_operator(const std::vector<double>& x,
                              std::vector<double>& ax)
{
    const int n = static_cast<int>(x.size());
    ax.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;
        dirichlet_coefficients(i, n, a, b, c);
        ax[i] = b * x[i];
        if (i > 0) ax[i] += a * x[i - 1];
        if (i + 1 < n) ax[i] += c * x[i + 1];
    }
}

void weighted_jacobi(const std::vector<double>& rhs,
                     std::vector<double>& x,
                     int sweeps)
{
    const int n = static_cast<int>(rhs.size());
    std::vector<double> next(n, 0.0);
    const double omega = 2.0 / 3.0;

    for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (int i = 0; i < n; ++i) {
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            dirichlet_coefficients(i, n, a, b, c);
            double r = rhs[i];
            if (i > 0) r -= a * x[i - 1];
            if (i + 1 < n) r -= c * x[i + 1];
            const double candidate = r / b;
            next[i] = (1.0 - omega) * x[i] + omega * candidate;
        }
        x.swap(next);
    }
}

void restrict_residual(const std::vector<double>& residual,
                       std::vector<double>& coarse)
{
    const int nf = static_cast<int>(residual.size());
    const int nc = std::max(1, nf / 2);
    coarse.assign(nc, 0.0);
    for (int ic = 0; ic < nc; ++ic) {
        const int i = std::min(2 * ic, nf - 1);
        double value = 0.5 * residual[i];
        if (i > 0) value += 0.25 * residual[i - 1];
        if (i + 1 < nf) value += 0.25 * residual[i + 1];
        coarse[ic] = value;
    }
}

void prolongate_and_add(const std::vector<double>& coarse,
                        std::vector<double>& fine)
{
    const int nf = static_cast<int>(fine.size());
    const int nc = static_cast<int>(coarse.size());
    for (int i = 0; i < nf; ++i) {
        const int ic = i / 2;
        if ((i & 1) == 0) {
            fine[i] += coarse[std::min(ic, nc - 1)];
        } else {
            const int il = std::min(ic, nc - 1);
            const int ir = std::min(ic + 1, nc - 1);
            fine[i] += 0.5 * (coarse[il] + coarse[ir]);
        }
    }
}

void multigrid_vcycle(const std::vector<double>& rhs,
                      std::vector<double>& x)
{
    const int n = static_cast<int>(rhs.size());
    if (n <= 32) {
        std::vector<double> a;
        std::vector<double> b;
        std::vector<double> c;
        std::vector<double> d;
        build_global_dirichlet_system(rhs, a, b, c, d);
        solve_tridiagonal(a, b, c, d, x);
        return;
    }

    weighted_jacobi(rhs, x, Param::poisson_multigrid_pre_smooth);

    std::vector<double> ax;
    apply_dirichlet_operator(x, ax);
    std::vector<double> residual(n, 0.0);
    for (int i = 0; i < n; ++i) {
        residual[i] = rhs[i] - ax[i];
    }

    std::vector<double> coarse_rhs;
    restrict_residual(residual, coarse_rhs);
    std::vector<double> coarse_error(coarse_rhs.size(), 0.0);
    multigrid_vcycle(coarse_rhs, coarse_error);
    prolongate_and_add(coarse_error, x);

    weighted_jacobi(rhs, x, Param::poisson_multigrid_post_smooth);
}

void solve_multigrid_system(const std::vector<double>& rhs,
                            std::vector<double>& phi)
{
    phi.assign(rhs.size(), 0.0);
    for (int iter = 0; iter < Param::poisson_multigrid_vcycles; ++iter) {
        multigrid_vcycle(rhs, phi);
    }
}

void solve_global_pcr(EMFields& fields,
                      const std::vector<int>& counts,
                      const std::vector<int>& displs,
                      int mpi_rank,
                      int nxl)
{
    gather_global_rhs(fields, counts, displs, nxl, fields.global_rhs);
    solve_parallel_cyclic_reduction_system(fields.global_rhs, fields.global_phi);
    copy_global_phi_to_local(fields, displs, mpi_rank, nxl, fields.global_phi);
}

void solve_global_multigrid(EMFields& fields,
                            const std::vector<int>& counts,
                            const std::vector<int>& displs,
                            int mpi_rank,
                            int nxl)
{
    gather_global_rhs(fields, counts, displs, nxl, fields.global_rhs);
    solve_multigrid_system(fields.global_rhs, fields.global_phi);
    copy_global_phi_to_local(fields, displs, mpi_rank, nxl, fields.global_phi);
}

void solve_distributed_tridiagonal(EMFields& fields,
                                   const std::vector<int>& counts,
                                   const std::vector<int>& displs,
                                   int mpi_rank,
                                   int mpi_size,
                                   int nxl)
{
    int local_ok = (nxl >= 2) ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) {
        solve_global_pcr(fields, counts, displs, mpi_rank, nxl);
        return;
    }

    const int ng = Param::Nghost;
    const int global_start = displs[mpi_rank];
    const double scale = fields.dx * fields.dx / Const::eps0;

    fields.tri_a.resize(nxl);
    fields.tri_b.resize(nxl);
    fields.tri_c.resize(nxl);
    fields.tri_rhs.resize(nxl);
    for (int i = 0; i < nxl; ++i) {
        const int ig = global_start + i;
        dirichlet_coefficients(ig, Param::nx, fields.tri_a[i],
                               fields.tri_b[i], fields.tri_c[i]);
        fields.tri_rhs[i] = scale * fields.rho[ng + i];
    }
    if (mpi_rank > 0) fields.tri_a[0] = 0.0;
    if (mpi_rank + 1 < mpi_size) fields.tri_c[nxl - 1] = 0.0;

    fields.tri_left_response.assign(nxl, 0.0);
    fields.tri_right_response.assign(nxl, 0.0);
    solve_tridiagonal_cached(fields.tri_a, fields.tri_b, fields.tri_c,
                             fields.tri_rhs, fields.tri_y,
                             fields.tri_work_a, fields.tri_work_b,
                             fields.tri_work_c, fields.tri_work_d);

    if (mpi_rank > 0) {
        fields.tri_rhs.assign(nxl, 0.0);
        fields.tri_rhs[0] = 1.0;
        solve_tridiagonal_cached(fields.tri_a, fields.tri_b, fields.tri_c,
                                 fields.tri_rhs, fields.tri_left_response,
                                 fields.tri_work_a, fields.tri_work_b,
                                 fields.tri_work_c, fields.tri_work_d);
    }
    if (mpi_rank + 1 < mpi_size) {
        fields.tri_rhs.assign(nxl, 0.0);
        fields.tri_rhs[nxl - 1] = 1.0;
        solve_tridiagonal_cached(fields.tri_a, fields.tri_b, fields.tri_c,
                                 fields.tri_rhs, fields.tri_right_response,
                                 fields.tri_work_a, fields.tri_work_b,
                                 fields.tri_work_c, fields.tri_work_d);
    }

    double local_interface[6] = {
        fields.tri_y.front(), fields.tri_y.back(),
        fields.tri_left_response.front(), fields.tri_left_response.back(),
        fields.tri_right_response.front(), fields.tri_right_response.back()
    };

    if (mpi_rank == 0) {
        fields.all_interfaces.resize(static_cast<size_t>(6) * mpi_size);
    }
    MPI_Gather(local_interface, 6, MPI_DOUBLE,
               mpi_rank == 0 ? fields.all_interfaces.data() : NULL, 6, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    fields.interface_values.assign(static_cast<size_t>(2) * mpi_size, 0.0);
    if (mpi_rank == 0) {
        const int nred = 2 * mpi_size;
        std::vector<double> mat(static_cast<size_t>(nred) * nred, 0.0);
        std::vector<double> reduced_rhs(nred, 0.0);

        for (int r = 0; r < mpi_size; ++r) {
            const double yf = fields.all_interfaces[static_cast<size_t>(6) * r + 0];
            const double yl = fields.all_interfaces[static_cast<size_t>(6) * r + 1];
            const double lf = fields.all_interfaces[static_cast<size_t>(6) * r + 2];
            const double ll = fields.all_interfaces[static_cast<size_t>(6) * r + 3];
            const double rf = fields.all_interfaces[static_cast<size_t>(6) * r + 4];
            const double rl = fields.all_interfaces[static_cast<size_t>(6) * r + 5];

            const int first = 2 * r;
            const int last = first + 1;
            mat[static_cast<size_t>(first) * nred + first] = 1.0;
            mat[static_cast<size_t>(last) * nred + last] = 1.0;
            reduced_rhs[first] = yf;
            reduced_rhs[last] = yl;

            if (r > 0) {
                const int prev_last = 2 * (r - 1) + 1;
                mat[static_cast<size_t>(first) * nred + prev_last] -= lf;
                mat[static_cast<size_t>(last) * nred + prev_last] -= ll;
            }
            if (r + 1 < mpi_size) {
                const int next_first = 2 * (r + 1);
                mat[static_cast<size_t>(first) * nred + next_first] -= rf;
                mat[static_cast<size_t>(last) * nred + next_first] -= rl;
            }
        }

        solve_dense(mat, reduced_rhs, nred, fields.interface_values);
    }

    MPI_Bcast(fields.interface_values.data(), static_cast<int>(fields.interface_values.size()),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

    const double left_value = (mpi_rank > 0)
                            ? fields.interface_values[static_cast<size_t>(2) * (mpi_rank - 1) + 1]
                            : 0.0;
    const double right_value = (mpi_rank + 1 < mpi_size)
                             ? fields.interface_values[static_cast<size_t>(2) * (mpi_rank + 1)]
                             : 0.0;
    for (int i = 0; i < nxl; ++i) {
        fields.phi[ng + i] = fields.tri_y[i]
                           + fields.tri_left_response[i] * left_value
                           + fields.tri_right_response[i] * right_value;
    }
}

void exchange_phi_edge_ghosts(EMFields& fields, int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;
    const int left = mpi_rank - 1;
    const int right = mpi_rank + 1;

    double send_left = fields.phi[ng];
    double send_right = fields.phi[ng + nxl - 1];
    double recv_left = 0.0;
    double recv_right = 0.0;

    MPI_Request reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(&send_left, 1, MPI_DOUBLE, left, 211, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, left, 212, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right, 1, MPI_DOUBLE, right, 212, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, right, 211, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) {
        MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
    }

    fields.phi[ng - 1] = (left >= 0) ? recv_left : 0.0;
    fields.phi[ng + nxl] = (right < mpi_size) ? recv_right : 0.0;
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
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;

    if (counts_mpi_size != mpi_size || counts_nx_local != nxl) {
        counts.resize(mpi_size);
        displs.resize(mpi_size);
        MPI_Allgather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) {
            displs[r] = displs[r - 1] + counts[r - 1];
        }
        counts_mpi_size = mpi_size;
        counts_nx_local = nxl;
    }

    switch (Param::poisson_solver) {
    case Param::POISSON_PARALLEL_CYCLIC_REDUCTION:
        solve_global_pcr(*this, counts, displs, mpi_rank, nxl);
        break;
    case Param::POISSON_MULTIGRID:
        solve_global_multigrid(*this, counts, displs, mpi_rank, nxl);
        break;
    case Param::POISSON_DISTRIBUTED_TRIDIAGONAL:
    default:
        solve_distributed_tridiagonal(*this, counts, displs, mpi_rank, mpi_size, nxl);
        break;
    }

    exchange_phi_edge_ghosts(*this, mpi_rank, mpi_size);

    const int global_start = displs[mpi_rank];
    for (int ix = 0; ix < nxl; ++ix) {
        const int i = ng + ix;
        const int ig = global_start + ix;
        if (ig == 0) {
            Ex[i] = -2.0 * phi[i] / dx;
        } else if (ig == Param::nx - 1) {
            Ex[i] = 2.0 * phi[i] / dx;
        } else {
            Ex[i] = -(phi[i + 1] - phi[i - 1]) / (2.0 * dx);
        }
    }

    exchange_ghosts(mpi_rank, mpi_size);
}

void EMFields::exchange_ghosts(int mpi_rank, int mpi_size)
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;
    int left = mpi_rank - 1;
    int right = mpi_rank + 1;

    if (send_left.size() != static_cast<size_t>(ng)) {
        send_left.resize(ng);
        send_right.resize(ng);
        recv_left.resize(ng);
        recv_right.resize(ng);
    }

    for (int g = 0; g < ng; ++g) {
        send_left[g] = Ex[ng + g];
        send_right[g] = Ex[ng + nxl - ng + g];
    }

    MPI_Request reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(send_left.data(), ng, MPI_DOUBLE, left, 201, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_left.data(), ng, MPI_DOUBLE, left, 202, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(send_right.data(), ng, MPI_DOUBLE, right, 202, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_right.data(), ng, MPI_DOUBLE, right, 201, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) {
        MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
    }

    for (int g = 0; g < ng; ++g) {
        Ex[g] = (left >= 0) ? recv_left[g] : Ex[ng];
        Ex[ng + nxl + g] = (right < mpi_size) ? recv_right[g] : Ex[ng + nxl - 1];
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
