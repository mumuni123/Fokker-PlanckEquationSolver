#include "vlasov.h"
#include "maxwell.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <omp.h>

namespace {
inline int mirror_mu_index(int imu)
{
    return Param::Nmu - 1 - imu;
}

void fill_closed_left_boundary(Species& sp, const SpatialGrid& sg)
{
    int ng = sg.nghost;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng - 1 - g;
        int ix_src = ng + g;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                sp.f[idx3(ix_dst, iv, imu)] =
                    sp.f[idx3(ix_src, iv, mirror_mu_index(imu))];
            }
        }
    }
}

void fill_closed_right_boundary(Species& sp, const SpatialGrid& sg)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng + nxl + g;
        int ix_src = ng + nxl - 1 - g;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                sp.f[idx3(ix_dst, iv, imu)] =
                    sp.f[idx3(ix_src, iv, mirror_mu_index(imu))];
            }
        }
    }
}

void fill_open_left_boundary(Species& sp, const SpatialGrid& sg)
{
    int ng = sg.nghost;
    std::fill(sp.f.begin(), sp.f.begin() + static_cast<size_t>(ng) * Param::Nvmu, 0.0);
}

void fill_open_right_boundary(Species& sp, const SpatialGrid& sg)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    size_t begin = static_cast<size_t>(ng + nxl) * Param::Nvmu;
    std::fill(sp.f.begin() + begin, sp.f.end(), 0.0);
}

void fill_global_left_boundary(Species& sp, const SpatialGrid& sg)
{
    if (sp.type == SpeciesType::BEAM) {
        fill_open_left_boundary(sp, sg);
    } else {
        fill_closed_left_boundary(sp, sg);
    }
}

void fill_global_right_boundary(Species& sp, const SpatialGrid& sg)
{
    if (sp.type == SpeciesType::BEAM) {
        fill_open_right_boundary(sp, sg);
    } else {
        fill_closed_right_boundary(sp, sg);
    }
}
}

void VlasovSolver::advect(Species& sp, const SpatialGrid& sg,
                          const EMFields& fields, double dt,
                          int mpi_rank, int mpi_size)
{
    advect_x(sp, sg, dt, mpi_rank, mpi_size);
    advect_v(sp, sg, fields, dt);
    advect_mu(sp, sg, fields, dt);
}

void VlasovSolver::advect_x(Species& sp, const SpatialGrid& sg, double dt,
                            int mpi_rank, int mpi_size)
{
    exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);

    int ng = sg.nghost;
    int nxl = sg.nx_local;
    double dxi = 1.0 / sg.dx;

    #pragma omp parallel for schedule(static)
    for (int ix = ng; ix < ng + nxl; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double v = sp.vgrid.v(iv);
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                // The Vlasov equation keeps vx explicitly; in spherical variables
                // vx = |v| cos(theta) = v * mu.
                double vel = v * sp.vgrid.mu(imu);
                double f0 = sp.f[idx3(ix, iv, imu)];
                double up = (vel >= 0.0)
                          ? sp.f[idx3(ix - 1, iv, imu)]
                          : sp.f[idx3(ix + 1, iv, imu)];
                double fval = f0 - dt * vel * dxi * (f0 - up);
                sp.f_tmp[idx3(ix, iv, imu)] = std::max(0.0, fval);
            }
        }
    }

    sp.f.swap(sp.f_tmp);
}

void VlasovSolver::advect_v(Species& sp, const SpatialGrid& sg,
                            const EMFields& fields, double dt)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double a = sp.charge * fields.Ex[ix_g] / sp.mass;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = sp.vgrid.mu(imu);
                // a_x * d/dvx = a_x * [mu d/dv + (1-mu^2)/v d/dmu].
                double av = a * mu;
                double f0 = sp.f[idx3(ix_g, iv, imu)];
                double up = 0.0;
                if (av >= 0.0) {
                    up = (iv > 0) ? sp.f[idx3(ix_g, iv - 1, imu)] : 0.0;
                } else {
                    up = (iv + 1 < Param::Nv) ? sp.f[idx3(ix_g, iv + 1, imu)] : 0.0;
                }
                double fval = f0 - dt * av / sp.vgrid.dv * (f0 - up);
                sp.f_tmp[idx3(ix_g, iv, imu)] = std::max(0.0, fval);
            }
        }
    }

    sp.f.swap(sp.f_tmp);
}

void VlasovSolver::advect_mu(Species& sp, const SpatialGrid& sg,
                             const EMFields& fields, double dt)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double a = sp.charge * fields.Ex[ix_g] / sp.mass;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = std::max(sp.vgrid.v(iv), Param::v_floor);

            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = sp.vgrid.mu(imu);
                double amu = a * (1.0 - mu * mu) / vv;
                double f0 = sp.f[idx3(ix_g, iv, imu)];
                double up = 0.0;
                if (amu >= 0.0) {
                    up = (imu > 0) ? sp.f[idx3(ix_g, iv, imu - 1)] : sp.f[idx3(ix_g, iv, imu)];
                } else {
                    up = (imu + 1 < Param::Nmu) ? sp.f[idx3(ix_g, iv, imu + 1)] : sp.f[idx3(ix_g, iv, imu)];
                }
                double fval = f0 - dt * amu / sp.vgrid.dmu * (f0 - up);
                sp.f_tmp[idx3(ix_g, iv, imu)] = std::max(0.0, fval);
            }
        }
    }

    sp.f.swap(sp.f_tmp);
}

void VlasovSolver::exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                                     int mpi_rank, int mpi_size)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    size_t slice_size = Param::Nvmu;

    int left_rank = mpi_rank - 1;
    int right_rank = mpi_rank + 1;

    size_t buffer_size = static_cast<size_t>(ng) * slice_size;
    if (send_left_.size() != buffer_size) {
        send_left_.resize(buffer_size);
        send_right_.resize(buffer_size);
        recv_left_.resize(buffer_size);
        recv_right_.resize(buffer_size);
    }

    std::memcpy(send_left_.data(),
                &sp.f[static_cast<size_t>(ng) * slice_size],
                buffer_size * sizeof(double));
    std::memcpy(send_right_.data(),
                &sp.f[static_cast<size_t>(ng + nxl - ng) * slice_size],
                buffer_size * sizeof(double));

    MPI_Request reqs[4];
    int nreq = 0;
    if (left_rank >= 0) {
        MPI_Isend(send_left_.data(), (int)buffer_size, MPI_DOUBLE,
                  left_rank, 101, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_left_.data(), (int)buffer_size, MPI_DOUBLE,
                  left_rank, 102, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (right_rank < mpi_size) {
        MPI_Isend(send_right_.data(), (int)buffer_size, MPI_DOUBLE,
                  right_rank, 102, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Irecv(recv_right_.data(), (int)buffer_size, MPI_DOUBLE,
                  right_rank, 101, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (nreq > 0) {
        MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
    }

    if (left_rank >= 0) {
        std::memcpy(&sp.f[0], recv_left_.data(), buffer_size * sizeof(double));
    } else {
        fill_global_left_boundary(sp, sg);
    }

    if (right_rank < mpi_size) {
        std::memcpy(&sp.f[static_cast<size_t>(ng + nxl) * slice_size],
                    recv_right_.data(), buffer_size * sizeof(double));
    } else {
        fill_global_right_boundary(sp, sg);
    }
}
