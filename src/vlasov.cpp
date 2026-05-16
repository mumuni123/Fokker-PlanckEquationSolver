#include "vlasov.h"
#include "maxwell.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <omp.h>

namespace {
double maxwellian_value(const Species& sp, double v)
{
    double norm = sp.density0 * std::pow(2.0 * Const::pi * sp.temperature / sp.mass, -1.5);
    double inv2vth2 = sp.mass / (2.0 * sp.temperature);
    return norm * std::exp(-v * v * inv2vth2);
}

void fill_background_left_boundary(Species& sp, const SpatialGrid& sg)
{
    // Open perturbations: outgoing vx < 0 uses zero-gradient extrapolation,
    // incoming vx > 0 is replenished by the unperturbed Maxwellian reservoir.
    int ng = sg.nghost;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng - 1 - g;
        int ix_src = ng;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double v = sp.vgrid.v(iv);
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = sp.vgrid.mu(imu);
                sp.f[idx3(ix_dst, iv, imu)] = (v * mu > 0.0)
                    ? maxwellian_value(sp, v)
                    : sp.f[idx3(ix_src, iv, imu)];
            }
        }
    }
}

void fill_background_right_boundary(Species& sp, const SpatialGrid& sg)
{
    // Open perturbations: outgoing vx > 0 uses zero-gradient extrapolation,
    // incoming vx < 0 is replenished by the unperturbed Maxwellian reservoir.
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng + nxl + g;
        int ix_src = ng + nxl - 1;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double v = sp.vgrid.v(iv);
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = sp.vgrid.mu(imu);
                sp.f[idx3(ix_dst, iv, imu)] = (v * mu < 0.0)
                    ? maxwellian_value(sp, v)
                    : sp.f[idx3(ix_src, iv, imu)];
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
        fill_background_left_boundary(sp, sg);
    }
}

void fill_global_right_boundary(Species& sp, const SpatialGrid& sg)
{
    if (sp.type == SpeciesType::BEAM) {
        fill_open_right_boundary(sp, sg);
    } else {
        fill_background_right_boundary(sp, sg);
    }
}

int substeps_from_cfl(double cfl)
{
    if (cfl <= Param::velocity_space_cfl) return 1;
    return static_cast<int>(std::ceil(cfl / Param::velocity_space_cfl));
}
}

VlasovSolver::VlasovSolver()
    : last_cfl_v_(0.0),
      last_cfl_mu_(0.0),
      last_nsub_v_(1),
      last_nsub_mu_(1)
{}

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
    double max_cfl = 0.0;
    double max_abs_mu = 0.0;
    for (int imu = 0; imu < Param::Nmu; ++imu) {
        max_abs_mu = std::max(max_abs_mu, std::fabs(sp.vgrid.mu(imu)));
    }

    #pragma omp parallel for schedule(static) reduction(max:max_cfl)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double a_abs = std::fabs(sp.charge * fields.Ex[ix_g] / sp.mass);
        max_cfl = std::max(max_cfl, a_abs * max_abs_mu * dt / sp.vgrid.dv);
    }

    last_nsub_v_ = substeps_from_cfl(max_cfl);
    last_cfl_v_ = max_cfl / last_nsub_v_;
    double dt_sub = dt / last_nsub_v_;

    for (int isub = 0; isub < last_nsub_v_; ++isub) {
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            int ix_g = ix + ng;
            double a = sp.charge * fields.Ex[ix_g] / sp.mass;

            for (int iv = 0; iv < Param::Nv; ++iv) {
                double v_center = std::max(sp.vgrid.v(iv), Param::v_floor);
                double jac = v_center * v_center;

                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    double mu = sp.vgrid.mu(imu);
                    double flux_left = 0.0;
                    double flux_right = 0.0;

                    for (int side = 0; side < 2; ++side) {
                        int iface = iv + side;
                        double v_face = sp.vgrid.v_face(iface);
                        double av = a * mu;
                        double f_up = 0.0;

                        if (av >= 0.0) {
                            if (iface > 0) {
                                f_up = sp.f[idx3(ix_g, iface - 1, imu)];
                            }
                        } else {
                            if (iface < Param::Nv) {
                                f_up = sp.f[idx3(ix_g, iface, imu)];
                            }
                        }

                        double flux = v_face * v_face * av * f_up;
                        if (side == 0) {
                            flux_left = flux;
                        } else {
                            flux_right = flux;
                        }
                    }

                    double g0 = jac * sp.f[idx3(ix_g, iv, imu)];
                    double gval = g0 - dt_sub / sp.vgrid.dv * (flux_right - flux_left);
                    sp.f_tmp[idx3(ix_g, iv, imu)] = std::max(0.0, gval / jac);
                }
            }
        }

        sp.f.swap(sp.f_tmp);
    }
}

void VlasovSolver::advect_mu(Species& sp, const SpatialGrid& sg,
                             const EMFields& fields, double dt)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    double max_cfl = 0.0;
    double max_mu_face_factor = 0.0;
    for (int imu = 0; imu <= Param::Nmu; ++imu) {
        double mu_face = sp.vgrid.mu_face(imu);
        max_mu_face_factor = std::max(max_mu_face_factor,
                                      1.0 - mu_face * mu_face);
    }

    #pragma omp parallel for schedule(static) reduction(max:max_cfl)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double a_abs = std::fabs(sp.charge * fields.Ex[ix_g] / sp.mass);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = std::max(sp.vgrid.v(iv), Param::v_floor);
            max_cfl = std::max(max_cfl,
                               a_abs * max_mu_face_factor * dt
                             / (vv * sp.vgrid.dmu));
        }
    }

    last_nsub_mu_ = substeps_from_cfl(max_cfl);
    last_cfl_mu_ = max_cfl / last_nsub_mu_;
    double dt_sub = dt / last_nsub_mu_;

    for (int isub = 0; isub < last_nsub_mu_; ++isub) {
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            int ix_g = ix + ng;
            double a = sp.charge * fields.Ex[ix_g] / sp.mass;

            for (int iv = 0; iv < Param::Nv; ++iv) {
                double vv = std::max(sp.vgrid.v(iv), Param::v_floor);
                double jac = vv * vv;

                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    double flux_left = 0.0;
                    double flux_right = 0.0;

                    for (int side = 0; side < 2; ++side) {
                        int iface = imu + side;
                        double mu_face = sp.vgrid.mu_face(iface);
                        double amu = a * (1.0 - mu_face * mu_face) / vv;
                        double f_up = 0.0;

                        if (amu >= 0.0) {
                            if (iface > 0) {
                                f_up = sp.f[idx3(ix_g, iv, iface - 1)];
                            }
                        } else {
                            if (iface < Param::Nmu) {
                                f_up = sp.f[idx3(ix_g, iv, iface)];
                            }
                        }

                        double flux = jac * amu * f_up;
                        if (side == 0) {
                            flux_left = flux;
                        } else {
                            flux_right = flux;
                        }
                    }

                    double g0 = jac * sp.f[idx3(ix_g, iv, imu)];
                    double gval = g0 - dt_sub / sp.vgrid.dmu * (flux_right - flux_left);
                    sp.f_tmp[idx3(ix_g, iv, imu)] = std::max(0.0, gval / jac);
                }
            }
        }

        sp.f.swap(sp.f_tmp);
    }
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
