#include "vlasov.h"
#include "maxwell.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <omp.h>

namespace {
void fill_background_left_boundary(Species& sp, const SpatialGrid& sg)
{
    // Open perturbations: outgoing vx < 0 uses zero-gradient extrapolation,
    // incoming vx > 0 is replenished by the return-current Maxwellian reservoir.
    int ng = sg.nghost;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng - 1 - g;
        int ix_src = ng;
        const size_t dst_base = static_cast<size_t>(ix_dst) * Param::Nvmu;
        const size_t src_base = static_cast<size_t>(ix_src) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const size_t row = static_cast<size_t>(iv) * Param::Nmu;
            const size_t dst_row = dst_base + row;
            const size_t src_row = src_base + row;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k = row + imu;
                const double f_reservoir = sp.maxwellian_left_boundary[k];
                sp.f[dst_row + imu] = (sp.vgrid.vx_cells[k] > 0.0)
                    ? f_reservoir
                    : sp.f[src_row + imu];
            }
        }
    }
}

void fill_background_right_boundary(Species& sp, const SpatialGrid& sg)
{
    // Open perturbations: outgoing vx > 0 uses zero-gradient extrapolation,
    // incoming vx < 0 is replenished by the return-current Maxwellian reservoir.
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng + nxl + g;
        int ix_src = ng + nxl - 1;
        const size_t dst_base = static_cast<size_t>(ix_dst) * Param::Nvmu;
        const size_t src_base = static_cast<size_t>(ix_src) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const size_t row = static_cast<size_t>(iv) * Param::Nmu;
            const size_t dst_row = dst_base + row;
            const size_t src_row = src_base + row;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k = row + imu;
                const double f_reservoir = sp.maxwellian_right_boundary[k];
                sp.f[dst_row + imu] = (sp.vgrid.vx_cells[k] < 0.0)
                    ? f_reservoir
                    : sp.f[src_row + imu];
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

int substeps_from_cfl(double cfl, double cfl_limit)
{
    if (cfl <= cfl_limit) return 1;
    return static_cast<int>(std::ceil(cfl / cfl_limit));
}

RemapStencil build_translation_stencil(double shift_over_dx)
{
    RemapStencil s;
    s.first = 0;
    s.count = 0;
    for (int i = 0; i < 6; ++i) s.weight[i] = 0.0;

    const double left = -shift_over_dx;
    const int first = static_cast<int>(std::floor(left));
    const double frac = left - first;
    const double eps = 1.0e-14;

    if (frac < eps) {
        s.first = first;
        s.count = 1;
        s.weight[0] = 1.0;
    } else if (1.0 - frac < eps) {
        s.first = first + 1;
        s.count = 1;
        s.weight[0] = 1.0;
    } else {
        s.first = first;
        s.count = 2;
        s.weight[0] = 1.0 - frac;
        s.weight[1] = frac;
    }
    return s;
}

}

VlasovSolver::VlasovSolver()
    : last_cfl_v_(0.0),
      last_cfl_mu_(0.0),
      last_loss_v_(0.0),
      last_loss_v_low_(0.0),
      last_loss_v_high_(0.0),
      last_loss_mu_(0.0),
      last_nsub_v_(1),
      last_nsub_mu_(1),
      step_diagnostics_enabled_(false)
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
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const double max_cfl = sp.vgrid.v_max * sp.vgrid.max_abs_mu * dt / sg.dx;
    const double cfl_limit = std::min(Param::semi_lagrangian_cfl,
                                      std::max(1.0, (double)sg.nghost - 0.5));
    const int nsub_x = substeps_from_cfl(max_cfl, cfl_limit);
    const double dt_sub = dt / nsub_x;

    x_stencil_.resize(Param::Nvmu);
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < static_cast<int>(Param::Nvmu); ++k) {
        x_stencil_[k] = build_translation_stencil(sp.vgrid.vx_cells[k] * dt_sub / sg.dx);
    }

    for (int isub = 0; isub < nsub_x; ++isub) {
        exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int ix = ng; ix < ng + nxl; ++ix) {
            for (int iv = 0; iv < Param::Nv; ++iv) {
                const size_t xbase = static_cast<size_t>(ix) * Param::Nvmu;
                const size_t row = static_cast<size_t>(iv) * Param::Nmu;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    // The Vlasov equation keeps vx explicitly; in spherical variables
                    // vx = |v| cos(theta) = v * mu.
                    const size_t k = row + imu;
                    const size_t offset = xbase + k;
                    const RemapStencil& st = x_stencil_[k];
                    double fval = 0.0;
                    for (int m = 0; m < st.count; ++m) {
                        const int donor_ix = ix + st.first + m;
                        if (donor_ix >= 0 && donor_ix < sg.nx_total) {
                            fval += st.weight[m] *
                                    sp.f[static_cast<size_t>(donor_ix) * Param::Nvmu + k];
                        }
                    }
                    sp.f_tmp[offset] = std::max(0.0, fval);
                }
            }
        }

        sp.f.swap(sp.f_tmp);
    }
}

void VlasovSolver::advect_v(Species& sp, const SpatialGrid& sg,
                            const EMFields& fields, double dt)
{
    const bool track_loss = step_diagnostics_enabled_ || Param::abort_on_vmax_loss;
    last_loss_v_ = 0.0;
    last_loss_v_low_ = 0.0;
    last_loss_v_high_ = 0.0;

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    double max_cfl = 0.0;
    const double cfl_scale = dt * sp.vgrid.inv_dv;

    #pragma omp parallel for schedule(static) reduction(max:max_cfl)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const double a_abs = std::fabs(sp.charge * fields.Ex[ix_g] / sp.mass);
        max_cfl = std::max(max_cfl, a_abs * cfl_scale);
    }

    last_nsub_v_ = substeps_from_cfl(max_cfl, Param::velocity_space_cfl);
    last_cfl_v_ = max_cfl / last_nsub_v_;
    const double dt_sub = dt / last_nsub_v_;

    for (int isub = 0; isub < last_nsub_v_; ++isub) {
        double n_before_sub = 0.0;
        double n_after_sub = 0.0;
        double loss_low_sub = 0.0;
        double loss_high_sub = 0.0;

        #pragma omp parallel for schedule(static) \
            reduction(+:n_before_sub,n_after_sub,loss_low_sub,loss_high_sub)
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
            const double dvx = (sp.charge * fields.Ex[ix_g] / sp.mass) * dt_sub;

            std::fill(sp.f_tmp.begin() + xbase,
                      sp.f_tmp.begin() + xbase + Param::Nvmu,
                      0.0);

            for (int iv = 0; iv < Param::Nv; ++iv) {
                const double v = sp.vgrid.v_cells[iv];
                const double v2 = v * v;
                const double shell = sp.vgrid.moment_weight[iv];
                const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const size_t offset = row + imu;
                    const double cell_mass = std::max(0.0, sp.f[offset]) * shell;
                    if (cell_mass == 0.0) continue;
                    if (track_loss) n_before_sub += cell_mass * sg.dx;

                    const double mu = sp.vgrid.mu_cells[imu];
                    const double vx_new = v * mu + dvx;
                    const double vperp2 = std::max(0.0, v2 * (1.0 - mu * mu));
                    const double v_new = std::sqrt(vx_new * vx_new + vperp2);

                    if (v_new >= sp.vgrid.v_max) {
                        if (track_loss) loss_high_sub += cell_mass * sg.dx;
                        continue;
                    }
                    if (v_new < sp.vgrid.v_min) {
                        if (track_loss) loss_low_sub += cell_mass * sg.dx;
                        continue;
                    }

                    double mu_new = 0.0;
                    if (v_new > Param::v_floor) {
                        mu_new = std::max(-1.0, std::min(1.0, vx_new / v_new));
                    }

                    double sv = (v_new - sp.vgrid.v_min) * sp.vgrid.inv_dv - 0.5;
                    int iv0 = static_cast<int>(std::floor(sv));
                    double wv1 = sv - iv0;
                    if (iv0 < 0) {
                        iv0 = 0;
                        wv1 = 0.0;
                    } else if (iv0 >= Param::Nv - 1) {
                        iv0 = Param::Nv - 1;
                        wv1 = 0.0;
                    }
                    const int iv1 = (iv0 + 1 < Param::Nv) ? iv0 + 1 : iv0;
                    const double wv0 = 1.0 - wv1;

                    double smu = (mu_new - sp.vgrid.mu_min) * sp.vgrid.inv_dmu - 0.5;
                    int imu0 = static_cast<int>(std::floor(smu));
                    double wm1 = smu - imu0;
                    if (imu0 < 0) {
                        imu0 = 0;
                        wm1 = 0.0;
                    } else if (imu0 >= Param::Nmu - 1) {
                        imu0 = Param::Nmu - 1;
                        wm1 = 0.0;
                    }
                    const int imu1 = (imu0 + 1 < Param::Nmu) ? imu0 + 1 : imu0;
                    const double wm0 = 1.0 - wm1;

                    sp.f_tmp[xbase + static_cast<size_t>(iv0) * Param::Nmu + imu0]
                        += cell_mass * wv0 * wm0;
                    if (imu1 != imu0) {
                        sp.f_tmp[xbase + static_cast<size_t>(iv0) * Param::Nmu + imu1]
                            += cell_mass * wv0 * wm1;
                    }
                    if (iv1 != iv0) {
                        sp.f_tmp[xbase + static_cast<size_t>(iv1) * Param::Nmu + imu0]
                            += cell_mass * wv1 * wm0;
                        if (imu1 != imu0) {
                            sp.f_tmp[xbase + static_cast<size_t>(iv1) * Param::Nmu + imu1]
                                += cell_mass * wv1 * wm1;
                        }
                    }
                }
            }

            for (int iv = 0; iv < Param::Nv; ++iv) {
                const double inv_shell = 1.0 / sp.vgrid.moment_weight[iv];
                const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const size_t offset = row + imu;
                    const double cell_mass = sp.f_tmp[offset];
                    if (track_loss) n_after_sub += cell_mass * sg.dx;
                    sp.f_tmp[offset] = std::max(0.0, cell_mass * inv_shell);
                }
            }
        }

        sp.f.swap(sp.f_tmp);
        if (track_loss) {
            last_loss_v_ += n_before_sub - n_after_sub;
            last_loss_v_low_ += loss_low_sub;
            last_loss_v_high_ += loss_high_sub;
        }
    }
}

void VlasovSolver::advect_mu(Species& sp, const SpatialGrid& sg,
                             const EMFields& fields, double dt)
{
    (void)sp;
    (void)sg;
    (void)fields;
    (void)dt;
    last_cfl_mu_ = 0.0;
    last_nsub_mu_ = 0;
    last_loss_mu_ = 0.0;
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
