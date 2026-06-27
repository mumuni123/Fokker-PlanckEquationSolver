#include "vlasov.h"
#include "maxwell.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <mpi.h>
#include <omp.h>

namespace {
void fill_periodic_ghosts_single_rank(Species& sp, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    if (nxl <= 0) return;

    const size_t slice_size = Param::Nvmu;
    for (int g = 0; g < ng; ++g) {
        const int left_src = ng + ((nxl - ng + g) % nxl);
        const int right_src = ng + (g % nxl);
        std::memcpy(&sp.f[static_cast<size_t>(g) * slice_size],
                    &sp.f[static_cast<size_t>(left_src) * slice_size],
                    slice_size * sizeof(double));
        std::memcpy(&sp.f[static_cast<size_t>(ng + nxl + g) * slice_size],
                    &sp.f[static_cast<size_t>(right_src) * slice_size],
                    slice_size * sizeof(double));
    }
}

int substeps_from_cfl(double cfl, double cfl_limit)
{
    if (cfl <= cfl_limit) return 1;
    return static_cast<int>(std::ceil(cfl / cfl_limit));
}

inline double mc_limited_slope(double fm, double f0, double fp)
{
    const double dl = f0 - fm;
    const double dr = fp - f0;
    if (dl * dr <= 0.0) return 0.0;

    const double centered = 0.5 * (dl + dr);
    const double sign = (centered >= 0.0) ? 1.0 : -1.0;
    return sign * std::min(std::fabs(centered),
                           std::min(2.0 * std::fabs(dl),
                                    2.0 * std::fabs(dr)));
}

inline double f_v_line(const Species& sp, size_t xbase, int iv, int imu)
{
    return std::max(0.0, sp.f[xbase + static_cast<size_t>(iv) * Param::Nmu
                              + static_cast<size_t>(imu)]);
}

inline double slope_v_line(const Species& sp, size_t xbase, int iv, int imu)
{
    if (iv <= 0 || iv >= Param::Nv - 1) return 0.0;
    const double fm = f_v_line(sp, xbase, iv - 1, imu);
    const double f0 = f_v_line(sp, xbase, iv, imu);
    const double fp = f_v_line(sp, xbase, iv + 1, imu);
    const double dl =
        (f0 - fm) * sp.vgrid.inv_v_center_dist[iv];
    const double dr =
        (fp - f0) * sp.vgrid.inv_v_center_dist[iv + 1];
    if (dl * dr <= 0.0) return 0.0;

    const double centered =
        (fp - fm) / (sp.vgrid.v_cells[iv + 1] - sp.vgrid.v_cells[iv - 1]);
    const double sign = (centered >= 0.0) ? 1.0 : -1.0;
    return sign * std::min(std::fabs(centered),
                           std::min(2.0 * std::fabs(dl),
                                    2.0 * std::fabs(dr)));
}

inline double reconstruct_v_face(const Species& sp, size_t xbase,
                                 int face, int imu, double speed)
{
    if (speed >= 0.0) {
        const int iv = face - 1;
        return std::max(0.0, f_v_line(sp, xbase, iv, imu)
                             + slope_v_line(sp, xbase, iv, imu)
                             * (sp.vgrid.v_faces[face]
                                - sp.vgrid.v_cells[iv]));
    }

    const int iv = face;
    return std::max(0.0, f_v_line(sp, xbase, iv, imu)
                         + slope_v_line(sp, xbase, iv, imu)
                         * (sp.vgrid.v_faces[face]
                            - sp.vgrid.v_cells[iv]));
}

inline double f_mu_line(const Species& sp, size_t row, int imu)
{
    return std::max(0.0, sp.f[row + static_cast<size_t>(imu)]);
}

inline double slope_mu_line(const Species& sp, size_t row, int imu)
{
    if (imu <= 0 || imu >= Param::Nmu - 1) return 0.0;
    return mc_limited_slope(f_mu_line(sp, row, imu - 1),
                            f_mu_line(sp, row, imu),
                            f_mu_line(sp, row, imu + 1));
}

inline double reconstruct_mu_face(const Species& sp, size_t row,
                                  int face, double speed)
{
    if (speed >= 0.0) {
        const int imu = face - 1;
        return std::max(0.0, f_mu_line(sp, row, imu)
                             + 0.5 * slope_mu_line(sp, row, imu));
    }

    const int imu = face;
    return std::max(0.0, f_mu_line(sp, row, imu)
                         - 0.5 * slope_mu_line(sp, row, imu));
}

inline double staggered_cell_ex(const EMFields& fields, int ix_local, int ix_g)
{
    const size_t left_face = static_cast<size_t>(ix_local);
    const size_t right_face = static_cast<size_t>(ix_local + 1);
    if (right_face < fields.Ex_face.size()) {
        return 0.5 * (fields.Ex_face[left_face] + fields.Ex_face[right_face]);
    }
    return fields.Ex[static_cast<size_t>(ix_g)];
}

void resize_or_zero(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    } else {
        std::fill(values.begin(), values.end(), 0.0);
    }
}

}

VlasovSolver::VlasovSolver()
    : last_cfl_v_(0.0),
      last_cfl_mu_(0.0),
      last_loss_v_(0.0),
      last_loss_v_low_(0.0),
      last_loss_v_high_(0.0),
      last_loss_v_high_local_max_(0.0),
      last_x_at_max_loss_v_high_(0.0),
      last_f_umax_at_max_loss_v_high_(0.0),
      last_integral_f_u_gt_8_at_max_loss_v_high_(0.0),
      last_loss_mu_(0.0),
      last_mass_error_v_(0.0),
      last_mass_error_mu_(0.0),
      last_momentum_delta_v_(0.0),
      last_momentum_delta_mu_(0.0),
      last_energy_delta_v_(0.0),
      last_energy_delta_mu_(0.0),
      last_nsub_v_(1),
      last_nsub_mu_(1),
      step_diagnostics_enabled_(false)
{}

void VlasovSolver::advect(Species& sp, const SpatialGrid& sg,
                          const EMFields& fields, double dt,
                          int mpi_rank, int mpi_size)
{
    advect_x(sp, sg, 0.5 * dt, mpi_rank, mpi_size);
    advect_v(sp, sg, fields, 0.5 * dt, mpi_rank, mpi_size);
    advect_mu(sp, sg, fields, dt);
    advect_v(sp, sg, fields, 0.5 * dt, mpi_rank, mpi_size);
    advect_x(sp, sg, 0.5 * dt, mpi_rank, mpi_size);
}

void VlasovSolver::advect_x(Species& sp, const SpatialGrid& sg, double dt,
                            int mpi_rank, int mpi_size, double time)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    if (dt <= 0.0 || nxl <= 0) return;
    const double max_cfl = sp.vgrid.max_speed * sp.vgrid.max_abs_mu * dt / sg.dx;
    const double cfl_limit = std::min(0.85, Param::semi_lagrangian_cfl);
    const int nsub_x = substeps_from_cfl(max_cfl, cfl_limit);
    const double dt_sub = dt / nsub_x;
    const double dt_dx = dt_sub / sg.dx;
    const double current_average_weight = dt_sub / dt;
    (void)time;

    resize_or_zero(sp.current_face_x, static_cast<size_t>(nxl + 1));
    resize_or_zero(sp.current_x, static_cast<size_t>(nxl));
    const size_t active_faces = static_cast<size_t>(nxl + 1) * Param::Nvmu;
    std::vector<double> f_old(sp.f.size(), 0.0);
    std::vector<double> flux_faces(active_faces, 0.0);
    const int max_midpoint_iters = 24;
    const double midpoint_tol = 1.0e-11;

    for (int isub = 0; isub < nsub_x; ++isub) {
        exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);
        f_old = sp.f;

        bool midpoint_converged = false;
        double global_midpoint_error = 0.0;
        for (int iter = 0; iter < max_midpoint_iters; ++iter) {
            exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);

            #pragma omp parallel for schedule(static)
            for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu);
                 ++k_int) {
                const size_t k = static_cast<size_t>(k_int);
                const double vx = sp.vgrid.vx_cells[k];
                for (int iface = 0; iface < nxl; ++iface) {
                    const int ix_left = ng + iface - 1;
                    const int ix_right = ng + iface;
                    const size_t left_offset =
                        static_cast<size_t>(ix_left) * Param::Nvmu + k;
                    const size_t right_offset =
                        static_cast<size_t>(ix_right) * Param::Nvmu + k;
                    const double f_left_mid =
                        0.5 * (f_old[left_offset] + sp.f[left_offset]);
                    const double f_right_mid =
                        0.5 * (f_old[right_offset] + sp.f[right_offset]);
                    flux_faces[static_cast<size_t>(iface) * Param::Nvmu + k] =
                        0.5 * vx * (f_left_mid + f_right_mid);
                }
            }

            if (mpi_size <= 1) {
                const size_t left_face = 0;
                const size_t right_face = static_cast<size_t>(nxl) * Param::Nvmu;
                std::copy(flux_faces.begin() + left_face,
                          flux_faces.begin() + left_face + Param::Nvmu,
                          flux_faces.begin() + right_face);
            } else {
                const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
                const int right_peer = (mpi_rank + 1) % mpi_size;
                const size_t right_face =
                    static_cast<size_t>(nxl) * Param::Nvmu;
                MPI_Sendrecv(flux_faces.data(), static_cast<int>(Param::Nvmu),
                             MPI_DOUBLE, left_peer, 731,
                             flux_faces.data() + right_face,
                             static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                             right_peer, 731, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
            }

            double local_num = 0.0;
            double local_den = 0.0;
            #pragma omp parallel for collapse(2) schedule(static) \
                reduction(max:local_num,local_den)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu);
                     ++k_int) {
                    const size_t k = static_cast<size_t>(k_int);
                    const size_t offset =
                        static_cast<size_t>(ng + ix) * Param::Nvmu + k;
                    const double updated =
                        f_old[offset]
                      - dt_dx *
                        (flux_faces[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
                       - flux_faces[static_cast<size_t>(ix) * Param::Nvmu + k]);
                    sp.f_tmp[offset] = updated;
                    local_num = std::max(local_num,
                                         std::fabs(updated - sp.f[offset]));
                    local_den = std::max(local_den,
                                         std::max(std::fabs(updated),
                                                  std::fabs(sp.f[offset])));
                }
            }

            const double local_midpoint_error =
                local_num / std::max(1.0, local_den);
            MPI_Allreduce(&local_midpoint_error, &global_midpoint_error, 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            sp.f.swap(sp.f_tmp);
            if (global_midpoint_error < midpoint_tol) {
                midpoint_converged = true;
                break;
            }
        }

        if (!midpoint_converged) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "ERROR: advect_x midpoint iteration failed to "
                             "converge at x-substep %d/%d; residual %.6e. "
                             "Reduce dt or tighten the spatial solve.\n",
                             isub + 1, nsub_x, global_midpoint_error);
            }
            MPI_Abort(MPI_COMM_WORLD, 7);
        }

        #pragma omp parallel
        {
            std::vector<double> current_face_local(
                static_cast<size_t>(nxl + 1), 0.0);
            #pragma omp for schedule(static)
            for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu);
                 ++k_int) {
                const size_t k = static_cast<size_t>(k_int);
                const double vx = sp.vgrid.vx_cells[k];
                if (vx == 0.0) continue;
                const double current_per_f =
                    sp.charge * sp.vgrid.current_weight[k] / vx;
                for (int iface = 0; iface <= nxl; ++iface) {
                    current_face_local[static_cast<size_t>(iface)] +=
                        current_average_weight * current_per_f *
                        flux_faces[static_cast<size_t>(iface) * Param::Nvmu + k];
                }
            }
            #pragma omp critical
            {
                for (int iface = 0; iface <= nxl; ++iface) {
                    sp.current_face_x[static_cast<size_t>(iface)] +=
                        current_face_local[static_cast<size_t>(iface)];
                }
            }
        }
    }

    if (mpi_size <= 1 && sp.current_face_x.size() >= static_cast<size_t>(nxl + 1)) {
        sp.current_face_x[static_cast<size_t>(nxl)] = sp.current_face_x[0];
    } else if (mpi_size > 1 &&
               sp.current_face_x.size() >= static_cast<size_t>(nxl + 1)) {
        const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
        const int right_peer = (mpi_rank + 1) % mpi_size;
        MPI_Sendrecv(sp.current_face_x.data(), 1, MPI_DOUBLE, left_peer, 732,
                     &sp.current_face_x[static_cast<size_t>(nxl)], 1,
                     MPI_DOUBLE, right_peer, 732, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
    }
    for (int ix = 0; ix < nxl; ++ix) {
        sp.current_x[static_cast<size_t>(ix)] =
            0.5 * (sp.current_face_x[static_cast<size_t>(ix)]
                 + sp.current_face_x[static_cast<size_t>(ix + 1)]);
    }
}

void VlasovSolver::advect_v(Species& sp, const SpatialGrid& sg,
                            const EMFields& fields, double dt,
                            int mpi_rank, int mpi_size)
{
    const bool track_step_diagnostics = step_diagnostics_enabled_;
    const bool track_vmax_loss = Param::abort_on_vmax_loss;
    last_loss_v_ = 0.0;
    last_loss_v_low_ = 0.0;
    last_loss_v_high_ = 0.0;
    last_loss_v_high_local_max_ = 0.0;
    last_x_at_max_loss_v_high_ = 0.0;
    last_f_umax_at_max_loss_v_high_ = 0.0;
    last_integral_f_u_gt_8_at_max_loss_v_high_ = 0.0;
    last_mass_error_v_ = 0.0;
    last_momentum_delta_v_ = 0.0;
    last_energy_delta_v_ = 0.0;

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    resize_or_zero(last_energy_current_cell_x_, static_cast<size_t>(nxl));
    resize_or_zero(last_energy_current_face_x_, static_cast<size_t>(nxl + 1));
    double max_cfl = 0.0;
    const bool track_high_loss_location = step_diagnostics_enabled_;
    std::vector<double> loss_high_by_x;
    if (track_high_loss_location) {
        loss_high_by_x.assign(static_cast<size_t>(nxl), 0.0);
    }

    #pragma omp parallel for collapse(2) schedule(static) reduction(max:max_cfl)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const int ix_g = ix + ng;
            const double ex_cell = staggered_cell_ex(fields, ix, ix_g);
            const double udot_abs =
                std::fabs(sp.charge * ex_cell / (sp.mass * Const::c));
            const double cfl =
                dt * udot_abs * sp.vgrid.max_abs_mu
                * sp.vgrid.inv_v_widths[iv];
            max_cfl = std::max(max_cfl, cfl);
        }
    }

    last_nsub_v_ = substeps_from_cfl(max_cfl, Param::velocity_space_cfl);
    last_cfl_v_ = max_cfl / last_nsub_v_;
    const double dt_sub = dt / last_nsub_v_;
    const double face_weight_base = 2.0 * Const::pi * sp.vgrid.dmu;
    std::array<double, Param::Nv> diag_px_base;
    std::array<double, Param::Nv> diag_ke;
    if (track_step_diagnostics) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            diag_px_base[static_cast<size_t>(iv)] =
                sp.mass * Const::c * sp.vgrid.v_cells[iv];
            diag_ke[static_cast<size_t>(iv)] =
                (sp.vgrid.gamma_cells[iv] - 1.0)
                * sp.mass * Const::c * Const::c;
        }
    }

    for (int isub = 0; isub < last_nsub_v_; ++isub) {
        double n_before_sub = 0.0;
        double n_after_sub = 0.0;
        double loss_low_sub = 0.0;
        double loss_high_sub = 0.0;
        double px_before_sub = 0.0;
        double px_after_sub = 0.0;
        double ke_flux_delta_sub = 0.0;

        #pragma omp parallel reduction(+:n_before_sub,n_after_sub,loss_low_sub,loss_high_sub, \
                                         px_before_sub,px_after_sub,ke_flux_delta_sub)
        {
            std::array<double, Param::Nv> mass;
            std::array<double, Param::Nv> new_mass;
            std::array<double, Param::Nv + 1> flux;
            std::array<double, Param::Nv> outgoing;
            std::array<double, Param::Nv> scale;

            #pragma omp for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                double energy_current_dt_sum = 0.0;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                const int ix_g = ix + ng;
                const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
                const double ex_cell = staggered_cell_ex(fields, ix, ix_g);
                const double accel_u =
                    sp.charge * ex_cell / (sp.mass * Const::c);
                const double vdot = accel_u * sp.vgrid.mu_cells[imu];

                flux.fill(0.0);
                outgoing.fill(0.0);
                scale.fill(1.0);

                for (int iv = 0; iv < Param::Nv; ++iv) {
                    const size_t offset = xbase + static_cast<size_t>(iv) * Param::Nmu
                                        + static_cast<size_t>(imu);
                    const double cell_mass =
                        std::max(0.0, sp.f[offset]) * sp.vgrid.moment_weight[iv];
                    mass[static_cast<size_t>(iv)] = cell_mass;
                    new_mass[static_cast<size_t>(iv)] = cell_mass;
                    if (track_step_diagnostics) {
                        const double px =
                            diag_px_base[static_cast<size_t>(iv)]
                            * sp.vgrid.mu_cells[imu];
                        n_before_sub += cell_mass * sg.dx;
                        px_before_sub += cell_mass * px * sg.dx;
                    }
                }

                for (int face = 1; face < Param::Nv; ++face) {
                    const double f_face =
                        reconstruct_v_face(sp, xbase, face, imu, vdot);
                    flux[static_cast<size_t>(face)] =
                        face_weight_base * sp.vgrid.v2_faces[face] * vdot * f_face;
                    const int donor =
                        (flux[static_cast<size_t>(face)] >= 0.0) ? face - 1 : face;
                    outgoing[static_cast<size_t>(donor)] +=
                        dt_sub * std::fabs(flux[static_cast<size_t>(face)]);
                }

                if (vdot > 0.0) {
                    const double f_face =
                        reconstruct_v_face(sp, xbase, Param::Nv, imu, vdot);
                    flux[static_cast<size_t>(Param::Nv)] =
                        face_weight_base * sp.vgrid.v2_faces[Param::Nv] * vdot * f_face;
                    outgoing[static_cast<size_t>(Param::Nv - 1)] +=
                        dt_sub * flux[static_cast<size_t>(Param::Nv)];
                }

                for (int iv = 0; iv < Param::Nv; ++iv) {
                    const double out = outgoing[static_cast<size_t>(iv)];
                    const double m = mass[static_cast<size_t>(iv)];
                    if (out > m && out > 0.0) {
                        scale[static_cast<size_t>(iv)] = m / out;
                    }
                }

                for (int face = 1; face <= Param::Nv; ++face) {
                    double& ff = flux[static_cast<size_t>(face)];
                    if (ff == 0.0) continue;
                    const int donor = (ff >= 0.0) ? face - 1 : face;
                    if (donor >= 0 && donor < Param::Nv) {
                        ff *= scale[static_cast<size_t>(donor)];
                    }
                }

                for (int face = 1; face < Param::Nv; ++face) {
                    const double ff = flux[static_cast<size_t>(face)];
                    if (ff != 0.0) {
                        const double du_center =
                            sp.vgrid.v_cells[face]
                            - sp.vgrid.v_cells[face - 1];
                        const double ke_jump =
                            sp.mass * Const::c
                            * sp.vgrid.chain_speed_faces[face]
                            * du_center;
                        if (track_step_diagnostics) {
                            ke_flux_delta_sub +=
                                dt_sub * ff * ke_jump * sg.dx;
                        }
                        if (ex_cell != 0.0) {
                            energy_current_dt_sum +=
                                dt_sub * ff * ke_jump / ex_cell;
                        }
                    }
                    new_mass[static_cast<size_t>(face - 1)] -= dt_sub * ff;
                    new_mass[static_cast<size_t>(face)] += dt_sub * ff;
                }

                if (flux[static_cast<size_t>(Param::Nv)] > 0.0) {
                    const double ff = flux[static_cast<size_t>(Param::Nv)];
                    new_mass[static_cast<size_t>(Param::Nv - 1)] -= dt_sub * ff;
                    if (ex_cell != 0.0) {
                        const double ke_high =
                            (sp.vgrid.gamma_cells[Param::Nv - 1] - 1.0)
                            * sp.mass * Const::c * Const::c;
                        energy_current_dt_sum -=
                            dt_sub * ff * ke_high / ex_cell;
                    }
                    if (track_step_diagnostics || track_vmax_loss) {
                        const double loss_high = dt_sub * ff * sg.dx;
                        loss_high_sub += loss_high;
                        if (track_step_diagnostics) {
                            const double ke_high =
                                diag_ke[static_cast<size_t>(Param::Nv - 1)];
                            ke_flux_delta_sub -= loss_high * ke_high;
                        }
                        if (track_high_loss_location) {
                            #pragma omp atomic
                            loss_high_by_x[static_cast<size_t>(ix)] +=
                                loss_high;
                        }
                    }
                }

                for (int iv = 0; iv < Param::Nv; ++iv) {
                    const size_t offset = xbase + static_cast<size_t>(iv) * Param::Nmu
                                        + static_cast<size_t>(imu);
                    const double cell_mass = new_mass[static_cast<size_t>(iv)];
                    if (track_step_diagnostics) {
                        const double px =
                            diag_px_base[static_cast<size_t>(iv)]
                            * sp.vgrid.mu_cells[imu];
                        n_after_sub += cell_mass * sg.dx;
                        px_after_sub += cell_mass * px * sg.dx;
                    }
                    sp.f_tmp[offset] = cell_mass * sp.vgrid.inv_moment_weight[iv];
                }
            }
                if (dt > 0.0) {
                    last_energy_current_cell_x_[static_cast<size_t>(ix)] +=
                        energy_current_dt_sum / dt;
                }
        }
        }

        sp.f.swap(sp.f_tmp);
        if (track_step_diagnostics) {
            last_loss_v_ += n_before_sub - n_after_sub;
            last_loss_v_low_ += loss_low_sub;
        }
        if (track_step_diagnostics || track_vmax_loss) {
            last_loss_v_high_ += loss_high_sub;
        }
        if (track_step_diagnostics) {
            last_mass_error_v_ +=
                n_after_sub + loss_low_sub + loss_high_sub - n_before_sub;
            last_momentum_delta_v_ += px_after_sub - px_before_sub;
            last_energy_delta_v_ += ke_flux_delta_sub;
        }
    }
    close_energy_current_faces(sg, mpi_rank, mpi_size);

    if (track_high_loss_location) {
        int ix_at_max = -1;
        for (int ix = 0; ix < nxl; ++ix) {
            const double loss = loss_high_by_x[static_cast<size_t>(ix)];
            if (loss > last_loss_v_high_local_max_) {
                last_loss_v_high_local_max_ = loss;
                ix_at_max = ix;
            }
        }
        if (ix_at_max >= 0) {
            const int ix_g = ix_at_max + ng;
            const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
            const int iv_umax = Param::Nv - 1;
            const size_t row_umax =
                xbase + static_cast<size_t>(iv_umax) * Param::Nmu;
            double f_umax_mu_integral = 0.0;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                f_umax_mu_integral +=
                    std::max(0.0, sp.f[row_umax + static_cast<size_t>(imu)])
                    * sp.vgrid.dmu;
            }

            double tail_density = 0.0;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                if (sp.vgrid.v_cells[iv] <= Param::diagnostic_tail_u_min) {
                    continue;
                }
                const double shell = sp.vgrid.moment_weight[iv];
                const size_t row =
                    xbase + static_cast<size_t>(iv) * Param::Nmu;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    tail_density +=
                        std::max(0.0, sp.f[row + static_cast<size_t>(imu)])
                        * shell;
                }
            }

            last_x_at_max_loss_v_high_ = sg.x(ix_g);
            last_f_umax_at_max_loss_v_high_ = 0.5 * f_umax_mu_integral;
            last_integral_f_u_gt_8_at_max_loss_v_high_ = tail_density;
        }
    }
}

void VlasovSolver::close_energy_current_faces(const SpatialGrid& sg,
                                              int mpi_rank, int mpi_size)
{
    const int nxl = sg.nx_local;
    if (nxl <= 0) return;
    resize_or_zero(last_energy_current_face_x_, static_cast<size_t>(nxl + 1));
    if (last_energy_current_cell_x_.size() < static_cast<size_t>(nxl)) {
        return;
    }

    if (mpi_size <= 1) {
        last_energy_current_face_x_[0] =
            0.5 * (last_energy_current_cell_x_[static_cast<size_t>(nxl - 1)]
                 + last_energy_current_cell_x_[0]);
        for (int iface = 1; iface < nxl; ++iface) {
            last_energy_current_face_x_[static_cast<size_t>(iface)] =
                0.5 * (last_energy_current_cell_x_[static_cast<size_t>(iface - 1)]
                     + last_energy_current_cell_x_[static_cast<size_t>(iface)]);
        }
        last_energy_current_face_x_[static_cast<size_t>(nxl)] =
            last_energy_current_face_x_[0];
        return;
    }

    const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right_peer = (mpi_rank + 1) % mpi_size;
    const double send_left_cell = last_energy_current_cell_x_[0];
    const double send_right_cell =
        last_energy_current_cell_x_[static_cast<size_t>(nxl - 1)];
    double recv_left_cell = 0.0;
    double recv_right_cell = 0.0;
    MPI_Sendrecv(&send_right_cell, 1, MPI_DOUBLE, right_peer, 621,
                 &recv_left_cell, 1, MPI_DOUBLE, left_peer, 621,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&send_left_cell, 1, MPI_DOUBLE, left_peer, 622,
                 &recv_right_cell, 1, MPI_DOUBLE, right_peer, 622,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    last_energy_current_face_x_[0] =
        0.5 * (recv_left_cell + last_energy_current_cell_x_[0]);
    for (int iface = 1; iface < nxl; ++iface) {
        last_energy_current_face_x_[static_cast<size_t>(iface)] =
            0.5 * (last_energy_current_cell_x_[static_cast<size_t>(iface - 1)]
                 + last_energy_current_cell_x_[static_cast<size_t>(iface)]);
    }
    last_energy_current_face_x_[static_cast<size_t>(nxl)] =
        0.5 * (last_energy_current_cell_x_[static_cast<size_t>(nxl - 1)]
             + recv_right_cell);
}

void VlasovSolver::advect_mu(Species& sp, const SpatialGrid& sg,
                             const EMFields& fields, double dt)
{
    const bool track_step_diagnostics = step_diagnostics_enabled_;
    last_loss_mu_ = 0.0;
    last_mass_error_mu_ = 0.0;
    last_momentum_delta_mu_ = 0.0;
    last_energy_delta_mu_ = 0.0;

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    double max_cfl = 0.0;

    #pragma omp parallel for collapse(2) schedule(static) reduction(max:max_cfl)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const int ix_g = ix + ng;
            const double ex_cell = staggered_cell_ex(fields, ix, ix_g);
            const double udot_abs =
                std::fabs(sp.charge * ex_cell / (sp.mass * Const::c));
            max_cfl = std::max(max_cfl,
                               dt * udot_abs * sp.vgrid.mu_cfl_factor[iv]);
        }
    }

    last_nsub_mu_ = substeps_from_cfl(max_cfl, Param::velocity_space_cfl);
    last_cfl_mu_ = max_cfl / last_nsub_mu_;
    const double dt_sub = dt / last_nsub_mu_;
    std::array<double, Param::Nv> diag_px_base;
    std::array<double, Param::Nv> diag_ke;
    if (track_step_diagnostics) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            diag_px_base[static_cast<size_t>(iv)] =
                sp.mass * Const::c * sp.vgrid.v_cells[iv];
            diag_ke[static_cast<size_t>(iv)] =
                (sp.vgrid.gamma_cells[iv] - 1.0)
                * sp.mass * Const::c * Const::c;
        }
    }

    for (int isub = 0; isub < last_nsub_mu_; ++isub) {
        double n_before_sub = 0.0;
        double n_after_sub = 0.0;
        double px_before_sub = 0.0;
        double px_after_sub = 0.0;
        double ke_before_sub = 0.0;
        double ke_after_sub = 0.0;

        #pragma omp parallel reduction(+:n_before_sub,n_after_sub, \
                                         px_before_sub,px_after_sub,ke_before_sub,ke_after_sub)
        {
            std::array<double, Param::Nmu> mass;
            std::array<double, Param::Nmu> new_mass;
            std::array<double, Param::Nmu + 1> flux;
            std::array<double, Param::Nmu> outgoing;
            std::array<double, Param::Nmu> scale;

            #pragma omp for collapse(2) schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int iv = 0; iv < Param::Nv; ++iv) {
                const int ix_g = ix + ng;
                const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
                const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
                const double ex_cell = staggered_cell_ex(fields, ix, ix_g);
                const double accel_u =
                    sp.charge * ex_cell / (sp.mass * Const::c);
                const double u = sp.vgrid.v_cells[iv];
                const double u_eff = std::max(u, Param::u_floor);
                const double mu_flux_scale = sp.vgrid.mu_flux_scale[iv];
                flux.fill(0.0);
                outgoing.fill(0.0);
                scale.fill(1.0);
                const double ke = track_step_diagnostics
                    ? diag_ke[static_cast<size_t>(iv)] : 0.0;
                const double px_base = track_step_diagnostics
                    ? diag_px_base[static_cast<size_t>(iv)] : 0.0;

                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const size_t offset = row + static_cast<size_t>(imu);
                    const double cell_mass =
                        std::max(0.0, sp.f[offset]) * sp.vgrid.moment_weight[iv];
                    mass[static_cast<size_t>(imu)] = cell_mass;
                    new_mass[static_cast<size_t>(imu)] = cell_mass;
                    if (track_step_diagnostics) {
                        n_before_sub += cell_mass * sg.dx;
                        px_before_sub +=
                            cell_mass * px_base * sp.vgrid.mu_cells[imu]
                            * sg.dx;
                        ke_before_sub += cell_mass * ke * sg.dx;
                    }
                }

                for (int face = 1; face < Param::Nmu; ++face) {
                    const double mudot =
                        accel_u * sp.vgrid.mu_face_factor[face] / u_eff;
                    const double f_face =
                        reconstruct_mu_face(sp, row, face, mudot);
                    flux[static_cast<size_t>(face)] =
                        mu_flux_scale * mudot * f_face;
                    const int donor =
                        (flux[static_cast<size_t>(face)] >= 0.0) ? face - 1 : face;
                    outgoing[static_cast<size_t>(donor)] +=
                        dt_sub * std::fabs(flux[static_cast<size_t>(face)]);
                }

                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double out = outgoing[static_cast<size_t>(imu)];
                    const double m = mass[static_cast<size_t>(imu)];
                    if (out > m && out > 0.0) {
                        scale[static_cast<size_t>(imu)] = m / out;
                    }
                }

                for (int face = 1; face < Param::Nmu; ++face) {
                    double& ff = flux[static_cast<size_t>(face)];
                    if (ff == 0.0) continue;
                    const int donor = (ff >= 0.0) ? face - 1 : face;
                    ff *= scale[static_cast<size_t>(donor)];
                    new_mass[static_cast<size_t>(face - 1)] -= dt_sub * ff;
                    new_mass[static_cast<size_t>(face)] += dt_sub * ff;
                }

                const double inv_shell = sp.vgrid.inv_moment_weight[iv];
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const size_t offset = row + static_cast<size_t>(imu);
                    const double cell_mass =
                        (new_mass[static_cast<size_t>(imu)] > 0.0)
                        ? new_mass[static_cast<size_t>(imu)] : 0.0;
                    if (track_step_diagnostics) {
                        n_after_sub += cell_mass * sg.dx;
                        px_after_sub +=
                            cell_mass * px_base * sp.vgrid.mu_cells[imu]
                            * sg.dx;
                        ke_after_sub += cell_mass * ke * sg.dx;
                    }
                    sp.f_tmp[offset] = cell_mass * inv_shell;
                }
            }
        }
        }

        sp.f.swap(sp.f_tmp);
        if (track_step_diagnostics) {
            last_loss_mu_ += n_before_sub - n_after_sub;
            last_mass_error_mu_ += n_after_sub - n_before_sub;
            last_momentum_delta_mu_ += px_after_sub - px_before_sub;
            last_energy_delta_mu_ += ke_after_sub - ke_before_sub;
        }
    }
}

void VlasovSolver::exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                                     int mpi_rank, int mpi_size)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;
    size_t slice_size = Param::Nvmu;
    int left_rank = (mpi_rank + mpi_size - 1) % mpi_size;
    int right_rank = (mpi_rank + 1) % mpi_size;

    if (mpi_size == 1) {
        fill_periodic_ghosts_single_rank(sp, sg);
        return;
    }

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
    MPI_Isend(send_left_.data(), (int)buffer_size, MPI_DOUBLE,
              left_rank, 101, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(recv_left_.data(), (int)buffer_size, MPI_DOUBLE,
              left_rank, 102, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Isend(send_right_.data(), (int)buffer_size, MPI_DOUBLE,
              right_rank, 102, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(recv_right_.data(), (int)buffer_size, MPI_DOUBLE,
              right_rank, 101, MPI_COMM_WORLD, &reqs[nreq++]);
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    std::memcpy(&sp.f[0], recv_left_.data(), buffer_size * sizeof(double));
    std::memcpy(&sp.f[static_cast<size_t>(ng + nxl) * slice_size],
                recv_right_.data(), buffer_size * sizeof(double));
}
