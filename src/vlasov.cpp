#include "vlasov.h"
#include "maxwell.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include <omp.h>

namespace {
void fill_left_physical_ghosts(Species& sp, const SpatialGrid& sg,
                               const std::vector<double>& incoming)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    if (nxl <= 0) return;
    const size_t source = static_cast<size_t>(ng) * Param::Nvmu;
    for (int g = 0; g < ng; ++g) {
        const size_t dst = static_cast<size_t>(ng - 1 - g) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx > 0.0) {
                sp.f[dst + k] =
                    (incoming.size() == Param::Nvmu) ? incoming[k] : 0.0;
            } else {
                sp.f[dst + k] = sp.f[source + k];
            }
        }
    }
}

void fill_right_physical_ghosts(Species& sp, const SpatialGrid& sg,
                                const std::vector<double>& incoming)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    if (nxl <= 0) return;
    const size_t source = static_cast<size_t>(ng + nxl - 1) * Param::Nvmu;
    for (int g = 0; g < ng; ++g) {
        const size_t dst = static_cast<size_t>(ng + nxl + g) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx < 0.0) {
                sp.f[dst + k] =
                    (incoming.size() == Param::Nvmu) ? incoming[k] : 0.0;
            } else {
                sp.f[dst + k] = sp.f[source + k];
            }
        }
    }
}

void local_boundary_outflow_moments(const Species& sp,
                                    const SpatialGrid& sg,
                                    int mpi_rank,
                                    int mpi_size,
                                    double& left_out,
                                    double& right_out,
                                    double& px_out,
                                    double& energy_out)
{
    left_out = 0.0;
    right_out = 0.0;
    px_out = 0.0;
    energy_out = 0.0;
    if (sp.type != SpeciesType::BACKGROUND_ELECTRON || sg.nx_local <= 0) {
        return;
    }

    const int ng = sg.nghost;
    if (mpi_rank == 0) {
        const size_t xbase = static_cast<size_t>(ng) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx >= 0.0) continue;
            const int iv = static_cast<int>(k / Param::Nmu);
            const int imu = static_cast<int>(k % Param::Nmu);
            const double flux =
                -vx * sp.f[xbase + k] * sp.vgrid.moment_weight[iv];
            const double u = sp.vgrid.v_cells[iv];
            const double px =
                sp.mass * Const::c * u * sp.vgrid.mu_cells[imu];
            const double ke =
                (gamma_from_u(u) - 1.0) * sp.mass * Const::c * Const::c;
            left_out += flux;
            px_out += flux * px;
            energy_out += flux * ke;
        }
    }

    if (mpi_rank == mpi_size - 1) {
        const size_t xbase =
            static_cast<size_t>(ng + sg.nx_local - 1) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx <= 0.0) continue;
            const int iv = static_cast<int>(k / Param::Nmu);
            const int imu = static_cast<int>(k % Param::Nmu);
            const double flux =
                vx * sp.f[xbase + k] * sp.vgrid.moment_weight[iv];
            const double u = sp.vgrid.v_cells[iv];
            const double px =
                sp.mass * Const::c * u * sp.vgrid.mu_cells[imu];
            const double ke =
                (gamma_from_u(u) - 1.0) * sp.mass * Const::c * Const::c;
            right_out += flux;
            px_out += flux * px;
            energy_out += flux * ke;
        }
    }
}

void local_boundary_outflow_fluxes(const Species& sp,
                                   const SpatialGrid& sg,
                                   int mpi_rank,
                                   int mpi_size,
                                   double& left_out,
                                   double& right_out)
{
    double px_out = 0.0;
    double energy_out = 0.0;
    local_boundary_outflow_moments(sp, sg, mpi_rank, mpi_size,
                                   left_out, right_out,
                                   px_out, energy_out);
}

void local_boundary_inflow_fluxes(const Species& sp,
                                  const SpatialGrid& sg,
                                  int mpi_rank,
                                  int mpi_size,
                                  double& left_in,
                                  double& right_in)
{
    left_in = 0.0;
    right_in = 0.0;
    if (sp.type != SpeciesType::BACKGROUND_ELECTRON || sg.nx_local <= 0) {
        return;
    }

    const int ng = sg.nghost;
    if (mpi_rank == 0) {
        const size_t xbase = static_cast<size_t>(ng - 1) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx <= 0.0) continue;
            const int iv = static_cast<int>(k / Param::Nmu);
            left_in += vx * sp.f[xbase + k] * sp.vgrid.moment_weight[iv];
        }
    }

    if (mpi_rank == mpi_size - 1) {
        const size_t xbase =
            static_cast<size_t>(ng + sg.nx_local) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = sp.vgrid.vx_cells[k];
            if (vx >= 0.0) continue;
            const int iv = static_cast<int>(k / Param::Nmu);
            right_in += -vx * sp.f[xbase + k] * sp.vgrid.moment_weight[iv];
        }
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

}

VlasovSolver::VlasovSolver()
    : cached_unit_flux_left_(0.0),
      cached_unit_flux_right_(0.0),
      cached_unit_flux_drift_left_(0.0),
      cached_unit_flux_drift_right_(0.0),
      cached_global_unit_flux_sum_(0.0),
      cached_global_unit_flux_drift_left_(0.0),
      cached_global_unit_flux_drift_right_(0.0),
      cached_global_unit_flux_mpi_size_(0),
      last_cfl_v_(0.0),
      last_cfl_mu_(0.0),
      last_loss_v_(0.0),
      last_loss_v_low_(0.0),
      last_loss_v_high_(0.0),
      last_loss_v_high_local_max_(0.0),
      last_x_at_max_loss_v_high_(0.0),
      last_f_umax_at_max_loss_v_high_(0.0),
      last_integral_f_u_gt_8_at_max_loss_v_high_(0.0),
      last_loss_mu_(0.0),
      last_loss_x_left_(0.0),
      last_loss_x_right_(0.0),
      last_inflow_x_left_(0.0),
      last_inflow_x_right_(0.0),
      last_loss_x_momentum_(0.0),
      last_loss_x_energy_(0.0),
      last_mass_error_v_(0.0),
      last_mass_error_mu_(0.0),
      last_momentum_delta_v_(0.0),
      last_momentum_delta_mu_(0.0),
      last_energy_delta_v_(0.0),
      last_energy_delta_mu_(0.0),
      last_nsub_v_(1),
      last_nsub_mu_(1),
      unit_flux_left_valid_(false),
      unit_flux_right_valid_(false),
      global_unit_flux_valid_(false),
      step_diagnostics_enabled_(false)
{}

void VlasovSolver::advect(Species& sp, const SpatialGrid& sg,
                          const EMFields& fields, double dt,
                          int mpi_rank, int mpi_size)
{
    advect_x(sp, sg, 0.5 * dt, mpi_rank, mpi_size);
    advect_v(sp, sg, fields, 0.5 * dt);
    advect_mu(sp, sg, fields, dt);
    advect_v(sp, sg, fields, 0.5 * dt);
    advect_x(sp, sg, 0.5 * dt, mpi_rank, mpi_size);
}

double VlasovSolver::cached_incoming_flux_per_density(const Species& sp,
                                                      bool left_boundary)
{
    const double drift = left_boundary
        ? sp.reservoir_drift_left
        : sp.reservoir_drift_right;
    bool& valid = left_boundary ? unit_flux_left_valid_ : unit_flux_right_valid_;
    double& cached_drift = left_boundary
        ? cached_unit_flux_drift_left_
        : cached_unit_flux_drift_right_;
    double& cached_flux = left_boundary
        ? cached_unit_flux_left_
        : cached_unit_flux_right_;

    if (valid && cached_drift == drift) return cached_flux;

    sp.fill_maxwellian_velocity_slice(unit_reservoir_, 1.0,
                                      sp.temperature, drift);
    double flux = 0.0;
    for (size_t k = 0; k < Param::Nvmu; ++k) {
        const double vx = sp.vgrid.vx_cells[k];
        const bool incoming = left_boundary ? (vx > 0.0) : (vx < 0.0);
        if (!incoming) continue;
        const int iv = static_cast<int>(k / Param::Nmu);
        const double speed = left_boundary ? vx : -vx;
        flux += speed * unit_reservoir_[k] * sp.vgrid.moment_weight[iv];
    }

    cached_drift = drift;
    cached_flux = flux;
    valid = true;
    return cached_flux;
}

double VlasovSolver::bounded_density_from_flux(const Species& sp,
                                               double target_flux,
                                               bool left_boundary)
{
    const double unit_flux =
        cached_incoming_flux_per_density(sp, left_boundary);
    if (!(unit_flux > 0.0) || !(target_flux > 0.0)) {
        return Param::background_reservoir_min_density_factor * sp.density0;
    }

    const double min_density =
        Param::background_reservoir_min_density_factor * sp.density0;
    const double max_density =
        Param::background_reservoir_max_density_factor * sp.density0;
    return std::max(min_density,
                    std::min(max_density, target_flux / unit_flux));
}

void VlasovSolver::update_dynamic_reservoir(Species& sp,
                                            const SpatialGrid& sg,
                                            double control_dt,
                                            int mpi_rank,
                                            int mpi_size)
{
    if (!Param::enable_dynamic_background_reservoir ||
        sp.type != SpeciesType::BACKGROUND_ELECTRON ||
        sg.nx_local <= 0 || !(control_dt > 0.0)) {
        return;
    }

    double local_left_out = 0.0;
    double local_right_out = 0.0;
    local_boundary_outflow_fluxes(sp, sg, mpi_rank, mpi_size,
                                  local_left_out, local_right_out);

    double global_values[2] = { 0.0, 0.0 };
    double local_values[2] = {
        local_left_out,
        local_right_out
    };
    MPI_Allreduce(local_values, global_values, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    const double target_density_left =
        bounded_density_from_flux(sp, global_values[0], true);
    const double target_density_right =
        bounded_density_from_flux(sp, global_values[1], false);
    const double relaxation =
        std::max(0.0,
                 std::min(1.0,
                          Param::background_reservoir_density_relaxation));

    sp.reservoir_density_left =
        sp.reservoir_density_left
        + relaxation * (target_density_left - sp.reservoir_density_left);
    sp.reservoir_density_right =
        sp.reservoir_density_right
        + relaxation * (target_density_right - sp.reservoir_density_right);
}

void VlasovSolver::configure_boundary_reservoir_inflow(Species& sp,
                                                       double target_number,
                                                       double dt,
                                                       int mpi_rank,
                                                       int mpi_size)
{
    const bool owns_left_boundary = (mpi_rank == 0);
    const bool owns_right_boundary = (mpi_rank == mpi_size - 1);

    if (!Param::enable_dynamic_background_reservoir ||
        sp.type != SpeciesType::BACKGROUND_ELECTRON ||
        !(target_number > 0.0) || !(dt > 0.0)) {
        if (owns_left_boundary && sp.reservoir_density_left != 0.0 &&
            reservoir_left_.size() == Param::Nvmu) {
            std::fill(reservoir_left_.begin(), reservoir_left_.end(), 0.0);
        }
        if (owns_right_boundary && sp.reservoir_density_right != 0.0 &&
            reservoir_right_.size() == Param::Nvmu) {
            std::fill(reservoir_right_.begin(), reservoir_right_.end(), 0.0);
        }
        sp.reservoir_density_left = 0.0;
        sp.reservoir_density_right = 0.0;
        return;
    }

    const bool global_unit_cache_current =
        global_unit_flux_valid_ &&
        cached_global_unit_flux_mpi_size_ == mpi_size &&
        cached_global_unit_flux_drift_left_ == sp.reservoir_drift_left &&
        cached_global_unit_flux_drift_right_ == sp.reservoir_drift_right;
    double global_unit_sum = cached_global_unit_flux_sum_;
    if (!global_unit_cache_current) {
        const double unit_left = owns_left_boundary
            ? cached_incoming_flux_per_density(sp, true) : 0.0;
        const double unit_right = owns_right_boundary
            ? cached_incoming_flux_per_density(sp, false) : 0.0;
        double local_unit_sum = unit_left + unit_right;
        MPI_Allreduce(&local_unit_sum, &global_unit_sum, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        cached_global_unit_flux_sum_ = global_unit_sum;
        cached_global_unit_flux_drift_left_ = sp.reservoir_drift_left;
        cached_global_unit_flux_drift_right_ = sp.reservoir_drift_right;
        cached_global_unit_flux_mpi_size_ = mpi_size;
        global_unit_flux_valid_ = true;
    }
    if (!(global_unit_sum > 0.0)) return;

    const double max_density =
        Param::background_reservoir_max_density_factor * sp.density0;
    const double target_density =
        std::min(max_density, target_number / (dt * global_unit_sum));

    sp.reservoir_density_left = owns_left_boundary ? target_density : 0.0;
    sp.reservoir_density_right = owns_right_boundary ? target_density : 0.0;

    if (owns_left_boundary) {
        sp.fill_maxwellian_velocity_slice(reservoir_left_,
                                          sp.reservoir_density_left,
                                          sp.temperature,
                                          sp.reservoir_drift_left);
    }
    if (owns_right_boundary) {
        sp.fill_maxwellian_velocity_slice(reservoir_right_,
                                          sp.reservoir_density_right,
                                          sp.temperature,
                                          sp.reservoir_drift_right);
    }
}

void VlasovSolver::advect_x(Species& sp, const SpatialGrid& sg, double dt,
                            int mpi_rank, int mpi_size)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const double max_cfl = sp.vgrid.max_speed * sp.vgrid.max_abs_mu * dt / sg.dx;
    const double cfl_limit = std::min(0.85, Param::semi_lagrangian_cfl);
    const int nsub_x = substeps_from_cfl(max_cfl, cfl_limit);
    const double dt_sub = dt / nsub_x;
    if (x_cfl_.size() != Param::Nvmu) {
        x_cfl_.resize(Param::Nvmu);
    }
    for (size_t k = 0; k < Param::Nvmu; ++k) {
        x_cfl_[k] = sp.vgrid.vx_cells[k] * dt_sub / sg.dx;
    }
    last_loss_x_left_ = 0.0;
    last_loss_x_right_ = 0.0;
    last_inflow_x_left_ = 0.0;
    last_inflow_x_right_ = 0.0;
    last_loss_x_momentum_ = 0.0;
    last_loss_x_energy_ = 0.0;

    for (int isub = 0; isub < nsub_x; ++isub) {
        exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);
        double out_left = 0.0;
        double out_right = 0.0;
        double px_out = 0.0;
        double energy_out = 0.0;
        double in_left = 0.0;
        double in_right = 0.0;
        local_boundary_inflow_fluxes(sp, sg, mpi_rank, mpi_size,
                                     in_left, in_right);
        local_boundary_outflow_moments(sp, sg, mpi_rank, mpi_size,
                                       out_left, out_right,
                                       px_out, energy_out);
        last_inflow_x_left_ += in_left * dt_sub;
        last_inflow_x_right_ += in_right * dt_sub;
        last_loss_x_left_ += out_left * dt_sub;
        last_loss_x_right_ += out_right * dt_sub;
        last_loss_x_momentum_ += px_out * dt_sub;
        last_loss_x_energy_ += energy_out * dt_sub;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int ix = ng; ix < ng + nxl; ++ix) {
            for (int iv = 0; iv < Param::Nv; ++iv) {
                const size_t xbase = static_cast<size_t>(ix) * Param::Nvmu;
                const size_t xbase_left = static_cast<size_t>(ix - 1) * Param::Nvmu;
                const size_t xbase_right = static_cast<size_t>(ix + 1) * Param::Nvmu;
                const size_t row = static_cast<size_t>(iv) * Param::Nmu;
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const size_t k = row + imu;
                    const size_t offset = xbase + k;
                    const double cfl = x_cfl_[k];
                    if (cfl >= 0.0) {
                        sp.f_tmp[offset] =
                            std::max(0.0, (1.0 - cfl) * sp.f[offset]
                                           + cfl * sp.f[xbase_left + k]);
                    } else {
                        const double a = -cfl;
                        sp.f_tmp[offset] =
                            std::max(0.0, (1.0 - a) * sp.f[offset]
                                           + a * sp.f[xbase_right + k]);
                    }
                }
            }
        }

        sp.f.swap(sp.f_tmp);
    }
}

void VlasovSolver::advect_v(Species& sp, const SpatialGrid& sg,
                            const EMFields& fields, double dt)
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
            const double udot_abs =
                std::fabs(sp.charge * fields.Ex[ix_g] / (sp.mass * Const::c));
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
        double ke_before_sub = 0.0;
        double ke_after_sub = 0.0;

        #pragma omp parallel reduction(+:n_before_sub,n_after_sub,loss_low_sub,loss_high_sub, \
                                         px_before_sub,px_after_sub,ke_before_sub,ke_after_sub)
        {
            std::array<double, Param::Nv> mass;
            std::array<double, Param::Nv> new_mass;
            std::array<double, Param::Nv + 1> flux;
            std::array<double, Param::Nv> outgoing;
            std::array<double, Param::Nv> scale;

            #pragma omp for collapse(2) schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                const int ix_g = ix + ng;
                const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
                const double accel_u =
                    sp.charge * fields.Ex[ix_g] / (sp.mass * Const::c);
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
                        const double ke = diag_ke[static_cast<size_t>(iv)];
                        n_before_sub += cell_mass * sg.dx;
                        px_before_sub += cell_mass * px * sg.dx;
                        ke_before_sub += cell_mass * ke * sg.dx;
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
                    new_mass[static_cast<size_t>(face - 1)] -= dt_sub * ff;
                    new_mass[static_cast<size_t>(face)] += dt_sub * ff;
                }

                if (flux[static_cast<size_t>(Param::Nv)] > 0.0) {
                    const double ff = flux[static_cast<size_t>(Param::Nv)];
                    new_mass[static_cast<size_t>(Param::Nv - 1)] -= dt_sub * ff;
                    if (track_step_diagnostics || track_vmax_loss) {
                        const double loss_high = dt_sub * ff * sg.dx;
                        loss_high_sub += loss_high;
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
                    const double cell_mass =
                        (new_mass[static_cast<size_t>(iv)] > 0.0)
                        ? new_mass[static_cast<size_t>(iv)] : 0.0;
                    if (track_step_diagnostics) {
                        const double px =
                            diag_px_base[static_cast<size_t>(iv)]
                            * sp.vgrid.mu_cells[imu];
                        const double ke = diag_ke[static_cast<size_t>(iv)];
                        n_after_sub += cell_mass * sg.dx;
                        px_after_sub += cell_mass * px * sg.dx;
                        ke_after_sub += cell_mass * ke * sg.dx;
                    }
                    sp.f_tmp[offset] = cell_mass * sp.vgrid.inv_moment_weight[iv];
                }
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
            last_energy_delta_v_ += ke_after_sub - ke_before_sub;
        }
    }

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
            const double udot_abs =
                std::fabs(sp.charge * fields.Ex[ix_g] / (sp.mass * Const::c));
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
                const double accel_u =
                    sp.charge * fields.Ex[ix_g] / (sp.mass * Const::c);
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
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (left_rank >= 0) {
        std::memcpy(&sp.f[0], recv_left_.data(), buffer_size * sizeof(double));
    } else {
        fill_left_physical_ghosts(sp, sg, reservoir_left_);
    }

    if (right_rank < mpi_size) {
        std::memcpy(&sp.f[static_cast<size_t>(ng + nxl) * slice_size],
                    recv_right_.data(), buffer_size * sizeof(double));
    } else {
        fill_right_physical_ghosts(sp, sg, reservoir_right_);
    }
}
