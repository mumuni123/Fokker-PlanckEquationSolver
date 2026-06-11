#ifndef VLASOV_H
#define VLASOV_H

#include "species.h"
#include "grid.h"
#include <vector>

struct EMFields;

struct RemapStencil {
    int first;
    int count;
    double weight[6];
};

class VlasovSolver {
public:
    VlasovSolver();

    void set_step_diagnostics_enabled(bool enabled) {
        step_diagnostics_enabled_ = enabled;
    }

    void advect(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                double dt, int mpi_rank, int mpi_size);
    void advect_x(Species& sp, const SpatialGrid& sg, double dt,
                  int mpi_rank, int mpi_size, double time = 0.0);
    void advect_v(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                  double dt);
    void advect_mu(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                   double dt);
    void update_upstream_phase_feedback(double time,
                                        double ex_left,
                                        double ne_left);

    double last_cfl_v() const { return last_cfl_v_; }
    double last_cfl_mu() const { return last_cfl_mu_; }
    int last_nsub_v() const { return last_nsub_v_; }
    int last_nsub_mu() const { return last_nsub_mu_; }
    double last_loss_v() const { return last_loss_v_; }
    double last_loss_v_low() const { return last_loss_v_low_; }
    double last_loss_v_high() const { return last_loss_v_high_; }
    double last_loss_v_high_local_max() const {
        return last_loss_v_high_local_max_;
    }
    double last_x_at_max_loss_v_high() const {
        return last_x_at_max_loss_v_high_;
    }
    double last_f_umax_at_max_loss_v_high() const {
        return last_f_umax_at_max_loss_v_high_;
    }
    double last_integral_f_u_gt_8_at_max_loss_v_high() const {
        return last_integral_f_u_gt_8_at_max_loss_v_high_;
    }
    double last_loss_mu() const { return last_loss_mu_; }
    double last_loss_x_left() const { return last_loss_x_left_; }
    double last_loss_x_right() const { return last_loss_x_right_; }
    double last_inflow_x_left() const { return last_inflow_x_left_; }
    double last_inflow_x_right() const { return last_inflow_x_right_; }
    double last_boundary_inflow() const {
        return last_inflow_x_left_ + last_inflow_x_right_;
    }
    double last_loss_x_momentum() const { return last_loss_x_momentum_; }
    double last_loss_x_energy() const { return last_loss_x_energy_; }
    double last_mass_error_v() const { return last_mass_error_v_; }
    double last_mass_error_mu() const { return last_mass_error_mu_; }
    double last_momentum_delta_v() const { return last_momentum_delta_v_; }
    double last_momentum_delta_mu() const { return last_momentum_delta_mu_; }
    double last_energy_delta_v() const { return last_energy_delta_v_; }
    double last_energy_delta_mu() const { return last_energy_delta_mu_; }

private:
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);
    void update_open_boundary_inflow(const Species& sp, double time,
                                     bool owns_left_boundary,
                                     bool owns_right_boundary);
    void ensure_upstream_basis(const Species& sp);
    void update_flux_balance(double in_left, double in_right,
                             double out_left, double out_right,
                             double dt_sub, int mpi_size);

    std::vector<double> send_left_;
    std::vector<double> send_right_;
    std::vector<double> recv_left_;
    std::vector<double> recv_right_;
    std::vector<double> upstream_left_;
    std::vector<double> upstream_right_;
    std::vector<double> upstream_base_;
    std::vector<double> upstream_current_shape_;
    std::vector<double> upstream_temperature_shape_;
    double upstream_basis_density_;
    double upstream_basis_temperature_;
    double upstream_basis_mass_;
    double upstream_left_cached_density_;
    double upstream_left_cached_temperature_;
    double upstream_left_cached_drift_;
    double upstream_right_cached_density_;
    double upstream_right_cached_temperature_;
    double upstream_right_cached_drift_;
    bool upstream_left_cache_valid_;
    bool upstream_right_cache_valid_;
    bool upstream_basis_valid_;
    double upstream_flux_in_left_avg_;
    double upstream_flux_in_right_avg_;
    double upstream_flux_out_avg_;
    double upstream_flux_correction_;
    double upstream_phase_feedback_;
    std::vector<RemapStencil> x_stencil_;
    std::vector<double> x_cfl_;

    double last_cfl_v_;
    double last_cfl_mu_;
    double last_loss_v_;
    double last_loss_v_low_;
    double last_loss_v_high_;
    double last_loss_v_high_local_max_;
    double last_x_at_max_loss_v_high_;
    double last_f_umax_at_max_loss_v_high_;
    double last_integral_f_u_gt_8_at_max_loss_v_high_;
    double last_loss_mu_;
    double last_loss_x_left_;
    double last_loss_x_right_;
    double last_inflow_x_left_;
    double last_inflow_x_right_;
    double last_loss_x_momentum_;
    double last_loss_x_energy_;
    double last_mass_error_v_;
    double last_mass_error_mu_;
    double last_momentum_delta_v_;
    double last_momentum_delta_mu_;
    double last_energy_delta_v_;
    double last_energy_delta_mu_;
    int last_nsub_v_;
    int last_nsub_mu_;
    bool step_diagnostics_enabled_;
};

#endif
