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
    void update_dynamic_reservoir(Species& sp,
                                  const SpatialGrid& sg,
                                  double control_dt,
                                  int mpi_rank,
                                  int mpi_size);

    void advect_x(Species& sp, const SpatialGrid& sg, double dt,
                  int mpi_rank, int mpi_size);
    void configure_boundary_reservoir_inflow(Species& sp,
                                             double target_number,
                                             double dt,
                                             int mpi_rank,
                                             int mpi_size);
    void advect_v(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                  double dt);
    void advect_mu(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                   double dt);

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
    double cached_incoming_flux_per_density(const Species& sp,
                                            bool left_boundary);
    double bounded_density_from_flux(const Species& sp,
                                     double target_flux,
                                     bool left_boundary);

    std::vector<double> send_left_;
    std::vector<double> send_right_;
    std::vector<double> recv_left_;
    std::vector<double> recv_right_;
    std::vector<double> reservoir_left_;
    std::vector<double> reservoir_right_;
    std::vector<double> unit_reservoir_;
    std::vector<RemapStencil> x_stencil_;
    std::vector<double> x_cfl_;

    double cached_unit_flux_left_;
    double cached_unit_flux_right_;
    double cached_unit_flux_drift_left_;
    double cached_unit_flux_drift_right_;
    double cached_global_unit_flux_sum_;
    double cached_global_unit_flux_drift_left_;
    double cached_global_unit_flux_drift_right_;
    int cached_global_unit_flux_mpi_size_;
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
    bool unit_flux_left_valid_;
    bool unit_flux_right_valid_;
    bool global_unit_flux_valid_;
    bool step_diagnostics_enabled_;
};

#endif
