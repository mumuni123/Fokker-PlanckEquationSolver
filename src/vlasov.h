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
                                  const EMFields& fields,
                                  double beam_injected_number,
                                  double beam_outflow_number,
                                  double control_dt,
                                  int mpi_rank,
                                  int mpi_size);

    void advect_x(Species& sp, const SpatialGrid& sg, double dt,
                  int mpi_rank, int mpi_size);
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
    double last_mass_error_v() const { return last_mass_error_v_; }
    double last_mass_error_mu() const { return last_mass_error_mu_; }
    double last_momentum_delta_v() const { return last_momentum_delta_v_; }
    double last_momentum_delta_mu() const { return last_momentum_delta_mu_; }
    double last_energy_delta_v() const { return last_energy_delta_v_; }
    double last_energy_delta_mu() const { return last_energy_delta_mu_; }

private:
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);
    void update_reservoir_cache(const Species& sp);
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

    double cached_reservoir_density_left_;
    double cached_reservoir_density_right_;
    double cached_reservoir_drift_left_;
    double cached_reservoir_drift_right_;
    double cached_unit_flux_left_;
    double cached_unit_flux_right_;
    double cached_unit_flux_drift_left_;
    double cached_unit_flux_drift_right_;
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
    double last_mass_error_v_;
    double last_mass_error_mu_;
    double last_momentum_delta_v_;
    double last_momentum_delta_mu_;
    double last_energy_delta_v_;
    double last_energy_delta_mu_;
    int last_nsub_v_;
    int last_nsub_mu_;
    bool reservoir_cache_valid_;
    bool unit_flux_left_valid_;
    bool unit_flux_right_valid_;
    bool step_diagnostics_enabled_;
};

#endif
