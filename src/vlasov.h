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
    double last_loss_mu() const { return last_loss_mu_; }

private:
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);

    std::vector<double> send_left_;
    std::vector<double> send_right_;
    std::vector<double> recv_left_;
    std::vector<double> recv_right_;
    std::vector<RemapStencil> x_stencil_;

    double last_cfl_v_;
    double last_cfl_mu_;
    double last_loss_v_;
    double last_loss_v_low_;
    double last_loss_v_high_;
    double last_loss_mu_;
    int last_nsub_v_;
    int last_nsub_mu_;
    bool step_diagnostics_enabled_;
};

#endif
