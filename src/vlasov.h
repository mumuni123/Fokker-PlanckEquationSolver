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

private:
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);

    std::vector<double> send_left_;
    std::vector<double> send_right_;
    std::vector<double> recv_left_;
    std::vector<double> recv_right_;
    std::vector<RemapStencil> x_stencil_;
    std::vector<RemapStencil> v_stencil_;
    std::vector<RemapStencil> mu_stencil_;

    double last_cfl_v_;
    double last_cfl_mu_;
    int last_nsub_v_;
    int last_nsub_mu_;
};

#endif
