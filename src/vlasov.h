#ifndef VLASOV_H
#define VLASOV_H

#include "species.h"
#include "grid.h"
#include <vector>

struct EMFields;

class VlasovSolver {
public:
    void advect(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                double dt, int mpi_rank, int mpi_size);

private:
    void advect_x(Species& sp, const SpatialGrid& sg, double dt,
                  int mpi_rank, int mpi_size);
    void advect_v(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                  double dt);
    void advect_mu(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                   double dt);
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);

    std::vector<double> send_left_;
    std::vector<double> send_right_;
    std::vector<double> recv_left_;
    std::vector<double> recv_right_;
};

#endif
