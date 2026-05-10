#ifndef MAXWELL_H
#define MAXWELL_H

#include "grid.h"
#include "parameters.h"
#include <vector>

struct EMFields {
    std::vector<double> Ex;
    std::vector<double> rho;

    int nx_total;
    double dx;

    void init(const SpatialGrid& sg);
    void zero_currents();
    void accumulate_moments(const class Species& sp);
    void set_charge_density(const class Species& electrons,
                            const std::vector<double>& beam_density);

    // Electrostatic field from Gauss law: dE/dx = rho / eps0.
    void solve_poisson(int mpi_rank, int mpi_size);
    void exchange_ghosts(int mpi_rank, int mpi_size);

    double total_energy() const;
};

#endif
