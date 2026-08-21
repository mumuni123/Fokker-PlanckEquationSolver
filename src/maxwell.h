#ifndef MAXWELL_H
#define MAXWELL_H

#include "grid.h"
#include <vector>

// Shared field storage for the open electrostatic solver.  Ex_face contains
// n_local+1 physical faces; the two global endpoints are independent.
struct EMFields {
    std::vector<double> Ex;
    std::vector<double> Ex_face;
    std::vector<double> phi;
    std::vector<double> rho;
    int nx_total;
    double dx;
    double last_gauss_residual_l1;
    double last_gauss_residual_linf;

    void init(const SpatialGrid& grid);
    void set_charge_density(const class Species& electrons,
                            const std::vector<double>& tail_density,
                            const std::vector<double>& beam_density,
                            const std::vector<double>& ion_density);
};

#endif
