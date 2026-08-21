#include "maxwell.h"
#include "species.h"

#include <algorithm>

void EMFields::init(const SpatialGrid& grid)
{
    nx_total = grid.nx_total;
    dx = grid.dx;
    Ex.assign(static_cast<size_t>(grid.nx_total), 0.0);
    Ex_face.assign(static_cast<size_t>(grid.nx_local + 1), 0.0);
    phi.assign(static_cast<size_t>(grid.nx_total), 0.0);
    rho.assign(static_cast<size_t>(grid.nx_total), 0.0);
    last_gauss_residual_l1 = 0.0;
    last_gauss_residual_linf = 0.0;
}

void EMFields::set_charge_density(const Species& electrons,
                                  const std::vector<double>& tail_density,
                                  const std::vector<double>& beam_density,
                                  const std::vector<double>& ion_density)
{
    const int ng = electrons.sgrid->nghost;
    const int nxl = electrons.sgrid->nx_local;
    std::fill(rho.begin(), rho.end(), 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        const double ni = ix < static_cast<int>(ion_density.size()) ? ion_density[ix] : 0.0;
        const double nt = ix < static_cast<int>(tail_density.size()) ? tail_density[ix] : 0.0;
        const double nb = ix < static_cast<int>(beam_density.size()) ? beam_density[ix] : 0.0;
        rho[ng + ix] = Const::qe * (ni - electrons.number_density[ix] - nt - nb);
    }
}
