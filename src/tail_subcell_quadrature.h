#ifndef TAIL_SUBCELL_QUADRATURE_H
#define TAIL_SUBCELL_QUADRATURE_H

#include "grid.h"

#include <vector>

// Read-only finite-volume representation audit for the bulk/tail interface.
// A Species entry is a cell-integrated mass, so a subcell node carries a
// fraction of that mass rather than a pointwise distribution value.
struct TailSubcellNode {
    double upar;
    double uperp;
    double mass_fraction;
    double kinetic_energy;
};

// Exact-in-u_perp cylindrical phase-volume fraction carried by one energy
// bin.  The remaining one-dimensional u_parallel integral is evaluated
// piecewise analytically after splitting at every branch-change point.
struct TailEnergyBinFraction {
    int bin;
    double mass_fraction;
};

struct TailSubcellSpectrum {
    std::vector<double> number;
    std::vector<double> energy;
    double input_number;
    double input_center_energy;
    double represented_number;
    double represented_energy;
    double number_residual;
    double energy_residual;
    int straddling_cell_count;
    double straddling_cell_mass;

    explicit TailSubcellSpectrum(size_t bin_count = 0)
        : number(bin_count, 0.0), energy(bin_count, 0.0),
          input_number(0.0), input_center_energy(0.0), represented_number(0.0),
          represented_energy(0.0), number_residual(0.0),
          energy_residual(0.0), straddling_cell_count(0),
          straddling_cell_mass(0.0)
    {}
};

class TailSubcellQuadrature {
public:
    // Deterministic tensor Gauss rule on one physical (u_parallel,u_perp)
    // cell.  The u_perp factor includes the cylindrical 2*pi*u_perp measure
    // through normalized positive quadrature weights.
    static std::vector<TailSubcellNode> nodes(
        const CylindricalVelocityGrid& grid, int iv, int imu,
        int points_per_dimension = 4);

    static int energy_bin(const std::vector<double>& edges, double energy);

    // Conservative cell-volume intersection with spherical energy shells.
    // Unlike point quadrature, this representation remains smooth when the
    // requested energy bins are narrower than the spacing between fixed
    // quadrature nodes.  Fractions are finite, non-negative and sum to one
    // when energy_edges cover the complete cell.
    static std::vector<TailEnergyBinFraction> energy_bin_fractions(
        const CylindricalVelocityGrid& grid, int iv, int imu,
        const std::vector<double>& energy_edges);

    static void accumulate_cell(const CylindricalVelocityGrid& grid,
                                int iv, int imu, double cell_mass,
                                const std::vector<double>& energy_edges,
                                double threshold_energy,
                                TailSubcellSpectrum& spectrum,
                                int points_per_dimension = 4);
};

#endif
