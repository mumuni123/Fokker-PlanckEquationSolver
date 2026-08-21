#ifndef OPEN_BOUNDARY_H
#define OPEN_BOUNDARY_H

#include "grid.h"
#include "species.h"

enum class BackgroundXBoundaryType {
    RESERVOIR,
    ABSORBING,
    PERIODIC
};

enum class PhysicalSide {
    LEFT,
    RIGHT
};

struct MaxwellianReservoir {
    double density;
    double temperature;
    double drift_vx;
};

struct OpenBackgroundBoundaryConfig {
    BackgroundXBoundaryType left_type;
    BackgroundXBoundaryType right_type;
    MaxwellianReservoir left_reservoir;
    MaxwellianReservoir right_reservoir;
};

// Exchanges internal ghost cells and creates either production open-boundary
// ghosts or an explicit periodic ghost wrap. The distribution is a
// cylindrical cell-integrated mass.
class OpenBackgroundBoundary {
public:
    OpenBackgroundBoundary();
    explicit OpenBackgroundBoundary(const OpenBackgroundBoundaryConfig& config);

    const OpenBackgroundBoundaryConfig& config() const { return config_; }
    // Short-stencil ghost filling: internal MPI ghosts plus physical
    // reservoir/absorbing ghosts for the current explicit-FCT path.
    void fill_ghosts(Species& electrons, const SpatialGrid& grid,
                     int mpi_rank, int mpi_size) const;

    // Remap boundary sampling (section 13.5).  Returns the incoming 1D line
    // density lambda = M / dx of velocity cell (j_upar, k_uperp) at the given
    // physical side; absorbing inflow returns 0.  Outflow characteristics must
    // never call this method: the domain reconstruction supplies them.
    double incoming_cell_average(PhysicalSide side, int j_upar, int k_uperp,
                                 double time, const Species& electrons) const;

    bool is_incoming(PhysicalSide side, double vx) const;

private:
    OpenBackgroundBoundaryConfig config_;
    // Discrete reservoir slices are precomputed once (lazily, on first use)
    // instead of recomputing exponentials per face per step.  The stored
    // values are line densities lambda = M/dx for each velocity cell.
    mutable bool reservoir_prepared_;
    mutable double prepared_dx_;
    mutable std::vector<double> left_reservoir_line_density_;
    mutable std::vector<double> right_reservoir_line_density_;
    void ensure_reservoir_prepared(const Species& electrons) const;
    void fill_physical_left(Species& electrons, const SpatialGrid& grid) const;
    void fill_physical_right(Species& electrons, const SpatialGrid& grid) const;
};

#endif
