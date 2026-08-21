#ifndef LOW_ENERGY_TAIL_AUDIT_H
#define LOW_ENERGY_TAIL_AUDIT_H

#include "background_tail_pic.h"
#include "grid.h"

#include <cstdint>
#include <vector>

struct LowEnergyTailThresholdStats {
    double threshold_mev;
    double number;
    double kinetic_energy;
    std::uint64_t macro_particles;
    double outside_shape_number;
    double outside_shape_energy;
    std::vector<double> cell_number;
    std::vector<double> cell_kinetic_energy;
    std::vector<double> cell_macro_supports;

    LowEnergyTailThresholdStats()
        : threshold_mev(0.0), number(0.0), kinetic_energy(0.0),
          macro_particles(0), outside_shape_number(0.0),
          outside_shape_energy(0.0)
    {}
};

struct LowEnergyTailAuditResult {
    bool finite;
    double total_number;
    double total_kinetic_energy;
    std::uint64_t total_macro_particles;
    std::vector<LowEnergyTailThresholdStats> thresholds;

    LowEnergyTailAuditResult()
        : finite(true), total_number(0.0), total_kinetic_energy(0.0),
          total_macro_particles(0)
    {}
};

// Read-only local contribution. Cell arrays use global cell indices so they
// can be summed directly with MPI_Allreduce by the checkpoint audit driver.
LowEnergyTailAuditResult audit_low_energy_tail_local(
    const std::vector<BackgroundTailParticle>& particles,
    const SpatialGrid& grid, const std::vector<double>& thresholds_mev);

#endif
