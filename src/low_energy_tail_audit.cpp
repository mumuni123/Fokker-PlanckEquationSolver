#include "low_energy_tail_audit.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

double kinetic_energy(const BackgroundTailParticle& particle)
{
    const double u2 = particle.ux * particle.ux +
                      particle.uy * particle.uy +
                      particle.uz * particle.uz;
    return (std::sqrt(1.0 + u2) - 1.0) *
           Const::me * Const::c * Const::c * particle.weight;
}

double kinetic_energy_mev(const BackgroundTailParticle& particle)
{
    return kinetic_energy(particle) /
           std::max(particle.weight, std::numeric_limits<double>::min()) /
           (1.0e6 * Const::eV);
}

void add_shape(int cell, double shape, double number, double energy,
               LowEnergyTailThresholdStats& stats, int nx)
{
    if (!(shape > 0.0)) return;
    if (cell >= 0 && cell < nx) {
        const size_t i = static_cast<size_t>(cell);
        stats.cell_number[i] += shape * number;
        stats.cell_kinetic_energy[i] += shape * energy;
        stats.cell_macro_supports[i] += shape;
    } else {
        stats.outside_shape_number += shape * number;
        stats.outside_shape_energy += shape * energy;
    }
}

} // namespace

LowEnergyTailAuditResult audit_low_energy_tail_local(
    const std::vector<BackgroundTailParticle>& particles,
    const SpatialGrid& grid, const std::vector<double>& thresholds_mev)
{
    LowEnergyTailAuditResult result;
    result.thresholds.resize(thresholds_mev.size());
    for (size_t k = 0; k < thresholds_mev.size(); ++k) {
        LowEnergyTailThresholdStats& stats = result.thresholds[k];
        stats.threshold_mev = thresholds_mev[k];
        stats.cell_number.assign(static_cast<size_t>(grid.nx_global), 0.0);
        stats.cell_kinetic_energy.assign(
            static_cast<size_t>(grid.nx_global), 0.0);
        stats.cell_macro_supports.assign(
            static_cast<size_t>(grid.nx_global), 0.0);
        if (!std::isfinite(stats.threshold_mev) ||
            !(stats.threshold_mev > 0.0)) result.finite = false;
    }

    for (size_t p = 0; p < particles.size(); ++p) {
        const BackgroundTailParticle& particle = particles[p];
        const double energy = kinetic_energy(particle);
        const double energy_mev = kinetic_energy_mev(particle);
        if (!std::isfinite(particle.x) || !std::isfinite(particle.ux) ||
            !std::isfinite(particle.uy) || !std::isfinite(particle.uz) ||
            !std::isfinite(particle.weight) || !(particle.weight >= 0.0) ||
            !std::isfinite(energy) || !std::isfinite(energy_mev)) {
            result.finite = false;
            continue;
        }
        result.total_number += particle.weight;
        result.total_kinetic_energy += energy;
        ++result.total_macro_particles;

        const double coordinate =
            (particle.x - grid.x_min) / grid.dx - 0.5;
        const int cell0 = static_cast<int>(std::floor(coordinate));
        const double fraction = coordinate - static_cast<double>(cell0);
        for (size_t k = 0; k < result.thresholds.size(); ++k) {
            LowEnergyTailThresholdStats& stats = result.thresholds[k];
            if (energy_mev > stats.threshold_mev) continue;
            stats.number += particle.weight;
            stats.kinetic_energy += energy;
            ++stats.macro_particles;
            add_shape(cell0, 1.0 - fraction, particle.weight, energy,
                      stats, grid.nx_global);
            add_shape(cell0 + 1, fraction, particle.weight, energy,
                      stats, grid.nx_global);
        }
    }
    return result;
}
