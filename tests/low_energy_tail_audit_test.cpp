#include "low_energy_tail_audit.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

double ux_for_mev(double mev)
{
    const double gamma = 1.0 + mev * 1.0e6 * Const::eV /
                                  (Const::me * Const::c * Const::c);
    return std::sqrt(gamma * gamma - 1.0);
}

BackgroundTailParticle particle(double x, double mev, double weight,
                                std::uint64_t id)
{
    BackgroundTailParticle p = {};
    p.x = x;
    p.ux = ux_for_mev(mev);
    p.weight = weight;
    p.id = id;
    return p;
}

} // namespace

int main()
{
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 4, 4.0);
    std::vector<BackgroundTailParticle> particles;
    particles.push_back(particle(0.5, 5.0, 2.0, 1));
    particles.push_back(particle(1.5, 5.6, 3.0, 2));
    particles.push_back(particle(2.5, 5.9, 5.0, 3));
    particles.push_back(particle(3.5, 6.2, 7.0, 4));
    const std::vector<BackgroundTailParticle> before = particles;
    std::vector<double> thresholds;
    thresholds.push_back(5.5);
    thresholds.push_back(5.75);
    thresholds.push_back(6.0);
    const LowEnergyTailAuditResult result =
        audit_low_energy_tail_local(particles, grid, thresholds);
    const bool unchanged = particles.size() == before.size() &&
        particles[0].x == before[0].x && particles[3].ux == before[3].ux;
    const bool pass = result.finite && unchanged &&
        result.total_macro_particles == 4 &&
        result.thresholds[0].macro_particles == 1 &&
        result.thresholds[1].macro_particles == 2 &&
        result.thresholds[2].macro_particles == 3 &&
        result.thresholds[0].number == 2.0 &&
        result.thresholds[1].number == 5.0 &&
        result.thresholds[2].number == 10.0 &&
        std::fabs(result.thresholds[2].cell_number[0] - 2.0) < 1.0e-14 &&
        std::fabs(result.thresholds[2].cell_number[1] - 3.0) < 1.0e-14 &&
        std::fabs(result.thresholds[2].cell_number[2] - 5.0) < 1.0e-14 &&
        result.thresholds[2].outside_shape_number == 0.0;
    std::cout << "finite=" << result.finite << "\n"
              << "read_only=" << unchanged << "\n"
              << "n_5p5=" << result.thresholds[0].number << "\n"
              << "n_5p75=" << result.thresholds[1].number << "\n"
              << "n_6p0=" << result.thresholds[2].number << "\n"
              << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
