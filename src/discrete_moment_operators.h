#ifndef DISCRETE_MOMENT_OPERATORS_H
#define DISCRETE_MOMENT_OPERATORS_H

#include "grid.h"

#include <algorithm>
#include <cmath>

namespace Stage5 {

// DiscreteMomentOperators: all energy and momentum moments use these same
// cell/face definitions.  The analytic x velocity remains the production
// transport velocity; energy-consistent cell speed is audit-only in stage 5.
inline double u_face_width(const CylindricalVelocityGrid& grid, int jf)
{
    return 0.5 * (grid.upar_widths[jf - 1] + grid.upar_widths[jf]);
}

inline double delta_energy(const CylindricalVelocityGrid& grid, int jf, int k)
{
    return grid.kinetic_energy[idx2(jf, k)] -
           grid.kinetic_energy[idx2(jf - 1, k)];
}

inline double energy_face_speed(const CylindricalVelocityGrid& grid,
                                double mass, int jf, int k)
{
    return delta_energy(grid, jf, k) /
           (mass * Const::c * u_face_width(grid, jf));
}

inline double energy_consistent_cell_speed_candidate(
    const CylindricalVelocityGrid& grid, double mass, int j, int k,
    double analytic_speed)
{
    if (j == 0 || j == Param::Nv - 1) return analytic_speed;
    return 0.5 * (energy_face_speed(grid, mass, j, k) +
                  energy_face_speed(grid, mass, j + 1, k));
}

// SharedBudgetFCT: solve min 1/2 sum_f (1-alpha_f)^2 subject to the one
// cell donor constraint sum_f b_f alpha_f <= available.  Each face alpha is
// independent, but all outgoing anti-fluxes from the same donor share the
// finite low-order mass budget.  The caller stores each alpha on its unique
// shared face, so it is applied identically to both cell updates.
inline void shared_budget_alphas(const double contribution[4],
                                 double available, double alpha[4])
{
    double total = 0.0;
    for (int f = 0; f < 4; ++f) {
        alpha[f] = 1.0;
        total += std::max(0.0, contribution[f]);
    }
    if (total <= std::max(0.0, available)) return;

    const double budget = std::max(0.0, available);
    // With four faces this KKT problem has an exact active-set solution.
    // Removing saturated faces in descending b order is equivalent to the
    // former 48-iteration bisection but avoids a hot inner-loop reduction.
    double b[4];
    bool active[4];
    for (int f = 0; f < 4; ++f) {
        b[f] = std::max(0.0, contribution[f]);
        active[f] = b[f] > 0.0;
    }
    double lambda = 0.0;
    for (;;) {
        double sum_b = 0.0;
        double sum_b2 = 0.0;
        int largest = -1;
        for (int f = 0; f < 4; ++f) {
            if (!active[f]) continue;
            sum_b += b[f];
            sum_b2 += b[f] * b[f];
            if (largest < 0 || b[f] > b[largest]) largest = f;
        }
        if (largest < 0 || sum_b2 == 0.0) break;
        lambda = std::max(0.0, (sum_b - budget) / sum_b2);
        if (lambda <= 1.0 / b[largest]) break;
        active[largest] = false;
    }
    for (int f = 0; f < 4; ++f) {
        alpha[f] = (b[f] > 0.0) ? std::max(0.0, 1.0 - lambda * b[f]) : 1.0;
    }
}

} // namespace Stage5

#endif
