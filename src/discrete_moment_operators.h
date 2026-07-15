#ifndef DISCRETE_MOMENT_OPERATORS_H
#define DISCRETE_MOMENT_OPERATORS_H

#include "grid.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Stage5 {

// DiscreteMomentOperators: all energy and momentum moments use these same
// cell/face definitions.  The analytic x velocity remains the production
// transport velocity; energy-consistent cell speed is audit-only in stage 5.
//
// Cylindrical u-face coefficient contract:
//   M[j,k] = f_3 * dx * du_parallel[j] * (2*pi*u_perp*du_perp)[k]
//   C_u[j+1/2,k] = M_donor / du_parallel_donor
//   Phi_u = a_x * C_u,  a_x = q E_x / (m c).
//
// Thus C_u carries dx and the perpendicular ring measure implicitly through
// M, but it contains neither a_x nor an additional 1/dx factor.  The energy
// current is J_E = q/(m c dx) sum(delta_K * C_u).  A substep accumulates
// h*J_E and divides by dt only after all substeps are complete.
enum CylindricalUFluxContract {
    CU_DIVIDES_BY_UPAR_DONOR_WIDTH = 1,
    CU_CONTAINS_DX_VIA_CELL_MASS = 1,
    CU_CONTAINS_UPERP_RING_VIA_CELL_MASS = 1,
    CU_CONTAINS_ACCELERATION = 0,
    JE_REQUIRES_CELL_DX_DIVISION = 1,
    JE_IS_SUBSTEP_TIME_AVERAGED = 1
};

// The donor width is selected by the caller from the upwind cell.  Keeping
// this primitive here prevents production and standalone audits from quietly
// diverging on nonuniform velocity grids.
inline double donor_cell_coefficient(double donor_mass, double donor_width)
{
    return donor_mass / donor_width;
}

inline double upar_center_distance(const CylindricalVelocityGrid& grid, int jf)
{
    return grid.upar_cells[jf] - grid.upar_cells[jf - 1];
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
           (mass * Const::c * upar_center_distance(grid, jf));
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
    // The high-order candidate remains the default.  Limit only when its
    // antidiffusive outflow exceeds the positive low-order mass budget by
    // more than the local floating-point summation envelope.
    const double budget_with_roundoff = std::max(0.0, available) +
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(total, std::fabs(available)));
    if (total <= budget_with_roundoff) return;

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
