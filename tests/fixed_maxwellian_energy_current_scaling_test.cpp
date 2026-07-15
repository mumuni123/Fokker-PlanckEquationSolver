#include "discrete_moment_operators.h"
#include "nonuniform_reconstruction.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

struct CaseResult {
    double field;
    double je_direct;
    double force_power;
    double current_scale;
};

CaseResult run_case(double field)
{
    CylindricalVelocityGrid grid;
    grid.init(Param::momentum_umax);
    const double dx = Param::dx;
    const double charge = -Const::qe;
    const double thermal_u = 0.08;
    std::vector<double> mass(Param::Nvmu, 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double upar = grid.upar_cells[j];
            const double uperp = grid.uperp_cells[k];
            const double fbar = std::exp(-0.5 *
                (upar * upar + uperp * uperp) / (thermal_u * thermal_u));
            mass[idx2(j, k)] = fbar * dx * grid.cell_phase_volume(j, k);
        }
    }

    double je_direct = 0.0;
    double current_scale = 0.0;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            current_scale += std::fabs(charge * grid.vx[idx2(j, k)] *
                                       mass[idx2(j, k)] / dx);
        }
    }
    for (int jf = 1; jf < Param::Nv; ++jf) {
        const int jl = jf - 1;
        const int jr = jf;
        const int jll = (jl == 0) ? jl : jl - 1;
        const int jrr = (jr + 1 == Param::Nv) ? jr : jr + 1;
        const double s_jll = (jl == 0)
            ? grid.upar_cells[jl] - grid.upar_widths[jl]
            : grid.upar_cells[jll];
        const double s_jrr = (jr + 1 == Param::Nv)
            ? grid.upar_cells[jr] + grid.upar_widths[jr]
            : grid.upar_cells[jrr];
        for (int k = 0; k < Param::Nmu; ++k) {
            const double area = grid.uperp_ring_areas[k];
            const auto fbar = [&](int j) {
                return mass[idx2(j, k)] /
                    (dx * grid.cell_phase_volume(j, k));
            };
            const NonuniformMuscl::FaceStates states =
                NonuniformMuscl::reconstruct_face(
                    fbar(jll), fbar(jl), fbar(jr), fbar(jrr),
                    s_jll, grid.upar_cells[jl], grid.upar_cells[jr], s_jrr,
                    grid.upar_faces[jf]);
            const double cu = NonuniformMuscl::upar_face_coefficient(
                NonuniformMuscl::centered_state(states), dx, area);
            je_direct += charge * Stage5::delta_energy(grid, jf, k) * cu /
                (Const::me * Const::c * dx);
        }
    }
    // The power is formed only after direct construction of J_E.  No P/E
    // division is used, so the zero-field limit is well-defined.
    CaseResult result;
    result.field = field;
    result.je_direct = je_direct;
    result.force_power = field * je_direct;
    result.current_scale = current_scale;
    return result;
}

} // namespace

int main()
{
    const double amplitudes[] = {1.0e3, 1.0, 1.0e-3, 1.0e-6};
    double worst_relative = 0.0;
    double sign_symmetry = 0.0;
    std::cout << std::scientific << std::setprecision(17)
              << "fixed_maxwellian_energy_current_scaling_test\n";
    for (size_t i = 0; i < sizeof(amplitudes) / sizeof(amplitudes[0]); ++i) {
        const CaseResult plus = run_case(amplitudes[i]);
        const CaseResult minus = run_case(-amplitudes[i]);
        const double scale = std::max(1.0e-300, plus.current_scale);
        worst_relative = std::max(worst_relative,
            std::max(std::fabs(plus.je_direct), std::fabs(minus.je_direct)) / scale);
        sign_symmetry = std::max(sign_symmetry,
            std::fabs(plus.je_direct - minus.je_direct) / scale);
        std::cout << "abs_E=" << amplitudes[i]
                  << " JE_plus=" << plus.je_direct
                  << " JE_minus=" << minus.je_direct
                  << " P_plus=" << plus.force_power
                  << " P_minus=" << minus.force_power
                  << " JE_over_current_scale=" <<
                     std::fabs(plus.je_direct) / scale << "\n";
    }
    // A stationary symmetric Maxwellian has zero current.  The tolerance is
    // quadrature/reconstruction roundoff relative to its unsigned current
    // scale, not a tolerance derived from E or from P/E.
    std::cout << "worst_relative_JE=" << worst_relative << "\n"
              << "signed_field_symmetry=" << sign_symmetry << "\n";
    return (std::isfinite(worst_relative) && std::isfinite(sign_symmetry) &&
            worst_relative <= 2.0e-12 && sign_symmetry <= 2.0e-14) ? 0 : 1;
}
