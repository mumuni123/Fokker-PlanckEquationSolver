#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

double relative_error(double residual, double lhs, double rhs)
{
    const double scale = std::max(std::numeric_limits<double>::min(),
        std::max(std::fabs(lhs), std::fabs(rhs)));
    return std::fabs(residual) / scale;
}

} // namespace

struct CaseResult {
    double fixed_ex;
    double dt;
    double max_force_cfl;
    double dke;
    double field_work;
    double boundary_energy;
    double energy_relative;
    double je_relative;
    double min_mass;
    double min_du;
    double max_du;
};

CaseResult run_case(double fixed_ex)
{
    // Section 7.1 baseline: one uniform x cell, fixed uniform E, no Beam,
    // no x transport, PPM, CTU, FCT, or Ampere feedback.  The update below
    // intentionally duplicates the production donor-cell u flux contract:
    // C_u = M_donor / du_donor, Phi_u = a_x * C_u.
    CylindricalVelocityGrid grid;
    grid.init(Param::momentum_umax);

    const double dx = Param::dx;
    const double charge = -Const::qe;
    const double mass = Const::me;
    const double acceleration = charge * fixed_ex / (mass * Const::c);

    double min_du = grid.upar_widths.front();
    for (size_t j = 1; j < grid.upar_widths.size(); ++j)
        min_du = std::min(min_du, grid.upar_widths[j]);
    const double dt = 0.20 * min_du / std::fabs(acceleration);

    std::vector<double> mass_old(Param::Nvmu, 0.0);
    std::vector<double> mass_new(Param::Nvmu, 0.0);
    std::vector<double> cu(static_cast<size_t>(Param::Nv + 1) * Param::Nmu,
                           0.0);
    std::vector<double> fu(static_cast<size_t>(Param::Nv + 1) * Param::Nmu,
                           0.0);
    const double uth = 0.08;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double up = grid.upar_cells[j];
            const double ut = grid.uperp_cells[k];
            const double f3 = std::exp(-0.5 * (up * up + ut * ut) /
                                       (uth * uth));
            mass_old[idx2(j, k)] =
                f3 * dx * grid.cell_phase_volume(j, k);
        }
    }

    for (int jf = 1; jf < Param::Nv; ++jf) {
        const int donor = (acceleration >= 0.0) ? jf - 1 : jf;
        const double donor_du = grid.upar_widths[donor];
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t face = static_cast<size_t>(jf) * Param::Nmu + k;
            cu[face] = Stage5::donor_cell_coefficient(
                mass_old[idx2(donor, k)], donor_du);
            fu[face] = acceleration * cu[face];
        }
    }
    if (acceleration < 0.0) {
        for (int k = 0; k < Param::Nmu; ++k) {
            cu[k] = Stage5::donor_cell_coefficient(
                mass_old[idx2(0, k)], grid.upar_widths[0]);
            fu[k] = acceleration * cu[k];
        }
    } else if (acceleration > 0.0) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t face = static_cast<size_t>(Param::Nv) * Param::Nmu + k;
            cu[face] = Stage5::donor_cell_coefficient(
                mass_old[idx2(Param::Nv - 1, k)],
                grid.upar_widths[Param::Nv - 1]);
            fu[face] = acceleration * cu[face];
        }
    }

    double dke = 0.0;
    double je_from_delta_energy = 0.0;
    double je_from_energy_speed = 0.0;
    double min_mass = std::numeric_limits<double>::infinity();
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t cell = idx2(j, k);
            const size_t lower = static_cast<size_t>(j) * Param::Nmu + k;
            const size_t upper = static_cast<size_t>(j + 1) * Param::Nmu + k;
            mass_new[cell] = mass_old[cell] - dt * (fu[upper] - fu[lower]);
            min_mass = std::min(min_mass, mass_new[cell]);
            dke += grid.kinetic_energy[cell] *
                (mass_new[cell] - mass_old[cell]);
        }
    }
    for (int jf = 1; jf < Param::Nv; ++jf) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t face = static_cast<size_t>(jf) * Param::Nmu + k;
            const double delta_k = Stage5::delta_energy(grid, jf, k);
            const double center_distance = Stage5::upar_center_distance(grid, jf);
            je_from_delta_energy +=
                charge * delta_k * cu[face] / (mass * Const::c * dx);
            const double v_energy =
                delta_k / (mass * Const::c * center_distance);
            je_from_energy_speed +=
                charge * v_energy * center_distance * cu[face] / dx;
        }
    }

    double boundary_energy = 0.0;
    for (int k = 0; k < Param::Nmu; ++k) {
        const double f_lo = fu[static_cast<size_t>(k)];
        const double f_hi =
            fu[static_cast<size_t>(Param::Nv) * Param::Nmu + k];
        boundary_energy += dt * (
            grid.kinetic_energy[idx2(0, k)] * f_lo -
            grid.kinetic_energy[idx2(Param::Nv - 1, k)] * f_hi);
    }

    const double field_work = dt * fixed_ex * je_from_delta_energy * dx;
    const double rhs = field_work + boundary_energy;
    const double residual = dke - rhs;
    const double relative = relative_error(residual, dke, rhs);
    const double je_difference =
        je_from_delta_energy - je_from_energy_speed;
    const double je_relative = relative_error(
        je_difference, je_from_delta_energy, je_from_energy_speed);

    double max_force_cfl = 0.0;
    double max_du = 0.0;
    for (size_t j = 0; j < grid.upar_widths.size(); ++j) {
        max_force_cfl = std::max(max_force_cfl,
            std::fabs(acceleration) * dt / grid.upar_widths[j]);
        max_du = std::max(max_du, grid.upar_widths[j]);
    }
    CaseResult result;
    result.fixed_ex = fixed_ex;
    result.dt = dt;
    result.max_force_cfl = max_force_cfl;
    result.dke = dke;
    result.field_work = field_work;
    result.boundary_energy = boundary_energy;
    result.energy_relative = relative;
    result.je_relative = je_relative;
    result.min_mass = min_mass;
    result.min_du = min_du;
    result.max_du = max_du;
    return result;
}

void print_case(const CaseResult& result)
{
    std::cout << "E_fixed=" << result.fixed_ex
              << " dt=" << result.dt
              << " max_force_cfl=" << result.max_force_cfl
              << " dK=" << result.dke
              << " dt_E_JE_dx=" << result.field_work
              << " B_u_K=" << result.boundary_energy
              << " energy_relative=" << result.energy_relative
              << " JE_relative_difference=" << result.je_relative
              << " min_mass_new=" << result.min_mass << "\n";
}

bool passes_case(const CaseResult& result)
{
    return std::isfinite(result.energy_relative) &&
        std::isfinite(result.je_relative) &&
        result.max_force_cfl <= 0.2 + 64.0 * std::numeric_limits<double>::epsilon() &&
        result.min_mass >= -1.0e-14 && result.energy_relative <= 1.0e-11 &&
        result.je_relative <= 1.0e-14;
}

int main()
{
    const CaseResult positive = run_case(1.0e10);
    const CaseResult negative = run_case(-1.0e10);
    const double min_du = std::min(positive.min_du, negative.min_du);
    const double max_du = std::max(positive.max_du, negative.max_du);
    std::cout << std::scientific << std::setprecision(16)
              << "fixed_field_low_order_test\n"
              << "upar_stretch=" << Param::momentum_upar_stretch
              << " min_du=" << min_du << " max_du=" << max_du
              << " du_ratio=" << max_du / min_du << "\n";
    print_case(positive);
    print_case(negative);
    return (passes_case(positive) && passes_case(negative)) ? 0 : 1;
}
