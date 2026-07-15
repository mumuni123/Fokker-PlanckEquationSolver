#include "grid.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

int main()
{
    CylindricalVelocityGrid grid;
    grid.init(Param::momentum_umax);

    bool monotone = true;
    bool positive_widths = true;
    double symmetry_linf = 0.0;
    for (int j = 0; j < Param::Nv; ++j) {
        monotone = monotone && grid.upar_faces[j + 1] > grid.upar_faces[j];
        positive_widths = positive_widths && grid.upar_widths[j] > 0.0 &&
            std::isfinite(grid.upar_widths[j]);
    }
    for (int j = 0; j <= Param::Nv; ++j)
        symmetry_linf = std::max(symmetry_linf,
            std::fabs(grid.upar_faces[j] + grid.upar_faces[Param::Nv - j]));

    double ring_sum = 0.0;
    for (int k = 0; k < Param::Nmu; ++k) {
        monotone = monotone && grid.uperp_faces[k + 1] > grid.uperp_faces[k];
        positive_widths = positive_widths && grid.uperp_widths[k] > 0.0 &&
            std::isfinite(grid.uperp_widths[k]);
        ring_sum += grid.uperp_ring_areas[k];
    }
    const double ring_target = Const::pi * Param::momentum_umax *
        Param::momentum_umax;
    const double ring_relative = std::fabs(ring_sum - ring_target) / ring_target;
    const double min_upar = *std::min_element(
        grid.upar_widths.begin(), grid.upar_widths.end());
    const double min_uperp = *std::min_element(
        grid.uperp_widths.begin(), grid.uperp_widths.end());
    const double symmetry_tolerance = 64.0 *
        std::numeric_limits<double>::epsilon() * Param::momentum_umax;
    const bool endpoints = grid.upar_faces.front() == -Param::momentum_umax &&
        grid.upar_faces.back() == Param::momentum_umax &&
        grid.upar_faces[Param::Nv / 2] == 0.0 &&
        grid.uperp_faces.front() == 0.0 &&
        grid.uperp_faces.back() == Param::momentum_umax;
    const bool pass = monotone && positive_widths && endpoints &&
        symmetry_linf <= symmetry_tolerance && ring_relative <= 1.0e-13 &&
        !grid.is_uniform();

    std::cout << std::scientific << std::setprecision(16)
              << "cylindrical_velocity_grid_test Nv=" << Param::Nv
              << " Nmu=" << Param::Nmu
              << " upar_stretch=" << Param::momentum_upar_stretch
              << " uperp_stretch=" << Param::momentum_uperp_stretch
              << " min_upar_width=" << min_upar
              << " min_uperp_width=" << min_uperp
              << " symmetry_linf=" << symmetry_linf
              << " ring_area_relative=" << ring_relative
              << " result=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
