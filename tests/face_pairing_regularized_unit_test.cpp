#include "regularized_face_pairing.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace {

std::vector<double> mode(int n, int wave_number)
{
    std::vector<double> value(static_cast<size_t>(n), 0.0);
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int i = 0; i < n; ++i)
        value[static_cast<size_t>(i)] = std::cos(
            two_pi * wave_number * i / static_cast<double>(n));
    return value;
}

double residual_norm(const std::vector<double>& residual,
                     const std::vector<double>& correction)
{
    std::vector<double> image;
    RegularizedFacePairing::apply_gstar(correction, image);
    for (size_t i = 0; i < image.size(); ++i) image[i] -= residual[i];
    return RegularizedFacePairing::l2_norm(image);
}

bool adjoint_case(int n)
{
    std::vector<double> face(static_cast<size_t>(n), 0.0);
    std::vector<double> cell(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        face[static_cast<size_t>(i)] =
            std::sin(0.37 * i) + 0.2 * std::cos(0.11 * i);
        cell[static_cast<size_t>(i)] =
            std::cos(0.23 * i) - 0.15 * std::sin(0.41 * i);
    }
    std::vector<double> g_face;
    std::vector<double> gstar_cell;
    RegularizedFacePairing::apply_g(face, g_face);
    RegularizedFacePairing::apply_gstar(cell, gstar_cell);
    long double left = 0.0L;
    long double right = 0.0L;
    for (int i = 0; i < n; ++i) {
        left += static_cast<long double>(g_face[static_cast<size_t>(i)]) *
                cell[static_cast<size_t>(i)];
        right += static_cast<long double>(face[static_cast<size_t>(i)]) *
                 gstar_cell[static_cast<size_t>(i)];
    }
    const double error = std::fabs(static_cast<double>(left - right));
    const double scale = std::max(
        1.0, std::max(std::fabs(static_cast<double>(left)),
                      std::fabs(static_cast<double>(right))));
    const bool pass = error <=
        256.0 * std::numeric_limits<double>::epsilon() * scale;
    std::printf("adjoint error=%.16e scale=%.16e PASS=%d\n",
                error, scale, pass ? 1 : 0);
    return pass;
}

bool solve_case(const char* name, const std::vector<double>& residual,
                double sigma_cutoff, double bound, bool expect_improvement,
                bool expect_unresolved)
{
    const size_t n = residual.size();
    std::vector<double> weights(n, 1.0);
    std::vector<double> lower(n, -bound);
    std::vector<double> upper(n, bound);
    RegularizedFacePairing::Config config;
    config.sigma_cutoff = sigma_cutoff;
    config.lambda = 1.0e-5;
    config.eta = 1.0e-10;
    config.max_iterations = 1000;
    config.relative_tolerance = 1.0e-11;
    RegularizedFacePairing::Diagnostics diagnostics;
    std::vector<double> correction;
    std::vector<double> unresolved;
    const bool solved = RegularizedFacePairing::solve(
        residual, weights, lower, upper, config, correction, unresolved,
        diagnostics);
    const double before = RegularizedFacePairing::l2_norm(residual);
    const double after = residual_norm(residual, correction);
    const bool improved = after < before;
    const bool unresolved_ok = expect_unresolved
        ? diagnostics.unresolved_mode_l2 > 0.9 * before
        : diagnostics.unresolved_mode_l2 <
              1.0e-10 * std::max(1.0, before);
    const double objective_sum =
        diagnostics.objective_residual +
        diagnostics.objective_smoothness +
        diagnostics.objective_amplitude;
    const bool objective_ok =
        std::isfinite(diagnostics.objective_total) &&
        diagnostics.objective_residual >= 0.0 &&
        diagnostics.objective_smoothness >= 0.0 &&
        diagnostics.objective_amplitude >= 0.0 &&
        std::fabs(diagnostics.objective_total - objective_sum) <=
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::fabs(objective_sum));
    const bool pass = solved && improved == expect_improvement &&
        unresolved_ok && diagnostics.correction_linf <=
            bound * (1.0 + 1.0e-12) && objective_ok &&
        diagnostics.nonzero_capacity_cells == static_cast<long long>(n) &&
        diagnostics.bound_saturated_cells ==
            diagnostics.trust_region_active_cells;
    std::printf(
        "%s solved=%d before=%.16e after=%.16e unresolved=%.16e "
        "objective=%.16e iterations=%d bound_active=%lld PASS=%d\n",
        name, solved ? 1 : 0, before, after,
        diagnostics.unresolved_mode_l2, diagnostics.objective_total,
        diagnostics.iterations,
        diagnostics.trust_region_active_cells, pass ? 1 : 0);
    return pass;
}

bool fallback_case(int n)
{
    const std::vector<double> residual = mode(n, n / 2);
    std::vector<double> weights(static_cast<size_t>(n), 1.0);
    std::vector<double> lower(static_cast<size_t>(n), -2.0);
    std::vector<double> upper(static_cast<size_t>(n), 2.0);
    const std::vector<double> baseline(static_cast<size_t>(n), 0.0);
    RegularizedFacePairing::Config config;
    config.sigma_cutoff = 1.0e-8;
    RegularizedFacePairing::Diagnostics diagnostics;
    std::vector<double> correction;
    std::vector<double> unresolved;
    const bool solved = RegularizedFacePairing::solve(
        residual, weights, lower, upper, config, correction, unresolved,
        diagnostics);
    const double before = residual_norm(residual, baseline);
    const double after = residual_norm(residual, correction);
    const bool accepted = solved && after < before;
    if (!accepted) correction = baseline;
    double restore_error = 0.0;
    for (size_t i = 0; i < baseline.size(); ++i)
        restore_error = std::max(
            restore_error, std::fabs(correction[i] - baseline[i]));
    const bool pass = solved && !accepted && restore_error == 0.0;
    std::printf(
        "monotone_fallback solved=%d accepted=%d before=%.16e "
        "after=%.16e restore_error=%.16e PASS=%d\n",
        solved ? 1 : 0, accepted ? 1 : 0, before, after, restore_error,
        pass ? 1 : 0);
    return pass;
}

}

int main()
{
    const int n = 64;
    bool pass = adjoint_case(n);
    pass = solve_case("constant", std::vector<double>(n, 1.0),
                      1.0e-8, 2.0, true, false) && pass;
    pass = solve_case("low_frequency", mode(n, 3), 1.0e-8, 2.0,
                      true, false) && pass;
    pass = solve_case("near_nyquist", mode(n, n / 2 - 1), 0.1, 2.0,
                      false, true) && pass;
    pass = solve_case("alternating_cutoff_1e-6", mode(n, n / 2), 1.0e-6,
                      2.0, false, true) && pass;
    pass = solve_case("alternating_cutoff_1e-8", mode(n, n / 2), 1.0e-8,
                      2.0, false, true) && pass;
    pass = solve_case("alternating_cutoff_1e-10", mode(n, n / 2), 1.0e-10,
                      2.0, false, true) && pass;
    pass = solve_case("capacity_saturated", mode(n, 2), 1.0e-8, 1.0e-3,
                      true, false) && pass;
    pass = fallback_case(n) && pass;
    std::printf("face_pairing_regularized_unit_test %s\n",
                pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
