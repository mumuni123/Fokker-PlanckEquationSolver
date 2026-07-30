#ifndef REGULARIZED_FACE_PAIRING_H
#define REGULARIZED_FACE_PAIRING_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace RegularizedFacePairing {

struct Config {
    double sigma_cutoff;
    double lambda;
    double eta;
    double trust_fraction;
    int max_iterations;
    double relative_tolerance;

    Config()
        : sigma_cutoff(1.0e-8), lambda(1.0e-3), eta(1.0e-8),
          trust_fraction(0.1), max_iterations(400),
          relative_tolerance(1.0e-10) {}
};

struct Diagnostics {
    int valid;
    int converged;
    int iterations;
    int unresolved_mode_count;
    long long capacity_active_cells;
    long long trust_region_active_cells;
    long long nonzero_capacity_cells;
    long long bound_saturated_cells;
    double unresolved_mode_l2;
    double representative_mode_l2;
    double projected_gradient_linf;
    double correction_l2;
    double correction_linf;
    double objective_residual;
    double objective_smoothness;
    double objective_amplitude;
    double objective_total;

    Diagnostics()
        : valid(1), converged(0), iterations(0),
          unresolved_mode_count(0), capacity_active_cells(0),
          trust_region_active_cells(0), nonzero_capacity_cells(0),
          bound_saturated_cells(0), unresolved_mode_l2(0.0),
          representative_mode_l2(0.0), projected_gradient_linf(0.0),
          correction_l2(0.0), correction_linf(0.0),
          objective_residual(0.0), objective_smoothness(0.0),
          objective_amplitude(0.0), objective_total(0.0) {}
};

inline void apply_gstar(const std::vector<double>& cell,
                        std::vector<double>& face)
{
    const size_t n = cell.size();
    face.assign(n, 0.0);
    if (n == 0) return;
    for (size_t i = 0; i < n; ++i)
        face[i] = 0.5 * (cell[(i + n - 1) % n] + cell[i]);
}

inline void apply_g(const std::vector<double>& face,
                    std::vector<double>& cell)
{
    const size_t n = face.size();
    cell.assign(n, 0.0);
    if (n == 0) return;
    for (size_t i = 0; i < n; ++i)
        cell[i] = 0.5 * (face[i] + face[(i + 1) % n]);
}

inline double l2_norm(const std::vector<double>& values,
                      const std::vector<double>* weights = 0)
{
    long double sum = 0.0L;
    for (size_t i = 0; i < values.size(); ++i) {
        const double weight = weights ? (*weights)[i] : 1.0;
        sum += static_cast<long double>(weight) * values[i] * values[i];
    }
    return std::sqrt(static_cast<double>(std::max(0.0L, sum)));
}

// The periodic arithmetic-average G* has singular values
// |cos(pi*k/N)|.  Only modes below the requested relative cutoff are
// projected out.  This is O(N*m_unresolved), not a dense pseudo-inverse.
inline bool split_representable_residual(
    const std::vector<double>& residual, double sigma_cutoff,
    std::vector<double>& representable, std::vector<double>& unresolved,
    Diagnostics& diagnostics)
{
    const size_t n = residual.size();
    representable = residual;
    unresolved.assign(n, 0.0);
    if (n == 0 || !(sigma_cutoff >= 0.0) ||
        !std::isfinite(sigma_cutoff)) {
        diagnostics.valid = 0;
        return false;
    }

    const double two_pi = 2.0 * std::acos(-1.0);
    for (size_t k = 0; k < n; ++k) {
        const double sigma = std::fabs(
            std::cos(std::acos(-1.0) * static_cast<double>(k) /
                     static_cast<double>(n)));
        if (!(sigma < sigma_cutoff)) continue;
        std::complex<long double> coefficient(0.0L, 0.0L);
        for (size_t i = 0; i < n; ++i) {
            const long double angle = -static_cast<long double>(two_pi) *
                static_cast<long double>(k * i) /
                static_cast<long double>(n);
            coefficient += static_cast<long double>(residual[i]) *
                std::complex<long double>(std::cos(angle), std::sin(angle));
        }
        coefficient /= static_cast<long double>(n);
        for (size_t i = 0; i < n; ++i) {
            const long double angle = static_cast<long double>(two_pi) *
                static_cast<long double>(k * i) /
                static_cast<long double>(n);
            unresolved[i] += static_cast<double>(std::real(
                coefficient *
                std::complex<long double>(std::cos(angle), std::sin(angle))));
        }
        ++diagnostics.unresolved_mode_count;
    }
    for (size_t i = 0; i < n; ++i)
        representable[i] -= unresolved[i];
    diagnostics.unresolved_mode_l2 = l2_norm(unresolved);
    diagnostics.representative_mode_l2 = l2_norm(representable);
    return true;
}

inline void apply_hessian(const std::vector<double>& value,
                          const std::vector<double>& weights,
                          double lambda, double eta,
                          std::vector<double>& output)
{
    std::vector<double> face;
    std::vector<double> weighted_face;
    apply_gstar(value, face);
    weighted_face.resize(face.size());
    for (size_t i = 0; i < face.size(); ++i)
        weighted_face[i] = weights[i] * face[i];
    apply_g(weighted_face, output);
    const size_t n = value.size();
    for (size_t i = 0; i < n; ++i) {
        const double laplacian =
            2.0 * value[i] - value[(i + n - 1) % n] -
            value[(i + 1) % n];
        output[i] += lambda * laplacian + eta * value[i];
    }
}

inline double projected_gradient_linf(
    const std::vector<double>& value, const std::vector<double>& gradient,
    const std::vector<double>& lower, const std::vector<double>& upper)
{
    double result = 0.0;
    for (size_t i = 0; i < value.size(); ++i) {
        double component = gradient[i];
        if ((value[i] <= lower[i] && component > 0.0) ||
            (value[i] >= upper[i] && component < 0.0))
            component = 0.0;
        result = std::max(result, std::fabs(component));
    }
    return result;
}

// Projected conjugate-gradient with restart at active-set changes.  The
// matrix is applied through G/G* and periodic stencils; no dense inverse is
// formed.  Bounds already include the requested trust fraction.
inline bool solve(const std::vector<double>& face_residual,
                  const std::vector<double>& face_weights,
                  const std::vector<double>& lower_bound,
                  const std::vector<double>& upper_bound,
                  const Config& config,
                  std::vector<double>& correction,
                  std::vector<double>& unresolved,
                  Diagnostics& diagnostics)
{
    const size_t n = face_residual.size();
    if (n == 0 || face_weights.size() != n ||
        lower_bound.size() != n || upper_bound.size() != n ||
        !(config.lambda >= 0.0) || !(config.eta > 0.0) ||
        !(config.trust_fraction > 0.0 &&
          config.trust_fraction <= 1.0) ||
        config.max_iterations < 1 ||
        !(config.relative_tolerance > 0.0)) {
        diagnostics.valid = 0;
        return false;
    }
    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(face_residual[i]) ||
            !std::isfinite(face_weights[i]) || face_weights[i] < 0.0 ||
            !std::isfinite(lower_bound[i]) ||
            !std::isfinite(upper_bound[i]) ||
            lower_bound[i] > upper_bound[i]) {
            diagnostics.valid = 0;
            return false;
        }
        if (lower_bound[i] < 0.0 || upper_bound[i] > 0.0)
            ++diagnostics.nonzero_capacity_cells;
    }

    std::vector<double> representative;
    if (!split_representable_residual(
            face_residual, config.sigma_cutoff, representative,
            unresolved, diagnostics))
        return false;

    std::vector<double> weighted_residual(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        weighted_residual[i] = face_weights[i] * representative[i];
    std::vector<double> rhs;
    apply_g(weighted_residual, rhs);

    correction.assign(n, 0.0);
    std::vector<double> hessian(n), gradient(n), direction(n), hdirection(n);
    apply_hessian(correction, face_weights, config.lambda, config.eta,
                  hessian);
    for (size_t i = 0; i < n; ++i) {
        gradient[i] = hessian[i] - rhs[i];
        direction[i] = -gradient[i];
    }
    const double rhs_scale = std::max(1.0, l2_norm(rhs));
    double old_projected_square = 0.0;
    for (size_t i = 0; i < n; ++i)
        old_projected_square += direction[i] * direction[i];

    for (int iteration = 0; iteration < config.max_iterations; ++iteration) {
        diagnostics.iterations = iteration + 1;
        diagnostics.projected_gradient_linf = projected_gradient_linf(
            correction, gradient, lower_bound, upper_bound);
        if (diagnostics.projected_gradient_linf <=
            config.relative_tolerance * rhs_scale) {
            diagnostics.converged = 1;
            break;
        }

        apply_hessian(direction, face_weights, config.lambda, config.eta,
                      hdirection);
        long double denominator = 0.0L;
        long double numerator = 0.0L;
        for (size_t i = 0; i < n; ++i) {
            denominator += static_cast<long double>(direction[i]) *
                           hdirection[i];
            numerator -= static_cast<long double>(gradient[i]) * direction[i];
        }
        if (!(denominator > 0.0L) || !std::isfinite(
                static_cast<double>(denominator))) {
            diagnostics.valid = 0;
            return false;
        }
        const double alpha = std::max(
            0.0, static_cast<double>(numerator / denominator));
        bool active_set_changed = false;
        for (size_t i = 0; i < n; ++i) {
            const double trial = correction[i] + alpha * direction[i];
            const double projected = std::max(
                lower_bound[i], std::min(upper_bound[i], trial));
            active_set_changed = active_set_changed || projected != trial;
            correction[i] = projected;
        }

        apply_hessian(correction, face_weights, config.lambda, config.eta,
                      hessian);
        std::vector<double> next_direction(n, 0.0);
        double projected_square = 0.0;
        for (size_t i = 0; i < n; ++i) {
            gradient[i] = hessian[i] - rhs[i];
            double pg = -gradient[i];
            if ((correction[i] <= lower_bound[i] && pg < 0.0) ||
                (correction[i] >= upper_bound[i] && pg > 0.0))
                pg = 0.0;
            next_direction[i] = pg;
            projected_square += pg * pg;
        }
        const double beta = active_set_changed || old_projected_square <= 0.0
            ? 0.0 : projected_square / old_projected_square;
        for (size_t i = 0; i < n; ++i)
            direction[i] = next_direction[i] + beta * direction[i];
        old_projected_square = projected_square;
    }

    diagnostics.correction_l2 = l2_norm(correction);
    for (size_t i = 0; i < n; ++i) {
        diagnostics.correction_linf = std::max(
            diagnostics.correction_linf, std::fabs(correction[i]));
        const double bound_scale = std::max(
            std::fabs(lower_bound[i]), std::fabs(upper_bound[i]));
        if (bound_scale > 0.0 &&
            (std::fabs(correction[i] - lower_bound[i]) <=
                 64.0 * std::numeric_limits<double>::epsilon() * bound_scale ||
             std::fabs(correction[i] - upper_bound[i]) <=
                 64.0 * std::numeric_limits<double>::epsilon() * bound_scale))
            ++diagnostics.trust_region_active_cells;
    }
    diagnostics.capacity_active_cells =
        diagnostics.nonzero_capacity_cells;
    diagnostics.bound_saturated_cells =
        diagnostics.trust_region_active_cells;

    std::vector<double> fitted_face;
    apply_gstar(correction, fitted_face);
    long double residual_objective = 0.0L;
    long double smoothness_objective = 0.0L;
    long double amplitude_objective = 0.0L;
    for (size_t i = 0; i < n; ++i) {
        const long double mismatch =
            static_cast<long double>(fitted_face[i]) -
            static_cast<long double>(representative[i]);
        residual_objective +=
            0.5L * static_cast<long double>(face_weights[i]) *
            mismatch * mismatch;
        const long double difference =
            static_cast<long double>(correction[i]) -
            static_cast<long double>(correction[(i + n - 1) % n]);
        smoothness_objective +=
            0.5L * static_cast<long double>(config.lambda) *
            difference * difference;
        const long double value = correction[i];
        amplitude_objective +=
            0.5L * static_cast<long double>(config.eta) * value * value;
    }
    diagnostics.objective_residual =
        static_cast<double>(residual_objective);
    diagnostics.objective_smoothness =
        static_cast<double>(smoothness_objective);
    diagnostics.objective_amplitude =
        static_cast<double>(amplitude_objective);
    diagnostics.objective_total =
        diagnostics.objective_residual +
        diagnostics.objective_smoothness +
        diagnostics.objective_amplitude;
    return diagnostics.valid != 0 && diagnostics.converged != 0;
}

} // namespace RegularizedFacePairing

#endif
