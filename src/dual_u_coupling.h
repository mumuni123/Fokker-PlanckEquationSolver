#ifndef DUAL_U_COUPLING_H
#define DUAL_U_COUPLING_H

#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace DualUCoupling {

enum FailureSubtype {
    FAILURE_NONE = 0,
    FAILURE_INPUT_CONTRACT = 1,
    FAILURE_GRAM_DEGENERATE = 2,
    FAILURE_CORRECTION_NONFINITE = 3,
    FAILURE_DUAL_CURRENT_NONFINITE = 4
};

struct Diagnostics {
    int valid;
    double target_replay_linf;
    double target_replay_scale;
    double legacy_current_linf;
    double dual_current_linf;
    double correction_l2;
    double correction_linf;
    long long corrected_cell_count;
    int failure_subtype;
    int failure_local_ix;
    double failure_target;
    double failure_replay;
    double failure_legacy;
    double failure_residual;
    double failure_denominator;
    double failure_maximum_coefficient;
    double failure_support_floor;
    double failure_scale;

    Diagnostics()
        : valid(1), target_replay_linf(0.0), target_replay_scale(0.0),
          legacy_current_linf(0.0), dual_current_linf(0.0),
          correction_l2(0.0), correction_linf(0.0),
          corrected_cell_count(0), failure_subtype(FAILURE_NONE),
          failure_local_ix(-1), failure_target(0.0), failure_replay(0.0),
          failure_legacy(0.0), failure_residual(0.0),
          failure_denominator(0.0), failure_maximum_coefficient(0.0),
          failure_support_floor(0.0), failure_scale(0.0) {}
};

// Fixed-state local rank-one update
//
//   A_new = A_legacy + b (w_N - w_E)^T / <l,b>,
//
// where l is the discrete kinetic-energy jump moment.  The implementation
// applies this operator to the supplied fixed midpoint state: target_current
// is the actual x-flux functional evaluated on that state, while
// legacy_coefficient is A_legacy M.  The returned coefficient is therefore a
// real shared u-face coefficient; multiplying it by acceleration updates M,
// and taking its kinetic-energy moment forms J_E.  No current or field is
// modified after transport.
inline Diagnostics apply_local_rank_one(
    const CylindricalVelocityGrid& grid, double charge, double mass,
    double dx, int nx_local, const std::vector<double>& target_current,
    const std::vector<double>& target_current_replay,
    const std::vector<double>& legacy_coefficient,
    std::vector<double>& dual_coefficient,
    std::vector<double>& legacy_current,
    std::vector<double>& dual_current)
{
    Diagnostics diagnostics;
    const size_t faces_per_x = static_cast<size_t>(Param::Nv + 1) * Param::Nmu;
    const size_t expected = static_cast<size_t>(nx_local) * faces_per_x;
    dual_coefficient = legacy_coefficient;
    legacy_current.assign(static_cast<size_t>(nx_local), 0.0);
    dual_current.assign(static_cast<size_t>(nx_local), 0.0);
    if (legacy_coefficient.size() != expected ||
        target_current.size() < static_cast<size_t>(nx_local) ||
        target_current_replay.size() < static_cast<size_t>(nx_local) ||
        !(dx > 0.0) || charge == 0.0 || !(mass > 0.0)) {
        diagnostics.valid = 0;
        diagnostics.failure_subtype = FAILURE_INPUT_CONTRACT;
        return diagnostics;
    }

    const double current_factor = charge / (mass * Const::c * dx);
    for (int ix = 0; ix < nx_local; ++ix) {
        double legacy = 0.0;
        double maximum_coefficient = 0.0;
        for (int jf = 1; jf < Param::Nv; ++jf) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = (static_cast<size_t>(ix) * (Param::Nv + 1) +
                                   static_cast<size_t>(jf)) * Param::Nmu + k;
                legacy += Stage5::delta_energy(grid, jf, k) *
                          legacy_coefficient[id];
                maximum_coefficient = std::max(maximum_coefficient,
                    std::fabs(legacy_coefficient[id]));
            }
        }
        legacy *= current_factor;
        legacy_current[static_cast<size_t>(ix)] = legacy;

        const double target = target_current[static_cast<size_t>(ix)];
        const double replay = target_current_replay[static_cast<size_t>(ix)];
        diagnostics.target_replay_linf = std::max(
            diagnostics.target_replay_linf, std::fabs(target - replay));
        diagnostics.target_replay_scale = std::max(
            diagnostics.target_replay_scale,
            std::max(std::fabs(target), std::fabs(replay)));
        diagnostics.legacy_current_linf = std::max(
            diagnostics.legacy_current_linf, std::fabs(target - legacy));

        // Freeze a compact support metric from the legacy face state.  A tiny
        // relative floor lets an empty internal face participate without
        // spreading the correction beyond this x cell or altering u bounds.
        const double support_floor = std::max(
            64.0 * std::numeric_limits<double>::denorm_min(),
            1.0e-14 * maximum_coefficient);
        // Positive Gram denominator for the energy-current functional.
        // a=current_factor*dK includes the complete physical scaling.  The
        // previous signed form was algebraically equivalent away from zero,
        // but compared a differently scaled denominator with a current and
        // could therefore reject a well-conditioned local solve.
        double denominator = 0.0;
        for (int jf = 1; jf < Param::Nv; ++jf) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = (static_cast<size_t>(ix) * (Param::Nv + 1) +
                                   static_cast<size_t>(jf)) * Param::Nmu + k;
                const double dke = Stage5::delta_energy(grid, jf, k);
                const double support = std::max(
                    std::fabs(legacy_coefficient[id]), support_floor);
                const double functional = current_factor * dke;
                denominator += functional * functional * support;
            }
        }
        const double residual = target - legacy;
        const double scale = std::max(1.0, std::max(std::fabs(target),
                                                   std::fabs(legacy)));
        if (!std::isfinite(denominator) || !std::isfinite(residual) ||
            !(denominator > 64.0 *
                std::numeric_limits<double>::denorm_min())) {
            if (std::fabs(residual) >
                4096.0 * std::numeric_limits<double>::epsilon() * scale) {
                diagnostics.valid = 0;
                if (diagnostics.failure_subtype == FAILURE_NONE) {
                    diagnostics.failure_subtype = FAILURE_GRAM_DEGENERATE;
                    diagnostics.failure_local_ix = ix;
                    diagnostics.failure_target = target;
                    diagnostics.failure_replay = replay;
                    diagnostics.failure_legacy = legacy;
                    diagnostics.failure_residual = residual;
                    diagnostics.failure_denominator = denominator;
                    diagnostics.failure_maximum_coefficient =
                        maximum_coefficient;
                    diagnostics.failure_support_floor = support_floor;
                    diagnostics.failure_scale = scale;
                }
            }
        } else {
            double correction_square = 0.0;
            double correction_linf = 0.0;
            for (int jf = 1; jf < Param::Nv; ++jf) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = (static_cast<size_t>(ix) *
                        (Param::Nv + 1) + static_cast<size_t>(jf)) *
                        Param::Nmu + k;
                    const double dke = Stage5::delta_energy(grid, jf, k);
                    const double support = std::max(
                        std::fabs(legacy_coefficient[id]), support_floor);
                    const double functional = current_factor * dke;
                    const double correction = residual * functional * support /
                                              denominator;
                    const double corrected = dual_coefficient[id] + correction;
                    if (!std::isfinite(correction) ||
                        !std::isfinite(corrected)) {
                        diagnostics.valid = 0;
                        if (diagnostics.failure_subtype == FAILURE_NONE) {
                            diagnostics.failure_subtype =
                                FAILURE_CORRECTION_NONFINITE;
                            diagnostics.failure_local_ix = ix;
                            diagnostics.failure_target = target;
                            diagnostics.failure_replay = replay;
                            diagnostics.failure_legacy = legacy;
                            diagnostics.failure_residual = residual;
                            diagnostics.failure_denominator = denominator;
                            diagnostics.failure_maximum_coefficient =
                                maximum_coefficient;
                            diagnostics.failure_support_floor = support_floor;
                            diagnostics.failure_scale = scale;
                        }
                        continue;
                    }
                    dual_coefficient[id] = corrected;
                    correction_square += correction * correction;
                    correction_linf = std::max(correction_linf,
                                                std::fabs(correction));
                }
            }
            diagnostics.correction_l2 += correction_square;
            diagnostics.correction_linf = std::max(
                diagnostics.correction_linf, correction_linf);
            if (correction_linf > 0.0) ++diagnostics.corrected_cell_count;
        }

        double dual = 0.0;
        for (int jf = 1; jf < Param::Nv; ++jf)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = (static_cast<size_t>(ix) * (Param::Nv + 1) +
                                   static_cast<size_t>(jf)) * Param::Nmu + k;
                dual += Stage5::delta_energy(grid, jf, k) *
                        dual_coefficient[id];
            }
        dual *= current_factor;
        dual_current[static_cast<size_t>(ix)] = dual;
        diagnostics.dual_current_linf = std::max(
            diagnostics.dual_current_linf, std::fabs(target - dual));
        if (!std::isfinite(dual)) {
            diagnostics.valid = 0;
            if (diagnostics.failure_subtype == FAILURE_NONE) {
                diagnostics.failure_subtype = FAILURE_DUAL_CURRENT_NONFINITE;
                diagnostics.failure_local_ix = ix;
                diagnostics.failure_target = target;
                diagnostics.failure_replay = replay;
                diagnostics.failure_legacy = legacy;
                diagnostics.failure_residual = residual;
                diagnostics.failure_denominator = denominator;
                diagnostics.failure_maximum_coefficient = maximum_coefficient;
                diagnostics.failure_support_floor = support_floor;
                diagnostics.failure_scale = scale;
            }
        }
    }
    diagnostics.correction_l2 = std::sqrt(diagnostics.correction_l2);
    return diagnostics;
}

} // namespace DualUCoupling

#endif
