#ifndef DUAL_U_COUPLING_H
#define DUAL_U_COUPLING_H

#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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

struct FinalLimitedDiagnostics {
    int valid;
    double target_linf;
    double residual_before_linf;
    double residual_after_linf;
    double minimum_scale;
    double correction_l2;
    double correction_linf;
    double candidate_min;
    long long corrected_cell_count;
    long long limited_cell_count;
    long long unresolved_cell_count;
    long long roundoff_zeroed_count;
    double roundoff_zeroed_mass;
    int failure_local_ix;

    FinalLimitedDiagnostics()
        : valid(1), target_linf(0.0), residual_before_linf(0.0),
          residual_after_linf(0.0), minimum_scale(1.0),
          correction_l2(0.0), correction_linf(0.0),
          candidate_min(std::numeric_limits<double>::infinity()),
          corrected_cell_count(0), limited_cell_count(0),
          unresolved_cell_count(0), roundoff_zeroed_count(0),
          roundoff_zeroed_mass(0.0), failure_local_ix(-1) {}
};

inline void merge_final_limited_diagnostics(
    FinalLimitedDiagnostics& target,
    const FinalLimitedDiagnostics& source)
{
    target.valid = std::min(target.valid, source.valid);
    target.target_linf = std::max(target.target_linf, source.target_linf);
    target.residual_before_linf = std::max(
        target.residual_before_linf, source.residual_before_linf);
    target.residual_after_linf = std::max(
        target.residual_after_linf, source.residual_after_linf);
    target.minimum_scale = std::min(target.minimum_scale,
                                    source.minimum_scale);
    target.correction_l2 += source.correction_l2;
    target.correction_linf = std::max(target.correction_linf,
                                      source.correction_linf);
    target.candidate_min = std::min(target.candidate_min,
                                    source.candidate_min);
    target.corrected_cell_count += source.corrected_cell_count;
    target.limited_cell_count += source.limited_cell_count;
    target.unresolved_cell_count += source.unresolved_cell_count;
    target.roundoff_zeroed_count += source.roundoff_zeroed_count;
    target.roundoff_zeroed_mass += source.roundoff_zeroed_mass;
    if (target.failure_local_ix < 0 && source.failure_local_ix >= 0)
        target.failure_local_ix = source.failure_local_ix;
}

inline void merge_diagnostics(Diagnostics& target,
                              const Diagnostics& source)
{
    target.target_replay_linf = std::max(target.target_replay_linf,
                                          source.target_replay_linf);
    target.target_replay_scale = std::max(target.target_replay_scale,
                                           source.target_replay_scale);
    target.legacy_current_linf = std::max(target.legacy_current_linf,
                                           source.legacy_current_linf);
    target.dual_current_linf = std::max(target.dual_current_linf,
                                         source.dual_current_linf);
    target.correction_l2 += source.correction_l2;
    target.correction_linf = std::max(target.correction_linf,
                                       source.correction_linf);
    target.corrected_cell_count += source.corrected_cell_count;
    if (!source.valid &&
        (target.valid || source.failure_local_ix < target.failure_local_ix)) {
        target.valid = 0;
        target.failure_subtype = source.failure_subtype;
        target.failure_local_ix = source.failure_local_ix;
        target.failure_target = source.failure_target;
        target.failure_replay = source.failure_replay;
        target.failure_legacy = source.failure_legacy;
        target.failure_residual = source.failure_residual;
        target.failure_denominator = source.failure_denominator;
        target.failure_maximum_coefficient =
            source.failure_maximum_coefficient;
        target.failure_support_floor = source.failure_support_floor;
        target.failure_scale = source.failure_scale;
    }
}

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
    const std::vector<double>& energy_current_weight,
    const std::vector<double>& legacy_coefficient,
    std::vector<double>& dual_coefficient,
    std::vector<double>& legacy_current,
    std::vector<double>& dual_current)
{
    Diagnostics diagnostics;
    const size_t faces_per_x = static_cast<size_t>(Param::Nv + 1) * Param::Nmu;
    const size_t expected = static_cast<size_t>(nx_local) * faces_per_x;
    if (legacy_coefficient.size() != expected ||
        target_current.size() < static_cast<size_t>(nx_local) ||
        target_current_replay.size() < static_cast<size_t>(nx_local) ||
        energy_current_weight.size() < Param::Nvmu ||
        !(dx > 0.0) || charge == 0.0 || !(mass > 0.0)) {
        diagnostics.valid = 0;
        diagnostics.failure_subtype = FAILURE_INPUT_CONTRACT;
        return diagnostics;
    }

    // Caller-owned work arrays are reused across midpoint iterations.  Keep
    // the copy explicit, rather than assigning temporary vectors, to avoid
    // repeated allocation and zero-fill on the production path.
    if (dual_coefficient.size() != expected) dual_coefficient.resize(expected);
    std::copy(legacy_coefficient.begin(), legacy_coefficient.end(),
              dual_coefficient.begin());
    const size_t cells = static_cast<size_t>(nx_local);
    if (legacy_current.size() != cells) legacy_current.resize(cells);
    if (dual_current.size() != cells) dual_current.resize(cells);
    std::fill(legacy_current.begin(), legacy_current.end(), 0.0);
    std::fill(dual_current.begin(), dual_current.end(), 0.0);

    const double current_factor = charge / (mass * Const::c * dx);
    int thread_count = 1;
#ifdef _OPENMP
    thread_count = omp_get_max_threads();
#endif
    std::vector<Diagnostics> thread_diagnostics(
        static_cast<size_t>(thread_count));

    #pragma omp parallel
    {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif
    Diagnostics& local = thread_diagnostics[static_cast<size_t>(thread_id)];
    #pragma omp for schedule(static)
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
        local.target_replay_linf = std::max(
            local.target_replay_linf, std::fabs(target - replay));
        local.target_replay_scale = std::max(
            local.target_replay_scale,
            std::max(std::fabs(target), std::fabs(replay)));
        local.legacy_current_linf = std::max(
            local.legacy_current_linf, std::fabs(target - legacy));

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
                const double support = std::max(
                    std::fabs(legacy_coefficient[id]), support_floor);
                const double functional = energy_current_weight[idx2(jf, k)];
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
                local.valid = 0;
                if (local.failure_subtype == FAILURE_NONE) {
                    local.failure_subtype = FAILURE_GRAM_DEGENERATE;
                    local.failure_local_ix = ix;
                    local.failure_target = target;
                    local.failure_replay = replay;
                    local.failure_legacy = legacy;
                    local.failure_residual = residual;
                    local.failure_denominator = denominator;
                    local.failure_maximum_coefficient =
                        maximum_coefficient;
                    local.failure_support_floor = support_floor;
                    local.failure_scale = scale;
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
                    const double support = std::max(
                        std::fabs(legacy_coefficient[id]), support_floor);
                    const double functional =
                        energy_current_weight[idx2(jf, k)];
                    const double correction = residual * functional * support /
                                              denominator;
                    const double corrected = dual_coefficient[id] + correction;
                    if (!std::isfinite(correction) ||
                        !std::isfinite(corrected)) {
                        local.valid = 0;
                        if (local.failure_subtype == FAILURE_NONE) {
                            local.failure_subtype =
                                FAILURE_CORRECTION_NONFINITE;
                            local.failure_local_ix = ix;
                            local.failure_target = target;
                            local.failure_replay = replay;
                            local.failure_legacy = legacy;
                            local.failure_residual = residual;
                            local.failure_denominator = denominator;
                            local.failure_maximum_coefficient =
                                maximum_coefficient;
                            local.failure_support_floor = support_floor;
                            local.failure_scale = scale;
                        }
                        continue;
                    }
                    dual_coefficient[id] = corrected;
                    correction_square += correction * correction;
                    correction_linf = std::max(correction_linf,
                                                std::fabs(correction));
                }
            }
            local.correction_l2 += correction_square;
            local.correction_linf = std::max(
                local.correction_linf, correction_linf);
            if (correction_linf > 0.0) ++local.corrected_cell_count;
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
        local.dual_current_linf = std::max(
            local.dual_current_linf, std::fabs(target - dual));
        if (!std::isfinite(dual)) {
            local.valid = 0;
            if (local.failure_subtype == FAILURE_NONE) {
                local.failure_subtype = FAILURE_DUAL_CURRENT_NONFINITE;
                local.failure_local_ix = ix;
                local.failure_target = target;
                local.failure_replay = replay;
                local.failure_legacy = legacy;
                local.failure_residual = residual;
                local.failure_denominator = denominator;
                local.failure_maximum_coefficient = maximum_coefficient;
                local.failure_support_floor = support_floor;
                local.failure_scale = scale;
            }
        }
    }
    } // OpenMP parallel region
    for (size_t thread = 0; thread < thread_diagnostics.size(); ++thread)
        merge_diagnostics(diagnostics, thread_diagnostics[thread]);
    diagnostics.correction_l2 = std::sqrt(diagnostics.correction_l2);
    return diagnostics;
}

// Match the final limited x-current with a conservative u-space transfer.
// Each populated donor cell selects at most one adjacent internal u face whose
// transfer changes the energy current in the requested direction.  The target
// residual is then distributed in one pass over the available donor capacity.
// This avoids a global positivity scale being pinned by empty velocity tails
// and avoids the repeated full-plane active-set scans used by the first 13.3
// implementation.  Velocity-boundary faces remain untouched.
const double kFinalPairingDonorRelativeFloor = 1.0e-12;

// Compute the current interval reachable by the same conservative donor
// transfers used by apply_final_limited_capacity_pairing().  The two bounds
// are independent one-sided capacities from the supplied candidate state;
// no mass, flux, or current is modified.
inline bool compute_final_pairing_current_bounds(
    int nx_local, double h,
    const std::vector<double>& acceleration,
    const std::vector<double>& energy_current_weight,
    const std::vector<double>& upar_widths,
    const std::vector<double>& candidate_mass,
    std::vector<double>& lower_delta_current,
    std::vector<double>& upper_delta_current)
{
    const size_t expected_cells =
        static_cast<size_t>(std::max(0, nx_local)) * Param::Nvmu;
    if (nx_local < 0 || !(h > 0.0) ||
        acceleration.size() < static_cast<size_t>(nx_local) ||
        energy_current_weight.size() < Param::Nvmu ||
        upar_widths.size() < static_cast<size_t>(Param::Nv) ||
        candidate_mass.size() < expected_cells)
        return false;

    lower_delta_current.assign(static_cast<size_t>(nx_local), 0.0);
    upper_delta_current.assign(static_cast<size_t>(nx_local), 0.0);
    int valid = 1;
    #pragma omp parallel for schedule(static) reduction(min:valid)
    for (int ix = 0; ix < nx_local; ++ix) {
        int cell_valid = 1;
        const size_t cell_base = static_cast<size_t>(ix) * Param::Nvmu;
        const double a = acceleration[static_cast<size_t>(ix)];
        if (!std::isfinite(a)) {
            valid = 0;
            continue;
        }
        if (a == 0.0) continue;

        double maximum_mass = 0.0;
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const double mass =
                    candidate_mass[cell_base + idx2(j, k)];
                if (!std::isfinite(mass)) {
                    cell_valid = 0;
                    continue;
                }
                maximum_mass = std::max(maximum_mass, std::max(0.0, mass));
            }
        if (!cell_valid) {
            valid = 0;
            continue;
        }
        const double donor_floor = std::max(
            64.0 * std::numeric_limits<double>::denorm_min(),
            kFinalPairingDonorRelativeFloor * maximum_mass);
        const double inverse_h_a = 1.0 / (h * a);
        long double positive_capacity = 0.0L;
        long double negative_capacity = 0.0L;
        for (int j = 0; j < Param::Nv; ++j) {
            if (!(upar_widths[static_cast<size_t>(j)] > 0.0) ||
                !std::isfinite(upar_widths[static_cast<size_t>(j)])) {
                cell_valid = 0;
                continue;
            }
            for (int k = 0; k < Param::Nmu; ++k) {
                const double available = std::max(
                    0.0, candidate_mass[cell_base + idx2(j, k)] -
                             donor_floor);
                if (!(available > 0.0)) continue;
                double best_positive = 0.0;
                double best_negative = 0.0;
                if (j > 0) {
                    const double gain =
                        -energy_current_weight[idx2(j, k)] * inverse_h_a;
                    if (gain > 0.0)
                        best_positive = std::max(best_positive, gain);
                    else
                        best_negative = std::min(best_negative, gain);
                }
                if (j + 1 < Param::Nv) {
                    const double gain =
                        energy_current_weight[idx2(j + 1, k)] * inverse_h_a;
                    if (gain > 0.0)
                        best_positive = std::max(best_positive, gain);
                    else
                        best_negative = std::min(best_negative, gain);
                }
                positive_capacity +=
                    static_cast<long double>(best_positive) * available;
                negative_capacity +=
                    static_cast<long double>(best_negative) * available;
            }
        }
        lower_delta_current[static_cast<size_t>(ix)] =
            static_cast<double>(negative_capacity);
        upper_delta_current[static_cast<size_t>(ix)] =
            static_cast<double>(positive_capacity);
        if (!std::isfinite(lower_delta_current[static_cast<size_t>(ix)]) ||
            !std::isfinite(upper_delta_current[static_cast<size_t>(ix)]))
            cell_valid = 0;
        valid = std::min(valid, cell_valid);
    }
    return valid != 0;
}

inline FinalLimitedDiagnostics apply_final_limited_capacity_pairing(
    int nx_local, double h,
    const std::vector<double>& target_current,
    const std::vector<double>& acceleration,
    const std::vector<double>& energy_current_weight,
    const std::vector<double>& upar_widths,
    std::vector<double>& candidate_mass,
    std::vector<double>& coefficient,
    std::vector<double>& flux)
{
    FinalLimitedDiagnostics diagnostics;
    const size_t faces_per_x =
        static_cast<size_t>(Param::Nv + 1) * Param::Nmu;
    const size_t expected_faces =
        static_cast<size_t>(nx_local) * faces_per_x;
    const size_t expected_cells =
        static_cast<size_t>(nx_local) * Param::Nvmu;
    if (nx_local < 0 || !(h > 0.0) ||
        target_current.size() < static_cast<size_t>(nx_local) ||
        acceleration.size() < static_cast<size_t>(nx_local) ||
        energy_current_weight.size() < Param::Nvmu ||
        upar_widths.size() < static_cast<size_t>(Param::Nv) ||
        candidate_mass.size() < expected_cells ||
        coefficient.size() != expected_faces ||
        flux.size() != expected_faces) {
        diagnostics.valid = 0;
        return diagnostics;
    }

    int thread_count = 1;
#ifdef _OPENMP
    thread_count = omp_get_max_threads();
#endif
    std::vector<FinalLimitedDiagnostics> thread_diagnostics(
        static_cast<size_t>(thread_count));

    #pragma omp parallel
    {
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif
    FinalLimitedDiagnostics& local =
        thread_diagnostics[static_cast<size_t>(thread_id)];
    // Reused by the OpenMP worker across solver calls.  Capturing the original
    // donor capacities prevents an incoming correction from being re-used as
    // outgoing mass later in the same simultaneous update.
    static thread_local std::vector<double> donor_capacity;
    if (donor_capacity.size() != Param::Nvmu)
        donor_capacity.assign(Param::Nvmu, 0.0);

    #pragma omp for schedule(static)
    for (int ix = 0; ix < nx_local; ++ix) {
        const size_t face_base = static_cast<size_t>(ix) * faces_per_x;
        const size_t cell_base = static_cast<size_t>(ix) * Param::Nvmu;
        double current = 0.0;
        for (int jf = 1; jf < Param::Nv; ++jf) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t local_face =
                    static_cast<size_t>(jf) * Param::Nmu + k;
                const double value = coefficient[face_base + local_face];
                current += energy_current_weight[idx2(jf, k)] * value;
            }
        }

        const double target = target_current[static_cast<size_t>(ix)];
        double residual = target - current;
        local.target_linf = std::max(local.target_linf, std::fabs(target));
        local.residual_before_linf = std::max(
            local.residual_before_linf, std::fabs(residual));
        const double residual_scale = std::max(
            1.0, std::max(std::fabs(target), std::fabs(current)));
        const double residual_floor =
            4096.0 * std::numeric_limits<double>::epsilon() * residual_scale;
        const double a = acceleration[static_cast<size_t>(ix)];
        if (!std::isfinite(a) || !std::isfinite(residual)) {
            local.valid = 0;
            if (local.failure_local_ix < 0) local.failure_local_ix = ix;
            continue;
        }

        double maximum_mass = 0.0;
        bool mass_finite = true;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double value =
                    candidate_mass[cell_base + idx2(j, k)];
                mass_finite = mass_finite && std::isfinite(value);
                maximum_mass = std::max(maximum_mass, std::max(0.0, value));
                local.candidate_min = std::min(local.candidate_min, value);
            }
        }
        if (!mass_finite) {
            local.valid = 0;
            if (local.failure_local_ix < 0) local.failure_local_ix = ix;
            continue;
        }
        if (std::fabs(residual) <= residual_floor || a == 0.0) {
            local.residual_after_linf = std::max(
                local.residual_after_linf, std::fabs(residual));
            if (a == 0.0 && std::fabs(residual) > residual_floor)
                ++local.unresolved_cell_count;
            continue;
        }

        // The correction support is a subset of the populated thermal body.
        // This threshold is only used by the pairing correction; it does not
        // clip or otherwise modify the evolved distribution.
        const double donor_floor = std::max(
            64.0 * std::numeric_limits<double>::denorm_min(),
            kFinalPairingDonorRelativeFloor * maximum_mass);
        const double requested = std::fabs(residual);
        const double requested_sign = residual > 0.0 ? 1.0 : -1.0;
        const double inverse_h_a = 1.0 / (h * a);

        // Return the adjacent face and signed current gain per unit donor mass.
        // A donor is used at most once, so all capacities are independent and
        // the simultaneous update is positivity preserving by construction.
        const auto donor_choice = [&](int j, int k, int& face,
                                      double& gain) {
            face = -1;
            gain = 0.0;
            if (j > 0) {
                const double lower_gain =
                    -energy_current_weight[idx2(j, k)] * inverse_h_a;
                if (requested_sign * lower_gain > 0.0) {
                    face = j;
                    gain = lower_gain;
                }
            }
            if (j + 1 < Param::Nv) {
                const double upper_gain =
                    energy_current_weight[idx2(j + 1, k)] * inverse_h_a;
                if (requested_sign * upper_gain > 0.0 &&
                    std::fabs(upper_gain) > std::fabs(gain)) {
                    face = j + 1;
                    gain = upper_gain;
                }
            }
        };

        long double capacity_current = 0.0L;
        for (int j = 0; j < Param::Nv; ++j) {
            const double donor_width =
                upar_widths[static_cast<size_t>(j)];
            if (!(donor_width > 0.0) || !std::isfinite(donor_width)) {
                local.valid = 0;
                if (local.failure_local_ix < 0) local.failure_local_ix = ix;
                continue;
            }
            for (int k = 0; k < Param::Nmu; ++k) {
                const double donor_mass =
                    candidate_mass[cell_base + idx2(j, k)];
                const double available = std::max(0.0, donor_mass - donor_floor);
                donor_capacity[idx2(j, k)] = available;
                if (!(available > 0.0)) continue;
                int face = -1;
                double gain = 0.0;
                donor_choice(j, k, face, gain);
                if (face >= 0 && std::isfinite(gain))
                    capacity_current +=
                        static_cast<long double>(std::fabs(gain)) *
                        static_cast<long double>(available);
            }
        }
        if (!local.valid) continue;

        const double capacity = static_cast<double>(capacity_current);
        if (!(capacity > 0.0) || !std::isfinite(capacity)) {
            local.residual_after_linf = std::max(
                local.residual_after_linf, requested);
            ++local.unresolved_cell_count;
            local.minimum_scale = 0.0;
            continue;
        }

        const double achieved_fraction = std::min(1.0, capacity / requested);
        // If all available mass is required, retain a roundoff margin at the
        // donor floor.  An unresolved residual is preferable to a negative
        // final candidate.
        double capacity_fraction = std::min(1.0, requested / capacity);
        if (capacity_fraction >= 1.0)
            capacity_fraction =
                1.0 - 128.0 * std::numeric_limits<double>::epsilon();
        local.minimum_scale = std::min(
            local.minimum_scale, achieved_fraction);

        long double applied_current_ld = 0.0L;
        double correction_square = 0.0;
        double correction_linf = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t donor_id = cell_base + idx2(j, k);
                const double available = donor_capacity[idx2(j, k)];
                if (!(available > 0.0)) continue;
                int face = -1;
                double gain = 0.0;
                donor_choice(j, k, face, gain);
                if (face < 0 || !std::isfinite(gain)) continue;

                const double moved_mass = capacity_fraction * available;
                if (!(moved_mass > 0.0)) continue;
                const int receiver_j = face == j ? j - 1 : j + 1;
                const double transfer =
                    receiver_j > j ? moved_mass : -moved_mass;
                const double delta_coefficient = transfer * inverse_h_a;
                const size_t face_id = face_base +
                    static_cast<size_t>(face) * Param::Nmu + k;
                coefficient[face_id] += delta_coefficient;
                flux[face_id] += a * delta_coefficient;
                candidate_mass[donor_id] -= moved_mass;
                candidate_mass[cell_base + idx2(receiver_j, k)] += moved_mass;
                const double current_change =
                    energy_current_weight[idx2(face, k)] * delta_coefficient;
                applied_current_ld += static_cast<long double>(current_change);
                correction_square += delta_coefficient * delta_coefficient;
                correction_linf = std::max(
                    correction_linf, std::fabs(delta_coefficient));
            }
        }

        residual -= static_cast<double>(applied_current_ld);

        local.residual_after_linf = std::max(
            local.residual_after_linf, std::fabs(residual));
        local.correction_l2 += correction_square;
        local.correction_linf = std::max(
            local.correction_linf, correction_linf);
        if (correction_linf > 0.0) ++local.corrected_cell_count;
        if (achieved_fraction < 1.0) ++local.limited_cell_count;
        if (std::fabs(residual) > residual_floor)
            ++local.unresolved_cell_count;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double value =
                    candidate_mass[cell_base + idx2(j, k)];
                if (!std::isfinite(value)) {
                    local.valid = 0;
                    if (local.failure_local_ix < 0)
                        local.failure_local_ix = ix;
                    continue;
                }
                if (value < 0.0) {
                    const double tolerance =
                        1024.0 * std::numeric_limits<double>::epsilon() *
                            std::max(maximum_mass, std::fabs(value)) +
                        64.0 * std::numeric_limits<double>::denorm_min();
                    if (std::fpclassify(value) == FP_SUBNORMAL ||
                        -value <= tolerance) {
                        candidate_mass[cell_base + idx2(j, k)] = 0.0;
                        ++local.roundoff_zeroed_count;
                        local.roundoff_zeroed_mass += -value;
                    } else {
                        local.valid = 0;
                        if (local.failure_local_ix < 0)
                            local.failure_local_ix = ix;
                    }
                }
                local.candidate_min = std::min(
                    local.candidate_min,
                    candidate_mass[cell_base + idx2(j, k)]);
            }
        }
    }
    } // OpenMP parallel region

    for (size_t thread = 0; thread < thread_diagnostics.size(); ++thread)
        merge_final_limited_diagnostics(
            diagnostics, thread_diagnostics[thread]);
    diagnostics.correction_l2 = std::sqrt(diagnostics.correction_l2);
    return diagnostics;
}

} // namespace DualUCoupling

#endif
