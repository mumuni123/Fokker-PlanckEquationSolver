#ifndef TAIL_MOMENT_CONSTRAINT_H
#define TAIL_MOMENT_CONSTRAINT_H

#include "parameters.h"

#include <cstddef>
#include <vector>

// Shared tail representation constraint module (sections 7.5, 7.10 and
// 14.5.1): the seven constrained moments N, Px, Jx, K, Pixx, Piperp and Xw
// are built with exactly one production formula, and the deterministic
// conservative support compression is implemented exactly once here.
// BulkTailConverter (six moments) and TailPopulationController (seven
// moments) both call into this module; no second copy of the formulas is
// allowed.
//
// Units follow section 4.1 (per unit transverse area):
//   n      : m^-2
//   px     : kg m s^-1 m^-2
//   jx     : A m^-1 (cell-integrated x current contribution)
//   ke     : J m^-2
//   pixx   : J m^-2 (parallel pressure contribution)
//   piperp : J m^-2 (perpendicular pressure contribution)
//   xw     : m^-1 (first spatial moment, section 7.10)
struct TailMoment7 {
    double n;
    double px;
    double jx;
    double ke;
    double pixx;
    double piperp;
    double xw;
    TailMoment7()
        : n(0.0), px(0.0), jx(0.0), ke(0.0), pixx(0.0), piperp(0.0),
          xw(0.0)
    {}
};

// Single production per-macro-particle moment formula (sections 4.1/7.6):
// weight w at position x with p/(m_e c) = (ux, uy, uz).
void tail_particle_moments(double weight, double x, double ux, double uy,
                           double uz, TailMoment7& m);

// Deterministic conservative support compression (section 7.5): given the
// per-unit-weight moment columns `cols` (cols[q] has `rows` entries, one per
// constrained moment), the current nonnegative `weights` and the reference
// moment vector `ref`, reduce the active supports to at most max_support
// while preserving ref to `tolerance` relative to |ref| per row.  On success
// the reduced weights are written into `weights` (non-active entries become
// zero) and true is returned.  On failure false is returned and `weights`
// is left unchanged (the caller keeps the original particles, section
// 7.10).
bool tail_compress_moment_supports(
    const std::vector<std::vector<double> >& cols,
    std::vector<double>& weights, const std::vector<double>& ref,
    size_t max_support, double tolerance);

// Deterministic non-negative moment fit used by the optional subcell tail
// loader.  `prior` supplies a geometric quadrature preference; feasibility
// is accepted only after the physical (unscaled) moment residual satisfies
// `tolerance`.  It never creates negative weights.
bool tail_solve_nonnegative_moment_weights(
    const std::vector<std::vector<double> >& cols,
    const std::vector<double>& ref, const std::vector<double>& prior,
    std::vector<double>& weights, double tolerance);

#endif
