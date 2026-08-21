#include "collision_coefficients.h"
#include "parameters.h"

#include <cmath>

CylindricalCollisionCoefficients ZeroCollisionCoefficients::evaluate(
    double, double, double, double, const LocalCollisionMoments&) const
{
    return { 0.0, 0.0, 0.0, 0.0, 0.0 };
}

std::string ZeroCollisionCoefficients::name() const { return "none"; }
CollisionCoefficientMode ZeroCollisionCoefficients::mode() const
{
    return CollisionCoefficientMode::NONE;
}

PrescribedCollisionCoefficients::PrescribedCollisionCoefficients(
    const CylindricalCollisionCoefficients& coefficients)
    : coefficients_(coefficients)
{}

CylindricalCollisionCoefficients PrescribedCollisionCoefficients::evaluate(
    double, double, double, double, const LocalCollisionMoments&) const
{
    return coefficients_;
}

std::string PrescribedCollisionCoefficients::name() const { return "prescribed"; }
CollisionCoefficientMode PrescribedCollisionCoefficients::mode() const
{
    return CollisionCoefficientMode::PRESCRIBED;
}

MomentClosureCollisionCoefficients::MomentClosureCollisionCoefficients(
    double coulomb_log)
    : coulomb_log_(coulomb_log)
{}

CylindricalCollisionCoefficients
MomentClosureCollisionCoefficients::evaluate(
    double, double u_parallel, double u_perp, double,
    const LocalCollisionMoments& moments) const
{
    CylindricalCollisionCoefficients c = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    const double n = moments.density;
    const double u_mean = moments.u_parallel_mean;
    const double ke_density = moments.kinetic_energy_density;
    if (!(n > 0.0) || !(ke_density > 0.0)) return c;
    // Non-relativistic effective temperature from the local moments:
    //   T_eff = (2/3)(KE_density/n - 0.5 m c^2 u_mean^2),
    //   u_th^2 = T_eff/(m c^2).
    const double random_ke = ke_density / n -
                             0.5 * Const::me * Const::c * Const::c *
                                 u_mean * u_mean;
    if (!(random_ke > 0.0)) return c;
    // Section 19.3: the operator equilibrium must be the discrete
    // Maxwellian whose moments reproduce the current distribution.  The
    // collision scan fills u_th2_closure from the grid self-consistency
    // table; without it the closure temperature chases the undercounted
    // grid moments and the bulk energy drains monotonically (no-beam 40fs
    // cooled 39%).
    const double u_th2 =
        (moments.u_th2_closure > 0.0)
            ? moments.u_th2_closure
            : (2.0 / 3.0) * random_ke / (Const::me * Const::c * Const::c);
    const double u_th = std::sqrt(u_th2);
    if (!(u_th > 0.0)) return c;
    // nu0 = n e^4 lnL / (4 pi eps0^2 m^2 c^3 u_th^3)  [1/s].
    const double nu0 =
        n * Const::qe * Const::qe * Const::qe * Const::qe *
        coulomb_log_ /
        (4.0 * Const::pi * Const::eps0 * Const::eps0 *
         Const::me * Const::me * Const::c * Const::c * Const::c *
         u_th2 * u_th);
    const double u_abs =
        std::sqrt(u_parallel * u_parallel + u_perp * u_perp);
    const double rate =
        nu0 / std::pow(1.0 + u_abs * u_abs / u_th2, 1.5);
    c.a_parallel = -rate * (u_parallel - u_mean);
    c.a_perp = -rate * u_perp;
    c.d_parallel_parallel = rate * u_th2;
    c.d_parallel_perp = 0.0;
    c.d_perp_perp = rate * u_th2;
    return c;
}

std::string MomentClosureCollisionCoefficients::name() const
{
    return "moment-closure";
}

CollisionCoefficientMode MomentClosureCollisionCoefficients::mode() const
{
    return CollisionCoefficientMode::MOMENT_CLOSURE;
}
