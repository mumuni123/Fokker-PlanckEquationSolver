#ifndef COLLISION_COEFFICIENTS_H
#define COLLISION_COEFFICIENTS_H

#include <memory>
#include <string>

// Section 10.2.1 coefficient-mode contract: the manifest must select one of
// these modes.  ZeroCollisionCoefficients is the collisionless path; the
// production bulk collision starts with moment_closure (H7) and the full
// self_consistent_landau is H8+.
enum class CollisionCoefficientMode {
    NONE,
    PRESCRIBED,
    MOMENT_CLOSURE,
    SELF_CONSISTENT_LANDAU
};

inline const char* collision_mode_name(CollisionCoefficientMode mode)
{
    switch (mode) {
        case CollisionCoefficientMode::NONE: return "none";
        case CollisionCoefficientMode::PRESCRIBED: return "prescribed";
        case CollisionCoefficientMode::MOMENT_CLOSURE: return "moment-closure";
        case CollisionCoefficientMode::SELF_CONSISTENT_LANDAU:
            return "self-consistent-landau";
    }
    return "unknown";
}

struct LocalCollisionMoments {
    double density;
    double u_parallel_mean;
    double kinetic_energy_density;
    // Section 19.3: closure temperature made self-consistent with the
    // discrete grid moments (u_th^2 of the operator equilibrium).  Filled
    // by the collision half-step scan via
    // CylindricalVelocityGrid::moment_closure_uth2_self_consistent; 0 means
    // "not set" and the moment-closure provider falls back to the
    // moments-derived value (gridless callers).
    double u_th2_closure;
    LocalCollisionMoments()
        : density(0.0), u_parallel_mean(0.0), kinetic_energy_density(0.0),
          u_th2_closure(0.0)
    {}
};

// Drift (a, per unit time in u = p/(m c)) and diffusion (d, per unit time
// in u^2) coefficients of the cylindrical FP operator (section 10.2).  The
// operator evaluates the provider at every velocity point and supports the
// non-diagonal d_parallel_perp term.
struct CylindricalCollisionCoefficients {
    double a_parallel;
    double a_perp;
    double d_parallel_parallel;
    double d_parallel_perp;
    double d_perp_perp;
};

class CollisionCoefficientProvider {
public:
    virtual ~CollisionCoefficientProvider() {}
    // Velocity-point evaluation: u_parallel/u_perp are the p/(m c) point
    // where the drift/diffusion coefficients are needed.  moments carries
    // the local cell moments for moment-closure modes.
    virtual CylindricalCollisionCoefficients evaluate(
        double x, double u_parallel, double u_perp, double time,
        const LocalCollisionMoments& moments) const = 0;
    virtual std::string name() const = 0;
    virtual CollisionCoefficientMode mode() const = 0;
};

class ZeroCollisionCoefficients : public CollisionCoefficientProvider {
public:
    CylindricalCollisionCoefficients evaluate(double, double, double, double,
                                               const LocalCollisionMoments&) const;
    std::string name() const;
    CollisionCoefficientMode mode() const;
};

class PrescribedCollisionCoefficients : public CollisionCoefficientProvider {
public:
    explicit PrescribedCollisionCoefficients(const CylindricalCollisionCoefficients& coefficients);
    CylindricalCollisionCoefficients evaluate(double, double, double, double,
                                               const LocalCollisionMoments&) const;
    std::string name() const;
    CollisionCoefficientMode mode() const;

private:
    CylindricalCollisionCoefficients coefficients_;
};

// Stage-H7 moment-closure provider (section 10.2.1 mode 2): the field is
// the local drifting Maxwellian built from the cell density, mean
// u_parallel and effective temperature.  The velocity-dependent collision
// rate mimics the Coulomb 1/u^3 scaling at high speed,
//   nu(u) = nu0 / (1 + (|u|/u_th)^2)^{3/2},
// with nu0 = n e^4 lnL / (4 pi eps0^2 m^2 c^3 u_th^3).  The drift and
// diffusion satisfy the Einstein relation with the drifting Maxwellian, so
// the equilibrium is exactly that Maxwellian:
//   A = -nu(u)(u - u_mean e_par),   D_parpar = D_perpperp = nu(u) u_th^2.
// This is an explicit closure approximation, not the full Landau operator
// (that is self_consistent_landau, H8+).
class MomentClosureCollisionCoefficients : public CollisionCoefficientProvider {
public:
    explicit MomentClosureCollisionCoefficients(double coulomb_log = 20.0);
    CylindricalCollisionCoefficients evaluate(
        double x, double u_parallel, double u_perp, double time,
        const LocalCollisionMoments& moments) const;
    std::string name() const;
    CollisionCoefficientMode mode() const;

private:
    double coulomb_log_;
};

#endif
