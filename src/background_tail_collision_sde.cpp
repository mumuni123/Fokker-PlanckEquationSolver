#include "background_tail_collision_sde.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace {

inline std::uint64_t mix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline double uniform01(std::uint64_t& state)
{
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return static_cast<double>(z >> 11) *
           (1.0 / 9007199254740992.0);
}

// Two independent standard normals via Box-Muller.
inline void gaussian_pair(std::uint64_t& state, double& g0, double& g1)
{
    const double u0 = std::max(1.0e-15, uniform01(state));
    const double u1 = uniform01(state);
    const double r = std::sqrt(-2.0 * std::log(u0));
    g0 = r * std::cos(2.0 * Const::pi * u1);
    g1 = r * std::sin(2.0 * Const::pi * u1);
}

// Build a square root in the natural (parallel, radial-perpendicular) basis
// and rotate it to Cartesian velocity coordinates.  The Cartesian tensor is
// rank deficient away from the axis because it has no azimuthal diffusion;
// applying an unpivoted 3x3 Cholesky directly to that tensor is numerically
// fragile when one Cartesian component of the radial unit vector is small.
bool cylindrical_diffusion_sqrt(
    const CylindricalCollisionCoefficients& c, const double u[3],
    double b[3][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) b[i][j] = 0.0;
    }

    const double dpp = c.d_parallel_parallel;
    const double dpr = c.d_parallel_perp;
    const double drr = c.d_perp_perp;
    if (!std::isfinite(dpp) || !std::isfinite(dpr) ||
        !std::isfinite(drr)) {
        return false;
    }
    const double scale = std::max(
        std::max(std::fabs(dpp), std::fabs(drr)), std::fabs(dpr));
    if (!(scale > 0.0)) return dpp == 0.0 && dpr == 0.0 && drr == 0.0;
    const double tol = 256.0 * std::numeric_limits<double>::epsilon() * scale;
    if (dpp < -tol || drr < -tol) return false;

    const double a = std::max(0.0, dpp);
    const double cperp = std::max(0.0, drr);
    const double determinant = a * cperp - dpr * dpr;
    if (determinant < -tol * scale) return false;

    const double u_perp = std::hypot(u[1], u[2]);
    if (u_perp <= 1.0e-12) {
        // At the axis the radial direction is undefined.  Symmetry requires
        // zero parallel-radial cross diffusion and an isotropic split of the
        // perpendicular trace between y and z.
        if (std::fabs(dpr) > tol) return false;
        b[0][0] = std::sqrt(a);
        b[1][1] = std::sqrt(0.5 * cperp);
        b[2][2] = std::sqrt(0.5 * cperp);
        return true;
    }

    // Symmetric 2x2 eigendecomposition of [[a,dpr],[dpr,cperp]].
    // Clamp only roundoff-scale negative eigenvalues; a genuinely indefinite
    // provider tensor remains a hard failure.
    const double trace = a + cperp;
    const double gap = std::hypot(a - cperp, 2.0 * dpr);
    double lambda_hi = 0.5 * (trace + gap);
    double lambda_lo = 0.5 * (trace - gap);
    if (lambda_lo < -tol) return false;
    lambda_hi = std::max(0.0, lambda_hi);
    lambda_lo = std::max(0.0, lambda_lo);
    const double theta = 0.5 * std::atan2(2.0 * dpr, a - cperp);
    const double ct = std::cos(theta);
    const double st = std::sin(theta);
    const double sh = std::sqrt(lambda_hi);
    const double sl = std::sqrt(lambda_lo);
    const double radial_y = u[1] / u_perp;
    const double radial_z = u[2] / u_perp;

    // Columns are independent Wiener increments.  B B^T reproduces the
    // cylindrical 2x2 tensor in the parallel/radial plane exactly.
    b[0][0] = ct * sh;
    b[0][1] = -st * sl;
    const double radial_col0 = st * sh;
    const double radial_col1 = ct * sl;
    b[1][0] = radial_y * radial_col0;
    b[1][1] = radial_y * radial_col1;
    b[2][0] = radial_z * radial_col0;
    b[2][1] = radial_z * radial_col1;
    return true;
}

// Cartesian drift/diffusion from the cylindrical coefficients at the
// particle velocity (u_par along x, perpendicular plane y-z).
void cylindrical_to_cartesian(const CylindricalCollisionCoefficients& c,
                              const double u[3], double a_out[3],
                              double d_out[3][3])
{
    for (int i = 0; i < 3; ++i) a_out[i] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) d_out[i][j] = 0.0;
    }
    const double u_perp = std::sqrt(u[1] * u[1] + u[2] * u[2]);
    a_out[0] = c.a_parallel;
    d_out[0][0] = c.d_parallel_parallel;
    if (u_perp > 1.0e-12) {
        const double e1 = u[1] / u_perp;
        const double e2 = u[2] / u_perp;
        a_out[1] = c.a_perp * e1;
        a_out[2] = c.a_perp * e2;
        d_out[0][1] = c.d_parallel_perp * e1;
        d_out[1][0] = d_out[0][1];
        d_out[0][2] = c.d_parallel_perp * e2;
        d_out[2][0] = d_out[0][2];
        d_out[1][1] = c.d_perp_perp * e1 * e1;
        d_out[1][2] = c.d_perp_perp * e1 * e2;
        d_out[2][1] = d_out[1][2];
        d_out[2][2] = c.d_perp_perp * e2 * e2;
    } else {
        // Axis (u_perp -> 0): the perpendicular diffusion is isotropic in
        // the (y,z) plane with the same trace as the u_perp > 0 case,
        // i.e. D_yy + D_zz = D_perp_perp.
        a_out[1] = 0.0;
        a_out[2] = 0.0;
        d_out[1][1] = 0.5 * c.d_perp_perp;
        d_out[2][2] = 0.5 * c.d_perp_perp;
    }
}

} // namespace

bool sde_collide(BackgroundTailPIC& tail, const SpatialGrid& grid,
                 const TailCollisionRequest& request,
                 const CollisionCoefficientProvider& provider,
                 const std::vector<LocalCollisionMoments>& cell_moments,
                 TailCollisionDiagnostics& diagnostics)
{
    diagnostics = TailCollisionDiagnostics();
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    if (!(request.dt > 0.0) || tail.particles.empty()) {
        diagnostics.success = true;
        return true;
    }
    diagnostics.particle_count_before = tail.particles.size();
    diagnostics.particle_count_attempted = tail.particles.size();
    diagnostics.particle_count_limit = tail.particles.size();
    if (cell_moments.size() <
        static_cast<size_t>(std::max(0, grid.nx_local))) {
        diagnostics.success = false;
        diagnostics.failure_reason =
            TailCollisionFailureReason::InvalidRequest;
        return false;
    }
    const int nxl = grid.nx_local;
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        const BackgroundTailParticle& pp = tail.particles[p];
        const double gamma =
            std::sqrt(1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
        diagnostics.number_before += pp.weight;
        diagnostics.px_before += Const::me * Const::c * pp.weight * pp.ux;
        diagnostics.py_before += Const::me * Const::c * pp.weight * pp.uy;
        diagnostics.pz_before += Const::me * Const::c * pp.weight * pp.uz;
        diagnostics.ke_before += Const::me * Const::c * Const::c *
                                 pp.weight * (gamma - 1.0);
    }

    bool ok = true;
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        BackgroundTailParticle& pp = tail.particles[p];
        const int cell = static_cast<int>(
            std::floor(pp.x / grid.dx)) - grid.ix_start;
        const LocalCollisionMoments& moments =
            (cell >= 0 && cell < nxl)
                ? cell_moments[static_cast<size_t>(cell)]
                : cell_moments.front();
        const double u[3] = { pp.ux, pp.uy, pp.uz };
        const CylindricalCollisionCoefficients c = provider.evaluate(
            grid.x(grid.nghost + std::max(0, std::min(cell, nxl - 1))),
            pp.ux, std::sqrt(pp.uy * pp.uy + pp.uz * pp.uz), 0.0,
            moments);
        double a[3];
        double d[3][3];
        cylindrical_to_cartesian(c, u, a, d);

        // Relativistic Ito/Jacobian correction (1/2) dD/d u: approximate the
        // derivative of the Cartesian D tensor with central differences of
        // the provider along u_parallel and u_perp.
        const double h = 1.0e-4 * std::max(1.0, std::sqrt(
            u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
        const double up = std::sqrt(pp.uy * pp.uy + pp.uz * pp.uz);
        const CylindricalCollisionCoefficients cp = provider.evaluate(
            grid.x(grid.nghost + std::max(0, std::min(cell, nxl - 1))),
            pp.ux + h, up, 0.0, moments);
        const CylindricalCollisionCoefficients cm = provider.evaluate(
            grid.x(grid.nghost + std::max(0, std::min(cell, nxl - 1))),
            pp.ux - h, up, 0.0, moments);
        const CylindricalCollisionCoefficients cpp = provider.evaluate(
            grid.x(grid.nghost + std::max(0, std::min(cell, nxl - 1))),
            pp.ux, up + h, 0.0, moments);
        const CylindricalCollisionCoefficients cmm = provider.evaluate(
            grid.x(grid.nghost + std::max(0, std::min(cell, nxl - 1))),
            pp.ux, up - h, 0.0, moments);
        const double dd_dd_par = (cp.d_parallel_parallel -
                                  cm.d_parallel_parallel) / (2.0 * h);
        const double dd_perp_perp_perp =
            (cpp.d_perp_perp - cmm.d_perp_perp) / (2.0 * h);
        const double dd_cross_perp =
            (cpp.d_parallel_perp - cmm.d_parallel_perp) / (2.0 * h);
        // Leading-order relativistic Ito/Jacobian correction from the
        // cylindrical divergence of D:
        //   A_ito_par  = (1/2) d D_parpar / d u_par
        //   A_ito_perp = (1/2)(d D_perpperp / d u_perp + D_perpperp/u_perp)
        // plus the cross-divergence projection; the D_perpperp/u_perp term
        // is the standard cylindrical spurious drift.
        a[0] += 0.5 * dd_dd_par;
        if (up > 1.0e-12) {
            const double e1 = pp.uy / up;
            const double e2 = pp.uz / up;
            const double a_perp_ito =
                0.5 * (dd_perp_perp_perp + (d[1][1] + d[2][2]) / up);
            a[1] += a_perp_ito * e1;
            a[2] += a_perp_ito * e2;
            a[0] += 0.5 * dd_cross_perp;
        }

        double b[3][3];
        if (!cylindrical_diffusion_sqrt(c, u, b)) {
            ok = false;
            break;
        }
        std::uint64_t state = mix64(
            request.rng_seed_base ^
            (static_cast<std::uint64_t>(request.accepted_step) << 32) ^
            (static_cast<std::uint64_t>(request.collision_half) << 40) ^
            (pp.id * 0x9e3779b97f4a7c15ULL));
        const double dt = request.dt;
        double du[3] = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < 3; ++i) du[i] = a[i] * dt;
        double noise[3];
        gaussian_pair(state, noise[0], noise[1]);
        {
            double extra0 = 0.0;
            double extra1 = 0.0;
            gaussian_pair(state, extra0, extra1);
            noise[2] = extra0;
        }
        const double sqrt_dt = std::sqrt(dt);
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int i = 0; i < 3; ++i) s += b[j][i] * noise[i];
            du[j] += s * sqrt_dt;
        }
        const double du_mag = std::sqrt(
            du[0] * du[0] + du[1] * du[1] + du[2] * du[2]);
        diagnostics.max_du = std::max(diagnostics.max_du, du_mag);
        pp.ux += du[0];
        pp.uy += du[1];
        pp.uz += du[2];
    }
    if (!ok) {
        diagnostics.success = false;
        diagnostics.failure_reason =
            TailCollisionFailureReason::DiffusionFactorization;
        return false;
    }
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        const BackgroundTailParticle& pp = tail.particles[p];
        const double gamma =
            std::sqrt(1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
        diagnostics.number_after += pp.weight;
        diagnostics.px_after += Const::me * Const::c * pp.weight * pp.ux;
        diagnostics.py_after += Const::me * Const::c * pp.weight * pp.uy;
        diagnostics.pz_after += Const::me * Const::c * pp.weight * pp.uz;
        diagnostics.ke_after += Const::me * Const::c * Const::c *
                                pp.weight * (gamma - 1.0);
    }
    diagnostics.success = true;
    diagnostics.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    return true;
}
