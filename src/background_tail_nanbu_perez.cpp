#include "background_tail_nanbu_perez.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace {

// Counter-based deterministic RNG (SplitMix64), keyed by a 64-bit counter
// that includes the particle id, accepted step, collision half, substep and
// pairing pass (section 10.3.1 item 2: never the vector order or MPI rank).
inline std::uint64_t mix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline double uniform01(std::uint64_t& state)
{
    // Advance the counter state (SplitMix64): uniform01 must consume the
    // state, otherwise every call returns the same value.
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return static_cast<double>(z >> 11) *
           (1.0 / 9007199254740992.0);   // (0,1)
}

// Relativistic boost of a 4-velocity u (p/(m c)) by beta (v/c).
inline void lorentz_boost(const double u[3], const double beta[3],
                          double u_out[3])
{
    const double b2 = beta[0] * beta[0] + beta[1] * beta[1] +
                      beta[2] * beta[2];
    if (b2 <= 0.0) {
        u_out[0] = u[0];
        u_out[1] = u[1];
        u_out[2] = u[2];
        return;
    }
    const double gamma = 1.0 / std::sqrt(1.0 - b2);
    const double gamma_u = std::sqrt(
        1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    const double bdotu = beta[0] * u[0] + beta[1] * u[1] + beta[2] * u[2];
    const double coeff = (gamma - 1.0) / b2 * bdotu - gamma * gamma_u;
    for (int i = 0; i < 3; ++i) {
        u_out[i] = u[i] + coeff * beta[i];
    }
}

// Same-species Coulomb scattering strength (section 10.3.1 item 3).  The
// Rutherford rate scales like n q^4 lnL / (g^3); written in the COM-frame
// relative momentum u_rel = 2 u* with the relativistic gamma this is
//   s12 = n e^4 lnL / (4 pi eps0^2 m^2 c^3) * (gamma_rel^3 / u_rel^3) dt,
// regularised at small relative momentum by the cold-plasma cutoff u_cold
// (a separate, config-controlled parameter, section 10.3.1 item 6).
inline double s12_value(double n_t, double gamma_rel, double u_star,
                        double dt, double coulomb_log,
                        double cold_velocity_cutoff)
{
    // gamma_rel is clamped to >= 1 by the caller; the gamma_rel == 1 case
    // (identical velocities) is exactly where the cold-plasma cutoff must
    // bound the rate from below (1/u_cold^3), not suppress it to zero.
    if (!(n_t > 0.0) || !(gamma_rel >= 1.0) || !(dt > 0.0)) return 0.0;
    const double u_rel = std::max(2.0 * u_star, cold_velocity_cutoff);
    const double gamma_over_u = gamma_rel / u_rel;
    return n_t * Const::qe * Const::qe * Const::qe * Const::qe *
           coulomb_log / (4.0 * Const::pi * Const::eps0 * Const::eps0 *
                          Const::me * Const::me * Const::c * Const::c *
                          Const::c) *
           gamma_over_u * gamma_over_u * gamma_over_u * dt;
}

// Rotate the unit direction n by cos_theta/phi and return the rotated
// direction.  The reference axis avoids the near-parallel case.
inline void rotate_direction(const double n[3], double cos_theta,
                             double sin_theta, double cos_phi,
                             double sin_phi, double n_out[3])
{
    double ref[3] = { 0.0, 0.0, 1.0 };
    if (std::fabs(n[2]) > 0.9) ref[2] = 0.0, ref[0] = 1.0;
    double a[3];
    a[0] = ref[1] * n[2] - ref[2] * n[1];
    a[1] = ref[2] * n[0] - ref[0] * n[2];
    a[2] = ref[0] * n[1] - ref[1] * n[0];
    const double an = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    if (an <= 0.0) {
        n_out[0] = n[0];
        n_out[1] = n[1];
        n_out[2] = n[2];
        return;
    }
    for (int i = 0; i < 3; ++i) a[i] /= an;
    double b[3];
    b[0] = n[1] * a[2] - n[2] * a[1];
    b[1] = n[2] * a[0] - n[0] * a[2];
    b[2] = n[0] * a[1] - n[1] * a[0];
    for (int i = 0; i < 3; ++i) {
        n_out[i] = cos_theta * n[i] +
                   sin_theta * (cos_phi * a[i] + sin_phi * b[i]);
    }
}

// Sentoku--Kemp unequal-weight correction (EPOCH collisions.F90,
// weighted_particles_correction).  The lighter macro particle represents the
// fully scattered share.  The heavier particle keeps a single support whose
// target energy is the weighted mixture of its unscattered and scattered
// states.  A random transverse component restores that energy without
// materialising the residual share as a permanent macro particle.
void weighted_particle_correction(double scattered_fraction,
                                  const double u_old[3],
                                  const double u_scattered[3],
                                  std::uint64_t& rng_state,
                                  double u_out[3])
{
    const double r = std::max(0.0, std::min(1.0, scattered_fraction));
    const double gamma_old = std::sqrt(
        1.0 + u_old[0] * u_old[0] + u_old[1] * u_old[1] +
        u_old[2] * u_old[2]);
    const double gamma_scattered = std::sqrt(
        1.0 + u_scattered[0] * u_scattered[0] +
        u_scattered[1] * u_scattered[1] +
        u_scattered[2] * u_scattered[2]);
    const double gamma_target =
        (1.0 - r) * gamma_old + r * gamma_scattered;
    for (int i = 0; i < 3; ++i) {
        u_out[i] = (1.0 - r) * u_old[i] + r * u_scattered[i];
    }
    const double u2 = u_out[0] * u_out[0] + u_out[1] * u_out[1] +
                      u_out[2] * u_out[2];
    const double gamma_mixed = std::sqrt(1.0 + u2);
    const double delta2 = std::max(
        0.0, gamma_target * gamma_target - gamma_mixed * gamma_mixed);
    if (!(delta2 > 0.0)) return;

    double axis[3] = { 1.0, 0.0, 0.0 };
    if (u2 > 1.0e-30) {
        const double inv_u = 1.0 / std::sqrt(u2);
        axis[0] = u_out[0] * inv_u;
        axis[1] = u_out[1] * inv_u;
        axis[2] = u_out[2] * inv_u;
    }
    // Pick the Cartesian direction least aligned with axis, then construct
    // an orthonormal transverse pair.
    int ref_index = 0;
    if (std::fabs(axis[1]) < std::fabs(axis[ref_index])) ref_index = 1;
    if (std::fabs(axis[2]) < std::fabs(axis[ref_index])) ref_index = 2;
    double ref[3] = { 0.0, 0.0, 0.0 };
    ref[ref_index] = 1.0;
    double e2[3] = {
        axis[1] * ref[2] - axis[2] * ref[1],
        axis[2] * ref[0] - axis[0] * ref[2],
        axis[0] * ref[1] - axis[1] * ref[0]
    };
    const double e2_norm = std::sqrt(
        e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
    for (int i = 0; i < 3; ++i) e2[i] /= e2_norm;
    const double e3[3] = {
        axis[1] * e2[2] - axis[2] * e2[1],
        axis[2] * e2[0] - axis[0] * e2[2],
        axis[0] * e2[1] - axis[1] * e2[0]
    };
    const double phi = 2.0 * Const::pi * uniform01(rng_state);
    const double delta = std::sqrt(delta2);
    const double cp = std::cos(phi);
    const double sp = std::sin(phi);
    for (int i = 0; i < 3; ++i) {
        u_out[i] += delta * (e2[i] * cp + e3[i] * sp);
    }
}

} // namespace

bool nanbu_perez_collide(BackgroundTailPIC& tail, const SpatialGrid& grid,
                         const TailCollisionRequest& request,
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
    const int nxl = grid.nx_local;
    const double inv_dx = 1.0 / grid.dx;

    // Per-cell physical density n_t = sum(w)/dx (m^-3, section 10.3.1).
    // Positions and weights are unchanged by the scattering, so the density
    // is constant over the collision substeps.
    std::vector<double> cell_density(static_cast<size_t>(nxl), 0.0);
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        const int cell = static_cast<int>(
            std::floor(tail.particles[p].x / grid.dx)) - grid.ix_start;
        if (cell >= 0 && cell < nxl) {
            cell_density[static_cast<size_t>(cell)] +=
                tail.particles[p].weight;
        }
    }
    for (int ix = 0; ix < nxl; ++ix) {
        cell_density[static_cast<size_t>(ix)] *= inv_dx;
    }

    // Substep bound: the largest s12 in a cell is bounded by the cold
    // cutoff (gamma_rel/u_rel <= max(1, 1/u_cold) since gamma/u -> 1/beta).
    const double rate_factor =
        std::max(1.0, 1.0 / request.cold_velocity_cutoff);
    const double rate_scale =
        Const::qe * Const::qe * Const::qe * Const::qe * request.coulomb_log /
        (4.0 * Const::pi * Const::eps0 * Const::eps0 * Const::me *
         Const::me * Const::c * Const::c * Const::c);
    double max_rate = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        max_rate = std::max(
            max_rate,
            cell_density[static_cast<size_t>(ix)] * rate_scale *
                rate_factor * rate_factor * rate_factor);
    }
    // The small-angle sampling formula stays accurate up to s ~ 0.5 (the
    // cosine is clamped to [-1,1]); 0.5 keeps the substep count bounded at
    // production density.
    const double s12_substep_target = 0.5;
    const int substeps = std::max(
        1, static_cast<int>(
               std::ceil(max_rate * request.dt / s12_substep_target)));
    if (request.max_substeps <= 0 || substeps > request.max_substeps) {
        diagnostics.collision_substeps = substeps;
        diagnostics.success = false;
        diagnostics.failure_reason =
            TailCollisionFailureReason::SubstepLimit;
        diagnostics.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
        return false;
    }
    const double dt_sub = request.dt / static_cast<double>(substeps);
    diagnostics.collision_substeps = substeps;

    // Diagnostics: totals before.
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
    const std::uint64_t particles_before = tail.particles.size();
    const std::uint64_t growth_budget = static_cast<std::uint64_t>(
        static_cast<double>(particles_before) *
        std::max(0.0, request.max_particle_growth));
    diagnostics.particle_count_limit = particles_before + growth_budget;

    // Collision substeps.  Each substep re-bins the CURRENT particle list
    // (never the main storage order), pairs deterministically and scatters
    // into a fresh list. The bounded unequal-weight correction keeps the
    // original particle IDs and never creates persistent residual particles.
    std::vector<BackgroundTailParticle> current;
    current.swap(tail.particles);
    bool failed = false;
    double pair_events = 0.0;
    for (int sub = 0; sub < substeps && !failed; ++sub) {
        std::vector<BackgroundTailParticle> next;
        next.reserve(current.size());
        std::vector<std::vector<size_t> > cell_bins(
            static_cast<size_t>(nxl));
        for (size_t p = 0; p < current.size(); ++p) {
            const int cell = static_cast<int>(
                std::floor(current[p].x / grid.dx)) - grid.ix_start;
            if (cell >= 0 && cell < nxl) {
                cell_bins[static_cast<size_t>(cell)].push_back(p);
            }
        }
        for (int ix = 0; ix < nxl; ++ix) {
            std::vector<size_t>& bin = cell_bins[static_cast<size_t>(ix)];
            if (bin.empty()) continue;
            const double n_t = cell_density[static_cast<size_t>(ix)];
            // Equal-strata pairs only identical weight layers.  The
            // virtual-split mode may pair arbitrary weights and applies the
            // bounded Sentoku--Kemp correction without creating particles.
            std::vector<size_t> order = bin;
            if (request.weight_mode == TailCollisionWeightMode::EqualStrata) {
                std::stable_sort(order.begin(), order.end(),
                    [&current](size_t lhs, size_t rhs) {
                        if (current[lhs].weight != current[rhs].weight) {
                            return current[lhs].weight < current[rhs].weight;
                        }
                        return current[lhs].id < current[rhs].id;
                    });
                size_t group_begin = 0;
                while (group_begin < order.size()) {
                    size_t group_end = group_begin + 1;
                    while (group_end < order.size() &&
                           current[order[group_end]].weight ==
                               current[order[group_begin]].weight) {
                        ++group_end;
                    }
                    std::uint64_t shuffle_state = mix64(
                        request.rng_seed_base ^
                        (static_cast<std::uint64_t>(request.accepted_step) << 32) ^
                        (static_cast<std::uint64_t>(request.collision_half) << 40) ^
                        (static_cast<std::uint64_t>(sub) << 48) ^
                        static_cast<std::uint64_t>(group_begin));
                    for (size_t i = group_end; i-- > group_begin + 1;) {
                        const size_t span = i - group_begin + 1;
                        const size_t j = group_begin + static_cast<size_t>(
                            uniform01(shuffle_state) *
                            static_cast<double>(span));
                        std::swap(order[i], order[std::min(j, i)]);
                    }
                    group_begin = group_end;
                }
            } else {
                std::uint64_t shuffle_state = mix64(
                    request.rng_seed_base ^
                    (static_cast<std::uint64_t>(request.accepted_step) << 32) ^
                    (static_cast<std::uint64_t>(request.collision_half) << 40) ^
                    (static_cast<std::uint64_t>(sub) << 48));
                for (size_t i = order.size(); i-- > 1;) {
                    const size_t j = static_cast<size_t>(
                        uniform01(shuffle_state) *
                        static_cast<double>(i + 1));
                    std::swap(order[i], order[std::min(j, i)]);
                }
            }
            // Deterministic rotation of the odd leftover (section 10.3.1).
            if (order.size() % 2 != 0) {
                const size_t shift = static_cast<size_t>(sub) % order.size();
                if (shift > 0) {
                    std::rotate(order.begin(), order.begin() + shift,
                                order.end());
                }
            }
            const size_t pair_count = order.size() / 2;
            for (size_t pa = 0; pa < pair_count; ++pa) {
                const size_t i1 = order[2 * pa];
                const size_t i2 = order[2 * pa + 1];
                const BackgroundTailParticle& p1 =
                    current[i1];
                const BackgroundTailParticle& p2 =
                    current[i2];
                if (p1.id == p2.id) {
                    next.push_back(p1);
                    next.push_back(p2);
                    continue;
                }
                const double w1 = p1.weight;
                const double w2 = p2.weight;
                if (w1 <= 0.0 || w2 <= 0.0 ||
                    (request.weight_mode ==
                         TailCollisionWeightMode::EqualStrata &&
                     w1 != w2)) {
                    next.push_back(p1);
                    next.push_back(p2);
                    continue;
                }
                if (request.weight_mode ==
                        TailCollisionWeightMode::VirtualSplit &&
                    w1 != w2) {
                    ++diagnostics.weight_split_count;
                }
                pair_events += 1.0;

                const double g1 = std::sqrt(
                    1.0 + p1.ux * p1.ux + p1.uy * p1.uy + p1.uz * p1.uz);
                const double g2 = std::sqrt(
                    1.0 + p2.ux * p2.ux + p2.uy * p2.uy + p2.uz * p2.uz);
                const double u1[3] = { p1.ux, p1.uy, p1.uz };
                const double u2[3] = { p2.ux, p2.uy, p2.uz };
                // COM boost and relative gamma.
                const double beta[3] = {
                    (u1[0] + u2[0]) / (g1 + g2),
                    (u1[1] + u2[1]) / (g1 + g2),
                    (u1[2] + u2[2]) / (g1 + g2)
                };
                const double beta1[3] = { u1[0] / g1, u1[1] / g1,
                                          u1[2] / g1 };
                const double beta2[3] = { u2[0] / g2, u2[1] / g2,
                                          u2[2] / g2 };
                const double b1dotb2 = beta1[0] * beta2[0] +
                                       beta1[1] * beta2[1] +
                                       beta1[2] * beta2[2];
                const double gamma_rel =
                    std::max(1.0, g1 * g2 * (1.0 - b1dotb2));
                // COM-frame magnitude of each particle's momentum: the
                // reduced gamma of the pair is
                //   gamma_star = sqrt((1 + gamma_rel)/2)
                // so |u*| = sqrt(gamma_star^2 - 1) = sqrt((gamma_rel-1)/2).
                const double u_star = std::sqrt(
                    std::max(0.0, 0.5 * (gamma_rel - 1.0)));
                const double s12 = s12_value(n_t, gamma_rel, u_star, dt_sub,
                                             request.coulomb_log,
                                             request.cold_velocity_cutoff);
                diagnostics.max_s12 = std::max(diagnostics.max_s12, s12);

                // Sample the deflection (small-angle Nanbu limit; substeps
                // keep s12 <= 0.5) and the azimuth.
                std::uint64_t pair_state = mix64(
                    request.rng_seed_base ^
                    (static_cast<std::uint64_t>(request.accepted_step) << 32) ^
                    (static_cast<std::uint64_t>(request.collision_half) << 40) ^
                    (static_cast<std::uint64_t>(sub) << 48) ^
                    (p1.id * 0x9e3779b97f4a7c15ULL) ^
                    (p2.id * 0xbf58476d1ce4e5b9ULL) ^
                    (static_cast<std::uint64_t>(pa) << 16));
                const double R = uniform01(pair_state);
                const double phi = 2.0 * Const::pi *
                                   uniform01(pair_state);
                double cos_theta = 1.0 - 2.0 * s12 * (-std::log(R));
                if (cos_theta < -1.0) cos_theta = -1.0;
                if (cos_theta > 1.0) cos_theta = 1.0;
                if (cos_theta < 0.5) ++diagnostics.large_angle_fraction;
                const double sin_theta = std::sqrt(
                    std::max(0.0, 1.0 - cos_theta * cos_theta));

                // Relative direction in the COM frame: the boost that makes
                // the total spatial momentum vanish is u' = u + ... - gamma
                // gamma_u beta with beta = +beta_com (the algebra fixes the
                // sign: p* = gamma(p - beta E/c) = 0 requires beta = beta_com).
                double u1s[3];
                lorentz_boost(u1, beta, u1s);
                const double ns = std::sqrt(u1s[0] * u1s[0] +
                                            u1s[1] * u1s[1] +
                                            u1s[2] * u1s[2]);
                double n[3] = { 1.0, 0.0, 0.0 };
                if (ns > 1.0e-15) {
                    n[0] = u1s[0] / ns;
                    n[1] = u1s[1] / ns;
                    n[2] = u1s[2] / ns;
                }
                double n_rot[3];
                rotate_direction(n, cos_theta, sin_theta, std::cos(phi),
                                 std::sin(phi), n_rot);
                const double u1s_new[3] = {
                    u_star * n_rot[0], u_star * n_rot[1], u_star * n_rot[2]
                };
                const double u2s_new[3] = {
                    -u_star * n_rot[0], -u_star * n_rot[1],
                    -u_star * n_rot[2]
                };
                double u1_new[3];
                double u2_new[3];
                // Boost back to the lab frame with -beta_com.
                const double neg_beta[3] = { -beta[0], -beta[1], -beta[2] };
                lorentz_boost(u1s_new, neg_beta, u1_new);
                lorentz_boost(u2s_new, neg_beta, u2_new);

                // Unequal-weight collision: retain exactly two persistent
                // macro particles.  The lighter support is fully scattered;
                // the heavier support receives the Sentoku--Kemp weighted
                // correction.  This preserves weighted energy eventwise and
                // weighted momentum in expectation, without the exponential
                // residual-particle growth of materialised virtual shares.
                if (request.weight_mode ==
                        TailCollisionWeightMode::VirtualSplit &&
                    w1 != w2) {
                    if (w1 > w2) {
                        weighted_particle_correction(
                            w2 / w1, u1, u1_new, pair_state, u1_new);
                    } else {
                        weighted_particle_correction(
                            w1 / w2, u2, u2_new, pair_state, u2_new);
                    }
                }
                BackgroundTailParticle n1 = p1;
                n1.weight = w1;
                n1.ux = u1_new[0];
                n1.uy = u1_new[1];
                n1.uz = u1_new[2];
                next.push_back(n1);
                BackgroundTailParticle n2 = p2;
                n2.weight = w2;
                n2.ux = u2_new[0];
                n2.uy = u2_new[1];
                n2.uz = u2_new[2];
                next.push_back(n2);
            }
            // The odd leftover keeps its momentum (no self-scatter).
            if (order.size() % 2 != 0 && !order.empty()) {
                next.push_back(current[order.back()]);
            }
        }
        // Growth budget (section 10.3.3): a soft failure must not silently
        // drop residual weight.
        if (next.size() > particles_before + growth_budget) {
            diagnostics.particle_count_attempted = next.size();
            diagnostics.failure_reason =
                TailCollisionFailureReason::ParticleGrowthBudget;
            failed = true;
            break;
        }
        current.swap(next);
    }
    if (failed) {
        // A soft failure must leave the trial particle list unchanged.
        tail.particles.swap(current);
        diagnostics.success = false;
        return false;
    }
    tail.particles.swap(current);
    diagnostics.particle_count_attempted = tail.particles.size();

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
    diagnostics.large_angle_fraction /= std::max(1.0, pair_events);
    diagnostics.success = true;
    diagnostics.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    return true;
}
