#include "hybrid_collision_step.h"

#include "background_tail_collision_sde.h"
#include "background_tail_nanbu_perez.h"

#include <algorithm>
#include <cmath>
#include <vector>

// Stage-H8 local field-particle reaction (section 10.7): the tail's actual
// random momentum/energy change in a cell (from C_tb) must appear as the
// opposite change in the bulk.  The bulk correction is deposited on a local
// three-cell u_parallel stencil (j, j+1, j+2 near the thermal peak) as a
// number-preserving linear source.  The weights are solved for the two
// constraints (P and K) and must stay nonnegative; otherwise the step fails
// (section 10.7 item 5: reduce substep or transactional failure).
namespace {

bool apply_bulk_reaction(Species& bulk, const SpatialGrid& grid, int cell,
                         double dp, double dk, double& applied_dp,
                         double& applied_dk)
{
    applied_dp = 0.0;
    applied_dk = 0.0;
    if (dp == 0.0 && dk == 0.0) return true;
    const int ng = grid.nghost;
    const int nv = Param::Nv;
    // Pick the (u_parallel, u_perp) cell with the largest local mass (the
    // folded Maxwellian vanishes on the axis, so the peak is in k>=1) and
    // use a three-cell stencil in u_parallel at that u_perp ring.
    int j_peak = 0;
    int k_peak = 0;
    double mass_peak = -1.0;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double m = bulk.f[idx3(ng + cell, j, k)];
            if (m > mass_peak) {
                mass_peak = m;
                j_peak = j;
                k_peak = k;
            }
        }
    }
    if (!(mass_peak > 0.0)) return false;
    const int j0 = std::max(0, std::min(nv - 3, j_peak));
    // Stencil cells (j0, j0+1, j0+2) with weights w0, w1, w2 that change
    // P and K while conserving N (w0+w1+w2 = 0).  Solve the 2x2 system for
    // (w1, w2) with w0 = -w1-w2:
    //   u1*w1 + u2*w2 + u0*w0 = dp / (m_e c)
    //   K1*w1 + K2*w2 + K0*w0 = dk
    // (in u-space: the per-cell mass weights are in m^-2; the reaction
    // source is added to f).
    const double u0 = bulk.cgrid.upar_cells[j0];
    const double u1 = bulk.cgrid.upar_cells[j0 + 1];
    const double u2 = bulk.cgrid.upar_cells[j0 + 2];
    const double ke0 = bulk.cgrid.kinetic_energy[idx2(j0, k_peak)];
    const double ke1 = bulk.cgrid.kinetic_energy[idx2(j0 + 1, k_peak)];
    const double ke2 = bulk.cgrid.kinetic_energy[idx2(j0 + 2, k_peak)];
    const double dp_u = dp / (Const::me * Const::c);
    const double a11 = u1 - u0;
    const double a12 = u2 - u0;
    const double a21 = ke1 - ke0;
    const double a22 = ke2 - ke0;
    const double det = a11 * a22 - a12 * a21;
    if (!(std::fabs(det) > 0.0)) return false;
    const double w1 = (dp_u * a22 - a12 * dk) / det;
    const double w2 = (a11 * dk - dp_u * a21) / det;
    const double w0 = -w1 - w2;
    // Nonnegativity: the correction must not drive any stencil cell below
    // zero (section 10.7 item 5: positive mass budget).
    const double m0 = bulk.f[idx3(ng + cell, j0, k_peak)];
    const double m1 = bulk.f[idx3(ng + cell, j0 + 1, k_peak)];
    const double m2 = bulk.f[idx3(ng + cell, j0 + 2, k_peak)];
    if (m0 + w0 < 0.0 || m1 + w1 < 0.0 || m2 + w2 < 0.0) {
        return false;
    }
    bulk.f[idx3(ng + cell, j0, k_peak)] = m0 + w0;
    bulk.f[idx3(ng + cell, j0 + 1, k_peak)] = m1 + w1;
    bulk.f[idx3(ng + cell, j0 + 2, k_peak)] = m2 + w2;
    applied_dp = Const::me * Const::c *
        (u0 * w0 + u1 * w1 + u2 * w2);
    applied_dk = ke0 * w0 + ke1 * w1 + ke2 * w2;
    return true;
}

} // namespace

HybridCollisionDiagnostics HybridCollisionStep::advance(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid, const HybridCollisionConfig& config)
{
    HybridCollisionDiagnostics result;
    if (config.bulk_provider == NULL || !(config.dt > 0.0)) {
        result.success = false;
        result.failure_reason = HybridCollisionFailureReason::InvalidConfig;
        return result;
    }
    // Snapshot the trials so a failed half-step leaves them unchanged.
    Species bulk_before = bulk_trial;
    BackgroundTailPIC tail_before = tail_trial;

    // 1) bulk--bulk FP update (the H7 operator) with the shared provider.
    if (config.pairs.bulk_bulk) {
        CylindricalFokkerPlanckCollision bulk_collision(
            *config.bulk_provider, CollisionIntegratorType::BACKWARD_EULER);
        bulk_collision.set_bulk_integrator(config.bulk_integrator);
        if (config.bulk_face_fluxes != NULL &&
            config.bulk_integrator == BulkCollisionIntegrator::CHANG_COOPER_FLUX) {
            result.bulk_diag = bulk_collision.apply_with_flux(
                bulk_trial, grid, config.time, config.dt,
                config.bulk_velocity_mask, config.bulk_partition,
                *config.bulk_face_fluxes);
        } else {
            result.bulk_diag = bulk_collision.apply(
                bulk_trial, grid, config.time, config.dt,
                config.bulk_velocity_mask);
        }
        if (!result.bulk_diag.success) {
            bulk_trial = bulk_before;
            tail_trial = tail_before;
            result.success = false;
            result.failure_reason = HybridCollisionFailureReason::BulkBulk;
            return result;
        }
        result.bulk_bulk_applied = true;
    } else {
        result.bulk_diag.success = true;
    }

    // Per-cell bulk moments for the tail--bulk coefficients: taken from the
    // PRE-step bulk state so bulk and tail use the same time layer as the
    // bulk FP update (section 10.4: never bulk half-step moments with tail
    // step-start moments).
    std::vector<LocalCollisionMoments> cell_moments(
        static_cast<size_t>(grid.nx_local));
    {
        const int ng = grid.nghost;
        const int nv = Param::Nv;
        const int nmu = Param::Nmu;
        for (int ix = 0; ix < grid.nx_local; ++ix) {
            const int sx = ng + ix;
            const size_t xbase = static_cast<size_t>(sx) *
                                 static_cast<size_t>(Param::Nvmu);
            double n = 0.0;
            double u_mean = 0.0;
            double ke = 0.0;
            for (int j = 0; j < nv; ++j) {
                const size_t row = xbase + static_cast<size_t>(j) * nmu;
                const size_t vrow = static_cast<size_t>(j) * nmu;
                for (int k = 0; k < nmu; ++k) {
                    if (config.bulk_velocity_mask != NULL &&
                        (*config.bulk_velocity_mask)[
                            static_cast<size_t>(j) * nmu +
                            static_cast<size_t>(k)] == 0) {
                        continue;
                    }
                    const double m = bulk_before.f[row + k];
                    n += m;
                    u_mean += m * bulk_before.cgrid.upar_cells[j];
                    ke += m * bulk_before.cgrid.kinetic_energy[vrow + k];
                }
            }
            LocalCollisionMoments& moments = cell_moments[static_cast<size_t>(ix)];
            moments.density = n / grid.dx;
            moments.u_parallel_mean =
                (n > 0.0) ? u_mean / std::max(n, 1.0e-300) : 0.0;
            moments.kinetic_energy_density = ke / grid.dx;
            // Section 19.3: same grid-consistent closure temperature as the
            // bulk FP scan, so the tail--bulk SDE coefficients use the same
            // equilibrium temperature as the bulk operator.
            if (config.bulk_provider->mode() ==
                CollisionCoefficientMode::MOMENT_CLOSURE) {
                const double random_ke =
                    moments.kinetic_energy_density /
                        std::max(moments.density, 1.0e-300) -
                    0.5 * Const::me * Const::c * Const::c *
                        moments.u_parallel_mean * moments.u_parallel_mean;
                if (random_ke > 0.0) {
                    const double uth2_mom =
                        (2.0 / 3.0) * random_ke /
                        (Const::me * Const::c * Const::c);
                    moments.u_th2_closure =
                        bulk_before.cgrid.moment_closure_uth2_self_consistent(
                            uth2_mom);
                }
            }
        }
    }

    // Shared request fields. Pair selection and backend selection are
    // intentionally independent: one physical half-step may invoke both
    // Nanbu--Perez for tail--tail and SDE for tail--bulk.
    TailCollisionRequest tail_request;
    tail_request.dt = config.dt;
    tail_request.accepted_step = config.accepted_step;
    tail_request.collision_half = config.collision_half;
    tail_request.coulomb_log = config.coulomb_log;
    tail_request.rng_seed_base = config.rng_seed_base;
    tail_request.mpi_rank = config.mpi_rank;
    tail_request.weight_mode = config.weight_mode;
    tail_request.max_substeps = config.max_substeps;
    tail_request.max_particle_growth = config.max_particle_growth;
    tail_request.pairs = config.pairs;

    // 2) tail--tail backend.
    if (config.pairs.tail_tail &&
        config.tail_tail_kernel != TailCollisionKernel::None) {
        tail_request.kernel = config.tail_tail_kernel;
        if (!nanbu_perez_collide(tail_trial, grid, tail_request,
                                 result.tail_tail_diag)) {
            bulk_trial = bulk_before;
            tail_trial = tail_before;
            result.success = false;
            result.failure_reason = HybridCollisionFailureReason::TailTail;
            return result;
        }
        result.tail_tail_applied = true;
        result.tail_tail_px_change =
            result.tail_tail_diag.px_after - result.tail_tail_diag.px_before;
        result.tail_tail_energy_change =
            result.tail_tail_diag.ke_after - result.tail_tail_diag.ke_before;
    }

    // 3) tail--bulk SDE (C_tb) with the per-cell field-particle reaction
    // (section 10.7): the tail's actual momentum/energy change in each cell
    // is balanced by the opposite bulk reaction.
    if (config.pairs.tail_bulk &&
        config.tail_bulk_kernel != TailCollisionKernel::None) {
        const BackgroundTailPIC tail_before_tail_bulk = tail_trial;
        double tail_px_before = 0.0;
        double tail_ke_before = 0.0;
        for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
            const BackgroundTailParticle& pp = tail_trial.particles[p];
            tail_px_before += Const::me * Const::c * pp.weight * pp.ux;
            const double gamma = std::sqrt(
                1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
            tail_ke_before += Const::me * Const::c * Const::c * pp.weight *
                              (gamma - 1.0);
        }
        tail_request.kernel = config.tail_bulk_kernel;
        if (!sde_collide(tail_trial, grid, tail_request,
                         *config.bulk_provider, cell_moments,
                         result.tail_bulk_diag)) {
            bulk_trial = bulk_before;
            tail_trial = tail_before;
            result.success = false;
            result.failure_reason = HybridCollisionFailureReason::TailBulk;
            return result;
        }
        result.tail_bulk_applied = true;
        double tail_px_after = 0.0;
        double tail_ke_after = 0.0;
        for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
            const BackgroundTailParticle& pp = tail_trial.particles[p];
            tail_px_after += Const::me * Const::c * pp.weight * pp.ux;
            const double gamma = std::sqrt(
                1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
            tail_ke_after += Const::me * Const::c * Const::c * pp.weight *
                             (gamma - 1.0);
        }
        const double tail_dp = tail_px_after - tail_px_before;
        const double tail_dk = tail_ke_after - tail_ke_before;
        result.tail_px_change = tail_dp;
        result.tail_energy_change = tail_dk;
        // C_tb may be run without C_bt only for an explicitly requested
        // trace/test-particle experiment. Production Coulomb enables both.
        if (!config.pairs.bulk_tail) {
            result.tail_balance_error = 1.0;
        } else {
        double bulk_reaction_dp = 0.0;
        double bulk_reaction_dk = 0.0;
        std::vector<double> tail_px_cell_before(
            static_cast<size_t>(grid.nx_local), 0.0);
        std::vector<double> tail_ke_cell_before(
            static_cast<size_t>(grid.nx_local), 0.0);
        std::vector<double> tail_px_cell_after(
            static_cast<size_t>(grid.nx_local), 0.0);
        std::vector<double> tail_ke_cell_after(
            static_cast<size_t>(grid.nx_local), 0.0);
        for (size_t p = 0; p < tail_before_tail_bulk.particles.size(); ++p) {
            const BackgroundTailParticle& pp =
                tail_before_tail_bulk.particles[p];
            const int cell = static_cast<int>(std::floor(pp.x / grid.dx)) -
                             grid.ix_start;
            if (cell < 0 || cell >= grid.nx_local) continue;
            const double gamma = std::sqrt(
                1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
            tail_px_cell_before[static_cast<size_t>(cell)] +=
                Const::me * Const::c * pp.weight * pp.ux;
            tail_ke_cell_before[static_cast<size_t>(cell)] +=
                Const::me * Const::c * Const::c * pp.weight * (gamma - 1.0);
        }
        for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
            const BackgroundTailParticle& pp = tail_trial.particles[p];
            const int cell = static_cast<int>(std::floor(pp.x / grid.dx)) -
                             grid.ix_start;
            if (cell < 0 || cell >= grid.nx_local) continue;
            const double gamma = std::sqrt(
                1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
            tail_px_cell_after[static_cast<size_t>(cell)] +=
                Const::me * Const::c * pp.weight * pp.ux;
            tail_ke_cell_after[static_cast<size_t>(cell)] +=
                Const::me * Const::c * Const::c * pp.weight * (gamma - 1.0);
        }
        for (int ix = 0; ix < grid.nx_local; ++ix) {
            const double delta_p =
                tail_px_cell_after[static_cast<size_t>(ix)] -
                tail_px_cell_before[static_cast<size_t>(ix)];
            const double delta_k =
                tail_ke_cell_after[static_cast<size_t>(ix)] -
                tail_ke_cell_before[static_cast<size_t>(ix)];
            if (delta_p == 0.0 && delta_k == 0.0) continue;
            double applied_dp = 0.0;
            double applied_dk = 0.0;
            if (!apply_bulk_reaction(bulk_trial, grid, ix, -delta_p,
                                     -delta_k, applied_dp, applied_dk)) {
                bulk_trial = bulk_before;
                tail_trial = tail_before;
                result.success = false;
                result.failure_reason =
                    HybridCollisionFailureReason::BulkReaction;
                result.failure_cell = ix;
                return result;
            }
            bulk_reaction_dp += applied_dp;
            bulk_reaction_dk += applied_dk;
            ++result.reaction_cells;
        }
        const double p_scale = std::max(
            1.0e-300, std::fabs(bulk_reaction_dp) + std::fabs(tail_dp));
        const double k_scale = std::max(
            1.0e-300, std::fabs(bulk_reaction_dk) + std::fabs(tail_dk));
        result.reaction_px_balance =
            std::fabs(bulk_reaction_dp + tail_dp) / p_scale;
        result.reaction_energy_balance =
            std::fabs(bulk_reaction_dk + tail_dk) / k_scale;
        result.bulk_reaction_px_change = bulk_reaction_dp;
        result.bulk_reaction_energy_change = bulk_reaction_dk;
        result.bulk_reaction_applied = result.reaction_cells > 0;
        }
    }

    // Combined ledger.
    // Combined ledger: the total bulk+tail energy change over the half-step
    // (equals the bulk FP reservoir exchange when the reaction balances).
    const double bulk_ke_before = bulk_before.total_kinetic_energy();
    const double bulk_ke_after = bulk_trial.total_kinetic_energy();
    double tail_ke_before = 0.0;
    double tail_ke_after = 0.0;
    for (size_t p = 0; p < tail_before.particles.size(); ++p) {
        const BackgroundTailParticle& pp = tail_before.particles[p];
        const double gamma = std::sqrt(
            1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
        tail_ke_before += Const::me * Const::c * Const::c * pp.weight *
                          (gamma - 1.0);
    }
    for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
        const BackgroundTailParticle& pp = tail_trial.particles[p];
        const double gamma = std::sqrt(
            1.0 + pp.ux * pp.ux + pp.uy * pp.uy + pp.uz * pp.uz);
        tail_ke_after += Const::me * Const::c * Const::c * pp.weight *
                         (gamma - 1.0);
    }
    result.combined_energy_change =
        (bulk_ke_after - bulk_ke_before) + (tail_ke_after - tail_ke_before);
    result.bulk_energy_change = bulk_ke_after - bulk_ke_before;
    result.bulk_bulk_energy_change =
        result.bulk_energy_change - result.bulk_reaction_energy_change;
    if (result.tail_bulk_applied && !config.pairs.bulk_tail) {
        result.tail_balance_error = 1.0;
    } else {
        result.tail_balance_error = result.reaction_energy_balance;
    }
    result.bulk_mass_change =
        bulk_trial.total_particle_number() - bulk_before.total_particle_number();
    result.success = true;
    return result;
}
