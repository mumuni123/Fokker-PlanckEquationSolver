#include "cylindrical_fp_collision.h"

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <limits>
#include <vector>

namespace {
bool finite_coefficients(const CylindricalCollisionCoefficients& c)
{
    return std::isfinite(c.a_parallel) && std::isfinite(c.a_perp) &&
           std::isfinite(c.d_parallel_parallel) &&
           std::isfinite(c.d_parallel_perp) && std::isfinite(c.d_perp_perp);
}

bool solve_tridiagonal(std::vector<double>& lower, std::vector<double>& diagonal,
                       std::vector<double>& upper, std::vector<double>& rhs,
                       size_t active_size)
{
    // The work arrays are sized to max(Nv,Nmu), but each sweep has its own
    // active length.  Solving the inactive tail contaminates the last active
    // row during back substitution and breaks the conservative face ledger.
    const size_t n = active_size;
    if (n == 0 || n > diagonal.size() || n > lower.size() ||
        n > upper.size() || n > rhs.size()) return false;
    for (size_t i = 1; i < n; ++i) {
        if (!(std::fabs(diagonal[i - 1]) > 0.0)) return false;
        const double ratio = lower[i] / diagonal[i - 1];
        diagonal[i] -= ratio * upper[i - 1];
        rhs[i] -= ratio * rhs[i - 1];
    }
    if (!(std::fabs(diagonal[n - 1]) > 0.0)) return false;
    rhs[n - 1] /= diagonal[n - 1];
    for (size_t i = n - 1; i-- > 0;) {
        if (!(std::fabs(diagonal[i]) > 0.0)) return false;
        rhs[i] = (rhs[i] - upper[i] * rhs[i + 1]) / diagonal[i];
    }
    return true;
}

// 3x3 linear solve with partial pivoting for the moment-closure invariant
// restoration.  The Gram matrix rows live on very different scales (mass,
// momentum, energy), so plain column pivoting with an absolute small-pivot
// test is used; the matrix is well-conditioned when the distribution is
// non-degenerate.
bool solve3x3(const double m[3][3], const double rhs[3],
              double& x0, double& x1, double& x2)
{
    double a[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) a[i][j] = m[i][j];
        a[i][3] = rhs[i];
    }
    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        double best = std::fabs(a[col][col]);
        for (int r = col + 1; r < 3; ++r) {
            const double v = std::fabs(a[r][col]);
            if (v > best) {
                best = v;
                pivot = r;
            }
        }
        if (!(best > 1.0e-300)) return false;
        if (pivot != col) {
            for (int j = 0; j < 4; ++j) std::swap(a[col][j], a[pivot][j]);
        }
        for (int r = col + 1; r < 3; ++r) {
            const double f = a[r][col] / a[col][col];
            if (f == 0.0) continue;
            for (int j = col; j < 4; ++j) a[r][j] -= f * a[col][j];
        }
    }
    x2 = a[2][3] / a[2][2];
    x1 = (a[1][3] - a[1][2] * x2) / a[1][1];
    x0 = (a[0][3] - a[0][1] * x1 - a[0][2] * x2) / a[0][0];
    return true;
}

double chang_cooper_delta(double peclet)
{
    if (std::fabs(peclet) < 1.0e-6) return 0.5 - peclet / 12.0;
    if (peclet > 50.0) return 1.0 / peclet;
    if (peclet < -50.0) return 1.0 + 1.0 / peclet;
    return 1.0 / peclet - 1.0 / std::expm1(peclet);
}

// Conservative explicit cross-diffusion step (d_parallel_perp, section
// 10.2 item 2).  The cross flux
//   F_par  = D_par_perp d f/d u_perp          (u_parallel direction)
//   F_perp = u_perp D_par_perp d f/d u_par    (cylindrical u_perp direction)
// is discretised on the same faces as the Chang-Cooper sweeps with zero
// flux at all domain boundaries, so the step conserves number exactly.
void apply_cross_diffusion_explicit(
    const Species& electrons, const SpatialGrid& grid,
    const LocalCollisionMoments& moments, double time,
    const CollisionCoefficientProvider& provider, double dt,
    int sx, std::vector<double>& stage,
    const std::vector<unsigned char>* bulk_mask,
    const HybridVelocityPartition* partition,
    std::vector<double>* cross_upar_flux,
    std::vector<double>* cross_uperp_flux)
{
    const int nv = Param::Nv;
    const int nmu = Param::Nmu;
    const CylindricalVelocityGrid& cgrid = electrons.cgrid;
    std::vector<double> g(static_cast<size_t>(nv), 0.0);
    std::vector<double> h(static_cast<size_t>(nmu), 0.0);
    std::vector<double> dpar(static_cast<size_t>(nv), 0.0);
    std::vector<double> dper(static_cast<size_t>(nmu), 0.0);
    std::vector<double> fpar(static_cast<size_t>(nv) + 1, 0.0);
    std::vector<double> fper(static_cast<size_t>(nmu) + 1, 0.0);
    const bool has_mask = bulk_mask != NULL;
    const std::vector<unsigned char>& mask = has_mask
        ? *bulk_mask : std::vector<unsigned char>();

    const auto bulk_cell = [&](int iv, int imu) {
        return !has_mask || mask[idx2(iv, imu)] != 0;
    };
    const auto upar_face_flux = [&](int face, int imu, double raw) {
        if (face <= 0 || face >= nv) return 0.0;
        const bool left_bulk = bulk_cell(face - 1, imu);
        const bool right_bulk = bulk_cell(face, imu);
        if (left_bulk && right_bulk) return raw;
        if (partition == NULL) return 0.0;
        const int interface_index =
            partition->upar_interface_index(face, imu);
        if (interface_index < 0) return 0.0;
        const BulkTailInterfaceFace& iface =
            partition->upar_interface_faces[
                static_cast<size_t>(interface_index)];
        if (iface.outward_sign > 0 && left_bulk)
            return std::max(raw, 0.0);
        if (iface.outward_sign < 0 && right_bulk)
            return std::min(raw, 0.0);
        return 0.0;
    };
    const auto uperp_face_flux = [&](int iv, int face, double raw) {
        if (face <= 0 || face >= nmu) return 0.0;
        const bool lower_bulk = bulk_cell(iv, face - 1);
        const bool upper_bulk = bulk_cell(iv, face);
        if (lower_bulk && upper_bulk) return raw;
        if (partition == NULL) return 0.0;
        const int interface_index =
            partition->uperp_interface_index(iv, face);
        if (interface_index < 0) return 0.0;
        const BulkTailInterfaceFace& iface =
            partition->uperp_interface_faces[
                static_cast<size_t>(interface_index)];
        if (iface.outward_sign > 0 && lower_bulk)
            return std::max(raw, 0.0);
        if (iface.outward_sign < 0 && upper_bulk)
            return std::min(raw, 0.0);
        return 0.0;
    };
    for (int imu = 0; imu < nmu; ++imu) {
        // d f/d u_perp (central, one-sided at the ends) and D_par_perp at
        // each u_parallel cell.
        for (int iv = 0; iv < nv; ++iv) {
            const double f0 = stage[idx2(iv, imu)] /
                std::max(grid.dx * cgrid.cell_phase_volume(iv, imu),
                         1.0e-300);
            const double fp = (imu + 1 < nmu)
                                  ? stage[idx2(iv, imu + 1)] /
                                    std::max(grid.dx * cgrid.cell_phase_volume(
                                                 iv, imu + 1), 1.0e-300)
                                  : f0;
            const double fm = (imu > 0)
                                  ? stage[idx2(iv, imu - 1)] /
                                    std::max(grid.dx * cgrid.cell_phase_volume(
                                                 iv, imu - 1), 1.0e-300)
                                  : f0;
            const double du = (imu + 1 < nmu && imu > 0)
                                  ? cgrid.uperp_cells[imu + 1] -
                                        cgrid.uperp_cells[imu - 1]
                                  : (imu + 1 < nmu)
                                        ? cgrid.uperp_cells[imu + 1] -
                                              cgrid.uperp_cells[imu]
                                        : cgrid.uperp_cells[imu] -
                                              cgrid.uperp_cells[imu - 1];
            g[static_cast<size_t>(iv)] = (fp - fm) / std::max(du, 1.0e-300);
            const LocalCollisionMoments m = moments;
            const CylindricalCollisionCoefficients c = provider.evaluate(
                grid.x(sx), cgrid.upar_cells[iv], cgrid.uperp_cells[imu],
                time, m);
            dpar[static_cast<size_t>(iv)] = c.d_parallel_perp;
        }
        fpar[0] = 0.0;
        fpar[static_cast<size_t>(nv)] = 0.0;
        for (int face = 1; face < nv; ++face) {
            const double raw =
                0.5 * (dpar[static_cast<size_t>(face - 1)] *
                           g[static_cast<size_t>(face - 1)] +
                       dpar[static_cast<size_t>(face)] *
                           g[static_cast<size_t>(face)]);
            fpar[static_cast<size_t>(face)] =
                upar_face_flux(face, imu, raw);
        }
        const double area_par =
            grid.dx * cgrid.uperp_ring_areas[imu];
        if (cross_upar_flux != NULL) {
            for (int face = 0; face <= nv; ++face) {
                (*cross_upar_flux)[static_cast<size_t>(face) *
                                   static_cast<size_t>(nmu) +
                                   static_cast<size_t>(imu)] +=
                    dt * area_par * fpar[static_cast<size_t>(face)];
            }
        }
        for (int iv = 0; iv < nv; ++iv) {
            if (!bulk_cell(iv, imu)) continue;
            stage[idx2(iv, imu)] -=
                dt * (fpar[static_cast<size_t>(iv + 1)] -
                      fpar[static_cast<size_t>(iv)]) *
                area_par;
        }
    }

        for (int iv = 0; iv < nv; ++iv) {
        for (int imu = 0; imu < nmu; ++imu) {
            const double f0 = stage[idx2(iv, imu)] /
                std::max(grid.dx * cgrid.cell_phase_volume(iv, imu),
                         1.0e-300);
            const double fp = (iv + 1 < nv)
                                  ? stage[idx2(iv + 1, imu)] /
                                    std::max(grid.dx * cgrid.cell_phase_volume(
                                                 iv + 1, imu), 1.0e-300)
                                  : f0;
            const double fm = (iv > 0)
                                  ? stage[idx2(iv - 1, imu)] /
                                    std::max(grid.dx * cgrid.cell_phase_volume(
                                                 iv - 1, imu), 1.0e-300)
                                  : f0;
            const double du = (iv + 1 < nv && iv > 0)
                                  ? cgrid.upar_cells[iv + 1] -
                                        cgrid.upar_cells[iv - 1]
                                  : (iv + 1 < nv)
                                        ? cgrid.upar_cells[iv + 1] -
                                              cgrid.upar_cells[iv]
                                        : cgrid.upar_cells[iv] -
                                              cgrid.upar_cells[iv - 1];
            h[static_cast<size_t>(imu)] = (fp - fm) / std::max(du, 1.0e-300);
            const LocalCollisionMoments m = moments;
            const CylindricalCollisionCoefficients c = provider.evaluate(
                grid.x(sx), cgrid.upar_cells[iv], cgrid.uperp_cells[imu],
                time, m);
            dper[static_cast<size_t>(imu)] = c.d_parallel_perp;
        }
        fper[0] = 0.0;
        fper[static_cast<size_t>(nmu)] = 0.0;
        for (int face = 1; face < nmu; ++face) {
            const double raw = cgrid.uperp_faces[face] *
                0.5 * (dper[static_cast<size_t>(face - 1)] *
                           h[static_cast<size_t>(face - 1)] +
                       dper[static_cast<size_t>(face)] *
                           h[static_cast<size_t>(face)]);
            fper[static_cast<size_t>(face)] =
                uperp_face_flux(iv, face, raw);
        }
        const double area_perp =
            grid.dx * 2.0 * Const::pi * cgrid.upar_widths[iv];
        if (cross_uperp_flux != NULL) {
            for (int face = 0; face <= nmu; ++face) {
                (*cross_uperp_flux)[static_cast<size_t>(iv) *
                                    static_cast<size_t>(nmu + 1) +
                                    static_cast<size_t>(face)] +=
                    dt * area_perp * fper[static_cast<size_t>(face)];
            }
        }
        for (int imu = 0; imu < nmu; ++imu) {
            if (!bulk_cell(iv, imu)) continue;
            stage[idx2(iv, imu)] -=
                dt * (fper[static_cast<size_t>(imu + 1)] -
                      fper[static_cast<size_t>(imu)]) *
                area_perp;
        }
    }
}
} // namespace

CylindricalFokkerPlanckCollision::CylindricalFokkerPlanckCollision(
    const CollisionCoefficientProvider& provider, CollisionIntegratorType type)
    : provider_(provider), type_(type),
      bulk_integrator_(BulkCollisionIntegrator::BGK_VALIDATION)
{}

bool CylindricalFokkerPlanckCollision::is_trivial() const
{
    return provider_.mode() == CollisionCoefficientMode::NONE;
}

CollisionDiagnostics CylindricalFokkerPlanckCollision::apply(
    Species& electrons, const SpatialGrid& grid, double time, double dt,
    const std::vector<unsigned char>* bulk_mask) const
{
    CollisionFaceFluxes unused;
    return apply_with_flux(electrons, grid, time, dt, bulk_mask, NULL,
                           unused);
}

CollisionDiagnostics CylindricalFokkerPlanckCollision::apply_with_flux(
    Species& electrons, const SpatialGrid& grid, double time, double dt,
    const std::vector<unsigned char>* bulk_mask,
    const HybridVelocityPartition* partition,
    CollisionFaceFluxes& fluxes) const
{
    CollisionDiagnostics result;
    const bool export_flux = partition != NULL &&
        bulk_integrator_ == BulkCollisionIntegrator::CHANG_COOPER_FLUX;
    // The flux path is a public transaction boundary: the caller may use the
    // returned face batch to create tail particles, so a late residual or
    // positivity failure must not leave a partially updated bulk state.  The
    // snapshot is only taken for the exporting path; the legacy validation
    // path already runs on the integrator's collision trial buffer.
    const std::vector<double> state_before_flux = export_flux
        ? electrons.f : std::vector<double>();
    const auto rollback_flux_transaction = [&]() {
        if (!export_flux) return;
        electrons.f = state_before_flux;
        fluxes.upar_flux.clear();
        fluxes.uperp_flux.clear();
        fluxes.cross_upar_flux.clear();
        fluxes.cross_uperp_flux.clear();
        fluxes.nx_local = 0;
        fluxes.nv = 0;
        fluxes.nmu = 0;
        result.transaction_rollback_count = 1;
    };
    if (export_flux) {
        fluxes.resize(grid.nx_local, Param::Nv, Param::Nmu);
    }
    if (!(dt > 0.0)) return result;
    if (is_trivial()) return result;

    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int nv = Param::Nv;
    const int nmu = Param::Nmu;
    const CylindricalVelocityGrid& cgrid = electrons.cgrid;

    const bool has_mask = bulk_mask != NULL;
    static const std::vector<unsigned char> no_mask;
    const std::vector<unsigned char>& mask =
        has_mask ? *bulk_mask : no_mask;

    // Per-cell local moments (section 10.2.1 mode 2 input): number, mean
    // u_parallel and kinetic-energy density.  Scanned once per collision
    // half-step; the provider is then evaluated at every velocity point with
    // these cell moments.
    std::vector<LocalCollisionMoments> cell_moments(
        static_cast<size_t>(nxl));
    int invalid = 0;
    int has_cross = 0;
    double max_rate = 0.0;
    #pragma omp parallel for reduction(|:invalid,has_cross) reduction(max:max_rate)
    for (int ix = 0; ix < nxl; ++ix) {
        const int sx = ng + ix;
        const size_t xbase = static_cast<size_t>(sx) *
                             static_cast<size_t>(Param::Nvmu);
        double n = 0.0;
        double u_mean = 0.0;
        double ke = 0.0;
            for (int iv = 0; iv < nv; ++iv) {
                const size_t row = xbase + static_cast<size_t>(iv) * nmu;
                const size_t vrow = static_cast<size_t>(iv) * nmu;
                for (int imu = 0; imu < nmu; ++imu) {
                    const double m = electrons.f[row + imu];
                    if (has_mask && mask[idx2(iv, imu)] == 0) {
                        continue;
                    }
                    // The physical distribution is non-negative; negative
                    // cells are remap artifacts and must not drive the
                    // closure moments.
                    if (m < 0.0) continue;
                    n += m;
                u_mean += m * cgrid.upar_cells[iv];
                ke += m * cgrid.kinetic_energy[vrow + imu];
            }
        }
        LocalCollisionMoments& moments = cell_moments[static_cast<size_t>(ix)];
        moments.density = n / grid.dx;
        moments.u_parallel_mean =
            (n > 0.0) ? u_mean / std::max(n, 1.0e-300) : 0.0;
        moments.kinetic_energy_density = ke / grid.dx;
        // Section 19.3 energy-conservation fix: make the moment-closure
        // temperature self-consistent with the discrete grid moments, so
        // the operator equilibrium is exactly the distribution with the
        // current temperature (otherwise the closure temperature chases the
        // undercounted quadrature moments and the bulk energy drains
        // monotonically; the H9 no-beam 40fs run cooled 39%).
        if (provider_.mode() == CollisionCoefficientMode::MOMENT_CLOSURE) {
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
                    cgrid.moment_closure_uth2_self_consistent(uth2_mom);
            }
        }

        // Validity + rate scan at the velocity cell centres.
        for (int iv = 0; iv < nv; ++iv) {
            for (int imu = 0; imu < nmu; ++imu) {
                if (has_mask &&
                    mask[idx2(iv, imu)] == 0) {
                    continue;
                }
                const CylindricalCollisionCoefficients c = provider_.evaluate(
                    grid.x(sx), cgrid.upar_cells[iv], cgrid.uperp_cells[imu],
                    time, moments);
                if (!finite_coefficients(c) ||
                    c.d_parallel_parallel < 0.0 || c.d_perp_perp < 0.0) {
                    invalid = 1;
                }
                if (std::fabs(c.d_parallel_perp) > 0.0) has_cross = 1;
                max_rate = std::max(
                    max_rate, std::fabs(c.a_parallel) + std::fabs(c.a_perp) +
                                  c.d_parallel_parallel + c.d_perp_perp);
            }
        }
    }
    if (invalid) {
        result.success = false;
        rollback_flux_transaction();
        return result;
    }
    if (max_rate == 0.0) return result;
    // Stage H7 (section 10.3.5): the physical e-e collision frequency at
    // production density can make max_rate*dt >> 1 in one Strang half-step.
    // Subcycle the implicit sweeps so every collision substep stays in the
    // Chang-Cooper positivity/stability regime; the coefficients (cell
    // moments) are fixed at the half-step start so bulk and tail use the
    // same time layer (section 10.4).
    const double substep_target = 0.1;
    const int collision_substeps = std::max(
        1, static_cast<int>(std::ceil(max_rate * dt / substep_target)));
    const double dt_sub = dt / static_cast<double>(collision_substeps);

    double local_mass_before = 0.0;
    double local_mass_after = 0.0;
    double local_mass_change = 0.0;
    double local_ke_before = 0.0;
    double local_ke_after = 0.0;
    double local_px_before = 0.0;
    double local_px_after = 0.0;
    double local_flux_residual_linf = 0.0;
    double local_flux_divergence_sum = 0.0;
    double local_interface_inward_clipped = 0.0;
    int solve_failed = 0;
    int first_failure_reported = 0;
    #pragma omp parallel reduction(+:local_mass_before,local_mass_after, \
                          local_mass_change, \
                         local_ke_before,local_ke_after, \
                                   local_px_before,local_px_after) \
                           reduction(max:local_flux_residual_linf) \
                           reduction(+:local_flux_divergence_sum) \
                          reduction(+:local_interface_inward_clipped) \
                          reduction(|:solve_failed)
    {
        std::vector<double> lower(std::max(nv, nmu), 0.0);
        std::vector<double> diagonal(std::max(nv, nmu), 0.0);
        std::vector<double> upper(std::max(nv, nmu), 0.0);
        std::vector<double> rhs(std::max(nv, nmu), 0.0);
        std::vector<double> stage(static_cast<size_t>(Param::Nvmu), 0.0);
        std::vector<double> cross_upar_local(
            static_cast<size_t>(nv + 1) * static_cast<size_t>(nmu), 0.0);
        std::vector<double> cross_uperp_local(
            static_cast<size_t>(nv) * static_cast<size_t>(nmu + 1), 0.0);
        #pragma omp for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            if (export_flux) {
                std::fill(cross_upar_local.begin(), cross_upar_local.end(), 0.0);
                std::fill(cross_uperp_local.begin(), cross_uperp_local.end(), 0.0);
            }
            const int sx = ng + ix;
            const LocalCollisionMoments& moments =
                cell_moments[static_cast<size_t>(ix)];
            double n0 = 0.0;
            double px_u0 = 0.0;
            double ke0 = 0.0;
            double fmax0 = 0.0;
            for (int iv = 0; iv < nv; ++iv) {
                for (int imu = 0; imu < nmu; ++imu) {
                    const double mass = electrons.f[idx3(sx, iv, imu)];
                    // Tail-owned (excluded) cells are not part of the bulk
                    // collision: keep them unchanged and out of the bulk
                    // moments and diagnostics.
                    if (has_mask && mask[idx2(iv, imu)] == 0) {
                        stage[idx2(iv, imu)] = mass;
                        continue;
                    }
                    // Keep the input as-is: the BGK relaxation below is
                    // positivity-preserving for non-negative input and
                    // damps (never amplifies) small numerical negatives
                    // from the remap, so no clipping is needed.  Negative
                    // cells are remap artifacts: they stay in stage for the
                    // BGK to damp but are excluded from the bulk moments,
                    // the invariant restoration targets and the
                    // diagnostics (consistent with the projection and the
                    // moments scan).
                    stage[idx2(iv, imu)] = mass;
                    // The conservation ledger is an audit of the actual
                    // state, so retain finite negative remap debt here.  The
                    // positivity gate below decides whether that debt is
                    // admissible; silently omitting it would manufacture a
                    // mass loss and make the face-flux closure appear broken.
                    local_mass_before += mass;
                    if (mass < 0.0) continue;
                    n0 += mass;
                    fmax0 = std::max(fmax0, mass);
                    px_u0 += mass * cgrid.upar_cells[iv];
                    ke0 += mass * cgrid.kinetic_energy[idx2(iv, imu)];
                    local_ke_before +=
                        mass * cgrid.kinetic_energy[idx2(iv, imu)];
                    local_px_before +=
                        mass * Const::me * Const::c *
                        cgrid.upar_cells[iv];
                }
            }

            // The remap can leave signed roundoff in exponentially small
            // Maxwellian tails.  Chang--Cooper legitimately redistributes
            // that inherited debt, but it is neither a physical negative
            // distribution nor a collision instability.  Remove only this
            // local, scale-relative residue from the collision trial; a
            // resolvable negative remains intact and still fails below.
            const double inherited_negative_floor =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max(fmax0, n0 / std::max(
                    1.0, static_cast<double>(Param::Nvmu)));
            if (inherited_negative_floor > 0.0) {
                for (int iv = 0; iv < nv; ++iv) {
                    for (int imu = 0; imu < nmu; ++imu) {
                        if (has_mask && mask[idx2(iv, imu)] == 0) continue;
                        double& value = stage[idx2(iv, imu)];
                        if (value < 0.0 &&
                            -value <= inherited_negative_floor) {
                            value = 0.0;
                        }
                    }
                }
            }

            if (provider_.mode() == CollisionCoefficientMode::MOMENT_CLOSURE &&
                bulk_integrator_ == BulkCollisionIntegrator::BGK_VALIDATION) {
                // Section 19.3: the moment-closure mode uses an implicit
                // BGK-style relaxation toward the local drifting Maxwellian
                // instead of the ADI Chang-Cooper sweeps.  The sweeps
                // oscillate and lose positivity on strongly non-Maxwellian
                // states (H9 beam-40fs failed with negative cells), while
                // the implicit relaxation is unconditionally positive,
                // oscillation-free and exact for the relaxation model:
                //   f_new = (f + dt*nu(u)*M) / (1 + dt*nu(u)).
                // M uses the self-consistent closure temperature
                // (u_th2_closure, section 19.3) and the cell moments.  The
                // exact number/momentum/energy invariants are restored by
                // the projection below.
                const double n_cell = moments.density * grid.dx;
                const double u_mean = moments.u_parallel_mean;
                const double random_ke =
                    moments.kinetic_energy_density /
                        std::max(moments.density, 1.0e-300) -
                    0.5 * Const::me * Const::c * Const::c * u_mean * u_mean;
                const double uth2 =
                    (moments.u_th2_closure > 0.0)
                        ? moments.u_th2_closure
                        : (2.0 / 3.0) * random_ke /
                              (Const::me * Const::c * Const::c);
                if (n_cell > 0.0 && uth2 > 0.0 && random_ke > 0.0) {
                    const double u_th = std::sqrt(uth2);
                    const double nu0 =
                        moments.density * Const::qe * Const::qe *
                        Const::qe * Const::qe *
                        20.0 /* Coulomb logarithm, matches the moment-closure
                                provider */ /
                        (4.0 * Const::pi * Const::eps0 * Const::eps0 *
                         Const::me * Const::me * Const::c * Const::c *
                         Const::c * uth2 * u_th);
                    double z = 0.0;
                    for (int iv = 0; iv < nv; ++iv) {
                        const double du = cgrid.upar_cells[iv] - u_mean;
                        for (int imu = 0; imu < nmu; ++imu) {
                            if (has_mask &&
                                mask[idx2(iv, imu)] == 0) continue;
                            const double u2 = du * du +
                                cgrid.uperp_cells[imu] *
                                    cgrid.uperp_cells[imu];
                            z += std::exp(-u2 * (0.5 / uth2)) *
                                 cgrid.cell_phase_volume(iv, imu);
                        }
                    }
                    if (z > 0.0) {
                        const double m_norm = n_cell / z;
                        for (int iv = 0; iv < nv; ++iv) {
                            const double du = cgrid.upar_cells[iv] - u_mean;
                            for (int imu = 0; imu < nmu; ++imu) {
                                if (has_mask &&
                                    mask[idx2(iv, imu)] == 0) continue;
                                const double up = cgrid.uperp_cells[imu];
                                const double u2 = du * du + up * up;
                                const double rate =
                                    nu0 / std::pow(1.0 + u2 / uth2, 1.5);
                                const double m_eq =
                                    m_norm * std::exp(-u2 * (0.5 / uth2)) *
                                    cgrid.cell_phase_volume(iv, imu);
                                stage[idx2(iv, imu)] =
                                    (stage[idx2(iv, imu)] + dt * rate * m_eq) /
                                    (1.0 + dt * rate);
                            }
                        }
                    }
                }
            } else {
                // Collision substeps: explicit conservative cross-diffusion
                // half, then the two implicit diagonal Chang-Cooper sweeps
                // with the velocity-point coefficients evaluated from the
                // fixed half-step cell moments.
                for (int sub = 0; sub < collision_substeps; ++sub) {
                    if (has_cross) {
                        apply_cross_diffusion_explicit(
                            electrons, grid, moments, time, provider_,
                            dt_sub, sx, stage,
                            has_mask ? &mask : NULL, partition,
                            export_flux ? &cross_upar_local : NULL,
                            export_flux ? &cross_uperp_local : NULL);
                    }

                    // Backward-Euler u_parallel sweep, F=0 at both outer
                    // faces.
                    for (int imu = 0; imu < nmu; ++imu) {
                    std::fill(lower.begin(), lower.begin() + nv, 0.0);
                    std::fill(diagonal.begin(), diagonal.begin() + nv, 1.0);
                    std::fill(upper.begin(), upper.begin() + nv, 0.0);
                    for (int iv = 0; iv < nv; ++iv) {
                        rhs[iv] = stage[idx2(iv, imu)];
                    }
                    for (int face = 1; face < nv; ++face) {
                        const int l = face - 1;
                        const int r = face;
                        const bool l_bulk = !has_mask ||
                            mask[idx2(l, imu)] != 0;
                        const bool r_bulk = !has_mask ||
                            mask[idx2(r, imu)] != 0;
                        const BulkTailInterfaceFace* interface_face = NULL;
                        if (export_flux && partition != NULL &&
                            static_cast<size_t>(face) *
                                    static_cast<size_t>(nmu) +
                                static_cast<size_t>(imu) <
                                partition->upar_interface_lookup.size()) {
                            const int interface_index =
                                partition->upar_interface_lookup[
                                    static_cast<size_t>(face) *
                                    static_cast<size_t>(nmu) +
                                    static_cast<size_t>(imu)];
                            if (interface_index >= 0) {
                                interface_face = &partition->upar_interface_faces[
                                    static_cast<size_t>(interface_index)];
                            }
                        }
                        if (!l_bulk || !r_bulk) {
                            if (interface_face == NULL) continue;
                        }
                        const double distance =
                            cgrid.upar_cells[r] - cgrid.upar_cells[l];
                        const double area =
                            grid.dx * cgrid.uperp_ring_areas[imu];
                        // Equilibrium-preserving evaluation point: the
                        // midpoint of the two adjacent cell centres, not the
                        // geometric face.  On the non-uniform production
                        // u_parallel grid (sinh stretch) the geometric face
                        // is not the midpoint of the cell centres, so the
                        // Chang-Cooper Peclet pe = a*distance/d would not
                        // equal the discrete Maxwellian log-ratio and the
                        // equilibrium would not be stationary, causing a
                        // spurious bulk energy drift.  The moment-closure
                        // a/d = -(u-u_mean)/u_th^2 is linear in u, so the
                        // midpoint makes pe exactly the equilibrium
                        // log-ratio on any grid.
                        const double u_mid =
                            0.5 * (cgrid.upar_cells[l] +
                                   cgrid.upar_cells[r]);
                        const CylindricalCollisionCoefficients c =
                            provider_.evaluate(
                                grid.x(sx), u_mid,
                                cgrid.uperp_cells[imu], time, moments);
                        const double a = c.a_parallel;
                        const double d = c.d_parallel_parallel;
                        const double pe =
                            d > 0.0 ? a * distance / d : 0.0;
                        // Upwind weight on the left cell.  The classic
                        // Chang-Cooper delta is 1/pe - 1/(e^pe-1); with the
                        // flux written as
                        //   F = a[(1-delta)f_l + delta f_r] - d(f_r-f_l)/du
                        // the equilibrium-preserving (positivity-preserving)
                        // left weight is the complement 1 - delta_cc.
                        const double upwind =
                            d > 0.0 ? chang_cooper_delta(pe)
                                    : (a >= 0.0 ? 0.0 : 1.0);
                        const double wl = 1.0 - upwind;
                        const double vl =
                            grid.dx * cgrid.cell_phase_volume(l, imu);
                        const double vr =
                            grid.dx * cgrid.cell_phase_volume(r, imu);
                        const double coeff_l =
                            area * (a * wl + d / distance) / vl;
                        const double coeff_r =
                            area * (a * (1.0 - wl) - d / distance) / vr;
                        if (l_bulk && r_bulk) {
                            diagonal[l] += dt_sub * coeff_l;
                            upper[l] += dt_sub * coeff_r;
                            lower[r] -= dt_sub * coeff_l;
                            diagonal[r] -= dt_sub * coeff_r;
                        } else if (interface_face != NULL) {
                            if (interface_face->outward_sign > 0 && l_bulk) {
                                // The tail-side density is fixed to zero and
                                // only the outward part is admitted.
                                if (coeff_l > 0.0)
                                    diagonal[l] += dt_sub * coeff_l;
                            } else if (interface_face->outward_sign < 0 &&
                                       r_bulk) {
                                if (coeff_r < 0.0)
                                    diagonal[r] -= dt_sub * coeff_r;
                            }
                        }
                    }
                        if (!solve_tridiagonal(
                                lower, diagonal, upper, rhs,
                                static_cast<size_t>(nv))) {
                            solve_failed = 1;
                            continue;
                        }
                        for (int iv = 0; iv < nv; ++iv) {
                            stage[idx2(iv, imu)] = rhs[iv];
                        }
                        if (export_flux) {
                            for (int face = 1; face < nv; ++face) {
                                const int l = face - 1;
                                const int r = face;
                                const bool l_bulk = !has_mask ||
                                    mask[idx2(l, imu)] != 0;
                                const bool r_bulk = !has_mask ||
                                    mask[idx2(r, imu)] != 0;
                                const BulkTailInterfaceFace* iface = NULL;
                                if (partition != NULL) {
                                    const int ii = partition->upar_interface_lookup[
                                        static_cast<size_t>(face) *
                                        static_cast<size_t>(nmu) +
                                        static_cast<size_t>(imu)];
                                    if (ii >= 0) iface =
                                        &partition->upar_interface_faces[
                                            static_cast<size_t>(ii)];
                                }
                                if (!l_bulk && !r_bulk) continue;
                                if (!l_bulk || !r_bulk) {
                                    if (iface == NULL) continue;
                                    const double distance =
                                        cgrid.upar_cells[r] -
                                        cgrid.upar_cells[l];
                                    const double area =
                                        grid.dx * cgrid.uperp_ring_areas[imu];
                                    const double u_mid =
                                        0.5 * (cgrid.upar_cells[l] +
                                               cgrid.upar_cells[r]);
                                    const CylindricalCollisionCoefficients c =
                                        provider_.evaluate(
                                            grid.x(sx), u_mid,
                                            cgrid.uperp_cells[imu], time,
                                            moments);
                                    const double pe = c.d_parallel_parallel > 0.0
                                        ? c.a_parallel * distance /
                                          c.d_parallel_parallel : 0.0;
                                    const double wl = c.d_parallel_parallel > 0.0
                                        ? 1.0 - chang_cooper_delta(pe)
                                        : (c.a_parallel >= 0.0 ? 1.0 : 0.0);
                                    const double vl = grid.dx *
                                        cgrid.cell_phase_volume(l, imu);
                                    const double vr = grid.dx *
                                        cgrid.cell_phase_volume(r, imu);
                                    double transfer = 0.0;
                                    if (iface->outward_sign > 0 && l_bulk) {
                                        const double coeff = area *
                                            (c.a_parallel * wl +
                                             c.d_parallel_parallel / distance) / vl;
                                        transfer = dt_sub * coeff *
                                            stage[idx2(l, imu)];
                                        if (transfer < 0.0) transfer = 0.0;
                                    } else if (iface->outward_sign < 0 && r_bulk) {
                                        const double coeff = area *
                                            (c.a_parallel * (1.0 - wl) -
                                             c.d_parallel_parallel / distance) / vr;
                                        transfer = dt_sub * coeff *
                                            stage[idx2(r, imu)];
                                        if (transfer > 0.0) {
                                            local_interface_inward_clipped +=
                                                transfer;
                                            transfer = 0.0;
                                        }
                                    }
                                    fluxes.upar_flux[fluxes.upar_index(
                                        ix, face, imu)] += transfer;
                                } else {
                                    const double distance =
                                        cgrid.upar_cells[r] -
                                        cgrid.upar_cells[l];
                                    const double area =
                                        grid.dx * cgrid.uperp_ring_areas[imu];
                                    const double u_mid =
                                        0.5 * (cgrid.upar_cells[l] +
                                               cgrid.upar_cells[r]);
                                    const CylindricalCollisionCoefficients c =
                                        provider_.evaluate(
                                            grid.x(sx), u_mid,
                                            cgrid.uperp_cells[imu], time,
                                            moments);
                                    const double pe = c.d_parallel_parallel > 0.0
                                        ? c.a_parallel * distance /
                                          c.d_parallel_parallel : 0.0;
                                    const double wl = c.d_parallel_parallel > 0.0
                                        ? 1.0 - chang_cooper_delta(pe)
                                        : (c.a_parallel >= 0.0 ? 1.0 : 0.0);
                                    const double vl = grid.dx *
                                        cgrid.cell_phase_volume(l, imu);
                                    const double vr = grid.dx *
                                        cgrid.cell_phase_volume(r, imu);
                                    const double f_l = stage[idx2(l, imu)] /
                                        std::max(vl, 1.0e-300);
                                    const double f_r = stage[idx2(r, imu)] /
                                        std::max(vr, 1.0e-300);
                                    const double face_flux = area *
                                        (c.a_parallel * (wl * f_l +
                                                         (1.0 - wl) * f_r) -
                                         c.d_parallel_parallel *
                                         (f_r - f_l) / distance);
                                    fluxes.upar_flux[fluxes.upar_index(
                                        ix, face, imu)] += dt_sub * face_flux;
                                }
                            }
                        }
                    }

                    // Backward-Euler u_perp sweep; axis and outer boundary
                    // have zero normal collision flux.
                    for (int iv = 0; iv < nv; ++iv) {
                    std::fill(lower.begin(), lower.begin() + nmu, 0.0);
                    std::fill(diagonal.begin(), diagonal.begin() + nmu, 1.0);
                    std::fill(upper.begin(), upper.begin() + nmu, 0.0);
                    for (int imu = 0; imu < nmu; ++imu) {
                        rhs[imu] = stage[idx2(iv, imu)];
                    }
                    for (int face = 1; face < nmu; ++face) {
                        const int l = face - 1;
                        const int r = face;
                        const bool l_bulk = !has_mask ||
                            mask[idx2(iv, l)] != 0;
                        const bool r_bulk = !has_mask ||
                            mask[idx2(iv, r)] != 0;
                        const BulkTailInterfaceFace* interface_face = NULL;
                        if (export_flux && partition != NULL &&
                            static_cast<size_t>(iv) *
                                    static_cast<size_t>(nmu + 1) +
                                static_cast<size_t>(face) <
                                partition->uperp_interface_lookup.size()) {
                            const int interface_index =
                                partition->uperp_interface_lookup[
                                    static_cast<size_t>(iv) *
                                    static_cast<size_t>(nmu + 1) +
                                    static_cast<size_t>(face)];
                            if (interface_index >= 0) {
                                interface_face = &partition->uperp_interface_faces[
                                    static_cast<size_t>(interface_index)];
                            }
                        }
                        if (!l_bulk || !r_bulk) {
                            if (interface_face == NULL) continue;
                        }
                        const double distance =
                            cgrid.uperp_cells[r] - cgrid.uperp_cells[l];
                        const double radius = cgrid.uperp_faces[face];
                        const double area =
                            grid.dx * 2.0 * Const::pi * radius *
                            cgrid.upar_widths[iv];
                        // Same equilibrium-preserving midpoint rule as the
                        // u_parallel sweep (the u_perp grid is also a
                        // non-uniform sinh stretch; the geometric face is
                        // not the midpoint of the adjacent cell centres).
                        const double uperp_mid =
                            0.5 * (cgrid.uperp_cells[l] +
                                   cgrid.uperp_cells[r]);
                        const CylindricalCollisionCoefficients c =
                            provider_.evaluate(
                                grid.x(sx), cgrid.upar_cells[iv],
                                uperp_mid, time, moments);
                        const double a = c.a_perp;
                        const double d = c.d_perp_perp;
                        const double pe =
                            d > 0.0 ? a * distance / d : 0.0;
                        // Same equilibrium-preserving upwind complement as
                        // the u_parallel sweep.
                        const double upwind =
                            d > 0.0 ? chang_cooper_delta(pe)
                                    : (a >= 0.0 ? 0.0 : 1.0);
                        const double wl = 1.0 - upwind;
                        const double vl =
                            grid.dx * cgrid.cell_phase_volume(iv, l);
                        const double vr =
                            grid.dx * cgrid.cell_phase_volume(iv, r);
                        const double coeff_l =
                            area * (a * wl + d / distance) / vl;
                        const double coeff_r =
                            area * (a * (1.0 - wl) - d / distance) / vr;
                        if (l_bulk && r_bulk) {
                            diagonal[l] += dt_sub * coeff_l;
                            upper[l] += dt_sub * coeff_r;
                            lower[r] -= dt_sub * coeff_l;
                            diagonal[r] -= dt_sub * coeff_r;
                        } else if (interface_face != NULL) {
                            if (interface_face->outward_sign > 0 && l_bulk) {
                                if (coeff_l > 0.0)
                                    diagonal[l] += dt_sub * coeff_l;
                            } else if (interface_face->outward_sign < 0 &&
                                       r_bulk) {
                                if (coeff_r < 0.0)
                                    diagonal[r] -= dt_sub * coeff_r;
                            }
                        }
                    }
                        if (!solve_tridiagonal(
                                lower, diagonal, upper, rhs,
                                static_cast<size_t>(nmu))) {
                            solve_failed = 1;
                            continue;
                        }
                        for (int imu = 0; imu < nmu; ++imu) {
                            stage[idx2(iv, imu)] = rhs[imu];
                        }
                        if (export_flux) {
                            for (int face = 1; face < nmu; ++face) {
                                const int l = face - 1;
                                const int r = face;
                                const bool l_bulk = !has_mask ||
                                    mask[idx2(iv, l)] != 0;
                                const bool r_bulk = !has_mask ||
                                    mask[idx2(iv, r)] != 0;
                                const BulkTailInterfaceFace* iface = NULL;
                                if (partition != NULL) {
                                    const int ii = partition->uperp_interface_lookup[
                                        static_cast<size_t>(iv) *
                                        static_cast<size_t>(nmu + 1) +
                                        static_cast<size_t>(face)];
                                    if (ii >= 0) iface =
                                        &partition->uperp_interface_faces[
                                            static_cast<size_t>(ii)];
                                }
                                if (!l_bulk && !r_bulk) continue;
                                const double distance =
                                    cgrid.uperp_cells[r] -
                                    cgrid.uperp_cells[l];
                                const double radius = cgrid.uperp_faces[face];
                                const double area =
                                    grid.dx * 2.0 * Const::pi * radius *
                                    cgrid.upar_widths[iv];
                                const double uperp_mid =
                                    0.5 * (cgrid.uperp_cells[l] +
                                           cgrid.uperp_cells[r]);
                                const CylindricalCollisionCoefficients c =
                                    provider_.evaluate(
                                        grid.x(sx), cgrid.upar_cells[iv],
                                        uperp_mid, time, moments);
                                const double pe = c.d_perp_perp > 0.0
                                    ? c.a_perp * distance / c.d_perp_perp : 0.0;
                                const double wl = c.d_perp_perp > 0.0
                                    ? 1.0 - chang_cooper_delta(pe)
                                    : (c.a_perp >= 0.0 ? 1.0 : 0.0);
                                const double vl = grid.dx *
                                    cgrid.cell_phase_volume(iv, l);
                                const double vr = grid.dx *
                                    cgrid.cell_phase_volume(iv, r);
                                double transfer = 0.0;
                                if (l_bulk && r_bulk) {
                                    const double f_l = stage[idx2(iv, l)] /
                                        std::max(vl, 1.0e-300);
                                    const double f_r = stage[idx2(iv, r)] /
                                        std::max(vr, 1.0e-300);
                                    const double face_flux = area *
                                        (c.a_perp * (wl * f_l +
                                                      (1.0 - wl) * f_r) -
                                         c.d_perp_perp *
                                         (f_r - f_l) / distance);
                                    transfer = dt_sub * face_flux;
                                } else if (iface != NULL &&
                                           iface->outward_sign > 0 && l_bulk) {
                                    const double coeff = area *
                                        (c.a_perp * wl +
                                         c.d_perp_perp / distance) / vl;
                                    transfer = dt_sub * coeff * stage[idx2(iv, l)];
                                    if (transfer < 0.0) transfer = 0.0;
                                } else if (iface != NULL &&
                                           iface->outward_sign < 0 && r_bulk) {
                                    const double coeff = area *
                                        (c.a_perp * (1.0 - wl) -
                                         c.d_perp_perp / distance) / vr;
                                    transfer = dt_sub * coeff * stage[idx2(iv, r)];
                                    if (transfer > 0.0) {
                                        local_interface_inward_clipped +=
                                            transfer;
                                        transfer = 0.0;
                                    }
                                }
                                fluxes.uperp_flux[fluxes.uperp_index(
                                    ix, iv, face)] += transfer;
                            }
                        }
                    }
                }
            }
            if (export_flux) {
                for (size_t q = 0; q < cross_upar_local.size(); ++q) {
                    fluxes.cross_upar_flux[
                        static_cast<size_t>(ix) * cross_upar_local.size() + q] +=
                        cross_upar_local[q];
                }
                for (size_t q = 0; q < cross_uperp_local.size(); ++q) {
                    fluxes.cross_uperp_flux[
                        static_cast<size_t>(ix) * cross_uperp_local.size() + q] +=
                        cross_uperp_local[q];
                }
            }
            // Moment-closure invariant restoration (section 19.3): the
            // drift-diffusion closure is not energy-conserving for
            // non-Maxwellian states (H9 beam-40fs showed a large negative
            // reservoir feeding the instability).  Restore the cell number,
            // parallel momentum and kinetic energy with a minimal
            // multiplicative weight w = 1 + a + b*u + c*K applied only to
            // cells with physical mass (excludes numerical far-tail noise).
            if (provider_.mode() == CollisionCoefficientMode::MOMENT_CLOSURE &&
                bulk_integrator_ == BulkCollisionIntegrator::BGK_VALIDATION) {
                double fmax_cell = 0.0;
                for (int iv = 0; iv < nv; ++iv) {
                    for (int imu = 0; imu < nmu; ++imu) {
                        if (has_mask && mask[idx2(iv, imu)] == 0) continue;
                        if (stage[idx2(iv, imu)] < 0.0) continue;
                        fmax_cell = std::max(
                            fmax_cell, stage[idx2(iv, imu)]);
                    }
                }
                if (fmax_cell > 0.0) {
                    const double fmin_include = 1.0e-30 * fmax_cell;
                    double msum = 0.0;
                    double usum = 0.0;
                    double ksum = 0.0;
                    double uusum = 0.0;
                    double uksum = 0.0;
                    double kksum = 0.0;
                    for (int iv = 0; iv < nv; ++iv) {
                        const double u = cgrid.upar_cells[iv];
                        for (int imu = 0; imu < nmu; ++imu) {
                            if (has_mask && mask[idx2(iv, imu)] == 0)
                                continue;
                            const double m = stage[idx2(iv, imu)];
                            if (m < 0.0 || m < fmin_include) continue;
                            const double k =
                                cgrid.kinetic_energy[idx2(iv, imu)];
                            msum += m;
                            usum += m * u;
                            ksum += m * k;
                            uusum += m * u * u;
                            uksum += m * u * k;
                            kksum += m * k * k;
                        }
                    }
                    if (msum > 0.0) {
                        const double mat[3][3] = {
                            { msum, usum, ksum },
                            { usum, uusum, uksum },
                            { ksum, uksum, kksum }
                        };
                        const double rhs[3] = {
                            n0 - msum, px_u0 - usum, ke0 - ksum
                        };
                        double a = 0.0, b = 0.0, c = 0.0;
                        if (solve3x3(mat, rhs, a, b, c)) {
                            for (int iv = 0; iv < nv; ++iv) {
                                const double u = cgrid.upar_cells[iv];
                                for (int imu = 0; imu < nmu; ++imu) {
                                    if (has_mask &&
                                        mask[idx2(iv, imu)] == 0) continue;
                                    if (stage[idx2(iv, imu)] < 0.0 ||
                                        stage[idx2(iv, imu)] < fmin_include)
                                        continue;
                                    const double k =
                                        cgrid.kinetic_energy[idx2(iv, imu)];
                                    const double w =
                                        1.0 + a + b * u + c * k;
                                    // The multiplicative weight must not
                                    // flip a positive cell negative: a
                                    // negative w at a high-energy cell would
                                    // create an unphysical negative (H9
                                    // beam run, u~12.3 cell).  Such cells
                                    // carry negligible mass, so skipping
                                    // them keeps the invariant restoration
                                    // exact to ~1e-30 relative while
                                    // preserving positivity.
                                    if (w <= 0.0) continue;
                                    stage[idx2(iv, imu)] *= w;
                                }
                            }
                        }
                    }
                }
            }
            for (int iv = 0; iv < nv; ++iv) {
                for (int imu = 0; imu < nmu; ++imu) {
                    if (has_mask && mask[idx2(iv, imu)] == 0) continue;
                    const double value = stage[idx2(iv, imu)];
                    const double original =
                        electrons.f[idx3(sx, iv, imu)];
                    // The positivity gate fires only when the collision
                    // itself turned a non-negative cell negative (or
                    // amplified an existing negative).  Pre-existing tiny
                    // negatives from the remap are not the collision's
                    // fault: the BGK relaxation only damps them (never
                    // amplifies), and the projection skips them.
                    const bool collision_made_negative =
                        original >= 0.0 &&
                        value < -1.0e-13 *
                                    std::max(1.0, std::fabs(original));
                    const bool amplified_negative =
                        original < -inherited_negative_floor &&
                        value < original * (1.0 + 1.0e-10);
                    if (!std::isfinite(value) ||
                        collision_made_negative || amplified_negative) {
                        int print_this_failure = 0;
                        #pragma omp atomic capture
                        { print_this_failure = first_failure_reported;
                          first_failure_reported = 1; }
                        if (print_this_failure == 0) {
                            std::fprintf(stderr,
                                "[collision-invalid] cell ix=%d iv=%d imu=%d "
                                "value=%.6e original=%.6e u=%.6e up=%.6e\n",
                                ix, iv, imu, value, original,
                                cgrid.upar_cells[iv],
                                cgrid.uperp_cells[imu]);
                        }
                        solve_failed = 1;
                    }
                    if (export_flux) {
                        const double up_right = fluxes.upar_flux[
                            fluxes.upar_index(ix, iv + 1, imu)] +
                            fluxes.cross_upar_flux[
                            fluxes.upar_index(ix, iv + 1, imu)];
                        const double up_left = fluxes.upar_flux[
                            fluxes.upar_index(ix, iv, imu)] +
                            fluxes.cross_upar_flux[
                            fluxes.upar_index(ix, iv, imu)];
                        const double pp_right = fluxes.uperp_flux[
                            fluxes.uperp_index(ix, iv, imu + 1)] +
                            fluxes.cross_uperp_flux[
                            fluxes.uperp_index(ix, iv, imu + 1)];
                        const double pp_left = fluxes.uperp_flux[
                            fluxes.uperp_index(ix, iv, imu)] +
                            fluxes.cross_uperp_flux[
                            fluxes.uperp_index(ix, iv, imu)];
                        const double residual = value - original +
                            up_right - up_left + pp_right - pp_left;
                        local_flux_residual_linf = std::max(
                            local_flux_residual_linf, std::fabs(residual));
                        local_flux_divergence_sum +=
                            up_right - up_left + pp_right - pp_left;
                    }
                    electrons.f[idx3(sx, iv, imu)] = value;
                    // Accumulate the change before the global reduction.  It
                    // avoids subtracting two O(N) totals and is therefore
                    // stable under OpenMP/MPI repartitioning.
                    local_mass_change += value - original;
                    local_mass_after += value;
                    if (value < 0.0) continue;
                    local_ke_after +=
                        value * cgrid.kinetic_energy[idx2(iv, imu)];
                    local_px_after +=
                        value * Const::me * Const::c * cgrid.upar_cells[iv];
                }
            }
        }
    }
    if (type_ == CollisionIntegratorType::TR_BDF2) result.order_reduced = true;
    double sums[6] = {
        local_mass_before, local_mass_after, local_ke_before, local_ke_after,
        local_px_before, local_px_after
    };
    double global[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(sums, global, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_mass_change = 0.0;
    MPI_Allreduce(&local_mass_change, &global_mass_change, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    int global_failed = 0;
    MPI_Allreduce(&solve_failed, &global_failed, 1, MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);
    result.mass_change = global_mass_change;
    // Section 10.2 item 6 / 11.2: energy and momentum exchanged with the
    // field (collision reservoir) are tracked explicitly.  Positive means
    // the bulk distribution lost energy/momentum to the reservoir.
    result.reservoir_energy_change = global[2] - global[3];
    result.reservoir_momentum_change = global[4] - global[5];
    double global_flux_residual_linf = 0.0;
    double global_cross_pair_residual_linf = 0.0;
    double global_flux_divergence_sum = 0.0;
    MPI_Allreduce(&local_flux_residual_linf, &global_flux_residual_linf, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (export_flux && partition != NULL) {
        // Every cross-diffusion face is stored once.  Summing its conservative
        // divergence over bulk-owned cells must therefore equal only the
        // cross flux that leaves through a bulk/tail interface; internal
        // bulk/bulk faces cancel algebraically.  This is an audit of the
        // shared-face topology, not a positivity or physics correction.
        double local_cross_divergence = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            for (int iv = 0; iv < nv; ++iv) {
                for (int imu = 0; imu < nmu; ++imu) {
                    if (has_mask && mask[idx2(iv, imu)] == 0) continue;
                    local_cross_divergence +=
                        fluxes.cross_upar_flux[fluxes.upar_index(
                            ix, iv + 1, imu)] -
                        fluxes.cross_upar_flux[fluxes.upar_index(
                            ix, iv, imu)] +
                        fluxes.cross_uperp_flux[fluxes.uperp_index(
                            ix, iv, imu + 1)] -
                        fluxes.cross_uperp_flux[fluxes.uperp_index(
                            ix, iv, imu)];
                }
            }
        }
        double local_cross_interface = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            for (size_t q = 0; q < partition->upar_interface_faces.size(); ++q) {
                const BulkTailInterfaceFace& face =
                    partition->upar_interface_faces[q];
                local_cross_interface +=
                    static_cast<double>(face.outward_sign) *
                    fluxes.cross_upar_flux[fluxes.upar_index(
                        ix, face.face_index, face.transverse_index)];
            }
            for (size_t q = 0; q < partition->uperp_interface_faces.size(); ++q) {
                const BulkTailInterfaceFace& face =
                    partition->uperp_interface_faces[q];
                local_cross_interface +=
                    static_cast<double>(face.outward_sign) *
                    fluxes.cross_uperp_flux[fluxes.uperp_index(
                        ix, face.transverse_index, face.face_index)];
            }
        }
        double cross_pair_local[2] = {local_cross_divergence,
                                      local_cross_interface};
        double cross_pair_global[2] = {0.0, 0.0};
        MPI_Allreduce(cross_pair_local, cross_pair_global, 2, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        global_cross_pair_residual_linf = std::fabs(
            cross_pair_global[0] - cross_pair_global[1]);
    }
    MPI_Allreduce(&local_flux_divergence_sum, &global_flux_divergence_sum, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double global_inward_clipped = 0.0;
    MPI_Allreduce(&local_interface_inward_clipped, &global_inward_clipped, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    result.implicit_flux_residual_linf = global_flux_residual_linf;
    result.cross_flux_pair_residual_linf = global_cross_pair_residual_linf;
    result.flux_divergence_sum = global_flux_divergence_sum;
    result.mass_flux_balance_residual = result.mass_change +
        result.flux_divergence_sum;
    result.interface_inward_clipped_number = global_inward_clipped;

    // Export a direct diagnostic from the final conservative face fluxes.
    // The integrator later replaces the number/energy with the exact
    // quadrature moments of the committed collision parcels; this fallback
    // keeps the standalone collision operator auditable and, importantly,
    // never reconstructs an interface source from a cell-state difference.
    if (export_flux && partition != NULL) {
        double local_interface_number = 0.0;
        double local_interface_energy = 0.0;
        double local_interface_faces = 0.0;
        const double mec2 = Const::me * Const::c * Const::c;
        for (int ix = 0; ix < nxl; ++ix) {
            for (size_t q = 0; q < partition->upar_interface_faces.size(); ++q) {
                const BulkTailInterfaceFace& face =
                    partition->upar_interface_faces[q];
                const double transfer =
                    fluxes.upar_flux[fluxes.upar_index(
                        ix, face.face_index, face.transverse_index)] +
                    fluxes.cross_upar_flux[fluxes.upar_index(
                        ix, face.face_index, face.transverse_index)];
                const double amount = face.outward_sign > 0
                    ? std::max(transfer, 0.0)
                    : std::max(-transfer, 0.0);
                if (amount <= 0.0) continue;
                const double up = cgrid.upar_faces[face.face_index];
                const double ut = cgrid.uperp_cells[face.transverse_index];
                local_interface_number += amount;
                local_interface_energy += amount * mec2 *
                    (std::sqrt(1.0 + up * up + ut * ut) - 1.0);
                local_interface_faces += 1.0;
            }
            for (size_t q = 0; q < partition->uperp_interface_faces.size(); ++q) {
                const BulkTailInterfaceFace& face =
                    partition->uperp_interface_faces[q];
                const double transfer =
                    fluxes.uperp_flux[fluxes.uperp_index(
                        ix, face.transverse_index, face.face_index)] +
                    fluxes.cross_uperp_flux[fluxes.uperp_index(
                        ix, face.transverse_index, face.face_index)];
                const double amount = face.outward_sign > 0
                    ? std::max(transfer, 0.0)
                    : std::max(-transfer, 0.0);
                if (amount <= 0.0) continue;
                const double up = cgrid.upar_cells[face.transverse_index];
                const double ut = cgrid.uperp_faces[face.face_index];
                local_interface_number += amount;
                local_interface_energy += amount * mec2 *
                    (std::sqrt(1.0 + up * up + ut * ut) - 1.0);
                local_interface_faces += 1.0;
            }
        }
        const double local_export[3] = {local_interface_number,
                                        local_interface_energy,
                                        local_interface_faces};
        double global_export[3] = {0.0, 0.0, 0.0};
        MPI_Allreduce(local_export, global_export, 3, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.interface_export_number = global_export[0];
        result.interface_export_energy = global_export[1];
        result.interface_parcel_count = static_cast<std::size_t>(
            std::llround(global_export[2]));
    }
    const double flux_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(global[0]),
                               std::fabs(global[1])));
    result.success = global_failed == 0 &&
        (!export_flux || (global_flux_residual_linf <= flux_tolerance &&
                          global_cross_pair_residual_linf <= flux_tolerance));
    if (!result.success) rollback_flux_transaction();
    return result;
}
