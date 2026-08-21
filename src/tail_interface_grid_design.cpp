#include "tail_interface_grid_design.h"

#include "tail_moment_constraint.h"
#include "tail_subcell_quadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace {

// Deterministic integration resolution: every baseline cell is subdivided
// into this many intervals for the monitor cumulative integral.  The faces
// are then inverted on the fine cumulative mesh, so the generated grid is a
// genuine equidistribution of the monitor rather than a copy of the
// baseline faces.
const int kMonitorSubintervalsPerCell = 512;

// "Low u_perp conversion cells" for the u* centroid: bins whose u_perp cell
// centre is within this many first-cell widths of the axis.  4 matches the
// upper bound of the sigma_perp scan (1, 2, 4 current low-u_perp cell
// widths); the rule is deterministic and documented, never a hard-coded iv.
const double kLowUperpCells = 4.0;

// Important cells for the important-cell maximum-relative-difference gate:
// mass >= 1e-6 of the total remapped mass.
const double kImportantCellFraction = 1.0e-6;

// Sparse support budget per cell (section 7.11.16C item 3).
const int kMaxSparseSupport = 7;

// The centre quartet creates 4 macro-particles per conversion cell.
const double kCenterQuartetPerCell = 4.0;

// Nominal representative axial field for the velocity-dt budget.  The
// absolute value cancels in dt_candidate/dt_current; it is only written so
// the per-constraint estimates are transparent.
const double kReferenceFieldVm = 1.0e12;

// Production collision coefficients: moment-closure Coulomb log and the
// subcycle target used by cylindrical_fp_collision.cpp.
const double kCoulombLog = 20.0;
const double kCollisionSubstepTarget = 0.1;

const double kTiny = 1.0e-300;

std::array<double, BULK_TAIL_MOMENT_COUNT> unit_moments(double upar,
                                                        double uperp)
{
    std::array<double, BULK_TAIL_MOMENT_COUNT> r = {};
    mass_cell_moments(1.0, upar, uperp, r[0], r[1], r[3], r[2], r[4], r[5]);
    return r;
}

void add_moments(std::array<double, BULK_TAIL_MOMENT_COUNT>& dst,
                 const std::array<double, BULK_TAIL_MOMENT_COUNT>& src,
                 double scale)
{
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) dst[m] += scale * src[m];
}

double relative_ratio(double a, double b)
{
    // r(a,b): a/b when b>0; 0 when both are zero; +inf when b==0 and a>0
    // (the offline prototype reports the undefined case as a metric, the
    // hard gates only compare defined ratios).
    if (b > 0.0) return a / b;
    if (a == 0.0) return 0.0;
    return std::numeric_limits<double>::infinity();
}

bool face_arrays_identical(const std::vector<double>& lhs,
                           const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (size_t i = 0; i < lhs.size(); ++i)
        if (lhs[i] != rhs[i]) return false;
    return true;
}

// Deterministic monitor equidistribution on one non-negative half of a
// symmetric domain.  baseline_half_faces are the baseline faces from the
// axis (index 0) to the outer endpoint (index n_cells).  The equidistribution
// runs in the baseline grid coordinate xi (baseline faces sit at uniform
// xi = j/n_cells), with the monitor evaluated at the physical midpoint of
// every fine subinterval.  Equal monitor mass in xi is inverted back to u
// through the baseline faces, so w == 1 reproduces the baseline faces
// exactly (Gx only moves faces near the threshold, never the whole grid),
// and w > 1 concentrates cells where the conversion threshold sits.  The
// result has n_cells+1 faces with faces[0]=0 and faces[n_cells]=outer
// endpoint exactly.
std::vector<double> invert_monitor_half(
    const std::vector<double>& baseline_half_faces, int n_cells,
    double monitor_amp, double monitor_sigma, double monitor_center)
{
    std::vector<double> faces(static_cast<size_t>(n_cells) + 1, 0.0);
    if (n_cells <= 0 || baseline_half_faces.size() !=
                            static_cast<size_t>(n_cells) + 1) {
        return faces;
    }
    // A zero amplitude must reproduce the baseline bit-for-bit (the
    // candidate grid is then the current grid in that direction).
    if (monitor_amp == 0.0) return baseline_half_faces;
    const int n_sub = n_cells * kMonitorSubintervalsPerCell;
    std::vector<double> cum(static_cast<size_t>(n_sub) + 1, 0.0);
    std::vector<double> xi_edge(static_cast<size_t>(n_sub) + 1, 0.0);
    size_t slot = 1;
    for (int j = 0; j < n_cells; ++j) {
        const double u_lo = baseline_half_faces[static_cast<size_t>(j)];
        const double u_hi = baseline_half_faces[static_cast<size_t>(j) + 1];
        const double dxi = 1.0 / static_cast<double>(n_cells);
        for (int s = 0; s < kMonitorSubintervalsPerCell; ++s) {
            const double xi_lo = (static_cast<double>(j) +
                                  static_cast<double>(s) /
                                      kMonitorSubintervalsPerCell) /
                                 static_cast<double>(n_cells);
            const double xi_hi = (static_cast<double>(j) +
                                  static_cast<double>(s + 1) /
                                      kMonitorSubintervalsPerCell) /
                                 static_cast<double>(n_cells);
            const double xi_mid = 0.5 * (xi_lo + xi_hi);
            // Physical midpoint inside the baseline cell (linear in xi).
            const double u_mid = u_lo + (u_hi - u_lo) *
                                             (xi_mid -
                                              static_cast<double>(j) /
                                                  n_cells) *
                                             static_cast<double>(n_cells);
            const double delta = u_mid - monitor_center;
            const double w =
                1.0 + monitor_amp *
                          std::exp(-(delta * delta) /
                                   std::max(monitor_sigma * monitor_sigma,
                                            kTiny));
            cum[slot] = cum[slot - 1] + w * dxi /
                                            kMonitorSubintervalsPerCell;
            xi_edge[slot] = xi_edge[slot - 1] +
                            dxi / kMonitorSubintervalsPerCell;
            ++slot;
        }
    }
    const double total = cum[n_sub];
    if (!(total > 0.0) || !std::isfinite(total)) return faces;
    faces[0] = baseline_half_faces[0];
    faces[static_cast<size_t>(n_cells)] =
        baseline_half_faces[static_cast<size_t>(n_cells)];
    for (int m = 1; m < n_cells; ++m) {
        const double target = total * static_cast<double>(m) / n_cells;
        // Binary search for the fine xi interval containing the target.
        int lo = 0;
        int hi = n_sub;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (cum[static_cast<size_t>(mid)] <= target) lo = mid;
            else hi = mid;
        }
        const double c_lo = cum[static_cast<size_t>(lo)];
        const double c_hi = cum[static_cast<size_t>(lo) + 1];
        const double span = std::max(c_hi - c_lo, kTiny);
        const double xi_m =
            xi_edge[static_cast<size_t>(lo)] +
            (target - c_lo) *
                (xi_edge[static_cast<size_t>(lo) + 1] -
                 xi_edge[static_cast<size_t>(lo)]) /
                span;
        // Map xi back to u through the baseline faces (linear in xi).
        int c = static_cast<int>(xi_m * static_cast<double>(n_cells));
        c = std::max(0, std::min(n_cells - 1, c));
        const double u_c = baseline_half_faces[static_cast<size_t>(c)];
        const double u_c1 =
            baseline_half_faces[static_cast<size_t>(c) + 1];
        faces[static_cast<size_t>(m)] =
            u_c + (xi_m * static_cast<double>(n_cells) -
                   static_cast<double>(c)) *
                      (u_c1 - u_c);
    }
    return faces;
}

// Compute the deterministic u* centroid from the aggregated histogram:
// mass-weighted |u_parallel| over the low-u_perp conversion bins.
double compute_u_star(const CylindricalVelocityGrid& baseline,
                      const std::vector<double>& agg, int nv, int nmu)
{
    const double uperp_cut = kLowUperpCells * baseline.uperp_widths[0];
    long double mass_sum = 0.0L;
    long double weighted = 0.0L;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = agg[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            if (baseline.uperp_cells[static_cast<size_t>(k)] <= uperp_cut) {
                mass_sum += m;
                weighted +=
                    static_cast<long double>(m) *
                    std::fabs(baseline.upar_cells[static_cast<size_t>(j)]);
            }
        }
    }
    if (!(mass_sum > 0.0L)) {
        // No low-u_perp conversion cell: fall back to the full histogram
        // centroid (the production request histogram only contains
        // conversion cells).
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const double m = agg[static_cast<size_t>(j) * nmu + k];
                if (!(m > 0.0)) continue;
                mass_sum += m;
                weighted +=
                    static_cast<long double>(m) *
                    std::fabs(baseline.upar_cells[static_cast<size_t>(j)]);
            }
        }
    }
    return (mass_sum > 0.0L)
               ? static_cast<double>(weighted / mass_sum)
               : 0.0;
}

// Build a CylindricalVelocityGrid carrying the candidate faces so the
// production geometry, kinetic-energy and subcell-quadrature formulas can
// be reused on the candidate (no production-default init is touched).
CylindricalVelocityGrid make_candidate_grid(
    const TailInterfaceGridCandidate& candidate, int nv, int nmu)
{
    CylindricalVelocityGrid grid;
    grid.upar_faces = candidate.upar_faces;
    // build_cell_geometry_and_moments / build_moment_closure_table assume
    // the geometry arrays are pre-sized exactly as CylindricalVelocityGrid
    // sizes them in init_grid (the production default init is never
    // touched; only the public fields are prepared here).
    grid.upar_cells.resize(static_cast<size_t>(nv));
    grid.upar_widths.resize(static_cast<size_t>(nv));
    grid.upar_center_distances.assign(static_cast<size_t>(nv) + 1, 0.0);
    grid.uperp_faces = candidate.uperp_faces;
    grid.uperp_cells.resize(static_cast<size_t>(nmu));
    grid.uperp_widths.resize(static_cast<size_t>(nmu));
    grid.uperp_ring_areas.resize(static_cast<size_t>(nmu));
    grid.kinetic_energy.resize(static_cast<size_t>(nv) * nmu);
    grid.vx.resize(static_cast<size_t>(nv) * nmu);
    grid.build_cell_geometry_and_moments(nv, nmu);
    grid.build_moment_closure_table(
        Param::temperature_e / (Const::me * Const::c * Const::c));
    return grid;
}

// State moments for the velocity-dt / operator-work model.  They are
// derived once from the input histogram (the same physical state on every
// grid); only the grid samples the coefficients at different points.
struct OfflineStateMoments {
    double density;
    double u_parallel_mean;
    double u_th2;
    double nu0;
    bool valid;
    OfflineStateMoments()
        : density(0.0), u_parallel_mean(0.0), u_th2(0.0), nu0(0.0),
          valid(false)
    {}
};

OfflineStateMoments state_moments_from_histogram(
    const CylindricalVelocityGrid& baseline,
    const std::vector<double>& agg, int nv, int nmu)
{
    OfflineStateMoments state;
    long double total_n = 0.0L;
    long double u_sum = 0.0L;
    long double ke_sum = 0.0L;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = agg[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            total_n += m;
            u_sum += static_cast<long double>(m) *
                     baseline.upar_cells[static_cast<size_t>(j)];
            ke_sum += static_cast<long double>(m) *
                      baseline.kinetic_energy[static_cast<size_t>(j) * nmu +
                                              k];
        }
    }
    if (!(total_n > 0.0L) || !(Param::dx > 0.0)) return state;
    state.density = static_cast<double>(total_n) / Param::dx;
    state.u_parallel_mean =
        static_cast<double>(u_sum / total_n);
    const double ke_per_n = static_cast<double>(ke_sum / total_n);
    const double random_ke =
        ke_per_n - 0.5 * Const::me * Const::c * Const::c *
                       state.u_parallel_mean * state.u_parallel_mean;
    state.u_th2 = (random_ke > 0.0)
                      ? (2.0 / 3.0) * random_ke /
                            (Const::me * Const::c * Const::c)
                      : Param::temperature_e /
                            (Const::me * Const::c * Const::c);
    if (!(state.u_th2 > 0.0)) return state;
    const double u_th = std::sqrt(state.u_th2);
    // nu0 = n e^4 lnL / (4 pi eps0^2 m^2 c^3 u_th^3)  [1/s]
    // (identical to MomentClosureCollisionCoefficients::evaluate).
    state.nu0 =
        state.density * Const::qe * Const::qe * Const::qe * Const::qe *
        kCoulombLog /
        (4.0 * Const::pi * Const::eps0 * Const::eps0 * Const::me * Const::me *
         Const::c * Const::c * Const::c * state.u_th2 * u_th);
    state.valid = std::isfinite(state.nu0) && state.nu0 > 0.0;
    return state;
}

double collision_rate(const OfflineStateMoments& state, double upar,
                      double uperp)
{
    if (!state.valid) return 0.0;
    const double u2 = upar * upar + uperp * uperp;
    return state.nu0 /
           std::pow(1.0 + u2 / std::max(state.u_th2, kTiny), 1.5);
}

// Per-grid velocity-step budget from the production operators:
//   A) T_u: PPM u_parallel force transport.  The departure point of each
//      interior face must stay inside its 4-cell reconstruction stencil
//      (2 cells each side), limited by semi_lagrangian_cfl = 2.5;
//   B) collision cross-diffusion: the explicit half is bounded by
//      velocity_space_cfl = 0.35 * min(du_par^2, du_perp^2) / D_cell with
//      the moment-closure diffusion D = 2*nu(u)*u_th^2 sampled at the cell
//      centres of the evaluated grid.
struct VelocityDtBudget {
    double upar_dt;
    double uperp_dt;
    double total_dt;
    VelocityDtBudget() : upar_dt(0.0), uperp_dt(0.0), total_dt(0.0) {}
};

VelocityDtBudget velocity_dt_budget(
    const CylindricalVelocityGrid& grid, const OfflineStateMoments& state,
    const std::vector<double>& masses)
{
    VelocityDtBudget budget;
    const int nv = static_cast<int>(grid.upar_widths.size());
    const int nmu = static_cast<int>(grid.uperp_widths.size());
    if (nv <= 0 || nmu <= 0) return budget;
    const double a_u_max = Const::qe * kReferenceFieldVm /
                           (Const::me * Const::c);
    if (a_u_max > 0.0) {
        double min_span = std::numeric_limits<double>::infinity();
        for (int f = 1; f < nv; ++f) {
            const double uf = grid.upar_faces[static_cast<size_t>(f)];
            const double left =
                uf - ((f >= 2)
                          ? grid.upar_faces[static_cast<size_t>(f) - 2]
                          : grid.upar_faces[0]);
            const double right =
                ((f + 2 <= nv)
                     ? grid.upar_faces[static_cast<size_t>(f) + 2]
                     : grid.upar_faces[static_cast<size_t>(nv)]) -
                uf;
            min_span = std::min(min_span, std::min(left, right));
        }
        if (std::isfinite(min_span) && min_span > 0.0) {
            budget.upar_dt =
                Param::semi_lagrangian_cfl * min_span / a_u_max;
        }
    }
    if (state.valid) {
        double min_diff = std::numeric_limits<double>::infinity();
        for (int j = 0; j < nv; ++j) {
            const double du2 =
                grid.upar_widths[static_cast<size_t>(j)] *
                grid.upar_widths[static_cast<size_t>(j)];
            for (int k = 0; k < nmu; ++k) {
                const double dp2 =
                    grid.uperp_widths[static_cast<size_t>(k)] *
                    grid.uperp_widths[static_cast<size_t>(k)];
                const double rate = collision_rate(
                    state, grid.upar_cells[static_cast<size_t>(j)],
                    grid.uperp_cells[static_cast<size_t>(k)]);
                const double diffusion = 2.0 * rate * state.u_th2;
                if (!(diffusion > 0.0)) continue;
                const double ratio =
                    std::min(du2, dp2) / diffusion;
                min_diff = std::min(min_diff, ratio);
            }
        }
        if (std::isfinite(min_diff) && min_diff > 0.0) {
            budget.uperp_dt = Param::velocity_space_cfl * min_diff;
        }
    }
    (void)masses;
    budget.total_dt = std::min(
        (budget.upar_dt > 0.0) ? budget.upar_dt
                               : std::numeric_limits<double>::infinity(),
        (budget.uperp_dt > 0.0) ? budget.uperp_dt
                                : std::numeric_limits<double>::infinity());
    if (!std::isfinite(budget.total_dt)) budget.total_dt = 0.0;
    return budget;
}

// Per-grid collision rate maximum over the cells carrying request mass
// (used only for the substep-budget ratio; the same physical state is
// sampled on every grid).
double max_collision_rate(const CylindricalVelocityGrid& grid,
                          const OfflineStateMoments& state,
                          const std::vector<double>& masses)
{
    double max_rate = 0.0;
    const int nv = static_cast<int>(grid.upar_widths.size());
    const int nmu = static_cast<int>(grid.uperp_widths.size());
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = masses[static_cast<size_t>(j) * nmu +
                                    static_cast<size_t>(k)];
            if (!(m > 0.0)) continue;
            const double upar = grid.upar_cells[static_cast<size_t>(j)];
            const double uperp = grid.uperp_cells[static_cast<size_t>(k)];
            const double rate = collision_rate(state, upar, uperp);
            const double a_par = rate * std::fabs(upar - state.u_parallel_mean);
            const double a_perp = rate * std::fabs(uperp);
            const double d = rate * state.u_th2;
            max_rate = std::max(max_rate,
                                a_par + a_perp + 2.0 * d);
        }
    }
    return max_rate;
}

bool parse_profile(const std::string& profile, bool& move_upar,
                   bool& move_uperp)
{
    if (profile == "G0") { move_upar = false; move_uperp = false; return true; }
    if (profile == "Gx") { move_upar = true; move_uperp = false; return true; }
    if (profile == "Gp") { move_upar = false; move_uperp = true; return true; }
    if (profile == "G2") { move_upar = true; move_uperp = true; return true; }
    return false;
}

} // namespace

TailInterfaceReplayResult::TailInterfaceReplayResult()
    : g0_identity_ok(false), number_residual(0.0), max_partition_error(0.0),
      negative_mass_cells(0), max_negative_mass(0.0),
      below_threshold_number_fraction(0.0), center_total(), volume_total(),
      delta_l1(), center_l1(), r_l1(), max_cell_relative(),
      important_max_cell_relative(), min_dupar(0.0), min_duperp(0.0),
      max_adjacent_width_ratio(1.0), estimated_cell_count_ratio(1.0),
      estimated_memory_ratio(1.0), estimated_velocity_dt_ratio(1.0),
      velocity_dt_upar_current(0.0), velocity_dt_upar_candidate(0.0),
      velocity_dt_uperp_current(0.0), velocity_dt_uperp_candidate(0.0),
      velocity_dt_total_current(0.0), velocity_dt_total_candidate(0.0),
      estimated_operator_work_ratio(1.0), scan_cost_ratio(1.0),
      pic_creation_ratio(1.0), collision_pair_ratio(1.0),
      collision_substep_ratio(1.0), center_target_feasible_count(0),
      center_target_failed_count(0), volume_self_feasible_count(0),
      volume_self_failed_count(0), volume_self_sparse_failed_count(0),
      max_sparse_support_count(0), estimated_created_macroparticles(0),
      estimated_particle_ratio_to_center_quartet(0.0), scale_floor_max(),
      candidate_status("RED"), status_reason("uncomputed")
{}

std::vector<double> tail_interface_aggregate_histogram(
    const std::vector<BulkTailVelocityBinAudit>& histogram, int nv, int nmu)
{
    std::vector<double> agg(static_cast<size_t>(nv) * nmu, 0.0);
    for (size_t b = 0; b < histogram.size(); ++b) {
        const BulkTailVelocityBinAudit& bin = histogram[b];
        if (bin.iv < 0 || bin.iv >= nv || bin.imu < 0 || bin.imu >= nmu)
            continue;
        agg[static_cast<size_t>(bin.iv) * nmu + static_cast<size_t>(bin.imu)] +=
            bin.request_number;
    }
    return agg;
}

TailInterfaceGridCandidate build_tail_interface_grid_candidate(
    const CylindricalVelocityGrid& baseline,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    const TailInterfaceGridDesignConfig& config,
    const std::string& profile)
{
    TailInterfaceGridCandidate candidate;
    bool move_upar = false;
    bool move_uperp = false;
    if (!parse_profile(profile, move_upar, move_uperp)) {
        candidate.invalid_reason = "unknown profile " + profile;
        return candidate;
    }
    const int nv = static_cast<int>(baseline.upar_faces.size()) - 1;
    const int nmu = static_cast<int>(baseline.uperp_faces.size()) - 1;
    if (nv <= 0 || nmu <= 0 || nv % 2 != 0) {
        candidate.invalid_reason = "invalid baseline grid";
        return candidate;
    }
    if (!(config.min_width_ratio > 0.0) ||
        !(config.max_adjacent_width_ratio >= 1.0)) {
        candidate.invalid_reason = "invalid width-constraint limits";
        return candidate;
    }
    if (move_upar && !(config.ax >= 0.0)) {
        candidate.invalid_reason = "invalid ax";
        return candidate;
    }
    if (move_uperp && !(config.aperp >= 0.0)) {
        candidate.invalid_reason = "invalid aperp";
        return candidate;
    }
    if (move_upar && !(config.sigma_x_cells > 0.0)) {
        candidate.invalid_reason = "invalid sigma_x_cells";
        return candidate;
    }
    if (move_uperp && !(config.sigma_perp_cells > 0.0)) {
        candidate.invalid_reason = "invalid sigma_perp_cells";
        return candidate;
    }

    const std::vector<double> agg =
        tail_interface_aggregate_histogram(histogram, nv, nmu);
    double total_mass = 0.0;
    for (size_t s = 0; s < agg.size(); ++s) total_mass += agg[s];
    if (!(total_mass > 0.0)) {
        candidate.invalid_reason = "empty request histogram";
        return candidate;
    }

    const int nv_half = nv / 2;
    std::vector<double> upar_half_baseline(static_cast<size_t>(nv_half) + 1,
                                           0.0);
    for (int j = 0; j <= nv_half; ++j)
        upar_half_baseline[static_cast<size_t>(j)] =
            baseline.upar_faces[static_cast<size_t>(nv_half) + j];
    std::vector<double> uperp_baseline(nmu + 1, 0.0);
    for (int k = 0; k <= nmu; ++k)
        uperp_baseline[static_cast<size_t>(k)] =
            baseline.uperp_faces[static_cast<size_t>(k)];

    candidate.grid_name = profile;
    std::ostringstream name_extra;
    if (move_upar) {
        name_extra << "_ax" << config.ax << "_sx" << config.sigma_x_cells;
    }
    if (move_uperp) {
        name_extra << "_ap" << config.aperp << "_sp" << config.sigma_perp_cells;
    }
    candidate.grid_name += name_extra.str();

    if (move_upar || move_uperp) {
        const double u_star =
            compute_u_star(baseline, agg, nv, nmu);
        // Threshold-adjacent u_parallel width: the baseline cell whose
        // centre is closest to u*.
        int j_star = nv_half;
        double best = std::numeric_limits<double>::infinity();
        for (int j = 0; j < nv; ++j) {
            const double d =
                std::fabs(baseline.upar_cells[static_cast<size_t>(j)] - u_star);
            if (d < best) { best = d; j_star = j; }
        }
        const double du_threshold =
            baseline.upar_widths[static_cast<size_t>(j_star)];
        const double du_perp_low = baseline.uperp_widths[0];

        if (move_upar) {
            const double sigma =
                config.sigma_x_cells * du_threshold;
            const std::vector<double> half =
                invert_monitor_half(upar_half_baseline, nv_half, config.ax,
                                    sigma, u_star);
            if (half.size() != static_cast<size_t>(nv_half) + 1 ||
                half.front() != 0.0 ||
                half.back() != upar_half_baseline.back()) {
                candidate.invalid_reason = "u_parallel face inversion failed";
                return candidate;
            }
            candidate.upar_faces.assign(static_cast<size_t>(nv) + 1, 0.0);
            for (int j = 0; j <= nv_half; ++j) {
                candidate.upar_faces[static_cast<size_t>(nv_half) + j] =
                    half[static_cast<size_t>(j)];
            }
            for (int j = 0; j < nv_half; ++j) {
                candidate.upar_faces[static_cast<size_t>(j)] =
                    -candidate.upar_faces[static_cast<size_t>(nv - j)];
            }
        } else {
            candidate.upar_faces = baseline.upar_faces;
        }
        if (move_uperp) {
            const double sigma = config.sigma_perp_cells * du_perp_low;
            const std::vector<double> perp =
                invert_monitor_half(uperp_baseline, nmu, config.aperp, sigma,
                                    0.0);
            if (perp.size() != static_cast<size_t>(nmu) + 1 ||
                perp.front() != 0.0 || perp.back() != uperp_baseline.back()) {
                candidate.invalid_reason = "u_perp face inversion failed";
                return candidate;
            }
            candidate.uperp_faces = perp;
        } else {
            candidate.uperp_faces = baseline.uperp_faces;
        }
    } else {
        candidate.upar_faces = baseline.upar_faces;
        candidate.uperp_faces = baseline.uperp_faces;
    }

    // Face sanity: exact endpoints, mirror symmetry, strict monotonicity.
    if (candidate.upar_faces.front() != baseline.upar_faces.front() ||
        candidate.upar_faces.back() != baseline.upar_faces.back() ||
        candidate.uperp_faces.front() != baseline.uperp_faces.front() ||
        candidate.uperp_faces.back() != baseline.uperp_faces.back()) {
        candidate.invalid_reason = "face endpoints changed";
        return candidate;
    }
    for (int j = 0; j < nv; ++j) {
        if (!(candidate.upar_faces[static_cast<size_t>(j) + 1] >
              candidate.upar_faces[static_cast<size_t>(j)])) {
            candidate.invalid_reason = "u_parallel faces not strictly monotone";
            return candidate;
        }
        const double mirror_error =
            std::fabs(candidate.upar_faces[static_cast<size_t>(j)] +
                      candidate.upar_faces[static_cast<size_t>(nv - j)]);
        const double scale = std::max(
            1.0, std::fabs(candidate.upar_faces[static_cast<size_t>(j)]));
        if (mirror_error > 64.0 * std::numeric_limits<double>::epsilon() *
                               scale) {
            candidate.invalid_reason = "u_parallel mirror symmetry broken";
            return candidate;
        }
    }
    if (candidate.upar_faces[static_cast<size_t>(nv_half)] != 0.0) {
        candidate.invalid_reason = "u_parallel=0 face not fixed";
        return candidate;
    }
    for (int k = 0; k < nmu; ++k) {
        if (!(candidate.uperp_faces[static_cast<size_t>(k) + 1] >
              candidate.uperp_faces[static_cast<size_t>(k)])) {
            candidate.invalid_reason = "u_perp faces not strictly monotone";
            return candidate;
        }
    }

    // Constraint checks: adjacent width ratio and minimum width (relative to
    // the baseline cell containing the new cell centre).  Violations mark
    // the candidate INVALID; nothing is clipped.
    double max_ratio = 1.0;
    for (int j = 1; j < nv; ++j) {
        const double wl =
            candidate.upar_faces[static_cast<size_t>(j)] -
            candidate.upar_faces[static_cast<size_t>(j) - 1];
        const double wr =
            candidate.upar_faces[static_cast<size_t>(j) + 1] -
            candidate.upar_faces[static_cast<size_t>(j)];
        if (wl <= 0.0 || wr <= 0.0) continue;
        max_ratio = std::max(max_ratio, std::max(wl / wr, wr / wl));
    }
    for (int k = 1; k < nmu; ++k) {
        const double wl =
            candidate.uperp_faces[static_cast<size_t>(k)] -
            candidate.uperp_faces[static_cast<size_t>(k) - 1];
        const double wr =
            candidate.uperp_faces[static_cast<size_t>(k) + 1] -
            candidate.uperp_faces[static_cast<size_t>(k)];
        if (wl <= 0.0 || wr <= 0.0) continue;
        max_ratio = std::max(max_ratio, std::max(wl / wr, wr / wl));
    }
    if (max_ratio > config.max_adjacent_width_ratio) {
        candidate.invalid_reason =
            "max adjacent width ratio " + std::to_string(max_ratio) +
            " exceeds limit " +
            std::to_string(config.max_adjacent_width_ratio);
        return candidate;
    }
    for (int j = 0; j < nv; ++j) {
        const double u_center = 0.5 *
            (candidate.upar_faces[static_cast<size_t>(j)] +
             candidate.upar_faces[static_cast<size_t>(j) + 1]);
        int c = 0;
        while (c + 1 < nv &&
               baseline.upar_faces[static_cast<size_t>(c) + 1] <= u_center) {
            ++c;
        }
        const double new_width =
            candidate.upar_faces[static_cast<size_t>(j) + 1] -
            candidate.upar_faces[static_cast<size_t>(j)];
        const double ref_width =
            baseline.upar_widths[static_cast<size_t>(c)];
        if (new_width < config.min_width_ratio * ref_width) {
            candidate.invalid_reason =
                "u_parallel cell " + std::to_string(j) +
                " below minimum width ratio";
            return candidate;
        }
    }
    for (int k = 0; k < nmu; ++k) {
        const double p_center = 0.5 *
            (candidate.uperp_faces[static_cast<size_t>(k)] +
             candidate.uperp_faces[static_cast<size_t>(k) + 1]);
        int c = 0;
        while (c + 1 < nmu &&
               baseline.uperp_faces[static_cast<size_t>(c) + 1] <= p_center) {
            ++c;
        }
        const double new_width =
            candidate.uperp_faces[static_cast<size_t>(k) + 1] -
            candidate.uperp_faces[static_cast<size_t>(k)];
        const double ref_width =
            baseline.uperp_widths[static_cast<size_t>(c)];
        if (new_width < config.min_width_ratio * ref_width) {
            candidate.invalid_reason =
                "u_perp cell " + std::to_string(k) +
                " below minimum width ratio";
            return candidate;
        }
    }
    candidate.valid = true;
    return candidate;
}

std::vector<double> tail_interface_remap_masses(
    const std::vector<double>& from_upar_faces,
    const std::vector<double>& from_uperp_faces,
    const std::vector<double>& to_upar_faces,
    const std::vector<double>& to_uperp_faces,
    const std::vector<double>& from_mass,
    double* max_partition_error_out)
{
    const int nv = static_cast<int>(from_upar_faces.size()) - 1;
    const int nmu = static_cast<int>(from_uperp_faces.size()) - 1;
    std::vector<double> result;
    if (nv <= 0 || nmu <= 0 ||
        from_upar_faces.size() != to_upar_faces.size() ||
        from_uperp_faces.size() != to_uperp_faces.size() ||
        from_mass.size() != static_cast<size_t>(nv) * nmu) {
        if (max_partition_error_out)
            *max_partition_error_out =
                std::numeric_limits<double>::infinity();
        return result;
    }
    result.assign(static_cast<size_t>(nv) * nmu, 0.0);
    double max_partition_error = 0.0;
    for (int j = 0; j < nv; ++j) {
        const double up_lo = from_upar_faces[static_cast<size_t>(j)];
        const double up_hi = from_upar_faces[static_cast<size_t>(j) + 1];
        const double up_width = up_hi - up_lo;
        if (!(up_width > 0.0)) continue;
        // New-cell overlap range for this old cell.
        int jp_lo = static_cast<int>(std::upper_bound(
            to_upar_faces.begin(), to_upar_faces.end(), up_lo) -
                                     to_upar_faces.begin());
        jp_lo = std::max(0, jp_lo - 1);
        int jp_hi = static_cast<int>(std::lower_bound(
            to_upar_faces.begin(), to_upar_faces.end(), up_hi) -
                                     to_upar_faces.begin());
        jp_hi = std::min(nv - 1, jp_hi);
        for (int k = 0; k < nmu; ++k) {
            const double m = from_mass[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            const double pp_lo = from_uperp_faces[static_cast<size_t>(k)];
            const double pp_hi = from_uperp_faces[static_cast<size_t>(k) + 1];
            const double pp2_range =
                pp_hi * pp_hi - pp_lo * pp_lo;
            if (!(pp2_range > 0.0)) continue;
            int kp_lo = static_cast<int>(std::upper_bound(
                to_uperp_faces.begin(), to_uperp_faces.end(), pp_lo) -
                                         to_uperp_faces.begin());
            kp_lo = std::max(0, kp_lo - 1);
            int kp_hi = static_cast<int>(std::lower_bound(
                to_uperp_faces.begin(), to_uperp_faces.end(), pp_hi) -
                                         to_uperp_faces.begin());
            kp_hi = std::min(nmu - 1, kp_hi);
            double theta_sum = 0.0;
            for (int jp = jp_lo; jp <= jp_hi; ++jp) {
                const double n_lo = std::max(
                    up_lo, to_upar_faces[static_cast<size_t>(jp)]);
                const double n_hi = std::min(
                    up_hi, to_upar_faces[static_cast<size_t>(jp) + 1]);
                if (!(n_hi > n_lo)) continue;
                const double du_overlap = n_hi - n_lo;
                const double theta_u = du_overlap / up_width;
                for (int kp = kp_lo; kp <= kp_hi; ++kp) {
                    const double np_lo = std::max(
                        pp_lo, to_uperp_faces[static_cast<size_t>(kp)]);
                    const double np_hi = std::min(
                        pp_hi, to_uperp_faces[static_cast<size_t>(kp) + 1]);
                    if (!(np_hi > np_lo)) continue;
                    const double overlap2 =
                        np_hi * np_hi - np_lo * np_lo;
                    const double theta_p = overlap2 / pp2_range;
                    const double theta = theta_u * theta_p;
                    if (!(theta > 0.0)) continue;
                    result[static_cast<size_t>(jp) * nmu + kp] += m * theta;
                    theta_sum += theta;
                }
            }
            max_partition_error = std::max(
                max_partition_error, std::fabs(theta_sum - 1.0));
        }
    }
    if (max_partition_error_out) *max_partition_error_out = max_partition_error;
    return result;
}

TailInterfaceReplayResult replay_tail_interface_histogram(
    const CylindricalVelocityGrid& baseline,
    const TailInterfaceGridCandidate& candidate,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    double conversion_energy)
{
    TailInterfaceReplayResult result;
    const int nv = static_cast<int>(baseline.upar_faces.size()) - 1;
    const int nmu = static_cast<int>(baseline.uperp_faces.size()) - 1;
    if (nv <= 0 || nmu <= 0 ||
        candidate.upar_faces.size() != baseline.upar_faces.size() ||
        candidate.uperp_faces.size() != baseline.uperp_faces.size() ||
        !(conversion_energy > 0.0)) {
        result.candidate_status = "INVALID";
        result.status_reason = "invalid replay inputs";
        return result;
    }
    const std::vector<double> agg =
        tail_interface_aggregate_histogram(histogram, nv, nmu);
    double input_total = 0.0;
    for (size_t s = 0; s < agg.size(); ++s) input_total += agg[s];
    if (!(input_total > 0.0)) {
        result.candidate_status = "INVALID";
        result.status_reason = "empty request histogram";
        return result;
    }
    const bool is_g0 = candidate.grid_name == "G0";

    // Conservative remap with the production overlap rule.
    double max_partition_error = 0.0;
    const std::vector<double> remapped = tail_interface_remap_masses(
        baseline.upar_faces, baseline.uperp_faces, candidate.upar_faces,
        candidate.uperp_faces, agg, &max_partition_error);
    result.max_partition_error = max_partition_error;
    if (remapped.size() != agg.size()) {
        result.candidate_status = "INVALID";
        result.status_reason = "remap produced no mass";
        return result;
    }
    double remapped_total = 0.0;
    for (size_t s = 0; s < remapped.size(); ++s) remapped_total += remapped[s];
    result.number_residual =
        std::fabs(remapped_total - input_total) /
        std::max(input_total, kTiny);
    for (size_t s = 0; s < remapped.size(); ++s) {
        if (remapped[s] < 0.0) {
            ++result.negative_mass_cells;
            result.max_negative_mass =
                std::min(result.max_negative_mass, remapped[s]);
        }
    }

    // G0 identity: the baseline faces/masses/moments must reconstruct the
    // input exactly; when it fails, every subsequent candidate is INVALID.
    if (is_g0) {
        result.g0_identity_ok =
            face_arrays_identical(candidate.upar_faces, baseline.upar_faces) &&
            face_arrays_identical(candidate.uperp_faces, baseline.uperp_faces) &&
            max_partition_error <= 1.0e-14 && result.negative_mass_cells == 0;
    }

    const CylindricalVelocityGrid cgrid =
        make_candidate_grid(candidate, nv, nmu);
    std::array<long double, BULK_TAIL_MOMENT_COUNT> center_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> volume_sum = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> delta_l1 = {};
    std::array<long double, BULK_TAIL_MOMENT_COUNT> center_l1 = {};
    std::array<double, BULK_TAIL_MOMENT_COUNT> max_cell_rel = {};
    std::array<double, BULK_TAIL_MOMENT_COUNT> important_max_rel = {};
    double below_number = 0.0;
    long long nonzero_cells = 0;
    std::vector<std::vector<double> > scale_floor_per_cell(
        static_cast<size_t>(BULK_TAIL_MOMENT_COUNT));

    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const double m = remapped[static_cast<size_t>(j) * nmu + k];
            if (!(m > 0.0)) continue;
            ++nonzero_cells;
            const double upar = cgrid.upar_cells[static_cast<size_t>(j)];
            const double uperp = cgrid.uperp_cells[static_cast<size_t>(k)];
            const std::array<double, BULK_TAIL_MOMENT_COUNT> center =
                unit_moments(upar, uperp);
            std::vector<TailSubcellNode> nodes =
                TailSubcellQuadrature::nodes(cgrid, j, k);
            std::array<double, BULK_TAIL_MOMENT_COUNT> volume = {};
            std::vector<std::vector<double> > cols;
            std::vector<double> prior;
            double fraction_sum = 0.0;
            for (size_t q = 0; q < nodes.size(); ++q) {
                const std::array<double, BULK_TAIL_MOMENT_COUNT> col =
                    unit_moments(nodes[q].upar, nodes[q].uperp);
                const double w = m * nodes[q].mass_fraction;
                add_moments(volume, col, w);
                fraction_sum += nodes[q].mass_fraction;
                cols.push_back(std::vector<double>(col.begin(), col.end()));
                prior.push_back(w);
            }
            if (std::fabs(fraction_sum - 1.0) > 1.0e-12) {
                result.candidate_status = "INVALID";
                result.status_reason = "subcell quadrature partition broken";
                return result;
            }
            std::array<double, BULK_TAIL_MOMENT_COUNT> center_cell = {};
            add_moments(center_cell, center, m);
            for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
                center_sum[static_cast<size_t>(x)] += center_cell[x];
                volume_sum[static_cast<size_t>(x)] += volume[x];
                const double delta = volume[x] - center_cell[x];
                delta_l1[static_cast<size_t>(x)] += std::fabs(delta);
                center_l1[static_cast<size_t>(x)] += std::fabs(center_cell[x]);
                const double rel = relative_ratio(
                    std::fabs(delta), std::fabs(center_cell[x]));
                max_cell_rel[static_cast<size_t>(x)] =
                    std::max(max_cell_rel[static_cast<size_t>(x)], rel);
                if (m >= kImportantCellFraction * input_total) {
                    important_max_rel[static_cast<size_t>(x)] =
                        std::max(important_max_rel[static_cast<size_t>(x)], rel);
                }
            }
            if (cgrid.kinetic_energy[static_cast<size_t>(j) * nmu + k] <
                conversion_energy) {
                below_number += m;
            }

            // Volume-self feasibility and deterministic sparse reduction.
            const std::vector<double> ref_volume(volume.begin(), volume.end());
            const std::vector<double> ref_center(center_cell.begin(),
                                                 center_cell.end());
            std::vector<double> center_weights;
            if (tail_solve_nonnegative_moment_weights(
                    cols, ref_center, prior, center_weights, 1.0e-10)) {
                ++result.center_target_feasible_count;
            } else {
                ++result.center_target_failed_count;
            }
            std::vector<double> self_weights;
            if (!tail_solve_nonnegative_moment_weights(
                    cols, ref_volume, prior, self_weights, 1.0e-10)) {
                ++result.volume_self_failed_count;
                ++result.volume_self_sparse_failed_count;
                continue;
            }
            ++result.volume_self_feasible_count;
            if (!tail_compress_moment_supports(cols, self_weights, ref_volume,
                                               static_cast<size_t>(kMaxSparseSupport),
                                               1.0e-10)) {
                ++result.volume_self_failed_count;
                ++result.volume_self_sparse_failed_count;
                continue;
            }
            // Scaled residual check with the per-cell scale floors
            // (section 7.11.16C item 3: never divide by a signed moment that
            // can vanish; floor = cell number * velocity/energy scale).
            const double gamma =
                std::sqrt(1.0 + upar * upar + uperp * uperp);
            const double u_abs = std::fabs(upar);
            std::array<double, BULK_TAIL_MOMENT_COUNT> scale_floor = {};
            scale_floor[0] = m;
            scale_floor[1] = m * Const::me * Const::c * u_abs;
            scale_floor[2] = m * Const::qe * Const::c * u_abs / gamma;
            scale_floor[3] = m * Const::me * Const::c * Const::c * (gamma - 1.0);
            scale_floor[4] =
                m * Const::me * Const::c * Const::c * u_abs * u_abs / gamma;
            scale_floor[5] =
                m * Const::me * Const::c * Const::c * uperp * uperp / gamma;
            for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
                scale_floor_per_cell[static_cast<size_t>(x)].push_back(
                    scale_floor[static_cast<size_t>(x)]);
            }
            std::array<long double, BULK_TAIL_MOMENT_COUNT> got = {};
            int support_count = 0;
            for (size_t q = 0; q < self_weights.size(); ++q) {
                if (!(self_weights[q] > 0.0)) continue;
                ++support_count;
                for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
                    got[static_cast<size_t>(x)] +=
                        static_cast<long double>(cols[q][static_cast<size_t>(x)]) *
                        self_weights[q];
                }
            }
            bool residual_ok = true;
            for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
                const double residual =
                    std::fabs(static_cast<double>(got[static_cast<size_t>(x)]) -
                              ref_volume[static_cast<size_t>(x)]);
                const double scaled =
                    residual /
                    std::max(std::fabs(ref_volume[static_cast<size_t>(x)]),
                             scale_floor[static_cast<size_t>(x)]);
                const double limit = (x == BULK_TAIL_MOMENT_N) ? 1.0e-12
                                                               : 1.0e-10;
                if (!(scaled <= limit)) residual_ok = false;
            }
            if (!residual_ok) {
                ++result.volume_self_failed_count;
                ++result.volume_self_sparse_failed_count;
                continue;
            }
            result.max_sparse_support_count =
                std::max(result.max_sparse_support_count, support_count);
            result.estimated_created_macroparticles += support_count;
        }
    }
    result.below_threshold_number_fraction =
        below_number / std::max(remapped_total, kTiny);
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
        result.center_total[static_cast<size_t>(x)] =
            static_cast<double>(center_sum[static_cast<size_t>(x)]);
        result.volume_total[static_cast<size_t>(x)] =
            static_cast<double>(volume_sum[static_cast<size_t>(x)]);
        result.delta_l1[static_cast<size_t>(x)] =
            static_cast<double>(delta_l1[static_cast<size_t>(x)]);
        result.center_l1[static_cast<size_t>(x)] =
            static_cast<double>(center_l1[static_cast<size_t>(x)]);
        result.r_l1[static_cast<size_t>(x)] = relative_ratio(
            result.delta_l1[static_cast<size_t>(x)],
            result.center_l1[static_cast<size_t>(x)]);
        result.max_cell_relative[static_cast<size_t>(x)] =
            max_cell_rel[static_cast<size_t>(x)];
        result.important_max_cell_relative[static_cast<size_t>(x)] =
            important_max_rel[static_cast<size_t>(x)];
        for (size_t s = 0; s < scale_floor_per_cell[static_cast<size_t>(x)].size();
             ++s) {
            result.scale_floor_max[static_cast<size_t>(x)] = std::max(
                result.scale_floor_max[static_cast<size_t>(x)],
                scale_floor_per_cell[static_cast<size_t>(x)][s]);
        }
    }
    result.estimated_particle_ratio_to_center_quartet =
        (nonzero_cells > 0)
            ? static_cast<double>(result.estimated_created_macroparticles) /
                  (kCenterQuartetPerCell * static_cast<double>(nonzero_cells))
            : 0.0;

    // Grid geometry metrics.
    result.min_dupar = *std::min_element(
        cgrid.upar_widths.begin(), cgrid.upar_widths.end());
    result.min_duperp = *std::min_element(
        cgrid.uperp_widths.begin(), cgrid.uperp_widths.end());
    double max_ratio = 1.0;
    for (int j = 1; j < nv; ++j) {
        const double wl = cgrid.upar_widths[static_cast<size_t>(j) - 1];
        const double wr = cgrid.upar_widths[static_cast<size_t>(j)];
        max_ratio = std::max(max_ratio, std::max(wl / wr, wr / wl));
    }
    for (int k = 1; k < nmu; ++k) {
        const double wl = cgrid.uperp_widths[static_cast<size_t>(k) - 1];
        const double wr = cgrid.uperp_widths[static_cast<size_t>(k)];
        max_ratio = std::max(max_ratio, std::max(wl / wr, wr / wl));
    }
    result.max_adjacent_width_ratio = max_ratio;
    result.estimated_cell_count_ratio = 1.0;
    result.estimated_memory_ratio = 1.0;

    // Velocity-step and operator-work budget (same physical state on both
    // grids; only the sampled grid differs).
    const OfflineStateMoments state = state_moments_from_histogram(
        baseline, agg, nv, nmu);
    const VelocityDtBudget current_budget =
        velocity_dt_budget(baseline, state, agg);
    const VelocityDtBudget candidate_budget =
        velocity_dt_budget(cgrid, state, remapped);
    result.velocity_dt_upar_current = current_budget.upar_dt;
    result.velocity_dt_upar_candidate = candidate_budget.upar_dt;
    result.velocity_dt_uperp_current = current_budget.uperp_dt;
    result.velocity_dt_uperp_candidate = candidate_budget.uperp_dt;
    result.velocity_dt_total_current = current_budget.total_dt;
    result.velocity_dt_total_candidate = candidate_budget.total_dt;
    result.estimated_velocity_dt_ratio =
        (current_budget.total_dt > 0.0)
            ? candidate_budget.total_dt / current_budget.total_dt
            : 0.0;

    const double current_rate =
        max_collision_rate(baseline, state, agg);
    const double candidate_rate =
        max_collision_rate(cgrid, state, remapped);
    const int current_substeps = std::max(
        1, static_cast<int>(std::ceil(
               current_rate * current_budget.total_dt /
               kCollisionSubstepTarget)));
    const int candidate_substeps = std::max(
        1, static_cast<int>(std::ceil(
               candidate_rate * current_budget.total_dt /
               kCollisionSubstepTarget)));
    result.collision_substep_ratio =
        static_cast<double>(candidate_substeps) /
        static_cast<double>(current_substeps);
    result.scan_cost_ratio = 1.0;
    result.pic_creation_ratio =
        (nonzero_cells > 0)
            ? static_cast<double>(result.estimated_created_macroparticles) /
                  (kCenterQuartetPerCell * static_cast<double>(nonzero_cells))
            : 1.0;
    // Tail collision pairs are formed per particle (one partner per
    // substep, pair_count = floor(n/2) per cell), so the pair budget scales
    // linearly with the created macro-particle count, not quadratically.
    result.collision_pair_ratio =
        result.pic_creation_ratio;
    // Production-step cost structure (section 7.11.16C): the bulk Vlasov
    // scan dominates the per-step wall time (diagnostics vlasov_s ~0.1 s vs
    // beam_s/collision_s ~1e-4..1e-3 s), so the fixed cell scan carries 0.7
    // of the weight and the PIC creation / collision budgets 0.15 each.
    // The contributions themselves are written to the result so the model
    // is transparent; only the final ratio is gated (<= 1.5).
    result.estimated_operator_work_ratio =
        0.7 * result.scan_cost_ratio +
        0.15 * result.pic_creation_ratio +
        0.15 * result.collision_pair_ratio * result.collision_substep_ratio;

    // Classification (section 7.11.16C offline gates).
    if (!candidate.valid) {
        result.candidate_status = "INVALID";
        result.status_reason = candidate.invalid_reason;
        return result;
    }
    if (is_g0 && !result.g0_identity_ok) {
        result.candidate_status = "INVALID";
        result.status_reason = "G0 identity check failed";
        return result;
    }
    const bool hard_ok =
        result.number_residual <= 1.0e-12 &&
        result.negative_mass_cells == 0 &&
        result.below_threshold_number_fraction <= 1.0e-6 &&
        result.volume_self_failed_count == 0 &&
        result.volume_self_sparse_failed_count == 0 &&
        result.max_sparse_support_count <= kMaxSparseSupport &&
        result.estimated_particle_ratio_to_center_quartet <= 2.0 &&
        result.estimated_cell_count_ratio == 1.0 &&
        result.estimated_memory_ratio == 1.0 &&
        result.estimated_velocity_dt_ratio >= 0.8 &&
        result.estimated_operator_work_ratio <= 1.5;
    if (!hard_ok) {
        result.candidate_status = "RED";
        result.status_reason = "hard offline gate failed";
        return result;
    }
    const bool moment_ok =
        result.r_l1[BULK_TAIL_MOMENT_JX] <= 1.0e-3 &&
        result.r_l1[BULK_TAIL_MOMENT_K] <= 1.0e-3 &&
        result.r_l1[BULK_TAIL_MOMENT_PIXX] <= 1.0e-3 &&
        result.r_l1[BULK_TAIL_MOMENT_PIPERP] <= 1.0e-3;
    bool important_ok = true;
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x)
        important_ok = important_ok &&
                       result.important_max_cell_relative[static_cast<size_t>(x)] <=
                           1.0e-2;
    if (moment_ok && important_ok) {
        result.candidate_status = "GREEN";
        result.status_reason = "all 16B offline gates pass";
        return result;
    }
    // Gray band: hard gates pass but the physics-moment gates stay in
    // (1e-3, 1e-2] (R_L1) or (1e-2, 1e-1] (important-cell relative diff).
    bool gray_ok = true;
    for (int x = 0; x < BULK_TAIL_MOMENT_COUNT; ++x) {
        gray_ok = gray_ok &&
                  result.r_l1[static_cast<size_t>(x)] <= 1.0e-2 &&
                  result.important_max_cell_relative[static_cast<size_t>(x)] <=
                      1.0e-1;
    }
    result.candidate_status = gray_ok ? "GRAY" : "RED";
    result.status_reason = gray_ok ? "physics-moment gates in gray band"
                                   : "physics-moment gate violation beyond gray band";
    return result;
}
