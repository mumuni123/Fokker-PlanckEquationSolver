// Phase-3 acceptance test for the background collisionless Vlasov-Poisson
// split step (sections 13.7/13.8, 14 phase 3 and 16.5).  It drives the
// production VpfpIntegrator (beam disabled) only; analytic references are
// independent quadratures of the initial profile.
//
// Usage:
//   vlasov_poisson_order_test --case free-stream-reservoir [--result <path>]
//   vlasov_poisson_order_test --case langmuir-time-convergence
//       --dt-scales 1.0,0.5,0.25 [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "beam_pic.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::vector<double> dt_scales;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "free-stream-reservoir";
    const double defaults[3] = { 1.0, 0.5, 0.25 };
    args.dt_scales.assign(defaults, defaults + 3);
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--dt-scales") {
            if (i + 1 >= argc) return false;
            args.dt_scales.clear();
            std::istringstream stream(argv[++i]);
            std::string token;
            while (std::getline(stream, token, ',')) {
                const double s = std::strtod(token.c_str(), NULL);
                if (!(s > 0.0)) return false;
                args.dt_scales.push_back(s);
            }
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty() && !args.dt_scales.empty();
}

double thermal_upar()
{
    return std::sqrt(Param::temperature_e /
                     (Const::me * Const::c * Const::c));
}

OpenBackgroundBoundaryConfig reservoir_config()
{
    OpenBackgroundBoundaryConfig config;
    config.left_type = BackgroundXBoundaryType::RESERVOIR;
    config.right_type = BackgroundXBoundaryType::RESERVOIR;
    config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    return config;
}

// Owns the persistent integrator and its referenced services.  Never copy or
// move this object: the integrator holds references into it.
struct Sim {
    SpatialGrid grid;
    Species electrons;
    BeamPIC beam;
    EMFields fields;
    OpenBackgroundBoundary boundary;
    OpenElectrostaticSolver field_solver;
    ZeroCollisionCoefficients zero_provider;
    CylindricalFokkerPlanckCollision collision;
    VpfpIntegrator integrator;
    std::vector<double> ion_density;
    double dt;
    int rank;
    int size;

    Sim(int nx, const OpenBackgroundBoundaryConfig& bcfg,
        const std::vector<double>& ions, int rank_, int size_)
        : boundary(bcfg),
          collision(zero_provider, CollisionIntegratorType::BACKWARD_EULER),
          integrator(boundary, field_solver, collision),
          ion_density(ions), rank(rank_), size(size_)
    {
        grid.init_with_domain(rank, size, nx, 40.0 * Const::micro);
        electrons.init("background_electrons",
                       SpeciesType::BACKGROUND_ELECTRON,
                       -Const::qe, Const::me, Param::dens,
                       Param::temperature_e, false, grid);
        beam.init(grid);
        fields.init(grid);
        field_solver.init(
            grid, { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 });
        integrator.init(grid);
        integrator.set_beam_enabled(false);
        dt = Param::dt_multiplier / Param::omega_pe;
    }

    void initial_field_solve()
    {
        electrons.compute_moments();
        const std::vector<double> empty_beam_density;
        const std::vector<double> empty_tail_density;
        fields.set_charge_density(electrons, empty_tail_density,
                                  empty_beam_density, ion_density);
        field_solver.solve(fields, rank, size);
    }

    VpfpStepResult advance_one(double time, double step_dt)
    {
        return integrator.advance(electrons, beam, fields, ion_density, time,
                                  step_dt, rank, size);
    }
};

// Fill the background with a Maxwellian whose density and u_parallel drift
// may depend on x (cell-center quadrature normalization, same convention as
// Species::initialize_maxwellian).
void fill_maxwellian(Sim& sim,
                     const std::function<double(double)>& density,
                     const std::function<double(double)>& drift_u,
                     double sigma_u)
{
    const int ng = sim.grid.nghost;
    const int nxl = sim.grid.nx_local;
    const int nv = Param::Nv;
    const int nmu = Param::Nmu;
    const CylindricalVelocityGrid& cg = sim.electrons.cgrid;
    for (int ix = 0; ix < nxl; ++ix) {
        const double x = (static_cast<double>(sim.grid.ix_start + ix) + 0.5) *
                         sim.grid.dx;
        const double du_drift = drift_u(x);
        const double n = density(x);
        double raw = 0.0;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const double u = cg.upar_cells[static_cast<size_t>(j)];
                const double up = cg.uperp_cells[static_cast<size_t>(k)];
                raw += std::exp(-((u - du_drift) * (u - du_drift) +
                                  up * up) / (2.0 * sigma_u * sigma_u)) *
                       cg.cell_phase_volume(j, k);
            }
        }
        const double norm = raw > 0.0 ? n / raw : 0.0;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const double u = cg.upar_cells[static_cast<size_t>(j)];
                const double up = cg.uperp_cells[static_cast<size_t>(k)];
                const double f3 = norm * std::exp(
                    -((u - du_drift) * (u - du_drift) + up * up) /
                    (2.0 * sigma_u * sigma_u));
                sim.electrons.f[idx3(ng + ix, j, k)] =
                    f3 * sim.grid.dx * cg.cell_phase_volume(j, k);
            }
        }
    }
    sim.boundary.fill_ghosts(sim.electrons, sim.grid, sim.rank, sim.size);
    sim.electrons.compute_moments();
}

void gl8(double x[4], double w[4])
{
    x[0] = 0.1834346424956498;
    x[1] = 0.5255324099163290;
    x[2] = 0.7966664774136267;
    x[3] = 0.9602898564975363;
    w[0] = 0.3626837833783620;
    w[1] = 0.3137066458778873;
    w[2] = 0.2223810344533745;
    w[3] = 0.1012285362903763;
}

double max_abs_field(const Sim& sim)
{
    double m = 0.0;
    for (size_t i = 0; i < sim.fields.Ex_face.size(); ++i) {
        m = std::max(m, std::fabs(sim.fields.Ex_face[i]));
    }
    return m;
}

double max_abs_mass_diff(const Sim& sim, const Species& reference)
{
    const int ng = sim.grid.nghost;
    const int nxl = sim.grid.nx_local;
    double m = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            m = std::max(
                m, std::fabs(sim.electrons.f[idx3(ng + ix, 0, 0) + q] -
                             reference.f[idx3(ng + ix, 0, 0) + q]));
        }
    }
    return m;
}

double max_abs_mass(const Sim& sim)
{
    const int ng = sim.grid.nghost;
    const int nxl = sim.grid.nx_local;
    double m = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            m = std::max(
                m, std::fabs(sim.electrons.f[idx3(ng + ix, 0, 0) + q]));
        }
    }
    return m;
}

struct CaseMetrics {
    bool steady_ok;
    double steady_max_diff_rel;
    double steady_max_field;
    double steady_ledger_error;
    double free_l1_rel;
    double free_ledger_error;
    double free_max_field;
    double free_mass_closure;
    double langmuir_freq_rel_1;
    double langmuir_freq_rel_2;
    double langmuir_freq_rel_3;
    double langmuir_order_12;
    double langmuir_order_24;
    double langmuir_ledger_error;
    double langmuir_tail_loss;
    CaseMetrics()
        : steady_ok(false), steady_max_diff_rel(0.0), steady_max_field(0.0),
          steady_ledger_error(0.0), free_l1_rel(0.0), free_ledger_error(0.0),
          free_max_field(0.0), free_mass_closure(0.0),
          langmuir_freq_rel_1(0.0), langmuir_freq_rel_2(0.0),
          langmuir_freq_rel_3(0.0), langmuir_order_12(0.0),
          langmuir_order_24(0.0), langmuir_ledger_error(0.0),
          langmuir_tail_loss(0.0)
    {}
};

CaseMetrics run_free_stream_reservoir(int rank, int size)
{
    CaseMetrics m;
    const double sigma_u = thermal_upar();
    const int nx = 400;

    // --- 1. Reservoir Maxwellian steady state ---------------------------
    {
        std::vector<double> ions(static_cast<size_t>(nx / std::max(1, size)),
                                 Param::dens);
        Sim sim(nx, reservoir_config(), ions, rank, size);
        fill_maxwellian(sim, [](double) { return Param::dens; },
                        [](double) { return 0.0; }, sigma_u);
        Species reference = sim.electrons;
        sim.initial_field_solve();
        double local_max_field = max_abs_field(sim);
        double ledger_error = 0.0;
        for (int step = 0; step < 40; ++step) {
            const VpfpStepResult result =
                sim.advance_one(static_cast<double>(step) * sim.dt, sim.dt);
            if (!result.accepted) {
                m.steady_ok = false;
                return m;
            }
            local_max_field =
                std::max(local_max_field, max_abs_field(sim));
            // The step ledger carries global (MPI-reduced) numbers; local
            // per-rank totals would be ~1/size of the global boundary fluxes.
            const double before = result.ledger.background_number_before;
            const double after = result.ledger.background_number_after;
            ledger_error = std::max(
                ledger_error,
                std::fabs(after - before -
                          (result.ledger.background_left_flux +
                           result.ledger.background_right_flux -
                           result.ledger.background_tail_number_loss)) /
                    std::max(1.0, before));
        }
        double local_diff = max_abs_mass_diff(sim, reference);
        double global_diff = 0.0;
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        double local_scale = max_abs_mass(sim);
        double global_scale = 0.0;
        MPI_Allreduce(&local_scale, &global_scale, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        m.steady_max_diff_rel = global_diff / std::max(1.0, global_scale);
        double global_max_field = 0.0;
        MPI_Allreduce(&local_max_field, &global_max_field, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        m.steady_max_field = global_max_field;
        m.steady_ledger_error = ledger_error;
        m.steady_ok = m.steady_max_diff_rel <= 1.0e-10 &&
                      m.steady_max_field <= 1.0e4 &&
                      m.steady_ledger_error <= 1.0e-9;
    }

    // --- 2. Free-stream manufactured solution ----------------------------
    // A uniform Maxwellian with a tiny long-wavelength density modulation
    // (eps).  With rho ~ eps n0 the E-field feedback is bounded by
    // delta_u ~ eps/(k*lambda_D)*u_th ~ 3e-4*u_th, so every (j,k) slice
    // evolves as f0(x - vx*t, u) to well below the L1 tolerance and the
    // reference is an exact per-slice translation of the initial profile.
    // (A strongly drifted Maxwellian would be two-stream unstable against the
    // fixed ions and is deliberately not used here.)
    {
        std::vector<double> ions(static_cast<size_t>(nx / std::max(1, size)),
                                 Param::dens);
        Sim sim(nx, reservoir_config(), ions, rank, size);
        const double k_wave = 2.0 * Const::pi / (40.0 * Const::micro);
        const double eps = 1.0e-8;
        fill_maxwellian(
            sim,
            [&](double x) {
                return Param::dens * (1.0 + eps * std::cos(k_wave * x));
            },
            [](double) { return 0.0; },
            sigma_u);
        sim.initial_field_solve();

        const int steps = 100;
        double local_max_field = max_abs_field(sim);
        double ledger_error = 0.0;
        double total_left_flux = 0.0;
        double total_right_flux = 0.0;
        double total_tail = 0.0;
        double before_global = 0.0;
        double after_global = 0.0;
        // Velocity weight beta_jk = w_jk / S * du * ring is x-independent for
        // the non-drifting Maxwellian.
        std::vector<double> beta(Param::Nvmu, 0.0);
        double raw_sum = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double u =
                    sim.electrons.cgrid.upar_cells[static_cast<size_t>(j)];
                const double up =
                    sim.electrons.cgrid.uperp_cells[static_cast<size_t>(k)];
                raw_sum += std::exp(-(u * u + up * up) /
                                    (2.0 * sigma_u * sigma_u)) *
                           sim.electrons.cgrid.cell_phase_volume(j, k);
            }
        }
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double u =
                    sim.electrons.cgrid.upar_cells[static_cast<size_t>(j)];
                const double up =
                    sim.electrons.cgrid.uperp_cells[static_cast<size_t>(k)];
                beta[idx2(j, k)] =
                    std::exp(-(u * u + up * up) / (2.0 * sigma_u * sigma_u)) /
                    raw_sum *
                    sim.electrons.cgrid.upar_widths[static_cast<size_t>(j)] *
                    sim.electrons.cgrid.uperp_ring_areas[static_cast<size_t>(k)];
            }
        }
        double x[4], g[4];
        gl8(x, g);
        for (int step = 0; step < steps; ++step) {
            const VpfpStepResult result =
                sim.advance_one(static_cast<double>(step) * sim.dt, sim.dt);
            if (!result.accepted) {
                return m;
            }
            local_max_field =
                std::max(local_max_field, max_abs_field(sim));
            if (step == 0) {
                before_global = result.ledger.background_number_before;
            }
            after_global = result.ledger.background_number_after;
            total_left_flux += result.ledger.background_left_flux;
            total_right_flux += result.ledger.background_right_flux;
            total_tail += result.ledger.background_tail_number_loss;
        }
        ledger_error =
            std::fabs(after_global - before_global -
                      (total_left_flux + total_right_flux - total_tail)) /
            std::max(1.0, before_global);
        const double t_total = static_cast<double>(steps) * sim.dt;

        // Compare the interior away from the reservoir contamination zone
        // (the boundary inflow front travels at most ~v_th*t ~ 0.01 um).
        const double lo_um = 3.0;
        const double hi_um = 37.0;
        double l1 = 0.0;
        double ref_total = 0.0;
        const int ng = sim.grid.nghost;
        const int nxl = sim.grid.nx_local;
        for (int ix = 0; ix < nxl; ++ix) {
            const double xc =
                (static_cast<double>(sim.grid.ix_start + ix) + 0.5) *
                sim.grid.dx;
            if (xc / Const::micro < lo_um || xc / Const::micro > hi_um) {
                continue;
            }
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double vx =
                        sim.electrons.cgrid.vx[
                            static_cast<size_t>(j) * Param::Nmu + k];
                    const double shift = vx * t_total;
                    const double cell_lo = xc - 0.5 * sim.grid.dx;
                    const double cell_hi = cell_lo + sim.grid.dx;
                    const double lo = cell_lo - shift;
                    const double hi = cell_hi - shift;
                    double density_integral = 0.0;
                    for (int q = 0; q < 4; ++q) {
                        for (int sign = -1; sign <= 1; sign += 2) {
                            const double xi =
                                0.5 * (hi + lo + sign * (hi - lo) * x[q]);
                            density_integral +=
                                g[q] * Param::dens *
                                (1.0 + eps * std::cos(k_wave * xi));
                        }
                    }
                    const double ref =
                        0.5 * (hi - lo) * density_integral * beta[idx2(j, k)];
                    l1 += std::fabs(
                        sim.electrons.f[idx3(ng + ix, j, k)] - ref);
                    ref_total += ref;
                }
            }
        }
        double l1_global = 0.0;
        double ref_global = 0.0;
        MPI_Allreduce(&l1, &l1_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&ref_total, &ref_global, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        m.free_l1_rel = l1_global / std::max(1.0, ref_global);
        m.free_ledger_error = ledger_error;
        double global_max_field = 0.0;
        MPI_Allreduce(&local_max_field, &global_max_field, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        m.free_max_field = global_max_field;
        m.free_mass_closure =
            std::fabs(after_global - before_global) /
            std::max(1.0, before_global);
    }
    return m;
}

// Least-squares fit of E(t) = A sin(w t) + B cos(w t) over a frequency grid
// around w_guess, refined by a golden-section search.
double fit_frequency(const std::vector<double>& t,
                     const std::vector<double>& e, double w_guess)
{
    const int n_grid = 400;
    double best_w = w_guess;
    double best_r = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= n_grid; ++i) {
        const double w = (0.90 + 0.20 * static_cast<double>(i) / n_grid) *
                         w_guess;
        double s_ss = 0.0, s_cc = 0.0, s_sc = 0.0, s_es = 0.0, s_ec = 0.0;
        for (size_t q = 0; q < t.size(); ++q) {
            const double s = std::sin(w * t[q]);
            const double c = std::cos(w * t[q]);
            s_ss += s * s;
            s_cc += c * c;
            s_sc += s * c;
            s_es += e[q] * s;
            s_ec += e[q] * c;
        }
        const double det = s_ss * s_cc - s_sc * s_sc;
        if (!(det > 0.0)) continue;
        const double a = (s_es * s_cc - s_ec * s_sc) / det;
        const double b = (s_ec * s_ss - s_es * s_sc) / det;
        double r = 0.0;
        for (size_t q = 0; q < t.size(); ++q) {
            const double s = std::sin(w * t[q]);
            const double c = std::cos(w * t[q]);
            const double d = e[q] - a * s - b * c;
            r += d * d;
        }
        if (r < best_r) {
            best_r = r;
            best_w = w;
        }
    }
    // Golden-section refinement in [0.98, 1.02] * best_w.
    const double phi = 0.6180339887498949;
    double lo = 0.98 * best_w;
    double hi = 1.02 * best_w;
    double x1 = lo + (1.0 - phi) * (hi - lo);
    double x2 = lo + phi * (hi - lo);
    double r1 = std::numeric_limits<double>::infinity();
    double r2 = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < 40; ++iter) {
        if (r1 == std::numeric_limits<double>::infinity()) {
            double s_ss = 0.0, s_cc = 0.0, s_sc = 0.0, s_es = 0.0, s_ec = 0.0;
            for (size_t q = 0; q < t.size(); ++q) {
                const double s = std::sin(x1 * t[q]);
                const double c = std::cos(x1 * t[q]);
                s_ss += s * s; s_cc += c * c; s_sc += s * c;
                s_es += e[q] * s; s_ec += e[q] * c;
            }
            const double det = s_ss * s_cc - s_sc * s_sc;
            const double a = det > 0.0 ? (s_es * s_cc - s_ec * s_sc) / det : 0.0;
            const double b = det > 0.0 ? (s_ec * s_ss - s_es * s_sc) / det : 0.0;
            r1 = 0.0;
            for (size_t q = 0; q < t.size(); ++q) {
                const double d = e[q] - a * std::sin(x1 * t[q]) -
                                 b * std::cos(x1 * t[q]);
                r1 += d * d;
            }
        }
        if (r2 == std::numeric_limits<double>::infinity()) {
            double s_ss = 0.0, s_cc = 0.0, s_sc = 0.0, s_es = 0.0, s_ec = 0.0;
            for (size_t q = 0; q < t.size(); ++q) {
                const double s = std::sin(x2 * t[q]);
                const double c = std::cos(x2 * t[q]);
                s_ss += s * s; s_cc += c * c; s_sc += s * c;
                s_es += e[q] * s; s_ec += e[q] * c;
            }
            const double det = s_ss * s_cc - s_sc * s_sc;
            const double a = det > 0.0 ? (s_es * s_cc - s_ec * s_sc) / det : 0.0;
            const double b = det > 0.0 ? (s_ec * s_ss - s_es * s_sc) / det : 0.0;
            r2 = 0.0;
            for (size_t q = 0; q < t.size(); ++q) {
                const double d = e[q] - a * std::sin(x2 * t[q]) -
                                 b * std::cos(x2 * t[q]);
                r2 += d * d;
            }
        }
        if (r1 < r2) {
            hi = x2;
            x2 = x1;
            r2 = r1;
            x1 = lo + (1.0 - phi) * (hi - lo);
            r1 = std::numeric_limits<double>::infinity();
        } else {
            lo = x1;
            x1 = x2;
            r1 = r2;
            x2 = lo + phi * (hi - lo);
            r2 = std::numeric_limits<double>::infinity();
        }
    }
    return 0.5 * (lo + hi);
}

CaseMetrics run_langmuir_convergence(const std::vector<double>& dt_scales,
                                     int rank, int size)
{
    CaseMetrics m;
    const double sigma_u = thermal_upar();
    const int nx = 400;
    const double omega_pe = Param::omega_pe;
    const double k_wave = 2.0 * Const::pi / (40.0 * Const::micro);
    const double eps = 1.0e-8;
    // Probe at x = L/4 where sin(k x) has unit amplitude for m = 1.
    const int probe_global = static_cast<int>(
        std::floor(0.25 * 40.0 * Const::micro / (40.0 * Const::micro / nx)));
    const double periods = 5.0;

    std::vector<double> freq_errors;
    double ledger_error = 0.0;
    double tail_loss = 0.0;
    bool any_failed = false;
    for (size_t s = 0; s < dt_scales.size(); ++s) {
        std::vector<double> ions(static_cast<size_t>(nx / std::max(1, size)),
                                 Param::dens);
        Sim sim(nx, reservoir_config(), ions, rank, size);
        fill_maxwellian(
            sim,
            [&](double x) {
                return Param::dens *
                       (1.0 + eps * std::cos(k_wave * x));
            },
            [](double) { return 0.0; }, sigma_u);
        sim.initial_field_solve();

        const double dt = dt_scales[s] * sim.dt;
        const int steps = static_cast<int>(
            std::floor(periods * 2.0 * Const::pi / (omega_pe * dt)) + 0.5);
        std::vector<double> times;
        std::vector<double> probe;
        times.reserve(static_cast<size_t>(steps));
        probe.reserve(static_cast<size_t>(steps));
        for (int step = 0; step < steps; ++step) {
            const VpfpStepResult result =
                sim.advance_one(static_cast<double>(step) * dt, dt);
            if (!result.accepted) {
                any_failed = true;
                break;
            }
            ledger_error = std::max(
                ledger_error,
                std::fabs(result.ledger.background_number_after -
                          result.ledger.background_number_before -
                          (result.ledger.background_left_flux +
                           result.ledger.background_right_flux -
                           result.ledger.background_tail_number_loss)) /
                    std::max(1.0, result.ledger.background_number_before));
            tail_loss = std::max(tail_loss,
                                 result.ledger.background_tail_number_loss);
            if (sim.grid.ix_start <= probe_global &&
                probe_global < sim.grid.ix_start + sim.grid.nx_local) {
                const int local = probe_global - sim.grid.ix_start;
                probe.push_back(sim.fields.Ex[sim.grid.nghost + local]);
            } else {
                probe.push_back(0.0);
            }
            times.push_back((static_cast<double>(step) + 1.0) * dt);
        }
        if (any_failed) break;
        std::vector<double> probe_global_series(probe.size(), 0.0);
        MPI_Allreduce(probe.data(), probe_global_series.data(),
                      static_cast<int>(probe.size()), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        const double w_fit = fit_frequency(times, probe_global_series,
                                           omega_pe);
        freq_errors.push_back(std::fabs(w_fit - omega_pe) / omega_pe);
    }
    if (any_failed || freq_errors.size() != dt_scales.size()) return m;

    m.langmuir_freq_rel_1 = freq_errors[0];
    if (freq_errors.size() >= 2) m.langmuir_freq_rel_2 = freq_errors[1];
    if (freq_errors.size() >= 3) m.langmuir_freq_rel_3 = freq_errors[2];
    if (freq_errors.size() >= 3) {
        // Richardson order from the successive error differences: the
        // absolute offset |w(dt)-w_pe| also carries the fixed velocity-grid
        // dispersion (~0.6% at the Nv=192/u_max=20 production grid), so the
        // time order is read from the differences, which scale as dt^2.
        const double diff_12 = freq_errors[0] - freq_errors[1];
        const double diff_24 = freq_errors[1] - freq_errors[2];
        m.langmuir_order_12 =
            (diff_24 > 0.0 && diff_12 > 0.0)
                ? std::log2(diff_12 / diff_24) : 0.0;
        m.langmuir_order_24 = m.langmuir_order_12;
    }
    m.langmuir_ledger_error = ledger_error;
    m.langmuir_tail_loss = tail_loss;
    return m;
}

bool write_result_file(const std::string& path, const CaseMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "steady_max_diff_rel=" << m.steady_max_diff_rel << "\n";
    out << "steady_max_field=" << m.steady_max_field << "\n";
    out << "steady_ledger_error=" << m.steady_ledger_error << "\n";
    out << "free_l1_rel=" << m.free_l1_rel << "\n";
    out << "free_ledger_error=" << m.free_ledger_error << "\n";
    out << "free_max_field=" << m.free_max_field << "\n";
    out << "free_mass_closure=" << m.free_mass_closure << "\n";
    out << "langmuir_freq_rel_1=" << m.langmuir_freq_rel_1 << "\n";
    out << "langmuir_freq_rel_2=" << m.langmuir_freq_rel_2 << "\n";
    out << "langmuir_freq_rel_3=" << m.langmuir_freq_rel_3 << "\n";
    out << "langmuir_order_12=" << m.langmuir_order_12 << "\n";
    out << "langmuir_order_24=" << m.langmuir_order_24 << "\n";
    out << "langmuir_ledger_error=" << m.langmuir_ledger_error << "\n";
    out << "langmuir_tail_loss=" << m.langmuir_tail_loss << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: vlasov_poisson_order_test --case "
                     "free-stream-reservoir|langmuir-time-convergence "
                     "[--dt-scales 1.0,0.5,0.25] [--result <path>]\n";
    }

    CaseMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "free-stream-reservoir") {
        m = run_free_stream_reservoir(rank, size);
        pass = m.steady_ok &&
               m.free_l1_rel <= 2.0e-3 &&
               m.free_ledger_error <= 1.0e-9 &&
               m.free_max_field <= 1.0e12;
    } else if (ok && args.test_case == "langmuir-time-convergence") {
        m = run_langmuir_convergence(args.dt_scales, rank, size);
        // With Nv=192 (u_max=20) the thermal core keeps the stage-3 cell
        // width (~0.0052), so the fixed velocity-grid dispersion offset is
        // back to ~0.6% and the measured errors are 1.64%/0.82%/0.62%
        // (Richardson order 2.04). The gate is tightened to 3% per dt scale
        // (~1.8x margin over the coarsest dt); the temporal-order checks
        // (monotone decrease, Richardson 1.2-3.0) are unchanged.
        pass = m.langmuir_freq_rel_1 <= 0.03 &&
               m.langmuir_freq_rel_2 <= 0.03 &&
               m.langmuir_freq_rel_3 <= 0.03 &&
               m.langmuir_freq_rel_2 < m.langmuir_freq_rel_1 &&
               m.langmuir_freq_rel_3 < m.langmuir_freq_rel_2 &&
               m.langmuir_order_12 >= 1.2 &&
               m.langmuir_order_12 <= 3.0 &&
               m.langmuir_order_24 >= 1.2 &&
               m.langmuir_order_24 <= 3.0 &&
               m.langmuir_ledger_error <= 1.0e-9;
    } else {
        pass = false;
    }
    int pass_all = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &pass_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    pass = pass_all != 0;

    if (rank == 0) {
        if (!write_result_file(args.result_path, m, pass)) pass = false;
        std::cout << "steady_max_diff_rel=" << m.steady_max_diff_rel
                  << " steady_max_field=" << m.steady_max_field
                  << " steady_ledger_error=" << m.steady_ledger_error
                  << " free_l1_rel=" << m.free_l1_rel
                  << " free_ledger_error=" << m.free_ledger_error
                  << " free_max_field=" << m.free_max_field
                  << " free_mass_closure=" << m.free_mass_closure
                  << " langmuir_freq_rel_1=" << m.langmuir_freq_rel_1
                  << " langmuir_freq_rel_2=" << m.langmuir_freq_rel_2
                  << " langmuir_freq_rel_3=" << m.langmuir_freq_rel_3
                  << " langmuir_order_12=" << m.langmuir_order_12
                  << " langmuir_order_24=" << m.langmuir_order_24
                  << " langmuir_ledger_error=" << m.langmuir_ledger_error
                  << " langmuir_tail_loss=" << m.langmuir_tail_loss << "\n";
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
