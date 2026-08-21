// Stage H7 acceptance for the bulk cylindrical Fokker-Planck collision
// module (sections 10.2, 14.11 and 19.3):
//   * number-conservation        : prescribed constant A/D conserves mass;
//   * drift-diffusion-moments    : first/second velocity moments follow the
//                                  input coefficients (A dt, 2 D dt);
//   * maxwellian-equilibrium     : a Maxwellian at the closure temperature
//                                  is stationary to discretization error;
//   * h-theorem                  : a non-Maxwellian relaxes with increasing
//                                  entropy;
//   * cross-diffusion-supported  : d_parallel_perp is accepted (no longer
//                                  rejected) and number is conserved;
//   * reservoir-accounting       : the tracked collision-reservoir energy
//                                  exactly equals the bulk kinetic-energy
//                                  change.
//
// Usage:
//   cylindrical_fp_collision_test [--case <case|all>] [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "all";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty();
}

struct Metrics {
    bool pass;
    double mass_change;
    double mass_scale;
    double drift_rel_error;
    double var_rel_error;
    double equilibrium_l2_rel;
    double entropy_change;
    double reservoir_rel_error;
    bool cross_accepted;
    Metrics()
        : pass(false), mass_change(0.0), mass_scale(1.0),
          drift_rel_error(0.0), var_rel_error(0.0),
          equilibrium_l2_rel(0.0), entropy_change(0.0),
          reservoir_rel_error(0.0), cross_accepted(false)
    {}
};

Species make_species(const SpatialGrid& grid)
{
    Species species;
    species.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                 -Const::qe, Const::me, Param::dens, Param::temperature_e,
                 false, grid);
    return species;
}

// Plain-sum velocity moments (the cylindrical mass representation folds the
// 2 pi u_perp weight into f).
void velocity_moments(const Species& species, const SpatialGrid& grid,
                      double& mean_u, double& variance_u,
                      double& mass_total, double& energy_total)
{
    const int ng = grid.nghost;
    mean_u = 0.0;
    variance_u = 0.0;
    mass_total = 0.0;
    energy_total = 0.0;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = ng + ix;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double m = species.f[idx3(sx, iv, imu)];
                const double u = species.cgrid.upar_cells[iv];
                mass_total += m;
                mean_u += m * u;
                energy_total +=
                    m * species.cgrid.kinetic_energy[idx2(iv, imu)];
            }
        }
    }
    if (mass_total > 0.0) {
        mean_u /= mass_total;
        for (int ix = 0; ix < grid.nx_local; ++ix) {
            const int sx = ng + ix;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double m = species.f[idx3(sx, iv, imu)];
                    const double u = species.cgrid.upar_cells[iv];
                    variance_u += m * (u - mean_u) * (u - mean_u);
                }
            }
        }
        variance_u /= mass_total;
    }
}

double entropy(const Species& species, const SpatialGrid& grid)
{
    const int ng = grid.nghost;
    double s = 0.0;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = ng + ix;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double m = species.f[idx3(sx, iv, imu)];
                if (m <= 0.0) continue;
                // The stored f is the folded cylindrical representation
                // f_tilde = 2 pi u_perp f_3D, so the correct entropy is
                //   S = -sum f_tilde ln f_tilde + sum f_tilde ln(2 pi u_perp).
                s -= m * std::log(m);
                const double ue = species.cgrid.uperp_cells[imu];
                if (ue > 0.0) s += m * std::log(2.0 * Const::pi * ue);
            }
        }
    }
    return s;
}

double l2_rel(const Species& a, const Species& b, const SpatialGrid& grid)
{
    const int ng = grid.nghost;
    double l2a = 0.0;
    double l2d = 0.0;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = ng + ix;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double va = a.f[idx3(sx, iv, imu)];
                const double vb = b.f[idx3(sx, iv, imu)];
                l2a += va * va;
                l2d += (va - vb) * (va - vb);
            }
        }
    }
    return std::sqrt(l2d) / std::max(1.0, std::sqrt(l2a));
}

void set_gaussian(Species& species, const SpatialGrid& grid,
                  double u0, double sigma_par, double sigma_perp,
                  double scale)
{
    const int ng = grid.nghost;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = ng + ix;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double up = species.cgrid.upar_cells[iv];
                const double ue = species.cgrid.uperp_cells[imu];
                // Folded cylindrical mass representation: f_tilde ~ u_perp
                // near the axis for a smooth 3-D distribution.  A nonzero
                // value at the axis is an unphysical ridge that the
                // cylindrical collision operator relaxes into negatives.
                species.f[idx3(sx, iv, imu)] =
                    scale * ue *
                    std::exp(-(up - u0) * (up - u0) /
                             (2.0 * sigma_par * sigma_par) -
                             ue * ue /
                             (2.0 * sigma_perp * sigma_perp));
            }
        }
    }
}

Metrics run_number_conservation()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species species = make_species(grid);
    set_gaussian(species, grid, 0.3, 0.15, 0.2, 1.0e18);
    const CylindricalCollisionCoefficients c =
        { 0.0, 0.0, 1.0e10, 0.0, 1.0e10 };
    const PrescribedCollisionCoefficients provider(c);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const CollisionDiagnostics d = collision.apply(species, grid, 0.0, 1.0e-16);
    m.mass_change = d.mass_change;
    double mean_u, variance_u, mass_total, energy_total;
    velocity_moments(species, grid, mean_u, variance_u, mass_total,
                     energy_total);
    m.mass_scale = mass_total;
    m.pass = d.success &&
             std::fabs(d.mass_change) <=
                 1.0e-10 * std::max(1.0, mass_total);
    return m;
}

Metrics run_drift_diffusion_moments()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species before = make_species(grid);
    set_gaussian(before, grid, 0.3, 0.15, 0.2, 1.0e18);
    Species after = before;
    const double A0 = 1.0e10;
    const double D0 = 1.0e10;
    const double dt = 1.0e-16;
    const CylindricalCollisionCoefficients c =
        { A0, 0.0, D0, 0.0, D0 };
    const PrescribedCollisionCoefficients provider(c);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const CollisionDiagnostics d = collision.apply(after, grid, 0.0, dt);
    double mean_b, var_b, mass_b, energy_b;
    double mean_a, var_a, mass_a, energy_a;
    velocity_moments(before, grid, mean_b, var_b, mass_b, energy_b);
    velocity_moments(after, grid, mean_a, var_a, mass_a, energy_a);
    const double expected_drift = A0 * dt;
    const double expected_var = 2.0 * D0 * dt;
    m.drift_rel_error =
        std::fabs((mean_a - mean_b) - expected_drift) /
        std::max(1.0e-30, std::fabs(expected_drift));
    m.var_rel_error =
        std::fabs((var_a - var_b) - expected_var) /
        std::max(1.0e-30, std::fabs(expected_var));
    m.pass = d.success && m.drift_rel_error < 0.05 &&
             m.var_rel_error < 0.05;
    return m;
}

Metrics run_maxwellian_equilibrium()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species before = make_species(grid);
    before.initialize_maxwellian();
    // Reduce the density so the relaxation rate nu0*dt stays small
    // (nu0 scales linearly with n): the closure then acts as a perturbation
    // and the Maxwellian must stay stationary to discretization error.
    const double scale = 1.0e-5;
    for (size_t i = 0; i < before.f.size(); ++i) before.f[i] *= scale;
    Species after = before;
    const MomentClosureCollisionCoefficients provider(20.0);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const CollisionDiagnostics d = collision.apply(after, grid, 0.0, 1.0e-16);
    m.equilibrium_l2_rel = l2_rel(before, after, grid);
    m.pass = d.success && m.equilibrium_l2_rel < 1.0e-6;
    return m;
}

Metrics run_h_theorem()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species species = make_species(grid);
    // Anisotropic initial distribution: u_parallel spread twice the
    // perpendicular spread (T_par > T_perp), which relaxes and isotropises.
    const int ng = grid.nghost;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = ng + ix;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double up = species.cgrid.upar_cells[iv];
                const double ue = species.cgrid.uperp_cells[imu];
                species.f[idx3(sx, iv, imu)] =
                    1.0e18 * ue *
                    std::exp(-up * up / (2.0 * 0.04) -
                             ue * ue / (2.0 * 0.01));
            }
        }
    }
    const double s_before = entropy(species, grid);
    const MomentClosureCollisionCoefficients provider(20.0);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const double dt = 1.0e-16;
    bool ok = true;
    for (int step = 0; step < 20; ++step) {
        const CollisionDiagnostics d = collision.apply(species, grid, 0.0, dt);
        ok = ok && d.success;
    }
    const double s_after = entropy(species, grid);
    m.entropy_change = s_after - s_before;
    m.pass = ok && m.entropy_change > 0.0;
    return m;
}

Metrics run_cross_diffusion()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species species = make_species(grid);
    set_gaussian(species, grid, 0.3, 0.15, 0.2, 1.0e18);
    double mean_u, variance_u, mass_total, energy_total;
    velocity_moments(species, grid, mean_u, variance_u, mass_total,
                     energy_total);
    const CylindricalCollisionCoefficients c =
        { 0.0, 0.0, 0.0, 1.0e9, 0.0 };
    const PrescribedCollisionCoefficients provider(c);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const CollisionDiagnostics d = collision.apply(species, grid, 0.0, 1.0e-16);
    double mean_u2, variance_u2, mass_total2, energy_total2;
    velocity_moments(species, grid, mean_u2, variance_u2, mass_total2,
                     energy_total2);
    m.cross_accepted = d.success && !d.unsupported_cross_diffusion;
    m.mass_change = d.mass_change;
    m.mass_scale = mass_total;
    m.pass = m.cross_accepted &&
             std::fabs(d.mass_change) <=
                 1.0e-10 * std::max(1.0, mass_total);
    return m;
}

Metrics run_reservoir_accounting()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species species = make_species(grid);
    set_gaussian(species, grid, 0.3, 0.15, 0.2, 1.0e18);
    const double ke_before = species.total_kinetic_energy();
    const MomentClosureCollisionCoefficients provider(20.0);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    const CollisionDiagnostics d = collision.apply(species, grid, 0.0, 1.0e-16);
    const double ke_after = species.total_kinetic_energy();
    m.reservoir_rel_error =
        std::fabs(d.reservoir_energy_change - (ke_before - ke_after)) /
        std::max(1.0, std::fabs(ke_before));
    m.pass = d.success && m.reservoir_rel_error < 1.0e-12;
    return m;
}

Metrics run_roundoff_negative_cleanup()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 16, 1.6 * Const::micro);
    Species species = make_species(grid);
    set_gaussian(species, grid, 0.0, 0.15, 0.2, 1.0e18);
    const int sx = grid.nghost;
    const int iv = Param::Nv / 2;
    const int imu = Param::Nmu - 1;
    const size_t tail = idx3(sx, iv, imu);
    species.f[tail] = -1.0e-108;
    const MomentClosureCollisionCoefficients provider(20.0);
    CylindricalFokkerPlanckCollision collision(provider,
        CollisionIntegratorType::BACKWARD_EULER);
    collision.set_bulk_integrator(BulkCollisionIntegrator::CHANG_COOPER_FLUX);
    const CollisionDiagnostics d = collision.apply(species, grid, 0.0, 1.0e-16);
    m.mass_change = d.mass_change;
    m.mass_scale = species.total_particle_number();
    m.pass = d.success && species.f[tail] >= 0.0;
    return m;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        // The test grid is a single-rank layout; the collision operator has
        // no cross-rank communication, so the module acceptance is single
        // rank.  Multi-rank collision-ledger consistency is covered by the
        // production collision short run (section 17.9).
        std::cerr << "cylindrical_fp_collision_test: single-rank test only\n";
        MPI_Finalize();
        return 2;
    }
    TestArgs args;
    const bool parsed = parse_args(argc, argv, args);
    if (!parsed) {
        std::cerr << "usage: cylindrical_fp_collision_test "
                     "[--case <case|all>] [--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;
    const char* names[7] = {
        "number-conservation",
        "drift-diffusion-moments",
        "maxwellian-equilibrium",
        "h-theorem",
        "cross-diffusion-supported",
        "reservoir-accounting",
        "roundoff-negative-cleanup"
    };
    Metrics (*runs[7])() = {
        &run_number_conservation,
        &run_drift_diffusion_moments,
        &run_maxwellian_equilibrium,
        &run_h_theorem,
        &run_cross_diffusion,
        &run_reservoir_accounting,
        &run_roundoff_negative_cleanup
    };
    for (int i = 0; i < 7; ++i) {
        const bool selected =
            args.test_case == "all" || args.test_case == names[i];
        if (!selected) continue;
        const Metrics m = runs[i]();
        if (!m.pass) pass = false;
        std::cout << names[i] << "=" << (m.pass ? 1 : 0)
                  << " mass_change=" << m.mass_change
                  << " drift_rel=" << m.drift_rel_error
                  << " var_rel=" << m.var_rel_error
                  << " eq_l2_rel=" << m.equilibrium_l2_rel
                  << " dS=" << m.entropy_change
                  << " reservoir_rel=" << m.reservoir_rel_error
                  << " cross=" << (m.cross_accepted ? 1 : 0) << "\n";
        if (!args.result_path.empty()) {
            std::ofstream out(args.result_path.c_str(), std::ios::app);
            if (out) {
                out << "case=" << names[i] << " pass=" << (m.pass ? 1 : 0)
                    << " mass_change=" << m.mass_change
                    << " drift_rel=" << m.drift_rel_error
                    << " var_rel=" << m.var_rel_error
                    << " eq_l2_rel=" << m.equilibrium_l2_rel
                    << " dS=" << m.entropy_change
                    << " reservoir_rel=" << m.reservoir_rel_error
                    << " cross=" << (m.cross_accepted ? 1 : 0) << "\n";
            }
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
