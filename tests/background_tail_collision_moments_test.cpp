// Stage H8 acceptance: Kramers-Moyal SDE weak moments (section 10.3.2).
// With a prescribed constant drift/diffusion, a large sample of particles at
// one fixed velocity point must show the first moment <du> = A dt and the
// second moments <du_i du_j> - <du_i><du_j> = D_ij dt to sampling accuracy,
// and dt/dt/2 must give the same first moment (linear scaling).
//
// Usage:
//   background_tail_collision_moments_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_collision.h"
#include "background_tail_collision_sde.h"
#include "grid.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string result_path;
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

struct Moments {
    double mean_ux;
    double var_ux;
    double perp_var;
};

Moments run_sde(BackgroundTailPIC& tail, const SpatialGrid& grid,
                const std::vector<LocalCollisionMoments>& moments,
                const CollisionCoefficientProvider& provider, double dt,
                std::uint64_t seed)
{
    TailCollisionRequest request;
    request.kernel = TailCollisionKernel::KramersMoyalSDE;
    request.dt = dt;
    request.accepted_step = 3;
    request.collision_half = 0;
    request.rng_seed_base = seed;
    request.mpi_rank = 0;
    TailCollisionDiagnostics diag;
    sde_collide(tail, grid, request, provider, moments, diag);
    double mean = 0.0;
    double var = 0.0;
    double perp = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        mean += p.ux;
        var += p.ux * p.ux;
        perp += p.uy * p.uy + p.uz * p.uz;
    }
    const double n = static_cast<double>(tail.particles.size());
    mean /= n;
    var = var / n - mean * mean;
    perp /= n;
    Moments m = { mean, var, perp };
    return m;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_collision_moments_test: single-rank "
                     "only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: background_tail_collision_moments_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    const double A0 = 1.0e12;
    const double D0 = 1.0e10;
    const CylindricalCollisionCoefficients c =
        { A0, 0.0, D0, 0.0, D0 };
    const PrescribedCollisionCoefficients provider(c);
    std::vector<LocalCollisionMoments> moments(
        static_cast<size_t>(grid.nx_local));
    const double dt = 1.0e-15;
    const double expected_drift = A0 * dt;
    const double expected_var = D0 * dt;

    const int n_particles = 100000;
    BackgroundTailPIC reference;
    reference.init(grid);
    for (int i = 0; i < n_particles; ++i) {
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        p.ux = 0.5;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0e18;
        p.id = static_cast<std::uint64_t>(i + 1);
        reference.particles.push_back(p);
    }
    const double u0 = 0.5;
    // Particles start on the axis (u_perp = 0): the 2 u_perp * du cross
    // terms vanish and the perpendicular diffusion D_perp_perp dt is
    // measured cleanly (no sampling-noise mask).
    const double perp0 = 0.0;
    BackgroundTailPIC full = reference;
    const Moments m_full = run_sde(full, grid, moments, provider, dt,
                                   0xabc);
    const double drift_err =
        std::fabs((m_full.mean_ux - u0) - expected_drift) /
        std::max(1.0e-30, expected_drift);
    const double var_err =
        std::fabs(m_full.var_ux - expected_var) /
        std::max(1.0e-30, expected_var);
    const double perp_err =
        std::fabs((m_full.perp_var - perp0) - expected_var) /
        std::max(1.0e-30, expected_var);
    const bool first_order =
        (m_full.mean_ux - u0) > 0.0 && drift_err < 0.2 && var_err < 0.2 &&
        perp_err < 0.2;
    std::cout << "sde-moments=" << (first_order ? 1 : 0)
              << " mean=" << m_full.mean_ux
              << " expected=" << expected_drift
              << " drift_err=" << drift_err
              << " var_err=" << var_err
              << " perp_err=" << perp_err << "\n";
    pass = pass && first_order;

    // dt/dt/2: two half steps give the same first moment.
    BackgroundTailPIC half;
    half.init(grid);
    for (int i = 0; i < n_particles; ++i) {
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        p.ux = 0.5;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0e18;
        p.id = static_cast<std::uint64_t>(i + 1);
        half.particles.push_back(p);
    }
    const Moments m_half1 = run_sde(half, grid, moments, provider,
                                    0.5 * dt, 0xdef);
    (void)m_half1;
    const Moments m_half2 = run_sde(half, grid, moments, provider,
                                    0.5 * dt, 0xdef + 1);
    // Two half-steps must give the same total first moment as one full step.
    const double half_drift = m_half2.mean_ux - u0;
    const double half_err =
        std::fabs(half_drift - expected_drift) /
        std::max(1.0e-30, expected_drift);
    const bool dt_conv = half_err < 0.2;
    std::cout << "sde-dt-half=" << (dt_conv ? 1 : 0)
              << " half_drift=" << half_drift
              << " half_err=" << half_err << "\n";
    pass = pass && dt_conv;

    // Regression for H9 step 413: away from the axis the Cartesian tensor is
    // rank deficient.  A nearly y-aligned radial direction must not make its
    // positive-semidefinite cylindrical diffusion fail factorization.
    BackgroundTailPIC near_axis;
    near_axis.init(grid);
    for (int i = 0; i < 4; ++i) {
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        p.ux = 12.0;
        p.uy = 1.0e-12;
        p.uz = 1.0;
        p.weight = 1.0;
        p.id = static_cast<std::uint64_t>(200001 + i);
        near_axis.particles.push_back(p);
    }
    TailCollisionRequest near_axis_request;
    near_axis_request.kernel = TailCollisionKernel::KramersMoyalSDE;
    near_axis_request.dt = dt;
    near_axis_request.accepted_step = 413;
    near_axis_request.collision_half = 0;
    near_axis_request.rng_seed_base = 0x413ULL;
    near_axis_request.mpi_rank = 3;
    TailCollisionDiagnostics near_axis_diag;
    const bool near_axis_ok = sde_collide(
        near_axis, grid, near_axis_request, provider, moments, near_axis_diag);
    bool near_axis_finite = near_axis_ok;
    for (size_t i = 0; i < near_axis.particles.size(); ++i) {
        const BackgroundTailParticle& p = near_axis.particles[i];
        near_axis_finite = near_axis_finite && std::isfinite(p.ux) &&
                           std::isfinite(p.uy) && std::isfinite(p.uz);
    }
    std::cout << "sde-rank-deficient=" << (near_axis_finite ? 1 : 0)
              << " failure_reason="
              << tail_collision_failure_reason_name(
                     near_axis_diag.failure_reason)
              << "\n";
    pass = pass && near_axis_finite;

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=background-tail-collision-moments pass="
                << (pass ? 1 : 0)
                << " drift_err=" << drift_err
                << " var_err=" << var_err
                << " perp_err=" << perp_err << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
