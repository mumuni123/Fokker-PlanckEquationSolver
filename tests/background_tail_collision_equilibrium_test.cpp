// Stage H8 acceptance: an isotropic equal-weight tail is a stationary state
// of the tail--tail Nanbu--Perez operator (Maxwellian equilibrium kept):
// after many collision calls the mean velocity stays ~0, the temperature
// <u^2> and the isotropy are preserved, and the weighted 3-momentum and
// relativistic energy are conserved.
//
// Usage:
//   background_tail_collision_equilibrium_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_collision.h"
#include "background_tail_nanbu_perez.h"
#include "grid.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
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

void totals(const BackgroundTailPIC& tail, double& px, double& py,
            double& pz, double& ke)
{
    px = py = pz = ke = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        px += Const::me * Const::c * p.weight * p.ux;
        py += Const::me * Const::c * p.weight * p.uy;
        pz += Const::me * Const::c * p.weight * p.uz;
        ke += Const::me * Const::c * Const::c * p.weight * (gamma - 1.0);
    }
}

void stats(const BackgroundTailPIC& tail, double& mean_ux, double& temp,
           double& anisotropy)
{
    double n = 0.0;
    double sx = 0.0;
    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        n += p.weight;
        sx += p.weight * p.ux;
        sxx += p.weight * p.ux * p.ux;
        syy += p.weight * p.uy * p.uy;
        szz += p.weight * p.uz * p.uz;
    }
    mean_ux = sx / n;
    temp = (sxx + syy + szz) / n;
    anisotropy = std::fabs(sxx / n - 0.5 * (syy + szz) / n) /
                 std::max(1.0e-30, temp);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_collision_equilibrium_test: "
                     "single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: background_tail_collision_equilibrium_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    std::uint64_t seed = 0x51a7;
    const double sigma = 0.6;
    for (int i = 0; i < 400; ++i) {
        std::uint64_t s = mix64(seed + static_cast<std::uint64_t>(i));
        const double u0 = 0.5 * std::sqrt(
            -2.0 * std::log(1.0e-15 + uniform01(s)));
        const double phi = 2.0 * Const::pi * uniform01(s);
        const double cz = 2.0 * uniform01(s) - 1.0;
        const double sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        p.ux = sigma * u0 * sz * std::cos(phi);
        p.uy = sigma * u0 * sz * std::sin(phi);
        p.uz = sigma * u0 * cz;
        p.weight = 1.0e20;
        p.id = static_cast<std::uint64_t>(i + 1);
        tail.particles.push_back(p);
    }
    double px0, py0, pz0, ke0;
    totals(tail, px0, py0, pz0, ke0);
    double mean0, temp0, aniso0;
    stats(tail, mean0, temp0, aniso0);

    TailCollisionRequest request;
    request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
    request.dt = 1.0e-13;
    request.coulomb_log = 20.0;
    request.rng_seed_base = 0x51a7;
    request.mpi_rank = 0;
    request.max_particle_growth = 0.0;
    bool ok = true;
    for (int step = 0; step < 50; ++step) {
        request.accepted_step = step + 1;
        request.collision_half = step % 2;
        TailCollisionDiagnostics diag;
        ok = ok && nanbu_perez_collide(tail, grid, request, diag);
    }
    double px1, py1, pz1, ke1;
    totals(tail, px1, py1, pz1, ke1);
    double mean1, temp1, aniso1;
    stats(tail, mean1, temp1, aniso1);
    const double p_scale = std::max(
        1.0, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
    const double p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
                          std::fabs(pz1 - pz0)) / p_scale;
    const double k_err = std::fabs(ke1 - ke0) / std::max(1.0, ke0);
    const double temp_err =
        std::fabs(temp1 - temp0) / std::max(1.0e-30, temp0);
    const bool equil = ok && p_err < 1.0e-12 && k_err < 1.0e-12 &&
                       temp_err < 0.02 && aniso1 < 0.05 &&
                       std::fabs(mean1) < 0.05;
    std::cout << "equilibrium=" << (equil ? 1 : 0)
              << " p_err=" << p_err << " k_err=" << k_err
              << " temp_err=" << temp_err
              << " aniso=" << aniso0 << "->" << aniso1
              << " mean=" << mean0 << "->" << mean1 << "\n";
    pass = pass && equil;

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=background-tail-collision-equilibrium pass="
                << (pass ? 1 : 0) << " temp_err=" << temp_err
                << " aniso_final=" << aniso1 << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
