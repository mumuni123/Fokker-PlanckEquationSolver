// Stage H8 acceptance (section 10.3.6 item 2): two equal-weight PIC test
// populations with different initial temperatures equilibrate through the
// tail--tail Nanbu--Perez kernel; the total weighted 3-momentum and
// relativistic energy are conserved and the final common temperature lies
// between the two initial values.
//
// Usage:
//   tail_collision_two_population_test [--result <path>]
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

double add_gaussian_population(BackgroundTailPIC& tail, double sigma,
                               int n, std::uint64_t seed_base,
                               std::uint64_t id_start)
{
    double u2_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        std::uint64_t s = mix64(seed_base + static_cast<std::uint64_t>(i));
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
        // Unique ids across populations (the collision pairs skip equal ids).
        p.id = id_start + static_cast<std::uint64_t>(i);
        tail.particles.push_back(p);
        u2_sum += p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;
    }
    return u2_sum / static_cast<double>(n);   // actual sample <u^2>
}

void population_temps(const BackgroundTailPIC& tail, int n_cold,
                      double& t_cold, double& t_hot)
{
    t_cold = t_hot = 0.0;
    double n_c = 0.0;
    double n_h = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double u2 = p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;
        // The collision substeps rebuild the particle list in cell order, so
        // the populations are tracked by their creation IDs (cold = 1..N).
        if (p.id <= static_cast<std::uint64_t>(n_cold)) {
            t_cold += u2;
            n_c += 1.0;
        } else {
            t_hot += u2;
            n_h += 1.0;
        }
    }
    t_cold /= std::max(1.0, n_c);
    t_hot /= std::max(1.0, n_h);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "tail_collision_two_population_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: tail_collision_two_population_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    const int n_cold = 200;
    const int n_hot = 200;
    const double t_cold0 =
        add_gaussian_population(tail, 0.3, n_cold, 0xaaaa, 1);
    const double t_hot0 =
        add_gaussian_population(tail, 1.2, n_hot, 0xbbbb, n_cold + 1);
    double px0, py0, pz0, ke0;
    px0 = py0 = pz0 = ke0 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        px0 += Const::me * Const::c * p.weight * p.ux;
        py0 += Const::me * Const::c * p.weight * p.uy;
        pz0 += Const::me * Const::c * p.weight * p.uz;
        ke0 += Const::me * Const::c * Const::c * p.weight * (gamma - 1.0);
    }

    TailCollisionRequest request;
    request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
    request.dt = 1.0e-13;
    request.coulomb_log = 20.0;
    request.rng_seed_base = 0xcccc;
    request.mpi_rank = 0;
    request.max_particle_growth = 0.0;
    bool ok = true;
    for (int step = 0; step < 80; ++step) {
        request.accepted_step = step + 1;
        request.collision_half = step % 2;
        TailCollisionDiagnostics diag;
        ok = ok && nanbu_perez_collide(tail, grid, request, diag);
    }
    double t_cold, t_hot;
    population_temps(tail, n_cold, t_cold, t_hot);
    double px1, py1, pz1, ke1;
    px1 = py1 = pz1 = ke1 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        px1 += Const::me * Const::c * p.weight * p.ux;
        py1 += Const::me * Const::c * p.weight * p.uy;
        pz1 += Const::me * Const::c * p.weight * p.uz;
        ke1 += Const::me * Const::c * Const::c * p.weight * (gamma - 1.0);
    }
    const double p_scale = std::max(
        1.0, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
    const double p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
                          std::fabs(pz1 - pz0)) / p_scale;
    const double k_err = std::fabs(ke1 - ke0) / std::max(1.0, ke0);
    const double t_common = 0.5 * (t_cold0 + t_hot0);
    const bool equil = ok && p_err < 1.0e-12 && k_err < 1.0e-12 &&
                       std::fabs(t_cold - t_hot) < 0.10 * (t_hot0 - t_cold0) &&
                       std::fabs(t_cold - t_common) < 0.10 * (t_hot0 - t_cold0) &&
                       std::fabs(t_hot - t_common) < 0.10 * (t_hot0 - t_cold0);
    std::cout << "two-population=" << (equil ? 1 : 0)
              << " t_cold=" << t_cold0 << "->" << t_cold
              << " t_hot=" << t_hot0 << "->" << t_hot
              << " common=" << t_common
              << " p_err=" << p_err << " k_err=" << k_err << "\n";
    pass = pass && equil;

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=tail-collision-two-population pass="
                << (pass ? 1 : 0) << " t_cold=" << t_cold
                << " t_hot=" << t_hot
                << " p_err=" << p_err << " k_err=" << k_err << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
