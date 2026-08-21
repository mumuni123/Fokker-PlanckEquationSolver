// Stage H8 acceptance: unequal macro-particle weights use the bounded
// Sentoku--Kemp correction.  The representation count must remain fixed,
// weighted relativistic energy is conserved eventwise, and momentum is
// conserved statistically without materialising residual-weight particles.
//
// Usage:
//   background_tail_collision_weight_test [--result <path>]
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

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_collision_weight_test: single-rank "
                     "only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: background_tail_collision_weight_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    {
        BackgroundTailParticle p1;
        p1.x = 0.4e-6;
        p1.ux = 10.0;
        p1.uy = 1.0;
        p1.uz = 0.5;
        p1.weight = 1.0e20;
        p1.id = 1;
        tail.particles.push_back(p1);
        BackgroundTailParticle p2;
        p2.x = 0.4e-6;
        p2.ux = -5.0;
        p2.uy = -2.0;
        p2.uz = 1.0;
        p2.weight = 3.0e20;
        p2.id = 2;
        tail.particles.push_back(p2);
    }
    double px0, py0, pz0, ke0;
    totals(tail, px0, py0, pz0, ke0);

    TailCollisionRequest request;
    request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
    request.weight_mode = TailCollisionWeightMode::VirtualSplit;
    // One substep: dt small enough that s12 stays below the substep bound.
    request.dt = 1.0e-16;
    request.accepted_step = 1;
    request.collision_half = 0;
    request.coulomb_log = 20.0;
    request.rng_seed_base = 0xfeed;
    request.mpi_rank = 0;
    request.max_particle_growth = 10.0;
    TailCollisionDiagnostics diag;
    const bool ok = nanbu_perez_collide(tail, grid, request, diag);

    double px1, py1, pz1, ke1;
    totals(tail, px1, py1, pz1, ke1);
    const double p_scale = std::max(
        1.0e-300, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
    const double p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
                          std::fabs(pz1 - pz0)) / p_scale;
    const double k_err = std::fabs(ke1 - ke0) /
                         std::max(1.0e-300, std::fabs(ke0));
    const bool split_ok = ok && diag.weight_split_count == 1 &&
                          tail.particles.size() == 2 &&
                          k_err < 5.0e-14;
    std::cout << "weight-correction=" << (split_ok ? 1 : 0)
              << " particles=" << tail.particles.size()
              << " split_count=" << diag.weight_split_count
              << " p_err=" << p_err << " k_err=" << k_err << "\n";
    pass = pass && split_ok;

    // Zero growth budget is now a valid production mode: unequal weights do
    // not increase the persistent macro-particle count.
    {
        BackgroundTailPIC tail2;
        tail2.init(grid);
        BackgroundTailParticle p1;
        p1.x = 0.4e-6;
        p1.ux = 10.0;
        p1.uy = 1.0;
        p1.uz = 0.5;
        p1.weight = 1.0e20;
        p1.id = 1;
        tail2.particles.push_back(p1);
        BackgroundTailParticle p2;
        p2.x = 0.4e-6;
        p2.ux = -5.0;
        p2.uy = -2.0;
        p2.uz = 1.0;
        p2.weight = 3.0e20;
        p2.id = 2;
        tail2.particles.push_back(p2);
        TailCollisionRequest r2 = request;
        r2.max_particle_growth = 0.0;   // no growth allowed
        TailCollisionDiagnostics d2;
        const bool bounded = nanbu_perez_collide(tail2, grid, r2, d2);
        std::cout << "bounded-growth=" << (bounded ? 1 : 0)
                   << " particles_after=" << tail2.particles.size() << "\n";
        pass = pass && bounded && tail2.particles.size() == 2 &&
               d2.particle_count_attempted == 2;
    }

    // Sentoku--Kemp conserves unequal-weight momentum statistically rather
    // than event by event. Put one identical pair in each independent cell
    // and verify that the random transverse corrections have no macroscopic
    // bias while eventwise energy conservation remains at roundoff.
    double ensemble_p_err = 0.0;
    double ensemble_k_err = 0.0;
    {
        const int pair_count = 1024;
        SpatialGrid ensemble_grid;
        ensemble_grid.init_with_domain(
            0, 1, pair_count, pair_count * 0.1 * Const::micro);
        BackgroundTailPIC ensemble;
        ensemble.init(ensemble_grid);
        ensemble.particles.reserve(static_cast<size_t>(2 * pair_count));
        for (int ix = 0; ix < pair_count; ++ix) {
            BackgroundTailParticle p1;
            p1.x = (static_cast<double>(ix) + 0.5) * 0.1e-6;
            p1.ux = 10.0;
            p1.uy = 1.0;
            p1.uz = 0.5;
            p1.weight = 1.0e20;
            p1.id = static_cast<std::uint64_t>(2 * ix + 1);
            ensemble.particles.push_back(p1);
            BackgroundTailParticle p2;
            p2.x = p1.x;
            p2.ux = -5.0;
            p2.uy = -2.0;
            p2.uz = 1.0;
            p2.weight = 3.0e20;
            p2.id = static_cast<std::uint64_t>(2 * ix + 2);
            ensemble.particles.push_back(p2);
        }
        double epx0, epy0, epz0, eke0;
        totals(ensemble, epx0, epy0, epz0, eke0);
        TailCollisionRequest er = request;
        er.max_particle_growth = 0.0;
        er.rng_seed_base = 0x12345678ULL;
        TailCollisionDiagnostics ed;
        const bool ensemble_ok =
            nanbu_perez_collide(ensemble, ensemble_grid, er, ed);
        double epx1, epy1, epz1, eke1;
        totals(ensemble, epx1, epy1, epz1, eke1);
        const double ep_scale = std::max(
            1.0e-300, std::fabs(epx0) + std::fabs(epy0) +
                        std::fabs(epz0));
        ensemble_p_err =
            (std::fabs(epx1 - epx0) + std::fabs(epy1 - epy0) +
             std::fabs(epz1 - epz0)) / ep_scale;
        ensemble_k_err = std::fabs(eke1 - eke0) /
                         std::max(1.0e-300, std::fabs(eke0));
        const bool ensemble_pass = ensemble_ok &&
            ensemble.particles.size() == static_cast<size_t>(2 * pair_count) &&
            ed.weight_split_count == static_cast<std::uint64_t>(pair_count) &&
            ensemble_p_err < 3.0e-3 && ensemble_k_err < 5.0e-13;
        std::cout << "ensemble-statistical=" << (ensemble_pass ? 1 : 0)
                  << " particles=" << ensemble.particles.size()
                  << " p_err=" << ensemble_p_err
                  << " k_err=" << ensemble_k_err << "\n";
        pass = pass && ensemble_pass;
    }

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=background-tail-collision-weight pass="
                << (pass ? 1 : 0) << " p_err=" << p_err
                << " k_err=" << k_err
                << " ensemble_p_err=" << ensemble_p_err
                << " ensemble_k_err=" << ensemble_k_err << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
