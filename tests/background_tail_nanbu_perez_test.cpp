// Stage H8 acceptance: tail--tail Nanbu--Perez equal-weight scattering
// (section 10.3.1).  A relativistic two-body COM update must conserve the
// weighted 3-momentum and relativistic energy to summation accuracy per
// event, keep the macro-particle count (equal weights -> no split), and a
// batch must show the same conservation globally.
//
// Usage:
//   background_tail_nanbu_perez_test [--result <path>]
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

void add_particle(BackgroundTailPIC& tail, double x, double ux, double uy,
                  double uz, double weight, std::uint64_t id)
{
    BackgroundTailParticle p;
    p.x = x;
    p.ux = ux;
    p.uy = uy;
    p.uz = uz;
    p.weight = weight;
    p.id = id;
    tail.particles.push_back(p);
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
        std::cerr << "background_tail_nanbu_perez_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: background_tail_nanbu_perez_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);

    // Single relativistic pair (equal weights).
    {
        BackgroundTailPIC tail;
        tail.init(grid);
        add_particle(tail, 0.4e-6, 12.0, 2.0, 1.0, 1.0e20, 1);
        add_particle(tail, 0.4e-6, -6.0, -3.0, 0.5, 1.0e20, 2);
        double px0, py0, pz0, ke0;
        totals(tail, px0, py0, pz0, ke0);
        TailCollisionRequest request;
        request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
        request.dt = 1.0e-16;
        request.accepted_step = 1;
        request.collision_half = 0;
        request.coulomb_log = 20.0;
        request.rng_seed_base = 0x12345678;
        request.mpi_rank = 0;
        request.max_particle_growth = 0.0;
        TailCollisionDiagnostics diag;
        const bool ok = nanbu_perez_collide(tail, grid, request, diag);
        double px1, py1, pz1, ke1;
        totals(tail, px1, py1, pz1, ke1);
        const double p_scale = std::max(
            1.0, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
        const double p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
                              std::fabs(pz1 - pz0)) / p_scale;
        const double k_err = std::fabs(ke1 - ke0) / std::max(1.0, ke0);
        const bool conserved = ok && p_err < 1.0e-12 && k_err < 1.0e-12 &&
                               diag.weight_split_count == 0 &&
                               tail.particles.size() == 2;
        std::cout << "single-pair=" << (conserved ? 1 : 0)
                  << " p_err=" << p_err << " k_err=" << k_err
                  << " split=" << diag.weight_split_count << "\n";
        pass = pass && conserved;
    }

    // Batch of 200 equal-weight particles: global conservation.
    {
        BackgroundTailPIC tail;
        tail.init(grid);
        for (int i = 0; i < 200; ++i) {
            const double ux = 4.0 * std::sin(0.31 * i) + 2.0;
            const double uy = 3.0 * std::cos(0.27 * i);
            const double uz = 2.0 * std::sin(0.13 * i + 1.0);
            add_particle(tail, 0.4e-6, ux, uy, uz, 1.0e19, 1000 + i);
        }
        double px0, py0, pz0, ke0;
        totals(tail, px0, py0, pz0, ke0);
        TailCollisionRequest request;
        request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
        request.dt = 1.0e-15;
        request.accepted_step = 2;
        request.collision_half = 1;
        request.coulomb_log = 20.0;
        request.rng_seed_base = 0xdeadbeef;
        request.mpi_rank = 0;
        request.max_particle_growth = 0.0;
        TailCollisionDiagnostics diag;
        const bool ok = nanbu_perez_collide(tail, grid, request, diag);
        double px1, py1, pz1, ke1;
        totals(tail, px1, py1, pz1, ke1);
        const double p_scale = std::max(
            1.0, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
        const double p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
                              std::fabs(pz1 - pz0)) / p_scale;
        const double k_err = std::fabs(ke1 - ke0) / std::max(1.0, ke0);
        const bool conserved = ok && p_err < 1.0e-12 && k_err < 1.0e-12 &&
                               tail.particles.size() == 200 &&
                               diag.collision_substeps >= 1;
        std::cout << "batch-200=" << (conserved ? 1 : 0)
                  << " p_err=" << p_err << " k_err=" << k_err
                  << " substeps=" << diag.collision_substeps
                  << " max_s12=" << diag.max_s12 << "\n";
        pass = pass && conserved;
    }

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=background-tail-nanbu-perez pass=" << (pass ? 1 : 0)
                << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
