// Stage H8 acceptance: HybridCollisionStep is transactional (section 10.4).
// A successful half-step mutates both trials; a failed half-step (an
// infeasible field-particle reaction on a massless bulk) leaves both the
// bulk trial and the tail trial bitwise unchanged.
//
// Usage:
//   hybrid_collision_transaction_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "grid.h"
#include "hybrid_collision_step.h"
#include "species.h"

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

bool particles_equal(const BackgroundTailPIC& a, const BackgroundTailPIC& b)
{
    if (a.particles.size() != b.particles.size()) return false;
    for (size_t i = 0; i < a.particles.size(); ++i) {
        if (a.particles[i].x != b.particles[i].x ||
            a.particles[i].ux != b.particles[i].ux ||
            a.particles[i].uy != b.particles[i].uy ||
            a.particles[i].uz != b.particles[i].uz ||
            a.particles[i].weight != b.particles[i].weight ||
            a.particles[i].id != b.particles[i].id) {
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "hybrid_collision_transaction_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: hybrid_collision_transaction_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    CylindricalCollisionCoefficients c =
        { 1.0e10, 0.0, 1.0e6, 0.0, 1.0e6 };
    PrescribedCollisionCoefficients provider(c);

    // Failure case: a massless bulk cannot absorb the reaction, so the
    // half-step must fail and leave both trials unchanged.
    {
        Species bulk;
        bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                  -Const::qe, Const::me, Param::dens, Param::temperature_e,
                  false, grid);
        for (size_t i = 0; i < bulk.f.size(); ++i) bulk.f[i] = 0.0;
        BackgroundTailPIC tail;
        tail.init(grid);
        for (int i = 0; i < 50; ++i) {
            BackgroundTailParticle p;
            p.x = 0.4e-6;
            p.ux = 0.5;
            p.uy = 0.1;
            p.uz = 0.05;
            p.weight = 1.0e14;
            p.id = static_cast<std::uint64_t>(i + 1);
            tail.particles.push_back(p);
        }
        const Species bulk_before = bulk;
        const BackgroundTailPIC tail_before = tail;
        HybridCollisionConfig config;
        config.bulk_provider = &provider;
    config.requested_kernel = TailCollisionKernel::KramersMoyalSDE;
    config.tail_bulk_kernel = TailCollisionKernel::KramersMoyalSDE;
    config.pairs.bulk_bulk = false;
    config.pairs.bulk_tail = true;
    config.pairs.tail_tail = false;
        config.pairs.tail_bulk = true;
        config.dt = 1.0e-12;   // large tail delta -> infeasible reaction
        config.accepted_step = 1;
        config.collision_half = 0;
        config.rng_seed_base = 0x1234;
        config.mpi_rank = 0;
        HybridCollisionStep step;
        const HybridCollisionDiagnostics d = step.advance(bulk, tail, grid,
                                                          config);
        bool bulk_same = bulk.f.size() == bulk_before.f.size();
        for (size_t i = 0; i < bulk.f.size() && bulk_same; ++i) {
            if (bulk.f[i] != bulk_before.f[i]) bulk_same = false;
        }
        const bool unchanged = !d.success && bulk_same &&
                               particles_equal(tail, tail_before);
        std::cout << "transaction-failure=" << (unchanged ? 1 : 0)
                  << " success=" << (d.success ? 1 : 0) << "\n";
        pass = pass && unchanged;
    }

    // Success case: a finite bulk absorbs the reaction; the trials change.
    {
        Species bulk;
        bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                  -Const::qe, Const::me, Param::dens, Param::temperature_e,
                  false, grid);
        bulk.initialize_maxwellian();
        const double scale = 1.0e-2;
        for (size_t i = 0; i < bulk.f.size(); ++i) bulk.f[i] *= scale;
        BackgroundTailPIC tail;
        tail.init(grid);
        for (int i = 0; i < 50; ++i) {
            BackgroundTailParticle p;
            p.x = 0.4e-6;
            p.ux = 0.5;
            p.uy = 0.1;
            p.uz = 0.05;
            p.weight = 1.0e14;
            p.id = static_cast<std::uint64_t>(i + 1);
            tail.particles.push_back(p);
        }
        HybridCollisionConfig config;
        config.bulk_provider = &provider;
    config.requested_kernel = TailCollisionKernel::KramersMoyalSDE;
    config.tail_bulk_kernel = TailCollisionKernel::KramersMoyalSDE;
    config.pairs.bulk_bulk = false;
    config.pairs.bulk_tail = true;
    config.pairs.tail_tail = false;
        config.pairs.tail_bulk = true;
        config.dt = 1.0e-16;
        config.accepted_step = 1;
        config.collision_half = 0;
        config.rng_seed_base = 0x5678;
        config.mpi_rank = 0;
        HybridCollisionStep step;
        const HybridCollisionDiagnostics d = step.advance(bulk, tail, grid,
                                                          config);
        const bool mutated = d.success && d.reaction_cells >= 1;
        std::cout << "transaction-success=" << (mutated ? 1 : 0)
                  << " reaction_cells=" << d.reaction_cells << "\n";
        pass = pass && mutated;
    }

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=hybrid-collision-transaction pass="
                << (pass ? 1 : 0) << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
