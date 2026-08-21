// Stage H8 acceptance (sections 10.6/10.7): the tail--bulk SDE (C_tb) plus
// the explicit field-particle reaction must make the bulk momentum/energy
// change cancel the tail's change per cell (pair balance), with the bulk
// number conserved.  The bulk-bulk and tail-tail pairs are disabled to
// isolate the reaction.
//
// Usage:
//   hybrid_collision_pair_balance_test [--result <path>]
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

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "hybrid_collision_pair_balance_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: hybrid_collision_pair_balance_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    bulk.initialize_maxwellian();
    const double scale = 1.0e-2;
    for (size_t i = 0; i < bulk.f.size(); ++i) bulk.f[i] *= scale;
    bulk.compute_moments();
    BackgroundTailPIC tail;
    tail.init(grid);
    for (int i = 0; i < 100; ++i) {
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        p.ux = 0.5;
        p.uy = 0.1;
        p.uz = 0.05;
        p.weight = 1.0e14;
        p.id = static_cast<std::uint64_t>(i + 1);
        tail.particles.push_back(p);
    }
    const Species bulk_initial = bulk;
    const BackgroundTailPIC tail_initial = tail;
    const double bulk_n0 = bulk.total_particle_number();
    double tail_px0 = 0.0;
    double tail_ke0 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        tail_px0 += Const::me * Const::c * p.weight * p.ux;
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        tail_ke0 += Const::me * Const::c * Const::c * p.weight *
                    (gamma - 1.0);
    }

    CylindricalCollisionCoefficients c =
        { 1.0e10, 0.0, 1.0e6, 0.0, 1.0e6 };
    PrescribedCollisionCoefficients provider(c);
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
    config.rng_seed_base = 0xabcd;
    config.mpi_rank = 0;
    HybridCollisionStep step;
    const HybridCollisionDiagnostics d = step.advance(bulk, tail, grid, config);

    const double bulk_n1 = bulk.total_particle_number();
    double tail_px1 = 0.0;
    double tail_ke1 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        tail_px1 += Const::me * Const::c * p.weight * p.ux;
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        tail_ke1 += Const::me * Const::c * Const::c * p.weight *
                    (gamma - 1.0);
    }
    const bool balance_ok =
        d.success && d.reaction_px_balance < 1.0e-12 &&
        d.reaction_energy_balance < 1.0e-12 &&
        std::fabs(bulk_n1 - bulk_n0) / std::max(1.0, bulk_n0) < 1.0e-12 &&
        std::fabs(d.bulk_mass_change) / std::max(1.0, bulk_n0) < 1.0e-12;
    std::cout << "pair-balance=" << (balance_ok ? 1 : 0)
              << " px_balance=" << d.reaction_px_balance
              << " ke_balance=" << d.reaction_energy_balance
              << " tail_dpx=" << (tail_px1 - tail_px0)
              << " tail_dke=" << (tail_ke1 - tail_ke0)
              << " mass_rel=" << (bulk_n1 - bulk_n0) /
                                     std::max(1.0, bulk_n0) << "\n";
    pass = pass && balance_ok;

    // Pair-mask isolation: C_tb without C_bt changes only the tail. This is
    // an explicit trace-particle mode and must not be confused with the
    // closed production pair.
    Species trace_bulk = bulk_initial;
    BackgroundTailPIC trace_tail = tail_initial;
    HybridCollisionConfig trace_config = config;
    trace_config.pairs.bulk_tail = false;
    const HybridCollisionDiagnostics trace_d =
        step.advance(trace_bulk, trace_tail, grid, trace_config);
    bool trace_bulk_unchanged =
        trace_bulk.f.size() == bulk_initial.f.size();
    for (size_t i = 0; trace_bulk_unchanged && i < trace_bulk.f.size(); ++i) {
        trace_bulk_unchanged = trace_bulk.f[i] == bulk_initial.f[i];
    }
    const bool trace_ok =
        trace_d.success && trace_d.tail_bulk_applied &&
        !trace_d.bulk_reaction_applied && trace_bulk_unchanged &&
        std::fabs(trace_d.tail_energy_change) > 0.0;
    std::cout << "tail-bulk-only=" << (trace_ok ? 1 : 0)
              << " tail_dke=" << trace_d.tail_energy_change
              << " bulk_unchanged=" << (trace_bulk_unchanged ? 1 : 0)
              << "\n";
    pass = pass && trace_ok;

    // Full production registry: the requested Coulomb mode must execute the
    // bulk FP, tail-tail Nanbu, tail-bulk SDE and bulk reaction in the same
    // transactional half-step.
    Species full_bulk = bulk_initial;
    BackgroundTailPIC full_tail = tail_initial;
    HybridCollisionConfig full_config = config;
    full_config.requested_kernel =
        TailCollisionKernel::CoulombLandauNanbuPerez;
    full_config.tail_tail_kernel =
        TailCollisionKernel::CoulombLandauNanbuPerez;
    full_config.tail_bulk_kernel = TailCollisionKernel::KramersMoyalSDE;
    full_config.weight_mode = TailCollisionWeightMode::VirtualSplit;
    full_config.max_substeps = 64;
    full_config.pairs.bulk_bulk = true;
    full_config.pairs.bulk_tail = true;
    full_config.pairs.tail_bulk = true;
    full_config.pairs.tail_tail = true;
    const HybridCollisionDiagnostics full_d =
        step.advance(full_bulk, full_tail, grid, full_config);
    const bool full_ok =
        full_d.success && full_d.bulk_bulk_applied &&
        full_d.tail_tail_applied && full_d.tail_bulk_applied &&
        full_d.bulk_reaction_applied && full_d.tail_tail_diag.success &&
        full_d.tail_bulk_diag.success &&
        full_d.reaction_px_balance < 1.0e-12 &&
        full_d.reaction_energy_balance < 1.0e-12;
    std::cout << "full-registry=" << (full_ok ? 1 : 0)
              << " bb=" << (full_d.bulk_bulk_applied ? 1 : 0)
              << " tt=" << (full_d.tail_tail_applied ? 1 : 0)
              << " tb=" << (full_d.tail_bulk_applied ? 1 : 0)
              << " bt=" << (full_d.bulk_reaction_applied ? 1 : 0)
              << "\n";
    pass = pass && full_ok;

    // All pairs disabled is an exact no-op.
    Species none_bulk = bulk_initial;
    BackgroundTailPIC none_tail = tail_initial;
    HybridCollisionConfig none_config = config;
    none_config.pairs.bulk_bulk = false;
    none_config.pairs.bulk_tail = false;
    none_config.pairs.tail_bulk = false;
    none_config.pairs.tail_tail = false;
    const HybridCollisionDiagnostics none_d =
        step.advance(none_bulk, none_tail, grid, none_config);
    bool none_unchanged = none_d.success &&
        none_bulk.f.size() == bulk_initial.f.size() &&
        none_tail.particles.size() == tail_initial.particles.size();
    for (size_t i = 0; none_unchanged && i < none_bulk.f.size(); ++i) {
        none_unchanged = none_bulk.f[i] == bulk_initial.f[i];
    }
    for (size_t i = 0; none_unchanged && i < none_tail.particles.size(); ++i) {
        const BackgroundTailParticle& lhs = none_tail.particles[i];
        const BackgroundTailParticle& rhs = tail_initial.particles[i];
        none_unchanged = lhs.x == rhs.x && lhs.ux == rhs.ux &&
                         lhs.uy == rhs.uy && lhs.uz == rhs.uz &&
                         lhs.weight == rhs.weight && lhs.id == rhs.id;
    }
    std::cout << "all-pairs-off=" << (none_unchanged ? 1 : 0) << "\n";
    pass = pass && none_unchanged;

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=hybrid-collision-pair-balance pass="
                << (pass ? 1 : 0)
                << " px_balance=" << d.reaction_px_balance
                << " ke_balance=" << d.reaction_energy_balance
                << " trace_ok=" << (trace_ok ? 1 : 0)
                << " full_registry_ok=" << (full_ok ? 1 : 0)
                << " all_pairs_off_ok=" << (none_unchanged ? 1 : 0)
                << "\n";
        }
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
