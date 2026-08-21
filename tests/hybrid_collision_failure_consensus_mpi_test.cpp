// Multi-rank regression for the production collision failure handshake.
// Rank 0 owns a deliberately dense tail particle that exceeds the configured
// Nanbu substep limit; all other ranks have an empty tail.  Every rank must
// return the same rejected VpfpStepResult instead of entering later MPI
// collectives with divergent control flow.

#include "background_tail_pic.h"
#include "beam_pic.h"
#include "bulk_tail_converter.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <cstdint>
#include <iostream>
#include <vector>

class ConversionConsensusTestAccess {
public:
    static bool synchronize(VpfpIntegrator& integrator, int local_reason,
                            const BulkTailConversionDiagnostics& diagnostics,
                            int rank, VpfpStepResult& result)
    {
        return integrator.synchronize_conversion_outcome(
            local_reason, diagnostics, rank, 10, result);
    }
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        if (rank == 0) {
            std::cerr << "hybrid_collision_failure_consensus_mpi_test: "
                         "single-rank smoke only; use at least 2 ranks for "
                         "the control-flow consensus gate\n";
        }
    }

    SpatialGrid grid;
    grid.init_with_domain(rank, size, 8 * size, 0.8 * Const::micro);
    Species electrons;
    electrons.init("background_electrons",
                   SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                   Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();

    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.right_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.left_reservoir =
        { Param::dens, Param::temperature_e, 0.0 };
    boundary_config.right_reservoir =
        { Param::dens, Param::temperature_e, 0.0 };
    OpenBackgroundBoundary boundary(boundary_config);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid,
                      { ElectrostaticBoundaryType::DIRICHLET_PHI,
                        0.0, 0.0, 0.0 });
    MomentClosureCollisionCoefficients provider(20.0);
    CylindricalFokkerPlanckCollision collision(
        provider, CollisionIntegratorType::BACKWARD_EULER);
    HybridVelocityPartition partition;
    partition.init(electrons.cgrid, 6.0, 1.0, 4, 4);
    BulkTailConverter converter;
    VpfpIntegrator integrator(boundary, field_solver, collision,
                              partition, converter, true);
    integrator.init(grid);
    integrator.set_beam_enabled(true);
    integrator.set_tail_collision(
        TailCollisionKernel::CoulombLandauNanbuPerez, 20.0,
        TailCollisionWeightMode::VirtualSplit, 1, 1.0);

    if (rank == 0) {
        BackgroundTailParticle particle;
        particle.x = grid.x(grid.nghost);
        particle.ux = 12.0;
        particle.uy = 0.1;
        particle.uz = 0.05;
        // n_tail ~ 1e31 m^-3 locally: enough to require more than one
        // Nanbu substep at the production dt, while avoiding integer
        // overflow in the required-substep diagnostic.
        particle.weight = 1.0e31 * grid.dx;
        particle.id = 1;
        integrator.tail_state().particles.push_back(particle);
    }

    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    std::vector<double> ion_density(static_cast<size_t>(grid.nx_local),
                                    Param::dens);
    const double dt = Param::dt_multiplier / Param::omega_pe;
    const VpfpStepResult result = integrator.advance(
        electrons, beam, fields, ion_density, 0.0, dt, rank, size);

    const int local_accepted = result.accepted ? 1 : 0;
    int accepted_min = 0;
    int accepted_max = 0;
    int failure_min = 0;
    int failure_max = 0;
    MPI_Allreduce(&local_accepted, &accepted_min, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_accepted, &accepted_max, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&result.failure_code, &failure_min, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&result.failure_code, &failure_max, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    const bool collision_consensus = accepted_min == 0 && accepted_max == 0 &&
                                     failure_min == 5 && failure_max == 5;

    // A rank-local conversion failure must be synchronized before peers can
    // enter the following tail migration collective.  Reason 2 is the
    // production nonfinite conversion classification; only rank 0 fails.
    BulkTailConversionDiagnostics conversion;
    conversion.finite = rank != 0;
    conversion.complete = rank != 0;
    conversion.conservative = rank != 0;
    conversion.fidelity_ok = rank != 0;
    VpfpStepResult conversion_result;
    const bool conversion_returned_ok =
        ConversionConsensusTestAccess::synchronize(
            integrator, rank == 0 ? 2 : 0, conversion, rank,
            conversion_result);
    const int local_conversion_ok = conversion_returned_ok ? 1 : 0;
    int conversion_ok_min = 0;
    int conversion_ok_max = 0;
    int conversion_failure_min = 0;
    int conversion_failure_max = 0;
    int conversion_rank_min = 0;
    int conversion_rank_max = 0;
    MPI_Allreduce(&local_conversion_ok, &conversion_ok_min, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_conversion_ok, &conversion_ok_max, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&conversion_result.failure_code, &conversion_failure_min,
                  1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&conversion_result.failure_code, &conversion_failure_max,
                  1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&conversion_result.failing_rank, &conversion_rank_min,
                  1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&conversion_result.failing_rank, &conversion_rank_max,
                  1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    const bool conversion_consensus =
        conversion_ok_min == 0 && conversion_ok_max == 0 &&
        conversion_failure_min == 10 && conversion_failure_max == 10 &&
        conversion_rank_min == 0 && conversion_rank_max == 0;
    const bool pass = collision_consensus && conversion_consensus;
    if (rank == 0) {
        std::cout << "accepted_min=" << accepted_min
                  << " accepted_max=" << accepted_max
                  << " failure_min=" << failure_min
                  << " failure_max=" << failure_max << "\n";
        std::cout << "conversion_ok_min=" << conversion_ok_min
                  << " conversion_ok_max=" << conversion_ok_max
                  << " conversion_failure_min=" << conversion_failure_min
                  << " conversion_failure_max=" << conversion_failure_max
                  << " conversion_rank_min=" << conversion_rank_min
                  << " conversion_rank_max=" << conversion_rank_max << "\n";
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
