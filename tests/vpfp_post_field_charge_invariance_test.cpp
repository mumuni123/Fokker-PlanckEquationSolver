// JC0 acceptance test (section 3): post-field charge invariance precheck gate.
//
// For a Bulk+Tail candidate that has completed the field-particle core, this
// test runs the production post-field operators -- the second collision half
// C2, the collision-face conversion and the H10 tail-to-bulk return -- through
// VpfpIntegrator::post_field_charge_invariance_transaction(), and verifies that
// the combined (bulk + tail) electron number is unchanged per cell, per rank
// and globally.  It covers no-collision, bulk-collision, hybrid-collision,
// return none/hysteretic, conversion-active and a combined case, plus the MPI
// multi-rank variant.
//
// No production collision, conversion or return formula is re-implemented
// here: every operator is the production implementation reached through the
// integrator's read-only JC0 hook.
//
// Usage:
//   vpfp_post_field_charge_invariance_test --case all   [--result <path>]
//   vpfp_post_field_charge_invariance_test --case mpi   [--result <path>]

#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "grid.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "tail_bulk_return_test_common.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
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
    return args.test_case == "all" || args.test_case == "mpi";
}

struct IntegratorStateSnapshot {
    long long step_count;
    VpfpCombinedChecksum checksum;
    VpfpTailCumulativeLedger tail_cum;
    double tail_particles_x_sum;
    double tail_particles_ux_sum;
    std::size_t tail_particle_count;
};

IntegratorStateSnapshot snapshot_integrator_state(VpfpIntegrator& integrator)
{
    IntegratorStateSnapshot s;
    s.step_count = integrator.step_count();
    s.checksum = integrator.combined_checksum();
    s.tail_cum = integrator.tail_cumulative();
    const BackgroundTailPIC& tail = integrator.tail_state();
    s.tail_particle_count = tail.particles.size();
    s.tail_particles_x_sum = 0.0;
    s.tail_particles_ux_sum = 0.0;
    for (const auto& p : tail.particles) {
        s.tail_particles_x_sum += p.x;
        s.tail_particles_ux_sum += p.ux;
    }
    return s;
}

bool compare_integrator_state(const IntegratorStateSnapshot& a,
                              const IntegratorStateSnapshot& b)
{
    if (a.step_count != b.step_count) return false;
    if (a.checksum.number != b.checksum.number) return false;
    if (a.checksum.kinetic_energy != b.checksum.kinetic_energy) return false;
    if (a.checksum.field_energy != b.checksum.field_energy) return false;
    if (a.tail_cum.conversion_number != b.tail_cum.conversion_number) return false;
    if (a.tail_cum.conversion_px != b.tail_cum.conversion_px) return false;
    if (a.tail_cum.conversion_energy != b.tail_cum.conversion_energy) return false;
    if (a.tail_cum.particles_created != b.tail_cum.particles_created) return false;
    if (a.tail_cum.return_number != b.tail_cum.return_number) return false;
    if (a.tail_cum.return_px != b.tail_cum.return_px) return false;
    if (a.tail_cum.return_energy != b.tail_cum.return_energy) return false;
    if (a.tail_cum.return_particles_removed != b.tail_cum.return_particles_removed) return false;
    if (a.tail_particle_count != b.tail_particle_count) return false;
    if (a.tail_particles_x_sum != b.tail_particles_x_sum) return false;
    if (a.tail_particles_ux_sum != b.tail_particles_ux_sum) return false;
    return true;
}

struct CaseOutcome {
    std::string name;
    bool operators_ok;
    bool c2_ok;
    bool conversion_ok;
    bool return_ok;
    int collision_half_calls;
    int conversion_calls;
    int return_calls;
    int expected_return_calls;
    bool stochastic_ran;
    bool expected_stochastic;
    double cell_residual_linf;
    double cell_tolerance;
    double rank_residual_linf;
    double rank_tolerance;
    double global_residual;
    double global_tolerance;
    double rho_rel_linf;
    double export_number;
    bool cell_pass;
    bool rank_pass;
    bool global_pass;
    bool rho_pass;
    bool integrator_state_unchanged;
    bool case_pass;
};

// Center of the rank's local physical cell range (ghost-independent).
double local_center_x(const SpatialGrid& grid)
{
    return grid.x_min +
        (grid.ix_start + 0.5 * static_cast<double>(grid.nx_local)) * grid.dx;
}

// For the conversion cases, seed a single bulk-owned velocity cell adjacent
// to a +u interface face so the drift/diffusion exports through exactly one
// collision face (four quadrature nodes -> one group, representable by the
// conversion compressor).  A production-density uniform bulk would export
// across every interface face, which the compressor cannot represent within
// max_supports and makes convert_flux_batch fail spuriously.
bool seed_conversion_bump(Species& bulk, const SpatialGrid& grid,
                          const HybridVelocityPartition& partition)
{
    int bump_iv = -1;
    int bump_imu = -1;
    for (size_t q = 0; q < partition.upar_interface_faces.size(); ++q) {
        const BulkTailInterfaceFace& face = partition.upar_interface_faces[q];
        if (face.outward_sign > 0 && face.bulk_iv >= 0 &&
            face.bulk_imu >= 0 && face.bulk_iv < Param::Nv &&
            face.bulk_imu < Param::Nmu) {
            bump_iv = face.bulk_iv;
            bump_imu = face.bulk_imu;
            break;
        }
    }
    if (bump_iv < 0) return false;
    std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
    for (int il = 0; il < grid.nx_local; ++il) {
        bulk.f[idx3(grid.nghost + il, bump_iv, bump_imu)] = 1.0e18;
    }
    bulk.compute_moments();
    return true;
}

// Test particle IDs live in a high-bit namespace disjoint from the production
// conversion ID scheme ((mpi_rank << 32) | per-rank counter).  This guarantees
// that conversion-created particles can never collide with the seeded test
// tail state (the production operator rejects duplicate IDs as a transaction
// failure).
std::uint64_t test_particle_id(int rank, std::uint64_t index)
{
    const std::uint64_t marker = 0x7E57CA11ULL;
    const std::uint64_t low = (static_cast<std::uint64_t>(rank) << 16) |
                              (index & 0xFFFFULL);
    return (marker << 32) | low;
}

// One ordinary in-domain tail particle (non-return, low energy).
void add_ordinary_tail_particle(BackgroundTailPIC& tail, double x,
                                int rank, std::uint64_t index)
{
    BackgroundTailParticle p;
    p.x = x;
    p.ux = 0.5;
    p.uy = 0.3;
    p.uz = 0.0;
    p.weight = 1.0e14;
    p.id = test_particle_id(rank, index);
    tail.particles.push_back(p);
}

// Return-eligible tail cloud (bulk-owned velocity cells below 5 MeV).
void add_return_cloud(BackgroundTailPIC& tail, const Species& bulk,
                      const HybridVelocityPartition& partition, double x,
                      int rank)
{
    const std::pair<int, int> slot =
        tail_return_test::safe_velocity_slot(bulk, partition);
    if (slot.first < 0) return;
    tail_return_test::add_representable_cloud(
        tail, bulk, partition, slot.first, slot.second, x, 1.0e14,
        test_particle_id(rank, 0), 0);
}

CaseOutcome run_case(int rank, int size, int case_id)
{
    CaseOutcome out;
    out.operators_ok = false;
    out.c2_ok = false;
    out.conversion_ok = false;
    out.return_ok = false;
    out.collision_half_calls = 0;
    out.conversion_calls = 0;
    out.return_calls = 0;
    out.expected_return_calls = 0;
    out.stochastic_ran = false;
    out.expected_stochastic = false;
    out.cell_residual_linf = 0.0;
    out.cell_tolerance = 0.0;
    out.rank_residual_linf = 0.0;
    out.rank_tolerance = 0.0;
    out.global_residual = 0.0;
    out.global_tolerance = 0.0;
    out.rho_rel_linf = 0.0;
    out.export_number = 0.0;
    out.cell_pass = false;
    out.rank_pass = false;
    out.global_pass = false;
    out.rho_pass = false;

    switch (case_id) {
        case 0: out.name = "no-collision"; break;
        case 1: out.name = "bulk-collision"; break;
        case 2: out.name = "hybrid-collision"; break;
        case 3: out.name = "return-none"; break;
        case 4: out.name = "return-hysteretic"; break;
        case 5: out.name = "conversion-active"; break;
        case 6: out.name = "combined"; break;
        default: out.name = "unknown"; return out;
    }

    const int nx_global = (size == 1) ? 32 : 8 * size;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, nx_global, Param::Lx);

    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);

    // Bulk state: a localized high-energy bump for the conversion cases (so
    // the collision-face export is a single, representable face); a physical
    // Maxwellian otherwise.
    if (case_id == 5 || case_id == 6) {
        if (!seed_conversion_bump(bulk, grid, partition)) {
            bulk.initialize_maxwellian();
        }
    } else {
        bulk.initialize_maxwellian();
    }

    // Tail particles.
    const double cx = local_center_x(grid);
    if (case_id == 2) {
        for (int n = 0; n < 3; ++n) {
            add_ordinary_tail_particle(
                tail, cx + 0.1 * static_cast<double>(n) * grid.dx,
                rank, static_cast<std::uint64_t>(n));
        }
    } else if (case_id == 3 || case_id == 4 || case_id == 6) {
        add_return_cloud(tail, bulk, partition, cx, rank);
    } else {
        add_ordinary_tail_particle(tail, cx, rank, 0);
    }

    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.right_type = BackgroundXBoundaryType::RESERVOIR;
    boundary_config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    boundary_config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    OpenBackgroundBoundary boundary(boundary_config);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid,
                      { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 });

    // Collision provider / operator selection.  The u_perp coefficients are
    // zero so the conversion bump only crosses the u_parallel interface face
    // (a single representable export) and never spreads across u_perp faces.
    const CylindricalCollisionCoefficients coeff = {
        1.0e8, 0.0, 1.0e8, 0.0, 0.0 };
    const PrescribedCollisionCoefficients prescribed_provider(coeff);
    const MomentClosureCollisionCoefficients moment_closure_provider(20.0);
    const ZeroCollisionCoefficients zero_provider;

    const bool zero_collision = (case_id == 0 || case_id == 3 || case_id == 4);
    const bool hybrid = (case_id == 2);
    out.expected_stochastic = hybrid;

    const CollisionCoefficientProvider& provider =
        zero_collision ? static_cast<const CollisionCoefficientProvider&>(zero_provider) :
        (hybrid ? static_cast<const CollisionCoefficientProvider&>(moment_closure_provider) :
                  static_cast<const CollisionCoefficientProvider&>(prescribed_provider));
    CylindricalFokkerPlanckCollision collision(
        provider, CollisionIntegratorType::BACKWARD_EULER);

    const bool conversion_active = (case_id == 5 || case_id == 6);
    if (conversion_active) {
        collision.set_bulk_integrator(BulkCollisionIntegrator::CHANG_COOPER_FLUX);
    }

    BulkTailConverter converter;
    VpfpIntegrator integrator(boundary, field_solver, collision,
                              partition, converter, true);
    integrator.init(grid);

    if (conversion_active) {
        integrator.set_tail_conversion_mode(TailConversionMode::FLUX_INTERFACE,
                                            4, 7, 0);
        integrator.set_collision_interface_exporting_absorbing(true);
        // Populate the bulk-owned velocity mask without enabling the tail
        // collision backend (hybrid stays off for the bulk-only case).
        integrator.set_tail_collision(TailCollisionKernel::None, 20.0,
                                      TailCollisionWeightMode::VirtualSplit,
                                      1, 1.0);
    }
    if (hybrid) {
        integrator.set_tail_collision(
            TailCollisionKernel::CoulombLandauNanbuPerez, 20.0,
            TailCollisionWeightMode::VirtualSplit, 1, 1.0);
    }
    if (case_id == 4 || case_id == 6) {
        TailBulkReturnConfig ret;
        ret.enabled = true;
        ret.return_energy_mev = 5.5;
        ret.residence_steps = 1;
        ret.max_stencil_radius = 3;
        ret.moment_tolerance = 1.0e-12;
        integrator.set_tail_bulk_return(ret);
        out.expected_return_calls = 1;
    }

    const double dt = 1.0e-15;
    IntegratorStateSnapshot snap_before = snapshot_integrator_state(integrator);

    VpfpPostFieldChargeInvarianceReport report;
    out.operators_ok = integrator.post_field_charge_invariance_transaction(
        bulk, tail, 0.0, dt, rank, size, report);

    IntegratorStateSnapshot snap_after = snapshot_integrator_state(integrator);

    out.c2_ok = report.c2_ok;
    out.conversion_ok = report.conversion_ok;
    out.return_ok = report.return_ok;
    out.collision_half_calls = report.collision_half_calls;
    out.conversion_calls = report.collision_face_conversion_calls;
    out.return_calls = report.tail_return_calls;
    out.stochastic_ran = report.tail_collision_stochastic_ran;
    out.cell_residual_linf = report.cell_number_residual_linf;
    out.cell_tolerance = report.cell_number_tolerance;
    out.rank_residual_linf = report.rank_number_residual_linf;
    out.rank_tolerance = report.rank_number_tolerance;
    out.global_residual = report.global_number_residual;
    out.global_tolerance = report.global_number_tolerance;
    out.rho_rel_linf = report.rho_before_after_relative_linf;
    out.export_number = report.collision_face_export_number;

    out.cell_pass = out.cell_residual_linf <= out.cell_tolerance;
    out.rank_pass = out.rank_residual_linf <= out.rank_tolerance;
    out.global_pass = out.global_residual <= out.global_tolerance;
    out.rho_pass = out.rho_rel_linf <= 1.0e-12;
    out.integrator_state_unchanged = compare_integrator_state(snap_before, snap_after);
    out.case_pass = out.cell_pass && out.rank_pass && out.global_pass &&
                    out.rho_pass && out.integrator_state_unchanged;
    return out;
}

struct GateMetrics {
    int case_count;
    double cell_number_residual_linf;
    double cell_number_tolerance;
    double rank_number_residual_linf;
    double rank_number_tolerance;
    double global_number_residual;
    double global_number_tolerance;
    double rho_before_after_relative_linf;
    int rng_advanced_exactly_once;
    int ledger_advanced_exactly_once;
    int operators_ok;
    int accepted_integrator_state_unchanged;
    int conversion_active_nonzero;
    int return_hysteretic_nonzero;
    int mpi_case_count;
};

GateMetrics run_all_cases(int rank, int size)
{
    GateMetrics m;
    m.case_count = 0;
    m.cell_number_residual_linf = 0.0;
    m.cell_number_tolerance = 0.0;
    m.rank_number_residual_linf = 0.0;
    m.rank_number_tolerance = 0.0;
    m.global_number_residual = 0.0;
    m.global_number_tolerance = 0.0;
    m.rho_before_after_relative_linf = 0.0;
    m.rng_advanced_exactly_once = 1;
    m.ledger_advanced_exactly_once = 1;
    m.operators_ok = 1;
    m.accepted_integrator_state_unchanged = 1;
    m.conversion_active_nonzero = 0;
    m.return_hysteretic_nonzero = 0;
    m.mpi_case_count = 0;

    const int case_count_limit = (size > 1) ? 7 : 7;
    for (int case_id = 0; case_id < case_count_limit; ++case_id) {
        const CaseOutcome out = run_case(rank, size, case_id);
        ++m.case_count;
        if (size > 1) ++m.mpi_case_count;
        if (rank == 0) {
            std::fprintf(stderr,
                "[jc0-case] rank=%d %s c2_ok=%d conversion_ok=%d return_ok=%d "
                "conversion_calls=%d return_calls=%d stochastic=%d "
                "cell_res=%.17g export=%.17g state_unchanged=%d case_pass=%d\n",
                rank, out.name.c_str(), out.c2_ok ? 1 : 0,
                out.conversion_ok ? 1 : 0, out.return_ok ? 1 : 0,
                out.conversion_calls, out.return_calls,
                out.stochastic_ran ? 1 : 0, out.cell_residual_linf,
                out.export_number, out.integrator_state_unchanged ? 1 : 0,
                out.case_pass ? 1 : 0);
        }
        m.cell_number_residual_linf = std::max(
            m.cell_number_residual_linf, out.cell_residual_linf);
        m.cell_number_tolerance = std::max(
            m.cell_number_tolerance, out.cell_tolerance);
        m.rank_number_residual_linf = std::max(
            m.rank_number_residual_linf, out.rank_residual_linf);
        m.rank_number_tolerance = std::max(
            m.rank_number_tolerance, out.rank_tolerance);
        m.global_number_residual = std::max(
            m.global_number_residual, out.global_residual);
        m.global_number_tolerance = std::max(
            m.global_number_tolerance, out.global_tolerance);
        m.rho_before_after_relative_linf = std::max(
            m.rho_before_after_relative_linf, out.rho_rel_linf);

        if (out.collision_half_calls != 1 ||
            out.stochastic_ran != out.expected_stochastic)
            m.rng_advanced_exactly_once = 0;
        if (out.conversion_calls > 1 ||
            out.return_calls != out.expected_return_calls)
            m.ledger_advanced_exactly_once = 0;
        if (!out.operators_ok) m.operators_ok = 0;
        if (!out.integrator_state_unchanged) m.accepted_integrator_state_unchanged = 0;
        if (out.export_number > 0.0) m.conversion_active_nonzero = 1;
        if (out.return_calls > 0) m.return_hysteretic_nonzero = 1;
    }
    return m;
}

bool write_result_file(const std::string& path, const GateMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "case_count=" << m.case_count << "\n";
    out << "cell_number_residual_linf=" << m.cell_number_residual_linf << "\n";
    out << "cell_number_tolerance=" << m.cell_number_tolerance << "\n";
    out << "rank_number_residual_linf=" << m.rank_number_residual_linf << "\n";
    out << "rank_number_tolerance=" << m.rank_number_tolerance << "\n";
    out << "global_number_residual=" << m.global_number_residual << "\n";
    out << "global_number_tolerance=" << m.global_number_tolerance << "\n";
    out << "rho_before_after_relative_linf="
        << m.rho_before_after_relative_linf << "\n";
    out << "rng_advanced_exactly_once=" << m.rng_advanced_exactly_once << "\n";
    out << "ledger_advanced_exactly_once=" << m.ledger_advanced_exactly_once
        << "\n";
    out << "accepted_integrator_state_unchanged="
        << m.accepted_integrator_state_unchanged << "\n";
    out << "conversion_active_nonzero=" << m.conversion_active_nonzero << "\n";
    out << "return_hysteretic_nonzero=" << m.return_hysteretic_nonzero << "\n";
    out << "mpi_case_count=" << m.mpi_case_count << "\n";
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
        std::cerr << "usage: vpfp_post_field_charge_invariance_test "
                     "--case all|mpi [--result <path>]\n";
    }

    GateMetrics m;
    if (ok) m = run_all_cases(rank, size);

    // The post-field charge-invariance gate (section 3.3): every cell, rank
    // and global residual stays inside the JC0 number tolerance
    // T_N = 4096 eps max(1, N_before, N_after), and the relative charge
    // density linf stays below 1e-12.
    const bool local_pass = ok &&
        m.cell_number_residual_linf <= m.cell_number_tolerance &&
        m.rank_number_residual_linf <= m.rank_number_tolerance &&
        m.global_number_residual <= m.global_number_tolerance &&
        m.rho_before_after_relative_linf <= 1.0e-12 &&
        m.rng_advanced_exactly_once == 1 &&
        m.ledger_advanced_exactly_once == 1 &&
        m.operators_ok == 1 &&
        m.accepted_integrator_state_unchanged == 1 &&
        m.conversion_active_nonzero == 1 &&
        m.return_hysteretic_nonzero == 1 &&
        m.case_count >= 7;
    int local_pass_int = local_pass ? 1 : 0;
    int global_pass_int = 0;
    MPI_Allreduce(&local_pass_int, &global_pass_int, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    const bool pass = global_pass_int != 0;

    if (rank == 0) {
        if (!write_result_file(args.result_path, m, pass)) {
            std::cerr << "cannot write result file\n";
        }
        std::cout << std::setprecision(17)
                  << "case_count=" << m.case_count
                  << " cell_number_residual_linf="
                  << m.cell_number_residual_linf
                  << " cell_number_tolerance=" << m.cell_number_tolerance
                  << " rank_number_residual_linf="
                  << m.rank_number_residual_linf
                  << " rank_number_tolerance=" << m.rank_number_tolerance
                  << " global_number_residual=" << m.global_number_residual
                  << " global_number_tolerance=" << m.global_number_tolerance
                  << " rho_before_after_relative_linf="
                  << m.rho_before_after_relative_linf
                  << " rng_advanced_exactly_once="
                  << m.rng_advanced_exactly_once
                  << " ledger_advanced_exactly_once="
                  << m.ledger_advanced_exactly_once
                  << " accepted_integrator_state_unchanged="
                  << m.accepted_integrator_state_unchanged
                  << " conversion_active_nonzero="
                  << m.conversion_active_nonzero
                  << " return_hysteretic_nonzero="
                  << m.return_hysteretic_nonzero
                  << " mpi_case_count=" << m.mpi_case_count << "\n";
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
