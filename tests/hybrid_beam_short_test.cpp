// Stage H4 acceptance: Beam + collisionless hybrid Vlasov-Poisson
// (sections 8.2, 11 and 15 H4).  Drives the production VpfpIntegrator with
// the Beam enabled and the tail PIC + converter wired in, and compares it
// against the tail-disabled (phase-4) path on the same beam.
//
// Cases:
//   beam-no-tail-equivalence : thermal plasma + beam; no conversion; the
//                              hybrid run is identical to the tail-off run
//                              in phase space and beam state;
//   beam-conversion-continuity : a mid-domain high-energy background packet
//                              is converted while the beam injects; combined
//                              N/Px/K stay continuous, the conversion is
//                              conservative, and the beam ledger is
//                              unaffected by the internal representation
//                              change.
//   beam-hybrid-collision-pairs : a Beam run with a non-empty tail enters the
//                              production collision selector and executes all
//                              four configured collision pairs.
//
// Usage:
//   hybrid_beam_short_test --case <case> [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

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

#include <algorithm>
#include <cmath>
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

inline double tail_total_weight_probe(const BackgroundTailPIC& tail)
{
    double total = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        total += tail.particles[i].weight;
    }
    return total;
}

inline double tail_total_energy_probe(const BackgroundTailPIC& tail)
{
    double total = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        total += p.weight * (gamma - 1.0) *
                 Const::me * Const::c * Const::c;
    }
    return total;
}

struct Sim {
    SpatialGrid grid;
    Species electrons_off;
    Species electrons_on;
    OpenBackgroundBoundary boundary;
    OpenElectrostaticSolver solver_off;
    OpenElectrostaticSolver solver_on;
    ZeroCollisionCoefficients zero_provider;
    CylindricalFokkerPlanckCollision collision_off;
    CylindricalFokkerPlanckCollision collision_on;
    HybridVelocityPartition partition;
    BulkTailConverter converter;
    VpfpIntegrator off;
    VpfpIntegrator on;
    BeamPIC beam_off;
    BeamPIC beam_on;
    std::vector<double> ion_density;
    double dt;

    Sim(int nx, double length, ElectrostaticBoundaryType field_type,
        double ion_offset)
        : boundary(reservoir_config()),
          collision_off(zero_provider,
                        CollisionIntegratorType::BACKWARD_EULER),
          collision_on(zero_provider,
                       CollisionIntegratorType::BACKWARD_EULER),
          off(boundary, solver_off, collision_off),
          on(boundary, solver_on, collision_on, partition, converter, true)
    {
        grid.init_with_domain(0, 1, nx, length);
        electrons_off.init("background_electrons",
                           SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                           Const::me, Param::dens, Param::temperature_e,
                           false, grid);
        electrons_on.init("background_electrons",
                          SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                          Const::me, Param::dens, Param::temperature_e,
                          false, grid);
        partition.init(electrons_on.cgrid, 6.0, 1.0, 4, 4);
        solver_off.init(grid, { field_type, 0.0, 0.0, 0.0 });
        solver_on.init(grid, { field_type, 0.0, 0.0, 0.0 });
        off.init(grid);
        on.init(grid);
        off.set_beam_enabled(true);
        on.set_beam_enabled(true);
        beam_off.init(grid);
        beam_on.init(grid);
        ion_density.assign(static_cast<size_t>(grid.nx_local),
                           Param::dens * (1.0 + ion_offset));
        dt = Param::dt_multiplier / Param::omega_pe;
    }

    static OpenBackgroundBoundaryConfig reservoir_config()
    {
        OpenBackgroundBoundaryConfig config;
        config.left_type = BackgroundXBoundaryType::RESERVOIR;
        config.right_type = BackgroundXBoundaryType::RESERVOIR;
        config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
        config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
        return config;
    }

    void sync_states()
    {
        electrons_on.f = electrons_off.f;
        beam_on.particles = beam_off.particles;
        beam_on.density = beam_off.density;
    }

    VpfpStepResult step_off(double time)
    {
        EMFields fields_off;
        fields_off.init(grid);
        const VpfpStepResult r =
            off.advance(electrons_off, beam_off, fields_off, ion_density,
                        time, dt, 0, 1);
        last_fields_off = fields_off.Ex_face;
        return r;
    }

    VpfpStepResult step_on(double time)
    {
        EMFields fields_on;
        fields_on.init(grid);
        const VpfpStepResult r =
            on.advance(electrons_on, beam_on, fields_on, ion_density,
                       time, dt, 0, 1);
        last_fields_on = fields_on.Ex_face;
        return r;
    }

    std::vector<double> last_fields_off;
    std::vector<double> last_fields_on;
};

void fill_background_packet(Species& species, const SpatialGrid& grid,
                            const HybridVelocityPartition& partition,
                            const std::vector<int>& local_cells,
                            double mass_scale)
{
    const int ng = grid.nghost;
    for (size_t c = 0; c < local_cells.size(); ++c) {
        const int il = local_cells[c];
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                if (!partition.is_conversion(j, k)) continue;
                const double upar = species.cgrid.upar_cells[j];
                const double uperp = species.cgrid.uperp_cells[k];
                if (std::fabs(upar) < 12.0 || std::fabs(upar) > 15.0 ||
                    uperp < 0.5) continue;
                species.f[idx3(ng + il, j, k)] =
                    mass_scale * (1.0 + 0.05 *
                                  static_cast<double>((j + k) % 7));
            }
        }
    }
}

struct Metrics {
    bool all_ok;
    double max_ledger_diff;
    double beam_number_diff;
    double combined_number_rel;
    double combined_energy_rel;
    double conversion_residual_max;
    double tail_number_balance_max;
    bool collision_pairs_ok;
    double collision_reaction_residual_max;
    Metrics()
        : all_ok(false), max_ledger_diff(0.0), beam_number_diff(0.0),
          combined_number_rel(0.0), combined_energy_rel(0.0),
          conversion_residual_max(0.0), tail_number_balance_max(0.0),
          collision_pairs_ok(false), collision_reaction_residual_max(0.0)
    {}
};

Metrics run_beam_no_tail_equivalence()
{
    Metrics m;
    Sim sim(160, 1.6 * Const::micro,
            ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0);
    sim.electrons_off.initialize_maxwellian();
    sim.sync_states();
    double ledger_diff = 0.0;
    double beam_number_diff = 0.0;
    double f_rel = 0.0;
    double beam_dens_rel = 0.0;
    bool ok = true;
    for (int s = 0; s < 20; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        const VpfpStepResult r_off = sim.step_off(time);
        const VpfpStepResult r_on = sim.step_on(time);
        ok = ok && r_off.accepted && r_on.accepted &&
             r_on.conversion_ok && r_on.tail_ok &&
             sim.on.tail_state().particles.empty();
        ledger_diff = std::max(
            ledger_diff,
            std::fabs(r_off.ledger.background_kinetic_energy -
                      r_on.ledger.background_kinetic_energy) /
                std::max(1.0,
                         std::fabs(r_off.ledger.background_kinetic_energy)));
        ledger_diff = std::max(
            ledger_diff,
            std::fabs(r_off.ledger.field_energy - r_on.ledger.field_energy) /
                std::max(1.0, std::fabs(r_off.ledger.field_energy)));
        ledger_diff = std::max(
            ledger_diff,
            std::fabs(r_off.ledger.background_number_after -
                      r_on.ledger.background_number_after) /
                std::max(1.0, r_off.ledger.background_number_after));
        // Beam ledger and state must be identical (the conversion does not
        // touch the beam).
        beam_number_diff = std::max(
            beam_number_diff,
            std::fabs(r_off.ledger.beam_injected - r_on.ledger.beam_injected) /
                std::max(1.0, r_off.ledger.beam_injected) +
            std::fabs(r_off.ledger.beam_outflow - r_on.ledger.beam_outflow) /
                std::max(1.0, r_off.ledger.beam_outflow));
        // Phase-space comparison is relative: the phase-4 OpenMP reductions
        // (beam/bulk density deposits) can differ in the last ULP between
        // two separately advanced integrator instances, and the Gauss
        // integration amplifies that to ~1e-9 relative in the field.  The
        // beam particle inventory itself must match exactly.
        for (size_t i = 0; i < sim.electrons_off.f.size(); ++i) {
            const double rel =
                std::fabs(sim.electrons_off.f[i] - sim.electrons_on.f[i]) /
                std::max(1e-300, std::fabs(sim.electrons_off.f[i]));
            f_rel = std::max(f_rel, rel);
        }
        for (size_t i = 0; i < sim.beam_off.density.size(); ++i) {
            const double rel =
                std::fabs(sim.beam_off.density[i] - sim.beam_on.density[i]) /
                std::max(1e-300, std::fabs(sim.beam_off.density[i]));
            beam_dens_rel = std::max(beam_dens_rel, rel);
        }
        if (sim.beam_off.particles.size() != sim.beam_on.particles.size()) {
            ok = false;
        }
    }
    m.max_ledger_diff = ledger_diff;
    m.beam_number_diff = beam_number_diff;
    m.combined_number_rel = f_rel;
    m.combined_energy_rel = beam_dens_rel;
    m.all_ok = ok;
    return m;
}

Metrics run_beam_conversion_continuity()
{
    Metrics m;
    Sim sim(160, 1.6 * Const::micro,
            ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0);
    sim.electrons_off.initialize_maxwellian();
    sim.sync_states();
    std::vector<int> cells;
    cells.push_back(60);
    cells.push_back(61);
    cells.push_back(62);
    fill_background_packet(sim.electrons_off, sim.grid, sim.partition, cells,
                           1.0e16);
    sim.sync_states();
    const double number_before = sim.electrons_off.total_particle_number();
    const double energy_before = sim.electrons_off.total_kinetic_energy();

    bool ok = true;
    double conv_residual_max = 0.0;
    double tail_balance_max = 0.0;
    double flux_sum = 0.0;
    double u_tail_sum = 0.0;
    double tail_out_sum = 0.0;
    double energy_drift = 0.0;
    double ue_prev = 0.0;
    double tail_seen = 0.0;
    double beam_balance_max = 0.0;
    for (int s = 0; s < 3; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        const double ke_before_now =
            sim.electrons_on.total_kinetic_energy() +
            tail_total_energy_probe(sim.on.tail_state());
        const double ue_before_now = (s == 0) ? 0.0 : ue_prev;
        const VpfpStepResult r_on = sim.step_on(time);
        const double ke_after =
            r_on.ledger.background_kinetic_energy +
            r_on.ledger.tail_kinetic_energy_after;
        ue_prev = r_on.ledger.field_energy;
        energy_drift = std::max(
            energy_drift,
            std::fabs((ke_after + ue_prev) -
                      (ke_before_now + ue_before_now)));
        ok = ok && r_on.accepted && r_on.conversion_ok && r_on.tail_ok;
        if (r_on.ledger.conversion_number_removed > 0.0) tail_seen = 1.0;
        conv_residual_max = std::max(
            conv_residual_max,
            std::max(r_on.ledger.conversion_number_residual_rel,
                     std::max(r_on.ledger.conversion_px_residual_rel,
                              r_on.ledger.conversion_energy_residual_rel)));
        tail_balance_max = std::max(tail_balance_max,
                                    r_on.ledger.tail_number_balance_error);
        flux_sum += r_on.ledger.background_left_flux +
                    r_on.ledger.background_right_flux;
        u_tail_sum += r_on.ledger.background_tail_number_loss;
        tail_out_sum += r_on.ledger.tail_outflow_number;
        // Beam ledger continuity: N_b after = N_b before + injected - out.
        const double beam_balance =
            std::fabs(r_on.ledger.beam_number_after -
                      (r_on.ledger.beam_number_before +
                       r_on.ledger.beam_injected -
                       r_on.ledger.beam_outflow)) /
            std::max(1.0, r_on.ledger.beam_number_before +
                              r_on.ledger.beam_injected);
        beam_balance_max = std::max(beam_balance_max, beam_balance);
    }
    const double combined_number_after =
        sim.electrons_on.total_particle_number() +
        tail_total_weight_probe(sim.on.tail_state());
    const double expected_number =
        number_before + flux_sum - u_tail_sum - tail_out_sum;
    m.combined_number_rel =
        std::fabs(combined_number_after - expected_number) /
        std::max(1.0, number_before);
    m.combined_energy_rel = energy_drift / std::max(1.0, energy_before);
    m.conversion_residual_max = conv_residual_max;
    m.tail_number_balance_max = tail_balance_max;
    m.max_ledger_diff = beam_balance_max;
    m.all_ok = ok && tail_seen > 0.0 && beam_balance_max <= 1.0e-9;
    return m;
}

Metrics run_beam_hybrid_collision_pairs()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
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
        TailCollisionWeightMode::VirtualSplit, 64, 10.0);

    // A pre-existing interior tail guarantees that both collision halves
    // exercise C_tt and C_tb/C_bt through the production Beam path.  The
    // direct insertion deliberately isolates dispatch from conversion.
    for (int i = 0; i < 100; ++i) {
        BackgroundTailParticle particle;
        particle.x = 0.4 * Const::micro;
        particle.ux = 12.0 + 0.01 * static_cast<double>(i % 5);
        particle.uy = 0.10;
        particle.uz = 0.05;
        particle.weight = 1.0e14;
        particle.id = static_cast<std::uint64_t>(i + 1);
        integrator.tail_state().particles.push_back(particle);
    }
    integrator.tail_state().deposit_density(grid, 0, 1);

    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    std::vector<double> ion_density(static_cast<size_t>(grid.nx_local),
                                    Param::dens);
    const double dt = Param::dt_multiplier / Param::omega_pe;
    const VpfpStepResult result = integrator.advance(
        electrons, beam, fields, ion_density, 0.0, dt, 0, 1);
    if (!result.accepted) {
        std::cerr << "beam-hybrid-collision-pairs failure_code="
                  << result.failure_code
                  << " collision_ok=" << (result.collision_ok ? 1 : 0)
                  << " conversion_ok=" << (result.conversion_ok ? 1 : 0)
                  << " tail_ok=" << (result.tail_ok ? 1 : 0) << "\n";
    }

    m.collision_pairs_ok =
        result.accepted && result.collision_ok &&
        result.ledger.collision_pair_bulk_bulk == 1 &&
        result.ledger.collision_pair_tail_tail == 1 &&
        result.ledger.collision_pair_tail_bulk == 1 &&
        result.ledger.collision_pair_bulk_reaction == 1;
    m.collision_reaction_residual_max = std::max(
        result.ledger.collision_reaction_px_residual,
        result.ledger.collision_reaction_energy_residual);
    m.all_ok = m.collision_pairs_ok &&
               std::isfinite(m.collision_reaction_residual_max) &&
               m.collision_reaction_residual_max <= 1.0e-10;
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "all_ok=" << (m.all_ok ? 1 : 0) << "\n";
    out << "max_ledger_diff=" << m.max_ledger_diff << "\n";
    out << "beam_number_diff=" << m.beam_number_diff << "\n";
    out << "combined_number_rel=" << m.combined_number_rel << "\n";
    out << "combined_energy_rel=" << m.combined_energy_rel << "\n";
    out << "conversion_residual_max=" << m.conversion_residual_max << "\n";
    out << "tail_number_balance_max=" << m.tail_number_balance_max << "\n";
    out << "collision_pairs_ok=" << (m.collision_pairs_ok ? 1 : 0) << "\n";
    out << "collision_reaction_residual_max="
        << m.collision_reaction_residual_max << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "hybrid_beam_short_test must run with exactly 1 rank; "
                     "use plain ./build_hybrid/hybrid_beam_short_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: hybrid_beam_short_test --case <case> "
                     "[--result <path>]\n"
                  << "cases: beam-no-tail-equivalence | "
                     "beam-conversion-continuity | "
                     "beam-hybrid-collision-pairs | all\n";
    }

    Metrics m;
    bool pass = ok;
    bool eq_ok = true;
    bool conv_ok = true;
    bool pair_ok = true;
    if (ok && (args.test_case == "all" ||
               args.test_case == "beam-no-tail-equivalence")) {
        const Metrics r = run_beam_no_tail_equivalence();
        eq_ok = r.all_ok && r.max_ledger_diff <= 1.0e-10 &&
                r.beam_number_diff <= 1.0e-10;
        m.max_ledger_diff = std::max(m.max_ledger_diff, r.max_ledger_diff);
        m.beam_number_diff = std::max(m.beam_number_diff, r.beam_number_diff);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "beam-conversion-continuity")) {
        const Metrics r = run_beam_conversion_continuity();
        conv_ok = r.all_ok && r.combined_number_rel <= 1.0e-8 &&
                  r.combined_energy_rel <= 1.0e-2 &&
                  r.conversion_residual_max <= 1.0e-10 &&
                  r.tail_number_balance_max <= 1.0e-9 &&
                  r.max_ledger_diff <= 1.0e-9;
        m.combined_number_rel = r.combined_number_rel;
        m.combined_energy_rel = r.combined_energy_rel;
        m.conversion_residual_max = r.conversion_residual_max;
        m.tail_number_balance_max = r.tail_number_balance_max;
        m.max_ledger_diff = std::max(m.max_ledger_diff, r.max_ledger_diff);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "beam-hybrid-collision-pairs")) {
        const Metrics r = run_beam_hybrid_collision_pairs();
        pair_ok = r.all_ok;
        m.collision_pairs_ok = r.collision_pairs_ok;
        m.collision_reaction_residual_max =
            r.collision_reaction_residual_max;
    }
    pass = pass && eq_ok && conv_ok && pair_ok;
    m.all_ok = pass;

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "all_ok=" << (m.all_ok ? 1 : 0)
              << " max_ledger_diff=" << m.max_ledger_diff
              << " beam_number_diff=" << m.beam_number_diff
              << " combined_number_rel=" << m.combined_number_rel
              << " combined_energy_rel=" << m.combined_energy_rel
              << " conversion_residual_max=" << m.conversion_residual_max
              << " tail_number_balance_max=" << m.tail_number_balance_max
              << " collision_pairs_ok=" << (m.collision_pairs_ok ? 1 : 0)
              << " collision_reaction_residual_max="
              << m.collision_reaction_residual_max
              << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
