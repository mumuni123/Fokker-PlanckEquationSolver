// Stage H3 acceptance: no-beam hybrid Vlasov-Poisson integration
// (sections 8.2, 9, 11 and 15 H3).  Drives the production VpfpIntegrator
// with the tail PIC + converter wired in (background-tail-mode pic) and
// compares it against the tail-disabled (phase-4) path on the same state.
//
// Cases:
//   maxwellian-no-tail            : a thermal Maxwellian produces no tail
//                                   and the hybrid step is bitwise equal to
//                                   the tail-off (phase-4) step;
//   packet-continuity             : a high-energy packet above the threshold
//                                   is converted; combined N/Px/K stay
//                                   continuous and match the Eulerian;
//   straddling-threshold          : a packet spanning the threshold (some
//                                   cells converted, some remaining Eulerian)
//                                   keeps the combined moments continuous;
//   langmuir-small-amplitude      : small-amplitude wave, no conversion;
//                                   hybrid == tail-off and U_E oscillates;
//   left-e-reservoir              : LEFT_E + reservoir equilibrium with no
//                                   beam; the field energy stays near zero
//                                   and the background is conserved.
//
// Usage:
//   hybrid_no_beam_test --case <case> [--result <path>]
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
#include <iomanip>
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
        off.set_beam_enabled(false);
        on.set_beam_enabled(false);
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
    }

    VpfpStepResult step_off(double time)
    {
        BeamPIC beam_off;
        beam_off.init(grid);
        EMFields fields_off;
        fields_off.init(grid);
        return off.advance(electrons_off, beam_off, fields_off, ion_density,
                           time, dt, 0, 1);
    }

    VpfpStepResult step_on(double time)
    {
        BeamPIC beam_on;
        beam_on.init(grid);
        EMFields fields_on;
        fields_on.init(grid);
        return on.advance(electrons_on, beam_on, fields_on, ion_density,
                          time, dt, 0, 1);
    }
};

// Fill a deterministic high-energy packet in the 2<|u_parallel|<=15,
// u_perp>0.5 cells of the given local cells.  only_above keeps only
// conversion cells, only_below only non-conversion cells, neither keeps all
// (which spans the threshold).  The outermost u_parallel cells are left
// empty on purpose: the Eulerian u-remap loses mass at the open velocity
// boundary, so a physical test packet must not start on the edge cells.
void fill_packet(Species& species, const SpatialGrid& grid,
                 const HybridVelocityPartition& partition,
                 const std::vector<int>& local_cells, double mass_scale,
                 bool only_above, bool only_below,
                 double upar_lo, double upar_hi)
{
    const int ng = grid.nghost;
    for (size_t c = 0; c < local_cells.size(); ++c) {
        const int il = local_cells[c];
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const bool conv = partition.is_conversion(j, k);
                if (only_above && !conv) continue;
                if (only_below && conv) continue;
                const double upar = species.cgrid.upar_cells[j];
                const double uperp = species.cgrid.uperp_cells[k];
                if (std::fabs(upar) < upar_lo || std::fabs(upar) > upar_hi ||
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
    double combined_number_rel;
    double combined_energy_rel;
    double conversion_residual_max;
    double tail_number_balance_max;
    double field_energy_rel_change;
    double left_e_max_field_energy;
    double left_e_number_rel;
    Metrics()
        : all_ok(false), max_ledger_diff(0.0), combined_number_rel(0.0),
          combined_energy_rel(0.0), conversion_residual_max(0.0),
          tail_number_balance_max(0.0), field_energy_rel_change(0.0),
          left_e_max_field_energy(0.0), left_e_number_rel(0.0)
    {}
};

Metrics run_maxwellian_no_tail()
{
    Metrics m;
    Sim sim(160, 1.6 * Const::micro,
            ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0);
    sim.electrons_off.initialize_maxwellian();
    sim.sync_states();
    // The phase-4 OpenMP reductions (total number / kinetic energy) can
    // differ in the last ULP between two separately advanced integrator
    // instances, so the ledger comparison is relative (~1e-16), while the
    // phase-space state comparison below stays bitwise.
    double ledger_diff = 0.0;
    bool ok = true;
    for (int s = 0; s < 5; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        const VpfpStepResult r_off = sim.step_off(time);
        const VpfpStepResult r_on = sim.step_on(time);
        ok = ok && r_off.accepted && r_on.accepted &&
             r_on.conversion_ok && r_on.tail_ok &&
             sim.on.tail_state().particles.empty() &&
             r_on.ledger.conversion_number_removed == 0.0;
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
        if (sim.electrons_off.f != sim.electrons_on.f) ok = false;
    }
    m.max_ledger_diff = ledger_diff;
    m.all_ok = ok;
    return m;
}

Metrics run_packet_continuity()
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
    fill_packet(sim.electrons_off, sim.grid, sim.partition, cells,
                1.0e16, true, false, 2.0, 15.0);
    sim.sync_states();

    const double number_before = sim.electrons_off.total_particle_number();
    const double energy_before = sim.electrons_off.total_kinetic_energy();
    bool ok = true;
    double conv_residual_max = 0.0;
    double tail_balance_max = 0.0;
    double remap_residual_max = 0.0;
    double flux_sum = 0.0;
    double u_tail_sum = 0.0;
    double tail_out_sum = 0.0;
    double energy_drift = 0.0;
    double ue_prev = 0.0;
    double tail_seen = 0.0;
    for (int s = 0; s < 2; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        // Total energy continuity: the packet's space charge exchanges KE
        // with the field energy, so track K_e + U_E across the step.  The
        // initial field is not solved in the test, so U_E starts at zero.
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
        if (r_on.ledger.conversion_number_removed > 0.0) tail_seen = 1.0;
        ok = ok && r_on.accepted && r_on.conversion_ok && r_on.tail_ok;
        conv_residual_max = std::max(
            conv_residual_max,
            std::max(r_on.ledger.conversion_number_residual_rel,
                     std::max(r_on.ledger.conversion_px_residual_rel,
                              r_on.ledger.conversion_energy_residual_rel)));
        tail_balance_max = std::max(tail_balance_max,
                                    r_on.ledger.tail_number_balance_error);
        remap_residual_max = std::max(remap_residual_max,
                                      r_on.ledger.remap_ledger_residual);
        flux_sum += r_on.ledger.background_left_flux +
                    r_on.ledger.background_right_flux;
        u_tail_sum += r_on.ledger.background_tail_number_loss;
        tail_out_sum += r_on.ledger.tail_outflow_number;
    }
    // Global combined N/K continuity of the hybrid run itself: the
    // conversion is an internal representation change, so bulk+tail must
    // reproduce the initial totals (the packet stays inside the spatial
    // domain and the velocity domain over this short run).
    const double combined_number_after =
        sim.electrons_on.total_particle_number() +
        tail_total_weight_probe(sim.on.tail_state());
    // Combined continuity including the reservoir exchange, the bulk
    // velocity-tail loss and the tail open outflows (section 11.1).
    const double expected_number =
        number_before + flux_sum - u_tail_sum - tail_out_sum;
    m.combined_number_rel =
        std::fabs(combined_number_after - expected_number) /
        std::max(1.0, number_before);
    m.combined_energy_rel =
        energy_drift / std::max(1.0, energy_before);
    m.conversion_residual_max = conv_residual_max;
    m.tail_number_balance_max = tail_balance_max;
    m.max_ledger_diff = remap_residual_max;
    m.all_ok = ok && tail_seen > 0.0;
    return m;
}

Metrics run_straddling_threshold()
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
    cells.push_back(63);
    cells.push_back(64);
    // All high-|u_parallel| cells: some above, some below the threshold.
    fill_packet(sim.electrons_off, sim.grid, sim.partition, cells,
                1.0e16, false, false, 12.0, 15.0);
    sim.sync_states();
    const double number_before = sim.electrons_off.total_particle_number();
    const double energy_before = sim.electrons_off.total_kinetic_energy();

    bool ok = true;
    double conv_residual_max = 0.0;
    double remap_residual_max = 0.0;
    double flux_sum = 0.0;
    double u_tail_sum = 0.0;
    double tail_out_sum = 0.0;
    double energy_drift = 0.0;
    double ue_prev = 0.0;
    double tail_seen = 0.0;
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
        remap_residual_max = std::max(remap_residual_max,
                                      r_on.ledger.remap_ledger_residual);
        flux_sum += r_on.ledger.background_left_flux +
                    r_on.ledger.background_right_flux;
        u_tail_sum += r_on.ledger.background_tail_number_loss;
        tail_out_sum += r_on.ledger.tail_outflow_number;
    }
    const double combined_number_after =
        sim.electrons_on.total_particle_number() +
        tail_total_weight_probe(sim.on.tail_state());
    const double expected_number =
        number_before + flux_sum - u_tail_sum - tail_out_sum;
    m.combined_number_rel =
        std::fabs(combined_number_after - expected_number) /
        std::max(1.0, number_before);
    m.combined_energy_rel =
        energy_drift / std::max(1.0, energy_before);
    m.conversion_residual_max = conv_residual_max;
    m.max_ledger_diff = remap_residual_max;
    m.all_ok = ok && tail_seen > 0.0;
    return m;
}

Metrics run_langmuir_small_amplitude()
{
    Metrics m;
    Sim sim(160, 1.6 * Const::micro,
            ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0);
    sim.electrons_off.initialize_maxwellian();
    const int ng = sim.grid.nghost;
    for (int il = 0; il < sim.grid.nx_local; ++il) {
        const double amp = 1.0e-3 * std::sin(
            2.0 * Const::pi * static_cast<double>(il) / 16.0);
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                sim.electrons_off.f[idx3(ng + il, j, k)] *= (1.0 + amp);
            }
        }
    }
    sim.sync_states();
    const int steps = 16;
    bool ok = true;
    double ledger_diff = 0.0;
    double ue0 = 0.0;
    double ue1 = 0.0;
    for (int s = 0; s < steps; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        const VpfpStepResult r_off = sim.step_off(time);
        const VpfpStepResult r_on = sim.step_on(time);
        ok = ok && r_off.accepted && r_on.accepted &&
             r_on.conversion_ok && r_on.tail_ok &&
             sim.on.tail_state().particles.empty();
        ledger_diff = std::max(
            ledger_diff,
            std::fabs(r_off.ledger.field_energy - r_on.ledger.field_energy) /
                std::max(1.0, std::fabs(r_off.ledger.field_energy)));
        if (s == 1) ue0 = r_on.ledger.field_energy;
        if (s == steps - 1) ue1 = r_on.ledger.field_energy;
    }
    m.field_energy_rel_change = std::fabs(ue1 - ue0) /
                                std::max(1.0, std::fabs(ue0) + 1.0);
    m.max_ledger_diff = ledger_diff;
    m.all_ok = ok && sim.electrons_off.f == sim.electrons_on.f;
    return m;
}

Metrics run_left_e_reservoir()
{
    Metrics m;
    Sim sim(160, 1.6 * Const::micro,
            ElectrostaticBoundaryType::LEFT_E, 0.0);
    sim.electrons_off.initialize_maxwellian();
    sim.sync_states();
    bool ok = true;
    double max_field_energy = 0.0;
    double number_rel = 0.0;
    for (int s = 0; s < 5; ++s) {
        const double time = static_cast<double>(s) * sim.dt;
        const VpfpStepResult r_off = sim.step_off(time);
        const VpfpStepResult r_on = sim.step_on(time);
        ok = ok && r_off.accepted && r_on.accepted &&
             r_on.conversion_ok && r_on.tail_ok &&
             sim.on.tail_state().particles.empty();
        max_field_energy = std::max(max_field_energy,
                                    r_on.ledger.field_energy);
        number_rel = std::max(
            number_rel,
            std::fabs(r_on.ledger.combined_number_after -
                      r_on.ledger.combined_number_before) /
                std::max(1.0, r_on.ledger.combined_number_before));
    }
    m.left_e_max_field_energy = max_field_energy;
    m.left_e_number_rel = number_rel;
    m.all_ok = ok;
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
    out << "combined_number_rel=" << m.combined_number_rel << "\n";
    out << "combined_energy_rel=" << m.combined_energy_rel << "\n";
    out << "conversion_residual_max=" << m.conversion_residual_max << "\n";
    out << "tail_number_balance_max=" << m.tail_number_balance_max << "\n";
    out << "field_energy_rel_change=" << m.field_energy_rel_change << "\n";
    out << "left_e_max_field_energy=" << m.left_e_max_field_energy << "\n";
    out << "left_e_number_rel=" << m.left_e_number_rel << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "hybrid_no_beam_test must run with exactly 1 rank; "
                     "use plain ./build_hybrid/hybrid_no_beam_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: hybrid_no_beam_test --case <case> "
                     "[--result <path>]\n"
                  << "cases: maxwellian-no-tail | packet-continuity | "
                     "straddling-threshold | langmuir-small-amplitude | "
                     "left-e-reservoir | all\n";
    }

    Metrics m;
    bool pass = ok;
    bool maxwell_ok = true;
    bool packet_ok = true;
    bool straddle_ok = true;
    bool langmuir_ok = true;
    bool left_e_ok = true;
    if (ok && (args.test_case == "all" ||
               args.test_case == "maxwellian-no-tail")) {
        const Metrics r = run_maxwellian_no_tail();
        maxwell_ok = r.all_ok && r.max_ledger_diff <= 1.0e-10;
        m.max_ledger_diff = std::max(m.max_ledger_diff, r.max_ledger_diff);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "packet-continuity")) {
        const Metrics r = run_packet_continuity();
        packet_ok = r.all_ok && r.combined_number_rel <= 1.0e-9 &&
                    r.combined_energy_rel <= 1.0e-2 &&
                    r.conversion_residual_max <= 1.0e-10 &&
                    r.tail_number_balance_max <= 1.0e-9;
        m.combined_number_rel = std::max(m.combined_number_rel,
                                         r.combined_number_rel);
        m.combined_energy_rel = std::max(m.combined_energy_rel,
                                         r.combined_energy_rel);
        m.conversion_residual_max = std::max(
            m.conversion_residual_max, r.conversion_residual_max);
        m.tail_number_balance_max = std::max(
            m.tail_number_balance_max, r.tail_number_balance_max);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "straddling-threshold")) {
        const Metrics r = run_straddling_threshold();
        straddle_ok = r.all_ok && r.combined_number_rel <= 1.0e-9 &&
                      r.combined_energy_rel <= 1.0e-2 &&
                      r.conversion_residual_max <= 1.0e-10;
        m.combined_number_rel = std::max(m.combined_number_rel,
                                         r.combined_number_rel);
        m.combined_energy_rel = std::max(m.combined_energy_rel,
                                         r.combined_energy_rel);
        m.conversion_residual_max = std::max(
            m.conversion_residual_max, r.conversion_residual_max);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "langmuir-small-amplitude")) {
        const Metrics r = run_langmuir_small_amplitude();
        langmuir_ok = r.all_ok && r.max_ledger_diff <= 1.0e-10;
        m.max_ledger_diff = std::max(m.max_ledger_diff, r.max_ledger_diff);
    }
    if (ok && (args.test_case == "all" ||
               args.test_case == "left-e-reservoir")) {
        const Metrics r = run_left_e_reservoir();
        left_e_ok = r.all_ok && r.left_e_number_rel <= 1.0e-9;
        m.left_e_max_field_energy = r.left_e_max_field_energy;
        m.left_e_number_rel = r.left_e_number_rel;
    }
    pass = pass && maxwell_ok && packet_ok && straddle_ok && langmuir_ok &&
           left_e_ok;
    m.all_ok = pass;

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "all_ok=" << (m.all_ok ? 1 : 0)
              << " max_ledger_diff=" << m.max_ledger_diff
              << " combined_number_rel=" << m.combined_number_rel
              << " combined_energy_rel=" << m.combined_energy_rel
              << " conversion_residual_max=" << m.conversion_residual_max
              << " tail_number_balance_max=" << m.tail_number_balance_max
              << " field_energy_rel_change=" << m.field_energy_rel_change
              << " left_e_max_field_energy=" << m.left_e_max_field_energy
              << " left_e_number_rel=" << m.left_e_number_rel << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
