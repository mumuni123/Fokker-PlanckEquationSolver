// Stage H1 acceptance for tail CIC deposition and step accounting
// (sections 6.3, 6.5, 6.6.2 and 15 H1):
//   * cic-charge-conservation: deposited density integrates to the total
//     weight and the cell weights sum to one;
//   * continuity-open-outflow-conversion: per-cell discrete continuity
//     closes to summation error with open outflow and an artificial
//     conversion source;
//   * trial-failure-safety: a rejected trial (NaN kick, invalid weight)
//     never mutates the accepted state.
//
// Usage:
//   background_tail_deposition_test --case <case> [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "grid.h"
#include "maxwell.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
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

unsigned long long lcg_state = 0x123456789abcdef0ULL;

double random_unit()
{
    lcg_state = lcg_state * 2862933555777941757ULL + 3037000493ULL;
    return static_cast<double>(lcg_state >> 11) * (1.0 / 9007199254740992.0);
}

struct Snapshot {
    std::vector<BackgroundTailParticle> particles;
    std::vector<double> density;
    std::uint64_t id_counter;
    TailOutflowLedger outflow;
    double truncation_shape_left;
    double truncation_shape_right;
};

Snapshot capture(const BackgroundTailPIC& tail)
{
    Snapshot s;
    s.particles = tail.particles;
    s.density = tail.density;
    s.id_counter = tail.id_counter();
    s.outflow = tail.outflow_ledger();
    s.truncation_shape_left = tail.truncation_shape_left();
    s.truncation_shape_right = tail.truncation_shape_right();
    return s;
}

bool identical(const Snapshot& a, const Snapshot& b)
{
    if (a.particles.size() != b.particles.size()) return false;
    for (size_t i = 0; i < a.particles.size(); ++i) {
        const BackgroundTailParticle& pa = a.particles[i];
        const BackgroundTailParticle& pb = b.particles[i];
        if (pa.x != pb.x || pa.ux != pb.ux || pa.uy != pb.uy ||
            pa.uz != pb.uz || pa.weight != pb.weight || pa.id != pb.id) {
            return false;
        }
    }
    if (a.density != b.density) return false;
    if (a.id_counter != b.id_counter) return false;
    if (a.outflow.left_number != b.outflow.left_number ||
        a.outflow.right_number != b.outflow.right_number ||
        a.outflow.left_px != b.outflow.left_px ||
        a.outflow.right_px != b.outflow.right_px ||
        a.outflow.left_kinetic_energy != b.outflow.left_kinetic_energy ||
        a.outflow.right_kinetic_energy != b.outflow.right_kinetic_energy) {
        return false;
    }
    if (a.truncation_shape_left != b.truncation_shape_left ||
        a.truncation_shape_right != b.truncation_shape_right) {
        return false;
    }
    return true;
}

struct Metrics {
    double cic_total_charge_rel;
    double cic_weight_balance_rel;
    double max_continuity_rel_l1;
    double max_continuity_rel_linf;
    double max_truncation_ledger_error;
    double max_number_balance_error;
    double source_sensitivity_ratio;
    double final_weight_balance_rel;
    bool trial_accepted_unchanged;
    bool trial_rejected;
    bool negative_weight_rejected;
    bool swap_roundtrip;
    Metrics()
        : cic_total_charge_rel(0.0), cic_weight_balance_rel(0.0),
          max_continuity_rel_l1(0.0), max_continuity_rel_linf(0.0),
          max_truncation_ledger_error(0.0), max_number_balance_error(0.0),
          source_sensitivity_ratio(0.0), final_weight_balance_rel(0.0),
          trial_accepted_unchanged(false), trial_rejected(false),
          negative_weight_rejected(false), swap_roundtrip(false)
    {}
};

Metrics run_cic_charge_conservation()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    const int n = 500;
    double total_weight = 0.0;
    for (int i = 0; i < n; ++i) {
        BackgroundTailParticle p;
        p.x = (0.2 + 0.6 * random_unit()) *
              (grid.dx * static_cast<double>(grid.nx_global));
        p.ux = (random_unit() - 0.5) * 4.0;
        p.uy = (random_unit() - 0.5) * 2.0;
        p.uz = (random_unit() - 0.5) * 2.0;
        p.weight = 1.0e20 * (0.5 + random_unit());
        p.id = tail.next_particle_id(0);
        tail.particles.push_back(p);
        total_weight += p.weight;
    }
    tail.deposit_density(grid, 0, 1);
    double deposited = 0.0;
    for (size_t i = 0; i < tail.density.size(); ++i) {
        deposited += tail.density[i] * grid.dx;
    }
    m.cic_total_charge_rel =
        std::fabs(deposited - total_weight) / std::max(1.0, total_weight);
    m.cic_weight_balance_rel =
        tail.density_weight_balance(grid) / std::max(1.0, total_weight);
    // Cell weights always sum to one (a CIC structural invariant).
    double cell_weight_sum_err = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const double x = random_unit() *
                         (grid.dx * static_cast<double>(grid.nx_global));
        const CellDepositWeights cw = ParticleShape1D::cell_weights(x, grid);
        cell_weight_sum_err = std::max(
            cell_weight_sum_err,
            std::fabs(cw.w0 + cw.w1 - 1.0));
    }
    if (cell_weight_sum_err > 1.0e-15) {
        m.cic_total_charge_rel = 1.0;  // force failure
    }
    return m;
}

Metrics run_continuity_open_outflow_conversion()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    const double dt = 1.0e-16;
    BackgroundTailPIC tail;
    tail.init(grid);

    const double w = 1.0e20;
    const double positions[5] = { 0.15, 1.85, 0.5, 1.2, 0.8 };
    const double ux[5] = { -3.0, 3.0, 1.0, -0.5, 2.0 };
    for (int i = 0; i < 5; ++i) {
        BackgroundTailParticle p;
        p.x = positions[i] * Const::micro;
        p.ux = ux[i];
        p.uy = 0.1 * static_cast<double>(i);
        p.uz = -0.1 * static_cast<double>(i);
        p.weight = w;
        p.id = tail.next_particle_id(0);
        tail.particles.push_back(p);
    }

    // Artificial conversion source density [m^-3 s^-1] in a mid-domain
    // block: positive into the tail.
    const size_t nxl = static_cast<size_t>(grid.nx_local);
    std::vector<double> source(nxl, 0.0);
    const double s0 = 1.0e33;
    for (size_t i = 45; i <= 65 && i < nxl; ++i) {
        const double d = (static_cast<double>(i) - 55.0) / 4.0;
        source[i] = s0 * std::exp(-d * d);
    }

    tail.deposit_density(grid, 0, 1);
    const int steps = 8;
    double total_outflow = 0.0;
    EMFields zero_fields;
    zero_fields.init(grid);
    for (int s = 0; s < steps; ++s) {
        tail.begin_step(grid, dt);
        tail.drift_half(grid, 0.5 * dt, 0, 1);
        tail.kick(grid, zero_fields, dt, 0, 1);
        tail.drift_half(grid, 0.5 * dt, 0, 1);
        tail.deposit_density(grid, 0, 1);
        tail.finalize_trajectory_current(grid, dt, 0, 1);
        // Materialize this step's artificial conversion source: one particle
        // per source cell at the cell centre with weight = source*dx*dt, so
        // the audited density change matches the source term.
        for (size_t i = 0; i < source.size(); ++i) {
            if (!(source[i] > 0.0)) continue;
            const double added_weight = source[i] * grid.dx * dt;
            BackgroundTailParticle sp;
            sp.x = (static_cast<double>(i) + 0.5) * grid.dx;
            sp.ux = 0.0;
            sp.uy = 0.0;
            sp.uz = 0.0;
            sp.weight = added_weight;
            sp.id = tail.next_particle_id(0);
            tail.particles.push_back(sp);
            tail.density[i] += source[i] * dt;
        }
        const TailContinuityResult r =
            tail.audit_continuity(grid, dt, source, 0, 1);
        m.max_continuity_rel_l1 =
            std::max(m.max_continuity_rel_l1, r.rel_l1);
        m.max_continuity_rel_linf =
            std::max(m.max_continuity_rel_linf, r.rel_linf);
        m.max_truncation_ledger_error = std::max(
            m.max_truncation_ledger_error, r.truncation_ledger_error);
        m.max_number_balance_error = std::max(
            m.max_number_balance_error, r.number_balance_error);
        total_outflow += tail.outflow_ledger().left_number +
                         tail.outflow_ledger().right_number;
    }
    // The audit must be sensitive to the source term: without C the residual
    // should be far larger than the closed residual.
    tail.begin_step(grid, dt);
    tail.drift_half(grid, 0.5 * dt, 0, 1);
    tail.kick(grid, zero_fields, dt, 0, 1);
    tail.drift_half(grid, 0.5 * dt, 0, 1);
    tail.deposit_density(grid, 0, 1);
    tail.finalize_trajectory_current(grid, dt, 0, 1);
    for (size_t i = 0; i < source.size(); ++i) {
        if (!(source[i] > 0.0)) continue;
        const double added_weight = source[i] * grid.dx * dt;
        BackgroundTailParticle sp;
        sp.x = (static_cast<double>(i) + 0.5) * grid.dx;
        sp.ux = 0.0;
        sp.uy = 0.0;
        sp.uz = 0.0;
        sp.weight = added_weight;
        sp.id = tail.next_particle_id(0);
        tail.particles.push_back(sp);
        tail.density[i] += source[i] * dt;
    }
    const TailContinuityResult r_closed =
        tail.audit_continuity(grid, dt, source, 0, 1);
    const TailContinuityResult r_without_source =
        tail.audit_continuity(grid, dt, std::vector<double>(), 0, 1);
    m.source_sensitivity_ratio =
        r_without_source.abs_l1 / std::max(1e-30, r_closed.abs_l1);

    const double total_weight = 5.0 * w;
    m.final_weight_balance_rel =
        tail.density_weight_balance(grid) / total_weight;
    (void)total_outflow;
    return m;
}

Metrics run_trial_failure_safety()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    const double dt = 1.0e-16;

    BackgroundTailPIC accepted;
    accepted.init(grid);
    BackgroundTailParticle p1;
    p1.x = 0.5 * Const::micro;
    p1.ux = 1.0;
    p1.uy = 0.2;
    p1.uz = 0.0;
    p1.weight = 1.0e20;
    p1.id = accepted.next_particle_id(0);
    accepted.particles.push_back(p1);
    BackgroundTailParticle p2;
    p2.x = 1.2 * Const::micro;
    p2.ux = -0.5;
    p2.uy = 0.0;
    p2.uz = 0.3;
    p2.weight = 2.0e20;
    p2.id = accepted.next_particle_id(1);
    accepted.particles.push_back(p2);
    accepted.deposit_density(grid, 0, 1);
    accepted.begin_step(grid, dt);
    const Snapshot before = capture(accepted);

    // NaN kick on the trial: the trial must be rejected and the accepted
    // state must be bitwise unchanged.
    BackgroundTailPIC trial = accepted;
    EMFields bad_fields;
    bad_fields.init(grid);
    std::fill(bad_fields.Ex_face.begin(), bad_fields.Ex_face.end(),
              std::numeric_limits<double>::quiet_NaN());
    trial.kick(grid, bad_fields, dt, 0, 1);
    m.trial_rejected = !trial.finite();
    m.trial_accepted_unchanged = identical(before, capture(accepted));

    // Invalid (non-positive) macro weight must be rejected by the guard.
    BackgroundTailPIC bad = accepted;
    BackgroundTailParticle bad_p;
    bad_p.x = 0.3 * Const::micro;
    bad_p.ux = 0.0;
    bad_p.uy = 0.0;
    bad_p.uz = 0.0;
    bad_p.weight = -1.0;
    bad_p.id = 99;
    bad.particles.push_back(bad_p);
    m.negative_weight_rejected = !bad.nonnegative_weights();

    // Transactional swap: the two objects exchange the complete state and a
    // second swap restores the accepted state.
    BackgroundTailPIC swapped = accepted;
    BackgroundTailPIC holder = accepted;
    holder.particles[0].ux += 0.5;  // distinct trial state
    const Snapshot holder_before = capture(holder);
    swapped.swap_state(holder);
    const Snapshot accepted_after_first = capture(swapped);
    const Snapshot holder_after_first = capture(holder);
    swapped.swap_state(holder);
    const Snapshot accepted_after_roundtrip = capture(swapped);
    const Snapshot holder_after_roundtrip = capture(holder);
    m.swap_roundtrip =
        identical(accepted_after_first, holder_before) &&
        identical(holder_after_first, before) &&
        identical(accepted_after_roundtrip, before) &&
        identical(holder_after_roundtrip, holder_before);
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "cic_total_charge_rel=" << m.cic_total_charge_rel << "\n";
    out << "cic_weight_balance_rel=" << m.cic_weight_balance_rel << "\n";
    out << "max_continuity_rel_l1=" << m.max_continuity_rel_l1 << "\n";
    out << "max_continuity_rel_linf=" << m.max_continuity_rel_linf << "\n";
    out << "max_truncation_ledger_error=" << m.max_truncation_ledger_error
        << "\n";
    out << "max_number_balance_error=" << m.max_number_balance_error << "\n";
    out << "source_sensitivity_ratio=" << m.source_sensitivity_ratio << "\n";
    out << "final_weight_balance_rel=" << m.final_weight_balance_rel << "\n";
    out << "trial_accepted_unchanged="
        << (m.trial_accepted_unchanged ? 1 : 0) << "\n";
    out << "trial_rejected=" << (m.trial_rejected ? 1 : 0) << "\n";
    out << "negative_weight_rejected="
        << (m.negative_weight_rejected ? 1 : 0) << "\n";
    out << "swap_roundtrip=" << (m.swap_roundtrip ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_deposition_test must run with exactly "
                     "1 rank; use plain ./build_hybrid/background_tail_"
                     "deposition_test (no yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: background_tail_deposition_test --case <case> "
                     "[--result <path>]\n"
                  << "cases: cic-charge-conservation | "
                     "continuity-open-outflow-conversion | "
                     "trial-failure-safety | all\n";
    }

    Metrics m;
    bool pass = ok;
    const bool do_cic =
        ok && (args.test_case == "all" ||
               args.test_case == "cic-charge-conservation");
    const bool do_cont =
        ok && (args.test_case == "all" ||
               args.test_case == "continuity-open-outflow-conversion");
    const bool do_trial =
        ok && (args.test_case == "all" ||
               args.test_case == "trial-failure-safety");
    if (do_cic) m = run_cic_charge_conservation();
    if (do_cont) {
        Metrics c = run_continuity_open_outflow_conversion();
        m.max_continuity_rel_l1 =
            std::max(m.max_continuity_rel_l1, c.max_continuity_rel_l1);
        m.max_continuity_rel_linf =
            std::max(m.max_continuity_rel_linf, c.max_continuity_rel_linf);
        m.max_truncation_ledger_error = std::max(
            m.max_truncation_ledger_error, c.max_truncation_ledger_error);
        m.max_number_balance_error = std::max(
            m.max_number_balance_error, c.max_number_balance_error);
        m.source_sensitivity_ratio = c.source_sensitivity_ratio;
        m.final_weight_balance_rel = c.final_weight_balance_rel;
    }
    if (do_trial) {
        Metrics t = run_trial_failure_safety();
        m.trial_accepted_unchanged = t.trial_accepted_unchanged;
        m.trial_rejected = t.trial_rejected;
        m.negative_weight_rejected = t.negative_weight_rejected;
        m.swap_roundtrip = t.swap_roundtrip;
    }

    const bool cic_ok =
        !do_cic ||
        (m.cic_total_charge_rel <= 1.0e-12 &&
         m.cic_weight_balance_rel <= 1.0e-12);
    const bool cont_ok =
        !do_cont ||
        (m.max_continuity_rel_l1 <= 1.0e-9 &&
         m.max_continuity_rel_linf <= 1.0e-9 &&
         m.max_truncation_ledger_error <= 1.0e-9 &&
         m.max_number_balance_error <= 1.0e-9 &&
         m.source_sensitivity_ratio >= 100.0 &&
         m.final_weight_balance_rel <= 1.0e-9);
    const bool trial_ok =
        !do_trial ||
        (m.trial_rejected && m.trial_accepted_unchanged &&
         m.negative_weight_rejected && m.swap_roundtrip);
    pass = pass && cic_ok && cont_ok && trial_ok;

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "cic_total_charge_rel=" << m.cic_total_charge_rel
              << " cic_weight_balance_rel=" << m.cic_weight_balance_rel
              << " max_continuity_rel_l1=" << m.max_continuity_rel_l1
              << " max_continuity_rel_linf=" << m.max_continuity_rel_linf
              << " max_truncation_ledger_error="
              << m.max_truncation_ledger_error
              << " max_number_balance_error=" << m.max_number_balance_error
              << " source_sensitivity_ratio=" << m.source_sensitivity_ratio
              << " final_weight_balance_rel=" << m.final_weight_balance_rel
              << " trial_accepted_unchanged="
              << (m.trial_accepted_unchanged ? 1 : 0)
              << " trial_rejected=" << (m.trial_rejected ? 1 : 0)
              << " negative_weight_rejected="
              << (m.negative_weight_rejected ? 1 : 0)
              << " swap_roundtrip=" << (m.swap_roundtrip ? 1 : 0) << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
