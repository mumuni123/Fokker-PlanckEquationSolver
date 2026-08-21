// Gate I acceptance test (section 4.7.3) for the four discrete pairing
// identities and the injected-fault negative test (section 4.7.5).  All
// production operators are driven directly (VlasovSplitStep, BeamPIC,
// BackgroundTailPIC, FieldParticlePowerAudit); no production formula is
// re-implemented here.
//
// Usage:
//   vpfp_field_particle_pairing_test --case all [--result <path>]
//   vpfp_field_particle_pairing_test --case injected-faults [--result <path>]

#include "conservative_ppm_remap.h"
#include "field_particle_power_audit.h"
#include "open_boundary.h"
#include "parameters.h"
#include "vlasov_split_step.h"
#include "vpfp_field_particle_pairing_test_support.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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

struct PairingMetrics {
    double continuity_residual_abs;
    double continuity_scale;
    double continuity_tolerance;
    bool continuity_pass;
    double poisson_transport_residual_abs;
    double poisson_transport_scale;
    bool poisson_transport_pass;
    double force_work_residual_abs;
    double force_work_scale;
    bool force_work_pass;
    double current_pair_residual_abs;
    double current_pair_scale;
    bool current_pair_pass;
    double poisson_identity_crosscheck;
    double conversion_transfer_residual_abs;
    double conversion_transfer_scale;
    double boundary_residual_abs;
    double boundary_scale;
    double full_residual_abs;
    double reconstructed_residual_abs;
    double reconstruction_mismatch;
    double reconstruction_tolerance;
    bool reconstruction_pass;
    double bulk_work_plus;
    double bulk_work_minus;
    double tail_cell_work_mismatch;
    double beam_cell_work_mismatch;
    double manufactured_face_work;
    double manufactured_cell_work;
    double manufactured_identity_error;
    bool manufactured_identity_pass;
    bool manufactured_endpoint_weight_pass;
    bool read_only_bitwise_equal;
    int negative_test_effective;
    int injected_detected_count;
    int injected_total_count;
    int injected_detected[6];
    int injected_first_bad[6];
    PairingMetrics()
        : continuity_residual_abs(0.0), continuity_scale(1.0),
          continuity_tolerance(0.0), continuity_pass(false),
          poisson_transport_residual_abs(0.0), poisson_transport_scale(1.0),
          poisson_transport_pass(false), force_work_residual_abs(0.0),
          force_work_scale(1.0), force_work_pass(false),
          current_pair_residual_abs(0.0), current_pair_scale(1.0),
          current_pair_pass(false), poisson_identity_crosscheck(0.0),
          conversion_transfer_residual_abs(0.0),
          conversion_transfer_scale(1.0), boundary_residual_abs(0.0),
          boundary_scale(1.0), full_residual_abs(0.0),
          reconstructed_residual_abs(0.0), reconstruction_mismatch(0.0),
          reconstruction_tolerance(0.0), reconstruction_pass(false),
          bulk_work_plus(0.0), bulk_work_minus(0.0),
          tail_cell_work_mismatch(0.0), beam_cell_work_mismatch(0.0),
          manufactured_face_work(0.0), manufactured_cell_work(0.0),
          manufactured_identity_error(0.0), manufactured_identity_pass(false),
          manufactured_endpoint_weight_pass(false),
          read_only_bitwise_equal(false), negative_test_effective(1),
          injected_detected_count(0), injected_total_count(0)
    {
        for (int i = 0; i < 6; ++i) {
            injected_detected[i] = 0;
            injected_first_bad[i] = -1;
        }
    }
};

// Manufactured shared-face/endpoint audit.  The production audit owns the
// gather and quadrature definitions; this test only supplies deterministic
// E_face and J_cell data and checks the resulting identity.
void run_manufactured_gstar_identity(PairingMetrics& m)
{
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 32, Param::Lx);
    FieldParticlePowerAudit audit;
    audit.init(grid);
    std::vector<double> e_face(static_cast<size_t>(grid.nx_local) + 1, 0.0);
    std::vector<double> j_cell(static_cast<size_t>(grid.nx_local), 0.0);
    for (size_t f = 0; f < e_face.size(); ++f) {
        const double x = static_cast<double>(f) / grid.nx_global;
        e_face[f] = 1.0e7 * (1.0 + 0.37 * std::sin(2.0 * Const::pi * x) +
                             0.11 * std::cos(6.0 * Const::pi * x));
    }
    for (size_t i = 0; i < j_cell.size(); ++i) {
        const double x = (static_cast<double>(i) + 0.5) / grid.nx_global;
        j_cell[i] = 2.0e5 * (1.0 - 0.23 * std::cos(4.0 * Const::pi * x));
    }
    const bool finite = audit.manufactured_gstar_identity(
        e_face, j_cell, m.manufactured_face_work,
        m.manufactured_cell_work, m.manufactured_identity_error);
    const double scale = std::max(
        1.0, std::max(std::fabs(m.manufactured_face_work),
                      std::fabs(m.manufactured_cell_work)));
    m.manufactured_identity_pass = finite &&
        m.manufactured_identity_error <= machine_scaled_tolerance(scale, scale);
    const double weight_tol = machine_scaled_tolerance(grid.dx, grid.dx);
    m.manufactured_endpoint_weight_pass =
        std::fabs(audit.face_quadrature_weight(0) - 0.5 * grid.dx) <=
            weight_tol &&
        std::fabs(audit.face_quadrature_weight(grid.nx_local) -
                  0.5 * grid.dx) <= weight_tol &&
        std::fabs(audit.face_quadrature_weight(grid.nx_local / 2) -
                  grid.dx) <= weight_tol;
}

// Case 1/10: bulk-only continuity through the production x half-steps, and
// the four-identity reconstruction of the pairing calculator on a smooth
// manufactured state.
void run_bulk_continuity_and_reconstruction(PairingMetrics& m)
{
    PairingTestState s = make_smooth_bulk_case(32, false, false);
    const SpatialGrid& grid = s.grid;
    const int nxl = grid.nx_local;
    const double dt = 1.0e-15;

    OpenBackgroundBoundaryConfig cfg;
    cfg.left_type = BackgroundXBoundaryType::ABSORBING;
    cfg.right_type = BackgroundXBoundaryType::ABSORBING;
    OpenBackgroundBoundary boundary(cfg);

    VlasovSplitStep vlasov;
    vlasov.init(grid, s.bulk_n, boundary);
    Species half;
    half.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
              Const::me, Param::dens, Param::temperature_e, false, grid);
    VlasovStepDiagnostics diag;
    XFaceTransportAudit x1;
    x1.enabled = true;
    x1.init(nxl);
    vlasov.first_x_half(s.bulk_n, half, 0.0, 0.5 * dt, diag, &x1);

    // Build a pairing workspace and run the calculator to exercise the four
    // identities directly (section 4.5).
    FieldParticlePowerAuditWorkspace ws;
    ws.enabled = true;
    ws.bulk_x1 = x1;
    ws.bulk_x2.init(nxl);
    ws.cell_work.init(nxl);
    ws.bulk_number_n.assign(static_cast<size_t>(nxl), 0.0);
    ws.bulk_number_np1.assign(static_cast<size_t>(nxl), 0.0);
    for (int i = 0; i < nxl; ++i) {
        double m_n = 0.0;
        double m_np1 = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                m_n += s.bulk_n.f[idx3(grid.nghost + i, j, k)];
                m_np1 += half.f[idx3(grid.nghost + i, j, k)];
            }
        }
        ws.bulk_number_n[static_cast<size_t>(i)] = m_n;
        ws.bulk_number_np1[static_cast<size_t>(i)] = m_np1;
    }
    ws.field_n_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.field_mid_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.field_np1_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.potential_pair_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.bulk_work_ledger = 0.0;
    ws.tail_work_ledger = 0.0;
    ws.beam_work_ledger = 0.0;

    FieldParticlePowerAudit calculator;
    calculator.init(grid);
    const FieldParticlePowerAuditResult r =
        calculator.finalize(ws, dt, -Const::qe);

    m.continuity_residual_abs = r.continuity_linf;
    m.continuity_scale = std::max(1.0, r.continuity_scale);
    m.continuity_tolerance =
        machine_scaled_tolerance(m.continuity_scale, m.continuity_scale);
    m.continuity_pass = r.continuity_pass;
    m.poisson_transport_residual_abs =
        std::fabs(r.poisson_transport_residual);
    m.poisson_transport_scale = std::max(1.0,
        std::fabs(r.poisson_transport_residual));
    m.poisson_transport_pass =
        std::isfinite(r.poisson_transport_residual) &&
        m.poisson_transport_residual_abs <= machine_scaled_tolerance(
            m.poisson_transport_scale, m.poisson_transport_scale);
    m.force_work_residual_abs = std::fabs(r.force_work_residual);
    m.force_work_scale = std::max(1.0, m.force_work_residual_abs);
    m.force_work_pass = std::isfinite(r.force_work_residual) &&
        m.force_work_residual_abs <= machine_scaled_tolerance(
            m.force_work_scale, m.force_work_scale);
    m.current_pair_residual_abs = std::fabs(r.current_pair_residual);
    m.current_pair_scale = std::max(1.0, m.current_pair_residual_abs);
    m.current_pair_pass = std::isfinite(r.current_pair_residual) &&
        m.current_pair_residual_abs <= machine_scaled_tolerance(
            m.current_pair_scale, m.current_pair_scale);
    m.poisson_identity_crosscheck = r.poisson_identity_crosscheck;
    m.conversion_transfer_residual_abs =
        std::fabs(r.conversion_transfer_residual);
    m.conversion_transfer_scale = std::max(1.0, m.conversion_transfer_residual_abs);
    m.boundary_residual_abs = std::fabs(r.boundary_residual);
    m.boundary_scale = std::max(1.0, m.boundary_residual_abs);
    m.full_residual_abs = std::fabs(r.full_residual);
    m.reconstructed_residual_abs = std::fabs(r.reconstructed_full_residual);
    m.reconstruction_mismatch = r.reconstruction_mismatch;
    m.reconstruction_tolerance = r.roundoff_tolerance;
    m.reconstruction_pass = r.full_residual_reconstruction_pass;
}

// Case 2: bulk force work with +/- E through advect_u_parallel; the per-cell
// work sums to the existing bulk upar face work and reverses sign with E.
void run_bulk_force(PairingMetrics& m)
{
    PairingTestState s = make_smooth_bulk_case(32, false, false);
    const SpatialGrid& grid = s.grid;
    const int nxl = grid.nx_local;
    const double dt = 1.0e-15;

    // A smooth-in-x, one-sided positive-u distribution gives an unambiguous
    // first-order work sign.  A zero-drift Maxwellian is even in u and its
    // acceleration energy is even in E, so it cannot test sign reversal.
    std::vector<double> drifting_slice(Param::Nvmu, 0.0);
    const int positive_u = 3 * Param::Nv / 4;
    const int low_uperp = 0;
    drifting_slice[static_cast<size_t>(positive_u) * Param::Nmu +
                   static_cast<size_t>(low_uperp)] = Param::dens * grid.dx;
    for (int ix = 0; ix < grid.nx_total; ++ix) {
        const double x = grid.x(ix);
        const double amp = 1.0 + 1.0e-3 *
            std::cos(2.0 * Const::pi * x / grid.length());
        const size_t base = static_cast<size_t>(ix) * Param::Nvmu;
        for (size_t q = 0; q < drifting_slice.size(); ++q)
            s.bulk_n.f[base + q] = amp * drifting_slice[q];
    }

    ConservativePpmRemap remap;
    remap.init(grid, s.velocity_grid);
    EMFields field;
    field.init(grid);
    std::fill(field.Ex.begin(), field.Ex.end(), 0.0);
    std::fill(field.Ex_face.begin(), field.Ex_face.end(), 0.0);

    for (int sign = 0; sign < 2; ++sign) {
        const double e_x = sign == 0 ? 1.0e12 : -1.0e12;
        std::fill(field.Ex.begin(), field.Ex.end(), e_x);
        std::fill(field.Ex_face.begin(), field.Ex_face.end(), e_x);
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        std::vector<double> cell_work(static_cast<size_t>(nxl), 0.0);
        const RemapDiagnostics d = remap.advect_u_parallel(
            s.bulk_n, out, field, dt, 0.0, NULL, NULL, 4, &cell_work);
        const double sum = std::accumulate(cell_work.begin(), cell_work.end(), 0.0);
        if (sign == 0) m.bulk_work_plus = sum;
        else m.bulk_work_minus = sum;
        (void)d;
    }
}

// Case 4/6: Tail and Beam kick cell-work sum matches the existing global kick
// work (section 4.4 item 3).
void run_pic_kick_work(PairingMetrics& m)
{
    const int nx = 32;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, nx, Param::Lx);
    const double dt = 1.0e-15;

    {
        BackgroundTailPIC tail;
        tail.init(grid);
        BackgroundTailParticle p = {};
        p.x = 0.5 * grid.length();
        p.ux = 0.5;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0e14;
        p.id = 1;
        tail.particles.push_back(p);
        EMFields field;
        field.init(grid);
        std::fill(field.Ex.begin(), field.Ex.end(), -1.0e12);
        std::fill(field.Ex_face.begin(), field.Ex_face.end(), -1.0e12);
        double work = 0.0;
        std::vector<double> cell_work(static_cast<size_t>(grid.nx_local), 0.0);
        tail.kick(grid, field, dt, 0, 1, &work, &cell_work);
        const double sum = std::accumulate(cell_work.begin(), cell_work.end(), 0.0);
        m.tail_cell_work_mismatch = std::fabs(sum - work);
    }
    {
        BeamPIC beam;
        beam.init(grid);
        BeamParticle p;
        p.x = 0.5 * grid.length();
        p.px = 0.5 * Const::me * Const::c;
        p.weight = 1.0e12;
        beam.particles.push_back(p);
        beam.begin_step(grid, dt);
        EMFields field;
        field.init(grid);
        std::fill(field.Ex.begin(), field.Ex.end(), -1.0e11);
        std::fill(field.Ex_face.begin(), field.Ex_face.end(), -1.0e11);
        BeamInjectionSchedule empty;
        BeamPIC work_beam;
        work_beam.init(grid);
        work_beam.begin_step(grid, dt);
        beam.predict_to_midpoint(empty, grid, field, 0.0, dt, 0, 1, work_beam);
        std::vector<double> cell_work(static_cast<size_t>(grid.nx_local), 0.0);
        work_beam.finish_from_midpoint(empty, grid, field, 0.0, dt, 0, 1,
                                       &cell_work);
        const double sum = std::accumulate(cell_work.begin(), cell_work.end(), 0.0);
        m.beam_cell_work_mismatch =
            std::fabs(sum - work_beam.last_field_work());
    }
}

// Case 10: read-only equivalence of the bulk x remap audit on/off.
void run_read_only_equivalence(PairingMetrics& m)
{
    PairingTestState a = make_smooth_bulk_case(32, false, false);
    PairingTestState b = make_smooth_bulk_case(32, false, false);
    const SpatialGrid& grid = a.grid;
    const double dt = 1.0e-15;
    OpenBackgroundBoundaryConfig cfg;
    cfg.left_type = BackgroundXBoundaryType::ABSORBING;
    cfg.right_type = BackgroundXBoundaryType::ABSORBING;
    OpenBackgroundBoundary boundary(cfg);

    ConservativePpmRemap remap;
    remap.init(grid, a.velocity_grid);
    Species out_a;
    Species out_b;
    out_a.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
               Const::me, Param::dens, Param::temperature_e, false, grid);
    out_b.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
               Const::me, Param::dens, Param::temperature_e, false, grid);
    XFaceTransportAudit audit;
    audit.enabled = true;
    audit.init(grid.nx_local);
    remap.advect_x(a.bulk_n, out_a, dt, 0.0, boundary, 0, 1, &audit);
    remap.advect_x(b.bulk_n, out_b, dt, 0.0, boundary, 0, 1, NULL);
    m.read_only_bitwise_equal = (out_a.f == out_b.f);
}

// Section 4.7.5: injected faults must each flip the corresponding gate.
// Faults are applied to an audit-input copy, never to production state.
bool run_injected_fault(PairingMetrics& m, int fault_id)
{
    PairingTestState s = make_smooth_bulk_case(32, false, false);
    const SpatialGrid& grid = s.grid;
    const int nxl = grid.nx_local;
    const double dt = 1.0e-15;

    FieldParticlePowerAuditWorkspace ws;
    ws.enabled = true;
    ws.bulk_x1.init(nxl);
    ws.bulk_x2.init(nxl);
    ws.cell_work.init(nxl);
    ws.bulk_number_n.assign(static_cast<size_t>(nxl), 1.0e18);
    ws.bulk_number_np1.assign(static_cast<size_t>(nxl), 1.0e18);
    ws.field_n_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.field_mid_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.field_np1_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.potential_pair_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    // Give the swept faces a small known pattern and construct the clean
    // endpoint state from that exact finite-volume update.
    for (int f = 0; f <= nxl; ++f) {
        ws.bulk_x1.bulk_number_swept_face[static_cast<size_t>(f)] =
            1.0e14 * static_cast<double>(f % 5);
    }
    ws.bulk_x2.bulk_number_swept_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    ws.bulk_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.tail_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.beam_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    // A right-boundary source exercises the independent source ledger.
    ws.bulk_source_number_cell.back() = 2.0e13;
    for (int i = 0; i < nxl; ++i) {
        ws.bulk_number_np1[static_cast<size_t>(i)] =
            ws.bulk_number_n[static_cast<size_t>(i)] +
            ws.bulk_x1.bulk_number_swept_face[static_cast<size_t>(i)] -
            ws.bulk_x1.bulk_number_swept_face[static_cast<size_t>(i + 1)] +
            ws.bulk_source_number_cell[static_cast<size_t>(i)];
    }
    // Nonzero, internally consistent cell/ledger work makes a misplaced cell
    // visible through current_pair_linf while preserving the clean total.
    ws.cell_work.bulk_delta_ke_cell[static_cast<size_t>(nxl / 2)] = 1.0e-6;
    ws.bulk_work_ledger = 1.0e-6;
    ws.electrode_work = 3.0e-6;
    ws.poisson_potential_charge_work = -3.0e-6;
    ws.poisson_identity_residual = 0.0;

    FieldParticlePowerAudit calculator;
    calculator.init(grid);
    FieldParticlePowerAuditResult baseline =
        calculator.finalize(ws, dt, -Const::qe);

    // Apply the fault to a copy of the raw audit inputs only.
    FieldParticlePowerAuditWorkspace bad = ws;
    if (fault_id == 1) {
        // Cycle-shift the J_charge face array by one cell.
        std::rotate(bad.bulk_x1.bulk_number_swept_face.begin(),
                    bad.bulk_x1.bulk_number_swept_face.begin() + 1,
                    bad.bulk_x1.bulk_number_swept_face.end());
    } else if (fault_id == 2) {
        // Flip the sign of one interior face.
        bad.bulk_x1.bulk_number_swept_face[1] *= -1.0;
    } else if (fault_id == 3) {
        // Count one shared/interior face twice.
        bad.bulk_x1.bulk_number_swept_face[static_cast<size_t>(nxl / 2)] *= 2.0;
    } else if (fault_id == 4) {
        // Drop the independent right-boundary source.
        bad.bulk_source_number_cell.back() = 0.0;
    } else if (fault_id == 5) {
        // Misplace J_work_bulk into an adjacent cell.
        const size_t from = static_cast<size_t>(nxl / 2);
        const size_t to = from + 1;
        const double moved = 0.5 * bad.cell_work.bulk_delta_ke_cell[from];
        bad.cell_work.bulk_delta_ke_cell[from] -= moved;
        bad.cell_work.bulk_delta_ke_cell[to] += moved;
    } else if (fault_id == 6) {
        // Add the electrode work a second time.
        bad.electrode_work *= 2.0;
    }
    const FieldParticlePowerAuditResult faulty =
        calculator.finalize(bad, dt, -Const::qe);

    ++m.injected_total_count;
    // A fault is "detected" when the reconstruction mismatch or continuity
    // linf differs from the clean baseline by more than the tolerance.
    const double delta_cont = std::fabs(
        faulty.continuity_linf - baseline.continuity_linf);
    const double delta_rec = std::fabs(
        faulty.reconstruction_mismatch - baseline.reconstruction_mismatch);
    const double delta_local = std::fabs(
        faulty.current_pair_linf - baseline.current_pair_linf);
    const double delta_identity = std::fabs(
        faulty.poisson_identity_crosscheck -
        baseline.poisson_identity_crosscheck);
    const double scale = std::max(1.0, baseline.reconstruction_scale);
    const double tol = machine_scaled_tolerance(scale, scale);
    const bool detected = delta_cont > tol || delta_rec > tol ||
                          delta_local > tol || delta_identity > tol;
    if (detected) ++m.injected_detected_count;
    if (fault_id >= 1 && fault_id <= 6) {
        m.injected_detected[fault_id - 1] = detected ? 1 : 0;
        m.injected_first_bad[fault_id - 1] = faulty.first_bad_index;
    }
    return detected;
}

bool write_result_file(const std::string& path, const PairingMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "continuity_residual_abs=" << m.continuity_residual_abs << "\n";
    out << "continuity_scale=" << m.continuity_scale << "\n";
    out << "continuity_tolerance=" << m.continuity_tolerance << "\n";
    out << "continuity_pass=" << (m.continuity_pass ? 1 : 0) << "\n";
    out << "poisson_transport_residual_abs="
        << m.poisson_transport_residual_abs << "\n";
    out << "poisson_transport_scale=" << m.poisson_transport_scale << "\n";
    out << "poisson_transport_pass=" << (m.poisson_transport_pass ? 1 : 0) << "\n";
    out << "force_work_residual_abs=" << m.force_work_residual_abs << "\n";
    out << "force_work_scale=" << m.force_work_scale << "\n";
    out << "force_work_pass=" << (m.force_work_pass ? 1 : 0) << "\n";
    out << "current_pair_residual_abs=" << m.current_pair_residual_abs << "\n";
    out << "current_pair_scale=" << m.current_pair_scale << "\n";
    out << "current_pair_pass=" << (m.current_pair_pass ? 1 : 0) << "\n";
    out << "poisson_identity_crosscheck="
        << m.poisson_identity_crosscheck << "\n";
    out << "conversion_transfer_residual_abs="
        << m.conversion_transfer_residual_abs << "\n";
    out << "conversion_transfer_scale=" << m.conversion_transfer_scale << "\n";
    out << "boundary_residual_abs=" << m.boundary_residual_abs << "\n";
    out << "boundary_scale=" << m.boundary_scale << "\n";
    out << "full_residual_abs=" << m.full_residual_abs << "\n";
    out << "reconstructed_residual_abs=" << m.reconstructed_residual_abs << "\n";
    out << "reconstruction_mismatch=" << m.reconstruction_mismatch << "\n";
    out << "reconstruction_tolerance=" << m.reconstruction_tolerance << "\n";
    out << "reconstruction_pass=" << (m.reconstruction_pass ? 1 : 0) << "\n";
    out << "bulk_work_plus=" << m.bulk_work_plus << "\n";
    out << "bulk_work_minus=" << m.bulk_work_minus << "\n";
    out << "tail_cell_work_mismatch=" << m.tail_cell_work_mismatch << "\n";
    out << "beam_cell_work_mismatch=" << m.beam_cell_work_mismatch << "\n";
    out << "manufactured_face_work=" << m.manufactured_face_work << "\n";
    out << "manufactured_cell_work=" << m.manufactured_cell_work << "\n";
    out << "manufactured_identity_error="
        << m.manufactured_identity_error << "\n";
    out << "manufactured_identity_pass="
        << (m.manufactured_identity_pass ? 1 : 0) << "\n";
    out << "manufactured_endpoint_weight_pass="
        << (m.manufactured_endpoint_weight_pass ? 1 : 0) << "\n";
    out << "read_only_bitwise_equal=" << (m.read_only_bitwise_equal ? 1 : 0) << "\n";
    out << "negative_test_effective=" << m.negative_test_effective << "\n";
    out << "injected_total_count=" << m.injected_total_count << "\n";
    out << "injected_detected_count=" << m.injected_detected_count << "\n";
    for (int i = 0; i < 6; ++i) {
        out << "fault_" << (i + 1) << "_detected="
            << m.injected_detected[i] << "\n";
        out << "fault_" << (i + 1) << "_first_bad_index="
            << m.injected_first_bad[i] << "\n";
    }
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "vpfp_field_particle_pairing_test must run with exactly 1 "
                     "rank; use the MPI variant for 2/5 ranks.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: vpfp_field_particle_pairing_test --case "
                     "all|injected-faults|fault-1..fault-6 "
                     "[--result <path>]\n";
    }

    PairingMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "all") {
        run_bulk_continuity_and_reconstruction(m);
        run_manufactured_gstar_identity(m);
        run_bulk_force(m);
        run_pic_kick_work(m);
        run_read_only_equivalence(m);
        const double work_scale = std::max(1.0,
            std::max(std::fabs(m.bulk_work_plus), std::fabs(m.bulk_work_minus)));
        pass = m.continuity_pass &&
               m.poisson_transport_pass &&
               m.force_work_pass &&
               m.current_pair_pass &&
               m.reconstruction_pass &&
               m.manufactured_identity_pass &&
               m.manufactured_endpoint_weight_pass &&
               // Electron charge is negative: +E decelerates the selected
               // positive-u population, while -E accelerates it.
               m.bulk_work_plus < 0.0 && m.bulk_work_minus > 0.0 &&
               m.tail_cell_work_mismatch <=
                   machine_scaled_tolerance(work_scale, work_scale) &&
               m.beam_cell_work_mismatch <=
                   machine_scaled_tolerance(work_scale, work_scale) &&
               m.read_only_bitwise_equal;
    } else if (ok && args.test_case == "injected-faults") {
        const int fault_ids[6] = { 1, 2, 3, 4, 5, 6 };
        for (int i = 0; i < 6; ++i) run_injected_fault(m, fault_ids[i]);
        m.negative_test_effective =
            (m.injected_detected_count == m.injected_total_count) ? 1 : 0;
        pass = m.negative_test_effective == 1;
    } else if (ok && args.test_case.size() == 7 &&
               args.test_case.compare(0, 6, "fault-") == 0 &&
               args.test_case[6] >= '1' && args.test_case[6] <= '6') {
        const int fault_id = args.test_case[6] - '0';
        pass = run_injected_fault(m, fault_id);
        m.negative_test_effective = pass ? 1 : 0;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << std::setprecision(17)
              << "continuity_residual_abs=" << m.continuity_residual_abs
              << " continuity_pass=" << (m.continuity_pass ? 1 : 0)
              << " reconstruction_mismatch=" << m.reconstruction_mismatch
              << " reconstruction_pass=" << (m.reconstruction_pass ? 1 : 0)
              << " bulk_work_plus=" << m.bulk_work_plus
              << " bulk_work_minus=" << m.bulk_work_minus
              << " tail_cell_work_mismatch=" << m.tail_cell_work_mismatch
              << " beam_cell_work_mismatch=" << m.beam_cell_work_mismatch
              << " manufactured_identity_error="
              << m.manufactured_identity_error
              << " manufactured_identity_pass="
              << (m.manufactured_identity_pass ? 1 : 0)
              << " manufactured_endpoint_weight_pass="
              << (m.manufactured_endpoint_weight_pass ? 1 : 0)
              << " read_only_bitwise_equal=" << (m.read_only_bitwise_equal ? 1 : 0)
              << " negative_test_effective=" << m.negative_test_effective << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
