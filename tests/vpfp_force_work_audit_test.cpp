// Gate C acceptance test (section 7.8) for the read-only discrete-force-work
// recording added to the production operators.  Every case drives the actual
// production implementation (ConservativePpmRemap::advect_u_parallel,
// BackgroundTailPIC::kick, BeamPIC::last_field_work) and compares the recorded
// work against an independent direct kinetic-energy difference.  No case
// re-implements a production flux formula.
//
// Usage:
//   vpfp_force_work_audit_test --case all [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "beam_pic.h"
#include "bulk_tail_flux_parcel.h"
#include "conservative_ppm_remap.h"
#include "grid.h"
#include "maxwell.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
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

// Machine-precision-scaled tolerance (never an arbitrary 1e-6).
double roundoff_tolerance(double scale)
{
    return 4096.0 * std::numeric_limits<double>::epsilon() *
           std::max(1.0, scale);
}

CylindricalVelocityGrid make_velocity_grid(int nv)
{
    CylindricalVelocityGrid grid;
    grid.init(Param::momentum_umax, nv, Param::Nmu);
    return grid;
}

Species make_packed_species(const SpatialGrid& grid, int nv)
{
    Species sp;
    sp.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
            Const::me, Param::dens, Param::temperature_e, false, grid);
    sp.f.assign(static_cast<size_t>(grid.nx_total) *
                static_cast<size_t>(nv) * Param::Nmu, 0.0);
    return sp;
}

size_t packed_index(int ix, int j, int k, int nv)
{
    return (static_cast<size_t>(ix) * nv + static_cast<size_t>(j)) *
               Param::Nmu + static_cast<size_t>(k);
}

// fbar is the point value of the phase-space density; fill_profile multiplies
// by the velocity-cell phase volume (du * ring) and dx to obtain the
// cell-integrated conservative mass stored in Species::f.
void fill_profile(Species& sp, const SpatialGrid& grid,
                  const CylindricalVelocityGrid& cg, int nv,
                  double amp, double center, double sigma)
{
    const int nmu = Param::Nmu;
    for (int ix = 0; ix < grid.nx_total; ++ix) {
        for (int j = 0; j < nv; ++j) {
            const double u = cg.upar_cells[static_cast<size_t>(j)];
            const double du = cg.upar_widths[static_cast<size_t>(j)];
            for (int k = 0; k < nmu; ++k) {
                const double up = cg.uperp_cells[static_cast<size_t>(k)];
                const double ring = cg.uperp_ring_areas[static_cast<size_t>(k)];
                const double fbar = amp * std::exp(
                    -((u - center) * (u - center) + up * up) /
                    (2.0 * sigma * sigma));
                sp.f[packed_index(ix, j, k, nv)] =
                    fbar * du * ring * grid.dx;
            }
        }
    }
}

// Uniform block over all velocity cells so the tail-owned high-energy cells
// carry mass (used by the flux-interface case).
void fill_uniform(Species& sp, const SpatialGrid& grid,
                  const CylindricalVelocityGrid& cg, int nv, double amp)
{
    const int nmu = Param::Nmu;
    for (int ix = 0; ix < grid.nx_total; ++ix) {
        for (int j = 0; j < nv; ++j) {
            const double du = cg.upar_widths[static_cast<size_t>(j)];
            for (int k = 0; k < nmu; ++k) {
                const double ring = cg.uperp_ring_areas[static_cast<size_t>(k)];
                sp.f[packed_index(ix, j, k, nv)] = amp * du * ring * grid.dx;
            }
        }
    }
}

// Block touching one u_parallel boundary, for the velocity-outflow sign cases.
void fill_boundary_block(Species& sp, const SpatialGrid& grid,
                         const CylindricalVelocityGrid& cg, int nv,
                         double lo, double hi, double amp)
{
    const int nmu = Param::Nmu;
    for (int ix = 0; ix < grid.nx_total; ++ix) {
        for (int j = 0; j < nv; ++j) {
            const double u = cg.upar_cells[static_cast<size_t>(j)];
            if (u < lo || u > hi) continue;
            const double du = cg.upar_widths[static_cast<size_t>(j)];
            for (int k = 0; k < nmu; ++k) {
                const double ring = cg.uperp_ring_areas[static_cast<size_t>(k)];
                sp.f[packed_index(ix, j, k, nv)] = amp * du * ring * grid.dx;
            }
        }
    }
}

double total_kinetic_energy(const Species& sp, const SpatialGrid& grid,
                            const CylindricalVelocityGrid& cg, int nv)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    double total = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                total += sp.f[packed_index(ng + ix, j, k, nv)] *
                         cg.kinetic_energy[static_cast<size_t>(j) * Param::Nmu + k];
            }
        }
    }
    return total;
}

EMFields make_constant_field(const SpatialGrid& grid, double e_x)
{
    EMFields fields;
    fields.init(grid);
    std::fill(fields.Ex.begin(), fields.Ex.end(), e_x);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), e_x);
    return fields;
}

// a_u = q E_x / (m c).  For background electrons q = -qe, so a positive E_x
// drives a_u negative (mass shifts toward -u).
double e_x_from_a_u(double a_u)
{
    return a_u * Const::me * Const::c / (-Const::qe);
}

struct IdentityMetrics {
    double residual_plus;
    double residual_minus;
    double scale_plus;
    double scale_minus;
    double left_boundary_energy;
    double right_boundary_energy;
    double interface_removed;
    double interface_residual;
    double interface_expected_removed;
    double tail_work_plus;
    double tail_work_minus;
    double tail_work_mismatch;
    double beam_work_mismatch;
    double beam_work_scale;
    double beam_work_zero_before_push;
    bool state_bitwise_equal_with_null;
    IdentityMetrics()
        : residual_plus(0.0), residual_minus(0.0), scale_plus(1.0),
          scale_minus(1.0), left_boundary_energy(0.0),
          right_boundary_energy(0.0), interface_removed(0.0),
          interface_residual(0.0), interface_expected_removed(0.0),
          tail_work_plus(0.0), tail_work_minus(0.0),
          tail_work_mismatch(0.0), beam_work_mismatch(0.0),
          beam_work_scale(1.0), beam_work_zero_before_push(0.0),
          state_bitwise_equal_with_null(false)
    {}
};

// Cases 1-3: non-uniform u_parallel grid, constant positive/negative field,
// no velocity-boundary outflow and no tail interface -> the discrete telescope
// identity R = dK - W_internal - W_left - W_right + K_interface_removed must
// close to machine precision.  A centered Gaussian never reaches a velocity
// boundary, so both boundary energies and the interface-removed energy are 0.
void run_bulk_identity(IdentityMetrics& m)
{
    const int nv = Param::Nv;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 40.0 * Const::micro);
    CylindricalVelocityGrid cg = make_velocity_grid(nv);
    ConservativePpmRemap remap;
    remap.init(grid, cg);
    const double amp = 1.0e20;
    const double sigma = 0.06;
    const double dt = 1.0e-15;
    const double a_u = 0.03 / dt;

    for (int sign = 0; sign < 2; ++sign) {
        const double a = sign == 0 ? a_u : -a_u;
        Species sp = make_packed_species(grid, nv);
        fill_profile(sp, grid, cg, nv, amp, 0.0, sigma);
        Species out = make_packed_species(grid, nv);
        const double ke_before = total_kinetic_energy(sp, grid, cg, nv);
        const RemapDiagnostics diag = remap.advect_u_parallel(
            sp, out, make_constant_field(grid, e_x_from_a_u(a)), dt, 0.0);
        const double ke_after = total_kinetic_energy(out, grid, cg, nv);
        const double dk = ke_after - ke_before;
        const double scale = std::max(
            1.0, std::max(std::fabs(dk),
                          std::max(std::fabs(diag.upar_internal_face_energy_transfer),
                                   std::max(
                                       std::fabs(diag.upar_left_velocity_boundary_energy),
                                       std::fabs(diag.upar_right_velocity_boundary_energy)))));
        if (sign == 0) {
            m.residual_plus = diag.upar_discrete_energy_identity_residual;
            m.scale_plus = scale;
        } else {
            m.residual_minus = diag.upar_discrete_energy_identity_residual;
            m.scale_minus = scale;
        }
    }
}

// Case 4: velocity-boundary outflow signs.  A block touching the left
// boundary pushed further left gives W_left < 0; a block touching the right
// boundary pushed further right gives W_right < 0.
void run_boundary_signs(IdentityMetrics& m)
{
    const int nv = Param::Nv;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 100, 40.0 * Const::micro);
    CylindricalVelocityGrid cg = make_velocity_grid(nv);
    ConservativePpmRemap remap;
    remap.init(grid, cg);
    const double amp = 1.0e20;
    const double umax = Param::momentum_upar_core_max;
    const double dt = 1.0e-15;

    {
        // Left outflow: a_u < 0 (positive E_x for electrons).
        Species sp = make_packed_species(grid, nv);
        fill_boundary_block(sp, grid, cg, nv, -umax, -umax + 2.0, amp);
        Species out = make_packed_species(grid, nv);
        const RemapDiagnostics diag = remap.advect_u_parallel(
            sp, out, make_constant_field(grid, e_x_from_a_u(-1.0 / dt)),
            dt, 0.0);
        m.left_boundary_energy = diag.upar_left_velocity_boundary_energy;
    }
    {
        // Right outflow: a_u > 0 (negative E_x for electrons).
        Species sp = make_packed_species(grid, nv);
        fill_boundary_block(sp, grid, cg, nv, umax - 2.0, umax, amp);
        Species out = make_packed_species(grid, nv);
        const RemapDiagnostics diag = remap.advect_u_parallel(
            sp, out, make_constant_field(grid, e_x_from_a_u(1.0 / dt)),
            dt, 0.0);
        m.right_boundary_energy = diag.upar_right_velocity_boundary_energy;
    }
}

// Case 5: flux-interface sink clears tail-owned Eulerian cells.  With a
// uniform profile and zero field the swept masses vanish, so the cleared
// energy is exactly sum_tail K_j * M_j and must be recorded exactly once:
// the identity closes with a single +K_interface_removed term.
void run_flux_interface(IdentityMetrics& m)
{
    const int nv = Param::Nv;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 40, 40.0 * Const::micro);
    CylindricalVelocityGrid cg = make_velocity_grid(nv);
    ConservativePpmRemap remap;
    remap.init(grid, cg);
    HybridVelocityPartition partition;
    partition.init(cg, 6.0, 1.0, 1, 1);
    const double amp = 1.0e18;

    Species sp = make_packed_species(grid, nv);
    fill_uniform(sp, grid, cg, nv, amp);
    Species out = make_packed_species(grid, nv);
    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    const RemapDiagnostics diag = remap.advect_u_parallel(
        sp, out, make_constant_field(grid, 0.0), 1.0e-15, 0.0,
        &partition, &batch, 4);

    // Independent reference: the energy cleared from tail-owned cells.
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int nmu = Param::Nmu;
    double expected = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (!partition.is_tail_owned(j, k)) continue;
                expected += sp.f[packed_index(ng + ix, j, k, nv)] *
                            cg.kinetic_energy[static_cast<size_t>(j) * nmu + k];
            }
        }
    }
    m.interface_removed = diag.upar_interface_energy_removed;
    m.interface_expected_removed = expected;
    m.interface_residual = diag.upar_discrete_energy_identity_residual;
}

double tail_particle_ke(const BackgroundTailParticle& p)
{
    const double gamma = std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy +
                                   p.uz * p.uz);
    return (gamma - 1.0) * Const::me * Const::c * Const::c;
}

// Case 6: Tail single-particle positive/negative-field kick.  The recorded
// rank-local work must equal weight * [K(u_after) - K(u_before)] directly.
void run_tail_kick(IdentityMetrics& m)
{
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 100, 40.0 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    const double weight = 1.0e14;
    const double dt = 1.0e-15;

    for (int sign = 0; sign < 2; ++sign) {
        BackgroundTailPIC t;
        t.init(grid);
        BackgroundTailParticle p = {};
        p.x = 20.0 * Const::micro;
        p.ux = 0.5;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = weight;
        p.id = 0;
        t.particles.push_back(p);
        const double ke_before = tail_particle_ke(t.particles[0]);
        // Electron: positive E_x decreases ux (work < 0), negative E_x
        // increases ux (work > 0).
        const double e_x = sign == 0 ? -1.0e12 : 1.0e12;
        double work = 0.0;
        t.kick(grid, make_constant_field(grid, e_x), dt, 0, 1, &work);
        const double ke_after = tail_particle_ke(t.particles[0]);
        const double expected = weight * (ke_after - ke_before);
        const double mismatch = std::fabs(work - expected);
        if (sign == 0) {
            m.tail_work_plus = work;
            m.tail_work_mismatch = std::max(m.tail_work_mismatch, mismatch);
        } else {
            m.tail_work_minus = work;
            m.tail_work_mismatch = std::max(m.tail_work_mismatch, mismatch);
        }
    }
}

// Case 8: with the nullable work pointer absent the kick must produce a
// bitwise-identical state (no extra sqrt, scan or state perturbation).
void run_diagnostic_off(IdentityMetrics& m)
{
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 100, 40.0 * Const::micro);
    BackgroundTailPIC a;
    BackgroundTailPIC b;
    a.init(grid);
    b.init(grid);
    for (int i = 0; i < 3; ++i) {
        BackgroundTailParticle p = {};
        p.x = (10.0 + 10.0 * i) * Const::micro;
        p.ux = 0.1 * (i + 1);
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0e14;
        p.id = static_cast<std::uint64_t>(i);
        a.particles.push_back(p);
        b.particles.push_back(p);
    }
    EMFields fields = make_constant_field(grid, 5.0e11);
    double work = 0.0;
    a.kick(grid, fields, 1.0e-15, 0, 1, NULL);
    b.kick(grid, fields, 1.0e-15, 0, 1, &work);
    m.state_bitwise_equal_with_null = true;
    if (a.particles.size() != b.particles.size()) {
        m.state_bitwise_equal_with_null = false;
        return;
    }
    for (size_t i = 0; i < a.particles.size(); ++i) {
        const BackgroundTailParticle& pa = a.particles[i];
        const BackgroundTailParticle& pb = b.particles[i];
        if (std::memcmp(&pa, &pb, sizeof(pa)) != 0) {
            m.state_bitwise_equal_with_null = false;
            return;
        }
    }
}

// Case 7: Beam last_field_work() must equal the direct kinetic-energy
// difference over the midpoint kick, and must read zero before the push.
void run_beam_work(IdentityMetrics& m)
{
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 2000, 200.0 * Const::micro);
    BeamPIC beam;
    beam.init(grid);
    BeamParticle particle;
    particle.x = 40.0 * Const::micro;
    particle.px = 0.5 * Const::me * Const::c;
    particle.weight = 1.0e12;
    beam.particles.push_back(particle);
    beam.begin_step(grid, 2.5e-15);
    m.beam_work_zero_before_push = beam.last_field_work();

    const double ke_before = beam.total_kinetic_energy();
    BeamPIC work;
    work.init(grid);
    work.begin_step(grid, 2.5e-15);
    const double e_field = -1.0e11;
    BeamInjectionSchedule empty;
    beam.predict_to_midpoint(empty, grid, make_constant_field(grid, e_field),
                             0.0, 2.5e-15, 0, 1, work);
    work.finish_from_midpoint(empty, grid, make_constant_field(grid, e_field),
                              0.0, 2.5e-15, 0, 1);
    const double ke_after = work.total_kinetic_energy();
    const double expected = ke_after - ke_before;
    m.beam_work_mismatch = std::fabs(work.last_field_work() - expected);
    m.beam_work_scale = std::max(1.0, std::fabs(expected));
}

bool write_result_file(const std::string& path, const IdentityMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "bulk_identity_residual_plus=" << m.residual_plus << "\n";
    out << "bulk_identity_scale_plus=" << m.scale_plus << "\n";
    out << "bulk_identity_residual_minus=" << m.residual_minus << "\n";
    out << "bulk_identity_scale_minus=" << m.scale_minus << "\n";
    out << "left_velocity_boundary_energy=" << m.left_boundary_energy << "\n";
    out << "right_velocity_boundary_energy=" << m.right_boundary_energy << "\n";
    out << "interface_removed=" << m.interface_removed << "\n";
    out << "interface_expected_removed=" << m.interface_expected_removed << "\n";
    out << "interface_identity_residual=" << m.interface_residual << "\n";
    out << "tail_work_plus=" << m.tail_work_plus << "\n";
    out << "tail_work_minus=" << m.tail_work_minus << "\n";
    out << "tail_work_mismatch=" << m.tail_work_mismatch << "\n";
    out << "beam_work_mismatch=" << m.beam_work_mismatch << "\n";
    out << "beam_work_zero_before_push=" << m.beam_work_zero_before_push << "\n";
    out << "state_bitwise_equal_with_null="
        << (m.state_bitwise_equal_with_null ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "vpfp_force_work_audit_test must run with exactly 1 "
                     "rank; use plain ./build/vpfp_force_work_audit_test "
                     "(no yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    const bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: vpfp_force_work_audit_test --case all "
                     "[--result <path>]\n";
    }

    IdentityMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "all") {
        run_bulk_identity(m);
        run_boundary_signs(m);
        run_flux_interface(m);
        run_tail_kick(m);
        run_beam_work(m);
        run_diagnostic_off(m);

        // Machine-precision-scaled telescope identity gates.
        pass = std::fabs(m.residual_plus) <= roundoff_tolerance(m.scale_plus) &&
               std::fabs(m.residual_minus) <= roundoff_tolerance(m.scale_minus) &&
               // Velocity-boundary outflow signs: energy leaves the bulk.
               m.left_boundary_energy < 0.0 &&
               m.right_boundary_energy < 0.0 &&
               // Interface removal recorded exactly once: identity closes and
               // the recorded value equals the independently summed cleared
               // tail-owned energy.
               m.interface_removed > 0.0 &&
               std::fabs(m.interface_removed - m.interface_expected_removed) <=
                   roundoff_tolerance(m.interface_expected_removed) &&
               std::fabs(m.interface_residual) <=
                   roundoff_tolerance(m.interface_expected_removed) &&
               // Tail kick work equals the direct KE difference.
               m.tail_work_plus > 0.0 &&
               m.tail_work_minus < 0.0 &&
               m.tail_work_mismatch <= roundoff_tolerance(
                   std::max(1.0, std::max(std::fabs(m.tail_work_plus),
                                          std::fabs(m.tail_work_minus)))) &&
               // Beam field work equals the direct KE difference.
               m.beam_work_zero_before_push == 0.0 &&
               m.beam_work_mismatch <= roundoff_tolerance(m.beam_work_scale) &&
               // Null diagnostic pointer leaves the pushed state bitwise
               // identical to the recording path.
               m.state_bitwise_equal_with_null;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << std::setprecision(17)
              << "bulk_identity_residual_plus=" << m.residual_plus
              << " bulk_identity_residual_minus=" << m.residual_minus
              << " left_velocity_boundary_energy=" << m.left_boundary_energy
              << " right_velocity_boundary_energy=" << m.right_boundary_energy
              << " interface_removed=" << m.interface_removed
              << " interface_expected_removed=" << m.interface_expected_removed
              << " tail_work_plus=" << m.tail_work_plus
              << " tail_work_minus=" << m.tail_work_minus
              << " tail_work_mismatch=" << m.tail_work_mismatch
              << " beam_work_mismatch=" << m.beam_work_mismatch
              << " state_bitwise_equal_with_null="
              << (m.state_bitwise_equal_with_null ? 1 : 0) << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
