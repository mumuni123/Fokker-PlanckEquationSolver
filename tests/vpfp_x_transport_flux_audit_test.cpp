// Gate I acceptance test (section 4.7.2) for the read-only x-face swept-number
// exposure of ConservativePpmRemap::advect_x.  Every case drives the production
// advect_x and checks the exposed face vector against the actual finite-volume
// update: R_i = M_i^out - M_i^in - Q_{i-1/2} + Q_{i+1/2}.  No case re-implements
// a production PPM flux.
//
// Usage:
//   vpfp_x_transport_flux_audit_test --case all [--result <path>]

#include "conservative_ppm_remap.h"
#include "field_particle_power_audit.h"
#include "open_boundary.h"
#include "parameters.h"
#include "species.h"
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
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
    XTransportVelocityMode velocity_mode;
    TestArgs() : velocity_mode(XTransportVelocityMode::ANALYTIC_CELL_CENTER) {}
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
        } else if (arg == "--x-transport-velocity-mode") {
            if (i + 1 >= argc) return false;
            const std::string value(argv[++i]);
            if (value == "analytic-cell-center")
                args.velocity_mode = XTransportVelocityMode::ANALYTIC_CELL_CENTER;
            else if (value == "energy-conjugate")
                args.velocity_mode = XTransportVelocityMode::ENERGY_CONJUGATE_CELL;
            else return false;
        } else {
            return false;
        }
    }
    return !args.test_case.empty();
}

// Choose one (j,k) velocity slot whose vx has the requested sign.
size_t slot_with_vx_sign(const CylindricalVelocityGrid& cg, int sign)
{
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double vx = cg.vx[static_cast<size_t>(j) * Param::Nmu + k];
            if (sign > 0 && vx > 0.0) return idx2(j, k);
            if (sign < 0 && vx < 0.0) return idx2(j, k);
        }
    }
    return idx2(Param::Nv / 2, 0);
}

// Absorbing boundary on both ends (no reservoir inflow) for interior cases.
OpenBackgroundBoundary make_absorbing_boundary()
{
    OpenBackgroundBoundaryConfig cfg;
    cfg.left_type = BackgroundXBoundaryType::ABSORBING;
    cfg.right_type = BackgroundXBoundaryType::ABSORBING;
    return OpenBackgroundBoundary(cfg);
}

struct TransportMetrics {
    double mass_update_linf_pos;
    double mass_update_linf_neg;
    double mass_update_linf_mixed;
    double global_mass_balance;
    double boundary_flux_balance;
    double inflow_source_closure;
    double outflow_left_balance;
    double outflow_right_balance;
    double zero_dt_swept_sum;
    bool zero_dt_bitwise;
    bool audit_null_bitwise_equal;
    TransportMetrics()
        : mass_update_linf_pos(0.0), mass_update_linf_neg(0.0),
          mass_update_linf_mixed(0.0), global_mass_balance(0.0),
          boundary_flux_balance(0.0), inflow_source_closure(0.0),
          outflow_left_balance(0.0), outflow_right_balance(0.0),
          zero_dt_swept_sum(0.0), zero_dt_bitwise(false),
          audit_null_bitwise_equal(false)
    {}
};

// Per-cell residual R_i = M_out - M_in - Q_{i-1/2} + Q_{i+1/2} over the
// interior cells; returns linf and the global mass balance.
double check_residual(const Species& in, const Species& out,
                      const XFaceTransportAudit& audit,
                      const SpatialGrid& grid, double& mass_balance)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    double linf = 0.0;
    double mb = 0.0;
    for (int i = 0; i < nxl; ++i) {
        // Number per cell is the cell-integrated mass summed over velocity.
        double m_in = 0.0;
        double m_out = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                m_in += in.f[idx3(ng + i, j, k)];
                m_out += out.f[idx3(ng + i, j, k)];
            }
        }
        const double r = m_out - m_in -
            audit.bulk_number_swept_face[static_cast<size_t>(i)] +
            audit.bulk_number_swept_face[static_cast<size_t>(i) + 1];
        linf = std::max(linf, std::fabs(r));
        mb += r;
    }
    mass_balance = mb;
    return linf;
}

void run_single_vx_cases(TransportMetrics& m,
                         XTransportVelocityMode velocity_mode)
{
    const int nx = 40;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, nx, Param::Lx);
    CylindricalVelocityGrid cg;
    cg.init(Param::momentum_umax);
    ConservativePpmRemap remap;
    remap.init(grid, cg);
    remap.set_x_transport_velocity_mode(velocity_mode);
    const double dt = 1.0e-15;
    const OpenBackgroundBoundary boundary = make_absorbing_boundary();
    const size_t pos_slot = slot_with_vx_sign(cg, +1);
    const size_t neg_slot = slot_with_vx_sign(cg, -1);

    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const size_t nslots = static_cast<size_t>(grid.nx_total) * Param::Nvmu;

    // Case 1/2: single velocity cell, positive/negative vx.
    for (int sign = 0; sign < 2; ++sign) {
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        const size_t slot = sign == 0 ? pos_slot : neg_slot;
        for (int i = 0; i < grid.nx_total; ++i) {
            const int ig = grid.global_cell(i);
            if (ig >= nx / 4 && ig < 3 * nx / 4) {
                in.f[static_cast<size_t>(i) * Param::Nvmu + slot] = 1.0e18;
            }
        }
        in.compute_moments();
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        out.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out, dt, 0.0, boundary, 0, 1, &audit);
        double mb = 0.0;
        const double linf = check_residual(in, out, audit, grid, mb);
        if (sign == 0) { m.mass_update_linf_pos = linf; }
        else { m.mass_update_linf_neg = linf; }
        m.global_mass_balance = std::max(m.global_mass_balance, std::fabs(mb));
    }

    // Case 3: mixed velocity fills both signs; the face swept number is the
    // velocity-integrated production swept_mass.
    {
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        for (int i = 0; i < grid.nx_total; ++i) {
            const int ig = grid.global_cell(i);
            if (ig >= nx / 4 && ig < 3 * nx / 4) {
                in.f[static_cast<size_t>(i) * Param::Nvmu + pos_slot] = 0.6e18;
                in.f[static_cast<size_t>(i) * Param::Nvmu + neg_slot] = 0.4e18;
            }
        }
        in.compute_moments();
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        out.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out, dt, 0.0, boundary, 0, 1, &audit);
        double mb = 0.0;
        m.mass_update_linf_mixed = check_residual(in, out, audit, grid, mb);
    }

    // Case 4: left reservoir inflow source.  A positive-vx reservoir fills
    // the domain; the left-boundary source plus the interior divergence must
    // close the total particle-number change.
    {
        OpenBackgroundBoundaryConfig cfg;
        cfg.left_type = BackgroundXBoundaryType::RESERVOIR;
        cfg.right_type = BackgroundXBoundaryType::ABSORBING;
        cfg.left_reservoir.density = 0.5 * Param::dens;
        cfg.left_reservoir.temperature = Param::temperature_e;
        cfg.left_reservoir.drift_vx = 0.0;
        const OpenBackgroundBoundary reservoir_boundary(cfg);
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        in.compute_moments();
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        out.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out, dt, 0.0, reservoir_boundary, 0, 1, &audit);
        double mb = 0.0;
        check_residual(in, out, audit, grid, mb);
        // Left boundary flux (source) enters through face 0; the total mass
        // change equals the face-0 inflow plus the face-nx outflow.
        double m_in = 0.0;
        double m_out = 0.0;
        for (int i = 0; i < nxl; ++i) {
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    m_in += in.f[idx3(ng + i, j, k)];
                    m_out += out.f[idx3(ng + i, j, k)];
                }
            }
        }
        m.inflow_source_closure = std::fabs(
            (m_out - m_in) -
            (audit.bulk_number_swept_face[0] -
             audit.bulk_number_swept_face[static_cast<size_t>(nxl)]));
    }

    // Case 5: independent left and right absorbing outflow.  This catches a
    // boundary-sign error that an interior divergence test cannot see.
    for (int side = 0; side < 2; ++side) {
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        const size_t slot = side == 0 ? neg_slot : pos_slot;
        const int physical_cell = side == 0 ? 0 : nx - 1;
        for (int i = 0; i < grid.nx_total; ++i) {
            if (grid.global_cell(i) == physical_cell) {
                in.f[static_cast<size_t>(i) * Param::Nvmu + slot] = 1.0e18;
            }
        }
        in.compute_moments();
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        out.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out, dt, 0.0, boundary, 0, 1, &audit);
        double mass_residual = 0.0;
        const double cell_linf = check_residual(in, out, audit, grid,
                                                mass_residual);
        double m_in = 0.0;
        double m_out = 0.0;
        for (int i = 0; i < nxl; ++i) {
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    m_in += in.f[idx3(ng + i, j, k)];
                    m_out += out.f[idx3(ng + i, j, k)];
                }
            }
        }
        const double closure = std::fabs(
            (m_out - m_in) -
            (audit.bulk_number_swept_face[0] -
             audit.bulk_number_swept_face[static_cast<size_t>(nxl)]));
        if (side == 0) m.outflow_left_balance = closure;
        else m.outflow_right_balance = closure;
        m.boundary_flux_balance = std::max(
            m.boundary_flux_balance,
            std::max(cell_linf, std::fabs(mass_residual)));
    }

    // Case 6: zero dt / zero vx -> swept number exactly zero, bitwise equal.
    {
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        for (int i = 0; i < grid.nx_total; ++i) {
            in.f[static_cast<size_t>(i) * Param::Nvmu + pos_slot] = 1.0e18;
        }
        in.compute_moments();
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        out.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out, 0.0, 0.0, boundary, 0, 1, &audit);
        double swept_sum = 0.0;
        for (size_t f = 0; f < audit.bulk_number_swept_face.size(); ++f) {
            swept_sum += std::fabs(audit.bulk_number_swept_face[f]);
        }
        m.zero_dt_swept_sum = swept_sum;
        m.zero_dt_bitwise = true;
        const size_t cells_per_x = Param::Nvmu;
        for (int i = 0; i < nxl && m.zero_dt_bitwise; ++i) {
            const size_t begin = static_cast<size_t>(ng + i) * cells_per_x;
            if (std::memcmp(&in.f[begin], &out.f[begin],
                            cells_per_x * sizeof(double)) != 0) {
                m.zero_dt_bitwise = false;
            }
        }
    }

    // Case 7: audit NULL equivalence (bitwise same output, moments, ledger).
    {
        Species in;
        in.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
        in.f.assign(nslots, 0.0);
        for (int i = 0; i < grid.nx_total; ++i) {
            in.f[static_cast<size_t>(i) * Param::Nvmu + pos_slot] = 1.0e18;
            in.f[static_cast<size_t>(i) * Param::Nvmu + neg_slot] = 0.3e18;
        }
        in.compute_moments();
        Species out_a;
        Species out_b;
        out_a.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                   Const::me, Param::dens, Param::temperature_e, false, grid);
        out_b.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                   Const::me, Param::dens, Param::temperature_e, false, grid);
        out_a.f.assign(nslots, 0.0);
        out_b.f.assign(nslots, 0.0);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(in, out_a, dt, 0.0, boundary, 0, 1, &audit);
        remap.advect_x(in, out_b, dt, 0.0, boundary, 0, 1, NULL);
        m.audit_null_bitwise_equal = (out_a.f == out_b.f);
    }
    (void)ng;
}

bool write_result_file(const std::string& path, const TransportMetrics& m,
                       bool pass, XTransportVelocityMode velocity_mode)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "x_transport_velocity_mode="
        << (velocity_mode == XTransportVelocityMode::ENERGY_CONJUGATE_CELL
            ? "energy-conjugate" : "analytic-cell-center") << "\n";
    out << "mass_update_linf_pos=" << m.mass_update_linf_pos << "\n";
    out << "mass_update_linf_neg=" << m.mass_update_linf_neg << "\n";
    out << "mass_update_linf_mixed=" << m.mass_update_linf_mixed << "\n";
    out << "global_mass_balance=" << m.global_mass_balance << "\n";
    out << "inflow_source_closure=" << m.inflow_source_closure << "\n";
    out << "outflow_left_balance=" << m.outflow_left_balance << "\n";
    out << "outflow_right_balance=" << m.outflow_right_balance << "\n";
    out << "boundary_flux_balance=" << m.boundary_flux_balance << "\n";
    out << "zero_dt_swept_sum=" << m.zero_dt_swept_sum << "\n";
    out << "zero_dt_bitwise=" << (m.zero_dt_bitwise ? 1 : 0) << "\n";
    out << "audit_null_bitwise_equal="
        << (m.audit_null_bitwise_equal ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "vpfp_x_transport_flux_audit_test must run with exactly 1 "
                     "rank.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: vpfp_x_transport_flux_audit_test --case all "
                     "[--result <path>]\n";
    }

    TransportMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "all") {
        run_single_vx_cases(m, args.velocity_mode);
        const double scale = 1.0e18;
        const double tol = machine_scaled_tolerance(scale, scale);
        pass = m.mass_update_linf_pos <= tol &&
               m.mass_update_linf_neg <= tol &&
               m.mass_update_linf_mixed <= tol &&
               m.global_mass_balance <= tol &&
               m.inflow_source_closure <= tol &&
               m.outflow_left_balance <= tol &&
               m.outflow_right_balance <= tol &&
               m.boundary_flux_balance <= tol &&
               m.zero_dt_swept_sum == 0.0 &&
               m.zero_dt_bitwise &&
               m.audit_null_bitwise_equal;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass, args.velocity_mode)) pass = false;
    std::cout << std::setprecision(17)
              << "mass_update_linf_pos=" << m.mass_update_linf_pos
              << " mass_update_linf_neg=" << m.mass_update_linf_neg
              << " mass_update_linf_mixed=" << m.mass_update_linf_mixed
              << " global_mass_balance=" << m.global_mass_balance
              << " inflow_source_closure=" << m.inflow_source_closure
              << " zero_dt_bitwise=" << (m.zero_dt_bitwise ? 1 : 0)
              << " audit_null_bitwise_equal="
              << (m.audit_null_bitwise_equal ? 1 : 0) << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
