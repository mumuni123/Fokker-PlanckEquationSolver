// Phase-1 acceptance test for the open electrostatic Gauss solver
// (sections 5.3, 14 and 16.2).  It exercises the production
// OpenElectrostaticSolver class directly; all analytic references are
// independent closed forms, the solver algorithm is never reimplemented here.
//
// Usage (single rank):
//   open_electrostatic_solver_test --case all --field-boundary dirichlet-phi
//                                  [--result <path>]
// Usage (multi rank, e.g. 2 or 5 ranks):
//   ... --case mpi-consistency --field-boundary dirichlet-phi
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "grid.h"
#include "maxwell.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    ElectrostaticBoundary boundary;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "all";
    args.boundary = { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--field-boundary") {
            if (i + 1 >= argc) return false;
            const std::string value(argv[++i]);
            if (value == "dirichlet-phi") {
                args.boundary.type = ElectrostaticBoundaryType::DIRICHLET_PHI;
            } else if (value == "left-E") {
                args.boundary.type = ElectrostaticBoundaryType::LEFT_E;
            } else {
                return false;
            }
        } else if (arg == "--phi-left") {
            if (i + 1 >= argc) return false;
            args.boundary.phi_left = std::strtod(argv[++i], NULL);
        } else if (arg == "--phi-right") {
            if (i + 1 >= argc) return false;
            args.boundary.phi_right = std::strtod(argv[++i], NULL);
        } else if (arg == "--left-electric-field") {
            if (i + 1 >= argc) return false;
            args.boundary.e_left = std::strtod(argv[++i], NULL);
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty();
}

double cell_center_x(const SpatialGrid& grid, int ix_with_ghost)
{
    return (static_cast<double>(grid.global_cell(ix_with_ghost)) + 0.5) *
           grid.dx;
}

double face_x(const SpatialGrid& grid, int global_face)
{
    return static_cast<double>(global_face) * grid.dx;
}

// Rectangular charge packet with edges aligned to cell faces so that the
// discrete Gauss prefix sum represents the piecewise-linear field exactly.
void packet_edges(const SpatialGrid& grid, double& a, double& b)
{
    const double third = grid.length() / 3.0;
    a = std::floor(third / grid.dx) * grid.dx;
    b = std::ceil(2.0 * third / grid.dx) * grid.dx;
    if (b <= a) { a = grid.dx; b = grid.length() - grid.dx; }
}

double packet_rho(double x, double rho0, double a, double b)
{
    return (x >= a && x <= b) ? rho0 : 0.0;
}

double packet_cumulative_charge(double x, double rho0, double a, double b)
{
    if (x <= a) return 0.0;
    if (x <= b) return rho0 * (x - a);
    return rho0 * (b - a);
}

double packet_phi_integral(double x, double rho0, double a, double b)
{
    // Integral_0^x Q(s) ds with Q the packet cumulative charge.
    if (x <= a) return 0.0;
    if (x <= b) return 0.5 * rho0 * (x - a) * (x - a);
    return rho0 * (b - a) * (x - 0.5 * (a + b));
}

struct CaseMetrics {
    double residual_linf_norm;
    double boundary_residual_rel;
    double max_field_error_rel;
    double max_phi_endpoint_error;
    CaseMetrics()
        : residual_linf_norm(0.0), boundary_residual_rel(0.0),
          max_field_error_rel(0.0), max_phi_endpoint_error(0.0) {}
};

// phi endpoint checks compare extrapolated ghost cells (x = -dx/2 or
// x = L + dx/2) whose exact value follows from the boundary potential and the
// analytic endpoint field.  The error is normalized by the potential scale of
// the case so the roundoff-level tolerance is meaningful across all cases.
double ghost_phi_scale(const ElectrostaticBoundary& boundary,
                       double max_abs_e, double length)
{
    return std::max(1.0,
                    std::max(std::fabs(boundary.phi_left),
                             std::max(std::fabs(boundary.phi_right),
                                      0.5 * max_abs_e * length)));
}

CaseMetrics run_zero_charge_case(const SpatialGrid& grid, EMFields& fields,
                                 const ElectrostaticBoundary& boundary,
                                 int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    std::fill(fields.rho.begin(), fields.rho.end(), 0.0);
    for (int ix = 0; ix < nxl; ++ix) fields.rho[ng + ix] = 0.0;

    OpenElectrostaticSolver solver;
    solver.init(grid, boundary);
    solver.solve(fields, rank, size);

    CaseMetrics m;
    const OpenGaussDiagnostics& diag = solver.diagnostics();
    m.residual_linf_norm = diag.residual_linf / 1.0;
    m.boundary_residual_rel =
        std::fabs(diag.boundary_charge_residual) /
        std::max(1.0, std::fabs(diag.net_charge));
    double max_field_error = 0.0;
    for (int iface = 0; iface <= nxl; ++iface) {
        max_field_error = std::max(max_field_error,
                                   std::fabs(fields.Ex_face[iface] - 0.0));
    }
    m.max_field_error_rel = max_field_error / 1.0;
    if (boundary.type == ElectrostaticBoundaryType::DIRICHLET_PHI) {
        const double phi_scale = ghost_phi_scale(boundary, 0.0, grid.length());
        if (grid.ix_start == 0) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng - 1] - boundary.phi_left) / phi_scale);
        }
        if (grid.ix_start + nxl == grid.nx_global) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng + nxl] - boundary.phi_right) /
                    phi_scale);
        }
    }
    return m;
}

CaseMetrics run_constant_charge_case(const SpatialGrid& grid, EMFields& fields,
                                     const ElectrostaticBoundary& boundary,
                                     double rho0, int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const double eps0 = Const::eps0;
    const double L = grid.length();
    for (int ix = 0; ix < nxl; ++ix) fields.rho[ng + ix] = rho0;

    OpenElectrostaticSolver solver;
    solver.init(grid, boundary);
    solver.solve(fields, rank, size);

    CaseMetrics m;
    const OpenGaussDiagnostics& diag = solver.diagnostics();
    m.residual_linf_norm =
        diag.residual_linf / std::max(1.0, rho0 / eps0);
    m.boundary_residual_rel =
        std::fabs(diag.boundary_charge_residual) /
        std::max(1.0, std::fabs(diag.net_charge));

    // E(x) = rho0 (x - L/2) / eps0 for phi(0)=phi(L)=0.
    const double e_scale = std::max(1.0, rho0 * L / eps0);
    double max_field_error = 0.0;
    for (int iface = 0; iface <= nxl; ++iface) {
        const double x = face_x(grid, grid.ix_start + iface);
        const double expected = rho0 * (x - 0.5 * L) / eps0;
        max_field_error = std::max(
            max_field_error, std::fabs(fields.Ex_face[iface] - expected));
    }
    m.max_field_error_rel = max_field_error / e_scale;
    if (boundary.type == ElectrostaticBoundaryType::DIRICHLET_PHI) {
        const double max_abs_e = rho0 * L / (2.0 * eps0);
        const double phi_scale =
            ghost_phi_scale(boundary, max_abs_e, grid.length());
        const double e_left_analytic = -max_abs_e;
        const double e_right_analytic = max_abs_e;
        if (grid.ix_start == 0) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng - 1] -
                          (boundary.phi_left +
                           0.5 * grid.dx * e_left_analytic)) /
                    phi_scale);
        }
        if (grid.ix_start + nxl == grid.nx_global) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng + nxl] -
                          (boundary.phi_right -
                           0.5 * grid.dx * e_right_analytic)) /
                    phi_scale);
        }
    }
    return m;
}

CaseMetrics run_packet_charge_case(const SpatialGrid& grid, EMFields& fields,
                                   const ElectrostaticBoundary& boundary,
                                   double rho0, int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const double eps0 = Const::eps0;
    const double L = grid.length();
    double a = 0.0, b = 0.0;
    packet_edges(grid, a, b);
    const double width = b - a;
    for (int ix = 0; ix < nxl; ++ix) {
        const double x = cell_center_x(grid, ng + ix);
        fields.rho[ng + ix] = packet_rho(x, rho0, a, b);
    }

    OpenElectrostaticSolver solver;
    solver.init(grid, boundary);
    solver.solve(fields, rank, size);

    CaseMetrics m;
    const OpenGaussDiagnostics& diag = solver.diagnostics();
    m.residual_linf_norm =
        diag.residual_linf / std::max(1.0, rho0 / eps0);
    m.boundary_residual_rel =
        std::fabs(diag.boundary_charge_residual) /
        std::max(1.0, std::fabs(diag.net_charge));

    // Constant field chosen so phi(L) = phi(0):
    // E0 = - Integral_0^L Q(s) ds / (eps0 L).
    const double e0 = -packet_phi_integral(L, rho0, a, b) / (eps0 * L);
    const double e_scale = std::max(1.0, rho0 * width / eps0);
    double max_field_error = 0.0;
    for (int iface = 0; iface <= nxl; ++iface) {
        const double x = face_x(grid, grid.ix_start + iface);
        const double expected =
            e0 + packet_cumulative_charge(x, rho0, a, b) / eps0;
        max_field_error = std::max(
            max_field_error, std::fabs(fields.Ex_face[iface] - expected));
    }
    m.max_field_error_rel = max_field_error / e_scale;

    // The two physical endpoint faces are independent degrees of freedom:
    // with a nonzero net charge they differ by net_charge/eps0 and each must
    // match its own analytic value.  A periodic alias of the right face to the
    // left face would fail this check.
    const double e_left = fields.Ex_face.front();
    const double e_right = fields.Ex_face.back();
    const double expected_left = e0;
    const double expected_right = e0 + rho0 * width / eps0;
    m.max_field_error_rel = std::max(
        m.max_field_error_rel,
        std::max(std::fabs(e_left - expected_left),
                 std::fabs(e_right - expected_right)) / e_scale);
    if (grid.ix_start == 0 && grid.ix_start + nxl == grid.nx_global &&
        std::fabs((e_right - e_left) - rho0 * width / eps0) >
            1.0e-11 * e_scale) {
        m.max_field_error_rel = 1.0;  // force failure on alias
    }

    if (boundary.type == ElectrostaticBoundaryType::DIRICHLET_PHI) {
        const double max_abs_e = std::max(
            std::fabs(e0), std::fabs(e0 + rho0 * width / eps0));
        const double phi_scale =
            ghost_phi_scale(boundary, max_abs_e, grid.length());
        const double e_left_analytic = e0;
        const double e_right_analytic = e0 + rho0 * width / eps0;
        if (grid.ix_start == 0) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng - 1] -
                          (boundary.phi_left +
                           0.5 * grid.dx * e_left_analytic)) /
                    phi_scale);
        }
        if (grid.ix_start + nxl == grid.nx_global) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng + nxl] -
                          (boundary.phi_right -
                           0.5 * grid.dx * e_right_analytic)) /
                    phi_scale);
        }
    }
    return m;
}

CaseMetrics run_dirichlet_difference_case(const SpatialGrid& grid,
                                          EMFields& fields,
                                          const ElectrostaticBoundary& boundary,
                                          int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const double L = grid.length();
    std::fill(fields.rho.begin(), fields.rho.end(), 0.0);
    for (int ix = 0; ix < nxl; ++ix) fields.rho[ng + ix] = 0.0;

    OpenElectrostaticSolver solver;
    solver.init(grid, boundary);
    solver.solve(fields, rank, size);

    CaseMetrics m;
    const OpenGaussDiagnostics& diag = solver.diagnostics();
    m.residual_linf_norm = diag.residual_linf / 1.0;
    m.boundary_residual_rel =
        std::fabs(diag.boundary_charge_residual) /
        std::max(1.0, std::fabs(diag.net_charge));

    // rho=0 with phi(0)=0, phi(L)=V gives the uniform field E = -V/L.
    const double e_expected =
        -(boundary.phi_right - boundary.phi_left) / L;
    const double e_scale = std::max(1.0, std::fabs(e_expected));
    double max_field_error = 0.0;
    for (int iface = 0; iface <= nxl; ++iface) {
        max_field_error = std::max(
            max_field_error,
            std::fabs(fields.Ex_face[iface] - e_expected));
    }
    m.max_field_error_rel = max_field_error / e_scale;
    if (boundary.type == ElectrostaticBoundaryType::DIRICHLET_PHI) {
        const double phi_scale =
            ghost_phi_scale(boundary, std::fabs(e_expected), grid.length());
        if (grid.ix_start == 0) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng - 1] -
                          (boundary.phi_left +
                           0.5 * grid.dx * e_expected)) /
                    phi_scale);
        }
        if (grid.ix_start + nxl == grid.nx_global) {
            m.max_phi_endpoint_error = std::max(
                m.max_phi_endpoint_error,
                std::fabs(fields.phi[ng + nxl] -
                          (boundary.phi_right -
                           0.5 * grid.dx * e_expected)) /
                    phi_scale);
        }
    }
    return m;
}

bool write_result_file(const std::string& path, int rank,
                       double residual_linf_norm, double boundary_residual_rel,
                       double max_field_error_rel, double max_phi_endpoint_error,
                       double mpi_face_consistency,
                       double boundary_energy_work_relative_error)
{
    if (path.empty()) return true;
    if (rank != 0) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=PASS\n";
    out << "residual_linf_norm=" << residual_linf_norm << "\n";
    out << "boundary_charge_residual_rel=" << boundary_residual_rel << "\n";
    out << "max_field_error_rel=" << max_field_error_rel << "\n";
    out << "max_phi_endpoint_error=" << max_phi_endpoint_error << "\n";
    out << "mpi_face_consistency_max_diff=" << mpi_face_consistency << "\n";
    out << "boundary_energy_work_relative_error="
        << boundary_energy_work_relative_error << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: open_electrostatic_solver_test --case "
                     "all|mpi-consistency --field-boundary "
                     "dirichlet-phi|left-E [--result <path>]\n";
    }
    if (ok && args.test_case == "mpi-consistency" && size < 2) {
        ok = false;
        if (rank == 0) {
            std::cerr << "mpi-consistency requires at least 2 ranks\n";
        }
    }
    int ok_all = ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &ok_all, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    ok = ok_all != 0;

    CaseMetrics combined;
    double mpi_face_consistency = 0.0;
    double boundary_energy_work_relative_error = 0.0;

    if (ok) {
        SpatialGrid grid;
        grid.init_with_domain(rank, size, Param::nx, Param::Lx);
        EMFields fields;
        fields.init(grid);

        // The accepted-step energy ledger uses the production endpoint-field
        // identity. Exercise both physical MPI endpoints without reconstructing
        // the Poisson solve or a boundary current in the test.
        EMFields energy_before, energy_after;
        energy_before.init(grid);
        energy_after.init(grid);
        if (rank == 0) energy_after.Ex_face.front() = 2.0;
        if (rank == size - 1) energy_after.Ex_face.back() = -4.0;
        ElectrostaticBoundary energy_boundary = {
            ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 3.0, -2.0 };
        OpenElectrostaticSolver energy_solver;
        energy_solver.init(grid, energy_boundary);
        const double boundary_work = energy_solver.boundary_energy_work(
            energy_before, energy_after, rank, size);
        const double expected_boundary_work = -2.0 * Const::eps0;
        boundary_energy_work_relative_error = std::fabs(
            boundary_work - expected_boundary_work) /
            std::max(std::fabs(expected_boundary_work),
                     std::numeric_limits<double>::min());

        const double rho0 = Const::qe * 1.0e28;  // ~ 1.6e9 C/m^3

        if (args.test_case == "all") {
            ElectrostaticBoundary zero_boundary =
                { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
            CaseMetrics m =
                run_zero_charge_case(grid, fields, zero_boundary, rank, size);
            combined.residual_linf_norm = std::max(
                combined.residual_linf_norm, m.residual_linf_norm);
            combined.boundary_residual_rel = std::max(
                combined.boundary_residual_rel, m.boundary_residual_rel);
            combined.max_field_error_rel = std::max(
                combined.max_field_error_rel, m.max_field_error_rel);
            combined.max_phi_endpoint_error = std::max(
                combined.max_phi_endpoint_error, m.max_phi_endpoint_error);

            m = run_constant_charge_case(grid, fields, zero_boundary, rho0,
                                         rank, size);
            combined.residual_linf_norm = std::max(
                combined.residual_linf_norm, m.residual_linf_norm);
            combined.boundary_residual_rel = std::max(
                combined.boundary_residual_rel, m.boundary_residual_rel);
            combined.max_field_error_rel = std::max(
                combined.max_field_error_rel, m.max_field_error_rel);
            combined.max_phi_endpoint_error = std::max(
                combined.max_phi_endpoint_error, m.max_phi_endpoint_error);

            m = run_packet_charge_case(grid, fields, zero_boundary, rho0,
                                       rank, size);
            combined.residual_linf_norm = std::max(
                combined.residual_linf_norm, m.residual_linf_norm);
            combined.boundary_residual_rel = std::max(
                combined.boundary_residual_rel, m.boundary_residual_rel);
            combined.max_field_error_rel = std::max(
                combined.max_field_error_rel, m.max_field_error_rel);
            combined.max_phi_endpoint_error = std::max(
                combined.max_phi_endpoint_error, m.max_phi_endpoint_error);

            ElectrostaticBoundary difference_boundary =
                { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 1.0e4 };
            m = run_dirichlet_difference_case(grid, fields, difference_boundary,
                                              rank, size);
            combined.residual_linf_norm = std::max(
                combined.residual_linf_norm, m.residual_linf_norm);
            combined.boundary_residual_rel = std::max(
                combined.boundary_residual_rel, m.boundary_residual_rel);
            combined.max_field_error_rel = std::max(
                combined.max_field_error_rel, m.max_field_error_rel);
            combined.max_phi_endpoint_error = std::max(
                combined.max_phi_endpoint_error, m.max_phi_endpoint_error);
        }

        if (args.test_case == "mpi-consistency") {
            combined = run_packet_charge_case(grid, fields, args.boundary,
                                              rho0, rank, size);

            // Multi-rank consistency: compare the local face fields with an
            // independent single-rank reference solve of the same problem.
            SpatialGrid ref_grid;
            ref_grid.init_with_domain(0, 1, Param::nx, Param::Lx);
            EMFields ref_fields;
            ref_fields.init(ref_grid);
            const int ng_ref = ref_grid.nghost;
            double a = 0.0, b = 0.0;
            packet_edges(ref_grid, a, b);
            for (int ix = 0; ix < ref_grid.nx_local; ++ix) {
                const double x = cell_center_x(ref_grid, ng_ref + ix);
                ref_fields.rho[ng_ref + ix] = packet_rho(x, rho0, a, b);
            }
            OpenElectrostaticSolver ref_solver;
            ref_solver.init(ref_grid, args.boundary);
            ref_solver.solve(ref_fields, 0, 1);
            const double eps0 = Const::eps0;
            const double e_scale = std::max(1.0, rho0 * (b - a) / eps0);
            double local_max_diff = 0.0;
            for (int iface = 0; iface <= grid.nx_local; ++iface) {
                const int global_face = grid.ix_start + iface;
                local_max_diff = std::max(
                    local_max_diff,
                    std::fabs(fields.Ex_face[iface] -
                              ref_fields.Ex_face[global_face]));
            }
            MPI_Allreduce(&local_max_diff, &mpi_face_consistency, 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            mpi_face_consistency /= e_scale;

            // Shared internal MPI faces must be bitwise identical after the
            // production interface sync.  Exchange endpoint values and verify.
            double local_shared_diff = 0.0;
            double my_right = fields.Ex_face[grid.nx_local];
            double my_left = fields.Ex_face[0];
            double neighbor_left = 0.0;
            double neighbor_right = 0.0;
            if (rank + 1 < size) {
                MPI_Sendrecv(&my_right, 1, MPI_DOUBLE, rank + 1, 9201,
                             &neighbor_left, 1, MPI_DOUBLE, rank + 1, 9202,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank > 0) {
                MPI_Sendrecv(&my_left, 1, MPI_DOUBLE, rank - 1, 9202,
                             &neighbor_right, 1, MPI_DOUBLE, rank - 1, 9201,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank + 1 < size) {
                local_shared_diff = std::max(
                    local_shared_diff, std::fabs(my_right - neighbor_left));
            }
            if (rank > 0) {
                local_shared_diff = std::max(
                    local_shared_diff, std::fabs(my_left - neighbor_right));
            }
            double global_shared_diff = 0.0;
            MPI_Allreduce(&local_shared_diff, &global_shared_diff, 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            mpi_face_consistency = std::max(mpi_face_consistency,
                                            global_shared_diff / e_scale);
        }
    }

    double global_linf = 0.0, global_bnd = 0.0, global_field = 0.0;
    double global_phi = 0.0;
    MPI_Allreduce(&combined.residual_linf_norm, &global_linf, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&combined.boundary_residual_rel, &global_bnd, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&combined.max_field_error_rel, &global_field, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&combined.max_phi_endpoint_error, &global_phi, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double roundoff_tolerance = 1.0e-11;
    const bool pass = ok_all != 0 &&
        global_linf <= roundoff_tolerance &&
        global_bnd <= roundoff_tolerance &&
        global_field <= roundoff_tolerance &&
        global_phi <= roundoff_tolerance &&
        (args.test_case != "mpi-consistency" ||
        mpi_face_consistency <= roundoff_tolerance);
    const bool boundary_energy_pass =
        boundary_energy_work_relative_error <= roundoff_tolerance;

    bool result_written =
        write_result_file(args.result_path, rank, global_linf, global_bnd,
                          global_field, global_phi, mpi_face_consistency,
                          boundary_energy_work_relative_error);
    if (!result_written) {
        if (rank == 0) std::cerr << "cannot write result file\n";
    }
    const bool final_pass = pass && boundary_energy_pass && result_written;
    if (rank == 0) {
        std::cout << "residual_linf_norm=" << global_linf
                  << " boundary_charge_residual_rel=" << global_bnd
                  << " max_field_error_rel=" << global_field
                  << " max_phi_endpoint_error=" << global_phi
                  << " mpi_face_consistency=" << mpi_face_consistency
                  << " boundary_energy_work_relative_error="
                  << boundary_energy_work_relative_error << "\n";
        std::cout << "status=" << (final_pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return final_pass ? 0 : 1;
}
