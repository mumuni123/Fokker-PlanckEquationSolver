// Phase-2A acceptance test for the independent conservative x remap
// (sections 13.6, 14 phase 2A and 16.3.1/16.3.4).  It drives the production
// ConservativePpmRemap class only; analytic references are independent.
//
// Usage:
//   conservative_x_remap_test --case all [--result <path>]
//   conservative_x_remap_test --case smooth-convergence
//       --resolutions 100,200,400,800 [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "conservative_ppm_remap.h"
#include "grid.h"
#include "open_boundary.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::vector<int> resolutions;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "all";
    const int defaults[4] = { 100, 200, 400, 800 };
    args.resolutions.assign(defaults, defaults + 4);
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--resolutions") {
            if (i + 1 >= argc) return false;
            args.resolutions.clear();
            std::istringstream stream(argv[++i]);
            std::string token;
            while (std::getline(stream, token, ',')) {
                const int r = std::atoi(token.c_str());
                if (r <= 0) return false;
                args.resolutions.push_back(r);
            }
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty() && !args.resolutions.empty();
}

Species make_species(const SpatialGrid& grid)
{
    Species sp;
    sp.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
            Const::me, Param::dens, Param::temperature_e, false, grid);
    return sp;
}

size_t find_slice_near(const Species& sp, double target_vx)
{
    size_t best = 0;
    double best_error = std::numeric_limits<double>::infinity();
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const size_t q = idx2(iv, imu);
            const double error = std::fabs(sp.cgrid.vx[q] - target_vx);
            if (error < best_error) {
                best_error = error;
                best = q;
            }
        }
    }
    return best;
}

void set_slice_profile(Species& sp, const SpatialGrid& grid, size_t q,
                       const std::function<double(double)>& lambda_x)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    std::fill(sp.f.begin(), sp.f.end(), 0.0);
    const int j = static_cast<int>(q) / Param::Nmu;
    const int k = static_cast<int>(q) % Param::Nmu;
    for (int i = 0; i < nxl; ++i) {
        const double x = (static_cast<double>(grid.ix_start + i) + 0.5) * grid.dx;
        sp.f[idx3(ng + i, j, k)] = lambda_x(x) * grid.dx;
    }
}

void fill_reservoir_state(Species& sp, const SpatialGrid& grid,
                          const OpenBackgroundBoundary& boundary)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double lambda = boundary.incoming_cell_average(
                PhysicalSide::LEFT, iv, imu, 0.0, sp);
            for (int i = 0; i < nxl; ++i) {
                sp.f[idx3(ng + i, iv, imu)] = lambda * grid.dx;
            }
        }
    }
}

double total_mass(const Species& sp, const SpatialGrid& grid)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    double total = 0.0;
    for (int i = 0; i < nxl; ++i) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                total += sp.f[idx3(ng + i, iv, imu)];
            }
        }
    }
    return total;
}

double max_abs_mass(const Species& sp, const SpatialGrid& grid)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    double max_mass = 0.0;
    for (int i = 0; i < nxl; ++i) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                max_mass = std::max(max_mass,
                                    std::fabs(sp.f[idx3(ng + i, iv, imu)]));
            }
        }
    }
    return max_mass;
}

double min_mass(const Species& sp, const SpatialGrid& grid)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    double min_mass = std::numeric_limits<double>::infinity();
    for (int i = 0; i < nxl; ++i) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                min_mass = std::min(min_mass, sp.f[idx3(ng + i, iv, imu)]);
            }
        }
    }
    return min_mass;
}

OpenBackgroundBoundaryConfig absorbing_config()
{
    OpenBackgroundBoundaryConfig config;
    config.left_type = BackgroundXBoundaryType::ABSORBING;
    config.right_type = BackgroundXBoundaryType::ABSORBING;
    config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    return config;
}

OpenBackgroundBoundaryConfig reservoir_config()
{
    OpenBackgroundBoundaryConfig config;
    config.left_type = BackgroundXBoundaryType::RESERVOIR;
    config.right_type = BackgroundXBoundaryType::RESERVOIR;
    config.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    return config;
}

bool nonnegativity_ok(const RemapDiagnostics& diag, const Species& sp,
                      const SpatialGrid& grid)
{
    const double scale = std::max(1.0, max_abs_mass(sp, grid));
    const double bound =
        128.0 * std::numeric_limits<double>::epsilon() * scale;
    return diag.minimum_cell_mass >= -bound && min_mass(sp, grid) >= -bound;
}

struct CaseMetrics {
    double max_cell_diff_rel;
    double mass_rel_error;
    double com_shift_error_um;
    double plateau_max_err_rel;
    double l1_rel_error;
    double mirror_max_diff_rel;
    double final_total_mass;
    double max_inflow;
    double ledger_rel_error;
    double min_cell_mass;
    double order_200_400;
    double order_400_800;
    bool nonnegativity_ok;
    CaseMetrics()
        : max_cell_diff_rel(0.0), mass_rel_error(0.0), com_shift_error_um(0.0),
          plateau_max_err_rel(0.0), l1_rel_error(0.0),
          mirror_max_diff_rel(0.0), final_total_mass(0.0), max_inflow(0.0),
          ledger_rel_error(0.0), min_cell_mass(0.0), order_200_400(0.0),
          order_400_800(0.0), nonnegativity_ok(false) {}
};

CaseMetrics run_constant_and_reservoir_steady_state()
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 40.0 * Const::micro);
    Species sp = make_species(grid);
    Species out = make_species(grid);
    OpenBackgroundBoundary boundary(reservoir_config());
    fill_reservoir_state(sp, grid, boundary);
    ConservativePpmRemap remap;
    remap.init(grid, sp.cgrid);
    const double dt = 6.0 * grid.dx / (0.4 * Const::c);
    const RemapDiagnostics diag = remap.advect_x(sp, out, dt, 0.0, boundary, 0, 1);

    double max_diff = 0.0;
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    for (int i = 0; i < nxl; ++i) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                max_diff = std::max(
                    max_diff, std::fabs(out.f[idx3(ng + i, iv, imu)] -
                                        sp.f[idx3(ng + i, iv, imu)]));
            }
        }
    }
    m.max_cell_diff_rel =
        max_diff / std::max(1.0, max_abs_mass(sp, grid));
    m.mass_rel_error =
        std::fabs(diag.number_after - diag.number_before) /
        std::max(1.0, diag.number_before);
    m.ledger_rel_error =
        std::fabs((diag.number_after - diag.number_before) -
                  (diag.inflow_number - diag.outflow_number)) /
        std::max(1.0, diag.number_before);
    m.min_cell_mass = diag.minimum_cell_mass;
    m.nonnegativity_ok = nonnegativity_ok(diag, out, grid);
    return m;
}

CaseMetrics run_square_wave()
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 400, 40.0 * Const::micro);
    Species sp = make_species(grid);
    Species out = make_species(grid);
    OpenBackgroundBoundary boundary(absorbing_config());
    const size_t q = find_slice_near(sp, 0.4 * Const::c);
    const int j = static_cast<int>(q) / Param::Nmu;
    const int k = static_cast<int>(q) % Param::Nmu;
    const double v = sp.cgrid.vx[q];
    const double lambda0 = 1.0e20;
    const double box_lo = 12.0 * Const::micro;
    const double box_hi = 24.0 * Const::micro;
    set_slice_profile(sp, grid, q, [&](double x) {
        return (x >= box_lo && x <= box_hi) ? lambda0 : 0.0;
    });

    ConservativePpmRemap remap;
    remap.init(grid, sp.cgrid);
    const double dt = 3.0 * grid.dx / v;
    const RemapDiagnostics diag = remap.advect_x(sp, out, dt, 0.0, boundary, 0, 1);

    const double shift = v * dt;
    // Number conservation uses the remap ledger (same-order global sums), not
    // independent naive summations which carry ~1e-13 relative noise.
    m.mass_rel_error =
        std::fabs(diag.number_after - diag.number_before) /
        std::max(1.0, diag.number_before);

    // Center of mass displacement.
    double com_before = 0.0, com_after = 0.0;
    double slice_mass_before = 0.0, slice_mass_after = 0.0;
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    for (int i = 0; i < nxl; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * grid.dx;
        const double mb = sp.f[idx3(ng + i, j, k)];
        const double ma = out.f[idx3(ng + i, j, k)];
        com_before += x * mb;
        com_after += x * ma;
        slice_mass_before += mb;
        slice_mass_after += ma;
    }
    m.com_shift_error_um =
        std::fabs((com_after / std::max(1e-300, slice_mass_after)) -
                  (com_before / std::max(1e-300, slice_mass_before)) - shift) /
        Const::micro;

    // Interior plateau (at least 4 cells from the translated edges) must be
    // reproduced exactly.
    const double new_lo = box_lo + shift;
    const double new_hi = box_hi + shift;
    double plateau_max_err = 0.0;
    for (int i = 4; i < nxl - 4; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * grid.dx;
        if (x > new_lo + 4.0 * grid.dx && x < new_hi - 4.0 * grid.dx) {
            plateau_max_err = std::max(
                plateau_max_err,
                std::fabs(out.f[idx3(ng + i, j, k)] - lambda0 * grid.dx));
        }
    }
    m.plateau_max_err_rel =
        plateau_max_err / std::max(1.0, lambda0 * grid.dx);
    m.min_cell_mass = diag.minimum_cell_mass;
    m.ledger_rel_error =
        std::fabs((diag.number_after - diag.number_before) -
                  (diag.inflow_number - diag.outflow_number)) /
        std::max(1.0, diag.number_before);
    m.nonnegativity_ok = nonnegativity_ok(diag, out, grid);
    return m;
}

CaseMetrics run_smooth_gaussian(bool convergence,
                                const std::vector<int>& resolutions,
                                double& order_200_400, double& order_400_800)
{
    CaseMetrics m;
    if (!convergence) {
        SpatialGrid grid;
        grid.init_with_domain(0, 1, 400, 40.0 * Const::micro);
        Species sp = make_species(grid);
        Species out = make_species(grid);
        OpenBackgroundBoundary boundary(absorbing_config());
        const size_t q = find_slice_near(sp, 0.3 * Const::c);
        const int j = static_cast<int>(q) / Param::Nmu;
        const int k = static_cast<int>(q) % Param::Nmu;
        const double v = sp.cgrid.vx[q];
        const double sigma = 3.0 * Const::micro;
        const double center = 20.0 * Const::micro;
        const double amplitude = 1.0e21;
        const double distance = 2.4 * Const::micro;
        set_slice_profile(sp, grid, q, [&](double x) {
            return amplitude *
                   std::exp(-(x - center) * (x - center) /
                            (2.0 * sigma * sigma));
        });
        ConservativePpmRemap remap;
        remap.init(grid, sp.cgrid);
        const double dt = distance / v;
        const RemapDiagnostics diag = remap.advect_x(sp, out, dt, 0.0, boundary, 0, 1);

        double l1 = 0.0;
        const int ng = grid.nghost;
        const int nxl = grid.nx_local;
        for (int i = 0; i < nxl; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * grid.dx;
            // Exact absorbing-boundary solution (zero where the profile would
            // enter from outside the domain).
            const double expected =
                (x >= distance)
                    ? amplitude * std::exp(-(x - distance - center) *
                                           (x - distance - center) /
                                           (2.0 * sigma * sigma)) * grid.dx
                    : 0.0;
            l1 += std::fabs(out.f[idx3(ng + i, j, k)] - expected);
        }
        m.l1_rel_error = l1 / std::max(1.0, total_mass(sp, grid));
        // The profile reaches the absorbing boundaries, so the physically
        // correct conservation statement is the ledger identity
        // after - before == inflow - outflow.
        m.mass_rel_error =
            std::fabs((diag.number_after - diag.number_before) -
                      (diag.inflow_number - diag.outflow_number)) /
            std::max(1.0, diag.number_before);
        m.min_cell_mass = diag.minimum_cell_mass;
        m.nonnegativity_ok = nonnegativity_ok(diag, out, grid);
        return m;
    }

    // Grid-convergence: fixed physical advection distance at each resolution.
    const double sigma = 3.0 * Const::micro;
    const double center = 20.0 * Const::micro;
    const double amplitude = 1.0e21;
    // Non-integer-cell shift: the swept mass then depends on the PPM
    // reconstruction, so the observed L1 order is meaningful (integer-cell
    // shifts are exact for conservative remaps and only measure roundoff).
    const double distance = 2.37 * Const::micro;
    std::vector<double> l1_errors;
    for (size_t r = 0; r < resolutions.size(); ++r) {
        SpatialGrid grid;
        grid.init_with_domain(0, 1, resolutions[r], 40.0 * Const::micro);
        Species sp = make_species(grid);
        Species out = make_species(grid);
        OpenBackgroundBoundary boundary(absorbing_config());
        const size_t q = find_slice_near(sp, 0.3 * Const::c);
        const int j = static_cast<int>(q) / Param::Nmu;
        const int k = static_cast<int>(q) % Param::Nmu;
        const double v = sp.cgrid.vx[q];
        set_slice_profile(sp, grid, q, [&](double x) {
            return amplitude *
                   std::exp(-(x - center) * (x - center) /
                            (2.0 * sigma * sigma));
        });
        ConservativePpmRemap remap;
        remap.init(grid, sp.cgrid);
        const double dt = distance / v;
        remap.advect_x(sp, out, dt, 0.0, boundary, 0, 1);

        double l1 = 0.0;
        const int ng = grid.nghost;
        const int nxl = grid.nx_local;
        for (int i = 0; i < nxl; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * grid.dx;
            // Exact absorbing-boundary solution: translated profile for
            // x >= v*dt, zero where the profile would come from outside.
            const double expected =
                (x >= distance)
                    ? amplitude * std::exp(-(x - distance - center) *
                                           (x - distance - center) /
                                           (2.0 * sigma * sigma)) * grid.dx
                    : 0.0;
            l1 += std::fabs(out.f[idx3(ng + i, j, k)] - expected);
        }
        l1_errors.push_back(l1 / std::max(1.0, total_mass(sp, grid)));
    }
    if (l1_errors.size() >= 4) {
        m.order_200_400 =
            std::log2(l1_errors[l1_errors.size() - 3] /
                      l1_errors[l1_errors.size() - 2]);
        m.order_400_800 =
            std::log2(l1_errors[l1_errors.size() - 2] /
                      l1_errors[l1_errors.size() - 1]);
        order_200_400 = m.order_200_400;
        order_400_800 = m.order_400_800;
    }
    return m;
}

CaseMetrics run_velocity_symmetry()
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 400, 40.0 * Const::micro);
    Species sp = make_species(grid);
    Species out = make_species(grid);
    OpenBackgroundBoundary boundary(absorbing_config());
    const size_t q_plus = find_slice_near(sp, 0.4 * Const::c);
    const size_t q_minus = find_slice_near(sp, -0.4 * Const::c);
    const int jp = static_cast<int>(q_plus) / Param::Nmu;
    const int kp = static_cast<int>(q_plus) % Param::Nmu;
    const int jm = static_cast<int>(q_minus) / Param::Nmu;
    const int km = static_cast<int>(q_minus) % Param::Nmu;
    const double sigma = 3.0 * Const::micro;
    const double center = 20.0 * Const::micro;
    const double amplitude = 1.0e21;
    const double distance = 2.4 * Const::micro;

    std::fill(sp.f.begin(), sp.f.end(), 0.0);
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    for (int i = 0; i < nxl; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * grid.dx;
        const double value = amplitude *
            std::exp(-(x - center) * (x - center) / (2.0 * sigma * sigma)) *
            grid.dx;
        sp.f[idx3(ng + i, jp, kp)] = value;
        sp.f[idx3(ng + i, jm, km)] = value;
    }
    ConservativePpmRemap remap;
    remap.init(grid, sp.cgrid);
    const double dt = distance / sp.cgrid.vx[q_plus];
    remap.advect_x(sp, out, dt, 0.0, boundary, 0, 1);

    double max_diff = 0.0;
    double max_mass = 0.0;
    for (int i = 0; i < nxl; ++i) {
        const double a = out.f[idx3(ng + i, jp, kp)];
        const double b = out.f[idx3(ng + nxl - 1 - i, jm, km)];
        max_diff = std::max(max_diff, std::fabs(a - b));
        max_mass = std::max(max_mass, std::fabs(a));
    }
    m.mirror_max_diff_rel = max_diff / std::max(1.0, max_mass);
    m.nonnegativity_ok = true;
    return m;
}

CaseMetrics run_absorbing_no_wrap()
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 100, 40.0 * Const::micro);
    Species sp = make_species(grid);
    OpenBackgroundBoundary boundary(absorbing_config());
    const size_t q = find_slice_near(sp, 0.8 * Const::c);
    const double v = sp.cgrid.vx[q];
    const double lambda0 = 1.0e20;
    set_slice_profile(sp, grid, q, [&](double x) {
        return (x >= 2.0 * Const::micro && x <= 6.0 * Const::micro)
                   ? lambda0 : 0.0;
    });

    ConservativePpmRemap remap;
    remap.init(grid, sp.cgrid);
    const double dt = 5.0 * grid.dx / v;
    double max_inflow = 0.0;
    double max_ledger_error = 0.0;
    const double initial_mass = total_mass(sp, grid);
    double before = initial_mass;
    bool nonneg = true;
    for (int step = 0; step < 25; ++step) {
        const RemapDiagnostics diag = remap.advect_x(sp, sp, dt, 0.0,
                                                     boundary, 0, 1);
        max_inflow = std::max(max_inflow, diag.inflow_number);
        const double after = total_mass(sp, grid);
        max_ledger_error = std::max(
            max_ledger_error,
            std::fabs((after - before) -
                      (diag.inflow_number - diag.outflow_number)) /
                std::max(1.0, before));
        before = after;
        m.min_cell_mass = std::min(m.min_cell_mass, diag.minimum_cell_mass);
        nonneg = nonneg && nonnegativity_ok(diag, sp, grid);
    }
    m.final_total_mass = total_mass(sp, grid) / std::max(1.0, initial_mass);
    m.max_inflow = max_inflow;
    m.ledger_rel_error = max_ledger_error;
    m.nonnegativity_ok = nonneg;
    return m;
}

bool write_result_file(const std::string& path, const CaseMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "max_cell_diff_rel=" << m.max_cell_diff_rel << "\n";
    out << "mass_rel_error=" << m.mass_rel_error << "\n";
    out << "com_shift_error_um=" << m.com_shift_error_um << "\n";
    out << "plateau_max_err_rel=" << m.plateau_max_err_rel << "\n";
    out << "l1_rel_error=" << m.l1_rel_error << "\n";
    out << "mirror_max_diff_rel=" << m.mirror_max_diff_rel << "\n";
    out << "final_total_mass=" << m.final_total_mass << "\n";
    out << "max_inflow=" << m.max_inflow << "\n";
    out << "ledger_rel_error=" << m.ledger_rel_error << "\n";
    out << "min_cell_mass=" << m.min_cell_mass << "\n";
    out << "order_200_400=" << m.order_200_400 << "\n";
    out << "order_400_800=" << m.order_400_800 << "\n";
    out << "nonnegativity_ok=" << (m.nonnegativity_ok ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        // The serial cases in this test assume one process holds the full
        // domain and advect with (rank=0, size=1).  Under a multi-rank
        // launcher the production MPI_COMM_WORLD reductions would sum
        // identical per-rank states, corrupting the ledgers.  The MPI
        // interface case lives in conservative_x_remap_mpi_test.
        std::cerr << "conservative_x_remap_test must run with exactly 1 rank; "
                     "use plain ./build/conservative_x_remap_test (no "
                     "yhrun/mpirun). MPI cases: conservative_x_remap_mpi_test.\n";
        MPI_Finalize();
        return 2;
    }
    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: conservative_x_remap_test --case "
                     "all|smooth-convergence [--resolutions ...] "
                     "[--result <path>]\n";
    }

    CaseMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "smooth-convergence") {
        double o2 = 0.0, o3 = 0.0;
        m = run_smooth_gaussian(true, args.resolutions, o2, o3);
        // PPM order approaches 3; the monotonicity limiter at the Gaussian
        // peak allows a local order reduction (section 15.5).
        pass = o3 >= 2.3 && o2 >= 2.2;
    } else if (ok && args.test_case == "all") {
        CaseMetrics a = run_constant_and_reservoir_steady_state();
        CaseMetrics b = run_square_wave();
        double o2 = 0.0, o3 = 0.0;
        CaseMetrics c = run_smooth_gaussian(false, args.resolutions, o2, o3);
        CaseMetrics d = run_velocity_symmetry();
        CaseMetrics e = run_absorbing_no_wrap();
        m = a;
        m.max_cell_diff_rel = std::max(a.max_cell_diff_rel, b.max_cell_diff_rel);
        m.mass_rel_error =
            std::max(std::max(a.mass_rel_error, b.mass_rel_error),
                     c.mass_rel_error);
        m.com_shift_error_um = b.com_shift_error_um;
        m.plateau_max_err_rel = b.plateau_max_err_rel;
        m.l1_rel_error = c.l1_rel_error;
        m.mirror_max_diff_rel = d.mirror_max_diff_rel;
        m.final_total_mass = e.final_total_mass;
        m.max_inflow = e.max_inflow;
        m.ledger_rel_error =
            std::max(std::max(a.ledger_rel_error, b.ledger_rel_error),
                     e.ledger_rel_error);
        m.min_cell_mass =
            std::min(std::min(a.min_cell_mass, b.min_cell_mass),
                     std::min(c.min_cell_mass, e.min_cell_mass));

        const double roundoff = 1.0e-10;
        pass = a.max_cell_diff_rel <= roundoff &&
               a.mass_rel_error <= 1.0e-12 &&
               b.mass_rel_error <= 1.0e-12 &&
               b.com_shift_error_um <= 1.0e-3 &&
               b.plateau_max_err_rel <= 1.0e-10 &&
               c.l1_rel_error <= 1.0e-3 &&
               c.mass_rel_error <= 1.0e-12 &&
               d.mirror_max_diff_rel <= 1.0e-8 &&
               e.final_total_mass <= 1.0e-10 &&
               e.max_inflow == 0.0 &&
               e.ledger_rel_error <= 1.0e-10 &&
               a.nonnegativity_ok && b.nonnegativity_ok &&
               c.nonnegativity_ok && e.nonnegativity_ok;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "max_cell_diff_rel=" << m.max_cell_diff_rel
              << " mass_rel_error=" << m.mass_rel_error
              << " com_shift_error_um=" << m.com_shift_error_um
              << " plateau_max_err_rel=" << m.plateau_max_err_rel
              << " l1_rel_error=" << m.l1_rel_error
              << " mirror_max_diff_rel=" << m.mirror_max_diff_rel
              << " final_total_mass=" << m.final_total_mass
              << " max_inflow=" << m.max_inflow
              << " ledger_rel_error=" << m.ledger_rel_error
              << " min_cell_mass=" << m.min_cell_mass
              << " order_200_400=" << m.order_200_400
              << " order_400_800=" << m.order_400_800 << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
