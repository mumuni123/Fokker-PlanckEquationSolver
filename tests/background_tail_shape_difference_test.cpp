// Stage H1 acceptance: shape-difference trajectory current (sections 6.3
// and 6.6.2).  Random CIC segments that stay inside one cell, cross one
// cell face, or cross multiple cell faces must satisfy
//   Delta(n dx) + D_x(number flux) dt = 0
// per cell to summation error, using only the production
// ChargeConservingTrajectory1D::deposit_segment and ParticleShape1D::cell_weights.
//
// Usage:
//   background_tail_shape_difference_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "grid.h"

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
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

unsigned long long lcg_state = 0xabcdef0123456789ULL;

double random_unit()
{
    lcg_state = lcg_state * 2862933555777941757ULL + 3037000493ULL;
    return static_cast<double>(lcg_state >> 11) * (1.0 / 9007199254740992.0);
}

struct Metrics {
    double max_residual_same_cell;
    double max_residual_one_face;
    double max_residual_multi_cell;
    double total_charge_error_rel;
    size_t checked_segments;
    Metrics()
        : max_residual_same_cell(0.0), max_residual_one_face(0.0),
          max_residual_multi_cell(0.0), total_charge_error_rel(0.0),
          checked_segments(0)
    {}
};

double cell_share(const CellDepositWeights& cw, int cell)
{
    if (cw.cell0 == cell) return cw.w0;
    if (cw.cell1 == cell) return cw.w1;
    return 0.0;
}

// Checks one segment class against the per-cell discrete continuity.
double check_class(const SpatialGrid& grid, double x0, double x1,
                   double weight, double dt, double& total_charge_error,
                   size_t& checked)
{
    const double inv_qe = 1.0 / Const::qe;
    FaceCurrentAccumulator acc;
    acc.init(grid);
    ChargeConservingTrajectory1D::deposit_segment(
        x0, x1, weight, dt, grid, acc);
    const CellDepositWeights w0 = ParticleShape1D::cell_weights(x0, grid);
    const CellDepositWeights w1 = ParticleShape1D::cell_weights(x1, grid);
    double max_residual = 0.0;
    double charge_change_total = 0.0;
    for (int i = 0; i < grid.nx_local; ++i) {
        const double dn = weight *
            (cell_share(w1, i) - cell_share(w0, i));
        const double div_j =
            (acc.current_face_x[static_cast<size_t>(i + 1)] -
             acc.current_face_x[static_cast<size_t>(i)]) * dt * inv_qe;
        const double residual = dn - div_j;
        max_residual = std::max(max_residual, std::fabs(residual));
        charge_change_total += dn;
    }
    total_charge_error += std::fabs(charge_change_total);
    ++checked;
    return max_residual;
}

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    const double dx = grid.dx;
    const double length = dx * static_cast<double>(grid.nx_global);
    const double lo = 0.75 * dx;
    const double hi = length - 0.75 * dx;
    const double weight = 1.0e20;
    const double dt = 1.0e-15;
    const int nseg = 200;

    for (int s = 0; s < nseg; ++s) {
        // Same-cell segment (does not cross a cell face).
        const int cell = 5 + static_cast<int>(random_unit() * 190.0);
        const double xc = (cell + 0.5) * dx;
        const double x0 = xc + (random_unit() - 0.5) * 0.2 * dx;
        const double x1 = x0 + (random_unit() - 0.5) * 0.2 * dx;
        double charge_error = 0.0;
        size_t checked = 0;
        m.max_residual_same_cell = std::max(
            m.max_residual_same_cell,
            check_class(grid, x0, x1, weight, dt, charge_error, checked));
        m.total_charge_error_rel += charge_error;
        m.checked_segments += checked;
    }
    for (int s = 0; s < nseg; ++s) {
        // Crosses exactly one cell face.
        const int face = 10 + static_cast<int>(random_unit() * 180.0);
        const double x0 = face * dx - (0.1 + 0.3 * random_unit()) * dx;
        const double x1 = face * dx + (0.1 + 0.3 * random_unit()) * dx;
        double charge_error = 0.0;
        size_t checked = 0;
        m.max_residual_one_face = std::max(
            m.max_residual_one_face,
            check_class(grid, x0, x1, weight, dt, charge_error, checked));
        m.total_charge_error_rel += charge_error;
        m.checked_segments += checked;
    }
    for (int s = 0; s < nseg; ++s) {
        // Crosses multiple cell faces.
        const double x0 = lo + random_unit() * (hi - lo);
        const double dir = (random_unit() < 0.5) ? -1.0 : 1.0;
        const double x1 = std::max(
            lo, std::min(hi, x0 + dir * (3.0 + 6.0 * random_unit()) * dx));
        double charge_error = 0.0;
        size_t checked = 0;
        m.max_residual_multi_cell = std::max(
            m.max_residual_multi_cell,
            check_class(grid, x0, x1, weight, dt, charge_error, checked));
        m.total_charge_error_rel += charge_error;
        m.checked_segments += checked;
    }
    // Segments stay inside the domain, so the total charge change over all
    // checked segments must vanish.
    m.total_charge_error_rel /=
        std::max(1.0, weight * static_cast<double>(m.checked_segments));
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "max_residual_same_cell=" << m.max_residual_same_cell << "\n";
    out << "max_residual_one_face=" << m.max_residual_one_face << "\n";
    out << "max_residual_multi_cell=" << m.max_residual_multi_cell << "\n";
    out << "total_charge_error_rel=" << m.total_charge_error_rel << "\n";
    out << "checked_segments=" << m.checked_segments << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_shape_difference_test must run with "
                     "exactly 1 rank; use plain ./build_hybrid/"
                     "background_tail_shape_difference_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: background_tail_shape_difference_test "
                     "[--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case();
    const double tol = 1.0e-10 * 1.0e20;  // 1e-10 relative to the weight
    const bool pass = ok &&
                      m.max_residual_same_cell <= tol &&
                      m.max_residual_one_face <= tol &&
                      m.max_residual_multi_cell <= tol &&
                      m.total_charge_error_rel <= 1.0e-12;
    if (!write_result_file(args.result_path, m, pass)) return 2;
    std::cout << "max_residual_same_cell=" << m.max_residual_same_cell
              << " max_residual_one_face=" << m.max_residual_one_face
              << " max_residual_multi_cell=" << m.max_residual_multi_cell
              << " total_charge_error_rel=" << m.total_charge_error_rel
              << " checked_segments=" << m.checked_segments << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
