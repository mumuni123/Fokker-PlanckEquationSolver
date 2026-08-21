// Phase-1 acceptance test for the fixed production spatial grid contract
// (sections 13.1, 13.2 and 16.2.1).  The production parameters must be
// 8000 cells / 40 um / 0.005 um, and SpatialGrid::init_with_domain() must
// produce the same open non-periodic partition for every rank count.
//
// Usage:
//   grid_contract_test --expected-nx 8000 --expected-length-um 40
//                      --expected-dx-um 0.005 [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "grid.h"
#include "parameters.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ExpectedGrid {
    int nx;
    double length_um;
    double dx_um;
    std::string result_path;
};

bool parse_args(int argc, char** argv, ExpectedGrid& expected)
{
    expected.nx = 8000;
    expected.length_um = 40.0;
    expected.dx_um = 0.005;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--expected-nx") {
            if (i + 1 >= argc) return false;
            expected.nx = std::atoi(argv[++i]);
        } else if (arg == "--expected-length-um") {
            if (i + 1 >= argc) return false;
            expected.length_um = std::strtod(argv[++i], NULL);
        } else if (arg == "--expected-dx-um") {
            if (i + 1 >= argc) return false;
            expected.dx_um = std::strtod(argv[++i], NULL);
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            expected.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return expected.nx > 0 && expected.length_um > 0.0 && expected.dx_um > 0.0;
}

bool near(double a, double b, double scale)
{
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, scale);
    return std::fabs(a - b) <= tolerance;
}

bool expect_throw(int rank, int nranks, int nx, double length)
{
    SpatialGrid grid;
    try {
        grid.init_with_domain(rank, nranks, nx, length);
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
    ExpectedGrid expected;
    bool ok = parse_args(argc, argv, expected);
    if (!ok) {
        std::cerr << "usage: grid_contract_test --expected-nx 8000 "
                     "--expected-length-um 40 --expected-dx-um 0.005 "
                     "[--result <path>]\n";
    }

    double max_partition_error = 0.0;
    if (ok) {
        const double length = expected.length_um * Const::micro;
        const double dx = expected.dx_um * Const::micro;
        const int nx = expected.nx;

        // 1. Production parameter contract (section 13.1).
        if (Param::nx != nx) {
            std::cerr << "Param::nx = " << Param::nx << " != expected "
                      << nx << "\n";
            ok = false;
        }
        if (!near(Param::Lx, length, length)) {
            std::cerr << "Param::Lx = " << Param::Lx << " != expected "
                      << length << "\n";
            ok = false;
        }
        if (!near(Param::dx, dx, dx)) {
            std::cerr << "Param::dx = " << Param::dx << " != expected "
                      << dx << "\n";
            ok = false;
        }

        // 2. round(L/dx) * dx must reproduce L to floating-point tolerance.
        const double cells_ratio = length / Param::dx;
        const double cells = std::floor(cells_ratio + 0.5);
        if (!near(cells * Param::dx, length, length)) {
            std::cerr << "round(L/dx)*dx does not reproduce L\n";
            ok = false;
        }

        // 3. Decomposition invariants for 1, 2 and 5 ranks.
        const int rank_counts[3] = { 1, 2, 5 };
        for (int c = 0; c < 3 && ok; ++c) {
            const int nranks = rank_counts[c];
            std::vector<SpatialGrid> grids(static_cast<size_t>(nranks));
            for (int r = 0; r < nranks; ++r) {
                grids[static_cast<size_t>(r)].init_with_domain(
                    r, nranks, Param::nx, Param::Lx);
            }
            int covered = 0;
            int previous_end = 0;
            for (int r = 0; r < nranks; ++r) {
                const SpatialGrid& g = grids[static_cast<size_t>(r)];
                if (g.nx_global != Param::nx || g.ix_start != previous_end ||
                    g.nx_local <= 0) {
                    ok = false;
                }
                if (!near(g.dx, Param::dx, Param::dx) ||
                    !near(g.length(), Param::Lx, Param::Lx)) {
                    ok = false;
                }
                if (!near(g.left_boundary(), 0.0, 1.0) ||
                    !near(g.right_boundary(), Param::Lx, Param::Lx)) {
                    ok = false;
                }
                if (g.owns_left_physical_boundary(r) != (r == 0) ||
                    g.owns_right_physical_boundary(r, nranks) !=
                        (r == nranks - 1)) {
                    ok = false;
                }
                for (int ix = 0; ix < g.nx_local; ++ix) {
                    if (g.global_cell(g.nghost + ix) != g.ix_start + ix) {
                        ok = false;
                    }
                    const double x_expected =
                        (static_cast<double>(g.ix_start + ix) + 0.5) * g.dx;
                    if (!near(g.x(g.nghost + ix), x_expected,
                              std::max(1.0, x_expected))) {
                        ok = false;
                    }
                }
                previous_end = g.ix_start + g.nx_local;
                covered += g.nx_local;
            }
            if (covered != Param::nx) ok = false;

            // Interface faces between adjacent ranks share one global face.
            for (int r = 0; r + 1 < nranks; ++r) {
                const SpatialGrid& left = grids[static_cast<size_t>(r)];
                const SpatialGrid& right = grids[static_cast<size_t>(r + 1)];
                max_partition_error = std::max(
                    max_partition_error,
                    std::fabs((left.ix_start + left.nx_local) - right.ix_start));
            }
            if (max_partition_error != 0.0) ok = false;
        }

        // 4. init_with_domain() validation (section 13.2).
        const double nan_length =
            std::numeric_limits<double>::quiet_NaN();
        if (!expect_throw(0, 1, Param::nx, -1.0) ||
            !expect_throw(0, 1, Param::nx, 0.0) ||
            !expect_throw(0, 1, Param::nx, nan_length) ||
            !expect_throw(0, 1, 0, Param::Lx) ||
            !expect_throw(0, 0, Param::nx, Param::Lx) ||
            !expect_throw(1, 1, Param::nx, Param::Lx) ||
            !expect_throw(0, 5, 4, Param::Lx)) {
            ok = false;
        }

        // 5. Runtime beam macro weight is derived from the actual grid.dx.
        SpatialGrid grid;
        grid.init_with_domain(0, 1, Param::nx, Param::Lx);
        const double expected_weight =
            Param::densb * grid.dx /
            static_cast<double>(Param::beam_macro_particles_per_cell);
        if (!near(beam_macro_weight(grid), expected_weight,
                  std::max(1.0, expected_weight))) {
            ok = false;
        }
    }

    if (!ok) std::cerr << "grid contract violation detected\n";
    if (!expected.result_path.empty()) {
        std::ofstream out(expected.result_path.c_str(), std::ios::trunc);
        if (out) {
            out << "status=" << (ok ? "PASS" : "FAIL") << "\n";
            out << "nx=" << Param::nx << "\n";
            out << "length_um=" << (Param::Lx / Const::micro) << "\n";
            out << "dx_um=" << (Param::dx / Const::micro) << "\n";
            out << "partition_overlap_error=" << max_partition_error << "\n";
            out.close();
        } else {
            ok = false;
        }
    }
    std::cout << "nx=" << Param::nx << " length_um=" << (Param::Lx / Const::micro)
              << " dx_um=" << (Param::dx / Const::micro)
              << " partition_error=" << max_partition_error << "\n";
    std::cout << "status=" << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 1;
}
