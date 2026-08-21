// Stage H1 acceptance: staggered face-field gather (section 6.6.1).  The
// production ParticleShape1D::shifted_face_weights applied to a nonlinear
// Ex_face must agree with the hand-computed shifted CIC interpolation and
// must use the face values directly (never a cell-averaged field).
//
// Usage:
//   background_tail_stagger_test [--result <path>]
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

struct Metrics {
    double max_gather_reference_diff_rel;
    double face_value_error;
    double non_cell_average_gap;
    bool reference_agrees;
    Metrics()
        : max_gather_reference_diff_rel(0.0), face_value_error(0.0),
          non_cell_average_gap(0.0), reference_agrees(false)
    {}
};

double gather(const FaceGatherWeights& w,
              const std::vector<double>& ex_face, int ix_start)
{
    const int il0 = w.face0 - ix_start;
    const int il1 = w.face1 - ix_start;
    return w.w0 * ex_face[static_cast<size_t>(il0)] +
           w.w1 * ex_face[static_cast<size_t>(il1)];
}

// Hand-computed shifted CIC reference (the same formula, written out as the
// independent reference for section 6.6.1).
double reference_gather(double x, double dx,
                        const std::vector<double>& ex_face)
{
    const double face_coordinate = x / dx;
    const int center = static_cast<int>(std::floor(face_coordinate + 0.5));
    const double shift = face_coordinate - static_cast<double>(center);
    double ex = (1.0 - std::fabs(shift)) *
                ex_face[static_cast<size_t>(center)];
    if (shift < 0.0) {
        ex += (-shift) * ex_face[static_cast<size_t>(center - 1)];
    } else if (shift > 0.0) {
        ex += shift * ex_face[static_cast<size_t>(center + 1)];
    }
    return ex;
}

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 400, 4.0 * Const::micro);
    const double dx = grid.dx;
    const double length = dx * static_cast<double>(grid.nx_global);

    // Strongly nonlinear face field.
    std::vector<double> ex_face(static_cast<size_t>(grid.nx_local) + 1);
    for (size_t j = 0; j < ex_face.size(); ++j) {
        const double xf = static_cast<double>(j) * dx / Const::micro;
        ex_face[j] = 1.0e12 * (1.0 + 0.4 * std::sin(2.0 * Const::pi * xf /
                                                     0.37) +
                                0.2 * std::cos(2.0 * Const::pi * xf / 0.11));
    }

    double max_rel = 0.0;
    const int nsamples = 2000;
    for (int s = 0; s < nsamples; ++s) {
        const double x = (static_cast<double>(s) + 0.5) / nsamples * length;
        const FaceGatherWeights gw =
            ParticleShape1D::shifted_face_weights(x, grid);
        const double g = gather(gw, ex_face, grid.ix_start);
        const double g_ref = reference_gather(x, dx, ex_face);
        const double scale = std::max(1.0, std::fabs(g_ref));
        max_rel = std::max(max_rel, std::fabs(g - g_ref) / scale);
    }
    m.max_gather_reference_diff_rel = max_rel;
    m.reference_agrees = max_rel <= 1.0e-12;

    // A particle exactly on a face must return that face value (not any
    // cell-centered combination).
    {
        std::vector<double> e(static_cast<size_t>(grid.nx_local) + 1, 0.0);
        e[1] = 10.0;
        e[2] = 0.0;
        const FaceGatherWeights gw_at_face =
            ParticleShape1D::shifted_face_weights(1.0 * dx, grid);
        const double g_face = gather(gw_at_face, e, grid.ix_start);
        m.face_value_error = std::fabs(g_face - 10.0);

        // At x = 1.25 dx the CIC face gather is 0.75*E[1] + 0.25*E[2] = 7.5;
        // a cell-averaged scheme would return 5.0.  This proves the gather
        // operates on the staggered face values directly.
        const FaceGatherWeights gw_quarter =
            ParticleShape1D::shifted_face_weights(1.25 * dx, grid);
        const double g_quarter = gather(gw_quarter, e, grid.ix_start);
        m.non_cell_average_gap = std::fabs(g_quarter - 7.5);
    }
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "max_gather_reference_diff_rel="
        << m.max_gather_reference_diff_rel << "\n";
    out << "face_value_error=" << m.face_value_error << "\n";
    out << "non_cell_average_gap=" << m.non_cell_average_gap << "\n";
    out << "reference_agrees=" << (m.reference_agrees ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_stagger_test must run with exactly "
                     "1 rank; use plain ./build_hybrid/"
                     "background_tail_stagger_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: background_tail_stagger_test [--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case();
    const bool pass = ok && m.reference_agrees &&
                      m.face_value_error <= 1.0e-12 &&
                      m.non_cell_average_gap <= 1.0e-12;
    if (!write_result_file(args.result_path, m, pass)) return 2;
    std::cout << "max_gather_reference_diff_rel="
              << m.max_gather_reference_diff_rel
              << " face_value_error=" << m.face_value_error
              << " non_cell_average_gap=" << m.non_cell_average_gap
              << " reference_agrees=" << (m.reference_agrees ? 1 : 0)
              << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
