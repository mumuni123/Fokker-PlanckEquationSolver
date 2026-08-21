// Phase-4 acceptance test for the symmetric Beam drift-kick-drift pusher
// (sections 7.2/13.9 and 16.6.1).  It drives the production BeamPIC
// predict_to_midpoint / finish_from_midpoint pair; analytic references are
// the exact constant-field relativistic single-particle solution.
//
// Usage:
//   beam_leapfrog_test --case zero-and-constant-field
//       --dt-scales 1.0,0.5,0.25 [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "parameters.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::vector<double> dt_scales;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "zero-and-constant-field";
    const double defaults[3] = { 1.0, 0.5, 0.25 };
    args.dt_scales.assign(defaults, defaults + 3);
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--dt-scales") {
            if (i + 1 >= argc) return false;
            args.dt_scales.clear();
            std::istringstream stream(argv[++i]);
            std::string token;
            while (std::getline(stream, token, ',')) {
                const double s = std::strtod(token.c_str(), NULL);
                if (!(s > 0.0)) return false;
                args.dt_scales.push_back(s);
            }
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty() && !args.dt_scales.empty();
}

double vx_from_px(double px)
{
    const double u = px / (Const::me * Const::c);
    return Const::c * u / std::sqrt(1.0 + u * u);
}

double exact_position(double x0, double p0, double e_field, double t)
{
    const double f = (-Const::qe) * e_field;
    if (f == 0.0) return x0 + vx_from_px(p0) * t;
    const double pt = p0 + f * t;
    const double u0 = p0 / (Const::me * Const::c);
    const double ut = pt / (Const::me * Const::c);
    return x0 + (Const::me * Const::c * Const::c / f) *
                (std::sqrt(1.0 + ut * ut) - std::sqrt(1.0 + u0 * u0));
}

double exact_momentum(double p0, double e_field, double t)
{
    return p0 + (-Const::qe) * e_field * t;
}

struct CaseMetrics {
    double pos_error_um_zero;
    double pos_error_um_1;
    double pos_error_um_2;
    double pos_error_um_3;
    double mom_error_rel_1;
    double mom_error_rel_2;
    double mom_error_rel_3;
    double order_pos_12;
    double order_pos_24;
    double order_mom_12;
    double order_mom_24;
    CaseMetrics()
        : pos_error_um_zero(0.0), pos_error_um_1(0.0), pos_error_um_2(0.0),
          pos_error_um_3(0.0), mom_error_rel_1(0.0), mom_error_rel_2(0.0),
          mom_error_rel_3(0.0), order_pos_12(0.0), order_pos_24(0.0),
          order_mom_12(0.0), order_mom_24(0.0)
    {}
};

CaseMetrics run_case(const std::vector<double>& dt_scales)
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 2000, 200.0 * Const::micro);
    const double dt0 = 2.5e-15;
    const double total_time = 1.0e-13;
    const double x0 = 20.0 * Const::micro;
    const double p0 = Const::me * Const::c;
    // Negative field accelerates the electron forward (right); the particle
    // stays inside the 200 um domain for the whole run.
    const double e_field = -1.0e12;

    // E = 0: the leapfrog is an exact drift (two half-drifts compose).
    {
        BeamPIC beam;
        beam.init(grid);
        BeamParticle particle = { x0, p0, 1.0 };
        beam.particles.push_back(particle);
        BeamPIC work;
        work.init(grid);
        EMFields fields;
        fields.init(grid);
        BeamInjectionSchedule empty;
        const int steps = static_cast<int>(total_time / dt0 + 0.5);
        for (int s = 0; s < steps; ++s) {
            beam.predict_to_midpoint(empty, grid, fields,
                                     static_cast<double>(s) * dt0, dt0,
                                     0, 1, work);
            work.finish_from_midpoint(empty, grid, fields,
                                      static_cast<double>(s) * dt0, dt0,
                                      0, 1);
            beam.swap_state(work);
        }
        const double expected = exact_position(x0, p0, 0.0, total_time);
        m.pos_error_um_zero =
            std::fabs(beam.particles[0].x - expected) / Const::micro;
    }

    // Constant field: run at each dt scale and read the second-order rate.
    std::vector<double> pos_errors;
    std::vector<double> mom_errors;
    for (size_t s = 0; s < dt_scales.size(); ++s) {
        const double dt = dt_scales[s] * dt0;
        const int steps = static_cast<int>(total_time / dt + 0.5);
        BeamPIC beam;
        beam.init(grid);
        BeamParticle particle = { x0, p0, 1.0 };
        beam.particles.push_back(particle);
        BeamPIC work;
        work.init(grid);
        EMFields fields;
        fields.init(grid);
        std::fill(fields.Ex.begin(), fields.Ex.end(), e_field);
        std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), e_field);
        BeamInjectionSchedule empty;
        for (int step = 0; step < steps; ++step) {
            beam.predict_to_midpoint(empty, grid, fields,
                                     static_cast<double>(step) * dt, dt,
                                     0, 1, work);
            work.finish_from_midpoint(empty, grid, fields,
                                      static_cast<double>(step) * dt, dt,
                                      0, 1);
            beam.swap_state(work);
        }
        const double x_expected = exact_position(x0, p0, e_field, total_time);
        const double p_expected = exact_momentum(p0, e_field, total_time);
        if (beam.particles.empty()) {
            pos_errors.push_back(1.0e30);
            mom_errors.push_back(1.0e30);
        } else {
            pos_errors.push_back(
                std::fabs(beam.particles[0].x - x_expected) / Const::micro);
            mom_errors.push_back(
                std::fabs(beam.particles[0].px - p_expected) /
                std::max(1.0, std::fabs(p_expected)));
        }
    }
    if (pos_errors.size() >= 3) {
        m.pos_error_um_1 = pos_errors[0];
        m.pos_error_um_2 = pos_errors[1];
        m.pos_error_um_3 = pos_errors[2];
        m.mom_error_rel_1 = mom_errors[0];
        m.mom_error_rel_2 = mom_errors[1];
        m.mom_error_rel_3 = mom_errors[2];
        m.order_pos_12 =
            std::log2(pos_errors[0] / std::max(1e-16, pos_errors[1]));
        m.order_pos_24 =
            std::log2(pos_errors[1] / std::max(1e-16, pos_errors[2]));
        m.order_mom_12 =
            std::log2(mom_errors[0] / std::max(1e-16, mom_errors[1]));
        m.order_mom_24 =
            std::log2(mom_errors[1] / std::max(1e-16, mom_errors[2]));
    }
    return m;
}

bool write_result_file(const std::string& path, const CaseMetrics& m,
                       bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "pos_error_um_zero=" << m.pos_error_um_zero << "\n";
    out << "pos_error_um_1=" << m.pos_error_um_1 << "\n";
    out << "pos_error_um_2=" << m.pos_error_um_2 << "\n";
    out << "pos_error_um_3=" << m.pos_error_um_3 << "\n";
    out << "mom_error_rel_1=" << m.mom_error_rel_1 << "\n";
    out << "mom_error_rel_2=" << m.mom_error_rel_2 << "\n";
    out << "mom_error_rel_3=" << m.mom_error_rel_3 << "\n";
    out << "order_pos_12=" << m.order_pos_12 << "\n";
    out << "order_pos_24=" << m.order_pos_24 << "\n";
    out << "order_mom_12=" << m.order_mom_12 << "\n";
    out << "order_mom_24=" << m.order_mom_24 << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "beam_leapfrog_test must run with exactly 1 rank; use "
                     "plain ./build/beam_leapfrog_test (no yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: beam_leapfrog_test --case "
                     "zero-and-constant-field [--dt-scales 1.0,0.5,0.25] "
                     "[--result <path>]\n";
    }

    CaseMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "zero-and-constant-field") {
        m = run_case(args.dt_scales);
        pass = m.pos_error_um_zero <= 1.0e-6 &&
               m.order_pos_12 >= 1.5 && m.order_pos_12 <= 2.6 &&
               m.order_pos_24 >= 1.5 && m.order_pos_24 <= 2.6 &&
               // The constant-field kick is exact, so the momentum error is
               // roundoff-level at every resolution.
               m.mom_error_rel_1 <= 1.0e-12 &&
               m.mom_error_rel_2 <= 1.0e-12 &&
               m.mom_error_rel_3 <= 1.0e-12;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "pos_error_um_zero=" << m.pos_error_um_zero
              << " pos_error_um_1=" << m.pos_error_um_1
              << " pos_error_um_2=" << m.pos_error_um_2
              << " pos_error_um_3=" << m.pos_error_um_3
              << " mom_error_rel_1=" << m.mom_error_rel_1
              << " mom_error_rel_2=" << m.mom_error_rel_2
              << " mom_error_rel_3=" << m.mom_error_rel_3
              << " order_pos_12=" << m.order_pos_12
              << " order_pos_24=" << m.order_pos_24
              << " order_mom_12=" << m.order_mom_12
              << " order_mom_24=" << m.order_mom_24 << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
