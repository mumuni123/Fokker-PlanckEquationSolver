// Stage H1 acceptance for the tail relativistic D-K-D pusher (sections 6.3,
// 6.6.6 and 15 H1):
//   * zero-field free drift is an exact analytic drift;
//   * constant-field momentum kick is exact and the position is second order
//     in time (dt, dt/2, dt/4);
//   * with B=0 the D-K-D result agrees with the full-relativistic Boris
//     electrostatic specialization (no magnetic rotation is needed).
//
// Usage:
//   background_tail_pusher_test --case zero-and-constant-field
//       [--dt-scales 1.0,0.5,0.25] [--result <path>]
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

double vx_from_u(double u)
{
    return Const::c * u / std::sqrt(1.0 + u * u);
}

double exact_u(double u0, double e_field, double t)
{
    const double force = -Const::qe * e_field;
    return u0 + force * t / (Const::me * Const::c);
}

double exact_position(double x0, double u0, double e_field, double t)
{
    const double force = -Const::qe * e_field;
    if (force == 0.0) return x0 + vx_from_u(u0) * t;
    const double ut = exact_u(u0, e_field, t);
    return x0 + (Const::me * Const::c * Const::c / force) *
                (std::sqrt(1.0 + ut * ut) - std::sqrt(1.0 + u0 * u0));
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
    // Electrostatic Boris specialization comparison.
    double boris_u_error_rel;
    double dkd_boris_diff_um_1;
    double dkd_boris_diff_um_3;
    double boris_pos_order;
    CaseMetrics()
        : pos_error_um_zero(0.0), pos_error_um_1(0.0), pos_error_um_2(0.0),
          pos_error_um_3(0.0), mom_error_rel_1(0.0), mom_error_rel_2(0.0),
          mom_error_rel_3(0.0), order_pos_12(0.0), order_pos_24(0.0),
          boris_u_error_rel(0.0), dkd_boris_diff_um_1(0.0),
          dkd_boris_diff_um_3(0.0), boris_pos_order(0.0)
    {}
};

void advance_dkd(BackgroundTailPIC& tail, const SpatialGrid& grid,
                 const EMFields& fields, double dt, int steps)
{
    for (int s = 0; s < steps; ++s) {
        tail.drift_half(grid, 0.5 * dt, 0, 1);
        tail.kick(grid, fields, dt, 0, 1);
        tail.drift_half(grid, 0.5 * dt, 0, 1);
    }
}

// Reference: full-relativistic Boris with B = 0 (only the electrostatic
// kick remains; there is no magnetic rotation to port).
void advance_boris_es(BackgroundTailPIC& tail, const SpatialGrid& grid,
                      const EMFields& fields, double dt, int steps)
{
    const double half_kick = 0.5 * Const::qe / (Const::me * Const::c);
    for (int s = 0; s < steps; ++s) {
        for (size_t i = 0; i < tail.particles.size(); ++i) {
            BackgroundTailParticle& p = tail.particles[i];
            const FaceGatherWeights gw =
                ParticleShape1D::shifted_face_weights(p.x, grid);
            const double ex =
                gw.w0 * fields.Ex_face[static_cast<size_t>(gw.face0)] +
                gw.w1 * fields.Ex_face[static_cast<size_t>(gw.face1)];
            // Boris: u^- = u^n - (q dt/(2 m c)) E; with B=0 the rotation is
            // the identity; u^+ = u^- - (q dt/(2 m c)) E.
            p.ux -= half_kick * ex * dt;
            p.ux -= half_kick * ex * dt;
            // Boris advances position with the new velocity over the full
            // step.
            const double gamma =
                std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
            p.x += Const::c * p.ux / gamma * dt;
        }
    }
}

CaseMetrics run_case(const std::vector<double>& dt_scales)
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 2000, 200.0 * Const::micro);
    const double dt0 = 2.5e-15;
    const double total_time = 1.0e-13;
    const double x0 = 20.0 * Const::micro;
    const double u0 = 1.0;
    const double e_field = -1.0e12;

    // E = 0: the symmetric D-K-D is an exact drift (two half-drifts compose).
    {
        BackgroundTailPIC tail;
        tail.init(grid);
        BackgroundTailParticle p;
        p.x = x0;
        p.ux = u0;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0;
        p.id = tail.next_particle_id(0);
        tail.particles.push_back(p);
        EMFields fields;
        fields.init(grid);
        const int steps = static_cast<int>(total_time / dt0 + 0.5);
        advance_dkd(tail, grid, fields, dt0, steps);
        const double expected = exact_position(x0, u0, 0.0, total_time);
        m.pos_error_um_zero =
            std::fabs(tail.particles[0].x - expected) / Const::micro;
    }

    // Constant field: run at each dt scale and read the second-order rate.
    std::vector<double> pos_errors;
    std::vector<double> mom_errors;
    for (size_t s = 0; s < dt_scales.size(); ++s) {
        const double dt = dt_scales[s] * dt0;
        const int steps = static_cast<int>(total_time / dt + 0.5);
        BackgroundTailPIC tail;
        tail.init(grid);
        BackgroundTailParticle p;
        p.x = x0;
        p.ux = u0;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0;
        p.id = tail.next_particle_id(0);
        tail.particles.push_back(p);
        EMFields fields;
        fields.init(grid);
        std::fill(fields.Ex.begin(), fields.Ex.end(), e_field);
        std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), e_field);
        advance_dkd(tail, grid, fields, dt, steps);
        const double x_expected = exact_position(x0, u0, e_field, total_time);
        const double u_expected = exact_u(u0, e_field, total_time);
        pos_errors.push_back(
            std::fabs(tail.particles[0].x - x_expected) / Const::micro);
        mom_errors.push_back(
            std::fabs(tail.particles[0].ux - u_expected) /
            std::max(1.0, std::fabs(u_expected)));
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
    }

    // Electrostatic Boris specialization (section 6.6.6): both schemes need
    // only the E kick; D-K-D is the second-order symmetric split while the
    // classical Boris position update is first order.  Both converge to the
    // same analytic solution and their difference shrinks with dt.
    {
        const double dt = dt_scales[0] * dt0;
        const int steps = static_cast<int>(total_time / dt + 0.5);
        BackgroundTailPIC dkd;
        dkd.init(grid);
        BackgroundTailParticle p;
        p.x = x0;
        p.ux = u0;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0;
        p.id = 1;
        dkd.particles.push_back(p);
        EMFields fields;
        fields.init(grid);
        std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), e_field);
        advance_dkd(dkd, grid, fields, dt, steps);
        BackgroundTailPIC boris = dkd;
        boris.particles[0].x = x0;
        boris.particles[0].ux = u0;
        boris.particles[0].uy = 0.0;
        boris.particles[0].uz = 0.0;
        advance_boris_es(boris, grid, fields, dt, steps);

        const double x_exact = exact_position(x0, u0, e_field, total_time);
        const double u_exact = exact_u(u0, e_field, total_time);
        m.dkd_boris_diff_um_1 =
            std::fabs(dkd.particles[0].x - boris.particles[0].x) /
            Const::micro;
        m.boris_u_error_rel =
            std::fabs(boris.particles[0].ux - u_exact) /
            std::max(1.0, std::fabs(u_exact));

        const double dt_fine = dt_scales[2] * dt0;
        const int steps_fine = static_cast<int>(total_time / dt_fine + 0.5);
        BackgroundTailPIC dkd_fine;
        dkd_fine.init(grid);
        p.x = x0;
        p.ux = u0;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = 1.0;
        p.id = 2;
        dkd_fine.particles.push_back(p);
        advance_dkd(dkd_fine, grid, fields, dt_fine, steps_fine);
        BackgroundTailPIC boris_fine = dkd_fine;
        boris_fine.particles[0].x = x0;
        boris_fine.particles[0].ux = u0;
        boris_fine.particles[0].uy = 0.0;
        boris_fine.particles[0].uz = 0.0;
        advance_boris_es(boris_fine, grid, fields, dt_fine, steps_fine);
        m.dkd_boris_diff_um_3 =
            std::fabs(dkd_fine.particles[0].x - boris_fine.particles[0].x) /
            Const::micro;
        // First-order Boris position error vs the analytic solution.
        const double boris_err_coarse =
            std::fabs(boris.particles[0].x - x_exact) / Const::micro;
        const double boris_err_fine =
            std::fabs(boris_fine.particles[0].x - x_exact) / Const::micro;
        m.boris_pos_order =
            std::log2(boris_err_coarse / std::max(1e-16, boris_err_fine));
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
    out << "boris_u_error_rel=" << m.boris_u_error_rel << "\n";
    out << "dkd_boris_diff_um_1=" << m.dkd_boris_diff_um_1 << "\n";
    out << "dkd_boris_diff_um_3=" << m.dkd_boris_diff_um_3 << "\n";
    out << "boris_pos_order=" << m.boris_pos_order << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "background_tail_pusher_test must run with exactly "
                     "1 rank; use plain ./build_hybrid/background_tail_"
                     "pusher_test (no yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: background_tail_pusher_test --case "
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
               // The constant-field kick is exact at every resolution.
               m.mom_error_rel_1 <= 1.0e-12 &&
               m.mom_error_rel_2 <= 1.0e-12 &&
               m.mom_error_rel_3 <= 1.0e-12 &&
               // Electrostatic Boris specialization: same momentum to
               // roundoff, first-order position convergence, and the D-K-D /
               // Boris difference shrinks at least ~3x over a 4x dt
               // refinement (both approaches converge to the same solution).
               m.boris_u_error_rel <= 1.0e-12 &&
               m.boris_pos_order >= 0.7 &&
               m.dkd_boris_diff_um_1 > 0.0 &&
               m.dkd_boris_diff_um_1 /
                   std::max(1e-16, m.dkd_boris_diff_um_3) >= 3.0;
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
              << " boris_u_error_rel=" << m.boris_u_error_rel
              << " dkd_boris_diff_um_1=" << m.dkd_boris_diff_um_1
              << " dkd_boris_diff_um_3=" << m.dkd_boris_diff_um_3
              << " boris_pos_order=" << m.boris_pos_order << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
