// Phase-4 acceptance test for in-step injection events (sections 13.9/16.6.2):
// a particle crossing before the midpoint enters the midpoint density at its
// real drift position, and one crossing after the midpoint is created in the
// second half and pushed for tau = t^{n+1} - t_cross.  Drives the production
// BeamPIC predict_to_midpoint / finish_from_midpoint pair.
//
// Usage:
//   beam_injection_event_test --case before-and-after-midpoint [--result <path>]
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
#include <string>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "before-and-after-midpoint";
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

double vx_from_px(double px)
{
    const double u = px / (Const::me * Const::c);
    return Const::c * u / std::sqrt(1.0 + u * u);
}

struct CaseMetrics {
    double midpoint_count_error;
    double midpoint_position_rel;
    double midpoint_density_rel;
    double final_count_error;
    double early_position_rel;
    double late_position_rel;
    double number_balance_rel;
    double zero_field_position_rel;
    CaseMetrics()
        : midpoint_count_error(0.0), midpoint_position_rel(0.0),
          midpoint_density_rel(0.0), final_count_error(0.0),
          early_position_rel(0.0), late_position_rel(0.0),
          number_balance_rel(0.0), zero_field_position_rel(0.0)
    {}
};

CaseMetrics run_case()
{
    CaseMetrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 400, 40.0 * Const::micro);
    // Deliberately larger than the production step so the injected particles
    // span several cells (the CIC deposit near x=0 otherwise leaks most of
    // the shape weight outside the open boundary, which is expected but not
    // what this case measures).
    const double dt = 1.0e-14;
    const double time = 0.0;
    const double t_mid = time + 0.5 * dt;
    const double p0 = 0.5 * Const::me * Const::c;
    const double v0 = vx_from_px(p0);
    const double weight = 2.0e14;
    const double e_field = -1.0e12;

    BeamInjectionSchedule schedule;
    BeamInjectionEvent early = { time + 0.25 * dt, p0, weight };
    BeamInjectionEvent late = { time + 0.75 * dt, p0, weight };
    schedule.events.push_back(early);
    schedule.events.push_back(late);

    // --- E = 0: pure drifts, exact positions -----------------------------
    {
        BeamPIC beam;
        beam.init(grid);
        BeamPIC work;
        work.init(grid);
        EMFields fields;
        fields.init(grid);
        beam.predict_to_midpoint(schedule, grid, fields, time, dt, 0, 1, work);
        // Midpoint state: only the early particle, at v0*(t_mid - t_cross).
        if (work.particles.size() == 1) {
            const double expected_mid = v0 * (t_mid - early.crossing_time);
            m.midpoint_position_rel =
                std::fabs(work.particles[0].x - expected_mid) /
                std::max(1e-30, expected_mid);
        } else {
            m.midpoint_count_error = 1.0;
        }
        work.deposit_density(grid, 0, 1);
        double mid_total = 0.0;
        for (size_t ix = 0; ix < work.density.size(); ++ix) {
            mid_total += work.density[ix] * grid.dx;
        }
        m.midpoint_density_rel =
            std::fabs(mid_total - weight) / std::max(1.0, weight);

        work.finish_from_midpoint(schedule, grid, fields, time, dt, 0, 1);
        if (work.particles.size() == 2) {
            const double expected_early = v0 * (time + dt - early.crossing_time);
            const double expected_late = v0 * (time + dt - late.crossing_time);
            // Sort by position to identify the two injected particles.
            double xa = work.particles[0].x;
            double xb = work.particles[1].x;
            if (xa > xb) std::swap(xa, xb);
            // The early particle (born first) ends further downstream, so the
            // smaller position is the late particle.
            m.early_position_rel =
                std::fabs(xb - expected_early) / std::max(1e-30, expected_early);
            m.late_position_rel =
                std::fabs(xa - expected_late) / std::max(1e-30, expected_late);
            m.zero_field_position_rel =
                std::max(m.early_position_rel, m.late_position_rel);
        } else {
            m.final_count_error = 1.0;
        }
        const double injected = work.last_injected_number();
        const double total = work.total_particle_number(grid);
        m.number_balance_rel =
            std::fabs(total - injected) / std::max(1.0, injected);
    }

    // --- Constant E: kicked positions ------------------------------------
    // Early particle: midpoint drift v0*(dt/4), then kick over
    // (t^{n+1}-t_cross) and a second drift of dt/2.  Late particle: created
    // in the second half, kicked and drifted over tau = dt/4.
    {
        BeamPIC beam;
        beam.init(grid);
        BeamPIC work;
        work.init(grid);
        EMFields fields;
        fields.init(grid);
        std::fill(fields.Ex.begin(), fields.Ex.end(), e_field);
        std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), e_field);
        beam.predict_to_midpoint(schedule, grid, fields, time, dt, 0, 1, work);
        work.finish_from_midpoint(schedule, grid, fields, time, dt, 0, 1);
        if (work.particles.size() == 2) {
            const double tau_early = time + dt - early.crossing_time;
            const double tau_late = time + dt - late.crossing_time;
            const double p_early = p0 + (-Const::qe) * e_field * tau_early;
            const double p_late = p0 + (-Const::qe) * e_field * tau_late;
            const double x_early_expected =
                v0 * (t_mid - early.crossing_time) +
                vx_from_px(p_early) * (time + dt - t_mid);
            const double x_late_expected = vx_from_px(p_late) * tau_late;
            double xa = work.particles[0].x;
            double xb = work.particles[1].x;
            if (xa > xb) std::swap(xa, xb);
            const double rel = std::max(
                std::fabs(xb - x_early_expected) /
                    std::max(1e-30, x_early_expected),
                std::fabs(xa - x_late_expected) /
                    std::max(1e-30, x_late_expected));
            m.early_position_rel = std::max(m.early_position_rel, rel);
        } else {
            m.final_count_error = std::max(m.final_count_error, 1.0);
        }
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
    out << "midpoint_count_error=" << m.midpoint_count_error << "\n";
    out << "midpoint_position_rel=" << m.midpoint_position_rel << "\n";
    out << "midpoint_density_rel=" << m.midpoint_density_rel << "\n";
    out << "final_count_error=" << m.final_count_error << "\n";
    out << "early_position_rel=" << m.early_position_rel << "\n";
    out << "late_position_rel=" << m.late_position_rel << "\n";
    out << "number_balance_rel=" << m.number_balance_rel << "\n";
    out << "zero_field_position_rel=" << m.zero_field_position_rel << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "beam_injection_event_test must run with exactly 1 rank; "
                     "use plain ./build/beam_injection_event_test (no "
                     "yhrun/mpirun).\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: beam_injection_event_test --case "
                     "before-and-after-midpoint [--result <path>]\n";
    }

    CaseMetrics m;
    bool pass = ok;
    if (ok && args.test_case == "before-and-after-midpoint") {
        m = run_case();
        pass = m.midpoint_count_error == 0.0 &&
               m.midpoint_position_rel <= 1.0e-12 &&
               m.midpoint_density_rel <= 1.0e-12 &&
               m.final_count_error == 0.0 &&
               m.zero_field_position_rel <= 1.0e-12 &&
               m.early_position_rel <= 1.0e-9 &&
               m.late_position_rel <= 1.0e-9 &&
               m.number_balance_rel <= 1.0e-12;
    } else {
        pass = false;
    }

    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "midpoint_count_error=" << m.midpoint_count_error
              << " midpoint_position_rel=" << m.midpoint_position_rel
              << " midpoint_density_rel=" << m.midpoint_density_rel
              << " final_count_error=" << m.final_count_error
              << " early_position_rel=" << m.early_position_rel
              << " late_position_rel=" << m.late_position_rel
              << " number_balance_rel=" << m.number_balance_rel
              << " zero_field_position_rel=" << m.zero_field_position_rel
              << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
