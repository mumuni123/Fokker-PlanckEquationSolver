#include "background_coupling_test_support.h"
#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    SpatialGrid sg;
    sg.init(rank, size);
    Species background;
    EMFields fields;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.05, 0.08, 1.0e3, 3, 0.0);
    const double dt = BackgroundCouplingTest::stable_dt(sg);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                size, dt);

    double local_max_speed_difference = 0.0;
    double local_l2_speed_difference = 0.0;
    long long local_speed_count = 0;
    for (int j = 1; j < Param::Nv - 1; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const double analytic = background.cgrid.vx[idx2(j, k)];
            const double discrete = Stage5::energy_consistent_cell_speed_candidate(
                background.cgrid, background.mass, j, k, analytic);
            const double difference = analytic - discrete;
            local_max_speed_difference = std::max(local_max_speed_difference,
                                                  std::fabs(difference));
            local_l2_speed_difference += difference * difference;
            ++local_speed_count;
        }
    }
    double global_max_speed_difference = 0.0;
    double global_l2_speed_difference = 0.0;
    long long global_speed_count = 0;
    MPI_Allreduce(&local_max_speed_difference, &global_max_speed_difference, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_l2_speed_difference, &global_l2_speed_difference, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_speed_count, &global_speed_count, 1,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    std::vector<double> jn_energy_speed(static_cast<size_t>(sg.nx_local + 1),
                                        0.0);
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        double gamma = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t q = idx2(j, k);
                const double analytic = background.cgrid.vx[q];
                const double discrete = Stage5::energy_consistent_cell_speed_candidate(
                    background.cgrid, background.mass, j, k, analytic);
                const double flux = bundle.fx_high[
                    (static_cast<size_t>(iface) * Param::Nv + j) * Param::Nmu + k];
                gamma += (std::fabs(analytic) > 1.0e-300)
                    ? flux * discrete / analytic : flux;
            }
        }
        jn_energy_speed[static_cast<size_t>(iface)] = background.charge * gamma;
    }
    if (sg.nx_local > 0)
        jn_energy_speed[static_cast<size_t>(sg.nx_local)] = jn_energy_speed[0];
    const BackgroundCouplingTest::Norms analytic_pair =
        BackgroundCouplingTest::face_difference_norms(
            bundle.jn_high, bundle.gstar_je_center, sg);
    const BackgroundCouplingTest::Norms discrete_pair =
        BackgroundCouplingTest::face_difference_norms(
            jn_energy_speed, bundle.gstar_je_center, sg);
    const BackgroundCouplingTest::Norms analytic_work =
        BackgroundCouplingTest::face_work_residual_norms(
            bundle.jn_high, bundle.gstar_je_center, fields, sg, dt);
    const BackgroundCouplingTest::Norms discrete_work =
        BackgroundCouplingTest::face_work_residual_norms(
            jn_energy_speed, bundle.gstar_je_center, fields, sg, dt);
    double local_r_analytic = 0.0;
    double local_r_discrete = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double e = fields.Ex_face[static_cast<size_t>(iface)];
        local_r_analytic += dt * sg.dx * e *
            (bundle.jn_high[static_cast<size_t>(iface)] -
             bundle.gstar_je_center[static_cast<size_t>(iface)]);
        local_r_discrete += dt * sg.dx * e *
            (jn_energy_speed[static_cast<size_t>(iface)] -
             bundle.gstar_je_center[static_cast<size_t>(iface)]);
    }
    double r_values[2] = {local_r_analytic, local_r_discrete};
    MPI_Allreduce(MPI_IN_PLACE, r_values, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double pair_l2_ratio = discrete_pair.l2 /
        std::max(1.0e-300, analytic_pair.l2);
    const double work_l2_ratio = discrete_work.l2 /
        std::max(1.0e-300, analytic_work.l2);
    const int discrete_velocity_improves = work_l2_ratio < 1.0;
    int operator_valid = bundle.state_advanced && bundle.finite &&
        std::isfinite(pair_l2_ratio) && std::isfinite(work_l2_ratio) &&
        global_speed_count > 0;
    MPI_Allreduce(MPI_IN_PLACE, &operator_valid, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    // A finite production evaluation can legitimately disprove the direct
    // velocity-replacement hypothesis.  Do not report that as a passing
    // physical conclusion merely because the operator itself completed.
    const int direct_velocity_replacement_hypothesis_pass =
        operator_valid && discrete_velocity_improves;
    if (rank == 0) {
        std::ofstream out("output/background_coupling_velocity_consistency.result");
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_velocity_consistency\n"
            << "production_kernel=1\nbeam_enabled=0\n"
            << "max_abs_vx_minus_vE=" << global_max_speed_difference << "\n"
            << "rms_vx_minus_vE="
            << std::sqrt(global_l2_speed_difference /
                         std::max(1LL, global_speed_count)) << "\n"
            << "JN_analytic_minus_GstarJE_L1=" << analytic_pair.l1 << "\n"
            << "JN_analytic_minus_GstarJE_L2=" << analytic_pair.l2 << "\n"
            << "JN_analytic_minus_GstarJE_Linf=" << analytic_pair.linf << "\n"
            << "JN_discrete_speed_minus_GstarJE_L1=" << discrete_pair.l1 << "\n"
            << "JN_discrete_speed_minus_GstarJE_L2=" << discrete_pair.l2 << "\n"
            << "JN_discrete_speed_minus_GstarJE_Linf=" << discrete_pair.linf << "\n"
            << "R_pair_analytic_L1=" << analytic_work.l1 << "\n"
            << "R_pair_analytic_L2=" << analytic_work.l2 << "\n"
            << "R_pair_analytic_Linf=" << analytic_work.linf << "\n"
            << "R_pair_discrete_speed_L1=" << discrete_work.l1 << "\n"
            << "R_pair_discrete_speed_L2=" << discrete_work.l2 << "\n"
            << "R_pair_discrete_speed_Linf=" << discrete_work.linf << "\n"
            << "discrete_pair_L2_ratio=" << pair_l2_ratio << "\n"
            << "discrete_work_L2_ratio=" << work_l2_ratio << "\n"
            << "discrete_velocity_improves=" << discrete_velocity_improves << "\n"
            << "operator_valid=" << operator_valid << "\n"
            << "direct_velocity_replacement_hypothesis_pass="
            << direct_velocity_replacement_hypothesis_pass << "\n"
            << "R_pair_analytic_velocity=" << r_values[0] << "\n"
            << "R_pair_discrete_energy_velocity=" << r_values[1] << "\n"
            << "passes=" << direct_velocity_replacement_hypothesis_pass << "\n";
        std::cout << "background_coupling_velocity_consistency_test passes="
                  << direct_velocity_replacement_hypothesis_pass
                  << " operator_valid=" << operator_valid << "\n";
    }
    MPI_Finalize();
    // A negative scientific result is a completed audit, not an MPI failure.
    // The formal conclusion is carried by the .result `passes` key.
    return operator_valid ? 0 : 1;
}
