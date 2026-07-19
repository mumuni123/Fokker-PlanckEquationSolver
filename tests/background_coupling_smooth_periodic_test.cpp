#include "background_coupling_test_support.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

using BackgroundCouplingTest::Norms;

namespace {

struct CaseResult {
    Norms current_norms;
    Norms work_norms;
    int state_advanced;
    int finite;
    int fct_active;
};

CaseResult run_case(const SpatialGrid& sg, int rank, int size,
                    double field_amplitude, double phase)
{
    Species background;
    EMFields fields;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.04, 0.06, field_amplitude, 3,
        phase);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(
            background, fields, sg, rank, size,
            BackgroundCouplingTest::stable_dt(sg));
    CaseResult result = {};
    result.current_norms = BackgroundCouplingTest::face_difference_norms(
        bundle.jn_high, bundle.gstar_je_center, sg);
    result.work_norms = BackgroundCouplingTest::face_work_residual_norms(
        bundle.jn_high, bundle.gstar_je_center, fields, sg,
        BackgroundCouplingTest::stable_dt(sg));
    result.state_advanced = bundle.state_advanced;
    result.finite = bundle.finite;
    result.fct_active = bundle.fct_active;
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SpatialGrid sg;
    sg.init(rank, size);
    const double phase_shift = 2.0 * Const::pi * 7.0 /
        static_cast<double>(sg.nx_global);
    const CaseResult full = run_case(sg, rank, size, 1.0e3, 0.0);
    const CaseResult half = run_case(sg, rank, size, 5.0e2, 0.0);
    const CaseResult shifted = run_case(sg, rank, size, 1.0e3, phase_shift);

    const double current_scale_l2 = half.current_norms.l2 > 0.0
        ? full.current_norms.l2 / half.current_norms.l2 : 0.0;
    const double work_scale_l1 = half.work_norms.l1 > 0.0
        ? full.work_norms.l1 / half.work_norms.l1 : 0.0;
    const double work_scale_l2 = half.work_norms.l2 > 0.0
        ? full.work_norms.l2 / half.work_norms.l2 : 0.0;
    const double work_scale_linf = half.work_norms.linf > 0.0
        ? full.work_norms.linf / half.work_norms.linf : 0.0;
    const double phase_relative = std::fabs(full.current_norms.l2 -
                                             shifted.current_norms.l2) /
        std::max(1.0e-300, std::max(full.current_norms.l2,
                                    shifted.current_norms.l2));
    int local_ok = full.state_advanced && half.state_advanced &&
        shifted.state_advanced && full.finite && half.finite && shifted.finite;
    MPI_Allreduce(MPI_IN_PLACE, &local_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    if (rank == 0) {
        std::ostringstream path;
        path << "output/background_coupling_smooth_periodic_Nx" << Param::nx
             << "_Nu" << Param::Nv << "_Nuperp" << Param::Nmu << ".result";
        std::ofstream out(path.str().c_str());
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_smooth_periodic\n"
            << "production_kernel=1\nbeam_enabled=0\n"
            << "nx=" << Param::nx << " Nu=" << Param::Nv
            << " Nuperp=" << Param::Nmu << "\n";
        BackgroundCouplingTest::write_key_value(log, "dt",
                                                 BackgroundCouplingTest::stable_dt(sg));
        BackgroundCouplingTest::write_key_value(log, "full_L1", full.current_norms.l1);
        BackgroundCouplingTest::write_key_value(log, "full_L2", full.current_norms.l2);
        BackgroundCouplingTest::write_key_value(log, "full_Linf", full.current_norms.linf);
        BackgroundCouplingTest::write_key_value(log, "half_L1", half.current_norms.l1);
        BackgroundCouplingTest::write_key_value(log, "half_L2", half.current_norms.l2);
        BackgroundCouplingTest::write_key_value(log, "half_Linf", half.current_norms.linf);
        BackgroundCouplingTest::write_key_value(log, "current_difference_scale_L2",
                                                 current_scale_l2);
        BackgroundCouplingTest::write_key_value(log, "full_work_L1", full.work_norms.l1);
        BackgroundCouplingTest::write_key_value(log, "full_work_L2", full.work_norms.l2);
        BackgroundCouplingTest::write_key_value(log, "full_work_Linf", full.work_norms.linf);
        BackgroundCouplingTest::write_key_value(log, "half_work_L1", half.work_norms.l1);
        BackgroundCouplingTest::write_key_value(log, "half_work_L2", half.work_norms.l2);
        BackgroundCouplingTest::write_key_value(log, "half_work_Linf", half.work_norms.linf);
        BackgroundCouplingTest::write_key_value(log, "field_work_scale_L1", work_scale_l1);
        BackgroundCouplingTest::write_key_value(log, "field_work_scale_L2", work_scale_l2);
        BackgroundCouplingTest::write_key_value(log, "field_work_scale_Linf", work_scale_linf);
        BackgroundCouplingTest::write_key_value(log, "phase_shift_relative_L2",
                                                 phase_relative);
        log << "full_fct_active=" << full.fct_active << "\n"
            << "half_fct_active=" << half.fct_active << "\n"
            << "shifted_fct_active=" << shifted.fct_active << "\n"
            << "passes=" << local_ok << "\n";
        std::cout << "background_coupling_smooth_periodic_test result="
                  << path.str() << " passes=" << local_ok << "\n";
    }
    MPI_Finalize();
    return local_ok ? 0 : 1;
}
