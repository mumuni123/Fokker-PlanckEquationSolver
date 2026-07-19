#include "background_coupling_test_support.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct CaseResult {
    BackgroundCouplingTest::Norms current;
    BackgroundCouplingTest::Norms work;
    std::vector<double> fx_high;
    std::vector<double> fu_center;
    std::vector<double> gstar_je_center;
    unsigned long long x_high_state_hash;
    unsigned long long end_field_hash;
    double signed_work;
    int state_advanced;
    int operator_failed;
    int outputs_finite;
    int fct_active;
};

CaseResult run_case(const Species& background, const EMFields& fields,
                    const SpatialGrid& sg, int rank, int size, double dt,
                    bool energy_consistent_x_speed)
{
    BackgroundCouplingTest::BundleOptions options;
    // Both arms use the same raw high-order production operator.  Only the
    // speed passed into x-face reconstruction/donor selection is changed.
    options.fct_enabled = false;
    options.allow_finite_negative_debt = true;
    options.energy_consistent_x_high_velocity = energy_consistent_x_speed;
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                size, dt, options);
    CaseResult result = {};
    const bool arrays_valid = bundle.jn_high.size() ==
            static_cast<size_t>(sg.nx_local + 1) &&
        bundle.gstar_je_center.size() == static_cast<size_t>(sg.nx_local + 1);
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    result.current = {invalid, invalid, invalid};
    result.work = {invalid, invalid, invalid};
    result.signed_work = invalid;
    if (arrays_valid) {
        result.current = BackgroundCouplingTest::face_difference_norms(
            bundle.jn_high, bundle.gstar_je_center, sg);
        result.work = BackgroundCouplingTest::face_work_residual_norms(
            bundle.jn_high, bundle.gstar_je_center, fields, sg, dt);
        result.signed_work = BackgroundCouplingTest::global_signed_face_work_residual(
            bundle.jn_high, bundle.gstar_je_center, fields, sg, dt);
    }
    result.fx_high = bundle.fx_high;
    result.fu_center = bundle.fu_center;
    result.gstar_je_center = bundle.gstar_je_center;
    result.x_high_state_hash = bundle.x_high_state_hash;
    result.end_field_hash = bundle.end_field_hash;
    result.state_advanced = bundle.state_advanced;
    result.operator_failed = bundle.operator_failed;
    result.outputs_finite = bundle.outputs_finite;
    result.fct_active = bundle.fct_active;
    return result;
}

struct Sample {
    const char* name;
    double field_amplitude;
    double phase;
};

void write_norms(std::ostream& out, const std::string& prefix,
                 const BackgroundCouplingTest::Norms& norms)
{
    out << prefix << "_L1=" << norms.l1 << "\n"
        << prefix << "_L2=" << norms.l2 << "\n"
        << prefix << "_Linf=" << norms.linf << "\n";
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
    const double phase_shift = 7.0 * 2.0 * Const::pi /
        static_cast<double>(sg.nx_global);
    const Sample samples[] = {
        {"plus_phase_0", 1.0e3, 0.0},
        {"minus_phase_0", -1.0e3, 0.0},
        {"plus_phase_7", 1.0e3, phase_shift},
        {"minus_phase_7", -1.0e3, phase_shift}
    };
    const double dt = BackgroundCouplingTest::stable_dt(sg);
    int operator_valid = 1;
    int experiment_control_valid = 1;
    int common_speed_improves_all_cases = 1;
    double common_speed_worst_l2_ratio = 0.0;

    std::ofstream out;
    if (rank == 0)
        out.open("output/background_coupling_common_speed_flux.result");
    std::ostream& log = out ? out : std::cout;
    if (rank == 0) {
        log << std::scientific << std::setprecision(17)
            << "test=background_coupling_common_speed_flux\n"
            << "production_kernel=1\nbeam_enabled=0\nfct_enabled=0\n"
            << "sample_count=4\n"
            << "dt=" << dt << "\n";
    }

    for (size_t sample_id = 0;
         sample_id < sizeof(samples) / sizeof(samples[0]); ++sample_id) {
        Species background;
        EMFields fields;
        BackgroundCouplingTest::initialize_periodic_state(
            background, fields, sg, rank, size, 0.05, 0.08,
            samples[sample_id].field_amplitude, 3, samples[sample_id].phase);
        const CaseResult analytic = run_case(background, fields, sg, rank,
                                             size, dt, false);
        const CaseResult common = run_case(background, fields, sg, rank,
                                           size, dt, true);
        const double fx_difference =
            BackgroundCouplingTest::global_vector_difference_linf(
                analytic.fx_high, common.fx_high);
        const double fx_scale = std::max(
            BackgroundCouplingTest::global_vector_abs_linf(analytic.fx_high),
            BackgroundCouplingTest::global_vector_abs_linf(common.fx_high));
        const double fu_difference =
            BackgroundCouplingTest::global_vector_difference_linf(
                analytic.fu_center, common.fu_center);
        const double fu_scale = std::max(
            BackgroundCouplingTest::global_vector_abs_linf(analytic.fu_center),
            BackgroundCouplingTest::global_vector_abs_linf(common.fu_center));
        const double gstar_difference =
            BackgroundCouplingTest::global_vector_difference_linf(
                analytic.gstar_je_center, common.gstar_je_center);
        const double gstar_scale = std::max(
            BackgroundCouplingTest::global_vector_abs_linf(analytic.gstar_je_center),
            BackgroundCouplingTest::global_vector_abs_linf(common.gstar_je_center));
        const double fx_relative = fx_difference / std::max(1.0e-300, fx_scale);
        const double fu_relative = fu_difference / std::max(1.0e-300, fu_scale);
        const double gstar_relative = gstar_difference /
            std::max(1.0e-300, gstar_scale);
        const double current_ratio = common.current.l2 /
            std::max(1.0e-300, analytic.current.l2);
        const double work_ratio = common.work.l2 /
            std::max(1.0e-300, analytic.work.l2);
        int sample_operator_valid = analytic.state_advanced &&
            common.state_advanced && !analytic.operator_failed &&
            !common.operator_failed && analytic.outputs_finite &&
            common.outputs_finite && analytic.fct_active == 0 &&
            common.fct_active == 0 && std::isfinite(analytic.current.l2) &&
            std::isfinite(common.current.l2) && std::isfinite(analytic.work.l2) &&
            std::isfinite(common.work.l2);
        const double roundoff = 4096.0 * std::numeric_limits<double>::epsilon();
        const int state_hash_equal = analytic.x_high_state_hash ==
            common.x_high_state_hash;
        const int field_hash_equal = analytic.end_field_hash == common.end_field_hash;
        int only_x_high_branch_changed = fx_relative > roundoff &&
            fu_relative <= roundoff && gstar_relative <= roundoff &&
            state_hash_equal && field_hash_equal;
        MPI_Allreduce(MPI_IN_PLACE, &sample_operator_valid, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &only_x_high_branch_changed, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);
        operator_valid = operator_valid && sample_operator_valid;
        experiment_control_valid = experiment_control_valid &&
            only_x_high_branch_changed;
        common_speed_improves_all_cases = common_speed_improves_all_cases &&
            work_ratio < 1.0;
        common_speed_worst_l2_ratio = std::max(common_speed_worst_l2_ratio,
                                                work_ratio);
        if (rank == 0) {
            const std::string prefix = std::string("sample_") + samples[sample_id].name;
            log << prefix << "_field_amplitude=" << samples[sample_id].field_amplitude << "\n"
                << prefix << "_phase=" << samples[sample_id].phase << "\n";
            write_norms(log, prefix + "_analytic_current", analytic.current);
            write_norms(log, prefix + "_common_current", common.current);
            write_norms(log, prefix + "_analytic_work", analytic.work);
            write_norms(log, prefix + "_common_work", common.work);
            log << prefix << "_analytic_signed_work=" << analytic.signed_work << "\n"
                << prefix << "_common_signed_work=" << common.signed_work << "\n"
                << prefix << "_current_L2_ratio=" << current_ratio << "\n"
                << prefix << "_work_L2_ratio=" << work_ratio << "\n"
                << prefix << "_fx_high_difference_relative_Linf=" << fx_relative << "\n"
                << prefix << "_fu_center_difference_relative_Linf=" << fu_relative << "\n"
                << prefix << "_GstarJE_center_difference_relative_Linf=" << gstar_relative << "\n"
                << prefix << "_state_hash_equal=" << state_hash_equal << "\n"
                << prefix << "_field_hash_equal=" << field_hash_equal << "\n"
                << prefix << "_only_x_high_branch_changed=" << only_x_high_branch_changed << "\n"
                << prefix << "_operator_valid=" << sample_operator_valid << "\n";
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &operator_valid, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &experiment_control_valid, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &common_speed_improves_all_cases, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &common_speed_worst_l2_ratio, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    const int passes = operator_valid && experiment_control_valid;
    const int hypothesis_pass = passes && common_speed_improves_all_cases;
    if (rank == 0) {
        log << "operator_valid=" << operator_valid << "\n"
            << "experiment_control_valid=" << experiment_control_valid << "\n"
            << "common_speed_improves_all_cases=" << common_speed_improves_all_cases << "\n"
            << "common_speed_worst_L2_ratio=" << common_speed_worst_l2_ratio << "\n"
            << "hypothesis_pass=" << hypothesis_pass << "\n"
            << "passes=" << passes << "\n";
        std::cout << "background_coupling_common_speed_flux_test passes="
                  << passes << " hypothesis_pass=" << hypothesis_pass << "\n";
    }
    MPI_Finalize();
    return passes ? 0 : 1;
}
