#include "background_coupling_test_support.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#ifndef FP_DIRECTIONAL_GRID_KIND
#define FP_DIRECTIONAL_GRID_KIND 0
#endif

namespace {

const char* refinement_kind()
{
#if FP_DIRECTIONAL_GRID_KIND == 1
    return "x";
#elif FP_DIRECTIONAL_GRID_KIND == 2
    return "u_parallel";
#elif FP_DIRECTIONAL_GRID_KIND == 3
    return "u_perp";
#else
    return "invalid";
#endif
}

double fixed_directional_dt(const SpatialGrid&)
{
#if FP_DIRECTIONAL_GRID_KIND == 1
    // Every x-resolution uses the finest x mesh time step.  This isolates
    // the x truncation trend from a simultaneous dt refinement.
    return 0.025 * (Param::Lx / 4000.0) / Const::c;
#else
    // Nx is fixed for both cylindrical velocity refinements.  The same
    // spatial-CFL-safe dt is therefore used at all u_parallel/u_perp levels;
    // the finest velocity grid still passes through the production CFL gate.
    return 0.025 * (Param::Lx / 2000.0) / Const::c;
#endif
}

const char* dt_policy()
{
#if FP_DIRECTIONAL_GRID_KIND == 1
    return "fixed_from_finest_x";
#else
    return "fixed_within_velocity_refinement";
#endif
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
    Species background;
    EMFields fields;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.04, 0.06, 1.0e3, 3, 0.0);

    // The directional study deliberately disables FCT for every resolution.
    // This keeps the comparison on the same high-order operator branch and
    // prevents resolution-dependent limiter activation from masking the
    // truncation trend being measured.
    BackgroundCouplingTest::BundleOptions options;
    options.fct_enabled = false;
    options.allow_finite_negative_debt = true;
    const double dt = fixed_directional_dt(sg);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                size, dt, options);
    const bool arrays_have_expected_size =
        bundle.jn_high.size() == static_cast<size_t>(sg.nx_local + 1) &&
        bundle.gstar_je_center.size() == static_cast<size_t>(sg.nx_local + 1);
    const int local_valid = bundle.state_advanced && !bundle.operator_failed &&
        bundle.outputs_finite && arrays_have_expected_size &&
        bundle.fct_active == 0;
    BackgroundCouplingTest::Norms current = {};
    BackgroundCouplingTest::Norms work = {};
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    current.l1 = current.l2 = current.linf = invalid;
    work.l1 = work.l2 = work.linf = invalid;
    if (local_valid) {
        current = BackgroundCouplingTest::face_difference_norms(
            bundle.jn_high, bundle.gstar_je_center, sg);
        work = BackgroundCouplingTest::face_work_residual_norms(
            bundle.jn_high, bundle.gstar_je_center, fields, sg, dt);
    }
    const int local_ok = local_valid && std::isfinite(current.l2) &&
        std::isfinite(work.l2) && current.l2 > 0.0 && work.l2 > 0.0;
    int global_ok = local_ok;
    MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    if (rank == 0) {
        std::ostringstream path;
        path << "output/background_coupling_directional_" << refinement_kind()
             << "_Nx" << Param::nx << "_Nu" << Param::Nv << "_Nuperp"
             << Param::Nmu << ".result";
        std::ofstream out(path.str().c_str());
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_directional_grid\n"
            << "production_kernel=1\nbeam_enabled=0\n"
            << "refinement_kind=" << refinement_kind() << "\n"
            << "nx=" << Param::nx << "\nNu=" << Param::Nv
            << "\nNuperp=" << Param::Nmu << "\n"
            << "dt=" << dt << "\n"
            << "dt_policy=" << dt_policy() << "\n"
            << "dt_reference_nx=" << (FP_DIRECTIONAL_GRID_KIND == 1 ? 4000 : 2000)
            << "\n"
            << "fct_enabled=0\n"
            << "fct_active=" << bundle.fct_active << "\n"
            << "state_advanced=" << bundle.state_advanced << "\n"
            << "operator_failed=" << bundle.operator_failed << "\n"
            << "outputs_finite=" << bundle.outputs_finite << "\n"
            << "failure_reason=" << bundle.failure_reason << "\n"
            << "failure_iteration=" << bundle.failure_iteration << "\n"
            << "failure_substep=" << bundle.failure_substep << "\n"
            << "low_order_candidate_min=" << bundle.low_order_candidate_min << "\n"
            << "low_order_negative_count=" << bundle.low_order_negative_count << "\n"
            << "low_order_negative_mass=" << bundle.low_order_negative_mass << "\n"
            << "final_candidate_min=" << bundle.final_candidate_min << "\n"
            << "limiter_active_fraction=" << bundle.limiter_active_fraction << "\n"
            << "limiter_min_alpha=" << bundle.limiter_min_alpha << "\n"
            << "current_L1=" << current.l1 << "\n"
            << "current_L2=" << current.l2 << "\n"
            << "current_Linf=" << current.linf << "\n"
            << "work_L1=" << work.l1 << "\n"
            << "work_L2=" << work.l2 << "\n"
            << "work_Linf=" << work.linf << "\n"
            << "audit_valid=" << global_ok << "\n"
            << "passes=" << global_ok << "\n";
        std::cout << "background_coupling_directional_grid_test result="
                  << path.str() << " passes=" << global_ok << "\n";
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
