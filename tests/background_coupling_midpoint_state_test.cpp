#include "background_coupling_test_support.h"

#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>

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
        background, fields, sg, rank, size, 0.06, 0.09, 1.5e3, 5, 0.3);
    Species guess_np1;
    EMFields fields_end_guess;
    BackgroundCouplingTest::make_distinct_endpoint_guess(
        background, fields, guess_np1, fields_end_guess, sg, rank, size);
    BackgroundCouplingTest::BundleOptions options;
    options.guess_np1 = &guess_np1;
    options.fields_end_guess = &fields_end_guess;
    // This is the existing test-only conservative anti-flux injection.  It
    // forces the production FCT path without changing production defaults.
    options.controlled_fct_injection = true;
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(
            background, fields, sg, rank, size,
            BackgroundCouplingTest::stable_dt(sg), options);

    const int high_state_shared =
        bundle.x_high_state_hash == bundle.u_high_state_hash ? 1 : 0;
    const int high_field_shared =
        bundle.x_high_field_hash == bundle.u_high_field_hash ? 1 : 0;
    const int low_state_shared =
        bundle.x_low_state_hash == bundle.u_low_state_hash ? 1 : 0;
    const int high_midpoint_labels = bundle.x_high_time_layer == 1 &&
        bundle.u_high_time_layer == 1;
    const int low_start_labels = bundle.x_low_time_layer == 0 &&
        bundle.u_low_time_layer == 0;
    const int low_u_midpoint_field = bundle.u_low_field_hash != 0ULL &&
        bundle.u_low_field_hash == bundle.u_high_field_hash;
    const int distinct_high_state = bundle.x_high_state_hash !=
        bundle.x_low_state_hash;
    const int distinct_midpoint_field = bundle.u_high_field_hash != 0ULL &&
        bundle.start_field_hash != bundle.end_field_hash &&
        bundle.u_high_field_hash != bundle.start_field_hash &&
        bundle.u_high_field_hash != bundle.end_field_hash;
    const BackgroundCouplingTest::Norms jn_final_minus_high =
        BackgroundCouplingTest::face_difference_norms(
            bundle.jn_final, bundle.jn_high, sg);
    const BackgroundCouplingTest::Norms je_final_minus_high =
        BackgroundCouplingTest::face_difference_norms(
            bundle.gstar_je_final, bundle.gstar_je_high, sg);
    const double fx_final_minus_high_linf =
        BackgroundCouplingTest::global_vector_difference_linf(
            bundle.fx_final, bundle.fx_high);
    const double fu_final_minus_high_linf =
        BackgroundCouplingTest::global_vector_difference_linf(
            bundle.fu_final, bundle.fu_high);
    const double moment_roundoff = 8192.0 * std::numeric_limits<double>::epsilon();
    const int final_flux_to_current_identity =
        bundle.final_flux_current_moment_audit_valid &&
        bundle.final_flux_current_moment_audit_finite &&
        bundle.final_flux_to_jn_linf <= moment_roundoff *
            std::max(1.0, bundle.final_flux_to_jn_scale) &&
        bundle.final_flux_to_je_linf <= moment_roundoff *
            std::max(1.0, bundle.final_flux_to_je_scale) &&
        bundle.final_flux_to_gstar_je_linf <= moment_roundoff *
            std::max(1.0, bundle.final_flux_to_gstar_je_scale);
    int substep_history_shared =
        bundle.x_low_state_hash_history.size() ==
        bundle.u_low_state_hash_history.size() &&
        bundle.x_high_state_hash_history.size() ==
        bundle.u_high_state_hash_history.size() &&
        bundle.u_high_state_hash_history.size() ==
        bundle.u_field_hash_history.size();
    for (size_t sub = 0; sub < bundle.x_low_state_hash_history.size() &&
         substep_history_shared; ++sub) {
        substep_history_shared =
            bundle.x_low_state_hash_history[sub] ==
                bundle.u_low_state_hash_history[sub] &&
            bundle.x_high_state_hash_history[sub] ==
                bundle.u_high_state_hash_history[sub];
    }
    int local_ok = bundle.state_advanced && bundle.finite &&
        high_state_shared && high_field_shared && low_state_shared &&
        high_midpoint_labels && low_start_labels && low_u_midpoint_field &&
        distinct_high_state && distinct_midpoint_field &&
        bundle.fct_active && bundle.limiter_min_alpha < 1.0 - 1.0e-14 &&
        fx_final_minus_high_linf > 0.0 && jn_final_minus_high.linf > 0.0 &&
        final_flux_to_current_identity;
    local_ok = local_ok && substep_history_shared;
    MPI_Allreduce(MPI_IN_PLACE, &local_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream out("output/background_coupling_midpoint_state.result");
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_midpoint_state\n"
            << "production_kernel=1\nbeam_enabled=0\n"
            << "x_low_state_hash=" << bundle.x_low_state_hash << "\n"
            << "u_low_state_hash=" << bundle.u_low_state_hash << "\n"
            << "x_high_state_hash=" << bundle.x_high_state_hash << "\n"
            << "u_high_state_hash=" << bundle.u_high_state_hash << "\n"
            << "x_low_field_hash=" << bundle.x_low_field_hash << "\n"
            << "u_low_field_hash=" << bundle.u_low_field_hash << "\n"
            << "x_high_field_hash=" << bundle.x_high_field_hash << "\n"
            << "u_high_field_hash=" << bundle.u_high_field_hash << "\n"
            << "start_field_hash=" << bundle.start_field_hash << "\n"
            << "end_field_hash=" << bundle.end_field_hash << "\n"
            << "x_low_time_layer=" << bundle.x_low_time_layer << "\n"
            << "u_low_time_layer=" << bundle.u_low_time_layer << "\n"
            << "x_high_time_layer=" << bundle.x_high_time_layer << "\n"
            << "u_high_time_layer=" << bundle.u_high_time_layer << "\n"
            << "high_state_shared=" << high_state_shared << "\n"
            << "high_field_shared=" << high_field_shared << "\n"
            << "low_state_shared=" << low_state_shared << "\n"
            << "low_u_midpoint_field=" << low_u_midpoint_field << "\n"
            << "distinct_high_state=" << distinct_high_state << "\n"
            << "distinct_midpoint_field=" << distinct_midpoint_field << "\n"
            << "limiter_active_fraction=" << bundle.limiter_active_fraction << "\n"
            << "limiter_min_alpha=" << bundle.limiter_min_alpha << "\n"
            << "fx_final_minus_high_Linf=" << fx_final_minus_high_linf << "\n"
            << "JN_final_minus_high_Linf=" << jn_final_minus_high.linf << "\n"
            << "GstarJE_final_minus_high_Linf=" << je_final_minus_high.linf << "\n"
            << "fu_final_minus_high_Linf=" << fu_final_minus_high_linf << "\n"
            << "final_flux_current_moment_audit_valid="
            << bundle.final_flux_current_moment_audit_valid << "\n"
            << "final_flux_current_moment_audit_finite="
            << bundle.final_flux_current_moment_audit_finite << "\n"
            << "final_flux_to_JN_Linf=" << bundle.final_flux_to_jn_linf << "\n"
            << "final_flux_to_JN_scale=" << bundle.final_flux_to_jn_scale << "\n"
            << "final_flux_to_JE_Linf=" << bundle.final_flux_to_je_linf << "\n"
            << "final_flux_to_JE_scale=" << bundle.final_flux_to_je_scale << "\n"
            << "final_flux_to_GstarJE_Linf="
            << bundle.final_flux_to_gstar_je_linf << "\n"
            << "final_flux_to_GstarJE_scale="
            << bundle.final_flux_to_gstar_je_scale << "\n"
            << "final_flux_current_moment_roundoff_bound=" << moment_roundoff << "\n"
            << "final_flux_to_current_identity="
            << final_flux_to_current_identity << "\n"
            << "substep_history_count="
            << bundle.x_low_state_hash_history.size() << "\n"
            << "substep_history_shared=" << substep_history_shared << "\n"
            << "fct_active=" << bundle.fct_active << "\n"
            << "passes=" << local_ok << "\n";
        for (size_t sub = 0; sub < bundle.x_low_state_hash_history.size(); ++sub) {
            log << "substep=" << sub
                << " x_low_state_hash=" << bundle.x_low_state_hash_history[sub]
                << " u_low_state_hash=" << bundle.u_low_state_hash_history[sub]
                << " x_high_state_hash=" << bundle.x_high_state_hash_history[sub]
                << " u_high_state_hash=" << bundle.u_high_state_hash_history[sub]
                << " u_field_hash=" << bundle.u_field_hash_history[sub]
                << "\n";
        }
        std::cout << "background_coupling_midpoint_state_test passes="
                  << local_ok << "\n";
    }
    MPI_Finalize();
    return local_ok ? 0 : 1;
}
