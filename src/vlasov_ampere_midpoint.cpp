#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <mpi.h>

namespace {

void fill_periodic_ghosts(Species& sp, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    for (int g = 0; g < ng; ++g) {
        const int left_source = ng + nxl - ng + g;
        const int right_source = ng + g;
        std::copy(sp.f.begin() + static_cast<size_t>(left_source) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(left_source + 1) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(g) * Param::Nvmu);
        std::copy(sp.f.begin() + static_cast<size_t>(right_source) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(right_source + 1) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(ng + nxl + g) * Param::Nvmu);
    }
}

}

VlasovAmpereMidpointSolver::VlasovAmpereMidpointSolver()
    : step_diagnostics_enabled_(false), accepted_energy_audit_enabled_(false),
      beam_enabled_(true),
      low_order_only_(false), nonuniform_high_order_enabled_(false),
      fct_enabled_(true), background_coupling_mode_(DUAL_U_COUPLING),
      face_pairing_mode_(FACE_PAIRING_CELL_BASELINE),
      face_pairing_sigma_cutoff_(1.0e-8),
      face_pairing_lambda_(1.0e-3), face_pairing_eta_(1.0e-8),
      face_pairing_trust_fraction_(0.1),
      face_pairing_correction_trust_fraction_(1.0),
      face_pairing_energy_pair_tolerance_(1.0e-8),
      face_pairing_energy_residual_fraction_(1.0),
      face_pairing_mass_relative_tolerance_(1.0e-10),
      face_pairing_f_residual_growth_tolerance_(0.1),
      capture_midpoint_input_(false),
      legacy_boundary_upwind_high_candidate_for_test_(false),
      energy_consistent_x_high_velocity_for_test_(false),
      max_midpoint_iterations_(20),
      midpoint_initial_guess_mode_(MIDPOINT_INITIAL_GUESS_NONE),
      midpoint_acceleration_mode_(MIDPOINT_ACCELERATION_NONE),
      anderson_depth_(3), acceleration_start_iter_(3),
      acceleration_accept_ratio_(0.95), acceleration_max_coefficient_(2.0),
      midpoint_iteration_trace_for_test_(false),
      progress_trace_start_fs_(-1.0), progress_trace_end_fs_(-1.0),
      fct_activation_audit_enabled_(false),
      controlled_fct_flux_injection_enabled_(false),
      controlled_u_fct_flux_injection_enabled_(false),
      allow_finite_negative_debt_for_test_(false),
      core_macro_debt_consecutive_steps_(0)
{}

void VlasovAmpereMidpointSolver::reset_current_diag(CurrentDiagnostics& diag) const
{
    diag = CurrentDiagnostics();
}

void VlasovAmpereMidpointSolver::reset_result(Result& result) const
{
    result = Result();
    result.operator_evaluations = 0;
    result.midpoint_initial_guess_mode =
        static_cast<int>(midpoint_initial_guess_mode_);
    result.midpoint_predictor_used = 0;
    result.midpoint_predictor_history_depth =
        midpoint_field_predictor_.history_depth();
    result.midpoint_acceleration_mode = MIDPOINT_ACCELERATION_NONE;
    result.acceleration_attempts = 0;
    result.acceleration_accepted = 0;
    result.acceleration_fallback_evaluations = 0;
    result.acceleration_rejected_residual = 0;
    result.acceleration_rejected_nonfinite = 0;
    result.acceleration_rejected_hard_failure = 0;
    result.acceleration_rejected_coefficient = 0;
    result.acceleration_history_resets = 0;
    result.max_residual_E = 0.0;
    result.max_residual_J_bkg = 0.0;
    result.max_residual_f = 0.0;
    result.final_dual_u_valid = 1;
    result.final_dual_u_target_linf = 0.0;
    result.final_dual_u_residual_before_linf = 0.0;
    result.final_dual_u_residual_after_linf = 0.0;
    result.final_dual_u_minimum_scale = 1.0;
    result.final_dual_u_correction_l2 = 0.0;
    result.final_dual_u_correction_linf = 0.0;
    result.final_dual_u_candidate_min =
        std::numeric_limits<double>::infinity();
    result.final_dual_u_corrected_cell_count = 0;
    result.final_dual_u_limited_cell_count = 0;
    result.final_dual_u_unresolved_cell_count = 0;
    result.face_pairing_attempted = 0;
    result.face_pairing_accepted = 0;
    result.face_pairing_fallback_to_cell_baseline = 0;
    result.face_pairing_solver_converged = 0;
    result.face_pairing_iterations = 0;
    result.face_pairing_unresolved_mode_count = 0;
    result.face_pairing_residual_before = 0.0;
    result.face_pairing_residual_after = 0.0;
    result.face_pairing_core_residual_before = 0.0;
    result.face_pairing_core_residual_after = 0.0;
    result.face_pairing_unresolved_mode_l2 = 0.0;
    result.face_pairing_correction_l2 = 0.0;
    result.face_pairing_correction_linf = 0.0;
    result.face_pairing_delta_ke = 0.0;
    result.face_pairing_delta_work = 0.0;
    result.face_pairing_candidate_min =
        std::numeric_limits<double>::infinity();
    result.face_pairing_mass_error = 0.0;
    result.face_pairing_mass_relative_error = 0.0;
    result.face_pairing_cell_mass_error_linf = 0.0;
    result.face_pairing_cell_mass_relative_linf = 0.0;
    result.face_pairing_energy_pair_error = 0.0;
    result.face_pairing_energy_pair_relative = 0.0;
    result.face_pairing_energy_residual_scale = 0.0;
    result.face_pairing_energy_residual_ratio = 0.0;
    result.face_pairing_correction_trust_limit = 0.0;
    result.face_pairing_correction_trust_ratio = 0.0;
    result.face_pairing_f_residual_relative_growth = 0.0;
    result.face_pairing_capacity_active_cells = 0;
    result.face_pairing_trust_region_active_cells = 0;
    result.face_pairing_candidate_valid = 0;
    result.face_pairing_candidate_residual_after =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_candidate_core_residual_after =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_candidate_delta_ke =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_candidate_delta_work =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_candidate_mass_error =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_candidate_min_before_fallback =
        std::numeric_limits<double>::quiet_NaN();
    result.face_pairing_requested_correction_l2 = 0.0;
    result.face_pairing_requested_correction_linf = 0.0;
    result.face_pairing_applied_correction_l2 = 0.0;
    result.face_pairing_applied_correction_linf = 0.0;
    result.face_pairing_nonzero_capacity_cells = 0;
    result.face_pairing_bound_saturated_cells = 0;
    result.face_pairing_objective_residual = 0.0;
    result.face_pairing_objective_smoothness = 0.0;
    result.face_pairing_objective_amplitude = 0.0;
    result.face_pairing_objective_total = 0.0;
    result.face_pairing_rejection_mask = 0u;
    result.face_pairing_pass_solver = 0;
    result.face_pairing_pass_apply = 0;
    result.face_pairing_pass_global_residual = 0;
    result.face_pairing_pass_core_residual = 0;
    result.face_pairing_pass_correction_trust = 0;
    result.face_pairing_pass_delta_ke = 0;
    result.face_pairing_pass_delta_work = 0;
    result.face_pairing_pass_candidate_min = 0;
    result.face_pairing_pass_mass = 0;
    result.face_pairing_pass_f_residual = 0;
    result.face_pairing_pass_energy_pair = 0;
    result.face_pairing_pass_energy_residual_scale = 0;
    result.fct_macro_budget_valid = 0;
    for (size_t i = 0; i < result.fct_macro_budget_x.size(); ++i) {
        FctMacroBudget* budgets[2] = {
            &result.fct_macro_budget_x[i], &result.fct_macro_budget_u[i]};
        for (int direction = 0; direction < 2; ++direction) {
            budgets[direction]->face_count = 0;
            budgets[direction]->active_face_count = 0;
            budgets[direction]->min_alpha = 1.0;
            budgets[direction]->delta_n = 0.0;
            budgets[direction]->delta_j = 0.0;
            budgets[direction]->delta_k = 0.0;
            budgets[direction]->e_dot_j = 0.0;
            budgets[direction]->r_fct_e = 0.0;
        }
    }
    result.limiter_min_alpha = 1.0;
    result.x_limiter_min_alpha = 1.0;
    result.u_limiter_min_alpha = 1.0;
    result.limiter_min_alpha_core = 1.0;
    result.limiter_min_alpha_boundary = 1.0;
    result.mu_low_u_alpha_min = 1.0;
    result.mu_low_u_alpha_min_boundary = 1.0;
    result.mu_low_u_alpha_min_core = 1.0;
    result.alpha_core_min = 1.0;
    result.alpha_boundary_min = 1.0;
    result.alpha_tail_min = 1.0;
    result.u_force_alpha_min = 1.0;
    result.low_u_average_subcycles = 1.0;
    result.low_u_max_subcycles = 1;
    result.boundary_force_nsub_max = 1;
    result.x_low_input_min_f = std::numeric_limits<double>::infinity();
    result.x_low_output_min_f = std::numeric_limits<double>::infinity();
    result.x_final_min_f = std::numeric_limits<double>::infinity();
    result.failure_low_min = std::numeric_limits<double>::infinity();
    result.failure_final_min = std::numeric_limits<double>::infinity();
    result.fct_high_candidate_min = std::numeric_limits<double>::infinity();
    result.fct_final_scratch_min = std::numeric_limits<double>::infinity();
    result.low_order_candidate_min = std::numeric_limits<double>::infinity();
    result.low_order_negative_count = 0;
    result.low_order_negative_mass = 0.0;
    result.low_order_roundoff_zeroed_count = 0;
    result.low_order_roundoff_zeroed_mass = 0.0;
    result.fct_roundoff_zeroed_count = 0;
    result.fct_roundoff_zeroed_mass = 0.0;
    result.accepted_energy_ledger = AcceptedEnergyLedger();
    result.accepted_energy_ledger.valid = 0;
    result.accepted_energy_ledger.strict_accepted = 0;
    result.accepted_energy_ledger.soft_accepted = 0;
    result.accepted_energy_ledger.transport_substeps = 0;
    result.beam_continuity_valid = 1;
    result.stage5_r_couple_centered = 0.0;
    result.stage5_r_couple_upwind_stabilization = 0.0;
    result.stage5_r_couple_fct_stabilization = 0.0;
    result.coupling_rj_global_sum = 0.0;
    result.coupling_rk_global_sum = 0.0;
    result.coupling_face_work_jn_global_sum = 0.0;
    result.coupling_rj_reconstruction_error = 0.0;
    result.coupling_rk_reconstruction_error = 0.0;
    result.coupling_face_work_jn_reconstruction_error = 0.0;
    result.coupling_beam_front_ix = -1;
    for (size_t region = 0; region < result.coupling_regions.size(); ++region) {
        result.coupling_regions[region].max_abs_jn_minus_gstar_je_face = -1;
    }
    result.fct_high_low_identity_worst_ix = -1;
    result.fct_high_low_identity_worst_iv = -1;
    result.fct_high_low_identity_worst_imu = -1;
    result.fct_donor_worst_ix = -1;
    result.fct_donor_worst_iv = -1;
    result.fct_donor_worst_imu = -1;
    result.fct_donor_beta_min = 1.0;
    result.failure_worst_ix = -1;
    result.failure_worst_iv = -1;
    result.failure_worst_imu = -1;
    result.failure_worst_is_core = -1;
    result.failure_worst_is_tail = -1;
    result.failure_low_work_input_min =
        std::numeric_limits<double>::infinity();
    result.accepted_transport.min_mass =
        std::numeric_limits<double>::infinity();
    result.accepted_transport.ix = -1;
    result.accepted_transport.iv = -1;
    result.accepted_transport.imu = -1;
    for (int f = 0; f < 4; ++f) {
        result.accepted_transport.alpha[f] = 1.0;
    }
    for (int d = 0; d < 3; ++d) {
        result.flux_pos[d].min_f_before = std::numeric_limits<double>::infinity();
        result.flux_pos[d].min_f_low = std::numeric_limits<double>::infinity();
        result.flux_pos[d].min_f_final = std::numeric_limits<double>::infinity();
        result.flux_pos[d].alpha_min = 1.0;
    }
}

void VlasovAmpereMidpointSolver::exchange_ghosts_x_persistent(
    Species& sp, const SpatialGrid& sg, int mpi_rank, int mpi_size) const
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t count = static_cast<size_t>(ng) * Param::Nvmu;
    if (mpi_size == 1) {
        fill_periodic_ghosts(sp, sg);
        return;
    }
    if (ghost_send_left_.size() != count) {
        ghost_send_left_.resize(count);
        ghost_send_right_.resize(count);
        ghost_recv_left_.resize(count);
        ghost_recv_right_.resize(count);
    }
    std::memcpy(ghost_send_left_.data(),
                sp.f.data() + static_cast<size_t>(ng) * Param::Nvmu,
                count * sizeof(double));
    std::memcpy(ghost_send_right_.data(),
                sp.f.data() + static_cast<size_t>(ng + nxl - ng) * Param::Nvmu,
                count * sizeof(double));
    const int left = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right = (mpi_rank + 1) % mpi_size;
    MPI_Request reqs[4];
    MPI_Isend(ghost_send_left_.data(), static_cast<int>(count), MPI_DOUBLE,
              left, 501, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(ghost_recv_left_.data(), static_cast<int>(count), MPI_DOUBLE,
              left, 502, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(ghost_send_right_.data(), static_cast<int>(count), MPI_DOUBLE,
              right, 502, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(ghost_recv_right_.data(), static_cast<int>(count), MPI_DOUBLE,
              right, 501, MPI_COMM_WORLD, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
    std::memcpy(sp.f.data(), ghost_recv_left_.data(), count * sizeof(double));
    std::memcpy(sp.f.data() + static_cast<size_t>(ng + nxl) * Param::Nvmu,
                ghost_recv_right_.data(), count * sizeof(double));
}

void VlasovAmpereMidpointSolver::set_midpoint_field(
    EMFields& fields_mid, const EMFields& fields_n,
    const std::vector<double>& ex_mid, const SpatialGrid& sg,
    int mpi_rank, int mpi_size) const
{
    fields_mid = fields_n;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        fields_mid.Ex_face[static_cast<size_t>(iface)] = ex_mid[static_cast<size_t>(iface)];
    }
    if (sg.nx_local > 0) fields_mid.Ex_face[static_cast<size_t>(sg.nx_local)] =
        fields_mid.Ex_face[0];
    fields_mid.sync_cell_ex_from_faces(mpi_rank, mpi_size);
}

double VlasovAmpereMidpointSolver::integrate_face_work(
    const std::vector<double>& current_face, const EMFields& fields_mid,
    const SpatialGrid& sg, double dt) const
{
    double work = 0.0;
    const int n = std::min(sg.nx_local, static_cast<int>(std::min(
        current_face.size(), fields_mid.Ex_face.size())));
    for (int iface = 0; iface < n; ++iface) {
        work += current_face[static_cast<size_t>(iface)] *
                fields_mid.Ex_face[static_cast<size_t>(iface)] * sg.dx;
    }
    return dt * work;
}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_background_and_fields(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size)
{
    if (!bkg_n.cylindrical_mass_representation) {
        Result failed;
        reset_result(failed);
        failed.failed = true;
        return failed;
    }
    if (!low_order_only_ && !bkg_n.cgrid.is_uniform() &&
        !nonuniform_high_order_enabled_) {
        Result failed;
        reset_result(failed);
        failed.failed = true;
        failed.state_advanced = 0;
        failed.failure_reason = 12;
        if (mpi_rank == 0) {
            std::cerr << "VlasovAmpereMidpointSolver configuration error: "
                      << "nonuniform cylindrical velocity grid requires "
                      << "set_nonuniform_high_order_enabled(true) when "
                      << "low_order_only=false\n";
        }
        return failed;
    }
    std::vector<double> owned_field(
        fields_n.Ex_face.begin(),
        fields_n.Ex_face.begin() +
            std::min(fields_n.Ex_face.size(),
                     static_cast<size_t>(std::max(0, sg.nx_local))));
    std::vector<double> predicted_e_end;
    const bool predictor_used = midpoint_field_predictor_.propose(
        owned_field, dt, predicted_e_end);
    Result result = evaluate_production_midpoint_operator(
        bkg_n, beam_n, fields_n, sg, dt, time, mpi_rank, mpi_size, 1,
        0, 0, 0, 0, false,
        predictor_used ? &predicted_e_end : 0);
    result.midpoint_predictor_used = predictor_used ? 1 : 0;
    if (result.state_advanced && !result.failed) {
        const bool core_macro_debt =
            result.neg_mass_core_fraction > 1.0e-6 ||
            result.neg_energy_core_fraction > 1.0e-6 ||
            result.neg_current_core_fraction > 1.0e-6;
        core_macro_debt_consecutive_steps_ = core_macro_debt
            ? core_macro_debt_consecutive_steps_ + 1 : 0;
        if (core_macro_debt_consecutive_steps_ >= 3) {
            result.failed = true;
            result.state_advanced = 0;
        }
    }
    const bool strict_accepted =
        result.state_advanced && !result.failed && result.converged &&
        !result.soft_unconverged;
    if (strict_accepted) {
        midpoint_field_predictor_.commit_strict(
            result.fields_np1.Ex_face,
            static_cast<size_t>(std::max(0, sg.nx_local)), dt);
    } else {
        midpoint_field_predictor_.clear();
    }
    result.midpoint_predictor_history_depth =
        midpoint_field_predictor_.history_depth();
    return result;
}

VlasovAmpereMidpointSolver::MidpointOperatorEvaluation
VlasovAmpereMidpointSolver::evaluate_fixed_midpoint_operator(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const Species& guess_np1, const EMFields& fields_end_guess,
    const std::vector<double>& fixed_j_beam_face_mid,
    const CouplingRegionLayout& coupling_layout,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size) const
{
    return evaluate_production_midpoint_operator(
        bkg_n, beam_n, fields_n, sg, dt, time, mpi_rank, mpi_size, 1,
        &guess_np1, &fields_end_guess, &fixed_j_beam_face_mid,
        &coupling_layout, true);
}

VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle
VlasovAmpereMidpointSolver::evaluate_background_coupling_flux_bundle(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const Species& guess_np1, const EMFields& fields_end_guess,
    const std::vector<double>& fixed_j_beam_face_mid,
    const CouplingRegionLayout& coupling_layout,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size) const
{
    const MidpointOperatorEvaluation evaluation = evaluate_fixed_midpoint_operator(
        bkg_n, beam_n, fields_n, guess_np1, fields_end_guess,
        fixed_j_beam_face_mid, coupling_layout, sg, dt, time, mpi_rank,
        mpi_size);
    BackgroundCouplingFluxBundle bundle;
    bundle.fx_low = evaluation.low_x_flux;
    bundle.fx_high = evaluation.high_x_flux;
    bundle.fx_final = evaluation.final_x_flux;
    bundle.fu_low = evaluation.low_u_flux;
    bundle.fu_center = evaluation.center_u_flux;
    bundle.fu_high = evaluation.high_u_flux;
    bundle.fu_final = evaluation.final_u_flux;
    bundle.fu_fct_limited = evaluation.fct_limited_u_flux;
    bundle.cu_low = evaluation.low_u_coefficient;
    bundle.cu_high = evaluation.high_u_coefficient;
    bundle.cu_final = evaluation.final_u_coefficient;
    bundle.cu_fct_limited = evaluation.fct_limited_u_coefficient;
    bundle.donor_beta = evaluation.fct_donor_beta;
    bundle.donor_low_mass = evaluation.fct_donor_low_mass;
    bundle.donor_limited_outflow = evaluation.fct_donor_limited_outflow;
    bundle.cu_center = evaluation.center_u_coefficient;
    bundle.cu_legacy_center = evaluation.legacy_center_u_coefficient;
    bundle.dual_target_jn_cell = evaluation.dual_target_jn_cell;
    bundle.dual_target_jn_replay_cell =
        evaluation.dual_target_jn_replay_cell;
    bundle.dual_legacy_je_cell = evaluation.dual_legacy_je_cell;
    bundle.dual_je_cell = evaluation.dual_je_cell;
    bundle.cu_reconstruction_mass = evaluation.center_u_reconstruction_mass;
    bundle.background_coupling_mode = evaluation.background_coupling_mode;
    bundle.dual_u_operator_valid = evaluation.dual_u_operator_valid;
    bundle.dual_u_target_replay_linf = evaluation.dual_u_target_replay_linf;
    bundle.dual_u_target_replay_scale = evaluation.dual_u_target_replay_scale;
    bundle.dual_u_legacy_operator_replay_linf =
        evaluation.dual_u_legacy_operator_replay_linf;
    bundle.dual_u_legacy_operator_replay_scale =
        evaluation.dual_u_legacy_operator_replay_scale;
    bundle.dual_u_legacy_current_linf =
        evaluation.dual_u_legacy_current_linf;
    bundle.dual_u_current_linf = evaluation.dual_u_current_linf;
    bundle.dual_u_correction_l2 = evaluation.dual_u_correction_l2;
    bundle.dual_u_correction_linf = evaluation.dual_u_correction_linf;
    bundle.dual_u_corrected_cell_count =
        evaluation.dual_u_corrected_cell_count;
    bundle.final_dual_u_valid = evaluation.final_dual_u_valid;
    bundle.final_dual_u_target_linf =
        evaluation.final_dual_u_target_linf;
    bundle.final_dual_u_residual_before_linf =
        evaluation.final_dual_u_residual_before_linf;
    bundle.final_dual_u_residual_after_linf =
        evaluation.final_dual_u_residual_after_linf;
    bundle.final_dual_u_minimum_scale =
        evaluation.final_dual_u_minimum_scale;
    bundle.final_dual_u_correction_l2 =
        evaluation.final_dual_u_correction_l2;
    bundle.final_dual_u_correction_linf =
        evaluation.final_dual_u_correction_linf;
    bundle.final_dual_u_candidate_min =
        evaluation.final_dual_u_candidate_min;
    bundle.final_dual_u_corrected_cell_count =
        evaluation.final_dual_u_corrected_cell_count;
    bundle.final_dual_u_limited_cell_count =
        evaluation.final_dual_u_limited_cell_count;
    bundle.final_dual_u_unresolved_cell_count =
        evaluation.final_dual_u_unresolved_cell_count;
    bundle.jn_low = evaluation.j_bkg_face_low_mid;
    bundle.jn_high = evaluation.j_bkg_face_high_mid;
    bundle.jn_final = evaluation.j_bkg_face_mid;
    bundle.je_low = evaluation.j_bkg_energy_low_cell_mid;
    bundle.je_center = evaluation.j_bkg_energy_center_cell_mid;
    bundle.je_high = evaluation.j_bkg_energy_high_cell_mid;
    bundle.je_final = evaluation.j_bkg_energy_cell_mid;
    bundle.gstar_je_low = evaluation.j_bkg_energy_low_debug_face;
    bundle.gstar_je_center = evaluation.j_bkg_energy_center_debug_face;
    bundle.gstar_je_high = evaluation.j_bkg_energy_high_debug_face;
    bundle.gstar_je_final = evaluation.j_bkg_energy_debug_face;
    bundle.final_state_mass.resize(static_cast<size_t>(sg.nx_local) *
                                   Param::Nvmu);
    if (evaluation.species_np1.f.size() >=
        static_cast<size_t>(sg.nx_total) * Param::Nvmu) {
        for (int ix = 0; ix < sg.nx_local; ++ix) {
            const size_t source = static_cast<size_t>(sg.nghost + ix) *
                                  Param::Nvmu;
            const size_t target = static_cast<size_t>(ix) * Param::Nvmu;
            std::copy(evaluation.species_np1.f.begin() + source,
                      evaluation.species_np1.f.begin() + source + Param::Nvmu,
                      bundle.final_state_mass.begin() + target);
        }
    }
    bundle.x_low_state_hash = evaluation.x_low_state_hash;
    bundle.u_low_state_hash = evaluation.u_low_state_hash;
    bundle.x_high_state_hash = evaluation.x_high_state_hash;
    bundle.u_high_state_hash = evaluation.u_high_state_hash;
    bundle.x_low_field_hash = evaluation.x_low_field_hash;
    bundle.u_low_field_hash = evaluation.u_low_field_hash;
    bundle.x_high_field_hash = evaluation.x_high_field_hash;
    bundle.u_high_field_hash = evaluation.u_high_field_hash;
    bundle.start_field_hash = evaluation.start_field_hash;
    bundle.end_field_hash = evaluation.end_field_hash;
    bundle.x_low_time_layer = evaluation.x_low_time_layer;
    bundle.u_low_time_layer = evaluation.u_low_time_layer;
    bundle.x_high_time_layer = evaluation.x_high_time_layer;
    bundle.u_high_time_layer = evaluation.u_high_time_layer;
    bundle.x_low_state_hash_history = evaluation.x_low_state_hash_history;
    bundle.u_low_state_hash_history = evaluation.u_low_state_hash_history;
    bundle.x_high_state_hash_history = evaluation.x_high_state_hash_history;
    bundle.u_high_state_hash_history = evaluation.u_high_state_hash_history;
    bundle.u_field_hash_history = evaluation.u_field_hash_history;
    bundle.coupling_substep_seam_audit =
        evaluation.coupling_substep_seam_audit;
    bundle.state_advanced = evaluation.state_advanced ? 1 : 0;
    bundle.operator_failed = evaluation.failed ? 1 : 0;
    bundle.failure_reason = evaluation.failure_reason;
    bundle.failure_iteration = evaluation.failure_iteration;
    bundle.failure_substep = evaluation.failure_substep;
    bundle.substeps_used = evaluation.substeps_used;
    bundle.low_order_candidate_min = evaluation.low_order_candidate_min;
    bundle.low_order_negative_count = evaluation.low_order_negative_count;
    bundle.low_order_negative_mass = evaluation.low_order_negative_mass;
    bundle.final_candidate_min = evaluation.failure_final_min;
    double local_final_min = std::numeric_limits<double>::infinity();
    if (evaluation.species_np1.f.size() >=
        static_cast<size_t>(sg.nx_total) * Param::Nvmu) {
        for (int ix = 0; ix < sg.nx_local; ++ix)
            for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k)
                    local_final_min = std::min(local_final_min,
                        evaluation.species_np1.f[idx3(sg.nghost + ix, j, k)]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &local_final_min, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    if (std::isfinite(local_final_min))
        bundle.final_candidate_min = local_final_min;

    const auto finite_vector = [](const std::vector<double>& values) {
        for (size_t p = 0; p < values.size(); ++p)
            if (!std::isfinite(values[p])) return false;
        return true;
    };
    bool outputs_finite = finite_vector(evaluation.species_np1.f) &&
        finite_vector(evaluation.fields_np1.Ex_face) &&
        finite_vector(evaluation.fields_np1.Ex) &&
        finite_vector(bundle.fx_high) && finite_vector(bundle.fu_center) &&
        finite_vector(bundle.jn_high) && finite_vector(bundle.je_center) &&
        finite_vector(bundle.gstar_je_center);
    int global_outputs_finite = outputs_finite ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_outputs_finite, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    bundle.outputs_finite = global_outputs_finite;
    // Keep the legacy field for existing tests, but it now means finite data,
    // not merely "the operator did not report failure".
    bundle.finite = bundle.outputs_finite;
    bundle.fct_active = evaluation.limiter_active_fraction > 0.0 ? 1 : 0;
    bundle.limiter_active_fraction = evaluation.limiter_active_fraction;
    bundle.limiter_min_alpha = evaluation.limiter_min_alpha;
    bundle.x_limiter_active_fraction =
        evaluation.x_limiter_active_fraction;
    bundle.x_limiter_min_alpha = evaluation.x_limiter_min_alpha;
    bundle.u_limiter_active_fraction =
        evaluation.u_limiter_active_fraction;
    bundle.u_limiter_min_alpha = evaluation.u_limiter_min_alpha;
    bundle.donor_beta_applied_count =
        evaluation.fct_donor_beta_applied_count;
    bundle.donor_beta_min = evaluation.fct_donor_beta_min;
    bundle.r_couple = evaluation.stage5_r_couple;
    bundle.r_fv = evaluation.stage5_r_fv;
    bundle.delta_ke_bkg = evaluation.delta_ke_bkg;
    bundle.stage5_u_energy_moment = evaluation.stage5_u_energy_moment;
    bundle.stage5_spatial_energy_boundary =
        evaluation.stage5_spatial_energy_boundary;
    bundle.mass_change = evaluation.stage2_mass_change;
    bundle.mass_scale = evaluation.stage2_mass_scale;
    bundle.mass_residual = evaluation.stage2_mass_residual;
    bundle.momentum_change = evaluation.stage2_momentum_change;
    bundle.momentum_scale = evaluation.stage2_momentum_scale;
    bundle.momentum_residual = evaluation.stage2_momentum_residual;
    bundle.u_boundary_particle = evaluation.u_boundary_particle_outflow;
    bundle.u_boundary_momentum = evaluation.u_boundary_momentum_outflow;
    bundle.u_boundary_energy = evaluation.u_boundary_energy_outflow;
    bundle.u_boundary_energy_lower = evaluation.u_boundary_energy_lower;
    bundle.u_boundary_energy_upper = evaluation.u_boundary_energy_upper;
    bundle.final_flux_current_moment_audit_valid =
        evaluation.final_flux_current_moment_audit_valid;
    bundle.final_flux_current_moment_audit_finite =
        evaluation.final_flux_current_moment_audit_finite;
    bundle.final_flux_to_jn_linf = evaluation.final_flux_to_jn_linf;
    bundle.final_flux_to_jn_scale = evaluation.final_flux_to_jn_scale;
    bundle.final_flux_to_je_linf = evaluation.final_flux_to_je_linf;
    bundle.final_flux_to_je_scale = evaluation.final_flux_to_je_scale;
    bundle.final_flux_to_gstar_je_linf =
        evaluation.final_flux_to_gstar_je_linf;
    bundle.final_flux_to_gstar_je_scale =
        evaluation.final_flux_to_gstar_je_scale;
    return bundle;
}
