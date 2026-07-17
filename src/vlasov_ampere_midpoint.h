#ifndef VLASOV_AMPERE_MIDPOINT_H
#define VLASOV_AMPERE_MIDPOINT_H

#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"

#include <array>
#include <string>
#include <vector>

class VlasovAmpereMidpointSolver {
public:
    // Dynamic part of the accepted regional closure partition.  A fixed
    // operator audit reuses this layout rather than inferring it from a
    // frozen Beam state/current.
    struct CouplingRegionLayout {
        int beam_front_ix;
        double wave_core_end_m;
    };

    // Serialized only for diagnostic-level=2 fixed-state audits.  It owns
    // complete endpoint objects; no solver scratch arrays are persisted.
    struct MidpointAuditState {
        long long step;
        double time_s;
        double dt_s;
        int substeps_used;
        int nonlinear_iterations;
        bool low_order_only;
        bool high_order_enabled;
        bool fct_enabled;
        std::string acceptance_type;
        Species bkg_n;
        Species guess_np1;
        Species operator_input_guess;
        EMFields fields_n;
        EMFields fields_end_guess;
        std::vector<double> j_beam_face_mid;
        std::vector<double> reference_jn_face;
        std::vector<double> reference_je_cell;
        std::vector<double> reference_gstar_je_face;
        // [JN(0), JN(L), G*JE(0), G*JE(L), pair(0), pair(L)].
        std::array<double, 6> periodic_seam_face_audit;
        CouplingRegionLayout coupling_layout;
        // The two global closure residuals belong to the same fixed-operator
        // reference as JN/JE/G*JE and must be reproduced independently.
        double reference_stage5_r_fv;
        double reference_stage5_r_couple;
    };
    enum NegativeDebtLevel {
        NEG_DEBT_OK = 0,
        NEG_DEBT_WARN = 1,
        NEG_DEBT_REPAIR = 2,
        NEG_DEBT_ABORT = 3
    };

    struct CurrentDiagnostics {
        double residual_if_charge;
        double residual_if_ampere;
        double e_dot_j_charge;
        double e_dot_j_energy;
        double e_dot_j_ampere;
        double max_abs_j_charge;
        double max_abs_j_energy;
        double max_abs_j_ampere;
        double max_abs_j_charge_minus_ampere;
        double max_abs_j_energy_minus_ampere;
    };

    // 7.1.6: per-direction flux-positivity diagnostics
    struct FluxPositivityDiag {
        double min_f_before;
        double min_f_low;
        double min_f_final;
        double low_order_failed_count;
        double alpha_active_fraction;
        double alpha_min;
        double alpha_core_fraction;
        double alpha_boundary_fraction;
        double negative_mass_prevented;
    };

    // 7.1.6: per-direction flux-defect diagnostics
    struct FluxDefectDiag {
        double mass_defect;
        double momentum_defect;
        double energy_defect;
        double boundary_mass_loss;
        double boundary_energy_loss;
    };

    struct CouplingRegionDiagnostics {
        double sum_rj;
        double sum_abs_rj;
        double sum_sq_rj;
        double linf_rj;
        double sum_rk;
        double sum_abs_rk;
        double sum_sq_rk;
        double linf_rk;
        double max_abs_jn_minus_gstar_je;
        int max_abs_jn_minus_gstar_je_face;
        double face_work_jn;
        double delta_ke_bkg;
        double fct_work;
        // Additive partition of sum_rj.  These are diagnostic layers only:
        // J_N remains the charge-conserving current used by Ampere.
        double r_couple_centered;
        double r_couple_upwind_stabilization;
        double r_couple_fct_stabilization;
        long long face_count;
        long long cell_count;
    };

    struct LowOrderFailureAudit {
        int valid;
        int rank;
        int ix;
        int iv;
        int imu;
        int region; // 0=core, 1=boundary
        double severity;
        double f_input;
        double f_after_x;
        double dx_div;
        double dmu_div_used;
        double du_div_low;
        double f_low;
        double left_lower_flux;
        double left_upper_flux;
        double right_lower_flux;
        double right_upper_flux;
        double left_lower_scale;
        double left_upper_scale;
        double right_lower_scale;
        double right_upper_scale;
        double left_lower_donor_f;
        double left_upper_donor_f;
        double right_lower_donor_f;
        double right_upper_donor_f;
        double lower_characteristic;
        double upper_characteristic;
        double moment_weight;
        double cell_budget;
        double low_order_failed_count;
        int left_lower_donor_index;
        int left_upper_donor_index;
        int right_lower_donor_index;
        int right_upper_donor_index;
    };

    struct AcceptedTransportSnapshot {
        int valid;
        int substep;
        double min_mass;
        int ix;
        int iv;
        int imu;
        double m_low;
        double m_candidate;
        double m_safe;
        double outflow;
        double inflow;
        double beta;
        double alpha[4];
        double a_limited[4];
        double a_safe[4];
        int mpi_interface;
        int k_zero;
        int next_step_donor;
    };

    struct Result {
        Species species_np1;
        BeamPIC beam_np1;
        EMFields fields_np1;
        std::vector<double> j_bkg_face_mid;
        // Fixed-state face-pairing layers.  The x high-order transport has no
        // separate centered variant, so j_bkg_face_center_mid aliases its
        // high-order layer; keeping it explicit makes the audit contract
        // complete without inventing a second transport formula.
        std::vector<double> j_bkg_face_low_mid;
        std::vector<double> j_bkg_face_center_mid;
        std::vector<double> j_bkg_face_high_mid;
        std::vector<double> j_beam_face_mid;
        std::vector<double> j_total_face_mid;
        std::vector<double> j_bkg_energy_debug_face;
        std::vector<double> j_bkg_energy_low_debug_face;
        std::vector<double> j_bkg_energy_center_debug_face;
        std::vector<double> j_bkg_energy_high_debug_face;
        // Cell-centered J_E assembled directly from final u-face energy
        // fluxes. It is never formed by dividing power by E.
        std::vector<double> j_bkg_energy_cell_mid;
        // Present for fixed midpoint-operator audits only.  These are the
        // exact final FV fluxes used to form the returned candidate state.
        std::vector<double> final_x_flux;
        std::vector<double> final_u_flux;
        Species operator_input_guess;
        // Exact endpoint field guess supplied to the final production-kernel
        // evaluation.  This is distinct from fields_np1, which is the result
        // of Ampere's update for that evaluation.
        EMFields operator_input_fields_end_guess;
        AcceptedTransportSnapshot accepted_transport;
        CurrentDiagnostics current_diag;
        double delta_ke_bkg;
        double delta_ke_beam;
        double field_work_bkg;
        double field_work_beam;
        double energy_residual_bkg;
        double energy_pair_residual_bkg;
        double beam_lag_energy_residual;
        double limiter_energy_defect_positive;
        double limiter_energy_defect_negative;
        double limiter_energy_defect_core;
        double limiter_energy_defect_boundary;
        double u_boundary_particle_outflow;
        double u_boundary_momentum_outflow;
        double u_boundary_energy_outflow;
        // Stage 5: finite-volume and field-current closure are deliberately
        // separate diagnostics.  No compatibility projection is active yet.
        double stage5_r_fv;
        double stage5_r_couple;
        std::array<CouplingRegionDiagnostics, 6> coupling_regions;
        // All three are global values formed by exactly one reduction of the
        // corresponding local regional terms.  Their reconstruction errors
        // audit the accepted finite-volume/Ampere state.
        double coupling_rj_global_sum;
        double coupling_rk_global_sum;
        double coupling_face_work_jn_global_sum;
        double coupling_rj_reconstruction_error;
        double coupling_rk_reconstruction_error;
        double coupling_face_work_jn_reconstruction_error;
        double coupling_rj_centered_global_sum;
        double coupling_rj_upwind_global_sum;
        double coupling_rj_fct_global_sum;
        double coupling_rj_centered_reconstruction_error;
        double coupling_rj_upwind_reconstruction_error;
        double coupling_rj_fct_reconstruction_error;
        int coupling_beam_front_ix;
        double coupling_wave_core_end_m;
        // [JN(0), JN(L), G*JE(0), G*JE(L), pair(0), pair(L)].
        std::array<double, 6> periodic_seam_face_audit;
        // R_couple decomposition.  The centered term is the physical
        // high-order u-force candidate; the other two terms expose the
        // boundary upwind and FCT safety contributions without changing J_N.
        double stage5_r_couple_centered;
        double stage5_r_couple_upwind_stabilization;
        double stage5_r_couple_fct_stabilization;
        double stage5_r_total;
        double stage5_r_physical_balance;
        double stage5_u_energy_moment;
        double stage5_u_momentum_moment;
        double stage5_spatial_energy_boundary;
        double stage5_spatial_momentum_boundary;
        double stage5_jn_minus_je_l2;
        double stage5_jn_minus_je_linf;
        double stage5_energy_speed_candidate_linf;
        double stage5_energy_speed_candidate_rel;
        double stage5_fct_budget_violation;
        int stage5_compatibility_enabled;
        int stage5_projection_skipped_zero_field;
        int stage5_projection_skipped_raw_residual;
        double stage5_correction_energy_target;
        double stage5_correction_energy_achieved;
        double stage5_correction_added_momentum_linf;
        double stage5_correction_l1;
        double stage5_correction_linf;
        double stage5_correction_relative;
        long long stage5_correction_active_bounds;
        int stage5_correction_infeasible;
        int stage5_correction_retry_count;
        // P1 FCT audit invariants.  These only audit the existing algorithm;
        // P2--P4 will change transfer units and the MPI owner protocol.
        double fct_high_low_identity_linf;
        double fct_high_low_identity_violation;
        double fct_high_low_identity_worst_residual;
        double fct_high_low_identity_worst_scale;
        double fct_high_low_identity_worst_relative;
        double fct_high_low_identity_ratio_linf;
        int fct_high_low_identity_worst_ix;
        int fct_high_low_identity_worst_iv;
        int fct_high_low_identity_worst_imu;
        double fct_donor_capacity_violation;
        double fct_donor_worst_m_low;
        double fct_donor_worst_ax_left;
        double fct_donor_worst_ax_right;
        double fct_donor_worst_au_lower;
        double fct_donor_worst_au_upper;
        double fct_donor_worst_alpha_x_left;
        double fct_donor_worst_alpha_x_right;
        double fct_donor_worst_alpha_u_lower;
        double fct_donor_worst_alpha_u_upper;
        double fct_donor_worst_outflow;
        double fct_donor_worst_scale;
        double fct_donor_worst_relative;
        int fct_donor_worst_ix;
        int fct_donor_worst_iv;
        int fct_donor_worst_imu;
        int fct_donor_roundoff_warning;
        long long fct_donor_beta_applied_count;
        double fct_donor_beta_min;
        double fct_interface_checksum_linf;
        double fct_interface_checksum_violation;
        double fct_low_order_tolerance_linf;
        double low_order_candidate_min;
        long long low_order_negative_count;
        double low_order_negative_mass;
        long long low_order_roundoff_zeroed_count;
        double low_order_roundoff_zeroed_mass;
        // Raw high-order candidate before FCT alpha/beta correction.
        double fct_high_candidate_min;
        double fct_high_candidate_donor_excess;
        long long fct_controlled_injection_count;
        double fct_final_scratch_min;
        double fct_final_tolerance_linf;
        // Final-scratch negative values set to +0.0 solely because they are
        // within the local floating-point summation bound.
        long long fct_roundoff_zeroed_count;
        double fct_roundoff_zeroed_mass;
        // 0=none, 1=CFL limit, 2=low-order positivity,
        // 3=FCT final positivity, 4=non-finite state,
        // 5=high-low identity, 6=donor capacity, 7=interface checksum,
        // 12=nonuniform high-order transport disabled by configuration.
        int failure_reason;
        int failure_iteration;
        int failure_substep;
        double failure_global_cfl;
        double failure_low_min;
        double failure_final_min;
        int failure_worst_ix;
        int failure_worst_iv;
        int failure_worst_imu;
        int failure_worst_is_core;
        int failure_worst_is_tail;
        double failure_low_m_in;
        double failure_low_m_low;
        double failure_low_transfer_x_left;
        double failure_low_transfer_x_right;
        double failure_low_transfer_u_lower;
        double failure_low_transfer_u_upper;
        double failure_low_out_x_left;
        double failure_low_out_x_right;
        double failure_low_out_u_lower;
        double failure_low_out_u_upper;
        double failure_low_in_x_left;
        double failure_low_in_x_right;
        double failure_low_in_u_lower;
        double failure_low_in_u_upper;
        double failure_low_cfl_x;
        double failure_low_cfl_u;
        double failure_low_work_input_min;
        int failure_low_on_mpi_interface;
        double continuity_residual_bkg;
        double beam_continuity_residual;
        double nonlinear_residual;
        double residual_E;
        double residual_f;
        double residual_J_bkg;
        double residual_J_beam;
        double limiter_active_fraction;
        double limiter_min_alpha;
        double limiter_active_fraction_core;
        double limiter_active_fraction_boundary;
        double limiter_min_alpha_core;
        double limiter_min_alpha_boundary;
        double limiter_energy_defect;
        double limiter_mass_defect;
        double limiter_momentum_defect;
        double x_limiter_energy_defect;
        double x_limiter_mass_defect;
        double mu_low_u_alpha_min;
        double mu_low_u_limiter_active_fraction;
        double mu_low_u_energy_delta;
        double mu_low_u_alpha_min_boundary;
        double mu_low_u_alpha_min_core;
        double mu_low_u_limiter_active_fraction_boundary;
        double mu_low_u_limiter_active_fraction_core;
        double mu_low_u_energy_delta_boundary;
        double mu_low_u_energy_delta_core;
        double mu_low_u_u_eff0;
        double mu_low_u_moment_weight0;
        double mu_low_u_mu_flux_scale0;
        double mu_low_u_half_dt_inv_shell0;
        double mu_low_u_dimless_scale0;
        double mu_low_u_endpoint_flux_max;
        double remap_active_fraction;
        long long remap_cell_count;
        double low_u_subcycle_active_fraction;
        double low_u_average_subcycles;
        int low_u_max_subcycles;
        // 7.1.1: unified flux diagnostics for Result
        double x_low_order_failed_count;
        double x_low_input_min_f;
        double x_low_max_cfl;
        double x_low_output_min_f;
        double x_low_failed_count;
        double x_low_input_neg_mass;
        double x_low_input_rel_neg;
        double x_low_output_rel_neg;
        double x_low_input_core_failed_count;
        double x_low_input_debt_accepted;
        int x_low_failure_kind;
        double x_final_min_f;
        double x_final_failed_count;
        double x_final_failed_max_debt;
        int x_final_worst_ix;
        int x_final_worst_iv;
        int x_final_worst_imu;
        int x_final_failure_region;
        double x_final_core_failed_count;
        double x_final_boundary_failed_count;
        // 7.1.6: per-direction flux diagnostics (0=x, 1=u, 2=mu)
        FluxPositivityDiag flux_pos[3];
        FluxDefectDiag      flux_defect[3];
        double region_u_limiter_energy_boundary[2];
        double region_u_limiter_energy_core[2];
        double region_abs_u_limiter_energy_boundary[2];
        double region_abs_u_limiter_energy_core[2];
        double region_limiter_active_fraction_boundary[2];
        double region_limiter_active_fraction_core[2];
        std::vector<double> stage_min_f;
        std::vector<double> stage_neg_mass;
        std::vector<long long> stage_neg_cell_count;
        std::vector<double> stage_low_u_neg_mass;
        std::vector<double> low_u_neg_added_by_div;
        std::vector<double> stage_core_low_u_min_f;
        std::vector<double> stage_min_f_core_by_u;
        std::vector<double> stage_neg_mass_core_by_u;
        std::vector<long long> stage_neg_cell_count_core_by_u;
        std::vector<double> stage_min_f_boundary_by_u;
        std::vector<double> stage_neg_mass_boundary_by_u;
        std::vector<long long> stage_neg_cell_count_boundary_by_u;
        double x_negative_mass_before_repair;
        double x_mass_added_by_positivity_repair;
        double floor_repair_mass;
        double floor_repair_energy;
        double floor_repair_core_fraction;
        int negative_debt_level;
        double neg_mass_boundary;
        double neg_mass_core;
        double neg_mass_tail;
        double neg_mass_total_guard;
        double neg_mass_core_fraction;
        double neg_mass_boundary_fraction;
        double neg_mass_tail_fraction;
        double neg_energy_core_abs;
        double neg_energy_core_fraction;
        double neg_energy_boundary_abs;
        double neg_energy_boundary_fraction;
        double neg_current_core_abs;
        double neg_current_core_fraction;
        double neg_debt_min_f_boundary;
        double neg_debt_min_f_core;
        double neg_debt_min_f_tail;
        long long neg_cell_boundary;
        long long neg_cell_core;
        long long neg_cell_tail;
        int debt_action;
        int limiter_reason;
        double alpha_core_min;
        double alpha_boundary_min;
        double alpha_tail_min;
        double limiter_modified_J_bkg_norm;
        double limiter_modified_J_bkg_boundary_norm;
        double limiter_modified_energy_norm;
        double low_u_mu_neg_mass_fraction;
        double low_u_mu_neg_energy_fraction;
        double low_u_mu_neg_current_fraction;
        double boundary_force_Cu_max;
        double boundary_force_Cmu_max;
        int boundary_force_nsub_max;
        long long boundary_force_remap_cell_count;
        double boundary_mu_low_L1_before;
        double boundary_mu_low_L1_after;
        double boundary_mu_high_L1_after;
        double J_bkg_neg_boundary;
        double delta_E_neg_boundary;
        double boundary_force_remap_mass_loss;
        double boundary_force_remap_energy_loss;
        double alpha_interface_BQ_min;
        double alpha_interface_QC_min;
        double interface_BQ_flux;
        double interface_BQ_high_correction;
        double interface_QC_flux_into_core;
        double interface_QC_high_correction_into_core;
        double boundary_energy_diagnostic_invalid;
        int trial_failure_downgraded;
        int accepted_with_negative_debt;
        int state_advanced;
        double positivity_energy_defect;
        double positivity_mass_defect;
        double u_force_alpha_min;
        double u_force_alpha_active_frac;
        int u_flux_audit_valid;
        int u_flux_audit_rank;
        int u_flux_audit_ix;
        int u_flux_audit_iv;
        int u_flux_audit_imu;
        double u_flux_audit_severity;
        double u_flux_audit_f0;
        double u_flux_audit_f_after_x;
        double u_flux_audit_f_low;
        double u_flux_audit_f_high;
        double u_flux_audit_alpha;
        double u_flux_audit_du_div_low;
        double u_flux_audit_du_div_high;
        double u_flux_audit_du_div_final;
        double u_flux_audit_updated;
        double u_flux_audit_final_xl_lo;
        double u_flux_audit_final_xl_hi;
        double u_flux_audit_final_xr_lo;
        double u_flux_audit_final_xr_hi;
        double u_flux_audit_lower_flux;
        double u_flux_audit_upper_flux;
        double u_flux_audit_lower_donor_f;
        double u_flux_audit_upper_donor_f;
        double u_flux_audit_cfl;
        double u_flux_audit_line_positive_peak;
        double u_flux_audit_line_negative_mass;
        double u_flux_audit_allowed_debt;
        int u_flux_audit_lower_donor_iv;
        int u_flux_audit_upper_donor_iv;
        int u_flux_audit_failure_kind;
        LowOrderFailureAudit u_low_failure_audit;
        LowOrderFailureAudit mu_low_failure_audit;
        int finite_failure_mask;
        long long u_low_order_failed_count;
        long long mu_low_order_failed_count;
        long long u_final_negative_hard_count;
        long long mu_final_negative_hard_count;
        double mu_final_failed_max_debt;
        int mu_final_worst_ix;
        int mu_final_worst_iv;
        int mu_final_worst_imu;
        int mu_final_failure_region;
        double mu_final_core_failed_count;
        double mu_final_boundary_failed_count;
        double mu_final_audit_f_base;
        double mu_final_audit_f_low;
        double mu_final_audit_f_high;
        double mu_final_audit_f_final;
        double mu_final_audit_f_floor;
        double mu_final_audit_A_left;
        double mu_final_audit_A_right;
        double mu_final_audit_P_minus;
        double mu_final_audit_Q_minus;
        double mu_final_audit_R_minus;
        double mu_final_audit_alpha_left_face;
        double mu_final_audit_alpha_right_face;
        double mu_final_audit_dt_inv_shell;
        double mu_final_audit_mu_dot_lower;
        double mu_final_audit_mu_dot_upper;
        int finite_stage_failure_valid;
        int finite_stage_failure_kind;
        int finite_stage_failure_rank;
        int finite_stage_failure_ix;
        int finite_stage_failure_iv;
        int finite_stage_failure_imu;
        double finite_stage_failure_severity;
        double finite_stage_failure_f_base;
        double finite_stage_failure_f_low;
        double finite_stage_failure_f_high;
        double finite_stage_failure_f_final;
        double finite_stage_failure_du_div;
        double finite_stage_failure_dmu_div;
        double finite_stage_failure_dx_div;
        double f_neg_min;
        double f_neg_ratio_max;
        double f_neg_mass_total;
        long long f_neg_cell_count;
        int f_neg_ix;
        int f_neg_iv;
        int f_neg_imu;
        int nonlinear_iterations;
        std::vector<double> midpoint_residual_e_history;
        std::vector<double> midpoint_residual_j_bkg_history;
        std::vector<double> midpoint_residual_j_beam_history;
        std::vector<double> midpoint_residual_f_history;
        bool converged;
        bool failed;
        bool soft_accepted;
        bool soft_unconverged;
        bool protected_converged;
        int transport_low_order_only;
        int substeps_used;
    };

    typedef Result MidpointOperatorEvaluation;

    VlasovAmpereMidpointSolver();

    void set_step_diagnostics_enabled(bool enabled) {
        step_diagnostics_enabled_ = enabled;
    }

    // Test harnesses can disable the open PIC source while retaining the
    // production background/Vlasov/Ampere implementation.  Normal runs keep
    // this enabled by default.
    void set_beam_enabled(bool enabled) { beam_enabled_ = enabled; }
    void set_low_order_only(bool enabled) { low_order_only_ = enabled; }
    // Explicit production configuration for MUSCL/CTU transport on a
    // nonuniform cylindrical velocity grid.
    void set_nonuniform_high_order_enabled(bool enabled) {
        nonuniform_high_order_enabled_ = enabled;
    }
    void set_fct_enabled(bool enabled) { fct_enabled_ = enabled; }
    void set_max_midpoint_iterations(int iterations) {
        max_midpoint_iterations_ = iterations > 0 ? iterations : 1;
    }
    void set_capture_midpoint_input(bool enabled) {
        capture_midpoint_input_ = enabled;
    }
    bool low_order_only() const { return low_order_only_; }
    bool nonuniform_high_order_enabled() const {
        return nonuniform_high_order_enabled_;
    }
    bool fct_enabled() const { return fct_enabled_; }
    int max_midpoint_iterations() const { return max_midpoint_iterations_; }
    // Test-only A/B control.  Production always uses the centered high-order
    // candidate at every x; enabling this restores the former boundary
    // upwind branch solely for comparison.
    void set_legacy_boundary_upwind_high_candidate_for_test(bool enabled) {
        legacy_boundary_upwind_high_candidate_for_test_ = enabled;
    }
    // Test-only convergence tracing does not alter the midpoint equation,
    // relaxation, fluxes, or production acceptance rules.
    void set_midpoint_iteration_trace_for_test(bool enabled) {
        midpoint_iteration_trace_for_test_ = enabled;
    }
    void set_midpoint_iteration_trace(bool enabled) {
        midpoint_iteration_trace_for_test_ = enabled;
    }
    // Raw high-candidate accounting is needed by the controlled FCT test,
    // but is intentionally off in production runs to avoid extra per-cell
    // long-double audit work.
    void set_fct_activation_audit_enabled(bool enabled) {
        fct_activation_audit_enabled_ = enabled;
    }
    // Test-only: inject one conservative anti-diffusive x-face correction
    // after the low-order state passes its transport gate.  Never enabled by
    // the production executable.
    void set_controlled_fct_flux_injection_enabled(bool enabled) {
        controlled_fct_flux_injection_enabled_ = enabled;
    }
    // Section-7.3 no-FCT verification may retain finite, non-macroscopic
    // negative debt so the raw high-order dynamics can be audited.  This is
    // disabled for the production solver and never relaxes NaN/Inf checks.
    void set_allow_finite_negative_debt_for_test(bool enabled) {
        allow_finite_negative_debt_for_test_ = enabled;
    }
    // Restart uses the exact production periodic background exchange.
    void synchronize_background_ghosts(Species& sp, const SpatialGrid& sg,
                                       int mpi_rank, int mpi_size) const {
        exchange_ghosts_x_persistent(sp, sg, mpi_rank, mpi_size);
    }

    Result advance_background_and_fields(const Species& bkg_n,
                                         const BeamPIC& beam_n,
                                         const EMFields& fields_n,
                                         const SpatialGrid& sg,
                                         double dt,
                                         double time,
                                         int mpi_rank,
                                         int mpi_size);
    // One production-kernel evaluation at an externally supplied Picard
    // endpoint.  It is intentionally side-effect free: callers receive a
    // Result and decide whether to advance physical state.
    MidpointOperatorEvaluation evaluate_fixed_midpoint_operator(
        const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
        const Species& guess_np1, const EMFields& fields_end_guess,
        const std::vector<double>& fixed_j_beam_face_mid,
        const CouplingRegionLayout& coupling_layout,
        const SpatialGrid& sg, double dt, double time, int mpi_rank,
        int mpi_size) const;

private:
    /*
     * Legacy FluxPack removed.  The cylindrical solver keeps its phase-space
     * fluxes as iteration-local storage so old (u,mu) arrays cannot be used
     * by a future caller accidentally.
     */
#if 0
    struct FluxPack {
        std::vector<double> x_high;
        std::vector<double> x_low;
        std::vector<double> x_final;
        std::vector<double> cell_alpha_u;   // 7.1.1: u-direction per-cell alpha
        std::vector<double> cell_alpha_x;   // 7.1.1: x-direction per-cell alpha
        std::vector<double> j_bkg_face;
        // 7.1.1: unified three-flux storage
        std::vector<double> u_high;
        std::vector<double> u_low;
        std::vector<double> u_final;
        std::vector<double> u_high_cell;
        std::vector<double> u_low_cell;
        std::vector<double> u_final_cell;
        std::vector<double> mu_high;
        std::vector<double> mu_low;
        std::vector<double> mu_final;
        std::vector<double> mu_high_cell;
        std::vector<double> mu_low_cell;
        std::vector<double> mu_final_cell;
        std::vector<double> cell_alpha_mu;
        double x_low_order_failed_count;
        double x_low_input_min_f;
        double x_low_max_cfl;
        double x_low_output_min_f;
        double x_low_failed_count;
        double x_low_input_neg_mass;
        double x_low_input_rel_neg;
        double x_low_output_rel_neg;
        double x_low_input_core_failed_count;
        double x_low_input_debt_accepted;
        int x_low_failure_kind;
        double x_final_min_f;
        double x_final_failed_count;
        double x_final_failed_max_debt;
        int x_final_worst_ix;
        int x_final_worst_iv;
        int x_final_worst_imu;
        int x_final_failure_region;
        double x_final_core_failed_count;
        double x_final_boundary_failed_count;
        // 7.1.6: per-direction flux diagnostics (0=x, 1=u, 2=mu)
        FluxPositivityDiag flux_pos[3];
        FluxDefectDiag      flux_defect[3];
        double limiter_active_fraction;
        double limiter_min_alpha;
        double limiter_active_fraction_core;
        double limiter_active_fraction_boundary;
        double limiter_min_alpha_core;
        double limiter_min_alpha_boundary;
        double limiter_energy_defect;
        double limiter_mass_defect;
        double limiter_momentum_defect;
        double x_limiter_energy_defect;
        double x_limiter_mass_defect;
        double mu_low_u_alpha_min;
        double mu_low_u_limiter_active_fraction;
        double mu_low_u_energy_delta;
        double mu_low_u_alpha_min_boundary;
        double mu_low_u_alpha_min_core;
        double mu_low_u_limiter_active_fraction_boundary;
        double mu_low_u_limiter_active_fraction_core;
        double mu_low_u_energy_delta_boundary;
        double mu_low_u_energy_delta_core;
        double mu_low_u_u_eff0;
        double mu_low_u_moment_weight0;
        double mu_low_u_mu_flux_scale0;
        double mu_low_u_half_dt_inv_shell0;
        double mu_low_u_dimless_scale0;
        double mu_low_u_endpoint_flux_max;
        double remap_active_fraction;
        long long remap_cell_count;
        double low_u_subcycle_active_fraction;
        double low_u_average_subcycles;
        int low_u_max_subcycles;
        double region_u_limiter_energy_boundary[2];
        double region_u_limiter_energy_core[2];
        double region_abs_u_limiter_energy_boundary[2];
        double region_abs_u_limiter_energy_core[2];
        double region_limiter_active_fraction_boundary[2];
        double region_limiter_active_fraction_core[2];
        std::vector<double> stage_min_f;
        std::vector<double> stage_neg_mass;
        std::vector<long long> stage_neg_cell_count;
        std::vector<double> stage_low_u_neg_mass;
        std::vector<double> low_u_neg_added_by_div;
        std::vector<double> stage_core_low_u_min_f;
        std::vector<double> stage_min_f_core_by_u;
        std::vector<double> stage_neg_mass_core_by_u;
        std::vector<long long> stage_neg_cell_count_core_by_u;
        std::vector<double> stage_min_f_boundary_by_u;
        std::vector<double> stage_neg_mass_boundary_by_u;
        std::vector<long long> stage_neg_cell_count_boundary_by_u;
        double negative_mass_before_repair;
        double mass_added_by_positivity_repair;
        double floor_repair_mass;
        double floor_repair_energy;
        double floor_repair_core_fraction;
        int negative_debt_level;
        double neg_mass_boundary;
        double neg_mass_core;
        double neg_mass_tail;
        double neg_mass_total_guard;
        double neg_mass_core_fraction;
        double neg_mass_boundary_fraction;
        double neg_mass_tail_fraction;
        double neg_energy_core_abs;
        double neg_energy_core_fraction;
        double neg_energy_boundary_abs;
        double neg_energy_boundary_fraction;
        double neg_current_core_abs;
        double neg_current_core_fraction;
        double neg_debt_min_f_boundary;
        double neg_debt_min_f_core;
        double neg_debt_min_f_tail;
        long long neg_cell_boundary;
        long long neg_cell_core;
        long long neg_cell_tail;
        int debt_action;
        int limiter_reason;
        double alpha_core_min;
        double alpha_boundary_min;
        double alpha_tail_min;
        double limiter_modified_J_bkg_norm;
        double limiter_modified_J_bkg_boundary_norm;
        double limiter_modified_energy_norm;
        double low_u_mu_neg_mass_fraction;
        double low_u_mu_neg_energy_fraction;
        double low_u_mu_neg_current_fraction;
        double boundary_force_Cu_max;
        double boundary_force_Cmu_max;
        int boundary_force_nsub_max;
        long long boundary_force_remap_cell_count;
        double boundary_mu_low_L1_before;
        double boundary_mu_low_L1_after;
        double boundary_mu_high_L1_after;
        double J_bkg_neg_boundary;
        double delta_E_neg_boundary;
        double boundary_force_remap_mass_loss;
        double boundary_force_remap_energy_loss;
        double alpha_interface_BQ_min;
        double alpha_interface_QC_min;
        double interface_BQ_flux;
        double interface_BQ_high_correction;
        double interface_QC_flux_into_core;
        double interface_QC_high_correction_into_core;
        double boundary_energy_diagnostic_invalid;
        int trial_failure_downgraded;
        int accepted_with_negative_debt;
        double positivity_energy_defect;
        double positivity_mass_defect;
        double u_force_alpha_min;
        double u_force_alpha_active_frac;
        int u_flux_audit_valid;
        int u_flux_audit_rank;
        int u_flux_audit_ix;
        int u_flux_audit_iv;
        int u_flux_audit_imu;
        double u_flux_audit_severity;
        double u_flux_audit_f0;
        double u_flux_audit_f_after_x;
        double u_flux_audit_f_low;
        double u_flux_audit_f_high;
        double u_flux_audit_alpha;
        double u_flux_audit_du_div_low;
        double u_flux_audit_du_div_high;
        double u_flux_audit_du_div_final;
        double u_flux_audit_updated;
        double u_flux_audit_final_xl_lo;
        double u_flux_audit_final_xl_hi;
        double u_flux_audit_final_xr_lo;
        double u_flux_audit_final_xr_hi;
        double u_flux_audit_lower_flux;
        double u_flux_audit_upper_flux;
        double u_flux_audit_lower_donor_f;
        double u_flux_audit_upper_donor_f;
        double u_flux_audit_cfl;
        double u_flux_audit_line_positive_peak;
        double u_flux_audit_line_negative_mass;
        double u_flux_audit_allowed_debt;
        int u_flux_audit_lower_donor_iv;
        int u_flux_audit_upper_donor_iv;
        int u_flux_audit_failure_kind;
        LowOrderFailureAudit u_low_failure_audit;
        LowOrderFailureAudit mu_low_failure_audit;
        int finite_failure_mask;
        long long u_low_order_failed_count;
        long long mu_low_order_failed_count;
        long long u_final_negative_hard_count;
        long long mu_final_negative_hard_count;
        double mu_final_failed_max_debt;
        int mu_final_worst_ix;
        int mu_final_worst_iv;
        int mu_final_worst_imu;
        int mu_final_failure_region;
        double mu_final_core_failed_count;
        double mu_final_boundary_failed_count;
        double mu_final_audit_f_base;
        double mu_final_audit_f_low;
        double mu_final_audit_f_high;
        double mu_final_audit_f_final;
        double mu_final_audit_f_floor;
        double mu_final_audit_A_left;
        double mu_final_audit_A_right;
        double mu_final_audit_P_minus;
        double mu_final_audit_Q_minus;
        double mu_final_audit_R_minus;
        double mu_final_audit_alpha_left_face;
        double mu_final_audit_alpha_right_face;
        double mu_final_audit_dt_inv_shell;
        double mu_final_audit_mu_dot_lower;
        double mu_final_audit_mu_dot_upper;
        int finite_stage_failure_valid;
        int finite_stage_failure_kind;
        int finite_stage_failure_rank;
        int finite_stage_failure_ix;
        int finite_stage_failure_iv;
        int finite_stage_failure_imu;
        double finite_stage_failure_severity;
        double finite_stage_failure_f_base;
        double finite_stage_failure_f_low;
        double finite_stage_failure_f_high;
        double finite_stage_failure_f_final;
        double finite_stage_failure_du_div;
        double finite_stage_failure_dmu_div;
        double finite_stage_failure_dx_div;
        double f_neg_min;
        double f_neg_ratio_max;
        double f_neg_mass_total;
        long long f_neg_cell_count;
        int f_neg_ix;
        int f_neg_iv;
        int f_neg_imu;
        double finite_flux_max_negative;
        double finite_flux_relative_negative;
        double finite_flux_updated;
        double finite_flux_f0;
        double finite_flux_dx_div;
        double finite_flux_du_div;
        double finite_flux_dmu_div;
        int finite_flux_rank;
        int finite_flux_ix;
        int finite_flux_iv;
        int finite_flux_imu;
        int finite_flux_has_failure;
    };
#endif

    void reset_result(Result& result) const;
    void reset_current_diag(CurrentDiagnostics& diag) const;
    // Shared production operator kernel.  The production step and later
    // fixed-midpoint audit both enter through this boundary; keeping the
    // original loop order here avoids a second implementation of FCT/JN/JE.
    Result evaluate_production_midpoint_operator(const Species& bkg_n,
                                                 const BeamPIC& beam_n,
                                                 const EMFields& fields_n,
                                                 const SpatialGrid& sg,
                                                 double dt,
                                                 double time,
                                                 int mpi_rank,
                                                 int mpi_size,
                                                 int substeps_used,
                                                 const Species* fixed_guess_np1 = 0,
                                                 const EMFields* fixed_fields_end = 0,
                                                 const std::vector<double>* fixed_j_beam_face_mid = 0,
                                                 const CouplingRegionLayout* fixed_coupling_layout = 0,
                                                 bool fixed_candidate = false) const;
    void set_midpoint_field(EMFields& fields_mid,
                            const EMFields& fields_n,
                            const std::vector<double>& ex_mid,
                            const SpatialGrid& sg,
                            int mpi_rank,
                            int mpi_size) const;
    double integrate_face_work(const std::vector<double>& current_face,
                               const EMFields& fields_mid,
                               const SpatialGrid& sg,
                               double dt) const;

    void exchange_ghosts_x_persistent(Species& sp, const SpatialGrid& sg,
                                      int mpi_rank, int mpi_size) const;

    bool step_diagnostics_enabled_;
    bool beam_enabled_;
    bool low_order_only_;
    bool nonuniform_high_order_enabled_;
    bool fct_enabled_;
    bool capture_midpoint_input_;
    bool legacy_boundary_upwind_high_candidate_for_test_;
    int max_midpoint_iterations_;
    bool midpoint_iteration_trace_for_test_;
    bool fct_activation_audit_enabled_;
    bool controlled_fct_flux_injection_enabled_;
    bool allow_finite_negative_debt_for_test_;
    mutable std::vector<double> ghost_send_left_;
    mutable std::vector<double> ghost_send_right_;
    mutable std::vector<double> ghost_recv_left_;
    mutable std::vector<double> ghost_recv_right_;
    int core_macro_debt_consecutive_steps_;
};

#endif
