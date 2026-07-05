#ifndef VLASOV_AMPERE_MIDPOINT_H
#define VLASOV_AMPERE_MIDPOINT_H

#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"

#include <vector>

class VlasovAmpereMidpointSolver {
public:
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

    struct Result {
        Species species_np1;
        BeamPIC beam_np1;
        EMFields fields_np1;
        std::vector<double> j_bkg_face_mid;
        std::vector<double> j_beam_face_mid;
        std::vector<double> j_total_face_mid;
        std::vector<double> j_bkg_energy_debug_face;
        CurrentDiagnostics current_diag;
        double delta_ke_bkg;
        double delta_ke_beam;
        double field_work_bkg;
        double field_work_beam;
        double energy_residual_bkg;
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
        LowOrderFailureAudit u_low_failure_audit;
        LowOrderFailureAudit mu_low_failure_audit;
        double f_neg_min;
        double f_neg_ratio_max;
        double f_neg_mass_total;
        long long f_neg_cell_count;
        int f_neg_ix;
        int f_neg_iv;
        int f_neg_imu;
        int nonlinear_iterations;
        bool converged;
        bool failed;
        bool soft_accepted;
        bool protected_converged;
        int substeps_used;
    };

    VlasovAmpereMidpointSolver();

    void set_step_diagnostics_enabled(bool enabled) {
        step_diagnostics_enabled_ = enabled;
    }

    Result advance_background_and_fields(const Species& bkg_n,
                                         const BeamPIC& beam_n,
                                         const EMFields& fields_n,
                                         const SpatialGrid& sg,
                                         double dt,
                                         double time,
                                         int mpi_rank,
                                         int mpi_size);

private:
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
        LowOrderFailureAudit u_low_failure_audit;
        LowOrderFailureAudit mu_low_failure_audit;
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

    void reset_result(Result& result) const;
    void reset_current_diag(CurrentDiagnostics& diag) const;
    Result advance_single_step(const Species& bkg_n,
                               const BeamPIC& beam_n,
                               const EMFields& fields_n,
                               const SpatialGrid& sg,
                               double dt,
                               double time,
                               int mpi_rank,
                               int mpi_size,
                               int substeps_used) const;
    Result advance_with_fixed_substeps(const Species& bkg_n,
                                       const BeamPIC& beam_n,
                                       const EMFields& fields_n,
                                       const SpatialGrid& sg,
                                       double dt,
                                       double time,
                                       int mpi_rank,
                                       int mpi_size,
                                       int substeps) const;
    void set_midpoint_field(EMFields& fields_mid,
                            const EMFields& fields_n,
                            const std::vector<double>& ex_mid,
                            const SpatialGrid& sg,
                            int mpi_rank,
                            int mpi_size) const;
    void compute_midpoint_fluxes(const Species& bkg_n,
                                 const Species& bkg_guess_np1,
                                 const EMFields& fields_mid,
                                 const SpatialGrid& sg,
                                 double dt,
                                 int mpi_rank,
                                 int mpi_size,
                                 Species& bkg_new,
                                 FluxPack& fluxes,
                                 const std::vector<double>* alpha_smooth_from,
                                 double alpha_smooth_beta,
                                 bool& finite) const;
    void compute_vlasov_midpoint_residual(const Species& bkg_n,
                                          const Species& bkg_guess_np1,
                                          const EMFields& fields_mid,
                                          const SpatialGrid& sg,
                                          double dt,
                                          int mpi_rank,
                                          int mpi_size,
                                          Species& bkg_new,
                                          FluxPack& fluxes,
                                          const std::vector<double>* alpha_smooth_from,
                                          double alpha_smooth_beta,
                                          bool& finite) const;
    void update_flux_current(const Species& sp,
                             const SpatialGrid& sg,
                             FluxPack& fluxes,
                             Species& bkg_new) const;
    double integrate_face_work(const std::vector<double>& current_face,
                               const EMFields& fields_mid,
                               const SpatialGrid& sg,
                               double dt) const;
    void build_current_diagnostics(const std::vector<double>& j_bkg_charge,
                                   const std::vector<double>& j_bkg_energy,
                                   const std::vector<double>& j_bkg_ampere,
                                   const EMFields& fields_mid,
                                   const SpatialGrid& sg,
                                   double local_delta_ke_bkg,
                                   double dt,
                                   CurrentDiagnostics& diag) const;
    bool check_finite_state(const Species& bkg,
                            const BeamPIC& beam,
                            const EMFields& fields,
                            const std::vector<double>& ex_mid_next,
                            const std::vector<double>& current_next) const;

    void exchange_ghosts_x_persistent(Species& sp, const SpatialGrid& sg,
                                      int mpi_rank, int mpi_size) const;

    bool step_diagnostics_enabled_;
    mutable std::vector<double> ghost_send_left_;
    mutable std::vector<double> ghost_send_right_;
    mutable std::vector<double> ghost_recv_left_;
    mutable std::vector<double> ghost_recv_right_;
};

#endif
