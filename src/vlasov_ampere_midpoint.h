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
        double region_u_limiter_energy_boundary[2];
        double region_u_limiter_energy_core[2];
        double region_abs_u_limiter_energy_boundary[2];
        double region_abs_u_limiter_energy_core[2];
        double region_limiter_active_fraction_boundary[2];
        double region_limiter_active_fraction_core[2];
        double x_negative_mass_before_repair;
        double x_mass_added_by_positivity_repair;
        double positivity_energy_defect;
        double positivity_mass_defect;
        double u_force_alpha_min;
        double u_force_alpha_active_frac;
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
        std::vector<double> x_final;
        std::vector<double> cell_alpha;
        std::vector<double> j_bkg_face;
        double limiter_active_fraction;
        double limiter_min_alpha;
        double limiter_active_fraction_core;
        double limiter_active_fraction_boundary;
        double limiter_min_alpha_core;
        double limiter_min_alpha_boundary;
        double limiter_energy_defect;
        double limiter_mass_defect;
        double limiter_momentum_defect;
        double region_u_limiter_energy_boundary[2];
        double region_u_limiter_energy_core[2];
        double region_abs_u_limiter_energy_boundary[2];
        double region_abs_u_limiter_energy_core[2];
        double region_limiter_active_fraction_boundary[2];
        double region_limiter_active_fraction_core[2];
        double negative_mass_before_repair;
        double mass_added_by_positivity_repair;
        double positivity_energy_defect;
        double positivity_mass_defect;
        double u_force_alpha_min;
        double u_force_alpha_active_frac;
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
                             const FluxPack& fluxes,
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
