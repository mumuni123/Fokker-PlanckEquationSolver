#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include <fstream>
#include <string>
#include <vector>

class BeamPIC;

class Diagnostics {
public:
    std::string output_dir;
    int snapshot_count;

    Diagnostics();

    void init(const std::string& dir, int mpi_rank,
              bool enable_debug_diagnostics,
              bool enable_step_diagnostics);

    void write_scalars(double time, int step,
                       const Species& electrons,
                       const BeamPIC& beam,
                       const EMFields& fields,
                       double cumulative_collision_energy_delta,
                       int mpi_rank, int mpi_size);

    void write_debug_state(int step, double time,
                           const std::string& stage,
                           const Species& electrons,
                           const BeamPIC& beam,
                           const EMFields& fields,
                           const SpatialGrid& sg,
                           int mpi_rank, int mpi_size,
                           double cfl_v = 0.0,
                           double cfl_mu = 0.0,
                           int nsub_v = 0,
                           int nsub_mu = 0);

    void write_step_diagnostics(int step, double time,
                                bool soft_unconverged,
                                const Species& electrons,
                                const BeamPIC& beam,
                                const EMFields& fields,
                                const SpatialGrid& sg,
                                int mpi_rank, int mpi_size,
                                int nsub_v1,
                                int nsub_mu1,
                                int nsub_v2,
                                int nsub_mu2,
                                double loss_v1,
                                double loss_mu1,
                                double loss_v2,
                                double loss_mu2,
                                double loss_v1_low,
                                double loss_v1_high,
                                double loss_v2_low,
                                double loss_v2_high,
                                double net_Nb_change,
                                double collision_energy_step,
                                double cumulative_collision_energy_delta,
                                double dke_bkg_step,
                                double dke_beam_push,
                                double dE_field_step,
                                 double global_W_bkg_E,
                                double W_beam_E,
                                double v_mass_error_step,
                                double mu_mass_error_step,
                                double v_momentum_delta_step,
                                double mu_momentum_delta_step,
                                double v_energy_delta_step,
                                double mu_energy_delta_step,
                                double E_src_in_step,
                                double E_src_out_step,
                                double E_balance_step,
                                double x_limiter_active_fraction,
                                double x_limiter_min_alpha,
                                 double global_bkg_energy_residual_step,
                                 double global_bkg_current_max_abs_charge,
                                 double global_bkg_current_max_abs_energy,
                                 double global_bkg_current_max_abs_ampere,
                                 double global_bkg_current_max_abs_charge_minus_ampere,
                                 double global_bkg_current_max_abs_energy_minus_ampere,
                                 double global_bkg_current_e_dot_charge,
                                 double global_bkg_current_e_dot_energy,
                                 double global_bkg_current_e_dot_ampere,
                                 double global_bkg_residual_if_charge_current,
                                 double global_bkg_residual_if_ampere_current,
                                double boundary_force_Cu_max,
                                double boundary_force_Cmu_max,
                                int boundary_force_nsub_max,
                                long long boundary_force_remap_cell_count,
                                double boundary_mu_low_L1_before,
                                double boundary_mu_low_L1_after,
                                double boundary_mu_high_L1_after,
                                double J_bkg_neg_boundary,
                                double delta_E_neg_boundary,
                                double boundary_force_remap_mass_loss,
                                double boundary_force_remap_energy_loss,
                                double alpha_interface_BQ_min,
                                double alpha_interface_QC_min,
                                double interface_BQ_flux,
                                double interface_BQ_high_correction,
                                double interface_QC_flux_into_core,
                                double interface_QC_high_correction_into_core,
                                double boundary_energy_diagnostic_invalid,
                                 int global_coupled_iter,
                                 double global_coupled_residual_E,
                                 double global_coupled_residual_J_bkg,
                                 double global_coupled_residual_J_beam,
                                 double global_jn_minus_gstar_je_linf,
                                 double global_stage5_r_fv,
                                 double global_stage5_r_couple,
                                 double local_max_loss_u_high,
                                double local_x_at_max_loss_u_high,
                                double local_f_u_max_x,
                                double local_integral_f_u_gt_8_x);

    void write_fields(double time,
                      const EMFields& fields,
                      const SpatialGrid& sg,
                      int mpi_rank, int mpi_size);

    void write_current_density(double time,
                               const Species& electrons,
                               const BeamPIC& beam,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size,
                               const std::vector<double>* bkg_energy_current_face = 0,
                               const std::vector<double>* bkg_ampere_current_face = 0);

    void write_px_distribution(double time,
                               const Species& sp,
                               int mpi_rank, int mpi_size);

    void write_density_profile(double time,
                               const Species& electrons,
                               const std::vector<double>& beam_density,
                               const std::vector<double>& ion_density_profile,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size);

    void write_electron_distribution(double time,
                                     const Species& electrons,
                                     const SpatialGrid& sg,
                                     int mpi_rank);

    void advance_snapshot();

private:
    std::ofstream scalar_file;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    std::ofstream debug_file;
#endif
    std::ofstream step_file;
    bool debug_enabled;
    bool step_enabled;
    bool has_energy_reference;
    double energy_reference;
    double initial_ke_per_particle_eV;
};

#endif
