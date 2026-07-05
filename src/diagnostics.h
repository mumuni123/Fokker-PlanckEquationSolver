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
                                double W_bkg_E,
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
                                double bkg_energy_residual_step,
                                double bkg_current_max_abs_charge,
                                double bkg_current_max_abs_energy,
                                double bkg_current_max_abs_ampere,
                                double bkg_current_max_abs_charge_minus_ampere,
                                double bkg_current_max_abs_energy_minus_ampere,
                                double bkg_current_e_dot_charge,
                                double bkg_current_e_dot_energy,
                                double bkg_current_e_dot_ampere,
                                double bkg_residual_if_charge_current,
                                double bkg_residual_if_ampere_current,
                                int coupled_iter,
                                double coupled_residual_E,
                                double coupled_residual_J_bkg,
                                double coupled_residual_J_beam,
                                double local_max_loss_u_high,
                                double local_x_at_max_loss_u_high,
                                double local_f_u_max_x,
                                double local_integral_f_u_gt_8_x);

    void write_bkg_stage_diagnostics(int step, double time,
                                     int coupled_iter,
                                     const std::string& stage,
                                     const Species& electrons,
                                     const SpatialGrid& sg,
                                     int mpi_rank, int mpi_size,
                                     double reference_total_mass_raw);

    void write_bkg_stage_negativity(
        int step, double time, int coupled_iter,
        const std::vector<double>& min_f,
        const std::vector<double>& neg_mass,
        const std::vector<long long>& neg_cell_count,
        const std::vector<double>& low_u_neg_mass,
        const std::vector<double>& core_low_u_min_f,
        int mpi_rank);

    void write_bkg_stage_by_u_diagnostics(
        int step, double time, int coupled_iter,
        const std::vector<double>& min_f_core_by_u,
        const std::vector<double>& neg_mass_core_by_u,
        const std::vector<long long>& neg_cell_count_core_by_u,
        const std::vector<double>& min_f_boundary_by_u,
        const std::vector<double>& neg_mass_boundary_by_u,
        const std::vector<long long>& neg_cell_count_boundary_by_u,
        int mpi_rank);

    void write_bkg_low_u_divergence_diagnostics(
        int step, double time, int coupled_iter,
        const std::vector<double>& low_u_neg_added_by_div,
        int mpi_rank);

    // 7.1.6: write per-direction flux-positivity and defect diagnostics
    void write_flux_positivity_diagnostics(
        int step, double time,
        const double flux_pos_min_f_before[3],
        const double flux_pos_min_f_low[3],
        const double flux_pos_min_f_final[3],
        const double flux_pos_low_order_failed[3],
        const double flux_pos_alpha_active[3],
        const double flux_pos_alpha_min[3],
        const double flux_pos_alpha_core[3],
        const double flux_pos_alpha_boundary[3],
        const double flux_pos_neg_mass_prevented[3],
        int mpi_rank);

    void write_stage_flux_defect_diagnostics(
        int step, double time,
        const double defect_mass[3],
        const double defect_momentum[3],
        const double defect_energy[3],
        const double defect_boundary_mass[3],
        const double defect_boundary_energy[3],
        int mpi_rank);

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
    std::ofstream bkg_stage_file;
    std::ofstream bkg_stage_by_u_file;
    std::ofstream bkg_low_u_divergence_file;
    bool debug_enabled;
    bool step_enabled;
    bool has_energy_reference;
    double energy_reference;
    double initial_ke_per_particle_eV;
};

#endif
