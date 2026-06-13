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
                                double local_max_loss_u_high,
                                double local_x_at_max_loss_u_high,
                                double local_f_u_max_x,
                                double local_integral_f_u_gt_8_x);

    void write_fields(double time,
                      const EMFields& fields,
                      const SpatialGrid& sg,
                      int mpi_rank, int mpi_size);

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
