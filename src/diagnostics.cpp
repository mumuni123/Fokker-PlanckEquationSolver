#include "diagnostics.h"
#include "beam_pic.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <sys/stat.h>

static void make_output_dir(const std::string& dir)
{
#ifdef _WIN32
    mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif
}

Diagnostics::Diagnostics()
    : snapshot_count(0),
      debug_enabled(false),
      step_enabled(false),
      has_energy_reference(false),
      has_energy_history_reference(false),
      has_last_energy_history_record(false),
      last_energy_history_step(0),
      last_energy_history_time(0.0),
      energy_reference(0.0),
      initial_background_ke(0.0),
      initial_field_energy(0.0),
      initial_background_plus_field_energy(0.0),
      initial_total_energy(0.0),
      initial_ke_per_particle_eV(0.0)
{}

namespace {
void compute_background_boundary_fluxes(const Species& electrons,
                                        const SpatialGrid& sg,
                                        int mpi_rank,
                                        int mpi_size,
                                        double values[6])
{
    for (int i = 0; i < 6; ++i) values[i] = 0.0;
    if (electrons.type != SpeciesType::BACKGROUND_ELECTRON ||
        sg.nx_local <= 0) {
        return;
    }

    const int ng = sg.nghost;
    if (mpi_rank == 0) {
        const size_t xbase = static_cast<size_t>(ng) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = electrons.vgrid.vx_cells[k];
            const int iv = static_cast<int>(k / Param::Nmu);
            const double weight = electrons.vgrid.moment_weight[iv];
            const double current_weight = electrons.vgrid.current_weight[k];
            if (vx > 0.0) {
                values[0] += vx * electrons.f[xbase + k] * weight;
                values[4] += electrons.charge * electrons.f[xbase + k]
                           * current_weight;
            } else if (vx < 0.0) {
                values[1] += -vx * electrons.f[xbase + k] * weight;
            }
        }
    }

    if (mpi_rank == mpi_size - 1) {
        const size_t xbase =
            static_cast<size_t>(ng + sg.nx_local - 1) * Param::Nvmu;
        for (size_t k = 0; k < Param::Nvmu; ++k) {
            const double vx = electrons.vgrid.vx_cells[k];
            const int iv = static_cast<int>(k / Param::Nmu);
            const double weight = electrons.vgrid.moment_weight[iv];
            const double current_weight = electrons.vgrid.current_weight[k];
            if (vx < 0.0) {
                values[2] += -vx * electrons.f[xbase + k] * weight;
                values[5] += electrons.charge * electrons.f[xbase + k]
                           * current_weight;
            } else if (vx > 0.0) {
                values[3] += vx * electrons.f[xbase + k] * weight;
            }
        }
    }
}

double compute_global_mean_rho(const EMFields& fields, const SpatialGrid& sg)
{
    double local_charge = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_charge += fields.rho[ix + sg.nghost] * sg.dx;
    }

    double global_charge = 0.0;
    MPI_Allreduce(&local_charge, &global_charge, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_charge / Param::Lx;
}

double gauss_charge_residual(const EMFields& fields,
                             const SpatialGrid& sg,
                             int ix,
                             double mean_rho)
{
    const int ix_g = ix + sg.nghost;
    double dEx_dx = 0.0;
    if (fields.Ex_face.size() >= static_cast<size_t>(sg.nx_local + 1)) {
        dEx_dx =
            (fields.Ex_face[static_cast<size_t>(ix + 1)]
           - fields.Ex_face[static_cast<size_t>(ix)]) / sg.dx;
    } else {
        dEx_dx =
            (fields.Ex[ix_g + 1] - fields.Ex[ix_g - 1]) / (2.0 * sg.dx);
    }
    return (Const::eps0 * dEx_dx - (fields.rho[ix_g] - mean_rho)) / Const::qe;
}

void gather_max_ex_location(double local_max_abs_Ex,
                            double local_x_at_max_abs_Ex,
                            int mpi_rank,
                            int mpi_size,
                            double& global_max_abs_Ex,
                            double& global_x_at_max_abs_Ex)
{
    double local_pair[2] = { local_max_abs_Ex, local_x_at_max_abs_Ex };
    std::vector<double> gathered;
    if (mpi_rank == 0) gathered.assign(static_cast<size_t>(2 * mpi_size), 0.0);
    MPI_Gather(local_pair, 2, MPI_DOUBLE,
               mpi_rank == 0 ? gathered.data() : static_cast<double*>(0),
               2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        global_max_abs_Ex = 0.0;
        global_x_at_max_abs_Ex = 0.0;
        for (int r = 0; r < mpi_size; ++r) {
            const double value = gathered[static_cast<size_t>(2 * r)];
            const double xpos = gathered[static_cast<size_t>(2 * r + 1)];
            if (value > global_max_abs_Ex) {
                global_max_abs_Ex = value;
                global_x_at_max_abs_Ex = xpos;
            }
        }
    }
}

void gather_max_loss_u_high_location(double local_max_loss,
                                     double local_x,
                                     double local_f_u_max_x,
                                     double local_integral_f_u_gt_8_x,
                                     int mpi_rank,
                                     int mpi_size,
                                     double& global_x,
                                     double& global_f_u_max_x,
                                     double& global_integral_f_u_gt_8_x)
{
    double local_values[4] = {
        local_max_loss,
        local_x,
        local_f_u_max_x,
        local_integral_f_u_gt_8_x
    };
    std::vector<double> gathered;
    if (mpi_rank == 0) gathered.assign(static_cast<size_t>(4 * mpi_size), 0.0);
    MPI_Gather(local_values, 4, MPI_DOUBLE,
               mpi_rank == 0 ? gathered.data() : static_cast<double*>(0),
               4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        double best_loss = -1.0;
        global_x = 0.0;
        global_f_u_max_x = 0.0;
        global_integral_f_u_gt_8_x = 0.0;
        for (int r = 0; r < mpi_size; ++r) {
            const size_t base = static_cast<size_t>(4 * r);
            if (gathered[base] > best_loss) {
                best_loss = gathered[base];
                global_x = gathered[base + 1];
                global_f_u_max_x = gathered[base + 2];
                global_integral_f_u_gt_8_x = gathered[base + 3];
            }
        }
    }
}

}

void Diagnostics::init(const std::string& dir, int mpi_rank,
                       bool enable_debug_diagnostics,
                       bool enable_step_diagnostics,
                       bool enable_accepted_energy_audit)
{
    output_dir = dir;
    snapshot_count = 0;
    step_enabled = enable_step_diagnostics;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    debug_enabled = enable_debug_diagnostics;
#else
    (void)enable_debug_diagnostics;
    debug_enabled = false;
#endif

    if (mpi_rank == 0) {
        make_output_dir(output_dir);
        scalar_file.open((output_dir + "/scalars.dat").c_str());
        scalar_file << "# step  time[fs]  N_bkg_e  KE_bkg_e[J/m2]  "
                    << "N_beam  KE_beam[J/m2]  E_field[J/m2]  "
                    << "E_total[J/m2]  E_beam_injected_cum[J/m2]  "
                    << "E_beam_outflow_cum[J/m2]  E_collision_cum[J/m2]  "
                    << "E_accounted[J/m2]  E_balance_error[J/m2]  "
                    << "Q_total[C/m2]  gauss_residual_int[m^-2]  "
                    << "gauss_residual_abs_max[m^-3]  "
                    << "gauss_residual_abs_max_over_n0  N_bkg_minus_n0  "
                    << "Gamma_bkg_in_left  Gamma_bkg_out_left  "
                    << "Gamma_bkg_in_right  Gamma_bkg_out_right  "
                    << "J_bkg_in_left[A/m2]  J_bkg_in_right[A/m2]  "
                    << "N_beam_in_left_step  N_beam_out_step  "
                    << "J_beam_in_left_step  J_beam_out_step  "
                    << "max_abs_Ex[V/m]  x_at_max_abs_Ex[m]\n";
        scalar_file << std::scientific << std::setprecision(8);
        energy_history_file.open(
            (output_dir + "/background_electron_energy_history.dat").c_str());
        energy_history_file
            << "# Energies are per transverse area. Background electron "
            << "quantities exclude Beam macroparticles.\n"
            << "# Deltas are relative to the first state written in this "
            << "output directory (initial state for a fresh run, restart "
            << "state for a continuation).\n"
            << "# step time_fs N_bkg_e_m-2 KE_bkg_e_J_m-2 "
            << "mean_KE_bkg_e_J mean_KE_bkg_e_eV "
            << "E_field_J_m-2 E_bkg_plus_field_J_m-2 "
            << "KE_beam_J_m-2 E_total_including_beam_J_m-2 "
            << "delta_KE_bkg_e_J_m-2 delta_E_field_J_m-2 "
            << "delta_E_bkg_plus_field_J_m-2 "
            << "delta_E_total_including_beam_J_m-2\n";
        energy_history_file << std::scientific << std::setprecision(16);

#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (debug_enabled) {
            debug_file.open((output_dir + "/debug_diagnostics.dat").c_str());
            debug_file << "# step  time[fs]  stage  max_abs_Ex[V/m]  N_bkg_e  "
                       << "N_beam_macro  N_beam_weighted  CFL_u  CFL_mu  "
                       << "nsub_u  nsub_mu\n";
            debug_file << std::scientific << std::setprecision(8);
        }
#endif
        if (step_enabled) {
            step_file.open((output_dir + "/step_diagnostics.dat").c_str());
            step_file << "# step  time[fs]  accepted  state_advanced  "
                      << "soft_unconverged  "
                      << "max_abs_Ex[V/m]  x_at_max_abs_Ex[m]  "
                      << "gauss_residual_int[m^-2]  gauss_residual_abs_max[m^-3]  "
                      << "gauss_residual_abs_max_over_n0  "
                      << "N_bkg_e  "
                      << "N_beam_macro  N_beam_weighted  "
                      << "beam_cont_abs_l1[A/m2]  "
                      << "beam_cont_abs_linf[A/m3]  "
                      << "beam_cont_rel_injection_l1  "
                      << "beam_cont_rel_injection_linf  "
                      << "beam_boundary_flux_relative_error  "
                      << "beam_trajectory_reconstruction_relative_error  "
                      << "nsub_u1  nsub_mu1  nsub_u2  nsub_mu2  "
                      << "loss_u1  loss_mu1  loss_u2  loss_mu2  "
                      << "loss_u1_low  loss_u1_high  "
                      << "loss_u2_low  loss_u2_high  "
                      << "net_Nb_change[m^-2]  "
                      << "KE_bkg_e[J/m2]  KE_beam[J/m2]  E_field[J/m2]  "
                      << "E_total[J/m2]  dKE_bkg[J/m2]  dKE_beam[J/m2]  "
                       << "dE_field[J/m2]  W_bkg_E_global[J/m2]  W_beam_E[J/m2]  "
                       << "bkg_energy_residual_global[J/m2]  "
                       << "max_abs_J_bkg_charge_global[A/m2]  "
                       << "max_abs_J_bkg_energy_diagnostic_global[A/m2]  "
                       << "max_abs_J_bkg_ampere_global[A/m2]  "
                       << "max_abs_J_bkg_charge_minus_ampere_global[A/m2]  "
                       << "max_abs_J_bkg_energy_minus_ampere_global[A/m2]  "
                       << "E_dot_J_bkg_charge_global[W/m2]  "
                       << "E_dot_J_bkg_energy_diagnostic_global[W/m2]  "
                       << "E_dot_J_bkg_ampere_global[W/m2]  "
                       << "residual_if_charge_current_global[J/m2]  "
                       << "residual_if_ampere_current_global[J/m2]  "
                      << "boundary_force_Cu_max  boundary_force_Cmu_max  "
                      << "boundary_force_nsub_max  "
                      << "boundary_force_remap_cell_count  "
                      << "boundary_mu_low_L1_before  "
                      << "boundary_mu_low_L1_after  "
                      << "boundary_mu_high_L1_after  "
                      << "J_bkg_neg_boundary[A/m2]  "
                      << "delta_E_neg_boundary[V/m]  "
                      << "boundary_force_remap_mass_loss[m^-2]  "
                      << "boundary_force_remap_energy_loss[J/m2]  "
                      << "alpha_interface_BQ_min  alpha_interface_QC_min  "
                      << "interface_BQ_flux  "
                      << "interface_BQ_high_correction  "
                      << "interface_QC_flux_into_core  "
                      << "interface_QC_high_correction_into_core  "
                      << "boundary_energy_diagnostic_invalid  "
                       << "coupled_iter_global  coupled_residual_E_global  "
                       << "coupled_residual_J_bkg_global  coupled_residual_J_beam_global  "
                       << "JN_minus_GstarJE_linf_global[A/m2]  "
                       << "stage5_R_FV_global[J/m2]  stage5_R_couple_global[J/m2]  "
                       << "x_limiter_active_fraction  x_limiter_min_alpha  "
                      << "u_mass_error  mu_mass_error  "
                      << "u_px_delta[kg/m/s/m2]  mu_px_delta[kg/m/s/m2]  "
                      << "u_energy_delta[J/m2]  mu_energy_delta[J/m2]  "
                      << "E_src_in[J/m2]  E_src_out[J/m2]  "
                      << "E_collision_step[J/m2]  "
                      << "E_balance_step[J/m2]  E_beam_injected_cum[J/m2]  "
                      << "E_beam_outflow_cum[J/m2]  E_collision_cum[J/m2]  "
                      << "initial_KE_per_particle_eV  "
                      << "effective_T_eV  "
                      << "x_at_max_loss_u_high[m]  "
                      << "f_u_max_x_mu_avg[u^-3_m^-3]  "
                      << "integral_f_u_gt_8_x[m^-3]  "
                      << "background_coupling_mode\n";
            step_file << std::scientific << std::setprecision(8);

        }
        if (enable_accepted_energy_audit) {
            accepted_energy_ledger_file.open(
                (output_dir + "/accepted_energy_ledger.dat").c_str());
            accepted_energy_ledger_file
                << "# physical_step time_fs dt_s strict soft iterations substeps "
                << "dKE_bkg_actual_J_m2 W_bkg_ampere_J_m2 "
                << "dKE_low_J_m2 dKE_high_J_m2 dKE_final_J_m2 "
                << "dKE_fct_x_J_m2 dKE_fct_u_J_m2 "
                << "dM_fct_x_m2 dM_fct_u_m2 "
                << "dPpar_fct_x_kg_m_s_m2 dPpar_fct_u_kg_m_s_m2 "
                << "B_upar_energy_lower_J_m2 B_upar_energy_upper_J_m2 "
                << "B_uperp_energy_upper_J_m2 "
                << "R_FV_J_m2 R_couple_J_m2 R_unexplained_J_m2 "
                << "R_reconstruction_error_J_m2 flux_telescope_error_J_m2 "
                << "R_total_step_J_m2 R_total_cumulative_J_m2\n";
            velocity_boundary_flux_ledger_file.open(
                (output_dir + "/velocity_boundary_flux_ledger.dat").c_str());
            velocity_boundary_flux_ledger_file
                << "# physical_step time_fs dt_s strict soft substeps "
                << "upar_number_lower_m2 upar_number_upper_m2 "
                << "upar_ppar_lower_kg_m_s_m2 upar_ppar_upper_kg_m_s_m2 "
                << "upar_energy_lower_J_m2 upar_energy_upper_J_m2 "
                << "uperp_number_upper_m2 uperp_ppar_upper_kg_m_s_m2 "
                << "uperp_energy_upper_J_m2\n";
            midpoint_acceptance_ledger_file.open(
                (output_dir + "/midpoint_acceptance_ledger.dat").c_str());
            midpoint_acceptance_ledger_file
                << "# physical_step time_fs iterations strict soft "
                << "residual_E residual_J_bkg residual_f "
                << "last_contraction_E last_contraction_J "
                << "x_fct_active u_fct_active x_min_alpha u_min_alpha\n";
            accepted_energy_ledger_file << std::scientific << std::setprecision(16);
            velocity_boundary_flux_ledger_file << std::scientific << std::setprecision(16);
            midpoint_acceptance_ledger_file << std::scientific << std::setprecision(16);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Diagnostics::write_accepted_energy_ledger(
    long long physical_step, double time, double dt,
    const VlasovAmpereMidpointSolver::Result& midpoint_result,
    double total_energy_residual_step,
    double total_energy_residual_cumulative,
    int mpi_rank)
{
    if (mpi_rank != 0 || !midpoint_result.accepted_energy_ledger.valid) return;
    const VlasovAmpereMidpointSolver::AcceptedEnergyLedger& ledger =
        midpoint_result.accepted_energy_ledger;
    const auto contraction = [](const std::vector<double>& values) {
        if (values.size() < 2) return 1.0;
        const double previous = values[values.size() - 2];
        const double current = values.back();
        if (previous == 0.0) return current == 0.0 ? 1.0 :
            std::numeric_limits<double>::infinity();
        return current / previous;
    };

    if (accepted_energy_ledger_file.is_open()) {
        accepted_energy_ledger_file
            << physical_step << " " << time / Const::femto << " " << dt << " "
            << ledger.strict_accepted << " " << ledger.soft_accepted << " "
            << midpoint_result.nonlinear_iterations << " "
            << ledger.transport_substeps << " "
            << ledger.delta_ke_bkg_actual << " " << ledger.work_ampere_bkg << " "
            << ledger.delta_ke_after_low << " " << ledger.delta_ke_after_high << " "
            << ledger.delta_ke_after_final << " "
            << ledger.delta_ke_fct_x << " " << ledger.delta_ke_fct_u << " "
            << ledger.delta_mass_fct_x << " " << ledger.delta_mass_fct_u << " "
            << ledger.delta_ppar_fct_x << " " << ledger.delta_ppar_fct_u << " "
            << ledger.upar_boundary_energy_lower << " "
            << ledger.upar_boundary_energy_upper << " "
            << ledger.uperp_boundary_energy_upper << " "
            << ledger.residual_fv << " " << ledger.residual_coupling << " "
            << ledger.residual_unexplained << " "
            << ledger.residual_reconstruction_error << " "
            << ledger.flux_telescope_error << " "
            << total_energy_residual_step << " "
            << total_energy_residual_cumulative << "\n";
        accepted_energy_ledger_file.flush();
    }
    if (velocity_boundary_flux_ledger_file.is_open()) {
        velocity_boundary_flux_ledger_file
            << physical_step << " " << time / Const::femto << " " << dt << " "
            << ledger.strict_accepted << " " << ledger.soft_accepted << " "
            << ledger.transport_substeps << " "
            << ledger.upar_boundary_number_lower << " "
            << ledger.upar_boundary_number_upper << " "
            << ledger.upar_boundary_ppar_lower << " "
            << ledger.upar_boundary_ppar_upper << " "
            << ledger.upar_boundary_energy_lower << " "
            << ledger.upar_boundary_energy_upper << " "
            << ledger.uperp_boundary_number_upper << " "
            << ledger.uperp_boundary_ppar_upper << " "
            << ledger.uperp_boundary_energy_upper << "\n";
        velocity_boundary_flux_ledger_file.flush();
    }
    if (midpoint_acceptance_ledger_file.is_open()) {
        midpoint_acceptance_ledger_file
            << physical_step << " " << time / Const::femto << " "
            << midpoint_result.nonlinear_iterations << " "
            << ledger.strict_accepted << " " << ledger.soft_accepted << " "
            << midpoint_result.residual_E << " "
            << midpoint_result.residual_J_bkg << " "
            << midpoint_result.residual_f << " "
            << contraction(midpoint_result.midpoint_residual_e_history) << " "
            << contraction(midpoint_result.midpoint_residual_j_bkg_history) << " "
            << midpoint_result.x_limiter_active_fraction << " "
            << midpoint_result.u_limiter_active_fraction << " "
            << midpoint_result.x_limiter_min_alpha << " "
            << midpoint_result.u_limiter_min_alpha << "\n";
        midpoint_acceptance_ledger_file.flush();
    }
}

void Diagnostics::write_scalars(double time, int step,
                                const Species& electrons,
                                const BeamPIC& beam,
                                const EMFields& fields,
                                double cumulative_collision_energy_delta,
                                int mpi_rank, int mpi_size)
{
    double local_bkg_number = 0.0;
    double local_bkg_ke = 0.0;
    electrons.total_particle_number_and_energy(local_bkg_number, local_bkg_ke);

    double local_q_total = 0.0;
    double local_charge_residual_int = 0.0;
    double local_charge_residual_abs_max = 0.0;
    double local_bkg_delta = 0.0;
    double local_max_abs_Ex = 0.0;
    double local_x_at_max_abs_Ex = 0.0;
    const int ng = electrons.sgrid->nghost;
    const double mean_rho =
        compute_global_mean_rho(fields, *electrons.sgrid);
    for (int ix = 0; ix < electrons.sgrid->nx_local; ++ix) {
        const double charge_residual =
            gauss_charge_residual(fields, *electrons.sgrid, ix, mean_rho);
        local_q_total += fields.rho[ix + ng] * electrons.sgrid->dx;
        local_charge_residual_int += charge_residual * electrons.sgrid->dx;
        local_charge_residual_abs_max =
            std::max(local_charge_residual_abs_max, std::fabs(charge_residual));
        local_bkg_delta +=
            (electrons.number_density[static_cast<size_t>(ix)] - Param::dens) *
            electrons.sgrid->dx;
        const double abs_ex = std::fabs(fields.Ex[ix + ng]);
        if (abs_ex > local_max_abs_Ex) {
            local_max_abs_Ex = abs_ex;
            local_x_at_max_abs_Ex = electrons.sgrid->x(ix + ng);
        }
    }

    double local_boundary[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double global_boundary[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    compute_background_boundary_fluxes(electrons, *electrons.sgrid,
                                       mpi_rank, mpi_size, local_boundary);

    double local_values[15] = {
        local_bkg_number,
        local_bkg_ke,
        beam.total_particle_number(*electrons.sgrid),
        beam.total_kinetic_energy(),
        fields.total_energy(),
        beam.cumulative_injected_energy(),
        beam.cumulative_outflow_energy(),
        cumulative_collision_energy_delta,
        local_q_total,
        local_charge_residual_int,
        local_bkg_delta,
        beam.last_injected_number(),
        beam.last_outflow_number(),
        beam.last_injected_current(),
        beam.last_outflow_current()
    };
    double global_values[15] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    double global_charge_residual_abs_max = 0.0;

    MPI_Reduce(local_values, global_values, 15, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    MPI_Reduce(&local_charge_residual_abs_max,
               &global_charge_residual_abs_max, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_boundary, global_boundary, 6, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);

    double global_max_abs_Ex = 0.0;
    double global_x_at_max_abs_Ex = 0.0;
    gather_max_ex_location(local_max_abs_Ex, local_x_at_max_abs_Ex,
                           mpi_rank, mpi_size,
                           global_max_abs_Ex, global_x_at_max_abs_Ex);

    if (mpi_rank == 0) {
        const double total_energy =
            global_values[1] + global_values[3] + global_values[4];
        const double background_plus_field_energy =
            global_values[1] + global_values[4];
        const double mean_background_ke =
            global_values[0] > 0.0
            ? global_values[1] / global_values[0] : 0.0;
        const double accounted_energy =
            total_energy - global_values[5] + global_values[6]
            - global_values[7];
        if (!has_energy_reference) {
            energy_reference = accounted_energy;
            has_energy_reference = true;
        }
        const double balance_error = accounted_energy - energy_reference;
        if (step == 0 && global_values[0] > 0.0) {
            initial_ke_per_particle_eV =
                global_values[1] / global_values[0] / Const::eV;
        }
        if (!has_energy_history_reference) {
            initial_background_ke = global_values[1];
            initial_field_energy = global_values[4];
            initial_background_plus_field_energy =
                background_plus_field_energy;
            initial_total_energy = total_energy;
            has_energy_history_reference = true;
        }

        scalar_file << step << "  "
                    << time / Const::femto << "  "
                    << global_values[0] << "  "
                    << global_values[1] << "  "
                    << global_values[2] << "  "
                    << global_values[3] << "  "
                    << global_values[4] << "  "
                    << total_energy << "  "
                    << global_values[5] << "  "
                    << global_values[6] << "  "
                    << global_values[7] << "  "
                    << accounted_energy << "  "
                    << balance_error << "  "
                    << global_values[8] << "  "
                    << global_values[9] << "  "
                    << global_charge_residual_abs_max << "  "
                    << global_charge_residual_abs_max / Param::dens << "  "
                    << global_values[10] << "  "
                    << global_boundary[0] << "  "
                    << global_boundary[1] << "  "
                    << global_boundary[2] << "  "
                    << global_boundary[3] << "  "
                    << global_boundary[4] << "  "
                    << global_boundary[5] << "  "
                    << global_values[11] << "  "
                    << global_values[12] << "  "
                    << global_values[13] << "  "
                    << global_values[14] << "  "
                    << global_max_abs_Ex << "  "
                    << global_x_at_max_abs_Ex << "\n";
        scalar_file.flush();

        if (!has_last_energy_history_record ||
            step != last_energy_history_step ||
            time != last_energy_history_time) {
            energy_history_file
                << step << " "
                << time / Const::femto << " "
                << global_values[0] << " "
                << global_values[1] << " "
                << mean_background_ke << " "
                << mean_background_ke / Const::eV << " "
                << global_values[4] << " "
                << background_plus_field_energy << " "
                << global_values[3] << " "
                << total_energy << " "
                << global_values[1] - initial_background_ke << " "
                << global_values[4] - initial_field_energy << " "
                << background_plus_field_energy -
                       initial_background_plus_field_energy << " "
                << total_energy - initial_total_energy << "\n";
            energy_history_file.flush();
            has_last_energy_history_record = true;
            last_energy_history_step = step;
            last_energy_history_time = time;
        }
    }
    (void)mpi_size;
}

void Diagnostics::write_debug_state(int step, double time,
                                    const std::string& stage,
                                    const Species& electrons,
                                    const BeamPIC& beam,
                                    const EMFields& fields,
                                    const SpatialGrid& sg,
                                    int mpi_rank, int mpi_size,
                                    double cfl_v,
                                    double cfl_mu,
                                    int nsub_v,
                                    int nsub_mu)
{
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (!debug_enabled) return;

    double local_max_abs_Ex = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_max_abs_Ex = std::max(local_max_abs_Ex,
                                    std::fabs(fields.Ex[ix + sg.nghost]));
    }

    double local_N_bkg_e = electrons.total_particle_number();
    double local_N_beam_macro = static_cast<double>(beam.particles.size());
    double local_N_beam_weighted = beam.total_particle_number(sg);

    double global_max_abs_Ex = 0.0;
    double global_N_bkg_e = 0.0;
    double global_N_beam_macro = 0.0;
    double global_N_beam_weighted = 0.0;
    double global_cfl_v = 0.0;
    double global_cfl_mu = 0.0;
    int global_nsub_v = 0;
    int global_nsub_mu = 0;

    MPI_Reduce(&local_max_abs_Ex, &global_max_abs_Ex, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_bkg_e, &global_N_bkg_e, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_macro, &global_N_beam_macro, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_weighted, &global_N_beam_weighted, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cfl_v, &global_cfl_v, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cfl_mu, &global_cfl_mu, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nsub_v, &global_nsub_v, 1,
               MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nsub_mu, &global_nsub_mu, 1,
               MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        debug_file << step << "  "
                   << time / Const::femto << "  "
                   << stage << "  "
                   << global_max_abs_Ex << "  "
                   << global_N_bkg_e << "  "
                   << global_N_beam_macro << "  "
                   << global_N_beam_weighted << "  "
                   << global_cfl_v << "  "
                   << global_cfl_mu << "  "
                   << global_nsub_v << "  "
                   << global_nsub_mu << "\n";
        debug_file.flush();
    }
#else
    (void)step;
    (void)time;
    (void)stage;
    (void)electrons;
    (void)beam;
    (void)fields;
    (void)sg;
    (void)mpi_rank;
    (void)mpi_size;
    (void)cfl_v;
    (void)cfl_mu;
    (void)nsub_v;
    (void)nsub_mu;
#endif
}

void Diagnostics::write_step_diagnostics(int step, double time,
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
                                         double local_integral_f_u_gt_8_x,
                                         int background_coupling_mode)
{
    if (!step_enabled) return;

    double local_max_abs_Ex = 0.0;
    double local_x_at_max_abs_Ex = 0.0;
    double local_charge_residual_int = 0.0;
    double local_charge_residual_abs_max = 0.0;
    const double mean_rho = compute_global_mean_rho(fields, sg);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int ix_g = ix + sg.nghost;
        const double charge_residual =
            gauss_charge_residual(fields, sg, ix, mean_rho);
        local_charge_residual_int += charge_residual * sg.dx;
        local_charge_residual_abs_max =
            std::max(local_charge_residual_abs_max, std::fabs(charge_residual));
        const double abs_ex = std::fabs(fields.Ex[ix_g]);
        if (abs_ex > local_max_abs_Ex) {
            local_max_abs_Ex = abs_ex;
            local_x_at_max_abs_Ex = sg.x(ix_g);
        }
    }

    const double local_N_bkg_e = electrons.total_particle_number();
    const double local_N_beam_macro = static_cast<double>(beam.particles.size());
    const double local_N_beam_weighted = beam.total_particle_number(sg);
    double local_beam_continuity_sum[2] = {
        beam.last_continuity_abs_l1_residual(),
        beam.last_continuity_l1_error()
    };
    double local_beam_continuity_max[4] = {
        beam.last_continuity_abs_linf_residual(),
        beam.last_continuity_linf_error(),
        beam.last_boundary_flux_error(),
        beam.last_trajectory_reconstruction_error()
    };
    double global_beam_continuity_sum[2] = { 0.0, 0.0 };
    double global_beam_continuity_max[4] = { 0.0, 0.0, 0.0, 0.0 };

    double local_losses[9] = {
        loss_v1, loss_mu1, loss_v2, loss_mu2,
        loss_v1_low, loss_v1_high, loss_v2_low, loss_v2_high,
        net_Nb_change
    };
    double global_losses[9] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const double local_bkg_ke = electrons.total_kinetic_energy();
    const double local_beam_ke = beam.total_kinetic_energy();
    const double local_field_energy = fields.total_energy();
    const double local_total_energy =
        local_bkg_ke + local_beam_ke + local_field_energy;
    // Entries sourced from Species/Beam/EMFields are rank-local and need one
    // SUM reduction.  Solver closure diagnostics are already global and are
    // deliberately kept out of this array.
    double local_energy[28] = {
        local_bkg_ke,
        local_beam_ke,
        local_field_energy,
        local_total_energy,
        dke_bkg_step,
        dke_beam_push,
        dE_field_step,
        0.0,
        W_beam_E,
        v_mass_error_step,
        mu_mass_error_step,
        v_momentum_delta_step,
        mu_momentum_delta_step,
        v_energy_delta_step,
        mu_energy_delta_step,
        E_src_in_step,
        E_src_out_step,
        collision_energy_step,
        E_balance_step,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        beam.cumulative_injected_energy(),
        beam.cumulative_outflow_energy(),
        cumulative_collision_energy_delta
    };
    double global_energy[28] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };
    // One collective for all rank-local SUM quantities.  The layout is kept
    // contiguous to reduce synchronization overhead without mixing in any
    // solver-global diagnostics.
    double local_sum_values[43] = {0.0};
    double global_sum_values[43] = {0.0};
    local_sum_values[0] = local_charge_residual_int;
    local_sum_values[1] = local_N_bkg_e;
    local_sum_values[2] = local_N_beam_macro;
    local_sum_values[3] = local_N_beam_weighted;
    local_sum_values[4] = local_beam_continuity_sum[0];
    local_sum_values[5] = local_beam_continuity_sum[1];
    for (int i = 0; i < 9; ++i)
        local_sum_values[6 + i] = local_losses[i];
    for (int i = 0; i < 28; ++i)
        local_sum_values[15 + i] = local_energy[i];
    const double global_solver_energy[7] = {
        global_W_bkg_E,
        global_bkg_energy_residual_step,
        global_bkg_current_e_dot_charge,
        global_bkg_current_e_dot_energy,
        global_bkg_current_e_dot_ampere,
        global_bkg_residual_if_charge_current,
        global_bkg_residual_if_ampere_current
    };
    double global_max_abs_Ex = 0.0;
    double global_x_at_max_abs_Ex = 0.0;
    double global_charge_residual_int = 0.0;
    double global_charge_residual_abs_max = 0.0;
    double global_N_bkg_e = 0.0;
    double global_N_beam_macro = 0.0;
    double global_N_beam_weighted = 0.0;
    const double global_bkg_current_max_values[5] = {
        global_bkg_current_max_abs_charge,
        global_bkg_current_max_abs_energy,
        global_bkg_current_max_abs_ampere,
        global_bkg_current_max_abs_charge_minus_ampere,
        global_bkg_current_max_abs_energy_minus_ampere
    };
    double global_x_at_max_loss_u_high = 0.0;
    double global_f_u_max_x = 0.0;
    double global_integral_f_u_gt_8_x = 0.0;
    int local_nsub[4] = { nsub_v1, nsub_mu1, nsub_v2, nsub_mu2 };
    int global_nsub[4] = { 0, 0, 0, 0 };

    gather_max_ex_location(local_max_abs_Ex, local_x_at_max_abs_Ex,
                           mpi_rank, mpi_size,
                           global_max_abs_Ex, global_x_at_max_abs_Ex);
    gather_max_loss_u_high_location(local_max_loss_u_high,
                                    local_x_at_max_loss_u_high,
                                    local_f_u_max_x,
                                    local_integral_f_u_gt_8_x,
                                    mpi_rank, mpi_size,
                                    global_x_at_max_loss_u_high,
                                    global_f_u_max_x,
                                    global_integral_f_u_gt_8_x);
    MPI_Reduce(&local_charge_residual_abs_max,
               &global_charge_residual_abs_max, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_beam_continuity_max, global_beam_continuity_max,
               4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_nsub, global_nsub, 4, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_sum_values, global_sum_values, 43, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        global_charge_residual_int = global_sum_values[0];
        global_N_bkg_e = global_sum_values[1];
        global_N_beam_macro = global_sum_values[2];
        global_N_beam_weighted = global_sum_values[3];
        global_beam_continuity_sum[0] = global_sum_values[4];
        global_beam_continuity_sum[1] = global_sum_values[5];
        for (int i = 0; i < 9; ++i)
            global_losses[i] = global_sum_values[6 + i];
        for (int i = 0; i < 28; ++i)
            global_energy[i] = global_sum_values[15 + i];
        const double current_ke_per_particle_eV =
            (global_N_bkg_e > 0.0)
            ? global_energy[0] / global_N_bkg_e / Const::eV
            : 0.0;
        const double effective_T_eV =
            (2.0 / 3.0) * current_ke_per_particle_eV;

        step_file << step << "  "
                  << time / Const::femto << "  "
                  << 1 << "  "
                  << 1 << "  "
                  << (soft_unconverged ? 1 : 0) << "  "
                  << global_max_abs_Ex << "  "
                  << global_x_at_max_abs_Ex << "  "
                  << global_charge_residual_int << "  "
                  << global_charge_residual_abs_max << "  "
                  << global_charge_residual_abs_max / Param::dens << "  "
                  << global_N_bkg_e << "  "
                  << global_N_beam_macro << "  "
                  << global_N_beam_weighted << "  "
                  << global_beam_continuity_sum[0] << "  "
                  << global_beam_continuity_max[0] << "  "
                  << global_beam_continuity_sum[1] << "  "
                  << global_beam_continuity_max[1] << "  "
                  << global_beam_continuity_max[2] << "  "
                  << global_beam_continuity_max[3] << "  "
                  << global_nsub[0] << "  "
                  << global_nsub[1] << "  "
                  << global_nsub[2] << "  "
                  << global_nsub[3] << "  "
                  << global_losses[0] << "  "
                  << global_losses[1] << "  "
                  << global_losses[2] << "  "
                  << global_losses[3] << "  "
                  << global_losses[4] << "  "
                  << global_losses[5] << "  "
                  << global_losses[6] << "  "
                  << global_losses[7] << "  "
                  << global_losses[8] << "  "
                  << global_energy[0] << "  "
                  << global_energy[1] << "  "
                  << global_energy[2] << "  "
                  << global_energy[3] << "  "
                  << global_energy[4] << "  "
                  << global_energy[5] << "  "
                  << global_energy[6] << "  "
                   << global_solver_energy[0] << "  "
                  << global_energy[8] << "  "
                   << global_solver_energy[1] << "  "
                  << global_bkg_current_max_values[0] << "  "
                  << global_bkg_current_max_values[1] << "  "
                  << global_bkg_current_max_values[2] << "  "
                  << global_bkg_current_max_values[3] << "  "
                  << global_bkg_current_max_values[4] << "  "
                   << global_solver_energy[2] << "  "
                   << global_solver_energy[3] << "  "
                   << global_solver_energy[4] << "  "
                   << global_solver_energy[5] << "  "
                   << global_solver_energy[6] << "  "
                  << boundary_force_Cu_max << "  "
                  << boundary_force_Cmu_max << "  "
                  << boundary_force_nsub_max << "  "
                  << boundary_force_remap_cell_count << "  "
                  << boundary_mu_low_L1_before << "  "
                  << boundary_mu_low_L1_after << "  "
                  << boundary_mu_high_L1_after << "  "
                  << J_bkg_neg_boundary << "  "
                  << delta_E_neg_boundary << "  "
                  << boundary_force_remap_mass_loss << "  "
                  << boundary_force_remap_energy_loss << "  "
                  << alpha_interface_BQ_min << "  "
                  << alpha_interface_QC_min << "  "
                  << interface_BQ_flux << "  "
                  << interface_BQ_high_correction << "  "
                  << interface_QC_flux_into_core << "  "
                  << interface_QC_high_correction_into_core << "  "
                  << boundary_energy_diagnostic_invalid << "  "
                   << global_coupled_iter << "  "
                   << global_coupled_residual_E << "  "
                   << global_coupled_residual_J_bkg << "  "
                   << global_coupled_residual_J_beam << "  "
                   << global_jn_minus_gstar_je_linf << "  "
                   << global_stage5_r_fv << "  "
                   << global_stage5_r_couple << "  "
                   << x_limiter_active_fraction << "  "
                  << x_limiter_min_alpha << "  "
                  << global_energy[9] << "  "
                  << global_energy[10] << "  "
                  << global_energy[11] << "  "
                  << global_energy[12] << "  "
                  << global_energy[13] << "  "
                  << global_energy[14] << "  "
                  << global_energy[15] << "  "
                  << global_energy[16] << "  "
                  << global_energy[17] << "  "
                  << global_energy[18] << "  "
                  << global_energy[25] << "  "
                  << global_energy[26] << "  "
                  << global_energy[27] << "  "
                  << initial_ke_per_particle_eV << "  "
                  << effective_T_eV << "  "
                  << global_x_at_max_loss_u_high << "  "
                  << global_f_u_max_x << "  "
                  << global_integral_f_u_gt_8_x << "  "
                  << background_coupling_mode << "\n";
        step_file.flush();
    }
}

void Diagnostics::write_fields(double time,
                               const EMFields& fields,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size)
{
    int ng = sg.nghost;
    int nxl = sg.nx_local;

    std::vector<double> local_Ex(nxl);
    std::vector<double> local_phi(nxl);
    for (int i = 0; i < nxl; ++i) local_Ex[i] = fields.Ex[i + ng];
    for (int i = 0; i < nxl; ++i) local_phi[i] = fields.phi[i + ng];

    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);
    MPI_Gather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) displs[r] = displs[r - 1] + counts[r - 1];
    }

    std::vector<double> global_Ex(sg.nx_global);
    std::vector<double> global_phi(sg.nx_global);
    MPI_Gatherv(local_Ex.data(), nxl, MPI_DOUBLE,
                global_Ex.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(local_phi.data(), nxl, MPI_DOUBLE,
                global_phi.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fields_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# x[um]  Ex[V/m]  phi[V]\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i < sg.nx_global; ++i) {
            out << (i + 0.5) * sg.dx / Const::micro << "  "
                << global_Ex[i] << "  "
                << global_phi[i] << "\n";
        }
    }

    const int local_face_count = nxl + ((mpi_rank == 0) ? 1 : 0);
    std::vector<double> local_Ex_face(static_cast<size_t>(local_face_count), 0.0);
    for (int i = 0; i < local_face_count; ++i) {
        const int local_face = i + ((mpi_rank == 0) ? 0 : 1);
        const size_t slot = static_cast<size_t>(local_face);
        local_Ex_face[static_cast<size_t>(i)] =
            (slot < fields.Ex_face.size()) ? fields.Ex_face[slot] : 0.0;
    }

    MPI_Gather(&local_face_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) displs[r] = displs[r - 1] + counts[r - 1];
    }

    std::vector<double> global_Ex_face(sg.nx_global + 1);
    MPI_Gatherv(local_Ex_face.data(), local_face_count, MPI_DOUBLE,
                global_Ex_face.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fields_face_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# time[fs] = " << time / Const::femto << "\n";
        out << "# x_face[um]  Ex_face[V/m]\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i <= sg.nx_global; ++i) {
            out << i * sg.dx / Const::micro << "  "
                << global_Ex_face[i] << "\n";
        }
    }
}

void Diagnostics::write_current_density(double time,
                                        const Species& electrons,
                                        const BeamPIC& beam,
                                        const SpatialGrid& sg,
                                        int mpi_rank, int mpi_size,
                                        const std::vector<double>* bkg_energy_current_face,
                                        const std::vector<double>* bkg_ampere_current_face,
                                        bool bkg_energy_current_valid)
{
    int nxl = sg.nx_local;
    const int local_face_count = nxl + ((mpi_rank == 0) ? 1 : 0);
    const int ncomp = 6;
    std::vector<double> local_J(
        static_cast<size_t>(local_face_count) * ncomp, 0.0);
    for (int i = 0; i < local_face_count; ++i) {
        const int local_face = i + ((mpi_rank == 0) ? 0 : 1);
        const size_t slot = static_cast<size_t>(local_face);
        const size_t base = static_cast<size_t>(i) * ncomp;
        const double J_bkg = (slot < electrons.current_face_x.size())
                           ? electrons.current_face_x[slot] : 0.0;
        const double J_bkg_energy =
            (bkg_energy_current_face && slot < bkg_energy_current_face->size())
            ? (*bkg_energy_current_face)[slot] : 0.0;
        const double J_bkg_ampere =
            (bkg_ampere_current_face && slot < bkg_ampere_current_face->size())
            ? (*bkg_ampere_current_face)[slot] : J_bkg;
        const double J_beam = (slot < beam.current_face_x.size())
                            ? beam.current_face_x[slot] : 0.0;
        local_J[base] = J_bkg;
        local_J[base + 1] = J_bkg_energy;
        local_J[base + 2] = J_bkg - J_bkg_ampere;
        local_J[base + 3] = J_bkg_ampere;
        local_J[base + 4] = J_beam;
        local_J[base + 5] = J_bkg_ampere + J_beam;
    }

    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);
    MPI_Gather(&local_face_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) displs[r] = displs[r - 1] + counts[r - 1];
        for (int r = 0; r < mpi_size; ++r) {
            counts[r] *= ncomp;
            displs[r] *= ncomp;
        }
    }

    std::vector<double> global_J(
        static_cast<size_t>(sg.nx_global + 1) * ncomp, 0.0);
    MPI_Gatherv(local_J.data(), local_face_count * ncomp, MPI_DOUBLE,
                global_J.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/current_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# time[fs] = " << time / Const::femto << "\n";
        out << "# J_bkg_energy_diagnostic_valid = "
            << (bkg_energy_current_valid ? 1 : 0)
            << " (0 means disabled diagnostic, not zero Ampere current)\n";
        out << "# x_face[um]  J_bkg_charge_face[A/m2]  "
            << "J_bkg_energy_diagnostic_face[A/m2]  "
            << "J_bkg_charge_minus_ampere[A/m2]  "
            << "J_bkg_ampere_face[A/m2]  J_beam_face[A/m2]  "
            << "J_total_charge_face[A/m2]  J_total_ampere_face[A/m2]\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i <= sg.nx_global; ++i) {
            const size_t base = static_cast<size_t>(i) * ncomp;
            out << i * sg.dx / Const::micro
                << "  " << global_J[base]
                << "  " << global_J[base + 1]
                << "  " << global_J[base + 2]
                << "  " << global_J[base + 3]
                << "  " << global_J[base + 4]
                << "  " << global_J[base] + global_J[base + 4]
                << "  " << global_J[base + 5] << "\n";
        }
    }
}

void Diagnostics::write_px_distribution(double time,
                                        const Species& sp,
                                        int mpi_rank, int mpi_size)
{
    int ng = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;
    (void)mpi_size;

    if (sp.cylindrical_mass_representation) {
        const int energy_bin_count = 256;
        std::vector<double> energy_edges(
            static_cast<size_t>(energy_bin_count + 1), 0.0);
        double minimum_positive_energy_eV =
            std::numeric_limits<double>::infinity();
        double maximum_energy_eV = 0.0;
        for (size_t velocity = 0;
             velocity < sp.cgrid.kinetic_energy.size(); ++velocity) {
            const double energy_eV =
                sp.cgrid.kinetic_energy[velocity] / Const::eV;
            if (energy_eV > 0.0)
                minimum_positive_energy_eV =
                    std::min(minimum_positive_energy_eV, energy_eV);
            maximum_energy_eV = std::max(maximum_energy_eV, energy_eV);
        }
        if (!std::isfinite(minimum_positive_energy_eV) ||
            !(maximum_energy_eV > 0.0)) {
            minimum_positive_energy_eV = 1.0e-12;
            maximum_energy_eV = 1.0;
        }
        const double logarithmic_minimum =
            std::max(0.5 * minimum_positive_energy_eV,
                     std::numeric_limits<double>::min());
        const double logarithmic_maximum =
            maximum_energy_eV *
            (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
        energy_edges[0] = 0.0;
        for (int edge = 1; edge <= energy_bin_count; ++edge) {
            const double fraction =
                static_cast<double>(edge - 1) /
                static_cast<double>(energy_bin_count - 1);
            energy_edges[static_cast<size_t>(edge)] =
                logarithmic_minimum *
                std::pow(logarithmic_maximum / logarithmic_minimum,
                         fraction);
        }
        std::vector<int> energy_bin_by_velocity(Param::Nvmu, 0);
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t velocity = idx2(j, k);
                const double energy_eV =
                    sp.cgrid.kinetic_energy[velocity] / Const::eV;
                const std::vector<double>::const_iterator upper =
                    std::upper_bound(energy_edges.begin(),
                                     energy_edges.end(), energy_eV);
                int bin = static_cast<int>(
                    upper - energy_edges.begin()) - 1;
                energy_bin_by_velocity[velocity] = std::max(
                    0, std::min(energy_bin_count - 1, bin));
            }
        }

        std::vector<double> local_parallel(
            static_cast<size_t>(Param::Nv), 0.0);
        std::vector<double> local_perpendicular(
            static_cast<size_t>(Param::Nmu), 0.0);
        std::vector<double> local_energy(
            static_cast<size_t>(energy_bin_count), 0.0);
        std::vector<double> local_energy_sum(
            static_cast<size_t>(energy_bin_count), 0.0);
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double particle_number =
                        sp.f[idx3(ix_g, j, k)];
                    local_parallel[static_cast<size_t>(j)] +=
                        particle_number;
                    local_perpendicular[static_cast<size_t>(k)] +=
                        particle_number;
                    const int bin =
                        energy_bin_by_velocity[idx2(j, k)];
                    local_energy[static_cast<size_t>(bin)] +=
                        particle_number;
                    local_energy_sum[static_cast<size_t>(bin)] +=
                        particle_number *
                        sp.cgrid.kinetic_energy[idx2(j, k)];
                }
            }
        }

        std::vector<double> global_parallel(
            static_cast<size_t>(Param::Nv), 0.0);
        std::vector<double> global_perpendicular(
            static_cast<size_t>(Param::Nmu), 0.0);
        std::vector<double> global_energy(
            static_cast<size_t>(energy_bin_count), 0.0);
        std::vector<double> global_energy_sum(
            static_cast<size_t>(energy_bin_count), 0.0);
        MPI_Reduce(local_parallel.data(), global_parallel.data(), Param::Nv,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_perpendicular.data(), global_perpendicular.data(),
                   Param::Nmu, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(local_energy.data(), global_energy.data(),
                   energy_bin_count, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(local_energy_sum.data(), global_energy_sum.data(),
                   energy_bin_count, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            double spectrum_number = 0.0;
            double spectrum_energy = 0.0;
            for (int bin = 0; bin < energy_bin_count; ++bin) {
                spectrum_number += global_energy[static_cast<size_t>(bin)];
                spectrum_energy +=
                    global_energy_sum[static_cast<size_t>(bin)];
            }

            std::ostringstream energy_name;
            energy_name << output_dir << "/energy_spectrum_" << sp.name
                        << "_" << std::setw(5) << std::setfill('0')
                        << snapshot_count << ".dat";
            std::ofstream energy_out(energy_name.str().c_str());
            energy_out << std::scientific << std::setprecision(16);
            energy_out << "# time_fs " << time / Const::femto << "\n"
                       << "# Background-electron kinetic-energy spectrum; "
                       << "N_bin integrates to N_bkg_e per transverse area.\n"
                       << "# integrated_N_m-2 " << spectrum_number << "\n"
                       << "# integrated_exact_kinetic_energy_J_m-2 "
                       << spectrum_energy << "\n"
                       << "# E_left_eV E_geometric_center_eV "
                       << "E_mean_in_bin_eV E_right_eV "
                       << "N_bin_m-2 dN_dE_m-2_eV-1 fraction\n";
            for (int bin = 0; bin < energy_bin_count; ++bin) {
                const double left =
                    energy_edges[static_cast<size_t>(bin)];
                const double right =
                    energy_edges[static_cast<size_t>(bin + 1)];
                const double center = left > 0.0
                    ? std::sqrt(left * right) : 0.5 * right;
                const double width = right - left;
                const double count =
                    global_energy[static_cast<size_t>(bin)];
                const double mean_energy_eV = count > 0.0
                    ? global_energy_sum[static_cast<size_t>(bin)] /
                      count / Const::eV : center;
                energy_out << left << " " << center << " "
                           << mean_energy_eV << " " << right << " "
                           << count << " "
                           << (width > 0.0 ? count / width : 0.0) << " "
                           << (spectrum_number > 0.0
                               ? count / spectrum_number : 0.0) << "\n";
            }

            std::ostringstream parallel_name;
            parallel_name << output_dir << "/momentum_parallel_" << sp.name
                          << "_" << std::setw(5) << std::setfill('0')
                          << snapshot_count << ".dat";
            std::ofstream parallel_out(parallel_name.str().c_str());
            parallel_out
                << "# time_fs " << time / Const::femto << "\n"
                << "# u_parallel=p_parallel/(m_e*c); N_bin integrates "
                << "to N_bkg_e per transverse area.\n"
                << "# u_left u_center u_right p_center_kg_m_s-1 "
                << "N_bin_m-2 dN_du_m-2 "
                << "dN_dp_m-2_per_kg_m_s\n";
            parallel_out << std::scientific << std::setprecision(16);
            for (int j = 0; j < Param::Nv; ++j) {
                const double left = sp.cgrid.upar_faces[j];
                const double right = sp.cgrid.upar_faces[j + 1];
                const double width = right - left;
                const double count =
                    global_parallel[static_cast<size_t>(j)];
                parallel_out
                    << left << " " << sp.cgrid.upar_cells[j] << " "
                    << right << " "
                    << sp.mass * Const::c * sp.cgrid.upar_cells[j] << " "
                    << count << " "
                    << (width > 0.0 ? count / width : 0.0) << " "
                    << (width > 0.0
                        ? count / (sp.mass * Const::c * width) : 0.0)
                    << "\n";
            }

            std::ostringstream perpendicular_name;
            perpendicular_name
                << output_dir << "/momentum_perpendicular_" << sp.name
                << "_" << std::setw(5) << std::setfill('0')
                << snapshot_count << ".dat";
            std::ofstream perpendicular_out(
                perpendicular_name.str().c_str());
            perpendicular_out
                << "# time_fs " << time / Const::femto << "\n"
                << "# u_perp=p_perp/(m_e*c); cylindrical ring measure is "
                << "already integrated into N_bin.\n"
                << "# u_left u_center u_right p_center_kg_m_s-1 "
                << "N_bin_m-2 dN_du_m-2 "
                << "dN_dp_m-2_per_kg_m_s\n";
            perpendicular_out << std::scientific
                              << std::setprecision(16);
            for (int k = 0; k < Param::Nmu; ++k) {
                const double left = sp.cgrid.uperp_faces[k];
                const double right = sp.cgrid.uperp_faces[k + 1];
                const double width = right - left;
                const double count =
                    global_perpendicular[static_cast<size_t>(k)];
                perpendicular_out
                    << left << " " << sp.cgrid.uperp_cells[k] << " "
                    << right << " "
                    << sp.mass * Const::c * sp.cgrid.uperp_cells[k] << " "
                    << count << " "
                    << (width > 0.0 ? count / width : 0.0) << " "
                    << (width > 0.0
                        ? count / (sp.mass * Const::c * width) : 0.0)
                    << "\n";
            }
        }
        return;
    }

    std::vector<double> local_Fu(Param::Nv, 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double u = sp.vgrid.v(iv);
            double sum = 0.0;
            size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                sum += sp.f[row + imu];
            }
            local_Fu[iv] += sum * 2.0 * Const::pi * u * u
                          * sp.vgrid.dmu * sp.sgrid->dx;
        }
    }

    std::vector<double> global_Fu(Param::Nv, 0.0);
    MPI_Reduce(local_Fu.data(), global_Fu.data(), Param::Nv,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fv_" << sp.name << "_"
              << std::setw(5) << std::setfill('0') << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# u[p/(m c)]  F(u)\n";
        out << std::scientific << std::setprecision(8);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            out << sp.vgrid.v(iv) << "  " << global_Fu[iv] << "\n";
        }
    }
}

void Diagnostics::write_density_profile(double time,
                                        const Species& electrons,
                                        const std::vector<double>& beam_density,
                                        const std::vector<double>& ion_density_profile,
                                        const SpatialGrid& sg,
                                        int mpi_rank, int mpi_size)
{
    int nxl = sg.nx_local;
    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);
    MPI_Gather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) displs[r] = displs[r - 1] + counts[r - 1];
    }

    std::vector<double> global_ne(sg.nx_global);
    std::vector<double> global_nb(sg.nx_global);
    std::vector<double> global_zni(sg.nx_global);
    MPI_Gatherv(electrons.number_density.data(), nxl, MPI_DOUBLE,
                global_ne.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(beam_density.data(), nxl, MPI_DOUBLE,
                global_nb.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(ion_density_profile.data(), nxl, MPI_DOUBLE,
                global_zni.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/density_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# time[fs] = " << time / Const::femto << "\n";
        out << "# x[um]  n_bkg_e[m^-3]  Zni_profile[m^-3]  n_beam[m^-3]\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i < sg.nx_global; ++i) {
            out << (i + 0.5) * sg.dx / Const::micro
                << "  " << global_ne[i]
                << "  " << global_zni[i]
                << "  " << global_nb[i] << "\n";
        }
    }
}

void Diagnostics::write_electron_distribution(double time,
                                              const Species& electrons,
                                              const SpatialGrid& sg,
                                              int mpi_rank)
{
    std::ostringstream fname;
    fname << output_dir << "/fe_" << electrons.name << "_"
          << std::setw(5) << std::setfill('0') << snapshot_count
          << "_rank" << std::setw(4) << std::setfill('0') << mpi_rank
          << ".dat";

    std::ofstream out(fname.str().c_str());
    out << "# time[fs] = " << time / Const::femto << "\n";
    out << (electrons.cylindrical_mass_representation
            ? "# x[um]  u_parallel[p/(m c)]  u_perp[p/(m c)]  f_e[u^-3 m^-3]\n"
            : "# x[um]  u[p/(m c)]  mu  f_e[u^-3 m^-3]\n");
    out << std::scientific << std::setprecision(8);

    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int ix_g = ix + ng;
        const double x_um = sg.x(ix_g) / Const::micro;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u = electrons.vgrid.v(iv);
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double coordinate_1 = electrons.cylindrical_mass_representation
                    ? electrons.cgrid.upar_cells[iv] : u;
                const double coordinate_2 = electrons.cylindrical_mass_representation
                    ? electrons.cgrid.uperp_cells[imu] : electrons.vgrid.mu(imu);
                out << x_um << "  " << coordinate_1 << "  " << coordinate_2
                    << "  " << electrons.distribution_value(ix_g, iv, imu)
                    << "\n";
            }
        }
    }
}

void Diagnostics::advance_snapshot()
{
    ++snapshot_count;
}
