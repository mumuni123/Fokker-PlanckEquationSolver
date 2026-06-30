#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"
#include "beam_pic.h"
#include "vlasov_ampere_midpoint.h"
#include "config.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <vector>

double compute_dt(const Species& electron, const SpatialGrid& sg)
{
    double dt_min = 0.4 * sg.dx / Const::c;
    double e_est = Param::densb * Const::qe * sg.dx / Const::eps0;
    double udot_est = std::abs(electron.charge) * e_est /
                    (electron.mass * Const::c);
    if (udot_est > 1.0e-30) {
        double min_du = electron.vgrid.v_widths.empty()
                      ? electron.vgrid.dv
                      : electron.vgrid.v_widths[0];
        for (size_t iv = 1; iv < electron.vgrid.v_widths.size(); ++iv) {
            min_du = std::min(min_du, electron.vgrid.v_widths[iv]);
        }
        dt_min = std::min(dt_min, 0.25 * min_du / udot_est);
    }
    dt_min *= Param::dt_multiplier;
    return std::min(dt_min, 0.01 * Const::femto);
}

const char* field_solver_name()
{
    return "face-centered periodic Vlasov-Ampere update";
}

std::vector<double> build_local_ion_density_profile(const SpatialGrid& sg)
{
    return std::vector<double>(static_cast<size_t>(sg.nx_local), Param::dens);
}

void sync_moments_and_charge(Species& electrons,
                             const BeamPIC& beam,
                             EMFields& fields,
                             const std::vector<double>& ion_density_profile,
                             bool& moments_current)
{
    if (!moments_current) {
        electrons.compute_moments();
        moments_current = true;
    }
    fields.set_charge_density(electrons, beam.density, ion_density_profile);
}

struct BackgroundCurrentDiagnostics {
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

void reset_background_current_diagnostics(BackgroundCurrentDiagnostics& diag)
{
    diag.residual_if_charge = 0.0;
    diag.residual_if_ampere = 0.0;
    diag.e_dot_j_charge = 0.0;
    diag.e_dot_j_energy = 0.0;
    diag.e_dot_j_ampere = 0.0;
    diag.max_abs_j_charge = 0.0;
    diag.max_abs_j_energy = 0.0;
    diag.max_abs_j_ampere = 0.0;
    diag.max_abs_j_charge_minus_ampere = 0.0;
    diag.max_abs_j_energy_minus_ampere = 0.0;
}

void write_snapshot(Diagnostics& diag,
                    double time,
                    const Species& bkg_e,
                    const BeamPIC& beam,
                    EMFields& fields,
                    const std::vector<double>& ion_density_profile,
                    const SpatialGrid& sgrid,
                    int mpi_rank,
                    int mpi_size,
                    bool write_full_fe,
                    const std::vector<double>* bkg_energy_current_face = 0,
                    const std::vector<double>* bkg_ampere_current_face = 0)
{
    fields.compute_potential(mpi_rank, mpi_size);
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_current_density(time, bkg_e, beam, sgrid, mpi_rank, mpi_size,
                               bkg_energy_current_face,
                               bkg_ampere_current_face);
    diag.write_density_profile(time, bkg_e, beam.density, ion_density_profile,
                               sgrid, mpi_rank, mpi_size);
    diag.write_px_distribution(time, bkg_e, mpi_rank, mpi_size);
    if (write_full_fe) {
        diag.write_electron_distribution(time, bkg_e, sgrid, mpi_rank);
    }

    if (mpi_rank == 0) {
        printf("  >> Snapshot %d written at t = %.4f fs\n",
               diag.snapshot_count, time / Const::femto);
    }
    diag.advance_snapshot();
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    RuntimeConfig config = load_runtime_config();

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Background-electron VFP + fixed ions + PIC beam solver\n");
        printf("  Spherical electron momentum grid: (u, mu), u = p / (m c)\n");
        printf("============================================================\n");
        printf("MPI ranks: %d\n", mpi_size);
        #pragma omp parallel
        {
            #pragma omp master
            printf("OpenMP threads per rank: %d\n", omp_get_num_threads());
        }
        printf("Spatial grid: nx = %d, dx = %.3e m, Lx = %.3e m\n",
               Param::nx, Param::dx, Param::Lx);
        printf("Density profile: uniform plasma over full domain, n0 = %.3e /m^3\n",
               Param::dens);
        printf("Electron momentum grid: Nu x Nmu = %d x %d, nonuniform u with %d cells below u = %.3f\n",
               Param::Nv, Param::Nmu,
               Param::momentum_refined_cells,
               Param::momentum_refined_u);
        printf("Electron momentum domain: 0 <= u <= %.3f, vx = c u mu / sqrt(1+u^2)\n",
               Param::momentum_umax);
        printf("Spatial boundary: periodic in x for background electrons and electrostatic field; beam is open\n");
        printf("Electrostatic update: face-centered dE/dt = -J_total/eps0; zero mode evolves explicitly\n");
        printf("Field solver: %s\n", field_solver_name());
        printf("Fixed ions: uniform Z*n_i = %.3e /m^3\n", Param::dens);
        printf("Background electrons: periodic Vlasov transport, T_e = %.1f eV\n",
               Param::temperature_e / Const::eV);
        printf("PIC beam: gamma*beta = %.2f, beta = %.4f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::betab, Param::densb);
        printf("Beam source: sampled left-boundary crossings at x = 0\n");
        printf("Beam injection: remaining-time push with charge-conserving face current\n");
        printf("Beam charge compensation source: OFF; background density perturbations are produced only by Poisson/Vlasov dynamics\n");
        printf("Beam boundary: particles crossing either global edge leave the domain\n");
        printf("Background boundary: periodic ghost cells\n");
        printf("Beam macro weight: %.6e particles/m^2\n", Param::beam_macro_weight);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        printf("Debug diagnostics: %s\n",
               config.enable_debug_diagnostics ? "ON" : "OFF");
#else
        printf("Debug diagnostics: compile-time disabled\n");
#endif
        printf("Full fe distribution output: %s\n",
               config.enable_full_fe_output ? "ON" : "OFF");
        printf("Step diagnostics: %s\n",
               config.enable_step_diagnostics ? "ON" : "OFF");
        if (config.enable_step_diagnostics) {
            printf("Step diagnostics interval: every %d steps\n",
                   config.step_diagnostics_interval);
        }
        printf("Progress trace: %s\n",
               config.enable_progress_trace ? "ON" : "OFF");
        printf("------------------------------------------------------------\n");
    }

    SpatialGrid sgrid;
    sgrid.init(mpi_rank, mpi_size);
    std::vector<double> ion_density_profile = build_local_ion_density_profile(sgrid);

    Species bkg_e;
    bkg_e.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
               -Const::qe, Const::me,
               Param::dens, Param::temperature_e, false, sgrid);
    bkg_e.initialize_maxwellian();

    BeamPIC beam;
    beam.init(sgrid);
    beam.deposit_density(sgrid, mpi_rank, mpi_size);

    EMFields fields;
    fields.init(sgrid);

    VlasovAmpereMidpointSolver midpoint_solver;
    midpoint_solver.set_step_diagnostics_enabled(false);
    CollisionOperator collision;
    Diagnostics diag;
    diag.init("output", mpi_rank, config.enable_debug_diagnostics,
              config.enable_step_diagnostics);

    double dt = compute_dt(bkg_e, sgrid);
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    int nsteps = static_cast<int>(std::ceil(Param::t_end / dt));

    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    double cumulative_collision_energy_delta = 0.0;
    bool moments_current = false;
    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    fields.solve_poisson(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(0.0, 0, bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    std::vector<double> latest_bkg_energy_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    std::vector<double> latest_bkg_ampere_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    write_snapshot(diag, 0.0, bkg_e, beam, fields, ion_density_profile,
                   sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                   &latest_bkg_energy_current_face,
                   &latest_bkg_ampere_current_face);

    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 1000;
    int last_snapshot_step = 0;
    double cumulative_bkg_energy_residual = 0.0;
    std::ofstream bkg_energy_monitor;
    if (mpi_rank == 0) {
        bkg_energy_monitor.open("output/bkg_energy_monitor.dat");
        bkg_energy_monitor
            << "# step  time[fs]  x_limiter_active_fraction  "
            << "x_limiter_min_alpha  dKE_bkg_plus_W_bkg_E[J/m2]  "
            << "relative_residual_step  cumulative_relative_residual  "
            << "x_limiter_active_fraction_core  "
            << "x_limiter_active_fraction_boundary  "
            << "x_limiter_min_alpha_core  "
            << "x_limiter_min_alpha_boundary  "
            << "x_negative_mass_before_repair[m^-2]  "
            << "x_mass_added_by_positivity_repair[m^-2]  "
            << "max_abs_J_bkg_charge[A/m2]  "
            << "max_abs_J_bkg_energy_diagnostic[A/m2]  "
            << "max_abs_J_bkg_ampere[A/m2]  "
            << "max_abs_J_bkg_charge_minus_ampere[A/m2]  "
            << "max_abs_J_bkg_energy_minus_ampere[A/m2]  "
            << "E_dot_J_bkg_charge[W/m2]  "
            << "E_dot_J_bkg_energy_diagnostic[W/m2]  "
            << "E_dot_J_bkg_ampere[W/m2]  "
            << "residual_if_charge_current[J/m2]  "
            << "residual_if_ampere_current[J/m2]  "
            << "positivity_mass_defect[m^-2]  "
            << "positivity_energy_defect[J/m2]  "
            << "u_limiter_mass_delta[m^-2]  "
            << "u_limiter_px_delta[kg/m/s/m2]  "
            << "u_limiter_energy_delta[J/m2]  "
            << "u_force_alpha_min  u_force_alpha_active_frac  "
            << "coupled_iter  coupled_residual_E  "
            << "coupled_residual_J_bkg  coupled_residual_J_beam  "
            << "coupled_residual_bkg_mass  "
            << "coupled_residual_beam_continuity\n";
        bkg_energy_monitor << std::scientific;

        std::ofstream f_neg_monitor;
        f_neg_monitor.open("output/f_negativity_monitor.dat");
        f_neg_monitor
            << "# step  time[fs]  min_f  neg_ratio_max  "
            << "neg_mass_total[m^-2]  neg_cell_count  "
            << "x_worst  u_worst  mu_worst\n";
        f_neg_monitor << std::scientific << std::setprecision(8);
        f_neg_monitor.close();
    }
    for (int step = 1; step <= nsteps; ++step) {
        double time = step * dt;
        int nsub_v1 = 0;
        int nsub_mu1 = 0;
        int nsub_v2 = 0;
        int nsub_mu2 = 0;
        double loss_v1 = 0.0;
        double loss_v1_low = 0.0;
        double loss_v1_high = 0.0;
        double loss_mu1 = 0.0;
        double loss_v2 = 0.0;
        double loss_v2_low = 0.0;
        double loss_v2_high = 0.0;
        double loss_mu2 = 0.0;
        double net_nb_change_step = 0.0;
        double collision_energy_step = 0.0;
        const bool collect_step_diagnostics =
            should_write_step_diagnostics(config, step);
        double dke_bkg_step = 0.0;
        double dke_beam_push = 0.0;
        double dE_field_step = 0.0;
        double W_bkg_E = 0.0;
        double W_beam_E = 0.0;
        double v_mass_error_step = 0.0;
        double mu_mass_error_step = 0.0;
        double v_momentum_delta_step = 0.0;
        double mu_momentum_delta_step = 0.0;
        double v_energy_delta_step = 0.0;
        double mu_energy_delta_step = 0.0;
        double E_src_in_step = 0.0;
        double E_src_out_step = 0.0;
        double E_balance_step = 0.0;
        double max_loss_u_high_step = 0.0;
        double x_at_max_loss_u_high_step = 0.0;
        double f_u_max_x_step = 0.0;
        double integral_f_u_gt_8_x_step = 0.0;
        double bkg_number_step_start = 0.0;
        double bkg_ke_step_start = 0.0;
        double global_bkg_ke_step_start = 0.0;
        double beam_ke_step_start = 0.0;
        double field_energy_step_start = 0.0;
        double x_limiter_active_fraction_step = 0.0;
        double x_limiter_min_alpha_step = 1.0;
        double x_limiter_active_fraction_core_step = 0.0;
        double x_limiter_active_fraction_boundary_step = 0.0;
        double x_limiter_min_alpha_core_step = 1.0;
        double x_limiter_min_alpha_boundary_step = 1.0;
        double x_negative_mass_before_repair_step = 0.0;
        double x_mass_added_by_positivity_repair_step = 0.0;
        double positivity_energy_defect_step = 0.0;
        double positivity_mass_defect_step = 0.0;
        double u_limiter_mass_delta_step = 0.0;
        double u_limiter_momentum_delta_step = 0.0;
        double u_limiter_energy_delta_step = 0.0;
        double u_force_alpha_min_step = 1.0;
        double u_force_alpha_active_frac_step = 0.0;
        double f_neg_min_step = 0.0;
        double f_neg_ratio_max_step = 0.0;
        double f_neg_mass_total_step = 0.0;
        long long f_neg_cell_count_step = 0;
        int f_neg_ix_step = -1;
        int f_neg_iv_step = -1;
        int f_neg_imu_step = -1;
        double local_bkg_energy_residual_step = 0.0;
        double bkg_energy_residual_step = 0.0;
        double bkg_energy_relative_residual_step = 0.0;
        int coupled_iter_step = 0;
        double coupled_residual_E_step = 0.0;
        double coupled_residual_J_bkg_step = 0.0;
        double coupled_residual_J_beam_step = 0.0;
        double coupled_residual_bkg_mass_step = 0.0;
        double coupled_residual_beam_continuity_step = 0.0;
        BackgroundCurrentDiagnostics local_bkg_current_diag_step;
        reset_background_current_diagnostics(local_bkg_current_diag_step);
        BackgroundCurrentDiagnostics global_bkg_current_diag_step;
        reset_background_current_diagnostics(global_bkg_current_diag_step);
        const Species bkg_step_start = bkg_e;
        const BeamPIC beam_step_start = beam;
        const EMFields fields_step_start = fields;
        bkg_e.total_particle_number_and_energy(bkg_number_step_start,
                                               bkg_ke_step_start);
        double global_bkg_start_values[2] = {
            bkg_number_step_start,
            bkg_ke_step_start
        };
        MPI_Allreduce(MPI_IN_PLACE, global_bkg_start_values, 2,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_bkg_ke_step_start = global_bkg_start_values[1];
        if (collect_step_diagnostics) {
            beam_ke_step_start = beam.total_kinetic_energy();
            field_energy_step_start = fields.total_energy();
        }

        trace_progress(config, mpi_rank, step,
                       "before coupled midpoint FV solve");
        midpoint_solver.set_step_diagnostics_enabled(collect_step_diagnostics);
        VlasovAmpereMidpointSolver::Result midpoint_result =
            midpoint_solver.advance_background_and_fields(
                bkg_step_start, beam_step_start, fields_step_start, sgrid,
                dt, time, mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step,
                       "after coupled midpoint FV solve");

        if (!midpoint_result.converged || midpoint_result.failed) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "WARNING: coupled midpoint FV solve failed at "
                             "step %d, t = %.6e s; residual %.6e, "
                             "field %.6e, f %.6e, J_bkg %.6e, J_beam %.6e, "
                             "bkg_mass %.6e, beam_continuity %.6e after "
                             "%d iterations, substeps %d; limiter active "
                             "%.6e, min_alpha %.6e, core_active %.6e, "
                             "boundary_active %.6e, core_min_alpha %.6e, "
                             "boundary_min_alpha %.6e. "
                             "Continuing with best available result.\n",
                             step, time,
                             midpoint_result.nonlinear_residual,
                             midpoint_result.residual_E,
                             midpoint_result.residual_f,
                             midpoint_result.residual_J_bkg,
                             midpoint_result.residual_J_beam,
                             midpoint_result.continuity_residual_bkg,
                             midpoint_result.beam_continuity_residual,
                             midpoint_result.nonlinear_iterations,
                             midpoint_result.substeps_used,
                             midpoint_result.limiter_active_fraction,
                             midpoint_result.limiter_min_alpha,
                             midpoint_result.limiter_active_fraction_core,
                             midpoint_result.limiter_active_fraction_boundary,
                             midpoint_result.limiter_min_alpha_core,
                             midpoint_result.limiter_min_alpha_boundary);
            }
            // Protection disabled for long-run test — continue with best result
            // MPI_Abort(MPI_COMM_WORLD, 8);
        }

        bkg_e = midpoint_result.species_np1;
        beam = midpoint_result.beam_np1;
        fields = midpoint_result.fields_np1;
        latest_bkg_energy_current_face =
            midpoint_result.j_bkg_energy_debug_face;
        moments_current = true;
        latest_bkg_ampere_current_face = midpoint_result.j_bkg_face_mid;
        local_bkg_current_diag_step.residual_if_charge =
            midpoint_result.current_diag.residual_if_charge;
        local_bkg_current_diag_step.residual_if_ampere =
            midpoint_result.current_diag.residual_if_ampere;
        local_bkg_current_diag_step.e_dot_j_charge =
            midpoint_result.current_diag.e_dot_j_charge;
        local_bkg_current_diag_step.e_dot_j_energy =
            midpoint_result.current_diag.e_dot_j_energy;
        local_bkg_current_diag_step.e_dot_j_ampere =
            midpoint_result.current_diag.e_dot_j_ampere;
        local_bkg_current_diag_step.max_abs_j_charge =
            midpoint_result.current_diag.max_abs_j_charge;
        local_bkg_current_diag_step.max_abs_j_energy =
            midpoint_result.current_diag.max_abs_j_energy;
        local_bkg_current_diag_step.max_abs_j_ampere =
            midpoint_result.current_diag.max_abs_j_ampere;
        local_bkg_current_diag_step.max_abs_j_charge_minus_ampere =
            midpoint_result.current_diag.max_abs_j_charge_minus_ampere;
        local_bkg_current_diag_step.max_abs_j_energy_minus_ampere =
            midpoint_result.current_diag.max_abs_j_energy_minus_ampere;
        W_bkg_E = midpoint_result.field_work_bkg;
        W_beam_E = midpoint_result.field_work_beam;
        dke_beam_push = midpoint_result.delta_ke_beam;
        x_limiter_active_fraction_step =
            midpoint_result.limiter_active_fraction;
        x_limiter_min_alpha_step = midpoint_result.limiter_min_alpha;
        x_limiter_active_fraction_core_step =
            midpoint_result.limiter_active_fraction_core;
        x_limiter_active_fraction_boundary_step =
            midpoint_result.limiter_active_fraction_boundary;
        x_limiter_min_alpha_core_step =
            midpoint_result.limiter_min_alpha_core;
        x_limiter_min_alpha_boundary_step =
            midpoint_result.limiter_min_alpha_boundary;
        x_negative_mass_before_repair_step =
            midpoint_result.x_negative_mass_before_repair;
        x_mass_added_by_positivity_repair_step =
            midpoint_result.x_mass_added_by_positivity_repair;
        positivity_energy_defect_step =
            midpoint_result.positivity_energy_defect;
        positivity_mass_defect_step =
            midpoint_result.positivity_mass_defect;
        u_limiter_mass_delta_step =
            midpoint_result.limiter_mass_defect;
        u_limiter_momentum_delta_step =
            midpoint_result.limiter_momentum_defect;
        u_limiter_energy_delta_step =
            midpoint_result.limiter_energy_defect;
        u_force_alpha_min_step =
            midpoint_result.u_force_alpha_min;
        u_force_alpha_active_frac_step =
            midpoint_result.u_force_alpha_active_frac;
        f_neg_min_step = midpoint_result.f_neg_min;
        f_neg_ratio_max_step = midpoint_result.f_neg_ratio_max;
        f_neg_mass_total_step = midpoint_result.f_neg_mass_total;
        f_neg_cell_count_step = midpoint_result.f_neg_cell_count;
        f_neg_ix_step = midpoint_result.f_neg_ix;
        f_neg_iv_step = midpoint_result.f_neg_iv;
        f_neg_imu_step = midpoint_result.f_neg_imu;
        coupled_iter_step = midpoint_result.nonlinear_iterations;
        coupled_residual_E_step = midpoint_result.residual_E;
        coupled_residual_J_bkg_step = midpoint_result.residual_J_bkg;
        coupled_residual_J_beam_step = midpoint_result.residual_J_beam;
        coupled_residual_bkg_mass_step =
            midpoint_result.continuity_residual_bkg;
        coupled_residual_beam_continuity_step =
            midpoint_result.beam_continuity_residual;

        const double bkg_ke_step_end_for_residual =
            bkg_ke_step_start + midpoint_result.delta_ke_bkg;
        local_bkg_energy_residual_step =
            (bkg_ke_step_end_for_residual - bkg_ke_step_start) + W_bkg_E;
        const double local_bkg_energy_values[8] = {
            bkg_ke_step_end_for_residual,
            W_bkg_E,
            local_bkg_energy_residual_step,
            local_bkg_current_diag_step.residual_if_charge,
            local_bkg_current_diag_step.residual_if_ampere,
            local_bkg_current_diag_step.e_dot_j_charge,
            local_bkg_current_diag_step.e_dot_j_energy,
            local_bkg_current_diag_step.e_dot_j_ampere
        };
        double global_bkg_energy_values[8] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        MPI_Allreduce(local_bkg_energy_values, global_bkg_energy_values, 8,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double local_current_max_values[5] = {
            local_bkg_current_diag_step.max_abs_j_charge,
            local_bkg_current_diag_step.max_abs_j_energy,
            local_bkg_current_diag_step.max_abs_j_ampere,
            local_bkg_current_diag_step.max_abs_j_charge_minus_ampere,
            local_bkg_current_diag_step.max_abs_j_energy_minus_ampere
        };
        double global_current_max_values[5] = {
            0.0, 0.0, 0.0, 0.0, 0.0
        };
        MPI_Allreduce(local_current_max_values, global_current_max_values, 5,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        global_bkg_current_diag_step.residual_if_charge =
            global_bkg_energy_values[3];
        global_bkg_current_diag_step.residual_if_ampere =
            global_bkg_energy_values[4];
        global_bkg_current_diag_step.e_dot_j_charge =
            global_bkg_energy_values[5];
        global_bkg_current_diag_step.e_dot_j_energy =
            global_bkg_energy_values[6];
        global_bkg_current_diag_step.e_dot_j_ampere =
            global_bkg_energy_values[7];
        global_bkg_current_diag_step.max_abs_j_charge =
            global_current_max_values[0];
        global_bkg_current_diag_step.max_abs_j_energy =
            global_current_max_values[1];
        global_bkg_current_diag_step.max_abs_j_ampere =
            global_current_max_values[2];
        global_bkg_current_diag_step.max_abs_j_charge_minus_ampere =
            global_current_max_values[3];
        global_bkg_current_diag_step.max_abs_j_energy_minus_ampere =
            global_current_max_values[4];
        const double global_dke_bkg_step =
            global_bkg_energy_values[0] - global_bkg_ke_step_start;
        const double global_W_bkg_E = global_bkg_energy_values[1];
        bkg_energy_residual_step = global_bkg_energy_values[2];
        const double bkg_energy_residual_den =
            std::max(std::max(std::fabs(global_dke_bkg_step),
                              std::fabs(global_W_bkg_E)),
                     std::max(1.0,
                              1.0e-12 *
                              std::fabs(global_bkg_ke_step_start)));
        bkg_energy_relative_residual_step =
            std::fabs(bkg_energy_residual_step) /
            bkg_energy_residual_den;
        cumulative_bkg_energy_residual +=
            bkg_energy_relative_residual_step;

        if (mpi_rank == 0 && step % 100 == 0) {
            bkg_energy_monitor << step << "  "
                               << time / Const::femto << "  "
                               << x_limiter_active_fraction_step << "  "
                               << x_limiter_min_alpha_step << "  "
                               << bkg_energy_residual_step << "  "
                               << bkg_energy_relative_residual_step << "  "
                               << cumulative_bkg_energy_residual << "  "
                               << x_limiter_active_fraction_core_step << "  "
                               << x_limiter_active_fraction_boundary_step << "  "
                               << x_limiter_min_alpha_core_step << "  "
                               << x_limiter_min_alpha_boundary_step << "  "
                               << x_negative_mass_before_repair_step << "  "
                               << x_mass_added_by_positivity_repair_step << "  "
                               << global_bkg_current_diag_step.max_abs_j_charge << "  "
                               << global_bkg_current_diag_step.max_abs_j_energy << "  "
                               << global_bkg_current_diag_step.max_abs_j_ampere << "  "
                               << global_bkg_current_diag_step.max_abs_j_charge_minus_ampere << "  "
                               << global_bkg_current_diag_step.max_abs_j_energy_minus_ampere << "  "
                               << global_bkg_current_diag_step.e_dot_j_charge << "  "
                               << global_bkg_current_diag_step.e_dot_j_energy << "  "
                               << global_bkg_current_diag_step.e_dot_j_ampere << "  "
                               << global_bkg_current_diag_step.residual_if_charge << "  "
                               << global_bkg_current_diag_step.residual_if_ampere << "  "
                               << positivity_mass_defect_step << "  "
                               << positivity_energy_defect_step << "  "
                               << u_limiter_mass_delta_step << "  "
                               << u_limiter_momentum_delta_step << "  "
                               << u_limiter_energy_delta_step << "  "
                               << u_force_alpha_min_step << "  "
                               << u_force_alpha_active_frac_step << "  "
                               << coupled_iter_step << "  "
                               << coupled_residual_E_step << "  "
                               << coupled_residual_J_bkg_step << "  "
                               << coupled_residual_J_beam_step << "  "
                               << coupled_residual_bkg_mass_step << "  "
                               << coupled_residual_beam_continuity_step << "\n";
            bkg_energy_monitor.flush();

            // f-negativity monitor: append every 100 steps
            std::ofstream f_neg_monitor;
            f_neg_monitor.open("output/f_negativity_monitor.dat",
                               std::ios::app);
            f_neg_monitor << step << "  "
                          << time / Const::femto << "  "
                          << f_neg_min_step << "  "
                          << f_neg_ratio_max_step << "  "
                          << f_neg_mass_total_step << "  "
                          << f_neg_cell_count_step << "  "
                          << f_neg_ix_step << "  "
                          << f_neg_iv_step << "  "
                          << f_neg_imu_step << "\n";
            f_neg_monitor.close();
        }
        if (bkg_e.collisions_enabled) {
            trace_progress(config, mpi_rank, step, "before collisions");
            collision_energy_step +=
                collision.apply(bkg_e, dt, Param::dens, Param::temperature_e,
                                Const::me, 1.0, 1.0);
            collision_energy_step +=
                collision.apply(bkg_e, dt, Param::dens / Param::Z_ion,
                                Param::temperature_i, Param::mass_ion,
                                (double)Param::Z_ion, 1.0);
            cumulative_collision_energy_delta += collision_energy_step;
            trace_progress(config, mpi_rank, step, "after collisions");
            moments_current = false;
        }

        trace_progress(config, mpi_rank, step, "before end sync");
        if (!moments_current) {
            bkg_e.compute_moments();
            moments_current = true;
        }
        fields.set_charge_density(bkg_e, beam.density, ion_density_profile);
        fields.update_gauss_residual_diagnostics(mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after end sync");

        net_nb_change_step = beam.last_injected_number()
                           - beam.last_outflow_number();

        if (collect_step_diagnostics) {
            const double bkg_ke_step_end = bkg_e.total_kinetic_energy();
            const double beam_ke_step_end = beam.total_kinetic_energy();
            const double field_energy_step_end = fields.total_energy();
            dke_bkg_step = bkg_ke_step_end - bkg_ke_step_start;
            dE_field_step = field_energy_step_end - field_energy_step_start;
            E_src_in_step = beam.last_injected_energy();
            E_src_out_step = beam.last_outflow_energy();
            const double total_energy_delta =
                (bkg_ke_step_end + beam_ke_step_end + field_energy_step_end) -
                (bkg_ke_step_start + beam_ke_step_start + field_energy_step_start);
            E_balance_step =
                total_energy_delta - E_src_in_step + E_src_out_step
                - collision_energy_step;

            diag.write_step_diagnostics(step, time, bkg_e, beam, fields,
                                        sgrid, mpi_rank, mpi_size,
                                        nsub_v1, nsub_mu1,
                                        nsub_v2, nsub_mu2,
                                        loss_v1, loss_mu1,
                                        loss_v2, loss_mu2,
                                        loss_v1_low, loss_v1_high,
                                        loss_v2_low, loss_v2_high,
                                        net_nb_change_step,
                                        collision_energy_step,
                                        cumulative_collision_energy_delta,
                                        dke_bkg_step, dke_beam_push,
                                        dE_field_step, W_bkg_E, W_beam_E,
                                        v_mass_error_step,
                                        mu_mass_error_step,
                                        v_momentum_delta_step,
                                        mu_momentum_delta_step,
                                        v_energy_delta_step,
                                        mu_energy_delta_step,
                                        E_src_in_step,
                                        E_src_out_step,
                                        E_balance_step,
                                        x_limiter_active_fraction_step,
                                        x_limiter_min_alpha_step,
                                        local_bkg_energy_residual_step,
                                        local_bkg_current_diag_step.max_abs_j_charge,
                                        local_bkg_current_diag_step.max_abs_j_energy,
                                        local_bkg_current_diag_step.max_abs_j_ampere,
                                        local_bkg_current_diag_step.max_abs_j_charge_minus_ampere,
                                        local_bkg_current_diag_step.max_abs_j_energy_minus_ampere,
                                        local_bkg_current_diag_step.e_dot_j_charge,
                                        local_bkg_current_diag_step.e_dot_j_energy,
                                        local_bkg_current_diag_step.e_dot_j_ampere,
                                        local_bkg_current_diag_step.residual_if_charge,
                                        local_bkg_current_diag_step.residual_if_ampere,
                                        coupled_iter_step,
                                        coupled_residual_E_step,
                                        coupled_residual_J_bkg_step,
                                        coupled_residual_J_beam_step,
                                        max_loss_u_high_step,
                                        x_at_max_loss_u_high_step,
                                        f_u_max_x_step,
                                        integral_f_u_gt_8_x_step);
        }

        if (step % stdout_freq == 0) {
            diag.write_scalars(time, step, bkg_e, beam, fields,
                               cumulative_collision_energy_delta,
                               mpi_rank, mpi_size);
            if (mpi_rank == 0) {
                printf("Step %d / %d, t = %.4f fs\n", step, nsteps, time / Const::femto);
            }
        }

        if (time >= next_snapshot) {
            write_snapshot(diag, time, bkg_e, beam, fields, ion_density_profile,
                           sgrid, mpi_rank, mpi_size,
                           config.enable_full_fe_output,
                           &latest_bkg_energy_current_face,
                           &latest_bkg_ampere_current_face);
            last_snapshot_step = step;
            next_snapshot += Param::dt_snapshot;
        }
    }

    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    fields.update_gauss_residual_diagnostics(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(nsteps, Param::t_end, "final", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(Param::t_end, nsteps, bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    if (last_snapshot_step != nsteps) {
        write_snapshot(diag, Param::t_end, bkg_e, beam, fields, ion_density_profile,
                       sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                       &latest_bkg_energy_current_face,
                       &latest_bkg_ampere_current_face);
    }

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Simulation complete: t = %.1f fs, %d steps\n",
               Param::t_end / Const::femto, nsteps);
        printf("============================================================\n");
    }

    MPI_Finalize();
    return 0;
}
