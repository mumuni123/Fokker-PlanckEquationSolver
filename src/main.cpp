#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "vlasov.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"
#include "beam_pic.h"
#include "config.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

void abort_if_vmax_loss(const VlasovSolver& vlasov,
                        int step,
                        double time,
                        const char* stage,
                        int mpi_rank)
{
    if (!Param::abort_on_vmax_loss) return;

    const double local_loss = vlasov.last_loss_v_high();
    double global_loss = 0.0;
    MPI_Allreduce(&local_loss, &global_loss, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    const double threshold =
        Param::umax_loss_abort_fraction * Param::dens * Param::plasma_length;
    if (!(global_loss <= threshold)) {
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "ERROR: background electron distribution reached umax "
                         "during %s at step %d, t = %.6e s. "
                         "loss_u_high = %.8e, threshold = %.8e. "
                         "Stopping to avoid amplifying nonphysical results.\n",
                         stage, step, time, global_loss, threshold);
        }
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
}

std::vector<double> copy_local_ex_face(const EMFields& fields,
                                       const SpatialGrid& sg)
{
    std::vector<double> ex_face(static_cast<size_t>(sg.nx_local), 0.0);
    const size_t n = std::min(ex_face.size(), fields.Ex_face.size());
    for (size_t iface = 0; iface < n; ++iface) {
        ex_face[iface] = fields.Ex_face[iface];
    }
    return ex_face;
}

double integrate_ampere_face_work(const std::vector<double>& current_face,
                                  const std::vector<double>& ex_face_before,
                                  const EMFields& fields,
                                  const SpatialGrid& sg,
                                  double dt)
{
    double work = 0.0;
    const size_t unique_faces = static_cast<size_t>(sg.nx_local);
    const size_t n = std::min(unique_faces,
                              std::min(current_face.size(),
                                       std::min(ex_face_before.size(),
                                                fields.Ex_face.size())));
    for (size_t iface = 0; iface < n; ++iface) {
        const double ex_mid =
            0.5 * (ex_face_before[iface] + fields.Ex_face[iface]);
        work -= current_face[iface] * ex_mid * sg.dx * dt;
    }
    return work;
}

double integrate_static_face_work(const std::vector<double>& current_face,
                                  const EMFields& fields,
                                  const SpatialGrid& sg,
                                  double dt)
{
    double work = 0.0;
    const size_t unique_faces = static_cast<size_t>(sg.nx_local);
    const size_t n =
        std::min(unique_faces, std::min(current_face.size(),
                                        fields.Ex_face.size()));
    for (size_t iface = 0; iface < n; ++iface) {
        work += current_face[iface] * fields.Ex_face[iface] * sg.dx * dt;
    }
    return work;
}

void set_midpoint_face_field(EMFields& fields,
                             const std::vector<double>& ex_old,
                             const std::vector<double>& ex_mid,
                             const SpatialGrid& sg,
                             int mpi_rank,
                             int mpi_size)
{
    const size_t n =
        std::min(static_cast<size_t>(sg.nx_local), fields.Ex_face.size());
    for (size_t iface = 0; iface < n; ++iface) {
        const double old_value = (iface < ex_old.size()) ? ex_old[iface] : 0.0;
        const double mid_value = (iface < ex_mid.size()) ? ex_mid[iface] : old_value;
        fields.Ex_face[iface] = mid_value;
    }
    if (fields.Ex_face.size() > static_cast<size_t>(sg.nx_local)) {
        fields.Ex_face[static_cast<size_t>(sg.nx_local)] =
            (sg.nx_local == 0) ? 0.0 : fields.Ex_face[0];
    }
    fields.sync_cell_ex_from_faces(mpi_rank, mpi_size);
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
                    bool write_full_fe)
{
    fields.compute_potential(mpi_rank, mpi_size);
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_current_density(time, bkg_e, beam, sgrid, mpi_rank, mpi_size);
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
        printf("Electrostatic update: dE/dt = -(J_total - <J_total>)/eps0 with <E> = 0\n");
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

    VlasovSolver vlasov;
    vlasov.set_step_diagnostics_enabled(false);
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
    write_snapshot(diag, 0.0, bkg_e, beam, fields, ion_density_profile,
                   sgrid, mpi_rank, mpi_size, config.enable_full_fe_output);

    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 1000;
    int last_snapshot_step = 0;
    for (int step = 1; step <= nsteps; ++step) {
        double time = step * dt;
        const double time_start = time - dt;
        const double time_center = time_start + 0.5 * dt;
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
        double beam_push_ke_before = 0.0;
        double bkg_ke_step_start = 0.0;
        double beam_ke_step_start = 0.0;
        double field_energy_step_start = 0.0;
        std::vector<double> ex_face_before_ampere;
        vlasov.set_step_diagnostics_enabled(collect_step_diagnostics);

        const Species bkg_step_start = bkg_e;
        const BeamPIC beam_step_start = beam;
        const EMFields fields_step_start = fields;
        const std::vector<double> ex_face_step_start =
            copy_local_ex_face(fields_step_start, sgrid);
        if (collect_step_diagnostics) {
            bkg_ke_step_start = bkg_e.total_kinetic_energy();
            beam_ke_step_start = beam.total_kinetic_energy();
            field_energy_step_start = fields.total_energy();
        }

        std::vector<double> ex_mid_trial = ex_face_step_start;
        std::vector<double> ex_mid_next = ex_face_step_start;
        Species accepted_bkg = bkg_step_start;
        BeamPIC accepted_beam = beam_step_start;
        EMFields accepted_fields = fields_step_start;
        const int max_coupled_iters = 30;
        const double field_tol = 1.0e-6;
        const double field_floor = 1.0;
        const double current_tol = 1.0e-5;
        double midpoint_omega = 0.5;
        double previous_normalized_error = -1.0;
        int residual_decrease_count = 0;
        std::vector<double> current_trial_prev(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        std::vector<double> current_trial_next(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        std::vector<double> current_bkg_trial_prev(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        std::vector<double> current_bkg_trial_next(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        std::vector<double> current_beam_trial_prev(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        std::vector<double> current_beam_trial_next(
            static_cast<size_t>(sgrid.nx_local), 0.0);
        bool have_previous_current_trial = false;
        bool coupled_converged = false;
        double final_coupled_error = 0.0;
        double final_field_error = 0.0;
        double final_current_error = 0.0;
        double final_max_delta_ex = 0.0;
        double final_e_abs_tol = 0.0;
        double final_j_abs_tol = 0.0;

        for (int coupled_iter = 0; coupled_iter < max_coupled_iters;
             ++coupled_iter) {
            bkg_e = bkg_step_start;
            beam = beam_step_start;
            fields = fields_step_start;
            set_midpoint_face_field(fields, ex_face_step_start, ex_mid_trial,
                                    sgrid, mpi_rank, mpi_size);

            nsub_v1 = nsub_mu1 = nsub_v2 = nsub_mu2 = 0;
            loss_v1 = loss_v1_low = loss_v1_high = 0.0;
            loss_mu1 = loss_v2 = loss_v2_low = loss_v2_high = loss_mu2 = 0.0;
            dke_beam_push = 0.0;
            W_bkg_E = 0.0;
            W_beam_E = 0.0;
            v_mass_error_step = 0.0;
            mu_mass_error_step = 0.0;
            v_momentum_delta_step = 0.0;
            mu_momentum_delta_step = 0.0;
            v_energy_delta_step = 0.0;
            mu_energy_delta_step = 0.0;
            max_loss_u_high_step = 0.0;
            x_at_max_loss_u_high_step = 0.0;
            f_u_max_x_step = 0.0;
            integral_f_u_gt_8_x_step = 0.0;

            beam.begin_step(sgrid, dt);

            trace_progress(config, mpi_rank, step, "midpoint before v half 1");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt,
                            mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "midpoint after v half 1");
            nsub_v1 = vlasov.last_nsub_v();
            loss_v1 = vlasov.last_loss_v();
            loss_v1_low = vlasov.last_loss_v_low();
            loss_v1_high = vlasov.last_loss_v_high();
            if (vlasov.last_loss_v_high_local_max() > max_loss_u_high_step) {
                max_loss_u_high_step = vlasov.last_loss_v_high_local_max();
                x_at_max_loss_u_high_step = vlasov.last_x_at_max_loss_v_high();
                f_u_max_x_step = vlasov.last_f_umax_at_max_loss_v_high();
                integral_f_u_gt_8_x_step =
                    vlasov.last_integral_f_u_gt_8_at_max_loss_v_high();
            }
            v_mass_error_step += vlasov.last_mass_error_v();
            v_momentum_delta_step += vlasov.last_momentum_delta_v();
            if (collect_step_diagnostics) {
                v_energy_delta_step +=
                    integrate_static_face_work(vlasov.last_energy_current_face_x(),
                                               fields, sgrid, 0.5 * dt);
            }
            abort_if_vmax_loss(vlasov, step, time, "midpoint_v_half_1",
                               mpi_rank);

            trace_progress(config, mpi_rank, step, "midpoint before mu half 1");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "midpoint after mu half 1");
            nsub_mu1 = vlasov.last_nsub_mu();
            loss_mu1 = vlasov.last_loss_mu();
            mu_mass_error_step += vlasov.last_mass_error_mu();
            mu_momentum_delta_step += vlasov.last_momentum_delta_mu();
            mu_energy_delta_step += vlasov.last_energy_delta_mu();

            trace_progress(config, mpi_rank, step, "midpoint before x full");
            vlasov.advect_x(bkg_e, sgrid, dt, mpi_rank, mpi_size, time_center);
            trace_progress(config, mpi_rank, step, "midpoint after x full");

            trace_progress(config, mpi_rank, step, "midpoint before mu half 2");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "midpoint after mu half 2");
            nsub_mu2 = vlasov.last_nsub_mu();
            loss_mu2 = vlasov.last_loss_mu();
            mu_mass_error_step += vlasov.last_mass_error_mu();
            mu_momentum_delta_step += vlasov.last_momentum_delta_mu();
            mu_energy_delta_step += vlasov.last_energy_delta_mu();

            trace_progress(config, mpi_rank, step, "midpoint before v half 2");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt,
                            mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "midpoint after v half 2");
            nsub_v2 = vlasov.last_nsub_v();
            loss_v2 = vlasov.last_loss_v();
            loss_v2_low = vlasov.last_loss_v_low();
            loss_v2_high = vlasov.last_loss_v_high();
            if (vlasov.last_loss_v_high_local_max() > max_loss_u_high_step) {
                max_loss_u_high_step = vlasov.last_loss_v_high_local_max();
                x_at_max_loss_u_high_step = vlasov.last_x_at_max_loss_v_high();
                f_u_max_x_step = vlasov.last_f_umax_at_max_loss_v_high();
                integral_f_u_gt_8_x_step =
                    vlasov.last_integral_f_u_gt_8_at_max_loss_v_high();
            }
            v_mass_error_step += vlasov.last_mass_error_v();
            v_momentum_delta_step += vlasov.last_momentum_delta_v();
            if (collect_step_diagnostics) {
                v_energy_delta_step +=
                    integrate_static_face_work(vlasov.last_energy_current_face_x(),
                                               fields, sgrid, 0.5 * dt);
            }
            abort_if_vmax_loss(vlasov, step, time, "midpoint_v_half_2",
                               mpi_rank);

            trace_progress(config, mpi_rank, step, "midpoint before beam");
            beam.inject(sgrid, fields, dt, time, mpi_rank, mpi_size);
            if (collect_step_diagnostics) {
                beam_push_ke_before = beam.total_kinetic_energy();
            }
            beam.push(sgrid, fields, dt, mpi_rank, mpi_size);
            if (collect_step_diagnostics) {
                dke_beam_push += beam.total_kinetic_energy() - beam_push_ke_before;
            }
            beam.deposit_density(sgrid, mpi_rank, mpi_size);
            beam.finalize_charge_conserving_current(sgrid, dt,
                                                    mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "midpoint after beam");
            W_beam_E = collect_step_diagnostics ? beam.last_field_work() : 0.0;

            bkg_e.compute_moments();
            fields = fields_step_start;
            if (collect_step_diagnostics) {
                ex_face_before_ampere = copy_local_ex_face(fields, sgrid);
            }
            fields.advance_ampere_face(bkg_e.current_face_x,
                                       beam.current_face_x,
                                       dt, mpi_rank, mpi_size);
            if (collect_step_diagnostics) {
                W_bkg_E =
                    integrate_ampere_face_work(bkg_e.current_face_x,
                                               ex_face_before_ampere,
                                               fields, sgrid, dt);
            }

            for (int iface = 0; iface < sgrid.nx_local; ++iface) {
                const size_t slot = static_cast<size_t>(iface);
                const double jbkg =
                    (slot < bkg_e.current_face_x.size())
                    ? bkg_e.current_face_x[slot] : 0.0;
                const double jbeam =
                    (slot < beam.current_face_x.size())
                    ? beam.current_face_x[slot] : 0.0;
                current_bkg_trial_next[slot] = jbkg;
                current_beam_trial_next[slot] = jbeam;
                current_trial_next[slot] = jbkg + jbeam;
            }

            const std::vector<double> ex_face_new =
                copy_local_ex_face(fields, sgrid);
            ex_mid_next = ex_face_new;
            const size_t nface =
                std::min(ex_mid_next.size(), ex_face_step_start.size());
            for (size_t iface = 0; iface < nface; ++iface) {
                ex_mid_next[iface] =
                    0.5 * (ex_face_step_start[iface] + ex_face_new[iface]);
            }

            double local_max_delta_ex = 0.0;
            double local_max_abs_ex = 0.0;
            const size_t field_faces =
                std::min(static_cast<size_t>(sgrid.nx_local),
                         std::min(ex_mid_next.size(), ex_mid_trial.size()));
            for (size_t iface = 0; iface < field_faces; ++iface) {
                local_max_delta_ex =
                    std::max(local_max_delta_ex,
                             std::fabs(ex_mid_next[iface]
                                     - ex_mid_trial[iface]));
                local_max_abs_ex =
                    std::max(local_max_abs_ex,
                             std::max(std::fabs(ex_mid_next[iface]),
                                      std::fabs(ex_mid_trial[iface])));
            }
            double local_max_delta_j_total = 0.0;
            double local_max_delta_j_bkg = 0.0;
            double local_max_delta_j_beam = 0.0;
            int local_max_delta_j_total_face = 0;
            int local_max_delta_j_bkg_face = 0;
            int local_max_delta_j_beam_face = 0;
            if (have_previous_current_trial) {
                for (size_t iface = 0; iface < current_trial_next.size();
                     ++iface) {
                    const double delta_total =
                        std::fabs(current_trial_next[iface]
                                - current_trial_prev[iface]);
                    const double delta_bkg =
                        std::fabs(current_bkg_trial_next[iface]
                                - current_bkg_trial_prev[iface]);
                    const double delta_beam =
                        std::fabs(current_beam_trial_next[iface]
                                - current_beam_trial_prev[iface]);
                    if (delta_total > local_max_delta_j_total) {
                        local_max_delta_j_total = delta_total;
                        local_max_delta_j_total_face =
                            static_cast<int>(iface);
                    }
                    if (delta_bkg > local_max_delta_j_bkg) {
                        local_max_delta_j_bkg = delta_bkg;
                        local_max_delta_j_bkg_face =
                            static_cast<int>(iface);
                    }
                    if (delta_beam > local_max_delta_j_beam) {
                        local_max_delta_j_beam = delta_beam;
                        local_max_delta_j_beam_face =
                            static_cast<int>(iface);
                    }
                }
            }
            double local_max_j_total = 0.0;
            for (size_t iface = 0; iface < current_trial_next.size(); ++iface) {
                local_max_j_total =
                    std::max(local_max_j_total,
                             std::fabs(current_trial_next[iface]));
            }
            double local_max_j_bkg = 0.0;
            for (size_t iface = 0; iface < bkg_e.current_face_x.size();
                 ++iface) {
                local_max_j_bkg =
                    std::max(local_max_j_bkg,
                             std::fabs(bkg_e.current_face_x[iface]));
            }
            double local_max_j_beam = 0.0;
            for (size_t iface = 0; iface < beam.current_face_x.size();
                 ++iface) {
                local_max_j_beam =
                    std::max(local_max_j_beam,
                             std::fabs(beam.current_face_x[iface]));
            }
            const double local_iter_errors[8] = {
                local_max_delta_ex,
                local_max_abs_ex,
                local_max_j_total,
                local_max_j_bkg,
                local_max_j_beam,
                local_max_delta_j_total,
                local_max_delta_j_bkg,
                local_max_delta_j_beam
            };
            double global_iter_errors[8] = {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            };
            MPI_Allreduce(local_iter_errors, global_iter_errors, 8,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            const double global_max_delta_ex = global_iter_errors[0];
            const double global_max_abs_ex = global_iter_errors[1];
            const double global_field_error =
                global_max_delta_ex /
                std::max(field_floor, global_max_abs_ex);
            const double global_current_scale =
                std::max(std::fabs(Param::jb), global_iter_errors[2]);
            const double global_j_abs_tol =
                current_tol * std::max(1.0, global_current_scale);
            const double global_e_abs_tol =
                std::max(field_tol * field_floor,
                         0.5 * dt / Const::eps0 * global_j_abs_tol);
            const double global_max_delta_j_total = global_iter_errors[5];
            const double global_j_total_error =
                global_max_delta_j_total / std::max(1.0, global_current_scale);
            const double global_j_bkg_error =
                global_iter_errors[6] / std::max(1.0, global_current_scale);
            const double global_j_beam_error =
                global_iter_errors[7] / std::max(1.0, global_current_scale);
            const double global_field_convergence_ratio =
                std::min(global_field_error / field_tol,
                         global_max_delta_ex /
                         std::max(field_tol * field_floor, global_e_abs_tol));
            const double global_normalized_error =
                have_previous_current_trial
                ? std::max(global_field_convergence_ratio,
                           global_max_delta_j_total /
                           std::max(1.0, global_j_abs_tol))
                : std::max(1.0, global_field_convergence_ratio);

            const int total_global_face =
                sgrid.ix_start + local_max_delta_j_total_face;
            const int bkg_global_face =
                sgrid.ix_start + local_max_delta_j_bkg_face;
            const int beam_global_face =
                sgrid.ix_start + local_max_delta_j_beam_face;
            const double local_location_values[12] = {
                local_max_delta_j_total,
                static_cast<double>(local_max_delta_j_total_face),
                static_cast<double>(total_global_face),
                total_global_face * sgrid.dx,
                local_max_delta_j_bkg,
                static_cast<double>(local_max_delta_j_bkg_face),
                static_cast<double>(bkg_global_face),
                bkg_global_face * sgrid.dx,
                local_max_delta_j_beam,
                static_cast<double>(local_max_delta_j_beam_face),
                static_cast<double>(beam_global_face),
                beam_global_face * sgrid.dx
            };
            std::vector<double> gathered_location_values;
            if (mpi_rank == 0) {
                gathered_location_values.assign(static_cast<size_t>(mpi_size) * 12,
                                                0.0);
            }
            MPI_Gather(local_location_values, 12, MPI_DOUBLE,
                       mpi_rank == 0 ? gathered_location_values.data() : NULL,
                       12, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (mpi_rank == 0) {
                double global_delta_j_total = -1.0;
                int global_delta_j_total_rank = 0;
                int global_delta_j_total_face = 0;
                int global_delta_j_total_global_face = 0;
                double global_delta_j_total_x = 0.0;
                double global_delta_j_bkg = -1.0;
                int global_delta_j_bkg_rank = 0;
                int global_delta_j_bkg_face = 0;
                int global_delta_j_bkg_global_face = 0;
                double global_delta_j_bkg_x = 0.0;
                double global_delta_j_beam = -1.0;
                int global_delta_j_beam_rank = 0;
                int global_delta_j_beam_face = 0;
                int global_delta_j_beam_global_face = 0;
                double global_delta_j_beam_x = 0.0;
                for (int rank = 0; rank < mpi_size; ++rank) {
                    const size_t base = static_cast<size_t>(rank) * 12;
                    if (gathered_location_values[base] >
                        global_delta_j_total) {
                        global_delta_j_total =
                            gathered_location_values[base];
                        global_delta_j_total_rank = rank;
                        global_delta_j_total_face =
                            static_cast<int>(gathered_location_values[base + 1]);
                        global_delta_j_total_global_face =
                            static_cast<int>(gathered_location_values[base + 2]);
                        global_delta_j_total_x =
                            gathered_location_values[base + 3];
                    }
                    if (gathered_location_values[base + 4] >
                        global_delta_j_bkg) {
                        global_delta_j_bkg =
                            gathered_location_values[base + 4];
                        global_delta_j_bkg_rank = rank;
                        global_delta_j_bkg_face =
                            static_cast<int>(gathered_location_values[base + 5]);
                        global_delta_j_bkg_global_face =
                            static_cast<int>(gathered_location_values[base + 6]);
                        global_delta_j_bkg_x =
                            gathered_location_values[base + 7];
                    }
                    if (gathered_location_values[base + 8] >
                        global_delta_j_beam) {
                        global_delta_j_beam =
                            gathered_location_values[base + 8];
                        global_delta_j_beam_rank = rank;
                        global_delta_j_beam_face =
                            static_cast<int>(gathered_location_values[base + 9]);
                        global_delta_j_beam_global_face =
                            static_cast<int>(gathered_location_values[base + 10]);
                        global_delta_j_beam_x =
                            gathered_location_values[base + 11];
                    }
                }
                std::printf("coupled_iter step=%d iter=%d "
                            "field_error=%.6e J_total_error=%.6e "
                            "J_bkg_error=%.6e J_beam_error=%.6e "
                            "E_abs_tol=%.6e max_delta_Ex=%.6e "
                            "J_abs_tol=%.6e current_scale=%.6e "
                            "max_abs_Ex_mid=%.6e max_abs_J_bkg=%.6e "
                            "max_abs_J_beam=%.6e "
                            "max_delta_J_total_face=%.6e "
                            "total_rank=%d total_face=%d "
                            "total_global_face=%d total_x=%.6e "
                            "max_delta_J_bkg_face=%.6e "
                            "bkg_rank=%d bkg_face=%d bkg_global_face=%d "
                            "bkg_x=%.6e "
                            "max_delta_J_beam_face=%.6e "
                            "beam_rank=%d beam_face=%d beam_global_face=%d "
                            "beam_x=%.6e omega=%.3f\n",
                            step, coupled_iter + 1,
                            global_field_error, global_j_total_error,
                            global_j_bkg_error, global_j_beam_error,
                            global_e_abs_tol, global_max_delta_ex,
                            global_j_abs_tol, global_current_scale,
                            global_iter_errors[1], global_iter_errors[3],
                            global_iter_errors[4],
                            global_delta_j_total,
                            global_delta_j_total_rank,
                            global_delta_j_total_face,
                            global_delta_j_total_global_face,
                            global_delta_j_total_x,
                            global_delta_j_bkg, global_delta_j_bkg_rank,
                            global_delta_j_bkg_face,
                            global_delta_j_bkg_global_face,
                            global_delta_j_bkg_x,
                            global_delta_j_beam, global_delta_j_beam_rank,
                            global_delta_j_beam_face,
                            global_delta_j_beam_global_face,
                            global_delta_j_beam_x, midpoint_omega);
                std::fflush(stdout);
            }

            final_coupled_error = global_normalized_error;
            final_field_error = global_field_error;
            final_current_error = global_j_total_error;
            final_max_delta_ex = global_max_delta_ex;
            final_e_abs_tol = global_e_abs_tol;
            final_j_abs_tol = global_j_abs_tol;
            if (have_previous_current_trial &&
                (global_field_error < field_tol ||
                 global_max_delta_ex < global_e_abs_tol) &&
                global_max_delta_j_total < global_j_abs_tol) {
                coupled_converged = true;
                accepted_bkg = bkg_e;
                accepted_beam = beam;
                accepted_fields = fields;
                break;
            }
            if (previous_normalized_error > 0.0) {
                if (global_normalized_error < 0.9 * previous_normalized_error) {
                    ++residual_decrease_count;
                    const double omega_cap =
                        (residual_decrease_count >= 2) ? 1.0 : 0.8;
                    midpoint_omega =
                        std::min(omega_cap, midpoint_omega + 0.1);
                } else if (global_normalized_error >
                           1.02 * previous_normalized_error) {
                    residual_decrease_count = 0;
                    midpoint_omega = std::max(0.3, 0.5 * midpoint_omega);
                } else {
                    residual_decrease_count = 0;
                    midpoint_omega = std::max(0.3, 0.8 * midpoint_omega);
                }
            }
            previous_normalized_error = global_normalized_error;
            const size_t relax_faces =
                std::min(ex_mid_trial.size(), ex_mid_next.size());
            for (size_t iface = 0; iface < relax_faces; ++iface) {
                ex_mid_trial[iface] =
                    (1.0 - midpoint_omega) * ex_mid_trial[iface]
                  + midpoint_omega * ex_mid_next[iface];
            }
            current_trial_prev.swap(current_trial_next);
            current_bkg_trial_prev.swap(current_bkg_trial_next);
            current_beam_trial_prev.swap(current_beam_trial_next);
            have_previous_current_trial = true;
        }

        if (!coupled_converged) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "ERROR: coupled midpoint iteration failed to "
                             "converge at step %d, t = %.6e s; normalized "
                             "residual %.6e, field_error %.6e, "
                             "max_delta_Ex %.6e, E_abs_tol %.6e, "
                             "current_error %.6e, J_abs_tol %.6e after %d "
                             "iterations. Reduce dt or increase the coupled "
                             "solve robustness.\n",
                             step, time, final_coupled_error,
                             final_field_error, final_max_delta_ex,
                             final_e_abs_tol, final_current_error,
                             final_j_abs_tol,
                             max_coupled_iters);
            }
            MPI_Abort(MPI_COMM_WORLD, 8);
        }

        bkg_e = accepted_bkg;
        beam = accepted_beam;
        fields = accepted_fields;
        moments_current = true;

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
                           sgrid, mpi_rank, mpi_size, config.enable_full_fe_output);
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
                       sgrid, mpi_rank, mpi_size, config.enable_full_fe_output);
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
