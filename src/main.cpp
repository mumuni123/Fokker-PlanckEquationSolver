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
    double a_est = std::abs(electron.charge) * e_est / electron.mass;
    if (a_est > 1.0e-30) {
        dt_min = std::min(dt_min, 0.25 * electron.vgrid.dv / a_est);
    }
    dt_min *= Param::dt_multiplier;
    return std::min(dt_min, 0.01 * Const::femto);
}

const char* poisson_solver_name()
{
    return "FFT periodic discrete Poisson with zero-mode removal";
}

std::vector<double> build_local_ion_density_profile(const SpatialGrid& sg)
{
    return std::vector<double>(static_cast<size_t>(sg.nx_local), Param::dens);
}

void sync_moments_and_fields(Species& electrons,
                             const BeamPIC& beam,
                             EMFields& fields,
                             const std::vector<double>& ion_density_profile,
                             int mpi_rank,
                             int mpi_size,
                             bool& moments_current)
{
    electrons.compute_moments();
    moments_current = true;
    fields.set_charge_density(electrons, beam.density, ion_density_profile);
    fields.solve_poisson(mpi_rank, mpi_size);
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
        Param::vmax_loss_abort_fraction * Param::dens * Param::plasma_length;
    if (!(global_loss <= threshold)) {
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "ERROR: background electron distribution reached vmax "
                         "during %s at step %d, t = %.6e s. "
                         "loss_v_high = %.8e, threshold = %.8e. "
                         "Stopping to avoid amplifying nonphysical results.\n",
                         stage, step, time, global_loss, threshold);
        }
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
}

std::vector<double> copy_local_ex(const EMFields& fields,
                                  const SpatialGrid& sg)
{
    std::vector<double> ex(static_cast<size_t>(sg.nx_local), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        ex[static_cast<size_t>(ix)] = fields.Ex[ix + sg.nghost];
    }
    return ex;
}

double integrate_current_work(const std::vector<double>& current_before,
                              const std::vector<double>& current_after,
                              const std::vector<double>& ex,
                              const SpatialGrid& sg,
                              double dt)
{
    double work = 0.0;
    const size_t n = std::min(current_before.size(),
                              std::min(current_after.size(), ex.size()));
    for (size_t ix = 0; ix < n; ++ix) {
        const double current_mid = 0.5 * (current_before[ix] + current_after[ix]);
        work += current_mid * ex[ix] * sg.dx * dt;
    }
    return work;
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
        printf("  Spherical electron velocity grid: (v, mu), mu = vx / |v|\n");
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
        printf("Electron velocity grid: Nv x Nmu = %d x %d\n", Param::Nv, Param::Nmu);
        printf("Electron velocity domain: 0 <= v <= %.1f v_th (cap %.3f c)\n",
               Param::Nsigma, Param::vmax_fraction_c);
        printf("Electrostatic boundary: periodic discrete FFT, net charge forced to zero before each solve\n");
        printf("Poisson solver: %s\n", poisson_solver_name());
        printf("Fixed ions: uniform Z*n_i = %.3e /m^3\n", Param::dens);
        printf("Background electrons: full-domain Maxwellian, T_e = %.1f eV, periodic ghosts\n",
               Param::temperature_e / Const::eV);
        printf("PIC beam: gamma*beta = %.2f, beta = %.4f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::betab, Param::densb);
        printf("Beam source: quiet-start internal flux plane at x = %.3f um\n",
               Param::beam_source_x_start / Const::micro);
        printf("Beam injection: charge-conserving path current, centered before Poisson\n");
        printf("Beam boundary: particles crossing the domain edge are deleted and counted in the energy ledger\n");
        printf("Return current: self-consistent background-electron response only\n");
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
    bool moments_current = true;
    sync_moments_and_fields(bkg_e, beam, fields, ion_density_profile,
                            mpi_rank, mpi_size, moments_current);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(0.0, 0, bkg_e, beam, fields,
                       cumulative_collision_energy_delta, mpi_rank, mpi_size);
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
        double loss_x1_left = 0.0;
        double loss_x1_right = 0.0;
        double loss_x2_left = 0.0;
        double loss_x2_right = 0.0;
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
        double beam_push_ke_before = 0.0;
        double bkg_ke_step_start = 0.0;
        double beam_ke_step_start = 0.0;
        double field_energy_step_start = 0.0;
        std::vector<double> bkg_current_start;
        std::vector<double> bkg_current_mid;
        std::vector<double> ex_step_start;
        std::vector<double> ex_center;
        vlasov.set_step_diagnostics_enabled(collect_step_diagnostics);

        beam.begin_step(sgrid, dt);

        if (collect_step_diagnostics) {
            bkg_ke_step_start = bkg_e.total_kinetic_energy();
            beam_ke_step_start = beam.total_kinetic_energy();
            field_energy_step_start = fields.total_energy();
            bkg_current_start = bkg_e.current_x;
            ex_step_start = copy_local_ex(fields, sgrid);
        }

        trace_progress(config, mpi_rank, step, "before x half 1");
        vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after x half 1");
        loss_x1_left = vlasov.last_loss_x_left();
        loss_x1_right = vlasov.last_loss_x_right();
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "x_half_1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);
        }
#endif

        trace_progress(config, mpi_rank, step, "before v half 1");
        vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
        trace_progress(config, mpi_rank, step, "after v half 1");
        nsub_v1 = vlasov.last_nsub_v();
        loss_v1 = vlasov.last_loss_v();
        loss_v1_low = vlasov.last_loss_v_low();
        loss_v1_high = vlasov.last_loss_v_high();
        v_mass_error_step += vlasov.last_mass_error_v();
        v_momentum_delta_step += vlasov.last_momentum_delta_v();
        v_energy_delta_step += vlasov.last_energy_delta_v();
        abort_if_vmax_loss(vlasov, step, time, "v_half_1", mpi_rank);
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "v_half_1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), 0.0,
                                   vlasov.last_nsub_v(), 0);
        }
#endif

        trace_progress(config, mpi_rank, step, "before mu half 1");
        vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
        trace_progress(config, mpi_rank, step, "after mu half 1");
        nsub_mu1 = vlasov.last_nsub_mu();
        loss_mu1 = vlasov.last_loss_mu();
        mu_mass_error_step += vlasov.last_mass_error_mu();
        mu_momentum_delta_step += vlasov.last_momentum_delta_mu();
        mu_energy_delta_step += vlasov.last_energy_delta_mu();
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "mu_half_1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   0.0, vlasov.last_cfl_mu(),
                                   0, vlasov.last_nsub_mu());
        }
#endif

        trace_progress(config, mpi_rank, step, "before beam half 1 push");
        if (collect_step_diagnostics) {
            beam_push_ke_before = beam.total_kinetic_energy();
        }
        beam.push(sgrid, fields, 0.5 * dt, mpi_rank, mpi_size);
        if (collect_step_diagnostics) {
            dke_beam_push += beam.total_kinetic_energy() - beam_push_ke_before;
        }
        trace_progress(config, mpi_rank, step, "after beam half 1 push");
        trace_progress(config, mpi_rank, step, "before beam half 1 inject");
        beam.inject(sgrid, fields, 0.5 * dt, time_center, mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after beam half 1 inject");
        trace_progress(config, mpi_rank, step, "before beam center deposit");
        beam.deposit_density(sgrid, mpi_rank, mpi_size);
        beam.finalize_charge_conserving_current(sgrid, 0.5 * dt,
                                                mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after beam center deposit");

        trace_progress(config, mpi_rank, step, "before center moments");
        bkg_e.compute_moments();
        moments_current = true;
        trace_progress(config, mpi_rank, step, "after center moments");
        if (collect_step_diagnostics) {
            W_bkg_E +=
                integrate_current_work(bkg_current_start, bkg_e.current_x,
                                       ex_step_start, sgrid, 0.5 * dt);
            bkg_current_mid = bkg_e.current_x;
        }
        trace_progress(config, mpi_rank, step, "before center Ex solve");
        fields.set_charge_density(bkg_e, beam.density, ion_density_profile);
        fields.solve_poisson(mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after center Ex solve");
        if (collect_step_diagnostics) {
            ex_center = copy_local_ex(fields, sgrid);
        }
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "center_solve_Ex",
                                   bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);
        }
#endif

        trace_progress(config, mpi_rank, step, "before beam half 2 push");
        if (collect_step_diagnostics) {
            beam_push_ke_before = beam.total_kinetic_energy();
        }
        beam.push(sgrid, fields, 0.5 * dt, mpi_rank, mpi_size);
        if (collect_step_diagnostics) {
            dke_beam_push += beam.total_kinetic_energy() - beam_push_ke_before;
        }
        trace_progress(config, mpi_rank, step, "after beam half 2 push");
        trace_progress(config, mpi_rank, step, "before beam half 2 inject");
        beam.inject(sgrid, fields, 0.5 * dt, time, mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after beam half 2 inject");
        trace_progress(config, mpi_rank, step, "before beam end deposit");
        beam.deposit_density(sgrid, mpi_rank, mpi_size);
        beam.finalize_charge_conserving_current(sgrid, dt,
                                                mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after beam end deposit");
        if (collect_step_diagnostics) {
            W_beam_E = beam.last_field_work();
        }

        trace_progress(config, mpi_rank, step, "before mu half 2");
        vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
        trace_progress(config, mpi_rank, step, "after mu half 2");
        nsub_mu2 = vlasov.last_nsub_mu();
        loss_mu2 = vlasov.last_loss_mu();
        mu_mass_error_step += vlasov.last_mass_error_mu();
        mu_momentum_delta_step += vlasov.last_momentum_delta_mu();
        mu_energy_delta_step += vlasov.last_energy_delta_mu();
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "mu_half_2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   0.0, vlasov.last_cfl_mu(),
                                   0, vlasov.last_nsub_mu());
        }
#endif

        trace_progress(config, mpi_rank, step, "before v half 2");
        vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
        trace_progress(config, mpi_rank, step, "after v half 2");
        nsub_v2 = vlasov.last_nsub_v();
        loss_v2 = vlasov.last_loss_v();
        loss_v2_low = vlasov.last_loss_v_low();
        loss_v2_high = vlasov.last_loss_v_high();
        v_mass_error_step += vlasov.last_mass_error_v();
        v_momentum_delta_step += vlasov.last_momentum_delta_v();
        v_energy_delta_step += vlasov.last_energy_delta_v();
        abort_if_vmax_loss(vlasov, step, time, "v_half_2", mpi_rank);
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "v_half_2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), 0.0,
                                   vlasov.last_nsub_v(), 0);
        }
#endif

        trace_progress(config, mpi_rank, step, "before x half 2");
        vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after x half 2");
        loss_x2_left = vlasov.last_loss_x_left();
        loss_x2_right = vlasov.last_loss_x_right();
        moments_current = false;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (config.enable_debug_diagnostics) {
            diag.write_debug_state(step, time, "x_half_2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);
        }
#endif

        if (collect_step_diagnostics) {
            bkg_e.compute_moments();
            W_bkg_E +=
                integrate_current_work(bkg_current_mid, bkg_e.current_x,
                                       ex_center, sgrid, 0.5 * dt);
            moments_current = true;
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
        sync_moments_and_fields(bkg_e, beam, fields, ion_density_profile,
                                mpi_rank, mpi_size, moments_current);
        trace_progress(config, mpi_rank, step, "after end sync");

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
                                        loss_x1_left, loss_x1_right,
                                        loss_x2_left, loss_x2_right,
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
                                        E_src_in_step, E_src_out_step,
                                        E_balance_step);
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

    sync_moments_and_fields(bkg_e, beam, fields, ion_density_profile,
                            mpi_rank, mpi_size, moments_current);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(nsteps, Param::t_end, "final", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(Param::t_end, nsteps, bkg_e, beam, fields,
                       cumulative_collision_energy_delta, mpi_rank, mpi_size);
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
