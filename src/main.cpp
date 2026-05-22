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
        dt_min = std::min(dt_min, 0.25 * electron.vgrid.dmu * electron.vgrid.dv / a_est);
    }
    dt_min *= Param::dt_multiplier;
    return std::min(dt_min, 0.01 * Const::femto);
}

const char* poisson_solver_name()
{
    switch (Param::poisson_solver) {
    case Param::POISSON_PARALLEL_CYCLIC_REDUCTION:
        return "parallel cyclic reduction";
    case Param::POISSON_MULTIGRID:
        return "multigrid";
    case Param::POISSON_DISTRIBUTED_TRIDIAGONAL:
    default:
        return "distributed tridiagonal";
    }
}

void write_snapshot(Diagnostics& diag,
                    double time,
                    const Species& bkg_e,
                    const BeamPIC& beam,
                    const EMFields& fields,
                    const SpatialGrid& sgrid,
                    int mpi_rank,
                    int mpi_size,
                    bool write_full_fe)
{
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_density_profile(time, bkg_e, beam.density, sgrid, mpi_rank, mpi_size);
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
        printf("Spatial grid: nx = %d, dx = %.3e m\n", Param::nx, Param::dx);
        printf("Electron velocity grid: Nv x Nmu = %d x %d\n", Param::Nv, Param::Nmu);
        printf("Electron velocity domain: 0 <= v <= %.1f v_th (cap 0.98 c)\n",
               Param::Nsigma);
        printf("Electrostatic boundary: phi(0) = phi(L) = 0\n");
        printf("Poisson solver: %s\n", poisson_solver_name());
        printf("Fixed ions: Z*n_i = %.3e /m^3\n", Param::dens);
        printf("Background electrons: n_e0 = %.3e /m^3, T_e = %.1f eV\n",
               Param::dens, Param::temperature_e / Const::eV);
        printf("PIC beam: gamma*beta = %.2f, beta = %.4f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::betab, Param::densb);
        printf("Beam macro weight: %.6e particles/m^2\n", Param::beam_macro_weight);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        printf("Debug diagnostics: %s\n",
               config.enable_debug_diagnostics ? "ON" : "OFF");
#else
        printf("Debug diagnostics: compile-time disabled\n");
#endif
        printf("Full fe distribution output: %s\n",
               config.enable_full_fe_output ? "ON" : "OFF");
        printf("Progress trace: %s\n",
               config.enable_progress_trace ? "ON" : "OFF");
        printf("------------------------------------------------------------\n");
    }

    SpatialGrid sgrid;
    sgrid.init(mpi_rank, mpi_size);

    Species bkg_e;
    bkg_e.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
               -Const::qe, Const::me,
               Param::dens, Param::temperature_e, false, sgrid);
    bkg_e.initialize_maxwellian();

    BeamPIC beam;
    beam.init(sgrid);
    beam.deposit_density(sgrid, mpi_rank, mpi_size);

    std::vector<Species*> electron_species;
    electron_species.push_back(&bkg_e);

    EMFields fields;
    fields.init(sgrid);

    VlasovSolver vlasov;
    CollisionOperator collision;
    Diagnostics diag;
    diag.init("output", mpi_rank, config.enable_debug_diagnostics);

    double dt = compute_dt(bkg_e, sgrid);
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    int nsteps = static_cast<int>(std::ceil(Param::t_end / dt));

    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    bkg_e.compute_moments();
    bool moments_current = true;
    fields.set_charge_density(bkg_e, beam.density);
    fields.solve_poisson(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(0.0, 0, electron_species, fields, mpi_rank, mpi_size);
    write_snapshot(diag, 0.0, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size,
                   config.enable_full_fe_output);

    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 1000;
    int last_snapshot_step = 0;

#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        for (int step = 1; step <= nsteps; ++step) {
            double time = step * dt;

            trace_progress(config, mpi_rank, step, "before x1");
            vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after x1");
            moments_current = false;
            diag.write_debug_state(step, time, "x1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "before v1");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after v1");
            moments_current = false;
            diag.write_debug_state(step, time, "v1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), 0.0,
                                   vlasov.last_nsub_v(), 0);
            trace_progress(config, mpi_rank, step, "before mu1");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after mu1");
            moments_current = false;
            diag.write_debug_state(step, time, "mu1", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), vlasov.last_cfl_mu(),
                                   vlasov.last_nsub_v(), vlasov.last_nsub_mu());

            trace_progress(config, mpi_rank, step, "before beam push");
            beam.push(sgrid, fields, dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after beam push");
            trace_progress(config, mpi_rank, step, "before beam inject");
            beam.inject(sgrid, dt, time, mpi_rank);
            trace_progress(config, mpi_rank, step, "after beam inject");
            trace_progress(config, mpi_rank, step, "before beam deposit");
            beam.deposit_density(sgrid, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after beam deposit");

            if (!moments_current) {
                trace_progress(config, mpi_rank, step, "before moments");
                bkg_e.compute_moments();
                trace_progress(config, mpi_rank, step, "after moments");
                moments_current = true;
            }
            trace_progress(config, mpi_rank, step, "before poisson");
            fields.set_charge_density(bkg_e, beam.density);
            fields.solve_poisson(mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after poisson");
            diag.write_debug_state(step, time, "solve_poisson", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);

            if (bkg_e.collisions_enabled) {
                trace_progress(config, mpi_rank, step, "before collisions");
                collision.apply(bkg_e, dt, Param::dens, Param::temperature_e,
                                Const::me, 1.0, 1.0);
                collision.apply(bkg_e, dt, Param::dens / Param::Z_ion,
                                Param::temperature_i, Param::mass_ion,
                                (double)Param::Z_ion, 1.0);
                trace_progress(config, mpi_rank, step, "after collisions");
                moments_current = false;
            }

            trace_progress(config, mpi_rank, step, "before x2");
            vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after x2");
            moments_current = false;
            diag.write_debug_state(step, time, "x2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "before v2");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after v2");
            moments_current = false;
            diag.write_debug_state(step, time, "v2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), 0.0,
                                   vlasov.last_nsub_v(), 0);
            trace_progress(config, mpi_rank, step, "before mu2");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after mu2");
            moments_current = false;
            diag.write_debug_state(step, time, "mu2", bkg_e, beam, fields,
                                   sgrid, mpi_rank, mpi_size,
                                   vlasov.last_cfl_v(), vlasov.last_cfl_mu(),
                                   vlasov.last_nsub_v(), vlasov.last_nsub_mu());

            if (step % stdout_freq == 0) {
                if (!moments_current) {
                    bkg_e.compute_moments();
                    moments_current = true;
                }
                diag.write_scalars(time, step, electron_species, fields, mpi_rank, mpi_size);
                if (mpi_rank == 0) {
                    printf("Step %d / %d, t = %.4f fs\n", step, nsteps, time / Const::femto);
                }
            }

            if (time >= next_snapshot) {
                if (!moments_current) {
                    bkg_e.compute_moments();
                    moments_current = true;
                }
                write_snapshot(diag, time, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size,
                               config.enable_full_fe_output);
                last_snapshot_step = step;
                next_snapshot += Param::dt_snapshot;
            }
        }
    } else
#endif
    {
        for (int step = 1; step <= nsteps; ++step) {
            double time = step * dt;

            trace_progress(config, mpi_rank, step, "before x1");
            vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after x1");
            moments_current = false;
            trace_progress(config, mpi_rank, step, "before v1");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after v1");
            moments_current = false;
            trace_progress(config, mpi_rank, step, "before mu1");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after mu1");
            moments_current = false;

            trace_progress(config, mpi_rank, step, "before beam push");
            beam.push(sgrid, fields, dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after beam push");
            trace_progress(config, mpi_rank, step, "before beam inject");
            beam.inject(sgrid, dt, time, mpi_rank);
            trace_progress(config, mpi_rank, step, "after beam inject");
            trace_progress(config, mpi_rank, step, "before beam deposit");
            beam.deposit_density(sgrid, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after beam deposit");

            if (!moments_current) {
                trace_progress(config, mpi_rank, step, "before moments");
                bkg_e.compute_moments();
                trace_progress(config, mpi_rank, step, "after moments");
                moments_current = true;
            }
            trace_progress(config, mpi_rank, step, "before poisson");
            fields.set_charge_density(bkg_e, beam.density);
            fields.solve_poisson(mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after poisson");

            if (bkg_e.collisions_enabled) {
                trace_progress(config, mpi_rank, step, "before collisions");
                collision.apply(bkg_e, dt, Param::dens, Param::temperature_e,
                                Const::me, 1.0, 1.0);
                collision.apply(bkg_e, dt, Param::dens / Param::Z_ion,
                                Param::temperature_i, Param::mass_ion,
                                (double)Param::Z_ion, 1.0);
                trace_progress(config, mpi_rank, step, "after collisions");
                moments_current = false;
            }

            trace_progress(config, mpi_rank, step, "before x2");
            vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
            trace_progress(config, mpi_rank, step, "after x2");
            moments_current = false;
            trace_progress(config, mpi_rank, step, "before v2");
            vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after v2");
            moments_current = false;
            trace_progress(config, mpi_rank, step, "before mu2");
            vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
            trace_progress(config, mpi_rank, step, "after mu2");
            moments_current = false;

            if (step % stdout_freq == 0) {
                if (!moments_current) {
                    bkg_e.compute_moments();
                    moments_current = true;
                }
                diag.write_scalars(time, step, electron_species, fields, mpi_rank, mpi_size);
                if (mpi_rank == 0) {
                    printf("Step %d / %d, t = %.4f fs\n", step, nsteps, time / Const::femto);
                }
            }

            if (time >= next_snapshot) {
                if (!moments_current) {
                    bkg_e.compute_moments();
                    moments_current = true;
                }
                write_snapshot(diag, time, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size,
                               config.enable_full_fe_output);
                last_snapshot_step = step;
                next_snapshot += Param::dt_snapshot;
            }
        }
    }

    if (!moments_current) {
        bkg_e.compute_moments();
        moments_current = true;
    }
    fields.set_charge_density(bkg_e, beam.density);
    fields.solve_poisson(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(nsteps, Param::t_end, "final", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(Param::t_end, nsteps, electron_species, fields, mpi_rank, mpi_size);
    if (last_snapshot_step != nsteps) {
        write_snapshot(diag, Param::t_end, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size,
                       config.enable_full_fe_output);
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
