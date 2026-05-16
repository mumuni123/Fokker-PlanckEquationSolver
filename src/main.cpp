#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "vlasov.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"
#include "beam_pic.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
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

bool parse_debug_diagnostics_flag(int argc, char** argv)
{
    bool enabled = Param::enable_debug_diagnostics;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--debug-diagnostics") == 0 ||
            std::strcmp(argv[i], "--debug") == 0) {
            enabled = true;
        } else if (std::strcmp(argv[i], "--no-debug-diagnostics") == 0 ||
                   std::strcmp(argv[i], "--no-debug") == 0) {
            enabled = false;
        }
    }
    return enabled;
}

void write_snapshot(Diagnostics& diag,
                    double time,
                    const Species& bkg_e,
                    const BeamPIC& beam,
                    const EMFields& fields,
                    const SpatialGrid& sgrid,
                    int mpi_rank,
                    int mpi_size)
{
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_density_profile(time, bkg_e, beam.density, sgrid, mpi_rank, mpi_size);
    diag.write_px_distribution(time, bkg_e, mpi_rank, mpi_size);
    diag.write_electron_distribution(time, bkg_e, sgrid, mpi_rank);

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

    bool enable_debug_diagnostics = parse_debug_diagnostics_flag(argc, argv);

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
        printf("Fixed ions: Z*n_i = %.3e /m^3\n", Param::dens);
        printf("Background electrons: n_e0 = %.3e /m^3, T_e = %.1f eV\n",
               Param::dens, Param::temperature_e / Const::eV);
        printf("PIC beam: gamma*beta = %.2f, beta = %.4f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::betab, Param::densb);
        printf("Beam macro weight: %.6e particles/m^2\n", Param::beam_macro_weight);
        printf("Debug diagnostics: %s\n", enable_debug_diagnostics ? "ON" : "OFF");
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
    diag.init("output", mpi_rank, enable_debug_diagnostics);

    double dt = compute_dt(bkg_e, sgrid);
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    int nsteps = static_cast<int>(std::ceil(Param::t_end / dt));

    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    bkg_e.compute_moments();
    fields.set_charge_density(bkg_e, beam.density);
    fields.solve_poisson(mpi_rank, mpi_size);
    diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                           sgrid, mpi_rank, mpi_size);
    diag.write_scalars(0.0, 0, electron_species, fields, mpi_rank, mpi_size);
    write_snapshot(diag, 0.0, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size);

    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 200;
    int last_snapshot_step = 0;

    for (int step = 1; step <= nsteps; ++step) {
        double time = step * dt;

        vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
        diag.write_debug_state(step, time, "x1", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
        vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
        diag.write_debug_state(step, time, "v1", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size,
                               vlasov.last_cfl_v(), 0.0,
                               vlasov.last_nsub_v(), 0);
        vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
        diag.write_debug_state(step, time, "mu1", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size,
                               vlasov.last_cfl_v(), vlasov.last_cfl_mu(),
                               vlasov.last_nsub_v(), vlasov.last_nsub_mu());

        beam.push(sgrid, fields, dt, mpi_rank, mpi_size);
        beam.inject(sgrid, dt, time, mpi_rank);
        beam.deposit_density(sgrid, mpi_rank, mpi_size);

        bkg_e.compute_moments();
        fields.set_charge_density(bkg_e, beam.density);
        fields.solve_poisson(mpi_rank, mpi_size);
        diag.write_debug_state(step, time, "solve_poisson", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);

        if (bkg_e.collisions_enabled) {
            collision.apply(bkg_e, dt, Param::dens, Param::temperature_e,
                            Const::me, 1.0, 1.0);
            collision.apply(bkg_e, dt, Param::dens / Param::Z_ion,
                            Param::temperature_i, Param::mass_ion,
                            (double)Param::Z_ion, 1.0);
        }

        vlasov.advect_x(bkg_e, sgrid, 0.5 * dt, mpi_rank, mpi_size);
        diag.write_debug_state(step, time, "x2", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
        vlasov.advect_v(bkg_e, sgrid, fields, 0.5 * dt);
        diag.write_debug_state(step, time, "v2", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size,
                               vlasov.last_cfl_v(), 0.0,
                               vlasov.last_nsub_v(), 0);
        vlasov.advect_mu(bkg_e, sgrid, fields, 0.5 * dt);
        diag.write_debug_state(step, time, "mu2", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size,
                               vlasov.last_cfl_v(), vlasov.last_cfl_mu(),
                               vlasov.last_nsub_v(), vlasov.last_nsub_mu());

        if (step % stdout_freq == 0) {
            bkg_e.compute_moments();
            diag.write_scalars(time, step, electron_species, fields, mpi_rank, mpi_size);
            if (mpi_rank == 0) {
                printf("Step %d / %d, t = %.4f fs\n", step, nsteps, time / Const::femto);
            }
        }

        if (time >= next_snapshot) {
            bkg_e.compute_moments();
            write_snapshot(diag, time, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size);
            last_snapshot_step = step;
            next_snapshot += Param::dt_snapshot;
        }
    }

    bkg_e.compute_moments();
    fields.set_charge_density(bkg_e, beam.density);
    fields.solve_poisson(mpi_rank, mpi_size);
    diag.write_debug_state(nsteps, Param::t_end, "final", bkg_e, beam, fields,
                           sgrid, mpi_rank, mpi_size);
    diag.write_scalars(Param::t_end, nsteps, electron_species, fields, mpi_rank, mpi_size);
    if (last_snapshot_step != nsteps) {
        write_snapshot(diag, Param::t_end, bkg_e, beam, fields, sgrid, mpi_rank, mpi_size);
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
