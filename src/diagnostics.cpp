#include "diagnostics.h"
#include "beam_pic.h"
#include "fft_utils.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iomanip>
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
      energy_reference(0.0)
{}

LowModeFractions::LowModeFractions()
{
    for (int i = 0; i < 3; ++i) {
        rho[i] = 0.0;
        Ex[i] = 0.0;
    }
}

namespace {
void fill_low_mode_fractions(std::vector<std::complex<double> >& spectrum,
                             double fractions[3])
{
    const int n = static_cast<int>(spectrum.size());
    if (n <= 1) return;

    fft_any(spectrum, false);

    double total_power = 0.0;
    for (int k = 1; k < n; ++k) {
        total_power += std::norm(spectrum[static_cast<size_t>(k)]);
    }
    if (!(total_power > 0.0)) return;

    for (int mode = 1; mode <= 3; ++mode) {
        if (mode >= n) break;
        double mode_power = std::norm(spectrum[static_cast<size_t>(mode)]);
        const int mirror = n - mode;
        if (mirror != mode && mirror > 0 && mirror < n) {
            mode_power += std::norm(spectrum[static_cast<size_t>(mirror)]);
        }
        fractions[mode - 1] = mode_power / total_power;
    }
}
}

LowModeFractions compute_low_mode_fractions(const EMFields& fields,
                                            const SpatialGrid& sg,
                                            int mpi_rank,
                                            int mpi_size)
{
    LowModeFractions result;

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    std::vector<double> local_rho(static_cast<size_t>(nxl), 0.0);
    std::vector<double> local_Ex(static_cast<size_t>(nxl), 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        local_rho[static_cast<size_t>(ix)] = fields.rho[ix + ng];
        local_Ex[static_cast<size_t>(ix)] = fields.Ex[ix + ng];
    }

    std::vector<int> counts(static_cast<size_t>(mpi_size), 0);
    std::vector<int> displs(static_cast<size_t>(mpi_size), 0);
    MPI_Gather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        for (int r = 1; r < mpi_size; ++r) {
            displs[static_cast<size_t>(r)] =
                displs[static_cast<size_t>(r - 1)] +
                counts[static_cast<size_t>(r - 1)];
        }
    }

    std::vector<double> global_rho;
    std::vector<double> global_Ex;
    if (mpi_rank == 0) {
        global_rho.assign(static_cast<size_t>(sg.nx_global), 0.0);
        global_Ex.assign(static_cast<size_t>(sg.nx_global), 0.0);
    }

    MPI_Gatherv(local_rho.data(), nxl, MPI_DOUBLE,
                mpi_rank == 0 ? global_rho.data() : static_cast<double*>(0),
                counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Gatherv(local_Ex.data(), nxl, MPI_DOUBLE,
                mpi_rank == 0 ? global_Ex.data() : static_cast<double*>(0),
                counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::vector<std::complex<double> > rho_spectrum(static_cast<size_t>(sg.nx_global));
        std::vector<std::complex<double> > ex_spectrum(static_cast<size_t>(sg.nx_global));
        for (int i = 0; i < sg.nx_global; ++i) {
            rho_spectrum[static_cast<size_t>(i)] =
                std::complex<double>(global_rho[static_cast<size_t>(i)], 0.0);
            ex_spectrum[static_cast<size_t>(i)] =
                std::complex<double>(global_Ex[static_cast<size_t>(i)], 0.0);
        }
        fill_low_mode_fractions(rho_spectrum, result.rho);
        fill_low_mode_fractions(ex_spectrum, result.Ex);
    }

    return result;
}

void Diagnostics::init(const std::string& dir, int mpi_rank,
                       bool enable_debug_diagnostics,
                       bool enable_step_diagnostics)
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
                    << "E_accounted[J/m2]  E_balance_error[J/m2]\n";
        scalar_file << std::scientific << std::setprecision(8);

#if FP_ENABLE_DEBUG_DIAGNOSTICS
        if (debug_enabled) {
            debug_file.open((output_dir + "/debug_diagnostics.dat").c_str());
            debug_file << "# step  time[fs]  stage  max_abs_Ex[V/m]  N_bkg_e  "
                       << "N_beam_macro  N_beam_weighted  CFL_v  CFL_mu  "
                       << "nsub_v  nsub_mu\n";
            debug_file << std::scientific << std::setprecision(8);
        }
#endif
        if (step_enabled) {
            step_file.open((output_dir + "/step_diagnostics.dat").c_str());
            step_file << "# step  time[fs]  max_abs_Ex[V/m]  N_bkg_e  "
                      << "N_beam_macro  N_beam_weighted  "
                      << "N_beam_source_step  N_beam_absorb_step  "
                      << "J_beam_source_int[A/m]  "
                      << "beam_cont_l1  beam_cont_linf  "
                      << "nsub_v1  nsub_mu1  nsub_v2  nsub_mu2  "
                      << "loss_v1  loss_mu1  loss_v2  loss_mu2  "
                      << "loss_v1_low  loss_v1_high  "
                      << "loss_v2_low  loss_v2_high  "
                      << "loss_x1_left  loss_x1_right  "
                      << "loss_x2_left  loss_x2_right  "
                      << "rho_k1_frac  rho_k2_frac  rho_k3_frac  "
                      << "Ex_k1_frac  Ex_k2_frac  Ex_k3_frac  "
                      << "KE_bkg_e[J/m2]  KE_beam[J/m2]  E_field[J/m2]  "
                      << "E_total[J/m2]  dKE_bkg[J/m2]  dKE_beam[J/m2]  "
                      << "dE_field[J/m2]  W_bkg_E[J/m2]  W_beam_E[J/m2]  "
                      << "v_mass_error  mu_mass_error  "
                      << "v_px_delta[kg/m/s/m2]  mu_px_delta[kg/m/s/m2]  "
                      << "v_energy_delta[J/m2]  mu_energy_delta[J/m2]  "
                      << "E_src_in[J/m2]  E_src_out[J/m2]  "
                      << "E_collision_step[J/m2]  E_balance_step[J/m2]  "
                      << "E_beam_injected_cum[J/m2]  "
                      << "E_beam_outflow_cum[J/m2]  E_collision_cum[J/m2]\n";
            step_file << std::scientific << std::setprecision(8);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
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

    double local_values[8] = {
        local_bkg_number,
        local_bkg_ke,
        beam.total_particle_number(*electrons.sgrid),
        beam.total_kinetic_energy(),
        fields.total_energy(),
        beam.cumulative_injected_energy(),
        beam.cumulative_outflow_energy(),
        cumulative_collision_energy_delta
    };
    double global_values[8] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    MPI_Reduce(local_values, global_values, 8, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        const double total_energy =
            global_values[1] + global_values[3] + global_values[4];
        const double accounted_energy =
            total_energy - global_values[5] + global_values[6] - global_values[7];
        if (!has_energy_reference) {
            energy_reference = accounted_energy;
            has_energy_reference = true;
        }
        const double balance_error = accounted_energy - energy_reference;

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
                    << balance_error << "\n";
        scalar_file.flush();
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
                                         double loss_x1_left,
                                         double loss_x1_right,
                                         double loss_x2_left,
                                         double loss_x2_right,
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
                                         double E_balance_step)
{
    if (!step_enabled) return;

    double local_max_abs_Ex = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_max_abs_Ex = std::max(local_max_abs_Ex,
                                    std::fabs(fields.Ex[ix + sg.nghost]));
    }

    const double local_N_bkg_e = electrons.total_particle_number();
    const double local_N_beam_macro = static_cast<double>(beam.particles.size());
    const double local_N_beam_weighted = beam.total_particle_number(sg);
    double local_beam_source[3] = { 0.0, 0.0, 0.0 };
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_beam_source[0] +=
            beam.source_density_delta[static_cast<size_t>(ix)] * sg.dx;
        local_beam_source[1] +=
            beam.absorber_density_delta[static_cast<size_t>(ix)] * sg.dx;
        local_beam_source[2] +=
            beam.source_current_x[static_cast<size_t>(ix)] * sg.dx;
    }
    double global_beam_source[3] = { 0.0, 0.0, 0.0 };
    double local_beam_continuity[2] = {
        beam.last_continuity_l1_error(),
        beam.last_continuity_linf_error()
    };
    double global_beam_continuity_l1 = 0.0;
    double global_beam_continuity_linf = 0.0;

    double local_losses[12] = {
        loss_v1, loss_mu1, loss_v2, loss_mu2,
        loss_v1_low, loss_v1_high, loss_v2_low, loss_v2_high,
        loss_x1_left, loss_x1_right, loss_x2_left, loss_x2_right
    };
    double global_losses[12] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const double local_bkg_ke = electrons.total_kinetic_energy();
    const double local_beam_ke = beam.total_kinetic_energy();
    const double local_field_energy = fields.total_energy();
    const double local_total_energy =
        local_bkg_ke + local_beam_ke + local_field_energy;
    double local_energy[22] = {
        local_bkg_ke,
        local_beam_ke,
        local_field_energy,
        local_total_energy,
        dke_bkg_step,
        dke_beam_push,
        dE_field_step,
        W_bkg_E,
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
        beam.cumulative_injected_energy(),
        beam.cumulative_outflow_energy(),
        cumulative_collision_energy_delta
    };
    double global_energy[22] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const LowModeFractions low_modes =
        compute_low_mode_fractions(fields, sg, mpi_rank, mpi_size);

    double global_max_abs_Ex = 0.0;
    double global_N_bkg_e = 0.0;
    double global_N_beam_macro = 0.0;
    double global_N_beam_weighted = 0.0;
    int local_nsub[4] = { nsub_v1, nsub_mu1, nsub_v2, nsub_mu2 };
    int global_nsub[4] = { 0, 0, 0, 0 };

    MPI_Reduce(&local_max_abs_Ex, &global_max_abs_Ex, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_bkg_e, &global_N_bkg_e, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_macro, &global_N_beam_macro, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_weighted, &global_N_beam_weighted, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_beam_source, global_beam_source, 3,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_beam_continuity[0], &global_beam_continuity_l1,
               1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_beam_continuity[1], &global_beam_continuity_linf,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_nsub, global_nsub, 4, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_losses, global_losses, 12, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_energy, global_energy, 22, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        step_file << step << "  "
                  << time / Const::femto << "  "
                  << global_max_abs_Ex << "  "
                  << global_N_bkg_e << "  "
                  << global_N_beam_macro << "  "
                  << global_N_beam_weighted << "  "
                  << global_beam_source[0] << "  "
                  << global_beam_source[1] << "  "
                  << global_beam_source[2] << "  "
                  << global_beam_continuity_l1 << "  "
                  << global_beam_continuity_linf << "  "
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
                  << global_losses[9] << "  "
                  << global_losses[10] << "  "
                  << global_losses[11] << "  "
                  << low_modes.rho[0] << "  "
                  << low_modes.rho[1] << "  "
                  << low_modes.rho[2] << "  "
                  << low_modes.Ex[0] << "  "
                  << low_modes.Ex[1] << "  "
                  << low_modes.Ex[2] << "  "
                  << global_energy[0] << "  "
                  << global_energy[1] << "  "
                  << global_energy[2] << "  "
                  << global_energy[3] << "  "
                  << global_energy[4] << "  "
                  << global_energy[5] << "  "
                  << global_energy[6] << "  "
                  << global_energy[7] << "  "
                  << global_energy[8] << "  "
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
                  << global_energy[19] << "  "
                  << global_energy[20] << "  "
                  << global_energy[21] << "\n";
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
}

void Diagnostics::write_px_distribution(double time,
                                        const Species& sp,
                                        int mpi_rank, int mpi_size)
{
    int ng = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;

    std::vector<double> local_Fv(Param::Nv, 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double v = sp.vgrid.v(iv);
            double sum = 0.0;
            size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                sum += sp.f[row + imu];
            }
            local_Fv[iv] += sum * 2.0 * Const::pi * v * v
                          * sp.vgrid.dmu * sp.sgrid->dx;
        }
    }

    std::vector<double> global_Fv(Param::Nv, 0.0);
    MPI_Reduce(local_Fv.data(), global_Fv.data(), Param::Nv,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fv_" << sp.name << "_"
              << std::setw(5) << std::setfill('0') << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# v[m/s]  F(v)\n";
        out << std::scientific << std::setprecision(8);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            out << sp.vgrid.v(iv) << "  " << global_Fv[iv] << "\n";
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
    out << "# x[um]  v[m/s]  mu  f_e[s^3/m^6]\n";
    out << std::scientific << std::setprecision(8);

    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int ix_g = ix + ng;
        const double x_um = sg.x(ix_g) / Const::micro;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double v = electrons.vgrid.v(iv);
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                out << x_um << "  "
                    << v << "  "
                    << electrons.vgrid.mu(imu) << "  "
                    << electrons.f[row + imu] << "\n";
            }
        }
    }
}

void Diagnostics::advance_snapshot()
{
    ++snapshot_count;
}
