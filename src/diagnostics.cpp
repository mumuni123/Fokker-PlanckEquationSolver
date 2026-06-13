#include "diagnostics.h"
#include "beam_pic.h"
#include <algorithm>
#include <cmath>
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
      energy_reference(0.0),
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
            if (vx > 0.0) {
                values[0] += vx * electrons.f[xbase + k] * weight;
                values[4] += electrons.charge * vx * electrons.f[xbase + k] * weight;
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
            if (vx < 0.0) {
                values[2] += -vx * electrons.f[xbase + k] * weight;
                values[5] += electrons.charge * vx * electrons.f[xbase + k] * weight;
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
    const double dEx_dx =
        (fields.Ex[ix_g + 1] - fields.Ex[ix_g - 1]) / (2.0 * sg.dx);
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
            step_file << "# step  time[fs]  max_abs_Ex[V/m]  x_at_max_abs_Ex[m]  "
                      << "gauss_residual_int[m^-2]  gauss_residual_abs_max[m^-3]  "
                      << "gauss_residual_abs_max_over_n0  "
                      << "N_bkg_e  "
                      << "N_beam_macro  N_beam_weighted  "
                      << "beam_cont_l1  beam_cont_linf  "
                      << "nsub_u1  nsub_mu1  nsub_u2  nsub_mu2  "
                      << "loss_u1  loss_mu1  loss_u2  loss_mu2  "
                      << "loss_u1_low  loss_u1_high  "
                      << "loss_u2_low  loss_u2_high  "
                      << "net_Nb_change[m^-2]  "
                      << "KE_bkg_e[J/m2]  KE_beam[J/m2]  E_field[J/m2]  "
                      << "E_total[J/m2]  dKE_bkg[J/m2]  dKE_beam[J/m2]  "
                      << "dE_field[J/m2]  W_bkg_E[J/m2]  W_beam_E[J/m2]  "
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
                      << "integral_f_u_gt_8_x[m^-3]\n";
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
                                         double local_integral_f_u_gt_8_x)
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
    double local_beam_continuity[2] = {
        beam.last_continuity_l1_error(),
        beam.last_continuity_linf_error()
    };
    double global_beam_continuity_l1 = 0.0;
    double global_beam_continuity_linf = 0.0;

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
    double global_max_abs_Ex = 0.0;
    double global_x_at_max_abs_Ex = 0.0;
    double global_charge_residual_int = 0.0;
    double global_charge_residual_abs_max = 0.0;
    double global_N_bkg_e = 0.0;
    double global_N_beam_macro = 0.0;
    double global_N_beam_weighted = 0.0;
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
    MPI_Reduce(&local_charge_residual_int, &global_charge_residual_int, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_charge_residual_abs_max,
               &global_charge_residual_abs_max, 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_bkg_e, &global_N_bkg_e, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_macro, &global_N_beam_macro, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_N_beam_weighted, &global_N_beam_weighted, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_beam_continuity[0], &global_beam_continuity_l1,
               1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_beam_continuity[1], &global_beam_continuity_linf,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_nsub, global_nsub, 4, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_losses, global_losses, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_energy, global_energy, 22, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        const double current_ke_per_particle_eV =
            (global_N_bkg_e > 0.0)
            ? global_energy[0] / global_N_bkg_e / Const::eV
            : 0.0;
        const double effective_T_eV =
            (2.0 / 3.0) * current_ke_per_particle_eV;

        step_file << step << "  "
                  << time / Const::femto << "  "
                  << global_max_abs_Ex << "  "
                  << global_x_at_max_abs_Ex << "  "
                  << global_charge_residual_int << "  "
                  << global_charge_residual_abs_max << "  "
                  << global_charge_residual_abs_max / Param::dens << "  "
                  << global_N_bkg_e << "  "
                  << global_N_beam_macro << "  "
                  << global_N_beam_weighted << "  "
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
                  << global_energy[21] << "  "
                  << initial_ke_per_particle_eV << "  "
                  << effective_T_eV << "  "
                  << global_x_at_max_loss_u_high << "  "
                  << global_f_u_max_x << "  "
                  << global_integral_f_u_gt_8_x << "\n";
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
    out << "# x[um]  u[p/(m c)]  mu  f_e[u^-3 m^-3]\n";
    out << std::scientific << std::setprecision(8);

    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int ix_g = ix + ng;
        const double x_um = sg.x(ix_g) / Const::micro;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u = electrons.vgrid.v(iv);
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                out << x_um << "  "
                    << u << "  "
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
