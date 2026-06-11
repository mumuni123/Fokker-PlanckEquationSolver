#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <cstddef>

#ifndef FP_ENABLE_DEBUG_DIAGNOSTICS
#define FP_ENABLE_DEBUG_DIAGNOSTICS 1
#endif

namespace Const {
    const double me   = 9.10938e-31;
    const double qe   = 1.60218e-19;
    const double c    = 2.99792e8;
    const double eps0 = 8.85419e-12;
    const double mu0  = 1.25664e-6;
    const double kB   = 1.38065e-23;
    const double pi   = 3.14159265358979323846;
    const double eV   = 1.60218e-19;
    const double micro = 1.0e-6;
    const double femto = 1.0e-15;
}

namespace Param {
    const double temperature_e = 100.0 * Const::eV;
    const double temperature_i = 10.0  * Const::eV;
    const double dens          = 1.2e29;
    const int    Z_ion         = 2;
    const double mass_ion_me   = 49572.0;
    const double mass_ion      = mass_ion_me * Const::me;

    const double jb        = 5.0e16;
    const double gambetab  = 8.60;
    const double gamb      = std::sqrt(1.0 + gambetab * gambetab);
    const double betab     = gambetab / gamb;
    const double densb     = jb / (Const::qe * betab * Const::c);
    const double beam_v0   = betab * Const::c;
    const double beam_p0   = gambetab * Const::me * Const::c;

    const double dx    = 0.002 * Const::micro;
    const double Lx    = 8.0 * Const::micro;
    const double plasma_length = Lx;
    const int    nx    = static_cast<int>(std::lround(Lx / dx));

    const double t_end         = 120.0 * Const::femto;
    const double t_inject_start = 0.0   * Const::femto;
    const double t_inject_end  = 25.0  * Const::femto;
    const double dt_multiplier = 0.5;
    const double dt_snapshot   = 0.6 * Const::femto;
    const bool   enable_debug_diagnostics = false;
    const bool   enable_full_fe_output = false;
    const bool   enable_beam_boundary_injection = true;
    const double velocity_space_cfl = 0.35;
    const double semi_lagrangian_cfl = 2.5;
    const bool   abort_on_vmax_loss = true;
    const double umax_loss_abort_fraction = 1.0e-12;
    const double vmax_loss_abort_fraction = umax_loss_abort_fraction;

    const int    beam_macro_particles_per_cell = 1000; // 1000
    const double beam_macro_weight = densb * dx / beam_macro_particles_per_cell;
    const double beam_source_x_start = 0.5 * Const::micro;
    const double beam_source_length  = 0.3 * Const::micro;
    const double beam_charge_compensation_alpha = 1.0;
    const double beam_charge_compensation_smoothing_skin_depths = 1.5;

    const double background_reservoir_drift_speed = jb / (Const::qe * dens);
    const double background_reservoir_u_return = -jb / (Const::qe * dens);
    const bool   enable_dynamic_background_reservoir = true;
    const double background_reservoir_density_relaxation = 0.25;
    const double background_reservoir_min_density_factor = 0.0;
    const double background_reservoir_max_density_factor = 8.0;
    const double background_reservoir_injection_skin_depths = 1.0;
    const double background_reservoir_source_length = 0.3 * Const::micro;
    const double background_reservoir_edge_exclusion = 0.0;
    const double background_reservoir_source_drift_speed = 0.0;
    const double background_reservoir_source_eta = 0.08;
    const double background_reservoir_feedback_beta = 0.05;
    const double background_reservoir_feedback_tau = 2.0 * Const::femto;

    const bool   poisson_remove_global_mean_charge = true;
    const double omega_pe = std::sqrt(dens * Const::qe * Const::qe /
                                      (Const::eps0 * Const::me));
    const double upstream_left_return_current_fraction = 1.0;
    const double upstream_left_drift_speed =
        -upstream_left_return_current_fraction * jb / (Const::qe * dens);
    const double upstream_left_density_wake_fraction = densb / dens;
    const double upstream_left_temperature_wake_fraction = 0.05;
    const double upstream_left_wake_frequency = omega_pe;
    const double upstream_left_wake_phase = 0.0;
    const double upstream_phase_lock_gain = 0.03;
    const double upstream_phase_lock_max_step = 0.02;
    const double upstream_phase_lock_ex_weight = 0.25;
    const double upstream_flux_balance_tau = 2.0 * Const::femto;
    const double upstream_flux_balance_gain = 0.25;
    const double upstream_flux_balance_min_scale = 0.05;
    const double upstream_flux_balance_max_scale = 64.0;
    const double upstream_right_density_wake_fraction = 0.0;
    const double upstream_right_temperature_wake_fraction = 0.0;
    const double upstream_right_drift_speed = 0.0;

    // Axisymmetric spherical momentum grid: (u, mu), u = p / (m c).
    // The distribution is normalized with d3u = 2*pi*u^2 du dmu.
    const int Nv  = 96;
    const int Nmu = 64;
    const size_t Nvmu = static_cast<size_t>(Nv) * Nmu;

    const double Nsigma = 80.0;
    const double momentum_umax = 10.0;
    const double momentum_refined_u = 0.2;
    const int    momentum_refined_cells = 32;
    const double diagnostic_tail_u_min = 8.0;
    const double vmax_fraction_c = 0.995;
    const int Nghost = 3;

    const double v_floor = 1.0e-12 * Const::c;
    const double u_floor = 1.0e-12;
}

#endif
