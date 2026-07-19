#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <cstddef>

#ifndef FP_ENABLE_DEBUG_DIAGNOSTICS
#define FP_ENABLE_DEBUG_DIAGNOSTICS 1
#endif

// Test-only overrides for the production section-7.2 velocity-grid
// convergence audit.  Normal solver builds do not define these macros and
// retain the established 96 x 64 momentum grid.
#ifndef FP_VELOCITY_GRID_NV
#define FP_VELOCITY_GRID_NV 96
#endif

#ifndef FP_VELOCITY_GRID_NMU
#define FP_VELOCITY_GRID_NMU 64
#endif

// Deprecated compatibility macro.  The cylindrical grid no longer has a
// piecewise "refined cells" region; all resolutions use one fixed smooth
// mapping and differ only in cell count.
#ifndef FP_VELOCITY_GRID_REFINED_CELLS
#define FP_VELOCITY_GRID_REFINED_CELLS 32
#endif

#ifndef FP_VELOCITY_GRID_UPAR_STRETCH
#define FP_VELOCITY_GRID_UPAR_STRETCH 6.2
#endif

#ifndef FP_VELOCITY_GRID_UPERP_STRETCH
#define FP_VELOCITY_GRID_UPERP_STRETCH 5.9
#endif

// Section-11 operator tests use small spatial-grid convergence builds.  The
// production build keeps this unset and therefore retains nx=4000 exactly.
#ifndef FP_SPATIAL_GRID_NX
#define FP_SPATIAL_GRID_NX 4000
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

    const double Lx    = 8.0 * Const::micro;
    const int    nx    = FP_SPATIAL_GRID_NX;
#if FP_SPATIAL_GRID_NX == 4000
    // Preserve the bit pattern used by production checkpoints created before
    // FP_SPATIAL_GRID_NX became configurable.  Although Lx / 4000 is
    // mathematically identical, it differs by one binary64 ULP and changes
    // dx-dependent Beam weights and checkpoint identity hashes.
    const double dx    = 0.002 * Const::micro;
#else
    // Operator convergence builds keep the same physical domain while
    // varying only the number of spatial cells.
    const double dx    = Lx / static_cast<double>(nx);
#endif
    const double plasma_length = Lx;

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

    const double omega_pe = std::sqrt(dens * Const::qe * Const::qe /
                                      (Const::eps0 * Const::me));

    // Axisymmetric spherical momentum grid: (u, mu), u = p / (m c).
    // The distribution is normalized with d3u = 2*pi*u^2 du dmu.
    const int Nv  = FP_VELOCITY_GRID_NV;
    const int Nmu = FP_VELOCITY_GRID_NMU;
    const size_t Nvmu = static_cast<size_t>(Nv) * Nmu;

    const double Nsigma = 80.0;
    const double momentum_umax = 10.0;
    const double momentum_refined_u = 0.2;
    // Deprecated: retained only so external builds that reference the old
    // name still compile.  It does not control CylindricalVelocityGrid.
    const int    momentum_refined_cells = FP_VELOCITY_GRID_REFINED_CELLS;
    const double momentum_upar_stretch = FP_VELOCITY_GRID_UPAR_STRETCH;
    const double momentum_uperp_stretch = FP_VELOCITY_GRID_UPERP_STRETCH;
    const double diagnostic_tail_u_min = 8.0;
    const double vmax_fraction_c = 0.995;
    const int Nghost = 3;

    const double v_floor = 1.0e-12 * Const::c;
    const double u_floor = 1.0e-12;
}

#endif
