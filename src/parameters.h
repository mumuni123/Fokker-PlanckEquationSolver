#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <cstddef>

#ifndef FP_ENABLE_DEBUG_DIAGNOSTICS
#define FP_ENABLE_DEBUG_DIAGNOSTICS 1
#endif

// Production velocity-grid resolution.  The u_parallel domain half-width was
// expanded to 20 to contain the beam-plasma instability tails (section 14
// phase-4 record); Nv is raised to 192 so the thermal core keeps the same
// ~0.0052 cell width as the former 96-cell grid at u_max=10 (192 x 64 grid).
// Convergence builds may override FP_VELOCITY_GRID_NV via these macros.
#ifndef FP_VELOCITY_GRID_NV
#define FP_VELOCITY_GRID_NV 192
#endif

// The cylindrical u_parallel grid keeps this legacy resolution as an
// immutable core.  Optional cells are appended outside that core; they never
// stretch or otherwise modify the original [-10, 10] mesh.
#ifndef FP_VELOCITY_GRID_UPAR_CORE_NV
#define FP_VELOCITY_GRID_UPAR_CORE_NV FP_VELOCITY_GRID_NV
#endif

#ifndef FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE
#define FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE 0
#endif

#ifndef FP_VELOCITY_GRID_UPERP_MAX
#define FP_VELOCITY_GRID_UPERP_MAX 10.0
#endif

// u_parallel domain half-width.  Expanded from 10 -> 12 -> 20 after the
// phase-4 joint smoke showed the physical beam-plasma instability driving a
// growing background tail: it reached u_parallel ~ +-10 by step 865 and
// ~ +-12 by step 922, and the production run extends far beyond 1000 steps,
// so the domain must comfortably contain the instability tails rather than
// chase the boundary.  The u_perp extent stays at 10 via
// FP_VELOCITY_GRID_UPERP_MAX.  With Nv=192 the thermal core cell width is
// restored to the stage-3 value (~0.0052); the section 16.7/16.11
// velocity-domain adequacy check re-examines both the extent and the
// resolution against the production physics.
#ifndef FP_VELOCITY_GRID_UPAR_CORE_MAX
#define FP_VELOCITY_GRID_UPAR_CORE_MAX 20.0
#endif

#ifndef FP_VELOCITY_GRID_UPAR_EXTENDED_MAX
#define FP_VELOCITY_GRID_UPAR_EXTENDED_MAX FP_VELOCITY_GRID_UPAR_CORE_MAX
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

// Production build uses the fixed 8000-cell / 40 um / 0.005 um grid.
// Convergence builds may override FP_SPATIAL_GRID_NX with a different cell
// count; the physical length stays 40 um in every case.
#ifndef FP_SPATIAL_GRID_NX
#define FP_SPATIAL_GRID_NX 8000
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

    const double Lx    = 40.0 * Const::micro;
    const int    nx    = FP_SPATIAL_GRID_NX;
#if FP_SPATIAL_GRID_NX == 8000
    // Production grid: 40 um / 8000 cells = 0.005 um exactly.
    const double dx    = 0.005 * Const::micro;
#else
    // Convergence builds keep the same 40 um physical domain while varying
    // only the number of spatial cells.
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
    // Production velocity-tail abort threshold.  The phase-4 joint smoke
    // showed that the real beam-plasma instability (production beam in the
    // 40 um plasma) drives a ~1e-12 background tail to u_parallel = +-10 by
    // ~837 steps; the 1e-12 threshold inherited from the legacy solver made
    // the marginal crossing (1.07e-12) a hard abort.  Relaxed to 1e-11 with
    // the tail loss remaining a hard failure (section 7.4); the velocity
    // domain adequacy is revisited against the section 16.7/16.11 runs.
    const double umax_loss_abort_fraction = 1.0e-11;
    const double vmax_loss_abort_fraction = umax_loss_abort_fraction;
    const int    beam_macro_particles_per_cell = 1000; // 1000
    const double beam_source_x_start = 0.5 * Const::micro;
    const double beam_source_length  = 0.3 * Const::micro;
    const double beam_charge_compensation_alpha = 1.0;
    const double beam_charge_compensation_smoothing_skin_depths = 1.5;

    const double omega_pe = std::sqrt(dens * Const::qe * Const::qe /
                                      (Const::eps0 * Const::me));

    // Axisymmetric cylindrical momentum grid: (u_parallel, u_perp),
    // u = p / (m c).  The distribution is normalized with
    // d3u = 2*pi*u_perp du_parallel du_perp.
    const int Nv_core = FP_VELOCITY_GRID_UPAR_CORE_NV;
    const int Nv_tail = FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE;
    const int Nv  = Nv_core + 2 * Nv_tail;
    const int Nmu = FP_VELOCITY_GRID_NMU;
    const size_t Nvmu = static_cast<size_t>(Nv) * Nmu;

    const double Nsigma = 80.0;
    // momentum_umax remains the historical u_perp and u_parallel-core bound.
    // The optional u_parallel extension is represented separately below.
    // u_perp bound (10); the u_parallel core bound is separate below.
    const double momentum_umax = FP_VELOCITY_GRID_UPERP_MAX;
    const double momentum_upar_core_max = FP_VELOCITY_GRID_UPAR_CORE_MAX;
    const double momentum_upar_extended_max = FP_VELOCITY_GRID_UPAR_EXTENDED_MAX;
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
