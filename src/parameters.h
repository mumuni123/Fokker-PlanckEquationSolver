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
    enum PoissonSolverKind {
        POISSON_DISTRIBUTED_TRIDIAGONAL = 0,
        POISSON_PARALLEL_CYCLIC_REDUCTION = 1,
        POISSON_MULTIGRID = 2
    };

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
    const double Lx    = 8.0   * Const::micro;
    const int    nx    = static_cast<int>(Lx / dx);

    const double t_end         = 6.0 * Const::femto;
    const double t_inject_end  = 25.0  * Const::femto;
    const double dt_multiplier = 0.5;
    const double dt_snapshot   = 0.6 * Const::femto;
    const bool   enable_debug_diagnostics = false;
    const bool   enable_full_fe_output = false;
    const double velocity_space_cfl = 0.45;
    const double semi_lagrangian_cfl = 2.5;
    const PoissonSolverKind poisson_solver = POISSON_DISTRIBUTED_TRIDIAGONAL;
    const int    poisson_multigrid_vcycles = 10;
    const int    poisson_multigrid_pre_smooth = 3;
    const int    poisson_multigrid_post_smooth = 3;
    const int    beam_macro_particles_per_cell = 1000;
    const double beam_macro_weight = densb * dx / beam_macro_particles_per_cell;

    // Axisymmetric spherical velocity grid: (v, mu), mu = cos(theta) = vx / |v|.
    // This is still a 3D velocity-space model after integrating over the azimuth.
    const int Nv  = 96;
    const int Nmu = 64;
    const size_t Nvmu = static_cast<size_t>(Nv) * Nmu;

    const double Nsigma = 12.0;
    const int Nghost = 3;

    const double v_floor = 1.0e-12 * Const::c;
}

#endif
