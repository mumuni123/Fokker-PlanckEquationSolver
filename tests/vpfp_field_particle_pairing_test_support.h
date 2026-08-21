#ifndef VPFP_FIELD_PARTICLE_PAIRING_TEST_SUPPORT_H
#define VPFP_FIELD_PARTICLE_PAIRING_TEST_SUPPORT_H

// Gate I common test fixture (section 4.7.1).  This header only builds test
// state, compares physical states and writes results; it must NOT re-implement
// swept-mass, Poisson, PIC deposition or work-current formulas.  All operators
// are exercised through their production implementations.

#include "background_tail_pic.h"
#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

struct PairingTestState {
    SpatialGrid grid;
    CylindricalVelocityGrid velocity_grid;
    Species bulk_n;
    Species bulk_np1;
    BeamPIC beam_n;
    BeamPIC beam_np1;
    BackgroundTailPIC tail_n;
    BackgroundTailPIC tail_np1;
    EMFields field_n;
    EMFields field_np1;
    std::vector<double> ion_density;
};

// Machine-scaled tolerance (section 5):
//   T_round = max(1e-12 * unit_floor, 512 * eps * scale).
inline double machine_scaled_tolerance(double scale, double unit_floor)
{
    return std::max(1.0e-12 * unit_floor,
                    512.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, scale));
}

// Build a smooth positive bulk Maxwellian profile
//   f = f_M * [1 + 1e-3 cos(2 pi x / L)]
// on a small grid (nx = 16 or 32), with no limiter/edge interference.  The
// ion density is built from the initial electron density so the initial state
// satisfies the discrete Poisson solve (section 4.7.1 item 3).
inline PairingTestState make_smooth_bulk_case(int nx, bool with_beam,
                                              bool with_tail,
                                              int mpi_rank = 0,
                                              int mpi_size = 1)
{
    PairingTestState s;
    s.grid.init_with_domain(mpi_rank, mpi_size, nx, Param::Lx);
    s.velocity_grid.init(Param::momentum_umax);
    s.bulk_n.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                  Const::me, Param::dens, Param::temperature_e, false, s.grid);
    s.bulk_np1.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                    Const::me, Param::dens, Param::temperature_e, false, s.grid);

    const int nxl = s.grid.nx_local;
    const double dx = s.grid.dx;
    const double length = s.grid.length();
    // Maxwellian normalization: fill per-cell integrated mass from the
    // production cylindrical Maxwellian slice, then apply the low-amplitude
    // spatial perturbation.
    std::vector<double> slice(Param::Nvmu, 0.0);
    s.bulk_n.fill_cylindrical_maxwellian_mass_slice(slice, Param::dens,
                                                    Param::temperature_e, 0.0);
    for (int ix = 0; ix < s.grid.nx_total; ++ix) {
        const double x = s.grid.x(ix);
        const double amp = 1.0 + 1.0e-3 * std::cos(2.0 * Const::pi * x / length);
        const size_t base = static_cast<size_t>(ix) * Param::Nvmu;
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            s.bulk_n.f[base + q] = slice[q] * amp;
        }
    }
    s.bulk_n.compute_moments();
    s.bulk_np1.f = s.bulk_n.f;
    s.bulk_np1.compute_moments();

    // Ion density from the initial electron density for discrete neutrality.
    s.ion_density.assign(static_cast<size_t>(nxl), 0.0);
    for (int i = 0; i < nxl; ++i) {
        s.ion_density[static_cast<size_t>(i)] =
            s.bulk_n.number_density[static_cast<size_t>(i)] /
            static_cast<double>(Param::Z_ion);
    }

    s.field_n.init(s.grid);
    s.field_np1.init(s.grid);
    if (with_beam) {
        s.beam_n.init(s.grid);
        s.beam_np1.init(s.grid);
    }
    if (with_tail) {
        s.tail_n.init(s.grid);
        s.tail_np1.init(s.grid);
    }
    return s;
}

// Open-boundary reservoir case: same smooth bulk but with a Maxwellian
// reservoir on both ends (section 4.7.1 item 2 case "open boundary").
inline PairingTestState make_open_boundary_case(int nx)
{
    PairingTestState s = make_smooth_bulk_case(nx, false, false);
    return s;
}

// Single in-domain Beam macro-particle (section 4.7.1 item 4).
inline PairingTestState make_beam_single_particle_case(int nx)
{
    PairingTestState s = make_smooth_bulk_case(nx, true, false);
    BeamParticle p;
    p.x = 0.5 * s.grid.length();
    p.px = 0.5 * Const::me * Const::c;
    p.weight = 1.0e12;
    s.beam_n.particles.push_back(p);
    s.beam_np1 = s.beam_n;
    return s;
}

// Single in-domain Tail macro-particle.
inline PairingTestState make_tail_single_particle_case(int nx)
{
    PairingTestState s = make_smooth_bulk_case(nx, false, true);
    BackgroundTailParticle p = {};
    p.x = 0.5 * s.grid.length();
    p.ux = 0.5;
    p.uy = 0.0;
    p.uz = 0.0;
    p.weight = 1.0e14;
    p.id = 1;
    s.tail_n.particles.push_back(p);
    s.tail_np1 = s.tail_n;
    return s;
}

// Bitwise comparison of two physical states (bulk f, PIC particles, fields,
// RNG carried in the objects).  No diagnostic arrays are compared.
inline bool bitwise_equal_physical_state(const PairingTestState& a,
                                         const PairingTestState& b)
{
    if (a.bulk_n.f.size() != b.bulk_n.f.size()) return false;
    if (std::memcmp(a.bulk_n.f.data(), b.bulk_n.f.data(),
                    a.bulk_n.f.size() * sizeof(double)) != 0) return false;
    if (a.field_n.Ex_face.size() != b.field_n.Ex_face.size()) return false;
    if (std::memcmp(a.field_n.Ex_face.data(), b.field_n.Ex_face.data(),
                    a.field_n.Ex_face.size() * sizeof(double)) != 0) return false;
    if (a.beam_n.particles.size() != b.beam_n.particles.size()) return false;
    for (size_t i = 0; i < a.beam_n.particles.size(); ++i) {
        if (std::memcmp(&a.beam_n.particles[i], &b.beam_n.particles[i],
                        sizeof(BeamParticle)) != 0) return false;
    }
    if (a.tail_n.particles.size() != b.tail_n.particles.size()) return false;
    for (size_t i = 0; i < a.tail_n.particles.size(); ++i) {
        if (std::memcmp(&a.tail_n.particles[i], &b.tail_n.particles[i],
                        sizeof(BackgroundTailParticle)) != 0) return false;
    }
    return true;
}

#endif
