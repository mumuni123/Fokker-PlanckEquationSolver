#include "species.h"
#include <algorithm>
#include <cmath>
#include <omp.h>

Species::Species()
    : charge(0.0), mass(0.0), density0(0.0), temperature(0.0),
      collisions_enabled(true), relativistic_push(false), sgrid(NULL)
{}

void Species::init(const std::string& n, SpeciesType t,
                   double q, double m,
                   double dens, double temp,
                   bool coll,
                   const SpatialGrid& sg)
{
    name = n;
    type = t;
    charge = q;
    mass = m;
    density0 = dens;
    temperature = temp;
    collisions_enabled = coll;
    relativistic_push = (type == SpeciesType::BEAM);
    sgrid = &sg;

    double vmax = 0.0;
    if (type == SpeciesType::BEAM) {
        vmax = 0.999999 * Const::c;
    } else {
        double vth = std::sqrt(temperature / mass);
        vmax = std::min(Param::Nsigma * vth, 0.98 * Const::c);
    }
    vgrid.init(vmax);
    precompute_maxwellian_boundary();

    f.assign(local_size(), 0.0);
    f_tmp.assign(local_size(), 0.0);
    number_density.assign(sgrid->nx_local, 0.0);
    charge_density.assign(sgrid->nx_local, 0.0);
    current_x.assign(sgrid->nx_local, 0.0);
}

void Species::precompute_maxwellian_boundary()
{
    maxwellian_boundary.assign(Param::Nv, 0.0);
    if (type == SpeciesType::BEAM) return;

    const double norm = density0 * std::pow(2.0 * Const::pi * temperature / mass, -1.5);
    const double inv2vth2 = mass / (2.0 * temperature);
    for (int iv = 0; iv < Param::Nv; ++iv) {
        const double vv = vgrid.v_cells[iv];
        maxwellian_boundary[iv] = norm * std::exp(-vv * vv * inv2vth2);
    }
}

void Species::initialize_maxwellian(double drift_vx)
{
    if (type == SpeciesType::BEAM) return;

    const double norm = density0 * std::pow(2.0 * Const::pi * temperature / mass, -1.5);
    const double inv2vth2 = mass / (2.0 * temperature);
    const int nxt = sgrid->nx_total;
    const bool zero_drift = std::fabs(drift_vx) == 0.0;

    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxt; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double vv = vgrid.v_cells[iv];
            const size_t row = static_cast<size_t>(ix) * Param::Nvmu
                             + static_cast<size_t>(iv) * Param::Nmu;
            if (zero_drift) {
                const double f0 = maxwellian_boundary[iv];
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    f[row + imu] = f0;
                }
                continue;
            }
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double mu = vgrid.mu_cells[imu];
                const double vx = vv * mu - drift_vx;
                const double vperp2 = vv * vv * (1.0 - mu * mu);
                const double v2 = vx * vx + vperp2;
                f[row + imu] = norm * std::exp(-v2 * inv2vth2);
            }
        }
    }
}

void Species::compute_moments()
{
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;

    std::fill(number_density.begin(), number_density.end(), 0.0);
    std::fill(charge_density.begin(), charge_density.end(), 0.0);
    std::fill(current_x.begin(), current_x.end(), 0.0);

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double n = 0.0;
        double jx_over_q = 0.0;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double w = vgrid.moment_weight[iv];
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            const size_t vrow = static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double fval = f[row + imu];
                n += fval * w;
                jx_over_q += fval * vgrid.current_weight[vrow + imu];
            }
        }
        number_density[ix] = n;
        charge_density[ix] = charge * n;
        current_x[ix] = charge * jx_over_q;
    }
}

double Species::total_particle_number() const
{
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double n = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double shell = vgrid.moment_weight[iv];
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                n += f[row + imu] * shell;
            }
        }
        total += n * sgrid->dx;
    }
    return total;
}

double Species::total_kinetic_energy() const
{
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double e = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double vv = vgrid.v_cells[iv];
            const double shell = vgrid.moment_weight[iv];
            const double ke = relativistic_push
                      ? (gamma_from_v(vv) - 1.0) * mass * Const::c * Const::c
                      : 0.5 * mass * vv * vv;
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                e += f[row + imu] * ke * shell;
            }
        }
        total += e * sgrid->dx;
    }
    return total;
}

void Species::total_particle_number_and_energy(double& number,
                                               double& kinetic_energy) const
{
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    double total_n = 0.0;
    double total_e = 0.0;

    #pragma omp parallel for reduction(+:total_n,total_e)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double n = 0.0;
        double e = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double vv = vgrid.v_cells[iv];
            const double shell = vgrid.moment_weight[iv];
            const double ke = relativistic_push
                      ? (gamma_from_v(vv) - 1.0) * mass * Const::c * Const::c
                      : 0.5 * mass * vv * vv;
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double weighted_f = f[row + imu] * shell;
                n += weighted_f;
                e += weighted_f * ke;
            }
        }
        total_n += n * sgrid->dx;
        total_e += e * sgrid->dx;
    }

    number = total_n;
    kinetic_energy = total_e;
}
