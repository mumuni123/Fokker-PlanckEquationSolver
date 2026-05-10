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
        vmax = std::min(Param::Nsigma * vth, 0.95 * Const::c);
    }
    vgrid.init(vmax);

    f.assign(local_size(), 0.0);
    f_tmp.assign(local_size(), 0.0);
    number_density.assign(sgrid->nx_local, 0.0);
    charge_density.assign(sgrid->nx_local, 0.0);
    current_x.assign(sgrid->nx_local, 0.0);
}

void Species::initialize_maxwellian(double drift_vx)
{
    if (type == SpeciesType::BEAM) return;

    double norm = density0 * std::pow(2.0 * Const::pi * temperature / mass, -1.5);
    double inv2vth2 = mass / (2.0 * temperature);
    int nxt = sgrid->nx_total;

    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxt; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = vgrid.v(iv);
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = vgrid.mu(imu);
                double vx = vv * mu - drift_vx;
                double vperp2 = vv * vv * (1.0 - mu * mu);
                double v2 = vx * vx + vperp2;
                f[idx3(ix, iv, imu)] = norm * std::exp(-v2 * inv2vth2);
            }
        }
    }
}

void Species::compute_moments()
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double weight0 = 2.0 * Const::pi * vgrid.dv * vgrid.dmu;

    std::fill(number_density.begin(), number_density.end(), 0.0);
    std::fill(charge_density.begin(), charge_density.end(), 0.0);
    std::fill(current_x.begin(), current_x.end(), 0.0);

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double n = 0.0;
        double jx_over_q = 0.0;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = vgrid.v(iv);
            double shell = weight0 * vv * vv;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double mu = vgrid.mu(imu);
                double fval = f[idx3(ix_g, iv, imu)];
                double w = shell;
                n += fval * w;
                jx_over_q += fval * (vv * mu) * w;
            }
        }
        number_density[ix] = n;
        charge_density[ix] = charge * n;
        current_x[ix] = charge * jx_over_q;
    }
}

double Species::total_particle_number() const
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double weight0 = 2.0 * Const::pi * vgrid.dv * vgrid.dmu;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double n = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = vgrid.v(iv);
            double shell = weight0 * vv * vv;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                n += f[idx3(ix_g, iv, imu)] * shell;
            }
        }
        total += n * sgrid->dx;
    }
    return total;
}

double Species::total_kinetic_energy() const
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double weight0 = 2.0 * Const::pi * vgrid.dv * vgrid.dmu;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double e = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            double vv = vgrid.v(iv);
            double shell = weight0 * vv * vv;
            double ke = relativistic_push
                      ? (gamma_from_v(vv) - 1.0) * mass * Const::c * Const::c
                      : 0.5 * mass * vv * vv;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                e += f[idx3(ix_g, iv, imu)] * ke * shell;
            }
        }
        total += e * sgrid->dx;
    }
    return total;
}
