#include "species.h"
#include <algorithm>
#include <cmath>
#include <omp.h>

namespace {
double maxwellian_raw_at(const Species& sp, int iv, int imu,
                         double drift_vx, double inv2vth2)
{
    const double vv = sp.vgrid.v_cells[iv];
    if (std::fabs(drift_vx) == 0.0) {
        return std::exp(-vv * vv * inv2vth2);
    }

    const double mu = sp.vgrid.mu_cells[imu];
    const double vx = vv * mu - drift_vx;
    const double vperp2 = vv * vv * (1.0 - mu * mu);
    return std::exp(-(vx * vx + vperp2) * inv2vth2);
}

double discrete_maxwellian_sum(const Species& sp,
                               double drift_vx,
                               double inv2vth2)
{
    double sum = 0.0;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        const double shell = sp.vgrid.moment_weight[iv];
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            sum += maxwellian_raw_at(sp, iv, imu, drift_vx, inv2vth2) * shell;
        }
    }
    return sum;
}

void fill_maxwellian_velocity_slice(const Species& sp,
                                    double drift_vx,
                                    std::vector<double>& out)
{
    out.assign(Param::Nvmu, 0.0);
    if (sp.type == SpeciesType::BEAM) return;

    const double inv2vth2 = sp.mass / (2.0 * sp.temperature);
    const double raw_sum = discrete_maxwellian_sum(sp, drift_vx, inv2vth2);
    const double norm = (raw_sum > 0.0) ? sp.density0 / raw_sum : 0.0;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        const size_t row = static_cast<size_t>(iv) * Param::Nmu;
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            out[row + imu] =
                norm * maxwellian_raw_at(sp, iv, imu, drift_vx, inv2vth2);
        }
    }
}
}

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
        vmax = std::min(Param::Nsigma * vth, Param::vmax_fraction_c * Const::c);
    }
    vgrid.init(vmax);
    set_maxwellian_boundary_drifts(0.0, 0.0);

    f.assign(local_size(), 0.0);
    f_tmp.assign(local_size(), 0.0);
    number_density.assign(sgrid->nx_local, 0.0);
    charge_density.assign(sgrid->nx_local, 0.0);
    current_x.assign(sgrid->nx_local, 0.0);
}

void Species::set_maxwellian_boundary_drifts(double left_drift_vx,
                                             double right_drift_vx)
{
    fill_maxwellian_velocity_slice(*this, left_drift_vx,
                                   maxwellian_left_boundary);
    fill_maxwellian_velocity_slice(*this, right_drift_vx,
                                   maxwellian_right_boundary);
}

void Species::initialize_maxwellian(double drift_vx)
{
    if (type == SpeciesType::BEAM) return;

    const double inv2vth2 = mass / (2.0 * temperature);
    const double raw_sum = discrete_maxwellian_sum(*this, drift_vx, inv2vth2);
    const double norm = (raw_sum > 0.0) ? density0 / raw_sum : 0.0;
    const int nxt = sgrid->nx_total;
    const bool zero_drift = std::fabs(drift_vx) == 0.0;

    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxt; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const size_t row = static_cast<size_t>(ix) * Param::Nvmu
                             + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double raw = zero_drift
                    ? maxwellian_raw_at(*this, iv, 0, 0.0, inv2vth2)
                    : maxwellian_raw_at(*this, iv, imu, drift_vx, inv2vth2);
                f[row + imu] = norm * raw;
            }
        }
    }
}

void Species::apply_vx_shift_profile(const std::vector<double>& delta_vx)
{
    if (type == SpeciesType::BEAM) return;
    if (delta_vx.empty()) return;

    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    if (delta_vx.size() < static_cast<size_t>(nxl)) return;

    #pragma omp parallel for schedule(static)
    for (int ix = 0; ix < nxl; ++ix) {
        const double shift = delta_vx[static_cast<size_t>(ix)];
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;

        if (std::fabs(shift) < 1.0e-30) {
            std::copy(f.begin() + xbase, f.begin() + xbase + Param::Nvmu,
                      f_tmp.begin() + xbase);
            continue;
        }

        std::fill(f_tmp.begin() + xbase,
                  f_tmp.begin() + xbase + Param::Nvmu,
                  0.0);

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double v = vgrid.v_cells[iv];
            const double v2 = v * v;
            const double shell = vgrid.moment_weight[iv];
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t offset = row + imu;
                const double cell_mass = std::max(0.0, f[offset]) * shell;
                if (cell_mass == 0.0) continue;

                const double mu = vgrid.mu_cells[imu];
                const double vx_new = v * mu + shift;
                const double vperp2 = std::max(0.0, v2 * (1.0 - mu * mu));
                const double v_new = std::sqrt(vx_new * vx_new + vperp2);
                if (v_new >= vgrid.v_max || v_new < vgrid.v_min) continue;

                double mu_new = 0.0;
                if (v_new > Param::v_floor) {
                    mu_new = std::max(-1.0, std::min(1.0, vx_new / v_new));
                }

                double sv = (v_new - vgrid.v_min) * vgrid.inv_dv - 0.5;
                int iv0 = static_cast<int>(std::floor(sv));
                double wv1 = sv - iv0;
                if (iv0 < 0) {
                    iv0 = 0;
                    wv1 = 0.0;
                } else if (iv0 >= Param::Nv - 1) {
                    iv0 = Param::Nv - 1;
                    wv1 = 0.0;
                }
                const int iv1 = (iv0 + 1 < Param::Nv) ? iv0 + 1 : iv0;
                const double wv0 = 1.0 - wv1;

                double smu = (mu_new - vgrid.mu_min) * vgrid.inv_dmu - 0.5;
                int imu0 = static_cast<int>(std::floor(smu));
                double wm1 = smu - imu0;
                if (imu0 < 0) {
                    imu0 = 0;
                    wm1 = 0.0;
                } else if (imu0 >= Param::Nmu - 1) {
                    imu0 = Param::Nmu - 1;
                    wm1 = 0.0;
                }
                const int imu1 = (imu0 + 1 < Param::Nmu) ? imu0 + 1 : imu0;
                const double wm0 = 1.0 - wm1;

                f_tmp[xbase + static_cast<size_t>(iv0) * Param::Nmu + imu0]
                    += cell_mass * wv0 * wm0;
                if (imu1 != imu0) {
                    f_tmp[xbase + static_cast<size_t>(iv0) * Param::Nmu + imu1]
                        += cell_mass * wv0 * wm1;
                }
                if (iv1 != iv0) {
                    f_tmp[xbase + static_cast<size_t>(iv1) * Param::Nmu + imu0]
                        += cell_mass * wv1 * wm0;
                    if (imu1 != imu0) {
                        f_tmp[xbase + static_cast<size_t>(iv1) * Param::Nmu + imu1]
                            += cell_mass * wv1 * wm1;
                    }
                }
            }
        }

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double inv_shell = 1.0 / vgrid.moment_weight[iv];
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t offset = row + imu;
                f_tmp[offset] = std::max(0.0, f_tmp[offset] * inv_shell);
            }
        }
    }

    for (int ix = 0; ix < nxl; ++ix) {
        const size_t xbase = static_cast<size_t>(ix + ng) * Param::Nvmu;
        std::copy(f_tmp.begin() + xbase, f_tmp.begin() + xbase + Param::Nvmu,
                  f.begin() + xbase);
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
