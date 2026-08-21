#include "species.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <omp.h>

namespace {
double maxwellian_raw_at(const Species& sp, int iv, int imu,
                         double drift_vx, double inv2uth2)
{
    const double u = sp.vgrid.v_cells[iv];
    if (std::fabs(drift_vx) == 0.0) {
        return std::exp(-u * u * inv2uth2);
    }

    const double mu = sp.vgrid.mu_cells[imu];
    const double drift_u = u_from_v(drift_vx);
    const double ux = u * mu - drift_u;
    const double uperp2 = u * u * (1.0 - mu * mu);
    return std::exp(-(ux * ux + uperp2) * inv2uth2);
}

double discrete_maxwellian_sum(const Species& sp,
                               double drift_vx,
                               double inv2uth2)
{
    double sum = 0.0;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        const double shell = sp.vgrid.moment_weight[iv];
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            sum += maxwellian_raw_at(sp, iv, imu, drift_vx, inv2uth2) * shell;
        }
    }
    return sum;
}

}

Species::Species()
    : charge(0.0), mass(0.0), density0(0.0), temperature(0.0),
      collisions_enabled(true), relativistic_push(false), sgrid(NULL),
      cylindrical_mass_representation(false)
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

    double umax = Param::momentum_umax;
    vgrid.init(umax);
    cylindrical_mass_representation =
        (type == SpeciesType::BACKGROUND_ELECTRON);
    if (cylindrical_mass_representation) cgrid.init(umax);

    f.assign(local_size(), 0.0);
    f_tmp.assign(local_size(), 0.0);
    number_density.assign(sgrid->nx_local, 0.0);
    charge_density.assign(sgrid->nx_local, 0.0);
    current_x.assign(sgrid->nx_local, 0.0);
    current_face_x.assign(sgrid->nx_local + 1, 0.0);
}

void Species::initialize_maxwellian(double drift_vx)
{
    if (type == SpeciesType::BEAM) return;

    if (cylindrical_mass_representation) {
        const double inv2uth2 = mass * Const::c * Const::c /
                                (2.0 * temperature);
        const double drift_u = u_from_v(drift_vx);
        double raw_number = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double up = cgrid.upar_cells[iv];
                const double ut = cgrid.uperp_cells[imu];
                raw_number += std::exp(-((up - drift_u) * (up - drift_u) +
                                         ut * ut) * inv2uth2) *
                              cgrid.cell_phase_volume(iv, imu);
            }
        }
        const double norm = (raw_number > 0.0) ? density0 / raw_number : 0.0;
        const int nxt = sgrid->nx_total;
        #pragma omp parallel for collapse(2)
        for (int ix = 0; ix < nxt; ++ix) {
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double up = cgrid.upar_cells[iv];
                    const double ut = cgrid.uperp_cells[imu];
                    const double f3 = norm *
                        std::exp(-((up - drift_u) * (up - drift_u) +
                                   ut * ut) * inv2uth2);
                    f[idx3(ix, iv, imu)] = f3 * sgrid->dx *
                        cgrid.cell_phase_volume(iv, imu);
                }
            }
        }
        f_tmp = f;
        return;
    }

    const double inv2uth2 = mass * Const::c * Const::c / (2.0 * temperature);
    const double raw_sum = discrete_maxwellian_sum(*this, drift_vx, inv2uth2);
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
                    ? maxwellian_raw_at(*this, iv, 0, 0.0, inv2uth2)
                    : maxwellian_raw_at(*this, iv, imu, drift_vx, inv2uth2);
                f[row + imu] = norm * raw;
            }
        }
    }
}

void Species::initialize_maxwellian_profile(const std::vector<double>& density_profile,
                                            double drift_vx)
{
    if (type == SpeciesType::BEAM) return;

    if (cylindrical_mass_representation) {
        const int ng = sgrid->nghost;
        const int nxl = sgrid->nx_local;
        if (density_profile.size() < static_cast<size_t>(nxl)) return;
        const double inv2uth2 = mass * Const::c * Const::c /
                                (2.0 * temperature);
        const double drift_u = u_from_v(drift_vx);
        double raw_number = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv)
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double up = cgrid.upar_cells[iv];
                const double ut = cgrid.uperp_cells[imu];
                raw_number += std::exp(-((up - drift_u) * (up - drift_u) +
                                         ut * ut) * inv2uth2) *
                              cgrid.cell_phase_volume(iv, imu);
            }
        #pragma omp parallel for
        for (int ix = 0; ix < nxl; ++ix) {
            const double norm = (raw_number > 0.0)
                ? std::max(0.0, density_profile[ix]) / raw_number : 0.0;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double up = cgrid.upar_cells[iv];
                    const double ut = cgrid.uperp_cells[imu];
                    f[idx3(ng + ix, iv, imu)] = norm *
                        std::exp(-((up - drift_u) * (up - drift_u) +
                                   ut * ut) * inv2uth2) * sgrid->dx *
                        cgrid.cell_phase_volume(iv, imu);
                }
            }
        }
        f_tmp = f;
        return;
    }

    std::fill(f.begin(), f.end(), 0.0);
    const double inv2uth2 = mass * Const::c * Const::c / (2.0 * temperature);
    const double raw_sum = discrete_maxwellian_sum(*this, drift_vx, inv2uth2);
    if (raw_sum <= 0.0) return;

    const bool zero_drift = std::fabs(drift_vx) == 0.0;
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    if (density_profile.size() < static_cast<size_t>(nxl)) return;

    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double density = std::max(0.0, density_profile[ix]);
            const double norm = density / raw_sum;
            const int ix_g = ix + ng;
            const size_t row = static_cast<size_t>(ix_g) * Param::Nvmu
                             + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double raw = zero_drift
                    ? maxwellian_raw_at(*this, iv, 0, 0.0, inv2uth2)
                    : maxwellian_raw_at(*this, iv, imu, drift_vx, inv2uth2);
                f[row + imu] = norm * raw;
            }
        }
    }
}

double Species::maxwellian_f_value(double density,
                                   double temp,
                                   double drift_vx,
                                   int iv,
                                   int imu) const
{
    if (type == SpeciesType::BEAM || iv < 0 || iv >= Param::Nv ||
        imu < 0 || imu >= Param::Nmu || !(temp > 0.0)) {
        return 0.0;
    }
    const double inv2uth2 = mass * Const::c * Const::c / (2.0 * temp);
    const double raw_sum = discrete_maxwellian_sum(*this, drift_vx, inv2uth2);
    if (!(raw_sum > 0.0)) return 0.0;
    return std::max(0.0, density) *
           maxwellian_raw_at(*this, iv, imu, drift_vx, inv2uth2) / raw_sum;
}

void Species::fill_maxwellian_velocity_slice(std::vector<double>& values,
                                             double density,
                                             double temp,
                                             double drift_vx) const
{
    values.assign(Param::Nvmu, 0.0);
    if (type == SpeciesType::BEAM || !(temp > 0.0)) return;

    const double inv2uth2 = mass * Const::c * Const::c / (2.0 * temp);
    const double raw_sum = discrete_maxwellian_sum(*this, drift_vx, inv2uth2);
    if (!(raw_sum > 0.0)) return;

    const double norm = std::max(0.0, density) / raw_sum;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        const size_t row = static_cast<size_t>(iv) * Param::Nmu;
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            values[row + static_cast<size_t>(imu)] =
                norm * maxwellian_raw_at(*this, iv, imu, drift_vx, inv2uth2);
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

    if (cylindrical_mass_representation) {
        #pragma omp parallel for
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            double n = 0.0;
            double gamma_x = 0.0;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double m = f[idx3(ix_g, iv, imu)];
                    n += m;
                    gamma_x += m * cgrid.vx[idx2(iv, imu)];
                }
            }
            number_density[ix] = n / sgrid->dx;
            charge_density[ix] = charge * number_density[ix];
            current_x[ix] = charge * gamma_x / sgrid->dx;
        }
        return;
    }

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

    if (cylindrical_mass_representation) {
        #pragma omp parallel for reduction(+:total)
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            for (int iv = 0; iv < Param::Nv; ++iv)
                for (int imu = 0; imu < Param::Nmu; ++imu)
                    total += f[idx3(ix_g, iv, imu)];
        }
        return total;
    }

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

    if (cylindrical_mass_representation) {
        #pragma omp parallel for reduction(+:total)
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            for (int iv = 0; iv < Param::Nv; ++iv)
                for (int imu = 0; imu < Param::Nmu; ++imu)
                    total += f[idx3(ix_g, iv, imu)] *
                             cgrid.kinetic_energy[idx2(iv, imu)];
        }
        return total;
    }

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double e = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double shell = vgrid.moment_weight[iv];
            const double ke =
                (vgrid.gamma_cells[iv] - 1.0) * mass * Const::c * Const::c;
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

    if (cylindrical_mass_representation) {
        #pragma omp parallel for reduction(+:total_n,total_e)
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            for (int iv = 0; iv < Param::Nv; ++iv) {
                for (int imu = 0; imu < Param::Nmu; ++imu) {
                    const double m = f[idx3(ix_g, iv, imu)];
                    total_n += m;
                    total_e += m * cgrid.kinetic_energy[idx2(iv, imu)];
                }
            }
        }
        number = total_n;
        kinetic_energy = total_e;
        return;
    }

    #pragma omp parallel for reduction(+:total_n,total_e)
    for (int ix = 0; ix < nxl; ++ix) {
        const int ix_g = ix + ng;
        const size_t xbase = static_cast<size_t>(ix_g) * Param::Nvmu;
        double n = 0.0;
        double e = 0.0;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double shell = vgrid.moment_weight[iv];
            const double ke =
                (vgrid.gamma_cells[iv] - 1.0) * mass * Const::c * Const::c;
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

void Species::fill_cylindrical_maxwellian_mass_slice(
    std::vector<double>& values, double density, double temp,
    double drift_vx) const
{
    values.assign(Param::Nvmu, 0.0);
    if (!cylindrical_mass_representation || !sgrid || !(temp > 0.0)) return;

    const double inv2uth2 = mass * Const::c * Const::c / (2.0 * temp);
    const double drift_u = u_from_v(drift_vx);
    double raw_number = 0.0;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double up = cgrid.upar_cells[iv];
            const double ut = cgrid.uperp_cells[imu];
            raw_number += std::exp(-((up - drift_u) * (up - drift_u) +
                                     ut * ut) * inv2uth2) *
                          cgrid.cell_phase_volume(iv, imu);
        }
    }
    if (!(raw_number > 0.0)) return;

    const double norm = std::max(0.0, density) / raw_number;
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double up = cgrid.upar_cells[iv];
            const double ut = cgrid.uperp_cells[imu];
            const double f3 = norm * std::exp(-((up - drift_u) * (up - drift_u) +
                                                ut * ut) * inv2uth2);
            values[idx2(iv, imu)] = f3 * sgrid->dx *
                                    cgrid.cell_phase_volume(iv, imu);
        }
    }
}

double Species::distribution_value(int ix_with_ghost, int iv, int imu) const
{
    const double stored = f[idx3(ix_with_ghost, iv, imu)];
    if (!cylindrical_mass_representation) return stored;
    const double volume = sgrid->dx * cgrid.cell_phase_volume(iv, imu);
    return (volume > 0.0) ? stored / volume : 0.0;
}

void Species::swap_state(Species& other)
{
    f.swap(other.f);
    f_tmp.swap(other.f_tmp);
    number_density.swap(other.number_density);
    charge_density.swap(other.charge_density);
    current_x.swap(other.current_x);
    current_face_x.swap(other.current_face_x);
}

SpeciesConversionResult Species::extract_conversion_masses(
    const std::vector<ConversionMassRequest>& requests)
{
    SpeciesConversionResult r;
    if (!cylindrical_mass_representation || !sgrid) return r;
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;

    // Validate every request before touching f (section 14.3).
    for (size_t q = 0; q < requests.size(); ++q) {
        const ConversionMassRequest& req = requests[q];
        const int il = req.ix_global - sgrid->ix_start;
        if (il < 0 || il >= nxl || req.iv < 0 || req.iv >= Param::Nv ||
            req.imu < 0 || req.imu >= Param::Nmu ||
            !std::isfinite(req.mass) || req.mass < 0.0) {
            return r;
        }
        const double stored = f[idx3(ng + il, req.iv, req.imu)];
        if (!std::isfinite(stored) || stored < 0.0) return r;
        const double tolerance =
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, stored);
        if (req.mass > stored + tolerance) return r;
    }

    // Apply: subtract exactly the requested masses and accumulate the
    // removed moments with the shared production formula.
    r.valid = true;
    r.complete = true;
    for (size_t q = 0; q < requests.size(); ++q) {
        const ConversionMassRequest& req = requests[q];
        const int il = req.ix_global - sgrid->ix_start;
        const size_t slot = idx3(ng + il, req.iv, req.imu);
        f[slot] -= req.mass;
        double number, px, energy, jx_dx, pixx_dx, piperp_dx;
        mass_cell_moments(req.mass, cgrid.upar_cells[req.iv],
                          cgrid.uperp_cells[req.imu], number, px, energy,
                          jx_dx, pixx_dx, piperp_dx);
        r.number_removed += number;
        r.px_removed += px;
        r.energy_removed += energy;
        r.jx_dx_removed += jx_dx;
        r.pixx_dx_removed += pixx_dx;
        r.piperp_dx_removed += piperp_dx;
        r.number_scale += std::fabs(number);
        r.px_scale += std::fabs(px);
        r.energy_scale += std::fabs(energy);
        r.jx_scale += std::fabs(jx_dx);
        r.pixx_scale += std::fabs(pixx_dx);
        r.piperp_scale += std::fabs(piperp_dx);
    }
    return r;
}

bool Species::validate_return_masses(
    const std::vector<ReturnMassRequest>& requests) const
{
    if (!cylindrical_mass_representation || !sgrid) return false;
    const int ng = sgrid->nghost;
    const int nxl = sgrid->nx_local;
    std::map<size_t, double> additions;
    for (size_t q = 0; q < requests.size(); ++q) {
        const ReturnMassRequest& req = requests[q];
        const int il = req.ix_global - sgrid->ix_start;
        if (il < 0 || il >= nxl || req.iv < 0 || req.iv >= Param::Nv ||
            req.imu < 0 || req.imu >= Param::Nmu ||
            !std::isfinite(req.mass) || req.mass < 0.0) return false;
        const size_t id = idx3(ng + il, req.iv, req.imu);
        const double next = additions[id] + req.mass;
        if (!std::isfinite(next)) return false;
        additions[id] = next;
    }
    for (std::map<size_t, double>::const_iterator it = additions.begin();
         it != additions.end(); ++it) {
        const double old = f[it->first];
        if (!std::isfinite(old) || !std::isfinite(old + it->second) ||
            old + it->second < 0.0) return false;
    }
    return true;
}

SpeciesReturnResult Species::add_return_masses(
    const std::vector<ReturnMassRequest>& requests)
{
    SpeciesReturnResult r;
    if (!validate_return_masses(requests)) return r;
    const int ng = sgrid->nghost;
    r.valid = true;
    r.complete = true;
    for (size_t q = 0; q < requests.size(); ++q) {
        const ReturnMassRequest& req = requests[q];
        const int il = req.ix_global - sgrid->ix_start;
        f[idx3(ng + il, req.iv, req.imu)] += req.mass;
        double n, px, e, jx, pixx, piperp;
        mass_cell_moments(req.mass, cgrid.upar_cells[req.iv],
                          cgrid.uperp_cells[req.imu], n, px, e, jx,
                          pixx, piperp);
        r.number_added += n; r.px_added += px; r.energy_added += e;
        r.jx_dx_added += jx; r.pixx_dx_added += pixx;
        r.piperp_dx_added += piperp;
        r.number_scale += std::fabs(n); r.px_scale += std::fabs(px);
        r.energy_scale += std::fabs(e); r.jx_scale += std::fabs(jx);
        r.pixx_scale += std::fabs(pixx); r.piperp_scale += std::fabs(piperp);
    }
    return r;
}
