#include "species.h"
#include <cmath>
#include <omp.h>

Species::Species()
    : charge(0), mass(0), density0(0), temperature(0),
      collisions_enabled(true), sgrid(NULL)
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
    sgrid = &sg;

    // 根据物种类型设置动量网格范围
    if (type == SpeciesType::BEAM) {
        // 束流：以漂移动量为中心，px 方向非对称
        double pd = Param::beam_drift_px;
        double pmax_px = 3.0 * pd;
        pgrid.init(-0.5 * pd, pmax_px,
                   -pd, pd,
                   -pd, pd);
    } else if (type == SpeciesType::BACKGROUND_ELECTRON) {
        // 背景电子：以 0 为中心，宽度为 Nsigma * 热动量
        double pth = std::sqrt(mass * temperature);
        double pmax = Param::Nsigma * pth;
        pgrid.init(-pmax, pmax, -pmax, pmax, -pmax, pmax);
    } else {
        // 背景离子：以 0 为中心
        double pth = std::sqrt(mass * temperature);
        double pmax = Param::Nsigma * pth;
        pgrid.init(-pmax, pmax, -pmax, pmax, -pmax, pmax);
    }

    // 分配数组
    size_t sz = local_size();
    f.assign(sz, 0.0);

    int nxl = sgrid->nx_local;
    number_density.assign(nxl, 0.0);
    charge_density.assign(nxl, 0.0);
    current_x.assign(nxl, 0.0);
    current_y.assign(nxl, 0.0);
    current_z.assign(nxl, 0.0);
}

void Species::initialize_maxwellian(double drift_px)
{
    if (type == SpeciesType::BEAM) {
        // 束流初始为空，在仿真过程中注入
        return;
    }

    // 非相对论麦克斯韦分布：f(p) = n / (2*pi*m*kT)^{3/2} * exp(-p^2/(2*m*kT))
    double norm = density0 * std::pow(2.0 * Const::pi * mass * temperature, -1.5);
    double inv2mkT = 1.0 / (2.0 * mass * temperature);
    double dp3 = pgrid.dpx * pgrid.dpy * pgrid.dpz;

    int ng = sgrid->nghost;
    int nxt = sgrid->nx_total;

    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxt; ++ix) {
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double ppx = pgrid.px(ipx) - drift_px;
            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                double ppy = pgrid.py(ipy);
                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    double ppz = pgrid.pz(ipz);
                    double p2 = ppx * ppx + ppy * ppy + ppz * ppz;
                    double val = norm * std::exp(-p2 * inv2mkT);
                    f[idx4(ix, ipx, ipy, ipz)] = val;
                }
            }
        }
    }
}

void Species::compute_moments()
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double dp3 = pgrid.dpx * pgrid.dpy * pgrid.dpz;

    // 先清零
    for (int i = 0; i < nxl; ++i) {
        number_density[i] = 0.0;
        charge_density[i] = 0.0;
        current_x[i] = 0.0;
        current_y[i] = 0.0;
        current_z[i] = 0.0;
    }

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng; // 含幽灵偏移的索引
        double local_n = 0.0;
        double local_jx = 0.0;
        double local_jy = 0.0;
        double local_jz = 0.0;

        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double ppx = pgrid.px(ipx);
            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                double ppy = pgrid.py(ipy);
                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    double ppz = pgrid.pz(ipz);
                    double fval = f[idx4(ix_g, ipx, ipy, ipz)];

                    local_n += fval;

                    double gam = gamma_from_p(ppx, ppy, ppz, mass);
                    double inv_gm = 1.0 / (gam * mass);
                    local_jx += fval * ppx * inv_gm;
                    local_jy += fval * ppy * inv_gm;
                    local_jz += fval * ppz * inv_gm;
                }
            }
        }
        number_density[ix] = local_n * dp3;
        charge_density[ix] = charge * local_n * dp3;
        current_x[ix] = charge * local_jx * dp3;
        current_y[ix] = charge * local_jy * dp3;
        current_z[ix] = charge * local_jz * dp3;
    }
}

double Species::total_particle_number() const
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double dp3 = pgrid.dpx * pgrid.dpy * pgrid.dpz;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double local_n = 0.0;
        for (size_t ip = 0; ip < Param::Np3; ++ip) {
            local_n += f[static_cast<size_t>(ix_g) * Param::Np3 + ip];
        }
        total += local_n * dp3 * sgrid->dx;
    }
    return total;
}

double Species::total_kinetic_energy() const
{
    int ng = sgrid->nghost;
    int nxl = sgrid->nx_local;
    double dp3 = pgrid.dpx * pgrid.dpy * pgrid.dpz;
    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double local_e = 0.0;
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double ppx = pgrid.px(ipx);
            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                double ppy = pgrid.py(ipy);
                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    double ppz = pgrid.pz(ipz);
                    double fval = f[idx4(ix_g, ipx, ipy, ipz)];
                    // 相对论动能：(gamma - 1)*m*c^2
                    double gam = gamma_from_p(ppx, ppy, ppz, mass);
                    local_e += fval * (gam - 1.0) * mass * Const::c * Const::c;
                }
            }
        }
        total += local_e * dp3 * sgrid->dx;
    }
    return total;
}
