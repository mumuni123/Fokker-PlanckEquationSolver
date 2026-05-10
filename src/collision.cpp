#include "collision.h"
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <vector>

double CollisionOperator::chandrasekhar_G(double x)
{
    if (std::abs(x) < 1.0e-6) {
        return 2.0 / (3.0 * std::sqrt(Const::pi));
    }
    double erfx = std::erf(x);
    double erfpx = (2.0 / std::sqrt(Const::pi)) * std::exp(-x * x);
    return (erfx - x * erfpx) / (2.0 * x * x);
}

double CollisionOperator::coulomb_log(double n_e, double T_e) const
{
    double lambda_D = std::sqrt(Const::eps0 * T_e / (n_e * Const::qe * Const::qe));
    double arg = 12.0 * Const::pi * n_e * lambda_D * lambda_D * lambda_D;
    double lnL = std::log(arg);
    return std::max(2.0, lnL);
}

CollisionOperator::CollisionRates
CollisionOperator::compute_rates(double v, double mass_test,
                                 double n_field, double T_field,
                                 double m_field, double lnLambda,
                                 double Z_test, double Z_field) const
{
    CollisionRates r;
    if (v < 1.0e-20) {
        r.nu_s = 0.0;
        r.nu_parallel = 0.0;
        r.nu_perp = 0.0;
        return r;
    }

    double v_th = std::sqrt(T_field / m_field);
    if (v_th < 1.0e-20) {
        r.nu_s = 0.0;
        r.nu_parallel = 0.0;
        r.nu_perp = 0.0;
        return r;
    }

    double x = v / (std::sqrt(2.0) * v_th);
    double qe4 = Const::qe * Const::qe * Const::qe * Const::qe;
    double zfac = Z_test * Z_test * Z_field * Z_field;
    double nu0 = n_field * zfac * qe4 * lnLambda
               / (4.0 * Const::pi * Const::eps0 * Const::eps0
                  * mass_test * mass_test * v * v * v);

    double Gx = chandrasekhar_G(x);
    double psi = 2.0 * x * x * Gx;
    double erfx = std::erf(x);

    r.nu_s = (1.0 + mass_test / m_field) * psi * nu0;
    if (x > 1.0e-6) {
        r.nu_parallel = psi / (x * x) * nu0;
        r.nu_perp = (erfx - psi) / (x * x) * nu0;
    } else {
        r.nu_parallel = (4.0 / (3.0 * std::sqrt(Const::pi))) * nu0;
        r.nu_perp = r.nu_parallel;
    }
    return r;
}

void CollisionOperator::apply(Species& sp, double dt,
                              double n_field, double T_field, double m_field,
                              double Z_field, double Z_test)
{
    double lnL = coulomb_log(n_field, T_field);
    int ng = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;

    std::vector<double> nu_perp(Param::Nv, 0.0);
    std::vector<double> v2_center(Param::Nv, 0.0);
    std::vector<double> flux_A_face(Param::Nv + 1, 0.0);
    std::vector<double> flux_D_face(Param::Nv + 1, 0.0);

    for (int iv = 0; iv < Param::Nv; ++iv) {
        double v = std::max(sp.vgrid.v(iv), Param::v_floor);
        CollisionRates rates = compute_rates(v, sp.mass, n_field, T_field,
                                             m_field, lnL, Z_test, Z_field);
        nu_perp[iv] = rates.nu_perp;
        v2_center[iv] = v * v;
    }

    for (int ivf = 1; ivf < Param::Nv; ++ivf) {
        double vf = std::max(sp.vgrid.v_face(ivf), Param::v_floor);
        CollisionRates rates = compute_rates(vf, sp.mass, n_field, T_field,
                                             m_field, lnL, Z_test, Z_field);
        double A = -rates.nu_s * vf;
        double D = rates.nu_parallel * vf * vf;
        flux_A_face[ivf] = vf * vf * A;
        flux_D_face[ivf] = vf * vf * D;
    }

    // Axisymmetric spherical form of the 3V Fokker-Planck collision operator:
    //
    // C[f] = -1/v^2 d/dv (v^2 A f)
    //        + 1/(2 v^2) d/dv (v^2 D_parallel df/dv)
    //        + D_perp/(2 v^2) d/dmu ((1-mu^2) df/dmu)
    //
    // with A = <Delta v_parallel>/dt = -nu_s v,
    // D_parallel = <Delta v_parallel^2>/dt = nu_parallel v^2,
    // D_perp = <Delta v_perp^2>/dt = nu_perp v^2.
    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            double v2 = v2_center[iv];

            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double f0 = sp.f[idx3(ix_g, iv, imu)];

                double flux_v_l = 0.0;
                double flux_v_r = 0.0;

                if (iv > 0) {
                    double fL = sp.f[idx3(ix_g, iv - 1, imu)];
                    double fR = f0;
                    double f_up = (flux_A_face[iv] >= 0.0) ? fL : fR;
                    flux_v_l = flux_A_face[iv] * f_up
                             - 0.5 * flux_D_face[iv] * (fR - fL) / sp.vgrid.dv;
                }

                if (iv + 1 < Param::Nv) {
                    double fL = f0;
                    double fR = sp.f[idx3(ix_g, iv + 1, imu)];
                    double f_up = (flux_A_face[iv + 1] >= 0.0) ? fL : fR;
                    flux_v_r = flux_A_face[iv + 1] * f_up
                             - 0.5 * flux_D_face[iv + 1] * (fR - fL) / sp.vgrid.dv;
                }

                double radial = -(flux_v_r - flux_v_l) / (v2 * sp.vgrid.dv);

                double flux_mu_l = 0.0;
                double flux_mu_r = 0.0;

                if (imu > 0) {
                    double muf = sp.vgrid.mu_face(imu);
                    double fL = sp.f[idx3(ix_g, iv, imu - 1)];
                    double fR = f0;
                    flux_mu_l = (1.0 - muf * muf) * (fR - fL) / sp.vgrid.dmu;
                }
                if (imu + 1 < Param::Nmu) {
                    double muf = sp.vgrid.mu_face(imu + 1);
                    double fL = f0;
                    double fR = sp.f[idx3(ix_g, iv, imu + 1)];
                    flux_mu_r = (1.0 - muf * muf) * (fR - fL) / sp.vgrid.dmu;
                }

                double pitch = 0.5 * nu_perp[iv]
                             * (flux_mu_r - flux_mu_l) / sp.vgrid.dmu;

                double val = f0 + dt * (radial + pitch);
                sp.f_tmp[idx3(ix_g, iv, imu)] = std::max(0.0, val);
            }
        }
    }

    sp.f.swap(sp.f_tmp);
}
