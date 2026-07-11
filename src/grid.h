#ifndef GRID_H
#define GRID_H

#include "parameters.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

struct VelocityGrid {
    double v_min, v_max, dv;
    double mu_min, mu_max, dmu;
    double inv_dv, inv_dmu;
    double max_abs_mu;
    double max_mu_face_factor;
    double max_inv_v;
    double max_speed;

    std::vector<double> v_cells;
    std::vector<double> v_faces;
    std::vector<double> v_widths;
    std::vector<double> inv_v_widths;
    std::vector<double> inv_v_center_dist;
    std::vector<double> mu_cells;
    std::vector<double> mu_faces;
    std::vector<double> v2_cells;
    std::vector<double> inv_v_cells;
    std::vector<double> inv_v2_cells;
    std::vector<double> v2_faces;
    std::vector<double> gamma_cells;
    std::vector<double> beta_cells;
    std::vector<double> speed_cells;
    std::vector<double> chain_speed_cells;
    std::vector<double> chain_speed_faces;
    std::vector<double> mu_face_factor;
    std::vector<double> moment_weight;
    std::vector<double> inv_moment_weight;
    std::vector<double> mu_flux_scale;
    std::vector<double> mu_cfl_factor;
    std::vector<double> vx_cells;
    std::vector<double> current_weight;

    void init(double umax)
    {
        v_min = 0.0;
        v_max = umax;
        dv = (v_max - v_min) / Param::Nv;
        inv_dv = 1.0 / dv;
        mu_min = -1.0;
        mu_max = 1.0;
        dmu = (mu_max - mu_min) / Param::Nmu;
        inv_dmu = 1.0 / dmu;

        v_cells.resize(Param::Nv);
        v_faces.resize(Param::Nv + 1);
        v_widths.resize(Param::Nv);
        inv_v_widths.resize(Param::Nv);
        inv_v_center_dist.resize(Param::Nv + 1);
        mu_cells.resize(Param::Nmu);
        mu_faces.resize(Param::Nmu + 1);
        v2_cells.resize(Param::Nv);
        inv_v_cells.resize(Param::Nv);
        inv_v2_cells.resize(Param::Nv);
        v2_faces.resize(Param::Nv + 1);
        gamma_cells.resize(Param::Nv);
        beta_cells.resize(Param::Nv);
        speed_cells.resize(Param::Nv);
        chain_speed_cells.resize(Param::Nv);
        chain_speed_faces.resize(Param::Nv + 1);
        mu_face_factor.resize(Param::Nmu + 1);
        moment_weight.resize(Param::Nv);
        inv_moment_weight.resize(Param::Nv);
        mu_flux_scale.resize(Param::Nv);
        mu_cfl_factor.resize(Param::Nv);
        vx_cells.resize(Param::Nvmu);
        current_weight.resize(Param::Nvmu);

        max_abs_mu = 0.0;
        max_mu_face_factor = 0.0;
        max_inv_v = 0.0;
        max_speed = 0.0;

        const double refined_u =
            std::min(Param::momentum_refined_u, v_max);
        const int refined_cells =
            std::max(1, std::min(Param::momentum_refined_cells, Param::Nv - 1));
        for (int iv = 0; iv <= Param::Nv; ++iv) {
            if (iv <= refined_cells) {
                v_faces[iv] =
                    refined_u * static_cast<double>(iv) / refined_cells;
            } else {
                const int coarse_cells = Param::Nv - refined_cells;
                const int coarse_iv = iv - refined_cells;
                v_faces[iv] =
                    refined_u + (v_max - refined_u)
                              * static_cast<double>(coarse_iv) / coarse_cells;
            }
        }
        v_faces[0] = v_min;
        v_faces[Param::Nv] = v_max;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u_left = v_faces[iv];
            const double u_right = v_faces[iv + 1];
            const double width = u_right - u_left;
            const double u = 0.5 * (u_left + u_right);
            const double u_eff = std::max(u, Param::u_floor);
            const double gamma = std::sqrt(1.0 + u * u);
            const double beta = u / gamma;
            const double gamma_left = std::sqrt(1.0 + u_left * u_left);
            const double gamma_right = std::sqrt(1.0 + u_right * u_right);
            v_cells[iv] = u;
            v_widths[iv] = width;
            inv_v_widths[iv] = 1.0 / width;
            v2_cells[iv] = u_eff * u_eff;
            inv_v_cells[iv] = 1.0 / u_eff;
            inv_v2_cells[iv] = 1.0 / (u_eff * u_eff);
            gamma_cells[iv] = gamma;
            beta_cells[iv] = beta;
            speed_cells[iv] = Const::c * beta;
            chain_speed_cells[iv] =
                (width > 0.0)
                ? Const::c * (gamma_right - gamma_left) / width
                : speed_cells[iv];
            moment_weight[iv] =
                2.0 * Const::pi * dmu
                * (u_right * u_right * u_right
                   - u_left * u_left * u_left) / 3.0;
            inv_moment_weight[iv] =
                (moment_weight[iv] > 0.0) ? 1.0 / moment_weight[iv] : 0.0;
            mu_flux_scale[iv] = moment_weight[iv] * inv_dmu;
            max_inv_v = std::max(max_inv_v, inv_v_cells[iv]);
            max_speed = std::max(max_speed, speed_cells[iv]);
        }
        inv_v_center_dist[0] = inv_v_widths[0];
        inv_v_center_dist[Param::Nv] = inv_v_widths[Param::Nv - 1];
        for (int iv = 1; iv < Param::Nv; ++iv) {
            inv_v_center_dist[iv] =
                1.0 / (v_cells[iv] - v_cells[iv - 1]);
        }
        for (int iv = 0; iv <= Param::Nv; ++iv) {
            v2_faces[iv] = v_faces[iv] * v_faces[iv];
        }
        chain_speed_faces[0] = chain_speed_cells[0];
        chain_speed_faces[Param::Nv] = chain_speed_cells[Param::Nv - 1];
        for (int iv = 1; iv < Param::Nv; ++iv) {
            const double du = v_cells[iv] - v_cells[iv - 1];
            chain_speed_faces[iv] =
                (du > 0.0)
                ? Const::c * (gamma_cells[iv] - gamma_cells[iv - 1]) / du
                : 0.5 * (speed_cells[iv] + speed_cells[iv - 1]);
        }
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double mu = mu_min + (imu + 0.5) * dmu;
            mu_cells[imu] = mu;
            max_abs_mu = std::max(max_abs_mu, std::fabs(mu));
        }
        for (int imu = 0; imu <= Param::Nmu; ++imu) {
            const double muf = mu_min + imu * dmu;
            mu_faces[imu] = muf;
            mu_face_factor[imu] = 1.0 - muf * muf;
            max_mu_face_factor = std::max(max_mu_face_factor, mu_face_factor[imu]);
        }
        for (int iv = 0; iv < Param::Nv; ++iv) {
            mu_cfl_factor[iv] =
                max_mu_face_factor * inv_v_cells[iv] * inv_dmu;
        }
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k = static_cast<size_t>(iv) * Param::Nmu + imu;
                vx_cells[k] = speed_cells[iv] * mu_cells[imu];
                current_weight[k] =
                    moment_weight[iv] * chain_speed_cells[iv] * mu_cells[imu];
            }
        }
    }

    double v(int i) const { return v_cells[i]; }
    double v_face(int i) const { return v_faces[i]; }
    double mu(int j) const { return mu_cells[j]; }
    double mu_face(int j) const { return mu_faces[j]; }
};

// Cylindrical momentum grid used by the collisionless background-electron
// Vlasov--Ampere path.  The two stored indices retain Param::Nv/Param::Nmu
// extents, but represent (u_parallel, u_perp), not the legacy (u, mu) grid.
// Values in Species::f are cell-integrated masses M on this grid.
struct CylindricalVelocityGrid {
    std::vector<double> upar_faces;
    std::vector<double> upar_cells;
    std::vector<double> upar_widths;
    std::vector<double> uperp_faces;
    std::vector<double> uperp_cells;
    std::vector<double> uperp_ring_areas;
    std::vector<double> kinetic_energy;
    std::vector<double> vx;

    void init(double umax)
    {
        upar_faces.resize(Param::Nv + 1);
        upar_cells.resize(Param::Nv);
        upar_widths.resize(Param::Nv);
        uperp_faces.resize(Param::Nmu + 1);
        uperp_cells.resize(Param::Nmu);
        uperp_ring_areas.resize(Param::Nmu);
        kinetic_energy.resize(Param::Nvmu);
        vx.resize(Param::Nvmu);

        // A symmetric parallel-momentum grid removes the u=0 directional
        // degeneracy of the former spherical representation.
        for (int j = 0; j <= Param::Nv; ++j) {
            upar_faces[j] = -umax + 2.0 * umax * j / Param::Nv;
        }
        for (int k = 0; k <= Param::Nmu; ++k) {
            uperp_faces[k] = umax * k / Param::Nmu;
        }
        for (int j = 0; j < Param::Nv; ++j) {
            upar_cells[j] = 0.5 * (upar_faces[j] + upar_faces[j + 1]);
            upar_widths[j] = upar_faces[j + 1] - upar_faces[j];
        }
        for (int k = 0; k < Param::Nmu; ++k) {
            const double lo = uperp_faces[k];
            const double hi = uperp_faces[k + 1];
            uperp_cells[k] = 0.5 * (lo + hi);
            uperp_ring_areas[k] = Const::pi * (hi * hi - lo * lo);
        }
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t slot = static_cast<size_t>(j) * Param::Nmu + k;
                const double gamma = std::sqrt(1.0 +
                    upar_cells[j] * upar_cells[j] +
                    uperp_cells[k] * uperp_cells[k]);
                kinetic_energy[slot] = Const::me * Const::c * Const::c *
                                       (gamma - 1.0);
                vx[slot] = Const::c * upar_cells[j] / gamma;
            }
        }
    }

    double cell_phase_volume(int j, int k) const
    {
        return upar_widths[j] * uperp_ring_areas[k];
    }
};

struct SpatialGrid {
    int nx_global;
    int nx_local;
    int ix_start;
    int nghost;
    int nx_total;
    double dx;
    double x_min;

    void init(int rank, int nranks) {
        nx_global = Param::nx;
        dx = Param::dx;
        x_min = 0.0;
        nghost = Param::Nghost;

        nx_local = nx_global / nranks;
        int remainder = nx_global % nranks;
        if (rank < remainder) {
            nx_local += 1;
            ix_start = rank * nx_local;
        } else {
            ix_start = remainder * (nx_local + 1) + (rank - remainder) * nx_local;
        }
        nx_total = nx_local + 2 * nghost;
    }

    double x(int i_local) const {
        int ig = ix_start + (i_local - nghost);
        return x_min + (ig + 0.5) * dx;
    }
};

inline double gamma_from_v(double v) {
    double beta = v / Const::c;
    if (beta > 0.999999999999) beta = 0.999999999999;
    return 1.0 / std::sqrt(1.0 - beta * beta);
}

inline double gamma_from_u(double u) {
    return std::sqrt(1.0 + u * u);
}

inline double speed_from_u(double u) {
    return Const::c * u / gamma_from_u(u);
}

inline double u_from_v(double v) {
    const double beta = std::max(-0.999999999999,
                                 std::min(0.999999999999, v / Const::c));
    return beta / std::sqrt(1.0 - beta * beta);
}

inline double momentum_from_v(double v, double mass) {
    return gamma_from_v(v) * mass * v;
}

inline size_t idx3(int ix, int iv, int imu) {
    return static_cast<size_t>(ix) * Param::Nvmu
         + static_cast<size_t>(iv) * Param::Nmu
         + static_cast<size_t>(imu);
}

inline size_t idx2(int iv, int imu) {
    return static_cast<size_t>(iv) * Param::Nmu
         + static_cast<size_t>(imu);
}

#endif
