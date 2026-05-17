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

    std::vector<double> v_cells;
    std::vector<double> v_faces;
    std::vector<double> mu_cells;
    std::vector<double> mu_faces;
    std::vector<double> v2_cells;
    std::vector<double> inv_v_cells;
    std::vector<double> inv_v2_cells;
    std::vector<double> v2_faces;
    std::vector<double> mu_face_factor;
    std::vector<double> moment_weight;
    std::vector<double> vx_cells;
    std::vector<double> current_weight;

    void init(double vmax)
    {
        v_min = 0.0;
        v_max = vmax;
        dv = (v_max - v_min) / Param::Nv;
        inv_dv = 1.0 / dv;
        mu_min = -1.0;
        mu_max = 1.0;
        dmu = (mu_max - mu_min) / Param::Nmu;
        inv_dmu = 1.0 / dmu;

        v_cells.resize(Param::Nv);
        v_faces.resize(Param::Nv + 1);
        mu_cells.resize(Param::Nmu);
        mu_faces.resize(Param::Nmu + 1);
        v2_cells.resize(Param::Nv);
        inv_v_cells.resize(Param::Nv);
        inv_v2_cells.resize(Param::Nv);
        v2_faces.resize(Param::Nv + 1);
        mu_face_factor.resize(Param::Nmu + 1);
        moment_weight.resize(Param::Nv);
        vx_cells.resize(Param::Nvmu);
        current_weight.resize(Param::Nvmu);

        max_abs_mu = 0.0;
        max_mu_face_factor = 0.0;
        max_inv_v = 0.0;

        const double weight0 = 2.0 * Const::pi * dv * dmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double v = v_min + (iv + 0.5) * dv;
            const double v_eff = std::max(v, Param::v_floor);
            v_cells[iv] = v;
            v2_cells[iv] = v_eff * v_eff;
            inv_v_cells[iv] = 1.0 / v_eff;
            inv_v2_cells[iv] = 1.0 / (v_eff * v_eff);
            moment_weight[iv] = weight0 * v * v;
            max_inv_v = std::max(max_inv_v, inv_v_cells[iv]);
        }
        for (int iv = 0; iv <= Param::Nv; ++iv) {
            const double vf = v_min + iv * dv;
            v_faces[iv] = vf;
            v2_faces[iv] = vf * vf;
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
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k = static_cast<size_t>(iv) * Param::Nmu + imu;
                vx_cells[k] = v_cells[iv] * mu_cells[imu];
                current_weight[k] = moment_weight[iv] * vx_cells[k];
            }
        }
    }

    double v(int i) const { return v_cells[i]; }
    double v_face(int i) const { return v_faces[i]; }
    double mu(int j) const { return mu_cells[j]; }
    double mu_face(int j) const { return mu_faces[j]; }
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
