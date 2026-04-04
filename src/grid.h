#ifndef GRID_H
#define GRID_H

#include "parameters.h"
#include <vector>
#include <cmath>
#include <cstddef>

// ============================================================================
// 单个物种的相空间网格
// ============================================================================
struct MomentumGrid {
    double px_min, px_max, dpx;
    double py_min, py_max, dpy;
    double pz_min, pz_max, dpz;

    void init(double pxmin, double pxmax, double pymin, double pymax,
              double pzmin, double pzmax)
    {
        px_min = pxmin; px_max = pxmax;
        py_min = pymin; py_max = pymax;
        pz_min = pzmin; pz_max = pzmax;
        dpx = (px_max - px_min) / Param::Npx;
        dpy = (py_max - py_min) / Param::Npy;
        dpz = (pz_max - pz_min) / Param::Npz;
    }

    // 单元中心动量值
    double px(int i) const { return px_min + (i + 0.5) * dpx; }
    double py(int j) const { return py_min + (j + 0.5) * dpy; }
    double pz(int k) const { return pz_min + (k + 0.5) * dpz; }

    // 单元面动量值
    double px_face(int i) const { return px_min + i * dpx; }
    double py_face(int j) const { return py_min + j * dpy; }
    double pz_face(int k) const { return pz_min + k * dpz; }
};

// ============================================================================
// 空间网格（1D，含 MPI 分解）
// ============================================================================
struct SpatialGrid {
    int nx_global;         // 全局单元总数
    int nx_local;          // 本 rank 的本地单元数
    int ix_start;          // 本地首单元的全局索引
    int nghost;            // 每侧幽灵单元数
    int nx_total;          // nx_local + 2*nghost
    double dx;
    double x_min;

    void init(int rank, int nranks) {
        nx_global = Param::nx;
        dx = Param::dx;
        x_min = 0.0;
        nghost = Param::Nghost;

        // 沿 x 方向进行简单均分
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

    // 单元中心 x 位置（局部索引，包含幽灵单元）
    double x(int i_local) const {
        int ig = ix_start + (i_local - nghost);
        return x_min + (ig + 0.5) * dx;
    }
};

// ============================================================================
// 相对论运动学辅助函数
// ============================================================================
inline double gamma_from_p(double px, double py, double pz, double mass) {
    double p2 = px * px + py * py + pz * pz;
    double mc = mass * Const::c;
    return std::sqrt(1.0 + p2 / (mc * mc));
}

inline double vx_from_p(double px, double py, double pz, double mass) {
    double gam = gamma_from_p(px, py, pz, mass);
    return px / (gam * mass);
}

inline double vy_from_p(double px, double py, double pz, double mass) {
    double gam = gamma_from_p(px, py, pz, mass);
    return py / (gam * mass);
}

inline double vz_from_p(double px, double py, double pz, double mass) {
    double gam = gamma_from_p(px, py, pz, mass);
    return pz / (gam * mass);
}

// ============================================================================
// 4D 相空间索引：f[ix][ipx][ipy][ipz] 按一维数组存储
// 布局：ix 变化最慢，ipz 变化最快
// ============================================================================
inline size_t idx4(int ix, int ipx, int ipy, int ipz) {
    return static_cast<size_t>(ix) * Param::Np3
         + static_cast<size_t>(ipx) * Param::Npy * Param::Npz
         + static_cast<size_t>(ipy) * Param::Npz
         + static_cast<size_t>(ipz);
}

// 仅动量 3D 索引
inline size_t idx3(int ipx, int ipy, int ipz) {
    return static_cast<size_t>(ipx) * Param::Npy * Param::Npz
         + static_cast<size_t>(ipy) * Param::Npz
         + static_cast<size_t>(ipz);
}

#endif // 头文件保护：GRID_H
