#ifndef VLASOV_H
#define VLASOV_H

#include "species.h"
#include "grid.h"
#include <vector>

// 前向声明
struct EMFields;

class VlasovSolver {
public:
    void advect(Species& sp, const SpatialGrid& sg, const EMFields& fields,
                double dt, int mpi_rank, int mpi_size);

private:
    // 维度分裂扫掠
    void advect_x(Species& sp, const SpatialGrid& sg, double dt,
                  int mpi_rank, int mpi_size);
    void advect_px(Species& sp, const SpatialGrid& sg,
                   const std::vector<double>& Ex,
                   const std::vector<double>& Ey,
                   const std::vector<double>& Ez,
                   const std::vector<double>& By,
                   const std::vector<double>& Bz,
                   double dt);
    void advect_py(Species& sp, const SpatialGrid& sg,
                   const std::vector<double>& Ex,
                   const std::vector<double>& Ey,
                   const std::vector<double>& Ez,
                   const std::vector<double>& By,
                   const std::vector<double>& Bz,
                   double dt);
    void advect_pz(Species& sp, const SpatialGrid& sg,
                   const std::vector<double>& Ex,
                   const std::vector<double>& Ey,
                   const std::vector<double>& Ez,
                   const std::vector<double>& By,
                   const std::vector<double>& Bz,
                   double dt);

    // WENO-5 重构：返回左右界面值
    // 输入：单元平均值 q[0..N-1]，两侧各有 3 个幽灵单元
    // 输出：内部单元的通量差分
    void weno5_sweep(const double* q, double* rhs, int N, int Ntot,
                     const double* vel, double dxi, int stride);

    // x 方向通过 MPI 交换幽灵单元
    void exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                           int mpi_rank, int mpi_size);
};

// 面 i+1/2 处左偏 WENO-5 重构
// 使用模板 q[i-2], q[i-1], q[i], q[i+1], q[i+2]
inline double weno5_left(double qm2, double qm1, double q0, double qp1, double qp2)
{
    // 三个候选多项式
    double p0 = ( 2.0*qm2 - 7.0*qm1 + 11.0*q0) / 6.0;
    double p1 = (-1.0*qm1 + 5.0*q0  +  2.0*qp1) / 6.0;
    double p2 = ( 2.0*q0  + 5.0*qp1  -  1.0*qp2) / 6.0;

    // 光滑性指标
    double b0 = (13.0/12.0) * (qm2 - 2.0*qm1 + q0)  * (qm2 - 2.0*qm1 + q0)
              + (1.0/4.0)   * (qm2 - 4.0*qm1 + 3.0*q0) * (qm2 - 4.0*qm1 + 3.0*q0);
    double b1 = (13.0/12.0) * (qm1 - 2.0*q0  + qp1) * (qm1 - 2.0*q0  + qp1)
              + (1.0/4.0)   * (qm1 - qp1) * (qm1 - qp1);
    double b2 = (13.0/12.0) * (q0  - 2.0*qp1 + qp2) * (q0  - 2.0*qp1 + qp2)
              + (1.0/4.0)   * (3.0*q0 - 4.0*qp1 + qp2) * (3.0*q0 - 4.0*qp1 + qp2);

    // 理想权重
    const double d0 = 1.0/10.0, d1 = 6.0/10.0, d2 = 3.0/10.0;
    const double eps = 1.0e-6;

    // 非线性权重
    double a0 = d0 / ((eps + b0) * (eps + b0));
    double a1 = d1 / ((eps + b1) * (eps + b1));
    double a2 = d2 / ((eps + b2) * (eps + b2));
    double asum = a0 + a1 + a2;

    return (a0 * p0 + a1 * p1 + a2 * p2) / asum;
}

// 面 i+1/2 处右偏 WENO-5 重构
// 右偏值使用左偏模板的镜像。
// 输入：q[i-1], q[i], q[i+1], q[i+2], q[i+3]（模板右移）
// 等价于对反序序列应用 weno5_left。
inline double weno5_right(double qm1, double q0, double qp1, double qp2, double qp3)
{
    // 三个候选多项式（左偏的镜像）
    double p0 = ( 2.0*qp3 - 7.0*qp2 + 11.0*qp1) / 6.0;
    double p1 = (-1.0*qp2 + 5.0*qp1 +  2.0*q0)  / 6.0;
    double p2 = ( 2.0*qp1 + 5.0*q0  -  1.0*qm1) / 6.0;

    // 光滑性指标（同一公式，索引镜像）
    double b0 = (13.0/12.0) * (qp3 - 2.0*qp2 + qp1) * (qp3 - 2.0*qp2 + qp1)
              + (1.0/4.0)   * (qp3 - 4.0*qp2 + 3.0*qp1) * (qp3 - 4.0*qp2 + 3.0*qp1);
    double b1 = (13.0/12.0) * (qp2 - 2.0*qp1 + q0) * (qp2 - 2.0*qp1 + q0)
              + (1.0/4.0)   * (qp2 - q0) * (qp2 - q0);
    double b2 = (13.0/12.0) * (qp1 - 2.0*q0  + qm1) * (qp1 - 2.0*q0  + qm1)
              + (1.0/4.0)   * (3.0*qp1 - 4.0*q0 + qm1) * (3.0*qp1 - 4.0*q0 + qm1);

    // 理想权重（与左偏相同：d0=1/10, d1=6/10, d2=3/10）
    const double d0 = 1.0/10.0, d1 = 6.0/10.0, d2 = 3.0/10.0;
    const double eps = 1.0e-6;

    double a0 = d0 / ((eps + b0) * (eps + b0));
    double a1 = d1 / ((eps + b1) * (eps + b1));
    double a2 = d2 / ((eps + b2) * (eps + b2));
    double asum = a0 + a1 + a2;

    return (a0 * p0 + a1 * p1 + a2 * p2) / asum;
}

#endif // 头文件保护：VLASOV_H
