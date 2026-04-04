#ifndef MAXWELL_H
#define MAXWELL_H

#include "parameters.h"
#include "grid.h"
#include <vector>

// ============================================================================
// 1D Yee 网格 FDTD 电磁场求解器
//
// 网格布局（1D）：
//   Ex[i], Ey[i], Ez[i] 定义在整数节点 i = 0, 1, ..., nx_total-1
//   By[i], Bz[i] 定义在半整数节点 i+1/2，并存放在索引 i
//     即 By[i] = By_{i+1/2}，其中 i = 0, ..., nx_total-2
//
// 在 1D 情况下（∂/∂y = ∂/∂z = 0）：
//   ∂By/∂t = +∂Ez/∂x       （符号来自 ∇×E）
//   ∂Bz/∂t = -∂Ey/∂x
//   ∂Ex/∂t = -Jx/ε₀
//   ∂Ey/∂t = -c²∂Bz/∂x - Jy/ε₀
//   ∂Ez/∂t = +c²∂By/∂x - Jz/ε₀
// ============================================================================

struct EMFields {
    // 局部网格上的场（包含幽灵单元）
    std::vector<double> Ex, Ey, Ez;  // 整数节点，长度 = nx_total
    std::vector<double> By, Bz;      // 半整数节点，长度 = nx_total

    // 总电荷/电流密度（所有物种求和）
    std::vector<double> rho;         // 整数节点处的电荷密度
    std::vector<double> Jx, Jy, Jz; // 整数节点处的电流密度

    int nx_total;
    double dx;

    void init(const SpatialGrid& sg);
    void zero_currents();
    void accumulate_moments(const class Species& sp);

    // Leapfrog 更新
    // B 位于半时间层：B^{n+1/2}
    // E 位于整数时间层：E^{n}
    void update_B(double dt);  // B^{n-1/2} -> B^{n+1/2}
    void update_E(double dt);  // E^{n} -> E^{n+1}，使用 B^{n+1/2} 与 J^{n+1/2}

    // 使用高斯定律修正 Ex（泊松修正）
    void gauss_correction();

    // 场量幽灵单元交换
    void exchange_ghosts(int mpi_rank, int mpi_size);

    // 总场能
    double total_energy() const;
};

#endif // 头文件保护：MAXWELL_H
