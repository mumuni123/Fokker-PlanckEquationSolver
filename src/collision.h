#ifndef COLLISION_H
#define COLLISION_H

#include "species.h"
#include "grid.h"
#include <vector>

// ============================================================================
// 线性化 Fokker-Planck 碰撞算子
//
// 仅作用于背景物种（根据 input.deck，束流碰撞关闭）。
// 使用 Chang-Cooper 离散以保证守恒性质。
//
// 动量空间中的 FP 方程：
//   ∂f/∂t = -∇_p · (f * D1) + ½ ∇_p ∇_p : (f * D2)
//
// 对于麦克斯韦场中的测试粒子：
//   D1_i = -ν_s * p_i / |p|    （减速拖曳）
//   D2_{ij} = ν_perp * (δ_{ij} - p_i*p_j/p²) + ν_parallel * p_i*p_j/p²
//
// 采用维度分裂：沿 px、py、pz 分别做 1D Chang-Cooper 扫掠，
// 使用拖曳与扩散的对角分量。
// ============================================================================

class CollisionOperator {
public:
    // 对物种 sp（测试粒子）施加 FP 碰撞，背景参数由给定密度和温度确定。
    // 仅修改内部单元（不含幽灵单元）。
    void apply(Species& sp, double dt,
               double n_field, double T_field, double m_field,
               double Z_field, double Z_test);

private:
    // 计算库仑对数
    double coulomb_log(double n_e, double T_e) const;

    // 测试粒子在动量 |p| 下的碰撞频率
    // 与温度 T、密度 n、质量 m 的麦克斯韦背景发生散射
    struct CollisionRates {
        double nu_s;        // 减速率
        double nu_parallel; // 平行扩散率
        double nu_perp;     // 垂直扩散率
    };

    CollisionRates compute_rates(double p_mag, double mass_test,
                                 double n_field, double T_field,
                                 double m_field, double lnLambda,
                                 double Z_test, double Z_field) const;

    // 沿某一动量方向执行 Chang-Cooper 1D 扫掠
    void chang_cooper_1d(double* f_line, int N, double dp,
                         const double* drag, const double* diff,
                         double dt) const;

    // 误差函数及相关函数
    static double erf_approx(double x);
    static double chandrasekhar_G(double x);
};

#endif // 头文件保护：COLLISION_H
