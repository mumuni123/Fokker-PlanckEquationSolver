#include "collision.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>

// ============================================================================
// 误差函数近似（C++11 提供 std::erf，这里保留封装）
// ============================================================================
double CollisionOperator::erf_approx(double x)
{
    return std::erf(x);
}

// ============================================================================
// Chandrasekhar 函数：G(x) = [erf(x) - x*erf'(x)] / (2x^2)
// 其中 erf'(x) = (2/sqrt(pi)) * exp(-x^2)
// ============================================================================
double CollisionOperator::chandrasekhar_G(double x)
{
    if (std::abs(x) < 1.0e-6) {
        // 小 x 的泰勒展开：G(x) ≈ 2/(3*sqrt(pi)) * (1 - 3x²/5)
        return 2.0 / (3.0 * std::sqrt(Const::pi));
    }
    double erfx = std::erf(x);
    double erfpx = (2.0 / std::sqrt(Const::pi)) * std::exp(-x * x);
    return (erfx - x * erfpx) / (2.0 * x * x);
}

// ============================================================================
// 库仑对数（自动计算，参考 NRL 手册）
// 对电子-电子：ln Λ ≈ 23.5 - ln(n_e^{1/2} * T_e^{-5/4})
// 简化形式：ln Λ ≈ log(12*pi*n*λ_D^3)
// ============================================================================
double CollisionOperator::coulomb_log(double n_e, double T_e) const
{
    // 德拜长度：λ_D = sqrt(ε₀ * kT / (n_e * e²))
    double lambda_D = std::sqrt(Const::eps0 * T_e / (n_e * Const::qe * Const::qe));
    double arg = 12.0 * Const::pi * n_e * lambda_D * lambda_D * lambda_D;
    double lnL = std::log(arg);
    if (lnL < 2.0) lnL = 2.0;  // 物理下限
    return lnL;
}

// ============================================================================
// 测试粒子（动量模长 |p|）与麦克斯韦背景粒子散射时的碰撞率
//
// NRL Plasma Formulary 记号：
//   x = v²/(2*v_th²)，其中 v_th = sqrt(kT/m_field)
//   ν₀ = n_f * Z_t² * Z_f² * e⁴ * lnΛ / (4π * ε₀² * m_t² * v³)
//
// 频率：
//   ν_s = (1 + m_t/m_f) * Ψ(x) * ν₀
//   ν_∥ = Ψ(x) / x * ν₀
//   ν_⊥ = [(1 - 1/(2x))Ψ(x) + Ψ'(x)] / x * ν₀
//
// 其中 Ψ(x) = G(√x) = Chandrasekhar 函数
// ============================================================================
CollisionOperator::CollisionRates
CollisionOperator::compute_rates(double p_mag, double mass_test,
                                 double n_field, double T_field,
                                 double m_field, double lnLambda,
                                 double Z_test, double Z_field) const
{
    CollisionRates rates;

    // 测试粒子速度（对背景采用非相对论近似）
    double v = p_mag / mass_test;
    if (v < 1.0e-20) {
        rates.nu_s = 0.0;
        rates.nu_parallel = 0.0;
        rates.nu_perp = 0.0;
        return rates;
    }

    // 背景粒子热速度
    double v_th = std::sqrt(T_field / m_field);
    if (v_th < 1.0e-20) {
        rates.nu_s = 0.0;
        rates.nu_parallel = 0.0;
        rates.nu_perp = 0.0;
        return rates;
    }

    double x = v / (std::sqrt(2.0) * v_th);  // 归一化速度

    // 基础碰撞频率
    double qe4 = Const::qe * Const::qe * Const::qe * Const::qe;
    double Zt2Zf2 = Z_test * Z_test * Z_field * Z_field;
    double eps02 = Const::eps0 * Const::eps0;
    double nu0 = n_field * Zt2Zf2 * qe4 * lnLambda
               / (4.0 * Const::pi * eps02 * mass_test * mass_test * v * v * v);

    // Chandrasekhar 函数及其导数
    double Gx = chandrasekhar_G(x);
    double erfx = std::erf(x);
    double erfpx = (2.0 / std::sqrt(Const::pi)) * std::exp(-x * x);

    // NRL 中的 Ψ(x) = erf(x) - (2x/√π)e^{-x²} = erf(x) - x*erfp(x)
    // 等价地，按 Chandrasekhar G 可写为 Ψ(x) = 2*x² * G(x)
    double Psi = 2.0 * x * x * Gx;

    // 减速率
    rates.nu_s = (1.0 + mass_test / m_field) * Psi * nu0;

    // 平行扩散率
    if (x > 1.0e-6) {
        rates.nu_parallel = Psi / (x * x) * nu0;
    } else {
        rates.nu_parallel = (4.0 / (3.0 * std::sqrt(Const::pi))) * nu0;
    }

    // 垂直扩散率
    // ν_⊥ = [(erfx - Ψ)/(2x²) + erfpx/(2x)] * ... 简化为：
    // ν_⊥ = [erf(x) - Psi] / x² * ν₀  （NRL 近似）
    if (x > 1.0e-6) {
        rates.nu_perp = (erfx - Psi) / (x * x) * nu0;
    } else {
        rates.nu_perp = (4.0 / (3.0 * std::sqrt(Const::pi))) * nu0;
    }

    return rates;
}

// ============================================================================
// Chang-Cooper 一维扫掠
//
// 求解：∂f/∂t = -∂/∂p [A(p)*f - ½ ∂/∂p(B(p)*f)]
//       = -∂J/∂p
// 其中通量 J_{j+1/2} = A_{j+1/2} * f_{j+1/2} - B_{j+1/2}/(2Δp) * (f_{j+1}-f_j)
//
// Chang-Cooper 插值：
//   f_{j+1/2} = (1-δ_j)*f_{j+1} + δ_j*f_j
//   其中 δ_j = 1/w - 1/(e^w - 1), w = A*Δp / (B/2)
//
// 为稳定性采用隐式（后向 Euler）求解：
//   f^{n+1}_j = f^n_j - (dt/dp) * (J^{n+1}_{j+1/2} - J^{n+1}_{j-1/2})
//
// 最终形成三对角线性系统。
// ============================================================================
void CollisionOperator::chang_cooper_1d(double* f_line, int N, double dp,
                                        const double* drag, const double* diff,
                                        double dt) const
{
    if (N < 3) return;

    // 在单元面计算 Chang-Cooper 权重
    // drag[j] 与 diff[j] 定义在单元中心，需插值到单元面
    std::vector<double> A_face(N + 1, 0.0);  // 面 j+1/2 上的 A
    std::vector<double> B_face(N + 1, 0.0);  // 面上的 B
    std::vector<double> delta(N + 1, 0.5);   // CC 权重

    for (int j = 0; j < N - 1; ++j) {
        A_face[j + 1] = 0.5 * (drag[j] + drag[j + 1]);
        B_face[j + 1] = 0.5 * (diff[j] + diff[j + 1]);
    }
    // 边界面：外推
    A_face[0] = drag[0];
    B_face[0] = diff[0];
    A_face[N] = drag[N - 1];
    B_face[N] = diff[N - 1];

    // 计算 Chang-Cooper 参数 δ
    for (int j = 0; j <= N; ++j) {
        double B_half = B_face[j] * 0.5;
        if (std::abs(B_half) < 1.0e-40) {
            delta[j] = 0.5;
            continue;
        }
        double w = A_face[j] * dp / B_half;
        if (std::abs(w) < 1.0e-8) {
            delta[j] = 0.5 - w / 12.0;  // 泰勒展开
        } else if (w > 500.0) {
            delta[j] = 0.0;
        } else if (w < -500.0) {
            delta[j] = 1.0;
        } else {
            delta[j] = 1.0 / w - 1.0 / (std::exp(w) - 1.0);
        }
        // 限幅
        if (delta[j] < 0.0) delta[j] = 0.0;
        if (delta[j] > 1.0) delta[j] = 1.0;
    }

    // 构建三对角系统：a_j f_{j-1} + b_j f_j + c_j f_{j+1} = f^n_j
    std::vector<double> a(N, 0.0), b(N, 0.0), cc(N, 0.0), rhs(N, 0.0);

    double dtdp = dt / dp;

    for (int j = 0; j < N; ++j) {
        // j+1/2 处通量系数：
        // J_{j+1/2} = A_{j+1/2} * [(1-δ_{j+1/2})*f_{j+1} + δ_{j+1/2}*f_j]
        //           - B_{j+1/2}/(2*dp) * (f_{j+1} - f_j)
        // J_{j+1/2} 中 f_j 的系数：
        double alpha_R = A_face[j + 1] * delta[j + 1] + B_face[j + 1] / (2.0 * dp);
        // J_{j+1/2} 中 f_{j+1} 的系数：
        double beta_R  = A_face[j + 1] * (1.0 - delta[j + 1]) - B_face[j + 1] / (2.0 * dp);

        // J_{j-1/2} 中 f_j 的系数：
        double beta_L  = A_face[j] * (1.0 - delta[j]) - B_face[j] / (2.0 * dp);
        // J_{j-1/2} 中 f_{j-1} 的系数：
        double alpha_L = A_face[j] * delta[j] + B_face[j] / (2.0 * dp);

        // f_j^{n+1} = f_j^n - dtdp * (J_{j+1/2} - J_{j-1/2})
        // -dtdp * J_{j+1/2} 贡献：-dtdp * [alpha_R * f_j + beta_R * f_{j+1}]
        // +dtdp * J_{j-1/2} 贡献：+dtdp * [alpha_L * f_{j-1} + beta_L * f_j]
        // 整理为隐式形式：
        // (1 + dtdp*(alpha_R - beta_L)) * f_j - dtdp*beta_R * f_{j+1} - dtdp*alpha_L * f_{j-1} = f_j^n

        a[j]   = -dtdp * alpha_L;        // 下对角（f_{j-1}）
        b[j]   = 1.0 + dtdp * (alpha_R - beta_L);  // 主对角（f_j）
        cc[j]  = -dtdp * beta_R;         // 上对角（f_{j+1}）
        rhs[j] = f_line[j];
    }

    // 边界条件：边界处 f = 0（零通量）
    // j=0: a[0] = 0（不存在 f_{-1}），A_face[0] 的边界项已处理
    a[0] = 0.0;
    cc[N - 1] = 0.0;

    // Thomas 算法求解三对角系统
    std::vector<double> cp(N, 0.0), dp_vec(N, 0.0);
    cp[0] = cc[0] / b[0];
    dp_vec[0] = rhs[0] / b[0];

    for (int j = 1; j < N; ++j) {
        double m = b[j] - a[j] * cp[j - 1];
        if (std::abs(m) < 1.0e-40) m = 1.0e-40;
        cp[j] = cc[j] / m;
        dp_vec[j] = (rhs[j] - a[j] * dp_vec[j - 1]) / m;
    }

    f_line[N - 1] = dp_vec[N - 1];
    for (int j = N - 2; j >= 0; --j) {
        f_line[j] = dp_vec[j] - cp[j] * f_line[j + 1];
    }

    // 正定性约束
    for (int j = 0; j < N; ++j) {
        if (f_line[j] < 0.0) f_line[j] = 0.0;
    }
}

// ============================================================================
// 对物种施加碰撞算子
// 维度分裂：沿 px、py、pz 分别执行 1D Chang-Cooper
// ============================================================================
void CollisionOperator::apply(Species& sp, double dt,
                              double n_field, double T_field, double m_field,
                              double Z_field, double Z_test)
{
    double lnL = coulomb_log(n_field, T_field);

    int ng  = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;

    // ---- 沿 px 扫掠 ----
    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int ipy = 0; ipy < Param::Npy; ++ipy) {
            int ix_g = ix + ng;
            double ppy = sp.pgrid.py(ipy);

            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = sp.pgrid.pz(ipz);

                // 提取沿 px 的 1D 线
                std::vector<double> f_line(Param::Npx);
                std::vector<double> drag_px(Param::Npx);
                std::vector<double> diff_px(Param::Npx);

                for (int ipx = 0; ipx < Param::Npx; ++ipx) {
                    f_line[ipx] = sp.f[idx4(ix_g, ipx, ipy, ipz)];

                    double ppx = sp.pgrid.px(ipx);
                    double p_mag = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz);

                    CollisionRates rates = compute_rates(
                        p_mag, sp.mass, n_field, T_field, m_field,
                        lnL, Z_test, Z_field);

                    // px 方向拖曳：D1_x = -ν_s * px/|p|
                    double p_hat_x = (p_mag > 1.0e-30) ? ppx / p_mag : 0.0;
                    drag_px[ipx] = -rates.nu_s * p_hat_x * p_mag;
                    // 实际拖曳系数 A = -ν_s * p（有量纲）
                    // 在 1D 投影下：A_px = -ν_s * px

                    // px 方向扩散：D2_xx = ν_perp*(1 - px²/p²) + ν_parallel*(px²/p²)
                    double cos2 = (p_mag > 1.0e-30) ? (ppx*ppx)/(p_mag*p_mag) : 0.0;
                    diff_px[ipx] = (rates.nu_perp * (1.0 - cos2) + rates.nu_parallel * cos2)
                                 * p_mag * p_mag;
                    // 乘以 p² 是因为 p 空间扩散系数量纲为 p²/time
                }

                chang_cooper_1d(f_line.data(), Param::Npx, sp.pgrid.dpx,
                                drag_px.data(), diff_px.data(), dt);

                // 回写结果
                for (int ipx = 0; ipx < Param::Npx; ++ipx) {
                    sp.f[idx4(ix_g, ipx, ipy, ipz)] = f_line[ipx];
                }
            }
        }
    }

    // ---- 沿 py 扫掠 ----
    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            int ix_g = ix + ng;
            double ppx = sp.pgrid.px(ipx);

            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = sp.pgrid.pz(ipz);

                std::vector<double> f_line(Param::Npy);
                std::vector<double> drag_py(Param::Npy);
                std::vector<double> diff_py(Param::Npy);

                for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                    f_line[ipy] = sp.f[idx4(ix_g, ipx, ipy, ipz)];

                    double ppy = sp.pgrid.py(ipy);
                    double p_mag = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz);

                    CollisionRates rates = compute_rates(
                        p_mag, sp.mass, n_field, T_field, m_field,
                        lnL, Z_test, Z_field);

                    double p_hat_y = (p_mag > 1.0e-30) ? ppy / p_mag : 0.0;
                    drag_py[ipy] = -rates.nu_s * p_hat_y * p_mag;

                    double cos2 = (p_mag > 1.0e-30) ? (ppy*ppy)/(p_mag*p_mag) : 0.0;
                    diff_py[ipy] = (rates.nu_perp * (1.0 - cos2) + rates.nu_parallel * cos2)
                                 * p_mag * p_mag;
                }

                chang_cooper_1d(f_line.data(), Param::Npy, sp.pgrid.dpy,
                                drag_py.data(), diff_py.data(), dt);

                for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                    sp.f[idx4(ix_g, ipx, ipy, ipz)] = f_line[ipy];
                }
            }
        }
    }

    // ---- 沿 pz 扫掠 ----
    #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            int ix_g = ix + ng;
            double ppx = sp.pgrid.px(ipx);

            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                double ppy = sp.pgrid.py(ipy);

                std::vector<double> f_line(Param::Npz);
                std::vector<double> drag_pz(Param::Npz);
                std::vector<double> diff_pz(Param::Npz);

                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    f_line[ipz] = sp.f[idx4(ix_g, ipx, ipy, ipz)];

                    double ppz = sp.pgrid.pz(ipz);
                    double p_mag = std::sqrt(ppx*ppx + ppy*ppy + ppz*ppz);

                    CollisionRates rates = compute_rates(
                        p_mag, sp.mass, n_field, T_field, m_field,
                        lnL, Z_test, Z_field);

                    double p_hat_z = (p_mag > 1.0e-30) ? ppz / p_mag : 0.0;
                    drag_pz[ipz] = -rates.nu_s * p_hat_z * p_mag;

                    double cos2 = (p_mag > 1.0e-30) ? (ppz*ppz)/(p_mag*p_mag) : 0.0;
                    diff_pz[ipz] = (rates.nu_perp * (1.0 - cos2) + rates.nu_parallel * cos2)
                                 * p_mag * p_mag;
                }

                chang_cooper_1d(f_line.data(), Param::Npz, sp.pgrid.dpz,
                                drag_pz.data(), diff_pz.data(), dt);

                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    sp.f[idx4(ix_g, ipx, ipy, ipz)] = f_line[ipz];
                }
            }
        }
    }
}
