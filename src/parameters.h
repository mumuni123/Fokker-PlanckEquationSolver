#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <cstddef>

// ============================================================================
// 物理常数（SI 单位）
// ============================================================================
namespace Const {
    const double me   = 9.10938e-31;   // 电子质量 [kg]
    const double qe   = 1.60218e-19;   // 元电荷 [C]
    const double c    = 2.99792e8;     // 光速 [m/s]
    const double eps0 = 8.85419e-12;   // 真空介电常数 [F/m]
    const double mu0  = 1.25664e-6;    // 真空磁导率 [H/m]
    const double kB   = 1.38065e-23;   // 玻尔兹曼常数 [J/K]
    const double pi   = 3.14159265358979323846;
    const double eV   = 1.60218e-19;   // 1 eV 对应焦耳
    const double micro = 1.0e-6;       // 1 微米
    const double femto = 1.0e-15;      // 1 飞秒
}

// ============================================================================
// 仿真参数（来自 EPOCH input.deck）
// ============================================================================
namespace Param {
    // --- 背景等离子体 ---
    const double temperature_e = 100.0 * Const::eV;  // 100 eV -> 焦耳
    const double temperature_i = 10.0  * Const::eV;   // 10 eV -> 焦耳
    const double dens          = 1.2e29;               // 背景电子密度 [/m^3]
    const int    Z_ion         = 2;                    // 离子电荷数
    const double mass_ion_me   = 49572.0;              // 离子质量（以电子质量为单位）
    const double mass_ion      = mass_ion_me * Const::me;

    // --- 束流 ---
    const double jb        = 5.0e16;                   // 束流电流密度 [A/m^2]
    const double gambetab  = 8.60;                     // 束流的 gamma*beta
    const double gamb      = std::sqrt(1.0 + gambetab * gambetab);
    const double betab     = gambetab / gamb;
    const double densb     = jb / (Const::qe * betab * Const::c); // 束流密度
    const double beam_drift_px = gambetab * Const::me * Const::c; // 束流漂移动量

    // --- 空间网格 ---
    const double dx    = 0.002 * Const::micro;         // 2 nm
    const double Lx    = 8.0   * Const::micro;         // 8 um
    const int    nx    = static_cast<int>(Lx / dx);     // 4000

    // --- 时间参数 ---
    const double t_end         = 120.0 * Const::femto;
    const double t_inject_end  = 25.0  * Const::femto;  // 束流注入持续时间
    const double dt_multiplier = 0.5;                    // CFL 安全系数
    const double dt_snapshot   = 0.6 * Const::femto;     // 输出间隔

    // --- 动量网格（每个维度） ---
    const int Npx = 32;
    const int Npy = 32;
    const int Npz = 32;
    const size_t Np3 = static_cast<size_t>(Npx) * Npy * Npz;  // 32768

    // 各物种的动量网格范围：
    // 束流：px in [-1*gambetab*me*c, 3*gambetab*me*c]，py/pz in [-gambetab*me*c, gambetab*me*c]
    // 背景电子：px,py,pz in [-Nσ*pth, Nσ*pth]，其中 pth = sqrt(me*kT)
    // 背景离子：px,py,pz in [-Nσ*pth, Nσ*pth]，其中 pth = sqrt(mi*kT_i)
    const double Nsigma = 8.0;  // 动量网格覆盖的热宽度倍数

    // WENO 幽灵单元数
    const int Nghost = 3;
}

#endif // 头文件保护：PARAMETERS_H
