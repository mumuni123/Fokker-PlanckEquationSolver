#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "vlasov.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

// ============================================================================
// 基于 CFL 条件计算全局时间步长
// WENO-5 + Euler 的 CFL 经验值约为 0.4
// ============================================================================
double compute_dt(const std::vector<Species*>& species, const SpatialGrid& sg)
{
    double dt_min = 1.0e10;

    for (size_t s = 0; s < species.size(); ++s) {
        const Species& sp = *species[s];

        // x 向对流 CFL：dt < CFL * dx / max(|vx|)
        // 对相对论粒子，max|vx| = c
        double dt_x = 0.4 * sg.dx / Const::c;
        dt_min = std::min(dt_min, dt_x);

        // p 向对流 CFL：dt < CFL * dp / max(|F|)
        // 由典型场强估计最大受力
        // 初始阶段采用保守估计，后续仿真中可再细化
        // 这里先用束流驱动回流电流的估算：
        // E ~ n_b * e * dx / eps0（非常粗略），但偏保守
        // 因此按热能尺度的分量给出估计
        double dp_min = std::min(sp.pgrid.dpx, std::min(sp.pgrid.dpy, sp.pgrid.dpz));
        // 粗略估计：F_max ~ e * E_max，E_max ~ dens * e / eps0 * dx
        // 该估计偏保守，因此再给一个相对宽松上界
        double E_est = Param::densb * Const::qe * sg.dx / Const::eps0;
        double F_max = Const::qe * E_est;
        if (F_max > 1.0e-30) {
            double dt_p = 0.4 * dp_min / F_max;
            dt_min = std::min(dt_min, dt_p);
        }
    }

    // 施加安全系数
    dt_min *= Param::dt_multiplier;

    // 限制到合理上限
    double dt_max = 0.01 * Const::femto;
    if (dt_min > dt_max) dt_min = dt_max;

    return dt_min;
}

// ============================================================================
// 束流注入：在 x = 0 处加入漂移（冷）束分布
// 注入导致密度变化率：dn/dt = n_beam * v_beam
// 因而每个时间步：Δf = n_beam * v_beam * dt / dx * delta(p - p_drift)
// 在离散网格上，delta 函数由束单元内的 1/dp³ 表示
// ============================================================================
void inject_beam(Species& beam, const SpatialGrid& sg, double dt, double time)
{
    if (time > Param::t_inject_end) return;

    // 仅在 x = 0 注入（rank 0 的第一个空间单元，
    // 即全局索引对应 x=0 的单元）
    int ng = sg.nghost;

    // 找到与首个空间单元（全局 ix=0）对应的本地单元
    // 在 MPI 分解下，由拥有 ix_start=0 的 rank 负责注入
    if (sg.ix_start != 0) return;

    int ix_inject = ng;  // 第一个内部单元

    // 找到最接近束流漂移动量的动量单元
    double px_drift = Param::beam_drift_px;
    int ipx_drift = static_cast<int>((px_drift - beam.pgrid.px_min) / beam.pgrid.dpx);
    if (ipx_drift < 0 || ipx_drift >= Param::Npx) return;

    // py=0, pz=0 对应的单元
    int ipy_center = Param::Npy / 2;
    int ipz_center = Param::Npz / 2;

    // 注入率：单位面积单位时间粒子数 = n_beam * v_beam
    // f 的单位为 [1/m³/(kg*m/s)³]
    // 每步 Δf：inject_rate * dt / dx / dp³
    // 其中 inject_rate [1/m² /s] = densb * betab * c
    double v_beam = Param::betab * Const::c;
    double inject_rate = Param::densb * v_beam;

    // 不直接使用 delta 函数，而在动量空间使用窄高斯
    // 宽度约 1 个网格单元，以提高数值平滑性
    double sigma_px = beam.pgrid.dpx;
    double sigma_py = beam.pgrid.dpy;
    double sigma_pz = beam.pgrid.dpz;
    double norm = inject_rate * dt / sg.dx;
    double norm_p = 1.0 / (std::pow(2.0 * Const::pi, 1.5) * sigma_px * sigma_py * sigma_pz);

    for (int ipx = 0; ipx < Param::Npx; ++ipx) {
        double ppx = beam.pgrid.px(ipx);
        double dpx2 = (ppx - px_drift) * (ppx - px_drift) / (2.0 * sigma_px * sigma_px);
        if (dpx2 > 20.0) continue;  // 跳过可忽略贡献

        for (int ipy = 0; ipy < Param::Npy; ++ipy) {
            double ppy = beam.pgrid.py(ipy);
            double dpy2 = ppy * ppy / (2.0 * sigma_py * sigma_py);
            if (dpx2 + dpy2 > 20.0) continue;

            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = beam.pgrid.pz(ipz);
                double dpz2 = ppz * ppz / (2.0 * sigma_pz * sigma_pz);
                double exponent = dpx2 + dpy2 + dpz2;
                if (exponent > 20.0) continue;

                double df = norm * norm_p * std::exp(-exponent);
                beam.f[idx4(ix_inject, ipx, ipy, ipz)] += df;
            }
        }
    }
}

// ============================================================================
// 主程序
// ============================================================================
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  1D3V Relativistic Vlasov-Fokker-Planck Solver\n");
        printf("  Matching EPOCH input.deck parameters\n");
        printf("============================================================\n");
        printf("MPI ranks: %d\n", mpi_size);
        #pragma omp parallel
        {
            #pragma omp master
            printf("OpenMP threads per rank: %d\n", omp_get_num_threads());
        }
        printf("Spatial grid: nx = %d, dx = %.3e m\n", Param::nx, Param::dx);
        printf("Momentum grid: %d x %d x %d = %lu cells\n",
               Param::Npx, Param::Npy, Param::Npz,
               (unsigned long)Param::Np3);
        printf("Background: n_e = %.3e /m^3, T_e = %.1f eV\n",
               Param::dens, Param::temperature_e / Const::eV);
        printf("Beam: gamma*beta = %.2f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::densb);
        printf("t_end = %.1f fs, injection duration = %.1f fs\n",
               Param::t_end / Const::femto,
               Param::t_inject_end / Const::femto);
        printf("------------------------------------------------------------\n");
    }

    // --- 初始化空间网格 ---
    SpatialGrid sgrid;
    sgrid.init(mpi_rank, mpi_size);

    if (mpi_rank == 0) {
        printf("Rank 0: nx_local = %d, ix_start = %d\n",
               sgrid.nx_local, sgrid.ix_start);
    }

    // --- 初始化物种 ---
    Species beam, bkg_e, bkg_i;

    beam.init("beam", SpeciesType::BEAM,
              -Const::qe, Const::me,
              Param::densb, 0.0,
              false,  // 束流碰撞关闭
              sgrid);

    bkg_e.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
               -Const::qe, Const::me,
               Param::dens, Param::temperature_e,
               true,   // 碰撞开启
               sgrid);

    bkg_i.init("bkg_i", SpeciesType::BACKGROUND_ION,
               Param::Z_ion * Const::qe, Param::mass_ion,
               Param::dens / Param::Z_ion, Param::temperature_i,
               true,
               sgrid);

    // 将背景物种初始化为麦克斯韦分布
    bkg_e.initialize_maxwellian(0.0);
    bkg_i.initialize_maxwellian(0.0);

    // 内存占用报告
    {
        double mem_per_species = (double)beam.local_size() * 8.0 / (1024.0*1024.0);
        if (mpi_rank == 0) {
            printf("Memory per species per rank: %.1f MB\n", mem_per_species);
            printf("Total memory estimate: %.1f MB per rank\n", mem_per_species * 3 * 3);
        }
    }

    // 为方便处理，保存物种指针数组
    std::vector<Species*> all_species;
    all_species.push_back(&beam);
    all_species.push_back(&bkg_e);
    all_species.push_back(&bkg_i);

    // --- 初始化电磁场 ---
    EMFields fields;
    fields.init(sgrid);

    // --- 求解器对象 ---
    VlasovSolver vlasov;
    CollisionOperator collision;
    Diagnostics diag;
    diag.init("output", mpi_rank);

    // --- 计算时间步长 ---
    double dt = compute_dt(all_species, sgrid);
    // 在所有 rank 间同步 dt
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    int nsteps = static_cast<int>(std::ceil(Param::t_end / dt));
    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    // --- t=0 时刻输出诊断 ---
    for (size_t s = 0; s < all_species.size(); ++s) {
        all_species[s]->compute_moments();
    }
    diag.write_scalars(0.0, 0, all_species, fields, mpi_rank, mpi_size);
    diag.write_density_profile(0.0, all_species, sgrid, mpi_rank, mpi_size);

    // 快照输出控制
    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 200;

    // ================================================================
    // 主时间循环：Strang 分裂
    //
    //   1. 半步 Vlasov 对流 (dt/2)
    //   2. 计算矩并更新电磁场 (dt)
    //   3. 一整步 FP 碰撞 (dt) [仅背景物种]
    //   4. 半步 Vlasov 对流 (dt/2)
    //   5. 束流注入
    // ================================================================
    for (int step = 1; step <= nsteps; ++step) {
        double time = step * dt;

        // ------ 步骤 1：半步 Vlasov 对流 ------
        for (size_t s = 0; s < all_species.size(); ++s) {
            vlasov.advect(*all_species[s], sgrid, fields, 0.5 * dt,
                          mpi_rank, mpi_size);
        }

        // ------ 步骤 2：计算矩并更新电磁场 ------
        // 从所有物种计算电荷密度与电流
        fields.zero_currents();
        for (size_t s = 0; s < all_species.size(); ++s) {
            all_species[s]->compute_moments();
            fields.accumulate_moments(*all_species[s]);
        }

        // 交换场的幽灵单元
        fields.exchange_ghosts(mpi_rank, mpi_size);

        // Maxwell 更新：先 B^{n-1/2} -> B^{n+1/2}，再 E^n -> E^{n+1}
        fields.update_B(dt);
        fields.update_E(dt);

        // 周期性高斯修正
        if (step % 10 == 0) {
            fields.gauss_correction();
        }

        fields.exchange_ghosts(mpi_rank, mpi_size);

        // ------ 步骤 3：整步 FP 碰撞（仅背景物种）------
        // bkg_e 与 bkg_e（自碰撞）和 bkg_i 发生碰撞
        if (bkg_e.collisions_enabled) {
            // 自碰撞：e-e
            collision.apply(bkg_e, dt,
                            Param::dens, Param::temperature_e, Const::me,
                            1.0, 1.0);  // Z_field=1，Z_test=1
            // e-i 碰撞
            collision.apply(bkg_e, dt,
                            Param::dens / Param::Z_ion, Param::temperature_i,
                            Param::mass_ion,
                            (double)Param::Z_ion, 1.0);
        }

        // bkg_i 与 bkg_i（自碰撞）和 bkg_e 发生碰撞
        if (bkg_i.collisions_enabled) {
            // 自碰撞：i-i
            collision.apply(bkg_i, dt,
                            Param::dens / Param::Z_ion, Param::temperature_i,
                            Param::mass_ion,
                            (double)Param::Z_ion, (double)Param::Z_ion);
            // i-e 碰撞
            collision.apply(bkg_i, dt,
                            Param::dens, Param::temperature_e, Const::me,
                            1.0, (double)Param::Z_ion);
        }

        // ------ 步骤 4：半步 Vlasov 对流 ------
        for (size_t s = 0; s < all_species.size(); ++s) {
            vlasov.advect(*all_species[s], sgrid, fields, 0.5 * dt,
                          mpi_rank, mpi_size);
        }

        // ------ 步骤 5：束流注入 ------
        inject_beam(beam, sgrid, dt, time);

        // ------ 诊断输出 ------
        if (step % stdout_freq == 0) {
            for (size_t s = 0; s < all_species.size(); ++s) {
                all_species[s]->compute_moments();
            }
            diag.write_scalars(time, step, all_species, fields, mpi_rank, mpi_size);

            if (mpi_rank == 0) {
                printf("Step %d / %d, t = %.4f fs\n",
                       step, nsteps, time / Const::femto);
            }
        }

        // 快照输出
        if (time >= next_snapshot) {
            for (size_t s = 0; s < all_species.size(); ++s) {
                all_species[s]->compute_moments();
            }
            diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
            diag.write_density_profile(time, all_species, sgrid, mpi_rank, mpi_size);

            for (size_t s = 0; s < all_species.size(); ++s) {
                diag.write_px_distribution(time, *all_species[s], mpi_rank, mpi_size);
            }

            if (mpi_rank == 0) {
                printf("  >> Snapshot %d written at t = %.4f fs\n",
                       diag.snapshot_count, time / Const::femto);
            }
            next_snapshot += Param::dt_snapshot;
        }
    }

    // 最终输出
    for (size_t s = 0; s < all_species.size(); ++s) {
        all_species[s]->compute_moments();
    }
    diag.write_scalars(Param::t_end, nsteps, all_species, fields, mpi_rank, mpi_size);
    diag.write_density_profile(Param::t_end, all_species, sgrid, mpi_rank, mpi_size);

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Simulation complete: t = %.1f fs, %d steps\n",
               Param::t_end / Const::femto, nsteps);
        printf("============================================================\n");
    }

    MPI_Finalize();
    return 0;
}
