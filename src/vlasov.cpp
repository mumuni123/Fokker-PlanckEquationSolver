#include "vlasov.h"
#include "maxwell.h"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <omp.h>
#include <mpi.h>

// ============================================================================
// 顶层对流推进：Strang 分裂维度扫掠
// 顺序：x -> px -> py -> pz（由主循环以半步调用）
// ============================================================================
void VlasovSolver::advect(Species& sp, const SpatialGrid& sg,
                          const EMFields& fields, double dt,
                          int mpi_rank, int mpi_size)
{
    advect_x(sp, sg, dt, mpi_rank, mpi_size);
    advect_px(sp, sg, fields.Ex, fields.Ey, fields.Ez, fields.By, fields.Bz, dt);
    advect_py(sp, sg, fields.Ex, fields.Ey, fields.Ez, fields.By, fields.Bz, dt);
    advect_pz(sp, sg, fields.Ex, fields.Ey, fields.Ez, fields.By, fields.Bz, dt);
}

// ============================================================================
// x 方向对流：∂f/∂t + vx(p) ∂f/∂x = 0
// 对每个动量单元 (ipx, ipy, ipz)，以速度 vx 沿 x 扫掠
// 使用 WENO-5 + 前向 Euler（主循环中的单子步）
// ============================================================================
void VlasovSolver::advect_x(Species& sp, const SpatialGrid& sg, double dt,
                            int mpi_rank, int mpi_size)
{
    // 先交换 x 方向幽灵单元
    exchange_ghosts_x(sp, sg, mpi_rank, mpi_size);

    int nxt = sg.nx_total;
    int ng  = sg.nghost;
    int nxl = sg.nx_local;
    double dxi = 1.0 / sg.dx;

    // 存放更新后 f 的临时数组
    std::vector<double> f_new(sp.f.size());
    std::copy(sp.f.begin(), sp.f.end(), f_new.begin());

    // 对每个动量单元执行 x 向 1D 对流
    #pragma omp parallel for collapse(2)
    for (int ipx = 0; ipx < Param::Npx; ++ipx) {
        for (int ipy = 0; ipy < Param::Npy; ++ipy) {
            double ppx = sp.pgrid.px(ipx);
            double ppy = sp.pgrid.py(ipy);

            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = sp.pgrid.pz(ipz);
                double vel = vx_from_p(ppx, ppy, ppz, sp.mass);

                // 提取该动量单元对应的 x 向 1D 切片
                // 包含幽灵单元
                for (int ix = ng; ix < ng + nxl; ++ix) {
                    // WENO-5 通量差分
                    double qm2 = sp.f[idx4(ix-2, ipx, ipy, ipz)];
                    double qm1 = sp.f[idx4(ix-1, ipx, ipy, ipz)];
                    double q0  = sp.f[idx4(ix,   ipx, ipy, ipz)];
                    double qp1 = sp.f[idx4(ix+1, ipx, ipy, ipz)];
                    double qp2 = sp.f[idx4(ix+2, ipx, ipy, ipz)];
                    double qp3 = (ix+3 < nxt) ? sp.f[idx4(ix+3, ipx, ipy, ipz)] : 0.0;
                    double qm3 = (ix-3 >= 0)  ? sp.f[idx4(ix-3, ipx, ipy, ipz)] : 0.0;

                    // 计算 ix-1/2 与 ix+1/2 两个面的通量
                    double flux_right, flux_left;

                    if (vel >= 0.0) {
                        // 从左侧迎风
                        // 面 ix+1/2：左偏，模板为 ix-2..ix+2
                        flux_right = vel * weno5_left(qm2, qm1, q0, qp1, qp2);
                        // 面 ix-1/2：左偏，模板为 ix-3..ix+1
                        flux_left = vel * weno5_left(qm3, qm2, qm1, q0, qp1);
                    } else {
                        // 从右侧迎风
                        // 面 ix+1/2：右偏
                        flux_right = vel * weno5_right(qm1, q0, qp1, qp2, qp3);
                        // 面 ix-1/2：右偏，使用平移后的模板
                        flux_left = vel * weno5_right(qm2, qm1, q0, qp1, qp2);
                    }

                    // 有限体积更新：df/dt = -(F_{i+1/2} - F_{i-1/2}) / dx
                    f_new[idx4(ix, ipx, ipy, ipz)] = q0 - dt * dxi * (flux_right - flux_left);

                    // 正定性约束
                    if (f_new[idx4(ix, ipx, ipy, ipz)] < 0.0) {
                        f_new[idx4(ix, ipx, ipy, ipz)] = 0.0;
                    }
                }
            }
        }
    }

    // 回写结果
    sp.f.swap(f_new);
}

// ============================================================================
// px 方向对流：∂f/∂t + Fx ∂f/∂px = 0
// 力项：Fx = q*(Ex + vy*Bz - vz*By)
// 对每个 (ix, ipy, ipz) 沿 px 方向按 Fx 扫掠
// ============================================================================
void VlasovSolver::advect_px(Species& sp, const SpatialGrid& sg,
                             const std::vector<double>& Ex,
                             const std::vector<double>& Ey,
                             const std::vector<double>& Ez,
                             const std::vector<double>& By,
                             const std::vector<double>& Bz,
                             double dt)
{
    int ng  = sg.nghost;
    int nxl = sg.nx_local;
    double dpxi = 1.0 / sp.pgrid.dpx;

    std::vector<double> f_new(sp.f.size());
    std::copy(sp.f.begin(), sp.f.end(), f_new.begin());

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double ex = Ex[ix_g];
        double ey = Ey[ix_g];
        double ez = Ez[ix_g];
        // By、Bz 定义在半节点，取平均到单元中心
        double by = 0.5 * (By[ix_g] + By[ix_g + 1]);
        double bz = 0.5 * (Bz[ix_g] + Bz[ix_g + 1]);

        for (int ipy = 0; ipy < Param::Npy; ++ipy) {
            double ppy = sp.pgrid.py(ipy);
            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = sp.pgrid.pz(ipz);

                // 沿 px 扫掠
                for (int ipx = 0; ipx < Param::Npx; ++ipx) {
                    double ppx = sp.pgrid.px(ipx);
                    double gam = gamma_from_p(ppx, ppy, ppz, sp.mass);
                    double inv_gm = 1.0 / (gam * sp.mass);
                    double vy_loc = ppy * inv_gm;
                    double vz_loc = ppz * inv_gm;

                    // 洛伦兹力的 px 分量
                    double Fx = sp.charge * (ex + vy_loc * bz - vz_loc * by);

                    // 沿 px 收集模板（边界使用零填充）
                    double q[Param::Npx + 6]; // 两侧各补 3 个幽灵单元
                    for (int s = -3; s < Param::Npx + 3; ++s) {
                        if (s >= 0 && s < Param::Npx) {
                            q[s + 3] = sp.f[idx4(ix_g, s, ipy, ipz)];
                        } else {
                            q[s + 3] = 0.0; // 动量边界取零
                        }
                    }

                    // ipx +/- 1/2 处的 WENO-5 通量
                    double flux_right, flux_left;
                    int si = ipx + 3; // 在补齐数组中的平移索引

                    if (Fx >= 0.0) {
                        flux_right = Fx * weno5_left(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                        flux_left  = Fx * weno5_left(q[si-3], q[si-2], q[si-1], q[si], q[si+1]);
                    } else {
                        flux_right = Fx * weno5_right(q[si-1], q[si], q[si+1], q[si+2], q[si+3]);
                        flux_left  = Fx * weno5_right(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                    }

                    double fval = q[si] - dt * dpxi * (flux_right - flux_left);
                    if (fval < 0.0) fval = 0.0;
                    f_new[idx4(ix_g, ipx, ipy, ipz)] = fval;
                }
            }
        }
    }
    sp.f.swap(f_new);
}

// ============================================================================
// py 方向对流：∂f/∂t + Fy ∂f/∂py = 0
// 力项：Fy = q*(Ey - vx*Bz)（在 1D 中 Bx=0）
// ============================================================================
void VlasovSolver::advect_py(Species& sp, const SpatialGrid& sg,
                             const std::vector<double>& Ex,
                             const std::vector<double>& Ey,
                             const std::vector<double>& Ez,
                             const std::vector<double>& By,
                             const std::vector<double>& Bz,
                             double dt)
{
    int ng  = sg.nghost;
    int nxl = sg.nx_local;
    double dpyi = 1.0 / sp.pgrid.dpy;

    std::vector<double> f_new(sp.f.size());
    std::copy(sp.f.begin(), sp.f.end(), f_new.begin());

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double ey_val = Ey[ix_g];
        double bz_val = 0.5 * (Bz[ix_g] + Bz[ix_g + 1]);

        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double ppx = sp.pgrid.px(ipx);
            for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                double ppz = sp.pgrid.pz(ipz);

                for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                    double ppy = sp.pgrid.py(ipy);
                    double gam = gamma_from_p(ppx, ppy, ppz, sp.mass);
                    double vx_loc = ppx / (gam * sp.mass);

                    // Fy = q*(Ey - vx*Bz)
                    double Fy = sp.charge * (ey_val - vx_loc * bz_val);

                    // 沿 py 收集模板
                    double q[Param::Npy + 6];
                    for (int s = -3; s < Param::Npy + 3; ++s) {
                        if (s >= 0 && s < Param::Npy) {
                            q[s + 3] = sp.f[idx4(ix_g, ipx, s, ipz)];
                        } else {
                            q[s + 3] = 0.0;
                        }
                    }

                    int si = ipy + 3;
                    double flux_right, flux_left;

                    if (Fy >= 0.0) {
                        flux_right = Fy * weno5_left(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                        flux_left  = Fy * weno5_left(q[si-3], q[si-2], q[si-1], q[si], q[si+1]);
                    } else {
                        flux_right = Fy * weno5_right(q[si-1], q[si], q[si+1], q[si+2], q[si+3]);
                        flux_left  = Fy * weno5_right(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                    }

                    double fval = q[si] - dt * dpyi * (flux_right - flux_left);
                    if (fval < 0.0) fval = 0.0;
                    f_new[idx4(ix_g, ipx, ipy, ipz)] = fval;
                }
            }
        }
    }
    sp.f.swap(f_new);
}

// ============================================================================
// pz 方向对流：∂f/∂t + Fz ∂f/∂pz = 0
// 力项：Fz = q*(Ez + vx*By)
// ============================================================================
void VlasovSolver::advect_pz(Species& sp, const SpatialGrid& sg,
                             const std::vector<double>& Ex,
                             const std::vector<double>& Ey,
                             const std::vector<double>& Ez,
                             const std::vector<double>& By,
                             const std::vector<double>& Bz,
                             double dt)
{
    int ng  = sg.nghost;
    int nxl = sg.nx_local;
    double dpzi = 1.0 / sp.pgrid.dpz;

    std::vector<double> f_new(sp.f.size());
    std::copy(sp.f.begin(), sp.f.end(), f_new.begin());

    #pragma omp parallel for
    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        double ez_val = Ez[ix_g];
        double by_val = 0.5 * (By[ix_g] + By[ix_g + 1]);

        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double ppx = sp.pgrid.px(ipx);
            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                double ppy = sp.pgrid.py(ipy);

                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    double ppz = sp.pgrid.pz(ipz);
                    double gam = gamma_from_p(ppx, ppy, ppz, sp.mass);
                    double vx_loc = ppx / (gam * sp.mass);

                    // Fz = q*(Ez + vx*By)
                    double Fz = sp.charge * (ez_val + vx_loc * by_val);

                    // 沿 pz 收集模板
                    double q[Param::Npz + 6];
                    for (int s = -3; s < Param::Npz + 3; ++s) {
                        if (s >= 0 && s < Param::Npz) {
                            q[s + 3] = sp.f[idx4(ix_g, ipx, ipy, s)];
                        } else {
                            q[s + 3] = 0.0;
                        }
                    }

                    int si = ipz + 3;
                    double flux_right, flux_left;

                    if (Fz >= 0.0) {
                        flux_right = Fz * weno5_left(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                        flux_left  = Fz * weno5_left(q[si-3], q[si-2], q[si-1], q[si], q[si+1]);
                    } else {
                        flux_right = Fz * weno5_right(q[si-1], q[si], q[si+1], q[si+2], q[si+3]);
                        flux_left  = Fz * weno5_right(q[si-2], q[si-1], q[si], q[si+1], q[si+2]);
                    }

                    double fval = q[si] - dt * dpzi * (flux_right - flux_left);
                    if (fval < 0.0) fval = 0.0;
                    f_new[idx4(ix_g, ipx, ipy, ipz)] = fval;
                }
            }
        }
    }
    sp.f.swap(f_new);
}

// ============================================================================
// x 方向 MPI 幽灵单元交换（周期边界）
// ============================================================================
void VlasovSolver::exchange_ghosts_x(Species& sp, const SpatialGrid& sg,
                                     int mpi_rank, int mpi_size)
{
    int ng  = sg.nghost;
    int nxl = sg.nx_local;
    int nxt = sg.nx_total;
    size_t slice_size = Param::Np3;  // 一个 x 切片的动量数据量

    // 确定相邻 rank（周期）
    int left_rank  = (mpi_rank - 1 + mpi_size) % mpi_size;
    int right_rank = (mpi_rank + 1) % mpi_size;

    // 发送/接收缓冲区
    std::vector<double> send_left(ng * slice_size);
    std::vector<double> send_right(ng * slice_size);
    std::vector<double> recv_left(ng * slice_size);
    std::vector<double> recv_right(ng * slice_size);

    // 打包发送缓冲区
    // 向左发送：前 ng 个内部单元 -> 左邻的右侧幽灵区
    for (int g = 0; g < ng; ++g) {
        int ix_src = ng + g;  // 前 ng 个内部单元
        std::memcpy(&send_left[g * slice_size],
                    &sp.f[ix_src * slice_size],
                    slice_size * sizeof(double));
    }
    // 向右发送：后 ng 个内部单元 -> 右邻的左侧幽灵区
    for (int g = 0; g < ng; ++g) {
        int ix_src = ng + nxl - ng + g;
        std::memcpy(&send_right[g * slice_size],
                    &sp.f[ix_src * slice_size],
                    slice_size * sizeof(double));
    }

    MPI_Request reqs[4];
    int tag_lr = 101, tag_rl = 102;

    MPI_Isend(send_left.data(),  (int)(ng * slice_size), MPI_DOUBLE,
              left_rank,  tag_lr, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(send_right.data(), (int)(ng * slice_size), MPI_DOUBLE,
              right_rank, tag_rl, MPI_COMM_WORLD, &reqs[1]);
    MPI_Irecv(recv_left.data(),  (int)(ng * slice_size), MPI_DOUBLE,
              left_rank,  tag_rl, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(recv_right.data(), (int)(ng * slice_size), MPI_DOUBLE,
              right_rank, tag_lr, MPI_COMM_WORLD, &reqs[3]);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // 解包接收缓冲区
    // 左侧幽灵单元（索引 0..ng-1）接收左邻数据
    for (int g = 0; g < ng; ++g) {
        std::memcpy(&sp.f[g * slice_size],
                    &recv_left[g * slice_size],
                    slice_size * sizeof(double));
    }
    // 右侧幽灵单元（索引 ng+nxl..nxt-1）接收右邻数据
    for (int g = 0; g < ng; ++g) {
        int ix_dst = ng + nxl + g;
        std::memcpy(&sp.f[ix_dst * slice_size],
                    &recv_right[g * slice_size],
                    slice_size * sizeof(double));
    }
}
