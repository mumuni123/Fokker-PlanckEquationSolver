#include "maxwell.h"
#include "species.h"
#include <cmath>
#include <cstring>
#include <mpi.h>

void EMFields::init(const SpatialGrid& sg)
{
    nx_total = sg.nx_total;
    dx = sg.dx;

    Ex.assign(nx_total, 0.0);
    Ey.assign(nx_total, 0.0);
    Ez.assign(nx_total, 0.0);
    By.assign(nx_total, 0.0);  // By_{i+1/2}，最后一个元素未使用
    Bz.assign(nx_total, 0.0);

    rho.assign(nx_total, 0.0);
    Jx.assign(nx_total, 0.0);
    Jy.assign(nx_total, 0.0);
    Jz.assign(nx_total, 0.0);
}

void EMFields::zero_currents()
{
    std::fill(rho.begin(), rho.end(), 0.0);
    std::fill(Jx.begin(), Jx.end(), 0.0);
    std::fill(Jy.begin(), Jy.end(), 0.0);
    std::fill(Jz.begin(), Jz.end(), 0.0);
}

void EMFields::accumulate_moments(const Species& sp)
{
    int ng = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;

    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        rho[ix_g] += sp.charge_density[ix];
        Jx[ix_g]  += sp.current_x[ix];
        Jy[ix_g]  += sp.current_y[ix];
        Jz[ix_g]  += sp.current_z[ix];
    }
}

// ============================================================================
// 更新 B：B^{n+1/2} = B^{n-1/2} + dt * curl(E^n)
// By_{i+1/2}^{n+1/2} = By_{i+1/2}^{n-1/2} + (dt/dx) * (Ez_{i+1} - Ez_i)
// Bz_{i+1/2}^{n+1/2} = Bz_{i+1/2}^{n-1/2} - (dt/dx) * (Ey_{i+1} - Ey_i)
// ============================================================================
void EMFields::update_B(double dt)
{
    double dtdx = dt / dx;
    // 更新内部与幽灵单元（半整数场的最后一个索引除外）
    for (int i = 0; i < nx_total - 1; ++i) {
        By[i] += dtdx * (Ez[i+1] - Ez[i]);
        Bz[i] -= dtdx * (Ey[i+1] - Ey[i]);
    }
}

// ============================================================================
// 更新 E：E^{n+1} = E^{n} + dt * (c^2 * curl(B^{n+1/2}) - J^{n+1/2}/eps0)
// Ex_i^{n+1} = Ex_i^{n} - dt * Jx_i / eps0
// Ey_i^{n+1} = Ey_i^{n} - c^2*(dt/dx)*(Bz_{i+1/2} - Bz_{i-1/2}) - dt*Jy_i/eps0
// Ez_i^{n+1} = Ez_i^{n} + c^2*(dt/dx)*(By_{i+1/2} - By_{i-1/2}) - dt*Jz_i/eps0
// ============================================================================
void EMFields::update_E(double dt)
{
    double c2 = Const::c * Const::c;
    double c2dtdx = c2 * dt / dx;
    double dt_eps = dt / Const::eps0;

    // Ex 更新（无空间导数，仅有电流源项）
    for (int i = 1; i < nx_total - 1; ++i) {
        Ex[i] -= dt_eps * Jx[i];
    }

    // Ey、Ez 更新（包含 curl B）
    for (int i = 1; i < nx_total - 1; ++i) {
        Ey[i] -= c2dtdx * (Bz[i] - Bz[i-1]) + dt_eps * Jy[i];
        Ez[i] += c2dtdx * (By[i] - By[i-1]) - dt_eps * Jz[i];
    }
}

// ============================================================================
// Ex 的高斯定律修正：∂Ex/∂x = ρ/ε₀
// 从左边界积分，得到与电荷密度一致的 Ex。
// 该步骤用于修正数值散度误差。
// ============================================================================
void EMFields::gauss_correction()
{
    // 简单积分：Ex[i+1] = Ex[i] + dx * rho[i] / eps0
    // 从最左侧内部单元开始
    int ng = Param::Nghost;
    for (int i = ng; i < nx_total - ng - 1; ++i) {
        Ex[i+1] = Ex[i] + dx * rho[i] / Const::eps0;
    }
}

// ============================================================================
// 电磁场幽灵单元交换（MPI 周期边界）
// ============================================================================
void EMFields::exchange_ghosts(int mpi_rank, int mpi_size)
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;

    int left  = (mpi_rank - 1 + mpi_size) % mpi_size;
    int right = (mpi_rank + 1) % mpi_size;

    // 每个场分量分别交换
    // 打包格式：ng 个单元的 [Ex, Ey, Ez, By, Bz]
    int buf_size = 5 * ng;
    std::vector<double> send_l(buf_size), send_r(buf_size);
    std::vector<double> recv_l(buf_size), recv_r(buf_size);

    // 打包左发送缓冲区（前 ng 个内部单元）
    for (int g = 0; g < ng; ++g) {
        int i = ng + g;
        send_l[5*g+0] = Ex[i];
        send_l[5*g+1] = Ey[i];
        send_l[5*g+2] = Ez[i];
        send_l[5*g+3] = By[i];
        send_l[5*g+4] = Bz[i];
    }
    // 打包右发送缓冲区（后 ng 个内部单元）
    for (int g = 0; g < ng; ++g) {
        int i = ng + nxl - ng + g;
        send_r[5*g+0] = Ex[i];
        send_r[5*g+1] = Ey[i];
        send_r[5*g+2] = Ez[i];
        send_r[5*g+3] = By[i];
        send_r[5*g+4] = Bz[i];
    }

    MPI_Request reqs[4];
    MPI_Isend(send_l.data(), buf_size, MPI_DOUBLE, left,  201, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(send_r.data(), buf_size, MPI_DOUBLE, right, 202, MPI_COMM_WORLD, &reqs[1]);
    MPI_Irecv(recv_l.data(), buf_size, MPI_DOUBLE, left,  202, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(recv_r.data(), buf_size, MPI_DOUBLE, right, 201, MPI_COMM_WORLD, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    // 解包左侧幽灵单元
    for (int g = 0; g < ng; ++g) {
        Ex[g] = recv_l[5*g+0];
        Ey[g] = recv_l[5*g+1];
        Ez[g] = recv_l[5*g+2];
        By[g] = recv_l[5*g+3];
        Bz[g] = recv_l[5*g+4];
    }
    // 解包右侧幽灵单元
    for (int g = 0; g < ng; ++g) {
        int i = ng + nxl + g;
        Ex[i] = recv_r[5*g+0];
        Ey[i] = recv_r[5*g+1];
        Ez[i] = recv_r[5*g+2];
        By[i] = recv_r[5*g+3];
        Bz[i] = recv_r[5*g+4];
    }
}

double EMFields::total_energy() const
{
    int ng = Param::Nghost;
    int nxl = nx_total - 2 * ng;
    double energy = 0.0;
    double c2 = Const::c * Const::c;

    for (int i = ng; i < ng + nxl; ++i) {
        // 电场能量：eps0/2 * (Ex^2 + Ey^2 + Ez^2)
        energy += 0.5 * Const::eps0 * (Ex[i]*Ex[i] + Ey[i]*Ey[i] + Ez[i]*Ez[i]);
        // 磁场能量：1/(2*mu0) * (By^2 + Bz^2)，定义在半整数节点
        if (i < ng + nxl - 1) {
            energy += 0.5 / Const::mu0 * (By[i]*By[i] + Bz[i]*Bz[i]);
        }
    }
    return energy * dx;  // 沿 x 积分
}
