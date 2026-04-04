#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "species.h"
#include "maxwell.h"
#include "grid.h"
#include <string>
#include <vector>
#include <fstream>

class Diagnostics {
public:
    std::string output_dir;
    int snapshot_count;

    Diagnostics();

    void init(const std::string& dir, int mpi_rank);

    // 将标量诊断量（能量、粒子数）写入 ASCII 文件
    void write_scalars(double time, int step,
                       const std::vector<Species*>& species,
                       const EMFields& fields,
                       int mpi_rank, int mpi_size);

    // 写出场快照（沿 x 的 Ex, Ey, Ez, By, Bz, rho, Jx）
    void write_fields(double time,
                      const EMFields& fields,
                      const SpatialGrid& sg,
                      int mpi_rank, int mpi_size);

    // 写出 1D 动量分布（对 x、py、pz 积分得到 f(px)）
    void write_px_distribution(double time,
                               const Species& sp,
                               int mpi_rank, int mpi_size);

    // 写出数密度剖面 n(x)
    void write_density_profile(double time,
                               const std::vector<Species*>& species,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size);

private:
    std::ofstream scalar_file;
};

#endif // 头文件保护：DIAGNOSTICS_H
