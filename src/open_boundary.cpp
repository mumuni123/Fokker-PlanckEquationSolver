#include "open_boundary.h"

#include <algorithm>
#include <mpi.h>

namespace {
void exchange_internal_ghosts(Species& electrons, const SpatialGrid& grid,
                              int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int slab = ng * static_cast<int>(Param::Nvmu);
    if (slab == 0 || size <= 1) return;

    const int left = rank > 0 ? rank - 1 : MPI_PROC_NULL;
    const int right = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
    double* data = electrons.f.data();
    MPI_Sendrecv(data + static_cast<size_t>(ng) * Param::Nvmu, slab, MPI_DOUBLE,
                 left, 7101,
                 data + static_cast<size_t>(ng + nxl) * Param::Nvmu, slab,
                 MPI_DOUBLE, right, 7101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(data + static_cast<size_t>(ng + nxl - ng) * Param::Nvmu, slab,
                 MPI_DOUBLE, right, 7102,
                 data, slab, MPI_DOUBLE, left, 7102,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void exchange_periodic_ghosts(Species& electrons, const SpatialGrid& grid,
                              int rank, int size)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int slab = ng * static_cast<int>(Param::Nvmu);
    if (slab == 0) return;
    const int left = rank > 0 ? rank - 1 : size - 1;
    const int right = rank + 1 < size ? rank + 1 : 0;
    double* data = electrons.f.data();
    MPI_Sendrecv(data + static_cast<size_t>(ng) * Param::Nvmu, slab,
                 MPI_DOUBLE, left, 7111,
                 data + static_cast<size_t>(ng + nxl) * Param::Nvmu, slab,
                 MPI_DOUBLE, right, 7111, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(data + static_cast<size_t>(ng + nxl - ng) * Param::Nvmu,
                 slab, MPI_DOUBLE, right, 7112,
                 data, slab, MPI_DOUBLE, left, 7112,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void copy_velocity_slice(std::vector<double>& dst, int dst_ix,
                         const std::vector<double>& src, int src_ix)
{
    const size_t count = Param::Nvmu;
    std::copy(src.begin() + static_cast<size_t>(src_ix) * count,
              src.begin() + static_cast<size_t>(src_ix + 1) * count,
              dst.begin() + static_cast<size_t>(dst_ix) * count);
}
}

OpenBackgroundBoundary::OpenBackgroundBoundary()
{
    config_.left_type = BackgroundXBoundaryType::RESERVOIR;
    config_.right_type = BackgroundXBoundaryType::RESERVOIR;
    config_.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    config_.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    reservoir_prepared_ = false;
    prepared_dx_ = 0.0;
}

OpenBackgroundBoundary::OpenBackgroundBoundary(const OpenBackgroundBoundaryConfig& config)
    : config_(config)
{
    reservoir_prepared_ = false;
    prepared_dx_ = 0.0;
}

void OpenBackgroundBoundary::ensure_reservoir_prepared(
    const Species& electrons) const
{
    if (!electrons.sgrid) return;
    const double dx = electrons.sgrid->dx;
    if (reservoir_prepared_ && prepared_dx_ == dx) return;

    std::vector<double> slice;
    if (config_.left_type == BackgroundXBoundaryType::RESERVOIR) {
        electrons.fill_cylindrical_maxwellian_mass_slice(
            slice, config_.left_reservoir.density,
            config_.left_reservoir.temperature, config_.left_reservoir.drift_vx);
        left_reservoir_line_density_.assign(Param::Nvmu, 0.0);
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            left_reservoir_line_density_[q] = slice[q] / dx;
        }
    } else {
        left_reservoir_line_density_.assign(Param::Nvmu, 0.0);
    }
    if (config_.right_type == BackgroundXBoundaryType::RESERVOIR) {
        electrons.fill_cylindrical_maxwellian_mass_slice(
            slice, config_.right_reservoir.density,
            config_.right_reservoir.temperature, config_.right_reservoir.drift_vx);
        right_reservoir_line_density_.assign(Param::Nvmu, 0.0);
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            right_reservoir_line_density_[q] = slice[q] / dx;
        }
    } else {
        right_reservoir_line_density_.assign(Param::Nvmu, 0.0);
    }
    prepared_dx_ = dx;
    reservoir_prepared_ = true;
}

double OpenBackgroundBoundary::incoming_cell_average(
    PhysicalSide side, int j_upar, int k_uperp, double time,
    const Species& electrons) const
{
    (void)time;
    const bool reservoir = (side == PhysicalSide::LEFT)
        ? config_.left_type == BackgroundXBoundaryType::RESERVOIR
        : config_.right_type == BackgroundXBoundaryType::RESERVOIR;
    if (!reservoir) return 0.0;
    ensure_reservoir_prepared(electrons);
    const size_t q = idx2(j_upar, k_uperp);
    return (side == PhysicalSide::LEFT)
               ? left_reservoir_line_density_[q]
               : right_reservoir_line_density_[q];
}

bool OpenBackgroundBoundary::is_incoming(PhysicalSide side, double vx) const
{
    return (side == PhysicalSide::LEFT) ? vx > 0.0 : vx < 0.0;
}

void OpenBackgroundBoundary::fill_ghosts(Species& electrons,
                                          const SpatialGrid& grid,
                                          int mpi_rank, int mpi_size) const
{
    if (config_.left_type == BackgroundXBoundaryType::PERIODIC &&
        config_.right_type == BackgroundXBoundaryType::PERIODIC) {
        exchange_periodic_ghosts(electrons, grid, mpi_rank, mpi_size);
        return;
    }
    exchange_internal_ghosts(electrons, grid, mpi_rank, mpi_size);
    if (grid.owns_left_physical_boundary(mpi_rank)) fill_physical_left(electrons, grid);
    if (grid.owns_right_physical_boundary(mpi_rank, mpi_size)) fill_physical_right(electrons, grid);
}

void OpenBackgroundBoundary::fill_physical_left(Species& electrons,
                                                 const SpatialGrid& grid) const
{
    ensure_reservoir_prepared(electrons);
    const int interior = grid.nghost;
    for (int g = 0; g < grid.nghost; ++g) {
        const int ghost_ix = grid.nghost - 1 - g;
        copy_velocity_slice(electrons.f, ghost_ix, electrons.f, interior);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t q = idx2(iv, imu);
                if (electrons.cgrid.vx[q] > 0.0) {
                    electrons.f[idx3(ghost_ix, iv, imu)] =
                        left_reservoir_line_density_[q] * grid.dx;
                }
            }
        }
    }
}

void OpenBackgroundBoundary::fill_physical_right(Species& electrons,
                                                  const SpatialGrid& grid) const
{
    ensure_reservoir_prepared(electrons);
    const int interior = grid.nghost + grid.nx_local - 1;
    for (int g = 0; g < grid.nghost; ++g) {
        const int ghost_ix = grid.nghost + grid.nx_local + g;
        copy_velocity_slice(electrons.f, ghost_ix, electrons.f, interior);
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t q = idx2(iv, imu);
                if (electrons.cgrid.vx[q] < 0.0) {
                    electrons.f[idx3(ghost_ix, iv, imu)] =
                        right_reservoir_line_density_[q] * grid.dx;
                }
            }
        }
    }
}
