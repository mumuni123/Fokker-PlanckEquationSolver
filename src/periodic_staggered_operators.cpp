#include "periodic_staggered_operators.h"

#include <algorithm>
#include <mpi.h>

namespace PeriodicStaggered {

void close_right_face_alias(std::vector<double>& face, int nxl, int block,
                            int mpi_rank, int mpi_size, int message_tag)
{
    if (nxl <= 0 || block <= 0 ||
        face.size() < static_cast<size_t>(nxl + 1) * block) return;
    const size_t first = 0;
    const size_t right_alias = static_cast<size_t>(nxl) * block;
    if (mpi_size == 1) {
        std::copy(face.begin() + first, face.begin() + first + block,
                  face.begin() + right_alias);
        return;
    }
    const int left = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right = (mpi_rank + 1) % mpi_size;
    MPI_Sendrecv(face.data() + first, block, MPI_DOUBLE, left, message_tag,
                 face.data() + right_alias, block, MPI_DOUBLE, right, message_tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void apply_face_to_cell_G(const std::vector<double>& face,
                          std::vector<double>& cell, int nxl)
{
    cell.assign(static_cast<size_t>(std::max(0, nxl)), 0.0);
    apply_face_to_cell_G(face, cell.data(), nxl);
}

void apply_face_to_cell_G(const std::vector<double>& face,
                          double* cell, int nxl)
{
    if (!cell || nxl <= 0 || face.size() < static_cast<size_t>(nxl + 1)) return;
    for (int ix = 0; ix < nxl; ++ix)
        cell[ix] = 0.5 * (
            face[static_cast<size_t>(ix)] +
            face[static_cast<size_t>(ix + 1)]);
}

void apply_cell_to_face_Gstar(const std::vector<double>& cell,
                              std::vector<double>& face, int nxl,
                              int mpi_rank, int mpi_size, int message_tag)
{
    face.assign(static_cast<size_t>(std::max(0, nxl) + 1), 0.0);
    if (nxl <= 0 || cell.size() < static_cast<size_t>(nxl)) return;
    double left_neighbor = 0.0;
    if (mpi_size == 1) {
        left_neighbor = cell[static_cast<size_t>(nxl - 1)];
    } else {
        const int left = (mpi_rank + mpi_size - 1) % mpi_size;
        const int right = (mpi_rank + 1) % mpi_size;
        const double last = cell[static_cast<size_t>(nxl - 1)];
        MPI_Sendrecv(&last, 1, MPI_DOUBLE, right, message_tag,
                     &left_neighbor, 1, MPI_DOUBLE, left, message_tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (int iface = 0; iface < nxl; ++iface) {
        const double left_cell = (iface == 0)
            ? left_neighbor : cell[static_cast<size_t>(iface - 1)];
        face[static_cast<size_t>(iface)] = 0.5 * (
            left_cell + cell[static_cast<size_t>(iface)]);
    }
    close_right_face_alias(face, nxl, 1, mpi_rank, mpi_size, message_tag + 1);
}

} // namespace PeriodicStaggered
