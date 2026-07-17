#include "periodic_staggered_operators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

namespace {

struct Partition {
    int count;
    int start;
};

Partition partition(int global_count, int rank, int size)
{
    const int base = global_count / size;
    const int remainder = global_count % size;
    Partition p;
    p.count = base + (rank < remainder ? 1 : 0);
    p.start = rank * base + (rank < remainder ? rank : remainder);
    return p;
}

double reproducible_sample(int index, unsigned salt)
{
    unsigned value = static_cast<unsigned>(index + 1) * 1103515245u + salt;
    value = value * 1664525u + 1013904223u;
    return static_cast<double>(value % 2000001u) / 1000000.0 - 1.0;
}

bool check_case(int global_cells, int rank, int size, int basis,
                double& worst_relative, int& completed_cases)
{
    const Partition p = partition(global_cells, rank, size);
    std::vector<double> face(static_cast<size_t>(p.count + 1), 0.0);
    std::vector<double> je(static_cast<size_t>(p.count), 0.0);
    for (int iface = 0; iface < p.count; ++iface) {
        const int global_face = p.start + iface;
        face[static_cast<size_t>(iface)] =
            basis == 1 ? (global_face == 0 ? 1.0 : 0.0) :
            reproducible_sample(global_face, 17u) +
            (global_face == 0 ? 2.0 : 0.0);
    }
    for (int ix = 0; ix < p.count; ++ix) {
        const int global_cell = p.start + ix;
        je[static_cast<size_t>(ix)] =
            basis == 2 ? (global_cell == global_cells - 1 ? 1.0 : 0.0) :
            reproducible_sample(global_cell, 97u);
    }
    PeriodicStaggered::close_right_face_alias(face, p.count, 1, rank, size,
                                              4100 + basis);

    std::vector<double> ge;
    std::vector<double> gstar_je;
    PeriodicStaggered::apply_face_to_cell_G(face, ge, p.count);
    PeriodicStaggered::apply_cell_to_face_Gstar(je, gstar_je, p.count,
                                                rank, size, 4200 + basis);
    double local_lhs = 0.0;
    double local_rhs = 0.0;
    for (int ix = 0; ix < p.count; ++ix)
        local_lhs += ge[static_cast<size_t>(ix)] * je[static_cast<size_t>(ix)];
    for (int iface = 0; iface < p.count; ++iface)
        local_rhs += face[static_cast<size_t>(iface)] *
                     gstar_je[static_cast<size_t>(iface)];

    double lhs = 0.0;
    double rhs = 0.0;
    int global_unique_face_count = 0;
    MPI_Allreduce(&local_lhs, &lhs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_rhs, &rhs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const int local_unique_face_count = p.count;
    MPI_Allreduce(&local_unique_face_count, &global_unique_face_count, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    const double scale = std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
    const double relative = std::fabs(lhs - rhs) / scale;
    worst_relative = std::max(worst_relative, relative);
    ++completed_cases;
    return global_unique_face_count == global_cells && relative <= 4096.0e-15;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Each rank owns at least two cells.  Both parity cases exercise the
    // periodic seam across the last-rank/rank-zero interface.
    const int odd_cells = 2 * size + 1;
    const int even_cells = 2 * size + 2;
    double local_worst_relative = 0.0;
    int local_cases = 0;
    int local_ok = 1;
    local_ok &= check_case(odd_cells, rank, size, 0, local_worst_relative, local_cases);
    local_ok &= check_case(even_cells, rank, size, 0, local_worst_relative, local_cases);
    local_ok &= check_case(odd_cells, rank, size, 1, local_worst_relative, local_cases);
    local_ok &= check_case(even_cells, rank, size, 2, local_worst_relative, local_cases);

    int global_ok = 0;
    double global_worst_relative = 0.0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    MPI_Allreduce(&local_worst_relative, &global_worst_relative, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream out("periodic_staggered_adjoint_test.result");
        out << std::setprecision(17)
            << "mpi_ranks " << size << "\n"
            << "cases 4\n"
            << "unique_faces_checked 1\n"
            << "seam_basis_checked 1\n"
            << "last_cell_basis_checked 1\n"
            << "worst_relative_residual " << global_worst_relative << "\n"
            << "passes " << global_ok << "\n";
        std::cout << "periodic_staggered_adjoint_test: passes=" << global_ok
                  << " worst_relative=" << global_worst_relative << std::endl;
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
