#include "grid.h"
#include "maxwell.h"
#include "parameters.h"
#include "periodic_staggered_operators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

namespace {

double face_value(int global_face, int component)
{
    return 0.125 * (global_face + 1) * (component + 1);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SpatialGrid sg;
    sg.init(rank, size);
    const int nxl = sg.nx_local;
    const int block = 3;
    std::vector<double> fx(static_cast<size_t>(nxl + 1) * block, 0.0);
    std::vector<double> jn(static_cast<size_t>(nxl + 1), 0.0);
    for (int iface = 0; iface < nxl; ++iface) {
        const int global_face = sg.ix_start + iface;
        for (int q = 0; q < block; ++q)
            fx[static_cast<size_t>(iface) * block + q] = face_value(global_face, q);
        jn[static_cast<size_t>(iface)] =
            fx[static_cast<size_t>(iface) * block] -
            fx[static_cast<size_t>(iface) * block + 1];
    }
    PeriodicStaggered::close_right_face_alias(fx, nxl, block, rank, size, 4300);
    PeriodicStaggered::close_right_face_alias(jn, nxl, 1, rank, size, 4310);

    const int alias_global_face = (sg.ix_start + nxl) % sg.nx_global;
    double local_alias_error = 0.0;
    for (int q = 0; q < block; ++q)
        local_alias_error = std::max(local_alias_error,
            std::fabs(fx[static_cast<size_t>(nxl) * block + q] -
                      face_value(alias_global_face, q)));
    local_alias_error = std::max(local_alias_error,
        std::fabs(jn[static_cast<size_t>(nxl)] -
                  (fx[static_cast<size_t>(nxl) * block] -
                   fx[static_cast<size_t>(nxl) * block + 1])));

    double local_jn_sum = 0.0;
    for (int iface = 0; iface < nxl; ++iface)
        local_jn_sum += jn[static_cast<size_t>(iface)];
    double global_jn_sum = 0.0;
    int global_owned_faces = 0;
    MPI_Allreduce(&local_jn_sum, &global_jn_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&nxl, &global_owned_faces, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    double expected_jn_sum = 0.0;
    for (int iface = 0; iface < sg.nx_global; ++iface)
        expected_jn_sum += face_value(iface, 0) - face_value(iface, 1);

    // Production output gathers owned faces only.  Reproduce that collection
    // here and verify that the periodic alias never becomes an extra row.
    std::vector<int> counts(static_cast<size_t>(size), 0);
    std::vector<int> displacements(static_cast<size_t>(size), 0);
    MPI_Allgather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int irank = 1; irank < size; ++irank)
        displacements[static_cast<size_t>(irank)] =
            displacements[static_cast<size_t>(irank - 1)] +
            counts[static_cast<size_t>(irank - 1)];
    std::vector<double> gathered_owned_faces(rank == 0
        ? static_cast<size_t>(global_owned_faces) : 0, 0.0);
    MPI_Gatherv(jn.data(), nxl, MPI_DOUBLE,
                rank == 0 ? gathered_owned_faces.data() : 0,
                counts.data(), displacements.data(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    double gathered_output_error = 0.0;
    if (rank == 0) {
        for (int iface = 0; iface < global_owned_faces; ++iface) {
            const double expected = face_value(iface, 0) - face_value(iface, 1);
            gathered_output_error = std::max(gathered_output_error,
                std::fabs(gathered_owned_faces[static_cast<size_t>(iface)] - expected));
        }
    }

    EMFields fields;
    fields.init(sg);
    const double dt = 1.0e-18;
    fields.advance_ampere_face_from_midpoint_current(jn, dt, rank, size);
    const double ampere_scale = -dt / Const::eps0;
    double local_ampere_error = 0.0;
    for (int iface = 0; iface < nxl; ++iface)
        local_ampere_error = std::max(local_ampere_error,
            std::fabs(fields.Ex_face[static_cast<size_t>(iface)] -
                      ampere_scale * jn[static_cast<size_t>(iface)]));
    double local_field_alias_error = std::fabs(
        fields.Ex_face[static_cast<size_t>(nxl)] - ampere_scale * jn[static_cast<size_t>(nxl)]);

    double global_alias_error = 0.0;
    double global_ampere_error = 0.0;
    double global_field_alias_error = 0.0;
    double global_gathered_output_error = 0.0;
    MPI_Allreduce(&local_alias_error, &global_alias_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_ampere_error, &global_ampere_error, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_field_alias_error, &global_field_alias_error, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&gathered_output_error, &global_gathered_output_error, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const bool passes = global_owned_faces == sg.nx_global &&
        std::fabs(global_jn_sum - expected_jn_sum) <= 1.0e-12 *
            std::max(1.0, std::fabs(expected_jn_sum)) &&
        global_alias_error <= 1.0e-14 && global_ampere_error <= 1.0e-14 &&
        global_field_alias_error <= 1.0e-14 &&
        global_gathered_output_error <= 1.0e-14;
    if (rank == 0) {
        std::ofstream out("periodic_x_flux_topology_test.result");
        out << std::setprecision(17)
            << "mpi_ranks " << size << "\n"
            << "unique_owned_face_count " << global_owned_faces << "\n"
            << "expected_unique_face_count " << sg.nx_global << "\n"
            << "alias_error " << global_alias_error << "\n"
            << "unique_jn_integral_error " << global_jn_sum - expected_jn_sum << "\n"
            << "ampere_owned_face_error " << global_ampere_error << "\n"
            << "ampere_alias_error " << global_field_alias_error << "\n"
            << "output_owned_face_error " << global_gathered_output_error << "\n"
            << "output_alias_counted_once "
            << (global_owned_faces == sg.nx_global ? 1 : 0) << "\n"
            << "passes " << (passes ? 1 : 0) << "\n";
        std::cout << "periodic_x_flux_topology_test: passes=" << passes
                  << " owned_faces=" << global_owned_faces << std::endl;
    }
    MPI_Finalize();
    return passes ? 0 : 1;
}
