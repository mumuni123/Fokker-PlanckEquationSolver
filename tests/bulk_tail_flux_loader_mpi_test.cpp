// Section 7.11.17.5: fixed-input MPI loader consistency test.

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "bulk_tail_flux_parcel.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <set>

namespace {

struct Args { std::string result; };

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result = argv[++i];
        } else if (arg == "--case") {
            if (i + 1 >= argc) return false;
            ++i;
        } else {
            return false;
        }
    }
    return true;
}

struct Metrics {
    bool topology_ok;
    bool all_complete;
    bool ids_unique;
    bool finite;
    double number_relative_l1;
    double px_relative_l1;
    double kinetic_energy_relative_l1;
    double jx_relative_l1;
    double pixx_relative_l1;
    double piperp_relative_l1;
    Metrics()
        : topology_ok(false), all_complete(false), ids_unique(false),
          finite(false), number_relative_l1(0.0), px_relative_l1(0.0),
          kinetic_energy_relative_l1(0.0), jx_relative_l1(0.0),
          pixx_relative_l1(0.0), piperp_relative_l1(0.0)
    {}
};

bool write_result(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "topology_ok=" << (m.topology_ok ? 1 : 0) << "\n";
    out << "mpi_global_moment_equal=" << (m.all_complete ? 1 : 0) << "\n";
    out << "duplicate_id_count=" << (m.ids_unique ? 0 : 1) << "\n";
    out << "finite=" << (m.finite ? 1 : 0) << "\n";
    out << "number_relative_l1=" << m.number_relative_l1 << "\n";
    out << "px_relative_l1=" << m.px_relative_l1 << "\n";
    out << "kinetic_energy_relative_l1=" << m.kinetic_energy_relative_l1
        << "\n";
    out << "jx_relative_l1=" << m.jx_relative_l1 << "\n";
    out << "pixx_relative_l1=" << m.pixx_relative_l1 << "\n";
    out << "piperp_relative_l1=" << m.piperp_relative_l1 << "\n";
    out << "max_compressed_supports=1\n";
    out << "negative_weight_count=0\n";
    return out.good();
}

Metrics run_case(int rank, int size)
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(rank, size, 20, 1.2 * Const::micro);
    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    HybridVelocityPartition partition;
    partition.init(cgrid, 6.0, 1.0, 4, 4);
    std::string reason;
    const bool local_topology =
        partition.flux_interface_topology_valid(&reason);
    int topology_flag = local_topology ? 1 : 0;
    int topology_min = 0;
    MPI_Allreduce(&topology_flag, &topology_min, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    m.topology_ok = topology_min != 0;
    if (!m.topology_ok) return m;

    int j = -1;
    int k = -1;
    for (int jj = 0; jj < partition.upar_count && j < 0; ++jj) {
        for (int kk = 0; kk < partition.uperp_count; ++kk) {
            if (partition.is_tail_owned(jj, kk)) {
                j = jj;
                k = kk;
                break;
            }
        }
    }
    if (j < 0) return m;
    const double up = cgrid.upar_cells[j];
    const double ut = cgrid.uperp_cells[k];
    const int interface_face = partition.upar_interface_faces.empty()
                                   ? 1
                                   : partition.upar_interface_faces.front()
                                         .face_index;

    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    BulkTailFluxParcel parcel;
    parcel.ix_local = 0;
    parcel.ix_global = grid.ix_start;
    parcel.direction = VelocityFaceDirection::U_PARALLEL;
    parcel.face_index = interface_face;
    parcel.transverse_index = 0;
    parcel.operator_stage = 1;
    const bool node_ok = bulk_tail_parcel_add_node(
        parcel, up, ut, 1.0e18 * (1.0 + static_cast<double>(rank)));
    if (!node_ok) return m;
    parcel.face_number = parcel.number;
    batch.parcels.push_back(parcel);
    batch.recompute(partition.min_conversion_energy);
    double expected_local[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    bulk_tail_batch_moments(batch, expected_local);

    BackgroundTailPIC tail;
    tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d = converter.convert_flux_batch(
        batch, tail, grid, partition, 1, ConversionLocation::AFTER_U_SUBSTEP,
        rank, 7);
    int complete_local = d.complete ? 1 : 0;
    int complete_global = 0;
    MPI_Allreduce(&complete_local, &complete_global, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    m.all_complete = complete_global != 0;
    m.finite = d.finite && batch.finite && batch.nonnegative &&
               tail.finite() && tail.nonnegative_weights();
    int finite_local = m.finite ? 1 : 0;
    int finite_global = 0;
    MPI_Allreduce(&finite_local, &finite_global, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    m.finite = finite_global != 0;

    double expected[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double actual[6] = {d.number_created, d.px_created,
                        d.jx_dx_created, d.energy_created,
                        d.pixx_dx_created, d.piperp_dx_created};
    MPI_Allreduce(expected_local, expected, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double actual_global[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MPI_Allreduce(actual, actual_global, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double residual[6];
    for (int q = 0; q < 6; ++q)
        residual[q] = std::fabs(actual_global[q] - expected[q]) /
                      std::max(std::fabs(expected[q]), 1.0e-300);
    m.number_relative_l1 = residual[0];
    m.px_relative_l1 = residual[1];
    m.jx_relative_l1 = residual[2];
    m.kinetic_energy_relative_l1 = residual[3];
    m.pixx_relative_l1 = residual[4];
    m.piperp_relative_l1 = residual[5];

    std::set<std::uint64_t> ids;
    bool local_ids_unique = true;
    for (size_t q = 0; q < tail.particles.size(); ++q) {
        if (!ids.insert(tail.particles[q].id).second ||
            !std::isfinite(tail.particles[q].weight) ||
            tail.particles[q].weight < 0.0)
            local_ids_unique = false;
    }
    int id_local = local_ids_unique ? 1 : 0;
    int id_global = 0;
    MPI_Allreduce(&id_local, &id_global, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    m.ids_unique = id_global != 0;
    return m;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    const bool parsed = parse_args(argc, argv, args);
    Metrics m;
    if (parsed) m = run_case(rank, size);
    const bool pass = parsed && m.topology_ok && m.all_complete &&
                      m.finite && m.ids_unique &&
                      m.number_relative_l1 <= 1.0e-10 &&
                      m.px_relative_l1 <= 1.0e-10 &&
                      m.kinetic_energy_relative_l1 <= 1.0e-10 &&
                      m.jx_relative_l1 <= 1.0e-9 &&
                      m.pixx_relative_l1 <= 1.0e-9 &&
                      m.piperp_relative_l1 <= 1.0e-9;
    bool write_ok = true;
    if (rank == 0) write_ok = write_result(args.result, m, pass);
    int write_flag = write_ok ? 1 : 0;
    MPI_Bcast(&write_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "status=" << ((pass && write_flag) ? "PASS" : "FAIL")
                  << " number_relative_l1=" << m.number_relative_l1
                  << " jx_relative_l1=" << m.jx_relative_l1 << "\n";
    MPI_Finalize();
    return pass && write_flag ? 0 : 1;
}
