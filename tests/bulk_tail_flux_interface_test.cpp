#include "grid.h"
#include "parameters.h"
#include "tail_subcell_quadrature.h"
#include "conservative_ppm_remap.h"
#include "species.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <set>
#include <string>
#include <vector>

namespace {
bool write_result(const std::string& path, bool ok,
                  const HybridVelocityPartition& p,
                  std::size_t hole_count, std::size_t axis_count,
                  std::size_t duplicate_count,
                  std::size_t negative_node_count,
                  std::size_t below_count, std::size_t supported_node_count,
                  double below_relative, bool null_interface_bitwise_equal)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str());
    if (!out) return false;
    out << "status=" << (ok ? "PASS" : "FAIL") << "\n"
        << "tail_cells=";
    std::size_t tail = 0;
    for (size_t i = 0; i < p.tail_owned_cell.size(); ++i) tail += p.tail_owned_cell[i] != 0;
    out << tail << "\n"
        << "upar_interface_faces=" << p.upar_interface_faces.size() << "\n"
        << "uperp_interface_faces=" << p.uperp_interface_faces.size() << "\n"
        << "interface_hole_count=" << hole_count << "\n"
        << "axis_interface_count=" << axis_count << "\n"
        << "interface_duplicate_count=" << duplicate_count << "\n"
        << "negative_node_count=" << negative_node_count << "\n"
        << "below_threshold_node_count=" << below_count << "\n"
        << "below_threshold_number_relative=" << below_relative << "\n"
        << "supported_node_count=" << supported_node_count << "\n"
        << "null_interface_bitwise_equal="
        << (null_interface_bitwise_equal ? 1 : 0) << "\n"
        << "topology_hash=" << p.topology_mask_hash() << "\n";
    return static_cast<bool>(out);
}
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string result_path;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--result" && i + 1 < argc)
            result_path = argv[++i];
    }
    bool ok = true;
    std::size_t hole_count = 0;
    std::size_t axis_count = 0;
    std::size_t below_count = 0;
    std::size_t duplicate_count = 0;
    std::size_t negative_node_count = 0;
    std::size_t supported_node_count = 0;
    double below_relative = 0.0;
    double total_node_weight = 0.0;
    bool null_interface_bitwise_equal = false;
    HybridVelocityPartition p;
    try {
        CylindricalVelocityGrid cg;
        cg.init(Param::momentum_umax, Param::Nv, Param::Nmu);
        p.init(cg, 6.0, 1.0, 4, 4);
        std::string reason;
        ok = p.flux_interface_topology_valid(&reason);
        for (int j = 0; j < p.upar_count; ++j) {
            for (int k = 0; k < p.uperp_count; ++k) {
                if (p.is_tail_owned(j, k) !=
                    p.is_tail_owned(p.upar_count - 1 - j, k)) ok = false;
            }
        }
        std::set<std::pair<int, int> > faces;
        for (size_t i = 0; i < p.upar_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = p.upar_interface_faces[i];
            if (!faces.insert(std::make_pair(f.face_index,
                                             f.transverse_index)).second) {
                ++duplicate_count;
                ok = false;
            }
            if (f.tail_iv < 0 || f.tail_iv >= p.upar_count ||
                f.tail_imu < 0 || f.tail_imu >= p.uperp_count) ok = false;
        }
        faces.clear();
        for (size_t i = 0; i < p.uperp_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = p.uperp_interface_faces[i];
            if (f.face_index <= 0 ||
                !faces.insert(std::make_pair(f.face_index,
                                              f.transverse_index)).second) {
                if (f.face_index > 0) ++duplicate_count;
                ok = false;
            }
            if (f.face_index <= 0) ++axis_count;
        }
        // Recheck the exterior-connected tail topology and all 4x4
        // threshold nodes independently of the production validator.
        const int nv = p.upar_count;
        const int nmu = p.uperp_count;
        std::vector<unsigned char> seen(
            static_cast<size_t>(nv) * static_cast<size_t>(nmu), 0);
        std::vector<std::pair<int, int> > stack;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (p.is_tail_owned(j, k) &&
                    (j == 0 || j == nv - 1 || k == nmu - 1))
                    stack.push_back(std::make_pair(j, k));
            }
        }
        const int dj[4] = {-1, 1, 0, 0};
        const int dk[4] = {0, 0, -1, 1};
        while (!stack.empty()) {
            const std::pair<int, int> cell = stack.back();
            stack.pop_back();
            const size_t slot = static_cast<size_t>(cell.first) * nmu +
                                static_cast<size_t>(cell.second);
            if (seen[slot]) continue;
            seen[slot] = 1;
            for (int d = 0; d < 4; ++d) {
                const int jn = cell.first + dj[d];
                const int kn = cell.second + dk[d];
                if (jn >= 0 && jn < nv && kn >= 0 && kn < nmu &&
                    p.is_tail_owned(jn, kn))
                    stack.push_back(std::make_pair(jn, kn));
            }
        }
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (p.is_tail_owned(j, k) &&
                    !seen[static_cast<size_t>(j) * nmu + k]) ++hole_count;
            }
        }
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (!p.is_tail_owned(j, k)) continue;
                const std::vector<TailSubcellNode> nodes =
                    TailSubcellQuadrature::nodes(cg, j, k, 4);
                for (size_t q = 0; q < nodes.size(); ++q) {
                    ++supported_node_count;
                    total_node_weight += nodes[q].mass_fraction;
                    if (!std::isfinite(nodes[q].upar) ||
                        !std::isfinite(nodes[q].uperp) ||
                        !std::isfinite(nodes[q].mass_fraction) ||
                        nodes[q].mass_fraction < 0.0) {
                        ++negative_node_count;
                        ok = false;
                    }
                    if (nodes[q].kinetic_energy < p.min_conversion_energy -
                        64.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, p.min_conversion_energy)) {
                        ++below_count;
                        below_relative += nodes[q].mass_fraction;
                    }
                }
            }
        }
        if (total_node_weight > 0.0)
            below_relative /= total_node_weight;
        // The NULL interface is a strict no-op contract.  Compare the
        // production remap with no interface object against the same remap
        // with the partition present but interface sink disabled.
        SpatialGrid test_grid;
        test_grid.init_with_domain(rank, 1, 4, 4.0 * Const::micro);
        Species test_input, test_null, test_interface;
        test_input.init("background", SpeciesType::BACKGROUND_ELECTRON,
                        -Const::qe, Const::me, Param::dens,
                        Param::temperature_e, false, test_grid);
        test_null = test_input;
        test_interface = test_input;
        for (size_t q = 0; q < test_input.f.size(); ++q)
            test_input.f[q] = 1.0e-12 * static_cast<double>(1 + q % 17);
        test_input.f_tmp = test_input.f;
        EMFields test_fields;
        test_fields.init(test_grid);
        ConservativePpmRemap test_remap;
        test_remap.init(test_grid, cg);
        const RemapDiagnostics dnull = test_remap.advect_u_parallel(
            test_input, test_null, test_fields, 1.0e-18, 0.0,
            NULL, NULL, 4);
        BulkTailFluxBatch no_sink;
        no_sink.apply_interface_sink = false;
        const RemapDiagnostics dinterface = test_remap.advect_u_parallel(
            test_input, test_interface, test_fields, 1.0e-18, 0.0,
            &p, &no_sink, 4);
        null_interface_bitwise_equal = dnull.finite && dinterface.finite &&
            test_null.f == test_interface.f;
        if (below_relative > 1.0e-12) ok = false;
        if (hole_count != 0 || axis_count != 0 || duplicate_count != 0 ||
            negative_node_count != 0 || !null_interface_bitwise_equal) {
            ok = false;
        }
    } catch (const std::exception& e) {
        ok = false;
        if (rank == 0) std::cerr << "bulk-tail interface test: " << e.what() << "\n";
    }
    int global_ok = ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        if (!write_result(result_path, global_ok != 0, p, hole_count,
                          axis_count, duplicate_count, negative_node_count,
                          below_count, supported_node_count, below_relative,
                          null_interface_bitwise_equal)) global_ok = 0;
        std::cout << "status=" << (global_ok ? "PASS" : "FAIL") << "\n";
    }
    MPI_Bcast(&global_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
