// Section 7.11.17.5: failed flux-parcel conversion must be transactional.

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "bulk_tail_flux_parcel.h"

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

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
    bool failure_rejected;
    bool state_unchanged;
    bool id_unchanged;
    bool ledgers_unchanged;
    bool finite_failure_classified;
    bool negative_failure_rejected;
    Metrics()
        : topology_ok(false), failure_rejected(false),
          state_unchanged(false), id_unchanged(false),
          ledgers_unchanged(false), finite_failure_classified(false),
          negative_failure_rejected(false)
    {}
};

bool snapshots_equal(const BackgroundTailStateSnapshot& a,
                     const BackgroundTailStateSnapshot& b)
{
    return a.particles.size() == b.particles.size() &&
           a.density == b.density && a.id_counter == b.id_counter &&
           a.outflow.left_number == b.outflow.left_number &&
           a.outflow.left_px == b.outflow.left_px &&
           a.outflow.left_kinetic_energy == b.outflow.left_kinetic_energy &&
           a.outflow.right_number == b.outflow.right_number &&
           a.outflow.right_px == b.outflow.right_px &&
           a.outflow.right_kinetic_energy == b.outflow.right_kinetic_energy &&
           a.truncation_shape_left == b.truncation_shape_left &&
           a.truncation_shape_right == b.truncation_shape_right &&
           a.deposit_shape_left == b.deposit_shape_left &&
           a.deposit_shape_right == b.deposit_shape_right &&
           a.collision_rng_seed == b.collision_rng_seed;
}

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 12, 1.2 * Const::micro);
    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    HybridVelocityPartition partition;
    partition.init(cgrid, 6.0, 1.0, 4, 4);
    std::string reason;
    m.topology_ok = partition.flux_interface_topology_valid(&reason);
    if (!m.topology_ok) return m;

    int selected_j = -1;
    int selected_k = -1;
    for (int j = 0; j < partition.upar_count && selected_j < 0; ++j) {
        for (int k = 0; k < partition.uperp_count; ++k) {
            if (partition.is_tail_owned(j, k)) {
                selected_j = j;
                selected_k = k;
                break;
            }
        }
    }
    if (selected_j < 0) return m;

    BackgroundTailPIC tail;
    tail.init(grid);
    BackgroundTailParticle p;
    p.x = 0.5 * grid.dx;
    p.ux = cgrid.upar_cells[selected_j];
    p.uy = cgrid.uperp_cells[selected_k];
    p.uz = 0.0;
    p.weight = 3.0e18;
    p.id = tail.next_particle_id(0);
    tail.particles.push_back(p);
    BackgroundTailStateSnapshot before;
    tail.export_accepted_state(before);

    BulkTailFluxBatch bad;
    bad.apply_interface_sink = true;
    BulkTailFluxParcel parcel;
    parcel.ix_local = 0;
    parcel.ix_global = 0;
    parcel.direction = VelocityFaceDirection::U_PARALLEL;
    parcel.face_index = partition.upar_interface_faces.empty()
                            ? 1 : partition.upar_interface_faces.front().face_index;
    parcel.transverse_index = 0;
    parcel.operator_stage = 1;
    FluxParcelNode invalid;
    invalid.upar = p.ux;
    invalid.uperp = p.uy;
    invalid.mass = std::numeric_limits<double>::quiet_NaN();
    parcel.nodes.push_back(invalid);
    bad.parcels.push_back(parcel);
    bad.recompute(partition.min_conversion_energy);

    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d = converter.convert_flux_batch(
        bad, tail, grid, partition, 2, ConversionLocation::AFTER_U_SUBSTEP,
        0, 7);
    BackgroundTailStateSnapshot after;
    tail.export_accepted_state(after);
    m.failure_rejected = !d.complete;
    m.finite_failure_classified = !d.finite && !d.conservative &&
                                  !d.fidelity_ok;
    m.state_unchanged = snapshots_equal(before, after) &&
                        before.particles.size() == after.particles.size();
    m.id_unchanged = before.id_counter == after.id_counter;
    m.ledgers_unchanged = before.outflow.left_number ==
                              after.outflow.left_number &&
                          before.outflow.right_number ==
                              after.outflow.right_number &&
                          before.truncation_shape_left ==
                              after.truncation_shape_left &&
                          before.deposit_shape_right ==
                              after.deposit_shape_right;

    // A negative parcel weight is a distinct hard failure from a non-finite
    // node and must also leave the same accepted state untouched.
    BulkTailFluxBatch negative;
    negative.apply_interface_sink = true;
    BulkTailFluxParcel negative_parcel = parcel;
    negative_parcel.nodes.clear();
    FluxParcelNode negative_node;
    negative_node.upar = p.ux;
    negative_node.uperp = p.uy;
    negative_node.mass = -1.0;
    negative_parcel.nodes.push_back(negative_node);
    negative.parcels.push_back(negative_parcel);
    negative.recompute(partition.min_conversion_energy);
    const BulkTailConversionDiagnostics negative_result =
        converter.convert_flux_batch(negative, tail, grid, partition, 3,
                                     ConversionLocation::AFTER_U_SUBSTEP, 0, 7);
    BackgroundTailStateSnapshot after_negative;
    tail.export_accepted_state(after_negative);
    m.negative_failure_rejected =
        !negative_result.complete && !negative_result.finite &&
        snapshots_equal(before, after_negative) &&
        before.id_counter == after_negative.id_counter;
    return m;
}

bool write_result(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "topology_ok=" << (m.topology_ok ? 1 : 0) << "\n";
    out << "failed_trial_rejected=" << (m.failure_rejected ? 1 : 0) << "\n";
    out << "failed_trial_state_unchanged=" << (m.state_unchanged ? 1 : 0)
        << "\n";
    out << "failed_trial_next_id_unchanged=" << (m.id_unchanged ? 1 : 0)
        << "\n";
    out << "failed_trial_ledgers_unchanged="
        << (m.ledgers_unchanged ? 1 : 0) << "\n";
    out << "finite_failure_classified="
        << (m.finite_failure_classified ? 1 : 0) << "\n";
    out << "negative_weight_failure_rejected="
        << (m.negative_failure_rejected ? 1 : 0) << "\n";
    return out.good();
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
    if (parsed && size == 1) m = run_case();
    const bool pass = parsed && size == 1 && m.topology_ok &&
                      m.failure_rejected && m.state_unchanged &&
                      m.id_unchanged && m.ledgers_unchanged &&
                      m.finite_failure_classified &&
                      m.negative_failure_rejected;
    bool write_ok = write_result(args.result, m, pass);
    if (rank == 0)
        std::cout << "status=" << ((pass && write_ok) ? "PASS" : "FAIL")
                  << " failed_trial_state_unchanged="
                  << (m.state_unchanged ? 1 : 0) << "\n";
    MPI_Finalize();
    return pass && write_ok ? 0 : 1;
}
