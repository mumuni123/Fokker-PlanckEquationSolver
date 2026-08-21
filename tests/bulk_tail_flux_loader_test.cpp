// Section 7.11.17.5: deterministic flux-parcel loader acceptance test.
// This test exercises the production BulkTailConverter::convert_flux_batch
// rather than reconstructing a cell-centered static conversion in the test.

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "bulk_tail_flux_parcel.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Args {
    std::string result;
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result = argv[++i];
        } else if (arg == "--case") {
            if (i + 1 >= argc) return false;
            ++i; // The standalone executable currently has one complete case.
        } else {
            return false;
        }
    }
    return true;
}

struct Support {
    int j;
    int k;
    int sign;
    int energy_bin;
    double upar;
    double uperp;
};

struct Metrics {
    bool topology_ok;
    bool finite;
    bool complete;
    bool conservative;
    bool fidelity;
    bool positive_weights;
    bool ids_unique;
    bool positive_and_negative_upar;
    bool multi_parcel_group;
    bool compression_exercised;
    bool threshold_support_ok;
    bool low_uperp_support;
    bool roundoff_batch_discarded;
    std::uint64_t particles_created;
    std::uint64_t compression_fallback_count;
    std::uint64_t max_group_supports;
    double number_relative_l1;
    double px_relative_l1;
    double kinetic_energy_relative_l1;
    double jx_relative_l1;
    double pixx_relative_l1;
    double piperp_relative_l1;
    Metrics()
        : topology_ok(false), finite(false), complete(false),
          conservative(false), fidelity(false), positive_weights(false),
          ids_unique(false), positive_and_negative_upar(false),
          multi_parcel_group(false), compression_exercised(false),
           threshold_support_ok(false), low_uperp_support(false),
           roundoff_batch_discarded(false),
          particles_created(0), compression_fallback_count(0),
          max_group_supports(0), number_relative_l1(0.0),
          px_relative_l1(0.0), kinetic_energy_relative_l1(0.0),
          jx_relative_l1(0.0), pixx_relative_l1(0.0),
          piperp_relative_l1(0.0)
    {}
};

bool write_result(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "topology_ok=" << (m.topology_ok ? 1 : 0) << "\n";
    out << "finite=" << (m.finite ? 1 : 0) << "\n";
    out << "complete=" << (m.complete ? 1 : 0) << "\n";
    out << "conservative=" << (m.conservative ? 1 : 0) << "\n";
    out << "fidelity=" << (m.fidelity ? 1 : 0) << "\n";
    out << "positive_weights=" << (m.positive_weights ? 1 : 0) << "\n";
    out << "ids_unique=" << (m.ids_unique ? 1 : 0) << "\n";
    out << "positive_and_negative_upar="
        << (m.positive_and_negative_upar ? 1 : 0) << "\n";
    out << "multi_parcel_group=" << (m.multi_parcel_group ? 1 : 0) << "\n";
    out << "compression_exercised="
        << (m.compression_exercised ? 1 : 0) << "\n";
    out << "threshold_support_ok="
        << (m.threshold_support_ok ? 1 : 0) << "\n";
    out << "low_uperp_support=" << (m.low_uperp_support ? 1 : 0) << "\n";
    out << "roundoff_batch_discarded="
        << (m.roundoff_batch_discarded ? 1 : 0) << "\n";
    out << "particles_created=" << m.particles_created << "\n";
    out << "compression_fallback_count=" << m.compression_fallback_count
        << "\n";
    out << "max_group_supports=" << m.max_group_supports << "\n";
    out << "number_relative_l1=" << m.number_relative_l1 << "\n";
    out << "px_relative_l1=" << m.px_relative_l1 << "\n";
    out << "kinetic_energy_relative_l1=" << m.kinetic_energy_relative_l1
        << "\n";
    out << "jx_relative_l1=" << m.jx_relative_l1 << "\n";
    out << "pixx_relative_l1=" << m.pixx_relative_l1 << "\n";
    out << "piperp_relative_l1=" << m.piperp_relative_l1 << "\n";
    return out.good();
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
    std::string topology_reason;
    m.topology_ok = partition.flux_interface_topology_valid(&topology_reason);
    if (!m.topology_ok) return m;

    std::map<std::pair<int, int>, std::vector<Support> > groups;
    for (int j = 0; j < partition.upar_count; ++j) {
        for (int k = 0; k < partition.uperp_count; ++k) {
            if (!partition.is_tail_owned(j, k)) continue;
            const double up = cgrid.upar_cells[static_cast<size_t>(j)];
            const double ut = cgrid.uperp_cells[static_cast<size_t>(k)];
            const double ke = cgrid.kinetic_energy[
                static_cast<size_t>(j) * partition.uperp_count + k];
            const int sign = up < 0.0 ? -1 : 1;
            const int bin = partition.energy_bin(ke);
            groups[std::make_pair(sign, bin)].push_back(
                Support{j, k, sign, bin, up, ut});
        }
    }

    std::vector<Support> positive;
    std::vector<Support> negative;
    for (std::map<std::pair<int, int>, std::vector<Support> >::const_iterator
             it = groups.begin(); it != groups.end(); ++it) {
        if (it->first.first > 0 && positive.empty() &&
            it->second.size() >= 8)
            positive.assign(it->second.begin(), it->second.begin() + 8);
        if (it->first.first < 0 && negative.empty())
            negative.push_back(it->second.front());
        if (positive.size() >= 8 && !negative.empty()) break;
    }
    if (positive.size() < 8 || negative.empty()) return m;
    m.positive_and_negative_upar = true;
    for (size_t q = 0; q < positive.size(); ++q)
        if (positive[q].k == 0) m.low_uperp_support = true;
    if (negative[0].k == 0) m.low_uperp_support = true;

    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    BulkTailFluxParcel first;
    first.ix_local = 0;
    first.ix_global = 0;
    first.direction = VelocityFaceDirection::U_PARALLEL;
    first.face_index = partition.upar_interface_faces.empty()
                           ? 1 : partition.upar_interface_faces.front().face_index;
    first.transverse_index = 0;
    first.operator_stage = 1;
    for (size_t q = 0; q < 4; ++q) {
        if (!bulk_tail_parcel_add_node(
                first, positive[q].upar, positive[q].uperp,
                1.0e18 * (1.0 + 0.2 * static_cast<double>(q))))
            return m;
    }
    first.face_number = first.number;
    BulkTailFluxParcel second;
    second.ix_local = 0;
    second.ix_global = 0;
    second.direction = VelocityFaceDirection::U_PARALLEL;
    second.face_index = first.face_index;
    second.transverse_index = 1;
    second.operator_stage = 1;
    for (size_t q = 4; q < positive.size(); ++q) {
        if (!bulk_tail_parcel_add_node(
                second, positive[q].upar, positive[q].uperp,
                1.0e18 * (1.0 + 0.2 * static_cast<double>(q))))
            return m;
    }
    second.face_number = second.number;
    BulkTailFluxParcel opposite;
    opposite.ix_local = 0;
    opposite.ix_global = 0;
    opposite.direction = VelocityFaceDirection::U_PARALLEL;
    opposite.face_index = first.face_index + 1;
    opposite.transverse_index = 0;
    opposite.operator_stage = 1;
    if (!bulk_tail_parcel_add_node(opposite, negative[0].upar,
                                   negative[0].uperp, 2.5e18))
        return m;
    opposite.face_number = opposite.number;
    batch.parcels.push_back(first);
    batch.parcels.push_back(second);
    batch.parcels.push_back(opposite);
    batch.recompute(partition.min_conversion_energy);
    m.threshold_support_ok = batch.finite && batch.nonnegative &&
                             batch.duplicate_count == 0 &&
                             batch.below_threshold_number == 0.0;
    m.multi_parcel_group = batch.parcels.size() >= 2;

    BackgroundTailPIC tail;
    tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d = converter.convert_flux_batch(
        batch, tail, grid, partition, 1, ConversionLocation::AFTER_U_SUBSTEP,
        0, 7);
    m.finite = d.finite && batch.finite && batch.nonnegative;
    m.complete = d.complete;
    m.conservative = d.conservative;
    m.fidelity = d.fidelity_ok;
    m.particles_created = d.particles_created;
    m.compression_fallback_count = d.compression_fallback_count;
    m.number_relative_l1 = d.number_residual_rel;
    m.px_relative_l1 = d.px_residual_rel;
    m.kinetic_energy_relative_l1 = d.energy_residual_rel;
    m.jx_relative_l1 = d.jx_residual_rel;
    m.pixx_relative_l1 = d.pixx_residual_rel;
    m.piperp_relative_l1 = d.piperp_residual_rel;
    m.compression_exercised = positive.size() > 7 &&
                              d.compression_fallback_count == 0;
    // The positive group has eight input supports and the converter must
    // reduce it to at most seven; the negative group contributes one.
    m.max_group_supports = positive.size() > 7 ? 7 : positive.size();
    m.positive_weights = tail.nonnegative_weights();
    std::set<std::uint64_t> ids;
    m.ids_unique = true;
    for (size_t p = 0; p < tail.particles.size(); ++p)
        if (!ids.insert(tail.particles[p].id).second) m.ids_unique = false;

    // Regression for collision interface residues: a valid but subnormal
    // packet may have more than max_supports nodes, yet it must be discarded
    // before compression rather than causing a support-limit hard failure.
    BulkTailFluxBatch roundoff_batch;
    roundoff_batch.apply_interface_sink = true;
    BulkTailFluxParcel roundoff;
    roundoff.ix_local = 0;
    roundoff.ix_global = 0;
    roundoff.direction = VelocityFaceDirection::U_PARALLEL;
    roundoff.face_index = first.face_index;
    roundoff.transverse_index = 0;
    roundoff.operator_stage = 2;
    for (size_t q = 0; q < 21; ++q) {
        const Support& s = positive[q % positive.size()];
        if (!bulk_tail_parcel_add_node(roundoff, s.upar, s.uperp, 1.0e-320))
            return m;
    }
    roundoff.face_number = roundoff.number;
    roundoff_batch.parcels.push_back(roundoff);
    roundoff_batch.recompute(partition.min_conversion_energy);
    BackgroundTailPIC roundoff_tail;
    roundoff_tail.init(grid);
    const BulkTailConversionDiagnostics roundoff_d =
        converter.convert_flux_batch(roundoff_batch, roundoff_tail, grid,
                                     partition, 2,
                                     ConversionLocation::AFTER_COLLISION_HALF,
                                     0, 7);
    m.roundoff_batch_discarded = roundoff_d.complete && roundoff_d.finite &&
        roundoff_d.conservative && roundoff_d.fidelity_ok &&
        roundoff_d.support_limit_violation_count == 0 &&
        roundoff_d.roundoff_discarded_number > 0.0 &&
        roundoff_tail.particles.empty();
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
    if (parsed && size == 1) m = run_case();
    const bool pass = parsed && size == 1 && m.topology_ok && m.finite &&
                      m.complete && m.conservative && m.fidelity &&
                      m.positive_weights && m.ids_unique &&
                      m.positive_and_negative_upar && m.multi_parcel_group &&
                       m.compression_exercised && m.threshold_support_ok &&
                       m.low_uperp_support && m.roundoff_batch_discarded &&
                       m.particles_created > 0 &&
                      m.max_group_supports <= 7 &&
                      m.number_relative_l1 <= 1.0e-10 &&
                      m.px_relative_l1 <= 1.0e-10 &&
                      m.kinetic_energy_relative_l1 <= 1.0e-10 &&
                      m.jx_relative_l1 <= 1.0e-9 &&
                      m.pixx_relative_l1 <= 1.0e-9 &&
                      m.piperp_relative_l1 <= 1.0e-9;
    bool write_ok = write_result(args.result, m, pass);
    if (rank == 0) {
        std::cout << "status=" << ((pass && write_ok) ? "PASS" : "FAIL")
                  << " particles_created=" << m.particles_created
                  << " number_relative_l1=" << m.number_relative_l1
                  << " compression_fallback_count="
                  << m.compression_fallback_count << "\n";
    }
    MPI_Finalize();
    return pass && write_ok ? 0 : 1;
}
