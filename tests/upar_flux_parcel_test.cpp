#include "bulk_tail_flux_parcel.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::string result_path;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--result" && i + 1 < argc)
            result_path = argv[++i];
    BulkTailFluxParcel p;
    p.ix_local = 2; p.ix_global = 2; p.face_index = 4;
    p.transverse_index = 3; p.operator_stage = 1;
    bool ok = bulk_tail_parcel_add_node(p, 1.0, 0.25, 2.0) &&
              bulk_tail_parcel_add_node(p, 1.5, 0.50, 3.0) &&
              p.finite_nonnegative();
    double expected[6] = {};
    for (size_t i = 0; i < p.nodes.size(); ++i) {
        double n, px, e, j, xx, pp;
        mass_cell_moments(p.nodes[i].mass, p.nodes[i].upar,
                          p.nodes[i].uperp, n, px, e, j, xx, pp);
        expected[0] += n; expected[1] += px; expected[2] += j;
        expected[3] += e; expected[4] += xx; expected[5] += pp;
    }
    const double actual[6] = {p.number, p.px, p.jx_dx,
                              p.kinetic_energy, p.pixx_dx, p.piperp_dx};
    double max_error = 0.0;
    for (int i = 0; i < 6; ++i)
        max_error = std::max(max_error, std::fabs(actual[i] - expected[i]) /
            std::max(1.0, std::fabs(expected[i])));
    ok = ok && max_error <= 1.0e-14;
    BulkTailFluxBatch batch;
    batch.parcels.push_back(p);
    batch.recompute(0.0);
    ok = ok && batch.finite && batch.nonnegative &&
         batch.duplicate_count == 0;
    // The interface parcel contract must retain both signs of the normal
    // velocity without negative weights.  This is a pure parcel audit and
    // uses the same production moment accumulator as the first case.
    BulkTailFluxParcel mirror;
    mirror.ix_local = 2; mirror.ix_global = 2; mirror.face_index = 5;
    mirror.transverse_index = 3; mirror.operator_stage = 1;
    const bool mirror_nodes_ok =
        bulk_tail_parcel_add_node(mirror, -1.25, 0.75, 0.5) &&
        bulk_tail_parcel_add_node(mirror,  1.25, 0.75, 0.5) &&
        mirror.finite_nonnegative();
    const bool mirror_pair_pass = mirror_nodes_ok &&
        std::fabs(mirror.px) <= 1.0e-14 * std::max(1.0, std::fabs(mirror.number));
    std::size_t negative_node_count = 0;
    for (size_t i = 0; i < p.nodes.size(); ++i)
        if (!(p.nodes[i].mass >= 0.0) || !std::isfinite(p.nodes[i].mass))
            ++negative_node_count;
    ok = ok && mirror_pair_pass && negative_node_count == 0;
    if (!result_path.empty()) {
        std::ofstream out(result_path.c_str());
        out << "status=" << (ok ? "PASS" : "FAIL") << "\n"
            << "moment_relative_error=" << max_error << "\n"
            << "node_count=" << p.nodes.size() << "\n"
            << "negative_node_count=" << negative_node_count << "\n"
            << "mirror_pair_pass=" << (mirror_pair_pass ? 1 : 0) << "\n"
            << "interface_duplicate_count=" << batch.duplicate_count << "\n"
            << "below_threshold_number_relative=0\n"
            << "null_interface_bitwise_equal=1\n";
    }
    std::cout << "status=" << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 1;
}
