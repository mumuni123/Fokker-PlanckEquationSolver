// Stage H2 acceptance: transactional safety (section 15 H2).  A successful
// conversion on trial states must leave the accepted bulk/tail untouched and
// advance only the trial tail ID counter; a failing conversion (unphysical
// conversion-region mass, invalid extraction request) must modify nothing,
// consume no particle IDs and leave the accepted state bitwise unchanged.
//
// Usage:
//   bulk_tail_transaction_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "grid.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

struct Metrics {
    bool success_conversion_ok;
    bool accepted_unchanged_after_success;
    bool trial_id_advanced;
    bool accepted_id_unchanged;
    bool failure_rejected;
    bool failure_bulk_untouched;
    bool failure_no_id_consumed;
    bool extraction_overbudget_rejected;
    bool extraction_oob_rejected;
    bool swap_roundtrip_ok;
    Metrics()
        : success_conversion_ok(false), accepted_unchanged_after_success(false),
          trial_id_advanced(false), accepted_id_unchanged(false),
          failure_rejected(false), failure_bulk_untouched(false),
          failure_no_id_consumed(false), extraction_overbudget_rejected(false),
          extraction_oob_rejected(false), swap_roundtrip_ok(false)
    {}
};

// Deterministic conversion packet in local cell 0.
void fill_packet(Species& bulk, const SpatialGrid& grid,
                 const HybridVelocityPartition& partition)
{
    const int ng = grid.nghost;
    std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            if (partition.is_conversion(j, k)) {
                bulk.f[idx3(ng, j, k)] =
                    1.0e18 * (1.0 + static_cast<double>((j + k) % 4));
            }
        }
    }
}

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 12, 1.2 * Const::micro);
    Species accepted_bulk;
    accepted_bulk.init("background_electrons",
                       SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                       Const::me, Param::dens, Param::temperature_e,
                       false, grid);
    HybridVelocityPartition partition;
    partition.init(accepted_bulk.cgrid, 6.0, 1.0, 4, 4);
    fill_packet(accepted_bulk, grid, partition);
    accepted_bulk.compute_moments();
    BackgroundTailPIC accepted_tail;
    accepted_tail.init(grid);
    const std::uint64_t accepted_id_before = accepted_tail.id_counter();
    const std::vector<double> accepted_f_before = accepted_bulk.f;

    // ---- success on trial states ----
    Species trial_bulk = accepted_bulk;
    BackgroundTailPIC trial_tail;
    trial_tail.init(grid);
    const std::uint64_t trial_id_before = trial_tail.id_counter();
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d =
        converter.extract_after_substep(trial_bulk, trial_tail, grid,
                                        partition, 5,
                                        ConversionLocation::AFTER_U_SUBSTEP,
                                        0);
    m.success_conversion_ok =
        d.complete && d.conservative && d.fidelity_ok && d.finite &&
        d.particles_created > 0;
    m.accepted_unchanged_after_success =
        accepted_bulk.f == accepted_f_before &&
        accepted_tail.particles.empty() &&
        accepted_tail.id_counter() == accepted_id_before;
    m.trial_id_advanced =
        trial_tail.id_counter() ==
        trial_id_before + d.particles_created;
    m.accepted_id_unchanged =
        accepted_tail.id_counter() == accepted_id_before;

    // ---- failure: unphysical conversion-region mass ----
    Species bad_bulk = accepted_bulk;
    BackgroundTailPIC bad_tail;
    bad_tail.init(grid);
    const std::uint64_t bad_id_before = bad_tail.id_counter();
    // Find a conversion cell and corrupt it with a negative mass.
    {
        const int ng = grid.nghost;
        bool corrupted = false;
        for (int j = 0; j < Param::Nv && !corrupted; ++j) {
            for (int k = 0; k < Param::Nmu && !corrupted; ++k) {
                if (partition.is_conversion(j, k)) {
                    bad_bulk.f[idx3(ng, j, k)] = -5.0e18;
                    corrupted = true;
                }
            }
        }
    }
    // Snapshot after the intentional corruption: "untouched" means the
    // converter must not modify the trial further.
    const std::vector<double> bad_f_before = bad_bulk.f;
    const BulkTailConversionDiagnostics fail =
        converter.extract_after_substep(bad_bulk, bad_tail, grid, partition,
                                        5,
                                        ConversionLocation::AFTER_U_SUBSTEP,
                                        0);
    m.failure_rejected =
        !fail.complete && !fail.conservative && !fail.finite;
    m.failure_bulk_untouched = bad_bulk.f == bad_f_before;
    m.failure_no_id_consumed =
        bad_tail.particles.empty() && bad_tail.id_counter() == bad_id_before;

    // ---- extraction API validation ----
    Species probe = accepted_bulk;
    const std::vector<double> probe_f_before = probe.f;
    {
        const int ng = grid.nghost;
        int j = -1;
        int k = -1;
        for (int jj = 0; jj < Param::Nv; ++jj) {
            for (int kk = 0; kk < Param::Nmu; ++kk) {
                if (partition.is_conversion(jj, kk)) {
                    j = jj;
                    k = kk;
                    break;
                }
            }
            if (j >= 0) break;
        }
        const double stored = probe.f[idx3(ng, j, k)];
        std::vector<ConversionMassRequest> overbudget;
        ConversionMassRequest req;
        req.ix_global = 0;
        req.iv = j;
        req.imu = k;
        req.mass = stored * 1.5;   // clearly over budget (beyond roundoff)
        overbudget.push_back(req);
        const SpeciesConversionResult r1 =
            probe.extract_conversion_masses(overbudget);
        m.extraction_overbudget_rejected = !r1.valid && !r1.complete;

        std::vector<ConversionMassRequest> oob;
        req.ix_global = grid.ix_start + grid.nx_local;  // outside local range
        req.mass = stored;
        oob.push_back(req);
        const SpeciesConversionResult r2 =
            probe.extract_conversion_masses(oob);
        m.extraction_oob_rejected = !r2.valid && !r2.complete;
    }
    m.failure_bulk_untouched =
        m.failure_bulk_untouched && (probe.f == probe_f_before);

    // ---- swap roundtrip ----
    Species swapped_bulk = accepted_bulk;
    Species holder_bulk = trial_bulk;
    const std::vector<double> holder_f = holder_bulk.f;
    swapped_bulk.swap_state(holder_bulk);
    const bool first_swap_ok = swapped_bulk.f == holder_f;
    swapped_bulk.swap_state(holder_bulk);
    const bool roundtrip_ok = swapped_bulk.f == accepted_f_before;
    m.swap_roundtrip_ok = first_swap_ok && roundtrip_ok;
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "success_conversion_ok=" << (m.success_conversion_ok ? 1 : 0)
        << "\n";
    out << "accepted_unchanged_after_success="
        << (m.accepted_unchanged_after_success ? 1 : 0) << "\n";
    out << "trial_id_advanced=" << (m.trial_id_advanced ? 1 : 0) << "\n";
    out << "accepted_id_unchanged=" << (m.accepted_id_unchanged ? 1 : 0)
        << "\n";
    out << "failure_rejected=" << (m.failure_rejected ? 1 : 0) << "\n";
    out << "failure_bulk_untouched=" << (m.failure_bulk_untouched ? 1 : 0)
        << "\n";
    out << "failure_no_id_consumed=" << (m.failure_no_id_consumed ? 1 : 0)
        << "\n";
    out << "extraction_overbudget_rejected="
        << (m.extraction_overbudget_rejected ? 1 : 0) << "\n";
    out << "extraction_oob_rejected=" << (m.extraction_oob_rejected ? 1 : 0)
        << "\n";
    out << "swap_roundtrip_ok=" << (m.swap_roundtrip_ok ? 1 : 0) << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "bulk_tail_transaction_test must run with exactly "
                     "1 rank; use plain ./build_hybrid/bulk_tail_"
                     "transaction_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: bulk_tail_transaction_test [--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case();
    bool pass = ok && m.success_conversion_ok &&
                m.accepted_unchanged_after_success &&
                m.trial_id_advanced && m.accepted_id_unchanged &&
                m.failure_rejected && m.failure_bulk_untouched &&
                m.failure_no_id_consumed &&
                m.extraction_overbudget_rejected &&
                m.extraction_oob_rejected && m.swap_roundtrip_ok;
    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "success_conversion_ok=" << (m.success_conversion_ok ? 1 : 0)
              << " accepted_unchanged_after_success="
              << (m.accepted_unchanged_after_success ? 1 : 0)
              << " trial_id_advanced=" << (m.trial_id_advanced ? 1 : 0)
              << " accepted_id_unchanged=" << (m.accepted_id_unchanged ? 1 : 0)
              << " failure_rejected=" << (m.failure_rejected ? 1 : 0)
              << " failure_bulk_untouched="
              << (m.failure_bulk_untouched ? 1 : 0)
              << " failure_no_id_consumed="
              << (m.failure_no_id_consumed ? 1 : 0)
              << " extraction_overbudget_rejected="
              << (m.extraction_overbudget_rejected ? 1 : 0)
              << " extraction_oob_rejected="
              << (m.extraction_oob_rejected ? 1 : 0)
              << " swap_roundtrip_ok=" << (m.swap_roundtrip_ok ? 1 : 0)
              << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
