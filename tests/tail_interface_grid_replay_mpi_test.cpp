// Section 7.11.16C / 16B: MPI rank-consistency test for the fixed-(192,64)
// offline replay.  The real velocity histogram rows are partitioned across
// ranks deterministically (row index % size), aggregated per rank, reduced
// with the production aggregate rule (per-(iv,imu) SUM), and every
// candidate's faces and summary metrics must be identical across rank
// partitions (bit-identical faces, reduced masses consistent with the
// single-rank reference within 1e-10).  Rank 0 writes the result file.
//
// Usage (same CLI as the serial replay test):
//   tail_interface_grid_replay_mpi_test
//     --input <hist> --profiles G0,Gx,Gp,G2
//     --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0
//     --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 [--result path]

#include "grid.h"
#include "tail_interface_grid_design.h"
#include "tail_interface_replay_common.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <mpi.h>

namespace {

struct Args {
    std::string input;
    std::string profiles;
    std::vector<double> ax_values;
    std::vector<double> aperp_values;
    std::vector<double> sigma_x_values;
    std::vector<double> sigma_perp_values;
    std::string result;
    Args() : profiles("G0,Gx,Gp,G2") {}
};

bool parse_double_list(const std::string& text, std::vector<double>& out)
{
    out.clear();
    std::istringstream iss(text);
    std::string token;
    while (std::getline(iss, token, ',')) {
        if (token.empty()) return false;
        char* end = NULL;
        const double value = std::strtod(token.c_str(), &end);
        if (!end || *end != '\0' || !std::isfinite(value)) return false;
        out.push_back(value);
    }
    return !out.empty();
}

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--input" && i + 1 < argc) args.input = argv[++i];
        else if (arg == "--profiles" && i + 1 < argc) args.profiles = argv[++i];
        else if (arg == "--ax-values" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], args.ax_values)) return false;
        } else if (arg == "--aperp-values" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], args.aperp_values)) return false;
        } else if (arg == "--sigma-x-cells" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], args.sigma_x_values)) return false;
        } else if (arg == "--sigma-perp-cells" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], args.sigma_perp_values))
                return false;
        } else if (arg == "--result" && i + 1 < argc) args.result = argv[++i];
        else return false;
    }
    return !args.input.empty() && !args.profiles.empty();
}

std::vector<std::string> split_profiles(const std::string& text)
{
    std::vector<std::string> out;
    std::istringstream iss(text);
    std::string token;
    while (std::getline(iss, token, ',')) {
        if (!token.empty()) out.push_back(token);
    }
    return out;
}

std::vector<BulkTailVelocityBinAudit> aggregate_rows(
    const std::vector<HistogramRow>& rows)
{
    std::vector<BulkTailVelocityBinAudit> bins;
    for (size_t r = 0; r < rows.size(); ++r) {
        bool merged = false;
        for (size_t q = 0; q < bins.size(); ++q) {
            if (bins[q].iv == rows[r].iv && bins[q].imu == rows[r].imu) {
                bins[q].request_number += rows[r].request_number;
                bins[q].request_cell_count +=
                    static_cast<unsigned long long>(rows[r].request_cell_count);
                merged = true;
                break;
            }
        }
        if (!merged) {
            BulkTailVelocityBinAudit bin;
            bin.iv = rows[r].iv;
            bin.imu = rows[r].imu;
            bin.request_cell_count =
                static_cast<unsigned long long>(rows[r].request_cell_count);
            bin.request_number = rows[r].request_number;
            bins.push_back(bin);
        }
    }
    std::sort(bins.begin(), bins.end(),
              [](const BulkTailVelocityBinAudit& a,
                 const BulkTailVelocityBinAudit& b) {
                  if (a.iv != b.iv) return a.iv < b.iv;
                  return a.imu < b.imu;
              });
    return bins;
}

// FNV-1a over the raw bytes of a face array (deterministic bit identity).
unsigned long long face_hash(const std::vector<double>& faces)
{
    unsigned long long hash = 0xcbf29ce484222325ULL;
    const unsigned long long prime = 0x100000001b3ULL;
    const unsigned char* bytes =
        reinterpret_cast<const unsigned char*>(faces.data());
    for (size_t i = 0; i < faces.size() * sizeof(double); ++i) {
        hash ^= bytes[i];
        hash *= prime;
    }
    return hash;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    if (!parse_args(argc, argv, args)) {
        if (rank == 0)
            std::cerr
                << "usage: tail_interface_grid_replay_mpi_test "
                   "--input <hist> --profiles G0,Gx,Gp,G2 "
                   "--ax-values a,b --aperp-values a,b "
                   "--sigma-x-cells a,b --sigma-perp-cells a,b "
                   "[--result path]\n";
        MPI_Finalize();
        return 2;
    }
    const std::vector<std::string> profiles = split_profiles(args.profiles);
    const std::vector<double> ax_values =
        args.ax_values.empty() ? std::vector<double>(1, 0.0) : args.ax_values;
    const std::vector<double> aperp_values =
        args.aperp_values.empty() ? std::vector<double>(1, 0.0)
                                  : args.aperp_values;
    const std::vector<double> sigma_x_values =
        args.sigma_x_values.empty() ? std::vector<double>(1, 1.0)
                                    : args.sigma_x_values;
    const std::vector<double> sigma_perp_values =
        args.sigma_perp_values.empty() ? std::vector<double>(1, 1.0)
                                       : args.sigma_perp_values;

    std::vector<HistogramRow> all_rows;
    std::string error;
    if (!read_velocity_histogram_rows(args.input, all_rows, error)) {
        if (rank == 0) std::cerr << error << "\n";
        MPI_Finalize();
        return 1;
    }
    // Single-rank reference from the full file (file order).
    const std::vector<BulkTailVelocityBinAudit> serial_bins =
        aggregate_rows(all_rows);
    // Deterministic row partition across ranks.
    std::vector<HistogramRow> local_rows;
    for (size_t r = 0; r < all_rows.size(); ++r) {
        if (static_cast<int>(r % static_cast<size_t>(size)) == rank)
            local_rows.push_back(all_rows[r]);
    }
    const std::vector<BulkTailVelocityBinAudit> local_bins =
        aggregate_rows(local_rows);

    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    const int nv = static_cast<int>(cgrid.upar_cells.size());
    const int nmu = static_cast<int>(cgrid.uperp_cells.size());
    const size_t nslots = static_cast<size_t>(nv) * nmu;

    // Global per-bin reduction (SUM over the row partition).
    std::vector<double> local_number(nslots, 0.0);
    std::vector<unsigned long long> local_count(nslots, 0ULL);
    for (size_t b = 0; b < local_bins.size(); ++b) {
        const size_t key = static_cast<size_t>(local_bins[b].iv) * nmu +
                           static_cast<size_t>(local_bins[b].imu);
        local_number[key] = local_bins[b].request_number;
        local_count[key] +=
            static_cast<unsigned long long>(local_bins[b].request_cell_count);
    }
    std::vector<double> global_number(nslots, 0.0);
    std::vector<unsigned long long> global_count(nslots, 0ULL);
    MPI_Allreduce(local_number.data(), global_number.data(),
                  static_cast<int>(nslots), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_count.data(), global_count.data(),
                  static_cast<int>(nslots), MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    std::vector<BulkTailVelocityBinAudit> global_bins;
    for (int j = 0; j < nv; ++j) {
        for (int k = 0; k < nmu; ++k) {
            const size_t key = static_cast<size_t>(j) * nmu + k;
            if (global_count[key] <= 0 || !(global_number[key] > 0.0)) continue;
            BulkTailVelocityBinAudit bin;
            bin.iv = j;
            bin.imu = k;
            bin.request_cell_count =
                static_cast<unsigned long long>(global_count[key]);
            bin.request_number = global_number[key];
            global_bins.push_back(bin);
        }
    }
    const std::vector<double> serial_agg =
        tail_interface_aggregate_histogram(serial_bins, nv, nmu);
    const std::vector<double> local_agg =
        tail_interface_aggregate_histogram(local_bins, nv, nmu);
    const std::vector<double> global_agg =
        tail_interface_aggregate_histogram(global_bins, nv, nmu);

    // Global bin consistency against the single-rank reference.
    bool bins_consistent = global_bins.size() == serial_bins.size();
    double bins_max_rel = 0.0;
    if (bins_consistent) {
        for (size_t s = 0; s < nslots; ++s) {
            const double scale = std::max(1.0, std::fabs(serial_agg[s]));
            bins_max_rel =
                std::max(bins_max_rel,
                         std::fabs(global_agg[s] - serial_agg[s]) / scale);
        }
        bins_consistent = bins_max_rel <= 1.0e-12;
    }

    bool program_ok = bins_consistent;
    bool g0_identity_failed = false;
    bool g0_seen = false;
    double serial_total = 0.0;
    for (size_t s = 0; s < serial_agg.size(); ++s) serial_total += serial_agg[s];
    std::ostringstream report;
    if (rank == 0) {
        // Event count from the unique (step, location) pairs.
        long long events = 0;
        long long last_step = -1, last_loc = -1;
        for (size_t r = 0; r < all_rows.size(); ++r) {
            if (r == 0 || all_rows[r].accepted_step != last_step ||
                all_rows[r].conversion_location != last_loc) {
                ++events;
                last_step = all_rows[r].accepted_step;
                last_loc = all_rows[r].conversion_location;
            }
        }
        report << "histogram_input=" << args.input
               << " histogram_events=" << events
               << " histogram_bins=" << static_cast<int>(serial_bins.size())
               << " histogram_number=" << std::setprecision(17)
               << serial_total << " ranks=" << size
               << " bins_max_rel=" << bins_max_rel << "\n";
    }

    for (size_t p = 0; p < profiles.size(); ++p) {
        const std::string& profile = profiles[p];
        size_t combo_count = 0;
        if (profile == "G0") combo_count = 1;
        else if (profile == "Gx") combo_count = ax_values.size() * sigma_x_values.size();
        else if (profile == "Gp") combo_count = aperp_values.size() * sigma_perp_values.size();
        else if (profile == "G2")
            combo_count = ax_values.size() * aperp_values.size() *
                          sigma_x_values.size() * sigma_perp_values.size();
        else {
            if (rank == 0) std::cerr << "unknown profile " << profile << "\n";
            MPI_Finalize();
            return 2;
        }
        for (size_t c = 0; c < combo_count; ++c) {
            size_t rem = c;
            double ax = 0.0, aperp = 0.0, sx = 1.0, sp = 1.0;
            if (profile == "G0") {
                ax = 0.0; aperp = 0.0; sx = 1.0; sp = 1.0;
            } else if (profile == "Gx") {
                ax = ax_values[rem % ax_values.size()]; rem /= ax_values.size();
                sx = sigma_x_values[rem % sigma_x_values.size()];
            } else if (profile == "Gp") {
                aperp = aperp_values[rem % aperp_values.size()];
                rem /= aperp_values.size();
                sp = sigma_perp_values[rem % sigma_perp_values.size()];
            } else {
                ax = ax_values[rem % ax_values.size()]; rem /= ax_values.size();
                aperp = aperp_values[rem % aperp_values.size()];
                rem /= aperp_values.size();
                sx = sigma_x_values[rem % sigma_x_values.size()];
                rem /= sigma_x_values.size();
                sp = sigma_perp_values[rem % sigma_perp_values.size()];
            }
            TailInterfaceGridDesignConfig config;
            config.ax = ax;
            config.aperp = aperp;
            config.sigma_x_cells = sx;
            config.sigma_perp_cells = sp;
            const TailInterfaceGridCandidate candidate =
                build_tail_interface_grid_candidate(cgrid, global_bins, config,
                                                    profile);
            if (profile == "G0") g0_seen = true;

            // Rank-local remap from the partitioned input, then the global
            // reduction; compare with the single-rank reference remap.
            double partition_error = 0.0;
            const std::vector<double> local_remap =
                tail_interface_remap_masses(
                    cgrid.upar_faces, cgrid.uperp_faces,
                    candidate.upar_faces, candidate.uperp_faces, local_agg,
                    &partition_error);
            std::vector<double> reduced_remap(nslots, 0.0);
            MPI_Allreduce(local_remap.data(), reduced_remap.data(),
                          static_cast<int>(nslots), MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);
            const std::vector<double> reference_remap =
                tail_interface_remap_masses(
                    cgrid.upar_faces, cgrid.uperp_faces,
                    candidate.upar_faces, candidate.uperp_faces, serial_agg,
                    NULL);
            double remap_max_rel = 0.0;
            double ref_max = 0.0;
            for (size_t s = 0; s < nslots; ++s)
                ref_max = std::max(ref_max, std::fabs(reference_remap[s]));
            for (size_t s = 0; s < nslots; ++s) {
                const double scale =
                    std::max(1.0e-30 * ref_max, std::fabs(reference_remap[s]));
                remap_max_rel = std::max(
                    remap_max_rel,
                    std::fabs(reduced_remap[s] - reference_remap[s]) / scale);
            }

            // Bit-identical faces across ranks.
            const unsigned long long upar_hash = face_hash(candidate.upar_faces);
            const unsigned long long uperp_hash = face_hash(candidate.uperp_faces);
            unsigned long long upar_min = 0, upar_max = 0;
            unsigned long long uperp_min = 0, uperp_max = 0;
            MPI_Allreduce(&upar_hash, &upar_min, 1, MPI_UNSIGNED_LONG_LONG,
                          MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&upar_hash, &upar_max, 1, MPI_UNSIGNED_LONG_LONG,
                          MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&uperp_hash, &uperp_min, 1, MPI_UNSIGNED_LONG_LONG,
                          MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&uperp_hash, &uperp_max, 1, MPI_UNSIGNED_LONG_LONG,
                          MPI_MAX, MPI_COMM_WORLD);
            const bool faces_identical =
                upar_min == upar_max && uperp_min == uperp_max;
            const bool combo_ok = faces_identical && remap_max_rel <= 1.0e-10;

            int local_ok = combo_ok ? 1 : 0;
            int all_ok = 0;
            MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN,
                          MPI_COMM_WORLD);
            program_ok = program_ok && all_ok == 1;

            if (rank == 0) {
                if (profile == "G0" && !g0_identity_failed) {
                    const TailInterfaceReplayResult result =
                        replay_tail_interface_histogram(
                            cgrid, candidate, global_bins, 6.0e6 * Const::eV);
                    program_ok = program_ok && candidate.valid &&
                                 result.g0_identity_ok &&
                                 result.max_partition_error <= 1.0e-14;
                    g0_identity_failed = !result.g0_identity_ok;
                    write_candidate_block(report, candidate.grid_name, result);
                } else if (g0_identity_failed) {
                    TailInterfaceReplayResult invalid;
                    invalid.candidate_status = "INVALID";
                    invalid.status_reason = "G0 identity failed";
                    write_candidate_block(report, candidate.grid_name, invalid);
                } else {
                    const TailInterfaceReplayResult result =
                        replay_tail_interface_histogram(
                            cgrid, candidate, global_bins, 6.0e6 * Const::eV);
                    if (result.max_partition_error > 1.0e-14)
                        program_ok = false;
                    write_candidate_block(report, candidate.grid_name, result);
                }
            }
        }
    }
    if (!g0_seen) {
        if (rank == 0) std::cerr << "profiles must contain G0 (identity gate)\n";
        MPI_Finalize();
        return 2;
    }
    int local_pass = program_ok ? 1 : 0;
    int all_pass = 0;
    MPI_Allreduce(&local_pass, &all_pass, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        report << "status=" << (all_pass ? "PASS" : "FAIL") << "\n";
        std::cout << report.str();
        if (!write_result_file(args.result, report.str())) {
            std::cerr << "cannot write result file: " << args.result << "\n";
            MPI_Finalize();
            return 1;
        }
    }
    MPI_Finalize();
    return all_pass ? 0 : 1;
}
