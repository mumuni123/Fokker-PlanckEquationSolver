// Section 7.11.16C / 16B: serial replay of the real request velocity
// histogram on the fixed-(192,64) candidate grids.  The file-end status
// only reports program correctness (histogram parsing, conservative
// overlap, G0 identity); each candidate carries its own
// candidate_status=GREEN|GRAY|RED|INVALID for the 16B offline gates.
//
// Usage:
//   tail_interface_grid_replay_test
//     --input <velocity_histogram.dat> --profiles G0,Gx,Gp,G2
//     --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0
//     --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 [--result path]

#include "grid.h"
#include "tail_interface_grid_design.h"
#include "tail_interface_replay_common.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

// One candidate block: build + replay + serialize.
void run_one(const CylindricalVelocityGrid& cgrid,
             const std::vector<BulkTailVelocityBinAudit>& bins,
             const std::string& profile, const TailInterfaceGridDesignConfig& config,
             bool g0_identity_failed, std::ostringstream& out,
             bool& program_ok)
{
    const TailInterfaceGridCandidate candidate =
        build_tail_interface_grid_candidate(cgrid, bins, config, profile);
    if (g0_identity_failed) {
        TailInterfaceReplayResult invalid;
        invalid.candidate_status = "INVALID";
        invalid.status_reason = "G0 identity failed";
        write_candidate_block(out, candidate.grid_name, invalid);
        return;
    }
    const TailInterfaceReplayResult result = replay_tail_interface_histogram(
        cgrid, candidate, bins, 6.0e6 * Const::eV);
    if (result.max_partition_error > 1.0e-14) program_ok = false;
    write_candidate_block(out, candidate.grid_name, result);
}

} // namespace

int main(int argc, char** argv)
{
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr
            << "usage: tail_interface_grid_replay_test --input <hist> "
               "--profiles G0,Gx,Gp,G2 --ax-values a,b --aperp-values a,b "
               "--sigma-x-cells a,b --sigma-perp-cells a,b [--result path]\n";
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

    HistogramFileData data;
    std::string error;
    if (!read_velocity_histogram(args.input, data, error)) {
        std::cerr << error << "\n";
        return 1;
    }
    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);

    std::ostringstream report;
    report << "histogram_input=" << args.input
           << " histogram_events=" << data.events
           << " histogram_cell_requests=" << data.cell_requests
           << " histogram_bins=" << data.bins_count
           << " histogram_number=" << std::setprecision(17)
           << data.total_number << "\n";

    // G0 first: identity failure invalidates every following candidate.
    bool g0_identity_failed = false;
    bool program_ok = true;
    bool g0_seen = false;
    for (size_t p = 0; p < profiles.size(); ++p) {
        const std::string& profile = profiles[p];
        if (profile == "G0") {
            g0_seen = true;
            TailInterfaceGridDesignConfig config;
            config.ax = 0.0;
            config.aperp = 0.0;
            const TailInterfaceGridCandidate candidate =
                build_tail_interface_grid_candidate(cgrid, data.bins, config,
                                                    "G0");
            const TailInterfaceReplayResult result =
                replay_tail_interface_histogram(cgrid, candidate, data.bins,
                                                6.0e6 * Const::eV);
            program_ok = program_ok && candidate.valid && result.g0_identity_ok;
            g0_identity_failed = !result.g0_identity_ok;
            write_candidate_block(report, candidate.grid_name, result);
            continue;
        }
        if (profile == "Gx") {
            for (size_t a = 0; a < ax_values.size(); ++a) {
                for (size_t s = 0; s < sigma_x_values.size(); ++s) {
                    TailInterfaceGridDesignConfig config;
                    config.ax = ax_values[a];
                    config.aperp = 0.0;
                    config.sigma_x_cells = sigma_x_values[s];
                    config.sigma_perp_cells = 1.0;
                    run_one(cgrid, data.bins, profile, config,
                            g0_identity_failed, report, program_ok);
                }
            }
            continue;
        }
        if (profile == "Gp") {
            for (size_t a = 0; a < aperp_values.size(); ++a) {
                for (size_t s = 0; s < sigma_perp_values.size(); ++s) {
                    TailInterfaceGridDesignConfig config;
                    config.ax = 0.0;
                    config.aperp = aperp_values[a];
                    config.sigma_x_cells = 1.0;
                    config.sigma_perp_cells = sigma_perp_values[s];
                    run_one(cgrid, data.bins, profile, config,
                            g0_identity_failed, report, program_ok);
                }
            }
            continue;
        }
        if (profile == "G2") {
            for (size_t a = 0; a < ax_values.size(); ++a) {
                for (size_t ap = 0; ap < aperp_values.size(); ++ap) {
                    for (size_t s = 0; s < sigma_x_values.size(); ++s) {
                        for (size_t sp = 0; sp < sigma_perp_values.size(); ++sp) {
                            TailInterfaceGridDesignConfig config;
                            config.ax = ax_values[a];
                            config.aperp = aperp_values[ap];
                            config.sigma_x_cells = sigma_x_values[s];
                            config.sigma_perp_cells = sigma_perp_values[sp];
                            run_one(cgrid, data.bins, profile, config,
                                    g0_identity_failed, report, program_ok);
                        }
                    }
                }
            }
            continue;
        }
        std::cerr << "unknown profile " << profile << "\n";
        return 2;
    }
    if (!g0_seen) {
        std::cerr << "profiles must contain G0 (identity gate)\n";
        return 2;
    }
    report << "status=" << (program_ok ? "PASS" : "FAIL") << "\n";
    std::cout << report.str();
    if (!write_result_file(args.result, report.str())) {
        std::cerr << "cannot write result file: " << args.result << "\n";
        return 1;
    }
    return program_ok ? 0 : 1;
}
