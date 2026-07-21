#ifndef RUNTIME_OPTIONS_H
#define RUNTIME_OPTIONS_H

#include <string>
#include <vector>

struct RuntimeOptions {
    std::string output_dir;
    std::string restart_dir;
    std::string operator_audit_dir;
    std::string beam_ledger_reference;
    std::vector<double> checkpoint_times_fs;
    double stop_time_fs;
    long long stop_after_accepted_steps;
    double dt_scale;
    int midpoint_max_iters;
    int diagnostic_level;
    // 0 = legacy charge/energy coupling; 1 = opt-in dual-u coupling.
    int background_coupling_mode;
    bool dump_final_midpoint;
    bool overwrite_output;
    bool checkpoint_enabled;
    bool restart_enabled;
    bool operator_audit_mode;
    bool beam_ledger_reference_enabled;
};

RuntimeOptions parse_runtime_options(int argc, char** argv, int mpi_rank,
                                     int mpi_size);
std::string output_path(const RuntimeOptions& options,
                        const std::string& filename);
bool prepare_output_directory(const RuntimeOptions& options, int mpi_rank,
                              int mpi_size, std::string& error);

#endif
