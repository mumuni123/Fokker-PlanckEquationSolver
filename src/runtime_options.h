#ifndef RUNTIME_OPTIONS_H
#define RUNTIME_OPTIONS_H

#include <string>
#include <vector>

enum BeamLedgerMode {
    BEAM_LEDGER_OFF = 0,
    BEAM_LEDGER_SUMMARY = 1,
    BEAM_LEDGER_FULL = 2
};

enum RuntimeMidpointAccelerationMode {
    RUNTIME_MIDPOINT_ACCELERATION_NONE = 0,
    RUNTIME_MIDPOINT_ACCELERATION_AITKEN = 1,
    RUNTIME_MIDPOINT_ACCELERATION_ANDERSON = 2
};

struct RuntimeOptions {
    std::string output_dir;
    std::string restart_dir;
    std::string operator_audit_dir;
    std::string beam_ledger_reference;
    std::vector<double> checkpoint_times_fs;
    double stop_time_fs;
    long long stop_after_accepted_steps;
    double dt_scale;
    double midpoint_trace_start_fs;
    double midpoint_trace_end_fs;
    int midpoint_max_iters;
    int diagnostic_level;
    RuntimeMidpointAccelerationMode midpoint_acceleration_mode;
    int anderson_depth;
    int acceleration_start_iter;
    double acceleration_accept_ratio;
    double acceleration_max_coefficient;
    BeamLedgerMode beam_ledger_mode;
    // 0 = explicit legacy fallback; 1 = default dual-u coupling.
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
