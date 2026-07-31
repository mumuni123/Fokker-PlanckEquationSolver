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

enum RuntimeMidpointInitialGuessMode {
    RUNTIME_MIDPOINT_INITIAL_GUESS_NONE = 0,
    RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR = 1
};

struct RuntimeOptions {
    std::string output_dir;
    std::string restart_dir;
    std::string operator_audit_dir;
    std::string beam_ledger_reference;
    std::string remap_checkpoint_source;
    std::string remap_checkpoint_destination;
    std::vector<double> checkpoint_times_fs;
    double stop_time_fs;
    long long stop_after_accepted_steps;
    double dt_scale;
    double midpoint_trace_start_fs;
    double midpoint_trace_end_fs;
    int midpoint_max_iters;
    int diagnostic_level;
    // Final accepted-state transport/energy ledger.  A negative window start
    // disables time filtering; cadence controls output/collection only.
    int accepted_energy_audit_cadence;
    double accepted_energy_audit_start_fs;
    double accepted_energy_audit_end_fs;
    // Calibration values for the transactional soft-candidate gate.  These
    // are runtime controls rather than hidden production constants.
    double soft_candidate_field_tolerance;
    double soft_candidate_current_tolerance;
    double soft_candidate_energy_p99_reference;
    double soft_candidate_energy_absolute_limit;
    RuntimeMidpointInitialGuessMode midpoint_initial_guess_mode;
    RuntimeMidpointAccelerationMode midpoint_acceleration_mode;
    int anderson_depth;
    int acceleration_start_iter;
    double acceleration_accept_ratio;
    double acceleration_max_coefficient;
    BeamLedgerMode beam_ledger_mode;
    // 0 = explicit legacy fallback; 1 = default dual-u coupling.
    int background_coupling_mode;
    // 0 = stable cell-capacity baseline; 1 = experimental regularized face
    // pairing with mandatory monotone fallback.
    int face_pairing_mode;
    double face_pairing_sigma_cutoff;
    double face_pairing_lambda;
    double face_pairing_eta;
    // This is only the fraction of the local dual-u capacity exposed to the
    // constrained solve.  Acceptance has separate amplitude/energy/state
    // tolerances below.
    double face_pairing_trust_fraction;
    double face_pairing_correction_trust_fraction;
    double face_pairing_energy_pair_tolerance;
    double face_pairing_energy_residual_fraction;
    double face_pairing_mass_relative_tolerance;
    double face_pairing_f_residual_growth_tolerance;
    bool dump_final_midpoint;
    bool overwrite_output;
    bool checkpoint_enabled;
    bool restart_enabled;
    bool operator_audit_mode;
    bool beam_ledger_reference_enabled;
    bool velocity_grid_remap_enabled;
};

RuntimeOptions parse_runtime_options(int argc, char** argv, int mpi_rank,
                                     int mpi_size);
std::string output_path(const RuntimeOptions& options,
                        const std::string& filename);
bool prepare_output_directory(const RuntimeOptions& options, int mpi_rank,
                              int mpi_size, std::string& error);

#endif
