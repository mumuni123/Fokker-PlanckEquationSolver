#include "runtime_options.h"

#include "parameters.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <mpi.h>
#include <sstream>
#include <sys/stat.h>
#if defined(_WIN32)
#include <direct.h>
#endif

namespace {
RuntimeOptions defaults()
{
    RuntimeOptions o;
    o.output_dir = "output";
    o.stop_time_fs = Param::t_end / Const::femto;
    o.stop_after_accepted_steps = -1;
    o.dt_scale = 1.0;
    o.midpoint_trace_start_fs = -1.0;
    o.midpoint_trace_end_fs = -1.0;
    o.midpoint_max_iters = 40;
    o.diagnostic_level = 1;
    o.accepted_energy_audit_cadence = 500;
    o.accepted_energy_audit_start_fs = -1.0;
    o.accepted_energy_audit_end_fs = -1.0;
    o.soft_candidate_field_tolerance = 3.0e-6;
    o.soft_candidate_current_tolerance = 3.0e-5;
    o.soft_candidate_energy_p99_reference = 1.0e-3;
    o.soft_candidate_energy_absolute_limit = 5.0e-3;
    o.midpoint_initial_guess_mode = RUNTIME_MIDPOINT_INITIAL_GUESS_NONE;
    o.midpoint_acceleration_mode = RUNTIME_MIDPOINT_ACCELERATION_NONE;
    o.anderson_depth = 3;
    o.acceleration_start_iter = 3;
    o.acceleration_accept_ratio = 0.95;
    o.acceleration_max_coefficient = 2.0;
    o.beam_ledger_mode = BEAM_LEDGER_SUMMARY;
    o.background_coupling_mode = 1;
    o.face_pairing_mode = 0;
    o.face_pairing_sigma_cutoff = 1.0e-8;
    o.face_pairing_lambda = 1.0e-3;
    o.face_pairing_eta = 1.0e-8;
    o.face_pairing_trust_fraction = 0.1;
    o.face_pairing_correction_trust_fraction = 1.0;
    o.face_pairing_energy_pair_tolerance = 1.0e-8;
    o.face_pairing_energy_residual_fraction = 1.0;
    o.face_pairing_mass_relative_tolerance = 1.0e-10;
    o.face_pairing_f_residual_growth_tolerance = 0.1;
    o.dump_final_midpoint = false;
    o.overwrite_output = false;
    o.checkpoint_enabled = false;
    o.restart_enabled = false;
    o.operator_audit_mode = false;
    o.beam_ledger_reference_enabled = false;
    o.velocity_grid_remap_enabled = false;
    return o;
}

bool parse_times(const char* text, std::vector<double>& values)
{
    std::stringstream stream(text ? text : "");
    std::string item;
    while (std::getline(stream, item, ',')) {
        char* end = 0;
        const double value = std::strtod(item.c_str(), &end);
        if (item.empty() || !end || *end != '\0' || !(value >= 0.0)) return false;
        values.push_back(value);
    }
    std::sort(values.begin(), values.end());
    return !values.empty();
}

void broadcast_string(std::string& value, int root)
{
    int size = static_cast<int>(value.size());
    MPI_Bcast(&size, 1, MPI_INT, root, MPI_COMM_WORLD);
    if (size < 0) size = 0;
    if (static_cast<int>(value.size()) != size) value.assign(static_cast<size_t>(size), '\0');
    if (size > 0) MPI_Bcast(&value[0], size, MPI_CHAR, root, MPI_COMM_WORLD);
}

bool directory_has_entries(const std::string& path)
{
    DIR* dir = opendir(path.c_str());
    if (!dir) return false;
    bool nonempty = false;
    for (dirent* entry = readdir(dir); entry; entry = readdir(dir)) {
        if (std::strcmp(entry->d_name, ".") && std::strcmp(entry->d_name, "..")) {
            nonempty = true;
            break;
        }
    }
    closedir(dir);
    return nonempty;
}
}

RuntimeOptions parse_runtime_options(int argc, char** argv, int mpi_rank,
                                     int mpi_size)
{
    RuntimeOptions options = defaults();
    std::string error;
    bool beam_ledger_mode_explicit = false;
    if (mpi_rank == 0) {
        for (int i = 1; i < argc && error.empty(); ++i) {
            const std::string arg(argv[i]);
            const bool has_value = i + 1 < argc;
            const char* value = has_value ? argv[i + 1] : 0;
            if (arg == "--dump-final-midpoint") options.dump_final_midpoint = true;
            else if (arg == "--overwrite-output") options.overwrite_output = true;
            else if (arg == "--restart" && has_value) {
                options.restart_dir = value; options.restart_enabled = true; ++i;
            } else if (arg == "--remap-checkpoint-velocity-grid" && i + 2 < argc) {
                options.remap_checkpoint_source = argv[++i];
                options.remap_checkpoint_destination = argv[++i];
                options.velocity_grid_remap_enabled = true;
            } else if (arg == "--remap-checkpoint-velocity-grid") {
                error = "--remap-checkpoint-velocity-grid requires <old_dir> <new_dir>";
            } else if (arg == "--operator-audit" && has_value) {
                options.operator_audit_dir = value; options.operator_audit_mode = true; ++i;
            } else if (arg == "--beam-ledger-reference" && has_value) {
                options.beam_ledger_reference = value;
                options.beam_ledger_reference_enabled = true;
                ++i;
            } else if (arg == "--beam-ledger-mode" && has_value) {
                const std::string mode(value);
                if (mode == "off") options.beam_ledger_mode = BEAM_LEDGER_OFF;
                else if (mode == "summary") options.beam_ledger_mode = BEAM_LEDGER_SUMMARY;
                else if (mode == "full") options.beam_ledger_mode = BEAM_LEDGER_FULL;
                else error = "invalid --beam-ledger-mode (use off, summary, or full)";
                beam_ledger_mode_explicit = true;
                ++i;
            } else if (arg == "--checkpoint-times" && has_value) {
                options.checkpoint_enabled = parse_times(value, options.checkpoint_times_fs); ++i;
                if (!options.checkpoint_enabled) error = "invalid --checkpoint-times";
            } else if (arg == "--stop-time-fs" && has_value) {
                options.stop_time_fs = std::strtod(value, 0); ++i;
            } else if (arg == "--stop-after-steps" && has_value) {
                options.stop_after_accepted_steps = std::strtoll(value, 0, 10); ++i;
            } else if (arg == "--dt-scale" && has_value) {
                options.dt_scale = std::strtod(value, 0); ++i;
            } else if (arg == "--midpoint-trace-window-fs" && has_value) {
                std::vector<double> window;
                if (!parse_times(value, window) || window.size() != 2) {
                    error = "invalid --midpoint-trace-window-fs (use start,end)";
                } else {
                    options.midpoint_trace_start_fs = window[0];
                    options.midpoint_trace_end_fs = window[1];
                }
                ++i;
            } else if (arg == "--midpoint-max-iters" && has_value) {
                options.midpoint_max_iters = std::atoi(value); ++i;
            } else if (arg == "--midpoint-initial-guess" && has_value) {
                const std::string mode(value);
                if (mode == "none") {
                    options.midpoint_initial_guess_mode =
                        RUNTIME_MIDPOINT_INITIAL_GUESS_NONE;
                } else if (mode == "field-linear") {
                    options.midpoint_initial_guess_mode =
                        RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR;
                } else {
                    error = "invalid --midpoint-initial-guess "
                            "(use none or field-linear)";
                }
                ++i;
            } else if (arg == "--midpoint-acceleration" && has_value) {
                const std::string mode(value);
                if (mode == "none") {
                    options.midpoint_acceleration_mode = RUNTIME_MIDPOINT_ACCELERATION_NONE;
                } else if (mode == "aitken") {
                    options.midpoint_acceleration_mode = RUNTIME_MIDPOINT_ACCELERATION_AITKEN;
                } else if (mode == "anderson") {
                    options.midpoint_acceleration_mode = RUNTIME_MIDPOINT_ACCELERATION_ANDERSON;
                } else {
                    error = "invalid --midpoint-acceleration (use none, aitken, or anderson)";
                }
                ++i;
            } else if (arg == "--anderson-depth" && has_value) {
                options.anderson_depth = std::atoi(value); ++i;
            } else if (arg == "--acceleration-start-iter" && has_value) {
                options.acceleration_start_iter = std::atoi(value); ++i;
            } else if (arg == "--acceleration-accept-ratio" && has_value) {
                options.acceleration_accept_ratio = std::strtod(value, 0); ++i;
            } else if (arg == "--acceleration-max-coefficient" && has_value) {
                options.acceleration_max_coefficient = std::strtod(value, 0); ++i;
            } else if (arg == "--diagnostic-level" && has_value) {
                options.diagnostic_level = std::atoi(value); ++i;
            } else if (arg == "--accepted-energy-audit-cadence" && has_value) {
                options.accepted_energy_audit_cadence = std::atoi(value); ++i;
            } else if (arg == "--accepted-energy-audit-window-fs" && has_value) {
                std::vector<double> window;
                if (!parse_times(value, window) || window.size() != 2) {
                    error = "invalid --accepted-energy-audit-window-fs (use start,end)";
                } else {
                    options.accepted_energy_audit_start_fs = window[0];
                    options.accepted_energy_audit_end_fs = window[1];
                }
                ++i;
            } else if (arg == "--soft-candidate-field-tol" && has_value) {
                options.soft_candidate_field_tolerance = std::strtod(value, 0); ++i;
            } else if (arg == "--soft-candidate-current-tol" && has_value) {
                options.soft_candidate_current_tolerance = std::strtod(value, 0); ++i;
            } else if (arg == "--soft-candidate-energy-p99" && has_value) {
                options.soft_candidate_energy_p99_reference = std::strtod(value, 0); ++i;
            } else if (arg == "--soft-candidate-energy-absolute-limit" && has_value) {
                options.soft_candidate_energy_absolute_limit = std::strtod(value, 0); ++i;
            } else if (arg == "--background-coupling-mode" && has_value) {
                const std::string mode(value);
                if (mode == "legacy") options.background_coupling_mode = 0;
                else if (mode == "dual-u" || mode == "dual_u")
                    options.background_coupling_mode = 1;
                else error = "invalid --background-coupling-mode (use legacy or dual-u)";
                ++i;
            } else if (arg == "--face-pairing-mode" && has_value) {
                const std::string mode(value);
                if (mode == "cell-baseline")
                    options.face_pairing_mode = 0;
                else if (mode == "regularized")
                    options.face_pairing_mode = 1;
                else
                    error = "invalid --face-pairing-mode "
                            "(use cell-baseline or regularized)";
                ++i;
            } else if (arg == "--face-pairing-sigma-cutoff" && has_value) {
                options.face_pairing_sigma_cutoff = std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-lambda" && has_value) {
                options.face_pairing_lambda = std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-eta" && has_value) {
                options.face_pairing_eta = std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-trust-fraction" && has_value) {
                options.face_pairing_trust_fraction = std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-correction-trust-fraction" &&
                       has_value) {
                options.face_pairing_correction_trust_fraction =
                    std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-energy-pair-tolerance" &&
                       has_value) {
                options.face_pairing_energy_pair_tolerance =
                    std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-energy-residual-fraction" &&
                       has_value) {
                options.face_pairing_energy_residual_fraction =
                    std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-mass-relative-tolerance" &&
                       has_value) {
                options.face_pairing_mass_relative_tolerance =
                    std::strtod(value, 0);
                ++i;
            } else if (arg == "--face-pairing-f-residual-growth-tolerance" &&
                       has_value) {
                options.face_pairing_f_residual_growth_tolerance =
                    std::strtod(value, 0);
                ++i;
            } else if (arg == "--output-dir" && has_value) {
                options.output_dir = value; ++i;
            } else error = "unknown or incomplete runtime option: " + arg;
        }
        if (options.dt_scale <= 0.0 || options.midpoint_max_iters < 1 ||
            options.diagnostic_level < 0 || options.diagnostic_level > 2 ||
            options.accepted_energy_audit_cadence < 1 ||
            !(options.soft_candidate_field_tolerance > 0.0) ||
            !(options.soft_candidate_current_tolerance > 0.0) ||
            !(options.soft_candidate_energy_p99_reference > 0.0) ||
            !(options.soft_candidate_energy_absolute_limit > 0.0) ||
            ((options.accepted_energy_audit_start_fs >= 0.0 ||
              options.accepted_energy_audit_end_fs >= 0.0) &&
             (!(options.accepted_energy_audit_start_fs >= 0.0) ||
              options.accepted_energy_audit_end_fs <
                  options.accepted_energy_audit_start_fs)) ||
            options.midpoint_acceleration_mode < RUNTIME_MIDPOINT_ACCELERATION_NONE ||
            options.midpoint_acceleration_mode > RUNTIME_MIDPOINT_ACCELERATION_ANDERSON ||
            options.midpoint_initial_guess_mode <
                RUNTIME_MIDPOINT_INITIAL_GUESS_NONE ||
            options.midpoint_initial_guess_mode >
                RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR ||
            (options.anderson_depth != 2 && options.anderson_depth != 3) ||
            options.acceleration_start_iter < 1 ||
            !(options.acceleration_accept_ratio > 0.0 &&
              options.acceleration_accept_ratio <= 1.0) ||
            !(options.acceleration_max_coefficient > 0.0) ||
            options.stop_time_fs < 0.0 || options.stop_after_accepted_steps < -1 ||
            ((options.midpoint_trace_start_fs >= 0.0 ||
              options.midpoint_trace_end_fs >= 0.0) &&
             (!(options.midpoint_trace_start_fs >= 0.0) ||
              options.midpoint_trace_end_fs < options.midpoint_trace_start_fs)) ||
            options.background_coupling_mode < 0 ||
            options.background_coupling_mode > 1 ||
            options.face_pairing_mode < 0 ||
            options.face_pairing_mode > 1 ||
            (options.face_pairing_mode == 1 &&
             options.background_coupling_mode != 1) ||
            !(options.face_pairing_sigma_cutoff >= 0.0 &&
              options.face_pairing_sigma_cutoff < 1.0) ||
            !(options.face_pairing_lambda >= 0.0) ||
            !(options.face_pairing_eta > 0.0) ||
            !(options.face_pairing_trust_fraction > 0.0 &&
              options.face_pairing_trust_fraction <= 1.0) ||
            !(options.face_pairing_correction_trust_fraction > 0.0) ||
            !(options.face_pairing_energy_pair_tolerance > 0.0) ||
            !(options.face_pairing_energy_residual_fraction > 0.0) ||
            !(options.face_pairing_mass_relative_tolerance > 0.0) ||
            !(options.face_pairing_f_residual_growth_tolerance >= 0.0) ||
            options.beam_ledger_mode < BEAM_LEDGER_OFF ||
            options.beam_ledger_mode > BEAM_LEDGER_FULL)
            error = "invalid runtime option value";
        if (!beam_ledger_mode_explicit) {
            options.beam_ledger_mode = options.diagnostic_level >= 2
                ? BEAM_LEDGER_FULL : BEAM_LEDGER_SUMMARY;
        }
        if (options.beam_ledger_reference_enabled) {
            options.beam_ledger_mode = BEAM_LEDGER_FULL;
        }
        if (options.operator_audit_mode &&
            options.midpoint_acceleration_mode !=
                RUNTIME_MIDPOINT_ACCELERATION_NONE)
            error = "--operator-audit requires --midpoint-acceleration none";
        else if (options.operator_audit_mode &&
                 (options.restart_enabled || options.checkpoint_enabled ||
                  options.stop_after_accepted_steps >= 0))
            error = "--operator-audit is exclusive with time advancement options";
        else if (options.velocity_grid_remap_enabled &&
                 (options.restart_enabled || options.checkpoint_enabled ||
                  options.operator_audit_mode ||
                  options.stop_after_accepted_steps >= 0))
            error = "--remap-checkpoint-velocity-grid is exclusive with time advancement options";
    }
    broadcast_string(error, 0);
    if (!error.empty()) {
        if (mpi_rank == 0) std::fprintf(stderr, "Runtime option error: %s\n", error.c_str());
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    broadcast_string(options.output_dir, 0);
    broadcast_string(options.restart_dir, 0);
    broadcast_string(options.operator_audit_dir, 0);
    broadcast_string(options.beam_ledger_reference, 0);
    broadcast_string(options.remap_checkpoint_source, 0);
    broadcast_string(options.remap_checkpoint_destination, 0);
    double numbers[21] = {
        options.stop_time_fs,
        options.dt_scale,
        options.midpoint_trace_start_fs,
        options.midpoint_trace_end_fs,
        options.accepted_energy_audit_start_fs,
        options.accepted_energy_audit_end_fs,
        options.acceleration_accept_ratio,
        options.acceleration_max_coefficient,
        options.face_pairing_sigma_cutoff,
        options.face_pairing_lambda,
        options.face_pairing_eta,
        options.face_pairing_trust_fraction,
        options.face_pairing_correction_trust_fraction,
        options.face_pairing_energy_pair_tolerance,
        options.face_pairing_energy_residual_fraction,
        options.face_pairing_mass_relative_tolerance,
        options.face_pairing_f_residual_growth_tolerance,
        options.soft_candidate_field_tolerance,
        options.soft_candidate_current_tolerance,
        options.soft_candidate_energy_p99_reference,
        options.soft_candidate_energy_absolute_limit
    };
    long long accepted = options.stop_after_accepted_steps;
    int ints[15] = {options.midpoint_max_iters, options.diagnostic_level,
                   options.accepted_energy_audit_cadence,
                   options.dump_final_midpoint ? 1 : 0, options.overwrite_output ? 1 : 0,
                   options.checkpoint_enabled ? 1 : 0, options.restart_enabled ? 1 : 0,
                   options.background_coupling_mode, static_cast<int>(options.beam_ledger_mode),
                   static_cast<int>(options.midpoint_acceleration_mode),
                   options.anderson_depth, options.acceleration_start_iter,
                   options.face_pairing_mode,
                   static_cast<int>(options.midpoint_initial_guess_mode),
                   options.velocity_grid_remap_enabled ? 1 : 0};
    MPI_Bcast(numbers, 21, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&accepted, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ints, 15, MPI_INT, 0, MPI_COMM_WORLD);
    int audit = options.operator_audit_mode ? 1 : 0;
    MPI_Bcast(&audit, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int beam_ledger_reference = options.beam_ledger_reference_enabled ? 1 : 0;
    MPI_Bcast(&beam_ledger_reference, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int time_count = static_cast<int>(options.checkpoint_times_fs.size());
    MPI_Bcast(&time_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) options.checkpoint_times_fs.resize(static_cast<size_t>(time_count));
    if (time_count) MPI_Bcast(&options.checkpoint_times_fs[0], time_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    options.stop_time_fs = numbers[0]; options.dt_scale = numbers[1];
    options.midpoint_trace_start_fs = numbers[2];
    options.midpoint_trace_end_fs = numbers[3];
    options.accepted_energy_audit_start_fs = numbers[4];
    options.accepted_energy_audit_end_fs = numbers[5];
    options.acceleration_accept_ratio = numbers[6];
    options.acceleration_max_coefficient = numbers[7];
    options.face_pairing_sigma_cutoff = numbers[8];
    options.face_pairing_lambda = numbers[9];
    options.face_pairing_eta = numbers[10];
    options.face_pairing_trust_fraction = numbers[11];
    options.face_pairing_correction_trust_fraction = numbers[12];
    options.face_pairing_energy_pair_tolerance = numbers[13];
    options.face_pairing_energy_residual_fraction = numbers[14];
    options.face_pairing_mass_relative_tolerance = numbers[15];
    options.face_pairing_f_residual_growth_tolerance = numbers[16];
    options.soft_candidate_field_tolerance = numbers[17];
    options.soft_candidate_current_tolerance = numbers[18];
    options.soft_candidate_energy_p99_reference = numbers[19];
    options.soft_candidate_energy_absolute_limit = numbers[20];
    options.stop_after_accepted_steps = accepted; options.midpoint_max_iters = ints[0];
    options.diagnostic_level = ints[1];
    options.accepted_energy_audit_cadence = ints[2];
    options.dump_final_midpoint = ints[3] != 0;
    options.overwrite_output = ints[4] != 0; options.checkpoint_enabled = ints[5] != 0;
    options.restart_enabled = ints[6] != 0; options.operator_audit_mode = audit != 0;
    options.background_coupling_mode = ints[7];
    options.beam_ledger_mode = static_cast<BeamLedgerMode>(ints[8]);
    options.midpoint_acceleration_mode =
        static_cast<RuntimeMidpointAccelerationMode>(ints[9]);
    options.anderson_depth = ints[10];
    options.acceleration_start_iter = ints[11];
    options.face_pairing_mode = ints[12];
    options.midpoint_initial_guess_mode =
        static_cast<RuntimeMidpointInitialGuessMode>(ints[13]);
    options.velocity_grid_remap_enabled = ints[14] != 0;
    options.beam_ledger_reference_enabled = beam_ledger_reference != 0;
    (void)mpi_size;
    return options;
}

std::string output_path(const RuntimeOptions& options, const std::string& filename)
{
    return options.output_dir + "/" + filename;
}

bool prepare_output_directory(const RuntimeOptions& options, int mpi_rank,
                              int mpi_size, std::string& error)
{
    int ok = 1;
    if (mpi_rank == 0) {
#if defined(_WIN32)
        const int mkdir_status = _mkdir(options.output_dir.c_str());
#else
        const int mkdir_status = mkdir(options.output_dir.c_str(), 0777);
#endif
        if (mkdir_status != 0 && errno != EEXIST) {
            error = "cannot create output directory: " + options.output_dir;
            ok = 0;
        } else if (!options.overwrite_output && directory_has_entries(options.output_dir)) {
            error = "output directory is nonempty; pass --overwrite-output or choose --output-dir";
            ok = 0;
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    broadcast_string(error, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    (void)mpi_size;
    return ok != 0;
}
