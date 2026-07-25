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
    o.midpoint_acceleration_mode = RUNTIME_MIDPOINT_ACCELERATION_NONE;
    o.anderson_depth = 3;
    o.acceleration_start_iter = 3;
    o.acceleration_accept_ratio = 0.95;
    o.acceleration_max_coefficient = 2.0;
    o.beam_ledger_mode = BEAM_LEDGER_SUMMARY;
    o.background_coupling_mode = 1;
    o.dump_final_midpoint = false;
    o.overwrite_output = false;
    o.checkpoint_enabled = false;
    o.restart_enabled = false;
    o.operator_audit_mode = false;
    o.beam_ledger_reference_enabled = false;
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
            } else if (arg == "--background-coupling-mode" && has_value) {
                const std::string mode(value);
                if (mode == "legacy") options.background_coupling_mode = 0;
                else if (mode == "dual-u" || mode == "dual_u")
                    options.background_coupling_mode = 1;
                else error = "invalid --background-coupling-mode (use legacy or dual-u)";
                ++i;
            } else if (arg == "--output-dir" && has_value) {
                options.output_dir = value; ++i;
            } else error = "unknown or incomplete runtime option: " + arg;
        }
        if (options.dt_scale <= 0.0 || options.midpoint_max_iters < 1 ||
            options.diagnostic_level < 0 || options.diagnostic_level > 2 ||
            options.midpoint_acceleration_mode < RUNTIME_MIDPOINT_ACCELERATION_NONE ||
            options.midpoint_acceleration_mode > RUNTIME_MIDPOINT_ACCELERATION_ANDERSON ||
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
    double numbers[6] = {
        options.stop_time_fs,
        options.dt_scale,
        options.midpoint_trace_start_fs,
        options.midpoint_trace_end_fs,
        options.acceleration_accept_ratio,
        options.acceleration_max_coefficient
    };
    long long accepted = options.stop_after_accepted_steps;
    int ints[11] = {options.midpoint_max_iters, options.diagnostic_level,
                   options.dump_final_midpoint ? 1 : 0, options.overwrite_output ? 1 : 0,
                   options.checkpoint_enabled ? 1 : 0, options.restart_enabled ? 1 : 0,
                   options.background_coupling_mode, static_cast<int>(options.beam_ledger_mode),
                   static_cast<int>(options.midpoint_acceleration_mode),
                   options.anderson_depth, options.acceleration_start_iter};
    MPI_Bcast(numbers, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&accepted, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ints, 11, MPI_INT, 0, MPI_COMM_WORLD);
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
    options.acceleration_accept_ratio = numbers[4];
    options.acceleration_max_coefficient = numbers[5];
    options.stop_after_accepted_steps = accepted; options.midpoint_max_iters = ints[0];
    options.diagnostic_level = ints[1]; options.dump_final_midpoint = ints[2] != 0;
    options.overwrite_output = ints[3] != 0; options.checkpoint_enabled = ints[4] != 0;
    options.restart_enabled = ints[5] != 0; options.operator_audit_mode = audit != 0;
    options.background_coupling_mode = ints[6];
    options.beam_ledger_mode = static_cast<BeamLedgerMode>(ints[7]);
    options.midpoint_acceleration_mode =
        static_cast<RuntimeMidpointAccelerationMode>(ints[8]);
    options.anderson_depth = ints[9];
    options.acceleration_start_iter = ints[10];
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
        if (mkdir(options.output_dir.c_str(), 0777) != 0 && errno != EEXIST) {
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
