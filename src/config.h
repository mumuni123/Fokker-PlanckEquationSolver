#ifndef CONFIG_H
#define CONFIG_H

namespace Config {
    // Set these values here and rebuild to control diagnostic modules.
    const bool enable_debug_diagnostics = false;
    const bool enable_step_diagnostics = true;
    const bool enable_full_fe_output = false;
    const bool enable_progress_trace = false;

    const int progress_trace_initial_steps = 3;
    const int progress_trace_interval = 100;
}

struct RuntimeConfig {
    bool enable_debug_diagnostics;
    bool enable_step_diagnostics;
    bool enable_full_fe_output;
    bool enable_progress_trace;
};

RuntimeConfig load_runtime_config();
bool should_trace_progress(const RuntimeConfig& config, int step);
void trace_progress(const RuntimeConfig& config, int mpi_rank, int step,
                    const char* stage);

#endif
