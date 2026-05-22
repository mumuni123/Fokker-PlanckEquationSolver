#include "config.h"
#include "parameters.h"

#include <cstdio>

RuntimeConfig load_runtime_config()
{
    RuntimeConfig config;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    config.enable_debug_diagnostics = Config::enable_debug_diagnostics;
#else
    config.enable_debug_diagnostics = false;
#endif
    config.enable_full_fe_output = Config::enable_full_fe_output;
    config.enable_progress_trace = Config::enable_progress_trace;
    return config;
}

bool should_trace_progress(const RuntimeConfig& config, int step)
{
    if (!config.enable_progress_trace) return false;
    if (step <= Config::progress_trace_initial_steps) return true;
    return Config::progress_trace_interval > 0 &&
           step % Config::progress_trace_interval == 0;
}

void trace_progress(const RuntimeConfig& config, int mpi_rank, int step,
                    const char* stage)
{
    if (should_trace_progress(config, step) && mpi_rank == 0) {
        printf("  [progress] step %d: %s\n", step, stage);
        fflush(stdout);
    }
}
