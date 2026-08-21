#include "vpfp_time_control.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::string result_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result" && i + 1 < argc) result_path = argv[++i];
    }

    const double stop = 1.02858508730862e-13;
    const double nominal_dt = 1.2792543654303769e-17;
    const double tolerance =
        VpfpTimeControl::roundoff_tolerance(stop, stop, nominal_dt);
    const double natural_endpoint = 1.0285850873086075e-13;
    const double remaining = stop - natural_endpoint;

    const bool tolerance_has_physical_scale =
        tolerance > 0.0 && tolerance < 1.0e-25;
    const bool suppresses_roundoff_tail_step =
        !VpfpTimeControl::should_advance(natural_endpoint, stop, nominal_dt);
    const bool preserves_full_step = VpfpTimeControl::should_advance(
        stop - nominal_dt, stop, nominal_dt);
    const bool reaches_rounded_target = VpfpTimeControl::target_reached(
        natural_endpoint, stop, nominal_dt);
    const bool does_not_reach_early_target = !VpfpTimeControl::target_reached(
        stop - nominal_dt, stop, nominal_dt);
    const bool pass = tolerance_has_physical_scale &&
                      remaining <= tolerance &&
                      suppresses_roundoff_tail_step &&
                      preserves_full_step &&
                      reaches_rounded_target &&
                      does_not_reach_early_target;

    std::ostream* output = &std::cout;
    std::ofstream file;
    if (!result_path.empty()) {
        file.open(result_path.c_str());
        if (!file) return 2;
        output = &file;
    }
    *output << "tolerance=" << tolerance << "\n"
            << "remaining=" << remaining << "\n"
            << "tolerance_has_physical_scale="
            << (tolerance_has_physical_scale ? 1 : 0) << "\n"
            << "suppresses_roundoff_tail_step="
            << (suppresses_roundoff_tail_step ? 1 : 0) << "\n"
            << "preserves_full_step=" << (preserves_full_step ? 1 : 0) << "\n"
            << "reaches_rounded_target="
            << (reaches_rounded_target ? 1 : 0) << "\n"
            << "does_not_reach_early_target="
            << (does_not_reach_early_target ? 1 : 0) << "\n"
            << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
