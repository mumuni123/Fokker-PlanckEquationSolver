#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <chrono>
#include <vector>

namespace {

std::map<std::string, double> read_result(const char* path)
{
    std::ifstream in(path);
    std::map<std::string, double> values;
    std::string line;
    while (std::getline(in, line)) {
        const std::string::size_type split = line.find('=');
        if (split == std::string::npos) continue;
        values[line.substr(0, split)] = std::strtod(
            line.substr(split + 1).c_str(), 0);
    }
    return values;
}

bool contains_required_keys(const std::map<std::string, double>& values)
{
    static const char* required[] = {
        "operator_valid", "outputs_finite", "nx", "Nu", "Nuperp",
        "wN_minus_wE_L2", "pair_work_L2", "tail_weight_relative_L1",
        "frozen_weight_fraction", "active_fully_covered_fraction",
        "gR_work_relative", "kink_active_mass_fraction",
        "kink_pair_work_fraction"
    };
    for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i)
        if (values.find(required[i]) == values.end()) return false;
    return true;
}

std::map<std::string, double> read_complete_result(const char* path)
{
    // Producers publish their compact files by atomic rename.  The retry is
    // still useful when an external launcher starts this reader concurrently.
    for (int attempt = 0; attempt < 240; ++attempt) {
        const std::map<std::string, double> values = read_result(path);
        if (contains_required_keys(values)) return values;
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
    return read_result(path);
}

double value(const std::map<std::string, double>& result, const char* key)
{
    const std::map<std::string, double>::const_iterator it = result.find(key);
    return it == result.end() ? std::numeric_limits<double>::quiet_NaN()
                              : it->second;
}

double observed_order(double coarse, double fine)
{
    return (std::isfinite(coarse) && std::isfinite(fine) && coarse > 0.0 &&
            fine > 0.0) ? std::log(coarse / fine) / std::log(2.0)
                         : std::numeric_limits<double>::quiet_NaN();
}

bool flag(const std::map<std::string, double>& result,
          const char* preferred, const char* fallback)
{
    const std::map<std::string, double>::const_iterator preferred_it =
        result.find(preferred);
    if (preferred_it != result.end()) return preferred_it->second == 1.0;
    return value(result, fallback) == 1.0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: " << argv[0]
                  << " coarse.result medium.result fine.result\n";
        return 2;
    }
    const std::map<std::string, double> coarse = read_complete_result(argv[1]);
    const std::map<std::string, double> medium = read_complete_result(argv[2]);
    const std::map<std::string, double> fine = read_complete_result(argv[3]);
    const double w[3] = {value(coarse, "wN_minus_wE_L2"),
                         value(medium, "wN_minus_wE_L2"),
                         value(fine, "wN_minus_wE_L2")};
    const double pair[3] = {value(coarse, "pair_work_L2"),
                            value(medium, "pair_work_L2"),
                            value(fine, "pair_work_L2")};
    const double tail[3] = {value(coarse, "tail_weight_relative_L1"),
                            value(medium, "tail_weight_relative_L1"),
                            value(fine, "tail_weight_relative_L1")};
    const double coverage[3] = {value(coarse, "frozen_weight_fraction"),
                                value(medium, "frozen_weight_fraction"),
                                value(fine, "frozen_weight_fraction")};
    const double active_coverage[3] = {
        value(coarse, "active_fully_covered_fraction"),
        value(medium, "active_fully_covered_fraction"),
        value(fine, "active_fully_covered_fraction")};
    const double work_relative[3] = {value(coarse, "gR_work_relative"),
                                     value(medium, "gR_work_relative"),
                                     value(fine, "gR_work_relative")};
    const double kink_mass[3] = {value(coarse, "kink_active_mass_fraction"),
                                 value(medium, "kink_active_mass_fraction"),
                                 value(fine, "kink_active_mass_fraction")};
    const double kink_work[3] = {value(coarse, "kink_pair_work_fraction"),
                                 value(medium, "kink_pair_work_fraction"),
                                 value(fine, "kink_pair_work_fraction")};
    const double active_coverage_threshold = 0.999;
    const bool direct_inputs_valid =
        flag(coarse, "direct_pair_sample_valid", "operator_valid") &&
        flag(medium, "direct_pair_sample_valid", "operator_valid") &&
        flag(fine, "direct_pair_sample_valid", "operator_valid") &&
        value(coarse, "outputs_finite") == 1.0 &&
        value(medium, "outputs_finite") == 1.0 &&
        value(fine, "outputs_finite") == 1.0 &&
        value(coarse, "Nuperp") == value(medium, "Nuperp") &&
        value(medium, "Nuperp") == value(fine, "Nuperp") &&
        value(coarse, "nx") == value(medium, "nx") &&
        value(medium, "nx") == value(fine, "nx");
    const bool jacobian_metrics_finite[3] = {
        std::isfinite(w[0]) && std::isfinite(work_relative[0]) &&
            std::isfinite(kink_mass[0]) && std::isfinite(kink_work[0]),
        std::isfinite(w[1]) && std::isfinite(work_relative[1]) &&
            std::isfinite(kink_mass[1]) && std::isfinite(kink_work[1]),
        std::isfinite(w[2]) && std::isfinite(work_relative[2]) &&
            std::isfinite(kink_mass[2]) && std::isfinite(kink_work[2])
    };
    const bool jacobian_level_valid[3] = {
        flag(coarse, "jacobian_transpose_valid", "frozen_jacobian_valid") &&
            jacobian_metrics_finite[0] &&
            active_coverage[0] >= active_coverage_threshold,
        flag(medium, "jacobian_transpose_valid", "frozen_jacobian_valid") &&
            jacobian_metrics_finite[1] &&
            active_coverage[1] >= active_coverage_threshold,
        flag(fine, "jacobian_transpose_valid", "frozen_jacobian_valid") &&
            jacobian_metrics_finite[2] &&
            active_coverage[2] >= active_coverage_threshold
    };
    const bool jacobian_inputs_valid = jacobian_level_valid[0] &&
        jacobian_level_valid[1] && jacobian_level_valid[2];
    const double w_order_cm = observed_order(w[0], w[1]);
    const double w_order_mf = observed_order(w[1], w[2]);
    const double pair_order_cm = observed_order(pair[0], pair[1]);
    const double pair_order_mf = observed_order(pair[1], pair[2]);
    const bool direct_norms_finite = std::isfinite(pair[0]) &&
        std::isfinite(pair[1]) && std::isfinite(pair[2]);
    const bool jacobian_norms_finite = jacobian_metrics_finite[0] &&
        jacobian_metrics_finite[1] && jacobian_metrics_finite[2];
    const bool direct_decreases = direct_norms_finite &&
        pair[2] < pair[1] && pair[1] < pair[0];
    const bool truncation_dominated = direct_inputs_valid &&
        direct_norms_finite &&
        pair_order_cm >= 1.5 && pair_order_mf >= 1.5 &&
        direct_decreases;
    const bool structural_plateau = direct_inputs_valid &&
        direct_norms_finite &&
        !truncation_dominated &&
        pair[2] / std::max(pair[0], 1.0e-300) > 0.5;
    const bool nonsmooth_limiter_dominated = jacobian_norms_finite &&
        (std::max(kink_mass[0], std::max(kink_mass[1], kink_mass[2])) > 1.0e-2 ||
         std::max(kink_work[0], std::max(kink_work[1], kink_work[2])) > 5.0e-2);
    const bool direct_pair_audit_inconclusive = !direct_inputs_valid ||
        !direct_norms_finite ||
        (!truncation_dominated && !structural_plateau);
    const bool direct_pair_convergence_valid = direct_inputs_valid &&
        direct_norms_finite && truncation_dominated;
    const bool jacobian_transpose_valid = jacobian_inputs_valid;
    const bool jacobian_transpose_incomplete = !jacobian_transpose_valid;

    std::ofstream out("output/background_coupling_upar_transpose_convergence_audit.result");
    std::ostream& log = out ? out : std::cout;
    log << std::scientific << std::setprecision(17)
        << "test=background_coupling_upar_transpose_convergence_audit\n"
        << "input_coarse=" << argv[1] << "\ninput_medium=" << argv[2]
        << "\ninput_fine=" << argv[3] << "\n"
        << "wN_minus_wE_L2_coarse=" << w[0] << "\n"
        << "wN_minus_wE_L2_medium=" << w[1] << "\n"
        << "wN_minus_wE_L2_fine=" << w[2] << "\n"
        << "pair_work_L2_coarse=" << pair[0] << "\n"
        << "pair_work_L2_medium=" << pair[1] << "\n"
        << "pair_work_L2_fine=" << pair[2] << "\n"
        << "tail_weight_relative_L1_coarse=" << tail[0] << "\n"
        << "tail_weight_relative_L1_medium=" << tail[1] << "\n"
        << "tail_weight_relative_L1_fine=" << tail[2] << "\n"
        << "frozen_weight_fraction_coarse=" << coverage[0] << "\n"
        << "frozen_weight_fraction_medium=" << coverage[1] << "\n"
        << "frozen_weight_fraction_fine=" << coverage[2] << "\n"
        << "active_fully_covered_fraction_coarse=" << active_coverage[0] << "\n"
        << "active_fully_covered_fraction_medium=" << active_coverage[1] << "\n"
        << "active_fully_covered_fraction_fine=" << active_coverage[2] << "\n"
        << "active_coverage_threshold=" << active_coverage_threshold << "\n"
        << "jacobian_level_valid_coarse=" << jacobian_level_valid[0] << "\n"
        << "jacobian_level_valid_medium=" << jacobian_level_valid[1] << "\n"
        << "jacobian_level_valid_fine=" << jacobian_level_valid[2] << "\n"
        << "gR_work_relative_coarse=" << work_relative[0] << "\n"
        << "gR_work_relative_medium=" << work_relative[1] << "\n"
        << "gR_work_relative_fine=" << work_relative[2] << "\n"
        << "kink_active_mass_fraction_coarse=" << kink_mass[0] << "\n"
        << "kink_active_mass_fraction_medium=" << kink_mass[1] << "\n"
        << "kink_active_mass_fraction_fine=" << kink_mass[2] << "\n"
        << "kink_pair_work_fraction_coarse=" << kink_work[0] << "\n"
        << "kink_pair_work_fraction_medium=" << kink_work[1] << "\n"
        << "kink_pair_work_fraction_fine=" << kink_work[2] << "\n"
        << "wN_minus_wE_order_coarse_to_medium=" << w_order_cm << "\n"
        << "wN_minus_wE_order_medium_to_fine=" << w_order_mf << "\n"
        << "pair_work_order_coarse_to_medium=" << pair_order_cm << "\n"
        << "pair_work_order_medium_to_fine=" << pair_order_mf << "\n"
        << "truncation_trend_detected=" << direct_decreases << "\n"
        << "truncation_dominated=" << truncation_dominated << "\n"
        << "structural_plateau=" << structural_plateau << "\n"
        << "nonsmooth_limiter_dominated=" << nonsmooth_limiter_dominated << "\n"
        << "direct_pair_audit_inconclusive="
        << direct_pair_audit_inconclusive << "\n"
        << "jacobian_transpose_incomplete="
        << jacobian_transpose_incomplete << "\n"
        << "direct_inputs_valid=" << direct_inputs_valid << "\n"
        << "direct_norms_finite=" << direct_norms_finite << "\n"
        << "jacobian_inputs_valid=" << jacobian_inputs_valid << "\n"
        << "jacobian_norms_finite=" << jacobian_norms_finite << "\n"
        << "direct_pair_convergence_valid="
        << direct_pair_convergence_valid << "\n"
        << "jacobian_transpose_valid=" << jacobian_transpose_valid << "\n";
    std::cout << "background_coupling_upar_transpose_convergence_audit result="
              << "output/background_coupling_upar_transpose_convergence_audit.result"
              << " direct_pair_convergence_valid="
              << direct_pair_convergence_valid
              << " jacobian_transpose_valid=" << jacobian_transpose_valid
              << "\n";
    return direct_pair_convergence_valid ? 0 : 1;
}
