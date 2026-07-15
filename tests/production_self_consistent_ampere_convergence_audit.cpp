#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const double kRoundoffFloor = 4096.0 * std::numeric_limits<double>::epsilon();

struct Result {
    int format_version;
    std::string test_run_id;
    int nv;
    int nmu;
    int nx;
    bool local_passes;
    double dt;
    double umax;
    double upar_stretch;
    double uperp_stretch;
    double min_upar_width;
    double min_uperp_width;
    double upar_width_at_drift;
    double initial_number_density;
    double initial_mean_vx;
    double initial_current_density;
    double initial_current_error;
    double initial_mean_vx_error;
    double cells_per_uth_upar;
    double cells_per_uth_uperp;
    int occupied_upar_cells;
    int occupied_uperp_cells;
    bool strict_dt;
    bool strict_dt_over_2;
    bool strict_dt_over_4;
    bool dt_state_advanced;
    bool dt_over_2_state_advanced;
    bool dt_over_4_state_advanced;
    bool dt_converged;
    bool dt_over_2_converged;
    bool dt_over_4_converged;
    bool dt_failed;
    bool dt_over_2_failed;
    bool dt_over_4_failed;
    double number_relative[3];
    double beam_current_linf[3];
    double limiter_active_fraction[3];
    double energy_relative;
    double current_pair_relative;
    double energy_scale;
    double r_fv;
    double r_couple;
    double direct_delta_k;
    double subtractive_delta_k;
    double delta_k_difference;
};

template <typename T>
bool parse_value(const std::map<std::string, std::string>& values,
                 const char* key, T& value)
{
    const std::map<std::string, std::string>::const_iterator found =
        values.find(key);
    if (found == values.end()) return false;
    std::istringstream input(found->second);
    input >> value;
    return static_cast<bool>(input);
}

bool read_result(const char* path, Result& result)
{
    std::ifstream in(path);
    if (!in) return false;
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(in, line)) {
        const std::string::size_type equal = line.find('=');
        if (equal != std::string::npos)
            values[line.substr(0, equal)] = line.substr(equal + 1);
    }

    int local_passes = 0;
    int strict_dt = 0;
    int strict_dt_over_2 = 0;
    int strict_dt_over_4 = 0;
    int dt_state_advanced = 0;
    int dt_over_2_state_advanced = 0;
    int dt_over_4_state_advanced = 0;
    int dt_converged = 0;
    int dt_over_2_converged = 0;
    int dt_over_4_converged = 0;
    int dt_failed = 0;
    int dt_over_2_failed = 0;
    int dt_over_4_failed = 0;
    const bool complete =
        parse_value(values, "format_version", result.format_version) &&
        parse_value(values, "test_run_id", result.test_run_id) &&
        parse_value(values, "velocity_nv", result.nv) &&
        parse_value(values, "velocity_nmu", result.nmu) &&
        parse_value(values, "spatial_nx", result.nx) &&
        parse_value(values, "local_passes", local_passes) &&
        parse_value(values, "dt", result.dt) &&
        parse_value(values, "velocity_umax", result.umax) &&
        parse_value(values, "velocity_upar_stretch", result.upar_stretch) &&
        parse_value(values, "velocity_uperp_stretch", result.uperp_stretch) &&
        parse_value(values, "min_upar_width", result.min_upar_width) &&
        parse_value(values, "min_uperp_width", result.min_uperp_width) &&
        parse_value(values, "upar_width_at_drift", result.upar_width_at_drift) &&
        parse_value(values, "initial_number_density", result.initial_number_density) &&
        parse_value(values, "initial_mean_vx", result.initial_mean_vx) &&
        parse_value(values, "initial_current_density", result.initial_current_density) &&
        parse_value(values, "initial_current_relative_error",
                    result.initial_current_error) &&
        parse_value(values, "initial_mean_vx_relative_error",
                    result.initial_mean_vx_error) &&
        parse_value(values, "cells_per_uth_upar_at_drift",
                    result.cells_per_uth_upar) &&
        parse_value(values, "cells_per_uth_uperp_at_axis",
                    result.cells_per_uth_uperp) &&
        parse_value(values, "occupied_upar_cells", result.occupied_upar_cells) &&
        parse_value(values, "occupied_uperp_cells", result.occupied_uperp_cells) &&
        parse_value(values, "dt_strict_accept", strict_dt) &&
        parse_value(values, "dt_over_2_strict_accept", strict_dt_over_2) &&
        parse_value(values, "dt_over_4_strict_accept", strict_dt_over_4) &&
        parse_value(values, "dt_state_advanced", dt_state_advanced) &&
        parse_value(values, "dt_over_2_state_advanced", dt_over_2_state_advanced) &&
        parse_value(values, "dt_over_4_state_advanced", dt_over_4_state_advanced) &&
        parse_value(values, "dt_converged", dt_converged) &&
        parse_value(values, "dt_over_2_converged", dt_over_2_converged) &&
        parse_value(values, "dt_over_4_converged", dt_over_4_converged) &&
        parse_value(values, "dt_failed", dt_failed) &&
        parse_value(values, "dt_over_2_failed", dt_over_2_failed) &&
        parse_value(values, "dt_over_4_failed", dt_over_4_failed) &&
        parse_value(values, "dt_number_relative", result.number_relative[0]) &&
        parse_value(values, "dt_over_2_number_relative", result.number_relative[1]) &&
        parse_value(values, "dt_over_4_number_relative", result.number_relative[2]) &&
        parse_value(values, "dt_beam_current_linf", result.beam_current_linf[0]) &&
        parse_value(values, "dt_over_2_beam_current_linf", result.beam_current_linf[1]) &&
        parse_value(values, "dt_over_4_beam_current_linf", result.beam_current_linf[2]) &&
        parse_value(values, "dt_limiter_active_fraction",
                    result.limiter_active_fraction[0]) &&
        parse_value(values, "dt_over_2_limiter_active_fraction",
                    result.limiter_active_fraction[1]) &&
        parse_value(values, "dt_over_4_limiter_active_fraction",
                    result.limiter_active_fraction[2]) &&
        parse_value(values, "dt_over_4_energy_relative", result.energy_relative) &&
        parse_value(values, "dt_over_4_j_pair_relative",
                    result.current_pair_relative) &&
        parse_value(values, "dt_over_4_energy_exchange_scale", result.energy_scale) &&
        parse_value(values, "dt_over_4_stage5_R_FV", result.r_fv) &&
        parse_value(values, "dt_over_4_stage5_R_couple", result.r_couple) &&
        parse_value(values, "direct_delta_K_bkg", result.direct_delta_k) &&
        parse_value(values, "subtractive_delta_K_bkg",
                    result.subtractive_delta_k) &&
        parse_value(values, "delta_K_method_difference", result.delta_k_difference);
    result.local_passes = local_passes != 0;
    result.strict_dt = strict_dt != 0;
    result.strict_dt_over_2 = strict_dt_over_2 != 0;
    result.strict_dt_over_4 = strict_dt_over_4 != 0;
    result.dt_state_advanced = dt_state_advanced != 0;
    result.dt_over_2_state_advanced = dt_over_2_state_advanced != 0;
    result.dt_over_4_state_advanced = dt_over_4_state_advanced != 0;
    result.dt_converged = dt_converged != 0;
    result.dt_over_2_converged = dt_over_2_converged != 0;
    result.dt_over_4_converged = dt_over_4_converged != 0;
    result.dt_failed = dt_failed != 0;
    result.dt_over_2_failed = dt_over_2_failed != 0;
    result.dt_over_4_failed = dt_over_4_failed != 0;
    return complete && result.format_version == 2 && !result.test_run_id.empty() &&
        std::isfinite(result.dt) && std::isfinite(result.umax) &&
        std::isfinite(result.upar_stretch) &&
        std::isfinite(result.uperp_stretch) &&
        std::isfinite(result.min_upar_width) &&
        std::isfinite(result.min_uperp_width) &&
        std::isfinite(result.upar_width_at_drift) &&
        std::isfinite(result.initial_number_density) &&
        std::isfinite(result.initial_mean_vx) &&
        std::isfinite(result.initial_current_density) &&
        std::isfinite(result.initial_current_error) &&
        std::isfinite(result.initial_mean_vx_error) &&
        std::isfinite(result.energy_relative) &&
        std::isfinite(result.current_pair_relative) &&
        std::isfinite(result.r_fv) && std::isfinite(result.r_couple) &&
        std::isfinite(result.number_relative[0]) &&
        std::isfinite(result.number_relative[1]) &&
        std::isfinite(result.number_relative[2]) &&
        std::isfinite(result.beam_current_linf[0]) &&
        std::isfinite(result.beam_current_linf[1]) &&
        std::isfinite(result.beam_current_linf[2]) &&
        std::isfinite(result.limiter_active_fraction[0]) &&
        std::isfinite(result.limiter_active_fraction[1]) &&
        std::isfinite(result.limiter_active_fraction[2]) &&
        std::isfinite(result.direct_delta_k) &&
        std::isfinite(result.subtractive_delta_k) &&
        std::isfinite(result.delta_k_difference);
}

bool same_mapping_value(double lhs, double rhs)
{
    return std::fabs(lhs - rhs) <= 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(lhs), std::fabs(rhs)));
}

bool same_positive_scale(double lhs, double rhs)
{
    return lhs > 0.0 && rhs > 0.0 &&
        std::fabs(lhs - rhs) <= 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(std::fabs(lhs), std::fabs(rhs));
}

bool width_halves(double coarse, double fine)
{
    const double ratio = fine / coarse;
    return ratio >= 0.45 && ratio <= 0.55;
}

double refinement_order(double coarse, double fine)
{
    if (coarse <= kRoundoffFloor && fine <= kRoundoffFloor)
        return std::numeric_limits<double>::infinity();
    if (!(coarse > 0.0) || !(fine > 0.0))
        return -std::numeric_limits<double>::infinity();
    return std::log(coarse / fine) / std::log(2.0);
}

bool decreases_without_rebound(double coarse, double base, double fine)
{
    if (coarse <= kRoundoffFloor && base <= kRoundoffFloor &&
        fine <= kRoundoffFloor) return true;
    return base < coarse && fine < base;
}

double normalized_residual(double residual, double scale)
{
    return std::fabs(residual) /
        std::max(std::numeric_limits<double>::min(), std::fabs(scale));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: " << argv[0]
                  << " coarse_result base_result fine_result\n";
        return 2;
    }

    std::vector<Result> results(3);
    for (int i = 0; i < 3; ++i) {
        if (!read_result(argv[i + 1], results[static_cast<size_t>(i)])) {
            std::cerr << "cannot read valid version-2 7.2 result: "
                      << argv[i + 1] << "\n";
            return 2;
        }
    }
    std::sort(results.begin(), results.end(),
        [](const Result& lhs, const Result& rhs) { return lhs.nv < rhs.nv; });

    const bool grids_are_nested =
        results[1].nv == 2 * results[0].nv &&
        results[2].nv == 2 * results[1].nv &&
        results[1].nmu == 2 * results[0].nmu &&
        results[2].nmu == 2 * results[1].nmu;
    const bool mapping_identical =
        same_mapping_value(results[0].umax, results[1].umax) &&
        same_mapping_value(results[1].umax, results[2].umax) &&
        same_mapping_value(results[0].upar_stretch, results[1].upar_stretch) &&
        same_mapping_value(results[1].upar_stretch, results[2].upar_stretch) &&
        same_mapping_value(results[0].uperp_stretch, results[1].uperp_stretch) &&
        same_mapping_value(results[1].uperp_stretch, results[2].uperp_stretch);
    const bool common_physical_dt =
        same_positive_scale(results[0].dt, results[1].dt) &&
        same_positive_scale(results[1].dt, results[2].dt) &&
        results[0].nx == results[1].nx && results[1].nx == results[2].nx;
    const bool test_batch_consistent =
        results[0].test_run_id == results[1].test_run_id &&
        results[1].test_run_id == results[2].test_run_id;
    const bool widths_converge =
        width_halves(results[0].min_upar_width, results[1].min_upar_width) &&
        width_halves(results[1].min_upar_width, results[2].min_upar_width) &&
        width_halves(results[0].min_uperp_width, results[1].min_uperp_width) &&
        width_halves(results[1].min_uperp_width, results[2].min_uperp_width);
    const bool initialization_converges = decreases_without_rebound(
        results[0].initial_current_error, results[1].initial_current_error,
        results[2].initial_current_error) && decreases_without_rebound(
        results[0].initial_mean_vx_error, results[1].initial_mean_vx_error,
        results[2].initial_mean_vx_error);
    const bool finest_initialization_resolved =
        results[2].initial_current_error < 2.0e-2 &&
        results[2].cells_per_uth_upar >= 2.0 &&
        results[2].cells_per_uth_uperp >= 2.0;
    const bool energy_grid_converges = decreases_without_rebound(
        results[0].energy_relative, results[1].energy_relative,
        results[2].energy_relative);
    const bool pair_grid_converges = decreases_without_rebound(
        results[0].current_pair_relative, results[1].current_pair_relative,
        results[2].current_pair_relative);

    double r_couple_relative[3];
    bool r_fv_roundoff = true;
    for (int i = 0; i < 3; ++i) {
        r_couple_relative[i] = normalized_residual(
            results[i].r_couple, results[i].energy_scale);
        const double r_fv_relative = normalized_residual(
            results[i].r_fv, results[i].energy_scale);
        const double accumulated_roundoff = 64.0 *
            std::numeric_limits<double>::epsilon() * results[i].nx *
            results[i].nv * results[i].nmu;
        r_fv_roundoff = r_fv_roundoff &&
            r_fv_relative <= accumulated_roundoff;
    }
    const bool r_couple_converges = decreases_without_rebound(
        r_couple_relative[0], r_couple_relative[1], r_couple_relative[2]);
    const bool strict_cases_accept = results[0].strict_dt &&
        results[0].strict_dt_over_2 && results[0].strict_dt_over_4 &&
        results[1].strict_dt && results[1].strict_dt_over_2 &&
        results[1].strict_dt_over_4 && results[2].strict_dt &&
        results[2].strict_dt_over_2 && results[2].strict_dt_over_4;
    const bool local_passes = results[0].local_passes &&
        results[1].local_passes && results[2].local_passes;
    const bool passes = grids_are_nested && mapping_identical && common_physical_dt &&
        test_batch_consistent && widths_converge &&
        local_passes && initialization_converges && finest_initialization_resolved &&
        energy_grid_converges && pair_grid_converges && r_fv_roundoff &&
        r_couple_converges && strict_cases_accept;

    const char* failure_class = "PASS";
    if (!test_batch_consistent)
        failure_class = "MIXED_TEST_BATCH";
    else if (!grids_are_nested || !mapping_identical || !common_physical_dt ||
             !widths_converge)
        failure_class = "INVALID_GRID_SUITE";
    else if (!strict_cases_accept || !local_passes)
        failure_class = "CASE_NOT_STRICTLY_ACCEPTED";
    else if (!initialization_converges || !finest_initialization_resolved)
        failure_class = "GRID_NOT_RESOLVED";
    else if (!energy_grid_converges || !pair_grid_converges ||
             !r_fv_roundoff || !r_couple_converges)
        failure_class = "LOW_ORDER_COUPLING_NOT_CONVERGED";

    std::cout << std::scientific << std::setprecision(16);
    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << "velocity_grid Nv=" << results[i].nv
                  << " Nmu=" << results[i].nmu
                  << " test_run_id=" << results[i].test_run_id
                  << " min_upar_width=" << results[i].min_upar_width
                  << " min_uperp_width=" << results[i].min_uperp_width
                  << " initial_number_density="
                  << results[i].initial_number_density
                  << " initial_mean_vx=" << results[i].initial_mean_vx
                  << " initial_current_density="
                  << results[i].initial_current_density
                  << " upar_width_at_drift="
                  << results[i].upar_width_at_drift
                  << " initial_current_relative_error="
                  << results[i].initial_current_error
                  << " initial_mean_vx_relative_error="
                  << results[i].initial_mean_vx_error
                  << " energy_relative_dt_over_4=" << results[i].energy_relative
                  << " JN_minus_GstarJE_relative_dt_over_4="
                  << results[i].current_pair_relative
                  << " R_FV_relative_dt_over_4="
                  << normalized_residual(results[i].r_fv, results[i].energy_scale)
                  << " R_couple_relative_dt_over_4=" << r_couple_relative[i]
                  << " delta_K_method_difference=" << results[i].delta_k_difference
                  << " strict_dt=" << results[i].strict_dt
                  << " strict_dt_over_2=" << results[i].strict_dt_over_2
                  << " strict_dt_over_4=" << results[i].strict_dt_over_4
                  << " dt_number_relative=" << results[i].number_relative[0]
                  << " dt_over_2_number_relative=" << results[i].number_relative[1]
                  << " dt_over_4_number_relative=" << results[i].number_relative[2]
                  << " dt_beam_current_linf=" << results[i].beam_current_linf[0]
                  << " dt_over_2_beam_current_linf=" << results[i].beam_current_linf[1]
                  << " dt_over_4_beam_current_linf=" << results[i].beam_current_linf[2]
                  << " dt_limiter_active_fraction="
                  << results[i].limiter_active_fraction[0]
                  << " dt_over_2_limiter_active_fraction="
                  << results[i].limiter_active_fraction[1]
                  << " dt_over_4_limiter_active_fraction="
                  << results[i].limiter_active_fraction[2]
                  << " local_passes=" << results[i].local_passes << "\n";
    }
    std::cout << "initial_current_order_coarse_to_base=" << refinement_order(
        results[0].initial_current_error, results[1].initial_current_error)
              << " initial_current_order_base_to_fine=" << refinement_order(
        results[1].initial_current_error, results[2].initial_current_error)
              << " energy_order_coarse_to_base=" << refinement_order(
        results[0].energy_relative, results[1].energy_relative)
              << " energy_order_base_to_fine=" << refinement_order(
        results[1].energy_relative, results[2].energy_relative)
              << " current_pair_order_coarse_to_base=" << refinement_order(
        results[0].current_pair_relative, results[1].current_pair_relative)
              << " current_pair_order_base_to_fine=" << refinement_order(
        results[1].current_pair_relative, results[2].current_pair_relative)
              << " mapping_identical=" << mapping_identical
              << " common_physical_dt=" << common_physical_dt
              << " test_batch_consistent=" << test_batch_consistent
              << " widths_converge=" << widths_converge
              << " initialization_converges=" << initialization_converges
              << " finest_initialization_resolved=" << finest_initialization_resolved
              << " energy_grid_converges=" << energy_grid_converges
              << " current_pair_grid_converges=" << pair_grid_converges
              << " R_FV_roundoff=" << r_fv_roundoff
              << " R_couple_converges=" << r_couple_converges
              << " strict_cases_accept=" << strict_cases_accept
              << " failure_class=" << failure_class
              << " result=" << (passes ? "PASS" : "FAIL") << "\n";
    return passes ? 0 : 1;
}
