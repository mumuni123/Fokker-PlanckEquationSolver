#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

namespace {

bool read_values(const char* path, std::map<std::string, double>& values,
                 std::string& revision)
{
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        const size_t equal = line.find('=');
        if (equal == std::string::npos) continue;
        if (line.substr(0, equal) == "git_revision")
            revision = line.substr(equal + 1);
        char* end = 0;
        const double value = std::strtod(line.c_str() + equal + 1, &end);
        if (end != line.c_str() + equal + 1)
            values[line.substr(0, equal)] = value;
    }
    const char* required[] = {"passes", "nx", "Nupar", "Nuperp",
        "velocity_grid_nonuniform", "JN_minus_GstarJE_L2",
        "legacy_JN_minus_GstarJE_L2", "mass_error_relative",
        "momentum_error_relative", "dual_cell_JE_minus_target_relative_L2",
        "projection_floor_L2", "projection_reconstruction_relative_L2",
        "uperp_grid_hash", "dual_decomposition_pass"};
    for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i)
        if (!values.count(required[i])) return false;
    return !revision.empty() && values.count("git_dirty");
}

std::string output_path(int argc, char** argv)
{
    for (int i = 4; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--output") return argv[i + 1];
    return "output/dual_u_stage2/manufactured_convergence/summary.result";
}

double order(double coarse, double fine, double nc, double nf)
{
    if (!(coarse > 0.0) || !(fine > 0.0) || !(nf > nc)) return 0.0;
    return std::log(coarse / fine) / std::log(nf / nc);
}

double relative_difference(double a, double b)
{
    const double scale = std::max(std::fabs(a), std::fabs(b));
    if (scale == 0.0) return 0.0;
    return std::fabs(a - b) / scale;
}

int parallel_process_count_hint()
{
    const char* variables[] = {"PMI_SIZE", "PMIX_SIZE",
        "OMPI_COMM_WORLD_SIZE", "MV2_COMM_WORLD_SIZE"};
    int maximum = 1;
    for (size_t i = 0; i < sizeof(variables) / sizeof(variables[0]); ++i) {
        const char* value = std::getenv(variables[i]);
        if (!value || !*value) continue;
        char* end = 0;
        const long parsed = std::strtol(value, &end, 10);
        if (end != value && parsed > maximum)
            maximum = static_cast<int>(parsed);
    }
    return maximum;
}

} // namespace

int main(int argc, char** argv)
{
    const int process_count_hint = parallel_process_count_hint();
    if (process_count_hint > 1) {
        std::cerr << "manufactured convergence aggregator must run as one "
                  << "process; detected launcher size=" << process_count_hint
                  << "\n";
        return 3;
    }
    if (argc < 4) {
        std::cerr << "usage: background_coupling_dual_u_convergence_audit "
                  << "coarse.result base.result fine.result "
                  << "[--output summary.result]\n";
        return 2;
    }
    std::map<std::string, double> coarse, base, fine;
    std::string coarse_revision, base_revision, fine_revision;
    if (!read_values(argv[1], coarse, coarse_revision) ||
        !read_values(argv[2], base, base_revision) ||
        !read_values(argv[3], fine, fine_revision)) {
        std::cerr << "unable to read stage-2 manufactured results\n";
        return 2;
    }
    const bool resolution_valid = coarse["Nupar"] == 48.0 &&
        base["Nupar"] == 96.0 && fine["Nupar"] == 192.0 &&
        coarse["nx"] == base["nx"] && base["nx"] == fine["nx"] &&
        coarse["Nuperp"] == base["Nuperp"] &&
        base["Nuperp"] == fine["Nuperp"];
    const bool inputs_pass = coarse["passes"] == 1.0 &&
        base["passes"] == 1.0 && fine["passes"] == 1.0;
    const bool nonuniform = coarse["velocity_grid_nonuniform"] == 1.0 &&
        base["velocity_grid_nonuniform"] == 1.0 &&
        fine["velocity_grid_nonuniform"] == 1.0;
    const bool transverse_grid_same = coarse["uperp_grid_hash"] ==
        base["uperp_grid_hash"] && base["uperp_grid_hash"] ==
        fine["uperp_grid_hash"];
    const bool same_build = coarse_revision == base_revision &&
        base_revision == fine_revision &&
        coarse["git_dirty"] == base["git_dirty"] &&
        base["git_dirty"] == fine["git_dirty"];
    const bool pair_improves =
        coarse["JN_minus_GstarJE_L2"] <=
            coarse["legacy_JN_minus_GstarJE_L2"] * (1.0 + 1.0e-10) &&
        base["JN_minus_GstarJE_L2"] <=
            base["legacy_JN_minus_GstarJE_L2"] * (1.0 + 1.0e-10) &&
        fine["JN_minus_GstarJE_L2"] <=
            fine["legacy_JN_minus_GstarJE_L2"] * (1.0 + 1.0e-10);
    const bool local_dual_closes =
        coarse["dual_cell_JE_minus_target_relative_L2"] <= 1.0e-10 &&
        base["dual_cell_JE_minus_target_relative_L2"] <= 1.0e-10 &&
        fine["dual_cell_JE_minus_target_relative_L2"] <= 1.0e-10 &&
        coarse["projection_reconstruction_relative_L2"] <= 1.0e-10 &&
        base["projection_reconstruction_relative_L2"] <= 1.0e-10 &&
        fine["projection_reconstruction_relative_L2"] <= 1.0e-10 &&
        coarse["dual_decomposition_pass"] == 1.0 &&
        base["dual_decomposition_pass"] == 1.0 &&
        fine["dual_decomposition_pass"] == 1.0;
    const double floor_match_coarse = relative_difference(
        coarse["JN_minus_GstarJE_L2"], coarse["projection_floor_L2"]);
    const double floor_match_base = relative_difference(
        base["JN_minus_GstarJE_L2"], base["projection_floor_L2"]);
    const double floor_match_fine = relative_difference(
        fine["JN_minus_GstarJE_L2"], fine["projection_floor_L2"]);
    const bool face_residual_is_projection_floor =
        floor_match_coarse <= 1.0e-10 && floor_match_base <= 1.0e-10 &&
        floor_match_fine <= 1.0e-10;
    const bool passes = resolution_valid && inputs_pass && nonuniform &&
                        transverse_grid_same && same_build && pair_improves &&
                        local_dual_closes && face_residual_is_projection_floor;
    std::ofstream out(output_path(argc, argv).c_str());
    std::ostream& log = out ? out : std::cout;
    log << std::scientific << std::setprecision(17)
        << "test=background_coupling_dual_u_convergence_audit\n"
        << "aggregator_process_count_hint=" << process_count_hint << "\n"
        << "resolution_sequence_valid=" << resolution_valid << "\n"
        << "inputs_pass=" << inputs_pass << "\n"
        << "nonuniform_grid_all_levels=" << nonuniform << "\n"
        << "transverse_grid_same=" << transverse_grid_same << "\n"
        << "same_git_build=" << same_build << "\n"
        << "git_revision=" << coarse_revision << "\n"
        << "dual_improves_legacy_all_levels=" << pair_improves << "\n"
        << "local_dual_closure_all_levels=" << local_dual_closes << "\n"
        << "face_residual_is_projection_floor="
        << face_residual_is_projection_floor << "\n"
        << "dual_L2_coarse=" << coarse["JN_minus_GstarJE_L2"] << "\n"
        << "dual_L2_base=" << base["JN_minus_GstarJE_L2"] << "\n"
        << "dual_L2_fine=" << fine["JN_minus_GstarJE_L2"] << "\n"
        << "projection_floor_L2_coarse=" << coarse["projection_floor_L2"]
        << "\nprojection_floor_L2_base=" << base["projection_floor_L2"]
        << "\nprojection_floor_L2_fine=" << fine["projection_floor_L2"]
        << "\nprojection_floor_match_relative_coarse=" << floor_match_coarse
        << "\nprojection_floor_match_relative_base=" << floor_match_base
        << "\nprojection_floor_match_relative_fine=" << floor_match_fine
        << "\nlocal_dual_relative_L2_coarse="
        << coarse["dual_cell_JE_minus_target_relative_L2"]
        << "\nlocal_dual_relative_L2_base="
        << base["dual_cell_JE_minus_target_relative_L2"]
        << "\nlocal_dual_relative_L2_fine="
        << fine["dual_cell_JE_minus_target_relative_L2"]
        << "\ndual_order_coarse_to_base=" << order(
            coarse["JN_minus_GstarJE_L2"], base["JN_minus_GstarJE_L2"],
            coarse["Nupar"], base["Nupar"]) << "\n"
        << "dual_order_base_to_fine=" << order(
            base["JN_minus_GstarJE_L2"], fine["JN_minus_GstarJE_L2"],
            base["Nupar"], fine["Nupar"]) << "\n"
        << "full_face_order_is_diagnostic_only=1\n"
        << "passes=" << passes << "\n";
    std::cout << "background_coupling_dual_u_convergence_audit passes="
              << passes << "\n";
    return passes ? 0 : 1;
}
