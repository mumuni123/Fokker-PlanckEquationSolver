#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

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
    return values.count("passes") && values.count("mpi_size") &&
           values.count("git_dirty") && !revision.empty();
}

std::string output_path(int argc, char** argv)
{
    for (int i = 3; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--output") return argv[i + 1];
    return "output/dual_u_stage2/mpi_consistency/summary.result";
}

double relative_difference(double a, double b)
{
    return std::fabs(a - b) /
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: background_coupling_dual_u_mpi_consistency_audit "
                  << "single.result multi.result [--output summary.result]\n";
        return 2;
    }
    std::map<std::string, double> single, multi;
    std::string single_revision, multi_revision;
    if (!read_values(argv[1], single, single_revision) ||
        !read_values(argv[2], multi, multi_revision)) {
        std::cerr << "unable to read MPI consistency results\n";
        return 2;
    }
    const char* relative_keys[] = {
        "JN_minus_GstarJE_L1", "JN_minus_GstarJE_L2",
        "JN_minus_GstarJE_Linf", "legacy_dual_JN_Linf",
        "legacy_dual_JE_Linf", "legacy_dual_state_relative_L1",
        "dual_cell_JE_minus_target_relative_L2",
        "projection_floor_L2", "projection_reconstruction_relative_L2"};
    const char* absolute_keys[] = {"JN_minus_GstarJE_signed", "R_FV",
        "R_couple", "u_boundary_work_left_plus_right"};
    double maximum = 0.0;
    std::string maximum_key = "none";
    bool fields_present = true;
    std::vector<std::pair<std::string, double> > differences;
    for (size_t i = 0; i < sizeof(relative_keys) /
                           sizeof(relative_keys[0]); ++i) {
        const char* key = relative_keys[i];
        if (!single.count(key) || !multi.count(key)) {
            fields_present = false;
            continue;
        }
        const double difference = relative_difference(single[key], multi[key]);
        differences.push_back(std::make_pair(std::string(key), difference));
        if (difference > maximum) {
            maximum = difference;
            maximum_key = key;
        }
    }
    double maximum_absolute_near_zero_difference = 0.0;
    std::string maximum_absolute_key = "none";
    for (size_t i = 0; i < sizeof(absolute_keys) /
                           sizeof(absolute_keys[0]); ++i) {
        const char* key = absolute_keys[i];
        if (!single.count(key) || !multi.count(key)) {
            fields_present = false;
            continue;
        }
        const double difference = std::fabs(single[key] - multi[key]);
        differences.push_back(std::make_pair(std::string(key), difference));
        if (difference > maximum_absolute_near_zero_difference) {
            maximum_absolute_near_zero_difference = difference;
            maximum_absolute_key = key;
        }
    }
    const bool conservation_fields_present =
        single.count("mass_error_relative") &&
        multi.count("mass_error_relative") &&
        single.count("momentum_error_relative") &&
        multi.count("momentum_error_relative");
    fields_present = fields_present && conservation_fields_present;
    const double maximum_mass_relative = conservation_fields_present ?
        std::max(single["mass_error_relative"],
                 multi["mass_error_relative"]) :
        std::numeric_limits<double>::infinity();
    const double maximum_momentum_relative = conservation_fields_present ?
        std::max(single["momentum_error_relative"],
                 multi["momentum_error_relative"]) :
        std::numeric_limits<double>::infinity();
    const bool conservation_valid = maximum_mass_relative <= 1.0e-10 &&
                                    maximum_momentum_relative <= 1.0e-9;
    const bool rank_layout_valid = single["mpi_size"] == 1.0 &&
                                   multi["mpi_size"] > 1.0;
    const bool inputs_pass = single["passes"] == 1.0 &&
                             multi["passes"] == 1.0;
    const bool same_build = single_revision == multi_revision &&
        single["git_dirty"] == multi["git_dirty"];
    const bool passes = fields_present && rank_layout_valid && inputs_pass &&
                        same_build && conservation_valid &&
                        maximum <= 1.0e-10 &&
                        maximum_absolute_near_zero_difference <= 1.0e-10;
    std::ofstream out(output_path(argc, argv).c_str());
    std::ostream& log = out ? out : std::cout;
    log << std::scientific << std::setprecision(17)
        << "test=background_coupling_dual_u_mpi_consistency_audit\n"
        << "single_mpi_size=" << single["mpi_size"] << "\n"
        << "multi_mpi_size=" << multi["mpi_size"] << "\n"
        << "fields_present=" << fields_present << "\n"
        << "same_git_build=" << same_build << "\n"
        << "git_revision=" << single_revision << "\n"
        << "conservation_valid=" << conservation_valid << "\n"
        << "maximum_mass_error_relative=" << maximum_mass_relative << "\n"
        << "maximum_momentum_error_relative=" << maximum_momentum_relative
        << "\n"
        << "maximum_relative_difference=" << maximum << "\n"
        << "maximum_relative_difference_key=" << maximum_key << "\n"
        << "maximum_absolute_near_zero_difference="
        << maximum_absolute_near_zero_difference << "\n"
        << "maximum_absolute_near_zero_difference_key="
        << maximum_absolute_key << "\n";
    for (size_t i = 0; i < differences.size(); ++i)
        log << "difference_" << differences[i].first << "="
            << differences[i].second << "\n";
    log
        << "passes=" << passes << "\n";
    std::cout << "background_coupling_dual_u_mpi_consistency_audit passes="
              << passes << "\n";
    return passes ? 0 : 1;
}
