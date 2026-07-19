#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>

namespace {

bool read_result(const char* path, std::map<std::string, double>& values)
{
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        const std::string::size_type equal = line.find('=');
        if (equal == std::string::npos) continue;
        char* end = 0;
        const double value = std::strtod(line.c_str() + equal + 1, &end);
        if (end != line.c_str() + equal + 1)
            values[line.substr(0, equal)] = value;
    }
    return values.count("nx") && values.count("Nu") &&
        values.count("Nuperp") && values.count("dt") &&
        values.count("current_L2") && values.count("work_L2") &&
        values.count("fct_enabled") && values.count("fct_active") &&
        values.count("state_advanced") && values.count("operator_failed") &&
        values.count("outputs_finite") && values.count("audit_valid") &&
        values.count("passes");
}

double observed_order(double coarse_error, double fine_error,
                      double coarse_h, double fine_h)
{
    if (coarse_error <= 0.0 || fine_error <= 0.0 || coarse_h <= fine_h)
        return 0.0;
    return std::log(coarse_error / fine_error) /
        std::log(coarse_h / fine_h);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "usage: background_coupling_directional_grid_audit "
                  << "kind coarse.result base.result fine.result\n";
        return 2;
    }
    std::map<std::string, double> coarse, base, fine;
    if (!read_result(argv[2], coarse) || !read_result(argv[3], base) ||
        !read_result(argv[4], fine)) {
        std::cerr << "unable to read directional grid result files\n";
        return 2;
    }
    const std::string kind(argv[1]);
    const char* resolution_key = kind == "x" ? "nx" :
        (kind == "u_parallel" ? "Nu" : "Nuperp");
    const double h_coarse = 1.0 / coarse[resolution_key];
    const double h_base = 1.0 / base[resolution_key];
    const double h_fine = 1.0 / fine[resolution_key];
    const int fct_consistent = coarse["fct_enabled"] == base["fct_enabled"] &&
        base["fct_enabled"] == fine["fct_enabled"] &&
        coarse["fct_active"] == base["fct_active"] &&
        base["fct_active"] == fine["fct_active"];
    const int fixed_other_grids = kind == "x"
        ? (coarse["Nu"] == base["Nu"] && base["Nu"] == fine["Nu"] &&
           coarse["Nuperp"] == base["Nuperp"] &&
           base["Nuperp"] == fine["Nuperp"])
        : (coarse["nx"] == base["nx"] && base["nx"] == fine["nx"] &&
           (kind != "u_parallel" || (coarse["Nuperp"] == base["Nuperp"] &&
                                 base["Nuperp"] == fine["Nuperp"])) &&
           (kind != "u_perp" || (coarse["Nu"] == base["Nu"] &&
                                  base["Nu"] == fine["Nu"])));
    const int inputs_valid = coarse["passes"] == 1.0 &&
        base["passes"] == 1.0 && fine["passes"] == 1.0 &&
        coarse["audit_valid"] == 1.0 && base["audit_valid"] == 1.0 &&
        fine["audit_valid"] == 1.0 &&
        coarse["state_advanced"] == 1.0 && base["state_advanced"] == 1.0 &&
        fine["state_advanced"] == 1.0 &&
        coarse["operator_failed"] == 0.0 && base["operator_failed"] == 0.0 &&
        fine["operator_failed"] == 0.0 &&
        coarse["outputs_finite"] == 1.0 && base["outputs_finite"] == 1.0 &&
        fine["outputs_finite"] == 1.0;
    const double eps = 64.0 * std::numeric_limits<double>::epsilon();
    const int dt_fixed = std::fabs(coarse["dt"] - base["dt"]) <=
            eps * std::max(1.0, std::max(std::fabs(coarse["dt"]), std::fabs(base["dt"]))) &&
        std::fabs(base["dt"] - fine["dt"]) <=
            eps * std::max(1.0, std::max(std::fabs(base["dt"]), std::fabs(fine["dt"])));
    const int finite_nonzero = std::isfinite(coarse["current_L2"]) &&
        std::isfinite(base["current_L2"]) && std::isfinite(fine["current_L2"]) &&
        std::isfinite(coarse["work_L2"]) && std::isfinite(base["work_L2"]) &&
        std::isfinite(fine["work_L2"]) && coarse["current_L2"] > 0.0 &&
        base["current_L2"] > 0.0 && fine["current_L2"] > 0.0 &&
        coarse["work_L2"] > 0.0 && base["work_L2"] > 0.0 && fine["work_L2"] > 0.0;
    const int monotone_current = base["current_L2"] < coarse["current_L2"] &&
        fine["current_L2"] < base["current_L2"];
    const int monotone_work = base["work_L2"] < coarse["work_L2"] &&
        fine["work_L2"] < base["work_L2"];
    const double current_order_fine = observed_order(
        base["current_L2"], fine["current_L2"], h_base, h_fine);
    const double work_order_fine = observed_order(
        base["work_L2"], fine["work_L2"], h_base, h_fine);
    const int asymptotic_order_established = monotone_current && monotone_work &&
        current_order_fine >= 0.5 && work_order_fine >= 0.5;
    const int passes = inputs_valid && fixed_other_grids && fct_consistent &&
        dt_fixed && finite_nonzero && monotone_current && monotone_work;
    std::ofstream out((std::string("output/background_coupling_directional_") +
                       kind + "_audit.result").c_str());
    std::ostream& log = out ? out : std::cout;
    log << "test=background_coupling_directional_grid_audit\n"
        << "refinement_kind=" << kind << "\n"
        << "audit_inputs_valid=" << inputs_valid << "\n"
        << "fct_state_consistent=" << fct_consistent << "\n"
        << "other_grids_fixed=" << fixed_other_grids << "\n"
        << "dt_fixed=" << dt_fixed << "\n"
        << "finite_nonzero=" << finite_nonzero << "\n"
        << "current_L2_order_coarse_to_base=" << observed_order(
            coarse["current_L2"], base["current_L2"], h_coarse, h_base) << "\n"
        << "current_L2_order_base_to_fine=" << current_order_fine << "\n"
        << "work_L2_order_coarse_to_base=" << observed_order(
            coarse["work_L2"], base["work_L2"], h_coarse, h_base) << "\n"
        << "work_L2_order_base_to_fine=" << work_order_fine << "\n"
        << "current_L2_monotone=" << monotone_current << "\n"
        << "work_L2_monotone=" << monotone_work << "\n"
        << "asymptotic_order_established=" << asymptotic_order_established << "\n"
        << "passes=" << passes << "\n";
    std::cout << "background_coupling_directional_grid_audit passes="
              << passes << "\n";
    return passes ? 0 : 1;
}
