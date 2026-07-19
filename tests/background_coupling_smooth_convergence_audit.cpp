#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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
        const std::string key = line.substr(0, equal);
        const char* begin = line.c_str() + equal + 1;
        char* end = 0;
        const double value = std::strtod(begin, &end);
        if (end != begin) values[key] = value;
    }
    return values.count("nx") != 0 && values.count("full_L2") != 0 &&
        values.count("field_work_scale_L2") != 0;
}

double observed_order(double coarse_error, double fine_error,
                      double coarse_dx, double fine_dx)
{
    if (coarse_error <= 0.0 || fine_error <= 0.0 || coarse_dx <= fine_dx)
        return 0.0;
    return std::log(coarse_error / fine_error) /
        std::log(coarse_dx / fine_dx);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: background_coupling_smooth_convergence_audit "
                  << "coarse.result base.result fine.result\n";
        return 2;
    }
    std::map<std::string, double> coarse, base, fine;
    if (!read_result(argv[1], coarse) || !read_result(argv[2], base) ||
        !read_result(argv[3], fine)) {
        std::cerr << "unable to read a smooth-periodic result file\n";
        return 2;
    }
    const double h_coarse = 1.0 / coarse["nx"];
    const double h_base = 1.0 / base["nx"];
    const double h_fine = 1.0 / fine["nx"];
    const double p_coarse_base = observed_order(coarse["full_L2"],
                                                base["full_L2"],
                                                h_coarse, h_base);
    const double p_base_fine = observed_order(base["full_L2"], fine["full_L2"],
                                              h_base, h_fine);
    // The accepted amplitude diagnostic is the local field-work residual
    // dt*dx*E*(JN-G*JE), not the bare current mismatch.  The latter can have
    // a field-independent transport component in a manufactured state.
    const bool amplitude_linear =
        std::fabs(base["field_work_scale_L2"] - 2.0) < 0.25;
    const bool phase_equivalent = base["phase_shift_relative_L2"] < 0.05;
    const bool errors_decrease = base["full_L2"] < coarse["full_L2"] &&
        fine["full_L2"] < base["full_L2"];
    const bool passes = amplitude_linear && phase_equivalent && errors_decrease;
    std::ofstream out("output/background_coupling_smooth_convergence_audit.result");
    std::ostream& log = out ? out : std::cout;
    log << std::scientific << std::setprecision(17)
        << "test=background_coupling_smooth_convergence_audit\n"
        << "coarse_L2=" << coarse["full_L2"] << "\n"
        << "base_L2=" << base["full_L2"] << "\n"
        << "fine_L2=" << fine["full_L2"] << "\n"
        << "observed_order_coarse_to_base=" << p_coarse_base << "\n"
        << "observed_order_base_to_fine=" << p_base_fine << "\n"
        << "base_current_difference_scale_L2="
        << base["current_difference_scale_L2"] << "\n"
        << "base_field_work_scale_L2=" << base["field_work_scale_L2"] << "\n"
        << "base_phase_shift_relative_L2=" << base["phase_shift_relative_L2"] << "\n"
        << "errors_decrease=" << errors_decrease << "\n"
        << "passes=" << passes << "\n";
    std::cout << "background_coupling_smooth_convergence_audit passes="
              << passes << "\n";
    return passes ? 0 : 1;
}
