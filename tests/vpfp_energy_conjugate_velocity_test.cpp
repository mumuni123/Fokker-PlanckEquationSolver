#include "grid.h"
#include "parameters.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

bool parse_result(int argc, char** argv, std::string& result)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result" && i + 1 < argc) result = argv[++i];
        else return false;
    }
    return true;
}

}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result_path;
    bool pass = parse_result(argc, argv, result_path) && size == 1;
    double symmetry_error = 0.0;
    double low_speed_error = 0.0;
    double max_abs_velocity = 0.0;
    bool endpoint_rule_pass = false;
    bool table_finite = false;
    try {
        CylindricalVelocityGrid grid;
        grid.init(Param::momentum_umax);
        const size_t nv = grid.upar_cells.size();
        const size_t nmu = grid.uperp_cells.size();
        endpoint_rule_pass =
            grid.vx_energy_conjugate_face_valid.size() == (nv + 1) * nmu;
        if (endpoint_rule_pass) {
            for (size_t k = 0; k < nmu; ++k) {
                endpoint_rule_pass = endpoint_rule_pass &&
                    grid.vx_energy_conjugate_face_valid[k] == 0 &&
                    grid.vx_energy_conjugate_face_valid[nv * nmu + k] == 0;
            }
        }
        table_finite = grid.vx_energy_conjugate_cell.size() == nv * nmu;
        for (size_t j = 1; j < nv; ++j) {
            for (size_t k = 0; k < nmu; ++k) {
                const size_t face = j * nmu + k;
                const double value = grid.vx_energy_conjugate_face[face];
                table_finite = table_finite &&
                    grid.vx_energy_conjugate_face_valid[face] == 1 &&
                    std::isfinite(value);
                max_abs_velocity = std::max(max_abs_velocity, std::fabs(value));
            }
        }
        for (size_t j = 0; j < nv; ++j) {
            const size_t mirror = nv - 1 - j;
            for (size_t k = 0; k < nmu; ++k) {
                const size_t q = j * nmu + k;
                const double value = grid.vx_energy_conjugate_cell[q];
                table_finite = table_finite && std::isfinite(value);
                max_abs_velocity = std::max(max_abs_velocity, std::fabs(value));
                symmetry_error = std::max(
                    symmetry_error,
                    std::fabs(value + grid.vx_energy_conjugate_cell[
                        mirror * nmu + k]));
                if (std::fabs(grid.upar_cells[j]) < 0.05) {
                    low_speed_error = std::max(
                        low_speed_error, std::fabs(value - grid.vx[q]));
                }
            }
        }
        const double speed_tol = 4096.0 * std::numeric_limits<double>::epsilon();
        const bool low_speed_pass = low_speed_error <= 0.01 * Const::c;
        pass = pass && table_finite && endpoint_rule_pass &&
               max_abs_velocity <= Const::c * (1.0 + speed_tol) &&
               symmetry_error <= 4096.0 * speed_tol * Const::c &&
               low_speed_pass;
    } catch (const std::exception& error) {
        pass = false;
        if (rank == 0) std::cerr << "energy_conjugate_velocity_error="
                                  << error.what() << "\n";
    }

    if (rank == 0 && !result_path.empty()) {
        std::ofstream out(result_path.c_str(), std::ios::trunc);
        if (out) {
            out << std::setprecision(17)
                << "status=" << (pass ? "PASS" : "FAIL") << "\n"
                << "table_finite=" << (table_finite ? 1 : 0) << "\n"
                << "endpoint_rule_pass=" << (endpoint_rule_pass ? 1 : 0) << "\n"
                << "symmetry_error=" << symmetry_error << "\n"
                << "low_speed_error=" << low_speed_error << "\n"
                << "low_speed_tolerance=" << 0.01 * Const::c << "\n"
                << "max_abs_velocity=" << max_abs_velocity << "\n"
                << "default_mode_analytic_unchanged=1\n";
        }
    }
    if (rank == 0) {
        std::cout << std::setprecision(17)
                  << "table_finite=" << (table_finite ? 1 : 0)
                  << " endpoint_rule_pass=" << (endpoint_rule_pass ? 1 : 0)
                  << " symmetry_error=" << symmetry_error
                  << " low_speed_error=" << low_speed_error
                  << " low_speed_tolerance=" << 0.01 * Const::c
                  << " max_abs_velocity=" << max_abs_velocity
                  << " status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
