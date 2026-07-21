#include "checkpoint.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <string>
#include <vector>

namespace {
std::string option_value(int argc, char** argv, const char* option)
{
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == option) return argv[i + 1];
    return std::string();
}

void ensure_directory(const std::string& directory)
{
#ifdef _WIN32
    const std::string command = "mkdir \"" + directory + "\" >nul 2>nul";
#else
    const std::string command = "mkdir -p \"" + directory + "\"";
#endif
    (void)std::system(command.c_str());
}

struct Norm {
    double diff2;
    double ref2;
    double diff_inf;
    Norm() : diff2(0.0), ref2(0.0), diff_inf(0.0) {}
    void add(double reference, double value) {
        const double delta = value - reference;
        diff2 += delta * delta;
        ref2 += reference * reference;
        diff_inf = std::max(diff_inf, std::fabs(delta));
    }
};

double relative_l2(const Norm& value)
{
    return std::sqrt(value.diff2 / std::max(std::numeric_limits<double>::min(),
                                             value.ref2));
}
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::string legacy_dir = option_value(argc, argv,
        "--legacy-midpoint-audit-dir");
    const std::string dual_dir = option_value(argc, argv,
        "--dual-midpoint-audit-dir");
    const std::string output_dir = option_value(argc, argv, "--output-dir");
    if (legacy_dir.empty() || dual_dir.empty() || output_dir.empty()) {
        if (rank == 0) std::cerr << "usage: " << argv[0]
            << " --legacy-midpoint-audit-dir DIR --dual-midpoint-audit-dir DIR"
            << " --output-dir DIR\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid sg;
    sg.init(rank, size);
    Species template_species;
    template_species.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
                          -Const::qe, Const::me, Param::dens,
                          Param::temperature_e, false, sg);
    EMFields template_fields;
    template_fields.init(sg);
    VlasovAmpereMidpointSolver::MidpointAuditState legacy = {};
    VlasovAmpereMidpointSolver::MidpointAuditState dual = {};
    std::string error;
    const bool legacy_ok = read_midpoint_audit_state(legacy_dir, legacy,
        template_species, template_fields, sg, rank, size, error);
    const std::string legacy_error = error;
    error.clear();
    const bool dual_ok = read_midpoint_audit_state(dual_dir, dual,
        template_species, template_fields, sg, rank, size, error);
    const std::string dual_error = error;
    const int local_ok = legacy_ok && dual_ok ? 1 : 0;
    int all_ok = 0;
    MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!all_ok) {
        if (rank == 0) std::cerr << "FATAL stage4 A/B dump read failed: legacy="
            << legacy_error << " dual=" << dual_error << "\n";
        MPI_Finalize();
        return 3;
    }

    legacy.guess_np1.compute_moments();
    dual.guess_np1.compute_moments();
    const int nxl = sg.nx_local;
    const int ng = sg.nghost;
    Norm field, density, current;
    std::vector<double> rows(static_cast<size_t>(nxl) * 10, 0.0);
    std::vector<double> u_rows(static_cast<size_t>(Param::Nv) * 5, 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        const double x = sg.x(ng + ix);
        const bool core = x >= 0.1 * Const::micro &&
            x <= Param::Lx - 0.1 * Const::micro;
        const double e_a = legacy.fields_np1.Ex[ng + ix];
        const double e_b = dual.fields_np1.Ex[ng + ix];
        const double n_a = legacy.guess_np1.number_density[ix];
        const double n_b = dual.guess_np1.number_density[ix];
        const double j_a = 0.5 * (legacy.reference_jn_face[static_cast<size_t>(ix)] +
                                  legacy.reference_jn_face[static_cast<size_t>(ix + 1)]);
        const double j_b = 0.5 * (dual.reference_jn_face[static_cast<size_t>(ix)] +
                                  dual.reference_jn_face[static_cast<size_t>(ix + 1)]);
        if (core) {
            field.add(e_a, e_b);
            density.add(n_a, n_b);
            current.add(j_a, j_b);
        }
        const size_t base = static_cast<size_t>(ix) * 10;
        rows[base] = sg.ix_start + ix;
        rows[base + 1] = sg.x(ng + ix);
        rows[base + 2] = e_a; rows[base + 3] = e_b;
        rows[base + 4] = e_b - e_a;
        rows[base + 5] = n_a; rows[base + 6] = n_b;
        rows[base + 7] = n_b - n_a;
        rows[base + 8] = j_a; rows[base + 9] = j_b - j_a;
        for (int j = 0; j < Param::Nv; ++j) {
            double ma = 0.0, mb = 0.0;
            for (int k = 0; k < Param::Nmu; ++k) {
                ma += legacy.guess_np1.f[idx3(ng + ix, j, k)];
                mb += dual.guess_np1.f[idx3(ng + ix, j, k)];
            }
            const size_t ub = static_cast<size_t>(j) * 5;
            u_rows[ub] += ma;
            u_rows[ub + 1] += mb;
            u_rows[ub + 2] += std::fabs(mb - ma);
            u_rows[ub + 3] += std::fabs(ma);
            u_rows[ub + 4] += std::fabs(mb);
        }
    }
    double norm_values[9] = {field.diff2, field.ref2, field.diff_inf,
        density.diff2, density.ref2, density.diff_inf,
        current.diff2, current.ref2, current.diff_inf};
    MPI_Allreduce(MPI_IN_PLACE, norm_values, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, norm_values + 2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, norm_values + 3, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, norm_values + 5, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, norm_values + 6, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, norm_values + 8, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, u_rows.data(), static_cast<int>(u_rows.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    std::vector<double> gathered_rows;
    if (rank == 0) gathered_rows.assign(static_cast<size_t>(Param::nx) * 10, 0.0);
    MPI_Gather(rows.data(), static_cast<int>(rows.size()), MPI_DOUBLE,
               rank == 0 ? gathered_rows.data() : 0,
               static_cast<int>(rows.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        ensure_directory(output_dir);
        std::ofstream result((output_dir + "/stage4_ab_comparison.result").c_str());
        std::ofstream by_x((output_dir + "/stage4_ab_by_x.dat").c_str());
        std::ofstream by_u((output_dir + "/stage4_ab_by_upar.dat").c_str());
        const double ex_relative = std::sqrt(norm_values[0] /
            std::max(std::numeric_limits<double>::min(), norm_values[1]));
        const double density_relative = std::sqrt(norm_values[3] /
            std::max(std::numeric_limits<double>::min(), norm_values[4]));
        const double current_relative = std::sqrt(norm_values[6] /
            std::max(std::numeric_limits<double>::min(), norm_values[7]));
        const bool legacy_strict = legacy.acceptance_type == "strict_converged";
        const bool dual_strict = dual.acceptance_type == "strict_converged";
        const bool legacy_soft = legacy.acceptance_type == "soft_accepted";
        const bool dual_soft = dual.acceptance_type == "soft_accepted";
        const bool core_state_within_gate = ex_relative <= 2.0e-2 &&
            density_relative <= 2.0e-2 && current_relative <= 5.0e-2;
        const bool final_iteration_not_increased =
            dual.nonlinear_iterations <= legacy.nonlinear_iterations;
        const bool mode_pair_valid =
            legacy.background_coupling_mode ==
                VlasovAmpereMidpointSolver::LEGACY_COUPLING &&
            dual.background_coupling_mode ==
                VlasovAmpereMidpointSolver::DUAL_U_COUPLING;
        const bool strict_endpoint_pass = mode_pair_valid &&
            core_state_within_gate && legacy_strict && dual_strict &&
            final_iteration_not_increased;
        const bool matched_soft_comparison_pass = mode_pair_valid &&
            core_state_within_gate && legacy_soft && dual_soft &&
            final_iteration_not_increased;
        const bool comparable_endpoint_pass = strict_endpoint_pass ||
            matched_soft_comparison_pass;
        result << std::setprecision(17)
            << "legacy_mode " << legacy.background_coupling_mode << "\n"
            << "dual_mode " << dual.background_coupling_mode << "\n"
            << "same_initial_time " << (legacy.time_s == dual.time_s ? 1 : 0) << "\n"
            << "legacy_nonlinear_iterations " << legacy.nonlinear_iterations << "\n"
            << "dual_nonlinear_iterations " << dual.nonlinear_iterations << "\n"
            << "legacy_acceptance " << legacy.acceptance_type << "\n"
            << "dual_acceptance " << dual.acceptance_type << "\n"
            << "core_Ex_relative_L2 " << ex_relative << "\n"
            << "core_n_relative_L2 " << density_relative << "\n"
            << "core_J_bkg_relative_L2 " << current_relative << "\n"
            << "Ex_difference_Linf " << norm_values[2] << "\n"
            << "n_difference_Linf " << norm_values[5] << "\n"
            << "J_bkg_difference_Linf " << norm_values[8] << "\n"
            << "legacy_limiter_active_fraction " << legacy.limiter_active_fraction << "\n"
            << "dual_limiter_active_fraction " << dual.limiter_active_fraction << "\n"
            << "legacy_limiter_active_fraction_core " << legacy.limiter_active_fraction_core << "\n"
            << "dual_limiter_active_fraction_core " << dual.limiter_active_fraction_core << "\n"
            << "legacy_limiter_active_fraction_boundary " << legacy.limiter_active_fraction_boundary << "\n"
            << "dual_limiter_active_fraction_boundary " << dual.limiter_active_fraction_boundary << "\n"
            << "legacy_x_limiter_active_fraction " << legacy.x_limiter_active_fraction << "\n"
            << "dual_x_limiter_active_fraction " << dual.x_limiter_active_fraction << "\n"
            << "legacy_u_limiter_active_fraction " << legacy.u_limiter_active_fraction << "\n"
            << "dual_u_limiter_active_fraction " << dual.u_limiter_active_fraction << "\n"
            << "legacy_limiter_min_alpha " << legacy.limiter_min_alpha << "\n"
            << "dual_limiter_min_alpha " << dual.limiter_min_alpha << "\n"
            << "core_state_gate_Ex_n_2pct_J_5pct "
            << (core_state_within_gate ? 1 : 0) << "\n"
            << "legacy_final_strict_converged " << (legacy_strict ? 1 : 0) << "\n"
            << "dual_final_strict_converged " << (dual_strict ? 1 : 0) << "\n"
            << "dual_final_iteration_not_increased "
            << (final_iteration_not_increased ? 1 : 0) << "\n"
            << "mode_pair_valid " << (mode_pair_valid ? 1 : 0) << "\n"
            << "matched_acceptance_class "
            << ((legacy_strict && dual_strict) || (legacy_soft && dual_soft) ? 1 : 0)
            << "\n"
            << "stage4_strict_endpoint_pass "
            << (strict_endpoint_pass ? 1 : 0) << "\n"
            << "stage4_matched_soft_comparison_pass "
            << (matched_soft_comparison_pass ? 1 : 0) << "\n"
            << "stage4_comparable_endpoint_pass "
            << (comparable_endpoint_pass ? 1 : 0) << "\n"
            << "stage4_requires_production_short_run_confirmation "
            << (matched_soft_comparison_pass ? 1 : 0) << "\n"
            << "stage4_final_endpoint_pass "
            << (strict_endpoint_pass ? 1 : 0) << "\n";
        by_x << "# global_ix x_m Ex_legacy Ex_dual delta_Ex n_legacy n_dual delta_n J_bkg_legacy delta_J_bkg\n";
        for (int ix = 0; ix < Param::nx; ++ix) {
            const size_t base = static_cast<size_t>(ix) * 10;
            for (int c = 0; c < 10; ++c)
                by_x << gathered_rows[base + c] << (c == 9 ? '\n' : ' ');
        }
        by_u << "# iu upar mass_legacy mass_dual abs_delta_mass abs_mass_legacy abs_mass_dual zone\n";
        for (int j = 0; j < Param::Nv; ++j) {
            const double u = legacy.guess_np1.cgrid.upar_cells[j];
            const char* zone = std::fabs(std::fabs(u) - 1.0) <= 0.15
                ? "near_abs_upar_1" :
                (std::fabs(u) <= 1.0 ? "thermal_body" : "outer");
            const size_t ub = static_cast<size_t>(j) * 5;
            by_u << j << ' ' << u << ' ' << u_rows[ub] << ' ' << u_rows[ub + 1]
                 << ' ' << u_rows[ub + 2] << ' ' << u_rows[ub + 3] << ' '
                 << u_rows[ub + 4] << ' ' << zone << '\n';
        }
    }
    MPI_Finalize();
    return 0;
}
