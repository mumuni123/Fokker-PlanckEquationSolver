#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace {

int create_directory(const char* path)
{
#ifdef _WIN32
    return mkdir(path);
#else
    return mkdir(path, 0777);
#endif
}

struct Options {
    std::string audit_dir;
    std::string output_dir;
    double sigma_cutoff;
    double lambda;
    double eta;
    double trust_fraction;
    double correction_trust_fraction;
    double energy_pair_tolerance;
    double energy_residual_fraction;
    double mass_relative_tolerance;
    double f_residual_growth_tolerance;

    Options()
        : sigma_cutoff(1.0e-6), lambda(1.0e-2), eta(1.0e-4),
          trust_fraction(0.1), correction_trust_fraction(1.0),
          energy_pair_tolerance(1.0e-8), energy_residual_fraction(1.0),
          mass_relative_tolerance(1.0e-10),
          f_residual_growth_tolerance(0.1) {}
};

bool parse_double(const char* text, double& value)
{
    char* end = 0;
    value = std::strtod(text, &end);
    return text && *text && end && *end == '\0' && std::isfinite(value);
}

bool parse_options(int argc, char** argv, Options& options,
                   std::string& error)
{
    for (int i = 1; i < argc && error.empty(); ++i) {
        const std::string arg(argv[i]);
        if (i + 1 >= argc) {
            error = "missing value for " + arg;
            break;
        }
        const char* value = argv[++i];
        if (arg == "--operator-audit") options.audit_dir = value;
        else if (arg == "--output-dir") options.output_dir = value;
        else if (arg == "--sigma-cutoff") {
            if (!parse_double(value, options.sigma_cutoff))
                error = "invalid --sigma-cutoff";
        } else if (arg == "--lambda") {
            if (!parse_double(value, options.lambda))
                error = "invalid --lambda";
        } else if (arg == "--eta") {
            if (!parse_double(value, options.eta))
                error = "invalid --eta";
        } else if (arg == "--trust-fraction") {
            if (!parse_double(value, options.trust_fraction))
                error = "invalid --trust-fraction";
        } else if (arg == "--correction-trust-fraction") {
            if (!parse_double(value, options.correction_trust_fraction))
                error = "invalid --correction-trust-fraction";
        } else if (arg == "--energy-pair-tolerance") {
            if (!parse_double(value, options.energy_pair_tolerance))
                error = "invalid --energy-pair-tolerance";
        } else if (arg == "--energy-residual-fraction") {
            if (!parse_double(value, options.energy_residual_fraction))
                error = "invalid --energy-residual-fraction";
        } else if (arg == "--mass-relative-tolerance") {
            if (!parse_double(value, options.mass_relative_tolerance))
                error = "invalid --mass-relative-tolerance";
        } else if (arg == "--f-residual-growth-tolerance") {
            if (!parse_double(value, options.f_residual_growth_tolerance))
                error = "invalid --f-residual-growth-tolerance";
        } else {
            error = "unknown option: " + arg;
        }
    }
    if (error.empty() && options.audit_dir.empty())
        error = "--operator-audit is required";
    if (error.empty() && options.output_dir.empty())
        error = "--output-dir is required";
    if (error.empty() &&
        (!(options.sigma_cutoff >= 0.0 && options.sigma_cutoff < 1.0) ||
         !(options.lambda >= 0.0) || !(options.eta > 0.0) ||
         !(options.trust_fraction > 0.0 &&
           options.trust_fraction <= 1.0)))
        error = "invalid regularized pairing parameter";
    if (error.empty() &&
        (!(options.correction_trust_fraction > 0.0) ||
         !(options.energy_pair_tolerance > 0.0) ||
         !(options.energy_residual_fraction > 0.0) ||
         !(options.mass_relative_tolerance > 0.0) ||
         !(options.f_residual_growth_tolerance >= 0.0)))
        error = "invalid regularized pairing acceptance parameter";
    return error.empty();
}

double vector_linf_difference(const std::vector<double>& left,
                              const std::vector<double>& right)
{
    if (left.size() != right.size())
        return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (size_t i = 0; i < left.size(); ++i)
        result = std::max(result, std::fabs(left[i] - right[i]));
    return result;
}

struct FaceNorms {
    double global_l2;
    double core_l2;
};

FaceNorms pairing_norms(
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation& result,
    const SpatialGrid& grid)
{
    long double global_square = 0.0L;
    long double core_square = 0.0L;
    const double boundary_width = 0.1 * Const::micro;
    for (int iface = 0; iface < grid.nx_local; ++iface) {
        const double residual =
            result.j_bkg_face_mid[static_cast<size_t>(iface)] -
            result.j_bkg_energy_debug_face[static_cast<size_t>(iface)];
        global_square += static_cast<long double>(residual) * residual;
        const double x = (grid.ix_start + iface) * grid.dx;
        if (x >= boundary_width && x <= Param::Lx - boundary_width)
            core_square += static_cast<long double>(residual) * residual;
    }
    double squares[2] = {
        static_cast<double>(global_square),
        static_cast<double>(core_square)};
    MPI_Allreduce(MPI_IN_PLACE, squares, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    FaceNorms norms = {
        std::sqrt(std::max(0.0, squares[0])),
        std::sqrt(std::max(0.0, squares[1]))};
    return norms;
}

double global_background_mass(
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation& result,
    const SpatialGrid&)
{
    double local = 0.0;
    double local_energy = 0.0;
    result.species_np1.total_particle_number_and_energy(
        local, local_energy);
    double global = local;
    MPI_Allreduce(MPI_IN_PLACE, &global, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return global;
}

void gather_unique_faces(const std::vector<double>& local, int local_count,
                         int rank, int size, std::vector<double>& global)
{
    std::vector<int> counts(static_cast<size_t>(size), 0);
    MPI_Gather(&local_count, 1, MPI_INT,
               rank == 0 ? counts.data() : 0, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    std::vector<int> displacements(static_cast<size_t>(size), 0);
    int total = 0;
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            displacements[static_cast<size_t>(r)] = total;
            total += counts[static_cast<size_t>(r)];
        }
        global.resize(static_cast<size_t>(total));
    }
    MPI_Gatherv(local.data(), local_count, MPI_DOUBLE,
                rank == 0 ? global.data() : 0,
                rank == 0 ? counts.data() : 0,
                rank == 0 ? displacements.data() : 0,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Options options;
    std::string error;
    int options_valid = 1;
    if (rank == 0)
        options_valid = parse_options(argc, argv, options, error) ? 1 : 0;
    MPI_Bcast(&options_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double parameters[9] = {
        options.sigma_cutoff, options.lambda, options.eta,
        options.trust_fraction, options.correction_trust_fraction,
        options.energy_pair_tolerance, options.energy_residual_fraction,
        options.mass_relative_tolerance,
        options.f_residual_growth_tolerance};
    MPI_Bcast(parameters, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int audit_length = rank == 0
        ? static_cast<int>(options.audit_dir.size()) : 0;
    int output_length = rank == 0
        ? static_cast<int>(options.output_dir.size()) : 0;
    MPI_Bcast(&audit_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&output_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        options.audit_dir.resize(static_cast<size_t>(audit_length));
        options.output_dir.resize(static_cast<size_t>(output_length));
    }
    if (audit_length)
        MPI_Bcast(&options.audit_dir[0], audit_length, MPI_CHAR, 0,
                  MPI_COMM_WORLD);
    if (output_length)
        MPI_Bcast(&options.output_dir[0], output_length, MPI_CHAR, 0,
                  MPI_COMM_WORLD);
    options.sigma_cutoff = parameters[0];
    options.lambda = parameters[1];
    options.eta = parameters[2];
    options.trust_fraction = parameters[3];
    options.correction_trust_fraction = parameters[4];
    options.energy_pair_tolerance = parameters[5];
    options.energy_residual_fraction = parameters[6];
    options.mass_relative_tolerance = parameters[7];
    options.f_residual_growth_tolerance = parameters[8];
    if (!options_valid) {
        if (rank == 0)
            std::fprintf(stderr, "face_pairing_checkpoint_audit: %s\n",
                         error.c_str());
        MPI_Finalize();
        return 2;
    }

    int output_ready = 1;
    if (rank == 0 && create_directory(options.output_dir.c_str()) != 0 &&
        errno != EEXIST) {
        std::fprintf(stderr, "cannot create output directory %s: %s\n",
                     options.output_dir.c_str(), std::strerror(errno));
        output_ready = 0;
    }
    MPI_Bcast(&output_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!output_ready) {
        MPI_Finalize();
        return 3;
    }

    SpatialGrid grid;
    grid.init(rank, size);
    Species species_template;
    species_template.init(
        "bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
        Param::dens, Param::temperature_e, false, grid);
    species_template.initialize_maxwellian();
    EMFields fields_template;
    fields_template.init(grid);
    VlasovAmpereMidpointSolver::MidpointAuditState state;
    if (!read_midpoint_audit_state(
            options.audit_dir, state, species_template, fields_template,
            grid, rank, size, error)) {
        if (rank == 0)
            std::fprintf(stderr, "midpoint audit read failed: %s\n",
                         error.c_str());
        MPI_Finalize();
        return 4;
    }

    BeamPIC beam;
    beam.init(grid);
    const auto configure_solver = [&](VlasovAmpereMidpointSolver& solver,
                                      bool regularized) {
        solver.set_low_order_only(state.low_order_only);
        solver.set_nonuniform_high_order_enabled(state.high_order_enabled);
        solver.set_fct_enabled(state.fct_enabled);
        solver.set_background_coupling_mode(
            VlasovAmpereMidpointSolver::DUAL_U_COUPLING);
        solver.set_face_pairing_mode(regularized
            ? VlasovAmpereMidpointSolver::FACE_PAIRING_REGULARIZED
            : VlasovAmpereMidpointSolver::FACE_PAIRING_CELL_BASELINE);
        solver.set_face_pairing_parameters(
            options.sigma_cutoff, options.lambda, options.eta,
            options.trust_fraction);
        solver.set_face_pairing_acceptance_parameters(
            options.correction_trust_fraction,
            options.energy_pair_tolerance,
            options.energy_residual_fraction,
            options.mass_relative_tolerance,
            options.f_residual_growth_tolerance);
        solver.set_max_midpoint_iterations(1);
    };

    VlasovAmpereMidpointSolver baseline_solver;
    VlasovAmpereMidpointSolver regularized_solver;
    configure_solver(baseline_solver, false);
    configure_solver(regularized_solver, true);
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation baseline =
        baseline_solver.evaluate_fixed_midpoint_operator(
            state.bkg_n, beam, state.fields_n, state.operator_input_guess,
            state.fields_end_guess, state.j_beam_face_mid,
            state.coupling_layout, grid, state.dt_s, state.time_s, rank,
            size);
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation regularized =
        regularized_solver.evaluate_fixed_midpoint_operator(
            state.bkg_n, beam, state.fields_n, state.operator_input_guess,
            state.fields_end_guess, state.j_beam_face_mid,
            state.coupling_layout, grid, state.dt_s, state.time_s, rank,
            size);

    const FaceNorms baseline_norms = pairing_norms(baseline, grid);
    const FaceNorms regularized_norms = pairing_norms(regularized, grid);
    double jn_difference = vector_linf_difference(
        baseline.j_bkg_face_mid, regularized.j_bkg_face_mid);
    double field_difference = vector_linf_difference(
        baseline.fields_np1.Ex_face, regularized.fields_np1.Ex_face);
    MPI_Allreduce(MPI_IN_PLACE, &jn_difference, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &field_difference, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double baseline_mass = global_background_mass(baseline, grid);
    const double regularized_mass = global_background_mass(regularized, grid);
    const double mass_difference = regularized_mass - baseline_mass;
    const double delta_ke =
        regularized.delta_ke_bkg - baseline.delta_ke_bkg;
    const double delta_work =
        regularized.field_work_bkg - baseline.field_work_bkg;

    const double jn_scale = std::max(
        1.0, baseline.current_diag.max_abs_j_charge);
    double field_scale = 1.0;
    for (size_t i = 0; i < baseline.fields_np1.Ex_face.size(); ++i)
        field_scale = std::max(
            field_scale, std::fabs(baseline.fields_np1.Ex_face[i]));
    MPI_Allreduce(MPI_IN_PLACE, &field_scale, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double mass_tolerance =
        4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::fabs(baseline_mass));
    const bool pass =
        baseline.state_advanced && !baseline.failed &&
        regularized.state_advanced && !regularized.failed &&
        regularized.face_pairing_attempted &&
        regularized.face_pairing_accepted &&
        !regularized.face_pairing_fallback_to_cell_baseline &&
        regularized_norms.global_l2 < baseline_norms.global_l2 &&
        regularized_norms.core_l2 < baseline_norms.core_l2 &&
        jn_difference <= 4096.0 *
            std::numeric_limits<double>::epsilon() * jn_scale &&
        field_difference <= 4096.0 *
            std::numeric_limits<double>::epsilon() * field_scale &&
        std::fabs(mass_difference) <= mass_tolerance;
    int global_pass = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_pass, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    std::vector<double> baseline_jn;
    std::vector<double> baseline_gstar;
    std::vector<double> regularized_gstar;
    gather_unique_faces(baseline.j_bkg_face_mid, grid.nx_local,
                        rank, size, baseline_jn);
    gather_unique_faces(baseline.j_bkg_energy_debug_face, grid.nx_local,
                        rank, size, baseline_gstar);
    gather_unique_faces(regularized.j_bkg_energy_debug_face, grid.nx_local,
                        rank, size, regularized_gstar);
    if (rank == 0) {
        std::ofstream summary(
            (options.output_dir +
             "/face_pairing_checkpoint_audit.result").c_str());
        summary << std::scientific << std::setprecision(16)
            << "sigma_cutoff " << options.sigma_cutoff << "\n"
            << "lambda " << options.lambda << "\n"
            << "eta " << options.eta << "\n"
            << "trust_fraction " << options.trust_fraction << "\n"
            << "correction_trust_fraction "
            << options.correction_trust_fraction << "\n"
            << "energy_pair_tolerance "
            << options.energy_pair_tolerance << "\n"
            << "energy_residual_fraction "
            << options.energy_residual_fraction << "\n"
            << "mass_relative_tolerance "
            << options.mass_relative_tolerance << "\n"
            << "f_residual_growth_tolerance "
            << options.f_residual_growth_tolerance << "\n"
            << "face_pairing_attempted "
            << regularized.face_pairing_attempted << "\n"
            << "face_pairing_accepted "
            << regularized.face_pairing_accepted << "\n"
            << "fallback_to_cell_baseline "
            << regularized.face_pairing_fallback_to_cell_baseline << "\n"
            << "solver_converged "
            << regularized.face_pairing_solver_converged << "\n"
            << "solver_iterations "
            << regularized.face_pairing_iterations << "\n"
            << "unresolved_mode_count "
            << regularized.face_pairing_unresolved_mode_count << "\n"
            << "unresolved_mode_l2 "
            << regularized.face_pairing_unresolved_mode_l2 << "\n"
            << "face_residual_before " << baseline_norms.global_l2 << "\n"
            << "face_residual_after " << regularized_norms.global_l2 << "\n"
            << "core_face_residual_before " << baseline_norms.core_l2 << "\n"
            << "core_face_residual_after " << regularized_norms.core_l2 << "\n"
            << "candidate_valid "
            << regularized.face_pairing_candidate_valid << "\n"
            << "candidate_face_residual_after "
            << regularized.face_pairing_candidate_residual_after << "\n"
            << "candidate_core_face_residual_after "
            << regularized.face_pairing_candidate_core_residual_after << "\n"
            << "candidate_delta_ke "
            << regularized.face_pairing_candidate_delta_ke << "\n"
            << "candidate_delta_work "
            << regularized.face_pairing_candidate_delta_work << "\n"
            << "candidate_mass_difference "
            << regularized.face_pairing_candidate_mass_error << "\n"
            << "candidate_min_before_fallback "
            << regularized.face_pairing_candidate_min_before_fallback << "\n"
            << "capacity_active_cells "
            << regularized.face_pairing_capacity_active_cells << "\n"
            << "trust_region_active_cells "
            << regularized.face_pairing_trust_region_active_cells << "\n"
            << "nonzero_capacity_cells "
            << regularized.face_pairing_nonzero_capacity_cells << "\n"
            << "bound_saturated_cells "
            << regularized.face_pairing_bound_saturated_cells << "\n"
            << "correction_l2 "
            << regularized.face_pairing_correction_l2 << "\n"
            << "correction_linf "
            << regularized.face_pairing_correction_linf << "\n"
            << "requested_correction_l2 "
            << regularized.face_pairing_requested_correction_l2 << "\n"
            << "requested_correction_linf "
            << regularized.face_pairing_requested_correction_linf << "\n"
            << "applied_correction_l2 "
            << regularized.face_pairing_applied_correction_l2 << "\n"
            << "applied_correction_linf "
            << regularized.face_pairing_applied_correction_linf << "\n"
            << "objective_residual "
            << regularized.face_pairing_objective_residual << "\n"
            << "objective_smoothness "
            << regularized.face_pairing_objective_smoothness << "\n"
            << "objective_amplitude "
            << regularized.face_pairing_objective_amplitude << "\n"
            << "objective_total "
            << regularized.face_pairing_objective_total << "\n"
            << "rejection_mask "
            << regularized.face_pairing_rejection_mask << "\n"
            << "pass_solver "
            << regularized.face_pairing_pass_solver << "\n"
            << "pass_apply "
            << regularized.face_pairing_pass_apply << "\n"
            << "pass_global_residual "
            << regularized.face_pairing_pass_global_residual << "\n"
            << "pass_core_residual "
            << regularized.face_pairing_pass_core_residual << "\n"
            << "pass_correction_trust "
            << regularized.face_pairing_pass_correction_trust << "\n"
            << "pass_delta_ke "
            << regularized.face_pairing_pass_delta_ke << "\n"
            << "pass_delta_work "
            << regularized.face_pairing_pass_delta_work << "\n"
            << "pass_candidate_min "
            << regularized.face_pairing_pass_candidate_min << "\n"
            << "pass_mass "
            << regularized.face_pairing_pass_mass << "\n"
            << "pass_f_residual "
            << regularized.face_pairing_pass_f_residual << "\n"
            << "candidate_min "
            << regularized.face_pairing_candidate_min << "\n"
            << "mass_difference " << mass_difference << "\n"
            << "candidate_mass_relative_error "
            << regularized.face_pairing_mass_relative_error << "\n"
            << "candidate_cell_mass_error_linf "
            << regularized.face_pairing_cell_mass_error_linf << "\n"
            << "candidate_cell_mass_relative_linf "
            << regularized.face_pairing_cell_mass_relative_linf << "\n"
            << "delta_ke " << delta_ke << "\n"
            << "delta_work " << delta_work << "\n"
            << "energy_pair_error "
            << regularized.face_pairing_energy_pair_error << "\n"
            << "energy_pair_relative "
            << regularized.face_pairing_energy_pair_relative << "\n"
            << "energy_residual_scale "
            << regularized.face_pairing_energy_residual_scale << "\n"
            << "energy_residual_ratio "
            << regularized.face_pairing_energy_residual_ratio << "\n"
            << "correction_trust_limit "
            << regularized.face_pairing_correction_trust_limit << "\n"
            << "correction_trust_ratio "
            << regularized.face_pairing_correction_trust_ratio << "\n"
            << "f_residual_relative_growth "
            << regularized.face_pairing_f_residual_relative_growth << "\n"
            << "JN_linf_difference " << jn_difference << "\n"
            << "Ex_face_linf_difference " << field_difference << "\n"
            << "baseline_R_FV " << baseline.stage5_r_fv << "\n"
            << "regularized_R_FV " << regularized.stage5_r_fv << "\n"
            << "baseline_R_couple " << baseline.stage5_r_couple << "\n"
            << "regularized_R_couple " << regularized.stage5_r_couple << "\n"
            << "pass_energy_pair "
            << regularized.face_pairing_pass_energy_pair << "\n"
            << "pass_energy_residual_scale "
            << regularized.face_pairing_pass_energy_residual_scale << "\n"
            << "pass_mass " << regularized.face_pairing_pass_mass << "\n"
            << "pass_f_residual "
            << regularized.face_pairing_pass_f_residual << "\n"
            << "passes " << global_pass << "\n";

        std::ofstream faces(
            (options.output_dir +
             "/face_pairing_checkpoint_faces.dat").c_str());
        faces << "# global_face x_face JN GstarJE_baseline "
                 "GstarJE_regularized residual_before residual_after\n"
              << std::scientific << std::setprecision(16);
        const size_t rows = std::min(
            baseline_jn.size(),
            std::min(baseline_gstar.size(), regularized_gstar.size()));
        for (size_t iface = 0; iface < rows; ++iface) {
            faces << iface << " " << iface * grid.dx << " "
                  << baseline_jn[iface] << " "
                  << baseline_gstar[iface] << " "
                  << regularized_gstar[iface] << " "
                  << baseline_jn[iface] - baseline_gstar[iface] << " "
                  << baseline_jn[iface] - regularized_gstar[iface] << "\n";
        }
    }

    MPI_Finalize();
    return global_pass ? 0 : 1;
}
