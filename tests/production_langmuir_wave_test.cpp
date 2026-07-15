#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#endif

#ifndef FP_TEST_RUN_ID
#define FP_TEST_RUN_ID "unidentified"
#endif

#ifndef FP_LANGMUIR_RESULT_DIR
#define FP_LANGMUIR_RESULT_DIR "output"
#endif

namespace {

const double kEpsilonDefault = 1.0e-4;
const int kModeDefault = 8;

struct ComplexMode {
    double re;
    double im;
};

struct Options {
    bool fct_enabled;
    double epsilon;
    int mode;
    double periods;
    double dt_scale;
    std::string series_path;
    std::string result_path;
};

bool create_directory_if_needed(const std::string& path)
{
    if (path.empty() || path == ".") return true;
#if defined(_WIN32)
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

bool ensure_parent_directory(const std::string& path)
{
    const std::string::size_type slash = path.find_last_of("/\\");
    if (slash == std::string::npos) return true;
    const std::string parent = path.substr(0, slash);
    const std::string::size_type parent_slash = parent.find_last_of("/\\");
    if (parent_slash != std::string::npos &&
        !ensure_parent_directory(parent))
        return false;
    return create_directory_if_needed(parent);
}

bool parse_double(const char* text, double& value)
{
    char* end = 0;
    value = std::strtod(text, &end);
    return end != text && *end == '\0' && std::isfinite(value);
}

bool parse_int(const char* text, int& value)
{
    char* end = 0;
    const long parsed = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || parsed < 1 ||
        parsed > std::numeric_limits<int>::max())
        return false;
    value = static_cast<int>(parsed);
    return true;
}

bool parse_options(int argc, char** argv, Options& options)
{
    options.fct_enabled = false;
    options.epsilon = kEpsilonDefault;
    options.mode = kModeDefault;
    options.periods = 0.5;
    options.dt_scale = 1.0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--transport" && i + 1 < argc) {
            const std::string mode(argv[++i]);
            if (mode == "high-no-fct") options.fct_enabled = false;
            else if (mode == "high-fct") options.fct_enabled = true;
            else return false;
        } else if (arg == "--run" && i + 1 < argc) {
            const std::string run(argv[++i]);
            if (run == "smoke") options.periods = 0.5;
            else if (run == "short") options.periods = 3.0;
            else if (run == "formal") options.periods = 8.0;
            else return false;
        } else if (arg == "--periods" && i + 1 < argc) {
            if (!parse_double(argv[++i], options.periods) || options.periods <= 0.0)
                return false;
        } else if (arg == "--epsilon" && i + 1 < argc) {
            if (!parse_double(argv[++i], options.epsilon) ||
                options.epsilon <= 0.0 || options.epsilon >= 1.0)
                return false;
        } else if (arg == "--mode" && i + 1 < argc) {
            if (!parse_int(argv[++i], options.mode) || options.mode >= Param::nx / 2)
                return false;
        } else if (arg == "--dt-scale" && i + 1 < argc) {
            if (!parse_double(argv[++i], options.dt_scale) || options.dt_scale <= 0.0)
                return false;
        } else if (arg == "--series" && i + 1 < argc) {
            options.series_path = argv[++i];
        } else if (arg == "--result" && i + 1 < argc) {
            options.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

double global_face_energy(const EMFields& fields, const SpatialGrid& sg)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double e = fields.Ex_face[static_cast<size_t>(iface)];
        local += 0.5 * Const::eps0 * e * e * sg.dx;
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

void global_background_moments(const Species& species,
                               double& number, double& kinetic_energy)
{
    double local_number = 0.0;
    double local_energy = 0.0;
    species.total_particle_number_and_energy(local_number, local_energy);
    double local[2] = {local_number, local_energy};
    double global[2] = {0.0, 0.0};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    number = global[0];
    kinetic_energy = global[1];
}

ComplexMode global_cell_mode(const std::vector<double>& values,
                             const SpatialGrid& sg, int mode)
{
    const int ng = sg.nghost;
    const double factor = 2.0 * Const::pi * static_cast<double>(mode) / Param::Lx;
    double local[2] = {0.0, 0.0};
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double phase = factor * sg.x(ng + ix);
        const double value = values[static_cast<size_t>(ix)];
        local[0] += value * std::cos(phase);
        local[1] -= value * std::sin(phase);
    }
    double global[2] = {0.0, 0.0};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ComplexMode result = {global[0] / Param::nx, global[1] / Param::nx};
    return result;
}

ComplexMode global_face_mode(const std::vector<double>& values,
                             const SpatialGrid& sg, int mode)
{
    const double factor = 2.0 * Const::pi * static_cast<double>(mode) / Param::Lx;
    double local[2] = {0.0, 0.0};
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double x = sg.x_min + (sg.ix_start + iface) * sg.dx;
        const double phase = factor * x;
        const double value = values[static_cast<size_t>(iface)];
        local[0] += value * std::cos(phase);
        local[1] -= value * std::sin(phase);
    }
    double global[2] = {0.0, 0.0};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ComplexMode result = {global[0] / Param::nx, global[1] / Param::nx};
    return result;
}

double mode_abs(const ComplexMode& value)
{
    return std::sqrt(value.re * value.re + value.im * value.im);
}

double global_linf_difference(const std::vector<double>& lhs,
                              const std::vector<double>& rhs,
                              const SpatialGrid& sg)
{
    double local = 0.0;
    const int nface = std::min(sg.nx_local,
        static_cast<int>(std::min(lhs.size(), rhs.size())));
    for (int iface = 0; iface < nface; ++iface)
        local = std::max(local, std::fabs(lhs[static_cast<size_t>(iface)] -
                                         rhs[static_cast<size_t>(iface)]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double global_linf(const std::vector<double>& values, const SpatialGrid& sg)
{
    double local = 0.0;
    const int nface = std::min(sg.nx_local, static_cast<int>(values.size()));
    for (int iface = 0; iface < nface; ++iface)
        local = std::max(local, std::fabs(values[static_cast<size_t>(iface)]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double global_distribution_min(const Species& species, const SpatialGrid& sg)
{
    double local = std::numeric_limits<double>::infinity();
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix)
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                local = std::min(local, species.f[idx3(ng + ix, j, k)]);
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global;
}

bool all_finite_state(const Species& species, const EMFields& fields,
                      const std::vector<double>& jn,
                      const std::vector<double>& je,
                      const SpatialGrid& sg)
{
    int local = 1;
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local && local; ++ix)
        for (int j = 0; j < Param::Nv && local; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                if (!std::isfinite(species.f[idx3(ng + ix, j, k)])) {
                    local = 0;
                    break;
                }
    for (int iface = 0; iface < sg.nx_local && local; ++iface) {
        const size_t slot = static_cast<size_t>(iface);
        if (!std::isfinite(fields.Ex_face[slot]) ||
            (slot < jn.size() && !std::isfinite(jn[slot])) ||
            (slot < je.size() && !std::isfinite(je[slot])))
            local = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return local != 0;
}

void initialize_langmuir_state(Species& background, EMFields& fields,
                               const SpatialGrid& sg, const Options& options,
                               int rank, int size)
{
    background.initialize_maxwellian();
    const int ng = sg.nghost;
    const double kx = 2.0 * Const::pi * options.mode / Param::Lx;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double profile = 1.0 + options.epsilon *
            std::cos(kx * sg.x(ng + ix));
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                background.f[idx3(ng + ix, j, k)] *= profile;
    }
    background.compute_moments();
    fields.init(sg);
    const std::vector<double> beam_density(static_cast<size_t>(sg.nx_local), 0.0);
    const std::vector<double> ion_density(static_cast<size_t>(sg.nx_local),
                                          Param::dens);
    fields.set_charge_density(background, beam_density, ion_density);
    // This is the one permitted Poisson operation: it builds a periodic,
    // discrete-Gauss-consistent initial face field.  All later steps use
    // only the production Yee/Ampere update.
    fields.solve_poisson(rank, size);
    fields.update_gauss_residual_diagnostics(rank, size);
}

double choose_dt(const Species& species, const EMFields& fields,
                 const SpatialGrid& sg, double scale)
{
    double local_max_e = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface)
        local_max_e = std::max(local_max_e,
            std::fabs(fields.Ex_face[static_cast<size_t>(iface)]));
    double max_e = 0.0;
    MPI_Allreduce(&local_max_e, &max_e, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const double min_du = *std::min_element(species.cgrid.upar_widths.begin(),
                                            species.cgrid.upar_widths.end());
    const double acceleration = Const::qe * max_e / (Const::me * Const::c);
    const double dt_x = 0.15 * sg.dx / Const::c;
    const double dt_u = acceleration > 0.0
        ? 0.15 * min_du / acceleration : dt_x;
    return scale * std::min(dt_x, dt_u);
}

void write_record(std::ofstream& out, double time, int step_index,
                  const ComplexMode& e_mode, const ComplexMode& n_mode,
                  const ComplexMode& jn_mode, const ComplexMode& je_mode,
                  double number_relative, double kinetic_energy,
                  double field_energy, double initial_total_energy,
                   const EMFields& fields,
                   const VlasovAmpereMidpointSolver::Result* step,
                   double high_candidate_min, double final_min,
                   double j_pair_l2, double j_pair_linf, double beam_current_linf,
                   int finite,
                  int strict_accept)
{
    const double total_energy = kinetic_energy + field_energy;
    const double energy_drift = (total_energy - initial_total_energy) /
        std::max(std::numeric_limits<double>::min(),
                 std::fabs(initial_total_energy));
    const int accepted = 1;
    out << std::scientific << std::setprecision(17)
        << time / Const::femto << " " << step_index << " "
        << accepted << " " << strict_accept << " "
        << e_mode.re << " " << e_mode.im << " " << mode_abs(e_mode) << " "
        << n_mode.re << " " << n_mode.im << " " << mode_abs(n_mode) << " "
        << jn_mode.re << " " << jn_mode.im << " " << mode_abs(jn_mode) << " "
        << je_mode.re << " " << je_mode.im << " " << mode_abs(je_mode) << " "
        << number_relative << " " << kinetic_energy << " " << field_energy
        << " " << total_energy << " " << energy_drift << " "
        << fields.last_gauss_residual_l1 << " " << fields.last_gauss_residual_linf
        << " " << (step ? step->stage5_r_fv : 0.0)
        << " " << (step ? step->stage5_r_couple : 0.0)
        << " " << (step ? step->stage5_r_couple_centered : 0.0)
        << " " << (step ? step->stage5_r_couple_upwind_stabilization : 0.0)
        << " " << (step ? step->stage5_r_couple_fct_stabilization : 0.0)
         << " " << j_pair_l2
         << " " << j_pair_linf
         << " " << beam_current_linf
        << " " << (step ? step->nonlinear_iterations : 0)
        << " " << (step ? step->residual_E : 0.0)
         << " " << (step ? step->residual_J_bkg : 0.0)
         << " " << (step ? step->residual_f : 0.0)
         << " " << high_candidate_min << " " << final_min
         << " " << (step ? step->low_order_candidate_min : 0.0)
         << " " << (step ? step->low_order_negative_count : 0)
         << " " << (step ? step->low_order_negative_mass : 0.0)
         << " " << (step ? step->low_order_roundoff_zeroed_count : 0)
         << " " << (step ? step->low_order_roundoff_zeroed_mass : 0.0)
         << " " << (step ? step->limiter_active_fraction : 0.0)
        << " " << (step ? step->limiter_active_fraction_core : 0.0)
        << " " << (step ? step->limiter_active_fraction_boundary : 0.0)
        << " " << (step ? step->limiter_min_alpha : 1.0)
        << " " << (step ? step->limiter_min_alpha_core : 1.0)
        << " " << (step ? step->limiter_min_alpha_boundary : 1.0)
        << " " << (step ? step->u_boundary_particle_outflow : 0.0)
        << " " << (step ? step->u_boundary_energy_outflow : 0.0)
        << " " << (step ? step->neg_mass_core : 0.0)
         << " " << (step ? step->neg_mass_boundary : 0.0)
         << " " << (step ? step->neg_current_core_abs : 0.0)
         << " " << (step ? step->neg_energy_core_abs : 0.0)
         << " " << (step ? step->neg_mass_total_guard : 0.0)
         << " " << (step ? step->negative_debt_level : 0)
         << " " << (step ? step->accepted_with_negative_debt : 0)
         << " " << finite
        << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Options options;
    if (!parse_options(argc, argv, options)) {
        if (rank == 0) {
            std::cerr << "usage: production_langmuir_wave_test "
                      << "[--transport high-no-fct|high-fct] "
                      << "[--run smoke|short|formal] [--periods P] "
                      << "[--epsilon E] [--mode M] [--dt-scale S] "
                      << "[--series path] [--result path]\n";
        }
        MPI_Finalize();
        return 2;
    }

    SpatialGrid sg;
    sg.init(rank, size);
    const std::string transport = options.fct_enabled ? "high_fct" : "high_no_fct";
    if (options.series_path.empty()) {
        std::ostringstream path;
        path << FP_LANGMUIR_RESULT_DIR << "/langmuir_m" << options.mode
             << "_eps" << std::scientific << std::setprecision(0)
             << options.epsilon << "_" << transport << ".dat";
        options.series_path = path.str();
    }
    if (options.result_path.empty())
        options.result_path = options.series_path + ".result";

    Species background;
    background.init("langmuir_background", SpeciesType::BACKGROUND_ELECTRON,
                    -Const::qe, Const::me, Param::dens, Param::temperature_e,
                    false, sg);
    EMFields fields;
    initialize_langmuir_state(background, fields, sg, options, rank, size);
    const double initial_gauss_linf = fields.last_gauss_residual_linf;
    double local_gauss_scale = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix)
        local_gauss_scale = std::max(local_gauss_scale,
            std::fabs(fields.rho[sg.nghost + ix]) / Const::eps0);
    double gauss_scale = 0.0;
    MPI_Allreduce(&local_gauss_scale, &gauss_scale, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double initial_gauss_relative = initial_gauss_linf /
        std::max(std::numeric_limits<double>::min(), gauss_scale);
    const double dt = choose_dt(background, fields, sg, options.dt_scale);
    const double plasma_period = 2.0 * Const::pi / Param::omega_pe;
    const double duration = options.periods * plasma_period;

    double initial_number = 0.0;
    double initial_ke = 0.0;
    global_background_moments(background, initial_number, initial_ke);
    const double initial_field_energy = global_face_energy(fields, sg);
    const double initial_total_energy = initial_ke + initial_field_energy;
    background.compute_moments();

    std::ofstream series;
    bool series_ok = true;
    if (rank == 0) {
        series_ok = ensure_parent_directory(options.series_path);
        if (series_ok) series.open(options.series_path.c_str());
        series_ok = series_ok && static_cast<bool>(series);
        if (series_ok) {
            series << "# time_fs step accepted strict_accept "
                   << "E_re E_im E_abs n_re n_im n_abs JN_re JN_im JN_abs "
                   << "JE_re JE_im JE_abs number_relative KE UE total_energy "
                    << "energy_drift gauss_l1 gauss_linf R_FV R_couple R_couple_centered "
                    << "R_couple_upwind_stabilization R_couple_fct_stabilization "
                    << "JN_minus_GstarJE_l2 JN_minus_GstarJE_linf beam_current_linf "
                    << "nonlinear_iterations residual_E "
                    << "residual_J residual_f high_candidate_min final_min "
                    << "low_order_min low_order_negative_count low_order_negative_mass "
                    << "low_order_roundoff_zeroed_count low_order_roundoff_zeroed_mass "
                   << "limiter_active limiter_active_core limiter_active_boundary "
                   << "min_alpha min_alpha_core min_alpha_boundary u_boundary_particles "
                    << "u_boundary_energy neg_mass_core neg_mass_boundary "
                    << "neg_current_core neg_energy_core neg_mass_total "
                    << "negative_debt_level finite_negative_debt_accepted finite\n";
        }
    }
    int series_ok_int = series_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &series_ok_int, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (series_ok_int == 0 || !(dt > 0.0) || !std::isfinite(dt)) {
        if (rank == 0) std::cerr << "Langmuir test initialization failed\n";
        MPI_Finalize();
        return 1;
    }

    const ComplexMode e0 = global_face_mode(fields.Ex_face, sg, options.mode);
    const ComplexMode n0 = global_cell_mode(background.number_density, sg, options.mode);
    const double initial_distribution_min = global_distribution_min(background, sg);
    if (rank == 0)
        write_record(series, 0.0, 0, e0, n0, ComplexMode{0.0, 0.0},
                     ComplexMode{0.0, 0.0}, 0.0, initial_ke,
                      initial_field_energy, initial_total_energy, fields, 0,
                      initial_distribution_min, initial_distribution_min,
                      0.0, 0.0, 0.0, 1, 1);

    VlasovAmpereMidpointSolver solver;
    solver.set_step_diagnostics_enabled(false);
    solver.set_beam_enabled(false);
    solver.set_low_order_only(false);
    solver.set_nonuniform_high_order_enabled(true);
    solver.set_fct_enabled(options.fct_enabled);
    solver.set_fct_activation_audit_enabled(options.fct_enabled);
    solver.set_allow_finite_negative_debt_for_test(!options.fct_enabled);
    solver.set_max_midpoint_iterations(40);
    solver.set_midpoint_iteration_trace_for_test(true);

    BeamPIC beam;
    double time = 0.0;
    int step_index = 0;
    int accepted_steps = 0;
    bool failed = false;
    std::string failure_reason = "none";
    VlasovAmpereMidpointSolver::Result failure_step;
    bool have_failure_step = false;
    int soft_unconverged_observed = 0;
    double max_rfv = 0.0;
    double max_rcouple = 0.0;
    double max_limiter_core = 0.0;
    double max_beam_current = 0.0;

    while (time + 0.5 * dt < duration) {
        const VlasovAmpereMidpointSolver::Result step =
            solver.advance_background_and_fields(background, beam, fields, sg,
                                                 dt, time, rank, size);
        ++step_index;
        const bool strict_local = step.state_advanced && step.converged &&
            !step.soft_unconverged && !step.failed;
        int strict_int = strict_local ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &strict_int, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if (strict_int == 0) {
            if (step.soft_unconverged) ++soft_unconverged_observed;
            failed = true;
            failure_step = step;
            have_failure_step = true;
            std::ostringstream reason;
            reason << "state_advanced=" << step.state_advanced
                   << ",converged=" << step.converged
                   << ",soft_unconverged=" << step.soft_unconverged
                   << ",failed=" << step.failed
                   << ",failure_reason=" << step.failure_reason;
            failure_reason = reason.str();
            break;
        }

        background = step.species_np1;
        fields = step.fields_np1;
        time += dt;
        ++accepted_steps;
        background.compute_moments();
        const std::vector<double> beam_density(static_cast<size_t>(sg.nx_local), 0.0);
        const std::vector<double> ion_density(static_cast<size_t>(sg.nx_local),
                                              Param::dens);
        fields.set_charge_density(background, beam_density, ion_density);
        fields.update_gauss_residual_diagnostics(rank, size);

        double number = 0.0;
        double kinetic_energy = 0.0;
        global_background_moments(background, number, kinetic_energy);
        const double field_energy = global_face_energy(fields, sg);
        const ComplexMode e_mode = global_face_mode(fields.Ex_face, sg, options.mode);
        const ComplexMode n_mode = global_cell_mode(background.number_density, sg, options.mode);
        const ComplexMode jn_mode = global_face_mode(step.j_bkg_face_mid, sg, options.mode);
        const ComplexMode je_mode = global_face_mode(step.j_bkg_energy_debug_face, sg, options.mode);
        const double final_min = global_distribution_min(background, sg);
        const double high_min = std::isfinite(step.fct_high_candidate_min)
            ? step.fct_high_candidate_min : final_min;
        const double j_pair_linf = global_linf_difference(
            step.j_bkg_face_mid, step.j_bkg_energy_debug_face, sg);
        const double beam_current_linf = global_linf(step.j_beam_face_mid, sg);
        const bool finite = all_finite_state(background, fields,
                                             step.j_bkg_face_mid,
                                             step.j_bkg_energy_debug_face, sg);
        if (!finite) {
            failed = true;
            failure_step = step;
            have_failure_step = true;
            failure_reason = "nonfinite accepted state";
            break;
        }
        max_rfv = std::max(max_rfv, std::fabs(step.stage5_r_fv));
        max_rcouple = std::max(max_rcouple, std::fabs(step.stage5_r_couple));
        max_limiter_core = std::max(max_limiter_core,
                                    step.limiter_active_fraction_core);
        max_beam_current = std::max(max_beam_current, beam_current_linf);
        if (rank == 0) {
            write_record(series, time, step_index, e_mode, n_mode, jn_mode,
                         je_mode, (number - initial_number) /
                         std::max(std::numeric_limits<double>::min(),
                                  std::fabs(initial_number)), kinetic_energy,
                          field_energy, initial_total_energy, fields, &step,
                          high_min, final_min, step.stage5_jn_minus_je_l2,
                          j_pair_linf, beam_current_linf, finite ? 1 : 0, 1);
        }
    }

    if (rank == 0 && series) series.close();
    int success = (!failed && accepted_steps > 0 &&
                   initial_gauss_relative <= 1.0e-10) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        const std::string trace_path =
            options.result_path + ".failure_midpoint_residuals.dat";
        if (have_failure_step &&
            !failure_step.midpoint_residual_e_history.empty() &&
            ensure_parent_directory(trace_path)) {
            std::ofstream trace(trace_path.c_str());
            if (trace) {
                trace << "# iteration residual_E\n";
                trace << std::scientific << std::setprecision(17);
                for (size_t i = 0;
                     i < failure_step.midpoint_residual_e_history.size(); ++i)
                    trace << (i + 1) << " "
                          << failure_step.midpoint_residual_e_history[i] << "\n";
            }
        }
        std::ofstream result;
        if (ensure_parent_directory(options.result_path))
            result.open(options.result_path.c_str());
        if (result) {
            result << std::scientific << std::setprecision(17)
                   << "format_version=1\n"
                   << "test=production_langmuir_wave\n"
                   << "test_run_id=" << FP_TEST_RUN_ID << "\n"
                   << "transport=" << transport << "\n"
                   << "beam_enabled=0\ncollision_enabled=0\n"
                   << "epsilon=" << options.epsilon << "\n"
                   << "mode=" << options.mode << "\n"
                   << "plasma_period=" << plasma_period << "\n"
                   << "target_periods=" << options.periods << "\n"
                    << "dt=" << dt << "\n"
                    << "accepted_steps=" << accepted_steps << "\n"
                    << "state_advanced=" << success << "\n"
                    << "strict_accept=" << success << "\n"
                    << "converged=" << success << "\n"
                    << "failed=" << (failed ? 1 : 0) << "\n"
                    << "initial_gauss_linf=" << initial_gauss_linf << "\n"
                   << "initial_gauss_relative=" << initial_gauss_relative << "\n"
                   << "max_abs_R_FV=" << max_rfv << "\n"
                    << "max_abs_R_couple=" << max_rcouple << "\n"
                    << "max_limiter_active_core=" << max_limiter_core << "\n"
                    << "max_abs_beam_current=" << max_beam_current << "\n"
                    << "midpoint_max_iterations=40\n"
                    << "soft_unconverged_observed=" << soft_unconverged_observed << "\n"
                    << "failure_midpoint_residual_trace=" << trace_path << "\n"
                    << "failure_midpoint_residual_count=" << (have_failure_step
                        ? failure_step.midpoint_residual_e_history.size() : 0) << "\n"
                    << "failure_final_min=" << (have_failure_step
                        ? failure_step.failure_final_min : 0.0) << "\n"
                    << "failure_worst_ix=" << (have_failure_step
                        ? failure_step.failure_worst_ix : -1) << "\n"
                    << "failure_worst_iv=" << (have_failure_step
                        ? failure_step.failure_worst_iv : -1) << "\n"
                    << "failure_worst_imu=" << (have_failure_step
                        ? failure_step.failure_worst_imu : -1) << "\n"
                    << "failure_worst_is_core=" << (have_failure_step
                        ? failure_step.failure_worst_is_core : -1) << "\n"
                    << "failure_worst_is_tail=" << (have_failure_step
                        ? failure_step.failure_worst_is_tail : -1) << "\n"
                    << "failure_neg_mass_total=" << (have_failure_step
                        ? failure_step.neg_mass_total_guard : 0.0) << "\n"
                    << "failure_neg_mass_core=" << (have_failure_step
                        ? failure_step.neg_mass_core : 0.0) << "\n"
                    << "failure_neg_mass_core_fraction=" << (have_failure_step
                        ? failure_step.neg_mass_core_fraction : 0.0) << "\n"
                    << "failure_neg_energy_core_fraction=" << (have_failure_step
                        ? failure_step.neg_energy_core_fraction : 0.0) << "\n"
                    << "failure_neg_current_core_fraction=" << (have_failure_step
                        ? failure_step.neg_current_core_fraction : 0.0) << "\n"
                    << "failure_negative_debt_level=" << (have_failure_step
                        ? failure_step.negative_debt_level : 0) << "\n"
                    << "failure_low_order_min=" << (have_failure_step
                        ? failure_step.low_order_candidate_min : 0.0) << "\n"
                    << "failure_low_order_negative_count=" << (have_failure_step
                        ? failure_step.low_order_negative_count : 0) << "\n"
                    << "failure_low_order_negative_mass=" << (have_failure_step
                        ? failure_step.low_order_negative_mass : 0.0) << "\n"
                    << "failure_low_order_roundoff_zeroed_count=" <<
                        (have_failure_step ? failure_step.low_order_roundoff_zeroed_count
                                           : 0) << "\n"
                    << "failure_low_order_roundoff_zeroed_mass=" <<
                        (have_failure_step ? failure_step.low_order_roundoff_zeroed_mass
                                           : 0.0) << "\n"
                    << "series_file=" << options.series_path << "\n"
                   << "failure_reason=" << failure_reason << "\n";
        }
        std::cout << "production_langmuir_wave_test transport=" << transport
                  << " accepted_steps=" << accepted_steps
                  << " success=" << success
                  << " series=" << options.series_path << "\n";
    }
    MPI_Finalize();
    return success ? 0 : 1;
}
