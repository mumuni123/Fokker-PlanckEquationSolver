#include "beam_pic.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

#ifndef FP_AMPERE_TEST_SUITE_FINEST_NV
#define FP_AMPERE_TEST_SUITE_FINEST_NV 192
#endif

#ifndef FP_TEST_RUN_ID
#define FP_TEST_RUN_ID "unidentified"
#endif

#ifndef FP_TEST_RESULT_DIR
#define FP_TEST_RESULT_DIR "output"
#endif

#ifndef FP_FCT_ACTIVATION_TEST
#define FP_FCT_ACTIVATION_TEST 0
#endif

const double kEnergyScaleFloor = 64.0 * std::numeric_limits<double>::min();
const double kTrendSlack = 1.05;
const double kDriftVx = 0.05 * Const::c;
const double kInitialEx = 2.0e3;

struct TransportMode {
    const char* name;
    bool low_order_only;
    bool fct_enabled;
    bool require_fct_activation;
};

const TransportMode kLowOrderMode = {"low_order", true, false, false};
const TransportMode kHighNoFctMode = {"high_no_fct", false, false, false};
const TransportMode kHighFctMode = {"high_fct", false, true, false};
const TransportMode kHighFctActivationMode = {
    "high_fct_activation", false, true, true};

const TransportMode* parse_transport_mode(const std::string& value)
{
    if (value == "low-order" || value == "low_order") return &kLowOrderMode;
    if (value == "high-no-fct" || value == "high_no_fct") return &kHighNoFctMode;
    if (value == "high-fct" || value == "high_fct") return &kHighFctMode;
    if (value == "high-fct-activation" || value == "high_fct_activation")
        return &kHighFctActivationMode;
    return 0;
}

struct InitializationAudit {
    double number_density;
    double mean_vx;
    double current_density;
    double target_mean_vx;
    double target_current_density;
    double mean_vx_relative_error;
    double current_relative_error;
    double min_upar_width;
    double min_uperp_width;
    double upar_width_at_drift;
    double cells_per_uth_upar_at_drift;
    double cells_per_uth_uperp_at_axis;
    int occupied_upar_cells;
    int occupied_uperp_cells;
};

struct ContinuousInitializationTarget {
    double mean_vx;
    double current_density;
};

struct Audit {
    double dt;
    double delta_k;
    double delta_k_subtractive;
    double delta_k_method_difference;
    double delta_u;
    double max_abs_jn;
    double max_abs_gstar_je;
    double e_dot_jn;
    double e_dot_gstar_je;
    double energy_exchange_scale;
    double energy_residual_abs;
    double energy_relative;
    double j_pair_relative;
    double number_relative;
    double beam_current_linf;
    double limiter_active_fraction;
    double limiter_min_alpha;
    double limiter_active_fraction_core;
    double limiter_active_fraction_boundary;
    double limiter_min_alpha_core;
    double limiter_min_alpha_boundary;
    double stage5_r_fv;
    double stage5_r_couple;
    double stage5_r_couple_centered;
    double stage5_r_couple_upwind_stabilization;
    double stage5_r_couple_fct_stabilization;
    double stage5_r_fv_relative;
    double fct_high_candidate_min;
    double fct_high_candidate_donor_excess;
    long long fct_controlled_injection_count;
    int state_advanced;
    int converged;
    int soft_unconverged;
    int failed;
    int low_order_only;
    int fct_enabled;
    int finite;
    InitializationAudit initial;
};

double gauss_legendre_node(int order, int index, double& weight)
{
    const double pi = std::acos(-1.0);
    double z = std::cos(pi * (static_cast<double>(index) + 0.75) /
                        (static_cast<double>(order) + 0.5));
    double z_previous = 0.0;
    double derivative = 0.0;
    do {
        double p0 = 1.0;
        double p1 = z;
        for (int n = 2; n <= order; ++n) {
            const double p2 = ((2.0 * n - 1.0) * z * p1 -
                               (n - 1.0) * p0) / n;
            p0 = p1;
            p1 = p2;
        }
        derivative = static_cast<double>(order) *
            (z * p1 - p0) / (z * z - 1.0);
        z_previous = z;
        z -= p1 / derivative;
    } while (std::fabs(z - z_previous) > 8.0 *
             std::numeric_limits<double>::epsilon());
    weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
    return z;
}

ContinuousInitializationTarget continuous_initialization_target(
    const Species& species)
{
    // This is deliberately independent of the velocity-cell quadrature used
    // by initialize_maxwellian().  The test compares the production grid
    // quadrature against the same continuous, truncated drifting Maxwellian.
    const int order = 96;
    const double uth = std::sqrt(species.temperature /
        (species.mass * Const::c * Const::c));
    const double drift_u = u_from_v(kDriftVx);
    const double extent = 12.0 * uth;
    const double upar_lo = std::max(-Param::momentum_umax, drift_u - extent);
    const double upar_hi = std::min( Param::momentum_umax, drift_u + extent);
    const double uperp_hi = std::min(Param::momentum_umax, extent);
    const double inv2uth2 = 1.0 / (2.0 * uth * uth);
    long double number = 0.0L;
    long double vx_moment = 0.0L;

    for (int a = 0; a < order; ++a) {
        double wa = 0.0;
        const double za = gauss_legendre_node(order, a, wa);
        const double upar = 0.5 * ((upar_hi - upar_lo) * za +
                                   (upar_hi + upar_lo));
        const double dupar_weight = 0.5 * (upar_hi - upar_lo) * wa;
        for (int b = 0; b < order; ++b) {
            double wb = 0.0;
            const double zb = gauss_legendre_node(order, b, wb);
            const double uperp = 0.5 * uperp_hi * (zb + 1.0);
            const double duperp_weight = 0.5 * uperp_hi * wb;
            const long double phase_weight = static_cast<long double>(
                2.0 * Const::pi * uperp * dupar_weight * duperp_weight);
            const long double raw = std::exp(-((upar - drift_u) *
                (upar - drift_u) + uperp * uperp) * inv2uth2);
            number += raw * phase_weight;
            const double gamma = std::sqrt(1.0 + upar * upar +
                                           uperp * uperp);
            vx_moment += raw * phase_weight * (Const::c * upar / gamma);
        }
    }

    ContinuousInitializationTarget target;
    target.mean_vx = static_cast<double>(vx_moment / number);
    target.current_density = species.charge * species.density0 * target.mean_vx;
    return target;
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

double global_background_delta_kinetic_energy(const Species& before,
                                               const Species& after,
                                               const SpatialGrid& sg)
{
    long double local = 0.0L;
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int cell_x = ng + ix;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t cell = idx3(cell_x, j, k);
                local += static_cast<long double>(
                    before.cgrid.kinetic_energy[idx2(j, k)]) *
                    (static_cast<long double>(after.f[cell]) -
                     static_cast<long double>(before.f[cell]));
            }
        }
    }
    long double global = 0.0L;
    MPI_Allreduce(&local, &global, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return static_cast<double>(global);
}

InitializationAudit audit_initial_distribution(const Species& species,
                                                const SpatialGrid& sg)
{
    long double local_moments[2] = {0.0L, 0.0L};
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int cell_x = ng + ix;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const long double mass = species.f[idx3(cell_x, j, k)];
                local_moments[0] += mass;
                local_moments[1] += mass * species.cgrid.vx[idx2(j, k)];
            }
        }
    }
    long double global_moments[2] = {0.0L, 0.0L};
    MPI_Allreduce(local_moments, global_moments, 2, MPI_LONG_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    std::vector<double> upar_marginal(Param::Nv, 0.0);
    std::vector<double> uperp_marginal(Param::Nmu, 0.0);
    if (sg.ix_start == 0 && sg.nx_local > 0) {
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const double mass = species.f[idx3(ng, j, k)];
                upar_marginal[j] += mass;
                uperp_marginal[k] += mass;
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, upar_marginal.data(), Param::Nv, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, uperp_marginal.data(), Param::Nmu, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double cell_mass = 0.0;
    for (int j = 0; j < Param::Nv; ++j) cell_mass += upar_marginal[j];
    const double occupied_threshold = 1.0e-12 * cell_mass;
    int occupied_upar = 0;
    int occupied_uperp = 0;
    for (int j = 0; j < Param::Nv; ++j)
        occupied_upar += upar_marginal[j] > occupied_threshold ? 1 : 0;
    for (int k = 0; k < Param::Nmu; ++k)
        occupied_uperp += uperp_marginal[k] > occupied_threshold ? 1 : 0;

    const double drift_u = u_from_v(kDriftVx);
    const double uth = std::sqrt(species.temperature /
        (species.mass * Const::c * Const::c));
    double width_at_drift = std::numeric_limits<double>::quiet_NaN();
    int cells_upar = 0;
    int cells_uperp = 0;
    const double upar_lo = drift_u - 0.5 * uth;
    const double upar_hi = drift_u + 0.5 * uth;
    for (int j = 0; j < Param::Nv; ++j) {
        if (species.cgrid.upar_faces[j] <= drift_u &&
            drift_u <= species.cgrid.upar_faces[j + 1])
            width_at_drift = species.cgrid.upar_widths[j];
        if (species.cgrid.upar_faces[j + 1] > upar_lo &&
            species.cgrid.upar_faces[j] < upar_hi) ++cells_upar;
    }
    for (int k = 0; k < Param::Nmu; ++k)
        if (species.cgrid.uperp_faces[k] < uth &&
            species.cgrid.uperp_faces[k + 1] > 0.0) ++cells_uperp;

    const double domain_length = sg.nx_global * sg.dx;
    const double total_mass = static_cast<double>(global_moments[0]);
    const double velocity_moment = static_cast<double>(global_moments[1]);
    InitializationAudit audit;
    audit.number_density = total_mass / domain_length;
    audit.mean_vx = velocity_moment / std::max(
        std::numeric_limits<double>::min(), total_mass);
    audit.current_density = species.charge * velocity_moment / domain_length;
    const ContinuousInitializationTarget target =
        continuous_initialization_target(species);
    audit.target_mean_vx = target.mean_vx;
    audit.target_current_density = target.current_density;
    audit.mean_vx_relative_error = std::fabs(audit.mean_vx - audit.target_mean_vx) /
        std::fabs(audit.target_mean_vx);
    audit.current_relative_error = std::fabs(
        audit.current_density - audit.target_current_density) /
        std::fabs(audit.target_current_density);
    audit.min_upar_width = *std::min_element(
        species.cgrid.upar_widths.begin(), species.cgrid.upar_widths.end());
    audit.min_uperp_width = *std::min_element(
        species.cgrid.uperp_widths.begin(), species.cgrid.uperp_widths.end());
    audit.upar_width_at_drift = width_at_drift;
    audit.cells_per_uth_upar_at_drift = cells_upar;
    audit.cells_per_uth_uperp_at_axis = cells_uperp;
    audit.occupied_upar_cells = occupied_upar;
    audit.occupied_uperp_cells = occupied_uperp;
    return audit;
}

double common_velocity_suite_dt(const SpatialGrid& sg)
{
    const int finest_nv = FP_AMPERE_TEST_SUITE_FINEST_NV;
    const double min_fine_du = Param::momentum_umax * std::sinh(
        2.0 * Param::momentum_upar_stretch / finest_nv) /
        std::sinh(Param::momentum_upar_stretch);
    const double acceleration = Const::qe * std::fabs(kInitialEx) /
        (Const::me * Const::c);
    const double dt_x = 0.20 * sg.dx / Const::c;
    const double dt_u = 0.20 * min_fine_du / acceleration;
    return std::min(dt_x, dt_u);
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

bool global_all_true(bool local_value)
{
    int value = local_value ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return value != 0;
}

bool create_directory_if_needed(const std::string& path)
{
    if (path.empty() || path == ".") return true;
#if defined(_WIN32)
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

bool ensure_result_parent_directory(const std::string& path)
{
    const std::string::size_type slash = path.find_last_of("/\\");
    if (slash == std::string::npos) return true;
    const std::string parent = path.substr(0, slash);
    const std::string::size_type parent_slash = parent.find_last_of("/\\");
    if (parent_slash != std::string::npos &&
        !create_directory_if_needed(parent.substr(0, parent_slash)))
        return false;
    return create_directory_if_needed(parent);
}

void run_case_checkpoint(const char* phase, double dt, int rank)
{
    // Every rank reaches both barriers in the same order.  This keeps the
    // diagnostic safe for MPI and identifies the last completed phase if a
    // rank faults inside initialization or the production advance.
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cerr << "production_7_2_checkpoint phase=" << phase
                  << " dt=" << std::scientific << std::setprecision(17)
                  << dt << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

double convergence_order(double coarse, double fine)
{
    if (coarse == 0.0 && fine == 0.0)
        return std::numeric_limits<double>::infinity();
    if (!(coarse > 0.0) || !(fine > 0.0))
        return -std::numeric_limits<double>::infinity();
    return std::log(coarse / fine) / std::log(2.0);
}

void apply_controlled_fct_profile(Species& species, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const double pi = std::acos(-1.0);
    const double drift_u = u_from_v(kDriftVx);
    int drift_j = 0;
    for (int j = 1; j < Param::Nv; ++j) {
        if (std::fabs(species.cgrid.upar_cells[j] - drift_u) <
            std::fabs(species.cgrid.upar_cells[drift_j] - drift_u))
            drift_j = j;
    }
    // Resolve the controlled velocity feature over a few *actual* cells at
    // the occupied drift, rather than placing it in an empty high-u tail.
    const double width = 2.5 * species.cgrid.upar_widths[drift_j];
    double number_before = 0.0;
    double ignored_energy = 0.0;
    global_background_moments(species, number_before, ignored_energy);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double phase = 2.0 * pi * sg.x(ng + ix) / Param::Lx;
        // Strictly positive, periodic and deliberately steep.  It is a
        // controlled numerical transport state, not an equilibrium claim.
        const double x_profile =
            1.0 + 0.98 * std::tanh(20.0 * std::sin(phase));
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
            {
                // Couple the sharp periodic x front to a narrow, positive
                // u_parallel band.  The distribution remains physical
                // (strictly positive), while the two-directional high-order
                // correction has a controlled opportunity to overdraw a
                // donor cell and therefore exercise FCT.
                const double du =
                    (species.cgrid.upar_cells[j] - drift_u) / width;
                const double u_profile =
                    0.20 + 1.80 * std::exp(-0.5 * du * du);
                species.f[static_cast<size_t>(ng + ix) * Param::Nvmu +
                          static_cast<size_t>(j) * Param::Nmu + k] *=
                    x_profile * u_profile;
            }
    }
    double number_after = 0.0;
    global_background_moments(species, number_after, ignored_energy);
    const double renormalize = (number_after > 0.0)
        ? number_before / number_after : 1.0;
    for (int ix = 0; ix < sg.nx_local; ++ix)
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                species.f[static_cast<size_t>(ng + ix) * Param::Nvmu +
                          static_cast<size_t>(j) * Param::Nmu + k] *=
                    renormalize;
}

Audit run_case(const SpatialGrid& sg, int rank, int size, double dt,
               const InitializationAudit& shared_initial,
               const TransportMode& transport)
{
    run_case_checkpoint("begin", dt, rank);
    Species background;
    background.init("production_7_2_background", SpeciesType::BACKGROUND_ELECTRON,
                    -Const::qe, Const::me, Param::dens, Param::temperature_e,
                    false, sg);
    run_case_checkpoint("background_initialized", dt, rank);
    background.initialize_maxwellian(kDriftVx);
    if (transport.require_fct_activation)
        apply_controlled_fct_profile(background, sg);
    run_case_checkpoint("background_maxwellian_initialized", dt, rank);
    const InitializationAudit& initial = shared_initial;
    run_case_checkpoint("initial_distribution_audited", dt, rank);

    // Beam is disabled for section 7.2.  Keep the default-constructed object
    // required by the production solver interface, but do not reserve the
    // real injector's macro-particle capacity on every rank.
    BeamPIC beam;
    EMFields fields;
    fields.init(sg);
    run_case_checkpoint("fields_initialized", dt, rank);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), kInitialEx);
    fields.sync_cell_ex_from_faces(rank, size);
    run_case_checkpoint("field_faces_synchronized", dt, rank);

    double n0 = 0.0;
    double k0 = 0.0;
    global_background_moments(background, n0, k0);
    const double u0 = global_face_energy(fields, sg);
    run_case_checkpoint("initial_global_moments_complete", dt, rank);

    VlasovAmpereMidpointSolver solver;
    solver.set_step_diagnostics_enabled(false);
    solver.set_beam_enabled(false);
    solver.set_low_order_only(transport.low_order_only);
    solver.set_fct_enabled(transport.fct_enabled);
    solver.set_fct_activation_audit_enabled(transport.require_fct_activation);
    solver.set_controlled_fct_flux_injection_enabled(
        transport.require_fct_activation);
    solver.set_nonuniform_high_order_enabled(!transport.low_order_only);
    run_case_checkpoint("before_production_advance", dt, rank);
    const VlasovAmpereMidpointSolver::Result step =
        solver.advance_background_and_fields(background, beam, fields, sg,
                                             dt, 0.0, rank, size);
    run_case_checkpoint("after_production_advance", dt, rank);

    double n1 = 0.0;
    double k1 = 0.0;
    global_background_moments(step.species_np1, n1, k1);
    const double u1 = global_face_energy(step.fields_np1, sg);
    run_case_checkpoint("final_global_moments_complete", dt, rank);
    const double delta_k_subtractive = k1 - k0;
    const double delta_k = global_background_delta_kinetic_energy(
        background, step.species_np1, sg);
    const double delta_k_method_difference = delta_k - delta_k_subtractive;
    const double delta_u = u1 - u0;
    const double energy_residual = delta_k + delta_u;
    // This is a relative closure defect on the actual one-step exchange
    // scale.  A 1 J/m^2 lower bound would turn small-step relative defects
    // into misleading absolute values.
    const double energy_scale = std::max(kEnergyScaleFloor,
        std::max(std::fabs(delta_k), std::fabs(delta_u)));
    const double jn_je_linf = global_linf_difference(step.j_bkg_face_mid,
                                                      step.j_bkg_energy_debug_face, sg);
    const double max_abs_jn = global_linf(step.j_bkg_face_mid, sg);
    const double max_abs_gstar_je = global_linf(step.j_bkg_energy_debug_face, sg);
    const double j_scale = std::max(1.0, std::max(
        max_abs_jn, max_abs_gstar_je));
    // These are the production midpoint pairings: <E_f^mid, J_N> and
    // <G E_f^mid, J_E>.  Reuse them instead of recomputing with E^(n+1).
    const double e_dot_jn = step.current_diag.e_dot_j_charge;
    const double e_dot_gstar_je = step.current_diag.e_dot_j_energy;
    const double beam_current_linf = global_linf(step.j_beam_face_mid, sg);

    bool local_finite = std::isfinite(delta_k) &&
        std::isfinite(delta_k_subtractive) &&
        std::isfinite(delta_k_method_difference) && std::isfinite(delta_u) &&
        std::isfinite(energy_residual) && std::isfinite(jn_je_linf) &&
        std::isfinite(max_abs_jn) && std::isfinite(max_abs_gstar_je) &&
        std::isfinite(e_dot_jn) && std::isfinite(e_dot_gstar_je) &&
        (!transport.require_fct_activation ||
         (std::isfinite(step.fct_high_candidate_min) &&
          std::isfinite(step.fct_high_candidate_donor_excess)));
    local_finite = local_finite &&
        std::isfinite(initial.number_density) &&
        std::isfinite(initial.mean_vx_relative_error) &&
        std::isfinite(initial.current_relative_error) &&
        std::isfinite(initial.min_upar_width) &&
        std::isfinite(initial.min_uperp_width) &&
        std::isfinite(initial.upar_width_at_drift);
    for (size_t i = 0; i < step.fields_np1.Ex_face.size(); ++i)
        local_finite = local_finite && std::isfinite(step.fields_np1.Ex_face[i]);

    Audit audit = Audit();
    audit.dt = dt;
    audit.delta_k = delta_k;
    audit.delta_k_subtractive = delta_k_subtractive;
    audit.delta_k_method_difference = delta_k_method_difference;
    audit.delta_u = delta_u;
    audit.max_abs_jn = max_abs_jn;
    audit.max_abs_gstar_je = max_abs_gstar_je;
    audit.e_dot_jn = e_dot_jn;
    audit.e_dot_gstar_je = e_dot_gstar_je;
    audit.energy_exchange_scale = energy_scale;
    audit.energy_residual_abs = std::fabs(energy_residual);
    audit.energy_relative = audit.energy_residual_abs / energy_scale;
    audit.j_pair_relative = jn_je_linf / j_scale;
    audit.number_relative = std::fabs(n1 - n0) /
        std::max(1.0, std::fabs(n0));
    audit.beam_current_linf = beam_current_linf;
    audit.limiter_active_fraction = step.limiter_active_fraction;
    audit.limiter_min_alpha = step.limiter_min_alpha;
    audit.limiter_active_fraction_core = step.limiter_active_fraction_core;
    audit.limiter_active_fraction_boundary = step.limiter_active_fraction_boundary;
    audit.limiter_min_alpha_core = step.limiter_min_alpha_core;
    audit.limiter_min_alpha_boundary = step.limiter_min_alpha_boundary;
    audit.stage5_r_fv = step.stage5_r_fv;
    audit.stage5_r_couple = step.stage5_r_couple;
    audit.stage5_r_couple_centered = step.stage5_r_couple_centered;
    audit.stage5_r_couple_upwind_stabilization =
        step.stage5_r_couple_upwind_stabilization;
    audit.stage5_r_couple_fct_stabilization =
        step.stage5_r_couple_fct_stabilization;
    audit.stage5_r_fv_relative = std::fabs(step.stage5_r_fv) /
        std::max(kEnergyScaleFloor, energy_scale);
    audit.fct_high_candidate_min = step.fct_high_candidate_min;
    audit.fct_high_candidate_donor_excess =
        step.fct_high_candidate_donor_excess;
    audit.fct_controlled_injection_count =
        step.fct_controlled_injection_count;
    audit.state_advanced = step.state_advanced ? 1 : 0;
    audit.converged = step.converged ? 1 : 0;
    audit.soft_unconverged = step.soft_unconverged ? 1 : 0;
    audit.failed = step.failed ? 1 : 0;
    audit.low_order_only = step.transport_low_order_only;
    audit.fct_enabled = transport.fct_enabled ? 1 : 0;
    audit.finite = global_all_true(local_finite) ? 1 : 0;
    audit.initial = initial;
    return audit;
}

void print_audit(const char* label, const Audit& audit)
{
    std::cout << label << " dt=" << audit.dt
              << " state_advanced=" << audit.state_advanced
              << " converged=" << audit.converged
              << " soft_unconverged=" << audit.soft_unconverged
              << " failed=" << audit.failed
              << " low_order_only=" << audit.low_order_only
              << " fct_enabled=" << audit.fct_enabled
              << " finite=" << audit.finite << "\n"
              << label << " delta_K_bkg=" << audit.delta_k
              << " subtractive_delta_K_bkg=" << audit.delta_k_subtractive
              << " delta_K_method_difference=" << audit.delta_k_method_difference
              << " delta_U_E=" << audit.delta_u
              << " max_abs_JN=" << audit.max_abs_jn
              << " max_abs_GstarJE=" << audit.max_abs_gstar_je
              << " E_dot_JN=" << audit.e_dot_jn
              << " E_dot_GstarJE=" << audit.e_dot_gstar_je
              << " energy_exchange_scale=" << audit.energy_exchange_scale
              << " abs_delta_K_plus_delta_U=" << audit.energy_residual_abs
              << " total_energy_relative=" << audit.energy_relative << "\n"
              << label << " JN_minus_GstarJE_relative=" << audit.j_pair_relative
              << " number_relative_change=" << audit.number_relative
              << " beam_current_linf=" << audit.beam_current_linf
              << " limiter_active_fraction=" << audit.limiter_active_fraction
              << " limiter_min_alpha=" << audit.limiter_min_alpha
              << " limiter_active_fraction_core=" << audit.limiter_active_fraction_core
              << " limiter_active_fraction_boundary=" << audit.limiter_active_fraction_boundary
              << " limiter_min_alpha_core=" << audit.limiter_min_alpha_core
              << " limiter_min_alpha_boundary=" << audit.limiter_min_alpha_boundary
              << " stage5_R_FV=" << audit.stage5_r_fv
              << " stage5_R_FV_relative=" << audit.stage5_r_fv_relative
              << " fct_high_candidate_min=" << audit.fct_high_candidate_min
              << " fct_high_candidate_donor_excess="
              << audit.fct_high_candidate_donor_excess
              << " fct_controlled_injection_count="
              << audit.fct_controlled_injection_count
              << " stage5_R_couple=" << audit.stage5_r_couple
              << " stage5_R_couple_centered=" << audit.stage5_r_couple_centered
              << " stage5_R_couple_upwind_stabilization="
              << audit.stage5_r_couple_upwind_stabilization
              << " stage5_R_couple_fct_stabilization="
              << audit.stage5_r_couple_fct_stabilization << "\n";
}

void print_initial_audit(const InitializationAudit& audit)
{
    std::cout << "initial_number_density=" << audit.number_density
              << " initial_mean_vx=" << audit.mean_vx
              << " initial_current_density=" << audit.current_density
              << " target_mean_vx=" << audit.target_mean_vx
              << " target_current_density=" << audit.target_current_density
              << " mean_vx_relative_error=" << audit.mean_vx_relative_error
              << " current_relative_error=" << audit.current_relative_error << "\n"
              << "min_upar_width=" << audit.min_upar_width
              << " min_uperp_width=" << audit.min_uperp_width
              << " upar_width_at_drift=" << audit.upar_width_at_drift
              << " cells_per_uth_upar_at_drift="
              << audit.cells_per_uth_upar_at_drift
              << " cells_per_uth_uperp_at_axis="
              << audit.cells_per_uth_uperp_at_axis
              << " occupied_upar_cells=" << audit.occupied_upar_cells
              << " occupied_uperp_cells=" << audit.occupied_uperp_cells << "\n";
}

bool strict_case_accepts(const Audit& audit, const TransportMode& transport)
{
    return audit.finite != 0 && audit.state_advanced != 0 &&
        audit.converged != 0 && audit.soft_unconverged == 0 &&
        audit.failed == 0 && audit.low_order_only == (transport.low_order_only ? 1 : 0) &&
        audit.fct_enabled == (transport.fct_enabled ? 1 : 0) &&
        audit.beam_current_linf == 0.0 &&
        audit.number_relative <= 1.0e-11 &&
        (!transport.fct_enabled || transport.require_fct_activation ||
         audit.limiter_active_fraction == 0.0) &&
        (!transport.fct_enabled || std::isfinite(audit.limiter_active_fraction)) &&
        (!transport.require_fct_activation ||
          (audit.limiter_active_fraction > 0.0 &&
          audit.limiter_min_alpha < 1.0 - 1.0e-14 &&
          audit.fct_controlled_injection_count > 0 &&
          (audit.fct_high_candidate_min < 0.0 ||
           audit.fct_high_candidate_donor_excess > 0.0) &&
          audit.stage5_r_fv_relative <= 1.0e-8)) &&
        audit.initial.min_upar_width > 0.0 &&
        audit.initial.min_uperp_width > 0.0 &&
        audit.initial.upar_width_at_drift > 0.0 &&
        audit.initial.occupied_upar_cells > 0 &&
        audit.initial.occupied_uperp_cells > 0;
}

bool write_result_file(const std::string& path, const Audit& full,
                       const Audit& half, const Audit& quarter,
                       const TransportMode& transport,
                       bool temporal_refinement_executed,
                       bool time_trend_compatible, bool local_passes)
{
    std::ofstream out(path.c_str());
    if (!out)
        return false;

    out << std::scientific << std::setprecision(17);
    out << "format_version=2\n";
    out << "transport_mode=" << transport.name << "\n";
    out << "transport_low_order_only=" << (transport.low_order_only ? 1 : 0) << "\n";
    out << "transport_fct_enabled=" << (transport.fct_enabled ? 1 : 0) << "\n";
    out << "transport_requires_fct_activation="
        << (transport.require_fct_activation ? 1 : 0) << "\n";
    out << "transport_final_flux_accounting_required="
        << (transport.require_fct_activation ? 1 : 0) << "\n";
    out << "temporal_refinement_executed="
        << (temporal_refinement_executed ? 1 : 0) << "\n";
    out << "test_run_id=" << FP_TEST_RUN_ID << "\n";
    out << "velocity_nv=" << Param::Nv << "\n";
    out << "velocity_nmu=" << Param::Nmu << "\n";
    out << "spatial_nx=" << Param::nx << "\n";
    out << "velocity_umax=" << Param::momentum_umax << "\n";
    out << "velocity_upar_stretch=" << Param::momentum_upar_stretch << "\n";
    out << "velocity_uperp_stretch=" << Param::momentum_uperp_stretch << "\n";
    out << "min_upar_width=" << full.initial.min_upar_width << "\n";
    out << "min_uperp_width=" << full.initial.min_uperp_width << "\n";
    out << "upar_width_at_drift=" << full.initial.upar_width_at_drift << "\n";
    out << "initial_number_density=" << full.initial.number_density << "\n";
    out << "initial_mean_vx=" << full.initial.mean_vx << "\n";
    out << "initial_current_density=" << full.initial.current_density << "\n";
    out << "target_mean_vx=" << full.initial.target_mean_vx << "\n";
    out << "target_current_density=" << full.initial.target_current_density << "\n";
    out << "initial_mean_vx_relative_error="
        << full.initial.mean_vx_relative_error << "\n";
    out << "initial_current_relative_error="
        << full.initial.current_relative_error << "\n";
    out << "cells_per_uth_upar_at_drift="
        << full.initial.cells_per_uth_upar_at_drift << "\n";
    out << "cells_per_uth_uperp_at_axis="
        << full.initial.cells_per_uth_uperp_at_axis << "\n";
    out << "occupied_upar_cells=" << full.initial.occupied_upar_cells << "\n";
    out << "occupied_uperp_cells=" << full.initial.occupied_uperp_cells << "\n";
    out << "dt=" << full.dt << "\n";
    out << "dt_over_2=" << half.dt << "\n";
    out << "dt_over_4=" << quarter.dt << "\n";
    out << "dt_energy_relative=" << full.energy_relative << "\n";
    out << "dt_over_2_energy_relative=" << half.energy_relative << "\n";
    out << "dt_over_4_energy_relative=" << quarter.energy_relative << "\n";
    out << "dt_energy_exchange_scale=" << full.energy_exchange_scale << "\n";
    out << "dt_over_2_energy_exchange_scale=" << half.energy_exchange_scale << "\n";
    out << "dt_over_4_energy_exchange_scale=" << quarter.energy_exchange_scale << "\n";
    out << "direct_delta_K_bkg=" << full.delta_k << "\n";
    out << "subtractive_delta_K_bkg=" << full.delta_k_subtractive << "\n";
    out << "delta_K_method_difference=" << full.delta_k_method_difference << "\n";
    out << "dt_over_2_direct_delta_K_bkg=" << half.delta_k << "\n";
    out << "dt_over_2_subtractive_delta_K_bkg=" << half.delta_k_subtractive << "\n";
    out << "dt_over_2_delta_K_method_difference="
        << half.delta_k_method_difference << "\n";
    out << "dt_over_4_direct_delta_K_bkg=" << quarter.delta_k << "\n";
    out << "dt_over_4_subtractive_delta_K_bkg=" << quarter.delta_k_subtractive << "\n";
    out << "dt_over_4_delta_K_method_difference="
        << quarter.delta_k_method_difference << "\n";
    out << "dt_delta_K_bkg=" << full.delta_k << "\n";
    out << "dt_over_2_delta_K_bkg=" << half.delta_k << "\n";
    out << "dt_over_4_delta_K_bkg=" << quarter.delta_k << "\n";
    out << "dt_delta_U_E=" << full.delta_u << "\n";
    out << "dt_over_2_delta_U_E=" << half.delta_u << "\n";
    out << "dt_over_4_delta_U_E=" << quarter.delta_u << "\n";
    out << "dt_max_abs_JN=" << full.max_abs_jn << "\n";
    out << "dt_over_2_max_abs_JN=" << half.max_abs_jn << "\n";
    out << "dt_over_4_max_abs_JN=" << quarter.max_abs_jn << "\n";
    out << "dt_max_abs_GstarJE=" << full.max_abs_gstar_je << "\n";
    out << "dt_over_2_max_abs_GstarJE=" << half.max_abs_gstar_je << "\n";
    out << "dt_over_4_max_abs_GstarJE=" << quarter.max_abs_gstar_je << "\n";
    out << "dt_E_dot_JN=" << full.e_dot_jn << "\n";
    out << "dt_over_2_E_dot_JN=" << half.e_dot_jn << "\n";
    out << "dt_over_4_E_dot_JN=" << quarter.e_dot_jn << "\n";
    out << "dt_E_dot_GstarJE=" << full.e_dot_gstar_je << "\n";
    out << "dt_over_2_E_dot_GstarJE=" << half.e_dot_gstar_je << "\n";
    out << "dt_over_4_E_dot_GstarJE=" << quarter.e_dot_gstar_je << "\n";
    out << "dt_j_pair_relative=" << full.j_pair_relative << "\n";
    out << "dt_over_2_j_pair_relative=" << half.j_pair_relative << "\n";
    out << "dt_over_4_j_pair_relative=" << quarter.j_pair_relative << "\n";
    out << "dt_energy_residual_abs=" << full.energy_residual_abs << "\n";
    out << "dt_over_2_energy_residual_abs=" << half.energy_residual_abs << "\n";
    out << "dt_over_4_energy_residual_abs=" << quarter.energy_residual_abs << "\n";
    out << "dt_stage5_R_FV=" << full.stage5_r_fv << "\n";
    out << "dt_over_2_stage5_R_FV=" << half.stage5_r_fv << "\n";
    out << "dt_over_4_stage5_R_FV=" << quarter.stage5_r_fv << "\n";
    out << "dt_stage5_R_FV_relative=" << full.stage5_r_fv_relative << "\n";
    out << "dt_over_2_stage5_R_FV_relative=" << half.stage5_r_fv_relative << "\n";
    out << "dt_over_4_stage5_R_FV_relative=" << quarter.stage5_r_fv_relative << "\n";
    out << "dt_fct_high_candidate_min=" << full.fct_high_candidate_min << "\n";
    out << "dt_over_2_fct_high_candidate_min=" << half.fct_high_candidate_min << "\n";
    out << "dt_over_4_fct_high_candidate_min=" << quarter.fct_high_candidate_min << "\n";
    out << "dt_fct_high_candidate_donor_excess="
        << full.fct_high_candidate_donor_excess << "\n";
    out << "dt_over_2_fct_high_candidate_donor_excess="
        << half.fct_high_candidate_donor_excess << "\n";
    out << "dt_over_4_fct_high_candidate_donor_excess="
        << quarter.fct_high_candidate_donor_excess << "\n";
    out << "dt_fct_controlled_injection_count="
        << full.fct_controlled_injection_count << "\n";
    out << "dt_fct_high_candidate_overshoot="
        << ((full.fct_high_candidate_min < 0.0 ||
             full.fct_high_candidate_donor_excess > 0.0) ? 1 : 0) << "\n";
    out << "dt_final_flux_accounting_pass="
        << ((std::isfinite(full.max_abs_jn) &&
             std::isfinite(full.max_abs_gstar_je) &&
             std::isfinite(full.e_dot_jn) &&
             std::isfinite(full.e_dot_gstar_je) &&
             full.stage5_r_fv_relative <= 1.0e-8) ? 1 : 0) << "\n";
    out << "dt_stage5_R_couple=" << full.stage5_r_couple << "\n";
    out << "dt_stage5_R_couple_centered=" << full.stage5_r_couple_centered << "\n";
    out << "dt_stage5_R_couple_upwind_stabilization="
        << full.stage5_r_couple_upwind_stabilization << "\n";
    out << "dt_stage5_R_couple_fct_stabilization="
        << full.stage5_r_couple_fct_stabilization << "\n";
    out << "dt_over_2_stage5_R_couple=" << half.stage5_r_couple << "\n";
    out << "dt_over_2_stage5_R_couple_centered="
        << half.stage5_r_couple_centered << "\n";
    out << "dt_over_2_stage5_R_couple_upwind_stabilization="
        << half.stage5_r_couple_upwind_stabilization << "\n";
    out << "dt_over_2_stage5_R_couple_fct_stabilization="
        << half.stage5_r_couple_fct_stabilization << "\n";
    out << "dt_over_4_stage5_R_couple=" << quarter.stage5_r_couple << "\n";
    out << "dt_over_4_stage5_R_couple_centered="
        << quarter.stage5_r_couple_centered << "\n";
    out << "dt_over_4_stage5_R_couple_upwind_stabilization="
        << quarter.stage5_r_couple_upwind_stabilization << "\n";
    out << "dt_over_4_stage5_R_couple_fct_stabilization="
        << quarter.stage5_r_couple_fct_stabilization << "\n";
    out << "dt_number_relative=" << full.number_relative << "\n";
    out << "dt_over_2_number_relative=" << half.number_relative << "\n";
    out << "dt_over_4_number_relative=" << quarter.number_relative << "\n";
    out << "dt_beam_current_linf=" << full.beam_current_linf << "\n";
    out << "dt_over_2_beam_current_linf=" << half.beam_current_linf << "\n";
    out << "dt_over_4_beam_current_linf=" << quarter.beam_current_linf << "\n";
    out << "dt_limiter_active_fraction=" << full.limiter_active_fraction << "\n";
    out << "dt_over_2_limiter_active_fraction=" << half.limiter_active_fraction << "\n";
    out << "dt_over_4_limiter_active_fraction="
        << quarter.limiter_active_fraction << "\n";
    out << "dt_limiter_min_alpha=" << full.limiter_min_alpha << "\n";
    out << "dt_over_2_limiter_min_alpha=" << half.limiter_min_alpha << "\n";
    out << "dt_over_4_limiter_min_alpha=" << quarter.limiter_min_alpha << "\n";
    out << "dt_limiter_active_fraction_core="
        << full.limiter_active_fraction_core << "\n";
    out << "dt_over_2_limiter_active_fraction_core="
        << half.limiter_active_fraction_core << "\n";
    out << "dt_over_4_limiter_active_fraction_core="
        << quarter.limiter_active_fraction_core << "\n";
    out << "dt_limiter_active_fraction_boundary="
        << full.limiter_active_fraction_boundary << "\n";
    out << "dt_over_2_limiter_active_fraction_boundary="
        << half.limiter_active_fraction_boundary << "\n";
    out << "dt_over_4_limiter_active_fraction_boundary="
        << quarter.limiter_active_fraction_boundary << "\n";
    out << "dt_limiter_min_alpha_core=" << full.limiter_min_alpha_core << "\n";
    out << "dt_over_2_limiter_min_alpha_core=" << half.limiter_min_alpha_core << "\n";
    out << "dt_over_4_limiter_min_alpha_core=" << quarter.limiter_min_alpha_core << "\n";
    out << "dt_limiter_min_alpha_boundary=" << full.limiter_min_alpha_boundary << "\n";
    out << "dt_over_2_limiter_min_alpha_boundary=" << half.limiter_min_alpha_boundary << "\n";
    out << "dt_over_4_limiter_min_alpha_boundary=" << quarter.limiter_min_alpha_boundary << "\n";
    out << "time_trend_compatible=" << (time_trend_compatible ? 1 : 0) << "\n";
    out << "dt_strict_accept=" << (strict_case_accepts(full, transport) ? 1 : 0) << "\n";
    out << "dt_over_2_strict_accept=" << (strict_case_accepts(half, transport) ? 1 : 0) << "\n";
    out << "dt_over_4_strict_accept=" << (strict_case_accepts(quarter, transport) ? 1 : 0) << "\n";
    out << "strict_cases_accept=" << (strict_case_accepts(full, transport) &&
        strict_case_accepts(half, transport) && strict_case_accepts(quarter, transport) ? 1 : 0) << "\n";
    out << "dt_state_advanced=" << full.state_advanced << "\n";
    out << "dt_over_2_state_advanced=" << half.state_advanced << "\n";
    out << "dt_over_4_state_advanced=" << quarter.state_advanced << "\n";
    out << "dt_converged=" << full.converged << "\n";
    out << "dt_over_2_converged=" << half.converged << "\n";
    out << "dt_over_4_converged=" << quarter.converged << "\n";
    out << "dt_soft_unconverged=" << full.soft_unconverged << "\n";
    out << "dt_over_2_soft_unconverged=" << half.soft_unconverged << "\n";
    out << "dt_over_4_soft_unconverged=" << quarter.soft_unconverged << "\n";
    out << "dt_failed=" << full.failed << "\n";
    out << "dt_over_2_failed=" << half.failed << "\n";
    out << "dt_over_4_failed=" << quarter.failed << "\n";
    out << "local_passes=" << (local_passes ? 1 : 0) << "\n";
    return static_cast<bool>(out);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Section 7.2 uses the production solver with Beam disabled.  The default
    // remains the low-order baseline; explicit CLI modes exercise the new
    // nonuniform MUSCL high-order path with or without its existing FCT.
    SpatialGrid sg;
    sg.init(rank, size);
    if (rank == 0) {
        const long double phase_cells = static_cast<long double>(sg.nx_total) *
            static_cast<long double>(Param::Nv) * static_cast<long double>(Param::Nmu);
        const long double one_array_mib = phase_cells * sizeof(double) /
            (1024.0L * 1024.0L);
        std::cerr << "production_7_2_start mpi_ranks=" << size
                  << " nx_local_rank0=" << sg.nx_local
                  << " Nv=" << Param::Nv << " Nmu=" << Param::Nmu
                  << " one_phase_space_array_MiB="
                  << static_cast<double>(one_array_mib) << "\n";
    }
    const TransportMode* transport = FP_FCT_ACTIVATION_TEST
        ? &kHighFctActivationMode : &kLowOrderMode;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--transport") {
            transport = parse_transport_mode(argv[++i]);
            if (!transport) {
                if (rank == 0)
                    std::cerr << "unknown --transport mode; use low-order, "
                              << "high-no-fct, high-fct, or high-fct-activation\n";
                MPI_Finalize();
                return 2;
            }
        }
    }

    // Production diagnostics already use output/.  Keep the section-7.2
    // machine-readable result there too, with a velocity-grid-specific name
    // so manually submitted coarse/base/fine jobs cannot overwrite each
    // other.  --result remains available for an explicit location.
    std::ostringstream default_result_path;
    default_result_path << FP_TEST_RESULT_DIR << "/production_ampere_7_2_Nv" << Param::Nv
                        << "_Nmu" << Param::Nmu;
    if (!transport->low_order_only)
        default_result_path << "_" << transport->name;
    default_result_path << ".result";
    std::string result_path = default_result_path.str();
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--result")
            result_path = argv[++i];
    }

    // Every velocity resolution uses this same physical dt, selected from
    // the finest-grid force CFL and the common x CFL.
    const double dt = common_velocity_suite_dt(sg);
    Species initial_reference;
    initial_reference.init("production_7_2_initial_reference",
                           SpeciesType::BACKGROUND_ELECTRON,
                           -Const::qe, Const::me, Param::dens,
                           Param::temperature_e, false, sg);
    initial_reference.initialize_maxwellian(kDriftVx);
    if (transport->require_fct_activation)
        apply_controlled_fct_profile(initial_reference, sg);
    const InitializationAudit shared_initial =
        audit_initial_distribution(initial_reference, sg);
    const Audit full = run_case(sg, rank, size, dt, shared_initial, *transport);
    // The controlled FCT regression is a single-step activation test.  Its
    // purpose is to demonstrate an actual raw-high overshoot and its final
    // flux accounting, not a temporal convergence sequence where smaller
    // steps may correctly require no limiting at all.
    const bool temporal_refinement_executed = !transport->require_fct_activation;
    const Audit half = temporal_refinement_executed
        ? run_case(sg, rank, size, 0.5 * dt, shared_initial, *transport) : full;
    const Audit quarter = temporal_refinement_executed
        ? run_case(sg, rank, size, 0.25 * dt, shared_initial, *transport) : full;

    // A fixed velocity-grid defect can make the normalized residual plateau
    // as dt decreases.  Require no material regression, rather than strict
    // decrease, and leave velocity-grid convergence to the suite audit.
    const bool time_trend_compatible = !temporal_refinement_executed ||
        (half.energy_residual_abs <= kTrendSlack * full.energy_residual_abs &&
         quarter.energy_residual_abs <= kTrendSlack * half.energy_residual_abs &&
         half.energy_relative <= kTrendSlack * full.energy_relative &&
         quarter.energy_relative <= kTrendSlack * half.energy_relative);
    const double energy_abs_order_1 = convergence_order(
        full.energy_residual_abs, half.energy_residual_abs);
    const double energy_abs_order_2 = convergence_order(
        half.energy_residual_abs, quarter.energy_residual_abs);
    const double energy_relative_order_1 = convergence_order(
        full.energy_relative, half.energy_relative);
    const double energy_relative_order_2 = convergence_order(
        half.energy_relative, quarter.energy_relative);
    const bool passes = global_all_true(strict_case_accepts(full, *transport) &&
        (!temporal_refinement_executed ||
         (strict_case_accepts(half, *transport) &&
          strict_case_accepts(quarter, *transport) && time_trend_compatible)));

    bool result_written = true;
    if (rank == 0 && !result_path.empty()) {
        result_written = ensure_result_parent_directory(result_path) &&
            write_result_file(result_path, full, half, quarter,
                              *transport, temporal_refinement_executed,
                              time_trend_compatible, passes);
    }
    result_written = global_all_true(result_written);

    if (rank == 0) {
        std::cout << std::scientific << std::setprecision(16)
                  << "production_self_consistent_ampere_test\n"
                  << "mpi_ranks=" << size << " beam_enabled=0 "
                  << "transport_mode=" << transport->name
                  << " low_order_only=" << (transport->low_order_only ? 1 : 0)
                  << " fct_enabled=" << (transport->fct_enabled ? 1 : 0)
                  << " nonuniform_muscl_enabled=" << (!transport->low_order_only ? 1 : 0)
                  << "\n"
                  << "production_solver=VlasovAmpereMidpointSolver"
                  << " velocity_nv=" << Param::Nv
                  << " velocity_nmu=" << Param::Nmu
                  << " velocity_umax=" << Param::momentum_umax
                  << " velocity_upar_stretch=" << Param::momentum_upar_stretch
                  << " velocity_uperp_stretch=" << Param::momentum_uperp_stretch
                  << " test_run_id=" << FP_TEST_RUN_ID
                  << " result_path=" << (result_path.empty() ? "none" : result_path) << "\n";
        print_initial_audit(full.initial);
        print_audit("dt", full);
        print_audit("dt_over_2", half);
        print_audit("dt_over_4", quarter);
        std::cout << "time_trend_compatible=" << time_trend_compatible
                  << " temporal_refinement_executed="
                  << temporal_refinement_executed
                  << " energy_abs_order_dt_to_dt_over_2=" << energy_abs_order_1
                  << " energy_abs_order_dt_over_2_to_dt_over_4=" << energy_abs_order_2
                  << " energy_relative_order_dt_to_dt_over_2=" << energy_relative_order_1
                  << " energy_relative_order_dt_over_2_to_dt_over_4=" << energy_relative_order_2
                  << " velocity_grid_suite_required=1"
                  << " local_result=" << (passes ? "PASS" : "FAIL")
                  << " result_file_written=" << result_written << "\n";
    }

    MPI_Finalize();
    return (passes && result_written) ? 0 : 1;
}
