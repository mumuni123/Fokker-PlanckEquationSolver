#include "joint_phase_space_midpoint.h"
#include "open_electrostatic_solver.h"
#include "maxwell.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

namespace {

size_t cell_index(int ix, int j, int k, int nv, int nmu)
{
    return (static_cast<size_t>(ix) * static_cast<size_t>(nv) +
            static_cast<size_t>(j)) * static_cast<size_t>(nmu) +
           static_cast<size_t>(k);
}

struct Scenario {
    const char* name;
    double field;
    int active_velocity_count;
};

struct ScenarioResult {
    std::string name;
    JointPhaseSpaceAuditResult audit;
    double mass_scale;
    double kinetic_scale;
    double poisson_scale;
    double g_gstar_scale;
    double potential_charge_abs_scale;
    double kinetic_roundoff_tolerance;
    double poisson_roundoff_tolerance;
    double g_gstar_roundoff_tolerance;
    double local_global_flux_relative_error;
    bool local_global_flux_pass;
    double u_energy_adjoint_relative_error;
    double u_energy_velocity_symmetry_error;
    double u_energy_velocity_max_abs;
    bool u_energy_adjoint_pass;
    double x_force_charge_adjoint_relative_error;
    bool x_force_charge_adjoint_pass;
    bool passed;
    ScenarioResult()
        : mass_scale(1.0), kinetic_scale(1.0), poisson_scale(1.0),
          g_gstar_scale(1.0), potential_charge_abs_scale(0.0),
          kinetic_roundoff_tolerance(0.0),
          poisson_roundoff_tolerance(0.0),
          g_gstar_roundoff_tolerance(0.0),
          local_global_flux_relative_error(0.0),
          local_global_flux_pass(false),
          u_energy_adjoint_relative_error(0.0),
          u_energy_velocity_symmetry_error(0.0),
          u_energy_velocity_max_abs(0.0),
          u_energy_adjoint_pass(false),
          x_force_charge_adjoint_relative_error(0.0),
          x_force_charge_adjoint_pass(false), passed(false) {}
};

struct UFluxGeometryResult {
    bool selected;
    bool passed;
    double relative_error_dx1;
    double relative_error_dx2;
    double dx1;
    double dx2;
    double ring_k1;
    double ring_k2;
    double cell_mass_relative_error;
    UFluxGeometryResult()
        : selected(false), passed(true), relative_error_dx1(0.0),
          relative_error_dx2(0.0), dx1(0.0), dx2(0.0), ring_k1(0.0),
          ring_k2(0.0), cell_mass_relative_error(0.0) {}
};

UFluxGeometryResult run_u_flux_geometry_test(int rank, int size)
{
    UFluxGeometryResult result;
    result.selected = true;
    if (size != 1) {
        result.passed = false;
        return result;
    }

    const double f0 = 3.0e20;
    const double field = 2.0e8;
    const int jface = 1;
    const int k1 = 1;
    const int k2 = 2;
    for (int case_id = 0; case_id < 2; ++case_id) {
        SpatialGrid sg;
        sg.init_with_domain(rank, size, 4,
                           case_id == 0 ? 4.0e-6 : 8.0e-6);
        CylindricalVelocityGrid vg;
        vg.init_grid(1.5, 32, 8, 32, 0, 1.5, 1.5, 2.0);
        const int nv = static_cast<int>(vg.upar_cells.size());
        const int nmu = static_cast<int>(vg.uperp_cells.size());
        const size_t count = static_cast<size_t>(sg.nx_global * nv * nmu);
        std::vector<double> mass(count, 0.0);
        for (int ix = 0; ix < sg.nx_global; ++ix)
            for (int j = 0; j < nv; ++j)
                for (int k = 0; k < nmu; ++k)
                    mass[cell_index(ix, j, k, nv, nmu)] =
                        f0 * sg.dx * vg.cell_phase_volume(j, k);
        for (int j = 0; j < nv; ++j)
            for (int k = 0; k < nmu; ++k) {
                const double expected_mass =
                    f0 * sg.dx * vg.cell_phase_volume(j, k);
                const double actual_mass = mass[cell_index(0, j, k, nv, nmu)];
                result.cell_mass_relative_error = std::max(
                    result.cell_mass_relative_error,
                    std::fabs(actual_mass - expected_mass) /
                        std::max(1.0, std::fabs(expected_mass)));
            }
        std::vector<double> e_cell(static_cast<size_t>(sg.nx_global), field);
        const JointPhaseSpaceFluxBundle bundle =
            JointPhaseSpaceMidpointOperator::build_periodic_center_flux(
                sg, vg, mass, e_cell, 1.0e-18);
        const int k = case_id == 0 ? k1 : k2;
        const double area = vg.uperp_ring_areas[static_cast<size_t>(k)];
        const double accel = field * (-Const::qe) /
            (Const::me * Const::c);
        const double expected = accel * f0 * sg.dx * area;
        const size_t flux_index =
            (static_cast<size_t>(jface) * static_cast<size_t>(nmu) +
             static_cast<size_t>(k));
        const double actual = bundle.u_flux_rate[flux_index];
        const double error = std::fabs(actual - expected) /
            std::max(1.0, std::fabs(expected));
        if (case_id == 0) {
            result.relative_error_dx1 = error;
            result.dx1 = sg.dx;
            result.ring_k1 = area;
        } else {
            result.relative_error_dx2 = error;
            result.dx2 = sg.dx;
            result.ring_k2 = area;
        }
    }
    result.passed = result.relative_error_dx1 <=
            4096.0 * std::numeric_limits<double>::epsilon() &&
        result.relative_error_dx2 <=
            4096.0 * std::numeric_limits<double>::epsilon() &&
        result.cell_mass_relative_error <=
            4096.0 * std::numeric_limits<double>::epsilon() &&
        result.dx1 != result.dx2 && result.ring_k1 != result.ring_k2;
    return result;
}

double face_inner_product(const std::vector<double>& e_face,
                          const std::vector<double>& current,
                          double dx)
{
    if (e_face.size() != current.size() || e_face.empty())
        return std::numeric_limits<double>::quiet_NaN();
    double value = 0.0;
    for (size_t i = 0; i < e_face.size(); ++i) {
        const double weight = (i == 0 || i + 1 == e_face.size()) ? 0.5 : 1.0;
        value += weight * e_face[i] * current[i];
    }
    return dx * value;
}

double face_abs_inner_product(const std::vector<double>& e_face,
                              const std::vector<double>& current,
                              double dx)
{
    if (e_face.size() != current.size() || e_face.empty())
        return std::numeric_limits<double>::quiet_NaN();
    double value = 0.0;
    for (size_t i = 0; i < e_face.size(); ++i) {
        const double weight = (i == 0 || i + 1 == e_face.size()) ? 0.5 : 1.0;
        value += weight * std::fabs(e_face[i] * current[i]);
    }
    return dx * value;
}

long double face_inner_product_extended(
    const std::vector<double>& e_face,
    const std::vector<double>& current, double dx)
{
    if (e_face.size() != current.size() || e_face.empty())
        return std::numeric_limits<long double>::quiet_NaN();
    long double value = 0.0L;
    for (size_t i = 0; i < e_face.size(); ++i) {
        const long double weight =
            (i == 0 || i + 1 == e_face.size()) ? 0.5L : 1.0L;
        value += weight * static_cast<long double>(e_face[i]) *
                 static_cast<long double>(current[i]);
    }
    return static_cast<long double>(dx) * value;
}

long double production_potential_work_extended(
    const EMFields& before, const EMFields& after,
    const std::vector<double>& rho_delta, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    long double work = 0.0L;
    for (int ix = 0; ix < nxl; ++ix) {
        const long double old_phi =
            static_cast<long double>(before.phi[static_cast<size_t>(ng + ix)]) +
            static_cast<long double>(sg.dx) *
            (static_cast<long double>(before.Ex_face[static_cast<size_t>(ix + 1)]) -
             static_cast<long double>(before.Ex_face[static_cast<size_t>(ix)])) /
            12.0L;
        const long double new_phi =
            static_cast<long double>(after.phi[static_cast<size_t>(ng + ix)]) +
            static_cast<long double>(sg.dx) *
            (static_cast<long double>(after.Ex_face[static_cast<size_t>(ix + 1)]) -
             static_cast<long double>(after.Ex_face[static_cast<size_t>(ix)])) /
            12.0L;
        work += 0.5L * (old_phi + new_phi) *
            static_cast<long double>(rho_delta[static_cast<size_t>(ng + ix)]) *
            static_cast<long double>(sg.dx);
    }
    return work;
}

struct PeriodicSeamAdjointResult {
    bool selected;
    bool helper_ok;
    bool endpoint_asymmetry_pass;
    bool weighted_adjoint_pass;
    bool prediction_nonzero_pass;
    bool prediction_identity_pass;
    bool passed;
    double weighted_adjoint_relative_error;
    double w_f_adjoint;
    double w_j;
    double w_f_naive;
    double old_naive_mismatch;
    double seam_predicted_residual;
    double prediction_error;
    double pairing_face_left;
    double pairing_face_right;
    double force_current_first_cell;
    double force_current_last_cell;
    PeriodicSeamAdjointResult()
        : selected(false), helper_ok(false), endpoint_asymmetry_pass(false),
          weighted_adjoint_pass(false), prediction_nonzero_pass(false),
          prediction_identity_pass(false), passed(false),
          weighted_adjoint_relative_error(0.0), w_f_adjoint(0.0), w_j(0.0),
          w_f_naive(0.0), old_naive_mismatch(0.0),
          seam_predicted_residual(0.0), prediction_error(0.0),
          pairing_face_left(0.0), pairing_face_right(0.0),
          force_current_first_cell(0.0), force_current_last_cell(0.0) {}
};

// J0-E2 (stage B4): isolated verification that the production periodic-seam
// weighted adjoint helper closes W_F == W_J globally, including the seam
// faces 0 and Nx that the pre-existing J0-E interior-face loop never visits.
PeriodicSeamAdjointResult run_periodic_seam_adjoint_test(int rank, int size)
{
    PeriodicSeamAdjointResult result;
    result.selected = true;
    const double eps_gate = 8192.0 * std::numeric_limits<double>::epsilon();

    SpatialGrid sg;
    sg.init_with_domain(rank, size, 4, 4.0e-6);
    CylindricalVelocityGrid vg;
    vg.init_grid(1.5, 32, 8, 32, 0, 1.5, 1.5, 2.0);
    const int nx = sg.nx_global;
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    const size_t count = static_cast<size_t>(nx * nv * nmu);
    const double u_extent = std::max(
        std::fabs(vg.upar_cells.front()), std::fabs(vg.upar_cells.back()));

    // B4.1: deterministic positive midpoint mass with no u_parallel symmetry
    // and a monotone x tilt, so the first-cell and last-cell force currents
    // cannot coincide.
    std::vector<double> m_mid(count, 0.0);
    for (int ix = 0; ix < nx; ++ix) {
        const double x_tilt =
            1.0 + 0.25 * ((static_cast<double>(ix) + 0.5) /
                          static_cast<double>(nx) - 0.5);
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const double u = vg.upar_cells[static_cast<size_t>(j)];
                const double up = vg.uperp_cells[static_cast<size_t>(k)];
                const double shape = std::exp(-(u * u + up * up) / 0.35);
                const double asymmetric_factor =
                    1.0 + 0.15 * u / u_extent + 0.03 * std::sin(
                        2.0 * Const::pi * static_cast<double>((ix + 1) *
                        (k + 1)) / static_cast<double>(nx * nmu));
                m_mid[cell_index(ix, j, k, nv, nmu)] =
                    1.0e20 * x_tilt * asymmetric_factor * shape *
                    vg.cell_phase_volume(j, k);
            }
        }
    }
    // Artificial non-periodic pairing face with distinct physical endpoints.
    std::vector<double> pairing_face(static_cast<size_t>(nx) + 1, 0.0);
    for (int f = 0; f <= nx; ++f) {
        pairing_face[static_cast<size_t>(f)] =
            2.0e8 * (1.0 + 0.5 * std::sin(
                2.0 * Const::pi * static_cast<double>(f) /
                static_cast<double>(nx))) +
            5.0e7 * static_cast<double>(f) / static_cast<double>(nx);
    }

    std::vector<double> e_cell(static_cast<size_t>(nx), 2.0e8);
    const double dt = 1.0e-18;
    // B4.3: production x bundle only; charge_current_face is never rebuilt.
    const JointPhaseSpaceFluxBundle bundle =
        JointPhaseSpaceMidpointOperator::build_periodic_center_flux(
            sg, vg, m_mid, e_cell, dt);
    // B4.2: force current from the production Hamiltonian velocity.
    const std::vector<double> hamiltonian_velocity =
        JointPhaseSpaceMidpointOperator::build_hamiltonian_velocity(vg);
    std::vector<double> force_current(static_cast<size_t>(nx), 0.0);
    for (int ix = 0; ix < nx; ++ix) {
        long double sum = 0.0L;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t q = static_cast<size_t>(j * nmu + k);
                sum += static_cast<long double>(hamiltonian_velocity[q]) *
                    static_cast<long double>(
                        m_mid[cell_index(ix, j, k, nv, nmu)]);
            }
        }
        force_current[static_cast<size_t>(ix)] = static_cast<double>(
            static_cast<long double>(-Const::qe) * sum /
            static_cast<long double>(sg.dx));
    }

    result.pairing_face_left = pairing_face.front();
    result.pairing_face_right = pairing_face.back();
    result.force_current_first_cell = force_current.front();
    result.force_current_last_cell = force_current.back();
    result.endpoint_asymmetry_pass =
        result.pairing_face_left != result.pairing_face_right &&
        result.force_current_first_cell != result.force_current_last_cell;

    // B4.4: production weighted-adjoint helper.
    std::vector<double> e_adjoint;
    result.helper_ok =
        JointPhaseSpaceMidpointOperator::build_periodic_x_adjoint_cell_field(
            sg, pairing_face, rank, size, e_adjoint);

    // B4.5: two independent work sums, no shared helper.
    long double w_f_adjoint = 0.0L;
    long double w_f_naive = 0.0L;
    for (int ix = 0; ix < nx; ++ix) {
        w_f_adjoint += static_cast<long double>(dt) *
            static_cast<long double>(sg.dx) *
            static_cast<long double>(e_adjoint[static_cast<size_t>(ix)]) *
            static_cast<long double>(force_current[static_cast<size_t>(ix)]);
        const double naive_field = 0.5 *
            (pairing_face[static_cast<size_t>(ix)] +
             pairing_face[static_cast<size_t>(ix + 1)]);
        w_f_naive += static_cast<long double>(dt) *
            static_cast<long double>(sg.dx) *
            static_cast<long double>(naive_field) *
            static_cast<long double>(force_current[static_cast<size_t>(ix)]);
    }
    long double w_j = 0.0L;
    for (int iface = 0; iface <= nx; ++iface) {
        const double weight =
            (iface == 0 || iface == nx) ? 0.5 : 1.0;
        w_j += static_cast<long double>(dt) *
            static_cast<long double>(weight * sg.dx) *
            static_cast<long double>(pairing_face[static_cast<size_t>(iface)]) *
            static_cast<long double>(bundle.charge_current_face[
                static_cast<size_t>(iface)]);
    }
    result.w_f_adjoint = static_cast<double>(w_f_adjoint);
    result.w_f_naive = static_cast<double>(w_f_naive);
    result.w_j = static_cast<double>(w_j);
    result.old_naive_mismatch =
        result.w_f_naive - result.w_j;
    // B4.6: analytic seam prediction of the old naive-gather mismatch.
    result.seam_predicted_residual = 0.25 * dt * sg.dx *
        (result.pairing_face_left - result.pairing_face_right) *
        (result.force_current_first_cell -
         result.force_current_last_cell);
    result.prediction_error = std::fabs(result.old_naive_mismatch -
                                        result.seam_predicted_residual);

    const double work_scale = std::max(
        1.0, std::max(std::fabs(result.w_f_adjoint),
                      std::fabs(result.w_j)));
    result.weighted_adjoint_relative_error =
        std::fabs(result.w_f_adjoint - result.w_j) / work_scale;
    result.weighted_adjoint_pass = result.helper_ok &&
        result.weighted_adjoint_relative_error <= eps_gate;
    const double naive_scale = std::max(
        1.0, std::max(std::fabs(result.w_f_naive),
                      std::fabs(result.w_j)));
    result.prediction_nonzero_pass =
        std::fabs(result.old_naive_mismatch) >
        eps_gate * naive_scale;
    result.prediction_identity_pass = result.prediction_error <=
        eps_gate * std::max(1.0, std::fabs(result.seam_predicted_residual));
    result.passed = result.helper_ok && result.endpoint_asymmetry_pass &&
        result.weighted_adjoint_pass && result.prediction_nonzero_pass &&
        result.prediction_identity_pass;
    return result;
}

struct A2CaseResult {
    bool ran;
    bool finite;
    bool poisson_scalar_identity_pass;
    bool charge_residual_projection_pass;
    bool poisson_current_prediction_pass;
    bool nontrivial_wc_pass;
    bool passed;
    double w_c;
    double r_p;
    double r_pj;
    double prediction_error;
    double tau_c;
    double tau_a;
    double mismatch_linf;
    double r_q_linf;
    double r_c_linf;
    double u_boundary_linf;
    double total_charge_change;
    A2CaseResult()
        : ran(false), finite(false),
          poisson_scalar_identity_pass(false),
          charge_residual_projection_pass(false),
          poisson_current_prediction_pass(false),
          nontrivial_wc_pass(false), passed(false),
          w_c(0.0), r_p(0.0), r_pj(0.0), prediction_error(0.0),
          tau_c(0.0), tau_a(0.0), mismatch_linf(0.0), r_q_linf(0.0),
          r_c_linf(0.0), u_boundary_linf(0.0),
          total_charge_change(0.0) {}
};

struct A2PairResult {
    bool selected;
    A2CaseResult zero;
    A2CaseResult injected;
    bool passed;
    A2PairResult() : selected(false), passed(false) {}
};

// Stage-A-S1 unit check (section 7B): the stable incremental charge
// assembly is algebraically identical to qe*(ni - n_e) but numerically
// superior for two near-neutral large totals.  This test verifies, with a
// fixed ion background and near-neutral electron states under positive,
// negative and tiny perturbations, that the production incremental form
// (per-velocity-cell long double differencing, one final rounding) matches
// an all-long-double high precision reference, while the legacy absolute
// form does not.
struct As1IncrementalRhoResult {
    bool selected;
    bool positive_pass;
    bool negative_pass;
    bool near_neutral_small_pass;
    bool absolute_form_degrades;
    bool passed;
    double positive_incremental_error;
    double positive_absolute_error;
    double negative_incremental_error;
    double near_neutral_incremental_error;
    double near_neutral_absolute_error;
    double roundoff_bound;
    As1IncrementalRhoResult()
        : selected(false), positive_pass(false), negative_pass(false),
          near_neutral_small_pass(false),
          absolute_form_degrades(false), passed(false),
          positive_incremental_error(0.0),
          positive_absolute_error(0.0),
          negative_incremental_error(0.0),
          near_neutral_incremental_error(0.0),
          near_neutral_absolute_error(0.0), roundoff_bound(0.0) {}
};

As1IncrementalRhoResult run_as1_incremental_rho_test(int rank, int size)
{
    (void)rank;
    (void)size;
    As1IncrementalRhoResult result;
    result.selected = true;
    const double eps_gate = 8192.0 *
        std::numeric_limits<double>::epsilon();

    // Fixed ion background and a near-neutral electron state: n_e^0 differs
    // from n_i by only 1e-13 relative, so rho^0 is a tiny difference of two
    // large numbers - exactly the regime where the absolute form degrades.
    const double ni = 3.0e20;
    const double dx = 4.0e-6;
    const int nq_local = 32 * 8;
    const long double neutral_offset = 1.0e-13L;
    const double m_cell =
        static_cast<double>(static_cast<long double>(ni) *
                            (1.0L + neutral_offset) *
                            static_cast<long double>(dx) /
                            static_cast<long double>(nq_local));
    const long double rho_old_ld =
        static_cast<long double>(Const::qe) *
        (static_cast<long double>(ni) -
         static_cast<long double>(m_cell) *
             static_cast<long double>(nq_local) /
             static_cast<long double>(dx));

    auto run_subcase = [&](double relative_perturbation)
        -> std::pair<double, double> {
        // candidate = m_old + perturbation applied per velocity cell.
        const double delta_m = m_cell * relative_perturbation;
        long double delta_number_ref = 0.0L;
        for (int q = 0; q < nq_local; ++q)
            delta_number_ref += static_cast<long double>(delta_m);
        const long double delta_rho_ref =
            -static_cast<long double>(Const::qe) * delta_number_ref /
            static_cast<long double>(dx);
        // High precision reference.
        const long double rho_reference = rho_old_ld + delta_rho_ref;
        // Production incremental form semantics: per-velocity-cell
        // long double differencing of double values, one final rounding,
        // added to the stored old rho in double.
        long double delta_number = 0.0L;
        for (int q = 0; q < nq_local; ++q)
            delta_number +=
                static_cast<long double>(m_cell + delta_m) -
                static_cast<long double>(m_cell);
        const long double delta_rho =
            -static_cast<long double>(Const::qe) * delta_number /
            static_cast<long double>(dx);
        const double rho_incremental = static_cast<double>(
            rho_old_ld + delta_rho);
        // Legacy absolute form semantics: plain double accumulation of the
        // candidate total, then qe*(ni - number/dx).
        double number = 0.0;
        for (int q = 0; q < nq_local; ++q) number += m_cell + delta_m;
        const double ni_double = static_cast<double>(
            static_cast<long double>(ni));
        const double rho_absolute =
            Const::qe * (ni_double - number / dx);

        const long double parent_scale =
            std::max(std::fabs(static_cast<long double>(Const::qe) *
                               static_cast<long double>(ni)),
                     static_cast<long double>(Const::qe) *
                         static_cast<long double>(m_cell + delta_m) *
                         static_cast<long double>(nq_local) /
                         static_cast<long double>(dx)) +
            std::fabs(rho_old_ld);
        result.roundoff_bound = std::max(result.roundoff_bound,
            static_cast<double>(32.0L *
                (static_cast<long double>(nq_local + 8) *
                 static_cast<long double>(
                     std::numeric_limits<double>::epsilon()) /
                 (1.0L - static_cast<long double>(nq_local + 8) *
                  static_cast<long double>(
                      std::numeric_limits<double>::epsilon()))) *
                std::max(1.0L, parent_scale)));
        return {std::fabs(static_cast<double>(
                    static_cast<long double>(rho_incremental) -
                    rho_reference)),
                std::fabs(static_cast<double>(
                    static_cast<long double>(rho_absolute) -
                    rho_reference))};
    };

    std::pair<double, double> pos = run_subcase(1.0e-6);
    result.positive_incremental_error = pos.first;
    result.positive_absolute_error = pos.second;
    std::pair<double, double> neg = run_subcase(-1.0e-6);
    result.negative_incremental_error = neg.first;
    std::pair<double, double> small = run_subcase(1.0e-15);
    result.near_neutral_incremental_error = small.first;
    result.near_neutral_absolute_error = small.second;

    const double scale = std::fabs(static_cast<double>(rho_old_ld));
    const double tol = eps_gate * std::max(scale, 1.0);
    result.positive_pass = pos.first <= tol;
    result.negative_pass = neg.first <= tol;
    result.near_neutral_small_pass = small.first <= tol;
    result.absolute_form_degrades =
        result.roundoff_bound > 0.0 &&
        pos.second > pos.first && small.second > small.first;
    result.passed = result.positive_pass && result.negative_pass &&
        result.near_neutral_small_pass &&
        result.absolute_form_degrades;
    return result;
}

// Stage-A2 (docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md
// section 6): direct production-object identity test.  Both cases call
// JointPhaseSpaceMidpointOperator::evaluate_local_residual(),
// OpenElectrostaticSolver::solve(), evaluate_work_identity() and
// build_potential_pairing_field(); no production formula is re-derived in
// this test.
A2CaseResult run_a2_case(bool inject_residual, int rank, int size)
{
    A2CaseResult result;
    result.ran = true;
    const double eps_gate = 8192.0 *
        std::numeric_limits<double>::epsilon();

    SpatialGrid sg;
    sg.init_with_domain(rank, size, 4, 4.0e-6);
    CylindricalVelocityGrid vg;
    vg.init_grid(1.5, 32, 8, 32, 0, 1.5, 1.5, 2.0);
    const int nx = sg.nx_global;
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    const size_t count = static_cast<size_t>(nx * nv * nmu);
    const double u_extent = std::max(
        std::fabs(vg.upar_cells.front()),
        std::fabs(vg.upar_cells.back()));

    // Deterministic positive midpoint mass with no u_parallel symmetry and a
    // monotone x tilt (same construction as the J0-E2 seam test).
    std::vector<double> m_mid(count, 0.0);
    for (int ix = 0; ix < nx; ++ix) {
        const double x_tilt =
            1.0 + 0.25 * ((static_cast<double>(ix) + 0.5) /
                          static_cast<double>(nx) - 0.5);
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const double u = vg.upar_cells[static_cast<size_t>(j)];
                const double up = vg.uperp_cells[static_cast<size_t>(k)];
                const double shape = std::exp(-(u * u + up * up) / 0.35);
                const double asymmetric_factor =
                    1.0 + 0.15 * u / u_extent + 0.03 * std::sin(
                        2.0 * Const::pi * static_cast<double>((ix + 1) *
                        (k + 1)) / static_cast<double>(nx * nmu));
                m_mid[cell_index(ix, j, k, nv, nmu)] =
                    1.0e20 * x_tilt * asymmetric_factor * shape *
                    vg.cell_phase_volume(j, k);
            }
        }
    }
    std::vector<double> e_cell(static_cast<size_t>(nx), 2.0e8);
    const double dt = 1.0e-18;
    const JointPhaseSpaceFluxBundle bundle =
        JointPhaseSpaceMidpointOperator::build_periodic_center_flux(
            sg, vg, m_mid, e_cell, dt);
    std::vector<double> m_old(count, 0.0);
    std::vector<double> m_new(count, 0.0);
    for (size_t i = 0; i < count; ++i) {
        m_old[i] = m_mid[i] - 0.5 * bundle.mass_delta_total[i];
        m_new[i] = m_mid[i] + 0.5 * bundle.mass_delta_total[i];
    }

    // Deterministic, x-asymmetric, strictly positive injection into the
    // candidate only, so the phase-space residual becomes nonzero while the
    // candidate stays positive and the total injected charge is controllable.
    std::vector<double> injection(count, 0.0);
    if (inject_residual) {
        long double local_total = 0.0L;
        for (int ix = 0; ix < nx; ++ix) {
            const double x_bump =
                1.0 + 0.5 * ((static_cast<double>(ix) + 0.5) /
                             static_cast<double>(nx)) +
                0.25 * std::sin(6.0 * Const::pi *
                                (static_cast<double>(ix) + 0.5) /
                                static_cast<double>(nx));
            for (int j = 0; j < nv; ++j) {
                for (int k = 0; k < nmu; ++k) {
                    const double u =
                        vg.upar_cells[static_cast<size_t>(j)];
                    const double up =
                        vg.uperp_cells[static_cast<size_t>(k)];
                    const double shape =
                        std::exp(-(u * u + up * up) / 0.35);
                    const size_t id =
                        cell_index(ix, j, k, nv, nmu);
                    injection[id] =
                        1.0e-4 * 1.0e20 * x_bump * shape *
                        vg.cell_phase_volume(j, k);
                    m_new[id] += injection[id];
                    local_total += static_cast<long double>(
                        injection[id]);
                }
            }
        }
        result.total_charge_change = static_cast<double>(local_total);
    }

    // Production residual evaluation on the final candidate.
    JointPhaseSpaceFluxBundle a2_bundle;
    std::vector<double> residual;
    double residual_linf = 0.0;
    double residual_scale = 0.0;
    const bool eval_ok =
        JointPhaseSpaceMidpointOperator::evaluate_local_residual(
            sg, vg, m_old, m_new, e_cell, dt, rank, size,
            a2_bundle, residual, residual_linf, residual_scale);

    // Production Poisson states: the base charge change follows the same
    // production face current used by the residual algebra (never rebuilt
    // from M differences); the injection adds its own exact charge moment so
    // that r^C matches r^Q by construction instead of by symmetry.
    ElectrostaticBoundary boundary;
    boundary.type = ElectrostaticBoundaryType::DIRICHLET_PHI;
    boundary.e_left = 0.0;
    boundary.phi_left = 0.0;
    boundary.phi_right = 0.0;
    OpenElectrostaticSolver poisson;
    poisson.init(sg, boundary);
    EMFields before;
    EMFields after;
    before.init(sg);
    after.init(sg);
    std::fill(before.rho.begin(), before.rho.end(), 0.0);
    std::fill(after.rho.begin(), after.rho.end(), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double charge_delta = -dt *
            (a2_bundle.charge_current_face[
                static_cast<size_t>(ix) + 1] -
             a2_bundle.charge_current_face[
                 static_cast<size_t>(ix)]) / sg.dx;
        if (inject_residual) {
            long double s_sum = 0.0L;
            for (int j = 0; j < nv; ++j)
                for (int k = 0; k < nmu; ++k)
                    s_sum += static_cast<long double>(
                        injection[cell_index(ix, j, k, nv, nmu)]);
            charge_delta += -Const::qe * static_cast<double>(s_sum) /
                sg.dx;
        }
        after.rho[static_cast<size_t>(sg.nghost + ix)] = charge_delta;
    }
    OpenGaussSolveOptions solve_options;
    solve_options.reconstruct_phi = true;
    solve_options.compute_l1 = true;
    solve_options.compute_boundary_audit = true;
    poisson.solve(before, rank, size, solve_options);
    poisson.solve(after, rank, size, solve_options);

    std::vector<double> rho_delta(before.rho.size(), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix)
        rho_delta[static_cast<size_t>(sg.nghost + ix)] =
            after.rho[static_cast<size_t>(sg.nghost + ix)] -
            before.rho[static_cast<size_t>(sg.nghost + ix)];
    const OpenPoissonWorkIdentity work = poisson.evaluate_work_identity(
        before, after, rho_delta, rank, size);
    std::vector<double> pairing_face;
    const bool pairing_ok = poisson.build_potential_pairing_field(
        before, after, pairing_face, rank, size);

    // Section 2/3 diagnostics built only from production outputs.
    long double local_r_q_linf = 0.0L;
    long double local_r_c_linf = 0.0L;
    long double local_mismatch_linf = 0.0L;
    long double local_ubnd_linf = 0.0L;
    long double local_max_drho_dx = 0.0L;
    long double local_max_dt_dj = 0.0L;
    long double w_c = 0.0L;
    long double abs_phi_rc_sum = 0.0L;
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        long double sum_residual = 0.0L;
        for (int q = 0; q < nv * nmu; ++q)
            sum_residual += static_cast<long double>(
                residual[static_cast<size_t>(ix) *
                         static_cast<size_t>(nv * nmu) +
                         static_cast<size_t>(q)]);
        const long double rQ =
            -static_cast<long double>(Const::qe) * sum_residual;
        const long double delta_rho_dx =
            static_cast<long double>(rho_delta[
                static_cast<size_t>(ng + ix)]) *
            static_cast<long double>(sg.dx);
        const long double current_div_dt =
            static_cast<long double>(dt) *
            static_cast<long double>(
                a2_bundle.charge_current_face[
                    static_cast<size_t>(ix) + 1] -
                a2_bundle.charge_current_face[
                    static_cast<size_t>(ix)]);
        const long double rC = delta_rho_dx + current_div_dt;
        long double u_boundary_sum = 0.0L;
        for (int k = 0; k < nmu; ++k) {
            const size_t bottom =
                static_cast<size_t>(ix) *
                static_cast<size_t>(nv + 1) *
                static_cast<size_t>(nmu) +
                static_cast<size_t>(k);
            const size_t top = bottom +
                static_cast<size_t>(nv) *
                static_cast<size_t>(nmu);
            u_boundary_sum += static_cast<long double>(
                a2_bundle.u_flux_rate[top] -
                a2_bundle.u_flux_rate[bottom]);
        }
        const long double r_ub =
            -static_cast<long double>(Const::qe) *
            static_cast<long double>(dt) * u_boundary_sum;
        // Section 2.3 verbatim: the same cell-average potential used by
        // evaluate_work_identity()/build_potential_pairing_field().
        const double old_phi_average =
            before.phi[static_cast<size_t>(ng + ix)] +
            sg.dx *
            (before.Ex_face[static_cast<size_t>(ix) + 1] -
             before.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const double new_phi_average =
            after.phi[static_cast<size_t>(ng + ix)] +
            sg.dx *
            (after.Ex_face[static_cast<size_t>(ix) + 1] -
             after.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const long double phi_bar = static_cast<long double>(
            0.5 * (old_phi_average + new_phi_average));
        local_r_q_linf = std::max(local_r_q_linf, std::fabs(rQ));
        local_r_c_linf = std::max(local_r_c_linf, std::fabs(rC));
        local_mismatch_linf = std::max(local_mismatch_linf,
            std::fabs(rC - rQ + r_ub));
        local_ubnd_linf = std::max(local_ubnd_linf,
                                   std::fabs(r_ub));
        local_max_drho_dx = std::max(local_max_drho_dx,
                                     std::fabs(delta_rho_dx));
        local_max_dt_dj = std::max(local_max_dt_dj,
                                   std::fabs(current_div_dt));
        const long double phi_rc = phi_bar * rC;
        w_c += phi_rc;
        abs_phi_rc_sum += std::fabs(phi_rc);
    }
    const double w_j = dt * face_inner_product(
        pairing_face, a2_bundle.charge_current_face, sg.dx);
    result.w_c = static_cast<double>(w_c);
    result.r_p = work.residual;
    result.r_pj = work.field_energy_change - work.electrode_work + w_j;
    const double predicted = result.r_p + result.w_c;
    result.prediction_error = result.r_pj - predicted;
    result.tau_c = eps_gate *
        std::max(1.0, std::max(
            static_cast<double>(local_max_drho_dx),
            std::max(static_cast<double>(local_max_dt_dj),
                     static_cast<double>(local_r_q_linf))));
    result.tau_a = eps_gate *
        std::max(1.0, std::max(std::fabs(result.r_pj),
                               std::max(std::fabs(result.r_p),
                                        std::max(std::fabs(result.w_c),
                                         static_cast<double>(
                                             abs_phi_rc_sum)))));
    result.mismatch_linf = static_cast<double>(local_mismatch_linf);
    result.r_q_linf = static_cast<double>(local_r_q_linf);
    result.r_c_linf = static_cast<double>(local_r_c_linf);
    result.u_boundary_linf = static_cast<double>(local_ubnd_linf);
    const double poisson_scale = std::max(1.0e-300,
        std::max(std::fabs(work.field_energy_before),
        std::max(std::fabs(work.field_energy_after),
        std::max(std::fabs(work.field_energy_change),
                 std::fabs(work.potential_charge_work)))));
    result.finite = eval_ok && pairing_ok && work.finite &&
        std::isfinite(result.w_c) && std::isfinite(result.r_p) &&
        std::isfinite(result.r_pj) &&
        std::isfinite(result.prediction_error) &&
        std::isfinite(result.tau_c) && std::isfinite(result.tau_a);
    result.poisson_scalar_identity_pass = result.finite &&
        std::fabs(work.residual) <= eps_gate * poisson_scale;
    result.charge_residual_projection_pass = result.finite &&
        result.mismatch_linf <= result.tau_c;
    result.poisson_current_prediction_pass = result.finite &&
        std::fabs(result.prediction_error) <= result.tau_a;
    result.nontrivial_wc_pass = !inject_residual ||
        std::fabs(result.w_c) > 100.0 * result.tau_a;
    result.passed = result.finite &&
        result.poisson_scalar_identity_pass &&
        result.charge_residual_projection_pass &&
        result.poisson_current_prediction_pass &&
        result.nontrivial_wc_pass;
    return result;
}

ScenarioResult run_scenario(const Scenario& scenario, int rank, int size)
{
    ScenarioResult result;
    result.name = scenario.name;

    SpatialGrid sg;
    sg.init_with_domain(rank, size, 4, 4.0e-6);
    CylindricalVelocityGrid vg;
    // Nv=6 is too coarse for the existing production grid contract: its
    // low-speed energy-conjugate validation rejects the central cells before
    // J0 can run.  Nv=32/Nmu=8 remains a small J0 grid while retaining the
    // production nonuniform geometry and passing that pre-existing check.
    vg.init_grid(1.5, 32, 8, 32, 0, 1.5, 1.5, 2.0);
    const int nx = sg.nx_global;
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    const size_t count = static_cast<size_t>(nx * nv * nmu);
    const double u_extent = std::max(
        std::fabs(vg.upar_cells.front()), std::fabs(vg.upar_cells.back()));
    std::vector<double> m_mid(count, 0.0);
    const int center_j = nv / 2;
    for (int ix = 0; ix < nx; ++ix) {
        const double x_factor = 1.0 + 0.15 *
            std::sin(2.0 * Const::pi * (static_cast<double>(ix) + 0.5) /
                     static_cast<double>(nx));
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                bool active = scenario.active_velocity_count == 0;
                if (scenario.active_velocity_count == 1)
                    active = j == center_j && k == 1;
                if (scenario.active_velocity_count == 2)
                    active = (j == center_j || j == center_j - 1) && k == 1;
                const double u = vg.upar_cells[static_cast<size_t>(j)];
                const double up = vg.uperp_cells[static_cast<size_t>(k)];
                const double shape = std::exp(-(u * u + up * up) / 0.35);
                // F5.3 requires a positive deterministic state without the
                // u_parallel symmetry of a Maxwellian.  Otherwise the net
                // electric work is near zero and roundoff in two summation
                // orders is incorrectly reported as an O(1) relative error.
                const double asymmetric_factor =
                    scenario.active_velocity_count == 0
                    ? 1.0 + 0.15 * u / u_extent + 0.03 * std::sin(
                        2.0 * Const::pi * static_cast<double>((ix + 1) *
                        (k + 1)) / static_cast<double>(nx * nmu))
                    : 1.0;
                m_mid[cell_index(ix, j, k, nv, nmu)] = active
                    ? 1.0e20 * x_factor * asymmetric_factor * shape *
                       vg.cell_phase_volume(j, k)
                    : 0.0;
            }
        }
    }
    std::vector<double> e_cell(static_cast<size_t>(nx), scenario.field);
    const double dt = 1.0e-18;
    const JointPhaseSpaceFluxBundle bundle =
        JointPhaseSpaceMidpointOperator::build_periodic_center_flux(
            sg, vg, m_mid, e_cell, dt);
    std::vector<double> m_old(count, 0.0);
    std::vector<double> m_new(count, 0.0);
    for (size_t i = 0; i < count; ++i) {
        m_old[i] = m_mid[i] - 0.5 * bundle.mass_delta_total[i];
        m_new[i] = m_mid[i] + 0.5 * bundle.mass_delta_total[i];
    }
    JointPhaseSpaceFluxBundle global_bundle;
    JointPhaseSpaceFluxBundle local_bundle;
    std::vector<double> global_residual;
    std::vector<double> local_residual;
    double global_residual_linf = 0.0;
    double global_residual_scale = 0.0;
    double local_residual_linf = 0.0;
    double local_residual_scale = 0.0;
    const bool global_ok =
        JointPhaseSpaceMidpointOperator::evaluate_residual(
            sg, vg, m_old, m_new, e_cell, dt, global_bundle,
            global_residual, global_residual_linf, global_residual_scale);
    const std::vector<double> e_cell_local = e_cell;
    const bool local_ok =
        JointPhaseSpaceMidpointOperator::evaluate_local_residual(
            sg, vg, m_old, m_new, e_cell_local, dt, rank, size,
            local_bundle, local_residual, local_residual_linf,
            local_residual_scale,
            // F4 compares the algebraic midpoint operator with the global
            // residual on the same signed manufactured Newton probe.  It is
            // not an accepted physical state, so use the same signed domain
            // as evaluate_residual().
            true);
    double local_global_error = 0.0;
    if (global_ok && local_ok &&
        global_bundle.x_flux_rate.size() == local_bundle.x_flux_rate.size() &&
        global_bundle.u_flux_rate.size() == local_bundle.u_flux_rate.size() &&
        global_bundle.mass_delta_x.size() == local_bundle.mass_delta_x.size() &&
        global_bundle.mass_delta_u.size() == local_bundle.mass_delta_u.size() &&
        global_bundle.mass_delta_total.size() == local_bundle.mass_delta_total.size() &&
        global_residual.size() == local_residual.size()) {
        for (size_t i = 0; i < global_bundle.x_flux_rate.size(); ++i)
            local_global_error = std::max(local_global_error,
                std::fabs(global_bundle.x_flux_rate[i] -
                          local_bundle.x_flux_rate[i]));
        for (size_t i = 0; i < global_bundle.u_flux_rate.size(); ++i)
            local_global_error = std::max(local_global_error,
                std::fabs(global_bundle.u_flux_rate[i] -
                          local_bundle.u_flux_rate[i]));
        for (size_t i = 0; i < global_bundle.mass_delta_x.size(); ++i) {
            local_global_error = std::max(local_global_error,
                std::fabs(global_bundle.mass_delta_x[i] -
                          local_bundle.mass_delta_x[i]));
            local_global_error = std::max(local_global_error,
                std::fabs(global_bundle.mass_delta_u[i] -
                          local_bundle.mass_delta_u[i]));
            local_global_error = std::max(local_global_error,
                std::fabs(global_bundle.mass_delta_total[i] -
                          local_bundle.mass_delta_total[i]));
            local_global_error = std::max(local_global_error,
                std::fabs(global_residual[i] - local_residual[i]));
        }
    } else {
        local_global_error = std::numeric_limits<double>::infinity();
    }
    result.local_global_flux_relative_error = local_global_error /
        std::max(1.0, std::max(global_residual_scale,
                               std::max(std::fabs(global_residual_linf),
                                        std::fabs(local_residual_linf))));
    result.local_global_flux_pass = global_ok && local_ok &&
        result.local_global_flux_relative_error <=
            4096.0 * std::numeric_limits<double>::epsilon();
    const std::vector<double> hamiltonian_velocity =
        JointPhaseSpaceMidpointOperator::build_hamiltonian_velocity(vg);
    long double u_work_from_flux = 0.0L;
    long double u_work_from_velocity = 0.0L;
    for (int ix = 0; ix < nx; ++ix) {
        for (int jf = 1; jf < nv; ++jf) {
            for (int k = 0; k < nmu; ++k) {
                const size_t left = static_cast<size_t>((jf - 1) * nmu + k);
                const size_t right = static_cast<size_t>(jf * nmu + k);
                const double delta_k = vg.kinetic_energy[right] -
                    vg.kinetic_energy[left];
                const size_t flux_id =
                    (static_cast<size_t>(ix) * static_cast<size_t>(nv + 1) +
                     static_cast<size_t>(jf)) * static_cast<size_t>(nmu) +
                    static_cast<size_t>(k);
                u_work_from_flux += static_cast<long double>(dt) *
                    static_cast<long double>(delta_k) *
                    static_cast<long double>(bundle.u_flux_rate[flux_id]);
            }
        }
        // build_hamiltonian_velocity() is the face-to-cell transpose already
        // summed over all adjacent u faces.  Its work form is therefore a
        // cell sum, not another u-face sum.
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t velocity_id =
                    static_cast<size_t>(j * nmu + k);
                const size_t mass_id = cell_index(ix, j, k, nv, nmu);
                u_work_from_velocity += static_cast<long double>(dt) *
                    static_cast<long double>(e_cell[static_cast<size_t>(ix)]) *
                    static_cast<long double>(-Const::qe) *
                    static_cast<long double>(hamiltonian_velocity[velocity_id]) *
                    static_cast<long double>(m_mid[mass_id]);
            }
        }
    }
    const long double u_work_scale = std::max(
        1.0e-300L, std::max(std::fabs(u_work_from_flux),
                            std::fabs(u_work_from_velocity)));
    result.u_energy_adjoint_relative_error = static_cast<double>(
        std::fabs(u_work_from_flux - u_work_from_velocity) / u_work_scale);
    for (int j = 0; j < nv; ++j) {
        const int mirror = nv - 1 - j;
        for (int k = 0; k < nmu; ++k) {
            const size_t id = static_cast<size_t>(j * nmu + k);
            const size_t mirror_id = static_cast<size_t>(mirror * nmu + k);
            result.u_energy_velocity_max_abs = std::max(
                result.u_energy_velocity_max_abs,
                std::fabs(hamiltonian_velocity[id]));
            result.u_energy_velocity_symmetry_error = std::max(
                result.u_energy_velocity_symmetry_error,
                std::fabs(hamiltonian_velocity[id] +
                          hamiltonian_velocity[mirror_id]));
        }
    }
    result.u_energy_adjoint_pass =
        result.u_energy_adjoint_relative_error <=
            4096.0 * std::numeric_limits<double>::epsilon() &&
        result.u_energy_velocity_symmetry_error <=
            4096.0 * std::numeric_limits<double>::epsilon() * Const::c &&
        result.u_energy_velocity_max_abs <= Const::c *
            (1.0 + 4096.0 * std::numeric_limits<double>::epsilon());
    std::vector<double> force_current(static_cast<size_t>(nx), 0.0);
    for (int ix = 0; ix < nx; ++ix) {
        long double sum = 0.0L;
        for (int j = 0; j < nv; ++j)
            for (int k = 0; k < nmu; ++k) {
                const size_t q = static_cast<size_t>(j * nmu + k);
                sum += static_cast<long double>(hamiltonian_velocity[q]) *
                    static_cast<long double>(m_mid[cell_index(ix, j, k, nv, nmu)]);
            }
        force_current[static_cast<size_t>(ix)] = static_cast<double>(
            static_cast<long double>(-Const::qe) * sum /
            static_cast<long double>(sg.dx));
    }
    double x_force_charge_error = 0.0;
    for (int iface = 1; iface < nx; ++iface) {
        const double expected = 0.5 *
            (force_current[static_cast<size_t>(iface - 1)] +
             force_current[static_cast<size_t>(iface)]);
        x_force_charge_error = std::max(
            x_force_charge_error,
            std::fabs(bundle.charge_current_face[static_cast<size_t>(iface)] -
                      expected) /
                std::max(1.0, std::fabs(expected)));
    }
    result.x_force_charge_adjoint_relative_error = x_force_charge_error;
    result.x_force_charge_adjoint_pass = x_force_charge_error <=
        4096.0 * std::numeric_limits<double>::epsilon();

    ElectrostaticBoundary boundary;
    boundary.type = ElectrostaticBoundaryType::DIRICHLET_PHI;
    boundary.e_left = 0.0;
    boundary.phi_left = 0.0;
    boundary.phi_right = 0.0;
    OpenElectrostaticSolver poisson;
    poisson.init(sg, boundary);
    EMFields before;
    EMFields after;
    before.init(sg);
    after.init(sg);
    // Keep the Poisson pairing audit independent from large-background
    // density addition.  A tiny manufactured flux delta must be representable
    // in the stored rho field; otherwise after.rho-before.rho can erase it in
    // double precision before the production G/G* helper sees it.
    std::fill(before.rho.begin(), before.rho.end(), 0.0);
    // For the G/G* audit, define the candidate charge change directly from
    // the same x-face charge current used by the common flux.  Reconstructing
    // rho_delta by summing M_new-M_old would mix large-mass cancellation and
    // the u-direction telescoping roundoff into the Poisson pairing audit.
    std::fill(after.rho.begin(), after.rho.end(), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double charge_delta = -dt *
            (bundle.charge_current_face[static_cast<size_t>(ix + 1)] -
             bundle.charge_current_face[static_cast<size_t>(ix)]) / sg.dx;
        after.rho[static_cast<size_t>(sg.nghost + ix)] =
            before.rho[static_cast<size_t>(sg.nghost + ix)] + charge_delta;
    }
    OpenGaussSolveOptions solve_options;
    solve_options.reconstruct_phi = true;
    solve_options.compute_l1 = true;
    solve_options.compute_boundary_audit = true;
    poisson.solve(before, rank, size, solve_options);
    poisson.solve(after, rank, size, solve_options);

    std::vector<double> rho_delta(before.rho.size(), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix)
        rho_delta[static_cast<size_t>(sg.nghost + ix)] =
            after.rho[static_cast<size_t>(sg.nghost + ix)] -
            before.rho[static_cast<size_t>(sg.nghost + ix)];
    const OpenPoissonWorkIdentity work = poisson.evaluate_work_identity(
        before, after, rho_delta, rank, size);
    std::vector<double> pairing_face;
    const bool pairing_ok = poisson.build_potential_pairing_field(
        before, after, pairing_face, rank, size);
    const double pairing_work = pairing_ok
        ? -dt * face_inner_product(pairing_face,
                                   bundle.charge_current_face, sg.dx)
        : std::numeric_limits<double>::quiet_NaN();
    const double pairing_abs_scale = pairing_ok
        ? dt * face_abs_inner_product(pairing_face,
                                      bundle.charge_current_face, sg.dx)
        : std::numeric_limits<double>::quiet_NaN();
    const long double pairing_work_extended = pairing_ok
        ? -static_cast<long double>(dt) *
          face_inner_product_extended(pairing_face,
                                      bundle.charge_current_face, sg.dx)
        : std::numeric_limits<long double>::quiet_NaN();
    const long double potential_work_extended = pairing_ok
        ? production_potential_work_extended(before, after, rho_delta, sg)
        : std::numeric_limits<long double>::quiet_NaN();
    // -dt <E_pair, J_charge> is the same signed potential-charge work used
    // by OpenElectrostaticSolver::evaluate_work_identity().  The J0 dual
    // residual is therefore their difference, not their sum.
    const double g_gstar_residual = pairing_ok
        ? static_cast<double>(pairing_work_extended - potential_work_extended)
        : std::numeric_limits<double>::quiet_NaN();

    result.audit = JointPhaseSpaceMidpointOperator::audit(
        sg, vg, m_old, m_new, bundle, work.residual, g_gstar_residual);
    double local_mass_scale = 1.0;
    for (size_t i = 0; i < count; ++i)
        local_mass_scale = std::max(local_mass_scale,
                                    std::fabs(m_old[i]) + std::fabs(m_new[i]));
    const double local_kinetic_scale = std::max(
        1.0e-300, result.audit.kinetic_absolute_work_scale);
    result.mass_scale = local_mass_scale;
    result.kinetic_scale = local_kinetic_scale;
    const double eps = 8192.0 * std::numeric_limits<double>::epsilon();
    // OpenElectrostaticSolver evaluates the production Poisson identity in
    // double precision.  Its residual is formed by subtracting field-energy
    // values that can be many orders larger than the step work, so the
    // roundoff floor must include both endpoint field energies, not only
    // |Delta U| or |potential work|.
    result.poisson_scale = std::max(1.0e-300,
        std::max(std::fabs(work.field_energy_before),
        std::max(std::fabs(work.field_energy_after),
        std::max(std::fabs(work.field_energy_change),
                 std::fabs(work.potential_charge_work)))));
    result.g_gstar_scale = std::max(1.0e-300,
        std::max(std::fabs(work.potential_charge_work),
                 std::max(std::fabs(pairing_work), pairing_abs_scale)));
    result.kinetic_roundoff_tolerance = eps * result.kinetic_scale;
    result.poisson_roundoff_tolerance = eps * result.poisson_scale;
    double max_pairing_phi = 0.0;
    double potential_charge_abs_scale = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double old_phi_average =
            before.phi[static_cast<size_t>(sg.nghost + ix)] + sg.dx *
            (before.Ex_face[static_cast<size_t>(ix + 1)] -
             before.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const double new_phi_average =
            after.phi[static_cast<size_t>(sg.nghost + ix)] + sg.dx *
            (after.Ex_face[static_cast<size_t>(ix + 1)] -
             after.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const double phi_average = 0.5 * (old_phi_average + new_phi_average);
        max_pairing_phi = std::max(max_pairing_phi, std::fabs(phi_average));
        potential_charge_abs_scale +=
            std::fabs(phi_average * rho_delta[static_cast<size_t>(sg.nghost + ix)]) *
            sg.dx;
    }
    result.potential_charge_abs_scale = potential_charge_abs_scale;
    // A rounded M_new-M_old enters rho_delta and is then multiplied by the
    // production cell-average potential.  Include that independently
    // bounded contribution in the G/G* audit floor; otherwise a large-M
    // manufactured case can fail even when the flux identity is exact.
    const double mass_roundoff_work_bound =
        std::fabs(Const::qe) * max_pairing_phi *
        result.audit.mass_roundoff_bound;
    result.g_gstar_roundoff_tolerance =
        eps * (result.g_gstar_scale + potential_charge_abs_scale) +
        mass_roundoff_work_bound;
    const double combined_roundoff_tolerance =
        result.kinetic_roundoff_tolerance +
        result.poisson_roundoff_tolerance;
    result.passed = result.local_global_flux_pass &&
        result.u_energy_adjoint_pass &&
        result.x_force_charge_adjoint_pass &&
        result.audit.finite && pairing_ok &&
        result.audit.mass_residual <=
            std::max(result.audit.mass_roundoff_bound, 1.0e-300) &&
        std::fabs(result.audit.kinetic_work_residual) <=
            result.kinetic_roundoff_tolerance &&
        std::fabs(result.audit.poisson_work_residual) <=
            result.poisson_roundoff_tolerance &&
        std::fabs(result.audit.g_gstar_residual) <=
            result.g_gstar_roundoff_tolerance &&
        std::fabs(result.audit.combined_energy_residual) <=
            combined_roundoff_tolerance &&
        result.audit.cell_volume_residual <= 1.0e-14 &&
        result.audit.hamiltonian_velocity_residual /
            (Const::me * Const::c * Const::c) <= 1.0e-14 &&
        result.audit.u_boundary_flux == 0.0;
    return result;
}

void write_a2_case(std::ostream& out, const char* tag,
                   const A2CaseResult& r)
{
    out << "a2_" << tag << "_ran=" << (r.ran ? 1 : 0) << "\n"
        << "a2_" << tag << "_finite=" << (r.finite ? 1 : 0) << "\n"
        << "a2_" << tag << "_poisson_scalar_identity_pass="
        << (r.poisson_scalar_identity_pass ? 1 : 0) << "\n"
        << "a2_" << tag << "_charge_residual_projection_pass="
        << (r.charge_residual_projection_pass ? 1 : 0) << "\n"
        << "a2_" << tag << "_poisson_current_prediction_pass="
        << (r.poisson_current_prediction_pass ? 1 : 0) << "\n"
        << "a2_" << tag << "_nontrivial_wc_pass="
        << (r.nontrivial_wc_pass ? 1 : 0) << "\n"
        << "a2_" << tag << "_w_c=" << r.w_c << "\n"
        << "a2_" << tag << "_r_p=" << r.r_p << "\n"
        << "a2_" << tag << "_r_pj=" << r.r_pj << "\n"
        << "a2_" << tag << "_prediction_error=" << r.prediction_error
        << "\n"
        << "a2_" << tag << "_tau_c=" << r.tau_c << "\n"
        << "a2_" << tag << "_tau_a=" << r.tau_a << "\n"
        << "a2_" << tag << "_mismatch_linf=" << r.mismatch_linf << "\n"
        << "a2_" << tag << "_r_q_linf=" << r.r_q_linf << "\n"
        << "a2_" << tag << "_r_c_linf=" << r.r_c_linf << "\n"
        << "a2_" << tag << "_u_boundary_linf=" << r.u_boundary_linf
        << "\n"
        << "a2_" << tag << "_total_charge_change="
        << r.total_charge_change << "\n"
        << "a2_" << tag << "_pass=" << (r.passed ? 1 : 0) << "\n";
}

void write_result(std::ostream& out,
                  const std::vector<ScenarioResult>& results,
                  const UFluxGeometryResult& geometry,
                  const PeriodicSeamAdjointResult* seam,
                  const A2PairResult* a2,
                  const As1IncrementalRhoResult* as1)
{
    double max_mass = 0.0;
    double max_kinetic = 0.0;
    double max_poisson = 0.0;
    double max_combined = 0.0;
    double max_pairing = 0.0;
    bool pass = !results.empty() && geometry.passed &&
        (!seam || !seam->selected || seam->passed) &&
        (!a2 || !a2->selected || a2->passed) &&
        (!as1 || !as1->selected || as1->passed);
    bool j0_c = !results.empty();
    bool j0_d = !results.empty();
    bool j0_e = !results.empty();
    bool j0_f = !results.empty();
    for (size_t i = 0; i < results.size(); ++i) {
        const ScenarioResult& r = results[i];
        pass = pass && r.passed;
        j0_c = j0_c && r.local_global_flux_pass;
        j0_d = j0_d && r.u_energy_adjoint_pass;
        j0_e = j0_e && r.x_force_charge_adjoint_pass;
        j0_f = j0_f && r.audit.finite &&
            std::fabs(r.audit.poisson_work_residual) <=
                r.poisson_roundoff_tolerance &&
            std::fabs(r.audit.g_gstar_residual) <=
                r.g_gstar_roundoff_tolerance &&
            std::fabs(r.audit.combined_energy_residual) <=
                r.kinetic_roundoff_tolerance +
                r.poisson_roundoff_tolerance &&
            r.audit.cell_volume_residual <= 1.0e-14 &&
            r.audit.u_boundary_flux == 0.0;
        max_mass = std::max(max_mass, r.audit.mass_residual);
        max_kinetic = std::max(max_kinetic,
                               std::fabs(r.audit.kinetic_work_residual));
        max_poisson = std::max(max_poisson,
                               std::fabs(r.audit.poisson_work_residual));
        max_combined = std::max(max_combined,
                                std::fabs(r.audit.combined_energy_residual));
        max_pairing = std::max(max_pairing,
                               std::fabs(r.audit.g_gstar_residual));
        out << "case_" << r.name << "_pass=" << (r.passed ? 1 : 0) << "\n"
            << "case_" << r.name << "_mass_residual="
            << r.audit.mass_residual << "\n"
            << "case_" << r.name << "_mass_roundoff_bound="
            << r.audit.mass_roundoff_bound << "\n"
            << "case_" << r.name << "_kinetic_work_residual="
            << r.audit.kinetic_work_residual << "\n"
            << "case_" << r.name << "_kinetic_absolute_work_scale="
            << r.audit.kinetic_absolute_work_scale << "\n"
            << "case_" << r.name << "_poisson_work_residual="
            << r.audit.poisson_work_residual << "\n"
            << "case_" << r.name << "_combined_energy_residual="
            << r.audit.combined_energy_residual << "\n"
            << "case_" << r.name << "_g_gstar_residual="
            << r.audit.g_gstar_residual << "\n"
            << "case_" << r.name << "_cell_volume_residual="
            << r.audit.cell_volume_residual << "\n"
            << "case_" << r.name << "_hamiltonian_velocity_residual="
            << r.audit.hamiltonian_velocity_residual << "\n"
            << "case_" << r.name << "_u_boundary_flux="
            << r.audit.u_boundary_flux << "\n"
            << "case_" << r.name << "_poisson_scale="
            << r.poisson_scale << "\n"
            << "case_" << r.name << "_poisson_roundoff_tolerance="
            << r.poisson_roundoff_tolerance << "\n"
            << "case_" << r.name << "_g_gstar_roundoff_tolerance="
            << r.g_gstar_roundoff_tolerance << "\n"
            << "case_" << r.name << "_local_global_flux_pass="
            << (r.local_global_flux_pass ? 1 : 0) << "\n"
            << "case_" << r.name << "_local_global_flux_relative_error="
            << r.local_global_flux_relative_error << "\n"
            << "case_" << r.name << "_u_energy_adjoint_pass="
            << (r.u_energy_adjoint_pass ? 1 : 0) << "\n"
            << "case_" << r.name << "_u_energy_adjoint_relative_error="
            << r.u_energy_adjoint_relative_error << "\n"
            << "case_" << r.name << "_u_energy_velocity_symmetry_error="
            << r.u_energy_velocity_symmetry_error << "\n"
            << "case_" << r.name << "_u_energy_velocity_max_abs="
            << r.u_energy_velocity_max_abs << "\n"
            << "case_" << r.name << "_x_force_charge_adjoint_pass="
            << (r.x_force_charge_adjoint_pass ? 1 : 0) << "\n"
            << "case_" << r.name << "_x_force_charge_adjoint_relative_error="
            << r.x_force_charge_adjoint_relative_error << "\n"
            << "case_" << r.name << "_potential_charge_abs_scale="
            << r.potential_charge_abs_scale << "\n";
    }
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n"
        << "j0_a_cell_mass_pass="
        << ((geometry.cell_mass_relative_error <=
             4096.0 * std::numeric_limits<double>::epsilon()) ? 1 : 0)
        << "\n"
        << "j0_b_u_flux_geometry_pass=" << (geometry.passed ? 1 : 0) << "\n"
        << "j0_c_midpoint_consistency_pass=" << (j0_c ? 1 : 0) << "\n"
        << "j0_d_u_work_force_adjoint_pass=" << (j0_d ? 1 : 0) << "\n"
        << "j0_e_x_current_force_adjoint_pass=" << (j0_e ? 1 : 0) << "\n"
        << "j0_f_poisson_pairing_pass=" << (j0_f ? 1 : 0) << "\n"
        << "j0_all_six_pass=" << (pass && j0_c && j0_d && j0_e && j0_f ? 1 : 0)
        << "\n"
        << "u_flux_geometry_selected=" << (geometry.selected ? 1 : 0) << "\n"
        << "u_flux_geometry_pass=" << (geometry.passed ? 1 : 0) << "\n"
        << "u_flux_geometry_dx1=" << geometry.dx1 << "\n"
        << "u_flux_geometry_dx2=" << geometry.dx2 << "\n"
        << "u_flux_geometry_ring_k1=" << geometry.ring_k1 << "\n"
        << "u_flux_geometry_ring_k2=" << geometry.ring_k2 << "\n"
        << "u_flux_geometry_relative_error_dx1="
        << geometry.relative_error_dx1 << "\n"
        << "u_flux_geometry_relative_error_dx2="
        << geometry.relative_error_dx2 << "\n"
        << "u_flux_geometry_cell_mass_relative_error="
        << geometry.cell_mass_relative_error << "\n"
        << "mass_residual=" << max_mass << "\n"
        << "kinetic_work_residual=" << max_kinetic << "\n"
        << "poisson_work_residual=" << max_poisson << "\n"
        << "combined_energy_residual=" << max_combined << "\n"
        << "g_gstar_residual=" << max_pairing << "\n"
        << "endpoint_face_weight_left=0.5\n"
        << "endpoint_face_weight_right=0.5\n"
        << "endpoint_face_weight_definition=dx_over_2\n"
        << "field_pairing_definition=production_OpenElectrostaticSolver_G_Gstar\n"
        << "x_flux_units=cell_integrated_mass_per_area_per_second\n"
        << "u_flux_units=cell_integrated_mass_per_area_per_second\n"
        << "charge_current_units=A_per_m2\n"
        << "poisson_rho_delta_source=x_charge_flux_divergence\n"
        << "u_boundary_condition=zero_inflow_zero_outflow_ledger\n"
        << "x_trace=periodic_center_trace\n"
        << "combined_energy_definition=kinetic_plus_poisson_identity_residual\n";
    if (seam && seam->selected) {
        out << "j0_e2_periodic_seam_weighted_adjoint_pass="
            << (seam->passed ? 1 : 0) << "\n"
            << "j0_e2_helper_ok=" << (seam->helper_ok ? 1 : 0) << "\n"
            << "j0_e2_endpoint_asymmetry_pass="
            << (seam->endpoint_asymmetry_pass ? 1 : 0) << "\n"
            << "j0_e2_weighted_adjoint_pass="
            << (seam->weighted_adjoint_pass ? 1 : 0) << "\n"
            << "j0_e2_prediction_nonzero_pass="
            << (seam->prediction_nonzero_pass ? 1 : 0) << "\n"
            << "j0_e2_prediction_identity_pass="
            << (seam->prediction_identity_pass ? 1 : 0) << "\n"
            << "j0_e2_weighted_adjoint_relative_error="
            << seam->weighted_adjoint_relative_error << "\n"
            << "j0_e2_w_f_adjoint=" << seam->w_f_adjoint << "\n"
            << "j0_e2_w_j=" << seam->w_j << "\n"
            << "j0_e2_w_f_naive=" << seam->w_f_naive << "\n"
            << "j0_e2_old_naive_mismatch=" << seam->old_naive_mismatch << "\n"
            << "j0_e2_seam_predicted_residual="
            << seam->seam_predicted_residual << "\n"
            << "j0_e2_prediction_error=" << seam->prediction_error << "\n"
            << "j0_e2_pairing_face_left=" << seam->pairing_face_left << "\n"
            << "j0_e2_pairing_face_right=" << seam->pairing_face_right << "\n"
            << "j0_e2_force_current_first_cell="
            << seam->force_current_first_cell << "\n"
             << "j0_e2_force_current_last_cell="
             << seam->force_current_last_cell << "\n";
    }
    if (a2 && a2->selected) {
        write_a2_case(out, "zero_residual", a2->zero);
        write_a2_case(out, "injected_residual", a2->injected);
        out << "a2_all_pass="
            << ((a2->zero.passed && a2->injected.passed) ? 1 : 0)
            << "\n";
    }
    if (as1 && as1->selected) {
        out << "as1_incremental_positive_pass="
            << (as1->positive_pass ? 1 : 0) << "\n"
            << "as1_incremental_negative_pass="
            << (as1->negative_pass ? 1 : 0) << "\n"
            << "as1_incremental_near_neutral_small_pass="
            << (as1->near_neutral_small_pass ? 1 : 0) << "\n"
            << "as1_absolute_form_degrades="
            << (as1->absolute_form_degrades ? 1 : 0) << "\n"
            << "as1_incremental_positive_error="
            << as1->positive_incremental_error << "\n"
            << "as1_absolute_positive_error="
            << as1->positive_absolute_error << "\n"
            << "as1_incremental_negative_error="
            << as1->negative_incremental_error << "\n"
            << "as1_incremental_near_neutral_error="
            << as1->near_neutral_incremental_error << "\n"
            << "as1_absolute_near_neutral_error="
            << as1->near_neutral_absolute_error << "\n"
            << "as1_roundoff_bound=" << as1->roundoff_bound << "\n"
            << "as1_incremental_rho_all_pass="
            << (as1->passed ? 1 : 0) << "\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string result_path;
    std::string test_case = "all";
    bool parsed = true;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) result_path = argv[++i];
        else parsed = false;
    }
    if (!parsed || (test_case != "all" && test_case != "positive-field" &&
                    test_case != "negative-field" && test_case != "single-velocity" &&
                    test_case != "two-velocity")) {
        if (rank == 0)
            std::cerr << "usage: joint_phase_space_midpoint_unit_test "
                         "--case all|positive-field|negative-field|single-velocity|two-velocity "
                         "[--result path]\n";
        MPI_Finalize();
        return 2;
    }
    if (size != 1) {
        if (rank == 0)
            std::cerr << "J0 unit test requires exactly one MPI rank\n";
        MPI_Finalize();
        return 2;
    }

    std::vector<Scenario> scenarios;
    if (test_case == "all" || test_case == "positive-field")
        scenarios.push_back({"positive_field", 2.0e8, 0});
    if (test_case == "all" || test_case == "negative-field")
        scenarios.push_back({"negative_field", -2.0e8, 0});
    if (test_case == "all" || test_case == "single-velocity")
        scenarios.push_back({"single_velocity", 2.0e8, 1});
    if (test_case == "all" || test_case == "two-velocity")
        scenarios.push_back({"two_velocity", -2.0e8, 2});

    std::vector<ScenarioResult> results;
    for (size_t i = 0; i < scenarios.size(); ++i)
        results.push_back(run_scenario(scenarios[i], rank, size));

    UFluxGeometryResult geometry;
    if (test_case == "all") geometry = run_u_flux_geometry_test(rank, size);

    PeriodicSeamAdjointResult seam_result;
    const bool seam_selected = test_case == "all";
    if (seam_selected)
        seam_result = run_periodic_seam_adjoint_test(rank, size);

    A2PairResult a2_result;
    const bool a2_selected = test_case == "all";
    if (a2_selected) {
        a2_result.selected = true;
        a2_result.zero =
            run_a2_case(/*inject_residual=*/false, rank, size);
        a2_result.injected =
            run_a2_case(/*inject_residual=*/true, rank, size);
        a2_result.passed =
            a2_result.zero.passed && a2_result.injected.passed;
    }

    As1IncrementalRhoResult as1_result;
    const bool as1_selected = test_case == "all";
    if (as1_selected)
        as1_result = run_as1_incremental_rho_test(rank, size);

    bool pass = true;
    for (size_t i = 0; i < results.size(); ++i) pass = pass && results[i].passed;
    if (!geometry.passed) pass = false;
    if (seam_selected && !seam_result.passed) pass = false;
    if (a2_selected && !a2_result.passed) pass = false;
    if (as1_selected && !as1_result.passed) pass = false;
    if (rank == 0) {
        std::cout << std::setprecision(17);
        write_result(std::cout, results, geometry,
                     seam_selected ? &seam_result : 0,
                     a2_selected ? &a2_result : 0,
                     as1_selected ? &as1_result : 0);
        if (!result_path.empty()) {
            std::ofstream out(result_path.c_str(), std::ios::trunc);
            if (!out) pass = false;
            else {
                out << std::setprecision(17);
                write_result(out, results, geometry,
                             seam_selected ? &seam_result : 0,
                             a2_selected ? &a2_result : 0,
                             as1_selected ? &as1_result : 0);
            }
        }
    }
    int local_pass = pass ? 1 : 0;
    int global_pass = 0;
    MPI_Allreduce(&local_pass, &global_pass, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Finalize();
    return global_pass ? 0 : 1;
}
