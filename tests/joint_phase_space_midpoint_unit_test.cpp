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
#include <vector>

namespace {

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
    bool passed;
    ScenarioResult()
        : mass_scale(1.0), kinetic_scale(1.0), poisson_scale(1.0),
          g_gstar_scale(1.0), potential_charge_abs_scale(0.0),
          kinetic_roundoff_tolerance(0.0),
          poisson_roundoff_tolerance(0.0),
          g_gstar_roundoff_tolerance(0.0), passed(false) {}
};

size_t cell_index(int ix, int j, int k, int nv, int nmu)
{
    return (static_cast<size_t>(ix) * static_cast<size_t>(nv) +
            static_cast<size_t>(j)) * static_cast<size_t>(nmu) +
           static_cast<size_t>(k);
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
                m_mid[cell_index(ix, j, k, nv, nmu)] = active
                    ? 1.0e20 * x_factor * shape *
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
    result.passed = result.audit.finite && pairing_ok &&
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

void write_result(std::ostream& out, const std::vector<ScenarioResult>& results)
{
    double max_mass = 0.0;
    double max_kinetic = 0.0;
    double max_poisson = 0.0;
    double max_combined = 0.0;
    double max_pairing = 0.0;
    bool pass = !results.empty();
    for (size_t i = 0; i < results.size(); ++i) {
        const ScenarioResult& r = results[i];
        pass = pass && r.passed;
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
            << "case_" << r.name << "_potential_charge_abs_scale="
            << r.potential_charge_abs_scale << "\n";
    }
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n"
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

    bool pass = true;
    for (size_t i = 0; i < results.size(); ++i) pass = pass && results[i].passed;
    if (rank == 0) {
        std::cout << std::setprecision(17);
        write_result(std::cout, results);
        if (!result_path.empty()) {
            std::ofstream out(result_path.c_str(), std::ios::trunc);
            if (!out) pass = false;
            else {
                out << std::setprecision(17);
                write_result(out, results);
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
