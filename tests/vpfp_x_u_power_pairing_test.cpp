// P3.0: read-only production Strang x/u power-pairing audit.
// The tested time layer is Tx(dt/2) -> Tu(E_mid,dt) -> Tx(dt/2).

#include "conservative_ppm_remap.h"
#include "field_particle_power_audit.h"
#include "open_boundary.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

const double kPairingTolerance = 1.0e-10;
const double kPairingEnergyFloor = 1.0e-30;

struct Args {
    std::string test_case;
    std::string result;
    XTransportVelocityMode x_transport_velocity_mode;
    Args() : test_case("all"),
              x_transport_velocity_mode(
                  XTransportVelocityMode::ANALYTIC_CELL_CENTER) {}
};

struct CaseSpec {
    const char* name;
    bool gradient;
    bool single_velocity;
    bool open_endpoint;
    bool baseline_case;
    bool limiter_case;
    double e_face;
};

struct CellPowerAudit {
    int global_cell;
    double e_cell;
    double nxhalf;
    double delta_ke_u;
    double j_x;
    double u_work_current;
    double x_flux_current;
    double signed_current_difference;
    double signed_power_difference;
    double v_u_power;
    double v_x_flux;
    double power_current_identity_error;
    bool excluded_small_field;
    bool excluded_small_number;
};

struct CaseMetrics {
    std::string name;
    double target_work;
    double actual_u_work;
    double absolute_error;
    double pairing_scale;
    double relative_error;
    double pairing_tolerance;
    double continuity_residual;
    double continuity_scale;
    double continuity_tolerance;
    bool finite;
    bool strang_sequence_valid;
    bool final_flux_current_valid;
    bool harness_pass;
    bool physical_pairing_pass;
    std::vector<CellPowerAudit> cell_audit;
    double max_power_current_identity_error;
    double max_power_current_identity_relative;
    int included_cell_count;
    int excluded_small_field_count;
    int excluded_small_number_count;
    bool power_current_identity_pass;
    double x1_linear_fraction;
    double x2_linear_fraction;
    double x1_constant_fraction;
    double x2_constant_fraction;

    CaseMetrics()
        : target_work(0.0), actual_u_work(0.0), absolute_error(0.0),
          pairing_scale(kPairingEnergyFloor), relative_error(0.0),
          pairing_tolerance(kPairingTolerance), continuity_residual(0.0),
          continuity_scale(1.0), continuity_tolerance(0.0), finite(true),
          strang_sequence_valid(false), final_flux_current_valid(false),
          harness_pass(false), physical_pairing_pass(false),
          max_power_current_identity_error(0.0),
          max_power_current_identity_relative(0.0),
          included_cell_count(0), excluded_small_field_count(0),
          excluded_small_number_count(0), power_current_identity_pass(false),
          x1_linear_fraction(0.0), x2_linear_fraction(0.0),
          x1_constant_fraction(0.0), x2_constant_fraction(0.0) {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) args.test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) args.result = argv[++i];
        else if (arg == "--x-transport-velocity-mode" && i + 1 < argc) {
            const std::string value(argv[++i]);
            if (value == "analytic-cell-center")
                args.x_transport_velocity_mode =
                    XTransportVelocityMode::ANALYTIC_CELL_CENTER;
            else if (value == "energy-conjugate")
                args.x_transport_velocity_mode =
                    XTransportVelocityMode::ENERGY_CONJUGATE_CELL;
            else return false;
        }
        else return false;
    }
    return args.test_case == "all" || args.test_case == "uniform" ||
           args.test_case == "gradient" || args.test_case == "mpi-shared-face";
}

size_t index3(int ix, int j, int k)
{
    return (static_cast<size_t>(ix) * Param::Nv + static_cast<size_t>(j)) *
           Param::Nmu + static_cast<size_t>(k);
}

void fill_species(Species& species, const SpatialGrid& grid,
                  const CylindricalVelocityGrid& velocity_grid,
                  bool gradient, bool single_velocity)
{
    const double sigma = 0.12;
    const double u_drift = 0.35;
    const double amplitude = 1.0e23;
    int single_j = 0;
    double single_distance = std::numeric_limits<double>::infinity();
    for (int j = 0; j < Param::Nv; ++j) {
        const double distance = std::fabs(
            velocity_grid.upar_cells[static_cast<size_t>(j)] - u_drift);
        if (distance < single_distance) {
            single_distance = distance;
            single_j = j;
        }
    }

    for (int ix = 0; ix < grid.nx_total; ++ix) {
        const int global_ix = grid.global_cell(ix);
        const double xi = (static_cast<double>(global_ix) + 0.5) /
                          static_cast<double>(grid.nx_global);
        const double x_factor = gradient
            ? 1.0 + 0.35 * std::sin(2.0 * Const::pi * xi) : 1.0;
        for (int j = 0; j < Param::Nv; ++j) {
            const double u = velocity_grid.upar_cells[static_cast<size_t>(j)];
            const double du = velocity_grid.upar_widths[static_cast<size_t>(j)];
            for (int k = 0; k < Param::Nmu; ++k) {
                const double up = velocity_grid.uperp_cells[static_cast<size_t>(k)];
                const double ring = velocity_grid.uperp_ring_areas[static_cast<size_t>(k)];
                const bool active = !single_velocity || j == single_j;
                const double f = active
                    ? amplitude * x_factor * std::exp(
                        -((u - u_drift) * (u - u_drift) + up * up) /
                         (2.0 * sigma * sigma))
                    : 0.0;
                species.f[index3(ix, j, k)] = f * grid.dx * du * ring;
            }
        }
    }
}

CaseMetrics run_case(const CaseSpec& spec, int rank, int mpi_size,
                     XTransportVelocityMode velocity_mode)
{
    CaseMetrics metrics;
    metrics.name = spec.name;

    SpatialGrid grid;
    grid.init_with_domain(rank, mpi_size, 64, 40.0 * Const::micro);
    CylindricalVelocityGrid velocity_grid;
    velocity_grid.init(Param::momentum_umax);
    ConservativePpmRemap remap;
    remap.init(grid, velocity_grid);
    remap.set_x_transport_velocity_mode(velocity_mode);

    Species input;
    input.init("bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
               Const::me, Param::dens, Param::temperature_e, false, grid);
    input.f.assign(static_cast<size_t>(grid.nx_total) * Param::Nvmu, 0.0);
    fill_species(input, grid, velocity_grid, spec.gradient,
                 spec.single_velocity);

    OpenBackgroundBoundaryConfig boundary_config;
    boundary_config.left_type = spec.open_endpoint
        ? BackgroundXBoundaryType::ABSORBING
        : BackgroundXBoundaryType::RESERVOIR;
    boundary_config.right_type = spec.open_endpoint
        ? BackgroundXBoundaryType::ABSORBING
        : BackgroundXBoundaryType::RESERVOIR;
    const OpenBackgroundBoundary boundary(boundary_config);

    Species x_half;
    x_half.init("bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, grid);
    x_half.f.assign(input.f.size(), 0.0);
    Species u_output;
    u_output.init("bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                  Const::me, Param::dens, Param::temperature_e, false, grid);
    u_output.f.assign(input.f.size(), 0.0);
    Species final_output;
    final_output.init("bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                      Const::me, Param::dens, Param::temperature_e, false, grid);
    final_output.f.assign(input.f.size(), 0.0);

    XFaceTransportAudit x_audit_1;
    x_audit_1.enabled = true;
    x_audit_1.init(grid.nx_local);
    XFaceTransportAudit x_audit_2;
    x_audit_2.enabled = true;
    x_audit_2.init(grid.nx_local);

    const double dt = 1.0e-18;
    const double half_dt = 0.5 * dt;
    const RemapDiagnostics x_diag_1 = remap.advect_x(
        input, x_half, half_dt, 0.0, boundary, rank, mpi_size, &x_audit_1);

    EMFields fields;
    fields.init(grid);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), spec.e_face);
    std::fill(fields.Ex.begin(), fields.Ex.end(), spec.e_face);

    std::vector<double> u_delta_ke;
    const RemapDiagnostics u_diag = remap.advect_u_parallel(
        x_half, u_output, fields, dt, half_dt, NULL, NULL, 4, &u_delta_ke);
    const RemapDiagnostics x_diag_2 = remap.advect_x(
        u_output, final_output, half_dt, dt, boundary, rank, mpi_size,
        &x_audit_2);

    const int nxl = grid.nx_local;
    double local_continuity = 0.0;
    double local_continuity_scale = 1.0;
    double local_target_work = 0.0;
    double local_actual_work = 0.0;
    for (int i = 0; i < nxl; ++i) {
        double before = 0.0;
        double after = 0.0;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                before += input.f[index3(grid.nghost + i, j, k)];
                after += final_output.f[index3(grid.nghost + i, j, k)];
            }
        }
        const double left = x_audit_1.bulk_number_swept_face[
                                static_cast<size_t>(i)] +
                            x_audit_2.bulk_number_swept_face[
                                static_cast<size_t>(i)];
        const double right = x_audit_1.bulk_number_swept_face[
                                 static_cast<size_t>(i + 1)] +
                             x_audit_2.bulk_number_swept_face[
                                 static_cast<size_t>(i + 1)];
        local_continuity = std::max(
            local_continuity, std::fabs(after - before - left + right));
        local_continuity_scale = std::max(
            local_continuity_scale, std::fabs(before) + std::fabs(after) +
                std::fabs(left) + std::fabs(right));

        const double j_left = input.charge * left / dt;
        const double j_right = input.charge * right / dt;
        const double j_cell = 0.5 * (j_left + j_right);
        local_target_work += dt * grid.dx * spec.e_face * j_cell;
        if (static_cast<size_t>(i) < u_delta_ke.size())
            local_actual_work += u_delta_ke[static_cast<size_t>(i)];
        else
            metrics.finite = false;
    }

    MPI_Allreduce(&local_continuity, &metrics.continuity_residual, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_continuity_scale, &metrics.continuity_scale, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_target_work, &metrics.target_work, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_actual_work, &metrics.actual_u_work, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const int local_finite = metrics.finite && x_diag_1.finite &&
        u_diag.finite && x_diag_2.finite &&
        std::isfinite(local_target_work) && std::isfinite(local_actual_work);
    int global_finite = 0;
    MPI_Allreduce(&local_finite, &global_finite, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    metrics.finite = global_finite != 0;
    metrics.strang_sequence_valid = metrics.finite;
    const int local_flux_valid =
        x_audit_1.bulk_number_swept_face.size() ==
            static_cast<size_t>(nxl) + 1 &&
        x_audit_2.bulk_number_swept_face.size() ==
            static_cast<size_t>(nxl) + 1;
    int global_flux_valid = local_flux_valid ? 1 : 0;
    MPI_Allreduce(&local_flux_valid, &global_flux_valid, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    metrics.final_flux_current_valid = global_flux_valid != 0;
    metrics.x1_linear_fraction = x_diag_1.limited_fraction;
    metrics.x2_linear_fraction = x_diag_2.limited_fraction;
    metrics.x1_constant_fraction = x_diag_1.constant_fraction;
    metrics.x2_constant_fraction = x_diag_2.constant_fraction;

    // P3-V.1R: signed, same-cell power/current audit.  Only cells whose field
    // is above the explicit relative floor enter the velocity auxiliary
    // statistics; every excluded cell is counted and reported below.
    double local_max_abs_e = 0.0;
    for (size_t f = 0; f < fields.Ex_face.size(); ++f)
        local_max_abs_e = std::max(local_max_abs_e,
                                   std::fabs(fields.Ex_face[f]));
    double global_max_abs_e = 0.0;
    MPI_Allreduce(&local_max_abs_e, &global_max_abs_e, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    const double field_relative_floor = 1.0e-12;
    const double field_threshold = field_relative_floor * global_max_abs_e;
    const double number_floor = 1.0e-300;
    double local_identity_error = 0.0;
    double local_identity_relative = 0.0;
    int local_included = 0;
    int local_excluded_field = 0;
    int local_excluded_number = 0;
    for (int i = 0; i < nxl; ++i) {
        CellPowerAudit cell = {};
        cell.global_cell = grid.ix_start + i;
        cell.e_cell = 0.5 * (fields.Ex_face[static_cast<size_t>(i)] +
                             fields.Ex_face[static_cast<size_t>(i + 1)]);
        cell.nxhalf = 0.0;
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                cell.nxhalf += x_half.f[index3(grid.nghost + i, j, k)];
        const double q_left =
            x_audit_1.bulk_number_swept_face[static_cast<size_t>(i)] +
            x_audit_2.bulk_number_swept_face[static_cast<size_t>(i)];
        const double q_right =
            x_audit_1.bulk_number_swept_face[static_cast<size_t>(i + 1)] +
            x_audit_2.bulk_number_swept_face[static_cast<size_t>(i + 1)];
        const double j_left = input.charge * q_left / dt;
        const double j_right = input.charge * q_right / dt;
        cell.delta_ke_u = i < static_cast<int>(u_delta_ke.size())
            ? u_delta_ke[static_cast<size_t>(i)] : 0.0;
        cell.j_x = 0.5 * (j_left + j_right);
        cell.x_flux_current = cell.j_x;
        cell.excluded_small_field =
            std::fabs(cell.e_cell) <= field_threshold;
        cell.excluded_small_number =
            std::fabs(cell.nxhalf) <= number_floor;
        if (cell.excluded_small_field) ++local_excluded_field;
        if (cell.excluded_small_number) ++local_excluded_number;
        if (!cell.excluded_small_field && !cell.excluded_small_number) {
            ++local_included;
            // Recompute the P3-V.1R identity in extended precision.  The
            // identity is an audit change of variables, so a double-valued
            // divide/multiply/subtract chain must not turn harmless
            // cancellation into a false harness failure.
            const long double dt_ld = static_cast<long double>(dt);
            const long double dx_ld = static_cast<long double>(grid.dx);
            const long double e_ld = static_cast<long double>(cell.e_cell);
            const long double dke_ld =
                static_cast<long double>(cell.delta_ke_u);
            const long double q_ld = static_cast<long double>(input.charge);
            const long double n_ld = static_cast<long double>(cell.nxhalf);
            const long double j_ld = static_cast<long double>(cell.j_x);
            const long double u_work_current_ld =
                dke_ld / (dt_ld * dx_ld * e_ld);
            const long double x_flux_current_ld =
                dx_ld * j_ld / (q_ld * n_ld);
            const long double signed_power_difference_ld =
                dke_ld - dt_ld * dx_ld * e_ld * j_ld;
            const long double v_u_power_ld =
                dke_ld / (dt_ld * q_ld * e_ld * n_ld);
            const long double v_x_flux_ld =
                dx_ld * j_ld / (q_ld * n_ld);
            const long double rhs_ld =
                dt_ld * q_ld * e_ld * n_ld *
                (v_u_power_ld - v_x_flux_ld);
            const long double identity_error_ld =
                signed_power_difference_ld - rhs_ld;
            const long double identity_scale_ld = std::max(
                static_cast<long double>(kPairingEnergyFloor),
                std::max(std::fabs(signed_power_difference_ld),
                         std::fabs(rhs_ld)));

            cell.u_work_current = static_cast<double>(u_work_current_ld);
            cell.x_flux_current = cell.j_x;
            cell.signed_current_difference =
                cell.u_work_current - cell.x_flux_current;
            cell.signed_power_difference =
                static_cast<double>(signed_power_difference_ld);
            cell.v_u_power = static_cast<double>(v_u_power_ld);
            cell.v_x_flux = static_cast<double>(v_x_flux_ld);
            cell.power_current_identity_error =
                static_cast<double>(identity_error_ld);
            local_identity_error = std::max(
                local_identity_error,
                static_cast<double>(std::fabs(identity_error_ld)));
            local_identity_relative = std::max(
                local_identity_relative,
                static_cast<double>(std::fabs(identity_error_ld) /
                                    identity_scale_ld));
        }
        metrics.cell_audit.push_back(cell);
    }
    MPI_Allreduce(&local_identity_error,
                  &metrics.max_power_current_identity_error, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_identity_relative,
                  &metrics.max_power_current_identity_relative, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_included, &metrics.included_cell_count, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_excluded_field, &metrics.excluded_small_field_count, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_excluded_number, &metrics.excluded_small_number_count, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    metrics.power_current_identity_pass =
        metrics.max_power_current_identity_relative <=
        4096.0 * std::numeric_limits<double>::epsilon();

    metrics.absolute_error = std::fabs(
        metrics.actual_u_work - metrics.target_work);
    metrics.pairing_scale = std::max(
        kPairingEnergyFloor,
        std::max(std::fabs(metrics.actual_u_work),
                 std::fabs(metrics.target_work)));
    metrics.relative_error = metrics.absolute_error / metrics.pairing_scale;
    metrics.continuity_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, metrics.continuity_scale);
    metrics.harness_pass = metrics.finite && metrics.strang_sequence_valid &&
        metrics.final_flux_current_valid &&
        metrics.continuity_residual <= metrics.continuity_tolerance &&
        metrics.power_current_identity_pass;
    metrics.physical_pairing_pass = metrics.harness_pass;
    return metrics;
}

void write_case(std::ofstream& out, const CaseMetrics& m)
{
    const std::string p = "case_" + m.name + "_";
    out << p << "target_work=" << m.target_work << "\n"
        << p << "actual_u_work=" << m.actual_u_work << "\n"
        << p << "absolute_error=" << m.absolute_error << "\n"
        << p << "pairing_scale=" << m.pairing_scale << "\n"
        << p << "relative_error=" << m.relative_error << "\n"
        << p << "pairing_tolerance=" << m.pairing_tolerance << "\n"
        << p << "pairing_energy_floor=" << kPairingEnergyFloor << "\n"
        << p << "continuity_residual=" << m.continuity_residual << "\n"
        << p << "continuity_scale=" << m.continuity_scale << "\n"
        << p << "continuity_tolerance=" << m.continuity_tolerance << "\n"
        << p << "finite=" << (m.finite ? 1 : 0) << "\n"
        << p << "strang_sequence_valid="
        << (m.strang_sequence_valid ? 1 : 0) << "\n"
        << p << "final_flux_current_valid="
        << (m.final_flux_current_valid ? 1 : 0) << "\n"
        << p << "harness_pass=" << (m.harness_pass ? 1 : 0) << "\n"
        << p << "physical_pairing_pass="
        << (m.physical_pairing_pass ? 1 : 0) << "\n"
        << p << "max_power_current_identity_error="
        << m.max_power_current_identity_error << "\n"
        << p << "max_power_current_identity_relative="
        << m.max_power_current_identity_relative << "\n"
        << p << "included_cell_count=" << m.included_cell_count << "\n"
        << p << "excluded_small_field_count="
        << m.excluded_small_field_count << "\n"
        << p << "excluded_small_number_count="
        << m.excluded_small_number_count << "\n"
        << p << "power_current_identity_pass="
        << (m.power_current_identity_pass ? 1 : 0) << "\n"
        << p << "x1_limiter_active_fraction=" << m.x1_linear_fraction << "\n"
        << p << "x2_limiter_active_fraction=" << m.x2_linear_fraction << "\n"
        << p << "x1_constant_fraction=" << m.x1_constant_fraction << "\n"
        << p << "x2_constant_fraction=" << m.x2_constant_fraction << "\n";
}

void write_cell_audit(std::ofstream& out, const CaseMetrics& m)
{
    if (out.tellp() == std::streampos(0)) {
        out << "# columns=case global_cell E_cell N_xhalf delta_ke_u J_x "
               "u_work_current_cell x_flux_current_cell signed_current_difference "
               "signed_power_difference v_u_power v_x_flux "
               "power_current_identity_error excluded_small_field "
               "excluded_small_number\n";
    }
    out << std::setprecision(17);
    for (size_t i = 0; i < m.cell_audit.size(); ++i) {
        const CellPowerAudit& c = m.cell_audit[i];
        out << m.name << " " << c.global_cell << " " << c.e_cell << " "
            << c.nxhalf << " " << c.delta_ke_u << " " << c.j_x << " "
            << c.u_work_current << " " << c.x_flux_current << " "
            << c.signed_current_difference << " "
            << c.signed_power_difference << " " << c.v_u_power << " "
            << c.v_x_flux << " " << c.power_current_identity_error << " "
            << (c.excluded_small_field ? 1 : 0) << " "
            << (c.excluded_small_number ? 1 : 0) << "\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    Args args;
    const bool parsed = parse_args(argc, argv, args);
    std::vector<CaseSpec> specs;
    if (parsed && args.test_case == "all" && mpi_size == 1) {
        specs.push_back({"uniform_drift_negative_field", false, false, true, true, false, -1.0e10});
        specs.push_back({"uniform_drift_positive_field", false, false, true, true, false, 1.0e10});
        specs.push_back({"single_velocity_negative_field", false, true, true, true, false, -1.0e10});
        specs.push_back({"single_velocity_positive_field", false, true, true, true, false, 1.0e10});
        specs.push_back({"single_velocity_gradient_negative_field", true, true, true, true, false, -1.0e10});
        specs.push_back({"single_velocity_gradient_positive_field", true, true, true, true, false, 1.0e10});
        specs.push_back({"gradient_limiter_negative_field", true, false, true, false, true, -1.0e10});
        specs.push_back({"open_endpoint_drift", true, false, true, false, false, -1.0e10});
    } else if (parsed && args.test_case == "uniform" && mpi_size == 1) {
        specs.push_back({"uniform_drift_negative_field", false, false, true, true, false, -1.0e10});
    } else if (parsed && args.test_case == "gradient" && mpi_size == 1) {
        specs.push_back({"gradient_limiter_negative_field", true, false, true, false, true, -1.0e10});
    } else if (parsed && args.test_case == "mpi-shared-face" &&
               (mpi_size == 2 || mpi_size == 5)) {
        specs.push_back({"mpi_shared_face_drift", true, false, true, false, false, -1.0e10});
    }

    bool harness_pass = parsed && !specs.empty();
    std::vector<CaseMetrics> results;
    for (size_t i = 0; i < specs.size(); ++i) {
        CaseMetrics result = run_case(
            specs[i], rank, mpi_size, args.x_transport_velocity_mode);
        harness_pass = harness_pass && result.harness_pass;
        result.physical_pairing_pass = result.harness_pass &&
            (!specs[i].baseline_case ||
             result.relative_error <= result.pairing_tolerance);
        results.push_back(result);
    }

    bool baseline_pass = true;
    double baseline_max_relative = 0.0;
    bool limiter_case_present = false;
    bool limiter_related_observable = false;
    double limiter_relative = 0.0;
    bool open_endpoint_pass = true;
    bool p3v1_power_identity_pass = true;
    int p3v1_single_velocity_cases = 0;
    for (size_t i = 0; i < specs.size(); ++i) {
        const CaseMetrics& result = results[i];
        if (specs[i].baseline_case) {
            baseline_pass = baseline_pass && result.physical_pairing_pass;
            baseline_max_relative = std::max(
                baseline_max_relative, result.relative_error);
        }
        if (specs[i].limiter_case) {
            limiter_case_present = true;
            limiter_relative = result.relative_error;
            limiter_related_observable =
                result.x1_linear_fraction > 0.0 ||
                result.x2_linear_fraction > 0.0;
        }
        if (specs[i].open_endpoint)
            open_endpoint_pass = open_endpoint_pass && result.harness_pass;
        if (specs[i].single_velocity && specs[i].baseline_case) {
            ++p3v1_single_velocity_cases;
            p3v1_power_identity_pass =
                p3v1_power_identity_pass &&
                result.power_current_identity_pass;
        }
    }
    const bool limiter_significant = limiter_case_present &&
        limiter_relative > std::max(kPairingTolerance,
                                    2.0 * baseline_max_relative);
    const bool physical_pass = harness_pass && baseline_pass &&
        limiter_case_present && limiter_significant &&
        limiter_related_observable && open_endpoint_pass;
    const bool full_physical_evaluation =
        mpi_size == 1 && args.test_case == "all";
    // P3-V.2 analytic-cell-center is a regression/control run.  Its known
    // P3.0 physical pairing defect must remain visible in the result, but it
    // must not make the control command fail before the energy-conjugate A/B
    // run can be executed.  The new mode retains the physical pairing gate.
    const bool analytic_mode =
        args.x_transport_velocity_mode ==
        XTransportVelocityMode::ANALYTIC_CELL_CENTER;
    const bool analytic_regression_pass = harness_pass &&
        p3v1_power_identity_pass && open_endpoint_pass;
    const bool command_status_pass = full_physical_evaluation
        ? (analytic_mode ? analytic_regression_pass : physical_pass)
        : harness_pass;

    if (!args.result.empty()) {
        std::string cell_path = args.result + ".cells";
        if (mpi_size > 1)
            cell_path += "_rank" + std::to_string(rank);
        std::ofstream cell_out(cell_path.c_str(), std::ios::trunc);
        if (cell_out) {
            for (size_t i = 0; i < results.size(); ++i)
                write_cell_audit(cell_out, results[i]);
        }
    }
    double uniform_relative = std::numeric_limits<double>::quiet_NaN();
    double gradient_relative = std::numeric_limits<double>::quiet_NaN();
    double uniform_continuity = std::numeric_limits<double>::quiet_NaN();
    double gradient_continuity = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].name == "uniform_drift_negative_field") {
            uniform_relative = results[i].relative_error;
            uniform_continuity = results[i].continuity_residual;
        }
        if (results[i].name == "gradient_limiter_negative_field") {
            gradient_relative = results[i].relative_error;
            gradient_continuity = results[i].continuity_residual;
        }
    }

    if (rank == 0 && !args.result.empty()) {
        std::ofstream out(args.result.c_str(), std::ios::trunc);
        if (out) {
            out << std::setprecision(17)
                << "status=" << (command_status_pass ? "PASS" : "FAIL") << "\n"
                << "uniform_relative=" << uniform_relative << "\n"
                << "gradient_relative=" << gradient_relative << "\n"
                << "uniform_continuity=" << uniform_continuity << "\n"
                << "gradient_continuity=" << gradient_continuity << "\n"
                << "harness_integrity_pass=" << (harness_pass ? 1 : 0) << "\n"
                << "physical_pairing_acceptance="
                << (full_physical_evaluation
                    ? (physical_pass ? "PASS" : "FAIL")
                    : "NOT_EVALUATED_PARTIAL_CASES") << "\n"
                << "physical_pairing_tolerance=" << kPairingTolerance << "\n"
                << "pairing_energy_floor=" << kPairingEnergyFloor << "\n"
                << "analytic_regression_acceptance="
                << (analytic_mode
                    ? (analytic_regression_pass ? "PASS" : "FAIL")
                    : "NOT_APPLICABLE") << "\n"
                << "nxhalf_definition=cell_integrated_number_m^-2\n"
                << "vx_flux_definition=Jx_times_dx_over_qN_xhalf\n"
                << "limiter_activity_definition=ppm_reconstruction_changed_cells\n"
                << "baseline_physical_pairing_pass=" << (baseline_pass ? 1 : 0) << "\n"
                << "baseline_max_relative=" << baseline_max_relative << "\n"
                << "limiter_case_present=" << (limiter_case_present ? 1 : 0) << "\n"
                << "limiter_relative=" << limiter_relative << "\n"
                << "limiter_significant=" << (limiter_significant ? 1 : 0) << "\n"
                << "limiter_related_observable="
                << (limiter_related_observable ? 1 : 0) << "\n"
                << "open_endpoint_harness_pass="
                << (open_endpoint_pass ? 1 : 0) << "\n"
                << "p3v1_single_velocity_case_count="
                << p3v1_single_velocity_cases << "\n"
                << "p3v1_power_identity_pass="
                << (p3v1_power_identity_pass ? 1 : 0) << "\n"
                << "time_layer_sequence=Tx_half_Tu_full_Tx_half\n"
                << "charge_current_source=Qx1_plus_Qx2_over_dt\n"
                << "x_transport_velocity_mode="
                << (args.x_transport_velocity_mode ==
                        XTransportVelocityMode::ENERGY_CONJUGATE_CELL
                    ? "energy-conjugate" : "analytic-cell-center") << "\n"
                << "upar_drift=0.35\n"
                << "mpi_size=" << mpi_size << "\n"
                << "p3_release_acceptance="
                << ((mpi_size == 1) ? "INCOMPLETE_MPI_CASE_NOT_RUN" :
                    "INCOMPLETE_SINGLE_RANK_CASES_NOT_RUN") << "\n";
            for (size_t i = 0; i < results.size(); ++i) write_case(out, results[i]);
        } else {
            harness_pass = false;
        }
    }

    double pass_value = command_status_pass ? 1.0 : 0.0;
    double global_pass = 0.0;
    MPI_Allreduce(&pass_value, &global_pass, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << std::setprecision(17)
                  << "harness_integrity_pass=" << (harness_pass ? 1 : 0)
                  << " physical_pairing_acceptance="
                  << (full_physical_evaluation
                      ? (physical_pass ? "PASS" : "FAIL")
                      : "NOT_EVALUATED_PARTIAL_CASES")
                  << " physical_pairing_tolerance=" << kPairingTolerance
                  << " pairing_energy_floor=" << kPairingEnergyFloor
                  << " analytic_regression_acceptance="
                  << (analytic_mode
                      ? (analytic_regression_pass ? "PASS" : "FAIL")
                      : "NOT_APPLICABLE")
                  << " nxhalf_definition=cell_integrated_number_m^-2"
                  << " vx_flux_definition=Jx_times_dx_over_qN_xhalf"
                  << " limiter_activity_definition=ppm_reconstruction_changed_cells"
                  << " baseline_physical_pairing_pass="
                  << (baseline_pass ? 1 : 0)
                  << " baseline_max_relative=" << baseline_max_relative
                  << " limiter_significant=" << (limiter_significant ? 1 : 0)
                  << " limiter_related_observable="
                  << (limiter_related_observable ? 1 : 0)
                  << " p3v1_power_identity_pass="
                  << (p3v1_power_identity_pass ? 1 : 0)
                  << " time_layer_sequence=Tx_half_Tu_full_Tx_half"
                  << " charge_current_source=Qx1_plus_Qx2_over_dt\n"
                  << " x_transport_velocity_mode="
                  << (args.x_transport_velocity_mode ==
                          XTransportVelocityMode::ENERGY_CONJUGATE_CELL
                      ? "energy-conjugate" : "analytic-cell-center")
                  << "status=" << (global_pass > 0.5 ? "PASS" : "FAIL")
                  << "\n";
    }
    MPI_Finalize();
    return global_pass > 0.5 ? 0 : 1;
}
