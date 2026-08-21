#include "beam_pic.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "parameters.h"
#include "species.h"
#include "vpfp_diagnostics.h"
#include "vpfp_integrator.h"
#include "vpfp_checkpoint.h"
#include "vpfp_time_control.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

namespace {
const char* bulk_collision_integrator_name(BulkCollisionIntegrator value)
{
    return value == BulkCollisionIntegrator::CHANG_COOPER_FLUX
        ? "chang-cooper-flux" : "bgk-validation";
}

const char* x_transport_velocity_mode_name(XTransportVelocityMode value)
{
    return value == XTransportVelocityMode::ENERGY_CONJUGATE_CELL
        ? "energy-conjugate" : "analytic-cell-center";
}

const char* background_phase_space_mode_name(BackgroundPhaseSpaceMode value)
{
    return value == BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY
        ? "joint-midpoint-energy" : "strang-ppm";
}

struct Options {
    double length;
    bool domain_length_explicitly_set;
    double dt_scale;
    double stop_time;
    int diagnostic_level;
    int diagnostic_interval;
    bool tail_stage_trace;
    std::string output_dir;
    ElectrostaticBoundary field_boundary;
    OpenBackgroundBoundaryConfig background_boundary;
    std::string collision_model;
    std::string collision_interface_mode;
    BulkCollisionIntegrator bulk_collision_integrator;
    std::string collision_coefficient_file;
    std::string checkpoint_dir;
    std::string restart_dir;
    int checkpoint_every;
    std::vector<double> checkpoint_times;
    std::vector<double> snapshot_times;
    int stop_after_steps;
    CollisionIntegratorType collision_integrator;
    bool beam_enabled;
    // Stage-H3 hybrid tail options (section 14.8).  Collision/population
    // options are parsed now and must be valid, but the corresponding
    // back-ends activate in H5/H8; unknown values fail immediately.
    std::string background_tail_mode;
    double tail_convert_energy_mev;
    int tail_conversion_upar_bins;
    int tail_conversion_energy_bins;
    std::uint64_t tail_max_particle_count;
    double tail_max_number_fraction;
    int tail_target_particles_per_bin;
    int tail_max_particles_per_bin;
    int tail_population_control_interval;
    double tail_max_weight_ratio;
    std::string tail_return_mode;
    double tail_return_energy_mev;
    int tail_return_residence_steps;
    int tail_return_max_stencil_radius;
    double tail_return_moment_tolerance;
    bool restart_allow_return_config_change;
    bool restart_allow_dt_scale_change;
    double restart_source_dt_scale;
    // JC4 (section 7.1): field-particle coupling CLI configuration.  The
    // production default is legacy; discrete-gradient is only activated by
    // explicit opt-in.
    std::string field_particle_coupling;
    XTransportVelocityMode x_transport_velocity_mode;
    BackgroundPhaseSpaceMode background_phase_space_mode;
    int field_particle_max_iters;
    double field_particle_relaxation;
    double field_particle_field_tol;
    double field_particle_pairing_tol;
    bool restart_allow_field_particle_coupling_change;
    std::string tail_collision_kernel;
    std::string tail_collision_weight_mode;
    int tail_collision_max_substeps;
    int tail_collision_max_particle_growth;
    bool tail_cell_moment_audit;
    int tail_cell_moment_audit_top_cells;
    TailConversionMode tail_conversion_mode;
    int tail_flux_quadrature_order;
    int tail_flux_max_supports;
    int tail_flux_max_created_particles_per_step;
    Options()
        : length(Param::Lx), domain_length_explicitly_set(false),
          dt_scale(1.0), stop_time(Param::t_end),
          diagnostic_level(1), diagnostic_interval(500),
          tail_stage_trace(false),
          output_dir("output_vpfp"),
          collision_model("none"), collision_interface_mode("none"),
          bulk_collision_integrator(BulkCollisionIntegrator::BGK_VALIDATION),
          checkpoint_every(0),
          stop_after_steps(0),
          collision_integrator(CollisionIntegratorType::BACKWARD_EULER),
          beam_enabled(true), background_tail_mode("off"),
          tail_convert_energy_mev(6.0), tail_conversion_upar_bins(4),
          tail_conversion_energy_bins(4), tail_max_particle_count(0),
          tail_max_number_fraction(0.0), tail_target_particles_per_bin(64),
          tail_max_particles_per_bin(1024),
          tail_population_control_interval(0), tail_max_weight_ratio(8.0),
          tail_return_mode("none"),
          tail_return_energy_mev(5.5), tail_return_residence_steps(8),
          tail_return_max_stencil_radius(3),
          tail_return_moment_tolerance(1.0e-12),
          restart_allow_return_config_change(false),
          restart_allow_dt_scale_change(false),
          restart_source_dt_scale(std::numeric_limits<double>::quiet_NaN()),
          field_particle_coupling("legacy"),
          x_transport_velocity_mode(XTransportVelocityMode::ANALYTIC_CELL_CENTER),
          background_phase_space_mode(BackgroundPhaseSpaceMode::STRANG_PPM),
          field_particle_max_iters(12),
          field_particle_relaxation(0.5), field_particle_field_tol(1.0e-8),
          field_particle_pairing_tol(1.0e-8),
          restart_allow_field_particle_coupling_change(false),
          tail_collision_kernel("none"), tail_collision_weight_mode(
              "equal-strata"), tail_collision_max_substeps(1024),
          tail_collision_max_particle_growth(0),
          tail_cell_moment_audit(false), tail_cell_moment_audit_top_cells(64),
          // Section 7.11.17.6: the validated flux-interface path is the
          // production default.  The static baseline remains available only
          // when explicitly requested for the 17B/17E A/B comparisons.
          tail_conversion_mode(TailConversionMode::FLUX_INTERFACE),
          tail_flux_quadrature_order(4), tail_flux_max_supports(7),
          tail_flux_max_created_particles_per_step(0)
    {
        field_boundary = { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
        background_boundary.left_type = BackgroundXBoundaryType::RESERVOIR;
        background_boundary.right_type = BackgroundXBoundaryType::RESERVOIR;
        background_boundary.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
        background_boundary.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    }
};

double parse_double(const char* value) { return std::strtod(value, NULL); }

bool require_value(int& i, int argc, char**, const char* option)
{
    if (i + 1 < argc) { ++i; return true; }
    std::cerr << "missing value for " << option << "\n";
    return false;
}

bool parse_time_list(const char* value, std::vector<double>& times_fs)
{
    std::istringstream stream(value);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) return false;
        char* end = NULL;
        const double t = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || !(t >= 0.0)) return false;
        times_fs.push_back(t);
    }
    std::sort(times_fs.begin(), times_fs.end());
    return !times_fs.empty();
}

bool parse_options(int argc, char** argv, Options& o)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--output-dir") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.output_dir = argv[i]; }
        else if (arg == "--domain-length-um") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.length = parse_double(argv[i]) * Const::micro;
            o.domain_length_explicitly_set = true;
        }
        else if (arg == "--dt-scale") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.dt_scale = parse_double(argv[i]); }
        else if (arg == "--stop-time-fs") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.stop_time = parse_double(argv[i]) * Const::femto; }
        else if (arg == "--stop-after-steps") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.stop_after_steps = std::atoi(argv[i]); }
        else if (arg == "--diagnostic-level") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.diagnostic_level = std::atoi(argv[i]); }
        else if (arg == "--diagnostic-interval") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.diagnostic_interval = std::atoi(argv[i]); }
        else if (arg == "--tail-stage-trace") { o.tail_stage_trace = true; }
        else if (arg == "--checkpoint-dir") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.checkpoint_dir = argv[i]; }
        else if (arg == "--restart-dir") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.restart_dir = argv[i]; }
        else if (arg == "--checkpoint-every") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.checkpoint_every = std::atoi(argv[i]); }
        else if (arg == "--checkpoint-times") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            if (!parse_time_list(argv[i], o.checkpoint_times)) {
                std::cerr << "invalid --checkpoint-times list: " << argv[i] << "\n";
                return false;
            }
        }
        else if (arg == "--snapshot-times") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            if (!parse_time_list(argv[i], o.snapshot_times)) {
                std::cerr << "invalid --snapshot-times list: " << argv[i] << "\n";
                return false;
            }
        }
        else if (arg == "--field-boundary") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value == "dirichlet-phi") {
                o.field_boundary.type = ElectrostaticBoundaryType::DIRICHLET_PHI;
            } else if (value == "left-E") {
                o.field_boundary.type = ElectrostaticBoundaryType::LEFT_E;
            } else {
                std::cerr << "unsupported --field-boundary: " << value
                          << " (expected dirichlet-phi or left-E)\n";
                return false;
            }
        }
        else if (arg == "--field-solver") { if (!require_value(i, argc, argv, arg.c_str())) return false; if (std::string(argv[i]) != "open-gauss") return false; }
        else if (arg == "--left-electric-field") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.field_boundary.e_left = parse_double(argv[i]); }
        else if (arg == "--phi-left") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.field_boundary.phi_left = parse_double(argv[i]); }
        else if (arg == "--phi-right") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.field_boundary.phi_right = parse_double(argv[i]); }
        else if (arg == "--background-x-boundary") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "reservoir" && value != "absorbing" &&
                value != "periodic") {
                std::cerr << "unsupported --background-x-boundary: " << value << "\n";
                return false;
            }
            const bool absorbing = value == "absorbing";
            const bool periodic = value == "periodic";
            o.background_boundary.left_type = periodic
                ? BackgroundXBoundaryType::PERIODIC
                : (absorbing ? BackgroundXBoundaryType::ABSORBING
                             : BackgroundXBoundaryType::RESERVOIR);
            o.background_boundary.right_type = o.background_boundary.left_type;
        }
        else if (arg == "--left-reservoir-density") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.left_reservoir.density = parse_double(argv[i]); }
        else if (arg == "--right-reservoir-density") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.right_reservoir.density = parse_double(argv[i]); }
        else if (arg == "--left-reservoir-temperature-ev") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.left_reservoir.temperature = parse_double(argv[i]) * Const::eV; }
        else if (arg == "--right-reservoir-temperature-ev") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.right_reservoir.temperature = parse_double(argv[i]) * Const::eV; }
        else if (arg == "--left-reservoir-drift") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.left_reservoir.drift_vx = parse_double(argv[i]); }
        else if (arg == "--right-reservoir-drift") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.background_boundary.right_reservoir.drift_vx = parse_double(argv[i]); }
        else if (arg == "--collision-model") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.collision_model = argv[i]; }
        else if (arg == "--collision-interface-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.collision_interface_mode = argv[i];
            if (o.collision_interface_mode != "none" &&
                o.collision_interface_mode != "zero-wall-validation" &&
                o.collision_interface_mode != "exporting-absorbing") {
                std::cerr << "--collision-interface-mode expects none, "
                             "zero-wall-validation, or exporting-absorbing\n";
                return false;
            }
        }
        else if (arg == "--bulk-collision-integrator") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value == "bgk-validation") {
                o.bulk_collision_integrator =
                    BulkCollisionIntegrator::BGK_VALIDATION;
            } else if (value == "chang-cooper-flux") {
                o.bulk_collision_integrator =
                    BulkCollisionIntegrator::CHANG_COOPER_FLUX;
            } else {
                std::cerr << "--bulk-collision-integrator expects "
                             "bgk-validation or chang-cooper-flux\n";
                return false;
            }
        }
        else if (arg == "--collision-integrator") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            if (std::string(argv[i]) == "tr-bdf2") {
                std::cerr << "--collision-integrator tr-bdf2 is unsupported: "
                          << "no production TR-BDF2 implementation exists yet\n";
                return false;
            }
            if (std::string(argv[i]) != "backward-euler") {
                std::cerr << "unsupported --collision-integrator: " << argv[i] << "\n";
                return false;
            }
            o.collision_integrator = CollisionIntegratorType::BACKWARD_EULER;
        }
        else if (arg == "--collision-coefficient-file") { if (!require_value(i, argc, argv, arg.c_str())) return false; o.collision_coefficient_file = argv[i]; }
        else if (arg == "--time-integrator") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            if (std::string(argv[i]) != "strang-csl") {
                std::cerr << "unsupported --time-integrator: " << argv[i]
                          << " (only strang-csl may enter production)\n";
                return false;
            }
        }
        else if (arg == "--beam-boundary") { if (!require_value(i, argc, argv, arg.c_str())) return false; if (std::string(argv[i]) != "open") return false; }
        else if (arg == "--beam-enabled") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "0" && value != "1") {
                std::cerr << "--beam-enabled expects 0 or 1\n";
                return false;
            }
            o.beam_enabled = value == "1";
        }
        else if (arg == "--background-tail-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "off" && value != "none" && value != "pic") {
                std::cerr << "--background-tail-mode expects none/off or pic\n";
                return false;
            }
            o.background_tail_mode = value == "none" ? "off" : value;
        }
        else if (arg == "--tail-convert-energy-mev") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_convert_energy_mev = parse_double(argv[i]);
            if (!(o.tail_convert_energy_mev > 0.0)) return false;
        }
        else if (arg == "--tail-conversion-upar-bins") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_conversion_upar_bins = std::atoi(argv[i]);
            if (!(o.tail_conversion_upar_bins > 0)) return false;
        }
        else if (arg == "--tail-conversion-energy-bins") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_conversion_energy_bins = std::atoi(argv[i]);
            if (!(o.tail_conversion_energy_bins > 0)) return false;
        }
        else if (arg == "--tail-cell-moment-audit") {
            o.tail_cell_moment_audit = true;
        }
        else if (arg == "--tail-cell-moment-audit-top-cells") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_cell_moment_audit_top_cells = std::atoi(argv[i]);
            if (!(o.tail_cell_moment_audit_top_cells > 0)) return false;
        }
        else if (arg == "--tail-conversion-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value == "static-cell") o.tail_conversion_mode = TailConversionMode::STATIC_CELL;
            else if (value == "flux-audit") o.tail_conversion_mode = TailConversionMode::FLUX_AUDIT;
            else if (value == "flux-interface") o.tail_conversion_mode = TailConversionMode::FLUX_INTERFACE;
            else {
                std::cerr << "--tail-conversion-mode expects static-cell, flux-audit or flux-interface\n";
                return false;
            }
        }
        else if (arg == "--tail-flux-quadrature-order") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_flux_quadrature_order = std::atoi(argv[i]);
            if (o.tail_flux_quadrature_order != 4 &&
                o.tail_flux_quadrature_order != 8) return false;
        }
        else if (arg == "--tail-flux-max-supports") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_flux_max_supports = std::atoi(argv[i]);
            if (o.tail_flux_max_supports <= 0) return false;
        }
        else if (arg == "--tail-flux-max-created-particles-per-step") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_flux_max_created_particles_per_step = std::atoi(argv[i]);
            if (o.tail_flux_max_created_particles_per_step < 0) return false;
        }
        else if (arg == "--tail-max-particle-count") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const long long value = std::atoll(argv[i]);
            if (value < 0) return false;
            o.tail_max_particle_count =
                static_cast<std::uint64_t>(value);
        }
        else if (arg == "--tail-max-number-fraction") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_max_number_fraction = parse_double(argv[i]);
            if (!(o.tail_max_number_fraction >= 0.0)) return false;
        }
        else if (arg == "--tail-target-particles-per-bin") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_target_particles_per_bin = std::atoi(argv[i]);
            if (!(o.tail_target_particles_per_bin > 0)) return false;
        }
        else if (arg == "--tail-max-particles-per-bin") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_max_particles_per_bin = std::atoi(argv[i]);
            if (!(o.tail_max_particles_per_bin > 0)) return false;
        }
        else if (arg == "--tail-population-control-interval") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_population_control_interval = std::atoi(argv[i]);
            if (!(o.tail_population_control_interval >= 0)) return false;
        }
        else if (arg == "--tail-max-weight-ratio") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_max_weight_ratio = parse_double(argv[i]);
            if (!(o.tail_max_weight_ratio > 1.0)) return false;
        }
        else if (arg == "--tail-return-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "none" && value != "hysteretic") {
                std::cerr << "--tail-return-mode expects none or hysteretic\n";
                return false;
            }
            o.tail_return_mode = value;
        }
        else if (arg == "--tail-return-energy-mev") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_return_energy_mev = parse_double(argv[i]);
            if (!(o.tail_return_energy_mev > 0.0)) return false;
        }
        else if (arg == "--tail-return-residence-steps") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_return_residence_steps = std::atoi(argv[i]);
            if (o.tail_return_residence_steps < 1) return false;
        }
        else if (arg == "--tail-return-max-stencil-radius") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_return_max_stencil_radius = std::atoi(argv[i]);
            if (o.tail_return_max_stencil_radius < 1 ||
                o.tail_return_max_stencil_radius > 3) return false;
        }
        else if (arg == "--tail-return-moment-tolerance") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_return_moment_tolerance = parse_double(argv[i]);
            if (!(o.tail_return_moment_tolerance > 0.0)) return false;
        }
        else if (arg == "--restart-allow-return-config-change") {
            o.restart_allow_return_config_change = true;
        }
        else if (arg == "--field-particle-coupling") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "legacy" && value != "discrete-gradient") {
                std::cerr << "--field-particle-coupling expects legacy or "
                             "discrete-gradient\n";
                return false;
            }
            o.field_particle_coupling = value;
        }
        else if (arg == "--x-transport-velocity-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value == "analytic-cell-center") {
                o.x_transport_velocity_mode =
                    XTransportVelocityMode::ANALYTIC_CELL_CENTER;
            } else if (value == "energy-conjugate") {
                o.x_transport_velocity_mode =
                    XTransportVelocityMode::ENERGY_CONJUGATE_CELL;
            } else {
                std::cerr << "--x-transport-velocity-mode expects "
                             "analytic-cell-center or energy-conjugate\n";
                return false;
            }
        }
        else if (arg == "--background-phase-space-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value == "strang-ppm")
                o.background_phase_space_mode = BackgroundPhaseSpaceMode::STRANG_PPM;
            else if (value == "joint-midpoint-energy")
                o.background_phase_space_mode = BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY;
            else {
                std::cerr << "--background-phase-space-mode expects strang-ppm "
                             "or joint-midpoint-energy\n";
                return false;
            }
        }
        else if (arg == "--field-particle-max-iters") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.field_particle_max_iters = std::atoi(argv[i]);
            if (o.field_particle_max_iters <= 0) {
                std::cerr << "--field-particle-max-iters must be > 0\n";
                return false;
            }
        }
        else if (arg == "--field-particle-relaxation") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.field_particle_relaxation = parse_double(argv[i]);
            if (!(o.field_particle_relaxation > 0.0) ||
                o.field_particle_relaxation > 1.0) {
                std::cerr << "--field-particle-relaxation must be in (0,1]\n";
                return false;
            }
        }
        else if (arg == "--field-particle-field-tol") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.field_particle_field_tol = parse_double(argv[i]);
            if (!std::isfinite(o.field_particle_field_tol) ||
                o.field_particle_field_tol <= 0.0) {
                std::cerr << "--field-particle-field-tol must be a finite "
                             "positive number\n";
                return false;
            }
        }
        else if (arg == "--field-particle-pairing-tol") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.field_particle_pairing_tol = parse_double(argv[i]);
            if (!std::isfinite(o.field_particle_pairing_tol) ||
                o.field_particle_pairing_tol <= 0.0) {
                std::cerr << "--field-particle-pairing-tol must be a finite "
                             "positive number\n";
                return false;
            }
        }
        else if (arg == "--restart-allow-field-particle-coupling-change") {
            o.restart_allow_field_particle_coupling_change = true;
        }
        else if (arg == "--restart-allow-dt-scale-change") {
            o.restart_allow_dt_scale_change = true;
        }
        else if (arg == "--restart-source-dt-scale") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.restart_source_dt_scale = parse_double(argv[i]);
        }
        else if (arg == "--tail-collision-kernel") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "none" && value != "coulomb-nanbu-perez" &&
                value != "kramers-moyal-sde") {
                std::cerr << "--tail-collision-kernel expects none, "
                             "coulomb-nanbu-perez or kramers-moyal-sde\n";
                return false;
            }
            o.tail_collision_kernel = value;
        }
        else if (arg == "--tail-collision-weight-mode") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const std::string value(argv[i]);
            if (value != "equal-strata" && value != "virtual-split") {
                std::cerr << "--tail-collision-weight-mode supports only "
                             "equal-strata or virtual-split (H8)\n";
                return false;
            }
            o.tail_collision_weight_mode = value;
        }
        else if (arg == "--tail-collision-max-substeps") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            o.tail_collision_max_substeps = std::atoi(argv[i]);
            if (!(o.tail_collision_max_substeps > 0)) return false;
        }
        else if (arg == "--tail-collision-max-particle-growth") {
            if (!require_value(i, argc, argv, arg.c_str())) return false;
            const long long value = std::atoll(argv[i]);
            if (value < 0) return false;
            o.tail_collision_max_particle_growth =
                static_cast<int>(value);
        }
        else { std::cerr << "unknown option: " << arg << "\n"; return false; }
    }
    const bool source_dt_scale_set =
        std::isfinite(o.restart_source_dt_scale) &&
        o.restart_source_dt_scale > 0.0;
    if (o.restart_allow_dt_scale_change != source_dt_scale_set) {
        std::cerr << "--restart-allow-dt-scale-change and "
                     "--restart-source-dt-scale must be supplied together\n";
        return false;
    }
    if (o.restart_allow_dt_scale_change && o.restart_dir.empty()) {
        std::cerr << "--restart-allow-dt-scale-change requires --restart-dir\n";
        return false;
    }
    // JC4 (section 7.4): field-particle coupling override legal combinations.
    // The override only permits a legacy checkpoint to start a
    // discrete-gradient A/B audit; it must not be usable without a restart or
    // while the requested mode is itself legacy.
    if (o.restart_allow_field_particle_coupling_change &&
        o.restart_dir.empty()) {
        std::cerr << "--restart-allow-field-particle-coupling-change "
                     "requires --restart-dir\n";
        return false;
    }
    if (o.background_phase_space_mode ==
            BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY &&
        (o.collision_model != "none" || o.background_tail_mode != "off" ||
         o.beam_enabled || o.field_particle_coupling != "legacy" ||
         o.background_boundary.left_type != BackgroundXBoundaryType::PERIODIC ||
         o.background_boundary.right_type != BackgroundXBoundaryType::PERIODIC)) {
        std::cerr << "joint-midpoint-energy requires collision-model none, "
                     "background-tail-mode off, beam-enabled 0, and legacy "
                     "field-particle coupling, and background-x-boundary periodic\n";
        return false;
    }
    if (o.restart_allow_field_particle_coupling_change &&
        o.field_particle_coupling != "discrete-gradient") {
        std::cerr << "--restart-allow-field-particle-coupling-change requires "
                     "--field-particle-coupling discrete-gradient\n";
        return false;
    }
    return o.length > 0.0 && o.dt_scale > 0.0 && o.stop_time >= 0.0 &&
           o.stop_after_steps >= 0 &&
           o.tail_max_particles_per_bin >= o.tail_target_particles_per_bin;
}

bool validate_domain_length(double length, int& nx_open)
{
    const double ratio = length / Param::dx;
    const double cells = std::floor(ratio + 0.5);
    if (!(cells >= 1.0) || cells > static_cast<double>(std::numeric_limits<int>::max())) {
        return false;
    }
    const double reconstructed = cells * Param::dx;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max(length, 1.0);
    if (std::fabs(reconstructed - length) > tolerance) return false;
    nx_open = static_cast<int>(cells);
    return true;
}

std::string boundary_type_name(ElectrostaticBoundaryType type)
{
    return type == ElectrostaticBoundaryType::DIRICHLET_PHI
               ? "dirichlet-phi" : "left-E";
}

std::string boundary_side_name(BackgroundXBoundaryType type)
{
    if (type == BackgroundXBoundaryType::PERIODIC) return "periodic";
    return type == BackgroundXBoundaryType::RESERVOIR ? "reservoir" : "absorbing";
}

void hash_u64(std::uint64_t& hash, std::uint64_t value)
{
    hash ^= value;
    hash *= 0x100000001b3ULL;
}

void hash_string(std::uint64_t& hash, const std::string& value)
{
    for (size_t i = 0; i < value.size(); ++i)
        hash_u64(hash, static_cast<unsigned char>(value[i]));
    hash_u64(hash, 0xffU);
}

void hash_double(std::uint64_t& hash, double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    hash_u64(hash, bits);
}

std::uint64_t physical_config_hash(
    const Options& options, const HybridVelocityPartition& partition,
    bool include_return_numeric_parameters = true,
    bool include_x_transport_velocity = true)
{
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    hash_u64(hash, partition.config_hash);
    hash_double(hash, options.length);
    hash_double(hash, options.dt_scale);
    hash_u64(hash, options.beam_enabled ? 1ULL : 0ULL);
    hash_string(hash, boundary_type_name(options.field_boundary.type));
    hash_double(hash, options.field_boundary.e_left);
    hash_double(hash, options.field_boundary.phi_left);
    hash_double(hash, options.field_boundary.phi_right);
    hash_string(hash, boundary_side_name(options.background_boundary.left_type));
    hash_string(hash, boundary_side_name(options.background_boundary.right_type));
    hash_double(hash, options.background_boundary.left_reservoir.density);
    hash_double(hash, options.background_boundary.left_reservoir.temperature);
    hash_double(hash, options.background_boundary.left_reservoir.drift_vx);
    hash_double(hash, options.background_boundary.right_reservoir.density);
    hash_double(hash, options.background_boundary.right_reservoir.temperature);
    hash_double(hash, options.background_boundary.right_reservoir.drift_vx);
    hash_string(hash, options.collision_model);
    hash_string(hash, options.collision_interface_mode);
    // Preserve the historical Strang physical hash exactly.  The explicit
    // J1 mode receives an additional mode token so its checkpoints cannot be
    // resumed under the default production path.
    if (options.background_phase_space_mode ==
        BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY) {
        hash_string(hash, background_phase_space_mode_name(
            options.background_phase_space_mode));
    }
    if (include_x_transport_velocity) {
        hash_string(hash, x_transport_velocity_mode_name(
            options.x_transport_velocity_mode));
        hash_u64(hash, 1ULL); // x transport velocity table schema v1
    }
    hash_string(hash, bulk_collision_integrator_name(
        options.bulk_collision_integrator));
    hash_string(hash, options.tail_return_mode);
    if (include_return_numeric_parameters) {
        hash_double(hash, options.tail_return_energy_mev);
        hash_u64(hash, static_cast<std::uint64_t>(
            options.tail_return_residence_steps));
        hash_u64(hash, static_cast<std::uint64_t>(
            options.tail_return_max_stencil_radius));
        hash_double(hash, options.tail_return_moment_tolerance);
    }
    hash_string(hash, tail_conversion_mode_name(options.tail_conversion_mode));
    hash_double(hash, options.tail_convert_energy_mev);
    hash_u64(hash, static_cast<std::uint64_t>(options.tail_conversion_upar_bins));
    hash_u64(hash, static_cast<std::uint64_t>(options.tail_conversion_energy_bins));
    hash_u64(hash, static_cast<std::uint64_t>(options.tail_flux_quadrature_order));
    hash_u64(hash, static_cast<std::uint64_t>(options.tail_flux_max_supports));
    hash_u64(hash, static_cast<std::uint64_t>(
        options.tail_flux_max_created_particles_per_step));
    return hash;
}

std::uint64_t diagnostic_config_hash(const Options& options)
{
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    hash_u64(hash, static_cast<std::uint64_t>(options.diagnostic_level));
    hash_u64(hash, static_cast<std::uint64_t>(options.diagnostic_interval));
    hash_u64(hash, options.tail_stage_trace ? 1ULL : 0ULL);
    return hash;
}

std::string checkpoint_time_directory(const std::string& output_dir,
                                      double target_time, double time, int step)
{
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer),
                  "%s/checkpoint_target%.9gfs_t%.9gfs_step%d",
                  output_dir.c_str(), target_time / Const::femto,
                  time / Const::femto, step);
    return std::string(buffer);
}

std::string failure_code_name(int code)
{
    switch (code) {
        case 1: return "invalid-dt-or-uninitialized";
        case 2: return "remap-non-finite";
        case 3: return "gauss-closure-or-non-finite";
        case 5: return "collision-failure";
        case 6: return "background-ledger-or-negative-mass";
        case 7: return "velocity-tail-loss-above-threshold";
        case 8: return "beam-number-balance";
        case 10: return "bulk-tail-conversion";
        case 12: return "collision-face-export-not-enabled";
        case 11: return "tail-resource-gate";
        // JC4 (section 7.4 item 8 / section 6.5): field-particle coupling
        // failure codes 201--208 must have explicit names, never "unknown".
        case 201: return "field_particle_trial_nonfinite";
        case 202: return "field_particle_poisson_failure";
        case 203: return "field_particle_pairing_field_failure";
        case 204: return "field_particle_work_identity_failure";
        case 205: return "field_particle_not_converged";
        case 206: return "post_field_charge_invariance_failure";
        case 207: return "field_particle_mpi_consensus_failure";
        case 208: return "accepted_field_poisson_mismatch";
        default: return "unknown";
    }
}

void write_startup_phase(const std::string& output_dir, int rank,
                         const char* phase)
{
    // Startup-only collective breadcrumbs locate hangs before the first
    // accepted step without changing the production time-stepping path.
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream out((output_dir + "/vpfp_startup_progress.dat").c_str(),
                          std::ios::app);
        if (out) {
            out << std::setprecision(17) << MPI_Wtime() << " " << phase << "\n";
            out.flush();
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Options options;
    const bool parsed = parse_options(argc, argv, options);
    int parsed_all = parsed ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &parsed_all, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if (!parsed_all) { MPI_Finalize(); return 2; }

    int nx_open = Param::nx;
    bool valid_domain = true;
    if (options.domain_length_explicitly_set) {
        valid_domain = validate_domain_length(options.length, nx_open);
    }
    int valid_domain_all = valid_domain ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &valid_domain_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!valid_domain_all) {
        if (rank == 0) {
            std::cerr << "--domain-length-um must be an integer multiple of the "
                      << "fixed production dx (" << Param::dx << " m)\n";
        }
        MPI_Finalize();
        return 2;
    }

    // Establish the requested output root before any large workspace
    // allocation or restart I/O.  A missing nested directory must never make
    // a production launch appear to run silently.
    VpfpDiagnostics diagnostics;
    if (!diagnostics.init(options.output_dir, rank, options.diagnostic_level,
                          options.diagnostic_interval)) {
        if (rank == 0) {
            std::cerr << "VPFP startup failed: cannot create or write output directory "
                      << options.output_dir << "\n";
        }
        MPI_Finalize();
        return 2;
    }
    write_startup_phase(options.output_dir, rank, "diagnostics_root_ready");

    SpatialGrid grid;
    grid.init_with_domain(rank, size, nx_open, options.length);
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    write_startup_phase(options.output_dir, rank, "grid_and_initial_state_ready");

    OpenBackgroundBoundary background_boundary(options.background_boundary);
    background_boundary.fill_ghosts(electrons, grid, rank, size);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid, options.field_boundary);
    write_startup_phase(options.output_dir, rank, "boundaries_and_field_ready");
    std::unique_ptr<CollisionCoefficientProvider> provider;
    if (options.collision_model == "none") provider.reset(new ZeroCollisionCoefficients());
    else if (options.collision_model == "moment-closure") {
        // Stage H7 moment-closure bulk collision (section 10.2.1 mode 2):
        // velocity-dependent coefficients from the local cell moments with
        // a fixed Coulomb logarithm (documented closure approximation).
        // Section 19.3: the closure temperature is made self-consistent
        // with the discrete grid moments inside the collision scan.
        provider.reset(new MomentClosureCollisionCoefficients(20.0));
    }
    else if (options.collision_model == "prescribed" && !options.collision_coefficient_file.empty()) {
        CylindricalCollisionCoefficients c = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        std::ifstream input(options.collision_coefficient_file.c_str());
        if (!(input >> c.a_parallel >> c.a_perp >> c.d_parallel_parallel >>
              c.d_parallel_perp >> c.d_perp_perp)) {
            if (rank == 0) std::cerr << "cannot read prescribed collision coefficients\n";
            MPI_Finalize(); return 2;
        }
        provider.reset(new PrescribedCollisionCoefficients(c));
    } else { if (rank == 0) std::cerr << "unsupported collision model or missing coefficient file\n"; MPI_Finalize(); return 2; }
    CylindricalFokkerPlanckCollision collision(*provider, options.collision_integrator);
    collision.set_bulk_integrator(options.bulk_collision_integrator);
    HybridVelocityPartition tail_partition;
    BulkTailConverter tail_converter;
    const bool tail_enabled = options.background_tail_mode == "pic";
    if (!tail_enabled && options.tail_return_mode != "none") {
        if (rank == 0) {
            std::cerr << "--tail-return-mode hysteretic requires --background-tail-mode pic\n";
        }
        MPI_Finalize();
        return 2;
    }
    // With tail disabled the conversion mode is a no-op.  Keeping the
    // selected mode in the manifest makes a tail-off production run explicit
    // without forcing an unrelated representation on the Eulerian-only path.
    if (tail_enabled) {
        if (options.tail_return_mode == "hysteretic" &&
            !(options.tail_return_energy_mev > 0.0 &&
              options.tail_return_energy_mev <
                  options.tail_convert_energy_mev)) {
            if (rank == 0) {
                std::cerr << "hysteretic tail return requires 0 < K_in < K_out\n";
            }
            MPI_Finalize();
            return 2;
        }
        try {
            tail_partition.init(electrons.cgrid,
                                options.tail_convert_energy_mev, 1.0,
                                options.tail_conversion_upar_bins,
                                options.tail_conversion_energy_bins);
        } catch (const std::runtime_error& error) {
            if (rank == 0) {
                std::cerr << "tail partition init failed: "
                          << error.what() << "\n";
            }
            MPI_Finalize();
            return 2;
        }
    if (tail_enabled &&
        options.tail_conversion_mode != TailConversionMode::STATIC_CELL) {
            std::string topology_error;
            if (!tail_partition.flux_interface_topology_valid(&topology_error)) {
                if (rank == 0) {
                    std::cerr << "tail flux interface topology invalid: "
                              << topology_error << "\n";
                }
                MPI_Finalize();
                return 2;
            }
            if (options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE &&
                options.collision_model != "none" &&
                options.collision_interface_mode != "zero-wall-validation" &&
                (options.bulk_collision_integrator !=
                     BulkCollisionIntegrator::CHANG_COOPER_FLUX ||
                 options.collision_interface_mode != "exporting-absorbing")) {
                if (rank == 0) {
                    std::cerr << "--tail-conversion-mode flux-interface requires "
                                 "--bulk-collision-integrator chang-cooper-flux "
                                 "and --collision-interface-mode "
                                 "exporting-absorbing for real collisions\n";
                }
                MPI_Finalize();
                return 2;
            }
        }
    }
    write_startup_phase(options.output_dir, rank, "tail_partition_ready");
    if (options.tail_cell_moment_audit &&
        tail_converter.subcell_loading_enabled()) {
        if (rank == 0) {
            std::cerr << "--tail-cell-moment-audit requires subcell loading disabled\n";
        }
        MPI_Finalize();
        return 2;
    }
    tail_converter.set_moment_audit_enabled(options.tail_cell_moment_audit,
                                            static_cast<size_t>(options.tail_cell_moment_audit_top_cells));
    VpfpIntegrator integrator(background_boundary, field_solver, collision,
                              tail_partition, tail_converter, tail_enabled);
    integrator.set_beam_enabled(options.beam_enabled);
    // Section 17.15.7.2: level 2 alone opts into the read-only stage-energy
    // audit.  Levels 0/1 retain the previous scan and collective profile.
    integrator.set_stage_energy_audit_enabled(options.diagnostic_level >= 2);
    // Gate I (section 4.6.1): the read-only field-particle power pairing
    // audit is also enabled only by diagnostic level 2.
    integrator.set_field_particle_power_audit_enabled(
        options.diagnostic_level >= 2);
    integrator.set_collision_interface_zero_wall_validation(
        options.collision_interface_mode == "zero-wall-validation");
    integrator.set_collision_interface_exporting_absorbing(
        options.collision_interface_mode == "exporting-absorbing");
    integrator.set_tail_conversion_mode(
        options.tail_conversion_mode, options.tail_flux_quadrature_order,
        static_cast<size_t>(options.tail_flux_max_supports),
        options.tail_flux_max_created_particles_per_step);
    TailBulkReturnConfig return_config;
    return_config.enabled = tail_enabled &&
        options.tail_return_mode == "hysteretic";
    return_config.return_energy_mev = options.tail_return_energy_mev;
    return_config.residence_steps = static_cast<std::uint32_t>(
        options.tail_return_residence_steps);
    return_config.max_stencil_radius = options.tail_return_max_stencil_radius;
    return_config.moment_tolerance = options.tail_return_moment_tolerance;
    integrator.set_tail_bulk_return(return_config);
    integrator.set_tail_limits(options.tail_max_particle_count,
                               options.tail_max_number_fraction);
    if (options.tail_collision_kernel == "coulomb-nanbu-perez") {
        integrator.set_tail_collision(
            TailCollisionKernel::CoulombLandauNanbuPerez, 20.0,
            options.tail_collision_weight_mode == "virtual-split"
                ? TailCollisionWeightMode::VirtualSplit
                : TailCollisionWeightMode::EqualStrata,
            options.tail_collision_max_substeps,
            static_cast<double>(options.tail_collision_max_particle_growth));
    } else if (options.tail_collision_kernel == "kramers-moyal-sde") {
        integrator.set_tail_collision(
            TailCollisionKernel::KramersMoyalSDE, 20.0,
            options.tail_collision_weight_mode == "virtual-split"
                ? TailCollisionWeightMode::VirtualSplit
                : TailCollisionWeightMode::EqualStrata,
            options.tail_collision_max_substeps,
            static_cast<double>(options.tail_collision_max_particle_growth));
    } else {
        integrator.set_tail_collision(
            TailCollisionKernel::None, 20.0,
            options.tail_collision_weight_mode == "virtual-split"
                ? TailCollisionWeightMode::VirtualSplit
                : TailCollisionWeightMode::EqualStrata,
            options.tail_collision_max_substeps,
            static_cast<double>(options.tail_collision_max_particle_growth));
    }
    TailPopulationController::Config population_config;
    population_config.enabled =
        tail_enabled && options.tail_population_control_interval > 0;
    population_config.control_interval =
        options.tail_population_control_interval;
    population_config.target_particles_per_phase_bin =
        options.tail_target_particles_per_bin;
    population_config.max_particles_per_phase_bin =
        options.tail_max_particles_per_bin;
    population_config.max_weight_ratio = options.tail_max_weight_ratio;
    integrator.set_tail_population_control(population_config);
    integrator.set_tail_stage_trace(options.tail_stage_trace,
                                    options.output_dir);
    // JC4 (section 7.1): parse the field-particle coupling configuration and
    // apply it to the integrator before init.  Production default is legacy.
    FieldParticleCouplingConfig coupling_config;
    coupling_config.mode =
        options.field_particle_coupling == "discrete-gradient"
            ? FieldParticleCouplingMode::DiscreteGradient
            : FieldParticleCouplingMode::Legacy;
    coupling_config.max_iterations = options.field_particle_max_iters;
    coupling_config.initial_relaxation = options.field_particle_relaxation;
    coupling_config.field_relative_tolerance = options.field_particle_field_tol;
    coupling_config.pairing_relative_tolerance =
        options.field_particle_pairing_tol;
    integrator.set_field_particle_coupling(coupling_config);
    integrator.set_x_transport_velocity_mode(
        options.x_transport_velocity_mode);
    integrator.set_background_phase_space_mode(
        options.background_phase_space_mode);
    integrator.init(grid);
    write_startup_phase(options.output_dir, rank, "integrator_workspace_ready");
    VpfpRunManifestConfig manifest_config;
    const std::uint64_t physical_hash =
        physical_config_hash(options, tail_partition);
    const std::uint64_t diagnostic_hash = diagnostic_config_hash(options);
    manifest_config.field_boundary =
        boundary_type_name(options.field_boundary.type);
    manifest_config.left_electric_field = options.field_boundary.e_left;
    manifest_config.phi_left = options.field_boundary.phi_left;
    manifest_config.phi_right = options.field_boundary.phi_right;
    manifest_config.background_boundary =
        boundary_side_name(options.background_boundary.left_type) + "/" +
        boundary_side_name(options.background_boundary.right_type);
    manifest_config.collision_model = provider->name();
    manifest_config.collision_interface_mode = options.collision_interface_mode;
    manifest_config.bulk_collision_integrator =
        bulk_collision_integrator_name(options.bulk_collision_integrator);
    manifest_config.collision_induced_conversion =
        options.collision_model != "none" &&
        options.collision_interface_mode != "zero-wall-validation";
    const HybridCollisionConfig& collision_config =
        integrator.tail_collision_config();
    manifest_config.requested_tail_collision_kernel =
        tail_collision_kernel_name(collision_config.requested_kernel);
    manifest_config.tail_tail_collision_backend =
        tail_collision_kernel_name(collision_config.tail_tail_kernel);
    manifest_config.tail_bulk_collision_backend =
        tail_collision_kernel_name(collision_config.tail_bulk_kernel);
    manifest_config.tail_collision_weight_mode =
        tail_collision_weight_mode_name(collision_config.weight_mode);
    manifest_config.tail_collision_max_substeps =
        collision_config.max_substeps;
    manifest_config.tail_collision_max_particle_growth =
        options.tail_collision_max_particle_growth;
    manifest_config.tail_cell_moment_audit = options.tail_cell_moment_audit;
    manifest_config.tail_cell_moment_audit_top_cells =
        options.tail_cell_moment_audit_top_cells;
    manifest_config.tail_subcell_loading = tail_converter.subcell_loading_enabled();
    manifest_config.tail_conversion_mode =
        tail_conversion_mode_name(options.tail_conversion_mode);
    manifest_config.tail_flux_quadrature_order = options.tail_flux_quadrature_order;
    manifest_config.tail_flux_max_supports = options.tail_flux_max_supports;
    manifest_config.tail_flux_max_created_particles_per_step =
        options.tail_flux_max_created_particles_per_step;
    manifest_config.physical_config_hash = physical_hash;
    manifest_config.diagnostic_config_hash = diagnostic_hash;
    manifest_config.interface_topology_hash = tail_partition.topology_mask_hash();
    manifest_config.interface_topology_version =
        tail_partition.interface_topology_version();
    manifest_config.interface_mask_hash =
        tail_partition.tail_owned_mask_hash();
    manifest_config.interface_face_list_hash =
        tail_partition.interface_face_list_hash();
    manifest_config.pair_bulk_bulk = collision_config.pairs.bulk_bulk;
    manifest_config.pair_bulk_tail = collision_config.pairs.bulk_tail;
    manifest_config.pair_tail_bulk = collision_config.pairs.tail_bulk;
    manifest_config.pair_tail_tail = collision_config.pairs.tail_tail;
    manifest_config.population_control_enabled = population_config.enabled;
    manifest_config.population_control_interval =
        population_config.control_interval;
    // JC4 (section 7.4 item 9 / section 7.2): persist coupling config in the
    // run manifest so restart can compare stored vs. requested values.
    manifest_config.coupling_mode = options.field_particle_coupling;
    manifest_config.coupling_max_iters = options.field_particle_max_iters;
    manifest_config.coupling_relaxation = options.field_particle_relaxation;
    manifest_config.coupling_field_tol = options.field_particle_field_tol;
    manifest_config.coupling_pairing_tol = options.field_particle_pairing_tol;
    manifest_config.x_transport_velocity_mode =
        x_transport_velocity_mode_name(options.x_transport_velocity_mode);
    manifest_config.x_transport_velocity_table_schema = 1;
    manifest_config.background_phase_space_mode =
        background_phase_space_mode_name(options.background_phase_space_mode);
    diagnostics.set_run_manifest_config(manifest_config);

    std::vector<double> ion_density(static_cast<size_t>(grid.nx_local), Param::dens);
    if (options.background_phase_space_mode ==
        BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY) {
        // J1 is a manufactured background-only closure test.  Use the actual
        // discrete Maxwellian density so the initial Poisson state is neutral
        // to the same quadrature as the joint operator.
        ion_density = electrons.number_density;
    }
    const double dt = options.dt_scale * Param::dt_multiplier / Param::omega_pe;
    double time = 0.0;
    int step = 0;
    BulkTailConversionDiagnostics initial_tail_conversion;
    write_startup_phase(options.output_dir, rank,
                        options.restart_dir.empty() ? "fresh_initialization_begin"
                                                    : "restart_read_begin");
    if (options.restart_dir.empty() && tail_enabled &&
        options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE) {
        const bool local_initial_ok = integrator.initialize_tail_from_bulk(
            electrons, rank, size, initial_tail_conversion);
        int initial_ok = local_initial_ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &initial_ok, 1, MPI_INT, MPI_LAND,
                      MPI_COMM_WORLD);
        if (!initial_ok) {
            if (rank == 0) {
                std::cerr << "initial flux-interface tail conversion failed"
                          << "; no state was advanced\n";
            }
            MPI_Finalize();
            return 3;
        }
    }
    if (!options.restart_dir.empty()) {
        VpfpCheckpointControl restored = {};
        VpfpCheckpointTailState restored_tail;
        std::string restart_error;
        if (!read_vpfp_checkpoint(options.restart_dir, restored, electrons, beam,
                                  fields, grid, tail_enabled ? &restored_tail : NULL,
                                  rank, size, restart_error)) {
            if (rank == 0) std::cerr << "VPFP restart failed: " << restart_error << "\n";
            MPI_Finalize(); return 2;
        }
        std::string stored_x_transport_velocity_mode;
        int stored_x_transport_velocity_table_schema = 1;
        // JC4 (section 7.5): read and validate the field-particle coupling
        // configuration from the checkpoint manifest.  Old checkpoints missing
        // these keys are treated as legacy defaults; new checkpoints must match
        // the requested coupling mode and numeric parameters.
        {
            std::string stored_mode;
            int stored_max_iters = 12;
            double stored_relaxation = 0.5;
            double stored_field_tol = 1.0e-8;
            double stored_pairing_tol = 1.0e-8;
            read_coupling_config_from_manifest(
                options.restart_dir, stored_mode, stored_max_iters,
                stored_relaxation, stored_field_tol, stored_pairing_tol,
                stored_x_transport_velocity_mode,
                stored_x_transport_velocity_table_schema);
            // §7.5 item 4: legacy->discrete-gradient requires explicit
            // --restart-allow-field-particle-coupling-change (already enforced
            // in parse_options).  Here we check the stored vs. requested mode
            // and numeric parameters.
            if (stored_mode != options.field_particle_coupling) {
                if (!options.restart_allow_field_particle_coupling_change) {
                    if (rank == 0) {
                        std::cerr << "VPFP restart failed: coupling mode mismatch"
                                  << " stored=" << stored_mode
                                  << " requested=" << options.field_particle_coupling
                                  << "; use --restart-allow-field-particle-coupling-change"
                                     " to override\n";
                    }
                    MPI_Finalize();
                    return 2;
                }
                // §7.5 item 5: override only allows coupling config differences;
                // the existing dt/grid/boundary/collision/Tail/Beam identity
                // checks still execute (handled by physical_config_hash below).
                if (rank == 0) {
                    std::cerr << "VPFP restart: accepted coupling mode switch"
                              << " stored=" << stored_mode
                              << " requested=" << options.field_particle_coupling
                              << "\n";
                }
            }
            std::string stored_phase_space_mode;
            read_background_phase_space_mode_from_manifest(
                options.restart_dir, stored_phase_space_mode);
            if (stored_phase_space_mode != background_phase_space_mode_name(
                    options.background_phase_space_mode)) {
                if (rank == 0)
                    std::cerr << "VPFP restart failed: background phase-space "
                                 "mode mismatch stored=" << stored_phase_space_mode
                              << " requested=" << background_phase_space_mode_name(
                                     options.background_phase_space_mode) << "\n";
                MPI_Finalize();
                return 2;
            }
            const std::string requested_x_mode =
                x_transport_velocity_mode_name(options.x_transport_velocity_mode);
            const bool legacy_x_velocity_metadata =
                stored_x_transport_velocity_table_schema == 0;
            const bool x_velocity_mismatch = legacy_x_velocity_metadata
                ? options.x_transport_velocity_mode !=
                    XTransportVelocityMode::ANALYTIC_CELL_CENTER
                : (stored_x_transport_velocity_mode != requested_x_mode ||
                   stored_x_transport_velocity_table_schema != 1);
            if (x_velocity_mismatch) {
                if (rank == 0) {
                    std::cerr << "VPFP restart rejected: x transport velocity "
                              << "mode/schema mismatch stored="
                              << stored_x_transport_velocity_mode << "/"
                              << stored_x_transport_velocity_table_schema
                              << " requested=" << requested_x_mode << "/1\n";
                }
                MPI_Finalize(); return 2;
            }
            // §7.5 item 6: numeric parameter mismatch is also an error unless
            // the override flag is set (same semantics as mode mismatch).
            if (stored_max_iters != options.field_particle_max_iters ||
                std::fabs(stored_relaxation - options.field_particle_relaxation) >
                    1.0e-12 ||
                std::fabs(stored_field_tol - options.field_particle_field_tol) >
                    1.0e-12 ||
                std::fabs(stored_pairing_tol - options.field_particle_pairing_tol) >
                    1.0e-12) {
                if (!options.restart_allow_field_particle_coupling_change) {
                    if (rank == 0) {
                        std::cerr << "VPFP restart failed: coupling numeric"
                                     " parameter mismatch\n";
                    }
                    MPI_Finalize();
                    return 2;
                }
            }
        }
        if (tail_enabled) {
            // Section 12.2/12.3: the checkpoint must carry the tail state
            // and match the current partition/controller configuration;
            // no silent upgrade or physical modification is allowed.
            const bool old_conversion_checkpoint =
                !restored_tail.config.conversion_metadata_present;
            const bool hash_mismatch =
                restored_tail.config.partition_config_hash !=
                tail_partition.config_hash;
            if (!restored_tail.present ||
                (hash_mismatch &&
                 !(old_conversion_checkpoint &&
                   options.tail_conversion_mode != TailConversionMode::FLUX_INTERFACE))) {
                if (rank == 0) {
                    std::cerr << "VPFP restart failed: checkpoint tail state "
                                 "or partition config mismatch\n";
                }
                MPI_Finalize();
                return 2;
            }
            const unsigned long long local_tail_particle_count =
                static_cast<unsigned long long>(restored_tail.tail.particles.size());
            unsigned long long global_tail_particle_count = 0;
            MPI_Allreduce(&local_tail_particle_count, &global_tail_particle_count,
                          1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            const bool preconversion_mode_switch =
                vpfp_preconversion_static_to_flux_restart_allowed(
                    restored_tail.config,
                    tail_conversion_mode_name(options.tail_conversion_mode),
                    static_cast<std::uint64_t>(global_tail_particle_count));
            if (!old_conversion_checkpoint) {
                VpfpCheckpointTailConfig expected_tail_config;
                expected_tail_config.physical_config_hash = physical_hash;
                expected_tail_config.partition_config_hash =
                    tail_partition.config_hash;
                expected_tail_config.conversion_energy_edges =
                    tail_partition.conversion_energy_edges;
                expected_tail_config.conversion_energy_edges_hash =
                    tail_partition.conversion_energy_edges_hash();
                expected_tail_config.convert_energy_mev =
                    options.tail_convert_energy_mev;
                expected_tail_config.return_mode = options.tail_return_mode;
                expected_tail_config.return_energy_mev =
                    options.tail_return_energy_mev;
                expected_tail_config.return_residence_steps =
                    options.tail_return_residence_steps;
                expected_tail_config.return_max_stencil_radius =
                    options.tail_return_max_stencil_radius;
                expected_tail_config.return_moment_tolerance =
                    options.tail_return_moment_tolerance;
                expected_tail_config.upar_bins =
                    options.tail_conversion_upar_bins;
                expected_tail_config.energy_bins =
                    options.tail_conversion_energy_bins;
                expected_tail_config.conversion_mode =
                    tail_conversion_mode_name(options.tail_conversion_mode);
                expected_tail_config.flux_quadrature_order =
                    options.tail_flux_quadrature_order;
                expected_tail_config.flux_max_supports =
                    options.tail_flux_max_supports;
                expected_tail_config.flux_max_created_particles_per_step =
                    options.tail_flux_max_created_particles_per_step;
                expected_tail_config.collision_interface_mode =
                    options.collision_interface_mode;
                expected_tail_config.bulk_collision_integrator =
                    bulk_collision_integrator_name(
                        options.bulk_collision_integrator);
                expected_tail_config.collision_induced_conversion =
                    options.collision_model != "none" &&
                    options.collision_interface_mode !=
                        "zero-wall-validation";
                expected_tail_config.interface_topology_metadata_present =
                    options.tail_conversion_mode ==
                    TailConversionMode::FLUX_INTERFACE;
                expected_tail_config.interface_topology_hash =
                    tail_partition.topology_mask_hash();
                expected_tail_config.interface_topology_detail_metadata_present =
                    options.tail_conversion_mode ==
                    TailConversionMode::FLUX_INTERFACE;
                expected_tail_config.interface_topology_version =
                    tail_partition.interface_topology_version();
                expected_tail_config.interface_mask_hash =
                    tail_partition.tail_owned_mask_hash();
                expected_tail_config.interface_face_list_hash =
                    tail_partition.interface_face_list_hash();
                Options stored_physical_options = options;
                bool rebuild_stored_physical_hash = false;
                bool include_stored_return_numeric_parameters = true;
                if (preconversion_mode_switch) {
                    // Preserve the stored static-cell identity while
                    // validating every other physical option against the
                    // requested flux-interface run.
                    stored_physical_options.tail_conversion_mode =
                        TailConversionMode::STATIC_CELL;
                    expected_tail_config.conversion_mode = "static-cell";
                    rebuild_stored_physical_hash = true;
                }
                if (options.restart_allow_return_config_change) {
                    // Rebuild the expected hash with exactly the stored H10
                    // settings.  This permits controlled none<->hysteretic
                    // A/B branches in either direction while proving that
                    // every non-return physical option is unchanged.
                    // Assigning the stored hash directly would incorrectly
                    // hide unrelated physics changes.
                    stored_physical_options.tail_return_mode =
                        restored_tail.config.return_mode;
                    stored_physical_options.tail_return_energy_mev =
                        restored_tail.config.return_energy_mev;
                    stored_physical_options.tail_return_residence_steps =
                        restored_tail.config.return_residence_steps;
                    stored_physical_options.tail_return_max_stencil_radius =
                        restored_tail.config.return_max_stencil_radius;
                    stored_physical_options.tail_return_moment_tolerance =
                        restored_tail.config.return_moment_tolerance;
                    // Checkpoints written before the H10 numeric return
                    // parameters were added hashed return_mode, but did not
                    // append energy/residence/radius/tolerance.  The v2
                    // reader represents those absent fields as exact zeros.
                    // Reproduce that historical hash layout instead of
                    // hashing four synthetic zero values.
                    const bool legacy_return_hash =
                        restored_tail.config.return_mode == "none" &&
                        restored_tail.config.return_energy_mev == 0.0 &&
                        restored_tail.config.return_residence_steps == 0 &&
                        restored_tail.config.return_max_stencil_radius == 0 &&
                        restored_tail.config.return_moment_tolerance == 0.0;
                    include_stored_return_numeric_parameters =
                        !legacy_return_hash;
                    rebuild_stored_physical_hash = true;
                }
                if (options.restart_allow_dt_scale_change) {
                    // Validate the checkpoint using the explicitly stated
                    // source scale before continuing with the requested one.
                    stored_physical_options.dt_scale =
                        options.restart_source_dt_scale;
                    rebuild_stored_physical_hash = true;
                }
                if (rebuild_stored_physical_hash) {
                    expected_tail_config.physical_config_hash =
                physical_config_hash(
                    stored_physical_options, tail_partition,
                    include_stored_return_numeric_parameters,
                    stored_x_transport_velocity_table_schema != 0);
                }
                std::string config_error;
                if (!validate_vpfp_checkpoint_tail_config(
                        restored_tail.config, expected_tail_config,
                        config_error)) {
                    if (rank == 0)
                        std::cerr << "VPFP restart failed: " << config_error
                                  << "\n";
                    MPI_Finalize();
                    return 2;
                }
            }
            if (old_conversion_checkpoint &&
                options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE) {
                if (rank == 0)
                    std::cerr << "VPFP restart failed: legacy checkpoint has no "
                                 "flux-interface topology metadata\n";
                MPI_Finalize();
                return 2;
            }
            if ((!preconversion_mode_switch &&
                 restored_tail.config.conversion_mode !=
                    tail_conversion_mode_name(options.tail_conversion_mode)) ||
                restored_tail.config.flux_quadrature_order !=
                    options.tail_flux_quadrature_order ||
                restored_tail.config.flux_max_supports !=
                    options.tail_flux_max_supports ||
                restored_tail.config.flux_max_created_particles_per_step !=
                    options.tail_flux_max_created_particles_per_step) {
                if (rank == 0) {
                    std::cerr << "VPFP restart failed: tail conversion "
                                 "mode/quadrature/support configuration mismatch\n";
                }
                MPI_Finalize();
                return 2;
            }
            if (!options.restart_allow_return_config_change &&
                (restored_tail.config.return_mode != options.tail_return_mode ||
                 restored_tail.config.return_energy_mev !=
                     options.tail_return_energy_mev ||
                 restored_tail.config.return_residence_steps !=
                     options.tail_return_residence_steps ||
                 restored_tail.config.return_max_stencil_radius !=
                     options.tail_return_max_stencil_radius ||
                 restored_tail.config.return_moment_tolerance !=
                     options.tail_return_moment_tolerance)) {
                if (rank == 0) {
                    std::cerr << "VPFP restart failed: tail return configuration mismatch\n";
                }
                MPI_Finalize();
                return 2;
            }
            if (preconversion_mode_switch && rank == 0) {
                std::cout << "VPFP restart: accepted zero-ledger "
                             "static-cell-to-flux-interface mode switch\n";
            }
            if (restored_tail.config.interface_topology_metadata_present &&
                restored_tail.config.interface_topology_hash !=
                tail_partition.topology_mask_hash()) {
                if (rank == 0) {
                    std::cerr << "VPFP restart failed: tail interface topology "
                                 "hash mismatch\n";
                }
                MPI_Finalize();
                return 2;
            }
            if (restored_tail.config.control_interval !=
                    options.tail_population_control_interval ||
                restored_tail.config.target_particles_per_phase_bin !=
                    options.tail_target_particles_per_bin ||
                restored_tail.config.max_particles_per_phase_bin !=
                    options.tail_max_particles_per_bin ||
                std::fabs(restored_tail.config.max_weight_ratio -
                          options.tail_max_weight_ratio) > 1.0e-12) {
                if (rank == 0) {
                    std::cerr << "VPFP restart failed: population-control "
                                 "configuration mismatch\n";
                }
                MPI_Finalize();
                return 2;
            }
            integrator.tail_state() = restored_tail.tail;
            VpfpTailCumulativeLedger cum;
            cum.conversion_number =
                restored_tail.config.conversion_cumulative_number;
            cum.conversion_px =
                restored_tail.config.conversion_cumulative_px;
            cum.conversion_energy =
                restored_tail.config.conversion_cumulative_energy;
            cum.particles_created =
                restored_tail.config.conversion_cumulative_particles_created;
            cum.outflow_number =
                restored_tail.config.tail_cumulative_outflow_number;
            cum.control_groups =
                restored_tail.config.control_cumulative_groups;
            cum.control_fallbacks =
                restored_tail.config.control_cumulative_fallbacks;
            cum.return_number = restored_tail.config.return_cumulative_number;
            cum.return_px = restored_tail.config.return_cumulative_px;
            cum.return_jx_dx = restored_tail.config.return_cumulative_jx_dx;
            cum.return_energy = restored_tail.config.return_cumulative_energy;
            cum.return_pixx_dx = restored_tail.config.return_cumulative_pixx_dx;
            cum.return_piperp_dx = restored_tail.config.return_cumulative_piperp_dx;
            cum.return_particles_removed =
                restored_tail.config.return_cumulative_particles_removed;
            cum.return_deferred_groups =
                restored_tail.config.return_cumulative_deferred_groups;
            integrator.restore_tail_cumulative(cum);
        }
        integrator.set_step_count(restored.step);
        time = restored.time;
        step = restored.step;
    } else {
        // Initialization order from section 13.3: deposit the initial Beam
        // density, compute the initial background moments, then solve the
        // initial Poisson field exactly once before the time loop.
        if (options.beam_enabled) beam.deposit_density(grid, rank, size);
        electrons.compute_moments();
        const std::vector<double> empty_beam_density;
        const std::vector<double> empty_tail_density;
        const std::vector<double>& initial_tail_density =
            tail_enabled ? integrator.tail_state().density : empty_tail_density;
        fields.set_charge_density(
            electrons, initial_tail_density,
            options.beam_enabled ? beam.density : empty_beam_density,
            ion_density);
        field_solver.solve(fields, rank, size);
    }

    if (!options.restart_dir.empty()) {
        background_boundary.fill_ghosts(electrons, grid, rank, size);
    }
    write_startup_phase(options.output_dir, rank,
                        options.restart_dir.empty() ? "fresh_initialization_complete"
                                                    : "restart_read_complete");

    // JC5/K1 (section 9.2/9.3): when restarting from a legacy checkpoint
    // into discrete-gradient mode, re-solve the initial Poisson field once
    // so the accepted field is self-consistent with the current solver.
    // A legacy checkpoint's stored Ex_face may not match a fresh solve under
    // the discrete-gradient bootstrap check (failure 208).  The non-restart
    // path already solves the initial field (section 13.3); this makes the
    // restart path behave identically.
    if (!options.restart_dir.empty() &&
        coupling_config.mode == FieldParticleCouplingMode::DiscreteGradient) {
        if (options.beam_enabled) beam.deposit_density(grid, rank, size);
        electrons.compute_moments();
        const std::vector<double> empty_beam_density;
        const std::vector<double> empty_tail_density;
        const std::vector<double>& tail_density =
            tail_enabled ? integrator.tail_state().density : empty_tail_density;
        fields.set_charge_density(
            electrons, tail_density,
            options.beam_enabled ? beam.density : empty_beam_density,
            ion_density);
        field_solver.solve(fields, rank, size);
    }

    if (rank == 0) {
        std::ostringstream log;
        log << "Open VPFP-PIC solver: "
            << "nx_global=" << grid.nx_global
            << " length_m=" << std::scientific << std::setprecision(1)
            << grid.length()
            << " dx_m=" << grid.dx
            << std::defaultfloat << std::setprecision(17)
            << " dt=" << dt
            << " field_boundary=" << boundary_type_name(options.field_boundary.type)
            << " background_boundary="
            << boundary_side_name(options.background_boundary.left_type) << "/"
            << boundary_side_name(options.background_boundary.right_type)
            << " collision_model=" << provider->name()
            << " collision_interface_mode=" << options.collision_interface_mode
            << " bulk_collision_integrator="
            << bulk_collision_integrator_name(options.bulk_collision_integrator)
            << " field_particle_coupling=" << options.field_particle_coupling
            << " field_particle_max_iters=" << options.field_particle_max_iters
            << " field_particle_relaxation=" << options.field_particle_relaxation
            << " field_particle_field_tol=" << options.field_particle_field_tol
            << " field_particle_pairing_tol="
            << options.field_particle_pairing_tol
            << " x_transport_velocity_mode="
            << x_transport_velocity_mode_name(options.x_transport_velocity_mode)
            << " x_transport_velocity_table_schema=1"
            << " restart_allow_field_particle_coupling_change="
            << (options.restart_allow_field_particle_coupling_change ? 1 : 0)
            << " collision_induced_conversion="
            << (options.collision_model != "none" &&
                options.collision_interface_mode != "zero-wall-validation" ? 1 : 0)
            << " background_phase_space_mode="
            << background_phase_space_mode_name(options.background_phase_space_mode)
            << " time_integrator="
            << background_phase_space_mode_name(options.background_phase_space_mode)
            << " remap_type=conservative-csl-split"
            << " beam_enabled=" << (options.beam_enabled ? 1 : 0)
            << " background_tail_mode=" << options.background_tail_mode
            << " tail_convert_energy_mev=" << options.tail_convert_energy_mev
            << " tail_conversion_bins=" << options.tail_conversion_upar_bins
            << "x" << options.tail_conversion_energy_bins
            << " tail_conversion_strategy="
            << bulk_tail_loading_policy_name(tail_converter.loading_policy())
            << " tail_subcell_loading="
            << (tail_converter.subcell_loading_enabled() ? 1 : 0)
            << " tail_conversion_mode="
            << tail_conversion_mode_name(options.tail_conversion_mode)
            << " observe_upar_flux="
            << ((options.tail_conversion_mode == TailConversionMode::FLUX_AUDIT ||
                 options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE) ? 1 : 0)
            << " apply_upar_sink="
            << (options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE ? 1 : 0)
            << " observe_collision_flux="
            << ((options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE &&
                 options.collision_interface_mode == "exporting-absorbing" &&
                 options.collision_model != "none" &&
                 options.bulk_collision_integrator ==
                     BulkCollisionIntegrator::CHANG_COOPER_FLUX) ? 1 : 0)
            << " tail_flux_quadrature_order="
            << options.tail_flux_quadrature_order
            << " tail_flux_max_supports=" << options.tail_flux_max_supports
             << " tail_flux_max_created_particles_per_step="
             << options.tail_flux_max_created_particles_per_step
             << " initial_tail_conversion_number="
             << initial_tail_conversion.number_created
             << " initial_tail_conversion_particles="
             << initial_tail_conversion.particles_created
             << " tail_conversion_strategy_stage=cell_volume_supported"
            << " tail_population_control_interval="
            << options.tail_population_control_interval
            << " tail_target_particles_per_bin="
            << options.tail_target_particles_per_bin
            << " tail_max_particles_per_bin="
            << options.tail_max_particles_per_bin
            << " tail_max_weight_ratio=" << options.tail_max_weight_ratio
            << " tail_return_mode=" << options.tail_return_mode
            << " tail_return_energy_mev=" << options.tail_return_energy_mev
            << " tail_return_residence_steps="
            << options.tail_return_residence_steps
            << " tail_return_max_stencil_radius="
            << options.tail_return_max_stencil_radius
            << " tail_return_moment_tolerance="
            << options.tail_return_moment_tolerance
            << " tail_return_projection_schema="
            << tail_bulk_return_projection_schema()
            << " restart_allow_return_config_change="
            << (options.restart_allow_return_config_change ? 1 : 0)
            << " restart_allow_dt_scale_change="
            << (options.restart_allow_dt_scale_change ? 1 : 0)
            << " restart_source_dt_scale="
            << options.restart_source_dt_scale
            << " tail_collision_kernel=" << options.tail_collision_kernel
            << " tail_tail_backend="
            << tail_collision_kernel_name(collision_config.tail_tail_kernel)
            << " tail_bulk_backend="
            << tail_collision_kernel_name(collision_config.tail_bulk_kernel)
            << " tail_collision_weight_mode="
            << tail_collision_weight_mode_name(collision_config.weight_mode)
            << " tail_collision_weight_algorithm="
            << (collision_config.weight_mode ==
                        TailCollisionWeightMode::VirtualSplit
                    ? "sentoku-kemp-bounded-v1"
                    : "equal-strata-exact-v1")
            << " tail_collision_max_substeps="
            << collision_config.max_substeps
            << " tail_collision_max_particle_growth="
            << options.tail_collision_max_particle_growth
            << " collision_pairs="
            << (collision_config.pairs.bulk_bulk ? 1 : 0)
            << (collision_config.pairs.bulk_tail ? 1 : 0)
            << (collision_config.pairs.tail_bulk ? 1 : 0)
            << (collision_config.pairs.tail_tail ? 1 : 0)
            << " checkpoint_schema=vpfp-open-v4"
            << " snapshot_schema=vpfp-open-csl-v1\n";
        std::cout << log.str();
        const double length_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                        std::max(options.length, 1.0);
        if (std::fabs(grid.length() - options.length) > length_tolerance) {
            std::cerr << "grid length does not close to 64*eps*L: "
                      << grid.length() << " vs " << options.length << "\n";
            MPI_Finalize();
            return 2;
        }
    }

    // Stage-H6: checkpoint tail-state builder (section 12.1).  Refreshed at
    // every checkpoint write (on demand) so the accepted cumulative ledgers
    // and combined checksums are current.  Never copied on steps that do not
    // write a checkpoint.
    VpfpCheckpointTailState checkpoint_tail_storage;
    if (tail_enabled) checkpoint_tail_storage.present = true;
    const VpfpCheckpointTailState* make_checkpoint_tail = NULL;
    auto build_checkpoint_tail = [&]() -> const VpfpCheckpointTailState*
    {
        if (!tail_enabled) return NULL;
        VpfpCheckpointTailConfig& c = checkpoint_tail_storage.config;
        c.partition_config_hash = tail_partition.config_hash;
        c.physical_config_hash = physical_hash;
        c.diagnostic_config_hash = diagnostic_hash;
        // Section 7.11.4 branch B: carry the explicit threshold-aware
        // energy edges and their hash into the checkpoint manifest so a
        // restart with a different grouping is rejected (the binary header
        // validation is via partition_config_hash, which folds the edges).
        c.conversion_energy_edges = tail_partition.conversion_energy_edges;
        c.conversion_energy_edges_hash =
            tail_partition.conversion_energy_edges_hash();
        c.convert_energy_mev = options.tail_convert_energy_mev;
        c.buffer_width_mev = 1.0;
        c.upar_bins = options.tail_conversion_upar_bins;
        c.energy_bins = options.tail_conversion_energy_bins;
        c.return_mode = options.tail_return_mode;
        c.return_energy_mev = options.tail_return_energy_mev;
        c.return_residence_steps = options.tail_return_residence_steps;
        c.return_max_stencil_radius = options.tail_return_max_stencil_radius;
        c.return_moment_tolerance = options.tail_return_moment_tolerance;
        c.collision_kernel = options.tail_collision_kernel;
        c.collision_weight_mode = options.tail_collision_weight_mode;
        c.collision_max_substeps = options.tail_collision_max_substeps;
        c.collision_max_particle_growth =
            options.tail_collision_max_particle_growth;
        c.population_control_enabled = population_config.enabled;
        c.control_interval = population_config.control_interval;
        c.target_particles_per_phase_bin =
            population_config.target_particles_per_phase_bin;
        c.max_particles_per_phase_bin =
            population_config.max_particles_per_phase_bin;
        c.max_weight_ratio = population_config.max_weight_ratio;
        c.max_support = population_config.max_support;
        c.conversion_mode = tail_conversion_mode_name(options.tail_conversion_mode);
        c.collision_interface_mode = options.collision_interface_mode;
        c.bulk_collision_integrator =
            bulk_collision_integrator_name(options.bulk_collision_integrator);
        c.collision_induced_conversion =
            options.collision_model != "none" &&
            options.collision_interface_mode != "zero-wall-validation";
        c.flux_quadrature_order = options.tail_flux_quadrature_order;
        c.flux_max_supports = options.tail_flux_max_supports;
        c.flux_max_created_particles_per_step =
            options.tail_flux_max_created_particles_per_step;
        c.interface_topology_hash = tail_partition.topology_mask_hash();
        c.interface_topology_version =
            tail_partition.interface_topology_version();
        c.interface_mask_hash = tail_partition.tail_owned_mask_hash();
        c.interface_face_list_hash =
            tail_partition.interface_face_list_hash();
        c.interface_topology_metadata_present =
            options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE;
        c.interface_topology_detail_metadata_present =
            options.tail_conversion_mode == TailConversionMode::FLUX_INTERFACE;
        const VpfpTailCumulativeLedger& cum = integrator.tail_cumulative();
        c.conversion_cumulative_number = cum.conversion_number;
        c.conversion_cumulative_px = cum.conversion_px;
        c.conversion_cumulative_energy = cum.conversion_energy;
        c.conversion_cumulative_particles_created = cum.particles_created;
        c.tail_cumulative_outflow_number = cum.outflow_number;
        c.control_cumulative_groups = cum.control_groups;
        c.control_cumulative_fallbacks = cum.control_fallbacks;
        c.return_cumulative_number = cum.return_number;
        c.return_cumulative_px = cum.return_px;
        c.return_cumulative_jx_dx = cum.return_jx_dx;
        c.return_cumulative_energy = cum.return_energy;
        c.return_cumulative_pixx_dx = cum.return_pixx_dx;
        c.return_cumulative_piperp_dx = cum.return_piperp_dx;
        c.return_cumulative_particles_removed = cum.return_particles_removed;
        c.return_cumulative_deferred_groups = cum.return_deferred_groups;
        const VpfpCombinedChecksum& chk = integrator.combined_checksum();
        c.combined_number = chk.number;
        c.combined_kinetic_energy = chk.kinetic_energy;
        c.combined_field_energy = chk.field_energy;
        checkpoint_tail_storage.tail = integrator.tail_state();
        return &checkpoint_tail_storage;
    };
    double last_diagnostics_wall = 0.0;

    size_t next_checkpoint_time = 0;
    size_t next_snapshot_time = 0;
    while (options.stop_after_steps <= 0 || step < options.stop_after_steps) {
        if (!VpfpTimeControl::should_advance(time, options.stop_time, dt)) break;
        const double step_dt = std::min(dt, options.stop_time - time);
        VpfpStepResult result = integrator.advance(electrons, beam, fields, ion_density,
                                                   time, step_dt, rank, size);
        if (!result.accepted) {
            diagnostics.write_failure(step, time, result, rank);
            if (rank == 0) {
                std::cerr << "VPFP step rejected, failure_code="
                          << result.failure_code << " ("
                          << failure_code_name(result.failure_code) << ")"
                          << " failure_stage=" << result.failure_stage
                          << " failing_rank=" << result.failing_rank
                          << " failing_ix=" << result.failing_ix
                          << " failing_iupar=" << result.failing_iupar
                          << " failing_iuperp=" << result.failing_iuperp
                          << " input_min=" << result.input_min
                          << " input_max=" << result.input_max
                          << " output_min=" << result.output_min
                          << " output_max=" << result.output_max
                          << " first_nonfinite_value="
                          << result.first_nonfinite_value
                          << " audit_valid=" << result.audit_valid
                          << " audit_failure_code="
                          << result.audit_failure_code << "\n";
            }
            MPI_Finalize();
            return 3;
        }
        time += step_dt;
        ++step;
        result.wall_diagnostics_seconds = last_diagnostics_wall;
        const std::chrono::steady_clock::time_point diag_begin =
            std::chrono::steady_clock::now();
        diagnostics.write_accepted_step(
            step, time, result, field_solver.diagnostics(), electrons,
            tail_enabled ? &integrator.tail_state() : NULL, grid, rank,
            result.pairing_audit_enabled ? &integrator.pairing_workspace()
                                         : NULL);
        diagnostics.write_population_control(step, time, result, rank);
        make_checkpoint_tail = NULL;
        if (!options.checkpoint_dir.empty() && options.checkpoint_every > 0 &&
            step % options.checkpoint_every == 0) {
            make_checkpoint_tail = build_checkpoint_tail();
            VpfpCheckpointControl checkpoint = { step, time, step_dt };
            VpfpCouplingManifestConfig coupling_manifest;
            coupling_manifest.mode = options.field_particle_coupling;
            coupling_manifest.max_iters = options.field_particle_max_iters;
            coupling_manifest.relaxation = options.field_particle_relaxation;
            coupling_manifest.field_tol = options.field_particle_field_tol;
            coupling_manifest.pairing_tol = options.field_particle_pairing_tol;
            coupling_manifest.background_phase_space_mode =
                background_phase_space_mode_name(options.background_phase_space_mode);
            coupling_manifest.x_transport_velocity_mode =
                x_transport_velocity_mode_name(options.x_transport_velocity_mode);
            coupling_manifest.x_transport_velocity_table_schema = 1;
            std::string checkpoint_error;
            if (!write_vpfp_checkpoint(options.checkpoint_dir, checkpoint, electrons, beam,
                                       fields, grid, options.field_boundary,
                                       options.background_boundary, provider->name(),
                                       make_checkpoint_tail, coupling_manifest,
                                       rank, size,
                                       checkpoint_error)) {
                if (rank == 0) std::cerr << "VPFP checkpoint failed: " << checkpoint_error << "\n";
                MPI_Finalize(); return 4;
            }
        }
        while (next_checkpoint_time < options.checkpoint_times.size()) {
            const double target =
                options.checkpoint_times[next_checkpoint_time] * Const::femto;
            if (!VpfpTimeControl::target_reached(time, target, dt)) break;
            const double tolerance = VpfpTimeControl::roundoff_tolerance(
                time, target, dt);
            // Save only on the first full accepted state at/after the target
            // physical time (the state whose step interval crosses it).
            if (time - step_dt < target + tolerance) {
                make_checkpoint_tail = build_checkpoint_tail();
                const std::string directory =
                    checkpoint_time_directory(options.output_dir, target,
                                              time, step);
                VpfpCheckpointControl checkpoint = { step, time, step_dt };
                VpfpCouplingManifestConfig coupling_manifest;
                coupling_manifest.mode = options.field_particle_coupling;
                coupling_manifest.max_iters = options.field_particle_max_iters;
                coupling_manifest.relaxation = options.field_particle_relaxation;
                coupling_manifest.field_tol = options.field_particle_field_tol;
                coupling_manifest.pairing_tol = options.field_particle_pairing_tol;
                coupling_manifest.background_phase_space_mode =
                    background_phase_space_mode_name(options.background_phase_space_mode);
                coupling_manifest.x_transport_velocity_mode =
                    x_transport_velocity_mode_name(options.x_transport_velocity_mode);
                coupling_manifest.x_transport_velocity_table_schema = 1;
                std::string checkpoint_error;
                if (!write_vpfp_checkpoint(directory, checkpoint, electrons, beam,
                                           fields, grid, options.field_boundary,
                                           options.background_boundary, provider->name(),
                                           make_checkpoint_tail, coupling_manifest,
                                           rank, size,
                                           checkpoint_error)) {
                    if (rank == 0) std::cerr << "VPFP time checkpoint failed: "
                                             << checkpoint_error << "\n";
                    MPI_Finalize(); return 4;
                }
            }
            ++next_checkpoint_time;
        }
        while (next_snapshot_time < options.snapshot_times.size()) {
            const double target =
                options.snapshot_times[next_snapshot_time] * Const::femto;
            if (!VpfpTimeControl::target_reached(time, target, dt)) break;
            const double tolerance = VpfpTimeControl::roundoff_tolerance(
                time, target, dt);
            // Snapshot the first full accepted state at/after the target
            // physical time (same crossing rule as checkpoint-times), so
            // dt-scale 1.0 and 0.5 runs land on the same physical times.
            if (time - step_dt < target + tolerance) {
                diagnostics.write_snapshot(step, time, electrons, beam,
                                           fields, grid,
                                           tail_enabled
                                               ? &integrator.tail_state() : NULL,
                                           tail_enabled ? &tail_partition : NULL,
                                           options.tail_convert_energy_mev,
                                           rank, size);
            }
            ++next_snapshot_time;
        }
        last_diagnostics_wall = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - diag_begin).count();
    }
    MPI_Finalize();
    return 0;
}
