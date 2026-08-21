// Gate F 10.5.2: fixed accepted-state dt/dt/2/dt/4 audit.  This executable
// reads a production checkpoint, advances only copies through the production
// field-particle split, and writes no checkpoint or production diagnostics.

#include "bulk_tail_converter.h"
#include "collision_coefficients.h"
#include "cylindrical_fp_collision.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "vpfp_checkpoint.h"
#include "vpfp_integrator.h"

#include <mpi.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

struct Options {
    std::string restart_dir;
    std::string output_dir;
    std::vector<double> scales;
    std::vector<int> intervals;
    double restart_source_dt_scale;
    std::string mode;
    Options()
        : scales(), intervals(), restart_source_dt_scale(1.0),
          mode("deterministic-field-particle")
    {}
};

struct IntervalRecord {
    double scale;
    int interval;
    double pair_signed;
    double pair_abs;
    double bulk_work;
    double tail_work;
    double beam_work;
    double field_energy_change;
    double electrode_work;
    double bulk_norm;
    double tail_norm;
    double beam_norm;
    double field_norm;
    double conversion_number;
    double conversion_px;
    double conversion_energy;
    unsigned long long conversion_event_count;
    unsigned long long tail_count;
    unsigned long long tail_rng_counter;
    unsigned long long beam_rng_counter;
    unsigned long long rng_hash;
    int accepted;
    int finite;
    int poisson_identity_pass;
    double poisson_residual;
    double poisson_bulk_work;
    double poisson_tail_work;
    double poisson_beam_work;
    double poisson_ion_work;
    double poisson_boundary_source_work;
    double poisson_component_reconstruction_error;
};

struct SourceOwnershipAudit {
    int count[7];
    double residual[7];
    SourceOwnershipAudit()
    {
        for (int i = 0; i < 7; ++i) {
            count[i] = 1;
            residual[i] = 0.0;
        }
    }
};

struct AuditRunStatus {
    bool same_initial_state;
    bool same_physical_window;
    bool source_ownership_valid;
    bool collision_frozen;
    bool conversion_frozen;
    bool h10_frozen;
    AuditRunStatus()
        : same_initial_state(true), same_physical_window(true),
          source_ownership_valid(true), collision_frozen(true),
          conversion_frozen(true), h10_frozen(true)
    {}
};

struct AuditBoundaryConfig {
    ElectrostaticBoundary field;
    OpenBackgroundBoundaryConfig background;
};

bool parse_double_list(const std::string& text, std::vector<double>& values)
{
    values.clear();
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        char* end = NULL;
        const double value = std::strtod(item.c_str(), &end);
        if (end == item.c_str() || *end != '\0' || !(value > 0.0)) return false;
        values.push_back(value);
    }
    return !values.empty();
}

bool parse_int_list(const std::string& text, std::vector<int>& values)
{
    values.clear();
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        char* end = NULL;
        const long value = std::strtol(item.c_str(), &end, 10);
        if (end == item.c_str() || *end != '\0' || value <= 0) return false;
        values.push_back(static_cast<int>(value));
    }
    return !values.empty();
}

bool parse_args(int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--restart-dir" && i + 1 < argc) {
            options.restart_dir = argv[++i];
        } else if (arg == "--restart-source-dt-scale" && i + 1 < argc) {
            options.restart_source_dt_scale = std::strtod(argv[++i], NULL);
        } else if (arg == "--mode" && i + 1 < argc) {
            options.mode = argv[++i];
        } else if (arg == "--dt-scales" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], options.scales)) return false;
        } else if (arg == "--common-dt-intervals" && i + 1 < argc) {
            if (!parse_int_list(argv[++i], options.intervals)) return false;
        } else if (arg == "--output-dir" && i + 1 < argc) {
            options.output_dir = argv[++i];
        } else {
            return false;
        }
    }
    return !options.restart_dir.empty() && !options.output_dir.empty() &&
        options.mode == "deterministic-field-particle" &&
        std::fabs(options.restart_source_dt_scale - 1.0) <= 1.0e-14 &&
        options.scales.size() == 3 && options.intervals.size() != 0;
}

AuditBoundaryConfig read_checkpoint_boundary(const std::string& directory,
                                             int rank)
{
    AuditBoundaryConfig config = {};
    config.field = { ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    config.background.left_type = BackgroundXBoundaryType::RESERVOIR;
    config.background.right_type = BackgroundXBoundaryType::RESERVOIR;
    config.background.left_reservoir =
        MaxwellianReservoir { Param::dens, Param::temperature_e, 0.0 };
    config.background.right_reservoir =
        MaxwellianReservoir { Param::dens, Param::temperature_e, 0.0 };
    double packed[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int types[3] = { static_cast<int>(config.field.type),
                     static_cast<int>(config.background.left_type),
                     static_cast<int>(config.background.right_type) };
    if (rank == 0) {
        std::ifstream input((directory + "/manifest.txt").c_str());
        std::string key;
        while (input >> key) {
            if (key == "field_boundary") input >> types[0];
            else if (key == "background_left_type") input >> types[1];
            else if (key == "background_right_type") input >> types[2];
            else if (key == "e_left") input >> packed[0];
            else if (key == "phi_left") input >> packed[1];
            else if (key == "phi_right") input >> packed[2];
            else { std::string ignored; input >> ignored; }
        }
    }
    MPI_Bcast(types, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(packed, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    config.field.type = static_cast<ElectrostaticBoundaryType>(types[0]);
    config.field.e_left = packed[0];
    config.field.phi_left = packed[1];
    config.field.phi_right = packed[2];
    config.background.left_type = static_cast<BackgroundXBoundaryType>(types[1]);
    config.background.right_type = static_cast<BackgroundXBoundaryType>(types[2]);
    return config;
}

bool ensure_directory(const std::string& path)
{
    if (path.empty()) return false;
    std::string current;
    for (size_t i = 0; i < path.size(); ++i) {
        current.push_back(path[i]);
        const bool boundary = path[i] == '/' || path[i] == '\\';
        if (!boundary || current.size() == 1) continue;
#ifdef _WIN32
        if (_mkdir(current.c_str()) != 0 && errno != EEXIST) return false;
#else
        if (mkdir(current.c_str(), 0777) != 0 && errno != EEXIST) return false;
#endif
    }
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

void initialize_progress_file(const std::string& output_dir, int rank,
                              double checkpoint_time, int checkpoint_step)
{
    if (rank != 0) return;
    if (!ensure_directory(output_dir)) {
        std::cerr << "cannot create audit output directory: " << output_dir << "\n";
        return;
    }
    std::ofstream progress((output_dir + "/audit_progress.dat").c_str());
    progress << std::setprecision(17)
             << "# Gate-F 10.5 fixed-state audit progress\n"
             << "checkpoint_time_s=" << checkpoint_time << "\n"
             << "checkpoint_step=" << checkpoint_step << "\n"
             << "dt_scale substep completed_substeps physical_time_s "
             << "completed_common_intervals status\n";
    progress.flush();
}

void append_progress(const std::string& output_dir, int rank, double scale,
                     int substep, int total_substeps, double physical_time,
                     int completed_intervals, const char* status)
{
    if (rank != 0) return;
    std::ofstream progress((output_dir + "/audit_progress.dat").c_str(),
                           std::ios::app);
    progress << std::setprecision(17) << scale << " " << substep << " "
             << total_substeps << " " << physical_time << " "
             << completed_intervals << " " << status << "\n";
    progress.flush();
}

std::uint64_t fnv_bytes(std::uint64_t hash, const void* data, size_t count)
{
    const unsigned char* bytes = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < count; ++i) {
        hash ^= static_cast<std::uint64_t>(bytes[i]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

template <typename T>
std::uint64_t hash_vector(std::uint64_t hash, const std::vector<T>& values)
{
    const std::uint64_t size = static_cast<std::uint64_t>(values.size());
    hash = fnv_bytes(hash, &size, sizeof(size));
    return values.empty() ? hash : fnv_bytes(hash, &values[0],
                                             values.size() * sizeof(T));
}

std::uint64_t state_hash(const Species& bulk, const BeamPIC& beam,
                         const EMFields& fields, const BackgroundTailPIC& tail)
{
    std::uint64_t hash = 1469598103934665603ULL;
    hash = hash_vector(hash, bulk.f);
    hash = hash_vector(hash, beam.particles);
    hash = hash_vector(hash, beam.density);
    const BeamPersistentState beam_state = beam.export_persistent_state();
    hash = fnv_bytes(hash, &beam_state, sizeof(beam_state));
    hash = hash_vector(hash, fields.Ex_face);
    hash = hash_vector(hash, fields.phi);
    hash = hash_vector(hash, fields.rho);
    BackgroundTailStateSnapshot snapshot;
    tail.export_accepted_state(snapshot);
    hash = hash_vector(hash, snapshot.particles);
    hash = fnv_bytes(hash, &snapshot.id_counter, sizeof(snapshot.id_counter));
    hash = fnv_bytes(hash, &snapshot.collision_rng_seed,
                     sizeof(snapshot.collision_rng_seed));
    return hash;
}

std::uint64_t rng_state_hash(const BeamPIC& beam,
                             const BackgroundTailPIC& tail)
{
    const BeamPersistentState beam_state = beam.export_persistent_state();
    BackgroundTailStateSnapshot tail_state;
    tail.export_accepted_state(tail_state);
    std::uint64_t hash = 1469598103934665603ULL;
    hash = fnv_bytes(hash, &beam_state.rng_state, sizeof(beam_state.rng_state));
    return fnv_bytes(hash, &tail_state.collision_rng_seed,
                     sizeof(tail_state.collision_rng_seed));
}

std::uint64_t checkpoint_ledger_hash(const VpfpCheckpointTailState& tail)
{
    const VpfpCheckpointTailConfig& c = tail.config;
    std::uint64_t hash = 1469598103934665603ULL;
    hash = fnv_bytes(hash, &c.conversion_cumulative_number,
                     sizeof(c.conversion_cumulative_number));
    hash = fnv_bytes(hash, &c.conversion_cumulative_px,
                     sizeof(c.conversion_cumulative_px));
    hash = fnv_bytes(hash, &c.conversion_cumulative_energy,
                     sizeof(c.conversion_cumulative_energy));
    hash = fnv_bytes(hash, &c.conversion_cumulative_particles_created,
                     sizeof(c.conversion_cumulative_particles_created));
    hash = fnv_bytes(hash, &c.return_cumulative_number,
                     sizeof(c.return_cumulative_number));
    return fnv_bytes(hash, &c.return_cumulative_energy,
                     sizeof(c.return_cumulative_energy));
}

double global_sum(double value)
{
    double global = 0.0;
    MPI_Allreduce(&value, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

unsigned long long global_sum_u64(unsigned long long value)
{
    unsigned long long global = 0;
    MPI_Allreduce(&value, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    return global;
}

unsigned long long global_xor_u64(unsigned long long value)
{
    unsigned long long global = 0;
    MPI_Allreduce(&value, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                  MPI_COMM_WORLD);
    return global;
}

double local_bulk_norm(const Species& bulk)
{
    double value = 0.0;
    for (size_t i = 0; i < bulk.f.size(); ++i) value += bulk.f[i] * bulk.f[i];
    return value;
}

double local_field_norm(const EMFields& fields)
{
    double value = 0.0;
    for (size_t i = 0; i < fields.Ex_face.size(); ++i)
        value += fields.Ex_face[i] * fields.Ex_face[i];
    return value;
}

double local_tail_norm(const BackgroundTailPIC& tail)
{
    double value = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i)
        value += tail.particles[i].weight * tail.particles[i].weight;
    return value;
}

double local_beam_norm(const BeamPIC& beam)
{
    double value = 0.0;
    for (size_t i = 0; i < beam.particles.size(); ++i)
        value += beam.particles[i].weight * beam.particles[i].weight;
    return value;
}

std::vector<double> component_charge(Species& bulk,
                                     const BackgroundTailPIC& tail,
                                     const BeamPIC& beam,
                                     const SpatialGrid& grid,
                                     int component)
{
    bulk.compute_moments();
    std::vector<double> charge(static_cast<size_t>(grid.nx_total), 0.0);
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const size_t cell = static_cast<size_t>(grid.nghost + ix);
        // All moment/deposit arrays contain physical local cells only.  The
        // output work array includes ghosts because evaluate_work_identity()
        // consumes the field layout.  Keep those two index spaces explicit.
        if (component == 0 &&
            static_cast<size_t>(ix) < bulk.charge_density.size())
            charge[cell] = bulk.charge_density[static_cast<size_t>(ix)];
        else if (component == 1 &&
                 static_cast<size_t>(ix) < tail.density.size())
            charge[cell] = -Const::qe * tail.density[static_cast<size_t>(ix)];
        else if (component == 2 &&
                 static_cast<size_t>(ix) < beam.density.size())
            charge[cell] = -Const::qe * beam.density[static_cast<size_t>(ix)];
    }
    return charge;
}

std::vector<double> subtract_charge(const std::vector<double>& after,
                                    const std::vector<double>& before)
{
    std::vector<double> delta(after.size(), 0.0);
    for (size_t i = 0; i < delta.size(); ++i) delta[i] = after[i] - before[i];
    return delta;
}

void stage_work(const VpfpStepResult& step, double& bulk, double& tail,
                double& beam)
{
    bulk = 0.0;
    tail = 0.0;
    beam = 0.0;
    for (int i = 0; i < step.stage_energy_count; ++i) {
        const VpfpStageEnergyRecord& record = step.stage_energy[i];
        if (record.stage_id == VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK) {
            bulk = record.bulk_upar_face_work +
                record.bulk_upar_velocity_boundary_work;
            tail = record.tail_kick_work;
            beam = record.beam_kick_work;
            return;
        }
    }
}

double source_value(const VpfpStageEnergyRecord& record, int source)
{
    switch (source) {
    case 0: return record.electrostatic_boundary_work;
    case 1: return record.background_left_inflow_energy -
                   record.background_left_outflow_energy +
                   record.background_right_inflow_energy -
                   record.background_right_outflow_energy;
    case 2: return record.beam_injected_energy - record.beam_outflow_energy;
    case 3: return record.bulk_upar_velocity_boundary_work;
    case 4: return record.conversion_energy;
    case 5: return record.tail_return_energy;
    case 6: return record.collision_reservoir_energy;
    default: return 0.0;
    }
}

bool source_ownership_ok(const VpfpStepResult& step,
                         SourceOwnershipAudit& audit)
{
    if (step.stage_energy_count != VPFP_STAGE_ENERGY_RECORD_COUNT) return false;
    bool valid = true;
    for (int source = 0; source != 7; ++source) {
        const double endpoint = source_value(step.stage_energy[
            VPFP_STAGE_ENERGY_RECORD_COUNT - 1], source) -
            source_value(step.stage_energy[0], source);
        long double increments = 0.0L;
        long double scale = 1.0L + std::fabs(endpoint);
        for (int i = 1; i < VPFP_STAGE_ENERGY_RECORD_COUNT; ++i) {
            const double delta = source_value(step.stage_energy[i], source) -
                source_value(step.stage_energy[i - 1], source);
            increments += static_cast<long double>(delta);
            scale += std::fabs(delta);
        }
        const double residual = static_cast<double>(increments) - endpoint;
        audit.residual[source] = std::max(audit.residual[source],
                                          std::fabs(residual));
        if (std::fabs(residual) > 8192.0 * std::numeric_limits<double>::epsilon() *
                static_cast<double>(scale)) {
            audit.count[source] = 0;
            valid = false;
        }
    }
    return valid;
}

VpfpTailCumulativeLedger checkpoint_tail_cumulative(
    const VpfpCheckpointTailState& tail)
{
    VpfpTailCumulativeLedger ledger;
    ledger.conversion_number = tail.config.conversion_cumulative_number;
    ledger.conversion_px = tail.config.conversion_cumulative_px;
    ledger.conversion_energy = tail.config.conversion_cumulative_energy;
    ledger.particles_created = tail.config.conversion_cumulative_particles_created;
    ledger.outflow_number = tail.config.tail_cumulative_outflow_number;
    ledger.control_groups = tail.config.control_cumulative_groups;
    ledger.control_fallbacks = tail.config.control_cumulative_fallbacks;
    ledger.return_number = tail.config.return_cumulative_number;
    ledger.return_px = tail.config.return_cumulative_px;
    ledger.return_jx_dx = tail.config.return_cumulative_jx_dx;
    ledger.return_energy = tail.config.return_cumulative_energy;
    ledger.return_pixx_dx = tail.config.return_cumulative_pixx_dx;
    ledger.return_piperp_dx = tail.config.return_cumulative_piperp_dx;
    ledger.return_particles_removed = tail.config.return_cumulative_particles_removed;
    ledger.return_deferred_groups = tail.config.return_cumulative_deferred_groups;
    return ledger;
}

bool run_scale(const SpatialGrid& grid, const VpfpCheckpointControl& control,
               const Species& bulk_initial, const BeamPIC& beam_initial,
               const EMFields& fields_initial,
               const VpfpCheckpointTailState& tail_initial,
               const AuditBoundaryConfig& boundary_config,
               double scale, const std::vector<int>& intervals,
               const std::string& output_dir, int rank, int size,
               std::vector<IntervalRecord>& output,
               SourceOwnershipAudit& ownership,
               AuditRunStatus& audit_status,
               std::uint64_t expected_initial_hash)
{
    Species bulk = bulk_initial;
    BeamPIC beam = beam_initial;
    EMFields fields = fields_initial;
    OpenBackgroundBoundary background_boundary(boundary_config.background);
    OpenElectrostaticSolver field_solver;
    field_solver.init(grid, boundary_config.field);
    ZeroCollisionCoefficients zero_collision;
    CylindricalFokkerPlanckCollision collision(
        zero_collision, CollisionIntegratorType::BACKWARD_EULER);

    HybridVelocityPartition partition;
    BulkTailConverter converter;
    const bool tail_on = tail_initial.present;
    if (tail_on) {
        partition.init(bulk.cgrid, tail_initial.config.convert_energy_mev,
                       std::max(1.0e-12, tail_initial.config.buffer_width_mev),
                       std::max(1, tail_initial.config.upar_bins),
                       std::max(1, tail_initial.config.energy_bins));
    }
    std::unique_ptr<VpfpIntegrator> integrator;
    if (tail_on) {
        integrator.reset(new VpfpIntegrator(background_boundary, field_solver,
                                             collision, partition, converter,
                                             true));
    } else {
        integrator.reset(new VpfpIntegrator(background_boundary, field_solver,
                                             collision));
    }
    integrator->init(grid);
    integrator->set_beam_enabled(true);
    integrator->set_stage_energy_audit_enabled(true);
    integrator->set_step_count(control.step);
    if (tail_on) {
        integrator->tail_state() = tail_initial.tail;
        integrator->restore_tail_cumulative(checkpoint_tail_cumulative(tail_initial));
    }
    audit_status.same_initial_state = audit_status.same_initial_state &&
        state_hash(bulk, beam, fields, integrator->tail_state()) == expected_initial_hash;

    const int max_interval = *std::max_element(intervals.begin(), intervals.end());
    const int steps_per_interval = static_cast<int>(std::llround(1.0 / scale));
    if (steps_per_interval <= 0 || std::fabs(scale * steps_per_interval - 1.0) > 1.0e-12)
        return false;
    const double dt = control.dt * scale;
    audit_status.same_physical_window = audit_status.same_physical_window &&
        std::fabs(dt * static_cast<double>(steps_per_interval) - control.dt) <=
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::fabs(control.dt));
    std::vector<double> ion_density(static_cast<size_t>(grid.nx_total), Param::dens);
    EMFields interval_start = fields;
    field_solver.reconstruct_phi(interval_start, rank, size);
    std::vector<double> bulk_charge_start = component_charge(
        bulk, integrator->tail_state(), beam, grid, 0);
    std::vector<double> tail_charge_start = component_charge(
        bulk, integrator->tail_state(), beam, grid, 1);
    std::vector<double> beam_charge_start = component_charge(
        bulk, integrator->tail_state(), beam, grid, 2);
    double pair_signed_sum = 0.0;
    double pair_abs_sum = 0.0;
    double bulk_work_sum = 0.0;
    double tail_work_sum = 0.0;
    double beam_work_sum = 0.0;
    double field_energy_change_sum = 0.0;
    double electrode_work_sum = 0.0;
    double conversion_number_sum = 0.0;
    double conversion_px_sum = 0.0;
    double conversion_energy_sum = 0.0;
    unsigned long long conversion_event_count_sum = 0;
    const int total_substeps = max_interval * steps_per_interval;
    for (int step_index = 1; step_index <= max_interval * steps_per_interval;
         ++step_index) {
        const double time = control.time + static_cast<double>(step_index - 1) * dt;
        const VpfpStepResult step = integrator->advance_fixed_state_field_particle_audit(
            bulk, beam, fields, ion_density, time, dt, rank, size);
        if (!step.accepted || !step.finite) {
            append_progress(output_dir, rank, scale, step_index, total_substeps,
                            time + dt, (step_index - 1) / steps_per_interval,
                            "FAILED");
            return false;
        }
        audit_status.source_ownership_valid =
            source_ownership_ok(step, ownership) && audit_status.source_ownership_valid;
        double step_bulk_work = 0.0;
        double step_tail_work = 0.0;
        double step_beam_work = 0.0;
        stage_work(step, step_bulk_work, step_tail_work, step_beam_work);
        const double step_pair = step.ledger.domain_energy_change -
            step.ledger.electrostatic_boundary_work -
            step.ledger.background_boundary_energy_net -
            step.ledger.beam_boundary_energy_net;
        pair_signed_sum += step_pair;
        pair_abs_sum += std::fabs(step_pair);
        bulk_work_sum += step_bulk_work;
        tail_work_sum += step_tail_work;
        beam_work_sum += step_beam_work;
        field_energy_change_sum += step.ledger.field_energy -
            step.ledger.field_energy_before;
        electrode_work_sum += step.ledger.electrostatic_boundary_work;
        conversion_number_sum += step.ledger.conversion_number_removed;
        conversion_px_sum += step.ledger.conversion_px_removed;
        conversion_energy_sum += step.ledger.conversion_energy_removed;
        conversion_event_count_sum +=
            static_cast<unsigned long long>(step.conversion_events.size());
        const double freeze_scale = std::max(1.0,
            std::fabs(step.ledger.domain_energy_change));
        const double freeze_tolerance = 8192.0 *
            std::numeric_limits<double>::epsilon() * freeze_scale;
        audit_status.collision_frozen = audit_status.collision_frozen &&
            std::fabs(step.ledger.collision_reservoir_energy) <= freeze_tolerance;
        audit_status.conversion_frozen = audit_status.conversion_frozen &&
            step.conversion_events.empty() &&
            std::fabs(step.ledger.conversion_number_removed) <= freeze_tolerance &&
            std::fabs(step.ledger.conversion_px_removed) <= freeze_tolerance &&
            std::fabs(step.ledger.conversion_energy_removed) <= freeze_tolerance;
        audit_status.h10_frozen = audit_status.h10_frozen &&
            std::fabs(step.tail_return.number) <= freeze_tolerance &&
            std::fabs(step.tail_return.px) <= freeze_tolerance &&
            std::fabs(step.tail_return.energy) <= freeze_tolerance;
        append_progress(output_dir, rank, scale, step_index, total_substeps,
                        time + dt, step_index / steps_per_interval,
                        "RUNNING");
        if (step_index % steps_per_interval != 0) continue;
        const int interval = step_index / steps_per_interval;
        const bool write_record = std::find(intervals.begin(), intervals.end(), interval) !=
            intervals.end();

        EMFields poisson_end = fields;
        field_solver.reconstruct_phi(poisson_end, rank, size);
        std::vector<double> rho_delta = subtract_charge(fields.rho, interval_start.rho);
        std::vector<double> bulk_charge_end = component_charge(
            bulk, integrator->tail_state(), beam, grid, 0);
        std::vector<double> tail_charge_end = component_charge(
            bulk, integrator->tail_state(), beam, grid, 1);
        std::vector<double> beam_charge_end = component_charge(
            bulk, integrator->tail_state(), beam, grid, 2);
        const std::vector<double> bulk_delta = subtract_charge(
            bulk_charge_end, bulk_charge_start);
        const std::vector<double> tail_delta = subtract_charge(
            tail_charge_end, tail_charge_start);
        const std::vector<double> beam_delta = subtract_charge(
            beam_charge_end, beam_charge_start);
        const std::vector<double> ion_delta(rho_delta.size(), 0.0);
        std::vector<double> boundary_delta = rho_delta;
        for (size_t i = 0; i < boundary_delta.size(); ++i) {
            boundary_delta[i] -= bulk_delta[i] + tail_delta[i] + beam_delta[i];
        }
        const OpenPoissonWorkIdentity poisson = field_solver.evaluate_work_identity(
            interval_start, poisson_end, rho_delta, rank, size);
        const OpenPoissonWorkIdentity poisson_bulk = field_solver.evaluate_work_identity(
            interval_start, poisson_end, bulk_delta, rank, size);
        const OpenPoissonWorkIdentity poisson_tail = field_solver.evaluate_work_identity(
            interval_start, poisson_end, tail_delta, rank, size);
        const OpenPoissonWorkIdentity poisson_beam = field_solver.evaluate_work_identity(
            interval_start, poisson_end, beam_delta, rank, size);
        const OpenPoissonWorkIdentity poisson_ion = field_solver.evaluate_work_identity(
            interval_start, poisson_end, ion_delta, rank, size);
        const OpenPoissonWorkIdentity poisson_boundary = field_solver.evaluate_work_identity(
            interval_start, poisson_end, boundary_delta, rank, size);
        const double poisson_tolerance = 8192.0 *
            std::numeric_limits<double>::epsilon() * poisson.scale;
        BackgroundTailStateSnapshot tail_snapshot;
        integrator->tail_state().export_accepted_state(tail_snapshot);
        const BeamPersistentState beam_snapshot = beam.export_persistent_state();
        IntervalRecord record = {};
        record.scale = scale;
        record.interval = interval;
        record.pair_signed = pair_signed_sum;
        record.pair_abs = pair_abs_sum;
        record.bulk_work = bulk_work_sum;
        record.tail_work = tail_work_sum;
        record.beam_work = beam_work_sum;
        record.field_energy_change = field_energy_change_sum;
        record.electrode_work = electrode_work_sum;
        record.bulk_norm = std::sqrt(global_sum(local_bulk_norm(bulk)));
        record.tail_norm = std::sqrt(global_sum(local_tail_norm(integrator->tail_state())));
        record.beam_norm = std::sqrt(global_sum(local_beam_norm(beam)));
        record.field_norm = std::sqrt(global_sum(local_field_norm(fields)));
        record.conversion_number = conversion_number_sum;
        record.conversion_px = conversion_px_sum;
        record.conversion_energy = conversion_energy_sum;
        record.conversion_event_count = conversion_event_count_sum;
        record.tail_count = global_sum_u64(
            static_cast<unsigned long long>(integrator->tail_state().particles.size()));
        record.tail_rng_counter = global_xor_u64(tail_snapshot.collision_rng_seed);
        record.beam_rng_counter = global_xor_u64(beam_snapshot.rng_state);
        record.rng_hash = global_xor_u64(fnv_bytes(tail_snapshot.collision_rng_seed,
                                                   &beam_snapshot.rng_state,
                                                   sizeof(beam_snapshot.rng_state)));
        record.accepted = step.accepted ? 1 : 0;
        record.finite = step.finite ? 1 : 0;
        const double poisson_component_error = poisson.potential_charge_work -
            (poisson_bulk.potential_charge_work + poisson_tail.potential_charge_work +
             poisson_beam.potential_charge_work + poisson_ion.potential_charge_work +
             poisson_boundary.potential_charge_work);
        record.poisson_identity_pass = poisson.finite &&
            std::fabs(poisson.residual) <= poisson_tolerance &&
            std::fabs(poisson_component_error) <= poisson_tolerance ? 1 : 0;
        record.poisson_residual = poisson.residual;
        record.poisson_bulk_work = poisson_bulk.potential_charge_work;
        record.poisson_tail_work = poisson_tail.potential_charge_work;
        record.poisson_beam_work = poisson_beam.potential_charge_work;
        record.poisson_ion_work = poisson_ion.potential_charge_work;
        record.poisson_boundary_source_work =
            poisson_boundary.potential_charge_work;
        record.poisson_component_reconstruction_error = poisson_component_error;
        if (write_record) output.push_back(record);
        interval_start = poisson_end;
        bulk_charge_start.swap(bulk_charge_end);
        tail_charge_start.swap(tail_charge_end);
        beam_charge_start.swap(beam_charge_end);
        pair_signed_sum = pair_abs_sum = bulk_work_sum = tail_work_sum =
            beam_work_sum = field_energy_change_sum = electrode_work_sum = 0.0;
        conversion_number_sum = conversion_px_sum = conversion_energy_sum = 0.0;
        conversion_event_count_sum = 0;
    }
    return true;
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
    bool valid = parse_args(argc, argv, options);
    SpatialGrid grid;
    if (valid) grid.init(rank, size);
    Species bulk;
    BeamPIC beam;
    EMFields fields;
    if (valid) {
        bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                  -Const::qe, Const::me, Param::dens, Param::temperature_e,
                  false, grid);
        beam.init(grid);
        fields.init(grid);
    }
    VpfpCheckpointControl control = {};
    VpfpCheckpointTailState tail;
    std::string checkpoint_error;
    if (valid) valid = read_vpfp_checkpoint(options.restart_dir, control, bulk,
                                             beam, fields, grid, &tail, rank,
                                             size, checkpoint_error);
    int global_valid = valid ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_valid, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (global_valid == 0) {
        if (rank == 0) std::cerr << "checkpoint read failed: " << checkpoint_error << "\n";
        MPI_Finalize();
        return 2;
    }

    // rho and Beam density are derived accepted-state data and are not part
    // of the checkpoint payload.  Production reconstructs them during the
    // next step; this audit needs the same reconstruction before defining the
    // left endpoint of its first Poisson work interval.  Otherwise interval 1
    // compares the restored E against an all-zero rho and reports a dt-
    // independent false Poisson residual.
    bulk.compute_moments();
    beam.deposit_density(grid, rank, size);
    const std::vector<double> empty_tail_density;
    const std::vector<double>& restored_tail_density =
        tail.present ? tail.tail.density : empty_tail_density;
    const std::vector<double> restored_ion_density(
        static_cast<size_t>(grid.nx_local), Param::dens);
    fields.set_charge_density(bulk, restored_tail_density, beam.density,
                              restored_ion_density);

    const std::uint64_t initial_hash = state_hash(bulk, beam, fields, tail.tail);
    const std::uint64_t initial_rng_hash = rng_state_hash(beam, tail.tail);
    const std::uint64_t initial_ledger_hash = checkpoint_ledger_hash(tail);
    const AuditBoundaryConfig boundary_config =
        read_checkpoint_boundary(options.restart_dir, rank);
    initialize_progress_file(options.output_dir, rank, control.time, control.step);
    std::vector<IntervalRecord> records;
    bool run_ok = true;
    SourceOwnershipAudit ownership;
    AuditRunStatus audit_status;
    for (size_t i = 0; i < options.scales.size(); ++i) {
        run_ok = run_scale(grid, control, bulk, beam, fields, tail,
                           boundary_config,
                           options.scales[i], options.intervals,
                           options.output_dir, rank, size,
                           records, ownership, audit_status,
                           initial_hash) && run_ok;
    }
    const std::uint64_t final_hash = state_hash(bulk, beam, fields, tail.tail);
    const std::uint64_t final_rng_hash = rng_state_hash(beam, tail.tail);
    const std::uint64_t final_ledger_hash = checkpoint_ledger_hash(tail);
    int state_unchanged = initial_hash == final_hash ? 1 : 0;
    int rng_unchanged = initial_rng_hash == final_rng_hash ? 1 : 0;
    int ledger_unchanged = initial_ledger_hash == final_ledger_hash ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &state_unchanged, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &rng_unchanged, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ledger_unchanged, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int global_run_ok = run_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_run_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int poisson_pass = 1;
    for (size_t i = 0; i < records.size(); ++i)
        poisson_pass = poisson_pass && records[i].poisson_identity_pass;
    MPI_Allreduce(MPI_IN_PLACE, &poisson_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int ownership_pass = audit_status.source_ownership_valid ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &ownership_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int same_initial_state = audit_status.same_initial_state ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &same_initial_state, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int same_physical_window = audit_status.same_physical_window ? 1 : 0;
    const int required_interval_count = static_cast<int>(options.intervals.size());
    const int expected_record_count = required_interval_count *
        static_cast<int>(options.scales.size());
    if (static_cast<int>(records.size()) != expected_record_count) {
        same_physical_window = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &same_physical_window, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int collision_frozen = audit_status.collision_frozen ? 1 : 0;
    int conversion_frozen = audit_status.conversion_frozen ? 1 : 0;
    int h10_frozen = audit_status.h10_frozen ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &collision_frozen, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &conversion_frozen, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &h10_frozen, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    double ownership_residual = 0.0;
    int ownership_count[7] = {};
    for (int source = 0; source < 7; ++source) {
        ownership_count[source] = ownership.count[source];
        ownership_residual = std::max(ownership_residual,
                                      ownership.residual[source]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &ownership_residual, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ownership_count, 7, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    if (rank == 0) {
        if (!ensure_directory(options.output_dir)) {
            std::cerr << "cannot create audit output directory: "
                      << options.output_dir << "\n";
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
        std::ofstream data((options.output_dir + "/local_defect_intervals.dat").c_str());
        data << std::setprecision(17)
             << "dt_scale interval field_particle_pair_residual_signed "
             << "field_particle_pair_residual_abs bulk_work tail_work beam_work "
             << "field_energy_change electrode_work bulk_state_norm tail_state_norm "
             << "beam_state_norm field_state_norm conversion_number conversion_px "
             << "conversion_energy conversion_event_count tail_particle_count "
             << "tail_rng_counter beam_rng_counter rng_hash accepted finite "
             << "poisson_identity_pass poisson_residual poisson_bulk_work "
             << "poisson_tail_work poisson_beam_work poisson_ion_work "
             << "poisson_boundary_source_work poisson_component_reconstruction_error\n";
        for (size_t i = 0; i < records.size(); ++i) {
            const IntervalRecord& r = records[i];
            data << r.scale << " " << r.interval << " " << r.pair_signed << " "
                 << r.pair_abs << " " << r.bulk_work << " " << r.tail_work << " "
                 << r.beam_work << " " << r.field_energy_change << " "
                 << r.electrode_work << " " << r.bulk_norm << " " << r.tail_norm
                 << " " << r.beam_norm << " " << r.field_norm << " "
                 << r.conversion_number << " " << r.conversion_px << " "
                 << r.conversion_energy << " " << r.conversion_event_count << " "
                 << r.tail_count << " " << r.tail_rng_counter << " "
                 << r.beam_rng_counter << " " << r.rng_hash << " "
                 << r.accepted << " " << r.finite << " "
                 << r.poisson_identity_pass << " " << r.poisson_residual << " "
                 << r.poisson_bulk_work << " " << r.poisson_tail_work << " "
                 << r.poisson_beam_work << " " << r.poisson_ion_work << " "
                 << r.poisson_boundary_source_work << " "
                 << r.poisson_component_reconstruction_error << "\n";
        }
        const bool unchanged = state_unchanged != 0;
        const bool comparable = global_run_ok != 0 && poisson_pass != 0 &&
            ownership_pass != 0 && unchanged && rng_unchanged != 0 &&
            ledger_unchanged != 0 && same_initial_state != 0 &&
            same_physical_window != 0 && collision_frozen != 0 &&
            conversion_frozen != 0 && h10_frozen != 0;
        append_progress(options.output_dir, rank, 0.0, 0, 0, control.time,
                        static_cast<int>(options.intervals.size()),
                        comparable ? "COMPLETE" : "FAILED_CONTRACT");
        std::ofstream summary((options.output_dir + "/audit_contract.result").c_str());
        summary << "checkpoint_derived_state_reconstructed=1\n"
                << "same_initial_state=" << same_initial_state << "\n"
                << "same_physical_window=" << same_physical_window << "\n"
                << "collision_frozen=" << collision_frozen << "\n"
                << "conversion_frozen=" << conversion_frozen << "\n"
                << "h10_frozen=" << h10_frozen << "\n"
                << "accepted_state_bitwise_equal_after_audit=" << (unchanged ? 1 : 0) << "\n"
                << "rng_bitwise_equal_after_audit=" << rng_unchanged << "\n"
                << "ledger_bitwise_equal_after_audit=" << ledger_unchanged << "\n"
                << "electrode_work_ownership_count=" << ownership_count[0] << "\n"
                << "background_boundary_energy_ownership_count=" << ownership_count[1] << "\n"
                << "beam_boundary_energy_ownership_count=" << ownership_count[2] << "\n"
                << "velocity_boundary_energy_ownership_count=" << ownership_count[3] << "\n"
                << "conversion_energy_ownership_count=" << ownership_count[4] << "\n"
                << "tail_return_energy_ownership_count=" << ownership_count[5] << "\n"
                << "collision_reservoir_ownership_count=" << ownership_count[6] << "\n"
                << "source_ownership_residual=" << ownership_residual << "\n"
                << "source_ownership_valid=" << ownership_pass << "\n"
                << "poisson_identity_pass=" << poisson_pass << "\n"
                << "status=" << (comparable ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return global_run_ok != 0 && poisson_pass != 0 && ownership_pass != 0 &&
        state_unchanged != 0 && rng_unchanged != 0 && ledger_unchanged != 0 &&
        same_initial_state != 0 &&
        same_physical_window != 0 && collision_frozen != 0 &&
        conversion_frozen != 0 && h10_frozen != 0 ? 0 : 1;
}
