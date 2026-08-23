#include "vpfp_integrator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
void capture_cell_integrated_number(const Species& species,
                                    std::vector<double>& number)
{
    const int nxl = species.sgrid->nx_local;
    const int ng = species.sgrid->nghost;
    number.assign(static_cast<size_t>(nxl), 0.0);
    if (species.cylindrical_mass_representation) {
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            double sum = 0.0;
            const int ixg = ng + ix;
            for (int iu = 0; iu < Param::Nv; ++iu) {
                for (int im = 0; im < Param::Nmu; ++im) {
                    sum += species.f[idx3(ixg, iu, im)];
                }
            }
            number[static_cast<size_t>(ix)] = sum;
        }
        return;
    }
    for (int ix = 0; ix < nxl; ++ix) {
        number[static_cast<size_t>(ix)] =
            species.number_density[static_cast<size_t>(ix)] * species.sgrid->dx;
    }
}

double global_sum(double local)
{
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

std::uint64_t global_sum_u64(std::uint64_t local)
{
    std::uint64_t global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    return global;
}

int global_sum_int(int local)
{
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

std::uint64_t global_max_u64(std::uint64_t local)
{
    std::uint64_t global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  MPI_COMM_WORLD);
    return global;
}

void global_max_doubles(const double local[7], double global[7])
{
    MPI_Allreduce(local, global, 7, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void record_remap_failure(VpfpStepResult& result,
                          const RemapDiagnostics& diag,
                          const char* stage)
{
    result.failure_stage = stage;
    result.failing_rank = diag.first_nonfinite_rank;
    result.failing_ix = diag.first_nonfinite_ix;
    result.failing_iupar = diag.first_nonfinite_iupar;
    result.failing_iuperp = diag.first_nonfinite_iuperp;
    result.input_min = diag.input_min_mass;
    result.input_max = diag.input_max_mass;
    result.output_min = diag.output_min_mass;
    result.output_max = diag.output_max_mass;
    result.first_nonfinite_value = diag.first_nonfinite_value;
}

void record_upar_audit(VpfpStepResult& result,
                       const RemapDiagnostics& diag)
{
    result.audit_valid = result.audit_valid && diag.audit_valid;
    if (!diag.audit_valid && result.audit_failure_code == 0)
        result.audit_failure_code = diag.audit_failure_code;
    result.audit_parcel_failure_reason = diag.audit_parcel_failure_reason;
    result.audit_parcel_failure_rank = diag.audit_parcel_failure_rank;
    result.audit_parcel_failure_ix = diag.audit_parcel_failure_ix;
    result.audit_parcel_failure_face = diag.audit_parcel_failure_face;
    result.audit_parcel_failure_iuperp = diag.audit_parcel_failure_iuperp;
    result.audit_parcel_failure_node_mass = diag.audit_parcel_failure_node_mass;
    result.audit_parcel_failure_target = diag.audit_parcel_failure_target;
    result.audit_parcel_failure_node_sum = diag.audit_parcel_failure_node_sum;
    result.audit_parcel_failure_scale = diag.audit_parcel_failure_scale;
}

enum ConversionFailureReason {
    ConversionFailureNone = 0,
    ConversionFailureAudit = 1,
    ConversionFailureNonfinite = 2,
    ConversionFailureSupportLimit = 3,
    ConversionFailureDuplicateId = 4,
    ConversionFailureFaceLedger = 5,
    ConversionFailureConservation = 6,
    ConversionFailureFidelity = 7,
    ConversionFailureIncomplete = 8,
    ConversionFailureParticleLimit = 9
};

const char* conversion_failure_reason_name(int reason)
{
    switch (reason) {
    case ConversionFailureAudit: return "audit-validation";
    case ConversionFailureNonfinite: return "nonfinite";
    case ConversionFailureSupportLimit: return "support-limit";
    case ConversionFailureDuplicateId: return "duplicate-id";
    case ConversionFailureFaceLedger: return "face-ledger";
    case ConversionFailureConservation: return "conservation";
    case ConversionFailureFidelity: return "fidelity";
    case ConversionFailureIncomplete: return "incomplete";
    case ConversionFailureParticleLimit: return "particle-limit";
    default: return "none";
    }
}
}

VpfpIntegrator::VpfpIntegrator(const OpenBackgroundBoundary& background_boundary,
                               OpenElectrostaticSolver& field_solver,
                               const CylindricalFokkerPlanckCollision& collision)
    : background_boundary_(background_boundary), field_solver_(field_solver),
      collision_(collision), partition_(NULL), converter_(NULL),
      background_tail_enabled_(false), fixed_state_field_particle_audit_mode_(false),
      tail_max_particles_(0),
      tail_max_number_fraction_(0.0),
      tail_collision_kernel_(TailCollisionKernel::None),
      hybrid_collision_active_(false),
      step_count_(0), beam_enabled_(true),
      tail_conversion_mode_(TailConversionMode::STATIC_CELL),
      collision_interface_zero_wall_validation_(false),
      collision_interface_exporting_absorbing_(false),
      tail_flux_quadrature_order_(4), tail_flux_max_supports_(7),
      tail_flux_max_created_particles_per_step_(0),
      tail_stage_trace_enabled_(false), stage_energy_audit_enabled_(false),
      field_particle_power_audit_enabled_(false),
      background_phase_space_mode_(BackgroundPhaseSpaceMode::STRANG_PPM),
      x_transport_velocity_mode_(XTransportVelocityMode::ANALYTIC_CELL_CENTER),
      energy_ledger_initialized_(false),
      accepted_field_energy_(0.0), accepted_background_kinetic_energy_(0.0),
      accepted_beam_kinetic_energy_(0.0), accepted_tail_kinetic_energy_(0.0),
      initialized_(false)
{
    field_particle_coupling_config_ = FieldParticleCouplingConfig();
}

VpfpIntegrator::VpfpIntegrator(
    const OpenBackgroundBoundary& background_boundary,
    OpenElectrostaticSolver& field_solver,
    const CylindricalFokkerPlanckCollision& collision,
    const HybridVelocityPartition& partition,
    BulkTailConverter& converter,
    bool background_tail_enabled)
    : background_boundary_(background_boundary), field_solver_(field_solver),
      collision_(collision), partition_(&partition), converter_(&converter),
      background_tail_enabled_(background_tail_enabled),
      fixed_state_field_particle_audit_mode_(false),
      tail_max_particles_(0), tail_max_number_fraction_(0.0),
      tail_collision_kernel_(TailCollisionKernel::None),
      hybrid_collision_active_(false),
      step_count_(0), beam_enabled_(true),
      tail_conversion_mode_(TailConversionMode::STATIC_CELL),
      collision_interface_zero_wall_validation_(false),
      collision_interface_exporting_absorbing_(false),
      tail_flux_quadrature_order_(4), tail_flux_max_supports_(7),
      tail_flux_max_created_particles_per_step_(0),
      tail_stage_trace_enabled_(false), stage_energy_audit_enabled_(false),
      field_particle_power_audit_enabled_(false),
      background_phase_space_mode_(BackgroundPhaseSpaceMode::STRANG_PPM),
      x_transport_velocity_mode_(XTransportVelocityMode::ANALYTIC_CELL_CENTER),
      energy_ledger_initialized_(false),
      accepted_field_energy_(0.0), accepted_background_kinetic_energy_(0.0),
      accepted_beam_kinetic_energy_(0.0), accepted_tail_kinetic_energy_(0.0),
      initialized_(false)
{
    field_particle_coupling_config_ = FieldParticleCouplingConfig();
}

void VpfpIntegrator::set_tail_stage_trace(bool enabled,
                                           const std::string& output_dir)
{
    tail_stage_trace_enabled_ = enabled;
    tail_stage_trace_output_dir_ = output_dir;
}

bool VpfpIntegrator::apply_tail_bulk_return(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    VpfpStepResult& result, int mpi_rank, int mpi_size)
{
    if (fixed_state_field_particle_audit_mode_ ||
        !background_tail_enabled_ || !tail_bulk_return_.config().enabled) {
        return true;
    }
    TailBulkReturnDiagnostics local;
    const bool local_ok = tail_bulk_return_.apply(
        bulk_trial, tail_trial, grid_, *partition_, step_count_, mpi_rank,
        mpi_size, local);
    int ok = local_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if (!ok) {
        result.failure_code = 10;
        result.failure_stage = "tail_bulk_return";
        result.conversion_ok = false;
        return false;
    }
    const unsigned long long local_counts[11] = {
        local.candidate_particles, local.resident_particles,
        local.attempted_groups, local.committed_groups,
        local.deferred_infeasible_groups, local.deferred_rank_boundary_groups,
        local.particles_removed, local.projection_invalid_input_cells,
        local.projection_insufficient_support_cells,
        local.projection_infeasible_invariant_cells,
        local.projection_representation_incompatible_cells };
    unsigned long long global_counts[11] = {};
    MPI_Allreduce(local_counts, global_counts, 11, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    // TailBulkReturn already performs the one required six-moment global
    // ledger reduction.  Re-summing those values here would multiply every
    // returned moment by the MPI size.  Only wall time is rank-local.
    double global_wall_seconds = 0.0;
    MPI_Allreduce(&local.wall_seconds, &global_wall_seconds, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    result.tail_return.candidate_particles = global_counts[0];
    result.tail_return.resident_particles = global_counts[1];
    result.tail_return.attempted_groups = global_counts[2];
    result.tail_return.committed_groups = global_counts[3];
    result.tail_return.deferred_infeasible_groups = global_counts[4];
    result.tail_return.deferred_rank_boundary_groups = global_counts[5];
    result.tail_return.particles_removed = global_counts[6];
    result.tail_return.projection_invalid_input_cells = global_counts[7];
    result.tail_return.projection_insufficient_support_cells = global_counts[8];
    result.tail_return.projection_infeasible_invariant_cells = global_counts[9];
    result.tail_return.projection_representation_incompatible_cells =
        global_counts[10];
    result.tail_return.number = local.number;
    result.tail_return.px = local.px;
    result.tail_return.jx_dx = local.jx_dx;
    result.tail_return.energy = local.energy;
    result.tail_return.pixx_dx = local.pixx_dx;
    result.tail_return.piperp_dx = local.piperp_dx;
    result.tail_return.wall_seconds = global_wall_seconds;
    result.tail_return.number_residual = local.number_residual;
    result.tail_return.px_residual = local.px_residual;
    result.tail_return.jx_residual = local.jx_residual;
    result.tail_return.energy_residual = local.energy_residual;
    result.tail_return.pixx_residual = local.pixx_residual;
    result.tail_return.piperp_residual = local.piperp_residual;
    result.tail_return.number_difference = local.number_difference;
    result.tail_return.px_difference = local.px_difference;
    result.tail_return.jx_difference = local.jx_difference;
    result.tail_return.energy_difference = local.energy_difference;
    result.tail_return.pixx_difference = local.pixx_difference;
    result.tail_return.piperp_difference = local.piperp_difference;
    result.tail_return.mpi_request_residual = local.mpi_request_residual;
    result.tail_return.returned_number_by_x = local.returned_number_by_x;
    result.tail_return.finite = local.finite;
    result.tail_return.committed = true;
    return true;
}

bool VpfpIntegrator::post_field_charge_invariance_transaction(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    double time, double dt, int mpi_rank, int mpi_size,
    VpfpPostFieldChargeInvarianceReport& report)
{
    report = VpfpPostFieldChargeInvarianceReport();
    const int nxl = grid_.nx_local;
    const double eps = std::numeric_limits<double>::epsilon();

    // Combined (bulk + tail) electron number per local cell BEFORE the
    // post-field operators.  Bulk uses the same integrated-number capture as
    // the Gate I pairing audit; the tail uses the freshly deposited PIC
    // density (deposit_density is idempotent w.r.t. the particle set).
    bulk_trial.compute_moments();
    tail_trial.deposit_density(grid_, mpi_rank, mpi_size);
    report.combined_number_before.assign(static_cast<size_t>(nxl), 0.0);
    {
        std::vector<double> bulk_number;
        capture_cell_integrated_number(bulk_trial, bulk_number);
        for (int ix = 0; ix < nxl; ++ix) {
            report.combined_number_before[static_cast<size_t>(ix)] =
                bulk_number[static_cast<size_t>(ix)] +
                tail_trial.density[static_cast<size_t>(ix)] * grid_.dx;
        }
    }

    // Same production post-field flags as advance_background /
    // advance_with_beam (section 7.11.17): collision-face export is only
    // observed under FLUX_INTERFACE + exporting-absorbing + a non-trivial
    // Chang-Cooper bulk integrator.
    const bool tail_on = background_tail_enabled_;
    const bool apply_upar_sink =
        tail_on && tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE;
    const bool observe_collision_flux =
        apply_upar_sink && collision_interface_exporting_absorbing_ &&
        !collision_.is_trivial() &&
        collision_.bulk_integrator() ==
            BulkCollisionIntegrator::CHANG_COOPER_FLUX;

    VpfpStepResult result = {};
    result.conversion_ok = true;
    // The collision-face conversion only runs when C2 exported interface
    // parcels.  When it does not run (no export, or conversion mode off), the
    // transaction is still charge-conserving, so conversion_ok defaults to
    // true exactly like the production result.conversion_ok and is cleared
    // only by an actual conversion failure below.
    report.conversion_ok = true;

    // Second Strang collision half C2 (production apply_collision_half with
    // collision_half == 1, exactly as in the accepted-step path).
    CollisionDiagnostics collision_diag = {};
    HybridCollisionDiagnostics second_hybrid;
    bool second_hybrid_valid = false;
    BulkTailFluxBatch second_collision_flux;
    report.collision_half_calls = 1;
    report.c2_ok = apply_collision_half(
        bulk_trial, tail_trial, time + dt, 0.5 * dt, 1, mpi_rank,
        collision_diag, second_hybrid, second_hybrid_valid,
        observe_collision_flux ? &second_collision_flux : NULL);
    report.tail_collision_stochastic_ran =
        second_hybrid_valid &&
        (second_hybrid.tail_tail_applied || second_hybrid.tail_bulk_applied);
    report.collision_face_export_number = collision_diag.interface_export_number;
    report.collision_reservoir_energy = collision_diag.reservoir_energy_change;
    stage_sources_.collision_reservoir_energy +=
        collision_diag.reservoir_energy_change;
    capture_dg_stage(VPFP_STAGE_COLLISION_HALF2, bulk_trial,
                     background_tail_enabled_ ? &tail_trial : NULL,
                     beam_enabled_ ? &beam_work_ : NULL, midpoint_fields_);

    // Collision-face conversion: the production convert_flux_batch on the C2
    // interface export.  interface_parcel_count is already global, so every
    // rank enters the transaction when any rank owns an export.
    if (report.c2_ok && observe_collision_flux &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        collision_diag.interface_parcel_count > 0) {
        second_collision_flux.recompute(partition_->min_conversion_energy);
        BulkTailConversionDiagnostics collision_conversion =
            converter_->convert_flux_batch(
                second_collision_flux, tail_trial, grid_, *partition_,
                static_cast<int>(step_count_),
                ConversionLocation::AFTER_COLLISION_HALF, mpi_rank,
                tail_flux_max_supports_);
        report.collision_face_conversion_calls = 1;
        report.conversion_ok = collision_conversion.complete &&
            collision_conversion.conservative &&
            collision_conversion.fidelity_ok && collision_conversion.finite;
        report.conversion_number_removed = collision_conversion.number_removed;
        report.conversion_energy_removed = collision_conversion.energy_removed;
        report.conversion_particles_created =
            collision_conversion.particles_created;
        // §15.13.4: carry per-cell C2 conversion diagnostics for Gate I.
        report.conversion_events.push_back(collision_conversion);
        int local_failure_reason = ConversionFailureNone;
        if (!report.conversion_ok) {
            if (!collision_conversion.finite)
                local_failure_reason = ConversionFailureNonfinite;
            else if (collision_conversion.support_limit_violation_count != 0)
                local_failure_reason = ConversionFailureSupportLimit;
            else if (collision_conversion.duplicate_id_count != 0)
                local_failure_reason = ConversionFailureDuplicateId;
            else if (collision_conversion.face_ledger_mismatch_count != 0)
                local_failure_reason = ConversionFailureFaceLedger;
            else if (!collision_conversion.conservative)
                local_failure_reason = ConversionFailureConservation;
            else if (!collision_conversion.fidelity_ok)
                local_failure_reason = ConversionFailureFidelity;
            else
                local_failure_reason = ConversionFailureIncomplete;
        }
        if (!synchronize_conversion_outcome(
                local_failure_reason, collision_conversion, mpi_rank, 8,
                result)) {
            report.conversion_ok = false;
        }
    }
    stage_sources_.conversion_energy += report.conversion_energy_removed;
    capture_dg_stage(VPFP_STAGE_CONVERSION_AFTER_COLLISION, bulk_trial,
                     background_tail_enabled_ ? &tail_trial : NULL,
                     beam_enabled_ ? &beam_work_ : NULL, midpoint_fields_);

    // H10 tail-to-bulk return (production apply_tail_bulk_return).
    report.tail_return_calls =
        (background_tail_enabled_ && tail_bulk_return_.config().enabled &&
         !fixed_state_field_particle_audit_mode_) ? 1 : 0;
    report.return_ok = apply_tail_bulk_return(
        bulk_trial, tail_trial, result, mpi_rank, mpi_size);
    report.tail_return = result.tail_return;
    capture_dg_stage(VPFP_STAGE_TAIL_BULK_RETURN, bulk_trial,
                     background_tail_enabled_ ? &tail_trial : NULL,
                     beam_enabled_ ? &beam_work_ : NULL, midpoint_fields_);

    // Combined electron number per local cell AFTER the post-field operators.
    bulk_trial.compute_moments();
    tail_trial.deposit_density(grid_, mpi_rank, mpi_size);
    report.combined_number_after.assign(static_cast<size_t>(nxl), 0.0);
    {
        std::vector<double> bulk_number;
        capture_cell_integrated_number(bulk_trial, bulk_number);
        for (int ix = 0; ix < nxl; ++ix) {
            report.combined_number_after[static_cast<size_t>(ix)] =
                bulk_number[static_cast<size_t>(ix)] +
                tail_trial.density[static_cast<size_t>(ix)] * grid_.dx;
        }
    }

    // Rank-local residuals, scales and totals.
    double cell_residual_linf = 0.0;
    double cell_scale = 1.0;
    double rho_scale = 1.0;
    double before_sum = 0.0;
    double after_sum = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        const double b = report.combined_number_before[static_cast<size_t>(ix)];
        const double a = report.combined_number_after[static_cast<size_t>(ix)];
        const double diff = std::fabs(a - b);
        cell_residual_linf = std::max(cell_residual_linf, diff);
        cell_scale = std::max(cell_scale, std::max(std::fabs(b), std::fabs(a)));
        rho_scale = std::max(rho_scale,
            std::max(std::fabs(b), std::fabs(a)) / grid_.dx);
        before_sum += b;
        after_sum += a;
    }
    report.rank_number_before = before_sum;
    report.rank_number_after = after_sum;
    const double rank_residual = std::fabs(after_sum - before_sum);
    const double rank_scale =
        std::max(1.0, std::max(std::fabs(before_sum), std::fabs(after_sum)));

    // One packed MAX reduction makes the report rank-consistent:
    //   [0] cell residual linf, [1] cell scale, [2] rank residual,
    //   [3] rank scale, [4] per-cell density residual linf, [5] rho scale.
    double local_metrics[6] = {
        cell_residual_linf, cell_scale, rank_residual, rank_scale,
        cell_residual_linf / grid_.dx, rho_scale
    };
    double global_metrics[6] = {};
    MPI_Allreduce(local_metrics, global_metrics, 6, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    report.cell_number_residual_linf = global_metrics[0];
    report.cell_number_tolerance = 4096.0 * eps * global_metrics[1];
    report.rank_number_residual_linf = global_metrics[2];
    report.rank_number_tolerance = 4096.0 * eps * global_metrics[3];
    report.rho_before_after_relative_linf =
        global_metrics[4] / std::max(1.0, global_metrics[5]);

    report.global_number_before = global_sum(report.rank_number_before);
    report.global_number_after = global_sum(report.rank_number_after);
    report.global_number_residual = std::fabs(
        report.global_number_after - report.global_number_before);
    report.global_number_tolerance =
        4096.0 * eps *
        std::max(1.0, std::max(std::fabs(report.global_number_before),
                               std::fabs(report.global_number_after)));

    return report.c2_ok && report.conversion_ok && report.return_ok;
}

void VpfpIntegrator::trace_tail_stage(int mpi_rank, double time,
                                      const char* stage) const
{
    if (!tail_stage_trace_enabled_) return;
    std::ostringstream path;
    path << tail_stage_trace_output_dir_ << "/tail_stage_trace_rank"
         << std::setfill('0') << std::setw(6) << mpi_rank << ".dat";
    std::ofstream out(path.str().c_str(), std::ios::app);
    if (!out) return;
    const double monotonic_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    out << std::setprecision(17) << step_count_ << " " << time << " "
        << mpi_rank << " " << stage << " "
        << tail_work_.particles.size() << " " << monotonic_s << "\n";
    out.flush();
}

bool VpfpIntegrator::append_collision_flux_parcels(
    const CollisionFaceFluxes& fluxes, BulkTailFluxBatch& batch,
    int collision_half, int mpi_rank) const
{
    if (partition_ == NULL || fluxes.nx_local != grid_.nx_local)
        return false;
    const CylindricalVelocityGrid& cgrid = state_np1_.cgrid;
    const int nq = tail_flux_quadrature_order_ >= 8 ? 8 : 4;
    static const double x4[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940523 };
    static const double w4[4] = {
         0.3478548451374539,  0.6521451548625461,
         0.6521451548625461,  0.3478548451374539 };
    static const double x8[8] = {
        -0.9602898564975363, -0.7966664774135913,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774135913,  0.9602898564975363 };
    static const double w8[8] = {
         0.1012285362903763,  0.2223810344533745,
         0.3137066458778873,  0.3626837833783620,
         0.3626837833783620,  0.3137066458778873,
         0.2223810344533745,  0.1012285362903763 };
    const double* gx = nq == 8 ? x8 : x4;
    const double* gw = nq == 8 ? w8 : w4;
    const int ng = grid_.nghost;
    const int stage = collision_half == 0 ? 2 : 3;
    // Collision face clipping can leave a positive roundoff residue on every
    // bulk-to-tail face.  Do not turn such residues into quadrature packets:
    // they are below one cell's representable particle budget and otherwise
    // create O(nx*Nv) useless parcels before the loader discards them.
    const double interface_roundoff_floor =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, Param::dens * grid_.dx);

    batch.apply_interface_sink = true;
    batch.finite = true;
    batch.nonnegative = true;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        const int ix_global = grid_.global_cell(ng + ix);
        for (size_t p = 0; p < partition_->upar_interface_faces.size(); ++p) {
            const BulkTailInterfaceFace& face =
                partition_->upar_interface_faces[p];
            const size_t q = fluxes.upar_index(ix, face.face_index,
                                               face.transverse_index);
            const double signed_transfer = fluxes.upar_flux[q] +
                fluxes.cross_upar_flux[q];
            const double amount = face.outward_sign > 0
                ? std::max(0.0, signed_transfer)
                : std::max(0.0, -signed_transfer);
            if (!(amount > interface_roundoff_floor)) continue;
            BulkTailFluxParcel parcel;
            parcel.ix_local = ix;
            parcel.ix_global = ix_global;
            parcel.direction = VelocityFaceDirection::U_PARALLEL;
            parcel.face_index = face.face_index;
            parcel.transverse_index = face.transverse_index;
            parcel.operator_stage = stage;
            parcel.face_number = amount;
            const double y0 = cgrid.uperp_faces[face.transverse_index] *
                cgrid.uperp_faces[face.transverse_index];
            const double y1 = cgrid.uperp_faces[face.transverse_index + 1] *
                cgrid.uperp_faces[face.transverse_index + 1];
            for (int g = 0; g < nq; ++g) {
                const double y = 0.5 * (y0 + y1) +
                    0.5 * (y1 - y0) * gx[g];
                const double weight = 0.5 * gw[g];
                const double uperp = std::sqrt(std::max(0.0, y));
                if (!bulk_tail_parcel_add_node(
                        parcel, cgrid.upar_faces[face.face_index], uperp,
                        amount * weight)) return false;
            }
            if (!parcel.finite_nonnegative()) return false;
            batch.parcels.push_back(parcel);
        }
        for (size_t p = 0; p < partition_->uperp_interface_faces.size(); ++p) {
            const BulkTailInterfaceFace& face =
                partition_->uperp_interface_faces[p];
            const size_t q = fluxes.uperp_index(ix, face.transverse_index,
                                                face.face_index);
            const double signed_transfer = fluxes.uperp_flux[q] +
                fluxes.cross_uperp_flux[q];
            const double amount = face.outward_sign > 0
                ? std::max(0.0, signed_transfer)
                : std::max(0.0, -signed_transfer);
            if (!(amount > interface_roundoff_floor)) continue;
            BulkTailFluxParcel parcel;
            parcel.ix_local = ix;
            parcel.ix_global = ix_global;
            parcel.direction = VelocityFaceDirection::U_PERP;
            parcel.face_index = face.face_index;
            parcel.transverse_index = face.transverse_index;
            parcel.operator_stage = stage;
            parcel.face_number = amount;
            const double u0 = cgrid.upar_faces[face.transverse_index];
            const double u1 = cgrid.upar_faces[face.transverse_index + 1];
            for (int g = 0; g < nq; ++g) {
                const double upar = 0.5 * (u0 + u1) +
                    0.5 * (u1 - u0) * gx[g];
                if (!bulk_tail_parcel_add_node(
                        parcel, upar,
                        cgrid.uperp_faces[face.face_index],
                        amount * 0.5 * gw[g])) return false;
            }
            if (!parcel.finite_nonnegative()) return false;
            batch.parcels.push_back(parcel);
        }
    }
    batch.recompute(partition_->min_conversion_energy);
    (void)mpi_rank;
    return batch.finite && batch.nonnegative && batch.duplicate_count == 0;
}

void VpfpIntegrator::append_flux_face_audit(
    const BulkTailFluxBatch& batch, VpfpStepResult& result) const
{
    result.flux_face_audit_count += batch.face_audit_count;
    result.flux_face_audit_face_abs_sum += batch.face_audit_face_abs_sum;
    result.flux_face_audit_parcel_abs_sum += batch.face_audit_parcel_abs_sum;
    result.flux_face_audit_abs_error_sum += batch.face_audit_abs_error_sum;
    if (batch.face_audit_max_valid &&
        (!result.flux_face_audit_max_valid ||
         batch.face_audit_max_relative > result.flux_face_audit_max_relative)) {
        result.flux_face_audit_max_valid = true;
        result.flux_face_audit_max_relative = batch.face_audit_max_relative;
        result.flux_face_audit_abs_at_max_relative =
            batch.face_audit_abs_at_max_relative;
        result.flux_face_audit_max = batch.face_audit_max;
    }
}

bool VpfpIntegrator::synchronize_conversion_outcome(
    int local_failure_reason,
    const BulkTailConversionDiagnostics& conversion,
    int mpi_rank,
    int failure_code,
    VpfpStepResult& result) const
{
    const int local_ok = local_failure_reason == ConversionFailureNone ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (global_ok != 0) return true;

    const int no_failure_rank = std::numeric_limits<int>::max();
    const int local_rank = local_ok != 0 ? no_failure_rank : mpi_rank;
    int first_failure_rank = no_failure_rank;
    MPI_Allreduce(&local_rank, &first_failure_rank, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    int local_detail[2] = { 0, 0 };
    unsigned long long local_counts[5] = { 0, 0, 0, 0, 0 };
    double local_residuals[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (mpi_rank == first_failure_rank) {
        local_detail[0] = local_failure_reason;
        local_detail[1] = conversion.finite ? 1 : 0;
        local_counts[0] = static_cast<unsigned long long>(
            conversion.particles_created);
        local_counts[1] = static_cast<unsigned long long>(
            conversion.compression_fallback_count);
        local_counts[2] = static_cast<unsigned long long>(
            conversion.support_limit_violation_count);
        local_counts[3] = static_cast<unsigned long long>(
            conversion.duplicate_id_count);
        local_counts[4] = static_cast<unsigned long long>(
            conversion.face_ledger_mismatch_count);
        local_residuals[0] = conversion.number_residual_rel;
        local_residuals[1] = conversion.px_residual_rel;
        local_residuals[2] = conversion.energy_residual_rel;
        local_residuals[3] = conversion.jx_residual_rel;
        local_residuals[4] = conversion.pixx_residual_rel;
        local_residuals[5] = conversion.piperp_residual_rel;
    }
    int global_detail[2] = { 0, 0 };
    unsigned long long global_counts[5] = { 0, 0, 0, 0, 0 };
    double global_residuals[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_detail, global_detail, 2, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_counts, global_counts, 5, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_residuals, global_residuals, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    result.conversion_ok = false;
    result.failure_code = failure_code;
    result.failure_stage = std::string("bulk_tail_conversion_") +
                           conversion_failure_reason_name(global_detail[0]);
    result.failing_rank = first_failure_rank;
    result.finite = global_detail[1] != 0;
    if (mpi_rank == 0) {
        std::fprintf(
            stderr,
            "[bulk-tail-conversion-fail] step=%lld failing_rank=%d "
            "reason=%s particles_created=%llu compression_fallbacks=%llu "
            "support_violations=%llu duplicate_ids=%llu "
            "face_ledger_mismatches=%llu residuals="
            "[%.6e,%.6e,%.6e,%.6e,%.6e,%.6e]\n",
            step_count_, first_failure_rank,
            conversion_failure_reason_name(global_detail[0]),
            global_counts[0], global_counts[1], global_counts[2],
            global_counts[3], global_counts[4], global_residuals[0],
            global_residuals[1], global_residuals[2], global_residuals[3],
            global_residuals[4], global_residuals[5]);
    }
    return false;
}

bool VpfpIntegrator::apply_upar_flux_conversion(
    BulkTailFluxBatch& exported_flux,
    const BulkTailFluxBatch& first_collision_flux,
    const CollisionDiagnostics& first_collision,
    bool flux_mode,
    VpfpStepResult& result,
    int mpi_rank,
    int failure_code)
{
    const std::chrono::steady_clock::time_point conversion_begin =
        std::chrono::steady_clock::now();
    BulkTailConversionDiagnostics conversion;
    int local_failure_reason = ConversionFailureNone;

    // Gate F deterministic transaction: preserve the actual Vlasov, field,
    // Tail and Beam maps while excluding representation changes.  Do not
    // fabricate a conversion ledger for the frozen operation.
    if (fixed_state_field_particle_audit_mode_) {
        conversion.complete = true;
        conversion.finite = true;
        conversion.conservative = true;
        conversion.fidelity_ok = true;
        result.conversion_ok = true;
        return true;
    }

    // FLUX_INTERFACE is a physical representation change, so an invalid
    // parcel ledger is a hard conversion failure there.  FLUX_AUDIT never
    // reaches this branch: it keeps the static-cell physical conversion and
    // leaves audit_valid as diagnostics only.
    if (flux_mode && tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !result.audit_valid) {
        result.failure_code = failure_code;
        result.failure_stage = "flux_interface_audit_validation";
        result.conversion_ok = false;
        result.finite = false;
        local_failure_reason = ConversionFailureAudit;
    } else if (flux_mode &&
               tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE) {
        // The Vlasov remap owns the final export batch.  Collision-half
        // parcels are appended only after that remap has finished, so both
        // sources enter one conversion transaction and one face ledger.
        if (!first_collision_flux.parcels.empty()) {
            exported_flux.parcels.insert(exported_flux.parcels.end(),
                                         first_collision_flux.parcels.begin(),
                                         first_collision_flux.parcels.end());
            exported_flux.recompute(partition_->min_conversion_energy);
        }
        append_flux_face_audit(exported_flux, result);
        result.flux_parcel_count += first_collision.interface_parcel_count;
        result.flux_export_number += first_collision.interface_export_number;
        result.flux_export_energy += first_collision.interface_export_energy;
        conversion = converter_->convert_flux_batch(
            exported_flux, tail_work_, grid_, *partition_,
            static_cast<int>(step_count_),
            ConversionLocation::AFTER_U_SUBSTEP, mpi_rank,
            tail_flux_max_supports_);
    } else {
        // static-cell remains the explicit A/B baseline.  flux-audit is
        // read-only: it records the production parcel batch but still uses
        // the original static extractor and must not mutate physical state.
        if (tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT &&
            flux_mode) {
            append_flux_face_audit(exported_flux, result);
        }
        conversion = converter_->extract_after_substep(
            state_u_full_, tail_work_, grid_, *partition_,
            static_cast<int>(step_count_),
            ConversionLocation::AFTER_U_SUBSTEP, mpi_rank);
        ++result.static_extractor_call_count;
    }

    result.wall_conversion_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - conversion_begin).count();
    const bool particle_limit_exceeded =
        flux_mode && tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        tail_flux_max_created_particles_per_step_ > 0 &&
        conversion.particles_created > static_cast<std::uint64_t>(
            tail_flux_max_created_particles_per_step_);
    result.conversion_ok = conversion.complete && conversion.conservative &&
        conversion.fidelity_ok && conversion.finite && !particle_limit_exceeded;

    result.ledger.conversion_number_removed = conversion.number_removed;
    result.ledger.conversion_px_removed = conversion.px_removed;
    result.ledger.conversion_energy_removed = conversion.energy_removed;
    result.ledger.conversion_particles_created = conversion.particles_created;
    result.flux_compression_fallback_count +=
        conversion.compression_fallback_count;
    result.flux_roundoff_discarded_number +=
        conversion.roundoff_discarded_number;
    result.flux_subcell_fallback_count += conversion.subcell_fallback_count;
    result.flux_support_limit_violation_count +=
        conversion.support_limit_violation_count;
    result.flux_duplicate_id_count += conversion.duplicate_id_count;
    result.flux_face_ledger_mismatch_count +=
        conversion.face_ledger_mismatch_count;
    result.flux_conversion_wall_seconds = result.wall_conversion_seconds;
    result.conversion_events.push_back(conversion);
    result.ledger.conversion_number_residual_rel =
        conversion.number_residual_rel;
    result.ledger.conversion_px_residual_rel = conversion.px_residual_rel;
    result.ledger.conversion_energy_residual_rel =
        conversion.energy_residual_rel;

    if (local_failure_reason == ConversionFailureNone &&
        !result.conversion_ok) {
        if (!conversion.finite) {
            local_failure_reason = ConversionFailureNonfinite;
        } else if (conversion.support_limit_violation_count != 0) {
            local_failure_reason = ConversionFailureSupportLimit;
        } else if (conversion.duplicate_id_count != 0) {
            local_failure_reason = ConversionFailureDuplicateId;
        } else if (conversion.face_ledger_mismatch_count != 0) {
            local_failure_reason = ConversionFailureFaceLedger;
        } else if (!conversion.conservative) {
            local_failure_reason = ConversionFailureConservation;
        } else if (!conversion.fidelity_ok) {
            local_failure_reason = ConversionFailureFidelity;
        } else if (!conversion.complete) {
            local_failure_reason = ConversionFailureIncomplete;
        } else {
            local_failure_reason = ConversionFailureParticleLimit;
        }
    }
    return synchronize_conversion_outcome(local_failure_reason, conversion,
                                          mpi_rank, failure_code, result);
}

void VpfpIntegrator::init(const SpatialGrid& grid)
{
    grid_ = grid;
    const std::string name = "background_electrons";
    const SpeciesType type = SpeciesType::BACKGROUND_ELECTRON;
    state_collision_trial_.init(name, type, -Const::qe, Const::me,
                                Param::dens, Param::temperature_e, false,
                                grid);
    state_x_half_.init(name, type, -Const::qe, Const::me, Param::dens,
                       Param::temperature_e, false, grid);
    state_u_full_.init(name, type, -Const::qe, Const::me, Param::dens,
                       Param::temperature_e, false, grid);
    state_np1_.init(name, type, -Const::qe, Const::me, Param::dens,
                    Param::temperature_e, false, grid);
    beam_work_.init(grid);
    tail_accepted_.init(grid);
    tail_work_.init(grid);
    midpoint_fields_.init(grid);
    final_fields_.init(grid);
    vlasov_.init(grid, state_x_half_, background_boundary_);
    vlasov_.set_x_transport_velocity_mode(x_transport_velocity_mode_);
    pairing_calculator_.init(grid);
    // JC1 (section 4.7 item 6): pre-allocate face/rho/residual arrays for
    // the discrete-gradient path.  These are sized once and reused.
    tail_field_trial_.init(grid);
    beam_field_trial_.init(grid);
    field_n_pairing_.init(grid);
    pairing_field_guess_.assign(grid.nx_local + 1, 0.0);
    pairing_field_map_.assign(grid.nx_local + 1, 0.0);
    pairing_field_residual_.assign(grid.nx_local + 1, 0.0);
    trial_force_fields_.init(grid);
    energy_ledger_initialized_ = false;
    initialized_ = true;
}

bool VpfpIntegrator::initialize_tail_from_bulk(
    Species& electrons, int mpi_rank, int mpi_size,
    BulkTailConversionDiagnostics& diagnostics)
{
    diagnostics = BulkTailConversionDiagnostics();
    if (!background_tail_enabled_ ||
        tail_conversion_mode_ != TailConversionMode::FLUX_INTERFACE) {
        diagnostics.complete = true;
        diagnostics.finite = true;
        diagnostics.conservative = true;
        diagnostics.fidelity_ok = true;
        return true;
    }
    if (!initialized_ || partition_ == NULL || converter_ == NULL)
        return false;
    diagnostics = converter_->convert_initial_tail_cells(
        electrons, tail_accepted_, grid_, *partition_, 0, mpi_rank,
        tail_flux_max_supports_);
    if (!diagnostics.complete || !diagnostics.finite ||
        !diagnostics.conservative || !diagnostics.fidelity_ok)
        return false;
    tail_accepted_.deposit_density(grid_, mpi_rank, mpi_size);
    return tail_accepted_.finite() && tail_accepted_.nonnegative_weights();
}

void VpfpIntegrator::set_beam_enabled(bool enabled)
{
    beam_enabled_ = enabled;
}

void VpfpIntegrator::set_tail_conversion_mode(
    TailConversionMode mode, int quadrature_order, size_t max_supports,
    int max_created_particles)
{
    tail_conversion_mode_ = mode;
    tail_flux_quadrature_order_ = quadrature_order >= 8 ? 8 : 4;
    tail_flux_max_supports_ = std::max<size_t>(1, max_supports);
    tail_flux_max_created_particles_per_step_ =
        std::max(0, max_created_particles);
}

void VpfpIntegrator::set_tail_limits(std::uint64_t max_particles,
                                     double max_number_fraction)
{
    tail_max_particles_ = max_particles;
    tail_max_number_fraction_ = max_number_fraction;
}

void VpfpIntegrator::set_tail_population_control(
    const TailPopulationController::Config& config)
{
    population_controller_.configure(config);
}

void VpfpIntegrator::set_tail_collision(TailCollisionKernel kernel,
                                        double coulomb_log,
                                        TailCollisionWeightMode weight_mode,
                                        int max_substeps,
                                        double max_particle_growth)
{
    tail_collision_kernel_ = kernel;
    hybrid_collision_active_ =
        background_tail_enabled_ && kernel != TailCollisionKernel::None;
    hybrid_collision_config_.max_particle_growth =
        std::max(0.0, max_particle_growth);
    hybrid_collision_config_.bulk_provider =
        const_cast<CollisionCoefficientProvider*>(&collision_.provider());
    hybrid_collision_config_.requested_kernel = kernel;
    hybrid_collision_config_.tail_tail_kernel = TailCollisionKernel::None;
    hybrid_collision_config_.tail_bulk_kernel = TailCollisionKernel::None;
    hybrid_collision_config_.weight_mode = weight_mode;
    hybrid_collision_config_.max_substeps = std::max(1, max_substeps);
    hybrid_collision_config_.pairs.bulk_bulk = !collision_.is_trivial();
    hybrid_collision_config_.pairs.bulk_tail = false;
    hybrid_collision_config_.pairs.tail_bulk = false;
    hybrid_collision_config_.pairs.tail_tail = false;
    if (kernel == TailCollisionKernel::CoulombLandauNanbuPerez) {
        hybrid_collision_config_.tail_tail_kernel =
            TailCollisionKernel::CoulombLandauNanbuPerez;
        hybrid_collision_config_.tail_bulk_kernel =
            TailCollisionKernel::KramersMoyalSDE;
        hybrid_collision_config_.pairs.tail_tail = true;
        hybrid_collision_config_.pairs.tail_bulk = true;
        hybrid_collision_config_.pairs.bulk_tail = true;
    } else if (kernel == TailCollisionKernel::KramersMoyalSDE) {
        hybrid_collision_config_.tail_bulk_kernel =
            TailCollisionKernel::KramersMoyalSDE;
        hybrid_collision_config_.pairs.tail_bulk = true;
        hybrid_collision_config_.pairs.bulk_tail = true;
    }
    hybrid_collision_config_.coulomb_log = coulomb_log;
    hybrid_collision_config_.rng_seed_base =
        static_cast<std::uint64_t>(coulomb_log * 1.0e6) ^ 0x5eed;
    // Section 17.7: only cells owned by the face-aligned bulk/tail
    // partition are excluded from the bulk collision operator.  The
    // threshold-intersecting buffer cells remain bulk-owned; using the old
    // center-cell conversion mask here would create an artificial wall
    // before the interface face.
    bulk_collision_mask_.clear();
    if (background_tail_enabled_ && partition_ != NULL &&
        partition_->bulk_owned_cell.size() ==
            static_cast<size_t>(Param::Nvmu)) {
        bulk_collision_mask_.assign(Param::Nvmu, 1);
        for (size_t slot = 0; slot < bulk_collision_mask_.size(); ++slot) {
            bulk_collision_mask_[slot] =
                partition_->bulk_owned_cell[slot] != 0 ? 1 : 0;
        }
    }
    hybrid_collision_config_.bulk_velocity_mask =
        bulk_collision_mask_.empty() ? NULL : &bulk_collision_mask_;
    hybrid_collision_config_.bulk_integrator = collision_.bulk_integrator();
    hybrid_collision_config_.bulk_partition = partition_;
}

void VpfpIntegrator::restore_tail_cumulative(
    const VpfpTailCumulativeLedger& ledger)
{
    tail_cumulative_ = ledger;
}

double VpfpIntegrator::tail_total_weight(const BackgroundTailPIC& tail) const
{
    double total = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        total += tail.particles[i].weight;
    }
    return total;
}

double VpfpIntegrator::tail_total_kinetic_energy(
    const BackgroundTailPIC& tail) const
{
    double total = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma =
            std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        total += p.weight * (gamma - 1.0) * Const::me * Const::c * Const::c;
    }
    return total;
}

bool VpfpIntegrator::finite_species(const Species& species) const
{
    for (size_t i = 0; i < species.f.size(); ++i) {
        if (!std::isfinite(species.f[i])) return false;
    }
    return true;
}

double VpfpIntegrator::field_energy(const EMFields& fields) const
{
    return global_sum(local_field_energy(fields));
}

double VpfpIntegrator::local_field_energy(const EMFields& fields) const
{
    double local = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        const double el = fields.Ex_face[static_cast<size_t>(ix)];
        const double er = fields.Ex_face[static_cast<size_t>(ix + 1)];
        local += Const::eps0 * grid_.dx *
                 (el * el + el * er + er * er) / 6.0;
    }
    return local;
}

void VpfpIntegrator::capture_stage_energy(
    VpfpStepResult& result, int stage_id, const Species& bulk,
    const BackgroundTailPIC* tail, const BeamPIC* beam,
    const EMFields& fields, const VpfpStageEnergyRecord& sources) const
{
    if (!stage_energy_audit_enabled_) return;
    result.stage_energy_audit_enabled = true;
    if (result.stage_energy_count >= VPFP_STAGE_ENERGY_RECORD_COUNT) {
        result.stage_energy_audit_valid = false;
        return;
    }

    VpfpStageEnergyRecord record = sources;
    record.stage_id = stage_id;
    // These are intentionally rank-local.  All eleven records are packed
    // into one collective at acceptance, avoiding per-stage synchronization.
    record.bulk_kinetic = bulk.total_kinetic_energy();
    record.tail_kinetic = tail ? tail_total_kinetic_energy(*tail) : 0.0;
    record.beam_kinetic = beam ? beam->total_kinetic_energy() : 0.0;
    record.field_energy = local_field_energy(fields);
    record.delta_bulk_kinetic = 0.0;
    record.delta_tail_kinetic = 0.0;
    record.delta_beam_kinetic = 0.0;
    record.delta_field_energy = 0.0;
    record.stage_balance = 0.0;
    result.stage_energy[result.stage_energy_count++] = record;
}

void VpfpIntegrator::finalize_stage_energy_audit(VpfpStepResult& result) const
{
    if (!stage_energy_audit_enabled_) return;

    enum { packed_values = 21 };
    double local[VPFP_STAGE_ENERGY_RECORD_COUNT * packed_values] = {};
    double global[VPFP_STAGE_ENERGY_RECORD_COUNT * packed_values] = {};
    const int local_count = result.stage_energy_count;
    int global_count = 0;
    int mpi_rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    if (global_count != mpi_size * VPFP_STAGE_ENERGY_RECORD_COUNT) {
        result.stage_energy_audit_valid = false;
        return;
    }

    for (int i = 0; i < VPFP_STAGE_ENERGY_RECORD_COUNT; ++i) {
        const VpfpStageEnergyRecord& r = result.stage_energy[i];
        const int o = i * packed_values;
        local[o + 0] = r.bulk_kinetic;
        local[o + 1] = r.tail_kinetic;
        local[o + 2] = r.beam_kinetic;
        local[o + 3] = r.field_energy;
        // ConservativePpmRemap and CylindricalFpCollision already return
        // globally reduced diagnostics on every rank. Contribute those
        // values once here; summing one identical copy per rank would inflate
        // the stage sources by mpi_size. Beam, tail-outflow and conversion
        // ledgers below remain rank-local and must still be summed.
        local[o + 4] = mpi_rank == 0 ? r.background_left_inflow_energy : 0.0;
        local[o + 5] = mpi_rank == 0 ? r.background_left_outflow_energy : 0.0;
        local[o + 6] = mpi_rank == 0 ? r.background_right_inflow_energy : 0.0;
        local[o + 7] = mpi_rank == 0 ? r.background_right_outflow_energy : 0.0;
        local[o + 8] = r.beam_injected_energy;
        local[o + 9] = r.beam_outflow_energy;
        local[o + 10] = r.tail_outflow_energy;
        local[o + 11] = mpi_rank == 0 ? r.collision_reservoir_energy : 0.0;
        local[o + 12] = r.conversion_energy;
        local[o + 13] = r.electrostatic_boundary_work;
        // tail_return_energy is already a globally reduced H10 ledger.
        local[o + 14] = 0.0;
        // Gate C discrete-work fields.  The bulk u_parallel quantities are
        // already global RemapDiagnostics; contribute one copy (rank 0 only)
        // so SUM preserves them.  Tail/Beam kick work stay rank-local.
        local[o + 15] = mpi_rank == 0 ? r.bulk_upar_face_work : 0.0;
        local[o + 16] = mpi_rank == 0 ? r.bulk_upar_velocity_boundary_work : 0.0;
        local[o + 17] = mpi_rank == 0 ? r.bulk_upar_interface_energy_removed : 0.0;
        local[o + 18] = mpi_rank == 0 ? r.bulk_upar_identity_residual : 0.0;
        local[o + 19] = r.tail_kick_work;
        local[o + 20] = r.beam_kick_work;
    }
    MPI_Allreduce(local, global, VPFP_STAGE_ENERGY_RECORD_COUNT * packed_values,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < VPFP_STAGE_ENERGY_RECORD_COUNT; ++i) {
        VpfpStageEnergyRecord& r = result.stage_energy[i];
        const int o = i * packed_values;
        r.bulk_kinetic = global[o + 0];
        r.tail_kinetic = global[o + 1];
        r.beam_kinetic = global[o + 2];
        r.field_energy = global[o + 3];
        r.background_left_inflow_energy = global[o + 4];
        r.background_left_outflow_energy = global[o + 5];
        r.background_right_inflow_energy = global[o + 6];
        r.background_right_outflow_energy = global[o + 7];
        r.beam_injected_energy = global[o + 8];
        r.beam_outflow_energy = global[o + 9];
        r.tail_outflow_energy = global[o + 10];
        r.collision_reservoir_energy = global[o + 11];
        r.conversion_energy = global[o + 12];
        r.electrostatic_boundary_work = global[o + 13];
        r.tail_return_energy = (i >= VPFP_STAGE_TAIL_BULK_RETURN)
            ? result.tail_return.energy : 0.0;
        r.bulk_upar_face_work = global[o + 15];
        r.bulk_upar_velocity_boundary_work = global[o + 16];
        r.bulk_upar_interface_energy_removed = global[o + 17];
        r.bulk_upar_identity_residual = global[o + 18];
        r.tail_kick_work = global[o + 19];
        r.beam_kick_work = global[o + 20];
        if (i == 0) continue;
        const VpfpStageEnergyRecord& p = result.stage_energy[i - 1];
        r.delta_bulk_kinetic = r.bulk_kinetic - p.bulk_kinetic;
        r.delta_tail_kinetic = r.tail_kinetic - p.tail_kinetic;
        r.delta_beam_kinetic = r.beam_kinetic - p.beam_kinetic;
        r.delta_field_energy = r.field_energy - p.field_energy;
        const double background_net =
            (r.background_left_inflow_energy - p.background_left_inflow_energy) -
            (r.background_left_outflow_energy - p.background_left_outflow_energy) +
            (r.background_right_inflow_energy - p.background_right_inflow_energy) -
            (r.background_right_outflow_energy - p.background_right_outflow_energy);
        const double beam_net = (r.beam_injected_energy - p.beam_injected_energy) -
            (r.beam_outflow_energy - p.beam_outflow_energy);
        const double tail_out = r.tail_outflow_energy - p.tail_outflow_energy;
        const double reservoir = r.collision_reservoir_energy -
            p.collision_reservoir_energy;
        const double boundary_work = r.electrostatic_boundary_work -
            p.electrostatic_boundary_work;
        r.stage_balance = r.delta_bulk_kinetic + r.delta_tail_kinetic +
            r.delta_beam_kinetic + r.delta_field_energy -
            (background_net + beam_net - tail_out - reservoir + boundary_work);
    }
    result.stage_energy_audit_valid = result.stage_energy_audit_valid &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].bulk_kinetic) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].field_energy) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].bulk_upar_face_work) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].bulk_upar_velocity_boundary_work) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].bulk_upar_interface_energy_removed) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].bulk_upar_identity_residual) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].tail_kick_work) &&
        std::isfinite(result.stage_energy[VPFP_STAGE_FINAL_POISSON].beam_kick_work);
}

void VpfpIntegrator::capture_dg_stage(
    int stage_id, const Species& bulk, const BackgroundTailPIC* tail,
    const BeamPIC* beam, const EMFields& fields)
{
    if (!stage_energy_audit_enabled_) return;
    capture_stage_energy(stage_energy_scratch_, stage_id, bulk, tail, beam,
                         fields, stage_sources_);
}

void VpfpIntegrator::finalize_field_particle_power_audit(
    VpfpStepResult& result,
    const Species& bulk_n, const Species& bulk_np1,
    const BeamPIC& beam_n, const BeamPIC& beam_np1,
    const BackgroundTailPIC* tail_n, const BackgroundTailPIC* tail_np1,
    const EMFields& field_n, const EMFields& field_np1,
    double bulk_work_global, double tail_work_local,
    double beam_work_local, double dt, int mpi_rank, int mpi_size)
{
    result.pairing_audit_enabled = false;
    result.pairing_audit = FieldParticlePowerAuditResult();
    if (!field_particle_power_audit_enabled_) return;
    result.pairing_audit_enabled = true;

    FieldParticlePowerAuditWorkspace& ws = pairing_workspace_;
    ws.enabled = true;
    const int nxl = grid_.nx_local;
    const int ng = grid_.nghost;

    // Accepted (n) and candidate (n+1) cell-integrated number (m^-2).
    capture_cell_integrated_number(bulk_n, ws.bulk_number_n);
    capture_cell_integrated_number(bulk_np1, ws.bulk_number_np1);
    if (tail_n != NULL && tail_np1 != NULL) {
        ws.tail_number_n.assign(static_cast<size_t>(nxl), 0.0);
        ws.tail_number_np1.assign(static_cast<size_t>(nxl), 0.0);
        for (int i = 0; i < nxl; ++i) {
            ws.tail_number_n[static_cast<size_t>(i)] =
                tail_n->density[static_cast<size_t>(i)] * grid_.dx;
            ws.tail_number_np1[static_cast<size_t>(i)] =
                tail_np1->density[static_cast<size_t>(i)] * grid_.dx;
        }
    } else {
        ws.tail_number_n.clear();
        ws.tail_number_np1.clear();
    }
    ws.beam_number_n.assign(static_cast<size_t>(nxl), 0.0);
    ws.beam_number_np1.assign(static_cast<size_t>(nxl), 0.0);
    for (int i = 0; i < nxl; ++i) {
        ws.beam_number_n[static_cast<size_t>(i)] =
            beam_n.density[static_cast<size_t>(i)] * grid_.dx;
        ws.beam_number_np1[static_cast<size_t>(i)] =
            beam_np1.density[static_cast<size_t>(i)] * grid_.dx;
    }

    // Field states and cell charge change (C/m^2).
    ws.field_n_ex_face = field_n.Ex_face;
    ws.field_mid_ex_face = midpoint_fields_.Ex_face;
    ws.field_np1_ex_face = field_np1.Ex_face;
    ws.cell_charge_change.assign(static_cast<size_t>(nxl), 0.0);
    for (int i = 0; i < nxl; ++i) {
        ws.cell_charge_change[static_cast<size_t>(i)] =
            (field_np1.rho[static_cast<size_t>(ng + i)] -
             field_n.rho[static_cast<size_t>(ng + i)]) * grid_.dx;
    }
    if (!field_solver_.build_potential_pairing_field(
            field_n, field_np1, ws.potential_pair_ex_face,
            mpi_rank, mpi_size)) {
        ws.potential_pair_ex_face.assign(
            static_cast<size_t>(nxl) + 1,
            std::numeric_limits<double>::quiet_NaN());
    }
    std::vector<double> rho_delta(field_n.rho.size(), 0.0);
    const size_t nrho = std::min(field_n.rho.size(), field_np1.rho.size());
    for (size_t i = 0; i < nrho; ++i) {
        rho_delta[i] = field_np1.rho[i] - field_n.rho[i];
    }
    const OpenPoissonWorkIdentity poisson_identity =
        field_solver_.evaluate_work_identity(
            field_n, field_np1, rho_delta, mpi_rank, mpi_size);
    ws.poisson_potential_charge_work = mpi_rank == 0
        ? poisson_identity.potential_charge_work : 0.0;
    ws.poisson_identity_residual = mpi_rank == 0
        ? poisson_identity.residual : 0.0;

    // Final trajectory current (A/m^2) for the charge-current pairing.
    if (tail_np1 != NULL) {
        ws.tail_trajectory.after_second_drift_current_face =
            tail_np1->current_face_x;
    } else {
        ws.tail_trajectory.after_second_drift_current_face.clear();
    }
    ws.beam_trajectory.after_second_drift_current_face =
        beam_np1.current_face_x;

    // The x-remap boundary fluxes are already in bulk_x1/bulk_x2.  Build the
    // remaining Bulk source from the accepted production stages on either
    // side of those remaps.  This records collision, u-remap, conversion and
    // return changes without reconstructing or modifying the x current.
    ws.bulk_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.tail_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.beam_source_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.bulk_conversion_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.bulk_return_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.bulk_other_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.tail_conversion_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.tail_return_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.tail_outflow_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    ws.beam_injection_number_cell.assign(static_cast<size_t>(nxl), 0.0);
    const bool bulk_stages_valid =
        ws.bulk_number_pre_x1.size() == static_cast<size_t>(nxl) &&
        ws.bulk_number_post_x1.size() == static_cast<size_t>(nxl) &&
        ws.bulk_number_pre_x2.size() == static_cast<size_t>(nxl) &&
        ws.bulk_number_post_x2.size() == static_cast<size_t>(nxl);
    if (bulk_stages_valid) {
        for (int i = 0; i < nxl; ++i) {
            const size_t k = static_cast<size_t>(i);
            ws.bulk_source_number_cell[k] =
                (ws.bulk_number_pre_x1[k] - ws.bulk_number_n[k]) +
                (ws.bulk_number_pre_x2[k] - ws.bulk_number_post_x1[k]) +
                (ws.bulk_number_np1[k] - ws.bulk_number_post_x2[k]);
        }
    }
    for (size_t event = 0; event < result.conversion_events.size(); ++event) {
        const BulkTailConversionDiagnostics& conversion =
            result.conversion_events[event];
        for (int i = 0; i < nxl; ++i) {
            if (static_cast<size_t>(i) <
                    conversion.removed_bulk_number_by_x.size()) {
                ws.bulk_conversion_number_cell[static_cast<size_t>(i)] -=
                    conversion.removed_bulk_number_by_x[static_cast<size_t>(i)];
            }
            if (static_cast<size_t>(i) <
                    conversion.created_tail_number_by_x.size()) {
                ws.tail_source_number_cell[static_cast<size_t>(i)] +=
                    conversion.created_tail_number_by_x[static_cast<size_t>(i)];
                ws.tail_conversion_number_cell[static_cast<size_t>(i)] +=
                    conversion.created_tail_number_by_x[static_cast<size_t>(i)];
            }
        }
    }
    const std::vector<double>& returned =
        result.tail_return.returned_number_by_x;
    for (int i = 0; i < nxl && static_cast<size_t>(i) < returned.size(); ++i) {
        ws.tail_source_number_cell[static_cast<size_t>(i)] -=
            returned[static_cast<size_t>(i)];
        ws.bulk_return_number_cell[static_cast<size_t>(i)] +=
            returned[static_cast<size_t>(i)];
        ws.tail_return_number_cell[static_cast<size_t>(i)] -=
            returned[static_cast<size_t>(i)];
    }
    for (int i = 0; i < nxl; ++i) {
        const size_t k = static_cast<size_t>(i);
        ws.bulk_other_number_cell[k] = ws.bulk_source_number_cell[k] -
            ws.bulk_conversion_number_cell[k] - ws.bulk_return_number_cell[k];
    }
    if (tail_np1 != NULL) {
        const std::vector<double>& removed =
            tail_np1->truncated_in_domain_number();
        for (int i = 0; i < nxl && static_cast<size_t>(i) < removed.size(); ++i) {
            ws.tail_source_number_cell[static_cast<size_t>(i)] -= removed[static_cast<size_t>(i)];
            ws.tail_outflow_number_cell[static_cast<size_t>(i)] -= removed[static_cast<size_t>(i)];
        }
    }
    const std::vector<double>& beam_source =
        beam_np1.boundary_source_density_delta();
    for (int i = 0; i < nxl && static_cast<size_t>(i) < beam_source.size(); ++i) {
        ws.beam_source_number_cell[static_cast<size_t>(i)] =
            beam_source[static_cast<size_t>(i)] * grid_.dx;
        ws.beam_injection_number_cell[static_cast<size_t>(i)] =
            beam_source[static_cast<size_t>(i)] * grid_.dx;
    }

    // Rank-local energy terms (section 4.6.1): the electrode work and the
    // already-global bulk work are contributed once by rank 0; the field
    // energy change, tail/beam work, conversion and boundary terms are
    // rank-local and are summed by the single reduction below.
    ws.field_energy_change =
        local_field_energy(field_np1) - local_field_energy(field_n);
    ws.electrode_work = (mpi_rank == 0)
        ? result.ledger.electrostatic_boundary_work : 0.0;
    ws.bulk_delta_ke = (mpi_rank == 0) ? bulk_work_global : 0.0;
    ws.tail_delta_ke = tail_work_local;
    ws.beam_delta_ke = beam_work_local;
    ws.bulk_work_ledger = (mpi_rank == 0) ? bulk_work_global : 0.0;
    ws.tail_work_ledger = tail_work_local;
    ws.beam_work_ledger = beam_work_local;
    // Representation changes are internal transfers.  Gate I therefore uses
    // only their created-minus-removed defect, never the transferred energy.
    double local_conversion_defect = 0.0;
    for (size_t event = 0; event < result.conversion_events.size(); ++event) {
        local_conversion_defect +=
            result.conversion_events[event].energy_created -
            result.conversion_events[event].energy_removed;
    }
    ws.conversion_energy_removed = local_conversion_defect;
    if (mpi_rank == 0)
        ws.conversion_energy_removed += result.tail_return.energy_difference;
    // R_fp is defined from force-only kinetic-energy changes.  Physical
    // particle energy carried through an open boundary is not part of that
    // quantity and must not be added as a residual.  Boundary charge sources
    // already enter the exact continuity/Poisson term above.  A future
    // independently measured boundary-ledger *defect* may be placed here;
    // the physical source itself may not.
    ws.boundary_energy_source = 0.0;

    // Rank-local calculator.  §15.13.4 step 3: the G* endpoint aliases need
    // the neighbouring rank's boundary cell force current.  finalize() is
    // rank-local (no MPI), so run it once to obtain j_force_cell, exchange
    // the two boundary cell values with the neighbours, then run it again so
    // the face dual region integrals use the exact G* endpoints.
    FieldParticlePowerAuditResult local =
        pairing_calculator_.finalize(ws, dt, -Const::qe);
    if (field_particle_power_audit_enabled_ && mpi_size > 1 &&
        local.j_force_cell.size() == static_cast<size_t>(grid_.nx_local)) {
        double send_left = local.j_force_cell[0];
        double send_right = local.j_force_cell[grid_.nx_local - 1];
        const int left = mpi_rank > 0 ? mpi_rank - 1 : MPI_PROC_NULL;
        const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : MPI_PROC_NULL;
        MPI_Sendrecv(&send_left, 1, MPI_DOUBLE, left, 7711,
                     &ws.gstar_right_neighbor_cell, 1, MPI_DOUBLE, right, 7711,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send_right, 1, MPI_DOUBLE, right, 7712,
                     &ws.gstar_left_neighbor_cell, 1, MPI_DOUBLE, left, 7712,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local = pairing_calculator_.finalize(ws, dt, -Const::qe);
    }

    // Single packed reduction: one SUM for the additive residuals/cell sums
    // and one MAX for the continuity linf and the all-finite flag.
    double sum_local[23];
    sum_local[0] = local.poisson_transport_residual;
    sum_local[1] = local.current_pair_residual;
    sum_local[2] = local.conversion_transfer_residual;
    sum_local[3] = local.boundary_residual;
    sum_local[4] = local.reconstructed_full_residual;
    sum_local[5] = local.full_residual;
    sum_local[6] = local.bulk_cell_work_sum;
    sum_local[7] = local.tail_cell_work_sum;
    sum_local[8] = local.beam_cell_work_sum;
    sum_local[9] = local.poisson_transport_residual_endpoint;
    sum_local[10] = local.poisson_transport_residual_midpoint;
    sum_local[11] = local.poisson_transport_residual_discrete_gradient;
    sum_local[12] = tail_work_local;
    sum_local[13] = beam_work_local;
    sum_local[14] = local.force_work_residual_bulk;
    sum_local[15] = local.force_work_residual_tail;
    sum_local[16] = local.force_work_residual_beam;
    sum_local[17] = local.poisson_identity_crosscheck;
    // §15.13.4 step 3: face dual audit region integrals (-dt<E_pair, dual>).
    sum_local[18] = local.dual_left5_integral;
    sum_local[19] = local.dual_core90_integral;
    sum_local[20] = local.dual_right5_integral;
    sum_local[21] = local.dual_in_domain_work;
    sum_local[22] = local.boundary_force_work;
    double sum_global[23] = {};
    MPI_Allreduce(sum_local, sum_global, 23, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const auto encode_bad = [&](int bad_index) {
        return bad_index >= 0
            ? -static_cast<double>(
                  mpi_rank * (grid_.nx_global + 1) + bad_index + 1)
            : -std::numeric_limits<double>::infinity();
    };
    double max_local[10] = {
        local.continuity_linf,
        local.continuity_linf_bulk,
        local.continuity_linf_tail,
        local.continuity_linf_beam,
        local.current_pair_linf,
        local.all_finite ? 0.0 : 1.0,
        local.continuity_scale,
        encode_bad(local.continuity_bulk.bad_global_index),
        encode_bad(local.continuity_tail.bad_global_index),
        encode_bad(local.continuity_beam.bad_global_index)
    };
    double max_global[10] = {};
    MPI_Allreduce(max_local, max_global, 10, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    FieldParticlePowerAuditResult& r = result.pairing_audit;
    r.poisson_transport_residual = sum_global[0];
    r.current_pair_residual = sum_global[1];
    r.conversion_transfer_residual = sum_global[2];
    r.boundary_residual = sum_global[3];
    r.reconstructed_full_residual = sum_global[4];
    r.full_residual = sum_global[5];
    r.bulk_cell_work_sum = sum_global[6];
    r.tail_cell_work_sum = sum_global[7];
    r.beam_cell_work_sum = sum_global[8];
    r.poisson_transport_residual_endpoint = sum_global[9];
    r.poisson_transport_residual_midpoint = sum_global[10];
    r.poisson_transport_residual_discrete_gradient = sum_global[11];
    const double tail_work_global = sum_global[12];
    const double beam_work_global = sum_global[13];
    r.continuity_linf = max_global[0];
    r.continuity_linf_bulk = max_global[1];
    r.continuity_linf_tail = max_global[2];
    r.continuity_linf_beam = max_global[3];
    r.current_pair_linf = max_global[4];
    r.all_finite = max_global[5] == 0.0;
    r.poisson_identity_crosscheck = sum_global[17];
    // §15.13.4 step 3: globalized face dual audit region integrals.
    r.dual_left5_integral = sum_global[18];
    r.dual_core90_integral = sum_global[19];
    r.dual_right5_integral = sum_global[20];
    r.dual_in_domain_work = sum_global[21];
    r.boundary_force_work = sum_global[22];
    r.dual_plus_boundary_work = r.dual_in_domain_work +
                                r.boundary_force_work;
    r.dual_reconstruction_error = r.dual_plus_boundary_work -
                                  r.current_pair_residual;
    r.dual_reconstruction_tolerance =
        FieldParticlePowerAudit::roundoff_tolerance(
            std::max(1.0, std::max(std::fabs(r.dual_plus_boundary_work),
                                   std::fabs(r.current_pair_residual))), 1.0);
    r.dual_reconstruction_pass =
        std::isfinite(r.dual_reconstruction_error) &&
        std::fabs(r.dual_reconstruction_error) <=
            r.dual_reconstruction_tolerance;
    r.all_finite = r.all_finite &&
        std::isfinite(r.dual_in_domain_work) &&
        std::isfinite(r.boundary_force_work) &&
        std::isfinite(r.dual_plus_boundary_work) &&
        std::isfinite(r.dual_reconstruction_error);
    // §15.13.4 step 3: the face dual audit profile needs the rank-local
    // per-face arrays, which are not part of the packed scalar reduction.
    // Copy them from the rank-local calculator result so the writer can emit
    // field_particle_power_dual_face_<step>_rank<rank>.dat on every rank.
    r.j_charge_face = local.j_charge_face;
    r.j_force_cell = local.j_force_cell;
    r.gstar_j_force_face = local.gstar_j_force_face;
    r.dual_face_residual = local.dual_face_residual;

    // Nonlinear quantities computed after the reduction.
    // Work-sum vs existing-ledger comparison (section 4.4 item 3).
    const double work_scale = std::max(
        1.0, std::max(std::fabs(bulk_work_global),
                      std::max(std::fabs(tail_work_global),
                               std::fabs(beam_work_global))));
    const double work_tol = FieldParticlePowerAudit::roundoff_tolerance(
        work_scale, work_scale);
    r.force_work_residual_bulk = sum_global[14];
    r.force_work_residual_tail = sum_global[15];
    r.force_work_residual_beam = sum_global[16];
    r.force_work_residual = r.force_work_residual_bulk +
                            r.force_work_residual_tail +
                            r.force_work_residual_beam;
    r.local_work_sum_matches_existing_ledger =
        std::fabs(r.force_work_residual_bulk) <= work_tol &&
        std::fabs(r.force_work_residual_tail) <= work_tol &&
        std::fabs(r.force_work_residual_beam) <= work_tol;

    // Reconstruct only after every signed term has been globalized.  Doing
    // this before the force-work reduction can hide a cell-ledger mismatch.
    r.reconstructed_full_residual =
        r.poisson_transport_residual + r.current_pair_residual +
        r.force_work_residual + r.conversion_transfer_residual +
        r.boundary_residual;
    r.reconstruction_mismatch =
        std::fabs(r.reconstructed_full_residual - r.full_residual);
    r.reconstruction_scale = std::max(
        1.0, std::max(std::fabs(r.reconstructed_full_residual),
                      std::fabs(r.full_residual)));
    r.roundoff_tolerance = FieldParticlePowerAudit::roundoff_tolerance(
        r.reconstruction_scale, r.reconstruction_scale);
    // Independent Gate-F Poisson identity crosscheck tolerance (section I6
    // fix).  The crosscheck residual is the roundoff of
    //   Delta U_E - W_electrode - W_potential_charge
    // over terms of magnitude poisson_identity_scale (~1e3-1e4 J/m^2), so it
    // must not be gated by the strongly-cancelled reconstruction tolerance.
    r.poisson_identity_scale = poisson_identity.scale;
    r.poisson_crosscheck_tolerance = std::max(
        4.0 * r.roundoff_tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() *
            r.poisson_identity_scale);
    r.full_residual_reconstruction_pass = r.all_finite &&
        std::isfinite(r.reconstruction_mismatch) &&
        r.reconstruction_mismatch <= r.roundoff_tolerance;

    // Continuity pass at the global scale.
    r.continuity_scale = max_global[6];
    r.continuity_pass = std::isfinite(r.continuity_linf) &&
        r.continuity_linf <= FieldParticlePowerAudit::roundoff_tolerance(
            r.continuity_scale, r.continuity_scale);

    // A-F candidate fractions (section 5.1 item 5).
    const double full_abs = std::fabs(r.full_residual);
    const double denom = full_abs > 0.0 ? full_abs : 1.0;
    r.transport_fraction = std::fabs(r.poisson_transport_residual) / denom;
    r.work_current_fraction = std::fabs(r.current_pair_residual) / denom;
    r.time_center_fraction =
        std::fabs(r.poisson_transport_residual_discrete_gradient -
                  r.poisson_transport_residual_midpoint) / denom;
    r.pic_fraction = std::fabs(r.force_work_residual) / denom;
    r.conversion_fraction = std::fabs(r.conversion_transfer_residual) / denom;
    r.boundary_fraction = std::fabs(r.boundary_residual) / denom;

    // Stable root-cause mask (nonzero candidate residuals).
    int mask = PAIRING_CAUSE_NONE;
    if (r.poisson_transport_residual != 0.0) mask |= PAIRING_CAUSE_TRANSPORT;
    if (r.current_pair_residual != 0.0) mask |= PAIRING_CAUSE_WORK;
    if (r.time_center_fraction != 0.0) mask |= PAIRING_CAUSE_TIME;
    if (r.force_work_residual != 0.0) mask |= PAIRING_CAUSE_PIC;
    if (r.conversion_transfer_residual != 0.0) mask |= PAIRING_CAUSE_CONVERSION;
    if (r.boundary_residual != 0.0) mask |= PAIRING_CAUSE_BOUNDARY;
    r.root_cause_mask = mask;
    r.validated_e_pair_candidate = EPAIR_DISCRETE_GRADIENT;
    // Decode the three per-component first-bad locations (section I2).  Each
    // component's deterministic global first-bad is the minimum (rank, index)
    // pair; the aggregate first_bad_rank/index is the minimum across the three.
    const auto decode_bad = [&](double value, int& bad_rank, int& bad_index) {
        bad_rank = -1;
        bad_index = -1;
        if (std::isfinite(value)) {
            const long long encoded =
                static_cast<long long>(-value) - 1;
            bad_rank = static_cast<int>(
                encoded / static_cast<long long>(grid_.nx_global + 1));
            bad_index = static_cast<int>(
                encoded % static_cast<long long>(grid_.nx_global + 1));
        }
    };
    r.continuity_bulk = local.continuity_bulk;
    r.continuity_tail = local.continuity_tail;
    r.continuity_beam = local.continuity_beam;
    decode_bad(max_global[7], r.continuity_bulk.first_bad_rank,
               r.continuity_bulk.bad_global_index);
    decode_bad(max_global[8], r.continuity_tail.first_bad_rank,
               r.continuity_tail.bad_global_index);
    decode_bad(max_global[9], r.continuity_beam.first_bad_rank,
               r.continuity_beam.bad_global_index);
    r.first_bad_rank = -1;
    r.first_bad_index = -1;
    const auto consider_global = [&](int bad_rank, int bad_index) {
        if (bad_rank < 0 || bad_index < 0) return;
        if (r.first_bad_rank < 0 ||
            bad_rank < r.first_bad_rank ||
            (bad_rank == r.first_bad_rank && bad_index < r.first_bad_index)) {
            r.first_bad_rank = bad_rank;
            r.first_bad_index = bad_index;
        }
    };
    consider_global(r.continuity_bulk.first_bad_rank,
                    r.continuity_bulk.bad_global_index);
    consider_global(r.continuity_tail.first_bad_rank,
                    r.continuity_tail.bad_global_index);
    consider_global(r.continuity_beam.first_bad_rank,
                    r.continuity_beam.bad_global_index);
    r.dt_s = dt;
    r.valid = r.all_finite;
}

void VpfpIntegrator::initialize_energy_ledger(
    const Species& electrons, const BeamPIC& beam,
    const EMFields& fields, bool tail_on)
{
    if (energy_ledger_initialized_) return;
    accepted_field_energy_ = field_energy(fields);
    accepted_background_kinetic_energy_ =
        global_sum(electrons.total_kinetic_energy());
    accepted_beam_kinetic_energy_ = beam_enabled_
        ? global_sum(beam.total_kinetic_energy()) : 0.0;
    accepted_tail_kinetic_energy_ = tail_on
        ? global_sum(tail_total_kinetic_energy(tail_accepted_)) : 0.0;
    energy_ledger_initialized_ = true;
}

void VpfpIntegrator::finalize_energy_ledger(VpfpStepResult& result)
{
    VpfpStepLedger& ledger = result.ledger;
    ledger.background_boundary_energy_net =
        ledger.background_left_inflow_energy -
        ledger.background_left_outflow_energy +
        ledger.background_right_inflow_energy -
        ledger.background_right_outflow_energy;
    ledger.beam_boundary_energy_net =
        ledger.beam_injected_energy - ledger.beam_outflow_energy;
    ledger.domain_energy_before =
        ledger.field_energy_before +
        ledger.background_kinetic_energy_before +
        ledger.beam_kinetic_energy_before +
        ledger.tail_kinetic_energy_before;
    ledger.domain_energy_after =
        ledger.field_energy + ledger.background_kinetic_energy +
        ledger.beam_kinetic_energy + ledger.tail_kinetic_energy_after;
    ledger.domain_energy_change =
        ledger.domain_energy_after - ledger.domain_energy_before;
    ledger.accounted_energy_source =
        ledger.background_boundary_energy_net +
        ledger.beam_boundary_energy_net - ledger.tail_outflow_energy -
        ledger.collision_reservoir_energy +
        ledger.electrostatic_boundary_work;
    ledger.energy_balance_residual =
        ledger.domain_energy_change - ledger.accounted_energy_source;
    ledger.energy_balance_relative = std::fabs(ledger.energy_balance_residual) /
        std::max(1.0, std::fabs(ledger.domain_energy_change) +
                      std::fabs(ledger.accounted_energy_source));

    accepted_field_energy_ = ledger.field_energy;
    accepted_background_kinetic_energy_ = ledger.background_kinetic_energy;
    accepted_beam_kinetic_energy_ = ledger.beam_kinetic_energy;
    accepted_tail_kinetic_energy_ = ledger.tail_kinetic_energy_after;
}

void VpfpIntegrator::swap_emfields(EMFields& a, EMFields& b)
{
    a.Ex.swap(b.Ex);
    a.Ex_face.swap(b.Ex_face);
    a.phi.swap(b.phi);
    a.rho.swap(b.rho);
    std::swap(a.nx_total, b.nx_total);
    std::swap(a.dx, b.dx);
    std::swap(a.last_gauss_residual_l1, b.last_gauss_residual_l1);
    std::swap(a.last_gauss_residual_linf, b.last_gauss_residual_linf);
}

bool VpfpIntegrator::apply_collision_half(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    double time, double half_dt, int collision_half, int mpi_rank,
    CollisionDiagnostics& bulk_diag,
    HybridCollisionDiagnostics& hybrid_diag, bool& hybrid_valid,
    BulkTailFluxBatch* collision_flux)
{
    hybrid_valid = false;
    bool local_ok = false;
    CollisionFaceFluxes face_flux;
    if (hybrid_collision_active_) {
        hybrid_collision_config_.dt = half_dt;
        hybrid_collision_config_.time = time;
        hybrid_collision_config_.accepted_step =
            static_cast<int>(step_count_);
        hybrid_collision_config_.collision_half = collision_half;
        hybrid_collision_config_.mpi_rank = mpi_rank;
        hybrid_collision_config_.bulk_face_fluxes = collision_flux != NULL
            ? &face_flux : NULL;
        hybrid_collision_config_.bulk_integrator = collision_.bulk_integrator();
        hybrid_collision_config_.bulk_partition = partition_;
        hybrid_diag = hybrid_collision_step_.advance(
            bulk_trial, tail_trial, grid_, hybrid_collision_config_);
        bulk_diag = hybrid_diag.bulk_diag;
        hybrid_valid = true;
        local_ok = hybrid_diag.success && bulk_diag.success;
        if (local_ok && collision_flux != NULL &&
            !append_collision_flux_parcels(face_flux, *collision_flux,
                                           collision_half, mpi_rank)) {
            local_ok = false;
            bulk_diag.success = false;
            bulk_diag.transaction_rollback_count += 1;
        }
    } else {
        if (collision_flux != NULL &&
            collision_.bulk_integrator() ==
                BulkCollisionIntegrator::CHANG_COOPER_FLUX) {
            bulk_diag = collision_.apply_with_flux(
                bulk_trial, grid_, time, half_dt,
                bulk_collision_mask_.empty() ? NULL : &bulk_collision_mask_,
                partition_, face_flux);
            if (bulk_diag.success &&
                !append_collision_flux_parcels(face_flux, *collision_flux,
                                               collision_half, mpi_rank)) {
                bulk_diag.success = false;
                bulk_diag.transaction_rollback_count += 1;
            }
        } else {
            bulk_diag = collision_.apply(
                bulk_trial, grid_, time, half_dt,
                bulk_collision_mask_.empty() ? NULL : &bulk_collision_mask_);
        }
        local_ok = bulk_diag.success;
    }

    // Export diagnostics are attached to the same collision transaction as
    // the face-flux batch.  Keep them global so the accepted-step ledger does
    // not depend on the MPI partition that happened to own an interface.
    if (collision_flux != NULL && local_ok) {
        double local_number = 0.0;
        double local_energy = 0.0;
        for (size_t p = 0; p < collision_flux->parcels.size(); ++p) {
            local_number += collision_flux->parcels[p].number;
            local_energy += collision_flux->parcels[p].kinetic_energy;
        }
        bulk_diag.interface_export_number = global_sum(local_number);
        bulk_diag.interface_export_energy = global_sum(local_energy);
        bulk_diag.interface_parcel_count = global_sum_u64(
            static_cast<std::uint64_t>(collision_flux->parcels.size()));
    }

    // Preserve the collisionless production fast path: the zero provider is
    // an exact no-op and cannot create rank-dependent failure.
    if (!hybrid_collision_active_ && collision_.is_trivial()) {
        return local_ok;
    }

    // Collision backends are rank-local.  No rank may return before every
    // rank has made the same accept/reject decision, otherwise a failed rank
    // enters MPI_Finalize while its peers continue into later collectives.
    int local_ok_int = local_ok ? 1 : 0;
    int global_ok_int = 0;
    MPI_Allreduce(&local_ok_int, &global_ok_int, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (global_ok_int != 0) return true;

    const int no_failure_rank = std::numeric_limits<int>::max();
    int local_failure_rank = local_ok ? no_failure_rank : mpi_rank;
    int first_failure_rank = no_failure_rank;
    MPI_Allreduce(&local_failure_rank, &first_failure_rank, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);

    int local_detail[4] = { 0, 0, 0, 0 };
    std::uint64_t local_counts[3] = { 0, 0, 0 };
    if (!local_ok && mpi_rank == first_failure_rank) {
        HybridCollisionFailureReason hybrid_reason =
            hybrid_diag.failure_reason;
        if (!hybrid_collision_active_) {
            hybrid_reason = HybridCollisionFailureReason::BulkBulk;
        }
        const TailCollisionDiagnostics* tail_diag = NULL;
        if (hybrid_reason == HybridCollisionFailureReason::TailTail) {
            tail_diag = &hybrid_diag.tail_tail_diag;
        } else if (hybrid_reason == HybridCollisionFailureReason::TailBulk) {
            tail_diag = &hybrid_diag.tail_bulk_diag;
        }
        local_detail[0] = static_cast<int>(hybrid_reason);
        local_detail[1] = tail_diag != NULL
            ? static_cast<int>(tail_diag->failure_reason) : 0;
        local_detail[2] = hybrid_diag.failure_cell;
        local_detail[3] = tail_diag != NULL
            ? tail_diag->collision_substeps : 0;
        if (tail_diag != NULL) {
            local_counts[0] = tail_diag->particle_count_before;
            local_counts[1] = tail_diag->particle_count_attempted;
            local_counts[2] = tail_diag->particle_count_limit;
        }
    }
    int global_detail[4] = {};
    std::uint64_t global_counts[3] = {};
    MPI_Allreduce(local_detail, global_detail, 4, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_counts, global_counts, 3,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        const HybridCollisionFailureReason hybrid_reason =
            static_cast<HybridCollisionFailureReason>(global_detail[0]);
        const TailCollisionFailureReason tail_reason =
            static_cast<TailCollisionFailureReason>(global_detail[1]);
        std::fprintf(
            stderr,
            "[hybrid-collision-fail] step=%lld half=%d failing_rank=%d "
            "hybrid_reason=%s tail_reason=%s local_cell=%d "
            "required_substeps=%d max_substeps=%d "
            "particles_before=%llu particles_attempted=%llu "
            "particles_limit=%llu\n",
            step_count_, collision_half, first_failure_rank,
            hybrid_collision_failure_reason_name(hybrid_reason),
            tail_collision_failure_reason_name(tail_reason),
            global_detail[2], global_detail[3],
            hybrid_collision_config_.max_substeps,
            static_cast<unsigned long long>(global_counts[0]),
            static_cast<unsigned long long>(global_counts[1]),
            static_cast<unsigned long long>(global_counts[2]));
    }
    bulk_diag.success = false;
    hybrid_diag.success = false;
    return false;
}

void VpfpIntegrator::finalize_collision_ledger(
    const CollisionDiagnostics& first_bulk,
    const CollisionDiagnostics& second_bulk,
    const HybridCollisionDiagnostics& first_hybrid,
    const HybridCollisionDiagnostics& second_hybrid,
    bool first_hybrid_valid, bool second_hybrid_valid,
    VpfpStepLedger& ledger) const
{
    ledger.collision_reservoir_energy =
        first_bulk.reservoir_energy_change +
        second_bulk.reservoir_energy_change;
    if (!first_hybrid_valid && !second_hybrid_valid) return;

    double local_sum[7] = {
        first_hybrid.bulk_bulk_energy_change +
            second_hybrid.bulk_bulk_energy_change,
        first_hybrid.tail_tail_px_change +
            second_hybrid.tail_tail_px_change,
        first_hybrid.tail_tail_energy_change +
            second_hybrid.tail_tail_energy_change,
        first_hybrid.tail_px_change + second_hybrid.tail_px_change,
        first_hybrid.tail_energy_change + second_hybrid.tail_energy_change,
        first_hybrid.bulk_reaction_px_change +
            second_hybrid.bulk_reaction_px_change,
        first_hybrid.bulk_reaction_energy_change +
            second_hybrid.bulk_reaction_energy_change
    };
    double global_values[7] = {};
    MPI_Allreduce(local_sum, global_values, 7, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    ledger.collision_bulk_bulk_energy_change = global_values[0];
    ledger.collision_tail_tail_px_change = global_values[1];
    ledger.collision_tail_tail_energy_change = global_values[2];
    ledger.collision_tail_bulk_px_change = global_values[3];
    ledger.collision_tail_bulk_energy_change = global_values[4];
    ledger.collision_bulk_reaction_px_change = global_values[5];
    ledger.collision_bulk_reaction_energy_change = global_values[6];

    double local_maximum[2] = {
        std::max(first_hybrid.reaction_px_balance,
                 second_hybrid.reaction_px_balance),
        std::max(first_hybrid.reaction_energy_balance,
                 second_hybrid.reaction_energy_balance)
    };
    double global_maximum[2] = {};
    MPI_Allreduce(local_maximum, global_maximum, 2, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    ledger.collision_reaction_px_residual = global_maximum[0];
    ledger.collision_reaction_energy_residual = global_maximum[1];

    int local_flags[4] = {
        (first_hybrid.bulk_bulk_applied ||
         second_hybrid.bulk_bulk_applied) ? 1 : 0,
        (first_hybrid.tail_tail_applied ||
         second_hybrid.tail_tail_applied) ? 1 : 0,
        (first_hybrid.tail_bulk_applied ||
         second_hybrid.tail_bulk_applied) ? 1 : 0,
        (first_hybrid.bulk_reaction_applied ||
         second_hybrid.bulk_reaction_applied) ? 1 : 0
    };
    int global_flags[4] = {};
    MPI_Allreduce(local_flags, global_flags, 4, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    ledger.collision_pair_bulk_bulk = global_flags[0];
    ledger.collision_pair_tail_tail = global_flags[1];
    ledger.collision_pair_tail_bulk = global_flags[2];
    ledger.collision_pair_bulk_reaction = global_flags[3];
}

void VpfpIntegrator::globalize_conversion_ledger(
    VpfpStepLedger& ledger) const
{
    double local_sums[3] = {
        ledger.conversion_number_removed,
        ledger.conversion_px_removed,
        ledger.conversion_energy_removed
    };
    double global_sums[3] = {};
    MPI_Allreduce(local_sums, global_sums, 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    ledger.conversion_number_removed = global_sums[0];
    ledger.conversion_px_removed = global_sums[1];
    ledger.conversion_energy_removed = global_sums[2];

    double local_residuals[3] = {
        ledger.conversion_number_residual_rel,
        ledger.conversion_px_residual_rel,
        ledger.conversion_energy_residual_rel
    };
    double global_residuals[3] = {};
    MPI_Allreduce(local_residuals, global_residuals, 3, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    ledger.conversion_number_residual_rel = global_residuals[0];
    ledger.conversion_px_residual_rel = global_residuals[1];
    ledger.conversion_energy_residual_rel = global_residuals[2];

    ledger.conversion_particles_created =
        global_sum_u64(ledger.conversion_particles_created);
}

void VpfpIntegrator::globalize_conversion_diagnostics(
    VpfpStepResult& result) const
{
    double local[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < result.conversion_events.size(); ++i) {
        const BulkTailConversionDiagnostics& event =
            result.conversion_events[i];
        local[0] = std::max(local[0], event.number_residual_rel);
        local[1] = std::max(local[1], event.px_residual_rel);
        local[2] = std::max(local[2], event.energy_residual_rel);
        local[3] = std::max(local[3], event.jx_residual_rel);
        local[4] = std::max(local[4], event.pixx_residual_rel);
        local[5] = std::max(local[5], event.piperp_residual_rel);
        local[6] = std::max(local[6], event.rho_l2_before_after);
    }
    // Diagnostics are written by rank 0, but conversion is strongly
    // load-imbalanced.  Report the slowest rank by piggybacking the wall time
    // on the existing max reduction; a rank-local value hid the step-970
    // compressor stall behind rank 0's ~10 us empty conversion.
    local[7] = result.wall_conversion_seconds;
    double global[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local, global, 8, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    result.conversion_number_residual = global[0];
    result.conversion_px_residual = global[1];
    result.conversion_energy_residual = global[2];
    result.conversion_jx_residual = global[3];
    result.conversion_pixx_residual = global[4];
    result.conversion_piperp_residual = global[5];
    result.conversion_rho_l2 = global[6];
    result.wall_conversion_seconds = global[7];
    result.flux_conversion_wall_seconds = global[7];
    unsigned long long local_contract[3] = {
        result.flux_support_limit_violation_count,
        result.flux_duplicate_id_count,
        result.flux_face_ledger_mismatch_count
    };
    unsigned long long global_contract[3] = { 0, 0, 0 };
    MPI_Allreduce(local_contract, global_contract, 3,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    result.flux_support_limit_violation_count = global_contract[0];
    result.flux_duplicate_id_count = global_contract[1];
    result.flux_face_ledger_mismatch_count = global_contract[2];

    unsigned long long local_face_count =
        static_cast<unsigned long long>(result.flux_face_audit_count);
    unsigned long long global_face_count = 0;
    MPI_Allreduce(&local_face_count, &global_face_count, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    double local_face_summary[3] = {
        result.flux_face_audit_face_abs_sum,
        result.flux_face_audit_parcel_abs_sum,
        result.flux_face_audit_abs_error_sum
    };
    double global_face_summary[3] = { 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_face_summary, global_face_summary, 3, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    result.flux_face_audit_count =
        static_cast<std::uint64_t>(global_face_count);
    result.flux_face_audit_face_abs_sum = global_face_summary[0];
    result.flux_face_audit_parcel_abs_sum = global_face_summary[1];
    result.flux_face_audit_abs_error_sum = global_face_summary[2];

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double local_max = result.flux_face_audit_max_valid
        ? result.flux_face_audit_max_relative : -1.0;
    double global_max = -1.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    int local_owner = result.flux_face_audit_max_valid &&
        local_max == global_max ? rank : size;
    int owner = size;
    MPI_Allreduce(&local_owner, &owner, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    result.flux_face_audit_max_valid = owner < size;
    if (result.flux_face_audit_max_valid) {
        int location[6] = { -1, -1, -1, -1, -1, 0 };
        double values[6] = { 0.0, 0.0, 0.0, 0.0,
                             std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN() };
        if (rank == owner) {
            location[0] = result.flux_face_audit_max.ix_global;
            location[1] = result.flux_face_audit_max.direction;
            location[2] = result.flux_face_audit_max.face_index;
            location[3] = result.flux_face_audit_max.transverse_index;
            location[4] = result.flux_face_audit_max.operator_stage;
            location[5] = result.flux_face_audit_max.node_failure_reason;
            values[0] = result.flux_face_audit_max.face_number;
            values[1] = result.flux_face_audit_max.parcel_number;
            values[2] = result.flux_face_audit_abs_at_max_relative;
            values[3] = global_max;
            values[4] = result.flux_face_audit_max.reconstructed_target;
            values[5] = result.flux_face_audit_max.node_sum;
        }
        MPI_Bcast(location, 6, MPI_INT, owner, MPI_COMM_WORLD);
        MPI_Bcast(values, 6, MPI_DOUBLE, owner, MPI_COMM_WORLD);
        result.flux_face_audit_max.ix_global = location[0];
        result.flux_face_audit_max.direction = location[1];
        result.flux_face_audit_max.face_index = location[2];
        result.flux_face_audit_max.transverse_index = location[3];
        result.flux_face_audit_max.operator_stage = location[4];
        result.flux_face_audit_max.node_failure_reason = location[5];
        result.flux_face_audit_max.face_number = values[0];
        result.flux_face_audit_max.parcel_number = values[1];
        result.flux_face_audit_abs_at_max_relative = values[2];
        result.flux_face_audit_max_relative = values[3];
        result.flux_face_audit_max.reconstructed_target = values[4];
        result.flux_face_audit_max.node_sum = values[5];
    }
}

VpfpStepResult VpfpIntegrator::advance_background(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density, double time, double dt,
    int mpi_rank, int mpi_size)
{
    VpfpStepResult result = {};
    result.stage_energy_audit_enabled = stage_energy_audit_enabled_;
    result.stage_energy_audit_valid = true;
    result.stage_energy_count = 0;
    result.tail_return_enabled =
        background_tail_enabled_ && tail_bulk_return_.config().enabled;
    result.split_used = false;
    result.finite = true;
    result.cfl_ok = dt > 0.0;
    result.gauss_ok = false;
    result.collision_ok = true;
    result.conversion_ok = true;
    result.tail_ok = true;
    result.failure_code = 0;
    result.failure_stage = "final_validation";
    // JC4 (section 7.7): legacy path has no discrete-gradient coupling.
    result.field_particle_coupling_enabled = false;
    result.post_field_charge_residual_linf = 0.0;
    result.failing_rank = -1;
    result.failing_ix = -1;
    result.failing_iupar = -1;
    result.failing_iuperp = -1;
    result.input_min = std::numeric_limits<double>::infinity();
    result.input_max = -std::numeric_limits<double>::infinity();
    result.output_min = std::numeric_limits<double>::infinity();
    result.output_max = -std::numeric_limits<double>::infinity();
    result.first_nonfinite_value = std::numeric_limits<double>::quiet_NaN();
    result.audit_valid = true;
    result.audit_failure_code = 0;
    result.audit_parcel_failure_reason =
        static_cast<int>(ParcelNodeFailureReason::None);
    result.audit_parcel_failure_rank = -1;
    result.audit_parcel_failure_ix = -1;
    result.audit_parcel_failure_face = -1;
    result.audit_parcel_failure_iuperp = -1;
    result.audit_parcel_failure_node_mass = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_target = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_node_sum = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_scale = std::numeric_limits<double>::quiet_NaN();
    result.audit_inplace_state_bitwise_equal = true;
    result.audit_inplace_rng_equal = true;
    result.audit_inplace_ledger_equal = true;
    result.tail_conversion_mode = static_cast<int>(tail_conversion_mode_);
    result.static_extractor_call_count = 0;
    result.flux_face_audit_count = 0;
    result.flux_face_audit_face_abs_sum = 0.0;
    result.flux_face_audit_parcel_abs_sum = 0.0;
    result.flux_face_audit_abs_error_sum = 0.0;
    result.flux_face_audit_max_relative = 0.0;
    result.flux_face_audit_abs_at_max_relative = 0.0;
    result.flux_face_audit_max_valid = false;
    if (!initialized_ || !result.cfl_ok) {
        result.failure_code = 1;
        return result;
    }
    (void)beam;
    const bool tail_on = background_tail_enabled_;
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        !collision_interface_zero_wall_validation_ &&
        !collision_interface_exporting_absorbing_) {
        // A non-trivial collision path must explicitly select either the
        // zero-wall validation contract or the production exporting-
        // absorbing contract.  Never infer an interface parcel from a cell
        // difference when neither contract was selected.
        result.failure_code = 12;
        result.collision_ok = false;
        result.conversion_ok = false;
        return result;
    }
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        collision_interface_exporting_absorbing_ &&
        collision_.bulk_integrator() !=
            BulkCollisionIntegrator::CHANG_COOPER_FLUX) {
        result.failure_code = 12;
        result.collision_ok = false;
        result.conversion_ok = false;
        return result;
    }
    ++step_count_;

    const std::chrono::steady_clock::time_point step_begin =
        std::chrono::steady_clock::now();
    result.ledger.background_number_before =
        global_sum(electrons.total_particle_number());
    result.ledger.beam_number_before = 0.0;
    if (tail_on) {
        // Trial tail: copy the accepted state, then snapshot its density and
        // reset the step ledgers (section 4.4).
        tail_work_ = tail_accepted_;
        tail_work_.begin_step(grid_, dt);
        result.ledger.tail_number_before =
            global_sum(tail_total_weight(tail_work_));
        result.ledger.tail_particle_count_before =
            global_sum_u64(tail_work_.particles.size());
        result.ledger.tail_kinetic_energy_before =
            global_sum(tail_total_kinetic_energy(tail_work_));
        result.ledger.combined_number_before =
            result.ledger.background_number_before +
            result.ledger.tail_number_before;
        trace_tail_stage(mpi_rank, time, "step_begin");
    }
    initialize_energy_ledger(electrons, beam, fields, tail_on);
    result.ledger.field_energy_before = accepted_field_energy_;
    result.ledger.background_kinetic_energy_before =
        accepted_background_kinetic_energy_;
    result.ledger.beam_kinetic_energy_before = accepted_beam_kinetic_energy_;
    if (!tail_on) result.ledger.tail_kinetic_energy_before = 0.0;
    result.conversion_ok = true;
    result.tail_ok = true;

    VpfpStageEnergyRecord stage_sources = {};
    const auto refresh_tail_stage_source = [&]() {
        stage_sources.tail_outflow_energy = tail_on
            ? tail_work_.outflow_ledger().left_kinetic_energy +
                  tail_work_.outflow_ledger().right_kinetic_energy
            : 0.0;
        stage_sources.conversion_energy = result.ledger.conversion_energy_removed;
    };
    const auto capture_stage = [&](int stage_id, const Species& bulk,
                                   const EMFields& stage_fields) {
        refresh_tail_stage_source();
        capture_stage_energy(result, stage_id, bulk,
                             tail_on ? &tail_work_ : NULL, NULL,
                             stage_fields, stage_sources);
    };
    capture_stage(VPFP_STAGE_ACCEPTED_N, electrons, fields);

    // Gate I (section 4.6.1): initialize the pairing workspace at the
    // accepted state (diagnostic level 2 only).  Level 0/1 leave it off so
    // no per-step scan or collective is added to production.
    pairing_workspace_.enabled = field_particle_power_audit_enabled_;
    if (field_particle_power_audit_enabled_) {
        pairing_workspace_.bulk_x1.init(grid_.nx_local);
        pairing_workspace_.bulk_x2.init(grid_.nx_local);
        pairing_workspace_.bulk_number_pre_x1.clear();
        pairing_workspace_.bulk_number_post_x1.clear();
        pairing_workspace_.bulk_number_pre_x2.clear();
        pairing_workspace_.bulk_number_post_x2.clear();
        pairing_workspace_.cell_work.init(grid_.nx_local);
        pairing_workspace_.tail_work_ledger = 0.0;
        pairing_workspace_.beam_work_ledger = 0.0;
    }
    XFaceTransportAudit* x1_audit = field_particle_power_audit_enabled_
        ? &pairing_workspace_.bulk_x1 : NULL;
    XFaceTransportAudit* x2_audit = field_particle_power_audit_enabled_
        ? &pairing_workspace_.bulk_x2 : NULL;
    std::vector<double>* bulk_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.bulk_delta_ke_cell : NULL;
    std::vector<double>* tail_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.tail_delta_ke_cell : NULL;
    double* tail_boundary_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.tail_delta_ke_boundary : NULL;
    double bulk_work_global = 0.0;
    double tail_work_local = 0.0;

    const bool observe_upar_flux = tail_on &&
        (tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT ||
         tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE);
    const bool apply_upar_sink = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE;
    const bool observe_collision_flux = apply_upar_sink &&
        collision_interface_exporting_absorbing_ &&
        !collision_.is_trivial() &&
        collision_.bulk_integrator() == BulkCollisionIntegrator::CHANG_COOPER_FLUX;
    // Keep the conversion helper's historical argument name local to the
    // u-parallel observation path only.  Collision-face export uses the
    // separate observe_collision_flux flag below.
    const bool flux_mode = observe_upar_flux;
    BulkTailFluxBatch exported_flux;
    BulkTailFluxBatch first_collision_flux;
    BulkTailFluxBatch second_collision_flux;
    // First Strang collision half-step C(dt/2).  With a real bulk collision
    // (stage H7) the collision acts on a trial copy of the accepted state,
    // so a rejected step never mutates the accepted bulk (section 14.7 item
    // 8).  collision=none stays a no-op without the per-step copy.
    const std::chrono::steady_clock::time_point first_collision_begin =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "collision_half1_begin");
    Species* collision_input = &electrons;
    CollisionDiagnostics first_collision = {};
    HybridCollisionDiagnostics first_hybrid;
    bool first_hybrid_valid = false;
    if (hybrid_collision_active_ || !collision_.is_trivial()) {
        state_collision_trial_ = electrons;
        collision_input = &state_collision_trial_;
    }
    const bool first_collision_ok = apply_collision_half(
        *collision_input, tail_work_, time, 0.5 * dt, 0, mpi_rank,
        first_collision, first_hybrid, first_hybrid_valid,
        observe_collision_flux ? &first_collision_flux : NULL);
    result.wall_collision_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - first_collision_begin).count();
    trace_tail_stage(mpi_rank, time, "collision_half1_end");
    if (!first_collision_ok) {
        result.failure_code = 5;
        result.collision_ok = false;
        return result;
    }
    stage_sources.collision_reservoir_energy +=
        first_collision.reservoir_energy_change;
    capture_stage(VPFP_STAGE_COLLISION_HALF1, *collision_input, fields);
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            *collision_input, pairing_workspace_.bulk_number_pre_x1);
    }
    if (observe_collision_flux) {
        result.collision_flux_half_export_number[0] =
            first_collision.interface_export_number;
        result.collision_flux_half_export_energy[0] =
            first_collision.interface_export_energy;
        result.collision_flux_half_implicit_residual_linf[0] =
            first_collision.implicit_flux_residual_linf;
        result.collision_flux_half_cross_pair_residual_linf[0] =
            first_collision.cross_flux_pair_residual_linf;
        result.collision_flux_half_inward_clipped_number[0] =
            first_collision.interface_inward_clipped_number;
        result.collision_flux_half_parcel_count[0] =
            first_collision.interface_parcel_count;
        result.collision_flux_half_rollback_count[0] =
            first_collision.transaction_rollback_count;
    }

    const std::chrono::steady_clock::time_point before_vlasov =
        std::chrono::steady_clock::now();
    VlasovStepDiagnostics vlasov_diag;
    exported_flux.apply_interface_sink = apply_upar_sink;
    trace_tail_stage(mpi_rank, time, "vlasov_x1_begin");
    if (!vlasov_.first_x_half(*collision_input, state_x_half_, time, 0.5 * dt,
                              vlasov_diag, x1_audit)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.x_first, "first_x");
        result.finite = false;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "vlasov_x1_end");
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_x_half_, pairing_workspace_.bulk_number_post_x1);
    }

    stage_sources.background_left_inflow_energy +=
        vlasov_diag.x_first.left_inflow_energy;
    stage_sources.background_left_outflow_energy +=
        vlasov_diag.x_first.left_outflow_energy;
    stage_sources.background_right_inflow_energy +=
        vlasov_diag.x_first.right_inflow_energy;
    stage_sources.background_right_outflow_energy +=
        vlasov_diag.x_first.right_outflow_energy;

    // Tail: first spatial half-drift to x^{n+1/2} with open-boundary
    // truncation and MPI migration (section 8.2 step 3), then the midpoint
    // density deposit (step 5).
    if (tail_on) {
        const std::chrono::steady_clock::time_point tail_push_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_drift1_begin");
        tail_work_.drift_half(grid_, 0.5 * dt, mpi_rank, mpi_size);
        if (field_particle_power_audit_enabled_) {
            tail_work_.finalize_trajectory_current(
                grid_, 0.5 * dt, mpi_rank, mpi_size);
            pairing_workspace_.tail_trajectory.after_first_drift_current_face =
                tail_work_.current_face_x;
        }
        result.wall_tail_push_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_push_begin).count();
        result.wall_tail_migrate_seconds +=
            tail_work_.last_migration_seconds();
        trace_tail_stage(mpi_rank, time, "tail_drift1_end");
        const std::chrono::steady_clock::time_point tail_dep_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_deposit_mid_begin");
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
        result.wall_tail_deposit_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_dep_begin).count();
        trace_tail_stage(mpi_rank, time, "tail_deposit_mid_end");
    }

    // Midpoint Poisson solve P[rho^{n+1/2}] on the cheap path (13.10):
    // no phi reconstruction, no L1/Linf/boundary audits.  The combined
    // background density is bulk (Eulerian) plus tail (PIC) (section 8.2
    // step 6 / section 9.1).
    capture_stage(VPFP_STAGE_X_HALF1, state_x_half_, fields);
    const std::vector<double> empty_beam_density;
    const std::vector<double> empty_tail_density;
    const std::vector<double>& tail_mid_density =
        tail_on ? tail_work_.density : empty_tail_density;
    midpoint_fields_.set_charge_density(state_x_half_, tail_mid_density,
                                        empty_beam_density, ion_density);
    OpenGaussSolveOptions midpoint_options;
    midpoint_options.reconstruct_phi = false;
    midpoint_options.compute_l1 = false;
    midpoint_options.compute_boundary_audit = false;
    const std::chrono::steady_clock::time_point before_field =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "field_mid_begin");
    field_solver_.solve(midpoint_fields_, mpi_rank, mpi_size, midpoint_options);
    result.wall_field_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before_field).count();
    result.mpi_collective_seconds +=
        field_solver_.diagnostics().mpi_collective_seconds;
    trace_tail_stage(mpi_rank, time, "field_mid_end");
    capture_stage(VPFP_STAGE_MIDPOINT_POISSON, state_x_half_, midpoint_fields_);

    // T_u(E^{n+1/2}, dt): non-uniform u_parallel conservative remap.
    const bool inplace_flux_audit = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT;
    trace_tail_stage(mpi_rank, time, "vlasov_u_begin");
    if (!vlasov_.u_full(state_x_half_, state_u_full_, midpoint_fields_,
                        time + 0.5 * dt, dt, vlasov_diag,
                        observe_upar_flux ? partition_ : NULL,
                        observe_upar_flux ? &exported_flux : NULL,
                        tail_flux_quadrature_order_, bulk_cell_work)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.u_full, "upar");
        record_upar_audit(result, vlasov_diag.u_full);
        result.finite = false;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "vlasov_u_end");
    if (inplace_flux_audit) {
        // u_full accepts its physical inputs by const reference and has no
        // access to the accepted tail, RNG or cumulative ledger.  Hashing
        // every full phase-space array before and after this read-only call
        // doubled memory traffic in flux-audit without strengthening this
        // ownership invariant.  The remap reports its own mutable-scratch
        // check below; checkpoint hashes remain the cross-run state check.
        result.audit_inplace_state_bitwise_equal =
            vlasov_diag.u_full.audit_physical_state_bitwise_equal;
        result.audit_inplace_rng_equal = true;
        result.audit_inplace_ledger_equal = true;
        if (!result.audit_inplace_state_bitwise_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 64;
        }
        if (!result.audit_inplace_rng_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 128;
        }
        if (!result.audit_inplace_ledger_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 256;
        }
    }
    record_upar_audit(result, vlasov_diag.u_full);
    // Gate C (sections 7.3/7.6): copy the bulk u_parallel discrete-work
    // values (already-global RemapDiagnostics) into the stage source record.
    // They are captured into every subsequent stage and contributed once by
    // rank 0 in finalize_stage_energy_audit.
    stage_sources.bulk_upar_face_work =
        vlasov_diag.u_full.upar_internal_face_energy_transfer;
    stage_sources.bulk_upar_velocity_boundary_work =
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy;
    stage_sources.bulk_upar_interface_energy_removed =
        vlasov_diag.u_full.upar_interface_energy_removed;
    stage_sources.bulk_upar_identity_residual =
        vlasov_diag.u_full.upar_discrete_energy_identity_residual;
    bulk_work_global =
        vlasov_diag.u_full.upar_internal_face_energy_transfer +
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy -
        vlasov_diag.u_full.upar_interface_energy_removed;
    pairing_workspace_.bulk_work_ledger = bulk_work_global;
    result.flux_parcel_count = vlasov_diag.u_full.interface_parcel_count;
    result.flux_node_count = vlasov_diag.u_full.interface_node_count;
    result.flux_duplicate_count = vlasov_diag.u_full.interface_duplicate_count;
    result.flux_face_export_number =
        vlasov_diag.u_full.interface_face_export_number;
    result.flux_parcel_number =
        vlasov_diag.u_full.interface_parcel_number;
    result.flux_export_number = vlasov_diag.u_full.interface_export_number;
    result.flux_export_energy = vlasov_diag.u_full.interface_export_energy;
    result.flux_below_threshold_number =
        vlasov_diag.u_full.interface_below_threshold_number;
    result.flux_roundoff_discarded_number =
        vlasov_diag.u_full.interface_roundoff_discarded_number;
    result.flux_quadrature_error_max =
        vlasov_diag.u_full.interface_quadrature_error_max;
    result.flux_tail_owned_bulk_residual =
        vlasov_diag.u_full.tail_owned_bulk_residual;
    result.flux_tail_owned_expected_transfer_number =
        vlasov_diag.u_full.tail_owned_expected_transfer_number;
    result.flux_tail_owned_roundoff_discarded_number =
        vlasov_diag.u_full.tail_owned_roundoff_discarded_number;
    result.flux_quadrature_order = tail_flux_quadrature_order_;
    result.flux_max_supports = static_cast<int>(tail_flux_max_supports_);
    result.flux_below_threshold_number =
        vlasov_diag.u_full.interface_below_threshold_number;
    result.flux_roundoff_discarded_number =
        vlasov_diag.u_full.interface_roundoff_discarded_number;
    if (apply_upar_sink) {
        const double tail_owned_floor = 4096.0 *
            std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::fabs(
                vlasov_diag.u_full.tail_owned_expected_transfer_number));
        if (!(vlasov_diag.u_full.tail_owned_bulk_residual <=
              tail_owned_floor)) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "[tail-owned-balance-fail] failure_code=6 step=%lld "
                             "unexplained=%.17g expected_transfer=%.17g "
                             "face_export=%.17g threshold=%.17g\n",
                             step_count_,
                             vlasov_diag.u_full.tail_owned_bulk_residual,
                             vlasov_diag.u_full.tail_owned_expected_transfer_number,
                             vlasov_diag.u_full.interface_face_export_number,
                             tail_owned_floor);
            }
            result.failure_code = 6;
            result.failure_stage = "upar_tail_interface";
            result.tail_ok = false;
            result.finite = false;
            return result;
        }
    }
    const double tail_fraction =
        vlasov_diag.u_full.tail_number_loss /
        std::max(1.0, vlasov_diag.u_full.number_before);
    if (tail_fraction > Param::umax_loss_abort_fraction) {
        // Velocity-domain under-resolution (sections 6.3/6.5/7.4): the tail
        // loss must stay below the production threshold, never be wrapped.
        result.failure_code = 7;
        result.finite = false;
        return result;
    }

    // Tail: kick with the same E^{n+1/2} at x^{n+1/2} (section 8.2 step 9),
    // then convert the post-T_u bulk mass above the threshold into tail
    // particles at the midpoint spatial positions (step 11).  The converter
    // only operates on the trial bulk/tail; a failed conversion leaves both
    // trial states untouched and rejects the step.
    if (!tail_on) {
        capture_stage(VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK, state_u_full_,
                      midpoint_fields_);
        capture_stage(VPFP_STAGE_CONVERSION_AFTER_FORCE, state_u_full_,
                      midpoint_fields_);
    }
    if (tail_on) {
        const std::chrono::steady_clock::time_point tail_kick_begin =
            std::chrono::steady_clock::now();
        double tail_kick_work_local = 0.0;
        tail_work_.kick(grid_, midpoint_fields_, dt, mpi_rank, mpi_size,
                        &tail_kick_work_local, tail_cell_work,
                        tail_boundary_work);
        tail_work_local = tail_kick_work_local;
        pairing_workspace_.tail_work_ledger = tail_kick_work_local;
        // Gate C (section 7.4/7.6): rank-local kick work, summed over ranks
        // in finalize_stage_energy_audit.
        stage_sources.tail_kick_work = tail_kick_work_local;
        result.wall_tail_push_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_kick_begin).count();
        capture_stage(VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK, state_u_full_,
                      midpoint_fields_);
        if (!apply_upar_flux_conversion(
                exported_flux, first_collision_flux, first_collision,
                flux_mode, result, mpi_rank, 8)) {
            return result;
        }
        capture_stage(VPFP_STAGE_CONVERSION_AFTER_FORCE, state_u_full_,
                      midpoint_fields_);
        // Tail: second spatial half-drift with MPI migration and open
        // deletion (step 12).  The newly converted particles take part in
        // this half-drift without any repeated kick.
        {
            const std::chrono::steady_clock::time_point tail_push2_begin =
                std::chrono::steady_clock::now();
            tail_work_.drift_half(grid_, 0.5 * dt, mpi_rank, mpi_size);
            result.wall_tail_push_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tail_push2_begin).count();
            result.wall_tail_migrate_seconds +=
                tail_work_.last_migration_seconds();
        }
        // Stage H5 (section 7.10): population control runs after the final
        // drift, so every particle is locally owned and in-domain and no
        // MPI-in-flight particle is merged.  Boundary cells are skipped
        // inside.  The final density deposit below then uses the
        // post-control representation, whose N/Xw moments (and therefore
        // the CIC density) are preserved to the compression tolerance.
        if (!fixed_state_field_particle_audit_mode_ &&
            population_controller_.enabled() &&
            population_controller_.active_step(
                static_cast<int>(step_count_))) {
            const TailPopulationController::Diagnostics diag =
                population_controller_.apply(
                    tail_work_, grid_, *partition_,
                    static_cast<int>(step_count_), mpi_rank);
            result.population_control_applied = diag.applied;
            result.population_control_local_groups =
                diag.groups_compressed + diag.groups_split;
            result.population_control_local_fallbacks =
                diag.compression_fallback_count;
            result.population_control_groups =
                global_sum_int(diag.groups_compressed +
                               diag.groups_split);
            result.population_control_fallbacks =
                global_sum_int(diag.compression_fallback_count);
            result.population_control_particles_before =
                global_sum_u64(diag.particles_before_local);
            result.population_control_particles_after =
                global_sum_u64(diag.particles_after_local);
            // The control report must be a global statement about the whole
            // run: groups/fallbacks are summed over ranks and the seven
            // moment residuals are worst-cased (max) over ranks, matching
            // the already-global particle counts.  Rank-local values would
            // hide compression that happens on other ranks (section 7.10
            // diagnostics).
            global_max_doubles(diag.max_residual,
                               result.population_control_max_residual);
        }
    }

    // Second T_x(dt/2) (section 8.2 step 13): operates on the post-conversion
    // bulk, so the remap ledger naturally includes the converted mass.
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_u_full_, pairing_workspace_.bulk_number_pre_x2);
    }
    if (!vlasov_.second_x_half(state_u_full_, state_np1_, time + 0.5 * dt,
                               0.5 * dt, vlasov_diag, x2_audit)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.x_second, "second_x");
        result.finite = false;
        return result;
    }
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_np1_, pairing_workspace_.bulk_number_post_x2);
    }
    stage_sources.background_left_inflow_energy +=
        vlasov_diag.x_second.left_inflow_energy;
    stage_sources.background_left_outflow_energy +=
        vlasov_diag.x_second.left_outflow_energy;
    stage_sources.background_right_inflow_energy +=
        vlasov_diag.x_second.right_inflow_energy;
    stage_sources.background_right_outflow_energy +=
        vlasov_diag.x_second.right_outflow_energy;
    capture_stage(VPFP_STAGE_X_HALF2, state_np1_, midpoint_fields_);
    result.wall_vlasov_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before_vlasov).count();

    // Second Strang collision half-step C(dt/2).
    const std::chrono::steady_clock::time_point second_collision_begin =
        std::chrono::steady_clock::now();
    CollisionDiagnostics collision_diag = {};
    HybridCollisionDiagnostics second_hybrid;
    bool second_hybrid_valid = false;
    const bool second_collision_ok = apply_collision_half(
        state_np1_, tail_work_, time + dt, 0.5 * dt, 1, mpi_rank,
        collision_diag, second_hybrid, second_hybrid_valid,
        observe_collision_flux ? &second_collision_flux : NULL);
    result.wall_collision_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - second_collision_begin).count();
    if (!second_collision_ok) {
        result.failure_code = 5;
        result.collision_ok = false;
        return result;
    }
    stage_sources.collision_reservoir_energy +=
        collision_diag.reservoir_energy_change;
    capture_stage(VPFP_STAGE_COLLISION_HALF2, state_np1_, midpoint_fields_);
    if (observe_collision_flux) {
        result.collision_flux_half_export_number[1] =
            collision_diag.interface_export_number;
        result.collision_flux_half_export_energy[1] =
            collision_diag.interface_export_energy;
        result.collision_flux_half_implicit_residual_linf[1] =
            collision_diag.implicit_flux_residual_linf;
        result.collision_flux_half_cross_pair_residual_linf[1] =
            collision_diag.cross_flux_pair_residual_linf;
        result.collision_flux_half_inward_clipped_number[1] =
            collision_diag.interface_inward_clipped_number;
        result.collision_flux_half_parcel_count[1] =
            collision_diag.interface_parcel_count;
        result.collision_flux_half_rollback_count[1] =
            collision_diag.transaction_rollback_count;
    }
    if (observe_collision_flux &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        collision_diag.interface_parcel_count > 0) {
        // interface_parcel_count is already global.  Therefore every rank
        // enters this block when any rank owns a second-half export; ranks
        // with an empty local batch still execute an empty transaction and
        // participate in the conversion outcome consensus.
        second_collision_flux.recompute(partition_->min_conversion_energy);
        append_flux_face_audit(second_collision_flux, result);
        const std::chrono::steady_clock::time_point conversion_begin =
            std::chrono::steady_clock::now();
        BulkTailConversionDiagnostics collision_conversion =
            converter_->convert_flux_batch(
                second_collision_flux, tail_work_, grid_, *partition_,
                static_cast<int>(step_count_),
                ConversionLocation::AFTER_COLLISION_HALF, mpi_rank,
                tail_flux_max_supports_);
        result.wall_conversion_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - conversion_begin).count();
        result.conversion_events.push_back(collision_conversion);
        // This second collision-half conversion is a separate accepted
        // event.  It must enter the same bulk/tail number and moment ledger
        // as the first event; otherwise the tail balance would silently omit
        // particles created after the second Strang half.
        result.ledger.conversion_number_removed +=
            collision_conversion.number_removed;
        result.ledger.conversion_px_removed +=
            collision_conversion.px_removed;
        result.ledger.conversion_energy_removed +=
            collision_conversion.energy_removed;
        result.ledger.conversion_particles_created +=
            collision_conversion.particles_created;
        result.flux_parcel_count += collision_diag.interface_parcel_count;
        result.flux_export_number += collision_diag.interface_export_number;
        result.flux_export_energy += collision_diag.interface_export_energy;
        result.flux_compression_fallback_count +=
            collision_conversion.compression_fallback_count;
        result.flux_roundoff_discarded_number +=
            collision_conversion.roundoff_discarded_number;
        result.flux_subcell_fallback_count +=
            collision_conversion.subcell_fallback_count;
        result.flux_support_limit_violation_count +=
            collision_conversion.support_limit_violation_count;
        result.flux_duplicate_id_count +=
            collision_conversion.duplicate_id_count;
        result.flux_face_ledger_mismatch_count +=
            collision_conversion.face_ledger_mismatch_count;
        result.ledger.conversion_number_residual_rel = std::max(
            result.ledger.conversion_number_residual_rel,
            collision_conversion.number_residual_rel);
        result.ledger.conversion_px_residual_rel = std::max(
            result.ledger.conversion_px_residual_rel,
            collision_conversion.px_residual_rel);
        result.ledger.conversion_energy_residual_rel = std::max(
            result.ledger.conversion_energy_residual_rel,
            collision_conversion.energy_residual_rel);
        const bool local_conversion_ok = collision_conversion.complete &&
            collision_conversion.conservative &&
            collision_conversion.fidelity_ok && collision_conversion.finite;
        int local_failure_reason = ConversionFailureNone;
        if (!local_conversion_ok) {
            if (!collision_conversion.finite)
                local_failure_reason = ConversionFailureNonfinite;
            else if (collision_conversion.support_limit_violation_count != 0)
                local_failure_reason = ConversionFailureSupportLimit;
            else if (collision_conversion.duplicate_id_count != 0)
                local_failure_reason = ConversionFailureDuplicateId;
            else if (collision_conversion.face_ledger_mismatch_count != 0)
                local_failure_reason = ConversionFailureFaceLedger;
            else if (!collision_conversion.conservative)
                local_failure_reason = ConversionFailureConservation;
            else if (!collision_conversion.fidelity_ok)
                local_failure_reason = ConversionFailureFidelity;
            else
                local_failure_reason = ConversionFailureIncomplete;
        }
        if (!synchronize_conversion_outcome(
                local_failure_reason, collision_conversion, mpi_rank, 8,
                result)) {
            return result;
        }
        result.conversion_ok = result.conversion_ok && local_conversion_ok;
    }
    if (observe_collision_flux) {
        // apply_collision_half() has already reduced parcel totals over MPI;
        // do not reduce these global diagnostics a second time.
        result.collision_flux_export_number =
            first_collision.interface_export_number +
            collision_diag.interface_export_number;
        result.collision_flux_export_energy =
            first_collision.interface_export_energy +
            collision_diag.interface_export_energy;
        result.collision_flux_implicit_residual_linf = std::max(
            first_collision.implicit_flux_residual_linf,
            collision_diag.implicit_flux_residual_linf);
        result.collision_flux_cross_pair_residual_linf = std::max(
            first_collision.cross_flux_pair_residual_linf,
            collision_diag.cross_flux_pair_residual_linf);
        result.collision_flux_inward_clipped_number =
            first_collision.interface_inward_clipped_number +
            collision_diag.interface_inward_clipped_number;
        result.collision_flux_parcel_count =
            first_collision.interface_parcel_count +
            collision_diag.interface_parcel_count;
        result.collision_flux_rollback_count =
            first_collision.transaction_rollback_count +
            collision_diag.transaction_rollback_count;
    }

    capture_stage(VPFP_STAGE_CONVERSION_AFTER_COLLISION, state_np1_,
                  midpoint_fields_);

    // H10 runs once after the second collision half and after any collision
    // interface export.  Its deletion changes only the end-of-step
    // representation; the trajectory current accumulated during drift stays
    // in tail_work_ and is finalized below.
    if (!apply_tail_bulk_return(state_np1_, tail_work_, result,
                                mpi_rank, mpi_size)) return result;
    capture_stage(VPFP_STAGE_TAIL_BULK_RETURN, state_np1_, midpoint_fields_);

    // Final accepted-state moments, density and Poisson solve with the full
    // Gauss audit (phi reconstruction stays off until an output step needs it).
    state_np1_.compute_moments();
    if (tail_on) {
        // Tail final density deposit and trajectory current finalisation
        // (section 9.2: the trajectory current is diagnostic-only).
        const std::chrono::steady_clock::time_point tail_dep2_begin =
            std::chrono::steady_clock::now();
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
        tail_work_.finalize_trajectory_current(grid_, dt, mpi_rank,
                                               mpi_size);
        result.wall_tail_deposit_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_dep2_begin).count();
    }
    const std::vector<double>& tail_final_density =
        tail_on ? tail_work_.density : empty_tail_density;
    final_fields_.set_charge_density(state_np1_, tail_final_density,
                                     empty_beam_density, ion_density);
    OpenGaussSolveOptions final_options;
    final_options.reconstruct_phi = false;
    const std::chrono::steady_clock::time_point final_field_begin =
        std::chrono::steady_clock::now();
    field_solver_.solve(final_fields_, mpi_rank, mpi_size, final_options);
    result.wall_field_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - final_field_begin).count();
    result.mpi_collective_seconds +=
        field_solver_.diagnostics().mpi_collective_seconds;

    // Hard validation (section 7.4): finite states/fields, closed Gauss
    // recursion, macroscopic ledger closure, roundoff-scale nonnegativity.
    result.finite = finite_species(state_np1_);
    for (size_t i = 0; i < final_fields_.Ex_face.size(); ++i) {
        result.finite = result.finite &&
                        std::isfinite(final_fields_.Ex_face[i]);
    }
    for (size_t i = 0; i < final_fields_.rho.size(); ++i) {
        result.finite = result.finite &&
                        std::isfinite(final_fields_.rho[i]);
    }
    int finite = result.finite ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &finite, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    result.finite = finite != 0;

    double local_rho_scale = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        local_rho_scale = std::max(
            local_rho_scale,
            std::fabs(final_fields_.rho[grid_.nghost + ix] / Const::eps0));
    }
    double rho_scale = 0.0;
    MPI_Allreduce(&local_rho_scale, &rho_scale, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double field_scale =
        std::max(std::fabs(final_fields_.Ex_face.front()),
                 std::fabs(final_fields_.Ex_face.back())) /
        (grid_.dx * static_cast<double>(grid_.nx_global));
    result.gauss_ok =
        std::isfinite(field_solver_.diagnostics().residual_linf) &&
        field_solver_.diagnostics().residual_linf <=
            1.0e-10 * std::max(1.0, rho_scale + field_scale);
    if (!result.finite || !result.gauss_ok) {
        result.failure_code = 3;
        record_remap_failure(result, vlasov_diag.x_second,
                             "final_validation");
        return result;
    }

    // Background global ledger: the remap diagnostics carry same-order global
    // sums, so
    //   after == before + x_in - x_out - u_tail_loss
    // must close to the deterministic-summation error scale.
    const double x_before = vlasov_diag.x_first.number_before;
    const double x_in = vlasov_diag.x_first.inflow_number +
                        vlasov_diag.x_second.inflow_number;
    const double x_out = vlasov_diag.x_first.outflow_number +
                         vlasov_diag.x_second.outflow_number;
    const double u_tail = vlasov_diag.u_full.tail_number_loss;
    // The conversion happens per rank on its local cells, while the remap
    // diagnostics and the tail/combined sums are global; the ledger must use
    // the MPI-summed conversion mass (section 11.1).
    const double conversion_removed =
        global_sum(result.ledger.conversion_number_removed);
    // The conversion removes bulk mass after T_u and before the second T_x,
    // so the remap ledger must account for it explicitly.
    const double expected_after =
        x_before + x_in - x_out - u_tail - conversion_removed +
        result.tail_return.number;
    // H10 return is applied after x_second has completed.  Consequently
    // x_second.number_after is the pre-return bulk number and cannot be used
    // to validate a ledger whose expected side includes tail_return.number.
    // Compare against the actual post-return trial state instead.
    const double actual_background_after =
        global_sum(state_np1_.total_particle_number());
    const double ledger_error =
        std::fabs(actual_background_after - expected_after) /
        std::max(1.0, x_before);
    result.ledger.remap_ledger_residual = ledger_error;
    result.ledger.background_left_flux =
        (vlasov_diag.x_first.left_inflow_number -
         vlasov_diag.x_first.left_outflow_number) +
        (vlasov_diag.x_second.left_inflow_number -
         vlasov_diag.x_second.left_outflow_number);
    result.ledger.background_right_flux =
        (vlasov_diag.x_first.right_inflow_number -
         vlasov_diag.x_first.right_outflow_number) +
        (vlasov_diag.x_second.right_inflow_number -
         vlasov_diag.x_second.right_outflow_number);
    result.ledger.background_left_inflow_energy =
        vlasov_diag.x_first.left_inflow_energy +
        vlasov_diag.x_second.left_inflow_energy;
    result.ledger.background_left_outflow_energy =
        vlasov_diag.x_first.left_outflow_energy +
        vlasov_diag.x_second.left_outflow_energy;
    result.ledger.background_right_inflow_energy =
        vlasov_diag.x_first.right_inflow_energy +
        vlasov_diag.x_second.right_inflow_energy;
    result.ledger.background_right_outflow_energy =
        vlasov_diag.x_first.right_outflow_energy +
        vlasov_diag.x_second.right_outflow_energy;
    result.ledger.background_tail_number_loss = u_tail;
    result.ledger.background_tail_energy_loss =
        vlasov_diag.u_full.tail_energy_loss;
    if (ledger_error > 1.0e-9) {
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "[ledger-fail] failure_code=6 step=%lld "
                         "remap_ledger_residual=%.17g threshold=1e-9 "
                         "x_before=%.17g x_in=%.17g x_out=%.17g "
                         "u_tail=%.17g conversion_removed=%.17g "
                         "tail_return=%.17g pre_return_after=%.17g "
                         "actual_after=%.17g expected_after=%.17g\n",
                         step_count_, ledger_error, x_before, x_in, x_out,
                         u_tail, conversion_removed,
                         result.tail_return.number,
                         vlasov_diag.x_second.number_after,
                         actual_background_after, expected_after);
        }
        result.failure_code = 6;
        result.failure_stage = "background_ledger_after_tail_return";
        return result;
    }

    double local_negative_mass = 0.0;
    double local_positive_mass = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double mass =
                    state_np1_.f[idx3(grid_.nghost + ix, iv, imu)];
                if (mass < 0.0) local_negative_mass -= mass;
                else local_positive_mass += mass;
            }
        }
    }
    double debt[2] = { local_negative_mass, local_positive_mass };
    double global_debt[2] = { 0.0, 0.0 };
    MPI_Allreduce(debt, global_debt, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (global_debt[0] > 1.0e-10 * std::max(1.0, global_debt[1])) {
        // The PPM positivity limiter only allows roundoff-scale debt; a
        // macroscopic deficit is a hard ledger failure, not a retry trigger.
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "[negmass-fail] failure_code=6 step=%lld "
                         "negative_mass=%.17g positive_mass=%.17g "
                         "ratio=%.17g threshold=1e-10\n",
                         step_count_, global_debt[0], global_debt[1],
                         global_debt[0] / std::max(1.0, global_debt[1]));
        }
        result.failure_code = 6;
        return result;
    }

    // Tail validation (section 11): finite nonnegative weights, no failed
    // migration, number balance against the conversion source and the open
    // outflows, and the section 11.5 resource gates.
    if (tail_on) {
        result.tail_ok =
            tail_work_.finite() && tail_work_.nonnegative_weights() &&
            !tail_work_.migration_failed();
        if (!result.tail_ok) {
            result.failure_code = 9;
            result.finite = false;
            return result;
        }
        const double tail_after_w =
            global_sum(tail_total_weight(tail_work_));
        const double tail_out =
            global_sum(tail_work_.outflow_ledger().left_number +
                       tail_work_.outflow_ledger().right_number);
        const double expected_tail_after =
            result.ledger.tail_number_before + conversion_removed - tail_out -
            result.tail_return.number;
        const double tail_balance_scale =
            std::max(1.0, result.ledger.tail_number_before +
                              conversion_removed + tail_out);
        result.ledger.tail_number_balance_error =
            std::fabs(tail_after_w - expected_tail_after) /
            tail_balance_scale;
        result.ledger.tail_outflow_number = tail_out;
        result.ledger.tail_outflow_energy =
            global_sum(tail_work_.outflow_ledger().left_kinetic_energy +
                       tail_work_.outflow_ledger().right_kinetic_energy);
        if (result.ledger.tail_number_balance_error > 1.0e-9) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "[tailbalance-fail] failure_code=6 step=%lld "
                             "tail_balance=%.17g threshold=1e-9 "
                             "tail_before=%.17g conversion=%.17g "
                             "tail_out=%.17g tail_after=%.17g "
                             "expected_after=%.17g scale=%.17g\n",
                             step_count_,
                             result.ledger.tail_number_balance_error,
                             result.ledger.tail_number_before,
                             conversion_removed, tail_out, tail_after_w,
                             expected_tail_after, tail_balance_scale);
            }
            result.failure_code = 6;
            result.failure_stage = "tail_ledger_after_tail_return";
            return result;
        }
        if (tail_max_particles_ > 0 &&
            tail_work_.particles.size() > tail_max_particles_) {
            result.failure_code = 9;
            result.finite = false;
            return result;
        }
        result.ledger.tail_number_after = tail_after_w;
        result.ledger.tail_kinetic_energy_after =
            global_sum(tail_total_kinetic_energy(tail_work_));
        // Compute the post-step combined number from the trial bulk before
        // the transactional swap (identical to the post-swap value).
        result.ledger.combined_number_after =
            actual_background_after + result.ledger.tail_number_after;
        if (tail_max_number_fraction_ > 0.0 &&
            result.ledger.tail_number_after >
                tail_max_number_fraction_ *
                    std::max(1.0, result.ledger.combined_number_after)) {
            result.failure_code = 9;
            result.finite = false;
            return result;
        }
    }

    result.ledger.electrostatic_boundary_work =
        field_solver_.boundary_energy_work(fields, final_fields_,
                                           mpi_rank, mpi_size);
    // boundary_energy_work is already global; contribute it on one rank so
    // the later packed stage-audit SUM preserves its value.
    if (mpi_rank == 0) {
        stage_sources.electrostatic_boundary_work =
            result.ledger.electrostatic_boundary_work;
    }
    capture_stage(VPFP_STAGE_FINAL_POISSON, state_np1_, final_fields_);

    // Transactional accept: swap state and fields, no full-array copy.
    electrons.swap_state(state_np1_);
    if (tail_on) {
        tail_accepted_.swap_state(tail_work_);
        result.ledger.tail_particle_count_after =
            global_sum_u64(tail_accepted_.particles.size());
    }
    swap_emfields(fields, final_fields_);
    result.ledger.background_number_after =
        global_sum(electrons.total_particle_number());
    result.ledger.beam_number_after = 0.0;
    result.ledger.beam_injected = 0.0;
    result.ledger.beam_outflow = 0.0;
    result.ledger.gauss_charge_residual =
        field_solver_.diagnostics().boundary_charge_residual;
    result.ledger.field_energy = field_energy(fields);
    result.ledger.background_kinetic_energy =
        global_sum(electrons.total_kinetic_energy());
    result.ledger.beam_kinetic_energy = 0.0;
    // Stage-H6 accepted cumulative ledgers and combined checksums (sections
    // 12.1/12.3).  Conversion/control counters are per-rank local; the
    // combined checksums are global.
    tail_cumulative_.conversion_number +=
        result.ledger.conversion_number_removed;
    tail_cumulative_.conversion_px += result.ledger.conversion_px_removed;
    tail_cumulative_.conversion_energy +=
        result.ledger.conversion_energy_removed;
    tail_cumulative_.return_number += result.tail_return.number;
    tail_cumulative_.return_px += result.tail_return.px;
    tail_cumulative_.return_jx_dx += result.tail_return.jx_dx;
    tail_cumulative_.return_energy += result.tail_return.energy;
    tail_cumulative_.return_pixx_dx += result.tail_return.pixx_dx;
    tail_cumulative_.return_piperp_dx += result.tail_return.piperp_dx;
    tail_cumulative_.return_particles_removed +=
        result.tail_return.particles_removed;
    tail_cumulative_.return_deferred_groups +=
        result.tail_return.deferred_infeasible_groups +
        result.tail_return.deferred_rank_boundary_groups;
    tail_cumulative_.particles_created +=
        result.ledger.conversion_particles_created;
    tail_cumulative_.control_groups += static_cast<std::uint64_t>(
        std::max(0, result.population_control_local_groups));
    tail_cumulative_.control_fallbacks += static_cast<std::uint64_t>(
        std::max(0, result.population_control_local_fallbacks));
    globalize_conversion_ledger(result.ledger);
    globalize_conversion_diagnostics(result);
    if (tail_on) {
        tail_cumulative_.outflow_number += static_cast<std::uint64_t>(
            tail_accepted_.outflow_ledger().left_number +
            tail_accepted_.outflow_ledger().right_number);
        combined_checksum_.number = result.ledger.combined_number_after;
        combined_checksum_.kinetic_energy =
            result.ledger.background_kinetic_energy +
            result.ledger.tail_kinetic_energy_after;
        result.tail_particles_local_max =
            global_max_u64(tail_accepted_.particles.size());
    } else {
        combined_checksum_.number = result.ledger.background_number_after;
        combined_checksum_.kinetic_energy =
            result.ledger.background_kinetic_energy;
        result.tail_particles_local_max = 0;
    }
    combined_checksum_.field_energy = result.ledger.field_energy;
    result.ledger.fct_energy_change = 0.0;
    finalize_collision_ledger(
        first_collision, collision_diag, first_hybrid, second_hybrid,
        first_hybrid_valid, second_hybrid_valid, result.ledger);
    finalize_energy_ledger(result);
    // Gate I is finalized only after all production ledgers have their final
    // global meaning.  The transactional swaps above leave the old accepted
    // state in the work objects, so both n and n+1 remain available without a
    // full-state copy.
    finalize_field_particle_power_audit(
        result, state_np1_, electrons, beam, beam,
        tail_on ? &tail_work_ : NULL, tail_on ? &tail_accepted_ : NULL,
        final_fields_, fields, bulk_work_global, tail_work_local,
        0.0, dt, mpi_rank, mpi_size);
    finalize_stage_energy_audit(result);
    result.wall_seconds_per_step = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - step_begin).count();
#ifndef _WIN32
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        result.max_rss_kib = usage.ru_maxrss;
    }
#endif
    // PPM limiter activity (section 6.4): mean over the three remap
    // substeps of the constant/linear reduced-order cell fractions.
    result.remap_constant_fraction =
        (vlasov_diag.x_first.constant_fraction +
         vlasov_diag.u_full.constant_fraction +
         vlasov_diag.x_second.constant_fraction) / 3.0;
    result.remap_linear_fraction =
        (vlasov_diag.x_first.linear_fraction +
         vlasov_diag.u_full.linear_fraction +
         vlasov_diag.x_second.linear_fraction) / 3.0;
    result.accepted = true;
    return result;
}

VpfpStepResult VpfpIntegrator::advance_joint_midpoint(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density, double time, double dt,
    int mpi_rank, int mpi_size)
{
    VpfpStepResult result = {};
    result.joint_midpoint_enabled = true;
    result.joint_midpoint_converged = false;
    result.joint_midpoint_iterations = 0;
    result.joint_midpoint_residual_linf = 0.0;
    result.joint_midpoint_poisson_residual_linf = 0.0;
    result.joint_midpoint_energy_residual = 0.0;
    result.joint_midpoint_pairing_field_built = false;
    result.joint_midpoint_pairing_face_left = 0.0;
    result.joint_midpoint_pairing_face_right = 0.0;
    result.joint_midpoint_force_current_first_cell = 0.0;
    result.joint_midpoint_force_current_last_cell = 0.0;
    result.joint_midpoint_naive_force_current_work = 0.0;
    result.joint_midpoint_seam_predicted_residual = 0.0;
    result.joint_midpoint_seam_prediction_error = 0.0;
    result.joint_midpoint_poisson_scalar_identity_residual = 0.0;
    result.joint_midpoint_continuity_charge_linf = 0.0;
    result.joint_midpoint_continuity_charge_l1 = 0.0;
    result.joint_midpoint_residual_charge_linf = 0.0;
    result.joint_midpoint_charge_projection_mismatch_linf = 0.0;
    result.joint_midpoint_u_boundary_charge_linf = 0.0;
    result.joint_midpoint_potential_weighted_continuity_defect = 0.0;
    result.joint_midpoint_poisson_current_predicted_residual = 0.0;
    result.joint_midpoint_poisson_current_prediction_error = 0.0;
    result.joint_midpoint_continuity_roundoff_bound = 0.0;
    result.joint_midpoint_prediction_roundoff_bound = 0.0;
    result.joint_midpoint_continuity_first_bad_global_ix = -1;
    result.joint_midpoint_density_assembly_mismatch_linf = 0.0;
    result.joint_midpoint_density_assembly_mismatch_l1 = 0.0;
    result.joint_midpoint_density_assembly_roundoff_bound = 0.0;
    result.joint_midpoint_mass_transport_charge_linf = 0.0;
    result.joint_midpoint_mass_transport_roundoff_bound = 0.0;
    result.joint_midpoint_transport_projection_mismatch_linf = 0.0;
    result.joint_midpoint_parent_charge_scale_max = 0.0;
    result.joint_midpoint_mass_delta_charge_linf = 0.0;
    result.joint_midpoint_density_assembly_first_bad_global_ix = -1;
    result.joint_midpoint_mass_transport_first_bad_global_ix = -1;
    result.joint_midpoint_potential_weighted_assembly_defect = 0.0;
    result.joint_midpoint_potential_weighted_transport_defect = 0.0;
    result.joint_midpoint_weighted_defect_reconstruction_error = 0.0;
    result.joint_midpoint_candidate_rho_incremental = 0.0;
    result.joint_midpoint_candidate_rho_absolute = 0.0;
    result.joint_midpoint_candidate_rho_form_difference = 0.0;
    result.joint_midpoint_candidate_rho_form_roundoff_bound = 0.0;
    result.joint_midpoint_poisson_identity_scale = 0.0;
    result.joint_midpoint_poisson_identity_roundoff_bound = 0.0;
    result.joint_midpoint_poisson_scalar_identity_pass = 0;
    result.joint_midpoint_poisson_identity_finite = 0;
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio = 0.0;
    result.joint_midpoint_poisson_identity_roundoff_bound_8192 = 0.0;
    result.joint_midpoint_poisson_identity_roundoff_bound_16384 = 0.0;
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio_8192 =
        0.0;
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio_16384 =
        0.0;
    result.joint_midpoint_poisson_scalar_identity_pass_8192 = 0;
    result.joint_midpoint_poisson_scalar_identity_pass_16384 = 0;
    result.joint_midpoint_poisson_term_abs_sum_energy_before = 0.0;
    result.joint_midpoint_poisson_term_abs_sum_energy_after = 0.0;
    result.joint_midpoint_poisson_term_abs_sum_potential_charge = 0.0;
    result.joint_midpoint_candidate_poisson_current_residual = 0.0;
    result.joint_midpoint_candidate_poisson_current_scale = 0.0;
    result.joint_midpoint_candidate_poisson_current_relative =
        std::numeric_limits<double>::infinity();
    result.joint_midpoint_pairing_tolerance = 0.0;
    result.joint_midpoint_candidate_poisson_scalar_residual = 0.0;
    result.joint_midpoint_candidate_weighted_continuity_defect = 0.0;
    result.joint_midpoint_phase_converged = 0;
    result.joint_midpoint_poisson_converged = 0;
    result.joint_midpoint_pairing_converged = 0;
    result.joint_midpoint_recent_pairing_relative_count = 0;
    for (double& value : result.joint_midpoint_recent_pairing_relative)
        value = std::numeric_limits<double>::infinity();
    result.accepted = false;
    result.finite = true;
    result.cfl_ok = dt > 0.0;
    result.gauss_ok = false;
    result.collision_ok = true;
    result.conversion_ok = true;
    result.tail_ok = true;
    result.failure_stage = "joint_midpoint_init";
    result.failure_code = 0;

    if (!initialized_ || !result.cfl_ok || beam_enabled_ ||
        background_tail_enabled_ ||
        field_particle_coupling_config_.mode != FieldParticleCouplingMode::Legacy) {
        result.failure_code = 70;
        result.failure_stage = "joint_midpoint_configuration";
        return result;
    }

    const int nq = Param::Nvmu;
    const size_t local_count = static_cast<size_t>(grid_.nx_local * nq);
    std::vector<double> m_old(local_count, 0.0);
    std::copy(electrons.f.begin() +
                  static_cast<std::ptrdiff_t>(grid_.nghost * nq),
              electrons.f.begin() + static_cast<std::ptrdiff_t>(
                  (grid_.nghost * nq) + local_count), m_old.begin());
    std::vector<double> candidate = m_old;
    EMFields candidate_fields;
    candidate_fields.init(grid_);
    JointPhaseSpaceFluxBundle accepted_bundle;
    std::vector<double> accepted_residual;
    std::vector<double> accepted_phi_residual;
    std::vector<double> accepted_e_local;
    double accepted_norm = std::numeric_limits<double>::infinity();
    double accepted_poisson_residual = std::numeric_limits<double>::infinity();
    const double pairing_tolerance = 1.0e-9;
    result.joint_midpoint_pairing_tolerance = pairing_tolerance;
    double accepted_pairing_residual = 0.0;
    double accepted_pairing_scale = 0.0;
    double accepted_pairing_relative =
        std::numeric_limits<double>::infinity();
    double accepted_scalar_residual = 0.0;
    double accepted_weighted_defect = 0.0;
    double log_pairing_residual = 0.0;
    double log_pairing_relative = std::numeric_limits<double>::infinity();
    double log_weighted_defect = 0.0;
    bool accepted_eval = false;

    const int max_iterations = 20;
    const double residual_tolerance = 1.0e-10;
    const double poisson_tolerance = 1.0e-8;
    const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
    auto dot_global = [&](const std::vector<double>& a,
                          const std::vector<double>& b) {
        long double local = 0.0L;
        for (size_t i = 0; i < a.size(); ++i)
            local += static_cast<long double>(a[i]) * b[i];
        double local_double = static_cast<double>(local);
        double global = 0.0;
        MPI_Allreduce(&local_double, &global, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        return global;
    };
    auto solve_small_normal = [&](const std::vector<std::vector<double>>& h,
                                  int m, double beta,
                                  std::vector<double>& y) {
        double a[4][4] = {};
        double b[4] = {};
        for (int p = 0; p < m; ++p) {
            for (int q = 0; q < m; ++q)
                for (int row = 0; row <= m; ++row)
                    a[p][q] += h[row][p] * h[row][q];
            b[p] = beta * h[0][p];
        }
        for (int p = 0; p < m; ++p) {
            int pivot = p;
            for (int row = p + 1; row < m; ++row)
                if (std::fabs(a[row][p]) > std::fabs(a[pivot][p])) pivot = row;
            if (std::fabs(a[pivot][p]) <= 1.0e-300) return false;
            if (pivot != p) {
                for (int q = p; q < m; ++q) std::swap(a[p][q], a[pivot][q]);
                std::swap(b[p], b[pivot]);
            }
            const double inv = 1.0 / a[p][p];
            for (int q = p; q < m; ++q) a[p][q] *= inv;
            b[p] *= inv;
            for (int row = 0; row < m; ++row) {
                if (row == p) continue;
                const double factor = a[row][p];
                for (int q = p; q < m; ++q) a[row][q] -= factor * a[p][q];
                b[row] -= factor * b[p];
            }
        }
        y.assign(b, b + m);
        return true;
    };
    auto record_iteration = [&](int iteration, int gmres_dimension,
                                double residual_linf, double phi_linf,
                                double alpha, const std::vector<double>& state,
                                int accepted, int failure_code) {
        double local_min = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < state.size(); ++i)
            local_min = std::min(local_min, state[i]);
        double global_min = 0.0;
        MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        VpfpJointPhaseSpaceIterationRecord a4_record;
        a4_record.iteration = iteration;
        a4_record.gmres_dimension = gmres_dimension;
        a4_record.residual_linf = residual_linf;
        a4_record.phi_residual_linf = phi_linf;
        a4_record.line_search_alpha = alpha;
        a4_record.trial_min_mass = global_min;
        a4_record.accepted = accepted;
        a4_record.failure_code = failure_code;
        a4_record.poisson_current_relative = log_pairing_relative;
        a4_record.poisson_current_residual = log_pairing_residual;
        a4_record.weighted_continuity_defect = log_weighted_defect;
        a4_record.phase_converged = residual_linf <= residual_tolerance;
        a4_record.poisson_converged = phi_linf <= poisson_tolerance;
        a4_record.pairing_converged =
            log_pairing_relative <= pairing_tolerance;
        if (result.joint_midpoint_recent_pairing_relative_count < 3) {
            result.joint_midpoint_recent_pairing_relative[
                result.joint_midpoint_recent_pairing_relative_count++] =
                log_pairing_relative;
        } else {
            result.joint_midpoint_recent_pairing_relative[0] =
                result.joint_midpoint_recent_pairing_relative[1];
            result.joint_midpoint_recent_pairing_relative[1] =
                result.joint_midpoint_recent_pairing_relative[2];
            result.joint_midpoint_recent_pairing_relative[2] =
                log_pairing_relative;
        }
        result.joint_midpoint_iterations_log.push_back(a4_record);
    };

    auto evaluate = [&](const std::vector<double>& state,
                        EMFields& eval_fields,
                        JointPhaseSpaceFluxBundle& bundle,
                        std::vector<double>& residual,
                        std::vector<double>& phi_residual_out,
                        double& normalized_norm,
                        double& poisson_residual,
                        double& poisson_current_residual,
                        double& poisson_current_scale,
                        double& poisson_current_relative,
                        double& poisson_scalar_residual,
                        double& weighted_continuity_defect,
                        bool& finite,
                        bool allow_negative_probe,
                        std::vector<double>& e_local_out) -> bool {
        bool local_ok = true;
        normalized_norm = 0.0;
        poisson_residual = 0.0;
        poisson_current_residual = 0.0;
        poisson_current_scale = 0.0;
        poisson_current_relative = std::numeric_limits<double>::infinity();
        poisson_scalar_residual = 0.0;
        weighted_continuity_defect = 0.0;
        try {
            if (state.size() != m_old.size()) local_ok = false;
            for (size_t i = 0; local_ok && i < state.size(); ++i)
                local_ok = std::isfinite(state[i]);
            if (local_ok) {
                std::fill(eval_fields.rho.begin(), eval_fields.rho.end(), 0.0);
                for (int ix = 0; ix < grid_.nx_local; ++ix) {
                    // Stage A-S1 (section 7B.3): stable incremental charge
                    // assembly.  In J1 the ion density is fixed and
                    // Beam/Tail/source are off, so this is algebraically
                    // identical to qe*(ni - n_e^{n+1}) but avoids subtracting
                    // two large near-neutral totals.  Only this J1 candidate
                    // evaluate path changes; EMFields::set_charge_density()
                    // and every other caller keep the production form.
                    const size_t base = static_cast<size_t>(ix) * nq;
                    long double delta_number = 0.0L;
                    for (int q = 0; q < nq; ++q) {
                        delta_number +=
                            static_cast<long double>(state[base + q]) -
                            static_cast<long double>(m_old[base + q]);
                    }
                    const long double delta_rho =
                        -static_cast<long double>(Const::qe) *
                        delta_number /
                        static_cast<long double>(grid_.dx);
                    eval_fields.rho[
                        static_cast<size_t>(grid_.nghost + ix)] =
                        fields.rho[
                            static_cast<size_t>(grid_.nghost + ix)] +
                        static_cast<double>(delta_rho);
                }
                OpenGaussSolveOptions options;
                options.reconstruct_phi = true;
                options.compute_l1 = true;
                options.compute_boundary_audit = true;
                field_solver_.solve(eval_fields, mpi_rank, mpi_size, options);
                poisson_residual = field_solver_.diagnostics().residual_linf;
                local_ok = std::isfinite(poisson_residual);
                std::vector<double> pairing_face;
                if (local_ok) {
                    const bool pairing_ok =
                        field_solver_.build_potential_pairing_field(
                            fields, eval_fields, pairing_face,
                            mpi_rank, mpi_size);
                    result.joint_midpoint_pairing_field_built =
                        result.joint_midpoint_pairing_field_built || pairing_ok;
                    local_ok = pairing_ok;
                }
                std::vector<double> e_pair_cell;
                if (local_ok && pairing_face.size() ==
                        static_cast<size_t>(grid_.nx_local + 1)) {
                    const bool adjoint_ok =
                        JointPhaseSpaceMidpointOperator::
                            build_periodic_x_adjoint_cell_field(
                                grid_, pairing_face, mpi_rank, mpi_size,
                                e_pair_cell);
                    local_ok = local_ok && adjoint_ok;
                } else {
                    local_ok = false;
                }
                std::vector<double> phi_local(static_cast<size_t>(grid_.nx_local), 0.0);
                for (int ix = 0; ix < grid_.nx_local; ++ix) {
                    const size_t id = static_cast<size_t>(grid_.nghost + ix);
                    phi_local[static_cast<size_t>(ix)] =
                        (eval_fields.Ex_face[static_cast<size_t>(ix + 1)] -
                         eval_fields.Ex_face[static_cast<size_t>(ix)]) / grid_.dx -
                        eval_fields.rho[id] / Const::eps0;
                }
                phi_residual_out.swap(phi_local);
                for (size_t i = 0; i < phi_residual_out.size(); ++i)
                    local_ok = local_ok && std::isfinite(phi_residual_out[i]);
                double residual_linf = 0.0;
                double residual_scale = 1.0;
                local_ok = local_ok &&
                    JointPhaseSpaceMidpointOperator::evaluate_local_residual(
                        grid_, electrons.cgrid, m_old, state, e_pair_cell, dt,
                        mpi_rank, mpi_size, bundle, residual, residual_linf,
                        residual_scale, allow_negative_probe);
                e_local_out = e_pair_cell;
                normalized_norm = residual_linf /
                    std::max(1.0, residual_scale);
                // The Gauss residual has units of rho/eps0.  Comparing its
                // absolute value to 1e-8 is invalid at plasma density: the
                // two O(rho/eps0) terms are formed in double precision and
                // leave an absolute cancellation residue even for the exact
                // production Poisson solution.  Use the same local physical
                // scale as the production Gauss gate and pass a dimensionless
                // residual to the J1 Newton convergence test.
                double phi_linf = 0.0;
                double phi_scale = 1.0e-30;
                for (int ix = 0; ix < grid_.nx_local; ++ix) {
                    const size_t face_left = static_cast<size_t>(ix);
                    const size_t face_right = face_left + 1;
                    const size_t cell = static_cast<size_t>(grid_.nghost + ix);
                    const double d_ex_dx =
                        (eval_fields.Ex_face[face_right] -
                         eval_fields.Ex_face[face_left]) / grid_.dx;
                    const double rho_eps0 = eval_fields.rho[cell] / Const::eps0;
                    phi_linf = std::max(phi_linf, std::fabs(d_ex_dx - rho_eps0));
                    phi_scale = std::max(phi_scale,
                        std::max(std::fabs(d_ex_dx), std::fabs(rho_eps0)));
                }
                poisson_residual = phi_linf / phi_scale;

                // A4 candidate metric: all quantities below are computed
                // from this candidate state, its solved fields, and this
                // candidate flux bundle.  No accepted/post-processed state
                // is consulted.
                int metric_ready = local_ok ? 1 : 0;
                MPI_Allreduce(MPI_IN_PLACE, &metric_ready, 1, MPI_INT,
                              MPI_MIN, MPI_COMM_WORLD);
                if (metric_ready) {
                    std::vector<double> rho_delta(eval_fields.rho.size(), 0.0);
                    for (int ix = 0; ix < grid_.nx_local; ++ix) {
                        rho_delta[static_cast<size_t>(grid_.nghost + ix)] =
                            eval_fields.rho[static_cast<size_t>(grid_.nghost + ix)] -
                            fields.rho[static_cast<size_t>(grid_.nghost + ix)];
                    }
                    const OpenPoissonWorkIdentity work =
                        field_solver_.evaluate_work_identity(
                            fields, eval_fields, rho_delta,
                            mpi_rank, mpi_size);
                    poisson_scalar_residual = work.residual;
                    long double local_charge_work = 0.0L;
                    long double local_weighted_defect = 0.0L;
                    for (int iface = 0; iface <= grid_.nx_local; ++iface) {
                        const int global_face = grid_.ix_start + iface;
                        if (iface == 0 && global_face != 0) continue;
                        const double weight =
                            (global_face == 0 || global_face == grid_.nx_global)
                            ? 0.5 : 1.0;
                        local_charge_work +=
                            static_cast<long double>(dt) *
                            static_cast<long double>(weight * grid_.dx) *
                            static_cast<long double>(pairing_face[
                                static_cast<size_t>(iface)]) *
                            static_cast<long double>(bundle.charge_current_face[
                                static_cast<size_t>(iface)]);
                    }
                    for (int ix = 0; ix < grid_.nx_local; ++ix) {
                        const size_t id = static_cast<size_t>(grid_.nghost + ix);
                        const long double r_c =
                            static_cast<long double>(rho_delta[id]) *
                                static_cast<long double>(grid_.dx) +
                            static_cast<long double>(dt) *
                                static_cast<long double>(
                                    bundle.charge_current_face[
                                        static_cast<size_t>(ix) + 1] -
                                    bundle.charge_current_face[
                                        static_cast<size_t>(ix)]);
                        const long double old_phi =
                            static_cast<long double>(fields.phi[id]) +
                            static_cast<long double>(grid_.dx) *
                                (static_cast<long double>(fields.Ex_face[
                                    static_cast<size_t>(ix) + 1]) -
                                 static_cast<long double>(fields.Ex_face[
                                    static_cast<size_t>(ix)])) / 12.0L;
                        const long double new_phi =
                            static_cast<long double>(eval_fields.phi[id]) +
                            static_cast<long double>(grid_.dx) *
                                (static_cast<long double>(eval_fields.Ex_face[
                                    static_cast<size_t>(ix) + 1]) -
                                 static_cast<long double>(eval_fields.Ex_face[
                                    static_cast<size_t>(ix)])) / 12.0L;
                        local_weighted_defect +=
                            0.5L * (old_phi + new_phi) * r_c;
                    }
                    double global_charge_work = 0.0;
                    const double local_charge_work_double =
                        static_cast<double>(local_charge_work);
                    MPI_Allreduce(&local_charge_work_double,
                                  &global_charge_work, 1, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    const double local_weighted_double =
                        static_cast<double>(local_weighted_defect);
                    double global_weighted_defect = 0.0;
                    MPI_Allreduce(&local_weighted_double,
                                  &global_weighted_defect, 1, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    poisson_current_residual =
                        work.field_energy_change - work.electrode_work +
                        global_charge_work;
                    poisson_current_scale = std::max(
                        1.0, std::max(std::fabs(work.field_energy_change),
                        std::max(std::fabs(work.electrode_work),
                                 std::fabs(global_charge_work))));
                    poisson_current_relative =
                        std::fabs(poisson_current_residual) /
                        poisson_current_scale;
                    weighted_continuity_defect = global_weighted_defect;
                    local_ok = work.finite &&
                        std::isfinite(poisson_scalar_residual) &&
                        std::isfinite(poisson_current_residual) &&
                        std::isfinite(poisson_current_scale) &&
                        std::isfinite(poisson_current_relative) &&
                        std::isfinite(weighted_continuity_defect);
                }
            }
        } catch (const std::exception&) {
            local_ok = false;
        }
        int global_ok = local_ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        finite = global_ok != 0;
        if (!finite) return false;
        double global_norm = 0.0;
        MPI_Allreduce(&normalized_norm, &global_norm, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        normalized_norm = global_norm;
        double global_poisson = 0.0;
        MPI_Allreduce(&poisson_residual, &global_poisson, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        poisson_residual = global_poisson;
        MPI_Allreduce(MPI_IN_PLACE, &poisson_current_residual, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &poisson_current_scale, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &poisson_current_relative, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &poisson_scalar_residual, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return true;
    };

    bool finite = false;
    if (!evaluate(candidate, candidate_fields, accepted_bundle, accepted_residual,
                  accepted_phi_residual,
                  accepted_norm, accepted_poisson_residual,
                  accepted_pairing_residual, accepted_pairing_scale,
                  accepted_pairing_relative, accepted_scalar_residual,
                  accepted_weighted_defect, finite,
                  // The initial candidate is the previously accepted state.
                  // Accepted J1 states may contain signed roundoff-level mass
                  // under the final code-76 tolerance, so the next step must
                  // evaluate that same signed residual domain.  Positivity is
                  // still enforced after convergence by the unchanged code-76
                  // acceptance gate below; no clipping or tolerance change is
                  // performed here.
                  true,
                  accepted_e_local)) {
        result.finite = false;
        result.failure_code = 71;
        result.failure_stage = "joint_midpoint_initial_residual";
        return result;
    }
    log_pairing_residual = accepted_pairing_residual;
    log_pairing_relative = accepted_pairing_relative;
    log_weighted_defect = accepted_weighted_defect;
    result.joint_midpoint_candidate_poisson_current_residual =
        accepted_pairing_residual;
    result.joint_midpoint_candidate_poisson_current_scale =
        accepted_pairing_scale;
    result.joint_midpoint_candidate_poisson_current_relative =
        accepted_pairing_relative;
    result.joint_midpoint_candidate_poisson_scalar_residual =
        accepted_scalar_residual;
    result.joint_midpoint_candidate_weighted_continuity_defect =
        accepted_weighted_defect;
    record_iteration(0, 0, accepted_norm, accepted_poisson_residual, 0.0,
                     candidate, 0, 0);
    // These values describe the last finite candidate even if Newton later
    // rejects every step.  They must not remain at their zero initialization
    // and make a line-search failure look like a Gauss failure.
    result.joint_midpoint_residual_linf = accepted_norm;
    result.joint_midpoint_poisson_residual_linf = accepted_poisson_residual;
    result.gauss_ok = accepted_poisson_residual <= poisson_tolerance;

    for (int iter = 0; iter < max_iterations; ++iter) {
        result.joint_midpoint_iterations = iter + 1;
        const bool accepted_phase_converged =
            accepted_norm <= residual_tolerance;
        const bool accepted_poisson_converged =
            accepted_poisson_residual <= poisson_tolerance;
        const bool accepted_pairing_converged =
            accepted_pairing_relative <= pairing_tolerance;
        result.joint_midpoint_phase_converged =
            accepted_phase_converged ? 1 : 0;
        result.joint_midpoint_poisson_converged =
            accepted_poisson_converged ? 1 : 0;
        result.joint_midpoint_pairing_converged =
            accepted_pairing_converged ? 1 : 0;
        if (accepted_phase_converged && accepted_poisson_converged &&
            accepted_pairing_converged) {
            accepted_eval = true;
            result.joint_midpoint_converged = true;
            log_pairing_residual = accepted_pairing_residual;
            log_pairing_relative = accepted_pairing_relative;
            log_weighted_defect = accepted_weighted_defect;
            record_iteration(iter + 1, 0, accepted_norm,
                             accepted_poisson_residual, 1.0, candidate, 1, 0);
            break;
        }

        std::vector<double> preconditioned_residual;
        if (!JointPhaseSpaceMidpointOperator::apply_local_block_diagonal_preconditioner(
                electrons.cgrid, accepted_e_local, grid_.dx, dt,
                accepted_residual, preconditioned_residual)) {
            result.failure_code = 72;
            result.failure_stage = "joint_midpoint_preconditioner";
            record_iteration(iter + 1, 0, accepted_norm,
                             accepted_poisson_residual, 0.0, candidate, 0, 72);
            return result;
        }
        // Newton requires J delta = -R.  GMRES solves the left-preconditioned
        // system P^{-1} J delta = -P^{-1} R, so the Krylov start vector must
        // carry the negative residual.  Using +P^{-1}R constructs an ascent
        // direction and makes every Armijo-style trial fail at iteration one.
        for (size_t i = 0; i < preconditioned_residual.size(); ++i)
            preconditioned_residual[i] = -preconditioned_residual[i];
        double local_state_norm = 0.0;
        for (size_t i = 0; i < candidate.size(); ++i)
            local_state_norm = std::max(local_state_norm, std::fabs(candidate[i]));
        double global_state_norm = 0.0;
        MPI_Allreduce(&local_state_norm, &global_state_norm, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        const int restart_dim = 4;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> hess(
            static_cast<size_t>(restart_dim + 1),
            std::vector<double>(static_cast<size_t>(restart_dim), 0.0));
        const double beta0 = std::sqrt(std::max(
            0.0, dot_global(preconditioned_residual, preconditioned_residual)));
        if (!(beta0 > 0.0) || !std::isfinite(beta0)) {
            result.failure_code = 72;
            result.failure_stage = "joint_midpoint_gmres_initial_vector";
            record_iteration(iter + 1, 0, accepted_norm,
                             accepted_poisson_residual, 0.0, candidate, 0, 72);
            return result;
        }
        std::vector<double> v0(preconditioned_residual.size(), 0.0);
        for (size_t i = 0; i < v0.size(); ++i) v0[i] = preconditioned_residual[i] / beta0;
        basis.push_back(v0);
        int gmres_dimension = 0;
        for (int k = 0; k < restart_dim; ++k) {
            std::vector<double> probe = candidate;
            const double fd = sqrt_eps * std::max(1.0, global_state_norm);
            for (size_t i = 0; i < probe.size(); ++i)
                probe[i] += fd * basis[static_cast<size_t>(k)][i];
            EMFields probe_fields;
            probe_fields.init(grid_);
            JointPhaseSpaceFluxBundle probe_bundle;
            std::vector<double> probe_residual;
            std::vector<double> probe_phi_residual;
            std::vector<double> probe_e_local;
            double probe_norm = 0.0;
            double probe_poisson = 0.0;
            double probe_pairing_residual = 0.0;
            double probe_pairing_scale = 0.0;
            double probe_pairing_relative =
                std::numeric_limits<double>::infinity();
            double probe_scalar_residual = 0.0;
            double probe_weighted_defect = 0.0;
            bool probe_finite = false;
            if (!evaluate(probe, probe_fields, probe_bundle, probe_residual,
                          probe_phi_residual, probe_norm, probe_poisson,
                          probe_pairing_residual, probe_pairing_scale,
                          probe_pairing_relative, probe_scalar_residual,
                          probe_weighted_defect,
                          probe_finite, true, probe_e_local)) {
                result.failure_code = 72;
                result.failure_stage = "joint_midpoint_jacobian_probe";
                record_iteration(iter + 1, k + 1, accepted_norm,
                                 accepted_poisson_residual, 0.0, candidate, 0, 72);
                return result;
            }
            std::vector<double> jv(probe_residual.size(), 0.0);
            for (size_t i = 0; i < jv.size(); ++i)
                jv[i] = (probe_residual[i] - accepted_residual[i]) / fd;
            std::vector<double> z;
            if (!JointPhaseSpaceMidpointOperator::apply_local_block_diagonal_preconditioner(
                    electrons.cgrid, accepted_e_local, grid_.dx, dt, jv, z)) {
                result.failure_code = 72;
                result.failure_stage = "joint_midpoint_gmres_preconditioner";
                record_iteration(iter + 1, k + 1, accepted_norm,
                                 accepted_poisson_residual, 0.0, candidate, 0, 72);
                return result;
            }
            for (int j = 0; j <= k; ++j) {
                hess[static_cast<size_t>(j)][static_cast<size_t>(k)] =
                    dot_global(z, basis[static_cast<size_t>(j)]);
                for (size_t i = 0; i < z.size(); ++i)
                    z[i] -= hess[static_cast<size_t>(j)][static_cast<size_t>(k)] *
                            basis[static_cast<size_t>(j)][i];
            }
            const double hnext = std::sqrt(std::max(0.0, dot_global(z, z)));
            hess[static_cast<size_t>(k + 1)][static_cast<size_t>(k)] = hnext;
            gmres_dimension = k + 1;
            if (hnext <= 1.0e-14 || k + 1 == restart_dim) break;
            for (size_t i = 0; i < z.size(); ++i) z[i] /= hnext;
            basis.push_back(z);
        }
        std::vector<double> y;
        if (!solve_small_normal(hess, gmres_dimension, beta0, y)) {
            result.failure_code = 72;
            result.failure_stage = "joint_midpoint_gmres_small_solve";
            record_iteration(iter + 1, gmres_dimension, accepted_norm,
                             accepted_poisson_residual, 0.0, candidate, 0, 72);
            return result;
        }
        std::vector<double> step(candidate.size(), 0.0);
        for (int j = 0; j < gmres_dimension; ++j)
            for (size_t i = 0; i < step.size(); ++i)
                step[i] += y[static_cast<size_t>(j)] * basis[static_cast<size_t>(j)][i];

        bool line_search_accepted = false;
        const bool pairing_search = accepted_phase_converged &&
            !accepted_pairing_converged;
        double lambda = 1.0;
        for (int ls = 0; ls < 12; ++ls) {
            std::vector<double> trial = candidate;
            for (size_t i = 0; i < trial.size(); ++i)
                trial[i] += lambda * step[i];
            EMFields trial_fields;
            trial_fields.init(grid_);
            JointPhaseSpaceFluxBundle trial_bundle;
            std::vector<double> trial_residual;
            std::vector<double> trial_phi_residual;
            std::vector<double> trial_e_local;
            double trial_norm = 0.0;
            double trial_poisson = 0.0;
            double trial_pairing_residual = 0.0;
            double trial_pairing_scale = 0.0;
            double trial_pairing_relative =
                std::numeric_limits<double>::infinity();
            double trial_scalar_residual = 0.0;
            double trial_weighted_defect = 0.0;
            bool trial_finite = false;
            if (evaluate(trial, trial_fields, trial_bundle, trial_residual,
                         trial_phi_residual,
                         trial_norm, trial_poisson,
                         trial_pairing_residual, trial_pairing_scale,
                         trial_pairing_relative, trial_scalar_residual,
                         trial_weighted_defect, trial_finite,
                         // A line-search candidate is an algebraic Newton
                         // trial.  It must use the same signed residual
                         // domain as the Jacobian probe; final positivity is
                         // enforced only after nonlinear convergence by the
                         // code-76 acceptance gate below.
                         true,
                         trial_e_local) &&
                trial_finite &&
                (pairing_search
                    ? (trial_norm <= residual_tolerance &&
                       trial_poisson <= poisson_tolerance &&
                       trial_pairing_relative < accepted_pairing_relative)
                    : (trial_norm < accepted_norm &&
                       trial_poisson <= std::max(
                           poisson_tolerance,
                           accepted_poisson_residual * 2.0)))) {
                candidate.swap(trial);
                candidate_fields = trial_fields;
                accepted_bundle = trial_bundle;
                accepted_residual.swap(trial_residual);
                accepted_phi_residual.swap(trial_phi_residual);
                accepted_e_local.swap(trial_e_local);
                accepted_norm = trial_norm;
                accepted_poisson_residual = trial_poisson;
                accepted_pairing_residual = trial_pairing_residual;
                accepted_pairing_scale = trial_pairing_scale;
                accepted_pairing_relative = trial_pairing_relative;
                accepted_scalar_residual = trial_scalar_residual;
                accepted_weighted_defect = trial_weighted_defect;
                result.joint_midpoint_candidate_poisson_current_residual =
                    accepted_pairing_residual;
                result.joint_midpoint_candidate_poisson_current_scale =
                    accepted_pairing_scale;
                result.joint_midpoint_candidate_poisson_current_relative =
                    accepted_pairing_relative;
                result.joint_midpoint_candidate_poisson_scalar_residual =
                    accepted_scalar_residual;
                result.joint_midpoint_candidate_weighted_continuity_defect =
                    accepted_weighted_defect;
                log_pairing_residual = accepted_pairing_residual;
                log_pairing_relative = accepted_pairing_relative;
                log_weighted_defect = accepted_weighted_defect;
                result.joint_midpoint_residual_linf = accepted_norm;
                result.joint_midpoint_poisson_residual_linf =
                    accepted_poisson_residual;
                result.gauss_ok = accepted_poisson_residual <=
                    poisson_tolerance;
                line_search_accepted = true;
                record_iteration(iter + 1, gmres_dimension, trial_norm,
                                 trial_poisson, lambda, candidate, 1, 0);
                break;
            }
            lambda *= 0.5;
        }
        // A tiny normal-equation GMRES solve can lose its descent property
        // when the Krylov columns are nearly dependent.  Before declaring a
        // nonlinear failure, try the already available left-preconditioned
        // Newton residual direction (-P^{-1}R) under the identical strict
        // line-search gates.  This is a solver safeguard only; it does not
        // alter the residual, fluxes, tolerances, or accepted-state rules.
        if (!line_search_accepted) {
            lambda = 1.0;
            for (int ls = 0; ls < 12; ++ls) {
                std::vector<double> trial = candidate;
                for (size_t i = 0; i < trial.size(); ++i)
                    trial[i] += lambda * preconditioned_residual[i];
                EMFields trial_fields;
                trial_fields.init(grid_);
                JointPhaseSpaceFluxBundle trial_bundle;
                std::vector<double> trial_residual;
                std::vector<double> trial_phi_residual;
                std::vector<double> trial_e_local;
                double trial_norm = 0.0;
                double trial_poisson = 0.0;
                double trial_pairing_residual = 0.0;
                double trial_pairing_scale = 0.0;
                double trial_pairing_relative =
                    std::numeric_limits<double>::infinity();
                double trial_scalar_residual = 0.0;
                double trial_weighted_defect = 0.0;
                bool trial_finite = false;
                if (evaluate(trial, trial_fields, trial_bundle,
                             trial_residual, trial_phi_residual,
                             trial_norm, trial_poisson,
                             trial_pairing_residual, trial_pairing_scale,
                             trial_pairing_relative, trial_scalar_residual,
                             trial_weighted_defect, trial_finite,
                             true, trial_e_local) && trial_finite &&
                    (pairing_search
                        ? (trial_norm <= residual_tolerance &&
                           trial_poisson <= poisson_tolerance &&
                           trial_pairing_relative < accepted_pairing_relative)
                        : (trial_norm < accepted_norm &&
                           trial_poisson <= std::max(
                               poisson_tolerance,
                               accepted_poisson_residual * 2.0)))) {
                    candidate.swap(trial);
                    candidate_fields = trial_fields;
                    accepted_bundle = trial_bundle;
                    accepted_residual.swap(trial_residual);
                    accepted_phi_residual.swap(trial_phi_residual);
                    accepted_e_local.swap(trial_e_local);
                    accepted_norm = trial_norm;
                    accepted_poisson_residual = trial_poisson;
                    accepted_pairing_residual = trial_pairing_residual;
                    accepted_pairing_scale = trial_pairing_scale;
                    accepted_pairing_relative = trial_pairing_relative;
                    accepted_scalar_residual = trial_scalar_residual;
                    accepted_weighted_defect = trial_weighted_defect;
                    result.joint_midpoint_candidate_poisson_current_residual =
                        accepted_pairing_residual;
                    result.joint_midpoint_candidate_poisson_current_scale =
                        accepted_pairing_scale;
                    result.joint_midpoint_candidate_poisson_current_relative =
                        accepted_pairing_relative;
                    result.joint_midpoint_candidate_poisson_scalar_residual =
                        accepted_scalar_residual;
                    result.joint_midpoint_candidate_weighted_continuity_defect =
                        accepted_weighted_defect;
                    log_pairing_residual = accepted_pairing_residual;
                    log_pairing_relative = accepted_pairing_relative;
                    log_weighted_defect = accepted_weighted_defect;
                    result.joint_midpoint_residual_linf = accepted_norm;
                    result.joint_midpoint_poisson_residual_linf =
                        accepted_poisson_residual;
                    result.gauss_ok = accepted_poisson_residual <=
                        poisson_tolerance;
                    line_search_accepted = true;
                    record_iteration(iter + 1, -1, trial_norm,
                                     trial_poisson, lambda, candidate, 1, 0);
                    break;
                }
                lambda *= 0.5;
            }
        }
        if (!line_search_accepted) {
            result.failure_code = pairing_search ? 74 : 73;
            result.failure_stage = pairing_search
                ? "joint_midpoint_poisson_current_not_converged"
                : "joint_midpoint_line_search";
            result.joint_midpoint_residual_linf = accepted_norm;
            result.joint_midpoint_poisson_residual_linf = accepted_poisson_residual;
            log_pairing_residual = accepted_pairing_residual;
            log_pairing_relative = accepted_pairing_relative;
            log_weighted_defect = accepted_weighted_defect;
            record_iteration(iter + 1, gmres_dimension, accepted_norm,
                             accepted_poisson_residual, 0.0, candidate, 0,
                             result.failure_code);
            return result;
        }
    }

    if (!accepted_eval && accepted_norm <= residual_tolerance &&
        accepted_poisson_residual <= poisson_tolerance &&
        accepted_pairing_relative <= pairing_tolerance) accepted_eval = true;
    result.joint_midpoint_residual_linf = accepted_norm;
    result.joint_midpoint_poisson_residual_linf = accepted_poisson_residual;
    result.joint_midpoint_candidate_poisson_current_residual =
        accepted_pairing_residual;
    result.joint_midpoint_candidate_poisson_current_scale =
        accepted_pairing_scale;
    result.joint_midpoint_candidate_poisson_current_relative =
        accepted_pairing_relative;
    result.joint_midpoint_candidate_poisson_scalar_residual =
        accepted_scalar_residual;
    result.joint_midpoint_candidate_weighted_continuity_defect =
        accepted_weighted_defect;
    result.joint_midpoint_phase_converged =
        accepted_norm <= residual_tolerance ? 1 : 0;
    result.joint_midpoint_poisson_converged =
        accepted_poisson_residual <= poisson_tolerance ? 1 : 0;
    result.joint_midpoint_pairing_converged =
        accepted_pairing_relative <= pairing_tolerance ? 1 : 0;
    if (!accepted_eval) {
        result.failure_code = 74;
        result.failure_stage =
            result.joint_midpoint_phase_converged &&
            result.joint_midpoint_poisson_converged &&
            !result.joint_midpoint_pairing_converged
            ? "joint_midpoint_poisson_current_not_converged"
            : "joint_midpoint_not_converged";
        log_pairing_residual = accepted_pairing_residual;
        log_pairing_relative = accepted_pairing_relative;
        log_weighted_defect = accepted_weighted_defect;
        record_iteration(result.joint_midpoint_iterations, 4,
                         accepted_norm, accepted_poisson_residual,
                         0.0, candidate, 0, result.failure_code);
        return result;
    }

    double local_min_candidate = std::numeric_limits<double>::infinity();
    double local_max_candidate = 0.0;
    std::size_t local_first_negative = candidate.size();
    for (std::size_t i = 0; i < candidate.size(); ++i) {
        local_min_candidate = std::min(local_min_candidate, candidate[i]);
        local_max_candidate = std::max(local_max_candidate,
                                       std::fabs(candidate[i]));
        if (candidate[i] < 0.0 && local_first_negative == candidate.size())
            local_first_negative = i;
    }
    double candidate_bounds[2] = { local_min_candidate, local_max_candidate };
    MPI_Allreduce(MPI_IN_PLACE, candidate_bounds, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, candidate_bounds + 1, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double negative_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, candidate_bounds[1]);
    result.joint_midpoint_min_mass = candidate_bounds[0];
    result.joint_midpoint_max_mass = candidate_bounds[1];
    if (candidate_bounds[0] < -negative_tolerance) {
        unsigned long long local_index =
            local_first_negative == candidate.size()
            ? std::numeric_limits<unsigned long long>::max()
            : static_cast<unsigned long long>(local_first_negative);
        unsigned long long global_index = 0;
        MPI_Allreduce(&local_index, &global_index, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        result.failure_code = 76;
        result.failure_stage = "joint_midpoint_negative_solution";
        result.input_min = candidate_bounds[0];
        result.input_max = candidate_bounds[1];
        if (global_index != std::numeric_limits<unsigned long long>::max()) {
            const std::size_t phase_size = static_cast<std::size_t>(Param::Nvmu);
            result.failing_ix = static_cast<int>(global_index / phase_size);
            const std::size_t velocity_index = global_index % phase_size;
            result.failing_iupar = static_cast<int>(velocity_index / Param::Nmu);
            result.failing_iuperp = static_cast<int>(velocity_index % Param::Nmu);
        }
        return result;
    }

    const double number_before = global_sum(electrons.total_particle_number());
    const double field_before = field_energy(fields);
    const double kinetic_before = global_sum(electrons.total_kinetic_energy());
    std::vector<double> rho_delta(fields.rho.size(), 0.0);
    for (size_t i = 0; i < rho_delta.size() && i < candidate_fields.rho.size(); ++i)
        rho_delta[i] = candidate_fields.rho[i] - fields.rho[i];
    const OpenPoissonWorkIdentity poisson_work =
        field_solver_.evaluate_work_identity(fields, candidate_fields,
                                             rho_delta, mpi_rank, mpi_size);
    std::vector<double> final_pairing_face;
    const bool final_pairing_ok =
        field_solver_.build_potential_pairing_field(
            fields, candidate_fields, final_pairing_face,
            mpi_rank, mpi_size);
    if (!final_pairing_ok || final_pairing_face.size() !=
            static_cast<size_t>(grid_.nx_local + 1)) {
        result.failure_code = 71;
        result.failure_stage = "joint_midpoint_final_pairing_field";
        return result;
    }
    result.joint_midpoint_pairing_field_built = true;
    double local_delta_ke_u = 0.0;
    double local_delta_ke_x = 0.0;
    long double local_u_face_work = 0.0L;
    long double local_force_current_work = 0.0L;
    long double local_charge_current_work = 0.0L;
    long double local_charge_current_work_interior = 0.0L;
    long double local_charge_current_work_endpoint = 0.0L;
    const int nmu = static_cast<int>(electrons.cgrid.uperp_cells.size());
    const int nupar = static_cast<int>(electrons.cgrid.upar_cells.size());
    const std::vector<double> hamiltonian_velocity =
        JointPhaseSpaceMidpointOperator::build_hamiltonian_velocity(
            electrons.cgrid);
    std::vector<double> final_pairing_cell;
    const bool final_adjoint_ok =
        JointPhaseSpaceMidpointOperator::
            build_periodic_x_adjoint_cell_field(
                grid_, final_pairing_face, mpi_rank, mpi_size,
                final_pairing_cell);
    if (!final_adjoint_ok) {
        result.failure_code = 71;
        result.failure_stage = "joint_midpoint_final_x_adjoint_field";
        return result;
    }
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        // accepted_bundle and candidate are rank-local slabs here too.
        const size_t base = static_cast<size_t>(ix) * nq;
        for (int q = 0; q < nq; ++q) {
            local_delta_ke_x +=
                electrons.cgrid.kinetic_energy[static_cast<size_t>(q)] *
                accepted_bundle.mass_delta_x[base + static_cast<size_t>(q)];
            local_delta_ke_u +=
                electrons.cgrid.kinetic_energy[static_cast<size_t>(q)] *
                accepted_bundle.mass_delta_u[base + static_cast<size_t>(q)];
            const double midpoint_mass = 0.5 *
                (m_old[base + static_cast<size_t>(q)] +
                 candidate[base + static_cast<size_t>(q)]);
            local_force_current_work += static_cast<long double>(dt) *
                static_cast<long double>(final_pairing_cell[
                    static_cast<size_t>(ix)]) *
                static_cast<long double>(-Const::qe) *
                static_cast<long double>(hamiltonian_velocity[
                    static_cast<size_t>(q)]) *
                static_cast<long double>(midpoint_mass);
        }
        for (int jf = 1; jf < nupar; ++jf) {
            for (int k = 0; k < nmu; ++k) {
                const size_t left = static_cast<size_t>((jf - 1) * nmu + k);
                const size_t right = static_cast<size_t>(jf * nmu + k);
                const double delta_k =
                    electrons.cgrid.kinetic_energy[right] -
                    electrons.cgrid.kinetic_energy[left];
                const size_t face_id =
                    (static_cast<size_t>(ix) *
                     static_cast<size_t>(nupar + 1) +
                     static_cast<size_t>(jf)) *
                    static_cast<size_t>(nmu) +
                    static_cast<size_t>(k);
                local_u_face_work += static_cast<long double>(dt) *
                    static_cast<long double>(delta_k) *
                    static_cast<long double>(accepted_bundle.u_flux_rate[face_id]);
            }
        }
    }
    for (int iface = 0; iface <= grid_.nx_local; ++iface) {
        const int global_face = grid_.ix_start + iface;
        // A shared MPI face is owned by the rank on its left.  Physical
        // endpoints retain the half-cell quadrature weights.
        if (iface == 0 && global_face != 0) continue;
        const double weight =
            (global_face == 0 || global_face == grid_.nx_global) ? 0.5 : 1.0;
        local_charge_current_work += static_cast<long double>(dt) *
            static_cast<long double>(weight * grid_.dx) *
            static_cast<long double>(final_pairing_face[
                static_cast<size_t>(iface)]) *
            static_cast<long double>(accepted_bundle.charge_current_face[
                static_cast<size_t>(iface)]);
        if (global_face == 0 || global_face == grid_.nx_global) {
            local_charge_current_work_endpoint += static_cast<long double>(dt) *
                static_cast<long double>(weight * grid_.dx) *
                static_cast<long double>(final_pairing_face[
                    static_cast<size_t>(iface)]) *
                static_cast<long double>(accepted_bundle.charge_current_face[
                    static_cast<size_t>(iface)]);
        } else {
            local_charge_current_work_interior += static_cast<long double>(dt) *
                static_cast<long double>(weight * grid_.dx) *
                static_cast<long double>(final_pairing_face[
                    static_cast<size_t>(iface)]) *
                static_cast<long double>(accepted_bundle.charge_current_face[
                    static_cast<size_t>(iface)]);
        }
    }
    result.joint_midpoint_delta_k_x = global_sum(local_delta_ke_x);
    const double delta_ke_u = global_sum(local_delta_ke_u);
    result.joint_midpoint_delta_k_u = delta_ke_u;
    double global_u_face_work = 0.0;
    double local_u_face_work_double = static_cast<double>(local_u_face_work);
    MPI_Allreduce(&local_u_face_work_double, &global_u_face_work, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    result.joint_midpoint_u_face_work = global_u_face_work;
    double global_force_current_work = 0.0;
    const double local_force_current_work_double =
        static_cast<double>(local_force_current_work);
    MPI_Allreduce(&local_force_current_work_double,
                  &global_force_current_work, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.joint_midpoint_force_current_work = global_force_current_work;
    double global_charge_current_work = 0.0;
    double local_charge_current_work_double =
        static_cast<double>(local_charge_current_work);
    MPI_Allreduce(&local_charge_current_work_double,
                  &global_charge_current_work, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.joint_midpoint_charge_current_work = global_charge_current_work;
    double global_charge_current_work_interior = 0.0;
    const double local_charge_current_work_interior_double =
        static_cast<double>(local_charge_current_work_interior);
    MPI_Allreduce(&local_charge_current_work_interior_double,
                  &global_charge_current_work_interior, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    result.joint_midpoint_charge_current_work_interior =
        global_charge_current_work_interior;
    double global_charge_current_work_endpoint = 0.0;
    const double local_charge_current_work_endpoint_double =
        static_cast<double>(local_charge_current_work_endpoint);
    MPI_Allreduce(&local_charge_current_work_endpoint_double,
                  &global_charge_current_work_endpoint, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    result.joint_midpoint_charge_current_work_endpoint =
        global_charge_current_work_endpoint;
    result.joint_midpoint_poisson_potential_charge_work =
        poisson_work.potential_charge_work;
    result.joint_midpoint_current_pair_residual =
        result.joint_midpoint_u_face_work -
        result.joint_midpoint_charge_current_work;
    result.joint_midpoint_force_charge_residual =
        result.joint_midpoint_force_current_work -
        result.joint_midpoint_charge_current_work;
    std::vector<double> force_current_cell(
        static_cast<size_t>(grid_.nx_local), 0.0);
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        long double sum = 0.0L;
        const size_t base = static_cast<size_t>(ix) * nq;
        for (int q = 0; q < nq; ++q) {
            const double midpoint_mass = 0.5 *
                (m_old[base + static_cast<size_t>(q)] +
                 candidate[base + static_cast<size_t>(q)]);
            sum += static_cast<long double>(
                       hamiltonian_velocity[static_cast<size_t>(q)]) *
                   static_cast<long double>(midpoint_mass);
        }
        force_current_cell[static_cast<size_t>(ix)] =
            (-Const::qe) * static_cast<double>(sum) / grid_.dx;
    }
    long double local_naive_force_current_work = 0.0L;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        const double naive_pairing_cell = 0.5 *
            (final_pairing_face[static_cast<size_t>(ix)] +
             final_pairing_face[static_cast<size_t>(ix + 1)]);
        local_naive_force_current_work +=
            static_cast<long double>(dt) *
            static_cast<long double>(grid_.dx) *
            static_cast<long double>(naive_pairing_cell) *
            static_cast<long double>(
                force_current_cell[static_cast<size_t>(ix)]);
    }
    double global_naive_force_current_work = 0.0;
    const double local_naive_force_current_work_double =
        static_cast<double>(local_naive_force_current_work);
    MPI_Allreduce(&local_naive_force_current_work_double,
                  &global_naive_force_current_work, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.joint_midpoint_naive_force_current_work =
        global_naive_force_current_work;
    double endpoint_local[4] = {0.0, 0.0, 0.0, 0.0};
    if (mpi_rank == 0) {
        endpoint_local[0] = final_pairing_face.front();
        endpoint_local[2] = force_current_cell.front();
    }
    if (mpi_rank == mpi_size - 1) {
        endpoint_local[1] = final_pairing_face.back();
        endpoint_local[3] = force_current_cell.back();
    }
    double endpoint_global[4] = {0.0, 0.0, 0.0, 0.0};
    MPI_Allreduce(endpoint_local, endpoint_global, 4, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.joint_midpoint_pairing_face_left = endpoint_global[0];
    result.joint_midpoint_pairing_face_right = endpoint_global[1];
    result.joint_midpoint_force_current_first_cell = endpoint_global[2];
    result.joint_midpoint_force_current_last_cell = endpoint_global[3];
    const double seam_pred =
        0.25 * dt * grid_.dx *
        (endpoint_global[0] - endpoint_global[1]) *
        (endpoint_global[2] - endpoint_global[3]);
    result.joint_midpoint_seam_predicted_residual = seam_pred;
    result.joint_midpoint_seam_prediction_error =
        result.joint_midpoint_force_charge_residual - seam_pred;
    result.joint_midpoint_poisson_transport_residual =
        poisson_work.field_energy_change - poisson_work.electrode_work +
        result.joint_midpoint_charge_current_work;
    result.joint_midpoint_field_energy_change =
        poisson_work.field_energy_change;
    result.joint_midpoint_electrode_work = poisson_work.electrode_work;
    result.joint_midpoint_domain_energy_change =
        result.joint_midpoint_delta_k_x +
        result.joint_midpoint_delta_k_u +
        poisson_work.field_energy_change;
    result.joint_midpoint_energy_residual =
        result.joint_midpoint_current_pair_residual +
        result.joint_midpoint_poisson_transport_residual;
    // Stage-A1 read-only residual decomposition (docs/VPFP_F10情形A_连续性
    // Poisson功配对严格修复实施方案.md sections 2 and 5).  Diagnostics only:
    // this block never modifies flux, residual, Newton, Poisson, the energy
    // gate, acceptance logic or dt.  It uses only the final accepted
    // candidate's accepted_residual, accepted_bundle, fields and
    // candidate_fields; no current is rebuilt and no extra time advance.
    result.joint_midpoint_poisson_scalar_identity_residual =
        poisson_work.residual;
    // Minimal read-only Poisson identity acceptance (production residual and
    // scale verbatim, established 8192*eps gate; no relaxation).  Diagnostic
    // only: it never alters the failure path.
    result.joint_midpoint_poisson_identity_scale = poisson_work.scale;
    // Stage A-FS-R1 (section 7C.10.6): dual Gate F sets; the legacy 8192
    // gate stays visible, the production alias maps to 16384.
    const double afsr1_eps = std::numeric_limits<double>::epsilon();
    result.joint_midpoint_poisson_identity_roundoff_bound_8192 =
        8192.0 * afsr1_eps * poisson_work.scale;
    result.joint_midpoint_poisson_identity_roundoff_bound_16384 =
        16384.0 * afsr1_eps * poisson_work.scale;
    result.joint_midpoint_poisson_identity_roundoff_bound =
        result.joint_midpoint_poisson_identity_roundoff_bound_8192;
    result.joint_midpoint_poisson_identity_finite =
        poisson_work.finite ? 1 : 0;
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio_8192 =
        poisson_work.finite &&
            result.joint_midpoint_poisson_identity_roundoff_bound_8192 > 0.0
        ? std::fabs(poisson_work.residual) /
              result.joint_midpoint_poisson_identity_roundoff_bound_8192
        : std::numeric_limits<double>::infinity();
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio_16384 =
        poisson_work.finite &&
            result.joint_midpoint_poisson_identity_roundoff_bound_16384 >
                0.0
        ? std::fabs(poisson_work.residual) /
              result.joint_midpoint_poisson_identity_roundoff_bound_16384
        : std::numeric_limits<double>::infinity();
    result.joint_midpoint_poisson_identity_residual_to_bound_ratio =
        result.joint_midpoint_poisson_identity_residual_to_bound_ratio_8192;
    result.joint_midpoint_poisson_scalar_identity_pass_8192 =
        poisson_work.finite &&
        result.joint_midpoint_poisson_identity_residual_to_bound_ratio_8192 <=
            1.0
        ? 1
        : 0;
    result.joint_midpoint_poisson_scalar_identity_pass_16384 =
        poisson_work.finite &&
        result.joint_midpoint_poisson_identity_residual_to_bound_ratio_16384 <=
            1.0
        ? 1
        : 0;
    // Compatibility alias maps to the production 16384 gate.
    result.joint_midpoint_poisson_scalar_identity_pass =
        result.joint_midpoint_poisson_scalar_identity_pass_16384;
    // Summation-error diagnosis for the three local integrals inside
    // evaluate_work_identity(): read-only absolute accumulation scales,
    // computed with the same per-cell formulas but |.| inside the sum.
    {
        const int as1b_ng = grid_.nghost;
        long double as1b_term0 = 0.0L;
        long double as1b_term1 = 0.0L;
        long double as1b_term2 = 0.0L;
        const long double eps0_ld =
            static_cast<long double>(Const::eps0);
        for (int ix = 0; ix < grid_.nx_local; ++ix) {
            const long double o_l = static_cast<long double>(
                fields.Ex_face[static_cast<size_t>(ix)]);
            const long double o_r = static_cast<long double>(
                fields.Ex_face[static_cast<size_t>(ix) + 1]);
            const long double n_l = static_cast<long double>(
                candidate_fields.Ex_face[static_cast<size_t>(ix)]);
            const long double n_r = static_cast<long double>(
                candidate_fields.Ex_face[static_cast<size_t>(ix) + 1]);
            as1b_term0 += std::fabs(eps0_ld *
                static_cast<long double>(grid_.dx) *
                (o_l * o_l + o_l * o_r + o_r * o_r) / 6.0L);
            as1b_term1 += std::fabs(eps0_ld *
                static_cast<long double>(grid_.dx) *
                (n_l * n_l + n_l * n_r + n_r * n_r) / 6.0L);
            const double old_phi_average =
                fields.phi[static_cast<size_t>(as1b_ng + ix)] +
                grid_.dx *
                (fields.Ex_face[static_cast<size_t>(ix) + 1] -
                 fields.Ex_face[static_cast<size_t>(ix)]) / 12.0;
            const double new_phi_average =
                candidate_fields.phi[
                    static_cast<size_t>(as1b_ng + ix)] +
                grid_.dx *
                (candidate_fields.Ex_face[
                     static_cast<size_t>(ix) + 1] -
                 candidate_fields.Ex_face[
                     static_cast<size_t>(ix)]) / 12.0;
            const double rho_delta_cell =
                candidate_fields.rho[
                    static_cast<size_t>(as1b_ng + ix)] -
                fields.rho[static_cast<size_t>(as1b_ng + ix)];
            as1b_term2 += std::fabs(
                0.5 * (old_phi_average + new_phi_average) *
                rho_delta_cell * grid_.dx);
        }
        double as1b_local[3] = {
            static_cast<double>(as1b_term0),
            static_cast<double>(as1b_term1),
            static_cast<double>(as1b_term2)};
        double as1b_global[3] = {0.0, 0.0, 0.0};
        MPI_Allreduce(as1b_local, as1b_global, 3, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.joint_midpoint_poisson_term_abs_sum_energy_before =
            as1b_global[0];
        result.joint_midpoint_poisson_term_abs_sum_energy_after =
            as1b_global[1];
        result.joint_midpoint_poisson_term_abs_sum_potential_charge =
            as1b_global[2];
    }
    const int a1_ng = grid_.nghost;
    std::vector<long double> a1_mismatch(
        static_cast<size_t>(grid_.nx_local), 0.0L);
    // Stage-A-S0 per-cell storage (section 7A): assembly/transport values,
    // per-cell assembly bound and transport scale.
    std::vector<long double> as0_r_assembly(
        static_cast<size_t>(grid_.nx_local), 0.0L);
    std::vector<long double> as0_tau_assembly(
        static_cast<size_t>(grid_.nx_local), 0.0L);
    std::vector<long double> as0_r_transport(
        static_cast<size_t>(grid_.nx_local), 0.0L);
    // gamma_m per section 7A.4: m = N_u*N_perp + 8.  If the denominator is
    // non-positive or non-finite the diagnostic bounds fail (infinite), so
    // every downstream classification gate reports not-closed.
    const long double as0_gamma_m_denom =
        1.0L - static_cast<long double>(nq + 8) *
        static_cast<long double>(std::numeric_limits<double>::epsilon());
    const bool as0_gamma_ok = as0_gamma_m_denom > 0.0L &&
        std::isfinite(static_cast<double>(as0_gamma_m_denom));
    const long double as0_gamma_m = as0_gamma_ok
        ? static_cast<long double>(nq + 8) *
          static_cast<long double>(
              std::numeric_limits<double>::epsilon()) /
          as0_gamma_m_denom
        : std::numeric_limits<long double>::infinity();
    long double a1_local_cont_linf = 0.0L;
    long double a1_local_cont_l1 = 0.0L;
    long double a1_local_rq_linf = 0.0L;
    long double a1_local_mismatch_linf = 0.0L;
    long double a1_local_ubnd_linf = 0.0L;
    long double a1_local_max_drho_dx = 0.0L;
    long double a1_local_max_dt_dj = 0.0L;
    long double a1_local_wc = 0.0L;
    long double a1_local_abs_phi_rc = 0.0L;
    // Stage-A-S0 local accumulators (section 7A.3/7A.4).
    long double as0_local_asm_linf = 0.0L;
    long double as0_local_asm_l1 = 0.0L;
    long double as0_local_transp_linf = 0.0L;
    long double as0_local_transp_proj_linf = 0.0L;
    long double as0_local_mass_delta_linf = 0.0L;
    long double as0_local_parent_max = 0.0L;
    long double as0_local_max_s_transport = 0.0L;
    long double as0_local_w_assembly = 0.0L;
    long double as0_local_w_transport = 0.0L;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        const size_t base = static_cast<size_t>(ix) * nq;
        long double sum_residual = 0.0L;
        // Section 7A.3: per-velocity-cell differencing first, then a single
        // long double accumulation; never two large totals subtracted.
        long double delta_number = 0.0L;
        long double sum_abs_old = 0.0L;
        long double sum_abs_new = 0.0L;
        long double sum_abs_delta_mass = 0.0L;
        for (int q = 0; q < nq; ++q) {
            const long double old_m = static_cast<long double>(
                m_old[base + static_cast<size_t>(q)]);
            const long double new_m = static_cast<long double>(
                candidate[base + static_cast<size_t>(q)]);
            sum_residual += static_cast<long double>(
                accepted_residual[base + static_cast<size_t>(q)]);
            delta_number += new_m - old_m;
            sum_abs_old += std::fabs(old_m);
            sum_abs_new += std::fabs(new_m);
            sum_abs_delta_mass += std::fabs(new_m - old_m);
        }
        const long double delta_q_mass =
            -static_cast<long double>(Const::qe) * delta_number;
        const long double rQ =
            -static_cast<long double>(Const::qe) * sum_residual;
        const size_t rho_id = static_cast<size_t>(a1_ng + ix);
        const long double delta_rho_dx =
            static_cast<long double>(
                candidate_fields.rho[rho_id] - fields.rho[rho_id]) *
            static_cast<long double>(grid_.dx);
        const long double current_div_dt =
            static_cast<long double>(dt) *
            static_cast<long double>(
                accepted_bundle.charge_current_face[
                    static_cast<size_t>(ix) + 1] -
                accepted_bundle.charge_current_face[
                    static_cast<size_t>(ix)]);
        const long double rC = delta_rho_dx + current_div_dt;
        long double u_boundary_sum = 0.0L;
        for (int k = 0; k < nmu; ++k) {
            const size_t bottom =
                static_cast<size_t>(ix) *
                static_cast<size_t>(nupar + 1) *
                static_cast<size_t>(nmu) +
                static_cast<size_t>(k);
            const size_t top = bottom +
                static_cast<size_t>(nupar) *
                static_cast<size_t>(nmu);
            u_boundary_sum += static_cast<long double>(
                accepted_bundle.u_flux_rate[top] -
                accepted_bundle.u_flux_rate[bottom]);
        }
        const long double r_ub =
            -static_cast<long double>(Const::qe) *
            static_cast<long double>(dt) * u_boundary_sum;
        // Same cell-average potential as evaluate_work_identity() and
        // build_potential_pairing_field(): section 2.3 verbatim.
        const double old_phi_average =
            fields.phi[static_cast<size_t>(a1_ng + ix)] +
            grid_.dx *
            (fields.Ex_face[static_cast<size_t>(ix) + 1] -
             fields.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const double new_phi_average =
            candidate_fields.phi[static_cast<size_t>(a1_ng + ix)] +
            grid_.dx *
            (candidate_fields.Ex_face[static_cast<size_t>(ix) + 1] -
             candidate_fields.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const long double phi_bar = static_cast<long double>(
            0.5 * (old_phi_average + new_phi_average));
        const long double mismatch = rC - rQ + r_ub;
        a1_mismatch[static_cast<size_t>(ix)] = mismatch;
        a1_local_cont_linf = std::max(a1_local_cont_linf,
                                      std::fabs(rC));
        a1_local_cont_l1 += std::fabs(rC);
        a1_local_rq_linf = std::max(a1_local_rq_linf, std::fabs(rQ));
        a1_local_mismatch_linf = std::max(a1_local_mismatch_linf,
                                          std::fabs(mismatch));
        a1_local_ubnd_linf = std::max(a1_local_ubnd_linf,
                                      std::fabs(r_ub));
        a1_local_max_drho_dx = std::max(a1_local_max_drho_dx,
                                        std::fabs(delta_rho_dx));
        a1_local_max_dt_dj = std::max(a1_local_max_dt_dj,
                                      std::fabs(current_div_dt));
        const long double phi_rc = phi_bar * rC;
        a1_local_wc += phi_rc;
        a1_local_abs_phi_rc += std::fabs(phi_rc);
        // Stage-A-S0 section 7A.3/7A.4 per-cell quantities.
        const long double r_assembly = delta_rho_dx - delta_q_mass;
        const long double r_transport = delta_q_mass + current_div_dt;
        const long double parent_charge =
            std::fabs(static_cast<long double>(Const::qe)) *
                (sum_abs_old + sum_abs_new) +
            std::max(
                std::fabs(static_cast<long double>(
                    fields.rho[rho_id])),
                std::fabs(static_cast<long double>(
                    candidate_fields.rho[rho_id]))) *
                static_cast<long double>(grid_.dx);
        const long double tau_i_assembly =
            32.0L * as0_gamma_m * std::max(1.0L, parent_charge);
        const long double abs_qe =
            std::fabs(static_cast<long double>(Const::qe));
        const long double s_transport =
            abs_qe * sum_abs_delta_mass +
            std::fabs(static_cast<long double>(dt) *
                      static_cast<long double>(
                          accepted_bundle.charge_current_face[
                              static_cast<size_t>(ix) + 1])) +
            std::fabs(static_cast<long double>(dt) *
                      static_cast<long double>(
                          accepted_bundle.charge_current_face[
                              static_cast<size_t>(ix)]));
        const long double transport_projection =
            r_transport - (rQ - r_ub);
        as0_r_assembly[static_cast<size_t>(ix)] = r_assembly;
        as0_tau_assembly[static_cast<size_t>(ix)] = tau_i_assembly;
        as0_r_transport[static_cast<size_t>(ix)] = r_transport;
        as0_local_asm_linf = std::max(as0_local_asm_linf,
                                      std::fabs(r_assembly));
        as0_local_asm_l1 += std::fabs(r_assembly);
        as0_local_transp_linf = std::max(as0_local_transp_linf,
                                         std::fabs(r_transport));
        as0_local_transp_proj_linf =
            std::max(as0_local_transp_proj_linf,
                     std::fabs(transport_projection));
        as0_local_mass_delta_linf = std::max(as0_local_mass_delta_linf,
                                             std::fabs(delta_q_mass));
        as0_local_parent_max = std::max(as0_local_parent_max,
                                        parent_charge);
        as0_local_max_s_transport = std::max(as0_local_max_s_transport,
                                             s_transport);
        const long double phi_asm = phi_bar * r_assembly;
        const long double phi_transp = phi_bar * r_transport;
        as0_local_w_assembly += phi_asm;
        as0_local_w_transport += phi_transp;
    }
    double a1_max_local[6] = {
        static_cast<double>(a1_local_cont_linf),
        static_cast<double>(a1_local_rq_linf),
        static_cast<double>(a1_local_mismatch_linf),
        static_cast<double>(a1_local_ubnd_linf),
        static_cast<double>(a1_local_max_drho_dx),
        static_cast<double>(a1_local_max_dt_dj)};
    double a1_max_global[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MPI_Allreduce(a1_max_local, a1_max_global, 6, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    double a1_sum_local[2] = {
        static_cast<double>(a1_local_wc),
        static_cast<double>(a1_local_abs_phi_rc)};
    double a1_sum_global[2] = {0.0, 0.0};
    MPI_Allreduce(a1_sum_local, a1_sum_global, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    result.joint_midpoint_continuity_charge_linf = a1_max_global[0];
    result.joint_midpoint_residual_charge_linf = a1_max_global[1];
    result.joint_midpoint_charge_projection_mismatch_linf =
        a1_max_global[2];
    result.joint_midpoint_u_boundary_charge_linf = a1_max_global[3];
    const double tau_c = 8192.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(a1_max_global[4],
                               std::max(a1_max_global[5],
                                        a1_max_global[0])));
    result.joint_midpoint_continuity_roundoff_bound = tau_c;
    // l1 uses its own SUM reduction (it was not part of a1_sum).
    {
        const double a1_l1_local =
            static_cast<double>(a1_local_cont_l1);
        double a1_l1_global = 0.0;
        MPI_Allreduce(&a1_l1_local, &a1_l1_global, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.joint_midpoint_continuity_charge_l1 = a1_l1_global;
    }
    result.joint_midpoint_potential_weighted_continuity_defect =
        a1_sum_global[0];
    const double r_pj_actual =
        result.joint_midpoint_poisson_transport_residual;
    const double predicted =
        poisson_work.residual + a1_sum_global[0];
    result.joint_midpoint_poisson_current_predicted_residual = predicted;
    result.joint_midpoint_poisson_current_prediction_error =
        r_pj_actual - predicted;
    const double tau_a = 8192.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(r_pj_actual),
                               std::max(std::fabs(poisson_work.residual),
                                        std::max(std::fabs(a1_sum_global[0]),
                                                 a1_sum_global[1]))));
    result.joint_midpoint_prediction_roundoff_bound = tau_a;
    unsigned long long a1_local_first_bad =
        std::numeric_limits<unsigned long long>::max();
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        if (std::fabs(static_cast<double>(
                a1_mismatch[static_cast<size_t>(ix)])) > tau_c) {
            a1_local_first_bad =
                static_cast<unsigned long long>(grid_.ix_start + ix);
            break;
        }
    }
    unsigned long long a1_global_first_bad = 0;
    MPI_Allreduce(&a1_local_first_bad, &a1_global_first_bad, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    result.joint_midpoint_continuity_first_bad_global_ix =
        a1_global_first_bad ==
            std::numeric_limits<unsigned long long>::max()
        ? -1
        : static_cast<int>(a1_global_first_bad);
    // Stage-A-S0 reductions and bounds (section 7A.4).  All ranks enter the
    // same collectives in the same order.
    double as0_max_local[5] = {
        static_cast<double>(as0_local_asm_linf),
        static_cast<double>(as0_local_transp_linf),
        static_cast<double>(as0_local_transp_proj_linf),
        static_cast<double>(as0_local_mass_delta_linf),
        static_cast<double>(as0_local_parent_max)};
    double as0_max_global[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    MPI_Allreduce(as0_max_local, as0_max_global, 5, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    double as0_sum_local[2] = {
        static_cast<double>(as0_local_asm_l1),
        static_cast<double>(as0_local_w_assembly)};
    double as0_sum_global[2] = {0.0, 0.0};
    MPI_Allreduce(as0_sum_local, as0_sum_global, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double as0_w_transport_local =
        static_cast<double>(as0_local_w_transport);
    double as0_w_transport_global = 0.0;
    MPI_Allreduce(&as0_w_transport_local, &as0_w_transport_global, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    result.joint_midpoint_density_assembly_mismatch_linf =
        as0_max_global[0];
    {
        const double as0_l1_local = as0_sum_global[0];
        double as0_l1_global = 0.0;
        MPI_Allreduce(&as0_l1_local, &as0_l1_global, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        result.joint_midpoint_density_assembly_mismatch_l1 =
            as0_l1_global;
    }
    result.joint_midpoint_mass_transport_charge_linf =
        as0_max_global[1];
    result.joint_midpoint_transport_projection_mismatch_linf =
        as0_max_global[2];
    result.joint_midpoint_mass_delta_charge_linf = as0_max_global[3];
    result.joint_midpoint_parent_charge_scale_max = as0_max_global[4];
    const double as0_tau_assembly_global = as0_gamma_ok
        ? 32.0 * static_cast<double>(as0_gamma_m) *
              std::max(1.0, as0_max_global[4])
        : std::numeric_limits<double>::infinity();
    result.joint_midpoint_density_assembly_roundoff_bound =
        as0_tau_assembly_global;
    // max_i S_i^transport is a per-cell maximum: reduce it too.
    {
        const double as0_strans_local =
            static_cast<double>(as0_local_max_s_transport);
        double as0_strans_global = 0.0;
        MPI_Allreduce(&as0_strans_local, &as0_strans_global, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        result.joint_midpoint_mass_transport_roundoff_bound =
            8192.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, as0_strans_global);
    }
    result.joint_midpoint_potential_weighted_assembly_defect =
        as0_sum_global[1];
    result.joint_midpoint_potential_weighted_transport_defect =
        as0_w_transport_global;
    result.joint_midpoint_weighted_defect_reconstruction_error =
        std::fabs(result.joint_midpoint_potential_weighted_continuity_defect -
                  result.joint_midpoint_potential_weighted_assembly_defect -
                  result.joint_midpoint_potential_weighted_transport_defect);
    if (as0_tau_assembly_global <=
            result.joint_midpoint_density_assembly_mismatch_linf) {
        for (int ix = 0; ix < grid_.nx_local; ++ix) {
            if (std::fabs(static_cast<double>(
                    as0_r_assembly[static_cast<size_t>(ix)])) >
                static_cast<double>(as0_tau_assembly[
                    static_cast<size_t>(ix)])) {
                result.joint_midpoint_density_assembly_first_bad_global_ix =
                    grid_.ix_start + ix;
                break;
            }
        }
    }
    {
        unsigned long long as0_local_first_transp =
            std::numeric_limits<unsigned long long>::max();
        const double tau_t =
            result.joint_midpoint_mass_transport_roundoff_bound;
        for (int ix = 0; ix < grid_.nx_local; ++ix) {
            if (std::fabs(static_cast<double>(
                    as0_r_transport[static_cast<size_t>(ix)])) > tau_t) {
                as0_local_first_transp =
                    static_cast<unsigned long long>(grid_.ix_start + ix);
                break;
            }
        }
        unsigned long long as0_global_first_transp = 0;
        MPI_Allreduce(&as0_local_first_transp, &as0_global_first_transp, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        result.joint_midpoint_mass_transport_first_bad_global_ix =
            as0_global_first_transp ==
                std::numeric_limits<unsigned long long>::max()
            ? -1
            : static_cast<int>(as0_global_first_transp);
    }
    // Stage-A-S1 (section 7B.4): read-only comparison between the accepted
    // incremental charge assembly and the legacy absolute form for the SAME
    // final candidate.  The absolute form never overwrites the incremental
    // result.
    {
        long double as1_local_incremental_linf = 0.0L;
        long double as1_local_absolute_linf = 0.0L;
        long double as1_local_diff_linf = 0.0L;
        long double as1_local_bound_max = 0.0L;
        const long double as1_gamma_m =
            static_cast<long double>(nq + 8) *
            static_cast<long double>(
                std::numeric_limits<double>::epsilon()) /
            as0_gamma_m_denom;
        for (int ix = 0; ix < grid_.nx_local; ++ix) {
            const size_t rho_id2 =
                static_cast<size_t>(a1_ng + ix);
            const size_t base2 = static_cast<size_t>(ix) * nq;
            const double rho_incremental =
                candidate_fields.rho[rho_id2];
            double number = 0.0;
            long double sum_abs_state = 0.0L;
            long double sum_abs_old2 = 0.0L;
            for (int q = 0; q < nq; ++q) {
                number += candidate[base2 + static_cast<size_t>(q)];
                sum_abs_state += std::fabs(static_cast<long double>(
                    candidate[base2 + static_cast<size_t>(q)]));
                sum_abs_old2 += std::fabs(static_cast<long double>(
                    m_old[base2 + static_cast<size_t>(q)]));
            }
            const double ni = ix < static_cast<int>(ion_density.size())
                ? ion_density[static_cast<size_t>(ix)]
                : 0.0;
            const double rho_absolute =
                Const::qe * (ni - number / grid_.dx);
            const double form_difference =
                rho_incremental - rho_absolute;
            const long double parent_scale =
                std::max(std::fabs(static_cast<long double>(Const::qe) *
                                   static_cast<long double>(ni)),
                         static_cast<long double>(Const::qe) *
                             (sum_abs_state + sum_abs_old2) /
                             static_cast<long double>(grid_.dx)) +
                std::fabs(static_cast<long double>(
                    fields.rho[rho_id2]));
            const long double bound_i =
                as0_gamma_ok
                ? 32.0L * as1_gamma_m * std::max(1.0L, parent_scale)
                : std::numeric_limits<long double>::infinity();
            as1_local_incremental_linf =
                std::max(as1_local_incremental_linf,
                         std::fabs(static_cast<long double>(
                             rho_incremental)));
            as1_local_absolute_linf =
                std::max(as1_local_absolute_linf,
                         std::fabs(static_cast<long double>(
                             rho_absolute)));
            as1_local_diff_linf = std::max(as1_local_diff_linf,
                std::fabs(static_cast<long double>(form_difference)));
            as1_local_bound_max = std::max(as1_local_bound_max,
                                           std::fabs(bound_i));
        }
        double as1_local[4] = {
            static_cast<double>(as1_local_incremental_linf),
            static_cast<double>(as1_local_absolute_linf),
            static_cast<double>(as1_local_diff_linf),
            static_cast<double>(as1_local_bound_max)};
        double as1_global[4] = {0.0, 0.0, 0.0, 0.0};
        MPI_Allreduce(as1_local, as1_global, 4, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        result.joint_midpoint_candidate_rho_incremental = as1_global[0];
        result.joint_midpoint_candidate_rho_absolute = as1_global[1];
        result.joint_midpoint_candidate_rho_form_difference =
            as1_global[2];
        result.joint_midpoint_candidate_rho_form_roundoff_bound =
            as1_global[3];
    }
    const double combined_energy_residual = delta_ke_u +
        poisson_work.field_energy_change - poisson_work.electrode_work;
    const double candidate_energy_scale = std::max(
        1.0e-300, std::max(std::fabs(delta_ke_u),
                            std::max(std::fabs(poisson_work.field_energy_change),
                                     std::fabs(poisson_work.electrode_work))));
    if (!poisson_work.finite ||
        std::fabs(combined_energy_residual) / candidate_energy_scale > 1.0e-8) {
        result.failure_code = 75;
        result.failure_stage = "joint_midpoint_energy_residual";
        return result;
    }
    std::copy(candidate.begin(), candidate.end(),
              electrons.f.begin() +
                  static_cast<std::ptrdiff_t>(grid_.nghost * nq));
    electrons.compute_moments();
    fields = candidate_fields;
    const double number_after = global_sum(electrons.total_particle_number());
    const double field_after = field_energy(fields);
    const double kinetic_after = global_sum(electrons.total_kinetic_energy());
    result.accepted = true;
    result.gauss_ok = accepted_poisson_residual <= poisson_tolerance;
    result.ledger.background_number_before = number_before;
    result.ledger.background_number_after = number_after;
    result.ledger.field_energy_before = field_before;
    result.ledger.field_energy = field_after;
    result.ledger.background_kinetic_energy_before = kinetic_before;
    result.ledger.background_kinetic_energy = kinetic_after;
    result.ledger.domain_energy_before = field_before + kinetic_before;
    result.ledger.domain_energy_after = field_after + kinetic_after;
    result.ledger.domain_energy_change =
        result.ledger.domain_energy_after - result.ledger.domain_energy_before;
    result.ledger.gauss_charge_residual =
        field_solver_.diagnostics().boundary_charge_residual;
    ++step_count_;
    (void)time;
    (void)beam;
    return result;
}

VpfpStepResult VpfpIntegrator::advance(Species& electrons, BeamPIC& beam,
                                       EMFields& fields,
                                       const std::vector<double>& ion_density,
                                       double time, double dt, int mpi_rank,
                                       int mpi_size)
{
    if (background_phase_space_mode_ ==
        BackgroundPhaseSpaceMode::JOINT_MIDPOINT_ENERGY) {
        return advance_joint_midpoint(electrons, beam, fields, ion_density,
                                      time, dt, mpi_rank, mpi_size);
    }
    // JC1 (section 4.7 item 7): mode dispatch.  Legacy continues to the
    // existing advance_background/advance_with_beam unchanged.
    if (field_particle_coupling_config_.mode ==
        FieldParticleCouplingMode::DiscreteGradient) {
        return advance_discrete_gradient(electrons, beam, fields,
                                         ion_density, time, dt,
                                         mpi_rank, mpi_size);
    }
    if (beam_enabled_) {
        return advance_with_beam(electrons, beam, fields, ion_density, time,
                                 dt, mpi_rank, mpi_size);
    }
    return advance_background(electrons, beam, fields, ion_density, time, dt,
                              mpi_rank, mpi_size);
}

VpfpStepResult VpfpIntegrator::advance_fixed_state_field_particle_audit(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density, double time, double dt,
    int mpi_rank, int mpi_size)
{
    // The caller creates a dedicated integrator and passes checkpoint copies.
    // This branch may therefore commit its own trial state between dt/2 and
    // dt/4 substeps, while no production integrator or checkpoint object is
    // touched.  Restoring tail_accepted_ here would incorrectly restart the
    // Tail pusher at every audit substep.
    const bool audit_mode_saved = fixed_state_field_particle_audit_mode_;

    fixed_state_field_particle_audit_mode_ = true;
    const VpfpStepResult result = advance(
        electrons, beam, fields, ion_density, time, dt, mpi_rank, mpi_size);

    fixed_state_field_particle_audit_mode_ = audit_mode_saved;
    return result;
}

// JC1 (section 4.7 item 8): stub that returns a clear "not implemented"
// failure code.  Must NOT call any incomplete new algorithm.
VpfpStepResult VpfpIntegrator::advance_discrete_gradient(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density, double time, double dt,
    int mpi_rank, int mpi_size)
{
    VpfpStepResult result = {};
    result.stage_energy_audit_enabled = stage_energy_audit_enabled_;
    result.stage_energy_audit_valid = true;
    result.stage_energy_count = 0;
    result.tail_return_enabled =
        background_tail_enabled_ && tail_bulk_return_.config().enabled;
    result.split_used = false;
    result.finite = true;
    result.cfl_ok = dt > 0.0;
    result.gauss_ok = false;
    result.collision_ok = true;
    result.conversion_ok = true;
    result.tail_ok = true;
    result.failure_code = 0;
    result.failure_stage = "advance_discrete_gradient";
    // Stage-energy audit scratch: the DG path fills the eleven records across
    // prepare (0-3), the accepted trial (4-6) and the post-field transaction
    // (7-10) instead of a single inline pass.  Reset it per step so a previous
    // step's trial records never leak into this step.
    stage_energy_scratch_ = VpfpStepResult();
    stage_energy_scratch_.stage_energy_audit_enabled =
        stage_energy_audit_enabled_;
    stage_energy_scratch_.stage_energy_audit_valid = true;
    stage_energy_scratch_.stage_energy_count = 0;
    stage_sources_ = VpfpStageEnergyRecord();
    stage_sources_frozen_ = VpfpStageEnergyRecord();
    // JC4 (section 7.7): always true for this path; legacy runs set it false.
    result.field_particle_coupling_enabled = true;
    result.post_field_charge_residual_linf = 0.0;

    if (!initialized_ || !result.cfl_ok) {
        result.failure_code = 1;
        result.failure_stage = "final_validation";
        return result;
    }

    // Full transaction: prepare once, solve, then commit atomically.  Any
    // failure restores the accepted step count and leaves every accepted
    // object, RNG and cumulative ledger unchanged.
    const long long step_before = step_count_;
    const bool tail_on = background_tail_enabled_;
    initialize_energy_ledger(electrons, beam, fields, tail_on);
    result.ledger.background_number_before =
        global_sum(electrons.total_particle_number());
    result.ledger.beam_number_before = beam_enabled_
        ? global_sum(beam.total_particle_number(grid_)) : 0.0;
    result.ledger.tail_number_before = tail_on
        ? global_sum(tail_total_weight(tail_accepted_)) : 0.0;
    result.ledger.tail_particle_count_before = tail_on
        ? global_sum_u64(tail_accepted_.particles.size()) : 0;
    result.ledger.combined_number_before =
        result.ledger.background_number_before +
        result.ledger.tail_number_before;
    result.ledger.field_energy_before = accepted_field_energy_;
    result.ledger.background_kinetic_energy_before =
        accepted_background_kinetic_energy_;
    result.ledger.beam_kinetic_energy_before =
        accepted_beam_kinetic_energy_;
    result.ledger.tail_kinetic_energy_before =
        accepted_tail_kinetic_energy_;

    FieldParticleFrozenState frozen;
    VpfpFailureInfo failure;
    if (!prepare_field_particle_frozen_state(
            electrons, beam, fields, ion_density, time, dt,
            mpi_rank, mpi_size, frozen, failure)) {
        step_count_ = step_before;
        result.accepted = false;
        result.failure_code = failure.code;
        result.failure_stage = failure.stage;
        result.failing_rank = failure.failing_rank;
        result.failing_ix = failure.failing_ix;
        return result;
    }
    stage_sources_frozen_ = stage_sources_;

    FieldParticleTrial accepted_trial;
    FieldParticleIterationDiagnostics diagnostics;
    if (!solve_field_particle_center(frozen, accepted_trial, diagnostics,
                                     failure)) {
        step_count_ = step_before;
        result.accepted = false;
        result.failure_code = failure.code;
        result.failure_stage = failure.stage;
        result.failing_rank = failure.failing_rank;
        result.failing_ix = failure.failing_ix;
        // §6.12: expose the Picard diagnostics on the failure path too so
        // the .result files show how many iterations ran and the residual
        // values that blocked convergence.
        result.field_particle_converged = diagnostics.converged;
        result.field_particle_iterations = diagnostics.iterations;
        result.field_particle_trial_evaluations = diagnostics.iterations;
        result.field_particle_relaxation = diagnostics.relaxation;
        result.field_particle_residual_l2 = diagnostics.field_residual_l2;
        result.field_particle_residual_linf = diagnostics.field_residual_linf;
        result.field_particle_pairing_residual = diagnostics.pairing_residual;
        result.post_field_charge_residual_linf = post_field_charge_residual_linf_;
        return result;
    }

    // Populate the accepted-step diagnostics from the exact converged trial.
    // This is bookkeeping only: no operator is replayed and no state is
    // altered before the atomic swaps below.
    const VlasovStepDiagnostics& vd = accepted_trial.vlasov_diagnostics;
    result.ledger.background_left_flux =
        (frozen.x1_diagnostics.x_first.left_inflow_number -
         frozen.x1_diagnostics.x_first.left_outflow_number) +
        (vd.x_second.left_inflow_number - vd.x_second.left_outflow_number);
    result.ledger.background_right_flux =
        (frozen.x1_diagnostics.x_first.right_inflow_number -
         frozen.x1_diagnostics.x_first.right_outflow_number) +
        (vd.x_second.right_inflow_number - vd.x_second.right_outflow_number);
    result.ledger.background_left_inflow_energy =
        frozen.x1_diagnostics.x_first.left_inflow_energy +
        vd.x_second.left_inflow_energy;
    result.ledger.background_left_outflow_energy =
        frozen.x1_diagnostics.x_first.left_outflow_energy +
        vd.x_second.left_outflow_energy;
    result.ledger.background_right_inflow_energy =
        frozen.x1_diagnostics.x_first.right_inflow_energy +
        vd.x_second.right_inflow_energy;
    result.ledger.background_right_outflow_energy =
        frozen.x1_diagnostics.x_first.right_outflow_energy +
        vd.x_second.right_outflow_energy;
    result.ledger.background_tail_number_loss = vd.u_full.tail_number_loss;
    result.ledger.background_tail_energy_loss = vd.u_full.tail_energy_loss;

    const double local_conversion_number =
        accepted_trial.conversion_number_removed +
        accepted_trial.post_conversion_number_removed;
    const double local_conversion_energy =
        accepted_trial.conversion_energy_removed +
        accepted_trial.post_conversion_energy_removed;
    const std::uint64_t local_conversion_particles =
        accepted_trial.conversion_particles_created +
        accepted_trial.post_conversion_particles_created;
    result.ledger.conversion_number_removed = local_conversion_number;
    result.ledger.conversion_energy_removed = local_conversion_energy;
    result.ledger.conversion_particles_created = local_conversion_particles;
    globalize_conversion_ledger(result.ledger);

    result.tail_return = accepted_trial.post_tail_return;
    result.ledger.collision_reservoir_energy =
        frozen.first_collision.reservoir_energy_change +
        accepted_trial.second_collision_reservoir_energy;
    result.ledger.electrostatic_boundary_work =
        field_solver_.boundary_energy_work(
            fields, accepted_trial.final_fields_trial, mpi_rank, mpi_size);

    result.ledger.background_number_after =
        global_sum(accepted_trial.bulk_trial.total_particle_number());
    result.ledger.beam_number_after = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.total_particle_number(grid_))
        : 0.0;
    result.ledger.beam_injected = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.last_injected_number()) : 0.0;
    result.ledger.beam_outflow = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.last_outflow_number()) : 0.0;
    result.ledger.beam_injected_energy = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.last_injected_energy()) : 0.0;
    result.ledger.beam_outflow_energy = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.last_outflow_energy()) : 0.0;
    result.ledger.tail_number_after = frozen.tail_on
        ? global_sum(tail_total_weight(accepted_trial.tail_trial)) : 0.0;
    result.ledger.tail_particle_count_after = frozen.tail_on
        ? global_sum_u64(accepted_trial.tail_trial.particles.size()) : 0;
    result.ledger.tail_kinetic_energy_after = frozen.tail_on
        ? global_sum(tail_total_kinetic_energy(accepted_trial.tail_trial))
        : 0.0;
    result.ledger.tail_outflow_number = frozen.tail_on
        ? global_sum(
            accepted_trial.tail_trial.outflow_ledger().left_number +
            accepted_trial.tail_trial.outflow_ledger().right_number) : 0.0;
    result.ledger.tail_outflow_energy = frozen.tail_on
        ? global_sum(
            accepted_trial.tail_trial.outflow_ledger().left_kinetic_energy +
            accepted_trial.tail_trial.outflow_ledger().right_kinetic_energy)
        : 0.0;
    result.ledger.combined_number_after =
        result.ledger.background_number_after +
        result.ledger.tail_number_after;
    if (frozen.tail_on) {
        const double expected_tail_after =
            result.ledger.tail_number_before +
            result.ledger.conversion_number_removed -
            result.ledger.tail_outflow_number - result.tail_return.number;
        result.ledger.tail_number_balance_error = std::fabs(
            result.ledger.tail_number_after - expected_tail_after) /
            std::max(1.0, result.ledger.tail_number_before +
                          result.ledger.conversion_number_removed +
                          result.ledger.tail_outflow_number);
    }
    const double x_before = frozen.x1_diagnostics.x_first.number_before;
    const double x_in = frozen.x1_diagnostics.x_first.inflow_number +
                        vd.x_second.inflow_number;
    const double x_out = frozen.x1_diagnostics.x_first.outflow_number +
                         vd.x_second.outflow_number;
    const double expected_bulk_after =
        x_before + x_in - x_out - vd.u_full.tail_number_loss -
        result.ledger.conversion_number_removed + result.tail_return.number;
    result.ledger.remap_ledger_residual = std::fabs(
        result.ledger.background_number_after - expected_bulk_after) /
        std::max(1.0, x_before);
    result.ledger.field_energy = field_energy(accepted_trial.final_fields_trial);
    // §16 条件3 核查：收敛 trial 实际用于推粒子的 E_pair 场（G_P(Phi^n,Phi^{n+1}),
    // 即 trial_force_fields_ 在收敛时的值）的全局场能。
    result.ledger.e_pair_field_energy = field_energy(trial_force_fields_);
    result.ledger.background_kinetic_energy =
        global_sum(accepted_trial.bulk_trial.total_kinetic_energy());
    result.ledger.beam_kinetic_energy = frozen.beam_on
        ? global_sum(accepted_trial.beam_trial.total_kinetic_energy()) : 0.0;
    result.ledger.gauss_charge_residual =
        field_solver_.diagnostics().boundary_charge_residual;
    result.gauss_ok = accepted_trial.final_poisson_pass;
    result.ledger.fct_energy_change = 0.0;
    finalize_energy_ledger(result);

    // §15.13.4: pass per-cell conversion diagnostics from the accepted DG
    // trial and the C2 post-field transaction to Gate I.  This is audit-only;
    // no physical state is modified.
    result.conversion_events = accepted_trial.conversion_events;
    result.conversion_events.insert(result.conversion_events.end(),
        accepted_trial.post_conversion_events.begin(),
        accepted_trial.post_conversion_events.end());

    // Gate I: the discrete-gradient branch must emit the same accepted-state
    // field-particle pairing ledger as the legacy branches.  Use the exact
    // converged trial before its transactional swap; no operator is replayed
    // and no field/current used by production is modified.
    finalize_field_particle_power_audit(
        result, electrons, accepted_trial.bulk_trial, beam,
        accepted_trial.beam_trial,
        frozen.tail_on ? &tail_accepted_ : NULL,
        frozen.tail_on ? &accepted_trial.tail_trial : NULL,
        fields, accepted_trial.final_fields_trial,
        accepted_trial.pairing_current_work,
        accepted_trial.tail_force_work, accepted_trial.beam_force_work,
        dt, mpi_rank, mpi_size);

    // Stage-energy audit: the accepted trial already left the exact
    // force/conversion/x2 records (stages 4-6) and the post-field transaction
    // recorded C2/conversion-after-collision/return (stages 7-9).  The final
    // Poisson (stage 10) reuses the converged trial fields after the
    // charge-invariant C2/return, matching the legacy FINAL_POISSON capture.
    if (stage_energy_audit_enabled_) {
        if (mpi_rank == 0) {
            stage_sources_.electrostatic_boundary_work =
                result.ledger.electrostatic_boundary_work;
        }
        capture_dg_stage(VPFP_STAGE_FINAL_POISSON, accepted_trial.bulk_trial,
                         frozen.tail_on ? &accepted_trial.tail_trial : NULL,
                         frozen.beam_on ? &accepted_trial.beam_trial : NULL,
                         accepted_trial.final_fields_trial);
        result.stage_energy_count = stage_energy_scratch_.stage_energy_count;
        result.stage_energy_audit_valid =
            stage_energy_scratch_.stage_energy_audit_valid;
        for (int i = 0; i < result.stage_energy_count; ++i) {
            result.stage_energy[i] = stage_energy_scratch_.stage_energy[i];
        }
        finalize_stage_energy_audit(result);
    }

    // JC3 (section 6.9): atomic commit in the documented order.
    electrons.swap_state(accepted_trial.bulk_trial);
    if (frozen.tail_on) tail_accepted_.swap_state(accepted_trial.tail_trial);
    if (frozen.beam_on) beam.swap_state(accepted_trial.beam_trial);
    swap_emfields(fields, accepted_trial.final_fields_trial);
    if (frozen.beam_on) beam.commit_injection_schedule(frozen.beam_schedule,
                                                       mpi_rank);

    // Accepted cumulative ledgers and combined checksums.
    tail_cumulative_.conversion_number += local_conversion_number;
    tail_cumulative_.conversion_energy += local_conversion_energy;
    tail_cumulative_.particles_created += local_conversion_particles;
    tail_cumulative_.return_number += result.tail_return.number;
    tail_cumulative_.return_px += result.tail_return.px;
    tail_cumulative_.return_jx_dx += result.tail_return.jx_dx;
    tail_cumulative_.return_energy += result.tail_return.energy;
    tail_cumulative_.return_pixx_dx += result.tail_return.pixx_dx;
    tail_cumulative_.return_piperp_dx += result.tail_return.piperp_dx;
    tail_cumulative_.return_particles_removed +=
        result.tail_return.particles_removed;
    tail_cumulative_.return_deferred_groups +=
        result.tail_return.deferred_infeasible_groups +
        result.tail_return.deferred_rank_boundary_groups;
    if (frozen.tail_on) {
        tail_cumulative_.outflow_number += static_cast<std::uint64_t>(
            tail_accepted_.outflow_ledger().left_number +
            tail_accepted_.outflow_ledger().right_number);
    }
    combined_checksum_.number = result.ledger.combined_number_after;
    combined_checksum_.kinetic_energy =
        result.ledger.background_kinetic_energy +
        result.ledger.tail_kinetic_energy_after;
    combined_checksum_.field_energy = result.ledger.field_energy;
    result.tail_particles_local_max = frozen.tail_on
        ? global_max_u64(tail_accepted_.particles.size()) : 0;

    // step_count_ stays at frozen.candidate_step (accepted step advanced by
    // prepare and retained by the successful commit).
    result.accepted = true;
    result.failure_code = 0;
    result.failure_stage = "final_validation";
    result.collision_ok = true;
    result.conversion_ok = true;
    result.tail_ok = true;
    result.finite = finite_species(electrons);
    // JC3 (section 6.12): accepted-step Picard diagnostics.
    result.field_particle_converged = diagnostics.converged;
    result.field_particle_iterations = diagnostics.iterations;
    result.field_particle_trial_evaluations = diagnostics.iterations;
    result.field_particle_relaxation = diagnostics.relaxation;
    result.field_particle_residual_l2 = diagnostics.field_residual_l2;
    result.field_particle_residual_linf = diagnostics.field_residual_linf;
    result.field_particle_pairing_residual = diagnostics.pairing_residual;
    result.post_field_charge_residual_linf = post_field_charge_residual_linf_;
    return result;
}

// JC2 (section 5.5.1): freeze the state for one physical step.
// Extracts blocks 1-10 from advance_with_beam() per §5.5.1.
bool VpfpIntegrator::prepare_field_particle_frozen_state(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density,
    double time, double dt, int mpi_rank, int mpi_size,
    FieldParticleFrozenState& frozen, VpfpFailureInfo& failure)
{
    failure = VpfpFailureInfo();

    // Gate-I workspace follows the same production x/u/Tail/Beam operators
    // as the discrete-gradient trial.  It is diagnostic-only and is reset
    // once per physical step; Picard trials overwrite only their trial-owned
    // u-force and second-x contributions.
    pairing_workspace_.enabled = field_particle_power_audit_enabled_;
    if (field_particle_power_audit_enabled_) {
        pairing_workspace_.bulk_x1.init(grid_.nx_local);
        pairing_workspace_.bulk_x2.init(grid_.nx_local);
        pairing_workspace_.bulk_number_pre_x1.clear();
        pairing_workspace_.bulk_number_post_x1.clear();
        pairing_workspace_.bulk_number_pre_x2.clear();
        pairing_workspace_.bulk_number_post_x2.clear();
        pairing_workspace_.cell_work.init(grid_.nx_local);
        pairing_workspace_.tail_work_ledger = 0.0;
        pairing_workspace_.beam_work_ledger = 0.0;
    }

    // Block 1: CFL and config gates
    if (!initialized_ || dt <= 0.0) {
        failure.code = 1;
        failure.stage = "frozen_cfl";
        return false;
    }
    const bool tail_on = background_tail_enabled_;
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        !collision_interface_zero_wall_validation_ &&
        !collision_interface_exporting_absorbing_) {
        failure.code = 12;
        failure.stage = "frozen_collision_interface";
        return false;
    }
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        collision_interface_exporting_absorbing_ &&
        collision_.bulk_integrator() !=
            BulkCollisionIntegrator::CHANG_COOPER_FLUX) {
        failure.code = 12;
        failure.stage = "frozen_collision_integrator";
        return false;
    }
    ++step_count_;

    // Block 2: Pre-step number/energy ledger (diagnostic only, not committed)
    // The actual ledger is computed in evaluate_field_particle_trial.

    // Block 3: beam_work_.begin_step() and generate_injection_schedule()
    beam_work_.density = beam.density;
    beam_work_.begin_step(grid_, dt);
    if (beam_enabled_) {
        frozen.beam_schedule =
            beam.generate_injection_schedule(grid_, time, dt, mpi_rank);
    }

    // Block 4: tail_work_ = tail_accepted_ and tail_work_.begin_step()
    if (tail_on) {
        tail_work_ = tail_accepted_;
        tail_work_.begin_step(grid_, dt);
    }
    stage_sources_.collision_reservoir_energy = 0.0;
    stage_sources_.conversion_energy = 0.0;
    stage_sources_.tail_outflow_energy = 0.0;
    stage_sources_.beam_injected_energy = 0.0;
    stage_sources_.beam_outflow_energy = 0.0;
    capture_dg_stage(VPFP_STAGE_ACCEPTED_N, electrons,
                     tail_on ? &tail_work_ : NULL, &beam, fields);

    // Block 5: First collision half-step C1
    Species* collision_input = &electrons;
    if (hybrid_collision_active_ || !collision_.is_trivial()) {
        state_collision_trial_ = electrons;
        collision_input = &state_collision_trial_;
    }
    const bool observe_collision_flux = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        collision_interface_exporting_absorbing_ &&
        !collision_.is_trivial() &&
        collision_.bulk_integrator() == BulkCollisionIntegrator::CHANG_COOPER_FLUX;
    if (!apply_collision_half(
            *collision_input, tail_work_, time, 0.5 * dt, 0, mpi_rank,
            frozen.first_collision, frozen.first_hybrid,
            frozen.first_hybrid_valid,
            observe_collision_flux ? &frozen.first_collision_flux : NULL)) {
        failure.code = 5;
        failure.stage = "frozen_collision_half1";
        return false;
    }
    stage_sources_.collision_reservoir_energy +=
        frozen.first_collision.reservoir_energy_change;
    capture_dg_stage(VPFP_STAGE_COLLISION_HALF1, *collision_input,
                     tail_on ? &tail_work_ : NULL, &beam, fields);

    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            *collision_input, pairing_workspace_.bulk_number_pre_x1);
    }

    // Block 6: vlasov_.first_x_half()
    VlasovStepDiagnostics vlasov_diag;
    if (!vlasov_.first_x_half(*collision_input, state_x_half_, time, 0.5 * dt,
                              vlasov_diag,
                              field_particle_power_audit_enabled_
                                  ? &pairing_workspace_.bulk_x1 : NULL)) {
        failure.code = 2;
        failure.stage = "frozen_first_x_half";
        return false;
    }
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_x_half_, pairing_workspace_.bulk_number_post_x1);
    }
    frozen.x1_diagnostics = vlasov_diag;
    stage_sources_.background_left_inflow_energy =
        vlasov_diag.x_first.left_inflow_energy;
    stage_sources_.background_left_outflow_energy =
        vlasov_diag.x_first.left_outflow_energy;
    stage_sources_.background_right_inflow_energy =
        vlasov_diag.x_first.right_inflow_energy;
    stage_sources_.background_right_outflow_energy =
        vlasov_diag.x_first.right_outflow_energy;

    // Block 7: tail_work_.drift_half() and midpoint density deposit
    if (tail_on) {
        tail_work_.drift_half(grid_, 0.5 * dt, mpi_rank, mpi_size);
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
    }

    // Block 8: beam.predict_to_midpoint() and midpoint density deposit.
    // Beam-off is a real fixture: it must not generate a schedule or advance
    // Beam state merely because the shared frozen path is reused.
    if (beam_enabled_) {
        beam.predict_to_midpoint(frozen.beam_schedule, grid_, fields, time, dt,
                                 mpi_rank, mpi_size, beam_work_);
        beam_work_.deposit_density(grid_, mpi_rank, mpi_size);
    }
    capture_dg_stage(VPFP_STAGE_X_HALF1, state_x_half_,
                     tail_on ? &tail_work_ : NULL,
                     beam_enabled_ ? &beam_work_ : &beam, fields);

    // Block 9: Legacy midpoint Poisson (only for first-round initial value)
    const std::vector<double> empty_tail_density;
    const std::vector<double>& tail_mid_density =
        tail_on ? tail_work_.density : empty_tail_density;
    midpoint_fields_.set_charge_density(state_x_half_, tail_mid_density,
                                        beam_work_.density, ion_density);
    OpenGaussSolveOptions midpoint_options;
    midpoint_options.reconstruct_phi = false;
    midpoint_options.compute_l1 = false;
    midpoint_options.compute_boundary_audit = false;
    field_solver_.solve(midpoint_fields_, mpi_rank, mpi_size, midpoint_options);
    capture_dg_stage(VPFP_STAGE_MIDPOINT_POISSON, state_x_half_,
                     tail_on ? &tail_work_ : NULL,
                     beam_enabled_ ? &beam_work_ : &beam, midpoint_fields_);

    // §5.6: n-layer potential bootstrap.
    // Copy accepted fields to field_n_pairing_, rebuild phi with full audit,
    // verify Ex_face matches accepted state within 1e-12 relative tolerance.
    field_n_pairing_ = fields;
    field_n_pairing_.rho = fields.rho;
    const std::vector<double> ex_face_before = field_n_pairing_.Ex_face;
    OpenGaussSolveOptions bootstrap_options;
    bootstrap_options.reconstruct_phi = true;
    bootstrap_options.compute_l1 = true;
    bootstrap_options.compute_boundary_audit = true;
    field_solver_.solve(field_n_pairing_, mpi_rank, mpi_size, bootstrap_options);
    double local_max_ex_accepted = 0.0;
    double local_max_ex_diff = 0.0;
    for (size_t i = 0; i < ex_face_before.size(); ++i) {
        local_max_ex_accepted = std::max(
            local_max_ex_accepted, std::fabs(ex_face_before[i]));
        local_max_ex_diff = std::max(local_max_ex_diff,
            std::fabs(field_n_pairing_.Ex_face[i] - ex_face_before[i]));
    }
    double local_bootstrap_scale[2] = {
        local_max_ex_accepted, local_max_ex_diff
    };
    double global_bootstrap_scale[2] = { 0.0, 0.0 };
    MPI_Allreduce(local_bootstrap_scale, global_bootstrap_scale, 2,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const double e_floor = 1.0e-30;
    const double r_bootstrap = global_bootstrap_scale[1] /
        std::max(e_floor, global_bootstrap_scale[0]);
    if (r_bootstrap > 1.0e-12) {
        failure.code = 208;
        failure.stage = "accepted_field_poisson_mismatch";
        return false;
    }

    // Block 10: Save frozen state pointers
    frozen.bulk_x_half = &state_x_half_;
    frozen.tail_midpoint = tail_on ? &tail_work_ : NULL;
    frozen.tail_work_snapshot = tail_on ? &tail_work_ : NULL;
    frozen.beam_midpoint = &beam_work_;
    frozen.beam_work_snapshot = &beam_work_;
    frozen.field_n_pairing = &field_n_pairing_;
    frozen.candidate_step = step_count_;
    frozen.time = time;
    frozen.dt = dt;
    frozen.beam_on = beam_enabled_;
    frozen.tail_on = tail_on;
    frozen.ion_density = ion_density;
    frozen.mpi_rank = mpi_rank;
    frozen.mpi_size = mpi_size;
    frozen.bootstrap_residual = r_bootstrap;
    frozen.bootstrap_pass = true;

    return true;
}

// JC2 (section 5.5.2): evaluate one deterministic trial from frozen state.
// §5.8: save/restore shared state for determinism across consecutive calls.
bool VpfpIntegrator::evaluate_field_particle_trial(
    const FieldParticleFrozenState& frozen,
    const std::vector<double>& pairing_field_guess,
    FieldParticleTrial& trial, VpfpFailureInfo& failure)
{
    failure = VpfpFailureInfo();
    const bool tail_on = frozen.tail_on;
    const int mpi_rank = frozen.mpi_rank;
    const int mpi_size = frozen.mpi_size;

    // §5.8: Save shared state for restore after trial (determinism)
    const Species saved_state_u_full = state_u_full_;
    const Species saved_state_np1 = state_np1_;
    const BackgroundTailPIC saved_tail_work = tail_work_;
    const BeamPIC saved_beam_work = beam_work_;

    // Stage-energy audit: each Picard trial re-runs the same sub-steps from
    // the same frozen state.  Reset the member scratch to the four stages
    // already captured by prepare (accepted_n .. midpoint_poisson) so only the
    // converged (accepted) trial's records survive.
    stage_energy_scratch_.stage_energy_count = 4;
    stage_sources_ = stage_sources_frozen_;

    // Copy frozen state to trial buffers (§5.5.2 preamble)
    trial.tail_trial = tail_work_;
    trial.beam_trial = beam_work_;
    trial.second_collision_reservoir_energy = 0.0;
    trial.post_conversion_number_removed = 0.0;
    trial.post_conversion_energy_removed = 0.0;
    trial.post_conversion_particles_created = 0;
    trial.post_tail_return = TailBulkReturnDiagnostics();
    trial.post_field_charge_residual_linf = 0.0;

    // Write pairing_field_guess into trial_force_fields_.Ex_face
    trial_force_fields_.Ex_face = pairing_field_guess;

    // Compute Ex and ghosts from Ex_face
    field_solver_.populate_electric_components_from_faces(
        trial_force_fields_, mpi_rank, mpi_size);
    trial.force_field_face = trial_force_fields_.Ex_face;

    // Block: vlasov_.u_full().  The JC trial must use the same threshold
    // interface contract as the legacy production path: u_full owns the
    // exported face batch and applies the interface sink before that exact
    // batch is consumed by the conversion transaction.
    VlasovStepDiagnostics vlasov_diag;
    std::vector<double>* bulk_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.bulk_delta_ke_cell : NULL;
    const bool observe_upar_flux = tail_on &&
        (tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT ||
         tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE);
    BulkTailFluxBatch exported_flux;
    exported_flux.apply_interface_sink = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE;
    if (!vlasov_.u_full(state_x_half_, state_u_full_, trial_force_fields_,
                        frozen.time + 0.5 * frozen.dt, frozen.dt, vlasov_diag,
                        observe_upar_flux ? partition_ : NULL,
                        observe_upar_flux ? &exported_flux : NULL,
                        tail_flux_quadrature_order_, bulk_cell_work)) {
        failure.code = 2;
        failure.stage = "trial_u_full";
        state_u_full_ = saved_state_u_full;
        state_np1_ = saved_state_np1;
        tail_work_ = saved_tail_work;
        beam_work_ = saved_beam_work;
        return false;
    }
    trial.bulk_force_work =
        state_u_full_.total_kinetic_energy() -
        state_x_half_.total_kinetic_energy();
    // §6.7: independent force-work computation from the production u_parallel
    // remap's Gate-C bulk work (already-global sums).  This is the same force
    // step, so at the self-consistent state it must equal the kinetic change
    // bulk_force_work to roundoff; the pairing residual R_W compares the two.
    trial.pairing_current_work =
        vlasov_diag.u_full.upar_internal_face_energy_transfer +
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy -
        vlasov_diag.u_full.upar_interface_energy_removed;

    // Block: beam_work_.finish_from_midpoint().  Beam-off is a genuine
    // frozen fixture and must not run a no-op finish path that can alter
    // trajectory ledgers.
    if (frozen.beam_on) {
        std::vector<double>* beam_cell_work = field_particle_power_audit_enabled_
            ? &pairing_workspace_.cell_work.beam_delta_ke_cell : NULL;
        double* beam_boundary_work = field_particle_power_audit_enabled_
            ? &pairing_workspace_.cell_work.beam_delta_ke_boundary : NULL;
        beam_work_.finish_from_midpoint(frozen.beam_schedule, grid_,
                                        trial_force_fields_, frozen.time,
                                        frozen.dt, mpi_rank, mpi_size,
                                        beam_cell_work, beam_boundary_work);
    }
    trial.beam_force_work = frozen.beam_on
        ? beam_work_.last_field_work() : 0.0;
    trial.beam_trial = beam_work_;

    // Block: tail_work_.kick()
    if (tail_on) {
        double tail_force_work = 0.0;
        std::vector<double>* tail_cell_work = field_particle_power_audit_enabled_
            ? &pairing_workspace_.cell_work.tail_delta_ke_cell : NULL;
        double* tail_boundary_work = field_particle_power_audit_enabled_
            ? &pairing_workspace_.cell_work.tail_delta_ke_boundary : NULL;
        tail_work_.kick(grid_, trial_force_fields_, frozen.dt,
                        mpi_rank, mpi_size, &tail_force_work, tail_cell_work,
                        tail_boundary_work);
        trial.tail_force_work = tail_force_work;
        trial.tail_trial = tail_work_;
    } else {
        trial.tail_force_work = 0.0;
    }
    stage_sources_.bulk_upar_face_work =
        vlasov_diag.u_full.upar_internal_face_energy_transfer;
    stage_sources_.bulk_upar_velocity_boundary_work =
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy;
    stage_sources_.bulk_upar_interface_energy_removed =
        vlasov_diag.u_full.upar_interface_energy_removed;
    stage_sources_.bulk_upar_identity_residual =
        vlasov_diag.u_full.upar_discrete_energy_identity_residual;
    stage_sources_.tail_kick_work = trial.tail_force_work;
    stage_sources_.beam_kick_work = trial.beam_force_work;
    // The stage-energy audit records the field energy present during the force
    // step.  Legacy uses midpoint_fields_ here (the field does not change
    // during u_full); the DG trial's iterated pairing field is only the
    // candidate, so recording it would make U_E jump to the candidate value
    // and then back at x_half2, fabricating a field_coupling gap.
    capture_dg_stage(VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK, state_u_full_,
                     tail_on ? &tail_work_ : NULL,
                     frozen.beam_on ? &beam_work_ : &trial.beam_trial,
                     midpoint_fields_);

    // Block: apply_upar_flux_conversion
    if (observe_upar_flux) {
        VpfpStepResult trial_result = {};
        // Aggregate initialization leaves bool fields false.  audit_valid
        // must describe the production u_full audit, otherwise every valid
        // FLUX_INTERFACE trial is rejected before conversion even when all
        // parcel residuals are exactly zero.
        trial_result.audit_valid = vlasov_diag.u_full.audit_valid;
        trial_result.finite = vlasov_diag.u_full.finite;
        trial_result.conversion_ok = true;
        if (!apply_upar_flux_conversion(
                exported_flux, frozen.first_collision_flux,
                frozen.first_collision, observe_upar_flux,
                trial_result, mpi_rank, 10)) {
            failure.code = trial_result.failure_code;
            failure.stage = "trial_conversion";
            state_u_full_ = saved_state_u_full;
            state_np1_ = saved_state_np1;
            tail_work_ = saved_tail_work;
            beam_work_ = saved_beam_work;
            return false;
        }
        trial.conversion_number_removed =
            trial_result.ledger.conversion_number_removed;
        trial.conversion_energy_removed =
            trial_result.ledger.conversion_energy_removed;
        trial.conversion_particles_created =
            trial_result.ledger.conversion_particles_created;
        // §15.13.4: carry the per-cell conversion diagnostics through the
        // DG trial so Gate I can audit bulk/tail conversion source.
        trial.conversion_events = trial_result.conversion_events;
    } else {
        trial.conversion_number_removed = 0.0;
        trial.conversion_energy_removed = 0.0;
        trial.conversion_particles_created = 0;
    }
    stage_sources_.conversion_energy = trial.conversion_energy_removed;
    capture_dg_stage(VPFP_STAGE_CONVERSION_AFTER_FORCE, state_u_full_,
                     tail_on ? &tail_work_ : NULL,
                     frozen.beam_on ? &beam_work_ : &trial.beam_trial,
                     midpoint_fields_);

    // Block: vlasov_.second_x_half()
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_u_full_, pairing_workspace_.bulk_number_pre_x2);
    }
    if (!vlasov_.second_x_half(state_u_full_, state_np1_,
                               frozen.time + 0.5 * frozen.dt,
                               0.5 * frozen.dt, vlasov_diag,
                               field_particle_power_audit_enabled_
                                   ? &pairing_workspace_.bulk_x2 : NULL)) {
        failure.code = 2;
        failure.stage = "trial_second_x_half";
        state_u_full_ = saved_state_u_full;
        state_np1_ = saved_state_np1;
        tail_work_ = saved_tail_work;
        beam_work_ = saved_beam_work;
        return false;
    }
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_np1_, pairing_workspace_.bulk_number_post_x2);
    }
    stage_sources_.background_left_inflow_energy +=
        vlasov_diag.x_second.left_inflow_energy;
    stage_sources_.background_left_outflow_energy +=
        vlasov_diag.x_second.left_outflow_energy;
    stage_sources_.background_right_inflow_energy +=
        vlasov_diag.x_second.right_inflow_energy;
    stage_sources_.background_right_outflow_energy +=
        vlasov_diag.x_second.right_outflow_energy;

    // Block: tail_work_.drift_half()
    if (tail_on) {
        tail_work_.drift_half(grid_, 0.5 * frozen.dt, mpi_rank, mpi_size);
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
        tail_work_.finalize_trajectory_current(grid_, frozen.dt,
                                                mpi_rank, mpi_size);
        trial.tail_trial = tail_work_;
    }
    if (frozen.beam_on) {
        beam_work_.deposit_density(grid_, mpi_rank, mpi_size);
        beam_work_.finalize_charge_conserving_current(
            grid_, frozen.dt, mpi_rank, mpi_size);
        trial.beam_trial = beam_work_;
    }
    capture_dg_stage(VPFP_STAGE_X_HALF2, state_np1_,
                     tail_on ? &tail_work_ : NULL,
                     frozen.beam_on ? &beam_work_ : &trial.beam_trial,
                     midpoint_fields_);

    // §5.7: Final candidate Poisson and pairing field construction.
    const std::vector<double> empty_tail_density_local;
    const std::vector<double>& tail_trial_density = tail_on
        ? tail_work_.density : empty_tail_density_local;
    trial.final_fields_trial = final_fields_;
    trial.final_fields_trial.set_charge_density(
        state_np1_, tail_trial_density, beam_work_.density, frozen.ion_density);
    OpenGaussSolveOptions final_options;
    final_options.reconstruct_phi = true;
    final_options.compute_l1 = true;
    final_options.compute_boundary_audit = true;
    field_solver_.solve(trial.final_fields_trial, mpi_rank, mpi_size,
                        final_options);
    std::vector<double> rho_delta = trial.final_fields_trial.rho;
    if (rho_delta.size() == field_n_pairing_.rho.size()) {
        for (size_t i = 0; i < rho_delta.size(); ++i)
            rho_delta[i] -= field_n_pairing_.rho[i];
    }
    const OpenPoissonWorkIdentity poisson_identity =
        field_solver_.evaluate_work_identity(
            field_n_pairing_, trial.final_fields_trial, rho_delta,
            mpi_rank, mpi_size);
    // §5.7 final Poisson pass: the natural solve-quality measure is the Gauss
    // residual |dEx/dx - rho/eps0|_inf, which is ~machine epsilon relative to
    // the charge scale by construction of the direct-integration solve.  The
    // tolerance scales with max|rho/eps0| (not an absolute 1.0 J/m^2 floor),
    // so a genuinely near-zero field passes instead of failing on roundoff.
    trial.final_poisson_residual_linf =
        field_solver_.diagnostics().residual_linf;
    double local_max_rho_eps0 = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        local_max_rho_eps0 = std::max(local_max_rho_eps0,
            std::fabs(trial.final_fields_trial.rho[
                static_cast<size_t>(grid_.nghost + ix)] / Const::eps0));
    }
    // diagnostics().residual_linf is already an MPI-global maximum.  Its
    // comparison scale must therefore also be global; a rank-local scale
    // makes some ranks reject while others enter the next pairing-field
    // collective, which deadlocks the physical step.
    double max_rho_eps0 = 0.0;
    MPI_Allreduce(&local_max_rho_eps0, &max_rho_eps0, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    trial.final_poisson_tolerance =
        4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(max_rho_eps0, 1.0e-30);
    trial.final_poisson_pass =
        !trial.final_fields_trial.Ex_face.empty() &&
        std::isfinite(trial.final_poisson_residual_linf) &&
        trial.final_poisson_residual_linf <= trial.final_poisson_tolerance;
    trial.pairing_field_build_pass = field_solver_.build_potential_pairing_field(
        field_n_pairing_, trial.final_fields_trial,
        trial.pairing_field_map, mpi_rank, mpi_size);

    // JC3 (section 6.10): test-only fault injection at the trial boundary.
    if (jc_fault_.fail_final_poisson) trial.final_poisson_pass = false;
    if (jc_fault_.fail_pairing_build) trial.pairing_field_build_pass = false;

    // §6.7: exchange-scale signed pairing quantities from the production
    // Poisson work identity.  pairing_current_work is dt <E_pair, J> and is
    // the same field guess the particles were actually pushed with.
    trial.potential_charge_work = poisson_identity.potential_charge_work;
    trial.field_energy_change = poisson_identity.field_energy_change;
    trial.electrode_work = poisson_identity.electrode_work;

    // Capture trial results before restoring shared state
    trial.bulk_trial = state_np1_;
    trial.vlasov_diagnostics = vlasov_diag;
    // §5.9: snapshot ledger values for deterministic comparison
    trial.work_pairing_current = trial.bulk_force_work +
        trial.tail_force_work + trial.beam_force_work;
    trial.total_force_work = trial.work_pairing_current;
    trial.rho_np1 = trial.final_fields_trial.rho;

    // §5.8: Restore shared state so next trial starts from same frozen state
    state_u_full_ = saved_state_u_full;
    state_np1_ = saved_state_np1;
    tail_work_ = saved_tail_work;
    beam_work_ = saved_beam_work;

    return true;
}

bool VpfpIntegrator::compute_field_work_residuals(
    const std::vector<double>& pairing_guess,
    const FieldParticleTrial& trial,
    double& field_l2, double& field_linf, double& pairing_relative,
    bool& all_finite, int mpi_rank, int mpi_size) const
{
    (void)mpi_size;
    const int nxl = grid_.nx_local;
    const double dx = grid_.dx;
    const std::size_t nface = static_cast<std::size_t>(nxl + 1);

    // A structurally invalid local map (e.g. build_potential_pairing_field
    // cleared it after a non-finite result) must still let every rank enter
    // the same fixed collective sequence; otherwise one rank early-returns
    // while the others block in MPI_Allreduce (deadlock).
    const bool local_valid =
        pairing_guess.size() == nface &&
        trial.pairing_field_map.size() == nface;

    double sum_r2 = 0.0;
    double sum_e2 = 0.0;
    double max_r = 0.0;
    double max_e = 0.0;
    bool local_finite = true;
    if (local_valid) {
        for (std::size_t f = 0; f < nface; ++f) {
            // Trapezoid face weights, consistent with the face field energy
            // norm (dx/2 at the two physical endpoints, dx in the interior).
            const double w = (f == 0 || f == nface - 1) ? 0.5 * dx : dx;
            const double r = trial.pairing_field_map[f] - pairing_guess[f];
            const double e = trial.pairing_field_map[f];
            sum_r2 += w * r * r;
            sum_e2 += w * e * e;
            max_r = std::max(max_r, std::fabs(r));
            max_e = std::max(max_e, std::fabs(e));
            if (!std::isfinite(r) || !std::isfinite(e)) local_finite = false;
        }
    } else {
        local_finite = false;
    }

    // JC3 (section 6.3): compare like with like.  total_force_work is
    // rank-local Bulk+Tail+Beam Delta K.  The Bulk Gate-C work stored in
    // pairing_current_work is already global, while Tail/Beam kick works are
    // rank-local.  Contribute the global Bulk term once (rank 0) and add the
    // local PIC terms before the SUM reduction; comparing total particle
    // work against Bulk-only work creates a nonzero residual plateau.
    const double local_pairing_work =
        (mpi_rank == 0 ? trial.pairing_current_work : 0.0) +
        trial.tail_force_work + trial.beam_force_work;
    double local_pack[4] = {
        sum_r2, sum_e2, trial.total_force_work, local_pairing_work
    };
    double global_pack[4] = { 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_pack, global_pack, 4, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double local_max[2] = { max_r, max_e };
    double global_max[2] = { 0.0, 0.0 };
    MPI_Allreduce(local_max, global_max, 2, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    int local_finite_i = local_finite ? 1 : 0;
    int global_finite = 0;
    MPI_Allreduce(&local_finite_i, &global_finite, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    all_finite = global_finite != 0;

    const double e_floor = 1.0;   // V/m absolute floor (physical units)
    const double k_floor = 1.0;   // J/m^2 absolute exchange-work floor
    field_l2 = std::sqrt(global_pack[0]) /
        std::max(e_floor, std::sqrt(global_pack[1]));
    field_linf = global_max[0] / std::max(e_floor, global_max[1]);
    const double dk = global_pack[2];
    const double wp = global_pack[3];
    const double denom = std::max(k_floor,
        std::max(std::fabs(dk), std::fabs(wp)));
    pairing_relative = std::fabs(dk - wp) / denom;
    return all_finite;
}

bool VpfpIntegrator::apply_post_field_once_and_validate_charge(
    const FieldParticleFrozenState& frozen,
    FieldParticleTrial& trial, VpfpFailureInfo& failure)
{
    // Restore the converged trial into the shared work buffers so the
    // production post-field operators run exactly once on the candidate.
    state_np1_ = trial.bulk_trial;
    tail_work_ = trial.tail_trial;
    beam_work_ = trial.beam_trial;
    final_fields_ = trial.final_fields_trial;

    // JC3 (section 6.10): test-only post-field charge failure injection.
    if (jc_fault_.fail_post_field_charge) {
        failure.code = 206;
        failure.stage = "post_field_charge_invariance_failure";
        failure.failing_rank = frozen.mpi_rank;
        return false;
    }

    VpfpPostFieldChargeInvarianceReport report;
    if (!post_field_charge_invariance_transaction(
            state_np1_, tail_work_, frozen.time, frozen.dt,
            frozen.mpi_rank, frozen.mpi_size, report)) {
        failure.code = 206;
        failure.stage = "post_field_charge_invariance_failure";
        failure.failing_rank = frozen.mpi_rank;
        return false;
    }
    // JC0 gate: per-cell number residual inside tolerance and relative rho
    // change below 1e-12.  Violating it means C2/return changed the charge.
    const bool charge_ok =
        report.cell_number_residual_linf <= report.cell_number_tolerance &&
        report.rho_before_after_relative_linf <= 1.0e-12;
    // JC4 (section 7.7): expose the post-field charge residual for
    // diagnostics; stored after the charge gate so it reflects the accepted
    // value.
    post_field_charge_residual_linf_ = report.rho_before_after_relative_linf;
    if (!charge_ok) {
        failure.code = 206;
        failure.stage = "post_field_charge_invariance_failure";
        failure.failing_rank = frozen.mpi_rank;
        return false;
    }

    // Candidate now holds the post-C2/return state; reuse the converged
    // final Poisson (§6.4 item 5: charge was invariant).
    trial.bulk_trial = state_np1_;
    trial.tail_trial = tail_work_;
    trial.final_fields_trial = final_fields_;
    trial.second_collision_reservoir_energy =
        report.collision_reservoir_energy;
    trial.post_conversion_number_removed =
        report.conversion_number_removed;
    trial.post_conversion_energy_removed =
        report.conversion_energy_removed;
    trial.post_conversion_particles_created =
        report.conversion_particles_created;
    trial.post_tail_return = report.tail_return;
    trial.post_field_charge_residual_linf =
        report.rho_before_after_relative_linf;
    // §15.13.4: carry per-cell C2 conversion diagnostics for Gate I.
    trial.post_conversion_events = report.conversion_events;
    return true;
}

bool VpfpIntegrator::solve_field_particle_center(
    const FieldParticleFrozenState& frozen,
    FieldParticleTrial& accepted_trial,
    FieldParticleIterationDiagnostics& diagnostics,
    VpfpFailureInfo& failure)
{
    failure = VpfpFailureInfo();
    diagnostics = FieldParticleIterationDiagnostics();
    const int mpi_rank = frozen.mpi_rank;
    const int mpi_size = frozen.mpi_size;

    // Step-count transaction (§4.9/§6.6): every trial uses the same candidate
    // step; on any failure restore the pre-prepare value so the accepted
    // step/time and all accepted state are unchanged.
    const long long candidate_step = frozen.candidate_step;
    const long long step_before = candidate_step - 1;
    step_count_ = candidate_step;

    const FieldParticleCouplingConfig& cfg = field_particle_coupling_config_;
    const double field_tol = cfg.field_relative_tolerance;
    const double pairing_tol = cfg.pairing_relative_tolerance;

    // Initial guess: the legacy midpoint Poisson field (§6.1), used only as
    // the first E_pair candidate.  A test-only sign override lets the
    // positive/negative signed-field cases exercise the full Picard path with
    // a signed initial force field.
    pairing_field_guess_ = midpoint_fields_.Ex_face;
    if (jc_fault_.initial_guess_sign != 0) {
        const double sign =
            jc_fault_.initial_guess_sign > 0 ? 1.0 : -1.0;
        for (std::size_t f = 0; f < pairing_field_guess_.size(); ++f)
            pairing_field_guess_[f] *= sign;
    }
    double omega = cfg.initial_relaxation;
    double residual_previous = std::numeric_limits<double>::infinity();
    int consecutive_growth = 0;
    bool converged = false;

    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        // JC3 (section 6.10): test-only NaN injection into the trial force
        // field at a chosen iteration on a chosen rank.
        if (jc_fault_.nan_inject_iteration >= 1 &&
            (iter + 1) >= jc_fault_.nan_inject_iteration &&
            (jc_fault_.nan_inject_rank == -1 ||
             jc_fault_.nan_inject_rank == mpi_rank) &&
            !pairing_field_guess_.empty()) {
            pairing_field_guess_[0] =
                std::numeric_limits<double>::quiet_NaN();
        }

        FieldParticleTrial trial;
        VpfpFailureInfo trial_failure;
        if (!evaluate_field_particle_trial(frozen, pairing_field_guess_,
                                           trial, trial_failure)) {
            failure = trial_failure;
            if (failure.code == 0) failure.code = 201;
            failure.failing_rank = mpi_rank;
            diagnostics.failure_code = failure.code;
            // §6.5 failure log: step, time, iteration, stage, failing rank.
            std::fprintf(stderr,
                "[jc3-picard-fail] step=%lld time=%.17g iter=%d "
                "failure_code=%d stage=%s failing_rank=%d\n",
                candidate_step, frozen.time, iter + 1, failure.code,
                failure.stage.c_str(), failure.failing_rank);
            step_count_ = step_before;
            return false;
        }

        double field_l2 = 0.0, field_linf = 0.0, pairing_relative = 0.0;
        bool all_finite = false;
        if (!compute_field_work_residuals(pairing_field_guess_, trial,
                                          field_l2, field_linf, pairing_relative,
                                          all_finite, mpi_rank, mpi_size) ||
            !all_finite) {
            failure.code = 201;
            failure.stage = "field_particle_trial_nonfinite";
            failure.failing_rank = mpi_rank;
            diagnostics.failure_code = 201;
            std::fprintf(stderr,
                "[jc3-picard-fail] step=%lld time=%.17g iter=%d "
                "failure_code=201 stage=nonfinite failing_rank=%d\n",
                candidate_step, frozen.time, iter + 1, mpi_rank);
            step_count_ = step_before;
            return false;
        }
        if (!trial.final_poisson_pass) {
            failure.code = 202;
            failure.stage = "field_particle_poisson_failure";
            failure.failing_rank = mpi_rank;
            diagnostics.failure_code = 202;
            std::fprintf(stderr,
                "[jc3-picard-fail] step=%lld time=%.17g iter=%d "
                "failure_code=202 stage=poisson failing_rank=%d "
                "gauss_residual=%.17g gauss_tolerance=%.17g\n",
                candidate_step, frozen.time, iter + 1, mpi_rank,
                trial.final_poisson_residual_linf,
                trial.final_poisson_tolerance);
            step_count_ = step_before;
            return false;
        }
        if (!trial.pairing_field_build_pass) {
            failure.code = 203;
            failure.stage = "field_particle_pairing_field_failure";
            failure.failing_rank = mpi_rank;
            diagnostics.failure_code = 203;
            std::fprintf(stderr,
                "[jc3-picard-fail] step=%lld time=%.17g iter=%d "
                "failure_code=203 stage=pairing_build failing_rank=%d\n",
                candidate_step, frozen.time, iter + 1, mpi_rank);
            step_count_ = step_before;
            return false;
        }

        diagnostics.field_residual_l2 = field_l2;
        diagnostics.field_residual_linf = field_linf;
        diagnostics.pairing_residual = pairing_relative;
        diagnostics.relaxation = omega;
        diagnostics.iterations = iter + 1;

        // §6.5 per-iteration log: step, time, iteration, relaxation, both
        // field residuals and the work pairing residual.
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                "[jc3-picard] step=%lld time=%.17g iter=%d omega=%.17g "
                "field_l2=%.17g field_linf=%.17g pairing=%.17g "
                "bulk_work=%.17g pairing_work=%.17g\n",
                candidate_step, frozen.time, iter + 1, omega,
                field_l2, field_linf, pairing_relative,
                trial.bulk_force_work, trial.pairing_current_work);
        }

        if (field_l2 <= field_tol && field_linf <= field_tol &&
            pairing_relative <= pairing_tol) {
            converged = true;
            diagnostics.iterations = iter + 1;
            diagnostics.converged = true;
            diagnostics.failure_code = 0;
            accepted_trial = trial;
            break;
        }

        if (field_l2 > 1.25 * residual_previous) ++consecutive_growth;
        else consecutive_growth = 0;
        if (consecutive_growth >= 2) {
            omega = std::max(cfg.minimum_relaxation, 0.5 * omega);
            consecutive_growth = 0;
        }

        for (std::size_t f = 0; f < pairing_field_guess_.size(); ++f) {
            pairing_field_guess_[f] += omega *
                (trial.pairing_field_map[f] - pairing_field_guess_[f]);
        }
        // Synchronize the shared physical MPI face exactly once, reusing the
        // project's face-owner convention (right rank owns the shared face).
        if (mpi_size > 1 && !pairing_field_guess_.empty()) {
            const int left = mpi_rank > 0 ? mpi_rank - 1 : MPI_PROC_NULL;
            const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : MPI_PROC_NULL;
            const std::size_t nxl = static_cast<std::size_t>(grid_.nx_local);
            double from_left = 0.0;
            const double to_right = pairing_field_guess_[nxl];
            MPI_Sendrecv(&to_right, 1, MPI_DOUBLE, right, 7901,
                         &from_left, 1, MPI_DOUBLE, left, 7901,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (left != MPI_PROC_NULL) pairing_field_guess_[0] = from_left;
        }
        residual_previous = field_l2;
    }

    if (!converged || jc_fault_.force_not_converged) {
        failure.code = 205;
        failure.stage = "field_particle_not_converged";
        failure.failing_rank = mpi_rank;
        diagnostics.failure_code = 205;
        diagnostics.converged = false;
        std::fprintf(stderr,
            "[jc3-picard-fail] step=%lld time=%.17g iterations=%d "
            "failure_code=205 stage=not_converged failing_rank=%d "
            "last_field_l2=%.17g last_field_linf=%.17g "
            "last_pairing=%.17g field_tol=%.17g pairing_tol=%.17g\n",
            candidate_step, frozen.time, diagnostics.iterations, mpi_rank,
            diagnostics.field_residual_l2, diagnostics.field_residual_linf,
            diagnostics.pairing_residual, field_tol, pairing_tol);
        step_count_ = step_before;
        return false;
    }

    // Post-field C2/return and charge validation (§6.4/§6.8).
    if (!apply_post_field_once_and_validate_charge(frozen, accepted_trial,
                                                   failure)) {
        if (failure.code == 0) failure.code = 206;
        diagnostics.failure_code = failure.code;
        std::fprintf(stderr,
            "[jc3-picard-fail] step=%lld time=%.17g iterations=%d "
            "failure_code=%d stage=%s failing_rank=%d\n",
            candidate_step, frozen.time, diagnostics.iterations,
            failure.code, failure.stage.c_str(), failure.failing_rank);
        step_count_ = step_before;
        return false;
    }

    // Success: step_count_ stays at candidate_step (committed by the caller).
    return true;
}

VpfpStepResult VpfpIntegrator::advance_with_beam(
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const std::vector<double>& ion_density, double time, double dt,
    int mpi_rank, int mpi_size)
{
    VpfpStepResult result = {};
    result.stage_energy_audit_enabled = stage_energy_audit_enabled_;
    result.stage_energy_audit_valid = true;
    result.stage_energy_count = 0;
    result.tail_return_enabled =
        background_tail_enabled_ && tail_bulk_return_.config().enabled;
    result.split_used = false;
    result.finite = true;
    result.cfl_ok = dt > 0.0;
    result.gauss_ok = false;
    result.collision_ok = true;
    result.conversion_ok = true;
    result.tail_ok = true;
    result.failure_code = 0;
    result.failure_stage = "final_validation";
    // JC4 (section 7.7): legacy path has no discrete-gradient coupling.
    result.field_particle_coupling_enabled = false;
    result.post_field_charge_residual_linf = 0.0;
    result.failing_rank = -1;
    result.failing_ix = -1;
    result.failing_iupar = -1;
    result.failing_iuperp = -1;
    result.input_min = std::numeric_limits<double>::infinity();
    result.input_max = -std::numeric_limits<double>::infinity();
    result.output_min = std::numeric_limits<double>::infinity();
    result.output_max = -std::numeric_limits<double>::infinity();
    result.first_nonfinite_value = std::numeric_limits<double>::quiet_NaN();
    result.audit_valid = true;
    result.audit_failure_code = 0;
    result.audit_parcel_failure_reason =
        static_cast<int>(ParcelNodeFailureReason::None);
    result.audit_parcel_failure_rank = -1;
    result.audit_parcel_failure_ix = -1;
    result.audit_parcel_failure_face = -1;
    result.audit_parcel_failure_iuperp = -1;
    result.audit_parcel_failure_node_mass = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_target = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_node_sum = std::numeric_limits<double>::quiet_NaN();
    result.audit_parcel_failure_scale = std::numeric_limits<double>::quiet_NaN();
    result.audit_inplace_state_bitwise_equal = true;
    result.audit_inplace_rng_equal = true;
    result.audit_inplace_ledger_equal = true;
    result.tail_conversion_mode = static_cast<int>(tail_conversion_mode_);
    result.static_extractor_call_count = 0;
    result.flux_face_audit_count = 0;
    result.flux_face_audit_face_abs_sum = 0.0;
    result.flux_face_audit_parcel_abs_sum = 0.0;
    result.flux_face_audit_abs_error_sum = 0.0;
    result.flux_face_audit_max_relative = 0.0;
    result.flux_face_audit_abs_at_max_relative = 0.0;
    result.flux_face_audit_max_valid = false;
    if (!initialized_ || !result.cfl_ok) {
        result.failure_code = 1;
        return result;
    }
    const bool tail_on = background_tail_enabled_;
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        !collision_interface_zero_wall_validation_ &&
        !collision_interface_exporting_absorbing_) {
        result.failure_code = 12;
        result.collision_ok = false;
        result.conversion_ok = false;
        return result;
    }
    if (tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        !collision_.is_trivial() &&
        collision_interface_exporting_absorbing_ &&
        collision_.bulk_integrator() !=
            BulkCollisionIntegrator::CHANG_COOPER_FLUX) {
        result.failure_code = 12;
        result.collision_ok = false;
        result.conversion_ok = false;
        return result;
    }
    ++step_count_;

    const std::chrono::steady_clock::time_point step_begin =
        std::chrono::steady_clock::now();
    result.ledger.background_number_before =
        global_sum(electrons.total_particle_number());
    result.ledger.beam_number_before =
        global_sum(beam.total_particle_number(grid_));

    // Working beam: snapshot the accepted density, reset the per-step
    // ledgers, then generate the step's injection event table once.
    beam_work_.density = beam.density;
    beam_work_.begin_step(grid_, dt);
    const BeamInjectionSchedule schedule =
        beam.generate_injection_schedule(grid_, time, dt, mpi_rank);
    if (tail_on) {
        // Trial tail: copy the accepted state, snapshot the density and
        // reset the step ledgers (section 4.4).
        tail_work_ = tail_accepted_;
        tail_work_.begin_step(grid_, dt);
        result.ledger.tail_number_before =
            global_sum(tail_total_weight(tail_work_));
        result.ledger.tail_particle_count_before =
            global_sum_u64(tail_work_.particles.size());
        result.ledger.tail_kinetic_energy_before =
            global_sum(tail_total_kinetic_energy(tail_work_));
        result.ledger.combined_number_before =
            result.ledger.background_number_before +
            result.ledger.tail_number_before;
        trace_tail_stage(mpi_rank, time, "step_begin");
    }
    initialize_energy_ledger(electrons, beam, fields, tail_on);
    result.ledger.field_energy_before = accepted_field_energy_;
    result.ledger.background_kinetic_energy_before =
        accepted_background_kinetic_energy_;
    result.ledger.beam_kinetic_energy_before = accepted_beam_kinetic_energy_;
    if (!tail_on) result.ledger.tail_kinetic_energy_before = 0.0;

    VpfpStageEnergyRecord stage_sources = {};
    const auto refresh_tail_stage_source = [&]() {
        stage_sources.tail_outflow_energy = tail_on
            ? tail_work_.outflow_ledger().left_kinetic_energy +
                  tail_work_.outflow_ledger().right_kinetic_energy
            : 0.0;
        stage_sources.conversion_energy = result.ledger.conversion_energy_removed;
        stage_sources.beam_injected_energy = beam_work_.last_injected_energy();
        stage_sources.beam_outflow_energy = beam_work_.last_outflow_energy();
    };
    const auto capture_stage = [&](int stage_id, const Species& bulk,
                                   const EMFields& stage_fields) {
        refresh_tail_stage_source();
        // Before predict_to_midpoint(), beam_work_ contains reset trial
        // ledgers rather than the accepted Beam state. Preserve the accepted
        // Beam kinetic energy through collision_half1; subsequent stages use
        // the predicted/pushed transactional work state.
        const BeamPIC* stage_beam = stage_id < VPFP_STAGE_X_HALF1
            ? &beam : &beam_work_;
        capture_stage_energy(result, stage_id, bulk,
                             tail_on ? &tail_work_ : NULL, stage_beam,
                             stage_fields, stage_sources);
    };
    capture_stage(VPFP_STAGE_ACCEPTED_N, electrons, fields);

    // Gate I (section 4.6.1): initialize the pairing workspace at the
    // accepted state (diagnostic level 2 only).  Level 0/1 leave it off so
    // no per-step scan or collective is added to production.
    pairing_workspace_.enabled = field_particle_power_audit_enabled_;
    if (field_particle_power_audit_enabled_) {
        pairing_workspace_.bulk_x1.init(grid_.nx_local);
        pairing_workspace_.bulk_x2.init(grid_.nx_local);
        pairing_workspace_.bulk_number_pre_x1.clear();
        pairing_workspace_.bulk_number_post_x1.clear();
        pairing_workspace_.bulk_number_pre_x2.clear();
        pairing_workspace_.bulk_number_post_x2.clear();
        pairing_workspace_.cell_work.init(grid_.nx_local);
        pairing_workspace_.tail_work_ledger = 0.0;
        pairing_workspace_.beam_work_ledger = 0.0;
    }
    XFaceTransportAudit* x1_audit = field_particle_power_audit_enabled_
        ? &pairing_workspace_.bulk_x1 : NULL;
    XFaceTransportAudit* x2_audit = field_particle_power_audit_enabled_
        ? &pairing_workspace_.bulk_x2 : NULL;
    std::vector<double>* bulk_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.bulk_delta_ke_cell : NULL;
    std::vector<double>* tail_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.tail_delta_ke_cell : NULL;
    double* tail_boundary_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.tail_delta_ke_boundary : NULL;
    double bulk_work_global = 0.0;
    double tail_work_local = 0.0;

    const bool observe_upar_flux = tail_on &&
        (tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT ||
         tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE);
    const bool apply_upar_sink = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE;
    const bool observe_collision_flux = apply_upar_sink &&
        collision_interface_exporting_absorbing_ &&
        !collision_.is_trivial() &&
        collision_.bulk_integrator() == BulkCollisionIntegrator::CHANG_COOPER_FLUX;
    const bool flux_mode = observe_upar_flux;
    BulkTailFluxBatch exported_flux;
    BulkTailFluxBatch first_collision_flux;
    BulkTailFluxBatch second_collision_flux;
    // First Strang collision half-step through the same selector used by the
    // no-Beam production path.  This must include all configured hybrid
    // pairs when a tail backend is active.
    const std::chrono::steady_clock::time_point first_collision_begin =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "collision_half1_begin");
    Species* collision_input = &electrons;
    CollisionDiagnostics first_collision = {};
    HybridCollisionDiagnostics first_hybrid;
    bool first_hybrid_valid = false;
    if (hybrid_collision_active_ || !collision_.is_trivial()) {
        state_collision_trial_ = electrons;
        collision_input = &state_collision_trial_;
    }
    const bool first_collision_ok = apply_collision_half(
        *collision_input, tail_work_, time, 0.5 * dt, 0, mpi_rank,
        first_collision, first_hybrid, first_hybrid_valid,
        observe_collision_flux ? &first_collision_flux : NULL);
    result.wall_collision_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - first_collision_begin).count();
    trace_tail_stage(mpi_rank, time, "collision_half1_end");
    if (!first_collision_ok) {
        result.failure_code = 5;
        result.collision_ok = false;
        return result;
    }
    stage_sources.collision_reservoir_energy +=
        first_collision.reservoir_energy_change;
    capture_stage(VPFP_STAGE_COLLISION_HALF1, *collision_input, fields);
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            *collision_input, pairing_workspace_.bulk_number_pre_x1);
    }
    if (observe_collision_flux) {
        result.collision_flux_half_export_number[0] =
            first_collision.interface_export_number;
        result.collision_flux_half_export_energy[0] =
            first_collision.interface_export_energy;
        result.collision_flux_half_implicit_residual_linf[0] =
            first_collision.implicit_flux_residual_linf;
        result.collision_flux_half_cross_pair_residual_linf[0] =
            first_collision.cross_flux_pair_residual_linf;
        result.collision_flux_half_inward_clipped_number[0] =
            first_collision.interface_inward_clipped_number;
        result.collision_flux_half_parcel_count[0] =
            first_collision.interface_parcel_count;
        result.collision_flux_half_rollback_count[0] =
            first_collision.transaction_rollback_count;
    }

    const std::chrono::steady_clock::time_point before_vlasov =
        std::chrono::steady_clock::now();
    VlasovStepDiagnostics vlasov_diag;
    exported_flux.apply_interface_sink = apply_upar_sink;
    trace_tail_stage(mpi_rank, time, "vlasov_x1_begin");
    if (!vlasov_.first_x_half(*collision_input, state_x_half_, time, 0.5 * dt,
                              vlasov_diag, x1_audit)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.x_first, "first_x");
        result.finite = false;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "vlasov_x1_end");
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_x_half_, pairing_workspace_.bulk_number_post_x1);
    }

    stage_sources.background_left_inflow_energy +=
        vlasov_diag.x_first.left_inflow_energy;
    stage_sources.background_left_outflow_energy +=
        vlasov_diag.x_first.left_outflow_energy;
    stage_sources.background_right_inflow_energy +=
        vlasov_diag.x_first.right_inflow_energy;
    stage_sources.background_right_outflow_energy +=
        vlasov_diag.x_first.right_outflow_energy;
    // Tail: first spatial half-drift to x^{n+1/2} with open-boundary
    // truncation and MPI migration (section 8.2 step 3), then the midpoint
    // density deposit (step 5).
    if (tail_on) {
        const std::chrono::steady_clock::time_point tail_push_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_drift1_begin");
        tail_work_.drift_half(grid_, 0.5 * dt, mpi_rank, mpi_size);
        result.wall_tail_push_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_push_begin).count();
        result.wall_tail_migrate_seconds +=
            tail_work_.last_migration_seconds();
        trace_tail_stage(mpi_rank, time, "tail_drift1_end");
        const std::chrono::steady_clock::time_point tail_dep_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_deposit_mid_begin");
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
        result.wall_tail_deposit_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_dep_begin).count();
        trace_tail_stage(mpi_rank, time, "tail_deposit_mid_end");
    }

    // Beam: pure first drift to x^{n+1/2} (predict_to_midpoint) and the
    // midpoint density deposit.
    const std::chrono::steady_clock::time_point before_beam =
        std::chrono::steady_clock::now();
    beam.predict_to_midpoint(schedule, grid_, fields, time, dt, mpi_rank,
                             mpi_size, beam_work_);
    if (field_particle_power_audit_enabled_) {
        beam_work_.snapshot_midpoint_trajectory_current(
            grid_, 0.5 * dt, mpi_rank, mpi_size,
            pairing_workspace_.beam_trajectory.after_first_drift_current_face);
    }
    beam_work_.deposit_density(grid_, mpi_rank, mpi_size);
    capture_stage(VPFP_STAGE_X_HALF1, state_x_half_, fields);

    // Midpoint Poisson solve P[rho^{n+1/2}] on the cheap path.
    const std::vector<double> empty_tail_density;
    const std::vector<double>& tail_mid_density =
        tail_on ? tail_work_.density : empty_tail_density;
    midpoint_fields_.set_charge_density(state_x_half_, tail_mid_density,
                                        beam_work_.density, ion_density);
    OpenGaussSolveOptions midpoint_options;
    midpoint_options.reconstruct_phi = false;
    midpoint_options.compute_l1 = false;
    midpoint_options.compute_boundary_audit = false;
    const std::chrono::steady_clock::time_point before_field =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "field_mid_begin");
    field_solver_.solve(midpoint_fields_, mpi_rank, mpi_size, midpoint_options);
    result.wall_field_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before_field).count();
    result.mpi_collective_seconds +=
        field_solver_.diagnostics().mpi_collective_seconds;
    trace_tail_stage(mpi_rank, time, "field_mid_end");
    capture_stage(VPFP_STAGE_MIDPOINT_POISSON, state_x_half_, midpoint_fields_);

    // Background T_u(E^{n+1/2}, dt) with the velocity-tail threshold.
    const bool inplace_flux_audit = tail_on &&
        tail_conversion_mode_ == TailConversionMode::FLUX_AUDIT;
    trace_tail_stage(mpi_rank, time, "vlasov_u_begin");
    if (!vlasov_.u_full(state_x_half_, state_u_full_, midpoint_fields_,
                        time + 0.5 * dt, dt, vlasov_diag,
                        observe_upar_flux ? partition_ : NULL,
                        observe_upar_flux ? &exported_flux : NULL,
                        tail_flux_quadrature_order_, bulk_cell_work)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.u_full, "upar");
        record_upar_audit(result, vlasov_diag.u_full);
        result.finite = false;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "vlasov_u_end");
    if (inplace_flux_audit) {
        // See the no-Beam branch above.  This operator receives Beam and
        // tail state only through const/read-only ownership boundaries, so
        // the per-step full-array digest is redundant with the remap scratch
        // invariant and final checkpoint comparison.
        result.audit_inplace_state_bitwise_equal =
            vlasov_diag.u_full.audit_physical_state_bitwise_equal;
        result.audit_inplace_rng_equal = true;
        result.audit_inplace_ledger_equal = true;
        if (!result.audit_inplace_state_bitwise_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 64;
        }
        if (!result.audit_inplace_rng_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 128;
        }
        if (!result.audit_inplace_ledger_equal) {
            vlasov_diag.u_full.audit_valid = false;
            vlasov_diag.u_full.audit_failure_code |= 256;
        }
    }
    record_upar_audit(result, vlasov_diag.u_full);
    // Gate C (sections 7.3/7.6): copy the bulk u_parallel discrete-work
    // values (already-global RemapDiagnostics) into the stage source record.
    // They are captured into every subsequent stage and contributed once by
    // rank 0 in finalize_stage_energy_audit.
    stage_sources.bulk_upar_face_work =
        vlasov_diag.u_full.upar_internal_face_energy_transfer;
    stage_sources.bulk_upar_velocity_boundary_work =
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy;
    stage_sources.bulk_upar_interface_energy_removed =
        vlasov_diag.u_full.upar_interface_energy_removed;
    stage_sources.bulk_upar_identity_residual =
        vlasov_diag.u_full.upar_discrete_energy_identity_residual;
    bulk_work_global =
        vlasov_diag.u_full.upar_internal_face_energy_transfer +
        vlasov_diag.u_full.upar_left_velocity_boundary_energy +
        vlasov_diag.u_full.upar_right_velocity_boundary_energy -
        vlasov_diag.u_full.upar_interface_energy_removed;
    pairing_workspace_.bulk_work_ledger = bulk_work_global;
    result.flux_parcel_count = vlasov_diag.u_full.interface_parcel_count;
    result.flux_node_count = vlasov_diag.u_full.interface_node_count;
    result.flux_duplicate_count = vlasov_diag.u_full.interface_duplicate_count;
    result.flux_face_export_number =
        vlasov_diag.u_full.interface_face_export_number;
    result.flux_parcel_number =
        vlasov_diag.u_full.interface_parcel_number;
    result.flux_export_number = vlasov_diag.u_full.interface_export_number;
    result.flux_export_energy = vlasov_diag.u_full.interface_export_energy;
    result.flux_below_threshold_number =
        vlasov_diag.u_full.interface_below_threshold_number;
    result.flux_quadrature_error_max =
        vlasov_diag.u_full.interface_quadrature_error_max;
    result.flux_tail_owned_bulk_residual =
        vlasov_diag.u_full.tail_owned_bulk_residual;
    result.flux_tail_owned_expected_transfer_number =
        vlasov_diag.u_full.tail_owned_expected_transfer_number;
    result.flux_tail_owned_roundoff_discarded_number =
        vlasov_diag.u_full.tail_owned_roundoff_discarded_number;
    result.flux_quadrature_order = tail_flux_quadrature_order_;
    result.flux_max_supports = static_cast<int>(tail_flux_max_supports_);
    result.flux_below_threshold_number =
        vlasov_diag.u_full.interface_below_threshold_number;
    if (apply_upar_sink) {
        const double tail_owned_floor = 4096.0 *
            std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::fabs(
                vlasov_diag.u_full.tail_owned_expected_transfer_number));
        if (!(vlasov_diag.u_full.tail_owned_bulk_residual <=
              tail_owned_floor)) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "[tail-owned-balance-fail] failure_code=6 step=%lld "
                             "unexplained=%.17g expected_transfer=%.17g "
                             "face_export=%.17g threshold=%.17g\n",
                             step_count_,
                             vlasov_diag.u_full.tail_owned_bulk_residual,
                             vlasov_diag.u_full.tail_owned_expected_transfer_number,
                             vlasov_diag.u_full.interface_face_export_number,
                             tail_owned_floor);
            }
            result.failure_code = 6;
            result.failure_stage = "upar_tail_interface";
            result.tail_ok = false;
            result.finite = false;
            return result;
        }
    }
    const double tail_fraction =
        vlasov_diag.u_full.tail_number_loss /
        std::max(1.0, vlasov_diag.u_full.number_before);
    if (tail_fraction > Param::umax_loss_abort_fraction) {
        double local_max_e = 0.0;
        double local_max_adu_dt = 0.0;
        for (int ix = 0; ix < grid_.nx_local; ++ix) {
            const double e = std::fabs(
                midpoint_fields_.Ex[grid_.nghost + ix]);
            local_max_e = std::max(local_max_e, e);
            local_max_adu_dt = std::max(
                local_max_adu_dt,
                std::fabs(electrons.charge) * e * dt /
                    (electrons.mass * Const::c));
        }
        std::fprintf(stderr,
                     "[tail-abort] tail=%.6g tail_fraction=%.6g "
                     "before=%.6g local_max_E=%.6g local_max_a_u_dt=%.6g\n",
                     vlasov_diag.u_full.tail_number_loss, tail_fraction,
                     vlasov_diag.u_full.number_before, local_max_e,
                     local_max_adu_dt);
        result.failure_code = 7;
        result.finite = false;
        return result;
    }

    // Beam: kick with E^{n+1/2} at x^{n+1/2}, second half-drift, late
    // injections, open-boundary removal and MPI migration.
    trace_tail_stage(mpi_rank, time, "beam_finish_begin");
    std::vector<double>* beam_cell_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.beam_delta_ke_cell : NULL;
    double* beam_boundary_work = field_particle_power_audit_enabled_
        ? &pairing_workspace_.cell_work.beam_delta_ke_boundary : NULL;
    beam_work_.finish_from_midpoint(schedule, grid_, midpoint_fields_, time,
                                    dt, mpi_rank, mpi_size, beam_cell_work,
                                    beam_boundary_work);
    trace_tail_stage(mpi_rank, time, "beam_finish_end");
    double beam_work_local = beam_work_.last_field_work();
    pairing_workspace_.beam_work_ledger = beam_work_local;
    // Gate C (section 7.5/7.6): rank-local Beam field work, summed over ranks
    // in finalize_stage_energy_audit.  Never re-push or re-trace the beam.
    stage_sources.beam_kick_work = beam_work_.last_field_work();
    if (!tail_on) {
        capture_stage(VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK, state_u_full_,
                      midpoint_fields_);
        capture_stage(VPFP_STAGE_CONVERSION_AFTER_FORCE, state_u_full_,
                      midpoint_fields_);
    }

    // Tail: kick with the same E^{n+1/2} at x^{n+1/2} (section 8.2 step 9),
    // then convert the post-T_u bulk mass above the threshold into tail
    // particles at the midpoint spatial positions (step 11), then the second
    // spatial half-drift (step 12).
    if (tail_on) {
        const std::chrono::steady_clock::time_point tail_kick_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_kick_begin");
        double tail_kick_work_local = 0.0;
        tail_work_.kick(grid_, midpoint_fields_, dt, mpi_rank, mpi_size,
                        &tail_kick_work_local, tail_cell_work,
                        tail_boundary_work);
        tail_work_local = tail_kick_work_local;
        pairing_workspace_.tail_work_ledger = tail_kick_work_local;
        // Gate C (section 7.4/7.6): rank-local kick work, summed over ranks
        // in finalize_stage_energy_audit.
        stage_sources.tail_kick_work = tail_kick_work_local;
        result.wall_tail_push_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_kick_begin).count();
        trace_tail_stage(mpi_rank, time, "tail_kick_end");
        capture_stage(VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK, state_u_full_,
                      midpoint_fields_);
        trace_tail_stage(mpi_rank, time, "conversion_begin");
        if (!apply_upar_flux_conversion(
                exported_flux, first_collision_flux, first_collision,
                flux_mode, result, mpi_rank, 10)) {
            return result;
        }
        trace_tail_stage(mpi_rank, time, "conversion_end");
        capture_stage(VPFP_STAGE_CONVERSION_AFTER_FORCE, state_u_full_,
                      midpoint_fields_);
        {
            const std::chrono::steady_clock::time_point tail_push2_begin =
                std::chrono::steady_clock::now();
            trace_tail_stage(mpi_rank, time, "tail_drift2_begin");
            tail_work_.drift_half(grid_, 0.5 * dt, mpi_rank, mpi_size);
            result.wall_tail_push_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tail_push2_begin).count();
            result.wall_tail_migrate_seconds +=
                tail_work_.last_migration_seconds();
            trace_tail_stage(mpi_rank, time, "tail_drift2_end");
        }
        // Stage H5 population control (section 7.10): post-final-drift, so
        // no MPI-in-flight particle is merged; boundary cells are skipped
        // inside.  The final density deposit uses the post-control
        // representation, whose N/Xw moments (and therefore the CIC
        // density) are preserved to the compression tolerance.
        if (!fixed_state_field_particle_audit_mode_ &&
            population_controller_.enabled() &&
            population_controller_.active_step(
                static_cast<int>(step_count_))) {
            const TailPopulationController::Diagnostics diag =
                population_controller_.apply(
                    tail_work_, grid_, *partition_,
                    static_cast<int>(step_count_), mpi_rank);
            result.population_control_applied = diag.applied;
            result.population_control_local_groups =
                diag.groups_compressed + diag.groups_split;
            result.population_control_local_fallbacks =
                diag.compression_fallback_count;
            result.population_control_groups =
                global_sum_int(diag.groups_compressed +
                               diag.groups_split);
            result.population_control_fallbacks =
                global_sum_int(diag.compression_fallback_count);
            result.population_control_particles_before =
                global_sum_u64(diag.particles_before_local);
            result.population_control_particles_after =
                global_sum_u64(diag.particles_after_local);
            // Global control report (sum over ranks for counts, max for the
            // moment residuals), see advance_background.
            global_max_doubles(diag.max_residual,
                               result.population_control_max_residual);
        }
    }

    // Background second T_x(dt/2).
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_u_full_, pairing_workspace_.bulk_number_pre_x2);
    }
    trace_tail_stage(mpi_rank, time, "vlasov_x2_begin");
    if (!vlasov_.second_x_half(state_u_full_, state_np1_, time + 0.5 * dt,
                               0.5 * dt, vlasov_diag, x2_audit)) {
        result.failure_code = 2;
        record_remap_failure(result, vlasov_diag.x_second, "second_x");
        result.finite = false;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "vlasov_x2_end");
    if (field_particle_power_audit_enabled_) {
        capture_cell_integrated_number(
            state_np1_, pairing_workspace_.bulk_number_post_x2);
    }
    stage_sources.background_left_inflow_energy +=
        vlasov_diag.x_second.left_inflow_energy;
    stage_sources.background_left_outflow_energy +=
        vlasov_diag.x_second.left_outflow_energy;
    stage_sources.background_right_inflow_energy +=
        vlasov_diag.x_second.right_inflow_energy;
    stage_sources.background_right_outflow_energy +=
        vlasov_diag.x_second.right_outflow_energy;
    capture_stage(VPFP_STAGE_X_HALF2, state_np1_, midpoint_fields_);
    result.wall_vlasov_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before_vlasov).count();

    // Second Strang collision half-step through the shared selector.
    const std::chrono::steady_clock::time_point second_collision_begin =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "collision_half2_begin");
    CollisionDiagnostics collision_diag = {};
    HybridCollisionDiagnostics second_hybrid;
    bool second_hybrid_valid = false;
    const bool second_collision_ok = apply_collision_half(
        state_np1_, tail_work_, time + dt, 0.5 * dt, 1, mpi_rank,
        collision_diag, second_hybrid, second_hybrid_valid,
        observe_collision_flux ? &second_collision_flux : NULL);
    result.wall_collision_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - second_collision_begin).count();
    trace_tail_stage(mpi_rank, time, "collision_half2_end");
    if (!second_collision_ok) {
        result.failure_code = 5;
        result.collision_ok = false;
        return result;
    }
    stage_sources.collision_reservoir_energy +=
        collision_diag.reservoir_energy_change;
    capture_stage(VPFP_STAGE_COLLISION_HALF2, state_np1_, midpoint_fields_);
    if (observe_collision_flux) {
        result.collision_flux_half_export_number[1] =
            collision_diag.interface_export_number;
        result.collision_flux_half_export_energy[1] =
            collision_diag.interface_export_energy;
        result.collision_flux_half_implicit_residual_linf[1] =
            collision_diag.implicit_flux_residual_linf;
        result.collision_flux_half_cross_pair_residual_linf[1] =
            collision_diag.cross_flux_pair_residual_linf;
        result.collision_flux_half_inward_clipped_number[1] =
            collision_diag.interface_inward_clipped_number;
        result.collision_flux_half_parcel_count[1] =
            collision_diag.interface_parcel_count;
        result.collision_flux_half_rollback_count[1] =
            collision_diag.transaction_rollback_count;
    }
    if (observe_collision_flux &&
        tail_conversion_mode_ == TailConversionMode::FLUX_INTERFACE &&
        collision_diag.interface_parcel_count > 0) {
        // interface_parcel_count is globally reduced by apply_collision_half().
        // All ranks must enter the conversion transaction when any rank owns
        // an interface parcel, otherwise a local failure can split control
        // flow before later collectives.
        second_collision_flux.recompute(partition_->min_conversion_energy);
        append_flux_face_audit(second_collision_flux, result);
        const std::chrono::steady_clock::time_point conversion_begin =
            std::chrono::steady_clock::now();
        BulkTailConversionDiagnostics collision_conversion =
            converter_->convert_flux_batch(
                second_collision_flux, tail_work_, grid_, *partition_,
                static_cast<int>(step_count_),
                ConversionLocation::AFTER_COLLISION_HALF, mpi_rank,
                tail_flux_max_supports_);
        result.wall_conversion_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - conversion_begin).count();
        result.conversion_events.push_back(collision_conversion);
        result.ledger.conversion_number_removed +=
            collision_conversion.number_removed;
        result.ledger.conversion_px_removed +=
            collision_conversion.px_removed;
        result.ledger.conversion_energy_removed +=
            collision_conversion.energy_removed;
        result.ledger.conversion_particles_created +=
            collision_conversion.particles_created;
        result.flux_parcel_count += collision_diag.interface_parcel_count;
        result.flux_export_number += collision_diag.interface_export_number;
        result.flux_export_energy += collision_diag.interface_export_energy;
        result.flux_compression_fallback_count +=
            collision_conversion.compression_fallback_count;
        result.flux_roundoff_discarded_number +=
            collision_conversion.roundoff_discarded_number;
        result.flux_subcell_fallback_count +=
            collision_conversion.subcell_fallback_count;
        result.flux_support_limit_violation_count +=
            collision_conversion.support_limit_violation_count;
        result.flux_duplicate_id_count +=
            collision_conversion.duplicate_id_count;
        result.flux_face_ledger_mismatch_count +=
            collision_conversion.face_ledger_mismatch_count;
        result.ledger.conversion_number_residual_rel = std::max(
            result.ledger.conversion_number_residual_rel,
            collision_conversion.number_residual_rel);
        result.ledger.conversion_px_residual_rel = std::max(
            result.ledger.conversion_px_residual_rel,
            collision_conversion.px_residual_rel);
        result.ledger.conversion_energy_residual_rel = std::max(
            result.ledger.conversion_energy_residual_rel,
            collision_conversion.energy_residual_rel);
        const bool local_conversion_ok = collision_conversion.complete &&
            collision_conversion.conservative &&
            collision_conversion.fidelity_ok && collision_conversion.finite;
        int local_failure_reason = ConversionFailureNone;
        if (!local_conversion_ok) {
            if (!collision_conversion.finite)
                local_failure_reason = ConversionFailureNonfinite;
            else if (collision_conversion.support_limit_violation_count != 0)
                local_failure_reason = ConversionFailureSupportLimit;
            else if (collision_conversion.duplicate_id_count != 0)
                local_failure_reason = ConversionFailureDuplicateId;
            else if (collision_conversion.face_ledger_mismatch_count != 0)
                local_failure_reason = ConversionFailureFaceLedger;
            else if (!collision_conversion.conservative)
                local_failure_reason = ConversionFailureConservation;
            else if (!collision_conversion.fidelity_ok)
                local_failure_reason = ConversionFailureFidelity;
            else
                local_failure_reason = ConversionFailureIncomplete;
        }
        if (!synchronize_conversion_outcome(
                local_failure_reason, collision_conversion, mpi_rank, 8,
                result)) {
            return result;
        }
        result.conversion_ok = result.conversion_ok && local_conversion_ok;
    }
    if (observe_collision_flux) {
        result.collision_flux_export_number =
            first_collision.interface_export_number +
            collision_diag.interface_export_number;
        result.collision_flux_export_energy =
            first_collision.interface_export_energy +
            collision_diag.interface_export_energy;
        result.collision_flux_implicit_residual_linf = std::max(
            first_collision.implicit_flux_residual_linf,
            collision_diag.implicit_flux_residual_linf);
        result.collision_flux_cross_pair_residual_linf = std::max(
            first_collision.cross_flux_pair_residual_linf,
            collision_diag.cross_flux_pair_residual_linf);
        result.collision_flux_inward_clipped_number =
            first_collision.interface_inward_clipped_number +
            collision_diag.interface_inward_clipped_number;
        result.collision_flux_parcel_count =
            first_collision.interface_parcel_count +
            collision_diag.interface_parcel_count;
        result.collision_flux_rollback_count =
            first_collision.transaction_rollback_count +
            collision_diag.transaction_rollback_count;
    }

    capture_stage(VPFP_STAGE_CONVERSION_AFTER_COLLISION, state_np1_,
                  midpoint_fields_);

    // H10 is deliberately after the second collision-half interface export
    // and before the final combined charge deposit/Poisson solve.
    if (!apply_tail_bulk_return(state_np1_, tail_work_, result,
                                mpi_rank, mpi_size)) return result;
    capture_stage(VPFP_STAGE_TAIL_BULK_RETURN, state_np1_, midpoint_fields_);

    // Final accepted-state moments, Beam density, continuity diagnostics and
    // the full-audit Poisson solve.
    state_np1_.compute_moments();
    beam_work_.deposit_density(grid_, mpi_rank, mpi_size);
    beam_work_.finalize_charge_conserving_current(grid_, dt, mpi_rank,
                                                  mpi_size);
    if (tail_on) {
        const std::chrono::steady_clock::time_point tail_dep2_begin =
            std::chrono::steady_clock::now();
        trace_tail_stage(mpi_rank, time, "tail_deposit_final_begin");
        tail_work_.deposit_density(grid_, mpi_rank, mpi_size);
        tail_work_.finalize_trajectory_current(grid_, dt, mpi_rank,
                                               mpi_size);
        result.wall_tail_deposit_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tail_dep2_begin).count();
        trace_tail_stage(mpi_rank, time, "tail_deposit_final_end");
    }
    result.wall_beam_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - before_beam).count();
    const std::vector<double>& tail_final_density =
        tail_on ? tail_work_.density : empty_tail_density;
    final_fields_.set_charge_density(state_np1_, tail_final_density,
                                     beam_work_.density, ion_density);
    OpenGaussSolveOptions final_options;
    final_options.reconstruct_phi = false;
    const std::chrono::steady_clock::time_point final_field_begin =
        std::chrono::steady_clock::now();
    trace_tail_stage(mpi_rank, time, "field_final_begin");
    field_solver_.solve(final_fields_, mpi_rank, mpi_size, final_options);
    result.wall_field_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - final_field_begin).count();
    result.mpi_collective_seconds +=
        field_solver_.diagnostics().mpi_collective_seconds;
    trace_tail_stage(mpi_rank, time, "field_final_end");

    // Hard validation (section 7.4).
    trace_tail_stage(mpi_rank, time, "final_validation_begin");
    result.finite = finite_species(state_np1_);
    for (size_t i = 0; i < beam_work_.particles.size(); ++i) {
        const BeamParticle& p = beam_work_.particles[i];
        result.finite = result.finite && std::isfinite(p.x) &&
                        std::isfinite(p.px) && std::isfinite(p.weight);
    }
    for (size_t i = 0; i < final_fields_.Ex_face.size(); ++i) {
        result.finite = result.finite &&
                        std::isfinite(final_fields_.Ex_face[i]);
    }
    for (size_t i = 0; i < final_fields_.rho.size(); ++i) {
        result.finite = result.finite &&
                        std::isfinite(final_fields_.rho[i]);
    }
    int finite = result.finite ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &finite, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    result.finite = finite != 0;

    double local_rho_scale = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        local_rho_scale = std::max(
            local_rho_scale,
            std::fabs(final_fields_.rho[grid_.nghost + ix] / Const::eps0));
    }
    double rho_scale = 0.0;
    MPI_Allreduce(&local_rho_scale, &rho_scale, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double field_scale =
        std::max(std::fabs(final_fields_.Ex_face.front()),
                 std::fabs(final_fields_.Ex_face.back())) /
        (grid_.dx * static_cast<double>(grid_.nx_global));
    result.gauss_ok =
        std::isfinite(field_solver_.diagnostics().residual_linf) &&
        field_solver_.diagnostics().residual_linf <=
            1.0e-10 * std::max(1.0, rho_scale + field_scale);
    if (!result.finite || !result.gauss_ok) {
        result.failure_code = 3;
        record_remap_failure(result, vlasov_diag.x_second,
                             "final_validation");
        return result;
    }

    // Background mass ledger (same-order global sums).
    const double x_before = vlasov_diag.x_first.number_before;
    const double x_in = vlasov_diag.x_first.inflow_number +
                        vlasov_diag.x_second.inflow_number;
    const double x_out = vlasov_diag.x_first.outflow_number +
                         vlasov_diag.x_second.outflow_number;
    const double u_tail = vlasov_diag.u_full.tail_number_loss;
    const double conversion_removed =
        global_sum(result.ledger.conversion_number_removed);
    const double expected_after =
        x_before + x_in - x_out - u_tail - conversion_removed +
        result.tail_return.number;
    // H10 return is applied after x_second.  Validate the ledger against the
    // post-return bulk trial state, not the pre-return remap diagnostic.
    const double actual_background_after =
        global_sum(state_np1_.total_particle_number());
    const double ledger_error =
        std::fabs(actual_background_after - expected_after) /
        std::max(1.0, x_before);
    result.ledger.remap_ledger_residual = ledger_error;
    result.ledger.background_left_flux =
        (vlasov_diag.x_first.left_inflow_number -
         vlasov_diag.x_first.left_outflow_number) +
        (vlasov_diag.x_second.left_inflow_number -
         vlasov_diag.x_second.left_outflow_number);
    result.ledger.background_right_flux =
        (vlasov_diag.x_first.right_inflow_number -
         vlasov_diag.x_first.right_outflow_number) +
        (vlasov_diag.x_second.right_inflow_number -
         vlasov_diag.x_second.right_outflow_number);
    result.ledger.background_left_inflow_energy =
        vlasov_diag.x_first.left_inflow_energy +
        vlasov_diag.x_second.left_inflow_energy;
    result.ledger.background_left_outflow_energy =
        vlasov_diag.x_first.left_outflow_energy +
        vlasov_diag.x_second.left_outflow_energy;
    result.ledger.background_right_inflow_energy =
        vlasov_diag.x_first.right_inflow_energy +
        vlasov_diag.x_second.right_inflow_energy;
    result.ledger.background_right_outflow_energy =
        vlasov_diag.x_first.right_outflow_energy +
        vlasov_diag.x_second.right_outflow_energy;
    result.ledger.background_tail_number_loss = u_tail;
    result.ledger.background_tail_energy_loss =
        vlasov_diag.u_full.tail_energy_loss;
    if (ledger_error > 1.0e-9) {
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "[ledger-fail] failure_code=6 step=%lld "
                         "remap_ledger_residual=%.17g threshold=1e-9 "
                         "x_before=%.17g x_in=%.17g x_out=%.17g "
                         "u_tail=%.17g conversion_removed=%.17g "
                         "tail_return=%.17g pre_return_after=%.17g "
                         "actual_after=%.17g expected_after=%.17g\n",
                         step_count_, ledger_error, x_before, x_in, x_out,
                         u_tail, conversion_removed,
                         result.tail_return.number,
                         vlasov_diag.x_second.number_after,
                         actual_background_after, expected_after);
        }
        result.failure_code = 6;
        result.failure_stage = "background_ledger_after_tail_return";
        return result;
    }

    // Beam number balance: N_b^{n+1} = N_b^n + N_in - N_out.  The injected
    // number lives only on rank 0 and the outflow is split over the ranks
    // that own the exiting particles, so both must be global sums; otherwise
    // one rank sees a spurious imbalance, rejects the step and breaks the
    // collective synchronization (PMIx fence timeout on the cluster).
    const double beam_after =
        global_sum(beam_work_.total_particle_number(grid_));
    const double beam_in = global_sum(beam_work_.last_injected_number());
    const double beam_out = global_sum(beam_work_.last_outflow_number());
    const double beam_balance =
        std::fabs(beam_after -
                  (result.ledger.beam_number_before + beam_in - beam_out)) /
        std::max(1.0, result.ledger.beam_number_before + beam_in);
    result.ledger.beam_injected = beam_in;
    result.ledger.beam_outflow = beam_out;
    result.ledger.beam_injected_energy =
        global_sum(beam_work_.last_injected_energy());
    result.ledger.beam_outflow_energy =
        global_sum(beam_work_.last_outflow_energy());
    if (beam_balance > 1.0e-9) {
        result.failure_code = 8;
        return result;
    }

    // Tail validation (section 11): finite nonnegative weights, no failed
    // migration, number balance against the conversion source and the open
    // outflows, and the section 11.5 resource gates.
    if (tail_on) {
        result.tail_ok =
            tail_work_.finite() && tail_work_.nonnegative_weights() &&
            !tail_work_.migration_failed();
        if (!result.tail_ok) {
            result.failure_code = 11;
            result.finite = false;
            return result;
        }
        const double tail_after_w =
            global_sum(tail_total_weight(tail_work_));
        const double tail_out =
            global_sum(tail_work_.outflow_ledger().left_number +
                       tail_work_.outflow_ledger().right_number);
        const double expected_tail_after =
            result.ledger.tail_number_before + conversion_removed - tail_out -
            result.tail_return.number;
        const double tail_balance_scale =
            std::max(1.0, result.ledger.tail_number_before +
                              conversion_removed + tail_out);
        result.ledger.tail_number_balance_error =
            std::fabs(tail_after_w - expected_tail_after) /
            tail_balance_scale;
        result.ledger.tail_outflow_number = tail_out;
        result.ledger.tail_outflow_energy =
            global_sum(tail_work_.outflow_ledger().left_kinetic_energy +
                       tail_work_.outflow_ledger().right_kinetic_energy);
        if (result.ledger.tail_number_balance_error > 1.0e-9) {
            if (mpi_rank == 0) {
                std::fprintf(stderr,
                             "[tailbalance-fail] failure_code=6 step=%lld "
                             "tail_balance=%.17g threshold=1e-9 "
                             "tail_before=%.17g conversion=%.17g "
                             "tail_out=%.17g tail_after=%.17g "
                             "expected_after=%.17g scale=%.17g\n",
                             step_count_,
                             result.ledger.tail_number_balance_error,
                             result.ledger.tail_number_before,
                             conversion_removed, tail_out, tail_after_w,
                             expected_tail_after, tail_balance_scale);
            }
            result.failure_code = 6;
            result.failure_stage = "tail_ledger_after_tail_return";
            return result;
        }
        if (tail_max_particles_ > 0 &&
            tail_work_.particles.size() > tail_max_particles_) {
            result.failure_code = 11;
            result.finite = false;
            return result;
        }
        result.ledger.tail_number_after = tail_after_w;
        result.ledger.tail_kinetic_energy_after =
            global_sum(tail_total_kinetic_energy(tail_work_));
        result.ledger.combined_number_after =
            actual_background_after + result.ledger.tail_number_after;
        if (tail_max_number_fraction_ > 0.0 &&
            result.ledger.tail_number_after >
                tail_max_number_fraction_ *
                    std::max(1.0, result.ledger.combined_number_after)) {
            result.failure_code = 11;
            result.finite = false;
            return result;
        }
    }

    // Roundoff-scale nonnegativity of the background.
    double local_negative_mass = 0.0;
    double local_positive_mass = 0.0;
    for (int ix = 0; ix < grid_.nx_local; ++ix) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double mass =
                    state_np1_.f[idx3(grid_.nghost + ix, iv, imu)];
                if (mass < 0.0) local_negative_mass -= mass;
                else local_positive_mass += mass;
            }
        }
    }
    double debt[2] = { local_negative_mass, local_positive_mass };
    double global_debt[2] = { 0.0, 0.0 };
    MPI_Allreduce(debt, global_debt, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (global_debt[0] > 1.0e-10 * std::max(1.0, global_debt[1])) {
        if (mpi_rank == 0) {
            std::fprintf(stderr,
                         "[negmass-fail] failure_code=6 step=%lld "
                         "negative_mass=%.17g positive_mass=%.17g "
                         "ratio=%.17g threshold=1e-10\n",
                         step_count_, global_debt[0], global_debt[1],
                         global_debt[0] / std::max(1.0, global_debt[1]));
        }
        result.failure_code = 6;
        return result;
    }
    trace_tail_stage(mpi_rank, time, "final_validation_end");

    result.ledger.electrostatic_boundary_work =
        field_solver_.boundary_energy_work(fields, final_fields_,
                                           mpi_rank, mpi_size);
    if (mpi_rank == 0) {
        stage_sources.electrostatic_boundary_work =
            result.ledger.electrostatic_boundary_work;
    }
    capture_stage(VPFP_STAGE_FINAL_POISSON, state_np1_, final_fields_);

    // Transactional accept.
    electrons.swap_state(state_np1_);
    beam.swap_state(beam_work_);
    if (tail_on) {
        tail_accepted_.swap_state(tail_work_);
        result.ledger.tail_particle_count_after =
            global_sum_u64(tail_accepted_.particles.size());
    }
    swap_emfields(fields, final_fields_);
    beam.commit_injection_schedule(schedule, mpi_rank);
    result.ledger.background_number_after =
        global_sum(electrons.total_particle_number());
    result.ledger.beam_number_after =
        global_sum(beam.total_particle_number(grid_));
    result.ledger.gauss_charge_residual =
        field_solver_.diagnostics().boundary_charge_residual;
    result.ledger.field_energy = field_energy(fields);
    result.ledger.background_kinetic_energy =
        global_sum(electrons.total_kinetic_energy());
    result.ledger.beam_kinetic_energy =
        global_sum(beam.total_kinetic_energy());
    // Stage-H6 accepted cumulative ledgers and combined checksums (sections
    // 12.1/12.3).
    tail_cumulative_.conversion_number +=
        result.ledger.conversion_number_removed;
    tail_cumulative_.conversion_px += result.ledger.conversion_px_removed;
    tail_cumulative_.conversion_energy +=
        result.ledger.conversion_energy_removed;
    tail_cumulative_.return_number += result.tail_return.number;
    tail_cumulative_.return_px += result.tail_return.px;
    tail_cumulative_.return_jx_dx += result.tail_return.jx_dx;
    tail_cumulative_.return_energy += result.tail_return.energy;
    tail_cumulative_.return_pixx_dx += result.tail_return.pixx_dx;
    tail_cumulative_.return_piperp_dx += result.tail_return.piperp_dx;
    tail_cumulative_.return_particles_removed +=
        result.tail_return.particles_removed;
    tail_cumulative_.return_deferred_groups +=
        result.tail_return.deferred_infeasible_groups +
        result.tail_return.deferred_rank_boundary_groups;
    tail_cumulative_.particles_created +=
        result.ledger.conversion_particles_created;
    tail_cumulative_.control_groups += static_cast<std::uint64_t>(
        std::max(0, result.population_control_local_groups));
    tail_cumulative_.control_fallbacks += static_cast<std::uint64_t>(
        std::max(0, result.population_control_local_fallbacks));
    globalize_conversion_ledger(result.ledger);
    globalize_conversion_diagnostics(result);
    if (tail_on) {
        tail_cumulative_.outflow_number += static_cast<std::uint64_t>(
            tail_accepted_.outflow_ledger().left_number +
            tail_accepted_.outflow_ledger().right_number);
        combined_checksum_.number = result.ledger.combined_number_after;
        combined_checksum_.kinetic_energy =
            result.ledger.background_kinetic_energy +
            result.ledger.tail_kinetic_energy_after;
        result.tail_particles_local_max =
            global_max_u64(tail_accepted_.particles.size());
    } else {
        combined_checksum_.number = result.ledger.background_number_after;
        combined_checksum_.kinetic_energy =
            result.ledger.background_kinetic_energy;
        result.tail_particles_local_max = 0;
    }
    combined_checksum_.field_energy = result.ledger.field_energy;
    result.ledger.fct_energy_change = 0.0;
    finalize_collision_ledger(
        first_collision, collision_diag, first_hybrid, second_hybrid,
        first_hybrid_valid, second_hybrid_valid, result.ledger);
    finalize_energy_ledger(result);
    finalize_field_particle_power_audit(
        result, state_np1_, electrons, beam_work_, beam,
        tail_on ? &tail_work_ : NULL, tail_on ? &tail_accepted_ : NULL,
        final_fields_, fields, bulk_work_global, tail_work_local,
        beam_work_local, dt, mpi_rank, mpi_size);
    finalize_stage_energy_audit(result);
    result.wall_seconds_per_step = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - step_begin).count();
    trace_tail_stage(mpi_rank, time, "step_accepted");
#ifndef _WIN32
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        result.max_rss_kib = usage.ru_maxrss;
    }
#endif
    // PPM limiter activity (section 6.4): mean over the three remap
    // substeps of the constant/linear reduced-order cell fractions.
    result.remap_constant_fraction =
        (vlasov_diag.x_first.constant_fraction +
         vlasov_diag.u_full.constant_fraction +
         vlasov_diag.x_second.constant_fraction) / 3.0;
    result.remap_linear_fraction =
        (vlasov_diag.x_first.linear_fraction +
         vlasov_diag.u_full.linear_fraction +
         vlasov_diag.x_second.linear_fraction) / 3.0;
    result.accepted = true;
    return result;
}
