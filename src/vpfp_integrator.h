#ifndef VPFP_INTEGRATOR_H
#define VPFP_INTEGRATOR_H

#include "background_tail_pic.h"
#include "beam_pic.h"
#include "bulk_tail_converter.h"
#include "cylindrical_fp_collision.h"
#include "field_particle_power_audit.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "joint_phase_space_midpoint.h"
#include "hybrid_collision_step.h"
#include "tail_population_controller.h"
#include "tail_bulk_return.h"
#include "vlasov_split_step.h"

#include <cstdint>
#include <string>
#include <vector>

struct VpfpStepLedger {
    double background_number_before;
    double background_number_after;
    double beam_number_before;
    double beam_number_after;
    double background_left_flux;
    double background_right_flux;
    // Accepted-step kinetic-energy transfer through each physical
    // background boundary. These are copied directly from the two x-remap
    // face-flux ledgers; no cell-state reconstruction is used.
    double background_left_inflow_energy;
    double background_left_outflow_energy;
    double background_right_inflow_energy;
    double background_right_outflow_energy;
    double beam_injected;
    double beam_outflow;
    double beam_injected_energy;
    double beam_outflow_energy;
    double gauss_charge_residual;
    double field_energy;
    double background_kinetic_energy;
    double beam_kinetic_energy;
    // Accepted-state energy ledger.  These fields diagnose one complete
    // physical step and never participate in acceptance or state updates.
    double field_energy_before;
    double background_kinetic_energy_before;
    double beam_kinetic_energy_before;
    double electrostatic_boundary_work;
    double background_boundary_energy_net;
    double beam_boundary_energy_net;
    double domain_energy_before;
    double domain_energy_after;
    double domain_energy_change;
    double accounted_energy_source;
    double energy_balance_residual;
    double energy_balance_relative;
    double fct_energy_change;
    double collision_reservoir_energy;
    double collision_bulk_bulk_energy_change;
    double collision_tail_tail_px_change;
    double collision_tail_tail_energy_change;
    double collision_tail_bulk_px_change;
    double collision_tail_bulk_energy_change;
    double collision_bulk_reaction_px_change;
    double collision_bulk_reaction_energy_change;
    double collision_reaction_px_residual;
    double collision_reaction_energy_residual;
    int collision_pair_bulk_bulk;
    int collision_pair_tail_tail;
    int collision_pair_tail_bulk;
    int collision_pair_bulk_reaction;
    double background_tail_number_loss;
    double background_tail_energy_loss;
    double remap_ledger_residual;
    // Hybrid background-tail fields (stage H3, section 8/11): the tail is
    // part of the same background electron species, so the combined ledger
    // must close with the reservoir, boundary outflows and the internal
    // bulk-to-tail conversion.
    double tail_number_before;
    double tail_number_after;
    double tail_kinetic_energy_before;
    double tail_kinetic_energy_after;
    double tail_outflow_number;
    double tail_outflow_energy;
    double conversion_number_removed;
    double conversion_px_removed;
    double conversion_energy_removed;
    std::uint64_t conversion_particles_created;
    double conversion_number_residual_rel;
    double conversion_px_residual_rel;
    double conversion_energy_residual_rel;
    double combined_number_before;
    double combined_number_after;
    double tail_number_balance_error;
    // Macro-particle tail counts (global sums over ranks).  Added in stage
    // H5 as the population-control prerequisite: the H5 convergence matrix
    // compares tail particle numbers between thresholds/bins/resolutions
    // and the controller off/on A/B (section 15 H5).
    std::uint64_t tail_particle_count_before;
    std::uint64_t tail_particle_count_after;
    int tail_conversion_mode;
    std::uint64_t flux_parcel_count;
    std::uint64_t flux_node_count;
    std::uint64_t flux_duplicate_count;
    double flux_export_number;
    double flux_export_energy;
    double flux_below_threshold_number;
    double flux_tail_owned_expected_transfer_number;
    double flux_tail_owned_roundoff_discarded_number;
    double flux_tail_owned_bulk_residual;
    int flux_quadrature_order;
    int flux_max_supports;
    int static_extractor_call_count;
    // §16 条件3 核查：收敛 trial 实际用于推粒子的 E_pair 场（G_P(Phi^n,Phi^{n+1})，
    // 即 trial_force_fields_ 在收敛时的值）的全局场能 U_E(E_pair)。 与
    // field_energy（最终 Poisson 场能）对比，可核查"粒子实际受的场"与
    // "场能账本对应的场"是否为同一离散对象。 仅 discrete-gradient 路径填充。
    double e_pair_field_energy;
};

// Section 17.15.7.2: fixed, accepted-step-only snapshots of the actual
// production trial states.  Values are globalized in one packed reduction
// immediately before acceptance; no per-stage MPI operation is performed.
enum VpfpStageEnergyId {
    VPFP_STAGE_ACCEPTED_N = 0,
    VPFP_STAGE_COLLISION_HALF1,
    VPFP_STAGE_X_HALF1,
    VPFP_STAGE_MIDPOINT_POISSON,
    VPFP_STAGE_U_FORCE_TAIL_BEAM_KICK,
    VPFP_STAGE_CONVERSION_AFTER_FORCE,
    VPFP_STAGE_X_HALF2,
    VPFP_STAGE_COLLISION_HALF2,
    VPFP_STAGE_CONVERSION_AFTER_COLLISION,
    VPFP_STAGE_TAIL_BULK_RETURN,
    VPFP_STAGE_FINAL_POISSON,
    VPFP_STAGE_ENERGY_RECORD_COUNT
};

struct VpfpStageEnergyRecord {
    int stage_id;
    double bulk_kinetic;
    double tail_kinetic;
    double beam_kinetic;
    double field_energy;
    double delta_bulk_kinetic;
    double delta_tail_kinetic;
    double delta_beam_kinetic;
    double delta_field_energy;
    double background_left_inflow_energy;
    double background_left_outflow_energy;
    double background_right_inflow_energy;
    double background_right_outflow_energy;
    double beam_injected_energy;
    double beam_outflow_energy;
    double tail_outflow_energy;
    double collision_reservoir_energy;
    double conversion_energy;
    double tail_return_energy;
    double electrostatic_boundary_work;
    double stage_balance;
    // Gate C (section 7.6): read-only discrete-work records.  The four bulk
    // u_parallel fields are already-global RemapDiagnostics contributed once
    // by rank 0 in the packed SUM; tail_kick_work and beam_kick_work are
    // rank-local and summed over all ranks.
    double bulk_upar_face_work;
    double bulk_upar_velocity_boundary_work;
    double bulk_upar_interface_energy_removed;
    double bulk_upar_identity_residual;
    double tail_kick_work;
    double beam_kick_work;
};

// JC1 (section 4.2): field-particle coupling mode selection.
enum class FieldParticleCouplingMode {
    Legacy,
    DiscreteGradient
};

// J1 joint phase-space midpoint mode.  The default keeps the existing
// collisionless Strang/PPM path byte-for-byte outside this explicit mode.
enum class BackgroundPhaseSpaceMode {
    STRANG_PPM = 0,
    JOINT_MIDPOINT_ENERGY = 1
};

// JC1 (section 4.2): Picard iteration configuration for the discrete-gradient
// field-particle coupling.  Production default is Legacy; DiscreteGradient is
// only activated by explicit CLI opt-in.
struct FieldParticleCouplingConfig {
    FieldParticleCouplingMode mode = FieldParticleCouplingMode::Legacy;
    int max_iterations = 12;
    double initial_relaxation = 0.5;
    double minimum_relaxation = 0.125;
    double maximum_relaxation = 1.0;
    double field_relative_tolerance = 1.0e-8;
    double pairing_relative_tolerance = 1.0e-8;
};

// JC1 (section 4.3): read-only snapshot of the accepted state at the start
// of a physical step.  Only scalars and const pointers — no owned arrays.
struct FieldParticleFrozenState {
    const Species* bulk_x_half;
    const BackgroundTailPIC* tail_midpoint;
    const BackgroundTailPIC* tail_work_snapshot;
    const BeamPIC* beam_midpoint;
    const BeamPIC* beam_work_snapshot;
    const EMFields* field_n_pairing;
    BeamInjectionSchedule beam_schedule;
    BulkTailFluxBatch first_collision_flux;
    CollisionDiagnostics first_collision;
    HybridCollisionDiagnostics first_hybrid;
    bool first_hybrid_valid;
    VlasovStepDiagnostics x1_diagnostics;
    long long candidate_step;
    double time;
    double dt;
    bool beam_on;
    bool tail_on;
    std::vector<double> ion_density;
    int mpi_rank;
    int mpi_size;
    double bootstrap_residual;
    bool bootstrap_pass;
};

// JC1 (section 4.3): independent trial buffers for one Picard iteration.
// Overwritten each iteration; never aliases accepted state.
struct FieldParticleTrial {
    Species bulk_trial;
    BackgroundTailPIC tail_trial;
    BeamPIC beam_trial;
    EMFields final_fields_trial;
    std::vector<double> rho_np1;
    std::vector<double> pairing_field_map;
    double pairing_residual;
    // §5.9: ledger snapshots for deterministic comparison
    double conversion_number_removed;
    double conversion_energy_removed;
    std::uint64_t conversion_particles_created;
    double work_pairing_current;
    double bulk_force_work;
    double tail_force_work;
    double beam_force_work;
    double total_force_work;
    std::vector<double> force_field_face;
    bool final_poisson_pass;
    bool pairing_field_build_pass;
    double final_poisson_residual_linf;
    double final_poisson_tolerance;
    // §6.7: independent force-work computation from the production u_parallel
    // remap's Gate-C bulk work (already-global sums).  It is the same force
    // step as bulk_force_work, so at the self-consistent state the pairing
    // residual R_W = |total_force_work - pairing_current_work|/scale is ~0.
    double pairing_current_work;
    double potential_charge_work;
    double field_energy_change;
    double electrode_work;
    // Accepted-step diagnostics captured from the same deterministic trial;
    // these never feed back into the field-particle map.
    VlasovStepDiagnostics vlasov_diagnostics;
    double second_collision_reservoir_energy;
    double post_conversion_number_removed;
    double post_conversion_energy_removed;
    std::uint64_t post_conversion_particles_created;
    TailBulkReturnDiagnostics post_tail_return;
    double post_field_charge_residual_linf;
    // §15.13.4: per-cell conversion diagnostics carried through the DG trial
    // so that advance_discrete_gradient() can pass them to Gate I.  These are
    // audit-only and never modify the physical state.
    std::vector<BulkTailConversionDiagnostics> conversion_events;
    std::vector<BulkTailConversionDiagnostics> post_conversion_events;
};

// JC1 (section 4.3): per-iteration diagnostics.
struct FieldParticleIterationDiagnostics {
    int iterations = 0;
    double relaxation = 0.0;
    double field_residual_l2 = 0.0;
    double field_residual_linf = 0.0;
    double pairing_residual = 0.0;
    bool converged = false;
    int failure_code = 0;
};

// JC1 (section 4.4): failure information passed through helper interfaces.
struct VpfpFailureInfo {
    int code = 0;
    std::string stage;
    int failing_rank = 0;
    int failing_ix = 0;
};

struct VpfpStepResult {
    bool accepted;
    bool split_used;
    bool finite;
    bool cfl_ok;
    bool gauss_ok;
    bool collision_ok;
    bool conversion_ok;
    bool tail_ok;
    int failure_code;
    std::string failure_stage;
    int failing_rank;
    int failing_ix;
    int failing_iupar;
    int failing_iuperp;
    double input_min;
    double input_max;
    double output_min;
    double output_max;
    double first_nonfinite_value;
    bool audit_valid;
    int audit_failure_code;
    int audit_parcel_failure_reason;
    int audit_parcel_failure_rank;
    int audit_parcel_failure_ix;
    int audit_parcel_failure_face;
    int audit_parcel_failure_iuperp;
    double audit_parcel_failure_node_mass;
    double audit_parcel_failure_target;
    double audit_parcel_failure_node_sum;
    double audit_parcel_failure_scale;
    bool audit_inplace_state_bitwise_equal;
    bool audit_inplace_rng_equal;
    bool audit_inplace_ledger_equal;
    // Mean PPM limiter activation over the three remap substeps
    // (T_x(dt/2), T_u(dt), T_x(dt/2)): fractions of cells reconstructed at
    // reduced (constant / linear) order (sections 6.4/12.3, stage-5 record).
    double remap_constant_fraction;
    double remap_linear_fraction;
    VpfpStepLedger ledger;
    double wall_vlasov_seconds;
    double wall_field_seconds;
    double wall_beam_seconds;
    double wall_collision_seconds;
    double wall_seconds_per_step;
    double mpi_collective_seconds;
    long max_rss_kib;
    // Stage-H6 per-step tail timings (section 13.3).  wall_diagnostics is
    // filled by the caller with the previous accepted step's diagnostics
    // wall time; the others are measured inside advance.
    double wall_tail_push_seconds;
    double wall_tail_deposit_seconds;
    double wall_tail_migrate_seconds;
    double wall_conversion_seconds;
    double wall_diagnostics_seconds;
    // JC3 (section 6.12): discrete-gradient Picard diagnostics on the
    // accepted step.  These never participate in acceptance.
    bool field_particle_coupling_enabled;
    bool field_particle_converged;
    int field_particle_iterations;
    int field_particle_trial_evaluations;
    double field_particle_relaxation;
    double field_particle_residual_l2;
    double field_particle_residual_linf;
    double field_particle_pairing_residual;
    // JC4 (section 7.7): post-field charge residual from the final
    // C2/return validation; zero for legacy runs that skip this check.
    double post_field_charge_residual_linf;
    // Stage-H6 (section 13.3): global maximum of the per-rank tail
    // macro-particle counts after the step.
    std::uint64_t tail_particles_local_max;
    // Section 7.11.17 accepted/trial flux-interface audit fields.  These are
    // result-level values because the accepted ledger remains reserved for
    // physical number/energy balances.
    int tail_conversion_mode;
    std::uint64_t flux_parcel_count;
    std::uint64_t flux_node_count;
    std::uint64_t flux_duplicate_count;
    double flux_face_export_number;
    double flux_parcel_number;
    double flux_export_number;
    double flux_export_energy;
    double collision_flux_export_number;
    double collision_flux_export_energy;
    double collision_flux_implicit_residual_linf;
    double collision_flux_cross_pair_residual_linf;
    double collision_flux_inward_clipped_number;
    std::uint64_t collision_flux_parcel_count;
    int collision_flux_rollback_count;
    // Per-collision-half accepted diagnostics.  The aggregate fields above
    // remain convenient for existing readers; these pairs preserve the
    // required AFTER_COLLISION_HALF_1/2 provenance.
    double collision_flux_half_export_number[2];
    double collision_flux_half_export_energy[2];
    double collision_flux_half_implicit_residual_linf[2];
    double collision_flux_half_cross_pair_residual_linf[2];
    double collision_flux_half_inward_clipped_number[2];
    std::uint64_t collision_flux_half_parcel_count[2];
    int collision_flux_half_rollback_count[2];
    double flux_below_threshold_number;
    double flux_roundoff_discarded_number;
    double flux_quadrature_error_max;
    double flux_tail_owned_expected_transfer_number;
    double flux_tail_owned_roundoff_discarded_number;
    double flux_tail_owned_bulk_residual;
    int flux_quadrature_order;
    int flux_max_supports;
    std::uint64_t flux_compression_fallback_count;
    std::uint64_t flux_subcell_fallback_count;
    std::uint64_t flux_support_limit_violation_count;
    std::uint64_t flux_duplicate_id_count;
    std::uint64_t flux_face_ledger_mismatch_count;
    double flux_conversion_wall_seconds;
    int static_extractor_call_count;
    // Read-only flux-audit summary.  Keep only scalar closures plus one
    // deterministic worst face; gathering every interface row is not a
    // production-safe diagnostic protocol.
    std::uint64_t flux_face_audit_count;
    double flux_face_audit_face_abs_sum;
    double flux_face_audit_parcel_abs_sum;
    double flux_face_audit_abs_error_sum;
    double flux_face_audit_max_relative;
    double flux_face_audit_abs_at_max_relative;
    bool flux_face_audit_max_valid;
    BulkTailFluxFaceAudit flux_face_audit_max;
    // Global maxima over all conversion events in this candidate step.  The
    // accepted-step diagnostics are written by rank 0, so these values must
    // be reduced before the writer is called.
    double conversion_number_residual;
    double conversion_px_residual;
    double conversion_energy_residual;
    double conversion_jx_residual;
    double conversion_pixx_residual;
    double conversion_piperp_residual;
    double conversion_rho_l2;
    // Every conversion event performed by this accepted-step candidate.
    // Keeping a vector prevents a later location from overwriting an earlier
    // source when the split sequence grows additional conversion hooks.
    std::vector<BulkTailConversionDiagnostics> conversion_events;
    // Stage-H5 population-control report (section 7.10); all zero when the
    // controller is disabled or the step is not a control step.
    bool population_control_applied;
    int population_control_groups;
    int population_control_fallbacks;
    // Per-rank-local control counters (accumulated into the accepted
    // cumulative ledger on accept; the global fields above are the report).
    int population_control_local_groups;
    int population_control_local_fallbacks;
    std::uint64_t population_control_particles_before;
    std::uint64_t population_control_particles_after;
    double population_control_max_residual[7];
    // H10 is reported for every accepted step when enabled, including the
    // residence-only steps before the first particle returns to bulk.
    bool tail_return_enabled;
    TailBulkReturnDiagnostics tail_return;
    bool stage_energy_audit_enabled;
    bool stage_energy_audit_valid;
    int stage_energy_count;
    VpfpStageEnergyRecord stage_energy[VPFP_STAGE_ENERGY_RECORD_COUNT];
    // Gate I (section 4.6.1): read-only field-particle power pairing audit,
    // filled only on accepted steps when diagnostic level 2 enabled the
    // pairing audit.  pairing_audit is the globalized result of the single
    // packed reduction.
    bool pairing_audit_enabled;
    FieldParticlePowerAuditResult pairing_audit;
    bool joint_midpoint_enabled;
    bool joint_midpoint_converged;
    int joint_midpoint_iterations;
    double joint_midpoint_residual_linf;
    double joint_midpoint_poisson_residual_linf;
    double joint_midpoint_energy_residual;
    bool joint_midpoint_pairing_field_built;
    double joint_midpoint_delta_k_x;
    double joint_midpoint_delta_k_u;
    double joint_midpoint_u_face_work;
    double joint_midpoint_force_current_work;
    double joint_midpoint_charge_current_work;
    double joint_midpoint_charge_current_work_interior;
    double joint_midpoint_charge_current_work_endpoint;
    double joint_midpoint_poisson_potential_charge_work;
    double joint_midpoint_poisson_transport_residual;
    double joint_midpoint_current_pair_residual;
    double joint_midpoint_force_charge_residual;
    // Stage-B1 seam diagnostics (docs/VPFP_F10情形B根因定位与严格修复实施方案.md
    // sections 9-15).  Diagnostics only: they never participate in flux,
    // residual, Newton, Poisson, energy gate or acceptance logic.
    double joint_midpoint_pairing_face_left;
    double joint_midpoint_pairing_face_right;
    double joint_midpoint_force_current_first_cell;
    double joint_midpoint_force_current_last_cell;
    double joint_midpoint_naive_force_current_work;
    double joint_midpoint_seam_predicted_residual;
    double joint_midpoint_seam_prediction_error;
    double joint_midpoint_domain_energy_change;
    double joint_midpoint_field_energy_change;
    double joint_midpoint_electrode_work;
    double joint_midpoint_min_mass;
    double joint_midpoint_max_mass;
    std::vector<JointPhaseSpaceIterationRecord> joint_midpoint_iterations_log;
};

// Stage-H6 accepted-level cumulative tail ledgers (per-rank local) and
// combined checksums (global) saved by the checkpoint (section 12.1) and
// restored on restart (section 12.3).  The conversion/control counters
// accumulate on every accepted step.
struct VpfpTailCumulativeLedger {
    double conversion_number;
    double conversion_px;
    double conversion_energy;
    std::uint64_t particles_created;
    std::uint64_t outflow_number;
    std::uint64_t control_groups;
    std::uint64_t control_fallbacks;
    double return_number;
    double return_px;
    double return_jx_dx;
    double return_energy;
    double return_pixx_dx;
    double return_piperp_dx;
    std::uint64_t return_particles_removed;
    std::uint64_t return_deferred_groups;
    VpfpTailCumulativeLedger()
        : conversion_number(0.0), conversion_px(0.0),
          conversion_energy(0.0), particles_created(0), outflow_number(0),
          control_groups(0), control_fallbacks(0), return_number(0.0),
          return_px(0.0), return_jx_dx(0.0), return_energy(0.0),
          return_pixx_dx(0.0), return_piperp_dx(0.0),
          return_particles_removed(0), return_deferred_groups(0)
    {}
};

struct VpfpCombinedChecksum {
    double number;
    double kinetic_energy;
    double field_energy;
    VpfpCombinedChecksum()
        : number(0.0), kinetic_energy(0.0), field_energy(0.0)
    {}
};

// JC0 (section 3): read-only post-field charge-invariance transaction report.
// The helper runs the production second collision half C2, the collision-face
// conversion and the H10 tail-to-bulk return on caller-owned trial states and
// reports the combined (bulk + tail) electron number per local cell before
// and after.  It never touches the integrator's accepted state, cumulative
// ledgers, particle-ID counter or collision RNG ownership.  Every residual and
// tolerance is reduced so the report is rank-consistent.
struct VpfpPostFieldChargeInvarianceReport {
    std::vector<double> combined_number_before;  // per local cell, integrated
    std::vector<double> combined_number_after;
    double rank_number_before;
    double rank_number_after;
    double global_number_before;
    double global_number_after;
    double cell_number_residual_linf;   // max over cells and ranks
    double cell_number_tolerance;
    double rank_number_residual_linf;   // max over ranks
    double rank_number_tolerance;
    double global_number_residual;
    double global_number_tolerance;
    double rho_before_after_relative_linf;
    int collision_half_calls;             // production C2 invocation count
    int collision_face_conversion_calls;  // 0 or 1
    int tail_return_calls;                // 0 or 1
    bool tail_collision_stochastic_ran;   // C2 ran a stochastic tail backend
    bool c2_ok;
    bool conversion_ok;
    bool return_ok;
    double collision_face_export_number;  // global interface export (C2)
    double collision_reservoir_energy;
    double conversion_number_removed;
    double conversion_energy_removed;
    std::uint64_t conversion_particles_created;
    TailBulkReturnDiagnostics tail_return;
    // §15.13.4: per-cell C2 conversion diagnostics for Gate I audit.
    std::vector<BulkTailConversionDiagnostics> conversion_events;
    VpfpPostFieldChargeInvarianceReport()
        : rank_number_before(0.0), rank_number_after(0.0),
          global_number_before(0.0), global_number_after(0.0),
          cell_number_residual_linf(0.0), cell_number_tolerance(0.0),
          rank_number_residual_linf(0.0), rank_number_tolerance(0.0),
          global_number_residual(0.0), global_number_tolerance(0.0),
          rho_before_after_relative_linf(0.0), collision_half_calls(0),
          collision_face_conversion_calls(0), tail_return_calls(0),
          tail_collision_stochastic_ran(false), c2_ok(false),
          conversion_ok(false), return_ok(false),
          collision_face_export_number(0.0),
          collision_reservoir_energy(0.0),
          conversion_number_removed(0.0), conversion_energy_removed(0.0),
          conversion_particles_created(0)
    {}
};

// Open-boundary VPFP path (sections 7.1, 13.8).  Each physical step is the
// symmetric split
//   C(dt/2) -> T_x(dt/2) -> P[rho^{n+1/2}] -> T_u(E^{n+1/2},dt) -> T_x(dt/2)
//   -> C(dt/2)
// for the background and the symmetric Beam leapfrog (sections 7.1/7.2/13.9).
// The accepted state is only mutated at the end through swap_state, so a
// rejected step never changes the background, Beam, fields or ledgers.
class ConversionConsensusTestAccess;

// JC3 (section 6.10): test-only fault injection.  Never exposed through the
// production CLI; only the named test-access friend may set the flags.
struct FieldParticleJcFaultConfig {
    bool fail_final_poisson;     // forces failure 202 at the next trial
    bool fail_pairing_build;     // forces failure 203 at the next trial
    bool fail_post_field_charge; // forces failure 206 after convergence
    bool force_not_converged;    // exhausts max_iterations and returns 205
    int nan_inject_iteration;    // >=1: inject NaN into guess at this iteration
    int nan_inject_rank;         // -1 = all ranks, else only this rank
    int initial_guess_sign;      // 0 = default; +1/-1 scales the initial E_pair guess
    FieldParticleJcFaultConfig()
        : fail_final_poisson(false), fail_pairing_build(false),
          fail_post_field_charge(false), force_not_converged(false),
          nan_inject_iteration(0), nan_inject_rank(-1),
          initial_guess_sign(0)
    {}
};

class VpfpIntegrator {
public:
    VpfpIntegrator(const OpenBackgroundBoundary& background_boundary,
                   OpenElectrostaticSolver& field_solver,
                   const CylindricalFokkerPlanckCollision& collision);
    // Stage-H3 hybrid constructor: wires the tail PIC state and the
    // conservative converter into the same Strang split (sections 8.2 and
    // 14.7).  With background_tail_enabled=false this constructor behaves
    // exactly like the three-argument one (phase-4 path).
    VpfpIntegrator(const OpenBackgroundBoundary& background_boundary,
                   OpenElectrostaticSolver& field_solver,
                   const CylindricalFokkerPlanckCollision& collision,
                   const HybridVelocityPartition& partition,
                   BulkTailConverter& converter,
                   bool background_tail_enabled);

    void init(const SpatialGrid& grid);

    // Fresh-start flux-interface initialization.  Move any material mass in
    // tail-owned Eulerian cells through the same deterministic parcel loader
    // used by production conversion before the first Poisson solve.
    bool initialize_tail_from_bulk(
        Species& electrons, int mpi_rank, int mpi_size,
        BulkTailConversionDiagnostics& diagnostics);

    void set_beam_enabled(bool enabled);
    void set_background_phase_space_mode(BackgroundPhaseSpaceMode mode)
    { background_phase_space_mode_ = mode; }
    BackgroundPhaseSpaceMode background_phase_space_mode() const
    { return background_phase_space_mode_; }
    void set_x_transport_velocity_mode(XTransportVelocityMode mode)
    {
        x_transport_velocity_mode_ = mode;
        if (initialized_) vlasov_.set_x_transport_velocity_mode(mode);
    }
    XTransportVelocityMode x_transport_velocity_mode() const
    { return x_transport_velocity_mode_; }
    void set_tail_conversion_mode(TailConversionMode mode,
                                  int quadrature_order = 4,
                                  size_t max_supports = 7,
                                  int max_created_particles = 0);
    void set_tail_bulk_return(const TailBulkReturnConfig& config)
    { tail_bulk_return_.set_config(config); }
    const TailBulkReturnConfig& tail_bulk_return_config() const
    { return tail_bulk_return_.config(); }
    // Collision-face parcel export is not enabled until the conservative
    // collision-face path is implemented.  The only permitted interim mode
    // is an explicit zero-wall validation run; keep this gate in the
    // integrator as well as in the CLI so library callers cannot bypass it.
    void set_collision_interface_zero_wall_validation(bool enabled)
    { collision_interface_zero_wall_validation_ = enabled; }
    void set_collision_interface_exporting_absorbing(bool enabled)
    { collision_interface_exporting_absorbing_ = enabled; }
    TailConversionMode tail_conversion_mode() const
    { return tail_conversion_mode_; }
    // Section 11.5 resource gates: max particle count (0 = unlimited) and
    // max tail number fraction of the combined background.
    void set_tail_limits(std::uint64_t max_particles,
                         double max_number_fraction);
    // Stage-H5 population controller (section 7.10): configures the
    // deterministic conservative compression / equal-weight splitting that
    // runs on the trial tail after the final spatial half-drift.
    void set_tail_population_control(
        const TailPopulationController::Config& config);
    // Stage-H8: enables the unified HybridCollisionStep for the two Strang
    // collision halves (sections 10.4 and 14.11.1).  With tail_kernel ==
    // None the H7 bulk-only path is used.
    void set_tail_collision(TailCollisionKernel kernel, double coulomb_log,
                            TailCollisionWeightMode weight_mode,
                            int max_substeps,
                            double max_particle_growth);
    const HybridCollisionConfig& tail_collision_config() const
    {
        return hybrid_collision_config_;
    }
    // Stage-H6 checkpoint hooks (section 12): accepted-level cumulative
    // ledgers and combined checksums, plus restart-only state restore.
    const VpfpTailCumulativeLedger& tail_cumulative() const
    {
        return tail_cumulative_;
    }
    void restore_tail_cumulative(const VpfpTailCumulativeLedger& ledger);
    const VpfpCombinedChecksum& combined_checksum() const
    {
        return combined_checksum_;
    }
    const TailPopulationController::Config& population_control_config() const
    {
        return population_controller_.config();
    }
    // Restart only: the checkpoint carries the accepted step; the next
    // advance must continue the same step counter (conversion accepted_step
    // and population-control cadence depend on it).
    void set_step_count(long long count) { step_count_ = count; }
    long long step_count() const { return step_count_; }
    // Explicit opt-in liveness trace for long hybrid runs.  It writes only
    // local progress markers and never participates in physics or MPI gates.
    void set_tail_stage_trace(bool enabled, const std::string& output_dir);
    // Read-only stage audit.  It is intentionally independent of every
    // acceptance gate and is enabled only by diagnostic level 2 in main.
    void set_stage_energy_audit_enabled(bool enabled)
    { stage_energy_audit_enabled_ = enabled; }
    // Gate I (section 4.6.1): read-only field-particle power pairing audit,
    // also enabled only by diagnostic level 2.  level 0/1 keep the audit off
    // so no per-step I/O, full scan or collective is added to production.
    void set_field_particle_power_audit_enabled(bool enabled)
    { field_particle_power_audit_enabled_ = enabled; }

    // JC1 (section 4.4): field-particle coupling configuration.
    void set_field_particle_coupling(const FieldParticleCouplingConfig& config)
    { field_particle_coupling_config_ = config; }
    const FieldParticleCouplingConfig& field_particle_coupling_config() const
    { return field_particle_coupling_config_; }
    bool field_particle_coupling_enabled() const
    { return field_particle_coupling_config_.mode != FieldParticleCouplingMode::Legacy; }

    // Gate F (section 10.5.2): execute the production field-particle split
    // on caller-owned copies while freezing collision-driven representation
    // changes, bulk-to-tail conversion and H10 return.  The caller must use
    // a dedicated integrator instance for the transaction; its Tail state is
    // intentionally committed between substeps. advance() retains production
    // semantics.
    VpfpStepResult advance_fixed_state_field_particle_audit(
        Species& electrons, BeamPIC& beam, EMFields& fields,
        const std::vector<double>& ion_density, double time, double dt,
        int mpi_rank, int mpi_size);

    // Accepted tail state (read/write for main diagnostics and future
    // checkpoint wiring).
    const BackgroundTailPIC& tail_state() const { return tail_accepted_; }
    BackgroundTailPIC& tail_state() { return tail_accepted_; }
    bool background_tail_enabled() const { return background_tail_enabled_; }
    // Gate I (section 4.6): read-only access to the per-step pairing
    // workspace so the diagnostics layer can write the per-rank profile.
    const FieldParticlePowerAuditWorkspace& pairing_workspace() const
    { return pairing_workspace_; }

    // JC0 (section 3): post-field charge-invariance precheck gate.  Runs the
    // production C2 collision half, the collision-face conversion and the H10
    // tail-to-bulk return on caller-owned trial states and reports the
    // combined per-cell/rank/global electron number before and after.  This
    // is a read-only test hook: it reuses the exact production operators
    // (apply_collision_half, convert_flux_batch, apply_tail_bulk_return) and
    // never mutates the integrator's accepted state, cumulative ledgers or
    // RNG ownership.  Collective across all ranks.
    bool post_field_charge_invariance_transaction(
        Species& bulk_trial, BackgroundTailPIC& tail_trial,
        double time, double dt, int mpi_rank, int mpi_size,
        VpfpPostFieldChargeInvarianceReport& report);

    VpfpStepResult advance(Species& electrons, BeamPIC& beam, EMFields& fields,
                           const std::vector<double>& ion_density, double time,
                           double dt, int mpi_rank, int mpi_size);

private:
    friend class ConversionConsensusTestAccess;
    friend class FieldParticleJcTestAccess;
    void set_jc_fault_config(const FieldParticleJcFaultConfig& config)
    { jc_fault_ = config; }
    const FieldParticleJcFaultConfig& jc_fault_config() const { return jc_fault_; }
    FieldParticleJcFaultConfig jc_fault_;
    BackgroundPhaseSpaceMode background_phase_space_mode_;
    VlasovSplitStep vlasov_;
    XTransportVelocityMode x_transport_velocity_mode_;
    // Stage H7 (section 14.7 item 8): trial buffer for the first Strang
    // collision half-step.  A real bulk collision acts on a copy of the
    // accepted species; a rejected step never mutates the accepted state.
    Species state_collision_trial_;
    Species state_x_half_;
    Species state_u_full_;
    Species state_np1_;
    BeamPIC beam_work_;
    EMFields midpoint_fields_;
    EMFields final_fields_;
    BackgroundTailPIC tail_accepted_;
    BackgroundTailPIC tail_work_;
    // JC1 (section 4.8): field-particle trial buffers.
    BackgroundTailPIC tail_field_trial_;
    BeamPIC beam_field_trial_;
    EMFields field_n_pairing_;
    std::vector<double> pairing_field_guess_;
    std::vector<double> pairing_field_map_;
    std::vector<double> pairing_field_residual_;
    EMFields trial_force_fields_;
    FieldParticleCouplingConfig field_particle_coupling_config_;

    const OpenBackgroundBoundary& background_boundary_;
    OpenElectrostaticSolver& field_solver_;
    const CylindricalFokkerPlanckCollision& collision_;
    const HybridVelocityPartition* partition_;
    BulkTailConverter* converter_;
    bool background_tail_enabled_;
    bool fixed_state_field_particle_audit_mode_;
    std::uint64_t tail_max_particles_;
    double tail_max_number_fraction_;
    TailPopulationController population_controller_;
    HybridCollisionStep hybrid_collision_step_;
    HybridCollisionConfig hybrid_collision_config_;
    TailCollisionKernel tail_collision_kernel_;
    bool hybrid_collision_active_;
    // Per-velocity-slot bulk-collision ownership mask (1 = bulk-owned,
    // 0 = tail-owned conversion region), section 19.3.  Built from the
    // partition once tail collision is wired; empty when the tail is
    // disabled.
    std::vector<unsigned char> bulk_collision_mask_;
    VpfpTailCumulativeLedger tail_cumulative_;
    VpfpCombinedChecksum combined_checksum_;
    long long step_count_;
    bool beam_enabled_;
    TailConversionMode tail_conversion_mode_;
    bool collision_interface_zero_wall_validation_;
    bool collision_interface_exporting_absorbing_;
    int tail_flux_quadrature_order_;
    size_t tail_flux_max_supports_;
    int tail_flux_max_created_particles_per_step_;
    TailBulkReturn tail_bulk_return_;
    bool tail_stage_trace_enabled_;
    std::string tail_stage_trace_output_dir_;
    bool stage_energy_audit_enabled_;
    // JC3 discrete-gradient path: the stage-energy audit is split across the
    // frozen-state prepare, the Picard trial evaluation and the post-field
    // transaction.  These members carry the records and the cumulative
    // sources through those helpers so the accepted step can be finalized in
    // advance_discrete_gradient() exactly like the legacy paths.
    VpfpStepResult stage_energy_scratch_;
    VpfpStageEnergyRecord stage_sources_;
    VpfpStageEnergyRecord stage_sources_frozen_;
    bool field_particle_power_audit_enabled_;
    // JC4 (section 7.7): post-field charge residual from the last
    // discrete-gradient advance; zero for legacy runs.
    double post_field_charge_residual_linf_;
    FieldParticlePowerAuditWorkspace pairing_workspace_;
    FieldParticlePowerAudit pairing_calculator_;
    bool energy_ledger_initialized_;
    double accepted_field_energy_;
    double accepted_background_kinetic_energy_;
    double accepted_beam_kinetic_energy_;
    double accepted_tail_kinetic_energy_;
    bool initialized_;
    SpatialGrid grid_;

    bool finite_species(const Species& species) const;
    double field_energy(const EMFields& fields) const;
    double local_field_energy(const EMFields& fields) const;
    void capture_stage_energy(VpfpStepResult& result, int stage_id,
                              const Species& bulk,
                              const BackgroundTailPIC* tail,
                              const BeamPIC* beam,
                              const EMFields& fields,
                              const VpfpStageEnergyRecord& sources) const;
    void finalize_stage_energy_audit(VpfpStepResult& result) const;
    // JC3 discrete-gradient path: capture one stage-energy record into
    // stage_energy_scratch_ using the cumulative stage_sources_.  Mirrors the
    // legacy capture_stage lambda but reads/writes the member scratch so the
    // record survives across the prepare/trial/post-field helper split.
    void capture_dg_stage(int stage_id, const Species& bulk,
                          const BackgroundTailPIC* tail,
                          const BeamPIC* beam,
                          const EMFields& fields);
    // Gate I (section 4.6.1): fill the pairing workspace from the accepted
    // (n) and candidate (n+1) states, run the rank-local calculator and
    // perform the single packed reduction into result.pairing_audit.  Called
    // once before the transactional swap, when both states coexist.
    void finalize_field_particle_power_audit(
        VpfpStepResult& result,
        const Species& bulk_n, const Species& bulk_np1,
        const BeamPIC& beam_n, const BeamPIC& beam_np1,
        const BackgroundTailPIC* tail_n, const BackgroundTailPIC* tail_np1,
        const EMFields& field_n, const EMFields& field_np1,
        double bulk_work_global, double tail_work_local,
        double beam_work_local, double dt, int mpi_rank, int mpi_size);
    void initialize_energy_ledger(const Species& electrons,
                                  const BeamPIC& beam,
                                  const EMFields& fields,
                                  bool tail_on);
    void finalize_energy_ledger(VpfpStepResult& result);
    void swap_emfields(EMFields& a, EMFields& b);
    VpfpStepResult advance_background(Species& electrons, BeamPIC& beam,
                                       EMFields& fields,
                                       const std::vector<double>& ion_density,
                                       double time, double dt, int mpi_rank,
                                       int mpi_size);
    VpfpStepResult advance_with_beam(Species& electrons, BeamPIC& beam,
                                      EMFields& fields,
                                      const std::vector<double>& ion_density,
                                      double time, double dt, int mpi_rank,
                                      int mpi_size);
    // JC1 (section 4.4): discrete-gradient field-particle coupling path.
    VpfpStepResult advance_discrete_gradient(
        Species& electrons, BeamPIC& beam, EMFields& fields,
        const std::vector<double>& ion_density, double time, double dt,
        int mpi_rank, int mpi_size);
    VpfpStepResult advance_joint_midpoint(
        Species& electrons, BeamPIC& beam, EMFields& fields,
        const std::vector<double>& ion_density, double time, double dt,
        int mpi_rank, int mpi_size);
public:
    // JC2 (section 5.5.1): freeze the state for one physical step.
    // Executes C1, first x-half, tail drift, beam predict, midpoint Poisson
    // exactly once. All subsequent trial evaluations reuse this frozen state.
    bool prepare_field_particle_frozen_state(
        Species& electrons, BeamPIC& beam, EMFields& fields,
        const std::vector<double>& ion_density,
        double time, double dt, int mpi_rank, int mpi_size,
        FieldParticleFrozenState& frozen, VpfpFailureInfo& failure);
    // JC2 (section 5.5.2): evaluate one deterministic trial from frozen state.
    bool evaluate_field_particle_trial(
        const FieldParticleFrozenState& frozen,
        const std::vector<double>& pairing_field_guess,
        FieldParticleTrial& trial, VpfpFailureInfo& failure);
private:
    bool solve_field_particle_center(
        const FieldParticleFrozenState& frozen,
        FieldParticleTrial& accepted_trial,
        FieldParticleIterationDiagnostics& diagnostics,
        VpfpFailureInfo& failure);
    // JC3 (section 6.3/6.7): compute the global field fixed-point residuals
    // (R_E_L2, R_E_Linf) and the work pairing residual R_W with the fixed
    // two-collective MPI consensus.  Returns false only when the arrays are
    // structurally invalid; all_finite reports the LAND non-finite gate.
    bool compute_field_work_residuals(
        const std::vector<double>& pairing_guess,
        const FieldParticleTrial& trial,
        double& field_l2, double& field_linf, double& pairing_relative,
        bool& all_finite, int mpi_rank, int mpi_size) const;
    // JC3 (section 6.4/6.8): run C2, collision-face conversion and H10 return
    // exactly once on the converged trial, then validate the per-cell
    // combined charge against the JC0 tolerance.  On success the trial holds
    // the post-field candidate state and the converged final Poisson is
    // reused.  On charge violation returns failure code 206.
    bool apply_post_field_once_and_validate_charge(
        const FieldParticleFrozenState& frozen,
        FieldParticleTrial& trial, VpfpFailureInfo& failure);
    // Shared collision path for both production advance variants.  The
    // caller owns the trial copy; this routine only selects and executes the
    // configured bulk-only or hybrid pair backend on that trial.
    bool apply_collision_half(
        Species& bulk_trial, BackgroundTailPIC& tail_trial,
        double time, double half_dt, int collision_half, int mpi_rank,
        CollisionDiagnostics& bulk_diag,
        HybridCollisionDiagnostics& hybrid_diag, bool& hybrid_valid,
        BulkTailFluxBatch* collision_flux);
    bool append_collision_flux_parcels(
        const CollisionFaceFluxes& fluxes,
        BulkTailFluxBatch& batch,
        int collision_half,
        int mpi_rank) const;
    // Shared bulk-to-tail conversion hook required by section 7.11.17.6.
    // Both Beam and no-Beam advances call this exact transaction after the
    // midpoint kick and before the second spatial half-drift.
    bool apply_upar_flux_conversion(
        BulkTailFluxBatch& exported_flux,
        const BulkTailFluxBatch& first_collision_flux,
        const CollisionDiagnostics& first_collision,
        bool flux_mode,
        VpfpStepResult& result,
        int mpi_rank,
        int failure_code);
    bool apply_tail_bulk_return(Species& bulk_trial,
                                BackgroundTailPIC& tail_trial,
                                VpfpStepResult& result,
                                int mpi_rank, int mpi_size);
    // Rank-local conversion is transactional, but the following tail drift
    // contains MPI collectives.  Every rank must agree on conversion success
    // before any rank may enter that drift.
    bool synchronize_conversion_outcome(
        int local_failure_reason,
        const BulkTailConversionDiagnostics& conversion,
        int mpi_rank,
        int failure_code,
        VpfpStepResult& result) const;
    // Convert the two half-step diagnostics into one accepted, globally
    // reduced ledger.  Keeping this in one routine prevents Beam and
    // no-Beam runs from reporting different collision semantics.
    void finalize_collision_ledger(
        const CollisionDiagnostics& first_bulk,
        const CollisionDiagnostics& second_bulk,
        const HybridCollisionDiagnostics& first_hybrid,
        const HybridCollisionDiagnostics& second_hybrid,
        bool first_hybrid_valid, bool second_hybrid_valid,
        VpfpStepLedger& ledger) const;
    // Conversion is performed locally, while the accepted step file is a
    // rank-0 global ledger.  Cumulative checkpoint counters remain local and
    // are updated before this reduction.
    void globalize_conversion_ledger(VpfpStepLedger& ledger) const;
    void globalize_conversion_diagnostics(VpfpStepResult& result) const;
    void append_flux_face_audit(const BulkTailFluxBatch& batch,
                                VpfpStepResult& result) const;
    double tail_total_weight(const BackgroundTailPIC& tail) const;
    double tail_total_kinetic_energy(const BackgroundTailPIC& tail) const;
    void trace_tail_stage(int mpi_rank, double time,
                          const char* stage) const;
};

#endif
