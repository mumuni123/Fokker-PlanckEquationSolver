#ifndef FIELD_PARTICLE_POWER_AUDIT_H
#define FIELD_PARTICLE_POWER_AUDIT_H

// Gate I (section 4 of VPFP场粒子离散功同源闭合_根因判别与修复实施方案.md):
// read-only, rank-local field-particle discrete-power pairing audit.
//
// This header carries the raw audit structures, the stable root-cause mask
// and the dimension contract of section 4.1, plus the rank-local calculator
// used by VpfpIntegrator to reconstruct the four identities of section 4.5.
// The calculator performs no MPI and no state mutation; every identity is a
// pure function of its inputs.  VpfpIntegrator owns the per-step workspace
// (section 3.3) and performs the single packed reduction of section 4.6.1.
//
// Dimension contract (section 4.1, never changed by any test):
//   *_number_swept_face   face   m^-2   NOT divided by dt
//   current_face_x        face   A/m^2  already a full-step average
//   *_delta_ke_cell       cell   J/m^2  NOT divided by dt
//   E_pair                face   V/m    n/a
//   pairing residual      global/region  J/m^2  NOT divided by dt
//   residual power        global/region  J/(m^2 s)  divided by dt

#include "grid.h"
#include "maxwell.h"

#include <cstdint>
#include <vector>

// Stable root-cause mask (section 4.1): result files write the integer, the
// Python analyzer maps the integer back to a name.  Never reorder the bits.
enum PairingRootCauseMask {
    PAIRING_CAUSE_NONE       = 0,
    PAIRING_CAUSE_TRANSPORT  = 1 << 0,  // A: Poisson/transport pairing
    PAIRING_CAUSE_WORK       = 1 << 1,  // B: work-current vs charge-current
    PAIRING_CAUSE_TIME       = 1 << 2,  // C: time-centering of E_pair
    PAIRING_CAUSE_PIC        = 1 << 3,  // D: PIC deposit/gather adjoint
    PAIRING_CAUSE_CONVERSION = 1 << 4,  // E: conversion transfer accounting
    PAIRING_CAUSE_BOUNDARY   = 1 << 5   // F: boundary-region localization
};

// Rank-local per-substep x-face transport audit (section 4.2).  The face
// vectors hold nx_local+1 faces; the two physical boundary faces are indices
// 0 and nx_local.  A positive entry is transport in the +x direction.
struct XFaceTransportAudit {
    bool enabled;
    double substep_dt;
    std::vector<double> bulk_number_swept_face;
    std::vector<double> tail_number_swept_face;
    std::vector<double> beam_number_swept_face;

    XFaceTransportAudit() : enabled(false), substep_dt(0.0) {}
    void init(int nx_local) {
        bulk_number_swept_face.assign(static_cast<size_t>(nx_local) + 1, 0.0);
        tail_number_swept_face.assign(static_cast<size_t>(nx_local) + 1, 0.0);
        beam_number_swept_face.assign(static_cast<size_t>(nx_local) + 1, 0.0);
    }
};

// Rank-local per-step cell force work (section 4.4), J/m^2, nx_local cells.
// The two PIC boundary scalars collect the out-of-domain CIC shares of the
// kick work so that sum(cell) + boundary equals the global kick ledger
// (section I3).
struct CellForceWorkAudit {
    bool enabled;
    std::vector<double> bulk_delta_ke_cell;
    std::vector<double> tail_delta_ke_cell;
    std::vector<double> beam_delta_ke_cell;
    double tail_delta_ke_boundary;
    double beam_delta_ke_boundary;

    CellForceWorkAudit()
        : enabled(false), tail_delta_ke_boundary(0.0),
          beam_delta_ke_boundary(0.0) {}
    void init(int nx_local) {
        bulk_delta_ke_cell.assign(static_cast<size_t>(nx_local), 0.0);
        tail_delta_ke_cell.assign(static_cast<size_t>(nx_local), 0.0);
        beam_delta_ke_cell.assign(static_cast<size_t>(nx_local), 0.0);
        tail_delta_ke_boundary = 0.0;
        beam_delta_ke_boundary = 0.0;
    }
};

// Read-only PIC trajectory stage snapshot (section 4.3).  The two face
// vectors are A/m^2 snapshots of the production trajectory accumulator.
struct PicTrajectoryStageAudit {
    std::vector<double> after_first_drift_current_face;
    std::vector<double> after_second_drift_current_face;
    double left_boundary_number;
    double right_boundary_number;

    PicTrajectoryStageAudit()
        : left_boundary_number(0.0), right_boundary_number(0.0) {}
};

// The three candidate E_pair face fields of section 4.5.  Only the candidate
// that passes the mathematical derivation AND the manufactured test is the
// reported primary; the others are comparison-only.
enum EPairCandidate {
    EPAIR_ENDPOINT_AVERAGE = 0,   // 0.5*(E_f^n + E_f^{n+1})
    EPAIR_PRODUCTION_MIDPOINT = 1, // midpoint_fields_.Ex_face
    EPAIR_DISCRETE_GRADIENT = 2    // from Poisson potential-charge work
};

// Rank-local per-step raw workspace (section 3.3).  Owned by VpfpIntegrator;
// the calculator only reads it.
struct FieldParticlePowerAuditWorkspace {
    bool enabled;

    // Two bulk x half-step swept numbers (never added before finalize).
    XFaceTransportAudit bulk_x1;
    XFaceTransportAudit bulk_x2;
    // First/second drift trajectory-current snapshots (Tail and Beam).
    PicTrajectoryStageAudit tail_trajectory;
    PicTrajectoryStageAudit beam_trajectory;
    // Per-cell force work (all three species).
    CellForceWorkAudit cell_work;

    // Accepted (n) and candidate (n+1) cell-integrated number per cell
    // (m^-2), used only by the continuity identity.
    std::vector<double> bulk_number_n;
    std::vector<double> bulk_number_np1;
    // Cell-integrated Bulk number at the four sides of the two x-remaps.
    // Their differences isolate every accepted non-x local operator without
    // reconstructing or changing the production x current:
    //   pre_x1-n, pre_x2-post_x1, np1-post_x2.
    std::vector<double> bulk_number_pre_x1;
    std::vector<double> bulk_number_post_x1;
    std::vector<double> bulk_number_pre_x2;
    std::vector<double> bulk_number_post_x2;
    std::vector<double> tail_number_n;
    std::vector<double> tail_number_np1;
    std::vector<double> beam_number_n;
    std::vector<double> beam_number_np1;

    // Cell-integrated non-transport sources per species (m^-2).  Boundary
    // injection/removal and bulk-tail transfer are spatially resolved; a
    // single scalar cannot represent independent left/right events.
    std::vector<double> bulk_source_number_cell;
    std::vector<double> tail_source_number_cell;
    std::vector<double> beam_source_number_cell;

    // Per-cell source components (m^-2), section I2 breakdown.  Signed
    // "positive = mass added to the component at this cell"; the net source
    // equals injection + outflow + conversion + return + other.  The Bulk
    // `other` term is measured from accepted non-x stage transitions
    // (collision and velocity-boundary change).  Physical x-boundary fluxes
    // remain exclusively in the swept-face/trajectory-current arrays.
    std::vector<double> bulk_conversion_number_cell;   // bulk->tail removal (<0)
    std::vector<double> bulk_return_number_cell;       // tail->bulk addition (>0)
    std::vector<double> bulk_other_number_cell;        // collision/u-boundary remainder
    std::vector<double> tail_conversion_number_cell;   // bulk->tail addition (>0)
    std::vector<double> tail_return_number_cell;       // tail->bulk removal (<0)
    std::vector<double> tail_outflow_number_cell;      // open-boundary truncation (<0)
    std::vector<double> beam_injection_number_cell;    // Beam boundary injection (>0)

    // Face field states for the three E_pair candidates (V/m, nx_local+1).
    std::vector<double> field_n_ex_face;
    std::vector<double> field_mid_ex_face;
    std::vector<double> field_np1_ex_face;
    // Exact face dual of the paired finite-volume potential, built by the
    // read-only OpenElectrostaticSolver helper.  It is not an alias of the
    // endpoint field and must be present for Gate I4.
    std::vector<double> potential_pair_ex_face;

    // Full-step cell charge change (C/m^2 = rho_np1*dx - rho_n*dx) used by
    // the Poisson-transport pair.  Length nx_local.
    std::vector<double> cell_charge_change;

    // Energy terms of the reconstructed residual (section 4.5 identity 4).
    // These are rank-local raw sums; the integrator globalizes before
    // finalize when required by section 4.6.1.
    double field_energy_change;   // Delta U_E (J/m^2)
    double electrode_work;        // W_electrode (J/m^2)
    double bulk_delta_ke;         // Delta K_bulk (J/m^2)
    double tail_delta_ke;         // Delta K_tail (J/m^2)
    double beam_delta_ke;         // Delta K_beam (J/m^2)
    // Existing global ledgers the local work sums must match (section 4.4).
    double bulk_work_ledger;      // bulk_upar_face_work (J/m^2)
    double tail_work_ledger;      // tail_kick_work (J/m^2)
    double beam_work_ledger;      // beam_kick_work (J/m^2)
    double conversion_energy_removed; // R_conversion source (J/m^2)
    double boundary_energy_source;    // R_boundary source (J/m^2)
    double poisson_potential_charge_work; // production Gate-F identity term
    double poisson_identity_residual;     // diagnostic cross-check

    // §15.13.4 step 3: G* endpoint neighbours.  The face-alias G* of the
    // cell-average gather needs the adjacent rank's boundary cell value at
    // shared faces: gstar_left_neighbor_cell is the left neighbour's last
    // local cell and gstar_right_neighbor_cell the right neighbour's first
    // local cell.  The integrator fills them via Sendrecv before finalize;
    // they are unused when the corresponding face is the global endpoint.
    double gstar_left_neighbor_cell;
    double gstar_right_neighbor_cell;

    FieldParticlePowerAuditWorkspace()
        : enabled(false), field_energy_change(0.0), electrode_work(0.0),
          bulk_delta_ke(0.0), tail_delta_ke(0.0), beam_delta_ke(0.0),
          bulk_work_ledger(0.0), tail_work_ledger(0.0),
          beam_work_ledger(0.0), conversion_energy_removed(0.0),
          boundary_energy_source(0.0), poisson_potential_charge_work(0.0),
          poisson_identity_residual(0.0),
          gstar_left_neighbor_cell(0.0), gstar_right_neighbor_cell(0.0) {}
};

// Per-component continuity composition of the first bad cell (section I2).
// All quantities are cell-integrated number [m^-2] except the residual, which
// is the signed closure error of the same identity.
struct ContinuityComponentBreakdown {
    double delta_n;             // N^{n+1} - N^n
    double left_face_swept;     // Q_{i-1/2} (+x positive)
    double right_face_swept;    // Q_{i+1/2} (+x positive)
    double injection_source;    // boundary/reservoir/Beam injection (added)
    double outflow_source;      // open-boundary outflow (removed)
    double conversion_source;   // bulk<->tail conversion (signed by component)
    double return_source;       // H10 return (signed by component)
    double other_source;        // accepted local non-x source not listed above
    double residual;            // delta_n + (right-left) - net_source
    int bad_global_index;       // rank-local first-bad cell (calculator)
    int first_bad_rank;         // global owner rank (integrator, after reduce)

    ContinuityComponentBreakdown()
        : delta_n(0.0), left_face_swept(0.0), right_face_swept(0.0),
          injection_source(0.0), outflow_source(0.0),
          conversion_source(0.0), return_source(0.0), other_source(0.0),
          residual(0.0),
          bad_global_index(-1), first_bad_rank(-1) {}
};

// Result of one accepted step's pairing audit (section 4.1 plus the
// per-species/per-candidate detail required by sections 4.6/4.7.6/5.1).
struct FieldParticlePowerAuditResult {
    // Core fields (section 4.1).  Default construction sets the booleans to
    // false and the indices to -1 (section 4.1: no NaN as "not run" marker).
    bool valid;
    int first_bad_rank;
    int first_bad_index;
    double continuity_linf;
    double poisson_transport_residual;
    double force_work_residual;
    double current_pair_residual;
    double conversion_transfer_residual;
    double boundary_residual;
    double reconstructed_full_residual;
    double reconstruction_mismatch;
    double poisson_identity_crosscheck;
    double current_pair_linf;

    // Roundoff tolerance and pass flags (section 4.5 item 6).  All finite is
    // checked before any ordinary SUM reduction (section 4.6.1 MPI note 4).
    double roundoff_tolerance;
    // Independent scale/tolerance for the Gate-F Poisson identity crosscheck.
    // The crosscheck is the residual of Delta U_E - W_elec - W_potential_charge
    // over terms of magnitude ~1e3-1e4 J/m^2, so the reconstruction
    // roundoff_tolerance (scaled by the strongly-cancelled full_residual) is
    // not the correct gate.  See the section I6 fix.
    double poisson_identity_scale;
    double poisson_crosscheck_tolerance;
    double continuity_scale;
    double reconstruction_scale;
    bool all_finite;
    bool continuity_pass;
    bool local_work_sum_matches_existing_ledger;
    bool full_residual_reconstruction_pass;

    // Per-species detail (section 4.6 fixed columns and 5.1 item 8).
    double continuity_linf_bulk;
    double continuity_linf_tail;
    double continuity_linf_beam;
    double force_work_residual_bulk;
    double force_work_residual_tail;
    double force_work_residual_beam;

    // Per-component first-bad-cell continuity composition (section I2).  These
    // are rank-local (the rank owning the first bad cell); the integrator
    // resolves the deterministic global first-bad rank/index separately.
    ContinuityComponentBreakdown continuity_bulk;
    ContinuityComponentBreakdown continuity_tail;
    ContinuityComponentBreakdown continuity_beam;

    // Three E_pair candidates (section 4.5), reported so the analyzer can
    // discriminate A vs C (section 5.1 item 7) without re-running the solver.
    int validated_e_pair_candidate;
    double poisson_transport_residual_endpoint;
    double poisson_transport_residual_midpoint;
    double poisson_transport_residual_discrete_gradient;

    // A-F candidate contribution fractions F_c = sum|R_c| / sum|R_fp|
    // (section 5.1 item 5), computed by the calculator.
    double transport_fraction;
    double work_current_fraction;
    double time_center_fraction;
    double pic_fraction;
    double conversion_fraction;
    double boundary_fraction;

    // Stable root-cause mask (section 4.1).
    int root_cause_mask;

    // Rank-local raw sums that the integrator globalizes in its single packed
    // reduction (section 4.6.1): the actual full field-particle residual and
    // the three per-cell work sums.  The work-sum-vs-ledger comparison and
    // the reconstruction mismatch are nonlinear and are computed by the
    // integrator AFTER the reduction, not inside the rank-local calculator.
    double full_residual;
    double bulk_cell_work_sum;
    double tail_cell_work_sum;
    double beam_cell_work_sum;

    // Full physical step length used to convert the raw swept numbers and
    // cell work into currents/powers (section 4.2/4.1 dimension contract).
    double dt_s;

    // §15.13.4 step 3: direct face dual audit.  J_charge_face is the
    // charge-conserving face current (A/m^2, nx_local+1); J_force_cell is the
    // velocity-force cell current (A/m^2, nx_local) reconstructed from the
    // per-cell force work; Gstar_J_force_face is the face-aliased force
    // current (arithmetic mean of the two neighbouring cells, the discrete
    // adjoint of the cell-average gather); dual_face_residual is
    // J_charge_face - Gstar_J_force_face per face.  The three region sums are
    // the -dt <E_pair, dual_face_residual> integrals over left5/core90/right5.
    // These are read-only diagnostics and never modify the physical state.
    std::vector<double> j_charge_face;
    std::vector<double> j_force_cell;
    std::vector<double> gstar_j_force_face;
    std::vector<double> dual_face_residual;
    double dual_left5_integral;
    double dual_core90_integral;
    double dual_right5_integral;
    // Scalar reconstruction of the current-pair residual (section 15.13.6):
    //   current_pair_residual = dual_in_domain_work + boundary_force_work.
    // The first term is the face dual of the in-domain cell force work; the
    // second contains only Tail/Beam force work deposited outside the domain.
    double dual_in_domain_work;
    double boundary_force_work;
    double dual_plus_boundary_work;
    double dual_reconstruction_error;
    double dual_reconstruction_tolerance;
    bool dual_reconstruction_pass;

    FieldParticlePowerAuditResult()
        : valid(false), first_bad_rank(-1), first_bad_index(-1),
          continuity_linf(0.0), poisson_transport_residual(0.0),
          force_work_residual(0.0), current_pair_residual(0.0),
          conversion_transfer_residual(0.0), boundary_residual(0.0),
          reconstructed_full_residual(0.0), reconstruction_mismatch(0.0),
          poisson_identity_crosscheck(0.0), current_pair_linf(0.0),
          roundoff_tolerance(0.0), poisson_identity_scale(0.0),
          poisson_crosscheck_tolerance(0.0), continuity_scale(0.0),
          reconstruction_scale(0.0), all_finite(true),
          continuity_pass(false),
          local_work_sum_matches_existing_ledger(false),
          full_residual_reconstruction_pass(false),
          continuity_linf_bulk(0.0), continuity_linf_tail(0.0),
          continuity_linf_beam(0.0), force_work_residual_bulk(0.0),
          force_work_residual_tail(0.0), force_work_residual_beam(0.0),
          validated_e_pair_candidate(EPAIR_DISCRETE_GRADIENT),
          poisson_transport_residual_endpoint(0.0),
          poisson_transport_residual_midpoint(0.0),
          poisson_transport_residual_discrete_gradient(0.0),
          transport_fraction(0.0), work_current_fraction(0.0),
          time_center_fraction(0.0), pic_fraction(0.0),
          conversion_fraction(0.0), boundary_fraction(0.0),
          root_cause_mask(PAIRING_CAUSE_NONE),
          full_residual(0.0), bulk_cell_work_sum(0.0),
          tail_cell_work_sum(0.0), beam_cell_work_sum(0.0),
          dt_s(0.0),
          dual_left5_integral(0.0), dual_core90_integral(0.0),
          dual_right5_integral(0.0), dual_in_domain_work(0.0),
          boundary_force_work(0.0), dual_plus_boundary_work(0.0),
          dual_reconstruction_error(0.0), dual_reconstruction_tolerance(0.0),
          dual_reconstruction_pass(false) {}
};

// Rank-local pairing calculator (section 4.5).  It owns the fixed region
// classification (left 5% / core 90% / right 5%, section 4.6) and exposes one
// finalize entry point that follows the fixed order of section 4.6.1.  No MPI
// is performed here; VpfpIntegrator performs the single packed reduction.
class FieldParticlePowerAudit {
public:
    FieldParticlePowerAudit();
    void init(const SpatialGrid& grid);

    // Fixed-finalize-order calculator (section 4.6.1 pseudocode, numerical
    // portion).  All inputs are rank-local; the caller must have already
    // globalized the energy/ledger scalars and filled the workspace.  The
    // returned result is rank-local before the caller's single Allreduce.
    FieldParticlePowerAuditResult finalize(
        const FieldParticlePowerAuditWorkspace& ws,
        double dt, double electron_charge) const;

    // Region id for a local cell: 0 = left 5%, 1 = core 90%, 2 = right 5%.
    int region_id(int local_cell) const;
    int region_count() const { return 3; }

    // Manufactured discrete-adjoint check used by the read-only audit test.
    // For a complete single-rank domain this verifies
    // <E,G*J>_f = <GE,J>_c with the same shared-face owner and endpoint
    // quadrature weights as the production audit.
    bool manufactured_gstar_identity(const std::vector<double>& e_face,
                                     const std::vector<double>& j_cell,
                                     double& face_work,
                                     double& cell_work,
                                     double& absolute_error) const;

    bool manufactured_gstar_identity(
        const std::vector<double>& e_face,
        const std::vector<double>& j_cell,
        double left_neighbor_cell, double right_neighbor_cell,
        int mpi_rank, int mpi_size,
        double& face_work, double& cell_work,
        double& absolute_error) const;

    double face_quadrature_weight(int local_face) const;

    // Machine-scaled roundoff tolerance (section 5):
    //   T_round = max(1e-12*S_unit, 512*eps_double*S).
    static double roundoff_tolerance(double scale, double unit_floor);

private:
    SpatialGrid grid_;
    int left_region_end_;    // exclusive
    int right_region_begin_; // inclusive

    // Section 4.5 identity 1 per species: R_i = M_i^{n+1} - M_i^n
    //   - Q_{i-1/2} + Q_{i+1/2} - S_boundary,i, returning the linf residual
    // and the absolute contribution sum.  S_boundary is distributed to the
    // physical boundary cell only.
    double evaluate_species_continuity(
        const std::vector<double>& number_n,
        const std::vector<double>& number_np1,
        const std::vector<double>& swept_face,
        const std::vector<double>& source_number_cell,
        double& residual_linf) const;

    // Section 4.5: charge current from a raw swept-number face vector.
    // J_f = charge * Q_f / dt_full (section 4.2; the charge is the code's
    // signed electron charge, -qe).
    void swept_number_to_current(const std::vector<double>& swept,
                                 double dt, double charge,
                                 std::vector<double>& current) const;

    // Section 4.5: discrete face inner product sum_i E_i * J_i * dx_i with
    // the shared-face quadrature weight already folded into dx_i by the
    // caller (open boundary faces use the boundary weight, section 4.5).
    double face_inner_product(const std::vector<double>& e_face,
                              const std::vector<double>& j_face) const;

    // Poisson-transport residual R_{P<->J} = (Delta U_E - W_electrode)
    //   + dt <E_pair, J_charge> for one candidate field.
    double poisson_transport_residual_candidate(
        double field_energy_change, double electrode_work, double dt,
        const std::vector<double>& e_pair,
        const std::vector<double>& j_charge) const;

    // Build the three candidate face fields (section 4.5) into out[3].
    void build_e_pair_candidates(const FieldParticlePowerAuditWorkspace& ws,
                                 std::vector<double> out[3]) const;

    // Classify the A-F root-cause mask (section 5) from the reconstructed
    // residual terms; returns the mask and fills the six fractions.
    int classify_root_cause_mask(const FieldParticlePowerAuditResult& r,
                                 double full_abs_sum) const;

    // §15.13.4 step 3: region id (0=left5, 1=core90, 2=right5) of a local
    // face, defined by its global face index.
    int region_id_for_face(int local_face) const;
};

#endif
