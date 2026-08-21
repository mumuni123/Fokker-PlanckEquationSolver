#ifndef TAIL_INTERFACE_GRID_DESIGN_H
#define TAIL_INTERFACE_GRID_DESIGN_H

// Section 7.11.16C / 16B: fixed-memory local velocity-grid offline
// prototype.  The module is deliberately independent: it never modifies
// CylindricalVelocityGrid::init() (the production default), never touches
// main_vpfp.cpp global CLI state and never alters Param::Nv/Param::Nmu.
// It builds the four fixed-cell candidates G0/Gx/Gp/G2 by deterministic
// monitor equidistribution on the baseline grid, conservatively remaps the
// real request histogram with the exact cylindrical overlap rule, and
// reports the representation, feasibility, memory and step-size budgets
// that the 16B offline gates consume.

#include "bulk_tail_moment_audit.h"
#include "grid.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

// Monitor-amplitude / width configuration for one candidate.  All values
// are CLI-provided in the scan; nothing is hard-coded into the production
// grid.  min_width_ratio and max_adjacent_width_ratio are the constraint
// limits applied after face generation (defaults 0.5 and 2.0).
struct TailInterfaceGridDesignConfig {
    double ax;
    double aperp;
    double sigma_x_cells;
    double sigma_perp_cells;
    double min_width_ratio;
    double max_adjacent_width_ratio;

    TailInterfaceGridDesignConfig()
        : ax(0.0), aperp(0.0), sigma_x_cells(1.0), sigma_perp_cells(1.0),
          min_width_ratio(0.5), max_adjacent_width_ratio(2.0)
    {}
};

// One generated candidate.  grid_name is "G0"/"Gx"/"Gp"/"G2" with the scan
// parameters appended for the monitor profiles; faces are strictly
// monotone, u_parallel mirror symmetric with u_parallel=0 fixed, u_perp=0
// fixed and both endpoints fixed.  A candidate that violates the adjacent
// width ratio or the minimum width constraint is marked valid=0 and is
// never clipped.
struct TailInterfaceGridCandidate {
    std::string grid_name;
    std::vector<double> upar_faces;
    std::vector<double> uperp_faces;
    bool valid;
    std::string invalid_reason;

    TailInterfaceGridCandidate()
        : valid(false), invalid_reason("unbuilt")
    {}
};

// Six-moment component order matches BulkTailMomentComponent
// (N, Px, Jx, K, Pixx, Piperp) so the replay reuses the production
// mass_cell_moments / TailSubcellQuadrature formulas.
struct TailInterfaceReplayResult {
    // Conservation and threshold identity.
    bool g0_identity_ok;
    double number_residual;                  // |sum' - sum| / max(sum, floor)
    double max_partition_error;              // max_jk |sum(theta)-1|
    long long negative_mass_cells;
    double max_negative_mass;
    double below_threshold_number_fraction;

    // Cumulative six-moment representation metrics over the remapped mass.
    std::array<double, BULK_TAIL_MOMENT_COUNT> center_total;
    std::array<double, BULK_TAIL_MOMENT_COUNT> volume_total;
    std::array<double, BULK_TAIL_MOMENT_COUNT> delta_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> center_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> r_l1;
    std::array<double, BULK_TAIL_MOMENT_COUNT> max_cell_relative;
    std::array<double, BULK_TAIL_MOMENT_COUNT> important_max_cell_relative;

    // Grid geometry and budget estimates.
    double min_dupar;
    double min_duperp;
    double max_adjacent_width_ratio;
    double estimated_cell_count_ratio;       // strictly 1
    double estimated_memory_ratio;           // strictly 1
    double estimated_velocity_dt_ratio;      // dt_candidate / dt_current
    double velocity_dt_upar_current;
    double velocity_dt_upar_candidate;
    double velocity_dt_uperp_current;
    double velocity_dt_uperp_candidate;
    double velocity_dt_total_current;
    double velocity_dt_total_candidate;
    double estimated_operator_work_ratio;
    double scan_cost_ratio;
    double pic_creation_ratio;
    double collision_pair_ratio;
    double collision_substep_ratio;

    // Volume-self feasibility and deterministic sparse support reduction.
    long long center_target_feasible_count;
    long long center_target_failed_count;
    long long volume_self_feasible_count;
    long long volume_self_failed_count;
    long long volume_self_sparse_failed_count;
    int max_sparse_support_count;
    long long estimated_created_macroparticles;
    double estimated_particle_ratio_to_center_quartet;
    std::array<double, BULK_TAIL_MOMENT_COUNT> scale_floor_max;

    // Classification: GREEN (all 16B gates pass), GRAY (hard gates pass but
    // the six-moment / important-cell physics gates stay in the gray band),
    // RED (a hard gate fails), INVALID (constraint-violated grid or the G0
    // identity check failed).
    std::string candidate_status;
    std::string status_reason;

    TailInterfaceReplayResult();
};

// Build one of the four candidates.  profile is "G0" (baseline faces),
// "Gx" (u_parallel monitor only, A_perp=0), "Gp" (u_perp monitor only,
// A_x=0) or "G2" (both).  u* is the mass-weighted |u_parallel| of the
// low-u_perp conversion cells of the real histogram (deterministic, never
// hard-coded to a specific iv).  Identical parameters and input must
// produce bit-identical faces.
TailInterfaceGridCandidate build_tail_interface_grid_candidate(
    const CylindricalVelocityGrid& baseline,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    const TailInterfaceGridDesignConfig& config,
    const std::string& profile);

// Replay the request histogram on the candidate grid: conservative
// cylindrical overlap remap, cumulative center/volume six moments, volume-
// self feasibility with the deterministic <=7-node sparse reduction, and
// the production velocity-dt / operator-work budget model.  conversion_energy
// (J) is the tail threshold K_out used for the below-threshold fraction.
TailInterfaceReplayResult replay_tail_interface_histogram(
    const CylindricalVelocityGrid& baseline,
    const TailInterfaceGridCandidate& candidate,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    double conversion_energy);

// Single production cylindrical conservative remap rule (section 7.11.16C):
// theta = (du_par overlap / from width) *
//         ((u_perp_hi^2 - u_perp_lo^2) / (from cell u_perp^2 range)).
// from_mass is indexed j*Nmu+k on the from faces; the result is indexed on
// the to faces (same cell count).  max_partition_error_out (optional)
// receives max_jk |sum_theta - 1|.
std::vector<double> tail_interface_remap_masses(
    const std::vector<double>& from_upar_faces,
    const std::vector<double>& from_uperp_faces,
    const std::vector<double>& to_upar_faces,
    const std::vector<double>& to_uperp_faces,
    const std::vector<double>& from_mass,
    double* max_partition_error_out = NULL);

// Aggregate the histogram rows by (iv, imu) into a j*Nmu+k cell-mass array.
// The per-event (accepted_step, iv, imu) rows of the audit file are summed
// into the cell mass the offline prototype consumes.
std::vector<double> tail_interface_aggregate_histogram(
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    int nv, int nmu);

#endif
