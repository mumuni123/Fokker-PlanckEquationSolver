#include "vpfp_diagnostics.h"

#include "beam_pic.h"
#include "bulk_tail_moment_audit_io.h"
#include "species.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <sstream>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {
const char* stage_energy_name(int stage_id)
{
    static const char* const names[VPFP_STAGE_ENERGY_RECORD_COUNT] = {
        "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
        "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
        "collision_half2", "conversion_after_collision", "tail_bulk_return",
        "final_poisson" };
    return stage_id >= 0 && stage_id < VPFP_STAGE_ENERGY_RECORD_COUNT
        ? names[stage_id] : "invalid";
}

bool create_directory_tree(const std::string& path)
{
    if (path.empty()) return false;
    std::string current;
    for (size_t i = 0; i < path.size(); ++i) {
        current.push_back(path[i]);
        const bool separator = path[i] == '/' || path[i] == '\\';
        if (!separator || current.size() == 1) continue;
#ifdef _WIN32
        if (_mkdir(current.c_str()) != 0 && errno != EEXIST) return false;
#else
        if (mkdir(current.c_str(), 0755) != 0 && errno != EEXIST) return false;
#endif
    }
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}
}

VpfpDiagnostics::VpfpDiagnostics()
    : level_(1), interval_(500), last_control_groups_(0),
      last_control_fallbacks_(0), moment_audit_headers_written_(false),
      moment_audit_hist_headers_written_(false),
      pairing_headers_written_(false),
      pairing_breakdown_headers_written_(false),
      field_particle_iteration_headers_written_(false)
{}

bool VpfpDiagnostics::init(const std::string& output_dir, int rank,
                           int level, int interval)
{
    output_dir_ = output_dir;
    level_ = level;
    interval_ = interval > 0 ? interval : 500;
    int root_ready = 1;
    if (rank == 0) {
        root_ready = create_directory_tree(output_dir_) ? 1 : 0;
        if (root_ready) {
            std::ofstream probe((output_dir_ + "/.vpfp_write_probe").c_str());
            root_ready = probe ? 1 : 0;
            probe.close();
            if (root_ready) std::remove((output_dir_ + "/.vpfp_write_probe").c_str());
        }
    }
    MPI_Bcast(&root_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!root_ready) return false;
    int headers_ready = 1;
    if (rank == 0) {
    std::ofstream out((output_dir_ + "/vpfp_step_diagnostics.dat").c_str());
    if (!out) headers_ready = 0;
    out << "step time_s accepted split remap_const_fraction "
        << "remap_linear_fraction gauss_linf gauss_charge_residual "
        << "N_e_before N_e_after N_b_before N_b_after "
        << "background_left_flux background_right_flux "
        << "background_left_inflow_energy background_left_outflow_energy "
        << "background_right_inflow_energy background_right_outflow_energy "
        << "beam_injected beam_outflow beam_injected_energy beam_outflow_energy "
        << "U_E K_e K_b fct_energy collision_reservoir "
        << "U_E_before K_e_before K_b_before K_tail_before "
        << "electrostatic_boundary_work background_boundary_energy_net "
        << "beam_boundary_energy_net domain_energy_before domain_energy_after "
        << "domain_energy_delta accounted_energy_source "
        << "energy_balance_residual energy_balance_relative "
        << "U_E_pair U_E_final "
        << "wall_s vlasov_s field_s beam_s collision_s mpi_collective_s "
        << "max_rss_kib tail_particle_count "
        << "wall_tail_push_s wall_tail_deposit_s wall_tail_migrate_s "
        << "wall_conversion_s wall_diagnostics_s tail_particles_local_max "
        << "P_bkg P_tail P_combined K_tail K_combined "
        << "N_tail_before N_tail_after N_combined_after tail_number_balance_error "
        << "tail_outflow_N tail_outflow_K tail_return_enabled "
        << "tail_return_candidate_particles tail_return_resident_particles "
        << "tail_return_attempted_groups tail_return_committed_groups "
        << "tail_return_deferred_infeasible_groups "
        << "tail_return_deferred_rank_boundary_groups "
        << "tail_return_projection_invalid_input_cells "
        << "tail_return_projection_insufficient_support_cells "
        << "tail_return_projection_infeasible_invariant_cells "
        << "tail_return_projection_representation_incompatible_cells "
        << "tail_return_particles_removed tail_return_N tail_return_Px "
        << "tail_return_Jx_dx tail_return_K tail_return_Pixx_dx "
        << "tail_return_Piperp_dx tail_return_N_residual "
        << "tail_return_Px_residual tail_return_Jx_residual "
        << "tail_return_K_residual tail_return_Pixx_residual "
        << "tail_return_Piperp_residual "
        << "tail_return_N_difference tail_return_Px_difference "
        << "tail_return_Jx_dx_difference tail_return_K_difference "
        << "tail_return_Pixx_dx_difference tail_return_Piperp_dx_difference "
        << "tail_return_mpi_request_residual "
        << "tail_return_wall_seconds "
        << "conversion_N conversion_Px conversion_K "
        << "conversion_N_residual conversion_Px_residual "
        << "conversion_K_residual velocity_boundary_loss_N "
        << "velocity_boundary_loss_K "
        << "Cbb_dK Ctt_dPx Ctt_dK Ctb_dPx Ctb_dK "
        << "Cbt_dPx Cbt_dK Cbt_P_residual Cbt_K_residual "
        << "pair_bb pair_tt pair_tb pair_bt "
        << "collision_flux_export_N collision_flux_export_K "
        << "collision_flux_implicit_residual collision_flux_cross_residual "
        << "collision_flux_inward_clipped_N collision_flux_parcel_count "
        << "collision_flux_rollback_count "
        << "x_transport_velocity_mode energy_conjugate_vx_min "
        << "energy_conjugate_vx_max energy_conjugate_vs_analytic_linf "
        << "energy_conjugate_vs_analytic_l2 "
        << "background_phase_space_mode joint_midpoint_converged "
        << "joint_midpoint_iterations joint_midpoint_residual_linf "
        << "joint_midpoint_poisson_residual_linf joint_midpoint_energy_residual\n";
    if (level_ >= 2) {
        std::ofstream conversion_out(
            (output_dir_ + "/conversion_source_accepted_steps.dat").c_str());
        conversion_out << "step time_s location bin_low_J bin_high_J "
                       << "pre_bulk_N pre_bulk_K removed_bulk_N removed_bulk_K "
                       << "created_tail_N created_tail_K "
                       << "accepted_tail_total_N accepted_tail_total_K\n";
        std::ofstream stage_out(
            (output_dir_ + "/vpfp_stage_energy_audit.dat").c_str());
        stage_out << "step time_s accepted audit_valid split failure_code "
                  << "energy_balance_residual "
                  << "stage_id stage_name "
                  << "K_bulk K_tail K_beam U_E "
                  << "dK_bulk dK_tail dK_beam dU_E "
                  << "Q_bkg_left_in Q_bkg_left_out Q_bkg_right_in Q_bkg_right_out "
                  << "Q_beam_in Q_beam_out Q_tail_out Q_collision_reservoir "
                  << "K_conversion K_tail_return W_electrostatic_boundary "
                  << "stage_balance "
                  << "bulk_upar_face_work bulk_upar_velocity_boundary_work "
                  << "bulk_upar_interface_energy_removed "
                  << "bulk_upar_identity_residual "
                  << "tail_kick_work beam_kick_work\n";
        // Create the scalar pairing ledger at startup.  Rows remain strictly
        // accepted-step-only, but an empty file no longer ambiguously means
        // either "first step not accepted yet" or "audit was disabled".
        std::ofstream pairing_out(
            (output_dir_ + "/field_particle_power_pairing.dat").c_str());
        pairing_out << "# schema=vpfp-field-particle-power-pairing-v3\n";
        pairing_out << "# columns=step time_s dt_s accepted "
            << "continuity_bulk continuity_tail continuity_beam "
            << "poisson_transport_residual "
            << "poisson_endpoint poisson_midpoint poisson_discrete_gradient "
            << "force_work_bulk force_work_tail force_work_beam "
            << "current_pair_residual conversion_residual boundary_residual "
            << "full_residual reconstructed_residual reconstruction_mismatch "
            << "poisson_identity_crosscheck poisson_identity_scale "
            << "poisson_crosscheck_tolerance current_pair_linf "
            << "roundoff_tolerance all_finite continuity_pass "
            << "local_work_ledger_pass reconstruction_pass root_cause_mask "
            << "first_bad_rank first_bad_index "
            << "transport_fraction work_current_fraction time_center_fraction "
            << "pic_fraction conversion_fraction boundary_fraction "
                 << "dual_left5_integral dual_core90_integral dual_right5_integral "
                 << "dual_in_domain_work boundary_force_work "
                 << "dual_plus_boundary_work dual_reconstruction_error "
                 << "dual_reconstruction_tolerance dual_reconstruction_pass "
                 << "dual_total_integral\n";
        pairing_headers_written_ = static_cast<bool>(pairing_out);
        if (!conversion_out || !stage_out || !pairing_out) headers_ready = 0;
    }
    std::ofstream flux_out(
        (output_dir_ + "/bulk_tail_flux_accepted_steps.dat").c_str());
    flux_out << "accepted_step time_s conversion_mode quadrature_order "
             << "max_supports parcel_count node_count exported_number exported_energy "
             << "face_export_number parcel_number below_threshold_number "
             << "roundoff_discarded_number duplicate_count "
             << "quadrature_error_max "
             << "tail_owned_expected_transfer_number "
             << "tail_owned_roundoff_discarded_number tail_owned_bulk_residual "
             << "static_extractor_call_count "
             << "particles_created compression_fallback_count "
             << "subcell_fallback_count support_limit_violation_count "
             << "duplicate_id_count face_ledger_mismatch_count "
             << "face_audit_count face_audit_face_abs_sum "
             << "face_audit_parcel_abs_sum face_audit_abs_error_sum "
             << "face_audit_max_relative "
             << "face_audit_abs_at_max_relative "
             << "conversion_wall_seconds "
             << "conversion_number conversion_energy "
             << "conversion_number_residual conversion_px_residual "
             << "conversion_energy_residual conversion_jx_residual "
             << "conversion_pixx_residual conversion_piperp_residual "
             << "conversion_rho_l2 "
             << "audit_parcel_failure_reason audit_parcel_failure_rank "
             << "audit_parcel_failure_ix audit_parcel_failure_face "
             << "audit_parcel_failure_iuperp audit_parcel_failure_node_mass "
             << "audit_parcel_failure_target audit_parcel_failure_node_sum "
             << "audit_parcel_failure_scale "
             << "audit_valid audit_failure_code "
             << "audit_inplace_state_bitwise_equal "
             << "audit_inplace_rng_equal audit_inplace_ledger_equal\n";
    std::ofstream flux_face_out(
        (output_dir_ + "/bulk_tail_flux_face_accepted_steps.dat").c_str());
    flux_face_out
        << "accepted_step time_s ix_global direction face_index transverse_index "
        << "operator_stage face_export_number parcel_number "
        << "absolute_error relative_error node_failure_reason "
        << "reconstructed_target node_sum\n";
    std::ofstream collision_flux_out(
        (output_dir_ + "/collision_flux_accepted_steps.dat").c_str());
    collision_flux_out
        << "accepted_step time_s exported_number exported_energy "
        << "implicit_flux_residual_linf cross_flux_pair_residual_linf "
        << "inward_clipped_number parcel_count rollback_count\n";
    std::ofstream collision_substep_out(
        (output_dir_ + "/collision_flux_substeps_accepted.dat").c_str());
    collision_substep_out
        << "accepted_step time_s collision_half exported_number "
        << "exported_energy implicit_flux_residual_linf "
        << "cross_flux_pair_residual_linf inward_clipped_number "
        << "parcel_count rollback_count\n";
    headers_ready = headers_ready && static_cast<bool>(out) && static_cast<bool>(flux_out) &&
        static_cast<bool>(flux_face_out) && static_cast<bool>(collision_flux_out) &&
        static_cast<bool>(collision_substep_out);
    }
    MPI_Bcast(&headers_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return headers_ready != 0;
}

void VpfpDiagnostics::write_accepted_step(int step, double time,
                                           const VpfpStepResult& r,
                                           const OpenGaussDiagnostics& g,
                                           const Species& electrons,
                                           const BackgroundTailPIC* tail,
                                           const SpatialGrid& grid,
                                           int rank,
                                           const FieldParticlePowerAuditWorkspace*
                                               pairing_ws)
{
    double energy_vx_min = std::numeric_limits<double>::quiet_NaN();
    double energy_vx_max = std::numeric_limits<double>::quiet_NaN();
    double energy_vx_linf = 0.0;
    long double energy_vx_l2_sum = 0.0L;
    size_t energy_vx_count = 0;
    const bool collect_velocity_table_diagnostic = level_ >= 2 ||
        (interval_ > 0 && step % interval_ == 0) || r.tail_return_enabled;
    if (collect_velocity_table_diagnostic) {
        const std::vector<double>& energy_vx =
            electrons.cgrid.vx_energy_conjugate_cell;
        const std::vector<double>& analytic_vx = electrons.cgrid.vx;
        energy_vx_count = std::min(energy_vx.size(), analytic_vx.size());
        if (energy_vx_count > 0) {
            energy_vx_min = std::numeric_limits<double>::infinity();
            energy_vx_max = -std::numeric_limits<double>::infinity();
        }
        for (size_t q = 0; q < energy_vx_count; ++q) {
            energy_vx_min = std::min(energy_vx_min, energy_vx[q]);
            energy_vx_max = std::max(energy_vx_max, energy_vx[q]);
            const double difference = energy_vx[q] - analytic_vx[q];
            energy_vx_linf = std::max(energy_vx_linf, std::fabs(difference));
            energy_vx_l2_sum += static_cast<long double>(difference) * difference;
        }
    }
    const double energy_vx_l2 = energy_vx_count == 0 ?
        std::numeric_limits<double>::quiet_NaN() :
        std::sqrt(static_cast<double>(energy_vx_l2_sum /
                                      static_cast<long double>(energy_vx_count)));
    // This audit has accepted-only semantics and must not inherit the normal
    // diagnostic interval.  Every rank enters its collectives for an
    // accepted event before rank 0 decides whether to write regular output.
    write_bulk_tail_moment_audit_accepted_step(step, time, r, rank);
    write_bulk_tail_flux_accepted_step(step, time, r, rank);
    // JC4 (section 7.3): field-particle iteration accepted-step diagnostic,
    // written for all diagnostic levels (0/1/2).  Rank 0 only.
    write_field_particle_iteration_accepted_step(step, time, r, rank);
    write_joint_midpoint_iterations(step, time, r, rank);
    if (level_ >= 2) {
        write_stage_energy_audit_accepted_step(step, time, r, rank);
        write_field_particle_power_pairing_accepted_step(
            step, time, r, rank, grid, pairing_ws);
    }
    if (level_ <= 1 && step % interval_ != 0 && !r.tail_return_enabled) return;
    double local_px[2] = { 0.0, 0.0 };
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int sx = grid.nghost + ix;
        const size_t xbase = static_cast<size_t>(sx) * Param::Nvmu;
        for (int j = 0; j < Param::Nv; ++j) {
            const double px_per_particle =
                Const::me * Const::c * electrons.cgrid.upar_cells[j];
            const size_t row = xbase + static_cast<size_t>(j) * Param::Nmu;
            for (int k = 0; k < Param::Nmu; ++k) {
                local_px[0] += electrons.f[row + static_cast<size_t>(k)] *
                               px_per_particle;
            }
        }
    }
    if (tail != NULL) {
        for (size_t p = 0; p < tail->particles.size(); ++p) {
            const BackgroundTailParticle& particle = tail->particles[p];
            local_px[1] += Const::me * Const::c * particle.weight *
                           particle.ux;
        }
    }
    double global_px[2] = {};
    MPI_Allreduce(local_px, global_px, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    // This routine contains collectives; all ranks must enter before the
    // rank-0-only step file return below.
    write_conversion_source_ledger(step, time, r, tail, rank);
    if (rank != 0) return;
    std::ofstream out((output_dir_ + "/vpfp_step_diagnostics.dat").c_str(),
                      std::ios::app);
    out << std::setprecision(17) << step << " " << time << " 1 "
        << (r.split_used ? 1 : 0) << " "
        << r.remap_constant_fraction << " " << r.remap_linear_fraction << " "
        << g.residual_linf << " "
        << g.boundary_charge_residual << " "
        << r.ledger.background_number_before << " " << r.ledger.background_number_after << " "
        << r.ledger.beam_number_before << " " << r.ledger.beam_number_after << " "
        << r.ledger.background_left_flux << " " << r.ledger.background_right_flux << " "
        << r.ledger.background_left_inflow_energy << " "
        << r.ledger.background_left_outflow_energy << " "
        << r.ledger.background_right_inflow_energy << " "
        << r.ledger.background_right_outflow_energy << " "
        << r.ledger.beam_injected << " " << r.ledger.beam_outflow << " "
        << r.ledger.beam_injected_energy << " "
        << r.ledger.beam_outflow_energy << " "
        << r.ledger.field_energy << " " << r.ledger.background_kinetic_energy << " "
        << r.ledger.beam_kinetic_energy << " " << r.ledger.fct_energy_change << " "
        << r.ledger.collision_reservoir_energy << " "
        << r.ledger.field_energy_before << " "
        << r.ledger.background_kinetic_energy_before << " "
        << r.ledger.beam_kinetic_energy_before << " "
        << r.ledger.tail_kinetic_energy_before << " "
        << r.ledger.electrostatic_boundary_work << " "
        << r.ledger.background_boundary_energy_net << " "
        << r.ledger.beam_boundary_energy_net << " "
        << r.ledger.domain_energy_before << " "
        << r.ledger.domain_energy_after << " "
        << r.ledger.domain_energy_change << " "
        << r.ledger.accounted_energy_source << " "
        << r.ledger.energy_balance_residual << " "
        << r.ledger.energy_balance_relative << " "
        << r.ledger.e_pair_field_energy << " "
        << r.ledger.field_energy << " "
        << r.wall_seconds_per_step << " "
        << r.wall_vlasov_seconds << " "
        << r.wall_field_seconds << " " << r.wall_beam_seconds << " "
        << r.wall_collision_seconds << " " << r.mpi_collective_seconds << " "
        << r.max_rss_kib << " " << r.ledger.tail_particle_count_after << " "
        << r.wall_tail_push_seconds << " " << r.wall_tail_deposit_seconds
        << " " << r.wall_tail_migrate_seconds << " "
        << r.wall_conversion_seconds << " " << r.wall_diagnostics_seconds
        << " " << r.tail_particles_local_max
        << " " << global_px[0] << " " << global_px[1]
        << " " << global_px[0] + global_px[1]
        << " " << r.ledger.tail_kinetic_energy_after
        << " " << r.ledger.background_kinetic_energy +
                       r.ledger.tail_kinetic_energy_after
        << " " << r.ledger.tail_number_before
        << " " << r.ledger.tail_number_after
        << " " << r.ledger.combined_number_after
        << " " << r.ledger.tail_number_balance_error
        << " " << r.ledger.tail_outflow_number
        << " " << r.ledger.tail_outflow_energy
        << " " << (r.tail_return_enabled ? 1 : 0)
        << " " << r.tail_return.candidate_particles
        << " " << r.tail_return.resident_particles
        << " " << r.tail_return.attempted_groups
        << " " << r.tail_return.committed_groups
        << " " << r.tail_return.deferred_infeasible_groups
        << " " << r.tail_return.deferred_rank_boundary_groups
        << " " << r.tail_return.projection_invalid_input_cells
        << " " << r.tail_return.projection_insufficient_support_cells
        << " " << r.tail_return.projection_infeasible_invariant_cells
        << " " << r.tail_return.projection_representation_incompatible_cells
        << " " << r.tail_return.particles_removed
        << " " << r.tail_return.number
        << " " << r.tail_return.px
        << " " << r.tail_return.jx_dx
        << " " << r.tail_return.energy
        << " " << r.tail_return.pixx_dx
        << " " << r.tail_return.piperp_dx
        << " " << r.tail_return.number_residual
        << " " << r.tail_return.px_residual
        << " " << r.tail_return.jx_residual
        << " " << r.tail_return.energy_residual
        << " " << r.tail_return.pixx_residual
        << " " << r.tail_return.piperp_residual
        << " " << r.tail_return.number_difference
        << " " << r.tail_return.px_difference
        << " " << r.tail_return.jx_difference
        << " " << r.tail_return.energy_difference
        << " " << r.tail_return.pixx_difference
        << " " << r.tail_return.piperp_difference
        << " " << r.tail_return.mpi_request_residual
        << " " << r.tail_return.wall_seconds
        << " " << r.ledger.conversion_number_removed
        << " " << r.ledger.conversion_px_removed
        << " " << r.ledger.conversion_energy_removed
        << " " << r.ledger.conversion_number_residual_rel
        << " " << r.ledger.conversion_px_residual_rel
        << " " << r.ledger.conversion_energy_residual_rel
        << " " << r.ledger.background_tail_number_loss
        << " " << r.ledger.background_tail_energy_loss
        << " " << r.ledger.collision_bulk_bulk_energy_change
        << " " << r.ledger.collision_tail_tail_px_change
        << " " << r.ledger.collision_tail_tail_energy_change
        << " " << r.ledger.collision_tail_bulk_px_change
        << " " << r.ledger.collision_tail_bulk_energy_change
        << " " << r.ledger.collision_bulk_reaction_px_change
        << " " << r.ledger.collision_bulk_reaction_energy_change
        << " " << r.ledger.collision_reaction_px_residual
        << " " << r.ledger.collision_reaction_energy_residual
        << " " << r.ledger.collision_pair_bulk_bulk
        << " " << r.ledger.collision_pair_tail_tail
        << " " << r.ledger.collision_pair_tail_bulk
        << " " << r.ledger.collision_pair_bulk_reaction
        << " " << r.collision_flux_export_number
        << " " << r.collision_flux_export_energy
        << " " << r.collision_flux_implicit_residual_linf
        << " " << r.collision_flux_cross_pair_residual_linf
        << " " << r.collision_flux_inward_clipped_number
        << " " << r.collision_flux_parcel_count
         << " " << r.collision_flux_rollback_count << " "
         << run_config_.x_transport_velocity_mode << " "
         << energy_vx_min << " " << energy_vx_max << " "
          << energy_vx_linf << " " << energy_vx_l2 << " "
          << run_config_.background_phase_space_mode << " "
          << (r.joint_midpoint_converged ? 1 : 0) << " "
          << r.joint_midpoint_iterations << " "
          << r.joint_midpoint_residual_linf << " "
          << r.joint_midpoint_poisson_residual_linf << " "
          << r.joint_midpoint_energy_residual << "\n";
    if (r.tail_conversion_mode !=
        static_cast<int>(TailConversionMode::STATIC_CELL)) {
        std::ofstream collision_flux(
            (output_dir_ + "/collision_flux_accepted_steps.dat").c_str(),
            std::ios::app);
        collision_flux << std::setprecision(17)
                       << step << " " << time << " "
                       << r.collision_flux_export_number << " "
                       << r.collision_flux_export_energy << " "
                       << r.collision_flux_implicit_residual_linf << " "
                       << r.collision_flux_cross_pair_residual_linf << " "
                       << r.collision_flux_inward_clipped_number << " "
                       << r.collision_flux_parcel_count << " "
                       << r.collision_flux_rollback_count << "\n";
        std::ofstream collision_substeps(
            (output_dir_ + "/collision_flux_substeps_accepted.dat").c_str(),
            std::ios::app);
        for (int half = 0; half < 2; ++half) {
            collision_substeps << std::setprecision(17)
                << step << " " << time << " " << (half + 1) << " "
                << r.collision_flux_half_export_number[half] << " "
                << r.collision_flux_half_export_energy[half] << " "
                << r.collision_flux_half_implicit_residual_linf[half] << " "
                << r.collision_flux_half_cross_pair_residual_linf[half] << " "
                << r.collision_flux_half_inward_clipped_number[half] << " "
                << r.collision_flux_half_parcel_count[half] << " "
                << r.collision_flux_half_rollback_count[half] << "\n";
        }
    }
    last_control_groups_ =
        static_cast<std::uint64_t>(std::max(0, r.population_control_groups));
    last_control_fallbacks_ = static_cast<std::uint64_t>(
        std::max(0, r.population_control_fallbacks));
}

void VpfpDiagnostics::write_stage_energy_audit_accepted_step(
    int step, double time, const VpfpStepResult& result, int rank)
{
    if (!result.accepted || !result.stage_energy_audit_enabled || rank != 0) {
        return;
    }
    std::ofstream out((output_dir_ + "/vpfp_stage_energy_audit.dat").c_str(),
                      std::ios::app);
    if (!out) return;
    out << std::setprecision(17);
    const int count = std::min(result.stage_energy_count,
                               static_cast<int>(VPFP_STAGE_ENERGY_RECORD_COUNT));
    for (int i = 0; i < count; ++i) {
        const VpfpStageEnergyRecord& r = result.stage_energy[i];
        out << step << " " << time << " 1 "
            << (result.stage_energy_audit_valid ? 1 : 0) << " "
            << (result.split_used ? 1 : 0) << " "
            << result.failure_code << " "
            << result.ledger.energy_balance_residual << " "
            << r.stage_id << " " << stage_energy_name(r.stage_id) << " "
            << r.bulk_kinetic << " " << r.tail_kinetic << " "
            << r.beam_kinetic << " " << r.field_energy << " "
            << r.delta_bulk_kinetic << " " << r.delta_tail_kinetic << " "
            << r.delta_beam_kinetic << " " << r.delta_field_energy << " "
            << r.background_left_inflow_energy << " "
            << r.background_left_outflow_energy << " "
            << r.background_right_inflow_energy << " "
            << r.background_right_outflow_energy << " "
            << r.beam_injected_energy << " " << r.beam_outflow_energy << " "
            << r.tail_outflow_energy << " " << r.collision_reservoir_energy << " "
            << r.conversion_energy << " " << r.tail_return_energy << " "
            << r.electrostatic_boundary_work << " " << r.stage_balance << " "
            << r.bulk_upar_face_work << " "
            << r.bulk_upar_velocity_boundary_work << " "
            << r.bulk_upar_interface_energy_removed << " "
            << r.bulk_upar_identity_residual << " "
            << r.tail_kick_work << " " << r.beam_kick_work << "\n";
    }
}

void VpfpDiagnostics::write_field_particle_power_pairing_accepted_step(
    int step, double time, const VpfpStepResult& result, int rank,
    const SpatialGrid& grid,
    const FieldParticlePowerAuditWorkspace* pairing_ws)
{
    // Accepted-only (section 4.6).  The scalar file is rank-0 only; the
    // per-rank profile is written by every rank on interval hits.
    if (!result.accepted || !result.pairing_audit_enabled) return;
    const FieldParticlePowerAuditResult& r = result.pairing_audit;

    if (rank == 0) {
        const std::string path =
            output_dir_ + "/field_particle_power_pairing.dat";
        std::ios::openmode mode =
            pairing_headers_written_ ? std::ios::app : std::ios::trunc;
        std::ofstream out(path.c_str(), mode);
        if (!out) return;
        out << std::setprecision(17);
        if (!pairing_headers_written_) {
            out << "# schema=vpfp-field-particle-power-pairing-v3\n";
            out << "# columns=step time_s dt_s accepted "
                << "continuity_bulk continuity_tail continuity_beam "
                << "poisson_transport_residual "
                << "poisson_endpoint poisson_midpoint poisson_discrete_gradient "
                << "force_work_bulk force_work_tail force_work_beam "
                << "current_pair_residual conversion_residual boundary_residual "
                << "full_residual reconstructed_residual reconstruction_mismatch "
                << "poisson_identity_crosscheck poisson_identity_scale "
                << "poisson_crosscheck_tolerance current_pair_linf "
                << "roundoff_tolerance all_finite continuity_pass "
                << "local_work_ledger_pass reconstruction_pass root_cause_mask "
                << "first_bad_rank first_bad_index "
                << "transport_fraction work_current_fraction time_center_fraction "
                << "pic_fraction conversion_fraction boundary_fraction "
                << "dual_left5_integral dual_core90_integral dual_right5_integral "
                 << "dual_in_domain_work boundary_force_work "
                 << "dual_plus_boundary_work dual_reconstruction_error "
                 << "dual_reconstruction_tolerance dual_reconstruction_pass "
                 << "dual_total_integral\n";
            pairing_headers_written_ = true;
        }
        out << step << " " << time << " " << r.dt_s << " "
            << (result.accepted ? 1 : 0) << " "
            << r.continuity_linf_bulk << " " << r.continuity_linf_tail << " "
            << r.continuity_linf_beam << " "
            << r.poisson_transport_residual << " "
            << r.poisson_transport_residual_endpoint << " "
            << r.poisson_transport_residual_midpoint << " "
            << r.poisson_transport_residual_discrete_gradient << " "
            << r.force_work_residual_bulk << " " << r.force_work_residual_tail << " "
            << r.force_work_residual_beam << " "
            << r.current_pair_residual << " " << r.conversion_transfer_residual << " "
            << r.boundary_residual << " "
            << r.full_residual << " " << r.reconstructed_full_residual << " "
            << r.reconstruction_mismatch << " "
            << r.poisson_identity_crosscheck << " " << r.poisson_identity_scale
            << " " << r.poisson_crosscheck_tolerance << " " << r.current_pair_linf
            << " "
            << r.roundoff_tolerance << " " << (r.all_finite ? 1 : 0) << " "
            << (r.continuity_pass ? 1 : 0) << " "
            << (r.local_work_sum_matches_existing_ledger ? 1 : 0) << " "
            << (r.full_residual_reconstruction_pass ? 1 : 0) << " "
            << r.root_cause_mask << " " << r.first_bad_rank << " "
            << r.first_bad_index << " "
            << r.transport_fraction << " " << r.work_current_fraction << " "
            << r.time_center_fraction << " " << r.pic_fraction << " "
            << r.conversion_fraction << " " << r.boundary_fraction << " "
             << r.dual_left5_integral << " " << r.dual_core90_integral << " "
             << r.dual_right5_integral << " "
                << r.dual_in_domain_work << " " << r.boundary_force_work << " "
             << r.dual_plus_boundary_work << " "
             << r.dual_reconstruction_error << " "
             << r.dual_reconstruction_tolerance << " "
             << (r.dual_reconstruction_pass ? 1 : 0) << " "
             << (r.dual_left5_integral + r.dual_core90_integral +
                 r.dual_right5_integral) << "\n";
    }

    // Per-component first-bad-cell continuity composition (section I2).  Each
    // component owns its own deterministic first-bad location; the owning rank
    // writes that component's breakdown to its per-rank file.
    {
        const std::string bpath = output_dir_ +
            "/field_particle_power_continuity_breakdown_rank" +
            std::to_string(rank) + ".dat";
        std::ios::openmode bmode = pairing_breakdown_headers_written_
            ? std::ios::app : std::ios::trunc;
        std::ofstream bout(bpath.c_str(), bmode);
        if (bout) {
            bout << std::setprecision(17);
            if (!pairing_breakdown_headers_written_) {
                bout << "# schema=vpfp-field-particle-power-continuity-breakdown-v1\n";
                bout << "# columns=step rank component first_bad_global_index "
                     << "delta_N left_face_swept right_face_swept "
                     << "injection_source outflow_source conversion_source "
                     << "return_source other_source residual\n";
                pairing_breakdown_headers_written_ = true;
            }
            const auto emit = [&](const char* name,
                                  const ContinuityComponentBreakdown& bd) {
                if (rank != bd.first_bad_rank || bd.bad_global_index < 0) return;
                bout << step << " " << rank << " " << name << " "
                     << bd.bad_global_index << " " << bd.delta_n << " "
                     << bd.left_face_swept << " " << bd.right_face_swept << " "
                     << bd.injection_source << " " << bd.outflow_source << " "
                     << bd.conversion_source << " " << bd.return_source << " "
                     << bd.other_source << " " << bd.residual << "\n";
            };
            emit("bulk", r.continuity_bulk);
            emit("tail", r.continuity_tail);
            emit("beam", r.continuity_beam);
        }
    }

    // Per-rank profile (section 4.6): level 2 + interval hit only.
    if (level_ < 2 || pairing_ws == NULL) return;
    if (interval_ > 0 && step % interval_ != 0) return;
    const int nxl = static_cast<int>(pairing_ws->bulk_number_n.size());
    if (nxl <= 0) return;
    std::ostringstream pname;
    pname << output_dir_ << "/field_particle_power_pairing_profile_"
          << step << "_rank" << rank << ".dat";
    std::ofstream pout(pname.str().c_str(), std::ios::trunc);
    if (!pout) return;
    pout << std::setprecision(17);
    pout << "# schema=vpfp-field-particle-power-pairing-profile-v2\n";
    pout << "# columns=global_cell x_m E_pair J_charge_bulk J_charge_tail "
         << "J_charge_beam bulk_force_power_density tail_force_power_density "
         << "beam_force_power_density charge_power_density force_power_density "
         << "pairing_residual_density region_id\n";
    const double inv_dt = result.pairing_audit.dt_s > 0.0
        ? 1.0 / result.pairing_audit.dt_s : 0.0;
    const double electron_charge = -Const::qe;
    for (int i = 0; i < nxl; ++i) {
        const double e_pair = 0.5 * (
            pairing_ws->potential_pair_ex_face[static_cast<size_t>(i)] +
            pairing_ws->potential_pair_ex_face[static_cast<size_t>(i + 1)]);
        const auto bulk_face_swept = [&](int f) {
            double swept = 0.0;
            if (f < static_cast<int>(
                    pairing_ws->bulk_x1.bulk_number_swept_face.size())) {
                swept += pairing_ws->bulk_x1.bulk_number_swept_face[
                    static_cast<size_t>(f)];
            }
            if (f < static_cast<int>(
                    pairing_ws->bulk_x2.bulk_number_swept_face.size())) {
                swept += pairing_ws->bulk_x2.bulk_number_swept_face[
                    static_cast<size_t>(f)];
            }
            return swept;
        };
        const auto pic_face_current = [](const std::vector<double>& current,
                                         int f) {
            return f < static_cast<int>(current.size())
                ? current[static_cast<size_t>(f)] : 0.0;
        };
        const double j_bulk = electron_charge * 0.5 *
            (bulk_face_swept(i) + bulk_face_swept(i + 1)) * inv_dt;
        const double j_tail = 0.5 * (
            pic_face_current(
                pairing_ws->tail_trajectory.after_second_drift_current_face, i) +
            pic_face_current(
                pairing_ws->tail_trajectory.after_second_drift_current_face, i + 1));
        const double j_beam = 0.5 * (
            pic_face_current(
                pairing_ws->beam_trajectory.after_second_drift_current_face, i) +
            pic_face_current(
                pairing_ws->beam_trajectory.after_second_drift_current_face, i + 1));
        const double inv_cell_measure = result.pairing_audit.dt_s > 0.0 &&
                                        grid.dx > 0.0
            ? 1.0 / (result.pairing_audit.dt_s * grid.dx) : 0.0;
        const double bulk_power = (i < nxl && i < static_cast<int>(
            pairing_ws->cell_work.bulk_delta_ke_cell.size()))
            ? pairing_ws->cell_work.bulk_delta_ke_cell[static_cast<size_t>(i)] *
                inv_cell_measure : 0.0;
        const double tail_power = (i < nxl && i < static_cast<int>(
            pairing_ws->cell_work.tail_delta_ke_cell.size()))
            ? pairing_ws->cell_work.tail_delta_ke_cell[static_cast<size_t>(i)] *
                inv_cell_measure : 0.0;
        const double beam_power = (i < nxl && i < static_cast<int>(
            pairing_ws->cell_work.beam_delta_ke_cell.size()))
            ? pairing_ws->cell_work.beam_delta_ke_cell[static_cast<size_t>(i)] *
                inv_cell_measure : 0.0;
        const double j_total = j_bulk + j_tail + j_beam;
        const double charge_power = e_pair * j_total;
        const double force_power = bulk_power + tail_power + beam_power;
        const int global_cell = grid.ix_start + i;
        int region = 1;
        if (global_cell < static_cast<int>(std::ceil(0.05 * grid.nx_global)))
            region = 0;
        else if (global_cell >=
                  static_cast<int>(std::floor(0.95 * grid.nx_global)))
            region = 2;
        pout << global_cell << " " << (global_cell + 0.5) * grid.dx << " "
             << e_pair << " " << j_bulk << " " << j_tail << " " << j_beam
             << " " << bulk_power << " " << tail_power << " " << beam_power
             << " " << charge_power << " " << force_power << " "
             << (force_power - charge_power) << " " << region << "\n";
    }

    // §15.13.4 step 3: per-rank direct face dual audit.  The face profile
    // carries the charge-conserving face current, the reconstructed force
    // cell current, its G* face alias and the per-face residual.  This is
    // read-only; nothing here feeds back into the production state.
    {
        const std::vector<double>& j_charge =
            result.pairing_audit.j_charge_face;
        const std::vector<double>& j_force =
            result.pairing_audit.j_force_cell;
        const std::vector<double>& gstar =
            result.pairing_audit.gstar_j_force_face;
        const std::vector<double>& dual =
            result.pairing_audit.dual_face_residual;
        const std::vector<double>& e_pair_face =
            pairing_ws->potential_pair_ex_face;
        const size_t nface = gstar.size();
        if (nface > 0 && j_charge.size() == nface &&
            dual.size() == nface && e_pair_face.size() == nface) {
            std::ostringstream fname;
            fname << output_dir_
                  << "/field_particle_power_dual_face_"
                  << step << "_rank" << rank << ".dat";
            std::ofstream fout(fname.str().c_str(), std::ios::trunc);
            if (fout) {
                fout << std::setprecision(17);
                fout << "# schema=vpfp-field-particle-power-dual-face-v2\n";
                fout << "# columns=global_face x_m E_pair J_charge_face "
                     << "J_force_cell_left J_force_cell_right "
                     << "Gstar_J_force_face J_charge_face_minus_Gstar_J_force_face "
                     << "quadrature_weight face_owner_rank face_is_owner "
                     << "physical_endpoint region_id\n";
                const size_t ncell = j_force.size();
                for (size_t f = 0; f < nface; ++f) {
                    const int global_face = grid.ix_start + static_cast<int>(f);
                    const double j_force_left = f > 0 && f - 1 < ncell
                        ? j_force[f - 1] : 0.0;
                    const double j_force_right = f < ncell
                        ? j_force[f] : 0.0;
                    int region = 1;
                    if (global_face <
                        static_cast<int>(std::ceil(0.05 * grid.nx_global)))
                        region = 0;
                    else if (global_face >=
                             static_cast<int>(std::floor(0.95 * grid.nx_global)))
                        region = 2;
                    const double quadrature_weight =
                        (global_face == 0 || global_face == grid.nx_global)
                        ? 0.5 * grid.dx
                        : ((f + 1 == nface &&
                            global_face < grid.nx_global) ? 0.0 : grid.dx);
                    const int owner_rank =
                        (f + 1 == nface && global_face < grid.nx_global)
                        ? rank + 1 : rank;
                    const int face_is_owner = owner_rank == rank ? 1 : 0;
                    const int physical_endpoint =
                        (global_face == 0 || global_face == grid.nx_global)
                        ? 1 : 0;
                    fout << global_face << " " << global_face * grid.dx << " "
                         << e_pair_face[f] << " " << j_charge[f] << " "
                         << j_force_left << " " << j_force_right << " "
                         << gstar[f] << " " << dual[f] << " "
                         << quadrature_weight << " " << owner_rank << " "
                         << face_is_owner << " " << physical_endpoint << " "
                         << region
                         << "\n";
                }
            }
        }
    }
}

// JC4 (section 7.3): field-particle iteration accepted-step diagnostic.
// Written by rank 0 only; header validated on restart-append.  The file
// records one row per accepted step with the final Picard iteration state.
void VpfpDiagnostics::write_field_particle_iteration_accepted_step(
    int step, double time, const VpfpStepResult& result, int rank)
{
    if (!result.accepted) return;
    if (rank != 0) return;
    // §7.3 rule: diagnostic-level=0 writes only summary counts (not per-
    // iteration), level=1 writes the final accepted-state line, level=2
    // writes the same accepted line.  The iteration-level detail is not
    // emitted here; it is in the stderr [jc3-picard] logs.  This file
    // always writes one row per accepted step regardless of level.
    const std::string path =
        output_dir_ + "/field_particle_iteration.dat";
    std::ios::openmode mode =
        field_particle_iteration_headers_written_ ? std::ios::app
                                                  : std::ios::trunc;
    std::ofstream out(path.c_str(), mode);
    if (!out) return;
    out << std::setprecision(17);
    if (!field_particle_iteration_headers_written_) {
        out << "# schema=vpfp-field-particle-iteration-v1\n"
            << "# columns=step time_s mode iterations converged relaxation "
               "field_residual_l2 field_residual_linf pairing_residual "
               "trial_evaluations post_field_charge_residual failure_code\n";
        field_particle_iteration_headers_written_ = true;
    }
    const char* mode_name = result.field_particle_coupling_enabled
                                ? "discrete-gradient" : "legacy";
    out << step << " " << time << " " << mode_name << " "
        << result.field_particle_iterations << " "
        << (result.field_particle_converged ? 1 : 0) << " "
        << result.field_particle_relaxation << " "
        << result.field_particle_residual_l2 << " "
        << result.field_particle_residual_linf << " "
        << result.field_particle_pairing_residual << " "
        << result.field_particle_trial_evaluations << " "
        << result.post_field_charge_residual_linf << " "
        << result.failure_code << "\n";
}

void VpfpDiagnostics::write_bulk_tail_flux_accepted_step(
    int step, double time, const VpfpStepResult& result, int rank)
{
    if (!result.accepted || result.tail_conversion_mode ==
        static_cast<int>(TailConversionMode::STATIC_CELL)) return;
    if (rank != 0) return;
    if (result.flux_face_audit_max_valid &&
        result.flux_face_audit_max_relative > 0.0) {
        const BulkTailFluxFaceAudit& a = result.flux_face_audit_max;
        std::ofstream face_out(
            (output_dir_ + "/bulk_tail_flux_face_accepted_steps.dat").c_str(),
            std::ios::app);
        if (face_out) {
            face_out << std::setprecision(17)
                     << step << " " << time << " "
                     << a.ix_global << " " << a.direction << " "
                     << a.face_index << " " << a.transverse_index << " "
                     << a.operator_stage << " " << a.face_number << " "
                     << a.parcel_number << " "
                     << result.flux_face_audit_abs_at_max_relative << " "
                     << result.flux_face_audit_max_relative << " "
                     << a.node_failure_reason << " "
                     << a.reconstructed_target << " " << a.node_sum << "\n";
        }
    }
    std::ofstream out((output_dir_ + "/bulk_tail_flux_accepted_steps.dat").c_str(),
                      std::ios::app);
    out << std::setprecision(17)
        << step << " " << time << " "
        << tail_conversion_mode_name(static_cast<TailConversionMode>(
               result.tail_conversion_mode)) << " "
        << result.flux_quadrature_order << " " << result.flux_max_supports << " "
        << result.flux_parcel_count << " " << result.flux_node_count << " "
        << result.flux_export_number << " " << result.flux_export_energy << " "
        << result.flux_face_export_number << " "
        << result.flux_parcel_number << " "
        << result.flux_below_threshold_number << " "
        << result.flux_roundoff_discarded_number << " "
        << result.flux_duplicate_count << " "
        << result.flux_quadrature_error_max << " "
        << result.flux_tail_owned_expected_transfer_number << " "
        << result.flux_tail_owned_roundoff_discarded_number << " "
        << result.flux_tail_owned_bulk_residual << " "
        << result.static_extractor_call_count << " "
        << result.ledger.conversion_particles_created << " "
        << result.flux_compression_fallback_count << " "
        << result.flux_subcell_fallback_count << " "
        << result.flux_support_limit_violation_count << " "
        << result.flux_duplicate_id_count << " "
        << result.flux_face_ledger_mismatch_count << " "
        << result.flux_face_audit_count << " "
        << result.flux_face_audit_face_abs_sum << " "
        << result.flux_face_audit_parcel_abs_sum << " "
        << result.flux_face_audit_abs_error_sum << " "
        << result.flux_face_audit_max_relative << " "
        << result.flux_face_audit_abs_at_max_relative << " "
        << result.flux_conversion_wall_seconds << " "
        << result.ledger.conversion_number_removed << " "
        << result.ledger.conversion_energy_removed << " "
        << result.conversion_number_residual << " "
        << result.conversion_px_residual << " "
        << result.conversion_energy_residual << " "
        << result.conversion_jx_residual << " "
        << result.conversion_pixx_residual << " "
        << result.conversion_piperp_residual << " "
        << result.conversion_rho_l2 << " "
        << result.audit_parcel_failure_reason << " "
        << result.audit_parcel_failure_rank << " "
        << result.audit_parcel_failure_ix << " "
        << result.audit_parcel_failure_face << " "
        << result.audit_parcel_failure_iuperp << " "
        << result.audit_parcel_failure_node_mass << " "
        << result.audit_parcel_failure_target << " "
        << result.audit_parcel_failure_node_sum << " "
        << result.audit_parcel_failure_scale << " "
        << (result.audit_valid ? 1 : 0) << " "
        << result.audit_failure_code << " "
        << (result.audit_inplace_state_bitwise_equal ? 1 : 0) << " "
        << (result.audit_inplace_rng_equal ? 1 : 0) << " "
        << (result.audit_inplace_ledger_equal ? 1 : 0) << "\n";
}

void VpfpDiagnostics::write_bulk_tail_moment_audit_accepted_step(
    int step, double time, const VpfpStepResult& result, int rank)
{
    if (!result.accepted || !run_config_.tail_cell_moment_audit) return;
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // Section 7.11.16B: every rank checks the event layout before entering
    // any per-event collective; a mismatch follows one failure path.
    if (bulk_tail_moment_audit_check_event_layout(
            result.conversion_events.size(), rank, mpi_size) !=
        BulkTailMomentAuditIoStatus::OK) {
        if (rank == 0) {
            std::ofstream out((output_dir_ +
                "/bulk_tail_moment_audit_accepted_steps.dat").c_str(),
                std::ios::app);
            if (out) out << "INVALID_EVENT_LAYOUT step=" << step << "\n";
        }
        return;
    }
    const size_t top_cell_limit = static_cast<size_t>(
        run_config_.tail_cell_moment_audit_top_cells);
    for (size_t event = 0; event < result.conversion_events.size(); ++event) {
        const MomentRepresentationAudit& a =
            result.conversion_events[event].moment_audit;
        if (!a.enabled) continue;
        BulkTailMomentAuditGlobal g;
        const BulkTailMomentAuditIoStatus status =
            bulk_tail_moment_audit_reduce_event(
                a, static_cast<int>(result.conversion_events[event].location),
                top_cell_limit, rank, mpi_size, g);
        if (status != BulkTailMomentAuditIoStatus::OK) {
            if (rank == 0) {
                std::ofstream out((output_dir_ +
                    "/bulk_tail_moment_audit_accepted_steps.dat").c_str(),
                    std::ios::app);
                if (out) {
                    out << "INVALID_EVENT_LAYOUT step=" << step
                        << " event=" << event << "\n";
                }
            }
            return;   // every rank takes the same failure path
        }
        if (!g.has_positive_requests) continue;   // no row for empty events
        if (rank != 0) continue;

        const int location =
            static_cast<int>(result.conversion_events[event].location);
        // Main accepted-step table (section 7.11.16B item 6: raw
        // center_l1/delta_l1 are written explicitly; relative columns use
        // the r rule with relative_defined flags).
        std::ofstream out((output_dir_ +
            "/bulk_tail_moment_audit_accepted_steps.dat").c_str(),
            std::ios::app);
        if (!moment_audit_headers_written_) {
            out << "accepted_step time_fs conversion_location "
                << "request_cell_count positive_request_cell_count "
                << "below_threshold_number_fraction "
                << "threshold_window_6p0_6p2_N "
                << "volume_target_feasible_count volume_target_failed_count "
                << "eligible_target_feasible_count "
                << "eligible_target_failed_count";
            for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
                out << " center_" << m << " volume_" << m
                    << " eligible_raw_" << m << " eligible_normalized_" << m
                    << " delta_signed_" << m << " center_l1_" << m
                    << " delta_l1_" << m << " rel_l1_" << m
                    << " rel_signed_" << m << " max_cell_rel_" << m
                    << " relative_defined_" << m;
            out << "\n";
            moment_audit_headers_written_ = true;
        }
        const double center_n = g.center[BULK_TAIL_MOMENT_N];
        const double below_threshold =
            center_n > 0.0
                ? std::max(0.0, 1.0 - g.eligible_raw[BULK_TAIL_MOMENT_N] /
                                        center_n)
                : 0.0;
        out << std::setprecision(17) << step << " " << time / Const::femto
            << " " << location
            << " " << g.request_cell_count << " "
            << g.positive_request_cell_count << " " << below_threshold
            << " " << g.threshold_window_number << " "
            << g.volume_target_feasible_count << " "
            << g.volume_target_failed_count << " "
            << g.eligible_target_feasible_count << " "
            << g.eligible_target_failed_count;
        for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
            const double rel_l1 =
                g.center_l1[m] > 0.0
                    ? g.delta_l1[m] / g.center_l1[m]
                    : (g.delta_l1[m] == 0.0
                           ? 0.0
                           : std::numeric_limits<double>::infinity());
            const double rel_signed =
                g.center_l1[m] > 0.0
                    ? std::fabs(g.delta_signed[m]) / g.center_l1[m]
                    : (g.delta_signed[m] == 0.0
                           ? 0.0
                           : std::numeric_limits<double>::infinity());
            out << " " << g.center[m] << " " << g.volume[m]
                << " " << g.eligible_raw[m] << " "
                << g.eligible_normalized[m] << " " << g.delta_signed[m]
                << " " << g.center_l1[m] << " " << g.delta_l1[m]
                << " " << rel_l1 << " " << rel_signed
                << " " << g.max_cell_relative[m]
                << " " << static_cast<int>(g.relative_defined[m]);
        }
        out << "\n";

        // Velocity histogram (section 7.11.16B item 7): only real accepted
        // events, one row per (iv, imu) with global count > 0, sorted by
        // (iv, imu).
        std::ofstream hist_out((output_dir_ +
            "/bulk_tail_moment_audit_velocity_histogram.dat").c_str(),
            std::ios::app);
        if (!moment_audit_hist_headers_written_) {
            hist_out << "accepted_step time_fs conversion_location iv imu "
                        "request_cell_count request_number\n";
            moment_audit_hist_headers_written_ = true;
        }
        for (size_t b = 0; b < g.velocity_bins.size(); ++b) {
            const BulkTailVelocityBinAudit& bin = g.velocity_bins[b];
            hist_out << std::setprecision(17) << step << " "
                     << time / Const::femto << " " << location << " "
                     << bin.iv << " " << bin.imu << " "
                     << bin.request_cell_count << " " << bin.request_number
                     << "\n";
        }

        // Global top cells (second global truncation already applied inside
        // the reduction).
        std::ofstream top_out((output_dir_ +
            "/bulk_tail_moment_audit_top_cells.dat").c_str(), std::ios::app);
        for (size_t q = 0; q < g.top_cells.size(); ++q) {
            const BulkTailMomentAuditTopCell& cell = g.top_cells[q];
            top_out << std::setprecision(17) << step << " "
                    << time / Const::femto << " " << cell.rank << " "
                    << cell.ix_global << " " << cell.iv << " " << cell.imu
                    << " " << cell.score;
            for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m)
                top_out << " " << cell.center[m] << " " << cell.volume[m];
            top_out << "\n";
        }
    }
}

void VpfpDiagnostics::write_conversion_source_ledger(
    int step, double time, const VpfpStepResult& result,
    const BackgroundTailPIC* tail, int rank)
{
    if (level_ < 2 || result.conversion_events.empty()) return;
    for (size_t event = 0; event < result.conversion_events.size(); ++event) {
        const BulkTailConversionDiagnostics& d = result.conversion_events[event];
        const size_t nb = d.removed_bulk_number_spectrum.size();
        if (d.conversion_source_energy_edges.size() != nb + 1) continue;
        std::vector<double> local(8 * nb, 0.0);
        for (size_t b = 0; b < nb; ++b) {
            local[b] = d.pre_extraction_bulk_number_spectrum[b];
            local[nb + b] = d.pre_extraction_bulk_energy_spectrum[b];
            local[2 * nb + b] = d.removed_bulk_number_spectrum[b];
            local[3 * nb + b] = d.removed_bulk_energy_spectrum[b];
            local[4 * nb + b] = d.created_tail_number_spectrum[b];
            local[5 * nb + b] = d.created_tail_energy_spectrum[b];
        }
        // Re-bin the actual final accepted tail state.  These columns are
        // total accepted-tail spectra, not aliases of this event's created
        // source; pre-existing particles are intentionally included.
        if (tail != NULL) {
            for (size_t p = 0; p < tail->particles.size(); ++p) {
                const BackgroundTailParticle& q = tail->particles[p];
                const double gamma = std::sqrt(
                    1.0 + q.ux * q.ux + q.uy * q.uy + q.uz * q.uz);
                const double ke = Const::me * Const::c * Const::c *
                                  (gamma - 1.0);
                const std::vector<double>::const_iterator upper =
                    std::upper_bound(d.conversion_source_energy_edges.begin(),
                                     d.conversion_source_energy_edges.end(), ke);
                if (upper == d.conversion_source_energy_edges.begin()) continue;
                const size_t b = static_cast<size_t>(
                    upper - d.conversion_source_energy_edges.begin() - 1);
                if (b >= nb) continue;
                local[6 * nb + b] += q.weight;
                local[7 * nb + b] += q.weight * ke;
            }
        }
        std::vector<double> global(local.size(), 0.0);
        MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Candidate subcell loads must leave an explicit local fallback
        // trail.  Each MPI rank writes only its own records, avoiding an
        // unsafe concurrent append while retaining the physical cell index.
        if (!d.subcell_fallbacks.empty()) {
            std::ostringstream path;
            path << output_dir_ << "/conversion_subcell_fallbacks_rank"
                 << rank << ".dat";
            std::ofstream fallback_out(path.str().c_str(), std::ios::app);
            for (size_t q = 0; q < d.subcell_fallbacks.size(); ++q) {
                const BulkTailConversionDiagnostics::SubcellFallbackRecord& r =
                    d.subcell_fallbacks[q];
                fallback_out << std::setprecision(17)
                             << step << " " << time << " "
                             << static_cast<int>(d.location) << " "
                             << r.ix_global << " " << r.iv << " " << r.imu
                             << " " << r.reason << " "
                             << r.fallback_particles << " "
                             << r.number_target << " " << r.px_target << " "
                             << r.jx_target << " " << r.energy_target << " "
                             << r.pixx_target << " " << r.piperp_target
                             << "\n";
            }
        }
        if (rank != 0) continue;
        std::ofstream out((output_dir_ + "/conversion_source_accepted_steps.dat").c_str(),
                          std::ios::app);
        for (size_t b = 0; b < nb; ++b) {
            out << std::setprecision(17) << step << " " << time << " "
                << static_cast<int>(d.location) << " "
                << d.conversion_source_energy_edges[b] << " "
                << d.conversion_source_energy_edges[b + 1] << " "
                << global[b] << " " << global[nb + b] << " "
                << global[2 * nb + b] << " " << global[3 * nb + b] << " "
                << global[4 * nb + b] << " " << global[5 * nb + b] << " "
                << global[6 * nb + b] << " " << global[7 * nb + b] << "\n";
        }
        last_conversion_source_edges_ = d.conversion_source_energy_edges;
        last_conversion_source_spectrum_.assign(global.begin() + 2 * nb,
                                                 global.begin() + 3 * nb);

    }
}

void VpfpDiagnostics::write_population_control(
    int step, double time, const VpfpStepResult& r, int rank)
{
    if (rank != 0 || !r.population_control_applied) return;
    const std::string path = output_dir_ + "/tail_population_control.dat";
    std::ofstream out(path.c_str(), std::ios::app);
    if (!out) return;
    out << std::setprecision(17) << step << " " << time << " "
        << "groups=" << r.population_control_groups
        << " fallbacks=" << r.population_control_fallbacks
        << " particles_before=" << r.population_control_particles_before
        << " particles_after=" << r.population_control_particles_after
        << " residual_n=" << r.population_control_max_residual[0]
        << " residual_px=" << r.population_control_max_residual[1]
        << " residual_jx=" << r.population_control_max_residual[2]
        << " residual_ke=" << r.population_control_max_residual[3]
        << " residual_pixx=" << r.population_control_max_residual[4]
        << " residual_piperp=" << r.population_control_max_residual[5]
        << " residual_xw=" << r.population_control_max_residual[6] << "\n";
}

void VpfpDiagnostics::write_failure(int step, double time,
                                    const VpfpStepResult& r, int rank)
{
    if (rank != 0) return;
    write_joint_midpoint_iterations(step, time, r, rank);
    const char* audit_reason = "none";
    if (r.audit_failure_code & 1) audit_reason = "parcel_nodes";
    else if (r.audit_failure_code & 2) audit_reason = "parcel_nonfinite";
    else if (r.audit_failure_code & 4) audit_reason = "parcel_negative";
    else if (r.audit_failure_code & 8) audit_reason = "duplicate_parcel";
    else if (r.audit_failure_code & 16) audit_reason = "quadrature";
    else if (r.audit_failure_code & 32) audit_reason = "observer_flux_mutation";
    else if (r.audit_failure_code & 64) audit_reason = "physical_state_mutation";
    else if (r.audit_failure_code & 128) audit_reason = "rng_mutation";
    else if (r.audit_failure_code & 256) audit_reason = "ledger_mutation";
    std::ofstream out((output_dir_ + "/vpfp_failure.dat").c_str(), std::ios::app);
    out << std::setprecision(17) << step << " " << time << " failure_code="
        << r.failure_code << " finite=" << r.finite << " cfl_ok=" << r.cfl_ok
        << " gauss_ok=" << r.gauss_ok << " collision_ok=" << r.collision_ok
        << " failure_stage=" << r.failure_stage
        << " failing_rank=" << r.failing_rank
        << " failing_ix=" << r.failing_ix
        << " failing_iupar=" << r.failing_iupar
        << " failing_iuperp=" << r.failing_iuperp
        << " input_min=" << r.input_min
        << " input_max=" << r.input_max
        << " output_min=" << r.output_min
        << " output_max=" << r.output_max
        << " first_nonfinite_value=" << r.first_nonfinite_value
        << " audit_valid=" << r.audit_valid
        << " audit_failure_code=" << r.audit_failure_code
        << " audit_failure_reason=" << audit_reason
        // JC4 (section 7.7): field-particle coupling diagnostics on failure.
        << " field_particle_coupling_enabled=" << (r.field_particle_coupling_enabled ? 1 : 0)
        << " field_particle_converged=" << (r.field_particle_converged ? 1 : 0)
        << " field_particle_iterations=" << r.field_particle_iterations
        << " field_particle_trial_evaluations=" << r.field_particle_trial_evaluations
        << " field_particle_relaxation=" << r.field_particle_relaxation
        << " field_particle_residual_l2=" << r.field_particle_residual_l2
        << " field_particle_residual_linf=" << r.field_particle_residual_linf
        << " field_particle_pairing_residual=" << r.field_particle_pairing_residual
        << " post_field_charge_residual_linf=" << r.post_field_charge_residual_linf
        << "\n";
    if (r.tail_conversion_mode !=
        static_cast<int>(TailConversionMode::STATIC_CELL)) {
        std::ofstream flux((output_dir_ +
            "/trial_bulk_tail_flux_failures.dat").c_str(), std::ios::app);
        flux << std::setprecision(17)
             << "step=" << step << " time_s=" << time
             << " failure_code=" << r.failure_code
             << " conversion_mode="
             << tail_conversion_mode_name(static_cast<TailConversionMode>(
                    r.tail_conversion_mode))
             << " parcel_count=" << r.flux_parcel_count
             << " node_count=" << r.flux_node_count
             << " exported_number=" << r.flux_export_number
             << " exported_energy=" << r.flux_export_energy
             << " face_export_number=" << r.flux_face_export_number
             << " parcel_number=" << r.flux_parcel_number
             << " duplicate_count=" << r.flux_duplicate_count
             << " quadrature_error_max=" << r.flux_quadrature_error_max
             << " compression_fallback_count="
             << r.flux_compression_fallback_count
             << " subcell_fallback_count="
             << r.flux_subcell_fallback_count
             << " support_limit_violation_count="
             << r.flux_support_limit_violation_count
             << " duplicate_id_count=" << r.flux_duplicate_id_count
             << " face_ledger_mismatch_count="
             << r.flux_face_ledger_mismatch_count
             << " conversion_wall_seconds="
             << r.flux_conversion_wall_seconds
             << " tail_owned_expected_transfer_number="
             << r.flux_tail_owned_expected_transfer_number
             << " tail_owned_roundoff_discarded_number="
             << r.flux_tail_owned_roundoff_discarded_number
             << " tail_owned_bulk_residual="
             << r.flux_tail_owned_bulk_residual << "\n";
    }
    if (!r.conversion_events.empty()) {
        const BulkTailConversionDiagnostics& c = r.conversion_events.back();
        out << "conversion_complete=" << c.complete
            << " conversion_conservative=" << c.conservative
            << " conversion_fidelity_ok=" << c.fidelity_ok
            << " conversion_finite=" << c.finite
            << " conversion_particles_created=" << c.particles_created
            << " conversion_number_residual=" << c.number_residual_rel
            << " conversion_px_residual=" << c.px_residual_rel
            << " conversion_energy_residual=" << c.energy_residual_rel
            << " conversion_jx_residual=" << c.jx_residual_rel
            << " conversion_pixx_residual=" << c.pixx_residual_rel
            << " conversion_piperp_residual=" << c.piperp_residual_rel
            << " conversion_support_limit_violations="
            << c.support_limit_violation_count
            << " conversion_duplicate_ids=" << c.duplicate_id_count
            << " conversion_face_ledger_mismatches="
            << c.face_ledger_mismatch_count << "\n";
    }
}

void VpfpDiagnostics::write_joint_midpoint_iterations(
    int step, double time, const VpfpStepResult& result, int rank)
{
    if (rank != 0 || !result.joint_midpoint_enabled ||
        result.joint_midpoint_iterations_log.empty()) return;
    const std::string path = output_dir_ + "/joint_midpoint_iteration.dat";
    bool empty = true;
    {
        std::ifstream probe(path.c_str(), std::ios::in | std::ios::binary);
        if (probe.good()) empty = probe.peek() == std::ifstream::traits_type::eof();
    }
    std::ofstream out(path.c_str(), std::ios::app);
    if (!out) return;
    if (empty) {
        out << "step time_s iteration gmres_dimension residual_linf "
               "phi_residual_linf line_search_alpha trial_min_mass accepted "
               "failure_code\n";
    }
    out << std::setprecision(17);
    for (size_t i = 0; i < result.joint_midpoint_iterations_log.size(); ++i) {
        const JointPhaseSpaceIterationRecord& record =
            result.joint_midpoint_iterations_log[i];
        out << step << " " << time << " " << record.iteration << " "
            << record.gmres_dimension << " " << record.residual_linf << " "
            << record.phi_residual_linf << " " << record.line_search_alpha << " "
            << record.trial_min_mass << " " << record.accepted << " "
            << record.failure_code << "\n";
    }
}

namespace {

std::string snapshot_directory(const std::string& output_dir,
                               int step, double time)
{
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer), "%s/snapshot_t%.9gfs_step%d",
                  output_dir.c_str(), time / Const::femto, step);
    return std::string(buffer);
}

bool make_directory_if_missing(const std::string& path)
{
#ifdef _WIN32
    if (_mkdir(path.c_str()) == 0) return true;
#else
    if (mkdir(path.c_str(), 0755) == 0) return true;
#endif
    // EEXIST (or any other failure) is treated as usable only if the caller
    // can still open files inside; opening below reports the real error.
    return true;
}

} // namespace

void VpfpDiagnostics::write_snapshot(int step, double time,
                                     const Species& electrons,
                                     const BeamPIC& beam,
                                     const EMFields& fields,
                                     const SpatialGrid& grid,
                                     const BackgroundTailPIC* tail,
                                     const HybridVelocityPartition* partition,
                                     double convert_energy_mev,
                                     int mpi_rank, int mpi_size)
{
    if (level_ < 1) return;
    const std::string directory =
        snapshot_directory(output_dir_, step, time);
    make_directory_if_missing(directory);
    const std::string complete_marker = directory + "/_COMPLETE";
    const std::string pending_marker = complete_marker + ".tmp";
    if (mpi_rank == 0) {
        std::remove(complete_marker.c_str());
        std::remove(pending_marker.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    write_fields(step, time, electrons, beam, fields, grid, mpi_rank);
    write_moments_and_spectra(step, time, electrons, beam, grid, tail,
                              mpi_rank);
    write_tail_stats(step, time, electrons, grid, tail, mpi_rank);
    write_threshold_interface(step, time, electrons, grid, tail, partition,
                              convert_energy_mev, mpi_rank);
    write_manifest(step, time, electrons, beam, grid, tail, partition,
                   convert_energy_mev, mpi_rank, mpi_size);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::ofstream complete(pending_marker.c_str());
        complete << "schema=vpfp-snapshot-complete-v1\n"
                 << "step=" << step << "\n"
                 << "time_s=" << std::setprecision(17) << time << "\n"
                 << "ranks=" << mpi_size << "\n";
        complete.close();
        if (complete) std::rename(pending_marker.c_str(),
                                  complete_marker.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void VpfpDiagnostics::write_fields(int step, double time,
                                   const Species& electrons,
                                   const BeamPIC& beam,
                                   const EMFields& fields,
                                   const SpatialGrid& grid, int mpi_rank)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer),
                  "%s/snapshot_t%.9gfs_step%d", output_dir_.c_str(),
                  time / Const::femto, step);
    const std::string directory(buffer);

    {
        std::ofstream out((directory + "/fields_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step
            << " ix_start=" << grid.ix_start << " nx_local=" << nxl << "\n";
        out << "x_um Ex_Vm phi_V rho_Cm3\n";
        for (int ix = 0; ix < nxl; ++ix) {
            const int ig = ix + ng;
            out << std::setprecision(17)
                << grid.x(ig) / Const::micro << " "
                << fields.Ex[ig] << " "
                << fields.phi[ig] << " "
                << fields.rho[ig] << "\n";
        }
    }

    {
        std::ofstream out((directory + "/density_background_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step
            << " ix_start=" << grid.ix_start << " nx_local=" << nxl << "\n";
        out << "x_um n_e_m3 rho_e_Cm3\n";
        for (int ix = 0; ix < nxl; ++ix) {
            out << std::setprecision(17)
                << grid.x(ix + ng) / Const::micro << " "
                << electrons.number_density[ix] << " "
                << electrons.charge_density[ix] << "\n";
        }
    }

    {
        std::ofstream out((directory + "/density_beam_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step
            << " ix_start=" << grid.ix_start << " nx_local=" << nxl << "\n";
        out << "x_um n_b_m3\n";
        const size_t beam_size =
            beam.density.size() >= static_cast<size_t>(grid.nx_local)
                ? static_cast<size_t>(grid.nx_local) : beam.density.size();
        for (int ix = 0; ix < nxl; ++ix) {
            const double nb =
                ix < static_cast<int>(beam_size) ? beam.density[ix] : 0.0;
            out << std::setprecision(17)
                << grid.x(ix + ng) / Const::micro << " " << nb << "\n";
        }
    }
}

void VpfpDiagnostics::write_moments_and_spectra(
    int step, double time, const Species& electrons, const BeamPIC& beam,
    const SpatialGrid& grid, const BackgroundTailPIC* tail, int mpi_rank)
{
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int nv = static_cast<int>(electrons.cgrid.upar_cells.size());
    const int nmu = static_cast<int>(electrons.cgrid.uperp_cells.size());
    if (nv <= 0 || nmu <= 0) return;
    const size_t nvmu = static_cast<size_t>(nv) * nmu;

    char buffer[256];
    std::snprintf(buffer, sizeof(buffer),
                  "%s/snapshot_t%.9gfs_step%d", output_dir_.c_str(),
                  time / Const::femto, step);
    const std::string directory(buffer);

    // Spatial moments (number, current, kinetic-energy density) plus the
    // background spectra.  Only the final accepted state is read; the scan
    // happens at snapshot cadence, never per step (section 13.13).
    const int nbin = 256;
    double emax = 0.0;
    for (size_t s = 0; s < nvmu; ++s)
        emax = std::max(emax, electrons.cgrid.kinetic_energy[s]);
    if (emax <= 0.0) emax = 1.0;
    // Log-spaced energy bins from well below the 100 eV thermal peak to the
    // grid maximum (10.5 MeV at u_max=22.4), so both the thermal core and
    // the instability tail stay resolved in the same spectrum.
    const double emin = 0.01 * Const::eV;
    const double log_range = std::log(emax / emin);
    std::vector<double> count_bg(static_cast<size_t>(nbin), 0.0);
    std::vector<double> count_tail(static_cast<size_t>(nbin), 0.0);
    std::vector<double> count_beam(static_cast<size_t>(nbin), 0.0);
    std::vector<double> upar_bg(static_cast<size_t>(nv), 0.0);
    std::vector<double> upar_tail(static_cast<size_t>(nv), 0.0);
    std::vector<double> uperp_bg(static_cast<size_t>(nmu), 0.0);

    {
        std::ofstream out((directory + "/distribution_moments_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step
            << " ix_start=" << grid.ix_start << " nx_local=" << nxl << "\n";
        out << "x_um n_e_m3 current_x_A_m2 ke_density_J_m3\n";
        for (int ix = 0; ix < nxl; ++ix) {
            const int ix_g = ix + ng;
            const size_t xbase = static_cast<size_t>(ix_g) * nvmu;
            double ke = 0.0;
            for (int j = 0; j < nv; ++j) {
                const size_t row = xbase + static_cast<size_t>(j) * nmu;
                for (int k = 0; k < nmu; ++k) {
                    const double m = electrons.f[row + k];
                    ke += m * electrons.cgrid.kinetic_energy[
                        static_cast<size_t>(j) * nmu + k];
                }
            }
            out << std::setprecision(17)
                << grid.x(ix_g) / Const::micro << " "
                << electrons.number_density[ix] << " "
                << electrons.current_x[ix] << " "
                << ke / grid.dx << "\n";
        }
    }

    for (int ix = 0; ix < nxl; ++ix) {
        const size_t xbase = static_cast<size_t>(ix + ng) * nvmu;
        for (int j = 0; j < nv; ++j) {
            const size_t row = xbase + static_cast<size_t>(j) * nmu;
            for (int k = 0; k < nmu; ++k) {
                const double m = electrons.f[row + k];
                const double e = electrons.cgrid.kinetic_energy[
                    static_cast<size_t>(j) * nmu + k];
                int bin = (e > emin)
                    ? static_cast<int>(std::log(e / emin) / log_range * nbin)
                    : 0;
                if (bin >= nbin) bin = nbin - 1;
                count_bg[static_cast<size_t>(bin)] += m;
                upar_bg[static_cast<size_t>(j)] += m;
                uperp_bg[static_cast<size_t>(k)] += m;
            }
        }
    }

    // Beam macro-particle spectra (energy and u_parallel projections).
    for (size_t p = 0; p < beam.particles.size(); ++p) {
        const double u = beam.particles[p].px /
                         (Const::me * Const::c);
        const double gamma = std::sqrt(1.0 + u * u);
        const double e = Const::me * Const::c * Const::c * (gamma - 1.0);
        int bin = (e > emin)
            ? static_cast<int>(std::log(e / emin) / log_range * nbin)
            : 0;
        if (bin >= nbin) bin = nbin - 1;
        count_beam[static_cast<size_t>(bin)] += beam.particles[p].weight;
    }

    // Tail macro-particle spectra on the same unified energy bins
    // (section 13.2: combined spectrum must use one binning), plus the
    // u_parallel projection.
    if (tail != NULL) {
        for (size_t p = 0; p < tail->particles.size(); ++p) {
            const BackgroundTailParticle& tp = tail->particles[p];
            const double gamma = std::sqrt(
                1.0 + tp.ux * tp.ux + tp.uy * tp.uy + tp.uz * tp.uz);
            const double e =
                Const::me * Const::c * Const::c * (gamma - 1.0);
            int bin = (e > emin)
                ? static_cast<int>(std::log(e / emin) / log_range * nbin)
                : 0;
            if (bin >= nbin) bin = nbin - 1;
            count_tail[static_cast<size_t>(bin)] += tp.weight;
            int j = 0;
            for (int jj = 1; jj < nv; ++jj) {
                if (tp.ux > electrons.cgrid.upar_cells[
                        static_cast<size_t>(jj)]) {
                    j = jj;
                }
            }
            if (tp.ux >= electrons.cgrid.upar_faces.front() &&
                tp.ux <= electrons.cgrid.upar_faces.back()) {
                upar_tail[static_cast<size_t>(j)] += tp.weight;
            }
        }
    }

    {
        std::ofstream out((directory + "/energy_spectrum_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step << " bins=" << nbin
            << " emin_eV=" << emin / Const::eV
            << " emax_eV=" << emax / Const::eV << "\n";
        out << "energy_eV count_bg count_tail count_beam count_combined "
            << "dN_dE_bg_per_eV dN_dE_tail_per_eV dN_dE_beam_per_eV\n";
        for (int b = 0; b < nbin; ++b) {
            const double lo = emin * std::exp(log_range *
                static_cast<double>(b) / nbin);
            const double hi = emin * std::exp(log_range *
                static_cast<double>(b + 1) / nbin);
            const double width_eV = (hi - lo) / Const::eV;
            const double center = std::sqrt(lo * hi) / Const::eV;
            const double combined = count_bg[static_cast<size_t>(b)] +
                                    count_tail[static_cast<size_t>(b)] +
                                    count_beam[static_cast<size_t>(b)];
            out << std::setprecision(17) << center << " "
                << count_bg[static_cast<size_t>(b)] << " "
                << count_tail[static_cast<size_t>(b)] << " "
                << count_beam[static_cast<size_t>(b)] << " "
                << combined << " "
                << count_bg[static_cast<size_t>(b)] /
                   std::max(width_eV, 1.0e-300) << " "
                << count_tail[static_cast<size_t>(b)] /
                   std::max(width_eV, 1.0e-300) << " "
                << count_beam[static_cast<size_t>(b)] /
                   std::max(width_eV, 1.0e-300) << "\n";
        }
    }

    {
        std::ofstream out((directory + "/momentum_distribution_upar_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step << "\n";
        out << "u_par dN_du_par_bg dN_du_par_tail dN_du_par_beam\n";
        std::vector<double> upar_beam(static_cast<size_t>(nv), 0.0);
        for (size_t p = 0; p < beam.particles.size(); ++p) {
            const double u = beam.particles[p].px /
                             (Const::me * Const::c);
            int j = 0;
            for (int jj = 1; jj < nv; ++jj) {
                if (u > electrons.cgrid.upar_cells[static_cast<size_t>(jj)])
                    j = jj;
            }
            if (u >= electrons.cgrid.upar_faces.front() &&
                u <= electrons.cgrid.upar_faces.back())
                upar_beam[static_cast<size_t>(j)] += beam.particles[p].weight;
        }
        for (int j = 0; j < nv; ++j) {
            const double width = electrons.cgrid.upar_widths[
                static_cast<size_t>(j)];
            out << std::setprecision(17)
                << electrons.cgrid.upar_cells[static_cast<size_t>(j)] << " "
                << upar_bg[static_cast<size_t>(j)] / width << " "
                << upar_tail[static_cast<size_t>(j)] / width << " "
                << upar_beam[static_cast<size_t>(j)] / width << "\n";
        }
    }

    {
        std::ofstream out((directory + "/momentum_distribution_uperp_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
            << time << " step=" << step << "\n";
        out << "u_perp dN_du_perp_bg\n";
        for (int k = 0; k < nmu; ++k) {
            const double width = electrons.cgrid.uperp_widths[
                static_cast<size_t>(k)];
            out << std::setprecision(17)
                << electrons.cgrid.uperp_cells[static_cast<size_t>(k)] << " "
                << uperp_bg[static_cast<size_t>(k)] / width << "\n";
        }
    }
}

void VpfpDiagnostics::write_tail_stats(
    int step, double time, const Species& electrons,
    const SpatialGrid& grid, const BackgroundTailPIC* tail, int mpi_rank)
{
    (void)electrons;
    const std::string directory =
        snapshot_directory(output_dir_, step, time);
    std::ofstream out((directory + "/tail_per_cell_stats_rank" +
                       std::to_string(mpi_rank) + ".dat").c_str());
    if (!out) return;
    out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
        << time << " step=" << step
        << " ix_start=" << grid.ix_start << " nx_local=" << grid.nx_local
        << "\n";
    if (tail == NULL) {
        out << "tail_mode off\n";
        return;
    }
    out << "x_um macro_particles weight_sum_m2 density_m3\n";
    std::vector<int> counts(static_cast<size_t>(grid.nx_local), 0);
    std::vector<double> weights(static_cast<size_t>(grid.nx_local), 0.0);
    for (size_t p = 0; p < tail->particles.size(); ++p) {
        const BackgroundTailParticle& tp = tail->particles[p];
        const int cell = static_cast<int>(std::floor(tp.x / grid.dx));
        const int il = cell - grid.ix_start;
        if (il >= 0 && il < grid.nx_local) {
            ++counts[static_cast<size_t>(il)];
            weights[static_cast<size_t>(il)] += tp.weight;
        }
    }
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        out << std::setprecision(17)
            << grid.x(ix + grid.nghost) / Const::micro << " "
            << counts[static_cast<size_t>(ix)] << " "
            << weights[static_cast<size_t>(ix)] << " "
            << tail->density[static_cast<size_t>(ix)] << "\n";
    }
    const BackgroundTailPIC::WeightStats ws = tail->weight_stats(grid);
    out << "# weight_stats min=" << ws.min << " max=" << ws.max
        << " mean=" << ws.mean << " std=" << ws.std
        << " max_single_charge_fraction=" << ws.max_single_charge_fraction
        << " macro_particle_count=" << ws.macro_particle_count
        << " per_cell_max_macro_count=" << ws.per_cell_max_macro_count
        << " density_noise_estimate=" << ws.density_noise_estimate << "\n";
}

void VpfpDiagnostics::write_threshold_interface(
    int step, double time, const Species& electrons,
    const SpatialGrid& grid, const BackgroundTailPIC* tail,
    const HybridVelocityPartition* partition, double convert_energy_mev,
    int mpi_rank)
{
    const std::string directory =
        snapshot_directory(output_dir_, step, time);
    std::ofstream out((directory + "/tail_threshold_interface_rank" +
                       std::to_string(mpi_rank) + ".dat").c_str());
    if (!out) return;
    out << "# schema=vpfp-open-csl-v1 time_s=" << std::setprecision(17)
        << time << " step=" << step
        << " ix_start=" << grid.ix_start << " nx_local=" << grid.nx_local
        << "\n";
    if (tail == NULL || partition == NULL) {
        out << "tail_mode off\n";
        return;
    }
    const int nv = Param::Nv;
    const int nmu = Param::Nmu;
    const double k_out = convert_energy_mev * 1.0e6 * Const::eV;
    const double delta = 0.2 * 1.0e6 * Const::eV;
    const int ng = grid.nghost;

    // Bulk conversion-region Jx/Pixx/Piperp and the four fine bins around
    // K_out (two below, two above; section 13.2 threshold-interface
    // report).  Snapshot cadence only, never per step.
    double bulk_jx = 0.0;
    double bulk_pixx = 0.0;
    double bulk_piperp = 0.0;
    double bulk_forbidden_number = 0.0;
    double bulk_forbidden_energy = 0.0;
    double bulk_n_bin[4] = { 0.0, 0.0, 0.0, 0.0 };
    double bulk_e_bin[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const size_t xbase = static_cast<size_t>(ng + ix) *
                             static_cast<size_t>(Param::Nvmu);
        for (int j = 0; j < nv; ++j) {
            const double upar = electrons.cgrid.upar_cells[j];
            for (int k = 0; k < nmu; ++k) {
                const double uperp = electrons.cgrid.uperp_cells[k];
                const double ke = electrons.cgrid.kinetic_energy[
                    static_cast<size_t>(j) * nmu + k];
                const double mass = electrons.f[xbase +
                    static_cast<size_t>(j) * nmu + k];
                if (partition->is_conversion(j, k)) {
                    bulk_forbidden_number += mass;
                    bulk_forbidden_energy += mass * ke;
                    const double gamma =
                        std::sqrt(1.0 + upar * upar + uperp * uperp);
                    bulk_jx += -Const::qe * Const::c * mass * upar / gamma;
                    bulk_pixx += Const::me * Const::c * Const::c * mass *
                                 upar * upar / gamma;
                    bulk_piperp += Const::me * Const::c * Const::c * mass *
                                   uperp * uperp / gamma;
                }
                const int b = (ke >= k_out - 2.0 * delta &&
                               ke < k_out + 2.0 * delta)
                                  ? static_cast<int>(
                                        (ke - (k_out - 2.0 * delta)) /
                                        (4.0 * delta) * 4.0)
                                  : -1;
                if (b >= 0 && b < 4) {
                    bulk_n_bin[b] += mass;
                    bulk_e_bin[b] += mass * ke;
                }
            }
        }
    }

    double tail_jx = 0.0;
    double tail_pixx = 0.0;
    double tail_piperp = 0.0;
    double tail_n_bin[4] = { 0.0, 0.0, 0.0, 0.0 };
    double tail_e_bin[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (size_t p = 0; p < tail->particles.size(); ++p) {
        const BackgroundTailParticle& tp = tail->particles[p];
        const double gamma = std::sqrt(
            1.0 + tp.ux * tp.ux + tp.uy * tp.uy + tp.uz * tp.uz);
        const double ke =
            Const::me * Const::c * Const::c * (gamma - 1.0);
        tail_jx += -Const::qe * Const::c * tp.weight * tp.ux / gamma;
        tail_pixx += Const::me * Const::c * Const::c * tp.weight *
                     tp.ux * tp.ux / gamma;
        tail_piperp += Const::me * Const::c * Const::c * tp.weight *
                       (tp.uy * tp.uy + tp.uz * tp.uz) / gamma;
        const int b = (ke >= k_out - 2.0 * delta && ke < k_out + 2.0 * delta)
                          ? static_cast<int>((ke - (k_out - 2.0 * delta)) /
                                             (4.0 * delta) * 4.0)
                          : -1;
        if (b >= 0 && b < 4) {
            tail_n_bin[b] += tp.weight;
            tail_e_bin[b] += tp.weight * ke;
        }
    }

    out << "K_out_mev " << std::setprecision(17) << convert_energy_mev
        << "\nquantity bulk tail\n"
        << "Jx_A_m2 " << bulk_jx << " " << tail_jx << "\n"
        << "Pixx_J_m2 " << bulk_pixx << " " << tail_pixx << "\n"
        << "Piperp_J_m2 " << bulk_piperp << " " << tail_piperp << "\n"
        << "forbidden_bulk_number_m2 " << bulk_forbidden_number << "\n"
        << "forbidden_bulk_energy_J_m2 " << bulk_forbidden_energy << "\n";
    const char* bin_names[4] = {
        "Kout_minus2", "Kout_minus1", "Kout_plus1", "Kout_plus2"
    };
    for (int b = 0; b < 4; ++b) {
        const double low_mev = convert_energy_mev - 0.4 + 0.2 * b;
        const double high_mev = low_mev + 0.2;
        const double width_mev = high_mev - low_mev;
        out << bin_names[b] << "_low_mev " << low_mev << "\n"
            << bin_names[b] << "_high_mev " << high_mev << "\n"
            << bin_names[b] << "_width_mev " << width_mev << "\n";
        out << bin_names[b] << "_count " << bulk_n_bin[b] << " "
            << tail_n_bin[b] << "\n"
            << bin_names[b] << "_dN_dK_per_mev "
            << bulk_n_bin[b] / width_mev << " "
            << tail_n_bin[b] / width_mev << "\n"
            << bin_names[b] << "_energy_J " << bulk_e_bin[b] << " "
            << tail_e_bin[b] << "\n";
    }

    // Conversion source spectrum of the last accepted step (section 13.2).
    out << "conversion_source_bins "
        << last_conversion_source_spectrum_.size() << "\n";
    for (size_t b = 0; b < last_conversion_source_spectrum_.size(); ++b) {
        out << "conversion_source_bin " << b << " "
            << (b < last_conversion_source_edges_.size()
                    ? last_conversion_source_edges_[b] : 0.0) << " "
            << (b + 1 < last_conversion_source_edges_.size()
                    ? last_conversion_source_edges_[b + 1] : 0.0) << " "
            << last_conversion_source_spectrum_[b] << "\n";
    }
    out << "control_groups " << last_control_groups_
        << " control_fallbacks " << last_control_fallbacks_ << "\n";
}

void VpfpDiagnostics::write_manifest(int step, double time,
                                     const Species& electrons,
                                     const BeamPIC& beam,
                                     const SpatialGrid& grid,
                                     const BackgroundTailPIC* tail,
                                     const HybridVelocityPartition* partition,
                                     double convert_energy_mev,
                                     int mpi_rank, int mpi_size)
{
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer),
                  "%s/snapshot_t%.9gfs_step%d", output_dir_.c_str(),
                  time / Const::femto, step);
    const std::string directory(buffer);

    {
        std::ofstream out((directory + "/manifest_rank" +
                           std::to_string(mpi_rank) + ".dat").c_str());
        out << "schema=vpfp-open-csl-v1\n"
            << "version=2\n"
            << "step=" << step << "\n"
            << "time_s=" << std::setprecision(17) << time << "\n"
            << "time_fs=" << time / Const::femto << "\n"
            << "rank=" << mpi_rank << " of " << mpi_size << "\n"
            << "ix_start=" << grid.ix_start << "\n"
            << "nx_local=" << grid.nx_local << "\n"
            << "nx_global=" << grid.nx_global << "\n"
            << "dx_m=" << grid.dx << "\n"
            << "nv=" << electrons.cgrid.upar_cells.size() << "\n"
            << "nmu=" << electrons.cgrid.uperp_cells.size() << "\n"
            << "upar_extended_max="
            << electrons.cgrid.upar_faces.back() << "\n"
            << "uperp_max=" << electrons.cgrid.uperp_faces.back() << "\n"
            << "beam_particles=" << beam.particles.size() << "\n"
            << "background_representation="
            << (tail != NULL ? "eulerian_bulk_plus_pic_tail" : "eulerian_only")
            << "\n"
            << "field_boundary=" << run_config_.field_boundary << "\n"
            << "left_electric_field=" << run_config_.left_electric_field
            << "\nphi_left=" << run_config_.phi_left
            << "\nphi_right=" << run_config_.phi_right << "\n"
            << "background_boundary=" << run_config_.background_boundary
            << "\ncollision_model=" << run_config_.collision_model << "\n"
            << "collision_interface_mode="
            << run_config_.collision_interface_mode << "\n"
            << "bulk_collision_integrator="
            << run_config_.bulk_collision_integrator << "\n"
            << "collision_induced_conversion="
            << (run_config_.collision_induced_conversion ? 1 : 0) << "\n"
            << "tail_collision_kernel_requested="
            << run_config_.requested_tail_collision_kernel << "\n"
            << "tail_tail_collision_backend="
            << run_config_.tail_tail_collision_backend << "\n"
            << "tail_bulk_collision_backend="
            << run_config_.tail_bulk_collision_backend << "\n"
            << "tail_collision_weight_mode="
            << run_config_.tail_collision_weight_mode << "\n"
            << "tail_collision_weight_algorithm="
            << (run_config_.tail_collision_weight_mode == "virtual-split"
                    ? "sentoku-kemp-bounded-v1"
                    : "equal-strata-exact-v1") << "\n"
            << "tail_collision_max_substeps="
            << run_config_.tail_collision_max_substeps << "\n"
            << "tail_collision_max_particle_growth="
            << run_config_.tail_collision_max_particle_growth << "\n"
            << "collision_pair_bulk_bulk="
            << (run_config_.pair_bulk_bulk ? 1 : 0) << "\n"
            << "collision_pair_bulk_tail="
            << (run_config_.pair_bulk_tail ? 1 : 0) << "\n"
            << "collision_pair_bulk_reaction="
            << (run_config_.pair_bulk_tail ? 1 : 0) << "\n"
            << "collision_pair_tail_bulk="
            << (run_config_.pair_tail_bulk ? 1 : 0) << "\n"
            << "collision_pair_tail_tail="
            << (run_config_.pair_tail_tail ? 1 : 0) << "\n"
            << "population_control_enabled="
            << (run_config_.population_control_enabled ? 1 : 0) << "\n"
            << "population_control_interval="
            << run_config_.population_control_interval << "\n"
            << "tail_cell_moment_audit="
            << (run_config_.tail_cell_moment_audit ? 1 : 0) << "\n"
        << "tail_subcell_loading="
        << (run_config_.tail_subcell_loading ? 1 : 0) << "\n"
        << "tail_conversion_mode=" << run_config_.tail_conversion_mode << "\n"
        << "physical_config_hash=" << run_config_.physical_config_hash << "\n"
        << "diagnostic_config_hash=" << run_config_.diagnostic_config_hash << "\n"
        << "x_transport_velocity_mode="
        << run_config_.x_transport_velocity_mode << "\n"
        << "x_transport_velocity_table_schema="
        << run_config_.x_transport_velocity_table_schema << "\n"
        << "tail_flux_quadrature_order=" << run_config_.tail_flux_quadrature_order << "\n"
        << "tail_flux_max_supports=" << run_config_.tail_flux_max_supports << "\n"
        << "tail_flux_max_created_particles_per_step="
        << run_config_.tail_flux_max_created_particles_per_step << "\n"
        << "interface_topology_hash=" << run_config_.interface_topology_hash << "\n"
        << "interface_topology_version="
        << run_config_.interface_topology_version << "\n"
        << "interface_mask_hash=" << run_config_.interface_mask_hash << "\n"
        << "interface_face_list_hash="
        << run_config_.interface_face_list_hash << "\n"
        // JC4 (section 7.2/7.5): field-particle coupling config in the run
        // manifest.  Matches the checkpoint manifest keys.
        << "coupling_mode=" << run_config_.coupling_mode << "\n"
        << "coupling_max_iters=" << run_config_.coupling_max_iters << "\n"
        << "coupling_relaxation=" << std::setprecision(17)
        << run_config_.coupling_relaxation << "\n"
        << "coupling_field_tol=" << run_config_.coupling_field_tol << "\n"
        << "coupling_pairing_tol=" << run_config_.coupling_pairing_tol << "\n"
        << "background_phase_space_mode="
        << run_config_.background_phase_space_mode << "\n"
        << "x_transport_velocity_mode="
        << run_config_.x_transport_velocity_mode << "\n"
        << "x_transport_velocity_table_schema="
        << run_config_.x_transport_velocity_table_schema << "\n";
        if (tail != NULL && partition != NULL) {
            out << "tail_particles=" << tail->particles.size() << "\n"
                << "tail_id_counter=" << tail->id_counter() << "\n"
                << "tail_convert_energy_mev=" << std::setprecision(17)
                << convert_energy_mev << "\n"
                << "partition_config_hash=" << partition->config_hash
                << "\nconversion_energy_edge_count="
                << partition->conversion_energy_edges.size()
                << "\nconversion_energy_edges=";
            for (size_t i = 0;
                 i < partition->conversion_energy_edges.size(); ++i) {
                out << (i == 0 ? "" : ",") << std::setprecision(17)
                    << partition->conversion_energy_edges[i];
            }
            out << "\nconversion_energy_edges_hash="
                << partition->conversion_energy_edges_hash() << "\n"
                << "tail_particle_schema=background_tail_particle_v1\n"
                << "tail_pusher_backend=relativistic_dkd_v1\n"
                << "tail_deposition_backend=cic_v1\n"
                << "tail_collision_backend="
                << run_config_.requested_tail_collision_kernel << "\n"
                << "control_groups=" << last_control_groups_ << "\n"
                << "control_fallbacks=" << last_control_fallbacks_ << "\n"
                << "conversion_source_bins="
                << last_conversion_source_spectrum_.size() << "\n";
        }
    }

    if (mpi_rank != 0) return;
    std::ofstream out((directory + "/manifest.dat").c_str());
    out << "schema=vpfp-open-csl-v1\n"
        << "version=2\n"
        << "step=" << step << "\n"
        << "time_s=" << std::setprecision(17) << time << "\n"
        << "time_fs=" << time / Const::femto << "\n"
        << "nx_global=" << grid.nx_global << "\n"
        << "dx_m=" << grid.dx << "\n"
        << "nv=" << electrons.cgrid.upar_cells.size() << "\n"
        << "nmu=" << electrons.cgrid.uperp_cells.size() << "\n"
        << "upar_extended_max=" << electrons.cgrid.upar_faces.back() << "\n"
        << "uperp_max=" << electrons.cgrid.uperp_faces.back() << "\n"
        << "ranks=" << mpi_size << "\n"
        << "field_boundary=" << run_config_.field_boundary << "\n"
        << "collision_model=" << run_config_.collision_model << "\n"
        << "physical_config_hash=" << run_config_.physical_config_hash << "\n"
        << "diagnostic_config_hash=" << run_config_.diagnostic_config_hash << "\n"
        << "collision_interface_mode="
        << run_config_.collision_interface_mode << "\n"
        << "bulk_collision_integrator="
        << run_config_.bulk_collision_integrator << "\n"
        << "collision_induced_conversion="
        << (run_config_.collision_induced_conversion ? 1 : 0) << "\n"
        << "tail_collision_kernel_requested="
        << run_config_.requested_tail_collision_kernel << "\n"
        << "tail_tail_collision_backend="
        << run_config_.tail_tail_collision_backend << "\n"
        << "tail_bulk_collision_backend="
        << run_config_.tail_bulk_collision_backend << "\n"
        << "tail_collision_weight_mode="
        << run_config_.tail_collision_weight_mode << "\n"
        << "tail_collision_weight_algorithm="
        << (run_config_.tail_collision_weight_mode == "virtual-split"
                ? "sentoku-kemp-bounded-v1"
                : "equal-strata-exact-v1") << "\n"
        << "tail_collision_max_substeps="
        << run_config_.tail_collision_max_substeps << "\n"
        << "tail_collision_max_particle_growth="
        << run_config_.tail_collision_max_particle_growth << "\n"
        << "collision_pair_bulk_bulk="
        << (run_config_.pair_bulk_bulk ? 1 : 0) << "\n"
        << "collision_pair_bulk_tail="
        << (run_config_.pair_bulk_tail ? 1 : 0) << "\n"
        << "collision_pair_bulk_reaction="
        << (run_config_.pair_bulk_tail ? 1 : 0) << "\n"
        << "collision_pair_tail_bulk="
        << (run_config_.pair_tail_bulk ? 1 : 0) << "\n"
        << "collision_pair_tail_tail="
        << (run_config_.pair_tail_tail ? 1 : 0) << "\n"
        << "population_control_enabled="
        << (run_config_.population_control_enabled ? 1 : 0) << "\n"
        << "population_control_interval="
        << run_config_.population_control_interval << "\n"
        << "files=fields,density_background,density_beam,"
           "distribution_moments,energy_spectrum,"
           "momentum_distribution_upar,momentum_distribution_uperp,"
           "tail_per_cell_stats,tail_threshold_interface\n"
        << "per_rank_manifest=manifest_rank{rank}.dat\n";
}
