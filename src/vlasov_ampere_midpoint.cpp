#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <mpi.h>

namespace {

void fill_periodic_ghosts(Species& sp, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    for (int g = 0; g < ng; ++g) {
        const int left_source = ng + nxl - ng + g;
        const int right_source = ng + g;
        std::copy(sp.f.begin() + static_cast<size_t>(left_source) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(left_source + 1) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(g) * Param::Nvmu);
        std::copy(sp.f.begin() + static_cast<size_t>(right_source) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(right_source + 1) * Param::Nvmu,
                  sp.f.begin() + static_cast<size_t>(ng + nxl + g) * Param::Nvmu);
    }
}

}

VlasovAmpereMidpointSolver::VlasovAmpereMidpointSolver()
    : step_diagnostics_enabled_(false), beam_enabled_(true),
      low_order_only_(false), nonuniform_high_order_enabled_(false),
      fct_enabled_(true), capture_midpoint_input_(false),
      legacy_boundary_upwind_high_candidate_for_test_(false),
      max_midpoint_iterations_(20), midpoint_iteration_trace_for_test_(false),
      fct_activation_audit_enabled_(false),
      controlled_fct_flux_injection_enabled_(false),
      allow_finite_negative_debt_for_test_(false),
      core_macro_debt_consecutive_steps_(0)
{}

void VlasovAmpereMidpointSolver::reset_current_diag(CurrentDiagnostics& diag) const
{
    diag = CurrentDiagnostics();
}

void VlasovAmpereMidpointSolver::reset_result(Result& result) const
{
    result = Result();
    result.limiter_min_alpha = 1.0;
    result.limiter_min_alpha_core = 1.0;
    result.limiter_min_alpha_boundary = 1.0;
    result.mu_low_u_alpha_min = 1.0;
    result.mu_low_u_alpha_min_boundary = 1.0;
    result.mu_low_u_alpha_min_core = 1.0;
    result.alpha_core_min = 1.0;
    result.alpha_boundary_min = 1.0;
    result.alpha_tail_min = 1.0;
    result.u_force_alpha_min = 1.0;
    result.low_u_average_subcycles = 1.0;
    result.low_u_max_subcycles = 1;
    result.boundary_force_nsub_max = 1;
    result.x_low_input_min_f = std::numeric_limits<double>::infinity();
    result.x_low_output_min_f = std::numeric_limits<double>::infinity();
    result.x_final_min_f = std::numeric_limits<double>::infinity();
    result.failure_low_min = std::numeric_limits<double>::infinity();
    result.failure_final_min = std::numeric_limits<double>::infinity();
    result.fct_high_candidate_min = std::numeric_limits<double>::infinity();
    result.fct_final_scratch_min = std::numeric_limits<double>::infinity();
    result.low_order_candidate_min = std::numeric_limits<double>::infinity();
    result.low_order_negative_count = 0;
    result.low_order_negative_mass = 0.0;
    result.low_order_roundoff_zeroed_count = 0;
    result.low_order_roundoff_zeroed_mass = 0.0;
    result.fct_roundoff_zeroed_count = 0;
    result.fct_roundoff_zeroed_mass = 0.0;
    result.stage5_r_couple_centered = 0.0;
    result.stage5_r_couple_upwind_stabilization = 0.0;
    result.stage5_r_couple_fct_stabilization = 0.0;
    result.coupling_rj_global_sum = 0.0;
    result.coupling_rk_global_sum = 0.0;
    result.coupling_face_work_jn_global_sum = 0.0;
    result.coupling_rj_reconstruction_error = 0.0;
    result.coupling_rk_reconstruction_error = 0.0;
    result.coupling_face_work_jn_reconstruction_error = 0.0;
    result.coupling_beam_front_ix = -1;
    for (size_t region = 0; region < result.coupling_regions.size(); ++region) {
        result.coupling_regions[region].max_abs_jn_minus_gstar_je_face = -1;
    }
    result.fct_high_low_identity_worst_ix = -1;
    result.fct_high_low_identity_worst_iv = -1;
    result.fct_high_low_identity_worst_imu = -1;
    result.fct_donor_worst_ix = -1;
    result.fct_donor_worst_iv = -1;
    result.fct_donor_worst_imu = -1;
    result.fct_donor_beta_min = 1.0;
    result.failure_worst_ix = -1;
    result.failure_worst_iv = -1;
    result.failure_worst_imu = -1;
    result.failure_worst_is_core = -1;
    result.failure_worst_is_tail = -1;
    result.failure_low_work_input_min =
        std::numeric_limits<double>::infinity();
    result.accepted_transport.min_mass =
        std::numeric_limits<double>::infinity();
    result.accepted_transport.ix = -1;
    result.accepted_transport.iv = -1;
    result.accepted_transport.imu = -1;
    for (int f = 0; f < 4; ++f) {
        result.accepted_transport.alpha[f] = 1.0;
    }
    for (int d = 0; d < 3; ++d) {
        result.flux_pos[d].min_f_before = std::numeric_limits<double>::infinity();
        result.flux_pos[d].min_f_low = std::numeric_limits<double>::infinity();
        result.flux_pos[d].min_f_final = std::numeric_limits<double>::infinity();
        result.flux_pos[d].alpha_min = 1.0;
    }
}

void VlasovAmpereMidpointSolver::exchange_ghosts_x_persistent(
    Species& sp, const SpatialGrid& sg, int mpi_rank, int mpi_size) const
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t count = static_cast<size_t>(ng) * Param::Nvmu;
    if (mpi_size == 1) {
        fill_periodic_ghosts(sp, sg);
        return;
    }
    if (ghost_send_left_.size() != count) {
        ghost_send_left_.resize(count);
        ghost_send_right_.resize(count);
        ghost_recv_left_.resize(count);
        ghost_recv_right_.resize(count);
    }
    std::memcpy(ghost_send_left_.data(),
                sp.f.data() + static_cast<size_t>(ng) * Param::Nvmu,
                count * sizeof(double));
    std::memcpy(ghost_send_right_.data(),
                sp.f.data() + static_cast<size_t>(ng + nxl - ng) * Param::Nvmu,
                count * sizeof(double));
    const int left = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right = (mpi_rank + 1) % mpi_size;
    MPI_Request reqs[4];
    MPI_Isend(ghost_send_left_.data(), static_cast<int>(count), MPI_DOUBLE,
              left, 501, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(ghost_recv_left_.data(), static_cast<int>(count), MPI_DOUBLE,
              left, 502, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(ghost_send_right_.data(), static_cast<int>(count), MPI_DOUBLE,
              right, 502, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(ghost_recv_right_.data(), static_cast<int>(count), MPI_DOUBLE,
              right, 501, MPI_COMM_WORLD, &reqs[3]);
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
    std::memcpy(sp.f.data(), ghost_recv_left_.data(), count * sizeof(double));
    std::memcpy(sp.f.data() + static_cast<size_t>(ng + nxl) * Param::Nvmu,
                ghost_recv_right_.data(), count * sizeof(double));
}

void VlasovAmpereMidpointSolver::set_midpoint_field(
    EMFields& fields_mid, const EMFields& fields_n,
    const std::vector<double>& ex_mid, const SpatialGrid& sg,
    int mpi_rank, int mpi_size) const
{
    fields_mid = fields_n;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        fields_mid.Ex_face[static_cast<size_t>(iface)] = ex_mid[static_cast<size_t>(iface)];
    }
    if (sg.nx_local > 0) fields_mid.Ex_face[static_cast<size_t>(sg.nx_local)] =
        fields_mid.Ex_face[0];
    fields_mid.sync_cell_ex_from_faces(mpi_rank, mpi_size);
}

double VlasovAmpereMidpointSolver::integrate_face_work(
    const std::vector<double>& current_face, const EMFields& fields_mid,
    const SpatialGrid& sg, double dt) const
{
    double work = 0.0;
    const int n = std::min(sg.nx_local, static_cast<int>(std::min(
        current_face.size(), fields_mid.Ex_face.size())));
    for (int iface = 0; iface < n; ++iface) {
        work += current_face[static_cast<size_t>(iface)] *
                fields_mid.Ex_face[static_cast<size_t>(iface)] * sg.dx;
    }
    return dt * work;
}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_background_and_fields(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size)
{
    if (!bkg_n.cylindrical_mass_representation) {
        Result failed;
        reset_result(failed);
        failed.failed = true;
        return failed;
    }
    if (!low_order_only_ && !bkg_n.cgrid.is_uniform() &&
        !nonuniform_high_order_enabled_) {
        Result failed;
        reset_result(failed);
        failed.failed = true;
        failed.state_advanced = 0;
        failed.failure_reason = 12;
        if (mpi_rank == 0) {
            std::cerr << "VlasovAmpereMidpointSolver configuration error: "
                      << "nonuniform cylindrical velocity grid requires "
                      << "set_nonuniform_high_order_enabled(true) when "
                      << "low_order_only=false\n";
        }
        return failed;
    }
    Result result = evaluate_production_midpoint_operator(
        bkg_n, beam_n, fields_n, sg, dt, time, mpi_rank, mpi_size, 1);
    if (result.state_advanced && !result.failed) {
        const bool core_macro_debt =
            result.neg_mass_core_fraction > 1.0e-6 ||
            result.neg_energy_core_fraction > 1.0e-6 ||
            result.neg_current_core_fraction > 1.0e-6;
        core_macro_debt_consecutive_steps_ = core_macro_debt
            ? core_macro_debt_consecutive_steps_ + 1 : 0;
        if (core_macro_debt_consecutive_steps_ >= 3) {
            result.failed = true;
            result.state_advanced = 0;
        }
    }
    return result;
}

VlasovAmpereMidpointSolver::MidpointOperatorEvaluation
VlasovAmpereMidpointSolver::evaluate_fixed_midpoint_operator(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const Species& guess_np1, const EMFields& fields_end_guess,
    const std::vector<double>& fixed_j_beam_face_mid,
    const CouplingRegionLayout& coupling_layout,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size) const
{
    return evaluate_production_midpoint_operator(
        bkg_n, beam_n, fields_n, sg, dt, time, mpi_rank, mpi_size, 1,
        &guess_np1, &fields_end_guess, &fixed_j_beam_face_mid,
        &coupling_layout, true);
}
