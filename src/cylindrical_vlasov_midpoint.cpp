#include "vlasov_ampere_midpoint.h"
#include "ppm_ctu_reconstruction.h"
#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <omp.h>

namespace {

inline size_t mass_index(int ix, int j, int k)
{
    return idx3(ix, j, k);
}

inline size_t xface_index(int iface, int j, int k)
{
    return (static_cast<size_t>(iface) * Param::Nv + j) * Param::Nmu + k;
}

inline size_t uface_index(int ix, int jface, int k)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1) + jface) *
           Param::Nmu + k;
}

inline double candidate_roundoff_tolerance(double m_low, double ax_left,
                                           double ax_right, double au_lower,
                                           double au_upper)
{
    // Local final-scratch summation bound.  The subnormal term is explicit:
    // in an empty tail cell all visible transport terms can underflow to zero,
    // while the reconstructed value can still be a few denormals negative.
    const double local_sum = std::fabs(m_low) + std::fabs(ax_left) +
        std::fabs(ax_right) + std::fabs(au_lower) + std::fabs(au_upper);
    return 64.0 * std::numeric_limits<double>::epsilon() * local_sum +
        64.0 * std::numeric_limits<double>::denorm_min();
}

inline bool normalize_roundoff_negative_candidate(double& candidate,
                                                  double roundoff_tolerance,
                                                  double& normalized_mass)
{
    if (!std::isfinite(candidate)) return false;
    // Canonicalize signed zero as well.  It carries no mass and must not
    // survive as a negative-looking scratch value in downstream diagnostics.
    if (candidate == 0.0) {
        candidate = 0.0;
        return false;
    }
    if (candidate > 0.0) return false;
    if (-candidate > roundoff_tolerance) return false;
    normalized_mass = -candidate;
    candidate = 0.0; // Canonical positive zero, not a global distribution clip.
    return true;
}

void close_periodic_face_blocks(std::vector<double>& faces, int nxl,
                                int block, int rank, int size, int tag)
{
    if (nxl <= 0 || block <= 0) return;
    const size_t first = 0;
    const size_t last = static_cast<size_t>(nxl) * block;
    if (size == 1) {
        std::copy(faces.begin() + first, faces.begin() + first + block,
                  faces.begin() + last);
        return;
    }
    const int left = (rank + size - 1) % size;
    const int right = (rank + 1) % size;
    MPI_Sendrecv(faces.data() + first, block, MPI_DOUBLE, left, tag,
                 faces.data() + last, block, MPI_DOUBLE, right, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void exchange_neighbor_scalars(const std::vector<double>& local, int rank,
                               int size, double& left_value,
                               double& right_value, int tag)
{
    if (local.empty()) {
        left_value = right_value = 1.0;
        return;
    }
    if (size == 1) {
        left_value = local.back();
        right_value = local.front();
        return;
    }
    const int left = (rank + size - 1) % size;
    const int right = (rank + 1) % size;
    MPI_Sendrecv(&local.front(), 1, MPI_DOUBLE, left, tag,
                 &right_value, 1, MPI_DOUBLE, right, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&local.back(), 1, MPI_DOUBLE, right, tag + 1,
                 &left_value, 1, MPI_DOUBLE, left, tag + 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

bool all_finite(const Species& sp, const EMFields& fields,
                const std::vector<double>& current)
{
    for (size_t i = 0; i < sp.f.size(); ++i)
        if (!std::isfinite(sp.f[i])) return false;
    for (size_t i = 0; i < fields.Ex_face.size(); ++i)
        if (!std::isfinite(fields.Ex_face[i])) return false;
    for (size_t i = 0; i < current.size(); ++i)
        if (!std::isfinite(current[i])) return false;
    return true;
}

}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_cylindrical_single_step(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size, int substeps_used) const
{
    Result result;
    reset_result(result);
    result.substeps_used = substeps_used;
    result.species_np1 = bkg_n;
    result.beam_np1 = beam_n;
    result.fields_np1 = fields_n;

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t nface = static_cast<size_t>(nxl + 1);
    const size_t xflux_size = nface * Param::Nvmu;
    const size_t uflux_size = static_cast<size_t>(nxl) *
                              (Param::Nv + 1) * Param::Nmu;
    const auto set_failure_location = [&](int global_linear) {
        result.failure_worst_ix = global_linear / Param::Nvmu;
        const int velocity_linear = global_linear % Param::Nvmu;
        result.failure_worst_iv = velocity_linear / Param::Nmu;
        result.failure_worst_imu = velocity_linear % Param::Nmu;
    };
    // Reused across nonlinear iterations/substeps.  These are the dominant
    // transient allocations in the new FV/FCT kernel.
    std::vector<double> fx_low(xflux_size, 0.0), fx_high(xflux_size, 0.0);
    std::vector<double> fu_low(uflux_size, 0.0), fu_high(uflux_size, 0.0);
    std::vector<double> cu_low(uflux_size, 0.0), cu_high(uflux_size, 0.0);
    const size_t cell_count = static_cast<size_t>(sg.nx_total) * Param::Nvmu;
    const size_t local_cell_count = static_cast<size_t>(nxl) * Param::Nvmu;
    // Full PPM/CTU work storage.  These arrays are reused for every nonlinear
    // iteration and substep; no per-face temporary allocation is performed.
    std::vector<PpmCtu::Parabola> ppm_x(cell_count), ppm_u(local_cell_count);
    std::vector<double> x_left_1d(xflux_size, 0.0), x_right_1d(xflux_size, 0.0);
    std::vector<double> u_left_1d(uflux_size, 0.0), u_right_1d(uflux_size, 0.0);
    std::vector<double> fx_provisional(xflux_size, 0.0);
    std::vector<double> fu_provisional(uflux_size, 0.0);
    std::vector<double> alpha_x_left(local_cell_count, 1.0);
    std::vector<double> alpha_x_right(local_cell_count, 1.0);
    std::vector<double> alpha_u_lower(local_cell_count, 1.0);
    std::vector<double> alpha_u_upper(local_cell_count, 1.0);
    std::vector<double> left_alpha_x_right(Param::Nvmu, 1.0);
    std::vector<double> fx_final(xflux_size, 0.0), fu_final(uflux_size, 0.0);
    std::vector<double> cu_final(uflux_size, 0.0);
    std::vector<double> candidate_mass(local_cell_count, 0.0);
    std::vector<double> donor_beta(local_cell_count, 1.0);
    Species low_state_buffer = bkg_n;
    std::vector<double> e_mid_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> next_e_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_jn_buffer(nface, 0.0);
    std::vector<double> integrated_je_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_cell_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> initial_ke_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> initial_p_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> final_ke_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> final_p_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> delta_ke_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> delta_p_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> local_energy_residual(static_cast<size_t>(nxl), 0.0);
    std::vector<double> local_momentum_residual(static_cast<size_t>(nxl), 0.0);
    std::vector<double> psi_k_x(nface, 0.0), psi_p_x(nface, 0.0);
    std::vector<double> u_energy_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> u_momentum_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> u_boundary_energy_cell(static_cast<size_t>(nxl), 0.0);
    std::vector<double> u_boundary_momentum_cell(static_cast<size_t>(nxl), 0.0);
    result.j_bkg_face_mid.assign(nface, 0.0);
    result.j_beam_face_mid.assign(nface, 0.0);
    result.j_total_face_mid.assign(nface, 0.0);
    result.j_bkg_energy_debug_face.assign(nface, 0.0);

    std::vector<double> e_end(static_cast<size_t>(nxl), 0.0);
    for (int iface = 0; iface < nxl; ++iface) {
        e_end[static_cast<size_t>(iface)] = fields_n.Ex_face[iface];
    }
    Species guess = bkg_n;
    Species last_work = bkg_n;
    EMFields last_fields = fields_n;
    std::vector<double> previous_j(static_cast<size_t>(nxl), 0.0);
    bool have_previous = false;

    // Stage 7 is deliberately not enabled: the beam is a controlled
    // predictor during phases 1--5, and its lag is kept explicit.
    BeamPIC beam_predictor = beam_n;
    EMFields beam_field = fields_n;
    beam_predictor.begin_step(sg, dt);
    beam_predictor.inject(sg, beam_field, dt, time, mpi_rank, mpi_size);
    const double beam_ke_before = beam_predictor.total_kinetic_energy();
    beam_predictor.push(sg, beam_field, dt, mpi_rank, mpi_size);
    const double beam_ke_after = beam_predictor.total_kinetic_energy();
    beam_predictor.deposit_density(sg, mpi_rank, mpi_size);
    beam_predictor.finalize_charge_conserving_current(sg, dt, mpi_rank,
                                                      mpi_size);
    std::vector<double> jbeam(nface, 0.0);
    for (size_t i = 0; i < std::min(jbeam.size(),
                                    beam_predictor.current_face_x.size()); ++i) {
        jbeam[i] = beam_predictor.current_face_x[i];
    }

    const int max_iters = 20;
    const double field_tol = 1.0e-6;
    const double current_tol = 1.0e-5;
    const double omega = 0.55;
    double initial_number = 0.0;
    double initial_ke = 0.0;
    bkg_n.total_particle_number_and_energy(initial_number, initial_ke);

    const auto write_stage5_diagnostics = [&](int accepted, int soft_unconverged) {
        static int stage5_write_count = 0;
        const bool flush_now = (++stage5_write_count % 20) == 0;
        if (mpi_rank == 0) {
            static std::ofstream global_file("output/background_closure_diagnostics.dat");
            if (global_file.tellp() == std::streampos(0)) {
                global_file << "# time_fs accepted soft R_FV R_couple R_total R_physical dK face_work "
                            << "u_energy u_momentum u_energy_outflow PsiK_boundary PsiP_boundary "
                            << "JN_minus_JE_L2 JN_minus_JE_Linf vEC_Linf vEC_rel "
                            << "fct_budget_violation correction_enabled skip_zero_field skip_raw "
                            << "correction_target correction_achieved correction_added_momentum_linf "
                            << "correction_l1 correction_linf correction_relative active_bounds infeasible retry_count "
                            << "fct_high_low_linf fct_high_low_violation "
                            << "fct_donor_capacity_violation fct_interface_checksum_linf "
                            << "fct_interface_checksum_violation fct_low_tolerance_linf "
                            << "fct_final_scratch_min fct_final_tolerance_linf "
                            << "fct_roundoff_zeroed_count fct_roundoff_zeroed_mass "
                            << "fct_id_R_worst fct_id_S_worst fct_id_R_over_S "
                            << "fct_id_abs_R_over_max1S "
                            << "fct_id_worst_ix fct_id_worst_iv fct_id_worst_imu "
                            << "donor_m_low donor_Ax_left donor_Ax_right "
                            << "donor_Au_lower donor_Au_upper "
                            << "donor_alpha_x_left donor_alpha_x_right "
                            << "donor_alpha_u_lower donor_alpha_u_upper "
                            << "donor_outflow donor_scale donor_relative "
                            << "donor_ix donor_iv donor_imu donor_roundoff_warning "
                            << "donor_beta_count donor_beta_min\n";
                global_file << std::scientific << std::setprecision(10);
            }
            global_file << time / Const::femto << " " << accepted << " " << soft_unconverged << " "
                        << result.stage5_r_fv << " " << result.stage5_r_couple << " "
                        << result.stage5_r_total << " " << result.stage5_r_physical_balance << " "
                        << result.delta_ke_bkg << " "
                        << result.field_work_bkg << " " << result.stage5_u_energy_moment << " "
                        << result.stage5_u_momentum_moment << " "
                        << result.u_boundary_energy_outflow << " "
                        << result.stage5_spatial_energy_boundary << " "
                        << result.stage5_spatial_momentum_boundary << " "
                        << result.stage5_jn_minus_je_l2 << " "
                        << result.stage5_jn_minus_je_linf << " "
                        << result.stage5_energy_speed_candidate_linf << " "
                        << result.stage5_energy_speed_candidate_rel << " "
                        << result.stage5_fct_budget_violation << " "
                        << result.stage5_compatibility_enabled << " "
                        << result.stage5_projection_skipped_zero_field << " "
                        << result.stage5_projection_skipped_raw_residual << " "
                        << result.stage5_correction_energy_target << " "
                        << result.stage5_correction_energy_achieved << " "
                        << result.stage5_correction_added_momentum_linf << " "
                        << result.stage5_correction_l1 << " "
                        << result.stage5_correction_linf << " "
                        << result.stage5_correction_relative << " "
                        << result.stage5_correction_active_bounds << " "
                        << result.stage5_correction_infeasible << " "
                        << result.stage5_correction_retry_count << " "
                        << result.fct_high_low_identity_linf << " "
                        << result.fct_high_low_identity_violation << " "
                        << result.fct_donor_capacity_violation << " "
                        << result.fct_interface_checksum_linf << " "
                        << result.fct_interface_checksum_violation << " "
                        << result.fct_low_order_tolerance_linf << " "
                        << result.fct_final_scratch_min << " "
                        << result.fct_final_tolerance_linf << " "
                        << result.fct_roundoff_zeroed_count << " "
                        << result.fct_roundoff_zeroed_mass << " "
                        << result.fct_high_low_identity_worst_residual << " "
                        << result.fct_high_low_identity_worst_scale << " "
                        << result.fct_high_low_identity_worst_relative << " "
                        << result.fct_high_low_identity_ratio_linf << " "
                        << result.fct_high_low_identity_worst_ix << " "
                        << result.fct_high_low_identity_worst_iv << " "
                        << result.fct_high_low_identity_worst_imu << " "
                        << result.fct_donor_worst_m_low << " "
                        << result.fct_donor_worst_ax_left << " "
                        << result.fct_donor_worst_ax_right << " "
                        << result.fct_donor_worst_au_lower << " "
                        << result.fct_donor_worst_au_upper << " "
                        << result.fct_donor_worst_alpha_x_left << " "
                        << result.fct_donor_worst_alpha_x_right << " "
                        << result.fct_donor_worst_alpha_u_lower << " "
                        << result.fct_donor_worst_alpha_u_upper << " "
                        << result.fct_donor_worst_outflow << " "
                        << result.fct_donor_worst_scale << " "
                        << result.fct_donor_worst_relative << " "
                        << result.fct_donor_worst_ix << " "
                        << result.fct_donor_worst_iv << " "
                        << result.fct_donor_worst_imu << " "
                        << result.fct_donor_roundoff_warning << " "
                        << result.fct_donor_beta_applied_count << " "
                        << result.fct_donor_beta_min << "\n";
            if (flush_now) global_file.flush();
        }

        // Keep the detailed closure ledger compact during long production
        // runs.  The global closure record remains per accepted step, while
        // this rank summary is sampled every 50 accepted states.
        const int local_summary_interval = 50;
        if (stage5_write_count % local_summary_interval != 0) return;

        int worst_energy_ix = -1;
        int worst_momentum_ix = -1;
        double worst_energy_abs = -1.0;
        double worst_momentum_abs = -1.0;
        double local_energy_l1 = 0.0;
        double local_momentum_l1 = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            const double energy_abs = std::fabs(local_energy_residual[ix]);
            const double momentum_abs = std::fabs(local_momentum_residual[ix]);
            local_energy_l1 += energy_abs;
            local_momentum_l1 += momentum_abs;
            if (energy_abs > worst_energy_abs) {
                worst_energy_abs = energy_abs;
                worst_energy_ix = sg.ix_start + ix;
            }
            if (momentum_abs > worst_momentum_abs) {
                worst_momentum_abs = momentum_abs;
                worst_momentum_ix = sg.ix_start + ix;
            }
        }
        const int values_per_rank = 7;
        const double local_summary[values_per_rank] = {
            static_cast<double>(mpi_rank),
            static_cast<double>(worst_energy_ix), worst_energy_abs,
            static_cast<double>(worst_momentum_ix), worst_momentum_abs,
            local_energy_l1, local_momentum_l1
        };
        std::vector<double> gathered_rows;
        if (mpi_rank == 0) {
            gathered_rows.resize(static_cast<size_t>(mpi_size) * values_per_rank);
        }
        MPI_Gather(local_summary, values_per_rank, MPI_DOUBLE,
                   mpi_rank == 0 ? gathered_rows.data() : 0,
                   values_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            static std::ofstream local_file("output/background_closure_rank_summary.dat");
            if (local_file.tellp() == std::streampos(0)) {
                local_file << "# time_fs accepted soft rank worst_energy_ix "
                           << "max_abs_local_energy worst_momentum_ix "
                           << "max_abs_local_momentum local_energy_L1 "
                           << "local_momentum_L1\n";
                local_file << std::scientific << std::setprecision(10);
            }
            for (size_t offset = 0; offset < gathered_rows.size(); offset += values_per_rank) {
                const double* row = gathered_rows.data() + offset;
                local_file << time / Const::femto << " " << accepted << " " << soft_unconverged << " ";
                for (int col = 0; col < values_per_rank; ++col)
                    local_file << row[col] << (col + 1 == values_per_rank ? '\n' : ' ');
            }
            if (flush_now) local_file.flush();
        }
    };

    for (int iter = 0; iter < max_iters; ++iter) {
        std::vector<double>& e_mid = e_mid_buffer;
        for (int iface = 0; iface < nxl; ++iface) {
            e_mid[iface] = 0.5 * (fields_n.Ex_face[iface] + e_end[iface]);
        }
        EMFields fields_mid;
        set_midpoint_field(fields_mid, fields_n, e_mid, sg, mpi_rank, mpi_size);

        double local_cfl = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            const double ex = fields_mid.Ex[ng + ix];
            const double a = bkg_n.charge * ex / (bkg_n.mass * Const::c);
            for (int j = 0; j < Param::Nv; ++j) {
                const double cu = dt * std::fabs(a) /
                    std::max(bkg_n.cgrid.upar_widths[j],
                             std::numeric_limits<double>::min());
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double cx = dt * std::fabs(bkg_n.cgrid.vx[idx2(j, k)]) /
                                      sg.dx;
                    local_cfl = std::max(local_cfl, cx + cu);
                }
            }
        }
        double global_cfl = local_cfl;
        MPI_Allreduce(MPI_IN_PLACE, &global_cfl, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        const int nsub = std::max(1, static_cast<int>(std::ceil(global_cfl / 0.8)));
        if (nsub > 8 || !std::isfinite(global_cfl)) {
            result.failed = true;
            result.state_advanced = 0;
            result.failure_reason = 1;
            result.failure_iteration = iter;
            result.failure_global_cfl = global_cfl;
            result.x_low_failure_kind = 3; // genuine combined CFL failure
            return result;
        }

        Species work = bkg_n;
        std::vector<double>& integrated_jn = integrated_jn_buffer;
        std::vector<double>& integrated_je_cell = integrated_je_buffer;
        std::fill(integrated_jn.begin(), integrated_jn.end(), 0.0);
        std::fill(integrated_je_cell.begin(), integrated_je_cell.end(), 0.0);
        std::fill(psi_k_x.begin(), psi_k_x.end(), 0.0);
        std::fill(psi_p_x.begin(), psi_p_x.end(), 0.0);
        std::fill(u_energy_cell.begin(), u_energy_cell.end(), 0.0);
        std::fill(u_momentum_cell.begin(), u_momentum_cell.end(), 0.0);
        std::fill(u_boundary_energy_cell.begin(), u_boundary_energy_cell.end(), 0.0);
        std::fill(u_boundary_momentum_cell.begin(), u_boundary_momentum_cell.end(), 0.0);
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            double ke = 0.0;
            double pp = 0.0;
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const double mass = bkg_n.f[mass_index(ng + ix, j, k)];
                ke += bkg_n.cgrid.kinetic_energy[idx2(j, k)] * mass;
                pp += bkg_n.mass * Const::c * bkg_n.cgrid.upar_cells[j] * mass;
            }
            initial_ke_cell[ix] = ke;
            initial_p_cell[ix] = pp;
        }
        double limiter_energy_change = 0.0;
        double limiter_energy_positive = 0.0;
        double limiter_energy_negative = 0.0;
        double limiter_energy_core = 0.0;
        double limiter_energy_boundary = 0.0;
        double u_boundary_energy = 0.0;
        double u_boundary_particle = 0.0;
        double u_boundary_momentum = 0.0;
        double limiter_faces = 0.0;
        double limiter_active = 0.0;
        double limiter_min = 1.0;
        double fct_budget_violation = 0.0;

        const double h = dt / nsub;
        const auto capture_accepted_transport =
            [&](const Species& committed) {
                double local_min = std::numeric_limits<double>::infinity();
                #pragma omp parallel for reduction(min:local_min)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k)
                            local_min = std::min(local_min,
                                committed.f[mass_index(ng + ix, j, k)]);
                struct MassMinLoc {
                    double value;
                    int rank;
                };
                MassMinLoc local_owner = {local_min, mpi_rank};
                MassMinLoc global_owner = {0.0, 0};
                MPI_Allreduce(&local_owner, &global_owner, 1, MPI_DOUBLE_INT,
                              MPI_MINLOC, MPI_COMM_WORLD);

                // global index, M_low, M_cand, M_safe, O, I, beta,
                // alpha[4], A_lim[4], A_safe[4], interface, k=0, donor.
                double record[22] = {0.0};
                if (mpi_rank == global_owner.rank) {
                    bool found = false;
                    for (int ix = 0; ix < nxl && !found; ++ix)
                        for (int j = 0; j < Param::Nv && !found; ++j)
                            for (int k = 0; k < Param::Nmu; ++k) {
                                const double m_safe = committed.f[
                                    mass_index(ng + ix, j, k)];
                                if (m_safe != global_owner.value) continue;
                                const size_t p = static_cast<size_t>(ix) *
                                    Param::Nvmu + idx2(j, k);
                                const double a_lim[4] = {
                                    h * (fx_final[xface_index(ix, j, k)] -
                                         fx_low[xface_index(ix, j, k)]),
                                    h * (fx_final[xface_index(ix + 1, j, k)] -
                                         fx_low[xface_index(ix + 1, j, k)]),
                                    h * (fu_final[uface_index(ix, j, k)] -
                                         fu_low[uface_index(ix, j, k)]),
                                    h * (fu_final[uface_index(ix, j + 1, k)] -
                                         fu_low[uface_index(ix, j + 1, k)])
                                };
                                const double a_raw[4] = {
                                    h * (fx_high[xface_index(ix, j, k)] -
                                         fx_low[xface_index(ix, j, k)]),
                                    h * (fx_high[xface_index(ix + 1, j, k)] -
                                         fx_low[xface_index(ix + 1, j, k)]),
                                    h * (fu_high[uface_index(ix, j, k)] -
                                         fu_low[uface_index(ix, j, k)]),
                                    h * (fu_high[uface_index(ix, j + 1, k)] -
                                         fu_low[uface_index(ix, j + 1, k)])
                                };
                                record[0] = static_cast<double>(
                                    (sg.ix_start + ix) * Param::Nvmu + idx2(j, k));
                                record[1] = low_state_buffer.f[
                                    mass_index(ng + ix, j, k)];
                                record[2] = candidate_mass[p];
                                record[3] = m_safe;
                                record[4] = std::max(0.0, -a_lim[0]) +
                                    std::max(0.0, a_lim[1]) +
                                    std::max(0.0, -a_lim[2]) +
                                    std::max(0.0, a_lim[3]);
                                record[5] = std::max(0.0, a_lim[0]) +
                                    std::max(0.0, -a_lim[1]) +
                                    std::max(0.0, a_lim[2]) +
                                    std::max(0.0, -a_lim[3]);
                                record[6] = donor_beta[p];
                                for (int f = 0; f < 4; ++f) {
                                    record[7 + f] = (a_raw[f] != 0.0)
                                        ? a_lim[f] / a_raw[f] : 1.0;
                                    record[11 + f] = a_lim[f];
                                    record[15 + f] = a_lim[f];
                                }
                                record[19] = (mpi_size > 1 &&
                                    (ix == 0 || ix == nxl - 1)) ? 1.0 : 0.0;
                                record[20] = (k == 0) ? 1.0 : 0.0;
                                record[21] = 1.0;
                                found = true;
                            }
                }
                MPI_Bcast(record, 22, MPI_DOUBLE, global_owner.rank,
                          MPI_COMM_WORLD);
                result.accepted_transport.valid = 1;
                result.accepted_transport.substep = nsub;
                result.accepted_transport.min_mass = global_owner.value;
                const int global_linear = static_cast<int>(record[0]);
                result.accepted_transport.ix =
                    global_linear / Param::Nvmu;
                const int velocity_linear = global_linear % Param::Nvmu;
                result.accepted_transport.iv =
                    velocity_linear / Param::Nmu;
                result.accepted_transport.imu =
                    velocity_linear % Param::Nmu;
                result.accepted_transport.m_low = record[1];
                result.accepted_transport.m_candidate = record[2];
                result.accepted_transport.m_safe = record[3];
                result.accepted_transport.outflow = record[4];
                result.accepted_transport.inflow = record[5];
                result.accepted_transport.beta = record[6];
                for (int f = 0; f < 4; ++f) {
                    result.accepted_transport.alpha[f] = record[7 + f];
                    result.accepted_transport.a_limited[f] = record[11 + f];
                    result.accepted_transport.a_safe[f] = record[15 + f];
                }
                result.accepted_transport.mpi_interface =
                    static_cast<int>(record[19]);
                result.accepted_transport.k_zero =
                    static_cast<int>(record[20]);
                result.accepted_transport.next_step_donor =
                    static_cast<int>(record[21]);
            };
        const auto transport_safe_step_gate =
            [&](Species& candidate_state, int completed_substeps) {
                // Synchronize first so the next physical step cannot obtain
                // an unchecked ghost donor.  Physical cells are then enough
                // for the global minimum because ghosts are copies of them.
                exchange_ghosts_x_persistent(candidate_state, sg, mpi_rank,
                                             mpi_size);
                double local_min = std::numeric_limits<double>::infinity();
                double local_nonfinite = 0.0;
                #pragma omp parallel for reduction(min:local_min) reduction(max:local_nonfinite)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double m = candidate_state.f[
                                mass_index(ng + ix, j, k)];
                            local_min = std::min(local_min, m);
                            if (!std::isfinite(m)) local_nonfinite = 1.0;
                        }
                struct GateMinLoc {
                    double value;
                    int rank;
                };
                GateMinLoc local_owner = {local_min, mpi_rank};
                GateMinLoc global_owner = {0.0, 0};
                MPI_Allreduce(&local_owner, &global_owner, 1, MPI_DOUBLE_INT,
                              MPI_MINLOC, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &local_nonfinite, 1, MPI_DOUBLE,
                              MPI_MAX, MPI_COMM_WORLD);
                if (local_nonfinite == 0.0 && global_owner.value >= 0.0)
                    return true;

                int global_linear = -1;
                if (mpi_rank == global_owner.rank) {
                    for (int ix = 0; ix < nxl && global_linear < 0; ++ix)
                        for (int j = 0; j < Param::Nv && global_linear < 0; ++j)
                            for (int k = 0; k < Param::Nmu; ++k)
                                if (candidate_state.f[mass_index(ng + ix, j, k)]
                                    == global_owner.value) {
                                    global_linear =
                                        (sg.ix_start + ix) * Param::Nvmu +
                                        idx2(j, k);
                                    break;
                                }
                }
                MPI_Bcast(&global_linear, 1, MPI_INT, global_owner.rank,
                          MPI_COMM_WORLD);
                result.failed = true;
                result.state_advanced = 0;
                result.failure_reason =
                    (local_nonfinite != 0.0) ? 4 : 3;
                result.failure_iteration = iter;
                result.failure_substep = completed_substeps;
                result.failure_global_cfl = global_cfl;
                result.failure_final_min = global_owner.value;
                if (global_linear >= 0)
                    set_failure_location(global_linear);
                return false;
            };
        for (int sub = 0; sub < nsub; ++sub) {
            if (sub == 0)
                exchange_ghosts_x_persistent(work, sg, mpi_rank, mpi_size);
            std::fill(fx_low.begin(), fx_low.end(), 0.0);
            std::fill(fx_high.begin(), fx_high.end(), 0.0);
            std::fill(fu_low.begin(), fu_low.end(), 0.0);
            std::fill(fu_high.begin(), fu_high.end(), 0.0);
            std::fill(cu_low.begin(), cu_low.end(), 0.0);
            std::fill(cu_high.begin(), cu_high.end(), 0.0);

            // Unsplit donor-cell baseline.  It shares the same start state,
            // substep, face topology, velocities and accelerations as PPM-CTU.
            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                const int il = ng + iface - 1;
                const int ir = ng + iface;
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const double v = bkg_n.cgrid.vx[idx2(j, k)];
                        const int donor = (v >= 0.0) ? il : ir;
                        fx_low[xface_index(iface, j, k)] =
                            v * work.f[mass_index(donor, j, k)] / sg.dx;
                    }
                }
            }
            close_periodic_face_blocks(fx_low, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 901);
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                const double a = bkg_n.charge * fields_mid.Ex[ng + ix] /
                                 (bkg_n.mass * Const::c);
                for (int jf = 1; jf < Param::Nv; ++jf) {
                    const double du = 0.5 * (bkg_n.cgrid.upar_widths[jf - 1] +
                                             bkg_n.cgrid.upar_widths[jf]);
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const int donor = (a >= 0.0) ? jf - 1 : jf;
                        const double coefficient =
                            work.f[mass_index(ng + ix, donor, k)] / du;
                        const size_t id = uface_index(ix, jf, k);
                        cu_low[id] = coefficient;
                        fu_low[id] = a * coefficient;
                    }
                }
            }

            // A--C: reconstruct complete PPM parabolas, then integrate their
            // swept half-step intervals.  The two directions use the same
            // substep-start state; neither is reconstructed from an updated
            // state of the other direction.
            #pragma omp parallel for schedule(static)
            for (int i = 2; i <= sg.nx_total - 3; ++i) {
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        ppm_x[mass_index(i, j, k)] = PpmCtu::reconstruct(
                            work.f[mass_index(i - 2, j, k)],
                            work.f[mass_index(i - 1, j, k)],
                            work.f[mass_index(i, j, k)],
                            work.f[mass_index(i + 1, j, k)],
                            work.f[mass_index(i + 2, j, k)]);
                    }
                }
            }
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) {
                    const int jm2 = std::max(0, j - 2);
                    const int jm1 = std::max(0, j - 1);
                    const int jp1 = std::min(Param::Nv - 1, j + 1);
                    const int jp2 = std::min(Param::Nv - 1, j + 2);
                    for (int k = 0; k < Param::Nmu; ++k) {
                        ppm_u[static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k)] =
                            PpmCtu::reconstruct(work.f[mass_index(ng + ix, jm2, k)],
                                               work.f[mass_index(ng + ix, jm1, k)],
                                               work.f[mass_index(ng + ix, j, k)],
                                               work.f[mass_index(ng + ix, jp1, k)],
                                               work.f[mass_index(ng + ix, jp2, k)]);
                    }
                }
            }

            std::fill(fx_provisional.begin(), fx_provisional.end(), 0.0);
            std::fill(fu_provisional.begin(), fu_provisional.end(), 0.0);
            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                const int il = ng + iface - 1;
                const int ir = ng + iface;
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                    const double v = bkg_n.cgrid.vx[idx2(j, k)];
                    const double sigma = 0.5 * h * std::fabs(v) / sg.dx;
                    const size_t id = xface_index(iface, j, k);
                    x_left_1d[id] = PpmCtu::trace_right(ppm_x[mass_index(il, j, k)], sigma);
                    x_right_1d[id] = PpmCtu::trace_left(ppm_x[mass_index(ir, j, k)], sigma);
                    fx_provisional[id] = v * PpmCtu::upwind_state(x_left_1d[id],
                                                           x_right_1d[id], v) / sg.dx;
                }
            }
            close_periodic_face_blocks(fx_provisional, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 902);

            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                const double a = bkg_n.charge * fields_mid.Ex[ng + ix] /
                                 (bkg_n.mass * Const::c);
                for (int jf = 1; jf < Param::Nv; ++jf) {
                    const int jl = jf - 1;
                    const int jr = jf;
                    const double du_l = bkg_n.cgrid.upar_widths[jl];
                    const double du_r = bkg_n.cgrid.upar_widths[jr];
                    const size_t id_base = static_cast<size_t>(ix) * Param::Nvmu;
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = uface_index(ix, jf, k);
                        u_left_1d[id] = PpmCtu::trace_right(
                            ppm_u[id_base + idx2(jl, k)], 0.5 * h * std::fabs(a) / du_l);
                        u_right_1d[id] = PpmCtu::trace_left(
                            ppm_u[id_base + idx2(jr, k)], 0.5 * h * std::fabs(a) / du_r);
                        const double state = PpmCtu::upwind_state(u_left_1d[id], u_right_1d[id], a);
                        const double du_face = 0.5 * (du_l + du_r);
                        cu_high[id] = state / du_face;
                        fu_provisional[id] = a * cu_high[id];
                    }
                }
            }

            // D--F: Colella transverse corrections use only provisional
            // normal fluxes.  A second Riemann solve produces the sole PPM-CTU
            // high-order candidate consumed by the shared two-dimensional FCT.
            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                const int ixl = (iface == 0) ? nxl - 1 : iface - 1;
                const int ixr = (iface == nxl) ? 0 : iface;
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = xface_index(iface, j, k);
                    const double v = bkg_n.cgrid.vx[idx2(j, k)];
                    const double left_star = x_left_1d[id] - 0.5 * h *
                        (fu_provisional[uface_index(ixl, j + 1, k)] -
                         fu_provisional[uface_index(ixl, j, k)]);
                    const double right_star = x_right_1d[id] - 0.5 * h *
                        (fu_provisional[uface_index(ixr, j + 1, k)] -
                         fu_provisional[uface_index(ixr, j, k)]);
                    fx_high[id] = v * PpmCtu::upwind_state(left_star, right_star, v) / sg.dx;
                }
            }
            close_periodic_face_blocks(fx_high, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 908);
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                const double a = bkg_n.charge * fields_mid.Ex[ng + ix] /
                                 (bkg_n.mass * Const::c);
                for (int jf = 1; jf < Param::Nv; ++jf) for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = uface_index(ix, jf, k);
                    const double left_star = u_left_1d[id] - 0.5 * h *
                        (fx_provisional[xface_index(ix + 1, jf - 1, k)] -
                         fx_provisional[xface_index(ix, jf - 1, k)]);
                    const double right_star = u_right_1d[id] - 0.5 * h *
                        (fx_provisional[xface_index(ix + 1, jf, k)] -
                         fx_provisional[xface_index(ix, jf, k)]);
                    const double du = 0.5 * (bkg_n.cgrid.upar_widths[jf - 1] +
                                             bkg_n.cgrid.upar_widths[jf]);
                    const double state = PpmCtu::upwind_state(left_star, right_star, a);
                    cu_high[id] = state / du;
                    fu_high[id] = a * cu_high[id];
                }
            }

            Species& low_state = low_state_buffer;
            low_state = work;
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        low_state.f[mass_index(ng + ix, j, k)] = work.f[
                            mass_index(ng + ix, j, k)] - h *
                            (fx_low[xface_index(ix + 1, j, k)] -
                             fx_low[xface_index(ix, j, k)] +
                             fu_low[uface_index(ix, j + 1, k)] -
                             fu_low[uface_index(ix, j, k)]);
                    }
                }
            }
            exchange_ghosts_x_persistent(low_state, sg, mpi_rank, mpi_size);

            double low_min = std::numeric_limits<double>::infinity();
            double low_scale = 0.0;
            double low_tolerance_linf = 0.0;
            double low_positivity_violation = 0.0;
            double low_nonfinite = 0.0;
            #pragma omp parallel for reduction(min:low_min) reduction(max:low_scale,low_tolerance_linf,low_positivity_violation,low_nonfinite)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double m = low_state.f[mass_index(ng + ix, j, k)];
                    const double transfer_scale = h * (
                        std::fabs(fx_low[xface_index(ix, j, k)]) +
                        std::fabs(fx_low[xface_index(ix + 1, j, k)]) +
                        std::fabs(fu_low[uface_index(ix, j, k)]) +
                        std::fabs(fu_low[uface_index(ix, j + 1, k)]));
                    const double tolerance = 128.0 *
                        std::numeric_limits<double>::epsilon() *
                        (std::fabs(work.f[mass_index(ng + ix, j, k)]) +
                         transfer_scale);
                    low_min = std::min(low_min, m);
                    low_scale = std::max(low_scale, std::fabs(m));
                    low_tolerance_linf = std::max(low_tolerance_linf, tolerance);
                    low_positivity_violation = std::max(
                        low_positivity_violation, -m - tolerance);
                    if (!std::isfinite(m)) low_nonfinite = 1.0;
                }
            double low_guard[5] = {low_min, low_scale, low_tolerance_linf,
                                   low_positivity_violation, low_nonfinite};
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[0], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[1], 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            const double mass_roundoff = low_guard[2];
            result.fct_low_order_tolerance_linf = std::max(
                result.fct_low_order_tolerance_linf, mass_roundoff);
            if (low_guard[4] != 0.0 || !std::isfinite(low_guard[0]) ||
                low_guard[3] > 0.0) {
                double work_input_min = std::numeric_limits<double>::infinity();
                #pragma omp parallel for reduction(min:work_input_min)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k)
                            work_input_min = std::min(work_input_min,
                                work.f[mass_index(ng + ix, j, k)]);
                MPI_Allreduce(MPI_IN_PLACE, &work_input_min, 1, MPI_DOUBLE,
                              MPI_MIN, MPI_COMM_WORLD);

                int local_owner = std::numeric_limits<int>::max();
                for (int ix = 0; ix < nxl && local_owner == std::numeric_limits<int>::max(); ++ix)
                    for (int j = 0; j < Param::Nv && local_owner == std::numeric_limits<int>::max(); ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double m =
                                low_state.f[mass_index(ng + ix, j, k)];
                            const bool target = (low_guard[4] != 0.0)
                                ? !std::isfinite(m) : (m == low_guard[0]);
                            if (target) local_owner = mpi_rank;
                        }
                MPI_Allreduce(MPI_IN_PLACE, &local_owner, 1, MPI_INT, MPI_MIN,
                              MPI_COMM_WORLD);

                // Only the rank that owns the failing cell fills this
                // packet; all ranks receive exactly the same failure data.
                double low_failure[19] = {0.0};
                low_failure[18] = work_input_min;
                if (mpi_rank == local_owner) {
                    bool found = false;
                    for (int ix = 0; ix < nxl && !found; ++ix)
                        for (int j = 0; j < Param::Nv && !found; ++j)
                            for (int k = 0; k < Param::Nmu; ++k) {
                                const double m_low =
                                    low_state.f[mass_index(ng + ix, j, k)];
                                const bool target = (low_guard[4] != 0.0)
                                    ? !std::isfinite(m_low)
                                    : (m_low == low_guard[0]);
                                if (!target) continue;
                                const double tx_left = h *
                                    fx_low[xface_index(ix, j, k)];
                                const double tx_right = h *
                                    fx_low[xface_index(ix + 1, j, k)];
                                const double tu_lower = h *
                                    fu_low[uface_index(ix, j, k)];
                                const double tu_upper = h *
                                    fu_low[uface_index(ix, j + 1, k)];
                                const double a = bkg_n.charge *
                                    fields_mid.Ex[ng + ix] /
                                    (bkg_n.mass * Const::c);
                                low_failure[0] = static_cast<double>(
                                    (sg.ix_start + ix) * Param::Nvmu + idx2(j, k));
                                low_failure[1] =
                                    work.f[mass_index(ng + ix, j, k)];
                                low_failure[2] = m_low;
                                low_failure[3] = tx_left;
                                low_failure[4] = tx_right;
                                low_failure[5] = tu_lower;
                                low_failure[6] = tu_upper;
                                low_failure[7] = std::max(0.0, -tx_left);
                                low_failure[8] = std::max(0.0, tx_right);
                                low_failure[9] = std::max(0.0, -tu_lower);
                                low_failure[10] = std::max(0.0, tu_upper);
                                low_failure[11] = std::max(0.0, tx_left);
                                low_failure[12] = std::max(0.0, -tx_right);
                                low_failure[13] = std::max(0.0, tu_lower);
                                low_failure[14] = std::max(0.0, -tu_upper);
                                low_failure[15] = h * std::fabs(
                                    bkg_n.cgrid.vx[idx2(j, k)]) / sg.dx;
                                low_failure[16] = h * std::fabs(a) /
                                    std::max(bkg_n.cgrid.upar_widths[j],
                                             std::numeric_limits<double>::min());
                                low_failure[17] = (mpi_size > 1 &&
                                    (ix == 0 || ix == nxl - 1)) ? 1.0 : 0.0;
                                found = true;
                            }
                }
                MPI_Bcast(low_failure, 19, MPI_DOUBLE, local_owner,
                          MPI_COMM_WORLD);
                result.failed = true;
                result.state_advanced = 0;
                result.failure_reason = (low_guard[4] == 0.0 &&
                                         std::isfinite(low_guard[0])) ? 2 : 4;
                result.failure_iteration = iter;
                result.failure_substep = sub;
                result.failure_global_cfl = global_cfl;
                result.failure_low_min = low_guard[0];
                result.failure_low_m_in = low_failure[1];
                result.failure_low_m_low = low_failure[2];
                result.failure_low_transfer_x_left = low_failure[3];
                result.failure_low_transfer_x_right = low_failure[4];
                result.failure_low_transfer_u_lower = low_failure[5];
                result.failure_low_transfer_u_upper = low_failure[6];
                result.failure_low_out_x_left = low_failure[7];
                result.failure_low_out_x_right = low_failure[8];
                result.failure_low_out_u_lower = low_failure[9];
                result.failure_low_out_u_upper = low_failure[10];
                result.failure_low_in_x_left = low_failure[11];
                result.failure_low_in_x_right = low_failure[12];
                result.failure_low_in_u_lower = low_failure[13];
                result.failure_low_in_u_upper = low_failure[14];
                result.failure_low_cfl_x = low_failure[15];
                result.failure_low_cfl_u = low_failure[16];
                result.failure_low_on_mpi_interface =
                    static_cast<int>(low_failure[17]);
                result.failure_low_work_input_min = low_failure[18];
                set_failure_location(static_cast<int>(low_failure[0]));
                return result;
            }

            std::fill(alpha_x_left.begin(), alpha_x_left.end(), 1.0);
            std::fill(alpha_x_right.begin(), alpha_x_right.end(), 1.0);
            std::fill(alpha_u_lower.begin(), alpha_u_lower.end(), 1.0);
            std::fill(alpha_u_upper.begin(), alpha_u_upper.end(), 1.0);
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                    const double ax_left = fx_high[xface_index(ix, j, k)] -
                                           fx_low[xface_index(ix, j, k)];
                    const double ax_right = fx_high[xface_index(ix + 1, j, k)] -
                                            fx_low[xface_index(ix + 1, j, k)];
                    const double au_lower = fu_high[uface_index(ix, j, k)] -
                                            fu_low[uface_index(ix, j, k)];
                    const double au_upper = fu_high[uface_index(ix, j + 1, k)] -
                                            fu_low[uface_index(ix, j + 1, k)];
                    const double contribution[4] = {
                        h * std::max(0.0, -ax_left),
                        h * std::max(0.0, ax_right),
                        h * std::max(0.0, -au_lower),
                        h * std::max(0.0, au_upper)
                    };
                    double alpha[4];
                    const size_t p = static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k);
                    Stage5::shared_budget_alphas(contribution,
                        low_state.f[mass_index(ng + ix, j, k)], alpha);
                    alpha_x_left[p] = alpha[0];
                    alpha_x_right[p] = alpha[1];
                    alpha_u_lower[p] = alpha[2];
                    alpha_u_upper[p] = alpha[3];
                }
            }
            #pragma omp parallel for reduction(max:fct_budget_violation) schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                    const double ax_left = fx_high[xface_index(ix, j, k)] -
                                           fx_low[xface_index(ix, j, k)];
                    const double ax_right = fx_high[xface_index(ix + 1, j, k)] -
                                            fx_low[xface_index(ix + 1, j, k)];
                    const double au_lower = fu_high[uface_index(ix, j, k)] -
                                            fu_low[uface_index(ix, j, k)];
                    const double au_upper = fu_high[uface_index(ix, j + 1, k)] -
                                            fu_low[uface_index(ix, j + 1, k)];
                    const size_t p = static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k);
                    const double used = h * (
                        std::max(0.0, -ax_left) * alpha_x_left[p] +
                        std::max(0.0, ax_right) * alpha_x_right[p] +
                        std::max(0.0, -au_lower) * alpha_u_lower[p] +
                        std::max(0.0, au_upper) * alpha_u_upper[p]);
                    fct_budget_violation = std::max(fct_budget_violation,
                        used - low_state.f[mass_index(ng + ix, j, k)]);
                }
            }

            // The alpha stored on a periodic shared x face is owned by its
            // donor.  Exchange only the left neighbour's right-face alpha.
            if (mpi_size == 1) {
                std::copy(alpha_x_right.begin() + static_cast<size_t>(nxl - 1) * Param::Nvmu,
                          alpha_x_right.begin() + static_cast<size_t>(nxl) * Param::Nvmu,
                          left_alpha_x_right.begin());
            } else {
                const int left = (mpi_rank + mpi_size - 1) % mpi_size;
                const int right = (mpi_rank + 1) % mpi_size;
                MPI_Sendrecv(alpha_x_right.data() + static_cast<size_t>(nxl - 1) * Param::Nvmu,
                             static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                             right, 904, left_alpha_x_right.data(), static_cast<int>(Param::Nvmu),
                             MPI_DOUBLE, left, 904, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            fx_final = fx_low;
            fu_final = fu_low;
            cu_final = cu_low;
            #pragma omp parallel for schedule(static) reduction(+:limiter_faces,limiter_active) reduction(min:limiter_min)
            for (int iface = 0; iface < nxl; ++iface) {
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = xface_index(iface, j, k);
                    const double anti = fx_high[id] - fx_low[id];
                    const size_t p = idx2(j, k);
                    const double alpha = (anti >= 0.0)
                        ? ((iface == 0) ? left_alpha_x_right[p] :
                           alpha_x_right[static_cast<size_t>(iface - 1) * Param::Nvmu + p])
                        : alpha_x_left[static_cast<size_t>(iface) * Param::Nvmu + p];
                    fx_final[id] += alpha * anti;
                    limiter_faces += 1.0;
                    if (alpha < 1.0 - 1.0e-14) limiter_active += 1.0;
                    limiter_min = std::min(limiter_min, alpha);
                }
            }
            close_periodic_face_blocks(fx_final, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 905);

            // P1 MPI-interface audit: the right ghost face must be the
            // right-neighbour's owned face-zero transfer, for low/high/final
            // fluxes alike.  This only verifies the current protocol.
            double interface_sent[3] = {0.0, 0.0, 0.0};
            double interface_expected[3] = {0.0, 0.0, 0.0};
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const double weight = 1.0 + 1.0e-3 * (idx2(j, k) + 1);
                interface_sent[0] += weight * h * fx_low[xface_index(0, j, k)];
                interface_sent[1] += weight * h * fx_high[xface_index(0, j, k)];
                interface_sent[2] += weight * h * fx_final[xface_index(0, j, k)];
                interface_expected[0] += weight * h * fx_low[xface_index(nxl, j, k)];
                interface_expected[1] += weight * h * fx_high[xface_index(nxl, j, k)];
                interface_expected[2] += weight * h * fx_final[xface_index(nxl, j, k)];
            }
            double interface_received[3] = {0.0, 0.0, 0.0};
            if (mpi_size == 1) {
                std::copy(interface_sent, interface_sent + 3, interface_received);
            } else {
                const int left = (mpi_rank + mpi_size - 1) % mpi_size;
                const int right = (mpi_rank + 1) % mpi_size;
                MPI_Sendrecv(interface_sent, 3, MPI_DOUBLE, left, 906,
                             interface_received, 3, MPI_DOUBLE, right, 906,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            double interface_checksum_linf = 0.0;
            double interface_checksum_violation = 0.0;
            for (int q = 0; q < 3; ++q) {
                const double difference = std::fabs(
                    interface_expected[q] - interface_received[q]);
                const double tolerance = 256.0 *
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::max(std::fabs(interface_expected[q]),
                                            std::fabs(interface_received[q])));
                interface_checksum_linf = std::max(interface_checksum_linf,
                                                    difference);
                interface_checksum_violation = std::max(
                    interface_checksum_violation, difference - tolerance);
            }
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) for (int jf = 1; jf < Param::Nv; ++jf)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = uface_index(ix, jf, k);
                    const double anti = fu_high[id] - fu_low[id];
                    const size_t base = static_cast<size_t>(ix) * Param::Nvmu;
                    const double alpha = (anti >= 0.0)
                        ? alpha_u_upper[base + idx2(jf - 1, k)]
                        : alpha_u_lower[base + idx2(jf, k)];
                    fu_final[id] += alpha * anti;
                    cu_final[id] += alpha * (cu_high[id] - cu_low[id]);
                }

            // P1: construct the final state in scratch from the low-order
            // state plus limited anti-diffusive mass transfers.  The four
            // transfer terms below are h * (Phi^H - Phi^L), hence have the
            // same units as the cell-integrated cylindrical mass.
            double candidate_min = std::numeric_limits<double>::infinity();
            double high_low_identity_linf = 0.0;
            double high_low_identity_violation = 0.0;
            double donor_capacity_violation = 0.0;
            double final_tolerance_linf = 0.0;
            double final_positivity_violation = 0.0;
            double candidate_nonfinite = 0.0;
            double local_identity_ratio = 0.0;
            double local_identity_residual = 0.0;
            double local_identity_scale = 0.0;
            int local_identity_global_linear = -1;
            double local_donor_ratio = 0.0;
            double local_donor_m_low = 0.0;
            double local_donor_transfer[4] = {0.0, 0.0, 0.0, 0.0};
            double local_donor_alpha[4] = {1.0, 1.0, 1.0, 1.0};
            double local_donor_outflow = 0.0;
            double local_donor_scale = 1.0;
            int local_donor_global_linear = -1;
            double donor_roundoff_warning = 0.0;
            double donor_beta_count = 0.0;
            double donor_beta_min = 1.0;
            double roundoff_normalized_count = 0.0;
            double roundoff_normalized_mass = 0.0;
            std::fill(donor_beta.begin(), donor_beta.end(), 1.0);
            #pragma omp parallel
            {
                double thread_identity_ratio = 0.0;
                double thread_identity_residual = 0.0;
                double thread_identity_scale = 0.0;
                int thread_identity_global_linear = -1;
                double thread_donor_ratio = 0.0;
                double thread_donor_m_low = 0.0;
                double thread_donor_transfer[4] = {0.0, 0.0, 0.0, 0.0};
                double thread_donor_alpha[4] = {1.0, 1.0, 1.0, 1.0};
                double thread_donor_outflow = 0.0;
                double thread_donor_scale = 1.0;
                int thread_donor_global_linear = -1;
            #pragma omp for schedule(static) reduction(min:candidate_min) reduction(max:high_low_identity_linf,high_low_identity_violation,final_tolerance_linf,final_positivity_violation,candidate_nonfinite,donor_roundoff_warning) reduction(+:donor_beta_count,roundoff_normalized_count,roundoff_normalized_mass) reduction(min:donor_beta_min)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t p = static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k);
                    const double m_low = low_state.f[mass_index(ng + ix, j, k)];
                    const double ax_left_raw = h * (
                        fx_high[xface_index(ix, j, k)] -
                        fx_low[xface_index(ix, j, k)]);
                    const double ax_right_raw = h * (
                        fx_high[xface_index(ix + 1, j, k)] -
                        fx_low[xface_index(ix + 1, j, k)]);
                    const double au_lower_raw = h * (
                        fu_high[uface_index(ix, j, k)] -
                        fu_low[uface_index(ix, j, k)]);
                    const double au_upper_raw = h * (
                        fu_high[uface_index(ix, j + 1, k)] -
                        fu_low[uface_index(ix, j + 1, k)]);
                    const double high_minus_low = -h * (
                        (fx_high[xface_index(ix + 1, j, k)] -
                         fx_high[xface_index(ix, j, k)] +
                         fu_high[uface_index(ix, j + 1, k)] -
                         fu_high[uface_index(ix, j, k)]) -
                        (fx_low[xface_index(ix + 1, j, k)] -
                         fx_low[xface_index(ix, j, k)] +
                         fu_low[uface_index(ix, j + 1, k)] -
                         fu_low[uface_index(ix, j, k)]));
                    const double transfer_divergence = ax_left_raw - ax_right_raw +
                        au_lower_raw - au_upper_raw;
                    const double identity_residual = high_minus_low -
                        transfer_divergence;
                    const double identity_error = std::fabs(identity_residual);
                    const double identity_scale = h * (
                        std::fabs(fx_high[xface_index(ix, j, k)]) +
                        std::fabs(fx_high[xface_index(ix + 1, j, k)]) +
                        std::fabs(fu_high[uface_index(ix, j, k)]) +
                        std::fabs(fu_high[uface_index(ix, j + 1, k)]) +
                        std::fabs(fx_low[xface_index(ix, j, k)]) +
                        std::fabs(fx_low[xface_index(ix + 1, j, k)]) +
                        std::fabs(fu_low[uface_index(ix, j, k)]) +
                        std::fabs(fu_low[uface_index(ix, j + 1, k)]));
                    const double identity_ratio = identity_error /
                        std::max(1.0, identity_scale);
                    const double identity_tolerance = 4096.0 *
                        std::numeric_limits<double>::epsilon();
                    high_low_identity_linf = std::max(high_low_identity_linf,
                                                      identity_error);
                    high_low_identity_violation = std::max(
                        high_low_identity_violation,
                        identity_ratio - identity_tolerance);
                    if (identity_ratio > thread_identity_ratio) {
                        thread_identity_ratio = identity_ratio;
                        thread_identity_residual = identity_residual;
                        thread_identity_scale = identity_scale;
                        thread_identity_global_linear =
                            (sg.ix_start + ix) * Param::Nvmu + idx2(j, k);
                    }

                    const double ax_left = h * (
                        fx_final[xface_index(ix, j, k)] -
                        fx_low[xface_index(ix, j, k)]);
                    const double ax_right = h * (
                        fx_final[xface_index(ix + 1, j, k)] -
                        fx_low[xface_index(ix + 1, j, k)]);
                    const double au_lower = h * (
                        fu_final[uface_index(ix, j, k)] -
                        fu_low[uface_index(ix, j, k)]);
                    const double au_upper = h * (
                        fu_final[uface_index(ix, j + 1, k)] -
                        fu_low[uface_index(ix, j + 1, k)]);
                    const long double candidate_ld =
                        static_cast<long double>(m_low) +
                        static_cast<long double>(ax_left) -
                        static_cast<long double>(ax_right) +
                        static_cast<long double>(au_lower) -
                        static_cast<long double>(au_upper);
                    double candidate =
                        static_cast<double>(candidate_ld);
                    const double final_tolerance = candidate_roundoff_tolerance(
                        m_low, ax_left, ax_right, au_lower, au_upper);
                    double normalized_mass = 0.0;
                    if (normalize_roundoff_negative_candidate(
                            candidate, final_tolerance, normalized_mass)) {
                        roundoff_normalized_count += 1.0;
                        roundoff_normalized_mass += normalized_mass;
                    }
                    candidate_mass[p] = candidate;
                    // Final donor closure uses only the low-order donor
                    // capacity.  It is computed for every cell, not only
                    // after an over-budget diagnostic.
                    const long double limited_outflow_ld =
                        std::max(0.0L, -static_cast<long double>(ax_left)) +
                        std::max(0.0L, static_cast<long double>(ax_right)) +
                        std::max(0.0L, -static_cast<long double>(au_lower)) +
                        std::max(0.0L, static_cast<long double>(au_upper));
                    const long double donor_scale_ld = std::max(1.0L,
                        std::max(static_cast<long double>(m_low),
                            std::fabs(static_cast<long double>(ax_left)) +
                            std::fabs(static_cast<long double>(ax_right)) +
                            std::fabs(static_cast<long double>(au_lower)) +
                            std::fabs(static_cast<long double>(au_upper))));
                    const long double q_safe_ld = std::max(0.0L,
                        static_cast<long double>(m_low)) *
                        (1.0L - 32.0L *
                         std::numeric_limits<double>::epsilon());
                    const long double donor_excess_ld =
                        limited_outflow_ld - q_safe_ld;
                    const double donor_relative =
                        (donor_excess_ld > 0.0L)
                        ? static_cast<double>(donor_excess_ld / donor_scale_ld)
                        : 0.0;
                    const double donor_roundoff_tolerance = 4096.0 *
                        std::numeric_limits<double>::epsilon();
                    const double face_transfer[4] = {
                        ax_left, ax_right, au_lower, au_upper};
                    const double face_raw_transfer[4] = {
                        ax_left_raw, ax_right_raw, au_lower_raw, au_upper_raw};
                    double face_alpha[4];
                    for (int f = 0; f < 4; ++f) {
                        face_alpha[f] = (face_raw_transfer[f] != 0.0)
                            ? face_transfer[f] / face_raw_transfer[f] : 1.0;
                    }
                    if (donor_relative > thread_donor_ratio) {
                        thread_donor_ratio = donor_relative;
                        thread_donor_m_low = m_low;
                        for (int f = 0; f < 4; ++f) {
                            thread_donor_transfer[f] = face_transfer[f];
                            thread_donor_alpha[f] = face_alpha[f];
                        }
                        thread_donor_outflow =
                            static_cast<double>(limited_outflow_ld);
                        thread_donor_scale = static_cast<double>(donor_scale_ld);
                        thread_donor_global_linear =
                            (sg.ix_start + ix) * Param::Nvmu + idx2(j, k);
                    }
                    if (donor_excess_ld > 0.0L &&
                        donor_relative <= donor_roundoff_tolerance) {
                        donor_roundoff_warning = 1.0;
                    }
                    const double beta = (limited_outflow_ld > 0.0L)
                        ? std::max(0.0, std::min(1.0,
                            static_cast<double>(q_safe_ld /
                                                limited_outflow_ld)))
                        : 1.0;
                    donor_beta[p] = beta;
                    if (beta < 1.0) {
                        donor_beta_count += 1.0;
                        donor_beta_min = std::min(donor_beta_min, beta);
                    }
                    final_tolerance_linf = std::max(final_tolerance_linf,
                                                    final_tolerance);
                    if (std::isfinite(candidate)) {
                        candidate_min = std::min(candidate_min, candidate);
                        final_positivity_violation = std::max(
                            final_positivity_violation,
                            -candidate);
                    } else {
                        candidate_nonfinite = 1.0;
                    }
                }
                #pragma omp critical(cylindrical_identity_worst)
                {
                    if (thread_identity_ratio > local_identity_ratio) {
                        local_identity_ratio = thread_identity_ratio;
                        local_identity_residual = thread_identity_residual;
                        local_identity_scale = thread_identity_scale;
                        local_identity_global_linear =
                            thread_identity_global_linear;
                    }
                }
                #pragma omp critical(cylindrical_donor_worst)
                {
                    if (thread_donor_ratio > local_donor_ratio) {
                        local_donor_ratio = thread_donor_ratio;
                        local_donor_m_low = thread_donor_m_low;
                        for (int f = 0; f < 4; ++f) {
                            local_donor_transfer[f] = thread_donor_transfer[f];
                            local_donor_alpha[f] = thread_donor_alpha[f];
                        }
                        local_donor_outflow = thread_donor_outflow;
                        local_donor_scale = thread_donor_scale;
                        local_donor_global_linear =
                            thread_donor_global_linear;
                    }
                }
            }

            struct IdentityMaxLoc {
                double ratio;
                int rank;
            };
            IdentityMaxLoc local_identity_owner = {
                local_identity_ratio, mpi_rank};
            IdentityMaxLoc global_identity_owner = {0.0, 0};
            MPI_Allreduce(&local_identity_owner, &global_identity_owner, 1,
                          MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
            if (global_identity_owner.ratio >
                result.fct_high_low_identity_ratio_linf) {
                double identity_metadata[3] = {0.0, 0.0, -1.0};
                if (mpi_rank == global_identity_owner.rank) {
                    identity_metadata[0] = local_identity_residual;
                    identity_metadata[1] = local_identity_scale;
                    identity_metadata[2] =
                        static_cast<double>(local_identity_global_linear);
                }
                MPI_Bcast(identity_metadata, 3, MPI_DOUBLE,
                          global_identity_owner.rank, MPI_COMM_WORLD);
                result.fct_high_low_identity_ratio_linf =
                    global_identity_owner.ratio;
                result.fct_high_low_identity_worst_residual =
                    identity_metadata[0];
                result.fct_high_low_identity_worst_scale =
                    identity_metadata[1];
                result.fct_high_low_identity_worst_relative =
                    (identity_metadata[1] != 0.0)
                    ? identity_metadata[0] / identity_metadata[1] : 0.0;
                const int global_linear =
                    static_cast<int>(identity_metadata[2]);
                if (global_linear >= 0) {
                    result.fct_high_low_identity_worst_ix =
                        global_linear / Param::Nvmu;
                    const int velocity_linear =
                        global_linear % Param::Nvmu;
                    result.fct_high_low_identity_worst_iv =
                        velocity_linear / Param::Nmu;
                    result.fct_high_low_identity_worst_imu =
                        velocity_linear % Param::Nmu;
                }
            }

            struct DonorMaxLoc {
                double ratio;
                int rank;
            };
            DonorMaxLoc local_donor_owner = {local_donor_ratio, mpi_rank};
            DonorMaxLoc global_donor_owner = {0.0, 0};
            MPI_Allreduce(&local_donor_owner, &global_donor_owner, 1,
                          MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
            if (global_donor_owner.ratio >
                result.fct_donor_worst_relative) {
                double donor_metadata[12] = {
                    0.0, 0.0, 0.0, 0.0, 0.0,
                    1.0, 1.0, 1.0, 1.0,
                    0.0, 1.0, -1.0};
                if (mpi_rank == global_donor_owner.rank) {
                    donor_metadata[0] = local_donor_m_low;
                    for (int f = 0; f < 4; ++f) {
                        donor_metadata[1 + f] = local_donor_transfer[f];
                        donor_metadata[5 + f] = local_donor_alpha[f];
                    }
                    donor_metadata[9] = local_donor_outflow;
                    donor_metadata[10] = local_donor_scale;
                    donor_metadata[11] =
                        static_cast<double>(local_donor_global_linear);
                }
                MPI_Bcast(donor_metadata, 12, MPI_DOUBLE,
                          global_donor_owner.rank, MPI_COMM_WORLD);
                result.fct_donor_worst_relative = global_donor_owner.ratio;
                result.fct_donor_worst_m_low = donor_metadata[0];
                result.fct_donor_worst_ax_left = donor_metadata[1];
                result.fct_donor_worst_ax_right = donor_metadata[2];
                result.fct_donor_worst_au_lower = donor_metadata[3];
                result.fct_donor_worst_au_upper = donor_metadata[4];
                result.fct_donor_worst_alpha_x_left = donor_metadata[5];
                result.fct_donor_worst_alpha_x_right = donor_metadata[6];
                result.fct_donor_worst_alpha_u_lower = donor_metadata[7];
                result.fct_donor_worst_alpha_u_upper = donor_metadata[8];
                result.fct_donor_worst_outflow = donor_metadata[9];
                result.fct_donor_worst_scale = donor_metadata[10];
                const int global_linear =
                    static_cast<int>(donor_metadata[11]);
                if (global_linear >= 0) {
                    result.fct_donor_worst_ix =
                        global_linear / Param::Nvmu;
                    const int velocity_linear =
                        global_linear % Param::Nvmu;
                    result.fct_donor_worst_iv =
                        velocity_linear / Param::Nmu;
                    result.fct_donor_worst_imu =
                        velocity_linear % Param::Nmu;
                }
            }

            double donor_control[3] = {
                donor_beta_count > 0.0 ? 1.0 : 0.0,
                donor_beta_count > 0.0 ? 1.0 - donor_beta_min : 0.0,
                donor_roundoff_warning
            };
            MPI_Allreduce(MPI_IN_PLACE, donor_control, 3, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            result.fct_donor_roundoff_warning = std::max(
                result.fct_donor_roundoff_warning,
                static_cast<int>(donor_control[2]));

            // A real donor excess is repaired locally by scaling only the
            // anti-diffusive outflow faces of the donor.  The low-order
            // state and every incoming correction remain untouched.
            if (donor_control[0] != 0.0) {
                double global_beta_count = donor_beta_count;
                MPI_Allreduce(MPI_IN_PLACE, &global_beta_count, 1, MPI_DOUBLE,
                              MPI_SUM, MPI_COMM_WORLD);
                result.fct_donor_beta_applied_count +=
                    static_cast<long long>(global_beta_count);
                result.fct_donor_beta_min = std::min(
                    result.fct_donor_beta_min, 1.0 - donor_control[1]);

                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t p = static_cast<size_t>(ix) *
                                Param::Nvmu + idx2(j, k);
                            const double beta = donor_beta[p];
                            if (beta >= 1.0) continue;
                            const double ax_left_raw = h * (
                                fx_high[xface_index(ix, j, k)] -
                                fx_low[xface_index(ix, j, k)]);
                            const double ax_right_raw = h * (
                                fx_high[xface_index(ix + 1, j, k)] -
                                fx_low[xface_index(ix + 1, j, k)]);
                            const double au_lower_raw = h * (
                                fu_high[uface_index(ix, j, k)] -
                                fu_low[uface_index(ix, j, k)]);
                            const double au_upper_raw = h * (
                                fu_high[uface_index(ix, j + 1, k)] -
                                fu_low[uface_index(ix, j + 1, k)]);
                            if (ax_left_raw < 0.0) alpha_x_left[p] *= beta;
                            if (ax_right_raw > 0.0) alpha_x_right[p] *= beta;
                            if (au_lower_raw < 0.0) alpha_u_lower[p] *= beta;
                            if (au_upper_raw > 0.0) alpha_u_upper[p] *= beta;
                        }

                if (mpi_size == 1) {
                    std::copy(alpha_x_right.begin() +
                                  static_cast<size_t>(nxl - 1) * Param::Nvmu,
                              alpha_x_right.begin() +
                                  static_cast<size_t>(nxl) * Param::Nvmu,
                              left_alpha_x_right.begin());
                } else {
                    const int left = (mpi_rank + mpi_size - 1) % mpi_size;
                    const int right = (mpi_rank + 1) % mpi_size;
                    MPI_Sendrecv(alpha_x_right.data() +
                                     static_cast<size_t>(nxl - 1) * Param::Nvmu,
                                 static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                                 right, 907, left_alpha_x_right.data(),
                                 static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                                 left, 907, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                fx_final = fx_low;
                #pragma omp parallel for schedule(static)
                for (int iface = 0; iface < nxl; ++iface)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t id = xface_index(iface, j, k);
                            const double anti = fx_high[id] - fx_low[id];
                            const size_t p = idx2(j, k);
                            const double alpha = (anti >= 0.0)
                                ? ((iface == 0) ? left_alpha_x_right[p] :
                                   alpha_x_right[static_cast<size_t>(iface - 1) *
                                                 Param::Nvmu + p])
                                : alpha_x_left[static_cast<size_t>(iface) *
                                               Param::Nvmu + p];
                            fx_final[id] += alpha * anti;
                        }
                close_periodic_face_blocks(fx_final, nxl, Param::Nvmu,
                                           mpi_rank, mpi_size, 905);

                fu_final = fu_low;
                cu_final = cu_low;
                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int jf = 1; jf < Param::Nv; ++jf)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t id = uface_index(ix, jf, k);
                            const double anti = fu_high[id] - fu_low[id];
                            const size_t base =
                                static_cast<size_t>(ix) * Param::Nvmu;
                            const double alpha = (anti >= 0.0)
                                ? alpha_u_upper[base + idx2(jf - 1, k)]
                                : alpha_u_lower[base + idx2(jf, k)];
                            fu_final[id] += alpha * anti;
                            cu_final[id] += alpha * (cu_high[id] - cu_low[id]);
                        }

                // Recheck the synchronized final transfers and the scratch
                // state after beta.  This is the only path that can retain a
                // donor-capacity hard failure.
                candidate_min = std::numeric_limits<double>::infinity();
                donor_capacity_violation = 0.0;
                final_tolerance_linf = 0.0;
                final_positivity_violation = 0.0;
                candidate_nonfinite = 0.0;
                roundoff_normalized_count = 0.0;
                roundoff_normalized_mass = 0.0;
                #pragma omp parallel for schedule(static) reduction(min:candidate_min) reduction(max:donor_capacity_violation,final_tolerance_linf,final_positivity_violation,candidate_nonfinite) reduction(+:roundoff_normalized_count,roundoff_normalized_mass)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t p = static_cast<size_t>(ix) *
                                Param::Nvmu + idx2(j, k);
                            const double m_low =
                                low_state.f[mass_index(ng + ix, j, k)];
                            const double ax_left = h * (
                                fx_final[xface_index(ix, j, k)] -
                                fx_low[xface_index(ix, j, k)]);
                            const double ax_right = h * (
                                fx_final[xface_index(ix + 1, j, k)] -
                                fx_low[xface_index(ix + 1, j, k)]);
                            const double au_lower = h * (
                                fu_final[uface_index(ix, j, k)] -
                                fu_low[uface_index(ix, j, k)]);
                            const double au_upper = h * (
                                fu_final[uface_index(ix, j + 1, k)] -
                                fu_low[uface_index(ix, j + 1, k)]);
                            const long double candidate_ld =
                                static_cast<long double>(m_low) +
                                static_cast<long double>(ax_left) -
                                static_cast<long double>(ax_right) +
                                static_cast<long double>(au_lower) -
                                static_cast<long double>(au_upper);
                            double candidate =
                                static_cast<double>(candidate_ld);
                            const long double outflow =
                                std::max(0.0L, -static_cast<long double>(ax_left)) +
                                std::max(0.0L, static_cast<long double>(ax_right)) +
                                std::max(0.0L, -static_cast<long double>(au_lower)) +
                                std::max(0.0L, static_cast<long double>(au_upper));
                            const long double scale = std::max(1.0L,
                                std::max(static_cast<long double>(m_low),
                                    std::fabs(static_cast<long double>(ax_left)) +
                                    std::fabs(static_cast<long double>(ax_right)) +
                                    std::fabs(static_cast<long double>(au_lower)) +
                                    std::fabs(static_cast<long double>(au_upper))));
                            const long double q_safe = std::max(0.0L,
                                static_cast<long double>(m_low)) *
                                (1.0L - 32.0L *
                                 std::numeric_limits<double>::epsilon());
                            const long double excess = outflow - q_safe;
                            const double relative = (excess > 0.0L)
                                ? static_cast<double>(excess / scale) : 0.0;
                            donor_capacity_violation = std::max(
                                donor_capacity_violation, relative -
                                4096.0 * std::numeric_limits<double>::epsilon());
                            const double tolerance =
                                candidate_roundoff_tolerance(
                                    m_low, ax_left, ax_right, au_lower,
                                    au_upper);
                            double normalized_mass = 0.0;
                            if (normalize_roundoff_negative_candidate(
                                    candidate, tolerance, normalized_mass)) {
                                roundoff_normalized_count += 1.0;
                                roundoff_normalized_mass += normalized_mass;
                            }
                            candidate_mass[p] = candidate;
                            final_tolerance_linf = std::max(
                                final_tolerance_linf, tolerance);
                            if (std::isfinite(candidate)) {
                                candidate_min = std::min(candidate_min,
                                                         candidate);
                                final_positivity_violation = std::max(
                                    final_positivity_violation,
                                    -candidate);
                            } else {
                                candidate_nonfinite = 1.0;
                            }
                        }

                double post_audit[4] = {
                    donor_capacity_violation, final_tolerance_linf,
                    final_positivity_violation, candidate_nonfinite};
                MPI_Allreduce(MPI_IN_PLACE, post_audit, 4, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &candidate_min, 1, MPI_DOUBLE,
                              MPI_MIN, MPI_COMM_WORLD);
                donor_capacity_violation = post_audit[0];
                final_tolerance_linf = post_audit[1];
                final_positivity_violation = post_audit[2];
                candidate_nonfinite = post_audit[3];

                double recheck_sent[3] = {0.0, 0.0, 0.0};
                double recheck_expected[3] = {0.0, 0.0, 0.0};
                for (int j = 0; j < Param::Nv; ++j)
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const double weight =
                            1.0 + 1.0e-3 * (idx2(j, k) + 1);
                        recheck_sent[0] += weight * h *
                            fx_low[xface_index(0, j, k)];
                        recheck_sent[1] += weight * h *
                            fx_high[xface_index(0, j, k)];
                        recheck_sent[2] += weight * h *
                            fx_final[xface_index(0, j, k)];
                        recheck_expected[0] += weight * h *
                            fx_low[xface_index(nxl, j, k)];
                        recheck_expected[1] += weight * h *
                            fx_high[xface_index(nxl, j, k)];
                        recheck_expected[2] += weight * h *
                            fx_final[xface_index(nxl, j, k)];
                    }
                double recheck_received[3] = {0.0, 0.0, 0.0};
                if (mpi_size == 1) {
                    std::copy(recheck_sent, recheck_sent + 3,
                              recheck_received);
                } else {
                    const int left = (mpi_rank + mpi_size - 1) % mpi_size;
                    const int right = (mpi_rank + 1) % mpi_size;
                    MPI_Sendrecv(recheck_sent, 3, MPI_DOUBLE, left, 909,
                                 recheck_received, 3, MPI_DOUBLE, right, 909,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                interface_checksum_linf = 0.0;
                interface_checksum_violation = 0.0;
                for (int q = 0; q < 3; ++q) {
                    const double difference = std::fabs(
                        recheck_expected[q] - recheck_received[q]);
                    const double tolerance = 256.0 *
                        std::numeric_limits<double>::epsilon() * std::max(
                            1.0, std::max(std::fabs(recheck_expected[q]),
                                          std::fabs(recheck_received[q])));
                    interface_checksum_linf = std::max(
                        interface_checksum_linf, difference);
                    interface_checksum_violation = std::max(
                        interface_checksum_violation, difference - tolerance);
                }
            }

            // The count and mass are accumulated from the final reconstructed
            // scratch state.  When donor beta was applied, the recheck above
            // reset these local accumulators, so pre-beta roundoff events are
            // not reported as committed normalizations.
            double roundoff_normalization[2] = {
                roundoff_normalized_count, roundoff_normalized_mass};
            MPI_Allreduce(MPI_IN_PLACE, roundoff_normalization, 2, MPI_DOUBLE,
                          MPI_SUM, MPI_COMM_WORLD);
            result.fct_roundoff_zeroed_count +=
                static_cast<long long>(roundoff_normalization[0]);
            result.fct_roundoff_zeroed_mass += roundoff_normalization[1];

            double invariant_max[8] = {
                high_low_identity_linf, high_low_identity_violation,
                donor_capacity_violation, interface_checksum_linf,
                interface_checksum_violation, final_tolerance_linf,
                final_positivity_violation, candidate_nonfinite
            };
            MPI_Allreduce(MPI_IN_PLACE, invariant_max, 8, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &candidate_min, 1, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            result.fct_high_low_identity_linf = std::max(
                result.fct_high_low_identity_linf, invariant_max[0]);
            result.fct_high_low_identity_violation = std::max(
                result.fct_high_low_identity_violation, invariant_max[1]);
            result.fct_donor_capacity_violation = std::max(
                result.fct_donor_capacity_violation, invariant_max[2]);
            result.fct_interface_checksum_linf = std::max(
                result.fct_interface_checksum_linf, invariant_max[3]);
            result.fct_interface_checksum_violation = std::max(
                result.fct_interface_checksum_violation, invariant_max[4]);
            result.fct_final_tolerance_linf = std::max(
                result.fct_final_tolerance_linf, invariant_max[5]);
            result.fct_final_scratch_min = std::min(
                result.fct_final_scratch_min, candidate_min);

            int invariant_failure_reason = 0;
            if (invariant_max[7] != 0.0) invariant_failure_reason = 4;
            else if (invariant_max[1] > 0.0) invariant_failure_reason = 5;
            else if (invariant_max[2] > 0.0) invariant_failure_reason = 6;
            else if (invariant_max[4] > 0.0) invariant_failure_reason = 7;
            else if (invariant_max[6] > 0.0) invariant_failure_reason = 3;
            if (invariant_failure_reason != 0) {
                int local_worst = std::numeric_limits<int>::max();
                if (candidate_min == std::numeric_limits<double>::infinity()) {
                    local_worst = std::numeric_limits<int>::max();
                } else {
                    for (int ix = 0; ix < nxl && local_worst == std::numeric_limits<int>::max(); ++ix)
                        for (int j = 0; j < Param::Nv && local_worst == std::numeric_limits<int>::max(); ++j)
                            for (int k = 0; k < Param::Nmu; ++k)
                                if (candidate_mass[static_cast<size_t>(ix) * Param::Nvmu +
                                                   idx2(j, k)] == candidate_min) {
                                    local_worst = (sg.ix_start + ix) * Param::Nvmu + idx2(j, k);
                                    break;
                                }
                }
                int global_worst = local_worst;
                MPI_Allreduce(MPI_IN_PLACE, &global_worst, 1, MPI_INT, MPI_MIN,
                              MPI_COMM_WORLD);
                result.failed = true;
                result.state_advanced = 0;
                result.failure_reason = invariant_failure_reason;
                result.failure_iteration = iter;
                result.failure_substep = sub;
                result.failure_global_cfl = global_cfl;
                result.failure_low_min = low_guard[0];
                result.failure_final_min = candidate_min;
                if (global_worst != std::numeric_limits<int>::max())
                    set_failure_location(global_worst);
                return result;
            }

            // Account for the final, possibly beta-corrected u transfer only
            // after it has passed the scratch audit.
            #pragma omp parallel for schedule(static) reduction(+:limiter_energy_change,limiter_energy_positive,limiter_energy_negative,limiter_energy_core,limiter_energy_boundary)
            for (int ix = 0; ix < nxl; ++ix)
                for (int jf = 1; jf < Param::Nv; ++jf)
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = uface_index(ix, jf, k);
                        const double dlim = h * (
                            bkg_n.cgrid.kinetic_energy[idx2(jf, k)] -
                            bkg_n.cgrid.kinetic_energy[idx2(jf - 1, k)]) *
                            (fu_final[id] - fu_high[id]);
                        limiter_energy_change += dlim;
                        if (dlim >= 0.0) limiter_energy_positive += dlim;
                        else limiter_energy_negative += dlim;
                        const double x = sg.x(ng + ix);
                        if (x >= 0.1 * Const::micro &&
                            x <= Param::Lx - 0.1 * Const::micro)
                            limiter_energy_core += dlim;
                        else
                            limiter_energy_boundary += dlim;
                    }

            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k)
                    work.f[mass_index(ng + ix, j, k)] =
                        candidate_mass[static_cast<size_t>(ix) * Param::Nvmu +
                                       idx2(j, k)];
            if (!transport_safe_step_gate(work, sub + 1))
                return result;

            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                double gamma = 0.0;
                double psi_k = 0.0;
                double psi_p = 0.0;
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k)
                {
                    const double flux = fx_final[xface_index(iface, j, k)];
                    gamma += flux;
                    psi_k += bkg_n.cgrid.kinetic_energy[idx2(j, k)] * flux;
                    psi_p += bkg_n.mass * Const::c * bkg_n.cgrid.upar_cells[j] * flux;
                }
                integrated_jn[iface] += h * bkg_n.charge * gamma;
                psi_k_x[iface] += h * psi_k;
                psi_p_x[iface] += h * psi_p;
            }
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                double je = 0.0;
                double u_energy = 0.0;
                double u_momentum = 0.0;
                for (int jf = 1; jf < Param::Nv; ++jf) for (int k = 0; k < Param::Nmu; ++k) {
                    const double dke = Stage5::delta_energy(bkg_n.cgrid, jf, k);
                    const double flux = fu_final[uface_index(ix, jf, k)];
                    je += dke * (bkg_n.charge / (bkg_n.mass * Const::c)) *
                          cu_final[uface_index(ix, jf, k)] / sg.dx;
                    u_energy += dke * flux;
                    u_momentum += bkg_n.mass * Const::c *
                        Stage5::u_face_width(bkg_n.cgrid, jf) * flux;
                }
                integrated_je_cell[ix] += h * je;
                u_energy_cell[ix] += h * u_energy;
                u_momentum_cell[ix] += h * u_momentum;
            }
            for (int k = 0; k < Param::Nmu; ++k) {
                const double ke_lo = bkg_n.cgrid.kinetic_energy[idx2(0, k)];
                const double ke_hi = bkg_n.cgrid.kinetic_energy[idx2(Param::Nv - 1, k)];
                for (int ix = 0; ix < nxl; ++ix) {
                    const double f_hi = fu_final[uface_index(ix, Param::Nv, k)];
                    const double f_lo = fu_final[uface_index(ix, 0, k)];
                    u_boundary_energy += h * (ke_hi * f_hi - ke_lo * f_lo);
                    u_boundary_particle += h * (f_hi - f_lo);
                    u_boundary_momentum += h * bkg_n.mass * Const::c *
                        (bkg_n.cgrid.upar_cells[Param::Nv - 1] * f_hi -
                         bkg_n.cgrid.upar_cells[0] * f_lo);
                    u_boundary_energy_cell[ix] += -h * (ke_hi * f_hi - ke_lo * f_lo);
                    u_boundary_momentum_cell[ix] += -h * bkg_n.mass * Const::c *
                        (bkg_n.cgrid.upar_cells[Param::Nv - 1] * f_hi -
                         bkg_n.cgrid.upar_cells[0] * f_lo);
                }
            }
        }

        close_periodic_face_blocks(psi_k_x, nxl, 1, mpi_rank, mpi_size, 909);
        close_periodic_face_blocks(psi_p_x, nxl, 1, mpi_rank, mpi_size, 910);
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            double ke = 0.0;
            double pp = 0.0;
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const double mass = work.f[mass_index(ng + ix, j, k)];
                ke += bkg_n.cgrid.kinetic_energy[idx2(j, k)] * mass;
                pp += bkg_n.mass * Const::c * bkg_n.cgrid.upar_cells[j] * mass;
            }
            final_ke_cell[ix] = ke;
            final_p_cell[ix] = pp;
            delta_ke_cell[ix] = ke - initial_ke_cell[ix];
            delta_p_cell[ix] = pp - initial_p_cell[ix];
            local_energy_residual[ix] = delta_ke_cell[ix] +
                psi_k_x[ix + 1] - psi_k_x[ix] - u_energy_cell[ix] -
                u_boundary_energy_cell[ix];
            local_momentum_residual[ix] = delta_p_cell[ix] +
                psi_p_x[ix + 1] - psi_p_x[ix] - u_momentum_cell[ix] -
                u_boundary_momentum_cell[ix];
        }

        for (size_t i = 0; i < nface; ++i) result.j_bkg_face_mid[i] = integrated_jn[i] / dt;
        std::vector<double>& je_cell = je_cell_buffer;
        for (int ix = 0; ix < nxl; ++ix) je_cell[ix] = integrated_je_cell[ix] / dt;
        double je_left_neighbor = 0.0;
        double je_right_neighbor = 0.0;
        exchange_neighbor_scalars(je_cell, mpi_rank, mpi_size,
                                  je_left_neighbor, je_right_neighbor, 906);
        (void)je_right_neighbor;
        for (int iface = 0; iface < nxl; ++iface) {
            const double je_left = (iface == 0) ? je_left_neighbor : je_cell[iface - 1];
            result.j_bkg_energy_debug_face[iface] = 0.5 * (je_left + je_cell[iface]);
            result.j_beam_face_mid[iface] = jbeam[iface];
            result.j_total_face_mid[iface] = result.j_bkg_face_mid[iface] + jbeam[iface];
        }
        if (nxl > 0) {
            result.j_beam_face_mid[nxl] = jbeam[nxl];
            result.j_total_face_mid[nxl] = result.j_total_face_mid[0];
        }
        close_periodic_face_blocks(result.j_bkg_energy_debug_face, nxl, 1,
                                   mpi_rank, mpi_size, 907);

        EMFields fields_new = fields_n;
        fields_new.advance_ampere_face_from_midpoint_current(result.j_total_face_mid,
                                                              dt, mpi_rank, mpi_size);
        last_fields = fields_new;
        std::vector<double>& next_e = next_e_buffer;
        double local_de = 0.0, local_e = 0.0, local_dj = 0.0, local_j = 0.0;
        for (int iface = 0; iface < nxl; ++iface) {
            next_e[iface] = fields_new.Ex_face[iface];
            local_de = std::max(local_de, std::fabs(next_e[iface] - e_end[iface]));
            local_e = std::max(local_e, std::fabs(next_e[iface]));
            local_j = std::max(local_j, std::fabs(result.j_bkg_face_mid[iface]));
            if (have_previous) local_dj = std::max(local_dj,
                std::fabs(result.j_bkg_face_mid[iface] - previous_j[iface]));
        }
        double norms[4] = {local_de, local_e, local_dj, local_j};
        MPI_Allreduce(MPI_IN_PLACE, norms, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        result.residual_E = norms[0] / std::max(1.0, norms[1]);
        result.residual_J_bkg = norms[2] / std::max(1.0, std::max(std::fabs(Param::jb), norms[3]));

        double local_df = 0.0, local_f = 0.0;
        for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t p = mass_index(ng + ix, j, k);
                local_df = std::max(local_df, std::fabs(work.f[p] - guess.f[p]));
                local_f = std::max(local_f, std::fabs(work.f[p]));
            }
        double f_norms[2] = {local_df, local_f};
        MPI_Allreduce(MPI_IN_PLACE, f_norms, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        result.residual_f = f_norms[0] / std::max(1.0, f_norms[1]);
        result.nonlinear_residual = std::max(result.residual_E / field_tol,
            result.residual_J_bkg / current_tol);

        work.compute_moments();
        work.current_face_x = result.j_bkg_face_mid;
        last_work = work;
        double local_cont = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            double mn = 0.0;
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k)
                mn += bkg_n.f[mass_index(ng + ix, j, k)];
            const double m1 = work.number_density[ix] * sg.dx;
            const double divj = (result.j_bkg_face_mid[ix + 1] -
                                 result.j_bkg_face_mid[ix]) / bkg_n.charge;
            local_cont = std::max(local_cont, std::fabs((m1 - mn) / dt + divj));
        }
        MPI_Allreduce(MPI_IN_PLACE, &local_cont, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        result.continuity_residual_bkg = local_cont;

        double local_cell_work = 0.0, local_face_work = 0.0;
        for (int ix = 0; ix < nxl; ++ix)
            local_cell_work += fields_mid.Ex[ng + ix] * je_cell[ix] * sg.dx;
        for (int iface = 0; iface < nxl; ++iface)
            local_face_work += e_mid[iface] * result.j_bkg_face_mid[iface] * sg.dx;
        double energy_pair[2] = {local_cell_work, local_face_work};
        MPI_Allreduce(MPI_IN_PLACE, energy_pair, 2, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.energy_residual_bkg = dt * (energy_pair[0] - energy_pair[1]);

        double local_stage5[6] = {0.0, 0.0, 0.0, 0.0,
                                  psi_k_x[nxl] - psi_k_x[0],
                                  psi_p_x[nxl] - psi_p_x[0]};
        for (int ix = 0; ix < nxl; ++ix) {
            local_stage5[0] += local_energy_residual[ix];
            local_stage5[1] += local_momentum_residual[ix];
            local_stage5[2] += u_energy_cell[ix];
            local_stage5[3] += u_momentum_cell[ix];
        }
        MPI_Allreduce(MPI_IN_PLACE, local_stage5, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        double u_boundary_global[3] = {u_boundary_particle, u_boundary_momentum,
                                       u_boundary_energy};
        MPI_Allreduce(MPI_IN_PLACE, u_boundary_global, 3, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        double global_fct_budget_violation = fct_budget_violation;
        MPI_Allreduce(MPI_IN_PLACE, &global_fct_budget_violation, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);

        // R_FV audits only the conservative finite-volume update and velocity
        // boundary bookkeeping.  R_couple independently audits the adjoint
        // pairing between the face Ampere current and the cell force current.
        result.stage5_r_fv = local_stage5[0];
        result.stage5_r_couple = dt * (energy_pair[1] - energy_pair[0]);
        // The documented stage-5 aggregate is the sum of the two independently
        // reported residuals.  Keep the algebraic dK-face-work-Bu balance as a
        // separate audit quantity because its sign convention differs.
        result.stage5_r_total = result.stage5_r_fv + result.stage5_r_couple;
        result.stage5_r_physical_balance = result.stage5_r_fv - result.stage5_r_couple;
        result.stage5_u_energy_moment = local_stage5[2];
        result.stage5_u_momentum_moment = local_stage5[3];
        result.stage5_spatial_energy_boundary = local_stage5[4];
        result.stage5_spatial_momentum_boundary = local_stage5[5];
        result.u_boundary_particle_outflow = u_boundary_global[0];
        result.u_boundary_momentum_outflow = u_boundary_global[1];
        result.u_boundary_energy_outflow = u_boundary_global[2];
        result.stage5_fct_budget_violation = global_fct_budget_violation;
        double local_jmax[4] = {0.0, 0.0, 0.0, 0.0};
        for (int iface = 0; iface < nxl; ++iface) {
            const double jn = result.j_bkg_face_mid[iface];
            const double je = result.j_bkg_energy_debug_face[iface];
            local_jmax[0] = std::max(local_jmax[0], std::fabs(jn));
            local_jmax[1] = std::max(local_jmax[1], std::fabs(je));
            local_jmax[2] = std::max(local_jmax[2], std::fabs(jn - je));
            local_jmax[3] += (jn - je) * (jn - je) * sg.dx;
        }
        double j_l2_local = local_jmax[3];
        MPI_Allreduce(MPI_IN_PLACE, local_jmax, 3, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &j_l2_local, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.stage5_jn_minus_je_l2 = std::sqrt(j_l2_local);
        result.stage5_jn_minus_je_linf = local_jmax[2];
        result.current_diag.max_abs_j_charge = local_jmax[0];
        result.current_diag.max_abs_j_energy = local_jmax[1];
        result.current_diag.max_abs_j_ampere = local_jmax[0];
        result.current_diag.max_abs_j_charge_minus_ampere = 0.0;
        result.current_diag.max_abs_j_energy_minus_ampere = local_jmax[2];
        result.current_diag.e_dot_j_charge = energy_pair[1];
        result.current_diag.e_dot_j_energy = energy_pair[0];
        result.current_diag.e_dot_j_ampere = energy_pair[1];
        result.current_diag.residual_if_charge = result.energy_residual_bkg;
        result.current_diag.residual_if_ampere = result.energy_residual_bkg;

        double speed_audit[2] = {0.0, 0.0};
        for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
            const double analytic = bkg_n.cgrid.vx[idx2(j, k)];
            const double candidate = Stage5::energy_consistent_cell_speed_candidate(
                bkg_n.cgrid, bkg_n.mass, j, k, analytic);
            speed_audit[0] = std::max(speed_audit[0], std::fabs(candidate - analytic));
            speed_audit[1] = std::max(speed_audit[1], std::fabs(analytic));
        }
        MPI_Allreduce(MPI_IN_PLACE, speed_audit, 2, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        result.stage5_energy_speed_candidate_linf = speed_audit[0];
        result.stage5_energy_speed_candidate_rel = speed_audit[0] /
            std::max(Const::c * 1.0e-14, speed_audit[1]);
        double max_mid_field = 0.0;
        for (int iface = 0; iface < nxl; ++iface)
            max_mid_field = std::max(max_mid_field, std::fabs(e_mid[iface]));
        MPI_Allreduce(MPI_IN_PLACE, &max_mid_field, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        const double raw_scale = std::max(std::fabs(local_stage5[2]), 1.0);
        result.stage5_compatibility_enabled = 0;
        result.stage5_projection_skipped_zero_field =
            (max_mid_field <= 1.0e-14) ? 1 : 0;
        result.stage5_projection_skipped_raw_residual =
            (std::fabs(result.stage5_r_couple) / raw_scale >= 1.0e-2) ? 1 : 0;
        // Stage 5 is audit-only.  Record the scalar target a later bounded
        // compatibility solve would receive, but never alter fu_final here.
        result.stage5_correction_energy_target = result.stage5_r_couple;
        result.stage5_correction_energy_achieved = 0.0;
        result.stage5_correction_added_momentum_linf = 0.0;
        result.stage5_correction_l1 = 0.0;
        result.stage5_correction_linf = 0.0;
        result.stage5_correction_relative = 0.0;
        result.stage5_correction_active_bounds = 0;
        result.stage5_correction_infeasible = 0;
        result.stage5_correction_retry_count = 0;
        double limiter_sum[8] = {limiter_energy_change, u_boundary_energy,
                                 limiter_faces, limiter_active,
                                 limiter_energy_positive, limiter_energy_negative,
                                 limiter_energy_core, limiter_energy_boundary};
        MPI_Allreduce(MPI_IN_PLACE, limiter_sum, 8, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        double limiter_alpha_global = limiter_min;
        MPI_Allreduce(MPI_IN_PLACE, &limiter_alpha_global, 1, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        result.limiter_energy_defect = limiter_sum[0];
        result.limiter_energy_defect_positive = limiter_sum[4];
        result.limiter_energy_defect_negative = limiter_sum[5];
        result.limiter_energy_defect_core = limiter_sum[6];
        result.limiter_energy_defect_boundary = limiter_sum[7];
        result.flux_defect[1].boundary_energy_loss = limiter_sum[1];
        result.limiter_active_fraction = (limiter_sum[2] > 0.0)
            ? limiter_sum[3] / limiter_sum[2] : 0.0;
        result.limiter_min_alpha = limiter_alpha_global;
        result.flux_pos[0].alpha_active_fraction = result.limiter_active_fraction;
        result.flux_pos[1].alpha_active_fraction = result.limiter_active_fraction;
        result.flux_pos[0].alpha_min = limiter_alpha_global;
        result.flux_pos[1].alpha_min = limiter_alpha_global;

        double final_number = 0.0, final_ke = 0.0;
        work.total_particle_number_and_energy(final_number, final_ke);
        result.delta_ke_bkg = final_ke - initial_ke;
        result.delta_ke_beam = beam_ke_after - beam_ke_before;
        result.field_work_bkg = -dt * energy_pair[1];
        result.field_work_beam = -integrate_face_work(jbeam, beam_field, sg, dt);
        result.beam_lag_energy_residual =
            -integrate_face_work(jbeam, fields_mid, sg, dt) - result.field_work_beam;
        result.energy_pair_residual_bkg = result.energy_residual_bkg;
        result.beam_continuity_residual = std::max(
            beam_predictor.last_continuity_linf_error(),
            beam_predictor.last_boundary_flux_error());
        result.nonlinear_iterations = iter + 1;

        if (!all_finite(work, fields_new, result.j_total_face_mid)) {
            result.failed = true;
            result.state_advanced = 0;
            result.failure_reason = 4;
            result.failure_iteration = iter;
            result.failure_global_cfl = global_cfl;
            return result;
        }
        if (result.residual_E < field_tol && result.residual_J_bkg < current_tol &&
            result.beam_continuity_residual < 1.0e-6) {
            if (!transport_safe_step_gate(work, nsub))
                return result;
            capture_accepted_transport(work);
            result.species_np1 = work;
            result.beam_np1 = beam_predictor;
            result.fields_np1 = fields_new;
            result.converged = true;
            result.state_advanced = 1;
            write_stage5_diagnostics(1, 0);
            return result;
        }

        for (int iface = 0; iface < nxl; ++iface) {
            e_end[iface] = (1.0 - omega) * e_end[iface] + omega * next_e[iface];
            previous_j[iface] = result.j_bkg_face_mid[iface];
        }
        have_previous = true;
        guess = work;
        if (iter + 1 == max_iters) {
            if (!transport_safe_step_gate(last_work, nsub))
                return result;
            capture_accepted_transport(last_work);
        }
    }

    result.species_np1 = last_work;
    result.beam_np1 = beam_predictor;
    result.fields_np1 = last_fields;
    result.soft_unconverged = true;
    result.soft_accepted = true;
    result.state_advanced = 1;
    write_stage5_diagnostics(1, 1);
    return result;
}
