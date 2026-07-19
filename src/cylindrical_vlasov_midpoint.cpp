#include "vlasov_ampere_midpoint.h"
#include "nonuniform_reconstruction.h"
#include "discrete_moment_operators.h"
#include "periodic_staggered_operators.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
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

const double kFctCoreBoundaryWidth = 0.1 * Const::micro;

inline bool fct_core_x(double x)
{
    return x >= kFctCoreBoundaryWidth &&
           x <= Param::Lx - kFctCoreBoundaryWidth;
}

void close_periodic_face_blocks(std::vector<double>& faces, int nxl,
                                int block, int rank, int size, int tag)
{
    PeriodicStaggered::close_right_face_alias(faces, nxl, block,
                                              rank, size, tag);
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

void low_order_solver_checkpoint(bool enabled, const char* phase, int mpi_rank)
{
    if (!enabled) return;
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::cerr << "production_7_2_solver_checkpoint phase=" << phase
                  << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

unsigned long long hash_physical_distribution(const Species& species,
                                              const SpatialGrid& sg,
                                              int mpi_rank, int mpi_size)
{
    const int ng = sg.nghost;
    unsigned long long local = 1469598103934665603ULL ^
        static_cast<unsigned long long>(sg.ix_start + 1);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const size_t base = static_cast<size_t>(ng + ix) * Param::Nvmu;
        for (size_t q = 0; q < Param::Nvmu; ++q) {
            unsigned long long bits = 0;
            std::memcpy(&bits, &species.f[base + q], sizeof(bits));
            local ^= bits + 0x9e3779b97f4a7c15ULL + (local << 6) +
                (local >> 2);
            local *= 1099511628211ULL;
        }
    }
    unsigned long long global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                  MPI_COMM_WORLD);
    (void)mpi_rank;
    (void)mpi_size;
    return global;
}

unsigned long long hash_face_values(const std::vector<double>& values,
                                    const SpatialGrid& sg)
{
    unsigned long long local = 1469598103934665603ULL ^
        static_cast<unsigned long long>(sg.ix_start + 1);
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        unsigned long long bits = 0;
        std::memcpy(&bits, &values[static_cast<size_t>(iface)],
                    sizeof(bits));
        local ^= bits + 0x9e3779b97f4a7c15ULL + (local << 6) +
            (local >> 2);
        local *= 1099511628211ULL;
    }
    unsigned long long global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                  MPI_COMM_WORLD);
    return global;
}

unsigned long long hash_midpoint_field(const EMFields& fields,
                                       const SpatialGrid& sg)
{
    return hash_face_values(fields.Ex_face, sg);
}

}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::evaluate_production_midpoint_operator(
    const Species& bkg_n, const BeamPIC& beam_n, const EMFields& fields_n,
    const SpatialGrid& sg, double dt, double time, int mpi_rank,
    int mpi_size, int substeps_used, const Species* fixed_guess_np1,
    const EMFields* fixed_fields_end,
    const std::vector<double>* fixed_j_beam_face_mid,
    const CouplingRegionLayout* fixed_coupling_layout,
    bool fixed_candidate) const
{
    low_order_solver_checkpoint(low_order_only_, "entry", mpi_rank);
    Result result;
    reset_result(result);
    result.transport_low_order_only = low_order_only_ ? 1 : 0;
    result.substeps_used = substeps_used;
    result.species_np1 = bkg_n;
    result.beam_np1 = beam_n;
    result.fields_np1 = fields_n;
    result.operator_input_fields_end_guess = fields_n;
    low_order_solver_checkpoint(low_order_only_, "result_initialized", mpi_rank);

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
        const double x = sg.x_min +
            (static_cast<double>(result.failure_worst_ix) + 0.5) * sg.dx;
        result.failure_worst_is_core = fct_core_x(x) ? 1 : 0;
        result.failure_worst_is_tail =
            (result.failure_worst_iv >= 3 * Param::Nv / 4 ||
             result.failure_worst_imu >= 3 * Param::Nmu / 4) ? 1 : 0;
    };
    // Reused across nonlinear iterations/substeps.  These are the dominant
    // transient allocations in the new FV/FCT kernel.
    std::vector<double> fx_low(xflux_size, 0.0);
    std::vector<double> fu_low(uflux_size, 0.0);
    std::vector<double> cu_low(uflux_size, 0.0);
    const size_t high_local_cell_count = low_order_only_
        ? 0 : static_cast<size_t>(nxl) * Param::Nvmu;
    const size_t high_xflux_size = low_order_only_ ? 0 : xflux_size;
    const size_t high_uflux_size = low_order_only_ ? 0 : uflux_size;
    std::vector<double> fx_high(high_xflux_size, 0.0);
    std::vector<double> fu_high(high_uflux_size, 0.0);
    std::vector<double> cu_high(high_uflux_size, 0.0);
    // Centered u-force candidate is kept separately from the selected high
    // flux so the closure audit can isolate boundary/upwind stabilization.
    std::vector<double> cu_high_center(high_uflux_size, 0.0);
    std::vector<double> inv_cell_volume(low_order_only_ ? 0 : Param::Nvmu, 0.0);
    if (!low_order_only_) {
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                inv_cell_volume[idx2(j, k)] = 1.0 /
                    (sg.dx * bkg_n.cgrid.cell_phase_volume(j, k));
            }
        }
    }
    std::vector<double> alpha_x_left(high_local_cell_count, 1.0);
    std::vector<double> alpha_x_right(high_local_cell_count, 1.0);
    std::vector<double> alpha_u_lower(high_local_cell_count, 1.0);
    std::vector<double> alpha_u_upper(high_local_cell_count, 1.0);
    std::vector<double> left_alpha_x_right(
        low_order_only_ ? 0 : Param::Nvmu, 1.0);
    std::vector<double> fx_final(high_xflux_size, 0.0),
        fu_final(high_uflux_size, 0.0);
    std::vector<double> cu_final(high_uflux_size, 0.0);
    std::vector<double> candidate_mass(high_local_cell_count, 0.0);
    std::vector<double> donor_beta(high_local_cell_count, 1.0);
    Species low_state_buffer = bkg_n;
    std::vector<double> e_mid_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> next_e_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_jn_buffer(nface, 0.0);
    std::vector<double> integrated_jn_low_buffer(nface, 0.0);
    std::vector<double> integrated_jn_high_buffer(nface, 0.0);
    std::vector<double> integrated_je_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_je_low_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_je_high_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_je_center_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_cell_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_low_cell_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_high_cell_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_center_cell_buffer(static_cast<size_t>(nxl), 0.0);
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
    const auto current_moment_from_x_flux = [&](const std::vector<double>& flux,
                                                 std::vector<double>& current) {
        current.assign(nface, 0.0);
        if (flux.size() < high_xflux_size) return;
        #pragma omp parallel for schedule(static)
        for (int iface = 0; iface <= nxl; ++iface) {
            double gamma = 0.0;
            for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k)
                    gamma += flux[xface_index(iface, j, k)];
            current[static_cast<size_t>(iface)] = bkg_n.charge * gamma;
        }
    };
    const auto current_moment_from_x_flux_by_u = [&]
        (const std::vector<double>& flux, std::vector<double>& current) {
        current.assign(nface * static_cast<size_t>(Param::Nv), 0.0);
        if (flux.size() < high_xflux_size) return;
        #pragma omp parallel for schedule(static)
        for (int iface = 0; iface <= nxl; ++iface) {
            for (int j = 0; j < Param::Nv; ++j) {
                double gamma = 0.0;
                for (int k = 0; k < Param::Nmu; ++k)
                    gamma += flux[xface_index(iface, j, k)];
                current[static_cast<size_t>(iface) * Param::Nv + j] =
                    bkg_n.charge * gamma;
            }
        }
    };
    result.j_bkg_face_mid.assign(nface, 0.0);
    result.j_bkg_face_low_mid.assign(nface, 0.0);
    result.j_bkg_face_center_mid.assign(nface, 0.0);
    result.j_bkg_face_high_mid.assign(nface, 0.0);
    result.j_beam_face_mid.assign(nface, 0.0);
    result.j_total_face_mid.assign(nface, 0.0);
    result.j_bkg_energy_debug_face.assign(nface, 0.0);
    result.j_bkg_energy_low_debug_face.assign(nface, 0.0);
    result.j_bkg_energy_center_debug_face.assign(nface, 0.0);
    result.j_bkg_energy_high_debug_face.assign(nface, 0.0);
    result.j_bkg_energy_cell_mid.assign(static_cast<size_t>(nxl), 0.0);

    // The production low-order Ampere benchmark deliberately omits all
    // high-order/FCT storage.  Checkpoints are test-only and use the same
    // collective sequence on every rank.
    low_order_solver_checkpoint(low_order_only_, "work_arrays_initialized",
                                mpi_rank);

    std::vector<double> e_end(static_cast<size_t>(nxl), 0.0);
    for (int iface = 0; iface < nxl; ++iface) {
        e_end[static_cast<size_t>(iface)] = fixed_fields_end
            ? fixed_fields_end->Ex_face[iface] : fields_n.Ex_face[iface];
    }
    Species guess = fixed_guess_np1 ? *fixed_guess_np1 : bkg_n;
    Species last_work = bkg_n;
    EMFields last_fields = fields_n;
    std::vector<double> previous_j(static_cast<size_t>(nxl), 0.0);
    bool have_previous = false;

    // Stage 7 is deliberately not enabled: the beam is a controlled
    // predictor during phases 1--5, and its lag is kept explicit.
    BeamPIC beam_predictor = beam_n;
    EMFields beam_field = fields_n;
    double beam_ke_before = 0.0;
    double beam_ke_after = 0.0;
    if (beam_enabled_ && !fixed_j_beam_face_mid) {
        beam_predictor.begin_step(sg, dt);
        beam_predictor.inject(sg, beam_field, dt, time, mpi_rank, mpi_size);
        beam_ke_before = beam_predictor.total_kinetic_energy();
        beam_predictor.push(sg, beam_field, dt, mpi_rank, mpi_size);
        beam_ke_after = beam_predictor.total_kinetic_energy();
        beam_predictor.deposit_density(sg, mpi_rank, mpi_size);
        beam_predictor.finalize_charge_conserving_current(sg, dt, mpi_rank,
                                                          mpi_size);
    }
    low_order_solver_checkpoint(low_order_only_, "beam_predictor_ready",
                                mpi_rank);
    std::vector<double> jbeam(nface, 0.0);
    if (fixed_j_beam_face_mid) {
        for (size_t i = 0; i < std::min(jbeam.size(),
                                        fixed_j_beam_face_mid->size()); ++i) {
            jbeam[i] = (*fixed_j_beam_face_mid)[i];
        }
    } else if (beam_enabled_) {
        for (size_t i = 0; i < std::min(jbeam.size(),
                                        beam_predictor.current_face_x.size()); ++i) {
            jbeam[i] = beam_predictor.current_face_x[i];
        }
    }

    const int max_iters = fixed_candidate ? 1 : max_midpoint_iterations_;
    if (midpoint_iteration_trace_for_test_) {
        result.midpoint_residual_e_history.reserve(
            static_cast<size_t>(max_iters));
        result.midpoint_residual_j_bkg_history.reserve(
            static_cast<size_t>(max_iters));
        result.midpoint_residual_j_beam_history.reserve(
            static_cast<size_t>(max_iters));
        result.midpoint_residual_f_history.reserve(
            static_cast<size_t>(max_iters));
    }
    const double field_tol = 1.0e-6;
    const double current_tol = 1.0e-5;
    const double omega = 0.55;
    for (int iter = 0; iter < max_iters; ++iter) {
        if (iter == 0)
            low_order_solver_checkpoint(low_order_only_, "iteration_begin", mpi_rank);
        // Preserve the exact endpoint guess used by this kernel call.  The
        // accepted field is an Ampere output and is not a valid replacement
        // for fixed-state operator replay.
        if (capture_midpoint_input_) {
            result.operator_input_guess = guess;
            result.operator_input_fields_end_guess = fields_n;
            for (int iface = 0; iface < nxl; ++iface) {
                result.operator_input_fields_end_guess.Ex_face[iface] = e_end[iface];
            }
            result.operator_input_fields_end_guess.sync_cell_ex_from_faces(
                mpi_rank, mpi_size);
        }
        std::vector<double>& e_mid = e_mid_buffer;
        for (int iface = 0; iface < nxl; ++iface) {
            e_mid[iface] = 0.5 * (fields_n.Ex_face[iface] + e_end[iface]);
        }
        EMFields fields_mid;
        set_midpoint_field(fields_mid, fields_n, e_mid, sg, mpi_rank, mpi_size);
        if (iter == 0)
            low_order_solver_checkpoint(low_order_only_, "midpoint_field_ready", mpi_rank);

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
        // Report the real force-subcycle count rather than only the outer
        // retry count supplied at entry.
        result.substeps_used = nsub;
        if (nsub > 8 || !std::isfinite(global_cfl)) {
            result.failed = true;
            result.state_advanced = 0;
            result.failure_reason = 1;
            result.failure_iteration = iter;
            result.failure_global_cfl = global_cfl;
            result.x_low_failure_kind = 3; // genuine combined CFL failure
            return result;
        }
        if (iter == 0)
            low_order_solver_checkpoint(low_order_only_, "cfl_ready", mpi_rank);

        // Every transport current is evaluated from the same Picard midpoint
        // state.  The updated work state remains the conservative endpoint
        // f^(n+1)=f^n-dt L(f^(n+1/2)).
        Species midpoint_state = bkg_n;
        #pragma omp parallel for schedule(static)
        for (long long p = 0;
             p < static_cast<long long>(midpoint_state.f.size()); ++p) {
            midpoint_state.f[static_cast<size_t>(p)] = 0.5 * (
                bkg_n.f[static_cast<size_t>(p)] +
                guess.f[static_cast<size_t>(p)]);
        }
        exchange_ghosts_x_persistent(midpoint_state, sg, mpi_rank, mpi_size);
        Species work = bkg_n;
        if (fixed_candidate) {
            // Section 11.5 audit: report exactly which state/field each
            // production transport layer reads.  The low donor base is
            // intentionally step-start; both high-order x and u candidates
            // read the same Picard midpoint state and midpoint field.
            result.x_low_state_hash = hash_physical_distribution(
                work, sg, mpi_rank, mpi_size);
            result.u_low_state_hash = result.x_low_state_hash;
            result.x_high_state_hash = hash_physical_distribution(
                midpoint_state, sg, mpi_rank, mpi_size);
            result.u_high_state_hash = result.x_high_state_hash;
            // x transport has no E input; the u donor and both high layers
            // read the same midpoint E.
            result.x_low_field_hash = 0ULL;
            result.u_low_field_hash = hash_midpoint_field(fields_mid, sg);
            result.x_high_field_hash = result.u_low_field_hash;
            result.u_high_field_hash = result.u_low_field_hash;
            result.start_field_hash = hash_midpoint_field(fields_n, sg);
            result.end_field_hash = hash_face_values(e_end, sg);
            result.x_low_time_layer = 0;
            result.u_low_time_layer = 0;
            result.x_high_time_layer = 1;
            result.u_high_time_layer = 1;
        }
        std::vector<double>& integrated_jn = integrated_jn_buffer;
        std::vector<double>& integrated_jn_low = integrated_jn_low_buffer;
        std::vector<double>& integrated_jn_high = integrated_jn_high_buffer;
        std::vector<double>& integrated_je_cell = integrated_je_buffer;
        std::vector<double>& integrated_je_low_cell = integrated_je_low_buffer;
        std::vector<double>& integrated_je_high_cell = integrated_je_high_buffer;
        std::vector<double>& integrated_je_center_cell = integrated_je_center_buffer;
        std::fill(integrated_jn.begin(), integrated_jn.end(), 0.0);
        std::fill(integrated_jn_low.begin(), integrated_jn_low.end(), 0.0);
        std::fill(integrated_jn_high.begin(), integrated_jn_high.end(), 0.0);
        std::fill(integrated_je_cell.begin(), integrated_je_cell.end(), 0.0);
        std::fill(integrated_je_low_cell.begin(), integrated_je_low_cell.end(), 0.0);
        std::fill(integrated_je_high_cell.begin(), integrated_je_high_cell.end(), 0.0);
        std::fill(integrated_je_center_cell.begin(), integrated_je_center_cell.end(), 0.0);
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
        double limiter_faces_core = 0.0;
        double limiter_faces_boundary = 0.0;
        double limiter_active_core = 0.0;
        double limiter_active_boundary = 0.0;
        double limiter_min_core = 1.0;
        double limiter_min_boundary = 1.0;
        double fct_budget_violation = 0.0;

        const double h = dt / nsub;
        const auto capture_accepted_transport =
            [&](const Species& committed) {
                // low_order_only_ intentionally leaves high/FCT work arrays
                // empty.  The accepted-state audit must therefore use the
                // same low-order layers as both its high and final views.
                const std::vector<double>& audit_fx_high = low_order_only_
                    ? fx_low : fx_high;
                const std::vector<double>& audit_fu_high = low_order_only_
                    ? fu_low : fu_high;
                const std::vector<double>& audit_fx_final = low_order_only_
                    ? fx_low : fx_final;
                const std::vector<double>& audit_fu_final = low_order_only_
                    ? fu_low : fu_final;
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
                                    h * (audit_fx_final[xface_index(ix, j, k)] -
                                         fx_low[xface_index(ix, j, k)]),
                                    h * (audit_fx_final[xface_index(ix + 1, j, k)] -
                                         fx_low[xface_index(ix + 1, j, k)]),
                                    h * (audit_fu_final[uface_index(ix, j, k)] -
                                         fu_low[uface_index(ix, j, k)]),
                                    h * (audit_fu_final[uface_index(ix, j + 1, k)] -
                                         fu_low[uface_index(ix, j + 1, k)])
                                };
                                const double a_raw[4] = {
                                    h * (audit_fx_high[xface_index(ix, j, k)] -
                                         fx_low[xface_index(ix, j, k)]),
                                    h * (audit_fx_high[xface_index(ix + 1, j, k)] -
                                         fx_low[xface_index(ix + 1, j, k)]),
                                    h * (audit_fu_high[uface_index(ix, j, k)] -
                                         fu_low[uface_index(ix, j, k)]),
                                    h * (audit_fu_high[uface_index(ix, j + 1, k)] -
                                         fu_low[uface_index(ix, j + 1, k)])
                                };
                                record[0] = static_cast<double>(
                                    (sg.ix_start + ix) * Param::Nvmu + idx2(j, k));
                                record[1] = low_state_buffer.f[
                                    mass_index(ng + ix, j, k)];
                                record[2] = low_order_only_ ? m_safe : candidate_mass[p];
                                record[3] = m_safe;
                                record[4] = std::max(0.0, -a_lim[0]) +
                                    std::max(0.0, a_lim[1]) +
                                    std::max(0.0, -a_lim[2]) +
                                    std::max(0.0, a_lim[3]);
                                record[5] = std::max(0.0, a_lim[0]) +
                                    std::max(0.0, -a_lim[1]) +
                                    std::max(0.0, a_lim[2]) +
                                    std::max(0.0, -a_lim[3]);
                                record[6] = low_order_only_ ? 1.0 : donor_beta[p];
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
                double local_negative_total = 0.0;
                double local_negative_mass_core = 0.0;
                double local_positive_mass_core = 0.0;
                double local_negative_energy_core = 0.0;
                double local_positive_energy_core = 0.0;
                double local_negative_current_core = 0.0;
                double local_positive_current_core = 0.0;
                long long local_negative_core_cells = 0;
                #pragma omp parallel for reduction(min:local_min) reduction(max:local_nonfinite) reduction(+:local_negative_total,local_negative_mass_core,local_positive_mass_core,local_negative_energy_core,local_positive_energy_core,local_negative_current_core,local_positive_current_core,local_negative_core_cells)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double m = candidate_state.f[
                                mass_index(ng + ix, j, k)];
                            local_min = std::min(local_min, m);
                            if (!std::isfinite(m)) local_nonfinite = 1.0;
                            const bool core = fct_core_x(sg.x(ng + ix));
                            const double abs_mass = std::fabs(m);
                            const double energy = bkg_n.cgrid.kinetic_energy[idx2(j, k)];
                            const double current_weight = std::fabs(
                                bkg_n.charge * bkg_n.cgrid.vx[idx2(j, k)]);
                            if (m < 0.0) {
                                local_negative_total += -m;
                                if (core) {
                                    local_negative_mass_core += -m;
                                    local_negative_energy_core += energy * abs_mass;
                                    local_negative_current_core += current_weight * abs_mass;
                                    local_negative_core_cells += 1;
                                }
                            } else if (core) {
                                local_positive_mass_core += m;
                                local_positive_energy_core += energy * m;
                                local_positive_current_core += current_weight * m;
                            }
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
                double debt_sums[7] = {
                    local_negative_total, local_negative_mass_core,
                    local_positive_mass_core, local_negative_energy_core,
                    local_positive_energy_core, local_negative_current_core,
                    local_positive_current_core};
                MPI_Allreduce(MPI_IN_PLACE, debt_sums, 7, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);
                long long global_negative_core_cells = local_negative_core_cells;
                MPI_Allreduce(MPI_IN_PLACE, &global_negative_core_cells, 1,
                              MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
                result.neg_mass_total_guard = debt_sums[0];
                result.neg_mass_core = debt_sums[1];
                result.neg_mass_core_fraction = debt_sums[1] /
                    std::max(std::numeric_limits<double>::min(), debt_sums[2]);
                result.neg_energy_core_abs = debt_sums[3];
                result.neg_energy_core_fraction = debt_sums[3] /
                    std::max(std::numeric_limits<double>::min(), debt_sums[4]);
                result.neg_current_core_abs = debt_sums[5];
                result.neg_current_core_fraction = debt_sums[5] /
                    std::max(std::numeric_limits<double>::min(), debt_sums[6]);
                result.neg_cell_core = global_negative_core_cells;
                const bool core_macro_debt =
                    result.neg_mass_core_fraction > 1.0e-6 ||
                    result.neg_energy_core_fraction > 1.0e-6 ||
                    result.neg_current_core_fraction > 1.0e-6;
                result.negative_debt_level = core_macro_debt ? NEG_DEBT_ABORT :
                    (debt_sums[0] > 0.0 ? NEG_DEBT_WARN : NEG_DEBT_OK);
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
                if (local_nonfinite == 0.0 &&
                    allow_finite_negative_debt_for_test_ && !core_macro_debt) {
                    result.trial_failure_downgraded = 1;
                    result.accepted_with_negative_debt = 1;
                    result.failure_final_min = global_owner.value;
                    if (global_linear >= 0)
                        set_failure_location(global_linear);
                    return true;
                }
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
            CouplingSubstepSeamAudit substep_audit = {};
            if (fixed_candidate) {
                substep_audit.substep = sub;
                substep_audit.dt_substep = h;
                substep_audit.e_face_mid = fields_mid.Ex_face;
            }
            if (iter == 0 && sub == 0)
                low_order_solver_checkpoint(low_order_only_, "substep_begin", mpi_rank);
            if (sub == 0)
                exchange_ghosts_x_persistent(work, sg, mpi_rank, mpi_size);
            if (fixed_candidate) {
                result.x_low_state_hash_history.push_back(
                    hash_physical_distribution(work, sg, mpi_rank, mpi_size));
                result.u_low_state_hash_history.push_back(
                    result.x_low_state_hash_history.back());
                result.x_high_state_hash_history.push_back(
                    hash_physical_distribution(midpoint_state, sg, mpi_rank,
                                               mpi_size));
                result.u_high_state_hash_history.push_back(
                    result.x_high_state_hash_history.back());
                result.u_field_hash_history.push_back(
                    hash_midpoint_field(fields_mid, sg));
            }
            std::fill(fx_low.begin(), fx_low.end(), 0.0);
            std::fill(fx_high.begin(), fx_high.end(), 0.0);
            std::fill(fu_low.begin(), fu_low.end(), 0.0);
            std::fill(fu_high.begin(), fu_high.end(), 0.0);
            std::fill(cu_low.begin(), cu_low.end(), 0.0);
            std::fill(cu_high.begin(), cu_high.end(), 0.0);
            std::fill(cu_high_center.begin(), cu_high_center.end(), 0.0);

            // Unsplit donor-cell FCT baseline.  Unlike the high-order
            // midpoint candidate, this is an explicit positive update of the
            // current substep state `work`; its donor mass and update state
            // are therefore identical.
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
            if (fixed_candidate)
                current_moment_from_x_flux(fx_low,
                    substep_audit.jn_low_pre_sync);
            close_periodic_face_blocks(fx_low, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 901);
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                const double a = bkg_n.charge * fields_mid.Ex[ng + ix] /
                                 (bkg_n.mass * Const::c);
                for (int jf = 1; jf < Param::Nv; ++jf) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const int donor = (a >= 0.0) ? jf - 1 : jf;
                        const double donor_du = bkg_n.cgrid.upar_widths[donor];
                        // cu_low stores C_u, not Phi_u.  Since f stores the
                        // cell mass M, this coefficient retains dx and the
                        // perpendicular ring measure through M/du_face.
                        const double coefficient = Stage5::donor_cell_coefficient(
                            work.f[mass_index(ng + ix, donor, k)], donor_du);
                        const size_t id = uface_index(ix, jf, k);
                        cu_low[id] = coefficient;
                        fu_low[id] = a * coefficient;
                    }
                }
                // Open velocity-domain boundaries: only outward donor flux
                // is present.  The boundary donor uses its own cell width,
                // and the existing boundary ledger records the resulting
                // particle, momentum and energy loss.
                if (a < 0.0) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = uface_index(ix, 0, k);
                        cu_low[id] = Stage5::donor_cell_coefficient(
                            work.f[mass_index(ng + ix, 0, k)],
                            bkg_n.cgrid.upar_widths[0]);
                        fu_low[id] = a * cu_low[id];
                    }
                } else if (a > 0.0) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = uface_index(ix, Param::Nv, k);
                        cu_low[id] = Stage5::donor_cell_coefficient(
                            work.f[mass_index(ng + ix, Param::Nv - 1, k)],
                            bkg_n.cgrid.upar_widths[Param::Nv - 1]);
                        fu_low[id] = a * cu_low[id];
                    }
                }
            }
            if (iter == 0 && sub == 0)
                low_order_solver_checkpoint(low_order_only_, "low_fluxes_ready", mpi_rank);

            if (low_order_only_) {
                // Test-only production mode: high-order reconstruction and
                // CTU are not evaluated.  Shared accounting below aliases
                // low fluxes as the high/final diagnostic layers.
            } else {
            // Reconstruct physical cell-average fbar rather than integrated
            // mass M.  Both directions read the same substep-start state and
            // write one MUSCL Riemann flux per shared face.
            const auto x_center = [&](int storage_ix) {
                return sg.x_min + (sg.ix_start + storage_ix - ng + 0.5) * sg.dx;
            };
            const auto fbar = [&](int storage_ix, int j, int k) {
                return midpoint_state.f[mass_index(storage_ix, j, k)] *
                    inv_cell_volume[idx2(j, k)];
            };

            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                const int il = ng + iface - 1;
                const int ir = ng + iface;
                const double s_face = sg.x_min + (sg.ix_start + iface) * sg.dx;
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const NonuniformMuscl::FaceStates states =
                            NonuniformMuscl::reconstruct_face(
                                fbar(il - 1, j, k), fbar(il, j, k),
                                fbar(ir, j, k), fbar(ir + 1, j, k),
                                x_center(il - 1), x_center(il),
                                x_center(ir), x_center(ir + 1), s_face);
                        const double analytic_speed =
                            bkg_n.cgrid.vx[idx2(j, k)];
                        const double transport_speed =
                            energy_consistent_x_high_velocity_for_test_
                            ? Stage5::energy_consistent_cell_speed_candidate(
                                bkg_n.cgrid, bkg_n.mass, j, k,
                                analytic_speed)
                            : analytic_speed;
                        const double state = NonuniformMuscl::upwind_state(
                            states, transport_speed);
                        fx_high[xface_index(iface, j, k)] =
                            transport_speed * state *
                            bkg_n.cgrid.upar_widths[j] *
                            bkg_n.cgrid.uperp_ring_areas[k];
                    }
                }
            }
            if (fixed_candidate)
                current_moment_from_x_flux(fx_high,
                    substep_audit.jn_high_pre_sync);
            close_periodic_face_blocks(fx_high, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 908);

            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                const int storage_ix = ng + ix;
                const double a = bkg_n.charge * fields_mid.Ex[storage_ix] /
                                 (bkg_n.mass * Const::c);
                for (int jf = 1; jf < Param::Nv; ++jf) {
                    const int jl = jf - 1;
                    const int jr = jf;
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const double area = bkg_n.cgrid.uperp_ring_areas[k];
                        const int jll = (jl == 0) ? jl : jl - 1;
                        const int jrr = (jr + 1 == Param::Nv) ? jr : jr + 1;
                        const double s_jll = (jl == 0)
                            ? bkg_n.cgrid.upar_cells[jl] - bkg_n.cgrid.upar_widths[jl]
                            : bkg_n.cgrid.upar_cells[jll];
                        const double s_jrr = (jr + 1 == Param::Nv)
                            ? bkg_n.cgrid.upar_cells[jr] + bkg_n.cgrid.upar_widths[jr]
                            : bkg_n.cgrid.upar_cells[jrr];
                        const NonuniformMuscl::FaceStates states =
                            NonuniformMuscl::reconstruct_face(
                                fbar(storage_ix, jll, k),
                                fbar(storage_ix, jl, k),
                                fbar(storage_ix, jr, k),
                                fbar(storage_ix, jrr, k),
                                s_jll,
                                bkg_n.cgrid.upar_cells[jl],
                                bkg_n.cgrid.upar_cells[jr],
                                s_jrr,
                                bkg_n.cgrid.upar_faces[jf]);
                        const size_t id = uface_index(ix, jf, k);
                        // Production high order is centered at every x.
                        // Donor-cell transport exists only in the separate
                        // FCT low-order base.  The former boundary-upwind
                        // high candidate is retained behind a test-only A/B
                        // switch and never selected by normal production.
                        const double centered =
                            NonuniformMuscl::centered_state(states);
                        cu_high_center[id] =
                            NonuniformMuscl::upar_face_coefficient(
                                centered, sg.dx, area);
                        cu_high[id] = cu_high_center[id];
                        if (legacy_boundary_upwind_high_candidate_for_test_ &&
                            !fct_core_x(sg.x(storage_ix))) {
                            const double upwind =
                                NonuniformMuscl::upwind_state(states, a);
                            cu_high[id] = NonuniformMuscl::upar_face_coefficient(
                                upwind, sg.dx, area);
                        }
                        fu_high[id] = NonuniformMuscl::upar_face_flux(a, cu_high[id]);
                    }
                }
                // Keep the physical velocity-domain outflow exactly shared
                // with the donor baseline; no one-sided MUSCL stencil exists
                // outside the truncated u_parallel domain.
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t lower = uface_index(ix, 0, k);
                    const size_t upper = uface_index(ix, Param::Nv, k);
                    cu_high[lower] = cu_low[lower];
                    cu_high[upper] = cu_low[upper];
                    fu_high[lower] = fu_low[lower];
                    fu_high[upper] = fu_low[upper];
                    cu_high_center[lower] = cu_low[lower];
                    cu_high_center[upper] = cu_low[upper];
                }
            }
            } // !low_order_only_: nonuniform MUSCL high-order candidate

            Species& low_state = low_state_buffer;
            low_state = work;
            const auto low_roundoff_bound = [](double m_old, double tx_left,
                                                double tx_right, double tu_lower,
                                                double tu_upper) {
                // Covers the old-mass remainder and every incoming/outgoing
                // transfer in the explicit donor-cell summation.  The
                // subnormal allowance is needed in empty velocity tails.
                const double scale = std::fabs(m_old) + std::fabs(tx_left) +
                    std::fabs(tx_right) + std::fabs(tu_lower) +
                    std::fabs(tu_upper);
                return 256.0 * std::numeric_limits<double>::epsilon() * scale +
                    64.0 * std::numeric_limits<double>::denorm_min();
            };
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const double m_old = work.f[mass_index(ng + ix, j, k)];
                        const double tx_left = h * fx_low[xface_index(ix, j, k)];
                        const double tx_right = h * fx_low[xface_index(ix + 1, j, k)];
                        const double tu_lower = h * fu_low[uface_index(ix, j, k)];
                        const double tu_upper = h * fu_low[uface_index(ix, j + 1, k)];
                        // Build the donor-cell state as a sum of nonnegative
                        // terms, rather than a cancellation of full face
                        // transfers.  With the enforced combined CFL, this is
                        // M_old*(1-C_out) plus the four incoming transfers.
                        const double incoming = std::max(0.0, tx_left) +
                            std::max(0.0, -tx_right) +
                            std::max(0.0, tu_lower) +
                            std::max(0.0, -tu_upper);
                        // outgoing/M_old is evaluated through the analytic
                        // donor CFL coefficients, avoiding a division by a
                        // tiny tail mass.
                        const double a = bkg_n.charge *
                            fields_mid.Ex[ng + ix] / (bkg_n.mass * Const::c);
                        const double c_out = h * (
                            std::fabs(bkg_n.cgrid.vx[idx2(j, k)]) / sg.dx +
                            std::fabs(a) / bkg_n.cgrid.upar_widths[j]);
                        low_state.f[mass_index(ng + ix, j, k)] =
                            m_old * (1.0 - c_out) + incoming;
                    }
                }
            }
            // The low-order state is the FCT positivity base, but only in
            // FCT/low-order modes.  Normalize locally bounded roundoff debt
            // before it enters a donor budget; no-FCT keeps this state intact
            // as a diagnostic-only comparison layer.
            double low_roundoff_zeroed_count = 0.0;
            double low_roundoff_zeroed_mass = 0.0;
            if (fct_enabled_ || low_order_only_) {
                #pragma omp parallel for schedule(static) reduction(+:low_roundoff_zeroed_count,low_roundoff_zeroed_mass)
                for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t state_id = mass_index(ng + ix, j, k);
                        double& m = low_state.f[state_id];
                        const double tolerance = low_roundoff_bound(
                            work.f[state_id],
                            h * fx_low[xface_index(ix, j, k)],
                            h * fx_low[xface_index(ix + 1, j, k)],
                            h * fu_low[uface_index(ix, j, k)],
                            h * fu_low[uface_index(ix, j + 1, k)]);
                        if (std::isfinite(m) && m < 0.0 && m >= -tolerance) {
                            low_roundoff_zeroed_count += 1.0;
                            low_roundoff_zeroed_mass += -m;
                            m = 0.0;
                        }
                    }
                MPI_Allreduce(MPI_IN_PLACE, &low_roundoff_zeroed_count, 1,
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &low_roundoff_zeroed_mass, 1,
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                result.low_order_roundoff_zeroed_count +=
                    static_cast<long long>(low_roundoff_zeroed_count);
                result.low_order_roundoff_zeroed_mass += low_roundoff_zeroed_mass;
            }
            exchange_ghosts_x_persistent(low_state, sg, mpi_rank, mpi_size);
            if (iter == 0 && sub == 0)
                low_order_solver_checkpoint(low_order_only_, "low_state_ghosts_ready", mpi_rank);

            double low_min = std::numeric_limits<double>::infinity();
            double low_scale = 0.0;
            double low_tolerance_linf = 0.0;
            double low_positivity_violation = 0.0;
            double low_nonfinite = 0.0;
            double low_negative_mass = 0.0;
            long long low_negative_count = 0;
            #pragma omp parallel for reduction(min:low_min) reduction(max:low_scale,low_tolerance_linf,low_positivity_violation,low_nonfinite) reduction(+:low_negative_mass,low_negative_count)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double m = low_state.f[mass_index(ng + ix, j, k)];
                    const double tolerance = low_roundoff_bound(
                        work.f[mass_index(ng + ix, j, k)],
                        h * fx_low[xface_index(ix, j, k)],
                        h * fx_low[xface_index(ix + 1, j, k)],
                        h * fu_low[uface_index(ix, j, k)],
                        h * fu_low[uface_index(ix, j + 1, k)]);
                    low_min = std::min(low_min, m);
                    low_scale = std::max(low_scale, std::fabs(m));
                    low_tolerance_linf = std::max(low_tolerance_linf, tolerance);
                    low_positivity_violation = std::max(
                        low_positivity_violation, -m - tolerance);
                    if (m < 0.0) {
                        low_negative_mass += -m;
                        low_negative_count += 1;
                    }
                    if (!std::isfinite(m)) low_nonfinite = 1.0;
                }
            double low_guard[5] = {low_min, low_scale, low_tolerance_linf,
                                   low_positivity_violation, low_nonfinite};
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[0], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[1], 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &low_negative_mass, 1, MPI_DOUBLE,
                          MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &low_negative_count, 1,
                          MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
            const double mass_roundoff = low_guard[2];
            result.fct_low_order_tolerance_linf = std::max(
                result.fct_low_order_tolerance_linf, mass_roundoff);
            result.low_order_candidate_min = std::min(
                result.low_order_candidate_min, low_guard[0]);
            result.low_order_negative_mass += low_negative_mass;
            result.low_order_negative_count += low_negative_count;
            const bool no_fct_diagnostic_only = !low_order_only_ && !fct_enabled_ &&
                allow_finite_negative_debt_for_test_;
            const bool low_order_is_final_or_fct_base = !no_fct_diagnostic_only;
            if (low_guard[4] != 0.0 || !std::isfinite(low_guard[0]) ||
                (low_order_is_final_or_fct_base && low_guard[3] > 0.0)) {
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

            if (controlled_fct_flux_injection_enabled_ && fct_enabled_) {
                double local_injection_count = 0.0;
                if (mpi_rank == 0 && nxl >= 3 && h > 0.0) {
                    // One owned, internal x face is sufficient for this test.
                    // The same added flux leaves its left cell and enters its
                    // right cell, so it is conservative before FCT limits it.
                    const int iface = nxl / 2;
                    const int donor_ix = iface - 1;
                    int donor_j = 0;
                    int donor_k = 0;
                    double donor_mass = -1.0;
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double mass = low_state.f[
                                mass_index(ng + donor_ix, j, k)];
                            if (mass > donor_mass) {
                                donor_mass = mass;
                                donor_j = j;
                                donor_k = k;
                            }
                        }
                    if (donor_mass > 0.0) {
                        const size_t face = xface_index(iface, donor_j, donor_k);
                        fx_high[face] += 1.25 * donor_mass / h;
                        local_injection_count = 1.0;
                    }
                }
                MPI_Allreduce(MPI_IN_PLACE, &local_injection_count, 1,
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                result.fct_controlled_injection_count +=
                    static_cast<long long>(local_injection_count);
            }

            if (low_order_only_) {
                // Do not enter FCT/donor-capacity processing in the low-order
                // production test.  Commit the already checked low-order
                // conservative state and retain identical low/high/final
                // diagnostic layers for the shared accounting below.
                work = low_state;
                if (!transport_safe_step_gate(work, sub + 1))
                    return result;
                if (iter == 0 && sub == 0)
                    low_order_solver_checkpoint(low_order_only_, "low_state_committed", mpi_rank);
                if (fixed_candidate)
                    substep_audit.jn_final_pre_sync =
                        substep_audit.jn_low_pre_sync;
            } else if (!fct_enabled_) {
                // High-order/no-FCT verification path.  The MUSCL face fluxes
                // are the final conservative fluxes: the same arrays update
                // M, construct J_N/J_E and supply every closure diagnostic.
                fx_final = fx_high;
                fu_final = fu_high;
                cu_final = cu_high;
                if (fixed_candidate)
                    substep_audit.jn_final_pre_sync =
                        substep_audit.jn_high_pre_sync;
                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix) {
                    for (int j = 0; j < Param::Nv; ++j) {
                        for (int k = 0; k < Param::Nmu; ++k) {
                            work.f[mass_index(ng + ix, j, k)] -= h *
                                (fx_final[xface_index(ix + 1, j, k)] -
                                 fx_final[xface_index(ix, j, k)] +
                                 fu_final[uface_index(ix, j + 1, k)] -
                                 fu_final[uface_index(ix, j, k)]);
                        }
                    }
                }
                exchange_ghosts_x_persistent(work, sg, mpi_rank, mpi_size);
                if (!transport_safe_step_gate(work, sub + 1))
                    return result;
            } else {
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
            // The FCT combination preserves the C_u contract exactly:
            // cu_final is a coefficient, while fu_final is its acceleration-
            // multiplied mass flux.  No energy-current formula is altered.
            cu_final = cu_low;
            #pragma omp parallel for schedule(static) reduction(+:limiter_faces,limiter_active,limiter_faces_core,limiter_faces_boundary,limiter_active_core,limiter_active_boundary) reduction(min:limiter_min,limiter_min_core,limiter_min_boundary)
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
                    const bool core_face = fct_core_x(
                        sg.x_min + (sg.ix_start + iface) * sg.dx);
                    if (alpha < 1.0 - 1.0e-14) limiter_active += 1.0;
                    limiter_min = std::min(limiter_min, alpha);
                    if (core_face) {
                        limiter_faces_core += 1.0;
                        if (alpha < 1.0 - 1.0e-14) limiter_active_core += 1.0;
                        limiter_min_core = std::min(limiter_min_core, alpha);
                    } else {
                        limiter_faces_boundary += 1.0;
                        if (alpha < 1.0 - 1.0e-14) limiter_active_boundary += 1.0;
                        limiter_min_boundary = std::min(limiter_min_boundary, alpha);
                    }
                }
            }
            if (fixed_candidate)
                current_moment_from_x_flux(fx_final,
                    substep_audit.jn_final_pre_sync);
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
            #pragma omp parallel for schedule(static) reduction(+:limiter_faces,limiter_active,limiter_faces_core,limiter_faces_boundary,limiter_active_core,limiter_active_boundary) reduction(min:limiter_min,limiter_min_core,limiter_min_boundary)
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
                    limiter_faces += 1.0;
                    if (alpha < 1.0 - 1.0e-14) limiter_active += 1.0;
                    limiter_min = std::min(limiter_min, alpha);
                    if (fct_core_x(sg.x(ng + ix))) {
                        limiter_faces_core += 1.0;
                        if (alpha < 1.0 - 1.0e-14) limiter_active_core += 1.0;
                        limiter_min_core = std::min(limiter_min_core, alpha);
                    } else {
                        limiter_faces_boundary += 1.0;
                        if (alpha < 1.0 - 1.0e-14) limiter_active_boundary += 1.0;
                        limiter_min_boundary = std::min(limiter_min_boundary, alpha);
                    }
                }

            // P1: construct the final state in scratch from the low-order
            // state plus limited anti-diffusive mass transfers.  The four
            // transfer terms below are h * (Phi^H - Phi^L), hence have the
            // same units as the cell-integrated cylindrical mass.
            double candidate_min = std::numeric_limits<double>::infinity();
            double high_candidate_min = std::numeric_limits<double>::infinity();
            double high_candidate_donor_excess = 0.0;
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
            #pragma omp for schedule(static) reduction(min:candidate_min,high_candidate_min) reduction(max:high_candidate_donor_excess,high_low_identity_linf,high_low_identity_violation,final_tolerance_linf,final_positivity_violation,candidate_nonfinite,donor_roundoff_warning) reduction(+:donor_beta_count,roundoff_normalized_count,roundoff_normalized_mass) reduction(min:donor_beta_min)
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
                    if (fct_activation_audit_enabled_) {
                        const double high_candidate = m_low + high_minus_low;
                        if (std::isfinite(high_candidate)) {
                            high_candidate_min = std::min(high_candidate_min,
                                                          high_candidate);
                        } else {
                            candidate_nonfinite = 1.0;
                        }
                    }
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
                    if (fct_activation_audit_enabled_) {
                        const long double high_outflow_ld =
                            std::max(0.0L, -static_cast<long double>(ax_left_raw)) +
                            std::max(0.0L, static_cast<long double>(ax_right_raw)) +
                            std::max(0.0L, -static_cast<long double>(au_lower_raw)) +
                            std::max(0.0L, static_cast<long double>(au_upper_raw));
                        high_candidate_donor_excess = std::max(
                            high_candidate_donor_excess,
                            static_cast<double>(std::max(0.0L, high_outflow_ld -
                                std::max(0.0L, static_cast<long double>(m_low)))));
                    }
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
                if (fixed_candidate)
                    current_moment_from_x_flux(fx_final,
                        substep_audit.jn_final_pre_sync);
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

            double invariant_max[9] = {
                high_low_identity_linf, high_low_identity_violation,
                donor_capacity_violation, interface_checksum_linf,
                interface_checksum_violation, final_tolerance_linf,
                final_positivity_violation, candidate_nonfinite,
                high_candidate_donor_excess
            };
            MPI_Allreduce(MPI_IN_PLACE, invariant_max, 9, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &candidate_min, 1, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &high_candidate_min, 1, MPI_DOUBLE,
                          MPI_MIN, MPI_COMM_WORLD);
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
            result.fct_high_candidate_min = std::min(
                result.fct_high_candidate_min, high_candidate_min);
            result.fct_high_candidate_donor_excess = std::max(
                result.fct_high_candidate_donor_excess, invariant_max[8]);

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
                        if (fct_core_x(x))
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
            } // high-order/FCT donor-capacity closure

            const std::vector<double>& cu_high_used = low_order_only_
                ? cu_low : cu_high;
            const std::vector<double>& cu_high_center_used = low_order_only_
                ? cu_low : cu_high_center;
            const std::vector<double>& fx_low_used = fx_low;
            const std::vector<double>& fx_high_used = low_order_only_
                ? fx_low : fx_high;
            const std::vector<double>& fx_final_used = low_order_only_
                ? fx_low : fx_final;
            const std::vector<double>& fu_final_used = low_order_only_
                ? fu_low : fu_final;
            const std::vector<double>& cu_final_used = low_order_only_
                ? cu_low : cu_final;
            std::vector<double> sub_je_low_cell;
            std::vector<double> sub_je_high_cell;
            std::vector<double> sub_je_final_cell;
            std::vector<double> sub_je_low_cell_by_u;
            std::vector<double> sub_je_high_cell_by_u;
            std::vector<double> sub_je_final_cell_by_u;
            if (fixed_candidate) {
                sub_je_low_cell.assign(static_cast<size_t>(nxl), 0.0);
                sub_je_high_cell.assign(static_cast<size_t>(nxl), 0.0);
                sub_je_final_cell.assign(static_cast<size_t>(nxl), 0.0);
                const size_t resolved_size = static_cast<size_t>(nxl) *
                    Param::Nv;
                sub_je_low_cell_by_u.assign(resolved_size, 0.0);
                sub_je_high_cell_by_u.assign(resolved_size, 0.0);
                sub_je_final_cell_by_u.assign(resolved_size, 0.0);
            }

            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                double gamma_low = 0.0;
                double gamma_high = 0.0;
                double gamma_final = 0.0;
                double psi_k = 0.0;
                double psi_p = 0.0;
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k)
                {
                    const size_t id = xface_index(iface, j, k);
                    gamma_low += fx_low_used[id];
                    gamma_high += fx_high_used[id];
                    const double flux_final = fx_final_used[id];
                    gamma_final += flux_final;
                    psi_k += bkg_n.cgrid.kinetic_energy[idx2(j, k)] *
                        flux_final;
                    psi_p += bkg_n.mass * Const::c *
                        bkg_n.cgrid.upar_cells[j] * flux_final;
                }
                integrated_jn_low[iface] += h * bkg_n.charge * gamma_low;
                integrated_jn_high[iface] += h * bkg_n.charge * gamma_high;
                integrated_jn[iface] += h * bkg_n.charge * gamma_final;
                psi_k_x[iface] += h * psi_k;
                psi_p_x[iface] += h * psi_p;
            }
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                double je_low = 0.0;
                double je_high = 0.0;
                double je_center = 0.0;
                double je_final = 0.0;
                double u_energy = 0.0;
                double u_momentum = 0.0;
                for (int jf = 1; jf < Param::Nv; ++jf) for (int k = 0; k < Param::Nmu; ++k) {
                    const double dke = Stage5::delta_energy(bkg_n.cgrid, jf, k);
                    const double center_distance = Stage5::upar_center_distance(
                        bkg_n.cgrid, jf);
                    const size_t id = uface_index(ix, jf, k);
                    const double coefficient_high = cu_high_used[id];
                    const double coefficient_center = cu_high_center_used[id];
                    const double coefficient_final = cu_final_used[id];
                    const double coefficient_low = cu_low[id];
                    const double flux_final = fu_final_used[id];
                    // J_E = q/(m c dx) sum(delta_K * C_u).  The dx division
                    // converts the cell-integrated C_u back to a current
                    // density; h-weighting and /dt occur outside this loop.
                    const double current_factor =
                        bkg_n.charge / (bkg_n.mass * Const::c * sg.dx);
                    je_low += dke * current_factor * coefficient_low;
                    je_high += dke * current_factor * coefficient_high;
                    je_center += dke * current_factor * coefficient_center;
                    je_final += dke * current_factor * coefficient_final;
                    if (fixed_candidate) {
                        const double low_contribution =
                            dke * current_factor * coefficient_low;
                        const double high_contribution =
                            dke * current_factor * coefficient_high;
                        const double final_contribution =
                            dke * current_factor * coefficient_final;
                        const size_t left = static_cast<size_t>(ix) *
                            Param::Nv + static_cast<size_t>(jf - 1);
                        const size_t right = left + 1;
                        sub_je_low_cell_by_u[left] += 0.5 * low_contribution;
                        sub_je_low_cell_by_u[right] += 0.5 * low_contribution;
                        sub_je_high_cell_by_u[left] += 0.5 * high_contribution;
                        sub_je_high_cell_by_u[right] += 0.5 * high_contribution;
                        sub_je_final_cell_by_u[left] += 0.5 * final_contribution;
                        sub_je_final_cell_by_u[right] += 0.5 * final_contribution;
                    }
                    u_energy += dke * flux_final;
                    u_momentum += bkg_n.mass * Const::c *
                        center_distance * flux_final;
                }
                integrated_je_low_cell[ix] += h * je_low;
                integrated_je_cell[ix] += h * je_final;
                integrated_je_high_cell[ix] += h * je_high;
                integrated_je_center_cell[ix] += h * je_center;
                u_energy_cell[ix] += h * u_energy;
                u_momentum_cell[ix] += h * u_momentum;
                if (fixed_candidate) {
                    sub_je_low_cell[static_cast<size_t>(ix)] = je_low;
                    sub_je_high_cell[static_cast<size_t>(ix)] = je_high;
                    sub_je_final_cell[static_cast<size_t>(ix)] = je_final;
                }
            }

            for (int k = 0; k < Param::Nmu; ++k) {
                const double ke_lo = bkg_n.cgrid.kinetic_energy[idx2(0, k)];
                const double ke_hi = bkg_n.cgrid.kinetic_energy[idx2(Param::Nv - 1, k)];
                for (int ix = 0; ix < nxl; ++ix) {
                    const double f_hi = fu_final_used[uface_index(ix, Param::Nv, k)];
                    const double f_lo = fu_final_used[uface_index(ix, 0, k)];
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
            if (fixed_candidate) {
                current_moment_from_x_flux(fx_low_used,
                    substep_audit.jn_low_post_sync);
                current_moment_from_x_flux(fx_high_used,
                    substep_audit.jn_high_post_sync);
                current_moment_from_x_flux(fx_final_used,
                    substep_audit.jn_final_post_sync);
                current_moment_from_x_flux_by_u(fx_low_used,
                    substep_audit.jn_low_by_u_post_sync);
                current_moment_from_x_flux_by_u(fx_high_used,
                    substep_audit.jn_high_by_u_post_sync);
                current_moment_from_x_flux_by_u(fx_final_used,
                    substep_audit.jn_final_by_u_post_sync);
                const int audit_tag = 1200 + 40 * sub;
                PeriodicStaggered::audit_cell_to_face_Gstar_sync(
                    sub_je_low_cell,
                    substep_audit.gstar_je_low_pre_sync,
                    substep_audit.gstar_je_low_post_sync,
                    nxl, mpi_rank, mpi_size, audit_tag);
                PeriodicStaggered::audit_cell_to_face_Gstar_sync(
                    sub_je_high_cell,
                    substep_audit.gstar_je_high_pre_sync,
                    substep_audit.gstar_je_high_post_sync,
                    nxl, mpi_rank, mpi_size, audit_tag + 4);
                PeriodicStaggered::audit_cell_to_face_Gstar_sync(
                    sub_je_final_cell,
                    substep_audit.gstar_je_final_pre_sync,
                    substep_audit.gstar_je_final_post_sync,
                    nxl, mpi_rank, mpi_size, audit_tag + 8);
                PeriodicStaggered::audit_cell_blocks_to_face_Gstar(
                    sub_je_low_cell_by_u,
                    substep_audit.gstar_je_low_by_u_post_sync,
                    nxl, Param::Nv, mpi_rank, mpi_size, audit_tag + 20);
                PeriodicStaggered::audit_cell_blocks_to_face_Gstar(
                    sub_je_high_cell_by_u,
                    substep_audit.gstar_je_high_by_u_post_sync,
                    nxl, Param::Nv, mpi_rank, mpi_size, audit_tag + 24);
                PeriodicStaggered::audit_cell_blocks_to_face_Gstar(
                    sub_je_final_cell_by_u,
                    substep_audit.gstar_je_final_by_u_post_sync,
                    nxl, Param::Nv, mpi_rank, mpi_size, audit_tag + 28);
                result.coupling_substep_seam_audit.push_back(substep_audit);
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

        for (size_t i = 0; i < nface; ++i) {
            result.j_bkg_face_low_mid[i] = integrated_jn_low[i] / dt;
            result.j_bkg_face_high_mid[i] = integrated_jn_high[i] / dt;
            // x transport has one production high-order layer.  Its explicit
            // centered audit view is therefore the same flux, rather than a
            // separately reconstructed current with different ownership.
            result.j_bkg_face_center_mid[i] = result.j_bkg_face_high_mid[i];
            result.j_bkg_face_mid[i] = integrated_jn[i] / dt;
        }
        std::vector<double>& je_cell = je_cell_buffer;
        std::vector<double>& je_low_cell = je_low_cell_buffer;
        for (int ix = 0; ix < nxl; ++ix) {
            je_cell[ix] = integrated_je_cell[ix] / dt;
            je_low_cell[ix] = integrated_je_low_cell[ix] / dt;
        }
        // These are audit-only layers.  They preserve the direct C_u -> J_E
        // construction and never replace the charge-conserving Ampere J_N.
        std::vector<double>& je_high_cell = je_high_cell_buffer;
        std::vector<double>& je_center_cell = je_center_cell_buffer;
        for (int ix = 0; ix < nxl; ++ix) {
            je_high_cell[ix] = integrated_je_high_cell[ix] / dt;
            je_center_cell[ix] = integrated_je_center_cell[ix] / dt;
        }
        result.j_bkg_energy_cell_mid = je_cell;
        if (fixed_candidate) {
            result.j_bkg_energy_low_cell_mid = je_low_cell;
            result.j_bkg_energy_center_cell_mid = je_center_cell;
            result.j_bkg_energy_high_cell_mid = je_high_cell;
        }
        const auto assemble_gstar_face = [&](const std::vector<double>& cell,
                                             std::vector<double>& face,
                                             int message_tag) {
            PeriodicStaggered::apply_cell_to_face_Gstar(
                cell, face, nxl, mpi_rank, mpi_size, message_tag);
        };
        assemble_gstar_face(je_low_cell, result.j_bkg_energy_low_debug_face, 906);
        assemble_gstar_face(je_center_cell, result.j_bkg_energy_center_debug_face, 908);
        assemble_gstar_face(je_high_cell, result.j_bkg_energy_high_debug_face, 910);
        assemble_gstar_face(je_cell, result.j_bkg_energy_debug_face, 912);
        if (fixed_candidate && substeps_used == 1) {
            // Independently re-form every final current moment from the same
            // final FV flux arrays.  This is intentionally an audit-only
            // calculation: J_N remains the production charge current.
            std::vector<double> direct_jn(nface, 0.0);
            std::vector<double> direct_je_cell(static_cast<size_t>(nxl), 0.0);
            std::vector<double> direct_gstar(nface, 0.0);
            const std::vector<double>& audit_fx = low_order_only_ ? fx_low : fx_final;
            const std::vector<double>& audit_cu = low_order_only_ ? cu_low : cu_final;
            for (int iface = 0; iface <= nxl; ++iface) {
                double gamma = 0.0;
                for (int j = 0; j < Param::Nv; ++j)
                    for (int k = 0; k < Param::Nmu; ++k)
                        gamma += audit_fx[xface_index(iface, j, k)];
                direct_jn[static_cast<size_t>(iface)] = bkg_n.charge * gamma;
            }
            for (int ix = 0; ix < nxl; ++ix) {
                double je = 0.0;
                for (int jf = 1; jf < Param::Nv; ++jf)
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const double current_factor = bkg_n.charge /
                            (bkg_n.mass * Const::c * sg.dx);
                        je += Stage5::delta_energy(bkg_n.cgrid, jf, k) *
                            current_factor * audit_cu[uface_index(ix, jf, k)];
                    }
                direct_je_cell[static_cast<size_t>(ix)] = je;
            }
            PeriodicStaggered::apply_cell_to_face_Gstar(
                direct_je_cell, direct_gstar, nxl, mpi_rank, mpi_size, 914);
            double local_audit[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
            for (int iface = 0; iface <= nxl; ++iface) {
                const size_t p = static_cast<size_t>(iface);
                local_audit[0] = std::max(local_audit[0], std::fabs(
                    direct_jn[p] - result.j_bkg_face_mid[p]));
                local_audit[1] = std::max(local_audit[1], std::max(
                    std::fabs(direct_jn[p]), std::fabs(result.j_bkg_face_mid[p])));
                local_audit[4] = std::max(local_audit[4], std::fabs(
                    direct_gstar[p] - result.j_bkg_energy_debug_face[p]));
                local_audit[5] = std::max(local_audit[5], std::max(
                    std::fabs(direct_gstar[p]),
                    std::fabs(result.j_bkg_energy_debug_face[p])));
                if (!std::isfinite(direct_jn[p]) || !std::isfinite(direct_gstar[p]))
                    local_audit[6] = 0.0;
            }
            for (int ix = 0; ix < nxl; ++ix) {
                const size_t p = static_cast<size_t>(ix);
                local_audit[2] = std::max(local_audit[2], std::fabs(
                    direct_je_cell[p] - je_cell[p]));
                local_audit[3] = std::max(local_audit[3], std::max(
                    std::fabs(direct_je_cell[p]), std::fabs(je_cell[p])));
                if (!std::isfinite(direct_je_cell[p])) local_audit[6] = 0.0;
            }
            MPI_Allreduce(MPI_IN_PLACE, local_audit, 6, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, local_audit + 6, 1, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            result.final_flux_current_moment_audit_valid = 1;
            result.final_flux_current_moment_audit_finite =
                local_audit[6] > 0.5 ? 1 : 0;
            result.final_flux_to_jn_linf = local_audit[0];
            result.final_flux_to_jn_scale = local_audit[1];
            result.final_flux_to_je_linf = local_audit[2];
            result.final_flux_to_je_scale = local_audit[3];
            result.final_flux_to_gstar_je_linf = local_audit[4];
            result.final_flux_to_gstar_je_scale = local_audit[5];
        }
        for (int iface = 0; iface < nxl; ++iface) {
            result.j_beam_face_mid[iface] = jbeam[iface];
            result.j_total_face_mid[iface] = result.j_bkg_face_mid[iface] + jbeam[iface];
        }
        if (nxl > 0) {
            result.j_beam_face_mid[nxl] = jbeam[nxl];
            result.j_total_face_mid[nxl] = result.j_total_face_mid[0];
        }
        if (iter == 0)
            low_order_solver_checkpoint(low_order_only_, "currents_ready", mpi_rank);
        // Keep the two representations of the global periodic seam separate
        // for fixed-state replay.  Face 0 is owned by rank 0; face nxl on the
        // final rank is its periodic alias at x = L.  This is audit-only: it
        // neither closes nor modifies either current.
        double seam_local[4] = {0.0, 0.0, 0.0, 0.0};
        if (nxl > 0 && mpi_rank == 0) {
            seam_local[0] = result.j_bkg_face_mid[0];
            seam_local[2] = result.j_bkg_energy_debug_face[0];
        }
        if (nxl > 0 && mpi_rank == mpi_size - 1) {
            seam_local[1] = result.j_bkg_face_mid[nxl];
            seam_local[3] = result.j_bkg_energy_debug_face[nxl];
        }
        double seam_global[4] = {0.0, 0.0, 0.0, 0.0};
        MPI_Allreduce(seam_local, seam_global, 4, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.periodic_seam_face_audit[0] = seam_global[0];
        result.periodic_seam_face_audit[1] = seam_global[1];
        result.periodic_seam_face_audit[2] = seam_global[2];
        result.periodic_seam_face_audit[3] = seam_global[3];
        result.periodic_seam_face_audit[4] = seam_global[0] - seam_global[2];
        result.periodic_seam_face_audit[5] = seam_global[1] - seam_global[3];

        // A fixed operator audit evaluates the FV/current operator at an
        // externally supplied endpoint. It must not perform a second Ampere
        // update or advance/reinject Beam particles.
        EMFields fields_new = (fixed_candidate && fixed_fields_end)
            ? *fixed_fields_end : fields_n;
        if (!fixed_candidate) {
            fields_new.advance_ampere_face_from_midpoint_current(
                result.j_total_face_mid, dt, mpi_rank, mpi_size);
        }
        if (iter == 0)
            low_order_solver_checkpoint(low_order_only_, "ampere_ready", mpi_rank);
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
        if (midpoint_iteration_trace_for_test_) {
            result.midpoint_residual_e_history.push_back(result.residual_E);
            result.midpoint_residual_j_bkg_history.push_back(result.residual_J_bkg);
            result.midpoint_residual_j_beam_history.push_back(result.residual_J_beam);
        }

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
        if (midpoint_iteration_trace_for_test_)
            result.midpoint_residual_f_history.push_back(result.residual_f);
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

        double local_cell_work = 0.0, local_high_cell_work = 0.0;
        double local_center_cell_work = 0.0, local_face_work = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            local_cell_work += fields_mid.Ex[ng + ix] * je_cell[ix] * sg.dx;
            local_high_cell_work += fields_mid.Ex[ng + ix] *
                je_high_cell[ix] * sg.dx;
            local_center_cell_work += fields_mid.Ex[ng + ix] *
                je_center_cell[ix] * sg.dx;
        }
        for (int iface = 0; iface < nxl; ++iface)
            local_face_work += e_mid[iface] * result.j_bkg_face_mid[iface] * sg.dx;
        double energy_pair[4] = {local_cell_work, local_face_work,
                                  local_high_cell_work,
                                  local_center_cell_work};
        MPI_Allreduce(MPI_IN_PLACE, energy_pair, 4, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        result.energy_residual_bkg = dt * (energy_pair[0] - energy_pair[1]);

        double local_stage5[7] = {0.0, 0.0, 0.0, 0.0,
                                  psi_k_x[nxl] - psi_k_x[0],
                                  psi_p_x[nxl] - psi_p_x[0], 0.0};
        for (int ix = 0; ix < nxl; ++ix) {
            local_stage5[0] += local_energy_residual[ix];
            local_stage5[1] += local_momentum_residual[ix];
            local_stage5[2] += u_energy_cell[ix];
            local_stage5[3] += u_momentum_cell[ix];
            local_stage5[6] += delta_ke_cell[ix];
        }
        MPI_Allreduce(MPI_IN_PLACE, local_stage5, 7, MPI_DOUBLE, MPI_SUM,
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
        result.stage5_r_couple_centered = dt *
            (energy_pair[1] - energy_pair[3]);
        result.stage5_r_couple_upwind_stabilization = dt *
            (energy_pair[3] - energy_pair[2]);
        result.stage5_r_couple_fct_stabilization = dt *
            (energy_pair[2] - energy_pair[0]);
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
        double limiter_sum[12] = {limiter_energy_change, u_boundary_energy,
                                  limiter_faces, limiter_active,
                                  limiter_energy_positive, limiter_energy_negative,
                                  limiter_energy_core, limiter_energy_boundary,
                                  limiter_faces_core, limiter_faces_boundary,
                                  limiter_active_core, limiter_active_boundary};
        MPI_Allreduce(MPI_IN_PLACE, limiter_sum, 12, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        double limiter_alpha_global[3] = {
            limiter_min, limiter_min_core, limiter_min_boundary};
        MPI_Allreduce(MPI_IN_PLACE, limiter_alpha_global, 3, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        result.limiter_energy_defect = limiter_sum[0];
        result.limiter_energy_defect_positive = limiter_sum[4];
        result.limiter_energy_defect_negative = limiter_sum[5];
        result.limiter_energy_defect_core = limiter_sum[6];
        result.limiter_energy_defect_boundary = limiter_sum[7];
        result.flux_defect[1].boundary_energy_loss = limiter_sum[1];
        result.limiter_active_fraction = (limiter_sum[2] > 0.0)
            ? limiter_sum[3] / limiter_sum[2] : 0.0;
        result.limiter_active_fraction_core = (limiter_sum[8] > 0.0)
            ? limiter_sum[10] / limiter_sum[8] : 0.0;
        result.limiter_active_fraction_boundary = (limiter_sum[9] > 0.0)
            ? limiter_sum[11] / limiter_sum[9] : 0.0;
        result.limiter_min_alpha = limiter_alpha_global[0];
        result.limiter_min_alpha_core = limiter_alpha_global[1];
        result.limiter_min_alpha_boundary = limiter_alpha_global[2];
        result.flux_pos[0].alpha_active_fraction = result.limiter_active_fraction;
        result.flux_pos[1].alpha_active_fraction = result.limiter_active_fraction;
        result.flux_pos[0].alpha_min = limiter_alpha_global[0];
        result.flux_pos[1].alpha_min = limiter_alpha_global[0];

        result.delta_ke_bkg = local_stage5[6];
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
        const auto finalize_coupling_regions = [&]() {
            if (fixed_coupling_layout) {
                result.coupling_beam_front_ix =
                    fixed_coupling_layout->beam_front_ix;
                result.coupling_wave_core_end_m =
                    fixed_coupling_layout->wave_core_end_m;
            } else {
                const double beam_threshold = 1.0e-6 * Param::densb;
                int local_beam_front = -1;
                for (int ix = 0; ix < nxl; ++ix) {
                    if (ix < static_cast<int>(beam_predictor.density.size()) &&
                        beam_predictor.density[ix] > beam_threshold) {
                        local_beam_front = std::max(local_beam_front,
                            sg.ix_start + ix);
                    }
                }
                MPI_Allreduce(&local_beam_front, &result.coupling_beam_front_ix, 1,
                              MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                const double front_x = result.coupling_beam_front_ix >= 0
                    ? (static_cast<double>(result.coupling_beam_front_ix) + 0.5) * sg.dx
                    : 0.8 * Const::micro;
                result.coupling_wave_core_end_m = std::max(0.8 * Const::micro,
                    std::min(7.2 * Const::micro, front_x + 0.2 * Const::micro));
            }
            const auto region_for_x = [&](double x) {
                if (x < 0.2 * Const::micro) return 0;       // B_left
                if (x < 0.8 * Const::micro) return 1;       // Q_left
                if (x <= result.coupling_wave_core_end_m) return 2; // wave_core
                if (x < 7.2 * Const::micro) return 3;       // quiet_right
                if (x < 7.8 * Const::micro) return 4;       // Q_right
                return 5;                                    // B_right
            };
            // Every entry below is accumulated from rank-local cells/faces,
            // then reduced exactly once after the loops.  Do not use any
            // Result field here: those fields are already global.
            const int values_per_region = 14;
            std::array<double, 6 * values_per_region> sums = {};
            std::array<double, 6> linf_rj = {};
            std::array<double, 6> linf_rk = {};
            std::array<double, 6> max_j_difference = {};
            std::array<int, 6> local_max_face;
            local_max_face.fill(std::numeric_limits<int>::max());
            for (int ix = 0; ix < nxl; ++ix) {
                const int region = region_for_x(sg.x(ng + ix));
                const size_t base = static_cast<size_t>(region * values_per_region);
                // This is the local finite-volume energy residual used to
                // form stage5_r_fv, including spatial and velocity-face
                // boundary bookkeeping.  Reusing it makes the regional sum
                // an exact decomposition of the stage-5 residual.
                const double rk = local_energy_residual[ix];
                sums[base + 3] += rk;
                sums[base + 4] += std::fabs(rk);
                sums[base + 5] += rk * rk;
                sums[base + 7] += delta_ke_cell[ix];
                sums[base + 8] += dt * fields_mid.Ex[ng + ix] *
                    (je_high_cell[ix] - je_cell[ix]) * sg.dx;
                sums[base + 10] += 1.0;
                linf_rk[region] = std::max(linf_rk[region], std::fabs(rk));
            }
            for (int iface = 0; iface < nxl; ++iface) {
                const double xface = sg.x_min +
                    static_cast<double>(sg.ix_start + iface) * sg.dx;
                const int region = region_for_x(xface);
                const size_t base = static_cast<size_t>(region * values_per_region);
                const double j_difference = result.j_bkg_face_mid[iface] -
                    result.j_bkg_energy_debug_face[iface];
                const double rj = dt * e_mid[iface] * j_difference * sg.dx;
                const double rj_centered = dt * e_mid[iface] *
                    (result.j_bkg_face_mid[iface] -
                     result.j_bkg_energy_center_debug_face[iface]) * sg.dx;
                const double rj_upwind = dt * e_mid[iface] *
                    (result.j_bkg_energy_center_debug_face[iface] -
                     result.j_bkg_energy_high_debug_face[iface]) * sg.dx;
                const double rj_fct = dt * e_mid[iface] *
                    (result.j_bkg_energy_high_debug_face[iface] -
                     result.j_bkg_energy_debug_face[iface]) * sg.dx;
                sums[base] += rj;
                sums[base + 1] += std::fabs(rj);
                sums[base + 2] += rj * rj;
                sums[base + 6] += e_mid[iface] *
                    result.j_bkg_face_mid[iface] * sg.dx;
                sums[base + 9] += 1.0;
                sums[base + 11] += rj_centered;
                sums[base + 12] += rj_upwind;
                sums[base + 13] += rj_fct;
                linf_rj[region] = std::max(linf_rj[region], std::fabs(rj));
                const double abs_j_difference = std::fabs(j_difference);
                if (abs_j_difference > max_j_difference[region]) {
                    max_j_difference[region] = abs_j_difference;
                    local_max_face[region] = sg.ix_start + iface;
                }
            }
            std::array<double, 6 * 3> maxima = {};
            for (int region = 0; region < 6; ++region) {
                maxima[3 * region] = linf_rj[region];
                maxima[3 * region + 1] = linf_rk[region];
                maxima[3 * region + 2] = max_j_difference[region];
            }
            MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(sums.size()),
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, maxima.data(), static_cast<int>(maxima.size()),
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            std::array<int, 6> max_face;
            for (int region = 0; region < 6; ++region) {
                if (max_j_difference[region] < maxima[3 * region + 2])
                    local_max_face[region] = std::numeric_limits<int>::max();
                max_face[region] = local_max_face[region];
            }
            MPI_Allreduce(MPI_IN_PLACE, max_face.data(), static_cast<int>(max_face.size()),
                          MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            result.coupling_rj_global_sum = 0.0;
            result.coupling_rk_global_sum = 0.0;
            result.coupling_face_work_jn_global_sum = 0.0;
            result.coupling_rj_centered_global_sum = 0.0;
            result.coupling_rj_upwind_global_sum = 0.0;
            result.coupling_rj_fct_global_sum = 0.0;
            for (int region = 0; region < 6; ++region) {
                const size_t base = static_cast<size_t>(region * values_per_region);
                CouplingRegionDiagnostics& output = result.coupling_regions[region];
                output.sum_rj = sums[base];
                output.sum_abs_rj = sums[base + 1];
                output.sum_sq_rj = sums[base + 2];
                output.linf_rj = maxima[3 * region];
                output.sum_rk = sums[base + 3];
                output.sum_abs_rk = sums[base + 4];
                output.sum_sq_rk = sums[base + 5];
                output.linf_rk = maxima[3 * region + 1];
                output.max_abs_jn_minus_gstar_je = maxima[3 * region + 2];
                output.max_abs_jn_minus_gstar_je_face =
                    max_face[region] == std::numeric_limits<int>::max() ? -1 : max_face[region];
                output.face_work_jn = sums[base + 6];
                output.delta_ke_bkg = sums[base + 7];
                output.fct_work = sums[base + 8];
                output.r_couple_centered = sums[base + 11];
                output.r_couple_upwind_stabilization = sums[base + 12];
                output.r_couple_fct_stabilization = sums[base + 13];
                output.face_count = static_cast<long long>(sums[base + 9]);
                output.cell_count = static_cast<long long>(sums[base + 10]);
                result.coupling_rj_global_sum += output.sum_rj;
                result.coupling_rk_global_sum += output.sum_rk;
                result.coupling_face_work_jn_global_sum += output.face_work_jn;
                result.coupling_rj_centered_global_sum += output.r_couple_centered;
                result.coupling_rj_upwind_global_sum +=
                    output.r_couple_upwind_stabilization;
                result.coupling_rj_fct_global_sum +=
                    output.r_couple_fct_stabilization;
            }
            // Automatic regional reconstruction identities.  stage5_r_* and
            // energy_pair[1] are already global; do not reduce them again.
            result.coupling_rj_reconstruction_error =
                result.coupling_rj_global_sum - result.stage5_r_couple;
            result.coupling_rk_reconstruction_error =
                result.coupling_rk_global_sum - result.stage5_r_fv;
            result.coupling_face_work_jn_reconstruction_error =
                result.coupling_face_work_jn_global_sum - energy_pair[1];
            result.coupling_rj_centered_reconstruction_error =
                result.coupling_rj_centered_global_sum -
                result.stage5_r_couple_centered;
            result.coupling_rj_upwind_reconstruction_error =
                result.coupling_rj_upwind_global_sum -
                result.stage5_r_couple_upwind_stabilization;
            result.coupling_rj_fct_reconstruction_error =
                result.coupling_rj_fct_global_sum -
                result.stage5_r_couple_fct_stabilization;
        };
        if (fixed_candidate) {
            capture_accepted_transport(work);
            result.species_np1 = work;
            result.beam_np1 = beam_n;
            result.fields_np1 = fields_new;
            result.converged = true;
            result.state_advanced = 1;
            // Preserve all production layers for the fixed-state operator
            // bundle.  Normal stepping does not retain these extra copies.
            result.low_x_flux = fx_low;
            result.high_x_flux = low_order_only_ ? fx_low : fx_high;
            result.final_x_flux = low_order_only_ ? fx_low : fx_final;
            result.low_u_flux = fu_low;
            result.center_u_flux = low_order_only_ ? fu_low :
                std::vector<double>(fu_high.size(), 0.0);
            result.center_u_coefficient = low_order_only_ ? cu_low :
                cu_high_center;
            result.center_u_reconstruction_mass = midpoint_state.f;
            if (!low_order_only_) {
                for (size_t p = 0; p < fu_high.size(); ++p) {
                    const double coefficient = cu_high_center[p];
                    const int ix = static_cast<int>(p /
                        (static_cast<size_t>(Param::Nv + 1) * Param::Nmu));
                    const double a = bkg_n.charge * fields_mid.Ex[ng + ix] /
                        (bkg_n.mass * Const::c);
                    result.center_u_flux[p] = a * coefficient;
                }
            }
            result.high_u_flux = low_order_only_ ? fu_low : fu_high;
            result.final_u_flux = low_order_only_ ? fu_low : fu_final;
            if (step_diagnostics_enabled_ || fixed_candidate)
                finalize_coupling_regions();
            return result;
        }
        if (have_previous && result.residual_E < field_tol &&
            result.residual_J_bkg < current_tol &&
            result.beam_continuity_residual < 1.0e-6) {
            if (!transport_safe_step_gate(work, nsub))
                return result;
            capture_accepted_transport(work);
            result.species_np1 = work;
            result.beam_np1 = beam_predictor;
            result.fields_np1 = fields_new;
            result.converged = true;
            result.state_advanced = 1;
            if (step_diagnostics_enabled_ || fixed_candidate)
                finalize_coupling_regions();
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
        if (step_diagnostics_enabled_ || fixed_candidate)
            finalize_coupling_regions();
        }
    }

    result.species_np1 = last_work;
    result.beam_np1 = beam_predictor;
    result.fields_np1 = last_fields;
    result.soft_unconverged = true;
    result.soft_accepted = true;
    result.state_advanced = 1;
    return result;
}
