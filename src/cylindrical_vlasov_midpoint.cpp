#include "vlasov_ampere_midpoint.h"
#include "dual_u_coupling.h"
#include "nonuniform_reconstruction.h"
#include "discrete_moment_operators.h"
#include "periodic_staggered_operators.h"
#include "regularized_face_pairing.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

inline double candidate_roundoff_tolerance(
    double m_old, const double low_transfer[4], double m_low,
    const double raw_antidiffusive_transfer[4],
    const double final_antidiffusive_transfer[4])
{
    // The final candidate inherits roundoff from three local constructions:
    // donor-cell low mass, raw high-low transfers, and alpha/beta-limited
    // final transfers.  Using only the last five terms underestimates the
    // bound after a nearly cancelled limiter reconstruction.
    double local_sum = std::fabs(m_old) + std::fabs(m_low);
    for (int face = 0; face < 4; ++face) {
        local_sum += std::fabs(low_transfer[face]);
        local_sum += std::fabs(raw_antidiffusive_transfer[face]);
        local_sum += std::fabs(final_antidiffusive_transfer[face]);
    }
    return 512.0 * std::numeric_limits<double>::epsilon() * local_sum +
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
    // A negative subnormal is an underflow residue, not a resolvable mass
    // deficit.  The flux-difference arithmetic can accumulate more than the
    // fixed denorm_min allowance before it underflows, so preserve the local
    // tolerance for normal values and canonicalize subnormals explicitly.
    const bool subnormal = std::fpclassify(candidate) == FP_SUBNORMAL;
    if (!subnormal && -candidate > roundoff_tolerance) return false;
    normalized_mass = -candidate;
    candidate = 0.0; // Canonical positive zero, not a global distribution clip.
    return true;
}

inline void allgather_spatial_cells(const std::vector<double>& local,
                                    int nx_global, int mpi_size,
                                    std::vector<double>& global)
{
    const int local_count = static_cast<int>(local.size());
    std::vector<int> counts(static_cast<size_t>(mpi_size), 0);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displacements(static_cast<size_t>(mpi_size), 0);
    int total = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
        displacements[static_cast<size_t>(rank)] = total;
        total += counts[static_cast<size_t>(rank)];
    }
    global.assign(static_cast<size_t>(std::max(0, total)), 0.0);
    MPI_Allgatherv(local.empty() ? 0 : local.data(), local_count, MPI_DOUBLE,
                   global.empty() ? 0 : global.data(), counts.data(),
                   displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    if (total != nx_global)
        global.clear();
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
    bool fixed_candidate,
    const std::vector<double>* initial_e_end) const
{
    low_order_solver_checkpoint(low_order_only_, "entry", mpi_rank);
    Result result;
    reset_result(result);
    const double time_fs = time / Const::femto;
    const bool live_trace =
        progress_trace_start_fs_ >= 0.0 &&
        time_fs >= progress_trace_start_fs_ &&
        time_fs <= progress_trace_end_fs_;
    const double live_trace_start = MPI_Wtime();
    const auto trace_live = [&](const char* stage, int iter, int sub) {
        if (live_trace && mpi_rank == 0) {
            std::printf(
                "[midpoint-live] t_fs=%.16e wall_s=%.6f iter=%d sub=%d "
                "stage=%s injection_active=%d beam_local_rank0=%llu\n",
                time_fs, MPI_Wtime() - live_trace_start, iter, sub, stage,
                time <= Param::t_inject_end ? 1 : 0,
                static_cast<unsigned long long>(beam_n.particles.size()));
            std::fflush(stdout);
        }
    };
    trace_live("operator_entry", -1, -1);
    result.background_coupling_mode =
        static_cast<int>(background_coupling_mode_);
    result.dual_u_operator_valid =
        background_coupling_mode_ == LEGACY_COUPLING ? 1 : 0;
    if (background_coupling_mode_ == DUAL_U_COUPLING && low_order_only_) {
        result.failed = true;
        result.state_advanced = 0;
        result.failure_reason = 13;
        if (mpi_rank == 0) {
            std::cerr << "dual_u_coupling configuration error: "
                      << "requires low_order_only=0; the explicit runtime "
                      << "dual-u experiment otherwise retains the production "
                      << "Beam, Ampere, and FCT paths\n";
        }
        return result;
    }
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
    // Reuse the dominant FV/FCT storage across physical steps as well as
    // nonlinear iterations.  Resizing retains capacity, while fill preserves
    // the exact initialization previously provided by vector(size, value).
    size_t workspace_slot = 0;
    const auto workspace =
        [&](size_t size, double value) -> std::vector<double>& {
            std::vector<double>& buffer =
                production_workspace_[workspace_slot++];
            buffer.resize(size);
            std::fill(buffer.begin(), buffer.end(), value);
            return buffer;
        };
    std::vector<double>& fx_low = workspace(xflux_size, 0.0);
    std::vector<double>& fu_low = workspace(uflux_size, 0.0);
    std::vector<double>& cu_low = workspace(uflux_size, 0.0);
    const size_t high_local_cell_count = low_order_only_
        ? 0 : static_cast<size_t>(nxl) * Param::Nvmu;
    const size_t high_xflux_size = low_order_only_ ? 0 : xflux_size;
    const size_t high_uflux_size = low_order_only_ ? 0 : uflux_size;
    std::vector<double>& fx_high = workspace(high_xflux_size, 0.0);
    std::vector<double>& fu_high = workspace(high_uflux_size, 0.0);
    std::vector<double>& cu_high = workspace(high_uflux_size, 0.0);
    // Centered u-force candidate is kept separately from the selected high
    // flux so the closure audit can isolate boundary/upwind stabilization.
    std::vector<double>& cu_high_center =
        workspace(high_uflux_size, 0.0);
    const bool dual_u_enabled =
        background_coupling_mode_ == DUAL_U_COUPLING;
    // Replaying the x and u reconstruction operators is a fixed-state/audit
    // check.  The rank-one production correction itself uses the original
    // production fluxes directly, so ordinary production iterations do not
    // need a second full phase-space reconstruction.
    const bool dual_u_replay_audit = fixed_candidate ||
        capture_midpoint_input_;
    // Worst-cell metadata is diagnostic-only.  The limiter decisions and
    // transport safety gates use their own reductions below, so ordinary
    // production steps do not need MAXLOC/Bcast traffic for this metadata.
    const bool detailed_operator_diagnostics =
        step_diagnostics_enabled_ || fixed_candidate;
    std::vector<double>& inv_cell_volume =
        workspace(low_order_only_ ? 0 : Param::Nvmu, 0.0);
    if (!low_order_only_) {
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                inv_cell_volume[idx2(j, k)] = 1.0 /
                    (sg.dx * bkg_n.cgrid.cell_phase_volume(j, k));
            }
        }
    }
    std::vector<double>& alpha_x_left =
        workspace(high_local_cell_count, 1.0);
    std::vector<double>& alpha_x_right =
        workspace(high_local_cell_count, 1.0);
    std::vector<double>& alpha_u_lower =
        workspace(high_local_cell_count, 1.0);
    std::vector<double>& alpha_u_upper =
        workspace(high_local_cell_count, 1.0);
    std::vector<double>& left_alpha_x_right =
        workspace(low_order_only_ ? 0 : Param::Nvmu, 1.0);
    std::vector<double>& fx_final = workspace(high_xflux_size, 0.0);
    std::vector<double>& fu_final = workspace(high_uflux_size, 0.0);
    std::vector<double>& cu_final = workspace(high_uflux_size, 0.0);
    std::vector<double>& candidate_mass =
        workspace(high_local_cell_count, 0.0);
    std::vector<double>& donor_beta =
        workspace(high_local_cell_count, 1.0);
    std::vector<double>& donor_low_mass =
        workspace(high_local_cell_count, 0.0);
    std::vector<double>& donor_limited_outflow =
        workspace(high_local_cell_count, 0.0);
    Species low_state_buffer = bkg_n;
    std::vector<double>& e_mid_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& next_e_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& integrated_jn_buffer = workspace(nface, 0.0);
    std::vector<double>& integrated_jn_low_buffer = workspace(nface, 0.0);
    std::vector<double>& integrated_jn_high_buffer = workspace(nface, 0.0);
    std::vector<double>& integrated_je_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& integrated_je_low_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& integrated_je_high_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& integrated_je_center_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& je_cell_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& je_low_cell_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& je_high_cell_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& je_center_cell_buffer =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& initial_ke_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& initial_p_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& final_ke_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& final_p_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& delta_ke_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& delta_p_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& local_energy_residual =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& local_momentum_residual =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& psi_k_x = workspace(nface, 0.0);
    std::vector<double>& psi_p_x = workspace(nface, 0.0);
    // Dual-u scratch is allocated once per outer solve and reused by every
    // midpoint iteration/substep.  Do not store these in Result unless an
    // explicit audit requests them.
    std::vector<double>& dual_jn_high_face =
        workspace(dual_u_enabled ? nface : 0, 0.0);
    std::vector<double>& dual_target_jn_cell = workspace(
        dual_u_enabled ? static_cast<size_t>(nxl) : 0, 0.0);
    std::vector<double>& dual_target_jn_replay_cell = workspace(
        dual_u_enabled ? static_cast<size_t>(nxl) : 0, 0.0);
    std::vector<double>& dual_replay_jn_face =
        workspace(dual_u_replay_audit ? nface : 0, 0.0);
    std::vector<double>& dual_legacy_je_cell = workspace(
        dual_u_enabled ? static_cast<size_t>(nxl) : 0, 0.0);
    std::vector<double>& dual_je_cell = workspace(
        dual_u_enabled ? static_cast<size_t>(nxl) : 0, 0.0);
    std::vector<double>& final_dual_acceleration = workspace(
        dual_u_enabled ? static_cast<size_t>(nxl) : 0, 0.0);
    std::vector<double>& u_energy_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& u_momentum_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& u_boundary_energy_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    std::vector<double>& u_boundary_momentum_cell =
        workspace(static_cast<size_t>(nxl), 0.0);
    // These u_parallel-face factors depend only on the immutable velocity
    // grid and species constants for this solve.  Reusing them removes the
    // same scalar divisions/products from every x cell and force substep.
    std::vector<double>& u_face_energy_delta =
        workspace(Param::Nvmu, 0.0);
    std::vector<double>& u_face_momentum_delta =
        workspace(Param::Nvmu, 0.0);
    std::vector<double>& u_face_energy_current_weight =
        workspace(Param::Nvmu, 0.0);
    // The x part of the combined CFL depends only on the immutable velocity
    // grid and this physical step's dt/dx.  Collapse the u_perp maximum once
    // here instead of rescanning every x-u_parallel-u_perp cell in every
    // nonlinear iteration.  The force CFL below is independent of u_perp,
    // so max_k(cx(j,k) + cu(ix,j)) = max_k(cx(j,k)) + cu(ix,j).
    std::vector<double>& max_x_cfl_by_upar =
        workspace(static_cast<size_t>(Param::Nv), 0.0);
    for (int j = 0; j < Param::Nv; ++j) {
        double maximum = 0.0;
        for (int k = 0; k < Param::Nmu; ++k) {
            maximum = std::max(
                maximum,
                dt * std::fabs(bkg_n.cgrid.vx[idx2(j, k)]) / sg.dx);
        }
        max_x_cfl_by_upar[static_cast<size_t>(j)] = maximum;
    }
    const double energy_current_factor =
        bkg_n.charge / (bkg_n.mass * Const::c * sg.dx);
    for (int jf = 1; jf < Param::Nv; ++jf) {
        const double momentum_delta = bkg_n.mass * Const::c *
            bkg_n.cgrid.upar_center_distances[jf];
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t face = idx2(jf, k);
            const double energy_delta =
                Stage5::delta_energy(bkg_n.cgrid, jf, k);
            u_face_energy_delta[face] = energy_delta;
            u_face_momentum_delta[face] = momentum_delta;
            u_face_energy_current_weight[face] =
                energy_delta * energy_current_factor;
        }
    }
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

    std::vector<double>& e_end =
        workspace(static_cast<size_t>(nxl), 0.0);
    for (int iface = 0; iface < nxl; ++iface) {
        e_end[static_cast<size_t>(iface)] = fixed_fields_end
            ? fixed_fields_end->Ex_face[iface]
            : initial_e_end &&
              initial_e_end->size() == static_cast<size_t>(nxl)
            ? (*initial_e_end)[static_cast<size_t>(iface)]
            : fields_n.Ex_face[iface];
    }
    Species guess = fixed_guess_np1 ? *fixed_guess_np1 : bkg_n;
    Species last_work = bkg_n;
    EMFields last_fields = fields_n;
    std::vector<double>& previous_j =
        workspace(static_cast<size_t>(nxl), 0.0);
    bool have_previous = false;

    // Acceleration is deliberately field-only.  The distribution guess keeps
    // the regular Picard update, so an attempted field extrapolation can be
    // rejected and re-evaluated with the exact same input Species state.
    const bool midpoint_acceleration_enabled =
        !fixed_candidate &&
        midpoint_acceleration_mode_ != MIDPOINT_ACCELERATION_NONE;
    result.midpoint_acceleration_mode = midpoint_acceleration_enabled
        ? static_cast<int>(midpoint_acceleration_mode_)
        : MIDPOINT_ACCELERATION_NONE;
    bool acceleration_trial_pending = false;
    bool fallback_plain_evaluation = false;
    double acceleration_residual_before = 0.0;
    double acceleration_trial_coefficient = 0.55;
    double aitken_previous_omega = 0.55;
    bool have_acceleration_residual = false;
    int anderson_history_count = 0;
    if (midpoint_acceleration_enabled) {
        const size_t face_count = static_cast<size_t>(nxl);
        const auto prepare_acceleration_vector =
            [face_count](std::vector<double>& values) {
                if (values.size() != face_count) values.assign(face_count, 0.0);
            };
        prepare_acceleration_vector(acceleration_re_previous_);
        prepare_acceleration_vector(acceleration_re_current_);
        prepare_acceleration_vector(acceleration_fallback_e_);
        for (int h = 0; h < 3; ++h) {
            prepare_acceleration_vector(acceleration_e_history_[h]);
            prepare_acceleration_vector(acceleration_re_history_[h]);
        }
    }
    const auto reject_accelerated_trial_before_residual = [&]() {
        if (!acceleration_trial_pending) return false;
        ++result.acceleration_rejected_hard_failure;
        ++result.acceleration_fallback_evaluations;
        ++result.acceleration_history_resets;
        if (midpoint_iteration_trace_for_test_ &&
            !result.midpoint_acceleration_status_history.empty()) {
            result.midpoint_acceleration_status_history.back() = 4;
        }
        e_end = acceleration_fallback_e_;
        acceleration_trial_pending = false;
        fallback_plain_evaluation = true;
        have_acceleration_residual = false;
        anderson_history_count = 0;
        aitken_previous_omega = 0.55;
        result.failed = false;
        result.state_advanced = 0;
        result.failure_reason = 0;
        return true;
    };

    // Stage 7 is deliberately not enabled: the beam is a controlled
    // predictor during phases 1--5, and its lag is kept explicit.
    trace_live("beam_copy_begin", -1, -1);
    BeamPIC beam_predictor = beam_n;
    trace_live("beam_copy_end", -1, -1);
    EMFields beam_field = fields_n;
    double beam_ke_before = 0.0;
    double beam_ke_after = 0.0;
    if (beam_enabled_ && !fixed_j_beam_face_mid) {
        trace_live("beam_begin_step_begin", -1, -1);
        beam_predictor.begin_step(sg, dt);
        trace_live("beam_inject_begin", -1, -1);
        beam_predictor.inject(sg, beam_field, dt, time, mpi_rank, mpi_size);
        trace_live("beam_inject_end", -1, -1);
        beam_ke_before = beam_predictor.total_kinetic_energy();
        trace_live("beam_push_begin", -1, -1);
        beam_predictor.push(sg, beam_field, dt, mpi_rank, mpi_size);
        trace_live("beam_push_end", -1, -1);
        beam_ke_after = beam_predictor.total_kinetic_energy();
        trace_live("beam_density_begin", -1, -1);
        beam_predictor.deposit_density(sg, mpi_rank, mpi_size);
        trace_live("beam_density_end", -1, -1);
        trace_live("beam_current_begin", -1, -1);
        beam_predictor.finalize_charge_conserving_current(sg, dt, mpi_rank,
                                                          mpi_size);
        trace_live("beam_current_end", -1, -1);
    }
    low_order_solver_checkpoint(low_order_only_, "beam_predictor_ready",
                                mpi_rank);
    std::vector<double>& jbeam = workspace(nface, 0.0);
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

    // Beam boundary diagnostics are assembled from an open particle topology.
    // In particular, the right-boundary contribution is owned by the last MPI
    // rank, so the raw BeamPIC diagnostic is not guaranteed to be identical on
    // every rank.  The value participates in the collective control-flow
    // decision below and therefore must be globalized before the Picard loop.
    const double local_beam_continuity_residual = std::max(
        beam_predictor.last_continuity_linf_error(),
        beam_predictor.last_boundary_flux_error());
    const double local_beam_boundary_source_residual =
        beam_predictor.last_boundary_source_error();
    const double local_beam_open_face_residual =
        beam_predictor.last_open_face_error();
    double global_beam_continuity_residual = local_beam_continuity_residual;
    MPI_Allreduce(MPI_IN_PLACE, &global_beam_continuity_residual, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double beam_boundary_components[2] = {
        local_beam_boundary_source_residual,
        local_beam_open_face_residual};
    MPI_Allreduce(MPI_IN_PLACE, beam_boundary_components, 2, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    int beam_diagnostics_finite =
        std::isfinite(local_beam_continuity_residual) &&
        std::isfinite(local_beam_boundary_source_residual) &&
        std::isfinite(local_beam_open_face_residual) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &beam_diagnostics_finite, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    result.beam_boundary_source_residual = beam_boundary_components[0];
    result.beam_open_face_residual = beam_boundary_components[1];
    result.beam_continuity_valid =
        beam_diagnostics_finite &&
        std::isfinite(global_beam_continuity_residual) &&
        global_beam_continuity_residual <= 1.0e-6 ? 1 : 0;
    if (live_trace) {
        double min_beam_continuity_residual = local_beam_continuity_residual;
        MPI_Allreduce(MPI_IN_PLACE, &min_beam_continuity_residual, 1,
                      MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            std::printf(
                "[midpoint-live] t_fs=%.16e wall_s=%.6f iter=-1 sub=-1 "
                "stage=beam_continuity_globalized local_min=%.16e "
                "global_max=%.16e source=%.16e open_face=%.16e valid=%d\n",
                time_fs, MPI_Wtime() - live_trace_start,
                min_beam_continuity_residual,
                global_beam_continuity_residual,
                result.beam_boundary_source_residual,
                result.beam_open_face_residual,
                result.beam_continuity_valid);
            std::fflush(stdout);
        }
    }
    if (!beam_diagnostics_finite ||
        !std::isfinite(global_beam_continuity_residual)) {
        result.failed = true;
        result.state_advanced = 0;
        result.failure_reason = 17;
        return result;
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
        result.midpoint_acceleration_omega_history.reserve(
            static_cast<size_t>(max_iters));
        result.midpoint_acceleration_residual_before_history.reserve(
            static_cast<size_t>(max_iters));
        result.midpoint_acceleration_status_history.reserve(
            static_cast<size_t>(max_iters));
    }
    const double field_tol = 1.0e-6;
    const double current_tol = 1.0e-5;
    const double omega = 0.55;
    enum AcceptedEnergyLedgerSlot {
        LEDGER_DKE_LOW = 0,
        LEDGER_DKE_HIGH,
        LEDGER_DKE_FINAL,
        LEDGER_DKE_FCT_X,
        LEDGER_DKE_FCT_U,
        LEDGER_DMASS_FCT_X,
        LEDGER_DMASS_FCT_U,
        LEDGER_DPPAR_FCT_X,
        LEDGER_DPPAR_FCT_U,
        LEDGER_BNUM_LOWER,
        LEDGER_BNUM_UPPER,
        LEDGER_BPPAR_LOWER,
        LEDGER_BPPAR_UPPER,
        LEDGER_BENERGY_LOWER,
        LEDGER_BENERGY_UPPER,
        LEDGER_SLOT_COUNT
    };
    std::array<long double, LEDGER_SLOT_COUNT> accepted_energy_ledger_local = {};
    long double accepted_energy_actual_delta_ke_local = 0.0L;
    for (int iter = 0; iter < max_iters; ++iter) {
midpoint_iteration_retry:
        ++result.operator_evaluations;
        trace_live("iteration_begin", iter + 1, -1);
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
                local_cfl = std::max(
                    local_cfl,
                    max_x_cfl_by_upar[static_cast<size_t>(j)] + cu);
            }
        }
        double global_cfl = local_cfl;
        MPI_Allreduce(MPI_IN_PLACE, &global_cfl, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        const int nsub = std::max(1, static_cast<int>(std::ceil(global_cfl / 0.8)));
        if (live_trace && mpi_rank == 0) {
            std::printf(
                "[midpoint-live] t_fs=%.16e wall_s=%.6f iter=%d sub=-1 "
                "stage=cfl_ready global_cfl=%.16e nsub=%d\n",
                time_fs, MPI_Wtime() - live_trace_start, iter + 1,
                global_cfl, nsub);
            std::fflush(stdout);
        }
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
            if (reject_accelerated_trial_before_residual())
                goto midpoint_iteration_retry;
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
        double u_boundary_energy_lower = 0.0;
        double u_boundary_energy_upper = 0.0;
        double u_boundary_particle = 0.0;
        double u_boundary_momentum = 0.0;
        std::array<long double, LEDGER_SLOT_COUNT> trial_energy_ledger_local = {};
        long double trial_actual_delta_ke_local = 0.0L;
        double limiter_faces = 0.0;
        double limiter_active = 0.0;
        double limiter_min = 1.0;
        double x_limiter_faces = 0.0;
        double x_limiter_active = 0.0;
        double x_limiter_min = 1.0;
        double u_limiter_faces = 0.0;
        double u_limiter_active = 0.0;
        double u_limiter_min = 1.0;
        double limiter_faces_core = 0.0;
        double limiter_faces_boundary = 0.0;
        double limiter_active_core = 0.0;
        double limiter_active_boundary = 0.0;
        double limiter_min_core = 1.0;
        double limiter_min_boundary = 1.0;
        double fct_budget_violation = 0.0;
        DualUCoupling::FinalLimitedDiagnostics
            final_dual_iteration_diagnostics;

        const double h = dt / nsub;
        const bool fct_macro_budget_enabled = step_diagnostics_enabled_ &&
            fct_enabled_ && !low_order_only_;
        std::array<FctMacroBudget, 6> fct_macro_budget_x_local;
        std::array<FctMacroBudget, 6> fct_macro_budget_u_local;
        const auto reset_fct_macro_budget = [](std::array<FctMacroBudget, 6>& budget) {
            for (size_t bin = 0; bin < budget.size(); ++bin) {
                budget[bin].face_count = 0;
                budget[bin].active_face_count = 0;
                budget[bin].min_alpha = 1.0;
                budget[bin].delta_n = 0.0;
                budget[bin].delta_j = 0.0;
                budget[bin].delta_k = 0.0;
                budget[bin].e_dot_j = 0.0;
                budget[bin].r_fct_e = 0.0;
            }
        };
        reset_fct_macro_budget(fct_macro_budget_x_local);
        reset_fct_macro_budget(fct_macro_budget_u_local);
        std::vector<double> fct_macro_fmax;
        if (fct_macro_budget_enabled) {
            fct_macro_fmax.assign(static_cast<size_t>(sg.nx_total), 0.0);
            for (int cell = 0; cell < sg.nx_total; ++cell) {
                double maximum = 0.0;
                for (int j = 0; j < Param::Nv; ++j)
                    for (int k = 0; k < Param::Nmu; ++k)
                        maximum = std::max(maximum,
                            midpoint_state.f[mass_index(cell, j, k)]);
                fct_macro_fmax[static_cast<size_t>(cell)] = maximum;
            }
        }
        const auto fct_macro_bin = [&](int cell, int j, int k) {
            int global_ix = sg.ix_start + cell - ng;
            global_ix %= sg.nx_global;
            if (global_ix < 0) global_ix += sg.nx_global;
            const double x = sg.x_min + (static_cast<double>(global_ix) + 0.5) * sg.dx;
            const double x_max = sg.x_min + sg.nx_global * sg.dx;
            const int x_region = x < sg.x_min + 0.2e-6 ? 0 :
                (x > x_max - 0.2e-6 ? 2 : 1);
            const double fmax = fct_macro_fmax[static_cast<size_t>(cell)];
            const bool velocity_core = fmax > 0.0 &&
                midpoint_state.f[mass_index(cell, j, k)] >= 1.0e-8 * fmax;
            return 2 * x_region + (velocity_core ? 0 : 1);
        };
        const auto add_fct_macro_cell = [&](std::array<FctMacroBudget, 6>& budgets,
                                            int cell, int j, int k,
                                            double delta_mass) {
            const int bin = fct_macro_bin(cell, j, k);
            FctMacroBudget& budget = budgets[
                static_cast<size_t>(bin)];
            budget.delta_n += delta_mass;
            budget.delta_j += bkg_n.charge * bkg_n.cgrid.vx[idx2(j, k)] *
                delta_mass;
            budget.delta_k += bkg_n.cgrid.kinetic_energy[idx2(j, k)] *
                delta_mass;
        };
        const auto accumulate_final_fct_macro_budget = [&]() {
            if (!fct_macro_budget_enabled) return;
            // Each owned x face is visited once.  Its conservative update is
            // split between its two neighbouring cells for region/velocity
            // attribution, while the global sum remains exact.
            for (int iface = 0; iface < nxl; ++iface) {
                const int left_cell = ng + iface - 1;
                const int right_cell = ng + iface;
                for (int j = 0; j < Param::Nv; ++j) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = xface_index(iface, j, k);
                        const double anti = fx_high[id] - fx_low[id];
                        if (anti == 0.0) continue;
                        const double alpha = (fx_final[id] - fx_low[id]) / anti;
                        const double delta_flux = fx_final[id] - fx_high[id];
                        const bool active = alpha < 1.0 - 1.0e-14;
                        const int face_bin = fct_macro_bin(right_cell, j, k);
                        FctMacroBudget& face_budget = fct_macro_budget_x_local[
                            static_cast<size_t>(face_bin)];
                        ++face_budget.face_count;
                        if (active) ++face_budget.active_face_count;
                        face_budget.min_alpha = std::min(face_budget.min_alpha,
                                                         alpha);
                        add_fct_macro_cell(fct_macro_budget_x_local,
                                           left_cell, j, k, -h * delta_flux);
                        add_fct_macro_cell(fct_macro_budget_x_local,
                                           right_cell, j, k, h * delta_flux);
                        const double face_work = h * sg.dx *
                            fields_mid.Ex_face[iface] * bkg_n.charge * delta_flux;
                        fct_macro_budget_x_local[static_cast<size_t>(
                            fct_macro_bin(left_cell, j, k))].e_dot_j +=
                            0.5 * face_work;
                        fct_macro_budget_x_local[static_cast<size_t>(
                            fct_macro_bin(right_cell, j, k))].e_dot_j +=
                            0.5 * face_work;
                    }
                }
            }
            for (int ix = 0; ix < nxl; ++ix) {
                const int cell = ng + ix;
                for (int jf = 1; jf < Param::Nv; ++jf) {
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = uface_index(ix, jf, k);
                        const double anti = fu_high[id] - fu_low[id];
                        if (anti == 0.0) continue;
                        const double alpha = (fu_final[id] - fu_low[id]) / anti;
                        const double delta_flux = fu_final[id] - fu_high[id];
                        FctMacroBudget& face_budget = fct_macro_budget_u_local[
                            static_cast<size_t>(fct_macro_bin(cell, jf - 1, k))];
                        ++face_budget.face_count;
                        if (alpha < 1.0 - 1.0e-14)
                            ++face_budget.active_face_count;
                        face_budget.min_alpha = std::min(face_budget.min_alpha,
                                                         alpha);
                        add_fct_macro_cell(fct_macro_budget_u_local,
                                           cell, jf - 1, k, -h * delta_flux);
                        add_fct_macro_cell(fct_macro_budget_u_local,
                                           cell, jf, k, h * delta_flux);
                        const double delta_je = u_face_energy_current_weight[
                            idx2(jf, k)] * (cu_final[id] - cu_high[id]);
                        const double face_work = h * sg.dx *
                            fields_mid.Ex[cell] * delta_je;
                        fct_macro_budget_u_local[static_cast<size_t>(
                            fct_macro_bin(cell, jf - 1, k))].e_dot_j +=
                            0.5 * face_work;
                        fct_macro_budget_u_local[static_cast<size_t>(
                            fct_macro_bin(cell, jf, k))].e_dot_j +=
                            0.5 * face_work;
                    }
                }
            }
        };
        const auto finalize_fct_macro_budget = [&]() {
            if (!fct_macro_budget_enabled) return;
            long long counts[24] = {0};
            double sums[48] = {0.0};
            double minima[12] = {
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            std::array<FctMacroBudget, 6>* local_budgets[2] = {
                &fct_macro_budget_x_local, &fct_macro_budget_u_local};
            for (int direction = 0; direction < 2; ++direction) {
                for (int bin = 0; bin < 6; ++bin) {
                    const int offset = 6 * direction + bin;
                    const FctMacroBudget& local = (*local_budgets[direction])[bin];
                    counts[2 * offset] = local.face_count;
                    counts[2 * offset + 1] = local.active_face_count;
                    sums[4 * offset] = local.delta_n;
                    sums[4 * offset + 1] = local.delta_j;
                    sums[4 * offset + 2] = local.delta_k;
                    sums[4 * offset + 3] = local.e_dot_j;
                    minima[offset] = local.min_alpha;
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, counts, 24, MPI_LONG_LONG_INT,
                          MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, sums, 48, MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, minima, 12, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            result.fct_macro_budget_valid = 1;
            std::array<FctMacroBudget, 6>* result_budgets[2] = {
                &result.fct_macro_budget_x, &result.fct_macro_budget_u};
            for (int direction = 0; direction < 2; ++direction) {
                for (int bin = 0; bin < 6; ++bin) {
                    const int offset = 6 * direction + bin;
                    FctMacroBudget& global = (*result_budgets[direction])[bin];
                    global.face_count = counts[2 * offset];
                    global.active_face_count = counts[2 * offset + 1];
                    global.min_alpha = minima[offset];
                    global.delta_n = sums[4 * offset];
                    global.delta_j = sums[4 * offset + 1];
                    global.delta_k = sums[4 * offset + 2];
                    global.e_dot_j = sums[4 * offset + 3];
                    global.r_fct_e = global.delta_k - global.e_dot_j;
                }
            }
        };
        const auto finalize_directional_limiter_statistics = [&]() {
            if (low_order_only_) {
                result.x_limiter_active_fraction = 0.0;
                result.x_limiter_min_alpha = 1.0;
                result.u_limiter_active_fraction = 0.0;
                result.u_limiter_min_alpha = 1.0;
                return;
            }
            double x_sum[2] = {x_limiter_faces, x_limiter_active};
            MPI_Allreduce(MPI_IN_PLACE, x_sum, 2, MPI_DOUBLE,
                          MPI_SUM, MPI_COMM_WORLD);
            double x_min = x_limiter_min;
            MPI_Allreduce(MPI_IN_PLACE, &x_min, 1, MPI_DOUBLE,
                          MPI_MIN, MPI_COMM_WORLD);
            result.x_limiter_active_fraction = x_sum[0] > 0.0
                ? x_sum[1] / x_sum[0] : 0.0;
            result.x_limiter_min_alpha = x_min;
            if (fixed_candidate || step_diagnostics_enabled_) {
                double u_sum[2] = {u_limiter_faces, u_limiter_active};
                MPI_Allreduce(MPI_IN_PLACE, u_sum, 2, MPI_DOUBLE,
                              MPI_SUM, MPI_COMM_WORLD);
                double u_min = u_limiter_min;
                MPI_Allreduce(MPI_IN_PLACE, &u_min, 1, MPI_DOUBLE,
                              MPI_MIN, MPI_COMM_WORLD);
                result.u_limiter_active_fraction = u_sum[0] > 0.0
                    ? u_sum[1] / u_sum[0] : 0.0;
                result.u_limiter_min_alpha = u_min;
            } else {
                result.u_limiter_active_fraction = 0.0;
                result.u_limiter_min_alpha = 1.0;
            }
        };
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
                // A final copy or ghost synchronization may preserve a signed
                // subnormal that was already below the FCT candidate roundoff
                // scale.  Canonicalize only this underflow residue on the
                // conservative low/FCT paths; ordinary finite negative mass
                // remains subject to the existing hard gate below.
                double local_subnormal_zeroed_count = 0.0;
                double local_subnormal_zeroed_mass = 0.0;
                if (low_order_only_ || fct_enabled_) {
                    #pragma omp parallel for reduction(+:local_subnormal_zeroed_count,local_subnormal_zeroed_mass) schedule(static)
                    for (int ix = 0; ix < nxl; ++ix)
                        for (int j = 0; j < Param::Nv; ++j)
                            for (int k = 0; k < Param::Nmu; ++k) {
                                double& m = candidate_state.f[
                                    mass_index(ng + ix, j, k)];
                                if (m < 0.0 && std::isfinite(m) &&
                                    std::fpclassify(m) == FP_SUBNORMAL) {
                                    local_subnormal_zeroed_count += 1.0;
                                    local_subnormal_zeroed_mass += -m;
                                    m = 0.0;
                                }
                            }
                    double subnormal_zeroed[2] = {
                        local_subnormal_zeroed_count,
                        local_subnormal_zeroed_mass};
                    MPI_Allreduce(MPI_IN_PLACE, subnormal_zeroed, 2,
                                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    result.fct_roundoff_zeroed_count +=
                        static_cast<long long>(subnormal_zeroed[0]);
                    result.fct_roundoff_zeroed_mass += subnormal_zeroed[1];
                }
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
            trace_live("transport_substep_begin", iter + 1, sub + 1);
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
            trace_live("low_fluxes_ready", iter + 1, sub + 1);

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

            if (dual_u_enabled) {
                // The correction uses the production x-current functional
                // directly.  Reconstruction replay is intentionally an
                // audit-only operation; it is not needed to advance dual-u.
                current_moment_from_x_flux(fx_high, dual_jn_high_face);
                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix) {
                    const double target = 0.5 * (
                        dual_jn_high_face[static_cast<size_t>(ix)] +
                        dual_jn_high_face[static_cast<size_t>(ix + 1)]);
                    dual_target_jn_cell[static_cast<size_t>(ix)] = target;
                    // In production this aliases the direct target.  The
                    // independently reconstructed value is filled below only
                    // for fixed-state or diagnostic-level-2 audit runs.
                    dual_target_jn_replay_cell[static_cast<size_t>(ix)] = target;
                }
                double u_replay_linf = 0.0;
                double u_replay_scale = 0.0;

                if (dual_u_replay_audit) {
                const auto replay_x_face_current = [&](int iface) {
                    const int il = ng + iface - 1;
                    const int ir = ng + iface;
                    const double s_face = sg.x_min +
                        (sg.ix_start + iface) * sg.dx;
                    double current = 0.0;
                    for (int j = 0; j < Param::Nv; ++j) {
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double speed = bkg_n.cgrid.vx[idx2(j, k)];
                            NonuniformMuscl::FrozenCenteredLinearization frozen;
                            if (speed == 0.0) {
                                frozen = NonuniformMuscl::frozen_centered_linearization(
                                    fbar(il - 1, j, k), fbar(il, j, k),
                                    fbar(ir, j, k), fbar(ir + 1, j, k),
                                    x_center(il - 1), x_center(il),
                                    x_center(ir), x_center(ir + 1), s_face);
                            } else {
                                frozen = NonuniformMuscl::frozen_upwind_linearization(
                                    fbar(il - 1, j, k), fbar(il, j, k),
                                    fbar(ir, j, k), fbar(ir + 1, j, k),
                                    x_center(il - 1), x_center(il),
                                    x_center(ir), x_center(ir + 1), s_face,
                                    speed > 0.0);
                            }
                            const double q[4] = {
                                fbar(il - 1, j, k), fbar(il, j, k),
                                fbar(ir, j, k), fbar(ir + 1, j, k)};
                            double face_state = 0.0;
                            for (int p = 0; p < 4; ++p)
                                face_state += frozen.coefficient[p] * q[p];
                            current += bkg_n.charge * speed * face_state *
                                bkg_n.cgrid.upar_widths[j] *
                                bkg_n.cgrid.uperp_ring_areas[k];
                        }
                    }
                    return current;
                };

                #pragma omp parallel for schedule(static)
                for (int iface = 0; iface <= nxl; ++iface)
                    dual_replay_jn_face[static_cast<size_t>(iface)] =
                        replay_x_face_current(iface);
                #pragma omp parallel for reduction(max:u_replay_linf,u_replay_scale) schedule(static)
                for (int ix = 0; ix < nxl; ++ix) {
                    dual_target_jn_replay_cell[static_cast<size_t>(ix)] =
                        0.5 * (dual_replay_jn_face[static_cast<size_t>(ix)] +
                               dual_replay_jn_face[static_cast<size_t>(ix + 1)]);
                    const int storage_ix = ng + ix;
                    for (int jf = 1; jf < Param::Nv; ++jf) {
                        const int jl = jf - 1;
                        const int jr = jf;
                        const int jll = jl == 0 ? jl : jl - 1;
                        const int jrr = jr + 1 == Param::Nv ? jr : jr + 1;
                        const double s_jll = jl == 0
                            ? bkg_n.cgrid.upar_cells[jl] -
                              bkg_n.cgrid.upar_widths[jl]
                            : bkg_n.cgrid.upar_cells[jll];
                        const double s_jrr = jr + 1 == Param::Nv
                            ? bkg_n.cgrid.upar_cells[jr] +
                              bkg_n.cgrid.upar_widths[jr]
                            : bkg_n.cgrid.upar_cells[jrr];
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const NonuniformMuscl::FrozenCenteredLinearization
                                frozen = NonuniformMuscl::frozen_centered_linearization(
                                    fbar(storage_ix, jll, k),
                                    fbar(storage_ix, jl, k),
                                    fbar(storage_ix, jr, k),
                                    fbar(storage_ix, jrr, k), s_jll,
                                    bkg_n.cgrid.upar_cells[jl],
                                    bkg_n.cgrid.upar_cells[jr], s_jrr,
                                    bkg_n.cgrid.upar_faces[jf]);
                            const double q[4] = {
                                fbar(storage_ix, jll, k),
                                fbar(storage_ix, jl, k),
                                fbar(storage_ix, jr, k),
                                fbar(storage_ix, jrr, k)};
                            double state = 0.0;
                            for (int p = 0; p < 4; ++p)
                                state += frozen.coefficient[p] * q[p];
                            const double replay =
                                NonuniformMuscl::upar_face_coefficient(state,
                                    sg.dx,
                                    bkg_n.cgrid.uperp_ring_areas[k]);
                            const size_t id = uface_index(ix, jf, k);
                            u_replay_linf = std::max(u_replay_linf,
                                std::fabs(replay - cu_high_center[id]));
                            u_replay_scale = std::max(u_replay_scale,
                                std::max(std::fabs(replay),
                                         std::fabs(cu_high_center[id])));
                        }
                    }
                }
                } // dual_u_replay_audit

                DualUCoupling::Diagnostics dual_diagnostics =
                    DualUCoupling::apply_local_rank_one(
                        bkg_n.cgrid, bkg_n.charge, bkg_n.mass, sg.dx, nxl,
                        dual_target_jn_cell,
                        dual_target_jn_replay_cell,
                        u_face_energy_current_weight,
                        cu_high_center, cu_high,
                        dual_legacy_je_cell, dual_je_cell);
                // A local dual failure must become one collective decision.
                // Returning on only the affected rank would leave the other
                // ranks in later collectives and eventually hang the job.
                const int local_failed_rank = dual_diagnostics.valid
                    ? std::numeric_limits<int>::max() : mpi_rank;
                int first_failed_rank = local_failed_rank;
                MPI_Allreduce(MPI_IN_PLACE, &first_failed_rank, 1, MPI_INT,
                              MPI_MIN, MPI_COMM_WORLD);
                int failure_subtype = DualUCoupling::FAILURE_NONE;
                int failure_global_ix = -1;
                double failure_record[8] = {0.0};
                if (first_failed_rank != std::numeric_limits<int>::max()) {
                    if (mpi_rank == first_failed_rank) {
                        failure_subtype = dual_diagnostics.failure_subtype;
                        failure_global_ix = dual_diagnostics.failure_local_ix < 0
                            ? -1 : sg.ix_start +
                                dual_diagnostics.failure_local_ix;
                        failure_record[0] = dual_diagnostics.failure_target;
                        failure_record[1] = dual_diagnostics.failure_replay;
                        failure_record[2] = dual_diagnostics.failure_legacy;
                        failure_record[3] = dual_diagnostics.failure_residual;
                        failure_record[4] = dual_diagnostics.failure_denominator;
                        failure_record[5] =
                            dual_diagnostics.failure_maximum_coefficient;
                        failure_record[6] =
                            dual_diagnostics.failure_support_floor;
                        failure_record[7] = dual_diagnostics.failure_scale;
                    }
                    MPI_Bcast(&failure_subtype, 1, MPI_INT, first_failed_rank,
                              MPI_COMM_WORLD);
                    MPI_Bcast(&failure_global_ix, 1, MPI_INT,
                              first_failed_rank, MPI_COMM_WORLD);
                    MPI_Bcast(failure_record, 8, MPI_DOUBLE,
                              first_failed_rank, MPI_COMM_WORLD);
                }
                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix) {
                    const double acceleration = bkg_n.charge *
                        fields_mid.Ex[ng + ix] / (bkg_n.mass * Const::c);
                    for (int jf = 0; jf <= Param::Nv; ++jf)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t id = uface_index(ix, jf, k);
                            fu_high[id] = acceleration * cu_high[id];
                        }
                }
                if (dual_u_replay_audit) {
                    result.legacy_center_u_coefficient = cu_high_center;
                    result.dual_target_jn_cell = dual_target_jn_cell;
                    result.dual_target_jn_replay_cell =
                        dual_target_jn_replay_cell;
                    result.dual_legacy_je_cell = dual_legacy_je_cell;
                    result.dual_je_cell = dual_je_cell;
                }
                result.dual_u_operator_valid =
                    first_failed_rank == std::numeric_limits<int>::max();
                result.dual_u_target_replay_linf =
                    dual_diagnostics.target_replay_linf;
                result.dual_u_target_replay_scale =
                    dual_diagnostics.target_replay_scale;
                result.dual_u_legacy_operator_replay_linf = u_replay_linf;
                result.dual_u_legacy_operator_replay_scale = u_replay_scale;
                result.dual_u_legacy_current_linf =
                    dual_diagnostics.legacy_current_linf;
                result.dual_u_current_linf =
                    dual_diagnostics.dual_current_linf;
                result.dual_u_correction_l2 =
                    dual_diagnostics.correction_l2;
                result.dual_u_correction_linf =
                    dual_diagnostics.correction_linf;
                result.dual_u_corrected_cell_count =
                    dual_diagnostics.corrected_cell_count;
                if (!result.dual_u_operator_valid) {
                    result.failed = true;
                    result.state_advanced = 0;
                    result.failure_reason = failure_subtype ==
                        DualUCoupling::FAILURE_GRAM_DEGENERATE ? 14 :
                        (failure_subtype ==
                         DualUCoupling::FAILURE_INPUT_CONTRACT ? 16 : 15);
                    result.failure_iteration = iter;
                    result.failure_substep = sub;
                    if (mpi_rank == 0) {
                        const char* subtype_name =
                            failure_subtype ==
                                DualUCoupling::FAILURE_GRAM_DEGENERATE
                            ? "gram_degenerate"
                            : (failure_subtype ==
                                DualUCoupling::FAILURE_CORRECTION_NONFINITE
                               ? "correction_nonfinite"
                               : (failure_subtype ==
                                    DualUCoupling::FAILURE_DUAL_CURRENT_NONFINITE
                                  ? "dual_current_nonfinite"
                                  : "input_contract"));
                        std::cerr << std::scientific
                                  << "dual_u_coupling_failure"
                                  << " subtype=" << subtype_name
                                  << " rank=" << first_failed_rank
                                  << " global_ix=" << failure_global_ix
                                  << " target=" << failure_record[0]
                                  << " replay=" << failure_record[1]
                                  << " legacy=" << failure_record[2]
                                  << " residual=" << failure_record[3]
                                  << " denominator=" << failure_record[4]
                                  << " maximum_coefficient="
                                  << failure_record[5]
                                  << " support_floor=" << failure_record[6]
                                  << " scale=" << failure_record[7]
                                  << std::endl;
                    }
                    if (reject_accelerated_trial_before_residual())
                        goto midpoint_iteration_retry;
                    return result;
                }
            }
            trace_live("high_fluxes_dual_u_ready", iter + 1, sub + 1);
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
                if (reject_accelerated_trial_before_residual())
                    goto midpoint_iteration_retry;
                return result;
            }
            trace_live("low_state_checked", iter + 1, sub + 1);

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

            if (controlled_u_fct_flux_injection_enabled_ && fct_enabled_) {
                double local_injection_count = 0.0;
                if (mpi_rank == 0 && nxl >= 1 && h > 0.0) {
                    const int ix = nxl / 2;
                    const double acceleration = bkg_n.charge *
                        fields_mid.Ex[ng + ix] / (bkg_n.mass * Const::c);
                    if (acceleration != 0.0 && std::isfinite(acceleration)) {
                        int donor_j = 1;
                        int donor_k = 0;
                        double donor_mass = -1.0;
                        for (int j = 1; j + 1 < Param::Nv; ++j)
                            for (int k = 0; k < Param::Nmu; ++k) {
                                const double mass = low_state.f[
                                    mass_index(ng + ix, j, k)];
                                if (mass > donor_mass) {
                                    donor_mass = mass;
                                    donor_j = j;
                                    donor_k = k;
                                }
                            }
                        if (donor_mass > 0.0) {
                            const int jf = acceleration > 0.0
                                ? donor_j + 1 : donor_j;
                            const size_t face = uface_index(ix, jf, donor_k);
                            const double added_flux = std::copysign(
                                1.25 * donor_mass / h, acceleration);
                            fu_high[face] += added_flux;
                            cu_high[face] += added_flux / acceleration;
                            local_injection_count = 1.0;
                        }
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
                if (!transport_safe_step_gate(work, sub + 1)) {
                    if (reject_accelerated_trial_before_residual())
                        goto midpoint_iteration_retry;
                    return result;
                }
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
                if (fixed_candidate) {
                    result.fct_limited_u_flux = fu_final;
                    result.fct_limited_u_coefficient = cu_final;
                }
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
                if (!transport_safe_step_gate(work, sub + 1)) {
                    if (reject_accelerated_trial_before_residual())
                        goto midpoint_iteration_retry;
                    return result;
                }
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
            for (int ix = 0; ix < nxl; ++ix) {
                for (int j = 0; j < Param::Nv; ++j) {
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
                    const double low_transfer[4] = {
                        h * fx_low[xface_index(ix, j, k)],
                        h * fx_low[xface_index(ix + 1, j, k)],
                        h * fu_low[uface_index(ix, j, k)],
                        h * fu_low[uface_index(ix, j + 1, k)]};
                    const double face_transfer[4] = {
                        ax_left, ax_right, au_lower, au_upper};
                    const double face_raw_transfer[4] = {
                        ax_left_raw, ax_right_raw, au_lower_raw, au_upper_raw};
                    const double final_tolerance = candidate_roundoff_tolerance(
                        work.f[mass_index(ng + ix, j, k)], low_transfer,
                        m_low, face_raw_transfer, face_transfer);
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
                    donor_low_mass[p] = m_low;
                    donor_limited_outflow[p] =
                        static_cast<double>(limited_outflow_ld);
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

            if (detailed_operator_diagnostics) {
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
                            const double low_transfer[4] = {
                                h * fx_low[xface_index(ix, j, k)],
                                h * fx_low[xface_index(ix + 1, j, k)],
                                h * fu_low[uface_index(ix, j, k)],
                                h * fu_low[uface_index(ix, j + 1, k)]};
                            const double raw_transfer[4] = {
                                h * (fx_high[xface_index(ix, j, k)] -
                                     fx_low[xface_index(ix, j, k)]),
                                h * (fx_high[xface_index(ix + 1, j, k)] -
                                     fx_low[xface_index(ix + 1, j, k)]),
                                h * (fu_high[uface_index(ix, j, k)] -
                                     fu_low[uface_index(ix, j, k)]),
                                h * (fu_high[uface_index(ix, j + 1, k)] -
                                     fu_low[uface_index(ix, j + 1, k)])};
                            const double final_transfer[4] = {
                                ax_left, ax_right, au_lower, au_upper};
                            const double tolerance =
                                candidate_roundoff_tolerance(
                                    work.f[mass_index(ng + ix, j, k)],
                                    low_transfer, m_low, raw_transfer,
                                    final_transfer);
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
            double invariant_min[2] = {
                candidate_min, high_candidate_min};
            MPI_Allreduce(MPI_IN_PLACE, invariant_min, 2, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            candidate_min = invariant_min[0];
            high_candidate_min = invariant_min[1];
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
                if (reject_accelerated_trial_before_residual())
                    goto midpoint_iteration_retry;
                return result;
            }

            // The FCT energy ledger and macro budget are defined before the
            // final dual correction so limiter effects remain distinct from
            // the pairing repair.
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
            accumulate_final_fct_macro_budget();
            if (fixed_candidate) {
                result.fct_limited_u_flux = fu_final;
                result.fct_limited_u_coefficient = cu_final;
            }

            if (dual_u_enabled) {
                current_moment_from_x_flux(fx_final, dual_jn_high_face);
                #pragma omp parallel for schedule(static)
                for (int ix = 0; ix < nxl; ++ix) {
                    dual_target_jn_cell[static_cast<size_t>(ix)] = 0.5 * (
                        dual_jn_high_face[static_cast<size_t>(ix)] +
                        dual_jn_high_face[static_cast<size_t>(ix + 1)]);
                    final_dual_acceleration[static_cast<size_t>(ix)] =
                        bkg_n.charge * fields_mid.Ex[ng + ix] /
                        (bkg_n.mass * Const::c);
                }

                double pre_baseline_ke = 0.0;
                double pre_baseline_work = 0.0;
                if (face_pairing_mode_ == FACE_PAIRING_REGULARIZED) {
                    #pragma omp parallel for schedule(static) reduction(+:pre_baseline_ke,pre_baseline_work)
                    for (int ix = 0; ix < nxl; ++ix) {
                        double je = 0.0;
                        for (int j = 0; j < Param::Nv; ++j)
                            for (int k = 0; k < Param::Nmu; ++k)
                                pre_baseline_ke +=
                                    candidate_mass[
                                        static_cast<size_t>(ix) *
                                            Param::Nvmu + idx2(j, k)] *
                                    bkg_n.cgrid.kinetic_energy[idx2(j, k)];
                        for (int jf = 1; jf < Param::Nv; ++jf)
                            for (int k = 0; k < Param::Nmu; ++k)
                                je += u_face_energy_current_weight[
                                    idx2(jf, k)] *
                                    cu_final[uface_index(ix, jf, k)];
                        pre_baseline_work += h * sg.dx *
                            fields_mid.Ex[ng + ix] * je;
                    }
                }

                DualUCoupling::FinalLimitedDiagnostics final_dual =
                    DualUCoupling::apply_final_limited_capacity_pairing(
                        nxl, h, dual_target_jn_cell,
                        final_dual_acceleration,
                        u_face_energy_current_weight,
                        bkg_n.cgrid.upar_widths,
                        candidate_mass, cu_final, fu_final);

                // The correction is local in x.  Trial iterations need only a
                // single collective validity decision; global norms and counts
                // are reduced once below, for the candidate that is accepted
                // or emitted.
                int final_dual_valid = final_dual.valid;
                MPI_Allreduce(MPI_IN_PLACE, &final_dual_valid, 1, MPI_INT,
                              MPI_MIN, MPI_COMM_WORLD);
                final_dual_iteration_diagnostics.valid =
                    std::min(final_dual_iteration_diagnostics.valid,
                             final_dual_valid);
                final_dual_iteration_diagnostics.target_linf = std::max(
                    final_dual_iteration_diagnostics.target_linf,
                    final_dual.target_linf);
                final_dual_iteration_diagnostics.residual_before_linf =
                    std::max(
                        final_dual_iteration_diagnostics.residual_before_linf,
                        final_dual.residual_before_linf);
                final_dual_iteration_diagnostics.residual_after_linf =
                    std::max(
                        final_dual_iteration_diagnostics.residual_after_linf,
                        final_dual.residual_after_linf);
                final_dual_iteration_diagnostics.minimum_scale = std::min(
                    final_dual_iteration_diagnostics.minimum_scale,
                    final_dual.minimum_scale);
                final_dual_iteration_diagnostics.correction_l2 =
                    std::hypot(
                        final_dual_iteration_diagnostics.correction_l2,
                        final_dual.correction_l2);
                final_dual_iteration_diagnostics.correction_linf = std::max(
                    final_dual_iteration_diagnostics.correction_linf,
                    final_dual.correction_linf);
                final_dual_iteration_diagnostics.candidate_min = std::min(
                    final_dual_iteration_diagnostics.candidate_min,
                    final_dual.candidate_min);
                final_dual_iteration_diagnostics.corrected_cell_count +=
                    final_dual.corrected_cell_count;
                final_dual_iteration_diagnostics.limited_cell_count +=
                    final_dual.limited_cell_count;
                final_dual_iteration_diagnostics.unresolved_cell_count +=
                    final_dual.unresolved_cell_count;
                final_dual_iteration_diagnostics.roundoff_zeroed_count +=
                    final_dual.roundoff_zeroed_count;
                final_dual_iteration_diagnostics.roundoff_zeroed_mass +=
                    final_dual.roundoff_zeroed_mass;

                if (!final_dual_valid ||
                    !std::isfinite(final_dual.candidate_min)) {
                    result.failed = true;
                    result.state_advanced = 0;
                    result.failure_reason = 15;
                    result.failure_iteration = iter;
                    result.failure_substep = sub;
                    result.failure_global_cfl = global_cfl;
                    if (reject_accelerated_trial_before_residual())
                        goto midpoint_iteration_retry;
                    return result;
                }

                if (face_pairing_mode_ == FACE_PAIRING_REGULARIZED) {
                    result.face_pairing_attempted = 1;

                    // This is the complete stable state to restore on every
                    // regularized-solver or trust-region rejection.
                    const std::vector<double> baseline_mass(candidate_mass);
                    const std::vector<double> baseline_coefficient(cu_final);
                    const std::vector<double> baseline_flux(fu_final);

                    std::vector<double> baseline_je(
                        static_cast<size_t>(nxl), 0.0);
                    double baseline_ke = 0.0;
                    double baseline_work = 0.0;
                    double baseline_candidate_min =
                        std::numeric_limits<double>::infinity();
                    double baseline_f_residual_square = 0.0;
                    #pragma omp parallel for schedule(static) reduction(+:baseline_ke,baseline_work,baseline_f_residual_square) reduction(min:baseline_candidate_min)
                    for (int ix = 0; ix < nxl; ++ix) {
                        double je = 0.0;
                        for (int j = 0; j < Param::Nv; ++j) {
                            for (int k = 0; k < Param::Nmu; ++k) {
                                const size_t local_cell =
                                    static_cast<size_t>(ix) * Param::Nvmu +
                                    idx2(j, k);
                                const double mass = baseline_mass[local_cell];
                                baseline_ke += mass *
                                    bkg_n.cgrid.kinetic_energy[idx2(j, k)];
                                baseline_candidate_min = std::min(
                                    baseline_candidate_min, mass);
                                const double difference = mass -
                                    guess.f[
                                        mass_index(ng + ix, j, k)];
                                baseline_f_residual_square +=
                                    difference * difference;
                            }
                        }
                        for (int jf = 1; jf < Param::Nv; ++jf)
                            for (int k = 0; k < Param::Nmu; ++k)
                                je += u_face_energy_current_weight[
                                    idx2(jf, k)] *
                                    baseline_coefficient[
                                        uface_index(ix, jf, k)];
                        baseline_je[static_cast<size_t>(ix)] = je;
                        baseline_work += h * sg.dx *
                            fields_mid.Ex[ng + ix] * je;
                    }

                    std::vector<double> baseline_gstar;
                    PeriodicStaggered::apply_cell_to_face_Gstar(
                        baseline_je, baseline_gstar, nxl, mpi_rank,
                        mpi_size, 972);
                    std::vector<double> local_face_residual(
                        static_cast<size_t>(nxl), 0.0);
                    for (int iface = 0; iface < nxl; ++iface)
                        local_face_residual[static_cast<size_t>(iface)] =
                            dual_jn_high_face[static_cast<size_t>(iface)] -
                            baseline_gstar[static_cast<size_t>(iface)];

                    std::vector<double> lower_capacity;
                    std::vector<double> upper_capacity;
                    int capacity_valid =
                        DualUCoupling::compute_final_pairing_current_bounds(
                            nxl, h, final_dual_acceleration,
                            u_face_energy_current_weight,
                            bkg_n.cgrid.upar_widths, baseline_mass,
                            lower_capacity, upper_capacity) ? 1 : 0;
                    MPI_Allreduce(MPI_IN_PLACE, &capacity_valid, 1, MPI_INT,
                                  MPI_MIN, MPI_COMM_WORLD);

                    std::vector<double> global_residual;
                    std::vector<double> global_lower;
                    std::vector<double> global_upper;
                    allgather_spatial_cells(
                        local_face_residual, sg.nx_global, mpi_size,
                        global_residual);
                    allgather_spatial_cells(
                        lower_capacity, sg.nx_global, mpi_size, global_lower);
                    allgather_spatial_cells(
                        upper_capacity, sg.nx_global, mpi_size, global_upper);
                    if (global_residual.size() !=
                            static_cast<size_t>(sg.nx_global) ||
                        global_lower.size() != global_residual.size() ||
                        global_upper.size() != global_residual.size())
                        capacity_valid = 0;

                    RegularizedFacePairing::Config pairing_config;
                    pairing_config.sigma_cutoff =
                        face_pairing_sigma_cutoff_;
                    pairing_config.lambda = face_pairing_lambda_;
                    pairing_config.eta = face_pairing_eta_;
                    pairing_config.trust_fraction =
                        face_pairing_trust_fraction_;
                    std::vector<double> global_weight(
                        global_residual.size(), 1.0);
                    for (size_t i = 0; i < global_lower.size(); ++i) {
                        global_lower[i] *= face_pairing_trust_fraction_;
                        global_upper[i] *= face_pairing_trust_fraction_;
                    }
                    std::vector<double> global_correction;
                    std::vector<double> unresolved_residual;
                    RegularizedFacePairing::Diagnostics pairing_diagnostics;
                    const bool solved = capacity_valid &&
                        RegularizedFacePairing::solve(
                            global_residual, global_weight, global_lower,
                            global_upper, pairing_config, global_correction,
                            unresolved_residual, pairing_diagnostics);

                    result.face_pairing_solver_converged =
                        pairing_diagnostics.converged;
                    result.face_pairing_iterations =
                        pairing_diagnostics.iterations;
                    result.face_pairing_unresolved_mode_count =
                        pairing_diagnostics.unresolved_mode_count;
                    result.face_pairing_unresolved_mode_l2 =
                        pairing_diagnostics.unresolved_mode_l2;
                    result.face_pairing_correction_l2 =
                        pairing_diagnostics.correction_l2;
                    result.face_pairing_correction_linf =
                        pairing_diagnostics.correction_linf;
                    result.face_pairing_requested_correction_l2 =
                        pairing_diagnostics.correction_l2;
                    result.face_pairing_requested_correction_linf =
                        pairing_diagnostics.correction_linf;
                    result.face_pairing_capacity_active_cells =
                        pairing_diagnostics.capacity_active_cells;
                    result.face_pairing_trust_region_active_cells =
                        pairing_diagnostics.trust_region_active_cells;
                    result.face_pairing_nonzero_capacity_cells =
                        pairing_diagnostics.nonzero_capacity_cells;
                    result.face_pairing_bound_saturated_cells =
                        pairing_diagnostics.bound_saturated_cells;
                    result.face_pairing_objective_residual =
                        pairing_diagnostics.objective_residual;
                    result.face_pairing_objective_smoothness =
                        pairing_diagnostics.objective_smoothness;
                    result.face_pairing_objective_amplitude =
                        pairing_diagnostics.objective_amplitude;
                    result.face_pairing_objective_total =
                        pairing_diagnostics.objective_total;
                    result.face_pairing_pass_solver = solved ? 1 : 0;

                    bool accept_pairing = solved;
                    if (accept_pairing) {
                        std::vector<double> requested_je(baseline_je);
                        for (int ix = 0; ix < nxl; ++ix)
                            requested_je[static_cast<size_t>(ix)] +=
                                global_correction[
                                    static_cast<size_t>(sg.ix_start + ix)];
                        const DualUCoupling::FinalLimitedDiagnostics applied =
                            DualUCoupling::apply_final_limited_capacity_pairing(
                                nxl, h, requested_je,
                                final_dual_acceleration,
                                u_face_energy_current_weight,
                                bkg_n.cgrid.upar_widths,
                                candidate_mass, cu_final, fu_final);
                        int applied_valid = applied.valid;
                        MPI_Allreduce(MPI_IN_PLACE, &applied_valid, 1,
                                      MPI_INT, MPI_MIN, MPI_COMM_WORLD);
                        accept_pairing = applied_valid != 0;
                    }
                    result.face_pairing_pass_apply =
                        accept_pairing ? 1 : 0;
                    result.face_pairing_candidate_valid =
                        accept_pairing ? 1 : 0;

                    std::vector<double> corrected_je(baseline_je);
                    double corrected_ke = baseline_ke;
                    double corrected_work = baseline_work;
                    double corrected_candidate_min =
                        baseline_candidate_min;
                    double corrected_f_residual_square =
                        baseline_f_residual_square;
                    long double local_mass_error_stable = 0.0L;
                    long double local_baseline_mass_stable = 0.0L;
                    double local_cell_mass_error_linf = 0.0;
                    double local_cell_mass_relative_linf = 0.0;
                    if (accept_pairing) {
                        std::fill(corrected_je.begin(), corrected_je.end(),
                                  0.0);
                        corrected_ke = 0.0;
                        corrected_work = 0.0;
                        corrected_candidate_min =
                            std::numeric_limits<double>::infinity();
                        corrected_f_residual_square = 0.0;
                        #pragma omp parallel for schedule(static) reduction(+:corrected_ke,corrected_work,corrected_f_residual_square,local_mass_error_stable,local_baseline_mass_stable) reduction(min:corrected_candidate_min) reduction(max:local_cell_mass_error_linf,local_cell_mass_relative_linf)
                        for (int ix = 0; ix < nxl; ++ix) {
                            double je = 0.0;
                            long double cell_mass_error = 0.0L;
                            long double cell_baseline_mass = 0.0L;
                            for (int j = 0; j < Param::Nv; ++j)
                                for (int k = 0; k < Param::Nmu; ++k) {
                                    const size_t local_cell =
                                        static_cast<size_t>(ix) *
                                            Param::Nvmu + idx2(j, k);
                                    const double mass =
                                        candidate_mass[local_cell];
                                    corrected_ke += mass *
                                        bkg_n.cgrid.kinetic_energy[
                                            idx2(j, k)];
                                    corrected_candidate_min = std::min(
                                        corrected_candidate_min, mass);
                                    const double mass_difference = mass -
                                        baseline_mass[local_cell];
                                    cell_mass_error +=
                                        static_cast<long double>(mass_difference);
                                    cell_baseline_mass += static_cast<long double>(
                                        baseline_mass[local_cell]);
                                    const double difference = mass -
                                        guess.f[
                                            mass_index(ng + ix, j, k)];
                                    corrected_f_residual_square +=
                                        difference * difference;
                                }
                            for (int jf = 1; jf < Param::Nv; ++jf)
                                for (int k = 0; k < Param::Nmu; ++k)
                                    je += u_face_energy_current_weight[
                                        idx2(jf, k)] *
                                        cu_final[uface_index(ix, jf, k)];
                            corrected_je[static_cast<size_t>(ix)] = je;
                            corrected_work += h * sg.dx *
                                fields_mid.Ex[ng + ix] * je;
                            local_mass_error_stable += cell_mass_error;
                            local_baseline_mass_stable += cell_baseline_mass;
                            const double cell_error = std::fabs(
                                static_cast<double>(cell_mass_error));
                            const double cell_scale = std::max(
                                std::fabs(static_cast<double>(cell_baseline_mass)),
                                std::numeric_limits<double>::min());
                            local_cell_mass_error_linf = std::max(
                                local_cell_mass_error_linf, cell_error);
                            local_cell_mass_relative_linf = std::max(
                                local_cell_mass_relative_linf,
                                cell_error / cell_scale);
                        }
                    }

                    long double global_mass_terms[2] = {
                        local_mass_error_stable, local_baseline_mass_stable};
                    MPI_Allreduce(MPI_IN_PLACE, global_mass_terms, 2,
                                  MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    const double global_mass_error =
                        static_cast<double>(global_mass_terms[0]);
                    const double global_mass_scale = std::max(
                        std::fabs(static_cast<double>(global_mass_terms[1])),
                        std::numeric_limits<double>::min());
                    const double global_mass_relative_error =
                        std::fabs(global_mass_error) / global_mass_scale;
                    MPI_Allreduce(MPI_IN_PLACE, &local_cell_mass_error_linf,
                                  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE,
                                  &local_cell_mass_relative_linf, 1,
                                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE, &corrected_candidate_min, 1,
                                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                    result.face_pairing_mass_error = global_mass_error;
                    result.face_pairing_mass_relative_error =
                        global_mass_relative_error;
                    result.face_pairing_cell_mass_error_linf =
                        local_cell_mass_error_linf;
                    result.face_pairing_cell_mass_relative_linf =
                        local_cell_mass_relative_linf;
                    result.face_pairing_candidate_min =
                        corrected_candidate_min;
                    if (result.face_pairing_candidate_valid) {
                        result.face_pairing_candidate_mass_error =
                            global_mass_error;
                        result.face_pairing_candidate_min_before_fallback =
                            corrected_candidate_min;
                    }

                    std::vector<double> corrected_gstar;
                    if (accept_pairing)
                        PeriodicStaggered::apply_cell_to_face_Gstar(
                            corrected_je, corrected_gstar, nxl, mpi_rank,
                            mpi_size, 974);
                    long double local_before = 0.0L;
                    long double local_after = 0.0L;
                    long double local_core_before = 0.0L;
                    long double local_core_after = 0.0L;
                    long double local_energy_residual_scale = 0.0L;
                    long double local_face_work_scale = 0.0L;
                    double local_j_scale = 0.0;
                    for (int iface = 0; iface < nxl; ++iface) {
                        const double before =
                            local_face_residual[
                                static_cast<size_t>(iface)];
                        const double after = accept_pairing
                            ? dual_jn_high_face[
                                  static_cast<size_t>(iface)] -
                                  corrected_gstar[
                                      static_cast<size_t>(iface)]
                            : before;
                        local_before +=
                            static_cast<long double>(before) * before;
                        local_after +=
                            static_cast<long double>(after) * after;
                        const double eface = fields_mid.Ex_face[
                            static_cast<size_t>(iface)];
                        local_energy_residual_scale +=
                            std::fabs(static_cast<long double>(h * sg.dx) *
                                      eface * before);
                        local_face_work_scale +=
                            std::fabs(static_cast<long double>(h * sg.dx) *
                                      eface * baseline_gstar[
                                          static_cast<size_t>(iface)]);
                        const int global_face = sg.ix_start + iface;
                        const double x = global_face * sg.dx;
                        if (x >= kFctCoreBoundaryWidth &&
                            x <= Param::Lx - kFctCoreBoundaryWidth) {
                            local_core_before +=
                                static_cast<long double>(before) * before;
                            local_core_after +=
                                static_cast<long double>(after) * after;
                        }
                        local_j_scale = std::max(
                            local_j_scale,
                            std::max(std::fabs(
                                dual_jn_high_face[
                                    static_cast<size_t>(iface)]),
                                     std::fabs(baseline_gstar[
                                         static_cast<size_t>(iface)])));
                    }
                    double norm_squares[6] = {
                        static_cast<double>(local_before),
                        static_cast<double>(local_after),
                        static_cast<double>(local_core_before),
                        static_cast<double>(local_core_after),
                        static_cast<double>(local_energy_residual_scale),
                        static_cast<double>(local_face_work_scale)};
                    MPI_Allreduce(MPI_IN_PLACE, norm_squares, 6, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE, &local_j_scale, 1,
                                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                    result.face_pairing_residual_before =
                        std::sqrt(std::max(0.0, norm_squares[0]));
                    result.face_pairing_residual_after =
                        std::sqrt(std::max(0.0, norm_squares[1]));
                    result.face_pairing_core_residual_before =
                        std::sqrt(std::max(0.0, norm_squares[2]));
                    result.face_pairing_core_residual_after =
                        std::sqrt(std::max(0.0, norm_squares[3]));
                    result.face_pairing_energy_residual_scale =
                        norm_squares[4];
                    if (result.face_pairing_candidate_valid) {
                        result.face_pairing_candidate_residual_after =
                            result.face_pairing_residual_after;
                        result.face_pairing_candidate_core_residual_after =
                            result.face_pairing_core_residual_after;
                    }

                    double energy_terms[4] = {
                        pre_baseline_ke, baseline_ke,
                        pre_baseline_work, baseline_work};
                    MPI_Allreduce(MPI_IN_PLACE, energy_terms, 4, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE, &baseline_candidate_min, 1,
                                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                    double corrected_terms[4] = {
                        corrected_ke, corrected_work,
                        baseline_f_residual_square,
                        corrected_f_residual_square};
                    MPI_Allreduce(MPI_IN_PLACE, corrected_terms, 4,
                                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    result.face_pairing_delta_ke =
                        corrected_terms[0] - energy_terms[1];
                    result.face_pairing_delta_work =
                        corrected_terms[1] - energy_terms[3];
                    if (result.face_pairing_candidate_valid) {
                        result.face_pairing_candidate_delta_ke =
                            result.face_pairing_delta_ke;
                        result.face_pairing_candidate_delta_work =
                            result.face_pairing_delta_work;
                    }
                    double actual_correction_norms[2] = {0.0, 0.0};
                    if (accept_pairing) {
                        for (int ix = 0; ix < nxl; ++ix) {
                            const double delta =
                                corrected_je[static_cast<size_t>(ix)] -
                                baseline_je[static_cast<size_t>(ix)];
                            actual_correction_norms[0] += delta * delta;
                            actual_correction_norms[1] = std::max(
                                actual_correction_norms[1],
                                std::fabs(delta));
                        }
                    }
                    MPI_Allreduce(MPI_IN_PLACE, actual_correction_norms, 1,
                                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE,
                                  actual_correction_norms + 1, 1,
                                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                    result.face_pairing_correction_l2 = std::sqrt(
                        std::max(0.0, actual_correction_norms[0]));
                    result.face_pairing_correction_linf =
                        actual_correction_norms[1];
                    result.face_pairing_applied_correction_l2 =
                        result.face_pairing_candidate_valid
                        ? result.face_pairing_correction_l2 : 0.0;
                    result.face_pairing_applied_correction_linf =
                        result.face_pairing_candidate_valid
                        ? result.face_pairing_correction_linf : 0.0;
                    const double correction_scale =
                        std::max(1.0, local_j_scale);
                    const double correction_trust_limit =
                        face_pairing_correction_trust_fraction_ *
                        correction_scale;
                    const double face_work_roundoff_floor = 4096.0 *
                        std::numeric_limits<double>::epsilon() *
                        std::max(1.0, norm_squares[5]);
                    const double energy_residual_scale = std::max(
                        result.face_pairing_energy_residual_scale,
                        face_work_roundoff_floor);
                    const double correction_energy_scale = std::max(
                        std::fabs(result.face_pairing_delta_ke),
                        std::fabs(result.face_pairing_delta_work));
                    const double energy_pair_scale = std::max(
                        correction_energy_scale, face_work_roundoff_floor);
                    const double energy_pair_error = std::fabs(
                        result.face_pairing_delta_ke -
                        result.face_pairing_delta_work);
                    const double baseline_f_residual = std::sqrt(
                        std::max(0.0, corrected_terms[2]));
                    const double corrected_f_residual = std::sqrt(
                        std::max(0.0, corrected_terms[3]));
                    const double f_residual_floor = 4096.0 *
                        std::numeric_limits<double>::epsilon() *
                        std::max(1.0, baseline_f_residual);
                    const double candidate_tolerance =
                        4096.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::fabs(
                            baseline_candidate_min));

                    const bool pass_global_residual =
                        result.face_pairing_candidate_valid &&
                        result.face_pairing_residual_after <
                            result.face_pairing_residual_before;
                    const bool pass_core_residual =
                        result.face_pairing_candidate_valid &&
                        result.face_pairing_core_residual_after <
                            result.face_pairing_core_residual_before;
                    const bool pass_correction_trust =
                        result.face_pairing_candidate_valid &&
                        result.face_pairing_correction_linf <=
                            correction_trust_limit;
                    const bool pass_energy_pair =
                        result.face_pairing_candidate_valid &&
                        energy_pair_error <=
                            face_pairing_energy_pair_tolerance_ *
                            energy_pair_scale;
                    const bool pass_energy_residual_scale =
                        result.face_pairing_candidate_valid &&
                        correction_energy_scale <=
                            face_pairing_energy_residual_fraction_ *
                            energy_residual_scale;
                    const bool pass_candidate_min =
                        result.face_pairing_candidate_valid &&
                        corrected_candidate_min + candidate_tolerance >=
                            baseline_candidate_min;
                    const bool pass_mass =
                        result.face_pairing_candidate_valid &&
                        global_mass_relative_error <=
                            face_pairing_mass_relative_tolerance_ &&
                        local_cell_mass_relative_linf <=
                            face_pairing_mass_relative_tolerance_;
                    const bool pass_f_residual =
                        result.face_pairing_candidate_valid &&
                        // The next Picard evaluation supplies the full
                        // nonlinear residual.  At this point bound the
                        // correction-induced growth against the same guess.
                        corrected_f_residual <= baseline_f_residual *
                            (1.0 +
                             face_pairing_f_residual_growth_tolerance_) +
                            f_residual_floor;
                    result.face_pairing_correction_trust_limit =
                        correction_trust_limit;
                    result.face_pairing_correction_trust_ratio =
                        result.face_pairing_correction_linf /
                        std::max(correction_trust_limit,
                                 std::numeric_limits<double>::min());
                    result.face_pairing_energy_pair_error = energy_pair_error;
                    result.face_pairing_energy_pair_relative =
                        energy_pair_error / energy_pair_scale;
                    result.face_pairing_energy_residual_ratio =
                        correction_energy_scale / energy_residual_scale;
                    result.face_pairing_f_residual_relative_growth =
                        (corrected_f_residual - baseline_f_residual) /
                        std::max(baseline_f_residual,
                                 std::numeric_limits<double>::min());
                    result.face_pairing_pass_global_residual =
                        pass_global_residual ? 1 : 0;
                    result.face_pairing_pass_core_residual =
                        pass_core_residual ? 1 : 0;
                    result.face_pairing_pass_correction_trust =
                        pass_correction_trust ? 1 : 0;
                    result.face_pairing_pass_delta_ke =
                        pass_energy_pair ? 1 : 0;
                    result.face_pairing_pass_delta_work =
                        pass_energy_residual_scale ? 1 : 0;
                    result.face_pairing_pass_candidate_min =
                        pass_candidate_min ? 1 : 0;
                    result.face_pairing_pass_mass = pass_mass ? 1 : 0;
                    result.face_pairing_pass_f_residual =
                        pass_f_residual ? 1 : 0;
                    result.face_pairing_pass_energy_pair =
                        pass_energy_pair ? 1 : 0;
                    result.face_pairing_pass_energy_residual_scale =
                        pass_energy_residual_scale ? 1 : 0;

                    unsigned int rejection_mask = 0u;
                    if (!result.face_pairing_pass_solver)
                        rejection_mask |= 1u << 0;
                    if (!result.face_pairing_pass_apply)
                        rejection_mask |= 1u << 1;
                    if (!pass_global_residual)
                        rejection_mask |= 1u << 2;
                    if (!pass_core_residual)
                        rejection_mask |= 1u << 3;
                    if (!pass_correction_trust)
                        rejection_mask |= 1u << 4;
                    if (!pass_energy_pair)
                        rejection_mask |= 1u << 5;
                    if (!pass_energy_residual_scale)
                        rejection_mask |= 1u << 6;
                    if (!pass_candidate_min)
                        rejection_mask |= 1u << 7;
                    if (!pass_mass)
                        rejection_mask |= 1u << 8;
                    if (!pass_f_residual)
                        rejection_mask |= 1u << 9;
                    result.face_pairing_rejection_mask = rejection_mask;

                    accept_pairing =
                        result.face_pairing_candidate_valid &&
                        rejection_mask == 0u;

                    int globally_accepted = accept_pairing ? 1 : 0;
                    MPI_Allreduce(MPI_IN_PLACE, &globally_accepted, 1,
                                  MPI_INT, MPI_MIN, MPI_COMM_WORLD);
                    accept_pairing = globally_accepted != 0;
                    if (accept_pairing) {
                        result.face_pairing_accepted = 1;
                    } else {
                        candidate_mass = baseline_mass;
                        cu_final = baseline_coefficient;
                        fu_final = baseline_flux;
                        result.face_pairing_fallback_to_cell_baseline = 1;
                    }
                }
                result.fct_final_scratch_min = std::min(
                    result.fct_final_scratch_min,
                    final_dual.candidate_min);
            }

            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k)
                    work.f[mass_index(ng + ix, j, k)] =
                        candidate_mass[static_cast<size_t>(ix) * Param::Nvmu +
                                       idx2(j, k)];
            if (!transport_safe_step_gate(work, sub + 1)) {
                if (reject_accelerated_trial_before_residual())
                    goto midpoint_iteration_retry;
                return result;
            }

            // Directional limiter statistics must describe the committed
            // fluxes after donor-beta closure, not the provisional FCT alpha
            // arrays.  This scan is intentionally local; its single global
            // reduction is deferred until a candidate is accepted.
            #pragma omp parallel for schedule(static) reduction(+:x_limiter_faces,x_limiter_active) reduction(min:x_limiter_min)
            for (int iface = 0; iface < nxl; ++iface)
                for (int j = 0; j < Param::Nv; ++j)
                    for (int k = 0; k < Param::Nmu; ++k) {
                        const size_t id = xface_index(iface, j, k);
                        const double anti = fx_high[id] - fx_low[id];
                        const double alpha = anti != 0.0
                            ? (fx_final[id] - fx_low[id]) / anti : 1.0;
                        x_limiter_faces += 1.0;
                        if (alpha < 1.0 - 1.0e-14)
                            x_limiter_active += 1.0;
                        x_limiter_min = std::min(x_limiter_min, alpha);
                    }
            if (fixed_candidate || step_diagnostics_enabled_) {
                #pragma omp parallel for schedule(static) reduction(+:u_limiter_faces,u_limiter_active) reduction(min:u_limiter_min)
                for (int ix = 0; ix < nxl; ++ix)
                    for (int jf = 1; jf < Param::Nv; ++jf)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t id = uface_index(ix, jf, k);
                            const double anti = fu_high[id] - fu_low[id];
                            const double alpha = anti != 0.0
                                ? (fu_final[id] - fu_low[id]) / anti : 1.0;
                            u_limiter_faces += 1.0;
                            if (alpha < 1.0 - 1.0e-14)
                                u_limiter_active += 1.0;
                            u_limiter_min = std::min(u_limiter_min, alpha);
                        }
            }
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
            const std::vector<double>& fu_high_used = low_order_only_
                ? fu_low : fu_high;
            const std::vector<double>& fu_final_used = low_order_only_
                ? fu_low : fu_final;
            const std::vector<double>& cu_final_used = low_order_only_
                ? cu_low : cu_final;

            if (accepted_energy_audit_enabled_) {
                // This is a flux-form diagnostic of this real transport
                // substep.  It intentionally uses no reconstructed state and
                // no operator replay: all three layers are the actual
                // low/high/final face fluxes just committed for this trial.
                #pragma omp parallel
                {
                    std::array<long double, LEDGER_SLOT_COUNT> partial = {};
                    #pragma omp for nowait schedule(static)
                    for (int ix = 0; ix < nxl; ++ix) {
                        for (int j = 0; j < Param::Nv; ++j) {
                            for (int k = 0; k < Param::Nmu; ++k) {
                                const double kinetic = bkg_n.cgrid.kinetic_energy[
                                    idx2(j, k)];
                                const double ppar = bkg_n.mass * Const::c *
                                    bkg_n.cgrid.upar_cells[j];
                                const double low_x = -h * (
                                    fx_low_used[xface_index(ix + 1, j, k)] -
                                    fx_low_used[xface_index(ix, j, k)]);
                                const double high_x = -h * (
                                    fx_high_used[xface_index(ix + 1, j, k)] -
                                    fx_high_used[xface_index(ix, j, k)]);
                                const double final_x = -h * (
                                    fx_final_used[xface_index(ix + 1, j, k)] -
                                    fx_final_used[xface_index(ix, j, k)]);
                                const double low_u = -h * (
                                    fu_low[uface_index(ix, j + 1, k)] -
                                    fu_low[uface_index(ix, j, k)]);
                                const double high_u = -h * (
                                    fu_high_used[uface_index(ix, j + 1, k)] -
                                    fu_high_used[uface_index(ix, j, k)]);
                                const double final_u = -h * (
                                    fu_final_used[uface_index(ix, j + 1, k)] -
                                    fu_final_used[uface_index(ix, j, k)]);
                                const double low_delta = low_x + low_u;
                                const double high_delta = high_x + high_u;
                                const double final_delta = final_x + final_u;
                                partial[LEDGER_DKE_LOW] +=
                                    static_cast<long double>(kinetic) * low_delta;
                                partial[LEDGER_DKE_HIGH] +=
                                    static_cast<long double>(kinetic) * high_delta;
                                partial[LEDGER_DKE_FINAL] +=
                                    static_cast<long double>(kinetic) * final_delta;
                                partial[LEDGER_DKE_FCT_X] +=
                                    static_cast<long double>(kinetic) *
                                    (final_x - high_x);
                                partial[LEDGER_DKE_FCT_U] +=
                                    static_cast<long double>(kinetic) *
                                    (final_u - high_u);
                                partial[LEDGER_DMASS_FCT_X] += final_x - high_x;
                                partial[LEDGER_DMASS_FCT_U] += final_u - high_u;
                                partial[LEDGER_DPPAR_FCT_X] +=
                                    static_cast<long double>(ppar) *
                                    (final_x - high_x);
                                partial[LEDGER_DPPAR_FCT_U] +=
                                    static_cast<long double>(ppar) *
                                    (final_u - high_u);
                            }
                        }
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const double flux_lower =
                                fu_final_used[uface_index(ix, 0, k)];
                            const double flux_upper = fu_final_used[
                                uface_index(ix, Param::Nv, k)];
                            const double ppar_lower = bkg_n.mass * Const::c *
                                bkg_n.cgrid.upar_cells[0];
                            const double ppar_upper = bkg_n.mass * Const::c *
                                bkg_n.cgrid.upar_cells[Param::Nv - 1];
                            partial[LEDGER_BNUM_LOWER] += -h * flux_lower;
                            partial[LEDGER_BNUM_UPPER] += h * flux_upper;
                            partial[LEDGER_BPPAR_LOWER] +=
                                -h * ppar_lower * flux_lower;
                            partial[LEDGER_BPPAR_UPPER] +=
                                h * ppar_upper * flux_upper;
                            partial[LEDGER_BENERGY_LOWER] += -h *
                                bkg_n.cgrid.kinetic_energy[idx2(0, k)] *
                                flux_lower;
                            partial[LEDGER_BENERGY_UPPER] += h *
                                bkg_n.cgrid.kinetic_energy[
                                    idx2(Param::Nv - 1, k)] * flux_upper;
                        }
                    }
                    #pragma omp critical(accepted_energy_ledger_accumulate)
                    {
                        for (int q = 0; q < LEDGER_SLOT_COUNT; ++q)
                            trial_energy_ledger_local[static_cast<size_t>(q)] +=
                                partial[static_cast<size_t>(q)];
                    }
                }
            }
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
                    const size_t id = uface_index(ix, jf, k);
                    const size_t face = idx2(jf, k);
                    const double dke = u_face_energy_delta[face];
                    const double energy_current_weight =
                        u_face_energy_current_weight[face];
                    const double coefficient_high = cu_high_used[id];
                    const double coefficient_center = cu_high_center_used[id];
                    const double coefficient_final = cu_final_used[id];
                    const double coefficient_low = cu_low[id];
                    const double flux_final = fu_final_used[id];
                    // J_E = q/(m c dx) sum(delta_K * C_u).  The dx division
                    // converts the cell-integrated C_u back to a current
                    // density; h-weighting and /dt occur outside this loop.
                    je_low += energy_current_weight * coefficient_low;
                    je_high += energy_current_weight * coefficient_high;
                    je_center += energy_current_weight * coefficient_center;
                    je_final += energy_current_weight * coefficient_final;
                    if (fixed_candidate) {
                        const double low_contribution =
                            energy_current_weight * coefficient_low;
                        const double high_contribution =
                            energy_current_weight * coefficient_high;
                        const double final_contribution =
                            energy_current_weight * coefficient_final;
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
                    u_momentum += u_face_momentum_delta[face] * flux_final;
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
                    u_boundary_energy_lower += -h * ke_lo * f_lo;
                    u_boundary_energy_upper += h * ke_hi * f_hi;
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
                // Diagnostic-only physical weighting of FCT activity.  This
                // reads the already-finalized low/high/final fluxes and never
                // feeds back into the transport update.
                substep_audit.fct_weighted_coverage.assign(2 * 3 * 6, 0.0);
                substep_audit.fct_weighted_counts.assign(2 * 3 * 3, 0.0);
                double local_peak = 0.0;
                for (int ix = 0; ix < nxl; ++ix)
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t p = static_cast<size_t>(ix) *
                                Param::Nvmu + idx2(j, k);
                            const double volume = sg.dx *
                                bkg_n.cgrid.cell_phase_volume(j, k);
                            if (volume > 0.0)
                                local_peak = std::max(local_peak,
                                    std::max(0.0, donor_low_mass[p]) / volume);
                        }
                double global_peak = 0.0;
                MPI_Allreduce(&local_peak, &global_peak, 1, MPI_DOUBLE,
                              MPI_MAX, MPI_COMM_WORLD);
                const double tail_floor = std::max(
                    std::numeric_limits<double>::min(), global_peak * 1.0e-8);
                const int seam_cells = std::max(1, static_cast<int>(std::ceil(
                    0.1 * Const::micro / sg.dx)));
                const double activity_tolerance = 1024.0 *
                    std::numeric_limits<double>::epsilon();
                for (int ix = 0; ix < nxl; ++ix) {
                    const int global_ix = sg.ix_start + ix;
                    const int x_region = (global_ix < seam_cells ||
                        global_ix >= sg.nx_global - seam_cells) ? 1 : 0;
                    for (int j = 0; j < Param::Nv; ++j)
                        for (int k = 0; k < Param::Nmu; ++k) {
                            const size_t p = static_cast<size_t>(ix) *
                                Param::Nvmu + idx2(j, k);
                            const double mass = std::max(0.0, donor_low_mass[p]);
                            const double volume = sg.dx *
                                bkg_n.cgrid.cell_phase_volume(j, k);
                            const double fbar = volume > 0.0 ? mass / volume : 0.0;
                            const double upar = bkg_n.cgrid.upar_cells[j];
                            const double uperp = bkg_n.cgrid.uperp_cells[k];
                            const double umag = std::sqrt(upar * upar + uperp * uperp);
                            const double unit_half_width = std::max(0.10,
                                std::max(bkg_n.cgrid.upar_widths[j],
                                         bkg_n.cgrid.uperp_widths[k]));
                            const int velocity_region =
                                std::fabs(umag - 1.0) <= unit_half_width ? 1 :
                                (fbar >= tail_floor ? 0 : 2);
                            bool active = false;
                            const size_t xfaces[2] = {
                                xface_index(ix, j, k), xface_index(ix + 1, j, k)};
                            const size_t ufaces[2] = {
                                uface_index(ix, j, k), uface_index(ix, j + 1, k)};
                            for (int face = 0; face < 2 && !active; ++face) {
                                const double dx_high = fx_high[xfaces[face]] -
                                    fx_low[xfaces[face]];
                                const double du_high = fu_high[ufaces[face]] -
                                    fu_low[ufaces[face]];
                                const bool x_changed = std::fabs(dx_high) >
                                    activity_tolerance * std::max(1.0,
                                        std::max(std::fabs(fx_high[xfaces[face]]),
                                                 std::fabs(fx_low[xfaces[face]]))) &&
                                    std::fabs(fx_final[xfaces[face]] -
                                              fx_high[xfaces[face]]) >
                                    activity_tolerance * std::max(1.0,
                                                                  std::fabs(dx_high));
                                const bool u_changed = std::fabs(du_high) >
                                    activity_tolerance * std::max(1.0,
                                        std::max(std::fabs(fu_high[ufaces[face]]),
                                                 std::fabs(fu_low[ufaces[face]]))) &&
                                    std::fabs(fu_final[ufaces[face]] -
                                              fu_high[ufaces[face]]) >
                                    activity_tolerance * std::max(1.0,
                                                                  std::fabs(du_high));
                                active = x_changed || u_changed;
                            }
                            const double gamma = std::sqrt(1.0 + umag * umag);
                            const double vx = Const::c * upar / gamma;
                            const double current_weight =
                                std::fabs(bkg_n.charge * vx) * mass;
                            const double energy_weight =
                                bkg_n.cgrid.kinetic_energy[idx2(j, k)] * mass;
                            const size_t bin = static_cast<size_t>(
                                x_region * 3 + velocity_region);
                            double* weights = &substep_audit.fct_weighted_coverage[
                                bin * 6];
                            double* counts = &substep_audit.fct_weighted_counts[
                                bin * 3];
                            counts[0] += 1.0;
                            weights[0] += mass;
                            weights[2] += current_weight;
                            weights[4] += energy_weight;
                            if (active) {
                                counts[1] += 1.0;
                                weights[1] += mass;
                                weights[3] += current_weight;
                                weights[5] += energy_weight;
                            }
                            if (donor_beta[p] < 1.0 - 1.0e-14 &&
                                donor_low_mass[p] > 0.0)
                                counts[2] += 1.0;
                        }
                }
                result.coupling_substep_seam_audit.push_back(substep_audit);
            }
            trace_live("transport_substep_end", iter + 1, sub + 1);
        }

        close_periodic_face_blocks(psi_k_x, nxl, 1, mpi_rank, mpi_size, 909);
        close_periodic_face_blocks(psi_p_x, nxl, 1, mpi_rank, mpi_size, 910);
        #pragma omp parallel for schedule(static)
        for (int ix = 0; ix < nxl; ++ix) {
            double ke = 0.0;
            double pp = 0.0;
            long double delta_ke = 0.0L;
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = mass_index(ng + ix, j, k);
                const double mass = work.f[id];
                ke += bkg_n.cgrid.kinetic_energy[idx2(j, k)] * mass;
                pp += bkg_n.mass * Const::c * bkg_n.cgrid.upar_cells[j] * mass;
                // Accumulate the small accepted-state increment directly.
                // Subtracting separately summed O(K) totals loses this signal
                // when the update is close to machine precision.
                delta_ke += static_cast<long double>(
                    bkg_n.cgrid.kinetic_energy[idx2(j, k)]) *
                    (static_cast<long double>(mass) -
                     static_cast<long double>(bkg_n.f[id]));
            }
            final_ke_cell[ix] = ke;
            final_p_cell[ix] = pp;
            delta_ke_cell[ix] = static_cast<double>(delta_ke);
            delta_p_cell[ix] = pp - initial_p_cell[ix];
            local_energy_residual[ix] = delta_ke_cell[ix] +
                psi_k_x[ix + 1] - psi_k_x[ix] - u_energy_cell[ix] -
                u_boundary_energy_cell[ix];
            local_momentum_residual[ix] = delta_p_cell[ix] +
                psi_p_x[ix + 1] - psi_p_x[ix] - u_momentum_cell[ix] -
                u_boundary_momentum_cell[ix];
        }
        if (accepted_energy_audit_enabled_) {
            // delta_ke_cell is formed from the actual endpoint state rather
            // than from an independently reconstructed flux state.  Keep the
            // local sum in long double so the accepted ledger can expose any
            // FV/state mismatch instead of hiding it in cancellation noise.
            for (int ix = 0; ix < nxl; ++ix) {
                trial_actual_delta_ke_local +=
                    static_cast<long double>(delta_ke_cell[ix]);
            }
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
                        je += u_face_energy_current_weight[idx2(jf, k)] *
                            audit_cu[uface_index(ix, jf, k)];
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
        // A fixed operator audit evaluates the FV/current operator at an
        // externally supplied endpoint. It must not perform a second Ampere
        // update or advance/reinject Beam particles.
        EMFields fields_new = (fixed_candidate && fixed_fields_end)
            ? *fixed_fields_end : fields_n;
        if (!fixed_candidate) {
            trace_live("ampere_begin", iter + 1, -1);
            fields_new.advance_ampere_face_from_midpoint_current(
                result.j_total_face_mid, dt, mpi_rank, mpi_size);
            trace_live("ampere_end", iter + 1, -1);
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
        double local_df = 0.0, local_f = 0.0;
        for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t p = mass_index(ng + ix, j, k);
                local_df = std::max(local_df, std::fabs(work.f[p] - guess.f[p]));
                local_f = std::max(local_f, std::fabs(work.f[p]));
            }
        // All six values use the same associative MAX reduction and are
        // independent.  Combining them removes one global synchronization
        // from every coupled midpoint iteration without changing any norm.
        double norms[6] = {local_de, local_e, local_dj, local_j,
                           local_df, local_f};
        MPI_Allreduce(MPI_IN_PLACE, norms, 6, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        result.residual_E = norms[0] / std::max(1.0, norms[1]);
        result.residual_J_bkg = norms[2] /
            std::max(1.0, std::max(std::fabs(Param::jb), norms[3]));
        result.residual_f = norms[4] / std::max(1.0, norms[5]);
        if (midpoint_iteration_trace_for_test_) {
            result.midpoint_residual_e_history.push_back(result.residual_E);
            result.midpoint_residual_j_bkg_history.push_back(result.residual_J_bkg);
            result.midpoint_residual_j_beam_history.push_back(result.residual_J_beam);
        }
        if (midpoint_iteration_trace_for_test_)
            result.midpoint_residual_f_history.push_back(result.residual_f);
        result.nonlinear_residual = std::max(result.residual_E / field_tol,
            result.residual_J_bkg / current_tol);
        if (midpoint_iteration_trace_for_test_ &&
            midpoint_acceleration_enabled) {
            result.midpoint_acceleration_omega_history.push_back(omega);
            result.midpoint_acceleration_residual_before_history.push_back(
                result.nonlinear_residual);
            result.midpoint_acceleration_status_history.push_back(0);
        }
        result.max_residual_E = std::max(result.max_residual_E,
            result.residual_E);
        result.max_residual_J_bkg = std::max(result.max_residual_J_bkg,
            result.residual_J_bkg);
        result.max_residual_f = std::max(result.max_residual_f,
            result.residual_f);

        result.beam_continuity_residual = global_beam_continuity_residual;
        if (live_trace && mpi_rank == 0) {
            std::printf(
                "[midpoint-live] t_fs=%.16e wall_s=%.6f iter=%d sub=-1 "
                "stage=residuals E=%.16e J_bkg=%.16e J_beam=%.16e "
                "f=%.16e beam_continuity=%.16e\n",
                time_fs, MPI_Wtime() - live_trace_start, iter + 1,
                result.residual_E, result.residual_J_bkg,
                result.residual_J_beam, result.residual_f,
                result.beam_continuity_residual);
            std::fflush(stdout);
        }
        trace_live("all_finite_begin", iter + 1, -1);
        int finite_candidate_flag =
            all_finite(work, fields_new, result.j_total_face_mid) ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &finite_candidate_flag, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);
        const bool finite_candidate = finite_candidate_flag != 0;
        trace_live("all_finite_end", iter + 1, -1);
        if (!finite_candidate) {
            if (acceleration_trial_pending) {
                ++result.acceleration_rejected_nonfinite;
                ++result.acceleration_fallback_evaluations;
                ++result.acceleration_history_resets;
                if (midpoint_iteration_trace_for_test_) {
                    result.midpoint_acceleration_omega_history.back() =
                        acceleration_trial_coefficient;
                    result.midpoint_acceleration_residual_before_history.back() =
                        acceleration_residual_before;
                    result.midpoint_acceleration_status_history.back() = 4;
                }
                e_end = acceleration_fallback_e_;
                acceleration_trial_pending = false;
                fallback_plain_evaluation = true;
                have_acceleration_residual = false;
                anderson_history_count = 0;
                aitken_previous_omega = omega;
                --iter;
                continue;
            }
            result.failed = true;
            result.state_advanced = 0;
            result.failure_reason = 4;
            result.failure_iteration = iter;
            result.failure_global_cfl = global_cfl;
            return result;
        }
        const bool strict_candidate = have_previous &&
            result.residual_E < field_tol &&
            result.residual_J_bkg < current_tol;
        if (acceleration_trial_pending) {
            const bool residual_improved = strict_candidate ||
                (std::isfinite(result.nonlinear_residual) &&
                 result.nonlinear_residual <= acceleration_accept_ratio_ *
                    acceleration_residual_before);
            if (!residual_improved) {
                ++result.acceleration_rejected_residual;
                ++result.acceleration_fallback_evaluations;
                ++result.acceleration_history_resets;
                if (midpoint_iteration_trace_for_test_) {
                    result.midpoint_acceleration_omega_history.back() =
                        acceleration_trial_coefficient;
                    result.midpoint_acceleration_residual_before_history.back() =
                        acceleration_residual_before;
                    result.midpoint_acceleration_status_history.back() = 3;
                }
                e_end = acceleration_fallback_e_;
                acceleration_trial_pending = false;
                fallback_plain_evaluation = true;
                have_acceleration_residual = false;
                anderson_history_count = 0;
                aitken_previous_omega = omega;
                --iter;
                continue;
            }
            ++result.acceleration_accepted;
            if (midpoint_iteration_trace_for_test_) {
                result.midpoint_acceleration_omega_history.back() =
                    acceleration_trial_coefficient;
                result.midpoint_acceleration_residual_before_history.back() =
                    acceleration_residual_before;
                result.midpoint_acceleration_status_history.back() = 2;
            }
            acceleration_trial_pending = false;
        }
        const bool final_candidate = fixed_candidate ||
            iter + 1 == max_iters;
        if (strict_candidate || final_candidate) {
            if (dual_u_enabled && fct_enabled_) {
                double final_dual_max[4] = {
                    final_dual_iteration_diagnostics.target_linf,
                    final_dual_iteration_diagnostics.residual_before_linf,
                    final_dual_iteration_diagnostics.residual_after_linf,
                    final_dual_iteration_diagnostics.correction_linf};
                MPI_Allreduce(MPI_IN_PLACE, final_dual_max, 4, MPI_DOUBLE,
                              MPI_MAX, MPI_COMM_WORLD);
                double final_dual_min[2] = {
                    final_dual_iteration_diagnostics.minimum_scale,
                    final_dual_iteration_diagnostics.candidate_min};
                MPI_Allreduce(MPI_IN_PLACE, final_dual_min, 2, MPI_DOUBLE,
                              MPI_MIN, MPI_COMM_WORLD);
                // Counts remain far below the exact-integer range of double,
                // so they can share the sum collective with the two norms.
                double final_dual_sum[6] = {
                    final_dual_iteration_diagnostics.correction_l2 *
                        final_dual_iteration_diagnostics.correction_l2,
                    final_dual_iteration_diagnostics.roundoff_zeroed_mass,
                    static_cast<double>(
                        final_dual_iteration_diagnostics.corrected_cell_count),
                    static_cast<double>(
                        final_dual_iteration_diagnostics.limited_cell_count),
                    static_cast<double>(
                        final_dual_iteration_diagnostics.unresolved_cell_count),
                    static_cast<double>(
                        final_dual_iteration_diagnostics.roundoff_zeroed_count)};
                MPI_Allreduce(MPI_IN_PLACE, final_dual_sum, 6, MPI_DOUBLE,
                              MPI_SUM, MPI_COMM_WORLD);

                result.final_dual_u_valid =
                    final_dual_iteration_diagnostics.valid;
                result.final_dual_u_target_linf =
                    final_dual_max[0];
                result.final_dual_u_residual_before_linf =
                    final_dual_max[1];
                result.final_dual_u_residual_after_linf =
                    final_dual_max[2];
                result.final_dual_u_minimum_scale =
                    final_dual_min[0];
                result.final_dual_u_correction_l2 =
                    std::sqrt(std::max(0.0, final_dual_sum[0]));
                result.final_dual_u_correction_linf =
                    final_dual_max[3];
                result.final_dual_u_candidate_min =
                    final_dual_min[1];
                result.final_dual_u_corrected_cell_count =
                    static_cast<long long>(std::llround(final_dual_sum[2]));
                result.final_dual_u_limited_cell_count =
                    static_cast<long long>(std::llround(final_dual_sum[3]));
                result.final_dual_u_unresolved_cell_count =
                    static_cast<long long>(std::llround(final_dual_sum[4]));
                result.fct_roundoff_zeroed_count +=
                    static_cast<long long>(std::llround(final_dual_sum[5]));
                result.fct_roundoff_zeroed_mass += final_dual_sum[1];
            }
            finalize_directional_limiter_statistics();
            finalize_fct_macro_budget();
        }
        // Emit expensive closure data once, for the candidate that will be
        // committed.  Earlier Picard trials need only transport, Ampere and
        // the convergence/finite checks above.
        const bool collect_expensive_iteration_audit = fixed_candidate ||
            ((step_diagnostics_enabled_ || accepted_energy_audit_enabled_) &&
             (strict_candidate || final_candidate));
        if (collect_expensive_iteration_audit) {
            trace_live("expensive_audit_begin", iter + 1, -1);
            assemble_gstar_face(je_low_cell,
                                result.j_bkg_energy_low_debug_face, 906);
            assemble_gstar_face(je_center_cell,
                                result.j_bkg_energy_center_debug_face, 908);
            assemble_gstar_face(je_high_cell,
                                result.j_bkg_energy_high_debug_face, 910);
            assemble_gstar_face(je_cell, result.j_bkg_energy_debug_face, 912);
            // Face 0 is owned by rank 0; the final rank's face nxl is its
            // periodic alias.  This is diagnostic-only topology accounting.
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
            trace_live("gstar_end", iter + 1, -1);
        }

        // Stage-5, limiter and regional ledgers are not part of the midpoint
        // acceptance test.  Evaluate them only for an emitted production
        // diagnostic or a fixed-state audit, rather than for every trial.
        // Kept in outer scope because accepted-state regional diagnostics
        // consume them below.  They remain zero on a non-diagnostic trial and
        // are never used by the acceptance decision.
        double energy_pair[4] = {0.0, 0.0, 0.0, 0.0};
        double local_stage5[7] = {0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0};
        work.current_face_x = result.j_bkg_face_mid;
        last_work = work;
        if (collect_expensive_iteration_audit) {
        work.compute_moments();
        trace_live("moments_end", iter + 1, -1);
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
        energy_pair[0] = local_cell_work;
        energy_pair[1] = local_face_work;
        energy_pair[2] = local_high_cell_work;
        energy_pair[3] = local_center_cell_work;

        local_stage5[4] = psi_k_x[nxl] - psi_k_x[0];
        local_stage5[5] = psi_p_x[nxl] - psi_p_x[0];
        double local_mass_change = 0.0;
        double local_mass_scale = 0.0;
        double local_momentum_scale = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            local_stage5[0] += local_energy_residual[ix];
            local_stage5[1] += local_momentum_residual[ix];
            local_stage5[2] += u_energy_cell[ix];
            local_stage5[3] += u_momentum_cell[ix];
            local_stage5[6] += delta_ke_cell[ix];
            for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    local_mass_change += work.f[
                        mass_index(ng + ix, j, k)] - bkg_n.f[
                        mass_index(ng + ix, j, k)];
                    local_mass_scale += std::fabs(work.f[
                        mass_index(ng + ix, j, k)]) + std::fabs(bkg_n.f[
                        mass_index(ng + ix, j, k)]);
                }
            local_momentum_scale += std::fabs(initial_p_cell[ix]) +
                                    std::fabs(final_p_cell[ix]);
        }
        double u_boundary_global[5] = {u_boundary_particle, u_boundary_momentum,
                                       u_boundary_energy,
                                       u_boundary_energy_lower,
                                       u_boundary_energy_upper};
        // All quantities below use the same communicator, datatype and SUM
        // operation.  A single packed collective preserves each scalar
        // reduction while removing five latency-dominated synchronization
        // points from diagnostic iterations.
        double stage5_sum[19];
        for (int i = 0; i < 4; ++i) stage5_sum[i] = energy_pair[i];
        for (int i = 0; i < 7; ++i) stage5_sum[4 + i] = local_stage5[i];
        stage5_sum[11] = local_mass_change;
        stage5_sum[12] = local_mass_scale;
        stage5_sum[13] = local_momentum_scale;
        for (int i = 0; i < 5; ++i)
            stage5_sum[14 + i] = u_boundary_global[i];
        MPI_Allreduce(MPI_IN_PLACE, stage5_sum, 19, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        for (int i = 0; i < 4; ++i) energy_pair[i] = stage5_sum[i];
        for (int i = 0; i < 7; ++i) local_stage5[i] = stage5_sum[4 + i];
        local_mass_change = stage5_sum[11];
        local_mass_scale = stage5_sum[12];
        local_momentum_scale = stage5_sum[13];
        for (int i = 0; i < 5; ++i)
            u_boundary_global[i] = stage5_sum[14 + i];
        result.energy_residual_bkg =
            dt * (energy_pair[0] - energy_pair[1]);
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
        result.u_boundary_energy_lower = u_boundary_global[3];
        result.u_boundary_energy_upper = u_boundary_global[4];
        result.stage2_mass_change = local_mass_change;
        result.stage2_mass_scale = local_mass_scale;
        result.stage2_mass_residual = result.stage2_mass_change +
            u_boundary_global[0];
        result.stage2_momentum_change = 0.0;
        result.stage2_momentum_scale = local_momentum_scale;
        for (int ix = 0; ix < nxl; ++ix)
            result.stage2_momentum_change += delta_p_cell[ix];
        MPI_Allreduce(MPI_IN_PLACE, &result.stage2_momentum_change, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        result.stage2_momentum_residual = local_stage5[1];
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
        trace_live("stage5_reduce_end", iter + 1, -1);
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
        trace_live("limiter_reduce_end", iter + 1, -1);
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
        result.field_work_bkg = -dt * energy_pair[1];
        result.energy_pair_residual_bkg = result.energy_residual_bkg;
        trace_live("expensive_audit_end", iter + 1, -1);
        } // collect_expensive_iteration_audit
        result.delta_ke_beam = beam_ke_after - beam_ke_before;
        result.field_work_beam = -integrate_face_work(jbeam, beam_field, sg, dt);
        result.beam_lag_energy_residual =
            -integrate_face_work(jbeam, fields_mid, sg, dt) - result.field_work_beam;
        result.nonlinear_iterations = iter + 1;
        const auto commit_accepted_energy_ledger =
            [&](const std::array<long double, LEDGER_SLOT_COUNT>& local_ledger,
                long double local_actual_delta_ke, bool strict) {
                if (!accepted_energy_audit_enabled_) return;

                // One packed SUM is deliberately deferred until this final
                // candidate.  Rejected Picard trials never enter the ledger.
                double packed[LEDGER_SLOT_COUNT + 2] = {};
                for (int q = 0; q < LEDGER_SLOT_COUNT; ++q) {
                    packed[q] = static_cast<double>(
                        local_ledger[static_cast<size_t>(q)]);
                }
                packed[LEDGER_SLOT_COUNT] =
                    static_cast<double>(local_actual_delta_ke);
                long double local_work = 0.0L;
                for (int iface = 0; iface < nxl; ++iface) {
                    local_work += -static_cast<long double>(dt) * sg.dx *
                        fields_mid.Ex_face[static_cast<size_t>(iface)] *
                        result.j_bkg_face_mid[static_cast<size_t>(iface)];
                }
                packed[LEDGER_SLOT_COUNT + 1] =
                    static_cast<double>(local_work);
                MPI_Allreduce(MPI_IN_PLACE, packed,
                              LEDGER_SLOT_COUNT + 2, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);

                AcceptedEnergyLedger& ledger = result.accepted_energy_ledger;
                ledger.valid = 1;
                ledger.strict_accepted = strict ? 1 : 0;
                ledger.soft_accepted = strict ? 0 : 1;
                ledger.transport_substeps = nsub;
                ledger.delta_ke_after_low = packed[LEDGER_DKE_LOW];
                ledger.delta_ke_after_high = packed[LEDGER_DKE_HIGH];
                ledger.delta_ke_after_final = packed[LEDGER_DKE_FINAL];
                ledger.delta_ke_fct_x = packed[LEDGER_DKE_FCT_X];
                ledger.delta_ke_fct_u = packed[LEDGER_DKE_FCT_U];
                ledger.delta_mass_fct_x = packed[LEDGER_DMASS_FCT_X];
                ledger.delta_mass_fct_u = packed[LEDGER_DMASS_FCT_U];
                ledger.delta_ppar_fct_x = packed[LEDGER_DPPAR_FCT_X];
                ledger.delta_ppar_fct_u = packed[LEDGER_DPPAR_FCT_U];
                ledger.upar_boundary_number_lower = packed[LEDGER_BNUM_LOWER];
                ledger.upar_boundary_number_upper = packed[LEDGER_BNUM_UPPER];
                ledger.upar_boundary_ppar_lower = packed[LEDGER_BPPAR_LOWER];
                ledger.upar_boundary_ppar_upper = packed[LEDGER_BPPAR_UPPER];
                ledger.upar_boundary_energy_lower =
                    packed[LEDGER_BENERGY_LOWER];
                ledger.upar_boundary_energy_upper =
                    packed[LEDGER_BENERGY_UPPER];
                ledger.uperp_boundary_number_upper = 0.0;
                ledger.uperp_boundary_ppar_upper = 0.0;
                ledger.uperp_boundary_energy_upper = 0.0;
                ledger.delta_ke_bkg_actual = packed[LEDGER_SLOT_COUNT];
                ledger.work_ampere_bkg = packed[LEDGER_SLOT_COUNT + 1];
                ledger.residual_fv = result.stage5_r_fv;
                ledger.residual_coupling = result.stage5_r_couple;
                ledger.residual_unexplained = ledger.delta_ke_bkg_actual +
                    ledger.work_ampere_bkg +
                    ledger.upar_boundary_energy_lower +
                    ledger.upar_boundary_energy_upper +
                    ledger.uperp_boundary_energy_upper;
                // Independent stage-5 components use the same final fluxes:
                // R_FV - R_couple must reconstruct the physical background
                // energy balance with the program's established signs.
                ledger.residual_reconstruction_error =
                    ledger.residual_unexplained -
                    (ledger.residual_fv - ledger.residual_coupling);
                ledger.flux_telescope_error = ledger.delta_ke_after_final -
                    (ledger.delta_ke_after_high + ledger.delta_ke_fct_x +
                     ledger.delta_ke_fct_u);
            };
        const auto finalize_coupling_regions = [&]() {
            trace_live("regional_audit_begin", iter + 1, -1);
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
            trace_live("regional_audit_end", iter + 1, -1);
        };
        if (strict_candidate || final_candidate) {
            accepted_energy_ledger_local = trial_energy_ledger_local;
            accepted_energy_actual_delta_ke_local =
                trial_actual_delta_ke_local;
            commit_accepted_energy_ledger(accepted_energy_ledger_local,
                                          accepted_energy_actual_delta_ke_local,
                                          strict_candidate || fixed_candidate);
        }
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
            if (dual_u_enabled)
                result.legacy_center_u_coefficient = cu_high_center;
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
            result.low_u_coefficient = cu_low;
            result.high_u_coefficient = low_order_only_ ? cu_low : cu_high;
            result.final_u_coefficient = low_order_only_ ? cu_low : cu_final;
            result.fct_donor_beta = donor_beta;
            result.fct_donor_low_mass = donor_low_mass;
            result.fct_donor_limited_outflow = donor_limited_outflow;
            if (step_diagnostics_enabled_ || fixed_candidate)
                finalize_coupling_regions();
            return result;
        }
        if (strict_candidate) {
            trace_live("transport_safe_gate_begin", iter + 1, -1);
            const bool transport_safe =
                transport_safe_step_gate(work, nsub);
            trace_live("transport_safe_gate_end", iter + 1, -1);
            if (!transport_safe) {
                if (reject_accelerated_trial_before_residual())
                    goto midpoint_iteration_retry;
                return result;
            }
            if (detailed_operator_diagnostics) {
                trace_live("accepted_transport_capture_begin", iter + 1, -1);
                capture_accepted_transport(work);
                trace_live("accepted_transport_capture_end", iter + 1, -1);
            }
            result.species_np1 = work;
            result.beam_np1 = beam_predictor;
            result.fields_np1 = fields_new;
            result.converged = true;
            result.state_advanced = 1;
            trace_live("strict_accept", iter + 1, -1);
            if (step_diagnostics_enabled_ || fixed_candidate)
                finalize_coupling_regions();
            return result;
        }

        bool accelerated_next_endpoint = false;
        bool coefficient_rejected = false;
        double next_acceleration_coefficient = omega;
        if (midpoint_acceleration_enabled) {
            for (int iface = 0; iface < nxl; ++iface) {
                acceleration_re_current_[iface] = next_e[iface] - e_end[iface];
                acceleration_fallback_e_[iface] = e_end[iface] +
                    omega * acceleration_re_current_[iface];
            }

            const bool can_attempt_acceleration =
                !fallback_plain_evaluation && iter + 1 < max_iters &&
                iter + 1 >= acceleration_start_iter_;
            if (can_attempt_acceleration &&
                midpoint_acceleration_mode_ == MIDPOINT_ACCELERATION_AITKEN &&
                have_acceleration_residual) {
                long double local_pair[2] = {0.0L, 0.0L};
                for (int iface = 0; iface < nxl; ++iface) {
                    const long double delta =
                        static_cast<long double>(acceleration_re_current_[iface]) -
                        static_cast<long double>(acceleration_re_previous_[iface]);
                    local_pair[0] += static_cast<long double>(
                        acceleration_re_previous_[iface]) * delta;
                    local_pair[1] += delta * delta;
                }
                double pair[2] = {static_cast<double>(local_pair[0]),
                                  static_cast<double>(local_pair[1])};
                MPI_Allreduce(MPI_IN_PLACE, pair, 2, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);
                const double denominator_floor = 128.0 *
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::fabs(pair[0]));
                double aitken_omega = std::numeric_limits<double>::quiet_NaN();
                if (std::isfinite(pair[0]) && std::isfinite(pair[1]) &&
                    pair[1] > denominator_floor) {
                    aitken_omega = -aitken_previous_omega * pair[0] / pair[1];
                }
                if (std::isfinite(aitken_omega) &&
                    aitken_omega >= 0.1 && aitken_omega <= 1.5) {
                    for (int iface = 0; iface < nxl; ++iface) {
                        e_end[iface] += aitken_omega *
                            acceleration_re_current_[iface];
                    }
                    accelerated_next_endpoint = true;
                    next_acceleration_coefficient = aitken_omega;
                    aitken_previous_omega = aitken_omega;
                } else {
                    coefficient_rejected = true;
                    ++result.acceleration_rejected_coefficient;
                    ++result.acceleration_history_resets;
                    have_acceleration_residual = false;
                    aitken_previous_omega = omega;
                }
            } else if (can_attempt_acceleration &&
                       midpoint_acceleration_mode_ == MIDPOINT_ACCELERATION_ANDERSON) {
                const int depth = anderson_depth_;
                const int old_count = anderson_history_count;
                for (int h = std::min(depth - 1, old_count); h > 0; --h) {
                    acceleration_e_history_[h] = acceleration_e_history_[h - 1];
                    acceleration_re_history_[h] = acceleration_re_history_[h - 1];
                }
                acceleration_e_history_[0] = e_end;
                acceleration_re_history_[0] = acceleration_re_current_;
                anderson_history_count = std::min(depth, old_count + 1);
                const int columns = anderson_history_count - 1;
                if (columns > 0) {
                    double local_field_scale = 1.0;
                    for (int iface = 0; iface < nxl; ++iface) {
                        local_field_scale = std::max(local_field_scale,
                            std::fabs(e_end[iface]));
                    }
                    double field_scale = local_field_scale;
                    MPI_Allreduce(MPI_IN_PLACE, &field_scale, 1, MPI_DOUBLE,
                                  MPI_MAX, MPI_COMM_WORLD);
                    const double inv_scale2 = 1.0 / (field_scale * field_scale);
                    long double local_ls[5] = {0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
                    for (int iface = 0; iface < nxl; ++iface) {
                        const long double d0 = static_cast<long double>(
                            acceleration_re_history_[0][iface] -
                            acceleration_re_history_[1][iface]);
                        local_ls[0] += d0 * d0;
                        local_ls[2] += d0 * static_cast<long double>(
                            acceleration_re_history_[0][iface]);
                        if (columns == 2) {
                            const long double d1 = static_cast<long double>(
                                acceleration_re_history_[1][iface] -
                                acceleration_re_history_[2][iface]);
                            local_ls[1] += d0 * d1;
                            local_ls[3] += d1 * d1;
                            local_ls[4] += d1 * static_cast<long double>(
                                acceleration_re_history_[0][iface]);
                        }
                    }
                    double ls[5] = {static_cast<double>(local_ls[0]) * inv_scale2,
                                    static_cast<double>(local_ls[1]) * inv_scale2,
                                    static_cast<double>(local_ls[2]) * inv_scale2,
                                    static_cast<double>(local_ls[3]) * inv_scale2,
                                    static_cast<double>(local_ls[4]) * inv_scale2};
                    MPI_Allreduce(MPI_IN_PLACE, ls, 5, MPI_DOUBLE, MPI_SUM,
                                  MPI_COMM_WORLD);
                    const double regularization = 1.0e-12 *
                        std::max(1.0, (ls[0] + (columns == 2 ? ls[3] : 0.0)) /
                            static_cast<double>(columns));
                    double gamma0 = 0.0, gamma1 = 0.0;
                    bool valid_system = std::isfinite(regularization);
                    if (columns == 1) {
                        const double a00 = ls[0] + regularization;
                        gamma0 = ls[2] / a00;
                        valid_system = valid_system && std::isfinite(a00) &&
                            a00 > 0.0 && std::isfinite(gamma0);
                    } else {
                        const double a00 = ls[0] + regularization;
                        const double a11 = ls[3] + regularization;
                        const double determinant = a00 * a11 - ls[1] * ls[1];
                        const double trace = a00 + a11;
                        const double discriminant = std::max(0.0,
                            trace * trace - 4.0 * determinant);
                        const double eigen_min = 0.5 * (trace - std::sqrt(discriminant));
                        const double eigen_max = 0.5 * (trace + std::sqrt(discriminant));
                        gamma0 = (ls[2] * a11 - ls[1] * ls[4]) / determinant;
                        gamma1 = (a00 * ls[4] - ls[1] * ls[2]) / determinant;
                        valid_system = valid_system && std::isfinite(determinant) &&
                            determinant > 0.0 && eigen_min > 0.0 &&
                            eigen_max / eigen_min <= 1.0e10 &&
                            std::isfinite(gamma0) && std::isfinite(gamma1);
                    }
                    const double coefficient_max = std::max(std::fabs(gamma0),
                        columns == 2 ? std::fabs(gamma1) : 0.0);
                    if (valid_system &&
                        coefficient_max <= acceleration_max_coefficient_) {
                        for (int iface = 0; iface < nxl; ++iface) {
                            double candidate = acceleration_fallback_e_[iface];
                            candidate -= gamma0 * ((acceleration_e_history_[0][iface] -
                                acceleration_e_history_[1][iface]) + omega *
                                (acceleration_re_history_[0][iface] -
                                 acceleration_re_history_[1][iface]));
                            if (columns == 2) {
                                candidate -= gamma1 * ((acceleration_e_history_[1][iface] -
                                    acceleration_e_history_[2][iface]) + omega *
                                    (acceleration_re_history_[1][iface] -
                                     acceleration_re_history_[2][iface]));
                            }
                            e_end[iface] = candidate;
                        }
                        accelerated_next_endpoint = true;
                        next_acceleration_coefficient = coefficient_max;
                    } else {
                        coefficient_rejected = true;
                        ++result.acceleration_rejected_coefficient;
                        ++result.acceleration_history_resets;
                        anderson_history_count = 0;
                        aitken_previous_omega = omega;
                    }
                }
            }

            if (!accelerated_next_endpoint) {
                e_end = acceleration_fallback_e_;
            }
            if (!coefficient_rejected) {
                acceleration_re_previous_ = acceleration_re_current_;
                have_acceleration_residual = true;
            } else {
                have_acceleration_residual = false;
            }
            if (accelerated_next_endpoint) {
                ++result.acceleration_attempts;
                acceleration_trial_pending = true;
                acceleration_residual_before = result.nonlinear_residual;
                acceleration_trial_coefficient = next_acceleration_coefficient;
            }
            if (midpoint_iteration_trace_for_test_) {
                result.midpoint_acceleration_omega_history.back() =
                    accelerated_next_endpoint ? next_acceleration_coefficient : omega;
                result.midpoint_acceleration_residual_before_history.back() =
                    result.nonlinear_residual;
                result.midpoint_acceleration_status_history.back() =
                    coefficient_rejected ? 6 :
                    (accelerated_next_endpoint ? 1 :
                     (fallback_plain_evaluation ? 5 : 0));
            }
            fallback_plain_evaluation = false;
        } else {
            for (int iface = 0; iface < nxl; ++iface) {
                e_end[iface] = (1.0 - omega) * e_end[iface] +
                    omega * next_e[iface];
            }
        }
        for (int iface = 0; iface < nxl; ++iface) {
            previous_j[iface] = result.j_bkg_face_mid[iface];
        }
        have_previous = true;
        guess = work;
        if (iter + 1 == max_iters) {
            if (!transport_safe_step_gate(last_work, nsub))
                return result;
            if (detailed_operator_diagnostics)
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
    trace_live("soft_accept", max_iters, -1);
    return result;
}
