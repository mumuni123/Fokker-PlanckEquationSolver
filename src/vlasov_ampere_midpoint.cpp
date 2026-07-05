#include "vlasov_ampere_midpoint.h"

#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <omp.h>

namespace {

const double NEG_TOL_SOFT   = 1.0e-8;
const double NEG_TOL_MEDIUM = 1.0e-7;
const double NEG_TOL_HARD   = 1.0e300;  // disabled for long-run test
const int BKG_STAGE_COUNT = 3;
const int BKG_DIV_COMPONENT_COUNT = 3;
const double CORE_DIAG_BOUNDARY_WIDTH = 0.2 * Const::micro;
const int LOW_U_MU_LIMIT_COUNT = 6;
const int LOW_U_ADAPTIVE_MAX_SUBCYCLES = 4;
const double LOW_U_ADAPTIVE_CFL_THRESHOLD = 0.25;
const double LOW_U_MU_ALPHA_SAFETY = 0.95;
const double LOW_U_MU_ALPHA_SMOOTH = 0.25;
const double LOW_U_MU_ENDPOINT_ALPHA_CAP = 0.5;
// 7.1.2: safety factor for FCT face-alpha to prevent limiter overshoot
const double X_FCT_SAFETY = 0.99;
const double X_LOW_ABS_F_TOL = 1.0e-5;
const double X_LOW_REL_F_TOL = 1.0e-10;
const double X_LOW_CFL_TOL = 1.0e-12;
const double X_LOW_INPUT_DEBT_REL_TOL = 1.0e-8;
const double X_LOW_INPUT_DEBT_NEG_MASS_FRAC = 1.0e-12;
const double X_LOW_INPUT_DEBT_CELL_FRAC = 1.0e-6;
const double X_LOW_INPUT_DEBT_OUTPUT_REL_TOL = 1.0e-8;

enum XLowFailureKind {
    X_LOW_OK = 0,
    X_LOW_INPUT_DEBT = 1,
    X_LOW_INPUT_BAD = 2,
    X_LOW_DONOR_BUG = 3,
    X_LOW_TRUE_CFL = 4
};

void resize_or_zero(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    } else {
        std::fill(values.begin(), values.end(), 0.0);
    }
}

void resize_without_fill(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    }
}

void fill_periodic_ghosts_single_rank(Species& sp, const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t slice = Param::Nvmu;
    if (nxl <= 0) return;

    for (int g = 0; g < ng; ++g) {
        const int left_src = ng + ((nxl - ng + g) % nxl);
        const int right_src = ng + (g % nxl);
        std::memcpy(&sp.f[static_cast<size_t>(g) * slice],
                    &sp.f[static_cast<size_t>(left_src) * slice],
                    slice * sizeof(double));
        std::memcpy(&sp.f[static_cast<size_t>(ng + nxl + g) * slice],
                    &sp.f[static_cast<size_t>(right_src) * slice],
                    slice * sizeof(double));
    }
}

void close_periodic_face_blocks(std::vector<double>& face_values,
                                int nxl,
                                int block,
                                int mpi_rank,
                                int mpi_size,
                                int tag)
{
    if (face_values.size() <
        static_cast<size_t>(nxl + 1) * static_cast<size_t>(block)) {
        return;
    }
    const size_t right = static_cast<size_t>(nxl) * block;
    if (mpi_size <= 1) {
        std::copy(face_values.begin(),
                  face_values.begin() + block,
                  face_values.begin() + right);
        return;
    }

    const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right_peer = (mpi_rank + 1) % mpi_size;
    MPI_Sendrecv(face_values.data(), block, MPI_DOUBLE, left_peer, tag,
                 face_values.data() + right, block, MPI_DOUBLE,
                 right_peer, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

inline size_t u_xface_index(int iface, int face, int imu)
{
    return (static_cast<size_t>(iface) * (Param::Nv + 1)
          + static_cast<size_t>(face)) * Param::Nmu
         + static_cast<size_t>(imu);
}

inline size_t u_cell_index(int ix, int face, int imu)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1)
          + static_cast<size_t>(face)) * Param::Nmu
         + static_cast<size_t>(imu);
}

inline size_t mu_xface_index(int iface, int iv, int face)
{
    return (static_cast<size_t>(iface) * Param::Nv
          + static_cast<size_t>(iv)) * (Param::Nmu + 1)
         + static_cast<size_t>(face);
}

double energy_compatible_chain_speed(const Species& sp, int iv)
{
    const double shell = sp.vgrid.moment_weight[iv];
    if (!(shell > 0.0)) return sp.vgrid.chain_speed_cells[iv];

    double numerator = 0.0;
    if (iv > 0) {
        const double dK_left =
            (sp.vgrid.gamma_cells[iv] - sp.vgrid.gamma_cells[iv - 1])
          * sp.mass * Const::c * Const::c;
        numerator += sp.vgrid.v2_faces[iv] * dK_left;
    }
    if (iv + 1 < Param::Nv) {
        const double dK_right =
            (sp.vgrid.gamma_cells[iv + 1] - sp.vgrid.gamma_cells[iv])
          * sp.mass * Const::c * Const::c;
        numerator += sp.vgrid.v2_faces[iv + 1] * dK_right;
    }

    const double weighted =
        Const::pi * sp.vgrid.dmu * numerator / (sp.mass * Const::c);
    return weighted / shell;
}

double shell_consistent_mu_u_eff(const VelocityGrid& vg, int iv)
{
    const double u_left = vg.v_faces[iv];
    const double u_right = vg.v_faces[iv + 1];
    const double int_u2 =
        (u_right * u_right * u_right - u_left * u_left * u_left) / 3.0;
    const double int_u =
        (u_right * u_right - u_left * u_left) / 2.0;
    if (int_u > 0.0 && int_u2 > 0.0) {
        return std::max(int_u2 / int_u, Param::u_floor);
    }
    return std::max(vg.v_cells[iv], Param::u_floor);
}

double smooth_low_u_mu_alpha(double alpha)
{
    alpha = std::max(0.0, std::min(1.0, alpha));
    if (alpha >= 1.0) return 1.0;
    const double safe = LOW_U_MU_ALPHA_SAFETY * alpha;
    return safe * (1.0 - LOW_U_MU_ALPHA_SMOOTH * (1.0 - safe));
}

int adaptive_low_u_subcycles(double out, double available)
{
    if (!(out > 0.0) || !std::isfinite(out)) return 1;
    const double cfl =
        out / std::max(available, std::numeric_limits<double>::min());
    if (cfl <= LOW_U_ADAPTIVE_CFL_THRESHOLD) return 1;
    const int cycles =
        static_cast<int>(std::ceil(cfl / LOW_U_ADAPTIVE_CFL_THRESHOLD));
    return std::max(1, std::min(LOW_U_ADAPTIVE_MAX_SUBCYCLES, cycles));
}

double max_abs_vector(const std::vector<double>& values, size_t n)
{
    double result = 0.0;
    for (size_t i = 0; i < n && i < values.size(); ++i) {
        result = std::max(result, std::fabs(values[i]));
    }
    return result;
}

void sync_cell_current_from_faces(Species& sp, int nxl)
{
    if (sp.current_x.size() != static_cast<size_t>(nxl)) {
        sp.current_x.assign(static_cast<size_t>(nxl), 0.0);
    }
    if (sp.current_face_x.size() < static_cast<size_t>(nxl + 1)) {
        return;
    }
    for (int ix = 0; ix < nxl; ++ix) {
        sp.current_x[static_cast<size_t>(ix)] =
            0.5 * (sp.current_face_x[static_cast<size_t>(ix)]
                 + sp.current_face_x[static_cast<size_t>(ix + 1)]);
    }
}

struct FiniteFluxCandidate {
    double severity;
    double max_negative;
    double relative_negative;
    double updated;
    double f0;
    double dx_div;
    double du_div;
    double dmu_div;
    int ix;
    int iv;
    int imu;
    int has_failure;
};

struct UFluxAuditCandidate {
    double severity;
    double f0;
    double f_low;
    double f_high;
    double alpha;
    double du_div_low;
    double du_div_high;
    int ix;
    int iv;
    int imu;
    int valid;
};

FiniteFluxCandidate empty_finite_flux_candidate()
{
    FiniteFluxCandidate info;
    info.severity = 0.0;
    info.max_negative = 0.0;
    info.relative_negative = 0.0;
    info.updated = 0.0;
    info.f0 = 0.0;
    info.dx_div = 0.0;
    info.du_div = 0.0;
    info.dmu_div = 0.0;
    info.ix = -1;
    info.iv = -1;
    info.imu = -1;
    info.has_failure = 0;
    return info;
}

}

void VlasovAmpereMidpointSolver::exchange_ghosts_x_persistent(
    Species& sp, const SpatialGrid& sg, int mpi_rank, int mpi_size) const
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t slice = Param::Nvmu;
    const size_t buffer_size = static_cast<size_t>(ng) * slice;
    if (mpi_size <= 1) {
        fill_periodic_ghosts_single_rank(sp, sg);
        return;
    }

    if (ghost_send_left_.size() != buffer_size) {
        ghost_send_left_.resize(buffer_size);
        ghost_send_right_.resize(buffer_size);
        ghost_recv_left_.resize(buffer_size);
        ghost_recv_right_.resize(buffer_size);
    }
    const int left_rank = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right_rank = (mpi_rank + 1) % mpi_size;
    std::memcpy(ghost_send_left_.data(),
                &sp.f[static_cast<size_t>(ng) * slice],
                buffer_size * sizeof(double));
    std::memcpy(ghost_send_right_.data(),
                &sp.f[static_cast<size_t>(ng + nxl - ng) * slice],
                buffer_size * sizeof(double));

    MPI_Request reqs[4];
    int nreq = 0;
    MPI_Isend(ghost_send_left_.data(), static_cast<int>(buffer_size),
              MPI_DOUBLE, left_rank, 501, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(ghost_recv_left_.data(), static_cast<int>(buffer_size),
              MPI_DOUBLE, left_rank, 502, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Isend(ghost_send_right_.data(), static_cast<int>(buffer_size),
              MPI_DOUBLE, right_rank, 502, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(ghost_recv_right_.data(), static_cast<int>(buffer_size),
              MPI_DOUBLE, right_rank, 501, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    std::memcpy(&sp.f[0], ghost_recv_left_.data(),
                buffer_size * sizeof(double));
    std::memcpy(&sp.f[static_cast<size_t>(ng + nxl) * slice],
                ghost_recv_right_.data(), buffer_size * sizeof(double));
}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_with_fixed_substeps(
    const Species& bkg_n,
    const BeamPIC& beam_n,
    const EMFields& fields_n,
    const SpatialGrid& sg,
    double dt,
    double time,
    int mpi_rank,
    int mpi_size,
    int substeps) const
{
    Result combined;
    reset_result(combined);
    combined.substeps_used = substeps;

    Species bkg_state = bkg_n;
    BeamPIC beam_state = beam_n;
    EMFields fields_state = fields_n;
    const double sub_dt = dt / static_cast<double>(substeps);
    const double t_start = time - dt;

    for (int isub = 0; isub < substeps; ++isub) {
        const double sub_time =
            t_start + static_cast<double>(isub + 1) * sub_dt;
        Result sub = advance_single_step(bkg_state, beam_state, fields_state,
                                         sg, sub_dt, sub_time, mpi_rank,
                                         mpi_size, substeps);
        if (!sub.converged || sub.failed) {
            sub.substeps_used = substeps;
            return sub;
        }

        if (isub == 0) {
            combined = sub;
            combined.substeps_used = substeps;
        } else {
            combined.species_np1 = sub.species_np1;
            combined.beam_np1 = sub.beam_np1;
            combined.fields_np1 = sub.fields_np1;
            combined.j_bkg_face_mid = sub.j_bkg_face_mid;
            combined.j_beam_face_mid = sub.j_beam_face_mid;
            combined.j_total_face_mid = sub.j_total_face_mid;
            combined.j_bkg_energy_debug_face =
                sub.j_bkg_energy_debug_face;

            combined.delta_ke_bkg += sub.delta_ke_bkg;
            combined.delta_ke_beam += sub.delta_ke_beam;
            combined.field_work_bkg += sub.field_work_bkg;
            combined.field_work_beam += sub.field_work_beam;
            combined.energy_residual_bkg += sub.energy_residual_bkg;
            combined.continuity_residual_bkg =
                std::max(combined.continuity_residual_bkg,
                         sub.continuity_residual_bkg);
            combined.beam_continuity_residual =
                std::max(combined.beam_continuity_residual,
                         sub.beam_continuity_residual);
            combined.nonlinear_residual =
                std::max(combined.nonlinear_residual,
                         sub.nonlinear_residual);
            combined.residual_E =
                std::max(combined.residual_E, sub.residual_E);
            combined.residual_f =
                std::max(combined.residual_f, sub.residual_f);
            combined.residual_J_bkg =
                std::max(combined.residual_J_bkg, sub.residual_J_bkg);
            combined.residual_J_beam =
                std::max(combined.residual_J_beam, sub.residual_J_beam);

            const double inv_count =
                1.0 / static_cast<double>(isub + 1);
            combined.limiter_active_fraction =
                (combined.limiter_active_fraction * isub
               + sub.limiter_active_fraction) * inv_count;
            combined.limiter_active_fraction_core =
                (combined.limiter_active_fraction_core * isub
               + sub.limiter_active_fraction_core) * inv_count;
            combined.limiter_active_fraction_boundary =
                (combined.limiter_active_fraction_boundary * isub
               + sub.limiter_active_fraction_boundary) * inv_count;
            combined.limiter_min_alpha =
                std::min(combined.limiter_min_alpha,
                         sub.limiter_min_alpha);
            combined.limiter_min_alpha_core =
                std::min(combined.limiter_min_alpha_core,
                         sub.limiter_min_alpha_core);
            combined.limiter_min_alpha_boundary =
                std::min(combined.limiter_min_alpha_boundary,
                         sub.limiter_min_alpha_boundary);
            combined.limiter_energy_defect += sub.limiter_energy_defect;
            combined.limiter_mass_defect += sub.limiter_mass_defect;
            combined.limiter_momentum_defect +=
                sub.limiter_momentum_defect;
            combined.x_limiter_energy_defect +=
                sub.x_limiter_energy_defect;
            combined.x_limiter_mass_defect +=
                sub.x_limiter_mass_defect;
            combined.mu_low_u_alpha_min =
                std::min(combined.mu_low_u_alpha_min,
                         sub.mu_low_u_alpha_min);
            combined.mu_low_u_limiter_active_fraction =
                std::max(combined.mu_low_u_limiter_active_fraction,
                         sub.mu_low_u_limiter_active_fraction);
            combined.mu_low_u_energy_delta += sub.mu_low_u_energy_delta;
            combined.mu_low_u_alpha_min_boundary =
                std::min(combined.mu_low_u_alpha_min_boundary,
                         sub.mu_low_u_alpha_min_boundary);
            combined.mu_low_u_alpha_min_core =
                std::min(combined.mu_low_u_alpha_min_core,
                         sub.mu_low_u_alpha_min_core);
            combined.mu_low_u_limiter_active_fraction_boundary =
                std::max(
                    combined.mu_low_u_limiter_active_fraction_boundary,
                    sub.mu_low_u_limiter_active_fraction_boundary);
            combined.mu_low_u_limiter_active_fraction_core =
                std::max(combined.mu_low_u_limiter_active_fraction_core,
                         sub.mu_low_u_limiter_active_fraction_core);
            combined.mu_low_u_energy_delta_boundary +=
                sub.mu_low_u_energy_delta_boundary;
            combined.mu_low_u_energy_delta_core +=
                sub.mu_low_u_energy_delta_core;
            combined.mu_low_u_u_eff0 = sub.mu_low_u_u_eff0;
            combined.mu_low_u_moment_weight0 = sub.mu_low_u_moment_weight0;
            combined.mu_low_u_mu_flux_scale0 =
                sub.mu_low_u_mu_flux_scale0;
            combined.mu_low_u_half_dt_inv_shell0 =
                sub.mu_low_u_half_dt_inv_shell0;
            combined.mu_low_u_dimless_scale0 =
                std::max(combined.mu_low_u_dimless_scale0,
                         sub.mu_low_u_dimless_scale0);
            combined.mu_low_u_endpoint_flux_max =
                std::max(combined.mu_low_u_endpoint_flux_max,
                         sub.mu_low_u_endpoint_flux_max);
            combined.remap_active_fraction =
                std::max(combined.remap_active_fraction,
                         sub.remap_active_fraction);
            combined.remap_cell_count += sub.remap_cell_count;
            combined.low_u_subcycle_active_fraction =
                std::max(combined.low_u_subcycle_active_fraction,
                         sub.low_u_subcycle_active_fraction);
            combined.low_u_average_subcycles =
                std::max(combined.low_u_average_subcycles,
                         sub.low_u_average_subcycles);
            combined.low_u_max_subcycles =
                std::max(combined.low_u_max_subcycles,
                         sub.low_u_max_subcycles);
            combined.x_low_order_failed_count +=
                sub.x_low_order_failed_count;
            combined.x_low_input_min_f =
                std::min(combined.x_low_input_min_f,
                         sub.x_low_input_min_f);
            combined.x_low_max_cfl =
                std::max(combined.x_low_max_cfl, sub.x_low_max_cfl);
            combined.x_low_output_min_f =
                std::min(combined.x_low_output_min_f,
                         sub.x_low_output_min_f);
            combined.x_low_failed_count += sub.x_low_failed_count;
            combined.x_low_input_neg_mass += sub.x_low_input_neg_mass;
            combined.x_low_input_rel_neg =
                std::max(combined.x_low_input_rel_neg,
                         sub.x_low_input_rel_neg);
            combined.x_low_output_rel_neg =
                std::max(combined.x_low_output_rel_neg,
                         sub.x_low_output_rel_neg);
            combined.x_low_input_core_failed_count +=
                sub.x_low_input_core_failed_count;
            combined.x_low_input_debt_accepted =
                std::max(combined.x_low_input_debt_accepted,
                         sub.x_low_input_debt_accepted);
            combined.x_low_failure_kind =
                std::max(combined.x_low_failure_kind,
                         sub.x_low_failure_kind);
            for (int ir = 0; ir < 2; ++ir) {
                combined.region_u_limiter_energy_boundary[ir] +=
                    sub.region_u_limiter_energy_boundary[ir];
                combined.region_u_limiter_energy_core[ir] +=
                    sub.region_u_limiter_energy_core[ir];
                combined.region_abs_u_limiter_energy_boundary[ir] +=
                    sub.region_abs_u_limiter_energy_boundary[ir];
                combined.region_abs_u_limiter_energy_core[ir] +=
                    sub.region_abs_u_limiter_energy_core[ir];
                combined.region_limiter_active_fraction_boundary[ir] =
                    std::max(
                        combined.region_limiter_active_fraction_boundary[ir],
                        sub.region_limiter_active_fraction_boundary[ir]);
                combined.region_limiter_active_fraction_core[ir] =
                    std::max(
                        combined.region_limiter_active_fraction_core[ir],
                        sub.region_limiter_active_fraction_core[ir]);
            }
            for (int istage = 0; istage < BKG_STAGE_COUNT; ++istage) {
                const size_t s = static_cast<size_t>(istage);
                combined.stage_min_f[s] =
                    std::min(combined.stage_min_f[s],
                             sub.stage_min_f[s]);
                combined.stage_neg_mass[s] += sub.stage_neg_mass[s];
                combined.stage_neg_cell_count[s] +=
                    sub.stage_neg_cell_count[s];
                combined.stage_low_u_neg_mass[s] +=
                    sub.stage_low_u_neg_mass[s];
                combined.stage_core_low_u_min_f[s] =
                    std::min(combined.stage_core_low_u_min_f[s],
                             sub.stage_core_low_u_min_f[s]);
            }
            for (int ic = 0; ic < BKG_DIV_COMPONENT_COUNT; ++ic) {
                const size_t c = static_cast<size_t>(ic);
                combined.low_u_neg_added_by_div[c] +=
                    sub.low_u_neg_added_by_div[c];
            }
            for (size_t i = 0; i < combined.stage_min_f_core_by_u.size();
                 ++i) {
                combined.stage_min_f_core_by_u[i] =
                    std::min(combined.stage_min_f_core_by_u[i],
                             sub.stage_min_f_core_by_u[i]);
                combined.stage_neg_mass_core_by_u[i] +=
                    sub.stage_neg_mass_core_by_u[i];
                combined.stage_neg_cell_count_core_by_u[i] +=
                    sub.stage_neg_cell_count_core_by_u[i];
                combined.stage_min_f_boundary_by_u[i] =
                    std::min(combined.stage_min_f_boundary_by_u[i],
                             sub.stage_min_f_boundary_by_u[i]);
                combined.stage_neg_mass_boundary_by_u[i] +=
                    sub.stage_neg_mass_boundary_by_u[i];
                combined.stage_neg_cell_count_boundary_by_u[i] +=
                    sub.stage_neg_cell_count_boundary_by_u[i];
            }
            combined.x_negative_mass_before_repair +=
                sub.x_negative_mass_before_repair;
            combined.x_mass_added_by_positivity_repair +=
                sub.x_mass_added_by_positivity_repair;
            combined.positivity_energy_defect +=
                sub.positivity_energy_defect;
            combined.positivity_mass_defect +=
                sub.positivity_mass_defect;
            combined.u_force_alpha_min =
                std::min(combined.u_force_alpha_min,
                         sub.u_force_alpha_min);
            combined.u_force_alpha_active_frac =
                (combined.u_force_alpha_active_frac * isub
               + sub.u_force_alpha_active_frac) * inv_count;
            if (sub.u_flux_audit_valid &&
                (!combined.u_flux_audit_valid ||
                 sub.u_flux_audit_severity >
                     combined.u_flux_audit_severity)) {
                combined.u_flux_audit_valid = sub.u_flux_audit_valid;
                combined.u_flux_audit_rank = sub.u_flux_audit_rank;
                combined.u_flux_audit_ix = sub.u_flux_audit_ix;
                combined.u_flux_audit_iv = sub.u_flux_audit_iv;
                combined.u_flux_audit_imu = sub.u_flux_audit_imu;
                combined.u_flux_audit_severity =
                    sub.u_flux_audit_severity;
                combined.u_flux_audit_f0 = sub.u_flux_audit_f0;
                combined.u_flux_audit_f_low = sub.u_flux_audit_f_low;
                combined.u_flux_audit_f_high = sub.u_flux_audit_f_high;
                combined.u_flux_audit_alpha = sub.u_flux_audit_alpha;
                combined.u_flux_audit_du_div_low =
                    sub.u_flux_audit_du_div_low;
                combined.u_flux_audit_du_div_high =
                    sub.u_flux_audit_du_div_high;
                combined.u_flux_audit_du_div_final =
                    sub.u_flux_audit_du_div_final;
                combined.u_flux_audit_updated =
                    sub.u_flux_audit_updated;
                combined.u_flux_audit_final_xl_lo =
                    sub.u_flux_audit_final_xl_lo;
                combined.u_flux_audit_final_xl_hi =
                    sub.u_flux_audit_final_xl_hi;
                combined.u_flux_audit_final_xr_lo =
                    sub.u_flux_audit_final_xr_lo;
                combined.u_flux_audit_final_xr_hi =
                    sub.u_flux_audit_final_xr_hi;
            }
            if (sub.u_low_failure_audit.valid &&
                (!combined.u_low_failure_audit.valid ||
                 sub.u_low_failure_audit.severity >
                     combined.u_low_failure_audit.severity)) {
                combined.u_low_failure_audit =
                    sub.u_low_failure_audit;
            }
            if (sub.mu_low_failure_audit.valid &&
                (!combined.mu_low_failure_audit.valid ||
                 sub.mu_low_failure_audit.severity >
                     combined.mu_low_failure_audit.severity)) {
                combined.mu_low_failure_audit =
                    sub.mu_low_failure_audit;
            }
            combined.f_neg_min =
                std::min(combined.f_neg_min, sub.f_neg_min);
            if (sub.f_neg_ratio_max > combined.f_neg_ratio_max) {
                combined.f_neg_ratio_max = sub.f_neg_ratio_max;
                combined.f_neg_ix = sub.f_neg_ix;
                combined.f_neg_iv = sub.f_neg_iv;
                combined.f_neg_imu = sub.f_neg_imu;
            }
            combined.f_neg_mass_total += sub.f_neg_mass_total;
            combined.f_neg_cell_count += sub.f_neg_cell_count;
            combined.nonlinear_iterations += sub.nonlinear_iterations;
            combined.soft_accepted =
                combined.soft_accepted || sub.soft_accepted;
            combined.protected_converged =
                combined.protected_converged && sub.protected_converged;

            combined.current_diag.residual_if_charge +=
                sub.current_diag.residual_if_charge;
            combined.current_diag.residual_if_ampere +=
                sub.current_diag.residual_if_ampere;
            combined.current_diag.e_dot_j_charge +=
                sub.current_diag.e_dot_j_charge;
            combined.current_diag.e_dot_j_energy +=
                sub.current_diag.e_dot_j_energy;
            combined.current_diag.e_dot_j_ampere +=
                sub.current_diag.e_dot_j_ampere;
            combined.current_diag.max_abs_j_charge =
                std::max(combined.current_diag.max_abs_j_charge,
                         sub.current_diag.max_abs_j_charge);
            combined.current_diag.max_abs_j_energy =
                std::max(combined.current_diag.max_abs_j_energy,
                         sub.current_diag.max_abs_j_energy);
            combined.current_diag.max_abs_j_ampere =
                std::max(combined.current_diag.max_abs_j_ampere,
                         sub.current_diag.max_abs_j_ampere);
            combined.current_diag.max_abs_j_charge_minus_ampere =
                std::max(combined.current_diag.max_abs_j_charge_minus_ampere,
                         sub.current_diag.max_abs_j_charge_minus_ampere);
            combined.current_diag.max_abs_j_energy_minus_ampere =
                std::max(combined.current_diag.max_abs_j_energy_minus_ampere,
                         sub.current_diag.max_abs_j_energy_minus_ampere);
        }

        bkg_state = sub.species_np1;
        beam_state = sub.beam_np1;
        fields_state = sub.fields_np1;
    }

    combined.converged = true;
    combined.failed = false;
    combined.substeps_used = substeps;
    return combined;
}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_background_and_fields(
    const Species& bkg_n,
    const BeamPIC& beam_n,
    const EMFields& fields_n,
    const SpatialGrid& sg,
    double dt,
    double time,
    int mpi_rank,
    int mpi_size)
{
    const int substep_trials[7] = {1, 2, 4, 8, 16, 32, 64};

    for (int trial = 0; trial < 7; ++trial) {
        const int substeps = substep_trials[trial];
        Result result =
            (substeps == 1)
            ? advance_single_step(bkg_n, beam_n, fields_n, sg, dt, time,
                                  mpi_rank, mpi_size, 1)
            : advance_with_fixed_substeps(bkg_n, beam_n, fields_n, sg, dt,
                                          time, mpi_rank, mpi_size,
                                          substeps);
        if (result.converged && !result.failed) {
            return result;
        }
        if (result.x_low_failure_kind == X_LOW_INPUT_BAD ||
            result.x_low_failure_kind == X_LOW_DONOR_BUG) {
            return result;
        }
        if (trial == 6) {
            result.converged = false;
            result.failed = true;
            return result;
        }
        if (mpi_rank == 0 && trial < 6) {
            std::fprintf(stderr,
                         "WARNING: energy-compatible midpoint solve failed "
                         "with %d substep(s) at t = %.6e s "
                         "(residual %.6e, f %.6e, field %.6e); retrying "
                         "with %d substeps.\n",
                         substeps, time, result.nonlinear_residual,
                         result.residual_f, result.residual_E,
                         substep_trials[trial + 1]);
        }
    }

    Result failed;
    reset_result(failed);
    failed.failed = true;
    return failed;
}

UFluxAuditCandidate empty_u_flux_audit_candidate()
{
    UFluxAuditCandidate info;
    info.severity = 0.0;
    info.f0 = 0.0;
    info.f_low = 0.0;
    info.f_high = 0.0;
    info.alpha = 1.0;
    info.du_div_low = 0.0;
    info.du_div_high = 0.0;
    info.ix = -1;
    info.iv = -1;
    info.imu = -1;
    info.valid = 0;
    return info;
}

VlasovAmpereMidpointSolver::LowOrderFailureAudit
empty_low_order_failure_audit()
{
    VlasovAmpereMidpointSolver::LowOrderFailureAudit info;
    info.valid = 0;
    info.rank = -1;
    info.ix = -1;
    info.iv = -1;
    info.imu = -1;
    info.region = 0;
    info.severity = 0.0;
    info.f_input = 0.0;
    info.f_after_x = 0.0;
    info.dx_div = 0.0;
    info.dmu_div_used = 0.0;
    info.du_div_low = 0.0;
    info.f_low = 0.0;
    info.left_lower_flux = 0.0;
    info.left_upper_flux = 0.0;
    info.right_lower_flux = 0.0;
    info.right_upper_flux = 0.0;
    info.left_lower_scale = 0.0;
    info.left_upper_scale = 0.0;
    info.right_lower_scale = 0.0;
    info.right_upper_scale = 0.0;
    info.left_lower_donor_f = 0.0;
    info.left_upper_donor_f = 0.0;
    info.right_lower_donor_f = 0.0;
    info.right_upper_donor_f = 0.0;
    info.lower_characteristic = 0.0;
    info.upper_characteristic = 0.0;
    info.moment_weight = 0.0;
    info.cell_budget = 0.0;
    info.low_order_failed_count = 0.0;
    info.left_lower_donor_index = -1;
    info.left_upper_donor_index = -1;
    info.right_lower_donor_index = -1;
    info.right_upper_donor_index = -1;
    return info;
}

VlasovAmpereMidpointSolver::VlasovAmpereMidpointSolver()
    : step_diagnostics_enabled_(false)
{}

void VlasovAmpereMidpointSolver::reset_current_diag(
    CurrentDiagnostics& diag) const
{
    diag.residual_if_charge = 0.0;
    diag.residual_if_ampere = 0.0;
    diag.e_dot_j_charge = 0.0;
    diag.e_dot_j_energy = 0.0;
    diag.e_dot_j_ampere = 0.0;
    diag.max_abs_j_charge = 0.0;
    diag.max_abs_j_energy = 0.0;
    diag.max_abs_j_ampere = 0.0;
    diag.max_abs_j_charge_minus_ampere = 0.0;
    diag.max_abs_j_energy_minus_ampere = 0.0;
}

void VlasovAmpereMidpointSolver::reset_result(Result& result) const
{
    reset_current_diag(result.current_diag);
    result.delta_ke_bkg = 0.0;
    result.delta_ke_beam = 0.0;
    result.field_work_bkg = 0.0;
    result.field_work_beam = 0.0;
    result.energy_residual_bkg = 0.0;
    result.continuity_residual_bkg = 0.0;
    result.beam_continuity_residual = 0.0;
    result.nonlinear_residual = 0.0;
    result.residual_E = 0.0;
    result.residual_f = 0.0;
    result.residual_J_bkg = 0.0;
    result.residual_J_beam = 0.0;
    result.limiter_active_fraction = 0.0;
    result.limiter_min_alpha = 1.0;
    result.limiter_active_fraction_core = 0.0;
    result.limiter_active_fraction_boundary = 0.0;
    result.limiter_min_alpha_core = 1.0;
    result.limiter_min_alpha_boundary = 1.0;
    result.limiter_energy_defect = 0.0;
    result.limiter_mass_defect = 0.0;
    result.limiter_momentum_defect = 0.0;
    result.x_limiter_energy_defect = 0.0;
    result.x_limiter_mass_defect = 0.0;
    result.mu_low_u_alpha_min = 1.0;
    result.mu_low_u_limiter_active_fraction = 0.0;
    result.mu_low_u_energy_delta = 0.0;
    result.mu_low_u_alpha_min_boundary = 1.0;
    result.mu_low_u_alpha_min_core = 1.0;
    result.mu_low_u_limiter_active_fraction_boundary = 0.0;
    result.mu_low_u_limiter_active_fraction_core = 0.0;
    result.mu_low_u_energy_delta_boundary = 0.0;
    result.mu_low_u_energy_delta_core = 0.0;
    result.mu_low_u_u_eff0 = 0.0;
    result.mu_low_u_moment_weight0 = 0.0;
    result.mu_low_u_mu_flux_scale0 = 0.0;
    result.mu_low_u_half_dt_inv_shell0 = 0.0;
    result.mu_low_u_dimless_scale0 = 0.0;
    result.mu_low_u_endpoint_flux_max = 0.0;
    result.remap_active_fraction = 0.0;
    result.remap_cell_count = 0;
    result.low_u_subcycle_active_fraction = 0.0;
    result.low_u_average_subcycles = 1.0;
    result.low_u_max_subcycles = 1;
    result.x_low_order_failed_count = 0.0;
    result.x_low_input_min_f = std::numeric_limits<double>::infinity();
    result.x_low_max_cfl = 0.0;
    result.x_low_output_min_f = std::numeric_limits<double>::infinity();
    result.x_low_failed_count = 0.0;
    result.x_low_input_neg_mass = 0.0;
    result.x_low_input_rel_neg = 0.0;
    result.x_low_output_rel_neg = 0.0;
    result.x_low_input_core_failed_count = 0.0;
    result.x_low_input_debt_accepted = 0.0;
    result.x_low_failure_kind = X_LOW_OK;
    // 7.1.6: reset per-direction flux diagnostics
    for (int d = 0; d < 3; ++d) {
        FluxPositivityDiag& fp = result.flux_pos[d];
        fp.min_f_before = std::numeric_limits<double>::infinity();
        fp.min_f_low    = std::numeric_limits<double>::infinity();
        fp.min_f_final  = std::numeric_limits<double>::infinity();
        fp.low_order_failed_count  = 0.0;
        fp.alpha_active_fraction   = 0.0;
        fp.alpha_min = 1.0;
        fp.alpha_core_fraction     = 0.0;
        fp.alpha_boundary_fraction = 0.0;
        fp.negative_mass_prevented = 0.0;
        FluxDefectDiag& fd = result.flux_defect[d];
        fd.mass_defect = 0.0;
        fd.momentum_defect = 0.0;
        fd.energy_defect = 0.0;
        fd.boundary_mass_loss = 0.0;
        fd.boundary_energy_loss = 0.0;
    }
    for (int ir = 0; ir < 2; ++ir) {
        result.region_u_limiter_energy_boundary[ir] = 0.0;
        result.region_u_limiter_energy_core[ir] = 0.0;
        result.region_abs_u_limiter_energy_boundary[ir] = 0.0;
        result.region_abs_u_limiter_energy_core[ir] = 0.0;
        result.region_limiter_active_fraction_boundary[ir] = 0.0;
        result.region_limiter_active_fraction_core[ir] = 0.0;
    }
    result.stage_min_f.assign(BKG_STAGE_COUNT,
                              std::numeric_limits<double>::infinity());
    result.stage_neg_mass.assign(BKG_STAGE_COUNT, 0.0);
    result.stage_neg_cell_count.assign(BKG_STAGE_COUNT, 0);
    result.stage_low_u_neg_mass.assign(BKG_STAGE_COUNT, 0.0);
    result.low_u_neg_added_by_div.assign(BKG_DIV_COMPONENT_COUNT, 0.0);
    result.stage_core_low_u_min_f.assign(
        BKG_STAGE_COUNT, std::numeric_limits<double>::infinity());
    result.stage_min_f_core_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv,
        std::numeric_limits<double>::infinity());
    result.stage_neg_mass_core_by_u.assign(BKG_STAGE_COUNT * Param::Nv, 0.0);
    result.stage_neg_cell_count_core_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0);
    result.stage_min_f_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv,
        std::numeric_limits<double>::infinity());
    result.stage_neg_mass_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0.0);
    result.stage_neg_cell_count_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0);
    result.x_negative_mass_before_repair = 0.0;
    result.x_mass_added_by_positivity_repair = 0.0;
    result.positivity_energy_defect = 0.0;
    result.positivity_mass_defect = 0.0;
    result.u_force_alpha_min = 1.0;
    result.u_force_alpha_active_frac = 0.0;
    result.u_flux_audit_valid = 0;
    result.u_flux_audit_rank = -1;
    result.u_flux_audit_ix = -1;
    result.u_flux_audit_iv = -1;
    result.u_flux_audit_imu = -1;
    result.u_flux_audit_severity = 0.0;
    result.u_flux_audit_f0 = 0.0;
    result.u_flux_audit_f_low = 0.0;
    result.u_flux_audit_f_high = 0.0;
    result.u_flux_audit_alpha = 1.0;
    result.u_flux_audit_du_div_low = 0.0;
    result.u_flux_audit_du_div_high = 0.0;
    result.u_flux_audit_du_div_final = 0.0;
    result.u_flux_audit_updated = 0.0;
    result.u_flux_audit_final_xl_lo = 0.0;
    result.u_flux_audit_final_xl_hi = 0.0;
    result.u_flux_audit_final_xr_lo = 0.0;
    result.u_flux_audit_final_xr_hi = 0.0;
    result.u_low_failure_audit = empty_low_order_failure_audit();
    result.mu_low_failure_audit = empty_low_order_failure_audit();
    result.f_neg_min = 0.0;
    result.f_neg_ratio_max = 0.0;
    result.f_neg_mass_total = 0.0;
    result.f_neg_cell_count = 0;
    result.f_neg_ix = -1;
    result.f_neg_iv = -1;
    result.f_neg_imu = -1;
    result.nonlinear_iterations = 0;
    result.converged = false;
    result.failed = false;
    result.soft_accepted = false;
    result.protected_converged = false;
    result.substeps_used = 1;
}

void VlasovAmpereMidpointSolver::set_midpoint_field(
    EMFields& fields_mid,
    const EMFields& fields_n,
    const std::vector<double>& ex_mid,
    const SpatialGrid& sg,
    int mpi_rank,
    int mpi_size) const
{
    fields_mid = fields_n;
    const size_t n =
        std::min(static_cast<size_t>(sg.nx_local), fields_mid.Ex_face.size());
    for (size_t iface = 0; iface < n; ++iface) {
        fields_mid.Ex_face[iface] =
            (iface < ex_mid.size()) ? ex_mid[iface] : fields_n.Ex_face[iface];
    }
    if (fields_mid.Ex_face.size() > static_cast<size_t>(sg.nx_local)) {
        fields_mid.Ex_face[static_cast<size_t>(sg.nx_local)] =
            (sg.nx_local > 0) ? fields_mid.Ex_face[0] : 0.0;
    }
    fields_mid.sync_cell_ex_from_faces(mpi_rank, mpi_size);
}

void VlasovAmpereMidpointSolver::compute_midpoint_fluxes(
    const Species& bkg_n,
    const Species& bkg_guess_np1,
    const EMFields& fields_mid,
    const SpatialGrid& sg,
    double dt,
    int mpi_rank,
    int mpi_size,
    Species& bkg_new,
    FluxPack& fluxes,
    const std::vector<double>* alpha_smooth_from,
    double alpha_smooth_beta,
    bool& finite) const
{
    compute_vlasov_midpoint_residual(bkg_n, bkg_guess_np1, fields_mid, sg,
                                     dt, mpi_rank, mpi_size, bkg_new, fluxes,
                                     alpha_smooth_from, alpha_smooth_beta,
                                     finite);
}

void VlasovAmpereMidpointSolver::compute_vlasov_midpoint_residual(
    const Species& bkg_n,
    const Species& bkg_guess_np1,
    const EMFields& fields_mid,
    const SpatialGrid& sg,
    double dt,
    int mpi_rank,
    int mpi_size,
    Species& bkg_new,
    FluxPack& fluxes,
    const std::vector<double>* alpha_smooth_from,
    double alpha_smooth_beta,
    bool& finite) const
{
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    const size_t cells = static_cast<size_t>(nxl) * Param::Nvmu;
    const size_t x_faces = static_cast<size_t>(nxl + 1) * Param::Nvmu;
    finite = true;
    (void)alpha_smooth_from;
    (void)alpha_smooth_beta;

    resize_or_zero(fluxes.x_high, x_faces);
    resize_without_fill(fluxes.x_final, x_faces);
    resize_or_zero(fluxes.cell_alpha_u, cells);
    resize_or_zero(fluxes.cell_alpha_x, cells);
    resize_or_zero(fluxes.j_bkg_face, static_cast<size_t>(nxl + 1));
    // 7.1.1: allocate unified u/mu flux storage
    resize_or_zero(fluxes.u_high,
                   static_cast<size_t>(nxl + 1) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.u_low,
                   static_cast<size_t>(nxl + 1) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.u_final,
                   static_cast<size_t>(nxl + 1) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.u_high_cell,
                   static_cast<size_t>(nxl) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.u_low_cell,
                   static_cast<size_t>(nxl) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.u_final_cell,
                   static_cast<size_t>(nxl) * (Param::Nv + 1) * Param::Nmu);
    resize_or_zero(fluxes.mu_high,
                   static_cast<size_t>(nxl + 1) * Param::Nv * (Param::Nmu + 1));
    resize_or_zero(fluxes.mu_low,
                   static_cast<size_t>(nxl + 1) * Param::Nv * (Param::Nmu + 1));
    resize_or_zero(fluxes.mu_final,
                   static_cast<size_t>(nxl + 1) * Param::Nv * (Param::Nmu + 1));
    resize_or_zero(fluxes.cell_alpha_mu, cells);
    fluxes.limiter_active_fraction = 0.0;
    fluxes.limiter_min_alpha = 1.0;
    fluxes.limiter_active_fraction_core = 0.0;
    fluxes.limiter_active_fraction_boundary = 0.0;
    fluxes.limiter_min_alpha_core = 1.0;
    fluxes.limiter_min_alpha_boundary = 1.0;
    fluxes.limiter_energy_defect = 0.0;
    fluxes.limiter_mass_defect = 0.0;
    fluxes.limiter_momentum_defect = 0.0;
    fluxes.x_limiter_energy_defect = 0.0;
    fluxes.x_limiter_mass_defect = 0.0;
    fluxes.x_low_order_failed_count = 0.0;
    fluxes.x_low_input_min_f = std::numeric_limits<double>::infinity();
    fluxes.x_low_max_cfl = 0.0;
    fluxes.x_low_output_min_f = std::numeric_limits<double>::infinity();
    fluxes.x_low_failed_count = 0.0;
    fluxes.x_low_input_neg_mass = 0.0;
    fluxes.x_low_input_rel_neg = 0.0;
    fluxes.x_low_output_rel_neg = 0.0;
    fluxes.x_low_input_core_failed_count = 0.0;
    fluxes.x_low_input_debt_accepted = 0.0;
    fluxes.x_low_failure_kind = X_LOW_OK;
    // 7.1.6: zero per-direction flux diagnostics
    for (int d = 0; d < 3; ++d) {
        FluxPositivityDiag& fp = fluxes.flux_pos[d];
        fp.min_f_before = std::numeric_limits<double>::infinity();
        fp.min_f_low    = std::numeric_limits<double>::infinity();
        fp.min_f_final  = std::numeric_limits<double>::infinity();
        fp.low_order_failed_count = 0.0;
        fp.alpha_active_fraction  = 0.0;
        fp.alpha_min = 1.0;
        fp.alpha_core_fraction    = 0.0;
        fp.alpha_boundary_fraction = 0.0;
        fp.negative_mass_prevented = 0.0;
        FluxDefectDiag& fd = fluxes.flux_defect[d];
        fd.mass_defect = 0.0;
        fd.momentum_defect = 0.0;
        fd.energy_defect = 0.0;
        fd.boundary_mass_loss = 0.0;
        fd.boundary_energy_loss = 0.0;
    }
    fluxes.mu_low_u_alpha_min = 1.0;
    fluxes.mu_low_u_limiter_active_fraction = 0.0;
    fluxes.mu_low_u_energy_delta = 0.0;
    fluxes.mu_low_u_alpha_min_boundary = 1.0;
    fluxes.mu_low_u_alpha_min_core = 1.0;
    fluxes.mu_low_u_limiter_active_fraction_boundary = 0.0;
    fluxes.mu_low_u_limiter_active_fraction_core = 0.0;
    fluxes.mu_low_u_energy_delta_boundary = 0.0;
    fluxes.mu_low_u_energy_delta_core = 0.0;
    fluxes.mu_low_u_u_eff0 = 0.0;
    fluxes.mu_low_u_moment_weight0 = 0.0;
    fluxes.mu_low_u_mu_flux_scale0 = 0.0;
    fluxes.mu_low_u_half_dt_inv_shell0 = 0.0;
    fluxes.mu_low_u_dimless_scale0 = 0.0;
    fluxes.mu_low_u_endpoint_flux_max = 0.0;
    fluxes.remap_active_fraction = 0.0;
    fluxes.remap_cell_count = 0;
    fluxes.low_u_subcycle_active_fraction = 0.0;
    fluxes.low_u_average_subcycles = 1.0;
    fluxes.low_u_max_subcycles = 1;
    for (int ir = 0; ir < 2; ++ir) {
        fluxes.region_u_limiter_energy_boundary[ir] = 0.0;
        fluxes.region_u_limiter_energy_core[ir] = 0.0;
        fluxes.region_abs_u_limiter_energy_boundary[ir] = 0.0;
        fluxes.region_abs_u_limiter_energy_core[ir] = 0.0;
        fluxes.region_limiter_active_fraction_boundary[ir] = 0.0;
        fluxes.region_limiter_active_fraction_core[ir] = 0.0;
    }
    fluxes.stage_min_f.assign(BKG_STAGE_COUNT,
                              std::numeric_limits<double>::infinity());
    fluxes.stage_neg_mass.assign(BKG_STAGE_COUNT, 0.0);
    fluxes.stage_neg_cell_count.assign(BKG_STAGE_COUNT, 0);
    fluxes.stage_low_u_neg_mass.assign(BKG_STAGE_COUNT, 0.0);
    fluxes.low_u_neg_added_by_div.assign(BKG_DIV_COMPONENT_COUNT, 0.0);
    fluxes.stage_core_low_u_min_f.assign(
        BKG_STAGE_COUNT, std::numeric_limits<double>::infinity());
    fluxes.stage_min_f_core_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv,
        std::numeric_limits<double>::infinity());
    fluxes.stage_neg_mass_core_by_u.assign(BKG_STAGE_COUNT * Param::Nv, 0.0);
    fluxes.stage_neg_cell_count_core_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0);
    fluxes.stage_min_f_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv,
        std::numeric_limits<double>::infinity());
    fluxes.stage_neg_mass_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0.0);
    fluxes.stage_neg_cell_count_boundary_by_u.assign(
        BKG_STAGE_COUNT * Param::Nv, 0);
    fluxes.negative_mass_before_repair = 0.0;
    fluxes.mass_added_by_positivity_repair = 0.0;
    fluxes.positivity_energy_defect = 0.0;
    fluxes.positivity_mass_defect = 0.0;
    fluxes.u_force_alpha_min = 1.0;
    fluxes.u_force_alpha_active_frac = 0.0;
    fluxes.u_flux_audit_valid = 0;
    fluxes.u_flux_audit_rank = -1;
    fluxes.u_flux_audit_ix = -1;
    fluxes.u_flux_audit_iv = -1;
    fluxes.u_flux_audit_imu = -1;
    fluxes.u_flux_audit_severity = 0.0;
    fluxes.u_flux_audit_f0 = 0.0;
    fluxes.u_flux_audit_f_low = 0.0;
    fluxes.u_flux_audit_f_high = 0.0;
    fluxes.u_flux_audit_alpha = 1.0;
    fluxes.u_flux_audit_du_div_low = 0.0;
    fluxes.u_flux_audit_du_div_high = 0.0;
    fluxes.u_flux_audit_du_div_final = 0.0;
    fluxes.u_flux_audit_updated = 0.0;
    fluxes.u_flux_audit_final_xl_lo = 0.0;
    fluxes.u_flux_audit_final_xl_hi = 0.0;
    fluxes.u_flux_audit_final_xr_lo = 0.0;
    fluxes.u_flux_audit_final_xr_hi = 0.0;
    fluxes.u_low_failure_audit = empty_low_order_failure_audit();
    fluxes.mu_low_failure_audit = empty_low_order_failure_audit();
    fluxes.f_neg_min = 0.0;
    fluxes.f_neg_ratio_max = 0.0;
    fluxes.f_neg_mass_total = 0.0;
    fluxes.f_neg_cell_count = 0;
    fluxes.f_neg_ix = -1;
    fluxes.f_neg_iv = -1;
    fluxes.f_neg_imu = -1;
    fluxes.finite_flux_max_negative = 0.0;
    fluxes.finite_flux_relative_negative = 0.0;
    fluxes.finite_flux_updated = 0.0;
    fluxes.finite_flux_f0 = 0.0;
    fluxes.finite_flux_dx_div = 0.0;
    fluxes.finite_flux_du_div = 0.0;
    fluxes.finite_flux_dmu_div = 0.0;
    fluxes.finite_flux_rank = -1;
    fluxes.finite_flux_ix = -1;
    fluxes.finite_flux_iv = -1;
    fluxes.finite_flux_imu = -1;
    fluxes.finite_flux_has_failure = 0;

    const double dt_dx = dt / sg.dx;
    double local_negative_mass = 0.0;
    const double eps_tol_base = 1.0e-12;
    const bool collect_stage_diagnostics = step_diagnostics_enabled_;

    Species f_mid = bkg_n;
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_negative_mass)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
            const int iv = k_int / Param::Nmu;
            const size_t offset =
                static_cast<size_t>(ng + ix) * Param::Nvmu
              + static_cast<size_t>(k_int);
            f_mid.f[offset] =
                0.5 * (bkg_n.f[offset] + bkg_guess_np1.f[offset]);
            const double f0 = bkg_n.f[offset];
            if (f0 < 0.0) {
                local_negative_mass +=
                    -f0 * bkg_n.vgrid.moment_weight[iv] * sg.dx;
            }
        }
    }
    exchange_ghosts_x_persistent(f_mid, sg, mpi_rank, mpi_size);

    #pragma omp parallel for schedule(static)
    for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
        const size_t k = static_cast<size_t>(k_int);
        const int iv = k_int / Param::Nmu;
        const int imu = k_int - iv * Param::Nmu;
        const double vx_star =
            energy_compatible_chain_speed(bkg_n, iv)
          * f_mid.vgrid.mu_cells[imu];
        for (int iface = 0; iface < nxl; ++iface) {
            const int ix_left = ng + iface - 1;
            const int ix_right = ng + iface;
            const size_t left =
                static_cast<size_t>(ix_left) * Param::Nvmu + k;
            const size_t right =
                static_cast<size_t>(ix_right) * Param::Nvmu + k;
            const double f_left = f_mid.f[left];
            const double f_right = f_mid.f[right];
            const size_t face_slot =
                static_cast<size_t>(iface) * Param::Nvmu + k;
            fluxes.x_high[face_slot] = vx_star * 0.5 * (f_left + f_right);
        }
    }
    close_periodic_face_blocks(fluxes.x_high, nxl,
                               static_cast<int>(Param::Nvmu),
                               mpi_rank, mpi_size, 504);

    std::vector<double> x_low_shell_max(Param::Nv, 0.0);
    {
        const int nthreads = std::max(1, omp_get_max_threads());
        std::vector<double> thread_shell_max(
            static_cast<size_t>(nthreads) * Param::Nv, 0.0);
        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            double* shell_max =
                &thread_shell_max[static_cast<size_t>(tid) * Param::Nv];
            #pragma omp for collapse(2) schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                for (int k_int = 0;
                     k_int < static_cast<int>(Param::Nvmu); ++k_int) {
                    const int iv = k_int / Param::Nmu;
                    const size_t k = static_cast<size_t>(k_int);
                    const size_t src =
                        static_cast<size_t>(ng + ix) * Param::Nvmu + k;
                    const double f = bkg_n.f[src];
                    if (std::isfinite(f) && f > shell_max[iv]) {
                        shell_max[iv] = f;
                    }
                }
            }
        }
        for (int tid = 0; tid < nthreads; ++tid) {
            const double* shell_max =
                &thread_shell_max[static_cast<size_t>(tid) * Param::Nv];
            for (int iv = 0; iv < Param::Nv; ++iv) {
                x_low_shell_max[iv] =
                    std::max(x_low_shell_max[iv], shell_max[iv]);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, x_low_shell_max.data(),
                      Param::Nv, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    // 7.1.2: x-direction donor-cell low-order flux.
    // Use the accepted f^n state on the upwind side of each face:
    //   F_{i+1/2} = max(vx,0) f_i + min(vx,0) f_{i+1}.
    // There is deliberately no max(0,f) clip here.  If this monotone
    // low-order update fails positivity, the caller must retry with a
    // smaller true x substep.
    resize_or_zero(fluxes.x_low, x_faces);
    double local_x_low_max_cfl = 0.0;
    #pragma omp parallel for schedule(static) \
        reduction(max:local_x_low_max_cfl)
    for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
        const size_t k = static_cast<size_t>(k_int);
        const int iv = k_int / Param::Nmu;
        const int imu = k_int - iv * Param::Nmu;
        const double vx_star =
            energy_compatible_chain_speed(bkg_n, iv)
          * f_mid.vgrid.mu_cells[imu];
        local_x_low_max_cfl =
            std::max(local_x_low_max_cfl, std::fabs(vx_star) * dt_dx);
        const double vx_pos = std::max(vx_star, 0.0);
        const double vx_neg = std::min(vx_star, 0.0);
        for (int iface = 0; iface < nxl; ++iface) {
            const int ix_left = ng + iface - 1;
            const int ix_right = ng + iface;
            const size_t left =
                static_cast<size_t>(ix_left) * Param::Nvmu + k;
            const size_t right =
                static_cast<size_t>(ix_right) * Param::Nvmu + k;
            const size_t face_slot =
                static_cast<size_t>(iface) * Param::Nvmu + k;
            fluxes.x_low[face_slot] =
                vx_pos * bkg_n.f[left] + vx_neg * bkg_n.f[right];
        }
    }
    close_periodic_face_blocks(fluxes.x_low, nxl,
                               static_cast<int>(Param::Nvmu),
                               mpi_rank, mpi_size, 509);

    // 7.1.2: x-direction FCT limiter with f_floor tolerance
    // f_floor = -eps * local_scale allows roundoff-level negatives
    // while catching true CFL violations when the low-order update itself
    // produces significant negative values.
    // 7.1.6: per-direction x diagnostics (collected in this loop)
    double local_x_min_f_before = std::numeric_limits<double>::infinity();
    double local_x_min_f_low    = std::numeric_limits<double>::infinity();
    double local_x_neg_mass_prevented = 0.0;
    std::fill(fluxes.cell_alpha_x.begin(), fluxes.cell_alpha_x.end(), 1.0);
    int local_x_any_limited = 0;
    double local_x_low_order_failed = 0.0;
    double local_x_low_input_failed = 0.0;
    double local_x_low_input_neg_mass = 0.0;
    double local_x_low_input_rel_neg = 0.0;
    double local_x_low_output_rel_neg = 0.0;
    double local_x_low_input_core_failed = 0.0;
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(max:local_x_any_limited) \
        reduction(+:local_x_low_order_failed) \
        reduction(+:local_x_low_input_failed) \
        reduction(+:local_x_low_input_neg_mass) \
        reduction(+:local_x_low_input_core_failed) \
        reduction(max:local_x_low_input_rel_neg) \
        reduction(max:local_x_low_output_rel_neg) \
        reduction(min:local_x_min_f_before) \
        reduction(min:local_x_min_f_low) \
        reduction(+:local_x_neg_mass_prevented)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
            const size_t k = static_cast<size_t>(k_int);
            const size_t slot = static_cast<size_t>(ix) * Param::Nvmu + k;
            const size_t src =
                static_cast<size_t>(ng + ix) * Param::Nvmu + k;
            const double f0 = bkg_n.f[src];
            // 7.1.6: track min f before any update
            if (f0 < local_x_min_f_before) local_x_min_f_before = f0;
            const double div_low = dt_dx *
                (fluxes.x_low[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_low[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double div_high = dt_dx *
                (fluxes.x_high[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_high[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double f_low = f0 - div_low;
            const double f_high = f0 - div_high;
            // 7.1.6: track min after low-order update
            if (f_low < local_x_min_f_low) local_x_min_f_low = f_low;
            const int iv = k_int / Param::Nmu;
            const double f_tol =
                std::max(X_LOW_ABS_F_TOL,
                         X_LOW_REL_F_TOL * x_low_shell_max[iv]);
            const double f_scale = std::max(1.0, x_low_shell_max[iv]);
            if (!std::isfinite(f_low) || !std::isfinite(f_high)) {
                fluxes.cell_alpha_x[slot] = 0.0;
                local_x_low_order_failed += 1.0;
                local_x_any_limited = 1;
                continue;
            }
            const double f_floor = -f_tol;
            if (f0 < f_floor) {
                local_x_low_input_failed += 1.0;
                local_x_low_input_neg_mass +=
                    (-f0) * bkg_n.vgrid.moment_weight[iv] * sg.dx;
                local_x_low_input_rel_neg =
                    std::max(local_x_low_input_rel_neg, (-f0) / f_scale);
                const double x_cell =
                    static_cast<double>(sg.ix_start + ix) * sg.dx;
                const bool in_core =
                    (x_cell >= CORE_DIAG_BOUNDARY_WIDTH) &&
                    (x_cell <= Param::Lx - CORE_DIAG_BOUNDARY_WIDTH);
                if (in_core) {
                    local_x_low_input_core_failed += 1.0;
                }
                if (f_low < 0.0) {
                    local_x_low_output_rel_neg =
                        std::max(local_x_low_output_rel_neg,
                                 (-f_low) / f_scale);
                }
                fluxes.cell_alpha_x[slot] = 0.0;
                local_x_any_limited = 1;
                continue;
            }
            if (f_low < f_floor) {
                // Low-order update itself is negative beyond the shell-scaled
                // tolerance.  The global classification below distinguishes
                // true CFL from donor/periodic-face errors.
                local_x_low_order_failed += 1.0;
                local_x_low_output_rel_neg =
                    std::max(local_x_low_output_rel_neg, (-f_low) / f_scale);
                fluxes.cell_alpha_x[slot] = 0.0;
                local_x_any_limited = 1;
                continue;
            }
            // Fast path: high-order update is safe (common case, ~75-85%
            // of cells even with active limiter).  Branch hint for GCC/Clang.
            if (__builtin_expect(f_high >= f_floor, 1)) continue;
            // 7.1.6: track negative mass prevented by the limiter
            // f_high < f_floor means the high-order flux would create
            // negative mass. The negative mass prevented is (-f_high) * weight.
            if (f_high < 0.0) {
                const int iv2 = k_int / Param::Nmu;
                local_x_neg_mass_prevented +=
                    (-f_high) * bkg_n.vgrid.moment_weight[iv2] * sg.dx;
            }
            const double antidiff_out = div_high - div_low;
            if (antidiff_out <= 0.0) {
                fluxes.cell_alpha_x[slot] = 0.0;
            } else {
                fluxes.cell_alpha_x[slot] =
                    std::max(0.0, std::min(1.0,
                             (f_low - f_floor) / antidiff_out));
            }
            local_x_any_limited = 1;
        }
    }

    int global_x_any_limited = local_x_any_limited;
    MPI_Allreduce(MPI_IN_PLACE, &global_x_any_limited, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    std::vector<double> x_alpha_left_ghost;
    if (global_x_any_limited != 0 && nxl > 0) {
        x_alpha_left_ghost.assign(Param::Nvmu, 1.0);
        if (mpi_size <= 1) {
            const size_t last_base =
                static_cast<size_t>(nxl - 1) * Param::Nvmu;
            std::copy(fluxes.cell_alpha_x.begin() + last_base,
                      fluxes.cell_alpha_x.begin() + last_base + Param::Nvmu,
                      x_alpha_left_ghost.begin());
        } else {
            std::vector<double> x_alpha_send_right(Param::Nvmu, 1.0);
            const size_t last_base =
                static_cast<size_t>(nxl - 1) * Param::Nvmu;
            std::copy(fluxes.cell_alpha_x.begin() + last_base,
                      fluxes.cell_alpha_x.begin() + last_base + Param::Nvmu,
                      x_alpha_send_right.begin());
            const int left_rank = (mpi_rank + mpi_size - 1) % mpi_size;
            const int right_rank = (mpi_rank + 1) % mpi_size;
            MPI_Sendrecv(x_alpha_send_right.data(),
                         static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                         right_rank, 510,
                         x_alpha_left_ghost.data(),
                         static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                         left_rank, 510, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
        }
    }

    long long local_x_active_boundary = 0;
    long long local_x_total_boundary = 0;
    long long local_x_active_core = 0;
    long long local_x_total_core = 0;
    double local_x_min_alpha_boundary = 1.0;
    double local_x_min_alpha_core = 1.0;
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_x_active_boundary,local_x_total_boundary, \
                    local_x_active_core,local_x_total_core) \
        reduction(min:local_x_min_alpha_boundary,local_x_min_alpha_core)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
            const size_t k = static_cast<size_t>(k_int);
            double alpha_face = 1.0;
            if (global_x_any_limited != 0) {
                const double* alpha_left =
                    (iface == 0)
                    ? x_alpha_left_ghost.data()
                    : &fluxes.cell_alpha_x[
                        static_cast<size_t>(iface - 1) * Param::Nvmu];
                const double* alpha_right =
                    (iface == nxl)
                    ? fluxes.cell_alpha_x.data()
                    : &fluxes.cell_alpha_x[
                        static_cast<size_t>(iface) * Param::Nvmu];
                alpha_face =
                    std::min(alpha_left[k], alpha_right[k]);
                // 7.1.2: safety factor to prevent limiter overshoot
                alpha_face = std::max(0.0, std::min(1.0,
                    alpha_face * X_FCT_SAFETY));
            }
            const size_t face_slot =
                static_cast<size_t>(iface) * Param::Nvmu + k;
            fluxes.x_final[face_slot] =
                fluxes.x_low[face_slot] +
                alpha_face * (fluxes.x_high[face_slot] - fluxes.x_low[face_slot]);

            const double x_face =
                static_cast<double>(sg.ix_start + iface) * sg.dx;
            const bool boundary =
                (x_face < CORE_DIAG_BOUNDARY_WIDTH) ||
                (x_face > Param::Lx - CORE_DIAG_BOUNDARY_WIDTH);
            if (boundary) {
                ++local_x_total_boundary;
                if (alpha_face < 0.999999) ++local_x_active_boundary;
                local_x_min_alpha_boundary =
                    std::min(local_x_min_alpha_boundary, alpha_face);
            } else {
                ++local_x_total_core;
                if (alpha_face < 0.999999) ++local_x_active_core;
                local_x_min_alpha_core =
                    std::min(local_x_min_alpha_core, alpha_face);
            }
        }
    }
    close_periodic_face_blocks(fluxes.x_final, nxl,
                               static_cast<int>(Param::Nvmu),
                               mpi_rank, mpi_size, 505);

    /*
     * Yee-face force flux. The same x-face midpoint distribution is used for:
     *   1. x transport,
     *   2. background current passed to Ampere,
     *   3. velocity-space force work.
     * Each cell receives half of the force residual from its left face and
     * half from its right face. This is the actual unified residual path.
     */
    std::vector<double> u_force_face(
        static_cast<size_t>(nxl + 1) * (Param::Nv + 1) * Param::Nmu, 0.0);
    std::vector<double> mu_force_face(
        static_cast<size_t>(nxl + 1) * Param::Nv * (Param::Nmu + 1), 0.0);
    double mu_u_eff[Param::Nv];
    for (int iv = 0; iv < Param::Nv; ++iv) {
        mu_u_eff[iv] = shell_consistent_mu_u_eff(bkg_n.vgrid, iv);
    }

    #pragma omp parallel for schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        const double ex_face =
            (static_cast<size_t>(iface) < fields_mid.Ex_face.size())
            ? fields_mid.Ex_face[static_cast<size_t>(iface)] : 0.0;
        const double accel_u =
            bkg_n.charge * ex_face / (bkg_n.mass * Const::c);
        const int ix_left = ng + iface - 1;
        const int ix_right = ng + iface;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u_eff = mu_u_eff[iv];
            for (int face = 1; face < Param::Nmu; ++face) {
                const double mu_dot =
                    accel_u * bkg_n.vgrid.mu_face_factor[face] / u_eff;
                const size_t left_lo =
                    static_cast<size_t>(ix_left) * Param::Nvmu
                  + static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(face - 1);
                const size_t left_hi =
                    static_cast<size_t>(ix_left) * Param::Nvmu
                  + static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(face);
                const size_t right_lo =
                    static_cast<size_t>(ix_right) * Param::Nvmu
                  + static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(face - 1);
                const size_t right_hi =
                    static_cast<size_t>(ix_right) * Param::Nvmu
                  + static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(face);
                const double f_lo =
                    0.5 * (f_mid.f[left_lo] + f_mid.f[right_lo]);
                const double f_hi =
                    0.5 * (f_mid.f[left_hi] + f_mid.f[right_hi]);
                const double scale = bkg_n.vgrid.mu_flux_scale[iv] * mu_dot;
                mu_force_face[mu_xface_index(iface, iv, face)] =
                    scale * 0.5 * (f_lo + f_hi);
            }
        }
    }
    close_periodic_face_blocks(mu_force_face, nxl,
                               Param::Nv * (Param::Nmu + 1),
                               mpi_rank, mpi_size, 507);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            mu_force_face[mu_xface_index(iface, iv, 0)] = 0.0;
            mu_force_face[mu_xface_index(iface, iv, Param::Nmu)] = 0.0;
        }
    }

    // 7.1.1: save pre-limiter (high-order) mu fluxes.  u high-order
    // fluxes are now built cell-locally below.
    {
        const size_t mu_size = static_cast<size_t>(nxl + 1) *
                               Param::Nv * (Param::Nmu + 1);
        std::copy(mu_force_face.begin(), mu_force_face.begin() + mu_size,
                  fluxes.mu_high.begin());
    }
    /*
     * 7.1.3 finite-volume u-boundary closure.
     * The u-force update is closed at u=0 and u=u_max in this model.  These
     * faces are physical velocity-space boundaries, not ghost cells, so both
     * the high-order candidate and the working flux are explicitly zeroed
     * before any positivity limiting or final update can see them.
     */
    #pragma omp parallel for collapse(2) schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            fluxes.u_high[u_xface_index(iface, 0, imu)] = 0.0;
            fluxes.u_high[u_xface_index(iface, Param::Nv, imu)] = 0.0;
            u_force_face[u_xface_index(iface, 0, imu)] = 0.0;
            u_force_face[u_xface_index(iface, Param::Nv, imu)] = 0.0;
        }
    }

    /*
     * 7.1.3 / reconstruction-plan step 5:
     * construct a truly cell-local donor-cell u_low operator.  The positivity
     * base is f_after_x in each fixed (ix, imu) velocity column; no x-face
     * left/right average is allowed to decide the donor mass.  The x-face
     * arrays are filled only as compatibility diagnostics after the cell-local
     * fluxes have been built.
     */
    #pragma omp parallel for collapse(2) schedule(static)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            double f_base[Param::Nv];
            double f_floor[Param::Nv];
            double raw_low[Param::Nv + 1];
            double raw_high[Param::Nv + 1];
            double scale_face[Param::Nv + 1];
            double low_out[Param::Nv];
            double theta[Param::Nv];
            for (int iv = 0; iv < Param::Nv; ++iv) {
                const size_t k =
                    static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(imu);
                const size_t src =
                    static_cast<size_t>(ng + ix) * Param::Nvmu + k;
                const double dx_div =
                    dt_dx *
                    (fluxes.x_final[
                         static_cast<size_t>(ix + 1) * Param::Nvmu + k]
                   - fluxes.x_final[
                         static_cast<size_t>(ix) * Param::Nvmu + k]);
                f_base[iv] = bkg_n.f[src] - dx_div;
                const double local_scale =
                    std::max(1.0, std::fabs(f_base[iv]));
                f_floor[iv] = -eps_tol_base * local_scale;
                low_out[iv] = 0.0;
                theta[iv] = 1.0;
            }
            for (int face = 0; face <= Param::Nv; ++face) {
                raw_low[face] = 0.0;
                raw_high[face] = 0.0;
                scale_face[face] = 0.0;
            }
            const double ex_left =
                (static_cast<size_t>(ix) < fields_mid.Ex_face.size())
                ? fields_mid.Ex_face[static_cast<size_t>(ix)] : 0.0;
            const double ex_right =
                (static_cast<size_t>(ix + 1) < fields_mid.Ex_face.size())
                ? fields_mid.Ex_face[static_cast<size_t>(ix + 1)] : ex_left;
            const double accel_u =
                bkg_n.charge * 0.5 * (ex_left + ex_right) /
                (bkg_n.mass * Const::c);
            const double u_dot_base =
                accel_u * bkg_n.vgrid.mu_cells[imu];
            for (int face = 1; face < Param::Nv; ++face) {
                const double scale =
                    2.0 * Const::pi * bkg_n.vgrid.dmu
                  * bkg_n.vgrid.v2_faces[face] * u_dot_base;
                const int donor_iv = (scale >= 0.0) ? face - 1 : face;
                scale_face[face] = scale;
                raw_low[face] = scale * f_base[donor_iv];
                raw_high[face] =
                    scale * 0.5 * (f_base[face - 1] + f_base[face]);
                const double dt_inv =
                    dt * bkg_n.vgrid.inv_moment_weight[donor_iv];
                const double donor_decrement =
                    (scale >= 0.0)
                    ? dt_inv * std::max(0.0, raw_low[face])
                    : dt_inv * std::max(0.0, -raw_low[face]);
                low_out[donor_iv] += donor_decrement;
            }
            for (int iv = 0; iv < Param::Nv; ++iv) {
                if (low_out[iv] > 0.0) {
                    const double budget =
                        std::max(0.0, f_base[iv] - f_floor[iv]);
                    theta[iv] =
                        std::max(0.0, std::min(1.0, budget / low_out[iv]));
                }
            }
            for (int face = 1; face < Param::Nv; ++face) {
                const int donor_iv =
                    (scale_face[face] >= 0.0) ? face - 1 : face;
                raw_low[face] *= theta[donor_iv];
            }
            for (int face = 0; face <= Param::Nv; ++face) {
                const size_t idx = u_cell_index(ix, face, imu);
                fluxes.u_low_cell[idx] = raw_low[face];
                fluxes.u_high_cell[idx] = raw_high[face];
                fluxes.u_final_cell[idx] = raw_high[face];
            }
        }
    }
    #pragma omp parallel for collapse(2) schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int face = 0; face <= Param::Nv; ++face) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double low_val = 0.0;
                double high_val = 0.0;
                if (nxl > 0) {
                    if (iface == 0) {
                        low_val = fluxes.u_low_cell[u_cell_index(0, face, imu)];
                        high_val = fluxes.u_high_cell[u_cell_index(0, face, imu)];
                    } else if (iface == nxl) {
                        low_val =
                            fluxes.u_low_cell[u_cell_index(nxl - 1, face, imu)];
                        high_val =
                            fluxes.u_high_cell[u_cell_index(nxl - 1, face, imu)];
                    } else {
                        low_val = 0.5 *
                            (fluxes.u_low_cell[
                                 u_cell_index(iface - 1, face, imu)]
                           + fluxes.u_low_cell[
                                 u_cell_index(iface, face, imu)]);
                        high_val = 0.5 *
                            (fluxes.u_high_cell[
                                 u_cell_index(iface - 1, face, imu)]
                           + fluxes.u_high_cell[
                                 u_cell_index(iface, face, imu)]);
                    }
                }
                fluxes.u_low[u_xface_index(iface, face, imu)] = low_val;
                fluxes.u_high[u_xface_index(iface, face, imu)] = high_val;
                u_force_face[u_xface_index(iface, face, imu)] = high_val;
            }
        }
    }

    /*
     * Positivity-constrained u-force FCT limiter.
     *
     * F_u_final = F_u_low + alpha_u_face * (F_u_high - F_u_low).
     * This is true flux-corrected transport: first form the low-order state,
     * then limit only antidiffusive fluxes that would drain a cell beyond its
     * remaining positive budget.  The accepted update and energy-defect
     * accounting use only F_u_final.
     */
    bkg_new = bkg_n;
    const int low_u_limit_count =
        std::min(LOW_U_MU_LIMIT_COUNT, Param::Nv);
    // Precompute dt*inv_shell/2 and other per-iv constants once
    double half_dt_inv_shell[Param::Nv];
    double ke_per_mass_arr[Param::Nv];
    for (int iv = 0; iv < Param::Nv; ++iv) {
        half_dt_inv_shell[iv] = 0.5 * dt * bkg_n.vgrid.inv_moment_weight[iv];
        ke_per_mass_arr[iv] =
            (bkg_n.vgrid.gamma_cells[iv] - 1.0)
          * bkg_n.mass * Const::c * Const::c;
    }
    fluxes.mu_low_u_u_eff0 = mu_u_eff[0];
    fluxes.mu_low_u_moment_weight0 = bkg_n.vgrid.moment_weight[0];
    fluxes.mu_low_u_mu_flux_scale0 = bkg_n.vgrid.mu_flux_scale[0];
    fluxes.mu_low_u_half_dt_inv_shell0 = half_dt_inv_shell[0];
    double local_mu_low_u_dimless_scale0 = 0.0;
    double local_mu_low_u_endpoint_flux_max = 0.0;
    #pragma omp parallel for schedule(static) \
        reduction(max:local_mu_low_u_dimless_scale0, \
                     local_mu_low_u_endpoint_flux_max)
    for (int iface = 0; iface <= nxl; ++iface) {
        const double ex_face =
            (static_cast<size_t>(iface) < fields_mid.Ex_face.size())
            ? fields_mid.Ex_face[static_cast<size_t>(iface)] : 0.0;
        const double accel_u =
            bkg_n.charge * ex_face / (bkg_n.mass * Const::c);
        for (int face = 1; face < Param::Nmu; ++face) {
            const double scale =
                std::fabs(half_dt_inv_shell[0]
                        * bkg_n.vgrid.mu_flux_scale[0]
                        * accel_u
                        * bkg_n.vgrid.mu_face_factor[face]
                        / mu_u_eff[0]);
            local_mu_low_u_dimless_scale0 =
                std::max(local_mu_low_u_dimless_scale0, scale);
        }
        local_mu_low_u_endpoint_flux_max =
            std::max(local_mu_low_u_endpoint_flux_max,
                     std::fabs(mu_force_face[
                         mu_xface_index(iface, 0, 0)]));
        local_mu_low_u_endpoint_flux_max =
            std::max(local_mu_low_u_endpoint_flux_max,
                     std::fabs(mu_force_face[
                         mu_xface_index(iface, 0, Param::Nmu)]));
    }
    MPI_Allreduce(&local_mu_low_u_dimless_scale0,
                  &fluxes.mu_low_u_dimless_scale0, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_mu_low_u_endpoint_flux_max,
                  &fluxes.mu_low_u_endpoint_flux_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    // ---- Pass 1: compute per-cell alpha and detect if limiting is needed ----
    std::vector<double>& alpha_cell = fluxes.cell_alpha_u;
    std::fill(alpha_cell.begin(), alpha_cell.end(), 1.0);
    int local_any_limited = 0;
    long long local_ec_negative_cells = 0;
    long long local_u_low_order_failed = 0;
    double local_u_min_f_low = std::numeric_limits<double>::infinity();
    const int u_audit_threads = std::max(1, omp_get_max_threads());
    std::vector<UFluxAuditCandidate> thread_u_audit(
        static_cast<size_t>(u_audit_threads),
        empty_u_flux_audit_candidate());
    std::vector<LowOrderFailureAudit> thread_u_low_failure(
        static_cast<size_t>(u_audit_threads),
        empty_low_order_failure_audit());
    auto u_low_scale_at = [&](int ix, int face, int imu) -> double {
        if (face <= 0 || face >= Param::Nv) return 0.0;
        const double ex_left =
            (static_cast<size_t>(ix) < fields_mid.Ex_face.size())
            ? fields_mid.Ex_face[static_cast<size_t>(ix)] : 0.0;
        const double ex_right =
            (static_cast<size_t>(ix + 1) < fields_mid.Ex_face.size())
            ? fields_mid.Ex_face[static_cast<size_t>(ix + 1)] : ex_left;
        const double accel_u =
            bkg_n.charge * 0.5 * (ex_left + ex_right) /
            (bkg_n.mass * Const::c);
        const double u_dot_base = accel_u * bkg_n.vgrid.mu_cells[imu];
        return 2.0 * Const::pi * bkg_n.vgrid.dmu
             * bkg_n.vgrid.v2_faces[face] * u_dot_base;
    };
    auto u_low_donor_iv_at = [&](double scale, int face) -> int {
        if (face <= 0 || face >= Param::Nv) return -1;
        return (scale >= 0.0) ? face - 1 : face;
    };
    auto u_low_donor_f_at =
        [&](int ix, int donor_iv, int imu) -> double {
            if (donor_iv < 0 || donor_iv >= Param::Nv) return 0.0;
            const size_t k =
                static_cast<size_t>(donor_iv) * Param::Nmu
              + static_cast<size_t>(imu);
            const size_t src =
                static_cast<size_t>(ng + ix) * Param::Nvmu + k;
            const double dx_div =
                dt_dx *
                (fluxes.x_final[
                     static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_final[
                     static_cast<size_t>(ix) * Param::Nvmu + k]);
            return bkg_n.f[src] - dx_div;
        };
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(max:local_any_limited) \
        reduction(+:local_ec_negative_cells,local_u_low_order_failed) \
        reduction(min:local_u_min_f_low)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
            const int iv = k_int / Param::Nmu;
            const int imu = k_int - iv * Param::Nmu;
            const size_t k = static_cast<size_t>(k_int);
            const size_t ix_k =
                static_cast<size_t>(ix) * Param::Nvmu + k;
            const size_t src =
                static_cast<size_t>(ng + ix) * Param::Nvmu + k;
            const double f0 = bkg_n.f[src];
            const double hdt_is = half_dt_inv_shell[iv];
            const double dx_div =
                dt_dx *
                (fluxes.x_final[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_final[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double f_after_x = f0 - dx_div;
            const double dt_inv_shell = 2.0 * hdt_is;
            const double u_div_low =
                fluxes.u_low_cell[u_cell_index(ix, iv + 1, imu)]
              - fluxes.u_low_cell[u_cell_index(ix, iv, imu)];
            const double du_div_low = dt_inv_shell * u_div_low;
            const double u_div_high =
                fluxes.u_high_cell[u_cell_index(ix, iv + 1, imu)]
              - fluxes.u_high_cell[u_cell_index(ix, iv, imu)];
            const double du_div_high = dt_inv_shell * u_div_high;
            const double local_scale = std::max(1.0, std::fabs(f_after_x));
            const double f_floor = -eps_tol_base * local_scale;
            const double f_low = f_after_x - du_div_low;
            const double f_high = f_after_x - du_div_high;
            local_u_min_f_low = std::min(local_u_min_f_low, f_low);
            double audit_alpha = 1.0;
            double audit_severity = 0.0;
            if (!std::isfinite(f_low) || !std::isfinite(f_high)) {
                alpha_cell[ix_k] = 0.0;
                ++local_u_low_order_failed;
                local_any_limited = 1;
                audit_alpha = 0.0;
                audit_severity = std::numeric_limits<double>::infinity();
                const int tid = omp_get_thread_num();
                UFluxAuditCandidate& audit =
                    thread_u_audit[static_cast<size_t>(tid)];
                if (!audit.valid || audit_severity > audit.severity) {
                    audit.valid = 1;
                    audit.severity = audit_severity;
                    audit.f0 = f0;
                    audit.f_low = f_low;
                    audit.f_high = f_high;
                    audit.alpha = audit_alpha;
                    audit.du_div_low = du_div_low;
                    audit.du_div_high = du_div_high;
                    audit.ix = sg.ix_start + ix;
                    audit.iv = iv;
                    audit.imu = imu;
                }
                const int low_tid = omp_get_thread_num();
                LowOrderFailureAudit& low_audit =
                    thread_u_low_failure[static_cast<size_t>(low_tid)];
                if (!low_audit.valid ||
                    audit_severity > low_audit.severity) {
                    const int lower_face = iv;
                    const int upper_face = iv + 1;
                    const double left_lower_scale =
                        u_low_scale_at(ix, lower_face, imu);
                    const double left_upper_scale =
                        u_low_scale_at(ix, upper_face, imu);
                    const double right_lower_scale =
                        left_lower_scale;
                    const double right_upper_scale =
                        left_upper_scale;
                    const int left_lower_donor =
                        u_low_donor_iv_at(left_lower_scale, lower_face);
                    const int left_upper_donor =
                        u_low_donor_iv_at(left_upper_scale, upper_face);
                    const int right_lower_donor =
                        u_low_donor_iv_at(right_lower_scale, lower_face);
                    const int right_upper_donor =
                        u_low_donor_iv_at(right_upper_scale, upper_face);
                    const double x_cell =
                        (static_cast<double>(sg.ix_start + ix) + 0.5)
                      * sg.dx;
                    low_audit.valid = 1;
                    low_audit.rank = mpi_rank;
                    low_audit.ix = sg.ix_start + ix;
                    low_audit.iv = iv;
                    low_audit.imu = imu;
                    low_audit.region =
                        (x_cell < CORE_DIAG_BOUNDARY_WIDTH ||
                         x_cell > Param::Lx - CORE_DIAG_BOUNDARY_WIDTH)
                        ? 1 : 0;
                    low_audit.severity = audit_severity;
                    low_audit.f_input = f0;
                    low_audit.f_after_x = f0 - dx_div;
                    low_audit.dx_div = dx_div;
                    low_audit.dmu_div_used = 0.0;
                    low_audit.du_div_low = du_div_low;
                    low_audit.f_low = f_low;
                    low_audit.left_lower_flux =
                        fluxes.u_low_cell[u_cell_index(ix, lower_face, imu)];
                    low_audit.left_upper_flux =
                        fluxes.u_low_cell[u_cell_index(ix, upper_face, imu)];
                    low_audit.right_lower_flux =
                        low_audit.left_lower_flux;
                    low_audit.right_upper_flux =
                        low_audit.left_upper_flux;
                    low_audit.left_lower_scale = left_lower_scale;
                    low_audit.left_upper_scale = left_upper_scale;
                    low_audit.right_lower_scale = right_lower_scale;
                    low_audit.right_upper_scale = right_upper_scale;
                    low_audit.left_lower_donor_f =
                        u_low_donor_f_at(ix, left_lower_donor, imu);
                    low_audit.left_upper_donor_f =
                        u_low_donor_f_at(ix, left_upper_donor, imu);
                    low_audit.right_lower_donor_f =
                        low_audit.left_lower_donor_f;
                    low_audit.right_upper_donor_f =
                        low_audit.left_upper_donor_f;
                    low_audit.lower_characteristic = 0.5 *
                        (left_lower_scale + right_lower_scale);
                    low_audit.upper_characteristic = 0.5 *
                        (left_upper_scale + right_upper_scale);
                    low_audit.moment_weight = bkg_n.vgrid.moment_weight[iv];
                    low_audit.cell_budget =
                        std::isfinite(f_low) ? f_low - f_floor : 0.0;
                    low_audit.left_lower_donor_index = left_lower_donor;
                    low_audit.left_upper_donor_index = left_upper_donor;
                    low_audit.right_lower_donor_index = right_lower_donor;
                    low_audit.right_upper_donor_index = right_upper_donor;
                }
                continue;
            }
            if (f_low < f_floor) {
                alpha_cell[ix_k] = 0.0;
                ++local_u_low_order_failed;
                local_any_limited = 1;
                audit_alpha = 0.0;
                audit_severity = std::max(0.0, -f_low);
                const int tid = omp_get_thread_num();
                UFluxAuditCandidate& audit =
                    thread_u_audit[static_cast<size_t>(tid)];
                if (!audit.valid || audit_severity > audit.severity) {
                    audit.valid = 1;
                    audit.severity = audit_severity;
                    audit.f0 = f0;
                    audit.f_low = f_low;
                    audit.f_high = f_high;
                    audit.alpha = audit_alpha;
                    audit.du_div_low = du_div_low;
                    audit.du_div_high = du_div_high;
                    audit.ix = sg.ix_start + ix;
                    audit.iv = iv;
                    audit.imu = imu;
                }
                const int low_tid = omp_get_thread_num();
                LowOrderFailureAudit& low_audit =
                    thread_u_low_failure[static_cast<size_t>(low_tid)];
                if (!low_audit.valid ||
                    audit_severity > low_audit.severity) {
                    const int lower_face = iv;
                    const int upper_face = iv + 1;
                    const double left_lower_scale =
                        u_low_scale_at(ix, lower_face, imu);
                    const double left_upper_scale =
                        u_low_scale_at(ix, upper_face, imu);
                    const double right_lower_scale =
                        left_lower_scale;
                    const double right_upper_scale =
                        left_upper_scale;
                    const int left_lower_donor =
                        u_low_donor_iv_at(left_lower_scale, lower_face);
                    const int left_upper_donor =
                        u_low_donor_iv_at(left_upper_scale, upper_face);
                    const int right_lower_donor =
                        u_low_donor_iv_at(right_lower_scale, lower_face);
                    const int right_upper_donor =
                        u_low_donor_iv_at(right_upper_scale, upper_face);
                    const double x_cell =
                        (static_cast<double>(sg.ix_start + ix) + 0.5)
                      * sg.dx;
                    low_audit.valid = 1;
                    low_audit.rank = mpi_rank;
                    low_audit.ix = sg.ix_start + ix;
                    low_audit.iv = iv;
                    low_audit.imu = imu;
                    low_audit.region =
                        (x_cell < CORE_DIAG_BOUNDARY_WIDTH ||
                         x_cell > Param::Lx - CORE_DIAG_BOUNDARY_WIDTH)
                        ? 1 : 0;
                    low_audit.severity = audit_severity;
                    low_audit.f_input = f0;
                    low_audit.f_after_x = f0 - dx_div;
                    low_audit.dx_div = dx_div;
                    low_audit.dmu_div_used = 0.0;
                    low_audit.du_div_low = du_div_low;
                    low_audit.f_low = f_low;
                    low_audit.left_lower_flux =
                        fluxes.u_low_cell[u_cell_index(ix, lower_face, imu)];
                    low_audit.left_upper_flux =
                        fluxes.u_low_cell[u_cell_index(ix, upper_face, imu)];
                    low_audit.right_lower_flux =
                        low_audit.left_lower_flux;
                    low_audit.right_upper_flux =
                        low_audit.left_upper_flux;
                    low_audit.left_lower_scale = left_lower_scale;
                    low_audit.left_upper_scale = left_upper_scale;
                    low_audit.right_lower_scale = right_lower_scale;
                    low_audit.right_upper_scale = right_upper_scale;
                    low_audit.left_lower_donor_f =
                        u_low_donor_f_at(ix, left_lower_donor, imu);
                    low_audit.left_upper_donor_f =
                        u_low_donor_f_at(ix, left_upper_donor, imu);
                    low_audit.right_lower_donor_f =
                        low_audit.left_lower_donor_f;
                    low_audit.right_upper_donor_f =
                        low_audit.left_upper_donor_f;
                    low_audit.lower_characteristic = 0.5 *
                        (left_lower_scale + right_lower_scale);
                    low_audit.upper_characteristic = 0.5 *
                        (left_upper_scale + right_upper_scale);
                    low_audit.moment_weight = bkg_n.vgrid.moment_weight[iv];
                    low_audit.cell_budget = f_low - f_floor;
                    low_audit.left_lower_donor_index = left_lower_donor;
                    low_audit.left_upper_donor_index = left_upper_donor;
                    low_audit.right_lower_donor_index = right_lower_donor;
                    low_audit.right_upper_donor_index = right_upper_donor;
                }
                continue;
            }
            const double a_lo =
                fluxes.u_high_cell[u_cell_index(ix, iv, imu)]
              - fluxes.u_low_cell[u_cell_index(ix, iv, imu)];
            const double a_hi =
                fluxes.u_high_cell[u_cell_index(ix, iv + 1, imu)]
              - fluxes.u_low_cell[u_cell_index(ix, iv + 1, imu)];
            const double antidiff_out =
                dt_inv_shell * (std::max(0.0,  a_hi)
                              + std::max(0.0, -a_lo));
            if (f_high < f_floor) ++local_ec_negative_cells;
            if (antidiff_out > 0.0) {
                const double budget = std::max(0.0, f_low - f_floor);
                const double alpha =
                    std::max(0.0, std::min(1.0, budget / antidiff_out));
                alpha_cell[ix_k] = alpha;
                audit_alpha = alpha;
                if (alpha < 1.0 - 1.0e-14) local_any_limited = 1;
            }
            audit_severity =
                std::max(0.0, -std::min(f_low, f_high));
            if (audit_severity > 0.0) {
                const int tid = omp_get_thread_num();
                UFluxAuditCandidate& audit =
                    thread_u_audit[static_cast<size_t>(tid)];
                if (!audit.valid || audit_severity > audit.severity) {
                    audit.valid = 1;
                    audit.severity = audit_severity;
                    audit.f0 = f0;
                    audit.f_low = f_low;
                    audit.f_high = f_high;
                    audit.alpha = audit_alpha;
                    audit.du_div_low = du_div_low;
                    audit.du_div_high = du_div_high;
                    audit.ix = sg.ix_start + ix;
                    audit.iv = iv;
                    audit.imu = imu;
                }
            }
        }
    }
    int global_any_limited = local_any_limited;
    MPI_Allreduce(MPI_IN_PLACE, &global_any_limited, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    long long global_u_low_order_failed_for_audit =
        local_u_low_order_failed;
    MPI_Allreduce(MPI_IN_PLACE, &global_u_low_order_failed_for_audit, 1,
                  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    (void)local_ec_negative_cells;
    UFluxAuditCandidate local_u_audit = empty_u_flux_audit_candidate();
    for (const UFluxAuditCandidate& audit : thread_u_audit) {
        if (audit.valid &&
            (!local_u_audit.valid ||
             audit.severity > local_u_audit.severity)) {
            local_u_audit = audit;
        }
    }
    struct {
        double value;
        int rank;
    } local_u_audit_loc, global_u_audit_loc;
    local_u_audit_loc.value =
        local_u_audit.valid ? local_u_audit.severity : 0.0;
    local_u_audit_loc.rank = mpi_rank;
    global_u_audit_loc.value = 0.0;
    global_u_audit_loc.rank = -1;
    MPI_Allreduce(&local_u_audit_loc, &global_u_audit_loc, 1,
                  MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    double u_audit_values[7] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    int u_audit_indices[4] = {-1, -1, -1, -1};
    if (local_u_audit.valid && mpi_rank == global_u_audit_loc.rank) {
        u_audit_values[0] = local_u_audit.severity;
        u_audit_values[1] = local_u_audit.f0;
        u_audit_values[2] = local_u_audit.f_low;
        u_audit_values[3] = local_u_audit.f_high;
        u_audit_values[4] = local_u_audit.alpha;
        u_audit_values[5] = local_u_audit.du_div_low;
        u_audit_values[6] = local_u_audit.du_div_high;
        u_audit_indices[0] = mpi_rank;
        u_audit_indices[1] = local_u_audit.ix;
        u_audit_indices[2] = local_u_audit.iv;
        u_audit_indices[3] = local_u_audit.imu;
    }
    if (global_u_audit_loc.value > 0.0 &&
        global_u_audit_loc.rank >= 0) {
        MPI_Bcast(u_audit_values, 7, MPI_DOUBLE,
                  global_u_audit_loc.rank, MPI_COMM_WORLD);
        MPI_Bcast(u_audit_indices, 4, MPI_INT,
                  global_u_audit_loc.rank, MPI_COMM_WORLD);
        fluxes.u_flux_audit_valid = 1;
        fluxes.u_flux_audit_rank = u_audit_indices[0];
        fluxes.u_flux_audit_ix = u_audit_indices[1];
        fluxes.u_flux_audit_iv = u_audit_indices[2];
        fluxes.u_flux_audit_imu = u_audit_indices[3];
        fluxes.u_flux_audit_severity = u_audit_values[0];
        fluxes.u_flux_audit_f0 = u_audit_values[1];
        fluxes.u_flux_audit_f_low = u_audit_values[2];
        fluxes.u_flux_audit_f_high = u_audit_values[3];
        fluxes.u_flux_audit_alpha = u_audit_values[4];
        fluxes.u_flux_audit_du_div_low = u_audit_values[5];
        fluxes.u_flux_audit_du_div_high = u_audit_values[6];
    }
    LowOrderFailureAudit local_u_low_audit =
        empty_low_order_failure_audit();
    for (const LowOrderFailureAudit& audit : thread_u_low_failure) {
        if (audit.valid &&
            (!local_u_low_audit.valid ||
             audit.severity > local_u_low_audit.severity)) {
            local_u_low_audit = audit;
        }
    }
    struct {
        double value;
        int rank;
    } local_u_low_loc, global_u_low_loc;
    local_u_low_loc.value =
        local_u_low_audit.valid ? local_u_low_audit.severity : 0.0;
    local_u_low_loc.rank = mpi_rank;
    global_u_low_loc.value = 0.0;
    global_u_low_loc.rank = -1;
    MPI_Allreduce(&local_u_low_loc, &global_u_low_loc, 1,
                  MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    double u_low_values[27] = {0.0};
    int u_low_indices[10] = {-1, -1, -1, -1, 0, -1, -1, -1, -1, 0};
    if (local_u_low_audit.valid &&
        mpi_rank == global_u_low_loc.rank) {
        u_low_values[0] = local_u_low_audit.severity;
        u_low_values[1] = local_u_low_audit.f_input;
        u_low_values[2] = local_u_low_audit.f_after_x;
        u_low_values[3] = local_u_low_audit.dx_div;
        u_low_values[4] = local_u_low_audit.dmu_div_used;
        u_low_values[5] = local_u_low_audit.du_div_low;
        u_low_values[6] = local_u_low_audit.f_low;
        u_low_values[7] = local_u_low_audit.left_lower_flux;
        u_low_values[8] = local_u_low_audit.left_upper_flux;
        u_low_values[9] = local_u_low_audit.right_lower_flux;
        u_low_values[10] = local_u_low_audit.right_upper_flux;
        u_low_values[11] = local_u_low_audit.left_lower_scale;
        u_low_values[12] = local_u_low_audit.left_upper_scale;
        u_low_values[13] = local_u_low_audit.right_lower_scale;
        u_low_values[14] = local_u_low_audit.right_upper_scale;
        u_low_values[15] = local_u_low_audit.left_lower_donor_f;
        u_low_values[16] = local_u_low_audit.left_upper_donor_f;
        u_low_values[17] = local_u_low_audit.right_lower_donor_f;
        u_low_values[18] = local_u_low_audit.right_upper_donor_f;
        u_low_values[19] = local_u_low_audit.lower_characteristic;
        u_low_values[20] = local_u_low_audit.upper_characteristic;
        u_low_values[21] = local_u_low_audit.moment_weight;
        u_low_values[22] = local_u_low_audit.cell_budget;
        u_low_values[23] =
            static_cast<double>(global_u_low_order_failed_for_audit);
        u_low_indices[0] = local_u_low_audit.rank;
        u_low_indices[1] = local_u_low_audit.ix;
        u_low_indices[2] = local_u_low_audit.iv;
        u_low_indices[3] = local_u_low_audit.imu;
        u_low_indices[4] = local_u_low_audit.region;
        u_low_indices[5] = local_u_low_audit.left_lower_donor_index;
        u_low_indices[6] = local_u_low_audit.left_upper_donor_index;
        u_low_indices[7] = local_u_low_audit.right_lower_donor_index;
        u_low_indices[8] = local_u_low_audit.right_upper_donor_index;
    }
    if (global_u_low_loc.value > 0.0 && global_u_low_loc.rank >= 0) {
        MPI_Bcast(u_low_values, 27, MPI_DOUBLE,
                  global_u_low_loc.rank, MPI_COMM_WORLD);
        MPI_Bcast(u_low_indices, 10, MPI_INT,
                  global_u_low_loc.rank, MPI_COMM_WORLD);
        LowOrderFailureAudit& audit = fluxes.u_low_failure_audit;
        audit.valid = 1;
        audit.rank = u_low_indices[0];
        audit.ix = u_low_indices[1];
        audit.iv = u_low_indices[2];
        audit.imu = u_low_indices[3];
        audit.region = u_low_indices[4];
        audit.severity = u_low_values[0];
        audit.f_input = u_low_values[1];
        audit.f_after_x = u_low_values[2];
        audit.dx_div = u_low_values[3];
        audit.dmu_div_used = u_low_values[4];
        audit.du_div_low = u_low_values[5];
        audit.f_low = u_low_values[6];
        audit.left_lower_flux = u_low_values[7];
        audit.left_upper_flux = u_low_values[8];
        audit.right_lower_flux = u_low_values[9];
        audit.right_upper_flux = u_low_values[10];
        audit.left_lower_scale = u_low_values[11];
        audit.left_upper_scale = u_low_values[12];
        audit.right_lower_scale = u_low_values[13];
        audit.right_upper_scale = u_low_values[14];
        audit.left_lower_donor_f = u_low_values[15];
        audit.left_upper_donor_f = u_low_values[16];
        audit.right_lower_donor_f = u_low_values[17];
        audit.right_upper_donor_f = u_low_values[18];
        audit.lower_characteristic = u_low_values[19];
        audit.upper_characteristic = u_low_values[20];
        audit.moment_weight = u_low_values[21];
        audit.cell_budget = u_low_values[22];
        audit.low_order_failed_count = u_low_values[23];
        audit.left_lower_donor_index = u_low_indices[5];
        audit.left_upper_donor_index = u_low_indices[6];
        audit.right_lower_donor_index = u_low_indices[7];
        audit.right_upper_donor_index = u_low_indices[8];
    }

    const bool need_face_limiting = (global_any_limited != 0);
    long long local_u_face_active = 0;
    long long local_u_face_total = 0;
    double local_u_alpha_min = 1.0;

    if (need_face_limiting) {
    // ---- Pass 2: cell-local face alpha and construct F_u_final_cell ----
    #pragma omp parallel for schedule(static) \
        reduction(+:local_u_face_active,local_u_face_total) \
        reduction(min:local_u_alpha_min)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int face = 0; face <= Param::Nv; ++face) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t face_idx = u_cell_index(ix, face, imu);
                const double antidiff =
                    fluxes.u_high_cell[face_idx] -
                    fluxes.u_low_cell[face_idx];
                double alpha_face = 1.0;
                if (face > 0 && face < Param::Nv && antidiff > 0.0) {
                    const size_t k_lo =
                        static_cast<size_t>(ix) * Param::Nvmu +
                        static_cast<size_t>(face - 1) * Param::Nmu
                      + static_cast<size_t>(imu);
                    alpha_face = alpha_cell[k_lo];
                } else if (face > 0 && face < Param::Nv && antidiff < 0.0) {
                    const size_t k_hi =
                        static_cast<size_t>(ix) * Param::Nvmu +
                        static_cast<size_t>(face) * Param::Nmu
                      + static_cast<size_t>(imu);
                    alpha_face = alpha_cell[k_hi];
                }
                alpha_face = std::max(0.0, std::min(1.0, alpha_face));
                if (face <= low_u_limit_count &&
                    alpha_face < 1.0 - 1.0e-12) {
                    alpha_face = smooth_low_u_mu_alpha(alpha_face);
                }
                fluxes.u_final_cell[face_idx] =
                    fluxes.u_low_cell[face_idx] +
                    alpha_face *
                    (fluxes.u_high_cell[face_idx] -
                     fluxes.u_low_cell[face_idx]);
                if (alpha_face < 0.999999) ++local_u_face_active;
                ++local_u_face_total;
                local_u_alpha_min =
                    std::min(local_u_alpha_min, alpha_face);
            }
        }
    }
    } // need_face_limiting

    // Map cell-local final u flux back to x-face storage for diagnostics only.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int face = 0; face <= Param::Nv; ++face) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                double final_val = 0.0;
                if (nxl > 0) {
                    if (iface == 0) {
                        final_val =
                            fluxes.u_final_cell[u_cell_index(0, face, imu)];
                    } else if (iface == nxl) {
                        final_val =
                            fluxes.u_final_cell[
                                u_cell_index(nxl - 1, face, imu)];
                    } else {
                        final_val = 0.5 *
                            (fluxes.u_final_cell[
                                 u_cell_index(iface - 1, face, imu)]
                           + fluxes.u_final_cell[
                                 u_cell_index(iface, face, imu)]);
                    }
                }
                fluxes.u_final[u_xface_index(iface, face, imu)] = final_val;
                u_force_face[u_xface_index(iface, face, imu)] = final_val;
            }
        }
    }

    std::vector<double> mu_force_low_u_ec(
        static_cast<size_t>(nxl + 1) * low_u_limit_count *
        static_cast<size_t>(Param::Nmu + 1), 0.0);
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int iv = 0; iv < low_u_limit_count; ++iv) {
            for (int face = 0; face <= Param::Nmu; ++face) {
                const size_t low_idx =
                    (static_cast<size_t>(iface) * low_u_limit_count
                   + static_cast<size_t>(iv)) *
                    static_cast<size_t>(Param::Nmu + 1)
                  + static_cast<size_t>(face);
                mu_force_low_u_ec[low_idx] =
                    mu_force_face[mu_xface_index(iface, iv, face)];
            }
        }
    }
    long long local_mu_low_u_face_active = 0;
    long long local_mu_low_u_face_total = 0;
    long long local_mu_low_u_face_active_boundary = 0;
    long long local_mu_low_u_face_total_boundary = 0;
    long long local_mu_low_u_face_active_core = 0;
    long long local_mu_low_u_face_total_core = 0;
    long long local_low_u_subcycle_active = 0;
    long long local_low_u_subcycle_total = 0;
    long long local_low_u_subcycle_sum = 0;
    int local_low_u_max_subcycles = 1;
    double local_mu_low_u_alpha_min = 1.0;
    double local_mu_low_u_alpha_min_boundary = 1.0;
    double local_mu_low_u_alpha_min_core = 1.0;
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_mu_low_u_face_active, \
                    local_mu_low_u_face_total, \
                    local_mu_low_u_face_active_boundary, \
                    local_mu_low_u_face_total_boundary, \
                    local_mu_low_u_face_active_core, \
                    local_mu_low_u_face_total_core, \
                    local_low_u_subcycle_active, \
                    local_low_u_subcycle_total, \
                    local_low_u_subcycle_sum) \
        reduction(min:local_mu_low_u_alpha_min, \
                      local_mu_low_u_alpha_min_boundary, \
                      local_mu_low_u_alpha_min_core) \
        reduction(max:local_low_u_max_subcycles)
    for (int iface = 0; iface < nxl; ++iface) {
        for (int iv = 0; iv < low_u_limit_count; ++iv) {
            const int ix_left = ng + iface - 1;
            const int ix_right = ng + iface;
            const double x_face =
                static_cast<double>(sg.ix_start + iface) * sg.dx;
            const bool boundary =
                (x_face < CORE_DIAG_BOUNDARY_WIDTH) ||
                (x_face > Param::Lx - CORE_DIAG_BOUNDARY_WIDTH);
            for (int face = 1; face < Param::Nmu; ++face) {
                const size_t face_idx = mu_xface_index(iface, iv, face);
                double& ff = mu_force_face[face_idx];
                double alpha_face = 1.0;
                if (ff != 0.0 && std::isfinite(ff)) {
                    const int donor_imu = (ff > 0.0) ? face - 1 : face;
                    const size_t left_donor =
                        static_cast<size_t>(ix_left) * Param::Nvmu
                      + static_cast<size_t>(iv) * Param::Nmu
                      + static_cast<size_t>(donor_imu);
                    const size_t right_donor =
                        static_cast<size_t>(ix_right) * Param::Nvmu
                      + static_cast<size_t>(iv) * Param::Nmu
                      + static_cast<size_t>(donor_imu);
                    const double donor_face_f =
                        0.5 * (f_mid.f[left_donor] + f_mid.f[right_donor]);
                    const double available = std::max(0.0, donor_face_f);
                    const double full_out =
                        half_dt_inv_shell[iv] * std::fabs(ff);
                    if (!std::isfinite(donor_face_f) ||
                        !std::isfinite(full_out)) {
                        alpha_face = 0.0;
                    } else if (full_out > 0.0) {
                        const int subcycles =
                            adaptive_low_u_subcycles(full_out, available);
                        ++local_low_u_subcycle_total;
                        local_low_u_subcycle_sum += subcycles;
                        local_low_u_max_subcycles =
                            std::max(local_low_u_max_subcycles, subcycles);
                        const bool adaptive_region = subcycles > 1;
                        if (adaptive_region) {
                            ++local_low_u_subcycle_active;
                        }
                        const double full_alpha =
                            std::min(1.0, available / full_out);
                        alpha_face = full_alpha;
                        if (adaptive_region &&
                            (face == 1 || face == Param::Nmu - 1)) {
                            alpha_face =
                                std::min(alpha_face,
                                         LOW_U_MU_ENDPOINT_ALPHA_CAP);
                        }
                        if (adaptive_region || alpha_face < 1.0 - 1.0e-12) {
                            alpha_face = smooth_low_u_mu_alpha(alpha_face);
                        }
                    }
                    alpha_face = std::max(0.0, std::min(1.0, alpha_face));
                    ff *= alpha_face;
                } else if (!std::isfinite(ff)) {
                    alpha_face = 0.0;
                    ff = 0.0;
                }
                if (alpha_face < 0.999999) {
                    ++local_mu_low_u_face_active;
                    if (boundary) {
                        ++local_mu_low_u_face_active_boundary;
                    } else {
                        ++local_mu_low_u_face_active_core;
                    }
                }
                ++local_mu_low_u_face_total;
                local_mu_low_u_alpha_min =
                    std::min(local_mu_low_u_alpha_min, alpha_face);
                if (boundary) {
                    ++local_mu_low_u_face_total_boundary;
                    local_mu_low_u_alpha_min_boundary =
                        std::min(local_mu_low_u_alpha_min_boundary,
                                 alpha_face);
                } else {
                    ++local_mu_low_u_face_total_core;
                    local_mu_low_u_alpha_min_core =
                        std::min(local_mu_low_u_alpha_min_core,
                                 alpha_face);
                }
            }
        }
    }
    close_periodic_face_blocks(mu_force_face, nxl,
                               Param::Nv * (Param::Nmu + 1),
                               mpi_rank, mpi_size, 507);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        for (int iv = 0; iv < Param::Nv; ++iv) {
            mu_force_face[mu_xface_index(iface, iv, 0)] = 0.0;
            mu_force_face[mu_xface_index(iface, iv, Param::Nmu)] = 0.0;
        }
    }
    // 7.1.1: save post-limiter (final) mu-flux; mu_low not yet implemented
    {
        const size_t mu_size = static_cast<size_t>(nxl + 1) *
                               Param::Nv * (Param::Nmu + 1);
        std::copy(mu_force_face.begin(), mu_force_face.begin() + mu_size,
                  fluxes.mu_final.begin());
        std::copy(fluxes.mu_final.begin(), fluxes.mu_final.begin() + mu_size,
                  fluxes.mu_low.begin());
    }
    {
        const int mu_audit_threads = std::max(1, omp_get_max_threads());
        std::vector<LowOrderFailureAudit> thread_mu_low_failure(
            static_cast<size_t>(mu_audit_threads),
            empty_low_order_failure_audit());
        long long local_mu_low_order_failed = 0;
        auto mu_dot_at = [&](int iface, int iv, int face) -> double {
            if (face <= 0 || face >= Param::Nmu) return 0.0;
            const double ex_face =
                (static_cast<size_t>(iface) < fields_mid.Ex_face.size())
                ? fields_mid.Ex_face[static_cast<size_t>(iface)] : 0.0;
            const double accel_u =
                bkg_n.charge * ex_face / (bkg_n.mass * Const::c);
            return accel_u * bkg_n.vgrid.mu_face_factor[face]
                 / mu_u_eff[iv];
        };
        auto mu_donor_imu_at = [&](double characteristic, int face) -> int {
            if (face <= 0 || face >= Param::Nmu) return -1;
            return (characteristic >= 0.0) ? face - 1 : face;
        };
        #pragma omp parallel for collapse(2) schedule(static) \
            reduction(+:local_mu_low_order_failed)
        for (int ix = 0; ix < nxl; ++ix) {
            for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu);
                 ++k_int) {
                const int iv = k_int / Param::Nmu;
                const int imu = k_int - iv * Param::Nmu;
                const size_t k = static_cast<size_t>(k_int);
                const size_t src =
                    static_cast<size_t>(ng + ix) * Param::Nvmu + k;
                const double hdt_is = half_dt_inv_shell[iv];
                const double dx_div =
                    dt_dx *
                    (fluxes.x_final[
                         static_cast<size_t>(ix + 1) * Param::Nvmu + k]
                   - fluxes.x_final[
                         static_cast<size_t>(ix) * Param::Nvmu + k]);
                const double u_div =
                    fluxes.u_final_cell[u_cell_index(ix, iv + 1, imu)]
                  - fluxes.u_final_cell[u_cell_index(ix, iv, imu)];
                const double du_div = (2.0 * hdt_is) * u_div;
                const double mu_left_lower =
                    fluxes.mu_low[mu_xface_index(ix, iv, imu)];
                const double mu_left_upper =
                    fluxes.mu_low[mu_xface_index(ix, iv, imu + 1)];
                const double mu_right_lower =
                    fluxes.mu_low[mu_xface_index(ix + 1, iv, imu)];
                const double mu_right_upper =
                    fluxes.mu_low[mu_xface_index(ix + 1, iv, imu + 1)];
                const double dmu_div_low =
                    hdt_is * ((mu_left_upper - mu_left_lower)
                            + (mu_right_upper - mu_right_lower));
                const double f_before_mu = bkg_n.f[src] - dx_div - du_div;
                const double f_mu_low = f_before_mu - dmu_div_low;
                const double local_scale =
                    std::max(1.0, std::fabs(f_before_mu));
                const double f_floor = -eps_tol_base * local_scale;
                bool failed = false;
                double severity = 0.0;
                if (!std::isfinite(f_mu_low)) {
                    failed = true;
                    severity = std::numeric_limits<double>::infinity();
                } else if (f_mu_low < f_floor) {
                    failed = true;
                    severity = -f_mu_low;
                }
                if (!failed) continue;
                ++local_mu_low_order_failed;
                const int tid = omp_get_thread_num();
                LowOrderFailureAudit& audit =
                    thread_mu_low_failure[static_cast<size_t>(tid)];
                if (!audit.valid || severity > audit.severity) {
                    const int lower_face = imu;
                    const int upper_face = imu + 1;
                    const double left_lower_mu_dot =
                        mu_dot_at(ix, iv, lower_face);
                    const double left_upper_mu_dot =
                        mu_dot_at(ix, iv, upper_face);
                    const double right_lower_mu_dot =
                        mu_dot_at(ix + 1, iv, lower_face);
                    const double right_upper_mu_dot =
                        mu_dot_at(ix + 1, iv, upper_face);
                    const double lower_char =
                        0.5 * (left_lower_mu_dot + right_lower_mu_dot);
                    const double upper_char =
                        0.5 * (left_upper_mu_dot + right_upper_mu_dot);
                    const double x_cell =
                        (static_cast<double>(sg.ix_start + ix) + 0.5)
                      * sg.dx;
                    audit.valid = 1;
                    audit.rank = mpi_rank;
                    audit.ix = sg.ix_start + ix;
                    audit.iv = iv;
                    audit.imu = imu;
                    audit.region =
                        (x_cell < CORE_DIAG_BOUNDARY_WIDTH ||
                         x_cell > Param::Lx - CORE_DIAG_BOUNDARY_WIDTH)
                        ? 1 : 0;
                    audit.severity = severity;
                    audit.f_input = f_before_mu;
                    audit.f_after_x = bkg_n.f[src] - dx_div;
                    audit.dx_div = dx_div;
                    audit.dmu_div_used = dmu_div_low;
                    audit.du_div_low = du_div;
                    audit.f_low = f_mu_low;
                    audit.left_lower_flux = mu_left_lower;
                    audit.left_upper_flux = mu_left_upper;
                    audit.right_lower_flux = mu_right_lower;
                    audit.right_upper_flux = mu_right_upper;
                    audit.left_lower_scale = left_lower_mu_dot;
                    audit.left_upper_scale = left_upper_mu_dot;
                    audit.right_lower_scale = right_lower_mu_dot;
                    audit.right_upper_scale = right_upper_mu_dot;
                    audit.left_lower_donor_f = 0.0;
                    audit.left_upper_donor_f = 0.0;
                    audit.right_lower_donor_f = 0.0;
                    audit.right_upper_donor_f = 0.0;
                    audit.lower_characteristic = lower_char;
                    audit.upper_characteristic = upper_char;
                    audit.moment_weight = bkg_n.vgrid.moment_weight[iv];
                    audit.cell_budget =
                        std::isfinite(f_mu_low) ? f_mu_low - f_floor : 0.0;
                    audit.left_lower_donor_index =
                        mu_donor_imu_at(left_lower_mu_dot, lower_face);
                    audit.left_upper_donor_index =
                        mu_donor_imu_at(left_upper_mu_dot, upper_face);
                    audit.right_lower_donor_index =
                        mu_donor_imu_at(right_lower_mu_dot, lower_face);
                    audit.right_upper_donor_index =
                        mu_donor_imu_at(right_upper_mu_dot, upper_face);
                }
            }
        }
        long long global_mu_low_order_failed =
            local_mu_low_order_failed;
        MPI_Allreduce(MPI_IN_PLACE, &global_mu_low_order_failed, 1,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        LowOrderFailureAudit local_mu_low_audit =
            empty_low_order_failure_audit();
        for (const LowOrderFailureAudit& audit : thread_mu_low_failure) {
            if (audit.valid &&
                (!local_mu_low_audit.valid ||
                 audit.severity > local_mu_low_audit.severity)) {
                local_mu_low_audit = audit;
            }
        }
        struct {
            double value;
            int rank;
        } local_mu_low_loc, global_mu_low_loc;
        local_mu_low_loc.value =
            local_mu_low_audit.valid ? local_mu_low_audit.severity : 0.0;
        local_mu_low_loc.rank = mpi_rank;
        global_mu_low_loc.value = 0.0;
        global_mu_low_loc.rank = -1;
        MPI_Allreduce(&local_mu_low_loc, &global_mu_low_loc, 1,
                      MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        double mu_low_values[27] = {0.0};
        int mu_low_indices[10] = {-1, -1, -1, -1, 0, -1, -1, -1, -1, 0};
        if (local_mu_low_audit.valid &&
            mpi_rank == global_mu_low_loc.rank) {
            mu_low_values[0] = local_mu_low_audit.severity;
            mu_low_values[1] = local_mu_low_audit.f_input;
            mu_low_values[2] = local_mu_low_audit.f_after_x;
            mu_low_values[3] = local_mu_low_audit.dx_div;
            mu_low_values[4] = local_mu_low_audit.dmu_div_used;
            mu_low_values[5] = local_mu_low_audit.du_div_low;
            mu_low_values[6] = local_mu_low_audit.f_low;
            mu_low_values[7] = local_mu_low_audit.left_lower_flux;
            mu_low_values[8] = local_mu_low_audit.left_upper_flux;
            mu_low_values[9] = local_mu_low_audit.right_lower_flux;
            mu_low_values[10] = local_mu_low_audit.right_upper_flux;
            mu_low_values[11] = local_mu_low_audit.left_lower_scale;
            mu_low_values[12] = local_mu_low_audit.left_upper_scale;
            mu_low_values[13] = local_mu_low_audit.right_lower_scale;
            mu_low_values[14] = local_mu_low_audit.right_upper_scale;
            mu_low_values[19] = local_mu_low_audit.lower_characteristic;
            mu_low_values[20] = local_mu_low_audit.upper_characteristic;
            mu_low_values[21] = local_mu_low_audit.moment_weight;
            mu_low_values[22] = local_mu_low_audit.cell_budget;
            mu_low_values[23] =
                static_cast<double>(global_mu_low_order_failed);
            mu_low_indices[0] = local_mu_low_audit.rank;
            mu_low_indices[1] = local_mu_low_audit.ix;
            mu_low_indices[2] = local_mu_low_audit.iv;
            mu_low_indices[3] = local_mu_low_audit.imu;
            mu_low_indices[4] = local_mu_low_audit.region;
            mu_low_indices[5] =
                local_mu_low_audit.left_lower_donor_index;
            mu_low_indices[6] =
                local_mu_low_audit.left_upper_donor_index;
            mu_low_indices[7] =
                local_mu_low_audit.right_lower_donor_index;
            mu_low_indices[8] =
                local_mu_low_audit.right_upper_donor_index;
        }
        if (global_mu_low_loc.value > 0.0 &&
            global_mu_low_loc.rank >= 0) {
            MPI_Bcast(mu_low_values, 27, MPI_DOUBLE,
                      global_mu_low_loc.rank, MPI_COMM_WORLD);
            MPI_Bcast(mu_low_indices, 10, MPI_INT,
                      global_mu_low_loc.rank, MPI_COMM_WORLD);
            LowOrderFailureAudit& audit = fluxes.mu_low_failure_audit;
            audit.valid = 1;
            audit.rank = mu_low_indices[0];
            audit.ix = mu_low_indices[1];
            audit.iv = mu_low_indices[2];
            audit.imu = mu_low_indices[3];
            audit.region = mu_low_indices[4];
            audit.severity = mu_low_values[0];
            audit.f_input = mu_low_values[1];
            audit.f_after_x = mu_low_values[2];
            audit.dx_div = mu_low_values[3];
            audit.dmu_div_used = mu_low_values[4];
            audit.du_div_low = mu_low_values[5];
            audit.f_low = mu_low_values[6];
            audit.left_lower_flux = mu_low_values[7];
            audit.left_upper_flux = mu_low_values[8];
            audit.right_lower_flux = mu_low_values[9];
            audit.right_upper_flux = mu_low_values[10];
            audit.left_lower_scale = mu_low_values[11];
            audit.left_upper_scale = mu_low_values[12];
            audit.right_lower_scale = mu_low_values[13];
            audit.right_upper_scale = mu_low_values[14];
            audit.lower_characteristic = mu_low_values[19];
            audit.upper_characteristic = mu_low_values[20];
            audit.moment_weight = mu_low_values[21];
            audit.cell_budget = mu_low_values[22];
            audit.low_order_failed_count = mu_low_values[23];
            audit.left_lower_donor_index = mu_low_indices[5];
            audit.left_upper_donor_index = mu_low_indices[6];
            audit.right_lower_donor_index = mu_low_indices[7];
            audit.right_upper_donor_index = mu_low_indices[8];
        }
    }
    if (fluxes.u_flux_audit_valid) {
        double final_values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const int audit_local_ix = fluxes.u_flux_audit_ix - sg.ix_start;
        const bool owns_audit =
            mpi_rank == fluxes.u_flux_audit_rank &&
            audit_local_ix >= 0 && audit_local_ix < nxl &&
            fluxes.u_flux_audit_iv >= 0 &&
            fluxes.u_flux_audit_iv < Param::Nv &&
            fluxes.u_flux_audit_imu >= 0 &&
            fluxes.u_flux_audit_imu < Param::Nmu;
        if (owns_audit) {
            const int ix = audit_local_ix;
            const int iv = fluxes.u_flux_audit_iv;
            const int imu = fluxes.u_flux_audit_imu;
            const size_t k =
                static_cast<size_t>(iv) * Param::Nmu
              + static_cast<size_t>(imu);
            const size_t dst =
                static_cast<size_t>(ng + ix) * Param::Nvmu + k;
            const double hdt_is = half_dt_inv_shell[iv];
            const double dx_div =
                dt_dx *
                (fluxes.x_final[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_final[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double u_xl_lo =
                fluxes.u_final_cell[u_cell_index(ix, iv, imu)];
            const double u_xl_hi =
                fluxes.u_final_cell[u_cell_index(ix, iv + 1, imu)];
            const double u_xr_lo = u_xl_lo;
            const double u_xr_hi = u_xl_hi;
            const double u_div_cell = u_xl_hi - u_xl_lo;
            const double du_div_final =
                (2.0 * hdt_is) * u_div_cell;
            const double mu_div_left =
                mu_force_face[mu_xface_index(ix, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix, iv, imu)];
            const double mu_div_right =
                mu_force_face[mu_xface_index(ix + 1, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix + 1, iv, imu)];
            const double dmu_div =
                hdt_is * (mu_div_left + mu_div_right);
            final_values[0] = du_div_final;
            final_values[1] =
                bkg_n.f[dst] - dx_div - du_div_final - dmu_div;
            final_values[2] = u_xl_lo;
            final_values[3] = u_xl_hi;
            final_values[4] = u_xr_lo;
            final_values[5] = u_xr_hi;
        }
        MPI_Bcast(final_values, 6, MPI_DOUBLE,
                  fluxes.u_flux_audit_rank, MPI_COMM_WORLD);
        fluxes.u_flux_audit_du_div_final = final_values[0];
        fluxes.u_flux_audit_updated = final_values[1];
        fluxes.u_flux_audit_final_xl_lo = final_values[2];
        fluxes.u_flux_audit_final_xl_hi = final_values[3];
        fluxes.u_flux_audit_final_xr_lo = final_values[4];
        fluxes.u_flux_audit_final_xr_hi = final_values[5];
    }
    {
        long long local_mu_counts[6] = {
            local_mu_low_u_face_active,
            local_mu_low_u_face_total,
            local_mu_low_u_face_active_boundary,
            local_mu_low_u_face_total_boundary,
            local_mu_low_u_face_active_core,
            local_mu_low_u_face_total_core
        };
        long long global_mu_counts[6] = {0, 0, 0, 0, 0, 0};
        MPI_Allreduce(local_mu_counts, global_mu_counts, 6,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        double local_mu_alpha[3] = {
            local_mu_low_u_alpha_min,
            local_mu_low_u_alpha_min_boundary,
            local_mu_low_u_alpha_min_core
        };
        double global_mu_alpha[3] = {1.0, 1.0, 1.0};
        MPI_Allreduce(local_mu_alpha, global_mu_alpha, 3, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        fluxes.mu_low_u_alpha_min = global_mu_alpha[0];
        fluxes.mu_low_u_alpha_min_boundary = global_mu_alpha[1];
        fluxes.mu_low_u_alpha_min_core = global_mu_alpha[2];
        fluxes.mu_low_u_limiter_active_fraction =
            (global_mu_counts[1] > 0)
            ? static_cast<double>(global_mu_counts[0]) /
              static_cast<double>(global_mu_counts[1])
            : 0.0;
        fluxes.mu_low_u_limiter_active_fraction_boundary =
            (global_mu_counts[3] > 0)
            ? static_cast<double>(global_mu_counts[2]) /
              static_cast<double>(global_mu_counts[3])
            : 0.0;
        fluxes.mu_low_u_limiter_active_fraction_core =
            (global_mu_counts[5] > 0)
            ? static_cast<double>(global_mu_counts[4]) /
              static_cast<double>(global_mu_counts[5])
            : 0.0;

        long long local_subcycle_counts[3] = {
            local_low_u_subcycle_active,
            local_low_u_subcycle_total,
            local_low_u_subcycle_sum
        };
        long long global_subcycle_counts[3] = {0, 0, 0};
        MPI_Allreduce(local_subcycle_counts, global_subcycle_counts, 3,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        int global_low_u_max_subcycles = local_low_u_max_subcycles;
        MPI_Allreduce(MPI_IN_PLACE, &global_low_u_max_subcycles, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
        fluxes.low_u_subcycle_active_fraction =
            (global_subcycle_counts[1] > 0)
            ? static_cast<double>(global_subcycle_counts[0]) /
              static_cast<double>(global_subcycle_counts[1])
            : 0.0;
        fluxes.low_u_average_subcycles =
            (global_subcycle_counts[1] > 0)
            ? static_cast<double>(global_subcycle_counts[2]) /
              static_cast<double>(global_subcycle_counts[1])
            : 1.0;
        fluxes.low_u_max_subcycles = global_low_u_max_subcycles;
    }

    // ---- Pass 3: final update (uses cell-local u_final_cell) ----
    int local_bad_update_centered = 0;
    int local_bad_negative_hard = 0;
    double local_neg_mass_defect = 0.0;
    double local_neg_energy_defect = 0.0;
    double local_limiter_mass_delta = 0.0;
    double local_limiter_momentum_delta = 0.0;
    double local_limiter_energy_delta = 0.0;
    double local_x_limiter_mass_delta = 0.0;
    double local_x_limiter_energy_delta = 0.0;
    double local_mu_low_u_energy_delta = 0.0;
    double local_mu_low_u_energy_delta_boundary = 0.0;
    double local_mu_low_u_energy_delta_core = 0.0;
    long long local_remap_cell_count = 0;
    long long local_remap_cell_total = 0;
    double local_limiter_energy_boundary_01 = 0.0;
    double local_limiter_energy_core_01 = 0.0;
    double local_abs_limiter_energy_boundary_01 = 0.0;
    double local_abs_limiter_energy_core_01 = 0.0;
    double local_limiter_energy_boundary_02 = 0.0;
    double local_limiter_energy_core_02 = 0.0;
    double local_abs_limiter_energy_boundary_02 = 0.0;
    double local_abs_limiter_energy_core_02 = 0.0;
    long long local_limiter_active_boundary_01 = 0;
    long long local_limiter_total_boundary_01 = 0;
    long long local_limiter_active_core_01 = 0;
    long long local_limiter_total_core_01 = 0;
    long long local_limiter_active_boundary_02 = 0;
    long long local_limiter_total_boundary_02 = 0;
    long long local_limiter_active_core_02 = 0;
    long long local_limiter_total_core_02 = 0;
    double local_f_neg_min = 0.0;
    double local_f_neg_ratio_max = 0.0;
    double local_f_neg_mass = 0.0;
    long long local_f_neg_count = 0;
    const int nthreads = std::max(1, omp_get_max_threads());
    std::vector<FiniteFluxCandidate> thread_flux_failures(
        static_cast<size_t>(nthreads), empty_finite_flux_candidate());
    struct NegCellInfo {
        double neg_ratio;
        double f_val;
        int ix;
        int iv;
        int imu;
        int valid;
    };
    std::vector<NegCellInfo> thread_worst_neg(
        static_cast<size_t>(nthreads), NegCellInfo{0.0, 0.0, -1, -1, -1, 0});
    struct StageThreadStats {
        std::vector<double> min_f;
        std::vector<double> neg_mass;
        std::vector<long long> neg_count;
        std::vector<double> low_u_neg_mass;
        std::vector<double> low_u_neg_added_by_div;
        std::vector<double> core_low_u_min_f;
        std::vector<double> min_f_core_by_u;
        std::vector<double> neg_mass_core_by_u;
        std::vector<long long> neg_count_core_by_u;
        std::vector<double> min_f_boundary_by_u;
        std::vector<double> neg_mass_boundary_by_u;
        std::vector<long long> neg_count_boundary_by_u;
    };
    std::vector<StageThreadStats> thread_stage(
        collect_stage_diagnostics ? static_cast<size_t>(nthreads) : 0);
    for (StageThreadStats& stats : thread_stage) {
        stats.min_f.assign(BKG_STAGE_COUNT,
                           std::numeric_limits<double>::infinity());
        stats.neg_mass.assign(BKG_STAGE_COUNT, 0.0);
        stats.neg_count.assign(BKG_STAGE_COUNT, 0);
        stats.low_u_neg_mass.assign(BKG_STAGE_COUNT, 0.0);
        stats.low_u_neg_added_by_div.assign(BKG_DIV_COMPONENT_COUNT, 0.0);
        stats.core_low_u_min_f.assign(
            BKG_STAGE_COUNT, std::numeric_limits<double>::infinity());
        stats.min_f_core_by_u.assign(
            BKG_STAGE_COUNT * Param::Nv,
            std::numeric_limits<double>::infinity());
        stats.neg_mass_core_by_u.assign(BKG_STAGE_COUNT * Param::Nv, 0.0);
        stats.neg_count_core_by_u.assign(BKG_STAGE_COUNT * Param::Nv, 0);
        stats.min_f_boundary_by_u.assign(
            BKG_STAGE_COUNT * Param::Nv,
            std::numeric_limits<double>::infinity());
        stats.neg_mass_boundary_by_u.assign(
            BKG_STAGE_COUNT * Param::Nv, 0.0);
        stats.neg_count_boundary_by_u.assign(
            BKG_STAGE_COUNT * Param::Nv, 0);
    }
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_neg_mass_defect,local_neg_energy_defect, \
                    local_limiter_mass_delta,local_limiter_momentum_delta, \
                    local_limiter_energy_delta, \
                    local_x_limiter_mass_delta, \
                    local_x_limiter_energy_delta, \
                    local_mu_low_u_energy_delta, \
                    local_mu_low_u_energy_delta_boundary, \
                    local_mu_low_u_energy_delta_core, \
                    local_remap_cell_count,local_remap_cell_total, \
                    local_limiter_energy_boundary_01, \
                    local_limiter_energy_core_01, \
                    local_abs_limiter_energy_boundary_01, \
                    local_abs_limiter_energy_core_01, \
                    local_limiter_energy_boundary_02, \
                    local_limiter_energy_core_02, \
                    local_abs_limiter_energy_boundary_02, \
                    local_abs_limiter_energy_core_02, \
                    local_limiter_active_boundary_01, \
                    local_limiter_total_boundary_01, \
                    local_limiter_active_core_01, \
                    local_limiter_total_core_01, \
                    local_limiter_active_boundary_02, \
                    local_limiter_total_boundary_02, \
                    local_limiter_active_core_02, \
                    local_limiter_total_core_02, \
                    local_f_neg_mass,local_f_neg_count) \
        reduction(min:local_f_neg_min) \
        reduction(max:local_f_neg_ratio_max)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int k_int = 0; k_int < static_cast<int>(Param::Nvmu); ++k_int) {
            const int iv = k_int / Param::Nmu;
            const int imu = k_int - iv * Param::Nmu;
            const size_t k = static_cast<size_t>(k_int);
            const size_t dst =
                static_cast<size_t>(ng + ix) * Param::Nvmu + k;
            const double hdt_is = half_dt_inv_shell[iv];
            const double dx_div =
                dt_dx *
                (fluxes.x_final[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_final[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double dx_div_high =
                dt_dx *
                (fluxes.x_high[static_cast<size_t>(ix + 1) * Param::Nvmu + k]
               - fluxes.x_high[static_cast<size_t>(ix) * Param::Nvmu + k]);
            const double u_div =
                fluxes.u_final_cell[u_cell_index(ix, iv + 1, imu)]
              - fluxes.u_final_cell[u_cell_index(ix, iv, imu)];
            const double du_div = (2.0 * hdt_is) * u_div;
            const double u_div_ec =
                fluxes.u_high_cell[u_cell_index(ix, iv + 1, imu)]
              - fluxes.u_high_cell[u_cell_index(ix, iv, imu)];
            const double du_div_ec = (2.0 * hdt_is) * u_div_ec;
            const double mu_div_left =
                mu_force_face[mu_xface_index(ix, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix, iv, imu)];
            const double mu_div_right =
                mu_force_face[mu_xface_index(ix + 1, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix + 1, iv, imu)];
            const double dmu_div = hdt_is * (mu_div_left + mu_div_right);
            double dmu_div_ec = dmu_div;
            if (iv < low_u_limit_count) {
                const size_t low_mu_left_base =
                    (static_cast<size_t>(ix) * low_u_limit_count
                   + static_cast<size_t>(iv)) *
                    static_cast<size_t>(Param::Nmu + 1);
                const size_t low_mu_right_base =
                    (static_cast<size_t>(ix + 1) * low_u_limit_count
                   + static_cast<size_t>(iv)) *
                    static_cast<size_t>(Param::Nmu + 1);
                const double mu_div_left_ec =
                    mu_force_low_u_ec[low_mu_left_base
                                      + static_cast<size_t>(imu + 1)]
                  - mu_force_low_u_ec[low_mu_left_base
                                      + static_cast<size_t>(imu)];
                const double mu_div_right_ec =
                    mu_force_low_u_ec[low_mu_right_base
                                      + static_cast<size_t>(imu + 1)]
                  - mu_force_low_u_ec[low_mu_right_base
                                      + static_cast<size_t>(imu)];
                dmu_div_ec = hdt_is * (mu_div_left_ec + mu_div_right_ec);
            }
            const double updated = bkg_n.f[dst] - dx_div - du_div - dmu_div;
            const double f0 = bkg_n.f[dst];
            const double local_scale = std::max(1.0, std::fabs(f0));
            const double cell_weight = bkg_n.vgrid.moment_weight[iv];
            const double ke_per_mass = ke_per_mass_arr[iv];
            const double x_cell =
                (static_cast<double>(sg.ix_start + ix) + 0.5) * sg.dx;
            const bool in_boundary_01 =
                (x_cell < 0.1 * Const::micro) ||
                (x_cell > Param::Lx - 0.1 * Const::micro);
            const bool in_boundary_02 =
                (x_cell < 0.2 * Const::micro) ||
                (x_cell > Param::Lx - 0.2 * Const::micro);
            if (collect_stage_diagnostics) {
                const double stage_f[BKG_STAGE_COUNT] = {
                    f0 - dx_div,
                    f0 - dx_div - du_div,
                    updated
                };
                StageThreadStats& stats =
                    thread_stage[static_cast<size_t>(omp_get_thread_num())];
                if (iv == 0) {
                    const double f0_neg = std::max(0.0, -f0);
                    const double after_dx_neg = std::max(0.0, -stage_f[0]);
                    const double after_du_neg = std::max(0.0, -stage_f[1]);
                    const double after_dmu_neg = std::max(0.0, -stage_f[2]);
                    const double neg_dx_added =
                        std::max(0.0, after_dx_neg - f0_neg);
                    const double neg_du_added =
                        std::max(0.0, after_du_neg - after_dx_neg);
                    const double neg_dmu_added =
                        std::max(0.0, after_dmu_neg - after_du_neg);
                    stats.low_u_neg_added_by_div[0] +=
                        neg_dx_added * cell_weight * sg.dx;
                    stats.low_u_neg_added_by_div[1] +=
                        neg_du_added * cell_weight * sg.dx;
                    stats.low_u_neg_added_by_div[2] +=
                        neg_dmu_added * cell_weight * sg.dx;
                }
                for (int istage = 0; istage < BKG_STAGE_COUNT; ++istage) {
                    const double fv = stage_f[istage];
                    const size_t by_u =
                        static_cast<size_t>(istage * Param::Nv + iv);
                    stats.min_f[static_cast<size_t>(istage)] =
                        std::min(stats.min_f[static_cast<size_t>(istage)], fv);
                    if (iv == 0 && fv < 0.0) {
                        stats.low_u_neg_mass[static_cast<size_t>(istage)] +=
                            (-fv) * cell_weight * sg.dx;
                    }
                    if (in_boundary_02) {
                        stats.min_f_boundary_by_u[by_u] =
                            std::min(stats.min_f_boundary_by_u[by_u], fv);
                    } else {
                        stats.core_low_u_min_f[static_cast<size_t>(istage)] =
                            (iv == 0)
                            ? std::min(
                                stats.core_low_u_min_f[
                                    static_cast<size_t>(istage)], fv)
                            : stats.core_low_u_min_f[
                                static_cast<size_t>(istage)];
                        stats.min_f_core_by_u[by_u] =
                            std::min(stats.min_f_core_by_u[by_u], fv);
                    }
                    if (fv < 0.0) {
                        const double neg_mass = -fv * cell_weight * sg.dx;
                        stats.neg_mass[static_cast<size_t>(istage)] +=
                            neg_mass;
                        ++stats.neg_count[static_cast<size_t>(istage)];
                        if (in_boundary_02) {
                            stats.neg_mass_boundary_by_u[by_u] += neg_mass;
                            ++stats.neg_count_boundary_by_u[by_u];
                        } else {
                            stats.neg_mass_core_by_u[by_u] += neg_mass;
                            ++stats.neg_count_core_by_u[by_u];
                        }
                    }
                }
            }
            const bool limiter_cell_active =
                fluxes.cell_alpha_u[static_cast<size_t>(ix) * Param::Nvmu + k]
                < 0.999999;
            if (in_boundary_01) {
                ++local_limiter_total_boundary_01;
                if (limiter_cell_active) ++local_limiter_active_boundary_01;
            } else {
                ++local_limiter_total_core_01;
                if (limiter_cell_active) ++local_limiter_active_core_01;
            }
            if (in_boundary_02) {
                ++local_limiter_total_boundary_02;
                if (limiter_cell_active) ++local_limiter_active_boundary_02;
            } else {
                ++local_limiter_total_core_02;
                if (limiter_cell_active) ++local_limiter_active_core_02;
            }
            const double limiter_delta_f = du_div_ec - du_div;
            const double x_limiter_delta_f = dx_div_high - dx_div;
            if (x_limiter_delta_f != 0.0) {
                const double x_limiter_delta_n =
                    x_limiter_delta_f * cell_weight * sg.dx;
                local_x_limiter_mass_delta += x_limiter_delta_n;
                local_x_limiter_energy_delta +=
                    x_limiter_delta_n * ke_per_mass;
            }
            const double mu_limiter_delta_f = dmu_div_ec - dmu_div;
            if (mu_limiter_delta_f != 0.0 && iv < low_u_limit_count) {
                const double mu_limiter_delta_e =
                    mu_limiter_delta_f * cell_weight * sg.dx * ke_per_mass;
                local_mu_low_u_energy_delta += mu_limiter_delta_e;
                if (in_boundary_02) {
                    local_mu_low_u_energy_delta_boundary +=
                        mu_limiter_delta_e;
                } else {
                    local_mu_low_u_energy_delta_core +=
                        mu_limiter_delta_e;
                }
            }
            if (limiter_delta_f != 0.0) {
                const double limiter_delta_n =
                    limiter_delta_f * cell_weight * sg.dx;
                const double limiter_delta_e =
                    limiter_delta_n * ke_per_mass;
                const double px =
                    bkg_n.mass * Const::c * bkg_n.vgrid.v_cells[iv] *
                    bkg_n.vgrid.mu_cells[imu];
                local_limiter_mass_delta += limiter_delta_n;
                local_limiter_momentum_delta += limiter_delta_n * px;
                local_limiter_energy_delta += limiter_delta_e;
                if (in_boundary_01) {
                    local_limiter_energy_boundary_01 += limiter_delta_e;
                    local_abs_limiter_energy_boundary_01 +=
                        std::fabs(limiter_delta_e);
                } else {
                    local_limiter_energy_core_01 += limiter_delta_e;
                    local_abs_limiter_energy_core_01 +=
                        std::fabs(limiter_delta_e);
                }
                if (in_boundary_02) {
                    local_limiter_energy_boundary_02 += limiter_delta_e;
                    local_abs_limiter_energy_boundary_02 +=
                        std::fabs(limiter_delta_e);
                } else {
                    local_limiter_energy_core_02 += limiter_delta_e;
                    local_abs_limiter_energy_core_02 +=
                        std::fabs(limiter_delta_e);
                }
            }

            bool record_failure = false;
            double severity = 0.0;
            double max_negative = 0.0;
            double relative_negative = 0.0;

            if (!std::isfinite(updated)) {
                #pragma omp atomic write
                local_bad_update_centered = 1;
                record_failure = true;
                severity = std::numeric_limits<double>::infinity();
                max_negative = std::numeric_limits<double>::infinity();
                relative_negative = std::numeric_limits<double>::infinity();
            } else if (updated < 0.0 && iv >= low_u_limit_count) {
                const double neg_ratio = -updated / local_scale;
                if (neg_ratio >= NEG_TOL_HARD) {
                    #pragma omp atomic write
                    local_bad_negative_hard = 1;
                    record_failure = true;
                    max_negative = -updated;
                    relative_negative = neg_ratio;
                    severity = max_negative;
                } else if (neg_ratio >= NEG_TOL_SOFT) {
                    const double mass_defect =
                        (-updated) * cell_weight * sg.dx;
                    local_neg_mass_defect += mass_defect;
                    local_neg_energy_defect += mass_defect * ke_per_mass;
                }
                // neg_ratio < NEG_TOL_SOFT: negligible roundoff, silent
            }

            // Track f-negativity for diagnostics (all levels)
            if (updated < 0.0 && iv >= low_u_limit_count) {
                local_f_neg_min = std::min(local_f_neg_min, updated);
                const double neg_ratio = (updated < 0.0)
                    ? (-updated) / local_scale : 0.0;
                local_f_neg_ratio_max =
                    std::max(local_f_neg_ratio_max, neg_ratio);
                const double neg_mass =
                    (-updated) * cell_weight * sg.dx;
                local_f_neg_mass += neg_mass;
                ++local_f_neg_count;
                const int tid = omp_get_thread_num();
                NegCellInfo& wn =
                    thread_worst_neg[static_cast<size_t>(tid)];
                if (!wn.valid || neg_ratio > wn.neg_ratio) {
                    wn.neg_ratio = neg_ratio;
                    wn.f_val = updated;
                    wn.ix = sg.ix_start + ix;
                    wn.iv = iv;
                    wn.imu = imu;
                    wn.valid = 1;
                }
            }

            if (record_failure) {
                const int tid = omp_get_thread_num();
                FiniteFluxCandidate& candidate =
                    thread_flux_failures[static_cast<size_t>(tid)];
                if (!candidate.has_failure ||
                    severity > candidate.severity) {
                    candidate.severity = severity;
                    candidate.max_negative = max_negative;
                    candidate.relative_negative = relative_negative;
                    candidate.updated = updated;
                    candidate.f0 = f0;
                    candidate.dx_div = dx_div;
                    candidate.du_div = du_div;
                    candidate.dmu_div = dmu_div;
                    candidate.ix = sg.ix_start + ix;
                    candidate.iv = iv;
                    candidate.imu = imu;
                    candidate.has_failure = 1;
                }
            }

            bkg_new.f[dst] = updated;
        }
    }

    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_mu_low_u_energy_delta, \
                    local_mu_low_u_energy_delta_boundary, \
                    local_mu_low_u_energy_delta_core, \
                    local_neg_mass_defect,local_neg_energy_defect, \
                    local_f_neg_mass,local_f_neg_count) \
        reduction(min:local_f_neg_min) \
        reduction(max:local_f_neg_ratio_max)
    for (int ix = 0; ix < nxl; ++ix) {
        for (int iv = 0; iv < low_u_limit_count; ++iv) {
            ++local_remap_cell_total;
            const double cell_weight = bkg_n.vgrid.moment_weight[iv];
            const double ke_per_mass = ke_per_mass_arr[iv];
            /*
             * Do not repair endpoint negatives with max(0,f) scaling.  7.1.3
             * requires positivity control to enter through conservative
             * fluxes.  Strong endpoint negatives remain in the existing
             * negative-mass diagnostics so the preceding u/mu flux limiter can
             * be fixed instead of hiding the defect with a post-update clip.
             */
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t dst =
                    static_cast<size_t>(ng + ix) * Param::Nvmu
                  + static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(imu);
                const double fv = bkg_new.f[dst];
                if (!std::isfinite(fv)) {
                    #pragma omp atomic write
                    local_bad_update_centered = 1;
                    continue;
                }
                if (fv >= 0.0) continue;

                const double local_scale = std::max(1.0, std::fabs(fv));
                const double neg_ratio = -fv / local_scale;
                if (neg_ratio >= NEG_TOL_SOFT) {
                    const double mass_defect =
                        (-fv) * cell_weight * sg.dx;
                    local_neg_mass_defect += mass_defect;
                    local_neg_energy_defect += mass_defect * ke_per_mass;
                }
                local_f_neg_min = std::min(local_f_neg_min, fv);
                local_f_neg_ratio_max =
                    std::max(local_f_neg_ratio_max, neg_ratio);
                local_f_neg_mass += (-fv) * cell_weight * sg.dx;
                ++local_f_neg_count;

                const int tid = omp_get_thread_num();
                NegCellInfo& wn =
                    thread_worst_neg[static_cast<size_t>(tid)];
                if (!wn.valid || neg_ratio > wn.neg_ratio) {
                    wn.neg_ratio = neg_ratio;
                    wn.f_val = fv;
                    wn.ix = sg.ix_start + ix;
                    wn.iv = iv;
                    wn.imu = imu;
                    wn.valid = 1;
                }
            }
        }
    }

    if (collect_stage_diagnostics) {
        std::vector<double> local_stage_min_f(
            BKG_STAGE_COUNT, std::numeric_limits<double>::infinity());
        std::vector<double> local_stage_neg_mass(BKG_STAGE_COUNT, 0.0);
        std::vector<long long> local_stage_neg_count(BKG_STAGE_COUNT, 0);
        std::vector<double> local_stage_low_u_neg_mass(BKG_STAGE_COUNT, 0.0);
        std::vector<double> local_low_u_neg_added_by_div(
            BKG_DIV_COMPONENT_COUNT, 0.0);
        std::vector<double> local_stage_core_low_u_min_f(
            BKG_STAGE_COUNT, std::numeric_limits<double>::infinity());
        std::vector<double> local_min_f_core_by_u(
            BKG_STAGE_COUNT * Param::Nv,
            std::numeric_limits<double>::infinity());
        std::vector<double> local_neg_mass_core_by_u(
            BKG_STAGE_COUNT * Param::Nv, 0.0);
        std::vector<long long> local_neg_count_core_by_u(
            BKG_STAGE_COUNT * Param::Nv, 0);
        std::vector<double> local_min_f_boundary_by_u(
            BKG_STAGE_COUNT * Param::Nv,
            std::numeric_limits<double>::infinity());
        std::vector<double> local_neg_mass_boundary_by_u(
            BKG_STAGE_COUNT * Param::Nv, 0.0);
        std::vector<long long> local_neg_count_boundary_by_u(
            BKG_STAGE_COUNT * Param::Nv, 0);

        for (const StageThreadStats& stats : thread_stage) {
            for (int istage = 0; istage < BKG_STAGE_COUNT; ++istage) {
                const size_t s = static_cast<size_t>(istage);
                local_stage_min_f[s] =
                    std::min(local_stage_min_f[s], stats.min_f[s]);
                local_stage_neg_mass[s] += stats.neg_mass[s];
                local_stage_neg_count[s] += stats.neg_count[s];
                local_stage_low_u_neg_mass[s] += stats.low_u_neg_mass[s];
                local_stage_core_low_u_min_f[s] =
                    std::min(local_stage_core_low_u_min_f[s],
                             stats.core_low_u_min_f[s]);
            }
            for (int ic = 0; ic < BKG_DIV_COMPONENT_COUNT; ++ic) {
                const size_t c = static_cast<size_t>(ic);
                local_low_u_neg_added_by_div[c] +=
                    stats.low_u_neg_added_by_div[c];
            }
            for (size_t i = 0; i < local_min_f_core_by_u.size(); ++i) {
                local_min_f_core_by_u[i] =
                    std::min(local_min_f_core_by_u[i],
                             stats.min_f_core_by_u[i]);
                local_neg_mass_core_by_u[i] += stats.neg_mass_core_by_u[i];
                local_neg_count_core_by_u[i] += stats.neg_count_core_by_u[i];
                local_min_f_boundary_by_u[i] =
                    std::min(local_min_f_boundary_by_u[i],
                             stats.min_f_boundary_by_u[i]);
                local_neg_mass_boundary_by_u[i] +=
                    stats.neg_mass_boundary_by_u[i];
                local_neg_count_boundary_by_u[i] +=
                    stats.neg_count_boundary_by_u[i];
            }
        }

        fluxes.stage_min_f = local_stage_min_f;
        fluxes.stage_neg_mass = local_stage_neg_mass;
        fluxes.stage_neg_cell_count = local_stage_neg_count;
        fluxes.stage_low_u_neg_mass = local_stage_low_u_neg_mass;
        fluxes.low_u_neg_added_by_div = local_low_u_neg_added_by_div;
        fluxes.stage_core_low_u_min_f = local_stage_core_low_u_min_f;
        fluxes.stage_min_f_core_by_u = local_min_f_core_by_u;
        fluxes.stage_neg_mass_core_by_u = local_neg_mass_core_by_u;
        fluxes.stage_neg_cell_count_core_by_u = local_neg_count_core_by_u;
        fluxes.stage_min_f_boundary_by_u = local_min_f_boundary_by_u;
        fluxes.stage_neg_mass_boundary_by_u = local_neg_mass_boundary_by_u;
        fluxes.stage_neg_cell_count_boundary_by_u =
            local_neg_count_boundary_by_u;
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_min_f.data(),
                      BKG_STAGE_COUNT, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_neg_mass.data(),
                      BKG_STAGE_COUNT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_neg_cell_count.data(),
                      BKG_STAGE_COUNT, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_low_u_neg_mass.data(),
                      BKG_STAGE_COUNT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.low_u_neg_added_by_div.data(),
                      BKG_DIV_COMPONENT_COUNT, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_core_low_u_min_f.data(),
                      BKG_STAGE_COUNT, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        const int by_u_count = BKG_STAGE_COUNT * Param::Nv;
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_min_f_core_by_u.data(),
                      by_u_count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_neg_mass_core_by_u.data(),
                      by_u_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,
                      fluxes.stage_neg_cell_count_core_by_u.data(),
                      by_u_count, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_min_f_boundary_by_u.data(),
                      by_u_count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, fluxes.stage_neg_mass_boundary_by_u.data(),
                      by_u_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,
                      fluxes.stage_neg_cell_count_boundary_by_u.data(),
                      by_u_count, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        for (double& value : fluxes.stage_min_f) {
            if (!std::isfinite(value)) value = 0.0;
        }
        for (double& value : fluxes.stage_core_low_u_min_f) {
            if (!std::isfinite(value)) value = 0.0;
        }
        for (double& value : fluxes.stage_min_f_core_by_u) {
            if (!std::isfinite(value)) value = 0.0;
        }
        for (double& value : fluxes.stage_min_f_boundary_by_u) {
            if (!std::isfinite(value)) value = 0.0;
        }
    }

    // 7.1.6: populate u-direction and mu-direction flux diagnostics
    // after stage_min_f is finalized via MPI reductions above.
    {
        double global_u_min_f_low = local_u_min_f_low;
        MPI_Allreduce(MPI_IN_PLACE, &global_u_min_f_low, 1,
                      MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        long long global_u_low_order_failed = local_u_low_order_failed;
        MPI_Allreduce(MPI_IN_PLACE, &global_u_low_order_failed, 1,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        // --- u-direction (index 1) ---
        FluxPositivityDiag& fu = fluxes.flux_pos[1];
        fu.min_f_before = fluxes.stage_min_f[0];  // after x, before u
        fu.min_f_low    =
            std::isfinite(global_u_min_f_low) ? global_u_min_f_low : 0.0;
        fu.min_f_final  = fluxes.stage_min_f[1];  // after x + u
        fu.low_order_failed_count =
            static_cast<double>(global_u_low_order_failed);
        fu.alpha_active_fraction   = fluxes.u_force_alpha_active_frac;
        fu.alpha_min               = fluxes.u_force_alpha_min;
        fu.alpha_core_fraction     = 0.0;
        fu.alpha_boundary_fraction = 0.0;
        fu.negative_mass_prevented = 0.0;
        FluxDefectDiag& du = fluxes.flux_defect[1];
        du.mass_defect       = fluxes.limiter_mass_defect;
        du.momentum_defect   = fluxes.limiter_momentum_defect;
        du.energy_defect     = fluxes.limiter_energy_defect;
        du.boundary_mass_loss   = 0.0;
        du.boundary_energy_loss = 0.0;

        // --- mu-direction (index 2) ---
        FluxPositivityDiag& fm = fluxes.flux_pos[2];
        fm.min_f_before = fluxes.stage_min_f[1];  // after x+u, before mu
        fm.min_f_low    = fluxes.stage_min_f[1];  // no mu_low yet → use final
        fm.min_f_final  = fluxes.stage_min_f[2];  // final
        fm.low_order_failed_count  = 0.0;  // not yet tracked for mu
        fm.alpha_active_fraction   = fluxes.mu_low_u_limiter_active_fraction;
        fm.alpha_min               = fluxes.mu_low_u_alpha_min;
        fm.alpha_core_fraction = fluxes.mu_low_u_limiter_active_fraction_core;
        fm.alpha_boundary_fraction =
            fluxes.mu_low_u_limiter_active_fraction_boundary;
        fm.negative_mass_prevented = 0.0;
        FluxDefectDiag& dm = fluxes.flux_defect[2];
        dm.mass_defect       = 0.0;
        dm.momentum_defect   = 0.0;
        dm.energy_defect     = fluxes.mu_low_u_energy_delta;
        dm.boundary_mass_loss   = 0.0;
        dm.boundary_energy_loss = 0.0;
    }

    // Reduce per-thread worst-negative
    NegCellInfo local_worst_neg = {0.0, 0.0, -1, -1, -1, 0};
    for (const NegCellInfo& wn : thread_worst_neg) {
        if (wn.valid &&
            (!local_worst_neg.valid ||
             wn.neg_ratio > local_worst_neg.neg_ratio)) {
            local_worst_neg = wn;
        }
    }

    FiniteFluxCandidate local_failure = empty_finite_flux_candidate();
    for (const FiniteFluxCandidate& candidate : thread_flux_failures) {
        if (candidate.has_failure &&
            (!local_failure.has_failure ||
             candidate.severity > local_failure.severity)) {
            local_failure = candidate;
        }
    }
    struct {
        double value;
        int rank;
    } local_worst, global_worst;
    local_worst.value =
        local_failure.has_failure ? local_failure.severity : 0.0;
    local_worst.rank = mpi_rank;
    global_worst.value = 0.0;
    global_worst.rank = -1;
    MPI_Allreduce(&local_worst, &global_worst, 1, MPI_DOUBLE_INT,
                  MPI_MAXLOC, MPI_COMM_WORLD);
    double failure_values[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int failure_indices[4] = {-1, -1, -1, -1};
    if (local_failure.has_failure && mpi_rank == global_worst.rank) {
        failure_values[0] = local_failure.max_negative;
        failure_values[1] = local_failure.relative_negative;
        failure_values[2] = local_failure.updated;
        failure_values[3] = local_failure.f0;
        failure_values[4] = local_failure.dx_div;
        failure_values[5] = local_failure.du_div;
        failure_values[6] = local_failure.dmu_div;
        failure_indices[0] = mpi_rank;
        failure_indices[1] = local_failure.ix;
        failure_indices[2] = local_failure.iv;
        failure_indices[3] = local_failure.imu;
    }
    if (global_worst.rank >= 0 && global_worst.value > 0.0) {
        MPI_Bcast(failure_values, 7, MPI_DOUBLE, global_worst.rank,
                  MPI_COMM_WORLD);
        MPI_Bcast(failure_indices, 4, MPI_INT, global_worst.rank,
                  MPI_COMM_WORLD);
        fluxes.finite_flux_has_failure = 1;
        fluxes.finite_flux_max_negative = failure_values[0];
        fluxes.finite_flux_relative_negative = failure_values[1];
        fluxes.finite_flux_updated = failure_values[2];
        fluxes.finite_flux_f0 = failure_values[3];
        fluxes.finite_flux_dx_div = failure_values[4];
        fluxes.finite_flux_du_div = failure_values[5];
        fluxes.finite_flux_dmu_div = failure_values[6];
        fluxes.finite_flux_rank = failure_indices[0];
        fluxes.finite_flux_ix = failure_indices[1];
        fluxes.finite_flux_iv = failure_indices[2];
        fluxes.finite_flux_imu = failure_indices[3];
    }
    if (local_bad_update_centered != 0) finite = false;
    // Negative-f protection remains disabled for long-run diagnostics.
    // if (local_bad_negative_hard != 0) finite = false;
    (void)local_bad_negative_hard;
    exchange_ghosts_x_persistent(bkg_new, sg, mpi_rank, mpi_size);
    update_flux_current(bkg_n, sg, fluxes, bkg_new);

    MPI_Allreduce(&local_negative_mass, &fluxes.negative_mass_before_repair,
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    {
        double local_pos[2] = {
            local_neg_mass_defect,
            local_neg_energy_defect
        };
        double global_pos[2] = { 0.0, 0.0 };
        MPI_Allreduce(local_pos, global_pos, 2, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        fluxes.mass_added_by_positivity_repair = 0.0;
        fluxes.positivity_mass_defect = global_pos[0];
        fluxes.positivity_energy_defect = global_pos[1];
    }
    {
        double local_limiter[3] = {
            local_limiter_mass_delta,
            local_limiter_momentum_delta,
            local_limiter_energy_delta
        };
        double global_limiter[3] = { 0.0, 0.0, 0.0 };
        MPI_Allreduce(local_limiter, global_limiter, 3, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        fluxes.limiter_mass_defect = global_limiter[0];
        fluxes.limiter_momentum_defect = global_limiter[1];
        fluxes.limiter_energy_defect = global_limiter[2];
        double local_x_limiter[2] = {
            local_x_limiter_mass_delta,
            local_x_limiter_energy_delta
        };
        double global_x_limiter[2] = {0.0, 0.0};
        MPI_Allreduce(local_x_limiter, global_x_limiter, 2,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        fluxes.x_limiter_mass_defect = global_x_limiter[0];
        fluxes.x_limiter_energy_defect = global_x_limiter[1];
        // 7.1.6: populate x-direction flux-positivity diagnostics
        {
            double local_x_low_counts[2] = {
                local_x_low_input_failed,
                local_x_low_order_failed
            };
            double global_x_low_counts[2] = {0.0, 0.0};
            MPI_Allreduce(local_x_low_counts, global_x_low_counts, 2,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double local_x_low_debt[4] = {
                local_x_low_input_neg_mass,
                local_x_low_input_rel_neg,
                local_x_low_output_rel_neg,
                local_x_low_input_core_failed
            };
            double global_x_low_debt[4] = {0.0, 0.0, 0.0, 0.0};
            MPI_Allreduce(&local_x_low_debt[0], &global_x_low_debt[0], 1,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_x_low_debt[1], &global_x_low_debt[1], 2,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&local_x_low_debt[3], &global_x_low_debt[3], 1,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double global_x_low_max_cfl = local_x_low_max_cfl;
            MPI_Allreduce(MPI_IN_PLACE, &global_x_low_max_cfl, 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            double local_x_pos[3] = {
                local_x_min_f_before,
                local_x_min_f_low,
                local_x_neg_mass_prevented
            };
            double global_x_pos[3] = {
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                0.0
            };
            MPI_Allreduce(local_x_pos, global_x_pos, 2, MPI_DOUBLE,
                          MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&local_x_pos[2], &global_x_pos[2], 1,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            /*
             * Classify x-low failures:
             * A: input state violates the shell-scaled tolerance.  Tiny,
             *    localized, non-core residual negatives are recorded as
             *    X_LOW_INPUT_DEBT and allowed to continue; larger debt stays
             *    a hard X_LOW_INPUT_BAD.
             * B: input is acceptable and CFL is monotone, so donor/periodic
             *    face construction is inconsistent.
             * C: input is acceptable and max CFL exceeds one; only this case
             *    should be retried with a smaller true x substep.
             */
            const double global_x_low_input_bad = global_x_low_counts[0];
            const double global_x_low_output_bad = global_x_low_counts[1];
            const double global_x_low_hard_count =
                global_x_low_input_bad + global_x_low_output_bad;
            fluxes.x_low_order_failed_count = global_x_low_hard_count;
            fluxes.x_low_failed_count = global_x_low_hard_count;
            fluxes.x_low_input_neg_mass = global_x_low_debt[0];
            fluxes.x_low_input_rel_neg = global_x_low_debt[1];
            fluxes.x_low_output_rel_neg = global_x_low_debt[2];
            fluxes.x_low_input_core_failed_count = global_x_low_debt[3];
            fluxes.x_low_input_debt_accepted = 0.0;
            if (global_x_low_input_bad > 0.0) {
                const double total_cells =
                    static_cast<double>(std::max(1, sg.nx_global))
                  * static_cast<double>(Param::Nvmu);
                const double soft_cell_limit =
                    std::max(16.0,
                             X_LOW_INPUT_DEBT_CELL_FRAC * total_cells);
                const double total_bkg_column =
                    std::max(1.0, Param::dens * Param::Lx);
                const double neg_mass_fraction =
                    global_x_low_debt[0] / total_bkg_column;
                const bool small_count =
                    global_x_low_input_bad <= soft_cell_limit;
                const bool small_relative =
                    global_x_low_debt[1] <= X_LOW_INPUT_DEBT_REL_TOL;
                const bool small_neg_mass =
                    neg_mass_fraction <= X_LOW_INPUT_DEBT_NEG_MASS_FRAC;
                const bool not_core =
                    global_x_low_debt[3] <= 0.0;
                const bool low_output_roundoff =
                    global_x_low_output_bad <= 0.0 &&
                    global_x_low_debt[2] <= X_LOW_INPUT_DEBT_OUTPUT_REL_TOL;
                if (small_count && small_relative && small_neg_mass &&
                    not_core && low_output_roundoff) {
                    fluxes.x_low_failure_kind = X_LOW_INPUT_DEBT;
                    fluxes.x_low_input_debt_accepted = 1.0;
                } else {
                    fluxes.x_low_failure_kind = X_LOW_INPUT_BAD;
                    finite = false;
                }
            } else if (global_x_low_output_bad > 0.0) {
                if (global_x_low_max_cfl <= 1.0 + X_LOW_CFL_TOL) {
                    fluxes.x_low_failure_kind = X_LOW_DONOR_BUG;
                } else {
                    fluxes.x_low_failure_kind = X_LOW_TRUE_CFL;
                }
                finite = false;
            } else {
                fluxes.x_low_failure_kind = X_LOW_OK;
            }
            /*
             * Keep the existing per-direction diagnostic fields populated
             * from the same classified quantities.
             */
            fluxes.x_low_input_min_f = global_x_pos[0];
            fluxes.x_low_max_cfl = global_x_low_max_cfl;
            fluxes.x_low_output_min_f = global_x_pos[1];
            FluxPositivityDiag& fx = fluxes.flux_pos[0];
            fx.min_f_before = global_x_pos[0];
            fx.min_f_low    = global_x_pos[1];
            // min_f_final from stage diagnostics (populated in Pass 3)
            fx.low_order_failed_count  = fluxes.x_low_order_failed_count;
            fx.alpha_active_fraction   = fluxes.limiter_active_fraction;
            fx.alpha_min               = fluxes.limiter_min_alpha;
            fx.alpha_core_fraction     = fluxes.limiter_active_fraction_core;
            fx.alpha_boundary_fraction = fluxes.limiter_active_fraction_boundary;
            fx.negative_mass_prevented = global_x_pos[2];
            FluxDefectDiag& dx = fluxes.flux_defect[0];
            dx.mass_defect     = fluxes.x_limiter_mass_defect;
            dx.momentum_defect = 0.0;  // x-direction has no momentum defect
            dx.energy_defect   = fluxes.x_limiter_energy_defect;
            dx.boundary_mass_loss = 0.0;
            dx.boundary_energy_loss = 0.0;
        }
        double local_mu_low_u_limiter[3] = {
            local_mu_low_u_energy_delta,
            local_mu_low_u_energy_delta_boundary,
            local_mu_low_u_energy_delta_core
        };
        double global_mu_low_u_limiter[3] = {0.0, 0.0, 0.0};
        MPI_Allreduce(local_mu_low_u_limiter, global_mu_low_u_limiter, 3,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        fluxes.mu_low_u_energy_delta = global_mu_low_u_limiter[0];
        fluxes.mu_low_u_energy_delta_boundary =
            global_mu_low_u_limiter[1];
        fluxes.mu_low_u_energy_delta_core = global_mu_low_u_limiter[2];

        long long local_remap_counts[2] = {
            local_remap_cell_count,
            local_remap_cell_total
        };
        long long global_remap_counts[2] = {0, 0};
        MPI_Allreduce(local_remap_counts, global_remap_counts, 2,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        fluxes.remap_cell_count = global_remap_counts[0];
        fluxes.remap_active_fraction =
            (global_remap_counts[1] > 0)
            ? static_cast<double>(global_remap_counts[0]) /
              static_cast<double>(global_remap_counts[1])
            : 0.0;
    }
    {
        double local_region_limiter[8] = {
            local_limiter_energy_boundary_01,
            local_limiter_energy_core_01,
            local_abs_limiter_energy_boundary_01,
            local_abs_limiter_energy_core_01,
            local_limiter_energy_boundary_02,
            local_limiter_energy_core_02,
            local_abs_limiter_energy_boundary_02,
            local_abs_limiter_energy_core_02
        };
        double global_region_limiter[8] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        MPI_Allreduce(local_region_limiter, global_region_limiter, 8,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        fluxes.region_u_limiter_energy_boundary[0] =
            global_region_limiter[0];
        fluxes.region_u_limiter_energy_core[0] = global_region_limiter[1];
        fluxes.region_abs_u_limiter_energy_boundary[0] =
            global_region_limiter[2];
        fluxes.region_abs_u_limiter_energy_core[0] =
            global_region_limiter[3];
        fluxes.region_u_limiter_energy_boundary[1] =
            global_region_limiter[4];
        fluxes.region_u_limiter_energy_core[1] = global_region_limiter[5];
        fluxes.region_abs_u_limiter_energy_boundary[1] =
            global_region_limiter[6];
        fluxes.region_abs_u_limiter_energy_core[1] =
            global_region_limiter[7];

        long long local_region_counts[8] = {
            local_limiter_active_boundary_01,
            local_limiter_total_boundary_01,
            local_limiter_active_core_01,
            local_limiter_total_core_01,
            local_limiter_active_boundary_02,
            local_limiter_total_boundary_02,
            local_limiter_active_core_02,
            local_limiter_total_core_02
        };
        long long global_region_counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        MPI_Allreduce(local_region_counts, global_region_counts, 8,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        fluxes.region_limiter_active_fraction_boundary[0] =
            (global_region_counts[1] > 0)
            ? static_cast<double>(global_region_counts[0]) /
              static_cast<double>(global_region_counts[1])
            : 0.0;
        fluxes.region_limiter_active_fraction_core[0] =
            (global_region_counts[3] > 0)
            ? static_cast<double>(global_region_counts[2]) /
              static_cast<double>(global_region_counts[3])
            : 0.0;
        fluxes.region_limiter_active_fraction_boundary[1] =
            (global_region_counts[5] > 0)
            ? static_cast<double>(global_region_counts[4]) /
              static_cast<double>(global_region_counts[5])
            : 0.0;
        fluxes.region_limiter_active_fraction_core[1] =
            (global_region_counts[7] > 0)
            ? static_cast<double>(global_region_counts[6]) /
              static_cast<double>(global_region_counts[7])
            : 0.0;
    }
    {
        long long local_x_counts[4] = {
            local_x_active_boundary,
            local_x_total_boundary,
            local_x_active_core,
            local_x_total_core
        };
        long long global_x_counts[4] = {0, 0, 0, 0};
        MPI_Allreduce(local_x_counts, global_x_counts, 4,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        double local_x_alpha[2] = {
            local_x_min_alpha_boundary,
            local_x_min_alpha_core
        };
        double global_x_alpha[2] = {1.0, 1.0};
        MPI_Allreduce(local_x_alpha, global_x_alpha, 2, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        const long long active_total =
            global_x_counts[0] + global_x_counts[2];
        const long long face_total =
            global_x_counts[1] + global_x_counts[3];
        fluxes.limiter_active_fraction =
            (face_total > 0)
            ? static_cast<double>(active_total) /
              static_cast<double>(face_total)
            : 0.0;
        fluxes.limiter_active_fraction_boundary =
            (global_x_counts[1] > 0)
            ? static_cast<double>(global_x_counts[0]) /
              static_cast<double>(global_x_counts[1])
            : 0.0;
        fluxes.limiter_active_fraction_core =
            (global_x_counts[3] > 0)
            ? static_cast<double>(global_x_counts[2]) /
              static_cast<double>(global_x_counts[3])
            : 0.0;
        fluxes.limiter_min_alpha_boundary = global_x_alpha[0];
        fluxes.limiter_min_alpha_core = global_x_alpha[1];
        fluxes.limiter_min_alpha =
            std::min(global_x_alpha[0], global_x_alpha[1]);
    }
    {
        long long global_u[2] = { local_u_face_active, local_u_face_total };
        MPI_Allreduce(MPI_IN_PLACE, global_u, 2, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        double global_u_alpha = local_u_alpha_min;
        MPI_Allreduce(MPI_IN_PLACE, &global_u_alpha, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        fluxes.u_force_alpha_min = global_u_alpha;
        if (global_u[1] > 0) {
            fluxes.u_force_alpha_active_frac =
                static_cast<double>(global_u[0]) /
                static_cast<double>(global_u[1]);
        }
    }
    {
        // MPI-reduce f-negativity stats
        double global_f_neg_min = local_f_neg_min;
        MPI_Allreduce(MPI_IN_PLACE, &global_f_neg_min, 1, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        fluxes.f_neg_min = global_f_neg_min;

        double global_neg_ratio_max = local_f_neg_ratio_max;
        MPI_Allreduce(MPI_IN_PLACE, &global_neg_ratio_max, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        fluxes.f_neg_ratio_max = global_neg_ratio_max;

        double local_fneg_mass = local_f_neg_mass;
        double global_fneg_mass = 0.0;
        MPI_Allreduce(&local_fneg_mass, &global_fneg_mass, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        fluxes.f_neg_mass_total = global_fneg_mass;

        long long local_fneg_count = local_f_neg_count;
        MPI_Allreduce(&local_fneg_count, &fluxes.f_neg_cell_count, 1,
                      MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        // Worst-negative location: MPI_MAXLOC on (neg_ratio, rank)
        struct { double val; int rank; } loc_worst, glob_worst;
        loc_worst.val = local_worst_neg.valid ? local_worst_neg.neg_ratio : -1.0;
        loc_worst.rank = mpi_rank;
        glob_worst.val = -1.0;
        glob_worst.rank = -1;
        MPI_Allreduce(&loc_worst, &glob_worst, 1, MPI_DOUBLE_INT,
                      MPI_MAXLOC, MPI_COMM_WORLD);
        int worst_indices[3] = {
            local_worst_neg.valid ? local_worst_neg.ix : -1,
            local_worst_neg.valid ? local_worst_neg.iv : -1,
            local_worst_neg.valid ? local_worst_neg.imu : -1
        };
        if (glob_worst.rank >= 0 && glob_worst.val > 0.0) {
            MPI_Bcast(worst_indices, 3, MPI_INT, glob_worst.rank,
                      MPI_COMM_WORLD);
            fluxes.f_neg_ix = worst_indices[0];
            fluxes.f_neg_iv = worst_indices[1];
            fluxes.f_neg_imu = worst_indices[2];
        }
    }
    std::fill(fluxes.cell_alpha_u.begin(), fluxes.cell_alpha_u.end(), 1.0);
    return;

}

void VlasovAmpereMidpointSolver::update_flux_current(
    const Species& sp,
    const SpatialGrid& sg,
    FluxPack& fluxes,
    Species& bkg_new) const
{
    const int nxl = sg.nx_local;
    resize_or_zero(bkg_new.current_face_x, static_cast<size_t>(nxl + 1));
    resize_or_zero(bkg_new.current_x, static_cast<size_t>(nxl));
    resize_or_zero(fluxes.j_bkg_face, static_cast<size_t>(nxl + 1));

    #pragma omp parallel for schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        double j_face = 0.0;
        const size_t face_base = static_cast<size_t>(iface) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double shell = sp.vgrid.moment_weight[iv];
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k =
                    static_cast<size_t>(iv) * Param::Nmu
                  + static_cast<size_t>(imu);
                j_face += sp.charge * shell *
                          fluxes.x_final[face_base + k];
            }
        }
        bkg_new.current_face_x[static_cast<size_t>(iface)] = j_face;
        fluxes.j_bkg_face[static_cast<size_t>(iface)] = j_face;
    }
    for (int ix = 0; ix < nxl; ++ix) {
        bkg_new.current_x[static_cast<size_t>(ix)] =
            0.5 * (bkg_new.current_face_x[static_cast<size_t>(ix)]
                 + bkg_new.current_face_x[static_cast<size_t>(ix + 1)]);
    }
}

double VlasovAmpereMidpointSolver::integrate_face_work(
    const std::vector<double>& current_face,
    const EMFields& fields_mid,
    const SpatialGrid& sg,
    double dt) const
{
    double work = 0.0;
    const size_t n = std::min(static_cast<size_t>(sg.nx_local),
                              std::min(current_face.size(),
                                       fields_mid.Ex_face.size()));
    for (size_t iface = 0; iface < n; ++iface) {
        work -= current_face[iface] * fields_mid.Ex_face[iface] * sg.dx * dt;
    }
    return work;
}

void VlasovAmpereMidpointSolver::build_current_diagnostics(
    const std::vector<double>& j_bkg_charge,
    const std::vector<double>& j_bkg_energy,
    const std::vector<double>& j_bkg_ampere,
    const EMFields& fields_mid,
    const SpatialGrid& sg,
    double local_delta_ke_bkg,
    double dt,
    CurrentDiagnostics& diag) const
{
    reset_current_diag(diag);
    const size_t n = std::min(static_cast<size_t>(sg.nx_local),
                              fields_mid.Ex_face.size());
    for (size_t iface = 0; iface < n; ++iface) {
        const double e = fields_mid.Ex_face[iface];
        const double j_charge =
            (iface < j_bkg_charge.size()) ? j_bkg_charge[iface] : 0.0;
        const double j_energy =
            (iface < j_bkg_energy.size()) ? j_bkg_energy[iface] : 0.0;
        const double j_ampere =
            (iface < j_bkg_ampere.size()) ? j_bkg_ampere[iface] : 0.0;
        diag.e_dot_j_charge += e * j_charge * sg.dx;
        diag.e_dot_j_energy += e * j_energy * sg.dx;
        diag.e_dot_j_ampere += e * j_ampere * sg.dx;
        diag.max_abs_j_charge =
            std::max(diag.max_abs_j_charge, std::fabs(j_charge));
        diag.max_abs_j_energy =
            std::max(diag.max_abs_j_energy, std::fabs(j_energy));
        diag.max_abs_j_ampere =
            std::max(diag.max_abs_j_ampere, std::fabs(j_ampere));
        diag.max_abs_j_charge_minus_ampere =
            std::max(diag.max_abs_j_charge_minus_ampere,
                     std::fabs(j_charge - j_ampere));
        diag.max_abs_j_energy_minus_ampere =
            std::max(diag.max_abs_j_energy_minus_ampere,
                     std::fabs(j_energy - j_ampere));
    }
    diag.residual_if_charge = local_delta_ke_bkg - dt * diag.e_dot_j_charge;
    diag.residual_if_ampere = local_delta_ke_bkg - dt * diag.e_dot_j_ampere;
}

bool VlasovAmpereMidpointSolver::check_finite_state(
    const Species& bkg,
    const BeamPIC& beam,
    const EMFields& fields,
    const std::vector<double>& ex_mid_next,
    const std::vector<double>& current_next) const
{
    for (size_t i = 0; i < bkg.f.size(); ++i) {
        if (!std::isfinite(bkg.f[i])) return false;
    }
    for (size_t i = 0; i < bkg.current_face_x.size(); ++i) {
        if (!std::isfinite(bkg.current_face_x[i])) return false;
    }
    for (size_t i = 0; i < beam.current_face_x.size(); ++i) {
        if (!std::isfinite(beam.current_face_x[i])) return false;
    }
    for (size_t i = 0; i < fields.Ex_face.size(); ++i) {
        if (!std::isfinite(fields.Ex_face[i])) return false;
    }
    for (size_t i = 0; i < ex_mid_next.size(); ++i) {
        if (!std::isfinite(ex_mid_next[i])) return false;
    }
    for (size_t i = 0; i < current_next.size(); ++i) {
        if (!std::isfinite(current_next[i])) return false;
    }
    return true;
}

VlasovAmpereMidpointSolver::Result
VlasovAmpereMidpointSolver::advance_single_step(
    const Species& bkg_n,
    const BeamPIC& beam_n,
    const EMFields& fields_n,
    const SpatialGrid& sg,
    double dt,
    double time,
    int mpi_rank,
    int mpi_size,
    int substeps_used) const
{
    Result result;
    reset_result(result);
    result.substeps_used = substeps_used;
    result.species_np1 = bkg_n;
    result.beam_np1 = beam_n;
    result.fields_np1 = fields_n;
    const int nxl = sg.nx_local;
    result.j_bkg_face_mid.assign(static_cast<size_t>(nxl + 1), 0.0);
    result.j_beam_face_mid.assign(static_cast<size_t>(nxl + 1), 0.0);
    result.j_total_face_mid.assign(static_cast<size_t>(nxl + 1), 0.0);
    result.j_bkg_energy_debug_face.assign(static_cast<size_t>(nxl + 1), 0.0);

    Species bkg_guess = bkg_n;
    std::vector<double> ex_mid_trial(
        fields_n.Ex_face.begin(),
        fields_n.Ex_face.begin() +
        std::min(fields_n.Ex_face.size(), static_cast<size_t>(nxl)));
    if (ex_mid_trial.size() != static_cast<size_t>(nxl)) {
        ex_mid_trial.assign(static_cast<size_t>(nxl), 0.0);
    }
    std::vector<double> ex_mid_next = ex_mid_trial;
    std::vector<double> j_total_prev(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_bkg_prev(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_beam_prev(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_total_next(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_bkg_next(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_beam_next(static_cast<size_t>(nxl), 0.0);
    std::vector<double> j_total_face(static_cast<size_t>(nxl + 1), 0.0);

    double local_n_start = 0.0;
    double local_ke_start = 0.0;
    bkg_n.total_particle_number_and_energy(local_n_start, local_ke_start);
    const int max_iters = 50;
    const double field_tol = 1.0e-6;
    const double current_tol = 1.0e-5;
    const double f_tol = 1.0e-4;
    const double omega_min = 0.05;
    double omega = 0.5;
    double previous_residual = -1.0;
    bool have_previous_current = false;

    auto compute_beam_trial =
        [&](const std::vector<double>& ex_mid,
            BeamPIC& beam_trial,
            double& beam_ke_before,
            double& beam_ke_after,
            double& beam_continuity,
            std::vector<double>& j_beam_face,
            std::vector<double>& j_beam_local) {
            EMFields fields_mid;
            set_midpoint_field(fields_mid, fields_n, ex_mid, sg,
                               mpi_rank, mpi_size);
            beam_trial = beam_n;
            beam_trial.begin_step(sg, dt);
            beam_trial.inject(sg, fields_mid, dt, time, mpi_rank, mpi_size);
            beam_ke_before = beam_trial.total_kinetic_energy();
            beam_trial.push(sg, fields_mid, dt, mpi_rank, mpi_size);
            beam_ke_after = beam_trial.total_kinetic_energy();
            beam_trial.deposit_density(sg, mpi_rank, mpi_size);
            beam_trial.finalize_charge_conserving_current(sg, dt, mpi_rank,
                                                          mpi_size);
            j_beam_face.assign(static_cast<size_t>(nxl + 1), 0.0);
            const size_t nface =
                std::min(j_beam_face.size(), beam_trial.current_face_x.size());
            std::copy(beam_trial.current_face_x.begin(),
                      beam_trial.current_face_x.begin() + nface,
                      j_beam_face.begin());
            j_beam_local.assign(static_cast<size_t>(nxl), 0.0);
            for (int iface = 0; iface < nxl; ++iface) {
                j_beam_local[static_cast<size_t>(iface)] =
                    j_beam_face[static_cast<size_t>(iface)];
            }
            beam_continuity =
                std::max(beam_trial.last_continuity_linf_error(),
                         beam_trial.last_boundary_flux_error());
        };

    BeamPIC beam_predictor;
    double beam_ke_before_predictor = 0.0;
    double beam_ke_after_predictor = 0.0;
    double beam_continuity_predictor = 0.0;
    std::vector<double> j_beam_predictor_face;
    std::vector<double> j_beam_predictor_local;
    compute_beam_trial(ex_mid_trial, beam_predictor,
                       beam_ke_before_predictor, beam_ke_after_predictor,
                       beam_continuity_predictor, j_beam_predictor_face,
                       j_beam_predictor_local);

    FluxPack fluxes;
    for (int iter = 0; iter < max_iters; ++iter) {
        EMFields fields_mid;
        set_midpoint_field(fields_mid, fields_n, ex_mid_trial, sg,
                           mpi_rank, mpi_size);

        Species bkg_new;
        bool finite_flux = true;
        compute_midpoint_fluxes(bkg_n, bkg_guess, fields_mid, sg, dt,
                                mpi_rank, mpi_size, bkg_new, fluxes, 0, 0.0,
                                finite_flux);
        bkg_new.compute_moments();
        sync_cell_current_from_faces(bkg_new, nxl);

        for (int iface = 0; iface <= nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            const double jb = (slot < bkg_new.current_face_x.size())
                            ? bkg_new.current_face_x[slot] : 0.0;
            const double jbeam = (slot < j_beam_predictor_face.size())
                               ? j_beam_predictor_face[slot] : 0.0;
            j_total_face[slot] = jb + jbeam;
        }

        EMFields fields_new = fields_n;
        fields_new.advance_ampere_face_from_midpoint_current(
            j_total_face, dt, mpi_rank, mpi_size);

        for (int iface = 0; iface < nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            j_bkg_next[slot] =
                (slot < bkg_new.current_face_x.size())
                ? bkg_new.current_face_x[slot] : 0.0;
            j_beam_next[slot] =
                (slot < j_beam_predictor_local.size())
                ? j_beam_predictor_local[slot] : 0.0;
            j_total_next[slot] = j_bkg_next[slot] + j_beam_next[slot];
            ex_mid_next[slot] =
                0.5 * (fields_n.Ex_face[slot] + fields_new.Ex_face[slot]);
        }

        double local_delta_ex = 0.0;
        double local_abs_ex = 0.0;
        double local_delta_f = 0.0;
        double local_abs_f = 0.0;
        #pragma omp parallel for schedule(static) \
            reduction(max:local_delta_f,local_abs_f)
        for (long long i = 0; i < static_cast<long long>(bkg_new.f.size());
             ++i) {
            const size_t slot = static_cast<size_t>(i);
            local_delta_f =
                std::max(local_delta_f,
                         std::fabs(bkg_new.f[slot] - bkg_guess.f[slot]));
            local_abs_f =
                std::max(local_abs_f,
                         std::max(std::fabs(bkg_new.f[slot]),
                                  std::fabs(bkg_guess.f[slot])));
        }
        for (int iface = 0; iface < nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            local_delta_ex =
                std::max(local_delta_ex,
                         std::fabs(ex_mid_next[slot] - ex_mid_trial[slot]));
            local_abs_ex =
                std::max(local_abs_ex,
                         std::max(std::fabs(ex_mid_next[slot]),
                                  std::fabs(ex_mid_trial[slot])));
        }

        double local_delta_j_total = 0.0;
        double local_delta_j_bkg = 0.0;
        double local_delta_j_beam = 0.0;
        if (have_previous_current) {
            for (int iface = 0; iface < nxl; ++iface) {
                const size_t slot = static_cast<size_t>(iface);
                local_delta_j_total =
                    std::max(local_delta_j_total,
                             std::fabs(j_total_next[slot]
                                     - j_total_prev[slot]));
                local_delta_j_bkg =
                    std::max(local_delta_j_bkg,
                             std::fabs(j_bkg_next[slot]
                                     - j_bkg_prev[slot]));
                local_delta_j_beam =
                    std::max(local_delta_j_beam,
                             std::fabs(j_beam_next[slot]
                                     - j_beam_prev[slot]));
            }
        }

        double local_end_values[2] = {0.0, 0.0};
        bkg_new.total_particle_number_and_energy(local_end_values[0],
                                                 local_end_values[1]);
        double global_mass_delta = local_end_values[0] - local_n_start;
        MPI_Allreduce(MPI_IN_PLACE, &global_mass_delta, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        const double local_beam_continuity = beam_continuity_predictor;
        const double local_max_j_total =
            max_abs_vector(j_total_next, static_cast<size_t>(nxl));
        const double local_max_j_bkg =
            max_abs_vector(j_bkg_next, static_cast<size_t>(nxl));
        const double local_max_j_beam =
            max_abs_vector(j_beam_next, static_cast<size_t>(nxl));
        const double local_errors[13] = {
            local_delta_ex,
            local_abs_ex,
            local_delta_f,
            local_abs_f,
            local_delta_j_total,
            local_delta_j_bkg,
            local_delta_j_beam,
            local_max_j_total,
            local_max_j_bkg,
            local_max_j_beam,
            std::fabs(global_mass_delta),
            local_beam_continuity,
            finite_flux ? 0.0 : 1.0
        };
        double global_errors[13] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        MPI_Allreduce(local_errors, global_errors, 13, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);

        const double global_field_error =
            global_errors[0] / std::max(1.0, global_errors[1]);
        const double global_f_error =
            global_errors[2] / std::max(1.0, global_errors[3]);
        const double global_current_scale =
            std::max(std::fabs(Param::jb), global_errors[7]);
        const double global_j_abs_tol =
            current_tol * std::max(1.0, global_current_scale);
        const double global_bkg_mass_tol =
            std::max(1.0, 1.0e-10 * Param::dens * Param::Lx);
        const double global_beam_cont_tol = 1.0e-6;
        const double normalized_residual =
            std::max(std::max(global_field_error / field_tol,
                              global_f_error / f_tol),
                     std::max(global_errors[4] /
                              std::max(1.0, global_j_abs_tol),
                              std::max(global_errors[10] /
                                       global_bkg_mass_tol,
                                       global_errors[11] /
                                       global_beam_cont_tol)));

        result.nonlinear_residual = normalized_residual;
        result.residual_E = global_field_error;
        result.residual_f = global_f_error;
        result.residual_J_bkg =
            global_errors[5] /
            std::max(1.0, std::max(std::fabs(Param::jb),
                                   global_errors[8]));
        result.residual_J_beam =
            global_errors[6] /
            std::max(1.0, std::max(std::fabs(Param::jb),
                                   global_errors[9]));
        result.continuity_residual_bkg =
            global_errors[10] / global_bkg_mass_tol;
        result.beam_continuity_residual = global_errors[11];
        result.nonlinear_iterations = iter + 1;
        result.limiter_active_fraction = fluxes.limiter_active_fraction;
        result.limiter_min_alpha = fluxes.limiter_min_alpha;
        result.limiter_active_fraction_core =
            fluxes.limiter_active_fraction_core;
        result.limiter_active_fraction_boundary =
            fluxes.limiter_active_fraction_boundary;
        result.limiter_min_alpha_core = fluxes.limiter_min_alpha_core;
        result.limiter_min_alpha_boundary =
            fluxes.limiter_min_alpha_boundary;
        result.limiter_energy_defect = fluxes.limiter_energy_defect;
        result.limiter_mass_defect = fluxes.limiter_mass_defect;
        result.limiter_momentum_defect =
            fluxes.limiter_momentum_defect;
        result.x_limiter_mass_defect = fluxes.x_limiter_mass_defect;
        result.x_limiter_energy_defect = fluxes.x_limiter_energy_defect;
        result.mu_low_u_alpha_min = fluxes.mu_low_u_alpha_min;
        result.mu_low_u_limiter_active_fraction =
            fluxes.mu_low_u_limiter_active_fraction;
        result.mu_low_u_energy_delta = fluxes.mu_low_u_energy_delta;
        result.mu_low_u_alpha_min_boundary =
            fluxes.mu_low_u_alpha_min_boundary;
        result.mu_low_u_alpha_min_core = fluxes.mu_low_u_alpha_min_core;
        result.mu_low_u_limiter_active_fraction_boundary =
            fluxes.mu_low_u_limiter_active_fraction_boundary;
        result.mu_low_u_limiter_active_fraction_core =
            fluxes.mu_low_u_limiter_active_fraction_core;
        result.mu_low_u_energy_delta_boundary =
            fluxes.mu_low_u_energy_delta_boundary;
        result.mu_low_u_energy_delta_core =
            fluxes.mu_low_u_energy_delta_core;
        result.mu_low_u_u_eff0 = fluxes.mu_low_u_u_eff0;
        result.mu_low_u_moment_weight0 = fluxes.mu_low_u_moment_weight0;
        result.mu_low_u_mu_flux_scale0 = fluxes.mu_low_u_mu_flux_scale0;
        result.mu_low_u_half_dt_inv_shell0 =
            fluxes.mu_low_u_half_dt_inv_shell0;
        result.mu_low_u_dimless_scale0 = fluxes.mu_low_u_dimless_scale0;
        result.mu_low_u_endpoint_flux_max =
            fluxes.mu_low_u_endpoint_flux_max;
        result.remap_active_fraction = fluxes.remap_active_fraction;
        result.remap_cell_count = fluxes.remap_cell_count;
        result.low_u_subcycle_active_fraction =
            fluxes.low_u_subcycle_active_fraction;
        result.low_u_average_subcycles = fluxes.low_u_average_subcycles;
        result.low_u_max_subcycles = fluxes.low_u_max_subcycles;
        for (int ir = 0; ir < 2; ++ir) {
            result.region_u_limiter_energy_boundary[ir] =
                fluxes.region_u_limiter_energy_boundary[ir];
            result.region_u_limiter_energy_core[ir] =
                fluxes.region_u_limiter_energy_core[ir];
            result.region_abs_u_limiter_energy_boundary[ir] =
                fluxes.region_abs_u_limiter_energy_boundary[ir];
            result.region_abs_u_limiter_energy_core[ir] =
                fluxes.region_abs_u_limiter_energy_core[ir];
            result.region_limiter_active_fraction_boundary[ir] =
                fluxes.region_limiter_active_fraction_boundary[ir];
            result.region_limiter_active_fraction_core[ir] =
                fluxes.region_limiter_active_fraction_core[ir];
        }
        result.stage_min_f = fluxes.stage_min_f;
        result.stage_neg_mass = fluxes.stage_neg_mass;
        result.stage_neg_cell_count = fluxes.stage_neg_cell_count;
        result.stage_low_u_neg_mass = fluxes.stage_low_u_neg_mass;
        result.low_u_neg_added_by_div = fluxes.low_u_neg_added_by_div;
        result.stage_core_low_u_min_f = fluxes.stage_core_low_u_min_f;
        result.stage_min_f_core_by_u = fluxes.stage_min_f_core_by_u;
        result.stage_neg_mass_core_by_u = fluxes.stage_neg_mass_core_by_u;
        result.stage_neg_cell_count_core_by_u =
            fluxes.stage_neg_cell_count_core_by_u;
        result.stage_min_f_boundary_by_u =
            fluxes.stage_min_f_boundary_by_u;
        result.stage_neg_mass_boundary_by_u =
            fluxes.stage_neg_mass_boundary_by_u;
        result.stage_neg_cell_count_boundary_by_u =
            fluxes.stage_neg_cell_count_boundary_by_u;
        result.x_negative_mass_before_repair =
            fluxes.negative_mass_before_repair;
        result.x_mass_added_by_positivity_repair =
            fluxes.mass_added_by_positivity_repair;
        // 7.1.2: x low-order flux failure count
        result.x_low_order_failed_count =
            fluxes.x_low_order_failed_count;
        result.x_low_input_min_f = fluxes.x_low_input_min_f;
        result.x_low_max_cfl = fluxes.x_low_max_cfl;
        result.x_low_output_min_f = fluxes.x_low_output_min_f;
        result.x_low_failed_count = fluxes.x_low_failed_count;
        result.x_low_input_neg_mass = fluxes.x_low_input_neg_mass;
        result.x_low_input_rel_neg = fluxes.x_low_input_rel_neg;
        result.x_low_output_rel_neg = fluxes.x_low_output_rel_neg;
        result.x_low_input_core_failed_count =
            fluxes.x_low_input_core_failed_count;
        result.x_low_input_debt_accepted =
            fluxes.x_low_input_debt_accepted;
        result.x_low_failure_kind = fluxes.x_low_failure_kind;
        // 7.1.6: copy per-direction flux diagnostics
        for (int d = 0; d < 3; ++d) {
            result.flux_pos[d] = fluxes.flux_pos[d];
            result.flux_defect[d] = fluxes.flux_defect[d];
        }
        result.positivity_energy_defect =
            fluxes.positivity_energy_defect;
        result.positivity_mass_defect =
            fluxes.positivity_mass_defect;
        result.u_force_alpha_min = fluxes.u_force_alpha_min;
        result.u_force_alpha_active_frac =
            fluxes.u_force_alpha_active_frac;
        result.u_flux_audit_valid = fluxes.u_flux_audit_valid;
        result.u_flux_audit_rank = fluxes.u_flux_audit_rank;
        result.u_flux_audit_ix = fluxes.u_flux_audit_ix;
        result.u_flux_audit_iv = fluxes.u_flux_audit_iv;
        result.u_flux_audit_imu = fluxes.u_flux_audit_imu;
        result.u_flux_audit_severity = fluxes.u_flux_audit_severity;
        result.u_flux_audit_f0 = fluxes.u_flux_audit_f0;
        result.u_flux_audit_f_low = fluxes.u_flux_audit_f_low;
        result.u_flux_audit_f_high = fluxes.u_flux_audit_f_high;
        result.u_flux_audit_alpha = fluxes.u_flux_audit_alpha;
        result.u_flux_audit_du_div_low = fluxes.u_flux_audit_du_div_low;
        result.u_flux_audit_du_div_high = fluxes.u_flux_audit_du_div_high;
        result.u_flux_audit_du_div_final =
            fluxes.u_flux_audit_du_div_final;
        result.u_flux_audit_updated = fluxes.u_flux_audit_updated;
        result.u_flux_audit_final_xl_lo =
            fluxes.u_flux_audit_final_xl_lo;
        result.u_flux_audit_final_xl_hi =
            fluxes.u_flux_audit_final_xl_hi;
        result.u_flux_audit_final_xr_lo =
            fluxes.u_flux_audit_final_xr_lo;
        result.u_flux_audit_final_xr_hi =
            fluxes.u_flux_audit_final_xr_hi;
        result.u_low_failure_audit = fluxes.u_low_failure_audit;
        result.mu_low_failure_audit = fluxes.mu_low_failure_audit;
        result.f_neg_min = fluxes.f_neg_min;
        result.f_neg_ratio_max = fluxes.f_neg_ratio_max;
        result.f_neg_mass_total = fluxes.f_neg_mass_total;
        result.f_neg_cell_count = fluxes.f_neg_cell_count;
        result.f_neg_ix = fluxes.f_neg_ix;
        result.f_neg_iv = fluxes.f_neg_iv;
        result.f_neg_imu = fluxes.f_neg_imu;

        const bool finite_state =
            global_errors[12] == 0.0 &&
            check_finite_state(bkg_new, beam_predictor, fields_new,
                               ex_mid_next, j_total_next);
        if (!finite_state) {
            if (mpi_rank == 0 && global_errors[12] != 0.0) {
                if (fluxes.x_low_failure_kind != X_LOW_OK) {
                    const char* kind = "unknown";
                    const char* action = "not retrying as a CFL failure";
                    if (fluxes.x_low_failure_kind == X_LOW_INPUT_BAD) {
                        kind = "input_state_bad";
                        action = "inspect the previous accepted u/mu/x stage";
                    } else if (fluxes.x_low_failure_kind == X_LOW_INPUT_DEBT) {
                        kind = "input_state_debt";
                        action = "accepted as a small residual negative debt";
                    } else if (fluxes.x_low_failure_kind == X_LOW_DONOR_BUG) {
                        kind = "donor_or_periodic_face_bug";
                        action = "fix x_low donor/periodic face construction";
                    } else if (fluxes.x_low_failure_kind == X_LOW_TRUE_CFL) {
                        kind = "true_cfl";
                        action = "retrying with a smaller true x substep";
                    }
                    std::fprintf(
                        stderr,
                        "ERROR: x low-order positivity failure classified "
                        "as %s: failed_cell_count=%.0f input_min_f=%.16e "
                        "max_cfl=%.16e output_min_f=%.16e "
                        "input_neg_mass=%.16e input_rel_neg=%.16e "
                        "output_rel_neg=%.16e input_core_failed_count=%.0f "
                        "debt_accepted=%.0f; %s.\n",
                        kind,
                        fluxes.x_low_failed_count,
                        fluxes.x_low_input_min_f,
                        fluxes.x_low_max_cfl,
                        fluxes.x_low_output_min_f,
                        fluxes.x_low_input_neg_mass,
                        fluxes.x_low_input_rel_neg,
                        fluxes.x_low_output_rel_neg,
                        fluxes.x_low_input_core_failed_count,
                        fluxes.x_low_input_debt_accepted,
                        action);
                } else if (fluxes.finite_flux_has_failure != 0) {
                    std::fprintf(
                        stderr,
                        "ERROR: finite_flux failure detail: "
                        "rank=%d ix=%d iv=%d imu=%d "
                        "max_negative=%.16e relative_negative=%.16e "
                        "updated=%.16e f0=%.16e "
                        "dx_div=%.16e du_div=%.16e dmu_div=%.16e\n",
                        fluxes.finite_flux_rank,
                        fluxes.finite_flux_ix,
                        fluxes.finite_flux_iv,
                        fluxes.finite_flux_imu,
                        fluxes.finite_flux_max_negative,
                        fluxes.finite_flux_relative_negative,
                        fluxes.finite_flux_updated,
                        fluxes.finite_flux_f0,
                        fluxes.finite_flux_dx_div,
                        fluxes.finite_flux_du_div,
                        fluxes.finite_flux_dmu_div);
                } else {
                    std::fprintf(
                        stderr,
                        "ERROR: finite_flux failure detail unavailable; "
                        "no negative/non-finite centered update was recorded\n");
                }
            }
            result.failed = true;
            result.converged = false;
            break;
        }

        const bool protected_quantities_converged =
            have_previous_current &&
            global_field_error < field_tol &&
            global_errors[4] < global_j_abs_tol &&
            global_errors[5] < global_j_abs_tol &&
            global_errors[6] < global_j_abs_tol &&
            global_errors[10] < global_bkg_mass_tol &&
            global_errors[11] < global_beam_cont_tol;
        result.protected_converged = protected_quantities_converged;
        const bool converged = protected_quantities_converged;
        const bool soft_accepted =
            protected_quantities_converged && global_f_error >= f_tol;

        if (converged) {
            BeamPIC beam_corrected;
            double beam_ke_before_corrected = 0.0;
            double beam_ke_after_corrected = 0.0;
            double beam_continuity_corrected = 0.0;
            std::vector<double> j_beam_corrected_face;
            std::vector<double> j_beam_corrected_local;
            compute_beam_trial(ex_mid_next, beam_corrected,
                               beam_ke_before_corrected,
                               beam_ke_after_corrected,
                               beam_continuity_corrected,
                               j_beam_corrected_face,
                               j_beam_corrected_local);
            MPI_Allreduce(MPI_IN_PLACE, &beam_continuity_corrected, 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            for (int iface = 0; iface <= nxl; ++iface) {
                const size_t slot = static_cast<size_t>(iface);
                const double jb =
                    (slot < bkg_new.current_face_x.size())
                    ? bkg_new.current_face_x[slot] : 0.0;
                const double jbeam =
                    (slot < j_beam_corrected_face.size())
                    ? j_beam_corrected_face[slot] : 0.0;
                j_total_face[slot] = jb + jbeam;
            }
            fields_new = fields_n;
            fields_new.advance_ampere_face_from_midpoint_current(
                j_total_face, dt, mpi_rank, mpi_size);
            EMFields fields_mid_corrected;
            set_midpoint_field(fields_mid_corrected, fields_n, ex_mid_next,
                               sg, mpi_rank, mpi_size);

            result.species_np1 = bkg_new;
            result.beam_np1 = beam_corrected;
            result.fields_np1 = fields_new;
            result.j_bkg_face_mid = bkg_new.current_face_x;
            result.j_beam_face_mid = j_beam_corrected_face;
            result.j_total_face_mid = j_total_face;
            result.j_bkg_energy_debug_face = bkg_new.current_face_x;
            result.delta_ke_bkg = local_end_values[1] - local_ke_start;
            result.delta_ke_beam =
                beam_ke_after_corrected - beam_ke_before_corrected;
            result.field_work_bkg =
                integrate_face_work(bkg_new.current_face_x,
                                    fields_mid_corrected, sg, dt);
            result.field_work_beam =
                integrate_face_work(j_beam_corrected_face,
                                    fields_mid_corrected, sg, dt);
            result.energy_residual_bkg =
                result.delta_ke_bkg + result.field_work_bkg;
            result.beam_continuity_residual = beam_continuity_corrected;
            build_current_diagnostics(result.j_bkg_face_mid,
                                      result.j_bkg_energy_debug_face,
                                      result.j_bkg_face_mid,
                                      fields_mid_corrected, sg,
                                      result.delta_ke_bkg, dt,
                                      result.current_diag);
            result.limiter_active_fraction =
                fluxes.limiter_active_fraction;
            result.limiter_min_alpha = fluxes.limiter_min_alpha;
            result.limiter_active_fraction_core =
                fluxes.limiter_active_fraction_core;
            result.limiter_active_fraction_boundary =
                fluxes.limiter_active_fraction_boundary;
            result.limiter_min_alpha_core =
                fluxes.limiter_min_alpha_core;
            result.limiter_min_alpha_boundary =
                fluxes.limiter_min_alpha_boundary;
            result.limiter_energy_defect = fluxes.limiter_energy_defect;
            result.limiter_mass_defect = fluxes.limiter_mass_defect;
            result.limiter_momentum_defect =
                fluxes.limiter_momentum_defect;
            result.x_limiter_mass_defect = fluxes.x_limiter_mass_defect;
            result.x_limiter_energy_defect =
                fluxes.x_limiter_energy_defect;
            result.mu_low_u_alpha_min = fluxes.mu_low_u_alpha_min;
            result.mu_low_u_limiter_active_fraction =
                fluxes.mu_low_u_limiter_active_fraction;
            result.mu_low_u_energy_delta = fluxes.mu_low_u_energy_delta;
            result.mu_low_u_alpha_min_boundary =
                fluxes.mu_low_u_alpha_min_boundary;
            result.mu_low_u_alpha_min_core =
                fluxes.mu_low_u_alpha_min_core;
            result.mu_low_u_limiter_active_fraction_boundary =
                fluxes.mu_low_u_limiter_active_fraction_boundary;
            result.mu_low_u_limiter_active_fraction_core =
                fluxes.mu_low_u_limiter_active_fraction_core;
            result.mu_low_u_energy_delta_boundary =
                fluxes.mu_low_u_energy_delta_boundary;
            result.mu_low_u_energy_delta_core =
                fluxes.mu_low_u_energy_delta_core;
            result.mu_low_u_u_eff0 = fluxes.mu_low_u_u_eff0;
            result.mu_low_u_moment_weight0 =
                fluxes.mu_low_u_moment_weight0;
            result.mu_low_u_mu_flux_scale0 =
                fluxes.mu_low_u_mu_flux_scale0;
            result.mu_low_u_half_dt_inv_shell0 =
                fluxes.mu_low_u_half_dt_inv_shell0;
            result.mu_low_u_dimless_scale0 =
                fluxes.mu_low_u_dimless_scale0;
            result.mu_low_u_endpoint_flux_max =
                fluxes.mu_low_u_endpoint_flux_max;
            result.remap_active_fraction = fluxes.remap_active_fraction;
            result.remap_cell_count = fluxes.remap_cell_count;
            result.low_u_subcycle_active_fraction =
                fluxes.low_u_subcycle_active_fraction;
            result.low_u_average_subcycles =
                fluxes.low_u_average_subcycles;
            result.low_u_max_subcycles = fluxes.low_u_max_subcycles;
            for (int ir = 0; ir < 2; ++ir) {
                result.region_u_limiter_energy_boundary[ir] =
                    fluxes.region_u_limiter_energy_boundary[ir];
                result.region_u_limiter_energy_core[ir] =
                    fluxes.region_u_limiter_energy_core[ir];
                result.region_abs_u_limiter_energy_boundary[ir] =
                    fluxes.region_abs_u_limiter_energy_boundary[ir];
                result.region_abs_u_limiter_energy_core[ir] =
                    fluxes.region_abs_u_limiter_energy_core[ir];
                result.region_limiter_active_fraction_boundary[ir] =
                    fluxes.region_limiter_active_fraction_boundary[ir];
                result.region_limiter_active_fraction_core[ir] =
                    fluxes.region_limiter_active_fraction_core[ir];
            }
            result.stage_min_f = fluxes.stage_min_f;
            result.stage_neg_mass = fluxes.stage_neg_mass;
            result.stage_neg_cell_count = fluxes.stage_neg_cell_count;
            result.stage_low_u_neg_mass = fluxes.stage_low_u_neg_mass;
            result.low_u_neg_added_by_div =
                fluxes.low_u_neg_added_by_div;
            result.stage_core_low_u_min_f = fluxes.stage_core_low_u_min_f;
            result.stage_min_f_core_by_u = fluxes.stage_min_f_core_by_u;
            result.stage_neg_mass_core_by_u =
                fluxes.stage_neg_mass_core_by_u;
            result.stage_neg_cell_count_core_by_u =
                fluxes.stage_neg_cell_count_core_by_u;
            result.stage_min_f_boundary_by_u =
                fluxes.stage_min_f_boundary_by_u;
            result.stage_neg_mass_boundary_by_u =
                fluxes.stage_neg_mass_boundary_by_u;
            result.stage_neg_cell_count_boundary_by_u =
                fluxes.stage_neg_cell_count_boundary_by_u;
            result.x_negative_mass_before_repair =
                fluxes.negative_mass_before_repair;
            result.x_mass_added_by_positivity_repair =
                fluxes.mass_added_by_positivity_repair;
            // 7.1.2: x low-order flux failure count
            result.x_low_order_failed_count =
                fluxes.x_low_order_failed_count;
            result.x_low_input_min_f = fluxes.x_low_input_min_f;
            result.x_low_max_cfl = fluxes.x_low_max_cfl;
            result.x_low_output_min_f = fluxes.x_low_output_min_f;
            result.x_low_failed_count = fluxes.x_low_failed_count;
            result.x_low_input_neg_mass = fluxes.x_low_input_neg_mass;
            result.x_low_input_rel_neg = fluxes.x_low_input_rel_neg;
            result.x_low_output_rel_neg = fluxes.x_low_output_rel_neg;
            result.x_low_input_core_failed_count =
                fluxes.x_low_input_core_failed_count;
            result.x_low_input_debt_accepted =
                fluxes.x_low_input_debt_accepted;
            result.x_low_failure_kind = fluxes.x_low_failure_kind;
            // 7.1.6: copy per-direction flux diagnostics
            for (int d = 0; d < 3; ++d) {
                result.flux_pos[d] = fluxes.flux_pos[d];
                result.flux_defect[d] = fluxes.flux_defect[d];
            }
            result.positivity_energy_defect =
                fluxes.positivity_energy_defect;
            result.positivity_mass_defect =
                fluxes.positivity_mass_defect;
            result.u_force_alpha_min = fluxes.u_force_alpha_min;
            result.u_force_alpha_active_frac =
                fluxes.u_force_alpha_active_frac;
            result.u_flux_audit_valid = fluxes.u_flux_audit_valid;
            result.u_flux_audit_rank = fluxes.u_flux_audit_rank;
            result.u_flux_audit_ix = fluxes.u_flux_audit_ix;
            result.u_flux_audit_iv = fluxes.u_flux_audit_iv;
            result.u_flux_audit_imu = fluxes.u_flux_audit_imu;
            result.u_flux_audit_severity = fluxes.u_flux_audit_severity;
            result.u_flux_audit_f0 = fluxes.u_flux_audit_f0;
            result.u_flux_audit_f_low = fluxes.u_flux_audit_f_low;
            result.u_flux_audit_f_high = fluxes.u_flux_audit_f_high;
            result.u_flux_audit_alpha = fluxes.u_flux_audit_alpha;
            result.u_flux_audit_du_div_low =
                fluxes.u_flux_audit_du_div_low;
            result.u_flux_audit_du_div_high =
                fluxes.u_flux_audit_du_div_high;
            result.u_flux_audit_du_div_final =
                fluxes.u_flux_audit_du_div_final;
            result.u_flux_audit_updated = fluxes.u_flux_audit_updated;
            result.u_flux_audit_final_xl_lo =
                fluxes.u_flux_audit_final_xl_lo;
            result.u_flux_audit_final_xl_hi =
                fluxes.u_flux_audit_final_xl_hi;
            result.u_flux_audit_final_xr_lo =
                fluxes.u_flux_audit_final_xr_lo;
            result.u_flux_audit_final_xr_hi =
                fluxes.u_flux_audit_final_xr_hi;
            result.u_low_failure_audit = fluxes.u_low_failure_audit;
            result.mu_low_failure_audit = fluxes.mu_low_failure_audit;
            result.f_neg_min = fluxes.f_neg_min;
            result.f_neg_ratio_max = fluxes.f_neg_ratio_max;
            result.f_neg_mass_total = fluxes.f_neg_mass_total;
            result.f_neg_cell_count = fluxes.f_neg_cell_count;
            result.f_neg_ix = fluxes.f_neg_ix;
            result.f_neg_iv = fluxes.f_neg_iv;
            result.f_neg_imu = fluxes.f_neg_imu;
            result.converged = true;
            result.soft_accepted = soft_accepted;
            break;
        }

        if (previous_residual > 0.0) {
            if (normalized_residual < 0.9 * previous_residual) {
                omega = std::min(1.0, omega + 0.1);
            } else if (normalized_residual > 1.02 * previous_residual) {
                const double ratio = normalized_residual /
                    std::max(previous_residual, 1.0e-300);
                if (ratio > 4.0) {
                    omega = 0.0625;
                } else if (ratio > 2.0) {
                    omega = 0.125;
                } else if (ratio > 1.25) {
                    omega = 0.25;
                } else {
                    omega = 0.5;
                }
                omega = std::max(omega_min, std::min(omega, 0.5));
            } else {
                omega = std::max(omega_min, 0.8 * omega);
            }
        }
        previous_residual = normalized_residual;

        for (int iface = 0; iface < nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            ex_mid_trial[slot] =
                (1.0 - omega) * ex_mid_trial[slot]
              + omega * ex_mid_next[slot];
        }
        for (size_t i = 0; i < bkg_guess.f.size() && i < bkg_new.f.size();
             ++i) {
            bkg_guess.f[i] =
                (1.0 - omega) * bkg_guess.f[i] + omega * bkg_new.f[i];
        }
        exchange_ghosts_x_persistent(bkg_guess, sg, mpi_rank, mpi_size);
        j_total_prev.swap(j_total_next);
        j_bkg_prev.swap(j_bkg_next);
        j_beam_prev.swap(j_beam_next);
        have_previous_current = true;
    }

    if (!result.converged && !result.failed) {
        result.failed = true;
    }
    return result;
}
