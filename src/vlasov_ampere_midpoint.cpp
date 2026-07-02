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
    Result best;
    reset_result(best);
    best.failed = true;
    best.nonlinear_residual = 1.0e300;

    for (int trial = 0; trial < 7; ++trial) {
        const int substeps = substep_trials[trial];
        Result result =
            (substeps == 1)
            ? advance_single_step(bkg_n, beam_n, fields_n, sg, dt, time,
                                  mpi_rank, mpi_size, 1)
            : advance_with_fixed_substeps(bkg_n, beam_n, fields_n, sg, dt,
                                          time, mpi_rank, mpi_size,
                                          substeps);
        if (result.nonlinear_residual < best.nonlinear_residual ||
            (best.failed && !result.failed)) {
            best = result;
        }
        if (result.converged && !result.failed) {
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

    return best;
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
    for (int ir = 0; ir < 2; ++ir) {
        result.region_u_limiter_energy_boundary[ir] = 0.0;
        result.region_u_limiter_energy_core[ir] = 0.0;
        result.region_abs_u_limiter_energy_boundary[ir] = 0.0;
        result.region_abs_u_limiter_energy_core[ir] = 0.0;
        result.region_limiter_active_fraction_boundary[ir] = 0.0;
        result.region_limiter_active_fraction_core[ir] = 0.0;
    }
    result.x_negative_mass_before_repair = 0.0;
    result.x_mass_added_by_positivity_repair = 0.0;
    result.positivity_energy_defect = 0.0;
    result.positivity_mass_defect = 0.0;
    result.u_force_alpha_min = 1.0;
    result.u_force_alpha_active_frac = 0.0;
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
    resize_or_zero(fluxes.cell_alpha, cells);
    resize_or_zero(fluxes.j_bkg_face, static_cast<size_t>(nxl + 1));
    fluxes.limiter_active_fraction = 0.0;
    fluxes.limiter_min_alpha = 1.0;
    fluxes.limiter_active_fraction_core = 0.0;
    fluxes.limiter_active_fraction_boundary = 0.0;
    fluxes.limiter_min_alpha_core = 1.0;
    fluxes.limiter_min_alpha_boundary = 1.0;
    fluxes.limiter_energy_defect = 0.0;
    fluxes.limiter_mass_defect = 0.0;
    fluxes.limiter_momentum_defect = 0.0;
    for (int ir = 0; ir < 2; ++ir) {
        fluxes.region_u_limiter_energy_boundary[ir] = 0.0;
        fluxes.region_u_limiter_energy_core[ir] = 0.0;
        fluxes.region_abs_u_limiter_energy_boundary[ir] = 0.0;
        fluxes.region_abs_u_limiter_energy_core[ir] = 0.0;
        fluxes.region_limiter_active_fraction_boundary[ir] = 0.0;
        fluxes.region_limiter_active_fraction_core[ir] = 0.0;
    }
    fluxes.negative_mass_before_repair = 0.0;
    fluxes.mass_added_by_positivity_repair = 0.0;
    fluxes.positivity_energy_defect = 0.0;
    fluxes.positivity_mass_defect = 0.0;
    fluxes.u_force_alpha_min = 1.0;
    fluxes.u_force_alpha_active_frac = 0.0;
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

    #pragma omp parallel for schedule(static)
    for (int iface = 0; iface <= nxl; ++iface) {
        const double ex_face =
            (static_cast<size_t>(iface) < fields_mid.Ex_face.size())
            ? fields_mid.Ex_face[static_cast<size_t>(iface)] : 0.0;
        const double accel_u =
            bkg_n.charge * ex_face / (bkg_n.mass * Const::c);
        const int ix_left = ng + iface - 1;
        const int ix_right = ng + iface;

        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double u_dot_base = accel_u * bkg_n.vgrid.mu_cells[imu];
            for (int face = 1; face < Param::Nv; ++face) {
                const size_t left_lo =
                    static_cast<size_t>(ix_left) * Param::Nvmu
                  + static_cast<size_t>(face - 1) * Param::Nmu
                  + static_cast<size_t>(imu);
                const size_t left_hi =
                    static_cast<size_t>(ix_left) * Param::Nvmu
                  + static_cast<size_t>(face) * Param::Nmu
                  + static_cast<size_t>(imu);
                const size_t right_lo =
                    static_cast<size_t>(ix_right) * Param::Nvmu
                  + static_cast<size_t>(face - 1) * Param::Nmu
                  + static_cast<size_t>(imu);
                const size_t right_hi =
                    static_cast<size_t>(ix_right) * Param::Nvmu
                  + static_cast<size_t>(face) * Param::Nmu
                  + static_cast<size_t>(imu);
                const double f_lo =
                    0.5 * (f_mid.f[left_lo] + f_mid.f[right_lo]);
                const double f_hi =
                    0.5 * (f_mid.f[left_hi] + f_mid.f[right_hi]);
                const double scale =
                    2.0 * Const::pi * bkg_n.vgrid.dmu
                  * bkg_n.vgrid.v2_faces[face] * u_dot_base;
                u_force_face[u_xface_index(iface, face, imu)] =
                    scale * 0.5 * (f_lo + f_hi);
            }
        }

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u_eff =
                std::max(bkg_n.vgrid.v_cells[iv], Param::u_floor);
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
    close_periodic_face_blocks(u_force_face, nxl,
                               (Param::Nv + 1) * Param::Nmu,
                               mpi_rank, mpi_size, 506);
    close_periodic_face_blocks(mu_force_face, nxl,
                               Param::Nv * (Param::Nmu + 1),
                               mpi_rank, mpi_size, 507);

    /*
     * Centered midpoint x-flux (no FCT limiter — the flux is the single source
     * for charge-conserving current and velocity-space work).
     */
    fluxes.x_final = fluxes.x_high;
    close_periodic_face_blocks(fluxes.x_final, nxl,
                               static_cast<int>(Param::Nvmu),
                               mpi_rank, mpi_size, 505);

    /*
     * Positivity-constrained u-force flux limiter.
     *
     * The centered midpoint u-flux can produce negative f when du_div is large.
     * Instead of clipping f to zero after the fact, we limit the u_force_face
     * fluxes proportionally so that the discrete update respects positivity.
     * The limiter is face-consistent: each u_face takes min(alpha_L, alpha_R)
     * from its two adjacent cells.
     */
    bkg_new = bkg_n;
    const double eps_tol_base = 1.0e-12;

    // Precompute dt*inv_shell/2 and other per-iv constants once
    double half_dt_inv_shell[Param::Nv];
    double ke_per_mass_arr[Param::Nv];
    for (int iv = 0; iv < Param::Nv; ++iv) {
        half_dt_inv_shell[iv] = 0.5 * dt * bkg_n.vgrid.inv_moment_weight[iv];
        ke_per_mass_arr[iv] =
            (bkg_n.vgrid.gamma_cells[iv] - 1.0)
          * bkg_n.mass * Const::c * Const::c;
    }

    // ---- Pass 1: compute per-cell alpha and detect if limiting is needed ----
    std::vector<double>& alpha_cell = fluxes.cell_alpha;
    std::fill(alpha_cell.begin(), alpha_cell.end(), 1.0);
    int local_any_limited = 0;
    long long local_ec_negative_cells = 0;
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(max:local_any_limited) reduction(+:local_ec_negative_cells)
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
            const double mu_div_left =
                mu_force_face[mu_xface_index(ix, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix, iv, imu)];
            const double mu_div_right =
                mu_force_face[mu_xface_index(ix + 1, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix + 1, iv, imu)];
            const double dmu_div = hdt_is * (mu_div_left + mu_div_right);
            const double contrib[4] = {
                +hdt_is * u_force_face[u_xface_index(ix, iv + 1, imu)],
                -hdt_is * u_force_face[u_xface_index(ix, iv, imu)],
                +hdt_is * u_force_face[u_xface_index(ix + 1, iv + 1, imu)],
                -hdt_is * u_force_face[u_xface_index(ix + 1, iv, imu)]
            };
            double du_div_out = 0.0;
            double du_div_ec = 0.0;
            for (int c = 0; c < 4; ++c) {
                if (contrib[c] > 0.0) du_div_out += contrib[c];
                du_div_ec += contrib[c];
            }
            const double f_ec = f0 - dx_div - dmu_div - du_div_ec;
            if (!std::isfinite(f_ec) || f_ec >= 0.0) continue;
            ++local_ec_negative_cells;
            if (du_div_out <= 0.0) continue;

            const double local_scale = std::max(1.0, std::fabs(f0));
            const double eps_tol = eps_tol_base * local_scale;
            const double needed_reduction = -f_ec + eps_tol;
            const double alpha =
                std::max(0.0, std::min(1.0,
                         1.0 - needed_reduction / du_div_out));
            alpha_cell[ix_k] = alpha;
            local_any_limited = 1;
        }
    }
    int global_any_limited = local_any_limited;
    MPI_Allreduce(MPI_IN_PLACE, &global_any_limited, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    (void)local_ec_negative_cells;

    const bool need_face_limiting = (global_any_limited != 0);
    std::vector<double> u_force_ec;
    if (need_face_limiting) {
        u_force_ec = u_force_face;
    }
    std::vector<double> alpha_left_ghost;
    if (need_face_limiting && nxl > 0) {
        alpha_left_ghost.assign(Param::Nvmu, 1.0);
        if (mpi_size <= 1) {
            const size_t last_base =
                static_cast<size_t>(nxl - 1) * Param::Nvmu;
            std::copy(alpha_cell.begin() + last_base,
                      alpha_cell.begin() + last_base + Param::Nvmu,
                      alpha_left_ghost.begin());
        } else {
            std::vector<double> alpha_send_right(Param::Nvmu, 1.0);
            const size_t last_base =
                static_cast<size_t>(nxl - 1) * Param::Nvmu;
            std::copy(alpha_cell.begin() + last_base,
                      alpha_cell.begin() + last_base + Param::Nvmu,
                      alpha_send_right.begin());
            const int left_rank = (mpi_rank + mpi_size - 1) % mpi_size;
            const int right_rank = (mpi_rank + 1) % mpi_size;
            MPI_Sendrecv(alpha_send_right.data(),
                         static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                         right_rank, 508,
                         alpha_left_ghost.data(),
                         static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                         left_rank, 508, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
        }
    }

    long long local_u_face_active = 0;
    long long local_u_face_total = 0;
    double local_u_alpha_min = 1.0;

    if (need_face_limiting) {
    // ---- Pass 2: face-level alpha and scale u_force_face ----
    #pragma omp parallel for schedule(static) \
        reduction(+:local_u_face_active,local_u_face_total) \
        reduction(min:local_u_alpha_min)
    for (int iface = 0; iface < nxl; ++iface) {
        const double* alpha_left =
            (iface == 0)
            ? alpha_left_ghost.data()
            : &alpha_cell[static_cast<size_t>(iface - 1) * Param::Nvmu];
        const double* alpha_right =
            &alpha_cell[static_cast<size_t>(iface) * Param::Nvmu];
        for (int face = 0; face <= Param::Nv; ++face) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t face_idx = u_xface_index(iface, face, imu);
                double& ff = u_force_face[face_idx];
                double alpha_face = 1.0;
                if (ff > 0.0 && face > 0) {
                    const size_t k =
                        static_cast<size_t>(face - 1) * Param::Nmu
                      + static_cast<size_t>(imu);
                    alpha_face = std::min(alpha_left[k], alpha_right[k]);
                } else if (ff < 0.0 && face < Param::Nv) {
                    const size_t k =
                        static_cast<size_t>(face) * Param::Nmu
                      + static_cast<size_t>(imu);
                    alpha_face = std::min(alpha_left[k], alpha_right[k]);
                }
                alpha_face = std::max(0.0, std::min(1.0, alpha_face));
                ff *= alpha_face;
                if (alpha_face < 0.999999) ++local_u_face_active;
                ++local_u_face_total;
                local_u_alpha_min =
                    std::min(local_u_alpha_min, alpha_face);
            }
        }
    }
    close_periodic_face_blocks(u_force_face, nxl,
                               (Param::Nv + 1) * Param::Nmu,
                               mpi_rank, mpi_size, 506);
    } // need_face_limiting

    // ---- Pass 3: final update (uses scaled u_force_face if limiting was needed) ----
    int local_bad_update_centered = 0;
    int local_bad_negative_hard = 0;
    double local_neg_mass_defect = 0.0;
    double local_neg_energy_defect = 0.0;
    double local_limiter_mass_delta = 0.0;
    double local_limiter_momentum_delta = 0.0;
    double local_limiter_energy_delta = 0.0;
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
    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:local_neg_mass_defect,local_neg_energy_defect, \
                    local_limiter_mass_delta,local_limiter_momentum_delta, \
                    local_limiter_energy_delta, \
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
            const double u_div_left =
                u_force_face[u_xface_index(ix, iv + 1, imu)]
              - u_force_face[u_xface_index(ix, iv, imu)];
            const double u_div_right =
                u_force_face[u_xface_index(ix + 1, iv + 1, imu)]
              - u_force_face[u_xface_index(ix + 1, iv, imu)];
            const double du_div = hdt_is * (u_div_left + u_div_right);
            double du_div_ec = du_div;
            if (need_face_limiting) {
                const double u_div_left_ec =
                    u_force_ec[u_xface_index(ix, iv + 1, imu)]
                  - u_force_ec[u_xface_index(ix, iv, imu)];
                const double u_div_right_ec =
                    u_force_ec[u_xface_index(ix + 1, iv + 1, imu)]
                  - u_force_ec[u_xface_index(ix + 1, iv, imu)];
                du_div_ec = hdt_is * (u_div_left_ec + u_div_right_ec);
            }
            const double mu_div_left =
                mu_force_face[mu_xface_index(ix, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix, iv, imu)];
            const double mu_div_right =
                mu_force_face[mu_xface_index(ix + 1, iv, imu + 1)]
              - mu_force_face[mu_xface_index(ix + 1, iv, imu)];
            const double dmu_div = hdt_is * (mu_div_left + mu_div_right);
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
            const bool limiter_cell_active =
                fluxes.cell_alpha[static_cast<size_t>(ix) * Param::Nvmu + k]
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
            } else if (updated < 0.0) {
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
            if (updated < 0.0) {
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
        fluxes.limiter_active_fraction = fluxes.u_force_alpha_active_frac;
        fluxes.limiter_min_alpha = fluxes.u_force_alpha_min;
        fluxes.limiter_active_fraction_core =
            fluxes.u_force_alpha_active_frac;
        fluxes.limiter_active_fraction_boundary = 0.0;
        fluxes.limiter_min_alpha_core = fluxes.u_force_alpha_min;
        fluxes.limiter_min_alpha_boundary = 1.0;
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
    std::fill(fluxes.cell_alpha.begin(), fluxes.cell_alpha.end(), 1.0);
    return;

}

void VlasovAmpereMidpointSolver::update_flux_current(
    const Species& sp,
    const SpatialGrid& sg,
    const FluxPack& fluxes,
    Species& bkg_new) const
{
    const int nxl = sg.nx_local;
    resize_or_zero(bkg_new.current_face_x, static_cast<size_t>(nxl + 1));
    resize_or_zero(bkg_new.current_x, static_cast<size_t>(nxl));

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
    const int max_iters = 80;
    const double field_tol = 1.0e-6;
    const double current_tol = 1.0e-5;
    const double f_tol = 1.0e-6;
    const double omega_min = 0.05;
    double omega = 0.5;
    double previous_residual = -1.0;
    bool have_previous_current = false;

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

        BeamPIC beam_trial = beam_n;
        beam_trial.begin_step(sg, dt);
        beam_trial.inject(sg, fields_mid, dt, time, mpi_rank, mpi_size);
        const double beam_ke_before = beam_trial.total_kinetic_energy();
        beam_trial.push(sg, fields_mid, dt, mpi_rank, mpi_size);
        const double beam_ke_after = beam_trial.total_kinetic_energy();
        beam_trial.deposit_density(sg, mpi_rank, mpi_size);
        beam_trial.finalize_charge_conserving_current(sg, dt, mpi_rank,
                                                      mpi_size);

        for (int iface = 0; iface <= nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            const double jb = (slot < bkg_new.current_face_x.size())
                            ? bkg_new.current_face_x[slot] : 0.0;
            const double jbeam = (slot < beam_trial.current_face_x.size())
                               ? beam_trial.current_face_x[slot] : 0.0;
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
                (slot < beam_trial.current_face_x.size())
                ? beam_trial.current_face_x[slot] : 0.0;
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
        const double local_beam_continuity =
            std::max(beam_trial.last_continuity_linf_error(),
                     beam_trial.last_boundary_flux_error());
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
        result.x_negative_mass_before_repair =
            fluxes.negative_mass_before_repair;
        result.x_mass_added_by_positivity_repair =
            fluxes.mass_added_by_positivity_repair;
        result.positivity_energy_defect =
            fluxes.positivity_energy_defect;
        result.positivity_mass_defect =
            fluxes.positivity_mass_defect;
        result.u_force_alpha_min = fluxes.u_force_alpha_min;
        result.u_force_alpha_active_frac =
            fluxes.u_force_alpha_active_frac;
        result.f_neg_min = fluxes.f_neg_min;
        result.f_neg_ratio_max = fluxes.f_neg_ratio_max;
        result.f_neg_mass_total = fluxes.f_neg_mass_total;
        result.f_neg_cell_count = fluxes.f_neg_cell_count;
        result.f_neg_ix = fluxes.f_neg_ix;
        result.f_neg_iv = fluxes.f_neg_iv;
        result.f_neg_imu = fluxes.f_neg_imu;

        const bool finite_state =
            global_errors[12] == 0.0 &&
            check_finite_state(bkg_new, beam_trial, fields_new, ex_mid_next,
                               j_total_next);
        if (!finite_state) {
            if (mpi_rank == 0 && global_errors[12] != 0.0) {
                if (fluxes.finite_flux_has_failure != 0) {
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
        const bool converged =
            protected_quantities_converged && global_f_error < f_tol;
        const bool soft_accepted = false;

        if (converged) {
            result.species_np1 = bkg_new;
            result.beam_np1 = beam_trial;
            result.fields_np1 = fields_new;
            result.j_bkg_face_mid = bkg_new.current_face_x;
            result.j_beam_face_mid = beam_trial.current_face_x;
            result.j_total_face_mid = j_total_face;
            result.j_bkg_energy_debug_face = bkg_new.current_face_x;
            result.delta_ke_bkg = local_end_values[1] - local_ke_start;
            result.delta_ke_beam = beam_ke_after - beam_ke_before;
            result.field_work_bkg =
                integrate_face_work(bkg_new.current_face_x, fields_mid, sg,
                                    dt);
            result.field_work_beam =
                integrate_face_work(beam_trial.current_face_x, fields_mid, sg,
                                    dt);
            result.energy_residual_bkg =
                result.delta_ke_bkg + result.field_work_bkg;
            build_current_diagnostics(result.j_bkg_face_mid,
                                      result.j_bkg_energy_debug_face,
                                      result.j_bkg_face_mid,
                                      fields_mid, sg, result.delta_ke_bkg,
                                      dt, result.current_diag);
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
            result.x_negative_mass_before_repair =
                fluxes.negative_mass_before_repair;
            result.x_mass_added_by_positivity_repair =
                fluxes.mass_added_by_positivity_repair;
            result.positivity_energy_defect =
                fluxes.positivity_energy_defect;
            result.positivity_mass_defect =
                fluxes.positivity_mass_defect;
            result.u_force_alpha_min = fluxes.u_force_alpha_min;
            result.u_force_alpha_active_frac =
                fluxes.u_force_alpha_active_frac;
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
