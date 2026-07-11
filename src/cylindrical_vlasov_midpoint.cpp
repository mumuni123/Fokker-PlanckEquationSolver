#include "vlasov_ampere_midpoint.h"
#include "ppm_ctu_reconstruction.h"

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
    std::vector<double> outgoing(static_cast<size_t>(nxl) * Param::Nvmu, 0.0);
    std::vector<double> donor_ratio(static_cast<size_t>(nxl) * Param::Nvmu, 1.0);
    std::vector<double> last_ratio(Param::Nvmu, 1.0);
    std::vector<double> left_ratio(Param::Nvmu, 1.0);
    std::vector<double> fx_final(xflux_size, 0.0), fu_final(uflux_size, 0.0);
    std::vector<double> cu_final(uflux_size, 0.0);
    Species low_state_buffer = bkg_n;
    std::vector<double> e_mid_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> next_e_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> integrated_jn_buffer(nface, 0.0);
    std::vector<double> integrated_je_buffer(static_cast<size_t>(nxl), 0.0);
    std::vector<double> je_cell_buffer(static_cast<size_t>(nxl), 0.0);
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
            result.x_low_failure_kind = 3; // genuine combined CFL failure
            return result;
        }

        Species work = bkg_n;
        std::vector<double>& integrated_jn = integrated_jn_buffer;
        std::vector<double>& integrated_je_cell = integrated_je_buffer;
        std::fill(integrated_jn.begin(), integrated_jn.end(), 0.0);
        std::fill(integrated_je_cell.begin(), integrated_je_cell.end(), 0.0);
        double limiter_energy_change = 0.0;
        double limiter_energy_positive = 0.0;
        double limiter_energy_negative = 0.0;
        double limiter_energy_core = 0.0;
        double limiter_energy_boundary = 0.0;
        double u_boundary_energy = 0.0;
        double limiter_faces = 0.0;
        double limiter_active = 0.0;
        double limiter_min = 1.0;

        const double h = dt / nsub;
        for (int sub = 0; sub < nsub; ++sub) {
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
            #pragma omp parallel for reduction(min:low_min) reduction(max:low_scale)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double m = low_state.f[mass_index(ng + ix, j, k)];
                    low_min = std::min(low_min, m);
                    low_scale = std::max(low_scale, std::fabs(m));
                }
            double low_guard[2] = {low_min, low_scale};
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[0], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &low_guard[1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            const double mass_roundoff = 128.0 * std::numeric_limits<double>::epsilon() *
                                         std::max(1.0, low_guard[1]);
            if (!std::isfinite(low_guard[0]) || low_guard[0] < -mass_roundoff) {
                result.failed = true;
                result.state_advanced = 0;
                return result;
            }

            std::fill(outgoing.begin(), outgoing.end(), 0.0);
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
                    outgoing[static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k)] = h *
                        (std::max(0.0, ax_right) + std::max(0.0, -ax_left) +
                         std::max(0.0, au_upper) + std::max(0.0, -au_lower));
                }
            }

            std::fill(donor_ratio.begin(), donor_ratio.end(), 1.0);
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t p = static_cast<size_t>(ix) * Param::Nvmu + idx2(j, k);
                    if (outgoing[p] > 0.0) donor_ratio[p] = std::max(0.0,
                        std::min(1.0, low_state.f[mass_index(ng + ix, j, k)] /
                                      outgoing[p]));
                }

            // MPI face owners use a single donor budget on each shared face.
            for (size_t p = 0; p < Param::Nvmu; ++p) {
                last_ratio[p] = donor_ratio[static_cast<size_t>(nxl - 1) * Param::Nvmu + p];
            }
            if (mpi_size == 1) {
                left_ratio = last_ratio;
            } else {
                const int left = (mpi_rank + mpi_size - 1) % mpi_size;
                const int right = (mpi_rank + 1) % mpi_size;
                MPI_Sendrecv(last_ratio.data(), static_cast<int>(Param::Nvmu), MPI_DOUBLE,
                             right, 904, left_ratio.data(), static_cast<int>(Param::Nvmu),
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
                    const int donor_ix = (anti >= 0.0) ? iface - 1 : iface;
                    const size_t p = idx2(j, k);
                    const double alpha = (donor_ix < 0) ? left_ratio[p] :
                        donor_ratio[static_cast<size_t>(donor_ix) * Param::Nvmu + p];
                    fx_final[id] += alpha * anti;
                    limiter_faces += 1.0;
                    if (alpha < 1.0 - 1.0e-14) limiter_active += 1.0;
                    limiter_min = std::min(limiter_min, alpha);
                }
            }
            close_periodic_face_blocks(fx_final, nxl, Param::Nvmu,
                                       mpi_rank, mpi_size, 905);
            #pragma omp parallel for schedule(static) reduction(+:limiter_energy_change,limiter_energy_positive,limiter_energy_negative,limiter_energy_core,limiter_energy_boundary)
            for (int ix = 0; ix < nxl; ++ix) for (int jf = 1; jf < Param::Nv; ++jf)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = uface_index(ix, jf, k);
                    const double anti = fu_high[id] - fu_low[id];
                    const int donor_j = (anti >= 0.0) ? jf - 1 : jf;
                    const double alpha = donor_ratio[
                        static_cast<size_t>(ix) * Param::Nvmu + idx2(donor_j, k)];
                    fu_final[id] += alpha * anti;
                    cu_final[id] += alpha * (cu_high[id] - cu_low[id]);
                    const double dlim = h *
                        (bkg_n.cgrid.kinetic_energy[idx2(jf, k)] -
                         bkg_n.cgrid.kinetic_energy[idx2(jf - 1, k)]) *
                        (fu_final[id] - fu_high[id]);
                    limiter_energy_change += dlim;
                    if (dlim >= 0.0) limiter_energy_positive += dlim;
                    else limiter_energy_negative += dlim;
                    const double x = sg.x(ng + ix);
                    if (x >= 0.1 * Const::micro && x <= Param::Lx - 0.1 * Const::micro)
                        limiter_energy_core += dlim;
                    else limiter_energy_boundary += dlim;
                }

            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k) {
                    work.f[mass_index(ng + ix, j, k)] -= h *
                        (fx_final[xface_index(ix + 1, j, k)] -
                         fx_final[xface_index(ix, j, k)] +
                         fu_final[uface_index(ix, j + 1, k)] -
                        fu_final[uface_index(ix, j, k)]);
                }

            double final_min = std::numeric_limits<double>::infinity();
            #pragma omp parallel for reduction(min:final_min)
            for (int ix = 0; ix < nxl; ++ix) for (int j = 0; j < Param::Nv; ++j)
                for (int k = 0; k < Param::Nmu; ++k)
                    final_min = std::min(final_min, work.f[mass_index(ng + ix, j, k)]);
            MPI_Allreduce(MPI_IN_PLACE, &final_min, 1, MPI_DOUBLE, MPI_MIN,
                          MPI_COMM_WORLD);
            if (!std::isfinite(final_min) || final_min < -mass_roundoff) {
                result.failed = true;
                result.state_advanced = 0;
                return result;
            }

            #pragma omp parallel for schedule(static)
            for (int iface = 0; iface <= nxl; ++iface) {
                double gamma = 0.0;
                for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k)
                    gamma += fx_final[xface_index(iface, j, k)];
                integrated_jn[iface] += h * bkg_n.charge * gamma;
            }
            #pragma omp parallel for schedule(static)
            for (int ix = 0; ix < nxl; ++ix) {
                double je = 0.0;
                for (int jf = 1; jf < Param::Nv; ++jf) for (int k = 0; k < Param::Nmu; ++k) {
                    const double dke = bkg_n.cgrid.kinetic_energy[idx2(jf, k)] -
                                       bkg_n.cgrid.kinetic_energy[idx2(jf - 1, k)];
                    je += dke * (bkg_n.charge / (bkg_n.mass * Const::c)) *
                          cu_final[uface_index(ix, jf, k)] / sg.dx;
                }
                integrated_je_cell[ix] += h * je;
            }
            for (int k = 0; k < Param::Nmu; ++k) {
                const double ke_lo = bkg_n.cgrid.kinetic_energy[idx2(0, k)];
                const double ke_hi = bkg_n.cgrid.kinetic_energy[idx2(Param::Nv - 1, k)];
                for (int ix = 0; ix < nxl; ++ix) {
                    u_boundary_energy += h * (ke_hi * fu_final[uface_index(ix, Param::Nv, k)] -
                                               ke_lo * fu_final[uface_index(ix, 0, k)]);
                }
            }
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
        double local_jmax[3] = {0.0, 0.0, 0.0};
        for (int iface = 0; iface < nxl; ++iface) {
            const double jn = result.j_bkg_face_mid[iface];
            const double je = result.j_bkg_energy_debug_face[iface];
            local_jmax[0] = std::max(local_jmax[0], std::fabs(jn));
            local_jmax[1] = std::max(local_jmax[1], std::fabs(je));
            local_jmax[2] = std::max(local_jmax[2], std::fabs(jn - je));
        }
        MPI_Allreduce(MPI_IN_PLACE, local_jmax, 3, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
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

        if (mpi_rank == 0 && (iter == 0 || result.residual_E < field_tol)) {
            static std::ofstream ledger("output/cylindrical_energy_ledger.dat");
            if (ledger.tellp() == std::streampos(0)) {
                ledger << "# time_fs iter dK_bkg dK_beam W_bkg W_beam R_pair "
                       << "Dlim_pos Dlim_neg Dlim_core Dlim_boundary R_beam_lag BuK\n";
                ledger << std::scientific << std::setprecision(10);
            }
            ledger << time / Const::femto << " " << (iter + 1) << " "
                   << result.delta_ke_bkg << " " << result.delta_ke_beam << " "
                   << result.field_work_bkg << " " << result.field_work_beam << " "
                   << result.energy_pair_residual_bkg << " "
                   << result.limiter_energy_defect_positive << " "
                   << result.limiter_energy_defect_negative << " "
                   << result.limiter_energy_defect_core << " "
                   << result.limiter_energy_defect_boundary << " "
                   << result.beam_lag_energy_residual << " "
                   << result.flux_defect[1].boundary_energy_loss << "\n";
            ledger.flush();
        }

        if (!all_finite(work, fields_new, result.j_total_face_mid)) {
            result.failed = true;
            result.state_advanced = 0;
            return result;
        }
        if (result.residual_E < field_tol && result.residual_J_bkg < current_tol &&
            result.beam_continuity_residual < 1.0e-6) {
            result.species_np1 = work;
            result.beam_np1 = beam_predictor;
            result.fields_np1 = fields_new;
            result.converged = true;
            result.state_advanced = 1;
            return result;
        }

        for (int iface = 0; iface < nxl; ++iface) {
            e_end[iface] = (1.0 - omega) * e_end[iface] + omega * next_e[iface];
            previous_j[iface] = result.j_bkg_face_mid[iface];
        }
        have_previous = true;
        guess = work;
    }

    result.species_np1 = last_work;
    result.beam_np1 = beam_predictor;
    result.fields_np1 = last_fields;
    result.soft_unconverged = true;
    result.soft_accepted = true;
    result.state_advanced = 1;
    return result;
}
