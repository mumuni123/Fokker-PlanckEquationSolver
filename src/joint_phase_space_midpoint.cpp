#include "joint_phase_space_midpoint.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mpi.h>
#include <stdexcept>
#include <string>

namespace {

void require_finite_vector(const std::vector<double>& values,
                           const char* name)
{
    for (size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i]))
            throw std::runtime_error(std::string(name) + " contains nonfinite value");
    }
}

} // namespace

std::size_t JointPhaseSpaceMidpointOperator::x_index(
    int iface, int j, int k, int nv, int nmu)
{
    return (static_cast<size_t>(iface) * static_cast<size_t>(nv) +
            static_cast<size_t>(j)) * static_cast<size_t>(nmu) +
           static_cast<size_t>(k);
}

std::size_t JointPhaseSpaceMidpointOperator::u_index(
    int ix, int jface, int k, int nv, int nmu)
{
    return (static_cast<size_t>(ix) * static_cast<size_t>(nv + 1) +
            static_cast<size_t>(jface)) * static_cast<size_t>(nmu) +
           static_cast<size_t>(k);
}

std::size_t JointPhaseSpaceMidpointOperator::cell_index(
    int ix, int j, int k, int nv, int nmu)
{
    return (static_cast<size_t>(ix) * static_cast<size_t>(nv) +
            static_cast<size_t>(j)) * static_cast<size_t>(nmu) +
           static_cast<size_t>(k);
}

std::vector<double> JointPhaseSpaceMidpointOperator::build_hamiltonian_velocity(
    const CylindricalVelocityGrid& vg)
{
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    if (nv < 2 || nmu <= 0 || vg.kinetic_energy.size() !=
        static_cast<size_t>(nv * nmu)) {
        throw std::runtime_error("J0 Hamiltonian velocity grid dimensions invalid");
    }

    // This velocity is not an analytic cell velocity.  It is the exact
    // algebraic transpose of the J1 centered u_parallel face trace under the
    // cell-integrated-mass representation.  Do not replace it with vg.vx,
    // vg.vx_energy_conjugate_cell, or a centered derivative on the
    // nonuniform grid.
    if (vg.upar_widths.size() != static_cast<size_t>(nv))
        throw std::runtime_error("J1 Hamiltonian velocity widths invalid");
    std::vector<double> result(static_cast<size_t>(nv * nmu), 0.0);
    for (int jf = 1; jf < nv; ++jf) {
        const int jl = jf - 1;
        const int jr = jf;
        const double du_left = vg.upar_widths[static_cast<size_t>(jl)];
        const double du_right = vg.upar_widths[static_cast<size_t>(jr)];
        if (!(du_left > 0.0) || !(du_right > 0.0) ||
            !std::isfinite(du_left) || !std::isfinite(du_right)) {
            throw std::runtime_error("J1 Hamiltonian velocity width invalid");
        }
        for (int k = 0; k < nmu; ++k) {
            const size_t left_id = static_cast<size_t>(jl * nmu + k);
            const size_t right_id = static_cast<size_t>(jr * nmu + k);
            const double delta_k = vg.kinetic_energy[right_id] -
                vg.kinetic_energy[left_id];
            result[left_id] += delta_k /
                (2.0 * Const::me * Const::c * du_left);
            result[right_id] += delta_k /
                (2.0 * Const::me * Const::c * du_right);
        }
    }
    for (size_t i = 0; i < result.size(); ++i) {
        if (!std::isfinite(result[i]) ||
            std::fabs(result[i]) > Const::c *
                (1.0 + 4096.0 * std::numeric_limits<double>::epsilon())) {
            throw std::runtime_error("J1 Hamiltonian velocity is invalid");
        }
    }
    return result;
}

JointPhaseSpaceFluxBundle
JointPhaseSpaceMidpointOperator::build_periodic_center_flux(
    const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
    const std::vector<double>& m_mid, const std::vector<double>& e_cell,
    double dt)
{
    const int nx = sg.nx_global;
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    if (nx <= 0 || nv < 2 || nmu <= 0 || !(sg.dx > 0.0) ||
        !std::isfinite(sg.dx) || !(dt > 0.0) || !std::isfinite(dt)) {
        throw std::runtime_error("J0 flux geometry or time step invalid");
    }
    const size_t cell_count = static_cast<size_t>(nx * nv * nmu);
    if (m_mid.size() != cell_count || e_cell.size() != static_cast<size_t>(nx)) {
        throw std::runtime_error("J0 flux input size mismatch");
    }
    require_finite_vector(m_mid, "J0 M_mid");
    require_finite_vector(e_cell, "J0 E_cell");

    const std::vector<double> hamiltonian_velocity =
        build_hamiltonian_velocity(vg);
    JointPhaseSpaceFluxBundle bundle;
    bundle.nx = nx;
    bundle.nv = nv;
    bundle.nmu = nmu;
    bundle.dt = dt;
    bundle.dx = sg.dx;
    bundle.x_flux_rate.assign(static_cast<size_t>((nx + 1) * nv * nmu), 0.0);
    bundle.u_flux_rate.assign(static_cast<size_t>(nx * (nv + 1) * nmu), 0.0);
    bundle.charge_current_face.assign(static_cast<size_t>(nx + 1), 0.0);
    bundle.energy_current_cell.assign(static_cast<size_t>(nx), 0.0);
    bundle.mass_delta_x.assign(cell_count, 0.0);
    bundle.mass_delta_u.assign(cell_count, 0.0);
    bundle.mass_delta_total.assign(cell_count, 0.0);
    bundle.kinetic_work_cell.assign(static_cast<size_t>(nx), 0.0);

    // The unique periodic face is shared by iface=0 and iface=nx.  This is a
    // test-only topology choice; J2 will replace it with production boundary
    // traces after J1 has passed.
    for (int iface = 0; iface <= nx; ++iface) {
        const int left_cell = iface == 0 ? nx - 1 : iface - 1;
        const int right_cell = iface == nx ? 0 : iface;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t left = cell_index(left_cell, j, k, nv, nmu);
                const size_t right = cell_index(right_cell, j, k, nv, nmu);
                const size_t q = static_cast<size_t>(j * nmu + k);
                const double trace = 0.5 * (m_mid[left] + m_mid[right]);
                // M is integrated over dx.  Dividing by dx converts it to the
                // x flux-rate density used by charge current and FV update.
                bundle.x_flux_rate[x_index(iface, j, k, nv, nmu)] =
                    hamiltonian_velocity[q] * trace / sg.dx;
            }
        }
    }

    for (int ix = 0; ix < nx; ++ix) {
        for (int jface = 1; jface < nv; ++jface) {
            const double du_center = vg.upar_cells[jface] -
                                     vg.upar_cells[jface - 1];
            const double a = e_cell[static_cast<size_t>(ix)] *
                (-Const::qe) / (Const::me * Const::c);
            for (int k = 0; k < nmu; ++k) {
                const size_t left = cell_index(ix, jface - 1, k, nv, nmu);
                const size_t right = cell_index(ix, jface, k, nv, nmu);
                const double f_trace = 0.5 *
                    (m_mid[left] / vg.upar_widths[jface - 1] +
                     m_mid[right] / vg.upar_widths[jface]);
                // Integrating f over dx and the perpendicular ring gives the
                // mass rate across a u face.  Endpoint faces stay exactly 0.
                bundle.u_flux_rate[u_index(ix, jface, k, nv, nmu)] =
                    a * f_trace;

                const size_t q = static_cast<size_t>(jface * nmu + k);
                const double delta_k = vg.kinetic_energy[q] -
                    vg.kinetic_energy[static_cast<size_t>((jface - 1) * nmu + k)];
                // The energy current is built from the actual FV face
                // difference, not from the cell Hamiltonian velocity used by
                // x transport.  This is the J0 discrete chain rule.
                const double v_energy_face = delta_k /
                    (Const::me * Const::c * du_center);
                bundle.energy_current_cell[static_cast<size_t>(ix)] +=
                    (-Const::qe) * v_energy_face * f_trace * du_center *
                    vg.uperp_ring_areas[k];
                // Accumulate the diagnostic work in extended precision.  The
                // stored flux remains double and is the production authority;
                // this only prevents a summation-order artifact in J0.
                const long double work =
                    static_cast<long double>(dt) *
                    static_cast<long double>(delta_k) *
                    static_cast<long double>(bundle.u_flux_rate[
                        u_index(ix, jface, k, nv, nmu)]);
                bundle.kinetic_work_cell[static_cast<size_t>(ix)] +=
                    static_cast<double>(work);
            }
        }
    }

    for (int iface = 0; iface <= nx; ++iface) {
        double current = 0.0;
        for (int j = 0; j < nv; ++j)
            for (int k = 0; k < nmu; ++k)
                current += (-Const::qe) *
                    bundle.x_flux_rate[x_index(iface, j, k, nv, nmu)];
        bundle.charge_current_face[static_cast<size_t>(iface)] = current;
    }

    for (int ix = 0; ix < nx; ++ix) {
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t id = cell_index(ix, j, k, nv, nmu);
                const double x_div =
                    bundle.x_flux_rate[x_index(ix + 1, j, k, nv, nmu)] -
                    bundle.x_flux_rate[x_index(ix, j, k, nv, nmu)];
                const double u_div =
                    bundle.u_flux_rate[u_index(ix, j + 1, k, nv, nmu)] -
                    bundle.u_flux_rate[u_index(ix, j, k, nv, nmu)];
                bundle.mass_delta_x[id] = -dt * x_div;
                bundle.mass_delta_u[id] = -dt * u_div;
                bundle.mass_delta_total[id] =
                    bundle.mass_delta_x[id] + bundle.mass_delta_u[id];
            }
        }
    }
    return bundle;
}

bool JointPhaseSpaceMidpointOperator::evaluate_residual(
    const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
    const std::vector<double>& m_old,
    const std::vector<double>& m_candidate,
    const std::vector<double>& e_cell, double dt,
    JointPhaseSpaceFluxBundle& bundle,
    std::vector<double>& residual,
    double& residual_linf,
    double& residual_scale)
{
    if (m_old.size() != m_candidate.size()) return false;
    std::vector<double> midpoint(m_old.size(), 0.0);
    for (size_t i = 0; i < midpoint.size(); ++i) {
        if (!std::isfinite(m_old[i]) || !std::isfinite(m_candidate[i]))
            return false;
        midpoint[i] = 0.5 * (m_old[i] + m_candidate[i]);
    }
    try {
        bundle = build_periodic_center_flux(sg, vg, midpoint, e_cell, dt);
    } catch (const std::exception&) {
        return false;
    }
    residual.assign(m_candidate.size(), 0.0);
    residual_linf = 0.0;
    residual_scale = 1.0;
    for (size_t i = 0; i < residual.size(); ++i) {
        residual[i] = m_candidate[i] - m_old[i] - bundle.mass_delta_total[i];
        if (!std::isfinite(residual[i])) return false;
        residual_linf = std::max(residual_linf, std::fabs(residual[i]));
        residual_scale = std::max(residual_scale,
            std::fabs(m_old[i]) + std::fabs(m_candidate[i]));
        // J1 evaluates the algebraic midpoint residual on signed Newton
        // trials.  Positivity is an acceptance property checked by the
        // integrator after convergence, not a residual-domain restriction:
        // exact zero Maxwellian tail cells would otherwise make every
        // fraction-to-boundary step identically zero before the center-trace
        // operator can be diagnosed.
    }
    return true;
}

bool JointPhaseSpaceMidpointOperator::apply_block_diagonal_preconditioner(
    const CylindricalVelocityGrid& vg, const std::vector<double>& e_cell,
    double dx, double dt, const std::vector<double>& residual,
    std::vector<double>& preconditioned)
{
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    if (nv <= 0 || nmu <= 0 || dx <= 0.0 || dt <= 0.0 || e_cell.empty() ||
        residual.size() != e_cell.size() * static_cast<size_t>(nv * nmu))
        return false;
    const std::vector<double> velocity = build_hamiltonian_velocity(vg);
    preconditioned.assign(residual.size(), 0.0);
    for (size_t ix = 0; ix < e_cell.size(); ++ix) {
        const double accel = std::fabs(Const::qe * e_cell[ix] /
                                       (Const::me * Const::c));
        for (int j = 0; j < nv; ++j) {
            const double du = std::max(vg.upar_widths[static_cast<size_t>(j)],
                                       std::numeric_limits<double>::min());
            for (int k = 0; k < nmu; ++k) {
                const size_t q = static_cast<size_t>(j * nmu + k);
                const size_t id = ix * static_cast<size_t>(nv * nmu) + q;
                const double diagonal = 1.0 + dt *
                    (std::fabs(velocity[q]) / dx + accel / du);
                if (!std::isfinite(diagonal) || diagonal <= 0.0 ||
                    !std::isfinite(residual[id])) return false;
                preconditioned[id] = residual[id] / diagonal;
            }
        }
    }
    return true;
}

bool JointPhaseSpaceMidpointOperator::apply_local_block_diagonal_preconditioner(
    const CylindricalVelocityGrid& vg, const std::vector<double>& e_cell,
    double dx, double dt, const std::vector<double>& residual,
    std::vector<double>& preconditioned)
{
    return apply_block_diagonal_preconditioner(
        vg, e_cell, dx, dt, residual, preconditioned);
}

bool JointPhaseSpaceMidpointOperator::evaluate_local_residual(
    const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
    const std::vector<double>& m_old_local,
    const std::vector<double>& m_candidate_local,
    const std::vector<double>& e_cell_local, double dt,
    int mpi_rank, int mpi_size,
    JointPhaseSpaceFluxBundle& bundle,
    std::vector<double>& residual,
    double& residual_linf,
    double& residual_scale,
    bool allow_negative_probe)
{
    const int nx = sg.nx_local;
    const int nv = static_cast<int>(vg.upar_cells.size());
    const int nmu = static_cast<int>(vg.uperp_cells.size());
    const size_t count = static_cast<size_t>(nx * nv * nmu);
    if (nx <= 0 || nv < 2 || nmu <= 0 || m_old_local.size() != count ||
        m_candidate_local.size() != count ||
        e_cell_local.size() != static_cast<size_t>(nx)) return false;
    for (size_t i = 0; i < m_old_local.size(); ++i)
        if (!std::isfinite(m_old_local[i]) || !std::isfinite(m_candidate_local[i]))
            return false;
    std::vector<double> m_mid_local(count, 0.0);
    for (size_t i = 0; i < count; ++i) {
        m_mid_local[i] = 0.5 *
            (m_old_local[i] + m_candidate_local[i]);
        if (!std::isfinite(m_mid_local[i])) return false;
    }
    for (size_t i = 0; i < e_cell_local.size(); ++i)
        if (!std::isfinite(e_cell_local[i])) return false;

    const std::vector<double> velocity = build_hamiltonian_velocity(vg);
    const int nq = nv * nmu;
    std::vector<double> left_ghost(static_cast<size_t>(nq), 0.0);
    std::vector<double> right_ghost(static_cast<size_t>(nq), 0.0);
    if (mpi_size == 1) {
        std::copy(m_mid_local.end() - nq, m_mid_local.end(),
                  left_ghost.begin());
        std::copy(m_mid_local.begin(), m_mid_local.begin() + nq,
                  right_ghost.begin());
    } else {
        const int left = mpi_rank > 0 ? mpi_rank - 1 : mpi_size - 1;
        const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : 0;
        MPI_Sendrecv(m_mid_local.data() +
                         static_cast<size_t>(nx - 1) * nq, nq, MPI_DOUBLE,
                     right, 7311, left_ghost.data(), nq, MPI_DOUBLE,
                     left, 7311, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(m_mid_local.data(), nq, MPI_DOUBLE, left, 7312,
                     right_ghost.data(), nq, MPI_DOUBLE, right, 7312,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    bundle.nx = nx;
    bundle.nv = nv;
    bundle.nmu = nmu;
    bundle.dt = dt;
    bundle.dx = sg.dx;
    bundle.x_flux_rate.assign(static_cast<size_t>((nx + 1) * nq), 0.0);
    bundle.u_flux_rate.assign(static_cast<size_t>(nx * (nv + 1) * nmu), 0.0);
    bundle.charge_current_face.assign(static_cast<size_t>(nx + 1), 0.0);
    bundle.energy_current_cell.assign(static_cast<size_t>(nx), 0.0);
    bundle.mass_delta_x.assign(count, 0.0);
    bundle.mass_delta_u.assign(count, 0.0);
    bundle.mass_delta_total.assign(count, 0.0);
    bundle.kinetic_work_cell.assign(static_cast<size_t>(nx), 0.0);

    for (int iface = 0; iface <= nx; ++iface) {
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const int q = j * nmu + k;
                const double left = iface == 0
                    ? left_ghost[static_cast<size_t>(q)]
                    : m_mid_local[static_cast<size_t>(iface - 1) * nq + q];
                const double right = iface == nx
                    ? right_ghost[static_cast<size_t>(q)]
                    : m_mid_local[static_cast<size_t>(iface) * nq + q];
                bundle.x_flux_rate[x_index(iface, j, k, nv, nmu)] =
                    velocity[static_cast<size_t>(q)] * 0.5 * (left + right) /
                    sg.dx;
            }
        }
    }
    for (int ix = 0; ix < nx; ++ix) {
        const double a = e_cell_local[static_cast<size_t>(ix)] *
            (-Const::qe) / (Const::me * Const::c);
        for (int jf = 1; jf < nv; ++jf) {
            const double du = vg.upar_cells[jf] - vg.upar_cells[jf - 1];
            for (int k = 0; k < nmu; ++k) {
                const size_t left = cell_index(ix, jf - 1, k, nv, nmu);
                const size_t right = cell_index(ix, jf, k, nv, nmu);
                const double ftrace = 0.5 *
                    (m_mid_local[left] / vg.upar_widths[jf - 1] +
                     m_mid_local[right] / vg.upar_widths[jf]);
                bundle.u_flux_rate[u_index(ix, jf, k, nv, nmu)] =
                    a * ftrace;
                // `left`/`right` above index the x-u distribution.  The
                // kinetic-energy table has velocity-slot layout only and
                // must never be indexed with the x-cell offset.
                const size_t velocity_left = static_cast<size_t>(
                    (jf - 1) * nmu + k);
                const size_t velocity_right = static_cast<size_t>(
                    jf * nmu + k);
                const double delta_k =
                    vg.kinetic_energy[velocity_right] -
                    vg.kinetic_energy[velocity_left];
                const double vface = delta_k / (Const::me * Const::c * du);
                bundle.energy_current_cell[static_cast<size_t>(ix)] +=
                    (-Const::qe) * vface * ftrace * du *
                    vg.uperp_ring_areas[k];
                bundle.kinetic_work_cell[static_cast<size_t>(ix)] +=
                    dt * delta_k * bundle.u_flux_rate[
                        u_index(ix, jf, k, nv, nmu)];
            }
        }
    }
    for (int iface = 0; iface <= nx; ++iface) {
        double current = 0.0;
        for (int q = 0; q < nq; ++q)
            current += (-Const::qe) * bundle.x_flux_rate[
                static_cast<size_t>(iface) * nq + q];
        bundle.charge_current_face[static_cast<size_t>(iface)] = current;
    }
    for (int ix = 0; ix < nx; ++ix) {
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t id = cell_index(ix, j, k, nv, nmu);
                const double xdiv = bundle.x_flux_rate[
                    x_index(ix + 1, j, k, nv, nmu)] - bundle.x_flux_rate[
                    x_index(ix, j, k, nv, nmu)];
                const double udiv = bundle.u_flux_rate[
                    u_index(ix, j + 1, k, nv, nmu)] - bundle.u_flux_rate[
                    u_index(ix, j, k, nv, nmu)];
                bundle.mass_delta_x[id] = -dt * xdiv;
                bundle.mass_delta_u[id] = -dt * udiv;
                bundle.mass_delta_total[id] =
                    bundle.mass_delta_x[id] + bundle.mass_delta_u[id];
            }
        }
    }
    residual.resize(count);
    residual_linf = 0.0;
    residual_scale = 1.0;
    for (size_t i = 0; i < count; ++i) {
        residual[i] = m_candidate_local[i] - m_old_local[i] -
                      bundle.mass_delta_total[i];
        if (!std::isfinite(residual[i])) return false;
        if (!allow_negative_probe && m_candidate_local[i] < 0.0) return false;
        residual_linf = std::max(residual_linf, std::fabs(residual[i]));
        residual_scale = std::max(residual_scale,
            std::fabs(m_old_local[i]) + std::fabs(m_candidate_local[i]));
    }
    return true;
}

JointPhaseSpaceAuditResult JointPhaseSpaceMidpointOperator::audit(
    const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
    const std::vector<double>& m_old, const std::vector<double>& m_new,
    const JointPhaseSpaceFluxBundle& bundle,
    double poisson_work_residual, double g_gstar_residual)
{
    JointPhaseSpaceAuditResult result;
    const size_t cell_count = static_cast<size_t>(bundle.nx * bundle.nv * bundle.nmu);
    const size_t velocity_count = static_cast<size_t>(bundle.nv * bundle.nmu);
    const size_t x_flux_count = static_cast<size_t>((bundle.nx + 1) *
                                                    bundle.nv * bundle.nmu);
    const size_t u_flux_count = static_cast<size_t>(bundle.nx *
                                                    (bundle.nv + 1) *
                                                    bundle.nmu);
    if (m_old.size() != cell_count || m_new.size() != cell_count ||
        bundle.mass_delta_x.size() != cell_count ||
        bundle.mass_delta_u.size() != cell_count ||
        bundle.mass_delta_total.size() != cell_count ||
        bundle.x_flux_rate.size() != x_flux_count ||
        bundle.u_flux_rate.size() != u_flux_count ||
        bundle.charge_current_face.size() !=
            static_cast<size_t>(bundle.nx + 1) ||
        bundle.energy_current_cell.size() != static_cast<size_t>(bundle.nx) ||
        bundle.kinetic_work_cell.size() != static_cast<size_t>(bundle.nx) ||
        bundle.nx <= 0 || bundle.nv < 2 || bundle.nmu <= 0 ||
        vg.upar_widths.size() != static_cast<size_t>(bundle.nv) ||
        vg.uperp_ring_areas.size() != static_cast<size_t>(bundle.nmu) ||
        vg.kinetic_energy.size() != velocity_count ||
        vg.upar_cells.size() != static_cast<size_t>(bundle.nv)) {
        return result;
    }
    double mass_residual = 0.0;
    double mass_roundoff_bound = 0.0;
    long double delta_ke = 0.0L;
    long double kinetic_work = 0.0L;
    long double kinetic_absolute_scale = 0.0L;
    double volume_error = 0.0;
    double hamiltonian_error = 0.0;
    double boundary_flux = 0.0;
    for (int ix = 0; ix < bundle.nx; ++ix) {
        for (int j = 0; j < bundle.nv; ++j) {
            for (int k = 0; k < bundle.nmu; ++k) {
                const size_t id = cell_index(ix, j, k, bundle.nv, bundle.nmu);
                if (!std::isfinite(m_old[id]) || !std::isfinite(m_new[id]) ||
                    !std::isfinite(bundle.mass_delta_x[id]) ||
                    !std::isfinite(bundle.mass_delta_u[id]) ||
                    !std::isfinite(bundle.mass_delta_total[id])) {
                    return result;
                }
                mass_residual = std::max(mass_residual,
                    std::fabs((m_new[id] - m_old[id]) -
                              bundle.mass_delta_total[id]));
                const double local_roundoff =
                    8192.0 * std::numeric_limits<double>::epsilon() *
                    (std::fabs(m_old[id]) + std::fabs(m_new[id]) +
                     std::fabs(bundle.mass_delta_total[id]));
                mass_roundoff_bound = std::max(mass_roundoff_bound,
                                               local_roundoff);
                delta_ke += static_cast<long double>(
                    vg.kinetic_energy[static_cast<size_t>(j * bundle.nmu + k)]) *
                    static_cast<long double>(bundle.mass_delta_u[id]);
                const double expected_volume = vg.upar_widths[j] *
                    vg.uperp_ring_areas[k];
                volume_error = std::max(volume_error,
                    std::fabs(vg.cell_phase_volume(j, k) - expected_volume));
            }
        }
        for (int jface = 1; jface < bundle.nv; ++jface) {
            for (int k = 0; k < bundle.nmu; ++k) {
                const std::size_t right = static_cast<std::size_t>(
                    jface * bundle.nmu + k);
                const std::size_t left = static_cast<std::size_t>(
                    (jface - 1) * bundle.nmu + k);
                const long double term = static_cast<long double>(bundle.dt) *
                    static_cast<long double>(vg.kinetic_energy[right] -
                                             vg.kinetic_energy[left]) *
                    static_cast<long double>(bundle.u_flux_rate[
                        u_index(ix, jface, k, bundle.nv, bundle.nmu)]);
                kinetic_work += term;
                kinetic_absolute_scale += std::fabs(term);
            }
        }
    }
    for (int jface = 1; jface < bundle.nv; ++jface) {
        const double du = vg.upar_cells[jface] - vg.upar_cells[jface - 1];
        for (int k = 0; k < bundle.nmu; ++k) {
            const size_t right = static_cast<size_t>(jface * bundle.nmu + k);
            const size_t left = static_cast<size_t>((jface - 1) * bundle.nmu + k);
            const double lhs = vg.kinetic_energy[right] - vg.kinetic_energy[left];
            const double v_energy_face = lhs /
                (Const::me * Const::c * du);
            const double rhs = Const::me * Const::c *
                (v_energy_face * du);
            hamiltonian_error = std::max(hamiltonian_error, std::fabs(lhs - rhs));
        }
    }
    for (int ix = 0; ix < bundle.nx; ++ix) {
        for (int k = 0; k < bundle.nmu; ++k) {
            boundary_flux = std::max(boundary_flux,
                std::max(std::fabs(bundle.u_flux_rate[
                    u_index(ix, 0, k, bundle.nv, bundle.nmu)]),
                         std::fabs(bundle.u_flux_rate[
                    u_index(ix, bundle.nv, k, bundle.nv, bundle.nmu)])));
        }
    }
    result.mass_residual = mass_residual;
    result.mass_roundoff_bound = mass_roundoff_bound;
    result.kinetic_work_residual = static_cast<double>(delta_ke - kinetic_work);
    result.kinetic_absolute_work_scale =
        static_cast<double>(kinetic_absolute_scale);
    result.poisson_work_residual = poisson_work_residual;
    result.combined_energy_residual =
        result.kinetic_work_residual + result.poisson_work_residual;
    result.g_gstar_residual = g_gstar_residual;
    result.cell_volume_residual = volume_error;
    result.hamiltonian_velocity_residual = hamiltonian_error;
    result.u_boundary_flux = boundary_flux;
    result.finite = std::isfinite(result.mass_residual) &&
        std::isfinite(result.kinetic_work_residual) &&
        std::isfinite(result.poisson_work_residual) &&
        std::isfinite(result.combined_energy_residual) &&
        std::isfinite(result.g_gstar_residual) &&
        std::isfinite(result.cell_volume_residual) &&
        std::isfinite(result.hamiltonian_velocity_residual) &&
        std::isfinite(result.u_boundary_flux);
    (void)sg;
    return result;
}

bool JointPhaseSpaceMidpointOperator::build_periodic_x_adjoint_cell_field(
    const SpatialGrid& sg,
    const std::vector<double>& pairing_face,
    int mpi_rank,
    int mpi_size,
    std::vector<double>& pairing_cell)
{
    // J1 TEST TOPOLOGY ONLY.
    //
    // The J1 x operator is periodic: global faces 0 and Nx represent
    // the same seam current, while the OpenElectrostaticSolver pairing
    // field retains two distinct non-periodic physical endpoint faces
    // with half-cell quadrature weights.
    //
    // Therefore the u-force cell field must be the exact weighted
    // transpose G* of the periodic centered x-current map G.
    //
    // Do not replace the first/last-cell formulas with
    // 0.5*(E_left_face + E_right_face).
    //
    // This is NOT the production open-boundary rule.
    // J2 must replace the periodic seam operator with the real
    // OpenBackgroundBoundary operator and derive its corresponding G*.
    pairing_cell.clear();
    if (!(sg.nx_local > 0) || !(sg.nx_global >= 2) ||
        pairing_face.size() != static_cast<size_t>(sg.nx_local) + 1 ||
        mpi_rank < 0 || mpi_rank >= mpi_size || mpi_size < 1) {
        return false;
    }
    for (size_t i = 0; i < pairing_face.size(); ++i) {
        if (!std::isfinite(pairing_face[i]))
            return false;
    }
    double endpoint_local[2] = {0.0, 0.0};
    if (mpi_rank == 0)
        endpoint_local[0] = pairing_face.front();
    if (mpi_rank == mpi_size - 1)
        endpoint_local[1] = pairing_face.back();
    double endpoint_global[2] = {0.0, 0.0};
    MPI_Allreduce(endpoint_local, endpoint_global, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double e_left = endpoint_global[0];
    const double e_right = endpoint_global[1];
    if (!std::isfinite(e_left) || !std::isfinite(e_right))
        return false;
    pairing_cell.assign(static_cast<size_t>(sg.nx_local), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int ig = sg.ix_start + ix;
        if (ig == 0) {
            pairing_cell[static_cast<size_t>(ix)] =
                0.5 * pairing_face[static_cast<size_t>(ix + 1)] +
                0.25 * (e_left + e_right);
        } else if (ig == sg.nx_global - 1) {
            pairing_cell[static_cast<size_t>(ix)] =
                0.5 * pairing_face[static_cast<size_t>(ix)] +
                0.25 * (e_left + e_right);
        } else {
            pairing_cell[static_cast<size_t>(ix)] =
                0.5 * (pairing_face[static_cast<size_t>(ix)] +
                       pairing_face[static_cast<size_t>(ix + 1)]);
        }
    }
    for (size_t i = 0; i < pairing_cell.size(); ++i) {
        if (!std::isfinite(pairing_cell[i])) {
            pairing_cell.clear();
            return false;
        }
    }
    return true;
}
