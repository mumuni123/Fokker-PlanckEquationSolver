#include "open_electrostatic_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <mpi.h>

OpenElectrostaticSolver::OpenElectrostaticSolver()
{
    boundary_.type = ElectrostaticBoundaryType::LEFT_E;
    boundary_.e_left = 0.0;
    boundary_.phi_left = 0.0;
    boundary_.phi_right = 0.0;
    diagnostics_ = { 0.0, 0.0, 0.0, 0.0, 0.0 };
}

void OpenElectrostaticSolver::init(const SpatialGrid& grid,
                                   const ElectrostaticBoundary& boundary)
{
    grid_ = grid;
    boundary_ = boundary;
}

void OpenElectrostaticSolver::sync_face_interfaces(EMFields& fields,
                                                    int rank, int size) const
{
    if (size <= 1) return;
    const int left = rank > 0 ? rank - 1 : MPI_PROC_NULL;
    const int right = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
    double from_left = 0.0;
    const double to_right = fields.Ex_face[static_cast<size_t>(grid_.nx_local)];
    MPI_Sendrecv(&to_right, 1, MPI_DOUBLE, right, 7201,
                 &from_left, 1, MPI_DOUBLE, left, 7201,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (left != MPI_PROC_NULL) fields.Ex_face[0] = from_left;
}

void OpenElectrostaticSolver::exchange_cell_ghosts(
    std::vector<double>& values, int rank, int size,
    double left_value, double right_value) const
{
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    if (size > 1) {
        const int left = rank > 0 ? rank - 1 : MPI_PROC_NULL;
        const int right = rank + 1 < size ? rank + 1 : MPI_PROC_NULL;
        MPI_Sendrecv(values.data() + ng, ng, MPI_DOUBLE, left, 7202,
                     values.data() + ng + nxl, ng, MPI_DOUBLE, right, 7202,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(values.data() + ng + nxl - ng, ng, MPI_DOUBLE, right, 7203,
                     values.data(), ng, MPI_DOUBLE, left, 7203,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank == 0) {
        std::fill(values.begin(), values.begin() + ng, left_value);
    }
    if (rank == size - 1) {
        std::fill(values.begin() + ng + nxl, values.end(), right_value);
    }
}

void OpenElectrostaticSolver::solve(EMFields& fields, int mpi_rank, int mpi_size)
{
    solve(fields, mpi_rank, mpi_size, OpenGaussSolveOptions());
}

void OpenElectrostaticSolver::solve(EMFields& fields, int mpi_rank, int mpi_size,
                                    const OpenGaussSolveOptions& options)
{
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    double local_charge = 0.0;
    for (int ix = 0; ix < nxl; ++ix) local_charge += fields.rho[ng + ix] * grid_.dx;

    double prefix_charge = 0.0;
    diagnostics_.mpi_collective_seconds = 0.0;
    std::chrono::steady_clock::time_point collective_begin = std::chrono::steady_clock::now();
    MPI_Exscan(&local_charge, &prefix_charge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    diagnostics_.mpi_collective_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - collective_begin).count();
    if (mpi_rank == 0) prefix_charge = 0.0;

    double e_left = boundary_.type == ElectrostaticBoundaryType::LEFT_E
                  ? boundary_.e_left : 0.0;
    fields.Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    fields.Ex_face[0] = e_left + prefix_charge / Const::eps0;
    for (int ix = 0; ix < nxl; ++ix) {
        fields.Ex_face[static_cast<size_t>(ix + 1)] =
            fields.Ex_face[static_cast<size_t>(ix)] +
            fields.rho[ng + ix] * grid_.dx / Const::eps0;
    }

    if (boundary_.type == ElectrostaticBoundaryType::DIRICHLET_PHI) {
        double local_integral = 0.0;
        for (int ix = 0; ix < nxl; ++ix) {
            local_integral += 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix + 1]) * grid_.dx;
        }
        double integral = 0.0;
        collective_begin = std::chrono::steady_clock::now();
        MPI_Allreduce(&local_integral, &integral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        diagnostics_.mpi_collective_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - collective_begin).count();
        const double length = grid_.dx * static_cast<double>(grid_.nx_global);
        e_left = -(boundary_.phi_right - boundary_.phi_left + integral) / length;
        for (size_t iface = 0; iface < fields.Ex_face.size(); ++iface) fields.Ex_face[iface] += e_left;
    }

    sync_face_interfaces(fields, mpi_rank, mpi_size);
    for (int ix = 0; ix < nxl; ++ix) {
        fields.Ex[ng + ix] = 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix + 1]);
    }
    exchange_cell_ghosts(fields.Ex, mpi_rank, mpi_size,
                         fields.Ex_face.front(), fields.Ex_face.back());
    if (options.reconstruct_phi) reconstruct_phi(fields, mpi_rank, mpi_size);

    // Detailed audits are optional so the per-step production midpoint solve
    // does not pay for L1/Linf scans, endpoint gathering or phi rebuilding.
    const double not_computed = std::numeric_limits<double>::quiet_NaN();
    double local_l1 = 0.0;
    double local_linf = 0.0;
    if (options.compute_l1) {
        for (int ix = 0; ix < nxl; ++ix) {
            const double residual =
                (fields.Ex_face[ix + 1] - fields.Ex_face[ix]) / grid_.dx -
                fields.rho[ng + ix] / Const::eps0;
            local_l1 += std::fabs(residual) * grid_.dx;
            local_linf = std::max(local_linf, std::fabs(residual));
        }
    }
    if (options.compute_l1 || options.compute_boundary_audit) {
        // One structured SUM reduction carries both the L1 integral and the
        // net charge; a separate MAX reduction carries Linf.  This avoids the
        // previous pattern of three consecutive scalar reductions.
        double totals[2] = { options.compute_l1 ? local_l1 : 0.0, local_charge };
        double globals[2] = { 0.0, 0.0 };
        collective_begin = std::chrono::steady_clock::now();
        MPI_Allreduce(totals, globals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (options.compute_l1) {
            double global_linf = 0.0;
            MPI_Allreduce(&local_linf, &global_linf, 1, MPI_DOUBLE, MPI_MAX,
                          MPI_COMM_WORLD);
            diagnostics_.residual_l1 = globals[0];
            diagnostics_.residual_linf = global_linf;
            fields.last_gauss_residual_l1 = globals[0];
            fields.last_gauss_residual_linf = global_linf;
        } else {
            diagnostics_.residual_l1 = not_computed;
            diagnostics_.residual_linf = not_computed;
        }
        diagnostics_.mpi_collective_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - collective_begin).count();
        diagnostics_.net_charge = globals[1];
    } else {
        diagnostics_.residual_l1 = not_computed;
        diagnostics_.residual_linf = not_computed;
        diagnostics_.net_charge = not_computed;
    }

    if (options.compute_boundary_audit) {
        double end_values[2] = { 0.0, 0.0 };
        if (mpi_rank == 0) end_values[0] = fields.Ex_face.front();
        if (mpi_rank == mpi_size - 1) end_values[1] = fields.Ex_face.back();
        collective_begin = std::chrono::steady_clock::now();
        MPI_Allreduce(MPI_IN_PLACE, end_values, 2, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        diagnostics_.mpi_collective_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - collective_begin).count();
        diagnostics_.boundary_charge_residual =
            Const::eps0 * (end_values[1] - end_values[0]) -
            diagnostics_.net_charge;
    } else {
        diagnostics_.boundary_charge_residual = not_computed;
    }
}

void OpenElectrostaticSolver::reconstruct_phi(EMFields& fields,
                                              int mpi_rank, int mpi_size)
{
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    double local_drop = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        local_drop -= 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix + 1]) * grid_.dx;
    }
    double prefix_drop = 0.0;
    MPI_Exscan(&local_drop, &prefix_drop, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) prefix_drop = 0.0;
    double phi_node = (boundary_.type == ElectrostaticBoundaryType::DIRICHLET_PHI)
                    ? boundary_.phi_left : 0.0;
    phi_node += prefix_drop;
    for (int ix = 0; ix < nxl; ++ix) {
        fields.phi[ng + ix] = phi_node - 0.25 *
            (fields.Ex_face[ix] + fields.Ex_face[ix + 1]) * grid_.dx;
        phi_node -= 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix + 1]) * grid_.dx;
    }
    exchange_cell_ghosts(fields.phi, mpi_rank, mpi_size,
                         fields.phi[ng], fields.phi[ng + nxl - 1]);
    if (mpi_rank == 0) {
        for (int g = 0; g < ng; ++g) {
            fields.phi[ng - 1 - g] =
                ((boundary_.type == ElectrostaticBoundaryType::DIRICHLET_PHI)
                     ? boundary_.phi_left : fields.phi[ng]) +
                (static_cast<double>(g) + 0.5) * grid_.dx * fields.Ex_face.front();
        }
    }
    if (mpi_rank == mpi_size - 1) {
        const double phi_right = phi_node;
        for (int g = 0; g < ng; ++g) {
            fields.phi[ng + nxl + g] = phi_right -
                (static_cast<double>(g) + 0.5) * grid_.dx * fields.Ex_face.back();
        }
    }
}

double OpenElectrostaticSolver::boundary_power(double left_current,
                                                double right_current) const
{
    if (boundary_.type != ElectrostaticBoundaryType::DIRICHLET_PHI) return 0.0;
    return boundary_.phi_left * left_current - boundary_.phi_right * right_current;
}

double OpenElectrostaticSolver::boundary_energy_work(
    const EMFields& before, const EMFields& after,
    int mpi_rank, int mpi_size) const
{
    if (boundary_.type != ElectrostaticBoundaryType::DIRICHLET_PHI) {
        return 0.0;
    }
    if (boundary_.phi_left == 0.0 && boundary_.phi_right == 0.0) {
        return 0.0;
    }

    double local = 0.0;
    if (mpi_rank == 0 && !before.Ex_face.empty() && !after.Ex_face.empty()) {
        local += boundary_.phi_left * Const::eps0 *
            (after.Ex_face.front() - before.Ex_face.front());
    }
    if (mpi_rank == mpi_size - 1 &&
        !before.Ex_face.empty() && !after.Ex_face.empty()) {
        local -= boundary_.phi_right * Const::eps0 *
            (after.Ex_face.back() - before.Ex_face.back());
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

OpenPoissonWorkIdentity OpenElectrostaticSolver::evaluate_work_identity(
    const EMFields& before, const EMFields& after,
    const std::vector<double>& rho_delta,
    int mpi_rank, int mpi_size) const
{
    OpenPoissonWorkIdentity result;
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    if (before.Ex_face.size() != static_cast<size_t>(nxl + 1) ||
        after.Ex_face.size() != static_cast<size_t>(nxl + 1) ||
        before.phi.size() < static_cast<size_t>(ng + nxl) ||
        after.phi.size() < static_cast<size_t>(ng + nxl) ||
        rho_delta.size() < static_cast<size_t>(ng + nxl)) {
        return result;
    }

    double local[3] = { 0.0, 0.0, 0.0 };
    for (int ix = 0; ix < nxl; ++ix) {
        const double old_left = before.Ex_face[static_cast<size_t>(ix)];
        const double old_right = before.Ex_face[static_cast<size_t>(ix + 1)];
        const double new_left = after.Ex_face[static_cast<size_t>(ix)];
        const double new_right = after.Ex_face[static_cast<size_t>(ix + 1)];
        local[0] += Const::eps0 * grid_.dx *
            (old_left * old_left + old_left * old_right +
             old_right * old_right) / 6.0;
        local[1] += Const::eps0 * grid_.dx *
            (new_left * new_left + new_left * new_right +
             new_right * new_right) / 6.0;
        // rho is a finite-volume cell average.  The production reconstruction
        // stores phi at the cell midpoint, while E is linear inside a cell and
        // the field energy above integrates E^2 exactly.  The matching dual
        // quantity is therefore the cell-average potential, not phi(x_i).
        // For linear E and quadratic phi,
        //   <phi>_i = phi(x_i) + dx * (E_R - E_L) / 12.
        const double old_phi_average =
            before.phi[static_cast<size_t>(ng + ix)] +
            grid_.dx * (old_right - old_left) / 12.0;
        const double new_phi_average =
            after.phi[static_cast<size_t>(ng + ix)] +
            grid_.dx * (new_right - new_left) / 12.0;
        local[2] += 0.5 * (old_phi_average + new_phi_average) *
            rho_delta[static_cast<size_t>(ng + ix)] * grid_.dx;
    }
    double global[3] = { 0.0, 0.0, 0.0 };
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    result.field_energy_before = global[0];
    result.field_energy_after = global[1];
    result.field_energy_change = global[1] - global[0];
    result.electrode_work = boundary_energy_work(before, after, mpi_rank,
                                                  mpi_size);
    result.potential_charge_work = global[2];
    result.residual = result.field_energy_change - result.electrode_work -
        result.potential_charge_work;
    result.scale = std::max(std::numeric_limits<double>::min(),
        std::max(std::fabs(result.field_energy_change),
        std::max(std::fabs(result.electrode_work),
                 std::fabs(result.potential_charge_work))));
    result.finite = std::isfinite(result.field_energy_before) &&
        std::isfinite(result.field_energy_after) &&
        std::isfinite(result.field_energy_change) &&
        std::isfinite(result.electrode_work) &&
        std::isfinite(result.potential_charge_work) &&
        std::isfinite(result.residual);
    return result;
}

bool OpenElectrostaticSolver::build_potential_pairing_field(
    const EMFields& before, const EMFields& after,
    std::vector<double>& pairing_face,
    int mpi_rank, int mpi_size) const
{
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    if (nxl <= 0 || before.Ex_face.size() != static_cast<size_t>(nxl + 1) ||
        after.Ex_face.size() != static_cast<size_t>(nxl + 1) ||
        before.phi.size() < static_cast<size_t>(ng + nxl) ||
        after.phi.size() < static_cast<size_t>(ng + nxl)) {
        pairing_face.clear();
        return false;
    }

    // This is the same finite-volume cell-average potential used in
    // evaluate_work_identity(), time-centered between the two Poisson states.
    std::vector<double> paired_phi(static_cast<size_t>(nxl), 0.0);
    for (int ix = 0; ix < nxl; ++ix) {
        const double old_average =
            before.phi[static_cast<size_t>(ng + ix)] + grid_.dx *
            (before.Ex_face[static_cast<size_t>(ix + 1)] -
             before.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        const double new_average =
            after.phi[static_cast<size_t>(ng + ix)] + grid_.dx *
            (after.Ex_face[static_cast<size_t>(ix + 1)] -
             after.Ex_face[static_cast<size_t>(ix)]) / 12.0;
        paired_phi[static_cast<size_t>(ix)] = 0.5 * (old_average + new_average);
    }

    double left_neighbor = 0.0;
    double right_neighbor = 0.0;
    const int left = mpi_rank > 0 ? mpi_rank - 1 : MPI_PROC_NULL;
    const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : MPI_PROC_NULL;
    MPI_Sendrecv(&paired_phi.front(), 1, MPI_DOUBLE, left, 7611,
                 &right_neighbor, 1, MPI_DOUBLE, right, 7611,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&paired_phi.back(), 1, MPI_DOUBLE, right, 7612,
                 &left_neighbor, 1, MPI_DOUBLE, left, 7612,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    pairing_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
    // With endpoint face weights dx/2, these endpoint values make
    //   -dt <E_pair,J> = sum_i phi_bar_i Delta(rho_i) dx
    // after applying the finite-volume continuity equation.
    pairing_face[0] = mpi_rank == 0
        ? -2.0 * paired_phi.front() / grid_.dx
        : (left_neighbor - paired_phi.front()) / grid_.dx;
    for (int f = 1; f < nxl; ++f) {
        pairing_face[static_cast<size_t>(f)] =
            (paired_phi[static_cast<size_t>(f - 1)] -
             paired_phi[static_cast<size_t>(f)]) / grid_.dx;
    }
    pairing_face[static_cast<size_t>(nxl)] = mpi_rank + 1 == mpi_size
        ? 2.0 * paired_phi.back() / grid_.dx
        : (paired_phi.back() - right_neighbor) / grid_.dx;

    for (size_t i = 0; i < pairing_face.size(); ++i) {
        if (!std::isfinite(pairing_face[i])) {
            pairing_face.clear();
            return false;
        }
    }
    return true;
}

// JC2 (section 5.5.2): read-only face-to-cell helper.
void OpenElectrostaticSolver::populate_electric_components_from_faces(
    EMFields& fields, int mpi_rank, int mpi_size) const
{
    sync_face_interfaces(fields, mpi_rank, mpi_size);
    const int ng = grid_.nghost;
    const int nxl = grid_.nx_local;
    for (int ix = 0; ix < nxl; ++ix) {
        fields.Ex[ng + ix] = 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix + 1]);
    }
    exchange_cell_ghosts(fields.Ex, mpi_rank, mpi_size,
                         fields.Ex_face.front(), fields.Ex_face.back());
}
