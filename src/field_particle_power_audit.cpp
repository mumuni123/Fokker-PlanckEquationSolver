#include "field_particle_power_audit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {
// Face quadrature weight of section 4.5 item 4: interior faces use the full
// cell width; the two open-domain endpoints use half the width (trapezoidal
// endpoint rule, never the interior half-weight rule).
double face_weight(const SpatialGrid& grid, int local_face)
{
    const int global_face = grid.ix_start + local_face;
    // A shared face is owned by the rank on its right.  Therefore a rank's
    // local right face is excluded unless it is the global physical endpoint.
    if (local_face == grid.nx_local && global_face < grid.nx_global) {
        return 0.0;
    }
    if (global_face == 0 || global_face == grid.nx_global) {
        return 0.5 * grid.dx;
    }
    return grid.dx;
}
} // namespace

FieldParticlePowerAudit::FieldParticlePowerAudit()
    : left_region_end_(0), right_region_begin_(0)
{}

void FieldParticlePowerAudit::init(const SpatialGrid& grid)
{
    grid_ = grid;
    const double n = static_cast<double>(grid.nx_global);
    left_region_end_ = static_cast<int>(std::ceil(0.05 * n));
    right_region_begin_ = static_cast<int>(std::floor(0.95 * n));
    if (left_region_end_ < 0) left_region_end_ = 0;
    if (right_region_begin_ <= left_region_end_) {
        right_region_begin_ = std::max(left_region_end_ + 1, grid.nx_global - 1);
    }
    if (right_region_begin_ > grid.nx_global) right_region_begin_ = grid.nx_global;
}

int FieldParticlePowerAudit::region_id(int local_cell) const
{
    const int ig = grid_.global_cell(grid_.nghost + local_cell);
    if (ig < left_region_end_) return 0;
    if (ig >= right_region_begin_) return 2;
    return 1;
}

int FieldParticlePowerAudit::region_id_for_face(int local_face) const
{
    const int ig = grid_.ix_start + local_face;
    if (ig < left_region_end_) return 0;
    if (ig >= right_region_begin_) return 2;
    return 1;
}

double FieldParticlePowerAudit::face_quadrature_weight(int local_face) const
{
    if (local_face < 0 || local_face > grid_.nx_local) return 0.0;
    return face_weight(grid_, local_face);
}

bool FieldParticlePowerAudit::manufactured_gstar_identity(
    const std::vector<double>& e_face,
    const std::vector<double>& j_cell,
    double& face_work, double& cell_work, double& absolute_error) const
{
    return manufactured_gstar_identity(
        e_face, j_cell, j_cell.empty() ? 0.0 : j_cell.front(),
        j_cell.empty() ? 0.0 : j_cell.back(), 0, 1,
        face_work, cell_work, absolute_error);
}

bool FieldParticlePowerAudit::manufactured_gstar_identity(
    const std::vector<double>& e_face,
    const std::vector<double>& j_cell,
    double left_neighbor_cell, double right_neighbor_cell,
    int mpi_rank, int mpi_size,
    double& face_work, double& cell_work, double& absolute_error) const
{
    const size_t nf = static_cast<size_t>(grid_.nx_local) + 1;
    const size_t nc = static_cast<size_t>(grid_.nx_local);
    face_work = 0.0;
    cell_work = 0.0;
    absolute_error = std::numeric_limits<double>::infinity();
    if (e_face.size() != nf || j_cell.size() != nc || nc == 0) return false;

    if (mpi_size <= 0 || mpi_rank < 0 || mpi_rank >= mpi_size) return false;

    // Shared-face neighbour values are supplied by the MPI test/integrator.
    // A rank's right shared face has zero quadrature weight, so the right rank
    // is the unique owner of that face in the global sum.
    long double lhs = 0.0L;
    long double rhs = 0.0L;
    for (size_t f = 0; f < nf; ++f) {
        const double gstar = f == 0
            ? (mpi_rank > 0
                ? 0.5 * (left_neighbor_cell + j_cell.front())
                : j_cell.front())
            : (f + 1 == nf
                ? (mpi_rank + 1 < mpi_size
                    ? 0.5 * (j_cell.back() + right_neighbor_cell)
                    : j_cell.back())
                : 0.5 * (j_cell[f - 1] + j_cell[f]));
        lhs += static_cast<long double>(e_face[f]) *
               static_cast<long double>(gstar) *
               static_cast<long double>(face_weight(
                   grid_, static_cast<int>(f)));
    }
    for (size_t i = 0; i < nc; ++i) {
        const double e_cell = 0.5 * (e_face[i] + e_face[i + 1]);
        rhs += static_cast<long double>(e_cell) *
               static_cast<long double>(j_cell[i]) *
               static_cast<long double>(grid_.dx);
    }
    face_work = static_cast<double>(lhs);
    cell_work = static_cast<double>(rhs);
    absolute_error = std::fabs(face_work - cell_work);
    return std::isfinite(face_work) && std::isfinite(cell_work) &&
           std::isfinite(absolute_error);
}

double FieldParticlePowerAudit::roundoff_tolerance(double scale,
                                                   double unit_floor)
{
    return std::max(1.0e-12 * unit_floor,
                    512.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, scale));
}

void FieldParticlePowerAudit::swept_number_to_current(
    const std::vector<double>& swept, double dt, double charge,
    std::vector<double>& current) const
{
    current.assign(swept.size(), 0.0);
    const double inv_dt = dt > 0.0 ? 1.0 / dt : 0.0;
    for (size_t f = 0; f < swept.size(); ++f) {
        current[f] = charge * swept[f] * inv_dt;
    }
}

double FieldParticlePowerAudit::face_inner_product(
    const std::vector<double>& e_face,
    const std::vector<double>& j_face) const
{
    const size_t n = std::min(e_face.size(), j_face.size());
    double sum = 0.0;
    for (size_t f = 0; f < n; ++f) {
        sum += e_face[f] * j_face[f] * face_weight(grid_, static_cast<int>(f));
    }
    return sum;
}

double FieldParticlePowerAudit::evaluate_species_continuity(
    const std::vector<double>& number_n,
    const std::vector<double>& number_np1,
    const std::vector<double>& swept_face,
    const std::vector<double>& source_number_cell,
    double& residual_linf) const
{
    const int nxl = grid_.nx_local;
    residual_linf = 0.0;
    if (number_n.size() != static_cast<size_t>(nxl) ||
        number_np1.size() != static_cast<size_t>(nxl) ||
        swept_face.size() != static_cast<size_t>(nxl) + 1) {
        residual_linf = std::numeric_limits<double>::infinity();
        return residual_linf;
    }
    // Section 4.7.2 residual form: R_i = M_i^{n+1} - M_i^n
    //   - Q_{i-1/2} + Q_{i+1/2}.  The physical-boundary flux is already part
    // of swept_face[0] / swept_face[nxl].  Explicit injection, conversion and
    // return are represented by their actual per-cell source arrays.
    double linv = 0.0;
    for (int i = 0; i < nxl; ++i) {
        const double s = source_number_cell.size() == static_cast<size_t>(nxl)
            ? source_number_cell[static_cast<size_t>(i)] : 0.0;
        const double r = number_np1[static_cast<size_t>(i)] -
                         number_n[static_cast<size_t>(i)] -
                         swept_face[static_cast<size_t>(i)] +
                         swept_face[static_cast<size_t>(i) + 1] - s;
        linv = std::max(linv, std::fabs(r));
    }
    residual_linf = linv;
    return linv;
}

double FieldParticlePowerAudit::poisson_transport_residual_candidate(
    double field_energy_change, double electrode_work, double dt,
    const std::vector<double>& e_pair,
    const std::vector<double>& j_charge) const
{
    // R_{P<->J} = Delta U_E - W_electrode + dt <E_pair, J_charge>.
    return field_energy_change - electrode_work +
           dt * face_inner_product(e_pair, j_charge);
}

void FieldParticlePowerAudit::build_e_pair_candidates(
    const FieldParticlePowerAuditWorkspace& ws,
    std::vector<double> out[3]) const
{
    const int nxl = grid_.nx_local;
    const size_t nf = static_cast<size_t>(nxl) + 1;
    out[0].assign(nf, 0.0); // endpoint average
    out[1].assign(nf, 0.0); // production midpoint
    out[2].assign(nf, 0.0); // exact dual from potential-charge work
    for (size_t f = 0; f < nf; ++f) {
        const double en = f < ws.field_n_ex_face.size() ? ws.field_n_ex_face[f] : 0.0;
        const double em = f < ws.field_mid_ex_face.size() ? ws.field_mid_ex_face[f] : 0.0;
        const double ep = f < ws.field_np1_ex_face.size() ? ws.field_np1_ex_face[f] : 0.0;
        out[0][f] = 0.5 * (en + ep);
        out[1][f] = em;
        out[2][f] = f < ws.potential_pair_ex_face.size()
            ? ws.potential_pair_ex_face[f]
            : std::numeric_limits<double>::quiet_NaN();
    }
    (void)0;
}

int FieldParticlePowerAudit::classify_root_cause_mask(
    const FieldParticlePowerAuditResult& r, double full_abs_sum) const
{
    // Section 5 root-cause gates use per-term explanation fractions
    // F_c = sum|R_c| / sum|R_fp|.  Thresholds are applied by the Python
    // analyzer across coarse/fine; the mask here is the per-step candidate
    // set for the analyzer's uniqueness decision (section 5.1).
    (void)full_abs_sum;
    int mask = PAIRING_CAUSE_NONE;
    if (r.poisson_transport_residual != 0.0) mask |= PAIRING_CAUSE_TRANSPORT;
    if (r.current_pair_residual != 0.0) mask |= PAIRING_CAUSE_WORK;
    if (r.time_center_fraction != 0.0) mask |= PAIRING_CAUSE_TIME;
    if (r.force_work_residual != 0.0) mask |= PAIRING_CAUSE_PIC;
    if (r.conversion_transfer_residual != 0.0) mask |= PAIRING_CAUSE_CONVERSION;
    if (r.boundary_residual != 0.0) mask |= PAIRING_CAUSE_BOUNDARY;
    return mask;
}

FieldParticlePowerAuditResult FieldParticlePowerAudit::finalize(
    const FieldParticlePowerAuditWorkspace& ws, double dt,
    double electron_charge) const
{
    FieldParticlePowerAuditResult r;
    r.valid = false;
    const int nxl = grid_.nx_local;
    const size_t nf = static_cast<size_t>(nxl) + 1;

    // 1. verify_workspace_shapes().
    const bool shapes_ok =
        ws.bulk_number_n.size() == static_cast<size_t>(nxl) &&
        ws.bulk_number_np1.size() == static_cast<size_t>(nxl) &&
        ws.bulk_x1.bulk_number_swept_face.size() == nf &&
        ws.bulk_x2.bulk_number_swept_face.size() == nf &&
        ws.field_n_ex_face.size() == nf &&
        ws.field_mid_ex_face.size() == nf &&
        ws.field_np1_ex_face.size() == nf &&
        ws.potential_pair_ex_face.size() == nf;
    if (!shapes_ok) {
        r.all_finite = false;
        return r;
    }

    // 2. combine_bulk_x_half_swept_numbers().
    std::vector<double> bulk_swept(nf, 0.0);
    for (size_t f = 0; f < nf; ++f) {
        bulk_swept[f] = ws.bulk_x1.bulk_number_swept_face[f] +
                        ws.bulk_x2.bulk_number_swept_face[f];
    }

    // 3/4. build species charge current from raw swept numbers and the
    // production trajectory-current snapshots (PIC current is already A/m^2).
    std::vector<double> j_bulk(nf, 0.0);
    swept_number_to_current(bulk_swept, dt, electron_charge, j_bulk);
    const std::vector<double>& j_tail =
        ws.tail_trajectory.after_second_drift_current_face;
    const std::vector<double>& j_beam =
        ws.beam_trajectory.after_second_drift_current_face;
    std::vector<double> j_charge(nf, 0.0);
    for (size_t f = 0; f < nf; ++f) {
        j_charge[f] = j_bulk[f];
        if (f < j_tail.size()) j_charge[f] += j_tail[f];
        if (f < j_beam.size()) j_charge[f] += j_beam[f];
    }

    // 5. evaluate_species_continuity_with_boundary_sources().
    const std::vector<double>& tail_n = ws.tail_number_n;
    const std::vector<double>& tail_np1 = ws.tail_number_np1;
    const std::vector<double>& beam_n = ws.beam_number_n;
    const std::vector<double>& beam_np1 = ws.beam_number_np1;
    // Tail/Beam swept number is reconstructed from the final trajectory
    // current: Q = J * dt / charge (charge already inside J).
    std::vector<double> tail_swept(nf, 0.0);
    std::vector<double> beam_swept(nf, 0.0);
    if (!j_tail.empty()) {
        for (size_t f = 0; f < nf && f < j_tail.size(); ++f) {
            tail_swept[f] = j_tail[f] * dt / electron_charge;
        }
    }
    if (!j_beam.empty()) {
        for (size_t f = 0; f < nf && f < j_beam.size(); ++f) {
            beam_swept[f] = j_beam[f] * dt / electron_charge;
        }
    }
    r.continuity_linf_bulk = 0.0;
    r.continuity_linf_tail = 0.0;
    r.continuity_linf_beam = 0.0;
    double scale_sum = 0.0;
    const auto continuity_scale_for = [](const std::vector<double>& n0,
                                         const std::vector<double>& n1,
                                         const std::vector<double>& swept,
                                         const std::vector<double>& source) {
        double scale = 0.0;
        for (size_t i = 0; i < n0.size() && i < n1.size(); ++i) {
            double cell = std::fabs(n0[i]) + std::fabs(n1[i]);
            if (i < swept.size()) cell += std::fabs(swept[i]);
            if (i + 1 < swept.size()) cell += std::fabs(swept[i + 1]);
            if (i < source.size()) cell += std::fabs(source[i]);
            scale = std::max(scale, cell);
        }
        return scale;
    };
    r.continuity_linf_bulk =
        std::fabs(evaluate_species_continuity(
            ws.bulk_number_n, ws.bulk_number_np1, bulk_swept,
            ws.bulk_source_number_cell, r.continuity_linf_bulk));
    scale_sum = std::max(scale_sum, continuity_scale_for(
        ws.bulk_number_n, ws.bulk_number_np1, bulk_swept,
        ws.bulk_source_number_cell));
    if (!tail_n.empty() && tail_n.size() == static_cast<size_t>(nxl) &&
        tail_np1.size() == static_cast<size_t>(nxl)) {
        r.continuity_linf_tail =
            std::fabs(evaluate_species_continuity(
                tail_n, tail_np1, tail_swept,
                ws.tail_source_number_cell, r.continuity_linf_tail));
        scale_sum = std::max(scale_sum, continuity_scale_for(
            tail_n, tail_np1, tail_swept, ws.tail_source_number_cell));
    }
    if (!beam_n.empty() && beam_n.size() == static_cast<size_t>(nxl) &&
        beam_np1.size() == static_cast<size_t>(nxl)) {
        r.continuity_linf_beam =
            std::fabs(evaluate_species_continuity(
                beam_n, beam_np1, beam_swept,
                ws.beam_source_number_cell, r.continuity_linf_beam));
        scale_sum = std::max(scale_sum, continuity_scale_for(
            beam_n, beam_np1, beam_swept, ws.beam_source_number_cell));
    }
    r.continuity_linf = std::max(
        r.continuity_linf_bulk,
        std::max(r.continuity_linf_tail, r.continuity_linf_beam));
    r.continuity_scale = scale_sum;
    r.continuity_pass = std::isfinite(r.continuity_linf) &&
        r.continuity_linf <= roundoff_tolerance(scale_sum, scale_sum);

    // Per-component first-bad-cell continuity composition (section I2).  The
    // four source fields are signed "positive = added to the component"; the
    // net source equals injection + outflow + conversion + return and is the
    // same quantity already used by evaluate_species_continuity.
    const double continuity_tol = roundoff_tolerance(scale_sum, scale_sum);
    const auto fill_breakdown =
        [&](const std::vector<double>& n0, const std::vector<double>& n1,
            const std::vector<double>& swept,
            const std::vector<double>& net_source,
            const std::vector<double>& injection,
            const std::vector<double>& outflow,
            const std::vector<double>& conversion,
            const std::vector<double>& retur,
            const std::vector<double>& other,
            ContinuityComponentBreakdown& bd, int& bad) {
        bad = -1;
        for (int i = 0; i < nxl; ++i) {
            if (static_cast<size_t>(i) >= n0.size() ||
                static_cast<size_t>(i) >= n1.size() ||
                static_cast<size_t>(i + 1) >= swept.size()) continue;
            const double s = net_source.size() == static_cast<size_t>(nxl)
                ? net_source[static_cast<size_t>(i)] : 0.0;
            const double delta = n1[static_cast<size_t>(i)] -
                                 n0[static_cast<size_t>(i)];
            const double left = swept[static_cast<size_t>(i)];
            const double right = swept[static_cast<size_t>(i + 1)];
            const double residual = delta + (right - left) - s;
            if (bad < 0 && std::fabs(residual) > continuity_tol) {
                bad = i;
                bd.delta_n = delta;
                bd.left_face_swept = left;
                bd.right_face_swept = right;
                bd.injection_source =
                    injection.size() == static_cast<size_t>(nxl)
                    ? injection[static_cast<size_t>(i)] : 0.0;
                bd.outflow_source =
                    outflow.size() == static_cast<size_t>(nxl)
                    ? outflow[static_cast<size_t>(i)] : 0.0;
                bd.conversion_source =
                    conversion.size() == static_cast<size_t>(nxl)
                    ? conversion[static_cast<size_t>(i)] : 0.0;
                bd.return_source =
                    retur.size() == static_cast<size_t>(nxl)
                    ? retur[static_cast<size_t>(i)] : 0.0;
                bd.other_source =
                    other.size() == static_cast<size_t>(nxl)
                    ? other[static_cast<size_t>(i)] : 0.0;
                bd.residual = residual;
                bd.bad_global_index = grid_.ix_start + i;
            }
        }
    };
    const std::vector<double> empty_source;
    int bad_bulk = -1;
    int bad_tail = -1;
    int bad_beam = -1;
    fill_breakdown(ws.bulk_number_n, ws.bulk_number_np1, bulk_swept,
                   ws.bulk_source_number_cell,
                   empty_source, empty_source,
                   ws.bulk_conversion_number_cell,
                   ws.bulk_return_number_cell,
                   ws.bulk_other_number_cell,
                   r.continuity_bulk, bad_bulk);
    if (!tail_n.empty() && tail_n.size() == static_cast<size_t>(nxl) &&
        tail_np1.size() == static_cast<size_t>(nxl)) {
        fill_breakdown(tail_n, tail_np1, tail_swept,
                       ws.tail_source_number_cell,
                       empty_source, ws.tail_outflow_number_cell,
                       ws.tail_conversion_number_cell,
                       ws.tail_return_number_cell,
                       empty_source,
                       r.continuity_tail, bad_tail);
    }
    if (!beam_n.empty() && beam_n.size() == static_cast<size_t>(nxl) &&
        beam_np1.size() == static_cast<size_t>(nxl)) {
        fill_breakdown(beam_n, beam_np1, beam_swept,
                       ws.beam_source_number_cell,
                       ws.beam_injection_number_cell, empty_source,
                       empty_source, empty_source, empty_source,
                       r.continuity_beam, bad_beam);
    }
    if (bad_bulk >= 0) r.first_bad_index = grid_.ix_start + bad_bulk;
    else if (bad_tail >= 0) r.first_bad_index = grid_.ix_start + bad_tail;
    else if (bad_beam >= 0) r.first_bad_index = grid_.ix_start + bad_beam;

    // 6. sum_local_force_work_by_species().  The rank-local sums are stored;
    // the work-sum-vs-ledger comparison is nonlinear across MPI and is
    // computed by the integrator after its single packed reduction.  The PIC
    // sums include the out-of-domain CIC shares (section I3) so they equal
    // the global kick ledger exactly.
    const double sum_bulk_work =
        std::accumulate(ws.cell_work.bulk_delta_ke_cell.begin(),
                        ws.cell_work.bulk_delta_ke_cell.end(), 0.0);
    const double sum_tail_work =
        std::accumulate(ws.cell_work.tail_delta_ke_cell.begin(),
                        ws.cell_work.tail_delta_ke_cell.end(), 0.0) +
        ws.cell_work.tail_delta_ke_boundary;
    const double sum_beam_work =
        std::accumulate(ws.cell_work.beam_delta_ke_cell.begin(),
                        ws.cell_work.beam_delta_ke_cell.end(), 0.0) +
        ws.cell_work.beam_delta_ke_boundary;
    r.bulk_cell_work_sum = sum_bulk_work;
    r.tail_cell_work_sum = sum_tail_work;
    r.beam_cell_work_sum = sum_beam_work;

    // 7. build E_pair candidates and Poisson-transport residuals.
    std::vector<double> e_pair[3];
    build_e_pair_candidates(ws, e_pair);
    r.poisson_transport_residual_endpoint =
        poisson_transport_residual_candidate(
            ws.field_energy_change, ws.electrode_work, dt, e_pair[0], j_charge);
    r.poisson_transport_residual_midpoint =
        poisson_transport_residual_candidate(
            ws.field_energy_change, ws.electrode_work, dt, e_pair[1], j_charge);
    // This candidate is additionally checked against the already validated
    // Gate-F scalar identity.  The face form and scalar potential form must
    // agree; neither may be replaced by the final Ex field.
    r.poisson_transport_residual_discrete_gradient =
        poisson_transport_residual_candidate(
            ws.field_energy_change, ws.electrode_work, dt, e_pair[2], j_charge);
    r.poisson_identity_crosscheck =
        (ws.field_energy_change - ws.electrode_work -
         ws.poisson_potential_charge_work) - ws.poisson_identity_residual;
    // Validated candidate (section 4.5): the discrete gradient is the
    // mathematically derived candidate; the manufactured test in the pairing
    // unit test confirms it before it is trusted.
    const std::vector<double>& primary_e_pair =
        e_pair[r.validated_e_pair_candidate];
    r.poisson_transport_residual =
        poisson_transport_residual_candidate(
            ws.field_energy_change, ws.electrode_work, dt,
            primary_e_pair, j_charge);

    // 8. Signed force-work identity.  The existing ledgers are the production
    // force-driven Delta K; the cell sums are the independently exposed local
    // accumulation of the same final kick/remap.  Keeping the sign is required
    // for the four-term reconstruction.
    r.force_work_residual_bulk = ws.bulk_work_ledger - sum_bulk_work;
    r.force_work_residual_tail = ws.tail_work_ledger - sum_tail_work;
    r.force_work_residual_beam = ws.beam_work_ledger - sum_beam_work;
    r.force_work_residual = r.force_work_residual_bulk +
                            r.force_work_residual_tail +
                            r.force_work_residual_beam;

    // 9. current-pair residual: total field work vs Poisson-transport pairing.
    // The total field work is the rank-local sum of the per-cell force work
    // (the work ledgers are the same quantity); the integrator separately
    // checks those cell sums against the existing global ledgers.
    const double total_field_work = sum_bulk_work + sum_tail_work +
                                    sum_beam_work;
    const double dt_pair_charge = dt * face_inner_product(primary_e_pair,
                                                          j_charge);
    r.current_pair_residual = total_field_work - dt_pair_charge;
    r.current_pair_linf = 0.0;
    int first_current_pair_bad = -1;
    for (int i = 0; i < nxl; ++i) {
        const double charge_work_cell = 0.5 * dt * grid_.dx *
            (primary_e_pair[static_cast<size_t>(i)] *
                 j_charge[static_cast<size_t>(i)] +
             primary_e_pair[static_cast<size_t>(i + 1)] *
                 j_charge[static_cast<size_t>(i + 1)]);
        double force_work_cell = 0.0;
        if (static_cast<size_t>(i) < ws.cell_work.bulk_delta_ke_cell.size())
            force_work_cell += ws.cell_work.bulk_delta_ke_cell[static_cast<size_t>(i)];
        if (static_cast<size_t>(i) < ws.cell_work.tail_delta_ke_cell.size())
            force_work_cell += ws.cell_work.tail_delta_ke_cell[static_cast<size_t>(i)];
        if (static_cast<size_t>(i) < ws.cell_work.beam_delta_ke_cell.size())
            force_work_cell += ws.cell_work.beam_delta_ke_cell[static_cast<size_t>(i)];
        const double local_mismatch =
            std::fabs(force_work_cell - charge_work_cell);
        r.current_pair_linf = std::max(r.current_pair_linf, local_mismatch);
        const double local_scale = std::max(
            1.0, std::fabs(force_work_cell) + std::fabs(charge_work_cell));
        if (first_current_pair_bad < 0 &&
            local_mismatch > roundoff_tolerance(local_scale, local_scale)) {
            first_current_pair_bad = i;
        }
    }
    if (r.first_bad_index < 0 && first_current_pair_bad >= 0)
        r.first_bad_index = grid_.ix_start + first_current_pair_bad;
    if (r.first_bad_index < 0 &&
        std::fabs(r.poisson_identity_crosscheck) > roundoff_tolerance(
            std::max(1.0, std::fabs(r.poisson_identity_crosscheck)), 1.0)) {
        r.first_bad_index = 0;
    }

    // 9b. §15.13.4 step 3: direct face dual audit.  Charge-conserving face
    // current J_charge_face is already built (step 3/4, A/m^2 on nx_local+1
    // faces).  The velocity-force cell current J_force_cell is reconstructed
    // from the per-cell force work so that
    //   dt <E_pair, G* J_force_face> = sum_i force_work_cell[i],
    // the discrete adjoint of the cell-average gather
    //   E_pair_cell[i] = 0.5 (E_pair[i] + E_pair[i+1]).
    // Then -dt <E_pair, J_charge_face - G* J_force_face> must reconstruct
    // current_pair_residual to roundoff if the two currents share the dual
    // relation.  Read-only: nothing here modifies the physical state.
    r.j_charge_face = j_charge;
    r.j_force_cell.assign(static_cast<size_t>(nxl), 0.0);
    r.gstar_j_force_face.assign(nf, 0.0);
    r.dual_face_residual.assign(nf, 0.0);
    r.dual_left5_integral = 0.0;
    r.dual_core90_integral = 0.0;
    r.dual_right5_integral = 0.0;
    {
        const double inv_dt = dt > 0.0 ? 1.0 / dt : 0.0;
        const double e_floor =
            grid_.dx * 1.0e-10;   // V/m guard against divide-by-tiny
        const auto cell_e_pair = [&](int i) {
            const double e = 0.5 * (primary_e_pair[static_cast<size_t>(i)] +
                                    primary_e_pair[static_cast<size_t>(i + 1)]);
            return std::fabs(e) > e_floor ? e : (e >= 0.0 ? e_floor : -e_floor);
        };
        for (int i = 0; i < nxl; ++i) {
            double force_work_cell = 0.0;
            if (static_cast<size_t>(i) < ws.cell_work.bulk_delta_ke_cell.size())
                force_work_cell += ws.cell_work.bulk_delta_ke_cell[static_cast<size_t>(i)];
            if (static_cast<size_t>(i) < ws.cell_work.tail_delta_ke_cell.size())
                force_work_cell += ws.cell_work.tail_delta_ke_cell[static_cast<size_t>(i)];
            if (static_cast<size_t>(i) < ws.cell_work.beam_delta_ke_cell.size())
                force_work_cell += ws.cell_work.beam_delta_ke_cell[static_cast<size_t>(i)];
            // J_force_cell[i] = force_work_cell[i] / (dt * dx * E_pair_cell[i]).
            r.j_force_cell[static_cast<size_t>(i)] =
                force_work_cell * inv_dt / (grid_.dx * cell_e_pair(i));
        }
        // G* (adjoint of the cell-average gather): the interior faces
        // average the two neighbouring cells.  The shared boundary faces use
        // the neighbouring rank's boundary cell (filled by the integrator via
        // Sendrecv); the global physical endpoints alias the boundary cell.
        r.gstar_j_force_face[0] = (grid_.ix_start == 0)
            ? r.j_force_cell[0]
            : 0.5 * (ws.gstar_left_neighbor_cell + r.j_force_cell[0]);
        r.gstar_j_force_face[nf - 1] =
            (grid_.ix_start + nxl == grid_.nx_global)
            ? r.j_force_cell[static_cast<size_t>(nxl - 1)]
            : 0.5 * (r.j_force_cell[static_cast<size_t>(nxl - 1)] +
                     ws.gstar_right_neighbor_cell);
        for (size_t f = 1; f + 1 < nf; ++f) {
            r.gstar_j_force_face[f] =
                0.5 * (r.j_force_cell[f - 1] + r.j_force_cell[f]);
        }
        for (size_t f = 0; f < nf; ++f) {
            r.dual_face_residual[f] = j_charge[f] - r.gstar_j_force_face[f];
            const double region = f < r.j_charge_face.size()
                ? static_cast<int>(region_id_for_face(f)) : 1;
            const double e_w = primary_e_pair[f] * face_weight(grid_, static_cast<int>(f));
            const double contrib = -dt * e_w * r.dual_face_residual[f];
            if (region == 0) r.dual_left5_integral += contrib;
            else if (region == 2) r.dual_right5_integral += contrib;
            else r.dual_core90_integral += contrib;
        }
    }

    // §15.13.6 item 2: separate the in-domain dual face work from the
    // out-of-domain Tail/Beam boundary force work.  The two terms must add to
    // the current-pair residual; neither term is used by the production step.
    r.dual_in_domain_work = r.dual_left5_integral +
                            r.dual_core90_integral +
                            r.dual_right5_integral;
    r.boundary_force_work = ws.cell_work.tail_delta_ke_boundary +
                            ws.cell_work.beam_delta_ke_boundary;
    r.dual_plus_boundary_work = r.dual_in_domain_work +
                                r.boundary_force_work;
    r.dual_reconstruction_error = r.dual_plus_boundary_work -
                                  r.current_pair_residual;
    r.dual_reconstruction_tolerance = roundoff_tolerance(
        std::max(1.0, std::max(std::fabs(r.dual_plus_boundary_work),
                               std::fabs(r.current_pair_residual))),
        1.0);
    r.dual_reconstruction_pass =
        std::isfinite(r.dual_reconstruction_error) &&
        std::fabs(r.dual_reconstruction_error) <=
            r.dual_reconstruction_tolerance;

    // 10. conversion and boundary terms.
    r.conversion_transfer_residual = ws.conversion_energy_removed;
    r.boundary_residual = ws.boundary_energy_source;

    // 11. reconstruct_full_field_particle_residual().  Both the reconstructed
    // sum and the actual full residual are stored rank-locally; the integrator
    // globalizes them and computes the nonlinear mismatch afterwards.
    r.reconstructed_full_residual =
        r.poisson_transport_residual + r.current_pair_residual +
        r.force_work_residual + r.conversion_transfer_residual +
        r.boundary_residual;
    const double total_force_delta_ke = ws.bulk_work_ledger +
                                        ws.tail_work_ledger +
                                        ws.beam_work_ledger;
    r.full_residual = ws.field_energy_change - ws.electrode_work +
                      total_force_delta_ke + ws.conversion_energy_removed +
                      ws.boundary_energy_source;

    // Complete the nonlinear acceptance fields for a one-rank/direct audit.
    // The production integrator recomputes these after its packed MPI
    // reductions, so these values are also valid for unit tests that invoke
    // the production calculator without an MPI wrapper.
    const double work_scale = std::max(
        1.0, std::fabs(ws.bulk_work_ledger) +
             std::fabs(ws.tail_work_ledger) +
             std::fabs(ws.beam_work_ledger) +
             std::fabs(sum_bulk_work) + std::fabs(sum_tail_work) +
             std::fabs(sum_beam_work));
    const double work_tol = roundoff_tolerance(work_scale, work_scale);
    r.local_work_sum_matches_existing_ledger =
        std::fabs(r.force_work_residual_bulk) <= work_tol &&
        std::fabs(r.force_work_residual_tail) <= work_tol &&
        std::fabs(r.force_work_residual_beam) <= work_tol;
    r.reconstruction_mismatch =
        std::fabs(r.reconstructed_full_residual - r.full_residual);
    r.reconstruction_scale = std::max(
        1.0, std::max(std::fabs(r.reconstructed_full_residual),
                      std::fabs(r.full_residual)));
    r.roundoff_tolerance = roundoff_tolerance(
        r.reconstruction_scale, r.reconstruction_scale);
    r.full_residual_reconstruction_pass =
        std::isfinite(r.reconstruction_mismatch) &&
        r.reconstruction_mismatch <= r.roundoff_tolerance;

    // 12. all-finite (before any SUM reduction, section 4.6.1 note 4).
    r.all_finite =
        std::isfinite(r.continuity_linf) &&
        std::isfinite(r.poisson_transport_residual) &&
        std::isfinite(r.force_work_residual) &&
        std::isfinite(r.current_pair_residual) &&
        std::isfinite(r.current_pair_linf) &&
        std::isfinite(r.poisson_identity_crosscheck) &&
        std::isfinite(r.conversion_transfer_residual) &&
        std::isfinite(r.boundary_residual) &&
        std::isfinite(r.reconstructed_full_residual) &&
        std::isfinite(r.full_residual) &&
        std::isfinite(r.bulk_cell_work_sum) &&
        std::isfinite(r.tail_cell_work_sum) &&
        std::isfinite(r.beam_cell_work_sum) &&
        std::isfinite(r.dual_in_domain_work) &&
        std::isfinite(r.boundary_force_work) &&
        std::isfinite(r.dual_plus_boundary_work) &&
        std::isfinite(r.dual_reconstruction_error);
    r.full_residual_reconstruction_pass =
        r.all_finite && r.full_residual_reconstruction_pass;

    // 13. A-F candidate fractions (section 5.1 item 5):
    //     F_c = sum|R_c| / sum|R_fp|.  Rank-local; the integrator re-derives
    //     these from the globalized residuals.
    const double full_abs = std::fabs(r.full_residual);
    const double denom = full_abs > 0.0 ? full_abs : 1.0;
    r.transport_fraction = std::fabs(r.poisson_transport_residual) / denom;
    r.work_current_fraction = std::fabs(r.current_pair_residual) / denom;
    const double time_center_diff =
        std::fabs(r.poisson_transport_residual_discrete_gradient -
                  r.poisson_transport_residual_midpoint);
    r.time_center_fraction = time_center_diff / denom;
    r.pic_fraction = std::fabs(r.force_work_residual) / denom;
    r.conversion_fraction = std::fabs(r.conversion_transfer_residual) / denom;
    r.boundary_fraction = std::fabs(r.boundary_residual) / denom;

    // 14. classify_regions_and_root_cause_mask().
    r.root_cause_mask = classify_root_cause_mask(r, full_abs);

    r.valid = r.all_finite;
    return r;
}
