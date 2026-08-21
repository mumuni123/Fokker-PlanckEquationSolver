// Stage H1 acceptance for tail open boundaries and MPI migration
// (sections 6.6.3, 6.6.4, 6.6.5 and 15 H1):
//   * migration-consistency: the same trajectory set run with 1, 2 and 5
//     ranks must produce the same global particles, density, trajectory
//     current samples and left/right outflow ledgers (checked against the
//     analytic E=0 solution inside the binary; the rank-0 result file makes
//     the 1/2/5-rank runs directly comparable);
//   * open-truncation: a particle starting near a physical boundary moves
//     out step by step; in-domain density, boundary-shape ledgers, truncated
//     trajectory current and the final outflow ledger stay continuous, with
//     no periodic wrap and no in-domain shape renormalization;
//   * cfl-contract: c*dt < local rank width is enforced at startup and the
//     multi-hop migration guard is armed.
//
// Usage:
//   background_tail_open_boundary_mpi_test --case <case> [--result <path>]
// Run with yhrun -n 1 / -n 2 / -n 5 (different --result paths).
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "grid.h"
#include "maxwell.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "migration-consistency";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty();
}

struct Metrics {
    bool all_steps_ok;
    bool ids_unique;
    double max_pos_error_m;
    double max_uy_uz_error;
    double max_continuity_rel_l1;
    double max_continuity_rel_linf;
    double max_truncation_ledger_error;
    double max_number_balance_error;
    double final_density_l2_rel;
    double density_weight_balance_rel;
    double final_total_weight;
    double left_outflow_number;
    double right_outflow_number;
    double left_outflow_px;
    double right_outflow_px;
    double left_outflow_ke;
    double right_outflow_ke;
    double current_face_0;
    double current_face_50;
    double current_face_100;
    double current_face_150;
    double current_face_200;
    double cfl_ok_safe;
    double cfl_bad_rejected;
    double boundary_flux_closure_rel;
    double boundary_ledger_evolution_ok;
    Metrics()
        : all_steps_ok(false), ids_unique(false), max_pos_error_m(0.0),
          max_uy_uz_error(0.0), max_continuity_rel_l1(0.0),
          max_continuity_rel_linf(0.0), max_truncation_ledger_error(0.0),
          max_number_balance_error(0.0), final_density_l2_rel(0.0),
          density_weight_balance_rel(0.0), final_total_weight(0.0),
          left_outflow_number(0.0), right_outflow_number(0.0),
          left_outflow_px(0.0), right_outflow_px(0.0),
          left_outflow_ke(0.0), right_outflow_ke(0.0),
          current_face_0(0.0), current_face_50(0.0),
          current_face_100(0.0), current_face_150(0.0),
          current_face_200(0.0), cfl_ok_safe(0.0),
          cfl_bad_rejected(0.0), boundary_flux_closure_rel(0.0),
          boundary_ledger_evolution_ok(0.0)
    {}
};

struct TestParticle {
    double x_um;
    double ux;
    double uy;
    double uz;
    double weight;
};

double vx_from_u(const TestParticle& p)
{
    const double gamma = std::sqrt(
        1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
    return Const::c * p.ux / gamma;
}

void add_particle_on_owner(BackgroundTailPIC& tail, const SpatialGrid& grid,
                           int rank, int mpi_size, const TestParticle& tp)
{
    const double x = tp.x_um * Const::micro;
    const double x_left = grid.ix_start * grid.dx;
    const double x_right = (grid.ix_start + grid.nx_local) * grid.dx;
    if (x >= x_left && x < x_right) {
        BackgroundTailParticle p;
        p.x = x;
        p.ux = tp.ux;
        p.uy = tp.uy;
        p.uz = tp.uz;
        p.weight = tp.weight;
        p.id = tail.next_particle_id(rank);
        tail.particles.push_back(p);
    }
    (void)mpi_size;
}

Metrics run_migration_consistency(int rank, int mpi_size)
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(rank, mpi_size, 200, 2.0 * Const::micro);
    const double dt = 1.0e-16;
    const int steps = 60;
    const double length = grid.dx * static_cast<double>(grid.nx_global);

    const TestParticle particles[6] = {
        { 0.10,  4.0, 0.0, 0.0, 1.0e20 },   // crosses all ranks rightward
        { 1.80, -2.0, 0.0, 0.0, 1.1e20 },   // crosses all ranks leftward
        { 0.05, -4.0, 0.0, 0.0, 1.2e20 },   // exits left
        { 1.95,  4.0, 0.0, 0.0, 1.3e20 },   // exits right
        { 0.35,  1.0, 2.0, 1.5, 1.4e20 },   // survives, uy/uz preserved
        { 1.70, -1.0, 0.0, 0.0, 1.5e20 }    // survives
    };

    BackgroundTailPIC tail;
    tail.init(grid);
    for (int i = 0; i < 6; ++i) {
        add_particle_on_owner(tail, grid, rank, mpi_size, particles[i]);
    }
    tail.deposit_density(grid, rank, mpi_size);

    EMFields zero_fields;
    zero_fields.init(grid);
    const double x_left = grid.ix_start * grid.dx;
    const double x_right = (grid.ix_start + grid.nx_local) * grid.dx;

    bool steps_ok = true;
    double local_out[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    for (int s = 0; s < steps; ++s) {
        tail.begin_step(grid, dt);
        tail.drift_half(grid, 0.5 * dt, rank, mpi_size);
        tail.kick(grid, zero_fields, dt, rank, mpi_size);
        tail.drift_half(grid, 0.5 * dt, rank, mpi_size);
        tail.deposit_density(grid, rank, mpi_size);
        tail.finalize_trajectory_current(grid, dt, rank, mpi_size);
        const TailContinuityResult r =
            tail.audit_continuity(grid, dt, std::vector<double>(), rank,
                                  mpi_size);
        m.max_continuity_rel_l1 = std::max(m.max_continuity_rel_l1, r.rel_l1);
        m.max_continuity_rel_linf =
            std::max(m.max_continuity_rel_linf, r.rel_linf);
        m.max_truncation_ledger_error = std::max(
            m.max_truncation_ledger_error, r.truncation_ledger_error);
        m.max_number_balance_error = std::max(
            m.max_number_balance_error, r.number_balance_error);
        local_out[0] += tail.outflow_ledger().left_number;
        local_out[1] += tail.outflow_ledger().right_number;
        local_out[2] += tail.outflow_ledger().left_px;
        local_out[3] += tail.outflow_ledger().right_px;
        local_out[4] += tail.outflow_ledger().left_kinetic_energy;
        local_out[5] += tail.outflow_ledger().right_kinetic_energy;
        if (tail.migration_failed() || !tail.finite() ||
            !tail.nonnegative_weights()) {
            steps_ok = false;
        }
        // Ownership invariant: every local particle must be inside the local
        // domain (rank boundaries are not physical boundaries).
        for (size_t i = 0; i < tail.particles.size(); ++i) {
            const double x = tail.particles[i].x;
            if (x < x_left - 1.0e-12 || x >= x_right + 1.0e-12) {
                steps_ok = false;
            }
        }
    }

    // Analytic expectation for E = 0: the drift is exact.
    const double T = static_cast<double>(steps) * dt;
    double max_pos_error = 0.0;
    double max_uy_uz_error = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const TestParticle* tp = NULL;
        for (int k = 0; k < 6; ++k) {
            if (p.weight == particles[k].weight) {
                tp = &particles[k];
                break;
            }
        }
        if (tp == NULL) {
            steps_ok = false;
            continue;
        }
        const double x_analytic = tp->x_um * Const::micro +
                                  vx_from_u(*tp) * T;
        max_pos_error = std::max(
            max_pos_error, std::fabs(p.x - x_analytic));
        max_uy_uz_error = std::max(
            max_uy_uz_error,
            std::max(std::fabs(p.uy - tp->uy), std::fabs(p.uz - tp->uz)));
    }
    m.max_pos_error_m = max_pos_error;
    m.max_uy_uz_error = max_uy_uz_error;

    // Global aggregates.
    double local_weight = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        local_weight += tail.particles[i].weight;
    }
    MPI_Allreduce(&local_weight, &m.final_total_weight, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double global_out[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_out, global_out, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    m.left_outflow_number = global_out[0];
    m.right_outflow_number = global_out[1];
    m.left_outflow_px = global_out[2];
    m.right_outflow_px = global_out[3];
    m.left_outflow_ke = global_out[4];
    m.right_outflow_ke = global_out[5];

    // ID uniqueness: high bits must carry the creating rank and low bits
    // must be unique within the rank.  Migrated particles legitimately
    // reside on ranks different from their creator (section 6.4), so the
    // check only verifies local uniqueness and a valid creator rank.
    std::set<std::uint64_t> ids;
    bool ids_ok = true;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const std::uint64_t id = tail.particles[i].id;
        if (ids.count(id) > 0) ids_ok = false;
        ids.insert(id);
        if ((id >> 32) >= static_cast<std::uint64_t>(mpi_size)) {
            ids_ok = false;
        }
    }
    m.ids_unique = ids_ok;
    m.all_steps_ok = steps_ok;

    // Density L2 against the analytic final state (all six particles known;
    // only survivors contribute).
    double local_sq = 0.0;
    double ref_sq = 0.0;
    for (int i = 0; i < grid.nx_local; ++i) {
        double n_analytic = 0.0;
        for (int k = 0; k < 6; ++k) {
            const double x_analytic = particles[k].x_um * Const::micro +
                                      vx_from_u(particles[k]) * T;
            if (x_analytic < 0.0 || x_analytic >= length) continue;
            const CellDepositWeights cw =
                ParticleShape1D::cell_weights(x_analytic, grid);
            const int ig = grid.ix_start + i;
            double share = 0.0;
            if (cw.cell0 == ig) share = cw.w0;
            else if (cw.cell1 == ig) share = cw.w1;
            n_analytic += particles[k].weight * share * (1.0 / grid.dx);
        }
        const double n_actual = tail.density[static_cast<size_t>(i)];
        const double d = n_actual - n_analytic;
        local_sq += d * d;
        ref_sq += n_analytic * n_analytic;
    }
    double global_sq = 0.0;
    double global_ref_sq = 0.0;
    MPI_Allreduce(&local_sq, &global_sq, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&ref_sq, &global_ref_sq, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    m.final_density_l2_rel =
        std::sqrt(global_sq) / std::max(1.0, std::sqrt(global_ref_sq));

    const double total_weight = 7.5e20;
    m.density_weight_balance_rel =
        tail.density_weight_balance(grid) / total_weight;

    // Sample the finalized trajectory current on global faces.
    {
        const int local_count = static_cast<int>(tail.current_face_x.size());
        std::vector<int> counts(static_cast<size_t>(mpi_size), 0);
        std::vector<int> displs(static_cast<size_t>(mpi_size), 0);
        MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
                   MPI_COMM_WORLD);
        std::vector<double> global_current;
        if (rank == 0) {
            int total = 0;
            for (int r = 0; r < mpi_size; ++r) {
                displs[static_cast<size_t>(r)] = total;
                total += counts[static_cast<size_t>(r)];
            }
            global_current.assign(static_cast<size_t>(total), 0.0);
        }
        MPI_Gatherv(tail.current_face_x.data(), local_count, MPI_DOUBLE,
                    global_current.data(), counts.data(), displs.data(),
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            // Each rank contributes nx_local+1 faces and the shared face is
            // counted twice, so the gathered index of global face g is
            // g + owner_rank(g) (with uniform nx_local).
            const int face_samples[5] = { 0, 50, 100, 150, 200 };
            double* samples[5] = {
                &m.current_face_0, &m.current_face_50,
                &m.current_face_100, &m.current_face_150,
                &m.current_face_200
            };
            for (int s = 0; s < 5; ++s) {
                const int g = face_samples[s];
                const int owner = g / grid.nx_local;
                const int index = g + owner;
                if (owner < mpi_size &&
                    index < static_cast<int>(global_current.size())) {
                    *(samples[s]) = global_current[static_cast<size_t>(index)];
                }
            }
        }
    }
    MPI_Bcast(&m.current_face_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m.current_face_50, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m.current_face_100, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m.current_face_150, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m.current_face_200, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return m;
}

Metrics run_open_truncation(int rank, int mpi_size)
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(rank, mpi_size, 200, 2.0 * Const::micro);
    const double dt = 1.0e-16;
    const int steps = 10;
    const double length = grid.dx * static_cast<double>(grid.nx_global);

    BackgroundTailPIC tail;
    tail.init(grid);
    const double w_left = 1.0e20;
    const double w_right = 2.0e20;
    if (rank == 0) {
        BackgroundTailParticle p;
        p.x = 0.25 * grid.dx;
        p.ux = -0.05;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = w_left;
        p.id = tail.next_particle_id(rank);
        tail.particles.push_back(p);
    }
    if (rank == mpi_size - 1) {
        BackgroundTailParticle p;
        p.x = length - 0.25 * grid.dx;
        p.ux = 0.05;
        p.uy = 0.0;
        p.uz = 0.0;
        p.weight = w_right;
        p.id = tail.next_particle_id(rank);
        tail.particles.push_back(p);
    }
    tail.deposit_density(grid, rank, mpi_size);

    EMFields zero_fields;
    zero_fields.init(grid);
    const double initial_ledger_left =
        (rank == 0) ? ParticleShape1D::cell_weights(
                          0.25 * grid.dx, grid).w0 * w_left : 0.0;
    const double initial_ledger_right =
        (rank == mpi_size - 1)
            ? ParticleShape1D::cell_weights(
                  length - 0.25 * grid.dx, grid).w1 * w_right : 0.0;

    double flux_left = 0.0;   // accumulated -J_face0*dt/qe (number, left)
    double flux_right = 0.0;  // accumulated +J_faceLast*dt/qe (number, right)
    double trunc_left_total = 0.0;
    double trunc_right_total = 0.0;
    double max_ledger_evolution_error = 0.0;
    double local_left_out = 0.0;
    double local_right_out = 0.0;
    bool ok = true;
    for (int s = 0; s < steps; ++s) {
        tail.begin_step(grid, dt);
        tail.drift_half(grid, 0.5 * dt, rank, mpi_size);
        tail.kick(grid, zero_fields, dt, rank, mpi_size);
        tail.drift_half(grid, 0.5 * dt, rank, mpi_size);
        tail.deposit_density(grid, rank, mpi_size);
        tail.finalize_trajectory_current(grid, dt, rank, mpi_size);
        const TailContinuityResult r =
            tail.audit_continuity(grid, dt, std::vector<double>(), rank,
                                  mpi_size);
        m.max_continuity_rel_l1 = std::max(m.max_continuity_rel_l1, r.rel_l1);
        m.max_continuity_rel_linf =
            std::max(m.max_continuity_rel_linf, r.rel_linf);
        m.max_truncation_ledger_error = std::max(
            m.max_truncation_ledger_error, r.truncation_ledger_error);
        m.max_number_balance_error = std::max(
            m.max_number_balance_error, r.number_balance_error);
        if (tail.migration_failed() || !tail.finite() ||
            !tail.nonnegative_weights()) {
            ok = false;
        }
        // While the boundary particle is still inside, the deposit ledger
        // must equal the outside CIC share (no in-domain renormalization):
        // density_weight_balance == 0 every step.
        const double balance_scale = w_left + w_right;
        max_ledger_evolution_error = std::max(
            max_ledger_evolution_error,
            tail.density_weight_balance(grid) / balance_scale);
        // Accumulate the boundary-face trajectory flux in number units.
        flux_left += -(tail.current_face_x[0] * dt) / Const::qe;
        const size_t last_face =
            static_cast<size_t>(tail.current_face_x.size() - 1);
        flux_right += -(tail.current_face_x[last_face] * dt) / Const::qe;
        trunc_left_total += tail.truncation_shape_left();
        trunc_right_total += tail.truncation_shape_right();
        local_left_out += tail.outflow_ledger().left_number;
        local_right_out += tail.outflow_ledger().right_number;
    }
    m.boundary_ledger_evolution_ok =
        (max_ledger_evolution_error <= 1.0e-12) ? 1.0 : 0.0;

    // The boundary-face current integrated over the run must equal the
    // change of the CIC share below the boundary face.  At the left boundary
    // that share IS the outside-share ledger; at the right boundary it is
    // the in-domain complement (hence the opposite sign).
    double local_flux_closure = 0.0;
    if (rank == 0) {
        local_flux_closure = std::max(
            local_flux_closure,
            std::fabs(flux_left + (trunc_left_total - initial_ledger_left)) /
                w_left);
    }
    if (rank == mpi_size - 1) {
        local_flux_closure = std::max(
            local_flux_closure,
            std::fabs(flux_right - (trunc_right_total - initial_ledger_right)) /
                w_right);
    }
    MPI_Allreduce(&local_flux_closure, &m.boundary_flux_closure_rel, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Global outflow: the left particle exits left, the right particle
    // exits right; no periodic wrap and no leftovers.
    MPI_Allreduce(&local_left_out, &m.left_outflow_number, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_right_out, &m.right_outflow_number, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    double local_weight = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        local_weight += tail.particles[i].weight;
    }
    MPI_Allreduce(&local_weight, &m.final_total_weight, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    m.all_steps_ok = ok && m.left_outflow_number == w_left &&
                     m.right_outflow_number == w_right &&
                     m.final_total_weight == 0.0;
    return m;
}

Metrics run_cfl_contract(int rank, int mpi_size)
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(rank, mpi_size, 200, 2.0 * Const::micro);
    BackgroundTailPIC tail;
    tail.init(grid);
    const double rank_width = grid.dx * static_cast<double>(grid.nx_local);
    const double dt_safe = 0.5 * rank_width / Const::c;
    const double dt_bad = 1.1 * rank_width / Const::c;
    m.cfl_ok_safe = tail.cfl_contract_ok(grid, dt_safe, mpi_size) ? 1.0 : 0.0;
    m.cfl_bad_rejected =
        tail.cfl_contract_ok(grid, dt_bad, mpi_size) ? 0.0 : 1.0;
    m.all_steps_ok = (m.cfl_ok_safe == 1.0 && m.cfl_bad_rejected == 1.0);
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << std::setprecision(17);
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "all_steps_ok=" << (m.all_steps_ok ? 1 : 0) << "\n";
    out << "ids_unique=" << (m.ids_unique ? 1 : 0) << "\n";
    out << "max_pos_error_m=" << m.max_pos_error_m << "\n";
    out << "max_uy_uz_error=" << m.max_uy_uz_error << "\n";
    out << "max_continuity_rel_l1=" << m.max_continuity_rel_l1 << "\n";
    out << "max_continuity_rel_linf=" << m.max_continuity_rel_linf << "\n";
    out << "max_truncation_ledger_error=" << m.max_truncation_ledger_error
        << "\n";
    out << "max_number_balance_error=" << m.max_number_balance_error << "\n";
    out << "final_density_l2_rel=" << m.final_density_l2_rel << "\n";
    out << "density_weight_balance_rel=" << m.density_weight_balance_rel
        << "\n";
    out << "final_total_weight=" << m.final_total_weight << "\n";
    out << "left_outflow_number=" << m.left_outflow_number << "\n";
    out << "right_outflow_number=" << m.right_outflow_number << "\n";
    out << "left_outflow_px=" << m.left_outflow_px << "\n";
    out << "right_outflow_px=" << m.right_outflow_px << "\n";
    out << "left_outflow_ke=" << m.left_outflow_ke << "\n";
    out << "right_outflow_ke=" << m.right_outflow_ke << "\n";
    out << "current_face_0=" << m.current_face_0 << "\n";
    out << "current_face_50=" << m.current_face_50 << "\n";
    out << "current_face_100=" << m.current_face_100 << "\n";
    out << "current_face_150=" << m.current_face_150 << "\n";
    out << "current_face_200=" << m.current_face_200 << "\n";
    out << "cfl_ok_safe=" << m.cfl_ok_safe << "\n";
    out << "cfl_bad_rejected=" << m.cfl_bad_rejected << "\n";
    out << "boundary_flux_closure_rel=" << m.boundary_flux_closure_rel << "\n";
    out << "boundary_ledger_evolution_ok=" << m.boundary_ledger_evolution_ok
        << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: background_tail_open_boundary_mpi_test --case "
                     "<case> [--result <path>]\n"
                  << "cases: migration-consistency | open-truncation | "
                     "cfl-contract\n";
    }

    Metrics m;
    bool pass = false;
    if (ok && args.test_case == "migration-consistency") {
        m = run_migration_consistency(rank, size);
        pass = m.all_steps_ok && m.ids_unique &&
               m.max_pos_error_m <= 1.0e-12 &&
               m.max_uy_uz_error <= 1.0e-12 &&
               m.max_continuity_rel_l1 <= 1.0e-9 &&
               m.max_continuity_rel_linf <= 1.0e-9 &&
               m.max_truncation_ledger_error <= 1.0e-9 &&
               m.max_number_balance_error <= 1.0e-9 &&
               m.final_density_l2_rel <= 1.0e-9 &&
               m.density_weight_balance_rel <= 1.0e-9 &&
               m.final_total_weight == 5.0e20 &&
               m.left_outflow_number == 1.2e20 &&
               m.right_outflow_number == 1.3e20;
    } else if (ok && args.test_case == "open-truncation") {
        m = run_open_truncation(rank, size);
        pass = m.all_steps_ok &&
               m.max_continuity_rel_l1 <= 1.0e-9 &&
               m.max_continuity_rel_linf <= 1.0e-9 &&
               m.max_truncation_ledger_error <= 1.0e-9 &&
               m.max_number_balance_error <= 1.0e-9 &&
               m.boundary_flux_closure_rel <= 1.0e-9 &&
               m.boundary_ledger_evolution_ok == 1.0 &&
               m.left_outflow_number == 1.0e20 &&
               m.right_outflow_number == 2.0e20 &&
               m.final_total_weight == 0.0;
    } else if (ok && args.test_case == "cfl-contract") {
        m = run_cfl_contract(rank, size);
        pass = m.all_steps_ok;
    }

    int pass_all = pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &pass_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    pass = pass_all != 0;

    if (rank == 0) {
        if (!write_result_file(args.result_path, m, pass)) pass = false;
        std::cout << std::setprecision(17);
        std::cout << "all_steps_ok=" << (m.all_steps_ok ? 1 : 0)
                  << " ids_unique=" << (m.ids_unique ? 1 : 0)
                  << " max_pos_error_m=" << m.max_pos_error_m
                  << " max_uy_uz_error=" << m.max_uy_uz_error
                  << " max_continuity_rel_l1=" << m.max_continuity_rel_l1
                  << " max_continuity_rel_linf=" << m.max_continuity_rel_linf
                  << " max_truncation_ledger_error="
                  << m.max_truncation_ledger_error
                  << " max_number_balance_error=" << m.max_number_balance_error
                  << " final_density_l2_rel=" << m.final_density_l2_rel
                  << " density_weight_balance_rel="
                  << m.density_weight_balance_rel
                  << " final_total_weight=" << m.final_total_weight
                  << " left_outflow_number=" << m.left_outflow_number
                  << " right_outflow_number=" << m.right_outflow_number
                  << " left_outflow_px=" << m.left_outflow_px
                  << " right_outflow_px=" << m.right_outflow_px
                  << " left_outflow_ke=" << m.left_outflow_ke
                  << " right_outflow_ke=" << m.right_outflow_ke
                  << " current_face_0=" << m.current_face_0
                  << " current_face_50=" << m.current_face_50
                  << " current_face_100=" << m.current_face_100
                  << " current_face_150=" << m.current_face_150
                  << " current_face_200=" << m.current_face_200
                  << " cfl_ok_safe=" << m.cfl_ok_safe
                  << " cfl_bad_rejected=" << m.cfl_bad_rejected
                  << " boundary_flux_closure_rel="
                  << m.boundary_flux_closure_rel
                  << " boundary_ledger_evolution_ok="
                  << m.boundary_ledger_evolution_ok << "\n";
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
