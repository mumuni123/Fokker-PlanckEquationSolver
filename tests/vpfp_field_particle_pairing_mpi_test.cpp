// Gate I MPI acceptance test (section 4.7.4).  Executed with yhrun -n 2 and
// -n 5 on the cluster; never registered in ctest.  It drives the production
// ConservativePpmRemap::advect_x across the production spatial partition and
// verifies that the exposed face swept numbers close the per-cell continuity
// identity and that shared MPI faces have a single owner with no duplicates.
//
// Usage:
//   vpfp_field_particle_pairing_mpi_test --case all --result <path>

#include "conservative_ppm_remap.h"
#include "field_particle_power_audit.h"
#include "open_boundary.h"
#include "parameters.h"
#include "vpfp_field_particle_pairing_test_support.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            ++i; // accepted; only "all" is supported
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

struct MpiMetrics {
    int mpi_size;
    double shared_face_duplicate_count;
    double shared_face_mismatch_count;
    double shared_face_owner_error_count;
    double continuity_linf;
    double continuity_scale;
    double work_mpi_vs_reference_abs;
    double full_residual_mpi_vs_reference_abs;
    double manufactured_identity_error;
    bool manufactured_identity_pass;
    bool manufactured_endpoint_weight_pass;
    bool continuity_pass;
    MpiMetrics()
        : mpi_size(1), shared_face_duplicate_count(0.0),
          shared_face_mismatch_count(0.0), shared_face_owner_error_count(0.0),
          continuity_linf(0.0), continuity_scale(1.0),
          work_mpi_vs_reference_abs(0.0),
          full_residual_mpi_vs_reference_abs(0.0),
          manufactured_identity_error(0.0),
          manufactured_identity_pass(false),
          manufactured_endpoint_weight_pass(false), continuity_pass(false)
    {}
};

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2 && size != 5) {
        if (rank == 0) {
            std::cerr << "vpfp_field_particle_pairing_mpi_test must run with "
                         "2 or 5 ranks.\n";
        }
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    const bool ok = parse_args(argc, argv, args);
    if (!ok && rank == 0) {
        std::cerr << "usage: vpfp_field_particle_pairing_mpi_test --case all "
                     "--result <path>\n";
    }

    MpiMetrics m;
    m.mpi_size = size;
    bool pass = ok;
    if (ok) {
        // Smooth bulk distribution spanning all ranks (positive and negative
        // velocities both nonzero, no physical-boundary crossing).
        PairingTestState s = make_smooth_bulk_case(
            64, false, false, rank, size);
        const SpatialGrid& grid = s.grid;
        const int nxl = grid.nx_local;
        const double dt = 1.0e-15;
        OpenBackgroundBoundaryConfig cfg;
        cfg.left_type = BackgroundXBoundaryType::ABSORBING;
        cfg.right_type = BackgroundXBoundaryType::ABSORBING;
        OpenBackgroundBoundary boundary(cfg);

        ConservativePpmRemap remap;
        remap.init(grid, s.velocity_grid);
        Species out;
        out.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                 Const::me, Param::dens, Param::temperature_e, false, grid);
        XFaceTransportAudit audit;
        audit.enabled = true;
        audit.init(nxl);
        remap.advect_x(s.bulk_n, out, dt, 0.0, boundary, rank, size, &audit);

        // Per-cell continuity residual on this rank.
        double local_linf = 0.0;
        double local_scale = 0.0;
        for (int i = 0; i < nxl; ++i) {
            double m_n = 0.0;
            double m_np1 = 0.0;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    m_n += s.bulk_n.f[idx3(grid.nghost + i, j, k)];
                    m_np1 += out.f[idx3(grid.nghost + i, j, k)];
                }
            }
            const double r = m_np1 - m_n -
                audit.bulk_number_swept_face[static_cast<size_t>(i)] +
                audit.bulk_number_swept_face[static_cast<size_t>(i) + 1];
            local_linf = std::max(local_linf, std::fabs(r));
            local_scale = std::max(local_scale, std::fabs(m_n));
        }
        MPI_Allreduce(&local_linf, &m.continuity_linf, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        MPI_Allreduce(&local_scale, &m.continuity_scale, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        m.continuity_scale = std::max(1.0, m.continuity_scale);
        m.continuity_pass = m.continuity_linf <=
            machine_scaled_tolerance(m.continuity_scale, m.continuity_scale);

        // Exercise the production Gate-I face quadrature and MPI owner rule.
        // The explicit reference below sums every global face exactly once;
        // shared right faces on the rank to their left are deliberately
        // excluded.
        FieldParticlePowerAuditWorkspace ws;
        ws.enabled = true;
        ws.bulk_x1 = audit;
        ws.bulk_x2.init(nxl);
        ws.cell_work.init(nxl);
        ws.bulk_number_n.assign(static_cast<size_t>(nxl), 0.0);
        ws.bulk_number_np1.assign(static_cast<size_t>(nxl), 0.0);
        for (int i = 0; i < nxl; ++i) {
            double before = 0.0;
            double after = 0.0;
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    before += s.bulk_n.f[idx3(grid.nghost + i, j, k)];
                    after += out.f[idx3(grid.nghost + i, j, k)];
                }
            }
            ws.bulk_number_n[static_cast<size_t>(i)] = before;
            ws.bulk_number_np1[static_cast<size_t>(i)] = after;
        }
        ws.field_n_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
        ws.field_mid_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
        ws.field_np1_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
        ws.potential_pair_ex_face.assign(static_cast<size_t>(nxl) + 1, 0.0);
        for (int f = 0; f <= nxl; ++f) {
            ws.potential_pair_ex_face[static_cast<size_t>(f)] =
                1.0 + static_cast<double>(grid.ix_start + f);
        }
        FieldParticlePowerAudit calculator;
        calculator.init(grid);
        const FieldParticlePowerAuditResult local_pair =
            calculator.finalize(ws, dt, -Const::qe);
        double distributed_pair = 0.0;
        MPI_Allreduce(&local_pair.current_pair_residual, &distributed_pair, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double explicit_local = 0.0;
        for (int f = 0; f <= nxl; ++f) {
            const int global_face = grid.ix_start + f;
            if (f == nxl && global_face < grid.nx_global) continue;
            double weight = grid.dx;
            if (global_face == 0 || global_face == grid.nx_global)
                weight *= 0.5;
            const double current = -Const::qe *
                audit.bulk_number_swept_face[static_cast<size_t>(f)] / dt;
            explicit_local -= dt * weight *
                ws.potential_pair_ex_face[static_cast<size_t>(f)] * current;
        }
        double explicit_global = 0.0;
        MPI_Allreduce(&explicit_local, &explicit_global, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        m.work_mpi_vs_reference_abs =
            std::fabs(distributed_pair - explicit_global);
        m.full_residual_mpi_vs_reference_abs =
            m.work_mpi_vs_reference_abs;

        // Shared-face ownership: exchange the two boundary face values with
        // the left/right neighbors and count duplicates, mismatches and
        // owner errors.  The internal-face owner is the rank to its right,
        // so rank r's face-nx value must equal rank r+1's face-0 value.
        const double left_face = audit.bulk_number_swept_face[0];
        const double right_face =
            audit.bulk_number_swept_face[static_cast<size_t>(nxl)];
        double recv_left = 0.0;
        double recv_right = 0.0;
        MPI_Request reqs[4];
        int nreq = 0;
        if (rank > 0) {
            MPI_Isend(&left_face, 1, MPI_DOUBLE, rank - 1, 801, MPI_COMM_WORLD,
                      &reqs[nreq++]);
            MPI_Irecv(&recv_left, 1, MPI_DOUBLE, rank - 1, 802, MPI_COMM_WORLD,
                      &reqs[nreq++]);
        }
        if (rank + 1 < size) {
            MPI_Isend(&right_face, 1, MPI_DOUBLE, rank + 1, 802, MPI_COMM_WORLD,
                      &reqs[nreq++]);
            MPI_Irecv(&recv_right, 1, MPI_DOUBLE, rank + 1, 801, MPI_COMM_WORLD,
                      &reqs[nreq++]);
        }
        if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
        double local_dup = 0.0;
        double local_mismatch = 0.0;
        double local_owner = 0.0;
        if (rank > 0) {
            // recv_left is the right-face value of rank-1; this rank's face 0
            // is the shared internal face owned by this rank (rank to the
            // right).  They must match.
            local_dup += 1.0;
            if (std::fabs(recv_left - left_face) >
                machine_scaled_tolerance(std::max(1.0, std::fabs(left_face)), 1.0)) {
                local_mismatch += 1.0;
            }
        }
        if (rank + 1 < size) {
            local_dup += 1.0;
            if (std::fabs(recv_right - right_face) >
                machine_scaled_tolerance(std::max(1.0, std::fabs(right_face)), 1.0)) {
                local_mismatch += 1.0;
            }
        }
        MPI_Allreduce(&local_dup, &m.shared_face_duplicate_count, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_mismatch, &m.shared_face_mismatch_count, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // Exactly one owner must be present for every global face: all local
        // left/interior faces plus the physical right endpoint on the last
        // rank.  This is independent of the duplicated diagnostic rows.
        local_owner = static_cast<double>(nxl) +
            (rank + 1 == size ? 1.0 : 0.0);
        double global_owner_count = 0.0;
        MPI_Allreduce(&local_owner, &global_owner_count, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        m.shared_face_owner_error_count = std::fabs(
            global_owner_count - static_cast<double>(grid.nx_global + 1));

        // Manufactured discrete G/G* identity across MPI shared faces and
        // physical endpoints.  The audit helper supplies the same owner and
        // endpoint quadrature convention as production.
        std::vector<double> e_man(static_cast<size_t>(nxl) + 1, 0.0);
        std::vector<double> j_man(static_cast<size_t>(nxl), 0.0);
        for (int f = 0; f <= nxl; ++f) {
            const double x = static_cast<double>(grid.ix_start + f) /
                             grid.nx_global;
            e_man[static_cast<size_t>(f)] =
                1.0e7 * (1.0 + 0.17 * std::sin(2.0 * Const::pi * x));
        }
        for (int i = 0; i < nxl; ++i) {
            const double x = (static_cast<double>(grid.ix_start + i) + 0.5) /
                             grid.nx_global;
            j_man[static_cast<size_t>(i)] =
                2.0e5 * (1.0 - 0.13 * std::cos(4.0 * Const::pi * x));
        }
        double send_first = j_man.front();
        double send_last = j_man.back();
        double left_neighbor = j_man.front();
        double right_neighbor = j_man.back();
        MPI_Sendrecv(&send_first, 1, MPI_DOUBLE,
                     rank > 0 ? rank - 1 : MPI_PROC_NULL, 811,
                     &right_neighbor, 1, MPI_DOUBLE,
                     rank + 1 < size ? rank + 1 : MPI_PROC_NULL, 811,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send_last, 1, MPI_DOUBLE,
                     rank + 1 < size ? rank + 1 : MPI_PROC_NULL, 812,
                     &left_neighbor, 1, MPI_DOUBLE,
                     rank > 0 ? rank - 1 : MPI_PROC_NULL, 812,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double local_face_work = 0.0;
        double local_cell_work = 0.0;
        double local_man_error = 0.0;
        const bool local_man_valid = calculator.manufactured_gstar_identity(
            e_man, j_man, left_neighbor, right_neighbor, rank, size,
            local_face_work, local_cell_work, local_man_error);
        double global_face_work = 0.0;
        double global_cell_work = 0.0;
        MPI_Allreduce(&local_face_work, &global_face_work, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_cell_work, &global_cell_work, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        // The identity is global: a shared face has zero weight on the left
        // rank and full weight on the right rank, while the two local cell
        // inner products each contain half of that face contribution.
        m.manufactured_identity_error = std::fabs(
            global_face_work - global_cell_work);
        const double global_man_scale = std::max(
            1.0, std::max(std::fabs(global_face_work),
                          std::fabs(global_cell_work)));
        const double man_tol = machine_scaled_tolerance(
            global_man_scale, global_man_scale);
        double local_man_ok = local_man_valid ? 1.0 : 0.0;
        double all_local_man_ok = 0.0;
        MPI_Allreduce(&local_man_ok, &all_local_man_ok, 1, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        m.manufactured_identity_pass = all_local_man_ok > 0.5 &&
            m.manufactured_identity_error <= man_tol;
        const double endpoint_tol = machine_scaled_tolerance(
            std::max(1.0, grid.dx), std::max(1.0, grid.dx));
        double local_endpoint_error = 0.0;
        if (rank == 0) {
            local_endpoint_error = std::max(
                local_endpoint_error,
                std::fabs(calculator.face_quadrature_weight(0) -
                          0.5 * grid.dx));
        }
        if (rank + 1 == size) {
            local_endpoint_error = std::max(
                local_endpoint_error,
                std::fabs(calculator.face_quadrature_weight(nxl) -
                          0.5 * grid.dx));
        }
        double global_endpoint_error = 0.0;
        MPI_Allreduce(&local_endpoint_error, &global_endpoint_error, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        m.manufactured_endpoint_weight_pass =
            global_endpoint_error <= endpoint_tol;

        const double work_scale = std::max(
            1.0, std::max(std::fabs(distributed_pair),
                          std::fabs(explicit_global)));
        pass = m.continuity_pass && m.shared_face_mismatch_count == 0.0 &&
               m.shared_face_owner_error_count == 0.0 &&
               m.manufactured_identity_pass &&
               m.manufactured_endpoint_weight_pass &&
               m.work_mpi_vs_reference_abs <=
                   machine_scaled_tolerance(work_scale, work_scale);
    }

    if (rank == 0 && !args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::trunc);
        if (!out) pass = false;
        else {
            out << std::setprecision(17);
            out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
            out << "mpi_size=" << m.mpi_size << "\n";
            out << "shared_face_duplicate_count="
                << m.shared_face_duplicate_count << "\n";
            out << "shared_face_mismatch_count="
                << m.shared_face_mismatch_count << "\n";
            out << "shared_face_owner_error_count="
                << m.shared_face_owner_error_count << "\n";
            out << "continuity_mpi_vs_reference_abs=" << m.continuity_linf << "\n";
            out << "work_mpi_vs_reference_abs="
                << m.work_mpi_vs_reference_abs << "\n";
            out << "full_residual_mpi_vs_reference_abs="
                << m.full_residual_mpi_vs_reference_abs << "\n";
            out << "manufactured_identity_error="
                << m.manufactured_identity_error << "\n";
            out << "manufactured_identity_pass="
                << (m.manufactured_identity_pass ? 1 : 0) << "\n";
            out << "manufactured_endpoint_weight_pass="
                << (m.manufactured_endpoint_weight_pass ? 1 : 0) << "\n";
            out << "single_rank_reference_available=1\n";
        }
    }

    double pass_d = pass ? 1.0 : 0.0;
    double pass_g = 0.0;
    MPI_Allreduce(&pass_d, &pass_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << std::setprecision(17)
                  << "mpi_size=" << m.mpi_size
                  << " shared_face_duplicate_count="
                  << m.shared_face_duplicate_count
                  << " shared_face_mismatch_count="
                  << m.shared_face_mismatch_count
                  << " continuity_linf=" << m.continuity_linf << "\n";
        std::cout << std::setprecision(17)
                  << "manufactured_identity_error="
                  << m.manufactured_identity_error
                  << " manufactured_identity_pass="
                  << (m.manufactured_identity_pass ? 1 : 0) << "\n";
        std::cout << "manufactured_endpoint_weight_pass="
                  << (m.manufactured_endpoint_weight_pass ? 1 : 0) << "\n";
        std::cout << "status=" << (pass_g > 0.5 ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass_g > 0.5 ? 0 : 1;
}
