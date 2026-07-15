#include "beam_pic.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <vector>

#ifndef FP_TEST_RUN_ID
#define FP_TEST_RUN_ID "unidentified"
#endif

namespace {

struct Record {
    const char* u_high_candidate_mode;
    double field;
    double dt;
    double jn_mean;
    double je_mean;
    double jn_minus_gstar_je_mean;
    double jn_linf;
    double je_linf;
    double jn_minus_gstar_je_linf;
    double r_couple;
    double r_couple_centered;
    double r_couple_upwind_stabilization;
    double r_couple_fct_stabilization;
    double limiter_active_fraction;
    double limiter_min_alpha;
    int fct_mixed_with_low_order;
    int force_substeps;
    int state_advanced;
    int converged;
    int soft_unconverged;
    int failed;
    int finite;
};

double global_face_mean(const std::vector<double>& values,
                        const SpatialGrid& sg)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface)
        local += values[static_cast<size_t>(iface)];
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global / static_cast<double>(sg.nx_global);
}

double global_face_linf(const std::vector<double>& values,
                        const SpatialGrid& sg)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface)
        local = std::max(local, std::fabs(values[static_cast<size_t>(iface)]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double global_face_difference_linf(const std::vector<double>& lhs,
                                   const std::vector<double>& rhs,
                                   const SpatialGrid& sg)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const size_t i = static_cast<size_t>(iface);
        local = std::max(local, std::fabs(lhs[i] - rhs[i]));
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

bool finite_step(const VlasovAmpereMidpointSolver::Result& step,
                 const SpatialGrid& sg)
{
    int local = 1;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const size_t i = static_cast<size_t>(iface);
        if (!std::isfinite(step.j_bkg_face_mid[i]) ||
            !std::isfinite(step.j_bkg_energy_debug_face[i]) ||
            !std::isfinite(step.fields_np1.Ex_face[i])) {
            local = 0;
            break;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return local != 0;
}

Record run_case(const SpatialGrid& sg, int rank, int size,
                double field_amplitude, double dt,
                bool legacy_boundary_upwind)
{
    Species background;
    background.init("energy_current_scaling_background",
                    SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
                    Param::dens, Param::temperature_e, false, sg);
    background.initialize_maxwellian();

    EMFields fields;
    fields.init(sg);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), field_amplitude);
    fields.sync_cell_ex_from_faces(rank, size);

    BeamPIC beam;
    VlasovAmpereMidpointSolver solver;
    solver.set_step_diagnostics_enabled(false);
    solver.set_beam_enabled(false);
    solver.set_low_order_only(false);
    solver.set_nonuniform_high_order_enabled(true);
    solver.set_fct_enabled(true);
    // This only selects a legacy branch inside the production solver.  The
    // test does not reconstruct, limit, or otherwise duplicate a transport
    // flux; the default branch is the production centered candidate.
    solver.set_legacy_boundary_upwind_high_candidate_for_test(
        legacy_boundary_upwind);
    solver.set_fct_activation_audit_enabled(false);
    const VlasovAmpereMidpointSolver::Result step =
        solver.advance_background_and_fields(background, beam, fields, sg,
                                             dt, 0.0, rank, size);

    Record record = {};
    record.u_high_candidate_mode = legacy_boundary_upwind
        ? "legacy_boundary_upwind_test_control" : "production_centered_all_x";
    record.field = field_amplitude;
    record.dt = dt;
    record.jn_mean = global_face_mean(step.j_bkg_face_mid, sg);
    record.je_mean = global_face_mean(step.j_bkg_energy_debug_face, sg);
    record.jn_minus_gstar_je_mean = record.jn_mean - record.je_mean;
    record.jn_linf = global_face_linf(step.j_bkg_face_mid, sg);
    record.je_linf = global_face_linf(step.j_bkg_energy_debug_face, sg);
    record.jn_minus_gstar_je_linf = global_face_difference_linf(
        step.j_bkg_face_mid, step.j_bkg_energy_debug_face, sg);
    record.r_couple = step.stage5_r_couple;
    record.r_couple_centered = step.stage5_r_couple_centered;
    record.r_couple_upwind_stabilization =
        step.stage5_r_couple_upwind_stabilization;
    record.r_couple_fct_stabilization = step.stage5_r_couple_fct_stabilization;
    record.limiter_active_fraction = step.limiter_active_fraction;
    record.limiter_min_alpha = step.limiter_min_alpha;
    record.fct_mixed_with_low_order =
        step.limiter_active_fraction > 0.0 ? 1 : 0;
    record.force_substeps = step.substeps_used;
    record.state_advanced = step.state_advanced ? 1 : 0;
    record.converged = step.converged ? 1 : 0;
    record.soft_unconverged = step.soft_unconverged ? 1 : 0;
    record.failed = step.failed ? 1 : 0;
    record.finite = finite_step(step, sg) ? 1 : 0;
    return record;
}

void write_record(std::ostream& out, const Record& r)
{
    out << "u_high_candidate_mode=" << r.u_high_candidate_mode
        << " E=" << r.field
        << " dt=" << r.dt
        << " JN_mean=" << r.jn_mean
        << " JE_mean=" << r.je_mean
        << " JN_minus_GstarJE_mean=" << r.jn_minus_gstar_je_mean
        << " JN_Linf=" << r.jn_linf
        << " JE_Linf=" << r.je_linf
        << " JN_minus_GstarJE_Linf=" << r.jn_minus_gstar_je_linf
        << " R_couple=" << r.r_couple
        << " R_couple_centered=" << r.r_couple_centered
        << " R_couple_upwind_stabilization="
        << r.r_couple_upwind_stabilization
        << " R_couple_fct_stabilization=" << r.r_couple_fct_stabilization
        << " limiter_active_fraction=" << r.limiter_active_fraction
        << " limiter_min_alpha=" << r.limiter_min_alpha
        << " fct_mixed_with_low_order=" << r.fct_mixed_with_low_order
        << " force_substeps=" << r.force_substeps
        << " state_advanced=" << r.state_advanced
        << " converged=" << r.converged
        << " soft_unconverged=" << r.soft_unconverged
        << " failed=" << r.failed
        << " finite=" << r.finite << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SpatialGrid sg;
    sg.init(rank, size);
    Species reference;
    reference.init("energy_current_scaling_reference",
                   SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
                   Param::dens, Param::temperature_e, false, sg);
    const double min_du = *std::min_element(reference.cgrid.upar_widths.begin(),
                                            reference.cgrid.upar_widths.end());
    const double e0 = 1.0e3;
    const double acceleration = Const::qe * e0 / (Const::me * Const::c);
    // One common step, chosen from the largest amplitude, makes the signed
    // amplitude comparison a true linear-response experiment.
    const double dt = 0.05 * std::min(0.8 * sg.dx / Const::c,
                                      0.8 * min_du / acceleration);
    const double amplitudes[] = {e0, -e0, 0.5 * e0, -0.5 * e0,
                                 0.25 * e0, -0.25 * e0};
    std::vector<Record> records;
    // A/B control lives only in the production solver option: A restores the
    // legacy boundary-upwind branch, while B is the production centered path.
    for (int candidate = 0; candidate < 2; ++candidate) {
        const bool legacy_boundary_upwind = candidate == 0;
        for (size_t i = 0; i < sizeof(amplitudes) / sizeof(amplitudes[0]); ++i)
            records.push_back(run_case(sg, rank, size, amplitudes[i], dt,
                                       legacy_boundary_upwind));
    }

    int local_ok = 1;
    for (size_t i = 0; i < records.size(); ++i) {
        const Record& r = records[i];
        if (!r.finite || !r.state_advanced || !r.converged ||
            r.soft_unconverged || r.failed || r.force_substeps != 1)
            local_ok = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &local_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream file("output/production_energy_current_scaling.result");
        std::cout << std::scientific << std::setprecision(17)
                  << "production_energy_current_scaling_test\n"
                  << "test_run_id=" << FP_TEST_RUN_ID << "\n"
                  << "beam_enabled=0\ntransport=high_fct\n"
                  << "ab_control=A:legacy_boundary_upwind_test_control,"
                  << "B:production_centered_all_x\n";
        if (file) {
            file << std::scientific << std::setprecision(17)
                 << "test=production_energy_current_scaling\n"
                 << "test_run_id=" << FP_TEST_RUN_ID << "\n"
                 << "beam_enabled=0\ntransport=high_fct\n"
                 << "ab_control=A:legacy_boundary_upwind_test_control,"
                 << "B:production_centered_all_x\n";
        }
        for (size_t i = 0; i < records.size(); ++i) {
            write_record(std::cout, records[i]);
            if (file) write_record(file, records[i]);
        }
        std::cout << "strict_solver_run=" << local_ok << "\n";
        if (file) file << "strict_solver_run=" << local_ok << "\n";
    }
    MPI_Finalize();
    return local_ok ? 0 : 1;
}
