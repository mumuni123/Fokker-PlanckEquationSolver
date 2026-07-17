#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <mpi.h>

namespace {
double max_difference(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double value = 0.0;
    for (size_t i = 0; i < a.size(); ++i) value = std::max(value, std::fabs(a[i] - b[i]));
    return value;
}
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
    SpatialGrid grid; grid.init(rank, size);
    Species bkg;
    bkg.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
             Param::dens, Param::temperature_e, false, grid);
    bkg.initialize_maxwellian();
    EMFields fields; fields.init(grid);
    for (int i = 0; i < grid.nx_local; ++i) {
        const double phase = 2.0 * Const::pi * (grid.ix_start + i + 0.5) / Param::nx;
        fields.Ex_face[static_cast<size_t>(i)] = 1.0e6 * std::cos(phase);
    }
    fields.sync_cell_ex_from_faces(rank, size);
    BeamPIC beam; beam.init(grid);
    VlasovAmpereMidpointSolver solver;
    solver.set_beam_enabled(false);
    solver.set_nonuniform_high_order_enabled(true);
    solver.set_fct_enabled(true);
    solver.set_max_midpoint_iterations(40);
    solver.synchronize_background_ghosts(bkg, grid, rank, size);
    const double dt = 0.05 * Param::dx / Const::c;
    const VlasovAmpereMidpointSolver::CouplingRegionLayout layout = {
        -1, 0.8 * Const::micro};
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation reference =
        solver.evaluate_fixed_midpoint_operator(bkg, beam, fields, bkg, fields,
            std::vector<double>(static_cast<size_t>(grid.nx_local + 1), 0.0),
            layout, grid, dt, 0.0, rank, size);
    VlasovAmpereMidpointSolver::MidpointAuditState written;
    written.step = 1; written.time_s = dt; written.dt_s = dt;
    written.substeps_used = reference.substeps_used;
    written.nonlinear_iterations = reference.nonlinear_iterations;
    written.low_order_only = solver.low_order_only();
    written.high_order_enabled = solver.nonuniform_high_order_enabled();
    written.fct_enabled = solver.fct_enabled();
    written.acceptance_type = reference.converged ? "strict" : "soft";
    written.bkg_n = bkg; written.guess_np1 = bkg;
    written.operator_input_guess = bkg;
    written.fields_n = fields; written.fields_end_guess = fields;
    written.coupling_layout = layout;
    written.j_beam_face_mid = reference.j_beam_face_mid;
    written.reference_jn_face = reference.j_bkg_face_mid;
    written.reference_je_cell = reference.j_bkg_energy_cell_mid;
    written.reference_gstar_je_face = reference.j_bkg_energy_debug_face;
    written.periodic_seam_face_audit = reference.periodic_seam_face_audit;
    written.reference_stage5_r_fv = reference.stage5_r_fv;
    written.reference_stage5_r_couple = reference.stage5_r_couple;
    std::string error;
    bool ok = reference.state_advanced && !reference.failed &&
        write_midpoint_audit_state("fixed_midpoint_operator_audit_tmp", written,
                                   grid, rank, size, error);
    VlasovAmpereMidpointSolver::MidpointAuditState read;
    ok = ok && read_midpoint_audit_state("fixed_midpoint_operator_audit_tmp", read,
                                         bkg, fields, grid, rank, size, error);
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation first =
        solver.evaluate_fixed_midpoint_operator(read.bkg_n, beam, read.fields_n,
            read.operator_input_guess, read.fields_end_guess, read.j_beam_face_mid,
            read.coupling_layout, grid, read.dt_s, read.time_s, rank, size);
    const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation second =
        solver.evaluate_fixed_midpoint_operator(read.bkg_n, beam, read.fields_n,
            read.operator_input_guess, read.fields_end_guess, read.j_beam_face_mid,
            read.coupling_layout, grid, read.dt_s, read.time_s, rank, size);
    double local_difference = std::max(
        max_difference(first.j_bkg_face_mid, read.reference_jn_face),
        max_difference(first.j_bkg_energy_cell_mid, read.reference_je_cell));
    local_difference = std::max(local_difference,
        max_difference(first.j_bkg_energy_debug_face,
                       read.reference_gstar_je_face));
    local_difference = std::max(local_difference,
        std::fabs(first.stage5_r_fv - read.reference_stage5_r_fv));
    local_difference = std::max(local_difference,
        std::fabs(first.stage5_r_couple - read.reference_stage5_r_couple));
    local_difference = std::max(local_difference,
        max_difference(first.j_beam_face_mid, read.j_beam_face_mid));
    local_difference = std::max(local_difference,
        max_difference(first.j_bkg_face_mid, second.j_bkg_face_mid));
    local_difference = std::max(local_difference,
        max_difference(first.j_bkg_energy_cell_mid,
                       second.j_bkg_energy_cell_mid));
    for (size_t i = 0; i < first.periodic_seam_face_audit.size(); ++i) {
        local_difference = std::max(local_difference,
            std::fabs(first.periodic_seam_face_audit[i] -
                      read.periodic_seam_face_audit[i]));
        local_difference = std::max(local_difference,
            std::fabs(first.periodic_seam_face_audit[i] -
                      second.periodic_seam_face_audit[i]));
    }
    if (first.coupling_beam_front_ix != read.coupling_layout.beam_front_ix) {
        local_difference = std::numeric_limits<double>::infinity();
    }
    local_difference = std::max(local_difference,
        std::fabs(first.coupling_wave_core_end_m -
                  read.coupling_layout.wave_core_end_m));
    MPI_Allreduce(MPI_IN_PLACE, &local_difference, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double reference_scale = std::max(1.0,
        std::max(max_difference(read.reference_jn_face,
                                std::vector<double>(read.reference_jn_face.size(), 0.0)),
                 max_difference(read.reference_je_cell,
                                std::vector<double>(read.reference_je_cell.size(), 0.0))));
    const double closure_scale = std::max(
        std::fabs(read.reference_stage5_r_fv),
        std::fabs(read.reference_stage5_r_couple));
    double periodic_seam_scale = 0.0;
    for (size_t i = 0; i < read.periodic_seam_face_audit.size(); ++i) {
        periodic_seam_scale = std::max(periodic_seam_scale,
            std::fabs(read.periodic_seam_face_audit[i]));
    }
    const double comparison_scale = std::max(std::max(reference_scale, closure_scale),
                                             periodic_seam_scale);
    const int local_ok = ok && first.state_advanced && second.state_advanced &&
                         !first.failed && !second.failed &&
                         local_difference <= 4096.0 * std::numeric_limits<double>::epsilon() *
                                             std::max(1.0, comparison_scale);
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf(
            "fixed_midpoint_operator_audit deterministic_difference=%.17e "
            "R_FV_reference=%.17e R_FV_difference=%.17e "
            "R_couple_reference=%.17e R_couple_difference=%.17e status=%s\n",
            local_difference, read.reference_stage5_r_fv,
            std::fabs(first.stage5_r_fv - read.reference_stage5_r_fv),
            read.reference_stage5_r_couple,
            std::fabs(first.stage5_r_couple - read.reference_stage5_r_couple),
            global_ok ? "PASS" : "FAIL");
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
