#include "background_coupling_test_support.h"
#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <vector>

namespace {

void log_phase(int rank, const char* phase, double elapsed)
{
    if (rank == 0) {
        std::cout << "dual_u_center_phase=" << phase
                  << " elapsed_s=" << std::scientific
                  << std::setprecision(6) << elapsed << std::endl;
    }
}

size_t uface_index(int ix, int jf, int k)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1) + jf) * Param::Nmu + k;
}

double global_max_difference(const std::vector<double>& a,
                             const std::vector<double>& b)
{
    double local = a.size() == b.size() ? 0.0 :
        std::numeric_limits<double>::infinity();
    if (a.size() == b.size())
        for (size_t i = 0; i < a.size(); ++i)
            local = std::max(local, std::fabs(a[i] - b[i]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double global_max_abs(const std::vector<double>& values)
{
    double local = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
        local = std::max(local, std::fabs(values[i]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const double test_start = MPI_Wtime();

    SpatialGrid sg;
    sg.init(rank, size);
    Species background;
    EMFields fields;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.04, 0.06, 1.0e3, 3, 0.0);
    const double dt = BackgroundCouplingTest::stable_dt(sg);
    log_phase(rank, "initialization_done", MPI_Wtime() - test_start);

    BackgroundCouplingTest::BundleOptions legacy_options;
    legacy_options.fct_enabled = false;
    legacy_options.allow_finite_negative_debt = true;
    legacy_options.coupling_mode =
        VlasovAmpereMidpointSolver::LEGACY_COUPLING;
    const double legacy_start = MPI_Wtime();
    log_phase(rank, "legacy_begin", legacy_start - test_start);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle legacy =
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                  size, dt, legacy_options);
    log_phase(rank, "legacy_end", MPI_Wtime() - legacy_start);

    BackgroundCouplingTest::BundleOptions dual_options = legacy_options;
    dual_options.coupling_mode = VlasovAmpereMidpointSolver::DUAL_U_COUPLING;
    const double dual_start = MPI_Wtime();
    log_phase(rank, "dual_begin", dual_start - test_start);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle dual =
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                  size, dt, dual_options);
    log_phase(rank, "dual_end", MPI_Wtime() - dual_start);

    const size_t expected_x_flux = static_cast<size_t>(sg.nx_local + 1) *
                                   Param::Nvmu;
    const size_t expected_u_flux = static_cast<size_t>(sg.nx_local) *
                                   (Param::Nv + 1) * Param::Nmu;
    const size_t expected_face = static_cast<size_t>(sg.nx_local + 1);
    const size_t expected_cell = static_cast<size_t>(sg.nx_local);
    const bool legacy_contract = legacy.state_advanced &&
        !legacy.operator_failed && legacy.fx_high.size() == expected_x_flux &&
        legacy.fu_high.size() == expected_u_flux &&
        legacy.jn_high.size() == expected_face &&
        legacy.gstar_je_high.size() == expected_face;
    const bool dual_contract = dual.state_advanced && !dual.operator_failed &&
        dual.fx_high.size() == expected_x_flux &&
        dual.fu_high.size() == expected_u_flux &&
        dual.cu_high.size() == expected_u_flux &&
        dual.cu_legacy_center.size() == expected_u_flux &&
        dual.jn_high.size() == expected_face &&
        dual.je_high.size() == expected_cell &&
        dual.gstar_je_high.size() == expected_face &&
        dual.dual_target_jn_cell.size() == expected_cell;
    int local_contract_ok = legacy_contract && dual_contract ? 1 : 0;
    int global_contract_ok = 0;
    MPI_Allreduce(&local_contract_ok, &global_contract_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    if (!global_contract_ok) {
        int local_reasons[2] = {legacy.failure_reason, dual.failure_reason};
        int global_reasons[2] = {0, 0};
        MPI_Allreduce(local_reasons, global_reasons, 2, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
        if (rank == 0) {
            std::cerr << "dual_u_center_bundle_contract_failed"
                      << " mpi_size=" << size
                      << " nx_global=" << sg.nx_global
                      << " nx_local_rank0=" << sg.nx_local
                      << " legacy_failure_reason=" << global_reasons[0]
                      << " dual_failure_reason=" << global_reasons[1]
                      << " legacy_state_advanced=" << legacy.state_advanced
                      << " dual_state_advanced=" << dual.state_advanced
                      << " legacy_sizes=(fx=" << legacy.fx_high.size()
                      << ",fu=" << legacy.fu_high.size()
                      << ",jn=" << legacy.jn_high.size() << ")"
                      << " dual_sizes=(fx=" << dual.fx_high.size()
                      << ",fu=" << dual.fu_high.size()
                      << ",cu=" << dual.cu_high.size()
                      << ",cu_legacy=" << dual.cu_legacy_center.size()
                      << ",jn=" << dual.jn_high.size()
                      << ",je=" << dual.je_high.size() << ")"
                      << " expected=(fx=" << expected_x_flux
                      << ",u=" << expected_u_flux
                      << ",face=" << expected_face
                      << ",cell=" << expected_cell << ")" << std::endl;
            std::ofstream out(
                "output/background_coupling_dual_u_center_test.result");
            if (out) {
                out << "test=background_coupling_dual_u_center_test\n"
                    << "bundle_contract_valid=0\n"
                    << "legacy_failure_reason=" << global_reasons[0] << "\n"
                    << "dual_failure_reason=" << global_reasons[1] << "\n"
                    << "passes=0\n";
            }
        }
        MPI_Finalize();
        return 2;
    }
    log_phase(rank, "bundle_contract_validated", MPI_Wtime() - test_start);

    const double fx_difference = global_max_difference(legacy.fx_high,
                                                        dual.fx_high);
    const double jn_difference = global_max_difference(legacy.jn_high,
                                                        dual.jn_high);
    const double cu_correction = global_max_difference(dual.cu_legacy_center,
                                                        dual.cu_high);
    const double target_scale = global_max_abs(dual.dual_target_jn_cell);
    const double target_replay_relative = dual.dual_u_target_replay_linf /
        std::max(1.0, dual.dual_u_target_replay_scale);
    const double u_replay_relative = dual.dual_u_legacy_operator_replay_linf /
        std::max(1.0, dual.dual_u_legacy_operator_replay_scale);
    const double dual_current_relative = dual.dual_u_current_linf /
        std::max(1.0, target_scale);

    double local_fu_contract = 0.0;
    double local_je_contract = 0.0;
    double local_boundary_coefficient_difference = 0.0;
    double local_correction_divergence = 0.0;
    double local_correction_flux_scale = 0.0;
    const double current_factor = background.charge /
        (background.mass * Const::c * sg.dx);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double acceleration = background.charge *
            fields.Ex[sg.nghost + ix] / (background.mass * Const::c);
        double direct_je = 0.0;
        for (int jf = 0; jf <= Param::Nv; ++jf) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = uface_index(ix, jf, k);
                local_fu_contract = std::max(local_fu_contract,
                    std::fabs(dual.fu_high[id] -
                              acceleration * dual.cu_high[id]));
                if (jf > 0 && jf < Param::Nv)
                    direct_je += current_factor *
                        Stage5::delta_energy(background.cgrid, jf, k) *
                        dual.cu_high[id];
            }
        }
        local_je_contract = std::max(local_je_contract,
            std::fabs(direct_je - dual.je_high[static_cast<size_t>(ix)]));
        for (int k = 0; k < Param::Nmu; ++k) {
            local_boundary_coefficient_difference = std::max(
                local_boundary_coefficient_difference,
                std::fabs(dual.cu_high[uface_index(ix, 0, k)] -
                          dual.cu_legacy_center[uface_index(ix, 0, k)]));
            local_boundary_coefficient_difference = std::max(
                local_boundary_coefficient_difference,
                std::fabs(dual.cu_high[uface_index(ix, Param::Nv, k)] -
                          dual.cu_legacy_center[
                              uface_index(ix, Param::Nv, k)]));
        }
        double telescoping = 0.0;
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const double lower = dual.fu_high[uface_index(ix, j, k)] -
                    legacy.fu_high[uface_index(ix, j, k)];
                const double upper = dual.fu_high[uface_index(ix, j + 1, k)] -
                    legacy.fu_high[uface_index(ix, j + 1, k)];
                telescoping += upper - lower;
                local_correction_flux_scale += std::fabs(lower) +
                                               std::fabs(upper);
            }
        local_correction_divergence = std::max(
            local_correction_divergence, std::fabs(telescoping));
    }
    double contracts[5] = {local_fu_contract, local_je_contract,
        local_boundary_coefficient_difference, local_correction_divergence,
        local_correction_flux_scale};
    MPI_Allreduce(MPI_IN_PLACE, contracts, 5, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    log_phase(rank, "contracts_done", MPI_Wtime() - test_start);

    const BackgroundCouplingTest::Norms legacy_pair =
        BackgroundCouplingTest::face_difference_norms(
            legacy.jn_high, legacy.gstar_je_high, sg);
    const BackgroundCouplingTest::Norms dual_pair =
        BackgroundCouplingTest::face_difference_norms(
            dual.jn_high, dual.gstar_je_high, sg);
    const double flux_scale = std::max(1.0, global_max_abs(legacy.fx_high));
    const double current_scale = std::max(1.0, global_max_abs(legacy.jn_high));
    const double coefficient_scale = std::max(1.0,
        global_max_abs(dual.cu_high));
    const double flux_tolerance = 256.0 *
        std::numeric_limits<double>::epsilon() * flux_scale;
    const double current_tolerance = 256.0 *
        std::numeric_limits<double>::epsilon() * current_scale;
    const double coefficient_tolerance = 2048.0 *
        std::numeric_limits<double>::epsilon() * coefficient_scale;
    const bool finite = legacy.outputs_finite && dual.outputs_finite;
    const bool passes = finite && legacy.state_advanced && dual.state_advanced &&
        !legacy.operator_failed && !dual.operator_failed &&
        dual.background_coupling_mode ==
            VlasovAmpereMidpointSolver::DUAL_U_COUPLING &&
        dual.dual_u_operator_valid && fx_difference <= flux_tolerance &&
        jn_difference <= current_tolerance && cu_correction > 0.0 &&
        target_replay_relative <= 1.0e-11 && u_replay_relative <= 1.0e-11 &&
        dual_current_relative <= 1.0e-11 &&
        contracts[0] <= coefficient_tolerance &&
        contracts[1] <= current_tolerance &&
        contracts[2] <= coefficient_tolerance &&
        contracts[3] <= 4096.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, contracts[4]) &&
        dual_pair.l2 < legacy_pair.l2;

    if (rank == 0) {
        std::ofstream out("output/background_coupling_dual_u_center_test.result");
        std::ostream& log = out ? out : std::cout;
        log << std::scientific << std::setprecision(17)
            << "test=background_coupling_dual_u_center_test\n"
            << "beam_enabled=0\nfct_enabled=0\nfixed_candidate=1\n"
            << "legacy_state_advanced=" << legacy.state_advanced << "\n"
            << "dual_state_advanced=" << dual.state_advanced << "\n"
            << "dual_u_operator_valid=" << dual.dual_u_operator_valid << "\n"
            << "fx_high_linf_difference=" << fx_difference << "\n"
            << "jn_high_linf_difference=" << jn_difference << "\n"
            << "cu_dual_correction_linf=" << cu_correction << "\n"
            << "wN_target_replay_relative=" << target_replay_relative << "\n"
            << "wE_legacy_replay_relative=" << u_replay_relative << "\n"
            << "dual_current_target_relative=" << dual_current_relative << "\n"
            << "cu_to_fu_contract_linf=" << contracts[0] << "\n"
            << "cu_to_je_contract_linf=" << contracts[1] << "\n"
            << "u_boundary_coefficient_difference=" << contracts[2] << "\n"
            << "u_correction_mass_divergence_linf=" << contracts[3] << "\n"
            << "u_correction_mass_divergence_relative=" << contracts[3] /
                std::max(1.0, contracts[4]) << "\n"
            << "legacy_pair_L1=" << legacy_pair.l1 << "\n"
            << "legacy_pair_L2=" << legacy_pair.l2 << "\n"
            << "legacy_pair_Linf=" << legacy_pair.linf << "\n"
            << "dual_pair_L1=" << dual_pair.l1 << "\n"
            << "dual_pair_L2=" << dual_pair.l2 << "\n"
            << "dual_pair_Linf=" << dual_pair.linf << "\n"
            << "pair_L2_reduction=" << legacy_pair.l2 /
                std::max(dual_pair.l2, 1.0e-300) << "\n"
            << "dual_u_corrected_cell_count="
            << dual.dual_u_corrected_cell_count << "\n"
            << "passes=" << passes << "\n";
        std::cout << "background_coupling_dual_u_center_test result="
                  << "output/background_coupling_dual_u_center_test.result"
                  << " passes=" << passes << "\n";
    }
    MPI_Finalize();
    return passes ? 0 : 1;
}
