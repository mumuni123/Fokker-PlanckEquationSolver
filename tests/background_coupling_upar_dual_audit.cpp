#include "background_coupling_test_support.h"
#include "discrete_moment_operators.h"
#include "nonuniform_reconstruction.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <climits>
#include <cstdio>

namespace {

size_t xface_index(int iface, int j, int k)
{
    return (static_cast<size_t>(iface) * Param::Nv + j) * Param::Nmu + k;
}

size_t uface_index(int ix, int jf, int k)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1) + jf) * Param::Nmu + k;
}

bool is_low_u(double u) { return std::fabs(u) <= 0.10; }
bool is_thermal_u(double u) { return std::fabs(u) <= 1.00; }

int zone_for_u(double u)
{
    if (is_low_u(u)) return 0;
    if (is_thermal_u(u)) return 1;
    return 2;
}

struct NormAccumulator {
    double signed_sum;
    double abs_sum;
    double square_sum;
    double linf;

    NormAccumulator() : signed_sum(0.0), abs_sum(0.0), square_sum(0.0),
                        linf(0.0) {}
    void add(double value) {
        signed_sum += value;
        abs_sum += std::fabs(value);
        square_sum += value * value;
        linf = std::max(linf, std::fabs(value));
    }
};

enum AuditFailureBit {
    U_REPLAY_FAILED = 0x01,
    X_REPLAY_FAILED = 0x02,
    JACOBIAN_ACTION_FAILED = 0x04,
    ADJOINT_FAILED = 0x08,
    U_COVERAGE_FAILED = 0x10,
    X_COVERAGE_FAILED = 0x20,
    ACTIVE_COVERAGE_FAILED = 0x40,
    NON_FINITE_FAILED = 0x80
};

struct MaxLocation {
    double value;
    int rank;
};

MaxLocation global_max_location(double local_value, int rank)
{
    struct {
        double value;
        int rank;
    } local = {local_value, rank}, global = {0.0, 0};
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                  MPI_COMM_WORLD);
    MaxLocation result = {global.value, global.rank};
    return result;
}

int first_failed_rank(int local_mask, int rank)
{
    const int candidate = local_mask == 0 ? INT_MAX : rank;
    int first = INT_MAX;
    MPI_Allreduce(&candidate, &first, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return first == INT_MAX ? -1 : first;
}

struct CenteredFaceAudit {
    double coefficient;
    double production_coefficient;
    int branch_signature;
    int stencil[4];
    NonuniformMuscl::FrozenCenteredLinearization frozen;
};

struct XFaceAudit {
    double flux;
    double frozen_flux;
    int branch_signature;
    int stencil[4];
    NonuniformMuscl::FrozenCenteredLinearization frozen;
};

XFaceAudit x_flux_from_mass_line(const CylindricalVelocityGrid& grid,
                                 const std::vector<double>& mass,
                                 int iface, int j, int k,
                                 const SpatialGrid& sg)
{
    const int ng = sg.nghost;
    const int il = ng + iface - 1;
    const int ir = ng + iface;
    const double area = grid.uperp_ring_areas[k];
    const double volume = sg.dx * grid.upar_widths[j] * area;
    const auto fbar = [&](int storage_ix) {
        return volume > 0.0 ? mass[static_cast<size_t>(storage_ix)] / volume
                            : 0.0;
    };
    const auto x_center = [&](int storage_ix) {
        return sg.x_min +
            (sg.ix_start + storage_ix - ng + 0.5) * sg.dx;
    };
    const double xface = sg.x_min + (sg.ix_start + iface) * sg.dx;
    const double q_im1 = fbar(il - 1);
    const double q_i = fbar(il);
    const double q_ip1 = fbar(ir);
    const double q_ip2 = fbar(ir + 1);
    const NonuniformMuscl::FaceStates states =
        NonuniformMuscl::reconstruct_face(
            q_im1, q_i, q_ip1, q_ip2,
            x_center(il - 1), x_center(il), x_center(ir),
            x_center(ir + 1), xface);
    const double speed = grid.vx[idx2(j, k)];
    XFaceAudit result = {};
    result.stencil[0] = il - 1;
    result.stencil[1] = il;
    result.stencil[2] = ir;
    result.stencil[3] = ir + 1;
    result.frozen = NonuniformMuscl::frozen_upwind_linearization(
        q_im1, q_i, q_ip1, q_ip2, x_center(il - 1), x_center(il),
        x_center(ir), x_center(ir + 1), xface, speed >= 0.0);
    result.flux = speed * NonuniformMuscl::upwind_state(states, speed) *
        grid.upar_widths[j] * area;
    const double frozen_state = result.frozen.coefficient[0] * q_im1 +
        result.frozen.coefficient[1] * q_i +
        result.frozen.coefficient[2] * q_ip1 +
        result.frozen.coefficient[3] * q_ip2;
    result.frozen_flux = speed * frozen_state * grid.upar_widths[j] * area;
    result.branch_signature = result.frozen.branch_signature;
    return result;
}

CenteredFaceAudit centered_coefficient_from_mass_line(
    const CylindricalVelocityGrid& grid, const std::vector<double>& mass,
    int jf, int k, double dx)
{
    const int jl = jf - 1;
    const int jr = jf;
    const int jll = (jl == 0) ? jl : jl - 1;
    const int jrr = (jr + 1 == Param::Nv) ? jr : jr + 1;
    const double area = grid.uperp_ring_areas[k];
    const auto fbar = [&](int j) {
        const double volume = dx * grid.upar_widths[j] * area;
        return volume > 0.0 ? mass[static_cast<size_t>(j)] / volume : 0.0;
    };
    const double s_jll = (jl == 0)
        ? grid.upar_cells[jl] - grid.upar_widths[jl]
        : grid.upar_cells[jll];
    const double s_jrr = (jr + 1 == Param::Nv)
        ? grid.upar_cells[jr] + grid.upar_widths[jr]
        : grid.upar_cells[jrr];
    const double q_jll = fbar(jll);
    const double q_jl = fbar(jl);
    const double q_jr = fbar(jr);
    const double q_jrr = fbar(jrr);
    const NonuniformMuscl::FaceStates states =
        NonuniformMuscl::reconstruct_face(q_jll, q_jl, q_jr, q_jrr,
                                          s_jll, grid.upar_cells[jl],
                                          grid.upar_cells[jr], s_jrr,
                                          grid.upar_faces[jf]);
    CenteredFaceAudit result = {};
    result.stencil[0] = jll;
    result.stencil[1] = jl;
    result.stencil[2] = jr;
    result.stencil[3] = jrr;
    result.frozen = NonuniformMuscl::frozen_centered_linearization(
        q_jll, q_jl, q_jr, q_jrr, s_jll, grid.upar_cells[jl],
        grid.upar_cells[jr], s_jrr, grid.upar_faces[jf]);
    const double frozen_state = result.frozen.coefficient[0] * q_jll +
        result.frozen.coefficient[1] * q_jl +
        result.frozen.coefficient[2] * q_jr +
        result.frozen.coefficient[3] * q_jrr;
    result.coefficient = NonuniformMuscl::upar_face_coefficient(
        frozen_state, dx, area);
    result.production_coefficient = NonuniformMuscl::upar_face_coefficient(
        NonuniformMuscl::centered_state(states), dx, area);
    result.branch_signature = result.frozen.branch_signature;
    return result;
}

void reduce_sum(std::vector<double>& values)
{
    if (!values.empty())
        MPI_Allreduce(MPI_IN_PLACE, values.data(),
                      static_cast<int>(values.size()), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
}

std::string option_value(int argc, char** argv, const char* name,
                         const std::string& fallback)
{
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == name) return argv[i + 1];
    return fallback;
}

bool has_option(int argc, char** argv, const char* name)
{
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == name) return true;
    return false;
}

void ensure_directory(const std::string& directory)
{
#ifdef _WIN32
    const std::string command = "mkdir \"" + directory + "\" >nul 2>nul";
#else
    const std::string command = "mkdir -p \"" + directory + "\"";
#endif
    (void)std::system(command.c_str());
}

const char* zone_name(int zone)
{
    static const char* const names[] = {
        "u_endpoint", "empty_tail", "high_u_active", "thermal_body", "low_u"
    };
    return names[zone];
}

int audit_zone(const CylindricalVelocityGrid& grid, int j, int k,
               double mass, double maximum_mass)
{
    (void)k;
    if (j == 0 || j + 1 == Param::Nv) return 0;
    if (maximum_mass > 0.0 && mass / maximum_mass <= 1.0e-16) return 1;
    const double u = std::fabs(grid.upar_cells[j]);
    if (u > 1.0) return 2;
    if (u > 0.1) return 3;
    return 4;
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
    Species background;
    EMFields fields;
    const bool controlled_fct = has_option(argc, argv, "--controlled-fct");
    const bool requested_fct = controlled_fct || has_option(argc, argv, "--fct");
    const bool high_only = has_option(argc, argv, "--high-only");
    const bool negative_field = has_option(argc, argv, "--negative-field");
    const bool quarter_phase = has_option(argc, argv, "--quarter-phase");
    const std::string midpoint_directory = option_value(argc, argv,
        "--midpoint-audit-dir", "");
    const double phase = quarter_phase ? 0.5 * Const::pi : 0.0;
    const double field_amplitude = negative_field ? -1.0e3 : 1.0e3;
    std::ostringstream default_directory_stream;
    default_directory_stream << "output/upar_dual_audit/manufactured/"
        << (controlled_fct ? "controlled_fct" :
            (requested_fct ? "final_fct" : "no_fct")) << "/"
        << (negative_field ? "Eminus" : "Eplus") << "_"
        << (quarter_phase ? "phase_quarter" : "phase0") << "/Nu"
        << std::setw(3) << std::setfill('0') << Param::Nv;
    const std::string output_directory = option_value(argc, argv,
        "--output-dir", default_directory_stream.str());
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.04, 0.06, field_amplitude, 3,
        phase);

    BackgroundCouplingTest::BundleOptions options;
    options.fct_enabled = requested_fct;
    options.controlled_fct_injection = controlled_fct;
    options.allow_finite_negative_debt = true;
    VlasovAmpereMidpointSolver::MidpointAuditState midpoint = {};
    bool checkpoint_loaded = false;
    std::string checkpoint_error;
    if (!midpoint_directory.empty()) {
        checkpoint_loaded = read_midpoint_audit_state(midpoint_directory,
            midpoint, background, fields, sg, rank, size, checkpoint_error);
        if (checkpoint_loaded) {
            background = midpoint.bkg_n;
            fields = midpoint.fields_n;
            options.guess_np1 = &midpoint.operator_input_guess;
            options.fields_end_guess = &midpoint.fields_end_guess;
            options.coupling_layout = &midpoint.coupling_layout;
            options.fct_enabled = high_only ? false :
                (requested_fct || midpoint.fct_enabled);
        }
    }
    if (!midpoint_directory.empty() && !checkpoint_loaded) {
        if (rank == 0) {
            std::cerr << "FATAL checkpoint_audit_read_failed: directory="
                      << midpoint_directory << " reason=" << checkpoint_error
                      << " expected_mpi_size=" << size
                      << " expected_grid=(nx=" << Param::nx
                      << ",Nu=" << Param::Nv << ",Nuperp=" << Param::Nmu
                      << ")\n";
        }
        MPI_Finalize();
        return 2;
    }
    const double dt = checkpoint_loaded ? midpoint.dt_s :
        BackgroundCouplingTest::stable_dt(sg);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        checkpoint_loaded || midpoint_directory.empty() ?
        BackgroundCouplingTest::evaluate_bundle(background, fields, sg, rank,
                                                size, dt, options) :
        VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle();

    const size_t nxf = static_cast<size_t>(sg.nx_local + 1) * Param::Nvmu;
    const size_t nuf = static_cast<size_t>(sg.nx_local) *
        (Param::Nv + 1) * Param::Nmu;
    const bool raw_flux_arrays_ok = bundle.fx_low.size() == nxf &&
        bundle.fx_high.size() == nxf && bundle.fx_final.size() == nxf &&
        bundle.fu_low.size() == nuf && bundle.fu_center.size() == nuf &&
        bundle.fu_high.size() == nuf && bundle.fu_final.size() == nuf &&
        bundle.cu_center.size() == nuf &&
        bundle.cu_reconstruction_mass.size() >=
            static_cast<size_t>(sg.nx_total) * Param::Nvmu;
    const size_t current_count = static_cast<size_t>(sg.nx_local);
    const bool current_arrays_ok = bundle.jn_low.size() >= current_count &&
        bundle.jn_high.size() >= current_count &&
        bundle.jn_final.size() >= current_count &&
        bundle.gstar_je_low.size() >= current_count &&
        bundle.gstar_je_center.size() >= current_count &&
        bundle.gstar_je_high.size() >= current_count &&
        bundle.gstar_je_final.size() >= current_count;
    int direct_work_valid = (midpoint_directory.empty() || checkpoint_loaded) &&
        bundle.state_advanced && !bundle.operator_failed &&
        bundle.outputs_finite && current_arrays_ok;
    MPI_Allreduce(MPI_IN_PLACE, &direct_work_valid, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    int single_substep_flux_valid = direct_work_valid && raw_flux_arrays_ok &&
        bundle.substeps_used == 1;
    MPI_Allreduce(MPI_IN_PLACE, &single_substep_flux_valid, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    const int jacobian_audit_applicable = single_substep_flux_valid;

    std::vector<double> jn_signed(Param::Nvmu, 0.0);
    std::vector<double> jn_abs(Param::Nvmu, 0.0);
    std::vector<double> je_signed(static_cast<size_t>(Param::Nv + 1) *
                                  Param::Nmu, 0.0);
    std::vector<double> je_abs(static_cast<size_t>(Param::Nv + 1) *
                               Param::Nmu, 0.0);
    std::vector<double> residual_jk(Param::Nvmu, 0.0);
    std::vector<double> residual_jk_abs(Param::Nvmu, 0.0);
    std::vector<double> mass_jk(Param::Nvmu, 0.0);
    std::vector<double> delta_k_lower(Param::Nvmu, 0.0);
    std::vector<double> delta_k_upper(Param::Nvmu, 0.0);
    double local_boundary_work = 0.0;

    int audit_points[7] = {};
    int audit_point_count = 0;
    if (single_substep_flux_valid) {
        const int ng = sg.nghost;
        const double energy_factor = background.charge /
            (background.mass * Const::c * sg.dx);
        for (int iface = 0; iface < sg.nx_local; ++iface) {
            const double eface = fields.Ex_face[static_cast<size_t>(iface)];
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const double current = background.charge *
                        bundle.fx_high[xface_index(iface, j, k)];
                    const size_t jk = idx2(j, k);
                    jn_signed[jk] += current;
                    jn_abs[jk] += std::fabs(current);
                    const double work = dt * sg.dx * eface * current;
                    residual_jk[jk] += work;
                    residual_jk_abs[jk] += std::fabs(work);
                }
            }
        }
        for (int ix = 0; ix < sg.nx_local; ++ix) {
            const double ecell = fields.Ex[ng + ix];
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k)
                    mass_jk[idx2(j, k)] += bundle.cu_reconstruction_mass[
                        idx3(ng + ix, j, k)];
            }
            for (int jf = 1; jf < Param::Nv; ++jf) {
                const int jl = jf - 1;
                const int jr = jf;
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t id = uface_index(ix, jf, k);
                    const double dke = Stage5::delta_energy(background.cgrid,
                                                            jf, k);
                    const double current = energy_factor * dke *
                        bundle.cu_center[id];
                    const size_t jf_k = static_cast<size_t>(jf) *
                        Param::Nmu + k;
                    je_signed[jf_k] += current;
                    je_abs[jf_k] += std::fabs(current);
                    const double work = dt * sg.dx * ecell * current;
                    const size_t left = idx2(jl, k);
                    const size_t right = idx2(jr, k);
                    residual_jk[left] -= 0.5 * work;
                    residual_jk[right] -= 0.5 * work;
                    residual_jk_abs[left] += 0.5 * std::fabs(work);
                    residual_jk_abs[right] += 0.5 * std::fabs(work);
                    delta_k_upper[left] = dke;
                    delta_k_lower[right] = dke;
                }
            }
            for (int k = 0; k < Param::Nmu; ++k) {
                const double lo = bundle.fu_center[uface_index(ix, 0, k)];
                const double hi = bundle.fu_center[
                    uface_index(ix, Param::Nv, k)];
                local_boundary_work += dt * (
                    background.cgrid.kinetic_energy[idx2(Param::Nv - 1, k)] * hi -
                    background.cgrid.kinetic_energy[idx2(0, k)] * lo);
            }
        }
    }
    reduce_sum(jn_signed);
    reduce_sum(jn_abs);
    reduce_sum(je_signed);
    reduce_sum(je_abs);
    reduce_sum(residual_jk);
    reduce_sum(residual_jk_abs);
    reduce_sum(mass_jk);
    MPI_Allreduce(MPI_IN_PLACE, &local_boundary_work, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    std::vector<double> w_n(Param::Nvmu, 0.0);
    std::vector<double> w_n_analytic(Param::Nvmu, 0.0);
    std::vector<double> w_e(Param::Nvmu, 0.0);
    std::vector<int> frozen_face_count(Param::Nvmu, 0);
    std::vector<int> candidate_adjacent_face_count(Param::Nvmu, 0);
    std::vector<int> x_frozen_face_count(Param::Nvmu, 0);
    std::vector<int> x_candidate_face_count(Param::Nvmu, 0);
    std::vector<int> u_kink_adjacent_face_count(Param::Nvmu, 0);
    std::vector<int> x_kink_adjacent_face_count(Param::Nvmu, 0);
    double nonlinear_response_linf = 0.0;
    double nonlinear_response_scale = 0.0;
    double frozen_mapping_linf = 0.0;
    double x_nonlinear_response_linf = 0.0;
    double x_nonlinear_response_scale = 0.0;
    double jacobian_action_linf = 0.0;
    double jacobian_action_scale = 0.0;
    double adjoint_lhs[8] = {};
    double adjoint_rhs[8] = {};
    double epsilon_scale_error[8][3] = {};
    double epsilon_scale_value[8][3] = {};
    double epsilon_scale_common_spread[8] = {};
    long long epsilon_scale_valid_count[8][3] = {};
    long long epsilon_scale_branch_crossed[3] = {};
    long long u_branch_crossed_face_count = 0;
    long long x_branch_crossed_face_count = 0;
    long long u_kink_one_sided_count = 0;
    long long x_kink_one_sided_count = 0;
    long long kink_one_sided_nonfinite_count = 0;
    double kink_one_sided_relative_error_min =
        std::numeric_limits<double>::infinity();
    double kink_one_sided_relative_error_max = 0.0;
    long long frozen_face_total = 0;
    long long candidate_face_total = 0;
    long long frozen_zero_slope_count = 0;
    long long x_frozen_face_total = 0;
    long long x_candidate_face_total = 0;
    const double derivative_eps = std::sqrt(
        std::numeric_limits<double>::epsilon());
    const double epsilon_factors[3] = {0.25, 1.0, 4.0};
    int audit_ix = 0;
    int audit_global_ix = -1;
    double audit_ecell = 0.0;
    if (single_substep_flux_valid) {
        const int interior_begin = std::min(2, std::max(0, sg.nx_local - 1));
        const int interior_end = std::max(interior_begin,
            sg.nx_local - 3);
        audit_ix = interior_begin;
        for (int candidate = interior_begin; candidate <= interior_end;
             ++candidate) {
            if (std::fabs(fields.Ex[sg.nghost + candidate]) >
                std::fabs(fields.Ex[sg.nghost + audit_ix]))
                audit_ix = candidate;
        }
        const int ix = audit_ix;
        const int storage_ix = sg.nghost + ix;
        audit_global_ix = sg.ix_start + ix;
        audit_ecell = fields.Ex[storage_ix];
        for (int k = 0; k < Param::Nmu; ++k) {
            std::vector<double> mass(static_cast<size_t>(Param::Nv), 0.0);
            double line_scale = 0.0;
            for (int j = 0; j < Param::Nv; ++j) {
                mass[static_cast<size_t>(j)] = bundle.cu_reconstruction_mass[
                    idx3(storage_ix, j, k)];
                line_scale = std::max(line_scale,
                    std::fabs(mass[static_cast<size_t>(j)]));
            }
            line_scale = std::max(line_scale, 1.0);
            std::vector<std::vector<double> > directions(8,
                std::vector<double>(static_cast<size_t>(Param::Nv), 0.0));
            for (int direction_id = 0; direction_id < 8; ++direction_id)
                for (int j = 0; j < Param::Nv; ++j)
                    directions[static_cast<size_t>(direction_id)]
                        [static_cast<size_t>(j)] = line_scale * std::sin(
                            (0.371 + 0.031 * direction_id) * (j + 1) +
                            (0.193 + 0.017 * direction_id) * (k + 1));
            for (int jf = 1; jf < Param::Nv; ++jf) {
                const size_t face = uface_index(ix, jf, k);
                const CenteredFaceAudit base =
                    centered_coefficient_from_mass_line(background.cgrid,
                                                        mass, jf, k, sg.dx);
                const double dke = Stage5::delta_energy(background.cgrid,
                                                        jf, k);
                // g_E is the gradient of the actual cell-field work
                // functional W_E, not the older velocity-only proxy.
                const double energy_factor = dt * audit_ecell *
                    background.charge * dke /
                    (background.mass * Const::c);
                nonlinear_response_linf = std::max(nonlinear_response_linf,
                    std::fabs(base.production_coefficient -
                              bundle.cu_center[face]));
                nonlinear_response_scale = std::max(nonlinear_response_scale,
                    std::fabs(bundle.cu_center[face]));
                frozen_mapping_linf = std::max(frozen_mapping_linf,
                    std::fabs(base.coefficient - base.production_coefficient));
                std::vector<double> derivative(static_cast<size_t>(Param::Nv),
                                               0.0);
                bool seen[Param::Nv] = {};
                for (int p = 0; p < 4; ++p) {
                    const int j = base.stencil[p];
                    derivative[static_cast<size_t>(j)] +=
                        base.frozen.coefficient[p] /
                        background.cgrid.upar_widths[j];
                    if (!seen[j]) {
                        seen[j] = true;
                        ++candidate_adjacent_face_count[idx2(j, k)];
                        ++frozen_face_count[idx2(j, k)];
                    }
                }
                ++candidate_face_total;
                ++frozen_face_total;
                if ((base.branch_signature & 0x3) == 0)
                    ++frozen_zero_slope_count;
                if (((base.branch_signature >> 2) & 0x3) == 0)
                    ++frozen_zero_slope_count;
                for (int j = 0; j < Param::Nv; ++j) {
                    w_e[idx2(j, k)] += energy_factor *
                        derivative[static_cast<size_t>(j)];
                }
                for (int direction_id = 0; direction_id < 8; ++direction_id) {
                    double matrix_derivative = 0.0;
                    double frozen_directional = 0.0;
                    for (int j = 0; j < Param::Nv; ++j)
                        matrix_derivative += derivative[static_cast<size_t>(j)] *
                            directions[static_cast<size_t>(direction_id)]
                                [static_cast<size_t>(j)];
                    for (int p = 0; p < 4; ++p) {
                        const int j = base.stencil[p];
                        const double volume = sg.dx *
                            background.cgrid.upar_widths[j] *
                            background.cgrid.uperp_ring_areas[k];
                        frozen_directional += base.frozen.coefficient[p] *
                            directions[static_cast<size_t>(direction_id)]
                                [static_cast<size_t>(j)] / volume;
                    }
                    frozen_directional *= sg.dx *
                        background.cgrid.uperp_ring_areas[k];
                    jacobian_action_linf = std::max(jacobian_action_linf,
                        std::fabs(frozen_directional - matrix_derivative));
                    jacobian_action_scale = std::max(jacobian_action_scale,
                        std::fabs(frozen_directional));
                    adjoint_lhs[direction_id] += energy_factor * frozen_directional;

                    // A frozen analytic branch is the primary Jacobian.  The
                    // three finite-difference scales audit its numerical
                    // stability and explicitly classify branch crossings.
                    CenteredFaceAudit forward[3];
                    CenteredFaceAudit backward[3];
                    double finite_difference[3] = {};
                    bool common_smooth_branch = true;
                    for (int scale_id = 0; scale_id < 3; ++scale_id) {
                        const double eps = epsilon_factors[scale_id] *
                            derivative_eps * std::max(1.0e-300, line_scale);
                        std::vector<double> plus = mass;
                        std::vector<double> minus = mass;
                        for (int q = 0; q < Param::Nv; ++q) {
                            const double perturb = eps *
                                directions[static_cast<size_t>(direction_id)]
                                    [static_cast<size_t>(q)] / line_scale;
                            plus[static_cast<size_t>(q)] += perturb;
                            minus[static_cast<size_t>(q)] -= perturb;
                        }
                        forward[scale_id] = centered_coefficient_from_mass_line(
                            background.cgrid, plus, jf, k, sg.dx);
                        backward[scale_id] = centered_coefficient_from_mass_line(
                            background.cgrid, minus, jf, k, sg.dx);
                        const bool smooth =
                            forward[scale_id].branch_signature == base.branch_signature &&
                            backward[scale_id].branch_signature == base.branch_signature;
                        if (!smooth) {
                            ++epsilon_scale_branch_crossed[scale_id];
                            common_smooth_branch = false;
                        }
                        finite_difference[scale_id] =
                            (forward[scale_id].production_coefficient -
                             backward[scale_id].production_coefficient) /
                            (2.0 * eps);
                        finite_difference[scale_id] *= line_scale;
                    }
                    if (!common_smooth_branch) {
                        ++u_branch_crossed_face_count;
                        for (int p = 0; p < 4; ++p)
                            ++u_kink_adjacent_face_count[idx2(base.stencil[p], k)];
                        // A limiter kink has no unique centered derivative.  Audit
                        // both one-sided derivatives instead of comparing it with
                        // the frozen smooth-branch Jacobian.
                        const int scale_id = 1;
                        const double eps = epsilon_factors[scale_id] *
                            derivative_eps * std::max(1.0e-300, line_scale);
                        const double forward_derivative =
                            (forward[scale_id].production_coefficient -
                             base.production_coefficient) / eps * line_scale;
                        const double backward_derivative =
                            (base.production_coefficient -
                             backward[scale_id].production_coefficient) / eps *
                            line_scale;
                        const double denominator = std::max(1.0e-300,
                            std::max(std::fabs(matrix_derivative),
                            std::max(std::fabs(forward_derivative),
                                     std::fabs(backward_derivative))));
                        if (std::isfinite(forward_derivative) &&
                            std::isfinite(backward_derivative)) {
                            const double forward_error = std::fabs(
                                forward_derivative - matrix_derivative) / denominator;
                            const double backward_error = std::fabs(
                                backward_derivative - matrix_derivative) / denominator;
                            kink_one_sided_relative_error_min = std::min(
                                kink_one_sided_relative_error_min,
                                std::min(forward_error, backward_error));
                            kink_one_sided_relative_error_max = std::max(
                                kink_one_sided_relative_error_max,
                                std::max(forward_error, backward_error));
                            ++u_kink_one_sided_count;
                        } else {
                            ++kink_one_sided_nonfinite_count;
                        }
                        continue;
                    }
                    double fd_min = finite_difference[0];
                    double fd_max = finite_difference[0];
                    for (int scale_id = 0; scale_id < 3; ++scale_id) {
                        const double scale = std::max(1.0e-300,
                            std::max(std::fabs(finite_difference[scale_id]),
                                     std::fabs(matrix_derivative)));
                        epsilon_scale_error[direction_id][scale_id] = std::max(
                            epsilon_scale_error[direction_id][scale_id],
                            std::fabs(finite_difference[scale_id] -
                                      matrix_derivative) / scale);
                        epsilon_scale_value[direction_id][scale_id] +=
                            finite_difference[scale_id];
                        ++epsilon_scale_valid_count[direction_id][scale_id];
                        fd_min = std::min(fd_min, finite_difference[scale_id]);
                        fd_max = std::max(fd_max, finite_difference[scale_id]);
                    }
                    epsilon_scale_common_spread[direction_id] = std::max(
                        epsilon_scale_common_spread[direction_id],
                        (fd_max - fd_min) / std::max(1.0e-300,
                            std::max(std::fabs(fd_min), std::fabs(fd_max))));
                }
            }
            for (int direction_id = 0; direction_id < 8; ++direction_id)
                for (int j = 0; j < Param::Nv; ++j)
                    adjoint_rhs[direction_id] += w_e[idx2(j, k)] *
                        directions[static_cast<size_t>(direction_id)]
                            [static_cast<size_t>(j)];
        }

        // Linearize the actual production x-MUSCL flux about the same
        // midpoint mass state.  The E_face/E_cell factor puts the face work
        // derivative into the same cell-field inner product as w_E.
        if (std::fabs(audit_ecell) > 64.0 *
            std::numeric_limits<double>::epsilon() * 1.0e3) {
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t jk = idx2(j, k);
                    w_n_analytic[jk] = background.cgrid.vx[jk] / sg.dx;
                    std::vector<double> mass_line(
                        static_cast<size_t>(sg.nx_total), 0.0);
                    double line_scale = 0.0;
                    for (int sx = 0; sx < sg.nx_total; ++sx) {
                        mass_line[static_cast<size_t>(sx)] =
                            bundle.cu_reconstruction_mass[idx3(sx, j, k)];
                        line_scale = std::max(line_scale,
                            std::fabs(mass_line[static_cast<size_t>(sx)]));
                    }
                    line_scale = std::max(line_scale, 1.0e-300);
                    const double epsilon = derivative_eps * std::max(
                        std::fabs(mass_line[static_cast<size_t>(storage_ix)]),
                        line_scale * 1.0e-12);
                    const int first_face = std::max(0, ix - 1);
                    const int last_face = std::min(sg.nx_local - 1, ix + 2);
                    for (int iface = first_face; iface <= last_face; ++iface) {
                        ++x_candidate_face_count[jk];
                        ++x_candidate_face_total;
                        const XFaceAudit base = x_flux_from_mass_line(
                            background.cgrid, mass_line, iface, j, k, sg);
                        const size_t face_id = xface_index(iface, j, k);
                        x_nonlinear_response_linf = std::max(
                            x_nonlinear_response_linf,
                            std::fabs(base.flux - bundle.fx_high[face_id]));
                        x_nonlinear_response_scale = std::max(
                            x_nonlinear_response_scale,
                            std::fabs(bundle.fx_high[face_id]));
                        std::vector<double> plus = mass_line;
                        std::vector<double> minus = mass_line;
                        plus[static_cast<size_t>(storage_ix)] += epsilon;
                        minus[static_cast<size_t>(storage_ix)] -= epsilon;
                        const XFaceAudit forward = x_flux_from_mass_line(
                            background.cgrid, plus, iface, j, k, sg);
                        const XFaceAudit backward = x_flux_from_mass_line(
                            background.cgrid, minus, iface, j, k, sg);
                        if (forward.branch_signature != base.branch_signature ||
                            backward.branch_signature != base.branch_signature) {
                            ++x_branch_crossed_face_count;
                            for (int p = 0; p < 4; ++p)
                                ++x_kink_adjacent_face_count[jk];
                            const double forward_derivative =
                                (forward.flux - base.flux) / epsilon;
                            const double backward_derivative =
                                (base.flux - backward.flux) / epsilon;
                            double frozen_derivative = 0.0;
                            for (int p = 0; p < 4; ++p) {
                                if (base.stencil[p] != storage_ix) continue;
                                frozen_derivative += background.cgrid.vx[jk] *
                                    base.frozen.coefficient[p] / sg.dx;
                            }
                            const double denominator = std::max(1.0e-300,
                                std::max(std::fabs(frozen_derivative),
                                std::max(std::fabs(forward_derivative),
                                         std::fabs(backward_derivative))));
                            if (std::isfinite(forward_derivative) &&
                                std::isfinite(backward_derivative)) {
                                const double forward_error = std::fabs(
                                    forward_derivative - frozen_derivative) /
                                    denominator;
                                const double backward_error = std::fabs(
                                    backward_derivative - frozen_derivative) /
                                    denominator;
                                kink_one_sided_relative_error_min = std::min(
                                    kink_one_sided_relative_error_min,
                                    std::min(forward_error, backward_error));
                                kink_one_sided_relative_error_max = std::max(
                                    kink_one_sided_relative_error_max,
                                    std::max(forward_error, backward_error));
                                ++x_kink_one_sided_count;
                            } else {
                                ++kink_one_sided_nonfinite_count;
                            }
                            continue;
                        }
                        ++x_frozen_face_count[jk];
                        ++x_frozen_face_total;
                        const double derivative = (forward.flux - backward.flux) /
                            (2.0 * epsilon);
                        double frozen_derivative = 0.0;
                        for (int p = 0; p < 4; ++p) {
                            if (base.stencil[p] != storage_ix) continue;
                            frozen_derivative += background.cgrid.vx[jk] *
                                base.frozen.coefficient[p] / sg.dx;
                        }
                        w_n[jk] += dt * sg.dx *
                            fields.Ex_face[static_cast<size_t>(iface)] *
                            background.charge * frozen_derivative;
                        w_n_analytic[jk] += dt * sg.dx *
                            fields.Ex_face[static_cast<size_t>(iface)] *
                            background.charge * derivative;
                    }
                }
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, w_n.data(), static_cast<int>(w_n.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, w_n_analytic.data(),
                  static_cast<int>(w_n_analytic.size()), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, w_e.data(), static_cast<int>(w_e.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, frozen_face_count.data(),
                  static_cast<int>(frozen_face_count.size()), MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, candidate_adjacent_face_count.data(),
                  static_cast<int>(candidate_adjacent_face_count.size()), MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, x_frozen_face_count.data(),
                  static_cast<int>(x_frozen_face_count.size()), MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, x_candidate_face_count.data(),
                  static_cast<int>(x_candidate_face_count.size()), MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, u_kink_adjacent_face_count.data(),
                  static_cast<int>(u_kink_adjacent_face_count.size()), MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, x_kink_adjacent_face_count.data(),
                  static_cast<int>(x_kink_adjacent_face_count.size()), MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);

    NormAccumulator zone[5];
    NormAccumulator all;
    NormAccumulator g_n_all;
    NormAccumulator g_e_all;
    NormAccumulator g_r_all;
    NormAccumulator support[3];
    NormAccumulator symmetry;
    double positive_mass = 0.0;
    double weighted_gr_square = 0.0;
    double weighted_gn_square = 0.0;
    double weighted_ge_square = 0.0;
    double weighted_gr_l1 = 0.0;
    double weighted_gn_l1 = 0.0;
    double weighted_ge_l1 = 0.0;
    double weighted_gr_signed = 0.0;
    double kink_active_mass = 0.0;
    double kink_pair_work = 0.0;
    double maximum_mass = 0.0;
    for (size_t id = 0; id < mass_jk.size(); ++id)
        maximum_mass = std::max(maximum_mass, std::max(0.0, mass_jk[id]));
    double tail_weight_difference = 0.0;
    double tail_weight_total = 0.0;
    double worst_abs = -1.0;
    int worst_j = -1;
    int worst_k = -1;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t id = idx2(j, k);
            const double residual = residual_jk[id];
            const double g_n = w_n[id];
            const double g_e = w_e[id];
            const double g_r = g_n - g_e;
            all.add(residual);
            zone[audit_zone(background.cgrid, j, k, mass_jk[id], maximum_mass)].add(residual);
            g_n_all.add(g_n);
            g_e_all.add(g_e);
            g_r_all.add(g_r);
            const double mass = std::max(0.0, mass_jk[id]);
            positive_mass += mass;
            weighted_gr_square += mass * g_r * g_r;
            weighted_gn_square += mass * g_n * g_n;
            weighted_ge_square += mass * g_e * g_e;
            weighted_gr_l1 += std::fabs(mass * g_r);
            weighted_gn_l1 += std::fabs(mass * g_n);
            weighted_ge_l1 += std::fabs(mass * g_e);
            weighted_gr_signed += mass * g_r;
            if (u_kink_adjacent_face_count[id] > 0 ||
                x_kink_adjacent_face_count[id] > 0) {
                kink_active_mass += mass;
                kink_pair_work += std::fabs(residual);
            }
            static const double support_cutoffs[] = {1.0e-8, 1.0e-12, 1.0e-16};
            for (int s = 0; s < 3; ++s)
                if (maximum_mass > 0.0 && mass / maximum_mass > support_cutoffs[s])
                    support[s].add(g_r);
            if (std::fabs(background.cgrid.upar_cells[j]) >=
                0.8 * Param::momentum_umax &&
                candidate_adjacent_face_count[id] > 0 &&
                frozen_face_count[id] == candidate_adjacent_face_count[id] &&
                x_candidate_face_count[id] > 0 &&
                x_frozen_face_count[id] == x_candidate_face_count[id]) {
                tail_weight_difference += std::fabs(w_n[id] - w_e[id]);
                tail_weight_total += std::fabs(w_n[id]);
            }
            if (std::fabs(residual) > worst_abs) {
                worst_abs = std::fabs(residual);
                worst_j = j;
                worst_k = k;
            }
            const int mirror = Param::Nv - 1 - j;
            symmetry.add(residual - residual_jk[idx2(mirror, k)]);
        }
    }
    const double gr_mass_l2 = positive_mass > 0.0 ?
        std::sqrt(weighted_gr_square / positive_mass) : 0.0;
    const double gn_mass_l2 = positive_mass > 0.0 ?
        std::sqrt(weighted_gn_square / positive_mass) : 0.0;
    const double ge_mass_l2 = positive_mass > 0.0 ?
        std::sqrt(weighted_ge_square / positive_mass) : 0.0;
    const double relative_gradient = gr_mass_l2 /
        std::max(1.0e-300, std::max(gn_mass_l2, ge_mass_l2));
    const double relative_work_gradient = weighted_gr_l1 /
        std::max(1.0e-300, std::max(weighted_gn_l1, weighted_ge_l1));
    double frozen_weight_cells = 0.0;
    double candidate_weight_cells = 0.0;
    double active_fully_covered_cells = 0.0;
    double active_candidate_cells = 0.0;
    double w_l1 = 0.0, w_l2 = 0.0, w_linf = 0.0;
    double analytic_w_l2 = 0.0;
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t id = idx2(j, k);
            ++candidate_weight_cells;
            const bool u_covered = candidate_adjacent_face_count[id] > 0 &&
                frozen_face_count[id] == candidate_adjacent_face_count[id];
            const bool x_covered = x_candidate_face_count[id] > 0 &&
                x_frozen_face_count[id] == x_candidate_face_count[id];
            if (!u_covered || !x_covered) continue;
            const double difference = w_n[id] - w_e[id];
            w_l1 += std::fabs(difference);
            w_l2 += difference * difference;
            w_linf = std::max(w_linf, std::fabs(difference));
            const double analytic_difference = w_n_analytic[id] - w_e[id];
            analytic_w_l2 += analytic_difference * analytic_difference;
            frozen_weight_cells += 1.0;
            if (maximum_mass > 0.0 && mass_jk[id] / maximum_mass > 1.0e-16)
                active_fully_covered_cells += 1.0;
        }
    }
    for (size_t id = 0; id < mass_jk.size(); ++id)
        if (maximum_mass > 0.0 && mass_jk[id] / maximum_mass > 1.0e-16)
            active_candidate_cells += 1.0;
    if (frozen_weight_cells > 0.0) {
        w_l1 /= frozen_weight_cells;
        w_l2 = std::sqrt(w_l2 / frozen_weight_cells);
        analytic_w_l2 = std::sqrt(analytic_w_l2 / frozen_weight_cells);
    } else {
        w_l1 = w_l2 = w_linf = analytic_w_l2 =
            std::numeric_limits<double>::quiet_NaN();
    }

    MPI_Allreduce(MPI_IN_PLACE, adjoint_lhs, 8, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, adjoint_rhs, 8, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    long long count_pack[5] = {frozen_face_total, candidate_face_total,
        x_frozen_face_total, x_candidate_face_total, frozen_zero_slope_count};
    MPI_Allreduce(MPI_IN_PLACE, count_pack, 5, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    frozen_face_total = count_pack[0];
    candidate_face_total = count_pack[1];
    x_frozen_face_total = count_pack[2];
    x_candidate_face_total = count_pack[3];
    frozen_zero_slope_count = count_pack[4];
    MPI_Allreduce(MPI_IN_PLACE, epsilon_scale_branch_crossed, 3,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, epsilon_scale_error, 24, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, epsilon_scale_value, 24, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, epsilon_scale_common_spread, 8, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, epsilon_scale_valid_count, 24,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &u_branch_crossed_face_count, 1,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &x_branch_crossed_face_count, 1,
                  MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    long long kink_counts[3] = {u_kink_one_sided_count,
        x_kink_one_sided_count, kink_one_sided_nonfinite_count};
    MPI_Allreduce(MPI_IN_PLACE, kink_counts, 3, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    u_kink_one_sided_count = kink_counts[0];
    x_kink_one_sided_count = kink_counts[1];
    kink_one_sided_nonfinite_count = kink_counts[2];
    MPI_Allreduce(MPI_IN_PLACE, &kink_one_sided_relative_error_min, 1,
                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &kink_one_sided_relative_error_max, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (u_kink_one_sided_count + x_kink_one_sided_count == 0)
        kink_one_sided_relative_error_min = 0.0;
    const MaxLocation global_u_replay = global_max_location(
        nonlinear_response_linf, rank);
    const MaxLocation global_frozen_mapping = global_max_location(
        frozen_mapping_linf, rank);
    const MaxLocation global_x_replay = global_max_location(
        x_nonlinear_response_linf, rank);
    const MaxLocation global_jacobian = global_max_location(
        jacobian_action_linf, rank);
    const MaxLocation global_u_scale = global_max_location(
        nonlinear_response_scale, rank);
    const MaxLocation global_x_scale = global_max_location(
        x_nonlinear_response_scale, rank);
    const MaxLocation global_jacobian_scale = global_max_location(
        jacobian_action_scale, rank);
    const double nonlinear_response_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, global_u_scale.value);
    const double jacobian_action_tolerance = 8192.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, global_jacobian_scale.value);
    const double x_nonlinear_response_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, global_x_scale.value);
    double adjoint_error = 0.0;
    double adjoint_scale = 1.0;
    double adjoint_lhs_worst = 0.0;
    double adjoint_rhs_worst = 0.0;
    for (int direction_id = 0; direction_id < 8; ++direction_id) {
        const double error = std::fabs(adjoint_lhs[direction_id] -
                                       adjoint_rhs[direction_id]);
        if (error >= adjoint_error) {
            adjoint_error = error;
            adjoint_lhs_worst = adjoint_lhs[direction_id];
            adjoint_rhs_worst = adjoint_rhs[direction_id];
        }
        adjoint_scale = std::max(adjoint_scale,
            std::max(std::fabs(adjoint_lhs[direction_id]),
                     std::fabs(adjoint_rhs[direction_id])));
    }
    const double adjoint_tolerance = 16384.0 *
        std::numeric_limits<double>::epsilon() * adjoint_scale;
    const MaxLocation global_adjoint = global_max_location(adjoint_error, rank);
    const double frozen_face_fraction = static_cast<double>(frozen_face_total) /
        std::max(1.0, static_cast<double>(candidate_face_total));
    const double x_frozen_face_fraction =
        static_cast<double>(x_frozen_face_total) /
        std::max(1.0, static_cast<double>(x_candidate_face_total));
    double coverage_counts[4] = {frozen_weight_cells, candidate_weight_cells,
        active_fully_covered_cells, active_candidate_cells};
    MPI_Allreduce(MPI_IN_PLACE, coverage_counts, 4, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    frozen_weight_cells = coverage_counts[0];
    candidate_weight_cells = coverage_counts[1];
    active_fully_covered_cells = coverage_counts[2];
    active_candidate_cells = coverage_counts[3];
    const double frozen_weight_fraction = frozen_weight_cells /
        std::max(1.0, candidate_weight_cells);
    const double active_fully_covered_fraction = active_fully_covered_cells /
        std::max(1.0, active_candidate_cells);
    const double kink_active_mass_fraction = kink_active_mass /
        std::max(1.0e-300, positive_mass);
    const double kink_pair_work_fraction = kink_pair_work /
        std::max(1.0e-300, all.abs_sum);
    const int nonsmooth_limiter_dominated = kink_active_mass_fraction > 1.0e-2 ||
        kink_pair_work_fraction > 5.0e-2;
    int epsilon_scale_stable = 1;
    double epsilon_scale_spread_max = 0.0;
    const double derivative_validation_tolerance = 16384.0 * derivative_eps;
    for (int direction_id = 0; direction_id < 8; ++direction_id) {
        bool all_scales_valid = true;
        double best_error = std::numeric_limits<double>::infinity();
        for (int scale_id = 0; scale_id < 3; ++scale_id) {
            all_scales_valid = all_scales_valid &&
                epsilon_scale_valid_count[direction_id][scale_id] > 0;
            best_error = std::min(best_error,
                epsilon_scale_error[direction_id][scale_id]);
        }
        if (!all_scales_valid) {
            epsilon_scale_stable = 0;
            continue;
        }
        const double spread = epsilon_scale_common_spread[direction_id];
        epsilon_scale_spread_max = std::max(epsilon_scale_spread_max, spread);
        if (!std::isfinite(spread) || !std::isfinite(best_error) ||
            spread > 2.0 * derivative_validation_tolerance ||
            best_error > derivative_validation_tolerance)
            epsilon_scale_stable = 0;
    }
    // Kink faces are not differentiable in the centered sense.  They count as
    // covered only when both finite one-sided responses were evaluated.
    const long long kink_face_count = u_branch_crossed_face_count +
        x_branch_crossed_face_count;
    const long long kink_one_sided_count = u_kink_one_sided_count +
        x_kink_one_sided_count;
    const int kink_one_sided_valid = kink_one_sided_nonfinite_count == 0 &&
        kink_one_sided_count == kink_face_count;
    const int coverage_valid = jacobian_audit_applicable &&
        kink_one_sided_valid;
    int local_failure_mask = 0;
    if (!direct_work_valid)
        local_failure_mask |= NON_FINITE_FAILED;
    if (jacobian_audit_applicable) {
        if (!std::isfinite(nonlinear_response_linf) ||
            !std::isfinite(x_nonlinear_response_linf) ||
            !std::isfinite(jacobian_action_linf) ||
            !std::isfinite(adjoint_error))
            local_failure_mask |= NON_FINITE_FAILED;
        if (nonlinear_response_linf > nonlinear_response_tolerance ||
            frozen_mapping_linf > nonlinear_response_tolerance)
            local_failure_mask |= U_REPLAY_FAILED;
        if (x_nonlinear_response_linf > x_nonlinear_response_tolerance)
            local_failure_mask |= X_REPLAY_FAILED;
        if (jacobian_action_linf > jacobian_action_tolerance ||
            !epsilon_scale_stable)
            local_failure_mask |= JACOBIAN_ACTION_FAILED;
        if (adjoint_error > adjoint_tolerance)
            local_failure_mask |= ADJOINT_FAILED;
        if (u_branch_crossed_face_count > 0 && !kink_one_sided_valid)
            local_failure_mask |= U_COVERAGE_FAILED;
        if (x_branch_crossed_face_count > 0 && !kink_one_sided_valid)
            local_failure_mask |= X_COVERAGE_FAILED;
        if (active_fully_covered_fraction < 0.999 && !kink_one_sided_valid)
            local_failure_mask |= ACTIVE_COVERAGE_FAILED;
    }
    int global_failure_mask = 0;
    MPI_Allreduce(&local_failure_mask, &global_failure_mask, 1, MPI_INT,
                  MPI_BOR, MPI_COMM_WORLD);
    const int failed_rank = first_failed_rank(local_failure_mask, rank);
    const int frozen_jacobian_valid = jacobian_audit_applicable &&
        global_failure_mask == 0;
    const int passes = direct_work_valid &&
        (!jacobian_audit_applicable || frozen_jacobian_valid);

    double pair_layers[4] = {};
    double fct_x_work_correction = 0.0;
    double fct_u_work_correction = 0.0;
    long long x_fct_active_face_count = 0;
    long long u_fct_active_face_count = 0;
    double min_alpha_x = 1.0;
    double min_alpha_u = 1.0;
    const double alpha_tolerance = 4096.0 *
        std::numeric_limits<double>::epsilon();
    if (direct_work_valid) for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double eface = fields.Ex_face[static_cast<size_t>(iface)];
        const size_t face = static_cast<size_t>(iface);
        const double layer_scale = dt * sg.dx * eface;
        const double r_low = layer_scale *
            (bundle.jn_low[face] - bundle.gstar_je_low[face]);
        const double r_high = layer_scale *
            (bundle.jn_high[face] - bundle.gstar_je_high[face]);
        const double r_final = layer_scale *
            (bundle.jn_final[face] - bundle.gstar_je_final[face]);
        pair_layers[0] += r_low;
        pair_layers[1] += r_high;
        pair_layers[2] += r_final - r_high;
        pair_layers[3] += r_final;
        fct_x_work_correction += layer_scale *
            (bundle.jn_final[face] - bundle.jn_high[face]);
        fct_u_work_correction -= layer_scale *
            (bundle.gstar_je_final[face] - bundle.gstar_je_high[face]);
        if (single_substep_flux_valid)
        for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
            const size_t id = xface_index(iface, j, k);
            const double low = bundle.fx_low[id];
            const double high = bundle.fx_high[id];
            const double final = bundle.fx_final[id];
            const double denominator = high - low;
            const double scale = std::max(1.0, std::max(std::fabs(low),
                std::max(std::fabs(high), std::fabs(final))));
            if (std::fabs(final - high) > alpha_tolerance * scale) {
                ++x_fct_active_face_count;
                if (std::fabs(denominator) > alpha_tolerance * scale)
                    min_alpha_x = std::min(min_alpha_x,
                        std::max(0.0, std::min(1.0, (final - low) / denominator)));
            }
        }
    }
    if (single_substep_flux_valid)
    for (int ix = 0; ix < sg.nx_local; ++ix)
        for (int jf = 0; jf <= Param::Nv; ++jf)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = uface_index(ix, jf, k);
                const double low = bundle.fu_low[id];
                const double high = bundle.fu_high[id];
                const double final = bundle.fu_final[id];
                const double denominator = high - low;
                const double scale = std::max(1.0, std::max(std::fabs(low),
                    std::max(std::fabs(high), std::fabs(final))));
                if (std::fabs(final - high) > alpha_tolerance * scale) {
                    ++u_fct_active_face_count;
                    if (std::fabs(denominator) > alpha_tolerance * scale)
                        min_alpha_u = std::min(min_alpha_u,
                            std::max(0.0, std::min(1.0, (final - low) / denominator)));
                }
            }
    MPI_Allreduce(MPI_IN_PLACE, pair_layers, 4, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &fct_x_work_correction, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &fct_u_work_correction, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &x_fct_active_face_count, 1, MPI_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &u_fct_active_face_count, 1, MPI_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &min_alpha_x, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &min_alpha_u, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);

    double work_n = 0.0;
    double work_e_cell = 0.0;
    if (direct_work_valid) for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double eface = fields.Ex_face[static_cast<size_t>(iface)];
        work_n += dt * sg.dx * eface * bundle.jn_high[static_cast<size_t>(iface)];
        work_e_cell += dt * sg.dx * eface *
            bundle.gstar_je_center[static_cast<size_t>(iface)];
    }
    MPI_Allreduce(MPI_IN_PLACE, &work_n, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &work_e_cell, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    std::vector<double> local_pair(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_pair_low(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_pair_high(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_pair_final(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_pair_x_fct(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_pair_u_fct(static_cast<size_t>(sg.nx_local), 0.0);
    if (direct_work_valid) {
        for (int ix = 0; ix < sg.nx_local; ++ix) {
            const double work = dt * sg.dx *
                fields.Ex_face[static_cast<size_t>(ix)];
            local_pair[static_cast<size_t>(ix)] = dt * sg.dx *
                fields.Ex_face[static_cast<size_t>(ix)] *
                (bundle.jn_high[static_cast<size_t>(ix)] -
                 bundle.gstar_je_center[static_cast<size_t>(ix)]);
            local_pair_low[static_cast<size_t>(ix)] = work *
                (bundle.jn_low[static_cast<size_t>(ix)] -
                 bundle.gstar_je_low[static_cast<size_t>(ix)]);
            local_pair_high[static_cast<size_t>(ix)] = work *
                (bundle.jn_high[static_cast<size_t>(ix)] -
                 bundle.gstar_je_high[static_cast<size_t>(ix)]);
            local_pair_final[static_cast<size_t>(ix)] = work *
                (bundle.jn_final[static_cast<size_t>(ix)] -
                 bundle.gstar_je_final[static_cast<size_t>(ix)]);
            local_pair_x_fct[static_cast<size_t>(ix)] = work *
                (bundle.jn_final[static_cast<size_t>(ix)] -
                 bundle.jn_high[static_cast<size_t>(ix)]);
            local_pair_u_fct[static_cast<size_t>(ix)] = -work *
                (bundle.gstar_je_final[static_cast<size_t>(ix)] -
                 bundle.gstar_je_high[static_cast<size_t>(ix)]);
        }
    }
    const std::vector<double> global_pair =
        BackgroundCouplingTest::gather_owned_faces(local_pair, sg, rank, size);
    const std::vector<double> global_pair_low =
        BackgroundCouplingTest::gather_owned_faces(local_pair_low, sg, rank, size);
    const std::vector<double> global_pair_high =
        BackgroundCouplingTest::gather_owned_faces(local_pair_high, sg, rank, size);
    const std::vector<double> global_pair_final =
        BackgroundCouplingTest::gather_owned_faces(local_pair_final, sg, rank, size);
    const std::vector<double> global_pair_x_fct =
        BackgroundCouplingTest::gather_owned_faces(local_pair_x_fct, sg, rank, size);
    const std::vector<double> global_pair_u_fct =
        BackgroundCouplingTest::gather_owned_faces(local_pair_u_fct, sg, rank, size);

    NormAccumulator spatial_pair_high;
    NormAccumulator spatial_pair_final;
    if (rank == 0 && direct_work_valid) {
        for (size_t ix = 0; ix < global_pair_high.size(); ++ix) {
            spatial_pair_high.add(global_pair_high[ix]);
            spatial_pair_final.add(global_pair_final[ix]);
        }
    }

    // A small, deterministic set of representative x locations is sufficient
    // for the expensive frozen-Jacobian audit.  Direct pair work remains
    // available at every face in by_x.dat.
    if (single_substep_flux_valid) {
    std::vector<double> local_e(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_de(static_cast<size_t>(sg.nx_local), 0.0);
    std::vector<double> local_density_gradient(static_cast<size_t>(sg.nx_local), 0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        local_e[static_cast<size_t>(ix)] = std::fabs(
            fields.Ex_face[static_cast<size_t>(ix)]);
        if (ix + 1 < sg.nx_local)
            local_de[static_cast<size_t>(ix)] = std::fabs(
                fields.Ex_face[static_cast<size_t>(ix + 1)] -
                fields.Ex_face[static_cast<size_t>(ix)]) / sg.dx;
        double left_density = 0.0;
        double right_density = 0.0;
        for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
            left_density += bundle.cu_reconstruction_mass[
                idx3(sg.nghost + ix - 1, j, k)];
            right_density += bundle.cu_reconstruction_mass[
                idx3(sg.nghost + ix + 1, j, k)];
        }
        local_density_gradient[static_cast<size_t>(ix)] =
            std::fabs(right_density - left_density) / (2.0 * sg.dx);
    }
    const std::vector<double> global_e =
        BackgroundCouplingTest::gather_owned_faces(local_e, sg, rank, size);
    const std::vector<double> global_de =
        BackgroundCouplingTest::gather_owned_faces(local_de, sg, rank, size);
    const std::vector<double> global_density_gradient =
        BackgroundCouplingTest::gather_owned_faces(local_density_gradient, sg, rank, size);
    if (rank == 0 && !global_pair.empty()) {
        const auto max_index = [](const std::vector<double>& values) {
            size_t best = 0;
            for (size_t i = 1; i < values.size(); ++i)
                if (values[i] > values[best]) best = i;
            return static_cast<int>(best);
        };
        const int proposed[7] = {1, std::max(1, Param::nx - 2),
            max_index(global_e), max_index(global_de),
            max_index(global_density_gradient),
            0,
            Param::nx / 2};
        // Replace the temporary zero vector selection with the max |pair| face.
        int pair_index = 0;
        for (size_t i = 1; i < global_pair.size(); ++i)
            if (std::fabs(global_pair[i]) > std::fabs(global_pair[pair_index]))
                pair_index = static_cast<int>(i);
        for (int q = 0; q < 7; ++q) {
            const int candidate = q == 5 ? pair_index : proposed[q];
            bool duplicate = false;
            for (int p = 0; p < audit_point_count; ++p)
                duplicate = duplicate || audit_points[p] == candidate;
            if (!duplicate) audit_points[audit_point_count++] = candidate;
        }
    }
    MPI_Bcast(&audit_point_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(audit_points, 7, MPI_INT, 0, MPI_COMM_WORLD);
    for (int point = 0; point < audit_point_count; ++point) {
        const int global_ix = audit_points[point];
        if (global_ix < sg.ix_start || global_ix >= sg.ix_start + sg.nx_local)
            continue;
        const int local_ix = global_ix - sg.ix_start;
        const int storage_ix = sg.nghost + local_ix;
        std::ostringstream point_directory;
        point_directory << output_directory << "/by_x_point/global_ix_"
                        << std::setw(4) << std::setfill('0') << global_ix;
        ensure_directory(point_directory.str());
        std::ofstream point_summary((point_directory.str() + "/summary.result").c_str());
        std::ofstream point_rows((point_directory.str() + "/by_u.dat").c_str());
        long long smooth_faces = 0, kink_faces = 0;
        long long one_sided_faces = 0;
        double max_jacobian_error = 0.0;
        double max_one_sided_error = 0.0;
        if (point_rows)
            point_rows << "j k branch_signature smooth frozen_direction "
                       << "finite_difference forward_difference backward_difference "
                       << "relative_error one_sided_relative_error\n";
        for (int k = 0; k < Param::Nmu; ++k) {
            std::vector<double> mass(static_cast<size_t>(Param::Nv), 0.0);
            double line_scale = 0.0;
            for (int j = 0; j < Param::Nv; ++j) {
                mass[static_cast<size_t>(j)] = bundle.cu_reconstruction_mass[
                    idx3(storage_ix, j, k)];
                line_scale = std::max(line_scale,
                    std::fabs(mass[static_cast<size_t>(j)]));
            }
            line_scale = std::max(line_scale, 1.0e-300);
            for (int jf = 1; jf < Param::Nv; ++jf) {
                const CenteredFaceAudit base = centered_coefficient_from_mass_line(
                    background.cgrid, mass, jf, k, sg.dx);
                std::vector<double> plus = mass;
                std::vector<double> minus = mass;
                const double eps = derivative_eps * line_scale;
                for (int j = 0; j < Param::Nv; ++j) {
                    const double direction = std::sin(0.371 * (j + 1) + 0.193 * (k + 1));
                    plus[static_cast<size_t>(j)] += eps * direction;
                    minus[static_cast<size_t>(j)] -= eps * direction;
                }
                const CenteredFaceAudit forward = centered_coefficient_from_mass_line(
                    background.cgrid, plus, jf, k, sg.dx);
                const CenteredFaceAudit backward = centered_coefficient_from_mass_line(
                    background.cgrid, minus, jf, k, sg.dx);
                const bool smooth = forward.branch_signature == base.branch_signature &&
                    backward.branch_signature == base.branch_signature;
                double frozen = 0.0;
                for (int p = 0; p < 4; ++p) {
                    const int j = base.stencil[p];
                    frozen += base.frozen.coefficient[p] /
                        background.cgrid.upar_widths[j] *
                        std::sin(0.371 * (j + 1) + 0.193 * (k + 1));
                }
                const double finite_difference = (forward.coefficient -
                    backward.coefficient) / (2.0 * eps);
                const double forward_difference =
                    (forward.production_coefficient -
                     base.production_coefficient) / eps;
                const double backward_difference =
                    (base.production_coefficient -
                     backward.production_coefficient) / eps;
                const double relative_error = smooth ? std::fabs(finite_difference - frozen) /
                    std::max(1.0e-300, std::max(std::fabs(finite_difference),
                                                   std::fabs(frozen))) : 0.0;
                const double one_sided_denominator = std::max(1.0e-300,
                    std::max(std::fabs(frozen),
                    std::max(std::fabs(forward_difference),
                             std::fabs(backward_difference))));
                const double one_sided_relative_error = std::max(
                    std::fabs(forward_difference - frozen),
                    std::fabs(backward_difference - frozen)) /
                    one_sided_denominator;
                if (smooth) {
                    ++smooth_faces;
                    max_jacobian_error = std::max(max_jacobian_error, relative_error);
                } else {
                    ++kink_faces;
                    if (std::isfinite(forward_difference) &&
                        std::isfinite(backward_difference)) {
                        ++one_sided_faces;
                        max_one_sided_error = std::max(max_one_sided_error,
                                                       one_sided_relative_error);
                    }
                }
                if (point_rows)
                    point_rows << jf << " " << k << " " << base.branch_signature << " "
                               << smooth << " " << frozen << " " << finite_difference
                               << " " << forward_difference << " "
                               << backward_difference << " " << relative_error << " "
                               << one_sided_relative_error << "\n";
            }
        }
        if (point_summary)
            point_summary << std::scientific << std::setprecision(17)
                << "global_ix=" << global_ix << "\nlocal_ix=" << local_ix << "\n"
                << "smooth_face_count=" << smooth_faces << "\n"
                << "kink_face_count=" << kink_faces << "\n"
                << "kink_one_sided_face_count=" << one_sided_faces << "\n"
                << "smooth_fraction=" << static_cast<double>(smooth_faces) /
                    std::max(1LL, smooth_faces + kink_faces) << "\n"
                << "jacobian_relative_error_linf=" << max_jacobian_error << "\n"
                << "kink_one_sided_relative_error_linf="
                << max_one_sided_error << "\n";
    }
    }
    BeamPIC hash_beam;
    const CheckpointStateHashes state_hashes = checkpoint_state_hashes(
        background, hash_beam, fields, sg, rank, size);
    const double unavailable = std::numeric_limits<double>::quiet_NaN();
    const double reported_w_l1 = single_substep_flux_valid ? w_l1 : unavailable;
    const double reported_w_l2 = single_substep_flux_valid ? w_l2 : unavailable;
    const double reported_w_linf = single_substep_flux_valid ? w_linf : unavailable;
    const double reported_pair_l1 = single_substep_flux_valid ?
        all.abs_sum : unavailable;
    const double reported_pair_l2 = single_substep_flux_valid ?
        std::sqrt(all.square_sum) : unavailable;
    const double reported_pair_linf = single_substep_flux_valid ?
        all.linf : unavailable;
    const double reported_pair_signed = single_substep_flux_valid ?
        all.signed_sum : unavailable;

    // Targeted periodic-seam audit.  Keep this deliberately narrow: six
    // physical faces around the seam plus the right alias of face zero.
    static const int seam_column_count = 13;
    const int seam_targets[7] = {0, 1, 2, Param::nx - 3,
                                 Param::nx - 2, Param::nx - 1, Param::nx};
    const int local_substep_count = static_cast<int>(
        bundle.coupling_substep_seam_audit.size());
    int substep_count_min = local_substep_count;
    int substep_count_max = local_substep_count;
    MPI_Allreduce(MPI_IN_PLACE, &substep_count_min, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &substep_count_max, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    int substep_seam_audit_valid = direct_work_valid &&
        substep_count_min == substep_count_max &&
        substep_count_min == bundle.substeps_used && substep_count_min > 0;
    const size_t local_face_count = static_cast<size_t>(sg.nx_local + 1);
    if (substep_seam_audit_valid) {
        for (int sub = 0; sub < local_substep_count; ++sub) {
            const VlasovAmpereMidpointSolver::CouplingSubstepSeamAudit& audit =
                bundle.coupling_substep_seam_audit[static_cast<size_t>(sub)];
            const std::vector<double>* arrays[] = {
                &audit.e_face_mid,
                &audit.jn_low_pre_sync, &audit.jn_high_pre_sync,
                &audit.jn_final_pre_sync, &audit.jn_low_post_sync,
                &audit.jn_high_post_sync, &audit.jn_final_post_sync,
                &audit.gstar_je_low_pre_sync, &audit.gstar_je_high_pre_sync,
                &audit.gstar_je_final_pre_sync,
                &audit.gstar_je_low_post_sync,
                &audit.gstar_je_high_post_sync,
                &audit.gstar_je_final_post_sync
            };
            for (int column = 0; column < seam_column_count; ++column)
                substep_seam_audit_valid = substep_seam_audit_valid &&
                    arrays[column]->size() >= local_face_count;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &substep_seam_audit_valid, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    std::vector<double> local_seam_values;
    std::vector<double> global_seam_values;
    if (substep_count_min > 0) {
        const size_t value_count = static_cast<size_t>(substep_count_min) *
            7 * seam_column_count;
        local_seam_values.assign(value_count, 0.0);
        global_seam_values.assign(value_count, 0.0);
    }
    if (substep_seam_audit_valid) {
        for (int sub = 0; sub < substep_count_min; ++sub) {
            const VlasovAmpereMidpointSolver::CouplingSubstepSeamAudit& audit =
                bundle.coupling_substep_seam_audit[static_cast<size_t>(sub)];
            const std::vector<double>* arrays[] = {
                &audit.e_face_mid,
                &audit.jn_low_pre_sync, &audit.jn_high_pre_sync,
                &audit.jn_final_pre_sync, &audit.jn_low_post_sync,
                &audit.jn_high_post_sync, &audit.jn_final_post_sync,
                &audit.gstar_je_low_pre_sync, &audit.gstar_je_high_pre_sync,
                &audit.gstar_je_final_pre_sync,
                &audit.gstar_je_low_post_sync,
                &audit.gstar_je_high_post_sync,
                &audit.gstar_je_final_post_sync
            };
            for (int target_id = 0; target_id < 7; ++target_id) {
                const int global_face = seam_targets[target_id];
                int local_face = -1;
                if (global_face < Param::nx && global_face >= sg.ix_start &&
                    global_face < sg.ix_start + sg.nx_local) {
                    local_face = global_face - sg.ix_start;
                } else if (global_face == Param::nx &&
                           sg.ix_start + sg.nx_local == Param::nx) {
                    local_face = sg.nx_local;
                }
                if (local_face < 0) continue;
                const size_t base = (static_cast<size_t>(sub) * 7 + target_id) *
                    seam_column_count;
                for (int column = 0; column < seam_column_count; ++column)
                    local_seam_values[base + column] =
                        (*arrays[column])[static_cast<size_t>(local_face)];
            }
        }
    }
    if (!local_seam_values.empty())
        MPI_Reduce(local_seam_values.data(), global_seam_values.data(),
                   static_cast<int>(local_seam_values.size()), MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);

    double seam_jn_sync_linf = 0.0;
    double seam_gstar_sync_linf = 0.0;
    double seam_owned_jn_sync_linf = 0.0;
    double seam_alias_jn_sync_linf = 0.0;
    double seam_owned_gstar_sync_linf = 0.0;
    double seam_alias_gstar_sync_linf = 0.0;
    double seam_pair_substep_jump_linf = 0.0;
    int seam_pair_jump_substep = -1;
    int seam_pair_jump_face = -1;
    if (rank == 0 && substep_seam_audit_valid) {
        ensure_directory(output_directory);
        std::ofstream seam((output_directory + "/seam_substeps.dat").c_str());
        if (seam) {
            seam << std::scientific << std::setprecision(17)
                 << "substep dt_substep global_face role E_mid "
                 << "JN_low_pre JN_low_post JN_high_pre JN_high_post "
                 << "JN_final_pre JN_final_post "
                 << "GstarJE_low_pre GstarJE_low_post "
                 << "GstarJE_high_pre GstarJE_high_post "
                 << "GstarJE_final_pre GstarJE_final_post "
                 << "pair_low_post pair_high_post pair_final_post "
                 << "x_FCT_correction_post u_FCT_correction_post "
                 << "JN_low_sync_delta JN_high_sync_delta JN_final_sync_delta "
                 << "GstarJE_low_sync_delta GstarJE_high_sync_delta "
                 << "GstarJE_final_sync_delta pair_final_work\n";
            for (int sub = 0; sub < substep_count_min; ++sub) {
                const double h = bundle.coupling_substep_seam_audit[
                    static_cast<size_t>(sub)].dt_substep;
                for (int target_id = 0; target_id < 7; ++target_id) {
                    const size_t base = (static_cast<size_t>(sub) * 7 +
                        target_id) * seam_column_count;
                    const double* v = global_seam_values.data() + base;
                    const double pair_low = v[4] - v[10];
                    const double pair_high = v[5] - v[11];
                    const double pair_final = v[6] - v[12];
                    const double x_fct = v[6] - v[5];
                    const double u_fct = v[12] - v[11];
                    const double jn_sync = std::max(
                        std::fabs(v[4] - v[1]),
                        std::max(std::fabs(v[5] - v[2]),
                                 std::fabs(v[6] - v[3])));
                    const double gstar_sync = std::max(
                        std::fabs(v[10] - v[7]),
                        std::max(std::fabs(v[11] - v[8]),
                                 std::fabs(v[12] - v[9])));
                    seam_jn_sync_linf = std::max(seam_jn_sync_linf, jn_sync);
                    seam_gstar_sync_linf = std::max(
                        seam_gstar_sync_linf, gstar_sync);
                    if (target_id == 6) {
                        seam_alias_jn_sync_linf = std::max(
                            seam_alias_jn_sync_linf, jn_sync);
                        seam_alias_gstar_sync_linf = std::max(
                            seam_alias_gstar_sync_linf, gstar_sync);
                    } else {
                        seam_owned_jn_sync_linf = std::max(
                            seam_owned_jn_sync_linf, jn_sync);
                        seam_owned_gstar_sync_linf = std::max(
                            seam_owned_gstar_sync_linf, gstar_sync);
                    }
                    if (sub > 0) {
                        const size_t previous = (static_cast<size_t>(sub - 1) *
                            7 + target_id) * seam_column_count;
                        const double* old = global_seam_values.data() + previous;
                        const double old_pair = old[6] - old[12];
                        const double jump = std::fabs(pair_final - old_pair);
                        if (jump > seam_pair_substep_jump_linf) {
                            seam_pair_substep_jump_linf = jump;
                            seam_pair_jump_substep = sub;
                            seam_pair_jump_face = seam_targets[target_id];
                        }
                    }
                    seam << sub << " " << h << " " << seam_targets[target_id]
                         << " " << (target_id == 6 ? "right_alias_of_0" : "owned")
                         << " " << v[0]
                         << " " << v[1] << " " << v[4]
                         << " " << v[2] << " " << v[5]
                         << " " << v[3] << " " << v[6]
                         << " " << v[7] << " " << v[10]
                         << " " << v[8] << " " << v[11]
                         << " " << v[9] << " " << v[12]
                         << " " << pair_low << " " << pair_high
                         << " " << pair_final << " " << x_fct
                         << " " << u_fct
                         << " " << v[4] - v[1]
                         << " " << v[5] - v[2]
                         << " " << v[6] - v[3]
                         << " " << v[10] - v[7]
                         << " " << v[11] - v[8]
                         << " " << v[12] - v[9]
                         << " " << h * sg.dx * v[0] * pair_final << "\n";
                }
            }
        }
    }

    // Resolve only the two dominant seam-adjacent physical faces by
    // u_parallel cell.  The six columns are J_N and G*J_E at low/high/final.
    static const int seam_upar_column_count = 6;
    const int seam_upar_targets[2] = {1, Param::nx - 1};
    int seam_upar_audit_valid = substep_seam_audit_valid;
    const size_t resolved_face_size = local_face_count * Param::Nv;
    if (seam_upar_audit_valid) {
        for (int sub = 0; sub < local_substep_count; ++sub) {
            const VlasovAmpereMidpointSolver::CouplingSubstepSeamAudit& audit =
                bundle.coupling_substep_seam_audit[static_cast<size_t>(sub)];
            const std::vector<double>* arrays[] = {
                &audit.jn_low_by_u_post_sync,
                &audit.jn_high_by_u_post_sync,
                &audit.jn_final_by_u_post_sync,
                &audit.gstar_je_low_by_u_post_sync,
                &audit.gstar_je_high_by_u_post_sync,
                &audit.gstar_je_final_by_u_post_sync
            };
            for (int column = 0; column < seam_upar_column_count; ++column)
                seam_upar_audit_valid = seam_upar_audit_valid &&
                    arrays[column]->size() >= resolved_face_size;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &seam_upar_audit_valid, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    std::vector<double> local_seam_upar_values;
    std::vector<double> global_seam_upar_values;
    if (substep_count_min > 0) {
        const size_t value_count = static_cast<size_t>(substep_count_min) *
            2 * Param::Nv * seam_upar_column_count;
        local_seam_upar_values.assign(value_count, 0.0);
        global_seam_upar_values.assign(value_count, 0.0);
    }
    if (seam_upar_audit_valid) {
        for (int sub = 0; sub < substep_count_min; ++sub) {
            const VlasovAmpereMidpointSolver::CouplingSubstepSeamAudit& audit =
                bundle.coupling_substep_seam_audit[static_cast<size_t>(sub)];
            const std::vector<double>* arrays[] = {
                &audit.jn_low_by_u_post_sync,
                &audit.jn_high_by_u_post_sync,
                &audit.jn_final_by_u_post_sync,
                &audit.gstar_je_low_by_u_post_sync,
                &audit.gstar_je_high_by_u_post_sync,
                &audit.gstar_je_final_by_u_post_sync
            };
            for (int target_id = 0; target_id < 2; ++target_id) {
                const int global_face = seam_upar_targets[target_id];
                if (global_face < sg.ix_start ||
                    global_face >= sg.ix_start + sg.nx_local) continue;
                const int local_face = global_face - sg.ix_start;
                for (int j = 0; j < Param::Nv; ++j) {
                    const size_t output_base =
                        ((static_cast<size_t>(sub) * 2 + target_id) *
                         Param::Nv + j) * seam_upar_column_count;
                    const size_t input = static_cast<size_t>(local_face) *
                        Param::Nv + j;
                    for (int column = 0; column < seam_upar_column_count;
                         ++column)
                        local_seam_upar_values[output_base + column] =
                            (*arrays[column])[input];
                }
            }
        }
    }
    if (!local_seam_upar_values.empty())
        MPI_Reduce(local_seam_upar_values.data(),
                   global_seam_upar_values.data(),
                   static_cast<int>(local_seam_upar_values.size()),
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double seam_upar_reconstruction_linf = 0.0;
    double seam_upar_reconstruction_scale = 0.0;
    double seam_upar_high_zone_work_l1[4] = {0.0, 0.0, 0.0, 0.0};
    double seam_upar_final_zone_work_l1[4] = {0.0, 0.0, 0.0, 0.0};
    int seam_upar_dominant_high_zone = -1;
    int seam_upar_dominant_final_zone = -1;
    if (rank == 0 && seam_upar_audit_valid) {
        const char* zone_names[4] = {"u_endpoint", "low_u",
                                     "thermal_body", "high_u"};
        std::ofstream raw((output_directory +
            "/seam_upar_substeps.dat").c_str());
        std::ofstream zones((output_directory +
            "/seam_upar_zones.dat").c_str());
        if (raw)
            raw << std::scientific << std::setprecision(17)
                << "substep dt_substep global_face j u_parallel zone "
                << "JN_low JN_high JN_final GstarJE_low GstarJE_high "
                << "GstarJE_final pair_low pair_high pair_final "
                << "x_FCT_correction u_FCT_correction pair_high_work "
                << "pair_final_work\n";
        if (zones)
            zones << std::scientific << std::setprecision(17)
                  << "substep dt_substep global_face zone "
                  << "JN_low JN_high JN_final GstarJE_low GstarJE_high "
                  << "GstarJE_final pair_low pair_high pair_final "
                  << "pair_high_work_signed pair_high_work_L1 "
                  << "pair_final_work_signed pair_final_work_L1\n";
        for (int sub = 0; sub < substep_count_min; ++sub) {
            const double h = bundle.coupling_substep_seam_audit[
                static_cast<size_t>(sub)].dt_substep;
            for (int target_id = 0; target_id < 2; ++target_id) {
                const int seam_target_id = target_id == 0 ? 1 : 5;
                const size_t seam_base =
                    (static_cast<size_t>(sub) * 7 + seam_target_id) *
                    seam_column_count;
                const double* seam_values =
                    global_seam_values.data() + seam_base;
                const double eface = seam_values[0];
                double reconstructed[6] = {0.0, 0.0, 0.0,
                                           0.0, 0.0, 0.0};
                double zone_sum[4][6] = {};
                double zone_high_signed[4] = {0.0, 0.0, 0.0, 0.0};
                double zone_high_l1[4] = {0.0, 0.0, 0.0, 0.0};
                double zone_final_signed[4] = {0.0, 0.0, 0.0, 0.0};
                double zone_final_l1[4] = {0.0, 0.0, 0.0, 0.0};
                for (int j = 0; j < Param::Nv; ++j) {
                    const size_t base =
                        ((static_cast<size_t>(sub) * 2 + target_id) *
                         Param::Nv + j) * seam_upar_column_count;
                    const double* v = global_seam_upar_values.data() + base;
                    const double u = background.cgrid.upar_cells[j];
                    int zone = 0;
                    if (j != 0 && j != Param::Nv - 1) {
                        const double abs_u = std::fabs(u);
                        zone = abs_u <= 0.1 ? 1 : (abs_u <= 1.0 ? 2 : 3);
                    }
                    for (int column = 0; column < 6; ++column) {
                        reconstructed[column] += v[column];
                        zone_sum[zone][column] += v[column];
                    }
                    const double pair_low = v[0] - v[3];
                    const double pair_high = v[1] - v[4];
                    const double pair_final = v[2] - v[5];
                    const double high_work = h * sg.dx * eface * pair_high;
                    const double final_work = h * sg.dx * eface * pair_final;
                    zone_high_signed[zone] += high_work;
                    zone_high_l1[zone] += std::fabs(high_work);
                    zone_final_signed[zone] += final_work;
                    zone_final_l1[zone] += std::fabs(final_work);
                    seam_upar_high_zone_work_l1[zone] += std::fabs(high_work);
                    seam_upar_final_zone_work_l1[zone] += std::fabs(final_work);
                    if (raw)
                        raw << sub << " " << h << " "
                            << seam_upar_targets[target_id] << " " << j
                            << " " << u << " " << zone_names[zone]
                            << " " << v[0] << " " << v[1] << " " << v[2]
                            << " " << v[3] << " " << v[4] << " " << v[5]
                            << " " << pair_low << " " << pair_high
                            << " " << pair_final << " " << v[2] - v[1]
                            << " " << v[5] - v[4] << " " << high_work
                            << " " << final_work << "\n";
                }
                const double aggregate[6] = {
                    seam_values[4], seam_values[5], seam_values[6],
                    seam_values[10], seam_values[11], seam_values[12]
                };
                for (int column = 0; column < 6; ++column)
                {
                    seam_upar_reconstruction_linf = std::max(
                        seam_upar_reconstruction_linf,
                        std::fabs(reconstructed[column] - aggregate[column]));
                    seam_upar_reconstruction_scale = std::max(
                        seam_upar_reconstruction_scale,
                        std::max(std::fabs(reconstructed[column]),
                                 std::fabs(aggregate[column])));
                }
                for (int zone = 0; zone < 4; ++zone) {
                    if (!zones) continue;
                    const double pair_low = zone_sum[zone][0] -
                        zone_sum[zone][3];
                    const double pair_high = zone_sum[zone][1] -
                        zone_sum[zone][4];
                    const double pair_final = zone_sum[zone][2] -
                        zone_sum[zone][5];
                    zones << sub << " " << h << " "
                          << seam_upar_targets[target_id] << " "
                          << zone_names[zone];
                    for (int column = 0; column < 6; ++column)
                        zones << " " << zone_sum[zone][column];
                    zones << " " << pair_low << " " << pair_high
                          << " " << pair_final
                          << " " << zone_high_signed[zone]
                          << " " << zone_high_l1[zone]
                          << " " << zone_final_signed[zone]
                          << " " << zone_final_l1[zone] << "\n";
                }
            }
        }
        for (int zone = 0; zone < 4; ++zone) {
            if (seam_upar_dominant_high_zone < 0 ||
                seam_upar_high_zone_work_l1[zone] >
                seam_upar_high_zone_work_l1[seam_upar_dominant_high_zone])
                seam_upar_dominant_high_zone = zone;
            if (seam_upar_dominant_final_zone < 0 ||
                seam_upar_final_zone_work_l1[zone] >
                seam_upar_final_zone_work_l1[seam_upar_dominant_final_zone])
                seam_upar_dominant_final_zone = zone;
        }
        const double reconstruction_tolerance = 1.0e-12 *
            std::max(1.0, seam_upar_reconstruction_scale);
        if (!std::isfinite(seam_upar_reconstruction_linf) ||
            seam_upar_reconstruction_linf > reconstruction_tolerance)
            seam_upar_audit_valid = 0;
    }
    MPI_Bcast(&seam_upar_audit_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
    const int audit_passes = passes &&
        (!checkpoint_loaded || seam_upar_audit_valid);

    if (rank == 0) {
        std::ostringstream stem;
#ifdef FP_UPAR_AUDIT_DETAILED
        stem << "background_coupling_upar_dual";
#else
        stem << "background_coupling_upar_transpose";
#endif
        stem << "_Nx" << Param::nx << "_Nu" << Param::Nv << "_Nuperp"
             << Param::Nmu;
        const std::string base = stem.str();
        ensure_directory(output_directory);
        std::ofstream out((output_directory + "/summary.result").c_str());
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_upar_dual_audit\n"
            << "production_kernel=1\nstate_source="
            << (checkpoint_loaded ? "midpoint_checkpoint" : "manufactured") << "\n"
            << "midpoint_audit_directory=" << midpoint_directory << "\n"
            << "midpoint_checkpoint_loaded=" << checkpoint_loaded << "\n"
            << "midpoint_checkpoint_error=" << checkpoint_error << "\n"
            << "state_hash_background=" << state_hashes.background << "\n"
            << "state_hash_field_faces=" << state_hashes.field_faces << "\n"
            << "state_hash_beam=" << state_hashes.beam << "\n"
            << "beam_enabled=0\nfct_enabled=" << options.fct_enabled << "\n"
            << "high_only=" << high_only << "\n"
            << "field_sign=" << (negative_field ? -1 : 1) << "\n"
            << "phase=" << phase << "\n"
            << "nx=" << Param::nx << "\nNu=" << Param::Nv
            << "\nNuperp=" << Param::Nmu << "\ndt=" << dt << "\n"
            << "operator_valid=" << direct_work_valid << "\n"
            << "direct_work_valid=" << direct_work_valid << "\n"
            << "single_substep_flux_valid=" << single_substep_flux_valid << "\n"
            << "velocity_resolved_audit_valid=" << single_substep_flux_valid << "\n"
            << "jacobian_audit_applicable=" << jacobian_audit_applicable << "\n"
            << "jacobian_audit_skipped=" << (!jacobian_audit_applicable) << "\n"
            << "state_advanced=" << bundle.state_advanced << "\n"
            << "operator_failed=" << bundle.operator_failed << "\n"
            << "outputs_finite=" << bundle.outputs_finite << "\n"
            << "substeps_used=" << bundle.substeps_used << "\n"
            << "substep_seam_audit_valid=" << substep_seam_audit_valid << "\n"
            << "substep_seam_audit_count=" << substep_count_min << "\n"
            << "seam_JN_sync_linf=" << seam_jn_sync_linf << "\n"
            << "seam_GstarJE_sync_linf=" << seam_gstar_sync_linf << "\n"
            << "seam_owned_JN_sync_linf=" << seam_owned_jn_sync_linf << "\n"
            << "seam_alias_JN_sync_linf=" << seam_alias_jn_sync_linf << "\n"
            << "seam_owned_GstarJE_sync_linf="
            << seam_owned_gstar_sync_linf << "\n"
            << "seam_alias_GstarJE_sync_linf="
            << seam_alias_gstar_sync_linf << "\n"
            << "seam_pair_substep_jump_linf=" << seam_pair_substep_jump_linf << "\n"
            << "seam_pair_jump_substep=" << seam_pair_jump_substep << "\n"
            << "seam_pair_jump_global_face=" << seam_pair_jump_face << "\n"
            << "seam_upar_audit_valid=" << seam_upar_audit_valid << "\n"
            << "seam_upar_reconstruction_linf="
            << seam_upar_reconstruction_linf << "\n"
            << "seam_upar_reconstruction_scale="
            << seam_upar_reconstruction_scale << "\n"
            << "seam_upar_reconstruction_relative="
            << seam_upar_reconstruction_linf /
                std::max(1.0, seam_upar_reconstruction_scale) << "\n"
            << "seam_upar_high_work_L1_u_endpoint="
            << seam_upar_high_zone_work_l1[0] << "\n"
            << "seam_upar_high_work_L1_low_u="
            << seam_upar_high_zone_work_l1[1] << "\n"
            << "seam_upar_high_work_L1_thermal_body="
            << seam_upar_high_zone_work_l1[2] << "\n"
            << "seam_upar_high_work_L1_high_u="
            << seam_upar_high_zone_work_l1[3] << "\n"
            << "seam_upar_final_work_L1_u_endpoint="
            << seam_upar_final_zone_work_l1[0] << "\n"
            << "seam_upar_final_work_L1_low_u="
            << seam_upar_final_zone_work_l1[1] << "\n"
            << "seam_upar_final_work_L1_thermal_body="
            << seam_upar_final_zone_work_l1[2] << "\n"
            << "seam_upar_final_work_L1_high_u="
            << seam_upar_final_zone_work_l1[3] << "\n"
            << "seam_upar_dominant_high_zone="
            << seam_upar_dominant_high_zone << "\n"
            << "seam_upar_dominant_final_zone="
            << seam_upar_dominant_final_zone << "\n"
            << "frozen_branch_jacobian=1\n"
            << "global_max_u_replay_error=" << global_u_replay.value << "\n"
            << "global_max_u_replay_rank=" << global_u_replay.rank << "\n"
            << "global_max_x_replay_error=" << global_x_replay.value << "\n"
            << "global_max_x_replay_rank=" << global_x_replay.rank << "\n"
            << "global_max_jacobian_error=" << global_jacobian.value << "\n"
            << "global_max_jacobian_rank=" << global_jacobian.rank << "\n"
            << "global_max_adjoint_error=" << global_adjoint.value << "\n"
            << "global_max_adjoint_rank=" << global_adjoint.rank << "\n"
            << "global_failure_mask=" << global_failure_mask << "\n"
            << "first_failed_rank=" << failed_rank << "\n"
            << "nonlinear_Au_response_linf=" << global_u_replay.value << "\n"
            << "nonlinear_Au_response_scale=" << global_u_scale.value << "\n"
            << "nonlinear_Au_response_tolerance=" << nonlinear_response_tolerance << "\n"
            << "frozen_branch_mapping_linf=" << global_frozen_mapping.value << "\n"
            << "x_nonlinear_response_linf=" << global_x_replay.value << "\n"
            << "x_nonlinear_response_scale=" << global_x_scale.value << "\n"
            << "x_nonlinear_response_tolerance="
            << x_nonlinear_response_tolerance << "\n"
            << "jacobian_action_linf=" << global_jacobian.value << "\n"
            << "jacobian_action_scale=" << global_jacobian_scale.value << "\n"
            << "jacobian_action_tolerance=" << jacobian_action_tolerance << "\n"
            << "adjoint_direction_count=8\n"
            << "adjoint_lhs=" << adjoint_lhs_worst << "\n"
            << "adjoint_rhs=" << adjoint_rhs_worst << "\n"
            << "adjoint_error=" << global_adjoint.value << "\n"
            << "adjoint_tolerance=" << adjoint_tolerance << "\n"
            << "frozen_face_count=" << frozen_face_total << "\n"
            << "candidate_face_count=" << candidate_face_total << "\n"
            << "frozen_face_fraction=" << frozen_face_fraction << "\n"
            << "frozen_zero_slope_count=" << frozen_zero_slope_count << "\n"
            << "x_frozen_face_count=" << x_frozen_face_total << "\n"
            << "x_candidate_face_count=" << x_candidate_face_total << "\n"
            << "x_frozen_face_fraction=" << x_frozen_face_fraction << "\n"
            << "frozen_weight_cell_count=" << frozen_weight_cells << "\n"
            << "candidate_weight_cell_count=" << candidate_weight_cells << "\n"
            << "frozen_weight_fraction=" << frozen_weight_fraction << "\n"
            << "active_candidate_weight_count=" << active_candidate_cells << "\n"
            << "active_fully_covered_count=" << active_fully_covered_cells << "\n"
            << "active_fully_covered_fraction="
            << active_fully_covered_fraction << "\n"
            << "coverage_valid=" << coverage_valid << "\n"
            << "epsilon_scale_stable=" << epsilon_scale_stable << "\n"
            << "epsilon_scale_spread_max=" << epsilon_scale_spread_max << "\n"
            << "derivative_validation_tolerance="
            << derivative_validation_tolerance << "\n"
            << "u_branch_crossed_face_count=" << u_branch_crossed_face_count << "\n"
            << "x_branch_crossed_face_count=" << x_branch_crossed_face_count << "\n"
            << "u_kink_one_sided_count=" << u_kink_one_sided_count << "\n"
            << "x_kink_one_sided_count=" << x_kink_one_sided_count << "\n"
            << "kink_one_sided_nonfinite_count="
            << kink_one_sided_nonfinite_count << "\n"
            << "kink_one_sided_valid=" << kink_one_sided_valid << "\n"
            << "audit_inconclusive=" <<
               (!direct_work_valid || !jacobian_audit_applicable ||
                !coverage_valid) << "\n"
            << "audit_local_ix=" << audit_ix << "\n"
            << "audit_global_ix=" << audit_global_ix << "\n"
            << "audit_ecell=" << audit_ecell << "\n"
            << "by_x_point_count=" << audit_point_count << "\n"
            << "frozen_jacobian_valid=" << frozen_jacobian_valid << "\n"
            << "wN_minus_wE_L1=" << reported_w_l1 << "\n"
            << "wN_minus_wE_L2=" << reported_w_l2 << "\n"
            << "wN_minus_wE_Linf=" << reported_w_linf << "\n"
            << "analytic_wN_minus_wE_L2=" << analytic_w_l2 << "\n"
            << "gN_all_L1=" << g_n_all.abs_sum << "\n"
            << "gN_all_L2=" << std::sqrt(g_n_all.square_sum) << "\n"
            << "gN_all_Linf=" << g_n_all.linf << "\n"
            << "gE_all_L1=" << g_e_all.abs_sum << "\n"
            << "gE_all_L2=" << std::sqrt(g_e_all.square_sum) << "\n"
            << "gE_all_Linf=" << g_e_all.linf << "\n"
            << "gR_all_L1=" << g_r_all.abs_sum << "\n"
            << "gR_all_L2=" << std::sqrt(g_r_all.square_sum) << "\n"
            << "gR_all_Linf=" << g_r_all.linf << "\n"
            << "gR_mass_weighted_L2=" << gr_mass_l2 << "\n"
            << "gN_mass_weighted_L2=" << gn_mass_l2 << "\n"
            << "gE_mass_weighted_L2=" << ge_mass_l2 << "\n"
            << "gR_relative=" << relative_gradient << "\n"
            << "gR_mass_weighted_L1=" << weighted_gr_l1 << "\n"
            << "gN_mass_weighted_L1=" << weighted_gn_l1 << "\n"
            << "gE_mass_weighted_L1=" << weighted_ge_l1 << "\n"
            << "gR_mass_weighted_signed=" << weighted_gr_signed << "\n"
            << "gR_work_relative=" << relative_work_gradient << "\n"
            << "W_N_signed=" << work_n << "\n"
            << "W_E_cell_signed=" << work_e_cell << "\n"
            << "B_u_K_signed=" << local_boundary_work << "\n"
            << "W_E_signed=" << work_e_cell + local_boundary_work << "\n"
            << "W_N_minus_W_E=" << work_n - work_e_cell - local_boundary_work
            << "\n"
            << "production_r_couple=" << bundle.r_couple << "\n"
            << "R_pair_low=" << pair_layers[0] << "\n"
            << "R_pair_high=" << pair_layers[1] << "\n"
            << "R_pair_FCT_correction=" << pair_layers[2] << "\n"
            << "R_pair_final=" << pair_layers[3] << "\n"
            << "R_pair_final_minus_production=" << pair_layers[3] - bundle.r_couple << "\n"
            << "spatial_pair_high_L1=" << spatial_pair_high.abs_sum << "\n"
            << "spatial_pair_high_L2=" << std::sqrt(spatial_pair_high.square_sum) << "\n"
            << "spatial_pair_high_Linf=" << spatial_pair_high.linf << "\n"
            << "spatial_pair_high_signed=" << spatial_pair_high.signed_sum << "\n"
            << "spatial_pair_final_L1=" << spatial_pair_final.abs_sum << "\n"
            << "spatial_pair_final_L2=" << std::sqrt(spatial_pair_final.square_sum) << "\n"
            << "spatial_pair_final_Linf=" << spatial_pair_final.linf << "\n"
            << "spatial_pair_final_signed=" << spatial_pair_final.signed_sum << "\n"
            << "fct_work_correction=" << pair_layers[2] << "\n"
            << "fct_x_work_correction=" << fct_x_work_correction << "\n"
            << "fct_u_work_correction=" << fct_u_work_correction << "\n"
            << "fct_face_activity_valid=" << single_substep_flux_valid << "\n"
            << "x_fct_active_face_count=" << x_fct_active_face_count << "\n"
            << "u_fct_active_face_count=" << u_fct_active_face_count << "\n"
            << "min_alpha_x=" << min_alpha_x << "\n"
            << "min_alpha_u=" << min_alpha_u << "\n"
            << "kink_active_mass_fraction=" << kink_active_mass_fraction << "\n"
            << "kink_pair_work_fraction=" << kink_pair_work_fraction << "\n"
            << "nonsmooth_limiter_dominated=" << nonsmooth_limiter_dominated << "\n"
            << "kink_one_sided_relative_error_min="
            << kink_one_sided_relative_error_min << "\n"
            << "kink_one_sided_relative_error_max="
            << kink_one_sided_relative_error_max << "\n"
            << "support_1e-8_gR_L1=" << support[0].abs_sum << "\n"
            << "support_1e-8_gR_L2=" << std::sqrt(support[0].square_sum) << "\n"
            << "support_1e-8_gR_Linf=" << support[0].linf << "\n"
            << "support_1e-12_gR_L1=" << support[1].abs_sum << "\n"
            << "support_1e-12_gR_L2=" << std::sqrt(support[1].square_sum) << "\n"
            << "support_1e-12_gR_Linf=" << support[1].linf << "\n"
            << "support_1e-16_gR_L1=" << support[2].abs_sum << "\n"
            << "support_1e-16_gR_L2=" << std::sqrt(support[2].square_sum) << "\n"
            << "support_1e-16_gR_Linf=" << support[2].linf << "\n"
            << "pair_work_L1=" << reported_pair_l1 << "\n"
            << "pair_work_L2=" << reported_pair_l2 << "\n"
            << "pair_work_Linf=" << reported_pair_linf << "\n"
            << "pair_work_signed=" << reported_pair_signed << "\n"
            << "tail_weight_difference_L1=" << tail_weight_difference << "\n"
            << "tail_weight_relative_L1=" << tail_weight_difference /
                   std::max(1.0, tail_weight_total) << "\n"
            << "u_boundary_energy_work=" << local_boundary_work << "\n"
            << "positive_negative_u_symmetry_L1=" << symmetry.abs_sum << "\n"
            << "positive_negative_u_symmetry_L2=" << std::sqrt(symmetry.square_sum)
            << "\npositive_negative_u_symmetry_Linf=" << symmetry.linf << "\n"
            << "worst_j=" << worst_j << "\nworst_k=" << worst_k << "\n"
            << "worst_u_parallel=" << (worst_j >= 0 ?
                background.cgrid.upar_cells[worst_j] : 0.0) << "\n"
            << "worst_u_perp=" << (worst_k >= 0 ?
                background.cgrid.uperp_cells[worst_k] : 0.0) << "\n"
            << "worst_delta_K_lower=" << (worst_j >= 0 && worst_k >= 0 ?
                delta_k_lower[idx2(worst_j, worst_k)] : 0.0) << "\n"
            << "worst_delta_K_upper=" << (worst_j >= 0 && worst_k >= 0 ?
                delta_k_upper[idx2(worst_j, worst_k)] : 0.0) << "\n"
            << "worst_u_grid_distance=" << (worst_j > 0 && worst_j < Param::Nv ?
                background.cgrid.upar_center_distances[worst_j] : 0.0) << "\n"
            << "worst_cylindrical_weight=" << (worst_k >= 0 ?
                background.cgrid.uperp_ring_areas[worst_k] : 0.0) << "\n";
        for (int point = 0; point < audit_point_count; ++point)
            log << "by_x_point_" << point << "_global_ix="
                << audit_points[point] << "\n";
        for (int z = 0; z < 5; ++z) {
            log << zone_name(z) << "_cell_count=";
            double count = 0.0;
            double mass = 0.0;
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = idx2(j, k);
                if (audit_zone(background.cgrid, j, k, mass_jk[id], maximum_mass) == z) {
                    count += 1.0;
                    mass += std::max(0.0, mass_jk[id]);
                }
            }
            log << count << "\n" << zone_name(z) << "_positive_mass_fraction="
                << mass / std::max(1.0e-300, positive_mass) << "\n";
            log << zone_name(z) << "_pair_work_signed=" << zone[z].signed_sum << "\n"
                << zone_name(z) << "_pair_work_L1=" << zone[z].abs_sum << "\n"
                << zone_name(z) << "_pair_work_L2=" << std::sqrt(zone[z].square_sum) << "\n"
                << zone_name(z) << "_pair_work_Linf=" << zone[z].linf << "\n"
                << zone_name(z) << "_pair_work_L1_fraction=" << zone[z].abs_sum /
                    std::max(1.0e-300, all.abs_sum) << "\n";
        }
        log << "passes=" << audit_passes << "\n";

        std::ofstream rows((output_directory + "/by_u.dat").c_str());
        if (rows) {
            rows << "# velocity_resolved_audit_valid="
                 << single_substep_flux_valid << "\n";
            rows << "j k u_parallel u_perp wN_production wN_analytic wE "
                 << "wN_minus_wE "
                 << "JN_signed JN_abs JE_lower_signed JE_upper_signed "
                 << "JE_lower_abs JE_upper_abs pair_work_signed pair_work_abs "
                 << "delta_K_lower delta_K_upper upar_width cylindrical_weight "
                 << "u_frozen_adjacent_face_count u_candidate_adjacent_face_count "
                 << "x_frozen_face_count x_candidate_face_count fully_covered\n";
            if (single_substep_flux_valid)
            for (int j = 0; j < Param::Nv; ++j) for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = idx2(j, k);
                rows << j << " " << k << " " << background.cgrid.upar_cells[j]
                     << " " << background.cgrid.uperp_cells[k] << " "
                     << w_n[id] << " " << w_n_analytic[id] << " " << w_e[id]
                     << " " << w_n[id] - w_e[id]
                     << " " << jn_signed[id]
                     << " " << jn_abs[id]
                     << " " << je_signed[static_cast<size_t>(j) * Param::Nmu + k]
                     << " " << je_signed[static_cast<size_t>(j + 1) * Param::Nmu + k]
                     << " " << je_abs[static_cast<size_t>(j) * Param::Nmu + k]
                     << " " << je_abs[static_cast<size_t>(j + 1) * Param::Nmu + k]
                     << " " << residual_jk[id] << " " << residual_jk_abs[id]
                     << " " << delta_k_lower[id] << " " << delta_k_upper[id]
                     << " " << background.cgrid.upar_widths[j] << " "
                     << background.cgrid.uperp_ring_areas[k] << " "
                     << frozen_face_count[id] << " "
                     << candidate_adjacent_face_count[id] << " "
                     << x_frozen_face_count[id] << " "
                     << x_candidate_face_count[id] << " "
                     << ((candidate_adjacent_face_count[id] > 0 &&
                          frozen_face_count[id] == candidate_adjacent_face_count[id] &&
                          x_candidate_face_count[id] > 0 &&
                          x_frozen_face_count[id] == x_candidate_face_count[id]) ? 1 : 0)
                     << "\n";
            }
        }
        std::ofstream coverage((output_directory + "/coverage.result").c_str());
        if (coverage) {
            coverage << std::scientific << std::setprecision(17)
                     << "jacobian_audit_applicable="
                     << jacobian_audit_applicable << "\n"
                     << "u_candidate_face_count=" << candidate_face_total << "\n"
                     << "u_frozen_face_count=" << frozen_face_total << "\n"
                     << "u_frozen_face_fraction=" << frozen_face_fraction << "\n"
                     << "x_candidate_face_count=" << x_candidate_face_total << "\n"
                     << "x_frozen_face_count=" << x_frozen_face_total << "\n"
                     << "x_frozen_face_fraction=" << x_frozen_face_fraction << "\n"
                     << "active_candidate_weight_count=" << active_candidate_cells << "\n"
                     << "active_fully_covered_count=" << active_fully_covered_cells << "\n"
                     << "active_fully_covered_fraction=" << active_fully_covered_fraction << "\n"
                     << "u_branch_crossed_face_count=" << u_branch_crossed_face_count << "\n"
                     << "x_branch_crossed_face_count=" << x_branch_crossed_face_count << "\n"
                     << "branch_crossed_face_count=" <<
                        u_branch_crossed_face_count + x_branch_crossed_face_count << "\n"
                     << "u_kink_one_sided_count=" << u_kink_one_sided_count << "\n"
                     << "x_kink_one_sided_count=" << x_kink_one_sided_count << "\n"
                     << "kink_one_sided_nonfinite_count="
                     << kink_one_sided_nonfinite_count << "\n"
                     << "kink_one_sided_valid=" << kink_one_sided_valid << "\n"
                     << "coverage_valid=" << coverage_valid << "\n";
        }
        std::ofstream directional((output_directory +
                                   "/directional_derivative.result").c_str());
        if (directional) {
            directional << std::scientific << std::setprecision(17)
                        << "jacobian_audit_applicable="
                        << jacobian_audit_applicable << "\n"
                        << "direction_count=8\n"
                        << "epsilon_relative=" << derivative_eps << "\n"
                        << "epsilon_scale_count=3\n"
                        << "derivative_method=frozen_branch_analytic\n"
                        << "frozen_mapping_linf=" << frozen_mapping_linf << "\n"
                        << "jacobian_action_linf=" << jacobian_action_linf << "\n"
                        << "jacobian_action_tolerance=" << jacobian_action_tolerance << "\n"
                        << "adjoint_error=" << global_adjoint.value << "\n"
                        << "adjoint_tolerance=" << adjoint_tolerance << "\n";
            for (int direction_id = 0; direction_id < 8; ++direction_id) {
                for (int scale_id = 0; scale_id < 3; ++scale_id) {
                    const double error = epsilon_scale_error[direction_id][scale_id];
                    directional << "direction_" << direction_id << "_scale_"
                                << (scale_id == 0 ? "025" :
                                    (scale_id == 1 ? "100" : "400"))
                                << "_error=" << error << "\n"
                                << "direction_" << direction_id << "_scale_"
                                << (scale_id == 0 ? "025" :
                                    (scale_id == 1 ? "100" : "400"))
                                << "_valid_count="
                                << epsilon_scale_valid_count[direction_id][scale_id]
                                << "\n";
                }
                directional << "direction_" << direction_id << "_scale_spread="
                            << epsilon_scale_common_spread[direction_id] << "\n";
            }
            directional << "epsilon_scale_spread_max="
                        << epsilon_scale_spread_max << "\n"
                        << "epsilon_scale_stable=" << epsilon_scale_stable << "\n"
                        << "derivative_validation_tolerance="
                        << derivative_validation_tolerance << "\n"
                        << "kink_one_sided_relative_error_min="
                        << kink_one_sided_relative_error_min << "\n"
                        << "kink_one_sided_relative_error_max="
                        << kink_one_sided_relative_error_max << "\n";
            for (int scale_id = 0; scale_id < 3; ++scale_id)
                directional << "epsilon_scale_"
                            << (scale_id == 0 ? "025" :
                                (scale_id == 1 ? "100" : "400"))
                            << "_branch_crossed="
                            << epsilon_scale_branch_crossed[scale_id] << "\n";
        }
        std::ofstream zones((output_directory + "/by_zone.dat").c_str());
        if (zones) {
            zones << "# velocity_resolved_audit_valid="
                  << single_substep_flux_valid << "\n";
            zones << "zone pair_work_signed pair_work_L1 pair_work_L2 pair_work_Linf\n";
            if (single_substep_flux_valid) for (int z = 0; z < 5; ++z)
                zones << zone_name(z) << " " << zone[z].signed_sum << " "
                      << zone[z].abs_sum << " " << std::sqrt(zone[z].square_sum)
                      << " " << zone[z].linf << "\n";
        }
        std::ofstream by_x((output_directory + "/by_x.dat").c_str());
        if (by_x) {
            by_x << "global_x_cell legacy_high_center_pair_work "
                 << "R_pair_low R_pair_high R_pair_final "
                 << "delta_R_x_FCT delta_R_u_FCT pair_work_L1\n";
            for (size_t ix = 0; ix < global_pair.size(); ++ix)
                by_x << ix << " " << global_pair[ix] << " "
                     << global_pair_low[ix] << " " << global_pair_high[ix]
                     << " " << global_pair_final[ix] << " "
                     << global_pair_x_fct[ix] << " "
                     << global_pair_u_fct[ix] << " "
                     << std::fabs(global_pair_final[ix]) << "\n";
        }
        // Only the transpose refinement producers own flat compatibility
        // aliases.  Detailed/checkpoint cases use their unique output
        // directories and must not overwrite one another.
#ifndef FP_UPAR_AUDIT_DETAILED
        const std::string legacy_path = "output/" + base + ".result";
        const std::string legacy_tmp = legacy_path + ".tmp";
        std::ofstream legacy(legacy_tmp.c_str());
        if (legacy) {
            legacy << "passes=" << audit_passes
                   << "\ncoverage_valid=" << coverage_valid
                   << "\noperator_valid=" << direct_work_valid
                   << "\ndirect_work_valid=" << direct_work_valid
                   << "\nsingle_substep_flux_valid=" << single_substep_flux_valid
                   << "\noutputs_finite=" << bundle.outputs_finite
                   << "\nfrozen_jacobian_valid=" << frozen_jacobian_valid
                   << "\nepsilon_scale_stable=" << epsilon_scale_stable
                   << "\nkink_one_sided_valid=" << kink_one_sided_valid
                   << "\nnx=" << Param::nx << "\nNu=" << Param::Nv
                   << "\nNuperp=" << Param::Nmu << "\n"
                   << "wN_minus_wE_L2=" << w_l2 << "\n"
                   << "pair_work_L2=" << std::sqrt(all.square_sum) << "\n"
                   << "tail_weight_relative_L1=" << tail_weight_difference /
                      std::max(1.0, tail_weight_total) << "\n"
                   << "frozen_weight_fraction=" << frozen_weight_fraction << "\n"
                   << "active_fully_covered_fraction="
                   << active_fully_covered_fraction << "\n"
                   << "gR_work_relative=" << relative_work_gradient << "\n"
                   << "kink_active_mass_fraction="
                   << kink_active_mass_fraction << "\n"
                   << "kink_pair_work_fraction="
                   << kink_pair_work_fraction << "\n";
            legacy.close();
            std::remove(legacy_path.c_str());
            if (std::rename(legacy_tmp.c_str(), legacy_path.c_str()) != 0)
                std::cerr << "WARNING: unable to publish audit result "
                          << legacy_path << "\n";
        }
#endif
        std::cout << "background_coupling_upar_dual_audit result="
                  << output_directory << "/summary.result passes="
                  << audit_passes << "\n";
    }

    MPI_Finalize();
    return audit_passes ? 0 : 1;
}
