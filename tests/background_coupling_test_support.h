#ifndef BACKGROUND_COUPLING_TEST_SUPPORT_H
#define BACKGROUND_COUPLING_TEST_SUPPORT_H

#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <string>
#include <vector>

namespace BackgroundCouplingTest {

struct Norms {
    double l1;
    double l2;
    double linf;
};

struct BundleOptions {
    const Species* guess_np1;
    const EMFields* fields_end_guess;
    bool controlled_fct_injection;
    bool fct_enabled;
    bool allow_finite_negative_debt;
    bool energy_consistent_x_high_velocity;
    const VlasovAmpereMidpointSolver::CouplingRegionLayout* coupling_layout;

    BundleOptions()
        : guess_np1(0), fields_end_guess(0),
          controlled_fct_injection(false), fct_enabled(true),
          allow_finite_negative_debt(false),
          energy_consistent_x_high_velocity(false), coupling_layout(0) {}
};

inline double wave_number(const SpatialGrid& sg, int mode)
{
    return 2.0 * Const::pi * static_cast<double>(mode) /
        (sg.dx * sg.nx_global);
}

inline void initialize_periodic_state(Species& background, EMFields& fields,
                                      const SpatialGrid& sg, int rank,
                                      int size, double density_amplitude,
                                      double velocity_amplitude,
                                      double field_amplitude,
                                      int mode, double phase)
{
    background.init("background_coupling_manufactured",
                    SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me,
                    Param::dens, Param::temperature_e, false, sg);
    background.initialize_maxwellian();
    const double kx = wave_number(sg, mode);
    const int ng = sg.nghost;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double x = sg.x(ng + ix);
        const double c = std::cos(kx * x + phase);
        const double s = std::sin(kx * x + phase);
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t q = idx2(j, k);
                const double factor = 1.0 + density_amplitude * c +
                    velocity_amplitude *
                    (background.cgrid.vx[q] / Const::c) * s;
                background.f[idx3(ng + ix, j, k)] *= factor;
            }
        }
    }
    VlasovAmpereMidpointSolver ghost_sync;
    ghost_sync.synchronize_background_ghosts(background, sg, rank, size);

    fields.init(sg);
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double xf = sg.x_min +
            (static_cast<double>(sg.ix_start + iface)) * sg.dx;
        fields.Ex_face[static_cast<size_t>(iface)] =
            field_amplitude * std::cos(kx * xf + phase);
    }
    fields.sync_cell_ex_from_faces(rank, size);
}

inline VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle
evaluate_bundle(const Species& background, const EMFields& fields,
                const SpatialGrid& sg, int rank, int size, double dt,
                const BundleOptions& options = BundleOptions())
{
    BeamPIC beam;
    VlasovAmpereMidpointSolver solver;
    solver.set_step_diagnostics_enabled(false);
    solver.set_beam_enabled(false);
    solver.set_low_order_only(false);
    solver.set_nonuniform_high_order_enabled(true);
    solver.set_fct_enabled(options.fct_enabled);
    solver.set_allow_finite_negative_debt_for_test(
        options.allow_finite_negative_debt);
    solver.set_energy_consistent_x_high_velocity_for_test(
        options.energy_consistent_x_high_velocity);
    solver.set_fct_activation_audit_enabled(options.controlled_fct_injection);
    solver.set_controlled_fct_flux_injection_enabled(
        options.controlled_fct_injection);
    const std::vector<double> no_beam(static_cast<size_t>(sg.nx_local + 1),
                                      0.0);
    VlasovAmpereMidpointSolver::CouplingRegionLayout default_layout = {};
    default_layout.beam_front_ix = -1;
    default_layout.wave_core_end_m = Param::Lx;
    const VlasovAmpereMidpointSolver::CouplingRegionLayout& layout =
        options.coupling_layout ? *options.coupling_layout : default_layout;
    return solver.evaluate_background_coupling_flux_bundle(
        background, beam, fields,
        options.guess_np1 ? *options.guess_np1 : background,
        options.fields_end_guess ? *options.fields_end_guess : fields,
        no_beam, layout, sg, dt, 0.0, rank, size);
}

inline Norms face_difference_norms(const std::vector<double>& left,
                                   const std::vector<double>& right,
                                   const SpatialGrid& sg)
{
    double local_l1 = 0.0;
    double local_l2 = 0.0;
    double local_linf = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double d = left[static_cast<size_t>(iface)] -
            right[static_cast<size_t>(iface)];
        local_l1 += std::fabs(d);
        local_l2 += d * d;
        local_linf = std::max(local_linf, std::fabs(d));
    }
    double values[3] = {local_l1, local_l2, local_linf};
    double global[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(values, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(values + 2, global + 2, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    Norms norms = {};
    norms.l1 = global[0] / static_cast<double>(sg.nx_global);
    norms.l2 = std::sqrt(global[1] / static_cast<double>(sg.nx_global));
    norms.linf = global[2];
    return norms;
}

inline double global_vector_difference_linf(const std::vector<double>& left,
                                            const std::vector<double>& right)
{
    double local = 0.0;
    const size_t count = std::min(left.size(), right.size());
    for (size_t p = 0; p < count; ++p)
        local = std::max(local, std::fabs(left[p] - right[p]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

inline double global_vector_abs_linf(const std::vector<double>& values)
{
    double local = 0.0;
    for (size_t p = 0; p < values.size(); ++p)
        local = std::max(local, std::fabs(values[p]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

inline double global_signed_face_work_residual(
    const std::vector<double>& jn, const std::vector<double>& gstar_je,
    const EMFields& fields, const SpatialGrid& sg, double dt)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface)
        local += dt * sg.dx * fields.Ex_face[static_cast<size_t>(iface)] *
            (jn[static_cast<size_t>(iface)] -
             gstar_je[static_cast<size_t>(iface)]);
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

inline double stable_dt(const SpatialGrid& sg)
{
    return 0.025 * sg.dx / Const::c;
}

inline Norms face_work_residual_norms(
    const std::vector<double>& jn, const std::vector<double>& gstar_je,
    const EMFields& fields, const SpatialGrid& sg, double dt)
{
    std::vector<double> local(static_cast<size_t>(sg.nx_local), 0.0);
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        local[static_cast<size_t>(iface)] = dt * sg.dx *
            fields.Ex_face[static_cast<size_t>(iface)] *
            (jn[static_cast<size_t>(iface)] -
             gstar_je[static_cast<size_t>(iface)]);
    }
    return face_difference_norms(local,
                                 std::vector<double>(local.size(), 0.0), sg);
}

inline std::vector<double> gather_owned_faces(const std::vector<double>& local,
                                              const SpatialGrid& sg,
                                              int rank, int size)
{
    std::vector<int> counts;
    std::vector<int> offsets;
    if (rank == 0) {
        counts.assign(static_cast<size_t>(size), 0);
        offsets.assign(static_cast<size_t>(size), 0);
    }
    const int nlocal = sg.nx_local;
    MPI_Gather(&nlocal, 1, MPI_INT, rank == 0 ? counts.data() : 0, 1,
               MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<double> gathered;
    if (rank == 0) {
        int total = 0;
        for (int r = 0; r < size; ++r) {
            offsets[static_cast<size_t>(r)] = total;
            total += counts[static_cast<size_t>(r)];
        }
        gathered.assign(static_cast<size_t>(total), 0.0);
    }
    MPI_Gatherv(local.data(), nlocal, MPI_DOUBLE,
                rank == 0 ? gathered.data() : 0,
                rank == 0 ? counts.data() : 0,
                rank == 0 ? offsets.data() : 0, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    return gathered;
}

inline void make_distinct_endpoint_guess(const Species& start,
                                         const EMFields& fields_start,
                                         Species& guess, EMFields& fields_end,
                                         const SpatialGrid& sg, int rank,
                                         int size)
{
    guess = start;
    fields_end = fields_start;
    const int ng = sg.nghost;
    const double kx = wave_number(sg, 5);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double x = sg.x(ng + ix);
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = idx3(ng + ix, j, k);
                guess.f[id] *= 1.0 + 0.08 * std::sin(kx * x) +
                    0.03 * start.cgrid.vx[idx2(j, k)] / Const::c;
            }
        }
    }
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double x = sg.x_min + (sg.ix_start + iface) * sg.dx;
        fields_end.Ex_face[static_cast<size_t>(iface)] +=
            0.35 * fields_start.Ex_face[static_cast<size_t>(iface)] +
            250.0 * std::sin(kx * x);
    }
    VlasovAmpereMidpointSolver ghost_sync;
    ghost_sync.synchronize_background_ghosts(guess, sg, rank, size);
    fields_end.sync_cell_ex_from_faces(rank, size);
}

inline void write_key_value(std::ostream& out, const char* key, double value)
{
    out << key << "=" << std::scientific << std::setprecision(17) << value
        << "\n";
}

} // namespace BackgroundCouplingTest

#endif
