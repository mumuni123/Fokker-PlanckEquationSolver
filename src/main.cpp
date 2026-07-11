#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"
#include "beam_pic.h"
#include "vlasov_ampere_midpoint.h"
#include "config.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <omp.h>
#include <vector>

double compute_dt(const Species& electron, const SpatialGrid& sg)
{
    double dt_min = 0.4 * sg.dx / Const::c;
    double e_est = Param::densb * Const::qe * sg.dx / Const::eps0;
    double udot_est = std::abs(electron.charge) * e_est /
                    (electron.mass * Const::c);
    if (udot_est > 1.0e-30) {
        double min_du = electron.vgrid.v_widths.empty()
                      ? electron.vgrid.dv
                      : electron.vgrid.v_widths[0];
        for (size_t iv = 1; iv < electron.vgrid.v_widths.size(); ++iv) {
            min_du = std::min(min_du, electron.vgrid.v_widths[iv]);
        }
        dt_min = std::min(dt_min, 0.25 * min_du / udot_est);
    }
    dt_min *= Param::dt_multiplier;
    return std::min(dt_min, 0.01 * Const::femto);
}

const char* field_solver_name()
{
    return "face-centered periodic Vlasov-Ampere update";
}

std::vector<double> build_local_ion_density_profile(const SpatialGrid& sg)
{
    return std::vector<double>(static_cast<size_t>(sg.nx_local), Param::dens);
}

void sync_moments_and_charge(Species& electrons,
                             const BeamPIC& beam,
                             EMFields& fields,
                             const std::vector<double>& ion_density_profile,
                             bool& moments_current)
{
    if (!moments_current) {
        electrons.compute_moments();
        moments_current = true;
    }
    fields.set_charge_density(electrons, beam.density, ion_density_profile);
}

struct BackgroundCurrentDiagnostics {
    double residual_if_charge;
    double residual_if_ampere;
    double e_dot_j_charge;
    double e_dot_j_energy;
    double e_dot_j_ampere;
    double max_abs_j_charge;
    double max_abs_j_energy;
    double max_abs_j_ampere;
    double max_abs_j_charge_minus_ampere;
    double max_abs_j_energy_minus_ampere;
};

struct BoundaryCoreDiagnostics {
    double boundary_width;
    double neg_mass_boundary;
    double neg_mass_core;
    long long neg_cell_count_boundary;
    long long neg_cell_count_core;
    double min_f_boundary;
    double min_f_core;
    double max_abs_J_bkg_boundary;
    double max_abs_J_bkg_core;
    double rms_J_bkg_boundary;
    double rms_J_bkg_core;
    double mean_abs_J_bkg_boundary;
    double mean_abs_J_bkg_core;
    double max_abs_Ex_boundary;
    double max_abs_Ex_core;
    double W_bkg_E_boundary;
    double W_bkg_E_core;
    double anomaly_inward_E;
    double anomaly_inward_J;
    double anomaly_inward_n;
    double anomaly_inward_Mneg;
};

struct FNegativitySnapshotDiagnostics {
    double min_f;
    double neg_ratio_max;
    double neg_mass_total;
    long long neg_cell_count;
    int x_worst;
    int u_worst;
    int mu_worst;
};

void reset_background_current_diagnostics(BackgroundCurrentDiagnostics& diag)
{
    diag.residual_if_charge = 0.0;
    diag.residual_if_ampere = 0.0;
    diag.e_dot_j_charge = 0.0;
    diag.e_dot_j_energy = 0.0;
    diag.e_dot_j_ampere = 0.0;
    diag.max_abs_j_charge = 0.0;
    diag.max_abs_j_energy = 0.0;
    diag.max_abs_j_ampere = 0.0;
    diag.max_abs_j_charge_minus_ampere = 0.0;
    diag.max_abs_j_energy_minus_ampere = 0.0;
}

void reset_boundary_core_diagnostics(BoundaryCoreDiagnostics& diag,
                                     double boundary_width)
{
    diag.boundary_width = boundary_width;
    diag.neg_mass_boundary = 0.0;
    diag.neg_mass_core = 0.0;
    diag.neg_cell_count_boundary = 0;
    diag.neg_cell_count_core = 0;
    diag.min_f_boundary = 0.0;
    diag.min_f_core = 0.0;
    diag.max_abs_J_bkg_boundary = 0.0;
    diag.max_abs_J_bkg_core = 0.0;
    diag.rms_J_bkg_boundary = 0.0;
    diag.rms_J_bkg_core = 0.0;
    diag.mean_abs_J_bkg_boundary = 0.0;
    diag.mean_abs_J_bkg_core = 0.0;
    diag.max_abs_Ex_boundary = 0.0;
    diag.max_abs_Ex_core = 0.0;
    diag.W_bkg_E_boundary = 0.0;
    diag.W_bkg_E_core = 0.0;
    diag.anomaly_inward_E = 0.0;
    diag.anomaly_inward_J = 0.0;
    diag.anomaly_inward_n = 0.0;
    diag.anomaly_inward_Mneg = 0.0;
}

void reset_f_negativity_snapshot_diagnostics(
    FNegativitySnapshotDiagnostics& diag)
{
    diag.min_f = 0.0;
    diag.neg_ratio_max = 0.0;
    diag.neg_mass_total = 0.0;
    diag.neg_cell_count = 0;
    diag.x_worst = -1;
    diag.u_worst = -1;
    diag.mu_worst = -1;
}

void compute_f_negativity_snapshot_diagnostics(
    const Species& electrons,
    const SpatialGrid& sg,
    FNegativitySnapshotDiagnostics& diag)
{
    reset_f_negativity_snapshot_diagnostics(diag);

    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    double local_min_f = std::numeric_limits<double>::infinity();
    double local_max_positive_f = 0.0;
    double local_neg_mass = 0.0;
    long long local_neg_count = 0;
    int local_min_indices[3] = { -1, -1, -1 };

    for (int ix = 0; ix < nxl; ++ix) {
        const size_t xbase =
            static_cast<size_t>(ng + ix) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double f = electrons.f[row + static_cast<size_t>(imu)];
                if (f < local_min_f) {
                    local_min_f = f;
                    local_min_indices[0] = sg.ix_start + ix;
                    local_min_indices[1] = iv;
                    local_min_indices[2] = imu;
                }
                if (f > 0.0) {
                    local_max_positive_f = std::max(local_max_positive_f, f);
                }
                if (f < 0.0) {
                    local_neg_mass +=
                        (-f) * electrons.vgrid.moment_weight[iv] * sg.dx;
                    ++local_neg_count;
                }
            }
        }
    }

    int mpi_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    struct { double val; int rank; } local_min_loc, global_min_loc;
    local_min_loc.val = local_min_f;
    local_min_loc.rank = mpi_rank;
    global_min_loc.val = std::numeric_limits<double>::infinity();
    global_min_loc.rank = -1;
    MPI_Allreduce(&local_min_loc, &global_min_loc, 1, MPI_DOUBLE_INT,
                  MPI_MINLOC, MPI_COMM_WORLD);
    double global_max_positive_f = local_max_positive_f;
    MPI_Allreduce(MPI_IN_PLACE, &global_max_positive_f, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    double global_neg_mass = 0.0;
    MPI_Allreduce(&local_neg_mass, &global_neg_mass, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    long long global_neg_count = 0;
    MPI_Allreduce(&local_neg_count, &global_neg_count, 1,
                  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    const double scale =
        std::max(global_max_positive_f,
                 std::numeric_limits<double>::min());
    if (global_min_loc.rank >= 0) {
        MPI_Bcast(local_min_indices, 3, MPI_INT, global_min_loc.rank,
                  MPI_COMM_WORLD);
    }

    diag.min_f = std::isfinite(global_min_loc.val) ? global_min_loc.val : 0.0;
    diag.neg_ratio_max = (diag.min_f < 0.0) ? -diag.min_f / scale : 0.0;
    diag.neg_mass_total = global_neg_mass;
    diag.neg_cell_count = global_neg_count;
    diag.x_worst = (diag.min_f < 0.0) ? local_min_indices[0] : -1;
    diag.u_worst = (diag.min_f < 0.0) ? local_min_indices[1] : -1;
    diag.mu_worst = (diag.min_f < 0.0) ? local_min_indices[2] : -1;
}

void compute_boundary_core_diagnostics(
    const Species& electrons,
    const EMFields& fields,
    const SpatialGrid& sg,
    double dt,
    double boundary_width,
    int mpi_rank,
    int mpi_size,
    BoundaryCoreDiagnostics& diag)
{
    reset_boundary_core_diagnostics(diag, boundary_width);
    const int ng = sg.nghost;
    const int nxl = sg.nx_local;
    double local_neg_mass_boundary = 0.0;
    double local_neg_mass_core = 0.0;
    long long local_neg_count_boundary = 0;
    long long local_neg_count_core = 0;
    double local_min_f_boundary = std::numeric_limits<double>::infinity();
    double local_min_f_core = std::numeric_limits<double>::infinity();

    for (int ix = 0; ix < nxl; ++ix) {
        const double x_cell =
            (static_cast<double>(sg.ix_start + ix) + 0.5) * sg.dx;
        const bool boundary =
            (x_cell < boundary_width) || (x_cell > Param::Lx - boundary_width);
        const size_t xbase =
            static_cast<size_t>(ng + ix) * Param::Nvmu;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double shell = electrons.vgrid.moment_weight[iv];
            const size_t row = xbase + static_cast<size_t>(iv) * Param::Nmu;
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const double f = electrons.f[row + static_cast<size_t>(imu)];
                if (boundary) {
                    local_min_f_boundary =
                        std::min(local_min_f_boundary, f);
                } else {
                    local_min_f_core = std::min(local_min_f_core, f);
                }
                if (f < 0.0) {
                    const double neg_mass = -f * shell * sg.dx;
                    if (boundary) {
                        local_neg_mass_boundary += neg_mass;
                        ++local_neg_count_boundary;
                    } else {
                        local_neg_mass_core += neg_mass;
                        ++local_neg_count_core;
                    }
                }
            }
        }
    }

    double local_max_j_boundary = 0.0;
    double local_max_j_core = 0.0;
    double local_sum_abs_j_boundary = 0.0;
    double local_sum_abs_j_core = 0.0;
    double local_sum_j2_boundary = 0.0;
    double local_sum_j2_core = 0.0;
    double local_max_ex_boundary = 0.0;
    double local_max_ex_core = 0.0;
    double local_work_boundary = 0.0;
    double local_work_core = 0.0;
    long long local_face_count_boundary = 0;
    long long local_face_count_core = 0;

    const int face_count =
        std::min(nxl, static_cast<int>(std::min(electrons.current_face_x.size(),
                                                fields.Ex_face.size())));
    for (int iface = 0; iface < face_count; ++iface) {
        const double x_face =
            static_cast<double>(sg.ix_start + iface) * sg.dx;
        const bool boundary =
            (x_face < boundary_width) || (x_face > Param::Lx - boundary_width);
        const double j = electrons.current_face_x[static_cast<size_t>(iface)];
        const double ex = fields.Ex_face[static_cast<size_t>(iface)];
        const double abs_j = std::fabs(j);
        const double abs_ex = std::fabs(ex);
        const double work = -j * ex * sg.dx * dt;
        if (boundary) {
            local_max_j_boundary = std::max(local_max_j_boundary, abs_j);
            local_sum_abs_j_boundary += abs_j;
            local_sum_j2_boundary += j * j;
            local_max_ex_boundary = std::max(local_max_ex_boundary, abs_ex);
            local_work_boundary += work;
            ++local_face_count_boundary;
        } else {
            local_max_j_core = std::max(local_max_j_core, abs_j);
            local_sum_abs_j_core += abs_j;
            local_sum_j2_core += j * j;
            local_max_ex_core = std::max(local_max_ex_core, abs_ex);
            local_work_core += work;
            ++local_face_count_core;
        }
    }

    double sums[10] = {
        local_neg_mass_boundary,
        local_neg_mass_core,
        local_sum_abs_j_boundary,
        local_sum_abs_j_core,
        local_sum_j2_boundary,
        local_sum_j2_core,
        local_work_boundary,
        local_work_core,
        0.0,
        0.0
    };
    MPI_Allreduce(MPI_IN_PLACE, sums, 10, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    double mins[2] = { local_min_f_boundary, local_min_f_core };
    MPI_Allreduce(MPI_IN_PLACE, mins, 2, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);

    double maxes[4] = {
        local_max_j_boundary,
        local_max_j_core,
        local_max_ex_boundary,
        local_max_ex_core
    };
    MPI_Allreduce(MPI_IN_PLACE, maxes, 4, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    long long counts[4] = {
        local_neg_count_boundary,
        local_neg_count_core,
        local_face_count_boundary,
        local_face_count_core
    };
    MPI_Allreduce(MPI_IN_PLACE, counts, 4, MPI_LONG_LONG_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    // Dynamic inward-propagation diagnostic.  Thresholds are tied to the
    // current core envelope, so the 0.2 um belt remains a test isolation
    // width rather than becoming a moving physical boundary.
    double local_extent[4] = {0.0,0.0,0.0,0.0};
    const double e_trigger=1.25*std::max(1.0,maxes[3]);
    const double j_trigger=1.25*std::max(std::fabs(Param::jb),maxes[1]);
    for (int ix=0;ix<nxl;++ix) {
        const double x=(static_cast<double>(sg.ix_start+ix)+0.5)*sg.dx;
        const double inward=std::min(x,Param::Lx-x);
        const size_t xb=static_cast<size_t>(ng+ix)*Param::Nvmu;
        double n=0.0,nneg=0.0;
        for (int iv=0;iv<Param::Nv;++iv) {
            const double w=electrons.vgrid.moment_weight[iv];
            for (int imu=0;imu<Param::Nmu;++imu) {
                const double f=electrons.f[xb+static_cast<size_t>(iv)*Param::Nmu+imu];
                n += f*w;
                if (f<0.0) nneg += (-f)*w;
            }
        }
        const double ex=(static_cast<size_t>(ix)<fields.Ex.size())
            ? std::fabs(fields.Ex[static_cast<size_t>(ix)]) : 0.0;
        const double j=(static_cast<size_t>(ix)<electrons.current_x.size())
            ? std::fabs(electrons.current_x[static_cast<size_t>(ix)]) : 0.0;
        if (ex>e_trigger) local_extent[0]=std::max(local_extent[0],inward);
        if (j>j_trigger) local_extent[1]=std::max(local_extent[1],inward);
        if (std::fabs(n-Param::dens)>1.0e-5*Param::dens)
            local_extent[2]=std::max(local_extent[2],inward);
        if (nneg>1.0e-8*std::max(1.0,std::fabs(n)))
            local_extent[3]=std::max(local_extent[3],inward);
    }
    MPI_Allreduce(MPI_IN_PLACE,local_extent,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    (void)mpi_rank;
    (void)mpi_size;
    diag.neg_mass_boundary = sums[0];
    diag.neg_mass_core = sums[1];
    diag.neg_cell_count_boundary = counts[0];
    diag.neg_cell_count_core = counts[1];
    diag.min_f_boundary =
        std::isfinite(mins[0]) ? mins[0] : 0.0;
    diag.min_f_core =
        std::isfinite(mins[1]) ? mins[1] : 0.0;
    diag.max_abs_J_bkg_boundary = maxes[0];
    diag.max_abs_J_bkg_core = maxes[1];
    diag.mean_abs_J_bkg_boundary =
        (counts[2] > 0) ? sums[2] / static_cast<double>(counts[2]) : 0.0;
    diag.mean_abs_J_bkg_core =
        (counts[3] > 0) ? sums[3] / static_cast<double>(counts[3]) : 0.0;
    diag.rms_J_bkg_boundary =
        (counts[2] > 0) ? std::sqrt(sums[4] / static_cast<double>(counts[2]))
                        : 0.0;
    diag.rms_J_bkg_core =
        (counts[3] > 0) ? std::sqrt(sums[5] / static_cast<double>(counts[3]))
                        : 0.0;
    diag.max_abs_Ex_boundary = maxes[2];
    diag.max_abs_Ex_core = maxes[3];
    diag.W_bkg_E_boundary = sums[6];
    diag.W_bkg_E_core = sums[7];
    diag.anomaly_inward_E=local_extent[0];
    diag.anomaly_inward_J=local_extent[1];
    diag.anomaly_inward_n=local_extent[2];
    diag.anomaly_inward_Mneg=local_extent[3];
}

const char* classify_boundary_core_state(double neg_mass_core_fraction,
                                         double neg_cell_core_fraction,
                                         double limiter_core_fraction,
                                         double j_core_to_boundary_ratio)
{
    if (neg_mass_core_fraction > 0.2 ||
        limiter_core_fraction > 0.2 ||
        j_core_to_boundary_ratio > 0.5) {
        return "core_contaminated";
    }
    if (neg_mass_core_fraction < 0.05 &&
        limiter_core_fraction < 0.05 &&
        j_core_to_boundary_ratio < 0.2) {
        return "boundary_only";
    }
    return "transition";
}

void write_snapshot(Diagnostics& diag,
                    double time,
                    const Species& bkg_e,
                    const BeamPIC& beam,
                    EMFields& fields,
                    const std::vector<double>& ion_density_profile,
                    const SpatialGrid& sgrid,
                    int mpi_rank,
                    int mpi_size,
                    bool write_full_fe,
                    const std::vector<double>* bkg_energy_current_face = 0,
                    const std::vector<double>* bkg_ampere_current_face = 0)
{
    fields.compute_potential(mpi_rank, mpi_size);
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_current_density(time, bkg_e, beam, sgrid, mpi_rank, mpi_size,
                               bkg_energy_current_face,
                               bkg_ampere_current_face);
    diag.write_density_profile(time, bkg_e, beam.density, ion_density_profile,
                               sgrid, mpi_rank, mpi_size);
    diag.write_px_distribution(time, bkg_e, mpi_rank, mpi_size);
    if (write_full_fe) {
        diag.write_electron_distribution(time, bkg_e, sgrid, mpi_rank);
    }

    if (mpi_rank == 0) {
        printf("  >> Snapshot %d written at t = %.4f fs\n",
               diag.snapshot_count, time / Const::femto);
    }
    diag.advance_snapshot();
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    RuntimeConfig config = load_runtime_config();

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Background-electron VFP + fixed ions + PIC beam solver\n");
        printf("  Spherical electron momentum grid: (u, mu), u = p / (m c)\n");
        printf("============================================================\n");
        printf("MPI ranks: %d\n", mpi_size);
        #pragma omp parallel
        {
            #pragma omp master
            printf("OpenMP threads per rank: %d\n", omp_get_num_threads());
        }
        printf("Spatial grid: nx = %d, dx = %.3e m, Lx = %.3e m\n",
               Param::nx, Param::dx, Param::Lx);
        printf("Density profile: uniform plasma over full domain, n0 = %.3e /m^3\n",
               Param::dens);
        printf("Electron momentum grid: Nu_parallel x Nu_perp = %d x %d, conservative cylindrical mass M\n",
               Param::Nv, Param::Nmu);
        printf("Electron momentum domain: |u_parallel| <= %.3f, 0 <= u_perp <= %.3f, vx = c u_parallel / gamma\n",
               Param::momentum_umax, Param::momentum_umax);
        printf("Spatial boundary: periodic in x for background electrons and electrostatic field; beam is open\n");
        printf("Electrostatic update: face-centered dE/dt = -J_total/eps0; zero mode evolves explicitly\n");
        printf("Field solver: %s\n", field_solver_name());
        printf("Fixed ions: uniform Z*n_i = %.3e /m^3\n", Param::dens);
        printf("Background electrons: periodic Vlasov transport, T_e = %.1f eV\n",
               Param::temperature_e / Const::eV);
        printf("PIC beam: gamma*beta = %.2f, beta = %.4f, n_b = %.3e /m^3\n",
               Param::gambetab, Param::betab, Param::densb);
        printf("Beam source: sampled left-boundary crossings at x = 0\n");
        printf("Beam injection: remaining-time push with charge-conserving face current\n");
        printf("Beam charge compensation source: OFF; background density perturbations are produced only by Poisson/Vlasov dynamics\n");
        printf("Beam boundary: particles crossing either global edge leave the domain\n");
        printf("Background boundary: periodic ghost cells\n");
        printf("Beam macro weight: %.6e particles/m^2\n", Param::beam_macro_weight);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
        printf("Debug diagnostics: %s\n",
               config.enable_debug_diagnostics ? "ON" : "OFF");
#else
        printf("Debug diagnostics: compile-time disabled\n");
#endif
        printf("Full fe distribution output: %s\n",
               config.enable_full_fe_output ? "ON" : "OFF");
        printf("Step diagnostics: %s\n",
               config.enable_step_diagnostics ? "ON" : "OFF");
        if (config.enable_step_diagnostics) {
            printf("Step diagnostics interval: every %d steps\n",
                   config.step_diagnostics_interval);
        }
        printf("Progress trace: %s\n",
               config.enable_progress_trace ? "ON" : "OFF");
        printf("------------------------------------------------------------\n");
    }

    SpatialGrid sgrid;
    sgrid.init(mpi_rank, mpi_size);
    std::vector<double> ion_density_profile = build_local_ion_density_profile(sgrid);

    Species bkg_e;
    bkg_e.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON,
               -Const::qe, Const::me,
               Param::dens, Param::temperature_e, false, sgrid);
    bkg_e.initialize_maxwellian();

    BeamPIC beam;
    beam.init(sgrid);
    beam.deposit_density(sgrid, mpi_rank, mpi_size);

    EMFields fields;
    fields.init(sgrid);

    VlasovAmpereMidpointSolver midpoint_solver;
    midpoint_solver.set_step_diagnostics_enabled(false);
    CollisionOperator collision;
    Diagnostics diag;
    diag.init("output", mpi_rank, config.enable_debug_diagnostics,
              config.enable_step_diagnostics);

    double dt = compute_dt(bkg_e, sgrid);
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    int nsteps = static_cast<int>(std::ceil(Param::t_end / dt));

    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    double cumulative_collision_energy_delta = 0.0;
    bool moments_current = false;
    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    fields.solve_poisson(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(0.0, 0, bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    std::vector<double> latest_bkg_energy_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    std::vector<double> latest_bkg_ampere_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    write_snapshot(diag, 0.0, bkg_e, beam, fields, ion_density_profile,
                   sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                   &latest_bkg_energy_current_face,
                   &latest_bkg_ampere_current_face);

    double next_snapshot = Param::dt_snapshot;
    int stdout_freq = 1000;
    int last_snapshot_step = 0;
    double cumulative_bkg_energy_residual = 0.0;
    std::ofstream bkg_energy_monitor;
    std::ofstream boundary_core_monitor;
    std::ofstream x_low_monitor;
    std::ofstream x_final_monitor;
    std::ofstream mu_final_monitor;
    std::ofstream u_flux_audit;
    std::ofstream u_origin_audit;
    std::ofstream u_low_failure_audit;
    std::ofstream mu_low_failure_audit;
    std::ofstream negative_debt_guard_monitor;
    std::ofstream trial_boundary_debt_monitor;
    if (mpi_rank == 0) {
        bkg_energy_monitor.open("output/bkg_energy_monitor.dat");
        bkg_energy_monitor
            << "# step  time[fs]  x_limiter_active_fraction  "
            << "x_limiter_min_alpha  dKE_bkg_plus_W_bkg_E[J/m2]  "
            << "relative_residual_step  cumulative_relative_residual  "
            << "x_limiter_active_fraction_core  "
            << "x_limiter_active_fraction_boundary  "
            << "x_limiter_min_alpha_core  "
            << "x_limiter_min_alpha_boundary  "
            << "x_negative_mass_before_repair[m^-2]  "
            << "x_mass_added_by_positivity_repair[m^-2]  "
            << "floor_repair_mass[m^-2]  "
            << "floor_repair_energy[J/m2]  "
            << "floor_repair_core_fraction  "
            << "max_abs_J_bkg_charge[A/m2]  "
            << "max_abs_J_bkg_energy_diagnostic[A/m2]  "
            << "max_abs_J_bkg_ampere[A/m2]  "
            << "max_abs_J_bkg_charge_minus_ampere[A/m2]  "
            << "max_abs_J_bkg_energy_minus_ampere[A/m2]  "
            << "E_dot_J_bkg_charge[W/m2]  "
            << "E_dot_J_bkg_energy_diagnostic[W/m2]  "
            << "E_dot_J_bkg_ampere[W/m2]  "
            << "residual_if_charge_current[J/m2]  "
            << "residual_if_ampere_current[J/m2]  "
            << "positivity_mass_defect[m^-2]  "
            << "positivity_energy_defect[J/m2]  "
            << "u_limiter_mass_delta[m^-2]  "
            << "u_limiter_px_delta[kg/m/s/m2]  "
            << "u_limiter_energy_delta[J/m2]  "
            << "x_limiter_mass_delta[m^-2]  "
            << "x_limiter_energy_delta[J/m2]  "
            << "mu_low_u_alpha_min  "
            << "mu_low_u_limiter_active_fraction  "
            << "mu_low_u_energy_delta[J/m2]  "
            << "mu_low_u_alpha_min_boundary  "
            << "mu_low_u_alpha_min_core  "
            << "mu_low_u_limiter_active_fraction_boundary  "
            << "mu_low_u_limiter_active_fraction_core  "
            << "mu_low_u_energy_delta_boundary[J/m2]  "
            << "mu_low_u_energy_delta_core[J/m2]  "
            << "mu_low_u_u_eff0  "
            << "mu_low_u_moment_weight0  "
            << "mu_low_u_mu_flux_scale0  "
            << "mu_low_u_half_dt_inv_shell0  "
            << "mu_low_u_dimless_scale0  "
            << "mu_low_u_endpoint_flux_max  "
            << "remap_active_fraction  "
            << "remap_cell_count  "
            << "low_u_subcycle_active_fraction  "
            << "low_u_average_subcycles  "
            << "low_u_max_subcycles  "
            << "u_force_alpha_min  u_force_alpha_active_frac  "
            << "coupled_iter  coupled_residual_E  "
            << "coupled_residual_J_bkg  coupled_residual_J_beam  "
            << "coupled_residual_bkg_mass  "
            << "coupled_residual_beam_continuity\n";
        bkg_energy_monitor << std::scientific;

        u_flux_audit.open("output/u_flux_audit.dat");
        u_flux_audit
            << "# step  time[fs]  valid  rank  ix_global  iv  imu  "
            << "severity  f0  f_low  f_high  alpha  "
            << "du_div_low  du_div_high  du_div_final  updated  "
            << "u_final_xleft_lower  u_final_xleft_upper  "
            << "u_final_xright_lower  u_final_xright_upper\n";
        u_flux_audit << std::scientific << std::setprecision(8);

        u_origin_audit.open("output/u_origin_audit.dat");
        u_origin_audit
            << "# step time[fs] rank ix_global iv imu region "
            << "f_raw f_after_x f_u_low f_u_high f_u_final "
            << "u_lower_flux u_upper_flux "
            << "u_lower_donor_iv u_upper_donor_iv "
            << "u_lower_donor_f u_upper_donor_f "
            << "C_u alpha_cell alpha_face "
            << "line_positive_peak line_negative_mass allowed_debt "
            << "failure_kind\n";
        u_origin_audit << std::scientific << std::setprecision(8);

        u_low_failure_audit.open("output/u_low_failure_audit.dat");
        u_low_failure_audit
            << "# step  time[fs]  rank  ix_global  iv  imu  region  "
            << "severity  f_input  f_after_x  dx_div  dmu_div_used  "
            << "du_div_low  f_low  "
            << "left_lower_u_low  left_upper_u_low  "
            << "right_lower_u_low  right_upper_u_low  "
            << "left_lower_scale  left_upper_scale  "
            << "right_lower_scale  right_upper_scale  "
            << "left_lower_donor_iv  left_upper_donor_iv  "
            << "right_lower_donor_iv  right_upper_donor_iv  "
            << "left_lower_donor_f  left_upper_donor_f  "
            << "right_lower_donor_f  right_upper_donor_f  "
            << "lower_characteristic  upper_characteristic  "
            << "moment_weight  cell_budget  "
            << "low_order_failed_count\n";
        u_low_failure_audit << std::scientific << std::setprecision(8);

        mu_low_failure_audit.open("output/mu_low_failure_audit.dat");
        mu_low_failure_audit
            << "# step  time[fs]  rank  ix_global  iv  imu  region  "
            << "severity  f_before_mu  f_after_x  dx_div  dmu_div_low  "
            << "du_div_used  f_mu_low  "
            << "left_lower_mu_low  left_upper_mu_low  "
            << "right_lower_mu_low  right_upper_mu_low  "
            << "left_lower_mu_dot  left_upper_mu_dot  "
            << "right_lower_mu_dot  right_upper_mu_dot  "
            << "left_lower_donor_imu  left_upper_donor_imu  "
            << "right_lower_donor_imu  right_upper_donor_imu  "
            << "left_lower_donor_f  left_upper_donor_f  "
            << "right_lower_donor_f  right_upper_donor_f  "
            << "lower_mu_dot_avg  upper_mu_dot_avg  "
            << "moment_weight  cell_budget  "
            << "low_order_failed_count\n";
        mu_low_failure_audit << std::scientific << std::setprecision(8);

        negative_debt_guard_monitor.open("output/negative_debt_guard.dat");
        negative_debt_guard_monitor
            << "# step time[fs] accepted long_run_diagnostic level "
            << "neg_mass_boundary neg_mass_core neg_mass_tail "
            << "neg_mass_core_fraction "
            << "neg_mass_boundary_fraction neg_mass_tail_fraction "
            << "neg_energy_core_abs neg_energy_core_fraction "
            << "neg_energy_boundary_abs neg_energy_boundary_fraction "
            << "neg_current_core_abs neg_current_core_fraction "
            << "min_f_boundary min_f_core min_f_tail "
            << "neg_cell_boundary neg_cell_core neg_cell_tail "
            << "debt_action limiter_reason "
            << "alpha_core_min alpha_boundary_min alpha_tail_min "
            << "limiter_modified_J_bkg_norm "
            << "limiter_modified_J_bkg_boundary_norm "
            << "limiter_modified_energy_norm "
            << "low_u_mu_neg_mass_fraction "
            << "low_u_mu_neg_energy_fraction "
            << "low_u_mu_neg_current_fraction "
            << "boundary_force_Cu_max boundary_force_Cmu_max "
            << "boundary_force_nsub_max "
            << "boundary_force_remap_cell_count "
            << "boundary_mu_low_L1_before boundary_mu_low_L1_after "
            << "boundary_mu_high_L1_after "
            << "J_bkg_neg_boundary delta_E_neg_boundary "
            << "boundary_force_remap_mass_loss "
            << "boundary_force_remap_energy_loss "
            << "alpha_interface_BQ_min alpha_interface_QC_min "
            << "interface_BQ_flux interface_BQ_high_correction "
            << "interface_QC_flux_into_core "
            << "interface_QC_high_correction_into_core "
            << "boundary_energy_diagnostic_invalid "
            << "trial_failure_downgraded accepted_with_negative_debt "
            << "state_advanced soft_unconverged converged failed "
            << "coupled_iter residual_E residual_J_bkg\n";
        negative_debt_guard_monitor << std::scientific
                                    << std::setprecision(8);

        trial_boundary_debt_monitor.open("output/trial_boundary_debt.dat");
        trial_boundary_debt_monitor
            << "# step time[fs] accepted state_advanced soft_unconverged "
            << "level neg_mass_boundary neg_mass_core neg_mass_tail "
            << "neg_mass_core_fraction "
            << "neg_mass_boundary_fraction neg_mass_tail_fraction "
            << "neg_energy_core_abs neg_energy_core_fraction "
            << "neg_energy_boundary_abs neg_energy_boundary_fraction "
            << "neg_current_core_abs neg_current_core_fraction "
            << "min_f_boundary min_f_core min_f_tail "
            << "neg_cell_boundary neg_cell_core neg_cell_tail "
            << "debt_action limiter_reason "
            << "alpha_core_min alpha_boundary_min alpha_tail_min "
            << "limiter_modified_J_bkg_norm "
            << "limiter_modified_J_bkg_boundary_norm "
            << "limiter_modified_energy_norm "
            << "low_u_mu_neg_mass_fraction "
            << "low_u_mu_neg_energy_fraction "
            << "low_u_mu_neg_current_fraction "
            << "boundary_force_Cu_max boundary_force_Cmu_max "
            << "boundary_force_nsub_max "
            << "boundary_force_remap_cell_count "
            << "boundary_mu_low_L1_before boundary_mu_low_L1_after "
            << "boundary_mu_high_L1_after "
            << "J_bkg_neg_boundary delta_E_neg_boundary "
            << "boundary_force_remap_mass_loss "
            << "boundary_force_remap_energy_loss "
            << "alpha_interface_BQ_min alpha_interface_QC_min "
            << "interface_BQ_flux interface_BQ_high_correction "
            << "interface_QC_flux_into_core "
            << "interface_QC_high_correction_into_core "
            << "boundary_energy_diagnostic_invalid "
            << "trial_failure_downgraded limiter_reason_repeat "
            << "coupled_iter residual_E residual_J_bkg failed\n";
        trial_boundary_debt_monitor << std::scientific
                                    << std::setprecision(8);

        std::ofstream f_neg_monitor;
        f_neg_monitor.open("output/f_negativity_monitor.dat");
        f_neg_monitor
            << "# accepted_snapshot: statistics scan final accepted bkg_e.f; "
            << "neg_ratio_max uses max positive f as scale\n"
            << "# step  time[fs]  accepted  state_advanced  "
            << "soft_unconverged  min_f  neg_ratio_max  "
            << "neg_mass_total[m^-2]  neg_cell_count  "
            << "x_worst  u_worst  mu_worst\n";
        f_neg_monitor << std::scientific << std::setprecision(8);
        f_neg_monitor.close();

        boundary_core_monitor.open("output/boundary_core_diagnostics.dat");
        boundary_core_monitor
            << "# step  time[fs]  accepted  state_advanced  "
            << "soft_unconverged  boundary_width[um]  "
            << "neg_mass_boundary[m^-2]  neg_mass_core[m^-2]  "
            << "neg_mass_core_fraction  "
            << "neg_cell_count_boundary  neg_cell_count_core  "
            << "neg_cell_count_core_fraction  "
            << "min_f_boundary  min_f_core  "
            << "u_limiter_energy_delta_boundary[J/m2]  "
            << "u_limiter_energy_delta_core[J/m2]  "
            << "abs_u_limiter_energy_delta_boundary[J/m2]  "
            << "abs_u_limiter_energy_delta_core[J/m2]  "
            << "u_limiter_energy_delta_core_fraction  "
            << "core_limiter_energy_relative_to_core_work  "
            << "x_limiter_active_fraction_boundary  "
            << "x_limiter_active_fraction_core  "
            << "x_limiter_min_alpha_boundary  "
            << "x_limiter_min_alpha_core  "
            << "max_abs_J_bkg_boundary[A/m2]  "
            << "max_abs_J_bkg_core[A/m2]  "
            << "J_bkg_core_to_boundary_ratio  "
            << "rms_J_bkg_boundary[A/m2]  rms_J_bkg_core[A/m2]  "
            << "mean_abs_J_bkg_boundary[A/m2]  "
            << "mean_abs_J_bkg_core[A/m2]  "
            << "max_abs_Ex_boundary[V/m]  max_abs_Ex_core[V/m]  "
            << "W_bkg_E_boundary[J/m2]  W_bkg_E_core[J/m2]  "
            << "anomaly_inward_E[um]  anomaly_inward_J[um]  "
            << "anomaly_inward_n[um]  anomaly_inward_Mneg[um]  "
            << "conclusion\n";
        boundary_core_monitor << std::scientific << std::setprecision(8);

        x_low_monitor.open("output/x_low_diagnostics.dat");
        x_low_monitor
            << "# step  time[fs]  x_low_input_min_f  x_low_max_cfl  "
            << "x_low_output_min_f  x_low_failed_count\n";
        x_low_monitor << std::scientific << std::setprecision(8);

        x_final_monitor.open("output/x_final_diagnostics.dat");
        x_final_monitor
            << "# step  time[fs]  x_final_min_f  "
            << "x_final_failed_count  x_final_failed_max_debt  "
            << "x_final_worst_ix  x_final_worst_iv  "
            << "x_final_worst_imu  x_final_failure_region  "
            << "x_final_core_failed_count  "
            << "x_final_boundary_failed_count\n";
        x_final_monitor << std::scientific << std::setprecision(8);

        mu_final_monitor.open("output/mu_final_diagnostics.dat");
        mu_final_monitor
            << "# step  time[fs]  mu_final_min_f  "
            << "mu_final_failed_count  mu_final_failed_max_debt  "
            << "mu_final_worst_ix  mu_final_worst_iv  "
            << "mu_final_worst_imu  mu_final_failure_region  "
            << "mu_final_core_failed_count  "
            << "mu_final_boundary_failed_count  "
            << "mu_alpha_active_fraction  mu_alpha_min  "
            << "mu_alpha_core_fraction  mu_alpha_boundary_fraction  "
            << "mu_alpha_core_min  mu_alpha_boundary_min  "
            << "mu_limiter_energy_delta[J/m2]  "
            << "mu_limiter_mass_delta[m^-2]\n";
        mu_final_monitor << std::scientific << std::setprecision(8);

    }
    bool last_accepted_stage_valid = false;
    double last_accepted_stage_min_f[3] = {0.0, 0.0, 0.0};
    double last_accepted_stage_neg_mass[3] = {0.0, 0.0, 0.0};
    long long last_accepted_stage_neg_cell_count[3] = {0, 0, 0};
    for (int step = 1; step <= nsteps; ++step) {
        double time = step * dt;
        int nsub_v1 = 0;
        int nsub_mu1 = 0;
        int nsub_v2 = 0;
        int nsub_mu2 = 0;
        double loss_v1 = 0.0;
        double loss_v1_low = 0.0;
        double loss_v1_high = 0.0;
        double loss_mu1 = 0.0;
        double loss_v2 = 0.0;
        double loss_v2_low = 0.0;
        double loss_v2_high = 0.0;
        double loss_mu2 = 0.0;
        double net_nb_change_step = 0.0;
        double collision_energy_step = 0.0;
        const bool collect_step_diagnostics =
            should_write_step_diagnostics(config, step);
        double dke_bkg_step = 0.0;
        double dke_beam_push = 0.0;
        double dE_field_step = 0.0;
        double W_bkg_E = 0.0;
        double W_beam_E = 0.0;
        double v_mass_error_step = 0.0;
        double mu_mass_error_step = 0.0;
        double v_momentum_delta_step = 0.0;
        double mu_momentum_delta_step = 0.0;
        double v_energy_delta_step = 0.0;
        double mu_energy_delta_step = 0.0;
        double E_src_in_step = 0.0;
        double E_src_out_step = 0.0;
        double E_balance_step = 0.0;
        double max_loss_u_high_step = 0.0;
        double x_at_max_loss_u_high_step = 0.0;
        double f_u_max_x_step = 0.0;
        double integral_f_u_gt_8_x_step = 0.0;
        double bkg_number_step_start = 0.0;
        double bkg_ke_step_start = 0.0;
        double global_bkg_ke_step_start = 0.0;
        double beam_ke_step_start = 0.0;
        double field_energy_step_start = 0.0;
        double x_limiter_active_fraction_step = 0.0;
        double x_limiter_min_alpha_step = 1.0;
        double x_limiter_active_fraction_core_step = 0.0;
        double x_limiter_active_fraction_boundary_step = 0.0;
        double x_limiter_min_alpha_core_step = 1.0;
        double x_limiter_min_alpha_boundary_step = 1.0;
        double x_negative_mass_before_repair_step = 0.0;
        double x_mass_added_by_positivity_repair_step = 0.0;
        double floor_repair_mass_step = 0.0;
        double floor_repair_energy_step = 0.0;
        double floor_repair_core_fraction_step = 0.0;
        double positivity_energy_defect_step = 0.0;
        double positivity_mass_defect_step = 0.0;
        double u_limiter_mass_delta_step = 0.0;
        double u_limiter_momentum_delta_step = 0.0;
        double u_limiter_energy_delta_step = 0.0;
        double x_limiter_mass_delta_step = 0.0;
        double x_limiter_energy_delta_step = 0.0;
        double x_low_input_min_f_step = 0.0;
        double x_low_max_cfl_step = 0.0;
        double x_low_output_min_f_step = 0.0;
        double x_low_failed_count_step = 0.0;
        double x_low_input_neg_mass_step = 0.0;
        double x_low_input_rel_neg_step = 0.0;
        double x_low_output_rel_neg_step = 0.0;
        double x_low_input_core_failed_count_step = 0.0;
        double x_low_input_debt_accepted_step = 0.0;
        int x_low_failure_kind_step = 0;
        double x_final_min_f_step = 0.0;
        double x_final_failed_count_step = 0.0;
        double x_final_failed_max_debt_step = 0.0;
        int x_final_worst_ix_step = -1;
        int x_final_worst_iv_step = -1;
        int x_final_worst_imu_step = -1;
        int x_final_failure_region_step = -1;
        double x_final_core_failed_count_step = 0.0;
        double x_final_boundary_failed_count_step = 0.0;
        double mu_final_min_f_step = 0.0;
        double mu_final_failed_count_step = 0.0;
        double mu_final_failed_max_debt_step = 0.0;
        int mu_final_worst_ix_step = -1;
        int mu_final_worst_iv_step = -1;
        int mu_final_worst_imu_step = -1;
        int mu_final_failure_region_step = -1;
        double mu_final_core_failed_count_step = 0.0;
        double mu_final_boundary_failed_count_step = 0.0;
        double mu_low_u_alpha_min_step = 1.0;
        double mu_low_u_limiter_active_fraction_step = 0.0;
        double mu_low_u_energy_delta_step = 0.0;
        double mu_low_u_alpha_min_boundary_step = 1.0;
        double mu_low_u_alpha_min_core_step = 1.0;
        double mu_low_u_limiter_active_fraction_boundary_step = 0.0;
        double mu_low_u_limiter_active_fraction_core_step = 0.0;
        double mu_low_u_energy_delta_boundary_step = 0.0;
        double mu_low_u_energy_delta_core_step = 0.0;
        double mu_low_u_u_eff0_step = 0.0;
        double mu_low_u_moment_weight0_step = 0.0;
        double mu_low_u_mu_flux_scale0_step = 0.0;
        double mu_low_u_half_dt_inv_shell0_step = 0.0;
        double mu_low_u_dimless_scale0_step = 0.0;
        double mu_low_u_endpoint_flux_max_step = 0.0;
        double remap_active_fraction_step = 0.0;
        long long remap_cell_count_step = 0;
        double low_u_subcycle_active_fraction_step = 0.0;
        double low_u_average_subcycles_step = 1.0;
        int low_u_max_subcycles_step = 1;
        double u_force_alpha_min_step = 1.0;
        double u_force_alpha_active_frac_step = 0.0;
        double local_bkg_energy_residual_step = 0.0;
        double bkg_energy_residual_step = 0.0;
        double bkg_energy_relative_residual_step = 0.0;
        int coupled_iter_step = 0;
        double coupled_residual_E_step = 0.0;
        double coupled_residual_J_bkg_step = 0.0;
        double coupled_residual_J_beam_step = 0.0;
        double coupled_residual_bkg_mass_step = 0.0;
        double coupled_residual_beam_continuity_step = 0.0;
        BackgroundCurrentDiagnostics local_bkg_current_diag_step;
        reset_background_current_diagnostics(local_bkg_current_diag_step);
        BackgroundCurrentDiagnostics global_bkg_current_diag_step;
        reset_background_current_diagnostics(global_bkg_current_diag_step);
        const Species bkg_step_start = bkg_e;
        const BeamPIC beam_step_start = beam;
        const EMFields fields_step_start = fields;
        bkg_e.total_particle_number_and_energy(bkg_number_step_start,
                                               bkg_ke_step_start);
        double global_bkg_start_values[2] = {
            bkg_number_step_start,
            bkg_ke_step_start
        };
        MPI_Allreduce(MPI_IN_PLACE, global_bkg_start_values, 2,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_bkg_ke_step_start = global_bkg_start_values[1];
        if (collect_step_diagnostics) {
            beam_ke_step_start = beam.total_kinetic_energy();
            field_energy_step_start = fields.total_energy();
        }

        trace_progress(config, mpi_rank, step,
                       "before coupled midpoint FV solve");
        midpoint_solver.set_step_diagnostics_enabled(collect_step_diagnostics);
        VlasovAmpereMidpointSolver::Result midpoint_result =
            midpoint_solver.advance_background_and_fields(
                bkg_step_start, beam_step_start, fields_step_start, sgrid,
                dt, time, mpi_rank, mpi_size);
        const bool negative_debt_attention =
            midpoint_result.negative_debt_level >=
            VlasovAmpereMidpointSolver::NEG_DEBT_WARN ||
            midpoint_result.accepted_with_negative_debt != 0 ||
            midpoint_result.trial_failure_downgraded != 0;
        const bool abnormal_midpoint =
            !midpoint_result.converged || midpoint_result.failed ||
            negative_debt_attention;
        trace_progress(config, mpi_rank, step,
                       "after coupled midpoint FV solve");

        auto write_low_order_failure_audit =
            [&](std::ofstream& out,
                const VlasovAmpereMidpointSolver::LowOrderFailureAudit&
                    audit) {
                if (!audit.valid) return;
                const char* region =
                    (audit.region != 0) ? "boundary" : "core";
                out << step << "  "
                    << time / Const::femto << "  "
                    << audit.rank << "  "
                    << audit.ix << "  "
                    << audit.iv << "  "
                    << audit.imu << "  "
                    << region << "  "
                    << audit.severity << "  "
                    << audit.f_input << "  "
                    << audit.f_after_x << "  "
                    << audit.dx_div << "  "
                    << audit.dmu_div_used << "  "
                    << audit.du_div_low << "  "
                    << audit.f_low << "  "
                    << audit.left_lower_flux << "  "
                    << audit.left_upper_flux << "  "
                    << audit.right_lower_flux << "  "
                    << audit.right_upper_flux << "  "
                    << audit.left_lower_scale << "  "
                    << audit.left_upper_scale << "  "
                    << audit.right_lower_scale << "  "
                    << audit.right_upper_scale << "  "
                    << audit.left_lower_donor_index << "  "
                    << audit.left_upper_donor_index << "  "
                    << audit.right_lower_donor_index << "  "
                    << audit.right_upper_donor_index << "  "
                    << audit.left_lower_donor_f << "  "
                    << audit.left_upper_donor_f << "  "
                    << audit.right_lower_donor_f << "  "
                    << audit.right_upper_donor_f << "  "
                    << audit.lower_characteristic << "  "
                    << audit.upper_characteristic << "  "
                    << audit.moment_weight << "  "
                    << audit.cell_budget << "  "
                    << audit.low_order_failed_count << "\n";
                out.flush();
            };

        if (mpi_rank == 0) {
            trial_boundary_debt_monitor
                << step << " "
                << time / Const::femto << " "
                << 0 << " "
                << midpoint_result.state_advanced << " "
                << midpoint_result.soft_unconverged << " "
                << midpoint_result.negative_debt_level << " "
                << midpoint_result.neg_mass_boundary << " "
                << midpoint_result.neg_mass_core << " "
                << midpoint_result.neg_mass_tail << " "
                << midpoint_result.neg_mass_core_fraction << " "
                << midpoint_result.neg_mass_boundary_fraction << " "
                << midpoint_result.neg_mass_tail_fraction << " "
                << midpoint_result.neg_energy_core_abs << " "
                << midpoint_result.neg_energy_core_fraction << " "
                << midpoint_result.neg_energy_boundary_abs << " "
                << midpoint_result.neg_energy_boundary_fraction << " "
                << midpoint_result.neg_current_core_abs << " "
                << midpoint_result.neg_current_core_fraction << " "
                << midpoint_result.neg_debt_min_f_boundary << " "
                << midpoint_result.neg_debt_min_f_core << " "
                << midpoint_result.neg_debt_min_f_tail << " "
                << midpoint_result.neg_cell_boundary << " "
                << midpoint_result.neg_cell_core << " "
                << midpoint_result.neg_cell_tail << " "
                << midpoint_result.debt_action << " "
                << midpoint_result.limiter_reason << " "
                << midpoint_result.alpha_core_min << " "
                << midpoint_result.alpha_boundary_min << " "
                << midpoint_result.alpha_tail_min << " "
                << midpoint_result.limiter_modified_J_bkg_norm << " "
                << midpoint_result.limiter_modified_J_bkg_boundary_norm << " "
                << midpoint_result.limiter_modified_energy_norm << " "
                << midpoint_result.low_u_mu_neg_mass_fraction << " "
                << midpoint_result.low_u_mu_neg_energy_fraction << " "
                << midpoint_result.low_u_mu_neg_current_fraction << " "
                << midpoint_result.boundary_force_Cu_max << " "
                << midpoint_result.boundary_force_Cmu_max << " "
                << midpoint_result.boundary_force_nsub_max << " "
                << midpoint_result.boundary_force_remap_cell_count << " "
                << midpoint_result.boundary_mu_low_L1_before << " "
                << midpoint_result.boundary_mu_low_L1_after << " "
                << midpoint_result.boundary_mu_high_L1_after << " "
                << midpoint_result.J_bkg_neg_boundary << " "
                << midpoint_result.delta_E_neg_boundary << " "
                << midpoint_result.boundary_force_remap_mass_loss << " "
                << midpoint_result.boundary_force_remap_energy_loss << " "
                << midpoint_result.alpha_interface_BQ_min << " "
                << midpoint_result.alpha_interface_QC_min << " "
                << midpoint_result.interface_BQ_flux << " "
                << midpoint_result.interface_BQ_high_correction << " "
                << midpoint_result.interface_QC_flux_into_core << " "
                << midpoint_result.interface_QC_high_correction_into_core << " "
                << midpoint_result.boundary_energy_diagnostic_invalid << " "
                << midpoint_result.trial_failure_downgraded << " "
                << midpoint_result.limiter_reason << " "
                << midpoint_result.nonlinear_iterations << " "
                << midpoint_result.residual_E << " "
                << midpoint_result.residual_J_bkg << " "
                << midpoint_result.failed << " "
                << "\n";
            if (abnormal_midpoint ||
                step + 1 == nsteps ||
                (config.enable_step_diagnostics &&
                 config.step_diagnostics_interval > 0 &&
                 step % config.step_diagnostics_interval == 0)) {
                trial_boundary_debt_monitor.flush();
            }
        }

        const bool midpoint_soft_accepted =
            midpoint_result.soft_unconverged &&
            midpoint_result.state_advanced != 0 &&
            !midpoint_result.failed;
        if ((!midpoint_result.converged && !midpoint_soft_accepted) ||
            midpoint_result.failed) {
            if (mpi_rank == 0) {
                x_low_monitor << step << "  "
                              << time / Const::femto << "  "
                              << midpoint_result.x_low_input_min_f << "  "
                              << midpoint_result.x_low_max_cfl << "  "
                              << midpoint_result.x_low_output_min_f << "  "
                              << midpoint_result.x_low_failed_count << "\n";
                x_low_monitor.flush();
                x_final_monitor << step << "  "
                                << time / Const::femto << "  "
                                << midpoint_result.x_final_min_f << "  "
                                << midpoint_result.x_final_failed_count << "  "
                                << midpoint_result.x_final_failed_max_debt << "  "
                                << midpoint_result.x_final_worst_ix << "  "
                                << midpoint_result.x_final_worst_iv << "  "
                                << midpoint_result.x_final_worst_imu << "  "
                                << midpoint_result.x_final_failure_region << "  "
                                << midpoint_result.x_final_core_failed_count << "  "
                                << midpoint_result.x_final_boundary_failed_count << "\n";
                x_final_monitor.flush();
                mu_final_monitor << step << "  "
                                 << time / Const::femto << "  "
                                 << midpoint_result.flux_pos[2].min_f_final << "  "
                                 << midpoint_result.mu_final_negative_hard_count << "  "
                                 << midpoint_result.mu_final_failed_max_debt << "  "
                                 << midpoint_result.mu_final_worst_ix << "  "
                                 << midpoint_result.mu_final_worst_iv << "  "
                                 << midpoint_result.mu_final_worst_imu << "  "
                                 << midpoint_result.mu_final_failure_region << "  "
                                 << midpoint_result.mu_final_core_failed_count << "  "
                                 << midpoint_result.mu_final_boundary_failed_count << "  "
                                 << midpoint_result.flux_pos[2].alpha_active_fraction << "  "
                                 << midpoint_result.flux_pos[2].alpha_min << "  "
                                 << midpoint_result.flux_pos[2].alpha_core_fraction << "  "
                                 << midpoint_result.flux_pos[2].alpha_boundary_fraction << "  "
                                 << midpoint_result.mu_low_u_alpha_min_core << "  "
                                 << midpoint_result.mu_low_u_alpha_min_boundary << "  "
                                 << midpoint_result.flux_defect[2].energy_defect << "  "
                                 << midpoint_result.flux_defect[2].mass_defect << "\n";
                mu_final_monitor.flush();
                write_low_order_failure_audit(
                    u_low_failure_audit,
                    midpoint_result.u_low_failure_audit);
                write_low_order_failure_audit(
                    mu_low_failure_audit,
                    midpoint_result.mu_low_failure_audit);
                std::fprintf(stderr,
                             "WARNING: coupled midpoint FV solve failed at "
                             "step %d, t = %.6e s; residual %.6e, "
                             "field %.6e, f %.6e, J_bkg %.6e, J_beam %.6e, "
                             "bkg_mass %.6e, beam_continuity %.6e after "
                             "%d iterations, substeps %d; limiter active "
                             "%.6e, min_alpha %.6e, core_active %.6e, "
                             "boundary_active %.6e, core_min_alpha %.6e, "
                             "boundary_min_alpha %.6e, "
                             "x_low_input_min_f %.6e, "
                             "x_low_max_cfl %.6e, "
                             "x_low_output_min_f %.6e, "
                             "x_low_failed_count %.0f, "
                             "x_low_input_neg_mass %.6e, "
                             "x_low_input_rel_neg %.6e, "
                             "x_low_output_rel_neg %.6e, "
                             "x_low_input_core_failed_count %.0f, "
                             "x_low_input_debt_accepted %.0f, "
                             "x_low_failure_kind %d, "
                             "x_final(min %.6e, failed %.0f, "
                             "max_debt %.6e, worst ix %d iv %d imu %d, "
                             "region %d, core_failed %.0f, "
                             "boundary_failed %.0f), "
                             "finite_failure_mask %d, "
                             "finite_counts(u_low %lld, u_final %lld, "
                             "mu_low %lld, mu_final %lld), "
                             "finite_stage(valid %d, kind %d, rank %d, "
                             "ix %d, iv %d, imu %d, severity %.6e, "
                             "f_base %.6e, f_low %.6e, f_high %.6e, "
                             "f_final %.6e, dx_div %.6e, du_div %.6e, "
                             "dmu_div %.6e). "
                             "previous accepted stage valid %d; "
                             "prev after_x(min %.6e, neg_mass %.6e, "
                             "neg_count %lld), "
                             "after_u(min %.6e, neg_mass %.6e, "
                             "neg_count %lld), "
                             "after_mu(min %.6e, neg_mass %.6e, "
                             "neg_count %lld). "
                             "Soft negative-debt failures are accepted only "
                             "when state_advanced=1; frozen states remain "
                             "hard failures.\n",
                             step, time,
                             midpoint_result.nonlinear_residual,
                             midpoint_result.residual_E,
                             midpoint_result.residual_f,
                             midpoint_result.residual_J_bkg,
                             midpoint_result.residual_J_beam,
                             midpoint_result.continuity_residual_bkg,
                             midpoint_result.beam_continuity_residual,
                             midpoint_result.nonlinear_iterations,
                             midpoint_result.substeps_used,
                             midpoint_result.limiter_active_fraction,
                             midpoint_result.limiter_min_alpha,
                             midpoint_result.limiter_active_fraction_core,
                             midpoint_result.limiter_active_fraction_boundary,
                             midpoint_result.limiter_min_alpha_core,
                             midpoint_result.limiter_min_alpha_boundary,
                             midpoint_result.x_low_input_min_f,
                             midpoint_result.x_low_max_cfl,
                             midpoint_result.x_low_output_min_f,
                             midpoint_result.x_low_failed_count,
                             midpoint_result.x_low_input_neg_mass,
                             midpoint_result.x_low_input_rel_neg,
                             midpoint_result.x_low_output_rel_neg,
                             midpoint_result.x_low_input_core_failed_count,
                             midpoint_result.x_low_input_debt_accepted,
                             midpoint_result.x_low_failure_kind,
                             midpoint_result.x_final_min_f,
                             midpoint_result.x_final_failed_count,
                             midpoint_result.x_final_failed_max_debt,
                             midpoint_result.x_final_worst_ix,
                             midpoint_result.x_final_worst_iv,
                             midpoint_result.x_final_worst_imu,
                             midpoint_result.x_final_failure_region,
                             midpoint_result.x_final_core_failed_count,
                             midpoint_result.x_final_boundary_failed_count,
                             midpoint_result.finite_failure_mask,
                             midpoint_result.u_low_order_failed_count,
                             midpoint_result.u_final_negative_hard_count,
                             midpoint_result.mu_low_order_failed_count,
                             midpoint_result.mu_final_negative_hard_count,
                             midpoint_result.finite_stage_failure_valid,
                             midpoint_result.finite_stage_failure_kind,
                             midpoint_result.finite_stage_failure_rank,
                             midpoint_result.finite_stage_failure_ix,
                             midpoint_result.finite_stage_failure_iv,
                             midpoint_result.finite_stage_failure_imu,
                             midpoint_result.finite_stage_failure_severity,
                             midpoint_result.finite_stage_failure_f_base,
                             midpoint_result.finite_stage_failure_f_low,
                             midpoint_result.finite_stage_failure_f_high,
                             midpoint_result.finite_stage_failure_f_final,
                             midpoint_result.finite_stage_failure_dx_div,
                             midpoint_result.finite_stage_failure_du_div,
                             midpoint_result.finite_stage_failure_dmu_div,
                             last_accepted_stage_valid ? 1 : 0,
                             last_accepted_stage_min_f[0],
                             last_accepted_stage_neg_mass[0],
                             last_accepted_stage_neg_cell_count[0],
                             last_accepted_stage_min_f[1],
                             last_accepted_stage_neg_mass[1],
                             last_accepted_stage_neg_cell_count[1],
                             last_accepted_stage_min_f[2],
                             last_accepted_stage_neg_mass[2],
                             last_accepted_stage_neg_cell_count[2]);
            }
            const bool hard_negative_debt =
                midpoint_result.negative_debt_level >=
                VlasovAmpereMidpointSolver::NEG_DEBT_ABORT;
            const bool frozen_output = midpoint_result.state_advanced == 0;
            const bool soft_negative_debt_failure =
                midpoint_result.trial_failure_downgraded != 0 &&
                !hard_negative_debt && !frozen_output;
            if (soft_negative_debt_failure) {
                if (mpi_rank == 0) {
                    std::fprintf(
                        stderr,
                        "WARNING: accepting downgraded midpoint failure at "
                        "step %d with negative_debt_level=%d and "
                        "state_advanced=%d.\n",
                        step, midpoint_result.negative_debt_level,
                        midpoint_result.state_advanced);
                }
            } else {
                if (mpi_rank == 0) {
                    std::fprintf(
                        stderr,
                        "ERROR: hard midpoint failure at step %d: "
                        "negative_debt_level=%d state_advanced=%d "
                        "trial_failure_downgraded=%d.\n",
                        step, midpoint_result.negative_debt_level,
                        midpoint_result.state_advanced,
                        midpoint_result.trial_failure_downgraded);
                }
                MPI_Abort(MPI_COMM_WORLD, 8);
                return 8;
            }
        }

        if (midpoint_result.state_advanced == 0) {
            if (mpi_rank == 0) {
                std::fprintf(
                    stderr,
                    "ERROR: midpoint solver returned without advancing the "
                    "state at step %d; aborting to prevent frozen output.\n",
                    step);
            }
            MPI_Abort(MPI_COMM_WORLD, 9);
            return 9;
        }

        bkg_e = midpoint_result.species_np1;
        beam = midpoint_result.beam_np1;
        fields = midpoint_result.fields_np1;
        if (mpi_rank == 0) {
            const int long_run_diagnostic =
                midpoint_result.soft_unconverged ? 1 : 0;
            negative_debt_guard_monitor
                << step << " "
                << time / Const::femto << " "
                << 1 << " "
                << long_run_diagnostic << " "
                << midpoint_result.negative_debt_level << " "
                << midpoint_result.neg_mass_boundary << " "
                << midpoint_result.neg_mass_core << " "
                << midpoint_result.neg_mass_tail << " "
                << midpoint_result.neg_mass_core_fraction << " "
                << midpoint_result.neg_mass_boundary_fraction << " "
                << midpoint_result.neg_mass_tail_fraction << " "
                << midpoint_result.neg_energy_core_abs << " "
                << midpoint_result.neg_energy_core_fraction << " "
                << midpoint_result.neg_energy_boundary_abs << " "
                << midpoint_result.neg_energy_boundary_fraction << " "
                << midpoint_result.neg_current_core_abs << " "
                << midpoint_result.neg_current_core_fraction << " "
                << midpoint_result.neg_debt_min_f_boundary << " "
                << midpoint_result.neg_debt_min_f_core << " "
                << midpoint_result.neg_debt_min_f_tail << " "
                << midpoint_result.neg_cell_boundary << " "
                << midpoint_result.neg_cell_core << " "
                << midpoint_result.neg_cell_tail << " "
                << midpoint_result.debt_action << " "
                << midpoint_result.limiter_reason << " "
                << midpoint_result.alpha_core_min << " "
                << midpoint_result.alpha_boundary_min << " "
                << midpoint_result.alpha_tail_min << " "
                << midpoint_result.limiter_modified_J_bkg_norm << " "
                << midpoint_result.limiter_modified_J_bkg_boundary_norm << " "
                << midpoint_result.limiter_modified_energy_norm << " "
                << midpoint_result.low_u_mu_neg_mass_fraction << " "
                << midpoint_result.low_u_mu_neg_energy_fraction << " "
                << midpoint_result.low_u_mu_neg_current_fraction << " "
                << midpoint_result.boundary_force_Cu_max << " "
                << midpoint_result.boundary_force_Cmu_max << " "
                << midpoint_result.boundary_force_nsub_max << " "
                << midpoint_result.boundary_force_remap_cell_count << " "
                << midpoint_result.boundary_mu_low_L1_before << " "
                << midpoint_result.boundary_mu_low_L1_after << " "
                << midpoint_result.boundary_mu_high_L1_after << " "
                << midpoint_result.J_bkg_neg_boundary << " "
                << midpoint_result.delta_E_neg_boundary << " "
                << midpoint_result.boundary_force_remap_mass_loss << " "
                << midpoint_result.boundary_force_remap_energy_loss << " "
                << midpoint_result.alpha_interface_BQ_min << " "
                << midpoint_result.alpha_interface_QC_min << " "
                << midpoint_result.interface_BQ_flux << " "
                << midpoint_result.interface_BQ_high_correction << " "
                << midpoint_result.interface_QC_flux_into_core << " "
                << midpoint_result.interface_QC_high_correction_into_core << " "
                << midpoint_result.boundary_energy_diagnostic_invalid << " "
                << midpoint_result.trial_failure_downgraded << " "
                << midpoint_result.accepted_with_negative_debt << " "
                << midpoint_result.state_advanced << " "
                << midpoint_result.soft_unconverged << " "
                << midpoint_result.converged << " "
                << midpoint_result.failed << " "
                << midpoint_result.nonlinear_iterations << " "
                << midpoint_result.residual_E << " "
                << midpoint_result.residual_J_bkg << "\n";
            if (abnormal_midpoint ||
                step + 1 == nsteps ||
                (config.enable_step_diagnostics &&
                 config.step_diagnostics_interval > 0 &&
                 step % config.step_diagnostics_interval == 0)) {
                negative_debt_guard_monitor.flush();
            }
        }
        latest_bkg_energy_current_face =
            midpoint_result.j_bkg_energy_debug_face;
        moments_current = true;
        latest_bkg_ampere_current_face = midpoint_result.j_bkg_face_mid;
        local_bkg_current_diag_step.residual_if_charge =
            midpoint_result.current_diag.residual_if_charge;
        local_bkg_current_diag_step.residual_if_ampere =
            midpoint_result.current_diag.residual_if_ampere;
        local_bkg_current_diag_step.e_dot_j_charge =
            midpoint_result.current_diag.e_dot_j_charge;
        local_bkg_current_diag_step.e_dot_j_energy =
            midpoint_result.current_diag.e_dot_j_energy;
        local_bkg_current_diag_step.e_dot_j_ampere =
            midpoint_result.current_diag.e_dot_j_ampere;
        local_bkg_current_diag_step.max_abs_j_charge =
            midpoint_result.current_diag.max_abs_j_charge;
        local_bkg_current_diag_step.max_abs_j_energy =
            midpoint_result.current_diag.max_abs_j_energy;
        local_bkg_current_diag_step.max_abs_j_ampere =
            midpoint_result.current_diag.max_abs_j_ampere;
        local_bkg_current_diag_step.max_abs_j_charge_minus_ampere =
            midpoint_result.current_diag.max_abs_j_charge_minus_ampere;
        local_bkg_current_diag_step.max_abs_j_energy_minus_ampere =
            midpoint_result.current_diag.max_abs_j_energy_minus_ampere;
        W_bkg_E = midpoint_result.field_work_bkg;
        W_beam_E = midpoint_result.field_work_beam;
        dke_beam_push = midpoint_result.delta_ke_beam;
        x_limiter_active_fraction_step =
            midpoint_result.limiter_active_fraction;
        x_limiter_min_alpha_step = midpoint_result.limiter_min_alpha;
        x_limiter_active_fraction_core_step =
            midpoint_result.limiter_active_fraction_core;
        x_limiter_active_fraction_boundary_step =
            midpoint_result.limiter_active_fraction_boundary;
        x_limiter_min_alpha_core_step =
            midpoint_result.limiter_min_alpha_core;
        x_limiter_min_alpha_boundary_step =
            midpoint_result.limiter_min_alpha_boundary;
        x_negative_mass_before_repair_step =
            midpoint_result.x_negative_mass_before_repair;
        x_mass_added_by_positivity_repair_step =
            midpoint_result.x_mass_added_by_positivity_repair;
        floor_repair_mass_step = midpoint_result.floor_repair_mass;
        floor_repair_energy_step = midpoint_result.floor_repair_energy;
        floor_repair_core_fraction_step =
            midpoint_result.floor_repair_core_fraction;
        positivity_energy_defect_step =
            midpoint_result.positivity_energy_defect;
        positivity_mass_defect_step =
            midpoint_result.positivity_mass_defect;
        u_limiter_mass_delta_step =
            midpoint_result.limiter_mass_defect;
        u_limiter_momentum_delta_step =
            midpoint_result.limiter_momentum_defect;
        u_limiter_energy_delta_step =
            midpoint_result.limiter_energy_defect;
        x_limiter_mass_delta_step =
            midpoint_result.x_limiter_mass_defect;
        x_limiter_energy_delta_step =
            midpoint_result.x_limiter_energy_defect;
        x_low_input_min_f_step = midpoint_result.x_low_input_min_f;
        x_low_max_cfl_step = midpoint_result.x_low_max_cfl;
        x_low_output_min_f_step = midpoint_result.x_low_output_min_f;
        x_low_failed_count_step = midpoint_result.x_low_failed_count;
        x_low_input_neg_mass_step = midpoint_result.x_low_input_neg_mass;
        x_low_input_rel_neg_step = midpoint_result.x_low_input_rel_neg;
        x_low_output_rel_neg_step = midpoint_result.x_low_output_rel_neg;
        x_low_input_core_failed_count_step =
            midpoint_result.x_low_input_core_failed_count;
        x_low_input_debt_accepted_step =
            midpoint_result.x_low_input_debt_accepted;
        x_low_failure_kind_step = midpoint_result.x_low_failure_kind;
        (void)x_low_failure_kind_step;
        x_final_min_f_step = midpoint_result.x_final_min_f;
        x_final_failed_count_step = midpoint_result.x_final_failed_count;
        x_final_failed_max_debt_step =
            midpoint_result.x_final_failed_max_debt;
        x_final_worst_ix_step = midpoint_result.x_final_worst_ix;
        x_final_worst_iv_step = midpoint_result.x_final_worst_iv;
        x_final_worst_imu_step = midpoint_result.x_final_worst_imu;
        x_final_failure_region_step =
            midpoint_result.x_final_failure_region;
        x_final_core_failed_count_step =
            midpoint_result.x_final_core_failed_count;
        x_final_boundary_failed_count_step =
            midpoint_result.x_final_boundary_failed_count;
        mu_final_min_f_step = midpoint_result.flux_pos[2].min_f_final;
        mu_final_failed_count_step =
            static_cast<double>(midpoint_result.mu_final_negative_hard_count);
        mu_final_failed_max_debt_step =
            midpoint_result.mu_final_failed_max_debt;
        mu_final_worst_ix_step = midpoint_result.mu_final_worst_ix;
        mu_final_worst_iv_step = midpoint_result.mu_final_worst_iv;
        mu_final_worst_imu_step = midpoint_result.mu_final_worst_imu;
        mu_final_failure_region_step =
            midpoint_result.mu_final_failure_region;
        mu_final_core_failed_count_step =
            midpoint_result.mu_final_core_failed_count;
        mu_final_boundary_failed_count_step =
            midpoint_result.mu_final_boundary_failed_count;
        if (mpi_rank == 0 && x_low_input_debt_accepted_step > 0.0) {
            std::fprintf(
                stderr,
                "WARNING: accepted x_low input negative debt at step %d, "
                "t = %.6e s; failed_cell_count %.0f, input_min_f %.6e, "
                "input_neg_mass %.6e, input_rel_neg %.6e, "
                "output_min_f %.6e, output_rel_neg %.6e, "
                "core_failed_count %.0f. Previous accepted stage valid %d: "
                "after_x(min %.6e, neg_mass %.6e, neg_count %lld), "
                "after_u(min %.6e, neg_mass %.6e, neg_count %lld), "
                "after_mu(min %.6e, neg_mass %.6e, neg_count %lld).\n",
                step, time, x_low_failed_count_step,
                x_low_input_min_f_step,
                x_low_input_neg_mass_step,
                x_low_input_rel_neg_step,
                x_low_output_min_f_step,
                x_low_output_rel_neg_step,
                x_low_input_core_failed_count_step,
                last_accepted_stage_valid ? 1 : 0,
                last_accepted_stage_min_f[0],
                last_accepted_stage_neg_mass[0],
                last_accepted_stage_neg_cell_count[0],
                last_accepted_stage_min_f[1],
                last_accepted_stage_neg_mass[1],
                last_accepted_stage_neg_cell_count[1],
                last_accepted_stage_min_f[2],
                last_accepted_stage_neg_mass[2],
                last_accepted_stage_neg_cell_count[2]);
        }
        mu_low_u_alpha_min_step =
            midpoint_result.mu_low_u_alpha_min;
        mu_low_u_limiter_active_fraction_step =
            midpoint_result.mu_low_u_limiter_active_fraction;
        mu_low_u_energy_delta_step =
            midpoint_result.mu_low_u_energy_delta;
        mu_low_u_alpha_min_boundary_step =
            midpoint_result.mu_low_u_alpha_min_boundary;
        mu_low_u_alpha_min_core_step =
            midpoint_result.mu_low_u_alpha_min_core;
        mu_low_u_limiter_active_fraction_boundary_step =
            midpoint_result.mu_low_u_limiter_active_fraction_boundary;
        mu_low_u_limiter_active_fraction_core_step =
            midpoint_result.mu_low_u_limiter_active_fraction_core;
        mu_low_u_energy_delta_boundary_step =
            midpoint_result.mu_low_u_energy_delta_boundary;
        mu_low_u_energy_delta_core_step =
            midpoint_result.mu_low_u_energy_delta_core;
        mu_low_u_u_eff0_step = midpoint_result.mu_low_u_u_eff0;
        mu_low_u_moment_weight0_step =
            midpoint_result.mu_low_u_moment_weight0;
        mu_low_u_mu_flux_scale0_step =
            midpoint_result.mu_low_u_mu_flux_scale0;
        mu_low_u_half_dt_inv_shell0_step =
            midpoint_result.mu_low_u_half_dt_inv_shell0;
        mu_low_u_dimless_scale0_step =
            midpoint_result.mu_low_u_dimless_scale0;
        mu_low_u_endpoint_flux_max_step =
            midpoint_result.mu_low_u_endpoint_flux_max;
        remap_active_fraction_step =
            midpoint_result.remap_active_fraction;
        remap_cell_count_step = midpoint_result.remap_cell_count;
        low_u_subcycle_active_fraction_step =
            midpoint_result.low_u_subcycle_active_fraction;
        low_u_average_subcycles_step =
            midpoint_result.low_u_average_subcycles;
        low_u_max_subcycles_step = midpoint_result.low_u_max_subcycles;
        u_force_alpha_min_step =
            midpoint_result.u_force_alpha_min;
        u_force_alpha_active_frac_step =
            midpoint_result.u_force_alpha_active_frac;
        coupled_iter_step = midpoint_result.nonlinear_iterations;
        coupled_residual_E_step = midpoint_result.residual_E;
        coupled_residual_J_bkg_step = midpoint_result.residual_J_bkg;
        coupled_residual_J_beam_step = midpoint_result.residual_J_beam;
        coupled_residual_bkg_mass_step =
            midpoint_result.continuity_residual_bkg;
        coupled_residual_beam_continuity_step =
            midpoint_result.beam_continuity_residual;

        const bool write_flux_audit =
            collect_step_diagnostics || abnormal_midpoint;
        if (mpi_rank == 0 && write_flux_audit &&
            midpoint_result.u_flux_audit_valid) {
            u_flux_audit << step << "  "
                         << time / Const::femto << "  "
                         << midpoint_result.u_flux_audit_valid << "  "
                         << midpoint_result.u_flux_audit_rank << "  "
                         << midpoint_result.u_flux_audit_ix << "  "
                         << midpoint_result.u_flux_audit_iv << "  "
                         << midpoint_result.u_flux_audit_imu << "  "
                         << midpoint_result.u_flux_audit_severity << "  "
                         << midpoint_result.u_flux_audit_f0 << "  "
                         << midpoint_result.u_flux_audit_f_low << "  "
                         << midpoint_result.u_flux_audit_f_high << "  "
                         << midpoint_result.u_flux_audit_alpha << "  "
                         << midpoint_result.u_flux_audit_du_div_low << "  "
                         << midpoint_result.u_flux_audit_du_div_high << "  "
                         << midpoint_result.u_flux_audit_du_div_final << "  "
                         << midpoint_result.u_flux_audit_updated << "  "
                         << midpoint_result.u_flux_audit_final_xl_lo << "  "
                         << midpoint_result.u_flux_audit_final_xl_hi << "  "
                         << midpoint_result.u_flux_audit_final_xr_lo << "  "
                         << midpoint_result.u_flux_audit_final_xr_hi << "\n";
            u_flux_audit.flush();
        }
        if (mpi_rank == 0 && midpoint_result.u_flux_audit_valid) {
            const double x_cell =
                (static_cast<double>(midpoint_result.u_flux_audit_ix) + 0.5) *
                sgrid.dx;
            const int region =
                (x_cell < 0.1 * Const::micro ||
                 x_cell > Param::Lx - 0.1 * Const::micro) ? 1 : 0;
            u_origin_audit << step << "  "
                           << time / Const::femto << "  "
                           << midpoint_result.u_flux_audit_rank << "  "
                           << midpoint_result.u_flux_audit_ix << "  "
                           << midpoint_result.u_flux_audit_iv << "  "
                           << midpoint_result.u_flux_audit_imu << "  "
                           << region << "  "
                           << midpoint_result.u_flux_audit_f0 << "  "
                           << midpoint_result.u_flux_audit_f_after_x << "  "
                           << midpoint_result.u_flux_audit_f_low << "  "
                           << midpoint_result.u_flux_audit_f_high << "  "
                           << midpoint_result.u_flux_audit_updated << "  "
                           << midpoint_result.u_flux_audit_lower_flux << "  "
                           << midpoint_result.u_flux_audit_upper_flux << "  "
                           << midpoint_result.u_flux_audit_lower_donor_iv
                           << "  "
                           << midpoint_result.u_flux_audit_upper_donor_iv
                           << "  "
                           << midpoint_result.u_flux_audit_lower_donor_f
                           << "  "
                           << midpoint_result.u_flux_audit_upper_donor_f
                           << "  "
                           << midpoint_result.u_flux_audit_cfl << "  "
                           << midpoint_result.u_flux_audit_alpha << "  "
                           << midpoint_result.u_flux_audit_alpha << "  "
                           << midpoint_result
                                  .u_flux_audit_line_positive_peak << "  "
                           << midpoint_result
                                  .u_flux_audit_line_negative_mass << "  "
                           << midpoint_result.u_flux_audit_allowed_debt
                           << "  "
                           << midpoint_result.u_flux_audit_failure_kind
                           << "\n";
            if (abnormal_midpoint) {
                u_origin_audit.flush();
            }
        }
        if (mpi_rank == 0 && abnormal_midpoint) {
            write_low_order_failure_audit(
                u_low_failure_audit,
                midpoint_result.u_low_failure_audit);
            write_low_order_failure_audit(
                mu_low_failure_audit,
                midpoint_result.mu_low_failure_audit);
        }

        const bool write_periodic_monitors =
            collect_step_diagnostics || abnormal_midpoint ||
            (x_low_input_debt_accepted_step > 0.0) ||
            (x_final_failed_count_step > 0.0) ||
            (mu_final_failed_count_step > 0.0);
        if (write_periodic_monitors) {
            const double bkg_ke_step_end_for_residual =
                bkg_ke_step_start + midpoint_result.delta_ke_bkg;
            local_bkg_energy_residual_step =
                (bkg_ke_step_end_for_residual - bkg_ke_step_start) + W_bkg_E;
            const double local_bkg_energy_values[8] = {
                bkg_ke_step_end_for_residual,
                W_bkg_E,
                local_bkg_energy_residual_step,
                local_bkg_current_diag_step.residual_if_charge,
                local_bkg_current_diag_step.residual_if_ampere,
                local_bkg_current_diag_step.e_dot_j_charge,
                local_bkg_current_diag_step.e_dot_j_energy,
                local_bkg_current_diag_step.e_dot_j_ampere
            };
            double global_bkg_energy_values[8] = {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            };
            MPI_Allreduce(local_bkg_energy_values, global_bkg_energy_values, 8,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            const double local_current_max_values[5] = {
                local_bkg_current_diag_step.max_abs_j_charge,
                local_bkg_current_diag_step.max_abs_j_energy,
                local_bkg_current_diag_step.max_abs_j_ampere,
                local_bkg_current_diag_step.max_abs_j_charge_minus_ampere,
                local_bkg_current_diag_step.max_abs_j_energy_minus_ampere
            };
            double global_current_max_values[5] = {
                0.0, 0.0, 0.0, 0.0, 0.0
            };
            MPI_Allreduce(local_current_max_values, global_current_max_values,
                          5, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            global_bkg_current_diag_step.residual_if_charge =
                global_bkg_energy_values[3];
            global_bkg_current_diag_step.residual_if_ampere =
                global_bkg_energy_values[4];
            global_bkg_current_diag_step.e_dot_j_charge =
                global_bkg_energy_values[5];
            global_bkg_current_diag_step.e_dot_j_energy =
                global_bkg_energy_values[6];
            global_bkg_current_diag_step.e_dot_j_ampere =
                global_bkg_energy_values[7];
            global_bkg_current_diag_step.max_abs_j_charge =
                global_current_max_values[0];
            global_bkg_current_diag_step.max_abs_j_energy =
                global_current_max_values[1];
            global_bkg_current_diag_step.max_abs_j_ampere =
                global_current_max_values[2];
            global_bkg_current_diag_step.max_abs_j_charge_minus_ampere =
                global_current_max_values[3];
            global_bkg_current_diag_step.max_abs_j_energy_minus_ampere =
                global_current_max_values[4];
            const double global_dke_bkg_step =
                global_bkg_energy_values[0] - global_bkg_ke_step_start;
            const double global_W_bkg_E = global_bkg_energy_values[1];
            bkg_energy_residual_step = global_bkg_energy_values[2];
            const double bkg_energy_residual_den =
                std::max(std::max(std::fabs(global_dke_bkg_step),
                                  std::fabs(global_W_bkg_E)),
                         std::max(1.0,
                                  1.0e-12 *
                                  std::fabs(global_bkg_ke_step_start)));
            bkg_energy_relative_residual_step =
                std::fabs(bkg_energy_residual_step) /
                bkg_energy_residual_den;
            cumulative_bkg_energy_residual +=
                bkg_energy_relative_residual_step;
        }
        if (mpi_rank == 0 && write_periodic_monitors) {
            bkg_energy_monitor << step << "  "
                               << time / Const::femto << "  "
                               << x_limiter_active_fraction_step << "  "
                               << x_limiter_min_alpha_step << "  "
                               << bkg_energy_residual_step << "  "
                               << bkg_energy_relative_residual_step << "  "
                               << cumulative_bkg_energy_residual << "  "
                               << x_limiter_active_fraction_core_step << "  "
                               << x_limiter_active_fraction_boundary_step << "  "
                               << x_limiter_min_alpha_core_step << "  "
                               << x_limiter_min_alpha_boundary_step << "  "
                               << x_negative_mass_before_repair_step << "  "
                               << x_mass_added_by_positivity_repair_step << "  "
                               << floor_repair_mass_step << "  "
                               << floor_repair_energy_step << "  "
                               << floor_repair_core_fraction_step << "  "
                               << global_bkg_current_diag_step.max_abs_j_charge << "  "
                               << global_bkg_current_diag_step.max_abs_j_energy << "  "
                               << global_bkg_current_diag_step.max_abs_j_ampere << "  "
                               << global_bkg_current_diag_step.max_abs_j_charge_minus_ampere << "  "
                               << global_bkg_current_diag_step.max_abs_j_energy_minus_ampere << "  "
                               << global_bkg_current_diag_step.e_dot_j_charge << "  "
                               << global_bkg_current_diag_step.e_dot_j_energy << "  "
                               << global_bkg_current_diag_step.e_dot_j_ampere << "  "
                               << global_bkg_current_diag_step.residual_if_charge << "  "
                               << global_bkg_current_diag_step.residual_if_ampere << "  "
                               << positivity_mass_defect_step << "  "
                               << positivity_energy_defect_step << "  "
                               << u_limiter_mass_delta_step << "  "
                               << u_limiter_momentum_delta_step << "  "
                               << u_limiter_energy_delta_step << "  "
                               << x_limiter_mass_delta_step << "  "
                               << x_limiter_energy_delta_step << "  "
                               << mu_low_u_alpha_min_step << "  "
                               << mu_low_u_limiter_active_fraction_step << "  "
                               << mu_low_u_energy_delta_step << "  "
                               << mu_low_u_alpha_min_boundary_step << "  "
                               << mu_low_u_alpha_min_core_step << "  "
                               << mu_low_u_limiter_active_fraction_boundary_step << "  "
                               << mu_low_u_limiter_active_fraction_core_step << "  "
                               << mu_low_u_energy_delta_boundary_step << "  "
                               << mu_low_u_energy_delta_core_step << "  "
                               << mu_low_u_u_eff0_step << "  "
                               << mu_low_u_moment_weight0_step << "  "
                               << mu_low_u_mu_flux_scale0_step << "  "
                               << mu_low_u_half_dt_inv_shell0_step << "  "
                               << mu_low_u_dimless_scale0_step << "  "
                               << mu_low_u_endpoint_flux_max_step << "  "
                               << remap_active_fraction_step << "  "
                               << remap_cell_count_step << "  "
                               << low_u_subcycle_active_fraction_step << "  "
                               << low_u_average_subcycles_step << "  "
                               << low_u_max_subcycles_step << "  "
                               << u_force_alpha_min_step << "  "
                               << u_force_alpha_active_frac_step << "  "
                               << coupled_iter_step << "  "
                               << coupled_residual_E_step << "  "
                               << coupled_residual_J_bkg_step << "  "
                               << coupled_residual_J_beam_step << "  "
                               << coupled_residual_bkg_mass_step << "  "
                               << coupled_residual_beam_continuity_step << "\n";
            bkg_energy_monitor.flush();
            x_low_monitor << step << "  "
                          << time / Const::femto << "  "
                          << x_low_input_min_f_step << "  "
                          << x_low_max_cfl_step << "  "
                          << x_low_output_min_f_step << "  "
                          << x_low_failed_count_step << "\n";
            x_low_monitor.flush();
            x_final_monitor << step << "  "
                            << time / Const::femto << "  "
                            << x_final_min_f_step << "  "
                            << x_final_failed_count_step << "  "
                            << x_final_failed_max_debt_step << "  "
                            << x_final_worst_ix_step << "  "
                            << x_final_worst_iv_step << "  "
                            << x_final_worst_imu_step << "  "
                            << x_final_failure_region_step << "  "
                            << x_final_core_failed_count_step << "  "
                            << x_final_boundary_failed_count_step << "\n";
            x_final_monitor.flush();
            mu_final_monitor << step << "  "
                             << time / Const::femto << "  "
                             << mu_final_min_f_step << "  "
                             << mu_final_failed_count_step << "  "
                             << mu_final_failed_max_debt_step << "  "
                             << mu_final_worst_ix_step << "  "
                             << mu_final_worst_iv_step << "  "
                             << mu_final_worst_imu_step << "  "
                             << mu_final_failure_region_step << "  "
                             << mu_final_core_failed_count_step << "  "
                             << mu_final_boundary_failed_count_step << "  "
                             << midpoint_result.flux_pos[2].alpha_active_fraction << "  "
                             << midpoint_result.flux_pos[2].alpha_min << "  "
                             << midpoint_result.flux_pos[2].alpha_core_fraction << "  "
                             << midpoint_result.flux_pos[2].alpha_boundary_fraction << "  "
                             << mu_low_u_alpha_min_core_step << "  "
                             << mu_low_u_alpha_min_boundary_step << "  "
                             << midpoint_result.flux_defect[2].energy_defect << "  "
                             << midpoint_result.flux_defect[2].mass_defect << "\n";
            mu_final_monitor.flush();
        }
        if (collect_step_diagnostics) {
            diag.write_bkg_stage_negativity(
                step, time, coupled_iter_step,
                midpoint_result.soft_unconverged,
                midpoint_result.stage_min_f,
                midpoint_result.stage_neg_mass,
                midpoint_result.stage_neg_cell_count,
                midpoint_result.stage_low_u_neg_mass,
                midpoint_result.stage_core_low_u_min_f,
                mpi_rank);
            diag.write_bkg_stage_by_u_diagnostics(
                step, time, coupled_iter_step,
                midpoint_result.soft_unconverged,
                midpoint_result.stage_min_f_core_by_u,
                midpoint_result.stage_neg_mass_core_by_u,
                midpoint_result.stage_neg_cell_count_core_by_u,
                midpoint_result.stage_min_f_boundary_by_u,
                midpoint_result.stage_neg_mass_boundary_by_u,
                midpoint_result.stage_neg_cell_count_boundary_by_u,
                mpi_rank);
            diag.write_bkg_low_u_divergence_diagnostics(
                step, time, coupled_iter_step,
                midpoint_result.low_u_neg_added_by_div,
                mpi_rank);

            // 7.1.6: per-direction flux diagnostics
            {
                double fp_min_before[3], fp_min_low[3], fp_min_final[3];
                double fp_low_failed[3], fp_alpha_active[3], fp_alpha_min[3];
                double fp_alpha_core[3], fp_alpha_boundary[3];
                double fp_neg_mass_prev[3];
                double fd_mass[3], fd_momentum[3], fd_energy[3];
                double fd_bound_mass[3], fd_bound_energy[3];
                for (int d = 0; d < 3; ++d) {
                    const auto& fp = midpoint_result.flux_pos[d];
                    fp_min_before[d]   = fp.min_f_before;
                    fp_min_low[d]      = fp.min_f_low;
                    fp_min_final[d]    = fp.min_f_final;
                    fp_low_failed[d]   = fp.low_order_failed_count;
                    fp_alpha_active[d] = fp.alpha_active_fraction;
                    fp_alpha_min[d]    = fp.alpha_min;
                    fp_alpha_core[d]   = fp.alpha_core_fraction;
                    fp_alpha_boundary[d] = fp.alpha_boundary_fraction;
                    fp_neg_mass_prev[d]  = fp.negative_mass_prevented;
                    const auto& fd = midpoint_result.flux_defect[d];
                    fd_mass[d]      = fd.mass_defect;
                    fd_momentum[d]  = fd.momentum_defect;
                    fd_energy[d]    = fd.energy_defect;
                    fd_bound_mass[d]  = fd.boundary_mass_loss;
                    fd_bound_energy[d] = fd.boundary_energy_loss;
                }
                diag.write_flux_positivity_diagnostics(
                    step, time,
                    fp_min_before, fp_min_low, fp_min_final,
                    fp_low_failed, fp_alpha_active, fp_alpha_min,
                    fp_alpha_core, fp_alpha_boundary, fp_neg_mass_prev,
                    mpi_rank);
                diag.write_stage_flux_defect_diagnostics(
                    step, time,
                    fd_mass, fd_momentum, fd_energy,
                    fd_bound_mass, fd_bound_energy,
                    mpi_rank);
            }

            FNegativitySnapshotDiagnostics f_neg_snapshot;
            compute_f_negativity_snapshot_diagnostics(
                bkg_e, sgrid, f_neg_snapshot);

            if (mpi_rank == 0) {
                // f-negativity monitor: final accepted snapshot at the
                // configured step-diagnostics cadence.
                std::ofstream f_neg_monitor;
                f_neg_monitor.open("output/f_negativity_monitor.dat",
                                   std::ios::app);
                f_neg_monitor << step << "  "
                              << time / Const::femto << "  "
                              << 1 << "  "
                              << midpoint_result.state_advanced << "  "
                              << midpoint_result.soft_unconverged << "  "
                              << f_neg_snapshot.min_f << "  "
                              << f_neg_snapshot.neg_ratio_max << "  "
                              << f_neg_snapshot.neg_mass_total << "  "
                              << f_neg_snapshot.neg_cell_count << "  "
                              << f_neg_snapshot.x_worst << "  "
                              << f_neg_snapshot.u_worst << "  "
                              << f_neg_snapshot.mu_worst << "\n";
                f_neg_monitor.close();
            }
        }
        if (collect_step_diagnostics || abnormal_midpoint) {
            const double region_widths[2] = {
                0.1 * Const::micro,
                0.2 * Const::micro
            };
            for (int ir = 0; ir < 2; ++ir) {
                BoundaryCoreDiagnostics region_diag;
                compute_boundary_core_diagnostics(
                    bkg_e, fields, sgrid, dt, region_widths[ir],
                    mpi_rank, mpi_size, region_diag);
                const double neg_mass_total_region =
                    region_diag.neg_mass_boundary + region_diag.neg_mass_core;
                const double neg_mass_core_fraction =
                    (neg_mass_total_region > 0.0)
                    ? region_diag.neg_mass_core / neg_mass_total_region
                    : 0.0;
                const long long neg_cell_total_region =
                    region_diag.neg_cell_count_boundary
                  + region_diag.neg_cell_count_core;
                const double neg_cell_core_fraction =
                    (neg_cell_total_region > 0)
                    ? static_cast<double>(region_diag.neg_cell_count_core) /
                      static_cast<double>(neg_cell_total_region)
                    : 0.0;
                const double abs_limiter_boundary =
                    midpoint_result
                        .region_abs_u_limiter_energy_boundary[ir];
                const double abs_limiter_core =
                    midpoint_result
                        .region_abs_u_limiter_energy_core[ir];
                const double abs_limiter_total =
                    abs_limiter_boundary + abs_limiter_core;
                const double limiter_core_fraction =
                    (abs_limiter_total > 0.0)
                    ? abs_limiter_core / abs_limiter_total
                    : 0.0;
                const double j_core_to_boundary_ratio =
                    (region_diag.max_abs_J_bkg_boundary > 0.0)
                    ? region_diag.max_abs_J_bkg_core /
                      region_diag.max_abs_J_bkg_boundary
                    : 0.0;
                const double core_limiter_relative_to_work =
                    abs_limiter_core /
                    std::max(std::fabs(region_diag.W_bkg_E_core), 1.0);
                const char* conclusion =
                    classify_boundary_core_state(
                        neg_mass_core_fraction,
                        neg_cell_core_fraction,
                        limiter_core_fraction,
                        j_core_to_boundary_ratio);
                if (mpi_rank == 0) {
                    boundary_core_monitor
                        << step << "  "
                        << time / Const::femto << "  "
                        << 1 << "  "
                        << midpoint_result.state_advanced << "  "
                        << midpoint_result.soft_unconverged << "  "
                        << region_widths[ir] / Const::micro << "  "
                        << region_diag.neg_mass_boundary << "  "
                        << region_diag.neg_mass_core << "  "
                        << neg_mass_core_fraction << "  "
                        << region_diag.neg_cell_count_boundary << "  "
                        << region_diag.neg_cell_count_core << "  "
                        << neg_cell_core_fraction << "  "
                        << region_diag.min_f_boundary << "  "
                        << region_diag.min_f_core << "  "
                        << midpoint_result
                               .region_u_limiter_energy_boundary[ir] << "  "
                        << midpoint_result
                               .region_u_limiter_energy_core[ir] << "  "
                        << abs_limiter_boundary << "  "
                        << abs_limiter_core << "  "
                        << limiter_core_fraction << "  "
                        << core_limiter_relative_to_work << "  "
                        << midpoint_result
                               .region_limiter_active_fraction_boundary[ir]
                        << "  "
                        << midpoint_result
                               .region_limiter_active_fraction_core[ir]
                        << "  "
                        << midpoint_result.limiter_min_alpha_boundary
                        << "  "
                        << midpoint_result.limiter_min_alpha_core << "  "
                        << region_diag.max_abs_J_bkg_boundary << "  "
                        << region_diag.max_abs_J_bkg_core << "  "
                        << j_core_to_boundary_ratio << "  "
                        << region_diag.rms_J_bkg_boundary << "  "
                        << region_diag.rms_J_bkg_core << "  "
                        << region_diag.mean_abs_J_bkg_boundary << "  "
                        << region_diag.mean_abs_J_bkg_core << "  "
                        << region_diag.max_abs_Ex_boundary << "  "
                        << region_diag.max_abs_Ex_core << "  "
                        << region_diag.W_bkg_E_boundary << "  "
                        << region_diag.W_bkg_E_core << "  "
                        << region_diag.anomaly_inward_E/Const::micro << "  "
                        << region_diag.anomaly_inward_J/Const::micro << "  "
                        << region_diag.anomaly_inward_n/Const::micro << "  "
                        << region_diag.anomaly_inward_Mneg/Const::micro << "  "
                        << conclusion << "\n";
                }
            }
            if (mpi_rank == 0) {
                boundary_core_monitor.flush();
            }
        }
        const bool midpoint_stage_values_are_valid =
            midpoint_result.stage_min_f.size() >= 3 &&
            midpoint_result.stage_neg_mass.size() >= 3 &&
            midpoint_result.stage_neg_cell_count.size() >= 3 &&
            std::isfinite(midpoint_result.stage_min_f[0]) &&
            std::isfinite(midpoint_result.stage_min_f[1]) &&
            std::isfinite(midpoint_result.stage_min_f[2]);
        if (midpoint_stage_values_are_valid) {
            for (int istage = 0; istage < 3; ++istage) {
                last_accepted_stage_min_f[istage] =
                    midpoint_result.stage_min_f[static_cast<size_t>(istage)];
                last_accepted_stage_neg_mass[istage] =
                    midpoint_result.stage_neg_mass[static_cast<size_t>(istage)];
                last_accepted_stage_neg_cell_count[istage] =
                    midpoint_result
                        .stage_neg_cell_count[static_cast<size_t>(istage)];
            }
            last_accepted_stage_valid = true;
        }
        if (bkg_e.collisions_enabled) {
            trace_progress(config, mpi_rank, step, "before collisions");
            collision_energy_step +=
                collision.apply(bkg_e, dt, Param::dens, Param::temperature_e,
                                Const::me, 1.0, 1.0);
            collision_energy_step +=
                collision.apply(bkg_e, dt, Param::dens / Param::Z_ion,
                                Param::temperature_i, Param::mass_ion,
                                (double)Param::Z_ion, 1.0);
            cumulative_collision_energy_delta += collision_energy_step;
            trace_progress(config, mpi_rank, step, "after collisions");
            moments_current = false;
        }

        trace_progress(config, mpi_rank, step, "before end sync");
        if (!moments_current) {
            bkg_e.compute_moments();
            moments_current = true;
        }
        fields.set_charge_density(bkg_e, beam.density, ion_density_profile);
        fields.update_gauss_residual_diagnostics(mpi_rank, mpi_size);
        trace_progress(config, mpi_rank, step, "after end sync");

        net_nb_change_step = beam.last_injected_number()
                           - beam.last_outflow_number();

        if (collect_step_diagnostics) {
            const double bkg_ke_step_end = bkg_e.total_kinetic_energy();
            const double beam_ke_step_end = beam.total_kinetic_energy();
            const double field_energy_step_end = fields.total_energy();
            dke_bkg_step = bkg_ke_step_end - bkg_ke_step_start;
            dE_field_step = field_energy_step_end - field_energy_step_start;
            E_src_in_step = beam.last_injected_energy();
            E_src_out_step = beam.last_outflow_energy();
            const double total_energy_delta =
                (bkg_ke_step_end + beam_ke_step_end + field_energy_step_end) -
                (bkg_ke_step_start + beam_ke_step_start + field_energy_step_start);
            E_balance_step =
                total_energy_delta - E_src_in_step + E_src_out_step
                - collision_energy_step;

            diag.write_step_diagnostics(step, time,
                                        midpoint_result.soft_unconverged,
                                        bkg_e, beam, fields,
                                        sgrid, mpi_rank, mpi_size,
                                        nsub_v1, nsub_mu1,
                                        nsub_v2, nsub_mu2,
                                        loss_v1, loss_mu1,
                                        loss_v2, loss_mu2,
                                        loss_v1_low, loss_v1_high,
                                        loss_v2_low, loss_v2_high,
                                        net_nb_change_step,
                                        collision_energy_step,
                                        cumulative_collision_energy_delta,
                                        dke_bkg_step, dke_beam_push,
                                        dE_field_step, W_bkg_E, W_beam_E,
                                        v_mass_error_step,
                                        mu_mass_error_step,
                                        v_momentum_delta_step,
                                        mu_momentum_delta_step,
                                        v_energy_delta_step,
                                        mu_energy_delta_step,
                                        E_src_in_step,
                                        E_src_out_step,
                                        E_balance_step,
                                        x_limiter_active_fraction_step,
                                        x_limiter_min_alpha_step,
                                        local_bkg_energy_residual_step,
                                        local_bkg_current_diag_step.max_abs_j_charge,
                                        local_bkg_current_diag_step.max_abs_j_energy,
                                        local_bkg_current_diag_step.max_abs_j_ampere,
                                        local_bkg_current_diag_step.max_abs_j_charge_minus_ampere,
                                        local_bkg_current_diag_step.max_abs_j_energy_minus_ampere,
                                        local_bkg_current_diag_step.e_dot_j_charge,
                                        local_bkg_current_diag_step.e_dot_j_energy,
                                        local_bkg_current_diag_step.e_dot_j_ampere,
                                        local_bkg_current_diag_step.residual_if_charge,
                                        local_bkg_current_diag_step.residual_if_ampere,
                                        midpoint_result.boundary_force_Cu_max,
                                        midpoint_result.boundary_force_Cmu_max,
                                        midpoint_result.boundary_force_nsub_max,
                                        midpoint_result.boundary_force_remap_cell_count,
                                        midpoint_result.boundary_mu_low_L1_before,
                                        midpoint_result.boundary_mu_low_L1_after,
                                        midpoint_result.boundary_mu_high_L1_after,
                                        midpoint_result.J_bkg_neg_boundary,
                                        midpoint_result.delta_E_neg_boundary,
                                        midpoint_result.boundary_force_remap_mass_loss,
                                        midpoint_result.boundary_force_remap_energy_loss,
                                        midpoint_result.alpha_interface_BQ_min,
                                        midpoint_result.alpha_interface_QC_min,
                                        midpoint_result.interface_BQ_flux,
                                        midpoint_result.interface_BQ_high_correction,
                                        midpoint_result.interface_QC_flux_into_core,
                                        midpoint_result.interface_QC_high_correction_into_core,
                                        midpoint_result.boundary_energy_diagnostic_invalid,
                                        coupled_iter_step,
                                        coupled_residual_E_step,
                                        coupled_residual_J_bkg_step,
                                        coupled_residual_J_beam_step,
                                        max_loss_u_high_step,
                                        x_at_max_loss_u_high_step,
                                        f_u_max_x_step,
                                        integral_f_u_gt_8_x_step);
        }

        if (step % stdout_freq == 0) {
            diag.write_scalars(time, step, bkg_e, beam, fields,
                               cumulative_collision_energy_delta,
                               mpi_rank, mpi_size);
            if (mpi_rank == 0) {
                printf("Step %d / %d, t = %.4f fs\n", step, nsteps, time / Const::femto);
            }
        }

        if (time >= next_snapshot) {
            write_snapshot(diag, time, bkg_e, beam, fields, ion_density_profile,
                           sgrid, mpi_rank, mpi_size,
                           config.enable_full_fe_output,
                           &latest_bkg_energy_current_face,
                           &latest_bkg_ampere_current_face);
            last_snapshot_step = step;
            next_snapshot += Param::dt_snapshot;
        }
    }

    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    fields.update_gauss_residual_diagnostics(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(nsteps, Param::t_end, "final", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(Param::t_end, nsteps, bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    if (last_snapshot_step != nsteps) {
        write_snapshot(diag, Param::t_end, bkg_e, beam, fields, ion_density_profile,
                       sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                       &latest_bkg_energy_current_face,
                       &latest_bkg_ampere_current_face);
    }

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Simulation complete: t = %.1f fs, %d steps\n",
               Param::t_end / Const::femto, nsteps);
        printf("============================================================\n");
    }

    MPI_Finalize();
    return 0;
}
