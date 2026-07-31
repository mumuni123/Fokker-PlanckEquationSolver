#include "parameters.h"
#include "grid.h"
#include "species.h"
#include "maxwell.h"
#include "collision.h"
#include "diagnostics.h"
#include "beam_pic.h"
#include "vlasov_ampere_midpoint.h"
#include "config.h"
#include "runtime_options.h"
#include "checkpoint.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <vector>
#if !defined(_WIN32)
#include <sys/resource.h>
#endif

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

bool read_beam_ledger_total(const std::string& filename, double& n_in_total,
                            double& current_impulse_total)
{
    std::ifstream input(filename.c_str());
    std::string key;
    double value = 0.0;
    bool have_n_in = false;
    bool have_impulse = false;
    while (input >> key >> value) {
        if (key == "N_in_total") {
            n_in_total = value;
            have_n_in = true;
        } else if (key == "injected_current_impulse_total") {
            current_impulse_total = value;
            have_impulse = true;
        }
    }
    return have_n_in && have_impulse && std::isfinite(n_in_total) &&
        std::isfinite(current_impulse_total);
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
                    const std::vector<double>* bkg_ampere_current_face = 0,
                    bool bkg_energy_current_valid = false)
{
    fields.compute_potential(mpi_rank, mpi_size);
    diag.write_fields(time, fields, sgrid, mpi_rank, mpi_size);
    diag.write_current_density(time, bkg_e, beam, sgrid, mpi_rank, mpi_size,
                               bkg_energy_current_face,
                               bkg_ampere_current_face,
                               bkg_energy_current_valid);
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

enum FixedPairTerm {
    FIXED_PAIR_LOW = 0,
    FIXED_PAIR_X_HIGH = 1,
    FIXED_PAIR_U_CENTER = 2,
    FIXED_PAIR_U_UPWIND = 3,
    FIXED_PAIR_X_FCT = 4,
    FIXED_PAIR_U_FCT = 5,
    FIXED_PAIR_RECONSTRUCTED_FINAL = 6,
    FIXED_PAIR_FINAL = 7,
    FIXED_PAIR_RECONSTRUCTION_ERROR = 8,
    FIXED_PAIR_TERM_COUNT = 9
};

struct FixedPairTermSummary {
    double signed_sum;
    double absolute_sum;
    double linf;
    int max_global_face;
};

// Gather only owned periodic faces (global 0..nx-1) and write the fixed-state
// pairing ledger. The x=L alias is intentionally excluded: it is audited in
// fixed_midpoint_operator_audit.result, but is not a second physical face.
struct FixedMidpointFacePairingAudit {
    int row_count;
    int expected_row_count;
    int unique_face_count;
    int complete_unique_coverage;
    int ledger_match;
    double r_couple_from_unique_faces;
    double reconstruction_error;
    std::array<FixedPairTermSummary, FIXED_PAIR_TERM_COUNT> global_terms;
    std::array<std::array<FixedPairTermSummary, FIXED_PAIR_TERM_COUNT>, 6>
        region_terms;
};

FixedMidpointFacePairingAudit write_fixed_midpoint_face_pairing(
    const std::string& filename, long long step, double time_s, double dt_s,
    const VlasovAmpereMidpointSolver::Result& result,
    const EMFields& fields_start, const EMFields& fields_end_guess,
    const SpatialGrid& sg, int mpi_rank, int mpi_size, bool append)
{
    FixedMidpointFacePairingAudit audit = {};
    audit.expected_row_count = sg.nx_global;
    const int nxl = sg.nx_local;
    // E_mid, JN_low/high/final, G*JE_low/center/high/final.  There is no
    // JN_center column: x transport has no independently constructed center
    // current, so exposing its high-order alias would be misleading.
    const int columns = 8;
    std::vector<int> local_faces(static_cast<size_t>(nxl), 0);
    std::vector<double> local_face_values(static_cast<size_t>(nxl * columns), 0.0);
    for (int iface = 0; iface < nxl; ++iface) {
        local_faces[static_cast<size_t>(iface)] = sg.ix_start + iface;
        const size_t base = static_cast<size_t>(iface * columns);
        local_face_values[base] = 0.5 * (fields_start.Ex_face[iface] +
                                         fields_end_guess.Ex_face[iface]);
        local_face_values[base + 1] = result.j_bkg_face_low_mid[iface];
        local_face_values[base + 2] = result.j_bkg_face_high_mid[iface];
        local_face_values[base + 3] = result.j_bkg_face_mid[iface];
        local_face_values[base + 4] = result.j_bkg_energy_low_debug_face[iface];
        local_face_values[base + 5] = result.j_bkg_energy_center_debug_face[iface];
        local_face_values[base + 6] = result.j_bkg_energy_high_debug_face[iface];
        local_face_values[base + 7] = result.j_bkg_energy_debug_face[iface];
    }

    std::vector<int> counts;
    if (mpi_rank == 0) counts.assign(static_cast<size_t>(mpi_size), 0);
    MPI_Gather(&nxl, 1, MPI_INT, mpi_rank == 0 ? counts.data() : 0, 1,
               MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> displacements;
    std::vector<int> value_counts;
    std::vector<int> value_displacements;
    std::vector<int> gathered_faces;
    std::vector<double> gathered_face_values;
    if (mpi_rank == 0) {
        displacements.assign(static_cast<size_t>(mpi_size), 0);
        value_counts.assign(static_cast<size_t>(mpi_size), 0);
        value_displacements.assign(static_cast<size_t>(mpi_size), 0);
        int total = 0;
        for (int rank = 0; rank < mpi_size; ++rank) {
            displacements[static_cast<size_t>(rank)] = total;
            value_displacements[static_cast<size_t>(rank)] = total * columns;
            value_counts[static_cast<size_t>(rank)] = counts[static_cast<size_t>(rank)] * columns;
            total += counts[static_cast<size_t>(rank)];
        }
        gathered_faces.assign(static_cast<size_t>(total), -1);
        gathered_face_values.assign(static_cast<size_t>(total * columns), 0.0);
    }
    MPI_Gatherv(local_faces.data(), nxl, MPI_INT,
                mpi_rank == 0 ? gathered_faces.data() : 0,
                mpi_rank == 0 ? counts.data() : 0,
                mpi_rank == 0 ? displacements.data() : 0,
                MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_face_values.data(), nxl * columns, MPI_DOUBLE,
                mpi_rank == 0 ? gathered_face_values.data() : 0,
                mpi_rank == 0 ? value_counts.data() : 0,
                mpi_rank == 0 ? value_displacements.data() : 0,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        double values[2] = {0.0, 0.0};
        int flags[5] = {0, 0, 0, 0, 0};
        MPI_Bcast(values, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(flags, 5, MPI_INT, 0, MPI_COMM_WORLD);
        audit.r_couple_from_unique_faces = values[0];
        audit.reconstruction_error = values[1];
        audit.row_count = flags[0];
        audit.expected_row_count = flags[1];
        audit.unique_face_count = flags[2];
        audit.complete_unique_coverage = flags[3] != 0;
        audit.ledger_match = flags[4] != 0;
        return audit;
    }

    const int nx_global = sg.nx_global;
    std::vector<double> face_values(static_cast<size_t>(nx_global * columns), 0.0);
    std::vector<int> owner(static_cast<size_t>(nx_global), -1);
    std::vector<int> local_face(static_cast<size_t>(nx_global), -1);
    std::vector<int> seen(static_cast<size_t>(nx_global), 0);
    for (int rank = 0; rank < mpi_size; ++rank) {
        const int begin = displacements[static_cast<size_t>(rank)];
        const int count = counts[static_cast<size_t>(rank)];
        for (int local = 0; local < count; ++local) {
            const int packed = begin + local;
            const int global_face = gathered_faces[static_cast<size_t>(packed)];
            if (global_face < 0 || global_face >= nx_global) continue;
            const size_t destination = static_cast<size_t>(global_face * columns);
            const size_t source = static_cast<size_t>(packed * columns);
            for (int column = 0; column < columns; ++column)
                face_values[destination + static_cast<size_t>(column)] =
                    gathered_face_values[source + static_cast<size_t>(column)];
            owner[static_cast<size_t>(global_face)] = rank;
            local_face[static_cast<size_t>(global_face)] = local;
            ++seen[static_cast<size_t>(global_face)];
        }
    }

    std::ofstream out(filename.c_str(), append ? std::ios::app : std::ios::out);
    if (!append) {
        out << "# unique periodic face ledger; x has no independently "
            << "constructed JN_center layer.\n";
        out << "# step time_fs global_face x_face E_mid JN_low JN_high JN_final "
            << "GstarJE_low GstarJE_center GstarJE_high GstarJE_final "
            << "r_low delta_r_x_high delta_r_u_center delta_r_u_upwind "
            << "delta_r_x_FCT delta_r_u_FCT r_reconstructed_final r_final "
            << "reconstruction_error owner_rank local_face is_global_seam "
            << "is_periodic_alias\n";
    }
    out << std::scientific << std::setprecision(17);
    const char* const term_names[FIXED_PAIR_TERM_COUNT] = {
        "R_pair_low", "delta_R_x_high", "delta_R_u_center",
        "delta_R_u_upwind", "delta_R_x_FCT", "delta_R_u_FCT",
        "R_pair_reconstructed_final", "R_pair_final",
        "R_pair_reconstruction_error"};
    const auto region_for_face = [&](double x_face) {
        if (x_face < 0.2 * Const::micro) return 0;
        if (x_face < 0.8 * Const::micro) return 1;
        if (x_face <= result.coupling_wave_core_end_m) return 2;
        if (x_face < 7.2 * Const::micro) return 3;
        if (x_face < 7.8 * Const::micro) return 4;
        return 5;
    };
    const auto accumulate_term = [](FixedPairTermSummary& summary,
                                    double value, int global_face) {
        summary.signed_sum += value;
        summary.absolute_sum += std::fabs(value);
        if (std::fabs(value) > summary.linf) {
            summary.linf = std::fabs(value);
            summary.max_global_face = global_face;
        }
    };
    for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
        audit.global_terms[static_cast<size_t>(term)].max_global_face = -1;
        for (int region = 0; region < 6; ++region)
            audit.region_terms[static_cast<size_t>(region)]
                [static_cast<size_t>(term)].max_global_face = -1;
    }
    for (int global_face = 0; global_face < nx_global; ++global_face) {
        const size_t base = static_cast<size_t>(global_face * columns);
        const double e_mid = face_values[base];
        const double jn_low = face_values[base + 1];
        const double jn_high = face_values[base + 2];
        const double jn_final = face_values[base + 3];
        const double gstar_low = face_values[base + 4];
        const double gstar_center = face_values[base + 5];
        const double gstar_high = face_values[base + 6];
        const double gstar_final = face_values[base + 7];
        const double factor = dt_s * sg.dx * e_mid;
        std::array<double, FIXED_PAIR_TERM_COUNT> terms = {};
        terms[FIXED_PAIR_LOW] = factor * (jn_low - gstar_low);
        terms[FIXED_PAIR_X_HIGH] = factor * (jn_high - jn_low);
        terms[FIXED_PAIR_U_CENTER] = -factor * (gstar_center - gstar_low);
        terms[FIXED_PAIR_U_UPWIND] = -factor * (gstar_high - gstar_center);
        terms[FIXED_PAIR_X_FCT] = factor * (jn_final - jn_high);
        terms[FIXED_PAIR_U_FCT] = -factor * (gstar_final - gstar_high);
        terms[FIXED_PAIR_RECONSTRUCTED_FINAL] = terms[FIXED_PAIR_LOW] +
            terms[FIXED_PAIR_X_HIGH] + terms[FIXED_PAIR_U_CENTER] +
            terms[FIXED_PAIR_U_UPWIND] + terms[FIXED_PAIR_X_FCT] +
            terms[FIXED_PAIR_U_FCT];
        terms[FIXED_PAIR_FINAL] = factor * (jn_final - gstar_final);
        terms[FIXED_PAIR_RECONSTRUCTION_ERROR] =
            terms[FIXED_PAIR_RECONSTRUCTED_FINAL] - terms[FIXED_PAIR_FINAL];
        const double x_face = sg.x_min + static_cast<double>(global_face) * sg.dx;
        const int region = region_for_face(x_face);
        for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
            accumulate_term(audit.global_terms[static_cast<size_t>(term)],
                            terms[static_cast<size_t>(term)], global_face);
            accumulate_term(audit.region_terms[static_cast<size_t>(region)]
                            [static_cast<size_t>(term)],
                            terms[static_cast<size_t>(term)], global_face);
        }
        out << step << " " << time_s / Const::femto << " " << global_face << " "
            << x_face << " " << e_mid << " "
            << jn_low << " " << jn_high << " " << jn_final << " "
            << gstar_low << " " << gstar_center << " " << gstar_high << " "
            << gstar_final << " " << terms[FIXED_PAIR_LOW] << " "
            << terms[FIXED_PAIR_X_HIGH] << " "
            << terms[FIXED_PAIR_U_CENTER] << " "
            << terms[FIXED_PAIR_U_UPWIND] << " "
            << terms[FIXED_PAIR_X_FCT] << " "
            << terms[FIXED_PAIR_U_FCT] << " "
            << terms[FIXED_PAIR_RECONSTRUCTED_FINAL] << " "
            << terms[FIXED_PAIR_FINAL] << " "
            << terms[FIXED_PAIR_RECONSTRUCTION_ERROR] << " "
            << owner[static_cast<size_t>(global_face)] << " "
            << local_face[static_cast<size_t>(global_face)] << " "
            << (global_face == 0 ? 1 : 0) << " 0\n";
    }
    const double r_couple_from_unique_faces =
        audit.global_terms[FIXED_PAIR_FINAL].signed_sum;
    out << "# R_couple_from_unique_faces " << r_couple_from_unique_faces
        << " stage5_R_couple " << result.stage5_r_couple
        << " reconstruction_error "
        << (r_couple_from_unique_faces - result.stage5_r_couple) << "\n";
    for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
        const FixedPairTermSummary& summary =
            audit.global_terms[static_cast<size_t>(term)];
        out << "# global " << term_names[term]
            << " signed_sum " << summary.signed_sum
            << " absolute_sum " << summary.absolute_sum
            << " Linf " << summary.linf
            << " max_global_face " << summary.max_global_face << "\n";
    }
    for (int region = 0; region < 6; ++region) {
        for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
            const FixedPairTermSummary& summary =
                audit.region_terms[static_cast<size_t>(region)]
                    [static_cast<size_t>(term)];
            out << "# region " << region << " " << term_names[term]
                << " signed_sum " << summary.signed_sum
                << " absolute_sum " << summary.absolute_sum
                << " Linf " << summary.linf
                << " max_global_face " << summary.max_global_face << "\n";
        }
    }
    audit.row_count = nx_global;
    audit.unique_face_count = 0;
    audit.complete_unique_coverage = 1;
    for (int global_face = 0; global_face < nx_global; ++global_face) {
        if (seen[static_cast<size_t>(global_face)] == 1) ++audit.unique_face_count;
        else audit.complete_unique_coverage = 0;
    }
    audit.r_couple_from_unique_faces = r_couple_from_unique_faces;
    audit.reconstruction_error = r_couple_from_unique_faces - result.stage5_r_couple;
    const double tolerance = 8192.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::fabs(result.stage5_r_couple));
    audit.ledger_match = audit.complete_unique_coverage &&
        audit.row_count == audit.expected_row_count &&
        std::isfinite(audit.reconstruction_error) &&
        std::fabs(audit.reconstruction_error) <= tolerance &&
        std::isfinite(audit.global_terms[FIXED_PAIR_RECONSTRUCTION_ERROR]
                      .signed_sum) &&
        std::fabs(audit.global_terms[FIXED_PAIR_RECONSTRUCTION_ERROR]
                  .signed_sum) <= tolerance;
    out << "# face_row_count " << audit.row_count
        << " expected_row_count " << audit.expected_row_count
        << " unique_face_count " << audit.unique_face_count
        << " complete_unique_coverage " << audit.complete_unique_coverage
        << " ledger_match " << audit.ledger_match << "\n";
    double values[2] = {audit.r_couple_from_unique_faces,
                        audit.reconstruction_error};
    int flags[5] = {audit.row_count, audit.expected_row_count,
                    audit.unique_face_count,
                    audit.complete_unique_coverage ? 1 : 0,
                    audit.ledger_match ? 1 : 0};
    MPI_Bcast(values, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(flags, 5, MPI_INT, 0, MPI_COMM_WORLD);
    return audit;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    RuntimeConfig config = load_runtime_config();
    RuntimeOptions runtime = parse_runtime_options(argc, argv, mpi_rank, mpi_size);
    std::string runtime_error;
    if (!runtime.velocity_grid_remap_enabled &&
        !prepare_output_directory(runtime, mpi_rank, mpi_size, runtime_error)) {
        if (mpi_rank == 0) std::fprintf(stderr, "Runtime setup error: %s\n", runtime_error.c_str());
        MPI_Finalize();
        return 2;
    }

    if (mpi_rank == 0) {
        printf("============================================================\n");
        printf("  Background-electron VFP + fixed ions + PIC beam solver\n");
        printf("  Cylindrical electron momentum grid: (u_parallel, u_perp)\n");
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
        printf("Electron momentum grid: Nu_parallel x Nu_perp = %d x %d (core=%d, tail/side=%d), conservative cylindrical mass M\n",
               Param::Nv, Param::Nmu, Param::Nv_core, Param::Nv_tail);
        printf("Electron momentum domain: core |u_parallel| <= %.3f, extended |u_parallel| <= %.3f, 0 <= u_perp <= %.3f\n",
               Param::momentum_upar_core_max, Param::momentum_upar_extended_max,
               Param::momentum_umax);
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
    BeamPIC beam;
    beam.init(sgrid);

    EMFields fields;
    fields.init(sgrid);

    if (mpi_rank == 0 && !runtime.velocity_grid_remap_enabled) {
        const CylindricalVelocityGrid::NestedGridAudit& grid_audit =
            bkg_e.cgrid.nested_audit;
        std::ofstream out(output_path(runtime, "velocity_grid_audit.result").c_str(),
                          std::ios::out | std::ios::trunc);
        out << std::scientific << std::setprecision(17)
            << "Nv_core " << Param::Nv_core << "\n"
            << "Nv_tail_per_side " << Param::Nv_tail << "\n"
            << "Nv_total " << Param::Nv << "\n"
            << "upar_core_max " << Param::momentum_upar_core_max << "\n"
            << "upar_extended_max " << Param::momentum_upar_extended_max << "\n"
            << "core_face_identity_linf "
            << grid_audit.core_face_identity_linf << "\n"
            << "core_cell_identity_linf "
            << grid_audit.core_cell_identity_linf << "\n"
            << "core_vx_identity_linf " << grid_audit.core_vx_identity_linf << "\n"
            << "core_kinetic_energy_identity_linf "
            << grid_audit.core_kinetic_energy_identity_linf << "\n"
            << "symmetry_linf " << grid_audit.symmetry_linf << "\n"
            << "max_adjacent_width_ratio "
            << grid_audit.max_adjacent_width_ratio << "\n"
            << "phase_volume_relative_error "
            << grid_audit.phase_volume_relative_error << "\n";
    }

    VlasovAmpereMidpointSolver midpoint_solver;
    midpoint_solver.set_step_diagnostics_enabled(false);
    midpoint_solver.set_accepted_energy_audit_enabled(false);
    midpoint_solver.set_low_order_only(false);
    midpoint_solver.set_nonuniform_high_order_enabled(true);
    midpoint_solver.set_fct_enabled(true);
    midpoint_solver.set_background_coupling_mode(
        runtime.background_coupling_mode == 1
        ? VlasovAmpereMidpointSolver::DUAL_U_COUPLING
        : VlasovAmpereMidpointSolver::LEGACY_COUPLING);
    midpoint_solver.set_face_pairing_mode(
        runtime.face_pairing_mode == 1
        ? VlasovAmpereMidpointSolver::FACE_PAIRING_REGULARIZED
        : VlasovAmpereMidpointSolver::FACE_PAIRING_CELL_BASELINE);
    midpoint_solver.set_face_pairing_parameters(
        runtime.face_pairing_sigma_cutoff, runtime.face_pairing_lambda,
        runtime.face_pairing_eta, runtime.face_pairing_trust_fraction);
    midpoint_solver.set_face_pairing_acceptance_parameters(
        runtime.face_pairing_correction_trust_fraction,
        runtime.face_pairing_energy_pair_tolerance,
        runtime.face_pairing_energy_residual_fraction,
        runtime.face_pairing_mass_relative_tolerance,
        runtime.face_pairing_f_residual_growth_tolerance);
    midpoint_solver.set_max_midpoint_iterations(runtime.midpoint_max_iters);
    midpoint_solver.set_midpoint_initial_guess_mode(
        runtime.midpoint_initial_guess_mode ==
                RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR
        ? VlasovAmpereMidpointSolver::MIDPOINT_INITIAL_GUESS_FIELD_LINEAR
        : VlasovAmpereMidpointSolver::MIDPOINT_INITIAL_GUESS_NONE);
    midpoint_solver.set_midpoint_acceleration_mode(
        runtime.midpoint_acceleration_mode == RUNTIME_MIDPOINT_ACCELERATION_AITKEN
        ? VlasovAmpereMidpointSolver::MIDPOINT_ACCELERATION_AITKEN
        : runtime.midpoint_acceleration_mode == RUNTIME_MIDPOINT_ACCELERATION_ANDERSON
        ? VlasovAmpereMidpointSolver::MIDPOINT_ACCELERATION_ANDERSON
        : VlasovAmpereMidpointSolver::MIDPOINT_ACCELERATION_NONE);
    midpoint_solver.set_anderson_depth(runtime.anderson_depth);
    midpoint_solver.set_acceleration_start_iter(runtime.acceleration_start_iter);
    midpoint_solver.set_acceleration_accept_ratio(runtime.acceleration_accept_ratio);
    midpoint_solver.set_acceleration_max_coefficient(runtime.acceleration_max_coefficient);
    midpoint_solver.set_capture_midpoint_input(runtime.diagnostic_level >= 2);
    midpoint_solver.set_midpoint_iteration_trace(runtime.diagnostic_level >= 2);
    midpoint_solver.set_progress_trace_window_fs(
        runtime.midpoint_trace_start_fs, runtime.midpoint_trace_end_fs);

    if (runtime.velocity_grid_remap_enabled) {
        CheckpointControlState control = {};
        CheckpointVelocityRemapAudit remap_audit = {};
        const bool read_ok = read_checkpoint(
            runtime.remap_checkpoint_source, control, bkg_e, beam, fields,
            sgrid, mpi_rank, mpi_size, runtime_error,
            midpoint_solver.low_order_only(),
            midpoint_solver.nonuniform_high_order_enabled(),
            midpoint_solver.fct_enabled(), true, &remap_audit);
        if (!read_ok) {
            if (mpi_rank == 0) {
                std::fprintf(stderr, "Velocity-grid remap input error: %s\n",
                             runtime_error.c_str());
            }
            MPI_Finalize();
            return 3;
        }
        const bool write_ok = write_checkpoint(
            runtime.remap_checkpoint_destination, control, bkg_e, beam, fields,
            sgrid, mpi_rank, mpi_size, runtime_error,
            midpoint_solver.low_order_only(),
            midpoint_solver.nonuniform_high_order_enabled(),
            midpoint_solver.fct_enabled(), &remap_audit);
        const bool audit_ok = write_ok && write_checkpoint_velocity_remap_audit(
            runtime.remap_checkpoint_destination, remap_audit,
            mpi_rank, mpi_size, runtime_error);
        if (!audit_ok) {
            if (mpi_rank == 0) {
                std::fprintf(stderr, "Velocity-grid remap output error: %s\n",
                             runtime_error.c_str());
            }
            MPI_Finalize();
            return 4;
        }
        if (mpi_rank == 0) {
            std::printf("Velocity-grid checkpoint remap completed: %s -> %s\n",
                        runtime.remap_checkpoint_source.c_str(),
                        runtime.remap_checkpoint_destination.c_str());
        }
        MPI_Finalize();
        return 0;
    }

    if (mpi_rank == 0) {
        const char* const beam_ledger_mode = runtime.beam_ledger_mode == BEAM_LEDGER_FULL
            ? "full" : runtime.beam_ledger_mode == BEAM_LEDGER_SUMMARY ? "summary" : "off";
        const char* const midpoint_acceleration =
            runtime.midpoint_acceleration_mode == RUNTIME_MIDPOINT_ACCELERATION_AITKEN
            ? "aitken" : runtime.midpoint_acceleration_mode ==
                RUNTIME_MIDPOINT_ACCELERATION_ANDERSON ? "anderson" : "none";
        const char* const midpoint_initial_guess =
            runtime.midpoint_initial_guess_mode ==
                    RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR
            ? "field-linear" : "none";
        printf("Transport configuration: low_order_only=%s, "
               "nonuniform_high_order=%s, FCT=%s, background_coupling=%s, "
               "midpoint_max_iters=%d, midpoint_initial_guess=%s, "
               "midpoint_acceleration=%s, "
               "anderson_depth=%d, acceleration_start_iter=%d, "
               "acceleration_accept_ratio=%.3f, acceleration_max_coefficient=%.3f, "
               "beam_ledger_mode=%s, face_pairing=%s, "
               "face_sigma_cutoff=%.3e, face_lambda=%.3e, "
               "face_eta=%.3e, face_capacity=%.3f, "
               "face_correction_trust=%.3f, face_energy_pair_tol=%.3e, "
               "face_energy_scale_fraction=%.3f, face_mass_tol=%.3e, "
               "face_f_growth_tol=%.3f, accepted_energy_audit_cadence=%d, "
               "accepted_energy_audit_window_fs=[%.6g,%.6g]\n",
               midpoint_solver.low_order_only() ? "ON" : "OFF",
               midpoint_solver.nonuniform_high_order_enabled() ? "ON" : "OFF",
               midpoint_solver.fct_enabled() ? "ON" : "OFF",
               midpoint_solver.background_coupling_mode() ==
                   VlasovAmpereMidpointSolver::DUAL_U_COUPLING
                   ? "dual_u" : "legacy",
                midpoint_solver.max_midpoint_iterations(),
                midpoint_initial_guess, midpoint_acceleration,
                runtime.anderson_depth, runtime.acceleration_start_iter,
                runtime.acceleration_accept_ratio, runtime.acceleration_max_coefficient,
                beam_ledger_mode,
                midpoint_solver.face_pairing_mode() ==
                    VlasovAmpereMidpointSolver::FACE_PAIRING_REGULARIZED
                    ? "regularized" : "cell-baseline",
                runtime.face_pairing_sigma_cutoff,
                runtime.face_pairing_lambda, runtime.face_pairing_eta,
                runtime.face_pairing_trust_fraction,
                runtime.face_pairing_correction_trust_fraction,
                runtime.face_pairing_energy_pair_tolerance,
                runtime.face_pairing_energy_residual_fraction,
                runtime.face_pairing_mass_relative_tolerance,
                runtime.face_pairing_f_residual_growth_tolerance,
                runtime.accepted_energy_audit_cadence,
                runtime.accepted_energy_audit_start_fs,
                runtime.accepted_energy_audit_end_fs);
        printf("Soft-candidate transaction: field_tol=%.3e, current_tol=%.3e, "
               "energy_p99_reference=%.3e, energy_absolute_limit=%.3e, "
               "predictor_retry_max=1, physical_split_depth_max=1\n",
               runtime.soft_candidate_field_tolerance,
               runtime.soft_candidate_current_tolerance,
               runtime.soft_candidate_energy_p99_reference,
               runtime.soft_candidate_energy_absolute_limit);
        std::ofstream pairing_configuration(
            output_path(runtime, "face_pairing_configuration.result").c_str());
        pairing_configuration << std::scientific << std::setprecision(16)
            << "mode "
            << (midpoint_solver.face_pairing_mode() ==
                    VlasovAmpereMidpointSolver::FACE_PAIRING_REGULARIZED
                    ? "regularized" : "cell-baseline") << "\n"
            << "sigma_cutoff " << runtime.face_pairing_sigma_cutoff << "\n"
            << "lambda " << runtime.face_pairing_lambda << "\n"
            << "eta " << runtime.face_pairing_eta << "\n"
            << "capacity_fraction "
            << runtime.face_pairing_trust_fraction << "\n"
            << "correction_trust_fraction "
            << runtime.face_pairing_correction_trust_fraction << "\n"
            << "energy_pair_tolerance "
            << runtime.face_pairing_energy_pair_tolerance << "\n"
            << "energy_residual_fraction "
            << runtime.face_pairing_energy_residual_fraction << "\n"
            << "mass_relative_tolerance "
            << runtime.face_pairing_mass_relative_tolerance << "\n"
            << "f_residual_growth_tolerance "
            << runtime.face_pairing_f_residual_growth_tolerance << "\n";
        std::ofstream accepted_energy_audit_configuration(
            output_path(runtime, "accepted_energy_audit_configuration.result").c_str());
        accepted_energy_audit_configuration << std::scientific << std::setprecision(16)
            << "cadence " << runtime.accepted_energy_audit_cadence << "\n"
            << "window_start_fs " << runtime.accepted_energy_audit_start_fs << "\n"
            << "window_end_fs " << runtime.accepted_energy_audit_end_fs << "\n";
    }

    if (runtime.operator_audit_mode) {
        VlasovAmpereMidpointSolver::MidpointAuditState audit_state;
        if (!read_midpoint_audit_state(runtime.operator_audit_dir, audit_state,
                                       bkg_e, fields, sgrid, mpi_rank, mpi_size,
                                       runtime_error)) {
            if (mpi_rank == 0) {
                std::fprintf(stderr, "Operator-audit input error: %s\n",
                             runtime_error.c_str());
            }
            MPI_Finalize();
            return 3;
        }
        midpoint_solver.set_low_order_only(audit_state.low_order_only);
        midpoint_solver.set_nonuniform_high_order_enabled(
            audit_state.high_order_enabled);
        midpoint_solver.set_fct_enabled(audit_state.fct_enabled);
        midpoint_solver.set_background_coupling_mode(
            audit_state.background_coupling_mode ==
            VlasovAmpereMidpointSolver::DUAL_U_COUPLING
            ? VlasovAmpereMidpointSolver::DUAL_U_COUPLING
            : VlasovAmpereMidpointSolver::LEGACY_COUPLING);
        if (mpi_rank == 0) {
            std::printf("Operator-audit background coupling: %s\n",
                        midpoint_solver.background_coupling_mode() ==
                        VlasovAmpereMidpointSolver::DUAL_U_COUPLING
                        ? "dual_u" : "legacy");
        }
        const VlasovAmpereMidpointSolver::MidpointOperatorEvaluation evaluated =
            midpoint_solver.evaluate_fixed_midpoint_operator(
                audit_state.bkg_n, beam, audit_state.fields_n,
                audit_state.operator_input_guess, audit_state.fields_end_guess,
                audit_state.j_beam_face_mid, audit_state.coupling_layout,
                sgrid, audit_state.dt_s,
                audit_state.time_s, mpi_rank, mpi_size);
        const auto max_difference = [](const std::vector<double>& a,
                                       const std::vector<double>& b) {
            if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
            double value = 0.0;
            for (size_t i = 0; i < a.size(); ++i)
                value = std::max(value, std::fabs(a[i] - b[i]));
            return value;
        };
        const auto max_magnitude = [](const std::vector<double>& values) {
            double value = 0.0;
            for (size_t i = 0; i < values.size(); ++i)
                value = std::max(value, std::fabs(values[i]));
            return value;
        };
        double differences[6] = {
            max_difference(evaluated.j_bkg_face_mid, audit_state.reference_jn_face),
            max_difference(evaluated.j_bkg_energy_cell_mid, audit_state.reference_je_cell),
            max_difference(evaluated.j_bkg_energy_debug_face,
                           audit_state.reference_gstar_je_face),
            max_difference(evaluated.j_beam_face_mid, audit_state.j_beam_face_mid),
            std::fabs(evaluated.stage5_r_fv - audit_state.reference_stage5_r_fv),
            std::fabs(evaluated.stage5_r_couple - audit_state.reference_stage5_r_couple)
        };
        MPI_Allreduce(MPI_IN_PLACE, differences, 6, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        double reference_scales[6] = {
            max_magnitude(audit_state.reference_jn_face),
            max_magnitude(audit_state.reference_je_cell),
            max_magnitude(audit_state.reference_gstar_je_face),
            max_magnitude(audit_state.j_beam_face_mid),
            std::fabs(audit_state.reference_stage5_r_fv),
            std::fabs(audit_state.reference_stage5_r_couple)
        };
        MPI_Allreduce(MPI_IN_PLACE, reference_scales, 6, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        double relative_differences[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        bool reference_match = true;
        for (int i = 0; i < 6; ++i) {
            relative_differences[i] = differences[i] /
                std::max(1.0, reference_scales[i]);
            reference_match = reference_match &&
                relative_differences[i] <= 4096.0 * std::numeric_limits<double>::epsilon();
        }
        double periodic_seam_difference = 0.0;
        double periodic_seam_scale = 0.0;
        for (size_t i = 0; i < evaluated.periodic_seam_face_audit.size(); ++i) {
            periodic_seam_difference = std::max(periodic_seam_difference,
                std::fabs(evaluated.periodic_seam_face_audit[i] -
                          audit_state.periodic_seam_face_audit[i]));
            periodic_seam_scale = std::max(periodic_seam_scale,
                std::max(std::fabs(evaluated.periodic_seam_face_audit[i]),
                         std::fabs(audit_state.periodic_seam_face_audit[i])));
        }
        const double periodic_seam_relative_difference = periodic_seam_difference /
            std::max(1.0, periodic_seam_scale);
        const bool periodic_seam_match = std::isfinite(periodic_seam_difference) &&
            periodic_seam_relative_difference <=
                4096.0 * std::numeric_limits<double>::epsilon();
        const bool layout_match =
            evaluated.coupling_beam_front_ix == audit_state.coupling_layout.beam_front_ix &&
            std::fabs(evaluated.coupling_wave_core_end_m -
                      audit_state.coupling_layout.wave_core_end_m) <=
                64.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, std::fabs(audit_state.coupling_layout.wave_core_end_m));
        const FixedMidpointFacePairingAudit face_pairing_audit =
            write_fixed_midpoint_face_pairing(
                output_path(runtime, "fixed_midpoint_face_pairing.dat"),
                audit_state.step, audit_state.time_s, audit_state.dt_s, evaluated,
                audit_state.fields_n, audit_state.fields_end_guess,
                sgrid, mpi_rank, mpi_size, false);
        const bool audit_ok = evaluated.state_advanced && !evaluated.failed &&
            std::isfinite(differences[0]) && std::isfinite(differences[1]) &&
            std::isfinite(differences[2]) && std::isfinite(differences[3]) &&
            std::isfinite(differences[4]) && std::isfinite(differences[5]) &&
            reference_match && layout_match && periodic_seam_match &&
            face_pairing_audit.ledger_match;
        if (mpi_rank == 0) {
            std::ofstream out(output_path(runtime, "fixed_midpoint_operator_audit.result").c_str());
            out << std::setprecision(17)
                << "audit_input " << runtime.operator_audit_dir << "\n"
                << "background_coupling_mode "
                << (audit_state.background_coupling_mode ==
                    VlasovAmpereMidpointSolver::DUAL_U_COUPLING
                    ? "dual_u" : "legacy") << "\n"
                << "step " << audit_state.step << "\n"
                << "time_s " << audit_state.time_s << "\n"
                << "dt_s " << audit_state.dt_s << "\n"
                << "state_advanced " << evaluated.state_advanced << "\n"
                << "failed " << evaluated.failed << "\n"
                << "j_n_max_abs_difference " << differences[0] << "\n"
                << "j_e_max_abs_difference " << differences[1] << "\n"
                << "gstar_j_e_max_abs_difference " << differences[2] << "\n"
                << "beam_current_max_abs_difference " << differences[3] << "\n"
                << "max_abs_recomputed_JN " << max_magnitude(evaluated.j_bkg_face_mid) << "\n"
                << "max_abs_recomputed_JE " << max_magnitude(evaluated.j_bkg_energy_cell_mid) << "\n"
                << "max_abs_recomputed_GstarJE " << max_magnitude(evaluated.j_bkg_energy_debug_face) << "\n"
                << "j_n_relative_difference " << relative_differences[0] << "\n"
                << "j_e_relative_difference " << relative_differences[1] << "\n"
                << "gstar_j_e_relative_difference " << relative_differences[2] << "\n"
                << "beam_current_relative_difference " << relative_differences[3] << "\n"
                << "stage5_R_FV_reference " << audit_state.reference_stage5_r_fv << "\n"
                << "stage5_R_FV_max_abs_difference " << differences[4] << "\n"
                << "stage5_R_FV_relative_difference " << relative_differences[4] << "\n"
                << "stage5_R_couple_reference " << audit_state.reference_stage5_r_couple << "\n"
                << "stage5_R_couple_max_abs_difference " << differences[5] << "\n"
                << "stage5_R_couple_relative_difference " << relative_differences[5] << "\n"
                << "face_pairing_row_count " << face_pairing_audit.row_count << "\n"
                << "face_pairing_expected_row_count "
                << face_pairing_audit.expected_row_count << "\n"
                << "face_pairing_unique_face_count "
                << face_pairing_audit.unique_face_count << "\n"
                << "face_pairing_complete_unique_coverage "
                << face_pairing_audit.complete_unique_coverage << "\n"
                << "face_pairing_R_couple_from_unique_faces "
                << face_pairing_audit.r_couple_from_unique_faces << "\n"
                << "face_pairing_R_couple_reconstruction_error "
                << face_pairing_audit.reconstruction_error << "\n"
                << "face_pairing_ledger_match "
                << face_pairing_audit.ledger_match << "\n"
                << "unique_face_count " << face_pairing_audit.unique_face_count << "\n"
                << "sum_of_face_rows "
                << face_pairing_audit.r_couple_from_unique_faces << "\n"
                << "global_reconstruction_error "
                << face_pairing_audit.reconstruction_error << "\n";
            const char* const pair_term_names[FIXED_PAIR_TERM_COUNT] = {
                "R_pair_low", "delta_R_x_high", "delta_R_u_center",
                "delta_R_u_upwind", "delta_R_x_FCT", "delta_R_u_FCT",
                "R_pair_reconstructed_final", "R_pair_final",
                "R_pair_reconstruction_error"};
            const double pair_final_signed =
                face_pairing_audit.global_terms[FIXED_PAIR_FINAL].signed_sum;
            const double pair_final_signed_scale =
                std::fabs(pair_final_signed) > std::numeric_limits<double>::min()
                ? pair_final_signed
                : std::numeric_limits<double>::min();
            const double pair_final_absolute_scale = std::max(
                std::fabs(pair_final_signed), std::numeric_limits<double>::min());
            for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
                const FixedPairTermSummary& summary =
                    face_pairing_audit.global_terms[static_cast<size_t>(term)];
                out << pair_term_names[term] << " " << summary.signed_sum << "\n"
                    << pair_term_names[term] << "_absolute_sum "
                    << summary.absolute_sum << "\n"
                    << pair_term_names[term] << "_Linf " << summary.linf << "\n"
                    << pair_term_names[term] << "_max_global_face "
                    << summary.max_global_face << "\n"
                    << pair_term_names[term] << "_fraction_of_final "
                    << summary.signed_sum / pair_final_signed_scale << "\n"
                    << pair_term_names[term] << "_absolute_fraction_of_final "
                    << summary.absolute_sum / pair_final_absolute_scale << "\n";
            }
            for (int region = 0; region < 6; ++region) {
                for (int term = 0; term < FIXED_PAIR_TERM_COUNT; ++term) {
                    const FixedPairTermSummary& summary =
                        face_pairing_audit.region_terms[static_cast<size_t>(region)]
                            [static_cast<size_t>(term)];
                    out << "pair_region " << region << " " << pair_term_names[term]
                        << " signed_sum " << summary.signed_sum
                        << " absolute_sum " << summary.absolute_sum
                        << " Linf " << summary.linf
                        << " max_global_face " << summary.max_global_face << "\n";
                }
            }
            out
                << "periodic_seam_JN_face0_reference "
                << audit_state.periodic_seam_face_audit[0] << "\n"
                << "periodic_seam_JN_face0_recomputed "
                << evaluated.periodic_seam_face_audit[0] << "\n"
                << "periodic_seam_JN_faceL_reference "
                << audit_state.periodic_seam_face_audit[1] << "\n"
                << "periodic_seam_JN_faceL_recomputed "
                << evaluated.periodic_seam_face_audit[1] << "\n"
                << "periodic_seam_GstarJE_face0_reference "
                << audit_state.periodic_seam_face_audit[2] << "\n"
                << "periodic_seam_GstarJE_face0_recomputed "
                << evaluated.periodic_seam_face_audit[2] << "\n"
                << "periodic_seam_GstarJE_faceL_reference "
                << audit_state.periodic_seam_face_audit[3] << "\n"
                << "periodic_seam_GstarJE_faceL_recomputed "
                << evaluated.periodic_seam_face_audit[3] << "\n"
                << "periodic_seam_JN_minus_GstarJE_face0_reference "
                << audit_state.periodic_seam_face_audit[4] << "\n"
                << "periodic_seam_JN_minus_GstarJE_face0_recomputed "
                << evaluated.periodic_seam_face_audit[4] << "\n"
                << "periodic_seam_JN_minus_GstarJE_faceL_reference "
                << audit_state.periodic_seam_face_audit[5] << "\n"
                << "periodic_seam_JN_minus_GstarJE_faceL_recomputed "
                << evaluated.periodic_seam_face_audit[5] << "\n"
                << "periodic_seam_JN_alias_error_reference "
                << (audit_state.periodic_seam_face_audit[1] -
                    audit_state.periodic_seam_face_audit[0]) << "\n"
                << "periodic_seam_JN_alias_error_recomputed "
                << (evaluated.periodic_seam_face_audit[1] -
                    evaluated.periodic_seam_face_audit[0]) << "\n"
                << "periodic_seam_GstarJE_alias_error_reference "
                << (audit_state.periodic_seam_face_audit[3] -
                    audit_state.periodic_seam_face_audit[2]) << "\n"
                << "periodic_seam_GstarJE_alias_error_recomputed "
                << (evaluated.periodic_seam_face_audit[3] -
                    evaluated.periodic_seam_face_audit[2]) << "\n"
                << "seam_JN_face0 " << evaluated.periodic_seam_face_audit[0] << "\n"
                << "seam_JN_aliasN " << evaluated.periodic_seam_face_audit[1] << "\n"
                << "seam_GstarJE_face0 " << evaluated.periodic_seam_face_audit[2] << "\n"
                << "seam_GstarJE_aliasN " << evaluated.periodic_seam_face_audit[3] << "\n"
                << "seam_deltaJ_face0 " << evaluated.periodic_seam_face_audit[4] << "\n"
                << "seam_deltaJ_aliasN " << evaluated.periodic_seam_face_audit[5] << "\n"
                << "seam_JN_alias_error "
                << (evaluated.periodic_seam_face_audit[1] -
                    evaluated.periodic_seam_face_audit[0]) << "\n"
                << "seam_GstarJE_alias_error "
                << (evaluated.periodic_seam_face_audit[3] -
                    evaluated.periodic_seam_face_audit[2]) << "\n"
                << "alias_error "
                << std::max(std::fabs(evaluated.periodic_seam_face_audit[1] -
                                      evaluated.periodic_seam_face_audit[0]),
                            std::fabs(evaluated.periodic_seam_face_audit[3] -
                                      evaluated.periodic_seam_face_audit[2])) << "\n"
                << "periodic_seam_max_abs_difference "
                << periodic_seam_difference << "\n"
                << "periodic_seam_relative_difference "
                << periodic_seam_relative_difference << "\n"
                << "periodic_seam_match " << (periodic_seam_match ? 1 : 0) << "\n"
                << "coupling_beam_front_ix_reference "
                << audit_state.coupling_layout.beam_front_ix << "\n"
                << "coupling_beam_front_ix_recomputed "
                << evaluated.coupling_beam_front_ix << "\n"
                << "coupling_beam_front_ix_replayed_from_B_dump 1\n"
                << "coupling_wave_core_end_m_reference "
                << audit_state.coupling_layout.wave_core_end_m << "\n"
                << "coupling_wave_core_end_m_recomputed "
                << evaluated.coupling_wave_core_end_m << "\n"
                << "coupling_layout_match " << (layout_match ? 1 : 0) << "\n"
                << "final_x_flux_count " << evaluated.final_x_flux.size() << "\n"
                << "final_u_flux_count " << evaluated.final_u_flux.size() << "\n"
                << "stage5_R_FV " << evaluated.stage5_r_fv << "\n"
                << "stage5_R_couple " << evaluated.stage5_r_couple << "\n"
                << "stage5_R_couple_centered "
                << evaluated.stage5_r_couple_centered << "\n"
                << "stage5_R_couple_upwind_stabilization "
                << evaluated.stage5_r_couple_upwind_stabilization << "\n"
                << "stage5_R_couple_fct_stabilization "
                << evaluated.stage5_r_couple_fct_stabilization << "\n"
                << "R_couple_decomposition_error "
                << (evaluated.stage5_r_couple -
                    (evaluated.stage5_r_couple_centered +
                     evaluated.stage5_r_couple_upwind_stabilization +
                     evaluated.stage5_r_couple_fct_stabilization)) << "\n"
                << "coupling_rJ_sum " << evaluated.coupling_rj_global_sum << "\n"
                << "coupling_rK_sum " << evaluated.coupling_rk_global_sum << "\n"
                << "coupling_face_work_sum "
                << evaluated.coupling_face_work_jn_global_sum << "\n"
                << "coupling_rJ_identity_error "
                << evaluated.coupling_rj_reconstruction_error << "\n"
                << "coupling_rK_identity_error "
                << evaluated.coupling_rk_reconstruction_error << "\n"
                << "coupling_face_work_identity_error "
                << evaluated.coupling_face_work_jn_reconstruction_error << "\n"
                << "coupling_R_couple_centered_sum "
                << evaluated.coupling_rj_centered_global_sum << "\n"
                << "coupling_R_couple_centered_identity_error "
                << evaluated.coupling_rj_centered_reconstruction_error << "\n"
                << "coupling_R_couple_upwind_sum "
                << evaluated.coupling_rj_upwind_global_sum << "\n"
                << "coupling_R_couple_upwind_identity_error "
                << evaluated.coupling_rj_upwind_reconstruction_error << "\n"
                << "coupling_R_couple_fct_sum "
                << evaluated.coupling_rj_fct_global_sum << "\n"
                << "coupling_R_couple_fct_identity_error "
                << evaluated.coupling_rj_fct_reconstruction_error << "\n"
                << "audit_finite " << (audit_ok ? 1 : 0) << "\n";
            for (size_t region = 0; region < evaluated.coupling_regions.size(); ++region) {
                const VlasovAmpereMidpointSolver::CouplingRegionDiagnostics& r =
                    evaluated.coupling_regions[region];
                out << "region " << region
                    << " sum_rJ " << r.sum_rj
                    << " sum_abs_rJ " << r.sum_abs_rj
                    << " Linf_rJ " << r.linf_rj
                    << " sum_rK " << r.sum_rk
                    << " sum_abs_rK " << r.sum_abs_rk
                    << " Linf_rK " << r.linf_rk
                    << " face_work_JN " << r.face_work_jn
                    << " delta_K_bkg " << r.delta_ke_bkg
                    << " fct_work " << r.fct_work
                    << " R_couple_centered " << r.r_couple_centered
                    << " R_couple_upwind_stabilization "
                    << r.r_couple_upwind_stabilization
                    << " R_couple_fct_stabilization "
                    << r.r_couple_fct_stabilization
                    << "\n";
            }
            std::printf("Fixed midpoint audit: JN=%e JE=%e GstarJE=%e Beam=%e\n",
                        differences[0], differences[1], differences[2], differences[3]);
        }
        MPI_Finalize();
        return audit_ok ? 0 : 4;
    }
    CollisionOperator collision;
    Diagnostics diag;
    diag.init(runtime.output_dir, mpi_rank, config.enable_debug_diagnostics,
              config.enable_step_diagnostics || runtime.diagnostic_level >= 2,
              runtime.accepted_energy_audit_cadence > 0);

    double dt = compute_dt(bkg_e, sgrid);
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt *= runtime.dt_scale;
    int nsteps = static_cast<int>(std::ceil(runtime.stop_time_fs * Const::femto / dt));

    if (mpi_rank == 0) {
        printf("Time step: dt = %.4e s (%.4f fs)\n", dt, dt / Const::femto);
        printf("Total steps: %d\n", nsteps);
        printf("============================================================\n");
    }

    double cumulative_collision_energy_delta = 0.0;
    double current_time = 0.0;
    long long restored_step = 0;
    double next_snapshot = Param::dt_snapshot;
    int last_snapshot_step = 0;
    if (runtime.restart_enabled) {
        CheckpointControlState checkpoint_control;
        if (!read_checkpoint(runtime.restart_dir, checkpoint_control, bkg_e, beam,
                             fields, sgrid, mpi_rank, mpi_size, runtime_error,
                             midpoint_solver.low_order_only(),
                             midpoint_solver.nonuniform_high_order_enabled(),
                             midpoint_solver.fct_enabled())) {
            if (mpi_rank == 0) std::fprintf(stderr, "Restart error: %s\n", runtime_error.c_str());
            MPI_Abort(MPI_COMM_WORLD, 11);
            return 11;
        }
        current_time = checkpoint_control.time_s;
        restored_step = checkpoint_control.step;
        dt = checkpoint_control.dt_s * runtime.dt_scale;
        next_snapshot = checkpoint_control.next_snapshot_s;
        last_snapshot_step = checkpoint_control.last_snapshot_step;
        cumulative_collision_energy_delta = checkpoint_control.cumulative_collision_energy_delta;
        midpoint_solver.synchronize_background_ghosts(bkg_e, sgrid, mpi_rank, mpi_size);
        bkg_e.compute_moments();
        beam.deposit_density(sgrid, mpi_rank, mpi_size);
        fields.sync_cell_ex_from_faces(mpi_rank, mpi_size);
        fields.set_charge_density(bkg_e, beam.density, ion_density_profile);
        const CheckpointStateHashes restored_hashes = checkpoint_state_hashes(
            bkg_e, beam, fields, sgrid, mpi_rank, mpi_size);
        CheckpointStateHashes reference_hashes = {0ULL, 0ULL, 0ULL};
        const bool have_reference_hashes = read_checkpoint_reference_hashes(
            runtime.restart_dir, reference_hashes, mpi_rank, mpi_size);
        double restart_n = 0.0, restart_ke = 0.0;
        bkg_e.total_particle_number_and_energy(restart_n, restart_ke);
        double restart_values[4] = {restart_n, restart_ke,
                                    beam.total_particle_number(sgrid),
                                    beam.total_kinetic_energy()};
        MPI_Allreduce(MPI_IN_PLACE, restart_values, 4, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        double restart_ex = 0.0;
        for (size_t iface = 0; iface < fields.Ex_face.size(); ++iface)
            restart_ex = std::max(restart_ex, std::fabs(fields.Ex_face[iface]));
        MPI_Allreduce(MPI_IN_PLACE, &restart_ex, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            std::ofstream restart_audit(output_path(runtime, "restart_audit.dat").c_str());
            restart_audit << "# step time_s dt_s N_bkg KE_bkg beam_weight KE_beam field_energy max_abs_Ex "
                          << "hash_bkg hash_field hash_beam ref_hash_bkg ref_hash_field ref_hash_beam hashes_match\n";
            restart_audit << std::setprecision(17) << restored_step << " "
                          << current_time << " " << dt << " " << restart_values[0]
                          << " " << restart_values[1] << " " << restart_values[2]
                          << " " << restart_values[3] << " " << fields.total_energy()
                          << " " << restart_ex << " " << restored_hashes.background
                          << " " << restored_hashes.field_faces << " " << restored_hashes.beam
                          << " " << reference_hashes.background << " "
                          << reference_hashes.field_faces << " " << reference_hashes.beam
                          << " " << (have_reference_hashes &&
                              restored_hashes.background == reference_hashes.background &&
                              restored_hashes.field_faces == reference_hashes.field_faces &&
                              restored_hashes.beam == reference_hashes.beam ? 1 : 0) << "\n";
        }
        if (mpi_rank == 0) std::printf("Restarted checkpoint: step=%lld t=%.6f fs\n", restored_step, current_time / Const::femto);
    } else {
        bkg_e.initialize_maxwellian();
        beam.deposit_density(sgrid, mpi_rank, mpi_size);
    }
    bool moments_current = false;
    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    if (!runtime.restart_enabled) fields.solve_poisson(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(0, 0.0, "initial", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(current_time, static_cast<int>(restored_step), bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    std::vector<double> latest_bkg_energy_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    std::vector<double> latest_bkg_ampere_current_face(
        static_cast<size_t>(sgrid.nx_local + 1), 0.0);
    bool latest_bkg_energy_current_valid = false;
    write_snapshot(diag, current_time, bkg_e, beam, fields, ion_density_profile,
                   sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                   &latest_bkg_energy_current_face,
                   &latest_bkg_ampere_current_face,
                   latest_bkg_energy_current_valid);

    int stdout_freq = 1000;
    double cumulative_bkg_energy_residual = 0.0;
    double cumulative_total_energy_residual_for_audit = 0.0;
    std::ofstream accepted_coupling_region_monitor;
    std::ofstream fct_macro_budget_monitor;
    std::ofstream high_order_candidate_monitor;
    std::ofstream final_dual_u_pairing_monitor;
    std::ofstream midpoint_iteration_residual_monitor;
    std::ofstream trial_retry_monitor;
    std::ofstream physical_interval_acceptance_monitor;
    std::ofstream beam_half_step_ledger;
    std::array<double, 6> beam_ledger_local_totals = {{0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0}};
    std::array<double, 6> beam_ledger_full_totals = {{0.0, 0.0, 0.0,
                                                       0.0, 0.0, 0.0}};
    std::array<double, 12> beam_reference_rows = {{0.0}};
    int beam_reference_row_count = 0;
    const double beam_ledger_start_time_s = current_time;
    const bool beam_ledger_enabled = runtime.beam_ledger_mode != BEAM_LEDGER_OFF;
    const bool write_beam_ledger = runtime.beam_ledger_mode == BEAM_LEDGER_FULL;
    const auto flush_full_beam_ledger = [&]() {
        if (mpi_rank == 0 && beam_half_step_ledger.is_open()) {
            beam_half_step_ledger.flush();
        }
    };
    if (mpi_rank == 0) {
        accepted_coupling_region_monitor.open(
            output_path(runtime, "accepted_coupling_residual_by_region.dat").c_str());
        accepted_coupling_region_monitor
            << "# step time_fs region beam_front_ix wave_core_end_m "
            << "sum_rJ sum_abs_rJ L2_rJ Linf_rJ "
            << "sum_rK sum_abs_rK L2_rK Linf_rK "
            << "max_abs_JN_minus_GstarJE max_difference_face "
            << "face_work_JN_global delta_K_bkg_global fct_work_global "
            << "face_count_global cell_count_global "
            << "stage5_R_couple_global region_sum_rJ_global "
            << "rJ_identity_error_global stage5_R_FV_global "
            << "region_sum_rK_global rK_identity_error_global "
            << "energy_pair_face_work_JN_global "
            << "region_sum_face_work_JN_global "
            << "face_work_JN_identity_error_global "
            << "R_couple_centered_region R_couple_upwind_region R_couple_fct_region "
            << "stage5_R_couple_centered_global region_R_couple_centered_sum "
            << "R_couple_centered_identity_error_global "
            << "stage5_R_couple_upwind_global region_R_couple_upwind_sum "
            << "R_couple_upwind_identity_error_global "
            << "stage5_R_couple_fct_global region_R_couple_fct_sum "
            << "R_couple_fct_identity_error_global "
            << "accepted state_advanced soft_unconverged\n";
        accepted_coupling_region_monitor << std::scientific
                                         << std::setprecision(10);

        fct_macro_budget_monitor.open(
            output_path(runtime, "fct_macro_budget_by_region.dat").c_str());
        fct_macro_budget_monitor
            << "# step time_fs accepted state_advanced soft_unconverged "
            << "flux_direction x_region velocity_region face_count active_face_count min_alpha "
            << "raw_abs_mass rejected_abs_mass rejected_mass_fraction "
            << "raw_abs_current rejected_abs_current rejected_current_fraction "
            << "raw_abs_energy rejected_abs_energy rejected_energy_fraction "
            << "delta_N_fct delta_J_fct delta_K_fct E_dot_J_fct R_fct_E\n";
        fct_macro_budget_monitor << std::scientific << std::setprecision(10);

        high_order_candidate_monitor.open(
            output_path(runtime, "high_order_candidate_negativity.dat").c_str());
        high_order_candidate_monitor
            << "# step time_fs accepted state_advanced soft_unconverged "
            << "reconstructed_face_count negative_reconstructed_face_count "
            << "negative_reconstructed_mass_weight "
            << "negative_reconstructed_current_weight "
            << "negative_reconstructed_energy_weight "
            << "high_candidate_negative_cell_count high_candidate_negative_mass "
            << "high_candidate_negative_current_weight "
            << "high_candidate_negative_energy_weight "
            << "low_order_negative_cell_count low_order_negative_mass\n";
        high_order_candidate_monitor << std::scientific << std::setprecision(10);

        final_dual_u_pairing_monitor.open(
            output_path(runtime, "final_dual_u_pairing.dat").c_str());
        final_dual_u_pairing_monitor
             << "# step time_fs accepted state_advanced soft_unconverged "
             << "valid target_linf residual_before_linf residual_after_linf "
             << "residual_after_relative_target residual_reduction_fraction "
             << "minimum_scale correction_l2 correction_linf candidate_min "
            << "corrected_cell_count limited_cell_count unresolved_cell_count "
            << "face_pairing_attempted face_pairing_accepted "
            << "face_pairing_fallback face_solver_converged "
            << "face_solver_iterations unresolved_mode_count "
            << "face_residual_before face_residual_after "
            << "face_core_residual_before face_core_residual_after "
            << "unresolved_mode_l2 face_correction_l2 "
            << "face_correction_linf capacity_active_cells "
            << "trust_region_active_cells delta_ke delta_work "
            << "face_candidate_min face_mass_error candidate_valid "
            << "face_mass_relative_error face_cell_mass_error_linf "
            << "face_cell_mass_relative_linf face_energy_pair_error "
            << "face_energy_pair_relative face_energy_residual_scale "
            << "face_energy_residual_ratio correction_trust_limit "
            << "correction_trust_ratio face_f_residual_relative_growth "
            << "candidate_face_residual_after "
            << "candidate_core_face_residual_after candidate_delta_ke "
            << "candidate_delta_work candidate_mass_error "
            << "candidate_min_before_fallback requested_correction_l2 "
            << "requested_correction_linf applied_correction_l2 "
            << "applied_correction_linf nonzero_capacity_cells "
            << "bound_saturated_cells objective_residual "
            << "objective_smoothness objective_amplitude objective_total "
            << "rejection_mask pass_solver pass_apply "
            << "pass_global_residual pass_core_residual "
            << "pass_correction_trust pass_energy_pair "
            << "pass_energy_residual_scale pass_candidate_min pass_mass "
            << "pass_f_residual\n";
        final_dual_u_pairing_monitor
            << std::scientific << std::setprecision(10);

        if (runtime.diagnostic_level >= 2) {
            midpoint_iteration_residual_monitor.open(
                output_path(runtime, "trial_midpoint_residuals.dat").c_str());
            midpoint_iteration_residual_monitor
                << "# step time_fs iteration residual_E residual_J_bkg "
                << "residual_E_contraction residual_J_bkg_contraction "
                << "record_interval residual_J_beam residual_f accepted state_advanced "
                << "converged soft_unconverged acceleration_coefficient "
                << "acceleration_residual_before acceleration_status\n";
            midpoint_iteration_residual_monitor << std::scientific
                                                << std::setprecision(17);
        }

        if (write_beam_ledger) {
            beam_half_step_ledger.open(
                output_path(runtime, "beam_half_step_ledger.dat").c_str());
            beam_half_step_ledger
                << "# accepted_step time_fs dt_s N_in N_out injected_energy "
                << "outflow_energy injected_current boundary_current_out "
                << "injected_current_impulse outflow_current_impulse "
                << "physical_substep\n";
            beam_half_step_ledger << std::scientific << std::setprecision(17);
        }

        std::ofstream f_neg_monitor;
        f_neg_monitor.open(output_path(runtime, "f_negativity_monitor.dat").c_str());
        f_neg_monitor
            << "# accepted_snapshot: statistics scan final accepted bkg_e.f; "
            << "neg_ratio_max uses max positive f as scale\n"
            << "# step  time[fs]  accepted  state_advanced  "
            << "soft_unconverged  min_f  neg_ratio_max  "
            << "neg_mass_total[m^-2]  neg_cell_count  "
            << "x_worst  u_worst  mu_worst\n";
        f_neg_monitor << std::scientific << std::setprecision(8);
        f_neg_monitor.close();

        trial_retry_monitor.open(
            output_path(runtime, "trial_retry_diagnostics.dat").c_str());
        trial_retry_monitor
            << "# physical_step time_fs trial retry_index split_depth "
            << "acceptance_kind residual_E residual_J_bkg residual_f "
            << "beam_continuity_valid transport_safe energy_eta "
            << "energy_residual energy_limit accepted reason\n";
        trial_retry_monitor << std::scientific << std::setprecision(17);

        physical_interval_acceptance_monitor.open(
            output_path(runtime, "physical_interval_acceptance.dat").c_str());
        physical_interval_acceptance_monitor
            << "# physical_step time_begin_fs time_end_fs dt_s accepted_substeps "
            << "strict soft retry_count split_count wall_time_retry_s "
            << "operator_evaluations_retry residual_E residual_J_bkg residual_f\n";
        physical_interval_acceptance_monitor << std::scientific
                                             << std::setprecision(17);

    }
    long long accepted_after_restart = 0;
    long long performance_total_nonlinear_iterations = 0;
    long long performance_total_operator_evaluations = 0;
    long long performance_midpoint_predictor_used_steps = 0;
    long long performance_strict_accepted_steps = 0;
    long long performance_soft_accepted_steps = 0;
    long long performance_acceleration_attempts = 0;
    long long performance_acceleration_accepted = 0;
    long long performance_acceleration_fallback_evaluations = 0;
    long long performance_acceleration_rejected_residual = 0;
    long long performance_acceleration_rejected_nonfinite = 0;
    long long performance_acceleration_rejected_hard_failure = 0;
    long long performance_acceleration_rejected_coefficient = 0;
    long long performance_acceleration_history_resets = 0;
    long long performance_face_pairing_attempted_steps = 0;
    long long performance_face_pairing_accepted_steps = 0;
    long long performance_face_pairing_fallback_steps = 0;
    int soft_core_macro_debt_consecutive_steps = 0;
    long long performance_retry_count = 0;
    long long performance_split_count = 0;
    long long performance_operator_evaluations_retry = 0;
    double performance_wall_time_retry = 0.0;
    double performance_final_residual_e = 0.0;
    double performance_final_residual_j_bkg = 0.0;
    double performance_final_residual_f = 0.0;
    double performance_max_residual_e = 0.0;
    double performance_max_residual_j_bkg = 0.0;
    double performance_max_residual_f = 0.0;
    double performance_limiter_active_sum = 0.0;
    double performance_limiter_active_min =
        std::numeric_limits<double>::infinity();
    double performance_limiter_active_max = 0.0;
    double performance_limiter_min_alpha = 1.0;
    double performance_x_limiter_active_sum = 0.0;
    double performance_x_limiter_active_min =
        std::numeric_limits<double>::infinity();
    double performance_x_limiter_active_max = 0.0;
    double performance_x_limiter_min_alpha = 1.0;
    const double performance_start_time_s = current_time;
    MPI_Barrier(MPI_COMM_WORLD);
    const double performance_wall_start = MPI_Wtime();
    size_t next_checkpoint = 0;
    while (next_checkpoint < runtime.checkpoint_times_fs.size() &&
           current_time / Const::femto >= runtime.checkpoint_times_fs[next_checkpoint]) ++next_checkpoint;
    bool fixed_midpoint_face_pairing_started = false;
    for (int step = 1; step <= nsteps; ++step) {
        const long long physical_step = restored_step + step;
        double time = current_time + dt;
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
        // Probe/audit runs need every accepted step, independent of the
        // production cadence, so their field/current closure can be replayed.
        const bool collect_step_diagnostics = runtime.diagnostic_level >= 2 ||
            should_write_step_diagnostics(config, step) ||
            time >= next_snapshot || step == nsteps;
        const double time_fs_for_audit = time / Const::femto;
        const bool accepted_energy_audit_in_window =
            runtime.accepted_energy_audit_start_fs < 0.0 ||
            (time_fs_for_audit >= runtime.accepted_energy_audit_start_fs &&
             time_fs_for_audit <= runtime.accepted_energy_audit_end_fs);
        const bool collect_accepted_energy_audit =
            accepted_energy_audit_in_window &&
            (physical_step % runtime.accepted_energy_audit_cadence == 0);
        const bool collect_energy_accounting =
            collect_step_diagnostics || collect_accepted_energy_audit;
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
        double local_bkg_energy_residual_step = 0.0;
        double bkg_energy_residual_step = 0.0;
        double bkg_energy_relative_residual_step = 0.0;
        int coupled_iter_step = 0;
        double coupled_residual_E_step = 0.0;
        double coupled_residual_J_bkg_step = 0.0;
        double coupled_residual_J_beam_step = 0.0;
        BackgroundCurrentDiagnostics local_bkg_current_diag_step;
        reset_background_current_diagnostics(local_bkg_current_diag_step);
        const Species bkg_step_start = bkg_e;
        const BeamPIC beam_step_start = beam;
        const EMFields fields_step_start = fields;
        if (collect_energy_accounting) {
            bkg_e.total_particle_number_and_energy(bkg_number_step_start,
                                                   bkg_ke_step_start);
            double global_bkg_start_values[2] = {
                bkg_number_step_start,
                bkg_ke_step_start
            };
            MPI_Allreduce(MPI_IN_PLACE, global_bkg_start_values, 2,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            global_bkg_ke_step_start = global_bkg_start_values[1];
            beam_ke_step_start = beam.total_kinetic_energy();
            field_energy_step_start = fields.total_energy();
        }

        trace_progress(config, mpi_rank, step,
                       "before coupled midpoint FV solve");
        const double time_fs_for_trace = time / Const::femto;
        const bool midpoint_live_trace =
            runtime.midpoint_trace_start_fs >= 0.0 &&
            time_fs_for_trace >= runtime.midpoint_trace_start_fs &&
            time_fs_for_trace <= runtime.midpoint_trace_end_fs;
        if (midpoint_live_trace && mpi_rank == 0) {
            std::printf(
                "[step-live] physical_step=%lld restart_relative_step=%d "
                "t_fs=%.16e stage=midpoint_begin\n",
                physical_step, step, time_fs_for_trace);
            std::fflush(stdout);
        }
        midpoint_solver.set_step_diagnostics_enabled(collect_step_diagnostics);
        midpoint_solver.set_accepted_energy_audit_enabled(
            collect_accepted_energy_audit);
        midpoint_solver.set_midpoint_iteration_trace(
            runtime.diagnostic_level >= 2 || collect_accepted_energy_audit);
        VlasovAmpereMidpointSolver::Result midpoint_result =
            midpoint_solver.advance_background_and_fields(
                bkg_step_start, beam_step_start, fields_step_start, sgrid,
                dt, time, mpi_rank, mpi_size);
        if (midpoint_live_trace && mpi_rank == 0) {
            std::printf(
                "[step-live] physical_step=%lld restart_relative_step=%d "
                "t_fs=%.16e stage=midpoint_end iterations=%d substeps=%d "
                "converged=%d soft=%d failed=%d state_advanced=%d\n",
                physical_step, step, time_fs_for_trace,
                midpoint_result.nonlinear_iterations,
                midpoint_result.substeps_used,
                midpoint_result.converged ? 1 : 0,
                midpoint_result.soft_unconverged ? 1 : 0,
                midpoint_result.failed ? 1 : 0,
                midpoint_result.state_advanced);
            std::fflush(stdout);
        }
        trace_progress(config, mpi_rank, step,
                       "after coupled midpoint FV solve");
        if (mpi_rank == 0 && runtime.diagnostic_level >= 2) {
            const size_t count = midpoint_result.midpoint_residual_e_history.size();
            double previous_record_e = 0.0;
            double previous_record_j_bkg = 0.0;
            size_t previous_record_iteration = 0;
            bool have_previous_record = false;
            for (size_t iter = 0; iter < count; ++iter) {
                // Level 2 is an explicit nonlinear-solver audit.  Preserve
                // every candidate so an accelerated field proposal and its
                // subsequent accept/reject evaluation stay adjacent.
                const bool record_iteration = true;
                if (!record_iteration) continue;
                const double residual_e = midpoint_result.midpoint_residual_e_history[iter];
                const double residual_j_bkg =
                    midpoint_result.midpoint_residual_j_bkg_history[iter];
                const double residual_e_contraction = have_previous_record
                    ? (previous_record_e != 0.0
                        ? residual_e / previous_record_e
                        : (residual_e == 0.0 ? 1.0
                                             : std::numeric_limits<double>::infinity()))
                    : 1.0;
                const double residual_j_bkg_contraction = have_previous_record
                    ? (previous_record_j_bkg != 0.0
                        ? residual_j_bkg / previous_record_j_bkg
                        : (residual_j_bkg == 0.0 ? 1.0
                                                 : std::numeric_limits<double>::infinity()))
                    : 1.0;
                const size_t record_interval = have_previous_record
                    ? iter + 1 - previous_record_iteration : 0;
                const double residual_j_beam =
                    iter < midpoint_result.midpoint_residual_j_beam_history.size()
                    ? midpoint_result.midpoint_residual_j_beam_history[iter]
                    : std::numeric_limits<double>::quiet_NaN();
                const double residual_f =
                    iter < midpoint_result.midpoint_residual_f_history.size()
                    ? midpoint_result.midpoint_residual_f_history[iter]
                    : std::numeric_limits<double>::quiet_NaN();
                const double acceleration_coefficient =
                    iter < midpoint_result.midpoint_acceleration_omega_history.size()
                    ? midpoint_result.midpoint_acceleration_omega_history[iter]
                    : std::numeric_limits<double>::quiet_NaN();
                const double acceleration_residual_before =
                    iter < midpoint_result.midpoint_acceleration_residual_before_history.size()
                    ? midpoint_result.midpoint_acceleration_residual_before_history[iter]
                    : std::numeric_limits<double>::quiet_NaN();
                const int acceleration_status =
                    iter < midpoint_result.midpoint_acceleration_status_history.size()
                    ? midpoint_result.midpoint_acceleration_status_history[iter]
                    : -1;
                midpoint_iteration_residual_monitor
                    << physical_step << " " << time / Const::femto << " "
                    << (iter + 1) << " "
                    << residual_e << " " << residual_j_bkg << " "
                    << residual_e_contraction << " "
                    << residual_j_bkg_contraction << " "
                    << record_interval << " " << residual_j_beam << " "
                    << residual_f << " "
                    << (midpoint_result.state_advanced && !midpoint_result.failed ? 1 : 0)
                    << " " << midpoint_result.state_advanced << " "
                    << midpoint_result.converged << " "
                    << midpoint_result.soft_unconverged << " "
                    << acceleration_coefficient << " "
                    << acceleration_residual_before << " "
                    << acceleration_status << "\n";
                previous_record_e = residual_e;
                previous_record_j_bkg = residual_j_bkg;
                previous_record_iteration = iter + 1;
                have_previous_record = true;
            }
            midpoint_iteration_residual_monitor.flush();
        }

        // A max-iteration result is a transaction candidate, not an accepted
        // state.  The compact energy account below is intentionally evaluated
        // only for that final candidate, never inside the Picard loop.
        double soft_energy_eta = std::numeric_limits<double>::infinity();
        double soft_energy_residual = std::numeric_limits<double>::quiet_NaN();
        const double soft_energy_limit = std::min(
            5.0 * runtime.soft_candidate_energy_p99_reference,
            runtime.soft_candidate_energy_absolute_limit);
        const auto residual_history_is_contracting =
            [](const std::vector<double>& history) {
                if (history.size() < 4) return false;
                const size_t begin = history.size() - 4;
                bool continuously_growing = true;
                double log_sum = 0.0;
                int ratio_count = 0;
                for (size_t i = begin + 1; i < history.size(); ++i) {
                    if (!std::isfinite(history[i]) ||
                        !std::isfinite(history[i - 1])) return false;
                    continuously_growing = continuously_growing &&
                        history[i] > history[i - 1];
                    if (history[i - 1] > 0.0 && history[i] >= 0.0) {
                        const double ratio = history[i] / history[i - 1];
                        if (!std::isfinite(ratio)) return false;
                        log_sum += std::log(std::max(ratio,
                            std::numeric_limits<double>::min()));
                        ++ratio_count;
                    }
                }
                if (continuously_growing || ratio_count == 0) return false;
                const double last_ratio = history.back() /
                    std::max(history[history.size() - 2],
                             std::numeric_limits<double>::min());
                return std::isfinite(last_ratio) && last_ratio <= 1.0 &&
                    std::exp(log_sum / static_cast<double>(ratio_count)) < 1.0;
            };
        const auto evaluate_soft_candidate =
            [&](VlasovAmpereMidpointSolver::Result& candidate,
                const Species& candidate_start_bkg,
                const BeamPIC& candidate_start_beam,
                const EMFields& candidate_start_fields,
                double candidate_time, double candidate_dt,
                int retry_index, int split_depth, const char* trial_name) {
                bool accepted = false;
                soft_energy_eta = std::numeric_limits<double>::infinity();
                soft_energy_residual = std::numeric_limits<double>::quiet_NaN();
                const bool residuals_ok =
                    std::isfinite(candidate.residual_E) &&
                    std::isfinite(candidate.residual_J_bkg) &&
                    std::isfinite(candidate.residual_f) &&
                    candidate.residual_E <= runtime.soft_candidate_field_tolerance &&
                    candidate.residual_J_bkg <= runtime.soft_candidate_current_tolerance;
                const bool histories_ok =
                    residual_history_is_contracting(
                        candidate.midpoint_residual_e_history) &&
                    residual_history_is_contracting(
                        candidate.midpoint_residual_j_bkg_history) &&
                    !candidate.midpoint_residual_f_history.empty() &&
                    std::isfinite(candidate.midpoint_residual_f_history.back()) &&
                    !((candidate.midpoint_residual_f_history.size() >= 4) &&
                      candidate.midpoint_residual_f_history.back() >
                          candidate.midpoint_residual_f_history[
                              candidate.midpoint_residual_f_history.size() - 2] &&
                      candidate.midpoint_residual_f_history[
                              candidate.midpoint_residual_f_history.size() - 2] >
                          candidate.midpoint_residual_f_history[
                              candidate.midpoint_residual_f_history.size() - 3] &&
                      candidate.midpoint_residual_f_history[
                              candidate.midpoint_residual_f_history.size() - 3] >
                          candidate.midpoint_residual_f_history[
                              candidate.midpoint_residual_f_history.size() - 4]);
                const bool basic_ok = !candidate.failed &&
                    candidate.acceptance_kind ==
                        VlasovAmpereMidpointSolver::SOFT_CANDIDATE &&
                    candidate.transport_safe != 0 &&
                    candidate.beam_continuity_valid != 0 &&
                    residuals_ok && histories_ok;
                if (basic_ok) {
                    double start_number = 0.0;
                    double start_bkg_ke = 0.0;
                    double end_number = 0.0;
                    double end_bkg_ke = 0.0;
                    candidate_start_bkg.total_particle_number_and_energy(
                        start_number, start_bkg_ke);
                    Species energy_candidate = candidate.species_np1;
                    double candidate_collision_energy = 0.0;
                    if (energy_candidate.collisions_enabled) {
                        candidate_collision_energy += collision.apply(
                            energy_candidate, candidate_dt, Param::dens,
                            Param::temperature_e, Const::me, 1.0, 1.0);
                        candidate_collision_energy += collision.apply(
                            energy_candidate, candidate_dt,
                            Param::dens / Param::Z_ion, Param::temperature_i,
                            Param::mass_ion, static_cast<double>(Param::Z_ion),
                            1.0);
                    }
                    energy_candidate.total_particle_number_and_energy(
                        end_number, end_bkg_ke);
                    double packed[9] = {
                        start_bkg_ke, end_bkg_ke,
                        candidate_start_beam.total_kinetic_energy(),
                        candidate.beam_np1.total_kinetic_energy(),
                        candidate_start_fields.total_energy(),
                        candidate.fields_np1.total_energy(),
                        candidate.beam_np1.last_injected_energy(),
                        candidate.beam_np1.last_outflow_energy(),
                        candidate_collision_energy
                    };
                    MPI_Allreduce(MPI_IN_PLACE, packed, 9, MPI_DOUBLE, MPI_SUM,
                                  MPI_COMM_WORLD);
                    const double delta_bkg = packed[1] - packed[0];
                    const double delta_beam = packed[3] - packed[2];
                    const double delta_field = packed[5] - packed[4];
                    soft_energy_residual = delta_bkg + delta_beam + delta_field -
                        packed[6] + packed[7] - packed[8];
                    const double energy_scale = std::fabs(delta_bkg) +
                        std::fabs(delta_beam) + std::fabs(delta_field) +
                        std::fabs(packed[6]) + std::fabs(packed[7]) +
                        std::fabs(packed[8]) + 1.0;
                    soft_energy_eta = std::fabs(soft_energy_residual) / energy_scale;
                    accepted = std::isfinite(soft_energy_eta) &&
                        soft_energy_eta <= soft_energy_limit;
                }
                if (mpi_rank == 0 && trial_retry_monitor.is_open()) {
                    trial_retry_monitor << physical_step << " "
                        << candidate_time / Const::femto << " " << trial_name << " "
                        << retry_index << " " << split_depth << " "
                        << static_cast<int>(candidate.acceptance_kind) << " "
                        << candidate.residual_E << " "
                        << candidate.residual_J_bkg << " "
                        << candidate.residual_f << " "
                        << candidate.beam_continuity_valid << " "
                        << candidate.transport_safe << " " << soft_energy_eta
                        << " " << soft_energy_residual << " "
                        << soft_energy_limit << " " << (accepted ? 1 : 0)
                        << " " << (basic_ok ? "energy_gate" : "candidate_gate")
                        << "\n";
                    trial_retry_monitor.flush();
                }
                if (accepted) {
                    candidate.state_advanced = 1;
                    candidate.soft_accepted = true;
                } else if (!candidate.failed) {
                    candidate.acceptance_kind =
                        VlasovAmpereMidpointSolver::RETRY_REQUIRED;
                }
                return accepted;
            };

        bool midpoint_soft_accepted = false;
        int retry_count_step = 0;
        int split_count_step = 0;
        int operator_evaluations_retry_step = 0;
        double wall_time_retry_step = 0.0;
        bool split_collisions_applied = false;
        double split_beam_injected_number = 0.0;
        double split_beam_outflow_number = 0.0;
        double split_beam_injected_energy = 0.0;
        double split_beam_outflow_energy = 0.0;
        double split_beam_injected_impulse = 0.0;
        double split_beam_outflow_impulse = 0.0;
        std::array<double, 12> split_beam_ledger_rows = {{0.0}};
        if (midpoint_result.acceptance_kind ==
            VlasovAmpereMidpointSolver::SOFT_CANDIDATE) {
            midpoint_soft_accepted = evaluate_soft_candidate(
                midpoint_result, bkg_step_start, beam_step_start,
                fields_step_start, time, dt, retry_count_step, 0, "initial");
            if (!midpoint_soft_accepted &&
                midpoint_result.midpoint_predictor_used != 0) {
                const double retry_wall_start = MPI_Wtime();
                const VlasovAmpereMidpointSolver::MidpointInitialGuessMode saved_mode =
                    midpoint_solver.midpoint_initial_guess_mode();
                midpoint_solver.set_midpoint_initial_guess_mode(
                    VlasovAmpereMidpointSolver::MIDPOINT_INITIAL_GUESS_NONE);
                ++retry_count_step;
                midpoint_result = midpoint_solver.advance_background_and_fields(
                    bkg_step_start, beam_step_start, fields_step_start, sgrid,
                    dt, time, mpi_rank, mpi_size);
                midpoint_solver.set_midpoint_initial_guess_mode(saved_mode);
                operator_evaluations_retry_step +=
                    midpoint_result.operator_evaluations;
                wall_time_retry_step += MPI_Wtime() - retry_wall_start;
                if (midpoint_result.acceptance_kind ==
                    VlasovAmpereMidpointSolver::SOFT_CANDIDATE) {
                    midpoint_soft_accepted = evaluate_soft_candidate(
                        midpoint_result, bkg_step_start, beam_step_start,
                        fields_step_start, time, dt, retry_count_step, 0,
                        "predictor_none");
                }
            }
        }
        if (!midpoint_soft_accepted && !midpoint_result.failed &&
            midpoint_result.acceptance_kind ==
                VlasovAmpereMidpointSolver::RETRY_REQUIRED) {
            // A split is two complete physical intervals.  Nothing below
            // touches bkg_e/beam/fields until both local candidates pass.
            const double split_wall_start = MPI_Wtime();
            const double half_dt = 0.5 * dt;
            const VlasovAmpereMidpointSolver::MidpointInitialGuessMode saved_mode =
                midpoint_solver.midpoint_initial_guess_mode();
            midpoint_solver.set_midpoint_initial_guess_mode(
                VlasovAmpereMidpointSolver::MIDPOINT_INITIAL_GUESS_NONE);
            VlasovAmpereMidpointSolver::Result first_half =
                midpoint_solver.advance_background_and_fields(
                    bkg_step_start, beam_step_start, fields_step_start, sgrid,
                    half_dt, current_time + half_dt, mpi_rank, mpi_size);
            bool first_ok = first_half.converged && !first_half.failed &&
                first_half.state_advanced != 0;
            if (!first_ok && first_half.acceptance_kind ==
                VlasovAmpereMidpointSolver::SOFT_CANDIDATE) {
                first_ok = evaluate_soft_candidate(
                    first_half, bkg_step_start, beam_step_start,
                    fields_step_start, current_time + half_dt, half_dt, retry_count_step,
                    1, "split_half_1");
            }
            double split_collision_energy = 0.0;
            VlasovAmpereMidpointSolver::Result second_half;
            bool second_ok = false;
            if (first_ok) {
                Species split_bkg_mid = first_half.species_np1;
                if (split_bkg_mid.collisions_enabled) {
                    split_collision_energy += collision.apply(
                        split_bkg_mid, half_dt, Param::dens,
                        Param::temperature_e, Const::me, 1.0, 1.0);
                    split_collision_energy += collision.apply(
                        split_bkg_mid, half_dt, Param::dens / Param::Z_ion,
                        Param::temperature_i, Param::mass_ion,
                        static_cast<double>(Param::Z_ion), 1.0);
                }
                second_half = midpoint_solver.advance_background_and_fields(
                    split_bkg_mid, first_half.beam_np1, first_half.fields_np1,
                    sgrid, half_dt, time, mpi_rank, mpi_size);
                second_ok = second_half.converged && !second_half.failed &&
                    second_half.state_advanced != 0;
                if (!second_ok && second_half.acceptance_kind ==
                    VlasovAmpereMidpointSolver::SOFT_CANDIDATE) {
                    second_ok = evaluate_soft_candidate(
                        second_half, split_bkg_mid, first_half.beam_np1,
                        first_half.fields_np1, time, half_dt, retry_count_step, 1,
                        "split_half_2");
                }
                if (second_ok && second_half.species_np1.collisions_enabled) {
                    split_collision_energy += collision.apply(
                        second_half.species_np1, half_dt, Param::dens,
                        Param::temperature_e, Const::me, 1.0, 1.0);
                    split_collision_energy += collision.apply(
                        second_half.species_np1, half_dt, Param::dens / Param::Z_ion,
                        Param::temperature_i, Param::mass_ion,
                        static_cast<double>(Param::Z_ion), 1.0);
                }
            }
            midpoint_solver.set_midpoint_initial_guess_mode(saved_mode);
            operator_evaluations_retry_step += first_half.operator_evaluations +
                (first_ok ? second_half.operator_evaluations : 0);
            wall_time_retry_step += MPI_Wtime() - split_wall_start;
            if (first_ok && second_ok) {
                split_beam_injected_number =
                    first_half.beam_np1.last_injected_number() +
                    second_half.beam_np1.last_injected_number();
                split_beam_outflow_number =
                    first_half.beam_np1.last_outflow_number() +
                    second_half.beam_np1.last_outflow_number();
                split_beam_injected_energy =
                    first_half.beam_np1.last_injected_energy() +
                    second_half.beam_np1.last_injected_energy();
                split_beam_outflow_energy =
                    first_half.beam_np1.last_outflow_energy() +
                    second_half.beam_np1.last_outflow_energy();
                split_beam_injected_impulse = half_dt *
                    (first_half.beam_np1.last_injected_current() +
                     second_half.beam_np1.last_injected_current());
                split_beam_outflow_impulse = half_dt *
                    (first_half.beam_np1.last_outflow_current() +
                     second_half.beam_np1.last_outflow_current());
                const BeamPIC* split_beams[2] = {
                    &first_half.beam_np1, &second_half.beam_np1};
                for (int half = 0; half < 2; ++half) {
                    const BeamPIC& split_beam = *split_beams[half];
                    const size_t row = static_cast<size_t>(6 * half);
                    split_beam_ledger_rows[row] =
                        split_beam.last_injected_number();
                    split_beam_ledger_rows[row + 1] =
                        split_beam.last_outflow_number();
                    split_beam_ledger_rows[row + 2] =
                        split_beam.last_injected_energy();
                    split_beam_ledger_rows[row + 3] =
                        split_beam.last_outflow_energy();
                    split_beam_ledger_rows[row + 4] = half_dt *
                        split_beam.last_injected_current();
                    split_beam_ledger_rows[row + 5] = half_dt *
                        split_beam.last_outflow_current();
                }
                midpoint_result = second_half;
                collision_energy_step = split_collision_energy;
                split_collisions_applied = true;
                split_count_step = 1;
                midpoint_soft_accepted = midpoint_result.soft_unconverged &&
                    midpoint_result.state_advanced != 0 && !midpoint_result.failed;
            } else {
                midpoint_result = first_ok ? second_half : first_half;
                midpoint_result.acceptance_kind =
                    VlasovAmpereMidpointSolver::RETRY_REQUIRED;
            }
        }
        if ((!midpoint_result.converged && !midpoint_soft_accepted) ||
            midpoint_result.failed) {
            const char* failure_reason = "NONE";
            switch (midpoint_result.failure_reason) {
            case 1: failure_reason = "CFL_LIMIT"; break;
            case 2: failure_reason = "LOW_ORDER_POSITIVITY"; break;
            case 3: failure_reason = "FCT_FINAL_POSITIVITY"; break;
            case 4: failure_reason = "NONFINITE_STATE"; break;
            case 5: failure_reason = "FCT_HIGH_LOW_IDENTITY"; break;
            case 6: failure_reason = "FCT_DONOR_CAPACITY"; break;
            case 7: failure_reason = "FCT_INTERFACE_CHECKSUM"; break;
            case 12: failure_reason = "NONUNIFORM_HIGH_ORDER_DISABLED"; break;
            case 17: failure_reason = "NONFINITE_BEAM_CONTINUITY"; break;
            default: break;
            }
            if (mpi_rank == 0) {
                if (midpoint_result.failure_reason != 0) {
                std::fprintf(stderr,
                             "EARLY_EXIT: step=%d reason=%s iter=%d substep=%d "
                             "global_cfl=%.16e low_min=%.16e final_min=%.16e "
                             "worst=(x=%d,j=%d,k=%d); "
                             "fct_identity=%.16e identity_violation=%.16e "
                             "donor_violation=%.16e interface_checksum=%.16e "
                             "interface_violation=%.16e low_tol=%.16e "
                             "final_tol=%.16e; "
                             "identity_worst=(x=%d,j=%d,k=%d,R=%.16e,S=%.16e,"
                             "R_over_S=%.16e,abs_R_over_max1S=%.16e); "
                             "donor_worst=(x=%d,j=%d,k=%d,m_low=%.16e,"
                             "A=[%.16e,%.16e,%.16e,%.16e],"
                             "alpha=[%.16e,%.16e,%.16e,%.16e],"
                             "outflow=%.16e,S=%.16e,relative=%.16e,"
                             "roundoff_warning=%d,beta_count=%lld,"
                             "beta_min=%.16e); "
                             "low_order=(M_in=%.16e,M_low=%.16e,"
                             "hPhi=[%.16e,%.16e,%.16e,%.16e],"
                             "out=[%.16e,%.16e,%.16e,%.16e],"
                             "in=[%.16e,%.16e,%.16e,%.16e],"
                             "CFL_x=%.16e,CFL_u=%.16e,"
                             "work_input_min=%.16e,mpi_interface=%d).\n",
                             step, failure_reason,
                             midpoint_result.failure_iteration,
                             midpoint_result.failure_substep,
                             midpoint_result.failure_global_cfl,
                             midpoint_result.failure_low_min,
                             midpoint_result.failure_final_min,
                             midpoint_result.failure_worst_ix,
                             midpoint_result.failure_worst_iv,
                             midpoint_result.failure_worst_imu,
                             midpoint_result.fct_high_low_identity_linf,
                             midpoint_result.fct_high_low_identity_violation,
                             midpoint_result.fct_donor_capacity_violation,
                             midpoint_result.fct_interface_checksum_linf,
                             midpoint_result.fct_interface_checksum_violation,
                             midpoint_result.fct_low_order_tolerance_linf,
                             midpoint_result.fct_final_tolerance_linf,
                             midpoint_result.fct_high_low_identity_worst_ix,
                             midpoint_result.fct_high_low_identity_worst_iv,
                             midpoint_result.fct_high_low_identity_worst_imu,
                             midpoint_result.fct_high_low_identity_worst_residual,
                             midpoint_result.fct_high_low_identity_worst_scale,
                             midpoint_result.fct_high_low_identity_worst_relative,
                             midpoint_result.fct_high_low_identity_ratio_linf,
                             midpoint_result.fct_donor_worst_ix,
                             midpoint_result.fct_donor_worst_iv,
                             midpoint_result.fct_donor_worst_imu,
                             midpoint_result.fct_donor_worst_m_low,
                             midpoint_result.fct_donor_worst_ax_left,
                             midpoint_result.fct_donor_worst_ax_right,
                             midpoint_result.fct_donor_worst_au_lower,
                             midpoint_result.fct_donor_worst_au_upper,
                             midpoint_result.fct_donor_worst_alpha_x_left,
                             midpoint_result.fct_donor_worst_alpha_x_right,
                             midpoint_result.fct_donor_worst_alpha_u_lower,
                             midpoint_result.fct_donor_worst_alpha_u_upper,
                             midpoint_result.fct_donor_worst_outflow,
                             midpoint_result.fct_donor_worst_scale,
                             midpoint_result.fct_donor_worst_relative,
                             midpoint_result.fct_donor_roundoff_warning,
                             midpoint_result.fct_donor_beta_applied_count,
                             midpoint_result.fct_donor_beta_min,
                             midpoint_result.failure_low_m_in,
                             midpoint_result.failure_low_m_low,
                             midpoint_result.failure_low_transfer_x_left,
                             midpoint_result.failure_low_transfer_x_right,
                             midpoint_result.failure_low_transfer_u_lower,
                             midpoint_result.failure_low_transfer_u_upper,
                             midpoint_result.failure_low_out_x_left,
                             midpoint_result.failure_low_out_x_right,
                             midpoint_result.failure_low_out_u_lower,
                             midpoint_result.failure_low_out_u_upper,
                             midpoint_result.failure_low_in_x_left,
                             midpoint_result.failure_low_in_x_right,
                             midpoint_result.failure_low_in_u_lower,
                             midpoint_result.failure_low_in_u_upper,
                             midpoint_result.failure_low_cfl_x,
                             midpoint_result.failure_low_cfl_u,
                             midpoint_result.failure_low_work_input_min,
                             midpoint_result.failure_low_on_mpi_interface);
                }
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
                              midpoint_result.finite_stage_failure_dmu_div);
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
                flush_full_beam_ledger();
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
            flush_full_beam_ledger();
            MPI_Abort(MPI_COMM_WORLD, 9);
            return 9;
        }
        if (midpoint_result.soft_unconverged) {
            const bool core_macro_debt =
                midpoint_result.neg_mass_core_fraction > 1.0e-6 ||
                midpoint_result.neg_energy_core_fraction > 1.0e-6 ||
                midpoint_result.neg_current_core_fraction > 1.0e-6;
            soft_core_macro_debt_consecutive_steps = core_macro_debt
                ? soft_core_macro_debt_consecutive_steps + 1 : 0;
            if (soft_core_macro_debt_consecutive_steps >= 3) {
                if (mpi_rank == 0) {
                    std::fprintf(stderr,
                        "ERROR: soft-candidate core macro debt exceeded the "
                        "three accepted-step limit at step %d.\n", step);
                }
                flush_full_beam_ledger();
                MPI_Abort(MPI_COMM_WORLD, 18);
                return 18;
            }
        } else {
            soft_core_macro_debt_consecutive_steps = 0;
        }

        bkg_e = midpoint_result.species_np1;
        beam = midpoint_result.beam_np1;
        fields = midpoint_result.fields_np1;
        const double accepted_beam_injected_number = split_count_step != 0
            ? split_beam_injected_number : beam.last_injected_number();
        const double accepted_beam_outflow_number = split_count_step != 0
            ? split_beam_outflow_number : beam.last_outflow_number();
        const double accepted_beam_injected_energy = split_count_step != 0
            ? split_beam_injected_energy : beam.last_injected_energy();
        const double accepted_beam_outflow_energy = split_count_step != 0
            ? split_beam_outflow_energy : beam.last_outflow_energy();
        const double accepted_beam_injected_impulse = split_count_step != 0
            ? split_beam_injected_impulse : beam.last_injected_current() * dt;
        const double accepted_beam_outflow_impulse = split_count_step != 0
            ? split_beam_outflow_impulse : beam.last_outflow_current() * dt;
        if (mpi_rank == 0 && physical_interval_acceptance_monitor.is_open()) {
            physical_interval_acceptance_monitor
                << physical_step << " " << current_time / Const::femto << " "
                << time / Const::femto << " " << dt << " "
                << (split_count_step != 0 ? 2 : 1) << " "
                << (midpoint_result.converged ? 1 : 0) << " "
                << (midpoint_result.soft_unconverged ? 1 : 0) << " "
                << retry_count_step << " " << split_count_step << " "
                << wall_time_retry_step << " "
                << operator_evaluations_retry_step << " "
                << midpoint_result.residual_E << " "
                << midpoint_result.residual_J_bkg << " "
                << midpoint_result.residual_f << "\n";
            physical_interval_acceptance_monitor.flush();
        }
        // The midpoint kernel only evaluates moments on diagnostic trials.
        // Refresh once for the accepted production state below, rather than
        // once for every Picard iterate.
        moments_current = false;
        if (runtime.diagnostic_level >= 2 && !midpoint_result.failed &&
            split_count_step == 0) {
            write_fixed_midpoint_face_pairing(
                output_path(runtime, "fixed_midpoint_face_pairing.dat"),
                physical_step, time, dt, midpoint_result, fields_step_start,
                midpoint_result.operator_input_fields_end_guess, sgrid,
                mpi_rank, mpi_size, fixed_midpoint_face_pairing_started);
            fixed_midpoint_face_pairing_started = true;
        }
        if (runtime.diagnostic_level >= 2 && runtime.dump_final_midpoint &&
            split_count_step == 0) {
            VlasovAmpereMidpointSolver::MidpointAuditState audit_state;
            audit_state.step = physical_step;
            audit_state.time_s = time;
            audit_state.dt_s = dt;
            audit_state.substeps_used = midpoint_result.substeps_used;
            audit_state.nonlinear_iterations = midpoint_result.nonlinear_iterations;
            audit_state.low_order_only = midpoint_solver.low_order_only();
            audit_state.high_order_enabled = midpoint_solver.nonuniform_high_order_enabled();
            audit_state.fct_enabled = midpoint_solver.fct_enabled();
            audit_state.background_coupling_mode =
                midpoint_result.background_coupling_mode;
            audit_state.acceptance_type = midpoint_result.converged
                ? "strict_converged"
                : (midpoint_result.soft_unconverged ? "soft_accepted" : "unconverged_best_trial");
            audit_state.bkg_n = bkg_step_start;
            audit_state.guess_np1 = bkg_e;
            audit_state.operator_input_guess = midpoint_result.operator_input_guess;
            audit_state.fields_n = fields_step_start;
            audit_state.fields_end_guess =
                midpoint_result.operator_input_fields_end_guess;
            audit_state.fields_np1 = fields;
            audit_state.coupling_layout.beam_front_ix =
                midpoint_result.coupling_beam_front_ix;
            audit_state.coupling_layout.wave_core_end_m =
                midpoint_result.coupling_wave_core_end_m;
            audit_state.j_beam_face_mid = midpoint_result.j_beam_face_mid;
            audit_state.reference_jn_face = midpoint_result.j_bkg_face_mid;
            // The existing energy-debug face quantity is the production
            // G*J_E reference. J_E cell storage is introduced with 10.1.
            audit_state.reference_je_cell = midpoint_result.j_bkg_energy_cell_mid;
            audit_state.reference_gstar_je_face = midpoint_result.j_bkg_energy_debug_face;
            audit_state.periodic_seam_face_audit =
                midpoint_result.periodic_seam_face_audit;
            audit_state.reference_stage5_r_fv = midpoint_result.stage5_r_fv;
            audit_state.reference_stage5_r_couple = midpoint_result.stage5_r_couple;
            audit_state.limiter_active_fraction = midpoint_result.limiter_active_fraction;
            audit_state.limiter_active_fraction_core =
                midpoint_result.limiter_active_fraction_core;
            audit_state.limiter_active_fraction_boundary =
                midpoint_result.limiter_active_fraction_boundary;
            audit_state.x_limiter_active_fraction =
                midpoint_result.x_limiter_active_fraction;
            audit_state.u_limiter_active_fraction =
                midpoint_result.u_limiter_active_fraction;
            audit_state.limiter_min_alpha = midpoint_result.limiter_min_alpha;
            if (!write_midpoint_audit_state(output_path(runtime, "final_midpoint"),
                                            audit_state, sgrid, mpi_rank,
                                            mpi_size, runtime_error)) {
                if (mpi_rank == 0) std::fprintf(stderr, "Midpoint audit error: %s\n", runtime_error.c_str());
                flush_full_beam_ledger();
                MPI_Abort(MPI_COMM_WORLD, 13);
                return 13;
            }
        }
        if (mpi_rank == 0 && collect_step_diagnostics &&
            midpoint_result.state_advanced != 0 && !midpoint_result.failed) {
            static const char* const region_names[6] = {
                "B_left", "Q_left", "wave_core", "quiet_right",
                "Q_right", "B_right"};
            for (int region = 0; region < 6; ++region) {
                const VlasovAmpereMidpointSolver::CouplingRegionDiagnostics& stats =
                    midpoint_result.coupling_regions[region];
                accepted_coupling_region_monitor
                    << step << " " << time / Const::femto << " "
                    << region_names[region] << " "
                    << midpoint_result.coupling_beam_front_ix << " "
                    << midpoint_result.coupling_wave_core_end_m << " "
                    << stats.sum_rj << " " << stats.sum_abs_rj << " "
                    << std::sqrt(stats.sum_sq_rj) << " " << stats.linf_rj << " "
                    << stats.sum_rk << " " << stats.sum_abs_rk << " "
                    << std::sqrt(stats.sum_sq_rk) << " " << stats.linf_rk << " "
                    << stats.max_abs_jn_minus_gstar_je << " "
                    << stats.max_abs_jn_minus_gstar_je_face << " "
                    << stats.face_work_jn << " " << stats.delta_ke_bkg << " "
                    << stats.fct_work << " " << stats.face_count << " "
                    << stats.cell_count << " "
                    << midpoint_result.stage5_r_couple << " "
                    << midpoint_result.coupling_rj_global_sum << " "
                    << midpoint_result.coupling_rj_reconstruction_error << " "
                    << midpoint_result.stage5_r_fv << " "
                    << midpoint_result.coupling_rk_global_sum << " "
                    << midpoint_result.coupling_rk_reconstruction_error << " "
                    << midpoint_result.current_diag.e_dot_j_charge << " "
                    << midpoint_result.coupling_face_work_jn_global_sum << " "
                    << midpoint_result
                           .coupling_face_work_jn_reconstruction_error << " "
                    << stats.r_couple_centered << " "
                    << stats.r_couple_upwind_stabilization << " "
                    << stats.r_couple_fct_stabilization << " "
                    << midpoint_result.stage5_r_couple_centered << " "
                    << midpoint_result.coupling_rj_centered_global_sum << " "
                    << midpoint_result.coupling_rj_centered_reconstruction_error << " "
                    << midpoint_result.stage5_r_couple_upwind_stabilization << " "
                    << midpoint_result.coupling_rj_upwind_global_sum << " "
                    << midpoint_result.coupling_rj_upwind_reconstruction_error << " "
                    << midpoint_result.stage5_r_couple_fct_stabilization << " "
                    << midpoint_result.coupling_rj_fct_global_sum << " "
                    << midpoint_result.coupling_rj_fct_reconstruction_error << " "
                    << 1 << " " << midpoint_result.state_advanced << " "
                    << midpoint_result.soft_unconverged << "\n";
            }
            accepted_coupling_region_monitor.flush();
        }
        if (mpi_rank == 0 && collect_step_diagnostics &&
            midpoint_result.state_advanced != 0 && !midpoint_result.failed &&
            midpoint_result.fct_macro_budget_valid != 0) {
            static const char* const x_region_names[3] = {
                "B_left", "core", "B_right"};
            static const char* const velocity_region_names[3] = {
                "velocity_core", "upar_8_to_10", "upar_tail_gt_10"};
            const std::array<VlasovAmpereMidpointSolver::FctMacroBudget,
                VlasovAmpereMidpointSolver::FCT_MACRO_REGION_COUNT>*
                direction_budgets[2] = {
                    &midpoint_result.fct_macro_budget_x,
                    &midpoint_result.fct_macro_budget_u};
            static const char* const direction_names[2] = {"x", "u"};
            for (int direction = 0; direction < 2; ++direction) {
                for (int x_region = 0; x_region < 3; ++x_region) {
                    for (int velocity_region = 0; velocity_region < 3;
                         ++velocity_region) {
                        const size_t bin = static_cast<size_t>(3 * x_region +
                                                               velocity_region);
                        const VlasovAmpereMidpointSolver::FctMacroBudget& budget =
                            (*direction_budgets[direction])[bin];
                        fct_macro_budget_monitor
                            << step << " " << time / Const::femto << " "
                            << 1 << " " << midpoint_result.state_advanced << " "
                            << midpoint_result.soft_unconverged << " "
                            << direction_names[direction] << " "
                            << x_region_names[x_region] << " "
                            << velocity_region_names[velocity_region] << " "
                            << budget.face_count << " "
                            << budget.active_face_count << " "
                            << budget.min_alpha << " "
                            << budget.raw_abs_mass << " "
                            << budget.rejected_abs_mass << " "
                            << (budget.raw_abs_mass > 0.0
                                ? budget.rejected_abs_mass / budget.raw_abs_mass
                                : 0.0) << " "
                            << budget.raw_abs_current << " "
                            << budget.rejected_abs_current << " "
                            << (budget.raw_abs_current > 0.0
                                ? budget.rejected_abs_current /
                                    budget.raw_abs_current : 0.0) << " "
                            << budget.raw_abs_energy << " "
                            << budget.rejected_abs_energy << " "
                            << (budget.raw_abs_energy > 0.0
                                ? budget.rejected_abs_energy /
                                    budget.raw_abs_energy : 0.0) << " "
                            << budget.delta_n << " " << budget.delta_j << " "
                            << budget.delta_k << " "
                            << budget.e_dot_j << " " << budget.r_fct_e << "\n";
                    }
                }
            }
            fct_macro_budget_monitor.flush();
        }
        if (mpi_rank == 0 && collect_step_diagnostics &&
            midpoint_result.state_advanced != 0 && !midpoint_result.failed) {
            const VlasovAmpereMidpointSolver::HighOrderCandidateAudit& audit =
                midpoint_result.high_order_candidate_audit;
            high_order_candidate_monitor
                << step << " " << time / Const::femto << " " << 1 << " "
                << midpoint_result.state_advanced << " "
                << midpoint_result.soft_unconverged << " "
                << audit.reconstructed_face_count << " "
                << audit.negative_reconstructed_face_count << " "
                << audit.negative_reconstructed_mass_weight << " "
                << audit.negative_reconstructed_current_weight << " "
                << audit.negative_reconstructed_energy_weight << " "
                << audit.high_candidate_negative_cell_count << " "
                << audit.high_candidate_negative_mass << " "
                << audit.high_candidate_negative_current_weight << " "
                << audit.high_candidate_negative_energy_weight << " "
                << audit.low_order_negative_cell_count << " "
                << audit.low_order_negative_mass << "\n";
            high_order_candidate_monitor.flush();
        }
        if (mpi_rank == 0 && collect_step_diagnostics &&
            midpoint_result.state_advanced != 0 && !midpoint_result.failed &&
            midpoint_result.background_coupling_mode ==
                VlasovAmpereMidpointSolver::DUAL_U_COUPLING) {
            const double dual_target_scale = std::max(
                std::numeric_limits<double>::min(),
                midpoint_result.final_dual_u_target_linf);
            const double dual_before_scale = std::max(
                std::numeric_limits<double>::min(),
                midpoint_result.final_dual_u_residual_before_linf);
            final_dual_u_pairing_monitor
                << step << " " << time / Const::femto << " "
                << 1 << " " << midpoint_result.state_advanced << " "
                << midpoint_result.soft_unconverged << " "
                << midpoint_result.final_dual_u_valid << " "
                << midpoint_result.final_dual_u_target_linf << " "
                << midpoint_result.final_dual_u_residual_before_linf << " "
                << midpoint_result.final_dual_u_residual_after_linf << " "
                << midpoint_result.final_dual_u_residual_after_linf /
                    dual_target_scale << " "
                << 1.0 -
                    midpoint_result.final_dual_u_residual_after_linf /
                    dual_before_scale << " "
                << midpoint_result.final_dual_u_minimum_scale << " "
                << midpoint_result.final_dual_u_correction_l2 << " "
                << midpoint_result.final_dual_u_correction_linf << " "
                << midpoint_result.final_dual_u_candidate_min << " "
                << midpoint_result.final_dual_u_corrected_cell_count << " "
                << midpoint_result.final_dual_u_limited_cell_count << " "
                << midpoint_result.final_dual_u_unresolved_cell_count << " "
                << midpoint_result.face_pairing_attempted << " "
                << midpoint_result.face_pairing_accepted << " "
                << midpoint_result.face_pairing_fallback_to_cell_baseline
                << " " << midpoint_result.face_pairing_solver_converged
                << " " << midpoint_result.face_pairing_iterations << " "
                << midpoint_result.face_pairing_unresolved_mode_count << " "
                << midpoint_result.face_pairing_residual_before << " "
                << midpoint_result.face_pairing_residual_after << " "
                << midpoint_result.face_pairing_core_residual_before << " "
                << midpoint_result.face_pairing_core_residual_after << " "
                << midpoint_result.face_pairing_unresolved_mode_l2 << " "
                << midpoint_result.face_pairing_correction_l2 << " "
                << midpoint_result.face_pairing_correction_linf << " "
                << midpoint_result.face_pairing_capacity_active_cells << " "
                << midpoint_result.face_pairing_trust_region_active_cells
                << " " << midpoint_result.face_pairing_delta_ke << " "
                << midpoint_result.face_pairing_delta_work << " "
                << midpoint_result.face_pairing_candidate_min << " "
                << midpoint_result.face_pairing_mass_error << " "
                << midpoint_result.face_pairing_candidate_valid << " "
                << midpoint_result.face_pairing_mass_relative_error << " "
                << midpoint_result.face_pairing_cell_mass_error_linf << " "
                << midpoint_result.face_pairing_cell_mass_relative_linf
                << " " << midpoint_result.face_pairing_energy_pair_error
                << " " << midpoint_result.face_pairing_energy_pair_relative
                << " " << midpoint_result.face_pairing_energy_residual_scale
                << " " << midpoint_result.face_pairing_energy_residual_ratio
                << " " << midpoint_result.face_pairing_correction_trust_limit
                << " " << midpoint_result.face_pairing_correction_trust_ratio
                << " " << midpoint_result.face_pairing_f_residual_relative_growth
                << " "
                << midpoint_result.face_pairing_candidate_residual_after
                << " "
                << midpoint_result.face_pairing_candidate_core_residual_after
                << " "
                << midpoint_result.face_pairing_candidate_delta_ke << " "
                << midpoint_result.face_pairing_candidate_delta_work << " "
                << midpoint_result.face_pairing_candidate_mass_error << " "
                << midpoint_result.face_pairing_candidate_min_before_fallback
                << " "
                << midpoint_result.face_pairing_requested_correction_l2
                << " "
                << midpoint_result.face_pairing_requested_correction_linf
                << " "
                << midpoint_result.face_pairing_applied_correction_l2 << " "
                << midpoint_result.face_pairing_applied_correction_linf
                << " "
                << midpoint_result.face_pairing_nonzero_capacity_cells
                << " "
                << midpoint_result.face_pairing_bound_saturated_cells << " "
                << midpoint_result.face_pairing_objective_residual << " "
                << midpoint_result.face_pairing_objective_smoothness << " "
                << midpoint_result.face_pairing_objective_amplitude << " "
                << midpoint_result.face_pairing_objective_total << " "
                << midpoint_result.face_pairing_rejection_mask << " "
                << midpoint_result.face_pairing_pass_solver << " "
                << midpoint_result.face_pairing_pass_apply << " "
                << midpoint_result.face_pairing_pass_global_residual << " "
                << midpoint_result.face_pairing_pass_core_residual << " "
                << midpoint_result.face_pairing_pass_correction_trust << " "
                << midpoint_result.face_pairing_pass_energy_pair << " "
                << midpoint_result.face_pairing_pass_energy_residual_scale
                << " "
                << midpoint_result.face_pairing_pass_candidate_min << " "
                << midpoint_result.face_pairing_pass_mass << " "
                << midpoint_result.face_pairing_pass_f_residual << "\n";
            final_dual_u_pairing_monitor.flush();
        }
        latest_bkg_energy_current_face =
            midpoint_result.j_bkg_energy_debug_face;
        latest_bkg_energy_current_valid =
            collect_step_diagnostics &&
            !latest_bkg_energy_current_face.empty();
        latest_bkg_ampere_current_face = midpoint_result.j_bkg_face_mid;
        // All current diagnostics returned by the midpoint solver are global
        // MPI values despite this legacy per-step storage name.
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
        coupled_iter_step = midpoint_result.nonlinear_iterations;
        coupled_residual_E_step = midpoint_result.residual_E;
        coupled_residual_J_bkg_step = midpoint_result.residual_J_bkg;
        coupled_residual_J_beam_step = midpoint_result.residual_J_beam;

        if (collect_energy_accounting) {
            // The midpoint solver already returns globally reduced current,
            // work, and kinetic-energy diagnostics.  Re-reducing these here
            // would multiply sums by the number of MPI ranks.
            const double global_dke_bkg_step = midpoint_result.delta_ke_bkg;
            const double global_W_bkg_E = midpoint_result.field_work_bkg;
            local_bkg_energy_residual_step =
                global_dke_bkg_step + global_W_bkg_E;
            bkg_energy_residual_step = local_bkg_energy_residual_step;
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
        if (collect_step_diagnostics) {
            FNegativitySnapshotDiagnostics f_neg_snapshot;
            compute_f_negativity_snapshot_diagnostics(
                bkg_e, sgrid, f_neg_snapshot);

            if (mpi_rank == 0) {
                // f-negativity monitor: final accepted snapshot at the
                // configured step-diagnostics cadence.
                std::ofstream f_neg_monitor;
                f_neg_monitor.open(output_path(runtime, "f_negativity_monitor.dat").c_str(),
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
        if (bkg_e.collisions_enabled && !split_collisions_applied) {
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
        if (split_collisions_applied) {
            // Both physical half steps applied their own collision update to
            // transaction-local states before the final commit above.
            cumulative_collision_energy_delta += collision_energy_step;
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

        net_nb_change_step = accepted_beam_injected_number -
            accepted_beam_outflow_number;

        if (beam_ledger_enabled) {
            // Beam counters are local because injection and outflow can be
            // owned by different ranks.  Summary mode accumulates locally and
            // performs one final reduction; only the explicit full audit mode
            // pays for a collective at every accepted step.
            double beam_ledger_local[6] = {
                accepted_beam_injected_number,
                accepted_beam_outflow_number,
                accepted_beam_injected_energy,
                accepted_beam_outflow_energy,
                accepted_beam_injected_impulse,
                accepted_beam_outflow_impulse
            };
            if (runtime.beam_ledger_mode == BEAM_LEDGER_SUMMARY) {
                for (int value = 0; value < 6; ++value) {
                    beam_ledger_local_totals[static_cast<size_t>(value)] +=
                        beam_ledger_local[value];
                }
            }
            if (write_beam_ledger) {
                double beam_ledger_global[12] = {0.0};
                const int ledger_rows = split_count_step != 0 ? 2 : 1;
                double beam_ledger_send[12] = {0.0};
                if (split_count_step != 0) {
                    for (int value = 0; value < 12; ++value)
                        beam_ledger_send[value] =
                            split_beam_ledger_rows[static_cast<size_t>(value)];
                } else {
                    for (int value = 0; value < 6; ++value)
                        beam_ledger_send[value] = beam_ledger_local[value];
                }
                MPI_Allreduce(beam_ledger_send, beam_ledger_global, 6 * ledger_rows,
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                if (mpi_rank == 0) {
                    for (int half = 0; half < ledger_rows; ++half) {
                        const size_t row = static_cast<size_t>(6 * half);
                        for (int value = 0; value < 6; ++value)
                            beam_ledger_full_totals[static_cast<size_t>(value)] +=
                                beam_ledger_global[row + static_cast<size_t>(value)];
                    }
                    if (runtime.beam_ledger_reference_enabled) {
                        for (int half = 0; half < ledger_rows &&
                             beam_reference_row_count < 2; ++half) {
                            const size_t destination = static_cast<size_t>(
                                6 * beam_reference_row_count);
                            const size_t source = static_cast<size_t>(6 * half);
                            for (int value = 0; value < 6; ++value)
                                beam_reference_rows[destination +
                                    static_cast<size_t>(value)] =
                                    beam_ledger_global[source +
                                        static_cast<size_t>(value)];
                            ++beam_reference_row_count;
                        }
                    }
                    for (int half = 0; half < ledger_rows; ++half) {
                        const size_t row = static_cast<size_t>(6 * half);
                        const double ledger_dt = split_count_step != 0 ?
                            0.5 * dt : dt;
                        const double ledger_time = split_count_step != 0 ?
                            current_time + (half + 1) * ledger_dt : time;
                        beam_half_step_ledger
                            << physical_step << " " << ledger_time / Const::femto
                            << " " << ledger_dt << " "
                            << beam_ledger_global[row] << " "
                            << beam_ledger_global[row + 1] << " "
                            << beam_ledger_global[row + 2] << " "
                            << beam_ledger_global[row + 3] << " "
                            << beam_ledger_global[row + 4] / ledger_dt << " "
                            << beam_ledger_global[row + 5] / ledger_dt << " "
                            << beam_ledger_global[row + 4] << " "
                            << beam_ledger_global[row + 5] << " " << half << "\n";
                    }
                    const bool checkpoint_due = runtime.checkpoint_enabled &&
                        next_checkpoint < runtime.checkpoint_times_fs.size() &&
                        time / Const::femto >= runtime.checkpoint_times_fs[next_checkpoint];
                    if (physical_step % 500 == 0 || time >= next_snapshot || checkpoint_due) {
                        flush_full_beam_ledger();
                    }
                }
            }
        }

        if (collect_energy_accounting) {
            const double bkg_ke_step_end = bkg_e.total_kinetic_energy();
            const double beam_ke_step_end = beam.total_kinetic_energy();
            const double field_energy_step_end = fields.total_energy();
            dke_bkg_step = bkg_ke_step_end - bkg_ke_step_start;
            dE_field_step = field_energy_step_end - field_energy_step_start;
            E_src_in_step = accepted_beam_injected_energy;
            E_src_out_step = accepted_beam_outflow_energy;
            const double total_energy_delta =
                (bkg_ke_step_end + beam_ke_step_end + field_energy_step_end) -
                (bkg_ke_step_start + beam_ke_step_start + field_energy_step_start);
            E_balance_step =
                total_energy_delta - E_src_in_step + E_src_out_step
                - collision_energy_step;

            // The quantities above are rank-local.  write_step_diagnostics()
            // performs its own packed reduction, but the accepted-energy
            // ledger is written directly by rank zero.  Reduce a separate
            // packed copy here so that its total-system account is global
            // without changing the established step-diagnostics protocol.
            double accepted_total_energy_terms[7] = {
                dke_bkg_step,
                beam_ke_step_end - beam_ke_step_start,
                dE_field_step,
                E_src_in_step,
                E_src_out_step,
                collision_energy_step,
                E_balance_step
            };
            if (collect_accepted_energy_audit) {
                MPI_Allreduce(MPI_IN_PLACE, accepted_total_energy_terms, 7,
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
            const double accepted_total_energy_residual =
                accepted_total_energy_terms[0] +
                accepted_total_energy_terms[1] +
                accepted_total_energy_terms[2] -
                accepted_total_energy_terms[3] +
                accepted_total_energy_terms[4] -
                accepted_total_energy_terms[5];
            const double accepted_total_energy_reconstruction_error =
                accepted_total_energy_residual -
                accepted_total_energy_terms[6];

            if (collect_step_diagnostics) {
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
                                         midpoint_result.stage5_jn_minus_je_linf,
                                         midpoint_result.stage5_r_fv,
                                         midpoint_result.stage5_r_couple,
                                         max_loss_u_high_step,
                                        x_at_max_loss_u_high_step,
                                        f_u_max_x_step,
                                        integral_f_u_gt_8_x_step,
                                        midpoint_result.background_coupling_mode);
            }
            if (collect_accepted_energy_audit &&
                midpoint_result.state_advanced != 0 && !midpoint_result.failed &&
                split_count_step == 0) {
                cumulative_total_energy_residual_for_audit +=
                    accepted_total_energy_residual;
                diag.write_accepted_energy_ledger(
                    physical_step, time, dt, midpoint_result,
                    accepted_total_energy_terms[0],
                    accepted_total_energy_terms[1],
                    accepted_total_energy_terms[2],
                    accepted_total_energy_terms[3],
                    accepted_total_energy_terms[4],
                    accepted_total_energy_terms[5],
                    accepted_total_energy_residual,
                    accepted_total_energy_reconstruction_error,
                    cumulative_total_energy_residual_for_audit,
                    mpi_rank);
            }
        }

        current_time = time;
        ++accepted_after_restart;
        performance_total_nonlinear_iterations +=
            midpoint_result.nonlinear_iterations;
        performance_total_operator_evaluations +=
            midpoint_result.operator_evaluations;
        performance_midpoint_predictor_used_steps +=
            midpoint_result.midpoint_predictor_used;
        performance_face_pairing_attempted_steps +=
            midpoint_result.face_pairing_attempted;
        performance_face_pairing_accepted_steps +=
            midpoint_result.face_pairing_accepted;
        performance_face_pairing_fallback_steps +=
            midpoint_result.face_pairing_fallback_to_cell_baseline;
        performance_retry_count += retry_count_step;
        performance_split_count += split_count_step;
        performance_operator_evaluations_retry +=
            operator_evaluations_retry_step;
        performance_wall_time_retry += wall_time_retry_step;
        performance_acceleration_attempts += midpoint_result.acceleration_attempts;
        performance_acceleration_accepted += midpoint_result.acceleration_accepted;
        performance_acceleration_fallback_evaluations +=
            midpoint_result.acceleration_fallback_evaluations;
        performance_acceleration_rejected_residual +=
            midpoint_result.acceleration_rejected_residual;
        performance_acceleration_rejected_nonfinite +=
            midpoint_result.acceleration_rejected_nonfinite;
        performance_acceleration_rejected_hard_failure +=
            midpoint_result.acceleration_rejected_hard_failure;
        performance_acceleration_rejected_coefficient +=
            midpoint_result.acceleration_rejected_coefficient;
        performance_acceleration_history_resets +=
            midpoint_result.acceleration_history_resets;
        performance_final_residual_e = midpoint_result.residual_E;
        performance_final_residual_j_bkg = midpoint_result.residual_J_bkg;
        performance_final_residual_f = midpoint_result.residual_f;
        performance_max_residual_e = std::max(performance_max_residual_e,
                                              midpoint_result.max_residual_E);
        performance_max_residual_j_bkg = std::max(performance_max_residual_j_bkg,
                                                  midpoint_result.max_residual_J_bkg);
        performance_max_residual_f = std::max(performance_max_residual_f,
                                              midpoint_result.max_residual_f);
        if (midpoint_result.soft_unconverged) {
            ++performance_soft_accepted_steps;
        } else {
            ++performance_strict_accepted_steps;
        }
        performance_limiter_active_sum +=
            midpoint_result.limiter_active_fraction;
        performance_limiter_active_min = std::min(
            performance_limiter_active_min,
            midpoint_result.limiter_active_fraction);
        performance_limiter_active_max = std::max(
            performance_limiter_active_max,
            midpoint_result.limiter_active_fraction);
        performance_limiter_min_alpha = std::min(
            performance_limiter_min_alpha,
            midpoint_result.limiter_min_alpha);
        performance_x_limiter_active_sum +=
            midpoint_result.x_limiter_active_fraction;
        performance_x_limiter_active_min = std::min(
            performance_x_limiter_active_min,
            midpoint_result.x_limiter_active_fraction);
        performance_x_limiter_active_max = std::max(
            performance_x_limiter_active_max,
            midpoint_result.x_limiter_active_fraction);
        performance_x_limiter_min_alpha = std::min(
            performance_x_limiter_min_alpha,
            midpoint_result.x_limiter_min_alpha);
        if (runtime.checkpoint_enabled && next_checkpoint < runtime.checkpoint_times_fs.size() &&
            current_time / Const::femto >= runtime.checkpoint_times_fs[next_checkpoint]) {
            flush_full_beam_ledger();
            std::ostringstream checkpoint_name_stream;
            checkpoint_name_stream << runtime.output_dir
                << "/checkpoints/t_" << std::fixed << std::setprecision(6)
                << std::setw(10) << std::setfill('0')
                << current_time / Const::femto
                << "fs_step_" << std::setw(8) << std::setfill('0')
                << physical_step;
            const std::string checkpoint_name = checkpoint_name_stream.str();
            CheckpointControlState control = { physical_step, current_time, dt, next_snapshot,
                                                last_snapshot_step, cumulative_collision_energy_delta };
            if (!write_checkpoint(checkpoint_name, control, bkg_e, beam, fields, sgrid,
                                  mpi_rank, mpi_size, runtime_error,
                                  midpoint_solver.low_order_only(),
                                  midpoint_solver.nonuniform_high_order_enabled(),
                                  midpoint_solver.fct_enabled())) {
                if (mpi_rank == 0) std::fprintf(stderr, "Checkpoint error: %s\n", runtime_error.c_str());
                flush_full_beam_ledger();
                MPI_Abort(MPI_COMM_WORLD, 12); return 12;
            }
            if (mpi_rank == 0)
                std::printf("Checkpoint written: %s\n",
                            checkpoint_name.c_str());
            ++next_checkpoint;
        }
        if (step % stdout_freq == 0) {
            diag.write_scalars(time, static_cast<int>(physical_step), bkg_e, beam, fields,
                               cumulative_collision_energy_delta,
                               mpi_rank, mpi_size);
            if (mpi_rank == 0) {
                printf("Step %lld, t = %.4f fs\n", physical_step, time / Const::femto);
            }
        }

        if (time >= next_snapshot) {
            write_snapshot(diag, time, bkg_e, beam, fields, ion_density_profile,
                           sgrid, mpi_rank, mpi_size,
                           config.enable_full_fe_output,
                           &latest_bkg_energy_current_face,
                           &latest_bkg_ampere_current_face,
                           latest_bkg_energy_current_valid);
            last_snapshot_step = static_cast<int>(physical_step);
            next_snapshot += Param::dt_snapshot;
        }
        if ((runtime.stop_after_accepted_steps >= 0 && accepted_after_restart >= runtime.stop_after_accepted_steps) ||
            current_time / Const::femto >= runtime.stop_time_fs) break;
    }

    sync_moments_and_charge(bkg_e, beam, fields, ion_density_profile,
                            moments_current);
    fields.update_gauss_residual_diagnostics(mpi_rank, mpi_size);
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    if (config.enable_debug_diagnostics) {
        diag.write_debug_state(static_cast<int>(restored_step + accepted_after_restart), current_time, "final", bkg_e, beam, fields,
                               sgrid, mpi_rank, mpi_size);
    }
#endif
    diag.write_scalars(current_time, static_cast<int>(restored_step + accepted_after_restart), bkg_e, beam, fields,
                       cumulative_collision_energy_delta,
                       mpi_rank, mpi_size);
    if (last_snapshot_step != restored_step + accepted_after_restart) {
        write_snapshot(diag, current_time, bkg_e, beam, fields, ion_density_profile,
                       sgrid, mpi_rank, mpi_size, config.enable_full_fe_output,
                       &latest_bkg_energy_current_face,
                       &latest_bkg_ampere_current_face,
                       latest_bkg_energy_current_valid);
    }

    const double performance_wall_local =
        MPI_Wtime() - performance_wall_start;
    double performance_wall_max = 0.0;
    MPI_Reduce(&performance_wall_local, &performance_wall_max, 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    unsigned long long local_max_rss_kib = 0ULL;
#if !defined(_WIN32)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        local_max_rss_kib =
            static_cast<unsigned long long>(usage.ru_maxrss);
    }
#endif
    unsigned long long global_max_rss_kib = 0ULL;
    MPI_Reduce(&local_max_rss_kib, &global_max_rss_kib, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    const CheckpointStateHashes final_state_hashes =
        checkpoint_state_hashes(bkg_e, beam, fields, sgrid,
                                mpi_rank, mpi_size);
    if (mpi_rank == 0) {
        const double accepted_scale =
            accepted_after_restart > 0
            ? static_cast<double>(accepted_after_restart) : 1.0;
        const double iteration_scale =
            performance_total_nonlinear_iterations > 0
            ? static_cast<double>(performance_total_nonlinear_iterations)
            : 1.0;
        const double operator_evaluation_scale =
            performance_total_operator_evaluations > 0
            ? static_cast<double>(performance_total_operator_evaluations)
            : 1.0;
        const double physical_time_fs = std::max(
            (current_time - performance_start_time_s) / Const::femto,
            std::numeric_limits<double>::min());
        std::ofstream performance_summary(
            output_path(runtime, "performance_summary.result").c_str());
        const char* const performance_beam_ledger_mode =
            runtime.beam_ledger_mode == BEAM_LEDGER_FULL ? "full" :
            runtime.beam_ledger_mode == BEAM_LEDGER_SUMMARY ? "summary" : "off";
        performance_summary << std::setprecision(17)
            << "status PASS\n"
            << "beam_ledger_mode " << performance_beam_ledger_mode << "\n"
            << "accepted_steps " << accepted_after_restart << "\n"
            << "total_nonlinear_iterations "
            << performance_total_nonlinear_iterations << "\n"
            << "total_operator_evaluations "
            << performance_total_operator_evaluations << "\n"
            << "mean_iterations_per_step "
            << performance_total_nonlinear_iterations / accepted_scale << "\n"
            << "mean_operator_evaluations_per_step "
            << performance_total_operator_evaluations / accepted_scale << "\n"
            << "retry_count " << performance_retry_count << "\n"
            << "split_count " << performance_split_count << "\n"
            << "wall_time_retry_s " << performance_wall_time_retry << "\n"
            << "operator_evaluations_retry "
            << performance_operator_evaluations_retry << "\n"
            << "wall_time_per_physical_fs_s "
            << performance_wall_max / physical_time_fs << "\n"
            << "soft_candidate_field_tolerance "
            << runtime.soft_candidate_field_tolerance << "\n"
            << "soft_candidate_current_tolerance "
            << runtime.soft_candidate_current_tolerance << "\n"
            << "soft_candidate_energy_p99_reference "
            << runtime.soft_candidate_energy_p99_reference << "\n"
            << "soft_candidate_energy_absolute_limit "
            << runtime.soft_candidate_energy_absolute_limit << "\n"
            << "midpoint_initial_guess_mode "
            << (runtime.midpoint_initial_guess_mode ==
                    RUNTIME_MIDPOINT_INITIAL_GUESS_FIELD_LINEAR
                ? "field-linear" : "none") << "\n"
            << "midpoint_predictor_used_steps "
            << performance_midpoint_predictor_used_steps << "\n"
            << "midpoint_predictor_use_fraction "
            << performance_midpoint_predictor_used_steps / accepted_scale
            << "\n"
            << "strict_accepted_steps "
            << performance_strict_accepted_steps << "\n"
            << "soft_accepted_steps "
            << performance_soft_accepted_steps << "\n"
            << "face_pairing_mode "
            << (runtime.face_pairing_mode == 1
                ? "regularized" : "cell-baseline") << "\n"
            << "face_pairing_sigma_cutoff "
            << runtime.face_pairing_sigma_cutoff << "\n"
            << "face_pairing_lambda " << runtime.face_pairing_lambda << "\n"
            << "face_pairing_eta " << runtime.face_pairing_eta << "\n"
            << "face_pairing_trust_fraction "
            << runtime.face_pairing_trust_fraction << "\n"
            << "face_pairing_correction_trust_fraction "
            << runtime.face_pairing_correction_trust_fraction << "\n"
            << "face_pairing_energy_pair_tolerance "
            << runtime.face_pairing_energy_pair_tolerance << "\n"
            << "face_pairing_energy_residual_fraction "
            << runtime.face_pairing_energy_residual_fraction << "\n"
            << "face_pairing_mass_relative_tolerance "
            << runtime.face_pairing_mass_relative_tolerance << "\n"
            << "face_pairing_f_residual_growth_tolerance "
            << runtime.face_pairing_f_residual_growth_tolerance << "\n"
            << "face_pairing_attempted_steps "
            << performance_face_pairing_attempted_steps << "\n"
            << "face_pairing_accepted_steps "
            << performance_face_pairing_accepted_steps << "\n"
            << "face_pairing_fallback_steps "
            << performance_face_pairing_fallback_steps << "\n"
            << "wall_seconds_internal " << performance_wall_max << "\n"
            << "wall_seconds_per_accepted_step "
            << performance_wall_max / accepted_scale << "\n"
            << "wall_seconds_per_nonlinear_iteration "
            << performance_wall_max / iteration_scale << "\n"
            << "wall_seconds_per_operator_evaluation "
            << performance_wall_max / operator_evaluation_scale << "\n"
            << "wall_seconds_per_physical_fs "
            << performance_wall_max / physical_time_fs << "\n"
            << "midpoint_acceleration_mode "
            << (runtime.midpoint_acceleration_mode == RUNTIME_MIDPOINT_ACCELERATION_AITKEN
                ? "aitken" : runtime.midpoint_acceleration_mode ==
                    RUNTIME_MIDPOINT_ACCELERATION_ANDERSON ? "anderson" : "none") << "\n"
            << "acceleration_attempts " << performance_acceleration_attempts << "\n"
            << "acceleration_accepted " << performance_acceleration_accepted << "\n"
            << "acceleration_fallback_evaluations "
            << performance_acceleration_fallback_evaluations << "\n"
            << "acceleration_rejected_residual "
            << performance_acceleration_rejected_residual << "\n"
            << "acceleration_rejected_nonfinite "
            << performance_acceleration_rejected_nonfinite << "\n"
            << "acceleration_rejected_hard_failure "
            << performance_acceleration_rejected_hard_failure << "\n"
            << "acceleration_rejected_coefficient "
            << performance_acceleration_rejected_coefficient << "\n"
            << "acceleration_history_resets " << performance_acceleration_history_resets << "\n"
            << "final_residual_E " << performance_final_residual_e << "\n"
            << "final_residual_J_bkg " << performance_final_residual_j_bkg << "\n"
            << "final_residual_f " << performance_final_residual_f << "\n"
            << "max_residual_E " << performance_max_residual_e << "\n"
            << "max_residual_J_bkg " << performance_max_residual_j_bkg << "\n"
            << "max_residual_f " << performance_max_residual_f << "\n"
            << "limiter_active_fraction_mean "
            << performance_limiter_active_sum / accepted_scale << "\n"
            << "limiter_active_fraction_min "
            << (accepted_after_restart > 0
                ? performance_limiter_active_min : 0.0) << "\n"
            << "limiter_active_fraction_max "
            << performance_limiter_active_max << "\n"
            << "limiter_min_alpha "
            << performance_limiter_min_alpha << "\n"
            << "x_limiter_active_fraction_mean "
            << performance_x_limiter_active_sum / accepted_scale << "\n"
            << "x_limiter_active_fraction_min "
            << (accepted_after_restart > 0
                ? performance_x_limiter_active_min : 0.0) << "\n"
            << "x_limiter_active_fraction_max "
            << performance_x_limiter_active_max << "\n"
            << "x_limiter_min_alpha "
            << performance_x_limiter_min_alpha << "\n"
            << "max_rss_per_rank_kib " << global_max_rss_kib << "\n"
            << "max_rss_available "
            << (global_max_rss_kib > 0ULL ? 1 : 0) << "\n"
            << "energy_diagnostic_current_valid "
            << (latest_bkg_energy_current_valid ? 1 : 0) << "\n"
            << "final_state_hash_background "
            << final_state_hashes.background << "\n"
            << "final_state_hash_field_faces "
            << final_state_hashes.field_faces << "\n"
            << "final_state_hash_beam "
            << final_state_hashes.beam << "\n";
    }

    std::array<double, 6> beam_ledger_global_totals = {{0.0, 0.0, 0.0,
                                                          0.0, 0.0, 0.0}};
    if (runtime.beam_ledger_mode == BEAM_LEDGER_SUMMARY) {
        MPI_Reduce(beam_ledger_local_totals.data(), beam_ledger_global_totals.data(), 6,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else if (runtime.beam_ledger_mode == BEAM_LEDGER_FULL && mpi_rank == 0) {
        beam_ledger_global_totals = beam_ledger_full_totals;
    }

    if (mpi_rank == 0) {
        if (write_beam_ledger) {
            flush_full_beam_ledger();
            beam_half_step_ledger.close();
        }
        if (beam_ledger_enabled) {
            const char* const mode_name = runtime.beam_ledger_mode == BEAM_LEDGER_FULL
                ? "full" : "summary";
            std::ofstream ledger_summary(
                output_path(runtime, "beam_ledger_summary.result").c_str());
            ledger_summary << std::setprecision(17)
                << "beam_ledger_mode " << mode_name << "\n"
                << "accepted_substeps " << accepted_after_restart << "\n"
                << "start_time_s " << beam_ledger_start_time_s << "\n"
                << "end_time_s " << current_time << "\n"
                << "dt_s " << dt << "\n"
                << "N_in_total " << beam_ledger_global_totals[0] << "\n"
                << "N_out_total " << beam_ledger_global_totals[1] << "\n"
                << "injected_energy_total " << beam_ledger_global_totals[2] << "\n"
                << "outflow_energy_total " << beam_ledger_global_totals[3] << "\n"
                << "injected_current_impulse_total " << beam_ledger_global_totals[4] << "\n"
                << "outflow_current_impulse_total " << beam_ledger_global_totals[5] << "\n";
        }
        if (runtime.beam_ledger_reference_enabled) {
            double reference_n_in = 0.0;
            double reference_current_impulse = 0.0;
            const bool have_reference = read_beam_ledger_total(
                runtime.beam_ledger_reference, reference_n_in, reference_current_impulse);
            const double n_in_difference = have_reference
                ? std::fabs(reference_n_in - beam_ledger_global_totals[0]) : 0.0;
            const double current_impulse_difference = have_reference
                ? std::fabs(reference_current_impulse - beam_ledger_global_totals[4]) : 0.0;
            const double n_in_scale = std::max(1.0, std::max(
                std::fabs(reference_n_in), std::fabs(beam_ledger_global_totals[0])));
            const double current_impulse_scale = std::max(1.0, std::max(
                std::fabs(reference_current_impulse), std::fabs(beam_ledger_global_totals[4])));
            const bool reference_comparison_pass = !have_reference ||
                (n_in_difference <= 4096.0 * std::numeric_limits<double>::epsilon() *
                 n_in_scale && current_impulse_difference <=
                 4096.0 * std::numeric_limits<double>::epsilon() * current_impulse_scale);
            const bool two_substep_ledger_complete = beam_reference_row_count == 2 &&
                accepted_after_restart == 2;
            std::ofstream ledger_result(
                output_path(runtime, "beam_dt_half_time_consistency.result").c_str());
            ledger_result << std::setprecision(17)
                << "accepted_substeps " << accepted_after_restart << "\n"
                << "dt_scale " << runtime.dt_scale << "\n"
                << "two_substep_ledger_complete "
                << (two_substep_ledger_complete ? 1 : 0) << "\n"
                << "N_in_total " << beam_ledger_global_totals[0] << "\n"
                << "N_out_total " << beam_ledger_global_totals[1] << "\n"
                << "injected_energy_total " << beam_ledger_global_totals[2] << "\n"
                << "outflow_energy_total " << beam_ledger_global_totals[3] << "\n"
                << "injected_current_impulse_total " << beam_ledger_global_totals[4] << "\n"
                << "outflow_current_impulse_total " << beam_ledger_global_totals[5] << "\n"
                << "reference_available " << (have_reference ? 1 : 0) << "\n"
                << "reference_comparison_pass " << (reference_comparison_pass ? 1 : 0) << "\n";
            for (int substep = 0; substep < beam_reference_row_count; ++substep) {
                const size_t base = static_cast<size_t>(6 * substep);
                ledger_result << "substep_" << (substep + 1) << "_N_in "
                    << beam_reference_rows[base] << "\n"
                    << "substep_" << (substep + 1) << "_N_out "
                    << beam_reference_rows[base + 1] << "\n"
                    << "substep_" << (substep + 1) << "_injected_energy "
                    << beam_reference_rows[base + 2] << "\n"
                    << "substep_" << (substep + 1) << "_outflow_energy "
                    << beam_reference_rows[base + 3] << "\n"
                    << "substep_" << (substep + 1) << "_injected_current_impulse "
                    << beam_reference_rows[base + 4] << "\n"
                    << "substep_" << (substep + 1) << "_outflow_current_impulse "
                    << beam_reference_rows[base + 5] << "\n";
            }
            if (have_reference) {
                ledger_result << "reference_N_in_total " << reference_n_in << "\n"
                    << "N_in_A_minus_D_sum " << n_in_difference << "\n"
                    << "reference_injected_current_impulse_total "
                    << reference_current_impulse << "\n"
                    << "injected_current_impulse_A_minus_D_sum "
                    << current_impulse_difference << "\n";
            }
            ledger_result.close();
            std::ifstream consistency_input(
                output_path(runtime, "beam_dt_half_time_consistency.result").c_str());
            std::ofstream legacy_summary(
                output_path(runtime, "beam_dt_half_time_summary.result").c_str());
            legacy_summary << consistency_input.rdbuf();
        }
        printf("============================================================\n");
        printf("  Simulation complete: t = %.1f fs, %d steps\n",
               current_time / Const::femto, static_cast<int>(restored_step + accepted_after_restart));
        printf("============================================================\n");
    }

    MPI_Finalize();
    return 0;
}
