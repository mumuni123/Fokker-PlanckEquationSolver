#include "background_coupling_test_support.h"
#include "discrete_moment_operators.h"
#include "periodic_staggered_operators.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#endif

#ifndef FP_DUAL_U_STAGE2_CASE
#define FP_DUAL_U_STAGE2_CASE 1
#endif
#ifndef FP_DUAL_U_STAGE2_GIT_REVISION
#define FP_DUAL_U_STAGE2_GIT_REVISION "unknown"
#endif
#ifndef FP_DUAL_U_STAGE2_GIT_DIRTY
#define FP_DUAL_U_STAGE2_GIT_DIRTY -1
#endif

namespace {

enum CaseKind {
    ZERO_FIELD = 1,
    FIELD_SYMMETRY = 2,
    PHASE_TRANSLATION = 3,
    NONUNIFORM_GRID = 4,
    MANUFACTURED = 5,
    ENDPOINT_FLUX = 6,
    MPI_CONSISTENCY = 7
};

struct PairMetrics {
    VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle legacy;
    VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle dual;
    BackgroundCouplingTest::Norms legacy_pair;
    BackgroundCouplingTest::Norms dual_pair;
    double legacy_signed;
    double dual_signed;
    double jn_difference;
    double je_difference;
    double state_difference_l1;
    BackgroundCouplingTest::Norms dual_cell_target;
    BackgroundCouplingTest::Norms dual_cell_target_scale;
    BackgroundCouplingTest::Norms projection_floor;
    BackgroundCouplingTest::Norms projection_reconstruction;
    BackgroundCouplingTest::Norms projection_target_scale;
    double dual_cell_target_signed;
    double projection_floor_signed;
    double projection_reconstruction_signed;
    int valid;
};

struct PhaseMetrics {
    double jn_relative_l2;
    double gstar_je_relative_l2;
    double residual_relative_l2;
    double state_relative_l2;

    PhaseMetrics()
        : jn_relative_l2(0.0), gstar_je_relative_l2(0.0),
          residual_relative_l2(0.0), state_relative_l2(0.0) {}

    double maximum() const
    {
        return std::max(std::max(jn_relative_l2, gstar_je_relative_l2),
                        std::max(residual_relative_l2, state_relative_l2));
    }
};

struct FieldSymmetryMetrics {
    double momentum_odd_absolute;
    double momentum_response_scale;
    double momentum_roundoff_scale;
    double momentum_odd_relative;
    double state_odd_relative;
    double u_flux_odd_relative;
    double jn_even_relative;
    double je_even_relative;

    FieldSymmetryMetrics()
        : momentum_odd_absolute(0.0), momentum_response_scale(0.0),
          momentum_roundoff_scale(0.0), momentum_odd_relative(0.0),
          state_odd_relative(0.0),
          u_flux_odd_relative(0.0), jn_even_relative(0.0),
          je_even_relative(0.0) {}

    double maximum() const
    {
        return std::max(std::max(momentum_odd_relative, state_odd_relative),
            std::max(u_flux_odd_relative,
                     std::max(jn_even_relative, je_even_relative)));
    }
};

struct ZeroFieldMetrics {
    double state_relative_l2;
    double u_flux_linf;
    double jn_linf;
    double cancellation_current_scale;
    double jn_relative_to_cancellation_scale;

    ZeroFieldMetrics()
        : state_relative_l2(0.0), u_flux_linf(0.0), jn_linf(0.0),
          cancellation_current_scale(0.0),
          jn_relative_to_cancellation_scale(0.0) {}
};

struct EndpointBalanceMetrics {
    double delta_kinetic_energy;
    double field_work;
    double boundary_work;
    double spatial_boundary_work;
    double residual;
    double scale;
    double relative;
    double r_fv_match_relative;
    double production_delta_kinetic_energy;
    double production_u_energy_moment;
    double production_r_fv_reconstructed;
    double production_r_fv_reconstruction_relative;
    double delta_kinetic_energy_match_relative;
    double u_energy_field_work_difference_absolute;
    double u_energy_field_work_scale;
    double u_energy_field_work_match_relative;
    double u_flux_coefficient_relation_linf;
    double u_flux_coefficient_relation_scale;
    double u_flux_coefficient_relation_relative;
    int final_coefficient_available;
    double reconstructed_u_energy_moment;
    double reconstructed_field_work;
    double production_u_energy_reconstruction_relative;
    double final_current_work_reconstruction_relative;
    double u_work_roundoff_bound;
    int u_work_roundoff_limited;
    double final_state_energy_gross_scale;
    int local_update_replay_available;
    double local_update_energy_defect_signed;
    double local_update_energy_defect_absolute;
    double local_update_energy_ulp_bound;
    double balance_roundoff_bound;
    int balance_roundoff_limited;
    double r_fv_match_difference_absolute;
    double r_fv_match_roundoff_bound;
    int r_fv_match_roundoff_limited;
    double spatial_boundary_relative;

    EndpointBalanceMetrics()
        : delta_kinetic_energy(0.0), field_work(0.0),
          boundary_work(0.0), spatial_boundary_work(0.0), residual(0.0),
          scale(0.0), relative(0.0), r_fv_match_relative(0.0),
          production_delta_kinetic_energy(0.0),
          production_u_energy_moment(0.0),
          production_r_fv_reconstructed(0.0),
          production_r_fv_reconstruction_relative(0.0),
          delta_kinetic_energy_match_relative(0.0),
          u_energy_field_work_difference_absolute(0.0),
          u_energy_field_work_scale(0.0),
          u_energy_field_work_match_relative(0.0),
          u_flux_coefficient_relation_linf(0.0),
          u_flux_coefficient_relation_scale(0.0),
          u_flux_coefficient_relation_relative(0.0),
          final_coefficient_available(0),
          reconstructed_u_energy_moment(0.0),
          reconstructed_field_work(0.0),
          production_u_energy_reconstruction_relative(0.0),
          final_current_work_reconstruction_relative(0.0),
          u_work_roundoff_bound(0.0), u_work_roundoff_limited(0),
          final_state_energy_gross_scale(0.0),
          local_update_replay_available(0),
          local_update_energy_defect_signed(0.0),
          local_update_energy_defect_absolute(0.0),
          local_update_energy_ulp_bound(0.0),
          balance_roundoff_bound(0.0), balance_roundoff_limited(0),
          r_fv_match_difference_absolute(0.0),
          r_fv_match_roundoff_bound(0.0),
          r_fv_match_roundoff_limited(0),
          spatial_boundary_relative(0.0) {}
};

const char* case_name()
{
    switch (FP_DUAL_U_STAGE2_CASE) {
    case ZERO_FIELD: return "zero_field";
    case FIELD_SYMMETRY: return "field_symmetry";
    case PHASE_TRANSLATION: return "phase_translation";
    case NONUNIFORM_GRID: return "nonuniform_grid";
    case MANUFACTURED: return "manufactured";
    case ENDPOINT_FLUX: return "endpoint_flux";
    case MPI_CONSISTENCY: return "mpi_consistency";
    default: return "unknown";
    }
}

bool make_directory_tree(const std::string& path)
{
    if (path.empty()) return false;
    std::string current;
    for (size_t i = 0; i < path.size(); ++i) {
        current.push_back(path[i]);
        const bool boundary = path[i] == '/' || path[i] == '\\' ||
                              i + 1 == path.size();
        if (!boundary) continue;
        while (!current.empty() &&
               (current[current.size() - 1] == '/' ||
                current[current.size() - 1] == '\\'))
            current.resize(current.size() - 1);
        if (current.empty()) continue;
#ifdef _WIN32
        const int created = _mkdir(current.c_str());
#else
        const int created = mkdir(current.c_str(), 0777);
#endif
        if (created != 0 && errno != EEXIST)
            return false;
        if (i + 1 < path.size()) current.push_back('/');
    }
    return true;
}

std::string output_directory(int argc, char** argv)
{
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--output-dir") return argv[i + 1];
    return std::string("output/dual_u_stage2/") + case_name();
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

BackgroundCouplingTest::Norms local_vector_difference_norms(
    const std::vector<double>& left, const std::vector<double>& right,
    int count, const SpatialGrid& sg)
{
    double local[3] = {0.0, 0.0, 0.0};
    if (left.size() < static_cast<size_t>(count) ||
        right.size() < static_cast<size_t>(count)) {
        local[0] = local[1] = local[2] =
            std::numeric_limits<double>::infinity();
    } else {
        for (int i = 0; i < count; ++i) {
            const double difference = left[static_cast<size_t>(i)] -
                                      right[static_cast<size_t>(i)];
            local[0] += std::fabs(difference);
            local[1] += difference * difference;
            local[2] = std::max(local[2], std::fabs(difference));
        }
    }
    double global[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local + 2, global + 2, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    BackgroundCouplingTest::Norms result = {};
    result.l1 = global[0] / static_cast<double>(sg.nx_global);
    result.l2 = std::sqrt(global[1] / static_cast<double>(sg.nx_global));
    result.linf = global[2];
    return result;
}

double local_vector_signed_difference(const std::vector<double>& left,
                                      const std::vector<double>& right,
                                      int count, const SpatialGrid& sg)
{
    double local = 0.0;
    for (int i = 0; i < count; ++i)
        local += (left[static_cast<size_t>(i)] -
                  right[static_cast<size_t>(i)]) * sg.dx;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

double global_linear_combination_relative_l2(
    const std::vector<double>& a, double ca,
    const std::vector<double>& b, double cb,
    const std::vector<double>& c, double cc)
{
    if (a.size() != b.size() || a.size() != c.size())
        return std::numeric_limits<double>::infinity();
    double local_error = 0.0;
    double local_scale_a = 0.0;
    double local_scale_b = 0.0;
    double local_scale_c = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double error = ca * a[i] + cb * b[i] + cc * c[i];
        local_error += error * error;
        local_scale_a += (ca * a[i]) * (ca * a[i]);
        local_scale_b += (cb * b[i]) * (cb * b[i]);
        local_scale_c += (cc * c[i]) * (cc * c[i]);
    }
    double sums[4] = {local_error, local_scale_a, local_scale_b,
                      local_scale_c};
    MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double scale = std::sqrt(std::max(sums[1],
                                  std::max(sums[2], sums[3])));
    if (scale == 0.0) return sums[0] == 0.0 ? 0.0 :
        std::numeric_limits<double>::infinity();
    return std::sqrt(sums[0]) / scale;
}

double global_odd_about_reference_relative_l2(
    const std::vector<double>& positive,
    const std::vector<double>& negative,
    const std::vector<double>& reference)
{
    if (positive.size() != negative.size() ||
        positive.size() != reference.size())
        return std::numeric_limits<double>::infinity();
    double local[3] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < positive.size(); ++i) {
        const double plus_response = positive[i] - reference[i];
        const double minus_response = negative[i] - reference[i];
        const double error = plus_response + minus_response;
        local[0] += error * error;
        local[1] += plus_response * plus_response;
        local[2] += minus_response * minus_response;
    }
    MPI_Allreduce(MPI_IN_PLACE, local, 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    const double scale = std::sqrt(std::max(local[1], local[2]));
    if (scale == 0.0) return local[0] == 0.0 ? 0.0 :
        std::numeric_limits<double>::infinity();
    return std::sqrt(local[0]) / scale;
}

std::vector<double> gather_rank_ordered(const std::vector<double>& local,
                                        int rank, int size)
{
    const int local_count = static_cast<int>(local.size());
    std::vector<int> counts(rank == 0 ? static_cast<size_t>(size) : 0, 0);
    MPI_Gather(&local_count, 1, MPI_INT,
               rank == 0 ? counts.data() : 0, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    std::vector<int> offsets(rank == 0 ? static_cast<size_t>(size) : 0, 0);
    std::vector<double> global;
    if (rank == 0) {
        int total = 0;
        for (int r = 0; r < size; ++r) {
            offsets[static_cast<size_t>(r)] = total;
            total += counts[static_cast<size_t>(r)];
        }
        global.assign(static_cast<size_t>(total), 0.0);
    }
    MPI_Gatherv(local.empty() ? 0 : local.data(), local_count, MPI_DOUBLE,
        rank == 0 ? global.data() : 0,
        rank == 0 ? counts.data() : 0,
        rank == 0 ? offsets.data() : 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return global;
}

double shifted_relative_l2(const std::vector<double>& base,
                           const std::vector<double>& shifted,
                           int nx, int block, int shift_cells,
                           int rank)
{
    double relative = 0.0;
    if (rank == 0) {
        if (base.size() != static_cast<size_t>(nx) * block ||
            shifted.size() != base.size()) {
            relative = std::numeric_limits<double>::infinity();
        } else {
            double error = 0.0;
            double base_scale = 0.0;
            double shifted_scale = 0.0;
            for (int ix = 0; ix < nx; ++ix) {
                const int shifted_ix = (ix - shift_cells + nx) % nx;
                for (int q = 0; q < block; ++q) {
                    const double lhs = base[static_cast<size_t>(ix) * block + q];
                    const double rhs = shifted[
                        static_cast<size_t>(shifted_ix) * block + q];
                    const double difference = lhs - rhs;
                    error += difference * difference;
                    base_scale += lhs * lhs;
                    shifted_scale += rhs * rhs;
                }
            }
            const double scale = std::sqrt(std::max(base_scale, shifted_scale));
            relative = scale == 0.0 ? (error == 0.0 ? 0.0 :
                std::numeric_limits<double>::infinity()) :
                std::sqrt(error) / scale;
        }
    }
    MPI_Bcast(&relative, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return relative;
}

double signed_pair(const std::vector<double>& left,
                   const std::vector<double>& right,
                   const SpatialGrid& sg)
{
    double local = 0.0;
    if (left.size() >= static_cast<size_t>(sg.nx_local) &&
        right.size() >= static_cast<size_t>(sg.nx_local))
        for (int i = 0; i < sg.nx_local; ++i)
            local += (left[static_cast<size_t>(i)] -
                      right[static_cast<size_t>(i)]) * sg.dx;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

void set_uniform_field(EMFields& fields, const SpatialGrid& sg,
                       int rank, int size, double value)
{
    for (int iface = 0; iface < sg.nx_local; ++iface)
        fields.Ex_face[static_cast<size_t>(iface)] = value;
    fields.sync_cell_ex_from_faces(rank, size);
}

void add_endpoint_population(Species& background, const SpatialGrid& sg,
                             int rank, int size)
{
    const int ng = sg.nghost;
    const int center = Param::Nv / 2;
    for (int ix = 0; ix < sg.nx_local; ++ix)
        for (int k = 0; k < Param::Nmu; ++k) {
            const double reference = background.f[idx3(ng + ix, center, k)];
            background.f[idx3(ng + ix, 0, k)] += 1.0e-8 * reference;
            background.f[idx3(ng + ix, Param::Nv - 1, k)] +=
                1.0e-8 * reference;
        }
    VlasovAmpereMidpointSolver sync;
    sync.synchronize_background_ghosts(background, sg, rank, size);
}

PairMetrics evaluate_pair(const Species& background, const EMFields& fields,
                          const SpatialGrid& sg, int rank, int size,
                          double dt)
{
    BackgroundCouplingTest::BundleOptions legacy_options;
    legacy_options.fct_enabled = false;
    legacy_options.allow_finite_negative_debt = true;
    legacy_options.coupling_mode =
        VlasovAmpereMidpointSolver::LEGACY_COUPLING;
    PairMetrics result;
    result.legacy = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, legacy_options);
    BackgroundCouplingTest::BundleOptions dual_options = legacy_options;
    dual_options.coupling_mode =
        VlasovAmpereMidpointSolver::DUAL_U_COUPLING;
    result.dual = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, dual_options);

    const size_t face_count = static_cast<size_t>(sg.nx_local + 1);
    const bool local_valid = result.legacy.state_advanced &&
        result.dual.state_advanced && !result.legacy.operator_failed &&
        !result.dual.operator_failed && result.legacy.outputs_finite &&
        result.dual.outputs_finite && result.dual.dual_u_operator_valid &&
        result.legacy.jn_high.size() == face_count &&
        result.dual.jn_high.size() == face_count &&
        result.legacy.gstar_je_high.size() == face_count &&
        result.dual.gstar_je_high.size() == face_count &&
        result.dual.dual_target_jn_cell.size() ==
            static_cast<size_t>(sg.nx_local) &&
        result.dual.dual_je_cell.size() ==
            static_cast<size_t>(sg.nx_local) &&
        result.dual.je_final.size() ==
            static_cast<size_t>(sg.nx_local) &&
        result.legacy.final_state_mass.size() ==
            static_cast<size_t>(sg.nx_local) * Param::Nvmu &&
        result.dual.final_state_mass.size() ==
            static_cast<size_t>(sg.nx_local) * Param::Nvmu;
    result.valid = local_valid ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &result.valid, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    result.legacy_pair = {0.0, 0.0, 0.0};
    result.dual_pair = {0.0, 0.0, 0.0};
    result.legacy_signed = 0.0;
    result.dual_signed = 0.0;
    result.jn_difference = std::numeric_limits<double>::infinity();
    result.je_difference = std::numeric_limits<double>::infinity();
    result.state_difference_l1 = std::numeric_limits<double>::infinity();
    result.dual_cell_target = {0.0, 0.0, 0.0};
    result.dual_cell_target_scale = {0.0, 0.0, 0.0};
    result.projection_floor = {0.0, 0.0, 0.0};
    result.projection_reconstruction = {0.0, 0.0, 0.0};
    result.projection_target_scale = {0.0, 0.0, 0.0};
    result.dual_cell_target_signed = 0.0;
    result.projection_floor_signed = 0.0;
    result.projection_reconstruction_signed = 0.0;
    if (!result.valid) return result;
    result.legacy_pair = BackgroundCouplingTest::face_difference_norms(
        result.legacy.jn_high, result.legacy.gstar_je_high, sg);
    result.dual_pair = BackgroundCouplingTest::face_difference_norms(
        result.dual.jn_high, result.dual.gstar_je_high, sg);
    result.legacy_signed = signed_pair(result.legacy.jn_high,
        result.legacy.gstar_je_high, sg);
    result.dual_signed = signed_pair(result.dual.jn_high,
        result.dual.gstar_je_high, sg);
    result.jn_difference = global_max_difference(result.legacy.jn_high,
                                                  result.dual.jn_high);
    result.je_difference = global_max_difference(
        result.legacy.gstar_je_high, result.dual.gstar_je_high);
    result.dual_cell_target = local_vector_difference_norms(
        result.dual.dual_je_cell, result.dual.dual_target_jn_cell,
        sg.nx_local, sg);
    result.dual_cell_target_scale = local_vector_difference_norms(
        result.dual.dual_target_jn_cell,
        std::vector<double>(static_cast<size_t>(sg.nx_local), 0.0),
        sg.nx_local, sg);
    result.dual_cell_target_signed = local_vector_signed_difference(
        result.dual.dual_je_cell, result.dual.dual_target_jn_cell,
        sg.nx_local, sg);
    std::vector<double> projected_target;
    PeriodicStaggered::apply_cell_to_face_Gstar(
        result.dual.dual_target_jn_cell, projected_target, sg.nx_local,
        rank, size, 48721);
    result.projection_floor = BackgroundCouplingTest::face_difference_norms(
        result.dual.jn_high, projected_target, sg);
    result.projection_target_scale =
        BackgroundCouplingTest::face_difference_norms(
            projected_target,
            std::vector<double>(projected_target.size(), 0.0), sg);
    result.projection_floor_signed = signed_pair(
        result.dual.jn_high, projected_target, sg);
    result.projection_reconstruction =
        BackgroundCouplingTest::face_difference_norms(
            result.dual.gstar_je_high, projected_target, sg);
    result.projection_reconstruction_signed = signed_pair(
        result.dual.gstar_je_high, projected_target, sg);
    double local_state_l1 = 0.0;
    if (result.legacy.final_state_mass.size() ==
        result.dual.final_state_mass.size())
        for (size_t i = 0; i < result.legacy.final_state_mass.size(); ++i)
            local_state_l1 += std::fabs(result.legacy.final_state_mass[i] -
                                        result.dual.final_state_mass[i]);
    else
        local_state_l1 = std::numeric_limits<double>::infinity();
    MPI_Allreduce(&local_state_l1, &result.state_difference_l1, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

unsigned long long grid_hash(const std::vector<double>& values)
{
    unsigned long long hash = 1469598103934665603ULL;
    for (size_t i = 0; i < values.size(); ++i) {
        unsigned long long bits = 0;
        std::memcpy(&bits, &values[i], sizeof(bits));
        hash ^= bits;
        hash *= 1099511628211ULL;
    }
    return hash;
}

double relative_error(double value, double scale)
{
    return std::fabs(value) / std::max(1.0, std::fabs(scale));
}

double relative_to_nonzero_scale(double value, double scale)
{
    if (scale == 0.0)
        return value == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    return std::fabs(value) / std::fabs(scale);
}

double rank_ordered_sum(double local, int rank, int size)
{
    std::vector<double> contributions(
        rank == 0 ? static_cast<size_t>(size) : 0, 0.0);
    MPI_Gather(&local, 1, MPI_DOUBLE,
               rank == 0 ? contributions.data() : 0, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    double result = 0.0;
    if (rank == 0) {
        long double ordered = 0.0L;
        for (int r = 0; r < size; ++r)
            ordered += static_cast<long double>(
                contributions[static_cast<size_t>(r)]);
        result = static_cast<double>(ordered);
    }
    MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return result;
}

long double rank_ordered_sum_long_double(long double local, int rank, int size)
{
    std::vector<long double> contributions(
        rank == 0 ? static_cast<size_t>(size) : 0, 0.0L);
    MPI_Gather(&local, 1, MPI_LONG_DOUBLE,
               rank == 0 ? contributions.data() : 0, 1, MPI_LONG_DOUBLE, 0,
               MPI_COMM_WORLD);
    long double result = 0.0L;
    if (rank == 0) {
        for (int r = 0; r < size; ++r)
            result += contributions[static_cast<size_t>(r)];
    }
    MPI_Bcast(&result, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    return result;
}

std::vector<double> physical_state(const Species& background,
                                   const SpatialGrid& sg)
{
    std::vector<double> state(static_cast<size_t>(sg.nx_local) * Param::Nvmu,
                              0.0);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const size_t source = static_cast<size_t>(sg.nghost + ix) * Param::Nvmu;
        const size_t target = static_cast<size_t>(ix) * Param::Nvmu;
        std::copy(background.f.begin() + source,
                  background.f.begin() + source + Param::Nvmu,
                  state.begin() + target);
    }
    return state;
}

double max_abs_owned_face(const std::vector<double>& values,
                          const SpatialGrid& sg)
{
    double local = 0.0;
    for (int iface = 0; iface < sg.nx_local; ++iface)
        local = std::max(local, std::fabs(values[static_cast<size_t>(iface)]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double cancellation_current_scale(const Species& background,
                                  const SpatialGrid& sg)
{
    double local = 0.0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        double scale = 0.0;
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                scale += std::fabs(background.charge *
                    background.cgrid.vx[idx2(j, k)] *
                    background.f[idx3(sg.nghost + ix, j, k)] / sg.dx);
        local = std::max(local, scale);
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

FieldSymmetryMetrics deterministic_field_symmetry(
    const PairMetrics& positive, const PairMetrics& negative,
    const PairMetrics& zero_reference, const Species& background,
    const SpatialGrid& sg, int rank, int size)
{
    FieldSymmetryMetrics metrics;
    long double local_plus = 0.0L;
    long double local_minus = 0.0L;
    long double local_odd = 0.0L;
    long double local_roundoff = 0.0L;
    for (int ix = 0; ix < sg.nx_local; ++ix)
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = (static_cast<size_t>(ix) * Param::Nv + j) *
                                  Param::Nmu + k;
                const long double weight = static_cast<long double>(
                    background.mass * Const::c *
                    background.cgrid.upar_cells[j]);
                const long double plus_response = weight *
                    static_cast<long double>(positive.dual.final_state_mass[id] -
                        zero_reference.dual.final_state_mass[id]);
                const long double minus_response = weight *
                    static_cast<long double>(negative.dual.final_state_mass[id] -
                        zero_reference.dual.final_state_mass[id]);
                local_plus += plus_response;
                local_minus += minus_response;
                local_odd += plus_response + minus_response;
                local_roundoff += std::fabs(plus_response) +
                                  std::fabs(minus_response);
            }
    const double plus = rank_ordered_sum(static_cast<double>(local_plus), rank,
                                         size);
    const double minus = rank_ordered_sum(static_cast<double>(local_minus), rank,
                                          size);
    const double odd = rank_ordered_sum(static_cast<double>(local_odd), rank,
                                        size);
    const double roundoff = rank_ordered_sum(
        static_cast<double>(local_roundoff), rank, size);
    metrics.momentum_odd_absolute = std::fabs(odd);
    metrics.momentum_response_scale = std::max(std::fabs(plus),
                                               std::fabs(minus));
    metrics.momentum_roundoff_scale = roundoff;
    const double denominator = std::max(metrics.momentum_response_scale,
        4096.0 * std::numeric_limits<double>::epsilon() * roundoff);
    metrics.momentum_odd_relative = denominator == 0.0
        ? (metrics.momentum_odd_absolute == 0.0 ? 0.0 :
           std::numeric_limits<double>::infinity())
        : metrics.momentum_odd_absolute / denominator;
    metrics.state_odd_relative = global_odd_about_reference_relative_l2(
        positive.dual.final_state_mass, negative.dual.final_state_mass,
        zero_reference.dual.final_state_mass);
    metrics.u_flux_odd_relative = global_odd_about_reference_relative_l2(
        positive.dual.fu_center, negative.dual.fu_center,
        zero_reference.dual.fu_center);
    metrics.jn_even_relative = global_linear_combination_relative_l2(
        positive.dual.jn_high, 1.0, negative.dual.jn_high, -1.0,
        zero_reference.dual.jn_high, 0.0);
    metrics.je_even_relative = global_linear_combination_relative_l2(
        positive.dual.gstar_je_center, 1.0,
        negative.dual.gstar_je_center, -1.0,
        zero_reference.dual.gstar_je_center, 0.0);
    return metrics;
}

EndpointBalanceMetrics endpoint_energy_balance(
    const Species& background, const EMFields& fields,
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle& bundle,
    const SpatialGrid& sg, double dt, int rank, int size)
{
    EndpointBalanceMetrics metrics;
    long double local_delta_k = 0.0L;
    long double local_field_work = 0.0L;
    long double local_gross_kinetic_scale = 0.0L;
    long double local_reconstructed_u_energy = 0.0L;
    long double local_reconstructed_field_work = 0.0L;
    long double local_u_work_abs_scale = 0.0L;
    long double local_update_energy_defect = 0.0L;
    long double local_update_energy_defect_abs = 0.0L;
    long double local_update_energy_ulp_bound = 0.0L;
    double local_flux_coefficient_linf = 0.0;
    double local_flux_coefficient_scale = 0.0;
    const size_t expected_u_size = static_cast<size_t>(sg.nx_local) *
        (Param::Nv + 1) * Param::Nmu;
    const bool have_final_coefficient =
        bundle.cu_final.size() == expected_u_size &&
        bundle.fu_final.size() == expected_u_size;
    const size_t expected_x_size = static_cast<size_t>(sg.nx_local + 1) *
        Param::Nv * Param::Nmu;
    const bool can_replay_local_update = bundle.substeps_used == 1 &&
        bundle.fx_final.size() == expected_x_size &&
        bundle.fu_final.size() == expected_u_size;
    metrics.final_coefficient_available = have_final_coefficient ? 1 : 0;
    metrics.local_update_replay_available = can_replay_local_update ? 1 : 0;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t local_id =
                    (static_cast<size_t>(ix) * Param::Nv + j) * Param::Nmu + k;
                const long double final_mass = static_cast<long double>(
                    bundle.final_state_mass[local_id]);
                const long double initial_mass = static_cast<long double>(
                    background.f[idx3(sg.nghost + ix, j, k)]);
                const long double kinetic_energy = static_cast<long double>(
                    background.cgrid.kinetic_energy[idx2(j, k)]);
                local_delta_k += kinetic_energy * (final_mass - initial_mass);
                local_gross_kinetic_scale += std::fabs(kinetic_energy) *
                    (std::fabs(final_mass) + std::fabs(initial_mass));
                if (can_replay_local_update) {
                    const size_t x_left =
                        (static_cast<size_t>(ix) * Param::Nv + j) *
                        Param::Nmu + k;
                    const size_t x_right =
                        (static_cast<size_t>(ix + 1) * Param::Nv + j) *
                        Param::Nmu + k;
                    const size_t u_lower =
                        (static_cast<size_t>(ix) * (Param::Nv + 1) + j) *
                        Param::Nmu + k;
                    const size_t u_upper = u_lower + Param::Nmu;
                    const long double divergence =
                        static_cast<long double>(bundle.fx_final[x_right]) -
                        static_cast<long double>(bundle.fx_final[x_left]) +
                        static_cast<long double>(bundle.fu_final[u_upper]) -
                        static_cast<long double>(bundle.fu_final[u_lower]);
                    const long double replayed_mass = initial_mass -
                        static_cast<long double>(dt) * divergence;
                    const long double energy_defect = kinetic_energy *
                        (final_mass - replayed_mass);
                    local_update_energy_defect += energy_defect;
                    local_update_energy_defect_abs += std::fabs(energy_defect);

                    const double final_value =
                        bundle.final_state_mass[local_id];
                    const double upward = std::nextafter(final_value,
                        std::numeric_limits<double>::infinity());
                    const double downward = std::nextafter(final_value,
                        -std::numeric_limits<double>::infinity());
                    const long double local_ulp = std::max(
                        std::fabs(static_cast<long double>(upward) - final_mass),
                        std::fabs(final_mass -
                                  static_cast<long double>(downward)));
                    local_update_energy_ulp_bound +=
                        0.5L * std::fabs(kinetic_energy) * local_ulp;
                }
            }
        // Endpoint closure must use the time-averaged current assembled from
        // the final u-face fluxes used by the complete FV update.  The
        // dual_je_cell layer belongs to the local frozen dual construction;
        // it is not the production work current accumulated over all split
        // substeps.
        local_field_work += static_cast<long double>(dt * sg.dx) *
            fields.Ex[sg.nghost + ix] *
            bundle.je_final[static_cast<size_t>(ix)];
        if (have_final_coefficient) {
            const long double acceleration = static_cast<long double>(
                background.charge) * fields.Ex[sg.nghost + ix] /
                (static_cast<long double>(background.mass) * Const::c);
            const long double current_factor = static_cast<long double>(
                background.charge) /
                (static_cast<long double>(background.mass) * Const::c * sg.dx);
            for (int jf = 1; jf < Param::Nv; ++jf)
                for (int k = 0; k < Param::Nmu; ++k) {
                    const size_t u_id =
                        (static_cast<size_t>(ix) * (Param::Nv + 1) + jf) *
                        Param::Nmu + k;
                    const long double dke = static_cast<long double>(
                        Stage5::delta_energy(background.cgrid, jf, k));
                    const long double coefficient = bundle.cu_final[u_id];
                    const long double flux = bundle.fu_final[u_id];
                    const long double expected_flux = acceleration * coefficient;
                    local_flux_coefficient_linf = std::max(
                        local_flux_coefficient_linf,
                        static_cast<double>(std::fabs(flux - expected_flux)));
                    local_flux_coefficient_scale = std::max(
                        local_flux_coefficient_scale,
                        static_cast<double>(std::max(std::fabs(flux),
                                                     std::fabs(expected_flux))));
                    const long double u_term = static_cast<long double>(dt) *
                        dke * flux;
                    const long double current_work_term =
                        static_cast<long double>(dt * sg.dx) *
                        fields.Ex[sg.nghost + ix] * dke * current_factor *
                        coefficient;
                    local_reconstructed_u_energy += u_term;
                    local_reconstructed_field_work += current_work_term;
                    local_u_work_abs_scale += std::fabs(u_term) +
                        std::fabs(current_work_term);
                }
        }
    }
    const long double global_delta_k = rank_ordered_sum_long_double(
        local_delta_k, rank, size);
    const long double global_field_work = rank_ordered_sum_long_double(
        local_field_work, rank, size);
    const long double global_gross_kinetic_scale = rank_ordered_sum_long_double(
        local_gross_kinetic_scale, rank, size);
    const long double global_reconstructed_u_energy =
        rank_ordered_sum_long_double(local_reconstructed_u_energy, rank, size);
    const long double global_reconstructed_field_work =
        rank_ordered_sum_long_double(local_reconstructed_field_work, rank, size);
    const long double global_u_work_abs_scale = rank_ordered_sum_long_double(
        local_u_work_abs_scale, rank, size);
    const long double global_update_energy_defect =
        rank_ordered_sum_long_double(local_update_energy_defect, rank, size);
    const long double global_update_energy_defect_abs =
        rank_ordered_sum_long_double(local_update_energy_defect_abs, rank, size);
    const long double global_update_energy_ulp_bound =
        rank_ordered_sum_long_double(local_update_energy_ulp_bound, rank, size);
    metrics.delta_kinetic_energy = static_cast<double>(global_delta_k);
    metrics.field_work = static_cast<double>(global_field_work);
    metrics.final_state_energy_gross_scale = static_cast<double>(
        global_gross_kinetic_scale);
    metrics.local_update_energy_defect_signed = static_cast<double>(
        global_update_energy_defect);
    metrics.local_update_energy_defect_absolute = static_cast<double>(
        global_update_energy_defect_abs);
    metrics.local_update_energy_ulp_bound = static_cast<double>(
        global_update_energy_ulp_bound);
    metrics.reconstructed_u_energy_moment = static_cast<double>(
        global_reconstructed_u_energy);
    metrics.reconstructed_field_work = static_cast<double>(
        global_reconstructed_field_work);
    double flux_relation[2] = {local_flux_coefficient_linf,
                               local_flux_coefficient_scale};
    MPI_Allreduce(MPI_IN_PLACE, flux_relation, 2, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    metrics.u_flux_coefficient_relation_linf = flux_relation[0];
    metrics.u_flux_coefficient_relation_scale = flux_relation[1];
    metrics.u_flux_coefficient_relation_relative = relative_to_nonzero_scale(
        flux_relation[0], flux_relation[1]);
    metrics.boundary_work = bundle.u_boundary_energy;
    metrics.spatial_boundary_work = bundle.stage5_spatial_energy_boundary;
    metrics.production_delta_kinetic_energy = bundle.delta_ke_bkg;
    metrics.production_u_energy_moment = bundle.stage5_u_energy_moment;

    // Independent physical balance: final-state kinetic-energy change plus
    // periodic x-boundary transport and u-boundary work, minus the work formed
    // from the complete time-averaged final u current.
    const long double balance_residual = global_delta_k +
        static_cast<long double>(metrics.spatial_boundary_work) +
        static_cast<long double>(metrics.boundary_work) - global_field_work;
    metrics.residual = static_cast<double>(balance_residual);
    metrics.scale = std::fabs(metrics.delta_kinetic_energy) +
                    std::fabs(metrics.spatial_boundary_work) +
                    std::fabs(metrics.boundary_work) +
                    std::fabs(metrics.field_work);
    metrics.relative = relative_to_nonzero_scale(metrics.residual,
                                                 metrics.scale);
    metrics.r_fv_match_relative = relative_to_nonzero_scale(
        metrics.residual - bundle.r_fv,
        std::max(metrics.scale, std::fabs(bundle.r_fv)));

    // Exact reconstruction of the production R_FV definition.  This uses
    // values exported by the production operator rather than replayed fluxes.
    metrics.production_r_fv_reconstructed =
        metrics.production_delta_kinetic_energy +
        metrics.spatial_boundary_work - metrics.production_u_energy_moment +
        metrics.boundary_work;
    const double production_scale =
        std::fabs(metrics.production_delta_kinetic_energy) +
        std::fabs(metrics.spatial_boundary_work) +
        std::fabs(metrics.production_u_energy_moment) +
        std::fabs(metrics.boundary_work);
    metrics.production_r_fv_reconstruction_relative =
        relative_to_nonzero_scale(
            metrics.production_r_fv_reconstructed - bundle.r_fv,
            std::max(production_scale, std::fabs(bundle.r_fv)));
    metrics.delta_kinetic_energy_match_relative = relative_to_nonzero_scale(
        metrics.delta_kinetic_energy -
            metrics.production_delta_kinetic_energy,
        std::max(std::fabs(metrics.delta_kinetic_energy),
                 std::fabs(metrics.production_delta_kinetic_energy)));
    metrics.u_energy_field_work_difference_absolute = std::fabs(
        metrics.production_u_energy_moment - metrics.field_work);
    // Both work terms can vanish by symmetry.  Normalizing their roundoff-level
    // difference by either near-zero term turns harmless cancellation into an
    // order-one error.  Use the physical energy exchanged in this endpoint
    // event instead; retain the absolute difference for transparent auditing.
    metrics.u_energy_field_work_scale = std::max(metrics.scale,
                                                  production_scale);
    metrics.u_energy_field_work_match_relative = relative_to_nonzero_scale(
        metrics.u_energy_field_work_difference_absolute,
        metrics.u_energy_field_work_scale);
    metrics.production_u_energy_reconstruction_relative =
        relative_to_nonzero_scale(
            metrics.production_u_energy_moment -
                metrics.reconstructed_u_energy_moment,
            metrics.u_energy_field_work_scale);
    metrics.final_current_work_reconstruction_relative =
        relative_to_nonzero_scale(
            metrics.field_work - metrics.reconstructed_field_work,
            metrics.u_energy_field_work_scale);
    const long double operation_count = static_cast<long double>(
        (Param::Nv - 1) * Param::Nmu + sg.nx_global + 32);
    metrics.u_work_roundoff_bound = static_cast<double>(
        8.0L * std::numeric_limits<double>::epsilon() * operation_count *
        global_u_work_abs_scale);
    metrics.u_work_roundoff_limited =
        metrics.u_energy_field_work_difference_absolute <=
        metrics.u_work_roundoff_bound ? 1 : 0;
    // Tight endpoint budget: replay the actual accepted one-substep FV update
    // cell by cell, then add only local storage ULP and the independently
    // measured current/work accumulation floor.  The gross background kinetic
    // energy is retained as a diagnostic but is no longer used as a bound.
    const long double long_double_sum_floor =
        8.0L * std::numeric_limits<long double>::epsilon() * operation_count *
        (global_update_energy_defect_abs + global_u_work_abs_scale);
    metrics.balance_roundoff_bound = can_replay_local_update
        ? static_cast<double>(global_update_energy_defect_abs +
              global_update_energy_ulp_bound +
              metrics.u_energy_field_work_difference_absolute +
              long_double_sum_floor)
        : 0.0;
    metrics.balance_roundoff_limited = std::fabs(metrics.residual) <=
        metrics.balance_roundoff_bound ? 1 : 0;
    metrics.r_fv_match_difference_absolute = std::fabs(
        metrics.residual - bundle.r_fv);
    const double delta_k_path_difference = std::fabs(
        metrics.delta_kinetic_energy -
        metrics.production_delta_kinetic_energy);
    const double production_reconstruction_difference = std::fabs(
        metrics.production_r_fv_reconstructed - bundle.r_fv);
    metrics.r_fv_match_roundoff_bound = delta_k_path_difference +
        metrics.u_energy_field_work_difference_absolute +
        production_reconstruction_difference +
        metrics.local_update_energy_ulp_bound;
    metrics.r_fv_match_roundoff_limited =
        metrics.r_fv_match_difference_absolute <=
        metrics.r_fv_match_roundoff_bound ? 1 : 0;
    metrics.spatial_boundary_relative = relative_to_nonzero_scale(
        metrics.spatial_boundary_work, production_scale);
    return metrics;
}

void write_metrics(std::ostream& out, const PairMetrics& metrics)
{
    out << std::scientific << std::setprecision(17)
        << "state_advanced=" << metrics.dual.state_advanced << "\n"
        << "operator_failed=" << metrics.dual.operator_failed << "\n"
        << "outputs_finite=" << metrics.dual.outputs_finite << "\n"
        << "dual_u_operator_valid=" << metrics.dual.dual_u_operator_valid << "\n"
        << "mass_error_absolute=" << std::fabs(metrics.dual.mass_residual) << "\n"
        << "mass_scale=" << metrics.dual.mass_scale << "\n"
        << "mass_error_relative=" << relative_error(
            metrics.dual.mass_residual, metrics.dual.mass_scale) << "\n"
        << "momentum_error_absolute="
        << std::fabs(metrics.dual.momentum_residual) << "\n"
        << "momentum_scale=" << metrics.dual.momentum_scale << "\n"
        << "momentum_error_relative=" << relative_error(
            metrics.dual.momentum_residual, metrics.dual.momentum_scale)
        << "\n"
        << "u_boundary_work_left_plus_right="
        << metrics.dual.u_boundary_energy << "\n"
        << "u_boundary_work_total=" << metrics.dual.u_boundary_energy << "\n"
        << "u_boundary_work_left=" << metrics.dual.u_boundary_energy_lower
        << "\nu_boundary_work_right="
        << metrics.dual.u_boundary_energy_upper << "\n"
        << "u_boundary_particle=" << metrics.dual.u_boundary_particle << "\n"
        << "u_boundary_momentum=" << metrics.dual.u_boundary_momentum << "\n"
        << "JN_minus_GstarJE_signed=" << metrics.dual_signed << "\n"
        << "JN_minus_GstarJE_L1=" << metrics.dual_pair.l1 << "\n"
        << "JN_minus_GstarJE_L2=" << metrics.dual_pair.l2 << "\n"
        << "JN_minus_GstarJE_Linf=" << metrics.dual_pair.linf << "\n"
        << "dual_cell_JE_minus_target_signed="
        << metrics.dual_cell_target_signed << "\n"
        << "dual_cell_JE_minus_target_L1=" << metrics.dual_cell_target.l1
        << "\ndual_cell_JE_minus_target_L2=" << metrics.dual_cell_target.l2
        << "\ndual_cell_JE_minus_target_Linf=" << metrics.dual_cell_target.linf
        << "\ndual_cell_JE_minus_target_relative_L2="
        << relative_to_nonzero_scale(metrics.dual_cell_target.l2,
                                     metrics.dual_cell_target_scale.l2)
        << "\nprojection_floor_signed=" << metrics.projection_floor_signed
        << "\nprojection_floor_L1=" << metrics.projection_floor.l1
        << "\nprojection_floor_L2=" << metrics.projection_floor.l2
        << "\nprojection_floor_Linf=" << metrics.projection_floor.linf
        << "\nprojection_reconstruction_signed="
        << metrics.projection_reconstruction_signed
        << "\nprojection_reconstruction_L1="
        << metrics.projection_reconstruction.l1
        << "\nprojection_reconstruction_L2="
        << metrics.projection_reconstruction.l2
        << "\nprojection_reconstruction_Linf="
        << metrics.projection_reconstruction.linf
        << "\nprojection_reconstruction_relative_L2="
        << relative_to_nonzero_scale(metrics.projection_reconstruction.l2,
                                     metrics.projection_target_scale.l2)
        << "\nlegacy_JN_minus_GstarJE_signed=" << metrics.legacy_signed << "\n"
        << "legacy_JN_minus_GstarJE_L1=" << metrics.legacy_pair.l1 << "\n"
        << "legacy_JN_minus_GstarJE_L2=" << metrics.legacy_pair.l2 << "\n"
        << "legacy_JN_minus_GstarJE_Linf=" << metrics.legacy_pair.linf << "\n"
        << "min_f_legacy=" << metrics.legacy.final_candidate_min << "\n"
        << "min_f_dual=" << metrics.dual.final_candidate_min << "\n"
        << "legacy_dual_JN_Linf=" << metrics.jn_difference << "\n"
        << "legacy_dual_JE_Linf=" << metrics.je_difference << "\n"
        << "legacy_dual_state_L1=" << metrics.state_difference_l1 << "\n"
        << "legacy_dual_state_relative_L1="
        << relative_to_nonzero_scale(metrics.state_difference_l1,
                                     metrics.dual.mass_scale)
        << "\nlegacy_dual_R_couple_difference="
        << metrics.dual.r_couple - metrics.legacy.r_couple << "\n"
        << "R_FV=" << metrics.dual.r_fv << "\n"
        << "R_couple=" << metrics.dual.r_couple << "\n";
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
    const bool smooth_case = FP_DUAL_U_STAGE2_CASE == PHASE_TRANSLATION ||
        FP_DUAL_U_STAGE2_CASE == NONUNIFORM_GRID ||
        FP_DUAL_U_STAGE2_CASE == MANUFACTURED ||
        FP_DUAL_U_STAGE2_CASE == MPI_CONSISTENCY;
    const int spatial_mode = FP_DUAL_U_STAGE2_CASE == PHASE_TRANSLATION ? 4 : 3;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size,
        smooth_case ? 0.04 : 0.0, smooth_case ? 0.06 : 0.0,
        smooth_case ? 1.0e3 : 0.0, spatial_mode, 0.0);
    if (FP_DUAL_U_STAGE2_CASE == ENDPOINT_FLUX) {
        add_endpoint_population(background, sg, rank, size);
        set_uniform_field(fields, sg, rank, size, 1.0e3);
    }
    const std::vector<double> initial_zero_state =
        FP_DUAL_U_STAGE2_CASE == ZERO_FIELD
        ? physical_state(background, sg) : std::vector<double>();
    const double dt = BackgroundCouplingTest::stable_dt(sg);
    PairMetrics primary = evaluate_pair(background, fields, sg, rank, size, dt);
    const PairMetrics zero_reference = primary;
    PairMetrics secondary = primary;
    FieldSymmetryMetrics field_symmetry;
    PhaseMetrics phase_metrics;
    ZeroFieldMetrics zero_field_metrics;
    EndpointBalanceMetrics endpoint_positive_balance;
    EndpointBalanceMetrics endpoint_negative_balance;
    double symmetry_error = 0.0;
    double phase_norm_error = 0.0;
    double endpoint_positive_ledger_relative = 0.0;
    double endpoint_negative_ledger_relative = 0.0;
    double endpoint_negative_work_left = 0.0;
    double endpoint_negative_work_right = 0.0;
    int secondary_valid = 1;

    if (FP_DUAL_U_STAGE2_CASE == ZERO_FIELD && primary.valid) {
        zero_field_metrics.state_relative_l2 =
            global_linear_combination_relative_l2(
                primary.dual.final_state_mass, 1.0, initial_zero_state, -1.0,
                std::vector<double>(initial_zero_state.size(), 0.0), 0.0);
        zero_field_metrics.u_flux_linf =
            BackgroundCouplingTest::global_vector_abs_linf(
                primary.dual.fu_center);
        zero_field_metrics.jn_linf = max_abs_owned_face(primary.dual.jn_high,
                                                        sg);
        zero_field_metrics.cancellation_current_scale =
            cancellation_current_scale(background, sg);
        zero_field_metrics.jn_relative_to_cancellation_scale =
            relative_to_nonzero_scale(
                zero_field_metrics.jn_linf,
                zero_field_metrics.cancellation_current_scale);
    }

    if (FP_DUAL_U_STAGE2_CASE == FIELD_SYMMETRY) {
        set_uniform_field(fields, sg, rank, size, 1.0e-3);
        const PairMetrics positive = evaluate_pair(background, fields, sg,
                                                   rank, size, dt);
        set_uniform_field(fields, sg, rank, size, -1.0e-3);
        const PairMetrics negative = evaluate_pair(background, fields, sg,
                                                   rank, size, dt);
        primary = positive;
        secondary = negative;
        secondary_valid = positive.valid && negative.valid;
        if (!secondary_valid) {
            symmetry_error = std::numeric_limits<double>::infinity();
        } else {
        field_symmetry = deterministic_field_symmetry(
            positive, negative, zero_reference, background, sg, rank, size);
        symmetry_error = field_symmetry.maximum();
        }
    } else if (FP_DUAL_U_STAGE2_CASE == PHASE_TRANSLATION) {
        Species shifted_background;
        EMFields shifted_fields;
        BackgroundCouplingTest::initialize_periodic_state(
            shifted_background, shifted_fields, sg, rank, size,
            0.04, 0.06, 1.0e3, spatial_mode, 0.5 * Const::pi);
        const PairMetrics shifted = evaluate_pair(shifted_background,
            shifted_fields, sg, rank, size, dt);
        secondary = shifted;
        secondary_valid = primary.valid && shifted.valid;
        if (!secondary_valid) {
            phase_norm_error = std::numeric_limits<double>::infinity();
        } else {
        const int shift_denominator = 4 * spatial_mode;
        const bool integer_shift = sg.nx_global % shift_denominator == 0;
        const int shift_cells = integer_shift ?
            sg.nx_global / shift_denominator : 0;
        std::vector<double> primary_jn(primary.dual.jn_high.begin(),
            primary.dual.jn_high.begin() + sg.nx_local);
        std::vector<double> shifted_jn(shifted.dual.jn_high.begin(),
            shifted.dual.jn_high.begin() + sg.nx_local);
        std::vector<double> primary_je(primary.dual.gstar_je_center.begin(),
            primary.dual.gstar_je_center.begin() + sg.nx_local);
        std::vector<double> shifted_je(shifted.dual.gstar_je_center.begin(),
            shifted.dual.gstar_je_center.begin() + sg.nx_local);
        std::vector<double> primary_residual(static_cast<size_t>(sg.nx_local));
        std::vector<double> shifted_residual(static_cast<size_t>(sg.nx_local));
        for (int ix = 0; ix < sg.nx_local; ++ix) {
            primary_residual[static_cast<size_t>(ix)] =
                primary_jn[static_cast<size_t>(ix)] -
                primary_je[static_cast<size_t>(ix)];
            shifted_residual[static_cast<size_t>(ix)] =
                shifted_jn[static_cast<size_t>(ix)] -
                shifted_je[static_cast<size_t>(ix)];
        }
        std::vector<double> primary_global = gather_rank_ordered(
            primary_jn, rank, size);
        std::vector<double> shifted_global = gather_rank_ordered(
            shifted_jn, rank, size);
        phase_metrics.jn_relative_l2 = shifted_relative_l2(
            primary_global, shifted_global, sg.nx_global, 1, shift_cells,
            rank);
        primary_global = gather_rank_ordered(primary_je, rank, size);
        shifted_global = gather_rank_ordered(shifted_je, rank, size);
        phase_metrics.gstar_je_relative_l2 = shifted_relative_l2(
            primary_global, shifted_global, sg.nx_global, 1, shift_cells,
            rank);
        primary_global = gather_rank_ordered(primary_residual, rank, size);
        shifted_global = gather_rank_ordered(shifted_residual, rank, size);
        phase_metrics.residual_relative_l2 = shifted_relative_l2(
            primary_global, shifted_global, sg.nx_global, 1, shift_cells,
            rank);
        primary_global = gather_rank_ordered(
            primary.dual.final_state_mass, rank, size);
        shifted_global = gather_rank_ordered(
            shifted.dual.final_state_mass, rank, size);
        phase_metrics.state_relative_l2 = shifted_relative_l2(
            primary_global, shifted_global, sg.nx_global, Param::Nvmu,
            shift_cells, rank);
        phase_norm_error = integer_shift ? phase_metrics.maximum() :
            std::numeric_limits<double>::infinity();
        }
    } else if (FP_DUAL_U_STAGE2_CASE == ENDPOINT_FLUX) {
        if (primary.valid)
            endpoint_positive_balance = endpoint_energy_balance(
                background, fields, primary.dual, sg, dt, rank, size);
        else {
            endpoint_positive_balance.relative =
                std::numeric_limits<double>::infinity();
            endpoint_positive_balance.r_fv_match_relative =
                std::numeric_limits<double>::infinity();
        }
        set_uniform_field(fields, sg, rank, size, -1.0e3);
        secondary = evaluate_pair(background, fields, sg, rank, size, dt);
        secondary_valid = primary.valid && secondary.valid;
        if (!secondary_valid) {
            endpoint_positive_ledger_relative =
                std::numeric_limits<double>::infinity();
            endpoint_negative_ledger_relative =
                std::numeric_limits<double>::infinity();
        } else {
        const double positive_scale = std::max(
            std::fabs(primary.dual.u_boundary_energy),
            std::fabs(primary.dual.u_boundary_energy_lower) +
            std::fabs(primary.dual.u_boundary_energy_upper));
        endpoint_positive_ledger_relative = relative_to_nonzero_scale(
            primary.dual.u_boundary_energy -
            primary.dual.u_boundary_energy_lower -
            primary.dual.u_boundary_energy_upper, positive_scale);
        const double negative_scale = std::max(
            std::fabs(secondary.dual.u_boundary_energy),
            std::fabs(secondary.dual.u_boundary_energy_lower) +
            std::fabs(secondary.dual.u_boundary_energy_upper));
        endpoint_negative_ledger_relative = relative_to_nonzero_scale(
            secondary.dual.u_boundary_energy -
            secondary.dual.u_boundary_energy_lower -
            secondary.dual.u_boundary_energy_upper, negative_scale);
        endpoint_negative_work_left = secondary.dual.u_boundary_energy_lower;
        endpoint_negative_work_right = secondary.dual.u_boundary_energy_upper;
        endpoint_negative_balance = endpoint_energy_balance(
            background, fields, secondary.dual, sg, dt, rank, size);
        }
    }

    const auto conservation_passes = [](const PairMetrics& metrics) {
        const double mass_scale = std::max(1.0, metrics.dual.mass_scale);
        const double momentum_scale = std::max(1.0,
                                                metrics.dual.momentum_scale);
        return std::fabs(metrics.dual.mass_residual) <= 1.0e-10 * mass_scale &&
            std::fabs(metrics.dual.momentum_residual) <=
                1.0e-9 * momentum_scale;
    };
    const auto pair_improves = [](const PairMetrics& metrics) {
        return metrics.dual_pair.l2 <=
            metrics.legacy_pair.l2 * (1.0 + 1.0e-10) + 1.0e-12;
    };
    const bool has_secondary_case =
        FP_DUAL_U_STAGE2_CASE == FIELD_SYMMETRY ||
        FP_DUAL_U_STAGE2_CASE == PHASE_TRANSLATION ||
        FP_DUAL_U_STAGE2_CASE == ENDPOINT_FLUX;
    const bool conservation_ok = conservation_passes(primary) &&
        (!has_secondary_case || conservation_passes(secondary));
    const bool pair_ok = pair_improves(primary) &&
        (!has_secondary_case || pair_improves(secondary));
    const auto decomposition_passes = [](const PairMetrics& metrics) {
        return relative_to_nonzero_scale(metrics.dual_cell_target.l2,
                   metrics.dual_cell_target_scale.l2) <= 1.0e-10 &&
            relative_to_nonzero_scale(metrics.projection_reconstruction.l2,
                   metrics.projection_target_scale.l2) <= 1.0e-10;
    };
    const bool dual_decomposition_applicable =
        FP_DUAL_U_STAGE2_CASE == PHASE_TRANSLATION ||
        FP_DUAL_U_STAGE2_CASE == NONUNIFORM_GRID ||
        FP_DUAL_U_STAGE2_CASE == MANUFACTURED ||
        FP_DUAL_U_STAGE2_CASE == MPI_CONSISTENCY;
    const bool dual_decomposition_ok = !dual_decomposition_applicable ||
        (decomposition_passes(primary) &&
         (!has_secondary_case || decomposition_passes(secondary)));
    const bool zero_field_ok = FP_DUAL_U_STAGE2_CASE != ZERO_FIELD ||
        (zero_field_metrics.state_relative_l2 <= 1.0e-14 &&
         zero_field_metrics.u_flux_linf == 0.0 &&
         zero_field_metrics.jn_relative_to_cancellation_scale <= 1.0e-12 &&
         primary.dual.u_boundary_energy == 0.0);
    const bool endpoint_ledger_ok = FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        (std::fabs(primary.dual.u_boundary_energy_lower) > 0.0 &&
         std::fabs(secondary.dual.u_boundary_energy_upper) > 0.0 &&
         endpoint_positive_ledger_relative <= 1.0e-12 &&
         endpoint_negative_ledger_relative <= 1.0e-12);
    const bool endpoint_production_reconstruction_ok =
        FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        (endpoint_positive_balance.production_r_fv_reconstruction_relative <=
             1.0e-10 &&
         endpoint_negative_balance.production_r_fv_reconstruction_relative <=
             1.0e-10);
    const bool endpoint_delta_ke_match_ok =
        FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        (endpoint_positive_balance.delta_kinetic_energy_match_relative <=
             1.0e-10 &&
         endpoint_negative_balance.delta_kinetic_energy_match_relative <=
             1.0e-10);
    const bool endpoint_u_work_match_ok =
        FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        (endpoint_positive_balance.final_coefficient_available == 1 &&
         endpoint_negative_balance.final_coefficient_available == 1 &&
         endpoint_positive_balance.u_flux_coefficient_relation_relative <=
             1.0e-12 &&
         endpoint_negative_balance.u_flux_coefficient_relation_relative <=
             1.0e-12 &&
         endpoint_positive_balance.production_u_energy_reconstruction_relative <=
             std::max(1.0e-12,
                 endpoint_positive_balance.u_work_roundoff_bound /
                 std::max(endpoint_positive_balance.u_energy_field_work_scale,
                          std::numeric_limits<double>::min())) &&
         endpoint_negative_balance.production_u_energy_reconstruction_relative <=
             std::max(1.0e-12,
                 endpoint_negative_balance.u_work_roundoff_bound /
                 std::max(endpoint_negative_balance.u_energy_field_work_scale,
                          std::numeric_limits<double>::min())) &&
         endpoint_positive_balance.final_current_work_reconstruction_relative <=
             std::max(1.0e-12,
                 endpoint_positive_balance.u_work_roundoff_bound /
                 std::max(endpoint_positive_balance.u_energy_field_work_scale,
                          std::numeric_limits<double>::min())) &&
         endpoint_negative_balance.final_current_work_reconstruction_relative <=
             std::max(1.0e-12,
                 endpoint_negative_balance.u_work_roundoff_bound /
                 std::max(endpoint_negative_balance.u_energy_field_work_scale,
                          std::numeric_limits<double>::min())) &&
         (endpoint_positive_balance.u_energy_field_work_match_relative <=
              1.0e-10 ||
          endpoint_positive_balance.u_work_roundoff_limited == 1) &&
         (endpoint_negative_balance.u_energy_field_work_match_relative <=
              1.0e-10 ||
          endpoint_negative_balance.u_work_roundoff_limited == 1));
    const bool endpoint_spatial_boundary_ok =
        FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        (endpoint_positive_balance.spatial_boundary_relative <= 1.0e-12 &&
         endpoint_negative_balance.spatial_boundary_relative <= 1.0e-12);
    const bool endpoint_balance_ok = FP_DUAL_U_STAGE2_CASE != ENDPOINT_FLUX ||
        ((endpoint_positive_balance.relative <= 1.0e-10 ||
          endpoint_positive_balance.balance_roundoff_limited == 1) &&
         (endpoint_negative_balance.relative <= 1.0e-10 ||
          endpoint_negative_balance.balance_roundoff_limited == 1) &&
         (endpoint_positive_balance.r_fv_match_relative <= 1.0e-10 ||
          endpoint_positive_balance.r_fv_match_roundoff_limited == 1) &&
         (endpoint_negative_balance.r_fv_match_relative <= 1.0e-10 ||
          endpoint_negative_balance.r_fv_match_roundoff_limited == 1) &&
         endpoint_production_reconstruction_ok &&
         endpoint_delta_ke_match_ok && endpoint_u_work_match_ok &&
         endpoint_spatial_boundary_ok);
    const bool endpoint_ok = endpoint_ledger_ok && endpoint_balance_ok;
    const bool symmetry_ok = FP_DUAL_U_STAGE2_CASE != FIELD_SYMMETRY ||
        symmetry_error <= 1.0e-8;
    const bool phase_ok = FP_DUAL_U_STAGE2_CASE != PHASE_TRANSLATION ||
        phase_norm_error <= 1.0e-8;
    const bool nonuniform_active_ok =
        FP_DUAL_U_STAGE2_CASE != NONUNIFORM_GRID ||
        primary.dual.dual_u_corrected_cell_count > 0;
    int passes = primary.valid && secondary_valid && conservation_ok &&
        pair_ok && dual_decomposition_ok && zero_field_ok && endpoint_ok &&
        symmetry_ok && phase_ok && nonuniform_active_ok;
    MPI_Allreduce(MPI_IN_PLACE, &passes, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    const std::string outdir = output_directory(argc, argv);
    if (rank == 0) {
        make_directory_tree(outdir);
        std::ofstream out((outdir + "/summary.result").c_str());
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_dual_u_stage2\n"
            << "case=" << case_name() << "\n"
            << "production_kernel=1\nbeam_enabled=0\nfct_enabled=0\n"
            << "nx=" << Param::nx << "\nNupar=" << Param::Nv
            << "\nNuperp=" << Param::Nmu << "\nmpi_size=" << size << "\n"
            << "git_revision=" << FP_DUAL_U_STAGE2_GIT_REVISION
            << "\ngit_dirty=" << FP_DUAL_U_STAGE2_GIT_DIRTY << "\n"
            << "velocity_grid_nonuniform=" << !background.cgrid.is_uniform()
            << "\nupar_grid_hash=" << grid_hash(background.cgrid.upar_widths)
            << "\nuperp_grid_hash=" << grid_hash(background.cgrid.uperp_widths)
            << "\n";
        write_metrics(log, primary);
        log << std::scientific << std::setprecision(17)
            << "positive_negative_field_symmetry_error=" << symmetry_error
            << "\nfield_symmetry_momentum_odd_absolute="
            << field_symmetry.momentum_odd_absolute
            << "\nfield_symmetry_momentum_response_scale="
            << field_symmetry.momentum_response_scale
            << "\nfield_symmetry_momentum_roundoff_scale="
            << field_symmetry.momentum_roundoff_scale
            << "\nfield_symmetry_momentum_odd_relative="
            << field_symmetry.momentum_odd_relative
            << "\nfield_symmetry_state_odd_relative="
            << field_symmetry.state_odd_relative
            << "\nfield_symmetry_u_flux_odd_relative="
            << field_symmetry.u_flux_odd_relative
            << "\nfield_symmetry_JN_even_relative="
            << field_symmetry.jn_even_relative
            << "\nfield_symmetry_GstarJE_even_relative="
            << field_symmetry.je_even_relative
            << "\nphase_translation_norm_error=" << phase_norm_error
            << "\nphase_translation_JN_relative_L2="
            << phase_metrics.jn_relative_l2
            << "\nphase_translation_GstarJE_relative_L2="
            << phase_metrics.gstar_je_relative_l2
            << "\nphase_translation_residual_relative_L2="
            << phase_metrics.residual_relative_l2
            << "\nphase_translation_state_relative_L2="
            << phase_metrics.state_relative_l2
            << "\nzero_field_state_relative_L2="
            << zero_field_metrics.state_relative_l2
            << "\nzero_field_u_flux_Linf="
            << zero_field_metrics.u_flux_linf
            << "\nzero_field_JN_Linf=" << zero_field_metrics.jn_linf
            << "\nzero_field_cancellation_current_scale="
            << zero_field_metrics.cancellation_current_scale
            << "\nzero_field_JN_relative_to_cancellation_scale="
            << zero_field_metrics.jn_relative_to_cancellation_scale
            << "\nendpoint_work_uses_final_current=1"
            << "\nendpoint_positive_ledger_relative="
            << endpoint_positive_ledger_relative
            << "\nendpoint_negative_ledger_relative="
            << endpoint_negative_ledger_relative
            << "\nendpoint_negative_work_left="
            << endpoint_negative_work_left
            << "\nendpoint_negative_work_right="
            << endpoint_negative_work_right
            << "\nendpoint_positive_delta_kinetic_energy="
            << endpoint_positive_balance.delta_kinetic_energy
            << "\nendpoint_positive_field_work="
            << endpoint_positive_balance.field_work
            << "\nendpoint_positive_boundary_work="
            << endpoint_positive_balance.boundary_work
            << "\nendpoint_positive_spatial_boundary_work="
            << endpoint_positive_balance.spatial_boundary_work
            << "\nendpoint_positive_balance_residual="
            << endpoint_positive_balance.residual
            << "\nendpoint_positive_balance_relative="
            << endpoint_positive_balance.relative
            << "\nendpoint_positive_r_fv_match_relative="
            << endpoint_positive_balance.r_fv_match_relative
            << "\nendpoint_positive_production_delta_kinetic_energy="
            << endpoint_positive_balance.production_delta_kinetic_energy
            << "\nendpoint_positive_production_u_energy_moment="
            << endpoint_positive_balance.production_u_energy_moment
            << "\nendpoint_positive_production_r_fv_reconstructed="
            << endpoint_positive_balance.production_r_fv_reconstructed
            << "\nendpoint_positive_production_r_fv_reconstruction_relative="
            << endpoint_positive_balance.production_r_fv_reconstruction_relative
            << "\nendpoint_positive_delta_kinetic_energy_match_relative="
            << endpoint_positive_balance.delta_kinetic_energy_match_relative
            << "\nendpoint_positive_u_energy_field_work_difference_absolute="
            << endpoint_positive_balance.u_energy_field_work_difference_absolute
            << "\nendpoint_positive_u_energy_field_work_scale="
            << endpoint_positive_balance.u_energy_field_work_scale
        << "\nendpoint_positive_u_energy_field_work_match_relative="
        << endpoint_positive_balance.u_energy_field_work_match_relative
        << "\nendpoint_positive_final_coefficient_available="
        << endpoint_positive_balance.final_coefficient_available
        << "\nendpoint_positive_u_flux_coefficient_relation_linf="
        << endpoint_positive_balance.u_flux_coefficient_relation_linf
        << "\nendpoint_positive_u_flux_coefficient_relation_scale="
        << endpoint_positive_balance.u_flux_coefficient_relation_scale
        << "\nendpoint_positive_u_flux_coefficient_relation_relative="
        << endpoint_positive_balance.u_flux_coefficient_relation_relative
        << "\nendpoint_positive_reconstructed_u_energy_moment="
        << endpoint_positive_balance.reconstructed_u_energy_moment
        << "\nendpoint_positive_reconstructed_field_work="
        << endpoint_positive_balance.reconstructed_field_work
        << "\nendpoint_positive_production_u_energy_reconstruction_relative="
        << endpoint_positive_balance.production_u_energy_reconstruction_relative
        << "\nendpoint_positive_final_current_work_reconstruction_relative="
        << endpoint_positive_balance.final_current_work_reconstruction_relative
        << "\nendpoint_positive_u_work_roundoff_bound="
        << endpoint_positive_balance.u_work_roundoff_bound
        << "\nendpoint_positive_u_work_roundoff_limited="
        << endpoint_positive_balance.u_work_roundoff_limited
        << "\nendpoint_positive_final_state_energy_gross_scale="
        << endpoint_positive_balance.final_state_energy_gross_scale
        << "\nendpoint_positive_local_update_replay_available="
        << endpoint_positive_balance.local_update_replay_available
        << "\nendpoint_positive_local_update_energy_defect_signed="
        << endpoint_positive_balance.local_update_energy_defect_signed
        << "\nendpoint_positive_local_update_energy_defect_absolute="
        << endpoint_positive_balance.local_update_energy_defect_absolute
        << "\nendpoint_positive_local_update_energy_ulp_bound="
        << endpoint_positive_balance.local_update_energy_ulp_bound
        << "\nendpoint_positive_balance_roundoff_bound="
        << endpoint_positive_balance.balance_roundoff_bound
        << "\nendpoint_positive_balance_roundoff_limited="
        << endpoint_positive_balance.balance_roundoff_limited
        << "\nendpoint_positive_r_fv_match_difference_absolute="
        << endpoint_positive_balance.r_fv_match_difference_absolute
        << "\nendpoint_positive_r_fv_match_roundoff_bound="
        << endpoint_positive_balance.r_fv_match_roundoff_bound
        << "\nendpoint_positive_r_fv_match_roundoff_limited="
        << endpoint_positive_balance.r_fv_match_roundoff_limited
        << "\nendpoint_positive_spatial_boundary_relative="
            << endpoint_positive_balance.spatial_boundary_relative
            << "\nendpoint_negative_delta_kinetic_energy="
            << endpoint_negative_balance.delta_kinetic_energy
            << "\nendpoint_negative_field_work="
            << endpoint_negative_balance.field_work
            << "\nendpoint_negative_boundary_work="
            << endpoint_negative_balance.boundary_work
            << "\nendpoint_negative_spatial_boundary_work="
            << endpoint_negative_balance.spatial_boundary_work
            << "\nendpoint_negative_balance_residual="
            << endpoint_negative_balance.residual
            << "\nendpoint_negative_balance_relative="
            << endpoint_negative_balance.relative
            << "\nendpoint_negative_r_fv_match_relative="
            << endpoint_negative_balance.r_fv_match_relative
            << "\nendpoint_negative_production_delta_kinetic_energy="
            << endpoint_negative_balance.production_delta_kinetic_energy
            << "\nendpoint_negative_production_u_energy_moment="
            << endpoint_negative_balance.production_u_energy_moment
            << "\nendpoint_negative_production_r_fv_reconstructed="
            << endpoint_negative_balance.production_r_fv_reconstructed
            << "\nendpoint_negative_production_r_fv_reconstruction_relative="
            << endpoint_negative_balance.production_r_fv_reconstruction_relative
            << "\nendpoint_negative_delta_kinetic_energy_match_relative="
            << endpoint_negative_balance.delta_kinetic_energy_match_relative
            << "\nendpoint_negative_u_energy_field_work_difference_absolute="
            << endpoint_negative_balance.u_energy_field_work_difference_absolute
            << "\nendpoint_negative_u_energy_field_work_scale="
            << endpoint_negative_balance.u_energy_field_work_scale
        << "\nendpoint_negative_u_energy_field_work_match_relative="
        << endpoint_negative_balance.u_energy_field_work_match_relative
        << "\nendpoint_negative_final_coefficient_available="
        << endpoint_negative_balance.final_coefficient_available
        << "\nendpoint_negative_u_flux_coefficient_relation_linf="
        << endpoint_negative_balance.u_flux_coefficient_relation_linf
        << "\nendpoint_negative_u_flux_coefficient_relation_scale="
        << endpoint_negative_balance.u_flux_coefficient_relation_scale
        << "\nendpoint_negative_u_flux_coefficient_relation_relative="
        << endpoint_negative_balance.u_flux_coefficient_relation_relative
        << "\nendpoint_negative_reconstructed_u_energy_moment="
        << endpoint_negative_balance.reconstructed_u_energy_moment
        << "\nendpoint_negative_reconstructed_field_work="
        << endpoint_negative_balance.reconstructed_field_work
        << "\nendpoint_negative_production_u_energy_reconstruction_relative="
        << endpoint_negative_balance.production_u_energy_reconstruction_relative
        << "\nendpoint_negative_final_current_work_reconstruction_relative="
        << endpoint_negative_balance.final_current_work_reconstruction_relative
        << "\nendpoint_negative_u_work_roundoff_bound="
        << endpoint_negative_balance.u_work_roundoff_bound
        << "\nendpoint_negative_u_work_roundoff_limited="
        << endpoint_negative_balance.u_work_roundoff_limited
        << "\nendpoint_negative_final_state_energy_gross_scale="
        << endpoint_negative_balance.final_state_energy_gross_scale
        << "\nendpoint_negative_local_update_replay_available="
        << endpoint_negative_balance.local_update_replay_available
        << "\nendpoint_negative_local_update_energy_defect_signed="
        << endpoint_negative_balance.local_update_energy_defect_signed
        << "\nendpoint_negative_local_update_energy_defect_absolute="
        << endpoint_negative_balance.local_update_energy_defect_absolute
        << "\nendpoint_negative_local_update_energy_ulp_bound="
        << endpoint_negative_balance.local_update_energy_ulp_bound
        << "\nendpoint_negative_balance_roundoff_bound="
        << endpoint_negative_balance.balance_roundoff_bound
        << "\nendpoint_negative_balance_roundoff_limited="
        << endpoint_negative_balance.balance_roundoff_limited
        << "\nendpoint_negative_r_fv_match_difference_absolute="
        << endpoint_negative_balance.r_fv_match_difference_absolute
        << "\nendpoint_negative_r_fv_match_roundoff_bound="
        << endpoint_negative_balance.r_fv_match_roundoff_bound
        << "\nendpoint_negative_r_fv_match_roundoff_limited="
        << endpoint_negative_balance.r_fv_match_roundoff_limited
        << "\nendpoint_negative_spatial_boundary_relative="
            << endpoint_negative_balance.spatial_boundary_relative
            << "\nconservation_pass=" << conservation_ok
            << "\npair_improvement_pass=" << pair_ok
            << "\ndual_decomposition_applicable="
            << dual_decomposition_applicable
            << "\ndual_decomposition_pass=" << dual_decomposition_ok
            << "\nzero_field_pass=" << zero_field_ok
            << "\nendpoint_ledger_pass=" << endpoint_ledger_ok
            << "\nendpoint_production_reconstruction_pass="
            << endpoint_production_reconstruction_ok
            << "\nendpoint_delta_ke_match_pass=" << endpoint_delta_ke_match_ok
            << "\nendpoint_u_work_match_pass=" << endpoint_u_work_match_ok
            << "\nendpoint_spatial_boundary_pass="
            << endpoint_spatial_boundary_ok
            << "\nendpoint_balance_pass=" << endpoint_balance_ok
            << "\nnonuniform_operator_active_pass=" << nonuniform_active_ok
            << "\npasses=" << passes << "\n";
        std::cout << "dual_u_stage2 case=" << case_name()
                  << " result=" << outdir << "/summary.result"
                  << " passes=" << passes << std::endl;
    }
    MPI_Finalize();
    return passes ? 0 : 1;
}
