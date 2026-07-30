#include "background_coupling_test_support.h"

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

namespace {

typedef VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle Bundle;

size_t uface_index(int ix, int jf, int k)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1) + jf) * Param::Nmu + k;
}

size_t xface_index(int iface, int j, int k)
{
    return (static_cast<size_t>(iface) * Param::Nv + j) * Param::Nmu + k;
}

bool make_directory(const std::string& path)
{
    if (path.empty() || path == ".") return true;
    std::string current;
    for (size_t i = 0; i < path.size(); ++i) {
        const char c = path[i];
        current += c;
        if (c != '/' && c != '\\' && i + 1 != path.size()) continue;
        while (!current.empty() &&
               (current[current.size() - 1] == '/' ||
                current[current.size() - 1] == '\\'))
            current.resize(current.size() - 1);
        if (current.empty()) continue;
#ifdef _WIN32
        const int rc = _mkdir(current.c_str());
#else
        const int rc = mkdir(current.c_str(), 0775);
#endif
        if (rc != 0 && errno != EEXIST) return false;
        current += '/';
    }
#ifdef _WIN32
    const int rc = _mkdir(path.c_str());
#else
    const int rc = mkdir(path.c_str(), 0775);
#endif
    return rc == 0 || errno == EEXIST;
}

double global_max_difference(const std::vector<double>& left,
                             const std::vector<double>& right)
{
    double local = left.size() == right.size() ? 0.0 :
        std::numeric_limits<double>::infinity();
    if (left.size() == right.size())
        for (size_t p = 0; p < left.size(); ++p)
            local = std::max(local, std::fabs(left[p] - right[p]));
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

bool valid_bundle(const Bundle& bundle, const SpatialGrid& sg)
{
    const size_t xsize = static_cast<size_t>(sg.nx_local + 1) * Param::Nvmu;
    const size_t usize = static_cast<size_t>(sg.nx_local) *
                         (Param::Nv + 1) * Param::Nmu;
    return bundle.state_advanced && !bundle.operator_failed &&
        bundle.outputs_finite && bundle.dual_u_operator_valid &&
        bundle.fx_low.size() == xsize && bundle.fx_high.size() == xsize &&
        bundle.fx_final.size() == xsize && bundle.fu_low.size() == usize &&
        bundle.fu_high.size() == usize && bundle.fu_final.size() == usize &&
        bundle.fu_fct_limited.size() == usize &&
        bundle.cu_low.size() == usize && bundle.cu_high.size() == usize &&
        bundle.cu_final.size() == usize &&
        bundle.cu_fct_limited.size() == usize &&
        bundle.donor_beta.size() == static_cast<size_t>(sg.nx_local) *
            Param::Nvmu &&
        bundle.donor_low_mass.size() == static_cast<size_t>(sg.nx_local) *
            Param::Nvmu &&
        bundle.donor_limited_outflow.size() == static_cast<size_t>(sg.nx_local) *
            Param::Nvmu;
}

struct AlphaAudit {
    double high_contract_linf;
    double high_contract_scale;
    double final_contract_linf;
    double final_contract_scale;
    double shared_alpha_linf;
    double alpha_range_violation;
    long long active_u_faces;
    long long comparable_u_faces;

    AlphaAudit()
        : high_contract_linf(0.0), high_contract_scale(0.0),
          final_contract_linf(0.0), final_contract_scale(0.0),
          shared_alpha_linf(0.0), alpha_range_violation(0.0),
          active_u_faces(0), comparable_u_faces(0) {}
};

enum VelocityRegion {
    THERMAL_BODY = 0,
    NEAR_UNIT_U = 1,
    LOW_DENSITY_TAIL = 2,
    VELOCITY_REGION_COUNT = 3
};

enum XRegion {
    CORE_REGION = 0,
    PERIODIC_SEAM_REGION = 1,
    X_REGION_COUNT = 2
};

struct LimiterBin {
    long long cells;
    long long active_cells;
    long long donor_beta_cells;
    long long physical_donor_beta_cells;
    long long beta_zero_cells;
    long double mass_total;
    long double mass_active;
    long double current_total;
    long double current_active;
    long double energy_total;
    long double energy_active;

    LimiterBin()
        : cells(0), active_cells(0), donor_beta_cells(0),
          physical_donor_beta_cells(0), beta_zero_cells(0),
          mass_total(0.0L), mass_active(0.0L), current_total(0.0L),
          current_active(0.0L), energy_total(0.0L), energy_active(0.0L) {}
};

struct BetaZeroRecord {
    int global_ix;
    int j;
    int k;
    double upar;
    double uperp;
    double beta;
    double m_low;
    double outflow;
};

const char* velocity_region_name(int region)
{
    switch (region) {
    case THERMAL_BODY: return "thermal_body";
    case NEAR_UNIT_U: return "near_abs_u_eq_1";
    default: return "low_density_tail";
    }
}

const char* x_region_name(int region)
{
    return region == PERIODIC_SEAM_REGION ? "periodic_seam" : "core";
}

bool cell_has_limiter(const Bundle& bundle, int ix, int j, int k)
{
    const double tolerance = 1024.0 * std::numeric_limits<double>::epsilon();
    const size_t p = static_cast<size_t>(ix) * Param::Nvmu +
        static_cast<size_t>(j) * Param::Nmu + k;
    const size_t xfaces[2] = {xface_index(ix, j, k),
                              xface_index(ix + 1, j, k)};
    const size_t ufaces[2] = {uface_index(ix, j, k),
                              uface_index(ix, j + 1, k)};
    for (int f = 0; f < 2; ++f) {
        const double xdelta = bundle.fx_high[xfaces[f]] - bundle.fx_low[xfaces[f]];
        if (std::fabs(xdelta) > tolerance * std::max(1.0,
                std::max(std::fabs(bundle.fx_high[xfaces[f]]),
                         std::fabs(bundle.fx_low[xfaces[f]]))) &&
            std::fabs(bundle.fx_final[xfaces[f]] - bundle.fx_high[xfaces[f]]) >
                tolerance * std::max(1.0, std::fabs(xdelta)))
            return true;
        const std::vector<double>& limited_u_flux =
            bundle.fu_fct_limited.size() == bundle.fu_final.size()
            ? bundle.fu_fct_limited : bundle.fu_final;
        const double udelta = bundle.fu_high[ufaces[f]] - bundle.fu_low[ufaces[f]];
        if (std::fabs(udelta) > tolerance * std::max(1.0,
                std::max(std::fabs(bundle.fu_high[ufaces[f]]),
                         std::fabs(bundle.fu_low[ufaces[f]]))) &&
            std::fabs(limited_u_flux[ufaces[f]] -
                      bundle.fu_high[ufaces[f]]) >
                tolerance * std::max(1.0, std::fabs(udelta)))
            return true;
    }
    (void)p;
    return false;
}

struct LimiterDistributionAudit {
    LimiterBin bins[X_REGION_COUNT][VELOCITY_REGION_COUNT];
    long long beta_zero_total;
    long long beta_zero_records_written;
    long long beta_zero_records_truncated;
    long long beta_zero_empty_low_mass;
    long long beta_zero_positive_low_mass;
    long long beta_zero_positive_outflow;
    std::vector<BetaZeroRecord> local_zero_records;

    LimiterDistributionAudit()
        : beta_zero_total(0), beta_zero_records_written(0),
          beta_zero_records_truncated(0), beta_zero_empty_low_mass(0),
          beta_zero_positive_low_mass(0), beta_zero_positive_outflow(0) {}
};

LimiterDistributionAudit audit_limiter_distribution(
    const Bundle& bundle, const Species& background, const SpatialGrid& sg)
{
    LimiterDistributionAudit audit;
    const size_t expected = static_cast<size_t>(sg.nx_local) * Param::Nvmu;
    if (bundle.donor_beta.size() != expected ||
        bundle.donor_low_mass.size() != expected ||
        bundle.donor_limited_outflow.size() != expected)
        return audit;

    double local_peak = 0.0;
    for (size_t p = 0; p < expected; ++p) {
        const int j = static_cast<int>((p % Param::Nvmu) / Param::Nmu);
        const int k = static_cast<int>(p % Param::Nmu);
        const double volume = sg.dx * background.cgrid.cell_phase_volume(j, k);
        if (volume > 0.0)
            local_peak = std::max(local_peak,
                std::fabs(bundle.donor_low_mass[p]) / volume);
    }
    double global_peak = 0.0;
    MPI_Allreduce(&local_peak, &global_peak, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const double tail_floor = std::max(std::numeric_limits<double>::min(),
        global_peak * 1.0e-8);
    const int seam_cells = std::max(1, static_cast<int>(std::ceil(
        0.1 * Const::micro / sg.dx)));
    const size_t per_rank_record_cap = 2048;

    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const int global_ix = sg.ix_start + ix;
        const int xr = (global_ix < seam_cells ||
                        global_ix >= sg.nx_global - seam_cells)
            ? PERIODIC_SEAM_REGION : CORE_REGION;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t p = static_cast<size_t>(ix) * Param::Nvmu +
                    static_cast<size_t>(j) * Param::Nmu + k;
                const double m = std::max(0.0, bundle.donor_low_mass[p]);
                const double volume = sg.dx * background.cgrid.cell_phase_volume(j, k);
                const double fbar = volume > 0.0 ? m / volume : 0.0;
                const double upar = background.cgrid.upar_cells[j];
                const double uperp = background.cgrid.uperp_cells[k];
                const double umag = std::sqrt(upar * upar + uperp * uperp);
                const double unit_half_width = std::max(0.10,
                    std::max(background.cgrid.upar_widths[j],
                             background.cgrid.uperp_widths[k]));
                const int vr = std::fabs(umag - 1.0) <= unit_half_width
                    ? NEAR_UNIT_U
                    : (fbar >= tail_floor ? THERMAL_BODY : LOW_DENSITY_TAIL);
                LimiterBin& bin = audit.bins[xr][vr];
                const bool active = cell_has_limiter(bundle, ix, j, k);
                const double beta = bundle.donor_beta[p];
                const double gamma = std::sqrt(1.0 + umag * umag);
                const double vx = Const::c * upar / gamma;
                const double energy = background.cgrid.kinetic_energy[
                    static_cast<size_t>(j) * Param::Nmu + k];
                const long double current_weight = std::fabs(background.charge * vx) * m;
                const long double energy_weight = energy * m;
                ++bin.cells;
                bin.mass_total += m;
                bin.current_total += current_weight;
                bin.energy_total += energy_weight;
                if (active) {
                    ++bin.active_cells;
                    bin.mass_active += m;
                    bin.current_active += current_weight;
                    bin.energy_active += energy_weight;
                }
                if (beta < 1.0 - 1.0e-14) {
                    ++bin.donor_beta_cells;
                    // A beta-zero empty cell is a numerical tail point with
                    // no donor mass.  Preserve it in the raw point count,
                    // but do not classify it as physical donor degradation.
                    if (bundle.donor_low_mass[p] > 0.0)
                        ++bin.physical_donor_beta_cells;
                }
                if (beta <= 0.0) {
                    ++bin.beta_zero_cells;
                    ++audit.beta_zero_total;
                    if (bundle.donor_low_mass[p] > 0.0)
                        ++audit.beta_zero_positive_low_mass;
                    else
                        ++audit.beta_zero_empty_low_mass;
                    if (bundle.donor_limited_outflow[p] > 0.0)
                        ++audit.beta_zero_positive_outflow;
                    if (audit.local_zero_records.size() < per_rank_record_cap) {
                        BetaZeroRecord record = {global_ix, j, k, upar, uperp,
                            beta, bundle.donor_low_mass[p],
                            bundle.donor_limited_outflow[p]};
                        audit.local_zero_records.push_back(record);
                    } else {
                        ++audit.beta_zero_records_truncated;
                    }
                }
            }
        }
    }
    return audit;
}

void reduce_limiter_distribution(LimiterDistributionAudit& audit)
{
    long long local_counts[X_REGION_COUNT * VELOCITY_REGION_COUNT * 5];
    long double local_weights[X_REGION_COUNT * VELOCITY_REGION_COUNT * 6];
    size_t ci = 0;
    size_t wi = 0;
    for (int xr = 0; xr < X_REGION_COUNT; ++xr) {
        for (int vr = 0; vr < VELOCITY_REGION_COUNT; ++vr) {
            const LimiterBin& bin = audit.bins[xr][vr];
            local_counts[ci++] = bin.cells;
            local_counts[ci++] = bin.active_cells;
            local_counts[ci++] = bin.donor_beta_cells;
            local_counts[ci++] = bin.physical_donor_beta_cells;
            local_counts[ci++] = bin.beta_zero_cells;
            local_weights[wi++] = bin.mass_total;
            local_weights[wi++] = bin.mass_active;
            local_weights[wi++] = bin.current_total;
            local_weights[wi++] = bin.current_active;
            local_weights[wi++] = bin.energy_total;
            local_weights[wi++] = bin.energy_active;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, local_counts,
                  X_REGION_COUNT * VELOCITY_REGION_COUNT * 5,
                  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, local_weights,
                  X_REGION_COUNT * VELOCITY_REGION_COUNT * 6,
                  MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ci = 0;
    wi = 0;
    for (int xr = 0; xr < X_REGION_COUNT; ++xr) {
        for (int vr = 0; vr < VELOCITY_REGION_COUNT; ++vr) {
            LimiterBin& bin = audit.bins[xr][vr];
            bin.cells = local_counts[ci++];
            bin.active_cells = local_counts[ci++];
            bin.donor_beta_cells = local_counts[ci++];
            bin.physical_donor_beta_cells = local_counts[ci++];
            bin.beta_zero_cells = local_counts[ci++];
            bin.mass_total = local_weights[wi++];
            bin.mass_active = local_weights[wi++];
            bin.current_total = local_weights[wi++];
            bin.current_active = local_weights[wi++];
            bin.energy_total = local_weights[wi++];
            bin.energy_active = local_weights[wi++];
        }
    }
    long long totals[6] = {audit.beta_zero_total,
        static_cast<long long>(audit.local_zero_records.size()),
        audit.beta_zero_records_truncated, audit.beta_zero_empty_low_mass,
        audit.beta_zero_positive_low_mass, audit.beta_zero_positive_outflow};
    MPI_Allreduce(MPI_IN_PLACE, totals, 6, MPI_LONG_LONG_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    audit.beta_zero_total = totals[0];
    audit.beta_zero_records_written = totals[1];
    audit.beta_zero_records_truncated = totals[2];
    audit.beta_zero_empty_low_mass = totals[3];
    audit.beta_zero_positive_low_mass = totals[4];
    audit.beta_zero_positive_outflow = totals[5];
}

void write_beta_zero_records(const std::string& output_dir,
                             const LimiterDistributionAudit& audit,
                             int rank, int size, const char* scenario)
{
    std::ostringstream local;
    local << std::scientific << std::setprecision(17);
    for (size_t i = 0; i < audit.local_zero_records.size(); ++i) {
        const BetaZeroRecord& r = audit.local_zero_records[i];
        local << scenario << " " << r.global_ix << " " << r.j << " " << r.k
              << " " << r.upar << " " << r.uperp << " " << r.beta
              << " " << r.m_low << " " << r.outflow << " "
              << (r.outflow > 0.0 ? r.m_low / r.outflow : 0.0) << "\n";
    }
    const std::string text = local.str();
    const int local_size = static_cast<int>(text.size());
    std::vector<int> sizes(rank == 0 ? size : 0, 0);
    MPI_Gather(&local_size, 1, MPI_INT, rank == 0 ? sizes.data() : 0, 1,
               MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displacements;
    std::vector<char> gathered;
    if (rank == 0) {
        displacements.resize(size, 0);
        int total = 0;
        for (int r = 0; r < size; ++r) {
            displacements[r] = total;
            total += sizes[r];
        }
        gathered.resize(static_cast<size_t>(total));
    }
    MPI_Gatherv(text.empty() ? 0 : const_cast<char*>(text.data()), local_size,
                MPI_CHAR, rank == 0 ? gathered.data() : 0,
                rank == 0 ? sizes.data() : 0,
                rank == 0 ? displacements.data() : 0, MPI_CHAR, 0,
                MPI_COMM_WORLD);
    if (rank == 0 && !gathered.empty()) {
        std::ofstream out((output_dir + "/limiter_beta_zero_cells.dat").c_str(),
                          std::ios::app);
        out.write(gathered.data(), static_cast<std::streamsize>(gathered.size()));
    }
}

void write_limiter_distribution(const std::string& output_dir,
                                const LimiterDistributionAudit& normal,
                                const LimiterDistributionAudit& controlled)
{
    std::ofstream out((output_dir + "/limiter_distribution.dat").c_str());
    out << "# scenario x_region velocity_region cells active_cells donor_beta_cells"
        << " physical_donor_beta_cells"
        << " beta_zero_cells mass_total mass_active mass_coverage"
        << " current_abs_total current_abs_active current_coverage"
        << " energy_total energy_active energy_coverage\n";
    const LimiterDistributionAudit* audits[2] = {&normal, &controlled};
    const char* names[2] = {"normal", "controlled"};
    out << std::scientific << std::setprecision(17);
    for (int scenario = 0; scenario < 2; ++scenario) {
        for (int xr = 0; xr < X_REGION_COUNT; ++xr) {
            for (int vr = 0; vr < VELOCITY_REGION_COUNT; ++vr) {
                const LimiterBin& bin = audits[scenario]->bins[xr][vr];
                const long double mass_coverage = bin.mass_active /
                    std::max(1.0L, bin.mass_total);
                const long double current_coverage = bin.current_active /
                    std::max(1.0L, bin.current_total);
                const long double energy_coverage = bin.energy_active /
                    std::max(1.0L, bin.energy_total);
                out << names[scenario] << " " << x_region_name(xr) << " "
                    << velocity_region_name(vr) << " " << bin.cells << " "
                    << bin.active_cells << " " << bin.donor_beta_cells << " "
                    << bin.physical_donor_beta_cells << " "
                    << bin.beta_zero_cells << " "
                    << static_cast<double>(bin.mass_total) << " "
                    << static_cast<double>(bin.mass_active) << " "
                    << static_cast<double>(mass_coverage) << " "
                    << static_cast<double>(bin.current_total) << " "
                    << static_cast<double>(bin.current_active) << " "
                    << static_cast<double>(current_coverage) << " "
                    << static_cast<double>(bin.energy_total) << " "
                    << static_cast<double>(bin.energy_active) << " "
                    << static_cast<double>(energy_coverage) << "\n";
            }
        }
    }
}

struct LimiterAcceptanceMetrics {
    long long unique_donor_cells;
    long long physical_unique_donor_cells;
    long long thermal_cells;
    long long thermal_physical_donor_cells;
    long double mass_total;
    long double mass_active;
    long double current_total;
    long double current_active;
    long double energy_total;
    long double energy_active;

    LimiterAcceptanceMetrics()
        : unique_donor_cells(0), physical_unique_donor_cells(0),
          thermal_cells(0), thermal_physical_donor_cells(0),
          mass_total(0.0L), mass_active(0.0L), current_total(0.0L),
          current_active(0.0L), energy_total(0.0L), energy_active(0.0L) {}
};

LimiterAcceptanceMetrics limiter_acceptance_metrics(
    const LimiterDistributionAudit& audit)
{
    LimiterAcceptanceMetrics metrics;
    for (int xr = 0; xr < X_REGION_COUNT; ++xr) {
        for (int vr = 0; vr < VELOCITY_REGION_COUNT; ++vr) {
            const LimiterBin& bin = audit.bins[xr][vr];
            metrics.unique_donor_cells += bin.donor_beta_cells;
            metrics.physical_unique_donor_cells +=
                bin.physical_donor_beta_cells;
            metrics.mass_total += bin.mass_total;
            metrics.mass_active += bin.mass_active;
            metrics.current_total += bin.current_total;
            metrics.current_active += bin.current_active;
            metrics.energy_total += bin.energy_total;
            metrics.energy_active += bin.energy_active;
            if (vr == THERMAL_BODY) {
                metrics.thermal_cells += bin.cells;
                metrics.thermal_physical_donor_cells +=
                    bin.physical_donor_beta_cells;
            }
        }
    }
    return metrics;
}

double normalized_coverage(long double active, long double total)
{
    return static_cast<double>(active / std::max(1.0L, total));
}

AlphaAudit audit_u_alpha(const Bundle& bundle, const Species& background,
                         const EMFields& fields, const SpatialGrid& sg)
{
    AlphaAudit local;
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        const double acceleration = background.charge *
            fields.Ex[sg.nghost + ix] / (background.mass * Const::c);
        for (int jf = 1; jf < Param::Nv; ++jf) {
            for (int k = 0; k < Param::Nmu; ++k) {
                const size_t id = uface_index(ix, jf, k);
                local.high_contract_linf = std::max(
                    local.high_contract_linf,
                    std::fabs(bundle.fu_high[id] -
                              acceleration * bundle.cu_high[id]));
                local.high_contract_scale = std::max(
                    local.high_contract_scale,
                    std::max(std::fabs(bundle.fu_high[id]),
                             std::fabs(acceleration * bundle.cu_high[id])));
                local.final_contract_linf = std::max(
                    local.final_contract_linf,
                    std::fabs(bundle.fu_final[id] -
                              acceleration * bundle.cu_final[id]));
                local.final_contract_scale = std::max(
                    local.final_contract_scale,
                    std::max(std::fabs(bundle.fu_final[id]),
                             std::fabs(acceleration * bundle.cu_final[id])));
                const double df = bundle.fu_high[id] - bundle.fu_low[id];
                const double dc = bundle.cu_high[id] - bundle.cu_low[id];
                const double fscale = std::max(1.0,
                    std::max(std::fabs(bundle.fu_high[id]),
                             std::fabs(bundle.fu_low[id])));
                const double cscale = std::max(1.0,
                    std::max(std::fabs(bundle.cu_high[id]),
                             std::fabs(bundle.cu_low[id])));
                if (std::fabs(df) <= 128.0 *
                        std::numeric_limits<double>::epsilon() * fscale ||
                    std::fabs(dc) <= 128.0 *
                        std::numeric_limits<double>::epsilon() * cscale)
                    continue;
                const double alpha_f = (bundle.fu_fct_limited[id] -
                    bundle.fu_low[id]) / df;
                const double alpha_c = (bundle.cu_fct_limited[id] -
                    bundle.cu_low[id]) / dc;
                local.shared_alpha_linf = std::max(local.shared_alpha_linf,
                    std::fabs(alpha_f - alpha_c));
                local.alpha_range_violation = std::max(
                    local.alpha_range_violation,
                    std::max(0.0, std::max(-alpha_f, alpha_f - 1.0)));
                ++local.comparable_u_faces;
                if (alpha_f < 1.0 - 1.0e-12) ++local.active_u_faces;
            }
        }
    }
    double maxima[6] = {local.high_contract_linf,
        local.high_contract_scale, local.final_contract_linf,
        local.final_contract_scale, local.shared_alpha_linf,
        local.alpha_range_violation};
    MPI_Allreduce(MPI_IN_PLACE, maxima, 6, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    local.high_contract_linf = maxima[0];
    local.high_contract_scale = maxima[1];
    local.final_contract_linf = maxima[2];
    local.final_contract_scale = maxima[3];
    local.shared_alpha_linf = maxima[4];
    local.alpha_range_violation = maxima[5];
    long long counts[2] = {local.active_u_faces, local.comparable_u_faces};
    MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_LONG_LONG_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    local.active_u_faces = counts[0];
    local.comparable_u_faces = counts[1];
    return local;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string output_dir = "output/dual_u_stage3/muscl_fct";
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--output-dir" && i + 1 < argc)
            output_dir = argv[++i];
    }
    int directory_ok = rank == 0 && make_directory(output_dir) ? 1 : 0;
    MPI_Bcast(&directory_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!directory_ok) {
        if (rank == 0) std::cerr << "cannot create " << output_dir << '\n';
        MPI_Finalize();
        return 2;
    }

    SpatialGrid sg;
    sg.init(rank, size);
    Species background;
    EMFields fields;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.08, 0.10, 2.0e3, 5, 0.21);
    const double dt = BackgroundCouplingTest::stable_dt(sg);

    BackgroundCouplingTest::BundleOptions legacy_high_options;
    legacy_high_options.fct_enabled = false;
    legacy_high_options.allow_finite_negative_debt = true;
    legacy_high_options.coupling_mode =
        VlasovAmpereMidpointSolver::LEGACY_COUPLING;
    BackgroundCouplingTest::BundleOptions dual_high_options =
        legacy_high_options;
    dual_high_options.coupling_mode =
        VlasovAmpereMidpointSolver::DUAL_U_COUPLING;

    BackgroundCouplingTest::BundleOptions legacy_fct_options;
    legacy_fct_options.fct_enabled = true;
    legacy_fct_options.coupling_mode =
        VlasovAmpereMidpointSolver::LEGACY_COUPLING;
    BackgroundCouplingTest::BundleOptions dual_fct_options =
        legacy_fct_options;
    dual_fct_options.coupling_mode =
        VlasovAmpereMidpointSolver::DUAL_U_COUPLING;
    BackgroundCouplingTest::BundleOptions controlled_options =
        dual_fct_options;
    controlled_options.controlled_u_fct_injection = true;

    const Bundle legacy_high = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, legacy_high_options);
    const Bundle dual_high = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, dual_high_options);
    const Bundle legacy_fct = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, legacy_fct_options);
    const Bundle dual_fct = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, dual_fct_options);
    const Bundle controlled_fct = BackgroundCouplingTest::evaluate_bundle(
        background, fields, sg, rank, size, dt, controlled_options);

    // These scans run only in the stage-3 fixed-candidate test.  They use the
    // exact production FCT outputs and do not feed back into the solver.
    LimiterDistributionAudit normal_distribution =
        audit_limiter_distribution(dual_fct, background, sg);
    LimiterDistributionAudit controlled_distribution =
        audit_limiter_distribution(controlled_fct, background, sg);
    reduce_limiter_distribution(normal_distribution);
    reduce_limiter_distribution(controlled_distribution);
    const auto donor_cells_in_velocity_region = [](const LimiterDistributionAudit& audit,
                                                     int velocity_region,
                                                     bool physical_only) {
        long long count = 0;
        for (int xr = 0; xr < X_REGION_COUNT; ++xr)
            count += physical_only
                ? audit.bins[xr][velocity_region].physical_donor_beta_cells
                : audit.bins[xr][velocity_region].donor_beta_cells;
        return count;
    };
    const long long normal_donor_tail_count =
        donor_cells_in_velocity_region(normal_distribution, LOW_DENSITY_TAIL,
                                       false);
    const long long controlled_donor_tail_count =
        donor_cells_in_velocity_region(controlled_distribution, LOW_DENSITY_TAIL,
                                       false);
    const long long normal_physical_donor_tail_count =
        donor_cells_in_velocity_region(normal_distribution, LOW_DENSITY_TAIL,
                                       true);
    const long long controlled_physical_donor_tail_count =
        donor_cells_in_velocity_region(controlled_distribution,
                                       LOW_DENSITY_TAIL, true);
    const LimiterAcceptanceMetrics normal_metrics =
        limiter_acceptance_metrics(normal_distribution);
    const LimiterAcceptanceMetrics controlled_metrics =
        limiter_acceptance_metrics(controlled_distribution);
    if (rank == 0) {
        std::ofstream beta_zero((output_dir + "/limiter_beta_zero_cells.dat").c_str());
        beta_zero << "# scenario global_ix iv iuperp upar uperp beta m_low"
                  << " limited_outflow m_low_over_outflow\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    write_beta_zero_records(output_dir, normal_distribution, rank, size,
                            "normal");
    write_beta_zero_records(output_dir, controlled_distribution, rank, size,
                            "controlled");
    if (rank == 0)
        write_limiter_distribution(output_dir, normal_distribution,
                                   controlled_distribution);

    const bool contracts = valid_bundle(dual_high, sg) &&
        valid_bundle(dual_fct, sg) && valid_bundle(controlled_fct, sg) &&
        legacy_high.state_advanced && !legacy_high.operator_failed &&
        legacy_fct.state_advanced && !legacy_fct.operator_failed;

    const double x_high_difference = global_max_difference(
        legacy_high.fx_high, dual_high.fx_high);
    const double jn_high_difference = global_max_difference(
        legacy_high.jn_high, dual_high.jn_high);
    const double low_x_difference = global_max_difference(
        legacy_fct.fx_low, dual_fct.fx_low);
    const double low_u_difference = global_max_difference(
        legacy_fct.fu_low, dual_fct.fu_low);
    const double low_cu_difference = global_max_difference(
        legacy_fct.cu_low, dual_fct.cu_low);
    const double target_replay_relative = dual_high.dual_u_target_replay_linf /
        std::max(1.0, dual_high.dual_u_target_replay_scale);
    const double operator_replay_relative =
        dual_high.dual_u_legacy_operator_replay_linf /
        std::max(1.0, dual_high.dual_u_legacy_operator_replay_scale);

    const BackgroundCouplingTest::Norms legacy_high_pair =
        BackgroundCouplingTest::face_difference_norms(
            legacy_high.jn_high, legacy_high.gstar_je_high, sg);
    const BackgroundCouplingTest::Norms dual_high_pair =
        BackgroundCouplingTest::face_difference_norms(
            dual_high.jn_high, dual_high.gstar_je_high, sg);
    const BackgroundCouplingTest::Norms legacy_final_pair =
        BackgroundCouplingTest::face_difference_norms(
            legacy_fct.jn_final, legacy_fct.gstar_je_final, sg);
    const BackgroundCouplingTest::Norms dual_final_pair =
        BackgroundCouplingTest::face_difference_norms(
            dual_fct.jn_final, dual_fct.gstar_je_final, sg);

    const AlphaAudit alpha = audit_u_alpha(
        controlled_fct, background, fields, sg);
    const double high_contract_relative = alpha.high_contract_linf /
        std::max(1.0, alpha.high_contract_scale);
    const double final_contract_relative = alpha.final_contract_linf /
        std::max(1.0, alpha.final_contract_scale);
    const double donor_degradation_fraction =
        static_cast<double>(controlled_fct.donor_beta_applied_count) /
        std::max(1.0, static_cast<double>(sg.nx_global) * Param::Nvmu);
    const double normal_donor_degradation_fraction =
        static_cast<double>(dual_fct.donor_beta_applied_count) /
        std::max(1.0, static_cast<double>(sg.nx_global) * Param::Nvmu);
    const double normal_thermal_physical_donor_coverage =
        static_cast<double>(normal_metrics.thermal_physical_donor_cells) /
        std::max(1LL, normal_metrics.thermal_cells);
    const double normal_mass_weighted_limiter_coverage = normalized_coverage(
        normal_metrics.mass_active, normal_metrics.mass_total);
    const double normal_current_weighted_limiter_coverage = normalized_coverage(
        normal_metrics.current_active, normal_metrics.current_total);
    const double normal_energy_weighted_limiter_coverage = normalized_coverage(
        normal_metrics.energy_active, normal_metrics.energy_total);
    const long long normal_trigger_minus_unique =
        dual_fct.donor_beta_applied_count - normal_metrics.unique_donor_cells;
    const bool normal_donor_count_semantics_ok =
        dual_fct.donor_beta_applied_count >= normal_metrics.unique_donor_cells;
    const double controlled_u_active_fraction =
        alpha.comparable_u_faces > 0 ?
        static_cast<double>(alpha.active_u_faces) /
            alpha.comparable_u_faces : 0.0;
    const double only_final_growth = dual_final_pair.l2 /
        std::max(dual_high_pair.l2, std::numeric_limits<double>::min());

    const bool high_ok = x_high_difference == 0.0 &&
        jn_high_difference == 0.0 &&
        dual_high_pair.l2 <= legacy_high_pair.l2 * (1.0 + 1.0e-12);
    const bool muscl_ok = target_replay_relative <= 1.0e-11 &&
        operator_replay_relative <= 1.0e-11 &&
        dual_high.dual_u_corrected_cell_count > 0;
    const bool low_ok = low_x_difference == 0.0 && low_u_difference == 0.0 &&
        low_cu_difference == 0.0;
    const bool fct_ok = controlled_fct.u_limiter_active_fraction > 0.0 &&
        controlled_fct.u_limiter_min_alpha < 1.0 &&
        alpha.active_u_faces > 0 &&
        alpha.shared_alpha_linf <= 1.0e-10 &&
        alpha.alpha_range_violation <= 1.0e-12;
    const bool final_moment_ok =
        controlled_fct.final_flux_current_moment_audit_valid &&
        controlled_fct.final_flux_current_moment_audit_finite &&
        controlled_fct.final_flux_to_jn_linf <= 1.0e-10 *
            std::max(1.0, controlled_fct.final_flux_to_jn_scale) &&
        controlled_fct.final_flux_to_je_linf <= 1.0e-10 *
            std::max(1.0, controlled_fct.final_flux_to_je_scale) &&
        high_contract_relative <= 4096.0 *
            std::numeric_limits<double>::epsilon() &&
        final_contract_relative <= 4096.0 *
            std::numeric_limits<double>::epsilon();
    // Point coverage remains a diagnostic: empty tail cells can be numerous
    // while carrying no mass, current, or energy.  The acceptance gate uses
    // the same 1e-2 scale for physical thermal donors and weighted impact.
    const double limiter_coverage_limit = 1.0e-2;
    const bool limiter_rate_ok = normal_donor_count_semantics_ok &&
        normal_thermal_physical_donor_coverage <= limiter_coverage_limit &&
        normal_mass_weighted_limiter_coverage <= limiter_coverage_limit &&
        normal_current_weighted_limiter_coverage <= limiter_coverage_limit &&
        normal_energy_weighted_limiter_coverage <= limiter_coverage_limit;
    const bool final_pair_ok = dual_final_pair.l2 <=
        legacy_final_pair.l2 * (1.0 + 1.0e-12);

    int passes = contracts && high_ok && muscl_ok && low_ok && fct_ok &&
        final_moment_ok && limiter_rate_ok && final_pair_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &passes, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream out((output_dir + "/summary.result").c_str());
        out << std::scientific << std::setprecision(17)
            << "test=background_coupling_dual_u_stage3\n"
            << "mpi_size=" << size << "\n"
            << "Nx=" << sg.nx_global << "\nNupar=" << Param::Nv
            << "\nNuperp=" << Param::Nmu << "\n"
            << "bundle_contract_pass=" << contracts << "\n"
            << "center_high_no_limiter_pass=" << high_ok << "\n"
            << "x_high_legacy_dual_linf=" << x_high_difference << "\n"
            << "jn_high_legacy_dual_linf=" << jn_high_difference << "\n"
            << "legacy_high_pair_L2=" << legacy_high_pair.l2 << "\n"
            << "dual_high_pair_L2=" << dual_high_pair.l2 << "\n"
            << "muscl_fixed_branch_pass=" << muscl_ok << "\n"
            << "target_replay_relative=" << target_replay_relative << "\n"
            << "operator_replay_relative=" << operator_replay_relative << "\n"
            << "dual_corrected_cell_count="
            << dual_high.dual_u_corrected_cell_count << "\n"
            << "low_order_base_pass=" << low_ok << "\n"
            << "fx_low_legacy_dual_linf=" << low_x_difference << "\n"
            << "fu_low_legacy_dual_linf=" << low_u_difference << "\n"
            << "cu_low_legacy_dual_linf=" << low_cu_difference << "\n"
            << "fct_shared_alpha_pass=" << fct_ok << "\n"
            << "controlled_u_limiter_active_fraction="
            << controlled_fct.u_limiter_active_fraction << "\n"
            << "controlled_u_limiter_min_alpha="
            << controlled_fct.u_limiter_min_alpha << "\n"
            << "shared_alpha_linf=" << alpha.shared_alpha_linf << "\n"
            << "shared_alpha_range_violation="
            << alpha.alpha_range_violation << "\n"
            << "shared_alpha_active_u_faces=" << alpha.active_u_faces << "\n"
            << "shared_alpha_comparable_u_faces="
            << alpha.comparable_u_faces << "\n"
            << "shared_alpha_active_fraction="
            << controlled_u_active_fraction << "\n"
            << "final_flux_moment_pass=" << final_moment_ok << "\n"
            << "high_fu_equals_a_cu_linf=" << alpha.high_contract_linf << "\n"
            << "high_fu_equals_a_cu_relative="
            << high_contract_relative << "\n"
            << "final_fu_equals_a_cu_linf=" << alpha.final_contract_linf << "\n"
            << "final_fu_equals_a_cu_relative="
            << final_contract_relative << "\n"
            << "final_flux_to_jn_linf="
            << controlled_fct.final_flux_to_jn_linf << "\n"
            << "final_flux_to_je_linf="
            << controlled_fct.final_flux_to_je_linf << "\n"
            << "limiter_rate_pass=" << limiter_rate_ok << "\n"
            << "normal_x_limiter_active_fraction="
            << dual_fct.x_limiter_active_fraction << "\n"
            << "normal_u_limiter_active_fraction="
            << dual_fct.u_limiter_active_fraction << "\n"
            << "normal_limiter_active_fraction="
            << dual_fct.limiter_active_fraction << "\n"
            << "normal_limiter_point_coverage_diagnostic="
            << dual_fct.limiter_active_fraction << "\n"
            << "normal_thermal_physical_donor_coverage="
            << normal_thermal_physical_donor_coverage << "\n"
            << "normal_mass_weighted_limiter_coverage="
            << normal_mass_weighted_limiter_coverage << "\n"
            << "normal_current_weighted_limiter_coverage="
            << normal_current_weighted_limiter_coverage << "\n"
            << "normal_energy_weighted_limiter_coverage="
            << normal_energy_weighted_limiter_coverage << "\n"
            << "limiter_weighted_coverage_limit="
            << limiter_coverage_limit << "\n"
            << "normal_donor_beta_cumulative_trigger_count="
            << dual_fct.donor_beta_applied_count << "\n"
            << "normal_donor_beta_unique_cell_count="
            << normal_metrics.unique_donor_cells << "\n"
            << "normal_physical_donor_beta_unique_cell_count="
            << normal_metrics.physical_unique_donor_cells << "\n"
            << "normal_donor_beta_trigger_minus_unique="
            << normal_trigger_minus_unique << "\n"
            << "normal_donor_beta_count_semantics_ok="
            << normal_donor_count_semantics_ok << "\n"
            << "normal_donor_beta_applied_count="
            << dual_fct.donor_beta_applied_count << "\n"
            << "normal_donor_beta_min=" << dual_fct.donor_beta_min << "\n"
            << "normal_donor_cell_degradation_fraction="
            << normal_donor_degradation_fraction << "\n"
            << "normal_donor_beta_tail_cell_count=" << normal_donor_tail_count
            << "\n"
            << "normal_donor_beta_tail_fraction="
            << static_cast<double>(normal_donor_tail_count) /
                std::max(1LL, normal_metrics.unique_donor_cells) << "\n"
            << "normal_physical_donor_beta_tail_cell_count="
            << normal_physical_donor_tail_count << "\n"
            << "normal_physical_donor_beta_tail_fraction="
            << static_cast<double>(normal_physical_donor_tail_count) /
                std::max(1LL, normal_metrics.physical_unique_donor_cells)
            << "\n"
            << "controlled_donor_beta_cumulative_trigger_count="
            << controlled_fct.donor_beta_applied_count << "\n"
            << "controlled_donor_beta_unique_cell_count="
            << controlled_metrics.unique_donor_cells << "\n"
            << "controlled_physical_donor_beta_unique_cell_count="
            << controlled_metrics.physical_unique_donor_cells << "\n"
            << "donor_beta_applied_count="
            << controlled_fct.donor_beta_applied_count << "\n"
            << "donor_beta_min=" << controlled_fct.donor_beta_min << "\n"
            << "donor_cell_degradation_fraction="
            << donor_degradation_fraction << "\n"
            << "controlled_donor_beta_tail_cell_count="
            << controlled_donor_tail_count << "\n"
            << "controlled_donor_beta_tail_fraction="
            << static_cast<double>(controlled_donor_tail_count) /
                std::max(1LL, controlled_metrics.unique_donor_cells) << "\n"
            << "controlled_physical_donor_beta_tail_cell_count="
            << controlled_physical_donor_tail_count << "\n"
            << "controlled_physical_donor_beta_tail_fraction="
            << static_cast<double>(controlled_physical_donor_tail_count) /
                std::max(1LL, controlled_metrics.physical_unique_donor_cells)
            << "\n"
            << "normal_beta_zero_total_count="
            << normal_distribution.beta_zero_total << "\n"
            << "normal_beta_zero_records_written="
            << normal_distribution.beta_zero_records_written << "\n"
            << "normal_beta_zero_records_truncated="
            << normal_distribution.beta_zero_records_truncated << "\n"
            << "normal_beta_zero_empty_low_mass_count="
            << normal_distribution.beta_zero_empty_low_mass << "\n"
            << "normal_beta_zero_positive_low_mass_count="
            << normal_distribution.beta_zero_positive_low_mass << "\n"
            << "normal_beta_zero_positive_outflow_count="
            << normal_distribution.beta_zero_positive_outflow << "\n"
            << "controlled_beta_zero_total_count="
            << controlled_distribution.beta_zero_total << "\n"
            << "controlled_beta_zero_records_written="
            << controlled_distribution.beta_zero_records_written << "\n"
            << "controlled_beta_zero_records_truncated="
            << controlled_distribution.beta_zero_records_truncated << "\n"
            << "controlled_beta_zero_empty_low_mass_count="
            << controlled_distribution.beta_zero_empty_low_mass << "\n"
            << "controlled_beta_zero_positive_low_mass_count="
            << controlled_distribution.beta_zero_positive_low_mass << "\n"
            << "controlled_beta_zero_positive_outflow_count="
            << controlled_distribution.beta_zero_positive_outflow << "\n"
            << "final_pair_pass=" << final_pair_ok << "\n"
            << "legacy_final_pair_L2=" << legacy_final_pair.l2 << "\n"
            << "dual_final_pair_L2=" << dual_final_pair.l2 << "\n"
            << "dual_only_final_pair_growth=" << only_final_growth << "\n"
            << "passes=" << passes << "\n";
        std::cout << "dual_u_stage3 result=" << output_dir
                  << "/summary.result passes=" << passes << std::endl;
    }

    MPI_Finalize();
    return passes ? 0 : 1;
}
