#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "species.h"
#include "tail_subcell_quadrature.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string checkpoint;
    std::string snapshot25;
    std::string snapshot40;
    std::string result;
    double kout_mev;
    double kin_mev;
    double project_fs;
    double core_noise_limit;
    std::vector<double> bin_widths_mev;
    Options()
        : kout_mev(6.0), kin_mev(0.0), project_fs(120.0),
          core_noise_limit(1.0e-3)
    {}
};

struct SnapshotStats {
    bool complete;
    double tail_number;
    double macro_count;
    SnapshotStats() : complete(false), tail_number(0.0), macro_count(0.0) {}
};

struct SpectrumStats {
    std::vector<double> values;
    std::vector<double> bulk_values;
    std::vector<double> tail_values;
    double total;
};

struct PerformanceStats {
    bool complete;
    double projected_wall_seconds;
    double projected_wall_seconds_per_step;
    double dt_fs;
    double fit_intercept_seconds;
    double fit_slope_seconds_per_particle;
    PerformanceStats()
        : complete(false),
          projected_wall_seconds(std::numeric_limits<double>::quiet_NaN()),
          projected_wall_seconds_per_step(
              std::numeric_limits<double>::quiet_NaN()),
          dt_fs(std::numeric_limits<double>::quiet_NaN()),
          fit_intercept_seconds(std::numeric_limits<double>::quiet_NaN()),
          fit_slope_seconds_per_particle(
              std::numeric_limits<double>::quiet_NaN())
    {}
};

struct CoreNoiseStats {
    bool complete;
    double relative_to_tail_envelope;
    double relative_to_background_rms;
    double tail_rms_relative_to_background_rms;
    double max_tail_relative_to_background;
    CoreNoiseStats()
        : complete(false),
          relative_to_tail_envelope(
              std::numeric_limits<double>::quiet_NaN()),
          relative_to_background_rms(
              std::numeric_limits<double>::quiet_NaN()),
          tail_rms_relative_to_background_rms(
              std::numeric_limits<double>::quiet_NaN()),
          max_tail_relative_to_background(
              std::numeric_limits<double>::quiet_NaN())
    {}
};

struct SpectrumGapStats {
    int deep_valley_count;
    double max_adjacent_jump_ratio;
    bool pass;
    SpectrumGapStats()
        : deep_valley_count(0), max_adjacent_jump_ratio(0.0), pass(false)
    {}
};

std::string parent_directory(const std::string& path)
{
    const std::string::size_type slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string(".") :
           path.substr(0, slash);
}

void globalize_spectrum(SpectrumStats& spectrum)
{
    if (!spectrum.values.empty()) {
        MPI_Allreduce(MPI_IN_PLACE, spectrum.values.data(),
                      static_cast<int>(spectrum.values.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, spectrum.bulk_values.data(),
                      static_cast<int>(spectrum.bulk_values.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, spectrum.tail_values.data(),
                      static_cast<int>(spectrum.tail_values.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
    }
    MPI_Allreduce(MPI_IN_PLACE, &spectrum.total, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
}

bool require_value(int& i, int argc, char**, const char* option)
{
    if (i + 1 < argc) return true;
    std::cerr << option << " requires a value\n";
    return false;
}

bool parse_double(const char* text, double& value)
{
    char* end = NULL;
    value = std::strtod(text, &end);
    return end != text && *end == '\0' && std::isfinite(value);
}

bool parse_widths(const std::string& text, std::vector<double>& widths)
{
    std::stringstream input(text);
    std::string token;
    widths.clear();
    while (std::getline(input, token, ',')) {
        double width = 0.0;
        if (!parse_double(token.c_str(), width) || !(width > 0.0))
            return false;
        widths.push_back(width);
    }
    if (widths.size() != 3) return false;
    std::sort(widths.begin(), widths.end());
    const double expected[3] = {0.1, 0.2, 0.4};
    for (size_t i = 0; i < widths.size(); ++i) {
        if (std::fabs(widths[i] - expected[i]) > 1.0e-12)
            return false;
    }
    return true;
}

bool parse_options(int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--checkpoint") {
            if (!require_value(i, argc, argv, "--checkpoint")) return false;
            options.checkpoint = argv[++i];
        } else if (arg == "--snapshot25") {
            if (!require_value(i, argc, argv, "--snapshot25")) return false;
            options.snapshot25 = argv[++i];
        } else if (arg == "--snapshot40") {
            if (!require_value(i, argc, argv, "--snapshot40")) return false;
            options.snapshot40 = argv[++i];
        } else if (arg == "--kout-mev") {
            if (!require_value(i, argc, argv, "--kout-mev") ||
                !parse_double(argv[++i], options.kout_mev)) return false;
        } else if (arg == "--kin-mev") {
            if (!require_value(i, argc, argv, "--kin-mev") ||
                !parse_double(argv[++i], options.kin_mev)) return false;
        } else if (arg == "--bin-widths-mev") {
            if (!require_value(i, argc, argv, "--bin-widths-mev") ||
                !parse_widths(argv[++i], options.bin_widths_mev)) return false;
        } else if (arg == "--project-to-fs") {
            if (!require_value(i, argc, argv, "--project-to-fs") ||
                !parse_double(argv[++i], options.project_fs)) return false;
        } else if (arg == "--core-noise-limit") {
            if (!require_value(i, argc, argv, "--core-noise-limit") ||
                !parse_double(argv[++i], options.core_noise_limit))
                return false;
        } else if (arg == "--result") {
            if (!require_value(i, argc, argv, "--result")) return false;
            options.result = argv[++i];
        } else {
            std::cerr << "unknown option: " << arg << "\n";
            return false;
        }
    }
    return !options.checkpoint.empty() && !options.snapshot25.empty() &&
           !options.snapshot40.empty() && !options.result.empty() &&
           options.kin_mev >= 0.0 && options.project_fs > 0.0 &&
           options.core_noise_limit > 0.0 &&
           options.bin_widths_mev.size() == 3;
}

std::string rank_path(const std::string& directory, int rank,
                      const char* prefix)
{
    std::ostringstream path;
    path << directory;
    if (!directory.empty() && directory[directory.size() - 1] != '/' &&
        directory[directory.size() - 1] != '\\') path << '/';
    path << prefix << rank << ".dat";
    return path.str();
}

SnapshotStats read_snapshot_stats(const std::string& directory, int ranks)
{
    SnapshotStats stats;
    if (ranks <= 0) return stats;
    stats.complete = true;
    for (int rank = 0; rank < ranks; ++rank) {
        std::ifstream input(rank_path(directory, rank,
                                      "tail_per_cell_stats_rank").c_str());
        if (!input) {
            stats.complete = false;
            continue;
        }
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            double x = 0.0, macro = 0.0, weight = 0.0, density = 0.0;
            if (row >> x >> macro >> weight >> density) {
                stats.macro_count += macro;
                stats.tail_number += weight;
            }
        }
    }
    return stats;
}

double tail_energy_mev(const BackgroundTailParticle& particle)
{
    const double gamma = std::sqrt(1.0 + particle.ux * particle.ux +
                                   particle.uy * particle.uy +
                                   particle.uz * particle.uz);
    return (gamma - 1.0) * Const::me * Const::c * Const::c /
           (1.0e6 * Const::eV);
}

void add_histogram(std::vector<double>& histogram, double energy_mev,
                   double weight, double width, double max_energy)
{
    if (!(weight >= 0.0) || !std::isfinite(weight) ||
        !std::isfinite(energy_mev) || energy_mev < 0.0) return;
    const int bin = static_cast<int>(std::floor(energy_mev / width));
    if (bin < 0 || static_cast<size_t>(bin) >= histogram.size()) return;
    histogram[static_cast<size_t>(bin)] += weight;
    (void)max_energy;
}

SpectrumStats build_spectrum(const Species& electrons,
                             const BackgroundTailPIC& tail,
                             const CylindricalVelocityGrid& cgrid,
                             double width, double max_energy)
{
    SpectrumStats result;
    const size_t bins = static_cast<size_t>(std::ceil(max_energy / width)) + 1;
    result.values.assign(bins, 0.0);
    result.bulk_values.assign(bins, 0.0);
    result.tail_values.assign(bins, 0.0);
    result.total = 0.0;
    const SpatialGrid& grid = *electrons.sgrid;
    const int nvmu = Param::Nvmu;
    std::vector<double> energy_edges(bins + 1, 0.0);
    for (size_t i = 0; i < energy_edges.size(); ++i)
        energy_edges[i] = static_cast<double>(i) * width * 1.0e6 * Const::eV;
    std::vector<std::vector<TailEnergyBinFraction> > fractions(
        static_cast<size_t>(nvmu));
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            fractions[static_cast<size_t>(iv) * Param::Nmu + imu] =
                TailSubcellQuadrature::energy_bin_fractions(
                    cgrid, iv, imu, energy_edges);
        }
    }
    for (int il = 0; il < grid.nx_local; ++il) {
        const int ix = grid.nghost + il;
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t slot = static_cast<size_t>(ix) * nvmu +
                                    static_cast<size_t>(iv) * Param::Nmu +
                                    static_cast<size_t>(imu);
                const double mass = electrons.f[slot];
                if (!std::isfinite(mass) || mass < 0.0) continue;
                const std::vector<TailEnergyBinFraction>& cell_fractions =
                    fractions[
                    static_cast<size_t>(iv) * Param::Nmu + imu];
                for (size_t q = 0; q < cell_fractions.size(); ++q) {
                    const size_t bin = static_cast<size_t>(
                        cell_fractions[q].bin);
                    if (bin >= bins) continue;
                    const double shell_mass =
                        mass * cell_fractions[q].mass_fraction;
                    result.values[bin] += shell_mass;
                    result.bulk_values[bin] += shell_mass;
                }
                result.total += mass;
            }
        }
    }
    for (size_t p = 0; p < tail.particles.size(); ++p) {
        const double energy = tail_energy_mev(tail.particles[p]);
        add_histogram(result.values, energy, tail.particles[p].weight,
                      width, max_energy);
        add_histogram(result.tail_values, energy, tail.particles[p].weight,
                      width, max_energy);
        result.total += tail.particles[p].weight;
    }
    return result;
}

SpectrumStats rebin_spectrum(const SpectrumStats& fine, double fine_width,
                             double coarse_width)
{
    SpectrumStats coarse;
    const size_t ratio = static_cast<size_t>(
        std::llround(coarse_width / fine_width));
    if (ratio == 0 ||
        std::fabs(static_cast<double>(ratio) * fine_width - coarse_width) >
            1.0e-12 * coarse_width) return coarse;
    const size_t bins = (fine.values.size() + ratio - 1) / ratio;
    coarse.values.assign(bins, 0.0);
    coarse.bulk_values.assign(bins, 0.0);
    coarse.tail_values.assign(bins, 0.0);
    coarse.total = fine.total;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const size_t bin = i / ratio;
        if (bin >= bins) continue;
        coarse.values[bin] += fine.values[i];
        coarse.bulk_values[bin] += fine.bulk_values[i];
        coarse.tail_values[bin] += fine.tail_values[i];
    }
    return coarse;
}

double reconstructed_l1(const SpectrumStats& fine,
                        const SpectrumStats& coarse, double fine_width,
                        double coarse_width, double window_lo,
                        double window_hi)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const double e0 = static_cast<double>(i) * fine_width;
        const double e1 = e0 + fine_width;
        if (e1 <= window_lo || e0 >= window_hi) continue;
        const int coarse_bin = static_cast<int>(std::floor(e0 / coarse_width));
        const double coarse_value =
            (coarse_bin >= 0 && static_cast<size_t>(coarse_bin) <
             coarse.values.size())
            ? coarse.values[static_cast<size_t>(coarse_bin)] *
              fine_width / coarse_width : 0.0;
        numerator += std::fabs(fine.values[i] - coarse_value);
        denominator += std::fabs(fine.values[i]);
    }
    return numerator / std::max(denominator, 1.0e-300);
}

double max_relative_bin_difference(const SpectrumStats& fine,
                                   const SpectrumStats& coarse,
                                   double fine_width, double coarse_width,
                                   double window_lo, double window_hi,
                                   double sparse_fraction,
                                   int& compared_bins,
                                   int& skipped_sparse_bins)
{
    double maximum = 0.0;
    double window_total = 0.0;
    compared_bins = 0;
    skipped_sparse_bins = 0;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const double e0 = static_cast<double>(i) * fine_width;
        const double e1 = e0 + fine_width;
        if (e1 > window_lo && e0 < window_hi)
            window_total += std::fabs(fine.values[i]);
    }
    const double sparse_floor = sparse_fraction * window_total;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const double e0 = static_cast<double>(i) * fine_width;
        const double e1 = e0 + fine_width;
        if (e1 <= window_lo || e0 >= window_hi) continue;
        const int coarse_bin = static_cast<int>(std::floor(e0 / coarse_width));
        const double coarse_value =
            (coarse_bin >= 0 && static_cast<size_t>(coarse_bin) <
             coarse.values.size())
            ? coarse.values[static_cast<size_t>(coarse_bin)] *
              fine_width / coarse_width : 0.0;
        if (std::fabs(fine.values[i]) <= sparse_floor &&
            std::fabs(coarse_value) <= sparse_floor) {
            ++skipped_sparse_bins;
            continue;
        }
        ++compared_bins;
        maximum = std::max(maximum,
            std::fabs(fine.values[i] - coarse_value) /
            std::max(std::max(std::fabs(fine.values[i]),
                              std::fabs(coarse_value)), sparse_floor));
    }
    return compared_bins > 0 ? maximum :
           std::numeric_limits<double>::quiet_NaN();
}

PerformanceStats read_performance_projection(
    const std::string& diagnostics_path, double project_fs,
    double macro_at_end, double projected_macro)
{
    PerformanceStats result;
    std::ifstream input(diagnostics_path.c_str());
    if (!input) return result;

    std::string header;
    if (!std::getline(input, header)) return result;
    std::istringstream names_in(header);
    std::vector<std::string> names;
    std::string name;
    while (names_in >> name) names.push_back(name);
    int step_col = -1, time_col = -1, accepted_col = -1;
    int wall_col = -1, macro_col = -1;
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == "step") step_col = static_cast<int>(i);
        else if (names[i] == "time_s") time_col = static_cast<int>(i);
        else if (names[i] == "accepted") accepted_col = static_cast<int>(i);
        else if (names[i] == "wall_s") wall_col = static_cast<int>(i);
        else if (names[i] == "tail_particle_count")
            macro_col = static_cast<int>(i);
    }
    if (step_col < 0 || time_col < 0 || accepted_col < 0 ||
        wall_col < 0 || macro_col < 0) return result;

    std::vector<double> times, walls, macros;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        std::vector<double> values;
        double value = 0.0;
        while (row >> value) values.push_back(value);
        const int needed = std::max(
            std::max(step_col, time_col),
            std::max(std::max(accepted_col, wall_col), macro_col));
        if (static_cast<int>(values.size()) <= needed ||
            values[static_cast<size_t>(accepted_col)] != 1.0 ||
            !std::isfinite(values[static_cast<size_t>(wall_col)]) ||
            !(values[static_cast<size_t>(wall_col)] >= 0.0)) continue;
        times.push_back(values[static_cast<size_t>(time_col)] / Const::femto);
        walls.push_back(values[static_cast<size_t>(wall_col)]);
        macros.push_back(values[static_cast<size_t>(macro_col)]);
    }
    if (times.size() < 3) return result;

    const size_t first = times.size() > 200 ? times.size() - 200 : 0;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    const double count = static_cast<double>(times.size() - first);
    for (size_t i = first; i < times.size(); ++i) {
        sx += macros[i];
        sy += walls[i];
        sxx += macros[i] * macros[i];
        sxy += macros[i] * walls[i];
    }
    const double denominator = count * sxx - sx * sx;
    double slope = denominator > 0.0 ?
        (count * sxy - sx * sy) / denominator : 0.0;
    slope = std::max(0.0, slope);
    double intercept = (sy - slope * sx) / count;
    intercept = std::max(0.0, intercept);

    std::vector<double> time_steps;
    for (size_t i = 1; i < times.size(); ++i) {
        const double dt = times[i] - times[i - 1];
        if (dt > 0.0 && std::isfinite(dt)) time_steps.push_back(dt);
    }
    if (time_steps.empty()) return result;
    std::sort(time_steps.begin(), time_steps.end());
    const double dt_fs = time_steps[time_steps.size() / 2];
    const double remaining_fs = std::max(0.0, project_fs - times.back());
    const double steps = remaining_fs / dt_fs;
    const double average_macro = 0.5 * (macro_at_end + projected_macro);
    const double projected_per_step = intercept + slope * average_macro;

    result.complete = std::isfinite(projected_per_step) &&
                      projected_per_step >= 0.0 && std::isfinite(steps);
    result.projected_wall_seconds_per_step = projected_per_step;
    result.projected_wall_seconds = projected_per_step * steps;
    result.dt_fs = dt_fs;
    result.fit_intercept_seconds = intercept;
    result.fit_slope_seconds_per_particle = slope;
    return result;
}

CoreNoiseStats core_high_frequency_noise(
    const std::vector<double>& tail_density,
    const std::vector<double>& background_density)
{
    CoreNoiseStats result;
    if (tail_density.size() < 32 ||
        tail_density.size() != background_density.size()) return result;
    const size_t begin = tail_density.size() / 10;
    const size_t end = tail_density.size() - begin;
    const int radius = 10;
    double residual_sq = 0.0, smooth_sq = 0.0;
    double tail_sq = 0.0, background_sq = 0.0;
    for (size_t i = begin; i < end; ++i) {
        const size_t lo = i > static_cast<size_t>(radius) ?
                          i - static_cast<size_t>(radius) : 0;
        const size_t hi = std::min(tail_density.size() - 1,
                                   i + static_cast<size_t>(radius));
        double smooth = 0.0;
        for (size_t j = lo; j <= hi; ++j) smooth += tail_density[j];
        smooth /= static_cast<double>(hi - lo + 1);
        const double residual = tail_density[i] - smooth;
        residual_sq += residual * residual;
        smooth_sq += smooth * smooth;
        tail_sq += tail_density[i] * tail_density[i];
        background_sq += background_density[i] * background_density[i];
    }
    const double background_rms = std::sqrt(background_sq /
        static_cast<double>(end - begin));
    double max_ratio = 0.0;
    const double background_floor = 1.0e-12 * background_rms;
    for (size_t i = begin; i < end; ++i) {
        if (std::fabs(background_density[i]) > background_floor) {
            max_ratio = std::max(max_ratio,
                std::fabs(tail_density[i] / background_density[i]));
        }
    }
    result.complete = std::isfinite(background_sq) && background_sq > 0.0;
    result.relative_to_tail_envelope = std::sqrt(residual_sq /
        std::max(smooth_sq, std::numeric_limits<double>::min()));
    result.relative_to_background_rms = std::sqrt(residual_sq /
        std::max(background_sq, std::numeric_limits<double>::min()));
    result.tail_rms_relative_to_background_rms = std::sqrt(tail_sq /
        std::max(background_sq, std::numeric_limits<double>::min()));
    result.max_tail_relative_to_background = max_ratio;
    return result;
}

SpectrumGapStats inspect_artificial_spectrum_gaps(
    const SpectrumStats& fine, double width, double window_lo,
    double window_hi)
{
    SpectrumGapStats result;
    double window_total = 0.0;
    double window_peak = 0.0;
    size_t first = fine.values.size();
    size_t last = 0;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const double center = (static_cast<double>(i) + 0.5) * width;
        if (center >= window_lo && center <= window_hi) {
            window_total += std::fabs(fine.values[i]);
            window_peak = std::max(window_peak, std::fabs(fine.values[i]));
            first = std::min(first, i);
            last = std::max(last, i);
        }
    }
    if (!(window_total > 0.0) || first >= fine.values.size() || first >= last)
        return result;
    const double sparse_floor = 1.0e-12 * window_total;
    const double valley_ceiling = std::max(sparse_floor, 0.01 * window_peak);
    for (size_t i = first; i < last; ++i) {
        const double a = std::fabs(fine.values[i]);
        const double b = std::fabs(fine.values[i + 1]);
        if (a > sparse_floor && b > sparse_floor) {
            result.max_adjacent_jump_ratio = std::max(
                result.max_adjacent_jump_ratio,
                std::max(a, b) / std::max(std::min(a, b), sparse_floor));
        }
    }
    size_t i = first;
    while (i <= last) {
        if (std::fabs(fine.values[i]) > valley_ceiling) {
            ++i;
            continue;
        }
        const size_t run_begin = i;
        while (i <= last &&
               std::fabs(fine.values[i]) <= valley_ceiling) ++i;
        const size_t run_end = i - 1;
        const bool bounded_left = run_begin > first &&
            std::fabs(fine.values[run_begin - 1]) > valley_ceiling;
        const bool bounded_right = run_end < last &&
            std::fabs(fine.values[run_end + 1]) > valley_ceiling;
        if (bounded_left && bounded_right) ++result.deep_valley_count;
    }
    result.pass = result.deep_valley_count == 0 &&
        result.max_adjacent_jump_ratio <= 100.0;
    return result;
}

void write_spectrum_detail(const Options& options, const SpectrumStats& fine,
                           const std::vector<SpectrumStats>& spectra,
                           double fine_width, int rank)
{
    if (rank != 0) return;
    const std::string path = options.result + ".spectrum.dat";
    std::ofstream out(path.c_str());
    if (!out) {
        std::cerr << "cannot write spectrum detail: " << path << "\n";
        return;
    }
    out << std::setprecision(17);
    double represented = 0.0, represented_bulk = 0.0, represented_tail = 0.0;
    for (size_t i = 0; i < fine.values.size(); ++i) {
        represented += fine.values[i];
        represented_bulk += fine.bulk_values[i];
        represented_tail += fine.tail_values[i];
    }
    out << "# total_input_number=" << fine.total
        << " represented_number=" << represented
        << " relative_number_residual="
        << std::fabs(represented - fine.total) /
           std::max(std::fabs(fine.total), 1.0) << "\n";
    out << "# represented_bulk_number=" << represented_bulk
        << " represented_tail_number=" << represented_tail << "\n";
    out << "energy_mev fine_combined_0p05 fine_bulk_0p05 fine_tail_0p05 "
           "rebin_0p1 rebin_0p2 rebin_0p4\n";
    for (size_t i = 0; i < fine.values.size(); ++i) {
        const double e0 = static_cast<double>(i) * fine_width;
        out << e0 + 0.5 * fine_width << " " << fine.values[i] << " "
            << fine.bulk_values[i] << " " << fine.tail_values[i];
        for (size_t s = 0; s < spectra.size(); ++s) {
            const double coarse_width = options.bin_widths_mev[s];
            const size_t coarse_bin = static_cast<size_t>(
                std::floor(e0 / coarse_width));
            const double reconstructed = coarse_bin < spectra[s].values.size()
                ? spectra[s].values[coarse_bin] * fine_width / coarse_width
                : 0.0;
            out << " " << reconstructed;
        }
        out << "\n";
    }
}

void write_result(const Options& options, const SnapshotStats& snapshot25,
                  const SnapshotStats& snapshot40, double tail_number,
                  double thermalized_number, double total_number,
                  const CoreNoiseStats& core_noise,
                  const std::vector<double>& l1,
                  double max_bin_relative, int compared_bins,
                  int skipped_sparse_bins,
                  const SpectrumGapStats& spectrum_gap,
                  double global_macro,
                  double max_local_macro, double projected_memory_global,
                  double projected_memory_max_rank,
                  const PerformanceStats& performance, bool checkpoint_ok,
                  bool residence_available, int rank)
{
    if (rank != 0) return;
    std::ofstream out(options.result.c_str());
    if (!out) {
        std::cerr << "cannot write result: " << options.result << "\n";
        return;
    }
    out << std::setprecision(17);
    const bool finite_metrics = checkpoint_ok && std::isfinite(tail_number) &&
        std::isfinite(total_number) && core_noise.complete;
    const bool rebin_shape_pass = finite_metrics && l1.size() == 3 &&
        std::isfinite(l1[0]) && std::isfinite(l1[1]) &&
        std::isfinite(l1[2]) &&
        l1[0] <= 0.05 && l1[1] <= 0.05 && l1[2] <= 0.05 &&
        compared_bins > 0 &&
        std::isfinite(max_bin_relative) && max_bin_relative <= 0.10;
    // Rebinning one checkpoint cannot establish simulation convergence.
    // It is retained only as a shape-sensitivity diagnostic.
    const int spectrum_gate_pass = -1;
    const bool artificial_gap_pass = spectrum_gap.pass;
    const double physical_growth_rate =
        (snapshot25.complete && snapshot40.complete && options.project_fs > 40.0)
        ? (snapshot40.tail_number - snapshot25.tail_number) / 15.0 :
          std::numeric_limits<double>::quiet_NaN();
    const double projected_physical_number = std::isfinite(physical_growth_rate)
        ? std::max(0.0, snapshot40.tail_number +
                   physical_growth_rate * (options.project_fs - 40.0))
        : std::numeric_limits<double>::quiet_NaN();
    const double macro_growth_rate =
        (snapshot25.complete && snapshot40.complete && options.project_fs > 40.0)
        ? (snapshot40.macro_count - snapshot25.macro_count) / 15.0 :
          std::numeric_limits<double>::quiet_NaN();
    const double projected_macro = std::isfinite(macro_growth_rate)
        ? std::max(0.0, snapshot40.macro_count +
                   macro_growth_rate * (options.project_fs - 40.0))
        : std::numeric_limits<double>::quiet_NaN();
    const double projected_local_macro =
        global_macro > 0.0 && std::isfinite(projected_macro)
        ? projected_macro * max_local_macro / global_macro :
          std::numeric_limits<double>::quiet_NaN();
    const double tail_fraction = total_number > 0.0
        ? tail_number / total_number : 0.0;
    const double thermalized_fraction = tail_number > 0.0
        ? thermalized_number / tail_number : 0.0;
    const double global_particle_budget = 1.2e7;
    const double local_particle_budget = 2.0e6;
    const bool particle_projection_available =
        std::isfinite(projected_macro) &&
        std::isfinite(projected_local_macro);
    const bool particle_budget_pass = particle_projection_available &&
        projected_macro <= global_particle_budget &&
        projected_local_macro <= local_particle_budget;
    const int population_controller_required =
        particle_projection_available ? (particle_budget_pass ? 0 : 1) : -1;
    const bool core_noise_pass = core_noise.complete &&
        core_noise.relative_to_background_rms <= options.core_noise_limit;
    const int h10_required = residence_available ?
        ((tail_fraction >= 0.01 || thermalized_fraction >= 0.10 ||
          !core_noise_pass) ? 1 : 0) : -1;
    // Independent conversion-bin runs and residence history are both absent
    // from this checkpoint-only invocation, so the combined decision remains
    // incomplete.  The independently evaluable sub-gates are still reported.
    const char* status = !checkpoint_ok ? "FAIL" :
        ((!residence_available || spectrum_gate_pass < 0) ? "INCOMPLETE" :
         ((!artificial_gap_pass || population_controller_required == 1 ||
           h10_required != 0) ? "FAIL" : "PASS"));
    out << "audit_schema=hybrid_checkpoint_gate_audit_v5\n";
    out << "status=" << status << "\n";
    out << "spectrum_convergence_available=0\n";
    out << "spectrum_gate_pass=" << spectrum_gate_pass << "\n";
    out << "spectrum_rebin_metrics_informational=1\n";
    out << "spectrum_rebin_shape_pass=" << (rebin_shape_pass ? 1 : 0)
        << "\n";
    out << "bin_width_0p1_pass=-1\n";
    out << "bin_width_0p2_pass=-1\n";
    out << "bin_width_0p4_pass=-1\n";
    out << "threshold_l1_0p1=" << (l1.size() > 0 ? l1[0] : std::numeric_limits<double>::quiet_NaN()) << "\n";
    out << "threshold_l1_0p2=" << (l1.size() > 1 ? l1[1] : std::numeric_limits<double>::quiet_NaN()) << "\n";
    out << "threshold_l1_0p4=" << (l1.size() > 2 ? l1[2] : std::numeric_limits<double>::quiet_NaN()) << "\n";
    out << "threshold_max_bin_relative=" << max_bin_relative << "\n";
    out << "threshold_compared_non_sparse_bin_count=" << compared_bins << "\n";
    out << "threshold_skipped_sparse_bin_count=" << skipped_sparse_bins << "\n";
    out << "threshold_reference_bin_width_mev="
        << 0.5 * *std::min_element(options.bin_widths_mev.begin(),
                                   options.bin_widths_mev.end()) << "\n";
    out << "threshold_artificial_gap_count="
        << spectrum_gap.deep_valley_count << "\n";
    out << "threshold_deep_valley_count="
        << spectrum_gap.deep_valley_count << "\n";
    out << "threshold_max_adjacent_jump_ratio="
        << spectrum_gap.max_adjacent_jump_ratio << "\n";
    out << "threshold_artificial_gap_pass="
        << (artificial_gap_pass ? 1 : 0) << "\n";
    out << "threshold_spectrum_file=" << options.result
        << ".spectrum.dat\n";
    out << "tail_physical_number_fraction=" << tail_fraction << "\n";
    out << "thermalized_tail_fraction=" << thermalized_fraction << "\n";
    out << "thermalized_tail_residence_available=" << (residence_available ? 1 : 0) << "\n";
    out << "core_tail_density_noise_relative="
        << core_noise.relative_to_background_rms << "\n";
    out << "core_tail_density_noise_relative_to_tail_envelope="
        << core_noise.relative_to_tail_envelope << "\n";
    out << "core_tail_density_rms_relative_to_background="
        << core_noise.tail_rms_relative_to_background_rms << "\n";
    out << "core_tail_density_max_relative_to_background="
        << core_noise.max_tail_relative_to_background << "\n";
    out << "core_tail_density_noise_limit=" << options.core_noise_limit
        << "\n";
    out << "core_tail_density_noise_pass=" << (core_noise_pass ? 1 : 0)
        << "\n";
    out << "tail_physical_number_projected_120fs="
        << projected_physical_number << "\n";
    out << "tail_particles_projected_120fs=" << projected_macro << "\n";
    out << "tail_particles_local_projected_120fs="
        << projected_local_macro << "\n";
    out << "tail_particles_checkpoint_global=" << global_macro << "\n";
    out << "tail_particles_checkpoint_local_max=" << max_local_macro << "\n";
    out << "tail_particles_budget_global=" << global_particle_budget << "\n";
    out << "tail_particles_budget_local=" << local_particle_budget << "\n";
    out << "memory_projected_global_bytes=" << projected_memory_global << "\n";
    out << "memory_projected_max_rank_bytes=" << projected_memory_max_rank << "\n";
    out << "wall_seconds_projected_40_to_120fs="
        << performance.projected_wall_seconds << "\n";
    out << "wall_seconds_per_step_projected="
        << performance.projected_wall_seconds_per_step << "\n";
    out << "performance_dt_fs=" << performance.dt_fs << "\n";
    out << "performance_fit_intercept_seconds="
        << performance.fit_intercept_seconds << "\n";
    out << "performance_fit_slope_seconds_per_particle="
        << performance.fit_slope_seconds_per_particle << "\n";
    out << "particle_projection_available="
        << (particle_projection_available ? 1 : 0) << "\n";
    out << "particle_budget_pass=" << (particle_budget_pass ? 1 : 0)
        << "\n";
    out << "performance_projection_available="
        << (performance.complete ? 1 : 0) << "\n";
    // The audit has no user-declared wall-time limit.  Therefore the
    // resource gate is the independently testable macro-particle budget;
    // wall time remains a reported projection, not a fabricated pass/fail.
    out << "resource_projection_available="
        << (particle_projection_available ? 1 : 0) << "\n";
    out << "resource_gate_pass=" << (particle_budget_pass ? 1 : 0) << "\n";
    out << "h10_required=" << h10_required << "\n";
    out << "population_controller_required="
        << population_controller_required << "\n";
    out << "checkpoint_audit_finite=" << (finite_metrics ? 1 : 0) << "\n";
    out << "checkpoint_tail_number=" << tail_number << "\n";
    out << "checkpoint_thermalized_number=" << thermalized_number << "\n";
    out << "snapshot25_complete=" << (snapshot25.complete ? 1 : 0) << "\n";
    out << "snapshot40_complete=" << (snapshot40.complete ? 1 : 0) << "\n";
    out << "snapshot25_tail_number=" << snapshot25.tail_number << "\n";
    out << "snapshot40_tail_number=" << snapshot40.tail_number << "\n";
    out << "residence_status="
        << (residence_available ? "available" :
            "unavailable_in_checkpoint_schema") << "\n";
    out << "decision_note="
        << ((!residence_available || spectrum_gate_pass < 0) ?
            "independent_spectrum_runs_and_residence_required" :
            "complete") << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Options options;
    const bool parsed = parse_options(argc, argv, options);
    int parsed_all = parsed ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &parsed_all, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!parsed_all) {
        if (rank == 0) std::cerr << "invalid hybrid checkpoint audit arguments\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    try {
        grid.init_with_domain(rank, size, Param::nx, Param::Lx);
    } catch (const std::exception& error) {
        if (rank == 0) std::cerr << "grid initialization failed: "
                                  << error.what() << "\n";
        MPI_Finalize();
        return 2;
    }
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    VpfpCheckpointControl control;
    VpfpCheckpointTailState tail_state;
    std::string error;
    bool local_ok = read_vpfp_checkpoint(
        options.checkpoint, control, electrons, beam, fields, grid,
        &tail_state, rank, size, error);
    int global_ok = local_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!global_ok || !local_ok || !tail_state.present) {
        if (rank == 0) {
            std::cerr << "checkpoint audit read failed";
            if (!local_ok) std::cerr << ": " << error;
            std::cerr << "\n";
        }
        SnapshotStats empty;
        write_result(options, empty, empty, 0.0, 0.0, 0.0,
                     CoreNoiseStats(),
                     std::vector<double>(3, std::numeric_limits<double>::quiet_NaN()),
                     std::numeric_limits<double>::quiet_NaN(), 0, 0,
                     SpectrumGapStats(), 0.0, 0.0,
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::quiet_NaN(),
                     PerformanceStats(), false, false, rank);
        MPI_Finalize();
        return 2;
    }

    double local_max_tail_energy = 0.0;
    for (size_t p = 0; p < tail_state.tail.particles.size(); ++p)
        local_max_tail_energy = std::max(
            local_max_tail_energy, tail_energy_mev(tail_state.tail.particles[p]));
    double max_tail_energy = 0.0;
    MPI_Allreduce(&local_max_tail_energy, &max_tail_energy, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    const double max_upar = std::max(
        std::fabs(electrons.cgrid.upar_faces.front()),
        std::fabs(electrons.cgrid.upar_faces.back()));
    const double max_uperp = std::max(
        std::fabs(electrons.cgrid.uperp_faces.front()),
        std::fabs(electrons.cgrid.uperp_faces.back()));
    const double bulk_grid_max_energy =
        (std::sqrt(1.0 + max_upar * max_upar + max_uperp * max_uperp) - 1.0) *
        Const::me * Const::c * Const::c / (1.0e6 * Const::eV);
    const double max_energy_mev = std::max(
        options.kout_mev + 2.0,
        1.001 * std::max(max_tail_energy, bulk_grid_max_energy));
    // Use an independently finer reference.  Reusing the 0.1 MeV candidate
    // as its own reference would make that convergence row pass trivially.
    const double fine_width = 0.5 *
        *std::min_element(options.bin_widths_mev.begin(),
                          options.bin_widths_mev.end());
    SpectrumStats fine = build_spectrum(electrons, tail_state.tail,
                                         electrons.cgrid, fine_width,
                                         max_energy_mev);
    std::vector<SpectrumStats> spectra;
    for (size_t i = 0; i < options.bin_widths_mev.size(); ++i)
        spectra.push_back(rebin_spectrum(
            fine, fine_width, options.bin_widths_mev[i]));
    // Both operands of every comparison must represent the same global
    // physical state.  The original audit reduced only the coarse spectra,
    // which compared one rank's fine spectrum with the global coarse one.
    globalize_spectrum(fine);
    for (size_t s = 0; s < spectra.size(); ++s)
        globalize_spectrum(spectra[s]);
    write_spectrum_detail(options, fine, spectra, fine_width, rank);
    std::vector<double> l1(3, 0.0);
    for (size_t i = 0; i < spectra.size(); ++i)
        l1[i] = reconstructed_l1(
            fine, spectra[i], fine_width, options.bin_widths_mev[i],
            std::max(0.0, options.kout_mev - 0.4), options.kout_mev + 0.4);
    double max_bin_relative = 0.0;
    int compared_bins = 0;
    int skipped_sparse_bins = 0;
    for (size_t i = 0; i < spectra.size(); ++i) {
        int compared = 0, skipped = 0;
        const double candidate = max_relative_bin_difference(
            fine, spectra[i], fine_width, options.bin_widths_mev[i],
            std::max(0.0, options.kout_mev - 0.4),
            options.kout_mev + 0.4, 1.0e-12, compared, skipped);
        if (std::isfinite(candidate))
            max_bin_relative = std::max(max_bin_relative, candidate);
        compared_bins += compared;
        skipped_sparse_bins += skipped;
    }
    const SpectrumGapStats spectrum_gap = inspect_artificial_spectrum_gaps(
        fine, fine_width, std::max(0.0, options.kout_mev - 0.4),
        options.kout_mev + 0.4);

    double local_tail_number = 0.0;
    double local_thermalized = 0.0;
    const double local_bulk_number = electrons.total_particle_number();
    std::vector<double> global_tail_density(
        static_cast<size_t>(grid.nx_global), 0.0);
    std::vector<double> global_background_density(
        static_cast<size_t>(grid.nx_global), 0.0);
    for (size_t p = 0; p < tail_state.tail.particles.size(); ++p) {
        local_tail_number += tail_state.tail.particles[p].weight;
        if (tail_energy_mev(tail_state.tail.particles[p]) <= options.kin_mev)
            local_thermalized += tail_state.tail.particles[p].weight;
    }
    for (int il = 0; il < grid.nx_local; ++il) {
        const int global_ix = grid.ix_start + il;
        if (global_ix >= 0 && global_ix < grid.nx_global &&
            static_cast<size_t>(il) < tail_state.tail.density.size()) {
            global_tail_density[static_cast<size_t>(global_ix)] =
                tail_state.tail.density[static_cast<size_t>(il)];
            global_background_density[static_cast<size_t>(global_ix)] =
                electrons.number_density[static_cast<size_t>(il)];
        }
    }
    double sums[3] = {local_tail_number, local_thermalized,
                      local_bulk_number};
    MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, global_tail_density.data(), grid.nx_global,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, global_background_density.data(),
                  grid.nx_global, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // Measure grid-scale structure relative to a 21-cell local envelope in
    // the central 80% of the domain.  The previous all-domain coefficient of
    // variation treated the physical Tail envelope as numerical noise.
    const CoreNoiseStats core_noise = core_high_frequency_noise(
        global_tail_density, global_background_density);
    const double local_base_bytes =
        static_cast<double>(electrons.f.size() * sizeof(double) +
                            fields.Ex_face.size() * sizeof(double) +
                            fields.Ex.size() * sizeof(double) +
                            fields.phi.size() * sizeof(double));
    double base_global_bytes = 0.0, base_max_bytes = 0.0;
    MPI_Allreduce(&local_base_bytes, &base_global_bytes, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_base_bytes, &base_max_bytes, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    double local_macro = static_cast<double>(tail_state.tail.particles.size());
    double global_macro = 0.0, max_macro = 0.0;
    MPI_Allreduce(&local_macro, &global_macro, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_macro, &max_macro, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    const SnapshotStats snapshot25 = rank == 0
        ? read_snapshot_stats(options.snapshot25, size) : SnapshotStats();
    const SnapshotStats snapshot40 = rank == 0
        ? read_snapshot_stats(options.snapshot40, size) : SnapshotStats();
    int snapshot_flags[2] = {snapshot25.complete ? 1 : 0,
                             snapshot40.complete ? 1 : 0};
    MPI_Bcast(snapshot_flags, 2, MPI_INT, 0, MPI_COMM_WORLD);
    SnapshotStats snapshot25_b = snapshot25;
    SnapshotStats snapshot40_b = snapshot40;
    double snapshot_values[4] = {snapshot25.tail_number,
                                 snapshot25.macro_count,
                                 snapshot40.tail_number,
                                 snapshot40.macro_count};
    MPI_Bcast(snapshot_values, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    snapshot25_b.complete = snapshot_flags[0] != 0;
    snapshot40_b.complete = snapshot_flags[1] != 0;
    snapshot25_b.tail_number = snapshot_values[0];
    snapshot25_b.macro_count = snapshot_values[1];
    snapshot40_b.tail_number = snapshot_values[2];
    snapshot40_b.macro_count = snapshot_values[3];
    const double macro_growth =
        (snapshot25_b.complete && snapshot40_b.complete)
        ? (snapshot40_b.macro_count - snapshot25_b.macro_count) / 15.0
        : std::numeric_limits<double>::quiet_NaN();
    const double projected_macro = std::isfinite(macro_growth)
        ? std::max(0.0, snapshot40_b.macro_count +
                   macro_growth * (options.project_fs - 40.0))
        : global_macro;
    const double projected_local_macro = global_macro > 0.0
        ? projected_macro * max_macro / global_macro : max_macro;
    const double projected_memory_global =
        base_global_bytes + projected_macro * sizeof(BackgroundTailParticle);
    const double projected_memory_max_rank =
        base_max_bytes + projected_local_macro * sizeof(BackgroundTailParticle);
    const std::string diagnostics_path =
        parent_directory(options.snapshot40) + "/vpfp_step_diagnostics.dat";
    const PerformanceStats performance = read_performance_projection(
        diagnostics_path, options.project_fs, snapshot40_b.macro_count,
        projected_macro);
    write_result(options, snapshot25_b, snapshot40_b, sums[0], sums[1],
                 sums[0] + sums[2], core_noise, l1, max_bin_relative,
                 compared_bins, skipped_sparse_bins, spectrum_gap,
                 global_macro, max_macro,
                 projected_memory_global, projected_memory_max_rank,
                 performance, true, false, rank);
    MPI_Finalize();
    return 0;
}
