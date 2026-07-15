#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Sample {
    double time_fs;
    double e_re;
    double e_im;
    double e_abs;
    double total_energy;
    double energy_drift;
    double r_fv;
    double r_couple;
    double j_pair_linf;
    double limiter_active;
    double limiter_active_core;
    double limiter_active_boundary;
    double min_alpha;
    int strict_accept;
};

struct Series {
    std::vector<Sample> samples;
};

bool read_series(const std::string& path, Series& series)
{
    std::ifstream input(path.c_str());
    if (!input) return false;
    std::map<std::string, size_t> columns;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::istringstream header(line.substr(1));
            std::string name;
            size_t index = 0;
            while (header >> name) columns[name] = index++;
            continue;
        }
        std::istringstream row(line);
        std::vector<double> values;
        double value = 0.0;
        while (row >> value) values.push_back(value);
        const char* required[] = {
            "time_fs", "E_re", "E_im", "E_abs", "strict_accept",
            "total_energy", "energy_drift", "R_FV", "R_couple",
            "JN_minus_GstarJE_linf", "limiter_active",
            "limiter_active_core", "limiter_active_boundary", "min_alpha"};
        bool present = true;
        for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i)
            present = present && columns.find(required[i]) != columns.end();
        if (!present) return false;
        size_t required_last = 0;
        for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i)
            required_last = std::max(required_last, columns[required[i]]);
        if (values.size() <= required_last) return false;
        Sample sample;
        sample.time_fs = values[columns["time_fs"]];
        sample.e_re = values[columns["E_re"]];
        sample.e_im = values[columns["E_im"]];
        sample.e_abs = values[columns["E_abs"]];
        sample.strict_accept = static_cast<int>(values[columns["strict_accept"]]);
        sample.total_energy = values[columns["total_energy"]];
        sample.energy_drift = values[columns["energy_drift"]];
        sample.r_fv = values[columns["R_FV"]];
        sample.r_couple = values[columns["R_couple"]];
        sample.j_pair_linf = values[columns["JN_minus_GstarJE_linf"]];
        sample.limiter_active = values[columns["limiter_active"]];
        sample.limiter_active_core = values[columns["limiter_active_core"]];
        sample.limiter_active_boundary = values[columns["limiter_active_boundary"]];
        sample.min_alpha = values[columns["min_alpha"]];
        if (!std::isfinite(sample.time_fs) || !std::isfinite(sample.e_re) ||
            !std::isfinite(sample.e_im) || !std::isfinite(sample.e_abs))
            return false;
        series.samples.push_back(sample);
    }
    return series.samples.size() >= 3;
}

double unwrap_near(double previous, double value)
{
    const double two_pi = 2.0 * Const::pi;
    while (value - previous > Const::pi) value -= two_pi;
    while (value - previous < -Const::pi) value += two_pi;
    return value;
}

double fitted_angular_frequency(const Series& series)
{
    if (series.samples.size() < 3) return std::numeric_limits<double>::quiet_NaN();
    double sum_t = 0.0;
    double sum_p = 0.0;
    double sum_tt = 0.0;
    double sum_tp = 0.0;
    double phase_previous = std::atan2(series.samples[0].e_im,
                                       series.samples[0].e_re);
    for (size_t i = 0; i < series.samples.size(); ++i) {
        double phase = std::atan2(series.samples[i].e_im, series.samples[i].e_re);
        if (i > 0) phase = unwrap_near(phase_previous, phase);
        phase_previous = phase;
        const double time = series.samples[i].time_fs * Const::femto;
        sum_t += time;
        sum_p += phase;
        sum_tt += time * time;
        sum_tp += time * phase;
    }
    const double n = static_cast<double>(series.samples.size());
    const double denominator = n * sum_tt - sum_t * sum_t;
    return denominator != 0.0
        ? (n * sum_tp - sum_t * sum_p) / denominator
        : std::numeric_limits<double>::quiet_NaN();
}

double max_abs_member(const Series& series, int member)
{
    double value = 0.0;
    for (size_t i = 0; i < series.samples.size(); ++i) {
        const Sample& sample = series.samples[i];
        double selected = 0.0;
        if (member == 0) selected = sample.energy_drift;
        else if (member == 1) selected = sample.r_fv;
        else if (member == 2) selected = sample.r_couple;
        else if (member == 3) selected = sample.j_pair_linf;
        else if (member == 4) selected = sample.limiter_active;
        else if (member == 5) selected = sample.limiter_active_core;
        else if (member == 6) selected = sample.limiter_active_boundary;
        value = std::max(value, std::fabs(selected));
    }
    return value;
}

double min_alpha(const Series& series)
{
    double value = 1.0;
    for (size_t i = 0; i < series.samples.size(); ++i)
        value = std::min(value, series.samples[i].min_alpha);
    return value;
}

bool write_audit(const std::string& output_path, const Series& series,
                 const Series* reference, int mode)
{
    std::ofstream output(output_path.c_str());
    if (!output) return false;
    const double omega_fit = fitted_angular_frequency(series);
    const double k = 2.0 * Const::pi * static_cast<double>(mode) / Param::Lx;
    const double vth = std::sqrt(Param::temperature_e / Const::me);
    const double omega_linear = std::sqrt(Param::omega_pe * Param::omega_pe +
                                          3.0 * k * k * vth * vth);
    const double omega_relative_error = std::fabs(std::fabs(omega_fit) - omega_linear) /
        std::max(std::numeric_limits<double>::min(), omega_linear);
    int strict_count = 0;
    for (size_t i = 0; i < series.samples.size(); ++i)
        strict_count += series.samples[i].strict_accept != 0 ? 1 : 0;

    output << std::scientific << std::setprecision(17)
           << "format_version=1\n"
           << "sample_count=" << series.samples.size() << "\n"
           << "strict_accept_count=" << strict_count << "\n"
           << "strict_accept_fraction=" << static_cast<double>(strict_count) /
              series.samples.size() << "\n"
           << "time_start_fs=" << series.samples.front().time_fs << "\n"
           << "time_end_fs=" << series.samples.back().time_fs << "\n"
           << "mode=" << mode << "\n"
           << "omega_fit=" << omega_fit << "\n"
           << "omega_linear_warm=" << omega_linear << "\n"
           << "omega_relative_error=" << omega_relative_error << "\n"
           << "max_abs_energy_drift=" << max_abs_member(series, 0) << "\n"
           << "max_abs_R_FV=" << max_abs_member(series, 1) << "\n"
           << "max_abs_R_couple=" << max_abs_member(series, 2) << "\n"
           << "max_JN_minus_GstarJE_linf=" << max_abs_member(series, 3) << "\n"
           << "max_limiter_active_fraction=" << max_abs_member(series, 4) << "\n"
           << "max_limiter_active_fraction_core=" << max_abs_member(series, 5) << "\n"
           << "max_limiter_active_fraction_boundary=" << max_abs_member(series, 6) << "\n"
           << "min_alpha=" << min_alpha(series) << "\n";
    if (reference) {
        const size_t count = std::min(series.samples.size(), reference->samples.size());
        double phase_difference_sum = 0.0;
        double amplitude_ratio_sum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            const Sample& lhs = series.samples[i];
            const Sample& rhs = reference->samples[i];
            const double lhs_phase = std::atan2(lhs.e_im, lhs.e_re);
            const double rhs_phase = std::atan2(rhs.e_im, rhs.e_re);
            phase_difference_sum += std::atan2(std::sin(lhs_phase - rhs_phase),
                                               std::cos(lhs_phase - rhs_phase));
            amplitude_ratio_sum += lhs.e_abs / std::max(
                std::numeric_limits<double>::min(), rhs.e_abs);
        }
        output << "reference_sample_count=" << count << "\n"
               << "mean_phase_difference=" << phase_difference_sum /
                  std::max<size_t>(1, count) << "\n"
               << "mean_amplitude_ratio=" << amplitude_ratio_sum /
                  std::max<size_t>(1, count) << "\n";
    }
    return static_cast<bool>(output);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: production_langmuir_wave_audit series.dat "
                  << "[--reference series.dat] [--output audit.result] [--mode M]\n";
        return 2;
    }
    const std::string input_path(argv[1]);
    std::string reference_path;
    std::string output_path = input_path + ".audit";
    int mode = 8;
    for (int i = 2; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--reference" && i + 1 < argc)
            reference_path = argv[++i];
        else if (arg == "--output" && i + 1 < argc)
            output_path = argv[++i];
        else if (arg == "--mode" && i + 1 < argc) {
            char* end = 0;
            const long parsed = std::strtol(argv[++i], &end, 10);
            if (end == argv[i] || *end != '\0' || parsed < 1 ||
                parsed > std::numeric_limits<int>::max())
                return 2;
            mode = static_cast<int>(parsed);
        }
        else return 2;
    }
    Series series;
    Series reference;
    if (!read_series(input_path, series) ||
        (!reference_path.empty() && !read_series(reference_path, reference))) {
        std::cerr << "failed to read Langmuir series\n";
        return 1;
    }
    if (!write_audit(output_path, series,
                     reference_path.empty() ? 0 : &reference, mode)) {
        std::cerr << "failed to write audit\n";
        return 1;
    }
    return 0;
}
