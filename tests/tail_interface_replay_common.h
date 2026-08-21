#ifndef TAIL_INTERFACE_REPLAY_COMMON_H
#define TAIL_INTERFACE_REPLAY_COMMON_H

// Shared reader / result writer for the 16B replay tests
// (tail_interface_grid_replay_test / tail_interface_grid_replay_mpi_test).
// The production module tail_interface_grid_design is the only place where
// the monitor integration, face inversion and cylindrical overlap formulas
// live; this header only parses the audit histogram file and formats the
// candidate result blocks.

#include "bulk_tail_moment_audit.h"
#include "grid.h"
#include "tail_interface_grid_design.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct HistogramFileData {
    // Aggregated per (iv, imu): request_number is the summed mass, the
    // request_cell_count is informational (sum over rows of the same bin).
    std::vector<BulkTailVelocityBinAudit> bins;
    long long events;          // unique (accepted_step, conversion_location)
    long long cell_requests;   // sum of request_cell_count over all rows
    int bins_count;            // unique (iv, imu) bins
    double total_number;

    HistogramFileData()
        : events(0), cell_requests(0), bins_count(0), total_number(0.0)
    {}
};

// Raw row of the audit velocity histogram (one accepted step / one
// (iv, imu) group).  Kept row-level so the MPI replay test can partition
// the input across ranks before aggregating.
struct HistogramRow {
    long long accepted_step;
    long long conversion_location;
    int iv;
    int imu;
    long long request_cell_count;
    double request_number;
    HistogramRow()
        : accepted_step(0), conversion_location(0), iv(-1), imu(-1),
          request_cell_count(0), request_number(0.0)
    {}
};

inline bool read_velocity_histogram_rows(const std::string& path,
                                         std::vector<HistogramRow>& rows,
                                         std::string& error)
{
    rows.clear();
    std::ifstream in(path.c_str());
    if (!in) {
        error = "cannot open histogram file " + path;
        return false;
    }
    std::string line;
    bool header_seen = false;
    std::vector<std::pair<std::string, size_t> > seen_keys;
    long long row = 0;
    while (std::getline(in, line)) {
        // Trim surrounding whitespace.
        size_t b = 0;
        while (b < line.size() &&
               std::isspace(static_cast<unsigned char>(line[b]))) ++b;
        size_t e = line.size();
        while (e > b && std::isspace(static_cast<unsigned char>(line[e - 1])))
            --e;
        if (e <= b) continue;
        line = line.substr(b, e - b);
        std::istringstream iss(line);
        if (!header_seen) {
            std::string col;
            std::vector<std::string> cols;
            while (iss >> col) cols.push_back(col);
            if (cols.size() != 7 ||
                cols[0] != "accepted_step" || cols[1] != "time_fs" ||
                cols[2] != "conversion_location" || cols[3] != "iv" ||
                cols[4] != "imu" || cols[5] != "request_cell_count" ||
                cols[6] != "request_number") {
                error = "invalid velocity histogram header";
                return false;
            }
            header_seen = true;
            continue;
        }
        ++row;
        long long step = 0, location = 0, iv = 0, imu = 0, count = 0;
        double time_fs = 0.0;
        double number = 0.0;
        if (!(iss >> step >> time_fs >> location >> iv >> imu >> count >>
              number) ||
            !(number > 0.0) || iv < 0 || imu < 0 || count <= 0) {
            error = "invalid velocity histogram row " + std::to_string(row);
            return false;
        }
        HistogramRow hrow;
        hrow.accepted_step = step;
        hrow.conversion_location = location;
        hrow.iv = static_cast<int>(iv);
        hrow.imu = static_cast<int>(imu);
        hrow.request_cell_count = count;
        hrow.request_number = number;
        rows.push_back(hrow);
    }
    if (!header_seen) {
        error = "empty histogram file (no header)";
        return false;
    }
    return true;
}

inline bool read_velocity_histogram(const std::string& path,
                                    HistogramFileData& data,
                                    std::string& error)
{
    data = HistogramFileData();
    std::vector<HistogramRow> rows;
    if (!read_velocity_histogram_rows(path, rows, error)) return false;
    std::vector<std::pair<std::string, size_t> > seen_keys;
    for (size_t r = 0; r < rows.size(); ++r) {
        const HistogramRow& hrow = rows[r];
        const std::string event_key =
            std::to_string(hrow.accepted_step) + ":" +
            std::to_string(hrow.conversion_location);
        if (seen_keys.empty() || event_key != seen_keys.back().first) {
            seen_keys.push_back(std::make_pair(event_key, 0ULL));
        }
        BulkTailVelocityBinAudit bin;
        bin.iv = hrow.iv;
        bin.imu = hrow.imu;
        bin.request_cell_count =
            static_cast<unsigned long long>(hrow.request_cell_count);
        bin.request_number = hrow.request_number;
        bool merged = false;
        for (size_t q = 0; q < data.bins.size(); ++q) {
            if (data.bins[q].iv == bin.iv && data.bins[q].imu == bin.imu) {
                data.bins[q].request_number += hrow.request_number;
                data.bins[q].request_cell_count += hrow.request_cell_count;
                merged = true;
                break;
            }
        }
        if (!merged) data.bins.push_back(bin);
        data.total_number += hrow.request_number;
        data.cell_requests += hrow.request_cell_count;
    }
    data.events = static_cast<long long>(seen_keys.size());
    data.bins_count = static_cast<int>(data.bins.size());
    std::sort(data.bins.begin(), data.bins.end(),
              [](const BulkTailVelocityBinAudit& a,
                 const BulkTailVelocityBinAudit& b) {
                  if (a.iv != b.iv) return a.iv < b.iv;
                  return a.imu < b.imu;
              });
    return true;
}

inline void write_candidate_block(std::ostream& out,
                                  const std::string& grid_name,
                                  const TailInterfaceReplayResult& r)
{
    const char* comp[] = {"N", "Px", "Jx", "K", "Pixx", "Piperp"};
    out << "grid_name=" << grid_name
        << " candidate_status=" << r.candidate_status
        << " status_reason=" << r.status_reason << "\n";
    out << "number_residual=" << std::setprecision(17) << r.number_residual
        << " max_partition_error=" << r.max_partition_error
        << " negative_mass=" << r.negative_mass_cells
        << " below_threshold_number_fraction="
        << r.below_threshold_number_fraction << "\n";
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        out << "R_L1_" << comp[m] << "=" << std::setprecision(17) << r.r_l1[m]
            << "\n";
    }
    for (int m = 0; m < BULK_TAIL_MOMENT_COUNT; ++m) {
        out << "max_cell_relative_" << comp[m] << "=" << std::setprecision(17)
            << r.max_cell_relative[m] << "\n";
        out << "important_max_cell_relative_" << comp[m] << "="
            << std::setprecision(17) << r.important_max_cell_relative[m]
            << "\n";
        out << "scale_floor_max_" << comp[m] << "=" << std::setprecision(17)
            << r.scale_floor_max[m] << "\n";
    }
    out << "min_dupar=" << std::setprecision(17) << r.min_dupar
        << " min_duperp=" << r.min_duperp
        << " max_adjacent_width_ratio=" << r.max_adjacent_width_ratio
        << " estimated_cell_count_ratio=" << r.estimated_cell_count_ratio
        << " estimated_memory_ratio=" << r.estimated_memory_ratio << "\n";
    out << "estimated_velocity_dt_ratio=" << std::setprecision(17)
        << r.estimated_velocity_dt_ratio << "\n";
    out << "velocity_dt_upar_current=" << r.velocity_dt_upar_current
        << " velocity_dt_upar_candidate=" << r.velocity_dt_upar_candidate
        << " velocity_dt_uperp_current=" << r.velocity_dt_uperp_current
        << " velocity_dt_uperp_candidate=" << r.velocity_dt_uperp_candidate
        << " velocity_dt_total_current=" << r.velocity_dt_total_current
        << " velocity_dt_total_candidate=" << r.velocity_dt_total_candidate
        << "\n";
    out << "estimated_operator_work_ratio=" << std::setprecision(17)
        << r.estimated_operator_work_ratio
        << " scan_cost_ratio=" << r.scan_cost_ratio
        << " pic_creation_ratio=" << r.pic_creation_ratio
        << " collision_pair_ratio=" << r.collision_pair_ratio
        << " collision_substep_ratio=" << r.collision_substep_ratio << "\n";
    out << "center_target_feasible_count=" << r.center_target_feasible_count
        << " center_target_failed_count=" << r.center_target_failed_count
        << " volume_self_feasible_count=" << r.volume_self_feasible_count
        << " volume_self_failed_count=" << r.volume_self_failed_count
        << " volume_self_sparse_failed_count="
        << r.volume_self_sparse_failed_count
        << " max_sparse_support_count=" << r.max_sparse_support_count
        << " estimated_created_macroparticles="
        << r.estimated_created_macroparticles
        << " estimated_particle_ratio_to_center_quartet="
        << r.estimated_particle_ratio_to_center_quartet << "\n";
}

// Write the full per-candidate result file for one run (used by the serial
// and MPI replay tests).  The caller supplies the top-level histogram
// identity lines and the candidate blocks.
inline bool write_result_file(const std::string& path,
                              const std::string& content)
{
    if (path.empty()) {
        std::cout << content;
        return true;
    }
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << content;
    return true;
}

#endif
