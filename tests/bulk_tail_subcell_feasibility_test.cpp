#include "grid.h"
#include "species.h"
#include "tail_moment_constraint.h"
#include "tail_subcell_quadrature.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Args {
    std::string test_case;
    std::string result;
    Args() : test_case("all") {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) args.test_case = argv[++i];
        else if (arg == "--result" && i + 1 < argc) args.result = argv[++i];
        else return false;
    }
    return args.test_case == "all" ||
           args.test_case == "self-moments" ||
           args.test_case == "center-augmented" ||
           args.test_case == "constraint-hierarchy";
}

std::vector<double> moment_column(double upar, double uperp)
{
    TailMoment7 m;
    tail_particle_moments(1.0, 0.0, upar, uperp, 0.0, m);
    std::vector<double> result(6, 0.0);
    result[0] = m.n;
    result[1] = m.px;
    result[2] = m.jx;
    result[3] = m.ke;
    result[4] = m.pixx;
    result[5] = m.piperp;
    return result;
}

std::vector<double> weighted_moments(
    const std::vector<std::vector<double> >& columns,
    const std::vector<double>& weights)
{
    const size_t rows = columns.empty() ? 0 : columns[0].size();
    std::vector<double> result(rows, 0.0);
    for (size_t q = 0; q < columns.size(); ++q) {
        for (size_t r = 0; r < rows; ++r) {
            result[r] += columns[q][r] * weights[q];
        }
    }
    return result;
}

double max_relative_residual(
    const std::vector<std::vector<double> >& columns,
    const std::vector<double>& weights,
    const std::vector<double>& reference)
{
    const std::vector<double> got = weighted_moments(columns, weights);
    double result = 0.0;
    for (size_t r = 0; r < reference.size(); ++r) {
        result = std::max(result,
            std::fabs(got[r] - reference[r]) /
            std::max(1.0, std::fabs(reference[r])));
    }
    return result;
}

double offcenter_fraction(const std::vector<double>& weights,
                          size_t subcell_count)
{
    double total = 0.0;
    double offcenter = 0.0;
    for (size_t q = 0; q < weights.size(); ++q) {
        total += std::max(0.0, weights[q]);
        if (q < subcell_count) offcenter += std::max(0.0, weights[q]);
    }
    return offcenter / std::max(1.0, total);
}

void select_rows(const std::vector<std::vector<double> >& full_columns,
                 const std::vector<double>& full_reference,
                 const std::vector<int>& rows,
                 std::vector<std::vector<double> >& columns,
                 std::vector<double>& reference)
{
    columns.assign(full_columns.size(), std::vector<double>(rows.size(), 0.0));
    reference.assign(rows.size(), 0.0);
    for (size_t r = 0; r < rows.size(); ++r) {
        reference[r] = full_reference[static_cast<size_t>(rows[r])];
        for (size_t q = 0; q < full_columns.size(); ++q) {
            columns[q][r] = full_columns[q][static_cast<size_t>(rows[r])];
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    Args args;
    if (ranks != 1 || !parse_args(argc, argv, args)) {
        std::cerr << "usage: bulk_tail_subcell_feasibility_test "
                     "--case all|self-moments|center-augmented|constraint-hierarchy "
                     "[--result path]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 1, Param::dx);
    Species bulk;
    bulk.init("subcell_feasibility", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);

    int iv = -1;
    int imu = -1;
    std::vector<TailSubcellNode> nodes;
    const double window_high =
        partition.min_conversion_energy + 0.2e6 * Const::eV;
    for (int k = 0; k < std::min(8, Param::Nmu) && iv < 0; ++k) {
        for (int j = 0; j < Param::Nv; ++j) {
            if (!partition.is_conversion(j, k)) continue;
            const std::vector<TailSubcellNode> candidate =
                TailSubcellQuadrature::nodes(bulk.cgrid, j, k);
            bool all_above = !candidate.empty();
            bool reaches_window = false;
            for (size_t q = 0; q < candidate.size(); ++q) {
                all_above = all_above &&
                    candidate[q].kinetic_energy >= partition.min_conversion_energy;
                reaches_window = reaches_window ||
                    candidate[q].kinetic_energy < window_high;
            }
            if (all_above && reaches_window) {
                iv = j;
                imu = k;
                nodes = candidate;
                break;
            }
        }
    }

    bool harness_pass = iv >= 0 && !nodes.empty();
    std::ostringstream report;
    report << std::setprecision(17);
    if (!harness_pass) {
        report << "status=FAIL reason=no-near-axis-candidate\n";
    } else {
        const double mass = 1.0e20;
        std::vector<std::vector<double> > subcell_columns;
        std::vector<double> subcell_prior;
        double min_energy = 1.0e300;
        double max_energy = 0.0;
        for (size_t q = 0; q < nodes.size(); ++q) {
            subcell_columns.push_back(
                moment_column(nodes[q].upar, nodes[q].uperp));
            subcell_prior.push_back(mass * nodes[q].mass_fraction);
            min_energy = std::min(min_energy, nodes[q].kinetic_energy);
            max_energy = std::max(max_energy, nodes[q].kinetic_energy);
        }
        const std::vector<double> subcell_target =
            weighted_moments(subcell_columns, subcell_prior);
        const std::vector<double> center_unit = moment_column(
            bulk.cgrid.upar_cells[iv], bulk.cgrid.uperp_cells[imu]);
        std::vector<double> center_target(center_unit.size(), 0.0);
        for (size_t r = 0; r < center_unit.size(); ++r) {
            center_target[r] = mass * center_unit[r];
        }

        report << "cell iv=" << iv << " imu=" << imu
               << " candidate_count=" << nodes.size()
               << " center_energy_mev="
               << bulk.cgrid.kinetic_energy[idx2(iv, imu)] /
                    (1.0e6 * Const::eV)
               << " candidate_energy_min_mev="
               << min_energy / (1.0e6 * Const::eV)
               << " candidate_energy_max_mev="
               << max_energy / (1.0e6 * Const::eV) << "\n";

        if (args.test_case == "all" || args.test_case == "self-moments") {
            std::vector<double> weights;
            const bool feasible = tail_solve_nonnegative_moment_weights(
                subcell_columns, subcell_target, subcell_prior, weights, 1.0e-10);
            const double residual = weights.empty() ? 1.0e300 :
                max_relative_residual(subcell_columns, weights, subcell_target);
            report << "self_moments feasible=" << (feasible ? 1 : 0)
                   << " max_relative_residual=" << residual << "\n";
            harness_pass = harness_pass && feasible && residual <= 1.0e-10;
        }

        if (args.test_case == "all" || args.test_case == "center-augmented") {
            std::vector<std::vector<double> > augmented = subcell_columns;
            augmented.push_back(center_unit);

            std::vector<double> center_seed(nodes.size() + 1, 0.0);
            center_seed.back() = mass;
            std::vector<double> center_weights;
            const bool center_seed_feasible =
                tail_solve_nonnegative_moment_weights(
                    augmented, center_target, center_seed, center_weights, 1.0e-10);
            const double center_seed_residual = center_weights.empty() ? 1.0e300 :
                max_relative_residual(augmented, center_weights, center_target);

            std::vector<double> spread_seed = subcell_prior;
            spread_seed.push_back(0.0);
            std::vector<double> spread_weights;
            const bool spread_seed_feasible =
                tail_solve_nonnegative_moment_weights(
                    augmented, center_target, spread_seed, spread_weights, 1.0e-10);
            const double spread_seed_residual = spread_weights.empty() ? 1.0e300 :
                max_relative_residual(augmented, spread_weights, center_target);
            const double spread_offcenter = spread_weights.empty() ? -1.0 :
                offcenter_fraction(spread_weights, nodes.size());

            report << "center_augmented_center_seed feasible="
                   << (center_seed_feasible ? 1 : 0)
                   << " max_relative_residual=" << center_seed_residual
                   << " offcenter_weight_fraction="
                   << (center_weights.empty() ? -1.0 :
                       offcenter_fraction(center_weights, nodes.size())) << "\n";
            report << "center_augmented_spread_seed feasible="
                   << (spread_seed_feasible ? 1 : 0)
                   << " max_relative_residual=" << spread_seed_residual
                   << " offcenter_weight_fraction=" << spread_offcenter << "\n";
            // The exact-center seed is a known feasible control.  The spread
            // result is diagnostic and may legitimately be infeasible.
            harness_pass = harness_pass && center_seed_feasible &&
                           center_seed_residual <= 1.0e-10;
        }

        if (args.test_case == "all" ||
            args.test_case == "constraint-hierarchy") {
            const int row_sets[][6] = {
                {0, -1, -1, -1, -1, -1},
                {0, 1, -1, -1, -1, -1},
                {0, 1, 3, -1, -1, -1},
                {0, 1, 3, 2, -1, -1},
                {0, 1, 3, 2, 4, -1},
                {0, 1, 3, 2, 4, 5}
            };
            const char* labels[] = {
                "N", "N_Px", "N_Px_K", "N_Px_K_Jx",
                "N_Px_K_Jx_Pixx", "N_Px_K_Jx_Pixx_Piperp"
            };
            for (int level = 0; level < 6; ++level) {
                std::vector<int> rows;
                for (int r = 0; r < 6 && row_sets[level][r] >= 0; ++r) {
                    rows.push_back(row_sets[level][r]);
                }
                std::vector<std::vector<double> > columns;
                std::vector<double> reference;
                select_rows(subcell_columns, center_target, rows,
                            columns, reference);
                std::vector<double> weights;
                const bool feasible = tail_solve_nonnegative_moment_weights(
                    columns, reference, subcell_prior, weights, 1.0e-10);
                const double residual = weights.empty() ? 1.0e300 :
                    max_relative_residual(columns, weights, reference);
                report << "hierarchy=" << labels[level]
                       << " feasible=" << (feasible ? 1 : 0)
                       << " max_relative_residual=" << residual << "\n";
            }
        }

        report << "status=" << (harness_pass ? "PASS" : "FAIL") << "\n";
    }

    std::cout << report.str();
    if (!args.result.empty()) {
        std::ofstream out(args.result.c_str(), std::ios::trunc);
        if (!out) {
            std::cerr << "cannot write result file\n";
            MPI_Finalize();
            return 3;
        }
        out << report.str();
    }
    MPI_Finalize();
    return harness_pass ? 0 : 1;
}
