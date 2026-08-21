// Stage H2 acceptance: conversion leaves the discrete space charge unchanged
// (sections 4.3, 7.7 and 15 H2).  A bulk state with a high-energy packet
// above the threshold is converted; the combined bulk+tail electron density
// and the resulting non-periodic Poisson field must match the pre-conversion
// state to deposition roundoff, and N/K must be conserved across the two
// representations.
//
// Usage:
//   bulk_tail_poisson_invariance_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "grid.h"
#include "maxwell.h"
#include "open_electrostatic_solver.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

struct Metrics {
    bool complete;
    bool conservative;
    bool fidelity_ok;
    bool finite;
    double jx_residual_rel;
    double pixx_residual_rel;
    double piperp_residual_rel;
    std::uint64_t fallback_count;
    double density_l2_rel;
    double density_linf_rel;
    double field_l2_rel;
    double field_linf_rel;
    double number_conservation_rel;
    double energy_conservation_rel;
    Metrics()
        : complete(false), conservative(false), fidelity_ok(false),
          finite(false), jx_residual_rel(0.0), pixx_residual_rel(0.0),
          piperp_residual_rel(0.0), fallback_count(0),
          density_l2_rel(0.0), density_linf_rel(0.0),
          field_l2_rel(0.0), field_linf_rel(0.0),
          number_conservation_rel(0.0), energy_conservation_rel(0.0)
    {}
};

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 12, 1.2 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);

    // Bulk state: a uniform core below the threshold plus a high-energy
    // packet in the conversion cells of local cells 2..5.
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
    for (int il = 0; il < nxl; ++il) {
        const int ixg = grid.ix_start + il;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                if (partition.is_conversion(j, k)) continue;
                bulk.f[idx3(ng + il, j, k)] =
                    1.0e18 * (1.0 + static_cast<double>(ixg % 2));
            }
        }
        if (il >= 2 && il <= 5) {
            for (int j = 0; j < Param::Nv; ++j) {
                for (int k = 0; k < Param::Nmu; ++k) {
                    if (!partition.is_conversion(j, k)) continue;
                    bulk.f[idx3(ng + il, j, k)] =
                        5.0e19 * (1.0 + static_cast<double>(ixg % 2));
                }
            }
        }
    }
    bulk.compute_moments();
    const double number_before = bulk.total_particle_number();
    const double energy_before = bulk.total_kinetic_energy();

    std::vector<double> ions(static_cast<size_t>(nxl), Param::dens);
    EMFields fields_before;
    fields_before.init(grid);
    for (int il = 0; il < nxl; ++il) {
        fields_before.rho[static_cast<size_t>(ng + il)] =
            Const::qe * (ions[static_cast<size_t>(il)] -
                         bulk.number_density[static_cast<size_t>(il)]);
    }
    OpenElectrostaticSolver solver;
    solver.init(grid, { ElectrostaticBoundaryType::LEFT_E, 0.0, 0.0, 0.0 });
    solver.solve(fields_before, 0, 1);

    // Convert on trial states.
    Species trial_bulk = bulk;
    BackgroundTailPIC trial_tail;
    trial_tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d =
        converter.extract_after_substep(trial_bulk, trial_tail, grid,
                                        partition, 3,
                                        ConversionLocation::AFTER_U_SUBSTEP,
                                        0);
    trial_tail.deposit_density(grid, 0, 1);
    m.complete = d.complete;
    m.conservative = d.conservative;
    m.fidelity_ok = d.fidelity_ok;
    m.finite = d.finite;
    m.jx_residual_rel = d.jx_residual_rel;
    m.pixx_residual_rel = d.pixx_residual_rel;
    m.piperp_residual_rel = d.piperp_residual_rel;
    m.fallback_count = d.compression_fallback_count;

    // Combined electron density must be unchanged (section 4.3 / 7.7).
    double l2_before = 0.0;
    double l2_diff = 0.0;
    double linf_diff = 0.0;
    for (int il = 0; il < nxl; ++il) {
        const double ne_before = bulk.number_density[static_cast<size_t>(il)];
        const double ne_after =
            trial_bulk.number_density[static_cast<size_t>(il)] +
            trial_tail.density[static_cast<size_t>(il)];
        const double diff = ne_before - ne_after;
        l2_before += ne_before * ne_before;
        l2_diff += diff * diff;
        linf_diff = std::max(linf_diff, std::fabs(diff));
    }
    m.density_l2_rel = std::sqrt(l2_diff) /
                       std::max(1.0, std::sqrt(l2_before));
    m.density_linf_rel = linf_diff / std::max(
        1.0, l2_before / static_cast<double>(nxl) *
             static_cast<double>(nxl) + 1.0);

    // The Poisson field must be unchanged because rho is unchanged.
    EMFields fields_after;
    fields_after.init(grid);
    for (int il = 0; il < nxl; ++il) {
        fields_after.rho[static_cast<size_t>(ng + il)] =
            Const::qe * (ions[static_cast<size_t>(il)] -
                         trial_bulk.number_density[static_cast<size_t>(il)] -
                         trial_tail.density[static_cast<size_t>(il)]);
    }
    solver.solve(fields_after, 0, 1);
    double e2_before = 0.0;
    double e2_diff = 0.0;
    double elinf_diff = 0.0;
    for (size_t f = 0; f < fields_before.Ex_face.size(); ++f) {
        const double eb = fields_before.Ex_face[f];
        const double ea = fields_after.Ex_face[f];
        const double diff = eb - ea;
        e2_before += eb * eb;
        e2_diff += diff * diff;
        elinf_diff = std::max(elinf_diff, std::fabs(diff));
    }
    m.field_l2_rel = std::sqrt(e2_diff) / std::max(1.0, std::sqrt(e2_before));
    m.field_linf_rel = elinf_diff / std::max(1.0, std::sqrt(e2_before) /
        std::max(1.0, static_cast<double>(fields_before.Ex_face.size())));

    // Combined N/K conservation across the two representations.
    double tail_weight = 0.0;
    for (size_t i = 0; i < trial_tail.particles.size(); ++i) {
        tail_weight += trial_tail.particles[i].weight;
    }
    const double number_after =
        trial_bulk.total_particle_number() + tail_weight;
    double tail_energy = 0.0;
    for (size_t i = 0; i < trial_tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = trial_tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        tail_energy += p.weight * (gamma - 1.0) *
                       Const::me * Const::c * Const::c;
    }
    m.number_conservation_rel =
        std::fabs(number_after - number_before) /
        std::max(1.0, number_before);
    m.energy_conservation_rel =
        std::fabs((trial_bulk.total_kinetic_energy() + tail_energy) -
                  energy_before) / std::max(1.0, energy_before);
    return m;
}

bool write_result_file(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "complete=" << (m.complete ? 1 : 0) << "\n";
    out << "conservative=" << (m.conservative ? 1 : 0) << "\n";
    out << "fidelity_ok=" << (m.fidelity_ok ? 1 : 0) << "\n";
    out << "finite=" << (m.finite ? 1 : 0) << "\n";
    out << "jx_residual_rel=" << m.jx_residual_rel << "\n";
    out << "pixx_residual_rel=" << m.pixx_residual_rel << "\n";
    out << "piperp_residual_rel=" << m.piperp_residual_rel << "\n";
    out << "fallback_count=" << m.fallback_count << "\n";
    out << "density_l2_rel=" << m.density_l2_rel << "\n";
    out << "density_linf_rel=" << m.density_linf_rel << "\n";
    out << "field_l2_rel=" << m.field_l2_rel << "\n";
    out << "field_linf_rel=" << m.field_linf_rel << "\n";
    out << "number_conservation_rel=" << m.number_conservation_rel << "\n";
    out << "energy_conservation_rel=" << m.energy_conservation_rel << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "bulk_tail_poisson_invariance_test must run with "
                     "exactly 1 rank; use plain ./build_hybrid/bulk_tail_"
                     "poisson_invariance_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: bulk_tail_poisson_invariance_test "
                     "[--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case();
    bool pass = ok && m.complete && m.conservative && m.fidelity_ok &&
                m.finite && m.density_l2_rel <= 1.0e-12 &&
                m.field_l2_rel <= 1.0e-12 &&
                m.number_conservation_rel <= 1.0e-12 &&
                m.energy_conservation_rel <= 1.0e-12;
    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "complete=" << (m.complete ? 1 : 0)
              << " conservative=" << (m.conservative ? 1 : 0)
              << " fidelity_ok=" << (m.fidelity_ok ? 1 : 0)
              << " finite=" << (m.finite ? 1 : 0)
              << " jx_residual_rel=" << m.jx_residual_rel
              << " pixx_residual_rel=" << m.pixx_residual_rel
              << " piperp_residual_rel=" << m.piperp_residual_rel
              << " fallback_count=" << m.fallback_count
              << " density_l2_rel=" << m.density_l2_rel
              << " density_linf_rel=" << m.density_linf_rel
              << " field_l2_rel=" << m.field_l2_rel
              << " field_linf_rel=" << m.field_linf_rel
              << " number_conservation_rel=" << m.number_conservation_rel
              << " energy_conservation_rel=" << m.energy_conservation_rel
              << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
