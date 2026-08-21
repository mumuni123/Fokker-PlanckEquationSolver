// Stage H2 acceptance: single-cell bulk-to-tail conversion (sections 7.4,
// 7.6, 7.7 and 15 H2).  One cell-integrated mass M at a single
// (u_parallel, u_perp) cell above the threshold must be removed from the
// bulk and replaced by exactly one azimuthal quartet of tail particles that
// preserves N/Px/K exactly and Jx/Pixx/Piperp to roundoff, with the CIC
// density unchanged.
//
// Usage:
//   bulk_tail_single_cell_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "grid.h"
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
    std::uint64_t particles_created;
    double number_residual_rel;
    double px_residual_rel;
    double energy_residual_rel;
    double jx_residual_rel;
    double pixx_residual_rel;
    double piperp_residual_rel;
    double rho_l2_before_after;
    double bulk_number_diff_rel;
    double bulk_energy_diff_rel;
    double tail_density_error_rel;
    double azimuth_structure_ok;
    Metrics()
        : complete(false), conservative(false), fidelity_ok(false),
          finite(false), particles_created(0), number_residual_rel(0.0),
          px_residual_rel(0.0), energy_residual_rel(0.0),
          jx_residual_rel(0.0), pixx_residual_rel(0.0),
          piperp_residual_rel(0.0), rho_l2_before_after(0.0),
          bulk_number_diff_rel(0.0), bulk_energy_diff_rel(0.0),
          tail_density_error_rel(0.0), azimuth_structure_ok(0.0)
    {}
};

bool find_conversion_cell(const HybridVelocityPartition& partition,
                          const CylindricalVelocityGrid& cgrid,
                          int& j_out, int& k_out)
{
    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            if (partition.is_conversion(j, k) &&
                std::fabs(cgrid.upar_cells[j]) > 1.0 &&
                cgrid.uperp_cells[k] > 1.0) {
                j_out = j;
                k_out = k;
                return true;
            }
        }
    }
    return false;
}

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

    int j = -1;
    int k = -1;
    if (!find_conversion_cell(partition, bulk.cgrid, j, k)) {
        std::cerr << "bulk_tail_single_cell_test: no conversion cell with "
                     "|upar|>1 and uperp>1 on this grid.\n";
        return m;
    }
    const double upar = bulk.cgrid.upar_cells[j];
    const double uperp = bulk.cgrid.uperp_cells[k];
    const double gamma = std::sqrt(1.0 + upar * upar + uperp * uperp);
    const double M0 = 1.5e20;
    const int ng = grid.nghost;
    const int ix_local = 0;
    bulk.f[idx3(ng + ix_local, j, k)] = M0;
    bulk.compute_moments();
    const double number_before = bulk.total_particle_number();
    const double energy_before = bulk.total_kinetic_energy();

    Species trial_bulk = bulk;
    BackgroundTailPIC trial_tail;
    trial_tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics d =
        converter.extract_after_substep(trial_bulk, trial_tail, grid,
                                        partition, 7,
                                        ConversionLocation::AFTER_U_SUBSTEP,
                                        0);

    m.complete = d.complete;
    m.conservative = d.conservative;
    m.fidelity_ok = d.fidelity_ok;
    m.finite = d.finite;
    m.particles_created = d.particles_created;
    m.number_residual_rel = d.number_residual_rel;
    m.px_residual_rel = d.px_residual_rel;
    m.energy_residual_rel = d.energy_residual_rel;
    m.jx_residual_rel = d.jx_residual_rel;
    m.pixx_residual_rel = d.pixx_residual_rel;
    m.piperp_residual_rel = d.piperp_residual_rel;
    m.rho_l2_before_after = d.rho_l2_before_after;

    // Bulk side: the mass is gone exactly and the moments drop by M0.
    m.bulk_number_diff_rel =
        std::fabs((number_before - trial_bulk.total_particle_number()) -
                  M0) / std::max(1.0, M0);
    const double expected_energy =
        M0 * (gamma - 1.0) * Const::me * Const::c * Const::c;
    m.bulk_energy_diff_rel =
        std::fabs((energy_before - trial_bulk.total_kinetic_energy()) -
                  expected_energy) / std::max(1.0, expected_energy);

    // Tail side: exactly one quartet at the cell centre with the reference
    // azimuthal structure.
    bool azimuth_ok =
        trial_tail.particles.size() == 4 && m.particles_created == 4;
    const double x_center = 0.5 * grid.dx;
    const double phi[4] = { 0.0, 0.5 * Const::pi, Const::pi, 1.5 * Const::pi };
    for (size_t p = 0; p < trial_tail.particles.size(); ++p) {
        const BackgroundTailParticle& tp = trial_tail.particles[p];
        if (std::fabs(tp.x - x_center) > 1.0e-15) azimuth_ok = false;
        if (std::fabs(tp.ux - upar) > 1.0e-14) azimuth_ok = false;
        if (std::fabs(tp.uy * tp.uy + tp.uz * tp.uz - uperp * uperp) >
            1.0e-10) azimuth_ok = false;
        if (std::fabs(tp.weight - 0.25 * M0) > 1.0e-12 * M0) azimuth_ok = false;
        const double tphi = std::atan2(tp.uz, tp.uy);
        double best = 1.0e30;
        for (int a = 0; a < 4; ++a) {
            double dphi = std::fabs(tphi - phi[a]);
            if (dphi > Const::pi) dphi = 2.0 * Const::pi - dphi;
            best = std::min(best, dphi);
        }
        if (best > 1.0e-6) azimuth_ok = false;
    }
    m.azimuth_structure_ok = azimuth_ok ? 1.0 : 0.0;

    // Density equivalence (section 7.7): the tail deposit reproduces the
    // removed density in the same cell.
    trial_tail.deposit_density(grid, 0, 1);
    const double expected_density = M0 / grid.dx;
    m.tail_density_error_rel =
        std::fabs(trial_tail.density[static_cast<size_t>(ix_local)] -
                  expected_density) /
        std::max(1.0, expected_density);
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
    out << "particles_created=" << m.particles_created << "\n";
    out << "number_residual_rel=" << m.number_residual_rel << "\n";
    out << "px_residual_rel=" << m.px_residual_rel << "\n";
    out << "energy_residual_rel=" << m.energy_residual_rel << "\n";
    out << "jx_residual_rel=" << m.jx_residual_rel << "\n";
    out << "pixx_residual_rel=" << m.pixx_residual_rel << "\n";
    out << "piperp_residual_rel=" << m.piperp_residual_rel << "\n";
    out << "rho_l2_before_after=" << m.rho_l2_before_after << "\n";
    out << "bulk_number_diff_rel=" << m.bulk_number_diff_rel << "\n";
    out << "bulk_energy_diff_rel=" << m.bulk_energy_diff_rel << "\n";
    out << "tail_density_error_rel=" << m.tail_density_error_rel << "\n";
    out << "azimuth_structure_ok=" << m.azimuth_structure_ok << "\n";
    return out.good();
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "bulk_tail_single_cell_test must run with exactly "
                     "1 rank; use plain ./build_hybrid/bulk_tail_"
                     "single_cell_test.\n";
        MPI_Finalize();
        return 2;
    }

    TestArgs args;
    bool ok = parse_args(argc, argv, args);
    if (!ok) {
        std::cerr << "usage: bulk_tail_single_cell_test [--result <path>]\n";
    }

    Metrics m;
    if (ok) m = run_case();
    bool pass = ok &&
                m.complete && m.conservative && m.fidelity_ok &&
                m.finite && m.particles_created == 4 &&
                m.number_residual_rel <= 1.0e-12 &&
                m.px_residual_rel <= 1.0e-12 &&
                m.energy_residual_rel <= 1.0e-12 &&
                m.jx_residual_rel <= 1.0e-12 &&
                m.pixx_residual_rel <= 1.0e-12 &&
                m.piperp_residual_rel <= 1.0e-12 &&
                m.rho_l2_before_after <= 1.0e-12 &&
                m.bulk_number_diff_rel <= 1.0e-12 &&
                m.bulk_energy_diff_rel <= 1.0e-12 &&
                m.tail_density_error_rel <= 1.0e-12 &&
                m.azimuth_structure_ok == 1.0;
    if (!write_result_file(args.result_path, m, pass)) pass = false;
    std::cout << "complete=" << (m.complete ? 1 : 0)
              << " conservative=" << (m.conservative ? 1 : 0)
              << " fidelity_ok=" << (m.fidelity_ok ? 1 : 0)
              << " finite=" << (m.finite ? 1 : 0)
              << " particles_created=" << m.particles_created
              << " number_residual_rel=" << m.number_residual_rel
              << " px_residual_rel=" << m.px_residual_rel
              << " energy_residual_rel=" << m.energy_residual_rel
              << " jx_residual_rel=" << m.jx_residual_rel
              << " pixx_residual_rel=" << m.pixx_residual_rel
              << " piperp_residual_rel=" << m.piperp_residual_rel
              << " rho_l2_before_after=" << m.rho_l2_before_after
              << " bulk_number_diff_rel=" << m.bulk_number_diff_rel
              << " bulk_energy_diff_rel=" << m.bulk_energy_diff_rel
              << " tail_density_error_rel=" << m.tail_density_error_rel
              << " azimuth_structure_ok=" << m.azimuth_structure_ok << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
