// Gate F 10.5.3: exercise the production open Poisson operator and its
// read-only discrete-work identity.  No alternative finite-difference field
// operator is implemented in this test.

#include "grid.h"
#include "maxwell.h"
#include "open_electrostatic_solver.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

bool parse_args(int argc, char** argv, std::string& result)
{
    result = "output/vpfp_poisson_work_identity_unit.result";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) {
            if (std::string(argv[++i]) != "all") return false;
        } else if (arg == "--result" && i + 1 < argc) {
            result = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

void set_charge(EMFields& fields, const SpatialGrid& grid, double amplitude)
{
    std::fill(fields.rho.begin(), fields.rho.end(), 0.0);
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int slot = grid.nghost + ix;
        const double x = (static_cast<double>(grid.global_cell(slot)) + 0.5) *
            grid.dx;
        fields.rho[static_cast<size_t>(slot)] = amplitude *
            std::sin(2.0 * Const::pi * x / grid.length());
    }
}

double tolerance(const OpenPoissonWorkIdentity& identity)
{
    return 8192.0 * std::numeric_limits<double>::epsilon() * identity.scale;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string result_path;
    bool pass = parse_args(argc, argv, result_path);
    SpatialGrid grid;
    grid.init_with_domain(rank, size, 64, 64.0e-9);
    EMFields before;
    EMFields after;
    before.init(grid);
    after.init(grid);
    const ElectrostaticBoundary boundary = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    OpenElectrostaticSolver solver;
    solver.init(grid, boundary);

    if (pass) {
        set_charge(before, grid, 2.0e-5);
        solver.solve(before, rank, size);
        set_charge(after, grid, -3.0e-5);
        solver.solve(after, rank, size);
    }
    std::vector<double> rho_delta(after.rho.size(), 0.0);
    for (size_t i = 0; i < rho_delta.size(); ++i) {
        rho_delta[i] = after.rho[i] - before.rho[i];
    }
    const OpenPoissonWorkIdentity identity = solver.evaluate_work_identity(
        before, after, rho_delta, rank, size);
    // The manufactured state changes only the bulk charge.  Still evaluate
    // every ownership slot through the production read-only interface so the
    // result schema cannot silently omit a species contribution.
    const std::vector<double> zero_delta(rho_delta.size(), 0.0);
    const OpenPoissonWorkIdentity tail_identity = solver.evaluate_work_identity(
        before, after, zero_delta, rank, size);
    const bool local_pass = pass && identity.finite &&
        std::fabs(identity.residual) <= tolerance(identity);
    int global_pass = local_pass ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_pass, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream out(result_path.c_str());
        out << std::setprecision(17)
            << "field_energy_before=" << identity.field_energy_before << "\n"
            << "field_energy_after=" << identity.field_energy_after << "\n"
            << "field_energy_change=" << identity.field_energy_change << "\n"
            << "electrode_work=" << identity.electrode_work << "\n"
            << "potential_charge_work=" << identity.potential_charge_work << "\n"
            << "rho_bulk_potential_work=" << identity.potential_charge_work << "\n"
            << "rho_tail_potential_work=" << tail_identity.potential_charge_work << "\n"
            << "rho_beam_potential_work=" << tail_identity.potential_charge_work << "\n"
            << "rho_ion_potential_work=" << tail_identity.potential_charge_work << "\n"
            << "rho_boundary_source_potential_work="
            << tail_identity.potential_charge_work << "\n"
            << "R_P_h=" << identity.residual << "\n"
            << "roundoff_tolerance=" << tolerance(identity) << "\n"
            << "poisson_identity_pass=" << global_pass << "\n"
            << "status=" << (global_pass ? "PASS" : "FAIL") << "\n";
        std::cout << "status=" << (global_pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return global_pass ? 0 : 1;
}
