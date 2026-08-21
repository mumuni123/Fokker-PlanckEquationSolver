#ifndef OPEN_ELECTROSTATIC_SOLVER_H
#define OPEN_ELECTROSTATIC_SOLVER_H

#include "grid.h"
#include "maxwell.h"

enum class ElectrostaticBoundaryType {
    LEFT_E,
    DIRICHLET_PHI
};

struct ElectrostaticBoundary {
    ElectrostaticBoundaryType type;
    double e_left;
    double phi_left;
    double phi_right;
};

// Production solves only build what the current step needs.  The midpoint
// Poisson solve skips phi reconstruction and the detailed L1/Linf/boundary
// audit; the final accepted state (or an output/checkpoint step) enables the
// full audit.  Defaults describe the full-audit mode.
struct OpenGaussSolveOptions {
    bool reconstruct_phi;
    bool compute_l1;
    bool compute_boundary_audit;
    OpenGaussSolveOptions()
        : reconstruct_phi(true), compute_l1(true), compute_boundary_audit(true) {}
};

struct OpenGaussDiagnostics {
    double residual_l1;
    double residual_linf;
    double boundary_charge_residual;
    double net_charge;
    double mpi_collective_seconds;
};

// Read-only Gate-F identity evaluated with the production face field energy,
// finite-volume cell-average potential and cell charge weights.  It is an
// audit result, never an input to the Poisson solve or a correction to its
// field.
struct OpenPoissonWorkIdentity {
    double field_energy_before;
    double field_energy_after;
    double field_energy_change;
    double electrode_work;
    double potential_charge_work;
    double residual;
    double scale;
    bool finite;
    OpenPoissonWorkIdentity()
        : field_energy_before(0.0), field_energy_after(0.0),
          field_energy_change(0.0), electrode_work(0.0),
          potential_charge_work(0.0), residual(0.0), scale(1.0),
          finite(false) {}
};

class OpenElectrostaticSolver {
public:
    OpenElectrostaticSolver();
    void init(const SpatialGrid& grid, const ElectrostaticBoundary& boundary);
    void solve(EMFields& fields, int mpi_rank, int mpi_size);
    void solve(EMFields& fields, int mpi_rank, int mpi_size,
               const OpenGaussSolveOptions& options);
    void reconstruct_phi(EMFields& fields, int mpi_rank, int mpi_size);
    double boundary_power(double left_current, double right_current) const;
    // Work per unit transverse area supplied by fixed-potential electrodes
    // between two accepted Poisson states.  For Dirichlet data this is
    // phi_L Delta(eps0 E_L) - phi_R Delta(eps0 E_R).  The production
    // configuration phi_L=phi_R=0 returns exactly zero without an MPI call.
    double boundary_energy_work(const EMFields& before,
                                const EMFields& after,
                                int mpi_rank, int mpi_size) const;
    OpenPoissonWorkIdentity evaluate_work_identity(
        const EMFields& before, const EMFields& after,
        const std::vector<double>& rho_delta,
        int mpi_rank, int mpi_size) const;
    // Gate I: construct the face field that is exactly dual to the
    // cell-average potential used by evaluate_work_identity().  This is a
    // read-only diagnostic helper; it never changes either field state.
    // Shared MPI faces contain the same value on both neighboring ranks.
    bool build_potential_pairing_field(
        const EMFields& before, const EMFields& after,
        std::vector<double>& pairing_face,
        int mpi_rank, int mpi_size) const;
    const OpenGaussDiagnostics& diagnostics() const { return diagnostics_; }
    // JC2 (section 5.5.2): read-only face-to-cell helper.  Given Ex_face,
    // computes Ex via cell averaging, syncs MPI face interfaces and exchanges
    // cell ghosts.  Does NOT call solve() or modify Ex_face.
    void populate_electric_components_from_faces(
        EMFields& fields, int mpi_rank, int mpi_size) const;

private:
    SpatialGrid grid_;
    ElectrostaticBoundary boundary_;
    OpenGaussDiagnostics diagnostics_;
    void sync_face_interfaces(EMFields& fields, int rank, int size) const;
    void exchange_cell_ghosts(std::vector<double>& values, int rank, int size,
                              double left_value, double right_value) const;
};

#endif
