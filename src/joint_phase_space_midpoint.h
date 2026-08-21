#ifndef JOINT_PHASE_SPACE_MIDPOINT_H
#define JOINT_PHASE_SPACE_MIDPOINT_H

#include "grid.h"

#include <cstddef>
#include <vector>

struct JointPhaseSpaceIterationRecord {
    int iteration;
    int gmres_dimension;
    double residual_linf;
    double phi_residual_linf;
    double line_search_alpha;
    double trial_min_mass;
    int accepted;
    int failure_code;
    JointPhaseSpaceIterationRecord()
        : iteration(0), gmres_dimension(0), residual_linf(0.0),
          phi_residual_linf(0.0), line_search_alpha(0.0),
          trial_min_mass(0.0), accepted(0), failure_code(0) {}
};

// J0 uses cell-integrated number mass M [m^-2].  The flux arrays below are
// rates [m^-2 s^-1]; multiplying a rate by dt gives the swept mass used by the
// finite-volume update.  The helper is deliberately independent of the
// production Strang/PPM/FCT path and is only the common-flux algebra needed by
// the J0 audit and the later joint midpoint operator.
struct JointPhaseSpaceFluxBundle {
    int nx;
    int nv;
    int nmu;
    double dt;
    double dx;

    // x_flux_rate[(iface * nv + j) * nmu + k].  iface=0 and iface=nx
    // represent the same periodic face in the J0 manufactured test.
    std::vector<double> x_flux_rate;
    // u_flux_rate[((ix * (nv + 1) + jface) * nmu) + k].  The two endpoint
    // faces are exactly zero (zero inflow/outflow ledger for J0).
    std::vector<double> u_flux_rate;
    std::vector<double> charge_current_face;
    std::vector<double> energy_current_cell;
    std::vector<double> mass_delta_x;
    std::vector<double> mass_delta_u;
    std::vector<double> mass_delta_total;
    std::vector<double> kinetic_work_cell;

    JointPhaseSpaceFluxBundle()
        : nx(0), nv(0), nmu(0), dt(0.0), dx(0.0) {}
};

struct JointPhaseSpaceAuditResult {
    bool finite;
    double mass_residual;
    double mass_roundoff_bound;
    double kinetic_work_residual;
    // Sum of absolute production u-face kinetic-work terms.  It is the
    // stable roundoff scale for delta-K versus the u-flux work identity.
    double kinetic_absolute_work_scale;
    double poisson_work_residual;
    double combined_energy_residual;
    double g_gstar_residual;
    double cell_volume_residual;
    double hamiltonian_velocity_residual;
    double u_boundary_flux;

    JointPhaseSpaceAuditResult()
        : finite(false), mass_residual(0.0), mass_roundoff_bound(0.0),
          kinetic_work_residual(0.0), kinetic_absolute_work_scale(0.0),
          poisson_work_residual(0.0), combined_energy_residual(0.0),
          g_gstar_residual(0.0), cell_volume_residual(0.0),
          hamiltonian_velocity_residual(0.0), u_boundary_flux(0.0) {}
};

class JointPhaseSpaceMidpointOperator {
public:
    // Construct the J0 common center-trace fluxes.  M_mid is indexed as
    // [ix * Nv * Nmu + j * Nmu + k] and contains cell-integrated mass.  E_cell
    // is the production face-to-cell gathered field on physical x cells.
    // The x direction is periodic only inside this manufactured J0 helper;
    // this avoids adding an open-boundary policy before J2.
    static JointPhaseSpaceFluxBundle build_periodic_center_flux(
        const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
        const std::vector<double>& m_mid,
        const std::vector<double>& e_cell, double dt);

    // Unique J1 FV residual evaluation for a candidate M^{n+1}.  The field
    // solve is supplied by the caller, while all midpoint mass/flux residual
    // algebra lives here and is never re-derived in the integrator.
    static bool evaluate_residual(
        const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
        const std::vector<double>& m_old,
        const std::vector<double>& m_candidate,
        const std::vector<double>& e_cell, double dt,
        JointPhaseSpaceFluxBundle& bundle,
        std::vector<double>& residual,
        double& residual_linf,
        double& residual_scale);

    static bool apply_block_diagonal_preconditioner(
        const CylindricalVelocityGrid& vg, const std::vector<double>& e_cell,
        double dx, double dt, const std::vector<double>& residual,
        std::vector<double>& preconditioned);

    static bool evaluate_local_residual(
        const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
        const std::vector<double>& m_old_local,
        const std::vector<double>& m_candidate_local,
        const std::vector<double>& e_cell_local, double dt,
        int mpi_rank, int mpi_size,
        JointPhaseSpaceFluxBundle& bundle,
        std::vector<double>& residual,
        double& residual_linf,
        double& residual_scale,
        bool allow_negative_probe = false);

    static bool apply_local_block_diagonal_preconditioner(
        const CylindricalVelocityGrid& vg, const std::vector<double>& e_cell,
        double dx, double dt, const std::vector<double>& residual,
        std::vector<double>& preconditioned);

    // Build the J0 Hamiltonian velocity from the production kinetic-energy
    // table.  Interior cells use the centered discrete derivative of K with
    // respect to p_parallel=m_e*c*u_parallel; the two endpoints use the
    // corresponding one-sided derivative.  This is not the rejected P3-V.2
    // face projection and is not an analytic-vx fallback.
    static std::vector<double> build_hamiltonian_velocity(
        const CylindricalVelocityGrid& vg);

    // Audit a manufactured midpoint update M_old -> M_new whose candidate
    // increment is expected to equal bundle.mass_delta_total.  The Poisson
    // and G/G* residuals are supplied by the production OpenElectrostaticSolver
    // helper in the J0 test; this class does not reimplement that operator.
    static JointPhaseSpaceAuditResult audit(
        const SpatialGrid& sg, const CylindricalVelocityGrid& vg,
        const std::vector<double>& m_old,
        const std::vector<double>& m_new,
        const JointPhaseSpaceFluxBundle& bundle,
        double poisson_work_residual, double g_gstar_residual);

private:
    static std::size_t x_index(int iface, int j, int k,
                               int nv, int nmu);
    static std::size_t u_index(int ix, int jface, int k,
                               int nv, int nmu);
    static std::size_t cell_index(int ix, int j, int k,
                                  int nv, int nmu);
};

#endif
