#ifndef CONSERVATIVE_PPM_REMAP_H
#define CONSERVATIVE_PPM_REMAP_H

#include "grid.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "species.h"
#include "bulk_tail_flux_parcel.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Gate I read-only x-face transport audit (section 4.2).  Forward-declared to
// keep the remap header independent of the audit calculator; the full type is
// included by conservative_ppm_remap.cpp.
struct XFaceTransportAudit;

// Detailed reason for a read-only u_parallel flux-parcel construction
// failure.  These codes refine RemapDiagnostics::audit_failure_code bit 1.
enum class ParcelNodeFailureReason {
    None = 0,
    NegativeQuadratureNode = 1,
    NonfiniteQuadratureNode = 2,
    NonpositiveNodeSum = 3,
    NonfiniteScale = 4,
    InvalidParcelMoments = 5,
    InvalidTransverseMeasure = 6
};

struct ParcelNodeFailure {
    ParcelNodeFailureReason reason;
    double node_mass;
    double target;
    double node_sum;
    double scale;

    ParcelNodeFailure()
        : reason(ParcelNodeFailureReason::None),
          node_mass(std::numeric_limits<double>::quiet_NaN()),
          target(std::numeric_limits<double>::quiet_NaN()),
          node_sum(std::numeric_limits<double>::quiet_NaN()),
          scale(std::numeric_limits<double>::quiet_NaN())
    {}
};

enum class XTransportVelocityMode {
    ANALYTIC_CELL_CENTER = 0,
    ENERGY_CONJUGATE_CELL = 1
};

// Per-call ledger of one conservative remap (section 13.6).  Numbers and
// energies are global sums; tail losses are zero for the x remap and reserved
// for the u_parallel remap.
struct RemapDiagnostics {
    bool finite;
    // Audit validity is deliberately separate from physical remap validity.
    // In FLUX_AUDIT the parcel/quadrature ledger is read-only diagnostics and
    // must not reject an otherwise finite physical remap.
    bool audit_valid;
    int audit_failure_code;
    // First (lowest-rank, then local traversal order) read-only parcel-node
    // construction failure.  This is diagnostic-only and does not alter the
    // physical u_parallel remap.
    int audit_parcel_failure_reason;
    int audit_parcel_failure_rank;
    int audit_parcel_failure_ix;
    int audit_parcel_failure_face;
    int audit_parcel_failure_iuperp;
    double audit_parcel_failure_node_mass;
    double audit_parcel_failure_target;
    double audit_parcel_failure_node_sum;
    double audit_parcel_failure_scale;
    // Set by the u_parallel remap after comparing its final swept-face
    // vector immediately before and after read-only parcel observation.
    bool audit_physical_state_bitwise_equal;
    double number_before;
    double number_after;
    double inflow_number;
    double outflow_number;
    double inflow_energy;
    double outflow_energy;
    // Per-physical-side boundary ledgers (phase-3 addition): the aggregate
    // inflow/outflow above merges both ends; the split quantities feed the
    // open-domain charge/energy ledger (sections 10.1/10.3).
    double left_inflow_number;
    double left_outflow_number;
    double right_inflow_number;
    double right_outflow_number;
    double left_inflow_energy;
    double left_outflow_energy;
    double right_inflow_energy;
    double right_outflow_energy;
    double tail_number_loss;
    double tail_energy_loss;
    double max_departure_cells;
    double constant_fraction;
    double linear_fraction;
    double limited_fraction;
    std::uint64_t limited_cell_count;
    double minimum_cell_mass;
    double interface_face_export_number;
    double interface_export_number;
    double interface_parcel_number;
    double interface_export_energy;
    std::uint64_t interface_parcel_count;
    std::uint64_t interface_node_count;
    std::uint64_t interface_duplicate_count;
    double interface_below_threshold_number;
    double interface_roundoff_discarded_number;
    double interface_quadrature_error_max;
    // Tail-owned cell mismatch below the same per-interface finite-volume
    // resolution used to suppress unrepresentable PIC sources.
    double tail_owned_roundoff_discarded_number;
    // Amount injected into tail-owned Eulerian cells by the same accepted
    // bulk-to-tail face transfer used to create PIC parcels.
    double tail_owned_expected_transfer_number;
    // Sum of per-connected-tail-region absolute signed-balance errors after
    // subtracting the expected interface transfer above.  Internal tail
    // faces telescope and are intentionally excluded.
    double tail_owned_bulk_residual;

    // Failure-localization data.  The extrema are over the physical cells of
    // this remap; indices identify the first non-finite value when present.
    double input_min_mass;
    double input_max_mass;
    double output_min_mass;
    double output_max_mass;
    double first_nonfinite_value;
    int first_nonfinite_rank;
    int first_nonfinite_ix;
    int first_nonfinite_iupar;
    int first_nonfinite_iuperp;

    // Gate C (section 7.3): bulk u_parallel discrete kinetic-energy identity.
    // These describe only the final upar_swept_ actually used by the m_new
    // update (i.e. after the flux-interface sink), never a re-reconstructed
    // higher-order candidate flux.  All are global sums inside advect_u_parallel.
    double upar_internal_face_energy_transfer;
    double upar_left_velocity_boundary_energy;
    double upar_right_velocity_boundary_energy;
    double upar_interface_energy_removed;
    double upar_discrete_energy_identity_residual;

    RemapDiagnostics()
        : finite(false), audit_valid(true), audit_failure_code(0),
          audit_parcel_failure_reason(
              static_cast<int>(ParcelNodeFailureReason::None)),
          audit_parcel_failure_rank(-1), audit_parcel_failure_ix(-1),
          audit_parcel_failure_face(-1), audit_parcel_failure_iuperp(-1),
          audit_parcel_failure_node_mass(std::numeric_limits<double>::quiet_NaN()),
          audit_parcel_failure_target(std::numeric_limits<double>::quiet_NaN()),
          audit_parcel_failure_node_sum(std::numeric_limits<double>::quiet_NaN()),
          audit_parcel_failure_scale(std::numeric_limits<double>::quiet_NaN()),
          audit_physical_state_bitwise_equal(true),
          number_before(0.0), number_after(0.0),
          inflow_number(0.0), outflow_number(0.0), inflow_energy(0.0),
          outflow_energy(0.0), left_inflow_number(0.0),
          left_outflow_number(0.0), right_inflow_number(0.0),
          right_outflow_number(0.0), left_inflow_energy(0.0),
          left_outflow_energy(0.0), right_inflow_energy(0.0),
          right_outflow_energy(0.0), tail_number_loss(0.0),
           tail_energy_loss(0.0), max_departure_cells(0.0),
           constant_fraction(0.0), linear_fraction(0.0),
           limited_fraction(0.0), limited_cell_count(0),
          minimum_cell_mass(0.0), interface_face_export_number(0.0),
          interface_export_number(0.0), interface_parcel_number(0.0),
          interface_export_energy(0.0), interface_parcel_count(0),
          interface_node_count(0), interface_duplicate_count(0),
          interface_below_threshold_number(0.0),
          interface_roundoff_discarded_number(0.0),
          interface_quadrature_error_max(0.0),
          tail_owned_roundoff_discarded_number(0.0),
          tail_owned_expected_transfer_number(0.0),
          tail_owned_bulk_residual(0.0),
          input_min_mass(std::numeric_limits<double>::infinity()),
          input_max_mass(-std::numeric_limits<double>::infinity()),
          output_min_mass(std::numeric_limits<double>::infinity()),
          output_max_mass(-std::numeric_limits<double>::infinity()),
          first_nonfinite_value(std::numeric_limits<double>::quiet_NaN()),
          first_nonfinite_rank(-1), first_nonfinite_ix(-1),
          first_nonfinite_iupar(-1), first_nonfinite_iuperp(-1),
          upar_internal_face_energy_transfer(0.0),
          upar_left_velocity_boundary_energy(0.0),
          upar_right_velocity_boundary_energy(0.0),
           upar_interface_energy_removed(0.0),
           upar_discrete_energy_identity_residual(0.0)
    {}
};

// Persistent scratch owned by ConservativePpmRemap (never allocated per
// velocity slice).  average/left_edge/right_edge/curvature/cumulative_mass
// are sized over one extended x slice (nx_local + 2*halo cells), swept_mass
// over the nx_local+1 local faces, halo over 2*halo*nvmu doubles.
struct PpmWorkspace {
    std::vector<double> average;
    std::vector<double> left_edge;
    std::vector<double> right_edge;
    std::vector<double> curvature;
    std::vector<double> cumulative_mass;
    std::vector<double> swept_mass;
    std::vector<double> halo;
};

// Conservative PPM semi-Lagrangian remap (section 13.6).  The x remap solves
// dt f + vx dx f = 0 per velocity cell; the u_parallel remap solves
// dt f + a_u du_parallel f = 0 per (x, u_perp) slice on the non-uniform
// u_parallel grid.  Both follow the grid passed to init(); the u_parallel
// path reads the velocity extent and cell geometry from that grid at runtime
// (needed by the 48/96/192-cell convergence test, section 16.4.2).
class ConservativePpmRemap {
public:
    ConservativePpmRemap();

    void set_x_transport_velocity_mode(XTransportVelocityMode mode);
    XTransportVelocityMode x_transport_velocity_mode() const
    { return x_transport_velocity_mode_; }

    void init(const SpatialGrid& grid, const CylindricalVelocityGrid& velocity_grid);

    // Solves dt f + vx dx f = 0 per velocity cell.  Reads cell-integrated mass
    // from input, writes output (in-place is safe), refreshes output ghosts
    // through the boundary object and recomputes output moments.  When the
    // Gate I audit pointer is non-NULL, the finally-applied face swept number
    // (m^-2, summed over velocity, not multiplied by charge, not divided by
    // dt) is accumulated into audit->bulk_number_swept_face; the caller owns
    // and clears the audit object (section 4.2).  A NULL audit adds no
    // allocation, no clearing and no reduction, so the physical state is
    // bitwise unchanged.
    RemapDiagnostics advect_x(
        const Species& input,
        Species& output,
        double dt,
        double time,
        const OpenBackgroundBoundary& boundary,
        int mpi_rank,
        int mpi_size,
        XFaceTransportAudit* audit = NULL);

    // Solves dt f + a_u du_parallel f = 0 with a_u = q E_x / (m c) on the
    // non-uniform u_parallel grid.  The departure coordinate is located with
    // the actual face coordinates (section 6.3); both velocity ends default
    // to zero inflow, and mass/energy leaving the domain enter the tail
    // ledger.  Every x column (interior and ghost) is transformed by the same
    // column operator so MPI ghost cells stay consistent; the caller is
    // responsible for refreshing physical-boundary ghosts and moments.
    // In-place (output == input) is safe.  When local_delta_ke_by_x is
    // non-NULL, the Gate I per-cell bulk discrete kinetic-energy change
    // (J/m^2, section 4.4) is accumulated into (*local_delta_ke_by_x)[ix]
    // (nx_local entries, owned by the caller); a NULL pointer adds no scan,
    // no sqrt and no reduction beyond the existing production path.
    RemapDiagnostics advect_u_parallel(
        const Species& input,
        Species& output,
        const EMFields& midpoint_fields,
        double dt,
        double time,
        const HybridVelocityPartition* partition = NULL,
        BulkTailFluxBatch* exported_flux = NULL,
        int quadrature_order = 4,
        std::vector<double>* local_delta_ke_by_x = NULL);

private:
    SpatialGrid grid_;
    const CylindricalVelocityGrid* cgrid_;
    PpmWorkspace work_;
    std::vector<double> send_left_halo_;
    std::vector<double> send_right_halo_;
    std::vector<double> right_boundary_flux_;
    std::vector<double> left_boundary_flux_received_;
    // Persistent u_parallel workspace (never allocated per slice).
    std::vector<double> upar_average_;
    std::vector<double> upar_left_edge_;
    std::vector<double> upar_right_edge_;
    std::vector<double> upar_curvature_;
    std::vector<double> upar_cumulative_;
    std::vector<double> upar_swept_;
    // Per-face 4-cell stencil for the non-uniform PPM interface value.
    std::vector<double> upar_face_coeffs_;
    std::vector<int> upar_stencil_base_;
    int upar_nv_;
    int upar_nmu_;
    double max_abs_vx_;
    int halo_width_;
    int nx_local_;
    int ng_;
    XTransportVelocityMode x_transport_velocity_mode_;

    double parabola_integral(int extended_cell, double s1, double s2) const;
    double swept_mass(int local_face, double vx, double dt) const;
    double transport_velocity(size_t q) const;
    void build_slice_reconstruction(long long& constant_cells,
                                    long long& linear_cells,
                                    long long& limited_cells);
    void exchange_halo(const Species& input, int mpi_rank, int mpi_size,
                       double dt);
    void exchange_boundary_fluxes(int mpi_rank, int mpi_size);
    void fill_physical_halo(PhysicalSide side, const Species& input,
                            int j_upar, int k_uperp, double vx, double time,
                            const OpenBackgroundBoundary& boundary);
    // u_parallel helpers (section 6.3): all geometry comes from cgrid_.
    size_t upar_index(int ix, int j_upar, int k_uperp) const;
    int locate_upar(double u) const;
    double upar_parabola_integral(int cell, double u1, double u2) const;
    double upar_swept_mass(int face, double a_u, double dt,
                           int& cells_spanned) const;
    double upar_profile_value(int cell, double u) const;
    bool append_upar_swept_nodes(int face, double a_u, double dt,
                                 int k_uperp, int quadrature_order,
                                 BulkTailFluxParcel& parcel,
                                 ParcelNodeFailure* failure) const;
    double upar_tail_energy(int face, double a_u, double dt,
                            int k_uperp, double particle_mass) const;
    void build_upar_reconstruction(long long& constant_cells,
                                   long long& linear_cells, bool count);
    void compute_upar_face_coefficients(const CylindricalVelocityGrid& grid);
};

#endif
