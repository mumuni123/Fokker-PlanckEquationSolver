#ifndef BACKGROUND_TAIL_PIC_H
#define BACKGROUND_TAIL_PIC_H

#include "grid.h"
#include "maxwell.h"

#include <cstddef>
#include <cstdint>
#include <vector>

// Background high-energy tail PIC infrastructure (section 6 of the open-Beam
// VPFP-PIC reconstruction plan).  H1 implements the standalone class only:
// it is NOT connected to the converter or the production integrator until
// H2/H3.
//
// Physical units follow section 4.1:
//   x        : metres
//   u        : p/(m_e c), dimensionless
//   weight   : electrons per unit transverse area [m^-2]
//   density  : number density [m^-3] (deposited weight divided by dx)
//   current_face_x : current density [A/m^2] = -q_e times number flux
//   outflow px / kinetic energy : per unit transverse area
//
// CIC is used for the first version (section 6.3): the density cell weights
// and the shifted face gather are numerically identical to the already
// verified BeamPIC deposit/gather formulas (constraint 6), so no Beam
// regression is required in H1.

struct BackgroundTailParticle {
    double x;       // m
    double ux;      // p_x/(m_e c)
    double uy;      // p_y/(m_e c)
    double uz;      // p_z/(m_e c)
    double weight;  // m^-2
    std::uint64_t id;
    // H10: incremented only after a complete accepted physical step while
    // the particle remains below K_in.  The in-class initializer keeps
    // legacy construction and migration buffers deterministic.
    std::uint32_t return_residence_steps = 0;
};

// CIC cell-deposit weights (identical to BeamPIC::add_shape_density /
// deposit_density conventions).  cell0 receives w0, cell1 receives w1;
// w0 + w1 == 1 and cell1 == cell0 + 1.
struct CellDepositWeights {
    int cell0;
    int cell1;
    double w0;
    double w1;
};

// Staggered face gather weights (identical to the BeamPIC shifted CIC
// gather).  Applied directly to Ex_face; never to a cell-averaged field.
struct FaceGatherWeights {
    int face0;
    int face1;
    double w0;
    double w1;
};

struct ParticleShape1D {
    static CellDepositWeights cell_weights(double x, const SpatialGrid& sg);
    static FaceGatherWeights shifted_face_weights(double x,
                                                  const SpatialGrid& sg);
};

// Face-resolved trajectory current accumulator (A/m^2).
struct FaceCurrentAccumulator {
    std::vector<double> current_face_x;
    void init(const SpatialGrid& sg)
    {
        current_face_x.assign(
            static_cast<size_t>(sg.nx_local) + 1, 0.0);
    }
};

// Charge-conserving trajectory deposition (section 6.3): the current is
// reconstructed from the CIC shape difference, satisfying
//   Delta(n dx)/dt + D_x(number flux) = 0
// per cell to summation error.  Multi-cell paths are split at every interior
// cell face; open-boundary outflow segments are truncated at the first
// physical-boundary crossing (constraint 4).
struct ChargeConservingTrajectory1D {
    // Core shape-difference rule: adds
    //   q_e * weight * (G_j(x1) - G_j(x0))
    // to charge_face[j], where G_j is the CIC cumulative weight strictly
    // below face j.  Only local faces [ix_start, ix_start+nx_local] are
    // written; rank-boundary face sums are closed in
    // BackgroundTailPIC::finalize_trajectory_current.
    static void deposit_charge_delta(double x0, double x1, double weight,
                                     const SpatialGrid& sg,
                                     std::vector<double>& charge_face);

    // Average current carried during the segment [x0,x1] of duration dt:
    //   J_face += q_e * weight * (G_face(x1) - G_face(x0)) / dt.
    static void deposit_segment(double x0, double x1, double weight,
                                double dt, const SpatialGrid& sg,
                                FaceCurrentAccumulator& current);
};

// Per-step open-boundary outflow accounting (section 6.2 / 2.3).
struct TailOutflowLedger {
    double left_number;
    double left_px;
    double left_kinetic_energy;
    double right_number;
    double right_px;
    double right_kinetic_energy;
    TailOutflowLedger()
        : left_number(0.0), left_px(0.0), left_kinetic_energy(0.0),
          right_number(0.0), right_px(0.0), right_kinetic_energy(0.0)
    {}
};

// Per-cell continuity and ledger closure result (section 6.6.2 / 6.6.3).
struct TailContinuityResult {
    double abs_l1;                  // sum |residual_i| * dx
    double abs_linf;
    double rel_l1;
    double rel_linf;
    double number_balance_error;    // |N1 - N0 - C + outflow| / scale
    double truncation_ledger_error; // |sum trunc_in_domain + trunc_ledger
                                    //  - outflow| / scale
    TailContinuityResult()
        : abs_l1(0.0), abs_linf(0.0), rel_l1(0.0), rel_linf(0.0),
          number_balance_error(0.0), truncation_ledger_error(0.0)
    {}
};

// Carrier used when a drift segment crosses an internal MPI rank boundary:
// the particle is placed at the shared edge and the remaining drift time is
// passed to the owning neighbor (section 6.3 constraint 5).
struct TailDriftCarrier {
    BackgroundTailParticle particle;
    double remaining_dt;
};

// Stage-H6 exact snapshot of the accepted tail state (section 12.1/12.3):
// particles, deposited density, ID counter, outflow/truncation/deposit
// ledgers, counter-based collision RNG key and moment stats.  export/
// import are bitwise round-trips used by the checkpoint writer/reader.
struct BackgroundTailStateSnapshot {
    std::vector<BackgroundTailParticle> particles;
    std::vector<double> density;
    std::uint64_t id_counter;
    TailOutflowLedger outflow;
    double truncation_shape_left;
    double truncation_shape_right;
    double deposit_shape_left;
    double deposit_shape_right;
    double deposit_shape_step_start_left;
    double deposit_shape_step_start_right;
    std::uint64_t collision_rng_seed;
    double max_abs_u;
    double max_kinetic_energy;
    BackgroundTailStateSnapshot()
        : id_counter(0), truncation_shape_left(0.0),
          truncation_shape_right(0.0), deposit_shape_left(0.0),
          deposit_shape_right(0.0), deposit_shape_step_start_left(0.0),
          deposit_shape_step_start_right(0.0), collision_rng_seed(0),
          max_abs_u(0.0), max_kinetic_energy(0.0)
    {}
};

// Background high-energy tail PIC state (section 6.2): local particles,
// cell-centered density, trajectory current, left/right outflow ledgers,
// boundary-shape ledgers, ID generation, transactional swap and diagnostics.
class BackgroundTailPIC {
public:
    std::vector<BackgroundTailParticle> particles;
    std::vector<double> density;           // local [nx_local] m^-3
    std::vector<double> current_face_x;    // local [nx_local+1] A/m^2
    std::vector<double> density_step_start; // accepted density at step start

    BackgroundTailPIC();

    void init(const SpatialGrid& sg);
    void clear();

    // Section 6.4: deterministic unique IDs.  High bits carry the creating
    // MPI rank, low bits a per-rank accepted-conversion counter.
    std::uint64_t next_particle_id(int mpi_rank);
    std::uint64_t id_counter() const { return id_counter_; }

    // Per-step bookkeeping: snapshot the accepted density and reset the step
    // ledgers, trajectory-charge accumulator and boundary-shape ledgers.
    void begin_step(const SpatialGrid& sg, double dt);

    // Relativistic D-K-D pieces (section 6.3):
    //   drift_half: x += (dt/2) c u_x/gamma, with open-boundary truncation,
    //               MPI migration and trajectory-current accumulation;
    //   kick:       u_x -= q_e dt E_x^{n+1/2}(x)/(m_e c); u_y/u_z unchanged.
    // A full step is drift_half(dt/2), kick(E^{n+1/2}, dt), drift_half(dt/2).
    void drift_half(const SpatialGrid& sg, double dt,
                    int mpi_rank, int mpi_size);
    // Gate C (section 7.4): optional rank-local kinetic-work output.  When
    // local_kinetic_work is non-NULL, the existing kick loop also accumulates
    //   sum_p weight_p * [K(u_p^after) - K(u_p^before)]
    // with the same relativistic kinetic energy used by
    // VpfpIntegrator::tail_total_kinetic_energy, using an OpenMP reduction.
    // A NULL pointer adds no per-particle sqrt, full-particle scan or MPI.
    // Gate I (section 4.4): when local_delta_ke_by_x is additionally
    // non-NULL, the same per-particle kick energy is distributed to the CIC
    // cells at the kick position x^{n+1/2} (nx_local entries, caller-owned).
    // Out-of-domain CIC shares are never renormalized back into the domain;
    // when local_delta_ke_boundary is non-NULL they are accumulated there so
    // that sum(cell) + boundary equals the same global work ledger (section
    // I3).
    void kick(const SpatialGrid& sg, const EMFields& field_mid, double dt,
              int mpi_rank, int mpi_size, double* local_kinetic_work = NULL,
              std::vector<double>* local_delta_ke_by_x = NULL,
              double* local_delta_ke_boundary = NULL);

    // CIC density deposit with rank-boundary guard exchange; physical-boundary
    // shape shares are accumulated into the deposit shape ledgers instead of
    // being renormalized into the domain (section 6.3 constraint 7).
    void deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size);

    // Turns the accumulated trajectory charge into the step-averaged face
    // current and closes the shared faces between MPI ranks.
    void finalize_trajectory_current(const SpatialGrid& sg, double dt,
                                     int mpi_rank, int mpi_size);

    // Section 6.6.5: production migration supports adjacent ranks only, so
    // c*dt must be strictly smaller than the local rank width.
    bool cfl_contract_ok(const SpatialGrid& sg, double dt,
                         int mpi_size) const;

    // Transactional accept (section 4.4): exchange the complete state.
    void swap_state(BackgroundTailPIC& other);

    // Stage H6: exact accepted-state serialization for the checkpoint
    // schema-v2 round trip (sections 12.1 and 12.3).  export_accepted_state
    // copies the full accepted representation; import_accepted_state
    // replaces this object (including the ID counter, outflow/truncation/
    // deposit ledgers and the counter-based collision RNG key) so a
    // restarted run continues with the same next particle IDs and ledgers.
    void export_accepted_state(BackgroundTailStateSnapshot& snapshot) const;
    void import_accepted_state(const BackgroundTailStateSnapshot& snapshot);

    bool finite() const;
    bool nonnegative_weights() const;
    bool migration_failed() const { return migration_failed_; }

    const TailOutflowLedger& outflow_ledger() const { return outflow_; }
    // Gate I read-only in-domain CIC shares removed when trajectories cross
    // an open physical boundary [m^-2 per local cell].
    const std::vector<double>& truncated_in_domain_number() const {
        return truncated_in_domain_;
    }
    // Outside-domain CIC shape shares (m^-2) accumulated per step:
    //   truncation ledgers: particles that exited a physical boundary;
    //   deposit ledgers:    in-domain particles whose shape straddles a
    //                       physical boundary (constraint 7).
    double truncation_shape_left() const { return truncation_shape_left_; }
    double truncation_shape_right() const { return truncation_shape_right_; }
    double deposit_shape_left() const { return deposit_shape_left_; }
    double deposit_shape_right() const { return deposit_shape_right_; }

    double max_abs_u() const { return max_abs_u_; }
    double max_kinetic_energy() const { return max_kinetic_energy_; }
    // Stage H6 (section 13.3): MPI-exchange wall time inside the last
    // drift_half migration pass (0 when no exchange happened).
    double last_migration_seconds() const { return last_migration_seconds_; }

    // Section 6.5 diagnostics.
    struct WeightStats {
        double min;
        double max;
        double mean;
        double std;
        double max_single_charge_fraction;
        size_t macro_particle_count;
        size_t per_cell_max_macro_count;
        double density_noise_estimate;
        WeightStats()
            : min(0.0), max(0.0), mean(0.0), std(0.0),
              max_single_charge_fraction(0.0), macro_particle_count(0),
              per_cell_max_macro_count(0), density_noise_estimate(0.0)
        {}
    };
    WeightStats weight_stats(
        const SpatialGrid& sg,
        const std::vector<double>* background_density = NULL) const;

    // Constraint 7 check: |sum_i n_i dx + deposit ledgers - total weight|.
    double density_weight_balance(const SpatialGrid& sg) const;

    // Section 6.6.2/6.6.3 audit.  Requires deposit_density() and
    // finalize_trajectory_current() to have been called for the final state.
    // conversion_source is an artificial conversion source density
    // [m^-3 s^-1] (empty vector = zero), positive into the tail.
    TailContinuityResult audit_continuity(
        const SpatialGrid& sg, double dt,
        const std::vector<double>& conversion_source,
        int mpi_rank, int mpi_size) const;

private:
    double step_dt_;
    TailOutflowLedger outflow_;
    double truncation_shape_left_;
    double truncation_shape_right_;
    double deposit_shape_left_;
    double deposit_shape_right_;
    // Deposit shape ledgers at the step start (captured by begin_step before
    // the per-step reset): needed by the number-balance audit when particles
    // straddle a physical boundary.
    double deposit_shape_step_start_left_;
    double deposit_shape_step_start_right_;
    double max_abs_u_;
    double max_kinetic_energy_;
    std::uint64_t id_counter_;
    // Reserved for the H8 counter-based collision RNG (section 6.4); H1 only
    // stores the seed and carries it through the transactional swap.
    std::uint64_t collision_rng_seed_;
    bool migration_failed_;
    double last_migration_seconds_;

    std::vector<double> trajectory_charge_face_x_; // C/m^2 accumulated
    std::vector<double> truncated_in_domain_;      // m^-2 per local cell

    std::vector<BackgroundTailParticle> keep_;
    std::vector<TailDriftCarrier> send_left_;
    std::vector<TailDriftCarrier> send_right_;
    std::vector<TailDriftCarrier> recv_left_;
    std::vector<TailDriftCarrier> recv_right_;

    // OpenMP per-thread working sets for the initial drift pass.
    std::vector<std::vector<BackgroundTailParticle> > thread_keep_;
    std::vector<std::vector<TailDriftCarrier> > thread_send_left_;
    std::vector<std::vector<TailDriftCarrier> > thread_send_right_;
    std::vector<std::vector<double> > thread_charge_;
    std::vector<std::vector<double> > thread_trunc_in_domain_;
    std::vector<double> thread_trunc_left_;
    std::vector<double> thread_trunc_right_;
    std::vector<TailOutflowLedger> thread_outflow_;

    struct DriftOutputs {
        std::vector<BackgroundTailParticle>* keep;
        std::vector<TailDriftCarrier>* send_left;
        std::vector<TailDriftCarrier>* send_right;
        std::vector<double>* charge_face;
        std::vector<double>* trunc_in_domain;
        double* trunc_left;
        double* trunc_right;
        TailOutflowLedger* outflow;
    };

    void process_drift_particle(const BackgroundTailParticle& p,
                                double remaining_dt, double vx,
                                double x_left, double x_right, double length,
                                const SpatialGrid& sg, int mpi_rank,
                                int mpi_size, DriftOutputs& out);
    void accumulate_segment(double x0, double x1, double weight,
                            double dt_seg, const SpatialGrid& sg,
                            std::vector<double>& charge_face);
    void add_removed_shape(const CellDepositWeights& cw, double weight,
                           const SpatialGrid& sg, DriftOutputs& out);
    void exchange_carriers(const SpatialGrid& sg, int mpi_rank,
                           int mpi_size);
    void refresh_moment_stats();
};

#endif
