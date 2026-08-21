#ifndef VPFP_DIAGNOSTICS_H
#define VPFP_DIAGNOSTICS_H

#include "background_tail_pic.h"
#include "grid.h"
#include "vpfp_integrator.h"

#include <string>
#include <vector>

struct VpfpRunManifestConfig {
    std::string field_boundary;
    double left_electric_field;
    double phi_left;
    double phi_right;
    std::string background_boundary;
    std::string collision_model;
    std::string collision_interface_mode;
    std::string bulk_collision_integrator;
    bool collision_induced_conversion;
    std::string requested_tail_collision_kernel;
    std::string tail_tail_collision_backend;
    std::string tail_bulk_collision_backend;
    std::string tail_collision_weight_mode;
    int tail_collision_max_substeps;
    int tail_collision_max_particle_growth;
    bool pair_bulk_bulk;
    bool pair_bulk_tail;
    bool pair_tail_bulk;
    bool pair_tail_tail;
    bool population_control_enabled;
    int population_control_interval;
    bool tail_cell_moment_audit;
    bool tail_subcell_loading;
    int tail_cell_moment_audit_top_cells;
    std::string tail_conversion_mode;
    int tail_flux_quadrature_order;
    int tail_flux_max_supports;
    int tail_flux_max_created_particles_per_step;
    std::uint64_t physical_config_hash;
    std::uint64_t diagnostic_config_hash;
    std::uint64_t interface_topology_hash;
    std::uint64_t interface_topology_version;
    std::uint64_t interface_mask_hash;
    std::uint64_t interface_face_list_hash;
    // JC4 (section 7.2/7.5): field-particle coupling configuration persisted
    // in both the run manifest and checkpoint manifest.  Old checkpoints
    // missing these fields default to legacy.
    std::string coupling_mode;
    int coupling_max_iters;
    double coupling_relaxation;
    double coupling_field_tol;
    double coupling_pairing_tol;
    std::string background_phase_space_mode;
    std::string x_transport_velocity_mode;
    int x_transport_velocity_table_schema;
    VpfpRunManifestConfig()
        : left_electric_field(0.0), phi_left(0.0), phi_right(0.0),
          bulk_collision_integrator("bgk-validation"),
          collision_induced_conversion(true), tail_collision_max_substeps(1024),
          tail_collision_max_particle_growth(0), pair_bulk_bulk(false),
          pair_bulk_tail(false), pair_tail_bulk(false), pair_tail_tail(false),
          population_control_enabled(false), population_control_interval(0),
          tail_cell_moment_audit(false), tail_subcell_loading(false),
          tail_cell_moment_audit_top_cells(64),
          tail_conversion_mode("static-cell"), tail_flux_quadrature_order(4),
          tail_flux_max_supports(7),
          tail_flux_max_created_particles_per_step(0),
          physical_config_hash(0), diagnostic_config_hash(0),
          interface_topology_hash(0), interface_topology_version(0),
          interface_mask_hash(0), interface_face_list_hash(0),
          coupling_mode("legacy"), coupling_max_iters(12),
          coupling_relaxation(0.5), coupling_field_tol(1.0e-8),
          coupling_pairing_tol(1.0e-8),
          background_phase_space_mode("strang-ppm"),
          x_transport_velocity_mode("analytic-cell-center"),
          x_transport_velocity_table_schema(1)
    {}
};

class VpfpDiagnostics {
public:
    VpfpDiagnostics();
    // Returns false when the shared output root or its mandatory step ledger
    // cannot be created.  Production must not silently run without output.
    bool init(const std::string& output_dir, int rank, int level, int interval);
    void set_run_manifest_config(const VpfpRunManifestConfig& config)
    {
        run_config_ = config;
    }
    void write_accepted_step(int step, double time, const VpfpStepResult& result,
                             const OpenGaussDiagnostics& gauss,
                             const Species& electrons,
                             const BackgroundTailPIC* tail,
                             const SpatialGrid& grid,
                             int mpi_rank,
                             const FieldParticlePowerAuditWorkspace*
                                 pairing_ws = NULL);
    void write_failure(int step, double time, const VpfpStepResult& result,
                       int mpi_rank);
    // Stage-H5 population-control report (section 7.10): written only on
    // control steps from rank 0.
    void write_population_control(int step, double time,
                                  const VpfpStepResult& result,
                                  int mpi_rank);
    // Physical-time snapshot of the accepted state (sections 11.2/13.13):
    // per-rank spatial and velocity profiles plus a rank-0 manifest.  Never
    // gathers the full 4-D distribution; each rank writes its own files.
    // Stage-H6 (section 13.2): tail/partition are NULL when the solver runs
    // tail-off; convert_energy_mev anchors the threshold-interface report.
    void write_snapshot(int step, double time,
                        const Species& electrons, const BeamPIC& beam,
                        const EMFields& fields, const SpatialGrid& grid,
                        const BackgroundTailPIC* tail,
                        const HybridVelocityPartition* partition,
                        double convert_energy_mev,
                        int mpi_rank, int mpi_size);

private:
    std::string output_dir_;
    int level_;
    int interval_;
    // Threshold source data from the last accepted event-bearing step.  It
    // is globally reduced before rank 0 writes it; level 2 additionally
    // writes an accepted-step event ledger.
    std::vector<double> last_conversion_source_edges_;
    std::vector<double> last_conversion_source_spectrum_;
    std::uint64_t last_control_groups_;
    std::uint64_t last_control_fallbacks_;
    bool moment_audit_headers_written_;
    bool moment_audit_hist_headers_written_;
    bool pairing_headers_written_;
    bool pairing_breakdown_headers_written_;
    // JC4 (section 7.3): field-particle iteration accepted-step diagnostic.
    bool field_particle_iteration_headers_written_;
    VpfpRunManifestConfig run_config_;

    void write_fields(int step, double time, const Species& electrons,
                      const BeamPIC& beam, const EMFields& fields,
                      const SpatialGrid& grid, int mpi_rank);
    void write_moments_and_spectra(int step, double time,
                                   const Species& electrons,
                                   const BeamPIC& beam,
                                   const SpatialGrid& grid,
                                   const BackgroundTailPIC* tail,
                                   int mpi_rank);
    void write_tail_stats(int step, double time,
                          const Species& electrons,
                          const SpatialGrid& grid,
                          const BackgroundTailPIC* tail, int mpi_rank);
    void write_threshold_interface(int step, double time,
                                   const Species& electrons,
                                   const SpatialGrid& grid,
                                   const BackgroundTailPIC* tail,
                                   const HybridVelocityPartition* partition,
                                   double convert_energy_mev, int mpi_rank);
    void write_conversion_source_ledger(int step, double time,
                                         const VpfpStepResult& result,
                                         const BackgroundTailPIC* tail,
                                         int mpi_rank);
    void write_bulk_tail_moment_audit_accepted_step(
        int step, double time, const VpfpStepResult& result, int mpi_rank);
    void write_bulk_tail_flux_accepted_step(
        int step, double time, const VpfpStepResult& result, int mpi_rank);
    void write_stage_energy_audit_accepted_step(
        int step, double time, const VpfpStepResult& result, int mpi_rank);
    // Gate I (section 4.6): accepted-step-only field-particle power pairing
    // scalar file (rank 0) and per-rank profile (level 2, interval hit).
    void write_field_particle_power_pairing_accepted_step(
        int step, double time, const VpfpStepResult& result, int mpi_rank,
        const SpatialGrid& grid,
        const FieldParticlePowerAuditWorkspace* pairing_ws);
    // JC4 (section 7.3): field-particle iteration accepted-step diagnostic
    // file.  Written by rank 0 only; header validated on restart-append.
    void write_field_particle_iteration_accepted_step(
        int step, double time, const VpfpStepResult& result, int mpi_rank);
    void write_manifest(int step, double time, const Species& electrons,
                        const BeamPIC& beam, const SpatialGrid& grid,
                        const BackgroundTailPIC* tail,
                        const HybridVelocityPartition* partition,
                        double convert_energy_mev,
                        int mpi_rank, int mpi_size);
    void write_joint_midpoint_iterations(
        int step, double time, const VpfpStepResult& result, int mpi_rank);
};

#endif
