#ifndef VPFP_CHECKPOINT_H
#define VPFP_CHECKPOINT_H

#include "background_tail_pic.h"
#include "beam_pic.h"
#include "maxwell.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "species.h"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

struct VpfpCheckpointControl {
    int step;
    double time;
    double dt;
};

// JC4 (section 7.2/7.5): field-particle coupling configuration written to the
// checkpoint manifest.  Independent of the tail state; passed as a separate
// parameter to write_vpfp_checkpoint so non-tail runs also persist it.
struct VpfpCouplingManifestConfig {
    std::string mode;
    int max_iters;
    double relaxation;
    double field_tol;
    double pairing_tol;
    std::string background_phase_space_mode;
    std::string x_transport_velocity_mode;
    int x_transport_velocity_table_schema;
    VpfpCouplingManifestConfig()
        : mode("legacy"), max_iters(12), relaxation(0.5),
          field_tol(1.0e-8), pairing_tol(1.0e-8),
          background_phase_space_mode("strang-ppm"),
          x_transport_velocity_mode("analytic-cell-center"),
          x_transport_velocity_table_schema(1)
    {}
};

// Stage-H6 tail configuration and accepted cumulative ledgers stored in the
// checkpoint (section 12.1).  Conversion/control cumulative counters are
// per-rank local (they feed the per-rank restart state); the combined
// checksums are global values saved at checkpoint time.
struct VpfpCheckpointTailConfig {
    std::uint64_t partition_config_hash;
    // Keep the physical representation/configuration identity separate from
    // diagnostic switches.  A read-only audit must not change the physical
    // restart identity, while the diagnostic identity remains available for
    // provenance.
    std::uint64_t physical_config_hash;
    std::uint64_t diagnostic_config_hash;
    // Section 7.11.4 branch B: explicit threshold-aware energy edges and
    // their FNV-1a hash, written to the checkpoint manifest (not part of the
    // binary header; the binary validation is via partition_config_hash,
    // which folds the edges).
    std::vector<double> conversion_energy_edges;
    std::uint64_t conversion_energy_edges_hash;
    double convert_energy_mev;
    double buffer_width_mev;
    int upar_bins;
    int energy_bins;
    std::string return_mode;
    double return_energy_mev;
    int return_residence_steps;
    int return_max_stencil_radius;
    double return_moment_tolerance;
    std::string collision_kernel;
    std::string collision_weight_mode;
    std::string collision_interface_mode;
    std::string bulk_collision_integrator;
    bool collision_induced_conversion;
    int collision_max_substeps;
    int collision_max_particle_growth;
    bool population_control_enabled;
    int control_interval;
    int target_particles_per_phase_bin;
    int max_particles_per_phase_bin;
    double max_weight_ratio;
    int max_support;
    std::string conversion_mode;
    int flux_quadrature_order;
    int flux_max_supports;
    int flux_max_created_particles_per_step;
    std::uint64_t interface_topology_hash;
    std::uint64_t interface_topology_version;
    std::uint64_t interface_mask_hash;
    std::uint64_t interface_face_list_hash;
    bool interface_topology_metadata_present;
    bool interface_topology_detail_metadata_present;
    bool conversion_metadata_present;
    double conversion_cumulative_number;
    double conversion_cumulative_px;
    double conversion_cumulative_energy;
    std::uint64_t conversion_cumulative_particles_created;
    std::uint64_t tail_cumulative_outflow_number;
    std::uint64_t control_cumulative_groups;
    std::uint64_t control_cumulative_fallbacks;
    double return_cumulative_number;
    double return_cumulative_px;
    double return_cumulative_jx_dx;
    double return_cumulative_energy;
    double return_cumulative_pixx_dx;
    double return_cumulative_piperp_dx;
    std::uint64_t return_cumulative_particles_removed;
    std::uint64_t return_cumulative_deferred_groups;
    double combined_number;
    double combined_kinetic_energy;
    double combined_field_energy;
    // JC4 (section 7.2/7.5): field-particle coupling configuration persisted
    // in the checkpoint manifest.  Old checkpoints missing these fields
    // default to legacy.
    std::string coupling_mode;
    int coupling_max_iters;
    double coupling_relaxation;
    double coupling_field_tol;
    double coupling_pairing_tol;
    VpfpCheckpointTailConfig()
        : partition_config_hash(0), physical_config_hash(0),
          diagnostic_config_hash(0), conversion_energy_edges_hash(0),
          convert_energy_mev(0.0), buffer_width_mev(0.0), upar_bins(0),
          energy_bins(0), return_mode("none"), return_energy_mev(0.0),
          return_residence_steps(0), return_max_stencil_radius(0),
          return_moment_tolerance(0.0),
          collision_kernel("none"),
          collision_weight_mode("equal-strata"),
          collision_interface_mode("none"),
          bulk_collision_integrator("bgk-validation"),
          collision_induced_conversion(true),
          collision_max_substeps(1024),
          collision_max_particle_growth(0), population_control_enabled(false),
          control_interval(0), target_particles_per_phase_bin(64),
          max_particles_per_phase_bin(1024), max_weight_ratio(8.0),
          max_support(7), conversion_mode("static-cell"),
          flux_quadrature_order(4), flux_max_supports(7),
          flux_max_created_particles_per_step(0),
          interface_topology_hash(0),
          interface_topology_version(0), interface_mask_hash(0),
          interface_face_list_hash(0),
          interface_topology_metadata_present(false),
          interface_topology_detail_metadata_present(false),
          conversion_metadata_present(false),
          conversion_cumulative_number(0.0),
          conversion_cumulative_px(0.0), conversion_cumulative_energy(0.0),
          conversion_cumulative_particles_created(0),
          tail_cumulative_outflow_number(0),
          control_cumulative_groups(0), control_cumulative_fallbacks(0),
          return_cumulative_number(0.0), return_cumulative_px(0.0),
          return_cumulative_jx_dx(0.0), return_cumulative_energy(0.0),
          return_cumulative_pixx_dx(0.0), return_cumulative_piperp_dx(0.0),
          return_cumulative_particles_removed(0),
          return_cumulative_deferred_groups(0),
          combined_number(0.0), combined_kinetic_energy(0.0),
          combined_field_energy(0.0),
          coupling_mode("legacy"), coupling_max_iters(12),
          coupling_relaxation(0.5), coupling_field_tol(1.0e-8),
          coupling_pairing_tol(1.0e-8)
    {}
};

// Stage-H6 accepted tail state carried by the checkpoint (section 12.1):
// the full BackgroundTailPIC accepted representation plus the
// configuration/ledger block.  present=false means the checkpoint (or the
// running solver) has no tail.
struct VpfpCheckpointTailState {
    bool present;
    VpfpCheckpointTailConfig config;
    BackgroundTailPIC tail;
    VpfpCheckpointTailState() : present(false) {}
};

// Narrow serialization seam used by checkpoint compatibility tests.  The
// production checkpoint reader calls the same decoder; versions 1/2 reset
// the H10 residence counter because those schemas did not store it.
bool read_vpfp_tail_particle_records(
    std::istream& input, int checkpoint_version, size_t count,
    std::vector<BackgroundTailParticle>& particles);

// Shared restart gate used by the production CLI and the 17D checkpoint
// tests.  A mismatch is a configuration error, not a reason to reinterpret
// the stored particle state under a different flux interface.
bool validate_vpfp_checkpoint_tail_config(
    const VpfpCheckpointTailConfig& stored,
    const VpfpCheckpointTailConfig& expected,
    std::string& error);

// The 17E A/B experiment may change the conversion implementation only at a
// checkpoint that has never transferred any bulk support to the PIC tail.
// This does not reinterpret an existing tail state.  All other mode changes
// remain incompatible restart requests.
bool vpfp_preconversion_static_to_flux_restart_allowed(
    const VpfpCheckpointTailConfig& stored,
    const std::string& requested_conversion_mode,
    std::uint64_t global_tail_particle_count);

// JC4 (section 7.5): read field-particle coupling configuration from a
// checkpoint manifest.  Returns true if the file was readable; the five
// output values default to legacy when the manifest is absent or lacks the
// coupling keys (old checkpoint).
bool read_coupling_config_from_manifest(
    const std::string& directory,
    std::string& coupling_mode, int& coupling_max_iters,
    double& coupling_relaxation, double& coupling_field_tol,
    double& coupling_pairing_tol,
    std::string& x_transport_velocity_mode,
    int& x_transport_velocity_table_schema);

bool read_background_phase_space_mode_from_manifest(
    const std::string& directory, std::string& mode);

// Schema-v4 checkpoint write (section 12/H10).  tail_state must be NULL when
// the solver runs with --background-tail-mode off and non-NULL (with
// present=true) otherwise.  Each rank writes its own binary file plus a
// rank-0 manifest; the manifest carries the section 12.1 configuration
// keys.
bool write_vpfp_checkpoint(
    const std::string& directory, const VpfpCheckpointControl& control,
    const Species& electrons, const BeamPIC& beam, const EMFields& fields,
    const SpatialGrid& grid, const ElectrostaticBoundary& field_boundary,
    const OpenBackgroundBoundaryConfig& background_boundary,
    const std::string& collision_model,
    const VpfpCheckpointTailState* tail_state,
    const VpfpCouplingManifestConfig& coupling_config,
    int mpi_rank, int mpi_size, std::string& error);

// Schema-v4 checkpoint read.  tail_state must be NULL for a tail-disabled
// solver; otherwise the restored tail state is returned (validated against
// the current grid).  Section 12.2: an old no-tail checkpoint is refused
// when the solver expects a tail, and a checkpoint containing a tail is
// refused when the solver runs tail-off; no silent physical modification is
// performed.
bool read_vpfp_checkpoint(
    const std::string& directory, VpfpCheckpointControl& control,
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const SpatialGrid& grid, VpfpCheckpointTailState* tail_state,
    int mpi_rank, int mpi_size, std::string& error);

#endif
