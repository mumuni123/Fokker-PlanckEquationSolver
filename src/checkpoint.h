#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "beam_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include "vlasov_ampere_midpoint.h"

#include <string>

struct CheckpointControlState {
    long long step;
    double time_s;
    double dt_s;
    double next_snapshot_s;
    int last_snapshot_step;
    double cumulative_collision_energy_delta;
};

struct RankCheckpointHeader {
    char magic[16];
    unsigned int version;
    unsigned int endian_tag;
    int mpi_rank;
    int mpi_size;
    int ix_start;
    int nx_local;
    int nx_total;
    int nghost;
    int nv;
    int nmu;
    long long step;
    double time_s;
    double dt_s;
    unsigned long long f_count;
    unsigned long long ex_face_count;
    unsigned long long beam_particle_count;
};

struct CheckpointStateHashes {
    unsigned long long background;
    unsigned long long field_faces;
    unsigned long long beam;
};

// Filled only by the explicit baseline-core -> nested-tail conversion mode.
// Values are global after MPI reduction and are persisted with the converted
// checkpoint so a zero-tail A/B restart is never mistaken for production data.
struct CheckpointVelocityRemapAudit {
    bool applied;
    double mass_before;
    double mass_after;
    double parallel_momentum_before;
    double parallel_momentum_after;
    double kinetic_energy_before;
    double kinetic_energy_after;
    unsigned long long core_hash_before;
    unsigned long long core_hash_after;
};

unsigned long long checkpoint_hash64(const void* data, size_t bytes,
                                     unsigned long long seed = 1469598103934665603ULL);
unsigned long long checkpoint_configuration_hash();
unsigned long long checkpoint_configuration_hash(bool low_order_only,
                                                 bool high_order_enabled,
                                                 bool fct_enabled);
unsigned long long checkpoint_velocity_grid_hash(const Species& bkg);
unsigned long long checkpoint_physics_parameter_hash();
CheckpointStateHashes checkpoint_state_hashes(const Species& bkg, const BeamPIC& beam,
                                              const EMFields& fields,
                                              const SpatialGrid& sg,
                                              int mpi_rank, int mpi_size);
bool read_checkpoint_reference_hashes(const std::string& directory,
                                      CheckpointStateHashes& hashes,
                                      int mpi_rank, int mpi_size);
bool write_checkpoint(const std::string& directory, const CheckpointControlState& control,
                      const Species& bkg, const BeamPIC& beam, const EMFields& fields,
                      const SpatialGrid& sg, int mpi_rank, int mpi_size,
                      std::string& error, bool low_order_only = false,
                       bool high_order_enabled = true, bool fct_enabled = true,
                       const CheckpointVelocityRemapAudit* remap_audit = 0);
bool read_checkpoint(const std::string& directory, CheckpointControlState& control,
                     Species& bkg, BeamPIC& beam, EMFields& fields,
                     const SpatialGrid& sg, int mpi_rank, int mpi_size,
                     std::string& error, bool low_order_only = false,
                      bool high_order_enabled = true, bool fct_enabled = true,
                      bool allow_velocity_grid_remap = false,
                      CheckpointVelocityRemapAudit* remap_audit = 0);
bool write_checkpoint_velocity_remap_audit(const std::string& directory,
                                           const CheckpointVelocityRemapAudit& audit,
                                           int mpi_rank, int mpi_size,
                                           std::string& error);
bool write_midpoint_audit_state(const std::string& directory,
                                const VlasovAmpereMidpointSolver::MidpointAuditState& state,
                                const SpatialGrid& sg, int mpi_rank, int mpi_size,
                                std::string& error);
bool read_midpoint_audit_state(const std::string& directory,
                               VlasovAmpereMidpointSolver::MidpointAuditState& state,
                               const Species& species_template,
                               const EMFields& fields_template,
                               const SpatialGrid& sg, int mpi_rank, int mpi_size,
                               std::string& error);

#endif
