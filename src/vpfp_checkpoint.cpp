#include "vpfp_checkpoint.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#ifdef _WIN32
#include <direct.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

namespace {
struct RankHeader {
    unsigned long long magic;
    int version;
    int nx_global;
    int nx_local;
    int nv;
    int nmu;
    int rank;
    int mpi_size;
    VpfpCheckpointControl control;
    BeamPersistentState beam_state;
    unsigned long long f_size;
    unsigned long long particle_count;
    unsigned long long ex_face_size;
    unsigned long long ex_size;
    unsigned long long phi_size;
    // --- schema-v2 tail block (section 12.1) ---
    int tail_present;
    unsigned long long tail_particle_count;
    unsigned long long tail_density_size;
    unsigned long long tail_id_counter;
    unsigned long long tail_collision_rng_seed;
    double tail_outflow_left_number;
    double tail_outflow_left_px;
    double tail_outflow_left_ke;
    double tail_outflow_right_number;
    double tail_outflow_right_px;
    double tail_outflow_right_ke;
    double tail_truncation_shape_left;
    double tail_truncation_shape_right;
    double tail_deposit_shape_left;
    double tail_deposit_shape_right;
    double tail_deposit_shape_step_start_left;
    double tail_deposit_shape_step_start_right;
    double tail_max_abs_u;
    double tail_max_kinetic_energy;
    unsigned long long partition_config_hash;
    double convert_energy_mev;
    double buffer_width_mev;
    int upar_bins;
    int energy_bins;
    int population_control_enabled;
    int control_interval;
    int target_particles_per_phase_bin;
    int max_particles_per_phase_bin;
    double max_weight_ratio;
    int max_support;
    char return_mode[16];
    double return_energy_mev;
    int return_residence_steps;
    int return_max_stencil_radius;
    double return_moment_tolerance;
    char collision_kernel[32];
    char collision_weight_mode[32];
    int collision_max_substeps;
    int collision_max_particle_growth;
    double conversion_cumulative_number;
    double conversion_cumulative_px;
    double conversion_cumulative_energy;
    unsigned long long conversion_cumulative_particles_created;
    unsigned long long tail_cumulative_outflow_number;
    unsigned long long control_cumulative_groups;
    unsigned long long control_cumulative_fallbacks;
    double return_cumulative_number;
    double return_cumulative_px;
    double return_cumulative_jx_dx;
    double return_cumulative_energy;
    double return_cumulative_pixx_dx;
    double return_cumulative_piperp_dx;
    unsigned long long return_cumulative_particles_removed;
    unsigned long long return_cumulative_deferred_groups;
    double combined_number;
    double combined_kinetic_energy;
    double combined_field_energy;
};

// Exact schema-v2 header layout, retained solely to load the H1--H9 tail
// checkpoints.  Do not use it for new writes.
struct RankHeaderV2 {
    unsigned long long magic; int version; int nx_global; int nx_local;
    int nv; int nmu; int rank; int mpi_size;
    VpfpCheckpointControl control; BeamPersistentState beam_state;
    unsigned long long f_size; unsigned long long particle_count;
    unsigned long long ex_face_size; unsigned long long ex_size;
    unsigned long long phi_size; int tail_present;
    unsigned long long tail_particle_count; unsigned long long tail_density_size;
    unsigned long long tail_id_counter; unsigned long long tail_collision_rng_seed;
    double tail_outflow_left_number; double tail_outflow_left_px;
    double tail_outflow_left_ke; double tail_outflow_right_number;
    double tail_outflow_right_px; double tail_outflow_right_ke;
    double tail_truncation_shape_left; double tail_truncation_shape_right;
    double tail_deposit_shape_left; double tail_deposit_shape_right;
    double tail_deposit_shape_step_start_left;
    double tail_deposit_shape_step_start_right; double tail_max_abs_u;
    double tail_max_kinetic_energy; unsigned long long partition_config_hash;
    double convert_energy_mev; double buffer_width_mev; int upar_bins;
    int energy_bins; int population_control_enabled; int control_interval;
    int target_particles_per_phase_bin; int max_particles_per_phase_bin;
    double max_weight_ratio; int max_support; char return_mode[16];
    char collision_kernel[32]; char collision_weight_mode[32];
    int collision_max_substeps; int collision_max_particle_growth;
    double conversion_cumulative_number; double conversion_cumulative_px;
    double conversion_cumulative_energy;
    unsigned long long conversion_cumulative_particles_created;
    unsigned long long tail_cumulative_outflow_number;
    unsigned long long control_cumulative_groups;
    unsigned long long control_cumulative_fallbacks; double combined_number;
    double combined_kinetic_energy; double combined_field_energy;
};

bool read_rank_header(std::istream& in, RankHeader& header)
{
    unsigned long long magic = 0;
    int version = 0;
    if (!in.read(reinterpret_cast<char*>(&magic), sizeof(magic)) ||
        !in.read(reinterpret_cast<char*>(&version), sizeof(version))) return false;
    in.clear(); in.seekg(0, std::ios::beg);
    if (version >= 3) {
        return static_cast<bool>(in.read(reinterpret_cast<char*>(&header),
                                         sizeof(header)));
    }
    if (version != 1 && version != 2) return false;
    RankHeaderV2 old = {};
    if (!in.read(reinterpret_cast<char*>(&old), sizeof(old))) return false;
    std::memset(&header, 0, sizeof(header));
    std::memcpy(&header, &old, offsetof(RankHeader, return_energy_mev));
    std::memcpy(&header.collision_kernel, &old.collision_kernel,
                sizeof(RankHeaderV2) - offsetof(RankHeaderV2, collision_kernel));
    header.return_energy_mev = 0.0;
    header.return_residence_steps = 0;
    header.return_max_stencil_radius = 0;
    header.return_moment_tolerance = 0.0;
    return magic != 0;
}

struct RankPayload {
    RankHeader header;
    std::vector<double> f;
    std::vector<BeamParticle> beam_particles;
    std::vector<BackgroundTailParticle> tail_particles;
    std::vector<double> tail_density;
    std::vector<double> ex_face;
    std::vector<double> ex;
    std::vector<double> phi;
};

// Schema-v4 serializes each v2 particle field independently.  Do not replace
// this with a record struct: compiler alignment after the uint32_t residence
// field would silently make the on-disk layout ABI-dependent again.
// Retain the short-lived v3 padded record only to read checkpoints written
// before the explicit-field correction.
struct TailParticleRecordV3Padded {
    double x;
    double ux;
    double uy;
    double uz;
    double weight;
    std::uint64_t id;
    std::uint32_t return_residence_steps;
};
struct TailParticleRecordV1 {
    double x;
    double ux;
    double uy;
    double uz;
    double weight;
    std::uint64_t id;
};
const unsigned long long checkpoint_magic = 0x565046504F50454EULL;
const int checkpoint_version = 4;

bool write_tail_particles(std::ostream& out,
                          const std::vector<BackgroundTailParticle>& particles)
{
    for (size_t i = 0; i < particles.size(); ++i) {
        const BackgroundTailParticle& p = particles[i];
        out.write(reinterpret_cast<const char*>(&p.x), sizeof(p.x));
        out.write(reinterpret_cast<const char*>(&p.ux), sizeof(p.ux));
        out.write(reinterpret_cast<const char*>(&p.uy), sizeof(p.uy));
        out.write(reinterpret_cast<const char*>(&p.uz), sizeof(p.uz));
        out.write(reinterpret_cast<const char*>(&p.weight), sizeof(p.weight));
        out.write(reinterpret_cast<const char*>(&p.id), sizeof(p.id));
        out.write(reinterpret_cast<const char*>(&p.return_residence_steps),
                  sizeof(p.return_residence_steps));
        if (!out) return false;
    }
    return true;
}

bool read_tail_particles(std::istream& in, int version, size_t count,
                         std::vector<BackgroundTailParticle>& particles)
{
    particles.resize(count);
    if (version >= 4) {
        for (size_t i = 0; i < count; ++i) {
            BackgroundTailParticle& p = particles[i];
            if (!in.read(reinterpret_cast<char*>(&p.x), sizeof(p.x)) ||
                !in.read(reinterpret_cast<char*>(&p.ux), sizeof(p.ux)) ||
                !in.read(reinterpret_cast<char*>(&p.uy), sizeof(p.uy)) ||
                !in.read(reinterpret_cast<char*>(&p.uz), sizeof(p.uz)) ||
                !in.read(reinterpret_cast<char*>(&p.weight), sizeof(p.weight)) ||
                !in.read(reinterpret_cast<char*>(&p.id), sizeof(p.id)) ||
                !in.read(reinterpret_cast<char*>(&p.return_residence_steps),
                         sizeof(p.return_residence_steps))) return false;
        }
    } else if (version == 3) {
        for (size_t i = 0; i < count; ++i) {
            TailParticleRecordV3Padded record;
            if (!in.read(reinterpret_cast<char*>(&record), sizeof(record)))
                return false;
            BackgroundTailParticle& p = particles[i];
            p.x = record.x; p.ux = record.ux; p.uy = record.uy;
            p.uz = record.uz; p.weight = record.weight; p.id = record.id;
            p.return_residence_steps = record.return_residence_steps;
        }
    } else {
        // v1/v2 used the old ABI-dependent structure.  Existing tail
        // checkpoints have no return history, so reset it by definition.
        for (size_t i = 0; i < count; ++i) {
            TailParticleRecordV1 record;
            if (!in.read(reinterpret_cast<char*>(&record), sizeof(record)))
                return false;
            BackgroundTailParticle& p = particles[i];
            p.x = record.x; p.ux = record.ux; p.uy = record.uy;
            p.uz = record.uz; p.weight = record.weight; p.id = record.id;
            p.return_residence_steps = 0;
        }
    }
    return true;
}

std::string rank_file(const std::string& directory, int rank)
{
    char buffer[64];
    std::sprintf(buffer, "/rank_%06d.bin", rank);
    return directory + buffer;
}

bool directory_exists(const std::string& path)
{
#ifdef _WIN32
    struct _stat status;
    return _stat(path.c_str(), &status) == 0 &&
           (status.st_mode & _S_IFDIR) != 0;
#else
    struct stat status;
    return stat(path.c_str(), &status) == 0 && S_ISDIR(status.st_mode);
#endif
}

bool create_one_directory(const std::string& path)
{
    if (path.empty() || directory_exists(path)) return true;
#ifdef _WIN32
    const int rc = _mkdir(path.c_str());
#else
    const int rc = mkdir(path.c_str(), 0755);
#endif
    return rc == 0 || (errno == EEXIST && directory_exists(path));
}

// C++11 has no std::filesystem.  Checkpoint paths are commonly nested below
// an output directory, so creating only the leaf silently fails when an
// intermediate component is absent.
bool create_directories(const std::string& path)
{
    if (path.empty()) return false;
    std::string normalized = path;
    for (size_t i = 0; i < normalized.size(); ++i) {
        if (normalized[i] == '\\') normalized[i] = '/';
    }

    for (size_t i = 1; i <= normalized.size(); ++i) {
        if (i != normalized.size() && normalized[i] != '/') continue;
        std::string component = normalized.substr(0, i);
        while (component.size() > 1 && component[component.size() - 1] == '/') {
            component.erase(component.size() - 1);
        }
#ifdef _WIN32
        if (component.size() == 2 && component[1] == ':') continue;
#endif
        if (!component.empty() && !create_one_directory(component)) return false;
    }
    return directory_exists(normalized);
}

void copy_string_into(const std::string& value, char* target, size_t size)
{
    const size_t n = std::min(value.size(), size - 1);
    std::memcpy(target, value.data(), n);
    target[n] = '\0';
}

std::string copy_string_out(const char* source, size_t size)
{
    char buffer[128];
    const size_t n = std::min(size, sizeof(buffer) - 1);
    std::memcpy(buffer, source, n);
    buffer[n] = '\0';
    return std::string(buffer);
}

bool read_tail_conversion_manifest(const std::string& directory,
                                   VpfpCheckpointTailConfig& config);

bool read_rank_payload(const std::string& path, RankPayload& payload,
                       std::string& error)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in || !read_rank_header(in, payload.header)) {
        error = "cannot read VPFP rank payload: " + path;
        return false;
    }
    const RankHeader& h = payload.header;
    if (h.magic != checkpoint_magic ||
        (h.version != 1 && h.version != 2 && h.version != 3 &&
         h.version != checkpoint_version) ||
        h.nx_global <= 0 || h.nx_local <= 0 || h.nv != Param::Nv ||
        h.nmu != Param::Nmu || h.f_size == 0 || h.ex_face_size == 0 ||
        h.ex_size == 0 || h.phi_size == 0) {
        error = "invalid VPFP rank payload header: " + path;
        return false;
    }
    payload.f.resize(static_cast<size_t>(h.f_size));
    payload.beam_particles.resize(static_cast<size_t>(h.particle_count));
    payload.ex_face.resize(static_cast<size_t>(h.ex_face_size));
    payload.ex.resize(static_cast<size_t>(h.ex_size));
    payload.phi.resize(static_cast<size_t>(h.phi_size));
    if (h.tail_present) {
        payload.tail_particles.resize(
            static_cast<size_t>(h.tail_particle_count));
        payload.tail_density.resize(
            static_cast<size_t>(h.tail_density_size));
    }
    in.read(reinterpret_cast<char*>(payload.f.data()),
            static_cast<std::streamsize>(payload.f.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(payload.beam_particles.data()),
            static_cast<std::streamsize>(
                payload.beam_particles.size() * sizeof(BeamParticle)));
    if (h.tail_present) {
        if (!read_tail_particles(in, h.version, payload.tail_particles.size(),
                                 payload.tail_particles)) {
            error = "truncated VPFP tail particle payload: " + path;
            return false;
        }
        in.read(reinterpret_cast<char*>(payload.tail_density.data()),
                static_cast<std::streamsize>(payload.tail_density.size() *
                                              sizeof(double)));
    }
    in.read(reinterpret_cast<char*>(payload.ex_face.data()),
            static_cast<std::streamsize>(payload.ex_face.size() *
                                         sizeof(double)));
    in.read(reinterpret_cast<char*>(payload.ex.data()),
            static_cast<std::streamsize>(payload.ex.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(payload.phi.data()),
            static_cast<std::streamsize>(payload.phi.size() * sizeof(double)));
    if (!in) {
        error = "truncated VPFP rank payload: " + path;
        return false;
    }
    return true;
}

bool read_repartitioned_checkpoint(
    const std::string& directory, VpfpCheckpointControl& control,
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const SpatialGrid& grid, VpfpCheckpointTailState* tail_state,
    int mpi_rank, int mpi_size, std::string& error)
{
    (void)mpi_size;
    RankPayload first;
    if (!read_rank_payload(rank_file(directory, 0), first, error)) return false;
    const int old_size = first.header.mpi_size;
    if (old_size <= 0 || first.header.nx_global != grid.nx_global) {
        error = "checkpoint global grid mismatch during repartition";
        return false;
    }
    if (tail_state != NULL && !first.header.tail_present) {
        error = "repartitioned checkpoint has no tail state";
        return false;
    }
    if (tail_state == NULL && first.header.tail_present) {
        error = "repartitioned checkpoint contains an unexpected tail";
        return false;
    }
    std::vector<RankPayload> payloads(static_cast<size_t>(old_size));
    payloads[0] = first;
    for (int old_rank = 1; old_rank < old_size; ++old_rank) {
        if (!read_rank_payload(rank_file(directory, old_rank),
                               payloads[static_cast<size_t>(old_rank)],
                               error)) return false;
    }
    const int nxg = grid.nx_global;
    const int nvmu = Param::Nvmu;
    const int ng = grid.nghost;
    const size_t global_f_size = static_cast<size_t>(nxg) *
                                 static_cast<size_t>(nvmu);
    std::vector<double> global_f(global_f_size, 0.0);
    std::vector<double> global_ex(static_cast<size_t>(nxg), 0.0);
    std::vector<double> global_phi(static_cast<size_t>(nxg), 0.0);
    std::vector<double> global_ex_face(static_cast<size_t>(nxg + 1), 0.0);
    std::vector<unsigned char> cell_seen(static_cast<size_t>(nxg), 0);
    std::vector<unsigned char> face_seen(static_cast<size_t>(nxg + 1), 0);
    std::vector<BeamParticle> all_beam;
    std::vector<BackgroundTailParticle> all_tail;
    std::vector<double> global_tail_density(static_cast<size_t>(nxg), 0.0);
    const double length = grid.length();
    for (int old_rank = 0; old_rank < old_size; ++old_rank) {
        const RankPayload& p = payloads[static_cast<size_t>(old_rank)];
        const RankHeader& h = p.header;
        if (h.mpi_size != old_size || h.rank != old_rank ||
            h.nx_global != nxg) {
            error = "inconsistent rank headers during repartition";
            return false;
        }
        SpatialGrid old_grid;
        old_grid.init_with_domain(old_rank, old_size, nxg, length);
        if (h.nx_local != old_grid.nx_local ||
            p.f.size() != static_cast<size_t>(old_grid.nx_total) *
                              static_cast<size_t>(nvmu)) {
            error = "old rank local array contract mismatch";
            return false;
        }
        for (int il = 0; il < old_grid.nx_local; ++il) {
            const int ig = old_grid.ix_start + il;
            if (cell_seen[static_cast<size_t>(ig)] != 0) {
                error = "duplicate physical cell in checkpoint";
                return false;
            }
            cell_seen[static_cast<size_t>(ig)] = 1;
            for (int slot = 0; slot < nvmu; ++slot) {
                global_f[static_cast<size_t>(ig) * nvmu + slot] =
                    p.f[static_cast<size_t>(old_grid.nghost + il) * nvmu + slot];
            }
            global_ex[static_cast<size_t>(ig)] =
                p.ex[static_cast<size_t>(old_grid.nghost + il)];
            global_phi[static_cast<size_t>(ig)] =
                p.phi[static_cast<size_t>(old_grid.nghost + il)];
            if (h.tail_present && p.tail_density.size() ==
                    static_cast<size_t>(old_grid.nx_local)) {
                global_tail_density[static_cast<size_t>(ig)] =
                    p.tail_density[static_cast<size_t>(il)];
            }
        }
        for (int lf = 0; lf <= old_grid.nx_local; ++lf) {
            const int gf = old_grid.ix_start + lf;
            if (gf < 0 || gf > nxg ||
                static_cast<size_t>(lf) >= p.ex_face.size()) {
                error = "invalid physical face in checkpoint";
                return false;
            }
            const double value = p.ex_face[static_cast<size_t>(lf)];
            if (face_seen[static_cast<size_t>(gf)] != 0 &&
                global_ex_face[static_cast<size_t>(gf)] != value) {
                error = "shared face values differ across checkpoint ranks";
                return false;
            }
            face_seen[static_cast<size_t>(gf)] = 1;
            global_ex_face[static_cast<size_t>(gf)] = value;
        }
        all_beam.insert(all_beam.end(), p.beam_particles.begin(),
                        p.beam_particles.end());
        if (h.tail_present) {
            all_tail.insert(all_tail.end(), p.tail_particles.begin(),
                            p.tail_particles.end());
        }
    }
    for (int ix = 0; ix < nxg; ++ix) {
        if (cell_seen[static_cast<size_t>(ix)] == 0 ||
            face_seen[static_cast<size_t>(ix)] == 0) {
            error = "checkpoint has a physical cell/face hole";
            return false;
        }
    }
    if (face_seen[static_cast<size_t>(nxg)] == 0) {
        error = "checkpoint has a final face hole";
        return false;
    }
    electrons.f.assign(electrons.f.size(), 0.0);
    fields.Ex.assign(fields.Ex.size(), 0.0);
    fields.phi.assign(fields.phi.size(), 0.0);
    fields.Ex_face.assign(fields.Ex_face.size(), 0.0);
    for (int il = 0; il < grid.nx_local; ++il) {
        const int ig = grid.ix_start + il;
        for (int slot = 0; slot < nvmu; ++slot) {
            electrons.f[static_cast<size_t>(ng + il) * nvmu + slot] =
                global_f[static_cast<size_t>(ig) * nvmu + slot];
        }
        fields.Ex[static_cast<size_t>(ng + il)] = global_ex[static_cast<size_t>(ig)];
        fields.phi[static_cast<size_t>(ng + il)] = global_phi[static_cast<size_t>(ig)];
    }
    for (int g = 0; g < grid.nx_total; ++g) {
        const int il = std::max(0, std::min(grid.nx_local - 1,
                                             g - ng));
        fields.Ex[static_cast<size_t>(g)] =
            fields.Ex[static_cast<size_t>(ng + il)];
        fields.phi[static_cast<size_t>(g)] =
            fields.phi[static_cast<size_t>(ng + il)];
    }
    for (int lf = 0; lf <= grid.nx_local; ++lf) {
        fields.Ex_face[static_cast<size_t>(lf)] =
            global_ex_face[static_cast<size_t>(grid.ix_start + lf)];
    }
    const int current_first = grid.ix_start;
    const int current_last = grid.ix_start + grid.nx_local;
    beam.particles.clear();
    for (size_t i = 0; i < all_beam.size(); ++i) {
        const double x = all_beam[i].x;
        int ig = static_cast<int>(std::floor(x / grid.dx));
        if (ig < 0) ig = 0;
        if (ig >= nxg) ig = nxg - 1;
        if (ig >= current_first && ig < current_last)
            beam.particles.push_back(all_beam[i]);
    }
    beam.import_persistent_state(payloads[0].header.beam_state, grid);
    if (tail_state != NULL) {
        BackgroundTailStateSnapshot snapshot;
        snapshot.density.assign(global_tail_density.begin() + grid.ix_start,
                                global_tail_density.begin() +
                                    grid.ix_start + grid.nx_local);
        std::uint64_t next_local_id = 0;
        for (size_t i = 0; i < all_tail.size(); ++i) {
            const double x = all_tail[i].x;
            int ig = static_cast<int>(std::floor(x / grid.dx));
            if (ig < 0) ig = 0;
            if (ig >= nxg) ig = nxg - 1;
            if (ig >= current_first && ig < current_last) {
                snapshot.particles.push_back(all_tail[i]);
                if ((all_tail[i].id >> 32) ==
                    static_cast<std::uint64_t>(static_cast<unsigned int>(mpi_rank))) {
                    const std::uint64_t candidate_id =
                        static_cast<std::uint64_t>(all_tail[i].id &
                                                   static_cast<std::uint64_t>(0xffffffffULL)) +
                        static_cast<std::uint64_t>(1);
                    next_local_id = std::max(
                        next_local_id, candidate_id);
                }
            }
        }
        snapshot.id_counter = next_local_id;
        snapshot.collision_rng_seed =
            payloads[0].header.tail_collision_rng_seed;
        snapshot.outflow.left_number =
            payloads[0].header.tail_outflow_left_number;
        snapshot.outflow.left_px = payloads[0].header.tail_outflow_left_px;
        snapshot.outflow.left_kinetic_energy =
            payloads[0].header.tail_outflow_left_ke;
        snapshot.outflow.right_number =
            payloads[0].header.tail_outflow_right_number;
        snapshot.outflow.right_px = payloads[0].header.tail_outflow_right_px;
        snapshot.outflow.right_kinetic_energy =
            payloads[0].header.tail_outflow_right_ke;
        snapshot.truncation_shape_left =
            payloads[0].header.tail_truncation_shape_left;
        snapshot.truncation_shape_right =
            payloads[0].header.tail_truncation_shape_right;
        snapshot.deposit_shape_left =
            payloads[0].header.tail_deposit_shape_left;
        snapshot.deposit_shape_right =
            payloads[0].header.tail_deposit_shape_right;
        snapshot.deposit_shape_step_start_left =
            payloads[0].header.tail_deposit_shape_step_start_left;
        snapshot.deposit_shape_step_start_right =
            payloads[0].header.tail_deposit_shape_step_start_right;
        snapshot.max_abs_u = payloads[0].header.tail_max_abs_u;
        snapshot.max_kinetic_energy = payloads[0].header.tail_max_kinetic_energy;
        tail_state->present = true;
        tail_state->tail.init(grid);
        tail_state->tail.import_accepted_state(snapshot);
        VpfpCheckpointTailConfig& c = tail_state->config;
        c.partition_config_hash = payloads[0].header.partition_config_hash;
        c.convert_energy_mev = payloads[0].header.convert_energy_mev;
        c.buffer_width_mev = payloads[0].header.buffer_width_mev;
        c.upar_bins = payloads[0].header.upar_bins;
        c.energy_bins = payloads[0].header.energy_bins;
        c.population_control_enabled =
            payloads[0].header.population_control_enabled != 0;
        c.control_interval = payloads[0].header.control_interval;
        c.target_particles_per_phase_bin =
            payloads[0].header.target_particles_per_phase_bin;
        c.max_particles_per_phase_bin =
            payloads[0].header.max_particles_per_phase_bin;
        c.max_weight_ratio = payloads[0].header.max_weight_ratio;
        c.max_support = payloads[0].header.max_support;
        c.return_mode = copy_string_out(payloads[0].header.return_mode,
                                        sizeof(payloads[0].header.return_mode));
        c.return_energy_mev = payloads[0].header.return_energy_mev;
        c.return_residence_steps =
            payloads[0].header.return_residence_steps;
        c.return_max_stencil_radius =
            payloads[0].header.return_max_stencil_radius;
        c.return_moment_tolerance =
            payloads[0].header.return_moment_tolerance;
        c.collision_kernel = copy_string_out(
            payloads[0].header.collision_kernel,
            sizeof(payloads[0].header.collision_kernel));
        c.collision_weight_mode = copy_string_out(
            payloads[0].header.collision_weight_mode,
            sizeof(payloads[0].header.collision_weight_mode));
        c.collision_max_substeps = payloads[0].header.collision_max_substeps;
        c.collision_max_particle_growth =
            payloads[0].header.collision_max_particle_growth;
        c.conversion_cumulative_number =
            payloads[0].header.conversion_cumulative_number;
        c.conversion_cumulative_px = payloads[0].header.conversion_cumulative_px;
        c.conversion_cumulative_energy =
            payloads[0].header.conversion_cumulative_energy;
        c.conversion_cumulative_particles_created = static_cast<std::uint64_t>(
            payloads[0].header.conversion_cumulative_particles_created);
        c.tail_cumulative_outflow_number = static_cast<std::uint64_t>(
            payloads[0].header.tail_cumulative_outflow_number);
        c.control_cumulative_groups = static_cast<std::uint64_t>(
            payloads[0].header.control_cumulative_groups);
        c.control_cumulative_fallbacks = static_cast<std::uint64_t>(
            payloads[0].header.control_cumulative_fallbacks);
        c.return_cumulative_number =
            payloads[0].header.return_cumulative_number;
        c.return_cumulative_px = payloads[0].header.return_cumulative_px;
        c.return_cumulative_jx_dx =
            payloads[0].header.return_cumulative_jx_dx;
        c.return_cumulative_energy =
            payloads[0].header.return_cumulative_energy;
        c.return_cumulative_pixx_dx =
            payloads[0].header.return_cumulative_pixx_dx;
        c.return_cumulative_piperp_dx =
            payloads[0].header.return_cumulative_piperp_dx;
        c.return_cumulative_particles_removed = static_cast<std::uint64_t>(
            payloads[0].header.return_cumulative_particles_removed);
        c.return_cumulative_deferred_groups = static_cast<std::uint64_t>(
            payloads[0].header.return_cumulative_deferred_groups);
        c.combined_number = payloads[0].header.combined_number;
        c.combined_kinetic_energy = payloads[0].header.combined_kinetic_energy;
        c.combined_field_energy = payloads[0].header.combined_field_energy;
        if (payloads[0].header.version == 1) {
            // The legacy schema has no conversion metadata.  It is only a
            // valid static-cell checkpoint; flux modes must be rejected by
            // the caller's configuration validation.
            c.conversion_mode = "static-cell";
            c.conversion_metadata_present = true;
        }
        if (!read_tail_conversion_manifest(directory, c)) {
            error = "cannot read VPFP conversion manifest during repartition";
            return false;
        }
    }
    electrons.f_tmp = electrons.f;
    electrons.compute_moments();
    control = payloads[0].header.control;
    return true;
}

bool read_tail_conversion_manifest(const std::string& directory,
                                   VpfpCheckpointTailConfig& config)
{
    std::ifstream in((directory + "/manifest.txt").c_str());
    if (!in) return true; // legacy rank files have no flux-mode metadata
    std::string key;
    std::size_t edge_count = 0;
    while (in >> key) {
        if (key == "tail_conversion_mode") {
            in >> config.conversion_mode;
            config.conversion_metadata_present = true;
        }
        else if (key == "physical_config_hash") {
            in >> config.physical_config_hash;
            config.conversion_metadata_present = true;
        }
        else if (key == "diagnostic_config_hash") {
            in >> config.diagnostic_config_hash;
        }
        else if (key == "tail_flux_quadrature_order") {
            in >> config.flux_quadrature_order;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_flux_max_supports") {
            in >> config.flux_max_supports;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_flux_max_created_particles_per_step")
        {
            in >> config.flux_max_created_particles_per_step;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_interface_topology_hash") {
            in >> config.interface_topology_hash;
            config.interface_topology_metadata_present = true;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_interface_topology_version") {
            in >> config.interface_topology_version;
            config.interface_topology_detail_metadata_present = true;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_interface_mask_hash") {
            in >> config.interface_mask_hash;
            config.interface_topology_detail_metadata_present = true;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_interface_face_list_hash") {
            in >> config.interface_face_list_hash;
            config.interface_topology_detail_metadata_present = true;
            config.conversion_metadata_present = true;
        }
        else if (key == "tail_convert_energy_mev") {
            in >> config.convert_energy_mev;
        }
        else if (key == "tail_return_mode") {
            in >> config.return_mode;
        }
        else if (key == "tail_return_energy_mev") {
            in >> config.return_energy_mev;
        }
        else if (key == "tail_return_residence_steps") {
            in >> config.return_residence_steps;
        }
        else if (key == "tail_return_max_stencil_radius") {
            in >> config.return_max_stencil_radius;
        }
        else if (key == "tail_return_moment_tolerance") {
            in >> config.return_moment_tolerance;
        }
        else if (key == "conversion_cumulative_number") {
            in >> config.conversion_cumulative_number;
        }
        else if (key == "conversion_cumulative_px") {
            in >> config.conversion_cumulative_px;
        }
        else if (key == "conversion_cumulative_energy") {
            in >> config.conversion_cumulative_energy;
        }
        else if (key == "conversion_cumulative_particles_created") {
            in >> config.conversion_cumulative_particles_created;
        }
        else if (key == "return_cumulative_number") {
            in >> config.return_cumulative_number;
        }
        else if (key == "return_cumulative_px") {
            in >> config.return_cumulative_px;
        }
        else if (key == "return_cumulative_jx_dx") {
            in >> config.return_cumulative_jx_dx;
        }
        else if (key == "return_cumulative_energy") {
            in >> config.return_cumulative_energy;
        }
        else if (key == "return_cumulative_pixx_dx") {
            in >> config.return_cumulative_pixx_dx;
        }
        else if (key == "return_cumulative_piperp_dx") {
            in >> config.return_cumulative_piperp_dx;
        }
        else if (key == "return_cumulative_particles_removed") {
            in >> config.return_cumulative_particles_removed;
        }
        else if (key == "return_cumulative_deferred_groups") {
            in >> config.return_cumulative_deferred_groups;
        }
        else if (key == "tail_conversion_bins") {
            std::string bins;
            in >> bins;
            const std::size_t separator = bins.find('x');
            if (separator != std::string::npos) {
                config.upar_bins = std::atoi(
                    bins.substr(0, separator).c_str());
                config.energy_bins = std::atoi(
                    bins.substr(separator + 1).c_str());
            }
        }
        else if (key == "conversion_energy_edge_count") {
            in >> edge_count;
        }
        else if (key == "conversion_energy_edges") {
            config.conversion_energy_edges.resize(edge_count);
            for (std::size_t i = 0; i < edge_count; ++i)
                in >> config.conversion_energy_edges[i];
        }
        else if (key == "conversion_energy_edges_hash") {
            in >> config.conversion_energy_edges_hash;
            config.conversion_metadata_present = true;
        }
        else if (key == "collision_interface_mode") {
            in >> config.collision_interface_mode;
        }
        else if (key == "bulk_collision_integrator") {
            in >> config.bulk_collision_integrator;
        }
        else if (key == "collision_induced_conversion") {
            int value = 1;
            in >> value;
            config.collision_induced_conversion = value != 0;
        }
        else {
            std::string ignored;
            std::getline(in, ignored);
        }
    }
    return static_cast<bool>(in) || in.eof();
}

bool fill_tail_header(const VpfpCheckpointTailState* tail_state,
                      RankHeader& header)
{
    if (tail_state == NULL) {
        header.tail_present = 0;
        return true;
    }
    if (!tail_state->present) return false;
    header.tail_present = 1;
    header.tail_particle_count =
        static_cast<unsigned long long>(tail_state->tail.particles.size());
    header.tail_density_size =
        static_cast<unsigned long long>(tail_state->tail.density.size());
    const BackgroundTailStateSnapshot snapshot;
    // export_accepted_state fills the snapshot; the header fields below are
    // copied from it.
    BackgroundTailStateSnapshot s = snapshot;
    tail_state->tail.export_accepted_state(s);
    header.tail_id_counter = static_cast<unsigned long long>(s.id_counter);
    header.tail_collision_rng_seed =
        static_cast<unsigned long long>(s.collision_rng_seed);
    header.tail_outflow_left_number = s.outflow.left_number;
    header.tail_outflow_left_px = s.outflow.left_px;
    header.tail_outflow_left_ke = s.outflow.left_kinetic_energy;
    header.tail_outflow_right_number = s.outflow.right_number;
    header.tail_outflow_right_px = s.outflow.right_px;
    header.tail_outflow_right_ke = s.outflow.right_kinetic_energy;
    header.tail_truncation_shape_left = s.truncation_shape_left;
    header.tail_truncation_shape_right = s.truncation_shape_right;
    header.tail_deposit_shape_left = s.deposit_shape_left;
    header.tail_deposit_shape_right = s.deposit_shape_right;
    header.tail_deposit_shape_step_start_left =
        s.deposit_shape_step_start_left;
    header.tail_deposit_shape_step_start_right =
        s.deposit_shape_step_start_right;
    header.tail_max_abs_u = s.max_abs_u;
    header.tail_max_kinetic_energy = s.max_kinetic_energy;

    const VpfpCheckpointTailConfig& c = tail_state->config;
    header.partition_config_hash =
        static_cast<unsigned long long>(c.partition_config_hash);
    header.convert_energy_mev = c.convert_energy_mev;
    header.buffer_width_mev = c.buffer_width_mev;
    header.upar_bins = c.upar_bins;
    header.energy_bins = c.energy_bins;
    header.population_control_enabled = c.population_control_enabled ? 1 : 0;
    header.control_interval = c.control_interval;
    header.target_particles_per_phase_bin = c.target_particles_per_phase_bin;
    header.max_particles_per_phase_bin = c.max_particles_per_phase_bin;
    header.max_weight_ratio = c.max_weight_ratio;
    header.max_support = c.max_support;
    copy_string_into(c.return_mode, header.return_mode,
                     sizeof(header.return_mode));
    header.return_energy_mev = c.return_energy_mev;
    header.return_residence_steps = c.return_residence_steps;
    header.return_max_stencil_radius = c.return_max_stencil_radius;
    header.return_moment_tolerance = c.return_moment_tolerance;
    copy_string_into(c.collision_kernel, header.collision_kernel,
                     sizeof(header.collision_kernel));
    copy_string_into(c.collision_weight_mode, header.collision_weight_mode,
                     sizeof(header.collision_weight_mode));
    header.collision_max_substeps = c.collision_max_substeps;
    header.collision_max_particle_growth = c.collision_max_particle_growth;
    header.conversion_cumulative_number = c.conversion_cumulative_number;
    header.conversion_cumulative_px = c.conversion_cumulative_px;
    header.conversion_cumulative_energy = c.conversion_cumulative_energy;
    header.conversion_cumulative_particles_created =
        static_cast<unsigned long long>(
            c.conversion_cumulative_particles_created);
    header.tail_cumulative_outflow_number =
        static_cast<unsigned long long>(c.tail_cumulative_outflow_number);
    header.control_cumulative_groups =
        static_cast<unsigned long long>(c.control_cumulative_groups);
    header.control_cumulative_fallbacks =
        static_cast<unsigned long long>(c.control_cumulative_fallbacks);
    header.return_cumulative_number = c.return_cumulative_number;
    header.return_cumulative_px = c.return_cumulative_px;
    header.return_cumulative_jx_dx = c.return_cumulative_jx_dx;
    header.return_cumulative_energy = c.return_cumulative_energy;
    header.return_cumulative_pixx_dx = c.return_cumulative_pixx_dx;
    header.return_cumulative_piperp_dx = c.return_cumulative_piperp_dx;
    header.return_cumulative_particles_removed =
        static_cast<unsigned long long>(c.return_cumulative_particles_removed);
    header.return_cumulative_deferred_groups =
        static_cast<unsigned long long>(c.return_cumulative_deferred_groups);
    header.combined_number = c.combined_number;
    header.combined_kinetic_energy = c.combined_kinetic_energy;
    header.combined_field_energy = c.combined_field_energy;
    return true;
}

} // namespace

bool read_vpfp_tail_particle_records(
    std::istream& input, int checkpoint_version, size_t count,
    std::vector<BackgroundTailParticle>& particles)
{
    return read_tail_particles(input, checkpoint_version, count, particles);
}

bool validate_vpfp_checkpoint_tail_config(
    const VpfpCheckpointTailConfig& stored,
    const VpfpCheckpointTailConfig& expected,
    std::string& error)
{
    if (!stored.conversion_metadata_present) {
        error = "checkpoint has no tail conversion metadata";
        return false;
    }
    if (stored.conversion_mode != expected.conversion_mode) {
        error = "tail conversion mode mismatch";
        return false;
    }
    if (stored.flux_quadrature_order != expected.flux_quadrature_order) {
        error = "tail flux quadrature order mismatch";
        return false;
    }
    if (stored.flux_max_supports != expected.flux_max_supports) {
        error = "tail flux support limit mismatch";
        return false;
    }
    if (stored.flux_max_created_particles_per_step !=
        expected.flux_max_created_particles_per_step) {
        error = "tail flux particle limit mismatch";
        return false;
    }
    if (expected.convert_energy_mev > 0.0 &&
        stored.convert_energy_mev != expected.convert_energy_mev) {
        error = "tail conversion energy threshold mismatch";
        return false;
    }
    if (expected.upar_bins > 0 && stored.upar_bins != expected.upar_bins) {
        error = "tail conversion u_parallel bin mismatch";
        return false;
    }
    if (expected.energy_bins > 0 &&
        stored.energy_bins != expected.energy_bins) {
        error = "tail conversion energy bin mismatch";
        return false;
    }
    if (expected.conversion_energy_edges_hash != 0 &&
        stored.conversion_energy_edges_hash !=
            expected.conversion_energy_edges_hash) {
        error = "tail conversion energy edge hash mismatch";
        return false;
    }
    if (!expected.conversion_energy_edges.empty() &&
        stored.conversion_energy_edges != expected.conversion_energy_edges) {
        error = "tail conversion energy edges mismatch";
        return false;
    }
    if (stored.collision_interface_mode != expected.collision_interface_mode ||
        stored.bulk_collision_integrator != expected.bulk_collision_integrator ||
        stored.collision_induced_conversion !=
            expected.collision_induced_conversion) {
        error = "collision interface configuration mismatch";
        return false;
    }
    if (stored.partition_config_hash != expected.partition_config_hash) {
        error = "tail partition configuration mismatch";
        return false;
    }
    if (expected.physical_config_hash != 0 &&
        stored.physical_config_hash != expected.physical_config_hash) {
        std::ostringstream message;
        message << "physical checkpoint configuration hash mismatch"
                << " stored=" << stored.physical_config_hash
                << " expected=" << expected.physical_config_hash;
        error = message.str();
        return false;
    }
    if (expected.interface_topology_metadata_present &&
        (!stored.interface_topology_metadata_present ||
         stored.interface_topology_hash != expected.interface_topology_hash)) {
        error = "tail interface topology mismatch";
        return false;
    }
    if (expected.interface_topology_detail_metadata_present &&
        (!stored.interface_topology_detail_metadata_present ||
         stored.interface_topology_version != expected.interface_topology_version ||
         stored.interface_mask_hash != expected.interface_mask_hash ||
         stored.interface_face_list_hash != expected.interface_face_list_hash)) {
        error = "tail interface topology detail mismatch";
        return false;
    }
    return true;
}

// JC4 (section 7.5): read field-particle coupling configuration from a
// checkpoint manifest.  Returns true if the file was readable; defaults to
// legacy when the manifest is absent or lacks the coupling keys (old
// checkpoint).
bool read_coupling_config_from_manifest(
    const std::string& directory,
    std::string& coupling_mode, int& coupling_max_iters,
    double& coupling_relaxation, double& coupling_field_tol,
    double& coupling_pairing_tol,
    std::string& x_transport_velocity_mode,
    int& x_transport_velocity_table_schema)
{
    coupling_mode = "legacy";
    coupling_max_iters = 12;
    coupling_relaxation = 0.5;
    coupling_field_tol = 1.0e-8;
    coupling_pairing_tol = 1.0e-8;
    x_transport_velocity_mode = "analytic-cell-center";
    // Old manifests predate P3-V.2.  Missing metadata means the analytic
    // default, but is distinguishable so its historical physical hash can be
    // validated without silently accepting energy-conjugate mode.
    x_transport_velocity_table_schema = 0;
    std::ifstream in((directory + "/manifest.txt").c_str());
    if (!in) return true; // no manifest → legacy defaults
    std::string key;
    while (in >> key) {
        if (key == "coupling_mode") in >> coupling_mode;
        else if (key == "coupling_max_iters") in >> coupling_max_iters;
        else if (key == "coupling_relaxation") in >> coupling_relaxation;
        else if (key == "coupling_field_tol") in >> coupling_field_tol;
        else if (key == "coupling_pairing_tol") in >> coupling_pairing_tol;
        else if (key == "x_transport_velocity_mode")
            in >> x_transport_velocity_mode;
        else if (key == "x_transport_velocity_table_schema")
            in >> x_transport_velocity_table_schema;
    }
    return true;
}

bool read_background_phase_space_mode_from_manifest(
    const std::string& directory, std::string& mode)
{
    mode = "strang-ppm";
    std::ifstream in((directory + "/manifest.txt").c_str());
    if (!in) return true;
    std::string key;
    while (in >> key) {
        if (key == "background_phase_space_mode") {
            in >> mode;
            return static_cast<bool>(in);
        }
        std::string ignored;
        std::getline(in, ignored);
    }
    return true;
}

bool vpfp_preconversion_static_to_flux_restart_allowed(
    const VpfpCheckpointTailConfig& stored,
    const std::string& requested_conversion_mode,
    std::uint64_t global_tail_particle_count)
{
    return stored.conversion_metadata_present &&
           stored.conversion_mode == "static-cell" &&
           requested_conversion_mode == "flux-interface" &&
           global_tail_particle_count == 0 &&
           stored.conversion_cumulative_number == 0.0 &&
           stored.conversion_cumulative_px == 0.0 &&
           stored.conversion_cumulative_energy == 0.0 &&
           stored.conversion_cumulative_particles_created == 0 &&
           stored.tail_cumulative_outflow_number == 0 &&
           stored.control_cumulative_groups == 0 &&
           stored.control_cumulative_fallbacks == 0;
}

bool write_vpfp_checkpoint(
    const std::string& directory, const VpfpCheckpointControl& control,
    const Species& electrons, const BeamPIC& beam, const EMFields& fields,
    const SpatialGrid& grid, const ElectrostaticBoundary& field_boundary,
    const OpenBackgroundBoundaryConfig& background_boundary,
    const std::string& collision_model,
    const VpfpCheckpointTailState* tail_state,
    const VpfpCouplingManifestConfig& coupling_config,
    int mpi_rank, int mpi_size, std::string& error)
{
    int directory_ok = create_directories(directory) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &directory_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!directory_ok) {
        error = "cannot create VPFP checkpoint directory: " + directory;
        return false;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    RankHeader header = {};
    header.magic = checkpoint_magic;
    header.version = checkpoint_version;
    header.nx_global = grid.nx_global;
    header.nx_local = grid.nx_local;
    header.nv = Param::Nv;
    header.nmu = Param::Nmu;
    header.rank = mpi_rank;
    header.mpi_size = mpi_size;
    header.control = control;
    header.beam_state = beam.export_persistent_state();
    header.f_size = electrons.f.size();
    header.particle_count = beam.particles.size();
    header.ex_face_size = fields.Ex_face.size();
    header.ex_size = fields.Ex.size();
    header.phi_size = fields.phi.size();
    if (!fill_tail_header(tail_state, header)) {
        error = "invalid tail checkpoint state (present=false)";
        return false;
    }
    const std::string final_rank_file = rank_file(directory, mpi_rank);
    const std::string temporary_rank_file = final_rank_file + ".tmp";
    std::ofstream out(temporary_rank_file.c_str(),
                      std::ios::binary | std::ios::trunc);
    if (!out) { error = "cannot create VPFP rank checkpoint"; return false; }
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));
    out.write(reinterpret_cast<const char*>(electrons.f.data()),
              static_cast<std::streamsize>(electrons.f.size() * sizeof(double)));
    out.write(reinterpret_cast<const char*>(beam.particles.data()),
              static_cast<std::streamsize>(beam.particles.size() * sizeof(BeamParticle)));
    if (header.tail_present) {
        if (!write_tail_particles(out, tail_state->tail.particles)) {
            error = "cannot write VPFP tail particle payload";
            return false;
        }
        out.write(reinterpret_cast<const char*>(tail_state->tail.density.data()),
                  static_cast<std::streamsize>(
                      tail_state->tail.density.size() * sizeof(double)));
    }
    out.write(reinterpret_cast<const char*>(fields.Ex_face.data()),
              static_cast<std::streamsize>(fields.Ex_face.size() * sizeof(double)));
    out.write(reinterpret_cast<const char*>(fields.Ex.data()),
              static_cast<std::streamsize>(fields.Ex.size() * sizeof(double)));
    out.write(reinterpret_cast<const char*>(fields.phi.data()),
              static_cast<std::streamsize>(fields.phi.size() * sizeof(double)));
    const bool stream_ok = static_cast<bool>(out);
    out.close();
    int ok = stream_ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if (!ok) { error = "VPFP rank checkpoint write failed"; return false; }
    std::remove(final_rank_file.c_str());
    int rank_commit_ok =
        std::rename(temporary_rank_file.c_str(), final_rank_file.c_str()) == 0 ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &rank_commit_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!rank_commit_ok) { error = "cannot commit VPFP rank checkpoint"; return false; }
    MPI_Barrier(MPI_COMM_WORLD);
    int manifest_ok_global = 1;
    if (mpi_rank == 0) {
        const std::string final_manifest = directory + "/manifest.txt";
        const std::string temporary_manifest = final_manifest + ".tmp";
        std::ofstream manifest(temporary_manifest.c_str(), std::ios::trunc);
        manifest << "solver_kind vpfp-open-v4\nversion " << checkpoint_version
                 << "\nmpi_size " << mpi_size << "\nnx_global " << grid.nx_global
                 << "\ndx " << std::setprecision(17) << grid.dx
                 << "\nfield_boundary " << static_cast<int>(field_boundary.type)
                 << "\ne_left " << field_boundary.e_left
                 << "\nphi_left " << field_boundary.phi_left
                 << "\nphi_right " << field_boundary.phi_right
                 << "\nbackground_left_type "
                 << static_cast<int>(background_boundary.left_type)
                 << "\nbackground_right_type "
                 << static_cast<int>(background_boundary.right_type)
                 << "\ncollision_model " << collision_model
                 << "\nstep " << control.step << "\ntime " << control.time
                 << "\ndt " << control.dt
                 // JC4 (section 7.2/7.5): field-particle coupling config in
                 // the common manifest section so it is written for ALL runs
                 // (tail and non-tail).  Old checkpoints missing these keys
                 // are parsed as legacy defaults.
                 << "\ncoupling_mode " << coupling_config.mode
                 << "\ncoupling_max_iters " << coupling_config.max_iters
                 << "\ncoupling_relaxation " << std::setprecision(17)
                 << coupling_config.relaxation
                  << "\ncoupling_field_tol " << coupling_config.field_tol
                  << "\ncoupling_pairing_tol " << coupling_config.pairing_tol
                  << "\nbackground_phase_space_mode "
                  << coupling_config.background_phase_space_mode
                  << "\nx_transport_velocity_mode "
                  << coupling_config.x_transport_velocity_mode
                  << "\nx_transport_velocity_table_schema "
                  << coupling_config.x_transport_velocity_table_schema;
        if (header.tail_present) {
            const VpfpCheckpointTailConfig& c = tail_state->config;
            const bool coulomb =
                c.collision_kernel == "coulomb-nanbu-perez";
            const bool sde = c.collision_kernel == "kramers-moyal-sde";
            const bool tail_bulk = coulomb || sde;
            manifest
                << "\nbackground_representation eulerian_bulk_plus_pic_tail"
                << "\ntail_return_mode " << c.return_mode
                << "\ntail_return_energy_mev " << c.return_energy_mev
                << "\ntail_return_residence_steps "
                << c.return_residence_steps
                << "\ntail_return_max_stencil_radius "
                << c.return_max_stencil_radius
                << "\ntail_return_moment_tolerance "
                << c.return_moment_tolerance
                << "\ntail_convert_energy_mev "
                << std::setprecision(17) << c.convert_energy_mev
                << "\ntail_conversion_bins " << c.upar_bins << "x"
                << c.energy_bins
                << "\ntail_conversion_mode " << c.conversion_mode
                // These state counters make a checkpoint independently
                // classifiable as pre- or post-conversion by offline tools.
                << "\nconversion_cumulative_number "
                << c.conversion_cumulative_number
                << "\nconversion_cumulative_px "
                << c.conversion_cumulative_px
                << "\nconversion_cumulative_energy "
                << c.conversion_cumulative_energy
                << "\nconversion_cumulative_particles_created "
                << c.conversion_cumulative_particles_created
                << "\nreturn_cumulative_number " << c.return_cumulative_number
                << "\nreturn_cumulative_px " << c.return_cumulative_px
                << "\nreturn_cumulative_jx_dx " << c.return_cumulative_jx_dx
                << "\nreturn_cumulative_energy " << c.return_cumulative_energy
                << "\nreturn_cumulative_pixx_dx " << c.return_cumulative_pixx_dx
                << "\nreturn_cumulative_piperp_dx " << c.return_cumulative_piperp_dx
                << "\nreturn_cumulative_particles_removed "
                << c.return_cumulative_particles_removed
                << "\nreturn_cumulative_deferred_groups "
                << c.return_cumulative_deferred_groups
                << "\nphysical_config_hash " << c.physical_config_hash
                << "\ndiagnostic_config_hash " << c.diagnostic_config_hash
                << "\ntail_flux_quadrature_order " << c.flux_quadrature_order
                << "\ntail_flux_max_supports " << c.flux_max_supports
                << "\ntail_flux_max_created_particles_per_step "
                << c.flux_max_created_particles_per_step
                << "\ntail_interface_topology_hash "
                << c.interface_topology_hash
                << "\ntail_interface_topology_version "
                << c.interface_topology_version
                << "\ntail_interface_mask_hash "
                << c.interface_mask_hash
                << "\ntail_interface_face_list_hash "
                << c.interface_face_list_hash
                << "\ncollision_interface_mode "
                << c.collision_interface_mode
                << "\nbulk_collision_integrator "
                << c.bulk_collision_integrator
                << "\ncollision_induced_conversion "
                << (c.collision_induced_conversion ? 1 : 0)
                << "\npartition_config_hash " << c.partition_config_hash
                << "\nconversion_energy_edge_count "
                << c.conversion_energy_edges.size()
                << "\nconversion_energy_edges";
            for (size_t i = 0; i < c.conversion_energy_edges.size(); ++i) {
                manifest << " " << std::setprecision(17)
                         << c.conversion_energy_edges[i];
            }
            manifest << "\nconversion_energy_edges_hash "
                     << c.conversion_energy_edges_hash
                << "\ntail_particle_schema background_tail_particle_v2"
                << "\ntail_pusher_backend relativistic_dkd_v1"
                << "\ntail_deposition_backend cic_v1"
                << "\ntail_collision_backend " << c.collision_kernel
                << "\ntail_tail_collision_backend "
                << (coulomb ? "coulomb-nanbu-perez" : "none")
                << "\ntail_bulk_collision_backend "
                << (tail_bulk ? "kramers-moyal-sde" : "none")
                << "\ntail_collision_weight_mode " << c.collision_weight_mode
                << "\ntail_collision_weight_algorithm "
                << (c.collision_weight_mode == "virtual-split"
                        ? "sentoku-kemp-bounded-v1"
                        : "equal-strata-exact-v1")
                << "\ntail_collision_max_substeps "
                << c.collision_max_substeps
                << "\ntail_collision_max_particle_growth "
                << c.collision_max_particle_growth
                << "\ncollision_pair_bulk_bulk "
                << (collision_model == "zero" || collision_model == "none"
                        ? 0 : 1)
                << "\ncollision_pair_bulk_tail " << (tail_bulk ? 1 : 0)
                << "\ncollision_pair_bulk_reaction "
                << (tail_bulk ? 1 : 0)
                << "\ncollision_pair_tail_bulk " << (tail_bulk ? 1 : 0)
                << "\ncollision_pair_tail_tail " << (coulomb ? 1 : 0)
                << "\npopulation_control_interval " << c.control_interval
                << "\npopulation_control_target_per_bin "
                << c.target_particles_per_phase_bin
                << "\npopulation_control_max_per_bin "
                << c.max_particles_per_phase_bin
                << "\npopulation_control_max_weight_ratio "
                << c.max_weight_ratio
                << "\ncombined_number " << c.combined_number
                << "\ncombined_kinetic_energy " << c.combined_kinetic_energy
                << "\ncombined_field_energy " << c.combined_field_energy;
        } else {
            manifest << "\nbackground_representation eulerian_only";
        }
        manifest << "\nix_start " << grid.ix_start << "\nnx_local "
                 << grid.nx_local << "\n";
        const bool manifest_ok = static_cast<bool>(manifest);
        manifest.close();
        std::remove(final_manifest.c_str());
        if (!manifest_ok ||
            std::rename(temporary_manifest.c_str(), final_manifest.c_str()) != 0)
            manifest_ok_global = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &manifest_ok_global, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (!manifest_ok_global) { error = "cannot commit VPFP checkpoint manifest"; return false; }
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
}

bool read_vpfp_checkpoint(
    const std::string& directory, VpfpCheckpointControl& control,
    Species& electrons, BeamPIC& beam, EMFields& fields,
    const SpatialGrid& grid, VpfpCheckpointTailState* tail_state,
    int mpi_rank, int mpi_size, std::string& error)
{
    std::ifstream in(rank_file(directory, mpi_rank).c_str(), std::ios::binary);
    RankHeader header = {};
    if (!in || !read_rank_header(in, header)) {
        error = "cannot read VPFP rank checkpoint";
        return false;
    }
    if (header.magic != checkpoint_magic ||
        (header.version != 1 && header.version != 2 && header.version != 3 &&
         header.version != checkpoint_version) ||
        header.nx_global != grid.nx_global || header.nv != Param::Nv ||
        header.nmu != Param::Nmu) {
        error = "VPFP checkpoint incompatible with solver kind/grid/MPI layout";
        return false;
    }
    if (header.mpi_size != mpi_size || header.nx_local != grid.nx_local) {
        in.close();
        return read_repartitioned_checkpoint(
            directory, control, electrons, beam, fields, grid, tail_state,
            mpi_rank, mpi_size, error);
    }
    if (header.rank != mpi_rank) {
        error = "VPFP checkpoint rank file does not match requested rank";
        return false;
    }
    if (header.f_size != electrons.f.size() ||
        header.ex_face_size != fields.Ex_face.size() ||
        header.ex_size != fields.Ex.size() ||
        header.phi_size != fields.phi.size()) {
        error = "VPFP checkpoint array contract mismatch";
        return false;
    }

    // Section 12.2: an old no-tail checkpoint is never silently upgraded to
    // a hybrid one, and a checkpoint carrying a tail is refused by a
    // tail-disabled solver.
    if (header.version == 1) {
        if (tail_state != NULL) {
            error = "old no-tail checkpoint cannot be upgraded to a hybrid "
                    "checkpoint (section 12.2)";
            return false;
        }
    } else if (header.tail_present) {
        if (tail_state == NULL) {
            error = "checkpoint contains a tail but the solver runs tail-off";
            return false;
        }
    } else {
        if (tail_state != NULL) {
            error = "checkpoint has no tail state (section 12.2: no silent "
                    "upgrade)";
            return false;
        }
    }

    beam.particles.resize(static_cast<size_t>(header.particle_count));
    in.read(reinterpret_cast<char*>(electrons.f.data()),
            static_cast<std::streamsize>(electrons.f.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(beam.particles.data()),
            static_cast<std::streamsize>(beam.particles.size() * sizeof(BeamParticle)));
    if (header.tail_present) {
        if (header.tail_density_size !=
            static_cast<unsigned long long>(grid.nx_local)) {
            error = "VPFP checkpoint tail density contract mismatch";
            return false;
        }
        BackgroundTailStateSnapshot snapshot;
        snapshot.particles.resize(static_cast<size_t>(header.tail_particle_count));
        snapshot.density.resize(static_cast<size_t>(header.tail_density_size));
        if (!read_tail_particles(in, header.version, snapshot.particles.size(),
                                 snapshot.particles)) {
            error = "truncated VPFP tail particle payload";
            return false;
        }
        in.read(reinterpret_cast<char*>(snapshot.density.data()),
                static_cast<std::streamsize>(
                    snapshot.density.size() * sizeof(double)));
        snapshot.id_counter = static_cast<std::uint64_t>(header.tail_id_counter);
        snapshot.collision_rng_seed =
            static_cast<std::uint64_t>(header.tail_collision_rng_seed);
        snapshot.outflow.left_number = header.tail_outflow_left_number;
        snapshot.outflow.left_px = header.tail_outflow_left_px;
        snapshot.outflow.left_kinetic_energy = header.tail_outflow_left_ke;
        snapshot.outflow.right_number = header.tail_outflow_right_number;
        snapshot.outflow.right_px = header.tail_outflow_right_px;
        snapshot.outflow.right_kinetic_energy = header.tail_outflow_right_ke;
        snapshot.truncation_shape_left = header.tail_truncation_shape_left;
        snapshot.truncation_shape_right = header.tail_truncation_shape_right;
        snapshot.deposit_shape_left = header.tail_deposit_shape_left;
        snapshot.deposit_shape_right = header.tail_deposit_shape_right;
        snapshot.deposit_shape_step_start_left =
            header.tail_deposit_shape_step_start_left;
        snapshot.deposit_shape_step_start_right =
            header.tail_deposit_shape_step_start_right;
        snapshot.max_abs_u = header.tail_max_abs_u;
        snapshot.max_kinetic_energy = header.tail_max_kinetic_energy;
        tail_state->present = true;
        tail_state->tail.init(grid);
        tail_state->tail.import_accepted_state(snapshot);
        VpfpCheckpointTailConfig& c = tail_state->config;
        c.partition_config_hash =
            static_cast<std::uint64_t>(header.partition_config_hash);
        c.convert_energy_mev = header.convert_energy_mev;
        c.buffer_width_mev = header.buffer_width_mev;
        c.upar_bins = header.upar_bins;
        c.energy_bins = header.energy_bins;
        c.population_control_enabled =
            header.population_control_enabled != 0;
        c.control_interval = header.control_interval;
        c.target_particles_per_phase_bin =
            header.target_particles_per_phase_bin;
        c.max_particles_per_phase_bin = header.max_particles_per_phase_bin;
        c.max_weight_ratio = header.max_weight_ratio;
        c.max_support = header.max_support;
        c.return_mode = copy_string_out(header.return_mode,
                                        sizeof(header.return_mode));
        c.return_energy_mev = header.return_energy_mev;
        c.return_residence_steps = header.return_residence_steps;
        c.return_max_stencil_radius = header.return_max_stencil_radius;
        c.return_moment_tolerance = header.return_moment_tolerance;
        c.collision_kernel = copy_string_out(header.collision_kernel,
                                             sizeof(header.collision_kernel));
        c.collision_weight_mode = copy_string_out(
            header.collision_weight_mode, sizeof(header.collision_weight_mode));
        c.collision_max_substeps = header.collision_max_substeps;
        c.collision_max_particle_growth = header.collision_max_particle_growth;
        c.conversion_cumulative_number = header.conversion_cumulative_number;
        c.conversion_cumulative_px = header.conversion_cumulative_px;
        c.conversion_cumulative_energy = header.conversion_cumulative_energy;
        c.conversion_cumulative_particles_created =
            static_cast<std::uint64_t>(
                header.conversion_cumulative_particles_created);
        c.tail_cumulative_outflow_number =
            static_cast<std::uint64_t>(header.tail_cumulative_outflow_number);
        c.control_cumulative_groups =
            static_cast<std::uint64_t>(header.control_cumulative_groups);
        c.control_cumulative_fallbacks =
            static_cast<std::uint64_t>(header.control_cumulative_fallbacks);
        c.return_cumulative_number = header.return_cumulative_number;
        c.return_cumulative_px = header.return_cumulative_px;
        c.return_cumulative_jx_dx = header.return_cumulative_jx_dx;
        c.return_cumulative_energy = header.return_cumulative_energy;
        c.return_cumulative_pixx_dx = header.return_cumulative_pixx_dx;
        c.return_cumulative_piperp_dx = header.return_cumulative_piperp_dx;
        c.return_cumulative_particles_removed = static_cast<std::uint64_t>(
            header.return_cumulative_particles_removed);
        c.return_cumulative_deferred_groups = static_cast<std::uint64_t>(
            header.return_cumulative_deferred_groups);
        c.combined_number = header.combined_number;
        c.combined_kinetic_energy = header.combined_kinetic_energy;
        c.combined_field_energy = header.combined_field_energy;
        if (!read_tail_conversion_manifest(directory, c)) {
            error = "cannot read VPFP checkpoint conversion manifest";
            return false;
        }
    }
    in.read(reinterpret_cast<char*>(fields.Ex_face.data()),
            static_cast<std::streamsize>(fields.Ex_face.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(fields.Ex.data()),
            static_cast<std::streamsize>(fields.Ex.size() * sizeof(double)));
    in.read(reinterpret_cast<char*>(fields.phi.data()),
            static_cast<std::streamsize>(fields.phi.size() * sizeof(double)));
    if (!in) { error = "truncated VPFP checkpoint"; return false; }
    beam.import_persistent_state(header.beam_state, grid);
    electrons.f_tmp = electrons.f;
    electrons.compute_moments();
    control = header.control;
    return true;
}
