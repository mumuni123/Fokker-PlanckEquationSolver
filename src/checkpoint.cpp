#include "checkpoint.h"

#include "parameters.h"

#include <cerrno>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

namespace {
const unsigned int kVersion = 1U;
const unsigned int kEndian = 0x01020304U;
const char kMagic[] = "FPCHKPT1";

bool make_dir(const std::string& path)
{
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}
std::string rank_path(const std::string& d, int r) { char b[64]; std::sprintf(b, "/rank_%06d.bin", r); return d + b; }
bool finite_vector(const std::vector<double>& values) { for (size_t i=0;i<values.size();++i) if (!std::isfinite(values[i])) return false; return true; }
bool finite_particles(const std::vector<BeamParticle>& particles) { for (size_t i=0;i<particles.size();++i) if (!std::isfinite(particles[i].x) || !std::isfinite(particles[i].px) || !std::isfinite(particles[i].weight)) return false; return true; }
bool finite_persistent_state(const BeamPersistentState& s)
{
    const double values[] = {s.injection_remainder, s.cumulative_injected_energy,
        s.cumulative_outflow_energy, s.last_injected_energy, s.last_outflow_energy,
        s.last_injected_number, s.last_outflow_number, s.last_injected_current,
        s.last_outflow_current, s.last_field_work, s.step_dt,
        s.step_signed_outflow_number, s.interval_injected_number,
        s.interval_left_outflow_signed_number, s.interval_right_outflow_number,
        s.interval_left_guard_path_number, s.interval_right_guard_path_number,
        s.last_continuity_abs_l1_residual, s.last_continuity_abs_linf_residual,
        s.last_continuity_l1_error, s.last_continuity_linf_error,
        s.last_boundary_flux_error, s.last_trajectory_reconstruction_error};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
        if (!std::isfinite(values[i])) return false;
    return true;
}
bool rename_atomic(const std::string& from, const std::string& to) { return std::rename(from.c_str(), to.c_str()) == 0; }
}

unsigned long long checkpoint_hash64(const void* data, size_t bytes, unsigned long long hash)
{
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < bytes; ++i) { hash ^= p[i]; hash *= 1099511628211ULL; }
    return hash;
}

unsigned long long checkpoint_configuration_hash()
{
    return checkpoint_configuration_hash(false, true, true);
}

unsigned long long checkpoint_configuration_hash(bool low_order_only,
                                                 bool high_order_enabled,
                                                 bool fct_enabled)
{
    unsigned long long h = 1469598103934665603ULL;
    // This is deliberately a physics/grid identity, not a build identity.
    // A restart must never silently mix a different plasma, beam, time, or
    // cylindrical velocity-grid definition with persisted state.
    const double values[] = {
        Const::qe, Const::me, Const::c, Const::eps0,
        Param::dens, Param::temperature_e, Param::temperature_i,
        Param::densb, Param::jb, Param::gambetab, Param::beam_macro_weight,
        Param::Lx, Param::dx, Param::t_inject_start, Param::t_inject_end,
        Param::dt_multiplier, Param::dt_snapshot, Param::velocity_space_cfl,
        Param::semi_lagrangian_cfl, Param::momentum_umax,
        Param::momentum_upar_stretch, Param::momentum_uperp_stretch,
        Param::v_floor, Param::u_floor
    };
    h = checkpoint_hash64(values, sizeof(values), h);
    const char topology[] = "background_periodic_x;field_periodic_face;beam_open_x";
    h = checkpoint_hash64(topology, sizeof(topology) - 1, h);
    const int dims[] = {Param::nx, Param::Nghost, Param::Nv, Param::Nmu,
                        Param::Z_ion, Param::beam_macro_particles_per_cell,
                        Param::enable_beam_boundary_injection ? 1 : 0,
                        Param::abort_on_vmax_loss ? 1 : 0,
                        low_order_only ? 1 : 0, high_order_enabled ? 1 : 0,
                        fct_enabled ? 1 : 0};
    return checkpoint_hash64(dims, sizeof(dims), h);
}

unsigned long long checkpoint_velocity_grid_hash(const Species& bkg)
{
    unsigned long long h = 1469598103934665603ULL;
    const CylindricalVelocityGrid& g = bkg.cgrid;
    const std::vector<double>* arrays[] = {
        &g.upar_faces, &g.upar_cells, &g.upar_widths,
        &g.uperp_faces, &g.uperp_cells, &g.uperp_widths,
        &g.uperp_ring_areas, &g.kinetic_energy, &g.vx
    };
    for (size_t i = 0; i < sizeof(arrays) / sizeof(arrays[0]); ++i) {
        const unsigned long long n = arrays[i]->size();
        h = checkpoint_hash64(&n, sizeof(n), h);
        if (n) h = checkpoint_hash64(arrays[i]->data(), n * sizeof(double), h);
    }
    return h;
}

unsigned long long checkpoint_physics_parameter_hash()
{
    return checkpoint_configuration_hash();
}

CheckpointStateHashes checkpoint_state_hashes(const Species& bkg, const BeamPIC& beam,
                                              const EMFields& fields,
                                              const SpatialGrid& sg,
                                              int rank, int size)
{
    std::vector<double> owned_f(static_cast<size_t>(sg.nx_local) * Param::Nvmu);
    for (int ix = 0; ix < sg.nx_local; ++ix) {
        std::copy(bkg.f.begin() + static_cast<size_t>(sg.nghost + ix) * Param::Nvmu,
                  bkg.f.begin() + static_cast<size_t>(sg.nghost + ix + 1) * Param::Nvmu,
                  owned_f.begin() + static_cast<size_t>(ix) * Param::Nvmu);
    }
    const std::vector<double> owned_faces(fields.Ex_face.begin(),
                                          fields.Ex_face.begin() + sg.nx_local);
    const BeamPersistentState persistent = beam.export_persistent_state();
    unsigned long long local[3] = {
        checkpoint_hash64(owned_f.data(), owned_f.size() * sizeof(double)),
        checkpoint_hash64(owned_faces.data(), owned_faces.size() * sizeof(double)),
        checkpoint_hash64(&persistent, sizeof(persistent))
    };
    if (!beam.particles.empty()) {
        local[2] = checkpoint_hash64(beam.particles.data(),
                                     beam.particles.size() * sizeof(BeamParticle),
                                     local[2]);
    }
    std::vector<unsigned long long> gathered(rank == 0 ? static_cast<size_t>(3 * size) : 0);
    MPI_Gather(local, 3, MPI_UNSIGNED_LONG_LONG, rank == 0 ? gathered.data() : 0,
               3, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    CheckpointStateHashes result = {1469598103934665603ULL,
                                    1469598103934665603ULL,
                                    1469598103934665603ULL};
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            result.background = checkpoint_hash64(&gathered[3 * r], sizeof(unsigned long long),
                                                  result.background);
            result.field_faces = checkpoint_hash64(&gathered[3 * r + 1], sizeof(unsigned long long),
                                                   result.field_faces);
            result.beam = checkpoint_hash64(&gathered[3 * r + 2], sizeof(unsigned long long),
                                            result.beam);
        }
    }
    unsigned long long packed[3] = {result.background, result.field_faces, result.beam};
    MPI_Bcast(packed, 3, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    result.background = packed[0]; result.field_faces = packed[1]; result.beam = packed[2];
    return result;
}

bool read_checkpoint_reference_hashes(const std::string& directory,
                                      CheckpointStateHashes& hashes,
                                      int rank, int size)
{
    int ok = 1;
    unsigned long long packed[3] = {0ULL, 0ULL, 0ULL};
    if (rank == 0) {
        std::ifstream manifest((directory + "/manifest.dat").c_str());
        std::string key;
        while (manifest >> key) {
            if (key == "state_hash_background") manifest >> packed[0];
            else if (key == "state_hash_field_faces") manifest >> packed[1];
            else if (key == "state_hash_beam") manifest >> packed[2];
            else { std::string ignored; std::getline(manifest, ignored); }
        }
        if (!manifest.eof() || packed[0] == 0ULL || packed[1] == 0ULL || packed[2] == 0ULL) ok = 0;
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(packed, 3, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    hashes.background = packed[0]; hashes.field_faces = packed[1]; hashes.beam = packed[2];
    (void)size;
    return ok != 0;
}

bool write_checkpoint(const std::string& directory, const CheckpointControlState& c,
                      const Species& bkg, const BeamPIC& beam, const EMFields& fields,
                      const SpatialGrid& sg, int rank, int size, std::string& error,
                      bool low_order_only, bool high_order_enabled, bool fct_enabled)
{
    const std::string::size_type slash = directory.find_last_of('/');
    const std::string parent = slash == std::string::npos ? std::string() : directory.substr(0, slash);
    int local_ok = (parent.empty() || make_dir(parent)) && make_dir(directory) ? 1 : 0;
    std::vector<double> f_owned(static_cast<size_t>(sg.nx_local) * Param::Nvmu);
    for (int ix=0; ix<sg.nx_local; ++ix) std::copy(bkg.f.begin()+static_cast<size_t>(sg.nghost+ix)*Param::Nvmu,
        bkg.f.begin()+static_cast<size_t>(sg.nghost+ix+1)*Param::Nvmu,
        f_owned.begin()+static_cast<size_t>(ix)*Param::Nvmu);
    std::vector<double> ex_owned(fields.Ex_face.begin(), fields.Ex_face.begin()+sg.nx_local);
    if (!finite_vector(f_owned) || !finite_vector(ex_owned) ||
        !finite_particles(beam.particles)) local_ok = 0;
    const std::string file = rank_path(directory, rank), tmp = file + ".tmp";
    unsigned long long bytes = 0, hash = 1469598103934665603ULL;
    if (local_ok) {
        RankCheckpointHeader h; std::memset(&h, 0, sizeof(h)); std::memcpy(h.magic,kMagic,sizeof(kMagic));
        h.version=kVersion; h.endian_tag=kEndian; h.mpi_rank=rank; h.mpi_size=size; h.ix_start=sg.ix_start;
        h.nx_local=sg.nx_local; h.nx_total=sg.nx_total; h.nghost=sg.nghost; h.nv=Param::Nv; h.nmu=Param::Nmu;
        h.step=c.step; h.time_s=c.time_s; h.dt_s=c.dt_s; h.f_count=f_owned.size(); h.ex_face_count=ex_owned.size(); h.beam_particle_count=beam.particles.size();
        std::ofstream out(tmp.c_str(), std::ios::binary|std::ios::trunc);
        if (!out) local_ok=0; else {
            const BeamPersistentState state=beam.export_persistent_state();
            if (!finite_persistent_state(state)) local_ok = 0;
            if (local_ok) {
                out.write(reinterpret_cast<const char*>(&h),sizeof(h)); out.write(reinterpret_cast<const char*>(&c),sizeof(c));
                out.write(reinterpret_cast<const char*>(&f_owned[0]),static_cast<std::streamsize>(f_owned.size()*sizeof(double)));
                out.write(reinterpret_cast<const char*>(&ex_owned[0]),static_cast<std::streamsize>(ex_owned.size()*sizeof(double)));
                out.write(reinterpret_cast<const char*>(&state),sizeof(state));
                if (!beam.particles.empty()) out.write(reinterpret_cast<const char*>(&beam.particles[0]),static_cast<std::streamsize>(beam.particles.size()*sizeof(BeamParticle)));
                if (!out) local_ok = 0;
            }
            out.close();
            if (local_ok) {
                std::ifstream in(tmp.c_str(),std::ios::binary); char block[4096]; while (in) { in.read(block,sizeof(block)); hash=checkpoint_hash64(block,static_cast<size_t>(in.gcount()),hash); bytes+=static_cast<unsigned long long>(in.gcount()); }
                if (!in.eof()) local_ok = 0;
                if (!rename_atomic(tmp,file)) local_ok=0;
            }
        }
    }
    int global_ok=0; MPI_Allreduce(&local_ok,&global_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    std::vector<unsigned long long> all_bytes(rank==0?static_cast<size_t>(size):0), all_hash(rank==0?static_cast<size_t>(size):0);
    MPI_Gather(&bytes,1,MPI_UNSIGNED_LONG_LONG,rank==0?&all_bytes[0]:0,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
    MPI_Gather(&hash,1,MPI_UNSIGNED_LONG_LONG,rank==0?&all_hash[0]:0,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
    if (rank == 0 && global_ok) {
        for (int r = 0; r < size; ++r) {
            std::ifstream verify(rank_path(directory, r).c_str(), std::ios::binary);
            verify.seekg(0, std::ios::end);
            if (!verify || verify.tellg() < 0 ||
                static_cast<unsigned long long>(verify.tellg()) != all_bytes[r]) {
                global_ok = 0;
                break;
            }
        }
    }
    MPI_Bcast(&global_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    const CheckpointStateHashes state_hashes =
        checkpoint_state_hashes(bkg, beam, fields, sg, rank, size);
    if (rank==0 && global_ok) {
        const std::string manifest=directory+"/manifest.dat", tmp_manifest=manifest+".tmp"; std::ofstream out(tmp_manifest.c_str());
        out<<"magic FPCHKPT1\nversion 1\nendianness "<<kEndian<<"\nsizeof_double "<<sizeof(double)<<"\nsizeof_uint64 "<<sizeof(unsigned long long)<<"\nmpi_size "<<size<<"\nnx "<<Param::nx<<"\nnghost "<<Param::Nghost<<"\nnv "<<Param::Nv<<"\nnmu "<<Param::Nmu
           <<"\nLx "<<std::setprecision(17)<<Param::Lx<<"\ndx "<<Param::dx<<"\nconfig_hash "<<checkpoint_configuration_hash(low_order_only, high_order_enabled, fct_enabled)
           <<"\nvelocity_grid_hash "<<checkpoint_velocity_grid_hash(bkg)
           <<"\nphysics_parameter_hash "<<checkpoint_physics_parameter_hash()
           <<"\nstate_hash_background "<<state_hashes.background
           <<"\nstate_hash_field_faces "<<state_hashes.field_faces
           <<"\nstate_hash_beam "<<state_hashes.beam
           <<"\nstep "<<c.step<<"\ntime_s "<<c.time_s<<"\ndt_s "<<c.dt_s<<"\nrank_files "<<size<<"\n";
        for (int r = 0; r < size; ++r) {
            out << "rank " << r << " " << all_bytes[r] << " "
                << all_hash[r] << "\n";
        }
        out.close();
        if (!out || !rename_atomic(tmp_manifest,manifest)) global_ok=0;
    }
    MPI_Bcast(&global_ok,1,MPI_INT,0,MPI_COMM_WORLD);
    if (!global_ok) error="checkpoint write failed";
    return global_ok != 0;
}

bool read_checkpoint(const std::string& directory, CheckpointControlState& c,
                     Species& bkg, BeamPIC& beam, EMFields& fields,
                     const SpatialGrid& sg, int rank, int size, std::string& error,
                     bool low_order_only, bool high_order_enabled, bool fct_enabled)
{
    enum CheckpointReadFailure {
        kManifestOpenOrMagic = 1 << 0,
        kMpiSizeMismatch = 1 << 1,
        kGridMismatch = 1 << 2,
        kConfigurationHashMismatch = 1 << 3,
        kVelocityGridHashMismatch = 1 << 4,
        kPhysicsHashMismatch = 1 << 5,
        kRankTableMissing = 1 << 6,
        kRankFileOpen = 1 << 7,
        kRankFileSize = 1 << 8,
        kRankFileChecksum = 1 << 9,
        kRankHeaderMismatch = 1 << 10,
        kControlStateRead = 1 << 11,
        kPayloadCountMismatch = 1 << 12,
        kPayloadReadOrNonfinite = 1 << 13
    };
    int manifest_failure_mask = 0;
    unsigned long long expected_bytes = 0, expected_hash = 0;
    if (rank == 0) {
        std::ifstream manifest((directory + "/manifest.dat").c_str());
        std::string key, magic;
        int manifest_size = -1, nx = -1, nv = -1, nmu = -1;
        unsigned long long config_hash = 0, velocity_hash = 0,
                           physics_hash = 0;
        if (!(manifest >> key >> magic) || key != "magic" || magic != "FPCHKPT1")
            manifest_failure_mask |= kManifestOpenOrMagic;
        while (manifest >> key) {
            if (key == "mpi_size") manifest >> manifest_size;
            else if (key == "nx") manifest >> nx;
            else if (key == "nv") manifest >> nv;
            else if (key == "nmu") manifest >> nmu;
            else if (key == "config_hash") manifest >> config_hash;
            else if (key == "velocity_grid_hash") manifest >> velocity_hash;
            else if (key == "physics_parameter_hash") manifest >> physics_hash;
            else if (key == "rank") {
                int entry_rank = -1; unsigned long long bytes = 0, hash = 0;
                manifest >> entry_rank >> bytes >> hash;
                if (entry_rank == rank) { expected_bytes = bytes; expected_hash = hash; }
            } else { std::string ignored; std::getline(manifest, ignored); }
        }
        const unsigned long long expected_config_hash = checkpoint_configuration_hash(
            low_order_only, high_order_enabled, fct_enabled);
        const unsigned long long expected_velocity_hash =
            checkpoint_velocity_grid_hash(bkg);
        const unsigned long long expected_physics_hash =
            checkpoint_physics_parameter_hash();
        if (manifest_size != size) manifest_failure_mask |= kMpiSizeMismatch;
        if (nx != Param::nx || nv != Param::Nv || nmu != Param::Nmu)
            manifest_failure_mask |= kGridMismatch;
        if (config_hash != expected_config_hash)
            manifest_failure_mask |= kConfigurationHashMismatch;
        if (velocity_hash != expected_velocity_hash)
            manifest_failure_mask |= kVelocityGridHashMismatch;
        if (physics_hash != expected_physics_hash)
            manifest_failure_mask |= kPhysicsHashMismatch;
        if (manifest_failure_mask != 0) {
            std::fprintf(stderr,
                "Checkpoint manifest mismatch: path=%s mask=0x%x "
                "mpi_size(saved=%d expected=%d) grid(saved=%d,%d,%d "
                "expected=%d,%d,%d) config_hash(saved=%llu expected=%llu) "
                "velocity_grid_hash(saved=%llu expected=%llu) "
                "physics_parameter_hash(saved=%llu expected=%llu)\n",
                directory.c_str(), manifest_failure_mask, manifest_size, size,
                nx, nv, nmu, Param::nx, Param::Nv, Param::Nmu,
                config_hash, expected_config_hash, velocity_hash,
                expected_velocity_hash, physics_hash, expected_physics_hash);
        }
    }
    MPI_Bcast(&manifest_failure_mask, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // The manifest lists every rank. Broadcast the complete compact table so
    // each rank can verify its own file before entering the solver.
    std::vector<unsigned long long> expected_table;
    if (rank == 0) {
        expected_table.assign(static_cast<size_t>(2 * size), 0ULL);
        std::ifstream manifest((directory + "/manifest.dat").c_str());
        std::string key;
        while (manifest >> key) {
            if (key == "rank") { int r; unsigned long long b, hsh; manifest >> r >> b >> hsh; if (r >= 0 && r < size) { expected_table[2*r]=b; expected_table[2*r+1]=hsh; } }
            else { std::string ignored; std::getline(manifest, ignored); }
        }
    } else expected_table.resize(static_cast<size_t>(2 * size));
    MPI_Bcast(&expected_table[0], 2 * size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    expected_bytes = expected_table[2 * rank]; expected_hash = expected_table[2 * rank + 1];
    int local_failure_mask = manifest_failure_mask;
    if (expected_bytes == 0 || expected_hash == 0)
        local_failure_mask |= kRankTableMissing;
    std::ifstream verify(rank_path(directory, rank).c_str(), std::ios::binary);
    unsigned long long actual_bytes = 0, actual_hash = 1469598103934665603ULL;
    char block[4096]; while (verify) { verify.read(block, sizeof(block)); actual_hash = checkpoint_hash64(block, static_cast<size_t>(verify.gcount()), actual_hash); actual_bytes += static_cast<unsigned long long>(verify.gcount()); }
    if (!verify.eof()) local_failure_mask |= kRankFileOpen;
    if (actual_bytes != expected_bytes) local_failure_mask |= kRankFileSize;
    if (actual_hash != expected_hash) local_failure_mask |= kRankFileChecksum;
    std::ifstream in(rank_path(directory,rank).c_str(),std::ios::binary); RankCheckpointHeader h;
    std::memset(&h, 0, sizeof(h));
    if (!in || !in.read(reinterpret_cast<char*>(&h),sizeof(h)) ||
        std::memcmp(h.magic,kMagic,sizeof(kMagic)) || h.version!=kVersion ||
        h.endian_tag!=kEndian || h.mpi_rank!=rank || h.mpi_size!=size ||
        h.ix_start!=sg.ix_start || h.nx_local!=sg.nx_local ||
        h.nghost!=sg.nghost || h.nv!=Param::Nv || h.nmu!=Param::Nmu)
        local_failure_mask |= kRankHeaderMismatch;
    if (local_failure_mask == 0 &&
        !in.read(reinterpret_cast<char*>(&c),sizeof(c)))
        local_failure_mask |= kControlStateRead;
    std::vector<double> f_owned;
    std::vector<double> ex_owned;
    BeamPersistentState state;
    if (local_failure_mask == 0) {
        if (h.f_count!=static_cast<unsigned long long>(sg.nx_local*Param::Nvmu) ||
            h.ex_face_count!=static_cast<unsigned long long>(sg.nx_local)) {
            local_failure_mask |= kPayloadCountMismatch;
        } else {
            f_owned.resize(static_cast<size_t>(h.f_count));
            ex_owned.resize(static_cast<size_t>(h.ex_face_count));
        }
    }
    if (local_failure_mask == 0) {
        in.read(reinterpret_cast<char*>(&f_owned[0]),static_cast<std::streamsize>(f_owned.size()*sizeof(double)));
        in.read(reinterpret_cast<char*>(&ex_owned[0]),static_cast<std::streamsize>(ex_owned.size()*sizeof(double)));
        in.read(reinterpret_cast<char*>(&state),sizeof(state));
        beam.particles.resize(static_cast<size_t>(h.beam_particle_count));
        if(!beam.particles.empty()) in.read(reinterpret_cast<char*>(&beam.particles[0]),static_cast<std::streamsize>(beam.particles.size()*sizeof(BeamParticle)));
        if(!in || !finite_vector(f_owned)||!finite_vector(ex_owned)||
           !finite_persistent_state(state)||!finite_particles(beam.particles))
            local_failure_mask |= kPayloadReadOrNonfinite;
    }
    std::vector<int> failure_masks(rank == 0 ? static_cast<size_t>(size) : 0);
    unsigned long long local_file_details[4] = {
        expected_bytes, actual_bytes, expected_hash, actual_hash
    };
    std::vector<unsigned long long> all_file_details(
        rank == 0 ? static_cast<size_t>(4 * size) : 0);
    MPI_Gather(&local_failure_mask, 1, MPI_INT,
               rank == 0 ? failure_masks.data() : 0, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    MPI_Gather(local_file_details, 4, MPI_UNSIGNED_LONG_LONG,
               rank == 0 ? all_file_details.data() : 0, 4,
               MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    int local_ok = local_failure_mask == 0 ? 1 : 0;
    int global_ok=0;
    MPI_Allreduce(&local_ok,&global_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(!global_ok){
        if (rank == 0) {
            int printed = 0;
            for (int r = 0; r < size && printed < 16; ++r) {
                if (failure_masks[static_cast<size_t>(r)] == 0) continue;
                std::fprintf(stderr,
                    "Checkpoint rank validation failure: rank=%d mask=0x%x "
                    "bytes(expected=%llu actual=%llu) "
                    "checksum(expected=%llu actual=%llu)\n",
                    r, failure_masks[static_cast<size_t>(r)],
                    all_file_details[static_cast<size_t>(4*r)],
                    all_file_details[static_cast<size_t>(4*r+1)],
                    all_file_details[static_cast<size_t>(4*r+2)],
                    all_file_details[static_cast<size_t>(4*r+3)]);
                ++printed;
            }
            if (printed == 16)
                std::fprintf(stderr, "Checkpoint validation: additional failing ranks omitted\n");
            std::fprintf(stderr,
                "Checkpoint failure mask legend: 0x1 manifest/magic, 0x2 mpi_size, "
                "0x4 grid, 0x8 config_hash, 0x10 velocity_grid_hash, "
                "0x20 physics_hash, 0x40 rank_table, 0x80 file_open, "
                "0x100 file_size, 0x200 checksum, 0x400 header, "
                "0x800 control, 0x1000 counts, 0x2000 payload/nonfinite\n");
        }
        error="checkpoint validation failed; see rank-0 diagnostics";
        return false;
    }
    for(int ix=0;ix<sg.nx_local;++ix) std::copy(f_owned.begin()+static_cast<size_t>(ix)*Param::Nvmu,f_owned.begin()+static_cast<size_t>(ix+1)*Param::Nvmu,bkg.f.begin()+static_cast<size_t>(sg.nghost+ix)*Param::Nvmu);
    std::copy(ex_owned.begin(),ex_owned.end(),fields.Ex_face.begin()); beam.import_persistent_state(state,sg); fields.sync_cell_ex_from_faces(rank,size);
    return true;
}

bool write_midpoint_audit_state(
    const std::string& directory,
    const VlasovAmpereMidpointSolver::MidpointAuditState& state,
    const SpatialGrid& sg, int rank, int size, std::string& error)
{
    const std::string::size_type slash = directory.find_last_of('/');
    const std::string parent = slash == std::string::npos ? std::string() : directory.substr(0, slash);
    int ok = (parent.empty() || make_dir(parent)) && make_dir(directory) ? 1 : 0;
    const std::string file = rank_path(directory, rank), tmp = file + ".tmp";
    unsigned long long bytes = 0ULL, hash = 1469598103934665603ULL;
    if (ok) {
        std::ofstream out(tmp.c_str(), std::ios::binary | std::ios::trunc);
        const unsigned long long count = static_cast<unsigned long long>(sg.nx_local) * Param::Nvmu;
        const unsigned long long face_count = static_cast<unsigned long long>(sg.nx_local);
        const int type_length = static_cast<int>(state.acceptance_type.size());
        out.write("FPMIDPT6", 8); out.write(reinterpret_cast<const char*>(&state.step), sizeof(state.step));
        out.write(reinterpret_cast<const char*>(&state.time_s), sizeof(state.time_s)); out.write(reinterpret_cast<const char*>(&state.dt_s), sizeof(state.dt_s));
        out.write(reinterpret_cast<const char*>(&state.substeps_used), sizeof(state.substeps_used)); out.write(reinterpret_cast<const char*>(&state.nonlinear_iterations), sizeof(state.nonlinear_iterations));
        const unsigned char switches[3] = { static_cast<unsigned char>(state.low_order_only), static_cast<unsigned char>(state.high_order_enabled), static_cast<unsigned char>(state.fct_enabled) };
        out.write(reinterpret_cast<const char*>(switches), sizeof(switches));
        out.write(reinterpret_cast<const char*>(&state.background_coupling_mode),
                  sizeof(state.background_coupling_mode));
        out.write(reinterpret_cast<const char*>(&type_length), sizeof(type_length)); out.write(state.acceptance_type.data(), type_length);
        out.write(reinterpret_cast<const char*>(&count), sizeof(count)); out.write(reinterpret_cast<const char*>(&face_count), sizeof(face_count));
        for (int ix=0; ix<sg.nx_local; ++ix) {
            const size_t base=static_cast<size_t>(sg.nghost+ix)*Param::Nvmu;
            out.write(reinterpret_cast<const char*>(&state.bkg_n.f[base]), Param::Nvmu*sizeof(double));
            out.write(reinterpret_cast<const char*>(&state.guess_np1.f[base]), Param::Nvmu*sizeof(double));
            out.write(reinterpret_cast<const char*>(&state.operator_input_guess.f[base]), Param::Nvmu*sizeof(double));
        }
        out.write(reinterpret_cast<const char*>(&state.fields_n.Ex_face[0]), static_cast<std::streamsize>(face_count*sizeof(double)));
        out.write(reinterpret_cast<const char*>(&state.fields_end_guess.Ex_face[0]), static_cast<std::streamsize>(face_count*sizeof(double)));
        out.write(reinterpret_cast<const char*>(&state.fields_np1.Ex_face[0]),
                  static_cast<std::streamsize>(face_count * sizeof(double)));
        out.write(reinterpret_cast<const char*>(&state.reference_stage5_r_fv),
                  sizeof(state.reference_stage5_r_fv));
        out.write(reinterpret_cast<const char*>(&state.reference_stage5_r_couple),
                  sizeof(state.reference_stage5_r_couple));
        const double limiter_values[6] = {
            state.limiter_active_fraction,
            state.limiter_active_fraction_core,
            state.limiter_active_fraction_boundary,
            state.x_limiter_active_fraction,
            state.u_limiter_active_fraction,
            state.limiter_min_alpha
        };
        out.write(reinterpret_cast<const char*>(limiter_values),
                  sizeof(limiter_values));
        out.write(reinterpret_cast<const char*>(&state.coupling_layout.beam_front_ix),
                  sizeof(state.coupling_layout.beam_front_ix));
        out.write(reinterpret_cast<const char*>(&state.coupling_layout.wave_core_end_m),
                  sizeof(state.coupling_layout.wave_core_end_m));
        out.write(reinterpret_cast<const char*>(state.periodic_seam_face_audit.data()),
                  static_cast<std::streamsize>(state.periodic_seam_face_audit.size() *
                                               sizeof(double)));
        const std::vector<double>* arrays[] = {&state.j_beam_face_mid, &state.reference_jn_face, &state.reference_je_cell, &state.reference_gstar_je_face};
        for (int i=0;i<4;++i) { const unsigned long long n=arrays[i]->size(); out.write(reinterpret_cast<const char*>(&n),sizeof(n)); if(n) out.write(reinterpret_cast<const char*>(&(*arrays[i])[0]),static_cast<std::streamsize>(n*sizeof(double))); }
        out.close();
        if (!out) ok = 0;
        if (ok) {
            std::ifstream verify(tmp.c_str(), std::ios::binary);
            char block[4096];
            while (verify) {
                verify.read(block, sizeof(block));
                const std::streamsize count_read = verify.gcount();
                hash = checkpoint_hash64(block, static_cast<size_t>(count_read), hash);
                bytes += static_cast<unsigned long long>(count_read);
            }
            if (!rename_atomic(tmp,file)) ok=0;
        }
    }
    int all_ok=0; MPI_Allreduce(&ok,&all_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    std::vector<unsigned long long> all_bytes(rank == 0 ? static_cast<size_t>(size) : 0),
        all_hash(rank == 0 ? static_cast<size_t>(size) : 0);
    MPI_Gather(&bytes, 1, MPI_UNSIGNED_LONG_LONG, rank == 0 ? &all_bytes[0] : 0,
               1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&hash, 1, MPI_UNSIGNED_LONG_LONG, rank == 0 ? &all_hash[0] : 0,
               1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    if(rank==0 && all_ok) {
        std::ofstream manifest((directory+"/manifest.dat.tmp").c_str());
        manifest<<"magic FPMIDPT6\nmpi_size "<<size<<"\nnx "<<Param::nx<<"\nnv "<<Param::Nv<<"\nnmu "<<Param::Nmu<<"\nconfig_hash "<<checkpoint_configuration_hash()<<"\nvelocity_grid_hash "<<checkpoint_velocity_grid_hash(state.bkg_n)<<"\nphysics_parameter_hash "<<checkpoint_physics_parameter_hash()<<"\nstep "<<state.step<<"\ntime_s "<<std::setprecision(17)<<state.time_s<<"\ndt_s "<<state.dt_s<<"\nbackground_coupling_mode "<<state.background_coupling_mode<<"\nacceptance "<<state.acceptance_type<<"\n";
        for (int r = 0; r < size; ++r) manifest << "rank " << r << " "
            << all_bytes[r] << " " << all_hash[r] << "\n";
        manifest.close();
        if(!manifest || !rename_atomic(directory+"/manifest.dat.tmp",directory+"/manifest.dat")) all_ok=0;
    }
    MPI_Bcast(&all_ok,1,MPI_INT,0,MPI_COMM_WORLD); if(!all_ok) error="midpoint audit write failed"; return all_ok!=0;
}

bool read_midpoint_audit_state(
    const std::string& directory,
    VlasovAmpereMidpointSolver::MidpointAuditState& state,
    const Species& species_template, const EMFields& fields_template,
    const SpatialGrid& sg, int rank, int size, std::string& error)
{
    int ok = 1;
    std::string manifest_error;
    std::vector<unsigned long long> audit_file_table;
    unsigned long long expected_velocity_hash = 0ULL;
    if (rank == 0) {
        std::ifstream manifest((directory + "/manifest.dat").c_str());
        std::string key, magic;
        int saved_size = -1, saved_nx = -1, saved_nv = -1, saved_nmu = -1;
        unsigned long long saved_hash = 0ULL, saved_velocity_hash = 0ULL,
                           saved_physics_hash = 0ULL;
        if (!(manifest >> key >> magic) || key != "magic" || magic != "FPMIDPT6") {
            ok = 0;
            manifest_error = "midpoint audit manifest is missing or has invalid magic";
        }
        while (manifest >> key) {
            if (key == "mpi_size") manifest >> saved_size;
            else if (key == "nx") manifest >> saved_nx;
            else if (key == "nv") manifest >> saved_nv;
            else if (key == "nmu") manifest >> saved_nmu;
            else if (key == "config_hash") manifest >> saved_hash;
            else if (key == "velocity_grid_hash") manifest >> saved_velocity_hash;
            else if (key == "physics_parameter_hash") manifest >> saved_physics_hash;
            else if (key == "rank") {
                int file_rank = -1; unsigned long long bytes = 0ULL, hash = 0ULL;
                manifest >> file_rank >> bytes >> hash;
                if (audit_file_table.empty()) audit_file_table.assign(static_cast<size_t>(2 * size), 0ULL);
                if (file_rank >= 0 && file_rank < size) {
                    audit_file_table[2 * file_rank] = bytes;
                    audit_file_table[2 * file_rank + 1] = hash;
                }
            }
            else { std::string ignored; std::getline(manifest, ignored); }
        }
        const unsigned long long expected_config_hash =
            checkpoint_configuration_hash();
        const unsigned long long expected_physics_hash =
            checkpoint_physics_parameter_hash();
        const bool size_mismatch = saved_size != size;
        const bool grid_mismatch = saved_nx != Param::nx ||
            saved_nv != Param::Nv || saved_nmu != Param::Nmu;
        const bool config_mismatch = saved_hash != expected_config_hash;
        const bool physics_mismatch = saved_physics_hash != expected_physics_hash;
        const bool table_mismatch = audit_file_table.size() !=
            static_cast<size_t>(2 * size);
        const bool velocity_hash_missing = saved_velocity_hash == 0ULL;
        if (size_mismatch || grid_mismatch || config_mismatch ||
            physics_mismatch || table_mismatch || velocity_hash_missing) {
            ok = 0;
            std::ostringstream detail;
            detail << "midpoint audit manifest is incompatible:"
                   << " mpi_size(saved=" << saved_size
                   << " expected=" << size << ")"
                   << " grid(saved=" << saved_nx << "," << saved_nv << ","
                   << saved_nmu << " expected=" << Param::nx << ","
                   << Param::Nv << "," << Param::Nmu << ")"
                   << " config_hash(saved=" << saved_hash
                   << " expected=" << expected_config_hash << ")"
                   << " physics_hash(saved=" << saved_physics_hash
                   << " expected=" << expected_physics_hash << ")"
                   << " velocity_hash=" << saved_velocity_hash
                   << " rank_table_entries=" << audit_file_table.size() / 2
                   << " expected_rank_entries=" << size;
            manifest_error = detail.str();
        }
        expected_velocity_hash = saved_velocity_hash;
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!ok) {
        error = (rank == 0 && !manifest_error.empty()) ? manifest_error :
            "midpoint audit manifest is incompatible on rank 0";
        return false;
    }
    if (rank != 0) audit_file_table.resize(static_cast<size_t>(2 * size), 0ULL);
    MPI_Bcast(&audit_file_table[0], 2 * size, MPI_UNSIGNED_LONG_LONG, 0,
              MPI_COMM_WORLD);
    std::ifstream verify(rank_path(directory, rank).c_str(), std::ios::binary);
    unsigned long long actual_bytes = 0ULL, actual_hash = 1469598103934665603ULL;
    char verify_block[4096];
    while (verify) {
        verify.read(verify_block, sizeof(verify_block));
        const std::streamsize count_read = verify.gcount();
        actual_hash = checkpoint_hash64(verify_block, static_cast<size_t>(count_read), actual_hash);
        actual_bytes += static_cast<unsigned long long>(count_read);
    }
    if (audit_file_table[2 * rank] == 0ULL ||
        actual_bytes != audit_file_table[2 * rank] ||
        actual_hash != audit_file_table[2 * rank + 1]) ok = 0;

    std::ifstream in(rank_path(directory, rank).c_str(), std::ios::binary);
    char magic[8] = {0};
    in.read(magic, sizeof(magic));
    if (!in || std::memcmp(magic, "FPMIDPT6", sizeof(magic)) != 0) ok = 0;
    unsigned char switches[3] = {0, 0, 0};
    int type_length = 0;
    unsigned long long count = 0, face_count = 0;
    if (ok) {
        in.read(reinterpret_cast<char*>(&state.step), sizeof(state.step));
        in.read(reinterpret_cast<char*>(&state.time_s), sizeof(state.time_s));
        in.read(reinterpret_cast<char*>(&state.dt_s), sizeof(state.dt_s));
        in.read(reinterpret_cast<char*>(&state.substeps_used), sizeof(state.substeps_used));
        in.read(reinterpret_cast<char*>(&state.nonlinear_iterations), sizeof(state.nonlinear_iterations));
        in.read(reinterpret_cast<char*>(switches), sizeof(switches));
        in.read(reinterpret_cast<char*>(&state.background_coupling_mode),
                sizeof(state.background_coupling_mode));
        in.read(reinterpret_cast<char*>(&type_length), sizeof(type_length));
        if (!in || type_length < 0 || type_length > 128 ||
            state.background_coupling_mode <
                VlasovAmpereMidpointSolver::LEGACY_COUPLING ||
            state.background_coupling_mode >
                VlasovAmpereMidpointSolver::DUAL_U_COUPLING) ok = 0;
    }
    if (ok) {
        state.acceptance_type.assign(static_cast<size_t>(type_length), '\0');
        if (type_length) in.read(&state.acceptance_type[0], type_length);
        in.read(reinterpret_cast<char*>(&count), sizeof(count));
        in.read(reinterpret_cast<char*>(&face_count), sizeof(face_count));
        if (!in || count != static_cast<unsigned long long>(sg.nx_local) * Param::Nvmu ||
            face_count != static_cast<unsigned long long>(sg.nx_local)) ok = 0;
    }
    state.low_order_only = switches[0] != 0;
    state.high_order_enabled = switches[1] != 0;
    state.fct_enabled = switches[2] != 0;
    state.bkg_n = species_template;
    state.guess_np1 = species_template;
    state.operator_input_guess = species_template;
    state.fields_n = fields_template;
    state.fields_end_guess = fields_template;
    state.fields_np1 = fields_template;
    for (int ix = 0; ok && ix < sg.nx_local; ++ix) {
        const size_t base = static_cast<size_t>(sg.nghost + ix) * Param::Nvmu;
        in.read(reinterpret_cast<char*>(&state.bkg_n.f[base]), Param::Nvmu * sizeof(double));
        in.read(reinterpret_cast<char*>(&state.guess_np1.f[base]), Param::Nvmu * sizeof(double));
        in.read(reinterpret_cast<char*>(&state.operator_input_guess.f[base]), Param::Nvmu * sizeof(double));
    }
    if (ok) {
        in.read(reinterpret_cast<char*>(&state.fields_n.Ex_face[0]),
                static_cast<std::streamsize>(face_count * sizeof(double)));
        in.read(reinterpret_cast<char*>(&state.fields_end_guess.Ex_face[0]),
                static_cast<std::streamsize>(face_count * sizeof(double)));
        in.read(reinterpret_cast<char*>(&state.fields_np1.Ex_face[0]),
                static_cast<std::streamsize>(face_count * sizeof(double)));
        in.read(reinterpret_cast<char*>(&state.reference_stage5_r_fv),
                sizeof(state.reference_stage5_r_fv));
        in.read(reinterpret_cast<char*>(&state.reference_stage5_r_couple),
                sizeof(state.reference_stage5_r_couple));
        double limiter_values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        in.read(reinterpret_cast<char*>(limiter_values), sizeof(limiter_values));
        state.limiter_active_fraction = limiter_values[0];
        state.limiter_active_fraction_core = limiter_values[1];
        state.limiter_active_fraction_boundary = limiter_values[2];
        state.x_limiter_active_fraction = limiter_values[3];
        state.u_limiter_active_fraction = limiter_values[4];
        state.limiter_min_alpha = limiter_values[5];
        in.read(reinterpret_cast<char*>(&state.coupling_layout.beam_front_ix),
                sizeof(state.coupling_layout.beam_front_ix));
        in.read(reinterpret_cast<char*>(&state.coupling_layout.wave_core_end_m),
                sizeof(state.coupling_layout.wave_core_end_m));
        in.read(reinterpret_cast<char*>(state.periodic_seam_face_audit.data()),
                static_cast<std::streamsize>(state.periodic_seam_face_audit.size() *
                                             sizeof(double)));
    }
    std::vector<double>* arrays[] = {&state.j_beam_face_mid, &state.reference_jn_face,
                                     &state.reference_je_cell, &state.reference_gstar_je_face};
    for (int i = 0; ok && i < 4; ++i) {
        unsigned long long n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!in || n > static_cast<unsigned long long>(2 * (sg.nx_local + 1))) { ok = 0; break; }
        arrays[i]->assign(static_cast<size_t>(n), 0.0);
        if (n) in.read(reinterpret_cast<char*>(&(*arrays[i])[0]),
                       static_cast<std::streamsize>(n * sizeof(double)));
        if (!in || !finite_vector(*arrays[i])) ok = 0;
    }
    bool finite_periodic_seam_audit = true;
    for (size_t i = 0; i < state.periodic_seam_face_audit.size(); ++i) {
        finite_periodic_seam_audit = finite_periodic_seam_audit &&
            std::isfinite(state.periodic_seam_face_audit[i]);
    }
    if (ok && (!finite_vector(state.bkg_n.f) || !finite_vector(state.guess_np1.f) ||
               !finite_vector(state.operator_input_guess.f) ||
               !finite_vector(state.fields_n.Ex_face) ||
               !finite_vector(state.fields_end_guess.Ex_face) ||
               !finite_vector(state.fields_np1.Ex_face) ||
                !std::isfinite(state.reference_stage5_r_fv) ||
               !std::isfinite(state.reference_stage5_r_couple) ||
               !std::isfinite(state.limiter_active_fraction) ||
               !std::isfinite(state.limiter_active_fraction_core) ||
               !std::isfinite(state.limiter_active_fraction_boundary) ||
               !std::isfinite(state.x_limiter_active_fraction) ||
               !std::isfinite(state.u_limiter_active_fraction) ||
               !std::isfinite(state.limiter_min_alpha) ||
                !std::isfinite(state.coupling_layout.wave_core_end_m) ||
                !finite_periodic_seam_audit ||
                !std::isfinite(state.time_s) || !std::isfinite(state.dt_s))) ok = 0;
    int global_ok = 0;
    MPI_Allreduce(&ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) { error = "midpoint audit rank file failed validation"; return false; }
    state.fields_n.sync_cell_ex_from_faces(rank, size);
    state.fields_end_guess.sync_cell_ex_from_faces(rank, size);
    state.fields_np1.sync_cell_ex_from_faces(rank, size);
    MPI_Bcast(&expected_velocity_hash, 1, MPI_UNSIGNED_LONG_LONG, 0,
              MPI_COMM_WORLD);
    if (expected_velocity_hash != checkpoint_velocity_grid_hash(state.bkg_n)) {
        error = "midpoint audit velocity grid hash mismatch";
        return false;
    }
    return true;
}
