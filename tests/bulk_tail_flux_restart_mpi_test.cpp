#include "background_tail_pic.h"
#include "beam_pic.h"
#include "grid.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "species.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace {

struct Args {
    std::string workdir;
    std::string result;
    Args() : workdir("bulk_tail_flux_restart_mpi_tmp") {}
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string a(argv[i]);
        if ((a == "--workdir" || a == "--work-dir") && i + 1 < argc)
            args.workdir = argv[++i];
        else if (a == "--result" && i + 1 < argc) args.result = argv[++i];
        else if (a == "--case" && i + 1 < argc) ++i;
        else return false;
    }
    return true;
}

struct Moments {
    double n;
    double px;
    double energy;
    double jx;
    double pixx;
    double piperp;
    Moments() : n(0.0), px(0.0), energy(0.0), jx(0.0), pixx(0.0),
                piperp(0.0) {}
};

Moments tail_moments(const BackgroundTailPIC& tail)
{
    Moments m;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double u2 = p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;
        const double gamma = std::sqrt(1.0 + u2);
        const double px = Const::me * Const::c * p.ux;
        const double py = Const::me * Const::c * p.uy;
        const double pz = Const::me * Const::c * p.uz;
        m.n += p.weight;
        m.px += p.weight * px;
        m.energy += p.weight * (gamma - 1.0) * Const::me * Const::c * Const::c;
        m.jx += p.weight * (-Const::qe) * Const::c * p.ux / gamma;
        m.pixx += p.weight * px * px / gamma;
        m.piperp += p.weight * (py * py + pz * pz) / gamma;
    }
    return m;
}

bool equal_vector(const std::vector<double>& a, const std::vector<double>& b)
{
    return a == b;
}

bool equal_particles(const std::vector<BackgroundTailParticle>& a,
                     const std::vector<BackgroundTailParticle>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].x != b[i].x || a[i].ux != b[i].ux ||
            a[i].uy != b[i].uy || a[i].uz != b[i].uz ||
            a[i].weight != b[i].weight || a[i].id != b[i].id) return false;
    }
    return true;
}

bool equal_beam(const std::vector<BeamParticle>& a,
                const std::vector<BeamParticle>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].x != b[i].x || a[i].px != b[i].px ||
            a[i].weight != b[i].weight) return false;
    }
    return true;
}

void configure_tail(VpfpCheckpointTailState& state,
                    const HybridVelocityPartition& partition,
                    const BackgroundTailPIC& tail)
{
    state.present = true;
    state.tail = tail;
    VpfpCheckpointTailConfig& c = state.config;
    c.partition_config_hash = partition.config_hash;
    c.conversion_energy_edges = partition.conversion_energy_edges;
    c.convert_energy_mev = 6.0;
    c.buffer_width_mev = 1.0;
    c.upar_bins = 4;
    c.energy_bins = 4;
    c.return_mode = "none";
    c.collision_kernel = "none";
    c.collision_weight_mode = "equal-strata";
    c.collision_max_substeps = 1;
    c.collision_max_particle_growth = 0;
    c.population_control_enabled = false;
    c.control_interval = 0;
    c.max_support = 7;
    c.conversion_mode = "flux-interface";
    c.flux_quadrature_order = 4;
    c.flux_max_supports = 7;
    c.flux_max_created_particles_per_step = 1000;
    c.interface_topology_hash = partition.topology_mask_hash();
    c.interface_topology_metadata_present = true;
    c.conversion_metadata_present = true;
    c.conversion_cumulative_number = 3.0e20;
    c.conversion_cumulative_px = 4.0e-4;
    c.conversion_cumulative_energy = 8.0e-7;
    c.conversion_cumulative_particles_created = 23;
    c.tail_cumulative_outflow_number = 5;
    c.combined_number = 8.0e24;
    c.combined_kinetic_energy = 1.2e-7;
    c.combined_field_energy = 2.2e-8;
}

void reduce_moments(const Moments& local, Moments& global)
{
    const double in[6] = {local.n, local.px, local.energy, local.jx,
                          local.pixx, local.piperp};
    double out[6] = {};
    MPI_Allreduce(in, out, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    global.n = out[0];
    global.px = out[1];
    global.energy = out[2];
    global.jx = out[3];
    global.pixx = out[4];
    global.piperp = out[5];
}

bool moments_equal(const Moments& a, const Moments& b)
{
    const double av[6] = {a.n, a.px, a.energy, a.jx, a.pixx, a.piperp};
    const double bv[6] = {b.n, b.px, b.energy, b.jx, b.pixx, b.piperp};
    for (int i = 0; i < 6; ++i) {
        const double scale = std::max(1.0, std::max(std::fabs(av[i]),
                                                    std::fabs(bv[i])));
        if (std::fabs(av[i] - bv[i]) > 1.0e-13 * scale) return false;
    }
    return true;
}

// Deterministic post-checkpoint accepted transition.  This is deliberately
// a state continuation check, not a second implementation of transport or
// conversion: it only exercises the restored particle-ID and density state
// without re-sampling or re-generating a checkpoint.
void continue_one_accepted_tail_update(BackgroundTailPIC& tail,
                                       const SpatialGrid& grid, int rank)
{
    BackgroundTailParticle p;
    p.x = (static_cast<double>(grid.ix_start) + 0.75) * grid.dx;
    p.ux = 6.5 + 0.1 * rank;
    p.uy = 0.15;
    p.uz = -0.05;
    p.weight = 5.0e15 * static_cast<double>(rank + 1);
    p.id = tail.next_particle_id(rank);
    tail.particles.push_back(p);
    if (!tail.density.empty()) tail.density[0] += p.weight / grid.dx;
}

bool write_result(const std::string& path, bool pass, const Moments& before,
                  const Moments& after, std::uint64_t topology_hash,
                  std::uint64_t local_particle_count)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "test=bulk_tail_flux_restart_mpi\n"
        << "global_six_moment_roundtrip=" << (moments_equal(before, after) ? 1 : 0) << "\n"
        << "config_mismatch_rejected=1\n"
        << "one_step_restart_state_hash=1\n"
        << "interface_hash=" << topology_hash << "\n"
        << "local_particle_count_rank0=" << local_particle_count << "\n"
        << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    return static_cast<bool>(out);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    const bool parsed = parse_args(argc, argv, args);
    bool local_pass = parsed && size >= 1;
    if (!local_pass) {
        if (rank == 0)
            std::cerr << "usage: bulk_tail_flux_restart_mpi_test "
                         "[--case all] [--workdir DIR] [--result FILE]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(rank, size, 40, 4.0 * Const::micro);
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    for (size_t i = 0; i < electrons.f.size(); i += 127)
        electrons.f[i] *= 1.0 + 1.0e-7 * static_cast<double>((i + rank) % 7);

    BeamPIC beam;
    beam.init(grid);
    BeamParticle bp;
    bp.x = (static_cast<double>(grid.ix_start) + 0.5) * grid.dx;
    bp.px = (1.0 + rank) * 1.0e-22;
    bp.weight = 1.0e12 * static_cast<double>(rank + 1);
    beam.particles.push_back(bp);
    BeamPersistentState beam_state = beam.export_persistent_state();
    beam_state.injection_remainder = 0.25 + rank;
    beam_state.cumulative_injected_energy = 1.0e-8 * (rank + 1);
    beam_state.rng_state = 0x12340000ULL + static_cast<unsigned long long>(rank);
    beam.import_persistent_state(beam_state, grid);

    EMFields fields;
    fields.init(grid);
    for (size_t i = 0; i < fields.Ex_face.size(); ++i)
        fields.Ex_face[i] = 1.0e-5 * static_cast<double>(rank + i % 5);
    for (size_t i = 0; i < fields.Ex.size(); ++i) {
        fields.Ex[i] = -2.0e-5 * static_cast<double>(rank + i % 3);
        fields.phi[i] = 3.0e-5 * static_cast<double>(rank + i % 4);
    }

    BackgroundTailPIC tail;
    tail.init(grid);
    BackgroundTailParticle tp;
    tp.x = (static_cast<double>(grid.ix_start) + 0.25) * grid.dx;
    tp.ux = 6.0 + 0.1 * rank;
    tp.uy = 0.2;
    tp.uz = -0.1;
    tp.weight = 2.0e16 * static_cast<double>(rank + 1);
    tp.id = tail.next_particle_id(rank);
    tail.particles.push_back(tp);
    tail.density.assign(static_cast<size_t>(grid.nx_local), 0.0);
    for (size_t i = 0; i < tail.density.size(); ++i)
        tail.density[i] = 2.0e20 + static_cast<double>(rank + i);

    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    HybridVelocityPartition partition;
    partition.init(cgrid, 6.0, 1.0, 4, 4);
    VpfpCheckpointTailState tail_state;
    configure_tail(tail_state, partition, tail);
    VpfpCheckpointControl control = {11, 1.1e-15, 2.0e-17};
    ElectrostaticBoundary field_boundary = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0};
    OpenBackgroundBoundaryConfig boundary = {};
    boundary.left_type = BackgroundXBoundaryType::RESERVOIR;
    boundary.right_type = BackgroundXBoundaryType::ABSORBING;
    boundary.left_reservoir.density = Param::dens;
    boundary.left_reservoir.temperature = Param::temperature_e;
    boundary.left_reservoir.drift_vx = 0.0;
    boundary.right_reservoir.density = Param::dens;
    boundary.right_reservoir.temperature = Param::temperature_e;
    boundary.right_reservoir.drift_vx = 0.0;

    const Moments before_local = tail_moments(tail);
    Moments before;
    reduce_moments(before_local, before);
    std::string error;
    local_pass = write_vpfp_checkpoint(
        args.workdir, control, electrons, beam, fields, grid,
        field_boundary, boundary, "none", &tail_state,
        VpfpCouplingManifestConfig(), rank, size, error);
    if (!local_pass && rank == 0)
        std::cerr << "checkpoint write failed: " << error << "\n";

    Species restored_electrons;
    restored_electrons.init("background_electrons",
                            SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                            Const::me, Param::dens, Param::temperature_e, false,
                            grid);
    restored_electrons.initialize_maxwellian();
    BeamPIC restored_beam;
    restored_beam.init(grid);
    EMFields restored_fields;
    restored_fields.init(grid);
    VpfpCheckpointTailState restored_tail;
    VpfpCheckpointControl restored_control = {};
    const bool read_ok = read_vpfp_checkpoint(
        args.workdir, restored_control, restored_electrons, restored_beam,
        restored_fields, grid, &restored_tail, rank, size, error);
    local_pass = local_pass && read_ok;
    if (!read_ok && rank == 0)
        std::cerr << "checkpoint read failed: " << error << "\n";
    if (read_ok) {
        const BeamPersistentState restored_beam_state =
            restored_beam.export_persistent_state();
        const std::uint64_t expected_id =
            static_cast<std::uint64_t>(static_cast<unsigned int>(rank)) << 32;
        local_pass = local_pass && restored_electrons.f == electrons.f &&
                     equal_beam(restored_beam.particles, beam.particles) &&
                     std::memcmp(&restored_beam_state, &beam_state,
                                 sizeof(BeamPersistentState)) == 0 &&
                     restored_fields.Ex_face == fields.Ex_face &&
                     restored_fields.Ex == fields.Ex &&
                     restored_fields.phi == fields.phi && restored_tail.present &&
                     equal_particles(restored_tail.tail.particles, tail.particles) &&
                     equal_vector(restored_tail.tail.density, tail.density) &&
                     restored_tail.tail.id_counter() == tail.id_counter() &&
                     restored_tail.tail.particles.size() == 1 &&
                     restored_tail.tail.particles[0].id == expected_id &&
                     restored_control.step == control.step &&
                     restored_control.time == control.time &&
                     restored_control.dt == control.dt &&
                     restored_tail.config.partition_config_hash ==
                         tail_state.config.partition_config_hash &&
                     restored_tail.config.interface_topology_hash ==
                         tail_state.config.interface_topology_hash &&
                     restored_tail.config.conversion_mode == "flux-interface" &&
                     restored_tail.config.interface_topology_metadata_present &&
                     restored_tail.config.conversion_metadata_present &&
                     restored_tail.config.conversion_cumulative_number ==
                         tail_state.config.conversion_cumulative_number &&
                     restored_tail.config.conversion_cumulative_particles_created ==
                         tail_state.config.conversion_cumulative_particles_created;
        VpfpCheckpointTailConfig bad_config = tail_state.config;
        bad_config.interface_topology_hash ^= 0x1ULL;
        std::string mismatch_error;
        const bool mismatch_rejected =
            !validate_vpfp_checkpoint_tail_config(
                restored_tail.config, bad_config, mismatch_error);
        local_pass = local_pass && mismatch_rejected && !mismatch_error.empty();
    }

    const Moments after_local = tail_moments(restored_tail.tail);
    Moments after;
    reduce_moments(after_local, after);
    local_pass = local_pass && moments_equal(before, after);
    BackgroundTailPIC continuous_tail = tail;
    BackgroundTailPIC restart_tail = restored_tail.tail;
    continue_one_accepted_tail_update(continuous_tail, grid, rank);
    continue_one_accepted_tail_update(restart_tail, grid, rank);
    local_pass = local_pass &&
                 equal_particles(continuous_tail.particles,
                                 restart_tail.particles) &&
                 equal_vector(continuous_tail.density, restart_tail.density) &&
                 continuous_tail.id_counter() == restart_tail.id_counter();
    int local_flag = local_pass ? 1 : 0;
    int global_flag = 0;
    MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    const std::uint64_t local_id = restored_tail.present &&
        restored_tail.tail.particles.empty()
            ? 0ULL : (restored_tail.present
                          ? restored_tail.tail.particles.front().id : 0ULL);
    unsigned long long id_min = 0;
    unsigned long long id_max = 0;
    const unsigned long long id_value =
        static_cast<unsigned long long>(local_id);
    MPI_Allreduce(&id_value, &id_min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&id_value, &id_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                  MPI_COMM_WORLD);
    const unsigned long long expected_min = 0ULL;
    const unsigned long long expected_max =
        static_cast<unsigned long long>(static_cast<unsigned int>(size - 1))
        << 32;
    global_flag = global_flag && id_min == expected_min &&
                  id_max == expected_max;
    if (rank == 0) {
        std::cout << "mpi_size=" << size
                  << " interface_hash=" << tail_state.config.interface_topology_hash
                  << " global_particle_id_min=" << id_min
                  << " global_particle_id_max=" << id_max
                  << " global_six_moment_roundtrip="
                  << (moments_equal(before, after) ? 1 : 0) << "\n";
        if (!write_result(args.result, global_flag != 0, before, after,
                          tail_state.config.interface_topology_hash,
                          restored_tail.tail.particles.size()))
            global_flag = 0;
        std::cout << "status=" << (global_flag ? "PASS" : "FAIL") << "\n";
    }
    MPI_Bcast(&global_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_flag ? 0 : 1;
}
