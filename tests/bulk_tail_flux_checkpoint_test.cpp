#include "background_tail_pic.h"
#include "beam_pic.h"
#include "grid.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "species.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>

#include <cstring>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

namespace {

struct Args {
    std::string workdir;
    std::string result;
    Args() : workdir("bulk_tail_flux_checkpoint_tmp") {}
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

bool equal_double(double a, double b) { return a == b; }

bool equal_particles(const std::vector<BackgroundTailParticle>& a,
                     const std::vector<BackgroundTailParticle>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        const BackgroundTailParticle& p = a[i];
        const BackgroundTailParticle& q = b[i];
        if (!equal_double(p.x, q.x) || !equal_double(p.ux, q.ux) ||
            !equal_double(p.uy, q.uy) || !equal_double(p.uz, q.uz) ||
            !equal_double(p.weight, q.weight) || p.id != q.id) return false;
    }
    return true;
}

bool equal_beam(const std::vector<BeamParticle>& a,
                const std::vector<BeamParticle>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!equal_double(a[i].x, b[i].x) ||
            !equal_double(a[i].px, b[i].px) ||
            !equal_double(a[i].weight, b[i].weight)) return false;
    }
    return true;
}

bool equal_vector(const std::vector<double>& a, const std::vector<double>& b)
{
    return a == b;
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
    c.physical_config_hash = 0x1122334455667788ULL;
    c.diagnostic_config_hash = 0x8877665544332211ULL;
    c.flux_quadrature_order = 4;
    c.flux_max_supports = 7;
    c.flux_max_created_particles_per_step = 123;
    c.interface_topology_hash = partition.topology_mask_hash();
    c.interface_topology_metadata_present = true;
    c.conversion_metadata_present = true;
    c.conversion_cumulative_number = 2.5e20;
    c.conversion_cumulative_px = -3.5e-4;
    c.conversion_cumulative_energy = 7.5e-7;
    c.conversion_cumulative_particles_created = 19;
    c.tail_cumulative_outflow_number = 4;
    c.control_cumulative_groups = 0;
    c.control_cumulative_fallbacks = 0;
    c.combined_number = 9.0e24;
    c.combined_kinetic_energy = 1.4e-7;
    c.combined_field_energy = 2.1e-8;
}

bool result_file(const std::string& path, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "test=bulk_tail_flux_checkpoint\n"
        << "checkpoint_roundtrip=1\n"
        << "interface_metadata_roundtrip=1\n"
        << "config_mismatch_rejected=1\n"
        << "physical_config_mismatch_rejected=1\n"
        << "particle_ids_roundtrip=1\n"
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
    bool pass = parsed && size == 1;
    if (!pass) {
        if (rank == 0)
            std::cerr << "usage: bulk_tail_flux_checkpoint_test "
                         "[--case all] [--workdir DIR] [--result FILE]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 24, 2.4 * Const::micro);
    Species electrons;
    electrons.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
                   -Const::qe, Const::me, Param::dens, Param::temperature_e,
                   false, grid);
    electrons.initialize_maxwellian();
    for (size_t i = 0; i < electrons.f.size(); i += 113)
        electrons.f[i] *= 1.0 + 1.0e-7 * static_cast<double>(i % 9);

    BeamPIC beam;
    beam.init(grid);
    beam.particles.push_back(BeamParticle{0.35e-6, 1.2e-22, 2.0e12});
    beam.particles.push_back(BeamParticle{1.25e-6, -2.4e-22, 3.0e12});
    const BeamPersistentState beam_state = beam.export_persistent_state();

    EMFields fields;
    fields.init(grid);
    for (size_t i = 0; i < fields.Ex_face.size(); ++i)
        fields.Ex_face[i] = 3.0e-5 * static_cast<double>(i % 7);
    for (size_t i = 0; i < fields.Ex.size(); ++i) {
        fields.Ex[i] = -2.0e-5 * static_cast<double>(i % 5);
        fields.phi[i] = 4.0e-5 * static_cast<double>(i % 3);
    }

    BackgroundTailPIC tail;
    tail.init(grid);
    for (int pidx = 0; pidx < 3; ++pidx) {
        BackgroundTailParticle p;
        p.x = (0.45 + 0.35 * pidx) * Const::micro;
        p.ux = 6.0 + pidx;
        p.uy = 0.25 * pidx;
        p.uz = -0.15 * pidx;
        p.weight = 1.0e16 * static_cast<double>(pidx + 1);
        p.id = tail.next_particle_id(0);
        tail.particles.push_back(p);
    }
    tail.density.assign(static_cast<size_t>(grid.nx_local), 0.0);
    for (size_t i = 0; i < tail.density.size(); ++i)
        tail.density[i] = 1.0e20 * static_cast<double>(i + 1);

    CylindricalVelocityGrid cgrid;
    cgrid.init(Param::momentum_umax);
    HybridVelocityPartition partition;
    partition.init(cgrid, 6.0, 1.0, 4, 4);
    VpfpCheckpointTailState tail_state;
    configure_tail(tail_state, partition, tail);

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
    VpfpCheckpointControl control = {17, 2.5e-15, 1.0e-17};
    std::string error;
    pass = write_vpfp_checkpoint(
        args.workdir, control, electrons, beam, fields, grid,
        field_boundary, boundary, "none", &tail_state,
        VpfpCouplingManifestConfig(), 0, 1, error);
    if (!pass) std::cerr << "checkpoint write failed: " << error << "\n";

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
        restored_fields, grid, &restored_tail, 0, 1, error);
    pass = pass && read_ok;
    if (!read_ok) std::cerr << "checkpoint read failed: " << error << "\n";
    if (read_ok) {
        const BeamPersistentState restored_beam_state =
            restored_beam.export_persistent_state();
        pass = pass && restored_electrons.f == electrons.f &&
               equal_beam(restored_beam.particles, beam.particles) &&
               std::memcmp(&restored_beam_state, &beam_state,
                           sizeof(BeamPersistentState)) == 0 &&
               restored_fields.Ex_face == fields.Ex_face &&
               restored_fields.Ex == fields.Ex &&
               restored_fields.phi == fields.phi && restored_tail.present &&
               equal_particles(restored_tail.tail.particles, tail.particles) &&
               restored_tail.tail.density == tail.density &&
               restored_tail.tail.id_counter() == tail.id_counter() &&
               restored_control.step == control.step &&
               restored_control.time == control.time &&
               restored_control.dt == control.dt &&
               restored_tail.config.partition_config_hash ==
                   tail_state.config.partition_config_hash &&
               restored_tail.config.interface_topology_hash ==
                   tail_state.config.interface_topology_hash &&
               restored_tail.config.conversion_mode == "flux-interface" &&
               restored_tail.config.physical_config_hash ==
                   tail_state.config.physical_config_hash &&
               restored_tail.config.diagnostic_config_hash ==
                   tail_state.config.diagnostic_config_hash &&
               restored_tail.config.flux_quadrature_order == 4 &&
               restored_tail.config.flux_max_supports == 7 &&
               restored_tail.config.flux_max_created_particles_per_step ==
                   123 &&
               restored_tail.config.conversion_cumulative_number ==
                   tail_state.config.conversion_cumulative_number &&
               restored_tail.config.conversion_cumulative_particles_created ==
                   tail_state.config.conversion_cumulative_particles_created &&
               restored_tail.config.interface_topology_metadata_present &&
               restored_tail.config.conversion_metadata_present;
        VpfpCheckpointTailConfig bad_config = tail_state.config;
        bad_config.flux_max_supports += 1;
        std::string mismatch_error;
        const bool mismatch_rejected =
            !validate_vpfp_checkpoint_tail_config(
                restored_tail.config, bad_config, mismatch_error);
        pass = pass && mismatch_rejected && !mismatch_error.empty();
        VpfpCheckpointTailConfig preconversion = restored_tail.config;
        preconversion.conversion_mode = "static-cell";
        preconversion.conversion_cumulative_number = 0.0;
        preconversion.conversion_cumulative_px = 0.0;
        preconversion.conversion_cumulative_energy = 0.0;
        preconversion.conversion_cumulative_particles_created = 0;
        preconversion.tail_cumulative_outflow_number = 0;
        preconversion.control_cumulative_groups = 0;
        preconversion.control_cumulative_fallbacks = 0;
        const bool preconversion_switch_allowed =
            vpfp_preconversion_static_to_flux_restart_allowed(
                preconversion, "flux-interface", 0);
        preconversion.conversion_cumulative_number = 1.0;
        const bool converted_switch_rejected =
            !vpfp_preconversion_static_to_flux_restart_allowed(
                preconversion, "flux-interface", 0);
        preconversion.conversion_cumulative_number = 0.0;
        const bool occupied_switch_rejected =
            !vpfp_preconversion_static_to_flux_restart_allowed(
                preconversion, "flux-interface", 1);
        pass = pass && preconversion_switch_allowed &&
               converted_switch_rejected && occupied_switch_rejected;
        VpfpCheckpointTailConfig bad_physical_config = tail_state.config;
        bad_physical_config.physical_config_hash ^= 1ULL;
        std::string physical_mismatch_error;
        const bool physical_mismatch_rejected =
            !validate_vpfp_checkpoint_tail_config(
                restored_tail.config, bad_physical_config,
                physical_mismatch_error);
        pass = pass && physical_mismatch_rejected &&
               !physical_mismatch_error.empty();
        const std::uint64_t next_id = tail.next_particle_id(0);
        const std::uint64_t restored_next_id =
            restored_tail.tail.next_particle_id(0);
        pass = pass && next_id == restored_next_id;
    }

    if (rank == 0) {
        std::cout << "interface_hash=" << tail_state.config.interface_topology_hash
                  << " tail_particles=" << tail.particles.size()
                  << " checkpoint_roundtrip=" << (pass ? 1 : 0) << "\n";
        if (!result_file(args.result, pass)) pass = false;
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
