#include "tail_bulk_return_test_common.h"

#include "beam_pic.h"
#include "open_boundary.h"
#include "open_electrostatic_solver.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

namespace {
struct Args { std::string workdir; std::string result; Args() : workdir("tail_return_checkpoint_tmp") {} };
bool parse(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        const std::string s(argv[i]);
        if (s == "--workdir" && i + 1 < argc) a.workdir = argv[++i];
        else if (s == "--result" && i + 1 < argc) a.result = argv[++i];
        else if (s == "--case" && i + 1 < argc &&
                 std::string(argv[i + 1]) == "all") ++i;
        else return false;
    }
    return true;
}
void remove_dir(const std::string& p) {
#ifdef _WIN32
    _rmdir(p.c_str());
#else
    rmdir(p.c_str());
#endif
}
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    bool pass = size == 1 && parse(argc, argv, args);
    SpatialGrid grid; Species bulk; HybridVelocityPartition partition;
    BackgroundTailPIC tail;
    tail_return_test::init_state(0, 1, grid, bulk, partition, tail);
    const std::pair<int, int> slot =
        tail_return_test::safe_velocity_slot(bulk, partition);
    if (slot.first < 0) pass = false;
    if (pass) tail.particles.push_back(tail_return_test::make_particle(
        bulk, slot.first, slot.second, 10.5 * grid.dx, 2.0e12, 17, 6));
    BeamPIC beam; beam.init(grid);
    EMFields fields; fields.init(grid);
    ElectrostaticBoundary field_bc = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    OpenBackgroundBoundaryConfig bkg_bc;
    bkg_bc.left_type = BackgroundXBoundaryType::RESERVOIR;
    bkg_bc.right_type = BackgroundXBoundaryType::RESERVOIR;
    bkg_bc.left_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    bkg_bc.right_reservoir = { Param::dens, Param::temperature_e, 0.0 };
    VpfpCheckpointTailState state;
    state.present = true; state.tail = tail;
    state.config.partition_config_hash = partition.config_hash;
    state.config.convert_energy_mev = 6.0;
    state.config.buffer_width_mev = 1.0;
    state.config.upar_bins = 4; state.config.energy_bins = 4;
    state.config.return_mode = "project";
    state.config.return_energy_mev = 5.5;
    state.config.return_residence_steps = 8;
    state.config.return_max_stencil_radius = 3;
    state.config.return_moment_tolerance = 1.0e-12;
    state.config.return_cumulative_number = 1.25e18;
    state.config.return_cumulative_px = -2.5e-4;
    state.config.return_cumulative_jx_dx = 3.75e7;
    state.config.return_cumulative_energy = 4.5e5;
    state.config.return_cumulative_pixx_dx = 5.25e4;
    state.config.return_cumulative_piperp_dx = 6.75e4;
    state.config.return_cumulative_particles_removed = 1234;
    state.config.return_cumulative_deferred_groups = 17;
    VpfpCheckpointControl control = { 77, 2.0e-15, 2.5e-17 };
    std::string error;
    pass = pass && write_vpfp_checkpoint(
        args.workdir, control, bulk, beam, fields, grid, field_bc, bkg_bc,
        "none", &state, VpfpCouplingManifestConfig(), 0, 1, error);
    Species got_bulk; got_bulk.init("background_electrons",
        SpeciesType::BACKGROUND_ELECTRON, -Const::qe, Const::me, Param::dens,
        Param::temperature_e, false, grid);
    BeamPIC got_beam; got_beam.init(grid);
    EMFields got_fields; got_fields.init(grid);
    VpfpCheckpointTailState got;
    VpfpCheckpointControl got_control = {};
    pass = pass && read_vpfp_checkpoint(args.workdir, got_control, got_bulk,
        got_beam, got_fields, grid, &got, 0, 1, error);
    const bool residence_equal = pass && got.present &&
        got.tail.particles.size() == 1 &&
        got.tail.particles[0].return_residence_steps == 6;
    const bool config_equal = pass && got.config.return_mode == "project" &&
        got.config.return_energy_mev == 5.5 &&
        got.config.return_residence_steps == 8 &&
        got.config.return_max_stencil_radius == 3 &&
        got.config.return_moment_tolerance == 1.0e-12 &&
        got.config.return_cumulative_number == 1.25e18 &&
        got.config.return_cumulative_px == -2.5e-4 &&
        got.config.return_cumulative_jx_dx == 3.75e7 &&
        got.config.return_cumulative_energy == 4.5e5 &&
        got.config.return_cumulative_pixx_dx == 5.25e4 &&
        got.config.return_cumulative_piperp_dx == 6.75e4 &&
        got.config.return_cumulative_particles_removed == 1234 &&
        got.config.return_cumulative_deferred_groups == 17;
    struct LegacyParticleRecord {
        double x, ux, uy, uz, weight;
        std::uint64_t id;
    } legacy = { 1.25e-7, 2.0, 0.5, -0.25, 9.0e11, 321 };
    std::stringstream legacy_stream(
        std::ios::in | std::ios::out | std::ios::binary);
    legacy_stream.write(reinterpret_cast<const char*>(&legacy), sizeof(legacy));
    legacy_stream.seekg(0);
    std::vector<BackgroundTailParticle> legacy_particles;
    const bool legacy_reset = read_vpfp_tail_particle_records(
        legacy_stream, 2, 1, legacy_particles) &&
        legacy_particles.size() == 1 &&
        legacy_particles[0].x == legacy.x &&
        legacy_particles[0].id == legacy.id &&
        legacy_particles[0].return_residence_steps == 0;
    pass = pass && residence_equal && config_equal &&
        legacy_reset &&
        got_control.step == control.step && got_control.time == control.time;
    std::remove((args.workdir + "/rank_000000.bin").c_str());
    std::remove((args.workdir + "/manifest.txt").c_str());
    remove_dir(args.workdir);
    if (rank == 0) {
        tail_return_test::write_result(args.result, {
            {"residence_roundtrip_equal", residence_equal ? 1.0 : 0.0},
            {"return_config_roundtrip_equal", config_equal ? 1.0 : 0.0},
            {"return_ledger_roundtrip_equal", config_equal ? 1.0 : 0.0},
            {"legacy_v1_residence_reset_covered", legacy_reset ? 1.0 : 0.0}}, pass);
        std::cout << "status=" << (pass ? "PASS" : "FAIL") << '\n';
        if (!pass && !error.empty()) std::cerr << error << '\n';
    }
    MPI_Finalize();
    return pass ? 0 : 1;
}
