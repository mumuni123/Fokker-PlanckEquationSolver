// Stage H5 acceptance for TailPopulationController (sections 7.10 and 15
// H5): deterministic conservative compression and equal-weight splitting on
// the tail trial representation, with the shared 7-moment constraint module.
//
// Cases:
//   below-target-no-op                 : count <= target -> no operation
//   above-max-conservative-compression : count > max -> <=7 supports,
//                                        N/Px/Jx/K/Pixx/Piperp/Xw preserved,
//                                        CIC density unchanged, deterministic
//   weight-ratio-split                 : heavy particle split into
//                                        near-equal children (moments exact)
//   boundary-cell-skip                 : first/last cells never merged
//   compression-failure-keeps-original : fallback keeps the original
//                                        particles bitwise
//
// Usage:
//   tail_population_controller_test [--case <case|all>] [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".

#include "background_tail_pic.h"
#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include "tail_moment_constraint.h"
#include "tail_population_controller.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TestArgs {
    std::string test_case;
    std::string result_path;
};

bool parse_args(int argc, char** argv, TestArgs& args)
{
    args.test_case = "all";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case") {
            if (i + 1 >= argc) return false;
            args.test_case = argv[++i];
        } else if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return !args.test_case.empty();
}

struct Metrics {
    bool pass;
    std::uint64_t particles_before;
    std::uint64_t particles_after;
    int groups_considered;
    int groups_below_target;
    int groups_compressed;
    int groups_split;
    int boundary_groups_skipped;
    int compression_fallback_count;
    int split_capped_count;
    double max_residual[7];
    double rho_l2_rel;
    bool bitwise_unchanged;
    Metrics()
        : pass(false), particles_before(0), particles_after(0),
          groups_considered(0), groups_below_target(0),
          groups_compressed(0), groups_split(0), boundary_groups_skipped(0),
          compression_fallback_count(0), split_capped_count(0),
          rho_l2_rel(0.0), bitwise_unchanged(true)
    {
        for (int i = 0; i < 7; ++i) max_residual[i] = 0.0;
    }
};

// Deterministic pseudo-random tail particle in `cell` (interior), all in the
// same phase-space bin for the 6 MeV / 4x4 partition.
void add_bin_particle(BackgroundTailPIC& tail, const SpatialGrid& grid,
                      int cell, int variant, double weight,
                      std::uint64_t id)
{
    BackgroundTailParticle p;
    const double xc = (static_cast<double>(cell) + 0.5) * grid.dx;
    // Offsets in [0.15, 0.35] dx keep every particle inside `cell` so the
    // whole bin maps to one spatial group (same CIC stencil).
    p.x = xc + (0.15 + 0.05 * static_cast<double>(variant % 5)) * grid.dx;
    p.ux = 7.0 + 0.01 * static_cast<double>(variant % 3);
    p.uy = 0.3 + 0.1 * static_cast<double>(variant % 4);
    p.uz = 0.2 + 0.05 * static_cast<double>(variant % 5);
    p.weight = weight;
    p.id = id;
    tail.particles.push_back(p);
}

void tail_moments_sum(const BackgroundTailPIC& tail, TailMoment7& m)
{
    m = TailMoment7();
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        TailMoment7 pm;
        tail_particle_moments(p.weight, p.x, p.ux, p.uy, p.uz, pm);
        m.n += pm.n;
        m.px += pm.px;
        m.jx += pm.jx;
        m.ke += pm.ke;
        m.pixx += pm.pixx;
        m.piperp += pm.piperp;
        m.xw += pm.xw;
    }
}

double rho_l2_rel(const BackgroundTailPIC& a, const BackgroundTailPIC& b)
{
    double l2a = 0.0;
    double l2d = 0.0;
    const size_t n = std::min(a.density.size(), b.density.size());
    for (size_t i = 0; i < n; ++i) {
        const double da = a.density[i];
        const double db = b.density[i];
        l2a += da * da;
        l2d += (da - db) * (da - db);
    }
    return std::sqrt(l2d) / std::max(1.0, std::sqrt(l2a));
}

bool particles_bitwise_equal(const BackgroundTailPIC& a,
                             const BackgroundTailPIC& b)
{
    if (a.particles.size() != b.particles.size()) return false;
    for (size_t i = 0; i < a.particles.size(); ++i) {
        const BackgroundTailParticle& pa = a.particles[i];
        const BackgroundTailParticle& pb = b.particles[i];
        if (pa.x != pb.x || pa.ux != pb.ux || pa.uy != pb.uy ||
            pa.uz != pb.uz || pa.weight != pb.weight || pa.id != pb.id) {
            return false;
        }
    }
    return true;
}

// Multiset comparison by ID: the controller may reorder particles when it
// groups them by (cell, phase bin); "keep the original particles" means the
// same particle set, not the same array order.
bool particles_multiset_equal(const BackgroundTailPIC& a,
                              const BackgroundTailPIC& b)
{
    if (a.particles.size() != b.particles.size()) return false;
    std::vector<BackgroundTailParticle> pa = a.particles;
    std::vector<BackgroundTailParticle> pb = b.particles;
    std::sort(pa.begin(), pa.end(),
              [](const BackgroundTailParticle& x,
                 const BackgroundTailParticle& y) { return x.id < y.id; });
    std::sort(pb.begin(), pb.end(),
              [](const BackgroundTailParticle& x,
                 const BackgroundTailParticle& y) { return x.id < y.id; });
    for (size_t i = 0; i < pa.size(); ++i) {
        if (pa[i].x != pb[i].x || pa[i].ux != pb[i].ux ||
            pa[i].uy != pb[i].uy || pa[i].uz != pb[i].uz ||
            pa[i].weight != pb[i].weight || pa[i].id != pb[i].id) {
            return false;
        }
    }
    return true;
}

Metrics run_below_target()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);
    for (int i = 0; i < 10; ++i) {
        add_bin_particle(tail, grid, 50, i,
                         1.0e14 + 1.0e12 * static_cast<double>(i),
                         static_cast<std::uint64_t>(100 + i));
    }
    const BackgroundTailPIC before = tail;
    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 16;
    config.max_particles_per_phase_bin = 32;
    config.max_weight_ratio = 8.0;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, 0);
    m.pass = d.applied && d.groups_considered == 1 &&
             d.groups_below_target == 1 && d.groups_compressed == 0 &&
             d.groups_split == 0 && d.compression_fallback_count == 0 &&
             particles_bitwise_equal(tail, before);
    m.particles_before = d.particles_before_local;
    m.particles_after = d.particles_after_local;
    m.groups_considered = d.groups_considered;
    m.groups_below_target = d.groups_below_target;
    return m;
}

Metrics run_above_max()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);
    for (int i = 0; i < 40; ++i) {
        add_bin_particle(tail, grid, 50, i,
                         1.0e14 + 5.0e12 * static_cast<double>(i % 7),
                         static_cast<std::uint64_t>(200 + i));
    }
    TailMoment7 ref;
    tail_moments_sum(tail, ref);
    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 8;
    config.max_particles_per_phase_bin = 16;
    config.max_weight_ratio = 8.0;
    config.max_support = 7;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, 0);
    TailMoment7 got;
    tail_moments_sum(tail, got);
    const double scales[7] = {
        std::max(1.0, std::fabs(ref.n)),
        std::max(1.0, std::fabs(ref.px)),
        std::max(1.0, std::fabs(ref.jx)),
        std::max(1.0, std::fabs(ref.ke)),
        std::max(1.0, std::fabs(ref.pixx)),
        std::max(1.0, std::fabs(ref.piperp)),
        std::max(1.0, std::fabs(ref.xw))
    };
    const double residuals[7] = {
        std::fabs(got.n - ref.n) / scales[0],
        std::fabs(got.px - ref.px) / scales[1],
        std::fabs(got.jx - ref.jx) / scales[2],
        std::fabs(got.ke - ref.ke) / scales[3],
        std::fabs(got.pixx - ref.pixx) / scales[4],
        std::fabs(got.piperp - ref.piperp) / scales[5],
        std::fabs(got.xw - ref.xw) / scales[6]
    };
    for (int i = 0; i < 7; ++i) m.max_residual[i] = residuals[i];

    // CIC density invariance (section 7.10: the same-cell N/Xw preservation
    // makes the two-cell CIC deposit unchanged to the compression error).
    BackgroundTailPIC tail_before;
    tail_before.init(grid);
    for (int i = 0; i < 40; ++i) {
        add_bin_particle(tail_before, grid, 50, i,
                         1.0e14 + 5.0e12 * static_cast<double>(i % 7),
                         static_cast<std::uint64_t>(200 + i));
    }
    tail_before.deposit_density(grid, 0, 1);
    tail.deposit_density(grid, 0, 1);
    m.rho_l2_rel = rho_l2_rel(tail_before, tail);

    // Determinism: re-apply on an identical input and compare bitwise.
    BackgroundTailPIC tail2;
    tail2.init(grid);
    for (int i = 0; i < 40; ++i) {
        add_bin_particle(tail2, grid, 50, i,
                         1.0e14 + 5.0e12 * static_cast<double>(i % 7),
                         static_cast<std::uint64_t>(200 + i));
    }
    controller.apply(tail2, grid, partition, 1, 0);
    const bool deterministic = particles_bitwise_equal(tail, tail2);
    const bool residual_ok =
        d.max_residual[0] <= 1.0e-9 &&
        d.max_residual[1] <= 1.0e-9 &&
        d.max_residual[2] <= 1.0e-9 &&
        d.max_residual[3] <= 1.0e-9 &&
        d.max_residual[4] <= 1.0e-9 &&
        d.max_residual[5] <= 1.0e-9 &&
        d.max_residual[6] <= 1.0e-9 &&
        residuals[0] <= 1.0e-9 && residuals[1] <= 1.0e-9 &&
        residuals[2] <= 1.0e-9 && residuals[3] <= 1.0e-9 &&
        residuals[4] <= 1.0e-9 && residuals[5] <= 1.0e-9 &&
        residuals[6] <= 1.0e-9;
    m.pass = d.applied && d.groups_compressed == 1 &&
             d.compression_fallback_count == 0 &&
             d.particles_after_local <= 16 && residual_ok &&
             m.rho_l2_rel <= 1.0e-8 && deterministic &&
             d.particles_after_local < d.particles_before_local;
    m.particles_before = d.particles_before_local;
    m.particles_after = d.particles_after_local;
    m.groups_considered = d.groups_considered;
    m.groups_compressed = d.groups_compressed;
    return m;
}

Metrics run_weight_ratio_split()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);
    // Five light (1e14) and one heavy (1e16): ratio 100 > 8, count 6 between
    // target (4) and max (16), so only the weight-ratio split applies.
    for (int i = 0; i < 5; ++i) {
        add_bin_particle(tail, grid, 50, i, 1.0e14,
                         static_cast<std::uint64_t>(300 + i));
    }
    add_bin_particle(tail, grid, 50, 7, 1.0e16, 310);
    TailMoment7 ref;
    tail_moments_sum(tail, ref);
    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 4;
    config.max_particles_per_phase_bin = 16;
    config.max_weight_ratio = 8.0;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, 0);
    TailMoment7 got;
    tail_moments_sum(tail, got);
    const double scales[7] = {
        std::max(1.0, std::fabs(ref.n)),
        std::max(1.0, std::fabs(ref.px)),
        std::max(1.0, std::fabs(ref.jx)),
        std::max(1.0, std::fabs(ref.ke)),
        std::max(1.0, std::fabs(ref.pixx)),
        std::max(1.0, std::fabs(ref.piperp)),
        std::max(1.0, std::fabs(ref.xw))
    };
    const double residuals[7] = {
        std::fabs(got.n - ref.n) / scales[0],
        std::fabs(got.px - ref.px) / scales[1],
        std::fabs(got.jx - ref.jx) / scales[2],
        std::fabs(got.ke - ref.ke) / scales[3],
        std::fabs(got.pixx - ref.pixx) / scales[4],
        std::fabs(got.piperp - ref.piperp) / scales[5],
        std::fabs(got.xw - ref.xw) / scales[6]
    };
    for (int i = 0; i < 7; ++i) m.max_residual[i] = residuals[i];

    // Weight ratio bound: max/min <= 8 (children near-equal).
    double w_min = 0.0;
    double w_max = 0.0;
    std::set<std::uint64_t> ids;
    bool ids_unique = true;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const double w = tail.particles[i].weight;
        if (i == 0 || w < w_min) w_min = w;
        if (i == 0 || w > w_max) w_max = w;
        if (!ids.insert(tail.particles[i].id).second) ids_unique = false;
    }
    const double ratio = w_min > 0.0 ? w_max / w_min : 0.0;
    const bool moments_exact =
        residuals[0] <= 1.0e-12 && residuals[1] <= 1.0e-12 &&
        residuals[2] <= 1.0e-12 && residuals[3] <= 1.0e-12 &&
        residuals[4] <= 1.0e-12 && residuals[5] <= 1.0e-12 &&
        residuals[6] <= 1.0e-12;
    m.pass = d.applied && d.groups_split == 1 &&
             d.groups_compressed == 0 && d.compression_fallback_count == 0 &&
             d.particles_after_local > d.particles_before_local &&
             ratio <= 8.0 + 1.0e-12 && moments_exact && ids_unique;
    m.particles_before = d.particles_before_local;
    m.particles_after = d.particles_after_local;
    m.groups_considered = d.groups_considered;
    m.groups_split = d.groups_split;
    return m;
}

Metrics run_boundary_cell_skip()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);
    for (int i = 0; i < 12; ++i) {
        add_bin_particle(tail, grid, 0, i, 1.0e14,
                         static_cast<std::uint64_t>(400 + i));
        add_bin_particle(tail, grid, 199, i, 1.0e14,
                         static_cast<std::uint64_t>(500 + i));
    }
    const BackgroundTailPIC before = tail;
    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 4;
    config.max_particles_per_phase_bin = 8;
    config.max_weight_ratio = 8.0;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, 0);
    const bool unchanged = particles_multiset_equal(tail, before);
    m.pass = d.applied && d.boundary_groups_skipped == 2 &&
             d.groups_compressed == 0 && d.groups_split == 0 &&
             unchanged;
    m.particles_before = d.particles_before_local;
    m.particles_after = d.particles_after_local;
    m.groups_considered = d.groups_considered;
    m.boundary_groups_skipped = d.boundary_groups_skipped;
    return m;
}

Metrics run_compression_failure()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 200, 2.0 * Const::micro);
    Species bulk;
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    BackgroundTailPIC tail;
    tail.init(grid);
    // Two non-proportional supports (different x): with max_support=1 the
    // shared 7-row compression has no null direction (2 columns, rank 2),
    // so it fails and the original particles must be kept untouched.
    add_bin_particle(tail, grid, 50, 0, 1.0e14, 600);
    add_bin_particle(tail, grid, 50, 2, 1.0e14, 601);
    const BackgroundTailPIC before = tail;
    TailPopulationController controller;
    TailPopulationController::Config config;
    config.enabled = true;
    config.control_interval = 1;
    config.target_particles_per_phase_bin = 1;
    config.max_particles_per_phase_bin = 1;
    config.max_weight_ratio = 8.0;
    config.max_support = 1;
    controller.configure(config);
    const TailPopulationController::Diagnostics d =
        controller.apply(tail, grid, partition, 1, 0);
    m.pass = d.applied && d.compression_fallback_count == 1 &&
             d.groups_compressed == 0 &&
             particles_multiset_equal(tail, before);
    m.particles_before = d.particles_before_local;
    m.particles_after = d.particles_after_local;
    m.groups_considered = d.groups_considered;
    m.compression_fallback_count = d.compression_fallback_count;
    return m;
}

void write_result(const std::string& path, const std::string& test_case,
                  const Metrics& m)
{
    if (path.empty()) return;
    std::ofstream out(path.c_str(), std::ios::app);
    if (!out) return;
    out << "case=" << test_case
        << " pass=" << (m.pass ? 1 : 0)
        << " particles_before=" << m.particles_before
        << " particles_after=" << m.particles_after
        << " groups_considered=" << m.groups_considered
        << " groups_below_target=" << m.groups_below_target
        << " groups_compressed=" << m.groups_compressed
        << " groups_split=" << m.groups_split
        << " boundary_groups_skipped=" << m.boundary_groups_skipped
        << " compression_fallback_count=" << m.compression_fallback_count
        << " split_capped_count=" << m.split_capped_count
        << " rho_l2_rel=" << std::setprecision(17) << m.rho_l2_rel
        << " residual_n=" << m.max_residual[0]
        << " residual_px=" << m.max_residual[1]
        << " residual_jx=" << m.max_residual[2]
        << " residual_ke=" << m.max_residual[3]
        << " residual_pixx=" << m.max_residual[4]
        << " residual_piperp=" << m.max_residual[5]
        << " residual_xw=" << m.max_residual[6] << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    TestArgs args;
    const bool parsed = parse_args(argc, argv, args);
    if (!parsed) {
        std::cerr << "usage: tail_population_controller_test "
                     "[--case <case|all>] [--result <path>]\n";
        MPI_Finalize();
        return 2;
    }
    bool pass = true;
    Metrics all_metrics[5];
    const char* names[5] = {
        "below-target-no-op",
        "above-max-conservative-compression",
        "weight-ratio-split",
        "boundary-cell-skip",
        "compression-failure-keeps-original"
    };
    Metrics (*runs[5])() = {
        &run_below_target,
        &run_above_max,
        &run_weight_ratio_split,
        &run_boundary_cell_skip,
        &run_compression_failure
    };
    for (int i = 0; i < 5; ++i) {
        const bool selected =
            args.test_case == "all" || args.test_case == names[i];
        if (!selected) continue;
        all_metrics[i] = runs[i]();
        if (!all_metrics[i].pass) pass = false;
        write_result(args.result_path, names[i], all_metrics[i]);
        std::cout << names[i] << "=" << (all_metrics[i].pass ? 1 : 0)
                  << " particles=" << all_metrics[i].particles_before
                  << "->" << all_metrics[i].particles_after
                  << " groups=" << all_metrics[i].groups_considered
                  << " below=" << all_metrics[i].groups_below_target
                  << " compressed=" << all_metrics[i].groups_compressed
                  << " split=" << all_metrics[i].groups_split
                  << " boundary_skipped="
                  << all_metrics[i].boundary_groups_skipped
                  << " fallbacks="
                  << all_metrics[i].compression_fallback_count
                  << " rho_l2_rel=" << std::setprecision(17)
                  << all_metrics[i].rho_l2_rel
                  << " res_n=" << all_metrics[i].max_residual[0]
                  << " res_xw=" << all_metrics[i].max_residual[6] << "\n";
    }
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return pass ? 0 : 1;
}
