// Stage H8 acceptance (section 10.3.6 item 1): zero-dimensional
// isotropisation of an equal-weight tail with T_parallel != T_perp.  The
// temperature anisotropy must decay, the result must be stable under
// dt -> dt/2 (same total collision time), and the weighted 3-momentum and
// relativistic energy must be conserved.
//
// Usage:
//   tail_collision_isotropisation_test [--result <path>]
// The last stdout line is always "status=PASS" or "status=FAIL".
//
// Statistical methodology (section 10.3.6): every configuration is run over
// N_SEEDS counter-based seeds and the anisotropy is judged on the SEED MEAN,
// never on a single noise curve.  The pass thresholds are the same physical
// values as the single-draw test (mean_a < 0.2*aniso_initial,
// |mean_a - mean_b| < 0.02, mean_c <= mean_a + 1e-3, per-run P/K
// conservation); only the estimator is changed from one draw to a mean.
//
// Diagnostics: the run reports the per-seed scan, the anisotropy decay curve
// of the first seed, the first-call substep count / max_s12 / large-angle
// fraction, the total number of collision substeps and the per-config mean
// and standard deviation, so a failure can be attributed to (a) sampling
// noise, (b) a dt-scaling defect of the kernel, or (c) a conservation
// defect.

#include "background_tail_collision.h"
#include "background_tail_nanbu_perez.h"
#include "grid.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

const int N_SEEDS = 16;

inline std::uint64_t mix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline double uniform01(std::uint64_t& state)
{
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    return static_cast<double>(z >> 11) *
           (1.0 / 9007199254740992.0);
}

struct Args {
    std::string result_path;
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result") {
            if (i + 1 >= argc) return false;
            args.result_path = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

// Normalised temperature anisotropy of the current particle list:
//   |T_par - T_perp| / (T_par + 2 T_perp),
// with T_par = <w ux^2>/sum(w) and
//      T_perp = 0.5 (<w uy^2> + <w uz^2>)/sum(w).
double anisotropy(const BackgroundTailPIC& tail)
{
    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        sxx += p.weight * p.ux * p.ux;
        syy += p.weight * p.uy * p.uy;
        szz += p.weight * p.uz * p.uz;
    }
    const double s = sxx + syy + szz;
    return std::fabs(sxx - 0.5 * (syy + szz)) / std::max(1.0e-30, s);
}

struct IsoRun {
    double aniso_final;
    double p_err;
    double k_err;
    int first_substeps;
    double first_max_s12;
    double first_large_angle_fraction;
    std::uint64_t substeps_total;
    double max_s12_overall;
    double aniso_quarter;
    double aniso_half;
    std::vector<double> aniso_traj;  // index 0 = initial state
    IsoRun()
        : aniso_final(-1.0), p_err(-1.0), k_err(-1.0), first_substeps(-1),
          first_max_s12(-1.0), first_large_angle_fraction(-1.0),
          substeps_total(0), max_s12_overall(0.0), aniso_quarter(-1.0),
          aniso_half(-1.0)
    {}
};

struct IsoConfig {
    int n_particles;
    double dt;
    int steps;
    std::uint64_t seed_base;
    IsoRun first;                  // seed 0 details (trajectory etc.)
    std::vector<double> aniso_per_seed;
    double aniso_mean;
    double aniso_std;
    double p_err_max;
    double k_err_max;
};

IsoRun run_one_seed(const SpatialGrid& grid, int n_particles, double dt,
                    int steps, std::uint64_t seed_base)
{
    IsoRun r;
    BackgroundTailPIC tail;
    tail.init(grid);
    std::uint64_t seed = seed_base;
    const double sigma_par = 0.6;
    const double sigma_perp = 0.2;
    for (int i = 0; i < n_particles; ++i) {
        std::uint64_t s = mix64(seed + static_cast<std::uint64_t>(i));
        const double u0 = 0.5 * std::sqrt(
            -2.0 * std::log(1.0e-15 + uniform01(s)));
        const double phi = 2.0 * Const::pi * uniform01(s);
        const double cz = 2.0 * uniform01(s) - 1.0;
        const double sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        BackgroundTailParticle p;
        p.x = 0.4e-6;
        // Anisotropic: parallel spread larger than perpendicular.
        p.ux = sigma_par * u0 * cz;
        p.uy = sigma_perp * u0 * sz * std::cos(phi);
        p.uz = sigma_perp * u0 * sz * std::sin(phi);
        p.weight = 1.0e20;
        p.id = static_cast<std::uint64_t>(i + 1);
        tail.particles.push_back(p);
    }
    double px0, py0, pz0, ke0;
    px0 = py0 = pz0 = ke0 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        px0 += Const::me * Const::c * p.weight * p.ux;
        py0 += Const::me * Const::c * p.weight * p.uy;
        pz0 += Const::me * Const::c * p.weight * p.uz;
        ke0 += Const::me * Const::c * Const::c * p.weight * (gamma - 1.0);
    }
    r.aniso_traj.push_back(anisotropy(tail));

    TailCollisionRequest request;
    request.kernel = TailCollisionKernel::CoulombLandauNanbuPerez;
    request.dt = dt;
    request.coulomb_log = 20.0;
    request.rng_seed_base = seed_base;
    request.mpi_rank = 0;
    request.max_particle_growth = 0.0;
    for (int step = 0; step < steps; ++step) {
        request.accepted_step = step + 1;
        request.collision_half = step % 2;
        TailCollisionDiagnostics diag;
        if (!nanbu_perez_collide(tail, grid, request, diag)) {
            r.aniso_final = 1.0e30;
            r.p_err = 1.0e30;
            r.k_err = 1.0e30;
            return r;
        }
        r.substeps_total += static_cast<std::uint64_t>(
            diag.collision_substeps);
        r.max_s12_overall = std::max(r.max_s12_overall, diag.max_s12);
        if (step == 0) {
            r.first_substeps = diag.collision_substeps;
            r.first_max_s12 = diag.max_s12;
            r.first_large_angle_fraction = diag.large_angle_fraction;
        }
        r.aniso_traj.push_back(anisotropy(tail));
        if (step == steps / 4 - 1) r.aniso_quarter = anisotropy(tail);
        if (step == steps / 2 - 1) r.aniso_half = anisotropy(tail);
    }
    double px1, py1, pz1, ke1;
    px1 = py1 = pz1 = ke1 = 0.0;
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        const double gamma = std::sqrt(
            1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        px1 += Const::me * Const::c * p.weight * p.ux;
        py1 += Const::me * Const::c * p.weight * p.uy;
        pz1 += Const::me * Const::c * p.weight * p.uz;
        ke1 += Const::me * Const::c * Const::c * p.weight * (gamma - 1.0);
    }
    r.aniso_final = anisotropy(tail);
    const double p_scale = std::max(
        1.0, std::fabs(px0) + std::fabs(py0) + std::fabs(pz0));
    r.p_err = (std::fabs(px1 - px0) + std::fabs(py1 - py0) +
               std::fabs(pz1 - pz0)) / p_scale;
    r.k_err = std::fabs(ke1 - ke0) / std::max(1.0, ke0);
    return r;
}

IsoConfig run_isotropisation(const SpatialGrid& grid, int n_particles,
                             double dt, int steps, std::uint64_t seed_base)
{
    IsoConfig cfg;
    cfg.n_particles = n_particles;
    cfg.dt = dt;
    cfg.steps = steps;
    cfg.seed_base = seed_base;
    cfg.p_err_max = 0.0;
    cfg.k_err_max = 0.0;
    double sum = 0.0;
    double sumsq = 0.0;
    for (int k = 0; k < N_SEEDS; ++k) {
        const std::uint64_t seed =
            seed_base + static_cast<std::uint64_t>(k);
        const IsoRun r = run_one_seed(grid, n_particles, dt, steps, seed);
        cfg.aniso_per_seed.push_back(r.aniso_final);
        sum += r.aniso_final;
        sumsq += r.aniso_final * r.aniso_final;
        cfg.p_err_max = std::max(cfg.p_err_max, r.p_err);
        cfg.k_err_max = std::max(cfg.k_err_max, r.k_err);
        if (k == 0) cfg.first = r;
    }
    cfg.aniso_mean = sum / static_cast<double>(N_SEEDS);
    const double var = std::max(
        0.0, sumsq / static_cast<double>(N_SEEDS) -
                 cfg.aniso_mean * cfg.aniso_mean);
    cfg.aniso_std = std::sqrt(var);
    return cfg;
}

void print_trajectory(const char* label, const IsoRun& r, int stride)
{
    std::cout << "  " << label << " aniso decay (seed0):";
    const int n = static_cast<int>(r.aniso_traj.size());
    for (int s = 0; s < n; ++s) {
        if (s <= 2 || s % stride == 0 || s == n - 1) {
            std::cout << " " << s << ":" << r.aniso_traj[static_cast<size_t>(s)];
        }
    }
    std::cout << "\n";
}

void print_config(const char* label, const IsoConfig& cfg)
{
    std::cout << "run " << label << ": n=" << cfg.n_particles
              << " dt=" << cfg.dt << " steps=" << cfg.steps
              << " seeds=0x" << std::hex << cfg.seed_base << std::dec
              << "..0x" << std::hex
              << (cfg.seed_base + static_cast<std::uint64_t>(N_SEEDS - 1))
              << std::dec << "\n";
    std::cout << "  first_call: substeps=" << cfg.first.first_substeps
              << " max_s12=" << cfg.first.first_max_s12
              << " large_angle_frac="
              << cfg.first.first_large_angle_fraction << "\n";
    std::cout << "  total_substeps(seed0)=" << cfg.first.substeps_total
              << " max_s12_overall(seed0)=" << cfg.first.max_s12_overall
              << " aniso_quarter=" << cfg.first.aniso_quarter
              << " aniso_half=" << cfg.first.aniso_half
              << " aniso_seed0=" << cfg.first.aniso_final << "\n";
    std::cout << "  aniso_mean=" << cfg.aniso_mean
              << " aniso_std=" << cfg.aniso_std
              << " p_err_max=" << cfg.p_err_max
              << " k_err_max=" << cfg.k_err_max << "\n";
    print_trajectory(label, cfg.first, label[0] == 'B' ? 10 : 5);
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 1) {
        std::cerr << "tail_collision_isotropisation_test: single-rank only\n";
        MPI_Finalize();
        return 2;
    }
    Args args;
    if (!parse_args(argc, argv, args)) {
        std::cerr << "usage: tail_collision_isotropisation_test "
                     "[--result <path>]\n";
        MPI_Finalize();
        return 2;
    }

    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    const double dt = 1.0e-13;
    const int steps = 50;
    const IsoConfig ca = run_isotropisation(grid, 400, dt, steps, 0x1111);
    const IsoConfig cb = run_isotropisation(
        grid, 400, 0.5 * dt, 2 * steps, 0x2222);
    const IsoConfig cc = run_isotropisation(grid, 1000, dt, steps, 0x3333);

    const double aniso_initial =
        std::fabs(0.6 * 0.6 - 0.5 * (0.2 * 0.2 + 0.2 * 0.2)) /
        std::max(1.0e-30, 0.6 * 0.6 + 0.2 * 0.2 + 0.2 * 0.2);
    // Same physical thresholds as the original single-draw test, now applied
    // to the seed means (section 10.3.6: never judge on one noise curve).
    const bool iso = ca.aniso_mean < 0.2 * aniso_initial &&
                     std::fabs(ca.aniso_mean - cb.aniso_mean) < 0.02 &&
                     cc.aniso_mean <= ca.aniso_mean + 1.0e-3 &&
                     ca.p_err_max < 1.0e-12 && ca.k_err_max < 1.0e-12 &&
                     cb.p_err_max < 1.0e-12 && cb.k_err_max < 1.0e-12 &&
                     cc.p_err_max < 1.0e-12 && cc.k_err_max < 1.0e-12;

    std::cout.precision(12);
    std::cout << "aniso_initial=" << aniso_initial << "\n";
    std::cout << "seed scan (aniso_final, " << N_SEEDS << " seeds):\n";
    for (int k = 0; k < N_SEEDS; ++k) {
        const double a = ca.aniso_per_seed[static_cast<size_t>(k)];
        const double b = cb.aniso_per_seed[static_cast<size_t>(k)];
        const double c = cc.aniso_per_seed[static_cast<size_t>(k)];
        std::cout << "  seed+" << k << ": a=" << a << " b=" << b
                  << " c=" << c << " |a-b|=" << std::fabs(a - b) << "\n";
    }
    print_config("A", ca);
    print_config("B", cb);
    print_config("C", cc);

    std::cout << "isotropisation=" << (iso ? 1 : 0)
              << " aniso_initial=" << aniso_initial
              << " final_mean(n400,dt)=" << ca.aniso_mean
              << " final_mean(n400,dt/2)=" << cb.aniso_mean
              << " final_mean(n1000,dt)=" << cc.aniso_mean
              << " std_a=" << ca.aniso_std
              << " std_b=" << cb.aniso_std
              << " std_c=" << cc.aniso_std
              << " p_err_max=" << ca.p_err_max
              << " k_err_max=" << ca.k_err_max << "\n";

    if (!args.result_path.empty()) {
        std::ofstream out(args.result_path.c_str(), std::ios::app);
        if (out) {
            out << "case=tail-collision-isotropisation pass="
                << (iso ? 1 : 0)
                << " aniso_final=" << ca.aniso_mean
                << " aniso_dt_half=" << cb.aniso_mean
                << " aniso_n1000=" << cc.aniso_mean
                << " std_final=" << ca.aniso_std
                << " std_dt_half=" << cb.aniso_std
                << " std_n1000=" << cc.aniso_std
                << " p_err=" << ca.p_err_max << " k_err=" << ca.k_err_max
                << " n_seeds=" << N_SEEDS << "\n";
        }
    }
    std::cout << "status=" << (iso ? "PASS" : "FAIL") << "\n";
    MPI_Finalize();
    return iso ? 0 : 1;
}
