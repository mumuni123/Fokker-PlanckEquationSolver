#ifndef TAIL_BULK_RETURN_TEST_COMMON_H
#define TAIL_BULK_RETURN_TEST_COMMON_H

#include "background_tail_pic.h"
#include "grid.h"
#include "species.h"
#include "tail_bulk_return.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace tail_return_test {

inline void init_state(int rank, int size, SpatialGrid& grid, Species& bulk,
                       HybridVelocityPartition& partition,
                       BackgroundTailPIC& tail, int nx = 40)
{
    grid.init_with_domain(rank, size, nx, static_cast<double>(nx) * 5.0e-9);
    bulk.init("background_electrons", SpeciesType::BACKGROUND_ELECTRON,
              -Const::qe, Const::me, Param::dens, Param::temperature_e,
              false, grid);
    std::fill(bulk.f.begin(), bulk.f.end(), 0.0);
    partition.init(bulk.cgrid, 6.0, 1.0, 4, 4);
    tail.init(grid);
}

inline std::pair<int, int> safe_velocity_slot(
    const Species& bulk, const HybridVelocityPartition& partition)
{
    const int jc = Param::Nv / 2;
    const int kc = Param::Nmu / 4;
    for (int radius = 0; radius < std::max(Param::Nv, Param::Nmu); ++radius) {
        for (int j = std::max(1, jc - radius);
             j <= std::min(Param::Nv - 2, jc + radius); ++j) {
            for (int k = std::max(1, kc - radius);
                 k <= std::min(Param::Nmu - 2, kc + radius); ++k) {
                const size_t id = idx2(j, k);
                if (partition.bulk_owned_cell[id] != 0 &&
                    partition.kinetic_energy[id] < 5.0e6 * Const::eV) {
                    return std::make_pair(j, k);
                }
            }
        }
    }
    return std::make_pair(-1, -1);
}

inline BackgroundTailParticle make_particle(
    const Species& bulk, int j, int k, double x, double weight,
    std::uint64_t id, std::uint32_t residence = 0)
{
    BackgroundTailParticle p;
    p.x = x;
    p.ux = bulk.cgrid.upar_cells[static_cast<size_t>(j)];
    p.uy = bulk.cgrid.uperp_cells[static_cast<size_t>(k)];
    p.uz = 0.0;
    p.weight = weight;
    p.id = id;
    p.return_residence_steps = residence;
    return p;
}

inline void add_representable_cloud(
    BackgroundTailPIC& tail, const Species& bulk,
    const HybridVelocityPartition& partition, int jc, int kc, double x,
    double weight, std::uint64_t id_base, std::uint32_t residence)
{
    std::uint64_t offset = 0;
    for (int dj = -1; dj <= 1; ++dj) {
        for (int dk = -1; dk <= 1; ++dk) {
            const int j = jc + dj;
            const int k = kc + dk;
            if (j < 0 || j >= Param::Nv || k < 0 || k >= Param::Nmu)
                continue;
            const size_t slot = idx2(j, k);
            if (partition.bulk_owned_cell[slot] == 0 ||
                partition.kinetic_energy[slot] >= 5.0e6 * Const::eV)
                continue;
            tail.particles.push_back(make_particle(
                bulk, j, k, x, weight,
                id_base + offset, residence));
            ++offset;
        }
    }
}

inline bool equal_particles(const std::vector<BackgroundTailParticle>& a,
                            const std::vector<BackgroundTailParticle>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].x != b[i].x || a[i].ux != b[i].ux ||
            a[i].uy != b[i].uy || a[i].uz != b[i].uz ||
            a[i].weight != b[i].weight || a[i].id != b[i].id ||
            a[i].return_residence_steps != b[i].return_residence_steps) {
            return false;
        }
    }
    return true;
}

inline bool equal_doubles(const std::vector<double>& a,
                          const std::vector<double>& b)
{
    return a.size() == b.size() &&
        std::equal(a.begin(), a.end(), b.begin());
}

inline double worst_residual(const TailBulkReturnDiagnostics& d)
{
    return std::max(d.number_residual,
           std::max(d.px_residual,
           std::max(d.jx_residual,
           std::max(d.energy_residual,
           std::max(d.pixx_residual, d.piperp_residual)))));
}

inline double invariant_residual(const TailBulkReturnDiagnostics& d)
{
    return std::max(d.number_residual,
           std::max(d.px_residual, d.energy_residual));
}

inline double representation_residual(const TailBulkReturnDiagnostics& d)
{
    return std::max(d.jx_residual,
           std::max(d.pixx_residual, d.piperp_residual));
}

inline void write_result(const std::string& path,
                         const std::vector<std::pair<std::string, double> >& v,
                         bool pass)
{
    if (path.empty()) return;
    std::ofstream out(path.c_str());
    out << std::setprecision(17);
    for (size_t i = 0; i < v.size(); ++i)
        out << v[i].first << '=' << v[i].second << '\n';
    out << "status=" << (pass ? "PASS" : "FAIL") << '\n';
}

inline bool parse_result_arg(int argc, char** argv, std::string& path)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result" && i + 1 < argc) path = argv[++i];
        else if (arg == "--case" && i + 1 < argc &&
                 std::string(argv[i + 1]) == "all") ++i;
        else return false;
    }
    return true;
}

} // namespace tail_return_test

#endif
