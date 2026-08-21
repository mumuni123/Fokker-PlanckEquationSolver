#include "tail_population_controller.h"
#include "tail_moment_constraint.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace {

struct GroupKey {
    int cell;        // global spatial cell
    int sign;        // u_parallel sign (+1 / -1)
    int upar_bin;
    int energy_bin;
    bool operator<(const GroupKey& other) const
    {
        if (cell != other.cell) return cell < other.cell;
        if (sign != other.sign) return sign < other.sign;
        if (upar_bin != other.upar_bin) return upar_bin < other.upar_bin;
        return energy_bin < other.energy_bin;
    }
};

bool is_boundary_cell(int cell_global, int nx_global)
{
    return cell_global == 0 || cell_global == nx_global - 1;
}

// Deterministic split factor: ceil(w / (max_weight_ratio * w_min)) bounded
// by k_split_cap.  Children share the parent position/velocity, so the
// total weight, position moment and velocity moments are unchanged exactly
// (section 7.10).
const int k_split_cap = 64;

} // namespace

TailPopulationController::TailPopulationController()
{}

void TailPopulationController::configure(const Config& config)
{
    config_ = config;
    if (!config_.enabled) return;
    if (config_.control_interval <= 0) config_.control_interval = 1;
    if (config_.target_particles_per_phase_bin <= 0) {
        config_.target_particles_per_phase_bin = 1;
    }
    if (config_.max_particles_per_phase_bin <
        config_.target_particles_per_phase_bin) {
        config_.max_particles_per_phase_bin =
            config_.target_particles_per_phase_bin;
    }
    if (!(config_.max_weight_ratio > 1.0)) config_.max_weight_ratio = 2.0;
    config_.max_support =
        std::max(1, std::min(7, config_.max_support));
}

bool TailPopulationController::active_step(int step) const
{
    if (!config_.enabled) return false;
    return step % config_.control_interval == 0;
}

TailPopulationController::Diagnostics TailPopulationController::apply(
    BackgroundTailPIC& tail_trial, const SpatialGrid& grid,
    const HybridVelocityPartition& partition, int step, int mpi_rank)
{
    Diagnostics d;
    if (!active_step(step)) return d;
    d.applied = true;
    d.particles_before_local = tail_trial.particles.size();
    if (tail_trial.particles.empty()) {
        d.particles_after_local = 0;
        return d;
    }
    // Section 6.4: split children are new particles and must receive
    // deterministic IDs that cannot collide with any existing particle.
    // The per-rank counter is therefore advanced past the largest low-32-bit
    // value in the local representation before any new ID is issued.  In
    // production the converter already maintains this invariant; the guard
    // makes the controller robust to any externally loaded particle set.
    {
        // GCC on Linux typedefs std::uint64_t as unsigned long while a
        // literal ULL is unsigned long long; keep the mask the same type as
        // max_low so std::max deduces a single template argument.
        const std::uint64_t low_mask = 0xffffffffULL;
        std::uint64_t max_low = 0;
        for (size_t i = 0; i < tail_trial.particles.size(); ++i) {
            max_low = std::max(max_low,
                               tail_trial.particles[i].id & low_mask);
        }
        while (tail_trial.id_counter() <= max_low) {
            tail_trial.next_particle_id(mpi_rank);
        }
    }

    // Group the whole trial representation (existing + newly converted) by
    // (global spatial cell, u_parallel sign, |u_parallel| bin, energy bin):
    // same physical cell, same CIC stencil, adjacent tail phase-space bins
    // (section 7.10).  Boundary cells are skipped below.
    std::map<GroupKey, std::vector<size_t> > groups;
    for (size_t i = 0; i < tail_trial.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail_trial.particles[i];
        const double upar = std::fabs(p.ux);
        const double gamma =
            std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        const double ke = Const::me * Const::c * Const::c * (gamma - 1.0);
        int cell = static_cast<int>(std::floor(p.x / grid.dx));
        if (cell < 0) cell = 0;
        if (cell >= grid.nx_global) cell = grid.nx_global - 1;
        GroupKey key;
        key.cell = cell;
        key.sign = (p.ux >= 0.0) ? 1 : -1;
        key.upar_bin = partition.upar_bin(upar);
        key.energy_bin = partition.energy_bin(ke);
        groups[key].push_back(i);
    }

    const double compression_tolerance = 1.0e-10;
    std::vector<BackgroundTailParticle> out;
    out.reserve(tail_trial.particles.size());
    for (std::map<GroupKey, std::vector<size_t> >::const_iterator
             git = groups.begin(); git != groups.end(); ++git) {
        const GroupKey& key = git->first;
        const std::vector<size_t>& indices = git->second;
        ++d.groups_considered;

        // Boundary cells (and therefore particles about to exit next step)
        // never take part in merging (section 7.10).
        if (is_boundary_cell(key.cell, grid.nx_global)) {
            ++d.boundary_groups_skipped;
            for (size_t c = 0; c < indices.size(); ++c) {
                out.push_back(tail_trial.particles[indices[c]]);
            }
            continue;
        }
        if (indices.size() <=
            static_cast<size_t>(config_.target_particles_per_phase_bin)) {
            ++d.groups_below_target;
            for (size_t c = 0; c < indices.size(); ++c) {
                out.push_back(tail_trial.particles[indices[c]]);
            }
            continue;
        }

        // Deterministic conservative compression (shared 7-moment module)
        // applies only above max_particles_per_phase_bin (section 7.10);
        // a group between target and max is still routed through the
        // weight-ratio split below.
        std::vector<BackgroundTailParticle> group_out;
        if (indices.size() >
            static_cast<size_t>(config_.max_particles_per_phase_bin)) {
            std::vector<std::vector<double> > cols;
            std::vector<double> w;
            std::vector<double> refv(7, 0.0);
            TailMoment7 ref;
            for (size_t c = 0; c < indices.size(); ++c) {
                const BackgroundTailParticle& p =
                    tail_trial.particles[indices[c]];
                TailMoment7 col;
                tail_particle_moments(1.0, p.x, p.ux, p.uy, p.uz, col);
                std::vector<double> colv(7, 0.0);
                colv[0] = col.n;
                colv[1] = col.px;
                colv[2] = col.jx;
                colv[3] = col.ke;
                colv[4] = col.pixx;
                colv[5] = col.piperp;
                colv[6] = col.xw;
                cols.push_back(colv);
                w.push_back(p.weight);
                TailMoment7 pm;
                tail_particle_moments(p.weight, p.x, p.ux, p.uy, p.uz, pm);
                ref.n += pm.n;
                ref.px += pm.px;
                ref.jx += pm.jx;
                ref.ke += pm.ke;
                ref.pixx += pm.pixx;
                ref.piperp += pm.piperp;
                ref.xw += pm.xw;
            }
            refv[0] = ref.n;
            refv[1] = ref.px;
            refv[2] = ref.jx;
            refv[3] = ref.ke;
            refv[4] = ref.pixx;
            refv[5] = ref.piperp;
            refv[6] = ref.xw;
            std::vector<double> w_out = w;
            const bool compressed = tail_compress_moment_supports(
                cols, w_out, refv,
                static_cast<size_t>(config_.max_support),
                compression_tolerance);
            if (!compressed) {
                // Keep the original particles; never accept a low-order
                // match.
                ++d.compression_fallback_count;
                for (size_t c = 0; c < indices.size(); ++c) {
                    out.push_back(tail_trial.particles[indices[c]]);
                }
                continue;
            }
            size_t support_count = 0;
            for (size_t q = 0; q < w_out.size(); ++q) {
                if (w_out[q] > 0.0) ++support_count;
            }
            if (support_count < indices.size()) ++d.groups_compressed;

            // Residual report (section 7.10: N/Px/Jx/K/Pixx/Piperp/Xw).
            TailMoment7 got;
            for (size_t q = 0; q < indices.size(); ++q) {
                if (!(w_out[q] > 0.0)) continue;
                const BackgroundTailParticle& p =
                    tail_trial.particles[indices[q]];
                TailMoment7 m;
                tail_particle_moments(w_out[q], p.x, p.ux, p.uy, p.uz, m);
                got.n += m.n;
                got.px += m.px;
                got.jx += m.jx;
                got.ke += m.ke;
                got.pixx += m.pixx;
                got.piperp += m.piperp;
                got.xw += m.xw;
            }
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
            for (int r = 0; r < 7; ++r) {
                d.max_residual[r] = std::max(d.max_residual[r], residuals[r]);
            }

            // Emit the merged supports (same positions/velocities, merged
            // weights; the kept particle ID stays with its support).
            for (size_t q = 0; q < indices.size(); ++q) {
                if (!(w_out[q] > 0.0)) continue;
                BackgroundTailParticle p =
                    tail_trial.particles[indices[q]];
                p.weight = w_out[q];
                group_out.push_back(p);
            }
        } else {
            for (size_t c = 0; c < indices.size(); ++c) {
                group_out.push_back(tail_trial.particles[indices[c]]);
            }
        }

        // Weight-ratio splitting: heavy macro particles are split into
        // near-equal children so future merging does not skew the noise.
        // Children keep the parent position/velocity; total weight, position
        // moment and velocity moments are unchanged exactly.
        double w_min = -1.0;
        for (size_t c = 0; c < group_out.size(); ++c) {
            const double w = group_out[c].weight;
            if (w_min < 0.0 || w < w_min) w_min = w;
        }
        if (config_.max_weight_ratio > 0.0 && w_min > 0.0) {
            bool split_any = false;
            std::vector<BackgroundTailParticle> split_out;
            for (size_t c = 0; c < group_out.size(); ++c) {
                const double w = group_out[c].weight;
                const double ratio = w / w_min;
                if (ratio <= config_.max_weight_ratio) {
                    split_out.push_back(group_out[c]);
                    continue;
                }
                int k = static_cast<int>(std::ceil(
                    ratio / config_.max_weight_ratio));
                if (k < 2) k = 2;
                if (k > k_split_cap) {
                    // Keep the parent unsplit: the split budget is bounded
                    // (section 18.3); the group is counted for the report.
                    ++d.split_capped_count;
                    split_out.push_back(group_out[c]);
                    continue;
                }
                split_any = true;
                const double child_weight = w / static_cast<double>(k);
                for (int a = 0; a < k; ++a) {
                    BackgroundTailParticle child = group_out[c];
                    child.weight = child_weight;
                    child.id = tail_trial.next_particle_id(mpi_rank);
                    // A split child is a new macro-particle representation.
                    // It must earn a fresh H10 residence history rather than
                    // inheriting a parent that may already be return-eligible.
                    child.return_residence_steps = 0;
                    split_out.push_back(child);
                }
            }
            if (split_any) ++d.groups_split;
            group_out.swap(split_out);
        }
        for (size_t c = 0; c < group_out.size(); ++c) {
            out.push_back(group_out[c]);
        }
    }
    tail_trial.particles.swap(out);
    d.particles_after_local = tail_trial.particles.size();
    return d;
}
