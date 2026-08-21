#include "bulk_tail_converter.h"
#include "tail_moment_constraint.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace {

const double k_hard_gate = 1.0e-10;   // N/Px/K relative residual gate
const double k_fidelity_gate = 1.0e-9; // Jx/Pixx/Piperp relative gate

double interface_roundoff_floor(const SpatialGrid& grid)
{
    // A flux batch is measured in particle number per transverse area.  Use
    // the same grid-scaled floor as static extraction: this is far below any
    // resolved conversion, while rejecting subnormal collision face residues.
    return 64.0 * std::numeric_limits<double>::epsilon() *
           std::max(1.0, Param::dens * grid.dx);
}

} // namespace

BulkTailConverter::BulkTailConverter()
    // Section 7.11.9 demonstrated geometric support, but the section 7.11.10
    // six-moment feasibility gate has not passed.  Keep the candidate out of
    // production; dedicated tests enable it explicitly.
    : loading_policy_(
          BulkTailLoadingPolicy::THRESHOLD_AWARE_COMPRESSION),
      subcell_loading_enabled_(false),
      moment_audit_enabled_(false),
      moment_audit_top_cell_count_(64),
      subcell_geometry_threshold_(0.0),
      subcell_geometry_valid_(false)
{}

void BulkTailConverter::prepare_subcell_geometry(
    const CylindricalVelocityGrid& cgrid, double min_conversion_energy)
{
    if (subcell_geometry_valid_ &&
        subcell_geometry_.size() == Param::Nvmu &&
        subcell_geometry_threshold_ == min_conversion_energy) {
        return;
    }

    subcell_geometry_.assign(Param::Nvmu, std::vector<SubcellGeometry>());
    for (int iv = 0; iv < Param::Nv; ++iv) {
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const size_t id = idx2(iv, imu);
            const std::vector<TailSubcellNode> nodes =
                TailSubcellQuadrature::nodes(cgrid, iv, imu);
            for (size_t q = 0; q < nodes.size(); ++q) {
                if (nodes[q].kinetic_energy < min_conversion_energy) {
                    continue;
                }
                SubcellGeometry geometry;
                geometry.upar = nodes[q].upar;
                geometry.uperp = nodes[q].uperp;
                geometry.mass_fraction = nodes[q].mass_fraction;
                mass_cell_moments(1.0, geometry.upar, geometry.uperp,
                                  geometry.column.n, geometry.column.px,
                                  geometry.column.ke, geometry.column.jx,
                                  geometry.column.pixx,
                                  geometry.column.piperp);
                subcell_geometry_[id].push_back(geometry);
            }
        }
    }
    subcell_geometry_threshold_ = min_conversion_energy;
    subcell_geometry_valid_ = true;
}

bool BulkTailConverter::compress_quartets(
    const std::vector<Moment6>& cols, std::vector<double>& w,
    const Moment6& ref, size_t max_support, double tolerance)
{
    const size_t n = w.size();
    if (n <= max_support) return true;
    std::vector<std::vector<double> > matrix(n, std::vector<double>(6, 0.0));
    std::vector<double> refv(6, 0.0);
    for (size_t q = 0; q < n; ++q) {
        matrix[q][0] = cols[q].n;
        matrix[q][1] = cols[q].px;
        matrix[q][2] = cols[q].jx;
        matrix[q][3] = cols[q].ke;
        matrix[q][4] = cols[q].pixx;
        matrix[q][5] = cols[q].piperp;
    }
    refv[0] = ref.n;
    refv[1] = ref.px;
    refv[2] = ref.jx;
    refv[3] = ref.ke;
    refv[4] = ref.pixx;
    refv[5] = ref.piperp;
    // Section 14.5.1: the compression math lives exactly once in the shared
    // TailMomentConstraint module; the converter only adapts its six-moment
    // layout.
    return tail_compress_moment_supports(matrix, w, refv, max_support,
                                         tolerance);
}

BulkTailConversionDiagnostics BulkTailConverter::extract_after_substep(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid, const HybridVelocityPartition& partition,
    int accepted_step, ConversionLocation location, int mpi_rank)
{
    BulkTailConversionDiagnostics d;
    d.accepted_step = accepted_step;
    d.location = location;
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const double inv_dx = 1.0 / grid.dx;
    const size_t nslots = static_cast<size_t>(Param::Nvmu);
    // Physically-anchored roundoff floor: after a real remap step the
    // conversion cells can carry ~1e-16 relative residues (positive or
    // negative).  Section 7.8 requires the remap's roundoff path to resolve
    // them before the converter; as long as the remap only tolerates them,
    // the converter treats masses within this floor as roundoff (skipped)
    // and hard-fails only on material negatives.
    const double roundoff_floor = interface_roundoff_floor(grid);
    // Stage-H6 conversion source spectrum (section 13.2): log-spaced bins
    // over the conversion energy window, filled below during the scan.
    const int source_bins = 64;
    const double ke_min = partition.min_conversion_energy;
    const double ke_max = partition.max_conversion_energy;
    const double source_log_span =
        (ke_max > ke_min) ? std::log(ke_max / ke_min) : 1.0;
    d.conversion_source_energy_edges.resize(static_cast<size_t>(source_bins + 1));
    for (int b = 0; b <= source_bins; ++b) {
        const double t = static_cast<double>(b) / source_bins;
        d.conversion_source_energy_edges[static_cast<size_t>(b)] =
            (ke_max > ke_min) ? ke_min * std::exp(source_log_span * t) : ke_min;
    }
    d.pre_extraction_bulk_number_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    d.pre_extraction_bulk_energy_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    d.removed_bulk_number_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    d.removed_bulk_energy_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    d.created_tail_number_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    d.created_tail_energy_spectrum.assign(static_cast<size_t>(source_bins), 0.0);
    if (subcell_loading_enabled_) {
        prepare_subcell_geometry(bulk_trial.cgrid,
                                 partition.min_conversion_energy);
    }

    // ---- scan conversion cells (sections 7.2/7.8) ----
    std::vector<ConversionMassRequest> requests;
    std::map<GroupKey, std::vector<CellMassRef> > groups;
    std::vector<double> conv_density_before(static_cast<size_t>(nxl), 0.0);
    for (int il = 0; il < nxl; ++il) {
        const int ix_global = grid.ix_start + il;
        const size_t xbase = static_cast<size_t>(ng + il) * nslots;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                if (!partition.is_conversion(j, k)) continue;
                const double mass =
                    bulk_trial.f[xbase + idx2(j, k)];
                if (!std::isfinite(mass) || mass < -roundoff_floor) {
                    // Section 7.8: unphysical bulk mass in the conversion
                    // region is a hard failure; nothing may be modified.
                    d.complete = false;
                    d.conservative = false;
                    d.finite = false;
                    return d;
                }
                if (mass <= roundoff_floor) continue;
                ConversionMassRequest req;
                req.ix_global = ix_global;
                req.iv = j;
                req.imu = k;
                req.mass = mass;
                requests.push_back(req);
                conv_density_before[static_cast<size_t>(il)] +=
                    mass * inv_dx;

                const double upar = bulk_trial.cgrid.upar_cells[j];
                const double ke = partition.kinetic_energy[idx2(j, k)];
                int source_bin =
                    (ke > ke_min && ke_max > ke_min)
                        ? static_cast<int>(std::log(ke / ke_min) /
                                           source_log_span * source_bins)
                        : 0;
                if (source_bin >= source_bins) source_bin = source_bins - 1;
                const size_t sb = static_cast<size_t>(source_bin);
                d.pre_extraction_bulk_number_spectrum[sb] += mass;
                d.pre_extraction_bulk_energy_spectrum[sb] += mass * ke;
                // Extraction removes exactly the scanned conversion mass.
                d.removed_bulk_number_spectrum[sb] += mass;
                d.removed_bulk_energy_spectrum[sb] += mass * ke;
                GroupKey key;
                key.ix_global = ix_global;
                key.sign = (upar >= 0.0) ? 1 : -1;
                key.upar_bin = partition.upar_bin(upar);
                // Section 7.11.3: the loading policy selects the energy
                // grouping.  The golden path uses the same grouping as the
                // current reference path (compression is skipped below, so
                // the grouping does not change the golden output); the
                // production path bins by the explicit
                // conversion_energy_edges (section 7.11.4 branch B) so no
                // compression group may span a threshold audit bin.
                key.energy_bin =
                    (loading_policy_ ==
                             BulkTailLoadingPolicy::THRESHOLD_AWARE_COMPRESSION)
                        ? partition.energy_bin_threshold_aware(ke)
                        : partition.energy_bin_uniform(ke);
                CellMassRef ref;
                ref.iv = j;
                ref.imu = k;
                ref.mass = mass;
                groups[key].push_back(ref);
            }
        }
    }
    if (requests.empty()) {
        d.complete = true;
        d.conservative = true;
        d.fidelity_ok = true;
        d.finite = true;
        return d;
    }

    // Requests are complete but the source mass has not yet been changed.
    // This read-only audit must never affect the conversion transaction.
    if (moment_audit_enabled_) {
        std::vector<BulkTailMomentAuditRequest> audit_requests;
        audit_requests.reserve(requests.size());
        for (size_t q = 0; q < requests.size(); ++q) {
            BulkTailMomentAuditRequest request;
            request.ix_global = requests[q].ix_global;
            request.iv = requests[q].iv;
            request.imu = requests[q].imu;
            request.mass = requests[q].mass;
            audit_requests.push_back(request);
        }
        d.moment_audit = bulk_tail_audit_conversion_requests(
            bulk_trial.cgrid, partition, audit_requests, mpi_rank,
            moment_audit_top_cell_count_);
    } else {
        d.moment_audit.enabled = false;
    }

    // ---- remove the mass from the bulk (single production extraction) ----
    const SpeciesConversionResult extraction =
        bulk_trial.extract_conversion_masses(requests);
    if (!extraction.valid || !extraction.complete) {
        d.complete = false;
        d.conservative = false;
        d.finite = true;
        return d;
    }
    d.number_removed = extraction.number_removed;
    d.px_removed = extraction.px_removed;
    d.energy_removed = extraction.energy_removed;
    d.jx_dx_removed = extraction.jx_dx_removed;
    d.pixx_dx_removed = extraction.pixx_dx_removed;
    d.piperp_dx_removed = extraction.piperp_dx_removed;
    d.number_scale = extraction.number_scale;
    d.px_scale = extraction.px_scale;
    d.energy_scale = extraction.energy_scale;
    d.jx_scale = extraction.jx_scale;
    d.pixx_scale = extraction.pixx_scale;
    d.piperp_scale = extraction.piperp_scale;

    // ---- materialise per group (reference quartet + compression) ----
    std::vector<double> conv_density_after(static_cast<size_t>(nxl), 0.0);
    for (typename std::map<GroupKey, std::vector<CellMassRef> >::iterator
             git = groups.begin(); git != groups.end(); ++git) {
        const GroupKey& key = git->first;
        const std::vector<CellMassRef>& cells = git->second;
        struct Support {
            double upar;
            double uperp;
            double prior_weight;
            Moment6 column;
        };
        std::vector<Support> supports;
        std::vector<Moment6> cols;
        std::vector<double> w;
        Moment6 ref;
        for (size_t c = 0; c < cells.size(); ++c) {
            const CellMassRef& cr = cells[c];
            const double upar = bulk_trial.cgrid.upar_cells[cr.iv];
            const double uperp = bulk_trial.cgrid.uperp_cells[cr.imu];
            double number, px, energy, jx_dx, pixx_dx, piperp_dx;
            mass_cell_moments(cr.mass, upar, uperp, number, px, energy,
                              jx_dx, pixx_dx, piperp_dx);
            ref.n += number;
            ref.px += px;
            ref.ke += energy;
            ref.jx += jx_dx;
            ref.pixx += pixx_dx;
            ref.piperp += piperp_dx;
            if (subcell_loading_enabled_) {
                const std::vector<SubcellGeometry>& geometry =
                    subcell_geometry_[idx2(cr.iv, cr.imu)];
                for (size_t q = 0; q < geometry.size(); ++q) {
                    Support support;
                    support.upar = geometry[q].upar;
                    support.uperp = geometry[q].uperp;
                    support.prior_weight = cr.mass * geometry[q].mass_fraction;
                    support.column = geometry[q].column;
                    supports.push_back(support);
                }
            }
            if (!subcell_loading_enabled_) {
                Support support;
                support.upar = upar;
                support.uperp = uperp;
                support.prior_weight = cr.mass;
                mass_cell_moments(1.0, support.upar, support.uperp,
                                  support.column.n, support.column.px,
                                  support.column.ke, support.column.jx,
                                  support.column.pixx,
                                  support.column.piperp);
                supports.push_back(support);
            }
        }

        for (size_t q = 0; q < supports.size(); ++q) {
            cols.push_back(supports[q].column);
            w.push_back(supports[q].prior_weight);
        }

        bool subcell_fit = !subcell_loading_enabled_;
        const bool subcell_support_missing = cols.empty();
        if (subcell_loading_enabled_ && !cols.empty()) {
            std::vector<std::vector<double> > matrix(cols.size(),
                                                       std::vector<double>(6));
            std::vector<double> refv(6);
            for (size_t q = 0; q < cols.size(); ++q) {
                matrix[q][0] = cols[q].n; matrix[q][1] = cols[q].px;
                matrix[q][2] = cols[q].jx; matrix[q][3] = cols[q].ke;
                matrix[q][4] = cols[q].pixx; matrix[q][5] = cols[q].piperp;
            }
            refv[0] = ref.n; refv[1] = ref.px; refv[2] = ref.jx;
            refv[3] = ref.ke; refv[4] = ref.pixx; refv[5] = ref.piperp;
            const std::vector<double> prior = w;
            subcell_fit = tail_solve_nonnegative_moment_weights(
                matrix, refv, prior, w, 1.0e-10);
        }
        if (!subcell_fit) {
            // Explicit, auditable fallback.  It is conservative but a
            // nonzero counter keeps the subcell-loading stage gate closed.
            supports.clear(); cols.clear(); w.clear();
            for (size_t c = 0; c < cells.size(); ++c) {
                const CellMassRef& cr = cells[c];
                Support support;
                support.upar = bulk_trial.cgrid.upar_cells[cr.iv];
                support.uperp = bulk_trial.cgrid.uperp_cells[cr.imu];
                support.prior_weight = cr.mass;
                mass_cell_moments(1.0, support.upar, support.uperp,
                                  support.column.n, support.column.px,
                                  support.column.ke, support.column.jx,
                                  support.column.pixx,
                                  support.column.piperp);
                supports.push_back(support);
                cols.push_back(support.column); w.push_back(cr.mass);
            }
            ++d.subcell_fallback_count;
            BulkTailConversionDiagnostics::SubcellFallbackRecord record;
            record.ix_global = key.ix_global;
            record.iv = cells.front().iv;
            record.imu = cells.front().imu;
            record.reason = subcell_support_missing ? 1 : 2;
            record.fallback_particles =
                static_cast<std::uint64_t>(4 * cells.size());
            record.number_target = ref.n;
            record.px_target = ref.px;
            record.jx_target = ref.jx;
            record.energy_target = ref.ke;
            record.pixx_target = ref.pixx;
            record.piperp_target = ref.piperp;
            d.subcell_fallbacks.push_back(record);
        }

        bool compressed = true;
        if (loading_policy_ !=
            BulkTailLoadingPolicy::GOLDEN_QUARTETS_NO_COMPRESSION) {
            if (cols.size() > 7) {
                compressed = compress_quartets(cols, w, ref, 7, 1.0e-10);
            }
            if (!compressed) {
                // Fallback: keep the golden reference quartets (section
                // 7.5).
                supports.clear(); cols.clear(); w.clear();
                for (size_t c = 0; c < cells.size(); ++c) {
                    const CellMassRef& cr = cells[c];
                    Support support;
                    support.upar = bulk_trial.cgrid.upar_cells[cr.iv];
                    support.uperp = bulk_trial.cgrid.uperp_cells[cr.imu];
                    support.prior_weight = cr.mass;
                    mass_cell_moments(1.0, support.upar, support.uperp,
                                      support.column.n, support.column.px,
                                      support.column.ke, support.column.jx,
                                      support.column.pixx,
                                      support.column.piperp);
                    supports.push_back(support);
                    cols.push_back(support.column); w.push_back(cr.mass);
                }
                ++d.compression_fallback_count;
                if (subcell_loading_enabled_) {
                    ++d.subcell_fallback_count;
                    BulkTailConversionDiagnostics::SubcellFallbackRecord record;
                    record.ix_global = key.ix_global;
                    record.iv = cells.front().iv;
                    record.imu = cells.front().imu;
                    record.reason = 3;
                    record.fallback_particles =
                        static_cast<std::uint64_t>(4 * cells.size());
                    record.number_target = ref.n;
                    record.px_target = ref.px;
                    record.jx_target = ref.jx;
                    record.energy_target = ref.ke;
                    record.pixx_target = ref.pixx;
                    record.piperp_target = ref.piperp;
                    d.subcell_fallbacks.push_back(record);
                }
            }
        }

        if (subcell_loading_enabled_ && subcell_fit && compressed) {
            d.subcell_cells_loaded += cells.size();
            d.subcell_support_count += supports.size();
        }

        const double x_center =
            (static_cast<double>(key.ix_global) + 0.5) * grid.dx;
        const int local_il = key.ix_global - grid.ix_start;
        for (size_t q = 0; q < cols.size(); ++q) {
            if (!(w[q] > 0.0)) continue;
            const double upar = supports[q].upar;
            const double uperp = supports[q].uperp;
            const double phi[4] = { 0.0, 0.5 * Const::pi,
                                    Const::pi, 1.5 * Const::pi };
            for (int a = 0; a < 4; ++a) {
                BackgroundTailParticle p;
                p.x = x_center;
                p.ux = upar;
                p.uy = uperp * std::cos(phi[a]);
                p.uz = uperp * std::sin(phi[a]);
                p.weight = w[q] * 0.25;
                p.id = tail_trial.next_particle_id(mpi_rank);
                p.return_residence_steps = 0;
                tail_trial.particles.push_back(p);
                ++d.particles_created;
                const double gamma = std::sqrt(
                    1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
                d.number_created += p.weight;
                d.px_created += Const::me * Const::c * p.weight * p.ux;
                d.energy_created += Const::me * Const::c * Const::c *
                                    p.weight * (gamma - 1.0);
                d.jx_dx_created +=
                    -Const::qe * Const::c * p.weight * p.ux / gamma;
                d.pixx_dx_created += Const::me * Const::c * Const::c *
                                     p.weight * p.ux * p.ux / gamma;
                d.piperp_dx_created += Const::me * Const::c * Const::c *
                    p.weight * (p.uy * p.uy + p.uz * p.uz) / gamma;
                const double pke = Const::me * Const::c * Const::c *
                                   (gamma - 1.0);
                int created_bin =
                    (pke > ke_min && ke_max > ke_min)
                        ? static_cast<int>(std::log(pke / ke_min) /
                                           source_log_span * source_bins)
                        : 0;
                if (created_bin < 0) created_bin = 0;
                if (created_bin >= source_bins) created_bin = source_bins - 1;
                d.created_tail_number_spectrum[static_cast<size_t>(created_bin)] +=
                    p.weight;
                d.created_tail_energy_spectrum[static_cast<size_t>(created_bin)] +=
                    p.weight * pke;
            }
            if (local_il >= 0 && local_il < nxl) {
                conv_density_after[static_cast<size_t>(local_il)] +=
                    w[q] * inv_dx;
            }
        }
    }

    // ---- post-conversion bulk scan (section 7.8) ----
    const double residue_tolerance = std::max(
        64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, d.number_removed),
        roundoff_floor);
    double residue_max = 0.0;
    for (int il = 0; il < nxl; ++il) {
        const size_t xbase = static_cast<size_t>(ng + il) * nslots;
        for (int j = 0; j < Param::Nv; ++j) {
            for (int k = 0; k < Param::Nmu; ++k) {
                if (!partition.is_conversion(j, k)) continue;
                residue_max = std::max(
                    residue_max,
                    std::fabs(bulk_trial.f[xbase + idx2(j, k)]));
            }
        }
    }
    d.complete = (residue_max <= residue_tolerance);

    // ---- refresh bulk moments (section 7.8) ----
    bulk_trial.compute_moments();

    // ---- density diagnostics (section 7.7) ----
    double l2_before = 0.0;
    double l2_diff = 0.0;
    double linf_diff = 0.0;
    for (int il = 0; il < nxl; ++il) {
        const double before = conv_density_before[static_cast<size_t>(il)];
        const double after = conv_density_after[static_cast<size_t>(il)];
        const double diff = before - after;
        l2_before += before * before;
        l2_diff += diff * diff;
        linf_diff = std::max(linf_diff, std::fabs(diff));
    }
    d.removed_bulk_number_by_x.assign(static_cast<size_t>(nxl), 0.0);
    d.created_tail_number_by_x.assign(static_cast<size_t>(nxl), 0.0);
    for (int il = 0; il < nxl; ++il) {
        d.removed_bulk_number_by_x[static_cast<size_t>(il)] =
            conv_density_before[static_cast<size_t>(il)] * grid.dx;
        d.created_tail_number_by_x[static_cast<size_t>(il)] =
            conv_density_after[static_cast<size_t>(il)] * grid.dx;
    }
    d.rho_l2_before_after =
        std::sqrt(l2_diff) / std::max(1.0, std::sqrt(l2_before));
    d.rho_linf_before_after = linf_diff;

    // ---- residual and gate flags ----
    // Normalise by the L1 scale of the removed contributions: the signed
    // px/jx totals can cancel between the u_parallel signs, so a naive
    // |signed total| denominator would amplify cancellation-level errors.
    d.number_residual_rel =
        std::fabs(d.number_created - d.number_removed) /
        std::max(1.0, d.number_scale);
    d.px_residual_rel =
        std::fabs(d.px_created - d.px_removed) /
        std::max(1.0, d.px_scale);
    d.energy_residual_rel =
        std::fabs(d.energy_created - d.energy_removed) /
        std::max(1.0, d.energy_scale);
    d.jx_residual_rel =
        std::fabs(d.jx_dx_created - d.jx_dx_removed) /
        std::max(1.0, d.jx_scale);
    d.pixx_residual_rel =
        std::fabs(d.pixx_dx_created - d.pixx_dx_removed) /
        std::max(1.0, d.pixx_scale);
    d.piperp_residual_rel =
        std::fabs(d.piperp_dx_created - d.piperp_dx_removed) /
        std::max(1.0, d.piperp_scale);
    d.conservative =
        d.number_residual_rel <= k_hard_gate &&
        d.px_residual_rel <= k_hard_gate &&
        d.energy_residual_rel <= k_hard_gate;
    d.fidelity_ok =
        d.jx_residual_rel <= k_fidelity_gate &&
        d.pixx_residual_rel <= k_fidelity_gate &&
        d.piperp_residual_rel <= k_fidelity_gate;
    d.finite = std::isfinite(d.number_created) &&
               std::isfinite(d.px_created) &&
               std::isfinite(d.energy_created) &&
               std::isfinite(d.jx_dx_created) &&
               std::isfinite(d.pixx_dx_created) &&
               std::isfinite(d.piperp_dx_created);
    return d;
}

BulkTailConversionDiagnostics BulkTailConverter::convert_flux_batch(
    const BulkTailFluxBatch& batch, BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid, const HybridVelocityPartition& partition,
    int accepted_step, ConversionLocation location, int mpi_rank,
    size_t max_supports)
{
    BulkTailConversionDiagnostics d;
    d.accepted_step = accepted_step;
    d.location = location;
    d.finite = batch.finite;
    d.conservative = false;
    d.fidelity_ok = false;
    d.complete = false;
    // Per-cell bulk removal / tail creation source (section I2): the flux-
    // interface sink clears the tail-owned Eulerian cells during
    // advect_u_parallel, so the continuity audit must account for that mass
    // here.  face_number is the authoritative bulk removal; the parcel node
    // sum is the created tail mass (equal by the face-ledger contract).
    d.removed_bulk_number_by_x.assign(static_cast<size_t>(grid.nx_local), 0.0);
    d.created_tail_number_by_x.assign(static_cast<size_t>(grid.nx_local), 0.0);
    // A duplicate key means that the same shared interface event would be
    // removed more than once.  The batch validator detects this before the
    // transactional loader is entered; never try to compensate by changing
    // the particle weights.
    if (!batch.finite || !batch.nonnegative || batch.duplicate_count != 0 ||
        max_supports == 0) return d;
    const double roundoff_floor = interface_roundoff_floor(grid);

    struct GroupKey {
        int ix;
        int stage;
        int sign;
        int energy_bin;
        bool operator<(const GroupKey& other) const
        {
            if (ix != other.ix) return ix < other.ix;
            if (stage != other.stage) return stage < other.stage;
            if (sign != other.sign) return sign < other.sign;
            return energy_bin < other.energy_bin;
        }
    };
    struct Node {
        double upar;
        double uperp;
        double mass;
    };
    std::map<GroupKey, std::vector<Node> > groups;
    Moment6 reference;
    double reference_scale[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::set<std::uint64_t> ids;
    for (size_t p = 0; p < batch.parcels.size(); ++p) {
        const BulkTailFluxParcel& parcel = batch.parcels[p];
        if (parcel.ix_local < 0 || parcel.ix_local >= grid.nx_local ||
            parcel.ix_global < 0 || parcel.ix_global >= grid.nx_global ||
            (parcel.direction != VelocityFaceDirection::U_PARALLEL &&
             parcel.direction != VelocityFaceDirection::U_PERP)) {
            d.finite = false;
            return d;
        }
        double parcel_number = 0.0;
        for (size_t q = 0; q < parcel.nodes.size(); ++q)
            parcel_number += parcel.nodes[q].mass;
        const double face_scale = std::max(
            1.0, std::max(std::fabs(parcel.face_number),
                          std::fabs(parcel_number)));
        if (!std::isfinite(parcel.face_number) || parcel.face_number < 0.0 ||
            std::fabs(parcel_number - parcel.face_number) >
                1024.0 * std::numeric_limits<double>::epsilon() * face_scale) {
            ++d.face_ledger_mismatch_count;
            d.finite = false;
            return d;
        }
        // Collision face fluxes can leave a positive denormal residue after
        // the one-way interface clipping.  It is below the representable
        // particle budget of one grid cell, and must not be expanded into a
        // many-support tail packet.  Validate the parcel ledger first, then
        // discard the whole packet so no partial moment is retained.
        if (parcel.face_number <= roundoff_floor) {
            d.roundoff_discarded_number += parcel.face_number;
            continue;
        }
        d.removed_bulk_number_by_x[static_cast<size_t>(parcel.ix_local)] +=
            parcel.face_number;
        d.created_tail_number_by_x[static_cast<size_t>(parcel.ix_local)] +=
            parcel_number;
        for (size_t q = 0; q < parcel.nodes.size(); ++q) {
            const FluxParcelNode& node = parcel.nodes[q];
            const double gamma = std::sqrt(1.0 + node.upar * node.upar +
                                            node.uperp * node.uperp);
            const double ke = Const::me * Const::c * Const::c * (gamma - 1.0);
            const double debt = 256.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, partition.min_conversion_energy);
            if (!std::isfinite(node.upar) || !std::isfinite(node.uperp) ||
                !std::isfinite(node.mass) || node.mass < 0.0 ||
                ke < partition.min_conversion_energy - debt) {
                d.finite = false;
                return d;
            }
            double n = 0.0, px = 0.0, energy = 0.0, jx = 0.0;
            double pixx = 0.0, piperp = 0.0;
            mass_cell_moments(node.mass, node.upar, node.uperp,
                              n, px, energy, jx, pixx, piperp);
            reference.n += n; reference.px += px; reference.ke += energy;
            reference.jx += jx; reference.pixx += pixx;
            reference.piperp += piperp;
            reference_scale[0] += std::fabs(n);
            reference_scale[1] += std::fabs(px);
            reference_scale[2] += std::fabs(jx);
            reference_scale[3] += std::fabs(energy);
            reference_scale[4] += std::fabs(pixx);
            reference_scale[5] += std::fabs(piperp);
            GroupKey key;
            key.ix = parcel.ix_global;
            key.stage = parcel.operator_stage;
            // The u_parallel sign is a physical grouping variable.  A
            // u_perp-face parcel has a nonnegative transverse face speed;
            // keep its two transverse quadrature nodes in one group.
            key.sign = parcel.direction == VelocityFaceDirection::U_PARALLEL &&
                       node.upar < 0.0 ? -1 : 1;
            key.energy_bin = partition.energy_bin(ke);
            groups[key].push_back(Node{node.upar, node.uperp, node.mass});
        }
    }
    d.number_removed = reference.n;
    d.px_removed = reference.px;
    d.energy_removed = reference.ke;
    d.jx_dx_removed = reference.jx;
    d.pixx_dx_removed = reference.pixx;
    d.piperp_dx_removed = reference.piperp;
    d.number_scale = reference_scale[0];
    d.px_scale = reference_scale[1];
    d.energy_scale = reference_scale[3];
    d.jx_scale = reference_scale[2];
    d.pixx_scale = reference_scale[4];
    d.piperp_scale = reference_scale[5];

    BackgroundTailPIC candidate = tail_trial;
    // Existing accepted particles are part of the transactional state.  The
    // generated-id check must therefore cover both the candidate additions
    // and the already accepted particles, not just ids created in this call.
    for (size_t p = 0; p < candidate.particles.size(); ++p) {
        if (!std::isfinite(candidate.particles[p].weight) ||
            !(candidate.particles[p].weight > 0.0) ||
            !ids.insert(candidate.particles[p].id).second) {
            d.finite = false;
            ++d.duplicate_id_count;
            return d;
        }
    }
    for (std::map<GroupKey, std::vector<Node> >::const_iterator it = groups.begin();
         it != groups.end(); ++it) {
        const std::vector<Node>& nodes = it->second;
        std::vector<Moment6> cols;
        std::vector<double> weights;
        cols.reserve(nodes.size());
        weights.reserve(nodes.size());
        Moment6 ref;
        for (size_t q = 0; q < nodes.size(); ++q) {
            Moment6 col;
            mass_cell_moments(1.0, nodes[q].upar, nodes[q].uperp,
                              col.n, col.px, col.ke, col.jx,
                              col.pixx, col.piperp);
            cols.push_back(col);
            weights.push_back(nodes[q].mass);
            ref.n += nodes[q].mass * col.n;
            ref.px += nodes[q].mass * col.px;
            ref.ke += nodes[q].mass * col.ke;
            ref.jx += nodes[q].mass * col.jx;
            ref.pixx += nodes[q].mass * col.pixx;
            ref.piperp += nodes[q].mass * col.piperp;
        }
        const bool compressed = compress_quartets(
            cols, weights, ref, max_supports, 1.0e-10);
        if (!compressed) ++d.compression_fallback_count;
        size_t active_supports = 0;
        for (size_t q = 0; q < weights.size(); ++q) {
            if (weights[q] > 0.0) ++active_supports;
        }
        if (active_supports > max_supports) {
            ++d.support_limit_violation_count;
            d.finite = false;
            return d;
        }
        const int local_ix = it->first.ix - grid.ix_start;
        if (local_ix < 0 || local_ix >= grid.nx_local) {
            d.finite = false;
            return d;
        }
        const double x = grid.x(grid.nghost + local_ix);
        for (size_t q = 0; q < weights.size(); ++q) {
            if (!(weights[q] >= 0.0) || !std::isfinite(weights[q]) ||
                weights[q] == 0.0) continue;
            const double phi[4] = { 0.0, 0.5 * Const::pi,
                                    Const::pi, 1.5 * Const::pi };
            for (int a = 0; a < 4; ++a) {
                BackgroundTailParticle p;
                p.x = x;
                p.ux = nodes[q].upar;
                p.uy = nodes[q].uperp * std::cos(phi[a]);
                p.uz = nodes[q].uperp * std::sin(phi[a]);
                p.weight = weights[q] / 4.0;
                if (!std::isfinite(p.weight) || p.weight < 0.0) {
                    d.finite = false;
                    return d;
                }
                p.id = candidate.next_particle_id(mpi_rank);
                p.return_residence_steps = 0;
                if (!ids.insert(p.id).second) {
                    d.finite = false;
                    ++d.duplicate_id_count;
                    return d;
                }
                candidate.particles.push_back(p);
                ++d.particles_created;
                d.number_created += p.weight;
                const double gamma = std::sqrt(1.0 + p.ux * p.ux +
                                                p.uy * p.uy + p.uz * p.uz);
                d.px_created += Const::me * Const::c * p.weight * p.ux;
                d.energy_created += Const::me * Const::c * Const::c *
                                     p.weight * (gamma - 1.0);
                d.jx_dx_created += -Const::qe * Const::c * p.weight * p.ux / gamma;
                d.pixx_dx_created += Const::me * Const::c * Const::c *
                                     p.weight * p.ux * p.ux / gamma;
                d.piperp_dx_created += Const::me * Const::c * Const::c *
                    p.weight * (p.uy * p.uy + p.uz * p.uz) / gamma;
            }
        }
    }
    const double scales[6] = { d.number_scale, d.px_scale, d.jx_scale,
                               d.energy_scale, d.pixx_scale, d.piperp_scale };
    const double deltas[6] = {
        d.number_created - d.number_removed,
        d.px_created - d.px_removed,
        d.jx_dx_created - d.jx_dx_removed,
        d.energy_created - d.energy_removed,
        d.pixx_dx_created - d.pixx_dx_removed,
        d.piperp_dx_created - d.piperp_dx_removed
    };
    double* residuals[6] = { &d.number_residual_rel, &d.px_residual_rel,
                             &d.jx_residual_rel, &d.energy_residual_rel,
                             &d.pixx_residual_rel, &d.piperp_residual_rel };
    for (int m = 0; m < 6; ++m) {
        *residuals[m] = std::fabs(deltas[m]) / std::max(scales[m], 1.0e-300);
    }
    d.finite = d.finite && candidate.finite() && candidate.nonnegative_weights();
    d.conservative = d.finite && d.number_residual_rel <= k_hard_gate &&
                     d.px_residual_rel <= k_hard_gate &&
                     d.energy_residual_rel <= k_hard_gate;
    d.fidelity_ok = d.finite && d.jx_residual_rel <= k_fidelity_gate &&
                    d.pixx_residual_rel <= k_fidelity_gate &&
                    d.piperp_residual_rel <= k_fidelity_gate;
    d.complete = d.finite && d.conservative && d.fidelity_ok &&
                 d.support_limit_violation_count == 0 &&
                 d.duplicate_id_count == 0 &&
                 d.face_ledger_mismatch_count == 0;
    if (d.complete) tail_trial.swap_state(candidate);
    return d;
}

BulkTailConversionDiagnostics BulkTailConverter::convert_initial_tail_cells(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid, const HybridVelocityPartition& partition,
    int accepted_step, int mpi_rank, size_t max_supports)
{
    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    std::vector<ConversionMassRequest> requests;
    const int ng = grid.nghost;
    const int nxl = grid.nx_local;
    const int nmu = partition.uperp_count;
    bool input_valid = true;
    for (int il = 0; il < nxl && input_valid; ++il) {
        const int ix_global = grid.global_cell(ng + il);
        for (int iv = 0; iv < partition.upar_count && input_valid; ++iv) {
            for (int imu = 0; imu < nmu; ++imu) {
                if (!partition.is_tail_owned(iv, imu)) continue;
                const size_t slot = idx3(ng + il, iv, imu);
                const double mass = bulk_trial.f[slot];
                if (!std::isfinite(mass)) {
                    input_valid = false;
                    break;
                }
                const double floor = 256.0 *
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::fabs(mass));
                if (mass < -floor) {
                    input_valid = false;
                    break;
                }
                if (!(mass > floor)) continue;

                BulkTailFluxParcel parcel;
                parcel.ix_local = il;
                parcel.ix_global = ix_global;
                parcel.direction = VelocityFaceDirection::U_PARALLEL;
                // Initialization is not an interface crossing event.  The
                // sentinel face is retained only for provenance; the loader
                // never uses it to construct a physical current.
                parcel.face_index = -1;
                parcel.transverse_index = imu;
                parcel.operator_stage = 0;
                parcel.face_number = mass;
                const std::vector<TailSubcellNode> nodes =
                    TailSubcellQuadrature::nodes(
                        bulk_trial.cgrid, iv, imu, 4);
                if (nodes.empty()) {
                    input_valid = false;
                    break;
                }
                for (size_t q = 0; q < nodes.size(); ++q) {
                    const double node_mass = mass * nodes[q].mass_fraction;
                    if (!bulk_tail_parcel_add_node(
                            parcel, nodes[q].upar, nodes[q].uperp,
                            node_mass)) {
                        input_valid = false;
                        break;
                    }
                }
                if (!input_valid) break;
                batch.parcels.push_back(parcel);
                ConversionMassRequest request;
                request.ix_global = ix_global;
                request.iv = iv;
                request.imu = imu;
                request.mass = mass;
                requests.push_back(request);
            }
        }
    }

    BulkTailConversionDiagnostics result;
    result.accepted_step = accepted_step;
    result.location = ConversionLocation::AFTER_U_SUBSTEP;
    if (!input_valid) return result;
    batch.recompute(partition.min_conversion_energy);
    if (!batch.finite || !batch.nonnegative || batch.duplicate_count != 0)
        return result;

    BackgroundTailPIC tail_candidate = tail_trial;
    result = convert_flux_batch(batch, tail_candidate, grid, partition,
                                accepted_step,
                                ConversionLocation::AFTER_U_SUBSTEP,
                                mpi_rank, max_supports);
    if (!result.complete) return result;

    Species bulk_candidate = bulk_trial;
    // Tail-owned cells must be empty after the one-time initialization
    // transfer.  Treat only the already-admitted roundoff debt as zero; any
    // material mass has an explicit parcel/request above.
    for (int il = 0; il < nxl; ++il) {
        for (int iv = 0; iv < partition.upar_count; ++iv) {
            for (int imu = 0; imu < nmu; ++imu) {
                if (!partition.is_tail_owned(iv, imu)) continue;
                const size_t slot = idx3(ng + il, iv, imu);
                const double value = bulk_candidate.f[slot];
                const double floor = 256.0 *
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::fabs(value));
                if (value < 0.0 && value >= -floor)
                    bulk_candidate.f[slot] = 0.0;
            }
        }
    }
    const SpeciesConversionResult removed =
        bulk_candidate.extract_conversion_masses(requests);
    if (!removed.valid || !removed.complete) {
        result.complete = false;
        result.conservative = false;
        result.fidelity_ok = false;
        result.finite = false;
        return result;
    }
    const double scale = std::max(1.0, std::fabs(result.number_removed));
    if (std::fabs(removed.number_removed - result.number_removed) >
        1024.0 * std::numeric_limits<double>::epsilon() * scale) {
        result.complete = false;
        result.conservative = false;
        result.finite = false;
        return result;
    }
    bulk_trial.swap_state(bulk_candidate);
    tail_trial.swap_state(tail_candidate);
    return result;
}
