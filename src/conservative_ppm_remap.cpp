#include "conservative_ppm_remap.h"
#include "field_particle_power_audit.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <stdexcept>

namespace {
inline size_t nvmu() { return static_cast<size_t>(Param::Nvmu); }

// Keep the flux-interface representation change on the same numerical floor
// as the static cell extractor.  A source below this scale cannot be
// distinguished from a finite-volume roundoff residue of n0*dx, so it must
// remain in the Eulerian representation rather than creating PIC parcels.
inline double interface_conversion_roundoff_floor(const SpatialGrid& grid)
{
    return 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, Param::dens * grid.dx);
}
}

ConservativePpmRemap::ConservativePpmRemap()
    : cgrid_(NULL), upar_nv_(0), upar_nmu_(0), max_abs_vx_(0.0),
      halo_width_(0), nx_local_(0), ng_(0),
      x_transport_velocity_mode_(XTransportVelocityMode::ANALYTIC_CELL_CENTER)
{}

double ConservativePpmRemap::transport_velocity(size_t q) const
{
    if (x_transport_velocity_mode_ ==
        XTransportVelocityMode::ENERGY_CONJUGATE_CELL) {
        if (q >= cgrid_->vx_energy_conjugate_cell.size())
            throw std::runtime_error("energy-conjugate transport velocity index out of range");
        const double value = cgrid_->vx_energy_conjugate_cell[q];
        if (!std::isfinite(value))
            throw std::runtime_error("nonfinite energy-conjugate transport velocity");
        return value;
    }
    return cgrid_->vx[q];
}

void ConservativePpmRemap::set_x_transport_velocity_mode(
    XTransportVelocityMode mode)
{
    x_transport_velocity_mode_ = mode;
    max_abs_vx_ = 0.0;
    const std::vector<double>& table =
        mode == XTransportVelocityMode::ENERGY_CONJUGATE_CELL
            ? cgrid_->vx_energy_conjugate_cell : cgrid_->vx;
    for (size_t q = 0; q < table.size(); ++q)
        max_abs_vx_ = std::max(max_abs_vx_, std::fabs(table[q]));
}

void ConservativePpmRemap::init(const SpatialGrid& grid,
                                const CylindricalVelocityGrid& velocity_grid)
{
    grid_ = grid;
    cgrid_ = &velocity_grid;
    nx_local_ = grid_.nx_local;
    ng_ = grid_.nghost;
    upar_nv_ = static_cast<int>(velocity_grid.upar_cells.size());
    upar_nmu_ = static_cast<int>(velocity_grid.uperp_cells.size());
    max_abs_vx_ = 0.0;
    for (int iv = 0; iv < upar_nv_; ++iv) {
        for (int imu = 0; imu < upar_nmu_; ++imu) {
            max_abs_vx_ = std::max(
                max_abs_vx_,
                std::fabs(velocity_grid.vx[static_cast<size_t>(iv) *
                                               upar_nmu_ + imu]));
        }
    }
    compute_upar_face_coefficients(velocity_grid);
}

void ConservativePpmRemap::compute_upar_face_coefficients(
    const CylindricalVelocityGrid& grid)
{
    const int nv = upar_nv_;
    upar_average_.assign(static_cast<size_t>(nv), 0.0);
    upar_left_edge_.assign(static_cast<size_t>(nv) + 1, 0.0);
    upar_right_edge_.assign(static_cast<size_t>(nv), 0.0);
    upar_curvature_.assign(static_cast<size_t>(nv), 0.0);
    upar_cumulative_.assign(static_cast<size_t>(nv) + 1, 0.0);
    upar_swept_.assign(static_cast<size_t>(nv) + 1, 0.0);
    upar_face_coeffs_.assign(static_cast<size_t>(nv + 1) * 4, 0.0);
    upar_stencil_base_.assign(nv + 1, 0);
    if (nv < 4) return;

    // Interior faces use the quartic-through-four-cell-averages interface
    // value: on a uniform grid this reduces exactly to the PPM formula
    // (7/12)(a_i + a_{i+1}) - (1/12)(a_{i-1} + a_{i+2}).  Domain-edge faces
    // are pinned to the adjacent cell average (same convention as the x
    // path), and the two faces next to the edge use a one-sided stencil.
    for (int f = 1; f < nv; ++f) {
        int s = f - 2;
        if (s < 0) s = 0;
        if (s > nv - 4) s = nv - 4;
        upar_stencil_base_[static_cast<size_t>(f)] = s;

        const double xf = grid.upar_faces[static_cast<size_t>(f)];
        double mat[4][4];
        double inv[4][4];
        for (int m = 0; m < 4; ++m) {
            const double l = grid.upar_faces[static_cast<size_t>(s + m)] - xf;
            const double r =
                grid.upar_faces[static_cast<size_t>(s + m + 1)] - xf;
            mat[m][0] = r - l;
            mat[m][1] = 0.5 * (r * r - l * l);
            mat[m][2] = (r * r * r - l * l * l) / 3.0;
            mat[m][3] = (r * r * r * r - l * l * l * l) / 4.0;
            for (int c = 0; c < 4; ++c) inv[m][c] = (m == c) ? 1.0 : 0.0;
        }
        for (int col = 0; col < 4; ++col) {
            int pivot = col;
            for (int row = col + 1; row < 4; ++row) {
                if (std::fabs(mat[row][col]) > std::fabs(mat[pivot][col])) {
                    pivot = row;
                }
            }
            if (std::fabs(mat[pivot][col]) <= 0.0) {
                throw std::runtime_error(
                    "singular non-uniform u_parallel face stencil");
            }
            if (pivot != col) {
                for (int c = 0; c < 4; ++c) {
                    std::swap(mat[pivot][c], mat[col][c]);
                    std::swap(inv[pivot][c], inv[col][c]);
                }
            }
            const double d = mat[col][col];
            for (int c = 0; c < 4; ++c) {
                mat[col][c] /= d;
                inv[col][c] /= d;
            }
            for (int row = 0; row < 4; ++row) {
                if (row == col) continue;
                const double factor = mat[row][col];
                if (factor == 0.0) continue;
                for (int c = 0; c < 4; ++c) {
                    mat[row][c] -= factor * mat[col][c];
                    inv[row][c] -= factor * inv[col][c];
                }
            }
        }
        // The interface value v = c[0] solves M c = b with
        // b[m] = w_{s+m} * a_{s+m}, so the coefficient on each cell average
        // is inv[0][m] times that cell's width.
        for (int m = 0; m < 4; ++m) {
            const double w = grid.upar_faces[static_cast<size_t>(s + m + 1)] -
                             grid.upar_faces[static_cast<size_t>(s + m)];
            upar_face_coeffs_[static_cast<size_t>(f) * 4 + m] =
                inv[0][m] * w;
        }
    }
    // Face 0: value = average[0]; face nv: value = average[nv-1].
    upar_face_coeffs_[0] = 1.0;
    upar_stencil_base_[0] = 0;
    upar_face_coeffs_[static_cast<size_t>(nv) * 4 + 3] = 1.0;
    upar_stencil_base_[static_cast<size_t>(nv)] = nv - 4;
}

size_t ConservativePpmRemap::upar_index(int ix, int j_upar,
                                        int k_uperp) const
{
    return (static_cast<size_t>(ix) * static_cast<size_t>(upar_nv_) +
            static_cast<size_t>(j_upar)) *
               static_cast<size_t>(upar_nmu_) + k_uperp;
}

int ConservativePpmRemap::locate_upar(double u) const
{
    const int nv = upar_nv_;
    if (u <= cgrid_->upar_faces[0]) return 0;
    if (u >= cgrid_->upar_faces[static_cast<size_t>(nv)]) return nv - 1;
    int lo = 0;
    int hi = nv;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (cgrid_->upar_faces[static_cast<size_t>(mid)] <= u) lo = mid;
        else hi = mid;
    }
    return lo;
}

double ConservativePpmRemap::upar_parabola_integral(int cell, double u1,
                                                    double u2) const
{
    if (!(u2 > u1)) return 0.0;
    const double w = cgrid_->upar_widths[static_cast<size_t>(cell)];
    const double cell_lo = cgrid_->upar_faces[static_cast<size_t>(cell)];
    const double L = upar_left_edge_[static_cast<size_t>(cell)];
    const double R = upar_right_edge_[static_cast<size_t>(cell)];
    const double cur = upar_curvature_[static_cast<size_t>(cell)];
    const double s1 = (u1 - cell_lo) / w;
    const double s2 = (u2 - cell_lo) / w;
    const double ds = s2 - s1;
    const double ds2 = s2 * s2 - s1 * s1;
    const double ds3 = s2 * s2 * s2 - s1 * s1 * s1;
    return w * (L * ds + (R - L) * 0.5 * ds2 +
                cur * (0.5 * ds2 - ds3 / 3.0));
}

double ConservativePpmRemap::upar_swept_mass(int face, double a_u, double dt,
                                             int& cells_spanned) const
{
    if (a_u == 0.0 || dt <= 0.0) {
        cells_spanned = 0;
        return 0.0;
    }
    const int nv = upar_nv_;
    const double u_f = cgrid_->upar_faces[static_cast<size_t>(face)];
    const double shift = a_u * dt;
    double lo = std::min(u_f, u_f - shift);
    double hi = std::max(u_f, u_f - shift);
    lo = std::max(lo, cgrid_->upar_faces[0]);
    hi = std::min(hi, cgrid_->upar_faces[static_cast<size_t>(nv)]);
    if (!(hi > lo)) {
        cells_spanned = 0;
        return 0.0;
    }
    const int c_lo = locate_upar(lo);
    const int c_hi = locate_upar(hi);
    cells_spanned = c_hi - c_lo + 1;
    double mass = 0.0;
    if (c_hi == c_lo) {
        mass = upar_parabola_integral(c_lo, lo, hi);
    } else {
        mass = upar_parabola_integral(
                   c_lo, lo, cgrid_->upar_faces[static_cast<size_t>(c_lo) + 1])
             + (upar_cumulative_[static_cast<size_t>(c_hi)] -
                upar_cumulative_[static_cast<size_t>(c_lo) + 1])
             + upar_parabola_integral(
                   c_hi, cgrid_->upar_faces[static_cast<size_t>(c_hi)], hi);
    }
    return (a_u > 0.0) ? mass : -mass;
}

double ConservativePpmRemap::upar_tail_energy(int face, double a_u, double dt,
                                              int k_uperp,
                                              double particle_mass) const
{
    if (a_u == 0.0 || dt <= 0.0) return 0.0;
    const int nv = upar_nv_;
    const double u_f = cgrid_->upar_faces[static_cast<size_t>(face)];
    const double shift = a_u * dt;
    double lo = std::min(u_f, u_f - shift);
    double hi = std::max(u_f, u_f - shift);
    lo = std::max(lo, cgrid_->upar_faces[0]);
    hi = std::min(hi, cgrid_->upar_faces[static_cast<size_t>(nv)]);
    if (!(hi > lo)) return 0.0;
    const int c_lo = locate_upar(lo);
    const int c_hi = locate_upar(hi);
    const double uperp = cgrid_->uperp_cells[static_cast<size_t>(k_uperp)];
    // 8-point Gauss-Legendre on [-1,1].
    static const double gl_x[4] = {
        0.1834346424956498, 0.5255324099163290,
        0.7966664774136267, 0.9602898564975363
    };
    static const double gl_w[4] = {
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763
    };
    double energy = 0.0;
    for (int c = c_lo; c <= c_hi; ++c) {
        const double piece_lo =
            std::max(lo, cgrid_->upar_faces[static_cast<size_t>(c)]);
        const double piece_hi =
            std::min(hi, cgrid_->upar_faces[static_cast<size_t>(c) + 1]);
        if (!(piece_hi > piece_lo)) continue;
        const double width = piece_hi - piece_lo;
        const double mid = 0.5 * (piece_lo + piece_hi);
        const double cell_lo = cgrid_->upar_faces[static_cast<size_t>(c)];
        const double inv_cell_w =
            1.0 / cgrid_->upar_widths[static_cast<size_t>(c)];
        const double L = upar_left_edge_[static_cast<size_t>(c)];
        const double R = upar_right_edge_[static_cast<size_t>(c)];
        const double cur = upar_curvature_[static_cast<size_t>(c)];
        for (int g = 0; g < 4; ++g) {
            for (int sign = -1; sign <= 1; sign += 2) {
                const double u = mid + sign * 0.5 * width * gl_x[g];
                const double s = (u - cell_lo) * inv_cell_w;
                const double mu =
                    L + (R - L) * s + cur * s * (1.0 - s);
                const double gamma =
                    std::sqrt(1.0 + u * u + uperp * uperp);
                const double ke =
                    particle_mass * Const::c * Const::c * (gamma - 1.0);
                energy += 0.5 * width * gl_w[g] * mu * ke;
            }
        }
    }
    return energy;
}

double ConservativePpmRemap::upar_profile_value(int cell, double u) const
{
    const double lo = cgrid_->upar_faces[static_cast<size_t>(cell)];
    const double width = cgrid_->upar_widths[static_cast<size_t>(cell)];
    if (!(width > 0.0)) return 0.0;
    const double s = std::max(0.0, std::min(1.0, (u - lo) / width));
    return upar_left_edge_[static_cast<size_t>(cell)] +
           (upar_right_edge_[static_cast<size_t>(cell)] -
            upar_left_edge_[static_cast<size_t>(cell)]) * s +
           upar_curvature_[static_cast<size_t>(cell)] * s * (1.0 - s);
}

bool ConservativePpmRemap::append_upar_swept_nodes(
    int face, double a_u, double dt, int k_uperp, int quadrature_order,
    BulkTailFluxParcel& parcel, ParcelNodeFailure* failure) const
{
    if (failure != NULL) *failure = ParcelNodeFailure();
    const auto fail = [failure](ParcelNodeFailureReason reason,
                                double node_mass, double target,
                                double node_sum, double scale) {
        if (failure != NULL) {
            failure->reason = reason;
            failure->node_mass = node_mass;
            failure->target = target;
            failure->node_sum = node_sum;
            failure->scale = scale;
        }
        return false;
    };
    if (a_u == 0.0 || dt <= 0.0) return true;
    const double uf = cgrid_->upar_faces[static_cast<size_t>(face)];
    const double shift = a_u * dt;
    const double lo = std::max(std::min(uf, uf - shift),
                               cgrid_->upar_faces.front());
    const double hi = std::min(std::max(uf, uf - shift),
                               cgrid_->upar_faces.back());
    if (!(hi > lo)) return true;
    static const double x4[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526
    };
    static const double w4[4] = {
        0.3478548451374539, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374539
    };
    static const double x8[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363
    };
    static const double w8[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763
    };
    const int nq = quadrature_order >= 8 ? 8 : 4;
    const double* gx = nq == 8 ? x8 : x4;
    const double* gw = nq == 8 ? w8 : w4;
    const double ut0 = cgrid_->uperp_faces[static_cast<size_t>(k_uperp)];
    const double ut1 = cgrid_->uperp_faces[static_cast<size_t>(k_uperp) + 1];
    const double ut_measure = ut1 * ut1 - ut0 * ut0;
    if (!(ut_measure > 0.0)) {
        return fail(ParcelNodeFailureReason::InvalidTransverseMeasure,
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN());
    }
    const int c0 = locate_upar(lo);
    const int c1 = locate_upar(hi);
    double node_sum = 0.0;
    for (int cell = c0; cell <= c1; ++cell) {
        const double a = std::max(lo, cgrid_->upar_faces[static_cast<size_t>(cell)]);
        const double b = std::min(hi, cgrid_->upar_faces[static_cast<size_t>(cell) + 1]);
        if (!(b > a)) continue;
        const double um = 0.5 * (a + b);
        const double uh = 0.5 * (b - a);
        const double tm = 0.5 * (ut0 + ut1);
        const double th = 0.5 * (ut1 - ut0);
        for (int gu = 0; gu < nq; ++gu) {
            const double up = um + uh * gx[gu];
            const double fbar = upar_profile_value(cell, up);
            const double cell_average =
                upar_average_[static_cast<size_t>(cell)];
            for (int gt = 0; gt < nq; ++gt) {
                const double ut = tm + th * gx[gt];
                // Integrate in s=u_perp^2.  After mapping the Gauss node
                // from [-1,1] to the physical u_perp interval, the
                // cylindrical measure contributes
                // (u_perp * (u1-u0) / (u1^2-u0^2)) dxi.  The width factor
                // is essential for the parcel mass contract even though
                // the later face-mass normalization cancels it for a single
                // transverse cell.
                const double wt = gw[gt] * (ut1 - ut0) * ut / ut_measure;
                const double mass = fbar * uh * gw[gu] * wt;
                // This tolerance must be relative to the reconstructed local
                // profile.  An absolute O(epsilon) floor treats every
                // negative node in a dilute tail as roundoff, even when it is
                // order-one relative to the local swept flux.
                const double profile_scale = std::max(
                    std::fabs(upar_average_[static_cast<size_t>(cell)]),
                    std::max(std::fabs(upar_left_edge_[static_cast<size_t>(cell)]),
                             std::max(std::fabs(upar_right_edge_[static_cast<size_t>(cell)]),
                                      std::fabs(upar_curvature_[static_cast<size_t>(cell)]))));
                const double node_scale =
                    std::fabs(uh * gw[gu] * wt) * profile_scale;
                const double tolerance = 256.0 *
                    std::numeric_limits<double>::epsilon() * node_scale;
                if (!std::isfinite(mass)) {
                    return fail(ParcelNodeFailureReason::NonfiniteQuadratureNode,
                                mass, std::numeric_limits<double>::quiet_NaN(),
                                node_sum, std::numeric_limits<double>::quiet_NaN());
                }
                if (mass < -tolerance) {
                    // The PIC representation has a non-negative-mass
                    // contract.  A PPM quadrature node below zero is not a
                    // particle support point, irrespective of whether it was
                    // caused by an inherited negative average or an interior
                    // reconstruction overshoot.  Drop only that support;
                    // the complete retained node set is normalized to the
                    // unchanged final FV face mass below.  This projection is
                    // confined to bulk-to-tail conversion and never changes
                    // the Eulerian Vlasov state or its conservative face flux.
                    (void)cell_average;
                    continue;
                }
                const double stored = mass < 0.0 ? 0.0 : mass;
                // The exported parcel is the material after the force
                // characteristic has crossed the shared face.  Keep the
                // swept mass from the final FV flux, but carry the node to
                // its arrival momentum so the converter audits the tail-side
                // phase-space support rather than the bulk departure cell.
                const double arrival_up = std::max(
                    cgrid_->upar_faces.front(), std::min(
                        cgrid_->upar_faces.back(), up + shift));
                parcel.nodes.push_back(
                    FluxParcelNode{arrival_up, ut, stored});
                node_sum += stored;
            }
        }
    }
    int spanned_dummy = 0;
    const double target = std::fabs(
        upar_swept_mass(face, a_u, dt, spanned_dummy));
    const double zero_support_limit = 64.0 *
        std::numeric_limits<double>::denorm_min();
    if (target <= zero_support_limit) {
        // No double-precision parcel mass can be represented at this face.
        // Treat it as zero support rather than reporting a failed positive
        // normalization with node_sum==0.  The caller records this as a
        // discarded roundoff face and does not alter the bulk state.
        parcel.nodes.clear();
        parcel.recompute_moments();
        if (failure != NULL) {
            failure->target = 0.0;
            failure->node_sum = 0.0;
        }
        return true;
    }
    if (!parcel.nodes.empty() && target > 0.0) {
        // The PPM swept integral is the conservation authority.  Quadrature
        // of the same reconstructed profile can differ by a few ulps (and,
        // near a limiter kink, by more than one ulp).  Rescale the complete
        // positive node set instead of correcting one node, which could
        // create a negative support and invalidate the whole transaction.
        if (!(node_sum > 0.0) || !std::isfinite(node_sum)) {
            return fail(ParcelNodeFailureReason::NonpositiveNodeSum,
                        std::numeric_limits<double>::quiet_NaN(), target,
                        node_sum, std::numeric_limits<double>::quiet_NaN());
        }
        const double scale = target / node_sum;
        if (!(scale >= 0.0) || !std::isfinite(scale)) {
            return fail(ParcelNodeFailureReason::NonfiniteScale,
                        std::numeric_limits<double>::quiet_NaN(), target,
                        node_sum, scale);
        }
        for (size_t q = 0; q < parcel.nodes.size(); ++q)
            parcel.nodes[q].mass *= scale;
    }
    parcel.recompute_moments();
    if (!parcel.finite_nonnegative()) {
        return fail(ParcelNodeFailureReason::InvalidParcelMoments,
                    std::numeric_limits<double>::quiet_NaN(), target,
                    node_sum, std::numeric_limits<double>::quiet_NaN());
    }
    if (failure != NULL) {
        // Successful construction still returns the independently
        // reconstructed target so the caller can distinguish an actual
        // positive swept transfer from a cancellation residue in the stored
        // face vector.
        failure->target = target;
        failure->node_sum = node_sum;
    }
    return true;
}

void ConservativePpmRemap::build_upar_reconstruction(
    long long& constant_cells, long long& linear_cells, bool count)
{
    const int nv = upar_nv_;
    for (int f = 0; f <= nv; ++f) {
        const int s = upar_stencil_base_[static_cast<size_t>(f)];
        const double* coeff =
            &upar_face_coeffs_[static_cast<size_t>(f) * 4];
        double v = 0.0;
        for (int m = 0; m < 4; ++m) {
            v += coeff[m] * upar_average_[static_cast<size_t>(s + m)];
        }
        upar_left_edge_[static_cast<size_t>(f)] = v;
    }

    const double eps = std::numeric_limits<double>::epsilon();
    for (int j = 0; j < nv; ++j) {
        const double a = upar_average_[static_cast<size_t>(j)];
        double L = upar_left_edge_[static_cast<size_t>(j)];
        double R = upar_left_edge_[static_cast<size_t>(j + 1)];

        // Colella-Woodward monotonicity constraint (no interior extrema).
        if ((R - a) * (a - L) <= 0.0) {
            L = a;
            R = a;
        } else {
            const double d = R - L;
            const double m = a - 0.5 * (L + R);
            if (d * m > d * d / 6.0) {
                R = 3.0 * a - 2.0 * L;
            } else if (-d * m > d * d / 6.0) {
                L = 3.0 * a - 2.0 * R;
            }
        }

        if (a < 0.0) {
            // Roundoff-scale negative input: record the debt and reduce the
            // reconstruction to constant in this cell only.
            L = a;
            R = a;
        } else {
            double cur = 6.0 * a - 3.0 * (L + R);
            double parabola_min = std::min(L, R);
            if (cur != 0.0) {
                const double s_star = (R - L + cur) / (2.0 * cur);
                if (s_star > 0.0 && s_star < 1.0) {
                    parabola_min = std::min(
                        parabola_min,
                        L + (R - L) * s_star + cur * s_star * (1.0 - s_star));
                }
            }
            // Positivity is a relative property of this velocity cell.  Do
            // not use an absolute O(epsilon) floor here: the physical tail
            // can be far below one in the stored distribution units.
            const double profile_scale = std::max(
                std::fabs(a), std::max(std::fabs(L), std::fabs(R)));
            const double pos_tol = 128.0 * eps * profile_scale;
            if (parabola_min < -pos_tol) {
                // Blend the entire PPM polynomial back to its cell average.
                // For p(s) with mean a, p_theta(s)=a+theta*(p(s)-a) has the
                // same finite-volume average for every theta.  Choosing theta
                // from the true polynomial minimum is therefore the least
                // dissipative scalar correction that makes every quadrature
                // node non-negative.  Limiting L/R independently does not
                // provide that guarantee because the curvature may still
                // create a negative interior extremum.
                const double denominator = a - parabola_min;
                const double theta = denominator > 0.0
                    ? std::max(0.0, std::min(1.0, a / denominator))
                    : 0.0;
                L = a + theta * (L - a);
                R = a + theta * (R - a);
            }
        }

        upar_left_edge_[static_cast<size_t>(j)] = L;
        upar_right_edge_[static_cast<size_t>(j)] = R;
        upar_curvature_[static_cast<size_t>(j)] =
            6.0 * a - 3.0 * (L + R);

        if (count) {
            const double cur = upar_curvature_[static_cast<size_t>(j)];
            const bool is_constant = (L == a && R == a);
            const bool is_linear =
                !is_constant &&
                std::fabs(cur) <= 128.0 * eps * std::max(1.0, std::fabs(a));
            if (is_constant) ++constant_cells;
            else if (is_linear) ++linear_cells;
        }
    }

    // Prefix masses (particles) for the full-cell part of departure
    // intervals on the non-uniform grid.
    upar_cumulative_[0] = 0.0;
    for (int j = 0; j < nv; ++j) {
        upar_cumulative_[static_cast<size_t>(j) + 1] =
            upar_cumulative_[static_cast<size_t>(j)] +
            upar_average_[static_cast<size_t>(j)] *
                cgrid_->upar_widths[static_cast<size_t>(j)];
    }
}

RemapDiagnostics ConservativePpmRemap::advect_u_parallel(
    const Species& input, Species& output, const EMFields& fields,
    double dt, double time, const HybridVelocityPartition* partition,
    BulkTailFluxBatch* exported_flux, int quadrature_order,
    std::vector<double>* local_delta_ke_by_x)
{
    RemapDiagnostics diag = {};
    diag.finite = true;
    diag.audit_valid = true;
    (void)time;
    if (local_delta_ke_by_x != NULL) {
        local_delta_ke_by_x->assign(static_cast<size_t>(nx_local_), 0.0);
    }
    const bool interface_mode = partition != NULL && exported_flux != NULL;
    const bool apply_interface_sink = interface_mode &&
        exported_flux->apply_interface_sink;
    if (interface_mode) exported_flux->clear();
    const int nxl = nx_local_;
    const int ng = ng_;
    const int nxt = grid_.nx_total;
    const int nv = upar_nv_;
    const int nmu = upar_nmu_;

    if (!(dt > 0.0) || nxl <= 0) {
        output.f = input.f;
        return diag;
    }

    double local_before = 0.0;
    double local_after = 0.0;
    double local_tail_number = 0.0;
    double local_tail_energy = 0.0;
    double local_min_mass = std::numeric_limits<double>::infinity();
    double local_input_min = std::numeric_limits<double>::infinity();
    double local_input_max = -std::numeric_limits<double>::infinity();
    double local_output_min = std::numeric_limits<double>::infinity();
    double local_output_max = -std::numeric_limits<double>::infinity();
    double local_max_departure = 0.0;
    long long constant_cells = 0;
    long long linear_cells = 0;
    bool local_finite = true;
    int local_bad_ix = -1;
    int local_bad_j = -1;
    int local_bad_k = -1;
    double local_bad_value = std::numeric_limits<double>::quiet_NaN();
    bool local_audit_valid = true;
    bool local_audit_state_equal = true;
    int local_audit_failure_code = 0;
    bool local_parcel_failure_seen = false;
    ParcelNodeFailure local_parcel_failure;
    int local_parcel_failure_ix = -1;
    int local_parcel_failure_face = -1;
    int local_parcel_failure_iuperp = -1;
    int mpi_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    double local_interface_face_number = 0.0;
    double local_interface_number = 0.0;
    double local_interface_energy = 0.0;
    double local_interface_roundoff_discarded = 0.0;
    double local_tail_owned_expected_transfer = 0.0;
    double local_tail_owned_roundoff_discarded = 0.0;
    double local_tail_owned_residual = 0.0;
    bool local_tail_component_failure_reported = false;
    std::uint64_t local_interface_parcels = 0;
    std::uint64_t local_interface_nodes = 0;
    // Gate C (section 7.3): bulk u_parallel discrete kinetic-energy identity.
    double local_delta_bulk_ke = 0.0;
    double local_internal_face_energy = 0.0;
    double local_left_face_energy = 0.0;
    double local_right_face_energy = 0.0;
    double local_interface_removed_energy = 0.0;

    // There may be many u_parallel faces, but for a fixed u_perp slice only
    // the first outward bulk-to-tail interface can contribute.  Pre-index it
    // once for each acceleration sign.  The previous audit path scanned all
    // (nv + 1) faces for every x/k slice just to discover this same topology.
    std::vector<int> first_positive_interface(static_cast<size_t>(nmu), -1);
    std::vector<int> first_negative_interface(static_cast<size_t>(nmu), -1);
    std::vector<std::vector<int> > interfaces_by_uperp(
        static_cast<size_t>(nmu));
    if (interface_mode) {
        for (size_t n = 0; n < partition->upar_interface_faces.size(); ++n) {
            const BulkTailInterfaceFace& iface =
                partition->upar_interface_faces[n];
            const int k = iface.transverse_index;
            if (k < 0 || k >= nmu) continue;
            interfaces_by_uperp[static_cast<size_t>(k)].push_back(
                static_cast<int>(n));
            if (iface.outward_sign > 0 &&
                (first_positive_interface[static_cast<size_t>(k)] < 0 ||
                 iface.face_index <
                    first_positive_interface[static_cast<size_t>(k)])) {
                first_positive_interface[static_cast<size_t>(k)] =
                    iface.face_index;
            }
            if (iface.outward_sign < 0 &&
                (first_negative_interface[static_cast<size_t>(k)] < 0 ||
                 iface.face_index >
                    first_negative_interface[static_cast<size_t>(k)])) {
                first_negative_interface[static_cast<size_t>(k)] =
                    iface.face_index;
            }
        }
    }

    for (int ix = 0; ix < nxt; ++ix) {
        const bool interior = ix >= ng && ix < ng + nxl;
        // JC stage 2/P2: use the same Yee-face field object that defines the
        // converged pairing map.  The old path read fields.Ex, which is a
        // cached cell array and can be a different time/source object when a
        // trial field is reconstructed.  Physical cells use the project's
        // canonical face-to-cell average; ghost cells retain their existing
        // halo value because they have no local physical face pair.
        double e_cell = fields.Ex[static_cast<size_t>(ix)];
        if (interior && fields.Ex_face.size() ==
                static_cast<size_t>(nxl + 1)) {
            const int face = ix - ng;
            e_cell = 0.5 *
                (fields.Ex_face[static_cast<size_t>(face)] +
                 fields.Ex_face[static_cast<size_t>(face + 1)]);
        }
        const double a_u = input.charge * e_cell /
                           (input.mass * Const::c);
        local_finite = local_finite && std::isfinite(a_u);
        if (!std::isfinite(a_u) && local_bad_ix < 0) {
            local_bad_ix = grid_.global_cell(ix);
            local_bad_j = -1;
            local_bad_k = -1;
            local_bad_value = a_u;
        }
        for (int k = 0; k < nmu; ++k) {
            // These identify the single one-way bulk-to-tail face that is
            // actually converted for this (x, u_perp) slice.  The matching
            // tail-owned cell receives exactly this transfer in the finite
            // volume update and must not be counted as unexplained residue.
            int selected_tail_cell = -1;
            double selected_tail_transfer = 0.0;
            for (int j = 0; j < nv; ++j) {
                const double m = input.f[upar_index(ix, j, k)];
                upar_average_[static_cast<size_t>(j)] =
                    m / cgrid_->upar_widths[static_cast<size_t>(j)];
                if (interior) {
                    local_before += m;
                    local_min_mass = std::min(local_min_mass, m);
                    local_input_min = std::min(local_input_min, m);
                    local_input_max = std::max(local_input_max, m);
                    local_finite = local_finite && std::isfinite(m);
                    if (!std::isfinite(m) && local_bad_ix < 0) {
                        local_bad_ix = grid_.global_cell(ix);
                        local_bad_j = j;
                        local_bad_k = k;
                        local_bad_value = m;
                    }
                }
            }
            build_upar_reconstruction(constant_cells, linear_cells, interior);
            for (int f = 0; f <= nv; ++f) {
                int spanned = 0;
                upar_swept_[static_cast<size_t>(f)] =
                    upar_swept_mass(f, a_u, dt, spanned);
                if (interior) {
                    local_max_departure = std::max(
                        local_max_departure, static_cast<double>(spanned));
                }
            }
            if (interface_mode && interior) {
                const int selected_face = a_u >= 0.0
                    ? first_positive_interface[static_cast<size_t>(k)]
                    : first_negative_interface[static_cast<size_t>(k)];
                const std::vector<int>& k_interfaces =
                    interfaces_by_uperp[static_cast<size_t>(k)];
                for (size_t interface_slot = 0;
                     interface_slot < k_interfaces.size(); ++interface_slot) {
                    const int face_id =
                        k_interfaces[interface_slot];
                    const BulkTailInterfaceFace& iface =
                        partition->upar_interface_faces[
                            static_cast<size_t>(face_id)];
                    const int f = iface.face_index;
                    const double audit_face_before =
                        upar_swept_[static_cast<size_t>(f)];
                    const bool read_only_face = !apply_interface_sink;
                    const auto verify_read_only_face = [&]() {
                        if (read_only_face &&
                            std::memcmp(&audit_face_before,
                                        &upar_swept_[static_cast<size_t>(f)],
                                        sizeof(double)) != 0) {
                            local_audit_state_equal = false;
                        }
                    };
                    const bool outward = a_u * static_cast<double>(iface.outward_sign) > 0.0;
                    if (!outward || f != selected_face) {
                        // In flux-interface mode the tail representation is
                        // one-way: a tail-owned cell is never an Eulerian
                        // donor back into bulk.  Explicitly zero both a
                        // tail-to-bulk face and every later interface after
                        // the first bulk-to-tail crossing.  Leave the audit
                        // path untouched so it remains bitwise read-only.
                        if (apply_interface_sink)
                            upar_swept_[static_cast<size_t>(f)] = 0.0;
                        verify_read_only_face();
                        continue;
                    }
                    const double swept = upar_swept_[static_cast<size_t>(f)];
                    // A non-negative distribution transported by a_u has a
                    // swept face mass with the same sign as a_u.  Taking an
                    // unconditional absolute value here turns an inherited
                    // negative-f numerical debt into a fictitious outward
                    // PIC source.  Such a debt is neither exported nor
                    // removed from the Eulerian state.
                    const double exported = a_u >= 0.0 ? swept : -swept;
                    if (!(exported > 0.0)) {
                        verify_read_only_face();
                        continue;
                    }
                    if (exported <= interface_conversion_roundoff_floor(grid_)) {
                        // Do not turn an unresolvable PPM face residue into a
                        // large collection of tail PIC supports.  In the
                        // production sink path, cancel the selected face as
                        // well: otherwise bulk would lose a quantity for
                        // which no tail parcel is created.  The audit path is
                        // read-only and records the same discarded amount.
                        local_interface_roundoff_discarded += exported;
                        if (apply_interface_sink)
                            upar_swept_[static_cast<size_t>(f)] = 0.0;
                        verify_read_only_face();
                        continue;
                    }
                    BulkTailFluxParcel parcel;
                    parcel.ix_local = ix - ng;
                    parcel.ix_global = grid_.global_cell(ix);
                    parcel.direction = VelocityFaceDirection::U_PARALLEL;
                    parcel.face_index = f;
                    parcel.transverse_index = k;
                    parcel.operator_stage = 1;
                    ParcelNodeFailure parcel_failure;
                    if (!append_upar_swept_nodes(f, a_u, dt, k,
                                                  quadrature_order, parcel,
                                                  &parcel_failure)) {
                        local_audit_valid = false;
                        local_audit_failure_code |= 1;
                        if (!local_parcel_failure_seen) {
                            local_parcel_failure_seen = true;
                            local_parcel_failure = parcel_failure;
                            local_parcel_failure_ix = grid_.global_cell(ix);
                            local_parcel_failure_face = f;
                            local_parcel_failure_iuperp = k;
                        }
                        // A failed quadrature construction is not a valid
                        // tail conversion.  In particular, do not append an
                        // empty parcel or consume this one-way interface;
                        // doing so would make the audit failure alter later
                        // production conversion topology.
                        verify_read_only_face();
                        continue;
                    }
                    parcel.audit_node_failure_reason = static_cast<int>(
                        parcel_failure.reason);
                    parcel.audit_reconstructed_target = parcel_failure.target;
                    parcel.audit_node_sum = parcel_failure.node_sum;
                    const bool zero_support_roundoff =
                        parcel_failure.reason == ParcelNodeFailureReason::None &&
                        parcel.number == 0.0 && parcel_failure.target == 0.0;
                    if (zero_support_roundoff) {
                        // The stored face residual is positive only because
                        // of cancellation or subnormal underflow, while the
                        // same final PPM profile has no representable
                        // quadrature support.  Do not emit a zero-weight
                        // parcel or let this residue block a later physical
                        // interface crossing.
                        local_interface_roundoff_discarded += exported;
                        if (apply_interface_sink)
                            upar_swept_[static_cast<size_t>(f)] = 0.0;
                        verify_read_only_face();
                        continue;
                    }
                    // Keep the originating final PPM face transfer separate
                    // from the quadrature reconstruction.  The parcel
                    // number is recomputed from its nodes and is an audit
                    // quantity; face_number is the actual conservative
                    // transfer used by the bulk update.
                    parcel.face_number = exported;
                    exported_flux->parcels.push_back(parcel);
                    local_interface_face_number += exported;
                    ++local_interface_parcels;
                    local_interface_nodes += parcel.nodes.size();
                    local_interface_number += parcel.number;
                    local_interface_energy += parcel.kinetic_energy;
                    if (apply_interface_sink) {
                        selected_tail_cell = a_u >= 0.0 ? f : f - 1;
                        selected_tail_transfer = exported;
                    }
                    verify_read_only_face();
                }
            }
            // Velocity-domain outflow ledger: both ends default to zero
            // inflow; mass (and its kinetic energy) leaving through the outer
            // faces is recorded instead of being wrapped to the other side.
            if (interior) {
                if (upar_swept_[0] < 0.0) {
                    local_tail_number += -upar_swept_[0];
                    local_tail_energy +=
                        upar_tail_energy(0, a_u, dt, k, input.mass);
                }
                if (upar_swept_[static_cast<size_t>(nv)] > 0.0) {
                    local_tail_number +=
                        upar_swept_[static_cast<size_t>(nv)];
                    local_tail_energy +=
                        upar_tail_energy(nv, a_u, dt, k, input.mass);
                }
            }
            // Conservative flux-difference update.
            // Validate each connected tail-owned interval by its signed
            // finite-volume balance.  Internal PPM faces telescope inside
            // the interval and must not be counted with a cell-wise L1 norm:
            // those faces only redistribute mass between Eulerian cells that
            // are all cleared in the one-way PIC representation.
            bool tail_component_active = false;
            int tail_component_start = -1;
            int tail_component_end = -1;
            double tail_component_number = 0.0;
            double tail_component_expected = 0.0;
            const auto finish_tail_component = [&]() {
                if (!tail_component_active) return;
                local_tail_owned_expected_transfer += tail_component_expected;
                const double unexplained = std::fabs(
                    tail_component_number - tail_component_expected);
                // Use the same resolution as the face-to-PIC conversion
                // gate.  A larger component balance is material and remains
                // a hard closure failure in VpfpIntegrator.
                if (unexplained <= interface_conversion_roundoff_floor(grid_)) {
                    local_tail_owned_roundoff_discarded += unexplained;
                } else {
                    local_tail_owned_residual += unexplained;
                    if (!local_tail_component_failure_reported) {
                        std::fprintf(stderr,
                                     "[tail-owned-component-fail] rank=%d "
                                     "ix=%d iuperp=%d jlo=%d jhi=%d "
                                     "component_number=%.17g expected_transfer=%.17g "
                                     "unexplained=%.17g floor=%.17g\n",
                                     mpi_rank, grid_.global_cell(ix), k,
                                     tail_component_start, tail_component_end,
                                     tail_component_number,
                                     tail_component_expected, unexplained,
                                     interface_conversion_roundoff_floor(grid_));
                        local_tail_component_failure_reported = true;
                    }
                }
                tail_component_active = false;
                tail_component_start = -1;
                tail_component_end = -1;
                tail_component_number = 0.0;
                tail_component_expected = 0.0;
            };
            for (int j = 0; j < nv; ++j) {
                // Capture the input cell mass before any in-place write so the
                // discrete kinetic-energy identity remains correct when
                // output == input.
                const double input_mass = input.f[upar_index(ix, j, k)];
                const double m_new =
                    input_mass
                    - upar_swept_[static_cast<size_t>(j) + 1]
                    + upar_swept_[static_cast<size_t>(j)];
                const bool tail_owned = apply_interface_sink &&
                    partition->is_tail_owned(j, k);
                if (tail_owned) {
                    if (interior) {
                        if (!tail_component_active) {
                            tail_component_active = true;
                            tail_component_start = j;
                        }
                        tail_component_end = j;
                        const double input_floor =
                            256.0 * std::numeric_limits<double>::epsilon() *
                            std::max(1.0, std::fabs(input_mass));
                        // Only physical x cells enter the conversion ledger.
                        // Guard cells are communication copies and have no
                        // independent face transfer; counting them duplicates
                        // the neighboring rank's tail balance.
                        if (!std::isfinite(input_mass) ||
                            std::fabs(input_mass) > input_floor)
                            local_finite = false;
                        const double expected_transfer =
                            j == selected_tail_cell ? selected_tail_transfer : 0.0;
                        tail_component_number += m_new;
                        tail_component_expected += expected_transfer;
                    }
                    output.f[upar_index(ix, j, k)] = 0.0;
                } else {
                    finish_tail_component();
                    output.f[upar_index(ix, j, k)] = m_new;
                }
                if (interior) {
                    const double stored = output.f[upar_index(ix, j, k)];
                    local_after += stored;
                    local_min_mass = std::min(local_min_mass, stored);
                    local_output_min = std::min(local_output_min, stored);
                    local_output_max = std::max(local_output_max, stored);
                    local_finite = local_finite && std::isfinite(stored);
                    if (!std::isfinite(stored) && local_bad_ix < 0) {
                        local_bad_ix = grid_.global_cell(ix);
                        local_bad_j = j;
                        local_bad_k = k;
                        local_bad_value = stored;
                    }
                    // Gate C (section 7.3): discrete kinetic-energy identity.
                    // Uses the final upar_swept_ actually applied to m_new and
                    // the production kinetic_energy cell moment.  Internal
                    // faces are indexed by the cell to their right; the two
                    // velocity-boundary faces are the j==0 left face and the
                    // j==nv-1 right face, matching the telescope of the flux
                    // update M_j^{n+1} = M_j - S_{j+1} + S_j.
                    const double ke = cgrid_->kinetic_energy[
                        static_cast<size_t>(j) *
                            static_cast<size_t>(nmu) + k];
                    local_delta_bulk_ke +=
                        ke * (stored - input_mass);
                    if (local_delta_ke_by_x != NULL) {
                        (*local_delta_ke_by_x)[static_cast<size_t>(ix - ng)] +=
                            ke * (stored - input_mass);
                    }
                    if (j == 0) {
                        local_left_face_energy += upar_swept_[0] * ke;
                    } else {
                        const double ke_left = cgrid_->kinetic_energy[
                            static_cast<size_t>(j - 1) *
                                static_cast<size_t>(nmu) + k];
                        local_internal_face_energy +=
                            upar_swept_[static_cast<size_t>(j)] *
                            (ke - ke_left);
                    }
                    if (j == nv - 1) {
                        local_right_face_energy +=
                            -upar_swept_[static_cast<size_t>(nv)] * ke;
                    }
                    if (tail_owned) {
                        // Energy cleared from a tail-owned Eulerian cell is a
                        // bulk-to-tail representation transfer, not ordinary
                        // internal face work; keep it separate.
                        local_interface_removed_energy += ke * m_new;
                    }
                }
            }
            finish_tail_component();
        }
    }

    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    const double local_values[13] = {
        local_before, local_after, 0.0, 0.0, 0.0, 0.0,
        local_tail_number, local_tail_energy,
        local_delta_bulk_ke, local_internal_face_energy,
        local_left_face_energy, local_right_face_energy,
        local_interface_removed_energy
    };
    double global_values[13] = {};
    MPI_Allreduce(local_values, global_values, 13, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_min_mass, &diag.minimum_cell_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    const double local_counts[2] = {
        static_cast<double>(constant_cells), static_cast<double>(linear_cells)
    };
    double global_counts[2] = { 0.0, 0.0 };
    MPI_Allreduce(local_counts, global_counts, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_departure, &diag.max_departure_cells, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    int finite_i = local_finite ? 1 : 0;
    int finite_g = 0;
    MPI_Allreduce(&finite_i, &finite_g, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    int audit_i = local_audit_valid ? 1 : 0;
    int audit_g = 0;
    MPI_Allreduce(&audit_i, &audit_g, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    int audit_code = 0;
    MPI_Allreduce(&local_audit_failure_code, &audit_code, 1, MPI_INT,
                  MPI_BOR, MPI_COMM_WORLD);
    int audit_state_i = local_audit_state_equal ? 1 : 0;
    int audit_state_g = 0;
    MPI_Allreduce(&audit_state_i, &audit_state_g, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    int parcel_failure_rank_local = local_parcel_failure_seen ? mpi_rank : mpi_size;
    int parcel_failure_rank = mpi_size;
    MPI_Allreduce(&parcel_failure_rank_local, &parcel_failure_rank, 1,
                  MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    int parcel_failure_data[4] = {
        static_cast<int>(ParcelNodeFailureReason::None), -1, -1, -1
    };
    double parcel_failure_values[4] = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    if (mpi_rank == parcel_failure_rank) {
        parcel_failure_data[0] = static_cast<int>(local_parcel_failure.reason);
        parcel_failure_data[1] = local_parcel_failure_ix;
        parcel_failure_data[2] = local_parcel_failure_face;
        parcel_failure_data[3] = local_parcel_failure_iuperp;
        parcel_failure_values[0] = local_parcel_failure.node_mass;
        parcel_failure_values[1] = local_parcel_failure.target;
        parcel_failure_values[2] = local_parcel_failure.node_sum;
        parcel_failure_values[3] = local_parcel_failure.scale;
    }
    const int parcel_failure_root =
        parcel_failure_rank < mpi_size ? parcel_failure_rank : 0;
    MPI_Bcast(parcel_failure_data, 4, MPI_INT, parcel_failure_root,
              MPI_COMM_WORLD);
    MPI_Bcast(parcel_failure_values, 4, MPI_DOUBLE, parcel_failure_root,
              MPI_COMM_WORLD);
    int bad_rank_local = local_bad_ix >= 0 ? mpi_rank : mpi_size;
    int bad_rank = mpi_size;
    MPI_Allreduce(&bad_rank_local, &bad_rank, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    int bad_data[3] = { -1, -1, -1 };
    double bad_value = std::numeric_limits<double>::quiet_NaN();
    if (mpi_rank == bad_rank) {
        bad_data[0] = local_bad_ix;
        bad_data[1] = local_bad_j;
        bad_data[2] = local_bad_k;
        bad_value = local_bad_value;
    }
    MPI_Bcast(bad_data, 3, MPI_INT,
              bad_rank < mpi_size ? bad_rank : 0, MPI_COMM_WORLD);
    MPI_Bcast(&bad_value, 1, MPI_DOUBLE,
              bad_rank < mpi_size ? bad_rank : 0, MPI_COMM_WORLD);

    diag.number_before = global_values[0];
    diag.number_after = global_values[1];
    diag.inflow_number = global_values[2];
    diag.outflow_number = global_values[3];
    diag.inflow_energy = global_values[4];
    diag.outflow_energy = global_values[5];
    diag.tail_number_loss = global_values[6];
    diag.tail_energy_loss = global_values[7];
    diag.upar_internal_face_energy_transfer = global_values[9];
    diag.upar_left_velocity_boundary_energy = global_values[10];
    diag.upar_right_velocity_boundary_energy = global_values[11];
    diag.upar_interface_energy_removed = global_values[12];
    diag.upar_discrete_energy_identity_residual =
        global_values[8] - global_values[9] - global_values[10] -
        global_values[11] + global_values[12];
    double local_interface_values[7] = {
        local_interface_face_number, local_interface_number,
        local_interface_energy, local_tail_owned_expected_transfer,
        local_tail_owned_roundoff_discarded, local_tail_owned_residual,
        local_interface_roundoff_discarded
    };
    double global_interface_values[7] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    MPI_Allreduce(local_interface_values, global_interface_values, 7,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::uint64_t local_interface_counts[2] = {
        local_interface_parcels, local_interface_nodes
    };
    std::uint64_t global_interface_counts[2] = { 0, 0 };
    MPI_Allreduce(local_interface_counts, global_interface_counts, 2,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    diag.interface_face_export_number = global_interface_values[0];
    diag.interface_parcel_number = global_interface_values[1];
    diag.interface_export_number = global_interface_values[1];
    diag.interface_export_energy = global_interface_values[2];
    diag.tail_owned_expected_transfer_number = global_interface_values[3];
    diag.tail_owned_roundoff_discarded_number = global_interface_values[4];
    diag.tail_owned_bulk_residual = global_interface_values[5];
    diag.interface_roundoff_discarded_number = global_interface_values[6];
    diag.interface_parcel_count = global_interface_counts[0];
    diag.interface_node_count = global_interface_counts[1];
    MPI_Allreduce(&local_input_min, &diag.input_min_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_input_max, &diag.input_max_mass, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_output_min, &diag.output_min_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_output_max, &diag.output_max_mass, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    diag.first_nonfinite_rank = bad_rank < mpi_size ? bad_rank : -1;
    diag.first_nonfinite_ix = bad_data[0];
    diag.first_nonfinite_iupar = bad_data[1];
    diag.first_nonfinite_iuperp = bad_data[2];
    diag.first_nonfinite_value = bad_value;
    diag.audit_valid = audit_g != 0;
    diag.audit_failure_code = audit_code;
    diag.audit_physical_state_bitwise_equal = audit_state_g != 0;
    if (parcel_failure_rank < mpi_size) {
        diag.audit_parcel_failure_reason = parcel_failure_data[0];
        diag.audit_parcel_failure_rank = parcel_failure_rank;
        diag.audit_parcel_failure_ix = parcel_failure_data[1];
        diag.audit_parcel_failure_face = parcel_failure_data[2];
        diag.audit_parcel_failure_iuperp = parcel_failure_data[3];
        diag.audit_parcel_failure_node_mass = parcel_failure_values[0];
        diag.audit_parcel_failure_target = parcel_failure_values[1];
        diag.audit_parcel_failure_node_sum = parcel_failure_values[2];
        diag.audit_parcel_failure_scale = parcel_failure_values[3];
    }
    if (interface_mode) {
        exported_flux->recompute(partition->min_conversion_energy);
        // The batch itself is rank-local.  Reduce its audit fields before
        // placing them in RemapDiagnostics; otherwise rank 0 would report
        // only the parcels owned by rank 0 while the export number above is
        // already global.
        double local_below_threshold =
            exported_flux->below_threshold_number;
        double global_below_threshold = 0.0;
        MPI_Allreduce(&local_below_threshold, &global_below_threshold, 1,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        unsigned long long local_duplicates = static_cast<unsigned long long>(
            exported_flux->duplicate_count);
        unsigned long long global_duplicates = 0;
        MPI_Allreduce(&local_duplicates, &global_duplicates, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        double local_quadrature_error = exported_flux->quadrature_error_max;
        double global_quadrature_error = 0.0;
        MPI_Allreduce(&local_quadrature_error, &global_quadrature_error, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        int local_batch_finite = exported_flux->finite ? 1 : 0;
        int global_batch_finite = 0;
        MPI_Allreduce(&local_batch_finite, &global_batch_finite, 1,
                      MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        int local_batch_nonnegative = exported_flux->nonnegative ? 1 : 0;
        int global_batch_nonnegative = 0;
        MPI_Allreduce(&local_batch_nonnegative, &global_batch_nonnegative, 1,
                      MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        diag.interface_below_threshold_number = global_below_threshold;
        diag.interface_quadrature_error_max = global_quadrature_error;
        diag.interface_duplicate_count = static_cast<std::uint64_t>(
            global_duplicates);
        if (global_batch_finite == 0) {
            diag.audit_valid = false;
            diag.audit_failure_code |= 2;
        }
        if (global_batch_nonnegative == 0) {
            diag.audit_valid = false;
            diag.audit_failure_code |= 4;
        }
        if (global_duplicates != 0) {
            diag.audit_valid = false;
            diag.audit_failure_code |= 8;
        }
        if (!std::isfinite(global_quadrature_error)) {
            diag.audit_valid = false;
            diag.audit_failure_code |= 16;
        }
    }
    const double total_cells =
        static_cast<double>(nxl) * static_cast<double>(nv) *
        static_cast<double>(nmu) * static_cast<double>(std::max(1, mpi_size));
    diag.constant_fraction =
        total_cells > 0.0 ? global_counts[0] / total_cells : 0.0;
    diag.linear_fraction =
        total_cells > 0.0 ? global_counts[1] / total_cells : 0.0;
    diag.finite = finite_g != 0;
    return diag;
}

double ConservativePpmRemap::parabola_integral(int extended_cell,
                                               double s1, double s2) const
{
    if (!(s2 > s1)) return 0.0;
    const double L = work_.left_edge[static_cast<size_t>(extended_cell)];
    const double R = work_.right_edge[static_cast<size_t>(extended_cell)];
    const double cur = work_.curvature[static_cast<size_t>(extended_cell)];
    const double ds = s2 - s1;
    const double ds2 = s2 * s2 - s1 * s1;
    const double ds3 = s2 * s2 * s2 - s1 * s1 * s1;
    // f(s) = L + (R-L) s + cur s(1-s); integral over [s1,s2] times dx.
    const double integral = L * ds
                          + (R - L) * 0.5 * ds2
                          + cur * (0.5 * ds2 - ds3 / 3.0);
    return integral * grid_.dx;
}

double ConservativePpmRemap::swept_mass(int local_face, double vx,
                                        double dt) const
{
    if (vx == 0.0 || dt <= 0.0) return 0.0;
    const int H = halo_width_;
    const int E = nx_local_ + 2 * H;
    const double x_face = static_cast<double>(H) + static_cast<double>(local_face);
    const double shift = vx * dt / grid_.dx;
    const double lo = std::min(x_face, x_face - shift);
    const double hi = std::max(x_face, x_face - shift);
    if (!(hi > lo)) return 0.0;

    int c_lo = static_cast<int>(std::floor(lo));
    int c_hi = static_cast<int>(std::ceil(hi)) - 1;
    if (c_lo < 0) c_lo = 0;
    if (c_hi >= E) c_hi = E - 1;
    double s_lo = lo - static_cast<double>(c_lo);
    double s_hi = hi - static_cast<double>(c_hi);
    s_lo = std::max(0.0, std::min(1.0, s_lo));
    s_hi = std::max(0.0, std::min(1.0, s_hi));

    double swept = 0.0;
    if (c_hi == c_lo) {
        swept = parabola_integral(c_lo, s_lo, s_hi);
    } else {
        swept = parabola_integral(c_lo, s_lo, 1.0)
              + (work_.cumulative_mass[static_cast<size_t>(c_hi)] -
                 work_.cumulative_mass[static_cast<size_t>(c_lo + 1)])
              + parabola_integral(c_hi, 0.0, s_hi);
    }
    return (vx > 0.0) ? swept : -swept;
}

void ConservativePpmRemap::fill_physical_halo(
    PhysicalSide side, const Species& input, int j_upar, int k_uperp,
    double vx, double time,
    const OpenBackgroundBoundary& boundary)
{
    const int H = halo_width_;
    const int nxl = nx_local_;
    // Physical halo: incoming reservoir/absorbing line density, or a
    // zero-gradient extension of the interior for outflow characteristics.
    if (side == PhysicalSide::LEFT) {
        if (boundary.is_incoming(PhysicalSide::LEFT, vx)) {
            const double lambda = boundary.incoming_cell_average(
                PhysicalSide::LEFT, j_upar, k_uperp, time, input);
            for (int h = 0; h < H; ++h) {
                work_.average[static_cast<size_t>(h)] = lambda;
            }
        } else {
            const double interior = work_.average[static_cast<size_t>(H)];
            for (int h = 0; h < H; ++h) {
                work_.average[static_cast<size_t>(h)] = interior;
            }
        }
    } else {
        if (boundary.is_incoming(PhysicalSide::RIGHT, vx)) {
            const double lambda = boundary.incoming_cell_average(
                PhysicalSide::RIGHT, j_upar, k_uperp, time, input);
            for (int h = 0; h < H; ++h) {
                work_.average[static_cast<size_t>(H + nxl + h)] = lambda;
            }
        } else {
            const double interior =
                work_.average[static_cast<size_t>(H + nxl - 1)];
            for (int h = 0; h < H; ++h) {
                work_.average[static_cast<size_t>(H + nxl + h)] = interior;
            }
        }
    }
}

void ConservativePpmRemap::build_slice_reconstruction(
    long long& constant_cells, long long& linear_cells,
    long long& limited_cells)
{
    const int H = halo_width_;
    const int nxl = nx_local_;
    const int E = nxl + 2 * H;
    const double dx = grid_.dx;

    // Unlimited 4th-order face values on the uniform x grid.
    work_.left_edge[0] = work_.average[0];
    for (int e = 1; e < E; ++e) {
        const int im2 = std::max(0, e - 2);
        const int ip1 = std::min(E - 1, e + 1);
        work_.left_edge[static_cast<size_t>(e)] =
            (7.0 / 12.0) * (work_.average[static_cast<size_t>(e - 1)] +
                            work_.average[static_cast<size_t>(e)])
          - (1.0 / 12.0) * (work_.average[static_cast<size_t>(im2)] +
                            work_.average[static_cast<size_t>(ip1)]);
    }
    work_.left_edge[static_cast<size_t>(E)] = work_.average[static_cast<size_t>(E - 1)];

    const double eps = std::numeric_limits<double>::epsilon();
    for (int c = 0; c < E; ++c) {
        const double a = work_.average[static_cast<size_t>(c)];
        double L = work_.left_edge[static_cast<size_t>(c)];
        double R = work_.left_edge[static_cast<size_t>(c + 1)];
        bool limited = false;

        // Colella-Woodward monotonicity constraint (no interior extrema).
        if ((R - a) * (a - L) <= 0.0) {
            const bool changed = (L != a) || (R != a);
            L = a;
            R = a;
            limited = limited || changed;
        } else {
            const double d = R - L;
            const double m = a - 0.5 * (L + R);
            if (d * m > d * d / 6.0) {
                const double new_r = 3.0 * a - 2.0 * L;
                limited = limited || (new_r != R);
                R = new_r;
            } else if (-d * m > d * d / 6.0) {
                const double new_l = 3.0 * a - 2.0 * R;
                limited = limited || (new_l != L);
                L = new_l;
            }
        }

        if (a < 0.0) {
            // Roundoff-scale negative input: record the debt and reduce the
            // reconstruction to constant in this cell only.
            const bool changed = (L != a) || (R != a);
            L = a;
            R = a;
            limited = limited || changed;
        } else {
            // Positivity of the reconstructed parabola; reduce PPM -> linear
            // locally if it dips below roundoff scale.
            double cur = 6.0 * a - 3.0 * (L + R);
            double parabola_min = std::min(L, R);
            if (cur != 0.0) {
                const double s_star = (R - L + cur) / (2.0 * cur);
                if (s_star > 0.0 && s_star < 1.0) {
                    parabola_min = std::min(
                        parabola_min,
                        L + (R - L) * s_star + cur * s_star * (1.0 - s_star));
                }
            }
            const double pos_tol =
                128.0 * eps * std::max(1.0, std::fabs(a));
            if (parabola_min < -pos_tol) {
                const double d =
                    std::min(std::min(std::fabs(L - a), std::fabs(R - a)), a);
                const double new_l = a - d;
                const double new_r = a + d;
                limited = limited || (new_l != L) || (new_r != R);
                L = new_l;
                R = new_r;
            }
        }

        work_.left_edge[static_cast<size_t>(c)] = L;
        work_.right_edge[static_cast<size_t>(c)] = R;
        work_.curvature[static_cast<size_t>(c)] = 6.0 * a - 3.0 * (L + R);

        if (c >= H && c < H + nxl) {
            const double cur = work_.curvature[static_cast<size_t>(c)];
            const bool is_constant = (L == a && R == a);
            const bool is_linear =
                !is_constant &&
                std::fabs(cur) <= 128.0 * eps * std::max(1.0, std::fabs(a));
            if (is_constant) ++constant_cells;
            else if (is_linear) ++linear_cells;
            if (limited) ++limited_cells;
        }
    }

    // Prefix cell masses for the full-cell part of departure intervals.
    work_.cumulative_mass[0] = 0.0;
    for (int c = 0; c < E; ++c) {
        work_.cumulative_mass[static_cast<size_t>(c + 1)] =
            work_.cumulative_mass[static_cast<size_t>(c)] +
            work_.average[static_cast<size_t>(c)] * dx;
    }
}

void ConservativePpmRemap::exchange_halo(const Species& input, int mpi_rank,
                                         int mpi_size, double dt)
{
    const int H = halo_width_;
    const int nxl = nx_local_;
    const int ng = ng_;
    const size_t slab = static_cast<size_t>(H) * nvmu();
    const size_t nv = static_cast<size_t>(Param::Nvmu);

    for (int h = 0; h < H; ++h) {
        for (size_t q = 0; q < nv; ++q) {
            const int j = static_cast<int>(q) / Param::Nmu;
            const int k = static_cast<int>(q) % Param::Nmu;
            send_left_halo_[static_cast<size_t>(h) * nv + q] =
                input.f[idx3(ng + h, j, k)];
            send_right_halo_[static_cast<size_t>(h) * nv + q] =
                input.f[idx3(ng + nxl - H + h, j, k)];
        }
    }

    if (mpi_size > 1) {
        const int left = mpi_rank > 0 ? mpi_rank - 1 : MPI_PROC_NULL;
        const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : MPI_PROC_NULL;
        MPI_Sendrecv(send_left_halo_.data(), static_cast<int>(slab), MPI_DOUBLE,
                     left, 7301, work_.halo.data(), static_cast<int>(slab),
                     MPI_DOUBLE, left, 7302, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_right_halo_.data(), static_cast<int>(slab), MPI_DOUBLE,
                     right, 7302, work_.halo.data() + slab,
                     static_cast<int>(slab), MPI_DOUBLE, right, 7301,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    (void)dt;
}

void ConservativePpmRemap::exchange_boundary_fluxes(int mpi_rank,
                                                    int mpi_size)
{
    if (mpi_size <= 1) return;
    const int left = mpi_rank > 0 ? mpi_rank - 1 : MPI_PROC_NULL;
    const int right = mpi_rank + 1 < mpi_size ? mpi_rank + 1 : MPI_PROC_NULL;
    MPI_Sendrecv(right_boundary_flux_.data(), static_cast<int>(nvmu()), MPI_DOUBLE,
                 right, 7303, left_boundary_flux_received_.data(),
                 static_cast<int>(nvmu()), MPI_DOUBLE, left, 7303,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

RemapDiagnostics ConservativePpmRemap::advect_x(
    const Species& input, Species& output, double dt, double time,
    const OpenBackgroundBoundary& boundary, int mpi_rank, int mpi_size,
    XFaceTransportAudit* audit)
{
    RemapDiagnostics diag = {};
    diag.finite = true;
    if (audit != NULL) {
        audit->substep_dt = dt;
        audit->bulk_number_swept_face.assign(
            static_cast<size_t>(nx_local_) + 1, 0.0);
        audit->tail_number_swept_face.assign(
            static_cast<size_t>(nx_local_) + 1, 0.0);
        audit->beam_number_swept_face.assign(
            static_cast<size_t>(nx_local_) + 1, 0.0);
    }
    const int nxl = nx_local_;
    const int ng = ng_;
    const int H = std::max(
        3, static_cast<int>(std::ceil(max_abs_vx_ * dt / grid_.dx)) + 2);
    halo_width_ = H;
    const int E = nxl + 2 * H;
    const double dx = grid_.dx;
    const size_t nv = nvmu();
    const size_t slab = static_cast<size_t>(H) * nv;

    if (!(dt > 0.0) || nxl <= 0) {
        output.f = input.f;
        boundary.fill_ghosts(output, grid_, mpi_rank, mpi_size);
        output.compute_moments();
        return diag;
    }

    work_.average.resize(static_cast<size_t>(E));
    work_.left_edge.resize(static_cast<size_t>(E) + 1);
    work_.right_edge.resize(static_cast<size_t>(E));
    work_.curvature.resize(static_cast<size_t>(E));
    work_.cumulative_mass.resize(static_cast<size_t>(E) + 1);
    work_.swept_mass.resize(static_cast<size_t>(nxl) + 1);
    work_.halo.resize(2 * slab);
    send_left_halo_.resize(slab);
    send_right_halo_.resize(slab);
    right_boundary_flux_.resize(nv);
    left_boundary_flux_received_.resize(nv);

    exchange_halo(input, mpi_rank, mpi_size, dt);

    double local_before = 0.0;
    double local_after = 0.0;
    double local_inflow_number = 0.0;
    double local_outflow_number = 0.0;
    double local_inflow_energy = 0.0;
    double local_outflow_energy = 0.0;
    double local_left_inflow_number = 0.0;
    double local_left_outflow_number = 0.0;
    double local_right_inflow_number = 0.0;
    double local_right_outflow_number = 0.0;
    double local_left_inflow_energy = 0.0;
    double local_left_outflow_energy = 0.0;
    double local_right_inflow_energy = 0.0;
    double local_right_outflow_energy = 0.0;
    double local_min_mass = std::numeric_limits<double>::infinity();
    double local_input_min = std::numeric_limits<double>::infinity();
    double local_input_max = -std::numeric_limits<double>::infinity();
    double local_output_min = std::numeric_limits<double>::infinity();
    double local_output_max = -std::numeric_limits<double>::infinity();
    long long constant_cells = 0;
    long long linear_cells = 0;
    long long limited_cells = 0;
    bool local_finite = true;
    int local_bad_ix = -1;
    int local_bad_j = -1;
    int local_bad_k = -1;
    double local_bad_value = std::numeric_limits<double>::quiet_NaN();

    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t q = idx2(j, k);
            const double vx = transport_velocity(q);
            const double ke = cgrid_->kinetic_energy[q];

            // Cell averages (line density lambda = M/dx) of the input.
            for (int i = 0; i < nxl; ++i) {
                const double m = input.f[idx3(ng + i, j, k)];
                work_.average[static_cast<size_t>(H + i)] = m / dx;
                local_before += m;
                local_min_mass = std::min(local_min_mass, m);
                local_input_min = std::min(local_input_min, m);
                local_input_max = std::max(local_input_max, m);
                local_finite = local_finite && std::isfinite(m);
                if (!std::isfinite(m) && local_bad_ix < 0) {
                    local_bad_ix = grid_.global_cell(ng + i);
                    local_bad_j = j;
                    local_bad_k = k;
                    local_bad_value = m;
                }
            }

            // Boundary halo: MPI-received cells or physical states.
            if (mpi_rank > 0) {
                for (int h = 0; h < H; ++h) {
                    work_.average[static_cast<size_t>(h)] =
                        work_.halo[static_cast<size_t>(h) * nv + q] / dx;
                }
            } else {
                fill_physical_halo(PhysicalSide::LEFT, input, j, k, vx, time,
                                   boundary);
            }
            if (mpi_rank + 1 < mpi_size) {
                for (int h = 0; h < H; ++h) {
                    work_.average[static_cast<size_t>(H + nxl + h)] =
                        work_.halo[slab + static_cast<size_t>(h) * nv + q] / dx;
                }
            } else {
                fill_physical_halo(PhysicalSide::RIGHT, input, j, k, vx, time,
                                   boundary);
            }

            build_slice_reconstruction(constant_cells, linear_cells,
                                       limited_cells);

            for (int f = 1; f <= nxl; ++f) {
                work_.swept_mass[static_cast<size_t>(f)] =
                    swept_mass(f, vx, dt);
            }
            if (mpi_rank == 0) {
                work_.swept_mass[0] = swept_mass(0, vx, dt);
                left_boundary_flux_received_[q] = work_.swept_mass[0];
            }
            if (mpi_rank + 1 < mpi_size) {
                right_boundary_flux_[q] =
                    work_.swept_mass[static_cast<size_t>(nxl)];
            }
            // Gate I (section 4.2): accumulate the finally-applied swept
            // number on the interior and right faces.  Face 0 is filled after
            // the boundary exchange below so shared MPI faces use the
            // authoritative exchanged value.
            if (audit != NULL) {
                for (int f = 1; f <= nxl; ++f) {
                    audit->bulk_number_swept_face[static_cast<size_t>(f)] +=
                        work_.swept_mass[static_cast<size_t>(f)];
                }
            }

            // Physical-boundary charge/energy ledger.
            if (mpi_rank == 0) {
                const double s0 = work_.swept_mass[0];
                if (s0 > 0.0) {
                    local_inflow_number += s0;
                    local_inflow_energy += s0 * ke;
                    local_left_inflow_number += s0;
                    local_left_inflow_energy += s0 * ke;
                } else {
                    local_outflow_number += -s0;
                    local_outflow_energy += -s0 * ke;
                    local_left_outflow_number += -s0;
                    local_left_outflow_energy += -s0 * ke;
                }
            }
            if (mpi_rank == mpi_size - 1) {
                const double sn = work_.swept_mass[static_cast<size_t>(nxl)];
                if (sn < 0.0) {
                    local_inflow_number += -sn;
                    local_inflow_energy += -sn * ke;
                    local_right_inflow_number += -sn;
                    local_right_inflow_energy += -sn * ke;
                } else {
                    local_outflow_number += sn;
                    local_outflow_energy += sn * ke;
                    local_right_outflow_number += sn;
                    local_right_outflow_energy += sn * ke;
                }
            }

            // Conservative update.  Cell 0 awaits its received left flux.
            for (int i = 1; i < nxl; ++i) {
                const double new_mass =
                    input.f[idx3(ng + i, j, k)]
                    - work_.swept_mass[static_cast<size_t>(i + 1)]
                    + work_.swept_mass[static_cast<size_t>(i)];
                output.f[idx3(ng + i, j, k)] = new_mass;
            }
            if (nxl > 0) {
                output.f[idx3(ng, j, k)] =
                    input.f[idx3(ng, j, k)]
                    - work_.swept_mass[1];
            }
        }
    }

    exchange_boundary_fluxes(mpi_rank, mpi_size);

    for (int j = 0; j < Param::Nv; ++j) {
        for (int k = 0; k < Param::Nmu; ++k) {
            const size_t q = idx2(j, k);
            const double left_flux = left_boundary_flux_received_[q];
            if (audit != NULL) {
                audit->bulk_number_swept_face[0] += left_flux;
            }
            if (nxl > 0) {
                const double m0 =
                    output.f[idx3(ng, j, k)] + left_flux;
                output.f[idx3(ng, j, k)] = m0;
                local_after += m0;
                local_min_mass = std::min(local_min_mass, m0);
                local_output_min = std::min(local_output_min, m0);
                local_output_max = std::max(local_output_max, m0);
                local_finite = local_finite && std::isfinite(m0);
                if (!std::isfinite(m0) && local_bad_ix < 0) {
                    local_bad_ix = grid_.global_cell(ng);
                    local_bad_j = j;
                    local_bad_k = k;
                    local_bad_value = m0;
                }
            }
            for (int i = 1; i < nxl; ++i) {
                const double m = output.f[idx3(ng + i, j, k)];
                local_after += m;
                local_min_mass = std::min(local_min_mass, m);
                local_output_min = std::min(local_output_min, m);
                local_output_max = std::max(local_output_max, m);
                local_finite = local_finite && std::isfinite(m);
                if (!std::isfinite(m) && local_bad_ix < 0) {
                    local_bad_ix = grid_.global_cell(ng + i);
                    local_bad_j = j;
                    local_bad_k = k;
                    local_bad_value = m;
                }
            }
        }
    }

    boundary.fill_ghosts(output, grid_, mpi_rank, mpi_size);
    output.compute_moments();

    const double local_values[14] = {
        local_before, local_after, local_inflow_number, local_outflow_number,
        local_inflow_energy, local_outflow_energy,
        local_left_inflow_number, local_left_outflow_number,
        local_right_inflow_number, local_right_outflow_number,
        local_left_inflow_energy, local_left_outflow_energy,
        local_right_inflow_energy, local_right_outflow_energy
    };
    double global_values[14] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0 };
    MPI_Allreduce(local_values, global_values, 14, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_min_mass, &diag.minimum_cell_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    const double local_counts[3] = {
        static_cast<double>(constant_cells), static_cast<double>(linear_cells),
        static_cast<double>(limited_cells)
    };
    double global_counts[3] = { 0.0, 0.0, 0.0 };
    MPI_Allreduce(local_counts, global_counts, 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    int finite_i = local_finite ? 1 : 0;
    int finite_g = 0;
    MPI_Allreduce(&finite_i, &finite_g, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    const int local_bad_rank = local_bad_ix >= 0 ? mpi_rank : mpi_size;
    int bad_rank = mpi_size;
    MPI_Allreduce(&local_bad_rank, &bad_rank, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    int bad_data[3] = { -1, -1, -1 };
    double bad_value = std::numeric_limits<double>::quiet_NaN();
    if (mpi_rank == bad_rank) {
        bad_data[0] = local_bad_ix;
        bad_data[1] = local_bad_j;
        bad_data[2] = local_bad_k;
        bad_value = local_bad_value;
    }
    MPI_Bcast(bad_data, 3, MPI_INT,
              bad_rank < mpi_size ? bad_rank : 0, MPI_COMM_WORLD);
    MPI_Bcast(&bad_value, 1, MPI_DOUBLE,
              bad_rank < mpi_size ? bad_rank : 0, MPI_COMM_WORLD);

    diag.number_before = global_values[0];
    diag.number_after = global_values[1];
    diag.inflow_number = global_values[2];
    diag.outflow_number = global_values[3];
    diag.inflow_energy = global_values[4];
    diag.outflow_energy = global_values[5];
    diag.left_inflow_number = global_values[6];
    diag.left_outflow_number = global_values[7];
    diag.right_inflow_number = global_values[8];
    diag.right_outflow_number = global_values[9];
    diag.left_inflow_energy = global_values[10];
    diag.left_outflow_energy = global_values[11];
    diag.right_inflow_energy = global_values[12];
    diag.right_outflow_energy = global_values[13];
    diag.tail_number_loss = 0.0;
    diag.tail_energy_loss = 0.0;
    diag.max_departure_cells =
        std::ceil(max_abs_vx_ * dt / grid_.dx);
    const double total_cells =
        static_cast<double>(nxl) * static_cast<double>(Param::Nvmu) *
        static_cast<double>(std::max(1, mpi_size));
    diag.constant_fraction =
        total_cells > 0.0 ? global_counts[0] / total_cells : 0.0;
    diag.linear_fraction =
        total_cells > 0.0 ? global_counts[1] / total_cells : 0.0;
    diag.limited_fraction =
        total_cells > 0.0 ? global_counts[2] / total_cells : 0.0;
    diag.limited_cell_count = static_cast<std::uint64_t>(global_counts[2]);
    MPI_Allreduce(&local_input_min, &diag.input_min_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_input_max, &diag.input_max_mass, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_output_min, &diag.output_min_mass, 1, MPI_DOUBLE,
                  MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_output_max, &diag.output_max_mass, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    diag.first_nonfinite_rank = bad_rank < mpi_size ? bad_rank : -1;
    diag.first_nonfinite_ix = bad_data[0];
    diag.first_nonfinite_iupar = bad_data[1];
    diag.first_nonfinite_iuperp = bad_data[2];
    diag.first_nonfinite_value = bad_value;
    diag.finite = finite_g != 0;
    return diag;
}
