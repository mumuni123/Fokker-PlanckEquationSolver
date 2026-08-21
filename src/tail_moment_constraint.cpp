#include "tail_moment_constraint.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Deterministic null vector of a rows x cols matrix in reduced row echelon
// form: partial pivoting by largest magnitude (ties by index), free column
// chosen by increasing index.  Returns false when no free column exists.
// Long-double arithmetic keeps |A z| near the working precision so the
// theta-amplified weight updates (section 7.5) do not drift the constraints.
bool null_vector_rref_ld(const std::vector<std::vector<long double> >& matrix,
                         long double pivot_tol, std::vector<long double>& z)
{
    const int rows = static_cast<int>(matrix.size());
    const int cols = static_cast<int>(matrix[0].size());
    if (rows <= 0 || cols <= rows) return false;
    std::vector<std::vector<long double> > m = matrix;
    std::vector<int> pivots;
    int row = 0;
    for (int col = 0; col < cols && row < rows; ++col) {
        int pr = row;
        long double pmax = std::fabs(m[static_cast<size_t>(row)]
                                          [static_cast<size_t>(col)]);
        for (int rr = row + 1; rr < rows; ++rr) {
            const long double v = std::fabs(m[static_cast<size_t>(rr)]
                                                [static_cast<size_t>(col)]);
            if (v > pmax) {
                pmax = v;
                pr = rr;
            }
        }
        if (!(pmax > pivot_tol)) continue;
        if (pr != row) std::swap(m[static_cast<size_t>(pr)],
                                 m[static_cast<size_t>(row)]);
        const long double inv = 1.0L / m[static_cast<size_t>(row)]
                                         [static_cast<size_t>(col)];
        for (int c = col; c < cols; ++c) {
            m[static_cast<size_t>(row)][static_cast<size_t>(c)] *= inv;
        }
        for (int rr = 0; rr < rows; ++rr) {
            if (rr == row) continue;
            const long double factor =
                m[static_cast<size_t>(rr)][static_cast<size_t>(col)];
            if (factor == 0.0L) continue;
            for (int c = col; c < cols; ++c) {
                m[static_cast<size_t>(rr)][static_cast<size_t>(c)] -=
                    factor * m[static_cast<size_t>(row)]
                                  [static_cast<size_t>(c)];
            }
        }
        pivots.push_back(col);
        ++row;
    }
    int free_col = -1;
    size_t pi = 0;
    for (int col = 0; col < cols; ++col) {
        if (pi < pivots.size() &&
            pivots[pi] == col) {
            ++pi;
            continue;
        }
        free_col = col;
        break;
    }
    if (free_col < 0) return false;
    z.assign(static_cast<size_t>(cols), 0.0L);
    z[static_cast<size_t>(free_col)] = 1.0L;
    for (size_t t = 0; t < pivots.size(); ++t) {
        z[static_cast<size_t>(pivots[t])] =
            -m[t][static_cast<size_t>(free_col)];
    }
    return true;
}

// Deterministic cyclic-Jacobi eigendecomposition of a symmetric n x n matrix
// (long double).  `a` is overwritten with the eigenvalues on the diagonal;
// `q` receives the orthonormal eigenvectors as columns.
void jacobi_symmetric_n_ld(std::vector<std::vector<long double> >& a,
                           std::vector<std::vector<long double> >& q, int n)
{
    q.assign(static_cast<size_t>(n),
             std::vector<long double>(static_cast<size_t>(n), 0.0L));
    for (int i = 0; i < n; ++i) {
        q[static_cast<size_t>(i)][static_cast<size_t>(i)] = 1.0L;
    }
    for (int sweep = 0; sweep < 400; ++sweep) {
        long double off = 0.0L;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                off += a[static_cast<size_t>(i)][static_cast<size_t>(j)] *
                       a[static_cast<size_t>(i)][static_cast<size_t>(j)];
            }
        }
        if (off < 1.0e-28L) break;
        for (int p = 0; p < n; ++p) {
            for (int r = p + 1; r < n; ++r) {
                const long double apq =
                    a[static_cast<size_t>(p)][static_cast<size_t>(r)];
                if (apq == 0.0L) continue;
                const long double theta =
                    (a[static_cast<size_t>(r)][static_cast<size_t>(r)] -
                     a[static_cast<size_t>(p)][static_cast<size_t>(p)]) /
                    (2.0L * apq);
                const long double t =
                    (theta >= 0.0L)
                        ? 1.0L / (theta + std::sqrt(1.0L + theta * theta))
                        : 1.0L / (theta - std::sqrt(1.0L + theta * theta));
                const long double c = 1.0L / std::sqrt(1.0L + t * t);
                const long double s = t * c;
                for (int k = 0; k < n; ++k) {
                    const long double akp =
                        a[static_cast<size_t>(k)][static_cast<size_t>(p)];
                    const long double akq =
                        a[static_cast<size_t>(k)][static_cast<size_t>(r)];
                    a[static_cast<size_t>(k)][static_cast<size_t>(p)] =
                        c * akp - s * akq;
                    a[static_cast<size_t>(k)][static_cast<size_t>(r)] =
                        s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const long double apk =
                        a[static_cast<size_t>(p)][static_cast<size_t>(k)];
                    const long double aqk =
                        a[static_cast<size_t>(r)][static_cast<size_t>(k)];
                    a[static_cast<size_t>(p)][static_cast<size_t>(k)] =
                        c * apk - s * aqk;
                    a[static_cast<size_t>(r)][static_cast<size_t>(k)] =
                        s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) {
                    const long double qkp =
                        q[static_cast<size_t>(k)][static_cast<size_t>(p)];
                    const long double qkq =
                        q[static_cast<size_t>(k)][static_cast<size_t>(r)];
                    q[static_cast<size_t>(k)][static_cast<size_t>(p)] =
                        c * qkp - s * qkq;
                    q[static_cast<size_t>(k)][static_cast<size_t>(r)] =
                        s * qkp + c * qkq;
                }
            }
        }
    }
}

// Deterministic Gaussian elimination (partial pivoting by magnitude, ties by
// index) for an n x n linear system in long double.  Returns false on a
// numerically singular matrix.
bool solve_nxn_ld(const std::vector<std::vector<long double> >& matrix,
                  const std::vector<long double>& rhs,
                  std::vector<long double>& x, int n)
{
    std::vector<std::vector<long double> > m = matrix;
    std::vector<long double> b = rhs;
    for (int col = 0; col < n; ++col) {
        int pr = col;
        long double pmax = std::fabs(m[static_cast<size_t>(col)]
                                         [static_cast<size_t>(col)]);
        for (int rr = col + 1; rr < n; ++rr) {
            const long double v = std::fabs(m[static_cast<size_t>(rr)]
                                                [static_cast<size_t>(col)]);
            if (v > pmax) {
                pmax = v;
                pr = rr;
            }
        }
        if (!(pmax > 1.0e-30L)) return false;
        if (pr != col) {
            std::swap(m[static_cast<size_t>(pr)],
                      m[static_cast<size_t>(col)]);
            std::swap(b[static_cast<size_t>(pr)],
                      b[static_cast<size_t>(col)]);
        }
        const long double inv = 1.0L / m[static_cast<size_t>(col)]
                                        [static_cast<size_t>(col)];
        for (int c = col; c < n; ++c) {
            m[static_cast<size_t>(col)][static_cast<size_t>(c)] *= inv;
        }
        b[static_cast<size_t>(col)] *= inv;
        for (int rr = 0; rr < n; ++rr) {
            if (rr == col) continue;
            const long double factor =
                m[static_cast<size_t>(rr)][static_cast<size_t>(col)];
            if (factor == 0.0L) continue;
            for (int c = col; c < n; ++c) {
                m[static_cast<size_t>(rr)][static_cast<size_t>(c)] -=
                    factor * m[static_cast<size_t>(col)]
                                  [static_cast<size_t>(c)];
            }
            b[static_cast<size_t>(rr)] -= factor * b[static_cast<size_t>(col)];
        }
    }
    x = b;
    return true;
}

} // namespace

void tail_particle_moments(double weight, double x, double ux, double uy,
                           double uz, TailMoment7& m)
{
    const double gamma = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
    m.n = weight;                                     // m^-2
    m.px = Const::me * Const::c * weight * ux;        // kg m s^-1 m^-2
    m.jx = -Const::qe * Const::c * weight * ux / gamma; // A m^-1
    m.ke = Const::me * Const::c * Const::c * weight *
           (gamma - 1.0);                             // J m^-2
    m.pixx = Const::me * Const::c * Const::c * weight *
             ux * ux / gamma;                         // J m^-2
    m.piperp = Const::me * Const::c * Const::c * weight *
               (uy * uy + uz * uz) / gamma;           // J m^-2
    m.xw = weight * x;                                // m^-1
}

bool tail_compress_moment_supports(
    const std::vector<std::vector<double> >& cols,
    std::vector<double>& weights, const std::vector<double>& ref,
    size_t max_support, double tolerance)
{
    const size_t n = weights.size();
    if (n <= max_support) return true;
    if (n == 0 || cols.size() != n || ref.empty()) return false;
    const int rows = static_cast<int>(cols[0].size());
    if (rows <= 0 || static_cast<int>(ref.size()) != rows) return false;
    for (size_t q = 0; q < n; ++q) {
        if (cols[q].size() != cols[0].size()) return false;
    }

    // Row scales for conditioning (section 7.5 step 1): each scaled row is
    // normalised by the candidate columns' own magnitude so every row has
    // O(1) entries (the per-unit-weight moments span many orders of
    // magnitude).  The null space of the scaled matrix equals the null space
    // of the unscaled one (diagonal scaling), and the final residual check
    // runs on the unscaled physical moments.
    std::vector<double> ref_abs(static_cast<size_t>(rows), 0.0);
    for (int r = 0; r < rows; ++r) {
        ref_abs[static_cast<size_t>(r)] = std::fabs(ref[static_cast<size_t>(r)]);
    }
    // The floor only prevents division by an exactly-zero row (e.g. a group
    // whose u_perp is zero everywhere).  It must stay far below the column
    // magnitudes so no constraint row is silently dropped.
    const double floor_scale = 1.0e-300;
    std::vector<double> col_max(static_cast<size_t>(rows), 0.0);
    for (size_t q = 0; q < n; ++q) {
        for (int r = 0; r < rows; ++r) {
            col_max[static_cast<size_t>(r)] = std::max(
                col_max[static_cast<size_t>(r)],
                std::fabs(cols[q][static_cast<size_t>(r)]));
        }
    }
    std::vector<double> scales(static_cast<size_t>(rows), floor_scale);
    for (int r = 0; r < rows; ++r) {
        scales[static_cast<size_t>(r)] =
            std::max(col_max[static_cast<size_t>(r)], floor_scale);
    }

    // Long-double working precision: the elimination weights can reach
    // ~1e19, so double roundoff of the null vector (|A z| ~ 1e-14) would
    // inject a large constraint drift per call.  Long double keeps the
    // drift below the final tolerance.
    std::vector<std::vector<long double> > a(
        static_cast<size_t>(rows),
        std::vector<long double>(n, 0.0L));
    for (size_t q = 0; q < n; ++q) {
        for (int r = 0; r < rows; ++r) {
            a[static_cast<size_t>(r)][q] =
                static_cast<long double>(cols[q][static_cast<size_t>(r)]) /
                scales[static_cast<size_t>(r)];
        }
    }

    // Whiten the constraint rows: the moment columns of a relativistic
    // group are near-degenerate (pixx/ke/jx scale like u_parallel for
    // u_perp << u_parallel), so the scaled matrix can be badly conditioned
    // and a null vector with |A z| ~ eps*cond gets amplified by theta.
    // Whitening via the Gram eigendecomposition makes the rows orthonormal
    // (condition 1); the physical back-transform only amplifies by the
    // largest singular value, keeping the per-call drift below the final
    // tolerance.
    std::vector<std::vector<long double> > wmat;
    {
        std::vector<std::vector<long double> > gram(
            static_cast<size_t>(rows),
            std::vector<long double>(static_cast<size_t>(rows), 0.0L));
        for (int r1 = 0; r1 < rows; ++r1) {
            for (int r2 = 0; r2 < rows; ++r2) {
                long double sum = 0.0L;
                for (size_t q = 0; q < n; ++q) {
                    sum += a[static_cast<size_t>(r1)][q] *
                           a[static_cast<size_t>(r2)][q];
                }
                gram[static_cast<size_t>(r1)][static_cast<size_t>(r2)] = sum;
            }
        }
        std::vector<std::vector<long double> > eigvec;
        jacobi_symmetric_n_ld(gram, eigvec, rows);
        // W = Lambda^{-1/2} Q^T  (rows of W are the whitening directions).
        wmat.assign(static_cast<size_t>(rows),
                    std::vector<long double>(static_cast<size_t>(rows), 0.0L));
        for (int r = 0; r < rows; ++r) {
            const long double lambda =
                std::max(gram[static_cast<size_t>(r)][static_cast<size_t>(r)],
                         1.0e-300L);
            const long double inv_sqrt = 1.0L / std::sqrt(lambda);
            for (int i = 0; i < rows; ++i) {
                wmat[static_cast<size_t>(r)][static_cast<size_t>(i)] =
                    inv_sqrt * eigvec[static_cast<size_t>(i)]
                                    [static_cast<size_t>(r)];
            }
        }
        std::vector<std::vector<long double> > b(
            static_cast<size_t>(rows),
            std::vector<long double>(n, 0.0L));
        for (int r = 0; r < rows; ++r) {
            for (size_t q = 0; q < n; ++q) {
                long double sum = 0.0L;
                for (int i = 0; i < rows; ++i) {
                    sum += wmat[static_cast<size_t>(r)][static_cast<size_t>(i)] *
                           a[static_cast<size_t>(i)][q];
                }
                b[static_cast<size_t>(r)][q] = sum;
            }
        }
        a.swap(b);
    }

    const long double pivot_tol = 1.0e-15L;
    const long double weight_floor =
        std::max(ref[0], 1.0) * 1.0e-13L + 1.0e-300L;
    std::vector<size_t> active;
    for (size_t q = 0; q < n; ++q) {
        if (weights[q] > 0.0) active.push_back(q);
    }
    std::vector<long double> wl(weights.begin(), weights.end());

    // Caratheodory reduction only needs rows + 1 active columns to obtain a
    // null direction.  The previous implementation rebuilt an RREF over the
    // entire active set while removing at least one support per iteration.
    // That is cubic in the group size and can leave one MPI rank spending
    // minutes in flux conversion while every other rank waits in the next
    // particle migration.  Restrict each elimination to a deterministic
    // block of rows + 1 columns.  The update is still an exact null-space
    // update of the full moment system (all other weights are unchanged), so
    // this changes only the work complexity, not the constrained moments.
    const size_t reduction_width = static_cast<size_t>(rows) + 1;
    size_t guard = n + 8;
    while (active.size() > max_support && guard-- > 0) {
        const size_t block_size = std::min(active.size(), reduction_width);
        if (block_size <= static_cast<size_t>(rows)) return false;
        const size_t block_begin = active.size() - block_size;
        std::vector<std::vector<long double> > sub(
            static_cast<size_t>(rows),
            std::vector<long double>(block_size, 0.0L));
        for (size_t c = 0; c < block_size; ++c) {
            for (int r = 0; r < rows; ++r) {
                sub[static_cast<size_t>(r)][c] =
                    a[static_cast<size_t>(r)][active[block_begin + c]];
            }
        }
        std::vector<long double> z;
        if (!null_vector_rref_ld(sub, pivot_tol, z)) return false;
        long double zmax = 0.0L;
        for (size_t c = 0; c < z.size(); ++c) {
            zmax = std::max(zmax, std::fabs(z[c]));
        }
        if (!(zmax > pivot_tol)) return false;
        bool has_positive = false;
        for (size_t c = 0; c < z.size(); ++c) {
            if (z[c] > 0.0L) has_positive = true;
        }
        if (!has_positive) {
            for (size_t c = 0; c < z.size(); ++c) z[c] = -z[c];
        }
        long double theta = -1.0L;
        for (size_t c = 0; c < z.size(); ++c) {
            if (z[c] > pivot_tol * zmax) {
                const long double ratio =
                    wl[active[block_begin + c]] / z[c];
                if (theta < 0.0L || ratio < theta) theta = ratio;
            }
        }
        if (!(theta > 0.0L)) return false;
        for (size_t c = 0; c < block_size; ++c) {
            wl[active[block_begin + c]] -= theta * z[c];
        }
        const size_t active_before = active.size();
        // Only the selected block changed.  Remove its exhausted entries in
        // place.  Selecting the block at the vector tail makes compaction
        // fixed-width as well, rather than shifting O(n) indices per step.
        size_t write = block_begin;
        for (size_t c = 0; c < block_size; ++c) {
            const size_t index = active[block_begin + c];
            if (wl[index] > weight_floor) {
                active[write++] = index;
            } else {
                wl[index] = 0.0L;
            }
        }
        active.resize(write);
        if (active.size() == active_before) return false;  // no progress
    }
    if (guard == 0 && active.size() > max_support) return false;

    // Polish (section 7.5 step 6): the sequential theta-elimination can
    // drift the constraints (theta can be many orders above the
    // null-vector residual), so its weights are used only to select the
    // support set.  The final weights are re-solved directly against the
    // moments through the whitened normal equations (condition ~1), then
    // projected into the nonnegative solution family along the null
    // direction.
    if (active.empty()) return false;
    {
        std::vector<long double> ref_scaled(static_cast<size_t>(rows), 0.0L);
        for (int r = 0; r < rows; ++r) {
            ref_scaled[static_cast<size_t>(r)] =
                static_cast<long double>(ref[static_cast<size_t>(r)]) /
                scales[static_cast<size_t>(r)];
        }
        std::vector<long double> target(static_cast<size_t>(rows), 0.0L);
        for (int r = 0; r < rows; ++r) {
            long double sum = 0.0L;
            for (int i = 0; i < rows; ++i) {
                sum += wmat[static_cast<size_t>(r)][static_cast<size_t>(i)] *
                       ref_scaled[static_cast<size_t>(i)];
            }
            target[static_cast<size_t>(r)] = sum;
        }
        std::vector<std::vector<long double> > c(
            static_cast<size_t>(rows),
            std::vector<long double>(static_cast<size_t>(rows), 0.0L));
        for (int r1 = 0; r1 < rows; ++r1) {
            for (int r2 = 0; r2 < rows; ++r2) {
                long double sum = 0.0L;
                for (size_t cc = 0; cc < active.size(); ++cc) {
                    sum += a[static_cast<size_t>(r1)][active[cc]] *
                           a[static_cast<size_t>(r2)][active[cc]];
                }
                c[static_cast<size_t>(r1)][static_cast<size_t>(r2)] = sum;
            }
        }
        std::vector<long double> y;
        bool polished = solve_nxn_ld(c, target, y, rows);
        std::vector<long double> x(active.size(), 0.0L);
        if (polished) {
            for (size_t cc = 0; cc < active.size(); ++cc) {
                long double sum = 0.0L;
                for (int r = 0; r < rows; ++r) {
                    sum += a[static_cast<size_t>(r)][active[cc]] *
                           y[static_cast<size_t>(r)];
                }
                x[cc] = sum;
            }
        }

        // Nonnegativity: if the min-norm solution has negative weights,
        // move along the (single) null direction into the feasible interval.
        const long double weight_tol =
            std::max(static_cast<long double>(ref[0]), 1.0L) * 1.0e-12L;
        bool feasible = polished;
        if (polished) {
            for (size_t cc = 0; cc < x.size(); ++cc) {
                if (x[cc] < -weight_tol) feasible = false;
            }
        }
        if (polished && !feasible) {
            std::vector<std::vector<long double> > sub(
                static_cast<size_t>(rows),
                std::vector<long double>(active.size(), 0.0L));
            for (size_t cc = 0; cc < active.size(); ++cc) {
                for (int r = 0; r < rows; ++r) {
                    sub[static_cast<size_t>(r)][cc] =
                        a[static_cast<size_t>(r)][active[cc]];
                }
            }
            std::vector<long double> z;
            if (null_vector_rref_ld(sub, 1.0e-15L, z)) {
                long double t_lo = -1.0e300L;
                long double t_hi = 1.0e300L;
                bool ok = true;
                for (size_t cc = 0; cc < z.size(); ++cc) {
                    if (z[cc] > 0.0L) {
                        t_lo = std::max(t_lo, -x[cc] / z[cc]);
                    } else if (z[cc] < 0.0L) {
                        t_hi = std::min(t_hi, -x[cc] / z[cc]);
                    } else if (x[cc] < 0.0L) {
                        ok = false;
                    }
                }
                if (ok && t_lo <= t_hi) {
                    const long double t_center = 0.5L * (t_lo + t_hi);
                    for (size_t cc = 0; cc < z.size(); ++cc) {
                        x[cc] += t_center * z[cc];
                    }
                    feasible = true;
                }
            }
        }
        // Rank-deficient groups are common when many face parcels reuse the
        // same quadrature nodes.  A singular full-row polish is not a failed
        // conservative reduction: the long-double null updates above already
        // preserve every represented row.  Fall back to those weights and
        // let the unscaled physical residual gate below make the decision.
        polished = polished && feasible;
        for (size_t q = 0; q < n; ++q) weights[q] = 0.0;
        for (size_t cc = 0; cc < active.size(); ++cc) {
            long double selected = polished ? x[cc] : wl[active[cc]];
            if (selected < 0.0L) selected = 0.0L;
            weights[active[cc]] = static_cast<double>(selected);
        }
    }

    // Unscaled physical residual check (section 7.5 step 6).
    std::vector<double> got(static_cast<size_t>(rows), 0.0);
    for (size_t q = 0; q < n; ++q) {
        if (!(weights[q] > 0.0)) continue;
        for (int r = 0; r < rows; ++r) {
            got[static_cast<size_t>(r)] +=
                cols[q][static_cast<size_t>(r)] * weights[q];
        }
    }
    for (int r = 0; r < rows; ++r) {
        const double residual =
            std::fabs(got[static_cast<size_t>(r)] -
                      ref[static_cast<size_t>(r)]);
        if (residual > tolerance * std::max(1.0, ref_abs[static_cast<size_t>(r)])) {
            return false;
        }
    }
    return true;
}

bool tail_solve_nonnegative_moment_weights(
    const std::vector<std::vector<double> >& cols,
    const std::vector<double>& ref, const std::vector<double>& prior,
    std::vector<double>& weights, double tolerance)
{
    const size_t n = cols.size();
    if (n == 0 || prior.size() != n || ref.empty()) return false;
    const size_t rows = ref.size();
    for (size_t q = 0; q < n; ++q)
        if (cols[q].size() != rows || !(prior[q] >= 0.0)) return false;
    const double mass_scale = std::max(1.0, std::fabs(ref[0]));
    std::vector<double> row_scale(rows, 1.0);
    for (size_t r = 0; r < rows; ++r) {
        double column_max = 0.0;
        for (size_t q = 0; q < n; ++q)
            column_max = std::max(column_max, std::fabs(cols[q][r]) * mass_scale);
        row_scale[r] = std::max(std::fabs(ref[r]), column_max);
    }
    std::vector<double> y(n, 0.0);
    for (size_t q = 0; q < n; ++q) y[q] = prior[q] / mass_scale;
    // A conservative bound for the gradient Lipschitz constant of A^T A.
    double lipschitz = 0.0;
    for (size_t q = 0; q < n; ++q) {
        double norm2 = 0.0;
        for (size_t r = 0; r < rows; ++r) {
            const double a = cols[q][r] * mass_scale / row_scale[r];
            norm2 += a * a;
        }
        lipschitz += norm2;
    }
    if (!(lipschitz > 0.0) || !std::isfinite(lipschitz)) return false;
    const double step = 0.9 / lipschitz;
    const double regularization = 1.0e-14;
    for (int iteration = 0; iteration < 12000; ++iteration) {
        std::vector<double> residual(rows, 0.0);
        for (size_t r = 0; r < rows; ++r) {
            for (size_t q = 0; q < n; ++q)
                residual[r] += cols[q][r] * mass_scale / row_scale[r] * y[q];
            residual[r] -= ref[r] / row_scale[r];
        }
        double residual_linf = 0.0;
        for (size_t r = 0; r < rows; ++r)
            residual_linf = std::max(residual_linf, std::fabs(residual[r]));
        if (residual_linf <= tolerance * 0.1) break;
        for (size_t q = 0; q < n; ++q) {
            double gradient = 0.0;
            for (size_t r = 0; r < rows; ++r)
                gradient += cols[q][r] * mass_scale / row_scale[r] * residual[r];
            gradient += regularization * (y[q] - prior[q] / mass_scale);
            y[q] = std::max(0.0, y[q] - step * gradient);
        }
    }
    weights.assign(n, 0.0);
    for (size_t q = 0; q < n; ++q) weights[q] = y[q] * mass_scale;
    for (size_t r = 0; r < rows; ++r) {
        long double value = 0.0L;
        for (size_t q = 0; q < n; ++q)
            value += static_cast<long double>(cols[q][r]) * weights[q];
        const double relative = std::fabs(static_cast<double>(value) - ref[r]) /
                                std::max(1.0, std::fabs(ref[r]));
        if (!std::isfinite(relative) || relative > tolerance) return false;
    }
    return true;
}
