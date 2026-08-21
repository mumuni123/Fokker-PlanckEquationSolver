#include "tail_bulk_return.h"

#include "tail_moment_constraint.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <mpi.h>
#include <set>

namespace {
// The ideal gate remains 0.5%.  On this fixed cell-centred velocity grid the
// nearest-cell representation of a continuous PIC cloud has an irreducible
// current/pressure error (about 3.3% for the audited threshold population).
// Permit no more than a small margin over that group's own native grid error,
// and never permit more than 5%.  This avoids both pathological all-defer
// behaviour and the old unconstrained ~10.8% pressure bias.
const double kRepresentationIdealTolerance = 5.0e-3;
const double kRepresentationNativeMargin = 1.05;
const double kRepresentationHardTolerance = 5.0e-2;


const int kMaxObservedVelocityAnchors = 24;

struct Moment6 {
    double v[6];
    Moment6() { for (int r = 0; r < 6; ++r) v[r] = 0.0; }
    void add(const Moment6& other)
    { for (int r = 0; r < 6; ++r) v[r] += other.v[r]; }
    void add_abs(const Moment6& other)
    { for (int r = 0; r < 6; ++r) v[r] += std::fabs(other.v[r]); }
};

struct GroupKey {
    int cell0;
    int cell1;
    GroupKey() : cell0(-1), cell1(-1) {}
    GroupKey(int left, int right) : cell0(left), cell1(right) {}
    bool operator<(const GroupKey& other) const
    {
        return cell0 != other.cell0 ? cell0 < other.cell0 : cell1 < other.cell1;
    }
};

struct RemoteGroupKey {
    int owner_rank;
    GroupKey group;
    RemoteGroupKey() : owner_rank(-1) {}
    RemoteGroupKey(int owner, const GroupKey& key)
        : owner_rank(owner), group(key) {}
    bool operator<(const RemoteGroupKey& other) const
    {
        if (owner_rank != other.owner_rank) return owner_rank < other.owner_rank;
        return group < other.group;
    }
};

struct GroupPlan {
    GroupKey key;
    std::vector<std::uint64_t> particle_ids;
    std::vector<ReturnMassRequest> local_requests;
    Moment6 cell_target[2];
    Moment6 cell_target_l1[2];
    Moment6 total_target;
    Moment6 total_target_l1;
    double guide_weight[2];
    double guide_upar[2];
    double guide_uperp[2];
    double guide_upar2[2];
    double guide_uperp2[2];
    double guide_upar_min[2];
    double guide_upar_max[2];
    double guide_uperp_min[2];
    double guide_uperp_max[2];
    std::map<std::pair<int, int>, double> observed_velocity_mass[2];
    int remote_rank;
    bool local_feasible;
    bool remote_response_received;
    bool remote_feasible;
    GroupPlan()
        : remote_rank(-1), local_feasible(true),
          remote_response_received(false), remote_feasible(false)
    {
        guide_weight[0] = guide_weight[1] = 0.0;
        guide_upar[0] = guide_upar[1] = 0.0;
        guide_uperp[0] = guide_uperp[1] = 0.0;
        guide_upar2[0] = guide_upar2[1] = 0.0;
        guide_uperp2[0] = guide_uperp2[1] = 0.0;
        guide_upar_min[0] = guide_upar_min[1] =
            std::numeric_limits<double>::infinity();
        guide_uperp_min[0] = guide_uperp_min[1] =
            std::numeric_limits<double>::infinity();
        guide_upar_max[0] = guide_upar_max[1] =
            -std::numeric_limits<double>::infinity();
        guide_uperp_max[0] = guide_uperp_max[1] =
            -std::numeric_limits<double>::infinity();
    }
};

// These are internal MPI records, never checkpoint records.  They contain
// one aggregate CIC group, so communication scales with boundary groups, not
// with the global number of tail macro-particles.
struct BoundaryGroupRequest {
    int owner_rank;
    int cell0;
    int cell1;
    int target_ix;
    double target[6];
    double target_l1[6];
    double guide_weight;
    double guide_upar;
    double guide_uperp;
    double guide_upar2;
    double guide_uperp2;
    double guide_upar_min;
    double guide_upar_max;
    double guide_uperp_min;
    double guide_uperp_max;
    int observed_velocity_count;
    int observed_velocity_j[kMaxObservedVelocityAnchors];
    int observed_velocity_k[kMaxObservedVelocityAnchors];
    double observed_velocity_mass[kMaxObservedVelocityAnchors];
};

struct BoundaryGroupResponse {
    int owner_rank;
    int cell0;
    int cell1;
    int responder_rank;
    int feasible;
};

struct BoundaryGroupCommit {
    int owner_rank;
    int cell0;
    int cell1;
    int commit;
};

int owner_of_cell(int ix, int nx_global, int mpi_size)
{
    const int base = nx_global / mpi_size;
    const int extra = nx_global % mpi_size;
    if (base <= 0 || ix < 0 || ix >= nx_global) return -1;
    const int wide = extra * (base + 1);
    if (ix < wide) return ix / (base + 1);
    return extra + (ix - wide) / base;
}

template <typename Record>
bool vector_fits_mpi_count(const std::vector<Record>& values)
{
    return values.size() <=
        static_cast<size_t>(INT_MAX / static_cast<int>(sizeof(Record)));
}

template <typename Record>
bool exchange_neighbor_records(const std::vector<Record>& send_left,
                               const std::vector<Record>& send_right,
                               std::vector<Record>& recv_left,
                               std::vector<Record>& recv_right,
                               int mpi_rank, int mpi_size, int tag_base)
{
    recv_left.clear();
    recv_right.clear();
    if (!vector_fits_mpi_count(send_left) || !vector_fits_mpi_count(send_right))
        return false;
    if (mpi_size == 1) return send_left.empty() && send_right.empty();

    const int left = mpi_rank - 1;
    const int right = mpi_rank + 1;
    const int send_left_count = static_cast<int>(send_left.size());
    const int send_right_count = static_cast<int>(send_right.size());
    int recv_left_count = 0;
    int recv_right_count = 0;
    MPI_Request count_requests[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(&send_left_count, 1, MPI_INT, left, tag_base,
                  MPI_COMM_WORLD, &count_requests[nreq++]);
        MPI_Irecv(&recv_left_count, 1, MPI_INT, left, tag_base + 1,
                  MPI_COMM_WORLD, &count_requests[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right_count, 1, MPI_INT, right, tag_base + 1,
                  MPI_COMM_WORLD, &count_requests[nreq++]);
        MPI_Irecv(&recv_right_count, 1, MPI_INT, right, tag_base,
                  MPI_COMM_WORLD, &count_requests[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, count_requests, MPI_STATUSES_IGNORE);
    if (recv_left_count < 0 || recv_right_count < 0) return false;

    recv_left.resize(static_cast<size_t>(recv_left_count));
    recv_right.resize(static_cast<size_t>(recv_right_count));
    const int send_left_bytes = send_left_count * static_cast<int>(sizeof(Record));
    const int send_right_bytes = send_right_count * static_cast<int>(sizeof(Record));
    const int recv_left_bytes = recv_left_count * static_cast<int>(sizeof(Record));
    const int recv_right_bytes = recv_right_count * static_cast<int>(sizeof(Record));
    MPI_Request data_requests[4];
    nreq = 0;
    if (left >= 0 && send_left_bytes > 0)
        MPI_Isend(send_left.data(), send_left_bytes, MPI_BYTE, left,
                  tag_base + 2, MPI_COMM_WORLD, &data_requests[nreq++]);
    if (left >= 0 && recv_left_bytes > 0)
        MPI_Irecv(recv_left.data(), recv_left_bytes, MPI_BYTE, left,
                  tag_base + 3, MPI_COMM_WORLD, &data_requests[nreq++]);
    if (right < mpi_size && send_right_bytes > 0)
        MPI_Isend(send_right.data(), send_right_bytes, MPI_BYTE, right,
                  tag_base + 3, MPI_COMM_WORLD, &data_requests[nreq++]);
    if (right < mpi_size && recv_right_bytes > 0)
        MPI_Irecv(recv_right.data(), recv_right_bytes, MPI_BYTE, right,
                  tag_base + 2, MPI_COMM_WORLD, &data_requests[nreq++]);
    if (nreq > 0) MPI_Waitall(nreq, data_requests, MPI_STATUSES_IGNORE);
    return true;
}

bool all_ranks_true(bool local)
{
    int value = local ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return value != 0;
}

double kinetic_energy(const BackgroundTailParticle& p)
{
    return Const::me * Const::c * Const::c *
        (std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz) - 1.0);
}

Moment6 particle_moments(const BackgroundTailParticle& p, double fraction)
{
    TailMoment7 m;
    tail_particle_moments(p.weight * fraction, p.x, p.ux, p.uy, p.uz, m);
    Moment6 out;
    out.v[0] = m.n; out.v[1] = m.px; out.v[2] = m.jx;
    out.v[3] = m.ke; out.v[4] = m.pixx; out.v[5] = m.piperp;
    return out;
}

Moment6 mass_moments(double mass, double upar, double uperp)
{
    double n, px, e, jx, pixx, piperp;
    mass_cell_moments(mass, upar, uperp, n, px, e, jx, pixx, piperp);
    Moment6 out;
    out.v[0] = n; out.v[1] = px; out.v[2] = jx;
    out.v[3] = e; out.v[4] = pixx; out.v[5] = piperp;
    return out;
}

double representation_relative_worst(const Moment6& candidate,
                                      const Moment6& target,
                                      const Moment6& target_l1)
{
    double worst = 0.0;
    const int rows[3] = { 2, 4, 5 };
    for (int q = 0; q < 3; ++q) {
        const int r = rows[q];
        const double scale = std::max(
            std::fabs(target.v[r]), target_l1.v[r]);
        if (!(scale > 0.0) || !std::isfinite(scale)) {
            if (candidate.v[r] != 0.0)
                return std::numeric_limits<double>::infinity();
            continue;
        }
        worst = std::max(worst,
            std::fabs(candidate.v[r] - target.v[r]) / scale);
    }
    return worst;
}

Moment6 rebuild_moments(const std::vector<std::vector<double> >& columns,
                        const std::vector<double>& weights)
{
    Moment6 rebuilt;
    if (columns.size() != weights.size()) {
        for (int r = 0; r < 6; ++r)
            rebuilt.v[r] = std::numeric_limits<double>::quiet_NaN();
        return rebuilt;
    }
    for (size_t q = 0; q < weights.size(); ++q)
        for (int r = 0; r < 6; ++r)
            rebuilt.v[r] += columns[q][static_cast<size_t>(r)] * weights[q];
    return rebuilt;
}

int nearest_index(const std::vector<double>& cells, double value)
{
    if (cells.empty() || !std::isfinite(value)) return -1;
    std::vector<double>::const_iterator it =
        std::lower_bound(cells.begin(), cells.end(), value);
    if (it == cells.begin()) return 0;
    if (it == cells.end()) return static_cast<int>(cells.size()) - 1;
    const int hi = static_cast<int>(it - cells.begin());
    return std::fabs(cells[static_cast<size_t>(hi)] - value) <
           std::fabs(cells[static_cast<size_t>(hi - 1)] - value) ? hi : hi - 1;
}

struct InvariantPoint {
    long double x;
    long double y;
    size_t column;
};

bool invariant_point_less(const InvariantPoint& a, const InvariantPoint& b)
{
    return a.x < b.x || (a.x == b.x &&
           (a.y < b.y || (a.y == b.y && a.column < b.column)));
}

long double invariant_cross(const InvariantPoint& a,
                            const InvariantPoint& b,
                            const InvariantPoint& c)
{
    return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
}

bool set_invariant_weights(
    const std::vector<std::vector<double> >& columns,
    const std::vector<double>& reference,
    const std::vector<std::pair<size_t, long double> >& fractions,
    std::vector<double>& weights, double tolerance)
{
    if (reference.size()!=3 || !(reference[0]>0.0)) return false;
    weights.assign(columns.size(),0.0);
    for (size_t q=0;q<fractions.size();++q) {
        const size_t i=fractions[q].first;
        const long double f=fractions[q].second;
        if (i>=columns.size() || columns[i].size()!=3 || f<-1.0e-14L ||
            !(columns[i][0]>0.0)) return false;
        if (f>0.0L) weights[i]+=static_cast<double>(
            f*static_cast<long double>(reference[0])/columns[i][0]);
    }
    for (size_t r=0;r<3;++r) {
        long double got=0.0L;
        for (size_t q=0;q<columns.size();++q)
            got+=static_cast<long double>(columns[q][r])*weights[q];
        const long double scale=std::max(
            1.0L,std::fabs(static_cast<long double>(reference[r])));
        if (std::fabs(got-reference[r])>
            static_cast<long double>(tolerance)*scale) return false;
    }
    return true;
}

// Deterministic exact fallback for the three representation invariants.
// With N factored out, the feasibility problem is two-dimensional in
// (u_parallel, gamma-1). First seek a local four-cell interpolation; if the
// sparse adaptive support is irregular, triangulate its convex hull.
bool solve_three_invariant_geometry(
    const std::vector<std::vector<double> >& columns,
    const std::vector<double>& reference,
    std::vector<double>& weights, double tolerance)
{
    if (reference.size()!=3 || !(reference[0]>0.0) ||
        !std::isfinite(reference[0])) return false;
    const long double mc=static_cast<long double>(Const::me)*Const::c;
    const long double mc2=mc*Const::c;
    const long double tx=reference[1]/
        (static_cast<long double>(reference[0])*mc);
    const long double ty=reference[2]/
        (static_cast<long double>(reference[0])*mc2);
    std::vector<InvariantPoint> points;
    for (size_t q=0;q<columns.size();++q) {
        if (columns[q].size()!=3 || !(columns[q][0]>0.0)) continue;
        InvariantPoint p;
        p.x=columns[q][1]/(static_cast<long double>(columns[q][0])*mc);
        p.y=columns[q][2]/(static_cast<long double>(columns[q][0])*mc2);
        p.column=q;
        if (std::isfinite(static_cast<double>(p.x)) &&
            std::isfinite(static_cast<double>(p.y))) points.push_back(p);
    }
    if (points.empty()) return false;
    std::sort(points.begin(),points.end(),invariant_point_less);
    const long double eps=1.0e-13L;

    // For each u_parallel line, find the two nearest energy levels bracketing
    // the target. Then choose the nearest lines bracketing target u_parallel.
    struct LineMix { long double x; size_t lo,hi; long double hi_fraction; };
    std::vector<LineMix> lines;
    for (size_t begin=0;begin<points.size();) {
        size_t end=begin+1;
        while (end<points.size() && points[end].x==points[begin].x) ++end;
        size_t hi=begin;
        while (hi<end && points[hi].y<ty) ++hi;
        if (hi<end && (hi>begin || std::fabs(points[hi].y-ty)<=eps)) {
            const size_t lo=(hi>begin)?hi-1:hi;
            if (points[lo].y<=ty+eps && points[hi].y>=ty-eps) {
                const long double dy=points[hi].y-points[lo].y;
                const long double hf=std::fabs(dy)<=eps?0.0L:
                    std::max(0.0L,std::min(1.0L,(ty-points[lo].y)/dy));
                LineMix line={points[begin].x,points[lo].column,
                              points[hi].column,hf};
                lines.push_back(line);
            }
        }
        begin=end;
    }
    if (!lines.empty()) {
        size_t upper=0;
        while (upper<lines.size() && lines[upper].x<tx) ++upper;
        if (upper<lines.size() &&
            (upper>0 || std::fabs(lines[upper].x-tx)<=eps)) {
            const size_t lower=(upper>0)?upper-1:upper;
            if (lines[lower].x<=tx+eps && lines[upper].x>=tx-eps) {
                const long double dx=lines[upper].x-lines[lower].x;
                const long double xf=std::fabs(dx)<=eps?0.0L:
                    std::max(0.0L,std::min(1.0L,(tx-lines[lower].x)/dx));
                std::vector<std::pair<size_t,long double> > f;
                const LineMix pair[2]={lines[lower],lines[upper]};
                const long double xw[2]={1.0L-xf,xf};
                for (int s=0;s<2;++s) {
                    f.push_back(std::make_pair(pair[s].lo,
                        xw[s]*(1.0L-pair[s].hi_fraction)));
                    f.push_back(std::make_pair(pair[s].hi,
                        xw[s]*pair[s].hi_fraction));
                }
                if (set_invariant_weights(columns,reference,f,weights,tolerance))
                    return true;
            }
        }
    }

    // Monotone-chain convex hull followed by deterministic fan
    // triangulation. Collinear interior points are intentionally discarded.
    std::vector<InvariantPoint> hull;
    for (size_t q=0;q<points.size();++q) {
        while (hull.size()>=2 && invariant_cross(
               hull[hull.size()-2],hull.back(),points[q])<=0.0L) hull.pop_back();
        hull.push_back(points[q]);
    }
    const size_t lower_size=hull.size();
    for (size_t q=points.size();q-- > 0;) {
        while (hull.size()>lower_size && invariant_cross(
               hull[hull.size()-2],hull.back(),points[q])<=0.0L) hull.pop_back();
        hull.push_back(points[q]);
    }
    if (hull.size()>1) hull.pop_back();
    if (hull.size()<3) return false;
    InvariantPoint target={tx,ty,0};
    for (size_t q=1;q+1<hull.size();++q) {
        const InvariantPoint& a=hull[0];
        const InvariantPoint& b=hull[q];
        const InvariantPoint& c=hull[q+1];
        const long double area=invariant_cross(a,b,c);
        if (!(area>0.0L)) continue;
        const long double fa=invariant_cross(target,b,c)/area;
        const long double fb=invariant_cross(a,target,c)/area;
        const long double fc=invariant_cross(a,b,target)/area;
        if (fa>=-eps && fb>=-eps && fc>=-eps) {
            std::vector<std::pair<size_t,long double> > f;
            f.push_back(std::make_pair(a.column,std::max(0.0L,fa)));
            f.push_back(std::make_pair(b.column,std::max(0.0L,fb)));
            f.push_back(std::make_pair(c.column,std::max(0.0L,fc)));
            const long double sum=f[0].second+f[1].second+f[2].second;
            if (sum>0.0L) for (size_t k=0;k<f.size();++k) f[k].second/=sum;
            if (set_invariant_weights(columns,reference,f,weights,tolerance))
                return true;
        }
    }
    return false;
}

void sample_observed_slots(
    const std::map<std::pair<int, int>, double>& source,
    const Species& bulk,
    std::vector<std::pair<int, int> >& sampled,
    std::vector<double>& sampled_mass)
{
    sampled.clear();
    sampled_mass.clear();
    if (source.size() <= static_cast<size_t>(kMaxObservedVelocityAnchors)) {
        for (std::map<std::pair<int, int>, double>::const_iterator it =
                 source.begin(); it != source.end(); ++it) {
            sampled.push_back(it->first);
            sampled_mass.push_back(it->second);
        }
        return;
    }
    // Compress the occupied velocity histogram itself.  The old equal-count
    // index buckets moved an entire bucket to one arbitrary cell centre and
    // could discard almost all perpendicular pressure while preserving N.
    // Caratheodory compression keeps all six cell-centred moments and produces
    // at most seven active supports, which fits comfortably in the fixed MPI
    // boundary record.
    std::vector<std::pair<int, int> > slots;
    std::vector<double> weights;
    std::vector<std::vector<double> > columns;
    std::vector<double> reference(6, 0.0);
    for (std::map<std::pair<int, int>, double>::const_iterator it =
             source.begin(); it != source.end(); ++it) {
        if (!(it->second > 0.0) || !std::isfinite(it->second)) continue;
        const int j = it->first.first;
        const int k = it->first.second;
        if (j < 0 || j >= Param::Nv || k < 0 || k >= Param::Nmu) continue;
        const Moment6 column = mass_moments(
            1.0, bulk.cgrid.upar_cells[static_cast<size_t>(j)],
            bulk.cgrid.uperp_cells[static_cast<size_t>(k)]);
        std::vector<double> values(6, 0.0);
        for (int r = 0; r < 6; ++r) {
            values[static_cast<size_t>(r)] = column.v[r];
            reference[static_cast<size_t>(r)] += column.v[r] * it->second;
        }
        slots.push_back(it->first);
        weights.push_back(it->second);
        columns.push_back(values);
    }
    const std::vector<double> original_weights = weights;
    const bool compressed = !weights.empty() &&
        tail_compress_moment_supports(
            columns, weights, reference,
            static_cast<size_t>(kMaxObservedVelocityAnchors), 1.0e-11);
    if (!compressed) weights = original_weights;

    std::vector<size_t> active;
    for (size_t q = 0; q < weights.size(); ++q)
        if (weights[q] > 0.0 && std::isfinite(weights[q])) active.push_back(q);
    if (active.size() > static_cast<size_t>(kMaxObservedVelocityAnchors)) {
        // This is only a defensive fallback for a failed compression.  Keep
        // the largest occupied supports rather than relocating bucket mass to
        // unrelated velocity cells.
        std::stable_sort(active.begin(), active.end(),
            [&weights](size_t a, size_t b) { return weights[a] > weights[b]; });
        active.resize(static_cast<size_t>(kMaxObservedVelocityAnchors));
        std::sort(active.begin(), active.end());
    }
    for (size_t q = 0; q < active.size(); ++q) {
        sampled.push_back(slots[active[q]]);
        sampled_mass.push_back(weights[active[q]]);
    }
}

void copy_observed_slots(
    const std::map<std::pair<int, int>, double>& source,
    std::vector<std::pair<int, int> >& slots,
    std::vector<double>& mass)
{
    slots.clear();
    mass.clear();
    slots.reserve(source.size());
    mass.reserve(source.size());
    for (std::map<std::pair<int, int>, double>::const_iterator it =
             source.begin(); it != source.end(); ++it) {
        if (!(it->second > 0.0) || !std::isfinite(it->second)) continue;
        slots.push_back(it->first);
        mass.push_back(it->second);
    }
}

// Build a representation-aware prior without turning Jx/Pixx/Piperp into
// hard constraints.  N/Px/K are polished by the exact invariant solve below;
// these iterations only choose, among feasible nonnegative representations,
// one that is closer to the PIC current and pressure tensor.
bool improve_representation_prior(
    const std::vector<std::vector<double> >& columns,
    const Moment6& target, const Moment6& target_l1,
    const std::vector<double>& input_prior, std::vector<double>& output_prior)
{
    const size_t n = columns.size();
    if (n == 0 || input_prior.size() != n) return false;
    const double mass_scale = std::max(1.0, std::fabs(target.v[0]));
    std::vector<double> row_scale(6, 1.0);
    for (int r = 0; r < 6; ++r) {
        double column_max = 0.0;
        for (size_t q = 0; q < n; ++q) {
            if (columns[q].size() != 6) return false;
            column_max = std::max(column_max,
                std::fabs(columns[q][static_cast<size_t>(r)]) * mass_scale);
        }
        const double target_scale = std::max(
            std::fabs(target.v[r]), target_l1.v[r]);
        // Moment rows carry different SI units. A dimensionless floor of 1
        // suppresses J/Pi rows whenever their cell-integrated values are
        // below one. Use the observed L1 scale, with only a roundoff-level
        // fallback tied to the available support.
        row_scale[static_cast<size_t>(r)] = std::max(
            target_scale,
            64.0 * std::numeric_limits<double>::epsilon() * column_max);
        if (!(row_scale[static_cast<size_t>(r)] > 0.0) ||
            !std::isfinite(row_scale[static_cast<size_t>(r)])) return false;
    }
    std::vector<double> y(n, 0.0);
    for (size_t q = 0; q < n; ++q) {
        if (!(input_prior[q] >= 0.0) || !std::isfinite(input_prior[q]))
            return false;
        y[q] = input_prior[q] / mass_scale;
    }
    double lipschitz = 0.0;
    for (size_t q = 0; q < n; ++q) {
        double norm2 = 0.0;
        for (int r = 0; r < 6; ++r) {
            const double a = columns[q][static_cast<size_t>(r)] * mass_scale /
                             row_scale[static_cast<size_t>(r)];
            norm2 += a * a;
        }
        lipschitz += norm2;
    }
    if (!(lipschitz > 0.0) || !std::isfinite(lipschitz)) return false;
    const double step = 0.9 / lipschitz;
    for (int iteration = 0; iteration < 256; ++iteration) {
        double residual[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        for (int r = 0; r < 6; ++r) {
            for (size_t q = 0; q < n; ++q)
                residual[r] += columns[q][static_cast<size_t>(r)] *
                    mass_scale / row_scale[static_cast<size_t>(r)] * y[q];
            residual[r] -= target.v[r] / row_scale[static_cast<size_t>(r)];
        }
        for (size_t q = 0; q < n; ++q) {
            double gradient = 0.0;
            for (int r = 0; r < 6; ++r)
                gradient += columns[q][static_cast<size_t>(r)] * mass_scale /
                    row_scale[static_cast<size_t>(r)] * residual[r];
            y[q] = std::max(0.0, y[q] - step * gradient);
        }
    }
    output_prior.assign(n, 0.0);
    double sum = 0.0;
    for (size_t q = 0; q < n; ++q) {
        output_prior[q] = y[q] * mass_scale;
        if (!std::isfinite(output_prior[q])) return false;
        sum += output_prior[q];
    }
    return sum > 0.0 && std::isfinite(sum);
}

bool solve_3x3(long double matrix[3][3], long double rhs[3],
               long double solution[3])
{
    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 3; ++row)
            if (std::fabs(matrix[row][col]) >
                std::fabs(matrix[pivot][col])) pivot = row;
        if (!(std::fabs(matrix[pivot][col]) > 1.0e-24L)) return false;
        if (pivot != col) {
            for (int q = col; q < 3; ++q)
                std::swap(matrix[col][q], matrix[pivot][q]);
            std::swap(rhs[col], rhs[pivot]);
        }
        const long double inverse = 1.0L / matrix[col][col];
        for (int q = col; q < 3; ++q) matrix[col][q] *= inverse;
        rhs[col] *= inverse;
        for (int row = 0; row < 3; ++row) {
            if (row == col) continue;
            const long double factor = matrix[row][col];
            for (int q = col; q < 3; ++q)
                matrix[row][q] -= factor * matrix[col][q];
            rhs[row] -= factor * rhs[col];
        }
    }
    for (int r = 0; r < 3; ++r) solution[r] = rhs[r];
    return true;
}

// Euclidean projection of a representation-aware prior onto the exact
// N/Px/K affine manifold with nonnegative weights.  Removing the most
// negative free variable is the standard active-set step for a projection
// onto an affine set intersected with the nonnegative orthant.
bool project_invariants_near_prior(
    const std::vector<std::vector<double> >& columns,
    const std::vector<double>& reference, const std::vector<double>& prior,
    std::vector<double>& weights, double tolerance)
{
    const size_t n = columns.size();
    if (n < 3 || reference.size() != 3 || prior.size() != n ||
        !(reference[0] > 0.0)) return false;
    const long double mass_scale = std::max(1.0, std::fabs(reference[0]));
    long double row_scale[3] = { 1.0L, 1.0L, 1.0L };
    for (int r = 0; r < 3; ++r) {
        long double largest = std::fabs(static_cast<long double>(reference[r]));
        for (size_t q = 0; q < n; ++q) {
            if (columns[q].size() != 3 || !(prior[q] >= 0.0) ||
                !std::isfinite(prior[q])) return false;
            largest = std::max(largest,
                std::fabs(static_cast<long double>(columns[q][r])) *
                mass_scale);
        }
        row_scale[r] = std::max(1.0L, largest);
    }
    std::vector<std::vector<long double> > a(
        3, std::vector<long double>(n, 0.0L));
    long double b[3];
    for (int r = 0; r < 3; ++r) {
        b[r] = static_cast<long double>(reference[r]) / row_scale[r];
        for (size_t q = 0; q < n; ++q)
            a[r][q] = static_cast<long double>(columns[q][r]) *
                      mass_scale / row_scale[r];
    }
    std::vector<long double> p(n, 0.0L);
    std::vector<long double> y(n, 0.0L);
    std::vector<unsigned char> free(n, 1);
    for (size_t q = 0; q < n; ++q)
        p[q] = static_cast<long double>(prior[q]) / mass_scale;

    for (size_t active_iteration = 0;
         active_iteration <= n - 3; ++active_iteration) {
        long double gram[3][3] = {{0.0L,0.0L,0.0L},
                                  {0.0L,0.0L,0.0L},
                                  {0.0L,0.0L,0.0L}};
        long double rhs[3] = { b[0], b[1], b[2] };
        size_t free_count = 0;
        for (size_t q = 0; q < n; ++q) {
            if (free[q] == 0) continue;
            ++free_count;
            for (int r = 0; r < 3; ++r) {
                rhs[r] -= a[r][q] * p[q];
                for (int s = 0; s < 3; ++s)
                    gram[r][s] += a[r][q] * a[s][q];
            }
        }
        if (free_count < 3) return false;
        long double lambda[3];
        if (!solve_3x3(gram, rhs, lambda)) return false;
        long double most_negative = 0.0L;
        size_t negative_index = n;
        for (size_t q = 0; q < n; ++q) {
            if (free[q] == 0) {
                y[q] = 0.0L;
                continue;
            }
            y[q] = p[q];
            for (int r = 0; r < 3; ++r) y[q] += a[r][q] * lambda[r];
            if (y[q] < most_negative) {
                most_negative = y[q];
                negative_index = q;
            }
        }
        if (negative_index != n) {
            free[negative_index] = 0;
            continue;
        }
        weights.assign(n, 0.0);
        for (size_t q = 0; q < n; ++q)
            weights[q] = static_cast<double>(std::max(0.0L, y[q]) *
                                              mass_scale);
        for (int r = 0; r < 3; ++r) {
            long double rebuilt = 0.0L;
            for (size_t q = 0; q < n; ++q)
                rebuilt += static_cast<long double>(columns[q][r]) *
                           weights[q];
            const long double scale = std::max(
                1.0L, std::fabs(static_cast<long double>(reference[r])));
            if (std::fabs(rebuilt - reference[r]) >
                static_cast<long double>(tolerance) * scale) return false;
        }
        return true;
    }
    return false;
}

bool solve_cell_projection(const Moment6& target, const Moment6& target_l1,
                            double guide_weight,
                            double guide_upar_sum, double guide_uperp_sum,
                            double guide_upar2_sum, double guide_uperp2_sum,
                            double guide_upar_min, double guide_upar_max,
                            double guide_uperp_min, double guide_uperp_max,
                            const std::vector<std::pair<int, int> >&
                                observed_velocity_slots,
                            const std::vector<double>& observed_velocity_mass,
                            bool preserve_observed_distribution,
                            const Species& bulk,
                           const HybridVelocityPartition& partition,
                           const TailBulkReturnConfig& config, int ix_global,
                           std::vector<ReturnMassRequest>& requests,
                           int& failure_kind)
{
    failure_kind=0;
    if (!(guide_weight > 0.0) || !std::isfinite(guide_weight)) {
        failure_kind=1; return false;
    }
    if (observed_velocity_slots.size() != observed_velocity_mass.size())
        { failure_kind=1; return false; }
    for (int r = 0; r < 6; ++r) if (!std::isfinite(target.v[r]))
        { failure_kind=1; return false; }
    const double mean_upar = guide_upar_sum / guide_weight;
    const double mean_uperp = guide_uperp_sum / guide_weight;
    const double sigma_upar = std::sqrt(std::max(
        0.0, guide_upar2_sum / guide_weight - mean_upar * mean_upar));
    const double sigma_uperp = std::sqrt(std::max(
        0.0, guide_uperp2_sum / guide_weight - mean_uperp * mean_uperp));
    if (!std::isfinite(mean_upar) || !std::isfinite(mean_uperp) ||
        !std::isfinite(sigma_upar) || !std::isfinite(sigma_uperp) ||
        !std::isfinite(guide_upar_min) || !std::isfinite(guide_upar_max) ||
        !std::isfinite(guide_uperp_min) || !std::isfinite(guide_uperp_max)) {
        failure_kind=1; return false;
    }
    std::vector<int> upar_anchors;
    std::vector<int> uperp_anchors;
    bool saw_representation_incompatibility = false;
    const double upar_values[5] = { guide_upar_min, mean_upar - sigma_upar,
        mean_upar, mean_upar + sigma_upar, guide_upar_max };
    const double uperp_values[5] = { guide_uperp_min,
        std::max(0.0, mean_uperp - sigma_uperp), mean_uperp,
        mean_uperp + sigma_uperp, guide_uperp_max };
    for (int q = 0; q < 5; ++q) {
        const int j = nearest_index(bulk.cgrid.upar_cells, upar_values[q]);
        const int k = nearest_index(bulk.cgrid.uperp_cells, uperp_values[q]);
        if (j >= 0) upar_anchors.push_back(j);
        if (k >= 0) uperp_anchors.push_back(k);
    }
    std::sort(upar_anchors.begin(), upar_anchors.end());
    upar_anchors.erase(std::unique(upar_anchors.begin(), upar_anchors.end()),
                       upar_anchors.end());
    std::sort(uperp_anchors.begin(), uperp_anchors.end());
    uperp_anchors.erase(std::unique(uperp_anchors.begin(), uperp_anchors.end()),
                        uperp_anchors.end());
    if (upar_anchors.empty() || uperp_anchors.empty()) return false;
    size_t max_columns=0;
    for (int radius = 1; radius <= config.max_stencil_radius; ++radius) {
        std::vector<std::vector<double> > columns;
        std::vector<std::vector<double> > full_columns;
        std::vector<double> prior;
        std::vector<std::pair<int, int> > slots;
        std::set<std::pair<int, int> > unique_slots;
        for (size_t q = 0; q < observed_velocity_slots.size(); ++q) {
            const int jc = observed_velocity_slots[q].first;
            const int kc = observed_velocity_slots[q].second;
            for (int j = std::max(0, jc - radius);
                 j <= std::min(Param::Nv - 1, jc + radius); ++j) {
                for (int k = std::max(0, kc - radius);
                     k <= std::min(Param::Nmu - 1, kc + radius); ++k) {
                    unique_slots.insert(std::make_pair(j, k));
                }
            }
        }
        for (size_t aj = 0; aj < upar_anchors.size(); ++aj) {
            for (size_t ak = 0; ak < uperp_anchors.size(); ++ak) {
                for (int j = std::max(0, upar_anchors[aj] - radius);
                     j <= std::min(Param::Nv - 1,
                                  upar_anchors[aj] + radius); ++j) {
                    for (int k = std::max(0, uperp_anchors[ak] - radius);
                         k <= std::min(Param::Nmu - 1,
                                      uperp_anchors[ak] + radius); ++k) {
                        unique_slots.insert(std::make_pair(j, k));
                    }
                }
            }
        }
        for (std::set<std::pair<int, int> >::const_iterator it =
                 unique_slots.begin(); it != unique_slots.end(); ++it) {
                const int j = it->first;
                const int k = it->second;
                const size_t id = idx2(j, k);
                if (partition.bulk_owned_cell.empty() ||
                    partition.bulk_owned_cell[id] == 0 ||
                    partition.kinetic_energy[id] >=
                        partition.min_conversion_energy) continue;
                const Moment6 column = mass_moments(
                    1.0, bulk.cgrid.upar_cells[j], bulk.cgrid.uperp_cells[k]);
                // The remap's exact invariants are number, parallel momentum
                // and kinetic energy.  Current and pressure are nonlinear
                // diagnostic moments of a continuous PIC velocity sample and
                // generally do not lie in the exact convex hull of a fixed
                // cell-centred Eulerian grid.
                std::vector<double> invariant_column(3);
                invariant_column[0] = column.v[0];
                invariant_column[1] = column.v[1];
                invariant_column[2] = column.v[3];
                columns.push_back(invariant_column);
                std::vector<double> full_column(6, 0.0);
                for (int r = 0; r < 6; ++r)
                    full_column[static_cast<size_t>(r)] = column.v[r];
                full_columns.push_back(full_column);
                slots.push_back(std::make_pair(j, k));
        }
        max_columns=std::max(max_columns,columns.size());
        if (columns.size() < 3) continue;
        // A moment-matched positive prior keeps the projected-gradient solve
        // inexpensive even though the adaptive support spans a broad cloud.
        // It is only an initial iterate; the exact N/Px/K invariant gate below
        // still decides whether this representation conversion is admissible.
        std::map<std::pair<int, int>, double> observed_prior;
        for (size_t q = 0; q < observed_velocity_slots.size(); ++q) {
            observed_prior[observed_velocity_slots[q]] +=
                observed_velocity_mass[q];
        }
        prior.assign(columns.size(), 0.0);
        double prior_sum = 0.0;
        for (size_t q = 0; q < slots.size(); ++q) {
            const int j = slots[q].first;
            const int k = slots[q].second;
            const double sj = std::max(
                sigma_upar, bulk.cgrid.upar_widths[static_cast<size_t>(j)]);
            const double sk = std::max(
                sigma_uperp, bulk.cgrid.uperp_widths[static_cast<size_t>(k)]);
            const double zj =
                (bulk.cgrid.upar_cells[static_cast<size_t>(j)] - mean_upar) / sj;
            const double zk =
                (bulk.cgrid.uperp_cells[static_cast<size_t>(k)] - mean_uperp) / sk;
            const std::pair<int, int> slot(j, k);
            const double observed = observed_prior[slot];
            prior[q] = observed + target.v[0] * 1.0e-14 *
                std::exp(-0.5 * (zj * zj + zk * zk));
            prior_sum += prior[q];
        }
        if (!(prior_sum > 0.0) || !std::isfinite(prior_sum)) {
            prior.assign(columns.size(), target.v[0] /
                                          static_cast<double>(columns.size()));
        } else {
            const double scale = target.v[0] / prior_sum;
            for (size_t q = 0; q < prior.size(); ++q) prior[q] *= scale;
        }
        // This is the error incurred by the native Eulerian representation:
        // each continuous PIC velocity is assigned to its nearest available
        // velocity cell and the cell masses are renormalized to exact N.  It
        // is a local discretization floor, not a tunable physics tolerance.
        const Moment6 nearest_cell_rebuilt =
            rebuild_moments(full_columns, prior);
        const double nearest_cell_representation_worst =
            representation_relative_worst(
                nearest_cell_rebuilt, target, target_l1);
        std::vector<double> reference(3);
        reference[0] = target.v[0];
        reference[1] = target.v[1];
        reference[2] = target.v[3];
        std::vector<double> representation_prior;
        if (!preserve_observed_distribution && improve_representation_prior(
                full_columns, target, target_l1, prior,
                representation_prior)) {
            prior.swap(representation_prior);
        }
        std::vector<double> weights;
        if (!project_invariants_near_prior(
                columns, reference, prior, weights, config.moment_tolerance) &&
            !solve_three_invariant_geometry(
                columns, reference, weights, config.moment_tolerance)) continue;
        Moment6 rebuilt;
        for (size_t q = 0; q < weights.size(); ++q) {
            if (!std::isfinite(weights[q]) || weights[q] < 0.0) return false;
            rebuilt.add(mass_moments(weights[q],
                                     bulk.cgrid.upar_cells[slots[q].first],
                                     bulk.cgrid.uperp_cells[slots[q].second]));
        }
        double worst = 0.0;
        const int invariant_rows[3] = { 0, 1, 3 };
        for (int q = 0; q < 3; ++q) {
            const int r = invariant_rows[q];
            // The H10 transaction gate is per returned CIC contribution.
            // Its denominator includes the pre-cancellation L1 amount, so a
            // nearly zero P_x target cannot produce a false failure.
            const double scale = std::max(
                1.0, std::max(std::fabs(target.v[r]), target_l1.v[r]));
            worst = std::max(worst,
                std::fabs(rebuilt.v[r] - target.v[r]) / scale);
        }
        if (worst > config.moment_tolerance) continue;
        const double representation_worst = representation_relative_worst(
            rebuilt, target, target_l1);
        const double representation_tolerance =
            tail_bulk_return_representation_tolerance(
                nearest_cell_representation_worst);
        // A representation-incompatible group remains in PIC and may be
        // reconsidered after further collisional relaxation.  Never accept a
        // large anisotropy/current change merely to reduce particle count.
        if (representation_worst > representation_tolerance) {
            saw_representation_incompatibility = true;
            continue;
        }
        for (size_t q = 0; q < weights.size(); ++q) {
            if (weights[q] == 0.0) continue;
            ReturnMassRequest req;
            req.ix_global = ix_global;
            req.iv = slots[q].first;
            req.imu = slots[q].second;
            req.mass = weights[q];
            requests.push_back(req);
        }
        return true;
    }
    failure_kind = max_columns < 3 ? 2 :
        (saw_representation_incompatibility ? 4 : 3);
    return false;
}

void update_residence(BackgroundTailParticle& particle,
                      const TailBulkReturnConfig& config)
{
    if (kinetic_energy(particle) < config.return_energy_mev * 1.0e6 * Const::qe) {
        if (particle.return_residence_steps !=
            std::numeric_limits<std::uint32_t>::max()) {
            ++particle.return_residence_steps;
        }
    } else {
        particle.return_residence_steps = 0;
    }
}

void copy_moment(const Moment6& source, double target[6])
{
    for (int r = 0; r < 6; ++r) target[r] = source.v[r];
}

Moment6 load_moment(const double source[6])
{
    Moment6 out;
    for (int r = 0; r < 6; ++r) out.v[r] = source[r];
    return out;
}

} // namespace

const char* tail_bulk_return_projection_schema()
{
    return "conservative_n_px_k_grid_aware_representation_gate_v8";
}

double tail_bulk_return_representation_tolerance(
    double nearest_cell_relative_residual)
{
    if (!(nearest_cell_relative_residual >= 0.0) ||
        !std::isfinite(nearest_cell_relative_residual))
        return kRepresentationIdealTolerance;
    return std::min(kRepresentationHardTolerance,
        std::max(kRepresentationIdealTolerance,
                 kRepresentationNativeMargin *
                     nearest_cell_relative_residual));
}

bool TailBulkReturn::apply(
    Species& bulk_trial, BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid, const HybridVelocityPartition& partition,
    long long accepted_step, int mpi_rank, int mpi_size,
    TailBulkReturnDiagnostics& diagnostics) const
{
    (void)accepted_step;
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    diagnostics = TailBulkReturnDiagnostics();
    if (!config_.enabled) return true;

    // Group every candidate which has the same CIC pair.  Each group is an
    // indivisible representation change: both cell projections and every
    // member deletion commit together or none of them do.
    std::map<GroupKey, GroupPlan> groups;
    std::set<std::uint64_t> local_candidate_ids;
    bool local_protocol_ok = true;
    std::vector<std::uint32_t> next_residence(
        tail_trial.particles.size(), 0);
    for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
        BackgroundTailParticle& particle = tail_trial.particles[p];
        BackgroundTailParticle residence_probe = particle;
        update_residence(residence_probe, config_);
        next_residence[p] = residence_probe.return_residence_steps;
        if (kinetic_energy(particle) >= config_.return_energy_mev *
                1.0e6 * Const::qe) continue;
        ++diagnostics.candidate_particles;
        if (next_residence[p] < config_.residence_steps) continue;
        ++diagnostics.resident_particles;
        if (!local_candidate_ids.insert(particle.id).second) {
            local_protocol_ok = false;
            continue;
        }
        const CellDepositWeights shape = ParticleShape1D::cell_weights(
            particle.x, grid);
        if (shape.cell0 < 0 || shape.cell1 >= grid.nx_global ||
            !(shape.w0 >= 0.0) || !(shape.w1 >= 0.0) ||
            !std::isfinite(shape.w0) || !std::isfinite(shape.w1)) {
            local_protocol_ok = false;
            continue;
        }
        const GroupKey key(shape.cell0, shape.cell1);
        GroupPlan& group = groups[key];
        group.key = key;
        group.particle_ids.push_back(particle.id);
        const int cells[2] = { shape.cell0, shape.cell1 };
        const double fractions[2] = { shape.w0, shape.w1 };
        const double uperp = std::sqrt(particle.uy * particle.uy +
                                       particle.uz * particle.uz);
        for (int side = 0; side < 2; ++side) {
            if (!(fractions[side] > 0.0)) continue;
            const Moment6 component = particle_moments(particle, fractions[side]);
            group.cell_target[side].add(component);
            group.cell_target_l1[side].add_abs(component);
            group.total_target.add(component);
            group.total_target_l1.add_abs(component);
            const double guide = particle.weight * fractions[side];
            group.guide_weight[side] += guide;
            group.guide_upar[side] += guide * particle.ux;
            group.guide_uperp[side] += guide * uperp;
            group.guide_upar2[side] += guide * particle.ux * particle.ux;
            group.guide_uperp2[side] += guide * uperp * uperp;
            group.guide_upar_min[side] =
                std::min(group.guide_upar_min[side], particle.ux);
            group.guide_upar_max[side] =
                std::max(group.guide_upar_max[side], particle.ux);
            group.guide_uperp_min[side] =
                std::min(group.guide_uperp_min[side], uperp);
            group.guide_uperp_max[side] =
                std::max(group.guide_uperp_max[side], uperp);
            const int observed_j = nearest_index(
                bulk_trial.cgrid.upar_cells, particle.ux);
            const int observed_k = nearest_index(
                bulk_trial.cgrid.uperp_cells, uperp);
            if (observed_j >= 0 && observed_k >= 0) {
                group.observed_velocity_mass[side][
                    std::make_pair(observed_j, observed_k)] += guide;
            }
            const int owner = owner_of_cell(cells[side], grid.nx_global, mpi_size);
            if (owner < 0 || (owner != mpi_rank &&
                              std::abs(owner - mpi_rank) != 1)) {
                local_protocol_ok = false;
            } else if (owner != mpi_rank) {
                if (group.remote_rank >= 0 && group.remote_rank != owner)
                    local_protocol_ok = false;
                group.remote_rank = owner;
            }
        }
    }
    diagnostics.attempted_groups =
        static_cast<std::uint64_t>(groups.size());
    if (!all_ranks_true(local_protocol_ok)) {
        diagnostics.finite = false;
        return false;
    }

    std::vector<BoundaryGroupRequest> request_left;
    std::vector<BoundaryGroupRequest> request_right;
    for (std::map<GroupKey, GroupPlan>::iterator it = groups.begin();
         it != groups.end(); ++it) {
        GroupPlan& group = it->second;
        for (int side = 0; side < 2; ++side) {
            // A particle centered exactly on a mesh point can have one zero
            // CIC component.  That component has no target moment and must
            // not create a spurious zero-moment remote projection request.
            if (!(group.guide_weight[side] > 0.0)) continue;
            const int cell = side == 0 ? group.key.cell0 : group.key.cell1;
            const int owner = owner_of_cell(cell, grid.nx_global, mpi_size);
            std::vector<std::pair<int, int> > observed_velocity_slots;
            std::vector<double> observed_velocity_mass;
            if (owner == mpi_rank) {
                // Local cells do not have a fixed-size MPI record, so retain
                // the complete occupied velocity histogram as the projection
                // prior. The solver may redistribute that prior to reduce the
                // J/Pi representation error while keeping N/Px/K exact.
                copy_observed_slots(group.observed_velocity_mass[side],
                                    observed_velocity_slots,
                                    observed_velocity_mass);
                int failure_kind=0;
                const bool ok = solve_cell_projection(
                    group.cell_target[side], group.cell_target_l1[side],
                    group.guide_weight[side],
                    group.guide_upar[side], group.guide_uperp[side],
                    group.guide_upar2[side], group.guide_uperp2[side],
                    group.guide_upar_min[side], group.guide_upar_max[side],
                    group.guide_uperp_min[side], group.guide_uperp_max[side],
                    observed_velocity_slots,
                    observed_velocity_mass,
                    false,
                    bulk_trial,
                    partition, config_, cell, group.local_requests,
                    failure_kind);
                if (!ok) {
                    if (failure_kind==1) ++diagnostics.projection_invalid_input_cells;
                    else if (failure_kind==2)
                        ++diagnostics.projection_insufficient_support_cells;
                    else if (failure_kind==4)
                        ++diagnostics.projection_representation_incompatible_cells;
                    else ++diagnostics.projection_infeasible_invariant_cells;
                }
                group.local_feasible = group.local_feasible && ok;
            } else {
                // Only the rank-boundary request must fit its fixed MPI
                // record. This affects at most the two decomposition seam
                // cells per rank, not the physical interior distribution.
                sample_observed_slots(group.observed_velocity_mass[side],
                                      bulk_trial,
                                      observed_velocity_slots,
                                      observed_velocity_mass);
                BoundaryGroupRequest request = {};
                request.owner_rank = mpi_rank;
                request.cell0 = group.key.cell0;
                request.cell1 = group.key.cell1;
                request.target_ix = cell;
                copy_moment(group.cell_target[side], request.target);
                copy_moment(group.cell_target_l1[side], request.target_l1);
                request.guide_weight = group.guide_weight[side];
                request.guide_upar = group.guide_upar[side];
                request.guide_uperp = group.guide_uperp[side];
                request.guide_upar2 = group.guide_upar2[side];
                request.guide_uperp2 = group.guide_uperp2[side];
                request.guide_upar_min = group.guide_upar_min[side];
                request.guide_upar_max = group.guide_upar_max[side];
                request.guide_uperp_min = group.guide_uperp_min[side];
                request.guide_uperp_max = group.guide_uperp_max[side];
                request.observed_velocity_count =
                    static_cast<int>(observed_velocity_slots.size());
                for (int q = 0; q < request.observed_velocity_count; ++q) {
                    request.observed_velocity_j[q] =
                        observed_velocity_slots[static_cast<size_t>(q)].first;
                    request.observed_velocity_k[q] =
                        observed_velocity_slots[static_cast<size_t>(q)].second;
                    request.observed_velocity_mass[q] =
                        observed_velocity_mass[static_cast<size_t>(q)];
                }
                if (owner < mpi_rank) request_left.push_back(request);
                else request_right.push_back(request);
            }
        }
    }
    if (!all_ranks_true(vector_fits_mpi_count(request_left) &&
                        vector_fits_mpi_count(request_right))) {
        diagnostics.finite = false;
        return false;
    }

    std::vector<BoundaryGroupRequest> received_request_left;
    std::vector<BoundaryGroupRequest> received_request_right;
    if (!exchange_neighbor_records(request_left, request_right,
                                   received_request_left, received_request_right,
                                   mpi_rank, mpi_size, 931)) {
        diagnostics.finite = false;
        return false;
    }

    std::map<RemoteGroupKey, std::vector<ReturnMassRequest> > remote_requests;
    std::vector<BoundaryGroupResponse> response_left;
    std::vector<BoundaryGroupResponse> response_right;
    bool remote_protocol_ok = true;
    const std::vector<BoundaryGroupRequest>* received[2] = {
        &received_request_left, &received_request_right };
    std::vector<BoundaryGroupResponse>* responses[2] = {
        &response_left, &response_right };
    const int source_rank[2] = { mpi_rank - 1, mpi_rank + 1 };
    for (int side = 0; side < 2; ++side) {
        for (size_t q = 0; q < received[side]->size(); ++q) {
            const BoundaryGroupRequest& request = (*received[side])[q];
            BoundaryGroupResponse response = {};
            response.owner_rank = request.owner_rank;
            response.cell0 = request.cell0;
            response.cell1 = request.cell1;
            response.responder_rank = mpi_rank;
            const GroupKey key(request.cell0, request.cell1);
            const RemoteGroupKey remote_key(request.owner_rank, key);
            const bool protocol_valid = request.owner_rank == source_rank[side] &&
                owner_of_cell(request.target_ix, grid.nx_global, mpi_size) ==
                    mpi_rank && request.target_ix >= request.cell0 &&
                request.target_ix <= request.cell1 &&
                remote_requests.find(remote_key) == remote_requests.end();
            std::vector<ReturnMassRequest> projected;
            bool feasible = protocol_valid;
            if (feasible) {
                std::vector<std::pair<int, int> > observed_velocity_slots;
                std::vector<double> observed_velocity_mass;
                if (request.observed_velocity_count < 0 ||
                    request.observed_velocity_count >
                        kMaxObservedVelocityAnchors) {
                    feasible = false;
                    remote_protocol_ok = false;
                }
                for (int r = 0; feasible &&
                     r < request.observed_velocity_count; ++r) {
                    const int j = request.observed_velocity_j[r];
                    const int k = request.observed_velocity_k[r];
                    if (j < 0 || j >= Param::Nv || k < 0 || k >= Param::Nmu) {
                        feasible = false;
                        remote_protocol_ok = false;
                        break;
                    }
                    observed_velocity_slots.push_back(std::make_pair(j, k));
                    const double mass = request.observed_velocity_mass[r];
                    if (!(mass >= 0.0) || !std::isfinite(mass)) {
                        feasible = false;
                        remote_protocol_ok = false;
                        break;
                    }
                    observed_velocity_mass.push_back(mass);
                }
                if (feasible) {
                int failure_kind=0;
                feasible = solve_cell_projection(
                    load_moment(request.target), load_moment(request.target_l1),
                    request.guide_weight,
                    request.guide_upar, request.guide_uperp,
                    request.guide_upar2, request.guide_uperp2,
                    request.guide_upar_min, request.guide_upar_max,
                    request.guide_uperp_min, request.guide_uperp_max,
                    observed_velocity_slots,
                    observed_velocity_mass,
                    false,
                    bulk_trial,
                    partition, config_, request.target_ix, projected,
                    failure_kind);
                if (!feasible) {
                    if (failure_kind==1) ++diagnostics.projection_invalid_input_cells;
                    else if (failure_kind==2)
                        ++diagnostics.projection_insufficient_support_cells;
                    else if (failure_kind==4)
                        ++diagnostics.projection_representation_incompatible_cells;
                    else ++diagnostics.projection_infeasible_invariant_cells;
                }
                }
            }
            response.feasible = feasible ? 1 : 0;
            if (!protocol_valid) {
                remote_protocol_ok = false;
            } else {
                // Keep an entry even for an infeasible projection.  The
                // owner will reply with commit=0, which is a normal deferred
                // group, not a missing MPI transaction.
                remote_requests[remote_key] = projected;
            }
            responses[side]->push_back(response);
        }
    }
    if (!all_ranks_true(remote_protocol_ok)) {
        diagnostics.finite = false;
        return false;
    }
    std::vector<BoundaryGroupResponse> received_response_left;
    std::vector<BoundaryGroupResponse> received_response_right;
    if (!all_ranks_true(vector_fits_mpi_count(response_left) &&
                        vector_fits_mpi_count(response_right))) {
        diagnostics.finite = false;
        return false;
    }
    if (!exchange_neighbor_records(response_left, response_right,
                                   received_response_left, received_response_right,
                                   mpi_rank, mpi_size, 941)) {
        diagnostics.finite = false;
        return false;
    }

    double local_request_residual = 0.0;
    bool owner_protocol_ok = true;
    const std::vector<BoundaryGroupResponse>* response_sets[2] = {
        &received_response_left, &received_response_right };
    for (int side = 0; side < 2; ++side) {
        const int responder = source_rank[side];
        for (size_t q = 0; q < response_sets[side]->size(); ++q) {
            const BoundaryGroupResponse& response = (*response_sets[side])[q];
            const GroupKey key(response.cell0, response.cell1);
            std::map<GroupKey, GroupPlan>::iterator group = groups.find(key);
            if (response.owner_rank != mpi_rank || group == groups.end() ||
                group->second.remote_rank != responder ||
                response.responder_rank != responder ||
                group->second.remote_response_received) {
                owner_protocol_ok = false;
                local_request_residual = 1.0;
                continue;
            }
            group->second.remote_response_received = true;
            group->second.remote_feasible = response.feasible != 0;
        }
    }

    std::vector<BoundaryGroupCommit> commit_left;
    std::vector<BoundaryGroupCommit> commit_right;
    for (std::map<GroupKey, GroupPlan>::iterator it = groups.begin();
         it != groups.end(); ++it) {
        GroupPlan& group = it->second;
        if (group.remote_rank < 0) continue;
        if (!group.remote_response_received) {
            owner_protocol_ok = false;
            local_request_residual = 1.0;
        }
        BoundaryGroupCommit commit = {};
        commit.owner_rank = mpi_rank;
        commit.cell0 = group.key.cell0;
        commit.cell1 = group.key.cell1;
        commit.commit = (group.local_feasible && group.remote_response_received &&
                         group.remote_feasible) ? 1 : 0;
        if (group.remote_rank < mpi_rank) commit_left.push_back(commit);
        else commit_right.push_back(commit);
    }
    if (!all_ranks_true(owner_protocol_ok)) {
        diagnostics.finite = false;
        return false;
    }
    std::vector<BoundaryGroupCommit> received_commit_left;
    std::vector<BoundaryGroupCommit> received_commit_right;
    if (!all_ranks_true(vector_fits_mpi_count(commit_left) &&
                        vector_fits_mpi_count(commit_right))) {
        diagnostics.finite = false;
        return false;
    }
    if (!exchange_neighbor_records(commit_left, commit_right,
                                   received_commit_left, received_commit_right,
                                   mpi_rank, mpi_size, 951)) {
        diagnostics.finite = false;
        return false;
    }

    std::vector<ReturnMassRequest> requests;
    std::set<std::uint64_t> remove_ids;
    Moment6 removed;
    Moment6 removed_l1;
    std::uint64_t local_committed_groups = 0;
    for (std::map<GroupKey, GroupPlan>::iterator it = groups.begin();
         it != groups.end(); ++it) {
        GroupPlan& group = it->second;
        const bool commit = group.remote_rank < 0
            ? group.local_feasible
            : group.local_feasible && group.remote_response_received &&
                group.remote_feasible;
        if (!commit) {
            ++diagnostics.deferred_infeasible_groups;
            continue;
        }
        requests.insert(requests.end(), group.local_requests.begin(),
                        group.local_requests.end());
        for (size_t q = 0; q < group.particle_ids.size(); ++q) {
            if (!remove_ids.insert(group.particle_ids[q]).second) {
                owner_protocol_ok = false;
                local_request_residual = 1.0;
            }
        }
        removed.add(group.total_target);
        removed_l1.add(group.total_target_l1);
        ++local_committed_groups;
    }

    const std::vector<BoundaryGroupCommit>* commit_sets[2] = {
        &received_commit_left, &received_commit_right };
    for (int side = 0; side < 2; ++side) {
        const int owner = source_rank[side];
        for (size_t q = 0; q < commit_sets[side]->size(); ++q) {
            const BoundaryGroupCommit& commit = (*commit_sets[side])[q];
            const RemoteGroupKey remote_key(commit.owner_rank,
                GroupKey(commit.cell0, commit.cell1));
            std::map<RemoteGroupKey, std::vector<ReturnMassRequest> >::iterator
                remote = remote_requests.find(remote_key);
            if (commit.owner_rank != owner || remote == remote_requests.end()) {
                owner_protocol_ok = false;
                local_request_residual = 1.0;
                continue;
            }
            if (commit.commit) {
                requests.insert(requests.end(), remote->second.begin(),
                                remote->second.end());
            }
            remote_requests.erase(remote);
        }
    }
    if (!remote_requests.empty()) {
        owner_protocol_ok = false;
        local_request_residual = 1.0;
    }
    double global_request_residual = 0.0;
    MPI_Allreduce(&local_request_residual, &global_request_residual, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    diagnostics.mpi_request_residual = global_request_residual;
    if (!all_ranks_true(owner_protocol_ok) ||
        global_request_residual > 1.0e-13) {
        diagnostics.finite = false;
        return false;
    }

    const bool local_valid = bulk_trial.validate_return_masses(requests);
    if (!all_ranks_true(local_valid)) {
        diagnostics.finite = false;
        return false;
    }

    // Close the representation transaction before modifying either trial
    // representation. The former ordering checked this ledger after adding
    // bulk mass and deleting Tail particles, so a residual-gate failure left
    // the caller's trial objects partially committed.
    Moment6 predicted_added;
    diagnostics.returned_number_by_x.assign(
        static_cast<size_t>(grid.nx_local), 0.0);
    for (size_t q = 0; q < requests.size(); ++q) {
        const ReturnMassRequest& request = requests[q];
        predicted_added.add(mass_moments(
            request.mass, bulk_trial.cgrid.upar_cells[request.iv],
            bulk_trial.cgrid.uperp_cells[request.imu]));
        if (request.ix_global >= grid.ix_start &&
            request.ix_global < grid.ix_start + grid.nx_local) {
            diagnostics.returned_number_by_x[
                static_cast<size_t>(request.ix_global - grid.ix_start)] +=
                request.mass;
        }
    }
    const double local_ledger[18] = {
        predicted_added.v[0], predicted_added.v[1], predicted_added.v[2],
        predicted_added.v[3], predicted_added.v[4], predicted_added.v[5],
        removed.v[0], removed.v[1], removed.v[2], removed.v[3],
        removed.v[4], removed.v[5],
        removed_l1.v[0], removed_l1.v[1], removed_l1.v[2],
        removed_l1.v[3], removed_l1.v[4], removed_l1.v[5] };
    double global_ledger[18] = { 0.0 };
    MPI_Allreduce(local_ledger, global_ledger, 18, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double* values[6] = { &diagnostics.number, &diagnostics.px,
        &diagnostics.jx_dx, &diagnostics.energy, &diagnostics.pixx_dx,
        &diagnostics.piperp_dx };
    double* residuals[6] = { &diagnostics.number_residual,
        &diagnostics.px_residual, &diagnostics.jx_residual,
        &diagnostics.energy_residual, &diagnostics.pixx_residual,
        &diagnostics.piperp_residual };
    double* differences[6] = { &diagnostics.number_difference,
        &diagnostics.px_difference, &diagnostics.jx_difference,
        &diagnostics.energy_difference, &diagnostics.pixx_difference,
        &diagnostics.piperp_difference };
    for (int r = 0; r < 6; ++r) {
        *values[r] = global_ledger[r];
        *differences[r] = global_ledger[r] - global_ledger[6 + r];
        const double l1 = std::max(
            std::max(std::fabs(global_ledger[6 + r]), global_ledger[12 + r]),
            std::numeric_limits<double>::min());
        *residuals[r] = std::fabs(global_ledger[r] - global_ledger[6 + r]) / l1;
    }
    // Exact transaction invariants for the PIC-to-Eulerian representation
    // change.  Jx and pressure residuals remain visible diagnostics and are
    // bounded by the velocity-grid convergence gates, but cannot be exact
    // constraints for arbitrary continuous particle velocities on a fixed
    // cell-centred grid.
    if (diagnostics.number_residual > config_.moment_tolerance ||
        diagnostics.px_residual > config_.moment_tolerance ||
        diagnostics.energy_residual > config_.moment_tolerance) {
        diagnostics.finite = false;
        return false;
    }

    const SpeciesReturnResult added = bulk_trial.add_return_masses(requests);
    if (!all_ranks_true(added.valid && added.complete)) {
        diagnostics.finite = false;
        return false;
    }

    std::vector<BackgroundTailParticle> kept;
    kept.reserve(tail_trial.particles.size() - remove_ids.size());
    for (size_t p = 0; p < tail_trial.particles.size(); ++p) {
        if (remove_ids.find(tail_trial.particles[p].id) == remove_ids.end()) {
            BackgroundTailParticle particle = tail_trial.particles[p];
            particle.return_residence_steps = next_residence[p];
            kept.push_back(particle);
        }
    }
    tail_trial.particles.swap(kept);

    diagnostics.committed_groups = local_committed_groups;
    diagnostics.particles_removed =
        static_cast<std::uint64_t>(remove_ids.size());
    diagnostics.committed = true;
    diagnostics.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    return true;
}
