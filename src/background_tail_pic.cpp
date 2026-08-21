#include "background_tail_pic.h"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace {

inline double tail_gamma(const BackgroundTailParticle& p)
{
    return std::sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
}

inline double tail_vx(const BackgroundTailParticle& p)
{
    return Const::c * p.ux / tail_gamma(p);
}

inline double tail_kinetic_energy_per_particle(
    const BackgroundTailParticle& p)
{
    return (tail_gamma(p) - 1.0) * Const::me * Const::c * Const::c;
}

inline int carrier_byte_count(int carrier_count)
{
    return static_cast<int>(static_cast<size_t>(carrier_count) *
                            sizeof(TailDriftCarrier));
}

} // namespace

CellDepositWeights ParticleShape1D::cell_weights(double x,
                                                 const SpatialGrid& sg)
{
    // Identical to the BeamPIC CIC density deposit (s = x/dx - 0.5).
    const double s = x * (1.0 / sg.dx) - 0.5;
    const int i0 = static_cast<int>(std::floor(s));
    const double frac = s - static_cast<double>(i0);
    CellDepositWeights w;
    w.cell0 = i0;
    w.cell1 = i0 + 1;
    w.w0 = 1.0 - frac;
    w.w1 = frac;
    return w;
}

FaceGatherWeights ParticleShape1D::shifted_face_weights(double x,
                                                        const SpatialGrid& sg)
{
    // Identical to the BeamPIC shifted CIC face gather: the nearest face
    // receives 1-|shift| and the adjacent face |shift|.  Applied directly to
    // Ex_face (section 6.3 constraint 1).
    const double face_coordinate = x * (1.0 / sg.dx);
    const int center_global =
        static_cast<int>(std::floor(face_coordinate + 0.5));
    const double shift = face_coordinate - static_cast<double>(center_global);
    FaceGatherWeights w;
    w.face0 = center_global;
    w.face1 = center_global;
    w.w0 = 1.0 - std::fabs(shift);
    w.w1 = 0.0;
    if (shift < 0.0) {
        w.face1 = center_global - 1;
        w.w1 = -shift;
    } else if (shift > 0.0) {
        w.face1 = center_global + 1;
        w.w1 = shift;
    }
    return w;
}

namespace {

// Shape-difference rule for one straight sub-segment [a,b]: adds
// q_e*weight*(G_j(b)-G_j(a)) to charge_face for every local face j, where
// G_j(x) is the CIC cumulative weight strictly below face j.
void deposit_sub_segment(double a, double b, double weight,
                         const SpatialGrid& sg,
                         std::vector<double>& charge_face)
{
    if (a == b) return;
    const CellDepositWeights wa = ParticleShape1D::cell_weights(a, sg);
    const CellDepositWeights wb = ParticleShape1D::cell_weights(b, sg);
    const int jmin = std::min(wa.cell0, wb.cell0) + 1;
    const int jmax = std::max(wa.cell1, wb.cell1) + 1;
    for (int j = jmin; j <= jmax; ++j) {
        const double ga =
            (j <= wa.cell0) ? 0.0
            : ((j == wa.cell0 + 1) ? wa.w0 : 1.0);
        const double gb =
            (j <= wb.cell0) ? 0.0
            : ((j == wb.cell0 + 1) ? wb.w0 : 1.0);
        const double dg = gb - ga;
        if (dg == 0.0) continue;
        const int il = j - sg.ix_start;
        if (il < 0 || il > sg.nx_local) continue;
        charge_face[static_cast<size_t>(il)] += Const::qe * weight * dg;
    }
}

} // namespace

void ChargeConservingTrajectory1D::deposit_charge_delta(
    double x0, double x1, double weight, const SpatialGrid& sg,
    std::vector<double>& charge_face)
{
    if (x0 == x1 || weight == 0.0) return;
    // Split multi-cell paths at every interior cell face (section 6.3
    // constraint 3).  For CIC the split is exact (the shape-difference rule
    // telescopes), and it keeps the structure required by higher-order
    // shapes later.
    const double inv_dx = 1.0 / sg.dx;
    const double eps = 1.0e-12 * std::max(
        1.0, std::fabs(x0 * inv_dx) + std::fabs(x1 * inv_dx) + 1.0);
    const long long max_segments =
        static_cast<long long>(sg.nx_global) + 8;
    double a = x0;
    long long segments = 0;
    while (a != x1 && segments < max_segments) {
        double b = x1;
        const double sa = a * inv_dx;
        if (x1 > a) {
            const double next = (std::floor(sa + eps) + 1.0) * sg.dx;
            if (next < x1 - eps * sg.dx) b = next;
        } else {
            const double next = std::floor(sa - eps) * sg.dx;
            if (next > x1 + eps * sg.dx) b = next;
        }
        if (b == a) b = x1;
        deposit_sub_segment(a, b, weight, sg, charge_face);
        a = b;
        ++segments;
    }
    if (a != x1) {
        // Pathological float accumulation guard: close the last piece.
        deposit_sub_segment(a, x1, weight, sg, charge_face);
    }
}

void ChargeConservingTrajectory1D::deposit_segment(
    double x0, double x1, double weight, double dt, const SpatialGrid& sg,
    FaceCurrentAccumulator& current)
{
    if (!(dt > 0.0) || x0 == x1) return;
    if (current.current_face_x.size() !=
        static_cast<size_t>(sg.nx_local) + 1) {
        current.init(sg);
    }
    const double inv_dt = 1.0 / dt;
    std::vector<double> scratch;
    scratch.assign(static_cast<size_t>(sg.nx_local) + 1, 0.0);
    deposit_charge_delta(x0, x1, weight, sg, scratch);
    for (size_t i = 0; i < scratch.size(); ++i) {
        current.current_face_x[i] += scratch[i] * inv_dt;
    }
}

BackgroundTailPIC::BackgroundTailPIC()
    : step_dt_(0.0),
      truncation_shape_left_(0.0),
      truncation_shape_right_(0.0),
      deposit_shape_left_(0.0),
      deposit_shape_right_(0.0),
      deposit_shape_step_start_left_(0.0),
      deposit_shape_step_start_right_(0.0),
      max_abs_u_(0.0),
      max_kinetic_energy_(0.0),
      id_counter_(0),
      collision_rng_seed_(0x9e3779b97f4a7c15ULL),
      migration_failed_(false),
      last_migration_seconds_(0.0)
{}

void BackgroundTailPIC::init(const SpatialGrid& sg)
{
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    particles.clear();
    density.assign(nxl, 0.0);
    current_face_x.assign(nxl + 1, 0.0);
    density_step_start.assign(nxl, 0.0);
    trajectory_charge_face_x_.assign(nxl + 1, 0.0);
    truncated_in_domain_.assign(nxl, 0.0);
    keep_.clear();
    send_left_.clear();
    send_right_.clear();
    recv_left_.clear();
    recv_right_.clear();
    step_dt_ = 0.0;
    outflow_ = TailOutflowLedger();
    truncation_shape_left_ = 0.0;
    truncation_shape_right_ = 0.0;
    deposit_shape_left_ = 0.0;
    deposit_shape_right_ = 0.0;
    deposit_shape_step_start_left_ = 0.0;
    deposit_shape_step_start_right_ = 0.0;
    max_abs_u_ = 0.0;
    max_kinetic_energy_ = 0.0;
    id_counter_ = 0;
    collision_rng_seed_ = 0x9e3779b97f4a7c15ULL;
    migration_failed_ = false;
    last_migration_seconds_ = 0.0;
}

void BackgroundTailPIC::clear()
{
    particles.clear();
    if (!density.empty()) std::fill(density.begin(), density.end(), 0.0);
    if (!current_face_x.empty())
        std::fill(current_face_x.begin(), current_face_x.end(), 0.0);
    if (!density_step_start.empty())
        std::fill(density_step_start.begin(), density_step_start.end(), 0.0);
    if (!trajectory_charge_face_x_.empty())
        std::fill(trajectory_charge_face_x_.begin(),
                  trajectory_charge_face_x_.end(), 0.0);
    if (!truncated_in_domain_.empty())
        std::fill(truncated_in_domain_.begin(),
                  truncated_in_domain_.end(), 0.0);
    step_dt_ = 0.0;
    outflow_ = TailOutflowLedger();
    truncation_shape_left_ = 0.0;
    truncation_shape_right_ = 0.0;
    deposit_shape_left_ = 0.0;
    deposit_shape_right_ = 0.0;
    deposit_shape_step_start_left_ = 0.0;
    deposit_shape_step_start_right_ = 0.0;
    max_abs_u_ = 0.0;
    max_kinetic_energy_ = 0.0;
    migration_failed_ = false;
    last_migration_seconds_ = 0.0;
}

std::uint64_t BackgroundTailPIC::next_particle_id(int mpi_rank)
{
    // Section 6.4: high bits carry the creating MPI rank, low bits the
    // per-rank accepted conversion counter.
    const std::uint64_t high =
        static_cast<std::uint64_t>(static_cast<unsigned int>(mpi_rank))
        << 32;
    const std::uint64_t low = id_counter_ & 0xFFFFFFFFu;
    ++id_counter_;
    return high | low;
}

void BackgroundTailPIC::begin_step(const SpatialGrid& sg, double dt)
{
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    density_step_start = density;
    // Capture the previous deposit's boundary shape ledgers (the step-start
    // outside-share of in-domain particles) before the per-step reset.
    deposit_shape_step_start_left_ = deposit_shape_left_;
    deposit_shape_step_start_right_ = deposit_shape_right_;
    step_dt_ = dt;
    if (trajectory_charge_face_x_.size() != nxl + 1)
        trajectory_charge_face_x_.assign(nxl + 1, 0.0);
    else
        std::fill(trajectory_charge_face_x_.begin(),
                  trajectory_charge_face_x_.end(), 0.0);
    if (truncated_in_domain_.size() != nxl)
        truncated_in_domain_.assign(nxl, 0.0);
    else
        std::fill(truncated_in_domain_.begin(),
                  truncated_in_domain_.end(), 0.0);
    if (current_face_x.size() != nxl + 1)
        current_face_x.assign(nxl + 1, 0.0);
    else
        std::fill(current_face_x.begin(), current_face_x.end(), 0.0);
    outflow_ = TailOutflowLedger();
    truncation_shape_left_ = 0.0;
    truncation_shape_right_ = 0.0;
    deposit_shape_left_ = 0.0;
    deposit_shape_right_ = 0.0;
    migration_failed_ = false;
}

void BackgroundTailPIC::accumulate_segment(
    double x0, double x1, double weight, double dt_seg,
    const SpatialGrid& sg, std::vector<double>& charge_face)
{
    if (!(dt_seg > 0.0) || x0 == x1) return;
    ChargeConservingTrajectory1D::deposit_charge_delta(
        x0, x1, weight, sg, charge_face);
}

void BackgroundTailPIC::add_removed_shape(
    const CellDepositWeights& cw, double weight, const SpatialGrid& sg,
    DriftOutputs& out)
{
    const int cells[2] = { cw.cell0, cw.cell1 };
    const double vals[2] = { cw.w0, cw.w1 };
    for (int k = 0; k < 2; ++k) {
        if (cells[k] < 0) {
            *out.trunc_left += weight * vals[k];
        } else if (cells[k] >= sg.nx_global) {
            *out.trunc_right += weight * vals[k];
        } else {
            const int il = cells[k] - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                (*out.trunc_in_domain)[static_cast<size_t>(il)] +=
                    weight * vals[k];
            }
        }
    }
}

void BackgroundTailPIC::process_drift_particle(
    const BackgroundTailParticle& p, double remaining_dt, double vx,
    double x_left, double x_right, double length, const SpatialGrid& sg,
    int mpi_rank, int mpi_size, DriftOutputs& out)
{
    double x0 = p.x;
    double t_rem = remaining_dt;
    while (t_rem > 0.0) {
        const double x1 = x0 + vx * t_rem;
        if (x1 < 0.0 || x1 >= length) {
            // Open-boundary truncation: deposit only up to the first
            // physical-boundary crossing (constraint 4).
            const double boundary = (x1 < 0.0) ? 0.0 : length;
            double t_cross = 0.0;
            if (vx != 0.0) t_cross = (boundary - x0) / vx;
            t_cross = std::max(0.0, std::min(t_rem, t_cross));
            const double x_cross = x0 + vx * t_cross;
            if (t_cross > 0.0) {
                accumulate_segment(x0, x_cross, p.weight, t_cross, sg,
                                   *out.charge_face);
            }
            const CellDepositWeights cw =
                ParticleShape1D::cell_weights(x_cross, sg);
            add_removed_shape(cw, p.weight, sg, out);
            const double ke = tail_kinetic_energy_per_particle(p);
            const double px = Const::me * Const::c * p.ux * p.weight;
            if (boundary == 0.0) {
                out.outflow->left_number += p.weight;
                out.outflow->left_px += px;
                out.outflow->left_kinetic_energy += ke;
            } else {
                out.outflow->right_number += p.weight;
                out.outflow->right_px += px;
                out.outflow->right_kinetic_energy += ke;
            }
            return;
        }
        if (x1 < x_left || x1 >= x_right) {
            // Internal MPI rank boundary: split the segment at the shared
            // edge and hand the particle plus the remaining drift time to
            // the owning neighbor (constraint 5).
            const double x_edge = (x1 < x_left) ? x_left : x_right;
            double t_edge = 0.0;
            if (vx != 0.0) t_edge = (x_edge - x0) / vx;
            t_edge = std::max(0.0, std::min(t_rem, t_edge));
            if (t_edge > 0.0) {
                accumulate_segment(x0, x_edge, p.weight, t_edge, sg,
                                   *out.charge_face);
            }
            TailDriftCarrier carrier;
            carrier.particle = p;
            carrier.particle.x = x_edge;
            carrier.remaining_dt = t_rem - t_edge;
            if (x_edge == x_left && mpi_rank > 0) {
                out.send_left->push_back(carrier);
            } else if (x_edge == x_right && mpi_rank + 1 < mpi_size) {
                out.send_right->push_back(carrier);
            } else {
                out.keep->push_back(carrier.particle);
            }
            return;
        }
        if (t_rem > 0.0) {
            accumulate_segment(x0, x1, p.weight, t_rem, sg,
                               *out.charge_face);
        }
        BackgroundTailParticle moved = p;
        moved.x = x1;
        out.keep->push_back(moved);
        return;
    }
    // Zero remaining drift: the particle stays where it is.
    out.keep->push_back(p);
}

void BackgroundTailPIC::exchange_carriers(const SpatialGrid& sg,
                                          int mpi_rank, int mpi_size)
{
    (void)sg;
    const std::chrono::steady_clock::time_point exchange_begin =
        std::chrono::steady_clock::now();
    recv_left_.clear();
    recv_right_.clear();
    if (mpi_size == 1) {
        send_left_.clear();
        send_right_.clear();
        last_migration_seconds_ = 0.0;
        return;
    }

    const int left = mpi_rank - 1;
    const int right = mpi_rank + 1;
    int send_left_count = static_cast<int>(send_left_.size());
    int send_right_count = static_cast<int>(send_right_.size());
    int recv_left_count = 0;
    int recv_right_count = 0;

    MPI_Request count_reqs[4];
    int nreq = 0;
    if (left >= 0) {
        MPI_Isend(&send_left_count, 1, MPI_INT, left, 601, MPI_COMM_WORLD,
                  &count_reqs[nreq++]);
        MPI_Irecv(&recv_left_count, 1, MPI_INT, left, 602, MPI_COMM_WORLD,
                  &count_reqs[nreq++]);
    }
    if (right < mpi_size) {
        MPI_Isend(&send_right_count, 1, MPI_INT, right, 602, MPI_COMM_WORLD,
                  &count_reqs[nreq++]);
        MPI_Irecv(&recv_right_count, 1, MPI_INT, right, 601, MPI_COMM_WORLD,
                  &count_reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, count_reqs, MPI_STATUSES_IGNORE);

    recv_left_.resize(static_cast<size_t>(recv_left_count));
    recv_right_.resize(static_cast<size_t>(recv_right_count));

    MPI_Request data_reqs[4];
    nreq = 0;
    if (left >= 0 && send_left_count > 0) {
        MPI_Isend(send_left_.data(), carrier_byte_count(send_left_count),
                  MPI_BYTE, left, 603, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (left >= 0 && recv_left_count > 0) {
        MPI_Irecv(recv_left_.data(), carrier_byte_count(recv_left_count),
                  MPI_BYTE, left, 604, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (right < mpi_size && send_right_count > 0) {
        MPI_Isend(send_right_.data(), carrier_byte_count(send_right_count),
                  MPI_BYTE, right, 604, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (right < mpi_size && recv_right_count > 0) {
        MPI_Irecv(recv_right_.data(), carrier_byte_count(recv_right_count),
                  MPI_BYTE, right, 603, MPI_COMM_WORLD, &data_reqs[nreq++]);
    }
    if (nreq > 0) MPI_Waitall(nreq, data_reqs, MPI_STATUSES_IGNORE);

    send_left_.clear();
    send_right_.clear();
    last_migration_seconds_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - exchange_begin).count();
}

void BackgroundTailPIC::export_accepted_state(
    BackgroundTailStateSnapshot& snapshot) const
{
    snapshot.particles = particles;
    snapshot.density = density;
    snapshot.id_counter = id_counter_;
    snapshot.outflow = outflow_;
    snapshot.truncation_shape_left = truncation_shape_left_;
    snapshot.truncation_shape_right = truncation_shape_right_;
    snapshot.deposit_shape_left = deposit_shape_left_;
    snapshot.deposit_shape_right = deposit_shape_right_;
    snapshot.deposit_shape_step_start_left = deposit_shape_step_start_left_;
    snapshot.deposit_shape_step_start_right = deposit_shape_step_start_right_;
    snapshot.collision_rng_seed = collision_rng_seed_;
    snapshot.max_abs_u = max_abs_u_;
    snapshot.max_kinetic_energy = max_kinetic_energy_;
}

void BackgroundTailPIC::import_accepted_state(
    const BackgroundTailStateSnapshot& snapshot)
{
    particles = snapshot.particles;
    density = snapshot.density;
    id_counter_ = snapshot.id_counter;
    outflow_ = snapshot.outflow;
    truncation_shape_left_ = snapshot.truncation_shape_left;
    truncation_shape_right_ = snapshot.truncation_shape_right;
    deposit_shape_left_ = snapshot.deposit_shape_left;
    deposit_shape_right_ = snapshot.deposit_shape_right;
    deposit_shape_step_start_left_ = snapshot.deposit_shape_step_start_left;
    deposit_shape_step_start_right_ = snapshot.deposit_shape_step_start_right;
    collision_rng_seed_ = snapshot.collision_rng_seed;
    max_abs_u_ = snapshot.max_abs_u;
    max_kinetic_energy_ = snapshot.max_kinetic_energy;
    migration_failed_ = false;
    last_migration_seconds_ = 0.0;
}

void BackgroundTailPIC::drift_half(const SpatialGrid& sg, double dt,
                                   int mpi_rank, int mpi_size)
{
    const double length = sg.dx * static_cast<double>(sg.nx_global);
    const double x_left = sg.ix_start * sg.dx;
    const double x_right = (sg.ix_start + sg.nx_local) * sg.dx;
    const int max_hops = std::max(1, mpi_size);

    keep_.clear();
    send_left_.clear();
    send_right_.clear();
    recv_left_.clear();
    recv_right_.clear();

    const long long np = static_cast<long long>(particles.size());
    const int nthreads = std::max(1, omp_get_max_threads());
    if (thread_keep_.size() != static_cast<size_t>(nthreads)) {
        thread_keep_.resize(nthreads);
        thread_send_left_.resize(nthreads);
        thread_send_right_.resize(nthreads);
        thread_charge_.resize(nthreads);
        thread_trunc_in_domain_.resize(nthreads);
        thread_trunc_left_.resize(nthreads);
        thread_trunc_right_.resize(nthreads);
        thread_outflow_.resize(nthreads);
    }
    const size_t charge_size = static_cast<size_t>(sg.nx_local) + 1;
    const size_t trunc_size = static_cast<size_t>(sg.nx_local);
    for (int t = 0; t < nthreads; ++t) {
        thread_keep_[t].clear();
        thread_send_left_[t].clear();
        thread_send_right_[t].clear();
        if (thread_charge_[t].size() != charge_size)
            thread_charge_[t].assign(charge_size, 0.0);
        else
            std::fill(thread_charge_[t].begin(), thread_charge_[t].end(),
                      0.0);
        if (thread_trunc_in_domain_[t].size() != trunc_size)
            thread_trunc_in_domain_[t].assign(trunc_size, 0.0);
        else
            std::fill(thread_trunc_in_domain_[t].begin(),
                      thread_trunc_in_domain_[t].end(), 0.0);
        thread_trunc_left_[t] = 0.0;
        thread_trunc_right_[t] = 0.0;
        thread_outflow_[t] = TailOutflowLedger();
    }

    // Initial pass over the local particles.
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        DriftOutputs out;
        out.keep = &thread_keep_[tid];
        out.send_left = &thread_send_left_[tid];
        out.send_right = &thread_send_right_[tid];
        out.charge_face = &thread_charge_[tid];
        out.trunc_in_domain = &thread_trunc_in_domain_[tid];
        out.trunc_left = &thread_trunc_left_[tid];
        out.trunc_right = &thread_trunc_right_[tid];
        out.outflow = &thread_outflow_[tid];
        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            const BackgroundTailParticle& p =
                particles[static_cast<size_t>(i)];
            process_drift_particle(p, dt, tail_vx(p), x_left, x_right,
                                   length, sg, mpi_rank, mpi_size, out);
        }
    }

    particles.clear();
    for (int t = 0; t < nthreads; ++t) {
        keep_.insert(keep_.end(), thread_keep_[t].begin(),
                     thread_keep_[t].end());
        send_left_.insert(send_left_.end(), thread_send_left_[t].begin(),
                          thread_send_left_[t].end());
        send_right_.insert(send_right_.end(), thread_send_right_[t].begin(),
                           thread_send_right_[t].end());
        for (size_t i = 0; i < charge_size; ++i) {
            trajectory_charge_face_x_[i] += thread_charge_[t][i];
        }
        for (size_t i = 0; i < trunc_size; ++i) {
            truncated_in_domain_[i] += thread_trunc_in_domain_[t][i];
        }
        truncation_shape_left_ += thread_trunc_left_[t];
        truncation_shape_right_ += thread_trunc_right_[t];
        outflow_.left_number += thread_outflow_[t].left_number;
        outflow_.left_px += thread_outflow_[t].left_px;
        outflow_.left_kinetic_energy +=
            thread_outflow_[t].left_kinetic_energy;
        outflow_.right_number += thread_outflow_[t].right_number;
        outflow_.right_px += thread_outflow_[t].right_px;
        outflow_.right_kinetic_energy +=
            thread_outflow_[t].right_kinetic_energy;
    }

    // Continuation loop: received carriers finish their remaining drift on
    // the owning rank.  With the CFL contract satisfied one round suffices;
    // the bounded loop is a deterministic multi-hop safety net (section
    // 6.6.5) and trips migration_failed_ before any state is accepted.
    int hop = 0;
    while (true) {
        exchange_carriers(sg, mpi_rank, mpi_size);
        const int any_recv =
            (!recv_left_.empty() || !recv_right_.empty()) ? 1 : 0;
        int global_any = 0;
        MPI_Allreduce(&any_recv, &global_any, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        if (global_any == 0) break;
        if (hop >= max_hops) {
            migration_failed_ = true;
            break;
        }
        ++hop;
        DriftOutputs out;
        out.keep = &keep_;
        out.send_left = &send_left_;
        out.send_right = &send_right_;
        out.charge_face = &trajectory_charge_face_x_;
        out.trunc_in_domain = &truncated_in_domain_;
        out.trunc_left = &truncation_shape_left_;
        out.trunc_right = &truncation_shape_right_;
        out.outflow = &outflow_;
        for (size_t i = 0; i < recv_left_.size(); ++i) {
            const TailDriftCarrier& c = recv_left_[i];
            process_drift_particle(c.particle, c.remaining_dt, tail_vx(c.particle),
                                   x_left, x_right, length, sg, mpi_rank,
                                   mpi_size, out);
        }
        for (size_t i = 0; i < recv_right_.size(); ++i) {
            const TailDriftCarrier& c = recv_right_[i];
            process_drift_particle(c.particle, c.remaining_dt, tail_vx(c.particle),
                                   x_left, x_right, length, sg, mpi_rank,
                                   mpi_size, out);
        }
        recv_left_.clear();
        recv_right_.clear();
    }

    particles.swap(keep_);
    keep_.clear();
    refresh_moment_stats();
}

void BackgroundTailPIC::kick(const SpatialGrid& sg,
                             const EMFields& field_mid, double dt,
                             int mpi_rank, int mpi_size,
                             double* local_kinetic_work,
                             std::vector<double>* local_delta_ke_by_x,
                             double* local_delta_ke_boundary)
{
    (void)mpi_rank;
    (void)mpi_size;
    if (local_delta_ke_by_x != NULL) {
        local_delta_ke_by_x->assign(static_cast<size_t>(sg.nx_local), 0.0);
    }
    const long long np = static_cast<long long>(particles.size());
    const double kick_coeff = Const::qe / (Const::me * Const::c);
    const int nx_local = sg.nx_local;
    double work_sum = 0.0;
    double boundary_work = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:work_sum) \
        reduction(+:boundary_work)
    for (long long i = 0; i < np; ++i) {
        BackgroundTailParticle& p = particles[static_cast<size_t>(i)];
        const FaceGatherWeights gw =
            ParticleShape1D::shifted_face_weights(p.x, sg);
        int il0 = gw.face0 - sg.ix_start;
        int il1 = gw.face1 - sg.ix_start;
        il0 = std::max(0, std::min(nx_local, il0));
        il1 = std::max(0, std::min(nx_local, il1));
        const double ex =
            gw.w0 * field_mid.Ex_face[static_cast<size_t>(il0)] +
            gw.w1 * field_mid.Ex_face[static_cast<size_t>(il1)];
        if (local_kinetic_work != NULL) {
            // Same relativistic kinetic energy as tail_total_kinetic_energy.
            const double ke_before = tail_kinetic_energy_per_particle(p);
            p.ux -= kick_coeff * ex * dt;
            const double ke_after = tail_kinetic_energy_per_particle(p);
            const double delta = p.weight * (ke_after - ke_before);
            work_sum += delta;
            if (local_delta_ke_by_x != NULL) {
                // Gate I (section 4.4): CIC cell attribution at the kick
                // position x^{n+1/2}.  Out-of-domain shares go to the boundary
                // work audit (section I3), never renormalized into the domain.
                const CellDepositWeights cw =
                    ParticleShape1D::cell_weights(p.x, sg);
                const int c0 = cw.cell0 - sg.ix_start;
                const int c1 = cw.cell1 - sg.ix_start;
                if (c0 >= 0 && c0 < nx_local) {
                    #pragma omp atomic
                    (*local_delta_ke_by_x)[static_cast<size_t>(c0)] +=
                        delta * cw.w0;
                } else {
                    boundary_work += delta * cw.w0;
                }
                if (c1 >= 0 && c1 < nx_local) {
                    #pragma omp atomic
                    (*local_delta_ke_by_x)[static_cast<size_t>(c1)] +=
                        delta * cw.w1;
                } else {
                    boundary_work += delta * cw.w1;
                }
            }
        } else {
            p.ux -= kick_coeff * ex * dt;
        }
    }
    if (local_kinetic_work != NULL) {
        *local_kinetic_work = work_sum;
    }
    if (local_delta_ke_boundary != NULL) {
        *local_delta_ke_boundary = boundary_work;
    }
    refresh_moment_stats();
}

void BackgroundTailPIC::deposit_density(const SpatialGrid& sg,
                                        int mpi_rank, int mpi_size)
{
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    std::fill(density.begin(), density.end(), 0.0);
    deposit_shape_left_ = 0.0;
    deposit_shape_right_ = 0.0;

    const int nthreads = std::max(1, omp_get_max_threads());
    std::vector<std::vector<double> > thread_density(nthreads);
    std::vector<double> thread_send_left(nthreads, 0.0);
    std::vector<double> thread_send_right(nthreads, 0.0);
    std::vector<double> thread_ledger_left(nthreads, 0.0);
    std::vector<double> thread_ledger_right(nthreads, 0.0);
    for (int t = 0; t < nthreads; ++t) {
        thread_density[t].assign(nxl, 0.0);
    }

    const double inv_dx = 1.0 / sg.dx;
    const long long np = static_cast<long long>(particles.size());
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        std::vector<double>& local = thread_density[tid];
        double local_send_left = 0.0;
        double local_send_right = 0.0;
        double local_ledger_left = 0.0;
        double local_ledger_right = 0.0;
        #pragma omp for schedule(static)
        for (long long i = 0; i < np; ++i) {
            const BackgroundTailParticle& p =
                particles[static_cast<size_t>(i)];
            const CellDepositWeights cw =
                ParticleShape1D::cell_weights(p.x, sg);
            const int cells[2] = { cw.cell0, cw.cell1 };
            const double vals[2] = {
                p.weight * cw.w0 * inv_dx,
                p.weight * cw.w1 * inv_dx
            };
            for (int k = 0; k < 2; ++k) {
                const int ig = cells[k];
                if (ig < 0) {
                    local_ledger_left += p.weight * ((k == 0) ? cw.w0 : cw.w1);
                } else if (ig >= sg.nx_global) {
                    local_ledger_right += p.weight * ((k == 0) ? cw.w0 : cw.w1);
                } else if (ig < sg.ix_start) {
                    local_send_left += vals[k];
                } else if (ig >= sg.ix_start + sg.nx_local) {
                    local_send_right += vals[k];
                } else {
                    local[static_cast<size_t>(ig - sg.ix_start)] += vals[k];
                }
            }
        }
        thread_send_left[tid] = local_send_left;
        thread_send_right[tid] = local_send_right;
        thread_ledger_left[tid] = local_ledger_left;
        thread_ledger_right[tid] = local_ledger_right;
    }

    for (int t = 0; t < nthreads; ++t) {
        for (size_t i = 0; i < nxl; ++i) {
            density[i] += thread_density[static_cast<size_t>(t)][i];
        }
        deposit_shape_left_ += thread_ledger_left[t];
        deposit_shape_right_ += thread_ledger_right[t];
    }

    double send_left_density = 0.0;
    double send_right_density = 0.0;
    for (int t = 0; t < nthreads; ++t) {
        send_left_density += thread_send_left[t];
        send_right_density += thread_send_right[t];
    }
    if (mpi_size > 1) {
        double recv_left_density = 0.0;
        double recv_right_density = 0.0;
        MPI_Request reqs[4];
        int nreq = 0;
        if (mpi_rank > 0) {
            MPI_Isend(&send_left_density, 1, MPI_DOUBLE, mpi_rank - 1, 701,
                      MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Irecv(&recv_left_density, 1, MPI_DOUBLE, mpi_rank - 1, 702,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (mpi_rank + 1 < mpi_size) {
            MPI_Isend(&send_right_density, 1, MPI_DOUBLE, mpi_rank + 1, 702,
                      MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Irecv(&recv_right_density, 1, MPI_DOUBLE, mpi_rank + 1, 701,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
        if (mpi_rank > 0 && !density.empty()) density[0] += recv_left_density;
        if (mpi_rank + 1 < mpi_size && !density.empty())
            density[static_cast<size_t>(sg.nx_local - 1)] +=
                recv_right_density;
    }
}

void BackgroundTailPIC::finalize_trajectory_current(
    const SpatialGrid& sg, double dt, int mpi_rank, int mpi_size)
{
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    if (!(dt > 0.0)) {
        if (current_face_x.size() != nxl + 1)
            current_face_x.assign(nxl + 1, 0.0);
        else
            std::fill(current_face_x.begin(), current_face_x.end(), 0.0);
        return;
    }
    const double inv_dt = 1.0 / dt;
    if (current_face_x.size() != nxl + 1)
        current_face_x.assign(nxl + 1, 0.0);
    for (size_t i = 0; i <= nxl; ++i) {
        current_face_x[i] = trajectory_charge_face_x_[i] * inv_dt;
    }
    if (mpi_size > 1) {
        // Close the shared faces: the flux across the rank boundary is the
        // sum of the contributions deposited by both ranks' segments.
        const double send_left_face = current_face_x[0];
        const double send_right_face = current_face_x[nxl];
        double recv_left_face = 0.0;
        double recv_right_face = 0.0;
        MPI_Request reqs[4];
        int nreq = 0;
        if (mpi_rank > 0) {
            MPI_Isend(&send_left_face, 1, MPI_DOUBLE, mpi_rank - 1, 711,
                      MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Irecv(&recv_left_face, 1, MPI_DOUBLE, mpi_rank - 1, 712,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (mpi_rank + 1 < mpi_size) {
            MPI_Isend(&send_right_face, 1, MPI_DOUBLE, mpi_rank + 1, 712,
                      MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Irecv(&recv_right_face, 1, MPI_DOUBLE, mpi_rank + 1, 711,
                      MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
        if (mpi_rank > 0) current_face_x[0] += recv_left_face;
        if (mpi_rank + 1 < mpi_size) current_face_x[nxl] += recv_right_face;
    }
}

bool BackgroundTailPIC::cfl_contract_ok(const SpatialGrid& sg, double dt,
                                        int mpi_size) const
{
    (void)mpi_size;
    const double rank_width = sg.dx * static_cast<double>(sg.nx_local);
    return dt > 0.0 && Const::c * dt < rank_width;
}

void BackgroundTailPIC::refresh_moment_stats()
{
    max_abs_u_ = 0.0;
    max_kinetic_energy_ = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        const BackgroundTailParticle& p = particles[i];
        const double abs_u =
            std::sqrt(p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
        max_abs_u_ = std::max(max_abs_u_, abs_u);
        max_kinetic_energy_ = std::max(
            max_kinetic_energy_, tail_kinetic_energy_per_particle(p));
    }
}

void BackgroundTailPIC::swap_state(BackgroundTailPIC& other)
{
    particles.swap(other.particles);
    density.swap(other.density);
    current_face_x.swap(other.current_face_x);
    density_step_start.swap(other.density_step_start);
    std::swap(step_dt_, other.step_dt_);
    std::swap(outflow_, other.outflow_);
    std::swap(truncation_shape_left_, other.truncation_shape_left_);
    std::swap(truncation_shape_right_, other.truncation_shape_right_);
    std::swap(deposit_shape_left_, other.deposit_shape_left_);
    std::swap(deposit_shape_right_, other.deposit_shape_right_);
    std::swap(deposit_shape_step_start_left_,
              other.deposit_shape_step_start_left_);
    std::swap(deposit_shape_step_start_right_,
              other.deposit_shape_step_start_right_);
    std::swap(max_abs_u_, other.max_abs_u_);
    std::swap(max_kinetic_energy_, other.max_kinetic_energy_);
    std::swap(id_counter_, other.id_counter_);
    std::swap(collision_rng_seed_, other.collision_rng_seed_);
    std::swap(migration_failed_, other.migration_failed_);
    trajectory_charge_face_x_.swap(other.trajectory_charge_face_x_);
    truncated_in_domain_.swap(other.truncated_in_domain_);
    keep_.swap(other.keep_);
    send_left_.swap(other.send_left_);
    send_right_.swap(other.send_right_);
    recv_left_.swap(other.recv_left_);
    recv_right_.swap(other.recv_right_);
    thread_keep_.swap(other.thread_keep_);
    thread_send_left_.swap(other.thread_send_left_);
    thread_send_right_.swap(other.thread_send_right_);
    thread_charge_.swap(other.thread_charge_);
    thread_trunc_in_domain_.swap(other.thread_trunc_in_domain_);
    thread_trunc_left_.swap(other.thread_trunc_left_);
    thread_trunc_right_.swap(other.thread_trunc_right_);
    thread_outflow_.swap(other.thread_outflow_);
}

bool BackgroundTailPIC::finite() const
{
    for (size_t i = 0; i < particles.size(); ++i) {
        const BackgroundTailParticle& p = particles[i];
        if (!std::isfinite(p.x) || !std::isfinite(p.ux) ||
            !std::isfinite(p.uy) || !std::isfinite(p.uz) ||
            !std::isfinite(p.weight)) {
            return false;
        }
    }
    for (size_t i = 0; i < density.size(); ++i)
        if (!std::isfinite(density[i])) return false;
    for (size_t i = 0; i < current_face_x.size(); ++i)
        if (!std::isfinite(current_face_x[i])) return false;
    for (size_t i = 0; i < trajectory_charge_face_x_.size(); ++i)
        if (!std::isfinite(trajectory_charge_face_x_[i])) return false;
    for (size_t i = 0; i < truncated_in_domain_.size(); ++i)
        if (!std::isfinite(truncated_in_domain_[i])) return false;
    if (!std::isfinite(outflow_.left_number) ||
        !std::isfinite(outflow_.right_number) ||
        !std::isfinite(truncation_shape_left_) ||
        !std::isfinite(truncation_shape_right_) ||
        !std::isfinite(deposit_shape_left_) ||
        !std::isfinite(deposit_shape_right_) ||
        !std::isfinite(deposit_shape_step_start_left_) ||
        !std::isfinite(deposit_shape_step_start_right_)) {
        return false;
    }
    return true;
}

bool BackgroundTailPIC::nonnegative_weights() const
{
    for (size_t i = 0; i < particles.size(); ++i) {
        if (!(particles[i].weight > 0.0)) return false;
    }
    return true;
}

double BackgroundTailPIC::density_weight_balance(
    const SpatialGrid& sg) const
{
    double in_domain = 0.0;
    for (size_t i = 0; i < density.size(); ++i) {
        in_domain += density[i] * sg.dx;
    }
    double total_weight = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        total_weight += particles[i].weight;
    }
    return std::fabs(in_domain + deposit_shape_left_ + deposit_shape_right_ -
                     total_weight);
}

TailContinuityResult BackgroundTailPIC::audit_continuity(
    const SpatialGrid& sg, double dt,
    const std::vector<double>& conversion_source,
    int mpi_rank, int mpi_size) const
{
    (void)mpi_rank;
    (void)mpi_size;
    TailContinuityResult r;
    const size_t nxl = static_cast<size_t>(sg.nx_local);
    const double dx = sg.dx;
    const double inv_qe = 1.0 / Const::qe;
    const bool has_source = (conversion_source.size() == nxl);

    double n0_total = deposit_shape_step_start_left_ +
                      deposit_shape_step_start_right_;
    double n1_total = deposit_shape_left_ + deposit_shape_right_;
    double c_total = 0.0;
    double abs_l1 = 0.0;
    double abs_linf = 0.0;
    double max_cell_scale = 1.0;
    for (size_t i = 0; i < nxl; ++i) {
        const double n0 = density_step_start[i];
        const double n1 = density[i];
        // Number-flux divergence from the trajectory current: F = -J/q_e.
        const double flux =
            -(current_face_x[i + 1] - current_face_x[i]) * inv_qe * dt;
        const double c = has_source ? conversion_source[i] * dx * dt : 0.0;
        const double removed = truncated_in_domain_[i];
        const double residual = (n1 - n0) * dx + flux - c + removed;
        n0_total += n0 * dx;
        n1_total += n1 * dx;
        c_total += c;
        abs_l1 += std::fabs(residual) * dx;
        abs_linf = std::max(abs_linf, std::fabs(residual));
        const double cell_scale =
            (n0 + n1) * dx + c + std::fabs(flux) + removed;
        max_cell_scale = std::max(max_cell_scale, cell_scale);
    }
    const double outflow_total =
        outflow_.left_number + outflow_.right_number;
    const double local_scale = std::max(
        1.0, n0_total + n1_total + c_total + outflow_total);
    r.abs_l1 = abs_l1;
    r.abs_linf = abs_linf;
    r.rel_l1 = abs_l1 / local_scale;
    r.rel_linf = abs_linf / max_cell_scale;
    // The number balance is a combined-system global invariant (section
    // 11.1): when a particle crosses an internal MPI rank boundary its
    // number leaves one rank's local density and enters the neighbor's, so
    // only the MPI sum of the signed local balances cancels.  The physical
    // outflow terms are already included (they only appear on the ranks
    // owning the physical boundary).
    const double local_balance =
        n1_total - n0_total - c_total + outflow_total;
    double global_balance = 0.0;
    MPI_Allreduce(&local_balance, &global_balance, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double global_scale = 0.0;
    MPI_Allreduce(&local_scale, &global_scale, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    r.number_balance_error =
        std::fabs(global_balance) / std::max(1.0, global_scale);

    double trunc_sum = 0.0;
    for (size_t i = 0; i < nxl; ++i) trunc_sum += truncated_in_domain_[i];
    const double trunc_total =
        trunc_sum + truncation_shape_left_ + truncation_shape_right_;
    r.truncation_ledger_error =
        std::fabs(trunc_total - outflow_total) /
        std::max(1.0, outflow_total);
    return r;
}

BackgroundTailPIC::WeightStats BackgroundTailPIC::weight_stats(
    const SpatialGrid& sg,
    const std::vector<double>* background_density) const
{
    WeightStats s;
    s.macro_particle_count = particles.size();
    if (particles.empty()) return s;

    double sum = 0.0;
    double sum2 = 0.0;
    s.min = particles[0].weight;
    s.max = particles[0].weight;
    std::vector<size_t> per_cell(static_cast<size_t>(sg.nx_local), 0);
    std::vector<double> cell_number(static_cast<size_t>(sg.nx_local), 0.0);
    for (size_t i = 0; i < particles.size(); ++i) {
        const BackgroundTailParticle& p = particles[i];
        sum += p.weight;
        sum2 += p.weight * p.weight;
        s.min = std::min(s.min, p.weight);
        s.max = std::max(s.max, p.weight);
        const CellDepositWeights cw =
            ParticleShape1D::cell_weights(p.x, sg);
        for (int k = 0; k < 2; ++k) {
            const int il = cw.cell0 + k - sg.ix_start;
            if (il >= 0 && il < sg.nx_local) {
                per_cell[static_cast<size_t>(il)] += 1;
                cell_number[static_cast<size_t>(il)] +=
                    p.weight * ((k == 0) ? cw.w0 : cw.w1);
            }
        }
    }
    s.mean = sum / static_cast<double>(particles.size());
    const double variance = std::max(
        0.0, sum2 / static_cast<double>(particles.size()) - s.mean * s.mean);
    s.std = std::sqrt(variance);
    for (size_t i = 0; i < per_cell.size(); ++i) {
        s.per_cell_max_macro_count =
            std::max(s.per_cell_max_macro_count, per_cell[i]);
    }

    double max_fraction = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        const BackgroundTailParticle& p = particles[i];
        const CellDepositWeights cw =
            ParticleShape1D::cell_weights(p.x, sg);
        int il = cw.cell0 - sg.ix_start;
        if (il < 0 || il >= sg.nx_local) il = cw.cell1 - sg.ix_start;
        if (il < 0 || il >= sg.nx_local) continue;
        const size_t slot = static_cast<size_t>(il);
        double cell_total = cell_number[slot];
        if (background_density != NULL &&
            slot < background_density->size()) {
            cell_total += (*background_density)[slot] * sg.dx;
        }
        const double denom = std::max(cell_total, p.weight * 1.0e-12);
        max_fraction = std::max(max_fraction, p.weight / denom);
    }
    s.max_single_charge_fraction = max_fraction;

    // Density noise estimate: relative L2 variation over occupied cells.
    double occupied_sum = 0.0;
    size_t occupied = 0;
    for (size_t i = 0; i < density.size(); ++i) {
        if (density[i] > 0.0) {
            occupied_sum += density[i];
            ++occupied;
        }
    }
    if (occupied > 0) {
        const double mean_density = occupied_sum /
            static_cast<double>(occupied);
        double sq = 0.0;
        for (size_t i = 0; i < density.size(); ++i) {
            if (density[i] > 0.0) {
                const double d = density[i] - mean_density;
                sq += d * d;
            }
        }
        s.density_noise_estimate =
            std::sqrt(sq / static_cast<double>(occupied)) / mean_density;
    }
    return s;
}
