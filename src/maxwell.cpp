#include "maxwell.h"
#include "species.h"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <mpi.h>
#include <vector>

namespace {
void resize_or_zero(std::vector<double>& values, size_t n)
{
    if (values.size() != n) {
        values.assign(n, 0.0);
    } else {
        std::fill(values.begin(), values.end(), 0.0);
    }
}

void prepare_mpi_layout(EMFields& fields, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;

    if (fields.counts_mpi_size == mpi_size &&
        fields.counts_nx_local == nxl &&
        fields.global_rhs.size() == static_cast<size_t>(Param::nx)) {
        return;
    }

    fields.counts_mpi_size = mpi_size;
    fields.counts_nx_local = nxl;
    fields.counts.assign(static_cast<size_t>(mpi_size), 0);
    fields.displs.assign(static_cast<size_t>(mpi_size), 0);
    fields.face_counts.assign(static_cast<size_t>(mpi_size), 0);
    fields.face_displs.assign(static_cast<size_t>(mpi_size), 0);
    MPI_Allgather(&nxl, 1, MPI_INT, fields.counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int r = 0; r < mpi_size; ++r) {
        fields.face_counts[static_cast<size_t>(r)] =
            fields.counts[static_cast<size_t>(r)] + 1;
    }
    for (int r = 1; r < mpi_size; ++r) {
        fields.displs[static_cast<size_t>(r)] =
            fields.displs[static_cast<size_t>(r - 1)] +
            fields.counts[static_cast<size_t>(r - 1)];
        fields.face_displs[static_cast<size_t>(r)] =
            fields.displs[static_cast<size_t>(r)];
    }

    fields.local_rhs.assign(static_cast<size_t>(nxl), 0.0);
    fields.local_face_rhs.assign(static_cast<size_t>(nxl + 1), 0.0);
    fields.global_rhs.assign(static_cast<size_t>(Param::nx), 0.0);
    fields.global_ex.assign(static_cast<size_t>(Param::nx), 0.0);
    fields.global_face.assign(static_cast<size_t>(Param::nx + 1), 0.0);
    fields.global_phi.assign(static_cast<size_t>(Param::nx), 0.0);
}

void gather_rho(EMFields& fields, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;
    prepare_mpi_layout(fields, mpi_size);

    for (int ix = 0; ix < nxl; ++ix) {
        fields.local_rhs[static_cast<size_t>(ix)] = fields.rho[ng + ix];
    }

    MPI_Gatherv(fields.local_rhs.data(), nxl, MPI_DOUBLE,
                fields.global_rhs.data(), fields.counts.data(),
                fields.displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void compute_gauss_field(EMFields& fields, bool compute_phi)
{
    const int n = Param::nx;
    if (n <= 0) return;

    double mean_rho = 0.0;
    for (int i = 0; i < n; ++i) {
        mean_rho += fields.global_rhs[static_cast<size_t>(i)];
    }
    mean_rho /= static_cast<double>(n);

    resize_or_zero(fields.global_ex, static_cast<size_t>(n));
    resize_or_zero(fields.all_interfaces, static_cast<size_t>(n + 1));

    const double scale = fields.dx / Const::eps0;
    for (int i = 0; i < n; ++i) {
        fields.all_interfaces[static_cast<size_t>(i + 1)] =
            fields.all_interfaces[static_cast<size_t>(i)] +
            (fields.global_rhs[static_cast<size_t>(i)] - mean_rho) * scale;
    }

    double mean_ex = 0.0;
    for (int i = 0; i < n; ++i) {
        fields.global_ex[static_cast<size_t>(i)] =
            0.5 * (fields.all_interfaces[static_cast<size_t>(i)] +
                   fields.all_interfaces[static_cast<size_t>(i + 1)]);
        mean_ex += fields.global_ex[static_cast<size_t>(i)];
    }
    mean_ex /= static_cast<double>(n);
    for (int i = 0; i < n; ++i) {
        fields.global_ex[static_cast<size_t>(i)] -= mean_ex;
    }
    for (int i = 0; i <= n; ++i) {
        fields.all_interfaces[static_cast<size_t>(i)] -= mean_ex;
    }

    if (compute_phi) {
        resize_or_zero(fields.global_phi, static_cast<size_t>(n));
        for (int i = 1; i < n; ++i) {
            fields.global_phi[static_cast<size_t>(i)] =
                fields.global_phi[static_cast<size_t>(i - 1)] -
                0.5 * (fields.global_ex[static_cast<size_t>(i - 1)] +
                       fields.global_ex[static_cast<size_t>(i)]) * fields.dx;
        }
        double mean_phi = 0.0;
        for (int i = 0; i < n; ++i) {
            mean_phi += fields.global_phi[static_cast<size_t>(i)];
        }
        mean_phi /= static_cast<double>(n);
        for (int i = 0; i < n; ++i) {
            fields.global_phi[static_cast<size_t>(i)] -= mean_phi;
        }
    }
}

void exchange_scalar_ghosts(EMFields& fields,
                            std::vector<double>& a,
                            int tag_base,
                            int mpi_rank,
                            int mpi_size);

void close_periodic_right_face(std::vector<double>& face,
                               int nxl,
                               int mpi_rank,
                               int mpi_size,
                               int tag)
{
    if (nxl <= 0 || face.size() < static_cast<size_t>(nxl + 1)) return;
    if (mpi_size == 1) {
        face[static_cast<size_t>(nxl)] = face[0];
        return;
    }

    const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right_peer = (mpi_rank + 1) % mpi_size;
    const double send_left_face = face[0];
    double recv_right_face = 0.0;
    MPI_Request reqs[2];
    MPI_Isend(&send_left_face, 1, MPI_DOUBLE, left_peer, tag,
              MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(&recv_right_face, 1, MPI_DOUBLE, right_peer, tag,
              MPI_COMM_WORLD, &reqs[1]);
    MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
    face[static_cast<size_t>(nxl)] = recv_right_face;
}

void update_cell_ex_from_faces(EMFields& fields,
                               int mpi_rank,
                               int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;
    if (fields.Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        fields.Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }
    for (int ix = 0; ix < nxl; ++ix) {
        fields.Ex[ng + ix] =
            0.5 * (fields.Ex_face[static_cast<size_t>(ix)]
                 + fields.Ex_face[static_cast<size_t>(ix + 1)]);
    }
    exchange_scalar_ghosts(fields, fields.Ex, 201, mpi_rank, mpi_size);
}

void exchange_scalar_ghosts(EMFields& fields,
                            std::vector<double>& a,
                            int tag_base,
                            int mpi_rank,
                            int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = fields.nx_total - 2 * ng;

    if (fields.send_left.size() != static_cast<size_t>(ng)) {
        fields.send_left.resize(ng);
        fields.send_right.resize(ng);
        fields.recv_left.resize(ng);
        fields.recv_right.resize(ng);
    }

    if (mpi_size == 1) {
        for (int g = 0; g < ng; ++g) {
            const int left_src = ng + ((nxl - ng + g) % nxl);
            const int right_src = ng + (g % nxl);
            a[g] = a[left_src];
            a[ng + nxl + g] = a[right_src];
        }
        return;
    }

    for (int g = 0; g < ng; ++g) {
        fields.send_left[g] = a[ng + g];
        fields.send_right[g] = a[ng + nxl - ng + g];
    }

    MPI_Request reqs[4];
    int nreq = 0;
    const int left_peer = (mpi_rank + mpi_size - 1) % mpi_size;
    const int right_peer = (mpi_rank + 1) % mpi_size;

    MPI_Isend(fields.send_left.data(), ng, MPI_DOUBLE, left_peer, tag_base,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(fields.recv_left.data(), ng, MPI_DOUBLE, left_peer, tag_base + 1,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Isend(fields.send_right.data(), ng, MPI_DOUBLE, right_peer, tag_base + 1,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(fields.recv_right.data(), ng, MPI_DOUBLE, right_peer, tag_base,
              MPI_COMM_WORLD, &reqs[nreq++]);
    if (nreq > 0) MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    for (int g = 0; g < ng; ++g) {
        a[g] = fields.recv_left[g];
        a[ng + nxl + g] = fields.recv_right[g];
    }
}
}

void EMFields::init(const SpatialGrid& sg)
{
    nx_total = sg.nx_total;
    counts_mpi_size = 0;
    counts_nx_local = -1;
    dx = sg.dx;
    Ex.assign(nx_total, 0.0);
    Ex_face.assign(sg.nx_local + 1, 0.0);
    phi.assign(nx_total, 0.0);
    rho.assign(nx_total, 0.0);
    send_left.assign(sg.nghost, 0.0);
    send_right.assign(sg.nghost, 0.0);
    recv_left.assign(sg.nghost, 0.0);
    recv_right.assign(sg.nghost, 0.0);
    last_gauss_residual_l1 = 0.0;
    last_gauss_residual_linf = 0.0;
}

void EMFields::zero_currents()
{
    std::fill(rho.begin(), rho.end(), 0.0);
}

void EMFields::accumulate_moments(const Species& sp)
{
    const int ng = sp.sgrid->nghost;
    const int nxl = sp.sgrid->nx_local;
    for (int ix = 0; ix < nxl; ++ix) {
        rho[ix + ng] += sp.charge_density[ix];
    }
}

void EMFields::set_charge_density(const Species& electrons,
                                  const std::vector<double>& beam_density,
                                  const std::vector<double>& ion_density_profile)
{
    const int ng = electrons.sgrid->nghost;
    const int nxl = electrons.sgrid->nx_local;
    std::fill(rho.begin(), rho.end(), 0.0);
    const bool inputs_sized =
        ion_density_profile.size() >= static_cast<size_t>(nxl) &&
        beam_density.size() >= static_cast<size_t>(nxl);
    if (inputs_sized) {
        for (int ix = 0; ix < nxl; ++ix) {
            const size_t slot = static_cast<size_t>(ix);
            rho[ix + ng] =
                Const::qe * (ion_density_profile[slot]
                           - electrons.number_density[slot]
                           - beam_density[slot]);
        }
    } else {
        for (int ix = 0; ix < nxl; ++ix) {
            const size_t slot = static_cast<size_t>(ix);
            const double zni = (slot < ion_density_profile.size())
                             ? ion_density_profile[slot] : 0.0;
            const double nb = (slot < beam_density.size())
                            ? beam_density[slot] : 0.0;
            rho[ix + ng] =
                Const::qe * (zni - electrons.number_density[slot] - nb);
        }
    }
}

void EMFields::solve_poisson(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    gather_rho(*this, mpi_size);

    if (mpi_rank == 0) {
        compute_gauss_field(*this, false);
    }

    MPI_Scatterv(all_interfaces.data(), face_counts.data(),
                 face_displs.data(), MPI_DOUBLE,
                 local_face_rhs.data(), nxl + 1, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    if (Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }
    for (int ix = 0; ix <= nxl; ++ix) {
        Ex_face[static_cast<size_t>(ix)] =
            local_face_rhs[static_cast<size_t>(ix)];
    }
    close_periodic_right_face(Ex_face, nxl, mpi_rank, mpi_size, 301);
    update_cell_ex_from_faces(*this, mpi_rank, mpi_size);
}

void EMFields::advance_ampere_face(const std::vector<double>& background_current_face,
                                   const std::vector<double>& open_beam_current_face,
                                   double dt,
                                   int mpi_rank,
                                   int mpi_size)
{
    const int nxl = nx_total - 2 * Param::Nghost;
    if (dt <= 0.0 || nxl <= 0) return;
    if (Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }
    const bool currents_sized =
        background_current_face.size() >= static_cast<size_t>(nxl) &&
        open_beam_current_face.size() >= static_cast<size_t>(nxl);
    const double ampere_scale = -dt / Const::eps0;

    if (currents_sized) {
        // Only unique periodic field faces are advanced. Beam guard/outside
        // currents have already been cleared and are never folded or wrapped.
        for (int iface = 0; iface < nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            Ex_face[slot] +=
                ampere_scale * (background_current_face[slot]
                              + open_beam_current_face[slot]);
        }
    } else {
        for (int iface = 0; iface < nxl; ++iface) {
            const size_t slot = static_cast<size_t>(iface);
            const double jb = (slot < background_current_face.size())
                            ? background_current_face[slot] : 0.0;
            const double jbeam = (slot < open_beam_current_face.size())
                                ? open_beam_current_face[slot] : 0.0;
            Ex_face[slot] +=
                ampere_scale * (jb + jbeam);
        }
    }

    // Periodicity is applied only to Ex after the local Ampere update.
    close_periodic_right_face(Ex_face, nxl, mpi_rank, mpi_size, 302);
    update_cell_ex_from_faces(*this, mpi_rank, mpi_size);
}

void EMFields::advance_ampere_face_from_midpoint_current(
    const std::vector<double>& total_current_face_mid,
    double dt,
    int mpi_rank,
    int mpi_size)
{
    const int nxl = nx_total - 2 * Param::Nghost;
    if (dt <= 0.0 || nxl <= 0) return;
    if (Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }

    const double ampere_scale = -dt / Const::eps0;
    for (int iface = 0; iface < nxl; ++iface) {
        const size_t slot = static_cast<size_t>(iface);
        const double j_mid =
            (slot < total_current_face_mid.size())
            ? total_current_face_mid[slot] : 0.0;
        Ex_face[slot] += ampere_scale * j_mid;
    }

    close_periodic_right_face(Ex_face, nxl, mpi_rank, mpi_size, 306);
    update_cell_ex_from_faces(*this, mpi_rank, mpi_size);
}

void EMFields::sync_cell_ex_from_faces(int mpi_rank, int mpi_size)
{
    const int nxl = nx_total - 2 * Param::Nghost;
    if (nxl <= 0) return;
    if (Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }
    close_periodic_right_face(Ex_face, nxl, mpi_rank, mpi_size, 304);
    update_cell_ex_from_faces(*this, mpi_rank, mpi_size);
}

void EMFields::update_gauss_residual_diagnostics(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    if (nxl <= 0) return;
    if (Ex_face.size() != static_cast<size_t>(nxl + 1)) {
        Ex_face.assign(static_cast<size_t>(nxl + 1), 0.0);
    }
    close_periodic_right_face(Ex_face, nxl, mpi_rank, mpi_size, 303);

    double local_charge = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        local_charge += rho[ng + ix] * dx;
    }
    double global_charge = 0.0;
    MPI_Allreduce(&local_charge, &global_charge, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double mean_rho = global_charge / Param::Lx;

    double local_l1 = 0.0;
    double local_linf = 0.0;
    for (int ix = 0; ix < nxl; ++ix) {
        const double div_e =
            (Ex_face[static_cast<size_t>(ix + 1)]
           - Ex_face[static_cast<size_t>(ix)]) / dx;
        const double residual =
            div_e - (rho[ng + ix] - mean_rho) / Const::eps0;
        const double abs_residual = std::fabs(residual);
        local_l1 += abs_residual * dx;
        local_linf = std::max(local_linf, abs_residual);
    }

    MPI_Allreduce(&local_l1, &last_gauss_residual_l1, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_linf, &last_gauss_residual_linf, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    (void)mpi_rank;
    (void)mpi_size;
}

void EMFields::compute_potential(int mpi_rank, int mpi_size)
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    gather_rho(*this, mpi_size);

    if (mpi_rank == 0) {
        compute_gauss_field(*this, true);
    }

    MPI_Scatterv(global_phi.data(), counts.data(), displs.data(), MPI_DOUBLE,
                 local_rhs.data(), nxl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int ix = 0; ix < nxl; ++ix) phi[ng + ix] = local_rhs[static_cast<size_t>(ix)];
    exchange_phi_ghosts(mpi_rank, mpi_size);
}

void EMFields::exchange_ex_ghosts(int mpi_rank, int mpi_size)
{
    exchange_scalar_ghosts(*this, Ex, 201, mpi_rank, mpi_size);
}

void EMFields::exchange_phi_ghosts(int mpi_rank, int mpi_size)
{
    exchange_scalar_ghosts(*this, phi, 211, mpi_rank, mpi_size);
}

double EMFields::total_energy() const
{
    const int ng = Param::Nghost;
    const int nxl = nx_total - 2 * ng;
    double energy = 0.0;
    for (int i = ng; i < ng + nxl; ++i) {
        energy += 0.5 * Const::eps0 * Ex[i] * Ex[i];
    }
    return energy * dx;
}
