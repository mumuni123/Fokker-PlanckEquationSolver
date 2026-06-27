#ifndef MAXWELL_H
#define MAXWELL_H

#include "grid.h"
#include "parameters.h"
#include <vector>

struct EMFields {
    std::vector<double> Ex;
    std::vector<double> Ex_face;
    std::vector<double> phi;
    std::vector<double> rho;
    std::vector<double> send_left;
    std::vector<double> send_right;
    std::vector<double> recv_left;
    std::vector<double> recv_right;

    std::vector<int> counts;
    std::vector<int> displs;
    std::vector<int> face_counts;
    std::vector<int> face_displs;
    std::vector<double> local_rhs;
    std::vector<double> local_face_rhs;
    std::vector<double> global_rhs;
    std::vector<double> global_ex;
    std::vector<double> global_face;
    std::vector<double> global_phi;
    std::vector<double> all_interfaces;

    int nx_total;
    int counts_mpi_size;
    int counts_nx_local;
    double dx;
    double last_gauss_residual_l1;
    double last_gauss_residual_linf;

    void init(const SpatialGrid& sg);
    void zero_currents();
    void accumulate_moments(const class Species& sp);
    void set_charge_density(const class Species& electrons,
                            const std::vector<double>& beam_density,
                            const std::vector<double>& ion_density_profile);

    // Periodic Gauss solve for initialization or low-frequency correction.
    void solve_poisson(int mpi_rank, int mpi_size);
    // Ex_face[nxl] is a periodic ghost. Beam face 0 and face nxl are open
    // boundary currents and are never averaged or periodically identified.
    void advance_ampere_face(const std::vector<double>& background_current_face,
                             const std::vector<double>& open_beam_current_face,
                             double dt,
                             int mpi_rank,
                             int mpi_size);
    void sync_cell_ex_from_faces(int mpi_rank, int mpi_size);
    void update_gauss_residual_diagnostics(int mpi_rank, int mpi_size);
    // Snapshot-only potential reconstruction.
    void compute_potential(int mpi_rank, int mpi_size);
    void exchange_ex_ghosts(int mpi_rank, int mpi_size);
    void exchange_phi_ghosts(int mpi_rank, int mpi_size);

    double total_energy() const;
};

#endif
