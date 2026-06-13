#ifndef MAXWELL_H
#define MAXWELL_H

#include "grid.h"
#include "parameters.h"
#include <vector>

struct EMFields {
    std::vector<double> Ex;
    std::vector<double> phi;
    std::vector<double> rho;
    std::vector<double> send_left;
    std::vector<double> send_right;
    std::vector<double> recv_left;
    std::vector<double> recv_right;

    std::vector<int> counts;
    std::vector<int> displs;
    std::vector<double> local_rhs;
    std::vector<double> global_rhs;
    std::vector<double> global_ex;
    std::vector<double> global_phi;
    std::vector<double> all_interfaces;

    int nx_total;
    int counts_mpi_size;
    int counts_nx_local;
    double dx;

    void init(const SpatialGrid& sg);
    void zero_currents();
    void accumulate_moments(const class Species& sp);
    void set_charge_density(const class Species& electrons,
                            const std::vector<double>& beam_density,
                            const std::vector<double>& ion_density_profile);

    // Periodic Gauss solve for initialization or low-frequency correction.
    void solve_poisson(int mpi_rank, int mpi_size);
    void advance_ampere(const std::vector<double>& background_current,
                        const std::vector<double>& beam_current,
                        double dt,
                        int mpi_rank,
                        int mpi_size);
    // Snapshot-only potential reconstruction.
    void compute_potential(int mpi_rank, int mpi_size);
    void exchange_ex_ghosts(int mpi_rank, int mpi_size);
    void exchange_phi_ghosts(int mpi_rank, int mpi_size);

    double total_energy() const;
};

#endif
