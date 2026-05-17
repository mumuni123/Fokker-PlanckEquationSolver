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
    std::vector<double> global_phi;
    std::vector<double> tri_a;
    std::vector<double> tri_b;
    std::vector<double> tri_c;
    std::vector<double> tri_rhs;
    std::vector<double> tri_y;
    std::vector<double> tri_left_response;
    std::vector<double> tri_right_response;
    std::vector<double> tri_work_a;
    std::vector<double> tri_work_b;
    std::vector<double> tri_work_c;
    std::vector<double> tri_work_d;
    std::vector<double> all_interfaces;
    std::vector<double> interface_values;

    int nx_total;
    int counts_mpi_size;
    int counts_nx_local;
    double dx;

    void init(const SpatialGrid& sg);
    void zero_currents();
    void accumulate_moments(const class Species& sp);
    void set_charge_density(const class Species& electrons,
                            const std::vector<double>& beam_density);

    // Electrostatic solve: phi(0) = phi(L) = 0, Ex = -dphi/dx.
    void solve_poisson(int mpi_rank, int mpi_size);
    void exchange_ghosts(int mpi_rank, int mpi_size);

    double total_energy() const;
};

#endif
