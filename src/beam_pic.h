#ifndef BEAM_PIC_H
#define BEAM_PIC_H

#include "grid.h"
#include "maxwell.h"
#include <vector>

struct BeamParticle {
    double x;
    double px;
    double weight;
};

class BeamPIC {
public:
    std::vector<BeamParticle> particles;
    std::vector<double> density;
    std::vector<double> current_x;

    BeamPIC();

    void init(const SpatialGrid& sg);
    void inject(const SpatialGrid& sg, double dt, double time, int mpi_rank);
    void push(const SpatialGrid& sg, const EMFields& fields, double dt,
              int mpi_rank, int mpi_size);
    void deposit_density(const SpatialGrid& sg, int mpi_rank, int mpi_size);

    double total_particle_number(const SpatialGrid& sg) const;
    double total_kinetic_energy() const;

private:
    double injection_remainder_;
    std::vector<BeamParticle> send_left_;
    std::vector<BeamParticle> send_right_;
    std::vector<BeamParticle> keep_;
    std::vector<BeamParticle> recv_left_;
    std::vector<BeamParticle> recv_right_;
    std::vector<std::vector<BeamParticle> > thread_keep_;
    std::vector<std::vector<BeamParticle> > thread_send_left_;
    std::vector<std::vector<BeamParticle> > thread_send_right_;
    std::vector<std::vector<double> > thread_density_;
    std::vector<std::vector<double> > thread_current_;
    std::vector<double> thread_send_left_density_;
    std::vector<double> thread_send_right_density_;
    std::vector<double> thread_send_left_current_;
    std::vector<double> thread_send_right_current_;

    void exchange_particles(const SpatialGrid& sg, int mpi_rank, int mpi_size);
};

#endif
