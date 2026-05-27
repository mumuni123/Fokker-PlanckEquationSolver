#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include <fstream>
#include <string>
#include <vector>

class BeamPIC;

class Diagnostics {
public:
    std::string output_dir;
    int snapshot_count;

    Diagnostics();

    void init(const std::string& dir, int mpi_rank,
              bool enable_debug_diagnostics,
              bool enable_step_diagnostics);

    void write_scalars(double time, int step,
                       const std::vector<Species*>& species,
                       const EMFields& fields,
                       int mpi_rank, int mpi_size);

    void write_debug_state(int step, double time,
                           const std::string& stage,
                           const Species& electrons,
                           const BeamPIC& beam,
                           const EMFields& fields,
                           const SpatialGrid& sg,
                           int mpi_rank, int mpi_size,
                           double cfl_v = 0.0,
                           double cfl_mu = 0.0,
                           int nsub_v = 0,
                           int nsub_mu = 0);

    void write_step_diagnostics(int step, double time,
                                const Species& electrons,
                                const BeamPIC& beam,
                                const EMFields& fields,
                                const SpatialGrid& sg,
                                int mpi_rank, int mpi_size,
                                int nsub_v1,
                                int nsub_mu1,
                                int nsub_v2,
                                int nsub_mu2,
                                double loss_v1,
                                double loss_mu1,
                                double loss_v2,
                                double loss_mu2,
                                double loss_v1_low,
                                double loss_v1_high,
                                double loss_v2_low,
                                double loss_v2_high);

    void write_fields(double time,
                      const EMFields& fields,
                      const SpatialGrid& sg,
                      int mpi_rank, int mpi_size);

    void write_px_distribution(double time,
                               const Species& sp,
                               int mpi_rank, int mpi_size);

    void write_density_profile(double time,
                               const Species& electrons,
                               const std::vector<double>& beam_density,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size);

    void write_electron_distribution(double time,
                                     const Species& electrons,
                                     const SpatialGrid& sg,
                                     int mpi_rank);

    void advance_snapshot();

private:
    std::ofstream scalar_file;
#if FP_ENABLE_DEBUG_DIAGNOSTICS
    std::ofstream debug_file;
#endif
    std::ofstream step_file;
    bool debug_enabled;
    bool step_enabled;
};

#endif
