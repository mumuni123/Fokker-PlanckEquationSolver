#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "grid.h"
#include "maxwell.h"
#include "species.h"
#include <fstream>
#include <string>
#include <vector>

class Diagnostics {
public:
    std::string output_dir;
    int snapshot_count;

    Diagnostics();

    void init(const std::string& dir, int mpi_rank);

    void write_scalars(double time, int step,
                       const std::vector<Species*>& species,
                       const EMFields& fields,
                       int mpi_rank, int mpi_size);

    void write_fields(double time,
                      const EMFields& fields,
                      const SpatialGrid& sg,
                      int mpi_rank, int mpi_size);

    void write_px_distribution(double time,
                               const Species& sp,
                               int mpi_rank, int mpi_size);

    void write_density_profile(double time,
                               const std::vector<Species*>& species,
                               const SpatialGrid& sg,
                               int mpi_rank, int mpi_size);

private:
    std::ofstream scalar_file;
};

#endif
