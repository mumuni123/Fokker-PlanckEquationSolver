#ifndef PERIODIC_STAGGERED_OPERATORS_H
#define PERIODIC_STAGGERED_OPERATORS_H

#include <vector>

// Production periodic Yee stagger operators.  Faces 0..nxl-1 are owned;
// face nxl is the right-neighbour alias and is never a second physical DOF.
namespace PeriodicStaggered {

void close_right_face_alias(std::vector<double>& face, int nxl, int block,
                            int mpi_rank, int mpi_size, int message_tag);

// G: unique face field to local cell-center field.
void apply_face_to_cell_G(const std::vector<double>& face,
                          std::vector<double>& cell, int nxl);
void apply_face_to_cell_G(const std::vector<double>& face,
                          double* cell, int nxl);

// G*: local cell-center quantity to face quantity, including the periodic
// right alias exchange.  The same definition is used by production J_E and
// the synthetic adjoint tests.
void apply_cell_to_face_Gstar(const std::vector<double>& cell,
                              std::vector<double>& face, int nxl,
                              int mpi_rank, int mpi_size, int message_tag);

} // namespace PeriodicStaggered

#endif
