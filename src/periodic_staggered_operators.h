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

// Audit-only form of G*.  `independent_face` evaluates both local endpoint
// faces from neighbouring cell values before the right face is replaced by
// its periodic alias; `closed_face` is the production periodic result.
void audit_cell_to_face_Gstar_sync(const std::vector<double>& cell,
                                   std::vector<double>& independent_face,
                                   std::vector<double>& closed_face, int nxl,
                                   int mpi_rank, int mpi_size,
                                   int message_tag);

// Audit-only block form of G*.  `cell` is cell-major with `block` values per
// spatial cell; `face` uses the corresponding local Yee-face layout and has
// its right periodic alias closed before return.
void audit_cell_blocks_to_face_Gstar(const std::vector<double>& cell,
                                     std::vector<double>& face, int nxl,
                                     int block, int mpi_rank, int mpi_size,
                                     int message_tag);

} // namespace PeriodicStaggered

#endif
