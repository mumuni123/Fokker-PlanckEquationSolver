#include "vlasov_split_step.h"

#include <mpi.h>

VlasovSplitStep::VlasovSplitStep()
    : boundary_(NULL), mpi_rank_(0), mpi_size_(1), initialized_(false)
{}

void VlasovSplitStep::init(const SpatialGrid& grid, const Species& prototype,
                           const OpenBackgroundBoundary& boundary)
{
    grid_ = grid;
    boundary_ = &boundary;
    remap_.init(grid, prototype.cgrid);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
    initialized_ = true;
}

bool VlasovSplitStep::first_x_half(const Species& state_n,
                                   Species& state_x_half, double time,
                                   double half_dt, VlasovStepDiagnostics& diag,
                                   XFaceTransportAudit* audit)
{
    if (!initialized_) return false;
    diag.x_first = remap_.advect_x(state_n, state_x_half, half_dt, time,
                                   *boundary_, mpi_rank_, mpi_size_, audit);
    return diag.x_first.finite;
}

bool VlasovSplitStep::u_full(const Species& state_x_half,
                             Species& state_u_full, const EMFields& field_mid,
                             double time_mid, double dt,
                             VlasovStepDiagnostics& diag,
                             const HybridVelocityPartition* partition,
                             BulkTailFluxBatch* exported_flux,
                             int quadrature_order,
                             std::vector<double>* local_delta_ke_by_x)
{
    if (!initialized_) return false;
    diag.u_full = remap_.advect_u_parallel(state_x_half, state_u_full,
                                           field_mid, dt, time_mid,
                                           partition, exported_flux,
                                           quadrature_order,
                                           local_delta_ke_by_x);
    if (!diag.u_full.audit_physical_state_bitwise_equal) {
        diag.u_full.audit_valid = false;
        diag.u_full.audit_failure_code |= 32;
    }
    return diag.u_full.finite;
}

bool VlasovSplitStep::second_x_half(const Species& state_u_full,
                                    Species& state_np1, double time_mid,
                                    double half_dt,
                                    VlasovStepDiagnostics& diag,
                                    XFaceTransportAudit* audit)
{
    if (!initialized_) return false;
    diag.x_second = remap_.advect_x(state_u_full, state_np1, half_dt,
                                    time_mid, *boundary_, mpi_rank_,
                                    mpi_size_, audit);
    return diag.x_second.finite;
}
