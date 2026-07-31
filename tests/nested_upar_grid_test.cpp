#include "grid.h"

#include <cmath>
#include <cstdio>
#include <limits>

int main()
{
    CylindricalVelocityGrid grid;
    try {
        grid.init(Param::momentum_umax);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "nested_upar_grid_test: %s\n", error.what());
        return 1;
    }

    const CylindricalVelocityGrid::NestedGridAudit& audit = grid.nested_audit;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, Param::momentum_upar_extended_max);
    const bool passes = audit.core_face_identity_linf <= tolerance &&
                        audit.core_cell_identity_linf <= tolerance &&
                        audit.core_vx_identity_linf <= tolerance * Const::c &&
                        audit.core_kinetic_energy_identity_linf <=
                            tolerance * Const::me * Const::c * Const::c &&
                        audit.symmetry_linf <= tolerance &&
                        audit.max_adjacent_width_ratio <= 1.25 + tolerance &&
                        audit.phase_volume_relative_error <= 1.0e-13;
    std::printf("Nv_core=%d Nv_tail=%d Nv_total=%d\n",
                Param::Nv_core, Param::Nv_tail, Param::Nv);
    std::printf("core_face_identity_linf=%.17e\n", audit.core_face_identity_linf);
    std::printf("core_cell_identity_linf=%.17e\n", audit.core_cell_identity_linf);
    std::printf("core_vx_identity_linf=%.17e\n", audit.core_vx_identity_linf);
    std::printf("core_kinetic_energy_identity_linf=%.17e\n",
                audit.core_kinetic_energy_identity_linf);
    std::printf("symmetry_linf=%.17e\n", audit.symmetry_linf);
    std::printf("max_adjacent_width_ratio=%.17e\n",
                audit.max_adjacent_width_ratio);
    std::printf("phase_volume_relative_error=%.17e\n",
                audit.phase_volume_relative_error);
    std::printf("passes=%d\n", passes ? 1 : 0);
    return passes ? 0 : 1;
}
