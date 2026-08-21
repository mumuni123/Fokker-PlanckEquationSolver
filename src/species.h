#ifndef SPECIES_H
#define SPECIES_H

#include "parameters.h"
#include "grid.h"
#include <cstdint>
#include <string>
#include <vector>

enum class SpeciesType {
    BEAM,
    BACKGROUND_ELECTRON,
    BACKGROUND_ION
};

// One cell-integrated bulk mass request for the conservative bulk-to-tail
// converter (sections 7.1 and 14.3).  ix_global is the global spatial cell
// index; iv/imu are the cylindrical (u_parallel, u_perp) cell indices; mass
// is the cell-integrated mass [m^-2].
struct ConversionMassRequest {
    int ix_global;
    int iv;
    int imu;
    double mass;
};

// Transactional inverse of ConversionMassRequest.  The mass is added to a
// locally owned Eulerian cell only after every request has passed validation.
struct ReturnMassRequest {
    int ix_global;
    int iv;
    int imu;
    double mass;
};

// Result of Species::extract_conversion_masses: the actually removed
// cell-integrated moments (section 7.6 units).  On any invalid request
// `valid`/`complete` are false and f is left untouched.
struct SpeciesConversionResult {
    bool valid;
    bool complete;
    double number_removed;
    double px_removed;
    double energy_removed;
    double jx_dx_removed;
    double pixx_dx_removed;
    double piperp_dx_removed;
    // L1 scales (sum of |per-request contributions|) for residual
    // normalisation: the signed totals of px/jx can cancel between positive
    // and negative u_parallel groups, so the fidelity/hard gates must be
    // measured against the scale of the removed contributions (section
    // 19.1), not the near-zero signed sums.
    double number_scale;
    double px_scale;
    double energy_scale;
    double jx_scale;
    double pixx_scale;
    double piperp_scale;
    SpeciesConversionResult()
        : valid(false), complete(false), number_removed(0.0),
          px_removed(0.0), energy_removed(0.0), jx_dx_removed(0.0),
          pixx_dx_removed(0.0), piperp_dx_removed(0.0),
          number_scale(0.0), px_scale(0.0), energy_scale(0.0),
          jx_scale(0.0), pixx_scale(0.0), piperp_scale(0.0)
    {}
};

struct SpeciesReturnResult {
    bool valid;
    bool complete;
    double number_added;
    double px_added;
    double energy_added;
    double jx_dx_added;
    double pixx_dx_added;
    double piperp_dx_added;
    double number_scale;
    double px_scale;
    double energy_scale;
    double jx_scale;
    double pixx_scale;
    double piperp_scale;
    SpeciesReturnResult()
        : valid(false), complete(false), number_added(0.0), px_added(0.0),
          energy_added(0.0), jx_dx_added(0.0), pixx_dx_added(0.0),
          piperp_dx_added(0.0), number_scale(0.0), px_scale(0.0),
          energy_scale(0.0), jx_scale(0.0), pixx_scale(0.0),
          piperp_scale(0.0)
    {}
};

class Species {
public:
    std::string name;
    SpeciesType type;
    double charge;
    double mass;
    double density0;
    double temperature;
    bool collisions_enabled;
    bool relativistic_push;

    VelocityGrid vgrid;
    CylindricalVelocityGrid cgrid;
    const SpatialGrid* sgrid;

    // Background electrons use the phase-1 cylindrical representation.
    // In that mode f stores the cell-integrated conservative mass M, not a
    // point value of the legacy spherical distribution.
    bool cylindrical_mass_representation;

    // Axisymmetric distribution f(x, u, mu), with d3u = 2*pi*u^2 du dmu.
    std::vector<double> f;
    std::vector<double> f_tmp;

    std::vector<double> number_density;
    std::vector<double> charge_density;
    std::vector<double> current_x;
    std::vector<double> current_face_x;

    Species();

    void init(const std::string& name, SpeciesType type,
              double charge, double mass,
              double density, double temperature,
              bool collisions,
              const SpatialGrid& sg);

    void initialize_maxwellian(double drift_vx = 0.0);
    void initialize_maxwellian_profile(const std::vector<double>& density_profile,
                                       double drift_vx = 0.0);
    double maxwellian_f_value(double density,
                              double temperature,
                              double drift_vx,
                              int iv,
                              int imu) const;
    void fill_maxwellian_velocity_slice(std::vector<double>& values,
                                        double density,
                                        double temperature,
                                        double drift_vx) const;
    // Returns one cylindrical conservative-mass slice M(u_parallel,u_perp)
    // for a single spatial cell.  Open-boundary code uses this routine for
    // reservoir ghosts so the Maxwellian normalization is not duplicated in
    // transport operators.
    void fill_cylindrical_maxwellian_mass_slice(std::vector<double>& values,
                                                double density,
                                                double temperature,
                                                double drift_vx) const;
    void compute_moments();
    double total_particle_number() const;
    double total_kinetic_energy() const;
    void total_particle_number_and_energy(double& number,
                                          double& kinetic_energy) const;
    double distribution_value(int ix_with_ghost, int iv, int imu) const;
    // Transactional accept: exchange the phase-space mass and moment arrays
    // with a same-layout working species without copying the full array.
    void swap_state(Species& other);

    // Safe batch conversion extraction (section 14.3): validates every
    // request first (locally owned cell, in-range velocity indices,
    // finite non-negative mass not exceeding the stored mass), then
    // subtracts exactly the requested masses from f.  The returned removed
    // moments are computed with the single shared mass_cell_moments formula.
    // On any invalid request nothing is modified.
    SpeciesConversionResult extract_conversion_masses(
        const std::vector<ConversionMassRequest>& requests);

    // H10 inverse transaction.  All requests are validated before f is
    // changed, so a rejected tail-to-bulk projection cannot partially alter
    // the trial representation.
    SpeciesReturnResult add_return_masses(
        const std::vector<ReturnMassRequest>& requests);
    bool validate_return_masses(
        const std::vector<ReturnMassRequest>& requests) const;

    size_t local_size() const {
        return static_cast<size_t>(sgrid->nx_total) * Param::Nvmu;
    }
};

#endif
