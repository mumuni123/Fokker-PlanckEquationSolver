#ifndef SPECIES_H
#define SPECIES_H

#include "parameters.h"
#include "grid.h"
#include <string>
#include <vector>

enum class SpeciesType {
    BEAM,
    BACKGROUND_ELECTRON,
    BACKGROUND_ION
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
    const SpatialGrid* sgrid;

    // Axisymmetric distribution f(x, v, mu), with d3v = 2*pi*v^2 dv dmu.
    std::vector<double> f;
    std::vector<double> f_tmp;

    std::vector<double> number_density;
    std::vector<double> charge_density;
    std::vector<double> current_x;
    std::vector<double> maxwellian_left_boundary;
    std::vector<double> maxwellian_right_boundary;

    Species();

    void init(const std::string& name, SpeciesType type,
              double charge, double mass,
              double density, double temperature,
              bool collisions,
              const SpatialGrid& sg);

    void initialize_maxwellian(double drift_vx = 0.0);
    void set_maxwellian_boundary_drifts(double left_drift_vx,
                                         double right_drift_vx);
    void apply_vx_shift_profile(const std::vector<double>& delta_vx);
    void compute_moments();
    double total_particle_number() const;
    double total_kinetic_energy() const;
    void total_particle_number_and_energy(double& number,
                                          double& kinetic_energy) const;

    size_t local_size() const {
        return static_cast<size_t>(sgrid->nx_total) * Param::Nvmu;
    }
};

#endif
