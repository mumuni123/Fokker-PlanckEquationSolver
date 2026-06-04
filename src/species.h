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
    double reservoir_density_left;
    double reservoir_density_right;
    double reservoir_drift_left;
    double reservoir_drift_right;
    bool collisions_enabled;
    bool relativistic_push;

    VelocityGrid vgrid;
    const SpatialGrid* sgrid;

    // Axisymmetric distribution f(x, u, mu), with d3u = 2*pi*u^2 du dmu.
    std::vector<double> f;
    std::vector<double> f_tmp;

    std::vector<double> number_density;
    std::vector<double> charge_density;
    std::vector<double> current_x;

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
