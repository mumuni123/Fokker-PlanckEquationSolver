#ifndef SPECIES_H
#define SPECIES_H

#include "parameters.h"
#include "grid.h"
#include <vector>
#include <string>

enum class SpeciesType {
    BEAM,
    BACKGROUND_ELECTRON,
    BACKGROUND_ION
};

class Species {
public:
    std::string name;
    SpeciesType type;
    double charge;    // 电荷量（库仑，带符号）
    double mass;      // 质量（kg）
    double density0;  // 参考数密度 [/m^3]
    double temperature; // 温度（焦耳）
    bool collisions_enabled;

    MomentumGrid pgrid;
    const SpatialGrid* sgrid;

    // 分布函数：f[ix_total * Np3]
    // ix 取值范围 [0, nx_total)，包含幽灵单元
    std::vector<double> f;

    // 空间网格上的宏观矩（局部区域，输出时不需要幽灵单元）
    std::vector<double> number_density;   // n(x)
    std::vector<double> charge_density;   // rho(x)
    std::vector<double> current_x;        // Jx(x)
    std::vector<double> current_y;        // Jy(x)
    std::vector<double> current_z;        // Jz(x)

    Species();

    void init(const std::string& name, SpeciesType type,
              double charge, double mass,
              double density, double temperature,
              bool collisions,
              const SpatialGrid& sg);

    void initialize_maxwellian(double drift_px = 0.0);
    void compute_moments();
    double total_particle_number() const;
    double total_kinetic_energy() const;

    size_t local_size() const {
        return static_cast<size_t>(sgrid->nx_total) * Param::Np3;
    }
};

#endif // 头文件保护：SPECIES_H
