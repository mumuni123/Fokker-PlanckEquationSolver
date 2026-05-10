#ifndef COLLISION_H
#define COLLISION_H

#include "species.h"

class CollisionOperator {
public:
    void apply(Species& sp, double dt,
               double n_field, double T_field, double m_field,
               double Z_field, double Z_test);

private:
    struct CollisionRates {
        double nu_s;
        double nu_parallel;
        double nu_perp;
    };

    double coulomb_log(double n_e, double T_e) const;
    CollisionRates compute_rates(double v, double mass_test,
                                 double n_field, double T_field,
                                 double m_field, double lnLambda,
                                 double Z_test, double Z_field) const;
    static double chandrasekhar_G(double x);
};

#endif
