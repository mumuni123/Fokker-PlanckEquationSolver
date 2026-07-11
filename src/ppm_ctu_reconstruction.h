#ifndef PPM_CTU_RECONSTRUCTION_H
#define PPM_CTU_RECONSTRUCTION_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace PpmCtu {

struct Parabola {
    double ql;
    double qr;
    double dq;
    double q6;
    double theta_positive;
};

inline double value(const Parabola& p, double xi)
{
    return p.ql + xi * (p.dq + p.q6 * (1.0 - xi));
}

inline Parabola reconstruct(double qmm, double qm, double q0,
                            double qp, double qpp)
{
    // The production cylindrical grid is uniform in x and u_parallel.  The
    // fourth-order interface estimates below are converted to a full PPM
    // polynomial and are never used directly as the final face state.
    Parabola p;
    p.theta_positive = 1.0;
    p.ql = (7.0 * (qm + q0) - qmm - qp) / 12.0;
    p.qr = (7.0 * (q0 + qp) - qm - qpp) / 12.0;
    p.ql = std::max(std::min(qm, q0), std::min(std::max(qm, q0), p.ql));
    p.qr = std::max(std::min(q0, qp), std::min(std::max(q0, qp), p.qr));

    p.dq = p.qr - p.ql;
    p.q6 = 6.0 * q0 - 3.0 * (p.ql + p.qr);
    if ((p.qr - q0) * (q0 - p.ql) <= 0.0) {
        p.ql = q0;
        p.qr = q0;
    } else if (p.dq * p.q6 > p.dq * p.dq) {
        p.ql = 3.0 * q0 - 2.0 * p.qr;
    } else if (p.dq * p.q6 < -p.dq * p.dq) {
        p.qr = 3.0 * q0 - 2.0 * p.ql;
    }
    p.dq = p.qr - p.ql;
    p.q6 = 6.0 * q0 - 3.0 * (p.ql + p.qr);

    double qmin = std::min(p.ql, p.qr);
    if (std::fabs(p.q6) > std::numeric_limits<double>::min()) {
        const double vertex = 0.5 * (p.dq + p.q6) / p.q6;
        if (vertex > 0.0 && vertex < 1.0)
            qmin = std::min(qmin, value(p, vertex));
    }
    if (qmin < 0.0 && q0 > 0.0) {
        p.theta_positive = std::min(1.0, q0 / (q0 - qmin));
        p.ql = q0 + p.theta_positive * (p.ql - q0);
        p.qr = q0 + p.theta_positive * (p.qr - q0);
        p.dq = p.qr - p.ql;
        p.q6 = 6.0 * q0 - 3.0 * (p.ql + p.qr);
    }
    return p;
}

inline double trace_right(const Parabola& p, double sigma)
{
    if (sigma <= 0.0) return p.qr;
    sigma = std::min(1.0, sigma);
    return p.qr - 0.5 * sigma *
        (p.dq - (1.0 - 2.0 * sigma / 3.0) * p.q6);
}

inline double trace_left(const Parabola& p, double sigma)
{
    if (sigma <= 0.0) return p.ql;
    sigma = std::min(1.0, sigma);
    return p.ql + 0.5 * sigma *
        (p.dq + (1.0 - 2.0 * sigma / 3.0) * p.q6);
}

inline double upwind_state(double left, double right, double speed)
{
    return (speed > 0.0) ? left : ((speed < 0.0) ? right : 0.5 * (left + right));
}

} // namespace PpmCtu

#endif
