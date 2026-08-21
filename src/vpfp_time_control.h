#ifndef VPFP_TIME_CONTROL_H
#define VPFP_TIME_CONTROL_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace VpfpTimeControl {

inline double roundoff_tolerance(double time,
                                 double target,
                                 double nominal_dt)
{
    const double scale = std::max(
        std::max(std::fabs(time), std::fabs(target)),
        std::max(std::fabs(nominal_dt), std::numeric_limits<double>::min()));
    return 64.0 * std::numeric_limits<double>::epsilon() * scale;
}

inline bool should_advance(double time, double stop_time, double nominal_dt)
{
    return stop_time - time > roundoff_tolerance(time, stop_time, nominal_dt);
}

inline bool target_reached(double time, double target, double nominal_dt)
{
    return time + roundoff_tolerance(time, target, nominal_dt) >= target;
}

} // namespace VpfpTimeControl

#endif
