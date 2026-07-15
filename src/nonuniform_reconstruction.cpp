#include "nonuniform_reconstruction.h"

#include <algorithm>
#include <cmath>

namespace {

double minmod3(double a, double b, double c)
{
    if ((a > 0.0 && b > 0.0 && c > 0.0) ||
        (a < 0.0 && b < 0.0 && c < 0.0)) {
        const double sign = (a > 0.0) ? 1.0 : -1.0;
        return sign * std::min(std::fabs(a), std::min(std::fabs(b), std::fabs(c)));
    }
    return 0.0;
}

double mc_slope(double q_left, double q_center, double q_right,
                double s_left, double s_center, double s_right)
{
    const double ds_left = s_center - s_left;
    const double ds_right = s_right - s_center;
    if (!(ds_left > 0.0) || !(ds_right > 0.0)) return 0.0;

    const double d_left = (q_center - q_left) / ds_left;
    const double d_right = (q_right - q_center) / ds_right;
    return minmod3(0.5 * (d_left + d_right), 2.0 * d_left, 2.0 * d_right);
}

} // namespace

namespace NonuniformMuscl {

FaceStates reconstruct_face(double q_im1, double q_i, double q_ip1,
                            double q_ip2, double s_im1, double s_i,
                            double s_ip1, double s_ip2, double s_face)
{
    const double slope_i = mc_slope(q_im1, q_i, q_ip1,
                                    s_im1, s_i, s_ip1);
    const double slope_ip1 = mc_slope(q_i, q_ip1, q_ip2,
                                      s_i, s_ip1, s_ip2);
    FaceStates result;
    result.left = q_i + slope_i * (s_face - s_i);
    result.right = q_ip1 + slope_ip1 * (s_face - s_ip1);
    return result;
}

double upwind_state(const FaceStates& states, double speed)
{
    return (speed > 0.0) ? states.left :
           ((speed < 0.0) ? states.right : 0.5 * (states.left + states.right));
}

double centered_state(const FaceStates& states)
{
    return 0.5 * (states.left + states.right);
}

double upar_face_coefficient(double fbar, double dx, double transverse_area)
{
    return fbar * dx * transverse_area;
}

double upar_face_flux(double acceleration, double coefficient)
{
    return acceleration * coefficient;
}

} // namespace NonuniformMuscl
