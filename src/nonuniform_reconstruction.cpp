#include "nonuniform_reconstruction.h"

#include <algorithm>
#include <cmath>

namespace {

double minmod3(double a, double b, double c, int* branch = 0)
{
    if ((a > 0.0 && b > 0.0 && c > 0.0) ||
        (a < 0.0 && b < 0.0 && c < 0.0)) {
        const double sign = (a > 0.0) ? 1.0 : -1.0;
        const double aa = std::fabs(a);
        const double bb = std::fabs(b);
        const double cc = std::fabs(c);
        if (branch) {
            *branch = (aa <= bb && aa <= cc) ? 1 :
                      ((bb <= cc) ? 2 : 3);
        }
        return sign * std::min(aa, std::min(bb, cc));
    }
    if (branch) *branch = 0;
    return 0.0;
}

double mc_slope(double q_left, double q_center, double q_right,
                double s_left, double s_center, double s_right,
                int* branch = 0)
{
    const double ds_left = s_center - s_left;
    const double ds_right = s_right - s_center;
    if (!(ds_left > 0.0) || !(ds_right > 0.0)) {
        if (branch) *branch = 0;
        return 0.0;
    }

    const double d_left = (q_center - q_left) / ds_left;
    const double d_right = (q_right - q_center) / ds_right;
    return minmod3(0.5 * (d_left + d_right), 2.0 * d_left, 2.0 * d_right,
                   branch);
}

void slope_linear_coefficients(int branch, double ds_left, double ds_right,
                               double (&coefficient)[3])
{
    coefficient[0] = coefficient[1] = coefficient[2] = 0.0;
    if (!(ds_left > 0.0) || !(ds_right > 0.0)) return;
    if (branch == 1) {
        coefficient[0] = -0.5 / ds_left;
        coefficient[1] = 0.5 / ds_left - 0.5 / ds_right;
        coefficient[2] = 0.5 / ds_right;
    } else if (branch == 2) {
        coefficient[0] = -2.0 / ds_left;
        coefficient[1] = 2.0 / ds_left;
    } else if (branch == 3) {
        coefficient[1] = -2.0 / ds_right;
        coefficient[2] = 2.0 / ds_right;
    }
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

int centered_state_branch_signature(double q_im1, double q_i, double q_ip1,
                                    double q_ip2, double s_im1, double s_i,
                                    double s_ip1, double s_ip2)
{
    int left_branch = 0;
    int right_branch = 0;
    mc_slope(q_im1, q_i, q_ip1, s_im1, s_i, s_ip1, &left_branch);
    mc_slope(q_i, q_ip1, q_ip2, s_i, s_ip1, s_ip2, &right_branch);
    return left_branch | (right_branch << 2);
}

FrozenCenteredLinearization frozen_centered_linearization(
    double q_im1, double q_i, double q_ip1, double q_ip2, double s_im1,
    double s_i, double s_ip1, double s_ip2, double s_face)
{
    int left_branch = 0;
    int right_branch = 0;
    mc_slope(q_im1, q_i, q_ip1, s_im1, s_i, s_ip1, &left_branch);
    mc_slope(q_i, q_ip1, q_ip2, s_i, s_ip1, s_ip2, &right_branch);
    double left_slope[3];
    double right_slope[3];
    slope_linear_coefficients(left_branch, s_i - s_im1, s_ip1 - s_i,
                              left_slope);
    slope_linear_coefficients(right_branch, s_ip1 - s_i, s_ip2 - s_ip1,
                              right_slope);
    FrozenCenteredLinearization result = {};
    const double left_offset = s_face - s_i;
    const double right_offset = s_face - s_ip1;
    result.coefficient[0] = 0.5 * left_offset * left_slope[0];
    result.coefficient[1] = 0.5 * (1.0 + left_offset * left_slope[1] +
                                   right_offset * right_slope[0]);
    result.coefficient[2] = 0.5 * (left_offset * left_slope[2] + 1.0 +
                                   right_offset * right_slope[1]);
    result.coefficient[3] = 0.5 * right_offset * right_slope[2];
    result.branch_signature = left_branch | (right_branch << 2);
    return result;
}

FrozenCenteredLinearization frozen_upwind_linearization(
    double q_im1, double q_i, double q_ip1, double q_ip2, double s_im1,
    double s_i, double s_ip1, double s_ip2, double s_face,
    bool use_left_state)
{
    int left_branch = 0;
    int right_branch = 0;
    mc_slope(q_im1, q_i, q_ip1, s_im1, s_i, s_ip1, &left_branch);
    mc_slope(q_i, q_ip1, q_ip2, s_i, s_ip1, s_ip2, &right_branch);
    double left_slope[3];
    double right_slope[3];
    slope_linear_coefficients(left_branch, s_i - s_im1, s_ip1 - s_i,
                              left_slope);
    slope_linear_coefficients(right_branch, s_ip1 - s_i, s_ip2 - s_ip1,
                              right_slope);
    FrozenCenteredLinearization result = {};
    if (use_left_state) {
        const double offset = s_face - s_i;
        result.coefficient[0] = offset * left_slope[0];
        result.coefficient[1] = 1.0 + offset * left_slope[1];
        result.coefficient[2] = offset * left_slope[2];
    } else {
        const double offset = s_face - s_ip1;
        result.coefficient[1] = offset * right_slope[0];
        result.coefficient[2] = 1.0 + offset * right_slope[1];
        result.coefficient[3] = offset * right_slope[2];
    }
    result.branch_signature = left_branch | (right_branch << 2);
    return result;
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
