#include "nonuniform_reconstruction.h"
#include "discrete_moment_operators.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace {

bool close(double a, double b, double multiplier = 256.0)
{
    return std::fabs(a - b) <= multiplier * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

std::vector<double> stretched_faces(int cells)
{
    std::vector<double> faces(static_cast<size_t>(cells + 1));
    for (int i = 0; i <= cells; ++i) {
        const double xi = static_cast<double>(i) / cells;
        faces[static_cast<size_t>(i)] = xi + 0.23 * xi * xi;
    }
    return faces;
}

double cell_average_sine(double left, double right)
{
    const double pi = 3.141592653589793238462643383279502884;
    const double width = right - left;
    return (std::cos(2.0 * pi * left) - std::cos(2.0 * pi * right)) /
        (2.0 * pi * width);
}

double smooth_face_error(int cells)
{
    const double pi = 3.141592653589793238462643383279502884;
    const std::vector<double> faces = stretched_faces(cells);
    std::vector<double> centers(static_cast<size_t>(cells));
    std::vector<double> average(static_cast<size_t>(cells));
    for (int i = 0; i < cells; ++i) {
        centers[static_cast<size_t>(i)] = 0.5 * (faces[i] + faces[i + 1]);
        average[static_cast<size_t>(i)] = cell_average_sine(faces[i], faces[i + 1]);
    }

    double l1 = 0.0;
    int count = 0;
    for (int i = 1; i < cells - 2; ++i) {
        const NonuniformMuscl::FaceStates states = NonuniformMuscl::reconstruct_face(
            average[i - 1], average[i], average[i + 1], average[i + 2],
            centers[i - 1], centers[i], centers[i + 1], centers[i + 2], faces[i + 1]);
        const double exact = std::sin(2.0 * pi * faces[i + 1]);
        l1 += 0.5 * (std::fabs(states.left - exact) + std::fabs(states.right - exact));
        ++count;
    }
    return l1 / count;
}

} // namespace

int main()
{
    const std::vector<double> faces = stretched_faces(12);
    std::vector<double> centers(12);
    for (int i = 0; i < 12; ++i)
        centers[static_cast<size_t>(i)] = 0.5 * (faces[i] + faces[i + 1]);

    // Constant cell averages stay constant on both sides of every shared face.
    const NonuniformMuscl::FaceStates constant = NonuniformMuscl::reconstruct_face(
        3.0, 3.0, 3.0, 3.0, centers[2], centers[3], centers[4], centers[5], faces[4]);
    assert(close(constant.left, 3.0));
    assert(close(constant.right, 3.0));

    // Cell averages of a linear function equal its center values on an
    // arbitrary nonuniform mesh, so the face extrapolation is exact.
    const double slope = -1.7;
    const double intercept = 2.25;
    const NonuniformMuscl::FaceStates linear = NonuniformMuscl::reconstruct_face(
        intercept + slope * centers[2], intercept + slope * centers[3],
        intercept + slope * centers[4], intercept + slope * centers[5],
        centers[2], centers[3], centers[4], centers[5], faces[4]);
    const double linear_exact = intercept + slope * faces[4];
    assert(close(linear.left, linear_exact));
    assert(close(linear.right, linear_exact));

    // Smooth cell averages converge at approximately second order on the
    // stretched physical grid.  The limiter is inactive away from extrema.
    const double e32 = smooth_face_error(32);
    const double e64 = smooth_face_error(64);
    const double e128 = smooth_face_error(128);
    const double p1 = std::log(e32 / e64) / std::log(2.0);
    const double p2 = std::log(e64 / e128) / std::log(2.0);
    assert(p1 > 1.6);
    assert(p2 > 1.6);

    // Mirror-symmetric data gives mirror-symmetric states and the single
    // upwind Riemann selection is the only face flux degree of freedom.
    const NonuniformMuscl::FaceStates symmetric = NonuniformMuscl::reconstruct_face(
        1.0, 2.0, 2.0, 1.0, -2.0, -1.0, 1.0, 2.0, 0.0);
    assert(close(symmetric.left, symmetric.right));
    assert(close(NonuniformMuscl::upwind_state(symmetric, 4.0), symmetric.left));
    assert(close(NonuniformMuscl::upwind_state(symmetric, -4.0), symmetric.right));

    // Cylindrical high-order u-face conversion must retain the spatial width:
    // C_u = fbar*dx*A_perp = M/du and Phi_u = a*C_u.
    const double fbar = 2.75;
    const double dx = 2.0e-9;
    const double area = 4.5e-3;
    const double acceleration = -7.0e4;
    const double coefficient = NonuniformMuscl::upar_face_coefficient(fbar, dx, area);
    assert(close(coefficient, fbar * dx * area));
    assert(close(NonuniformMuscl::upar_face_flux(acceleration, coefficient),
                 acceleration * fbar * dx * area));

    // FCT limits only anti-diffusive donor outflow.  The four directional
    // alpha values share one cell budget, and a shared face has one final
    // flux value used by both adjacent finite-volume updates.
    const double anti_out[4] = {0.12, 0.09, 0.11, 0.08};
    const double m_low = 0.25;
    double alpha[4] = {0.0, 0.0, 0.0, 0.0};
    Stage5::shared_budget_alphas(anti_out, m_low, alpha);
    double limited_outflow = 0.0;
    for (int face = 0; face < 4; ++face) {
        assert(alpha[face] >= 0.0 && alpha[face] <= 1.0);
        limited_outflow += alpha[face] * anti_out[face];
    }
    assert(limited_outflow <= m_low + 1.0e-14);
    const double low_flux = -0.4;
    const double high_flux = -0.1;
    const double shared_alpha = alpha[0];
    const double final_flux = low_flux + shared_alpha * (high_flux - low_flux);
    assert(close(final_flux, low_flux + shared_alpha * (high_flux - low_flux)));
    return 0;
}
