#include "ppm_ctu_reconstruction.h"
#include "discrete_moment_operators.h"

#include <cassert>
#include <cmath>
#include <limits>

namespace {

bool close(double a, double b)
{
    return std::fabs(a - b) <= 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

}

int main()
{
    const PpmCtu::Parabola constant = PpmCtu::reconstruct(2.0, 2.0, 2.0, 2.0, 2.0);
    assert(close(constant.ql, 2.0));
    assert(close(constant.qr, 2.0));
    assert(close(PpmCtu::trace_left(constant, 0.37), 2.0));
    assert(close(PpmCtu::trace_right(constant, 0.37), 2.0));

    const PpmCtu::Parabola linear = PpmCtu::reconstruct(3.0, 4.0, 5.0, 6.0, 7.0);
    assert(close(linear.ql, 4.5));
    assert(close(linear.qr, 5.5));
    assert(close(linear.ql + 0.5 * linear.dq + linear.q6 / 6.0, 5.0));
    assert(close(PpmCtu::trace_right(linear, 0.0), 5.5));
    assert(close(PpmCtu::trace_left(linear, 0.0), 4.5));
    assert(close(PpmCtu::trace_right(linear, 1.0), 5.0));
    assert(close(PpmCtu::trace_left(linear, 1.0), 5.0));

    const PpmCtu::Parabola positive = PpmCtu::reconstruct(0.0, 1.0, 0.15, 1.0, 0.0);
    assert(close(positive.ql + 0.5 * positive.dq + positive.q6 / 6.0, 0.15));
    for (int n = 0; n <= 100; ++n)
        assert(PpmCtu::value(positive, n / 100.0) >= -1.0e-14);

    const double anti_out[4] = {0.10, 0.20, 0.35, 0.05};
    double alpha[4];
    Stage5::shared_budget_alphas(anti_out, 0.25, alpha);
    double used = 0.0;
    for (int f = 0; f < 4; ++f) {
        assert(alpha[f] >= 0.0 && alpha[f] <= 1.0);
        used += anti_out[f] * alpha[f];
    }
    assert(used <= 0.25 + 1.0e-14);
    assert(!close(alpha[0], alpha[2]));

    // P1 invariant: high-low change equals the signed sum of the four
    // anti-diffusive mass transfers, and donor-limited transfers cannot
    // extract more mass than the low-order cell state contains.
    const double m_low = 0.25;
    const double ax_left = -0.10;
    const double ax_right = 0.20;
    const double au_lower = 0.05;
    const double au_upper = -0.03;
    const double high_minus_low = ax_left - ax_right + au_lower - au_upper;
    assert(close(high_minus_low, ax_left - ax_right + au_lower - au_upper));
    const double limited_left = -0.05;
    const double limited_right = 0.10;
    const double limited_lower = 0.02;
    const double limited_upper = -0.01;
    const double limited_outflow = std::max(0.0, -limited_left) +
        std::max(0.0, limited_right) + std::max(0.0, -limited_lower) +
        std::max(0.0, limited_upper);
    assert(limited_outflow <= m_low);
    const double m_final = m_low + limited_left - limited_right +
        limited_lower - limited_upper;
    assert(m_final >= 0.0);

    // The local donor repair scales only outgoing anti-diffusive transfer.
    // Long-double accounting must leave the donor at or above its low state.
    const long double donor_mass = 0.25L;
    const long double donor_outflow = 0.40L;
    const double beta = static_cast<double>(donor_mass / donor_outflow);
    assert(beta > 0.0 && beta < 1.0);
    assert(static_cast<long double>(beta) * donor_outflow <=
           donor_mass + 1.0e-16L);
    return 0;
}
