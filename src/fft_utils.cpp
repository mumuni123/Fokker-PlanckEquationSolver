#include "fft_utils.h"
#include "parameters.h"
#include <algorithm>
#include <cmath>

namespace {
bool is_power_of_two(size_t n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

size_t next_power_of_two(size_t n)
{
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void fft_radix2(std::vector<std::complex<double> >& a, bool inverse)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        const double angle = 2.0 * Const::pi / static_cast<double>(len)
                           * (inverse ? 1.0 : -1.0);
        const std::complex<double> wlen(std::cos(angle), std::sin(angle));
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (size_t j = 0; j < len / 2; ++j) {
                const std::complex<double> u = a[i + j];
                const std::complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (size_t i = 0; i < n; ++i) a[i] *= inv_n;
    }
}

void bluestein_forward(std::vector<std::complex<double> >& a)
{
    static int plan_n = 0;
    static size_t plan_m = 0;
    static std::vector<std::complex<double> > chirp;
    static std::vector<std::complex<double> > kernel_fft;
    static std::vector<std::complex<double> > work;

    const int n = static_cast<int>(a.size());
    if (n <= 1) return;

    if (plan_n != n) {
        plan_n = n;
        plan_m = next_power_of_two(static_cast<size_t>(2 * n - 1));
        chirp.assign(static_cast<size_t>(n), std::complex<double>(0.0, 0.0));
        kernel_fft.assign(plan_m, std::complex<double>(0.0, 0.0));
        work.assign(plan_m, std::complex<double>(0.0, 0.0));

        for (int i = 0; i < n; ++i) {
            const double ii = static_cast<double>(i);
            const double angle = Const::pi * ii * ii / static_cast<double>(n);
            chirp[static_cast<size_t>(i)] =
                std::complex<double>(std::cos(-angle), std::sin(-angle));
            const std::complex<double> kernel_value =
                std::conj(chirp[static_cast<size_t>(i)]);
            kernel_fft[static_cast<size_t>(i)] = kernel_value;
            if (i > 0) kernel_fft[plan_m - static_cast<size_t>(i)] = kernel_value;
        }
        fft_radix2(kernel_fft, false);
    }

    std::fill(work.begin(), work.end(), std::complex<double>(0.0, 0.0));
    for (int i = 0; i < n; ++i) {
        work[static_cast<size_t>(i)] =
            a[static_cast<size_t>(i)] * chirp[static_cast<size_t>(i)];
    }

    fft_radix2(work, false);
    for (size_t i = 0; i < plan_m; ++i) work[i] *= kernel_fft[i];
    fft_radix2(work, true);

    for (int i = 0; i < n; ++i) {
        a[static_cast<size_t>(i)] =
            work[static_cast<size_t>(i)] * chirp[static_cast<size_t>(i)];
    }
}
}

void fft_any(std::vector<std::complex<double> >& a, bool inverse)
{
    if (is_power_of_two(a.size())) {
        fft_radix2(a, inverse);
        return;
    }

    if (inverse) {
        for (size_t i = 0; i < a.size(); ++i) a[i] = std::conj(a[i]);
        bluestein_forward(a);
        const double inv_n = 1.0 / static_cast<double>(a.size());
        for (size_t i = 0; i < a.size(); ++i) a[i] = std::conj(a[i]) * inv_n;
    } else {
        bluestein_forward(a);
    }
}
