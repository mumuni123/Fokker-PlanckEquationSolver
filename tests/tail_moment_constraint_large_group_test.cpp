#include "tail_moment_constraint.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::vector<double> moment_column(double upar, double uperp)
{
    const double gamma = std::sqrt(1.0 + upar * upar + uperp * uperp);
    std::vector<double> c(6, 0.0);
    c[0] = 1.0;
    c[1] = Const::me * Const::c * upar;
    c[2] = -Const::qe * Const::c * upar / gamma;
    c[3] = Const::me * Const::c * Const::c * (gamma - 1.0);
    c[4] = Const::me * Const::c * Const::c * upar * upar / gamma;
    c[5] = Const::me * Const::c * Const::c * uperp * uperp / gamma;
    return c;
}

}  // namespace

int main()
{
    const size_t support_count = 4096;
    const size_t max_support = 7;
    std::vector<std::vector<double> > columns;
    std::vector<double> weights;
    std::vector<double> original_weights;
    std::vector<double> reference(6, 0.0);
    columns.reserve(support_count);
    weights.reserve(support_count);

    for (size_t q = 0; q < support_count; ++q) {
        const double s = static_cast<double>(q) /
                         static_cast<double>(support_count - 1);
        // Production flux parcels repeat the same quadrature coordinates
        // across many faces.  Keep a long duplicate-heavy prefix so the
        // regression also rejects reducers that only work for uniformly
        // diverse support ordering.
        const bool repeated_prefix = q < 3 * support_count / 4;
        const double sr = repeated_prefix
                              ? 0.15 * static_cast<double>(q % 4)
                              : s;
        const double upar = 12.0 + 5.0 * sr + 0.07 * std::sin(31.0 * sr);
        const double uperp = 0.2 + 2.0 * sr + 0.03 * std::cos(17.0 * sr);
        const double weight = 1.0e8 * (1.0 + 0.25 * std::sin(13.0 * s));
        columns.push_back(moment_column(upar, uperp));
        weights.push_back(weight);
        original_weights.push_back(weight);
        for (size_t r = 0; r < reference.size(); ++r)
            reference[r] += weight * columns.back()[r];
    }

    const bool compressed = tail_compress_moment_supports(
        columns, weights, reference, max_support, 1.0e-10);
    size_t active = 0;
    double max_relative_residual = 0.0;
    bool nonnegative = true;
    for (size_t q = 0; q < weights.size(); ++q) {
        nonnegative = nonnegative && std::isfinite(weights[q]) &&
                      weights[q] >= 0.0;
        if (weights[q] > 0.0) ++active;
    }
    for (size_t r = 0; r < reference.size(); ++r) {
        double got = 0.0;
        for (size_t q = 0; q < weights.size(); ++q)
            got += weights[q] * columns[q][r];
        max_relative_residual = std::max(
            max_relative_residual,
            std::fabs(got - reference[r]) /
                std::max(1.0, std::fabs(reference[r])));
    }

    // The same shared reducer is used by the population controller with an
    // additional Xw constraint.  Exercise that seven-row path as well even
    // though the production stall originated in the six-row converter.
    std::vector<std::vector<double> > columns7 = columns;
    std::vector<double> reference7 = reference;
    reference7.push_back(0.0);
    for (size_t q = 0; q < support_count; ++q) {
        const double x = 1.0e-6 * static_cast<double>(q % 17);
        columns7[q].push_back(x);
        reference7[6] += original_weights[q] * x;
    }
    std::vector<double> weights7 = original_weights;
    const bool compressed7 = tail_compress_moment_supports(
        columns7, weights7, reference7, max_support, 1.0e-10);
    size_t active7 = 0;
    double max_relative_residual7 = 0.0;
    bool nonnegative7 = true;
    for (size_t q = 0; q < weights7.size(); ++q) {
        nonnegative7 = nonnegative7 && std::isfinite(weights7[q]) &&
                       weights7[q] >= 0.0;
        if (weights7[q] > 0.0) ++active7;
    }
    for (size_t r = 0; r < reference7.size(); ++r) {
        double got = 0.0;
        for (size_t q = 0; q < weights7.size(); ++q)
            got += weights7[q] * columns7[q][r];
        max_relative_residual7 = std::max(
            max_relative_residual7,
            std::fabs(got - reference7[r]) /
                std::max(1.0, std::fabs(reference7[r])));
    }

    // A completely repeated quadrature group has rank one.  It is a valid
    // conservative input and must not be rejected merely because the optional
    // full-row polish matrix is singular.
    std::vector<std::vector<double> > rank_deficient_columns(
        support_count, moment_column(14.0, 0.75));
    std::vector<double> rank_deficient_weights(support_count, 2.0e7);
    std::vector<double> rank_deficient_reference(6, 0.0);
    for (size_t q = 0; q < support_count; ++q)
        for (size_t r = 0; r < rank_deficient_reference.size(); ++r)
            rank_deficient_reference[r] += rank_deficient_weights[q] *
                                            rank_deficient_columns[q][r];
    const bool rank_deficient_compressed = tail_compress_moment_supports(
        rank_deficient_columns, rank_deficient_weights,
        rank_deficient_reference, max_support, 1.0e-10);
    size_t rank_deficient_active = 0;
    double rank_deficient_residual = 0.0;
    for (size_t q = 0; q < rank_deficient_weights.size(); ++q)
        if (rank_deficient_weights[q] > 0.0) ++rank_deficient_active;
    for (size_t r = 0; r < rank_deficient_reference.size(); ++r) {
        double got = 0.0;
        for (size_t q = 0; q < rank_deficient_weights.size(); ++q)
            got += rank_deficient_weights[q] * rank_deficient_columns[q][r];
        rank_deficient_residual = std::max(
            rank_deficient_residual,
            std::fabs(got - rank_deficient_reference[r]) /
                std::max(1.0, std::fabs(rank_deficient_reference[r])));
    }

    const bool pass = compressed && nonnegative && active <= max_support &&
                      max_relative_residual <= 1.0e-10 && compressed7 &&
                      nonnegative7 && active7 <= max_support &&
                      max_relative_residual7 <= 1.0e-10 &&
                      rank_deficient_compressed &&
                      rank_deficient_active <= max_support &&
                      rank_deficient_residual <= 1.0e-10;
    std::printf("compressed=%d\n", compressed ? 1 : 0);
    std::printf("nonnegative=%d\n", nonnegative ? 1 : 0);
    std::printf("active_supports=%zu\n", active);
    std::printf("max_relative_residual=%.17g\n", max_relative_residual);
    std::printf("seven_moment_compressed=%d\n", compressed7 ? 1 : 0);
    std::printf("seven_moment_nonnegative=%d\n", nonnegative7 ? 1 : 0);
    std::printf("seven_moment_active_supports=%zu\n", active7);
    std::printf("seven_moment_max_relative_residual=%.17g\n",
                max_relative_residual7);
    std::printf("rank_deficient_compressed=%d\n",
                rank_deficient_compressed ? 1 : 0);
    std::printf("rank_deficient_active_supports=%zu\n",
                rank_deficient_active);
    std::printf("rank_deficient_max_relative_residual=%.17g\n",
                rank_deficient_residual);
    std::printf("status=%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
