#include "midpoint_field_predictor.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {
bool close(double a, double b)
{
    return std::fabs(a - b) <=
        16.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}
}

int main()
{
    MidpointFieldPredictor predictor;
    std::vector<double> predicted;
    const std::vector<double> e1 = {2.0, -4.0};
    const std::vector<double> e2 = {3.0, -3.0};

    predictor.commit_strict(e1, e1.size(), 0.5);
    if (predictor.propose(e1, 0.5, predicted)) return 1;

    predictor.set_enabled(true);
    predictor.commit_strict(e1, e1.size(), 0.5);
    if (predictor.propose(e1, 0.5, predicted)) return 2;
    predictor.commit_strict(e2, e2.size(), 0.5);
    if (!predictor.propose(e2, 0.5, predicted)) return 3;
    if (predicted.size() != 2 ||
        !close(predicted[0], 4.0) || !close(predicted[1], -2.0)) return 4;

    // A changed step invalidates the two-level history.
    if (predictor.propose(e2, 0.25, predicted)) return 5;
    if (predictor.history_depth() != 0) return 6;

    // A dt change while only one endpoint is available must also restart
    // history instead of pairing endpoints from unlike step sizes.
    predictor.commit_strict(e1, e1.size(), 0.5);
    predictor.commit_strict(e2, e2.size(), 0.25);
    if (predictor.history_depth() != 1) return 11;
    if (predictor.propose(e2, 0.25, predicted)) return 12;
    predictor.clear();

    // The pointwise extrapolation is bounded by half the local field scale.
    const std::vector<double> e3 = {0.0, 0.0};
    const std::vector<double> e4 = {10.0, 0.0};
    predictor.commit_strict(e3, e3.size(), 0.5);
    predictor.commit_strict(e4, e4.size(), 0.5);
    if (!predictor.propose(e4, 0.5, predicted)) return 7;
    if (!close(predicted[0], 15.0) || !close(predicted[1], 0.0)) return 8;

    // A state discontinuity (including an external restart/load) invalidates
    // history instead of applying a stale proposal.
    const std::vector<double> different = {9.0, 0.0};
    if (predictor.propose(different, 0.5, predicted)) return 9;
    if (predictor.history_depth() != 0) return 10;

    std::cout << "midpoint_field_predictor_test: PASS\n";
    return 0;
}
