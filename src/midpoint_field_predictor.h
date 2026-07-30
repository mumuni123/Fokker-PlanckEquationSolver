#ifndef MIDPOINT_FIELD_PREDICTOR_H
#define MIDPOINT_FIELD_PREDICTOR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// History-only initial guess for the nonlinear midpoint solve.  It never
// changes an accepted state or a convergence test: callers may use the
// proposal as the first endpoint guess and still solve the original fixed
// point to the original tolerances.
class MidpointFieldPredictor {
public:
    MidpointFieldPredictor()
        : enabled_(false), history_depth_(0), latest_dt_(0.0)
    {}

    void set_enabled(bool enabled)
    {
        if (enabled_ != enabled) clear();
        enabled_ = enabled;
    }

    bool enabled() const { return enabled_; }
    int history_depth() const { return history_depth_; }

    void clear()
    {
        previous_.clear();
        latest_.clear();
        history_depth_ = 0;
        latest_dt_ = 0.0;
    }

    bool propose(const std::vector<double>& field_n, double dt,
                 std::vector<double>& field_predicted)
    {
        field_predicted.clear();
        if (!enabled_ || history_depth_ < 2 || !(dt > 0.0) ||
            dt != latest_dt_ || field_n.size() != latest_.size() ||
            previous_.size() != latest_.size()) {
            if (enabled_ && history_depth_ >= 2 &&
                (!(dt > 0.0) || dt != latest_dt_ ||
                 field_n.size() != latest_.size())) {
                clear();
            }
            return false;
        }

        double scale = 0.0;
        for (size_t i = 0; i < latest_.size(); ++i) {
            if (!std::isfinite(field_n[i]) ||
                !std::isfinite(latest_[i]) ||
                !std::isfinite(previous_[i]) ||
                field_n[i] != latest_[i]) {
                clear();
                return false;
            }
            scale = std::max(scale, std::fabs(latest_[i]));
            scale = std::max(scale, std::fabs(previous_[i]));
        }
        scale = std::max(scale, std::numeric_limits<double>::min());

        const double ratio = std::max(
            0.0, std::min(1.25, dt / latest_dt_));
        field_predicted.resize(latest_.size());
        for (size_t i = 0; i < latest_.size(); ++i) {
            const double raw_increment =
                ratio * (latest_[i] - previous_[i]);
            const double increment_limit =
                0.5 * std::max(scale, std::fabs(latest_[i]));
            const double increment = std::max(
                -increment_limit,
                std::min(increment_limit, raw_increment));
            const double predicted = latest_[i] + increment;
            if (!std::isfinite(predicted)) {
                field_predicted.clear();
                clear();
                return false;
            }
            field_predicted[i] = predicted;
        }
        return true;
    }

    void commit_strict(const std::vector<double>& accepted_field,
                       size_t owned_size, double dt)
    {
        if (!enabled_) return;
        if (!(dt > 0.0) || accepted_field.size() < owned_size) {
            clear();
            return;
        }
        // Do not combine endpoints produced by different physical step
        // sizes.  After a dt change, collect two fresh strict endpoints.
        if (history_depth_ > 0 && dt != latest_dt_) clear();
        std::vector<double> next(owned_size, 0.0);
        for (size_t i = 0; i < owned_size; ++i) {
            if (!std::isfinite(accepted_field[i])) {
                clear();
                return;
            }
            next[i] = accepted_field[i];
        }
        if (history_depth_ > 0 && latest_.size() != owned_size) clear();
        if (history_depth_ > 0) previous_.swap(latest_);
        latest_.swap(next);
        history_depth_ = std::min(2, history_depth_ + 1);
        latest_dt_ = dt;
    }

private:
    bool enabled_;
    int history_depth_;
    double latest_dt_;
    std::vector<double> previous_;
    std::vector<double> latest_;
};

#endif
