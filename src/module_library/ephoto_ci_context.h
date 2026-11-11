#ifndef EPHOTO_CI_CONTEXT_H
#define EPHOTO_CI_CONTEXT_H

#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include "ePhoto_assim.h"

// Owned by one leaf conductance solve. Input files and biochemical settings
// remain fixed for this scope; no results are shared across leaves or timesteps.
struct ephoto_ci_context {
    ephoto_state biochemical_state;
    // May be disabled for a matched cold-start comparison.
    bool reuse_biochemical_state = true;
    double previous_ci = std::numeric_limits<double>::quiet_NaN();
    double previous_slope = std::numeric_limits<double>::quiet_NaN();

    template <typename Evaluate>
    ephoto_outputs evaluate(double temperature, double ppfd, double ci,
                            double experiment, Evaluate&& calculate) {
        // NaN cannot be used safely in an ordered cache key. Preserve the
        // evaluator's existing handling of invalid inputs instead.
        if (!std::isfinite(temperature) || !std::isfinite(ppfd) ||
            !std::isfinite(ci) || !std::isfinite(experiment)) {
            return calculate();
        }
        const auto key = std::make_tuple(temperature, ppfd, ci, experiment);
        const auto found = results.find(key);
        if (found != results.end()) return found->second;
        const auto value = calculate();
        // Bound storage even when a difficult conductance solve tries many
        // temperatures. Per-Ci-solve residual/failure tracking is independent.
        if (results.size() >= 512) results.clear();
        results.emplace(key, value);
        return value;
    }

private:
    std::map<std::tuple<double, double, double, double>, ephoto_outputs> results;
};

#endif
