#ifndef EPHOTO_CI_BRACKET_H
#define EPHOTO_CI_BRACKET_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "../math/roots/onedim/roots.h"

// Expand around a useful estimate instead of evaluating the near-zero bound
// unconditionally. A nonfinite residual marks an invalid ODE evaluation.
struct ephoto_ci_bracket {
    double lower;
    double upper;
};

template <typename Function>
ephoto_ci_bracket bracket_ephoto_ci(Function& residual, double initial,
                                     double minimum, double maximum,
                                     double initial_radius = 0.0) {
    if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum >= maximum ||
        !std::isfinite(initial) || !std::isfinite(initial_radius) || initial_radius < 0.0) {
        throw std::runtime_error("invalid Ci search interval");
    }
    initial = std::max(minimum, std::min(initial, maximum));
    root_finding::graph_t lower{initial, residual(initial)};
    if (!std::isfinite(lower.y)) {
        throw std::runtime_error("initial Ci evaluation failed");
    }
    auto upper = lower;
    double lowerLimit = minimum;
    double upperLimit = maximum;
    bool lowerBlocked = false;
    bool upperBlocked = false;
    double radius = initial_radius > 0.0 ? initial_radius :
        std::max(0.05 * std::abs(initial), 1e-6);
    const bool lowerFirst = lower.y > 0.0;

    // If a trial fails, approach it from the valid side. Never construct a
    // bracket across a known failed evaluation or invent its residual's sign.
    auto expand = [&](root_finding::graph_t& edge, double target,
                      double& limit, bool& blocked) {
        if (blocked || target == edge.x) return;
        for (int refinement = 0; refinement < 6; ++refinement) {
            const double value = residual(target);
            if (std::isfinite(value)) {
                edge = {target, value};
                return;
            }
            limit = target;
            target = edge.x + 0.5 * (target - edge.x);
            if (target == edge.x) break;
        }
        blocked = true;
    };
    for (int expansion = 0; expansion < 32; ++expansion) {
        for (int side = 0; side < 2; ++side) {
            if (std::abs(lower.y) <= 1e-8) return {lower.x, lower.x};
            if (std::abs(upper.y) <= 1e-8) return {upper.x, upper.x};
            if (root_finding::opposite_signs(lower.y, upper.y)) {
                return {lower.x, upper.x};
            }
            const bool moveLower = side == 0 ? lowerFirst : !lowerFirst;
            if (moveLower) {
                expand(lower, std::max(lowerLimit, initial - radius), lowerLimit, lowerBlocked);
            } else {
                expand(upper, std::min(upperLimit, initial + radius), upperLimit, upperBlocked);
            }
        }
        radius = std::min(radius * 2.0, maximum - minimum);
    }
    if (std::abs(lower.y) <= 1e-8) return {lower.x, lower.x};
    if (std::abs(upper.y) <= 1e-8) return {upper.x, upper.x};
    if (root_finding::opposite_signs(lower.y, upper.y)) return {lower.x, upper.x};
    throw std::runtime_error("no valid sign-changing Ci bracket found");
}

#endif
