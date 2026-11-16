#ifndef EPHOTO_CI_LOCAL_H
#define EPHOTO_CI_LOCAL_H

#include "ephoto_ci_bracket.h"
#include "../math/roots/onedim/roots.h"

// Exceptional recovery only: search within 0.1% of the selected Ci (at least
// 0.05 ppm). This is a branch-local numerical window, not a physical bound.
template <typename Function>
root_finding::result_t refine_local_ephoto_ci(Function& residual, double anchor,
    double minimum, double maximum, double absolute_tolerance, double relative_tolerance)
{
    if (!std::isfinite(anchor) || anchor < minimum || anchor > maximum)
        throw std::runtime_error("invalid local Ci anchor");
    const double window = std::max(.05, .001 * std::abs(anchor));
    const double lower_limit = std::max(minimum, anchor - window);
    const double upper_limit = std::min(maximum, anchor + window);
    auto bracket = bracket_ephoto_ci(residual, anchor, lower_limit, upper_limit,
        std::min(window, std::max(1e-6, .0001 * std::abs(anchor))), absolute_tolerance, false, 8);
    double left = bracket.lower, right = bracket.upper;
    double fl = residual(left), fr = residual(right);
    if (left == right) return {left, fl, 0, root_finding::Flag::residual_zero};
    for (std::size_t iteration = 0; iteration < 100; ++iteration) {
        const bool choose_left = std::abs(fl) < std::abs(fr);
        const double best = choose_left ? left : right;
        const double value = choose_left ? fl : fr;
        const double width_tolerance = std::max(absolute_tolerance,
            relative_tolerance * std::min(std::abs(left), std::abs(right)));
        if (value == 0.0 || (iteration > 0 && right-left <= width_tolerance && std::abs(value) <= absolute_tolerance))
            return {best, value, iteration, root_finding::Flag::residual_zero};
        // Three additional decimal refinements are enough to diagnose an
        // unresolved local residual; do not chase a jump to machine precision.
        if (right-left <= width_tolerance*.001) break;
        const double middle = left + .5 * (right-left);
        if (middle == left || middle == right) break;
        const double fm = residual(middle);
        if (!std::isfinite(fm))
            throw std::runtime_error("invalid evaluation inside local Ci bracket");
        if (fm == 0.0) return {middle, fm, iteration+1, root_finding::Flag::residual_zero};
        if (root_finding::same_signs(fl, fm)) { left=middle; fl=fm; }
        else { right=middle; fr=fm; }
    }
    throw std::runtime_error("local Ci bracket did not satisfy the residual tolerance");
}
#endif
