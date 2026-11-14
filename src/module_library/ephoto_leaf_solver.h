#ifndef EPHOTO_LEAF_SOLVER_H
#define EPHOTO_LEAF_SOLVER_H

#include <cstdlib>
#include <limits>
#include <iomanip>
#include <sstream>
#include <map>
#include <tuple>
#include <utility>
#include <string>
#include <vector>
#include "ephoto_ci_context.h"
#include "c3_temperature_response.h"
#include "ephoto_ci_bracket.h"
#include "../math/roots/onedim/fixed_point.h"
#include "../math/roots/onedim/dekker.h"
#include "../math/roots/onedim/illinois.h"

// Process options support matched validation; numerical histories are module-owned.
inline bool ephoto_leaf_option(const char* name)
{
    const char* value = std::getenv(name);
    return !value || std::string(value) != "0";
}
struct ephoto_leaf_hint {
    double ci = std::numeric_limits<double>::quiet_NaN();
    double gs = std::numeric_limits<double>::quiet_NaN();
    double temperature = 0.0, ppfd = 0.0, rh = 0.0;
    std::vector<double> configuration;
    bool unsteady = false;
    bool nearby(const std::vector<double>& config, double t, double light, double humidity) const
    {
        return configuration == config && ppfd >= 5.0 && light >= 5.0 &&
               std::isfinite(t) && std::isfinite(light) && std::abs(t - temperature) <= 5.0 &&
               std::abs(rh - humidity) <= .2 && std::abs(light - ppfd) <= std::max(20.0, .25 * std::abs(ppfd));
    }
    bool usable(const std::vector<double>& config, double t, double light, double humidity, double minimum_gs) const
    {
        return !unsteady && nearby(config, t, light, humidity) && std::isfinite(ci) && ci > 0.0 &&
               std::isfinite(gs) && gs >= minimum_gs && gs <= 5.0;
    }
};
struct ephoto_canopy_history {
    std::vector<ephoto_leaf_hint> leaves;
    int layers = 0;
    double lai = 0.0;
};

inline std::vector<double> ephoto_hint_configuration(
    const c3_temperature_response_parameters& tr, std::initializer_list<double> parameters)
{
    std::vector<double> result(parameters);
    result.insert(result.end(), {tr.gm_Ha, tr.gm_Hd, tr.gm_S, tr.Gstar_Ea, tr.Jmax_Ea,
                                 tr.Kc_Ea, tr.Ko_Ea, tr.phi_PSII_0, tr.phi_PSII_1, tr.phi_PSII_2, tr.RL_Ea,
                                 tr.theta_0, tr.theta_1, tr.theta_2, tr.Tp_Ha, tr.Tp_Hd, tr.Tp_S, tr.Vcmax_Ea});
    return result;
}

inline bool ephoto_gs_close(double a, double b, double factor = 1.0)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a - b) <= factor * std::max(1e-3, 1e-3 * std::min(std::abs(a), std::abs(b)));
}

// The outer root finder requires repeatable f(gs). Inner Ci hints and state
// seeds can otherwise select a different result when the same gs is revisited.
// Restore the complete leaf evaluation, including all side-effect outputs.
// Accuracy and cold-start generation separate coarse work and cold retries.
template <typename Evaluate, typename Energy, typename Photo>
auto memoize_ephoto_leaf_trial(Evaluate&& evaluate, ephoto_ci_context& context, Energy& energy, Photo& photo, bool use_ephoto)
{
    using Key = std::tuple<double, double, double, std::size_t>;
    using Snapshot = std::tuple<Energy, Photo, double, double>;
    return [&, use_ephoto, evaluate = std::forward<Evaluate>(evaluate),
            trials = std::map<Key, Snapshot>{}](double gs) mutable {
        if (!use_ephoto || !std::isfinite(gs)) return evaluate(gs);
        const Key key{gs, context.ci_abs_tolerance, context.ci_rel_tolerance, context.cold_start_generation};
        const auto found = trials.find(key);
        if (found != trials.end()) {
            energy = std::get<0>(found->second);
            photo = std::get<1>(found->second);
            context.previous_ci = std::get<2>(found->second);
            context.previous_slope = std::get<3>(found->second);
            if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_trial_cache_hits;
            return photo.Gs;
        }
        const double gs_calculated = evaluate(gs);
        if (std::isfinite(gs_calculated))
            trials.emplace(key, Snapshot{energy, photo, context.previous_ci, context.previous_slope});
        return gs_calculated;
    };
}

// Same fixed-point/bracketed strategy and final tolerances as the existing
// conductance solver. Coarse Ci work never establishes final convergence.
template <typename Evaluate>
root_finding::result_t solve_ephoto_leaf_gs(Evaluate&& raw, ephoto_ci_context& context, double initial_gs, double minimum_gs, double hinted_gs, const std::string& label, bool use_ephoto, bool allow_acceleration = true)
{
    const bool adaptive = allow_acceleration && use_ephoto && ephoto_leaf_option("EPHOTO_ADAPTIVE_CI");
    bool strict = !adaptive;
    context.guard_acceleration = allow_acceleration && use_ephoto &&
                                 (adaptive || std::isfinite(hinted_gs));
    auto require_strict = [&]() {
        strict = true;
        context.ci_abs_tolerance = ephoto_ci_abs_tolerance;
        context.ci_rel_tolerance = ephoto_ci_rel_tolerance;
    };
    auto calculate = [&](double gs) {
        context.ci_abs_tolerance = strict ? ephoto_ci_abs_tolerance : .05;
        context.ci_rel_tolerance = strict ? ephoto_ci_rel_tolerance : 1e-4;
        double next = raw(gs);
        if (!strict && ephoto_gs_close(gs, next, 5.0)) {
            require_strict();
            if (ephoto_runtime().enabled()) ++ephoto_runtime().strict_refinements;
            next = raw(gs);
        }
        return next;
    };
    auto bracketed = [&](double estimate, double interval_scale = 1.0) {
        require_strict();
        auto residual = [&](double gs) {
            const double value = calculate(gs) - gs;
            if (!std::isfinite(value)) throw std::runtime_error(label + ": nonfinite conductance residual");
            return value;
        };
        const double lower = std::max(1e-8, minimum_gs);
        double upper = std::max(1.0, 2.0 * std::max(estimate, lower));
        const double f_lower = residual(lower);
        double f_upper = residual(upper);
        for (int i = 0; root_finding::same_signs(f_lower, f_upper) && i < 10; ++i) {
            upper *= 2.0;
            f_upper = residual(upper);
        }
        if (root_finding::same_signs(f_lower, f_upper))
            throw std::runtime_error(label + ": could not bracket conductance");
        // A narrow bracket alone need not imply a small conductance residual
        // on a steep response. Refinement must tighten the search, rather than
        // repeating the same convergence thresholds. Bisection guarantees
        // interval progress on these exceptional refinement passes.
        if (interval_scale < 1.0) {
            double left = lower, right = upper, left_residual = f_lower;
            double best = std::abs(f_lower) < std::abs(f_upper) ? lower : upper;
            double best_residual = std::abs(f_lower) < std::abs(f_upper) ? f_lower : f_upper;
            const double tolerance = 1e-3 * interval_scale;
            for (std::size_t iteration = 0; iteration < 100; ++iteration) {
                if (std::abs(best_residual) <= tolerance)
                    return root_finding::result_t{best, best_residual, iteration, root_finding::Flag::residual_zero};
                const double middle = left + .5 * (right - left);
                if (middle == left || middle == right ||
                    right - left <= std::max(tolerance, tolerance * std::min(std::abs(left), std::abs(right))))
                    return root_finding::result_t{best, best_residual, iteration, root_finding::Flag::bracket_width_zero};
                const double value = residual(middle);
                if (std::abs(value) < std::abs(best_residual)) {
                    best = middle;
                    best_residual = value;
                }
                if (root_finding::same_signs(left_residual, value)) {
                    left = middle;
                    left_residual = value;
                } else
                    right = middle;
            }
            throw std::runtime_error(label + ": conductance refinement exhausted its iteration limit");
        }
        root_finding::dekker solver{100, 1e-3, 1e-3};
        auto result = solver.solve(residual, lower, upper);
        if (!root_finding::is_successful(result.flag)) {
            root_finding::illinois fallback{200, 1e-3, 1e-3};
            result = fallback.solve(residual, lower, upper);
        }
        if (!root_finding::is_successful(result.flag))
            throw std::runtime_error(label + ": conductance solve failed: " + root_finding::flag_message(result.flag));
        return result;
    };
    auto fixed_then_bracketed = [&](double estimate) {
        root_finding::fixed_point solver{10, 1e-3, 1e-3};
        auto result = solver.solve(calculate, estimate);
        require_strict();
        if (!root_finding::is_successful(result.flag)) result = bracketed(estimate);
        return result;
    };
    const bool use_hint = allow_acceleration && use_ephoto && std::isfinite(hinted_gs);
    if (use_hint && ephoto_runtime().enabled()) ++ephoto_runtime().weather_hints_used;
    root_finding::result_t result;
    try {
        result = fixed_then_bracketed(use_hint ? hinted_gs : initial_gs);
    } catch (const ephoto_unsteady_acceleration&) {
        context.reset_for_cold_start();
        if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_cold_restarts;
        return solve_ephoto_leaf_gs(raw, context, initial_gs, minimum_gs,
                                    std::numeric_limits<double>::quiet_NaN(), label, use_ephoto, false);
    } catch (const std::runtime_error&) {
        if (!adaptive && !use_hint) throw;
        // A tentative approximation or weather guess must not prevent the
        // original, fully accurate cold-start strategy from being attempted.
        require_strict();
        context.reset_for_cold_start();
        if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_cold_restarts;
        result = fixed_then_bracketed(initial_gs);
    }
    require_strict();
    for (int refinement = 0; refinement < 3; ++refinement) {
        const double root_residual = result.residual;
        const double synchronized = calculate(result.root);
        if (!use_ephoto) return result;
        if (!ephoto_gs_close(synchronized, result.root)) {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(17) << "ePhoto conductance verification: refinement=" << refinement
                       << " gs_trial=" << result.root << " gs_calculated=" << synchronized
                       << " solver_residual=" << root_residual << " synchronized_residual=" << synchronized - result.root
                       << " root_flag=" << root_finding::flag_message(result.flag) << "\n";
            std::cerr << diagnostic.str() << std::flush;
        }
        result.residual = synchronized - result.root;
        if (ephoto_gs_close(synchronized, result.root)) return result;
        if (ephoto_runtime().enabled()) ++ephoto_runtime().strict_refinements;
        result = bracketed(result.root, std::pow(.1, refinement + 1));
    }
    // Ensure outputs and residuals refer to the last candidate, even on failure.
    const double synchronized = calculate(result.root);
    result.residual = synchronized - result.root;
    if (ephoto_gs_close(synchronized, result.root)) return result;
    if (adaptive || use_hint) {
        // A failed final verification also gets the original strict cold solve.
        context.reset_for_cold_start();
        if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_cold_restarts;
        return solve_ephoto_leaf_gs(raw, context, initial_gs, minimum_gs,
                                    std::numeric_limits<double>::quiet_NaN(), label, use_ephoto, false);
    }
    std::ostringstream message;
    message << std::setprecision(17) << label << ": strict final conductance residual exceeds tolerance"
            << " gs_trial=" << result.root << " gs_calculated=" << synchronized
            << " residual=" << result.residual << " initial_gs=" << initial_gs
            << " minimum_gs=" << minimum_gs << " root_flag=" << root_finding::flag_message(result.flag);
    throw std::runtime_error(message.str());
}

#endif
