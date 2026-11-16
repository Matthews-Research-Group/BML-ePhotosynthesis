#ifndef EPHOTO_CI_CONTEXT_H
#define EPHOTO_CI_CONTEXT_H

#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include "ePhoto_assim.h"
#include "ephoto_runtime_diagnostics.h"
#include "ephoto_ci_bracket.h"

struct ephoto_unsteady_acceleration {};

// Owned by one leaf conductance solve. Input files and biochemical settings
// remain fixed for this scope; no results are shared across leaves or timesteps.
struct ephoto_ci_context {
    ephoto_state biochemical_state;
    bool guard_acceleration = false;
    std::size_t cold_start_generation = 0;
    double ci_abs_tolerance = ephoto_ci_abs_tolerance;
    double ci_rel_tolerance = ephoto_ci_rel_tolerance;
    // May be disabled for a matched cold-start comparison.
    bool reuse_biochemical_state = true;
    double previous_ci = std::numeric_limits<double>::quiet_NaN();
    double previous_slope = std::numeric_limits<double>::quiet_NaN();
    bool local_ci_refinement = false;
    double local_ci_anchor = std::numeric_limits<double>::quiet_NaN();
    bool joint_evaluation = false;
    double joint_ci = std::numeric_limits<double>::quiet_NaN();
    double last_gs_trial = std::numeric_limits<double>::quiet_NaN();

    void begin_joint_evaluation() {
        // Ordinary results may have used eligible biochemical seeds. Mixing
        // those transient endpoints with cold candidates changes the residual
        // function. Joint candidates must all use the same initialization.
        local_results.clear();
        joint_evaluation=true;
        guard_acceleration=false;
        reuse_biochemical_state=false;
    }

    void begin_local_ci_refinement(double anchor) {
        local_results=results;
        local_ci_anchor=anchor;
        local_ci_refinement=true;
    }

    void reset_for_cold_start() {
        ++cold_start_generation;
        results.clear();
        local_results.clear();
        biochemical_state=ephoto_state{};
        previous_ci=std::numeric_limits<double>::quiet_NaN();
        previous_slope=std::numeric_limits<double>::quiet_NaN();
        guard_acceleration=false;
        local_ci_refinement=false;
        local_ci_anchor=std::numeric_limits<double>::quiet_NaN();
        joint_evaluation=false;
        joint_ci=std::numeric_limits<double>::quiet_NaN();
    }

    ~ephoto_ci_context() {
        if (diagnostics_used) ++ephoto_runtime().leaf_completed;
    }

    template <typename Evaluate>
    ephoto_outputs evaluate(double temperature, double ppfd, double ci,
                            double experiment, Evaluate&& calculate) {
        auto& diagnostics = ephoto_runtime();
        if (diagnostics.enabled()) {
            if (!diagnostics_used) {
                diagnostics_used = true;
                ++diagnostics.leaf_started;
            }
            ++diagnostics.cache_requests;
        }
        // NaN cannot be used safely in an ordered cache key. Preserve the
        // evaluator's existing handling of invalid inputs instead.
        if (!std::isfinite(temperature) || !std::isfinite(ppfd) ||
            !std::isfinite(ci) || !std::isfinite(experiment)) {
            return calculate();
        }
        const auto key = std::make_tuple(temperature, ppfd, ci, experiment);
        auto& cache = (local_ci_refinement || joint_evaluation) ? local_results : results;
        const auto found = cache.find(key);
        if (found != cache.end()) {
            if (diagnostics.enabled()) ++diagnostics.cache_hits;
            return found->second;
        }
        const auto value = calculate();
        // Bound storage even when a difficult conductance solve tries many
        // temperatures. Per-Ci-solve residual/failure tracking is independent.
        if (cache.size() >= 512) cache.clear();
        cache.emplace(key, value);
        return value;
    }

private:
    bool diagnostics_used = false;
    std::map<std::tuple<double, double, double, double>, ephoto_outputs> results;
    std::map<std::tuple<double, double, double, double>, ephoto_outputs> local_results;
};

#endif
