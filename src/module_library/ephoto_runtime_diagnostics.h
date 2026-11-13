#ifndef EPHOTO_RUNTIME_DIAGNOSTICS_H
#define EPHOTO_RUNTIME_DIAGNOSTICS_H

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>

// Process-local counters only. No solver options, inputs, or results are changed.
// Set EPHOTO_RUNTIME_DIAGNOSTICS to the reporting interval in seconds; 0 disables.
struct ephoto_runtime_diagnostics {
    using clock = std::chrono::steady_clock;
    double interval = 0.0;
    clock::time_point started = clock::now(), last_report = started;
    bool announced = false;
    std::uint64_t ode_started = 0, ode_completed = 0, ode_failed = 0, slow_calls = 0;
    std::uint64_t cache_requests = 0, cache_hits = 0;
    std::uint64_t ci_started = 0, ci_completed = 0, ci_failed = 0;
    std::uint64_t leaf_started = 0, leaf_completed = 0;
    std::uint64_t seed_available = 0, seed_requested = 0, seed_used = 0, seed_saved = 0;
    std::uint64_t coarse_ci_solves=0, strict_ci_solves=0, strict_refinements=0;
    std::uint64_t weather_hints_used=0, leaf_cold_restarts=0, unsteady_acceleration_rejections=0;
    std::uint64_t internal_steps = 0, rhs_evaluations = 0;
    double ode_seconds = 0.0, ci_seconds = 0.0, max_ode_seconds = 0.0;
    double last_ode_seconds = 0.0, profile = 0.0, ci = 0.0, ppfd = 0.0, temperature = 0.0;
    double assimilation = 0.0, penalty = 0.0;
    double last_profile = 0.0, last_ci = 0.0, last_ppfd = 0.0, last_temperature = 0.0;
    long last_steps = 0, last_rhs = 0;

    ephoto_runtime_diagnostics() {
        if (const char* setting = std::getenv("EPHOTO_RUNTIME_DIAGNOSTICS")) {
            char* end = nullptr;
            const double parsed = std::strtod(setting, &end);
            if (end != setting && *end == '\0' && std::isfinite(parsed) && parsed > 0.0)
                interval = parsed;
        }
    }
    bool enabled() const { return interval > 0.0; }

    void report(bool force = false) {
        if (!enabled()) return;
        const auto now = clock::now();
        if (!announced) {
            std::cerr << "ePhoto runtime diagnostics: version=20261112 interval_seconds=" << interval
                      << " totals=cumulative ode_scope=adapter ci_scope=inclusive"
                      << " seed_policy=eligible_nearby end_time=5000\n" << std::flush;
            announced = true;
        }
        if (!force && std::chrono::duration<double>(now - last_report).count() < interval) return;
        last_report = now;
        std::ostringstream line;
        line << std::setprecision(10)
             << "ePhoto runtime summary: wall_seconds=" << std::chrono::duration<double>(now - started).count()
             << " ode_started=" << ode_started << " ode_completed=" << ode_completed
             << " ode_failed=" << ode_failed << " ode_seconds=" << ode_seconds
             << " max_ode_seconds=" << max_ode_seconds << " slow_calls_gt5s=" << slow_calls
             << " internal_steps=" << internal_steps << " rhs_evaluations=" << rhs_evaluations
             << " ci_started=" << ci_started << " ci_completed=" << ci_completed
             << " ci_failed=" << ci_failed << " ci_seconds=" << ci_seconds
             << " leaf_started=" << leaf_started << " leaf_completed=" << leaf_completed
             << " cache_requests=" << cache_requests << " cache_hits=" << cache_hits
             << " seed_available=" << seed_available << " seed_requested=" << seed_requested
             << " seed_used=" << seed_used << " seed_saved=" << seed_saved
             << " coarse_ci_solves=" << coarse_ci_solves << " strict_ci_solves=" << strict_ci_solves
             << " strict_refinements=" << strict_refinements << " weather_hints_used=" << weather_hints_used
             << " leaf_cold_restarts=" << leaf_cold_restarts
             << " unsteady_acceleration_rejections=" << unsteady_acceleration_rejections
             << " last_ode_seconds=" << last_ode_seconds << " last_steps=" << last_steps
             << " last_rhs=" << last_rhs << " last_A=" << assimilation << " last_penalty=" << penalty
             << " last_exp_id=" << last_profile << " last_Ci=" << last_ci
             << " last_PPFD=" << last_ppfd << " last_Tleaf=" << last_temperature
             << " active_exp_id=" << profile << " active_Ci=" << ci
             << " active_PPFD=" << ppfd << " active_Tleaf=" << temperature << '\n';
        std::cerr << line.str() << std::flush;
    }
};

inline ephoto_runtime_diagnostics& ephoto_runtime() {
    static ephoto_runtime_diagnostics diagnostics;
    return diagnostics;
}

// RAII counts failed/exceptional Ci solves as well as successful ones.
struct ephoto_ci_diagnostic_scope {
    ephoto_runtime_diagnostics& diagnostics = ephoto_runtime();
    ephoto_runtime_diagnostics::clock::time_point started;
    int exceptions = std::uncaught_exceptions();
    ephoto_ci_diagnostic_scope() {
        if (diagnostics.enabled()) {
            started = ephoto_runtime_diagnostics::clock::now();
            ++diagnostics.ci_started;
        }
    }
    ~ephoto_ci_diagnostic_scope() {
        if (diagnostics.enabled()) {
            diagnostics.ci_seconds += std::chrono::duration<double>(
                ephoto_runtime_diagnostics::clock::now() - started).count();
            ++diagnostics.ci_completed;
            if (std::uncaught_exceptions() > exceptions) ++diagnostics.ci_failed;
            diagnostics.report();
        }
    }
};

#endif
