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
#include "ephoto_joint_solver.h"
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
// Accuracy, cold-start generation and local-recovery mode separate trials.
template <typename Evaluate, typename Energy, typename Photo>
auto memoize_ephoto_leaf_trial(Evaluate&& evaluate, ephoto_ci_context& context, Energy& energy, Photo& photo, bool use_ephoto)
{
    using Key = std::tuple<double, double, double, std::size_t, bool>;
    using Snapshot = std::tuple<Energy, Photo, double, double>;
    return [&, use_ephoto, evaluate = std::forward<Evaluate>(evaluate),
            trials = std::map<Key, Snapshot>{}](double gs) mutable {
        if (!use_ephoto || !std::isfinite(gs)) return evaluate(gs);
        const Key key{gs, context.ci_abs_tolerance, context.ci_rel_tolerance, context.cold_start_generation,
                      context.local_ci_refinement};
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

// Exceptional local recovery preserves ordinary caches and policies on failure.
template <typename Evaluate>
root_finding::result_t recover_local_ephoto_leaf(Evaluate& raw, ephoto_ci_context& context,
    double candidate, double minimum_gs, const std::string& label)
{
    const double anchor = context.previous_ci;
    const double slope = context.previous_slope;
    const bool reuse = context.reuse_biochemical_state;
    const bool guard = context.guard_acceleration;
    context.begin_local_ci_refinement(anchor);
    context.guard_acceleration = false;
    context.reuse_biochemical_state = false;
    std::cerr << "ePhoto local Ci recovery: gs=" << candidate << " anchor_Ci=" << anchor << '\n' << std::flush;
    try {
        std::string invalid_reason;
        auto residual = [&](double gs) {
            try { return raw(gs)-gs; }
            catch (const std::runtime_error& error) {
                invalid_reason=error.what();
                return std::numeric_limits<double>::quiet_NaN();
            }
        };
        const double window = std::max(.005, .02*std::abs(candidate));
        const auto bracket=[&]() {
            try {
                return bracket_ephoto_ci(residual, candidate,
                    std::max(minimum_gs, candidate-window), candidate+window,
                    std::max(1e-4, window/8), 1e-3, true, 8);
            } catch (const std::runtime_error& error) {
                throw std::runtime_error(label+": "+error.what()+
                    (invalid_reason.empty() ? std::string{} : ": "+invalid_reason));
            }
        }();
        root_finding::result_t result;
        if (bracket.lower==bracket.upper)
            result={bracket.lower,residual(bracket.lower),0,root_finding::Flag::residual_zero};
        else {
            root_finding::dekker solver{100,1e-4,1e-4};
            result=solver.solve(residual,bracket.lower,bracket.upper);
        }
        const double synchronized=raw(result.root);
        result.residual=synchronized-result.root;
        if (!root_finding::is_successful(result.flag) || !ephoto_gs_close(synchronized,result.root))
            throw std::runtime_error("local conductance recovery did not satisfy the final residual" +
                (invalid_reason.empty() ? std::string{} : ": "+invalid_reason));
        std::cerr << "ePhoto local Ci recovery accepted: residual=" << result.residual << '\n' << std::flush;
        return result;
    } catch (...) {
        context.local_ci_refinement=false;
        context.local_ci_anchor=std::numeric_limits<double>::quiet_NaN();
        context.previous_ci=anchor;
        context.previous_slope=slope;
        context.reuse_biochemical_state=reuse;
        context.guard_acceleration=guard;
        throw;
    }
}

// Internal control signal: only a failed strict verification can request the
// experimental early handoff. It never establishes convergence.
struct ephoto_early_joint : std::runtime_error {
    ephoto_early_joint() : std::runtime_error("strict cold conductance verification requested early joint recovery") {}
};

// Ordinary fixed-point/bracketed strategy and final acceptance tolerances.
// Coarse Ci work never establishes final convergence.
template <typename Evaluate>
root_finding::result_t solve_ephoto_leaf_gs(Evaluate&& raw, ephoto_ci_context& context, double initial_gs, double minimum_gs, double hinted_gs, const std::string& label, bool use_ephoto, bool allow_acceleration = true, bool allow_early_joint = true)
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
        context.last_gs_trial=gs;
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
    // The acceleration guard can reject a trial during final verification
    // or refinement as well as during the initial root search. All stages
    // must route that internal signal to the strict cold retry.
    try {
        bool local_recovery_attempted=false;
        std::string local_recovery_error;
        root_finding::result_t result;
        auto try_local_recovery = [&](double synchronized, bool final_check = false) {
            if (std::getenv("EPHOTO_JOINT_CI_GS") && ephoto_leaf_option("EPHOTO_JOINT_CI_GS")) return false;
            const double window=std::max(.005,.02*std::abs(result.root));
            // A correction outside this window needs the ordinary global
            // refinement first. Do not spend local ODE calls on that search.
            if (adaptive || use_hint || local_recovery_attempted ||
                !std::isfinite(context.previous_ci) || !std::isfinite(synchronized) ||
                std::abs(synchronized-result.root)>window ||
                !std::getenv("EPHOTO_LOCAL_CI_REFINEMENT") ||
                !ephoto_leaf_option("EPHOTO_LOCAL_CI_REFINEMENT")) return false;
            if (!final_check && !ephoto_gs_close(synchronized,result.root,5.0)) return false;
            local_recovery_attempted=true;
            try {
                result=recover_local_ephoto_leaf(raw,context,result.root,minimum_gs,label);
                return true;
            } catch (const std::runtime_error& error) {
                local_recovery_error=error.what();
                std::cerr << "ePhoto local Ci recovery failed: " << local_recovery_error << '\n' << std::flush;
                return false;
            }
        };
        try {
            result = fixed_then_bracketed(use_hint ? hinted_gs : initial_gs);
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
            // Keep the normal solve first. Once its fully strict, unhinted
            // candidate fails verification, avoid repeating global refinement
            // before trying the coupled residual function. Opt in explicitly.
            if (allow_early_joint && !adaptive && !use_hint &&
                std::getenv("EPHOTO_EARLY_JOINT") && ephoto_leaf_option("EPHOTO_EARLY_JOINT") &&
                std::getenv("EPHOTO_JOINT_CI_GS") && ephoto_leaf_option("EPHOTO_JOINT_CI_GS"))
                throw ephoto_early_joint{};
            if (try_local_recovery(synchronized)) return result;
            if (ephoto_runtime().enabled()) ++ephoto_runtime().strict_refinements;
            result = bracketed(result.root, std::pow(.1, refinement + 1));
        }
        // Ensure outputs and residuals refer to the last candidate, even on failure.
        const double synchronized = calculate(result.root);
        result.residual = synchronized - result.root;
        if (ephoto_gs_close(synchronized, result.root)) return result;
        if (try_local_recovery(synchronized,true)) return result;
        if (adaptive || use_hint) {
            // A failed final verification also gets the original strict cold solve.
            context.reset_for_cold_start();
            if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_cold_restarts;
            return solve_ephoto_leaf_gs(raw, context, initial_gs, minimum_gs,
                                        std::numeric_limits<double>::quiet_NaN(), label, use_ephoto, false, allow_early_joint);
        }
        std::ostringstream message;
        message << std::setprecision(17) << label << ": strict final conductance residual exceeds tolerance"
                << " gs_trial=" << result.root << " gs_calculated=" << synchronized
                << " residual=" << result.residual << " initial_gs=" << initial_gs
                << " minimum_gs=" << minimum_gs << " root_flag=" << root_finding::flag_message(result.flag);
        if (!local_recovery_error.empty()) message << " local_Ci_recovery=" << local_recovery_error;
        throw std::runtime_error(message.str());
    } catch (const ephoto_unsteady_acceleration&) {
        context.reset_for_cold_start();
        if (ephoto_runtime().enabled()) ++ephoto_runtime().leaf_cold_restarts;
        return solve_ephoto_leaf_gs(raw, context, initial_gs, minimum_gs,
                                    std::numeric_limits<double>::quiet_NaN(), label, use_ephoto, false, allow_early_joint);
    }
}

template<class Raw,class Energy,class Photo>
root_finding::result_t solve_joint_ephoto_leaf(Raw& raw,ephoto_ci_context& context,
    Energy& energy,Photo& photo,double ci,double gs,double minimum_gs,double maximum_ci)
{
    const bool reuse=context.reuse_biochemical_state, guard=context.guard_acceleration;
    context.begin_joint_evaluation();
    using Snapshot=std::tuple<Energy,Photo>;
    std::map<std::pair<double,double>,Snapshot> cache;
    auto evaluate=[&](double trial_ci,double trial_gs) {
        context.joint_ci=trial_ci;
        const auto key=std::make_pair(trial_ci,trial_gs);
        const auto found=cache.find(key);
        if(found!=cache.end()) {
            energy=std::get<0>(found->second); photo=std::get<1>(found->second);
        } else {
            raw(trial_gs);
            cache.emplace(key,Snapshot{energy,photo});
        }
        context.previous_ci=photo.Ci;
        context.last_gs_trial=trial_gs;
        return std::array<double,2>{{photo.residual,photo.Gs-trial_gs}};
    };
    auto& diagnostics=ephoto_runtime();
    const auto before=diagnostics.ode_completed;
    std::ostringstream started;
    started<<std::setprecision(17)<<"ePhoto joint Ci-Gs fallback started: search_version="<<ephoto_joint_search_version
        <<" initialization=cold Ci="<<ci<<" gs="<<gs<<'\n';
    std::cerr<<started.str()<<std::flush;
    try {
        const auto result=solve_ephoto_joint(evaluate,ci,gs,minimum_gs,maximum_ci);
        const auto verified=evaluate(result.ci,result.gs);
        if(!ephoto_joint_accepted(result.gs,verified))
            throw std::runtime_error("joint Ci-Gs final verification failed");
        std::ostringstream accepted;
        accepted<<std::setprecision(17)<<"ePhoto joint Ci-Gs accepted: Ci="<<result.ci
            <<" gs="<<result.gs<<" Ci_residual="<<verified[0]<<" Gs_residual="<<verified[1]
            <<" evaluations="<<result.evaluations;
        if(diagnostics.enabled()) accepted<<" ode_calls="<<diagnostics.ode_completed-before;
        accepted<<'\n';
        std::cerr<<accepted.str()<<std::flush;
        return {result.gs,verified[1],result.evaluations,root_finding::Flag::residual_zero};
    } catch(...) {
        context.joint_evaluation=false;
        context.joint_ci=NAN;
        context.reuse_biochemical_state=reuse; context.guard_acceleration=guard;
        throw;
    }
}

template<class Nested,class Raw,class Energy,class Photo>
root_finding::result_t solve_ephoto_leaf_with_joint(Nested& nested,Raw& raw,
    ephoto_ci_context& context,Energy& energy,Photo& photo,double initial_gs,
    double minimum_gs,double hinted_gs,double maximum_ci,const std::string& label,
    bool use_ephoto,bool allow_acceleration=true,bool allow_early_joint=true)
{
    const auto started=ephoto_runtime_diagnostics::clock::now();
    const auto ode_before=ephoto_runtime().ode_completed;
    std::size_t nested_calls=0;
    auto counted_nested=[&](double gs) {++nested_calls;return nested(gs);};
    auto report=[&](const char* trigger) {
        std::ostringstream message;
        message<<std::setprecision(10)<<"ePhoto joint handoff: trigger="<<trigger
            <<" nested_calls="<<nested_calls<<" nested_seconds="
            <<std::chrono::duration<double>(ephoto_runtime_diagnostics::clock::now()-started).count();
        if(ephoto_runtime().enabled()) message<<" nested_ode_calls="<<ephoto_runtime().ode_completed-ode_before;
        message<<'\n';
        std::cerr<<message.str()<<std::flush;
    };
    try {
        return solve_ephoto_leaf_gs(counted_nested,context,initial_gs,minimum_gs,hinted_gs,label,use_ephoto,allow_acceleration,allow_early_joint);
    } catch(const ephoto_early_joint&) {
        report("strict_verification");
        try {
            return solve_joint_ephoto_leaf(raw,context,energy,photo,photo.Ci,
                context.last_gs_trial,minimum_gs,maximum_ci);
        } catch(const std::runtime_error& error) {
            std::cerr<<"ePhoto early joint failed; resuming ordinary cold solve: "<<error.what()<<'\n'<<std::flush;
            // A failed speculative recovery cannot remove the original path.
            // A new generation also isolates all memoized ordinary snapshots.
            context.reset_for_cold_start();
            return solve_ephoto_leaf_with_joint(nested,raw,context,energy,photo,
                initial_gs,minimum_gs,NAN,maximum_ci,label,use_ephoto,false,false);
        }
    } catch(const std::runtime_error& error) {
        const std::string original=error.what();
        if(!use_ephoto || !std::getenv("EPHOTO_JOINT_CI_GS") || !ephoto_leaf_option("EPHOTO_JOINT_CI_GS") ||
           original.find("strict final conductance residual")==std::string::npos) throw;
        report("exhausted_refinement");
        try {
            return solve_joint_ephoto_leaf(raw,context,energy,photo,photo.Ci,
                context.last_gs_trial,minimum_gs,maximum_ci);
        } catch(const std::runtime_error& joint_error) {
            throw std::runtime_error(original+"; joint fallback failed: "+joint_error.what());
        }
    }
}

#endif
