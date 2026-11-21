#ifndef EPHOTO_REDUCED_SOLVER_H
#define EPHOTO_REDUCED_SOLVER_H
#include <functional>
#include <map>
#include <set>
#include <vector>
#include "ephoto_joint_solver.h"

// For a trial Gs, return the Ci and net assimilation required by the positive
// Ball-Berry branch and diffusion. The caller supplies the physical projection.
using ephoto_reduced_projection=std::function<std::array<double,2>(double)>;
constexpr int ephoto_reduced_search_version=1;
struct ephoto_reduced_unsupported : std::runtime_error {
    ephoto_reduced_unsupported() : std::runtime_error("reduced Ci-Gs requires positive net assimilation") {}
};

// Scalar search along the Ball-Berry/diffusion manifold. A scalar sign change
// guides search only; both original actual residuals establish acceptance.
template<class Evaluate>
ephoto_joint_result solve_ephoto_reduced(Evaluate&& evaluate,
    const ephoto_reduced_projection& project,double initial_gs,double minimum_gs,
    double maximum_ci,std::size_t maximum_evaluations=40)
{
    if(!project || !std::isfinite(initial_gs) || !std::isfinite(minimum_gs) ||
       minimum_gs<0. || !std::isfinite(maximum_ci) || maximum_ci<=1e-6 || maximum_evaluations==0)
        throw std::runtime_error("invalid reduced Ci-Gs inputs");
    initial_gs=std::max(minimum_gs,initial_gs);
    const double lower=minimum_gs+1e-7, upper=std::max(1.,2*initial_gs);
    initial_gs=std::max(lower,std::min(initial_gs,upper));
    std::map<double,double> points;
    std::size_t calls=0;
    bool success=false;
    ephoto_joint_result result{};
    auto sample=[&](double gs) {
        gs=std::max(lower,std::min(gs,upper));
        const auto found=points.find(gs);
        if(found!=points.end()) return found->second;
        if(calls>=maximum_evaluations) throw std::runtime_error("reduced Ci-Gs evaluation limit reached");
        ++calls;
        double f=NAN;
        try {
            const auto required=project(gs);
            if(std::isfinite(required[0]) && required[0]>=1e-6 && required[0]<=maximum_ci &&
               std::isfinite(required[1]) && required[1]>0.) {
                const auto actual=evaluate(required[0],gs);
                if(std::isfinite(actual[2]) && actual[2]<=0.) throw ephoto_reduced_unsupported{};
                if(std::isfinite(actual[0]) && std::isfinite(actual[1]) && std::isfinite(actual[2])) {
                    const std::array<double,2> residual{{actual[0],actual[1]}};
                    f=actual[2]-required[1];
                    if(ephoto_joint_accepted(gs,residual)) {
                        success=true;result={required[0],gs,residual,calls};
                    }
                }
            }
        } catch(const ephoto_reduced_unsupported&) {throw;}
        catch(const std::runtime_error&) {}
        points.emplace(gs,f);
        return f;
    };
    sample(initial_gs);
    std::set<std::pair<double,double>> attempted;
    for(double radius : {.0025,.005,.01,.02,.04,.08,.16,.32,.64,upper-lower}) {
        if(success) return result;
        sample(initial_gs-radius);
        if(success) return result;
        sample(initial_gs+radius);
        if(success) return result;
        // Snapshot brackets before inserting refinement points. Invalid samples
        // remain in the ordered map and separate disconnected valid intervals.
        std::vector<std::pair<double,double>> brackets;
        for(auto left=points.begin(),right=std::next(left);right!=points.end();++left,++right)
            if(std::isfinite(left->second) && std::isfinite(right->second) &&
               std::signbit(left->second)!=std::signbit(right->second))
                brackets.emplace_back(left->first,right->first);
        for(const auto& bracket:brackets) {
            if(!attempted.insert(bracket).second) continue;
            double a=bracket.first,b=bracket.second,fa=points.at(a),fb=points.at(b);
            for(int iteration=0;iteration<24 && b-a>1e-9;++iteration) {
                const double fraction=fa/(fa-fb);
                const double next=a+(b-a)*(fraction>.1 && fraction<.9 ? fraction : .5);
                const double f=sample(next);
                if(success) return result;
                if(!std::isfinite(f)) break;
                if(std::signbit(fa)!=std::signbit(f)) {b=next;fb=f;}
                else {a=next;fa=f;}
            }
        }
    }
    if(success) return result;
    throw std::runtime_error("reduced Ci-Gs search did not satisfy both actual residuals");
}
#endif
