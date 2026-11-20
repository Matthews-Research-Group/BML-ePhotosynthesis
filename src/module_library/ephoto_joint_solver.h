#ifndef EPHOTO_JOINT_SOLVER_H
#define EPHOTO_JOINT_SOLVER_H
#include <array>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <iostream>

constexpr int ephoto_joint_search_version=4;
constexpr std::size_t ephoto_joint_max_evaluations=200;

struct ephoto_joint_result {
    double ci, gs;
    std::array<double,2> residual;
    std::size_t evaluations;
};
inline bool ephoto_joint_accepted(double gs,const std::array<double,2>& r) {
    return std::isfinite(r[0]) && std::isfinite(r[1]) && std::abs(r[0])<=1e-3 &&
        std::abs(r[1])<=std::max(1e-3,1e-3*std::min(std::abs(gs),std::abs(gs+r[1])));
}

// Damped two-variable Newton search, with finite-difference Jacobians and
// regularized least-squares steps when the local Jacobian is ill-conditioned.
// Coordinates are Ci (ppm) and Gs/0.01 to keep step sizes comparable.
template<class Evaluate>
ephoto_joint_result solve_ephoto_joint_once(Evaluate&& evaluate,double ci,double gs,
    double minimum_gs,double maximum_ci,std::size_t maximum_evaluations=120)
{
    using Pair=std::array<double,2>;
    Pair x{{std::max(1e-6,std::min(ci,maximum_ci)),std::max(minimum_gs,gs)/.01}};
    const double upper_gs=std::max(1.,2*gs);
    std::size_t calls=0;
    auto trial=[&](Pair p) {
        if(calls>=maximum_evaluations) throw std::runtime_error("joint Ci-Gs evaluation limit reached");
        ++calls;
        try {
            auto r=evaluate(p[0],p[1]*.01);
            if(std::isfinite(r[0]) && std::isfinite(r[1])) return r;
        } catch(const std::runtime_error&) {}
        return Pair{{NAN,NAN}};
    };
    auto merit=[](Pair r) { return std::hypot(r[0],r[1]); };
    auto valid=[](Pair r) { return std::isfinite(r[0]) && std::isfinite(r[1]); };
    auto inside=[&](Pair p) { return p[0]>=1e-6 && p[0]<=maximum_ci &&
        p[1]*.01>=minimum_gs && p[1]*.01<=upper_gs; };
    Pair f=trial(x);
    if(!valid(f)) throw std::runtime_error("invalid initial joint Ci-Gs evaluation");
    for(std::size_t iteration=0;iteration<40;++iteration) {
        if(ephoto_joint_accepted(x[1]*.01,f)) return {x[0],x[1]*.01,f,calls};
        bool advanced=false;
        // Larger difference steps are useful when ODE error obscures a slope.
        for(double h : {.01,.1}) {
            double j[2][2]; bool jacobian=true;
            for(int axis=0;axis<2;++axis) {
                Pair probe=x; probe[axis]+=h;
                Pair fp=inside(probe) ? trial(probe) : Pair{{NAN,NAN}};
                double step=h;
                if(!valid(fp)) {
                    probe=x; probe[axis]-=h; step=-h;
                    fp=inside(probe) ? trial(probe) : Pair{{NAN,NAN}};
                }
                if(!valid(fp)) { jacobian=false; break; }
                for(int k=0;k<2;++k) j[k][axis]=(fp[k]-f[k])/step;
            }
            if(!jacobian) continue;
            for(double regularization : {0.,1e-6,.001,.1}) {
                Pair d;
                const double determinant=j[0][0]*j[1][1]-j[0][1]*j[1][0];
                if(regularization==0.) {
                    if(std::abs(determinant)<1e-14) continue;
                    d={{(-f[0]*j[1][1]+j[0][1]*f[1])/determinant,
                        (-j[0][0]*f[1]+f[0]*j[1][0])/determinant}};
                } else {
                    const double a=j[0][0]*j[0][0]+j[1][0]*j[1][0]+regularization;
                    const double b=j[0][0]*j[0][1]+j[1][0]*j[1][1];
                    const double c=j[0][1]*j[0][1]+j[1][1]*j[1][1]+regularization;
                    const double u=-(j[0][0]*f[0]+j[1][0]*f[1]);
                    const double v=-(j[0][1]*f[0]+j[1][1]*f[1]);
                    d={{(c*u-b*v)/(a*c-b*b),(a*v-b*u)/(a*c-b*b)}};
                }
                const double scale=std::max({1.,std::abs(d[0])/2.,std::abs(d[1])/2.});
                d[0]/=scale; d[1]/=scale;
                for(double alpha=1.;alpha>=1./64;alpha*=.5) {
                    Pair next{{x[0]+alpha*d[0],x[1]+alpha*d[1]}};
                    if(!inside(next)) continue;
                    Pair fn=trial(next);
                    if(valid(fn) && (ephoto_joint_accepted(next[1]*.01,fn) || merit(fn)<merit(f)*(1.-1e-4*alpha))) {
                        x=next; f=fn; advanced=true; break;
                    }
                }
                if(advanced) break;
            }
            if(advanced) break;
        }
        if(!advanced) throw std::runtime_error("joint Ci-Gs search could not reduce both residuals");
    }
    throw std::runtime_error("joint Ci-Gs iteration limit reached");
}

// Multiple nearby basins can contain distinct roots. Reserve work for new
// starting pairs instead of spending the entire budget at a stalled point.
template<class Evaluate>
ephoto_joint_result solve_ephoto_joint(Evaluate&& evaluate,double ci,double gs,
    double minimum_gs,double maximum_ci,std::size_t maximum_evaluations=ephoto_joint_max_evaluations)
{
    using Pair=std::array<double,2>;
    if(!std::isfinite(ci) || !std::isfinite(gs) || !std::isfinite(minimum_gs) ||
       !std::isfinite(maximum_ci) || minimum_gs<0. || maximum_ci<=1e-6 || maximum_evaluations==0)
        throw std::runtime_error("invalid joint Ci-Gs search inputs");
    ci=std::max(1e-6,std::min(ci,maximum_ci));
    gs=std::max(minimum_gs,gs);
    const double upper=std::max(1.,2*gs);
    std::size_t calls=0;
    Pair first{{NAN,NAN}},best_x{{ci,gs}},best_f{{NAN,NAN}};
    double best_merit=std::numeric_limits<double>::infinity();
    bool have_solution=false;
    ephoto_joint_result sampled{};
    auto shared=[&](double c,double g) {
        if(calls>=maximum_evaluations) throw std::runtime_error("joint total evaluation limit reached");
        ++calls;
        const auto f=evaluate(c,g);
        if(calls==1) first=f;
        const double score=std::hypot(f[0],f[1]);
        if(std::isfinite(score) && score<best_merit) {best_merit=score;best_x={{c,g}};best_f=f;}
        if(ephoto_joint_accepted(g,f) && !have_solution) {
            have_solution=true;sampled={c,g,f,calls};
        }
        return f;
    };
    std::string reason;
    auto attempt=[&](double c,double g,std::size_t budget,ephoto_joint_result& result) {
        try {
            result=solve_ephoto_joint_once(shared,c,g,minimum_gs,maximum_ci,budget);
            result.evaluations=calls;
            return true;
        } catch(const std::runtime_error& error) {reason=error.what();}
        if(have_solution) {result=sampled;result.evaluations=calls;return true;}
        return false;
    };
    ephoto_joint_result result;
    const std::size_t initial_budget=std::min<std::size_t>(120,maximum_evaluations);
    if(attempt(ci,gs,std::min<std::size_t>(40,initial_budget),result)) return result;
    const double direction=std::isfinite(first[1]) && first[1]<0. ? -1. : 1.;
    for(int restart=0;restart<4 && calls<maximum_evaluations;++restart) {
        const double sign=restart%2==0 ? direction : -direction;
        const double radius=restart<2 ? 1. : 2.;
        const double c=std::max(1e-6,std::min(ci+sign*radius*.25,maximum_ci));
        const double g=std::max(minimum_gs,std::min(gs+sign*radius*.01,upper));
        // Preserve the first three starts exactly. Wider starts only use the
        // extra budget after the original recovery sequence has failed.
        const std::size_t remaining=(restart<2 ? initial_budget : maximum_evaluations)-calls;
        const std::size_t budget=restart==0 ? std::max<std::size_t>(1,remaining/2) :
            restart==1 ? remaining : std::min<std::size_t>(40,remaining);
        std::ostringstream message;
        message<<std::setprecision(17)<<"ePhoto joint Ci-Gs restart: attempt="<<restart+2
            <<" Ci="<<c<<" gs="<<g<<" budget="<<budget<<" used="<<calls<<'\n';
        std::cerr<<message.str()<<std::flush;
        if(attempt(c,g,budget,result)) return result;
    }
    std::ostringstream message;
    message<<std::setprecision(17)<<"joint Ci-Gs recovery failed: evaluations="<<calls
        <<" best_Ci="<<best_x[0]<<" best_Gs="<<best_x[1]
        <<" Ci_residual="<<best_f[0]<<" Gs_residual="<<best_f[1]<<" reason="<<reason;
    std::cerr<<message.str()<<'\n'<<std::flush;
    throw std::runtime_error(message.str());
}
#endif
