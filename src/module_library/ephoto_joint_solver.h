#ifndef EPHOTO_JOINT_SOLVER_H
#define EPHOTO_JOINT_SOLVER_H
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
ephoto_joint_result solve_ephoto_joint(Evaluate&& evaluate,double ci,double gs,
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
#endif
