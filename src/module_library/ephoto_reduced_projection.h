#ifndef EPHOTO_REDUCED_PROJECTION_H
#define EPHOTO_REDUCED_PROJECTION_H
#include "ephoto_reduced_solver.h"
#include "leaf_energy_balance.h"
#include "water_and_air_properties.h"
#include "conductance_helpers.h"
#include "ball_berry_gs.h"
#include "../framework/constants.h"

inline ephoto_reduced_projection make_ephoto_reduced_projection(double longwave,double shortwave,
    double pressure,double air_temperature,double canopy_gbw,double width,double rh,double wind,
    double ca,double adjusted_b0,double adjusted_b1)
{
    if(!std::isfinite(adjusted_b1) || adjusted_b1<=0. || !std::isfinite(adjusted_b0) ||
       adjusted_b0<0. || !std::isfinite(ca) || ca<=0.) return {};
    return [=](double gs) {
        using physical_constants::dr_boundary;
        using physical_constants::dr_stomata;
        const auto energy=leaf_energy_balance(longwave,shortwave,pressure,air_temperature,
            canopy_gbw,width,rh,gs,wind);
        const double temperature=air_temperature+energy.Deltat,gbw=energy.gbw_molar;
        if(!std::isfinite(gbw) || gbw<=0. || !std::isfinite(gs) || gs<=adjusted_b0)
            throw std::runtime_error("invalid reduced Ci-Gs conductance");
        const double hs=std::min(1.,(gs+gbw*rh*saturation_vapor_pressure(air_temperature)/
            saturation_vapor_pressure(temperature))/(gs+gbw));
        const double wanted=ca/(adjusted_b1*hs/(gs-adjusted_b0)+dr_boundary/gbw);
        const double gt=sequential_conductance({gbw/dr_boundary,gs/dr_stomata});
        const double ci=ca-wanted/gt;
        const auto identity=ball_berry_gs(wanted*1e-6,ca*1e-6,rh,adjusted_b0,
            adjusted_b1,gbw,temperature,air_temperature);
        if(!std::isfinite(identity.gsw) || std::abs(identity.gsw-gs)>1e-10)
            throw std::runtime_error("reduced Ball-Berry inverse identity failed");
        return std::array<double,2>{{ci,wanted}};
    };
}
#endif
