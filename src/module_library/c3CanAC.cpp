#include <algorithm>
#include <cmath>
#include "../framework/constants.h"
#include "../math/quadrature/quad.h"
#include "../math/roots/onedim/dekker.h"
#include "../math/roots/onedim/fixed_point.h"
#include "../math/roots/onedim/illinois.h"
#include "c3photo.h"
#include "core/photosynthesis.h"
#include "leaf_energy_balance.h"
#include "core/atmosphere_light_scattering.h"
#include "respiration.h"
#include "c3CanAC.h"

canopy_photosynthesis_outputs c3CanAC(
    c3_temperature_response_parameters const tr_param,
    double const absorbed_longwave,
    double const ambient_temperature,
    double const atmospheric_pressure,
    double const atmospheric_scattering,
    double const atmospheric_transmittance,
    double const b0,
    double const b1,
    double const beta_PSII,
    double const Catm,
    double const chil,
    double const cosine_zenith_angle,
    double const electrons_per_carboxylation,
    double const electrons_per_oxygenation,
    double const exp_id,
    int const c3_model_type,
    double const gbw_canopy,
    double const gm_at_25,
    double const growth_respiration_fraction,
    double const Gs_min,
    double const Gstar_at_25,
    double const heightf,
    double const Jmax_at_25,
    double const k_diffuse,
    double const Kc_at_25,
    double const Ko_at_25,
    double const kpLN,
    double const LAI,
    double const leaf_reflectance_nir,
    double const leaf_reflectance_par,
    double const leaf_transmittance_nir,
    double const leaf_transmittance_par,
    double const leaf_width,
    double const leafN,
    double const lnb0,
    double const lnb1,
    double const o2,
    double const par_energy_content,
    double const par_energy_fraction,
    double const RH,
    double const RL_at_25,
    double const solarR,
    double const StomataWS,
    double const Tp_at_25,
    double Vcmax_at_25,
    double const WindSpeed,
    double const WindSpeedHeight,
    int const lnfun,
    int const nlayers,
    ephoto_canopy_history* history)
{
    core::atmosphere_light_scattering const light_model(
        cosine_zenith_angle,
        atmospheric_pressure,
        atmospheric_transmittance,
        atmospheric_scattering);

    core::canopy_light::parameters params = {
        chil,
        cosine_zenith_angle,
        heightf,
        k_diffuse,
        LAI,
        leaf_reflectance_nir,
        leaf_reflectance_par,
        leaf_transmittance_nir,
        leaf_transmittance_par,
        par_energy_content,
        par_energy_fraction};
    core::canopy_light light_dist =
        core::canopy_light::from_solar(solarR, light_model, params);

    const auto configuration=ephoto_hint_configuration(tr_param, {exp_id,double(c3_model_type),
        b0,b1,Gs_min,Catm,RL_at_25,gm_at_25,Gstar_at_25,Kc_at_25,Ko_at_25,
        Vcmax_at_25,Jmax_at_25,Tp_at_25,StomataWS,leaf_width,kpLN,leafN,double(lnfun)});
    const bool history_compatible=history && history->layers==nlayers &&
        std::abs(history->lai-LAI)<=std::max(.1,.2*std::abs(history->lai));
    std::vector<ephoto_leaf_hint> next_hints;
    std::size_t hint_slot=0;

    auto leaf_photo =
        [&](double iabs,
            double j_shortwave,
            double layer_wind_speed,
            double layer_leafN) -> core::leaf_assim {
        // Fixed Gauss-Legendre traversal visits each quadrature point's
        // sunlit then shaded leaf in the same order across weather rows.
        const std::size_t slot=hint_slot++;
        next_hints.emplace_back();
        const bool active_ephoto=c3_model_type==2 && iabs/.85>=5.0;
        const double gs_minimum=StomataWS*b0+Gs_min*(1.0-StomataWS);
        const ephoto_leaf_hint* hint=nullptr;
        if(active_ephoto && ephoto_leaf_option("EPHOTO_WEATHER_HINTS") && history_compatible &&
           slot<history->leaves.size() && history->leaves[slot].usable(
               configuration,ambient_temperature,iabs/.85,RH,gs_minimum))
            hint=&history->leaves[slot];
        const bool known_unsteady=active_ephoto && history_compatible && slot<history->leaves.size() &&
            history->leaves[slot].unsteady && history->leaves[slot].nearby(
                configuration,ambient_temperature,iabs/.85,RH);
        double const effective_Vcmax =
            (lnfun != 0) ? layer_leafN * lnb1 + lnb0 : Vcmax_at_25;
        double constexpr gbw_guess = 1.2;

        // Use FvCB for the initial conductance guess, then use the selected
        // biochemical model for the coupled iteration.
        double const gsw_estimate =
            c3photoC(
                tr_param,
                iabs,
                ambient_temperature,
                ambient_temperature,
                RH,
                gm_at_25,
                Gstar_at_25,
                Kc_at_25,
                Ko_at_25,
                effective_Vcmax,
                Jmax_at_25,
                Tp_at_25,
                RL_at_25,
                b0,
                b1,
                Gs_min,
                Catm,
                atmospheric_pressure,
                o2,
                StomataWS,
                electrons_per_carboxylation,
                electrons_per_oxygenation,
                beta_PSII,
                gbw_guess,
                exp_id,
                1)
                .Gs;

        energy_balance_outputs et;
        photosynthesis_outputs photo;
        ephoto_ci_context ci_context;
        if(hint) ci_context.previous_ci=hint->ci;

        auto evaluate_gs = [&](double current_gs) {
            et = leaf_energy_balance(
                absorbed_longwave,
                j_shortwave,
                atmospheric_pressure,
                ambient_temperature,
                gbw_canopy,
                leaf_width,
                RH,
                current_gs,
                layer_wind_speed);

            double const leaf_temperature =
                ambient_temperature + et.Deltat;

            photo = c3photoC(
                tr_param,
                iabs,
                leaf_temperature,
                ambient_temperature,
                RH,
                gm_at_25,
                Gstar_at_25,
                Kc_at_25,
                Ko_at_25,
                effective_Vcmax,
                Jmax_at_25,
                Tp_at_25,
                RL_at_25,
                b0,
                b1,
                Gs_min,
                Catm,
                atmospheric_pressure,
                o2,
                StomataWS,
                electrons_per_carboxylation,
                electrons_per_oxygenation,
                beta_PSII,
                et.gbw_molar,
                exp_id,
                c3_model_type,
                &ci_context);

            return photo.Gs;
        };

        const auto result=solve_ephoto_leaf_gs(evaluate_gs,ci_context,gsw_estimate,
            gs_minimum,hint ? hint->gs : std::numeric_limits<double>::quiet_NaN(),
            "c3Canopy",active_ephoto,!known_unsteady);
        if(active_ephoto) {
            next_hints[slot]={NAN,NAN,ambient_temperature,iabs/.85,RH,configuration,photo.penalty>=.01};
            if(photo.penalty<.01 && std::isfinite(photo.Ci) && std::abs(photo.residual)<=ephoto_ci_abs_tolerance) {
                next_hints[slot].ci=photo.Ci;
                next_hints[slot].gs=result.root;
            }
        }

        double constexpr transpiration_conversion =
            physical_constants::molar_mass_of_water * 36;

        return core::leaf_assim{
            photo.Assim,
            photo.Gs,
            et.EPenman,
            et.EPriestly,
            photo.GrossAssim,
            photo.RL,
            photo.Rp,
            et.TransR * transpiration_conversion};
    };

    core::canopy_integrand integrand(
        leaf_photo,
        light_dist,
        kpLN,
        leafN,
        WindSpeed);

    core::leaf_assim const canopy =
        quadrature::gauss_legendre<2, core::leaf_assim>(
            integrand,
            0.0,
            LAI,
            nlayers);

    if(history) {
        history->leaves=std::move(next_hints);
        history->layers=nlayers;
        history->lai=LAI;
    }

    double const whole_plant_gr =
        growth_resp(canopy.assim, growth_respiration_fraction);

    return canopy_photosynthesis_outputs{
        canopy.assim - whole_plant_gr,
        canopy.stomatal_vapor_conductance,
        canopy.penman,
        canopy.priestly,
        canopy.carboxylation,
        canopy.leaf_respiration,
        canopy.photorespiration,
        canopy.transpiration,
        whole_plant_gr};
}
