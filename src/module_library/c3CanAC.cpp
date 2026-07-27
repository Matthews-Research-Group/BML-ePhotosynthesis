#include "../framework/constants.h"
#include "../math/quadrature/quad.h"
#include "../math/roots/onedim/fixed_point.h"
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
    int const nlayers)
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

    root_finding::fixed_point solver(50, 1e-3, 1e-3);

    auto leaf_photo =
        [&](double iabs,
            double j_shortwave,
            double layer_wind_speed,
            double layer_leafN) -> core::leaf_assim {
        double const effective_Vcmax =
            (lnfun != 0) ? layer_leafN * lnb1 + lnb0 : Vcmax_at_25;
        double constexpr gbw_guess = 1.2;

        // As in BML's c3_leaf_photosynthesis module, use FvCB for the initial
        // conductance guess and ePhotosynthesis for the coupled iteration.
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

        auto gs_func = [&](double current_gs) {
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
                2);

            return photo.Gs;
        };

        root_finding::result_t const result =
            solver.solve(gs_func, gsw_estimate);

        if (!root_finding::is_successful(result.flag)) {
            throw std::runtime_error(
                "c3Canopy solver reports failed convergence. Termination flag:\n    " +
                root_finding::flag_message(result.flag));
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
