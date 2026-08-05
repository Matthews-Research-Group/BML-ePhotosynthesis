#include "c3_temperature_response.h"
#include "c3CanAC.h"
#include "c3_canopy.h"
#include <stdexcept>

using BMLePhoto::c3_canopy;

string_vector c3_canopy::get_inputs()
{
    return {
        "absorbed_longwave",
        "atmospheric_pressure",
        "atmospheric_scattering",
        "atmospheric_transmittance",
        "b0",
        "b1",
        "beta_PSII",
        "Catm",
        "c3_model_type",
        "chil",
        "cosine_zenith_angle",
        "electrons_per_carboxylation",
        "electrons_per_oxygenation",
        "exp_id",
        "gbw_canopy",
        "gm_at_25",
        "gm_Ha",
        "gm_Hd",
        "gm_S",
        "growth_respiration_fraction",
        "Gs_min",
        "Gstar_at_25",
        "Gstar_Ea",
        "heightf",
        "Jmax_at_25",
        "Jmax_Ea",
        "k_diffuse",
        "Kc_at_25",
        "Kc_Ea",
        "Ko_at_25",
        "Ko_Ea",
        "kpLN",
        "lai",
        "leaf_reflectance_nir",
        "leaf_reflectance_par",
        "leaf_transmittance_nir",
        "leaf_transmittance_par",
        "LeafN",
        "leafwidth",
        "lnb0",
        "lnb1",
        "lnfun",
        "nlayers",
        "O2",
        "par_energy_content",
        "par_energy_fraction",
        "phi_PSII_0",
        "phi_PSII_1",
        "phi_PSII_2",
        "rh",
        "RL_at_25",
        "RL_Ea",
        "solar",
        "StomataWS",
        "temp",
        "theta_0",
        "theta_1",
        "theta_2",
        "Tp_at_25",
        "Tp_Ha",
        "Tp_Hd",
        "Tp_S",
        "Vcmax_at_25",
        "Vcmax_Ea",
        "windspeed",
        "windspeed_height"};
}

string_vector c3_canopy::get_outputs()
{
    return {
        "canopy_assimilation_molar_flux",
        "canopy_conductance",
        "canopy_gross_assimilation_molar_flux",
        "canopy_non_photorespiratory_CO2_release_molar_flux",
        "canopy_photorespiration_molar_flux",
        "canopy_transpiration_rate",
        "whole_plant_growth_respiration_molar_flux"};
}

void c3_canopy::do_operation() const
{
    if (c3_model_type != 1.0 && c3_model_type != 2.0) {
        throw std::invalid_argument(
            "c3_model_type must be 1 (FvCB) or 2 (ePhotosynthesis)");
    }

    c3_temperature_response_parameters const tr_param{
        gm_Ha,
        gm_Hd,
        gm_S,
        Gstar_Ea,
        Jmax_Ea,
        Kc_Ea,
        Ko_Ea,
        phi_PSII_0,
        phi_PSII_1,
        phi_PSII_2,
        RL_Ea,
        theta_0,
        theta_1,
        theta_2,
        Tp_Ha,
        Tp_Hd,
        Tp_S,
        Vcmax_Ea};

    canopy_photosynthesis_outputs const result = c3CanAC(
        tr_param,
        absorbed_longwave,
        temp,
        atmospheric_pressure,
        atmospheric_scattering,
        atmospheric_transmittance,
        b0,
        b1,
        beta_PSII,
        Catm,
        chil,
        cosine_zenith_angle,
        electrons_per_carboxylation,
        electrons_per_oxygenation,
        exp_id,
        static_cast<int>(c3_model_type),
        gbw_canopy,
        gm_at_25,
        growth_respiration_fraction,
        Gs_min,
        Gstar_at_25,
        heightf,
        Jmax_at_25,
        k_diffuse,
        Kc_at_25,
        Ko_at_25,
        kpLN,
        lai,
        leaf_reflectance_nir,
        leaf_reflectance_par,
        leaf_transmittance_nir,
        leaf_transmittance_par,
        leafwidth,
        LeafN,
        lnb0,
        lnb1,
        O2,
        par_energy_content,
        par_energy_fraction,
        rh,
        RL_at_25,
        solar,
        StomataWS,
        Tp_at_25,
        Vcmax_at_25,
        windspeed,
        windspeed_height,
        static_cast<int>(lnfun),
        static_cast<int>(nlayers));

    update(canopy_assimilation_molar_flux_op, result.Assim);
    update(canopy_conductance_op, result.canopy_conductance);
    update(canopy_gross_assimilation_molar_flux_op, result.GrossAssim);
    update(canopy_non_photorespiratory_CO2_release_rate_op, result.RL);
    update(canopy_photorespiration_molar_flux_op, result.Rp);
    update(canopy_transpiration_rate_op, result.Trans);
    update(whole_plant_growth_respiration_molar_flux_op, result.whole_plant_gr);
}
