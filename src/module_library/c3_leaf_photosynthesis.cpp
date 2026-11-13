#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "../math/roots/onedim/dekker.h"
#include "../math/roots/onedim/fixed_point.h"
#include "../math/roots/onedim/illinois.h"
#include "c3_temperature_response.h"           // for c3_temperature_response_parameters
#include "c3photo.h"                           // for c3photoC
#include "leaf_energy_balance.h"               // for leaf_energy_balance
#include "c3_leaf_photosynthesis.h"

using BMLePhoto::c3_leaf_photosynthesis;

string_vector c3_leaf_photosynthesis::get_inputs()
{
    return {
        "absorbed_longwave",            // J / (m^2 leaf) / s
        "absorbed_ppfd",                // micromol / (m^2 leaf) / s
        "absorbed_shortwave",           // J / (m^2 leaf) / s
        "atmospheric_pressure",         // Pa
        "b0",                           // mol / m^2 / s
        "b1",                           // dimensionless
        "beta_PSII",                    // dimensionless (fraction of absorbed light that reaches photosystem II)
        "Catm",                         // micromol / mol
        "c3_model_type",                // 1: FvCB; 2: ePhotosynthesis
        "electrons_per_carboxylation",  // electron / carboxylation
        "electrons_per_oxygenation",    // electron / oxygenation
        "gbw_canopy",                   // m / s
        "gm_at_25",                     // mol / m^2 / s / Pa
        "gm_Ha",                        // J / mol
        "gm_Hd",                        // J / mol
        "gm_S",                         // J / K / mol
        "Gs_min",                       // mol / m^2 / s
        "Gstar_at_25",                  // micromol / mol
        "Gstar_Ea",                     // J / mol
        "height",                       // m
        "Jmax_at_25",                   // micromol / m^2 / s
        "Jmax_Ea",                      // J / mol
        "Kc_at_25",                     // micromol / mol
        "Kc_Ea",                        // J / mol
        "Ko_at_25",                     // mmol / mol
        "Ko_Ea",                        // J / mol
        "leafwidth",                    // m
        "O2",                           // mmol / mol
        "phi_PSII_0",                   // dimensionless
        "phi_PSII_1",                   // (degrees C)^(-1)
        "phi_PSII_2",                   // (degrees C)^(-2)
        "rh",                           // dimensionless
        "RL_at_25",                     // micromol / m^2 / s
        "RL_Ea",                        // J / mol
        "StomataWS",                    // dimensionless
        "temp",                         // degrees C
        "theta_0",                      // dimensionless
        "theta_1",                      // (degrees C)^(-1)
        "theta_2",                      // (degrees C)^(-2)
        "Tp_at_25",                     // micromol / m^2 / s
        "Tp_Ha",                        // J / mol
        "Tp_Hd",                        // J / mol
        "Tp_S",                         // J / K / mol
        "Vcmax_at_25",                  // micromol / m^2 / s
        "Vcmax_Ea",                     // J / mol
        "windspeed",                    // m / s
        "exp_id"                        // dimensionless
    };
}

string_vector c3_leaf_photosynthesis::get_outputs()
{
    return {
        "Assim",             // micromol / m^2 /s
        "Cc",                // micromol / mol
        "Ci",                // micromol / mol
        "Cs",                // micromol / m^2 / s
        "EPenman",           // mmol / m^2 / s
        "EPriestly",         // mmol / m^2 / s
        "gbw",               // mol / m^2 / s
        "GrossAssim",        // micromol / m^2 /s
        "Gs",                // mol / m^2 / s
        "leaf_temperature",  // degrees C
        "RHs",               // dimensionless from Pa / Pa
        "RH_canopy",         // dimensionless
        "RL",                // micromol / m^2 / s
        "Rp",                // micromol / m^2 / s
        "TransR",            // mmol / m^2 / s
        "iteration_C3_Gs",   // not a physical quantity
        "residual_C3_Gs"     // mol / m^2 / s
    };
}

void c3_leaf_photosynthesis::do_operation() const
{
    if (c3_model_type != 1.0 && c3_model_type != 2.0) {
        throw std::invalid_argument(
            "c3_model_type must be 1 (FvCB) or 2 (ePhotosynthesis)");
    }
    int const model_type = static_cast<int>(c3_model_type);

    // Combine temperature response parameters
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

    // Make an initial guess for boundary layer conductance
    double const gbw_guess{1.2};  // mol / m^2 / s

    // Get an initial estimate of stomatal conductance, assuming the leaf is at
    // air temperature
    double const initial_stomatal_conductance =
        c3photoC(
            tr_param, absorbed_ppfd, ambient_temperature, ambient_temperature,
            rh, gm_at_25, Gstar_at_25, Kc_at_25, Ko_at_25, Vcmax_at_25,
            Jmax_at_25, Tp_at_25, RL_at_25, b0, b1, Gs_min, Catm,
            atmospheric_pressure, O2, StomataWS, electrons_per_carboxylation,
            electrons_per_oxygenation, beta_PSII, gbw_guess, exp_id, 1)
            .Gs;  // mol / m^2 / s

    const bool active_ephoto=model_type==2 && absorbed_ppfd/.85>=5.0;
    const double gs_minimum=StomataWS*b0+Gs_min*(1.0-StomataWS);
    const auto configuration=ephoto_hint_configuration(tr_param, {exp_id,double(model_type),
        b0,b1,Gs_min,Catm,RL_at_25,gm_at_25,Gstar_at_25,Kc_at_25,Ko_at_25,
        Vcmax_at_25,Jmax_at_25,Tp_at_25,StomataWS,leafwidth});
    const bool use_hint=active_ephoto && ephoto_leaf_option("EPHOTO_WEATHER_HINTS") &&
        previous_leaf_hint.usable(configuration,ambient_temperature,absorbed_ppfd/.85,rh,gs_minimum);
    const bool known_unsteady=active_ephoto && previous_leaf_hint.unsteady &&
        previous_leaf_hint.nearby(configuration,ambient_temperature,absorbed_ppfd/.85,rh);
    photosynthesis_outputs photo;
    energy_balance_outputs et;
    ephoto_ci_context ci_context;
    if(use_hint) ci_context.previous_ci=previous_leaf_hint.ci;

    auto evaluate_gs = [=, &photo, &et, &ci_context](double current_gs) {
        // 2. Solve Energy Balance with current g_s
        et = leaf_energy_balance(
            absorbed_longwave,
            absorbed_shortwave,
            atmospheric_pressure,
            ambient_temperature,
            gbw_canopy,
            leafwidth,
            rh,
            current_gs,
            windspeed);

        double current_Tleaf = ambient_temperature + et.Deltat;  // degrees C

        // 3. Recalculate g_s with current Tleaf
        photo =
            c3photoC(
                tr_param, absorbed_ppfd, current_Tleaf, ambient_temperature,
                rh, gm_at_25, Gstar_at_25, Kc_at_25, Ko_at_25, Vcmax_at_25,
                Jmax_at_25, Tp_at_25, RL_at_25, b0, b1, Gs_min, Catm,
                atmospheric_pressure, O2, StomataWS, electrons_per_carboxylation,
                electrons_per_oxygenation, beta_PSII, et.gbw_molar, exp_id,
                model_type, &ci_context);

        return photo.Gs;
    };

    const auto result=solve_ephoto_leaf_gs(evaluate_gs,ci_context,initial_stomatal_conductance,
        gs_minimum,use_hint ? previous_leaf_hint.gs : std::numeric_limits<double>::quiet_NaN(),
        "c3_leaf_photosynthesis",active_ephoto,!known_unsteady);
    if(active_ephoto) {
        previous_leaf_hint={NAN,NAN,ambient_temperature,absorbed_ppfd/.85,rh,configuration,photo.penalty>=.01};
        if(photo.penalty<.01 && std::isfinite(photo.Ci) && std::abs(photo.residual)<=ephoto_ci_abs_tolerance) {
            previous_leaf_hint.ci=photo.Ci;
            previous_leaf_hint.gs=result.root;
        }
    } else {
        previous_leaf_hint=ephoto_leaf_hint{};
    }

    // Update the outputs
    update(Assim_op, photo.Assim);
    update(Cc_op, photo.Cc);
    update(Ci_op, photo.Ci);
    update(Cs_op, photo.Cs);
    update(EPenman_op, et.EPenman);
    update(EPriestly_op, et.EPriestly);
    update(gbw_op, et.gbw_molar);
    update(GrossAssim_op, photo.GrossAssim);
    update(Gs_op, photo.Gs);
    update(leaf_temperature_op, ambient_temperature + et.Deltat);
    update(RHs_op, photo.RHs);
    update(RH_canopy_op, et.RH_canopy);
    update(RL_op, photo.RL);
    update(Rp_op, photo.Rp);
    update(TransR_op, et.TransR);
    update(iteration_C3_Gs_op, result.iteration);
    update(residual_C3_Gs_op, result.residual);
}
