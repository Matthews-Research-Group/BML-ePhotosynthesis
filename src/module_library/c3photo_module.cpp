#include <stdexcept>
#include "c3_temperature_response.h"
#include "c3photo.h"
#include "c3photo_module.h"

using BMLePhoto::c3photo_module;

c3photo_module::c3photo_module(
    state_map const& input_quantities,
    state_map* output_quantities)
    : direct_module{},
      absorbed_ppfd{get_input(input_quantities, "absorbed_ppfd")},
      atmospheric_pressure{get_input(input_quantities, "atmospheric_pressure")},
      b0{get_input(input_quantities, "b0")},
      b1{get_input(input_quantities, "b1")},
      beta_PSII{get_input(input_quantities, "beta_PSII")},
      Catm{get_input(input_quantities, "Catm")},
      c3_model_type{get_input(input_quantities, "c3_model_type")},
      electrons_per_carboxylation{get_input(input_quantities, "electrons_per_carboxylation")},
      electrons_per_oxygenation{get_input(input_quantities, "electrons_per_oxygenation")},
      exp_id{get_input(input_quantities, "exp_id")},
      gbw{get_input(input_quantities, "gbw")},
      gm_at_25{get_input(input_quantities, "gm_at_25")},
      gm_Ha{get_input(input_quantities, "gm_Ha")},
      gm_Hd{get_input(input_quantities, "gm_Hd")},
      gm_S{get_input(input_quantities, "gm_S")},
      Gs_min{get_input(input_quantities, "Gs_min")},
      Gstar_at_25{get_input(input_quantities, "Gstar_at_25")},
      Gstar_Ea{get_input(input_quantities, "Gstar_Ea")},
      Jmax_at_25{get_input(input_quantities, "Jmax_at_25")},
      Jmax_Ea{get_input(input_quantities, "Jmax_Ea")},
      Kc_at_25{get_input(input_quantities, "Kc_at_25")},
      Kc_Ea{get_input(input_quantities, "Kc_Ea")},
      Ko_at_25{get_input(input_quantities, "Ko_at_25")},
      Ko_Ea{get_input(input_quantities, "Ko_Ea")},
      O2{get_input(input_quantities, "O2")},
      phi_PSII_0{get_input(input_quantities, "phi_PSII_0")},
      phi_PSII_1{get_input(input_quantities, "phi_PSII_1")},
      phi_PSII_2{get_input(input_quantities, "phi_PSII_2")},
      rh{get_input(input_quantities, "rh")},
      RL_at_25{get_input(input_quantities, "RL_at_25")},
      RL_Ea{get_input(input_quantities, "RL_Ea")},
      StomataWS{get_input(input_quantities, "StomataWS")},
      temp{get_input(input_quantities, "temp")},
      theta_0{get_input(input_quantities, "theta_0")},
      theta_1{get_input(input_quantities, "theta_1")},
      theta_2{get_input(input_quantities, "theta_2")},
      Tleaf{get_input(input_quantities, "Tleaf")},
      Tp_at_25{get_input(input_quantities, "Tp_at_25")},
      Tp_Ha{get_input(input_quantities, "Tp_Ha")},
      Tp_Hd{get_input(input_quantities, "Tp_Hd")},
      Tp_S{get_input(input_quantities, "Tp_S")},
      Vcmax_at_25{get_input(input_quantities, "Vcmax_at_25")},
      Vcmax_Ea{get_input(input_quantities, "Vcmax_Ea")},
      Assim_op{get_op(output_quantities, "Assim")},
      Assim_conductance_op{get_op(output_quantities, "Assim_conductance")},
      Cc_op{get_op(output_quantities, "Cc")},
      Ci_op{get_op(output_quantities, "Ci")},
      Cs_op{get_op(output_quantities, "Cs")},
      GrossAssim_op{get_op(output_quantities, "GrossAssim")},
      Gs_op{get_op(output_quantities, "Gs")},
      iteration_C3_Assim_op{get_op(output_quantities, "iteration_C3_Assim")},
      penalty_op{get_op(output_quantities, "penalty")},
      residual_C3_Assim_op{get_op(output_quantities, "residual_C3_Assim")},
      RHs_op{get_op(output_quantities, "RHs")},
      RL_op{get_op(output_quantities, "RL")},
      Rp_op{get_op(output_quantities, "Rp")}
{
}

string_vector c3photo_module::get_inputs()
{
    return {
        "absorbed_ppfd",
        "atmospheric_pressure",
        "b0",
        "b1",
        "beta_PSII",
        "Catm",
        "c3_model_type",
        "electrons_per_carboxylation",
        "electrons_per_oxygenation",
        "exp_id",
        "gbw",
        "gm_at_25",
        "gm_Ha",
        "gm_Hd",
        "gm_S",
        "Gs_min",
        "Gstar_at_25",
        "Gstar_Ea",
        "Jmax_at_25",
        "Jmax_Ea",
        "Kc_at_25",
        "Kc_Ea",
        "Ko_at_25",
        "Ko_Ea",
        "O2",
        "phi_PSII_0",
        "phi_PSII_1",
        "phi_PSII_2",
        "rh",
        "RL_at_25",
        "RL_Ea",
        "StomataWS",
        "temp",
        "theta_0",
        "theta_1",
        "theta_2",
        "Tleaf",
        "Tp_at_25",
        "Tp_Ha",
        "Tp_Hd",
        "Tp_S",
        "Vcmax_at_25",
        "Vcmax_Ea"};
}

string_vector c3photo_module::get_outputs()
{
    return {
        "Assim",
        "Assim_conductance",
        "Cc",
        "Ci",
        "Cs",
        "GrossAssim",
        "Gs",
        "iteration_C3_Assim",
        "penalty",
        "residual_C3_Assim",
        "RHs",
        "RL",
        "Rp"};
}

void c3photo_module::do_operation() const
{
    if (c3_model_type != 1.0 && c3_model_type != 2.0) {
        throw std::invalid_argument(
            "c3_model_type must be 1 (FvCB) or 2 (ePhotosynthesis)");
    }

    c3_temperature_response_parameters const temperature_parameters{
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

    photosynthesis_outputs const result = c3photoC(
        temperature_parameters,
        absorbed_ppfd,
        Tleaf,
        temp,
        rh,
        gm_at_25,
        Gstar_at_25,
        Kc_at_25,
        Ko_at_25,
        Vcmax_at_25,
        Jmax_at_25,
        Tp_at_25,
        RL_at_25,
        b0,
        b1,
        Gs_min,
        Catm,
        atmospheric_pressure,
        O2,
        StomataWS,
        electrons_per_carboxylation,
        electrons_per_oxygenation,
        beta_PSII,
        gbw,
        exp_id,
        static_cast<int>(c3_model_type));

    update(Assim_op, result.Assim);
    update(Assim_conductance_op, result.Assim_conductance);
    update(Cc_op, result.Cc);
    update(Ci_op, result.Ci);
    update(Cs_op, result.Cs);
    update(GrossAssim_op, result.GrossAssim);
    update(Gs_op, result.Gs);
    update(iteration_C3_Assim_op, result.iteration);
    update(penalty_op, result.penalty);
    update(residual_C3_Assim_op, result.residual);
    update(RHs_op, result.RHs);
    update(RL_op, result.RL);
    update(Rp_op, result.Rp);
}
