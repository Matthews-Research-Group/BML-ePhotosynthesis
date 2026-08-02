#include <algorithm>
#include <cmath>
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

    photosynthesis_outputs photo;
    energy_balance_outputs et;

    auto calculate_gs = [=, &photo, &et](double current_gs) {
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
                electrons_per_oxygenation, beta_PSII, et.gbw_molar, exp_id, 2);

        return photo.Gs;
    };

    // Use a short fixed-point solve as a fast path for easy, contractive
    // cases. If it does not converge promptly, switch to bracketed methods.
    int constexpr fast_path_iterations = 10;
    root_finding::fixed_point fixed_point_solver(
        fast_path_iterations,
        1e-3,
        1e-3);
    root_finding::result_t result =
        fixed_point_solver.solve(calculate_gs, initial_stomatal_conductance);

    if (!root_finding::is_successful(result.flag)) {
        // Solve the self-consistency equation calculated_gs(gs) - gs = 0.
        // Unlike direct fixed-point iteration, the bracketed methods do not
        // require calculated_gs to be a contraction near the solution.
        auto gs_residual = [&](double current_gs) {
            double const residual = calculate_gs(current_gs) - current_gs;
            if (!std::isfinite(residual)) {
                throw std::runtime_error(
                    "c3_leaf_photosynthesis conductance residual is not "
                    "finite at gs = " +
                    std::to_string(current_gs) +
                    ", absorbed PPFD = " + std::to_string(absorbed_ppfd) +
                    ", wind speed = " + std::to_string(windspeed) +
                    ", canopy height = " + std::to_string(height));
            }
            return residual;
        };

        // Ball-Berry conductance cannot be less than its
        // water-stress-adjusted intercept, so this provides a physical lower
        // bound.
        double const b0_adjusted =
            StomataWS * b0 + Gs_min * (1.0 - StomataWS);
        double const gs_lower = std::max(1e-8, b0_adjusted);
        double gs_upper = std::max(
            1.0,
            2.0 * std::max(initial_stomatal_conductance, gs_lower));

        double const residual_lower = gs_residual(gs_lower);
        double residual_upper = gs_residual(gs_upper);

        int bracket_expansions = 0;
        int constexpr max_bracket_expansions = 10;
        while (root_finding::same_signs(residual_lower, residual_upper) &&
               bracket_expansions < max_bracket_expansions) {
            gs_upper *= 2.0;
            residual_upper = gs_residual(gs_upper);
            ++bracket_expansions;
        }

        if (root_finding::same_signs(residual_lower, residual_upper)) {
            throw std::runtime_error(
                "c3_leaf_photosynthesis could not bracket a self-consistent "
                "stomatal conductance after the fixed-point fast path failed: "
                "lower gs = " +
                std::to_string(gs_lower) +
                ", lower residual = " + std::to_string(residual_lower) +
                ", upper gs = " + std::to_string(gs_upper) +
                ", upper residual = " + std::to_string(residual_upper) +
                ", absorbed PPFD = " + std::to_string(absorbed_ppfd) +
                ", wind speed = " + std::to_string(windspeed) +
                ", canopy height = " + std::to_string(height));
        }

        root_finding::dekker dekker_solver(100, 1e-3, 1e-3);
        result = dekker_solver.solve(gs_residual, gs_lower, gs_upper);

        if (!root_finding::is_successful(result.flag)) {
            root_finding::illinois illinois_solver(200, 1e-3, 1e-3);
            result = illinois_solver.solve(gs_residual, gs_lower, gs_upper);
        }

        if (!root_finding::is_successful(result.flag)) {
            throw std::runtime_error(
                "c3_leaf_photosynthesis conductance solvers failed after the "
                "fixed-point fast path. Termination flag:\n    " +
                root_finding::flag_message(result.flag) +
                "\n    lower gs = " + std::to_string(gs_lower) +
                ", upper gs = " + std::to_string(gs_upper) +
                ", absorbed PPFD = " + std::to_string(absorbed_ppfd) +
                ", wind speed = " + std::to_string(windspeed) +
                ", canopy height = " + std::to_string(height));
        }
    }

    // Ensure all side-effect outputs correspond to the accepted conductance.
    calculate_gs(result.root);

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
