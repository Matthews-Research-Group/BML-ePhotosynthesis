#ifndef C3PHOTO_H
#define C3PHOTO_H

#include "photosynthesis_outputs.h"   // for photosynthesis_outputs
#include "c3_temperature_response.h"  // for c3_temperature_response_parameters

photosynthesis_outputs c3photoC(
    c3_temperature_response_parameters const tr_param,
    double const absorbed_ppfd,
    double const Tleaf,
    double const Tambient,
    double const RH,
    double const Vcmax_at_25,
    double const Jmax_at_25,
    double const TPU_rate_max,
    double const RL_at_25,
    double const bb0,
    double const bb1,
    double const Gs_min,
    double Ca,
    double const AP,
    double const O2,
    double const StomWS,
    double const electrons_per_carboxylation,
    double const electrons_per_oxygenation,
    double const beta_PSII,
    double const gbw,
    double const exp_id,
    int    const model_type);

double solc(double LeafT);
double solo(double LeafT);

double Vcmax_multiplier(double T_kelvin);

#endif
