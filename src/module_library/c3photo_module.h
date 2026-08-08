#ifndef C3PHOTO_MODULE_H
#define C3PHOTO_MODULE_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace BMLePhoto
{
/**
 * @class c3photo_module
 *
 * @brief Exposes `c3photoC` as a BioCro direct module.
 *
 * This module calculates leaf-level C3 photosynthesis at a prescribed leaf
 * temperature and absorbed PPFD. The biochemical model is selected by
 * `c3_model_type`: 1 for FvCB and 2 for ePhotosynthesis.
 */
class c3photo_module : public direct_module
{
   public:
    c3photo_module(
        state_map const& input_quantities,
        state_map* output_quantities);

    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "c3photo"; }

   private:
    double const& absorbed_ppfd;
    double const& atmospheric_pressure;
    double const& b0;
    double const& b1;
    double const& beta_PSII;
    double const& Catm;
    double const& c3_model_type;
    double const& electrons_per_carboxylation;
    double const& electrons_per_oxygenation;
    double const& exp_id;
    double const& gbw;
    double const& gm_at_25;
    double const& gm_Ha;
    double const& gm_Hd;
    double const& gm_S;
    double const& Gs_min;
    double const& Gstar_at_25;
    double const& Gstar_Ea;
    double const& Jmax_at_25;
    double const& Jmax_Ea;
    double const& Kc_at_25;
    double const& Kc_Ea;
    double const& Ko_at_25;
    double const& Ko_Ea;
    double const& O2;
    double const& phi_PSII_0;
    double const& phi_PSII_1;
    double const& phi_PSII_2;
    double const& rh;
    double const& RL_at_25;
    double const& RL_Ea;
    double const& StomataWS;
    double const& temp;
    double const& theta_0;
    double const& theta_1;
    double const& theta_2;
    double const& Tleaf;
    double const& Tp_at_25;
    double const& Tp_Ha;
    double const& Tp_Hd;
    double const& Tp_S;
    double const& Vcmax_at_25;
    double const& Vcmax_Ea;

    double* Assim_op;
    double* Assim_conductance_op;
    double* Cc_op;
    double* Ci_op;
    double* Cs_op;
    double* GrossAssim_op;
    double* Gs_op;
    double* iteration_C3_Assim_op;
    double* penalty_op;
    double* residual_C3_Assim_op;
    double* RHs_op;
    double* RL_op;
    double* Rp_op;

    void do_operation() const;
};

}  // namespace BMLePhoto

#endif
