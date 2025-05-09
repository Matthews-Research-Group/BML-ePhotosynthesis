/* This file will contain functions which are common to several */
/* routines in the BioCro package. These are functions needed */
/* internally. The normal user will not need them */

#include <stdexcept>
#include <string>
#include <cmath>
#include "BioCro.h"
#include "../framework/constants.h"  // for pi, e, atmospheric_pressure_at_sea_level,
                           // ideal_gas_constant, molar_mass_of_water,
                           // stefan_boltzmann, celsius_to_kelvin

double poisson_density(int x, double lambda)
{
    // The probability density for the Poisson distribution is
    // e^-lambda * lambda^x / x!
    // The factorial term produces numbers too large to hold, so perform the calculation in the log domain.
    // n! can be estimated by the approximation sqrt(2 * pi * x) * (x / e)^x.
    using math_constants::pi;
    using math_constants::e;

    // Stirling's approximation for n!:
    double factorial_x = sqrt(2 * pi * x) * pow((x / e), x);

    double log_result = -lambda + x * log(lambda) - log(factorial_x);

    return exp(log_result);
}

/**
 *  @brief Calculates some properties of the "light macro environment," i.e.,
 *  the light just above the top of the canopy.
 *
 *  @param [in] cosine_zenith_angle The cosine of the solar zenith angle; when
 *      the Sun is directly overhead, the angle is 0 and its cosine is 1; when
 *      the Sun is at the horizon, the angle is 90 and its cosine is 0
 *
 *  @param [in] atmospheric_pressure Local atmospheric pressure (in Pa)
 *
 *  The return value from this function is a customized structure
 *  (`light_model`) with the following members:
 *  - `light_model.direct_irradiance_fraction`: The fraction of irradiance that
 *    is direct radiation (dimensionless)
 *  - `light_model.diffuse_irradiance_fraction`: The fraction of irradiance
 *    that is diffuse radiation (dimensionless)
 *
 *  The basis for this function is given in chapter 11 of Norman and Campbell,
 *  _An Introduction to Environmental Biophysics_, 2nd edition.
 */
Light_model lightME(double cosine_zenith_angle, double atmospheric_pressure)
{
    double direct_irradiance_transmittance;
    double diffuse_irradiance_transmittance;

    if (cosine_zenith_angle <= 0) { // Check if the Sun is at or below the horizon. If it is, directly set direct_irradiance_transmittance and diffuse_irradiance_transmittance because the equations below would otherwise produce erroneous results. The Sun is at or below the horizon when cosine_zenith_angle is less than or equal to 0.
        direct_irradiance_transmittance = 0;
        diffuse_irradiance_transmittance = 1;
    } else { // If the Sun is above the horizon, calculate diffuse and direct irradiance from the Sun's angle to the ground and the path through the atmosphere.
        constexpr double atmospheric_transmittance = 0.85; // dimensionless.
        const double pressure_ratio = atmospheric_pressure / physical_constants::atmospheric_pressure_at_sea_level; // dimensionless.
        constexpr double proportion_of_irradiance_scattered = 0.3; // dimensionless.

        direct_irradiance_transmittance = pow(atmospheric_transmittance, (pressure_ratio / cosine_zenith_angle)); // Modified from equation 11.11 of Norman and Campbell.
        diffuse_irradiance_transmittance = proportion_of_irradiance_scattered * (1 - direct_irradiance_transmittance) * cosine_zenith_angle; // Modified from equation 11.13 of Norman and Campbell.
    }

    Light_model light_model;
    light_model.direct_irradiance_fraction = direct_irradiance_transmittance / (direct_irradiance_transmittance + diffuse_irradiance_transmittance); // dimensionless.
    light_model.diffuse_irradiance_fraction = diffuse_irradiance_transmittance / (direct_irradiance_transmittance + diffuse_irradiance_transmittance); // dimensionless.

    return light_model;
}

/**
 *  @brief Computes total absorbed shortwave radiation from the
 *  photosynthetically active photon flux density (PPFD) incident on a leaf
 *
 *  @param [in] incident_ppfd Photosynthetically active photon flux density
 *              (PPFD) incident on a leaf expressed in micromol / m^2 / s
 *
 *  @param [in] par_energy_content The energy content of PPFD expressed in J /
 *              micromol
 *
 *  @param [in] par_energy_fraction The fraction of total shortwave energy
 *              contained in the PAR band expressed as a real number between 0
 *              and 1
 *
 *  @param [in] leaf_reflectance The fractional amount of shortwave radiation
 *              reflected by the leaf (weighted across all shortwave radiation)
 *
 *  @param [in] leaf_transmittance The fractional amount of shortwave radiation
 *              transmitted through the leaf (weighted across all shortwave
 *              radiation)
 *
 *  @return The total shortwave radiation absorbed by the leaf expressed in
 *          J / m^2 / s
 *
 *  The total absorbed shortwave radiation is determined using the following
 *  steps:
 *  - determine the incident photosynthetically active radiation (PAR) using the
 *    energy content of PAR
 *  - determine the incident near-infrared radiation (NIR) from the incident PAR
 *    using the PAR energy fraction
 *  - determine the total incident radiation by adding the incident PAR and NIR
 *  - determine the total absorbed radiation using the leaf's reflection and
 *    transmission coefficients
 */
double absorbed_shortwave_from_incident_ppfd(
    double incident_ppfd,        // micromol / m^2 / s
    double par_energy_content,   // J / micromol
    double par_energy_fraction,  // dimensionless
    double leaf_reflectance,     // dimensionless
    double leaf_transmittance    // dimensionless
)
{
    double incident_par = incident_ppfd * par_energy_content;  // J / m^2 / s

    double incident_nir = incident_par *
                          (1 - par_energy_fraction) /
                          par_energy_fraction;  // J / m^2 /s

    double incident_shortwave = incident_par + incident_nir;  // J / m^2 / s

    return incident_shortwave *
           (1 - leaf_reflectance - leaf_transmittance) /
           (1 - leaf_transmittance);  // J / m^2 / s
}

/**
 *  @brief Computes an n-layered light profile from the direct light, diffuse
 *  light, leaf area index, solar zenith angle, and other parameters.
 *
 *  @param [in] ambient_ppfd_beam Photosynthetically active photon flux density
 *              (PPFD) for beam light passing through a surface perpendicular
 *              to the beam direction at the top of the canopy; this represents
 *              direct sunlight for a plant in a field
 *              (micromol / (m^2 beam) / s)
 *
 *  @param [in] ambient_ppfd_diffuse Photosynthetically active photon flux
 *              density (PPFD) for diffuse light at the top of the canopy; this
 *              represents diffuse light scattered out of the solar beam by the
 *              Earth's atmosphere for a plant in a field; as a diffuse flux
 *              density, this represents the flux through any surface
 *              (micromol / m^2 / s)
 *
 *  @param [in] lai Leaf area index (LAI) of the entire canopy, which represents
 *              the leaf area per unit of ground area (dimensionless from m^2
 *              leaf / m^2 ground)
 *
 *  @param [in] nlayers Integer number of layers in the canopy
 *
 *  @param [in] cosine_zenith_angle Cosine of the solar zenith angle
 *              (dimensionless)
 *
 *  @param [in] kd Extinction coefficient for diffuse light (dimensionless)
 *
 *  @param [in] chil Ratio of average projected areas of canopy elements on
 *              horizontal surfaces; for a spherical leaf distribution,
 *              `chil = 0`; for a vertical leaf distribution, `chil = 1`; for a
 *              horizontal leaf distribution, `chil` approaches infinity
 *              (dimensionless from m^2 / m^2)
 *
 *  @param [in] absorptivity The leaf absorptivity on a quantum basis
 *              (dimensionless from mol / mol)
 *
 *  @param [in] heightf Leaf area density, i.e., LAI per height of canopy (m^-1
 *              from m^2 leaf / m^2 ground / m height)
 *
 *  @return An n-layered light profile representing quantities within
 *          the canopy, including several photon flux densities and
 *          the relative fractions of shaded and sunlit leaves
 */
Light_profile sunML(
    double ambient_ppfd_beam,     // micromol / (m^2 beam) / s
    double ambient_ppfd_diffuse,  // micromol / m^2 / s
    double lai,                   // dimensionless from m^2 / m^2
    int nlayers,                  // dimensionless
    double cosine_zenith_angle,   // dimensionless
    double kd,                    // dimensionless
    double chil,                  // dimensionless from m^2 / m^2
    double absorptivity,          // dimensionless from mol / mol
    double heightf,               // m^-1 from m^2 leaf / m^2 ground / m height
    double par_energy_content,    // J / micromol
    double par_energy_fraction,   // dimensionless
    double leaf_transmittance,    // dimensionless
    double leaf_reflectance       // dimensionless
)
{
    if (nlayers < 1 || nlayers > MAXLAY) {
        throw std::out_of_range("nlayers must be at least 1 but no more than MAXLAY.");
    }
    if (cosine_zenith_angle > 1 || cosine_zenith_angle < -1) {
        throw std::out_of_range("cosine_zenith_angle must be between -1 and 1.");
    }
    if (kd > 1 || kd < 0) {
        throw std::out_of_range("kd must be between 0 and 1.");
    }
    if (chil < 0) {
        throw std::out_of_range("chil must be non-negative.");
    }
    if (absorptivity > 1 || absorptivity < 0) {
        throw std::out_of_range("absorptivity must be between 0 and 1.");
    }
    if (heightf <= 0) {
        throw std::out_of_range("heightf must greater than zero.");
    }

    // Calculate the leaf shape factor for an ellipsoidal leaf angle
    // distribution using the equation from page 251 of Campbell & Norman
    // (1998). We will use this value as `k`, the canopy extinction coefficient
    // for photosynthetically active radiation throughout the canopy. This
    // quantity represents the ratio of horizontal area to total area for leaves
    // in the canopy and is therefore dimensionless from
    // (m^2 ground) / (m^2 leaf).
    double zenith_angle = acos(cosine_zenith_angle);  // radians
    double k0 = sqrt(pow(chil, 2) + pow(tan(zenith_angle), 2));
    double k1 = chil + 1.744 * pow((chil + 1.182), -0.733);
    double k = k0 / k1;  // dimensionless

    double lai_per_layer = lai / nlayers;

    // Calculate the fraction of direct radiation that passes through the canopy
    // using Equation 15.1. Note that this is equivalent to the fraction of
    // ground area below the canopy that is exposed to direct sunlight. Note
    // that if the sun is at or below the horizon, no part of the soil is
    // sunlit; this corresponds to the case where cosine_zenith_angle is close
    // to or below zero.
    double canopy_direct_transmission_fraction =
        cosine_zenith_angle <= 1E-10 ? 0.0 : exp(-k * lai); // dimensionless

    // Calculate the ambient direct PPFD through a surface parallel to the ground
    const double ambient_ppfd_beam_ground = ambient_ppfd_beam * cosine_zenith_angle;  // micromol / (m^2 ground) / s

    // Calculate the ambient direct PPFD through a unit area of leaf surface
    double ambient_ppfd_beam_leaf = ambient_ppfd_beam_ground * k;  // micromol / (m^2 leaf) / s

    // Start to fill in the light profile values
    Light_profile light_profile;
    light_profile.canopy_direct_transmission_fraction = canopy_direct_transmission_fraction;

    // Fill in the layer-dependent light profile values
    for (int i = 0; i < nlayers; ++i) {
        // Get the cumulative LAI for this layer, which represents the total
        // leaf area above this layer
        const double cumulative_lai = lai_per_layer * (i + 0.5);

        // Calculate the amount of PPFD scattered out of the direct beam using
        // Equations 15.6 and 15.1 from Campbell & Norman (1998), following
        // example 15.2. This is a diffuse flux density representing the flux
        // through any surface.
        const double scattered_ppfd =
            ambient_ppfd_beam_ground * (exp(-k * sqrt(absorptivity) * cumulative_lai) -
                                        exp(-k * cumulative_lai));  // micromol / m^2 / s

        // Calculate the total flux of diffuse photosynthetically active light
        // in this layer by combining the scattered PPFD with the ambient
        // diffuse PPFD. Here we use Equation 15.6 with `alpha` = 1 and
        // `kbe(phi)` = kd.
        double diffuse_ppfd =
            ambient_ppfd_diffuse * exp(-kd * cumulative_lai) + scattered_ppfd;  // micromol / m^2 / s

        // Calculate the fraction of sunlit and shaded leaves in this canopy
        // layer using Equation 15.21.
        const double Ls = (1 - exp(-k * lai_per_layer)) * exp(-k * cumulative_lai) / k;  // dimensionless
        double sunlit_fraction = Ls / lai_per_layer;                                     // dimensionless
        double shaded_fraction = 1 - sunlit_fraction;                                    // dimensionless

        // Calculate an "average" incident PPFD for the sunlit and shaded leaves
        // that doesn't seem to be based on a formula from Campbell & Norman
        // (1998). It's interpreted as a flux density through a unit of leaf
        // area, but that may not be correct.
        double average_ppfd =
            (sunlit_fraction * (ambient_ppfd_beam_leaf + diffuse_ppfd) + shaded_fraction * diffuse_ppfd) *
            (1 - exp(-k * lai_per_layer)) / k;  // micromol / (m^2 leaf) / s

        // For values of cosine_zenith_angle close to or less than 0, in place
        // of the calculations above, we want to use the limits of the above
        // expressions as cosine_zenith_angle approaches 0 from the right:
        if (cosine_zenith_angle <= 1E-10) {
            ambient_ppfd_beam_leaf = ambient_ppfd_beam / k1;
            diffuse_ppfd = ambient_ppfd_diffuse * exp(-kd * cumulative_lai);
            sunlit_fraction = 0;
            shaded_fraction = 1;
            average_ppfd = 0;
        }

        // Store these values of incident PPFD
        light_profile.sunlit_incident_ppfd[i] = ambient_ppfd_beam_leaf + diffuse_ppfd;  // micromole / (m^2 leaf) / s
        light_profile.incident_ppfd_scattered[i] = scattered_ppfd;                      // micromole / m^2 / s
        light_profile.shaded_incident_ppfd[i] = diffuse_ppfd;                           // micromole / (m^2 leaf) / s
        light_profile.average_incident_ppfd[i] = average_ppfd;                          // micromole / (m^2 leaf) / s
        light_profile.sunlit_fraction[i] = sunlit_fraction;                             // dimensionless from m^2 / m^2
        light_profile.shaded_fraction[i] = shaded_fraction;                             // dimensionless from m^2 / m^2
        light_profile.height[i] = (lai - cumulative_lai) / heightf;                     // m

        // We also need to determine the total amount of absorbed solar energy
        // for sunlit and shaded leaves
        light_profile.sunlit_absorbed_shortwave[i] =
            absorbed_shortwave_from_incident_ppfd(
                ambient_ppfd_beam_leaf + diffuse_ppfd,
                par_energy_content,
                par_energy_fraction,
                leaf_reflectance,
                leaf_transmittance);  // J / (m^2 leaf) / s

        light_profile.shaded_absorbed_shortwave[i] =
            absorbed_shortwave_from_incident_ppfd(
                diffuse_ppfd,
                par_energy_content,
                par_energy_fraction,
                leaf_reflectance,
                leaf_transmittance);  // J / (m^2 leaf) / s

        light_profile.average_absorbed_shortwave[i] =
            absorbed_shortwave_from_incident_ppfd(
                average_ppfd,
                par_energy_content,
                par_energy_fraction,
                leaf_reflectance,
                leaf_transmittance);  // J / (m^2 leaf) / s
    }
    return light_profile;
}


/* Additional Functions needed for EvapoTrans */


/**
 * @brief Wind profile function
 *
 * Preconditions:
 *     `WindSpeed` is non-negative.
 *     `LAI` is non-negative
 *     `nlayers` is at least 1 and at most MAXLAY.
 *     `wind_speed_profile` is an array of at least size `nlayers`.
 */
void WINDprof(double WindSpeed, double LAI, int nlayers,
              double* wind_speed_profile)
{
    constexpr double k = 0.7;
    double LI = LAI / nlayers;

    for (int i = 0; i < nlayers; ++i)
    {
        double CumLAI = LI * (i + 1);
        wind_speed_profile[i] = WindSpeed * exp(-k * (CumLAI - LI));
    }
}

/**
 * @brief Calculates a relative humidity profile throughout a multilayer
 * canopy.
 *
 * @param[in] RH relative humidity just above the canopy `(0 <= RH <= 1)`
 *
 * @param[in] nlayers number of layers in the canopy `(1 <= nlayers <= MAXLAY)`
 *
 * @param[out] relative_humidity_profile array of relative humidity values
 * expressed as fractions between 0 and 1, where the value at index `i`
 * represents relative humidity at the bottom of canopy layer `i` and `i = 0`
 * corresponds to the top canopy layer.
 *
 * One can derive an expression for the relative humidity (`h`) throughout a
 * plant canopy by making the following assumptions:
 * - `h` at the top of the canopy is the same as the ambient value (`h = h0`)
 * - `h` at the bottom of the canopy is one hundred percent (`h = 1`)
 * - `h` follows an exponential profile throughout the canopy
 *
 * To enforce the exponential profile, we can write
 *
 * `h(x) = A * exp(B * x)`              [Equation (1)]
 *
 * where `x` is a normalized expression of depth within the canopy (`x = 0` at
 * the canopy top and `x = 1` at the bottom) and the values of `A` and `B` are
 * yet to be determined.
 *
 * To enforce the `h` value at the canopy top, we require `h0 = h(0) = A`; in
 * other words, `A = h0`. To enforce the `h` value at the canopy bottom, we
 * require `1 = h(1) = h0 * exp(B)`; in other words, `B = -ln(h0)`. Putting it
 * all together, we see that under these assumptions, `h` throughout the canopy
 * is given by
 *
 * `h(x) = h0 * exp(-ln(h0) * x)`       [Equation (2)]
 *
 * If we additionally assume that `h0` is close to 1, we can simplify
 * `B = -ln(h0)` by taking just the linear part of the Taylor series for
 * `-ln(h0)` centered at `h0 = 1`, which is `B = -(h0 - 1)`. With this
 * modification, we have
 *
 * `h(x) = h0 * exp((1 - h0) * x)`      [Equation (3)]
 *
 * When `h0 = 1`, Equations (2) and (3) are in perfect agreement, reducing to
 * `h(x) = 1`. However, as `h0` deviates further from `h0 = 1`, the
 * linearization for `B` becomes less accurate and the two versions diverge,
 * especially deeper in the canopy. For example, when `h0 = 0.7`, Equation (2)
 * becomes `h(x) = 0.7 * exp(0.357 * x)` while Equation (3) becomes
 * `h(x) = 0.7 * exp(0.300 * x)`. Both versions agree at the top of the canopy
 * where `x = 0` and `h(x) = 1`, but are different at the bottom: 1.00 vs. 0.94.
 * For lower `h0` values, the difference at the canopy bottom becomes even more
 * pronounced, significantly violating one of the original assumptions. However,
 * there isn't a strong scientific justification for assuming humidity is 1 at
 * the bottom of *every* canopy, so the error due to the approximation for `B`
 * is deemed to be acceptable.
 *
 * ---
 *
 * In BioCro, we divide the canopy into equally sized layers; i.e., the interval
 * `0 <= x <= 1` is divided into `n` segments of length `1 / n` by `n + 1`
 * boundaries occurring at `x = 0`, `x = 1 / n`, `x = 2 / n`, ..., `x = 1`. Here
 * we wish to find the `h` value at the bottom of each layer, so we use Equation
 * (3) with `x = 1 / n`, `x = 2 / n`, ..., `x = 1`.
 *
 * In the code below, the `RH` input argument corresponds to `h0`, the
 * exponential growth constant `kh` corresponds to `B`, and `nlayers`
 * corresponds to `n`.
 *
 * ---
 *
 * Note 1: the explanation for this code was "reverse engineered" by EBL from
 * some crytpic documentation found in earlier versions of this function and
 * from Stephen Humphries's thesis, which is not available online:
 *
 * Humphries, S. "Will mechanistically rich models provide us with new insights
 * into the response of plant production to climate change?: development and
 * experiments with WIMOVAC: (Windows Intuitive Model of Vegetation response
 * to Atmosphere & Climate Change)" (University of Essex, 2002).
 *
 * ---
 *
 * Note 2: Equation (2) can be rewritten by noting that `h0 = exp(ln(h0))`. With
 * this replacement, the equation becomes
 *
 * `h(x) = exp(ln(h0) - ln(h0) * x) = exp(-B * (1 - x))`
 *
 * The equation takes this form in Humphries's thesis.
 *
 */
void RHprof(double RH, int nlayers, double* relative_humidity_profile)
{
    if (RH > 1 || RH < 0) {
        throw std::out_of_range("RH must be between 0 and 1.");
    }
    if (nlayers < 1 || nlayers > MAXLAY) {
        throw std::out_of_range("nlayers must be at least 1 but no more than MAXLAY.");
    }

    const double kh = 1 - RH;

    for (int i = 0; i < nlayers; ++i)
    {
        double j = i + 1;  // Explicitly make j a double so that j / nlayers isn't truncated.
        double temp_rh = RH * exp(kh * (j / nlayers));
        relative_humidity_profile[i] = temp_rh;
    }
}

void LNprof(double LeafN, double LAI, int nlayers, double kpLN, double* leafN_profile)
{
    double LI = LAI / nlayers;
    for(int i = 0; i < nlayers; ++i)
    {
        double CumLAI = LI * (i + 1);
        leafN_profile[i] = LeafN * exp(-kpLN * (CumLAI - LI));
    }
}

/**
 *  @brief Determines the density of dry air from the air temperature.
 *
 *  This function appears to be derived from fitting a linear equation to the
 *  values of `rho` in Table 14.3 on pg. 408 in Thornley and Johnson (1990):
 *
 *  | T (degrees C) | `rho` (kg / m^3) |
 *  | :-----------: | :--------------: |
 *  |  0            | 1.29             |
 *  |  5            | 1.27             |
 *  | 10            | 1.25             |
 *  | 15            | 1.23             |
 *  | 20            | 1.20             |
 *  | 25            | 1.18             |
 *  | 30            | 1.16             |
 *  | 35            | 1.15             |
 *  | 40            | 1.13             |
 *
 *  Thornley, J.H.M. and Johnson, I.R. (1990) Plant and Crop Modelling. A
 *  Mathematical Approach to Plant and Crop Physiology.
 *
 *  @param [in] air_temperature Air temperature in degrees C
 *
 *  @return Density of dry air in kg / m^3
 */
double TempToDdryA(
    double air_temperature  // degrees C
)
{
    return 1.295163636 + -0.004258182 * air_temperature;  // kg / m^3
}

/**
 *  @brief Calculate the density of dry air from temperature and pressure using
 *  the ideal gas law.
 *
 *  See https://en.wikipedia.org/wiki/Density_of_air for details.
 *
 *  @param [in] air_temperature Air temperature in kelvin
 *
 *  @param [in] air_pressure Air pressure in Pa
 *
 *  @return Density of dry air in kg / m^3
 */
double dry_air_density(
    const double air_temperature,  // kelvin
    const double air_pressure      // Pa
)
{
    constexpr double r_specific = physical_constants::ideal_gas_constant /
        physical_constants::molar_mass_of_dry_air; // J / kg / K

    return air_pressure / r_specific / air_temperature;  // kg / m^3
}

/**
 *  @brief Determine the latent heat of vaporization for water from its
 *  temperature.
 *
 *  This function appears to be derived from fitting a linear equation to the
 *  values of `lambda` in Table 14.3 on pg. 408 in Thornley and Johnson (1990):
 *
 *  | T (degrees C) | `lambda` (J / kg) |
 *  | :-----------: | :---------------: |
 *  |  0            | 2.50              |
 *  |  5            | 2.49              |
 *  | 10            | 2.48              |
 *  | 15            | 2.47              |
 *  | 20            | 2.45              |
 *  | 25            | 2.44              |
 *  | 30            | 2.43              |
 *  | 35            | 2.42              |
 *  | 40            | 2.41              |
 *
 *  Thornley, J.H.M. and Johnson, I.R. (1990) Plant and Crop Modelling. A
 *  Mathematical Approach to Plant and Crop Physiology.
 *
 *  @param [in] temperature Water temperature in degrees C
 *
 *  @return Latent heat of vaporization for water in J / kg
 */
double TempToLHV(
    double temperature  // degrees C
)
{
    return 2501000 + -2372.727 * temperature;  // J / kg.
}

/**
 *  @brief Determine the latent heat of vaporization of water from its
 *  temperature.
 *
 *  Here we use Equation 8 from
 *  [Henderson-Sellers, B. "A new formula for latent heat of vaporization of water as a function of
 *  temperature" Quarterly Journal of the Royal Meteorological Society 110, 1186–1190 (1984)]
 *  (https://doi.org/10.1002/qj.49711046626)
 *
 *  @param [in] temperature Water temperature in kelvin
 *
 *  @return Latent heat of vaporization for water in J / kg
 */
double water_latent_heat_of_vaporization_henderson(
    double temperature)  // kelvin
{
    return 1.91846e6 * pow(temperature / (temperature - 33.91), 2);  // J / kg
}

/**
 *  @brief Determine the derivative of saturation water vapor pressure with
 *  respect to temperature at a particular value of air temperature.
 *
 *  This function appears to be derived from fitting a quadratic function to the
 *  values of `s` in Table 14.3 on pg. 408 in Thornley and Johnson (1990):
 *
 *  | T (degrees C) | `s` (10^(-3) kg / m^3 / K) |
 *  | :-----------: | :------------------------: |
 *  |  0            | 0.33                       |
 *  |  5            | 0.45                       |
 *  | 10            | 0.60                       |
 *  | 15            | 0.78                       |
 *  | 20            | 1.01                       |
 *  | 25            | 1.30                       |
 *  | 30            | 1.65                       |
 *  | 35            | 2.07                       |
 *  | 40            | 2.57                       |
 *
 *  Thornley, J.H.M. and Johnson, I.R. (1990) Plant and Crop Modelling. A
 *  Mathematical Approach to Plant and Crop Physiology.
 *
 *  @param [in] air_temperature Air temperature in degrees C
 *
 *  @return Derivative of saturation water vapor pressure with respect to
 *  temperature in kg / m^3 / K (equivalent to Pa / K)
 */
double TempToSFS(
    double air_temperature  // degrees C
)
{
    return (0.338376068 + 0.011435897 * air_temperature + 0.001111111 *
            pow(air_temperature, 2)) * 1e-3;  //  kg / m^3 / K
}

/**
 * @brief Determine saturation water vapor pressure (Pa) from air temperature
 * (degrees C) using the Arden Buck equation.
 *
 * Equations of this form were used by Arden Buck to model saturation water
 * vapor pressure in 1981: [A. L. Buck, J. Appl. Meteor. 20, 1527–1532 (1981)]
 * (https://doi.org/10.1175/1520-0450(1981)020%3C1527:NEFCVP%3E2.0.CO;2)
 *
 * In 1996, some of the fitting coefficients were updated. Here we use these values,
 * found in [a Buck hygrometer manual]
 * (http://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf),
 * and also displayed on the [Wikipedia page for the Arden Buck equation]
 * (https://en.wikipedia.org/wiki/Arden_Buck_equation).
 *
 * We use the values for vapor over water (rather than ice) and disregard the
 * "enhancement factor" since we are only concerned with the pressure at
 * saturation.
 *
 *  @param [in] air_temperature Air temperature in degrees C
 *
 *  @return Saturation water vapor pressure in Pa
 */
double saturation_vapor_pressure(
    double air_temperature  // degrees C
)
{
    double a = (18.678 - air_temperature / 234.5) * air_temperature;
    double b = 257.14 + air_temperature;
    return 611.21 * exp(a / b);  // Pa
}

struct ET_Str EvapoTrans2(
    double absorbed_shortwave_radiation_et,  // J / m^2 / s (used to calculate evapotranspiration rate)
    double absorbed_shortwave_radiation_lt,  // J / m^2 / s (used to calculate leaf temperature)
    double airTemp,                          // degrees C
    double RH,                               // dimensionless from Pa / Pa
    double WindSpeed,                        // m / s
    double stomatal_conductance,             // mmol / m^2 / s
    double leaf_width,                       // meter
    double specific_heat_of_air,             // J / kg / K
    double minimum_gbw,                      // mol / m^2 / s
    int eteq                                 // unitless parameter
)
{
    const double DdryA = TempToDdryA(airTemp);               // kg / m^3. Density of dry air.,
    const double LHV = TempToLHV(airTemp);                   // J / kg
    const double SlopeFS = TempToSFS(airTemp);               // kg / m^3 / K
    const double SWVP = saturation_vapor_pressure(airTemp);  // Pa.

    // TODO: This is for about 20 degrees C at 100000 Pa. Change it to use the
    // model state. (1 * R * temperature) / pressure
    double constexpr volume_of_one_mole_of_air = 24.39e-3;  // m^3 / mol

    double minimum_gbw_in_m_per_s = minimum_gbw * volume_of_one_mole_of_air;  // m / s

    if (stomatal_conductance <= 0) {
        throw std::range_error("Thrown in EvapoTrans2: stomatal conductance is not positive.");
    }

    double conductance_in_m_per_s = stomatal_conductance * 1e-3 * volume_of_one_mole_of_air;  // m / s

    if (RH > 1) {
        throw std::range_error("Thrown in EvapoTrans2: RH (relative humidity) is greater than 1.");
    }

    // Convert from vapor pressure to vapor density using the ideal gas law.
    // This is approximately right for temperatures what won't kill plants.
    const double SWVC =
        SWVP / physical_constants::ideal_gas_constant /
        (airTemp + conversion_constants::celsius_to_kelvin) * physical_constants::molar_mass_of_water;  // kg / m^3

    if (SWVC < 0) {
        throw std::range_error("Thrown in EvapoTrans2: SWVC is less than 0.");
    }

    const double PsycParam = DdryA * specific_heat_of_air / LHV;  // kg / m^3 / K

    const double vapor_density_deficit = SWVC * (1 - RH);  // kg / m^3

    const double ActualVaporPressure = RH * SWVP;  // Pa

    /* This is the original from WIMOVAC*/
    double Deltat = 0.01;  // degrees C
    double ga;
    double rlc; /* Long wave radiation for iterative calculation */
    {
        double ChangeInLeafTemp = 10.0;  // degrees C
        double Counter = 0;
        do {
            ga = leaf_boundary_layer_conductance_nikolov(
                WindSpeed, leaf_width, airTemp, Deltat, conductance_in_m_per_s,
                ActualVaporPressure, minimum_gbw_in_m_per_s);  // m / s

            /* In WIMOVAC, ga was added to the canopy conductance */
            /* ga = (ga * gbcW)/(ga + gbcW); */

            double OldDeltaT = Deltat;

            rlc = 4 * physical_constants::stefan_boltzmann * pow(conversion_constants::celsius_to_kelvin + airTemp, 3) * Deltat;  // W / m^2

            /* rlc = net long wave radiation emittted per second
             *     = radiation emitted per second - radiation absorbed per second
             *     = sigma * (Tair + deltaT)^4 - sigma * Tair^4
             *
             * To make it a linear function of deltaT, do a Taylor series about
             * deltaT = 0 and keep only the zero and first order terms.
             *
             * rlc = sigma * Tair^4 + deltaT * (4 * sigma * Tair^3) - sigma * Tair^4
             *     = 4 * sigma * Tair^3 * deltaT
             *
             * where 4 * sigma * Tair^3 is the derivative of
             * sigma * (Tair + deltaT)^4 evaluated at deltaT = 0
             */

            const double PhiN2 = absorbed_shortwave_radiation_lt - rlc;  // W / m^2

            /* This equation is from Thornley and Johnson pg. 418 */
            const double TopValue = PhiN2 * (1 / ga + 1 / conductance_in_m_per_s) - LHV * vapor_density_deficit;  // J / m^3
            const double BottomValue = LHV * (SlopeFS + PsycParam * (1 + ga / conductance_in_m_per_s));           // J / m^3 / K
            Deltat = fmin(fmax(TopValue / BottomValue, -10), 10);                                                 // kelvin. Confine Deltat to the interval [-10, 10]:

            ChangeInLeafTemp = fabs(OldDeltaT - Deltat);  // kelvin
        } while ((++Counter <= 10) && (ChangeInLeafTemp > 0.5));
    }

    /* Net radiation */
    const double PhiN = fmax(0, absorbed_shortwave_radiation_et - rlc);  // W / m^2

    const double penman_monteith =
        (SlopeFS * PhiN + LHV * PsycParam * ga * vapor_density_deficit) /
        (LHV * (SlopeFS + PsycParam * (1 + ga / conductance_in_m_per_s)));  // kg / m^2 / s.  Thornley and Johnson. 1990. Plant and Crop Modeling. Equation 14.4k. Page 408.

    const double EPen =
        (SlopeFS * PhiN + LHV * PsycParam * ga * vapor_density_deficit) /
        (LHV * (SlopeFS + PsycParam));  // kg / m^2 / s

    const double EPries = 1.26 * SlopeFS * PhiN / (LHV * (SlopeFS + PsycParam));  // kg / m^2 / s

    /* Choose equation to report */
    double TransR;
    switch (eteq) {
        case 1:
            TransR = EPen;
            break;
        case 2:
            TransR = EPries;
            break;
        default:
            TransR = penman_monteith;
            break;
    }

    // TransR has units of kg / m^2 / s. Convert to mmol / m^2 / s using the
    // molar mass of water (in kg / mol) and noting that 1e3 mmol = 1 mol
    double cf = 1e3 / physical_constants::molar_mass_of_water;  // mmol / kg for water

    struct ET_Str et_results;
    et_results.TransR = TransR * cf;                                         // mmol / m^2 / s
    et_results.EPenman = EPen * cf;                                          // mmol / m^2 / s
    et_results.EPriestly = EPries * cf;                                      // mmol / m^2 / s
    et_results.Deltat = Deltat;                                              // degrees C
    et_results.boundary_layer_conductance = ga / volume_of_one_mole_of_air;  // mol / m^2 / s

    return et_results;
}

/**
 *  @brief Caluclates the conductance for water vapor flow between the leaf
 *  surface and the atmosphere (AKA, the boundary layer conductance) using a
 *  model described in Nikolov, Massman, and Schoettle (1995).
 *
 *  In this model, two types of gas flow are considered: "forced" flow driven
 *  by wind-created eddy currents and "free" flow driven by temperature-related
 *  buoyancy effects. The overall conductance is determined to be the larger of
 *  the free and forced conductances.
 *
 *  In this function, we use equations 29, 33, 34, and 35 to calculate boundary
 *  layer conductance. This is the same approach taken in the `MLcan` model of
 *  Drewry et al. (2010).
 *
 *  In this model, the minimum possible boundary layer conductance that could
 *  occur is zero. This would happen if wind speed is zero and the air and leaf
 *  temperatures are the same. In realistic field conditions, boundary layer
 *  conductance can never truly be zero. To accomodate this, an option is
 *  provided for setting a minimum value for the boundary layer counductance.
 *
 *  References:
 *
 *  - [Nikolov, N. T., Massman, W. J. & Schoettle, A. W. "Coupling biochemical and biophysical processes at the
 *    leaf level: an equilibrium photosynthesis model for leaves of C3 plants" Ecological Modelling 80, 205–235 (1995)]
 *    (https://doi.org/10.1016/0304-3800(94)00072-P)
 *
 *  - [Drewry, D. T. et al. "Ecohydrological responses of dense canopies to environmental variability: 1. Interplay between
 *    vertical structure and photosynthetic pathway" Journal of Geophysical Research: Biogeosciences 115, (2010)]
 *    (https://doi.org/10.1029/2010JG001340)
 *
 *  @param [in] windspeed The wind speed in m / s
 *
 *  @param [in] leafwidth The characteristic leaf dimension in m
 *
 *  @param [in] air_temperature The air temperature in degrees C
 *
 *  @param [in] delta_t The temperature difference between the leaf and air in
 *              degrees C
 *
 *  @param [in] stomcond The stomatal conductance in m / s
 *
 *  @param [in] water_vapor_pressure The partial pressure of water vapor in the
 *              atmosphere in Pa
 *
 *  @param [in] minimum_gbw The lowest possible value for boundary layer
 *              conductance in m / s that should be returned
 *
 *  @return The boundary layer conductance in m / s
 */
double leaf_boundary_layer_conductance_nikolov(
    double windspeed,             // m / s
    double leafwidth,             // m
    double air_temperature,       // degrees C
    double delta_t,               // degrees C
    double stomcond,              // m / s
    double water_vapor_pressure,  // Pa
    double minimum_gbw            // m / s
)
{
    constexpr double p = physical_constants::atmospheric_pressure_at_sea_level;  // Pa

    double leaftemp = air_temperature + delta_t;                             // degrees C
    double gsv = stomcond;                                                   // m / s
    double Tak = air_temperature + conversion_constants::celsius_to_kelvin;  // K
    double Tlk = leaftemp + conversion_constants::celsius_to_kelvin;         // K
    double ea = water_vapor_pressure;                                        // Pa
    double lw = leafwidth;                                                   // m

    double esTl = saturation_vapor_pressure(leaftemp);  // Pa.

    // Forced convection
    constexpr double cf = 1.6361e-3;  // TODO: Nikolov et. al equation 29 use cf = 4.322e-3, not cf = 1.6e-3 as is used here.

    double gbv_forced = cf * pow(Tak, 0.56) * pow((Tak + 120) * ((windspeed / lw) / p), 0.5);  // m / s.

    // Free convection
    double gbv_free = gbv_forced;
    double eb = (gsv * esTl + gbv_free * ea) / (gsv + gbv_free);  // Pa. Eq 35

    double Tvdiff = (Tlk / (1 - 0.378 * eb / p)) - (Tak / (1 - 0.378 * ea / p));  // kelvin. It is also degrees C since it is a temperature difference. Eq. 34

    if (Tvdiff < 0) Tvdiff = -Tvdiff;

    gbv_free = cf * pow(Tlk, 0.56) * pow((Tlk + 120) / p, 0.5) * pow(Tvdiff / lw, 0.25);  // m / s. Eq. 33

    // Overall conductance
    double gbv = std::max(gbv_forced, gbv_free);  // m / s

    // Apply the minimum
    return std::max(gbv, minimum_gbw);  // m / s
}

/**
 *  @brief Caluclates the conductance for water vapor flow between the leaf
 *  surface and the atmosphere (AKA, the boundary layer conductance) using a
 *  model described in Thornley and Johnson (1990).
 *
 *  This model considers gas flow due to wind-driven eddy currents. Here, the
 *  conductance is calculated using Equation 14.9 from pages 414 - 416 of the
 *  Thornley textbook. Unfortunately, an electronic version of this reference is
 *  not available.
 *
 *  In this model, the minimum possible boundary layer conductance that could
 *  occur is zero, which would correspond to zero wind speed or canopy height.
 *  In realistic field conditions, boundary layer conductance can never truly be
 *  zero. To accomodate this, an option is provided for setting a minimum value
 *  for the boundary layer counductance.
 *
 *  This model contains two singularities, which occur when either of the
 *  following conditions are met:
 *
 *  - `WindSpeedHeight + Zeta - d = 0`, which is equivalent to `CanopyHeight =
 *     WindSpeedHeight / (dCoef - ZetaCoef) = WindSpeedHeight * 1.34`
 *
 *  - `WindSpeedHeight + ZetaM - d = 0`, which is equivalent to `CanopyHeight =
 *     WindSpeedHeight / (dCoef - ZetaMCoef) = WindSpeedHeight * 1.56`
 *
 *  So, as the canopy height approaches or exceeds the height at which wind
 *  speed was measured, the calculated boundary layer conductance becomes
 *  unbounded. For even larger canopy heights, the conductance eventually begins
 *  to decrease. For tall crops, this is a severe limitation to this model. Here
 *  we address this issue by limiting the canopy height to
 *  `0.98 * WindSpeedHeight`.
 *
 *  References:
 *
 *  - Thornley, J. H. M. & Johnson, I. R. "Plant and Crop Modelling: A
 *    Mathematical Approach to Plant and Crop Physiology" (1990)
 *
 *  @param [in] CanopyHeight The height of the canopy above the ground in m
 *
 *  @param [in] WindSpeed The wind speed in m / s as measured above the canopy
 *              at a reference height of five meters
 *
 *  @param [in] minimum_gbw The lowest possible value for boundary layer
 *              conductance in m / s that should be returned
 *
 *  @param [in] WindSpeedHeight The height in m at which the wind speed was
 *              measured
 *
 *  @return The boundary layer conductance in m / s
 */
double leaf_boundary_layer_conductance_thornley(
    double CanopyHeight,    // m
    double WindSpeed,       // m / s
    double minimum_gbw,     // m / s
    double WindSpeedHeight  // m
)
{
    // Define constants used in the model
    constexpr double kappa = 0.41;         // dimensionless. von Karmon's constant. Thornley and Johnson pgs 414 and 416.
    constexpr double ZetaCoef = 0.026;     // dimensionless, Thornley and Johnson 1990, Eq. 14.9o
    constexpr double ZetaMCoef = 0.13;     // dimensionless, Thornley and Johnson 1990, Eq. 14.9o
    constexpr double dCoef = 0.77;         // dimensionless, Thornley and Johnson 1990, Eq. 14.9o.
                                           // In the original text this value is reported as 0.64.
                                           // In the 2000 reprinting of this text, the authors state that this value should be 0.77.
                                           // See "Errata to the 2000 printing" on the page after the preface of the 2000 reprinting of the 1990 text.

    // Apply the height limit
    CanopyHeight = std::min(CanopyHeight, 0.98 * WindSpeedHeight);  // meters

    // Calculate terms that depend on the canopy height
    const double Zeta = ZetaCoef * CanopyHeight;    // meters
    const double Zetam = ZetaMCoef * CanopyHeight;  // meters
    const double d = dCoef * CanopyHeight;          // meters

    // Calculate the boundary layer conductance `ga` according to Thornley and
    // Johnson Eq. 14.9n, pg. 416
    const double ga0 = pow(kappa, 2) * WindSpeed;                   // m / s
    const double ga1 = log((WindSpeedHeight + Zeta - d) / Zeta);    // dimensionless
    const double ga2 = log((WindSpeedHeight + Zetam - d) / Zetam);  // dimensionless
    const double gbv = ga0 / (ga1 * ga2);                           // m / s

    // Apply the minimum
    return std::max(gbv, minimum_gbw);  // m / s
}

/* Soil Evaporation Function */
/* Variables I need */
/* LAI = Leaf Area Index */
/* k = extinction coefficient */
/* ppfd = Photosynthetic photon flux density */
/* awc, wiltp, fieldc = available water content, wilting point and field capacity */
/* winds = wind speed */

double SoilEvapo(
    double LAI, double k, double air_temperature, double ppfd,
    double soil_water_content, double fieldc, double wiltp, double winds,
    double RelH, double rsec, double soil_clod_size, double soil_reflectance,
    double soil_transmission, double specific_heat_of_air,
    double par_energy_content)
{
    int method = 1;
    /* A simple way of calculating the proportion of the soil that is hit by direct radiation. */
    double soil_area_sunlit_fraction = exp(-k * LAI);  // dimensionless.

    /* For now the temperature of the soil will be the same as the air.
       At a later time this can be made more accurate. I looked at the equations for this and the issue is that it is strongly dependent on
       depth. Since the soil model now has a single layer, this cannot be implemented correctly at the moment.  */
    double SoilTemp = air_temperature;

    /* From Campbell and Norman. Environmental Biophysics. */
    /* If relative available water content is */
    double rawc = (soil_water_content - wiltp) / (fieldc - wiltp);  // dimensionless. relative available water content.

    /* Campbell and Norman. Environmental Physics, page 142 */
    double maximum_uptake_rate = 1 - pow((1 + 1.3 * rawc), -5);  // dimensionless
    /* This is a useful idea because dry soils evaporate little water when dry*/

    /* Total Radiation */
    /* At the moment soil evaporation is grossly overestimated. In WIMOVAC the light reaching the last layer of leaves is used. Here instead
       of calculating this again, I will for now assume a 10% as a rough estimate. Note that I could maybe get this since layIdir and
       layIDiff in sunML are external variables.  Rprintf("ppfd %.5f",layIdir[0],"\n"); Update: 03-13-2009. I tried printing this
       value but it is still too high and will likely overestimate soil evaporation. However, this is still a work in progress.
       */
    ppfd *= rsec; /* Radiation soil evaporation coefficient */

    double TotalRadiation = ppfd * par_energy_content;

    double DdryA = TempToDdryA(air_temperature);
    double LHV = TempToLHV(air_temperature);  // J / kg
    double SlopeFS = TempToSFS(air_temperature);
    double SWVC = saturation_vapor_pressure(air_temperature) * 1e-5;

    double PsycParam = (DdryA * specific_heat_of_air) / LHV;
    double vapor_density_deficit = SWVC * (1 - RelH / 100);

    double BoundaryLayerThickness = 4e-3 * sqrt(soil_clod_size / winds);
    double DiffCoef = 2.126e-5 + 1.48e-7 * SoilTemp;
    double SoilBoundaryLayer = DiffCoef / BoundaryLayerThickness;

    double Ja = 2 * TotalRadiation * ((1 - soil_reflectance - soil_transmission) / (1 - soil_transmission));

    double rlc = 4 * physical_constants::stefan_boltzmann * pow((conversion_constants::celsius_to_kelvin + SoilTemp), 3) * 0.005;
    /* the last term should be the difference between air temperature and soil. This is not actually calculated at the moment. Since this is
       mostly relevant to the first soil layer where the temperatures are similar. I will leave it like this for now. */

    double PhiN = Ja - rlc; /* Calculate the net radiation balance*/
    if (PhiN < 0) PhiN = 1e-7;

    double Evaporation = 0.0;
    if (method == 0) {
        /* Priestly-Taylor */
        Evaporation = 1.26 * (SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam));  // kg / m^2 / s.
    } else {
        /* Penman-Monteith */
        Evaporation = (SlopeFS * PhiN + LHV * PsycParam * SoilBoundaryLayer * vapor_density_deficit) / (LHV * (SlopeFS + PsycParam));  // kg / m^2 / s.
    }

    Evaporation *= soil_area_sunlit_fraction * maximum_uptake_rate;  // kg / m^2 / s.
    if (Evaporation < 0) Evaporation = 1e-6;  // Prevent odd values at very low light levels.

    return Evaporation;  // kg / m^2 / s.
}

double compute_wsPhoto(int wsFun, double fieldc, double wiltp, double phi1, double awc) {
    // Three different type of equations for modeling the effect of
    // water stress on vmax and leaf area expansion.  The equation for
    // leaf area expansion is more severe than the one for vmax. */
    double wsPhoto;
    switch (wsFun) {
    case 0: { /* linear */
        double slp = 1 / (fieldc - wiltp);
        double intcpt = 1 - fieldc * slp;
        wsPhoto = slp * awc + intcpt;
        break;
    }
    case 1: {
        double phi10 = (fieldc + wiltp) / 2;
        wsPhoto = 1 / (1 + exp((phi10 - awc) / phi1));
        break;
    }
    case 2: {
        double slp = (1 - wiltp) / (fieldc - wiltp);
        double intcpt = 1 - fieldc * slp;
        double theta = slp * awc + intcpt;
        wsPhoto = (1 - exp(-2.5 * (theta - wiltp)/(1 - wiltp))) / (1 - exp(-2.5));
        break;
    }
    case 3:
        wsPhoto = 1;
        break;
    default:
        wsPhoto = 0;
    }

    // wsPhoto can be mathematically lower than zero in some cases but
    // I should prevent that:
    if (wsPhoto <= 0) {
        wsPhoto = 1e-10;
    }

    // Apparently wsPhoto can be greater than 1.
    if (wsPhoto > 1) wsPhoto = 1;

    return wsPhoto;
}

/* This is a new function that attempts to keep a water budget and then
   calcualte an empirical coefficient that reduces the specific leaf area.
   This results from the general idea that water stress reduces first the
   rate of leaf expansion. */

/* This is meant to be a simple function that calculates a
   simple empirical coefficient that reduces specifi leaf area
   according to the water stress of the plant. This is done
   for now, with a very simple empirical approach. */

struct ws_str watstr(double precipit, double evapo, double cws, double soildepth, double fieldc,
                     double wiltp, double soil_saturation_capacity, double sand,
                     double Ks, double air_entry, double b)
{
    constexpr double g = 9.8; // m / s^2  ##  http://en.wikipedia.org/wiki/Standard_gravity

    double precipM = precipit * 1e-3; /* convert precip in mm to m*/
    double soil_water_fraction = precipM / soildepth + cws;  // m^3 m^-3

    double runoff = 0.0;
    double Nleach = 0.0;  /* Nleach is the NO3 leached. */

    if (soil_water_fraction > soil_saturation_capacity) {
        runoff = (soil_water_fraction - soil_saturation_capacity) * soildepth; /* This is in meters */
        /* Here runoff is interpreted as water content exceeding saturation level */
        /* Need to convert to units used in the Parton et al 1988 paper. */
        /* The data come in mm/hr and need to be in cm/month */
        Nleach = runoff / (1e3 * physical_constants::molar_mass_of_water) * (0.2 + 0.7 * sand);
        soil_water_fraction = soil_saturation_capacity;
    }

    /* The density of water is 0.9982 Mg / m^3 at 20 degrees Celsius. */
    /* evapo is demanded water (Mg / ha) */
    double npaw = soil_water_fraction - wiltp - evapo / 0.9982 / 1e4 / soildepth;  // fraction of saturation.

    /* If demand exceeds supply, the crop is getting close to wilting point and transpiration
       will be over estimated. In this one layer model though, the crop is practically dead. */
    if (npaw < 0) npaw = 0.0;
    double awc = npaw + wiltp;

    /* Calculating the soil water potential based on equations from Norman and Campbell */
    /* tmp.psim = soTexS.air_entry * pow((awc/soTexS.fieldc*1.1),-soTexS.b); */
    /* New version of the soil water potential is based on
     * "Dynamic Simulation of Water Deficit Effects upon Maize Yield" R. F. Grant Agricultural Systems. 33(1990) 13-39. */
    ws_str tmp;
    tmp.psim = -exp(log(0.033) + (log(fieldc) - log(awc)) / (log(fieldc) - log(wiltp)) * (log(1.5) - log(0.033))) * 1e3; /* This last term converts from MPa to kPa */

    /* This is drainage */
    if (awc > fieldc) {
        double K_psim = Ks * pow((air_entry / tmp.psim), 2 + 3 / b); /* This is hydraulic conductivity */

        double J_w = -K_psim * (-tmp.psim / (soildepth * 0.5)) - g * K_psim; /*  Campbell, pg 129 do not ignore the graviational effect. I multiply soil depth by 0.5 to calculate the average depth */
        double drainage = J_w * 3600 * 0.9982 * 1e-3; /* This is flow in m^3 / (m^2 * hr). */
        awc = awc + drainage / soildepth;
    }

    tmp.awc = awc;
    tmp.runoff = runoff;
    tmp.Nleach = Nleach;
    return tmp;
}

/* Function to simulate the multilayer behavior of soil water. In the
   future this could be coupled with Campbell (BASIC) ideas to
   esitmate water potential. */
struct soilML_str soilML(double precipit, double transp, double *cws, double soildepth, double *depths,
        double soil_field_capacity, double soil_wilting_point, double soil_saturation_capacity, double soil_air_entry, double soil_saturated_conductivity,
        double soil_b_coefficient, double soil_sand_content, double phi1, double phi2, int wsFun,
        int layers, double rootDB, double LAI, double k, double AirTemp,
        double IRad, double winds, double RelH, int hydrDist, double rfl,
        double rsec, double rsdf, double soil_clod_size, double soil_reflectance, double soil_transmission,
        double specific_heat_of_air, double par_energy_content)
{
    constexpr double g = 9.8; /* m / s-2  ##  http://en.wikipedia.org/wiki/Standard_gravity */

    soilML_str return_value;

    /* Crude empirical relationship between root biomass and rooting depth*/
    double rootDepth = fmin(rootDB * rsdf, soildepth);

    rd_str root_distribution = rootDist(layers, rootDepth, &depths[0], rfl);

    /* unit conversion for precip */
    double oldWaterIn = 0.0;
    double waterIn = precipit * 1e-3; /* convert precip in mm to m */

    double drainage = 0.0;
    double wsPhotoCol = 0.0;
    double LeafWSCol = 0.0;
    double Sevap = 0.0;
    double oldEvapoTra = 0.0;

    for (int i = layers - 1; i >= 0; --i) { /* The counter, i, decreases because I increase the water content due to precipitation in the last layer first*/

        double layerDepth = depths[i+1] - depths[i];  /* This supports unequal depths. */

        if (hydrDist > 0) {
            /* For this section see Campbell and Norman "Environmental BioPhysics" Chapter 9*/
            /* First compute the matric potential */
            double psim1 = soil_air_entry * pow((cws[i]/soil_saturation_capacity), -soil_b_coefficient); /* This is matric potential of current layer */
            double dPsim;
            if (i > 0) {
                double psim2 = soil_air_entry * pow((cws[i-1]/soil_saturation_capacity), -soil_b_coefficient); /* This is matric potential of next layer */
                dPsim = psim1 - psim2;
                /* The substraction is from the layer i - (i-1). If this last term is positive then it will move upwards. If it is negative it will move downwards. Presumably this term is almost always positive. */
            } else {
                dPsim = 0.0;
            }

            double K_psim = soil_saturated_conductivity * pow((soil_air_entry/psim1), 2+3/soil_b_coefficient); /* This is hydraulic conductivity */
            double J_w = K_psim * (dPsim/layerDepth) - g * K_psim; /*  Campbell, pg 129 do not ignore the graviational effect*/
            /* Notice that K_psim is positive because my reference system is reversed */
            /* This last result should be in kg/(m2 * s) */
            J_w *= 3600 * 0.9882 * 1e-3; /* This is flow in m3 / (m^2 * hr). */

            if (i == (layers-1) && J_w < 0) {
                /* cws[i] = cws[i] + J_w / layerDepth; Although this should be done it drains the last layer too much.*/
                drainage += J_w;
            } else {
                    cws[i] = cws[i] - J_w / layerDepth;
                if (i > 0) {
                    cws[i - 1] = cws[i - 1] + J_w / layerDepth;
                }
            }
        }

        if (cws[i] > soil_saturation_capacity) cws[i] = soil_saturation_capacity;
        if (cws[i] < soil_wilting_point) cws[i] = soil_wilting_point;
        // Here is a convention aw is available water in volume and awc
        // is available water content as a fraction of the soil section being investigated.
        double aw = cws[i] * layerDepth;
        /* Available water (for this layer) is the current water status times the layer depth */

        if (waterIn > 0) {
            /* There is some rain. Need to add it.*/
            aw += waterIn / layers + oldWaterIn; /* They are both in meters so it works */
            /* Adding the same amount to water to each layer */
            /* In case there is overflow */
            double diffw = soil_field_capacity * layerDepth - aw;

            if (diffw < 0) {
                /* This means that precipitation exceeded the capacity of the first layer */
                /* Save this amount of water for the next layer */
                oldWaterIn = -diffw;
                aw = soil_field_capacity * layerDepth;
            } else {
                oldWaterIn = 0.0;
            }
        }

        /* Root Biomass */
        double rootATdepth = rootDB * root_distribution.rootDist[i];
        return_value.rootDist[i] = rootATdepth;
        /* Plant available water is only between current water status and permanent wilting point */
        /* Plant available water */
        double pawha = (aw - soil_wilting_point * layerDepth) * 1e4;

        if (pawha < 0) pawha = 0;

        double Ctransp = 0.0;
        double EvapoTra = 0.0;
        double Newpawha;

        if (i == 0) {
            /* Only the first layer is affected by soil evaporation */
            double awc2 = aw / layerDepth;
            /* SoilEvapo function needs soil water content  */
            Sevap = SoilEvapo(
                LAI, k, AirTemp, IRad, awc2, soil_field_capacity,
                soil_wilting_point, winds, RelH, rsec, soil_clod_size,
                soil_reflectance, soil_transmission, specific_heat_of_air,
                par_energy_content) * 3600 * 1e-3 * 10000;  // Mg / ha / hr. 3600 s / hr * 1e-3 Mg / kg * 10000 m^2 / ha.
            /* I assume that crop transpiration is distributed simlarly to
               root density.  In other words the crop takes up water proportionally
               to the amount of root in each respective layer.*/
            Ctransp = transp * root_distribution.rootDist[0];
            EvapoTra = Ctransp + Sevap;
            constexpr double density_of_water_at_20_celcius = 0.9982;  // Mg m^-3.
            Newpawha = pawha - EvapoTra / density_of_water_at_20_celcius;
            /* The first term in the rhs pawha is the m3 of water available in this layer.
               EvapoTra is the Mg H2O ha-1 of transpired and evaporated water. 1/0.9882 converts from Mg to m3 */
        } else {
            Ctransp = transp * root_distribution.rootDist[i];
            EvapoTra = Ctransp;
            Newpawha = pawha - (EvapoTra + oldEvapoTra);
        }

        if (Newpawha < 0) {
            /* If the Demand is not satisfied by this layer. This will be stored and added to subsequent layers*/
            oldEvapoTra = -Newpawha;
            aw = soil_wilting_point * layerDepth;
        }

        double awc = Newpawha / 1e4 / layerDepth + soil_wilting_point;

        /* This might look like a weird place to populate the structure, but is more convenient*/
        return_value.cws[i] = awc;

        // To-do: Replace this block with a call to compute_wsPhoto.
    /* three different type of equations for modeling the effect of water stress on vmax and leaf area expansion.
       The equation for leaf area expansion is more severe than the one for vmax. */
        double wsPhoto = 0.0;
        double slp = 0.0;
        double intcpt = 0.0;
        double theta = 0.0;

        if (wsFun == 0) { /* linear */
            slp = 1/(soil_field_capacity - soil_wilting_point);
            intcpt = 1 - soil_field_capacity * slp;
            wsPhoto = slp * awc + intcpt;
        } else if (wsFun == 1) {
            double phi10 = (soil_field_capacity + soil_wilting_point)/2;
            wsPhoto = 1/(1 + exp((phi10 - awc)/ phi1));
        } else if (wsFun == 2) {
            slp = (1 - soil_wilting_point)/(soil_field_capacity - soil_wilting_point);
            intcpt = 1 - soil_field_capacity * slp;
            theta = slp * awc + intcpt;
            wsPhoto = (1 - exp(-2.5 * (theta - soil_wilting_point)/(1 - soil_wilting_point))) / (1 - exp(-2.5));
        } else if (wsFun == 3) {
            wsPhoto = 1;
        }

        if (wsPhoto <= 0 )
            wsPhoto = 1e-10; /* This can be mathematically lower than zero in some cases but I should prevent that. */

        wsPhotoCol += wsPhoto;

        double LeafWS = pow(awc, phi2) * 1/pow(soil_field_capacity, phi2);
        if (wsFun == 3) {
            LeafWS = 1;
        }
        LeafWSCol += LeafWS;

    }

    if (waterIn > 0) {
        drainage = waterIn;
        /* Need to convert to units used in the Parton et al 1988 paper. */
        /* The data comes in mm/hr and it needs to be in cm/month */
        return_value.Nleach = drainage * 0.1 * (1/24 * 30) / (1e3 * physical_constants::molar_mass_of_water * (0.2 + 0.7 * soil_sand_content));
    }
    else {
        return_value.Nleach = 0.0;
    }


    return_value.rcoefPhoto = (wsPhotoCol/layers);
    return_value.drainage = drainage;
    return_value.rcoefSpleaf = (LeafWSCol/layers);
    return_value.SoilEvapo = Sevap;

    return return_value;
}

/**
 *  @brief Subtracts respiratory losses from a carbon production rate
 *
 *  @param [in] base_rate The base rate of carbon production that does not
 *                        include respiratory losses. Any units are acceptable,
 *                        e.g. mol / m^2 / s or Mg / ha / hour.
 *
 *  @param [in] mrc Maintenance respiration coefficient (dimensionless)
 *
 *  @param [in] temp Temperature (degrees C)
 *
 *  @return A modified rate having the same units as `base_rate` where
 *          respiratory losses have been subtracted
 *
 *  ### Model overview
 *
 *  The idea here is that the net carbon assimilation rate (`A_n`) is given by
 *
 *  > `A_n = A_g - R` (1)
 *
 *  where `A_g` is the gross assimilation rate or the production rate and `R` is
 *  the respiration rate. A further assumption is that respiration is
 *  proportional to production, i.e.
 *
 *  > `R = R_c * A_g` (2)
 *
 *  where `R_c` is a dimensionless coefficient that lies on the interval [0,1].
 *  Combining equations (1) and (2), we can express the net assimilation rate
 *  as a function of the gross assimilation rate and the respiration
 *  coefficient:
 *
 *  > `A_n = A_g - R_c * A_g =  A_g * (1 - R_c)` (3)
 *
 *  In this function, the temperature dependence of the respiration coefficient
 *  is modeled using a simple "Q10" method where `Q10 = 2` and the base
 *  temperature is 0 degrees C, i.e.,
 *
 *  > `R_c = R_c_0 * 2^(T / 10)` (4)
 *
 *  where `T` is the temperature in Celsius. Combining equations (3) and (4), we
 *  finally arrive at a relationship between net assimilation, gross
 *  assimilation, and temperature:
 *
 *  > `A_n = A_g * (1 - R_c_0 * 2^(T / 10))` (5)
 *
 *  In the code below, `base_rate` represents the gross assimilation rate `A_g`,
 *  `mrc` represents the maintenance respiration coefficient `R_c_0`, and `temp`
 *  represents the temperature `T`.
 *
 *  ### Sources
 *
 *  In Stephen Humphries's thesis, he describes the respiration model in the
 *  following way: "The respiration model of McCree (1970) modified according to
 *  Penning de Vries (1972) and Thornley (1970) is used here to predict plant
 *  respiration... The respiration associated with each plant structure is
 *  modified according to the temperature of the structure using a Q10
 *  approximation with a value of 2 as described by Spain and Keen (1992)."
 *
 *  I (EBL) believe these references to be to the following documents:
 *
 *  - [de Vries, F. P. "Respiration and growth" in "Crop processes in controlled
 *    environments" 327–347 (Academic Press, 1972)]
 *    (https://library.wur.nl/WebQuery/wurpubs/fulltext/218533)
 *
 *  - [Thornley, J. H. M. "Respiration, Growth and Maintenance in Plants" Nature
 *    227, 304–305 (1970)](https://doi.org/10.1038/227304b0)
 *
 *  - McCree, K. "An equation for the rate of respiration of whiteclover plants
 *    grown under controlled conditions" in: "Prediction and measurement
 *    of photosynthetic productivity" Wageningen: Centre for Agricultural
 *    Publishing and Documentation, 221-229 (1970)
 *
 *  - Keen, R. E. & Spain, J. D. "Computer simulation in biology. A BASIC
 *    introduction" (1992)
 *
 *  I believe the McCree, Thornley, and de Vries papers describe the
 *  relationship between gross assimilation and respiration, while the Keen &
 *  Spain book describes the temperature dependence. Unfortunately, I can't find
 *  an online version of that book so I can't be sure.
 *
 *  For some general discussions about respiration, see the following two
 *  sources:
 *
 *  - [Amthor, J. S. "The McCree–de Wit–Penning de Vries–Thornley Respiration
 *    Paradigms: 30 Years Later" Ann Bot 86, 1–20 (2000)]
 *    (https://doi.org/10.1006/anbo.2000.1175)
 *
 *  - [Amthor, J. S. "The role of maintenance respiration in plant growth."
 *    Plant, Cell & Environment 7, 561–569 (1984)]
 *    (https://doi.org/10.1111/1365-3040.ep11591833)
 *
 *  The Humphries thesis is also unfortunately not available online:
 *
 *  Humphries, S. "Will mechanistically rich models provide us with new insights
 *  into the response of plant production to climate change?: development and
 *  experiments with WIMOVAC: (Windows Intuitive Model of Vegetation response
 *  to Atmosphere & Climate Change)" (University of Essex, 2002)
 */
double resp(double base_rate, double mrc, double temp) {
    double ans = base_rate *  (1 - (mrc * pow(2, (temp / 10.0))));

    if (ans <0) ans = 0;

    return ans;

}

struct seqRD_str seqRootDepth(double to, int lengthOut ) {
    double by = to / lengthOut;

    struct seqRD_str result;
    for (int i = 0; i <= lengthOut; ++i) {
        result.rootDepths[i] = i * by;
    }
    return result;
}

struct rd_str rootDist(int n_layers, double rootDepth, double *depths, double rfl)
{
    /*
     * Calculate the fraction of total root mass for each layer in `depths` assuming the mass
     * is follows a Poisson distribution along the depth.
     *
     * Returns a struct with one array of size `n_layers`.
     * Each element in the array is the fraction of total root mass in that layer.
     * The sum of all elements of the result equals 1.
     */

    double layerDepth = 0.0;
    double CumLayerDepth = 0.0;
    int CumRootDist = 1;
    double rootDist[n_layers];
    double cumulative_a = 0.0;

    for (int i = 0; i < n_layers; ++i) {
        if (i == 0) {
            layerDepth = depths[1];
        } else {
            layerDepth = depths[i] - depths[i-1];
        }

        CumLayerDepth += layerDepth;

        if (rootDepth > CumLayerDepth) {
            CumRootDist++;
        }
    }

    for (int j = 0; j < n_layers; ++j) {
        if (j < CumRootDist) {
            double a = poisson_density(j + 1, (double)CumRootDist * rfl);
            rootDist[j] = a;
            cumulative_a += a;
        } else {
            rootDist[j] = 0;
        }
    }

    struct rd_str result;
    for (int k = 0; k < n_layers; ++k) {
        result.rootDist[k] = rootDist[k] / cumulative_a;
    }
    return  result;
}
