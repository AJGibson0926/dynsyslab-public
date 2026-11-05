#pragma once

#include "SystemSupport.hpp"

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {

struct BubbleMedium {
    double density{998.0};
    double viscosity{1.0e-3};
    double surfaceTension{0.0725};
    double vaporPressure{2330.0};
    double kappa{4.0 / 3.0};
};

inline BubbleMedium bubbleMedium(SystemConfiguration& configuration,
                                 BubbleMedium fallback = {}) {
    auto medium = fallback;
    const auto& preset = configuration.preset;
    if (preset == "air_water_0c")
        medium = {999.8, 1.792e-3, 0.0756, 611.0, 1.4};
    else if (preset == "air_water_37c")
        medium = {993.3, 0.691e-3, 0.0696, 6280.0, 1.4};
    else if (preset == "air_saline_37c")
        medium = {1004.6, 0.75e-3, 0.0730, 6280.0, 1.4};
    else if (preset == "air_blood_37c")
        medium = {1025.0, 3.5e-3, 0.055, 6280.0, 1.4};
    else if (preset == "air_glycerol_20c")
        medium = {1261.0, 1.49, 0.0634, 0.01, 1.4};
    else if (preset == "air_glycerol_water_50pct_20c")
        medium = {1126.0, 0.006, 0.067, 1800.0, 1.4};
    else if (preset == "air_silicone_oil_100cst_25c")
        medium = {965.0, 0.0965, 0.020, 1.0, 1.4};
    else if (preset == "co2_water_20c")
        medium = {998.2, 1.002e-3, 0.0728, 2339.0, 1.30};
    else if (preset == "co2_honey_25c")
        medium = {1420.0, 10.0, 0.055, 3200.0, 1.30};
    else if (preset == "marmottant_matlab_reference")
        medium = {1000.0, 1.0e-3, 0.0725, 0.0, 1.095};
    medium.density = resolvedParameter(configuration, "rho_kg_m3", medium.density);
    medium.viscosity = resolvedParameter(configuration, "mu_pa_s", medium.viscosity);
    medium.surfaceTension = resolvedParameter(configuration, "sigma_n_m", medium.surfaceTension);
    medium.vaporPressure = resolvedParameter(configuration, "vapor_pressure_pa", medium.vaporPressure);
    medium.kappa = resolvedParameter(configuration, "kappa", medium.kappa);
    return medium;
}

struct BubbleScaling {
    double radius{};
    double ambientPressure{};
    double gasPressure{};
    double naturalFrequency{};
    double timeScale{};
    double pressureScale{};
};

inline BubbleScaling bubbleScaling(const BubbleMedium& medium, double radius,
                                   double ambientPressure, double frequencyFactor = 1.0) {
    const double gasPressure = ambientPressure - medium.vaporPressure
        + 2.0 * medium.surfaceTension / radius;
    const double stiffness = 3.0 * medium.kappa * gasPressure
        - 2.0 * medium.surfaceTension / radius;
    if (!(radius > 0.0) || !(ambientPressure > 0.0) || !(stiffness > 0.0))
        throw std::invalid_argument("bubble scaling requires positive radius, pressure, and stiffness");
    constexpr double twoPi = 6.283185307179586476925286766559;
    const double frequency = frequencyFactor * std::sqrt(stiffness)
        / (twoPi * radius * std::sqrt(medium.density));
    const double timeScale = 1.0 / frequency;
    return {radius, ambientPressure, gasPressure, frequency, timeScale,
            medium.density * std::pow(radius / timeScale, 2.0)};
}

std::unique_ptr<ODESystem> makeRPE(SystemConfiguration configuration);
std::unique_ptr<ODESystem> makeMarmottant(const SystemConfiguration& configuration);
std::unique_ptr<ODESystem> makeReddySzeri(const SystemConfiguration& configuration);

}  // namespace dynsyslab::core::detail
