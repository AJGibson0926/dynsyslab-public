#include <algorithm>
#include <cmath>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <stdexcept>
namespace dynsyslab::core {
double phaseSpaceForcingPeriod(const ODESystem &system, double requested) {
    const auto &c = system.configuration();
    constexpr double tau = 6.283185307179586476925286766559;
    std::vector<double> periods;
    if (!std::isfinite(c.driver.amplitude) || !std::isfinite(c.driver.frequencyHz) ||
        !std::isfinite(c.driver.phaseRadians) || !std::isfinite(requested) || requested < 0)
        throw std::invalid_argument("Periodic forcing requires finite driver values and period");
    if (c.driver.amplitude != 0 && c.driver.frequencyHz != 0) {
        const auto coordinates = system.coordinates();
        const double scale =
            coordinates.physicalTimeAvailable ? coordinates.secondsPerModelTime : 1.;
        periods.push_back(1 / std::abs(c.driver.frequencyHz * scale));
    }
    if (c.identifier == "reddy_szeri") {
        const auto &g = c.dimensionlessGroups;
        if (c.parameters.at("translation_enabled") != 0 && g.at("angular_frequency") != 0 &&
            (g.at("pressure_amplitude") != 0 || g.at("pressure_derivative_amplitude") != 0))
            periods.push_back(tau / std::abs(g.at("angular_frequency")));
    }
    if (periods.empty())
        throw std::invalid_argument("Periodic sidecar requires an active, known sinusoidal clock");
    auto common = [&](double t) {
        if (!std::isfinite(t) || t <= 0)
            return false;
        for (double p : periods) {
            const double cycles = t / p;
            if (!std::isfinite(p) || p <= 0 || cycles < .5 || cycles > 100000 ||
                std::abs(cycles - std::round(cycles)) > 1e-10 * std::max(1., cycles))
                return false;
        }
        return true;
    };
    if (requested > 0) {
        if (!common(requested))
            throw std::invalid_argument(
                "Forcing period must contain an integer number of every active forcing cycle");
        return requested;
    }
    const double longest = *std::max_element(periods.begin(), periods.end());
    for (int multiple = 1; multiple <= 64; ++multiple)
        if (common(multiple * longest))
            return multiple * longest;
    throw std::invalid_argument("No common forcing period found within 64 cycles; specify a "
                                "verified common period or use a quasiperiodic analysis");
}
} // namespace dynsyslab::core
