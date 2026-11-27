#include <dynsyslab/core/Targets.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<double> HarmonicStateTarget::evaluate(double time) const {
    if (amplitudes.empty() || amplitudes.size() != angularFrequencies.size()
        || (!phases.empty() && phases.size() != amplitudes.size()))
        throw std::logic_error("harmonic target parameters are inconsistent");
    double state = offset, velocity = 0.0;
    for (std::size_t index = 0; index < amplitudes.size(); ++index) {
        const double angle = angularFrequencies[index] * time
            + (phases.empty() ? 0.0 : phases[index]);
        state += amplitudes[index] * std::sin(angle);
        velocity += amplitudes[index] * angularFrequencies[index] * std::cos(angle);
    }
    return {state, velocity};
}

std::pair<std::vector<double>, std::vector<double>> HarmonicStateTarget::conservativeEnvelope() const {
    double radius = 0.0, velocity = 0.0;
    for (std::size_t index = 0; index < amplitudes.size(); ++index) {
        radius += std::abs(amplitudes[index]);
        velocity += std::abs(amplitudes[index] * angularFrequencies[index]);
    }
    return {{offset - radius, -velocity}, {offset + radius, velocity}};
}

const std::vector<HarmonicStateTarget>& rpeControlTargets() {
    constexpr double pi = 3.141592653589793238462643383279502884;
    static const std::vector<HarmonicStateTarget> targets{
        {"rpe_multisine_equilibrium_small", 1.0, {0.06,0.08}, {4.0,3.0}, {0.0,pi}},
        {"rpe_multisine_equilibrium_moderate_identification", 1.0,
            {0.09,0.12}, {4.0,3.0}, {0.0,pi}},
        {"rpe_multisine_offset_3p5_periodic", 3.5,
            {1.5,1.0/3.0,1.0}, {pi,2*pi,4*pi}, {}},
        {"rpe_multisine_offset_4_mixed", 4.0, {1.5,1.0}, {4.0,3.0}, {}},
        {"rpe_multisine_offset_4_periodic", 4.0,
            {1.5,1.0/3.0,1.0}, {pi,2*pi,4*pi}, {}},
        {"rpe_multisine_offset_2p85_mixed", 2.85, {1.5,1.0}, {4.0,3.0}, {}},
        {"rpe_single_sine_offset_2_experimental", 2.0, {1.5}, {4.59333}, {}}
    };
    return targets;
}

const HarmonicStateTarget& rpeControlTarget(const std::string& name) {
    for (const auto& target : rpeControlTargets()) if (target.name == name) return target;
    throw std::out_of_range("unknown canonical RPE control target: " + name);
}

}  // namespace dynsyslab::core
