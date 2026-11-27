#include "PrescribedControl.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> values(const json::Value& source) {
    std::vector<double> result;
    for (const auto& item : source.array()) result.push_back(item.number());
    return result;
}

std::vector<double> vectorOrScalar(const json::Value& settings,
                                   const char* plural, const char* singular,
                                   double fallback) {
    auto result = values(settings.at(plural));
    if (result.empty()) result.push_back(settings.at(singular).number(fallback));
    return result;
}

void requireFinite(const std::vector<double>& items, const char* name) {
    for (const double item : items)
        if (!std::isfinite(item))
            throw std::invalid_argument(std::string(name) + " must be finite");
}

std::uint64_t seed(const json::Value& source) {
    const double value = source.number();
    if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value
        || value > 9007199254740991.0)
        throw std::invalid_argument(
            "controller seed must be an exactly representable nonnegative integer");
    return static_cast<std::uint64_t>(value);
}

std::vector<double> window(const json::Value& settings,
                           const dynsyslab::core::CoordinateScaling& scaling) {
    bool physical = !settings.at("active_time_window_seconds").isNull();
    bool model = !settings.at("model_active_time_window").isNull();
    if (physical && model)
        throw std::invalid_argument("active window cannot mix physical and model time");
    auto result = values(physical ? settings.at("active_time_window_seconds")
                                  : settings.at("model_active_time_window"));
    if (!physical && !model) {
        result = values(settings.at("active_time_window"));
        physical = scaling.physicalTimeAvailable;
    }
    if (result.empty()) return result;
    if (result.size() != 2 || !std::isfinite(result[0]) || !std::isfinite(result[1])
        || result[0] > result[1])
        throw std::invalid_argument("active time window must be [start, stop]");
    if (physical) {
        result[0] = scaling.modelTime(result[0]);
        result[1] = scaling.modelTime(result[1]);
    }
    return result;
}

double modelAmplitude(const json::Value& settings,
                      const dynsyslab::core::CoordinateScaling& scaling,
                      double fallback, bool& physical) {
    if (!settings.at("model_amplitude").isNull()) {
        physical = false;
        return settings.at("model_amplitude").number(fallback);
    }
    physical = scaling.physicalInputAvailable;
    const double value = settings.at("amplitude").number(fallback);
    return physical ? scaling.modelInput(value) : value;
}

PrescribedControl waveform(const json::Value& settings,
                           const dynsyslab::core::CoordinateScaling& scaling) {
    const bool hasPhysicalFrequency = !settings.at("frequencies_hz").isNull()
        || !settings.at("frequency_hz").isNull();
    const bool hasModelFrequency = !settings.at("frequencies").isNull()
        || !settings.at("frequency").isNull();
    if (hasPhysicalFrequency && hasModelFrequency)
        throw std::invalid_argument("waveform cannot mix physical and model frequencies");
    const bool physicalFrequency = hasPhysicalFrequency
        || (!hasModelFrequency && scaling.physicalTimeAvailable);
    auto frequencies = physicalFrequency
        ? vectorOrScalar(settings, "frequencies_hz", "frequency_hz", 1.0)
        : vectorOrScalar(settings, "frequencies", "frequency", 1.0);

    const bool hasModelAmplitudes = !settings.at("model_amplitudes").isNull()
        || !settings.at("model_amplitude").isNull();
    const bool hasPhysicalAmplitudes = !settings.at("amplitudes").isNull()
        || !settings.at("amplitude").isNull();
    if (hasModelAmplitudes && hasPhysicalAmplitudes)
        throw std::invalid_argument("waveform cannot mix physical and model amplitudes");
    const bool physicalAmplitude = !hasModelAmplitudes && scaling.physicalInputAvailable;
    auto amplitudes = hasModelAmplitudes
        ? vectorOrScalar(settings, "model_amplitudes", "model_amplitude", 1.0)
        : vectorOrScalar(settings, "amplitudes", "amplitude", 1.0);
    if (frequencies.size() != amplitudes.size())
        throw std::invalid_argument("waveform frequencies and amplitudes must have equal size");
    if (physicalAmplitude)
        for (double& amplitude : amplitudes) amplitude = scaling.modelInput(amplitude);

    auto phases = values(settings.at("phases_rad"));
    if (phases.empty()) phases.assign(frequencies.size(),
        settings.at("phase_radians").number(settings.at("phase").number()));
    if (phases.size() != frequencies.size())
        throw std::invalid_argument("waveform phases must match its frequencies");
    requireFinite(frequencies, "waveform frequencies");
    requireFinite(amplitudes, "waveform amplitudes");
    requireFinite(phases, "waveform phases");

    const bool modelOffset = !settings.at("model_offset").isNull();
    const bool physicalOffset = !settings.at("offset").isNull()
        && scaling.physicalInputAvailable;
    const double rawOffset = modelOffset ? settings.at("model_offset").number()
                                         : settings.at("offset").number();
    const double offset = !modelOffset && scaling.physicalInputAvailable
        ? scaling.modelInput(rawOffset) : rawOffset;
    if (!std::isfinite(offset)) throw std::invalid_argument("waveform offset must be finite");
    constexpr double twoPi = 6.283185307179586476925286766559;
    PrescribedControl result;
    result.law = [frequencies = std::move(frequencies), amplitudes = std::move(amplitudes),
                  phases = std::move(phases), offset, scaling, physicalFrequency]
        (std::size_t, double time, const double*) {
        const double clock = physicalFrequency ? scaling.physicalTime(time) : time;
        double value = offset;
        for (std::size_t index = 0; index < frequencies.size(); ++index)
            value += amplitudes[index]
                * std::sin(twoPi * frequencies[index] * clock + phases[index]);
        return value;
    };
    result.activeWindow = window(settings, scaling);
    const bool physicalInput = physicalAmplitude || physicalOffset;
    result.settingsCoordinates = physicalFrequency == physicalInput
        ? (physicalFrequency ? "physical" : "model") : "mixed";
    return result;
}

PrescribedControl randomPulse(const json::Value& settings,
                              const dynsyslab::core::CoordinateScaling& scaling) {
    bool physicalAmplitude = false;
    const double amplitude = modelAmplitude(settings, scaling, 1.0, physicalAmplitude);
    const bool explicitPhysical = !settings.at("pulse_duration_seconds").isNull();
    const bool explicitModel = !settings.at("model_pulse_duration").isNull();
    if (explicitPhysical && explicitModel)
        throw std::invalid_argument("pulse duration cannot mix physical and model time");
    bool physicalDuration = explicitPhysical;
    double duration{};
    if (explicitPhysical) duration = settings.at("pulse_duration_seconds").number();
    else if (explicitModel) duration = settings.at("model_pulse_duration").number();
    else {
        duration = settings.at("pulse_duration").number(
            settings.at("interval").number(0.1));
        physicalDuration = scaling.physicalTimeAvailable;
    }
    if (physicalDuration) duration = scaling.modelTime(duration);
    const auto randomSeed = seed(settings.at("seed"));
    if (!std::isfinite(amplitude) || !std::isfinite(duration) || duration <= 0.0)
        throw std::invalid_argument(
            "random pulse amplitude must be finite and duration positive");
    PrescribedControl result;
    result.law = [amplitude, duration, randomSeed](std::size_t trajectory, double time,
                                                   const double*) {
        std::uint64_t value = randomSeed ^ (trajectory + 0x9e3779b97f4a7c15ull)
            ^ (static_cast<std::uint64_t>(std::floor(time / duration))
                * 0xbf58476d1ce4e5b9ull);
        value ^= value >> 30; value *= 0xbf58476d1ce4e5b9ull;
        value ^= value >> 27; value *= 0x94d049bb133111ebull;
        value ^= value >> 31;
        return amplitude * (2.0 * (value >> 11) / 9007199254740992.0 - 1.0);
    };
    result.activeWindow = window(settings, scaling);
    result.settingsCoordinates = physicalAmplitude == physicalDuration
        ? (physicalAmplitude ? "physical" : "model") : "mixed";
    return result;
}

}  // namespace

PrescribedControl buildPrescribedControl(
    const std::string& kind, const json::Value& settings,
    const dynsyslab::core::ODESystem& system) {
    const auto scaling = system.coordinates();
    if (kind == "waveform") return waveform(settings, scaling);
    if (kind == "random_pulse") return randomPulse(settings, scaling);
    auto result = buildDesignedInput(kind, settings, scaling);
    result.activeWindow = window(settings, scaling);
    return result;
}

bool isPrescribedControl(const std::string& kind) {
    return kind == "waveform" || kind == "random_pulse" || kind == "tabulated"
        || kind == "command_tape" || kind == "telegraph_multisine" || kind == "symmetric_pulse" || kind == "pulse_branches";
}

}  // namespace dynsyslab::workflow::detail
