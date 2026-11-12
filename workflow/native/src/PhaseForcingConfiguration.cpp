#include "PhaseForcingConfiguration.hpp"
#include "PrescribedControl.hpp"
#include <stdexcept>
namespace dynsyslab::workflow::detail {
namespace {
double single(const json::Value &s, const char *plural, const char *singular, double fallback) {
    const auto &items = s.at(plural).array();
    if (items.empty())
        return s.at(singular).number(fallback);
    if (items.size() != 1)
        throw std::invalid_argument(
            "Periodic archive sidecar currently requires a single sinusoidal controller tone");
    return items[0].number();
}
} // namespace
dynsyslab::core::SinusoidalDriver
phaseDriverFromController(const json::Value &controller, const std::string &evaluation,
                          const dynsyslab::core::ODESystem &system) {
    if (controller.at("kind").string() != "waveform" ||
        (evaluation != "continuous_rhs" && evaluation != "automatic" && !evaluation.empty()))
        throw std::invalid_argument(
            "Periodic archive sidecar requires a continuous waveform controller");
    const auto &s = controller.at("settings");
    // Use the production controller validator before translating its explicit functional form.
    (void)buildPrescribedControl("waveform", s, system);
    for (const char *key :
         {"active_time_window", "model_active_time_window", "active_time_window_seconds"})
        if (!s.at(key).array().empty())
            throw std::invalid_argument("A time-windowed waveform is not globally periodic");
    if (s.at("offset").number() != 0 || s.at("model_offset").number() != 0)
        throw std::invalid_argument(
            "Periodic archive sidecar does not yet support a waveform offset");
    const auto coordinates = system.coordinates();
    const bool modelAmplitude =
        !s.at("model_amplitudes").isNull() || !s.at("model_amplitude").isNull();
    double amplitude = modelAmplitude ? single(s, "model_amplitudes", "model_amplitude", 1.)
                                      : single(s, "amplitudes", "amplitude", 1.);
    if (modelAmplitude && coordinates.physicalInputAvailable)
        amplitude = coordinates.physicalInput(amplitude);
    const bool physicalFrequency = !s.at("frequencies_hz").isNull() ||
                                   !s.at("frequency_hz").isNull() ||
                                   ((s.at("frequencies").isNull() && s.at("frequency").isNull()) &&
                                    coordinates.physicalTimeAvailable);
    double frequency = physicalFrequency ? single(s, "frequencies_hz", "frequency_hz", 1.)
                                         : single(s, "frequencies", "frequency", 1.);
    if (!physicalFrequency && coordinates.physicalTimeAvailable)
        frequency /= coordinates.secondsPerModelTime;
    const double phase = single(s, "phases_rad", "phase_radians", s.at("phase").number());
    // PrescribedControl uses sin; System::SinusoidalDriver uses cos.
    return {amplitude, frequency, phase - 1.5707963267948966192313216916398};
}
} // namespace dynsyslab::workflow::detail
