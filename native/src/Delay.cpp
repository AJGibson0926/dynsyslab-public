#include <dynsyslab/core/Delay.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

double parameter(const SystemConfiguration& configuration, const char* name,
                 double fallback) {
    const auto found = configuration.parameters.find(name);
    return found == configuration.parameters.end() ? fallback : found->second;
}

std::size_t roundedCount(double value, const char* label) {
    if (!std::isfinite(value) || value < 0.0
        || value > static_cast<double>((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument(std::string(label) + " is not addressable");
    return static_cast<std::size_t>(std::llround(value));
}

std::vector<SystemDescriptor> delayCatalog() {
    return {mackeyGlassDescriptor()};
}

const bool registered = [] {
    registerSystemCatalogProvider(&delayCatalog);
    return true;
}();

}  // namespace

const SystemDescriptor& mackeyGlassDescriptor() {
    static const SystemDescriptor value = [] {
      SystemDescriptor result{
        "mackey_glass", "Mackey-Glass delay system", 1, "dimensional",
        {"x"}, {"x"}, {{"classic_chaos", "Classic chaotic delay dynamics"}},
        "dde"};
      result.research = mackeyGlassKnowledge();
      result.analyticSolutions = mackeyGlassAnalyticSolutions;
      return result;
    }();
    return value;
}

SystemConfiguration resolveMackeyGlass(const SystemConfiguration& requested) {
    SystemConfiguration result = requested;
    result.identifier = "mackey_glass";
    if (result.preset.empty()) result.preset = "classic_chaos";
    const std::pair<const char*, double> defaults[]{
        {"beta", 0.2}, {"gamma", 0.1}, {"delay", 17.0},
        {"exponent", 10.0}, {"history_value", 1.2}};
    for (const auto& item : defaults)
        if (result.parameters.find(item.first) == result.parameters.end())
            result.parameters[item.first] = item.second;
    result.parameterUnits = {{"beta", "1/time"}, {"gamma", "1/time"},
        {"delay", "time"}, {"exponent", "dimensionless"},
        {"history_value", "state"}};
    result.equationsVersion = "dynsyslab.mackey_glass/1";
    result.initialState = {result.parameters.at("history_value")};
    return result;
}

DelaySolution integrateMackeyGlass(
    const SystemConfiguration& requested, const std::vector<double>& historyValues,
    const DelayIntegrationOptions& options, const ExecutionContext* execution) {
    static_cast<void>(registered);
    const auto configuration = resolveMackeyGlass(requested);
    const double beta = parameter(configuration, "beta", 0.2);
    const double gamma = parameter(configuration, "gamma", 0.1);
    const double delay = parameter(configuration, "delay", 17.0);
    const double exponent = parameter(configuration, "exponent", 10.0);
    if (historyValues.empty()) throw std::invalid_argument("delay solve needs history values");
    if (!std::isfinite(options.initialTime) || !std::isfinite(options.finalTime)
        || options.finalTime <= options.initialTime || !std::isfinite(options.step)
        || options.step <= 0.0 || !std::isfinite(beta) || !std::isfinite(gamma)
        || !std::isfinite(delay) || delay < options.step || !std::isfinite(exponent))
        throw std::invalid_argument("invalid Mackey-Glass parameters or integration interval");
    const std::size_t steps = roundedCount(
        (options.finalTime - options.initialTime) / options.step, "delay step count");
    const std::size_t lag = roundedCount(delay / options.step, "delay lag count");
    if (steps > (std::numeric_limits<std::size_t>::max)() - lag - 1)
        throw ExecutionLimitExceeded("delay history buffer is not addressable");
    const std::size_t samples = steps + 1;
    if (historyValues.size() > (std::numeric_limits<std::size_t>::max)() / samples)
        throw ExecutionLimitExceeded("delay solution buffer is not addressable");
    if (execution) {
        const long double values = static_cast<long double>(historyValues.size())
            * static_cast<long double>(steps + lag + 1)
            + static_cast<long double>(samples)
            + static_cast<long double>(historyValues.size()) * samples;
        if (values > static_cast<long double>((std::numeric_limits<std::size_t>::max)())
            / sizeof(double))
            throw ExecutionLimitExceeded("delay solution buffer is not addressable");
        execution->requireMemory(static_cast<std::size_t>(values) * sizeof(double),
                                 "delay solution buffers");
    }
    DelaySolution result;
    result.trajectoryCount = historyValues.size();
    result.sampleCount = samples;
    result.acceptedSteps = steps * historyValues.size();
    result.time.resize(samples);
    result.values.resize(historyValues.size() * samples);
    for (std::size_t sample = 0; sample < samples; ++sample)
        result.time[sample] = options.initialTime + sample * options.step;
    for (std::size_t trajectory = 0; trajectory < historyValues.size(); ++trajectory) {
        if (!std::isfinite(historyValues[trajectory]))
            throw std::invalid_argument("Mackey-Glass history value must be finite");
        std::vector<double> history(steps + lag + 1, historyValues[trajectory]);
        for (std::size_t step = 0; step < steps; ++step) {
            if (execution) execution->checkpoint();
            const std::size_t active = lag + step;
            const double delayed = history[active - lag];
            history[active + 1] = history[active] + options.step
                * (beta * delayed / (1.0 + std::pow(delayed, exponent))
                   - gamma * history[active]);
            if (options.checkpoint && options.checkpointSteps && (step + 1) % options.checkpointSteps == 0
                && step + 1 < steps)
                options.checkpoint(trajectory, step + 2, result.time.data(), history.data() + lag);
        }
        std::copy(history.begin() + static_cast<std::ptrdiff_t>(lag), history.end(),
                  result.values.begin() + static_cast<std::ptrdiff_t>(trajectory * samples));
        if (execution) execution->report({"dde", "Integrated delay trajectory",
            static_cast<double>(trajectory + 1) / historyValues.size(),
            trajectory + 1, historyValues.size()});
    }
    return result;
}

}  // namespace dynsyslab::core
