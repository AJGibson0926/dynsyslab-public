#include "FourierPayload.hpp"
#include "DisplayCoordinates.hpp"
#include "TrajectoryPresentation.hpp"

#include <dynsyslab/core/Extrema.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dsl {
namespace {

using dynsyslab::core::FourierOptions;
using dynsyslab::core::FourierTransform;
using dynsyslab::core::FourierWindow;
using dynsyslab::core::TimeSeries;
using dynsyslab::core::TimeSeriesArchiveSelection;

struct DisplaySeries {
    TimeSeries values;
    DisplayCoordinateTransform coordinates;
};

DisplaySeries loadDisplaySeries(const dynsyslab::core::NpzRecord& archive,
                                const json::Value& metadata, bool control) {
    // Prefer the input applied at every controller step, on the controller clock.
    const auto rows = [&](const char* name) {
        return archive.contains(name) && !archive.at(name).shape.empty() ? archive.at(name).shape[0] : 0;
    };
    const bool controllerClock = control && rows("applied_input") > 1
        && rows("applied_input") == rows("estimator_time");
    const TimeSeriesArchiveSelection selection = controllerClock
        ? TimeSeriesArchiveSelection{{"estimator_time"}, {"applied_input"}}
        : control
        ? TimeSeriesArchiveSelection{
            {"control_time", "time", "controlled_time", "free_time"},
            {"control_values", "controls", "controlled_controls",
             "controlled_controls_sample_aligned"}}
        : TimeSeriesArchiveSelection{
            {"time", "controlled_time", "control_time", "free_time"},
            {"values", "states", "controlled_states", "free_training_values"}};
    auto source = dynsyslab::core::loadTimeSeries(archive, selection);
    auto coordinates = displayCoordinateTransform(metadata, source.componentCount());
    if (control) {
        if (coordinates.timeMultiplier == 1.0)
            return {std::move(source), std::move(coordinates)};
        std::vector<double> time(source.sampleCount());
        std::vector<double> values(source.valuesData(),
            source.valuesData() + source.sampleCount() * source.componentCount());
        for (std::size_t sample = 0; sample < source.sampleCount(); ++sample)
            time[sample] = source.time(sample) * coordinates.timeMultiplier;
        return {TimeSeries(std::move(time), std::move(values), source.componentCount()),
                std::move(coordinates)};
    }
    const bool identityStates = std::all_of(
        coordinates.stateMultipliers.begin(), coordinates.stateMultipliers.end(),
        [](double value) { return value == 1.0; });
    if (coordinates.timeMultiplier == 1.0 && identityStates)
        return {std::move(source), std::move(coordinates)};
    std::vector<double> time(source.sampleCount());
    std::vector<double> values(source.sampleCount() * source.componentCount());
    for (std::size_t sample = 0; sample < source.sampleCount(); ++sample) {
        time[sample] = source.time(sample) * coordinates.timeMultiplier;
        for (std::size_t component = 0; component < source.componentCount(); ++component) {
            const double multiplier = component < coordinates.stateMultipliers.size()
                ? coordinates.stateMultipliers[component] : 1.0;
            values[sample * source.componentCount() + component]
                = source.value(sample, component) * multiplier;
        }
    }
    return {TimeSeries(std::move(time), std::move(values), source.componentCount()),
            std::move(coordinates)};
}

std::vector<std::size_t> displayBins(const FourierTransform& transform,
                                     std::size_t component) {
    // Peak-preserving thinning only beyond 65536 bins, so an --x-range zoom keeps
    // full resolution (a 0.001 tau controller clock gives about 30000 bins).
    constexpr std::size_t limit = 65536;
    if (transform.binCount <= limit) {
        std::vector<std::size_t> result(transform.binCount);
        for (std::size_t index = 0; index < result.size(); ++index) result[index] = index;
        return result;
    }
    const std::size_t block = (transform.binCount - 1 + limit - 2) / (limit - 1);
    std::vector<std::size_t> result{0};
    for (std::size_t first = 1; first < transform.binCount; first += block) {
        const std::size_t last = std::min(transform.binCount, first + block);
        std::size_t peak = first;
        for (std::size_t bin = first + 1; bin < last; ++bin)
            if (transform.amplitude(bin, component) > transform.amplitude(peak, component))
                peak = bin;
        result.push_back(peak);
    }
    if (result.back() != transform.binCount - 1) result.push_back(transform.binCount - 1);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

Trace amplitudeTrace(const FourierTransform& transform, std::size_t component) {
    Trace result; result.label = "Fourier amplitude"; result.role = "spectrum";
    const auto bins = displayBins(transform, component);
    result.x.reserve(bins.size()); result.y.reserve(bins.size());
    for (const auto bin : bins) {
        result.x.push_back(transform.frequencies[bin]);
        result.y.push_back(transform.amplitude(bin, component));
    }
    return result;
}

}  // namespace

Payload buildFourierPayload(const std::filesystem::path& source,
                            const dynsyslab::core::NpzRecord& archive,
                            const json::Value& metadata) {
    auto displayed = loadDisplaySeries(archive, metadata, false);
    FourierOptions options;
    options.oneSided = true;
    options.removeMean = true;
    options.window = FourierWindow::hann;
    const auto transform = displayed.values.fourier().transform(options);
    Payload payload{source.u8string(), "fourier",
        metadata.at("label").string(source.stem().u8string()), 2};
    if (!displayed.coordinates.warning.empty())
        payload.warnings.push_back(displayed.coordinates.warning);
    const std::string frequencyLabel = displayed.coordinates.coordinates == "dimensionless"
        ? "frequency (cycles / tau)"
        : metadata.at("time_scale_to_seconds").number(0.0) > 0 ? "frequency (Hz)" : "frequency (cycles / time)";
    for (std::size_t component = 0;
         component < displayed.values.componentCount(); ++component) {
        const std::string name = component < displayed.coordinates.stateLabels.size()
            ? displayed.coordinates.stateLabels[component]
            : "state " + std::to_string(component);
        Figure figure{"fourier", name + " Fourier spectrum",
                      frequencyLabel, "amplitude", false};
        figure.valueScale = "linear";  // log on request (--scale log or the Log toggle)
        figure.traces.push_back(amplitudeTrace(*transform, component));
        // Open on the dominant band (99.5% of the amplitude, at least 1.25x the peak).
        if (const auto band = dynsyslab::core::dominantFrequencyRange(figure.traces.back().x, figure.traces.back().y);
            band.size() == 2 && band[1] > band[0]) figure.focusX = std::make_pair(band[0], band[1]);
        payload.figures.push_back(std::move(figure));
    }
    decorateTrajectory(payload, metadata, displayed.coordinates, {});
    return payload;
}

Payload buildControlFourierPayload(
    const std::filesystem::path& source,
    const dynsyslab::core::NpzRecord& archive,
    const json::Value& metadata) {
    auto displayed = loadDisplaySeries(archive, metadata, true);
    FourierOptions options;
    options.oneSided = true;
    options.removeMean = true;
    options.window = FourierWindow::hann;
    const auto transform = displayed.values.fourier().transform(options);
    Payload payload{source.u8string(), "control_fourier",
        metadata.at("label").string(source.stem().u8string()), 2};
    if (!displayed.coordinates.warning.empty())
        payload.warnings.push_back(displayed.coordinates.warning);
    const std::string controlCoordinates = metadata.at("control_coordinate_system")
        .string(metadata.at("integration_form").string(
            metadata.at("coordinate_system").string()));
    // Some retained control archives predate integration_form but do carry the
    // authoritative control-coordinate tag.  The controller and dynamics share
    // the same time scaling, so that tag is sufficient to label its spectrum.
    const bool dimensionlessTime = displayed.coordinates.coordinates == "dimensionless" ||
        controlCoordinates == "dimensionless";
    const std::string frequencyLabel = dimensionlessTime
        ? "frequency (cycles / tau)"
        : metadata.at("time_scale_to_seconds").number(0.0) > 0 ? "frequency (Hz)" : "frequency (cycles / time)";
    const std::string amplitudeLabel = controlCoordinates == "dimensionless"
        ? "control amplitude (dimensionless)" : "control amplitude";
    for (std::size_t component = 0;
         component < displayed.values.componentCount(); ++component) {
        const std::string name = "u" + std::to_string(component);
        Figure figure{"control_fourier", name + " Fourier spectrum",
                      frequencyLabel, amplitudeLabel, false};
        figure.valueScale = "linear";
        auto trace = amplitudeTrace(*transform, component);
        trace.label = name + " amplitude";
        trace.role = "control_spectrum";
        figure.traces.push_back(std::move(trace));
        payload.figures.push_back(std::move(figure));
    }
    decorateTrajectory(payload, metadata, displayed.coordinates, {});
    return payload;
}

}  // namespace dsl
