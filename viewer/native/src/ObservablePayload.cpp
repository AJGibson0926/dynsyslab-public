#include "ObservablePayload.hpp"
#include "ArchiveSystem.hpp"
#include "DisplayCoordinates.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/MathText.hpp"

#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace dsl {
namespace {

using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;

const NpyArray* states(const NpzRecord& record) {
    for (const char* name : {"values", "states", "controlled_states",
                             "free_training_values"})
        if (record.contains(name) && record.at(name).shape.size() == 2)
            return &record.at(name);
    return nullptr;
}


std::vector<std::size_t> sampleIndices(std::size_t count) {
    constexpr std::size_t limit = 100000;
    const std::size_t shown = std::min(count, limit);
    std::vector<std::size_t> result(shown);
    for (std::size_t index = 0; index < shown; ++index)
        result[index] = shown < 2 ? 0 : index * (count - 1) / (shown - 1);
    return result;
}

std::vector<double> times(const NpzRecord& record, std::size_t samples,
                          const std::vector<std::size_t>& indices) {
    std::vector<double> source(samples);
    for (std::size_t index = 0; index < samples; ++index) source[index] = static_cast<double>(index);
    for (const char* name : {"time", "controlled_time", "control_time", "free_time"})
        if (record.contains(name) && record.at(name).shape.size() == 1
            && record.at(name).elementCount() >= samples) {
            source = record.at(name).doubles(); break;
        }
    std::vector<double> result; result.reserve(indices.size());
    for (const auto index : indices) result.push_back(source[index]);
    return result;
}

std::string label(const dynsyslab::core::Observable& observable) {
    auto result = observable.descriptor().label;
    if (!observable.descriptor().units.empty() && observable.descriptor().units != "dimensionless") result += " (" + observable.descriptor().units + ')';
    return result;
}

Payload timePayload(const std::filesystem::path& source, const std::string& view,
                    const NpzRecord& record, const json::Value& metadata,
                    const dynsyslab::core::ODESystem& system, bool tolerateInvalid = false) {
    const auto* stateArray = states(record);
    if (!stateArray || !stateArray->shape[0] || stateArray->shape[1] != system.descriptor().stateDimension)
        throw std::runtime_error("observable time view has incompatible state data");
    const auto dimension = stateArray->shape[1];
    const auto indices = sampleIndices(stateArray->shape[0]);
    const auto current = displayCoordinateTransform(metadata, dimension);
    const auto alternate = displayCoordinateTransform(metadata, dimension,
        current.coordinates == "dimensionless" ? "dimensional" : "dimensionless");
    const auto time = transformDisplayTime(times(record, stateArray->shape[0], indices), current);
    const auto catalog = system.observableCatalog();
    Payload payload{source.u8string(), view, metadata.at("label").string(source.stem().u8string()), 2};
    for (const auto& name : catalog.names()) {
        const auto& observable = catalog.at(name);
        const auto& descriptor = observable.descriptor();
        // Catalog evaluators consume the same native state vector as this system's dynamics.
        const auto nativeCoordinates = system.descriptor().preferredCoordinates;
        const auto native = displayCoordinateTransform(metadata, dimension, nativeCoordinates);
        if (native.coordinates != nativeCoordinates)
            throw std::runtime_error("Observable input coordinates are unavailable: " + name);
        const auto values = transformDisplayStates(*stateArray, native).doubles();
        const auto parts = descriptor.complexValued ? 2u : 1u;
        const auto outputs = observable.outputDimension();
        std::vector<std::vector<double>> evaluated(outputs * parts);
        for (auto& output : evaluated) output.reserve(indices.size());
        for (const auto index : indices) {
            std::vector<std::complex<double>> row;
            try { row = observable.evaluate(values.data() + index * dimension); }
            catch (...) {
                if (!tolerateInvalid) throw;
                const auto invalid = std::numeric_limits<double>::quiet_NaN();
                row.assign(outputs, {invalid, invalid});
            }
            for (std::size_t component = 0; component < outputs; ++component) {
                evaluated[component * parts].push_back(row[component].real());
                if (parts == 2) evaluated[component * parts + 1].push_back(row[component].imag());
            }
        }
        for (std::size_t component = 0; component < outputs; ++component) {
            const auto outputName = component < descriptor.outputNames.size()
                ? descriptor.outputNames[component] : std::to_string(component + 1);
            // Deduplicate only factory-proven state projections, never names or sampled curves.
            const auto identity = observable.stateComponentIndex(component);
            if (outputs > 1 && identity && !descriptor.complexValued) {
                bool covered = false;
                for (const auto& scalarName : catalog.names()) {
                    const auto& scalar = catalog.at(scalarName);
                    if (scalar.outputDimension() == 1 && !scalar.descriptor().complexValued
                        && scalar.stateComponentIndex(0) == identity
                        && scalar.descriptor().stateCoordinateSystem == descriptor.stateCoordinateSystem
                        && scalar.descriptor().units == descriptor.units) { covered = true; break; }
                }
                if (covered) continue;
            }
            auto title = scientificMathLabel(label(observable));
            if (outputs > 1) title += " — " + scientificMathLabel(outputName);
            const bool stateOutput = identity.has_value();
            const auto stateIndex = identity.value_or(0);
            if (stateOutput && outputs > 1)
                title = scientificMathLabel(label(observable)) + " — " + current.stateLabels[stateIndex];
            for (unsigned part = 0; part < parts; ++part) {
                const auto plotTitle = parts == 2 ? title + (part ? " — imaginary" : " — real") : title;
                const auto yLabel = stateOutput ? current.stateLabels[stateIndex] : plotTitle;
                Trace trace; trace.label = plotTitle; trace.role = "data"; trace.x = time;
                trace.y = std::move(evaluated[component * parts + part]);
                if (stateOutput) for (auto& value : trace.y)
                    value *= current.stateMultipliers[stateIndex] / native.stateMultipliers[stateIndex];
                Figure figure{"observable_time", plotTitle, current.timeLabel, yLabel, false};
                figure.observableIdentity = dynsyslab::core::equationsFingerprint(system.configuration())
                    + ":" + name + ":" + std::to_string(component) + ":" + std::to_string(part);
                figure.coordinates = current.coordinates;
                figure.dashboardPage = static_cast<int>(payload.figures.size() / 4);
                figure.traces.push_back(std::move(trace));
                if (alternate.coordinates != current.coordinates) {
                    Figure::CoordinateAxes axes;
                    axes.name = alternate.coordinates;
                    axes.xScale = axes.timeScale = alternate.timeMultiplier / current.timeMultiplier;
                    axes.yScale = stateOutput ? alternate.stateMultipliers[stateIndex] / current.stateMultipliers[stateIndex] : 1.;
                    axes.xLabel = alternate.timeLabel;
                    axes.yLabel = stateOutput ? alternate.stateLabels[stateIndex] : yLabel;
                    axes.traceLabels = {stateOutput ? (outputs > 1 ? scientificMathLabel(label(observable)) + " — " : "")
                        + alternate.stateLabels[stateIndex] : plotTitle};
                    figure.alternateCoordinates = std::move(axes);
                }
                payload.figures.push_back(std::move(figure));
            }
        }
    }
    if (payload.figures.empty()) throw std::runtime_error("system has no observables");
    for (std::size_t page = 0; page * 4 < payload.figures.size(); ++page)
        payload.pageLabels.push_back("Observables " + std::to_string(page * 4 + 1) + "–"
            + std::to_string(std::min(payload.figures.size(), (page + 1) * 4)));
    return payload;
}

Payload heatmapPayload(const std::filesystem::path& source,
                       const NpzRecord& record, const json::Value& metadata,
                       const dynsyslab::core::ODESystem& system,
                       const std::optional<PhaseSpaceDomain>& requested) {
    if (system.descriptor().stateDimension != 2)
        throw std::runtime_error("observable heat maps require a two-state dynamical system");
    auto phase = buildNativePayload(source, "phase", record, metadata).figures.front();
    const auto current = displayCoordinateTransform(metadata, 2);
    auto domain = requested.value_or(PhaseSpaceDomain{automaticAxisLimits(phase), phase.coordinates});
    if (!std::isfinite(domain.bounds.xMin) || !std::isfinite(domain.bounds.xMax)
        || !std::isfinite(domain.bounds.yMin) || !std::isfinite(domain.bounds.yMax)
        || domain.bounds.xMax <= domain.bounds.xMin || domain.bounds.yMax <= domain.bounds.yMin)
        throw std::invalid_argument("Phase-space bounds must be finite and increasing");
    const auto display = displayCoordinateTransform(metadata, 2,
        domain.coordinates.empty() ? current.coordinates : domain.coordinates);
    if (!domain.coordinates.empty() && display.coordinates != domain.coordinates)
        throw std::runtime_error("The phase-space coordinate scale is unavailable for these observables");
    const double xScale = display.stateMultipliers[0] / current.stateMultipliers[0];
    const double yScale = display.stateMultipliers[1] / current.stateMultipliers[1];
    constexpr std::size_t count = 101;
    std::vector<double> x(count), y(count), grid;
    for (std::size_t i = 0; i < count; ++i) {
        x[i] = (domain.bounds.xMin + (domain.bounds.xMax-domain.bounds.xMin)*i/(count-1)) / xScale;
        y[i] = (domain.bounds.yMin + (domain.bounds.yMax-domain.bounds.yMin)*i/(count-1)) / yScale;
    }
    for (const auto yi : y) for (const auto xi : x) {
        // Store the grid in archive coordinates, exactly as the time-series evaluator expects.
        grid.push_back(xi / current.stateMultipliers[0]);
        grid.push_back(yi / current.stateMultipliers[1]);
    }
    NpzRecord gridRecord;
    gridRecord.arrays["values"] = NpyArray::fromDoubles({count*count, 2}, grid);
    auto payload = timePayload(source, "observable_heatmap", gridRecord, metadata, system, true);
    for (auto& figure : payload.figures) {
        auto& trace = figure.traces.front();
        trace.z = std::move(trace.y); trace.x = x; trace.y = y;
        trace.zShape = {static_cast<int>(count), static_cast<int>(count)};
        trace.heatmap = true; trace.encoding = "scalar";
        figure.kind = "observable_heatmap"; figure.square = false;
        figure.colorLabel = figure.yLabel;
        figure.xLabel = current.stateLabels[0]; figure.yLabel = current.stateLabels[1];
        figure.references.clear();
        if (figure.alternateCoordinates) {
            auto& axes = *figure.alternateCoordinates;
            axes.colorLabel = axes.yLabel; axes.zScale = axes.yScale;
            const auto alternate = displayCoordinateTransform(metadata, 2, axes.name);
            axes.xScale = alternate.stateMultipliers[0] / current.stateMultipliers[0];
            axes.yScale = alternate.stateMultipliers[1] / current.stateMultipliers[1];
            axes.xLabel = alternate.stateLabels[0]; axes.yLabel = alternate.stateLabels[1];
        }
        changeDisplayCoordinates(figure, display.coordinates);
    }
    return payload;
}

}  // namespace

Payload buildObservablePayload(const std::filesystem::path& source, const std::string& view,
    const NpzRecord& record, const json::Value& metadata,
    const std::optional<PhaseSpaceDomain>& phaseDomain) {
    auto system = dynsyslab::core::createSystem(archiveSystemConfiguration(metadata));
    auto resolvedMetadata = metadata;
    resolvedMetadata.object()["system"] = system->descriptor().identifier;
    if (metadata.at("integration_form").isNull() && metadata.at("coordinate_system").isNull())
        resolvedMetadata.object()["integration_form"] = system->descriptor().preferredCoordinates;
    if (view == "observable_heatmap") return heatmapPayload(source, record, resolvedMetadata, *system, phaseDomain);
    auto result = timePayload(source, view, record, resolvedMetadata, *system);
    const auto* primary = states(record);
    for (const auto& names : {std::vector<const char*>{"control_target_values", "target_values"},
                             std::vector<const char*>{"true_states", "truth_values"}}) {
        const NpyArray* overlay = nullptr;
        for (const auto* name : names) if (record.contains(name)) { overlay = &record.at(name); break; }
        if (!overlay || overlay->shape != primary->shape) continue;
        if (std::string(names.front()) == "control_target_values"
            && metadata.at("control_target_space").string() == "observable") continue;
        NpzRecord related;
        related.arrays["values"] = *overlay;
        for (const auto* name : {"time", "controlled_time", "control_time", "free_time"})
            if (record.contains(name)) related.arrays[name] = record.at(name);
        auto evaluated = timePayload(source, view, related, resolvedMetadata, *system);
        const bool target = std::string(names.front()) == "control_target_values";
        for (std::size_t i = 0; i < result.figures.size(); ++i) {
            auto trace = std::move(evaluated.figures[i].traces.front());
            trace.label = target ? "target" : "true state";
            trace.role = target ? "control_target" : "truth";
            trace.color = target ? "#ef4444" : "#94a3b8"; trace.lineStyle = "dash";
            if (result.figures[i].alternateCoordinates)
                result.figures[i].alternateCoordinates->traceLabels.push_back(trace.label);
            result.figures[i].traces.push_back(std::move(trace));
        }
    }
    // Match stored coordinate identities; a state reference is not an observable objective.
    if (record.contains("control_observable_target_values")
        && record.contains("control_observable_target_names")) {
        const auto& target = record.at("control_observable_target_values");
        std::vector<std::string> names; std::string name;
        std::istringstream encoded(record.at("control_observable_target_names").string());
        while (std::getline(encoded,name)) names.push_back(name);
        if (target.shape.size() != 2 || target.shape[0] != primary->shape[0]
            || target.shape[1] != names.size())
            throw std::runtime_error("Observable target dimensions do not match coordinate identities");
        const auto samples = sampleIndices(primary->shape[0]);
        const auto values = target.doubles();
        const auto fingerprint = dynsyslab::core::equationsFingerprint(system->configuration());
        for (std::size_t column=0; column<names.size(); ++column)
            for (auto& figure : result.figures) {
                if (figure.observableIdentity != fingerprint+":"+names[column]+":0:0") continue;
                Trace trace; trace.label="target"; trace.role="control_target";
                trace.color="#ef4444"; trace.lineStyle="dash";
                trace.x=figure.traces.front().x;
                for (const auto index:samples) trace.y.push_back(values[index*names.size()+column]);
                figure.traces.push_back(std::move(trace));
                if(figure.alternateCoordinates)
                    figure.alternateCoordinates->traceLabels.push_back("target");
            }
    }
    if (view == "observable_workspace" && system->descriptor().stateDimension == 2) {
        auto maps = heatmapPayload(source, record, resolvedMetadata, *system, phaseDomain);
        if (!maps.figures.empty() && maps.figures.front().coordinates != result.figures.front().coordinates)
            for (auto& figure : result.figures) changeDisplayCoordinates(figure, maps.figures.front().coordinates);
        std::vector<Figure> paired;
        for (std::size_t i = 0; i < result.figures.size(); ++i) {
            paired.push_back(std::move(result.figures[i]));
            paired.push_back(std::move(maps.figures.at(i)));
        }
        result.figures = std::move(paired); result.pageLabels.clear();
        for (std::size_t i = 0; i < result.figures.size(); ++i) {
            result.figures[i].dashboardPage = static_cast<int>(i / 4);
            if (i % 4 == 0) result.pageLabels.push_back("Observables " + std::to_string(i/2+1)
                + "–" + std::to_string(std::min(result.figures.size(), i+4)/2));
        }
    }
    if (view == "observable_workspace") {
        result.pageLabels.clear();
        for (auto& figure : result.figures) figure.dashboardPage = 0;
    }
    return result;
}

}  // namespace dsl
