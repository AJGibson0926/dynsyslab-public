#include "ArchiveSystem.hpp"
#include "DisplayCoordinates.hpp"
#include "HierarchyPayload.hpp"
#include "PhaseContext.hpp"
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/core/Archive.hpp>
#include <limits>
#include <set>
#include <stdexcept>
namespace dsl {
namespace {
json::Value numbers(const std::vector<double> &values) {
    json::Value::Array a;
    for (double x : values)
        a.emplace_back(x);
    return a;
}
bool hierarchy(const CatalogItem &item) {
    return item.kind == "batch" || item.kind == "cluster" || item.kind == "system_run" ||
           item.kind == "forcing_campaign";
}
} // namespace
json::Value Backend::queryPhaseSpace(const CatalogItem &item, const std::string &mode,
                                     std::optional<PhaseSpaceDomain> domain,
                                     const json::Value &options,
                                     const std::atomic_bool *cancelled) const {
    if (mode != "lookup" && mode != "compute")
        throw std::invalid_argument("Phase-space mode must be lookup or compute");
    auto paths = hierarchy(item)
                     ? hierarchyTrajectorySources(item, config_.projectRoot, config_.archiveRoot)
                     : std::vector<std::filesystem::path>{resolveSource(item)};
    json::Value::Array results;
    std::set<std::string> identities;
    for (const auto &source : paths) {
        if (cancelled && cancelled->load())
            throw std::runtime_error("Phase-space request cancelled");
        if (source.extension() != L".npz")
            continue;
        const auto index = dynsyslab::core::readNpzIndex(source);
        auto metadata = json::parse(index.metadataJson);
        auto system = dynsyslab::core::createSystem(archiveSystemConfiguration(metadata));
        const auto &resolved = system->configuration();
        const auto n = system->descriptor().stateDimension;
        const auto fingerprint = metadata.at("equations_fingerprint").string();
        if (!fingerprint.empty() && fingerprint != dynsyslab::core::equationsFingerprint(resolved))
            throw std::runtime_error("Archived equations differ from the installed system; phase "
                                     "context was not substituted");
        json::Value::Object parameters;
        for (const auto &p : resolved.parameters)
            parameters[p.first] = p.second;
        json::Value spec = json::Value::Object{{"identifier", resolved.identifier},
                                               {"preset", resolved.preset},
                                               {"parameters", parameters},
                                               {"driver", json::Value::Object{{"amplitude", resolved.driver.amplitude},
                                                 {"frequency_hz", resolved.driver.frequencyHz}, {"phase_radians", resolved.driver.phaseRadians}}}};
        auto queryOptions = options.isNull() ? json::Value(json::Value::Object{}) : options;
        if (queryOptions.at("lower").isNull() || queryOptions.at("upper").isNull()) {
            if (n < 2)
                throw std::runtime_error("Phase context requires at least two state components");
            const auto native =
                displayCoordinateTransform(metadata, n, system->descriptor().preferredCoordinates);
            if (native.coordinates != system->descriptor().preferredCoordinates)
                throw std::runtime_error("Native phase coordinates are unavailable");
            std::vector<double> low(n), high(n);
            std::optional<dynsyslab::core::NpzRecord> record;
            if (n > 2 || !domain) {
                record = dynsyslab::core::readNpz(source);
                const dynsyslab::core::NpyArray *state = nullptr;
                for (const char *name :
                     {"values", "states", "controlled_states", "free_training_values"})
                    if (record->contains(name) && record->at(name).shape.size() == 2 &&
                        record->at(name).shape[1] == n) {
                        state = &record->at(name);
                        break;
                    }
                if (!state)
                    throw std::runtime_error("Phase context requires archived state trajectories");
                const auto values = transformDisplayStates(*state, native).doubles();
                std::fill(low.begin(), low.end(), std::numeric_limits<double>::infinity());
                std::fill(high.begin(), high.end(), -std::numeric_limits<double>::infinity());
                for (size_t i = 0; i < values.size(); ++i)
                    if (std::isfinite(values[i])) {
                        low[i % n] = std::min(low[i % n], values[i]);
                        high[i % n] = std::max(high[i % n], values[i]);
                    }
                for (size_t j = 0; j < n; ++j) {
                    if (high[j] == low[j]) {
                        low[j] -= .5;
                        high[j] += .5;
                    }
                    const double pad = .045 * (high[j] - low[j]);
                    low[j] -= pad;
                    high[j] += pad;
                }
            }
            auto displayed = domain;
            if (!displayed) {
                const auto phase =
                    buildNativePayload(source, "phase", *record, metadata).figures.front();
                displayed = PhaseSpaceDomain{automaticAxisLimits(phase), phase.coordinates};
            }
            const auto display = displayCoordinateTransform(metadata, n, displayed->coordinates);
            const double sx = display.stateMultipliers[0] / native.stateMultipliers[0],
                         sy = display.stateMultipliers[1] / native.stateMultipliers[1];
            low[0] = displayed->bounds.xMin / sx;
            high[0] = displayed->bounds.xMax / sx;
            low[1] = displayed->bounds.yMin / sy;
            high[1] = displayed->bounds.yMax / sy;
            queryOptions.object()["lower"] = numbers(low);
            queryOptions.object()["upper"] = numbers(high);
        }
        const auto &controller = metadata.at("controller_configuration");
        const auto &settings = controller.at("settings");
        const bool modelAmplitude = !settings.at("model_amplitudes").isNull() || !settings.at("model_amplitude").isNull();
        const auto &amplitudes = settings.at(modelAmplitude ? "model_amplitudes" : "amplitudes").array();
        const bool nonzeroAmplitude = amplitudes.empty()
            ? settings.at(modelAmplitude ? "model_amplitude" : "amplitude").number(1.) != 0
            : std::any_of(amplitudes.begin(), amplitudes.end(), [](const auto &v) { return v.number() != 0; });
        const bool waveform = controller.at("kind").string() == "waveform" &&
            (nonzeroAmplitude || settings.at("offset").number() != 0 || settings.at("model_offset").number() != 0);
        bool internalClock = false;
        if (resolved.identifier == "reddy_szeri") {
            const auto &g = resolved.dimensionlessGroups;
            internalClock = resolved.parameters.at("translation_enabled") != 0 && g.at("angular_frequency") != 0 &&
                (g.at("pressure_amplitude") != 0 || g.at("pressure_derivative_amplitude") != 0);
        }
        if (queryOptions.at("dynamics").isNull() && (waveform || internalClock || resolved.driver.amplitude != 0))
            queryOptions.object()["dynamics"] = "periodic-forcing";
        // A higher-dimensional projection has no unique vector field without a specified slice.
        if (n > 2 && queryOptions.at("slice_state").isNull() &&
            queryOptions.at("vector_field").isNull())
            queryOptions.object()["vector_field"] = false;
        const auto identity = json::dump(spec) + json::dump(queryOptions) + json::dump(controller) + metadata.at("control_evaluation").string();
        if (!identities.insert(identity).second)
            continue;
        json::Value request = json::Value::Object{{"schema", "dynsyslab.phase-space-query/v1"},
                                                  {"mode", mode},
                                                  {"system", spec},
                                                  {"options", queryOptions}};
        if (waveform && queryOptions.at("dynamics").string() == "periodic-forcing") {
            request.object()["controller"] = controller;
            request.object()["control_evaluation"] = metadata.at("control_evaluation");
        }
        const auto result = callPhaseSidecar(request, config_.archiveRoot, cancelled);
        results.emplace_back(json::Value::Object{{"source", source.u8string()},
                                                 {"display_metadata", metadata},
                                                 {"request", request},
                                                 {"result", result}});
    }
    return json::Value::Object{
        {"schema", "dynsyslab.viewer.phase-space/v1"}, {"mode", mode}, {"results", results}};
}
Payload Backend::applyPhaseSpace(Payload payload, const std::vector<CatalogItem> &items,
                                 const std::string &mode, std::optional<PhaseSpaceDomain> domain,
                                 const std::atomic_bool *cancelled) const {
    if (std::none_of(payload.figures.begin(), payload.figures.end(),
                     [](const auto &f) { return f.kind == "phase"; }))
        return payload;
    // Remove an earlier context before changing its domain; it must not influence new plot bounds.
    for (auto &f : payload.figures)
        if (f.kind == "phase") {
            for (size_t i = f.traces.size(); i-- > 0;)
                if (f.traces[i].role.rfind("phase_", 0) == 0) {
                    f.traces.erase(f.traces.begin() + i);
                    if (f.alternateCoordinates && i < f.alternateCoordinates->traceLabels.size())
                        f.alternateCoordinates->traceLabels.erase(
                            f.alternateCoordinates->traceLabels.begin() + i);
                }
        }
    if (!domain) {
        for (const auto &f : payload.figures)
            if (f.kind == "phase") {
                domain = PhaseSpaceDomain{automaticAxisLimits(f), f.coordinates};
                break;
            }
    }
    std::set<std::string> seen;
    for (const auto &item : items)
        try {
            const auto bundle = queryPhaseSpace(item, mode, domain, {}, cancelled);
            for (const auto &entry : bundle.at("results").array()) {
                const auto &result = entry.at("result");
                const auto key = result.at("cache").at("key").string();
                if (!seen.insert(key).second)
                    continue;
                if (result.at("cache").at("status").string() == "miss")
                    continue;
                for (auto &figure : payload.figures)
                    appendPhaseContext(figure, entry.at("display_metadata"), result);
                for (const auto &diagnostic : result.at("diagnostics").array()) {
                    auto text = "Phase context: " + diagnostic.string();
                    if (std::find(payload.warnings.begin(), payload.warnings.end(), text) ==
                        payload.warnings.end())
                        payload.warnings.push_back(text);
                }
            }
        } catch (const std::exception &error) {
            if (mode == "compute")
                throw;
            payload.warnings.push_back(std::string("Phase context unavailable: ") + error.what());
        }
    return payload;
}
} // namespace dsl
