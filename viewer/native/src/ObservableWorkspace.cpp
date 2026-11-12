#include "ObservablePayload.hpp"
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "HierarchyPayload.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace dsl {
void appendObservablePayload(Payload& target, Payload source, const std::string& label,
                             const std::string& color, bool dashed) {
    for (auto& figure : source.figures) {
        auto match = std::find_if(target.figures.begin(), target.figures.end(), [&](const Figure& other) {
            if (figure.observableIdentity.empty() || other.observableIdentity != figure.observableIdentity) return false;
            if (other.kind != figure.kind || other.title != figure.title || other.xLabel != figure.xLabel || other.yLabel != figure.yLabel
                || other.coordinates != figure.coordinates || bool(other.alternateCoordinates) != bool(figure.alternateCoordinates)) return false;
            if (figure.kind == "observable_heatmap") {
                const auto& a = other.traces.front(); const auto& b = figure.traces.front();
                if (a.x != b.x || a.y != b.y || a.z.size() != b.z.size()) return false;
                // Same defined observable/system/index on the same grid; sampled values are not an identity test.
            }
            if (!figure.alternateCoordinates) return true;
            const auto& a = *other.alternateCoordinates; const auto& b = *figure.alternateCoordinates;
            return a.name == b.name && a.xScale == b.xScale && a.yScale == b.yScale && a.zScale == b.zScale;
        });
        for (auto& trace : figure.traces) {
            trace.label = trace.role == "data" ? label : label + " — " + trace.label;
            if (!color.empty() && trace.role == "data") trace.color = color;
            if (dashed && trace.role == "data") trace.lineStyle = "dash";
        }
        if (figure.alternateCoordinates) {
            figure.alternateCoordinates->traceLabels.clear();
            for (const auto& trace : figure.traces) figure.alternateCoordinates->traceLabels.push_back(trace.label);
        }
        if (match == target.figures.end()) {
            auto position = target.figures.end();
            if (figure.kind == "observable_heatmap")
                for (auto i = target.figures.begin(); i != target.figures.end(); ++i)
                    if (i->title == figure.title) position = std::next(i);
            target.figures.insert(position, std::move(figure));
        } else if (figure.kind == "observable_heatmap") {
            for (const auto& ref : figure.references)
                if (std::none_of(match->references.begin(),match->references.end(),[&](const auto& other) {
                    return other.x == ref.x && other.y == ref.y && other.equilibrium == ref.equilibrium;
                })) match->references.push_back(ref);
        } else {
            match->traces.insert(match->traces.end(), figure.traces.begin(), figure.traces.end());
            match->references.insert(match->references.end(), figure.references.begin(), figure.references.end());
            if (match->alternateCoordinates)
                match->alternateCoordinates->traceLabels.insert(match->alternateCoordinates->traceLabels.end(),
                    figure.alternateCoordinates->traceLabels.begin(), figure.alternateCoordinates->traceLabels.end());
        }
    }
    target.warnings.insert(target.warnings.end(), source.warnings.begin(), source.warnings.end());
    target.pageLabels.clear();
    for (std::size_t i = 0; i < target.figures.size(); ++i) {
        target.figures[i].dashboardPage = 0;
    }
}

bool supportsObservableWorkspace(const CatalogItem& item) {
    return item.kind == "trajectory" || item.kind == "controlled_trajectory"
        || item.kind == "batch" || item.kind == "system_run" || item.kind == "cluster"
        || item.kind == "forcing_campaign" || std::find(item.views.begin(), item.views.end(),
            "observable_time") != item.views.end();
}

Payload Backend::loadObservables(const std::vector<CatalogItem>& items,
    std::optional<PhaseSpaceDomain> phaseDomain) const {
    if (items.empty()) throw std::runtime_error("Select trajectories on Main before adding observables");
    // A shared phase domain also gives compared/grouped trajectories the same map grid.
    const auto hierarchy = [](const CatalogItem& item) {
        return item.kind == "batch" || item.kind == "system_run" || item.kind == "cluster" || item.kind == "forcing_campaign";
    };
    if (!phaseDomain) {
        Figure combined;
        for (const auto& item : items) {
            try {
                const auto phase = loadPayload(item, hierarchy(item) ? "ensemble_phase" : "phase").figures.front();
                if (combined.traces.empty()) combined = phase;
                else if (phase.coordinates == combined.coordinates) {
                    combined.traces.insert(combined.traces.end(), phase.traces.begin(), phase.traces.end());
                    combined.references.insert(combined.references.end(), phase.references.begin(), phase.references.end());
                }
            } catch (const std::exception&) { /* Higher-dimensional systems can still supply time series. */ }
        }
        if (!combined.traces.empty()) phaseDomain = PhaseSpaceDomain{automaticAxisLimits(combined), combined.coordinates};
    }
    Payload result{"", "observable_workspace", "Observables", 2};
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (!supportsObservableWorkspace(items[i]))
            throw std::runtime_error("This Main object has no system observable trajectories: " + items[i].label);
        Payload payload;
        if (hierarchy(items[i])) payload = buildHierarchyPayload(items[i], "observable_workspace",
            config_.projectRoot, config_.archiveRoot, phaseDomain);
        else {
            const auto source = resolveSource(items[i]);
            const auto record = dynsyslab::core::readNpz(source);
            const auto metadata = json::parse(record.metadataJson());
            payload = buildObservablePayload(source, "observable_workspace", record, metadata, phaseDomain);
        }
        if (!result.source.empty()) result.source += '\n';
        result.source += payload.source;
        // Hierarchy payloads already label and color their individual trajectories.
        if (items.size() == 1) return payload;
        appendObservablePayload(result, payload, items[i].label, {}, i != 0);
    }
    return result;
}
}
