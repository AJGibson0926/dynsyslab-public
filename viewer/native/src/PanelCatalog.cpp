#include "dynsyslab_viewer/PanelCatalog.hpp"

#include <algorithm>
#include <set>

namespace dsl {
namespace {

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

std::string family(const std::string& route, const Figure& figure) {
    if (route == "control_fourier" || figure.kind == "control_fourier")
        return "Control spectra";
    if (figure.title == "Control history") return "Control";
    if (route == "fourier" || figure.kind == "fourier") return "Fourier spectra";
    if (figure.kind == "time") return "Time series";
    if (figure.kind == "phase" || figure.kind == "three_d") return "State space";
    if (route.find("observable") != std::string::npos) return "Observables";
    if (figure.kind.find("koopman") != std::string::npos ||
        route.find("krom") != std::string::npos) return "Koopman";
    if (figure.kind.find("error") != std::string::npos) return "Error maps";
    return "Other scientific plots";
}

std::string fingerprint(const Figure& figure) {
    std::string result = figure.kind + '\n' + figure.title + '\n' +
        figure.xLabel + '\n' + figure.yLabel;
    for (const auto& trace : figure.traces)
        result += '\n' + trace.label + ':' + trace.role;
    return result;
}

}  // namespace

std::vector<std::string> panelSourceViews(const CatalogItem& item) {
    static const std::vector<std::string> preferred{
        "time", "phase", "three_d", "fourier", "control_fourier",
        "observable_time", "observable_heatmap", "field", "lattice",
        "ensemble_control", "ensemble_time", "ensemble_phase",
        "campaign_summary", "campaign_bifurcation",
        "koopman_spectrum", "krom_dashboard", "control_map_error",
        "observability", "koopman_validation", "eigenvalue_optimization",
        "krom_control_campaign", "data_overview"};
    std::vector<std::string> result;
    for (const auto& view : preferred)
        if (contains(item.views, view)) result.push_back(view);
    for (const auto& view : item.views) {
        if (view == "open" || view == "control" || view == "time_frequency")
            continue;
        if (!contains(result, view)) result.push_back(view);
    }
    // Fourier views are derived lazily from the decoded time/control arrays.
    // They must not depend on a catalogue entry written before those Viewer
    // routes existed.  Probe both for any trajectory-like archive; the native
    // payload builder rejects an unavailable control history without affecting
    // the remaining choices.
    const bool trajectoryLike = contains(item.views, "data_overview") ||
        contains(item.views, "time") || contains(item.views, "control") ||
        contains(item.views, "time_frequency");
    // Historical catalogs may advertise only the dashboard. Probe its component
    // routes through the same payload builder for UI and background inspection.
    if (trajectoryLike) for (const auto* view : {"time", "phase", "three_d", "observable_time"})
        if (!contains(result,view)) result.push_back(view);
    if (trajectoryLike && !contains(result, "fourier"))
        result.push_back("fourier");
    if (trajectoryLike && !contains(result, "control_fourier"))
        result.push_back("control_fourier");
    if (!contains(item.views, "time") && contains(item.views, "control"))
        result.insert(result.begin(), "control");
    return result;
}

std::vector<PanelChoice> composePanelChoices(
    std::vector<std::pair<std::string, Payload>> payloads) {
    std::vector<PanelChoice> result;
    std::set<std::string> seen;
    for (auto& entry : payloads) {
        for (auto& figure : entry.second.figures) {
            const std::string key = fingerprint(figure);
            if (!seen.insert(key).second) continue;
            PanelChoice choice;
            choice.route = entry.first;
            choice.family = family(entry.first, figure);
            choice.label = figure.title.empty() ? figure.yLabel : figure.title;
            choice.detail = figure.yLabel + "  vs  " + figure.xLabel;
            choice.figure = std::move(figure);
            choice.figure.dashboardPage = -1;
            choice.figure.gridSpan = 1;
            choice.figure.rowSpan = 1;
            result.push_back(std::move(choice));
        }
    }
    std::stable_sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        return left.family == right.family ? left.label < right.label
                                           : left.family < right.family;
    });
    return result;
}

}  // namespace dsl
