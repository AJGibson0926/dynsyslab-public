#include "dynsyslab_viewer/ComparisonPayload.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {

constexpr std::size_t kFiguresPerPage = 4;
constexpr const char* kSharedTargetColor = "#202020";
const std::string kSeparator = " — ";  // run label separator (UTF-8: five bytes)

using Pages = std::vector<std::vector<const Figure*>>;

Pages pages(const Payload& payload) {
    int greatest = -1;
    for (const auto& figure : payload.figures)
        greatest = std::max(greatest, figure.dashboardPage);
    if (greatest >= 0) {
        Pages result(static_cast<std::size_t>(greatest + 1));
        for (const auto& figure : payload.figures)
            if (figure.dashboardPage >= 0)
                result[static_cast<std::size_t>(figure.dashboardPage)].push_back(&figure);
        return result;
    }
    Pages result((payload.figures.size() + kFiguresPerPage - 1) / kFiguresPerPage);
    for (std::size_t index = 0; index < payload.figures.size(); ++index)
        result[index / kFiguresPerPage].push_back(&payload.figures[index]);
    return result;
}

bool overlayable(const Figure& first, const Figure& second) {
    if (first.kind != second.kind || first.xLabel != second.xLabel
        || first.yLabel != second.yLabel) return false;
    if (first.coordinates != second.coordinates ||
        bool(first.alternateCoordinates) != bool(second.alternateCoordinates)) return false;
    if (first.alternateCoordinates) {
        const auto& a = *first.alternateCoordinates;
        const auto& b = *second.alternateCoordinates;
        if (a.name != b.name || a.xScale != b.xScale || a.yScale != b.yScale
            || a.zScale != b.zScale || a.timeScale != b.timeScale) return false;
    }
    return std::none_of(first.traces.begin(), first.traces.end(),
               [](const Trace& trace) { return trace.heatmap; })
        && std::none_of(second.traces.begin(), second.traces.end(),
               [](const Trace& trace) { return trace.heatmap; });
}

bool sharedRole(const Trace& trace) { return trace.role == "control_target" || trace.role == "truth"; }

bool sameCurve(const Trace& a, const Trace& b) {
    if (a.x.size() != b.x.size() || a.y.size() != b.y.size()) return false;
    const auto close = [](const std::vector<double>& p, const std::vector<double>& q) {
        for (std::size_t i = 0; i < p.size(); ++i)
            if (std::abs(p[i] - q[i]) > 1e-12 * (1.0 + std::abs(p[i]))) return false;
        return true;
    };
    return close(a.x, b.x) && close(a.y, b.y);
}

void appendTraces(Figure& target, const Figure& source, const ComparisonEntry& entry,
                  std::size_t index, bool colored) {
    if (target.alternateCoordinates && source.alternateCoordinates)
        for (const auto& traceLabel : source.alternateCoordinates->traceLabels)
            target.alternateCoordinates->traceLabels.push_back(entry.label + " — " + traceLabel);
    for (auto trace : source.traces) {
        trace.label = entry.label + " — " + trace.label;
        const bool own = !sharedRole(trace) && !trace.heatmap
            && trace.encoding != "complex-domain-color";
        if (colored && own) {
            trace.color = entry.color;
            if (trace.encoding == "time") { trace.encoding.clear(); trace.colorValues.clear(); }
            if (trace.role != "estimate") trace.lineStyle = "solid";
        } else if (!colored && index > 0 && trace.lineStyle == "solid" && !trace.marker
                   && trace.encoding != "complex-domain-color") {
            trace.lineStyle = "dash";
        }
        target.traces.push_back(std::move(trace));
    }
}

// Identical targets (or true states) are drawn once: unprefixed when every run that has
// such a trace shares it (runs without one, e.g. an open-loop reference, do not count),
// otherwise prefixed with the runs that share them. In colored comparisons targets are
// dark grey, the most widely shared one solid and the others dashed.
void mergeSharedTraces(Figure& figure, std::size_t, bool colored) {
    std::vector<bool> drop(figure.traces.size(), false);
    std::map<std::string, std::set<std::string>> holders;  // role -> runs with such a trace
    for (const auto& trace : figure.traces) if (sharedRole(trace)) {
        const auto dash = trace.label.find(kSeparator);
        holders[trace.role].insert(dash == std::string::npos ? std::string{} : trace.label.substr(0, dash));
    }
    std::vector<std::pair<std::size_t, std::size_t>> groups;  // kept trace, copies
    for (std::size_t i = 0; i < figure.traces.size(); ++i) {
        if (drop[i] || !sharedRole(figure.traces[i])) continue;
        std::vector<std::size_t> copies{i};
        for (std::size_t j = i + 1; j < figure.traces.size(); ++j)
            if (!drop[j] && figure.traces[j].role == figure.traces[i].role
                && sameCurve(figure.traces[i], figure.traces[j])) copies.push_back(j);
        const auto split = [](const std::string& label) {
            const auto dash = label.find(kSeparator);
            return dash == std::string::npos ? std::make_pair(std::string{}, label)
                                             : std::make_pair(label.substr(0, dash), label.substr(dash + kSeparator.size()));
        };
        auto& kept = figure.traces[i];
        const auto name = split(kept.label).second;
        const auto runs = holders[kept.role].size();
        if (copies.size() == runs && runs > 1) kept.label = name;
        else if (copies.size() > 1) {
            std::string owners;
            for (const auto copy : copies) owners += (owners.empty() ? "" : ", ") + split(figure.traces[copy].label).first;
            kept.label = owners + " — " + name;
        }
        for (std::size_t k = 1; k < copies.size(); ++k) drop[copies[k]] = true;
        groups.emplace_back(i, copies.size());
    }
    if (colored) {
        std::size_t widest = 0;
        for (const auto& group : groups) widest = std::max(widest, group.second);
        bool solidTaken = false;
        for (const auto& [index, copies] : groups) {
            auto& trace = figure.traces[index];
            if (trace.role != "control_target") continue;
            trace.color = kSharedTargetColor; trace.lineWidthPoints = 2.4;
            const bool solid = copies == widest && !solidTaken;
            trace.lineStyle = solid ? "solid" : "dash"; solidTaken = solidTaken || solid;
        }
    }
    std::vector<Trace> kept;
    for (std::size_t i = 0; i < figure.traces.size(); ++i)
        if (!drop[i]) kept.push_back(std::move(figure.traces[i]));
    // References are drawn first so the runs stay on top of them.
    std::stable_partition(kept.begin(), kept.end(), [](const Trace& t) { return sharedRole(t); });
    figure.traces = std::move(kept);
}

std::string joinedLabels(const std::vector<ComparisonEntry>& entries) {
    std::string result;
    for (std::size_t i = 0; i < entries.size(); ++i) result += (i ? " vs " : "") + entries[i].label;
    return result;
}

Figure overlay(const std::vector<const Figure*>& figures, const std::vector<ComparisonEntry>& entries,
               int page, bool colored) {
    Figure result = *figures.front();
    result.title = joinedLabels(entries) + " — " + figures.front()->title;
    result.titleColor.clear(); result.dashboardPage = page;
    result.traces.clear(); result.inset.reset();
    if (result.alternateCoordinates) result.alternateCoordinates->traceLabels.clear();
    for (std::size_t k = 1; k < figures.size(); ++k)
        result.references.insert(result.references.end(), figures[k]->references.begin(),
                                 figures[k]->references.end());
    for (std::size_t k = 0; k < figures.size(); ++k)
        appendTraces(result, *figures[k], entries[k], k, colored);
    mergeSharedTraces(result, figures.size(), colored);
    if (colored && result.colorLabel.size()) result.colorLabel.clear();
    return result;
}

std::string pageLabel(const Payload& payload, std::size_t page) {
    if (page < payload.pageLabels.size()) return payload.pageLabels[page];
    return "Page " + std::to_string(page + 1);
}

}  // namespace

Payload comparePayloads(const std::vector<ComparisonEntry>& entries) {
    if (entries.size() < 2) throw std::invalid_argument("A comparison needs at least two runs");
    const bool colored = std::all_of(entries.begin(), entries.end(),
        [](const ComparisonEntry& entry) { return !entry.color.empty(); });
    std::vector<Pages> runPages;
    std::size_t count = 0;
    std::string sources;
    for (const auto& entry : entries) {
        runPages.push_back(pages(entry.payload));
        count = std::max(count, runPages.back().size());
        sources += (sources.empty() ? "" : "\n") + entry.payload.source;
    }
    Payload result{sources, "comparison", joinedLabels(entries), 2};
    for (const auto& entry : entries)
        for (const auto& warning : entry.payload.warnings)
            result.warnings.push_back(entry.label + ": " + warning);
    for (std::size_t page = 0; page < count; ++page) {
        std::string label;
        bool same = true;
        for (std::size_t k = 0; k < entries.size(); ++k) {
            const auto name = pageLabel(entries[k].payload, page);
            if (k == 0) label = name; else if (name != label) same = false;
        }
        if (!same) {
            label.clear();
            for (std::size_t k = 0; k < entries.size(); ++k)
                label += (k ? " | " : "") + pageLabel(entries[k].payload, page);
        }
        result.pageLabels.push_back(label);
        std::size_t figuresOnPage = 0;
        for (const auto& run : runPages) if (page < run.size()) figuresOnPage = std::max(figuresOnPage, run[page].size());
        for (std::size_t index = 0; index < figuresOnPage; ++index) {
            std::vector<const Figure*> figures;
            std::vector<ComparisonEntry> present;
            for (std::size_t k = 0; k < entries.size(); ++k)
                if (page < runPages[k].size() && index < runPages[k][page].size()) {
                    figures.push_back(runPages[k][page][index]);
                    present.push_back({entries[k].label, {}, entries[k].color});
                }
            const bool all = figures.size() == entries.size() && std::all_of(figures.begin() + 1, figures.end(),
                [&](const Figure* figure) { return overlayable(*figures.front(), *figure); });
            if (all) { result.figures.push_back(overlay(figures, present, static_cast<int>(page), colored)); continue; }
            for (std::size_t k = 0; k < figures.size(); ++k) {
                Figure figure = *figures[k];
                figure.title = present[k].label + " — " + figure.title;
                figure.dashboardPage = static_cast<int>(page); figure.gridSpan = 1;
                result.figures.push_back(std::move(figure));
            }
        }
    }
    return result;
}

Payload comparePayload(const Payload& first, const std::string& firstLabel,
                       const Payload& second, const std::string& secondLabel) {
    return comparePayloads({{firstLabel, first, {}}, {secondLabel, second, {}}});
}

}  // namespace dsl
