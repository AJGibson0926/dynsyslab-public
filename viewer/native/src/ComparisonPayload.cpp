#include "dynsyslab_viewer/ComparisonPayload.hpp"

#include <algorithm>

namespace dsl {
namespace {

constexpr std::size_t kFiguresPerPage = 4;

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

void appendTraces(Figure& target, const Figure& source,
                  const std::string& label, bool alternate) {
    if (target.alternateCoordinates && source.alternateCoordinates)
        for (const auto& traceLabel : source.alternateCoordinates->traceLabels)
            target.alternateCoordinates->traceLabels.push_back(label + " — " + traceLabel);
    for (auto trace : source.traces) {
        trace.label = label + " — " + trace.label;
        if (alternate && trace.lineStyle == "solid" && !trace.marker
            && trace.encoding != "complex-domain-color") trace.lineStyle = "dash";
        target.traces.push_back(std::move(trace));
    }
}

Figure overlay(const Figure& first, const std::string& firstLabel,
               const Figure& second, const std::string& secondLabel, int page) {
    Figure result = first;
    result.title = firstLabel + " vs " + secondLabel + " — " + first.title;
    result.titleColor.clear(); result.dashboardPage = page;
    result.traces.clear(); result.inset.reset();
    if (result.alternateCoordinates) result.alternateCoordinates->traceLabels.clear();
    result.references.insert(result.references.end(), second.references.begin(), second.references.end());
    appendTraces(result, first, firstLabel, false);
    appendTraces(result, second, secondLabel, true);
    return result;
}

std::string pageLabel(const Payload& payload, std::size_t page) {
    if (page < payload.pageLabels.size()) return payload.pageLabels[page];
    return "Page " + std::to_string(page + 1);
}

}  // namespace

Payload comparePayload(const Payload& first, const std::string& firstLabel,
                       const Payload& second, const std::string& secondLabel) {
    const auto firstPages = pages(first), secondPages = pages(second);
    const std::size_t count = std::max(firstPages.size(), secondPages.size());
    Payload result{first.source + "\n" + second.source, "comparison",
                   firstLabel + " vs " + secondLabel, 2};
    for (const auto& warning : first.warnings)
        result.warnings.push_back(firstLabel + ": " + warning);
    for (const auto& warning : second.warnings)
        result.warnings.push_back(secondLabel + ": " + warning);
    const std::vector<const Figure*> empty;
    for (std::size_t page = 0; page < count; ++page) {
        const auto& left = page < firstPages.size() ? firstPages[page] : empty;
        const auto& right = page < secondPages.size() ? secondPages[page] : empty;
        const std::string leftPage = pageLabel(first, page);
        const std::string rightPage = pageLabel(second, page);
        result.pageLabels.push_back(leftPage == rightPage ? leftPage
            : leftPage + " | " + rightPage);
        const std::size_t paired = std::min(left.size(), right.size());
        for (std::size_t index = 0; index < paired; ++index) {
            if (overlayable(*left[index], *right[index])) {
                result.figures.push_back(overlay(*left[index], firstLabel,
                    *right[index], secondLabel, static_cast<int>(page)));
            } else {
                Figure firstFigure = *left[index], secondFigure = *right[index];
                firstFigure.title = firstLabel + " — " + firstFigure.title;
                secondFigure.title = secondLabel + " — " + secondFigure.title;
                firstFigure.dashboardPage = secondFigure.dashboardPage = static_cast<int>(page);
                firstFigure.gridSpan = secondFigure.gridSpan = 1;
                result.figures.push_back(std::move(firstFigure));
                result.figures.push_back(std::move(secondFigure));
            }
        }
        for (std::size_t index = paired; index < left.size(); ++index) {
            Figure figure = *left[index]; figure.title = firstLabel + " — " + figure.title;
            figure.dashboardPage = static_cast<int>(page); figure.gridSpan = 1;
            result.figures.push_back(std::move(figure));
        }
        for (std::size_t index = paired; index < right.size(); ++index) {
            Figure figure = *right[index]; figure.title = secondLabel + " — " + figure.title;
            figure.dashboardPage = static_cast<int>(page); figure.gridSpan = 1;
            result.figures.push_back(std::move(figure));
        }
    }
    return result;
}

}  // namespace dsl
