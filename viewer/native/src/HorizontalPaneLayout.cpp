#include "dynsyslab_viewer/HorizontalPaneLayout.hpp"

#include <algorithm>

namespace dsl {
namespace {

int bounded(int value, int minimum, int maximum) {
    return std::clamp(value, minimum, std::max(minimum, maximum));
}

}  // namespace

HorizontalPaneLayout resolveHorizontalPaneLayout(
    const HorizontalPaneConstraints& constraints) {
    HorizontalPaneLayout result;
    const int width = std::max(0, constraints.totalWidth);
    const int margin = std::max(0, constraints.margin);
    const int gap = std::max(0, constraints.gap);
    const int navigation = constraints.desiredNavigationWidth > 0
        ? constraints.desiredNavigationWidth : constraints.defaultNavigationWidth;
    const int inspector = constraints.desiredInspectorWidth > 0
        ? constraints.desiredInspectorWidth : constraints.defaultInspectorWidth;

    const int threePaneMinimum = 2 * margin + 2 * gap +
        constraints.minimumNavigationWidth + constraints.minimumCenterWidth +
        constraints.minimumInspectorWidth;
    result.inspectorVisible = width >= constraints.inspectorThreshold &&
        width >= threePaneMinimum;

    if (result.inspectorVisible) {
        const int available = std::max(0, width - 2 * margin - 2 * gap);
        const int navigationMaximum = std::min(constraints.maximumNavigationWidth,
            available - constraints.minimumCenterWidth - constraints.minimumInspectorWidth);
        result.navigationWidth = bounded(navigation,
            constraints.minimumNavigationWidth, navigationMaximum);
        const int inspectorMaximum = std::min(constraints.maximumInspectorWidth,
            available - constraints.minimumCenterWidth - result.navigationWidth);
        result.inspectorWidth = bounded(inspector,
            constraints.minimumInspectorWidth, inspectorMaximum);
    } else {
        const int available = std::max(0, width - 2 * margin - gap);
        result.navigationWidth = bounded(navigation,
            constraints.minimumNavigationWidth,
            std::min(constraints.maximumNavigationWidth,
                available - constraints.minimumCenterWidth));
    }

    result.centerX = margin + result.navigationWidth + gap;
    result.centerWidth = std::max(0, width - result.centerX - margin -
        (result.inspectorVisible ? result.inspectorWidth + gap : 0));
    result.inspectorX = result.centerX + result.centerWidth + gap;
    result.leftSplitterCenter = margin + result.navigationWidth + gap / 2;
    result.rightSplitterCenter = result.centerX + result.centerWidth + gap / 2;
    return result;
}

int navigationWidthFromSplitter(int splitterCenterX, int margin, int gap) {
    return splitterCenterX - margin - gap / 2;
}

int inspectorWidthFromSplitter(
    int splitterCenterX, int totalWidth, int margin, int gap) {
    return totalWidth - margin - gap / 2 - splitterCenterX;
}

}  // namespace dsl
