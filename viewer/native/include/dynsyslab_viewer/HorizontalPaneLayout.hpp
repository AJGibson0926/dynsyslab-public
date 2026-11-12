#pragma once

namespace dsl {

struct HorizontalPaneConstraints {
    int totalWidth = 0;
    int margin = 12;
    int gap = 12;
    int desiredNavigationWidth = 0;
    int desiredInspectorWidth = 0;
    int defaultNavigationWidth = 285;
    int defaultInspectorWidth = 300;
    int minimumNavigationWidth = 210;
    int maximumNavigationWidth = 480;
    int minimumCenterWidth = 420;
    int minimumInspectorWidth = 250;
    int maximumInspectorWidth = 480;
    int inspectorThreshold = 1180;
};

struct HorizontalPaneLayout {
    int navigationWidth = 0;
    int centerX = 0;
    int centerWidth = 0;
    int inspectorX = 0;
    int inspectorWidth = 0;
    int leftSplitterCenter = 0;
    int rightSplitterCenter = 0;
    bool inspectorVisible = false;
};

HorizontalPaneLayout resolveHorizontalPaneLayout(
    const HorizontalPaneConstraints& constraints);

int navigationWidthFromSplitter(int splitterCenterX, int margin, int gap);
int inspectorWidthFromSplitter(
    int splitterCenterX, int totalWidth, int margin, int gap);

}  // namespace dsl
