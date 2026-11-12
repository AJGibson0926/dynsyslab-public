#include "dynsyslab_viewer/HorizontalPaneLayout.hpp"

#include <cassert>

int main() {
    dsl::HorizontalPaneConstraints constraints;
    constraints.totalWidth = 1480;
    constraints.desiredNavigationWidth = 285;
    constraints.desiredInspectorWidth = 300;
    constraints.maximumNavigationWidth = 500;
    constraints.maximumInspectorWidth = 520;
    const auto standard = dsl::resolveHorizontalPaneLayout(constraints);
    assert(standard.inspectorVisible);
    assert(standard.navigationWidth == 285);
    assert(standard.inspectorWidth == 300);
    assert(standard.centerWidth == 847);
    assert(dsl::navigationWidthFromSplitter(
        standard.leftSplitterCenter, constraints.margin, constraints.gap) == 285);
    assert(dsl::inspectorWidthFromSplitter(standard.rightSplitterCenter,
        constraints.totalWidth, constraints.margin, constraints.gap) == 300);
    const int stackedLeftWidth = standard.navigationWidth + constraints.gap +
        standard.inspectorWidth;
    const int shiftedPlotX = constraints.margin + stackedLeftWidth + constraints.gap;
    assert(stackedLeftWidth == 597);
    assert(shiftedPlotX == 621);
    assert(constraints.totalWidth - shiftedPlotX - constraints.margin ==
        standard.centerWidth);

    constraints.desiredNavigationWidth = 999;
    constraints.desiredInspectorWidth = 999;
    const auto bounded = dsl::resolveHorizontalPaneLayout(constraints);
    assert(bounded.navigationWidth == 500);
    assert(bounded.inspectorWidth == 512);
    assert(bounded.centerWidth == constraints.minimumCenterWidth);

    constraints.totalWidth = 1000;
    const auto compact = dsl::resolveHorizontalPaneLayout(constraints);
    assert(!compact.inspectorVisible);
    assert(compact.navigationWidth == 500);
    assert(compact.centerWidth == 464);

    constraints.totalWidth = 1200;
    constraints.desiredNavigationWidth = 300;
    constraints.desiredInspectorWidth = 310;
    constraints.minimumCenterWidth = 520;
    constraints.inspectorThreshold = 1200;
    const auto compute = dsl::resolveHorizontalPaneLayout(constraints);
    assert(compute.inspectorVisible);
    assert(compute.navigationWidth == 300);
    assert(compute.inspectorWidth == 310);
    assert(compute.centerWidth == 542);
    return 0;
}
