#pragma once
// Vector (SVG) and self-contained interactive HTML export of Viewer payload pages.
#include "dynsyslab_viewer/Model.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dsl {

// A trace in the figure's drawing plane: (x, y), or (x, y, z) projected orthographically at the
// figure's camera (azimuth, elevation; default 42 and 24 degrees) for three-state views.
struct Projected { std::vector<double> x, y; };
[[nodiscard]] Projected projectedTrace(const Figure& figure, const Trace& trace);
// A math-text plot label (_{sub}, ^{sup}, \frac{a}{b}, Greek names) as plain Unicode text.
[[nodiscard]] std::string plainLabel(const std::string& text);
// Figures on a dashboard page (all figures when the payload has no pages), and the page count.
[[nodiscard]] std::vector<std::size_t> pageFigures(const Payload& payload, int page);
[[nodiscard]] int pageCount(const Payload& payload);
// One figure drawn into the rectangle (SVG fragment, no <svg> root).
[[nodiscard]] std::string figureSvgBody(const Figure& figure, double left, double top, double width, double height);
// A page as a standalone SVG document, figures tiled payload.columns per row.
[[nodiscard]] std::string payloadSvg(const Payload& payload, int page, int width, int height);
// Every page of the payload in one self-contained HTML file: the traces as JSON and a small
// script that draws them on a canvas with hover readout, wheel zoom, drag pan, double-click reset
// and page tabs (PlotHtml.cpp). Opens offline in any browser.
[[nodiscard]] std::string payloadHtml(const Payload& payload, const std::string& title);

}  // namespace dsl
