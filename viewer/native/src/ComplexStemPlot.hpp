#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <array>
#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace dsl {

// Standard Viewer plot: stems rising from points of the complex plane.
// x = Re(base), y = Im(base), z = height. The figure is a "three_d" panel with
// independent axis scaling that opens orthographic at azimuth -45, elevation 45
// (plane tilted 45 degrees out of the page, height up); dragging rotates it and
// reset returns to that view.
// An additional set of stems drawn over the main series (for example the modes a
// reduced model kept, drawn over the parent spectrum).
struct StemSeries {
    std::string label;
    std::string color;
    std::vector<std::complex<double>> bases;
    std::vector<double> heights;         // same length as bases
};

struct ComplexStemPlot {
    std::string title;
    std::string baseSymbol = "z";        // axis labels become Re(...) and Im(...)
    std::string heightLabel = "height";
    std::string seriesLabel = "stems";
    std::string color = "#0f766e";
    std::vector<std::complex<double>> bases;
    std::vector<double> heights;         // same length as bases
    bool unitCircle = false;             // dashed |z|=1 reference on the base plane
    bool imaginaryAxis = false;          // dashed Re(z)=0 reference (stability boundary of rates)
    // Optional window {Re min, Re max, Im min, Im max}: reference curves (unit circle,
    // Re = 0) are drawn only inside it, so a zoomed plane is not rescaled by them.
    std::optional<std::array<double, 4>> referenceWindow;
    // Store log10(height) in the data itself (stems from the decade below the smallest
    // height). Normally off: heights stay linear and the Viewer's value-scale control
    // switches the height axis between linear and log ("auto" picks log when the
    // heights span three or more decades).
    bool logHeights = false;
    std::vector<StemSeries> overlays;    // drawn after (on top of) the main stems
    // Curves over the base plane (height as a function of a path of bases, e.g. a density around the
    // unit circle), drawn as lines in order, after the stems.
    std::vector<StemSeries> curves;
};

Figure complexPlaneStemFigure(const ComplexStemPlot& plot);

}  // namespace dsl
