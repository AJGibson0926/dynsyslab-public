#include "ComplexHeatmap.hpp"

#include <algorithm>
#include <cmath>

namespace dsl {

bool validComplexDisplay(const std::string& display) {
    return display == "real-imaginary" || display == "magnitude-phase" || display == "domain-coloring";
}

double complexComponent(std::complex<double> value, const std::string& display, unsigned part) {
    return display == "real-imaginary" ? (part ? value.imag() : value.real()) : (part ? std::arg(value) : std::abs(value));
}

std::string complexComponentLabel(const std::string& name, const std::string& display, unsigned part) {
    if (display == "real-imaginary") return std::string(part ? "Im(" : "Re(") + name + ")";
    return part ? "arg(" + name + ") (rad)" : "|" + name + "|";
}

void encodeComplexHeatmap(Figure& figure, const std::string& display) {
    for (auto& trace : figure.traces) {
        if (!trace.heatmap || trace.complexValues.empty() || trace.complexPart < 0) continue;
        const auto part = static_cast<unsigned>(trace.complexPart);
        const auto& name = trace.complexName;
        trace.z.clear();
        trace.magnitudes.clear();
        if (display == "domain-coloring" && part == 0) {
            double scale = 0;
            for (const auto v : trace.complexValues) if (std::isfinite(std::abs(v))) scale = std::max(scale, std::abs(v));
            if (!(scale > 0)) scale = 1;
            for (const auto v : trace.complexValues) {
                trace.z.push_back(std::arg(v));
                trace.magnitudes.push_back(std::abs(v) / scale);
            }
            trace.encoding = "domain-color";
            figure.title = name + " — domain coloring";
            figure.colorLabel = "phase (rad); brightness = magnitude";
            continue;
        }
        // Split panels; in domain-coloring mode the second panel shows the magnitude.
        const auto shown = display == "domain-coloring" ? std::string("magnitude-phase") : display;
        const unsigned shownPart = display == "domain-coloring" ? 0u : part;
        for (const auto v : trace.complexValues) trace.z.push_back(complexComponent(v, shown, shownPart));
        trace.encoding = "scalar";
        figure.title = complexComponentLabel(name, shown, shownPart);
        figure.colorLabel = figure.title;
    }
}

}  // namespace dsl
