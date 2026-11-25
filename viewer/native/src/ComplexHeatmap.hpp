#pragma once
// Complex values on a 2-D domain as heat maps. A complex scalar on a grid is drawn as a pair of
// panels whose encoding can be switched without re-evaluating the grid (the values stay on the
// trace, Trace::complexValues):
//   domain-coloring   panel 0: hue = arg, brightness = |g| (DomainColor.hpp); panel 1: |g|
//   real-imaginary    panel 0: Re g;  panel 1: Im g
//   magnitude-phase   panel 0: |g|;   panel 1: arg g
// Domain coloring is the default for complex scalars on 2-D domains (Wegert, Visual Complex
// Functions, Birkhauser 2012).
#include "dynsyslab_viewer/Model.hpp"

#include <complex>
#include <string>

namespace dsl {

inline constexpr const char* kDefaultComplexDisplay = "domain-coloring";

[[nodiscard]] bool validComplexDisplay(const std::string& display);
// Scalar shown in a split panel (part 0 or 1) for real-imaginary or magnitude-phase.
[[nodiscard]] double complexComponent(std::complex<double> value, const std::string& display, unsigned part);
[[nodiscard]] std::string complexComponentLabel(const std::string& name, const std::string& display, unsigned part);
// Re-encode every heat-map trace that carries complex values for `display` (titles, colour labels,
// z, magnitudes, encoding). `name` is the function's display name, e.g. "\\psi_{3}".
void encodeComplexHeatmap(Figure& figure, const std::string& display);

}  // namespace dsl
