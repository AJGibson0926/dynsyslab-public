#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <string>
#include <vector>

namespace dsl {

// One run in a comparison. A non-empty color (#RRGGBB) recolors that run's own
// curves (solution, control, spectrum, estimate) and keeps them solid; estimates
// stay dashed. Without colors every run after the first is dashed instead.
struct ComparisonEntry {
    std::string label;
    Payload payload;
    std::string color;
};

// Combines compatible Viewer payloads under one synchronized page model. Figures at
// the same page position are overlaid when their kinds and axes agree, otherwise they
// are placed side by side. Targets identical in every run are drawn once, unprefixed.
[[nodiscard]] Payload comparePayloads(const std::vector<ComparisonEntry>& entries);

[[nodiscard]] Payload comparePayload(const Payload& first, const std::string& firstLabel,
                                     const Payload& second, const std::string& secondLabel);

}  // namespace dsl
