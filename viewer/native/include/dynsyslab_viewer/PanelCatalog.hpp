#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <string>
#include <utility>
#include <vector>

namespace dsl {

struct PanelChoice {
    std::string route;
    std::string family;
    std::string label;
    std::string detail;
    Figure figure;
};

// Orders the smallest useful route set for discovering panels in one archive
// object. Compound routes are excluded because their constituent plots are
// already represented independently.
std::vector<std::string> panelSourceViews(const CatalogItem& item);

// Flattens independently loaded payloads and removes duplicate scientific
// panels. This is intentionally independent of Win32 and archive I/O.
std::vector<PanelChoice> composePanelChoices(
    std::vector<std::pair<std::string, Payload>> payloads);

}  // namespace dsl
