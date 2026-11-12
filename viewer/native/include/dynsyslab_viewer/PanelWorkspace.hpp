#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace dsl {

struct PanelSpan {
    int columns = 1;
    int rows = 1;
};

struct PanelGridPlacement {
    std::size_t panel = 0;
    int column = 0;
    int row = 0;
    int columnSpan = 1;
    int rowSpan = 1;
};

// Chooses a compact grid while preserving useful scientific plot aspect ratios.
int adaptivePanelColumns(std::size_t panelCount, int pixelWidth, int pixelHeight);

// First-fit, stable packing. Enlarged panels reserve cells and later panels
// move around them; no scientific object or archive is changed.
std::vector<PanelGridPlacement> packPanelGrid(
    const std::vector<PanelSpan>& panels, int columns);

int packedPanelRows(const std::vector<PanelGridPlacement>& placements);

template <class T>
bool movePanel(std::vector<T>& panels, std::size_t from, std::size_t to) {
    if (from >= panels.size() || to >= panels.size() || from == to) return false;
    T moving = std::move(panels[from]);
    panels.erase(panels.begin() + static_cast<std::ptrdiff_t>(from));
    panels.insert(panels.begin() + static_cast<std::ptrdiff_t>(to), std::move(moving));
    return true;
}

}  // namespace dsl
