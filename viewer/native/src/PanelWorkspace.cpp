#include "dynsyslab_viewer/PanelWorkspace.hpp"

#include <algorithm>
#include <cmath>

namespace dsl {

int adaptivePanelColumns(std::size_t count, int width, int height) {
    if (count <= 1) return 1;
    const double aspect = static_cast<double>(std::max(width, 1)) /
        static_cast<double>(std::max(height, 1));
    const int geometric = static_cast<int>(std::lround(
        std::sqrt(static_cast<double>(count) * aspect / 1.45)));
    const int byWidth = std::max(1, (std::max(width, 1) + 16) / 340);
    if (count == 2 && byWidth >= 2 && aspect >= 1.05) return 2;
    return std::clamp(geometric, 1,
        std::min({4, byWidth, static_cast<int>(count)}));
}

std::vector<PanelGridPlacement> packPanelGrid(
    const std::vector<PanelSpan>& panels, int requestedColumns) {
    const int columns = std::max(1, requestedColumns);
    std::vector<std::vector<bool>> occupied;
    auto ensureRows = [&](int count) {
        while (static_cast<int>(occupied.size()) < count)
            occupied.emplace_back(static_cast<std::size_t>(columns), false);
    };
    auto fits = [&](int row, int column, int rows, int span) {
        ensureRows(row + rows);
        if (column + span > columns) return false;
        for (int y = row; y < row + rows; ++y)
            for (int x = column; x < column + span; ++x)
                if (occupied[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)])
                    return false;
        return true;
    };

    std::vector<PanelGridPlacement> result;
    result.reserve(panels.size());
    for (std::size_t index = 0; index < panels.size(); ++index) {
        const int span = std::clamp(panels[index].columns, 1, columns);
        const int rows = std::clamp(panels[index].rows, 1, 8);
        int selectedRow = 0, selectedColumn = 0;
        bool placed = false;
        for (int row = 0; !placed; ++row) {
            for (int column = 0; column + span <= columns; ++column) {
                if (!fits(row, column, rows, span)) continue;
                selectedRow = row; selectedColumn = column; placed = true; break;
            }
        }
        for (int row = selectedRow; row < selectedRow + rows; ++row)
            for (int column = selectedColumn; column < selectedColumn + span; ++column)
                occupied[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] = true;
        result.push_back({index, selectedColumn, selectedRow, span, rows});
    }
    return result;
}

int packedPanelRows(const std::vector<PanelGridPlacement>& placements) {
    int rows = 0;
    for (const auto& item : placements) rows = std::max(rows, item.row + item.rowSpan);
    return std::max(rows, 1);
}

}  // namespace dsl
