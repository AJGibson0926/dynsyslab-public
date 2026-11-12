#include "dynsyslab_viewer/PanelWorkspace.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    try {
        require(dsl::adaptivePanelColumns(1, 1400, 800) == 1,
            "one panel must occupy one adaptive column");
        require(dsl::adaptivePanelColumns(6, 1400, 800) >= 2,
            "six panels did not adapt to a multi-column workspace");
        require(dsl::adaptivePanelColumns(2, 1200, 700) == 2,
            "two panels did not share a wide workspace");
        require(dsl::adaptivePanelColumns(4, 300, 900) == 1,
            "narrow workspace did not collapse to one column");

        const std::vector<dsl::PanelSpan> spans{{2, 1}, {1, 2}, {1, 1}, {2, 1}};
        const auto layout = dsl::packPanelGrid(spans, 2);
        require(layout.size() == spans.size(), "panel pack lost an item");
        require(layout[0].column == 0 && layout[0].row == 0 &&
            layout[0].columnSpan == 2, "full-row panel placement");
        require(layout[1].column == 0 && layout[1].row == 1 &&
            layout[1].rowSpan == 2, "column-spanning panel placement");
        require(layout[2].column == 1 && layout[2].row == 1,
            "neighbor did not flow beside tall panel");
        require(layout[3].row == 3 && dsl::packedPanelRows(layout) == 4,
            "collision did not displace the following full-row panel");

        std::vector<int> order{0, 1, 2, 3};
        require(dsl::movePanel(order, 0, 2), "panel move was rejected");
        require(order == std::vector<int>({1, 2, 0, 3}),
            "panel move did not preserve relative order");
        std::cout << "dynsyslab_panel_workspace_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_panel_workspace_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
