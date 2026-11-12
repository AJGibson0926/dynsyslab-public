#include "dynsyslab_viewer/PanelCatalog.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dsl::Figure figure(std::string kind, std::string title) {
    return {std::move(kind), std::move(title), "x", "y", false};
}

}  // namespace

int main() {
    try {
        dsl::CatalogItem item;
        item.views = {"data_overview", "control", "time", "phase", "fourier",
                      "control_fourier", "time_frequency"};
        const auto routes = dsl::panelSourceViews(item);
        require(routes == std::vector<std::string>({"time", "phase", "fourier",
            "control_fourier", "data_overview", "three_d", "observable_time"}), "panel route policy");

        dsl::CatalogItem staleTrajectory;
        staleTrajectory.views = {"data_overview"};
        const auto derivedRoutes = dsl::panelSourceViews(staleTrajectory);
        for (const auto* route : {"time", "phase", "three_d", "observable_time"})
            require(std::find(derivedRoutes.begin(),derivedRoutes.end(),route) != derivedRoutes.end(),
                "dashboard component routes must be available to background inspection");
        require(std::find(derivedRoutes.begin(), derivedRoutes.end(), "fourier") !=
                    derivedRoutes.end(),
                "state Fourier route must not depend on fresh catalogue metadata");
        require(std::find(derivedRoutes.begin(), derivedRoutes.end(),
                          "control_fourier") != derivedRoutes.end(),
                "control Fourier route must be probed from native record data");

        dsl::Payload time; time.figures.push_back(figure("time", "radius versus time"));
        dsl::Payload overview;
        overview.figures.push_back(figure("time", "radius versus time"));
        overview.figures.push_back(figure("phase", "Phase space"));
        dsl::Payload spectrum;
        spectrum.figures.push_back(figure("control_fourier", "u0 Fourier spectrum"));
        auto choices = dsl::composePanelChoices({
            {"time", std::move(time)}, {"data_overview", std::move(overview)},
            {"control_fourier", std::move(spectrum)}});
        require(choices.size() == 3, "duplicate panels were not removed");
        require(std::any_of(choices.begin(), choices.end(), [](const auto& choice) {
            return choice.family == "Control spectra" && choice.figure.gridSpan == 1;
        }), "control Fourier choice classification");
        std::cout << "dynsyslab_panel_catalog_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_panel_catalog_tests: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
