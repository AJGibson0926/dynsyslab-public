#include "dynsyslab_viewer/CatalogTreeSort.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dsl::CatalogItem item(std::string id, std::string kind, std::string label,
                      std::string parent, double modified,
                      std::string source, double bytes) {
    dsl::CatalogItem result;
    result.id = std::move(id);
    result.kind = std::move(kind);
    result.label = std::move(label);
    result.parentId = std::move(parent);
    result.modified = modified;
    result.source = std::move(source);
    result.bytes = bytes;
    return result;
}

}  // namespace

int main() {
    try {
        dsl::Catalog catalog;
        catalog.items.push_back(item("a", "presentation_group", "Branch A", {}, 0, {}, 0));
        catalog.items.push_back(item("a-old", "trajectory", "IC #10", "a", 10,
            "archives/shared.npz", 100));
        catalog.items.push_back(item("a-new", "trajectory", "IC #2", "a", 20,
            "ARCHIVES\\shared.npz", 100));
        catalog.items.push_back(item("a-control", "control_attempt", "Attempt #1", "a", 15,
            "archives/control.npz", 300));
        catalog.items.push_back(item("b", "presentation_group", "Branch B", {}, 0, {}, 0));
        catalog.items.push_back(item("b-new", "trajectory", "IC #1", "b", 30,
            "archives/large.npz", 1000));
        catalog.items[1].properties["trajectory count"] = "2";
        catalog.items[2].properties["trajectory count"] = "10";
        catalog.items[3].properties["trajectory count"] = "5";

        const auto date = dsl::orderedCatalogTree(catalog, dsl::CatalogTreeSort::date);
        require(date.roots.size() == 2 && date.roots[0] == 0 && date.roots[1] == 4 &&
                date.children[0].size() == 3 && date.children[0][0] == 2 &&
                date.children[0][1] == 3 && date.children[0][2] == 1,
                "date sort changed a branch or failed to order leaves newest-first");
        require(date.aggregates[0].newest == 20 &&
                date.aggregates[0].trajectoryCount == 17 &&
                std::abs(date.aggregates[0].totalBytes - 400) < 0.5,
                "subtree date/count/size aggregation is incorrect or double-counts a source");

        const auto trajectories = dsl::orderedCatalogTree(
            catalog, dsl::CatalogTreeSort::trajectoryCount);
        require(trajectories.roots[0] == 0 && trajectories.roots[1] == 4 &&
                trajectories.children[0][0] == 2 && trajectories.children[0][1] == 3 &&
                trajectories.children[0][2] == 1,
                "trajectory-count sort changed a branch or is not highest-first at leaves");

        const auto size = dsl::orderedCatalogTree(catalog, dsl::CatalogTreeSort::totalBytes);
        require(size.roots[0] == 0 && size.roots[1] == 4 &&
                size.children[0][0] == 3 && size.children[0][1] == 2 &&
                size.children[0][2] == 1,
                "total-size leaf sort, branch order, or natural tie-break is incorrect");

        dsl::Catalog startup;
        auto other = item("other", "trajectory", "Data #1", {}, 100,
            "archives/other.npz", 500);
        other.views = {"phase"};
        startup.items.push_back(std::move(other));
        startup.items.push_back(item("training", "presentation_group",
            "Training data", {}, 0, {}, 0));
        auto olderTraining = item("training-old", "trajectory", "Training set #1",
            "training", 10, "archives/training-old.npz", 100);
        olderTraining.views = {"phase"};
        olderTraining.properties["physics"] = "Air in water · 20 C";
        startup.items.push_back(std::move(olderTraining));
        auto newerTraining = item("training-new", "trajectory", "Training set #2",
            "training", 20, "archives/training-new.npz", 200);
        newerTraining.views = {"phase"};
        newerTraining.properties["physics"] = "General parameters";
        startup.items.push_back(std::move(newerTraining));
        auto staleMetadata = item("stale", "trajectory", "Data #2", {}, 200,
            "archives/stale.npz", 900);
        staleMetadata.category = "Data";
        staleMetadata.properties["archive category"] = "Training data";
        staleMetadata.views = {"phase"};
        startup.items.push_back(std::move(staleMetadata));
        require(dsl::isTrainingDataItem(startup, 2) &&
                dsl::isTrainingDataItem(startup, 3) &&
                !dsl::isTrainingDataItem(startup, 0) &&
                !dsl::isTrainingDataItem(startup, 4) &&
                dsl::preferredInitialCatalogItem(startup,
                    dsl::CatalogTreeSort::date) == 2,
                "Viewer startup did not prefer a contextualized training set branch");

        dsl::Catalog withoutTraining;
        withoutTraining.items.push_back(std::move(startup.items.front()));
        require(dsl::preferredInitialCatalogItem(withoutTraining,
                    dsl::CatalogTreeSort::date) == -1,
                "a catalog without training data invented a training selection");
        std::cout << "catalog tree sort tests passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
