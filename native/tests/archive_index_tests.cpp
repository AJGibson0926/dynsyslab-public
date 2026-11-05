#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
}

int main() {
    const auto root = std::filesystem::temp_directory_path() /
        "dynsyslab-archive-index-tests";
    std::error_code ignored; std::filesystem::remove_all(root, ignored);
    const auto source = root / "solutions" / "rpe" / "trajectory.npz";
    dynsyslab::core::NpzRecord record;
    record.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({2}, {0.0, 1.0});
    record.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {2, 2}, {1.0, 0.0, 1.1, 0.2});
    record.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"system":"rayleigh-plesset","integration_form":"dimensionless",
            "controller":{"controller":"koopman-lqr","control_representation":"bilinear",
                "target":"rpe_equilibrium","target_space":"state",
                "rom":{"source_archive":"krom.npz"},"large_gain_matrix":[1,2,3]},
            "system_parameters":{"preset":{"key":"air-water-20c",
                "description":"Air in ultrapure water","gas":"air","temperature_c":20}}})");
    dynsyslab::core::writeNpzAtomic(source, record);

    const auto first = dynsyslab::core::reconcileArchiveIndex(root);
    require(first.discovered == 1 && first.indexed == 1, "initial indexing");
    const auto firstGeneration = dynsyslab::core::archiveIndexGeneration(root);
    require(firstGeneration == 1, "initial archive generation");
    const auto entries = dynsyslab::core::loadArchiveIndex(root);
    require(entries.size() == 1 && entries[0].sha256.size() == 64,
            "indexed digest and pointer");
    require(entries[0].summaryJson.find("rayleigh-plesset") != std::string::npos,
            "indexed compact summary");
    require(entries[0].summaryJson.find("rpe_equilibrium") != std::string::npos &&
            entries[0].summaryJson.find("Air in ultrapure water") != std::string::npos &&
            entries[0].summaryJson.find("large_gain_matrix") == std::string::npos,
            "compact summary lost grouping semantics or retained detailed solver data");
    require(dynsyslab::core::storeArchiveCatalogSnapshot(
        root, "archive-index-test", 7, firstGeneration,
        R"({"ready":true})"), "catalog snapshot storage");
    const auto snapshot = dynsyslab::core::loadArchiveCatalogSnapshot(
        root, "archive-index-test", 7);
    require(snapshot && snapshot->find("ready") != std::string::npos,
            "current catalog snapshot");
    const auto second = dynsyslab::core::reconcileArchiveIndex(root);
    require(second.unchanged == 1 && second.indexed == 0,
            "unchanged archive is not reopened or rehashed");
    require(dynsyslab::core::archiveIndexGeneration(root) == firstGeneration &&
            dynsyslab::core::loadArchiveCatalogSnapshot(
                root, "archive-index-test", 7).has_value(),
            "unchanged reconciliation preserves the catalog snapshot");

    std::filesystem::remove(source);
    const auto third = dynsyslab::core::reconcileArchiveIndex(root);
    require(third.removed == 1 && dynsyslab::core::loadArchiveIndex(root).empty(),
            "deleted pointer removal");
    require(dynsyslab::core::archiveIndexGeneration(root) == firstGeneration + 1 &&
            !dynsyslab::core::loadArchiveCatalogSnapshot(
                root, "archive-index-test", 7).has_value(),
            "archive changes invalidate derived catalog snapshots");
    std::filesystem::remove_all(root, ignored);
    std::cout << "archive index tests passed\n";
    return 0;
}
