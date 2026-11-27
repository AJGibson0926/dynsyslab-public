#include "../src/CatalogPresentation.hpp"
#include "../src/ArchiveCatalog.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dsl::CatalogItem item(std::string id, std::string kind, std::string label,
                      std::string parent = {}) {
    dsl::CatalogItem result;
    result.id = std::move(id);
    result.kind = std::move(kind);
    result.category = "Archived solves";
    result.label = std::move(label);
    result.parentId = std::move(parent);
    result.system = "Rayleigh-Plesset";
    result.source = "archives/" + result.id + ".npz";
    result.properties["physics preset"] = "Air in ultrapure water";
    result.properties["temperature"] = "20 C";
    return result;
}

std::vector<std::string> ancestry(const dsl::Catalog& catalog,
                                  const std::string& id) {
    std::map<std::string, const dsl::CatalogItem*> byId;
    for (const auto& candidate : catalog.items) byId[candidate.id] = &candidate;
    std::vector<std::string> result;
    auto found = byId.find(id);
    while (found != byId.end() && !found->second->parentId.empty()) {
        found = byId.find(found->second->parentId);
        if (found != byId.end()) result.push_back(found->second->label);
    }
    return result;
}

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

void controlAttempt(dsl::Catalog& catalog, const std::string& id,
                    const std::string& krom, const std::string& representation) {
    auto attempt = item(id, "control_attempt",
        "rpe-approximate-bpsi-klqr-persistent-data-" + representation);
    attempt.properties["controller"] = "koopman-lqr";
    attempt.properties["control target"] = "rpe_multisine_equilibrium_identification";
    attempt.properties["control representation"] = representation;
    attempt.relations.push_back({krom, "controlled by KROM", "long KROM archive label"});
    catalog.items.push_back(std::move(attempt));
}

}  // namespace

int main() {
    try {
        dsl::CatalogItem dmdc;
        dsl::enrichArchiveCatalogItem(dmdc, {{"eigenvalues", "operator", "control_matrix",
            "projection_basis", "eigenfunction_coefficients"},
            R"({"method":"dmdc","snapshot_array":"measured_latents","schema":"dynsyslab.derived/v1"})"});
        require(dmdc.kind=="koopman_spectrum"&&dmdc.category=="Koopman spectra"
            &&dmdc.properties.at("snapshot array")=="measured_latents"
            &&dmdc.subtitle.find("not a control-ready KROM")!=std::string::npos,
            "DMDc is missing or misclassified as a control-ready KROM");
        dsl::CatalogItem edmdc;
        dsl::enrichArchiveCatalogItem(edmdc, {{"eigenvalues", "operator", "control_matrix",
            "projection_basis", "eigenfunction_coefficients", "polynomial_exponents", "state_center", "state_scale"},
            R"({"method":"edmdc","snapshot_array":"measured_latents","schema":"dynsyslab.derived/v1"})"});
        require(edmdc.kind=="koopman_spectrum"&&edmdc.properties.at("eigenfunction coordinates").find("Polynomials")!=std::string::npos
            &&edmdc.subtitle.find("not a control-ready KROM")!=std::string::npos,
            "EDMDc polynomial model misclassified");
        dsl::CatalogItem legacy;
        legacy.id = "legacy-control";
        legacy.kind = "trajectory";
        legacy.category = "Training data";
        dsl::enrichArchiveCatalogItem(legacy, {{"time", "values", "control_values"},
            R"({"label":"legacy-rpe-klqr-bilinear-long-name","system":"rayleigh-plesset",
                "forcing_type":"forced",
                "controller":{"controller":"koopman-lqr","control_representation":"bilinear",
                    "target":"rpe_equilibrium","target_space":"state",
                    "rom":{"source_archive":"C:/archives/krom-a.npz"}},
                "system_parameters":{"preset":{"key":"air-water-20c",
                    "description":"Air in ultrapure water","gas":"air","temperature_c":20}}})"});
        require(legacy.kind == "control_attempt" &&
                legacy.properties.at("physics preset") == "Air in ultrapure water" &&
                legacy.properties.at("forcing type") == "forced" &&
                legacy.properties.at("control representation") == "bilinear" &&
                legacy.properties.at("control target") == "rpe_equilibrium" &&
                !legacy.relations.empty(),
                "legacy archive grouping metadata was not recovered from its compact summary");

        dsl::Catalog catalog;
        catalog.projectName = "fixture";

        auto krom1 = item("krom-a", "krom",
            "rpe-krom-source-with-solver-and-control-details-hdmd2000");
        krom1.properties.erase("physics preset");
        krom1.properties.erase("temperature");
        krom1.properties["created UTC"] = "2026-01-01";
        krom1.relations.push_back({"training", "derived from", "training data"});
        auto krom2 = item("krom-b", "krom",
            "rpe-krom-source-with-solver-and-control-details-hdmd500");
        krom2.properties.erase("physics preset");
        krom2.properties.erase("temperature");
        krom2.properties["created UTC"] = "2026-01-02";
        krom2.relations.push_back({"training", "derived from", "training data"});
        catalog.items.push_back(std::move(krom1));
        catalog.items.push_back(std::move(krom2));
        controlAttempt(catalog, "attempt-a", "krom-a", "constant");
        controlAttempt(catalog, "attempt-b", "krom-a", "bilinear");
        controlAttempt(catalog, "attempt-c", "krom-b", "bilinear");

        auto persistent = item("bpsi-data-a", "control_data",
            "rpe-bpsi-persistent-excitation-full-solver-configuration");
        persistent.relations.push_back({"krom-a", "uses KROM", "source KROM"});
        auto symmetric = item("bpsi-data-b", "control_data",
            "rpe-bpsi-symmetric-pulses-full-solver-configuration");
        symmetric.relations.push_back({"krom-a", "uses KROM", "source KROM"});
        catalog.items.push_back(std::move(persistent));
        catalog.items.push_back(std::move(symmetric));

        catalog.items.push_back(item("batch", "batch",
            "rpe-periodically-forced-high-accuracy-bdf-run"));
        catalog.items.push_back(item("system-run", "system_run",
            "Rayleigh-Plesset air-water exact parameter listing", "batch"));
        catalog.items.push_back(item("cluster", "cluster",
            "initial-condition-cluster-near-equilibrium", "system-run"));
        catalog.items.push_back(item("trajectory", "trajectory",
            "trajectory-rpe-periodic-forcing-rtol-1e-11", "cluster"));
        auto training = item("training", "trajectory",
            "rpe-training-data-full-solver-configuration");
        training.properties["forcing type"] = "free";
        catalog.items.push_back(std::move(training));
        dsl::CatalogItem timeGrid;
        timeGrid.id = "time-grid";
        timeGrid.kind = "trajectory";
        timeGrid.category = "Training data";
        timeGrid.label = std::string(64, 'a');
        timeGrid.source = "archives/time_grids/" + timeGrid.label + ".npz";
        catalog.items.push_back(std::move(timeGrid));
        auto error = item("bpsi-error", "control_data", "control-map-errors");
        error.source = "archives/experiments/study/attempt-9/control-map-errors.npz";
        error.relations.push_back({"krom-a", "uses KROM", "source KROM"});
        catalog.items.push_back(std::move(error));
        auto forcedTraining = item("forced-training", "trajectory",
            "rpe-periodic-training-data");
        forcedTraining.properties["forcing type"] = "forced";
        forcedTraining.relations.push_back(
            {"krom-a", "produced derived result", "trained KROM"});
        catalog.items.push_back(std::move(forcedTraining));

        dsl::applyCatalogPresentation(catalog);
        const auto firstSize = catalog.items.size();
        require(firstSize > 14, "presentation hierarchy was not created");
        require(catalog.items[0].label == "KROM #1" &&
                catalog.items[1].label == "KROM #2",
                "KROMs do not have compact stable ordering");
        require(catalog.items[0].properties.at("archive label").find("solver") !=
                std::string::npos, "original archive label was not preserved");
        require(catalog.items[2].label == "Attempt #1" &&
                catalog.items[2].label.size() < 24 &&
                catalog.items[2].system == "Bubble · Rayleigh-Plesset",
                "control attempt label is not compact");
        require(catalog.items[5].label == "Dataset #1" &&
                catalog.items[6].label == "Dataset #1" &&
                catalog.items[7].label == "Run #1" &&
                catalog.items[8].label == "Parameter set #1" &&
                catalog.items[9].label == "IC set #1" &&
                catalog.items[10].label == "IC #1",
                "terminal records retain contextual descriptors in their labels");

        const auto firstKrom = ancestry(catalog, "krom-a");
        require(contains(firstKrom, "Models & assets") &&
                contains(firstKrom, "H-DMD 2000"),
                "KROM model was not grouped with model assets and method context");
        const auto trainingTree = ancestry(catalog, "training");
        require(catalog.items[11].label == "Training set #1" &&
                contains(trainingTree, "Training data") &&
                contains(trainingTree, "Free"),
                "a free KROM source was not classified and named as training data");
        const auto forcedTrainingTree = ancestry(catalog, "forced-training");
        require(catalog.items[14].label == "Training set #1" &&
                contains(forcedTrainingTree, "Training data") &&
                contains(forcedTrainingTree, "Forced") &&
                !contains(forcedTrainingTree, "Free"),
                "forced and free KROM training data were not separated");
        const auto unspecifiedTrainingTree = ancestry(catalog, "time-grid");
        require(catalog.items[12].label == "Training set #1" &&
                contains(unspecifiedTrainingTree, "Training data") &&
                contains(unspecifiedTrainingTree, "Unspecified"),
                "training-category records without forcing metadata were not consistently named");

        const auto firstControl = ancestry(catalog, "attempt-a");
        const auto secondKromControl = ancestry(catalog, "attempt-c");
        require(contains(firstControl, "Air in water · 20 C") &&
                contains(firstControl, "Control & estimation") &&
                contains(firstControl, "KLQR controls") &&
                contains(firstControl, "KROM #1 · H-DMD 2000") &&
                contains(firstControl, "Target · EC") &&
                contains(firstControl, "Constant fit"),
                "KLQR hierarchy lacks physics, KROM, target, or strategy grouping");
        require(contains(secondKromControl, "KROM #2 · H-DMD 500") &&
                !contains(secondKromControl, "KROM #1 · H-DMD 2000"),
                "same target under a second KROM was merged into the first KROM");
        const auto persistentTree = ancestry(catalog, "bpsi-data-a");
        const auto symmetricTree = ancestry(catalog, "bpsi-data-b");
        require(contains(persistentTree, "Models & assets") &&
                contains(persistentTree, "Bψ identification") &&
                contains(persistentTree, "Persistent-excitation data") &&
                contains(symmetricTree, "Symmetric-pulse data"),
                "Bψ comparison datasets are not grouped by evidence type");
        const auto errorTree = ancestry(catalog, "bpsi-error");
        require(contains(errorTree, "Search attempt #9") &&
                contains(errorTree, "Error fields") &&
                catalog.items[12].system == "Other" &&
                catalog.items[13].label == "Error #1",
                "attempt products or unknown-system records are not compactly grouped");
        require(catalog.items[2].relations.front().label == catalog.items[0].label,
                "relation labels did not follow compact presentation names");
        for (const auto& entry : catalog.items)
            if (entry.kind != "presentation_group") {
                require(entry.category == "Data", "scientific items did not share Data root");
                require(entry.label.size() <= 20 && entry.label.find(" · ") == std::string::npos,
                        "a data-record label retained contextual detail");
            }

        dsl::applyCatalogPresentation(catalog);
        require(catalog.items.size() == firstSize,
                "presentation transform is not idempotent");
        std::cout << "catalog presentation tests passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
