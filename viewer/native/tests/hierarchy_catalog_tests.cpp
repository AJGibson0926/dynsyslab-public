#include "dynsyslab_viewer/Model.hpp"
#include "../src/HierarchyPayload.hpp"
#include "../src/HierarchyCatalog.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

dynsyslab::core::BatchRecord batchFor(const std::filesystem::path& solution) {
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "duffing";
    configuration.parameters = {{"damping", 0.02}};
    dynsyslab::core::TrajectoryRecord trajectory;
    trajectory.id = "trajectory";
    trajectory.label = "portable trajectory";
    trajectory.initialState = {1.0, 0.0};
    trajectory.solutionArchive = solution;
    trajectory.tags = {"free"};
    dynsyslab::core::ClusterRecord cluster;
    cluster.id = "cluster";
    cluster.trajectories.push_back(trajectory);
    dynsyslab::core::SystemRunRecord system;
    system.id = "system";
    system.systemIdentifier = "duffing";
    system.systemLabel = "Duffing oscillator";
    system.systemFamily = "ode";
    system.preset = "default";
    system.presetLabel = "Default";
    system.parameterSetKey = dynsyslab::core::parameterSetKey(configuration);
    system.coordinateSystem = "dimensionless";
    system.resolvedParameters = configuration.parameters;
    system.parameterUnits = {{"damping", "dimensionless"}};
    system.dimensionlessGroups = {{"damping_ratio", 0.01}};
    system.dimensionlessGroupUnits = {{"damping_ratio", "dimensionless"}};
    system.equationsVersion = "duffing/test";
    system.equationsFingerprint = system.parameterSetKey;
    system.stateNames = {"x", "velocity"};
    system.clusters.push_back(cluster);
    dynsyslab::core::BatchRecord batch;
    batch.id = "batch";
    batch.label = "portable batch";
    batch.dataRole = "campaign";
    batch.campaignId = "campaign";
    batch.campaignKind = "frequency_sweep";
    batch.tags = {"free", "campaign"};
    batch.systems.push_back(system);
    return batch;
}

void writeSolution(const std::filesystem::path& path) {
    std::filesystem::create_directories(path.parent_path());
    dynsyslab::core::NpzRecord record;
    record.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({2}, {0.0, 1.0});
    record.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {2, 2}, {1.0, 0.0, 0.9, -0.1});
    record.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"format":"dynsyslab.solution/v1","system":"duffing"})");
    dynsyslab::core::writeNpzAtomic(path, record);
}

void portableHierarchyReference() {
    const auto oldWorkingDirectory = std::filesystem::current_path();
    const auto root = std::filesystem::temp_directory_path() /
        "dynsyslab-native-portable-hierarchy";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    const auto archiveRoot = root / "archives";
    const auto solution = archiveRoot / "solutions" / "duffing" / "trajectory.npz";
    const auto manifest = archiveRoot / "hierarchies" / "duffing" / "batch.hierarchy.json";
    writeSolution(solution);
    const auto draft = archiveRoot / "tentative-job.json";
    std::filesystem::create_directories(draft.parent_path());
    { std::ofstream output(draft);
      output << R"({"schema":"dynsyslab.compute-job/v1","kind":"solve_set","label":"tentative","batches":[]})"; }
    const auto progress = archiveRoot / "job_manager" / "progress" / "tentative.json";
    std::filesystem::create_directories(progress.parent_path());
    { std::ofstream output(progress); output << R"({"status":"queued"})"; }
    const auto genericRecord = archiveRoot / "legacy-record.json";
    { std::ofstream output(genericRecord); output << R"({"status":"complete"})"; }
    const auto legacyReport = root / "viewer_outputs" / "legacy-report.html";
    std::filesystem::create_directories(legacyReport.parent_path());
    { std::ofstream output(legacyReport); output << "<html><body>legacy</body></html>"; }
    std::filesystem::current_path(root.parent_path());
    const auto relativeSolution = std::filesystem::relative(solution, root.parent_path());
    dynsyslab::core::writeHierarchyAtomic(manifest, batchFor(relativeSolution));
    std::filesystem::current_path(oldWorkingDirectory);
    const auto encoded = dynsyslab::core::json::load(manifest.wstring())
        .at("systems").array().front().at("clusters").array().front()
        .at("trajectories").array().front().at("solution_archive").string();
    require(encoded.rfind("..", 0) == 0,
            "new hierarchy did not persist a portable manifest-relative archive path");
    const auto legacy = dsl::resolveHierarchySource(root, archiveRoot, manifest,
        std::filesystem::relative(solution, root).generic_u8string());
    require(std::filesystem::equivalent(legacy, solution),
        "legacy project-relative hierarchy reference is no longer readable");
    const auto catalog = dsl::Catalog::scan(root, archiveRoot);
    require(std::none_of(catalog.items.begin(), catalog.items.end(), [&](const auto& item) {
        const auto source = std::filesystem::u8path(item.source);
        return source == draft || source == progress || source == genericRecord ||
            source == legacyReport;
    }), "non-renderable records or reports leaked into the scientific data tree");
    const int index = catalog.indexOf("trajectory");
    require(index >= 0, "hierarchy trajectory is absent from native catalog");
    require(std::filesystem::equivalent(
        std::filesystem::u8path(catalog.items[static_cast<std::size_t>(index)].source),
        solution), "portable hierarchy reference did not resolve to its solution");
    const int batchIndex = catalog.indexOf("batch");
    require(batchIndex >= 0, "batch container is absent from native catalog");
    const auto& batch = catalog.items[static_cast<std::size_t>(batchIndex)];
    require(batch.kind == "batch" && batch.category == "Archived solves" &&
        batch.parentId.empty(), "completed batch is not an archive-tree root");
    require(std::find(batch.views.begin(), batch.views.end(), "ensemble_phase") !=
        batch.views.end(), "batch does not advertise an aggregate phase view");
    const auto payload = dsl::buildHierarchyPayload(
        batch, "ensemble_phase", root, archiveRoot);
    require(payload.figures.size() == 1 && !payload.figures.front().square &&
        payload.figures.front().traces.size() == 1,
        "native hierarchy aggregate phase payload is incomplete");
    const auto observables = dsl::buildHierarchyPayload(batch, "observable_time", root, archiveRoot);
    require(!observables.figures.empty(), "hierarchy observable workspace is empty");
    for (const auto& figure : observables.figures)
        require(figure.traces.size() == payload.figures.front().traces.size()
            && figure.traces.front().label == payload.figures.front().traces.front().label
            && figure.traces.front().color == payload.figures.front().traces.front().color,
            "hierarchy observable traces differ from Main trajectories");
    const int systemIndex = catalog.indexOf("system");
    require(systemIndex >= 0, "completed system run is absent from archive tree");
    const auto& system = catalog.items[static_cast<std::size_t>(systemIndex)];
    require(system.kind == "system_run" && system.parentId == batch.id &&
        system.properties.count("parameter · damping") == 1 &&
        system.properties.count("group · damping_ratio") == 1,
        "system run does not expose its exact parameters beneath the batch");
    const int clusterIndex = catalog.indexOf("cluster");
    require(clusterIndex >= 0 &&
        catalog.items[static_cast<std::size_t>(clusterIndex)].parentId == system.id,
        "cluster is not nested beneath its system run");
    require(catalog.items[static_cast<std::size_t>(index)].parentId == "cluster",
        "trajectory is not nested beneath its cluster");
    const auto& trajectory = catalog.items[static_cast<std::size_t>(index)];
    require(trajectory.properties.count("IC") == 1 &&
            trajectory.properties.count("initial state") == 0 &&
            trajectory.bytes > 0.0 && trajectory.modified != 0.0,
        "trajectory does not expose IC naming or sortable archive metrics");
    require(std::find(trajectory.views.begin(), trajectory.views.end(),
        "time_frequency") != trajectory.views.end(),
        "ODE trajectory does not advertise combined time/frequency analysis");
    require(std::none_of(catalog.items.begin(), catalog.items.end(), [](const auto& item) {
        return item.kind == "system_configuration" || item.kind == "forcing_campaign";
    }), "parameter sets or campaigns were inserted as solve-tree hierarchy levels");
    require(std::find(batch.tags.begin(), batch.tags.end(), "free") !=
        batch.tags.end() && std::find(batch.tags.begin(), batch.tags.end(), "forced") ==
        batch.tags.end(), "free campaign modifier was changed to forced");
    require(std::find(batch.views.begin(), batch.views.end(), "campaign_summary") !=
        batch.views.end(), "campaign-tagged batch does not advertise its response view");
    const auto summary = dsl::buildHierarchyPayload(
        batch, "campaign_summary", root, archiveRoot);
    require(summary.figures.size() == 2,
        "native campaign response payload is incomplete");
    const auto bifurcation = dsl::buildHierarchyPayload(
        batch, "campaign_bifurcation", root, archiveRoot);
    require(bifurcation.figures.size() == 1 &&
        !bifurcation.figures.front().traces.front().x.empty(),
        "native campaign bifurcation response is incomplete");
    std::filesystem::remove_all(root, ignored);
}

void legacyHierarchyReference() {
    const auto root = std::filesystem::temp_directory_path() /
        "dynsyslab-native-legacy-hierarchy";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    const auto archiveRoot = root / "archives";
    const auto solution = archiveRoot / "solutions" / "duffing" / "legacy.npz";
    const auto manifest = archiveRoot / "hierarchies" / "legacy.npz";
    writeSolution(solution);
    dynsyslab::core::NpzRecord record;
    const std::string metadata =
        R"({"format":"dynsyslab.hierarchy/v1","job":{"job_id":"legacy-job","label":"legacy batch","tags":{}},"systems":[{"label":"legacy Duffing","dynamical_system":{"name":"duffing","preset":{"description":"reference"},"driver":null},"clusters":[{"label":"legacy cluster","trajectories":[{"label":"legacy trajectory","initial_condition":[1,0],"metadata":{"status":"complete","solver_metadata":{"integration_form":"dimensionless"}},"solution_archive":")" +
        solution.generic_u8string() + R"("}]}]}]})";
    record.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(metadata);
    dynsyslab::core::writeNpzAtomic(manifest, record);

    const auto catalog = dsl::Catalog::scan(root, archiveRoot);
    const int batchIndex = catalog.indexOf("legacy-batch-legacy-job");
    require(batchIndex >= 0, "legacy completed batch is absent from archive tree");
    const auto& batch = catalog.items[static_cast<std::size_t>(batchIndex)];
    require(batch.category == "Archived solves" && batch.parentId.empty(),
        "legacy batch is not an archive-tree root");
    const int systemIndex = catalog.indexOf("legacy-batch-legacy-job-system-0000");
    require(systemIndex >= 0 &&
        catalog.items[static_cast<std::size_t>(systemIndex)].parentId == batch.id,
        "legacy system is not nested beneath its batch");
    const std::string clusterId = "legacy-batch-legacy-job-system-0000-cluster-0000";
    const int clusterIndex = catalog.indexOf(clusterId);
    require(clusterIndex >= 0 && catalog.items[static_cast<std::size_t>(clusterIndex)]
        .parentId == catalog.items[static_cast<std::size_t>(systemIndex)].id,
        "legacy cluster is not nested beneath its system");
    const int trajectoryIndex = catalog.indexOf(clusterId + "-trajectory-0000");
    require(trajectoryIndex >= 0 && catalog.items[static_cast<std::size_t>(trajectoryIndex)]
        .parentId == clusterId, "legacy trajectory is not nested beneath its cluster");
    require(std::none_of(catalog.items.begin(), catalog.items.end(), [&](const auto& item) {
        return item.kind == "trajectory" &&
            std::filesystem::u8path(item.source) == manifest;
    }), "legacy hierarchy container leaked into the tree as a trajectory");
    const auto update = dynsyslab::core::reconcileArchiveIndex(archiveRoot);
    require(update.indexed == 2, "legacy hierarchy index fixture was not published");
    const auto indexed = dsl::Catalog::indexed(root, archiveRoot);
    require(indexed.indexOf("legacy-batch-legacy-job") >= 0 &&
            indexed.indexOf(clusterId + "-trajectory-0000") >= 0,
            "indexed legacy hierarchy summary does not reproduce the archive tree");
    const auto snapshot = dynsyslab::core::loadArchiveCatalogSnapshot(
        archiveRoot, "dynsyslab.viewer.native", 11);
    require(snapshot.has_value(), "indexed catalog was not retained for startup reuse");
    const auto cached = dsl::Catalog::indexed(root, archiveRoot);
    require(cached.items.size() == indexed.items.size() &&
            cached.indexOf(clusterId + "-trajectory-0000") >= 0,
            "retained catalog snapshot does not reproduce indexed contents");
    const auto payload = dsl::buildHierarchyPayload(
        batch, "ensemble_phase", root, archiveRoot);
    require(payload.figures.size() == 1 && payload.figures.front().traces.size() == 1,
        "legacy hierarchy aggregate phase payload is incomplete");
    std::filesystem::remove_all(root, ignored);
}

}  // namespace

int main() {
    try {
        require(dsl::canonicalSystemLabel("rayleigh-plesset") == "Rayleigh-Plesset" &&
                dsl::canonicalSystemLabel("Rayleigh-Plesset") == "Rayleigh-Plesset" &&
                dsl::canonicalSystemLabel("RPE") == "Rayleigh-Plesset",
                "Rayleigh-Plesset aliases do not share one catalog identity");
        portableHierarchyReference(); legacyHierarchyReference();
    }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
