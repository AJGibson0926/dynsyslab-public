#include "../src/ModeWeightPayload.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/ComparisonPayload.hpp"
#include "dynsyslab_viewer/RecordWindow.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void latticeView() {
    const auto path=std::filesystem::temp_directory_path()/"dynsyslab-native-lattice-view-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["lattice_initial"]=dynsyslab::core::NpyArray::fromInt64s({2,2},{1,-1,-1,1});
    archive.arrays["lattice_final"]=dynsyslab::core::NpyArray::fromInt64s({2,2},{1,1,1,1});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"Ising lattice","system":"ising","system_family":"lattice"})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    const auto payload=dsl::buildNativePayload(path,"lattice");
    require(payload.figures.size()==2 && payload.figures[0].square
        && payload.figures[1].traces[0].z.back()==1.0,"native lattice payload");
    std::error_code ignored;std::filesystem::remove(path,ignored);
}

void decompositionViews() {
    const auto root=std::filesystem::temp_directory_path();
    const auto pcaPath=root/"dynsyslab-native-pca-payload-test.npz";
    dynsyslab::core::NpzRecord pca;
    pca.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({3},{0.0,.1,.2});
    pca.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles({2},{4.0,1.0});
    pca.arrays["explained_variance_ratio"]=
        dynsyslab::core::NpyArray::fromDoubles({2},{.8,.2});
    pca.arrays["components"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{1.0,0.0,0.0,1.0});
    pca.arrays["scores"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{0.0,1.0,1.0,0.0,2.0,-1.0});
    pca.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"native PCA"})");
    dynsyslab::core::writeNpzAtomic(pcaPath,pca);
    const auto pcaPayload=dsl::buildNativePayload(pcaPath,"data_overview");
    require(pcaPayload.figures.size()==4,"PCA payload figure count");
    require(pcaPayload.figures[2].traces[0].heatmap,
            "PCA component heatmap");
    const auto robustPath=root/"dynsyslab-native-robust-pca-payload-test.npz";
    dynsyslab::core::NpzRecord robust;
    robust.arrays["time"]=pca.arrays["time"];
    robust.arrays["low_rank"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{1.0,2.0,2.0,4.0,3.0,6.0});
    robust.arrays["sparse"]=dynsyslab::core::NpyArray::fromDoubles(
        {3,2},{0.0,0.0,0.0,5.0,0.0,0.0});
    robust.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"native robust PCA"})");
    dynsyslab::core::writeNpzAtomic(robustPath,robust);
    const auto robustPayload=dsl::buildNativePayload(robustPath,"data_overview");
    require(robustPayload.figures.size()==4,"robust PCA payload figure count");
    require(robustPayload.figures[1].traces[0].role=="sparse",
            "robust PCA sparse trace role");
    std::error_code ignored; std::filesystem::remove(pcaPath,ignored);
    std::filesystem::remove(robustPath,ignored);
}

void freeCampaignCatalog() {
    const auto root = std::filesystem::temp_directory_path()
        / "dynsyslab-native-free-campaign-catalog";
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    const auto archiveRoot = root / "archives";
    const auto solution = archiveRoot / "solutions" / "ising" / "sample.npz";
    std::filesystem::create_directories(solution.parent_path());
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles({1}, {0.0});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles({1, 2}, {0.0, -2.0});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString("{}");
    dynsyslab::core::writeNpzAtomic(solution, archive);
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier = "ising";
    configuration.parameters = {{"temperature", 2.0}};
    dynsyslab::core::BatchRecord batch;
    batch.id = "batch"; batch.label = "temperature sweep";
    batch.createdUtc = "now"; batch.dataRole = "campaign";
    batch.campaignId = "campaign"; batch.campaignKind = "temperature_sweep";
    batch.tags = {"free", "campaign"};
    dynsyslab::core::SystemRunRecord system;
    system.id = "system"; system.systemIdentifier = "ising";
    system.systemLabel = "2-D Ising model"; system.parameterSetKey =
        dynsyslab::core::parameterSetKey(configuration);
    system.coordinateSystem = "discrete"; system.stateNames = {"m", "energy"};
    dynsyslab::core::ClusterRecord cluster; cluster.id = "cluster";
    dynsyslab::core::TrajectoryRecord trajectory;
    trajectory.id = "trajectory"; trajectory.label = "sample";
    trajectory.solutionArchive = solution; trajectory.tags = {"free", "campaign"};
    cluster.trajectories.push_back(trajectory);
    system.clusters.push_back(cluster); batch.systems.push_back(system);
    dynsyslab::core::writeHierarchyAtomic(
        archiveRoot / "hierarchies" / "ising" / "campaign.hierarchy.json", batch);
    const auto catalog = dsl::Catalog::scan(root, archiveRoot);
    const int index = catalog.indexOf("trajectory");
    require(index >= 0, "free campaign trajectory catalog entry");
    const auto& item = catalog.items[static_cast<std::size_t>(index)];
    require(std::find(item.tags.begin(), item.tags.end(), "free") != item.tags.end(),
            "free campaign retains free modifier");
    require(std::find(item.tags.begin(), item.tags.end(), "forced") == item.tags.end(),
            "free campaign is not mislabeled forced");
    require(std::find(item.views.begin(), item.views.end(), "control") == item.views.end(),
            "free campaign does not advertise control history");
    std::filesystem::remove_all(root, ignored);
}

void metadataDrivenCatalogRelations() {
    const auto root=std::filesystem::temp_directory_path()/
        "dynsyslab-native-metadata-catalog";
    std::error_code ignored; std::filesystem::remove_all(root,ignored);
    const auto archiveRoot=root/"archives";
    const auto source=archiveRoot/"solutions"/"duffing"/"source.npz";
    const auto derived=archiveRoot/"derived"/"koopman-krom-extract"/"opaque-name.npz";
    const auto spectrum=archiveRoot/"derived"/"koopman-hankel-dmd"/"spectrum.npz";
    const auto expansion=archiveRoot/"solutions"/"duffing"/"domain-expansion.npz";
    std::filesystem::create_directories(source.parent_path());
    std::filesystem::create_directories(derived.parent_path());
    std::filesystem::create_directories(spectrum.parent_path());
    dynsyslab::core::NpzRecord trajectory;
    trajectory.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({2},{0.0,1.0});
    trajectory.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.0,1.0,1.0});
    trajectory.arrays["control_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,1},{0.0,0.25});
    trajectory.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"schema":"dynsyslab.solution/v1","label":"Duffing source","system":"duffing","system_label":"Duffing oscillator"})");
    dynsyslab::core::writeNpzAtomic(source,trajectory);
    dynsyslab::core::NpzRecord krom;
    krom.arrays["continuous_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {1},{{-0.1,1.0}});
    krom.arrays["discrete_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {1},{{0.9,0.1}});
    krom.arrays["coefficients"]=dynsyslab::core::NpyArray::fromComplex(
        {1,1},{{1.0,0.0}});
    krom.arrays["polynomial_exponents"]=dynsyslab::core::NpyArray::fromDoubles(
        {1,2},{1.0,0.0});
    krom.arrays["state_center"]=dynsyslab::core::NpyArray::fromDoubles({2},{0.0,0.0});
    krom.arrays["state_scale"]=dynsyslab::core::NpyArray::fromDoubles({2},{1.0,1.0});
    krom.arrays["training_states"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.0,1.0,1.0});
    krom.arrays["control_error_axis0"]=dynsyslab::core::NpyArray::fromDoubles(
        {2},{-1.0,1.0});
    krom.arrays["control_error_axis1"]=dynsyslab::core::NpyArray::fromDoubles(
        {2},{-1.0,1.0});
    krom.arrays["constant_relative_error"]=dynsyslab::core::NpyArray::fromDoubles(
        {2,2},{0.0,0.1,0.2,0.3});
    const std::string metadata="{\"schema\":\"dynsyslab.derived/v1\","
        "\"label\":\"Independent minimal KROM\",\"method\":\"krom_extract\","
        "\"system\":\"duffing\",\"system_label\":\"Duffing oscillator\","
        "\"source_archive\":\""+source.generic_u8string()+"\","
        "\"selected_modes\":[7],\"control_maps\":[\"exact\",\"bilinear\"],"
        "\"independent_archive\":true,\"schema\":\"dynsyslab.krom/v1\","
        "\"krom\":{\"parent\":null},\"eigenfunction_representation\":\"polynomial\","
        "\"data_object\":{\"schema\":\"dynsyslab.data-object/v1\",\"kind\":\"krom\"}}";
    krom.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(metadata);
    dynsyslab::core::writeNpzAtomic(derived,krom);
    dynsyslab::core::NpzRecord expanded=trajectory;
    expanded.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        "{\"schema\":\"dynsyslab.solution/v1\",\"label\":\"Domain expansion\","
        "\"system\":\"duffing\",\"training_preflight\":{\"enabled\":true,"
        "\"source_krom\":\""+derived.generic_u8string()+"\"}}");
    dynsyslab::core::writeNpzAtomic(expansion,expanded);
    dynsyslab::core::NpzRecord hdmd;
    hdmd.arrays["continuous_eigenvalues"]=krom.arrays["continuous_eigenvalues"];
    hdmd.arrays["eigenfunction_values"]=dynsyslab::core::NpyArray::fromComplex(
        {2,1},{{1.0,0.0},{0.9,0.1}});
    hdmd.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"schema":"dynsyslab.derived/v1","label":"H-DMD spectrum","method":"hdmd","system":"duffing"})");
    dynsyslab::core::writeNpzAtomic(spectrum,hdmd);
    const auto catalog=dsl::Catalog::scan(root,archiveRoot);
    const dsl::CatalogItem* sourceItem=nullptr; const dsl::CatalogItem* kromItem=nullptr;
    const dsl::CatalogItem* spectrumItem=nullptr;
    const dsl::CatalogItem* expansionItem=nullptr;
    for (const auto& item : catalog.items) {
        if (item.label=="Duffing source") sourceItem=&item;
        if (item.label=="Independent minimal KROM") kromItem=&item;
        if (item.label=="H-DMD spectrum") spectrumItem=&item;
        if (item.label=="Domain expansion") expansionItem=&item;
    }
    require(sourceItem && kromItem,"metadata-driven catalog objects are missing");
    require(std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "observable_time")!=sourceItem->views.end()
        && std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "observable_heatmap")!=sourceItem->views.end(),
        "two-state trajectory did not advertise native observable views");
    require(std::find(sourceItem->views.begin(),sourceItem->views.end(),
        "control_fourier")!=sourceItem->views.end(),
        "controlled trajectory did not advertise control spectra");
    require(spectrumItem && spectrumItem->kind=="koopman_spectrum"
        && spectrumItem->category=="Koopman spectra",
        "H-DMD spectrum was incorrectly exposed as a control-ready KROM");
    require(kromItem->kind=="krom" && kromItem->category=="KROMs",
        "metadata-driven KROM classification");
    require(!kromItem->views.empty() && kromItem->views.front()=="krom_dashboard",
        "renderable KROM does not default to the standard dashboard");
    require(std::find(kromItem->views.begin(),kromItem->views.end(),
        "control_map_error")!=kromItem->views.end(),"KROM error-map view discovery");
    require(kromItem->properties.at("source mode indices")=="7"
        && kromItem->properties.at("parent required at runtime")=="no",
        "independent KROM properties");
    require(std::any_of(kromItem->relations.begin(),kromItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==sourceItem->id
            && relation.kind=="extracted from";}),
        "archive lineage relation was not resolved");
    require(std::any_of(sourceItem->relations.begin(),sourceItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==kromItem->id
            && relation.kind=="has extracted KROM";}),
        "archive reverse lineage relation was not generated");
    require(expansionItem && std::any_of(expansionItem->relations.begin(),
        expansionItem->relations.end(),[&](const dsl::Relation& relation){
            return relation.target==kromItem->id
                && relation.kind=="expands training domain of";}),
        "training-domain source KROM relation was not generated");
    require(std::any_of(kromItem->relations.begin(),kromItem->relations.end(),
        [&](const dsl::Relation& relation){return relation.target==expansionItem->id
            && relation.kind=="has domain-expansion data";}),
        "training-domain reverse relation was not generated");
    const auto update = dynsyslab::core::reconcileArchiveIndex(archiveRoot);
    const auto indexedCatalog = dsl::Catalog::indexed(root, archiveRoot);
    require(update.indexed >= 4 && indexedCatalog.indexOf(kromItem->id) >= 0 &&
            indexedCatalog.indexOf(sourceItem->id) >= 0,
            "central archive index did not reproduce catalog metadata and lineage");
    std::filesystem::remove_all(root,ignored);
}

void pdeFieldView() {
    const auto path = std::filesystem::temp_directory_path()
        / "dynsyslab-native-pde-field-test.npz";
    dynsyslab::core::NpzRecord archive;
    archive.arrays["time"] = dynsyslab::core::NpyArray::fromDoubles(
        {3}, {0.0, 0.1, 0.2});
    archive.arrays["space"] = dynsyslab::core::NpyArray::fromDoubles(
        {4}, {0.0, 0.25, 0.5, 0.75});
    archive.arrays["values"] = dynsyslab::core::NpyArray::fromDoubles(
        {3, 4}, {0, 1, 0, -1, 0, 0.8, 0, -0.8, 0, 0.6, 0, -0.6});
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        R"({"label":"heat field","system_family":"pde"})");
    dynsyslab::core::writeNpzAtomic(path, archive);
    const auto payload = dsl::buildNativePayload(path, "field");
    require(payload.figures.size() == 2, "PDE field and profile figures");
    require(payload.figures[0].traces[0].heatmap,
            "PDE space-time heatmap");
    require(payload.figures[0].traces[0].zShape == std::vector<int>({3, 4}),
            "PDE heatmap shape");
    require(payload.figures[1].traces[0].y.back() == -0.6,
            "PDE final profile values");
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
}

}  // namespace

int main() {
    try {
        latticeView(); decompositionViews(); freeCampaignCatalog(); metadataDrivenCatalogRelations(); pdeFieldView();
        std::cout << "dynsyslab_native_payload_catalog_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_native_payload_catalog_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
