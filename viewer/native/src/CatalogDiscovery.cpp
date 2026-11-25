#include "dynsyslab_viewer/Model.hpp"
#include "ArchiveCatalog.hpp"
#include "HierarchyCatalog.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <set>
#include <sstream>

namespace dsl {
namespace {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char item) {
        return static_cast<char>(std::tolower(item));
    });
    return value;
}

std::string systemToken(std::string value) {
    value = lower(std::move(value));
    std::string result;
    for (const unsigned char item : value)
        if (std::isalnum(item)) result.push_back(static_cast<char>(item));
    return result;
}

void canonicalize(Catalog& catalog) {
    for (auto& item : catalog.items) {
        const auto previous = item.system;
        item.system = canonicalSystemLabel(previous);
        for (auto& tag : item.tags)
            if (systemToken(tag) == systemToken(previous)) tag = item.system;
        if (item.kind == "system_run" &&
            systemToken(item.label) == systemToken(previous)) item.label = item.system;
    }
}

Catalog baseCatalog(const std::filesystem::path& projectRoot) {
    Catalog result;
    result.projectName = "DynSysLab";
    result.projectRoot = projectRoot.generic_u8string();
    result.viewLabels = {{"open", "Open"}, {"data_overview", "Data overview"},
        {"koopman_spectrum", "Spectrum"}, {"mode_weights", "Mode weights"}, {"time", "Time series"},
        {"krom_dashboard", "KROM dashboard"}, {"phase", "Phase space"},
        {"field", "Space-time field"}, {"three_d", "Three-state view"},
        {"fourier", "Frequency spectra"},
        {"control_fourier", "Control frequency spectra"},
        {"estimator", "Estimator measurement"},
        {"time_frequency", "Time + frequency"},
        {"ensemble_time", "Ensemble time series"},
        {"ensemble_phase", "Ensemble phase space"},
        {"ensemble_control", "Ensemble controls"},
        {"campaign_summary", "Campaign response"},
        {"campaign_bifurcation", "Bifurcation response"},
        {"lattice", "Lattice states"}, {"control", "Control"},
        {"observable_time", "Observable time series"},
        {"observable_heatmap", "Observable phase maps"},
        {"domain_field", "Function on its domain"},
        {"observable_phase", "Observable on the phase plane"},
        {"animation", "Animation"},
        {"control_map_error", "Control-map error"},
        {"observability", "Observability"},
        {"hankel_spectrum", "Hankel singular values"},
        {"gramian_spectrum", "Gramians"},
        {"frequency_response", "Frequency response"},
        {"pole_zero", "Poles and zeros"},
        {"bode", "Bode plot"},
        {"nyquist", "Nyquist plot"},
        {"loop_bode", "Loop Bode plot"},
        {"sensitivity", "Sensitivity"},
        {"koopman_validation", "KROM validation"},
        {"eigenvalue_optimization", "Eigenvalue optimization"},
        {"krom_control_campaign", "KROM control campaign"}};
    return result;
}

bool operationalOrDraft(const std::filesystem::path& path,
                        const std::filesystem::path& archiveRoot) {
    const auto relative = path.lexically_relative(archiveRoot);
    if (!relative.empty() && *relative.begin() == L"job_manager") return true;
    if (lower(path.extension().u8string()) != ".json") return false;
    try {
        return json::load(path.wstring()).at("schema").string()
            == "dynsyslab.compute-job/v1";
    } catch (...) { return false; }
}

std::string stableIdentifier(const std::filesystem::path& root,
                             const std::filesystem::path& path) {
    const auto absoluteRoot = std::filesystem::absolute(root).lexically_normal();
    const auto absolutePath = std::filesystem::absolute(path).lexically_normal();
    const auto relative = absolutePath.lexically_relative(absoluteRoot);
    const std::string key = lower(relative.empty() || *relative.begin() == L".."
        ? absolutePath.generic_u8string() : relative.generic_u8string());
    std::uint64_t hash = 1469598103934665603ull;
    for (const unsigned char value : key) { hash ^= value; hash *= 1099511628211ull; }
    std::ostringstream result;
    result << std::hex << std::setw(16) << std::setfill('0') << hash;
    return result.str();
}

std::string systemFromPath(const std::filesystem::path& path) {
    const std::string value = lower(path.generic_u8string());
    const std::pair<const char*, const char*> systems[] = {
        {"reddy", "Reddy-Szeri"}, {"marmottant", "Marmottant"},
        {"duffing", "Duffing"}, {"rpe", "Rayleigh-Plesset"},
        {"rayleigh", "Rayleigh-Plesset"}, {"lorenz", "Lorenz"},
        {"rossler", "Rossler"}, {"van_der_pol", "Van der Pol"}};
    for (const auto& entry : systems)
        if (value.find(entry.first) != std::string::npos) return entry.second;
    return {};
}

CatalogItem scannedItem(const std::filesystem::path& root,
                        const std::filesystem::path& path,
                        bool readFileMetadata = true) {
    CatalogItem item;
    item.id = stableIdentifier(root, path); item.source = path.generic_u8string();
    item.label = path.stem().u8string(); item.system = systemFromPath(path);
    if (readFileMetadata) {
        std::error_code error;
        item.bytes = static_cast<double>(std::filesystem::file_size(path, error));
        if (error) item.bytes = 0.0;
        error.clear(); const auto modified = std::filesystem::last_write_time(path, error);
        item.modified = error ? 0.0 : static_cast<double>(
            modified.time_since_epoch().count());
    }
    const auto extension = lower(path.extension().u8string());
    const auto key = lower(item.label);
    item.kind = "trajectory"; item.category = "Training data";
    item.subtitle = "Native archive"; item.views = {"data_overview"};
    if (extension == ".html" || extension == ".pdf") {
        item.kind = "report"; item.category = "Reports";
        item.subtitle = extension == ".pdf" ? "PDF Viewer artifact" : "HTML Viewer artifact";
        item.views = {"open"};
    } else if (extension == ".json") {
        item.kind = "operation_record"; item.category = "Operations";
        item.subtitle = "JSON job, campaign, or provenance record"; item.views = {"open"};
    } else if (key.find("observability") != std::string::npos) {
        item.kind = "observability"; item.category = "Estimation";
        item.views = {"observability"};
    } else if (key.find("pca") != std::string::npos ||
               key.find("pod") != std::string::npos ||
               key.find("svd") != std::string::npos) {
        item.kind = "decomposition"; item.category = "Post-processing";
        item.subtitle = "Native matrix decomposition";
    } else if (key.find("krom") != std::string::npos ||
               key.find("sindy") != std::string::npos ||
               key.find("hdmd") != std::string::npos ||
               key.find("koopman") != std::string::npos) {
        item.kind = "krom"; item.category = "KROMs";
        item.views = {"koopman_spectrum", "data_overview"};
        if (key.find("control-regression") != std::string::npos)
            item.views.push_back("control_map_error");
    } else if (key.find("control") != std::string::npos) {
        item.kind = "control_data"; item.category = "Control data";
    }
    if (!extension.empty()) item.tags.push_back(extension.substr(1));
    if (!item.system.empty()) item.tags.push_back(item.system);
    item.properties.emplace("format", extension.empty() ? "file" : extension.substr(1));
    item.properties.emplace("native catalog", "filesystem metadata");
    item.lineageId = item.id;
    return item;
}

void finish(Catalog& catalog, const std::filesystem::path& projectRoot) {
    resolveArchiveCatalogRelations(catalog, projectRoot);
    std::sort(catalog.items.begin(), catalog.items.end(),
        [](const CatalogItem& left, const CatalogItem& right) {
            return left.modified != right.modified
                ? left.modified > right.modified : left.label < right.label;
        });
    canonicalize(catalog);
}

}  // namespace

Catalog Catalog::scan(const std::filesystem::path& projectRoot,
                      const std::filesystem::path& archiveRoot) {
    auto result = baseCatalog(projectRoot);
    std::set<std::string> hierarchySources;
    appendHierarchyCatalog(result, projectRoot, archiveRoot, hierarchySources);
    auto scan = [&](const std::filesystem::path& base, bool recursive) {
        if (!std::filesystem::is_directory(base)) return;
        auto accept = [&](const std::filesystem::directory_entry& entry) {
            if (!entry.is_regular_file() ||
                (base == archiveRoot && operationalOrDraft(entry.path(), archiveRoot)) ||
                hierarchySources.count(normalizedSource(entry.path()))) return;
            if (lower(entry.path().extension().u8string()) == ".npz")
                result.items.push_back(scannedItem(projectRoot, entry.path()));
        };
        if (recursive)
            for (const auto& entry : std::filesystem::recursive_directory_iterator(
                base, std::filesystem::directory_options::skip_permission_denied)) accept(entry);
        else for (const auto& entry : std::filesystem::directory_iterator(base)) accept(entry);
    };
    scan(archiveRoot, true); scan(projectRoot / L"viewer_outputs", false);
    for (auto& item : result.items) {
        auto source = std::filesystem::u8path(item.source);
        if (source.is_relative()) source = projectRoot / source;
        if (source.extension() == L".npz" && std::filesystem::is_regular_file(source))
            enrichArchiveCatalogItem(item, source);
    }
    finish(result, projectRoot);
    return result;
}

Catalog Catalog::indexed(const std::filesystem::path& projectRoot,
                         const std::filesystem::path& archiveRoot) {
    constexpr std::uint32_t snapshotVersion = 11;  // 11: index summary v5; 10: observable series; 9: animation, domain fields
    constexpr const char* consumer = "dynsyslab.viewer.native";
    if (!std::filesystem::is_regular_file(
            dynsyslab::core::archiveIndexPath(archiveRoot)))
        (void)dynsyslab::core::reconcileArchiveIndex(archiveRoot);
    if (const auto snapshot = dynsyslab::core::loadArchiveCatalogSnapshot(
            archiveRoot, consumer, snapshotVersion)) {
        try {
            auto cached = Catalog::fromJson(*snapshot);
            if (std::filesystem::absolute(std::filesystem::u8path(cached.projectRoot))
                    .lexically_normal() ==
                std::filesystem::absolute(projectRoot).lexically_normal()) return cached;
        } catch (...) {}
    }
    const auto generation = dynsyslab::core::archiveIndexGeneration(archiveRoot);
    const auto entries = dynsyslab::core::loadArchiveIndex(archiveRoot);
    if (entries.empty()) return {};
    auto result = baseCatalog(projectRoot);
    std::set<std::string> hierarchySources;
    appendHierarchyCatalog(result, projectRoot, archiveRoot, hierarchySources, &entries);
    for (const auto& entry : entries) {
        if (hierarchySources.count(normalizedSource(entry.path))) continue;
        auto item = scannedItem(projectRoot, entry.path, false);
        item.bytes = static_cast<double>(entry.sizeBytes);
        item.modified = static_cast<double>(entry.modifiedTicks);
        enrichArchiveCatalogItem(item, {entry.arrays, entry.summaryJson});
        result.items.push_back(std::move(item));
    }
    finish(result, projectRoot);
    try {
        (void)dynsyslab::core::storeArchiveCatalogSnapshot(
            archiveRoot, consumer, snapshotVersion, generation, result.toJson());
    } catch (...) {}
    return result;
}

}  // namespace dsl
