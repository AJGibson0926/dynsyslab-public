// Storage-location registry, verified copies and campaign relocation (core/StorageLocations.hpp).
#include <dynsyslab/core/StorageLocations.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeGrid.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::string utc(const char* format) {
    const std::time_t raw = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm value{};
#ifdef _WIN32
    gmtime_s(&value, &raw);
#else
    gmtime_r(&raw, &value);
#endif
    std::ostringstream out;
    out << std::put_time(&value, format);
    return out.str();
}

void writeJsonAtomic(const std::filesystem::path& path, const json::Value& value) {
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream || !(stream << json::dump(value) << '\n') || !stream.flush())
            throw std::runtime_error("unable to write " + path.u8string());
    }
    std::filesystem::rename(temporary, path);
}

std::string generic(const std::filesystem::path& path) { return path.generic_u8string(); }

std::filesystem::path normalized(const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
}

// `path` relative to `root` when it lies inside it; empty otherwise.
std::filesystem::path inside(const std::filesystem::path& path, const std::filesystem::path& root) {
    const auto relative = normalized(path).lexically_relative(normalized(root));
    if (relative.empty() || *relative.begin() == "..") return {};
    return relative;
}

bool isFile(const std::filesystem::path& path) {
    std::error_code ignored;
    return std::filesystem::is_regular_file(path, ignored);
}

json::Value locationValue(const StorageLocation& location) {
    return json::Value(json::Value::Object{{"name", json::Value(location.name)},
        {"root", json::Value(generic(location.root))}, {"writable", json::Value(location.writable)},
        {"purpose", json::Value(location.purpose)}});
}

json::Value copyValue(const VerifiedCopy& copy) {
    json::Value::Object value{{"source", json::Value(generic(copy.source))},
        {"destination", json::Value(generic(copy.destination))}, {"sha256", json::Value(copy.sha256)},
        {"bytes", json::Value(static_cast<double>(copy.bytes))}, {"verified", json::Value(true)},
        {"source_removed", json::Value(copy.sourceRemoved)},
        {"already_at_destination", json::Value(copy.alreadyAtDestination)},
        {"migrated_utc", json::Value(copy.copiedUtc)}};
    if (!copy.relativePath.empty()) value["relative_archive_path"] = json::Value(copy.relativePath);
    return json::Value(std::move(value));
}

json::Value readRegistry(const std::filesystem::path& archiveRoot) {
    const auto path = storageRegistryPath(archiveRoot);
    if (!isFile(path)) return json::Value(json::Value::Object{{"schema", json::Value(kStorageRegistrySchema)}});
    auto registry = json::load(path.wstring());
    if (registry.at("schema").string() != kStorageRegistrySchema)
        throw std::runtime_error(path.u8string() + " is not a " + kStorageRegistrySchema + " registry");
    return registry;
}

}  // namespace

std::filesystem::path storageRegistryPath(const std::filesystem::path& archiveRoot) {
    return archiveRoot / L"management" / L"storage_locations.json";
}

std::vector<StorageLocation> storageLocations(const std::filesystem::path& archiveRoot) {
    std::vector<StorageLocation> result{{"primary", normalized(archiveRoot), true, "solutions"}};
    const auto registry = readRegistry(archiveRoot);
    for (const auto& [name, entry] : registry.at("locations").object()) {
        if (name == "primary") continue;
        result.push_back({name, std::filesystem::u8path(entry.at("root").string()),
                          entry.at("writable").boolean(true), entry.at("purpose").string("solutions")});
    }
    return result;
}

void registerStorageLocation(const std::filesystem::path& archiveRoot, const StorageLocation& location) {
    if (location.name.empty() || location.name == "primary")
        throw std::invalid_argument("a storage location needs a name other than 'primary'");
    if (location.root.empty()) throw std::invalid_argument("a storage location needs a root folder");
    StorageLocation stored = location;
    stored.root = normalized(location.root);
    std::filesystem::create_directories(stored.root);
    auto registry = readRegistry(archiveRoot);
    registry.at("locations").at("primary") = locationValue({"primary", normalized(archiveRoot), true, "solutions"});
    registry.at("locations").at(stored.name) = locationValue(stored);
    registry.at("updated_utc") = json::Value(utc("%Y-%m-%dT%H:%M:%SZ"));
    writeJsonAtomic(storageRegistryPath(archiveRoot), registry);
}

VerifiedCopy copyVerified(const std::filesystem::path& sourceArgument, const std::filesystem::path& folder,
                          bool removeSource) {
    const auto source = normalized(sourceArgument);
    if (!isFile(source)) throw std::invalid_argument("copy source is not a file: " + source.u8string());
    if (source.filename().u8string().find(".partial") != std::string::npos)
        throw std::invalid_argument("active partial archives may not be copied: " + source.u8string());
    const auto target = normalized(folder) / source.filename();
    if (target == source) throw std::invalid_argument("copy destination is the source itself");
    VerifiedCopy copy;
    copy.source = source; copy.destination = target;
    copy.sha256 = sha256File(source); copy.bytes = std::filesystem::file_size(source);
    if (isFile(target)) {
        if (std::filesystem::file_size(target) != copy.bytes || sha256File(target) != copy.sha256)
            throw std::runtime_error("a different file already exists at " + target.u8string());
    } else {
        std::filesystem::create_directories(target.parent_path());
        const auto temporary = target.parent_path() / (target.filename().wstring() + L".copy.partial");
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        try {
            std::filesystem::copy_file(source, temporary);
            if (std::filesystem::file_size(temporary) != copy.bytes || sha256File(temporary) != copy.sha256)
                throw std::runtime_error("copy verification failed; the source is kept");
            std::filesystem::rename(temporary, target);
        } catch (...) {
            std::filesystem::remove(temporary, ignored);
            throw;
        }
    }
    if (removeSource) { std::filesystem::remove(source); copy.sourceRemoved = true; }
    copy.copiedUtc = utc("%Y-%m-%dT%H:%M:%SZ");
    return copy;
}

std::vector<std::filesystem::path> campaignArchives(const std::filesystem::path& archiveRoot,
                                                    const std::string& campaignId) {
    if (campaignId.empty()) throw std::invalid_argument("a campaign id is required");
    std::vector<std::filesystem::path> result;
    const auto add = [&](const std::filesystem::path& path) { if (isFile(path)) result.push_back(normalized(path)); };
    std::error_code ignored;
    if (std::filesystem::is_directory(archiveRoot / L"campaigns", ignored))
        for (const auto& entry : std::filesystem::directory_iterator(archiveRoot / L"campaigns")) {
            if (!entry.is_regular_file() || entry.path().extension() != ".json") continue;
            try {
                const auto manifest = json::load(entry.path().wstring());
                if (manifest.at("id").string(manifest.at("campaign_id").string()) != campaignId) continue;
                add(entry.path());
                for (const auto& point : manifest.at("points").array()) {
                    add(std::filesystem::u8path(point.at("hierarchy_manifest").string()));
                    for (const auto& archive : point.at("solution_archives").array())
                        add(std::filesystem::u8path(archive.string()));
                }
            } catch (const std::exception&) {
                continue;  // not a campaign manifest
            }
        }
    ArchiveAuditScope scope;
    scope.modifiedWithinDays = -1.0;
    scope.sections = {"solutions"};
    scope.campaignId = campaignId;
    scope.includeDerived = true;
    for (const auto& path : scopedArchives(archiveRoot, scope)) add(path);
    // Shared time grids the archives reference travel with them (TimeGrid.hpp).
    const auto archives = result;
    for (const auto& path : archives) {
        if (path.extension() != ".npz") continue;
        try {
            const auto metadata = readNpzIndexStored(path).metadataJson;
            for (const char* key : {"time_grid_ref", "raw_time_grid_ref"})
                if (const auto grid = resolveTimeGrid(path, metadata, key); !grid.empty()) add(grid);
        } catch (const std::exception&) {
            continue;  // an unreadable archive is still copied; its grid cannot be named
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::filesystem::path createRelocationPlan(const std::filesystem::path& archiveRoot, const std::string& campaignId,
                                           const std::string& destination) {
    if (campaignId.empty()) throw std::invalid_argument("a campaign id is required");
    const auto locations = storageLocations(archiveRoot);
    const auto found = std::find_if(locations.begin(), locations.end(),
                                    [&](const auto& location) { return location.name == destination; });
    if (found == locations.end() || !found->writable || found->name == "primary")
        throw std::invalid_argument("destination must be a registered writable storage location: " + destination);
    std::string safe = campaignId;
    for (auto& character : safe)
        if (!std::isalnum(static_cast<unsigned char>(character)) && character != '-' && character != '_') character = '-';
    const auto path = archiveRoot / L"management" / L"relocations"
        / std::filesystem::u8path("campaign-relocation__" + safe + "__" + utc("%Y%m%dT%H%M%SZ") + ".json");
    writeJsonAtomic(path, json::Value(json::Value::Object{{"schema", json::Value(kCampaignRelocationSchema)},
        {"created_utc", json::Value(utc("%Y-%m-%dT%H:%M:%SZ"))}, {"campaign_id", json::Value(campaignId)},
        {"source_root", json::Value(generic(normalized(archiveRoot)))},
        {"destination_location", json::Value(destination)},
        {"destination_root", json::Value(generic(found->root))}, {"state", json::Value("planned")},
        {"copies", json::Value(json::Value::Array{})}, {"pending_job_ids", json::Value(json::Value::Array{})},
        {"finalized", json::Value(false)}}));
    return path;
}

VerifiedCopy appendRelocationCopy(const std::filesystem::path& planPath, const std::filesystem::path& sourceArgument) {
    auto plan = json::load(planPath.wstring());
    if (plan.at("schema").string() != kCampaignRelocationSchema)
        throw std::invalid_argument(planPath.u8string() + " is not a campaign relocation plan");
    if (plan.at("finalized").boolean()) throw std::runtime_error("the relocation plan is finalized");
    const auto source = normalized(sourceArgument);
    const auto sourceRoot = std::filesystem::u8path(plan.at("source_root").string());
    const auto destinationRoot = std::filesystem::u8path(plan.at("destination_root").string());
    VerifiedCopy copy;
    if (!inside(source, destinationRoot).empty()) {  // written straight to the destination
        if (!isFile(source)) throw std::invalid_argument("copy source is not a file: " + source.u8string());
        copy.source = source; copy.destination = source; copy.alreadyAtDestination = true;
        copy.sha256 = sha256File(source); copy.bytes = std::filesystem::file_size(source);
        copy.copiedUtc = utc("%Y-%m-%dT%H:%M:%SZ");
    } else {
        // Keep the archive layout: <source root>/solutions/x/y.npz -> <destination>/solutions/x/y.npz.
        auto relative = inside(source, sourceRoot);
        if (relative.empty()) relative = std::filesystem::path(L"external-results") / source.filename();
        copy = copyVerified(source, destinationRoot / relative.parent_path(), false);
        copy.relativePath = generic(relative);
    }
    plan.at("state") = json::Value("copying");
    plan.at("copies").array().push_back(copyValue(copy));
    writeJsonAtomic(planPath, plan);
    return copy;
}

std::string storageLocationsJson(const std::vector<StorageLocation>& locations) {
    json::Value::Array values;
    for (const auto& location : locations) values.push_back(locationValue(location));
    return json::dump(json::Value(json::Value::Object{{"schema", json::Value(kStorageRegistrySchema)},
                                                      {"locations", json::Value(std::move(values))}}));
}

std::string verifiedCopyJson(const VerifiedCopy& copy) { return json::dump(copyValue(copy)); }

}  // namespace dynsyslab::core
