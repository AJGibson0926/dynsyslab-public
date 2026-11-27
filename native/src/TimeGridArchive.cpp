// Reading archives whose time arrays live in a shared grid (core/TimeGrid.hpp).
#include <dynsyslab/core/TimeGrid.hpp>
#include <dynsyslab/core/Json.hpp>

#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
namespace {

constexpr std::pair<const char*, const char*> kReferences[] = {
    {"time_grid_ref", "time"}, {"raw_time_grid_ref", "raw_time"}};

bool mentionsGrid(const std::string& metadataJson) {
    return metadataJson.find("time_grid_ref\"") != std::string::npos;
}

std::filesystem::path safeReference(const std::string& reference) {
    const auto path = std::filesystem::u8path(reference);
    if (reference.empty() || path.is_absolute() || path.has_root_name())
        throw std::runtime_error("time-grid reference must be relative: " + reference);
    for (const auto& part : path)
        if (part == "..") throw std::runtime_error("time-grid reference may not leave its root: " + reference);
    return path;
}

bool isFile(const std::filesystem::path& path) {
    std::error_code ignored;
    return std::filesystem::is_regular_file(path, ignored);
}

}  // namespace

std::filesystem::path resolveTimeGrid(const std::filesystem::path& archive, const std::string& metadataJson,
                                      const std::string& referenceKey) {
    if (!mentionsGrid(metadataJson)) return {};
    const auto metadata = json::parse(metadataJson);
    const auto& reference = metadata.at(referenceKey);
    if (reference.isNull()) return {};
    const auto relative = safeReference(reference.string());
    for (const char* rootKey : {"time_grid_root", "archive_root"}) {
        const auto root = metadata.at(rootKey).string();
        if (!root.empty() && isFile(std::filesystem::u8path(root) / relative))
            return std::filesystem::u8path(root) / relative;
    }
    std::error_code error;
    auto folder = std::filesystem::absolute(archive, error).parent_path();
    for (; !folder.empty(); folder = folder.parent_path()) {
        if (isFile(folder / relative)) return folder / relative;
        if (folder == folder.parent_path()) break;
    }
    throw std::runtime_error(archive.u8string() + " references a shared time grid that is unavailable: "
                             + reference.string());
}

void attachTimeGrids(const std::filesystem::path& archive, NpzRecord& record) {
    if (!record.contains("metadata_json")) return;
    const auto text = record.metadataJson();
    if (!mentionsGrid(text)) return;
    const auto metadata = json::parse(text);
    for (const auto& [key, name] : kReferences) {
        if (record.contains(name)) continue;
        const auto grid = resolveTimeGrid(archive, text, key);
        if (grid.empty()) continue;
        auto stored = readNpzStored(grid);
        if (!stored.contains("time")) throw std::runtime_error(grid.u8string() + " has no time array");
        auto time = std::move(stored.arrays.at("time"));
        const double expected = metadata.at(std::string(name) == "time" ? "time_grid_samples"
                                                                         : "raw_time_grid_samples").number(-1.0);
        const char* values = std::string(name) == "time" ? "values" : "raw_values";
        const bool countMismatch = expected >= 0.0 && static_cast<double>(time.elementCount()) != expected;
        const bool rowMismatch = record.contains(values) && !record.at(values).shape.empty()
            && record.at(values).shape.front() != time.elementCount();
        if (countMismatch || rowMismatch)
            throw std::runtime_error(archive.u8string() + ": shared time grid " + grid.u8string()
                                     + " does not match the archived samples");
        record.arrays[name] = std::move(time);
    }
}

void addTimeGridNames(NpzIndex& index) {
    if (!mentionsGrid(index.metadataJson)) return;
    const auto metadata = json::parse(index.metadataJson);
    for (const auto& [key, name] : kReferences)
        if (!metadata.at(key).isNull() && !index.contains(name)) index.arrays.emplace_back(name);
}

}  // namespace dynsyslab::core
