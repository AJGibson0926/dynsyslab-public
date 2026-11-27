// Lossless archive storage maintenance (core/ArchiveMaintenance.hpp).
#include <dynsyslab/core/ArchiveMaintenance.hpp>

#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/TimeGrid.hpp>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dynsyslab::core {
namespace {

std::string isoNow() {
    const std::time_t raw = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

bool sameArray(const NpyArray& a, const NpyArray& b) {
    return a.descriptor == b.descriptor && a.shape == b.shape && a.data == b.data;
}

// Every array of `expected` except `skip` is present and identical in `actual`, and `actual`
// has no other arrays (except `skip`).
bool sameArrays(const NpzRecord& expected, const NpzRecord& actual, const std::string& skip = {}) {
    std::size_t compared = 0;
    for (const auto& [name, array] : expected.arrays) {
        if (name == skip) continue;
        if (!actual.contains(name) || !sameArray(array, actual.at(name))) return false;
        ++compared;
    }
    std::size_t other = 0;
    for (const auto& item : actual.arrays) other += item.first == skip ? 0 : 1;
    return compared == other;
}

std::uintmax_t fileBytes(const std::filesystem::path& path) {
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : size;
}

// Writes `record` beside `path`, checks it with `verify`, then replaces `path` with it.
template <typename Verify>
void replaceVerified(const std::filesystem::path& path, const NpzRecord& record,
                     ArchiveCompression compression, Verify verify) {
    const auto candidate = path.parent_path()
        / (path.stem().wstring() + L".rewrite" + path.extension().wstring());
    std::error_code ignored;
    std::filesystem::remove(candidate, ignored);
    try {
        writeNpz(candidate, record, {compression, false, false});
        if (!verify(candidate)) throw std::runtime_error("verification failed; the original is unchanged");
        std::filesystem::rename(candidate, path);
    } catch (...) {
        std::filesystem::remove(candidate, ignored);
        throw;
    }
}

std::string timeKey(const NpyArray& time) {
    Sha256Stream digest;
    digest.update(time.descriptor);
    for (const auto extent : time.shape) digest.update(std::to_string(extent) + ",");
    digest.update(time.data.data(), time.data.size());
    return digest.finish();
}

json::Value metadataOf(const NpzRecord& record) {
    const auto text = record.metadataJson();
    return text.empty() ? json::Value(json::Value::Object{}) : json::parse(text);
}

bool protectedKey(const std::string& key) {
    return key == "recorded_archive" || key == "run_provenance" || key == "metadata_edits"
        || key == "metadata_revision" || key.rfind("time_grid_", 0) == 0
        || key.rfind("raw_time_grid_", 0) == 0;
}

// Writes the shared grid file once; an existing one must hold the same samples.
void ensureGrid(const std::filesystem::path& grid, const NpyArray& time, const std::string& key) {
    std::error_code missing;
    if (!std::filesystem::is_regular_file(grid, missing)) {
        NpzRecord record;
        record.arrays["time"] = time;
        json::Value::Object meta{{"schema", json::Value(kTimeGridSchema)}, {"sha256", json::Value(key)},
            {"samples", json::Value(static_cast<double>(time.elementCount()))}};
        record.arrays["metadata_json"] = NpyArray::fromString(json::dump(json::Value(std::move(meta))));
        writeNpz(grid, record, {ArchiveCompression::deflate, false, false});
    }
    const auto check = readNpzStored(grid);
    if (!check.contains("time") || !sameArray(check.at("time"), time))
        throw std::runtime_error("existing shared grid " + grid.u8string() + " differs");
}

}  // namespace

ArchiveMaintenanceResult recompressArchives(const std::vector<std::filesystem::path>& archives,
                                            ArchiveCompression compression, bool apply) {
    ArchiveMaintenanceResult result{"dynsyslab.archive-recompress/v1", {}, !apply,
                                    archiveCompressionName(compression), {}, {}};
    for (const auto& path : archives) {
        ArchiveRewrite item{path, fileBytes(path), 0, false, {}};
        try {
            if (apply) {
                const auto stored = readNpzStored(path);
                replaceVerified(path, stored, compression, [&](const std::filesystem::path& candidate) {
                    return sameArrays(stored, readNpzStored(candidate));
                });
                item.rewritten = true;
            }
            item.bytesAfter = fileBytes(path);
        } catch (const std::exception& error) {
            item.error = error.what();
            item.bytesAfter = item.bytesBefore;
        }
        result.files.push_back(std::move(item));
    }
    return result;
}

ArchiveMaintenanceResult shareTimeGrids(const std::filesystem::path& rootArgument,
                                        const std::vector<std::filesystem::path>& archives, bool apply) {
    const auto root = std::filesystem::absolute(rootArgument).lexically_normal();
    ArchiveMaintenanceResult result{"dynsyslab.time-grid-sharing/v1", root, !apply,
                                    archiveCompressionName(archiveCompression()), {}, {}};
    std::map<std::string, std::vector<std::filesystem::path>> groups;
    for (const auto& path : archives) {
        try {
            const auto stored = readNpzStored(path);
            if (!stored.contains("time") || isRecordedArchive(stored.metadataJson())
                    || !metadataOf(stored).at("time_grid_ref").isNull()) continue;
            groups[timeKey(stored.at("time"))].push_back(path);
        } catch (const std::exception& error) {
            result.files.push_back({path, fileBytes(path), fileBytes(path), false, error.what()});
        }
    }
    for (const auto& [key, members] : groups) {
        if (members.size() < 2) continue;
        const auto reference = std::string(kTimeGridFolder) + "/" + key + ".npz";
        const auto grid = root / std::filesystem::u8path(reference);
        result.grids.push_back(grid);
        for (const auto& path : members) {
            ArchiveRewrite item{path, fileBytes(path), fileBytes(path), false, {}};
            try {
                if (apply) {
                    auto stored = readNpzStored(path);
                    if (!stored.contains("time") || timeKey(stored.at("time")) != key)
                        throw std::runtime_error("the archive changed during sharing");
                    ensureGrid(grid, stored.at("time"), key);
                    const auto original = stored;
                    const auto samples = static_cast<double>(stored.at("time").elementCount());
                    stored.arrays.erase("time");
                    auto metadata = metadataOf(stored);
                    metadata.at("time_grid_ref") = json::Value(reference);
                    metadata.at("time_grid_root") = json::Value(root.generic_u8string());
                    metadata.at("time_grid_shared") = json::Value(true);
                    metadata.at("time_grid_samples") = json::Value(samples);
                    stored.arrays["metadata_json"] = NpyArray::fromString(json::dump(metadata));
                    replaceVerified(path, stored, archiveCompression(), [&](const std::filesystem::path& candidate) {
                        return sameArrays(original, readNpz(candidate), "metadata_json");
                    });
                    item.rewritten = true;
                    item.bytesAfter = fileBytes(path);
                }
            } catch (const std::exception& error) {
                item.error = error.what();
            }
            result.files.push_back(std::move(item));
        }
    }
    return result;
}

std::string archiveMaintenanceJson(const ArchiveMaintenanceResult& result) {
    json::Value::Array files, grids;
    std::uintmax_t before{}, after{};
    std::size_t rewritten{}, failed{};
    for (const auto& item : result.files) {
        json::Value::Object entry{{"path", json::Value(item.path.generic_u8string())},
            {"bytes_before", json::Value(static_cast<double>(item.bytesBefore))},
            {"bytes_after", json::Value(static_cast<double>(item.bytesAfter))},
            {"rewritten", json::Value(item.rewritten)}};
        if (!item.error.empty()) {
            entry["error"] = json::Value(item.error);
            ++failed;
        }
        files.emplace_back(std::move(entry));
        before += item.bytesBefore;
        after += item.bytesAfter;
        rewritten += item.rewritten ? 1 : 0;
    }
    for (const auto& grid : result.grids) grids.emplace_back(grid.generic_u8string());
    json::Value::Object out{{"schema", json::Value(result.schema)}, {"dry_run", json::Value(result.dryRun)},
        {"compression", json::Value(result.compression)}, {"files", json::Value(std::move(files))},
        {"rewritten", json::Value(static_cast<double>(rewritten))},
        {"failed", json::Value(static_cast<double>(failed))},
        {"bytes_before", json::Value(static_cast<double>(before))},
        {"bytes_after", json::Value(static_cast<double>(after))}};
    if (!result.root.empty()) out["root"] = json::Value(result.root.generic_u8string());
    if (result.schema == "dynsyslab.time-grid-sharing/v1") out["grids"] = json::Value(std::move(grids));
    return json::dump(json::Value(std::move(out)));
}

std::filesystem::path editArchiveMetadata(const std::filesystem::path& source,
                                          const std::filesystem::path& outputArgument,
                                          const std::string& patchJson) {
    const auto output = prepareOutputFile(outputArgument, "edit-metadata output");
    std::error_code ignored;
    if (std::filesystem::exists(output, ignored))
        throw std::invalid_argument("edit-metadata output already exists (archives are immutable): "
                                    + output.u8string());
    const auto patch = json::parse(patchJson);
    if (!std::holds_alternative<json::Value::Object>(patch.storage()))
        throw std::invalid_argument("the metadata patch must be a JSON object");
    auto record = readNpzStored(source);
    if (!record.contains("metadata_json"))
        throw std::invalid_argument("the source archive has no metadata_json");
    auto metadata = metadataOf(record);
    json::Value::Array keys;
    for (const auto& [key, value] : patch.object()) {
        if (protectedKey(key))
            throw std::invalid_argument("metadata key '" + key
                                        + "' defines how the archive is read and cannot be edited");
        if (value.isNull()) metadata.object().erase(key);
        else metadata.at(key) = value;
        keys.emplace_back(key);
    }
    json::Value::Array edits;
    if (const auto& previous = metadata.at("metadata_edits");
            std::holds_alternative<json::Value::Array>(previous.storage()))
        edits = previous.array();
    edits.emplace_back(json::Value::Object{
        {"source", json::Value(std::filesystem::absolute(source).generic_u8string())},
        {"source_sha256", json::Value(sha256File(source))}, {"keys", json::Value(std::move(keys))},
        {"edited_utc", json::Value(isoNow())}});
    metadata.at("metadata_edits") = json::Value(std::move(edits));
    metadata.at("metadata_revision") = json::Value(metadata.at("metadata_revision").number(0.0) + 1.0);
    record.arrays["metadata_json"] = NpyArray::fromString(json::dump(metadata));
    writeNpz(output, record, {archiveCompression(), false, false});
    return output;
}

}  // namespace dynsyslab::core
