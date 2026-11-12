#include <dynsyslab/core/ArchiveDeduplication.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
namespace {
namespace json = dynsyslab::core::json;

std::string timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &raw);
#else
    gmtime_r(&raw, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y%m%dT%H%M%SZ");
    return out.str();
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to inspect JSON provenance");
    std::ostringstream out; out << stream.rdbuf(); return out.str();
}

bool containsAny(const std::string& text, const std::vector<std::string>& tokens) {
    return std::any_of(tokens.begin(), tokens.end(), [&](const auto& token) {
        return !token.empty() && text.find(token) != std::string::npos;
    });
}

struct ReferenceCorpus {
    std::vector<std::pair<std::filesystem::path, std::string>> documents;
    bool incomplete{};
};

ReferenceCorpus loadReferences(const std::filesystem::path& root) {
    ReferenceCorpus result;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file()) continue;
        const auto extension = entry.path().extension().u8string();
        if (extension != ".json" && extension != ".npz") continue;
        try {
            const auto text = extension == ".json"
                ? readText(entry.path()) : readNpz(entry.path()).metadataJson();
            result.documents.emplace_back(entry.path(), text);
        } catch (...) {
            // Cleanup must become more conservative, never less, when a
            // provenance-bearing artifact cannot be inspected.
            result.incomplete = true;
        }
    }
    return result;
}

std::vector<std::string> referenceTokens(const std::filesystem::path& root,
                                         const std::filesystem::path& path) {
    std::vector<std::string> result{path.generic_u8string(), path.filename().u8string()};
    std::error_code error;
    const auto absolute = std::filesystem::absolute(path, error);
    if (!error) result.push_back(absolute.generic_u8string());
    const auto relative = std::filesystem::relative(path, root, error);
    if (!error) result.push_back(relative.generic_u8string());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool referenced(const ReferenceCorpus& corpus, const std::filesystem::path& root,
                const std::filesystem::path& candidate) {
    if (corpus.incomplete) return true;
    const auto tokens = referenceTokens(root, candidate);
    return std::any_of(corpus.documents.begin(), corpus.documents.end(),
        [&](const auto& document) {
            return document.first != candidate && containsAny(document.second, tokens);
        });
}

std::vector<std::filesystem::path> newestFirst(
        std::vector<std::filesystem::path> paths) {
    std::sort(paths.begin(), paths.end(), [](const auto& left, const auto& right) {
        const auto leftTime = std::filesystem::last_write_time(left);
        const auto rightTime = std::filesystem::last_write_time(right);
        if (leftTime != rightTime) return leftTime > rightTime;
        return left.generic_u8string() < right.generic_u8string();
    });
    return paths;
}

json::Value pathsJson(const std::vector<std::filesystem::path>& paths) {
    json::Value::Array result;
    for (const auto& path : paths) result.emplace_back(path.generic_u8string());
    return json::Value(std::move(result));
}

void writeManifest(const ArchiveDeduplicationResult& result) {
    std::filesystem::create_directories(result.trashRoot);
    const auto destination = result.trashRoot / "deduplication.json";
    const auto temporary = destination.parent_path() / "deduplication.json.partial";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("unable to write deduplication recovery manifest");
    out << archiveDeduplicationJson(result); out.close();
    std::filesystem::rename(temporary, destination);
}

}  // namespace

ArchiveDeduplicationResult deduplicateArchives(
        const std::filesystem::path& root, const ArchiveAuditScope& scope,
        bool apply) {
    ArchiveDeduplicationResult result;
    result.root = std::filesystem::absolute(root);
    result.dryRun = !apply;
    const auto audit = auditArchives(result.root, scope);
    const auto trash = result.root / "management" / "trash"
        / ("deduplicate-" + timestamp());
    const auto references = audit.exactDuplicateGroups.empty()
        ? ReferenceCorpus{} : loadReferences(result.root);
    for (const auto& group : audit.exactDuplicateGroups) {
        const auto ordered = newestFirst(group.files);
        if (ordered.empty()) continue;
        result.canonicalFiles.push_back(ordered.front());
        for (std::size_t index = 1; index < ordered.size(); ++index) {
            const auto& candidate = ordered[index];
            if (referenced(references, result.root, candidate)) {
                result.skipped.push_back({ordered.front(), candidate,
                    "referenced by archive provenance or unreadable provenance"});
                continue;
            }
            result.candidates.push_back(candidate);
            result.reclaimableBytes += std::filesystem::file_size(candidate);
        }
    }
    if (!apply || result.candidates.empty()) return result;
    result.trashRoot = trash;
    for (const auto& candidate : result.candidates) {
        std::error_code error;
        const auto relative = std::filesystem::relative(candidate, result.root, error);
        if (error || relative.empty() || relative.is_absolute()
                || *relative.begin() == "..")
            throw std::runtime_error("deduplication candidate escaped the archive root");
        const auto destination = trash / relative;
        if (std::filesystem::exists(destination))
            throw std::runtime_error("deduplication trash destination already exists");
    }
    std::filesystem::create_directories(trash);
    for (const auto& candidate : result.candidates) {
        const auto relative = std::filesystem::relative(candidate, result.root);
        const auto destination = trash / relative;
        std::filesystem::create_directories(destination.parent_path());
        const auto bytes = std::filesystem::file_size(candidate);
        std::filesystem::rename(candidate, destination);
        result.movedFiles.push_back(destination); result.movedBytes += bytes;
    }
    writeManifest(result);
    return result;
}

std::string archiveDeduplicationJson(const ArchiveDeduplicationResult& result) {
    json::Value::Array skipped;
    for (const auto& item : result.skipped) {
        json::Value::Object value;
        value["canonical"] = item.canonical.generic_u8string();
        value["candidate"] = item.candidate.generic_u8string();
        value["reason"] = item.reason; skipped.emplace_back(std::move(value));
    }
    json::Value::Object value;
    value["schema"] = "dynsyslab.archive-deduplication/v1";
    value["root"] = result.root.generic_u8string();
    value["dry_run"] = result.dryRun;
    value["recoverable"] = true;
    value["trash_root"] = result.trashRoot.empty()
        ? json::Value{} : json::Value(result.trashRoot.generic_u8string());
    value["canonical_files"] = pathsJson(result.canonicalFiles);
    value["candidates"] = pathsJson(result.candidates);
    value["moved_files"] = pathsJson(result.movedFiles);
    value["skipped_groups"] = json::Value(std::move(skipped));
    value["reclaimable_bytes"] = static_cast<double>(result.reclaimableBytes);
    value["moved_bytes"] = static_cast<double>(result.movedBytes);
    return json::dump(json::Value(std::move(value)));
}

}  // namespace dynsyslab::core
