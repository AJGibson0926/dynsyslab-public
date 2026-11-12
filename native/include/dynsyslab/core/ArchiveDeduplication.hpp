#pragma once

#include <dynsyslab/core/ArchiveAudit.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ArchiveDeduplicationSkip {
    std::filesystem::path canonical;
    std::filesystem::path candidate;
    std::string reason;
};

struct ArchiveDeduplicationResult {
    std::filesystem::path root;
    std::filesystem::path trashRoot;
    bool dryRun{true};
    std::vector<std::filesystem::path> canonicalFiles;
    std::vector<std::filesystem::path> candidates;
    std::vector<std::filesystem::path> movedFiles;
    std::vector<ArchiveDeduplicationSkip> skipped;
    std::uintmax_t reclaimableBytes{};
    std::uintmax_t movedBytes{};
};

// Conservatively select exact duplicate solutions. The default is read-only.
// Apply mode moves only unreferenced duplicates into archive-local trash so
// every mutation remains recoverable; it never unlinks or rewrites provenance.
[[nodiscard]] ArchiveDeduplicationResult deduplicateArchives(
    const std::filesystem::path& root, const ArchiveAuditScope& scope = {},
    bool apply = false);

[[nodiscard]] std::string archiveDeduplicationJson(
    const ArchiveDeduplicationResult& result);

}  // namespace dynsyslab::core
