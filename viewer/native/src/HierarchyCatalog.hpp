#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include <dynsyslab/core/ArchiveIndex.hpp>

#include <filesystem>
#include <set>
#include <string>

namespace dsl {

void appendHierarchyCatalog(Catalog& catalog, const std::filesystem::path& projectRoot,
                            const std::filesystem::path& archiveRoot,
                            std::set<std::string>& ownedSources,
                            const std::vector<dynsyslab::core::ArchiveIndexEntry>*
                                indexEntries = nullptr);


[[nodiscard]] std::filesystem::path resolveHierarchySource(
    const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot,
    const std::filesystem::path& manifest,
    const std::string& source);

[[nodiscard]] std::string normalizedSource(const std::filesystem::path& path);

}  // namespace dsl
