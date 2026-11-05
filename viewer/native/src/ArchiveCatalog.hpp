#pragma once

#include <dynsyslab_viewer/Model.hpp>
#include <dynsyslab/core/Archive.hpp>

#include <filesystem>

namespace dsl {

void enrichArchiveCatalogItem(CatalogItem& item,
                              const std::filesystem::path& source);
void enrichArchiveCatalogItem(CatalogItem& item,
                              const dynsyslab::core::NpzIndex& archive);
void resolveArchiveCatalogRelations(Catalog& catalog,
                                    const std::filesystem::path& projectRoot);

}  // namespace dsl
