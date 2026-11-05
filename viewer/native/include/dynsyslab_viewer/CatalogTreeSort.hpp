#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <cstddef>
#include <vector>

namespace dsl {

enum class CatalogTreeSort { date, trajectoryCount, totalBytes };

struct CatalogTreeAggregate {
    double newest = 0.0;
    std::size_t trajectoryCount = 0;
    double totalBytes = 0.0;
};

struct OrderedCatalogTree {
    std::vector<std::vector<int>> children;
    std::vector<int> roots;
    std::vector<CatalogTreeAggregate> aggregates;
};

// Computes subtree metrics once, de-duplicates archive sources for byte totals,
// and orders only terminal data records. Structural branches retain catalog order.
OrderedCatalogTree orderedCatalogTree(const Catalog& catalog,
                                      CatalogTreeSort sort);

// Identifies records presented beneath a Training data branch and chooses the
// most useful initial record according to the active leaf ordering. Returns -1
// when the catalog contains no selectable training record.
bool isTrainingDataItem(const Catalog& catalog, int index);
int preferredInitialCatalogItem(const Catalog& catalog,
                                CatalogTreeSort sort);

}  // namespace dsl
