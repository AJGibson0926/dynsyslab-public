#pragma once

#include <dynsyslab_viewer/Model.hpp>

namespace dsl {

// Builds a compact, derived navigation hierarchy. Archive identities, sources,
// relationships, and original labels remain available on the catalog items.
void applyCatalogPresentation(Catalog& catalog);

}  // namespace dsl
