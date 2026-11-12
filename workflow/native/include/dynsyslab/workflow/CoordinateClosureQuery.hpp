#pragma once

#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow {

// Read-only bounded numerical search. The result never asserts a global proof.
[[nodiscard]] json::Value evaluateCoordinateClosureQuery(const json::Value& request);

}  // namespace dynsyslab::workflow
