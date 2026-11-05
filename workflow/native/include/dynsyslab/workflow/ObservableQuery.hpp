#pragma once
#include <dynsyslab/workflow/Json.hpp>
namespace dynsyslab::workflow {
// Read-only expression/derivative inspection in native system coordinates.
[[nodiscard]] json::Value evaluateObservableQuery(const json::Value& request);
}
