#pragma once
#include <dynsyslab/core/Json.hpp>
namespace dynsyslab::math {
// Stateless, versioned adapter. Request schema dynsyslab.math-query/v1;
// response dynsyslab.mathematics/v1. No file, archive or solver ownership.
core::json::Value evaluateQuery(const core::json::Value& request);
}
