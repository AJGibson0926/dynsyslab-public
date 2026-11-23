#pragma once
#include <dynsyslab/core/Json.hpp>
#include <optional>
#include <string>
namespace dynsyslab::math {
// JSON adapter for the data-driven resolvent tools (docs/architecture/KOOPMAN_RESOLVENT_TOOLS.md).
// Returns nullopt for an operation it does not own; throws std::invalid_argument on bad input.
std::optional<core::json::Value> evaluateResolventDataQuery(const std::string& op,const core::json::Value& q);
}
