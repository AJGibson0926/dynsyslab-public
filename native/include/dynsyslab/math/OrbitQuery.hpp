#pragma once
#include <dynsyslab/core/Json.hpp>
#include <optional>
#include <string>
namespace dynsyslab::math {
// JSON operations of the orbit-family tools; std::nullopt when op is not one of them.
// orbit_family, damped_averages, melnikov_homoclinic, capture_probabilities,
// modal_gramian, orbit_gramian, damped_orbit_gramian, jordan_response.
std::optional<core::json::Value> evaluateOrbitQuery(const std::string& op,const core::json::Value& q);
}
