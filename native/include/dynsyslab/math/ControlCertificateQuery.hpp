#pragma once
#include <dynsyslab/core/Json.hpp>

#include <optional>
#include <string>

namespace dynsyslab::math {
// Operations sdp_feasibility, bilinear_fit, residual_sector and certified_design
// of dynsyslab.math-query/v1. Returns nullopt for any other operation, so the
// main dispatcher can chain it. Guide: docs/architecture/CERTIFIED_BILINEAR_CONTROL_TOOLS.md
[[nodiscard]] std::optional<core::json::Value> evaluateControlCertificateQuery(
    const std::string& op, const core::json::Value& q);
}  // namespace dynsyslab::math
