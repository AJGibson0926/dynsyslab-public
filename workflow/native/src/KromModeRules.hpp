#pragma once
// krom_extract mode-selection rules state_observability, nearest_frequency and validation_score
// (KromModeRules.cpp).
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

[[nodiscard]] bool isExtendedModeRule(const std::string& rule);
// Mode indices into `model`, conjugates completed; `report` receives the rule's evidence.
[[nodiscard]] std::vector<std::size_t> extendedModeRule(const json::Value& rule,
    const dynsyslab::core::KoopmanModel& model, const dynsyslab::core::NpzRecord& parentArchive,
    json::Value& report);

}  // namespace dynsyslab::workflow::postprocess_detail
