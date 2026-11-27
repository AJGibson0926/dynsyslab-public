#pragma once
// input_map fit "symmetric_pulse" (InputMapPulses.cpp).
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

[[nodiscard]] json::Value::Object symmetricPulseRegression(dynsyslab::core::KoopmanModel& model,
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::string& representation, const std::string& mapName);

}  // namespace dynsyslab::workflow::postprocess_detail
