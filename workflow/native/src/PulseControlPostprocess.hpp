#pragma once
// pulse_response and control_map_audit (PulseControlPostprocess.cpp, docs/architecture/KROM_TOOLS.md).
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runPulseResponse(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label);
std::filesystem::path runControlMapAudit(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                                         const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
