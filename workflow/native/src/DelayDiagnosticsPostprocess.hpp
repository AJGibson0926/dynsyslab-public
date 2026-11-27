#pragma once
// Delay-tap stages (DelayDiagnosticsPostprocess.cpp): delay_pairs, delay_consistency.
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runDelayPairs(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runDelayConsistency(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                          const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
