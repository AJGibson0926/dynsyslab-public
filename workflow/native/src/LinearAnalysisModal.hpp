#pragma once
// linear_analysis "modal-gramian" and "klqe-assessment" (LinearAnalysisModal.cpp).
#include "PostprocessInternal.hpp"

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runModalGramian(const std::filesystem::path& sourcePath, const SourceData& source,
                                      const json::Value& settings, const std::filesystem::path& root,
                                      const std::string& label);
std::filesystem::path runKlqeAssessment(const std::filesystem::path& sourcePath, const SourceData& source,
                                        const json::Value& settings, const std::filesystem::path& root,
                                        const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
