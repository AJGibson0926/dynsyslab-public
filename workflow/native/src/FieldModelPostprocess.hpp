#pragma once
// field_model: data-driven vector fields of every field-model kind (FieldModelPostprocess.cpp,
// docs/architecture/FIELD_MODELS.md).
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runFieldModel(const std::vector<std::filesystem::path>& sources, const dynsyslab::core::json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
