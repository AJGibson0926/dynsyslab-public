#pragma once
// Pooled K-ROM tool methods: control_error_compare, krom_rollout, pod_galerkin, krom_refit,
// reconstruction_transfer (KromTools.cpp).
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

[[nodiscard]] bool isKromToolMethod(const std::string& name);
[[nodiscard]] std::filesystem::path runKromTool(const std::string& name,
    const std::vector<std::filesystem::path>& sources, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
