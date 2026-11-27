#pragma once
// Versioned starter job documents (`dynsyslab-cli template solve|postprocess|estimate`).
// They validate as dynsyslab.compute-job/v1; postprocess and estimate templates still need
// their sources and KROM filled in before they can run.

#include <dynsyslab/workflow/Json.hpp>

#include <string>
#include <vector>

namespace dynsyslab::workflow {

[[nodiscard]] std::vector<std::string> jobTemplateNames();
[[nodiscard]] json::Value jobTemplate(const std::string& name);

}  // namespace dynsyslab::workflow
