#pragma once

#include <dynsyslab/core/Extension.hpp>
#include <dynsyslab/extension_api.h>

#include <vector>

namespace dynsyslab::core::detail {

// The delay, field and lattice models an extension declares (empty for extensions built
// against the ODE-only v1 header, whose struct_size ends at execute_json).
[[nodiscard]] std::vector<ExtensionModel> extensionModels(const dsl_extension_v1& extension);

}  // namespace dynsyslab::core::detail
