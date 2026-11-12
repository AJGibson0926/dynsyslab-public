#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <string>

namespace dsl {

// Combines two compatible Viewer payloads under one synchronized page model.
[[nodiscard]] Payload comparePayload(const Payload& first, const std::string& firstLabel,
                                     const Payload& second, const std::string& secondLabel);

}  // namespace dsl
