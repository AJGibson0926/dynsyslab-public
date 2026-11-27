#pragma once

#include <dynsyslab/core/Archive.hpp>

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

[[nodiscard]] std::vector<std::size_t> inheritedSourceModes(
    const dynsyslab::core::NpzRecord& parent,
    const std::vector<std::size_t>& selectedModes,
    std::size_t outputModes);

}  // namespace dynsyslab::workflow::postprocess_detail
