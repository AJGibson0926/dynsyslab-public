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

// Copies and consistently reduces the evidence needed to inspect a KROM
// without reopening the KROM from which it was extracted or regressed.
void inheritKromEvidence(
    const dynsyslab::core::NpzRecord& parent,
    const std::vector<std::size_t>& selectedModes,
    std::map<std::string, dynsyslab::core::NpyArray>& output);

}  // namespace dynsyslab::workflow::postprocess_detail
