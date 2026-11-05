#pragma once

#include <string>
#include <utility>
#include <vector>

namespace dsl::compute {

// Replaces top-level settings with typed JSON fragments while preserving
// unrelated advanced settings. Malformed source text is retained so the
// canonical workflow validation gate can report it to the user.
std::string mergeSettingsObject(
    std::string object,
    const std::vector<std::pair<std::string, std::string>>& replacements);

}  // namespace dsl::compute
