#pragma once
#include <cstddef>
#include <string>
#include <vector>
namespace dynsyslab::core {
struct SystemDescriptor;
// Expand component labels for a resolved state dimension without changing state identities.
[[nodiscard]] std::vector<std::string> systemStateLabels(const SystemDescriptor &,
                                                         std::size_t dimension);
} // namespace dynsyslab::core
