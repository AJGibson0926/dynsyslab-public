#pragma once

#include <cstddef>
#include <filesystem>
#include <string>

namespace dynsyslab::workflow::detail {

[[nodiscard]] std::string jsonEscape(const std::string& value);
[[nodiscard]] std::string timestamp();
[[nodiscard]] std::string safeName(std::string value,
                                   const std::string& fallback = "native-result");
void atomicProgress(const std::filesystem::path& path, double fraction,
                    std::size_t complete, std::size_t planned,
                    const std::string& message);

}  // namespace dynsyslab::workflow::detail
