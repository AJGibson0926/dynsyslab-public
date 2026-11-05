#pragma once

#include <string>

namespace dynsyslab::core {

[[nodiscard]] std::wstring wide(const std::string& value);
[[nodiscard]] std::string narrow(const std::wstring& value);

}  // namespace dynsyslab::core
