#pragma once

#include <windows.h>

#include <filesystem>
#include <string>

namespace dsl {

inline constexpr wchar_t kUiSupportVersion[] = L"0.1.0";

void showHelpAbout(HWND parent, const std::wstring& application,
                   const std::filesystem::path& projectRoot,
                   const std::filesystem::path& archiveRoot);

}  // namespace dsl
