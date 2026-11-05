#include <dynsyslab/core/Text.hpp>

#ifdef _WIN32
#include <windows.h>
#else
#include <codecvt>
#include <locale>
#endif

namespace dynsyslab::core {

std::wstring wide(const std::string& value) {
    if (value.empty()) return {};
#ifdef _WIN32
    const int count = MultiByteToWideChar(CP_UTF8, 0, value.data(),
                                           static_cast<int>(value.size()), nullptr, 0);
    std::wstring result(static_cast<std::size_t>(count), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, value.data(), static_cast<int>(value.size()),
                        result.data(), count);
    return result;
#else
    return std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>>{}.from_bytes(value);
#endif
}

std::string narrow(const std::wstring& value) {
    if (value.empty()) return {};
#ifdef _WIN32
    const int count = WideCharToMultiByte(CP_UTF8, 0, value.data(),
                                          static_cast<int>(value.size()), nullptr, 0,
                                          nullptr, nullptr);
    std::string result(static_cast<std::size_t>(count), '\0');
    WideCharToMultiByte(CP_UTF8, 0, value.data(), static_cast<int>(value.size()),
                        result.data(), count, nullptr, nullptr);
    return result;
#else
    return std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>>{}.to_bytes(value);
#endif
}

}  // namespace dynsyslab::core
