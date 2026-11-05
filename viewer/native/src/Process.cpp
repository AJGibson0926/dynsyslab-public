#include "dynsyslab_viewer/Process.hpp"

#include <windows.h>
#include <shellapi.h>

#include <sstream>

namespace dsl {
namespace {

std::wstring quote(const std::wstring& value) {
    if (value.find_first_of(L" \t\"") == std::wstring::npos) return value;
    std::wstring result = L"\"";
    unsigned slashes = 0;
    for (const wchar_t current : value) {
        if (current == L'\\') {
            ++slashes;
        } else if (current == L'"') {
            result.append(slashes * 2 + 1, L'\\');
            result.push_back(L'"');
            slashes = 0;
        } else {
            result.append(slashes, L'\\');
            slashes = 0;
            result.push_back(current);
        }
    }
    result.append(slashes * 2, L'\\');
    result.push_back(L'"');
    return result;
}

std::wstring systemMessage(DWORD code) {
    wchar_t* buffer = nullptr;
    const DWORD size = FormatMessageW(FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, nullptr, code, 0,
        reinterpret_cast<wchar_t*>(&buffer), 0, nullptr);
    std::wstring result = size && buffer ? std::wstring(buffer, size) : L"Unknown Windows error";
    if (buffer) LocalFree(buffer);
    return result;
}

}  // namespace

ProcessResult runHidden(const std::vector<std::wstring>& arguments,
                        const std::wstring& workingDirectory) {
    if (arguments.empty()) return {ERROR_INVALID_PARAMETER, L"No process was specified"};
    std::wostringstream command;
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        if (index) command << L' ';
        command << quote(arguments[index]);
    }
    std::wstring mutableCommand = command.str();
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESHOWWINDOW;
    startup.wShowWindow = SW_HIDE;
    PROCESS_INFORMATION process{};
    const BOOL created = CreateProcessW(nullptr, mutableCommand.data(), nullptr, nullptr,
        FALSE, CREATE_NO_WINDOW, nullptr, workingDirectory.c_str(), &startup, &process);
    if (!created) {
        const DWORD code = GetLastError();
        return {code, systemMessage(code)};
    }
    WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exitCode = 1;
    GetExitCodeProcess(process.hProcess, &exitCode);
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return {exitCode, exitCode ? L"Viewer data bridge exited with code " +
        std::to_wstring(exitCode) : L""};
}

ProcessResult launchDetached(const std::vector<std::wstring>& arguments,
                             const std::wstring& workingDirectory) {
    if (arguments.empty()) return {ERROR_INVALID_PARAMETER, L"No process was specified"};
    std::wostringstream command;
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        if (index) command << L' ';
        command << quote(arguments[index]);
    }
    std::wstring mutableCommand = command.str();
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process{};
    const BOOL created = CreateProcessW(arguments.front().c_str(), mutableCommand.data(),
        nullptr, nullptr, FALSE, CREATE_NEW_PROCESS_GROUP, nullptr,
        workingDirectory.empty() ? nullptr : workingDirectory.c_str(), &startup, &process);
    if (!created) {
        const DWORD code = GetLastError();
        return {code, systemMessage(code)};
    }
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    return {};
}

bool openDocument(const std::wstring& path) {
    return reinterpret_cast<INT_PTR>(ShellExecuteW(nullptr, L"open", path.c_str(),
        nullptr, nullptr, SW_SHOWNORMAL)) > 32;
}

}  // namespace dsl
