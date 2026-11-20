#include "CliProcess.hpp"

#include <thread>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::cli {

std::wstring quoteArgument(const std::wstring& value) {
    std::wstring out = L"\"";
    std::size_t slashes = 0;
    for (const wchar_t c : value) {
        if (c == L'\\') { ++slashes; continue; }
        if (c == L'"') out.append(slashes * 2 + 1, L'\\');
        else out.append(slashes, L'\\');
        slashes = 0;
        out.push_back(c);
    }
    out.append(slashes * 2, L'\\');
    return out + L"\"";
}

#ifdef _WIN32
namespace {
void drain(HANDLE pipe, std::string& into) {
    char buffer[4096]; DWORD read = 0;
    while (ReadFile(pipe, buffer, sizeof(buffer), &read, nullptr) && read) into.append(buffer, read);
}
}  // namespace

ProcessResult runProcess(const std::filesystem::path& program, const std::vector<std::wstring>& args,
                         unsigned long timeoutMs) {
    ProcessResult result;
    SECURITY_ATTRIBUTES inherit{sizeof(inherit), nullptr, TRUE};
    HANDLE outRead{}, outWrite{}, errRead{}, errWrite{};
    if (!CreatePipe(&outRead, &outWrite, &inherit, 0) || !CreatePipe(&errRead, &errWrite, &inherit, 0)) return result;
    SetHandleInformation(outRead, HANDLE_FLAG_INHERIT, 0);
    SetHandleInformation(errRead, HANDLE_FLAG_INHERIT, 0);
    HANDLE nul = CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &inherit, OPEN_EXISTING, 0, nullptr);
    std::wstring command = quoteArgument(program.wstring());
    for (const auto& arg : args) command += L" " + quoteArgument(arg);
    STARTUPINFOW startup{}; startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = nul; startup.hStdOutput = outWrite; startup.hStdError = errWrite;
    PROCESS_INFORMATION process{};
    result.started = CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW,
                                    nullptr, nullptr, &startup, &process) != 0;
    CloseHandle(outWrite); CloseHandle(errWrite);
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    if (!result.started) { CloseHandle(outRead); CloseHandle(errRead); return result; }
    std::thread errors([&] { drain(errRead, result.err); });
    drain(outRead, result.out);
    if (WaitForSingleObject(process.hProcess, timeoutMs ? timeoutMs : INFINITE) != WAIT_OBJECT_0)
        TerminateProcess(process.hProcess, 1);
    errors.join();
    DWORD code = 1; GetExitCodeProcess(process.hProcess, &code);
    result.exitCode = static_cast<int>(code);
    CloseHandle(process.hThread); CloseHandle(process.hProcess);
    CloseHandle(outRead); CloseHandle(errRead);
    return result;
}

std::filesystem::path findOnPath(const std::wstring& name) {
    std::vector<wchar_t> buffer(32768);
    const DWORD size = SearchPathW(nullptr, name.c_str(), nullptr, static_cast<DWORD>(buffer.size()), buffer.data(), nullptr);
    return size && size < buffer.size() ? std::filesystem::path(std::wstring(buffer.data(), size)) : std::filesystem::path{};
}
#else
ProcessResult runProcess(const std::filesystem::path&, const std::vector<std::wstring>&, unsigned long) { return {}; }
std::filesystem::path findOnPath(const std::wstring&) { return {}; }
#endif

}  // namespace dynsyslab::cli
