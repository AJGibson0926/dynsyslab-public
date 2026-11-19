#include "CliViewerFigures.hpp"

#include <dynsyslab/core/LongPath.hpp>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::cli {
namespace {

std::wstring quoted(const std::wstring& value) {  // CommandLineToArgvW quoting
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

}  // namespace

bool renderFigure(const std::filesystem::path& viewerCli, const FigureRequest& request,
                  const std::filesystem::path& cache, std::string& reason) {
    std::error_code ignored;
    if (!std::filesystem::is_regular_file(viewerCli, ignored)) {
        reason = "Viewer CLI not found next to dynsyslab-cli";
        return false;
    }
    std::filesystem::create_directories(core::extendedLengthPath(request.output.parent_path()), ignored);
#ifdef _WIN32
    std::wstring command = quoted(viewerCli.wstring()) + L" render --archives "
        + quoted(request.archive.parent_path().wstring()) + L" --cache " + quoted(cache.wstring())
        + L" --source " + quoted(request.archive.wstring()) + L" --view " + quoted(std::wstring(request.view.begin(), request.view.end()))
        + L" --width 1200 --height 800 --output " + quoted(request.output.wstring());
    for (const auto& option : request.extra) command += L" " + quoted(option);
    SECURITY_ATTRIBUTES inherit{sizeof(inherit), nullptr, TRUE};
    HANDLE nul = CreateFileW(L"NUL", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                             &inherit, OPEN_EXISTING, 0, nullptr);
    STARTUPINFOW startup{}; startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = startup.hStdOutput = startup.hStdError = nul;
    PROCESS_INFORMATION process{};
    const BOOL started = CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW,
                                        nullptr, nullptr, &startup, &process);
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    if (!started) { reason = "unable to start the Viewer CLI"; return false; }
    const DWORD waited = WaitForSingleObject(process.hProcess, 60000);
    DWORD code = 1;
    if (waited == WAIT_OBJECT_0) GetExitCodeProcess(process.hProcess, &code);
    else TerminateProcess(process.hProcess, 1);
    CloseHandle(process.hThread); CloseHandle(process.hProcess);
    if (waited != WAIT_OBJECT_0) { reason = "the Viewer CLI timed out"; return false; }
    if (code != 0 || !std::filesystem::is_regular_file(core::extendedLengthPath(request.output), ignored)) {
        reason = "the Viewer CLI could not render this view";
        return false;
    }
    return true;
#else
    reason = "figures need the Windows Viewer CLI";
    return false;
#endif
}

}  // namespace dynsyslab::cli
