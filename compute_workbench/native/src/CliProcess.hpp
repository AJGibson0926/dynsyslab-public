#pragma once
// Child processes for the CLI front ends (mcp, agent-setup): run a program without a window,
// capture its standard output and error, and return its exit code.
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::cli {

struct ProcessResult {
    int exitCode{-1};
    std::string out, err;
    bool started{};
};

// Quotes one argument the way CommandLineToArgvW parses it.
[[nodiscard]] std::wstring quoteArgument(const std::wstring& value);

// Runs `program args...`; stdin is empty. timeoutMs = 0 waits indefinitely.
[[nodiscard]] ProcessResult runProcess(const std::filesystem::path& program, const std::vector<std::wstring>& args,
                                       unsigned long timeoutMs = 0);

// First match of `name` (e.g. L"codex.exe") on PATH; empty when absent.
[[nodiscard]] std::filesystem::path findOnPath(const std::wstring& name);

}  // namespace dynsyslab::cli
