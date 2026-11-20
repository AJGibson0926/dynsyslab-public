#include "CliResources.hpp"
#include "CliRuntimePaths.hpp"

#include <cstring>
#include <fstream>
#include <utility>

namespace dynsyslab::cli {
namespace fs = std::filesystem;

Resources locateResources() {
    const auto exe = cliExecutablePath();
    const auto viewer = exe.parent_path() / L"DynSysLabViewerCli.exe";
    const auto root = exe.parent_path().parent_path() / L"resources";
    std::error_code ignored;
    if (fs::is_regular_file(root / L"examples" / L"quickstart" / L"quickstart.json", ignored))
        return {root / L"examples" / L"quickstart", root / L"docs" / L"AGENT_WELCOME.md", viewer};
    if (const auto checkout = sourceCheckoutContaining(exe); !checkout.empty())
        return {checkout / L"examples" / L"native" / L"quickstart", checkout / L"docs" / L"AGENT_WELCOME.md", viewer};
    return {root / L"examples" / L"quickstart", root / L"docs" / L"AGENT_WELCOME.md", viewer};
}

std::string asciiText(std::string text) {
    const std::pair<const char*, const char*> folds[]{{"\xE2\x80\x93", "-"}, {"\xE2\x80\x94", "-"},
        {"\xE2\x80\xA6", "..."}, {"\xE2\x86\x92", "->"}, {"\xE2\x80\x98", "'"}, {"\xE2\x80\x99", "'"},
        {"\xE2\x80\x9C", "\""}, {"\xE2\x80\x9D", "\""}, {"\xC3\xA9", "e"}};
    for (const auto& [from, to] : folds)
        for (auto at = text.find(from); at != std::string::npos; at = text.find(from, at + std::strlen(to)))
            text.replace(at, std::strlen(from), to);
    return text;
}

Welcome readWelcome(const fs::path& path) {
    Welcome out;
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        out.overview = "# Welcome to DynSysLab\n\nA native suite for simulating, analysing, modelling and "
                       "controlling dynamical systems.\n";
        return out;
    }
    std::string line; bool overview = true, asks = false;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("## ", 0) == 0) {
            if (line.rfind("## Get started", 0) == 0) overview = false;
            asks = line.rfind("## Things you can ask", 0) == 0;
            if (!overview) continue;
        }
        if (asks && line.rfind("- ", 0) == 0) out.next.push_back(asciiText(line.substr(2)));
        if (overview && line != "---") out.overview += asciiText(line) + "\n";
    }
    while (out.overview.size() > 1 && out.overview.substr(out.overview.size() - 2) == "\n\n") out.overview.pop_back();
    return out;
}

std::string shownPath(const fs::path& path) { return path.lexically_normal().make_preferred().u8string(); }

}  // namespace dynsyslab::cli
