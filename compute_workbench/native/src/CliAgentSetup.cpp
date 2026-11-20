#include "CliAgentSetup.hpp"
#include "CliHub.hpp"
#include "CliProcess.hpp"
#include "CliResources.hpp"
#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/Version.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

namespace dynsyslab::cli {
namespace {
namespace fs = std::filesystem;
namespace json = core::json;
constexpr const char* kMarker = "<!-- dynsyslab-workspace: written by dynsyslab-cli agent-setup; "
                                "delete this line to keep your own edits -->";

struct Step { std::string agent, action, detail; bool ok{true}; };

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
void writeFile(const fs::path& path, const std::string& text) {
    std::ofstream(core::prepareOutputFile(path, "agent configuration"), std::ios::binary) << text;
}

std::string pretty(const json::Value& value, int depth = 0) {  // settings files stay readable
    const std::string pad(static_cast<std::size_t>(depth + 1) * 2, ' '), end(static_cast<std::size_t>(depth) * 2, ' ');
    if (std::holds_alternative<json::Value::Object>(value.storage())) {
        const auto& object = value.object();
        if (object.empty()) return "{}";
        std::string out = "{\n"; bool first = true;
        for (const auto& [key, item] : object) {
            out += (first ? "" : ",\n") + pad + json::dump(json::Value(key)) + ": " + pretty(item, depth + 1); first = false;
        }
        return out + "\n" + end + "}";
    }
    if (std::holds_alternative<json::Value::Array>(value.storage())) {
        const auto& array = value.array();
        if (array.empty()) return "[]";
        std::string out = "[\n";
        for (std::size_t i = 0; i < array.size(); ++i) out += (i ? ",\n" : "") + pad + pretty(array[i], depth + 1);
        return out + "\n" + end + "]";
    }
    return json::dump(value);
}

fs::path claudeCommand() {
    if (auto found = findOnPath(L"claude.exe"); !found.empty()) return found;
    fs::path best; fs::file_time_type newest{}; std::error_code error;  // the desktop app's bundled copy
    const auto root = userDataRoot().parent_path() / L"Claude-Data";
    for (fs::recursive_directory_iterator it(root, fs::directory_options::skip_permission_denied, error), end;
         !error && it != end; it.increment(error)) {
        if (it.depth() > 4) { it.disable_recursion_pending(); continue; }
        const auto& path = it->path();
        if (path.filename() == L"claude.exe" && path.wstring().find(L"claude-code") != std::wstring::npos) {
            const auto time = fs::last_write_time(path, error);
            if (best.empty() || time > newest) { best = path; newest = time; }
        }
    }
    return best;
}

fs::path codexCommand() {
    if (auto found = findOnPath(L"codex.exe"); !found.empty()) return found;
    const auto bundled = userDataRoot().parent_path() / L"Programs" / L"OpenAI" / L"Codex" / L"bin" / L"codex.exe";
    std::error_code ignored;
    return fs::is_regular_file(bundled, ignored) ? bundled : fs::path{};
}

void registerMcp(const fs::path& agent, std::vector<std::wstring> scope, const char* name, bool remove, std::vector<Step>& steps) {
    auto removeArgs = std::vector<std::wstring>{L"mcp", L"remove"};
    removeArgs.insert(removeArgs.end(), scope.begin(), scope.end()); removeArgs.push_back(L"dynsyslab");
    (void)runProcess(agent, removeArgs, 60000);  // absent is fine
    if (remove) { steps.push_back({name, "removed MCP server", "dynsyslab"}); return; }
    auto addArgs = std::vector<std::wstring>{L"mcp", L"add"};
    addArgs.insert(addArgs.end(), scope.begin(), scope.end());
    addArgs.insert(addArgs.end(), {L"dynsyslab", L"--", cliExecutablePath().wstring(), L"mcp"});
    const auto run = runProcess(agent, addArgs, 60000);
    steps.push_back({name, "registered MCP server", run.exitCode == 0 ? "dynsyslab -> dynsyslab-cli mcp"
        : "failed: " + run.err + run.out, run.exitCode == 0});
}

bool ourHook(const json::Value& entry) {
    if (!std::holds_alternative<json::Value::Object>(entry.storage()) || entry.at("hooks").isNull()) return false;
    for (const auto& hook : entry.at("hooks").array()) {
        const auto command = hook.at("command").string();
        if (command.find("dynsyslab-cli") != std::string::npos && command.find("hub --hook") != std::string::npos) return true;
    }
    return false;
}

void claudeHook(const fs::path& home, bool remove, std::vector<Step>& steps) {
    const auto file = home / L".claude" / L"settings.json";
    std::error_code ignored;
    json::Value settings = json::Value::Object{};
    if (fs::is_regular_file(file, ignored)) settings = json::load(file.wstring());
    auto& root = settings.object();
    json::Value::Array kept;
    if (root.count("hooks") && root["hooks"].object().count("SessionStart"))
        for (const auto& entry : root["hooks"].object()["SessionStart"].array()) if (!ourHook(entry)) kept.push_back(entry);
    if (!remove) {
        const auto command = quoteArgument(cliExecutablePath().wstring()) + L" hub --hook";
        kept.push_back(json::Value::Object{{"matcher", "startup"}, {"hooks", json::Value::Array{json::Value::Object{
            {"type", "command"}, {"command", core::narrow(command)}, {"timeout", 30.0}}}}});
    }
    if (!root.count("hooks")) root["hooks"] = json::Value::Object{};
    if (kept.empty()) root["hooks"].object().erase("SessionStart");
    else root["hooks"].object()["SessionStart"] = json::Value(std::move(kept));
    if (root["hooks"].object().empty()) root.erase("hooks");
    const auto backup = fs::path(file.wstring() + L".dynsyslab-backup");
    if (fs::is_regular_file(file, ignored) && !fs::exists(backup, ignored)) fs::copy_file(file, backup, ignored);
    writeFile(file, pretty(settings) + "\n");
    steps.push_back({"Claude Code", remove ? "removed start-up hook" : "added start-up hook", shownPath(file)});
}

void codexToml(const fs::path& home, bool remove, std::vector<Step>& steps) {
    const auto file = home / L".codex" / L"config.toml";
    std::istringstream in(readFile(file));
    std::string line, out; bool skipping = false;
    while (std::getline(in, line)) {
        const auto trimmed = line.substr(0, line.find_last_not_of("\r \t") + 1);
        if (!trimmed.empty() && trimmed.front() == '[')
            skipping = trimmed == "[mcp_servers.dynsyslab]" || trimmed.rfind("[mcp_servers.dynsyslab.", 0) == 0;
        if (!skipping) out += line + "\n";
    }
    if (!remove) {
        if (!out.empty() && out.substr(out.size() - std::min<std::size_t>(out.size(), 2)) != "\n\n") out += "\n";
        out += "[mcp_servers.dynsyslab]\ncommand = '" + core::narrow(cliExecutablePath().wstring()) + "'\nargs = [\"mcp\"]\n";
    }
    writeFile(file, out);
    steps.push_back({"Codex", remove ? "removed MCP server" : "registered MCP server", shownPath(file)});
}

void workspaceFiles(const fs::path& workspace, bool remove, std::vector<Step>& steps) {
    std::error_code ignored;
    const auto archives = defaultArchiveRoot().path;
    const std::string agents = std::string(kMarker) + "\n# DynSysLab workspace\n\n"
        "This folder is the user's DynSysLab workspace (DynSysLab " + std::string(core::version) + ", installed with "
        "`" + shownPath(cliExecutablePath()) + "`).\n\nAt the start of every session here, open the DynSysLab Hub: call the "
        "`hub` tool of the `dynsyslab` MCP server (or run `dynsyslab-cli hub`), show its greeting, overview and numbered "
        "options in the chat, and ask the user what they would like to do. Do the work with the `dynsyslab` tools or "
        "`dynsyslab-cli`; every run generates its own data in `" + shownPath(archives) + "`. Start the desktop "
        "applications (Solver/Post-processor, Viewer) only when the user asks.\n";
    const std::pair<const wchar_t*, std::string> files[]{{L"AGENTS.md", agents},
        {L"CLAUDE.md", std::string(kMarker) + "\n@AGENTS.md\n"}, {kWorkspaceMarker, "DynSysLab workspace\n"}};
    for (const auto& [name, contents] : files) {
        const auto path = workspace / name;
        const bool exists = fs::exists(path, ignored);
        const bool ours = !exists || name == std::wstring(kWorkspaceMarker) || readFile(path).rfind(kMarker, 0) == 0;
        if (!ours) { steps.push_back({"Workspace", "kept your file", shownPath(path)}); continue; }
        if (remove) { if (exists) fs::remove(path, ignored); continue; }
        writeFile(path, contents);
    }
    steps.push_back({"Workspace", remove ? "removed DynSysLab files" : "ready", shownPath(workspace)});
}

}  // namespace

int agentSetup(const AgentSetupOptions& options) {
    std::vector<Step> steps;
    std::error_code ignored;
    // %USERPROFILE%: the parent of the user's local application-data folder's "AppData".
    const auto home = !options.home.empty() ? options.home : userDataRoot().parent_path().parent_path().parent_path();
    const auto claude = options.useAgentCommands ? claudeCommand() : fs::path{};
    if (!claude.empty()) registerMcp(claude, {L"--scope", L"user"}, "Claude Code", options.remove, steps);
    if (!claude.empty() || fs::is_directory(home / L".claude", ignored)) {
        try { claudeHook(home, options.remove, steps); }
        catch (const std::exception& error) { steps.push_back({"Claude Code", "start-up hook not changed", error.what(), false}); }
    }
    const auto codex = options.useAgentCommands ? codexCommand() : fs::path{};
    if (!codex.empty()) registerMcp(codex, {}, "Codex", options.remove, steps);
    else if (fs::is_directory(home / L".codex", ignored)) codexToml(home, options.remove, steps);
    if (steps.empty()) steps.push_back({"Agents", "none found", "install Claude Code or Codex, then run dynsyslab-cli agent-setup"});
    const auto workspace = !options.workspace.empty() ? options.workspace : installedProjectRoot();
    if (!workspace.empty()) workspaceFiles(core::prepareOutputDirectory(workspace, "--workspace"), options.remove, steps);
    bool ok = true;
    if (options.format == "json") {
        json::Value::Array items;
        for (const auto& step : steps) {
            items.push_back(json::Value::Object{{"agent", step.agent}, {"action", step.action}, {"detail", step.detail}, {"ok", step.ok}});
            ok = ok && step.ok;
        }
        std::cout << json::dump(json::Value::Object{{"schema", "dynsyslab.agent-setup/v1"}, {"removed", options.remove},
            {"steps", json::Value(std::move(items))}}) << '\n';
    } else {
        std::cout << (options.remove ? "DynSysLab agent connection removed:\n" : "DynSysLab agent connection:\n");
        for (const auto& step : steps) { std::cout << "- " << step.agent << ": " << step.action << " (" << step.detail << ")\n"; ok = ok && step.ok; }
        if (!options.remove)
            std::cout << "New Claude Code and Codex sessions open with the DynSysLab Hub. Undo with: dynsyslab-cli agent-setup --remove\n";
    }
    return ok ? 0 : 1;
}

}  // namespace dynsyslab::cli
