#include "CliHub.hpp"
#include "CliResources.hpp"
#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/Version.hpp>

#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

namespace dynsyslab::cli {
namespace {
namespace fs = std::filesystem;
namespace json = core::json;

struct Option { const char* id; const char* title; const char* detail; };
constexpr Option kOptions[]{
    {"quickstart", "Run the quickstart", "simulate a Duffing double well and the chaotic Lorenz system, then show the results and figures"},
    {"simulate", "Simulate a system", "pick a system (Duffing, Lorenz, pendula, bubble dynamics/RPE, Hodgkin-Huxley, ...), set parameters, initial conditions and forcing"},
    {"analyse", "Analyse a run", "Fourier spectra, phase portraits, SVD/POD and comparisons of runs in your archive"},
    {"model", "Build a Koopman model (K-ROM)", "fit DMD, H-DMD, SINDy or spectral Koopman models to simulated data and test their predictions"},
    {"control", "Estimate and control", "design KLQE/KLQR, LQR, MPC or PID, check observability and loop margins, simulate the closed loop"},
    {"apps", "Open the desktop apps", "the Solver/Post-processor and the Viewer"},
    {"other", "Something else", "describe what you want to do"}};

double now() { return static_cast<double>(std::time(nullptr)); }

struct State { std::string welcomedVersion; double firstLogin{}, lastLogin{}; double logins{}; };

fs::path statePath() { return userDataRoot() / L"hub.json"; }

State loadState() {
    State state;
    try {
        const auto value = json::load(statePath().wstring());
        state.welcomedVersion = value.at("welcomed_version").string();
        state.firstLogin = value.at("first_login_unix").number(0);
        state.lastLogin = value.at("last_login_unix").number(0);
        state.logins = value.at("logins").number(0);
    } catch (...) {}
    return state;
}

void saveState(const State& state) {
    const auto path = core::prepareOutputFile(statePath(), "hub state");
    std::ofstream(path, std::ios::binary) << json::dump(json::Value::Object{{"schema", "dynsyslab.hub-state/v1"},
        {"welcomed_version", state.welcomedVersion}, {"first_login_unix", state.firstLogin},
        {"last_login_unix", state.lastLogin}, {"logins", state.logins}}) << '\n';
}

// First login after an install; it stays "first" for 30 minutes so a hook and an MCP call in the
// same session agree.
bool firstLogin(const State& state) {
    return state.welcomedVersion != core::version || now() - state.firstLogin < 1800.0;
}

std::size_t runCount(const fs::path& archiveRoot) {
    std::size_t count = 0; std::error_code error;
    for (fs::recursive_directory_iterator it(archiveRoot / L"solutions", error), end; !error && it != end; it.increment(error))
        if (it->path().extension() == L".npz") ++count;
    return count;
}

bool inWorkspace(fs::path folder) {
    std::error_code ignored;
    for (; !folder.empty(); folder = folder.parent_path()) {
        if (fs::is_regular_file(folder / kWorkspaceMarker, ignored)) return true;
        if (folder == folder.parent_path()) break;
    }
    return false;
}

std::string greeting(bool first, const fs::path& archiveRoot) {
    std::ostringstream out;
    const auto runs = runCount(archiveRoot);
    if (first) {
        out << readWelcome(locateResources().welcome).overview << "\n\n## Your DynSysLab\n\n"
            << "- Version " << core::version << ", running on this computer; nothing is sent elsewhere.\n"
            << "- Results are saved in `" << shownPath(archiveRoot) << "`"
            << (runs ? " (" + std::to_string(runs) + " runs so far)." : " (empty: every example generates its own data).") << "\n"
            << "- The desktop apps (Solver/Post-processor, Viewer) are optional; they open only when you ask.\n";
    } else {
        out << "# Welcome back to DynSysLab\n\nDynSysLab " << core::version << ". Results are in `"
            << shownPath(archiveRoot) << "`" << (runs ? " (" + std::to_string(runs) + " runs)." : ".") << "\n";
    }
    out << "\n## What would you like to do?\n\n";
    int n = 0;
    for (const auto& option : kOptions) out << ++n << ". **" << option.title << "**: " << option.detail << ".\n";
    return out.str();
}

std::string hookCwd() {
    try {
        const std::string input{std::istreambuf_iterator<char>(std::cin), std::istreambuf_iterator<char>()};
        if (!input.empty()) return json::parse(input).at("cwd").string();
    } catch (...) {}
    return {};
}

}  // namespace

int hub(const HubOptions& options) {
    if (options.format != "markdown" && options.format != "json")
        throw std::invalid_argument("hub --format must be markdown or json");
    auto state = loadState();
    const bool first = firstLogin(state);
    if (options.hook) {
        const auto cwd = hookCwd();
        const bool workspace = inWorkspace(cwd.empty() ? fs::current_path() : fs::u8path(cwd));
        if (!first && !workspace) return 0;  // unrelated project after the first login: stay quiet
    }
    if (!options.peek) {
        if (state.welcomedVersion != core::version) { state.welcomedVersion = core::version; state.firstLogin = now(); }
        state.lastLogin = now(); state.logins += 1;
        try { saveState(state); } catch (...) {}  // the greeting never fails on a read-only profile
    }
    const auto text = asciiText(greeting(first, options.archiveRoot));
    if (options.hook) {
        const std::string context = "DynSysLab Hub (start-up hook installed with DynSysLab): begin your first reply "
            "by greeting the user with the following, keeping its overview and numbered options, then ask what they "
            "would like to do. Use the dynsyslab MCP tools or dynsyslab-cli for the work, and start the desktop "
            "applications only when the user asks.\n\n" + text;
        std::cout << json::dump(json::Value::Object{{"hookSpecificOutput", json::Value::Object{
            {"hookEventName", "SessionStart"}, {"additionalContext", context}}}}) << '\n';
        return 0;
    }
    if (options.format == "json") {
        json::Value::Array items;
        for (const auto& option : kOptions)
            items.push_back(json::Value::Object{{"id", option.id}, {"title", option.title}, {"detail", option.detail}});
        std::cout << json::dump(json::Value::Object{{"schema", "dynsyslab.hub/v1"}, {"version", std::string(core::version)},
            {"first_login", first}, {"greeting_markdown", text}, {"options", json::Value(std::move(items))},
            {"archive_root", shownPath(options.archiveRoot)}}) << '\n';
    } else {
        std::cout << text;
    }
    return 0;
}

}  // namespace dynsyslab::cli
