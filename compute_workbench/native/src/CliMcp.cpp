#include "CliMcp.hpp"
#include "CliProcess.hpp"
#include "CliResources.hpp"
#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/Version.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace dynsyslab::cli {
namespace {
namespace fs = std::filesystem;
namespace json = core::json;
using Object = json::Value::Object;

const char* kInstructions =
    "DynSysLab is a native suite, installed on this computer, for simulating, analysing, modelling (Koopman "
    "reduced-order models) and controlling dynamical systems. At the start of a conversation about DynSysLab, or "
    "in the user's DynSysLab workspace folder, call the `hub` tool first, show the user its greeting (including "
    "the overview when it has one) and its numbered options in the chat, and ask what they would like to do. "
    "Then use the other tools for their requests. Every run generates its own data in the user's DynSysLab "
    "archive folder. Start the desktop applications only when the user asks.";

json::Value schema(Object properties = {}, json::Value::Array required = {}) {
    return Object{{"type", "object"}, {"properties", json::Value(std::move(properties))}, {"required", json::Value(std::move(required))}};
}
json::Value text(const std::string& description) { return Object{{"type", "string"}, {"description", description}}; }

json::Value tool(const char* name, const char* description, json::Value input) {
    return Object{{"name", name}, {"description", description}, {"inputSchema", std::move(input)}};
}

json::Value tools() {
    const auto job = [] { return schema({{"job_path", text("Path of a dynsyslab.compute-job/v1 JSON file")},
        {"job", Object{{"type", "object"}, {"description", "The job itself, instead of job_path"}}},
        {"archive_root", text("Archive folder (default: the installed archive folder)")}}); };
    return json::Value::Array{
        tool("hub", "Start here: the DynSysLab Hub greeting. On the first login after an install it includes a "
             "welcome and an overview of the suite; it always ends with numbered options. Show it to the user and "
             "ask what they would like to do.", schema()),
        tool("quickstart", "Run the bundled examples (Duffing double well, Lorenz 1963): results table, archive paths "
             "and figure paths, as Markdown.", schema()),
        tool("list_systems", "List the dynamical systems DynSysLab can simulate, with presets and parameters.", schema()),
        tool("doctor", "Check the installation and report its folders.", schema()),
        tool("validate_job", "Validate a compute job.", job()),
        tool("plan_job", "Show the stages a compute job or pipeline would run.", job()),
        tool("run_job", "Run a compute job or pipeline; returns its archives as JSON.", job()),
        tool("fourier", "Fourier spectrum of one state of a trajectory archive.",
             schema({{"archive", text("Trajectory .npz archive")}, {"component", Object{{"type", "integer"}}}}, {"archive"})),
        tool("render_figure", "Render a Viewer figure of an archive to PNG and return the image.",
             schema({{"archive", text("Archive .npz in the archive folder")},
                     {"view", text("Viewer route: time, phase, fourier, control_fourier, estimator, ...")},
                     {"output", text("PNG path (default: next to the archive folder in viewer_outputs)")}}, {"archive", "view"})),
    };
}

std::wstring w(const std::string& value) { return core::wide(value); }

std::string base64(const std::string& bytes) {
    static const char* table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out; std::size_t i = 0;
    for (; i + 2 < bytes.size(); i += 3) {
        const unsigned v = (unsigned char)bytes[i] << 16 | (unsigned char)bytes[i + 1] << 8 | (unsigned char)bytes[i + 2];
        out += {table[v >> 18], table[v >> 12 & 63], table[v >> 6 & 63], table[v & 63]};
    }
    if (i + 1 == bytes.size()) { const unsigned v = (unsigned char)bytes[i] << 16; out += {table[v >> 18], table[v >> 12 & 63], '=', '='}; }
    else if (i + 2 == bytes.size()) { const unsigned v = (unsigned char)bytes[i] << 16 | (unsigned char)bytes[i + 1] << 8;
        out += {table[v >> 18], table[v >> 12 & 63], table[v >> 6 & 63], '='}; }
    return out;
}

json::Value result(const std::string& body, bool error = false) {
    return Object{{"content", json::Value::Array{Object{{"type", "text"}, {"text", body}}}}, {"isError", error}};
}

json::Value cli(std::vector<std::wstring> args) {
    const auto run = runProcess(cliExecutablePath(), args);
    if (!run.started) return result("unable to start dynsyslab-cli", true);
    return result(run.exitCode == 0 ? run.out : run.out + run.err, run.exitCode != 0);
}

json::Value jobCommand(const char* command, const json::Value& args) {
    std::wstring path = w(args.at("job_path").string());
    if (path.empty()) {
        if (args.at("job").isNull()) return result("give job_path or job", true);
        const auto file = fs::temp_directory_path() / (L"dynsyslab-mcp-job-" +
            std::to_wstring(std::chrono::steady_clock::now().time_since_epoch().count()) + L".json");
        std::ofstream(file, std::ios::binary) << json::dump(args.at("job"));
        path = file.wstring();
    }
    std::vector<std::wstring> argv{w(command), path};
    if (const auto root = args.at("archive_root").string(); !root.empty()) argv.insert(argv.end(), {L"--archive-root", w(root)});
    return cli(argv);
}

json::Value renderFigure(const json::Value& args) {
    const auto archive = fs::u8path(args.at("archive").string());
    auto root = defaultArchiveRoot().path;
    auto output = fs::u8path(args.at("output").string());
    if (output.empty()) output = root.parent_path() / L"viewer_outputs" / L"mcp" /
        (archive.stem().wstring() + L"-" + w(args.at("view").string()) + L".png");
    const auto run = runProcess(locateResources().viewer, {L"render", L"--archives", root.wstring(), L"--source",
        archive.wstring(), L"--view", w(args.at("view").string()), L"--width", L"1200", L"--height", L"800", L"--output", output.wstring()});
    if (!run.started || run.exitCode != 0) return result("figure not rendered: " + run.err, true);
    std::ifstream png(output, std::ios::binary);
    const std::string bytes{std::istreambuf_iterator<char>(png), std::istreambuf_iterator<char>()};
    return Object{{"content", json::Value::Array{Object{{"type", "image"}, {"data", base64(bytes)}, {"mimeType", "image/png"}},
        Object{{"type", "text"}, {"text", "Saved to " + shownPath(output)}}}}, {"isError", false}};
}

json::Value callTool(const std::string& name, const json::Value& args) {
    if (name == "hub") return cli({L"hub"});
    if (name == "quickstart") return cli({L"quickstart"});
    if (name == "list_systems") return cli({L"systems"});
    if (name == "doctor") return cli({L"doctor"});
    if (name == "validate_job") return jobCommand("validate", args);
    if (name == "plan_job") return jobCommand("plan", args);
    if (name == "run_job") return jobCommand("run", args);
    if (name == "fourier") return cli({L"fourier", w(args.at("archive").string()),
        L"--component", std::to_wstring(static_cast<long long>(args.at("component").number(0)))});
    if (name == "render_figure") return renderFigure(args);
    return result("unknown tool: " + name, true);
}

json::Value handle(const std::string& method, const json::Value& params) {
    if (method == "initialize") {
        const auto requested = params.at("protocolVersion").string();
        const bool known = requested == "2024-11-05" || requested == "2025-03-26" || requested == "2025-06-18";
        return Object{{"protocolVersion", known ? requested : std::string("2025-06-18")},
            {"capabilities", Object{{"tools", Object{}}, {"prompts", Object{}}}},
            {"serverInfo", Object{{"name", "dynsyslab"}, {"version", std::string(core::version)}}},
            {"instructions", kInstructions}};
    }
    if (method == "ping") return Object{};
    if (method == "tools/list") return Object{{"tools", tools()}};
    if (method == "tools/call") return callTool(params.at("name").string(), params.at("arguments"));
    if (method == "prompts/list")
        return Object{{"prompts", json::Value::Array{Object{{"name", "hub"},
            {"description", "Open the DynSysLab Hub: greeting, first-login overview and what you can do"}}}}};
    if (method == "prompts/get") {
        const auto run = runProcess(cliExecutablePath(), {L"hub"});
        return Object{{"messages", json::Value::Array{Object{{"role", "user"}, {"content", Object{{"type", "text"},
            {"text", "Open the DynSysLab Hub: greet me with this, keep the numbered options, and ask what I would "
                     "like to do.\n\n" + run.out}}}}}}};
    }
    throw std::invalid_argument("method not found: " + method);
}

}  // namespace

int mcpServer() {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY); _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string line;
    while (std::getline(std::cin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("\xEF\xBB\xBF", 0) == 0) line.erase(0, 3);  // a byte-order mark (PowerShell pipes)
        if (line.find_first_not_of(" \t") == std::string::npos) continue;
        json::Value request;
        try { request = json::parse(line); }
        catch (const std::exception& error) {
            std::cout << json::dump(Object{{"jsonrpc", "2.0"}, {"id", json::Value()},
                {"error", Object{{"code", -32700.0}, {"message", error.what()}}}}) << '\n' << std::flush;
            continue;
        }
        const auto& id = request.at("id");
        if (id.isNull()) continue;  // notifications (notifications/initialized, cancelled) need no reply
        Object reply{{"jsonrpc", "2.0"}, {"id", id}};
        try { reply["result"] = handle(request.at("method").string(), request.at("params")); }
        catch (const std::invalid_argument& error) {
            reply["error"] = Object{{"code", -32601.0}, {"message", error.what()}};
        } catch (const std::exception& error) {
            reply["error"] = Object{{"code", -32603.0}, {"message", error.what()}};
        }
        std::cout << json::dump(json::Value(std::move(reply))) << '\n' << std::flush;
    }
    return 0;
}

}  // namespace dynsyslab::cli
