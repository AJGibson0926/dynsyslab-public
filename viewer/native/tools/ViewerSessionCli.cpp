#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
struct Options {
    std::map<std::wstring,std::wstring> values;
    bool samples = false;
    std::wstring get(const wchar_t* key) const {
        const auto found = values.find(key);
        return found == values.end() ? L"" : found->second;
    }
    std::string text(const wchar_t* key) const { return dsl::narrow(get(key)); }
    int number(const wchar_t* key,int fallback) const {
        const auto value = get(key); if (value.empty()) return fallback;
        std::size_t end = 0; const int result = std::stoi(value,&end);
        if (end != value.size() || result < 0) throw std::runtime_error("Invalid nonnegative integer: " + dsl::narrow(value));
        return result;
    }
};
Options parse(int argc,wchar_t** argv) {
    const std::set<std::wstring> keys{L"--project",L"--archives",L"--cache",L"--catalog",L"--focus",
        L"--id",L"--source",L"--view",L"--coordinates",L"--search",L"--parent",L"--output",
        L"--mode",L"--options",L"--width",L"--height",L"--page",L"--figure",L"--workspace",L"--compare-id",L"--modes",L"--complex-display",L"--scale",L"--scroll-y",L"--benchmark-scroll"};
    Options result;
    for (int i = 2; i < argc; ++i) {
        const std::wstring key = argv[i];
        if (key == L"--json") continue;
        if (key == L"--samples") { result.samples = true; continue; }
        if (!keys.count(key)) throw std::runtime_error("Unknown option: " + dsl::narrow(key));
        if (++i >= argc) throw std::runtime_error("Missing value: " + dsl::narrow(key));
        if (!result.values.emplace(key,argv[i]).second) throw std::runtime_error("Duplicate option: " + dsl::narrow(key));
    }
    return result;
}
}

void viewerSessionUsage() {
    std::cout
        << "  tree [--search text] [--parent id]\n"
        << "  inspect (--id id | --source path) [--view name] [--coordinates dimensional|dimensionless] [--samples]\n"
        << "  panels (--id id | --source path)\n"
        << "  phase-space (--id id | --source path) [--mode lookup|compute] [--options options.json]\n"
        << "  render (--id id | --source path) --output plot.png [--view name] [--coordinates name]\n"
        << "         [--width 1100 --height 760] [--page 0 | --figure 0]\n"
        << "  Optional render profiling: --benchmark-scroll 60 (invisible redraws)\n"
        << "  Values: --scale auto|linear|log|symlog; scrolling workspaces: --scroll-y pixels\n"
        << "  Workspace: --workspace main|observables [--compare-id id] (inspect/render)\n"
        << "  KROM: --modes 0,1|all|none --complex-display real-imaginary|magnitude-phase|domain-coloring\n"
        << "  Common: --project path --archives path --cache path --catalog path --focus path\n"
        << "  JSON output by default; tree/inspect/panels accept --output file.json.\n"
        << "  The full catalog is used by default; --focus explicitly applies a temporary Mask.\n"
        << "  All commands run without visible windows.\n";
}

int viewerSessionCommand(int argc,wchar_t** argv) {
    try {
        const auto args = parse(argc,argv);
        const auto paths = dsl::discoverRuntimePaths(L"Viewer",args.get(L"--project"),
            args.get(L"--archives"),args.get(L"--cache"));
        dsl::BackendConfig config;
        config.projectRoot = paths.projectRoot; config.archiveRoot = paths.archiveRoot;
        config.cache = paths.cacheRoot; config.catalog = args.get(L"--catalog");
        config.archiveFocus = args.get(L"--focus");
        dsl::ViewerSession session(config);
        const std::wstring command = argv[1];
        dsl::json::Value result;
        if (command == L"tree") {
            std::optional<std::string> parent;
            if (args.values.count(L"--parent")) parent = args.text(L"--parent");
            result = session.tree(args.text(L"--search"),parent);
        } else {
            if (!args.get(L"--id").empty() && !args.get(L"--source").empty())
                throw std::runtime_error("Select with either --id or --source");
            if (!args.get(L"--id").empty()) session.select(args.text(L"--id"));
            else if (!args.get(L"--source").empty()) session.selectSource(args.get(L"--source"));
            else throw std::runtime_error("Select a data-tree item with --id or --source");
            if (command == L"phase-space") {
                const auto options = args.get(L"--options").empty() ? dsl::json::Value{} : dsl::json::load(args.get(L"--options"));
                result = session.phaseSpace(args.get(L"--mode").empty() ? "lookup" : args.text(L"--mode"), options);
            } else if (command == L"panels") {
                dsl::json::Value::Array panels;
                for (const auto& panel : session.panels())
                    panels.emplace_back(dsl::json::Value::Object{{"route",panel.route},{"family",panel.family},
                        {"label",panel.label},{"detail",panel.detail}});
                result = dsl::json::Value::Object{{"schema","dynsyslab.viewer.panels/v1"},
                    {"selected_id",session.selected().id},{"panels",std::move(panels)}};
            } else {
                const auto workspace = args.get(L"--workspace").empty() ? "main" : args.text(L"--workspace");
                dsl::KromWorkspaceOptions krom;
                if (args.values.count(L"--complex-display")) krom.complexDisplay = args.text(L"--complex-display");
                if (args.values.count(L"--modes")) {
                    const auto text = args.text(L"--modes");
                    krom.allModes = text == "all";
                    if (!krom.allModes) {
                        krom.modes = std::vector<std::size_t>{};
                        if (text != "none") {
                            if (text.empty() || text.back() == ',') throw std::runtime_error("Invalid KROM mode list");
                            std::stringstream stream(text); std::string token;
                            while (std::getline(stream,token,',')) {
                                if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos)
                                    throw std::runtime_error("Modes must be retained indices separated by commas");
                                krom.modes->push_back(std::stoull(token));
                            }
                        }
                    }
                }
                if (session.selected().kind != "krom" && (args.values.count(L"--modes") || args.values.count(L"--complex-display")))
                    throw std::runtime_error("KROM options require a KROM data-tree item");
                auto payload = session.workspace(workspace, args.text(L"--view"),
                    args.text(L"--coordinates"), args.text(L"--compare-id"),krom);
                if (args.values.count(L"--scale")) dsl::setValueScale(payload,args.text(L"--scale"));
                if (command == L"render") {
                    if (args.get(L"--output").empty()) throw std::runtime_error("render requires --output plot.png");
                    dsl::PlotRenderOptions options;
                    options.width = args.number(L"--width",1100); options.height = args.number(L"--height",760);
                    options.page = args.number(L"--page",0);
                    options.scrollY = args.number(L"--scroll-y",0);
                    options.benchmarkScrollSteps = args.number(L"--benchmark-scroll",0);
                    if (args.values.count(L"--figure")) {
                        if (args.values.count(L"--page")) throw std::runtime_error("Use --page or --figure, not both");
                        options.figure = static_cast<std::size_t>(args.number(L"--figure",0));
                    }
                    result = dsl::renderPlotPng(payload,args.get(L"--output"),options);
                } else result = dsl::inspectPlot(payload,args.samples);
                result.object()["workspace"] = workspace;
                result.object()["selected_id"] = session.selected().id;
                result.object()["focus_label"] = session.catalog().focusLabel;
            }
        }
        const auto encoded = dsl::json::dump(result);
        if (command != L"render" && !args.get(L"--output").empty()) {
            std::ofstream output(std::filesystem::path(args.get(L"--output")),std::ios::binary);
            output << encoded << '\n';
            if (!output) throw std::runtime_error("Unable to write JSON output");
        } else std::cout << encoded << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << dsl::json::dump(dsl::json::Value::Object{
            {"schema","dynsyslab.viewer.error/v1"},{"error",error.what()}}) << '\n';
        return 1;
    }
}
