#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/TimeSeriesOverlay.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"
#include <dynsyslab/core/OutputPath.hpp>
#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
struct Options {
    dsl::TimeSeriesOverlay overlay;
    std::map<std::wstring,std::wstring> values;
    // Repeatable archive-view options, in order: further runs, then one label/color per run.
    std::vector<std::string> compareIds, labels, colors;
    std::vector<std::wstring> compareSources;
    bool samples = false;
    std::wstring get(const wchar_t* key) const {
        const auto found = values.find(key);
        return found == values.end() ? L"" : found->second;
    }
    std::string text(const wchar_t* key) const { return dsl::narrow(get(key)); }
    std::optional<std::pair<double,double>> range(const wchar_t* key) const {
        const auto value = get(key); if (value.empty()) return {};
        const auto colon = value.find(L':');
        if (colon == std::wstring::npos) throw std::runtime_error("Expected a range low:high, got " + dsl::narrow(value));
        std::size_t a = 0, b = 0;
        const auto low = std::stod(value.substr(0, colon), &a), high = std::stod(value.substr(colon + 1), &b);
        if (a != colon || b != value.size() - colon - 1 || !(high > low))
            throw std::runtime_error("Expected a range low:high with high > low, got " + dsl::narrow(value));
        return std::make_pair(low, high);
    }
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
        L"--mode",L"--options",L"--width",L"--height",L"--page",L"--figure",L"--workspace",L"--modes",L"--complex-display",L"--scale",L"--scroll-y",L"--benchmark-scroll",
        L"--panel",L"--start",L"--stop",L"--title",L"--x-range",L"--y-range"};
    Options result;
    const bool series=std::wstring(argv[1])==L"time-series";
    const std::set<std::wstring> seriesKeys{L"--source",L"--output",L"--width",L"--height",L"--page",
        L"--title",L"--time-key",L"--time-label",L"--start",L"--stop",L"--columns",L"--value-label",L"--font-scale",L"--x-label",
        L"--frequency-scale",L"--scale",L"--x-range",L"--y-range",L"--stem-baseline"};
    for (int i = 2; i < argc; ++i) {
        const std::wstring key = argv[i];
        if (series && key==L"--line-width") {
            if (result.overlay.curves.empty() || ++i>=argc)
                throw std::runtime_error("Use --line-width after a --curve");
            const std::wstring text=argv[i]; std::size_t used=0;
            const double width=std::stod(text,&used);
            if (used!=text.size() || !std::isfinite(width) || width<=0 || width>32)
                throw std::runtime_error("Expected a positive line width in points up to 32");
            result.overlay.curves.back().lineWidthPoints=width;
            continue;
        }
        if (series && key==L"--line-style") {
            if (result.overlay.curves.empty() || ++i>=argc)
                throw std::runtime_error("Use --line-style after a --curve");
            const auto style=dsl::narrow(argv[i]);
            if (style!="solid" && style!="dash" && style!="time-dash" && style!="sparse-dash"
                && style!="dot" && style!="dash-dot" && style!="stem")
                throw std::runtime_error("Expected solid, dash, time-dash, sparse-dash, dot, dash-dot, or stem");
            result.overlay.curves.back().lineStyle=style;
            continue;
        }
        if (series && key==L"--color") {
            if (result.overlay.curves.empty() || ++i>=argc)
                throw std::runtime_error("Use --color #RRGGBB after a --curve");
            // One colour, or a comma-separated colour per component (panel).
            const auto list=dsl::narrow(argv[i]);
            std::vector<std::string> colors;
            for (std::size_t from=0;from<=list.size();) {
                const auto to=list.find(',',from);
                colors.push_back(list.substr(from,to==std::string::npos?std::string::npos:to-from));
                if (to==std::string::npos) break; from=to+1;
            }
            for (const auto& color:colors)
                if (color.size()!=7 || color[0]!='#' || color.find_first_not_of("0123456789abcdefABCDEF",1)!=std::string::npos)
                    throw std::runtime_error("Expected hex curve colors such as #64189B or #D8242A,#004DB8");
            result.overlay.curves.back().color=colors.front();
            if (colors.size()>1) result.overlay.curves.back().componentColors=colors;
            continue;
        }
        if (series && (key==L"--curve" || key==L"--axis")) {
            if (++i>=argc) throw std::runtime_error("Missing value for "+dsl::narrow(key));
            const auto text=dsl::narrow(argv[i]);
            if (key==L"--axis") result.overlay.componentLabels.push_back(text);
            else {
                const auto split=text.find('=');
                const auto name=text.substr(0,split);
                if (name.empty()) throw std::runtime_error("Empty curve key");
                result.overlay.curves.push_back({name,split==std::string::npos?name:text.substr(split+1),{},
                    result.overlay.curves.size()%2==1?"dash":"solid"});
            }
            continue;
        }
        if (!series && (key==L"--compare-id" || key==L"--compare-source" || key==L"--label" || key==L"--color")) {
            if (++i>=argc) throw std::runtime_error("Missing value: "+dsl::narrow(key));
            if (key==L"--compare-source") { result.compareSources.push_back(argv[i]); continue; }
            const auto text=dsl::narrow(argv[i]);
            if (key==L"--color" && (text.size()!=7 || text[0]!='#'
                || text.find_first_not_of("0123456789abcdefABCDEF",1)!=std::string::npos))
                throw std::runtime_error("Expected a hex run color such as #004DB8");
            (key==L"--compare-id"?result.compareIds:key==L"--label"?result.labels:result.colors).push_back(text);
            continue;
        }
        if (series && key==L"--spectrum") { result.overlay.spectrum=true; continue; }
        if (series && key==L"--integer-x-ticks") { result.overlay.integerTicks=true; continue; }
        if (key == L"--json") continue;
        if (key == L"--samples") { result.samples = true; continue; }
        if (!(series?seriesKeys:keys).count(key)) throw std::runtime_error("Unknown option: " + dsl::narrow(key));
        if (++i >= argc) throw std::runtime_error("Missing value: " + dsl::narrow(key));
        if (!result.values.emplace(key,argv[i]).second) throw std::runtime_error("Duplicate option: " + dsl::narrow(key));
    }
    return result;
}
}

void viewerSessionUsage() {
    std::cout
        << "  Archive views (see docs/architecture/VIEWER_CLI.md):\n"
        << "  panels (--id id | --source path)              routes and panel labels of one archive\n"
        << "  render (--id id | --source path) --output plot.png [--view route|label] [--panel label]\n"
        << "         [--start t --stop t] [--coordinates dimensional|dimensionless] [--scale auto|linear|log|symlog]\n"
        << "         [--compare-source path | --compare-id id ...] [--label text ...] [--color #RRGGBB ...]\n"
        << "         [--title text] [--x-range low:high] [--y-range low:high] [--width 1100 --height 760] [--page 0 | --figure 0]\n"
        << "  inspect: same selection as render, JSON instead of PNG [--samples]\n"
        << "      Views: time, phase, fourier, control_fourier, estimator, ... or a label from panels.\n"
        << "      --start/--stop: archive clock units; windows every route (spectra too).\n"
        << "      Comparisons overlay runs in order (selected first); --label/--color: one per run.\n"
        << "  time-series --source data.npz --curve key|xkey,ykey[=label] [--color #RRGGBB[,#RRGGBB per panel]] [--line-width points] [--line-style solid|dash|time-dash|sparse-dash|dot|dash-dot|stem] ... --output plot.png\n"
        << "      [--axis component-label ...] [--title text] [--time-key time] [--time-label Time] [--start t --stop t] [--columns n] [--value-label text] [--x-label text] [--font-scale x] [--width n --height n --page n]\n"
        << "      [--spectrum [--frequency-scale k]] one-sided Hann amplitude spectra of the curves over [start, stop], frequency x k\n"
        << "      [--scale auto|linear|log|symlog] [--x-range low:high] [--y-range low:high]\n"
        << "      [--stem-baseline v] baseline of stem curves (default 0; a decade below the smallest value with --scale log)\n"
        << "      [--integer-x-ticks] ticks and grid at integer abscissas (an index axis)\n"
        << "  tree [--search text] [--parent id]\n"
        << "  phase-space (--id id | --source path) [--mode lookup|compute] [--options options.json]\n"
        << "  Workspace: --workspace main|observables (inspect/render)\n"
        << "  KROM: --modes 0,1|all|none --complex-display real-imaginary|magnitude-phase|domain-coloring\n"
        << "  Optional render profiling: --benchmark-scroll 60 (invisible redraws); scrolling workspaces: --scroll-y pixels\n"
        << "  Common: --project path --archives path --cache path --catalog path --focus path\n"
        << "  JSON output by default; tree/inspect/panels accept --output file.json.\n"
        << "  The full catalog is used by default; --focus explicitly applies a temporary Mask.\n"
        << "  All commands run without visible windows.\n";
}

int viewerSessionCommand(int argc,wchar_t** argv) {
    try {
        const auto args = parse(argc,argv);
        if (std::wstring(argv[1])==L"time-series") {
            auto request=args.overlay;
            const auto source=args.get(L"--source"),output=args.get(L"--output");
            if (source.empty() || output.empty()) throw std::runtime_error("time-series requires --source and --output");
            if (args.values.count(L"--time-key")) request.timeKey=args.text(L"--time-key");
            if (args.values.count(L"--time-label")) request.timeLabel=args.text(L"--time-label");
            request.title=args.text(L"--title");
            const auto limit=[&](const wchar_t* flag)->std::optional<double> {
                if (!args.values.count(flag)) return {};
                const auto text=args.get(flag); std::size_t used=0;
                const auto value=std::stod(text,&used);
                if (used!=text.size()) throw std::invalid_argument("Invalid time limit: "+dsl::narrow(text));
                return value;
            };
            request.start=limit(L"--start"); request.stop=limit(L"--stop");
            if (args.values.count(L"--value-label")) request.valueLabel=args.text(L"--value-label");
            if (args.values.count(L"--x-label")) request.xLabel=args.text(L"--x-label");
            if (args.values.count(L"--columns")) request.columns=static_cast<std::size_t>(args.number(L"--columns",0));
            if (args.values.count(L"--frequency-scale")) {
                if (!request.spectrum) throw std::runtime_error("--frequency-scale requires --spectrum");
                const auto text=args.get(L"--frequency-scale"); std::size_t used=0;
                request.frequencyScale=std::stod(text,&used);
                if (used!=text.size() || !(request.frequencyScale>0)) throw std::invalid_argument("Invalid frequency scale: "+dsl::narrow(text));
            }
            if (args.values.count(L"--stem-baseline")) {
                const auto text=args.get(L"--stem-baseline"); std::size_t used=0;
                request.stemBaseline=std::stod(text,&used);
                if (used!=text.size() || !std::isfinite(*request.stemBaseline)) throw std::invalid_argument("Invalid stem baseline: "+dsl::narrow(text));
            }
            request.logBaseline=args.values.count(L"--scale") && args.text(L"--scale")=="log";
            dsl::PlotRenderOptions options;
            options.width=args.number(L"--width",1100); options.height=args.number(L"--height",760);
            options.page=args.number(L"--page",0); options.toolbar=false;
            if (args.values.count(L"--font-scale")) {
                const auto text=args.get(L"--font-scale"); std::size_t used=0;
                options.fontScale=std::stod(text,&used);
                if (used!=text.size()) throw std::invalid_argument("Invalid font scale: "+dsl::narrow(text));
            }
            options.xRange=args.range(L"--x-range"); options.yRange=args.range(L"--y-range");
            auto payload=dsl::timeSeriesOverlay(source,request,4ULL*options.width);
            if (args.values.count(L"--scale")) dsl::setValueScale(payload,args.text(L"--scale"));
            std::cout << dsl::json::dump(dsl::renderPlotPng(payload,
                dynsyslab::core::prepareOutputFile(output,"--output"),options)) << '\n';
            return 0;
        }
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
                dsl::ViewRequest request;
                request.workspace = args.get(L"--workspace").empty() ? "main" : args.text(L"--workspace");
                const auto& workspace = request.workspace;
                request.view = args.text(L"--view"); request.panel = args.text(L"--panel");
                request.coordinates = args.text(L"--coordinates");
                request.compareIds = args.compareIds; request.labels = args.labels; request.colors = args.colors;
                for (const auto& source : args.compareSources) request.compareSources.emplace_back(source);
                const auto limit=[&](const wchar_t* flag)->std::optional<double> {
                    if (!args.values.count(flag)) return {};
                    const auto text=args.get(flag); std::size_t used=0;
                    const auto value=std::stod(text,&used);
                    if (used!=text.size()) throw std::invalid_argument("Invalid time limit: "+dsl::narrow(text));
                    return value;
                };
                request.window.start = limit(L"--start"); request.window.stop = limit(L"--stop");
                auto& krom = request.krom;
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
                auto payload = session.workspace(request);
                if (args.values.count(L"--title"))
                    for (auto& figure : payload.figures)
                        figure.title = payload.figures.size() == 1 ? args.text(L"--title")
                                                                     : args.text(L"--title") + " | " + figure.title;
                if (args.values.count(L"--scale")) dsl::setValueScale(payload,args.text(L"--scale"));
                if (command == L"render") {
                    if (args.get(L"--output").empty()) throw std::runtime_error("render requires --output plot.png");
                    dsl::PlotRenderOptions options;
                    options.width = args.number(L"--width",1100); options.height = args.number(L"--height",760);
                    options.page = args.number(L"--page",0);
                    options.scrollY = args.number(L"--scroll-y",0);
                    options.benchmarkScrollSteps = args.number(L"--benchmark-scroll",0);
                    options.xRange = args.range(L"--x-range"); options.yRange = args.range(L"--y-range");
                    if (args.values.count(L"--figure")) {
                        if (args.values.count(L"--page")) throw std::runtime_error("Use --page or --figure, not both");
                        options.figure = static_cast<std::size_t>(args.number(L"--figure",0));
                    }
                    const auto output = args.get(L"--output");
                    result = dsl::renderPlotPng(payload, output.empty() ? std::filesystem::path(output)
                        : dynsyslab::core::prepareOutputFile(output,"--output"), options);
                } else result = dsl::inspectPlot(payload,args.samples);
                result.object()["workspace"] = workspace;
                result.object()["selected_id"] = session.selected().id;
                result.object()["focus_label"] = session.catalog().focusLabel;
            }
        }
        const auto encoded = dsl::json::dump(result);
        if (command != L"render" && !args.get(L"--output").empty()) {
            std::ofstream output(dynsyslab::core::prepareOutputFile(args.get(L"--output"),"--output"),std::ios::binary);
            output << encoded << '\n';
            if (!output) throw std::runtime_error("Unable to write JSON output");
        } else std::cout << encoded << '\n';
        return 0;
    } catch (const dynsyslab::core::OutputPathError& error) {
        std::cerr << dsl::json::dump(dsl::json::Value::Object{{"schema","dynsyslab.viewer.error/v1"},
            {"error",error.what()},{"option",error.option()},{"path",error.path().u8string()}}) << '\n';
        return 1;
    } catch (const std::exception& error) {
        std::cerr << dsl::json::dump(dsl::json::Value::Object{
            {"schema","dynsyslab.viewer.error/v1"},{"error",error.what()}}) << '\n';
        return 1;
    }
}
