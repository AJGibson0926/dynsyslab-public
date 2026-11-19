#include "CliQuickstart.hpp"
#include "CliRuntimePaths.hpp"
#include "CliViewerFigures.hpp"

#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>
#include <dynsyslab/core/Version.hpp>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace dynsyslab::cli {
namespace {
namespace fs = std::filesystem;
namespace json = core::json;

struct Resources { fs::path examples, welcome, viewer; };

Resources locate() {
    const auto exe = cliExecutablePath();
    if (const auto checkout = sourceCheckoutContaining(exe); !checkout.empty())
        return {checkout / L"examples" / L"native" / L"quickstart", checkout / L"docs" / L"AGENT_WELCOME.md",
                exe.parent_path() / L"DynSysLabViewerCli.exe"};
    const auto root = exe.parent_path().parent_path() / L"resources";
    return {root / L"examples" / L"quickstart", root / L"docs" / L"AGENT_WELCOME.md",
            exe.parent_path() / L"DynSysLabViewerCli.exe"};
}

std::string ascii(std::string text) {  // typographic punctuation survives any console code page
    const std::pair<const char*, const char*> folds[]{{"\xE2\x80\x93", "-"}, {"\xE2\x80\x94", "-"},
        {"\xE2\x80\xA6", "..."}, {"\xE2\x86\x92", "->"}, {"\xE2\x80\x98", "'"}, {"\xE2\x80\x99", "'"},
        {"\xE2\x80\x9C", "\""}, {"\xE2\x80\x9D", "\""}, {"\xC3\xA9", "e"}};
    for (const auto& [from, to] : folds)
        for (auto at = text.find(from); at != std::string::npos; at = text.find(from, at + std::strlen(to)))
            text.replace(at, std::strlen(from), to);
    return text;
}

struct Welcome { std::string overview; std::vector<std::string> next; };

// Overview = the document up to "## Get started"; next steps = bullets under "## Things you can ask".
Welcome welcome(const fs::path& path) {
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
        if (asks && line.rfind("- ", 0) == 0) out.next.push_back(ascii(line.substr(2)));
        if (overview && line != "---") out.overview += ascii(line) + "\n";
    }
    while (out.overview.size() > 1 && out.overview.substr(out.overview.size() - 2) == "\n\n") out.overview.pop_back();
    return out;
}

std::string shown(const fs::path& path) { return path.lexically_normal().make_preferred().u8string(); }

std::string number(double value, int digits = 4) {
    std::ostringstream out; out << std::setprecision(digits) << value; return out.str();
}

struct Figure { std::string title, path, error; };
struct Run {
    std::string title, about, system, preset, archive, finalState, fourierState;
    std::size_t samples{}, accepted{}, rejected{};
    double start{}, stop{}, dominant{};
    std::vector<Figure> figures;
};

Run runExample(const json::Value& example, const fs::path& examples, const fs::path& archiveRoot) {
    const auto job = json::load((examples / core::wide(example.at("job").string())).wstring());
    const auto result = workflow::executeSolve(job, archiveRoot, {});
    if (result.archives.empty()) throw std::runtime_error("quickstart example produced no archive");
    Run run;
    run.title = example.at("title").string(); run.about = example.at("about").string();
    run.system = result.system; run.preset = job.at("system").at("preset").string();
    run.archive = shown(result.archives.front());
    run.samples = result.samplesPerTrajectory; run.accepted = result.acceptedSteps; run.rejected = result.rejectedSteps;
    const auto record = core::readNpz(result.archives.front());
    const auto series = core::loadTimeSeries(record, {{"time"}, {"values"}});
    core::FourierOptions fourier; fourier.removeMean = true; fourier.window = core::FourierWindow::hann;
    const auto state = static_cast<std::size_t>(example.at("fourier_state").number(0));
    run.fourierState = example.at("fourier_state_name").string("state " + std::to_string(state));
    run.dominant = series.fourier().transform(fourier)->dominantFrequency(state);
    const auto time = record.at("time").doubles(); const auto& values = record.at("values");
    run.start = time.front(); run.stop = time.back();
    const auto data = values.doubles(); const std::size_t width = values.shape.size() > 1 ? values.shape[1] : 1;
    run.finalState = "(";
    for (std::size_t k = 0; k < width; ++k) run.finalState += (k ? ", " : "") + number(data[data.size() - width + k]);
    run.finalState += ")";
    return run;
}

std::string uniqueFolder(const fs::path& base, const std::string& stem) {
    std::error_code ignored;
    for (int n = 1; n < 10000; ++n) {
        const auto name = stem + (n > 1 ? "-" + std::to_string(n) : "");
        if (!fs::exists(core::extendedLengthPath(base / core::wide(name)), ignored)) return name;
    }
    throw std::runtime_error("unable to allocate a quickstart figure folder");
}

void markdown(const Welcome& w, const std::vector<Run>& runs, const QuickstartOptions& options,
              const fs::path& figures, double seconds) {
    std::cout << w.overview << "\n## Quickstart results\n\n"
              << "Ran " << runs.size() << " bundled examples in " << number(seconds, 3) << " s. Results are in `"
              << shown(options.archiveRoot) << "`.\n\n"
              << "| Example | System (preset) | Time span | Samples | Steps accepted / rejected | "
                 "Dominant frequency | Final state |\n|---|---|---|---|---|---|---|\n";
    for (const auto& r : runs)
        std::cout << "| " << r.title << " | " << r.system << " (" << r.preset << ") | " << number(r.start) << " to "
                  << number(r.stop) << " | " << r.samples << " | " << r.accepted << " / " << r.rejected << " | "
                  << number(r.dominant) << " cycles per time unit (" << r.fourierState << ") | " << r.finalState << " |\n";
    std::cout << "\n";
    for (const auto& r : runs) std::cout << "- **" << r.title << "**: " << r.about << " Archive: `" << r.archive << "`\n";
    std::cout << "\n## Figures\n\nRendered by the DynSysLab Viewer CLI into `" << shown(figures) << "`:\n\n";
    for (const auto& r : runs)
        for (const auto& f : r.figures)
            std::cout << "- " << f.title << ": " << (f.error.empty() ? "`" + f.path + "`" : "not rendered (" + f.error + ")") << "\n";
    std::cout << "\n## Things you can ask next\n\n";
    for (const auto& item : w.next) std::cout << "- " << item << "\n";
    std::cout << "- \"Open these runs in the DynSysLab Viewer.\"\n\nThe desktop applications were not started; "
                 "open DynSysLab or the Viewer from the Start menu when you want them.\n";
}

void jsonReport(const Welcome& w, const std::vector<Run>& runs, const QuickstartOptions& options,
                const fs::path& figures, double seconds) {
    json::Value::Array items, next;
    for (const auto& r : runs) {
        json::Value::Array figs;
        for (const auto& f : r.figures)
            figs.push_back(json::Value::Object{{"title", f.title}, {"path", f.path}, {"rendered", f.error.empty()},
                                               {"error", f.error}});
        items.push_back(json::Value::Object{{"title", r.title}, {"about", r.about}, {"system", r.system},
            {"preset", r.preset}, {"archive", r.archive}, {"samples", double(r.samples)},
            {"accepted_steps", double(r.accepted)}, {"rejected_steps", double(r.rejected)},
            {"time_start", r.start}, {"time_stop", r.stop}, {"dominant_frequency", r.dominant}, {"fourier_state", r.fourierState},
            {"final_state", r.finalState}, {"figures", json::Value(std::move(figs))}});
    }
    for (const auto& item : w.next) next.push_back(item);
    std::cout << json::dump(json::Value::Object{{"schema", "dynsyslab.quickstart/v1"},
        {"version", std::string(core::version)}, {"welcome_markdown", w.overview},
        {"archive_root", shown(options.archiveRoot)}, {"archive_root_source", options.archiveRootSource},
        {"figure_folder", shown(figures)}, {"seconds", seconds}, {"runs", json::Value(std::move(items))},
        {"next", json::Value(std::move(next))}, {"desktop_applications_started", false}}) << '\n';
}

}  // namespace

int quickstart(const QuickstartOptions& options) {
    if (options.format != "markdown" && options.format != "json")
        throw std::invalid_argument("quickstart --format must be markdown or json");
    const auto started = std::chrono::steady_clock::now();
    const auto resources = locate();
    const auto index = json::load((resources.examples / L"quickstart.json").wstring());
    std::vector<Run> runs;
    for (const auto& example : index.at("examples").array()) runs.push_back(runExample(example, resources.examples, options.archiveRoot));
    const auto stem = [&] { const auto name = fs::path(core::wide(runs.front().archive)).stem().u8string();
                            const auto at = name.rfind("__"); return "quickstart-" + (at == std::string::npos ? name : name.substr(at + 2)); }();
    const auto base = options.archiveRoot.parent_path() / L"viewer_outputs";
    const auto figures = base / core::wide(uniqueFolder(base, stem));
    const auto cache = figures / L".viewer-cache";
    for (std::size_t k = 0; k < runs.size(); ++k)
        for (const auto& spec : index.at("examples").array()[k].at("figures").array()) {
            FigureRequest request{fs::path(core::wide(runs[k].archive)), spec.at("view").string(), {},
                                  figures / core::wide(spec.at("name").string() + ".png")};
            if (!spec.at("options").isNull())
                for (const auto& option : spec.at("options").array()) request.extra.push_back(core::wide(option.string()));
            Figure figure{spec.at("title").string(), shown(request.output), {}};
            if (!renderFigure(resources.viewer, request, cache, figure.error)) figure.path.clear();
            runs[k].figures.push_back(std::move(figure));
        }
    std::error_code ignored; fs::remove_all(core::extendedLengthPath(cache), ignored);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    const auto w = welcome(resources.welcome);
    if (options.format == "json") jsonReport(w, runs, options, figures, seconds);
    else markdown(w, runs, options, figures, seconds);
    return 0;
}

}  // namespace dynsyslab::cli
