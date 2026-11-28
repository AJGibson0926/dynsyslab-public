// Job kind "viewer_figure": one PNG rendered headlessly by the compiled Viewer CLI from an
// archive reference and complete figure settings.
#include "ViewerFigureJob.hpp"

#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow::detail {
namespace {

std::string numberText(double value) {
    char text[40];
    std::snprintf(text, sizeof(text), "%.17g", value);
    return text;
}

void textOption(std::vector<std::string>& args, const json::Value& figure, const char* key, const char* option) {
    if (!figure.at(key).isNull()) { args.emplace_back(option); args.push_back(figure.at(key).string()); }
}

void numberOption(std::vector<std::string>& args, const json::Value& figure, const char* key, const char* option) {
    if (!figure.at(key).isNull()) { args.emplace_back(option); args.push_back(numberText(figure.at(key).number())); }
}

void rangeOption(std::vector<std::string>& args, const json::Value& figure, const char* key, const char* option) {
    const auto& range = figure.at(key);
    if (range.isNull()) return;
    if (range.array().size() != 2) throw std::invalid_argument(std::string("viewer_figure: ") + key + " needs [low, high]");
    args.emplace_back(option);
    args.push_back(numberText(range.array()[0].number()) + ":" + numberText(range.array()[1].number()));
}

// "style":{"font_scale":F,"line_width_scale":S,"default_line_width":D}. A curve width is scaled
// by S and rounded to 0.01; a curve with no width takes D. A shared style is one job parameter.
double styledLineWidth(const json::Value& curve, const json::Value& style) {
    if (curve.at("line_width").isNull()) return style.at("default_line_width").number(0.0);
    const double scaled = curve.at("line_width").number() * style.at("line_width_scale").number(1.0);
    return style.at("line_width_scale").isNull() ? scaled : std::round(scaled * 100.0) / 100.0;
}

std::string checkedOutputName(const json::Value& job) {
    const auto name = job.at("output").string();
    if (name.empty() || name.find_first_of("/\\:") != std::string::npos || name.find("..") != std::string::npos
        || name.size() < 5 || name.compare(name.size() - 4, 4, ".png") != 0)
        throw std::invalid_argument("viewer_figure: output must be a plain file name ending in .png");
    return name;
}

#ifdef _WIN32
std::wstring quoted(const std::wstring& value) {  // CommandLineToArgvW rules
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

std::filesystem::path viewerCli() {
    std::wstring value(32768, L'\0');
    const DWORD length = GetModuleFileNameW(nullptr, value.data(), static_cast<DWORD>(value.size()));
    if (!length || length >= value.size()) throw std::runtime_error("viewer_figure: cannot locate this executable");
    value.resize(length);
    const auto path = std::filesystem::path(value).parent_path() / L"DynSysLabViewerCli.exe";
    if (!std::filesystem::is_regular_file(path))
        throw std::runtime_error("viewer_figure: DynSysLabViewerCli.exe is not next to this executable");
    return path;
}

std::string runViewer(const std::vector<std::string>& args) {
    std::wstring command = quoted(viewerCli().wstring());
    for (const auto& item : args) command += L" " + quoted(core::wide(item));
    SECURITY_ATTRIBUTES inherit{sizeof(inherit), nullptr, TRUE};
    HANDLE readPipe = nullptr, writePipe = nullptr;
    if (!CreatePipe(&readPipe, &writePipe, &inherit, 0)) throw std::runtime_error("viewer_figure: unable to create a pipe");
    SetHandleInformation(readPipe, HANDLE_FLAG_INHERIT, 0);
    HANDLE nul = CreateFileW(L"NUL", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                             &inherit, OPEN_EXISTING, 0, nullptr);
    STARTUPINFOW startup{}; startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdInput = startup.hStdError = nul;
    startup.hStdOutput = writePipe;
    PROCESS_INFORMATION process{};
    const BOOL started = CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW,
                                        nullptr, nullptr, &startup, &process);
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    CloseHandle(writePipe);
    if (!started) { CloseHandle(readPipe); throw std::runtime_error("viewer_figure: unable to start the Viewer CLI"); }
    std::string output; char buffer[4096]; DWORD count = 0;
    while (ReadFile(readPipe, buffer, sizeof(buffer), &count, nullptr) && count) output.append(buffer, count);
    CloseHandle(readPipe);
    const DWORD waited = WaitForSingleObject(process.hProcess, 300000);
    DWORD code = 1;
    if (waited == WAIT_OBJECT_0) GetExitCodeProcess(process.hProcess, &code);
    else TerminateProcess(process.hProcess, 1);
    CloseHandle(process.hThread); CloseHandle(process.hProcess);
    if (waited != WAIT_OBJECT_0) throw std::runtime_error("viewer_figure: the Viewer CLI timed out");
    if (code != 0) throw std::runtime_error("viewer_figure: the Viewer CLI failed (exit " + std::to_string(code) + ")");
    return output;
}
#else
std::string runViewer(const std::vector<std::string>&) {
    throw std::runtime_error("viewer_figure needs the Windows Viewer CLI");
}
#endif

}  // namespace

std::vector<std::string> viewerFigureArguments(const json::Value& job, const std::filesystem::path& output) {
    const auto& figure = job.at("figure");
    const auto view = figure.at("view").string("time-series");
    if (job.at("source").string().empty()) throw std::invalid_argument("viewer_figure needs a source archive");
    if (view != "time-series") {  // an archive view: the Viewer's own route for the source (render --view)
        std::vector<std::string> args{"render", "--source", job.at("source").string(), "--view", view,
                                      "--output", output.u8string()};
        if (!job.at("viewer_root").string().empty()) {  // the render's own data tree (executeViewerFigure)
            const auto root = job.at("viewer_root").string();
            args.insert(args.end(), {"--archives", root, "--cache", root + "/cache"});
        }
        numberOption(args, figure, "page", "--page");
        numberOption(args, figure, "figure", "--figure");
        textOption(args, figure, "panel", "--panel");
        textOption(args, figure, "title", "--title");
        textOption(args, figure, "scale", "--scale");                        // auto | linear | log | symlog
        textOption(args, figure, "legend", "--legend");
        textOption(args, figure, "coordinates", "--coordinates");
        textOption(args, figure, "modes", "--modes");
        textOption(args, figure, "complex_display", "--complex-display");
        numberOption(args, figure, "width", "--width");
        numberOption(args, figure, "height", "--height");
        rangeOption(args, figure, "x_range", "--x-range");
        rangeOption(args, figure, "y_range", "--y-range");
        return args;
    }
    if (figure.at("curves").array().empty()) throw std::invalid_argument("viewer_figure needs at least one curve");
    std::vector<std::string> args{view};
    if (!figure.at("style").at("font_scale").isNull()) {
        args.emplace_back("--font-scale"); args.push_back(numberText(figure.at("style").at("font_scale").number()));
    } else numberOption(args, figure, "font_scale", "--font-scale");
    args.insert(args.end(), {"--source", job.at("source").string(), "--output", output.u8string()});
    for (const auto& curve : figure.at("curves").array()) {
        args.emplace_back("--curve");
        args.push_back(curve.at("series").string() + "=" + curve.at("label").string());
        textOption(args, curve, "color", "--color");
        const auto width = styledLineWidth(curve, figure.at("style"));
        if (width > 0.0) { args.emplace_back("--line-width"); args.push_back(numberText(width)); }
        textOption(args, curve, "line_style", "--line-style");
        const auto side = curve.at("axis").string("left");  // "right": its own fitted right ordinate
        if (side != "left" && side != "right") throw std::invalid_argument("viewer_figure: curve axis must be left or right");
        if (side == "right") { args.emplace_back("--curve-axis"); args.emplace_back("right"); }
    }
    for (const auto& axis : figure.at("axes").array()) { args.emplace_back("--axis"); args.push_back(axis.string()); }
    textOption(args, figure, "x_label", "--x-label");
    textOption(args, figure, "value_label", "--value-label");
    textOption(args, figure, "time_label", "--time-label");
    textOption(args, figure, "time_key", "--time-key");
    textOption(args, figure, "scale", "--scale");
    textOption(args, figure, "title", "--title");
    numberOption(args, figure, "columns", "--columns");
    textOption(args, figure, "legend", "--legend");                      // "auto" | "top"
    if (figure.at("stack_axes").boolean(false)) args.emplace_back("--stack-axes");
    if (figure.at("spectrum").boolean(false)) args.emplace_back("--spectrum");  // Hann amplitude spectra
    numberOption(args, figure, "frequency_scale", "--frequency-scale");
    if (!figure.at("right_axis").isNull()) {  // plotted = offset + scale * right-axis value
        textOption(args, figure.at("right_axis"), "label", "--right-axis-label");
        numberOption(args, figure.at("right_axis"), "scale", "--right-axis-scale");
        numberOption(args, figure.at("right_axis"), "offset", "--right-axis-offset");
    }
    numberOption(args, figure, "width", "--width");
    numberOption(args, figure, "height", "--height");
    numberOption(args, figure, "start", "--start");
    numberOption(args, figure, "stop", "--stop");
    rangeOption(args, figure, "x_range", "--x-range");
    rangeOption(args, figure, "y_range", "--y-range");
    return args;
}

json::Value planViewerFigure(const json::Value& job) {
    (void)checkedOutputName(job);
    if (job.at("figure").at("view").string("time-series") == "time-series" && job.at("figure").at("curves").array().empty()) throw std::invalid_argument("viewer_figure needs at least one curve");
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.workflow-plan/v1")},
                                           {"kind", json::Value("viewer_figure")}, {"status", json::Value("planned")},
                                           {"output", job.at("output")}});
}

RunResult executeViewerFigure(const json::Value& job, const std::filesystem::path& archiveRoot) {
    const auto name = checkedOutputName(job);
    const auto recorded = job.at("output_directory").string();
    const auto directory = recorded.empty() ? archiveRoot / L"figures" : std::filesystem::u8path(recorded);
    std::filesystem::create_directories(directory);
    const auto output = directory / std::filesystem::u8path(name);
    std::error_code ignored;
    std::filesystem::remove(output, ignored);  // a figure job always renders afresh
    auto resolved = job;
    const auto& limits = job.at("axis_limits").at("reference");
    if (!limits.isNull()) {  // axis limits of a declared reference render, computed now
        json::Value reference(json::Value::Object{{"source", limits.at("source")}, {"figure", limits.at("figure")}});
        reference.at("figure").at("style") = job.at("figure").at("style");
        const auto scratch = directory / std::filesystem::u8path(name + ".reference.png");
        const auto report = json::parse(runViewer(viewerFigureArguments(reference, scratch)));
        std::filesystem::remove(scratch, ignored);
        const auto& axes = report.at("axes").array();
        if (axes.empty()) throw std::runtime_error("viewer_figure: the reference render reported no axes");
        const auto panel = static_cast<std::size_t>(limits.at("axis").number(0.0));
        const auto& axis = axes.at(panel);
        resolved.at("figure").at("x_range") = json::Value(json::Value::Array{axis.at("x_min"), axis.at("x_max")});
        resolved.at("figure").at("y_range") = json::Value(json::Value::Array{axis.at("y_min"), axis.at("y_max")});
    }
    // An archive view renders the source in a data tree of its own: a hard link to it (a copy across
    // volumes), indexed there and removed afterwards, so nothing is written beside the archives.
    const auto tree = directory / std::filesystem::u8path(name + ".viewer");
    if (resolved.at("figure").at("view").string("time-series") != "time-series") {
        std::filesystem::remove_all(core::extendedLengthPath(tree), ignored);
        std::filesystem::create_directories(tree);
        const auto from = std::filesystem::u8path(resolved.at("source").string());
        const auto linked = tree / from.filename();
        std::error_code failed;
        std::filesystem::create_hard_link(core::extendedLengthPath(from), core::extendedLengthPath(linked), failed);
        if (failed) std::filesystem::copy_file(core::extendedLengthPath(from), core::extendedLengthPath(linked));
        resolved.at("source") = json::Value(linked.generic_u8string());
        resolved.at("viewer_root") = json::Value(tree.generic_u8string());
    }
    const auto scratch = directory / std::filesystem::u8path(name + ".curves.npz");
    if (hasCurveSources(resolved)) resolved = assembleCurveSources(resolved, scratch);
    resolved = fitOrdinate(resolved);
    try { (void)runViewer(viewerFigureArguments(resolved, output)); }
    catch (...) { std::filesystem::remove(scratch, ignored); std::filesystem::remove_all(core::extendedLengthPath(tree), ignored); throw; }
    std::filesystem::remove(scratch, ignored);
    std::filesystem::remove_all(core::extendedLengthPath(tree), ignored);
    if (!std::filesystem::is_regular_file(output)) throw std::runtime_error("viewer_figure: no PNG was written");
    return {"viewer_figure", L"Figure written:\r\n" + output.wstring(), {output}, {}, {"viewer_figure"}};
}

}  // namespace dynsyslab::workflow::detail
