#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/ViewerWindow.hpp"
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Plot3D.hpp"

#include <windows.h>

#include <cmath>
#include <iostream>
#include <string>

void viewerSessionUsage();
int viewerSessionCommand(int argc, wchar_t** argv);
int viewerReportCommand(int argc, wchar_t** argv);

namespace {

std::string utf8(const std::wstring& value) {
    if (value.empty()) return {};
    const int count = WideCharToMultiByte(CP_UTF8, 0, value.data(),
        static_cast<int>(value.size()), nullptr, 0, nullptr, nullptr);
    std::string result(static_cast<std::size_t>(count), '\0');
    WideCharToMultiByte(CP_UTF8, 0, value.data(), static_cast<int>(value.size()),
        result.data(), count, nullptr, nullptr);
    return result;
}

std::string jsonEscaped(const std::string& value) {
    std::string result;
    for (const char item : value) {
        if (item == '"' || item == '\\') result += '\\';
        if (item == '\n') { result += "\\n"; continue; }
        if (item == '\r') { result += "\\r"; continue; }
        if (item == '\t') { result += "\\t"; continue; }
        result += item;
    }
    return result;
}

const char* pointerName(dsl::PlotPointerKind value) {
    if (value == dsl::PlotPointerKind::curve) return "curve_pointer";
    if (value == dsl::PlotPointerKind::pinned_probe) return "pinned_probe";
    return "background_move";
}

const char* roleName(dsl::MathTextRole value) {
    if (value == dsl::MathTextRole::subscript) return "subscript";
    if (value == dsl::MathTextRole::superscript) return "superscript";
    if (value == dsl::MathTextRole::fraction) return "fraction";
    return "normal";
}

bool approximately(double first, double second) {
    return std::abs(first - second) < 1.0e-8;
}

int verifyInteractions(bool json) {
    dsl::Figure figure;
    dsl::Trace curve; curve.label = "\\psi_{2}";
    curve.x = {0.0, 10.0}; curve.y = {0.0, 10.0};
    figure.traces.push_back(std::move(curve));
    dsl::PlotInteractionSession session(figure,
        {0.0, 0.0, 100.0, 100.0, 0.0, 10.0, 0.0, 10.0});

    session.move(50.0, 50.0);
    const auto hover = session.state();
    session.press(50.0, 50.0);
    const bool pinned = session.state().pinned.has_value();
    session.drag(75.0, 25.0); session.release();
    const auto dragged = session.state();
    session.move(10.0, 20.0);
    const auto background = session.state();
    session.press(10.0, 20.0); session.drag(20.0, 30.0); session.release();
    const auto panned = session.state();
    const auto math = dsl::parseMathText(L"$\\psi_{2}^{*}$");
    dsl::Plot3DCamera camera;
    const auto initialCamera = camera.state();
    const dsl::Plot3DBounds bounds{{-1,-1,-1},{1,1,1}};
    const dsl::Plot3DViewport viewport{0,0,640,480};
    const auto projected = camera.project({1,0,0},bounds,viewport);
    camera.orbit(40,-20); camera.pan(20,-10,viewport); camera.zoomBy(1.2);
    const auto movedCamera = camera.state();
    const auto rotated = camera.project({1,0,0},bounds,viewport);
    camera.reset();

    const bool passed = hover.pointer == dsl::PlotPointerKind::curve
        && hover.hover && approximately(hover.hover->x, 5.0) && pinned
        && dragged.pinned && approximately(dragged.pinned->x, 7.5)
        && background.pointer == dsl::PlotPointerKind::background
        && approximately(panned.transform.xMin, -1.0)
        && approximately(panned.transform.yMin, 1.0)
        && math.size() == 3 && math[1].role == dsl::MathTextRole::subscript
        && math[2].role == dsl::MathTextRole::superscript
        && projected.visible && rotated.visible
        && (!approximately(projected.x,rotated.x) || !approximately(projected.y,rotated.y))
        && movedCamera.zoom > initialCamera.zoom
        && camera.state() == initialCamera;
    if (json) {
        const double hoverX = hover.hover ? hover.hover->x : 0.0;
        const double hoverY = hover.hover ? hover.hover->y : 0.0;
        const double draggedX = dragged.pinned ? dragged.pinned->x : 0.0;
        const double draggedY = dragged.pinned ? dragged.pinned->y : 0.0;
        std::cout << "{\"schema\":\"dynsyslab.viewer.interaction-check/v1\",";
        std::cout << "\"passed\":" << (passed ? "true" : "false") << ',';
        std::cout << "\"hover_cursor\":\"" << pointerName(hover.pointer) << "\",";
        std::cout << "\"hover_coordinates\":[" << hoverX << ',' << hoverY << "],";
        std::cout << "\"pinned\":" << (pinned ? "true" : "false") << ',';
        std::cout << "\"dragged_coordinates\":[" << draggedX << ',' << draggedY << "],";
        std::cout << "\"background_cursor\":\"" << pointerName(background.pointer) << "\",";
        std::cout << "\"panned_limits\":[" << panned.transform.xMin << ','
                  << panned.transform.xMax << ',' << panned.transform.yMin << ','
                  << panned.transform.yMax << "],\"math_span_count\":" << math.size()
                  << ",\"three_d\":{\"projection\":\"perspective\",\"orbit\":true,"
                     "\"pan\":true,\"zoom\":true,\"reset\":true}}\n";
    } else {
        std::cout << (passed ? "PASS" : "FAIL")
                  << ": curve hover, pin/drag, background pan, math legend, and true 3D camera\n";
    }
    return passed ? 0 : 1;
}

void printMathSpans(const std::vector<dsl::MathTextSpan>& spans) {
    std::cout << '[';
    for (std::size_t i = 0; i < spans.size(); ++i) {
        if (i) std::cout << ',';
        const auto& span = spans[i];
        std::cout << "{\"text\":\"" << jsonEscaped(utf8(span.text))
                  << "\",\"role\":\"" << roleName(span.role) << "\"";
        if (!span.numerator.empty()) {
            std::cout << (span.role == dsl::MathTextRole::fraction
                ? ",\"numerator\":" : ",\"children\":");
            printMathSpans(span.numerator);
        }
        if (span.role == dsl::MathTextRole::fraction) {
            std::cout << ",\"denominator\":"; printMathSpans(span.denominator);
        }
        std::cout << '}';
    }
    std::cout << ']';
}

int parseMath(int argc, wchar_t** argv) {
    if (argc < 3) { std::cerr << "parse-math requires an expression\n"; return 2; }
    const auto spans = dsl::parseMathText(argv[2]);
    std::cout << "{\"schema\":\"dynsyslab.viewer.math-text/v1\",\"spans\":";
    printMathSpans(spans);
    std::cout << "}\n"; return 0;
}

void usage() {
    std::cout << "DynSysLabViewerCli\n"
              << "  verify-interactions [--json]\n"
              << "  parse-math <TeX-like-label>\n";
    std::cout << "  windows  (read-only status of running Viewer trees)\n";
    std::cout << "  reload <process-id>  (refresh that Viewer without activation)\n";
    std::cout << "  report SPEC.json --output report.html|report.pdf  (dynsyslab.report/v1 multi-page report)\n";
    viewerSessionUsage();
}

}  // namespace

int wmain(int argc, wchar_t** argv) {
    if (argc < 2 || std::wstring(argv[1]) == L"--help") { usage(); return 0; }
    const std::wstring command = argv[1];
    if (command == L"windows") {
        std::cout << dsl::json::dump(dsl::runningViewerWindows()) << '\n'; return 0;
    }
    if (command == L"reload" && argc == 3) {
        wchar_t* end = nullptr; const auto pid = std::wcstoul(argv[2],&end,10);
        if (!pid || !end || *end) return 2;
        const bool requested = dsl::reloadRunningViewer(pid);
        std::cout << dsl::json::dump(dsl::json::Value::Object{{"schema","dynsyslab.viewer.reload/v1"},
            {"process_id",static_cast<double>(pid)},{"requested",requested}}) << '\n';
        return requested ? 0 : 1;
    }
    if (command == L"verify-interactions")
        return verifyInteractions(argc > 2 && std::wstring(argv[2]) == L"--json");
    if (command == L"parse-math") return parseMath(argc, argv);
    if (command == L"report") return viewerReportCommand(argc, argv);
    if (command == L"time-series" || command == L"tree" || command == L"inspect" || command == L"panels" || command == L"render" || command == L"animate" || command == L"phase-space")
        return viewerSessionCommand(argc, argv);
    usage(); return 2;
}
