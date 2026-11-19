#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"
#include <cmath>
#include <limits>
#include <gdiplus.h>
#include <stdexcept>
#include <algorithm>
#include <chrono>

namespace dsl {
namespace {
struct GraphicsRuntime {
    ULONG_PTR token = 0;
    GraphicsRuntime() {
        Gdiplus::GdiplusStartupInput input;
        if (Gdiplus::GdiplusStartup(&token,&input,nullptr) != Gdiplus::Ok)
            throw std::runtime_error("Unable to initialize background plot rendering");
    }
    ~GraphicsRuntime() { if (token) Gdiplus::GdiplusShutdown(token); }
};
struct InvisibleSurface {
    HWND window = CreateWindowExW(WS_EX_NOACTIVATE | WS_EX_TOOLWINDOW,L"STATIC",L"",
        WS_POPUP,0,0,1,1,nullptr,nullptr,GetModuleHandleW(nullptr),nullptr);
    ~InvisibleSurface() { if (window) DestroyWindow(window); }
};
void savePng(Gdiplus::Bitmap& bitmap, const std::filesystem::path& path) {
    UINT count = 0, bytes = 0;
    Gdiplus::GetImageEncodersSize(&count,&bytes);
    std::vector<BYTE> buffer(bytes);
    auto* encoders = reinterpret_cast<Gdiplus::ImageCodecInfo*>(buffer.data());
    Gdiplus::GetImageEncoders(count,bytes,encoders);
    for (UINT i = 0; i < count; ++i) if (wcscmp(encoders[i].MimeType,L"image/png") == 0) {
        if (bitmap.Save(path.c_str(),&encoders[i].Clsid) != Gdiplus::Ok)
            throw std::runtime_error("Unable to write plot PNG: " + path.u8string());
        return;
    }
    throw std::runtime_error("PNG encoder is unavailable");
}
}
json::Value renderPlotPng(const Payload& input, const std::filesystem::path& output,
                          const PlotRenderOptions& options) {
    if (options.width < 320 || options.height < 240 || options.width > 4096 || options.height > 4096)
        throw std::runtime_error("Render dimensions must be 320..4096 by 240..4096");
    if (options.benchmarkScrollSteps < 0 || options.benchmarkScrollSteps > 240)
        throw std::runtime_error("Scroll profiling requires 0..240 frames");
    if (output.extension() != L".png") throw std::runtime_error("Render output must have a .png extension");
    auto payload = std::make_shared<Payload>(input);
    if (options.figure) {
        if (*options.figure >= payload->figures.size()) throw std::runtime_error("Figure index is out of range");
        payload->figures = {payload->figures[*options.figure]};
        payload->figures[0].dashboardPage = -1;
        payload->figures[0].gridSpan = payload->figures[0].rowSpan = 1;
        payload->pageLabels.clear(); payload->layout.clear(); payload->view = "single_figure";
    }
    if (payload->figures.empty()) throw std::runtime_error("There are no figures to render");
    if (!(options.fontScale >= 0.5 && options.fontScale <= 4.0))
        throw std::runtime_error("Font scale must be 0.5..4");
    GraphicsRuntime runtime;
    theme::Fonts fonts(options.fontScale);
    InvisibleSurface surface;
    PlotCanvas canvas; PlotDragToolbar toolbar;
    if (!surface.window || !PlotCanvas::registerClass(GetModuleHandleW(nullptr))
        || !canvas.create(surface.window,fonts)
        || (options.toolbar && !toolbar.create(surface.window,canvas,fonts)))
        throw std::runtime_error("Unable to create invisible plot surface");
    canvas.setPayload(payload);
    if (options.page < 0 || options.page >= canvas.pageCount()) throw std::runtime_error("Page index is out of range");
    canvas.setPage(options.page);
    const int top = options.toolbar ? toolbar.rowCount()*40 : 0;
    MoveWindow(surface.window,0,0,options.width,options.height,FALSE);
    MoveWindow(canvas.handle(),0,top,options.width,options.height-top,FALSE);
    if (options.toolbar) {
        MoveWindow(toolbar.handle(),0,0,options.width,top,FALSE);
        toolbar.updateCoordinates();
    }
    canvas.setScrollPosition(options.scrollY);
    if (IsWindowVisible(surface.window)) throw std::runtime_error("Background render surface must stay invisible");
    Gdiplus::Bitmap bitmap(options.width,options.height,PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);
    HDC dc = graphics.GetHDC();
    if (options.toolbar) SendMessageW(toolbar.handle(),WM_PRINT,reinterpret_cast<WPARAM>(dc),PRF_CLIENT | PRF_CHILDREN);
    SetViewportOrgEx(dc,0,top,nullptr);
    const auto redraw = [&] {
        const auto start = std::chrono::steady_clock::now();
        SendMessageW(canvas.handle(),WM_PRINTCLIENT,reinterpret_cast<WPARAM>(dc),PRF_CLIENT);
        return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    };
    if (options.xRange || options.yRange) {
        redraw();  // lays out the page so its visible panels are known
        const auto visible = canvas.visibleFigureIndices();
        for (std::size_t panel = 0; panel < visible.size(); ++panel) {
            const auto& figure = payload->figures[visible[panel]];
            const bool logarithmic = ordinateScale(figure).kind == ValueScaleKind::logarithmic;
            constexpr double inf = std::numeric_limits<double>::infinity();
            double x0 = inf, x1 = -inf, y0 = inf, y1 = -inf;
            for (const auto& trace : figure.traces) for (std::size_t i = 0; i < trace.x.size() && i < trace.y.size(); ++i) {
                const double x = trace.x[i], y = trace.y[i];
                if (!std::isfinite(x) || !std::isfinite(y)) continue;
                x0 = std::min(x0, x); x1 = std::max(x1, x);
                if (options.xRange && (x < options.xRange->first || x > options.xRange->second)) continue;
                if (logarithmic && y <= 0) continue;
                y0 = std::min(y0, y); y1 = std::max(y1, y);
            }
            AxisLimits limits;
            if (options.xRange) { limits.xMin = options.xRange->first; limits.xMax = options.xRange->second; }
            else { const double pad = 0.02 * (x1 - x0); limits.xMin = x0 - pad; limits.xMax = x1 + pad; }
            if (options.yRange) { limits.yMin = options.yRange->first; limits.yMax = options.yRange->second; }
            else if (logarithmic) { limits.yMin = y0 / 1.5; limits.yMax = y1 * 1.5; }
            else { const double pad = 0.05 * std::max(y1 - y0, 1e-12 * std::max(1.0, std::abs(y1))); limits.yMin = y0 - pad; limits.yMax = y1 + pad; }
            if (!canvas.setAxisLimits(panel, limits))
                throw std::runtime_error("Invalid axis range for panel " + std::to_string(panel)
                    + (logarithmic ? " (a logarithmic ordinate needs positive bounds)" : ""));
        }
    }
    const double initialMilliseconds = redraw();
    std::vector<double> timings;
    const int initialPosition = canvas.scrollPosition();
    for (int i = 0; i < options.benchmarkScrollSteps; ++i) {
        canvas.setScrollPosition(initialPosition+(i%20)*4);
        timings.push_back(redraw());
    }
    if (!timings.empty()) { canvas.setScrollPosition(initialPosition); redraw(); }
    const auto statistics = canvas.renderStatistics();
    graphics.ReleaseHDC(dc);
    std::sort(timings.begin(),timings.end());
    const auto percentile = [&](double fraction) {
        return timings.empty() ? 0. : timings[static_cast<std::size_t>((timings.size()-1)*fraction)];
    };
    if (!output.parent_path().empty()) std::filesystem::create_directories(output.parent_path());
    savePng(bitmap,output);
    json::Value::Array visible;
    for (const auto i : canvas.visibleFigureIndices()) visible.emplace_back(static_cast<double>(i));
    json::Value::Array axes;
    for (const auto& axis : canvas.currentAxisPanels())
        axes.emplace_back(json::Value::Object{{"title",narrow(axis.title)},
            {"x_min",axis.automatic.xMin},{"x_max",axis.automatic.xMax},
            {"y_min",axis.automatic.yMin},{"y_max",axis.automatic.yMax}});
    return json::Value::Object{{"schema","dynsyslab.viewer.render/v1"},{"output",std::filesystem::absolute(output).u8string()},
        {"width",static_cast<double>(options.width)},{"height",static_cast<double>(options.height)},
        {"page",static_cast<double>(options.page)},{"page_count",static_cast<double>(canvas.pageCount())},
        {"scroll_y",static_cast<double>(canvas.scrollPosition())},{"scroll_maximum",static_cast<double>(canvas.scrollMaximum())},
        {"render_statistics",json::Value::Object{{"panel_renders",static_cast<double>(statistics.panelRenders)},
            {"cache_hits",static_cast<double>(statistics.cacheHits)},{"cache_bytes",static_cast<double>(statistics.cacheBytes)},
            {"initial_frame_ms",initialMilliseconds},{"scroll_frames",static_cast<double>(timings.size())},
            {"scroll_median_ms",percentile(.5)},{"scroll_p95_ms",percentile(.95)}}},
        {"visible_figures",std::move(visible)},{"renderer","Viewer PlotCanvas"},{"visible_windows",0.0},{"axes",std::move(axes)}};
}
}
