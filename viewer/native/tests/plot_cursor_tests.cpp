#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void move(HWND window, POINT point) {
    SendMessageW(window, WM_MOUSEMOVE, 0, MAKELPARAM(point.x, point.y));
}

void drag(HWND window, POINT from, POINT to) {
    SendMessageW(window, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(from.x, from.y));
    move(window, to);
    SendMessageW(window, WM_LBUTTONUP, 0, MAKELPARAM(to.x, to.y));
    SendMessageW(window, WM_PAINT, 0, 0);
}

void checkCursorSequence(HWND window, POINT point, LPCWSTR cursor);

void checkDragModes(dsl::PlotCanvas& canvas, dsl::PlotDragToolbar& toolbar,
                    const dsl::Figure& figure) {
    const auto window = canvas.handle();
    const auto arrange = FindWindowExW(toolbar.handle(), nullptr, L"BUTTON", L"Move plots");
    const auto pan = FindWindowExW(toolbar.handle(), nullptr, L"BUTTON", L"Pan curves");
    require(arrange && pan, "plot drag toolbar buttons are missing");
    SendMessageW(arrange, BM_CLICK, 0, 0);
    require(SendMessageW(arrange, BM_GETCHECK, 0, 0) == BST_CHECKED
        && SendMessageW(pan, BM_GETCHECK, 0, 0) == BST_UNCHECKED,
        "Move plots must be the only selected tool");

    auto payload = std::make_shared<dsl::Payload>();
    payload->figures = {figure, figure};
    payload->figures[0].title = "First";
    payload->figures[1].title = "Second";
    canvas.setPayload(payload);
    canvas.setAxisLimits(0, dsl::AxisLimits{-3.0, 1.0, -1.0, 1.0});
    canvas.setAxisLimits(1, dsl::AxisLimits{-2.0, 2.0, -2.0, 2.0});
    SendMessageW(window, WM_PAINT, 0, 0);
    drag(window, {200, 150}, {600, 150});
    require(payload->figures[0].title == "Second" && payload->figures[1].title == "First",
            "Move plots did not rearrange panels from the background");
    auto axes = canvas.currentAxisPanels();
    require(axes[0].manual && axes[1].manual && axes[0].manual->xMin == -2.0
        && axes[1].manual->xMin == -3.0, "rearranging lost the panels' axis limits");

    SendMessageW(pan, BM_CLICK, 0, 0);
    require(canvas.backgroundDragMode() == dsl::PlotBackgroundDragMode::pan_curves
        && SendMessageW(pan, BM_GETCHECK, 0, 0) == BST_CHECKED
        && SendMessageW(arrange, BM_GETCHECK, 0, 0) == BST_UNCHECKED,
        "Pan curves did not exclusively activate its mode");
    SendMessageW(pan, BM_CLICK, 0, 0);
    require(SendMessageW(pan, BM_GETCHECK, 0, 0) == BST_CHECKED,
        "clicking the active tool must keep one mode selected");
    drag(window, {200, 150}, {600, 160});
    axes = canvas.currentAxisPanels();
    require(payload->figures[0].title == "Second" && axes[0].manual->xMin < -2.0
        && axes[1].manual->xMin == -3.0, "Pan curves moved a panel or the wrong axes");
    const auto pannedMinimum = axes[0].manual->xMin;

    // Borders retain their panel-reorder gesture even while Pan curves is active.
    drag(window, {10, 150}, {600, 150});
    axes = canvas.currentAxisPanels();
    require(payload->figures[0].title == "First" && axes[0].manual->xMin == -3.0
        && axes[1].manual->xMin == pannedMinimum,
        "border drag did not preserve panel order and panned limits");
    require(payload->figures[0].traces[0].x == figure.traces[0].x
        && payload->figures[1].traces[0].y == figure.traces[0].y,
        "display gestures mutated scientific samples");

    SendMessageW(window, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(200, 150));
    SendMessageW(arrange, BM_CLICK, 0, 0);
    require(GetCapture() != window, "switching tools left a plot drag captured");
    require(canvas.backgroundDragMode() == dsl::PlotBackgroundDragMode::rearrange_panels,
        "toolbar did not restore Move plots");
    checkCursorSequence(window, {380, 20}, IDC_ARROW);
    SendMessageW(window, WM_LBUTTONDOWN, 0, MAKELPARAM(380, 20));
    SendMessageW(window, WM_LBUTTONUP, 0, MAKELPARAM(380, 20));
    require(payload->figures.size() == 1 && payload->figures[0].title == "Second"
        && canvas.currentAxisPanels()[0].manual->xMin == pannedMinimum,
        "close button did not remove just its panel and preserve remaining limits");
}

void renderToolbar(HWND toolbar, const std::filesystem::path& path, int width) {
    const int height = dsl::theme::metrics(toolbar).toolbarHeight
        + dsl::theme::metrics(toolbar).compactGap;
    MoveWindow(toolbar, 0, 0, width, height, FALSE);
    Gdiplus::Bitmap bitmap(width, height, PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);
    const auto target = graphics.GetHDC();
    SendMessageW(toolbar, WM_PRINT, reinterpret_cast<WPARAM>(target), PRF_CLIENT | PRF_CHILDREN);
    graphics.ReleaseHDC(target);
    UINT count = 0, bytes = 0;
    Gdiplus::GetImageEncodersSize(&count, &bytes);
    std::vector<BYTE> storage(bytes);
    auto* encoders = reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    Gdiplus::GetImageEncoders(count, bytes, encoders);
    for (UINT index = 0; index < count; ++index) {
        if (wcscmp(encoders[index].MimeType, L"image/png") == 0) {
            require(bitmap.Save(path.c_str(), &encoders[index].Clsid) == Gdiplus::Ok,
                "toolbar preview could not be saved");
            return;
        }
    }
    throw std::runtime_error("PNG encoder is unavailable");
}

void checkCursorSequence(HWND window, POINT point, LPCWSTR cursor) {
    const auto expected = LoadCursorW(nullptr, cursor);
    move(window, point);
    require(GetCursor() == expected, "incorrect cursor after mouse movement");
    for (int step = 0; step < 16; ++step) {
        // Windows requests the cursor before delivering the next mouse move.
        // Check this intermediate state: checking only after the move hides flicker.
        const auto handled = SendMessageW(window, WM_SETCURSOR,
            reinterpret_cast<WPARAM>(window), MAKELPARAM(HTCLIENT, WM_MOUSEMOVE));
        require(GetCursor() == expected, "cursor changed between mouse moves");
        require(handled == TRUE, "canvas did not handle its client cursor");
        move(window, {point.x + step % 2, point.y});
        require(GetCursor() == expected, "cursor changed within the same region");
    }
}

}  // namespace

int wmain(int argc, wchar_t** argv) {
    ULONG_PTR token = 0;
    HWND parent = nullptr;
    const auto previousCursor = GetCursor();
    try {
        Gdiplus::GdiplusStartupInput input;
        require(Gdiplus::GdiplusStartup(&token, &input, nullptr) == Gdiplus::Ok,
                "GDI+ startup failed");
        const auto instance = GetModuleHandleW(nullptr);
        WNDCLASSW fixture{};
        fixture.lpfnWndProc = DefWindowProcW;
        fixture.hInstance = instance;
        fixture.hCursor = LoadCursorW(nullptr, IDC_ARROW);
        fixture.lpszClassName = L"DynSysLabPlotCursorTestFixture";
        require(RegisterClassW(&fixture) != 0, "fixture class registration failed");
        parent = CreateWindowExW(0, fixture.lpszClassName, L"", WS_OVERLAPPEDWINDOW,
            0, 0, 800, 600, nullptr, nullptr, instance, nullptr);
        require(parent != nullptr, "fixture window creation failed");
        dsl::theme::Fonts fonts;
        dsl::PlotCanvas canvas;
        require(dsl::PlotCanvas::registerClass(instance) && canvas.create(parent, fonts),
                "canvas creation failed");
        dsl::PlotDragToolbar toolbar;
        require(toolbar.create(parent, canvas, fonts), "plot drag toolbar creation failed");
        require(canvas.backgroundDragMode() == dsl::PlotBackgroundDragMode::rearrange_panels,
            "plot backgrounds must default to rearranging panels");
        const auto window = canvas.handle();
        MoveWindow(window, 0, 0, 800, 600, FALSE);
        auto payload = std::make_shared<dsl::Payload>();
        dsl::Figure figure;
        figure.kind = "time";
        figure.title = "Cursor regression";
        dsl::Trace trace;
        trace.x = {-1.0, 1.0};
        trace.y = {0.0, 0.0};
        figure.traces.push_back(trace);
        payload->figures.push_back(figure);
        canvas.setPayload(payload);
        canvas.setAxisLimits(0, dsl::AxisLimits{-1.0, 1.0, -1.0, 1.0});
        // Paint the real canvas into its hidden window to populate hit regions.
        SendMessageW(window, WM_PAINT, 0, 0);

        checkCursorSequence(window, {400, 150}, IDC_SIZEALL);
        checkCursorSequence(window, {2, 2}, IDC_ARROW);
        checkCursorSequence(window, {400, 20}, IDC_SIZEALL);
        checkCursorSequence(window, {780, 20}, IDC_ARROW);
        checkCursorSequence(window, {780, 580}, IDC_SIZENWSE);
        checkCursorSequence(window, {400, 150}, IDC_SIZEALL);

        // Find the rendered horizontal curve without coupling to axis margins.
        POINT curve{400, 200};
        for (; curve.y < 400; ++curve.y) {
            move(window, curve);
            if (GetCursor() == LoadCursorW(nullptr, IDC_ARROW)) break;
        }
        require(curve.y < 400, "curve hover did not select the pointer cursor");
        checkCursorSequence(window, curve, IDC_ARROW);
        checkCursorSequence(window, {400, 150}, IDC_SIZEALL);

        // Captured drags receive mouse moves without WM_SETCURSOR.
        toolbar.setMode(dsl::PlotBackgroundDragMode::pan_curves);
        checkCursorSequence(window, {780, 20}, IDC_ARROW);
        SendMessageW(window, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(400, 150));
        move(window, {405, 155});
        require(GetCursor() == LoadCursorW(nullptr, IDC_SIZEALL),
                "pan drag lost the move cursor");
        require(canvas.currentAxisPanels().front().manual->xMin < -1.0,
                "background drag no longer pans the plot");
        SendMessageW(window, WM_LBUTTONUP, 0, MAKELPARAM(405, 155));
        checkCursorSequence(window, {2, 2}, IDC_ARROW);

        // Non-client cursor handling must still reach Windows.
        SetCursor(LoadCursorW(nullptr, IDC_SIZEALL));
        SendMessageW(window, WM_SETCURSOR, reinterpret_cast<WPARAM>(window),
                     MAKELPARAM(HTCAPTION, WM_MOUSEMOVE));
        require(GetCursor() == LoadCursorW(nullptr, IDC_ARROW),
                "non-client cursor handling was swallowed");

        checkDragModes(canvas, toolbar, figure);
        if (argc > 1) {
            const std::filesystem::path directory(argv[1]);
            std::filesystem::create_directories(directory);
            for (int width : {420, 760}) {
                toolbar.setMode(dsl::PlotBackgroundDragMode::rearrange_panels);
                renderToolbar(toolbar.handle(), directory /
                    ("move-plots-" + std::to_string(width) + ".png"), width);
                toolbar.setMode(dsl::PlotBackgroundDragMode::pan_curves);
                renderToolbar(toolbar.handle(), directory /
                    ("pan-curves-" + std::to_string(width) + ".png"), width);
            }
        }

        DestroyWindow(parent); parent = nullptr;
        SetCursor(previousCursor);
        Gdiplus::GdiplusShutdown(token); token = 0;
        std::cout << "dynsyslab_plot_cursor_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (parent) DestroyWindow(parent);
        SetCursor(previousCursor);
        if (token) Gdiplus::GdiplusShutdown(token);
        std::cerr << "dynsyslab_plot_cursor_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
