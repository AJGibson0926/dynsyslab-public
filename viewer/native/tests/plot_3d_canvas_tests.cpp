#include "dynsyslab_viewer/PlotCanvas.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void drag(HWND window, UINT down, WPARAM flags, POINT from, POINT to, UINT up) {
    SendMessageW(window, down, flags, MAKELPARAM(from.x, from.y));
    SendMessageW(window, WM_MOUSEMOVE, flags, MAKELPARAM(to.x, to.y));
    SendMessageW(window, up, 0, MAKELPARAM(to.x, to.y));
}

}  // namespace

int wmain() {
    ULONG_PTR token = 0;
    HWND parent = nullptr;
    try {
        Gdiplus::GdiplusStartupInput input;
        require(Gdiplus::GdiplusStartup(&token, &input, nullptr) == Gdiplus::Ok,
                "GDI+ startup failed");
        const auto instance = GetModuleHandleW(nullptr);
        WNDCLASSW fixture{};
        fixture.lpfnWndProc = DefWindowProcW;
        fixture.hInstance = instance;
        fixture.lpszClassName = L"DynSysLabPlot3DTestFixture";
        require(RegisterClassW(&fixture) != 0, "fixture class registration failed");
        parent = CreateWindowExW(0, fixture.lpszClassName, L"", WS_OVERLAPPEDWINDOW,
            0, 0, 900, 700, nullptr, nullptr, instance, nullptr);
        require(parent != nullptr, "fixture window creation failed");
        dsl::theme::Fonts fonts;
        dsl::PlotCanvas canvas;
        require(dsl::PlotCanvas::registerClass(instance) && canvas.create(parent, fonts),
                "3D canvas creation failed");
        MoveWindow(canvas.handle(), 0, 0, 900, 700, FALSE);

        auto payload = std::make_shared<dsl::Payload>();
        dsl::Figure figure;
        figure.kind = "three_d"; figure.title = "True 3D regression";
        figure.xLabel = "x"; figure.yLabel = "y"; figure.zLabel = "z";
        dsl::Trace trace;
        for (int index = 0; index <= 120; ++index) {
            const double t = index * 0.08;
            trace.x.push_back(std::cos(t)); trace.y.push_back(std::sin(t));
            trace.z.push_back(index / 60.0 - 1.0); trace.colorValues.push_back(t);
        }
        trace.encoding = "time"; figure.traces.push_back(std::move(trace));
        payload->figures.push_back(std::move(figure));
        canvas.setPayload(payload);
        SendMessageW(canvas.handle(), WM_PAINT, 0, 0);
        auto state = canvas.captureState();
        require(state.cameras3D.size() == 1, "paint did not initialize 3D camera state");
        const auto initial = state.cameras3D.at(0);

        drag(canvas.handle(), WM_LBUTTONDOWN, MK_LBUTTON,
             {380, 360}, {440, 325}, WM_LBUTTONUP);
        state = canvas.captureState();
        require(state.cameras3D.at(0).azimuthDegrees != initial.azimuthDegrees
            && state.cameras3D.at(0).elevationDegrees != initial.elevationDegrees,
            "left drag did not orbit the true 3D plot");

        const auto beforePan = state.cameras3D.at(0);
        drag(canvas.handle(), WM_LBUTTONDOWN, MK_LBUTTON | MK_CONTROL,
             {380, 360}, {420, 330}, WM_LBUTTONUP);
        state = canvas.captureState();
        require(state.cameras3D.at(0).panX != beforePan.panX
            && state.cameras3D.at(0).panY != beforePan.panY,
            "Ctrl+drag did not move the 3D camera frame");

        const auto beforeZoom = state.cameras3D.at(0).zoom;
        drag(canvas.handle(), WM_RBUTTONDOWN, MK_RBUTTON,
             {380, 360}, {380, 320}, WM_RBUTTONUP);
        state = canvas.captureState();
        require(state.cameras3D.at(0).zoom > beforeZoom,
                "right drag did not zoom the 3D camera");
        const auto draggedState = state.cameras3D.at(0);

        SendMessageW(canvas.handle(), WM_KEYDOWN, L'P', 0);
        state = canvas.captureState();
        require(state.cameras3D.at(0).projection == dsl::Plot3DProjection::orthographic,
                "P did not toggle the 3D projection");
        SendMessageW(canvas.handle(), WM_KEYDOWN, VK_HOME, 0);
        state = canvas.captureState();
        require(state.cameras3D.at(0) == dsl::Plot3DCameraState{},
                "Home did not reset the 3D camera");

        state.cameras3D[0] = draggedState;
        canvas.restoreState(std::move(state));
        require(canvas.captureState().cameras3D.at(0) == draggedState,
                "3D camera was not preserved across canvas state restoration");
        require(payload->figures[0].traces[0].z.front() == -1.0,
                "3D interaction mutated scientific samples");

        DestroyWindow(parent); parent = nullptr;
        Gdiplus::GdiplusShutdown(token); token = 0;
        std::cout << "dynsyslab_plot_3d_canvas_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        if (parent) DestroyWindow(parent);
        if (token) Gdiplus::GdiplusShutdown(token);
        std::cerr << "dynsyslab_plot_3d_canvas_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
