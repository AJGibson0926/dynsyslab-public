#include "dynsyslab_viewer/VerticalSplitter.hpp"

#include "dynsyslab_viewer/Theme.hpp"

#include <algorithm>

namespace dsl {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabVerticalSplitter";

}  // namespace

bool VerticalSplitter::registerClass(HINSTANCE instance) {
    WNDCLASSEXW type{sizeof(type)};
    type.lpfnWndProc = procedure;
    type.hInstance = instance;
    type.hCursor = LoadCursorW(nullptr, IDC_SIZEWE);
    type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

bool VerticalSplitter::create(HWND parent, int controlId) {
    const HINSTANCE instance = reinterpret_cast<HINSTANCE>(
        GetWindowLongPtrW(parent, GWLP_HINSTANCE));
    if (!registerClass(instance)) return false;
    window_ = CreateWindowExW(0, kClassName, L"", WS_CHILD | WS_VISIBLE |
        WS_CLIPSIBLINGS, 0, 0, 0, 0, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(controlId)), instance, this);
    return window_ != nullptr;
}

void VerticalSplitter::move(int x, int y, int width, int height) {
    MoveWindow(window_, x, y, std::max(1, width), std::max(1, height), TRUE);
}

void VerticalSplitter::show(bool visible) {
    ShowWindow(window_, visible ? SW_SHOW : SW_HIDE);
}

void VerticalSplitter::reportPosition() const {
    POINT cursor{};
    if (!GetCursorPos(&cursor)) return;
    HWND parent = GetParent(window_);
    ScreenToClient(parent, &cursor);
    SendMessageW(parent, kVerticalSplitterMovedMessage,
        static_cast<WPARAM>(GetDlgCtrlID(window_)), static_cast<LPARAM>(cursor.x));
}

LRESULT VerticalSplitter::message(UINT messageValue, WPARAM wparam, LPARAM lparam) {
    if (messageValue == WM_SETCURSOR) {
        SetCursor(LoadCursorW(nullptr, IDC_SIZEWE));
        return TRUE;
    }
    if (messageValue == WM_LBUTTONDOWN) {
        dragging_ = true;
        SetCapture(window_);
        reportPosition();
        InvalidateRect(window_, nullptr, FALSE);
        return 0;
    }
    if (messageValue == WM_MOUSEMOVE) {
        if (!hovered_) {
            hovered_ = true;
            TRACKMOUSEEVENT tracking{sizeof(tracking), TME_LEAVE, window_, 0};
            TrackMouseEvent(&tracking);
            InvalidateRect(window_, nullptr, FALSE);
        }
        if (dragging_) reportPosition();
        return 0;
    }
    if (messageValue == WM_MOUSELEAVE) {
        hovered_ = false;
        InvalidateRect(window_, nullptr, FALSE);
        return 0;
    }
    if (messageValue == WM_LBUTTONUP) {
        if (dragging_) reportPosition();
        dragging_ = false;
        if (GetCapture() == window_) ReleaseCapture();
        InvalidateRect(window_, nullptr, FALSE);
        return 0;
    }
    if (messageValue == WM_CAPTURECHANGED) {
        dragging_ = false;
        InvalidateRect(window_, nullptr, FALSE);
        return 0;
    }
    if (messageValue == WM_PAINT) {
        PAINTSTRUCT state{};
        HDC target = BeginPaint(window_, &state);
        RECT area{};
        GetClientRect(window_, &area);
        HBRUSH background = CreateSolidBrush(theme::background);
        FillRect(target, &area, background);
        DeleteObject(background);
        const int center = (area.left + area.right) / 2;
        RECT grip{center - 1, area.top, center + 1, area.bottom};
        HBRUSH line = CreateSolidBrush(
            hovered_ || dragging_ ? theme::accent : theme::border);
        FillRect(target, &grip, line);
        DeleteObject(line);
        EndPaint(window_, &state);
        return 0;
    }
    const LRESULT result = DefWindowProcW(window_, messageValue, wparam, lparam);
    if (messageValue == WM_NCDESTROY) window_ = nullptr;
    return result;
}

LRESULT CALLBACK VerticalSplitter::procedure(
    HWND window, UINT messageValue, WPARAM wparam, LPARAM lparam) {
    auto* self = reinterpret_cast<VerticalSplitter*>(
        GetWindowLongPtrW(window, GWLP_USERDATA));
    if (messageValue == WM_NCCREATE) {
        self = static_cast<VerticalSplitter*>(
            reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        self->window_ = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    }
    return self ? self->message(messageValue, wparam, lparam) :
        DefWindowProcW(window, messageValue, wparam, lparam);
}

}  // namespace dsl
