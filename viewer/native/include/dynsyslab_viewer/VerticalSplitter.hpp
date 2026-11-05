#pragma once

#include <windows.h>

namespace dsl {

inline constexpr UINT kVerticalSplitterMovedMessage = WM_APP + 0x4d;

class VerticalSplitter {
public:
    bool create(HWND parent, int controlId);
    void move(int x, int y, int width, int height);
    void show(bool visible);
    HWND handle() const { return window_; }

private:
    HWND window_ = nullptr;
    bool dragging_ = false;
    bool hovered_ = false;

    void reportPosition() const;
    LRESULT message(UINT message, WPARAM wparam, LPARAM lparam);
    static bool registerClass(HINSTANCE instance);
    static LRESULT CALLBACK procedure(
        HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

}  // namespace dsl
