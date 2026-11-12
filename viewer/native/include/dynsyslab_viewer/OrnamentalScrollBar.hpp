#pragma once

#include <windows.h>

namespace dsl {

// Client-drawn vertical scrollbar shared by DynSysLab scientific canvases.
// It avoids platform scrollbar chrome while retaining wheel, track, and drag
// interaction. The owner remains responsible for invalidating after changes.
class OrnamentalScrollBar {
public:
    static int preferredWidth(HWND owner);
    void configure(int contentHeight, int viewportHeight);
    int position() const { return position_; }
    bool visible() const { return maximum_ > 0; }
    bool contains(HWND owner, const RECT& client, int x, int y) const;
    bool buttonDown(HWND owner, const RECT& client, int x, int y);
    bool buttonUp(HWND owner, const RECT& client, int x, int y);
    bool mouseMove(HWND owner, const RECT& client, int x, int y);
    bool scrollBy(int pixels);
    void cancelDrag(HWND owner);
    void paint(HWND owner, HDC dc, const RECT& client) const;
private:
    int contentHeight_{};
    int viewportHeight_{};
    int maximum_{};
    int position_{};
    int dragOffset_{};
    bool dragging_{};
    RECT barRect(HWND owner, const RECT& client) const;
    RECT trackRect(HWND owner, const RECT& client) const;
    RECT thumbRect(HWND owner, const RECT& client) const;
};

}  // namespace dsl
