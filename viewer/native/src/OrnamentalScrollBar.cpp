#include "dynsyslab_viewer/OrnamentalScrollBar.hpp"

#include "dynsyslab_viewer/Theme.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>

namespace dsl {
namespace {

int scaled(HWND window, int value) {
    const UINT dpi = GetDpiForWindow(window);
    return MulDiv(value, dpi ? static_cast<int>(dpi) : 96, 96);
}

void fill(HDC dc, const RECT& area, COLORREF color) {
    HBRUSH brush = CreateSolidBrush(color);
    FillRect(dc, &area, brush); DeleteObject(brush);
}

void frame(HDC dc, const RECT& area, COLORREF color) {
    HBRUSH brush = CreateSolidBrush(color);
    FrameRect(dc, &area, brush); DeleteObject(brush);
}

void diamond(HDC dc, int centerX, int centerY, int radius) {
    POINT points[]{{centerX, centerY - radius}, {centerX + radius, centerY},
        {centerX, centerY + radius}, {centerX - radius, centerY},
        {centerX, centerY - radius}};
    Polyline(dc, points, static_cast<int>(std::size(points)));
}

}  // namespace

int OrnamentalScrollBar::preferredWidth(HWND owner) { return scaled(owner, 20); }

void OrnamentalScrollBar::configure(int contentHeight, int viewportHeight) {
    contentHeight_ = std::max(0, contentHeight);
    viewportHeight_ = std::max(0, viewportHeight);
    maximum_ = std::max(0, contentHeight_ - viewportHeight_);
    position_ = std::clamp(position_, 0, maximum_);
    if (!visible()) dragging_ = false;
}

RECT OrnamentalScrollBar::barRect(HWND owner, const RECT& client) const {
    return {client.right - preferredWidth(owner), client.top,
        client.right, client.bottom};
}

RECT OrnamentalScrollBar::trackRect(HWND owner, const RECT& client) const {
    const RECT bar = barRect(owner, client);
    const int inset = scaled(owner, 5), endInset = scaled(owner, 8);
    return {bar.left + inset, bar.top + endInset,
        bar.right - inset, bar.bottom - endInset};
}

RECT OrnamentalScrollBar::thumbRect(HWND owner, const RECT& client) const {
    RECT track = trackRect(owner, client);
    const int trackHeight = std::max(0L, track.bottom - track.top);
    const int minimum = scaled(owner, 38);
    const int proportional = contentHeight_ > 0
        ? static_cast<int>(std::llround(static_cast<double>(trackHeight)
            * viewportHeight_ / contentHeight_)) : trackHeight;
    const int thumbHeight = std::clamp(proportional,
        std::min(minimum, trackHeight), trackHeight);
    const int travel = trackHeight - thumbHeight;
    const int offset = maximum_ > 0
        ? static_cast<int>(std::llround(static_cast<double>(travel)
            * position_ / maximum_)) : 0;
    track.top += offset; track.bottom = track.top + thumbHeight;
    InflateRect(&track, scaled(owner, 2), 0);
    return track;
}

bool OrnamentalScrollBar::contains(HWND owner, const RECT& client,
                                   int x, int y) const {
    const RECT bar = barRect(owner, client);
    const POINT point{x, y}; return PtInRect(&bar, point) != FALSE;
}

bool OrnamentalScrollBar::scrollBy(int pixels) {
    const int next = std::clamp(position_ + pixels, 0, maximum_);
    if (next == position_) return false;
    position_ = next; return true;
}

bool OrnamentalScrollBar::buttonDown(HWND owner, const RECT& client,
                                     int x, int y) {
    if (!contains(owner, client, x, y)) return false;
    if (!visible()) return true;
    const RECT thumb = thumbRect(owner, client);
    const POINT point{x, y};
    if (PtInRect(&thumb, point)) {
        dragging_ = true; dragOffset_ = y - thumb.top; SetCapture(owner);
    } else {
        scrollBy(y < thumb.top ? -viewportHeight_ : viewportHeight_);
    }
    return true;
}

bool OrnamentalScrollBar::buttonUp(HWND owner, const RECT& client, int x, int y) {
    const bool handled = dragging_ || contains(owner, client, x, y);
    cancelDrag(owner); return handled;
}

bool OrnamentalScrollBar::mouseMove(HWND owner, const RECT& client, int, int y) {
    if (!dragging_ || !visible()) return false;
    const RECT track = trackRect(owner, client);
    const RECT thumb = thumbRect(owner, client);
    const int thumbHeight = thumb.bottom - thumb.top;
    const int travel = std::max(0L, track.bottom - track.top - thumbHeight);
    const int thumbTop = std::clamp(y - dragOffset_,
        static_cast<int>(track.top), static_cast<int>(track.top) + travel);
    position_ = travel > 0 ? static_cast<int>(std::llround(
        static_cast<double>(thumbTop - track.top) * maximum_ / travel)) : 0;
    return true;
}

void OrnamentalScrollBar::cancelDrag(HWND owner) {
    if (!dragging_) return;
    dragging_ = false; if (GetCapture() == owner) ReleaseCapture();
}

void OrnamentalScrollBar::paint(HWND owner, HDC dc, const RECT& client) const {
    const RECT bar = barRect(owner, client);
    fill(dc, bar, theme::panel);
    HPEN edge = CreatePen(PS_SOLID, 1, theme::border);
    HGDIOBJ previous = SelectObject(dc, edge);
    MoveToEx(dc, bar.left, bar.top, nullptr); LineTo(dc, bar.left, bar.bottom);
    SelectObject(dc, previous); DeleteObject(edge);
    if (!visible()) return;
    const RECT track = trackRect(owner, client);
    fill(dc, track, theme::controlBackground); frame(dc, track, theme::borderStrong);
    const int centerX = (track.left + track.right) / 2;
    HPEN detail = CreatePen(PS_SOLID, 1, theme::accent);
    previous = SelectObject(dc, detail);
    diamond(dc, centerX, track.top - scaled(owner, 3), scaled(owner, 2));
    diamond(dc, centerX, track.bottom + scaled(owner, 3), scaled(owner, 2));
    const RECT thumb = thumbRect(owner, client);
    fill(dc, thumb, theme::selection); frame(dc, thumb, theme::accent);
    RECT inner = thumb; InflateRect(&inner, -scaled(owner, 2), -scaled(owner, 2));
    fill(dc, inner, theme::panelRaised);
    const int middle = (inner.top + inner.bottom) / 2;
    MoveToEx(dc, centerX, inner.top + scaled(owner, 3), nullptr);
    LineTo(dc, centerX, inner.bottom - scaled(owner, 3));
    diamond(dc, centerX, middle, scaled(owner, 3));
    const int curl = scaled(owner, 4), lift = scaled(owner, 7);
    POINT upper[]{{centerX, middle - scaled(owner, 2)},
        {centerX - curl, middle - lift}, {centerX - curl, middle - lift},
        {centerX, middle - 2 * lift}, {centerX + curl, middle - lift},
        {centerX + curl, middle - lift}, {centerX, middle - scaled(owner, 2)}};
    PolyBezier(dc, upper, static_cast<DWORD>(std::size(upper)));
    POINT lower[]{{centerX, middle + scaled(owner, 2)},
        {centerX - curl, middle + lift}, {centerX - curl, middle + lift},
        {centerX, middle + 2 * lift}, {centerX + curl, middle + lift},
        {centerX + curl, middle + lift}, {centerX, middle + scaled(owner, 2)}};
    PolyBezier(dc, lower, static_cast<DWORD>(std::size(lower)));
    SelectObject(dc, previous); DeleteObject(detail);
}

}  // namespace dsl
