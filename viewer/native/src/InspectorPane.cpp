#include "dynsyslab_viewer/InspectorPane.hpp"

#include <gdiplus.h>
#include <windowsx.h>

#include <algorithm>

namespace dsl {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabInspectorPane";

void fill(Gdiplus::Graphics& graphics, const Gdiplus::RectF& area, COLORREF color,
          float radius = 0.0f) {
    Gdiplus::SolidBrush brush(theme::color(color));
    if (radius <= 0.0f) graphics.FillRectangle(&brush, area);
    else {
        Gdiplus::GraphicsPath path;
        path.AddArc(area.X, area.Y, radius, radius, 180, 90);
        path.AddArc(area.GetRight() - radius, area.Y, radius, radius, 270, 90);
        path.AddArc(area.GetRight() - radius, area.GetBottom() - radius, radius, radius, 0, 90);
        path.AddArc(area.X, area.GetBottom() - radius, radius, radius, 90, 90);
        path.CloseFigure();
        graphics.FillPath(&brush, &path);
    }
}

void text(Gdiplus::Graphics& graphics, const std::wstring& value, const Gdiplus::RectF& area,
          HFONT font, COLORREF color, Gdiplus::StringAlignment alignment = Gdiplus::StringAlignmentNear) {
    HDC dc = graphics.GetHDC();
    Gdiplus::Font display(dc, font);
    graphics.ReleaseHDC(dc);
    Gdiplus::SolidBrush brush(theme::color(color));
    Gdiplus::StringFormat format;
    format.SetAlignment(alignment);
    format.SetLineAlignment(Gdiplus::StringAlignmentNear);
    format.SetTrimming(Gdiplus::StringTrimmingEllipsisCharacter);
    graphics.DrawString(value.c_str(), -1, &display, area, &format, &brush);
}

}  // namespace

bool InspectorPane::registerClass(HINSTANCE instance) {
    WNDCLASSEXW type{sizeof(type)};
    type.style = CS_HREDRAW | CS_VREDRAW;
    type.lpfnWndProc = procedure;
    type.hInstance = instance;
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.lpszClassName = kClassName;
    type.hbrBackground = CreateSolidBrush(theme::panel);
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

bool InspectorPane::create(HWND parent, const theme::Fonts& fonts) {
    parent_ = parent;
    fonts_ = &fonts;
    window_ = CreateWindowExW(0, kClassName, L"", WS_CHILD | WS_VISIBLE |
        WS_CLIPSIBLINGS,
        0, 0, 0, 0, parent, nullptr, GetModuleHandleW(nullptr), this);
    return window_ != nullptr;
}

void InspectorPane::setItem(const Catalog* catalog, const CatalogItem* item) {
    catalog_ = catalog;
    item_ = item;
    links_.clear();
    scroll_ = 0;
    InvalidateRect(window_, nullptr, FALSE);
}

void InspectorPane::paint(HDC target) {
    RECT client{};
    GetClientRect(window_, &client);
    const float scale = static_cast<float>(GetDpiForWindow(window_)) / 96.0f;
    const auto px = [scale](float logical) { return logical * scale; };
    const float inset = px(18.0f);
    const float width = static_cast<float>(client.right);
    Gdiplus::Graphics graphics(target);
    graphics.SetSmoothingMode(Gdiplus::SmoothingModeAntiAlias);
    graphics.Clear(theme::color(theme::panel));
    links_.clear();
    if (!item_ || !fonts_) {
        text(graphics, L"Select an object to inspect its details and lineage.",
             Gdiplus::RectF(px(20), px(28), width - px(40), px(80)),
             fonts_ ? fonts_->body : static_cast<HFONT>(GetStockObject(DEFAULT_GUI_FONT)), theme::muted);
        return;
    }
    float y = px(20.0f) - static_cast<float>(scroll_);
    fill(graphics, Gdiplus::RectF(inset, y, width - 2 * inset, px(28)),
         theme::accent, px(12));
    text(graphics, wide(item_->kind), Gdiplus::RectF(
         px(28), y + px(3), width - px(56), px(23)),
         fonts_->caption, theme::panelRaised);
    y += px(42);
    if (!item_->subtitle.empty()) {
        text(graphics, wide(item_->subtitle), Gdiplus::RectF(
             inset, y, width - 2 * inset, px(34)),
             fonts_->body, theme::muted);
        y += px(43);
    }

    float tagX = inset;
    for (const std::string& tag : item_->tags) {
        const float tagWidth = std::min(width - 2 * inset,
            px(22.0f + 6.3f * static_cast<float>(tag.size())));
        if (tagX + tagWidth > width - inset) {
            tagX = inset;
            y += px(34);
        }
        fill(graphics, Gdiplus::RectF(tagX, y, tagWidth, px(28)),
             theme::panelRaised, px(9));
        text(graphics, wide(tag), Gdiplus::RectF(
             tagX + px(9), y + px(3), tagWidth - px(18), px(23)),
             fonts_->caption, theme::text);
        tagX += tagWidth + px(7);
    }
    y += px(44);
    text(graphics, L"Properties", Gdiplus::RectF(
         inset, y, width - 2 * inset, px(24)),
         fonts_->caption, theme::accentStrong);
    y += px(31);
    for (const auto& property : item_->properties) {
        text(graphics, wide(property.first), Gdiplus::RectF(
             inset, y, width * 0.38f - inset, px(30)),
             fonts_->body, theme::muted);
        text(graphics, wide(property.second), Gdiplus::RectF(
             width * 0.40f, y, width * 0.56f, px(30)),
             fonts_->body, theme::text);
        y += px(34);
    }
    y += px(14);
    text(graphics, L"Available views", Gdiplus::RectF(
         inset, y, width - 2 * inset, px(24)),
         fonts_->caption, theme::accentStrong);
    y += px(31);
    for (const std::string& view : item_->views) {
        const auto found = catalog_ ? catalog_->viewLabels.find(view) :
                                     std::map<std::string, std::string>::const_iterator{};
        const std::string label = catalog_ && found != catalog_->viewLabels.end() ? found->second : view;
        text(graphics, L"•  " + wide(label), Gdiplus::RectF(
             px(22), y, width - px(44), px(30)),
             fonts_->body, theme::text);
        y += px(34);
    }
    if (!item_->relations.empty()) {
        y += px(17);
        text(graphics, L"Related objects", Gdiplus::RectF(
             inset, y, width - 2 * inset, px(24)),
             fonts_->caption, theme::accentStrong);
        y += px(32);
        for (std::size_t index = 0; index < item_->relations.size(); ++index) {
            const Relation& relation = item_->relations[index];
            const Gdiplus::RectF area(inset, y, width - 2 * inset, px(62));
            fill(graphics, area, theme::panelRaised, px(8));
            text(graphics, wide(relation.kind), Gdiplus::RectF(
                 px(28), y + px(6), width - px(56), px(22)),
                 fonts_->caption, theme::accentStrong);
            text(graphics, wide(relation.label), Gdiplus::RectF(
                 px(28), y + px(29), width - px(56), px(28)),
                 fonts_->body, theme::text);
            links_.push_back({RECT{static_cast<LONG>(inset), static_cast<LONG>(y),
                client.right - static_cast<LONG>(inset),
                static_cast<LONG>(y + px(62))}, index});
            y += px(70);
        }
    }
    y += px(20);
    text(graphics, L"Source", Gdiplus::RectF(
         inset, y, width - 2 * inset, px(24)),
         fonts_->caption, theme::accentStrong);
    y += px(29);
    text(graphics, wide(item_->source), Gdiplus::RectF(
         inset, y, width - 2 * inset, px(72)),
         fonts_->caption, theme::muted);
    contentHeight_ = static_cast<int>(y + px(88) + scroll_);
}

void InspectorPane::mouseWheel(short delta) {
    RECT client{}; GetClientRect(window_, &client);
    scroll_ = std::clamp(scroll_ - static_cast<int>(delta / 2), 0,
                         std::max(0, contentHeight_ - static_cast<int>(client.bottom)));
    InvalidateRect(window_, nullptr, FALSE);
}

void InspectorPane::click(int x, int y) {
    POINT point{x, y};
    for (const Link& link : links_) {
        if (PtInRect(&link.bounds, point)) {
            PostMessageW(parent_, kInspectorRelationMessage,
                         static_cast<WPARAM>(link.relationIndex), 0);
            break;
        }
    }
}

LRESULT CALLBACK InspectorPane::procedure(HWND window, UINT message,
                                           WPARAM wparam, LPARAM lparam) {
    auto* self = reinterpret_cast<InspectorPane*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        self = static_cast<InspectorPane*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        self->window_ = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    }
    if (!self) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_PRINTCLIENT) { self->paint(reinterpret_cast<HDC>(wparam)); return 0; }
    if (message == WM_PAINT) {
        PAINTSTRUCT paint{}; HDC target = BeginPaint(window, &paint);
        self->paint(target); EndPaint(window, &paint); return 0;
    }
    if (message == WM_MOUSEWHEEL) { self->mouseWheel(GET_WHEEL_DELTA_WPARAM(wparam)); return 0; }
    if (message == WM_LBUTTONUP) { self->click(GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam)); return 0; }
    if (message == WM_ERASEBKGND) return 1;
    return DefWindowProcW(window, message, wparam, lparam);
}

}  // namespace dsl
