#include "PlotPanelChrome.hpp"

#include <cmath>

namespace dsl {

RECT panelRectangle(const Gdiplus::RectF& area) {
    return {static_cast<LONG>(std::lround(area.X)),
        static_cast<LONG>(std::lround(area.Y)),
        static_cast<LONG>(std::lround(area.GetRight())),
        static_cast<LONG>(std::lround(area.GetBottom()))};
}

void drawPanelChrome(Gdiplus::Graphics& graphics, const Gdiplus::RectF& card,
                     bool highlighted, const theme::Fonts& fonts) {
    const BYTE alpha = highlighted ? 245 : 145;
    Gdiplus::SolidBrush glyph(theme::color(
        highlighted ? theme::accent : theme::plotMuted, alpha));
    const float startX = card.X + 10.0f;
    for (int column = 0; column < 2; ++column)
        for (int row = 0; row < 3; ++row)
            graphics.FillEllipse(&glyph, startX + column * 5.0f,
                card.Y + 8.0f + row * 5.0f, 2.4f, 2.4f);

    HDC dc = graphics.GetHDC(); Gdiplus::Font display(dc, fonts.body);
    graphics.ReleaseHDC(dc);
    Gdiplus::StringFormat format;
    format.SetAlignment(Gdiplus::StringAlignmentCenter);
    format.SetLineAlignment(Gdiplus::StringAlignmentCenter);
    const wchar_t close[] = L"×";
    graphics.DrawString(close, 1, &display,
        Gdiplus::RectF(card.GetRight() - 31.0f, card.Y + 1.0f, 25.0f, 27.0f),
        &format, &glyph);

    Gdiplus::Pen grip(theme::color(
        highlighted ? theme::accent : theme::plotMuted, alpha), 1.1f);
    for (int offset = 0; offset < 3; ++offset) {
        const float distance = 6.0f + offset * 4.0f;
        graphics.DrawLine(&grip, card.GetRight() - distance, card.GetBottom() - 3.0f,
            card.GetRight() - 3.0f, card.GetBottom() - distance);
    }
}

}  // namespace dsl
