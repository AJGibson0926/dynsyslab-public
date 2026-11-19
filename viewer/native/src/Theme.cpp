#include "dynsyslab_viewer/Theme.hpp"

#include "dynsyslab_viewer/OrnamentalComboBox.hpp"
#include "TrajectoryColors.hpp"

#include <dwmapi.h>
#include <uxtheme.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>

namespace dsl::theme {
namespace {

HFONT makeFont(int size, int weight, const wchar_t* face = L"Segoe UI Variable Text") {
    HDC screen = GetDC(nullptr);
    const int height = -MulDiv(size, GetDeviceCaps(screen, LOGPIXELSY), 72);
    ReleaseDC(nullptr, screen);
    return CreateFontW(height, 0, 0, 0, weight, FALSE, FALSE, FALSE,
        DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
        CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, face);
}

int hexDigit(char value) {
    if (value >= '0' && value <= '9') return value - '0';
    if (value >= 'a' && value <= 'f') return value - 'a' + 10;
    if (value >= 'A' && value <= 'F') return value - 'A' + 10;
    return 0;
}

Gdiplus::Color hsv(double hue, double saturation, double light, BYTE alpha) {
    hue -= std::floor(hue);
    const double chroma = (1.0 - std::abs(2.0 * light - 1.0)) * saturation;
    const double section = hue * 6.0;
    const double x = chroma * (1.0 - std::abs(std::fmod(section, 2.0) - 1.0));
    double red = 0, green = 0, blueValue = 0;
    const int region = static_cast<int>(section) % 6;
    if (region == 0) { red = chroma; green = x; }
    else if (region == 1) { red = x; green = chroma; }
    else if (region == 2) { green = chroma; blueValue = x; }
    else if (region == 3) { green = x; blueValue = chroma; }
    else if (region == 4) { red = x; blueValue = chroma; }
    else { red = chroma; blueValue = x; }
    const double match = light - chroma / 2.0;
    return Gdiplus::Color(alpha,
        static_cast<BYTE>(255 * std::clamp(red + match, 0.0, 1.0)),
        static_cast<BYTE>(255 * std::clamp(green + match, 0.0, 1.0)),
        static_cast<BYTE>(255 * std::clamp(blueValue + match, 0.0, 1.0)));
}

}  // namespace

Fonts::Fonts(double scale) {
    const auto size = [scale](int points) { return std::max(1, static_cast<int>(std::lround(points * scale))); };
    title = makeFont(size(19), FW_SEMIBOLD, L"Bahnschrift SemiBold");
    heading = makeFont(size(11), FW_SEMIBOLD);
    body = makeFont(size(10), FW_NORMAL);
    caption = makeFont(size(9), FW_NORMAL);
    mono = makeFont(size(9), FW_NORMAL, L"Cascadia Mono");
}

Fonts::~Fonts() {
    for (HFONT font : {title, heading, body, caption, mono}) if (font) DeleteObject(font);
}

void configureTitleBar(HWND window) {
    const BOOL dark = FALSE;
    DwmSetWindowAttribute(window, 20, &dark, sizeof(dark));
    const COLORREF caption = background;
    const COLORREF captionText = text;
    DwmSetWindowAttribute(window, 35, &caption, sizeof(caption));
    DwmSetWindowAttribute(window, 36, &captionText, sizeof(captionText));
}

void styleControl(HWND control, HFONT font) {
    SendMessageW(control, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    SetWindowTheme(control, L"Explorer", nullptr);
    wchar_t className[32]{};
    GetClassNameW(control, className, static_cast<int>(std::size(className)));
    if (lstrcmpiW(className, L"ComboBox") == 0)
        dsl::styleOrnamentalComboBox(control);
}

int scalePixels(HWND window, int logicalPixels) {
    const UINT dpi = window ? GetDpiForWindow(window) : 96;
    return MulDiv(logicalPixels, dpi ? static_cast<int>(dpi) : 96, 96);
}

Metrics metrics(HWND window) {
    return {scalePixels(window, 12), scalePixels(window, 12),
        scalePixels(window, 8), scalePixels(window, 56),
        scalePixels(window, 34)};
}

Gdiplus::Color color(COLORREF value, BYTE alpha) {
    return Gdiplus::Color(alpha, GetRValue(value), GetGValue(value), GetBValue(value));
}

Gdiplus::Color parseColor(const std::string& value, std::size_t fallbackIndex, BYTE alpha) {
    static constexpr std::array<COLORREF, 8> palette = {
        RGB(59, 191, 184), RGB(78, 139, 221), RGB(232, 116, 91), RGB(179, 132, 224),
        RGB(233, 186, 81), RGB(100, 190, 122), RGB(230, 111, 153), RGB(143, 170, 184)
    };
    if (value.size() == 7 && value.front() == '#') {
        const BYTE red = static_cast<BYTE>(hexDigit(value[1]) * 16 + hexDigit(value[2]));
        const BYTE green = static_cast<BYTE>(hexDigit(value[3]) * 16 + hexDigit(value[4]));
        const BYTE blueValue = static_cast<BYTE>(hexDigit(value[5]) * 16 + hexDigit(value[6]));
        return Gdiplus::Color(alpha, red, green, blueValue);
    }
    return color(palette[fallbackIndex % palette.size()], alpha);
}

Gdiplus::Color plotColor(const std::string& value, std::size_t fallbackIndex) {
    const auto c = parseColor(value,fallbackIndex);
    const double luminance = .2126*c.GetR()+.7152*c.GetG()+.0722*c.GetB();
    const double factor = luminance > 145 ? 145/luminance : 1;
    return Gdiplus::Color(c.GetA(),static_cast<BYTE>(c.GetR()*factor),
        static_cast<BYTE>(c.GetG()*factor),static_cast<BYTE>(c.GetB()*factor));
}
Gdiplus::Color phaseColor(double phase, double magnitude) {
    const double normalized = (phase + 3.141592653589793) / (2.0 * 3.141592653589793);
    const double light = 0.36 + 0.24 * std::clamp(magnitude, 0.0, 1.0);
    return hsv(normalized, 0.82, light, 235);
}

Gdiplus::Color heatColor(double normalized) {
    if (!std::isfinite(normalized)) return color(panel);
    const auto index = std::min(255, static_cast<int>(std::clamp(normalized, 0., 1.) * 256));
    const auto rgb = trajectoryJet[static_cast<std::size_t>(index)];
    return {255, static_cast<BYTE>(rgb >> 16), static_cast<BYTE>(rgb >> 8), static_cast<BYTE>(rgb)};
}

}  // namespace dsl::theme
