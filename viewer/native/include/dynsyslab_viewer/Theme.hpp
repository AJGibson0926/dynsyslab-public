#pragma once

#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>

#include <string>

namespace dsl::theme {

constexpr COLORREF plotBackground = RGB(255,255,255);
constexpr COLORREF plotText = RGB(17,24,39);
constexpr COLORREF plotMuted = RGB(71,85,105);
constexpr COLORREF plotBorder = RGB(166,177,189);
constexpr COLORREF plotGrid = RGB(215,221,228);
// Shared light scientific-workbench tokens. Scientific plot palettes remain
// independent so a presentation refactor cannot alter plotted meaning.
constexpr COLORREF background = RGB(244, 248, 251);     // Lab paper.
constexpr COLORREF panel = RGB(232, 241, 246);          // Instrument mist.
constexpr COLORREF panelRaised = RGB(255, 255, 255);    // Working surface.
constexpr COLORREF controlBackground = RGB(238, 244, 248);
constexpr COLORREF selection = RGB(214, 244, 240);      // Trace selection.
constexpr COLORREF border = RGB(200, 216, 226);
constexpr COLORREF borderStrong = RGB(151, 176, 191);
constexpr COLORREF text = RGB(19, 33, 43);              // Graphite ink.
constexpr COLORREF muted = RGB(82, 104, 121);           // Supporting slate.
constexpr COLORREF tertiary = RGB(97, 119, 134);
constexpr COLORREF disabled = RGB(139, 154, 164);
constexpr COLORREF accent = RGB(11, 125, 117);           // Oscilloscope teal.
constexpr COLORREF accentStrong = RGB(11, 94, 88);       // Teal text on light surfaces.
constexpr COLORREF blue = RGB(11, 111, 191);             // Focus cobalt.
constexpr COLORREF focus = blue;
constexpr COLORREF warning = RGB(163, 79, 0);
constexpr COLORREF danger = RGB(199, 62, 69);

struct Fonts {
    HFONT title = nullptr;
    HFONT heading = nullptr;
    HFONT body = nullptr;
    HFONT caption = nullptr;
    HFONT mono = nullptr;

    Fonts();
    ~Fonts();
    Fonts(const Fonts&) = delete;
    Fonts& operator=(const Fonts&) = delete;
};

struct Metrics {
    int margin = 0;
    int gap = 0;
    int compactGap = 0;
    int headerHeight = 0;
    int toolbarHeight = 0;
};

void configureTitleBar(HWND window);
void styleControl(HWND control, HFONT font);
// Theme a native push-like radio/check button while retaining its checked state.
void styleToggleButton(HWND control, HFONT font);
int scalePixels(HWND window, int logicalPixels);
Metrics metrics(HWND window);
Gdiplus::Color color(COLORREF value, BYTE alpha = 255);
Gdiplus::Color parseColor(const std::string& value, std::size_t fallbackIndex = 0,
                          BYTE alpha = 255);
// Keep ordinary line colors legible on white while preserving their hue.
Gdiplus::Color plotColor(const std::string& value, std::size_t fallbackIndex = 0);
Gdiplus::Color phaseColor(double phase, double magnitude = 1.0);
Gdiplus::Color heatColor(double normalized);

}  // namespace dsl::theme
