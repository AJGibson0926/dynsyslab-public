#include "dynsyslab_compute/WaveformPreview.hpp"

#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabWaveformPreview";
constexpr double kTwoPi = 6.283185307179586476925286766559;

struct State {
    HFONT body = nullptr;
    HFONT caption = nullptr;
    std::vector<double> frequencies{1.0};
    std::vector<double> amplitudes{1.0};
    std::vector<double> phases{0.0};
    double offset{};
    double start{};
    double stop{10.0};
    bool angular{};
    bool physical{};
    std::wstring error;
};

std::vector<double> numbers(const dsl::json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::vector<double> pluralOrScalar(const dsl::json::Value& settings,
                                   const char* plural, const char* singular,
                                   double fallback) {
    auto result = numbers(settings.at(plural));
    if (result.empty()) result.push_back(settings.at(singular).number(fallback));
    return result;
}

double sample(const State& state, double time) {
    double result = state.offset;
    for (std::size_t index = 0; index < state.frequencies.size(); ++index) {
        const double frequency = state.angular
            ? state.frequencies[index] : kTwoPi * state.frequencies[index];
        result += state.amplitudes[index]
            * std::sin(frequency * time + state.phases[index]);
    }
    return result;
}

void line(HDC dc, int x1, int y1, int x2, int y2) {
    MoveToEx(dc, x1, y1, nullptr); LineTo(dc, x2, y2);
}

void paint(HWND window, State& state, HDC dc) {
    RECT area{}; GetClientRect(window, &area);
    HBRUSH background = CreateSolidBrush(dsl::theme::panel);
    FillRect(dc, &area, background); DeleteObject(background);
    HBRUSH border = CreateSolidBrush(dsl::theme::border);
    FrameRect(dc, &area, border); DeleteObject(border);
    SetBkMode(dc, TRANSPARENT);
    SelectObject(dc, state.caption); SetTextColor(dc, dsl::theme::muted);
    RECT title{9, 4, area.right - 7, 23};
    const wchar_t* heading = state.physical
        ? L"PLANNED INPUT  •  physical time" : L"PLANNED INPUT  •  model time";
    DrawTextW(dc, heading, -1, &title, DT_LEFT | DT_SINGLELINE | DT_VCENTER);
    if (!state.error.empty()) {
        SelectObject(dc, state.body); SetTextColor(dc, dsl::theme::warning);
        RECT message{10, 29, area.right - 10, area.bottom - 8};
        DrawTextW(dc, state.error.c_str(), -1, &message,
                  DT_LEFT | DT_WORDBREAK | DT_NOPREFIX);
        return;
    }
    RECT plot{42, 28, std::max(48L, area.right - 12),
              std::max(35L, area.bottom - 27)};
    if (plot.right <= plot.left || plot.bottom <= plot.top) return;
    constexpr int count = 400;
    std::vector<double> values(count);
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (int index = 0; index < count; ++index) {
        const double time = state.start + (state.stop - state.start)
            * static_cast<double>(index) / static_cast<double>(count - 1);
        values[index] = sample(state, time);
        low = std::min(low, values[index]); high = std::max(high, values[index]);
    }
    double span = high - low;
    if (!(span > 0.0)) span = std::max(1.0, std::abs(high));
    low -= 0.08 * span; high += 0.08 * span;
    HPEN grid = CreatePen(PS_SOLID, 1, dsl::theme::border);
    HGDIOBJ previous = SelectObject(dc, grid);
    for (int tick = 0; tick <= 4; ++tick) {
        const int x = plot.left + (plot.right - plot.left) * tick / 4;
        const int y = plot.top + (plot.bottom - plot.top) * tick / 4;
        line(dc, x, plot.top, x, plot.bottom); line(dc, plot.left, y, plot.right, y);
    }
    SelectObject(dc, previous); DeleteObject(grid);
    HPEN curve = CreatePen(PS_SOLID, 2, dsl::theme::accent);
    previous = SelectObject(dc, curve);
    for (int index = 0; index < count; ++index) {
        const int x = plot.left + (plot.right - plot.left) * index / (count - 1);
        const int y = plot.bottom - static_cast<int>((values[index] - low)
            / (high - low) * static_cast<double>(plot.bottom - plot.top));
        if (!index) MoveToEx(dc, x, y, nullptr); else LineTo(dc, x, y);
    }
    SelectObject(dc, previous); DeleteObject(curve);
    SelectObject(dc, state.caption); SetTextColor(dc, dsl::theme::muted);
    RECT xLabel{plot.left, plot.bottom + 2, plot.right, area.bottom - 2};
    DrawTextW(dc, state.physical ? L"time (s)" : L"model time", -1, &xLabel,
              DT_CENTER | DT_SINGLELINE | DT_VCENTER);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<State*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_PAINT) {
        PAINTSTRUCT paintState{}; HDC dc = BeginPaint(window, &paintState);
        paint(window, *state, dc); EndPaint(window, &paintState); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_NCDESTROY) {
        delete state; SetWindowLongPtrW(window, GWLP_USERDATA, 0);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerClass() {
    WNDCLASSEXW type{sizeof(type)}; type.lpfnWndProc = procedure;
    type.hInstance = GetModuleHandleW(nullptr); type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(dsl::theme::panel); type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

HWND makeWaveformPreview(HWND parent, const dsl::theme::Fonts& fonts) {
    if (!registerClass()) return nullptr;
    auto* state = new State; state->body = fonts.body; state->caption = fonts.caption;
    HWND result = CreateWindowExW(0, kClassName, L"", WS_CHILD | WS_VISIBLE,
        0, 0, 0, 0, parent, nullptr, GetModuleHandleW(nullptr), state);
    if (!result) delete state;
    return result;
}

void setWaveformPreview(HWND preview, const std::string& encoded,
                        double startTime, double stopTime) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(preview, GWLP_USERDATA));
    if (!state) return;
    state->error.clear(); state->start = startTime; state->stop = stopTime;
    if (!std::isfinite(startTime) || !std::isfinite(stopTime) || stopTime <= startTime) {
        state->start = 0.0; state->stop = 10.0;
    }
    try {
        const auto settings = dsl::json::parse(encoded.empty() ? "{}" : encoded);
        state->physical = !settings.at("frequencies_hz").isNull()
            || !settings.at("frequency_hz").isNull();
        state->angular = !settings.at("angular_frequencies").isNull();
        state->frequencies = state->angular
            ? numbers(settings.at("angular_frequencies"))
            : state->physical
                ? pluralOrScalar(settings, "frequencies_hz", "frequency_hz", 1.0)
                : pluralOrScalar(settings, "frequencies", "frequency", 1.0);
        const bool model = !settings.at("model_amplitudes").isNull()
            || !settings.at("model_amplitude").isNull();
        state->amplitudes = model
            ? pluralOrScalar(settings, "model_amplitudes", "model_amplitude", 1.0)
            : pluralOrScalar(settings, "amplitudes", "amplitude", 1.0);
        state->phases = numbers(settings.at("phases_rad"));
        if (state->phases.empty()) state->phases.assign(state->frequencies.size(), 0.0);
        state->offset = settings.at(model ? "model_offset" : "offset").number();
        if (state->frequencies.size() != state->amplitudes.size()
            || state->frequencies.size() != state->phases.size())
            throw std::invalid_argument("Frequency, amplitude, and phase counts must match.");
    } catch (const std::exception& error) {
        state->error = dynsyslab::core::wide(error.what());
    }
    InvalidateRect(preview, nullptr, FALSE);
}

}  // namespace dsl::compute
