#include "dynsyslab_viewer/PlotDragToolbar.hpp"

#include <algorithm>

namespace dsl {
namespace {
constexpr wchar_t kClassName[] = L"DynSysLabPlotDragToolbar";
constexpr int kArrangeId = 1;
constexpr int kPanId = 2;
constexpr int kDimensionalId = 3;
constexpr int kDimensionlessId = 4;
}

bool PlotDragToolbar::create(HWND parent, PlotCanvas& canvas, const theme::Fonts& fonts) {
    canvas_ = &canvas; fonts_ = &fonts;
    const auto instance = GetModuleHandleW(nullptr);
    WNDCLASSEXW type{sizeof(type)};
    type.lpfnWndProc = procedure;
    type.hInstance = instance;
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.lpszClassName = kClassName;
    if (!RegisterClassExW(&type) && GetLastError() != ERROR_CLASS_ALREADY_EXISTS)
        return false;
    window_ = CreateWindowExW(WS_EX_CONTROLPARENT, kClassName, L"Plot drag mode",
        WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
        0, 0, 0, 0, parent, nullptr, instance, this);
    return window_ != nullptr;
}

void PlotDragToolbar::setMode(PlotBackgroundDragMode mode) {
    canvas_->setBackgroundDragMode(mode);
    const bool arrange = mode == PlotBackgroundDragMode::rearrange_panels;
    SendMessageW(arrange_, BM_SETCHECK, arrange ? BST_CHECKED : BST_UNCHECKED, 0);
    SendMessageW(pan_, BM_SETCHECK, arrange ? BST_UNCHECKED : BST_CHECKED, 0);

}

void PlotDragToolbar::updateCoordinates() {
    if (!canvas_ || !dimensional_) return;
    const char* scaleNames[] = {"auto","linear","log","symlog"};
    for (int i = 0; i < 4; ++i)
        SendMessageW(scales_[i],BM_SETCHECK,canvas_->valueScaling() == scaleNames[i] ? BST_CHECKED : BST_UNCHECKED,0);
    const auto current = canvas_->displayCoordinates();
    const bool krom = canvas_->isKrom();
    for (const auto h : {realImaginary_,magnitudePhase_,domainColoring_,allModes_,noModes_})
        ShowWindow(h,krom ? SW_SHOWNA : SW_HIDE);
    const auto display = canvas_->complexDisplay();
    SendMessageW(realImaginary_,BM_SETCHECK,display == "real-imaginary" ? BST_CHECKED : BST_UNCHECKED,0);
    SendMessageW(magnitudePhase_,BM_SETCHECK,display == "magnitude-phase" ? BST_CHECKED : BST_UNCHECKED,0);
    SendMessageW(domainColoring_,BM_SETCHECK,display == "domain-coloring" ? BST_CHECKED : BST_UNCHECKED,0);
    EnableWindow(dimensional_, canvas_->canDisplayCoordinates("dimensional"));
    EnableWindow(dimensionless_, canvas_->canDisplayCoordinates("dimensionless"));
    SendMessageW(dimensional_, BM_SETCHECK, current == "dimensional" ? BST_CHECKED : BST_UNCHECKED, 0);
    SendMessageW(dimensionless_, BM_SETCHECK, current == "dimensionless" ? BST_CHECKED : BST_UNCHECKED, 0);
}

void PlotDragToolbar::layout() {
    RECT area{}; GetClientRect(window_, &area);
    const auto ui = theme::metrics(window_);
    const int gap = ui.compactGap;
    const int width = std::min(theme::scalePixels(window_, 128),
        std::max(1, static_cast<int>(area.right - gap * 4) / 4));
    const int height = std::max(1, static_cast<int>(area.bottom) / rowCount() - gap);
    MoveWindow(arrange_, 0, 0, width, height, TRUE);
    MoveWindow(pan_, width + gap, 0, width, height, TRUE);
    MoveWindow(dimensionless_, (width + gap) * 2 + gap, 0, width, height, TRUE);
    MoveWindow(dimensional_, (width + gap) * 3 + gap, 0, width, height, TRUE);
    const int modeWidth = std::min(theme::scalePixels(window_,164),std::max(1,static_cast<int>(area.right-gap*5)/4));
    const int scaleY = (rowCount()-1)*(height+gap);
    for (int i = 0; i < 4; ++i) MoveWindow(scales_[i],i*(modeWidth+gap),scaleY,modeWidth,height,TRUE);
    const int y = height+gap;
    MoveWindow(realImaginary_,0,y,modeWidth,height,TRUE);
    MoveWindow(magnitudePhase_,modeWidth+gap,y,modeWidth,height,TRUE);
    MoveWindow(domainColoring_,2*(modeWidth+gap),y,modeWidth,height,TRUE);
    MoveWindow(allModes_,3*(modeWidth+gap),y,modeWidth/2,height,TRUE);
    MoveWindow(noModes_,3*(modeWidth+gap)+modeWidth/2+gap,y,modeWidth/2,height,TRUE);
}

void PlotDragToolbar::paint(HDC target) {
    RECT area{}; GetClientRect(window_, &area);
    const auto brush = CreateSolidBrush(theme::background);
    FillRect(target, &area, brush); DeleteObject(brush);
}

LRESULT PlotDragToolbar::message(UINT messageValue, WPARAM wparam, LPARAM lparam) {
    if (messageValue == WM_CREATE) {
        auto button = [&](int id, const wchar_t* label, bool first, bool pushLike = true) {
            const auto result = CreateWindowExW(0, L"BUTTON", label,
                WS_CHILD | WS_VISIBLE | WS_TABSTOP
                    | (pushLike ? BS_RADIOBUTTON | BS_PUSHLIKE : BS_AUTORADIOBUTTON)
                    | (first ? WS_GROUP : 0),
                0, 0, 0, 0, window_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
                GetModuleHandleW(nullptr), nullptr);
            if (pushLike) theme::styleToggleButton(result, fonts_->body);
            else theme::styleControl(result, fonts_->body);
            return result;
        };
        arrange_ = button(kArrangeId, L"Move plots", true);
        pan_ = button(kPanId, L"Pan curves", false);
        dimensionless_ = button(kDimensionlessId, L"Dimensionless", true, false);
        dimensional_ = button(kDimensionalId, L"Dimensional", false, false);
        realImaginary_ = button(5,L"Real / Imaginary",true,false);
        magnitudePhase_ = button(6,L"Magnitude / Phase",false,false);
        domainColoring_ = button(7,L"Domain coloring",false,false);
        allModes_ = button(8,L"All",true);
        noModes_ = button(9,L"None",false);
        for (const auto h : {allModes_,noModes_}) SetWindowLongPtrW(h,GWL_STYLE,
            (GetWindowLongPtrW(h,GWL_STYLE) & ~BS_TYPEMASK) | BS_PUSHBUTTON);
        const wchar_t* scales[] = {L"Auto scale",L"Linear",L"Log",L"Symmetric log"};
        for (int i = 0; i < 4; ++i) scales_[i] = button(10+i,scales[i],i == 0,false);
        if (!arrange_ || !pan_ || !dimensionless_ || !dimensional_) return -1;
        updateCoordinates();
        setMode(canvas_->backgroundDragMode());
        return 0;
    }
    if (messageValue == WM_COMMAND && HIWORD(wparam) == BN_CLICKED) {
        if (LOWORD(wparam) == kArrangeId) setMode(PlotBackgroundDragMode::rearrange_panels);
        if (LOWORD(wparam) == kPanId) setMode(PlotBackgroundDragMode::pan_curves);
        if (LOWORD(wparam) == kDimensionalId) canvas_->setDisplayCoordinates("dimensional");
        if (LOWORD(wparam) == kDimensionlessId) canvas_->setDisplayCoordinates("dimensionless");
        if (LOWORD(wparam) == 5) canvas_->setComplexDisplay("real-imaginary");
        if (LOWORD(wparam) == 6) canvas_->setComplexDisplay("magnitude-phase");
        if (LOWORD(wparam) == 7) canvas_->setComplexDisplay("domain-coloring");
        if (LOWORD(wparam) == 8) canvas_->selectAllKromModes(true);
        if (LOWORD(wparam) == 9) canvas_->selectAllKromModes(false);
        if (LOWORD(wparam) >= 10 && LOWORD(wparam) <= 13) {
            const char* scales[] = {"auto","linear","log","symlog"};
            canvas_->setValueScaling(scales[LOWORD(wparam)-10]);
        }
        updateCoordinates();
        return 0;
    }
    if (messageValue == WM_SIZE) { layout(); return 0; }
    if (messageValue == WM_PAINT) {
        PAINTSTRUCT state{}; const auto target = BeginPaint(window_, &state);
        paint(target); EndPaint(window_, &state); return 0;
    }
    if (messageValue == WM_PRINTCLIENT) { paint(reinterpret_cast<HDC>(wparam)); return 0; }
    if (messageValue == WM_ERASEBKGND) return 1;
    if (messageValue == WM_CTLCOLORSTATIC || messageValue == WM_CTLCOLORBTN) {
        const auto target = reinterpret_cast<HDC>(wparam);
        SetTextColor(target, theme::muted); SetBkColor(target, theme::background);
        SetDCBrushColor(target, theme::background);
        return reinterpret_cast<LRESULT>(GetStockObject(DC_BRUSH));
    }
    return DefWindowProcW(window_, messageValue, wparam, lparam);
}

LRESULT CALLBACK PlotDragToolbar::procedure(HWND window, UINT message,
                                           WPARAM wparam, LPARAM lparam) {
    auto* self = reinterpret_cast<PlotDragToolbar*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        self = static_cast<PlotDragToolbar*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        self->window_ = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    }
    return self ? self->message(message, wparam, lparam)
                : DefWindowProcW(window, message, wparam, lparam);
}

}  // namespace dsl
