#include "dynsyslab_viewer/Theme.hpp"

#include <commctrl.h>
#include <iostream>
#include <windows.h>

namespace {

LRESULT CALLBACK fixtureProcedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    return DefWindowProcW(window, message, wparam, lparam);
}

bool require(bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    return condition;
}

}  // namespace

int main() {
    INITCOMMONCONTROLSEX controls{sizeof(controls), ICC_STANDARD_CLASSES};
    InitCommonControlsEx(&controls);
    WNDCLASSW fixture{};
    fixture.lpfnWndProc = fixtureProcedure;
    fixture.hInstance = GetModuleHandleW(nullptr);
    fixture.lpszClassName = L"DynSysLabComboTestFixture";
    RegisterClassW(&fixture);
    HWND parent = CreateWindowExW(WS_EX_COMPOSITED, fixture.lpszClassName, L"",
        WS_OVERLAPPEDWINDOW,
        0, 0, 320, 200, nullptr, nullptr, GetModuleHandleW(nullptr), nullptr);
    HWND combo = CreateWindowExW(WS_EX_CLIENTEDGE, L"COMBOBOX", L"",
        WS_CHILD | CBS_DROPDOWNLIST | WS_VSCROLL,
        10, 10, 240, 160, parent, nullptr, GetModuleHandleW(nullptr), nullptr);
    if (!require(parent && combo, "native ComboBox fixture could not be created"))
        return 1;

    dsl::theme::Fonts fonts;
    dsl::theme::styleControl(combo, fonts.body);
    bool ok = require((GetWindowLongPtrW(combo, GWL_EXSTYLE) & WS_EX_CLIENTEDGE) == 0,
                      "ornamental styling did not remove grey client-edge chrome");
    SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(L"RPE"));
    SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(L"Duffing"));
    SendMessageW(combo, CB_SETCURSEL, 1, 0);
    ok &= require(SendMessageW(combo, CB_GETCURSEL, 0, 0) == 1,
                  "selection messages no longer reach the native ComboBox");

    HDC screen = GetDC(nullptr);
    HDC target = CreateCompatibleDC(screen);
    HBITMAP bitmap = CreateCompatibleBitmap(screen, 240, 28);
    HGDIOBJ previous = SelectObject(target, bitmap);
    SendMessageW(combo, WM_PRINTCLIENT, reinterpret_cast<WPARAM>(target),
                 PRF_CLIENT);
    ok &= require(GetPixel(target, 1, 1) != CLR_INVALID,
                  "ornamental ComboBox did not render to an off-screen target");
    SelectObject(target, previous); DeleteObject(bitmap);
    DeleteDC(target); ReleaseDC(nullptr, screen);

    SendMessageW(combo, CB_SHOWDROPDOWN, TRUE, 0);
    ok &= require(SendMessageW(combo, CB_GETDROPPEDSTATE, 0, 0) != FALSE,
                  "native dropdown behavior was blocked by ornamental styling");
    SendMessageW(combo, CB_SHOWDROPDOWN, FALSE, 0);
    DestroyWindow(parent);
    return ok ? 0 : 1;
}
