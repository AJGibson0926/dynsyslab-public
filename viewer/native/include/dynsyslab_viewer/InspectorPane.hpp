#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <vector>

namespace dsl {

constexpr UINT kInspectorRelationMessage = WM_APP + 31;

class InspectorPane {
public:
    static bool registerClass(HINSTANCE instance);
    bool create(HWND parent, const theme::Fonts& fonts);
    void setItem(const Catalog* catalog, const CatalogItem* item);
    HWND handle() const { return window_; }

private:
    struct Link { RECT bounds{}; std::size_t relationIndex = 0; };
    HWND window_ = nullptr;
    HWND parent_ = nullptr;
    const theme::Fonts* fonts_ = nullptr;
    const Catalog* catalog_ = nullptr;
    const CatalogItem* item_ = nullptr;
    std::vector<Link> links_;
    int scroll_ = 0;
    int contentHeight_ = 0;

    void paint(HDC target);
    void mouseWheel(short delta);
    void click(int x, int y);
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

}  // namespace dsl
