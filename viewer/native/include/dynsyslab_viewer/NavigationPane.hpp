#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/CatalogTreeSort.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <string>

namespace dsl {

class NavigationPane {
public:
    bool create(HWND parent, int controlId, const theme::Fonts& fonts);
    void setCatalog(const Catalog* catalog);
    void setFilter(std::wstring filter);
    void setSort(CatalogTreeSort sort);
    int selectedIndex() const;
    void selectIndex(int index);
    void resetHorizontalScroll();
    std::size_t nodeCount() const;
    HWND handle() const { return window_; }

private:
    HWND window_ = nullptr;
    const Catalog* catalog_ = nullptr;
    std::wstring filter_;
    CatalogTreeSort sort_ = CatalogTreeSort::date;
    void rebuild();
};

}  // namespace dsl
