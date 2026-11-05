#include "dynsyslab_viewer/NavigationPane.hpp"

#include <commctrl.h>

#include <algorithm>
#include <cwctype>
#include <functional>
#include <map>

namespace dsl {
namespace {

std::wstring lowered(std::wstring value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](wchar_t item) { return std::towlower(item); });
    return value;
}

HTREEITEM insert(HWND tree, HTREEITEM parent, const std::wstring& label,
                 LPARAM data = 0, bool bold = false) {
    TVINSERTSTRUCTW item{};
    item.hParent = parent;
    item.hInsertAfter = TVI_LAST;
    item.item.mask = TVIF_TEXT | TVIF_PARAM | TVIF_STATE;
    item.item.pszText = const_cast<wchar_t*>(label.c_str());
    item.item.lParam = data;
    item.item.stateMask = TVIS_BOLD;
    item.item.state = bold ? TVIS_BOLD : 0;
    return TreeView_InsertItem(tree, &item);
}

HTREEITEM findByData(HWND tree, HTREEITEM item, LPARAM target) {
    while (item) {
        TVITEMW data{};
        data.mask = TVIF_PARAM;
        data.hItem = item;
        TreeView_GetItem(tree, &data);
        if (data.lParam == target) return item;
        if (HTREEITEM found = findByData(tree, TreeView_GetChild(tree, item), target))
            return found;
        item = TreeView_GetNextSibling(tree, item);
    }
    return nullptr;
}

}  // namespace

bool NavigationPane::create(HWND parent, int controlId, const theme::Fonts& fonts) {
    window_ = CreateWindowExW(0, WC_TREEVIEWW, L"", WS_CHILD | WS_VISIBLE |
        WS_TABSTOP |
        WS_CLIPSIBLINGS |
        TVS_HASBUTTONS | TVS_HASLINES | TVS_LINESATROOT | TVS_SHOWSELALWAYS |
        TVS_FULLROWSELECT, 0, 0, 0, 0, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(controlId)),
        GetModuleHandleW(nullptr), nullptr);
    if (!window_) return false;
    theme::styleControl(window_, fonts.body);
    TreeView_SetBkColor(window_, theme::panel);
    TreeView_SetTextColor(window_, theme::text);
    TreeView_SetLineColor(window_, theme::border);
    TreeView_SetItemHeight(window_, theme::scalePixels(parent, 28));
    return true;
}

void NavigationPane::setCatalog(const Catalog* catalog) {
    catalog_ = catalog;
    rebuild();
}

void NavigationPane::setFilter(std::wstring filter) {
    filter_ = lowered(std::move(filter));
    rebuild();
}

void NavigationPane::setSort(CatalogTreeSort sort) {
    if (sort_ == sort) return;
    sort_ = sort;
    rebuild();
}

void NavigationPane::rebuild() {
    if (!window_) return;
    TreeView_DeleteAllItems(window_);
    if (!catalog_) return;
    const auto orderedTree = orderedCatalogTree(*catalog_, sort_);
    const auto& children = orderedTree.children;
    const auto& roots = orderedTree.roots;

    auto selfMatches = [&](int index) {
        const CatalogItem& item = catalog_->items[static_cast<std::size_t>(index)];
        std::string text = item.label + " " + item.system + " " + item.category + " " +
            item.role + " " + item.source;
        for (const std::string& tag : item.tags) text += " " + tag;
        for (const auto& property : item.properties)
            text += " " + property.first + " " + property.second;
        return filter_.empty() || lowered(wide(text)).find(filter_) != std::wstring::npos;
    };
    std::function<bool(int)> subtreeMatches = [&](int index) {
        if (selfMatches(index)) return true;
        return std::any_of(children[static_cast<std::size_t>(index)].begin(),
                           children[static_cast<std::size_t>(index)].end(), subtreeMatches);
    };
    std::function<int(int)> visibleMembers = [&](int index) {
        if (!subtreeMatches(index)) return 0;
        int count = 1;
        for (int child : children[static_cast<std::size_t>(index)])
            count += visibleMembers(child);
        return count;
    };
    std::function<int(int)> visibleRecords = [&](int index) {
        if (!subtreeMatches(index)) return 0;
        int count = catalog_->items[static_cast<std::size_t>(index)].kind ==
            "presentation_group" ? 0 : 1;
        for (int child : children[static_cast<std::size_t>(index)])
            count += visibleRecords(child);
        return count;
    };

    std::map<std::wstring, std::vector<int>> categories;
    for (int index : roots) {
        if (subtreeMatches(index))
            categories[wide(catalog_->items[static_cast<std::size_t>(index)].category)].push_back(index);
    }
    for (const auto& category : categories) {
        int members = 0;
        for (int index : category.second) members += visibleMembers(index);
        std::wstring label;
        if (category.first == L"Data") {
            int records = 0;
            for (int index : category.second) records += visibleRecords(index);
            label = L"Data  ·  " + std::to_wstring(records) + L" records";
        } else if (category.first == L"Archived solves") {
            label = category.first + L"  " + std::to_wstring(category.second.size()) +
                (category.second.size() == 1 ? L" batch" : L" batches") + L"  ·  " +
                std::to_wstring(members - static_cast<int>(category.second.size())) +
                L" archived objects";
        } else if (category.first == L"Systems") {
            label = category.first + L"  " + std::to_wstring(category.second.size()) +
                L" physical realizations  ·  " + std::to_wstring(members -
                static_cast<int>(category.second.size())) + L" linked objects";
        } else if (category.first == L"Training data") {
            int datasets = 0;
            for (int index : category.second) {
                const auto found = catalog_->items[static_cast<std::size_t>(index)].properties.find("dataset count");
                if (found != catalog_->items[static_cast<std::size_t>(index)].properties.end())
                    try { datasets += std::stoi(found->second); } catch (...) {}
            }
            label = category.first + L"  " + std::to_wstring(category.second.size()) +
                L" parameter sets  ·  " + std::to_wstring(datasets) + L" datasets";
        } else {
            label = category.first + L"  " + std::to_wstring(category.second.size());
            if (members > static_cast<int>(category.second.size()))
                label += L" models  ·  " + std::to_wstring(members - category.second.size()) + L" variants";
        }
        HTREEITEM root = insert(window_, TVI_ROOT, label, 0, true);
        std::function<void(HTREEITEM, int)> addItem = [&](HTREEITEM owner, int index) {
            const CatalogItem& item = catalog_->items[static_cast<std::size_t>(index)];
            std::wstring itemLabel = wide(item.label);
            const auto workflowCount = item.properties.find("workflow object count");
            const auto datasetCount = item.properties.find("dataset count");
            const auto itemCount = item.properties.find("item count");
            if (item.kind == "presentation_group" && itemCount != item.properties.end())
                itemLabel += L"  ·  " + wide(itemCount->second);
            else if (workflowCount != item.properties.end())
                itemLabel += L"  ·  " + wide(workflowCount->second) + L" objects";
            else if (datasetCount != item.properties.end())
                itemLabel += L"  ·  " + wide(datasetCount->second) + L" datasets";
            else if (item.variantCount > 0)
                itemLabel += L"  ·  " + std::to_wstring(item.variantCount) +
                    (item.variantCount == 1 ? L" variant" : L" variants");
            const bool group = item.kind == "presentation_group";
            HTREEITEM node = insert(window_, owner, itemLabel,
                group ? 0 : static_cast<LPARAM>(index + 1), group || item.primary);
            for (int child : children[static_cast<std::size_t>(index)])
                if (subtreeMatches(child)) addItem(node, child);
            if (!filter_.empty()) TreeView_Expand(window_, node, TVE_EXPAND);
        };
        if (category.first == L"Archived solves") {
            for (int index : category.second) addItem(root, index);
            TreeView_Expand(window_, root, TVE_EXPAND);
            continue;
        }
        std::map<std::wstring, std::vector<int>> systems;
        for (int index : category.second) {
            const std::wstring system = wide(catalog_->items[static_cast<std::size_t>(index)].system);
            systems[system.empty() ? L"Other" : system].push_back(index);
        }
        for (const auto& system : systems) {
            HTREEITEM parent = root;
            if (systems.size() > 1 || system.first != L"Other") {
                parent = insert(window_, root, system.first, 0, true);
            }
            for (int index : system.second) {
                addItem(parent, index);
            }
            if (!filter_.empty() && parent != root) TreeView_Expand(window_, parent, TVE_EXPAND);
        }
        if (!filter_.empty()) TreeView_Expand(window_, root, TVE_EXPAND);
    }
}

int NavigationPane::selectedIndex() const {
    if (!window_) return -1;
    const HTREEITEM selected = TreeView_GetSelection(window_);
    if (!selected) return -1;
    TVITEMW item{};
    item.mask = TVIF_PARAM;
    item.hItem = selected;
    TreeView_GetItem(window_, &item);
    return item.lParam > 0 ? static_cast<int>(item.lParam - 1) : -1;
}

void NavigationPane::selectIndex(int index) {
    if (!window_ || index < 0) return;
    if (HTREEITEM item = findByData(window_, TreeView_GetRoot(window_), index + 1)) {
        TreeView_SelectItem(window_, item);
        TreeView_EnsureVisible(window_, item);
    }
}

void NavigationPane::resetHorizontalScroll() {
    if (!window_) return;
    SetScrollPos(window_, SB_HORZ, 0, TRUE);
    SendMessageW(window_, WM_HSCROLL, MAKEWPARAM(SB_LEFT, 0), 0);
}

std::size_t NavigationPane::nodeCount() const {
    return window_ ? static_cast<std::size_t>(TreeView_GetCount(window_)) : 0;
}

}  // namespace dsl
