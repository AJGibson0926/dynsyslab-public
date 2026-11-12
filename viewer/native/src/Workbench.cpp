#include "dynsyslab_viewer/Workbench.hpp"
#include "AxisRangeDialog.hpp"
#include "PanelPickerDialog.hpp"
#include "dynsyslab_viewer/Help.hpp"
#include "dynsyslab_viewer/FolderTab.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"
#include "dynsyslab_viewer/StartupDiagnostics.hpp"
#include "dynsyslab_resource.h"
#include <commctrl.h>
#include <gdiplus.h>
#include <shellapi.h>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <stdexcept>
namespace dsl {
namespace {
constexpr wchar_t kClassName[] = L"DynSysLabViewerWorkbench";
constexpr int kNavigationId = 100;
constexpr int kSearchId = 101;
constexpr int kRefreshId = 102;
constexpr int kPreviousId = 103;
constexpr int kNextId = 104;
constexpr int kHelpId = 105;
constexpr int kCopyId = 106;
constexpr int kSaveId = 107;
constexpr int kComputeId = 108;
constexpr int kPageSelectId = 109;
constexpr int kAxesId = 110;
constexpr int kCompareId = 111;
constexpr int kLeftSplitterId = 112;
constexpr int kAddPanelId = 113;
constexpr int kMainTabId = 114;
constexpr int kObservablesTabId = 115;
constexpr int kPhaseContextId = 116;
constexpr int kSortDateId = 117;
constexpr int kSortTrajectoriesId = 118;
constexpr int kSortSizeId = 119;
constexpr int kViewSelectId = 200;
constexpr UINT_PTR kTaskTimer = 301;
HWND button(HWND parent, int id, const wchar_t* label, HFONT font) {
    HWND result = CreateWindowExW(0, L"BUTTON", label, WS_CHILD | WS_VISIBLE |
        WS_TABSTOP | BS_OWNERDRAW, 0, 0, 0, 0, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(result, font);
    return result;
}
std::wstring buttonText(HWND window) {
    const int count = GetWindowTextLengthW(window);
    std::wstring result(static_cast<std::size_t>(count + 1), L'\0');
    GetWindowTextW(window, result.data(), count + 1);
    result.resize(static_cast<std::size_t>(count));
    return result;
}
void addTooltip(HWND tooltip, HWND control, const wchar_t* text) {
    TTTOOLINFOW info{};
    info.cbSize = TTTOOLINFOW_V1_SIZE;
    info.uFlags = TTF_SUBCLASS;
    info.hwnd = control;
    info.uId = 1;
    GetClientRect(control, &info.rect);
    info.lpszText = const_cast<wchar_t*>(text);
    SendMessageW(tooltip, TTM_ADDTOOLW, 0, reinterpret_cast<LPARAM>(&info));
}

}  // namespace

Workbench::Workbench(WorkbenchArguments arguments)
    : arguments_(std::move(arguments)), backend_(arguments_.backend) {}

bool Workbench::registerClass(HINSTANCE instance) {
    WNDCLASSEXW type{sizeof(type)};
    type.style = CS_HREDRAW | CS_VREDRAW;
    type.lpfnWndProc = procedure;
    type.hInstance = instance;
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hIcon = LoadIconW(instance, MAKEINTRESOURCEW(IDI_DYNSYSLAB_APP));
    type.hIconSm = static_cast<HICON>(LoadImageW(instance,
        MAKEINTRESOURCEW(IDI_DYNSYSLAB_APP), IMAGE_ICON,
        GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), LR_DEFAULTCOLOR));
    type.hbrBackground = CreateSolidBrush(theme::background);
    type.lpszClassName = kClassName;
    return (RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS) &&
           PlotCanvas::registerClass(instance) && InspectorPane::registerClass(instance);
}

int Workbench::run(HINSTANCE instance, int /*showCommand*/) {
    window_ = CreateWindowExW(WS_EX_COMPOSITED, kClassName, L"DynSysLab Viewer", WS_OVERLAPPEDWINDOW |
        WS_CLIPCHILDREN, CW_USEDEFAULT, CW_USEDEFAULT, 1480, 920, nullptr, nullptr,
        instance, this);
    if (!window_) throw std::runtime_error("Unable to create Viewer workbench window");
    theme::configureTitleBar(window_);
    if (!arguments_.backgroundCheck) {
        ShowWindow(window_, SW_SHOWMAXIMIZED);
        UpdateWindow(window_);
    } else {
        RECT workspace{};
        if (SystemParametersInfoW(SPI_GETWORKAREA,0,&workspace,0))
            SetWindowPos(window_,nullptr,workspace.left,workspace.top,
                workspace.right-workspace.left,workspace.bottom-workspace.top,SWP_NOACTIVATE|SWP_NOZORDER);
        layout();
        SetTimer(window_,302,30000,nullptr);
    }
    MSG event{};
    while (GetMessageW(&event, nullptr, 0, 0) > 0) {
        if (event.message == WM_MOUSEWHEEL && WindowFromPoint(event.pt) == plots_.handle())
            event.hwnd = plots_.handle();
        TranslateMessage(&event); DispatchMessageW(&event);
    }
    return static_cast<int>(event.wParam);
}

bool Workbench::createChildren() {
    search_ = CreateWindowExW(WS_EX_CLIENTEDGE, L"EDIT", L"", WS_CHILD | WS_VISIBLE |
        WS_TABSTOP |
        ES_AUTOHSCROLL, 0, 0, 0, 0, window_,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kSearchId)),
        GetModuleHandleW(nullptr), nullptr);
    SendMessageW(search_, EM_SETCUEBANNER, TRUE,
        reinterpret_cast<LPARAM>(L"Search archive tree…"));
    theme::styleControl(search_, fonts_.body);
    refresh_ = button(window_, kRefreshId, L"Refresh archive", fonts_.heading);
    sortDate_ = button(window_, kSortDateId, L"Sort by date", fonts_.heading);
    sortTrajectories_ = button(window_, kSortTrajectoriesId,
        L"Sort by trajectory count", fonts_.heading);
    sortSize_ = button(window_, kSortSizeId, L"Sort by total size", fonts_.heading);
    sortTooltip_ = CreateWindowExW(WS_EX_TOPMOST, TOOLTIPS_CLASSW, nullptr,
        WS_POPUP | TTS_ALWAYSTIP | TTS_NOPREFIX, CW_USEDEFAULT, CW_USEDEFAULT,
        CW_USEDEFAULT, CW_USEDEFAULT, window_, nullptr, GetModuleHandleW(nullptr), nullptr);
    addTooltip(sortTooltip_, sortDate_, L"Sort leaves by date (newest first)");
    addTooltip(sortTooltip_, sortTrajectories_,
        L"Sort leaves by number of trajectories (highest first)");
    addTooltip(sortTooltip_, sortSize_,
        L"Sort leaves by total data size (largest first)");
    addTooltip(sortTooltip_, refresh_, L"Refresh the archive index");
    previous_ = button(window_, kPreviousId, L"‹", fonts_.heading);
    next_ = button(window_, kNextId, L"›", fonts_.heading);
    viewSelect_ = CreateWindowExW(0, L"COMBOBOX", L"", WS_CHILD | WS_VISIBLE |
        WS_TABSTOP |
        CBS_DROPDOWNLIST | WS_VSCROLL, 0, 0, 0, 0, window_,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kViewSelectId)),
        GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(viewSelect_, fonts_.caption);
    pageSelect_ = CreateWindowExW(0, L"COMBOBOX", L"", WS_CHILD | WS_VISIBLE |
        WS_TABSTOP |
        CBS_DROPDOWNLIST | WS_VSCROLL, 0, 0, 0, 0, window_,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(kPageSelectId)),
        GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(pageSelect_, fonts_.caption);
    copy_ = button(window_, kCopyId, L"Copy plot", fonts_.body);
    save_ = button(window_, kSaveId, L"Save PNG/PDF", fonts_.body);
    axes_ = button(window_, kAxesId, L"Axes", fonts_.body);
    addPanel_ = button(window_, kAddPanelId, L"+ Panel", fonts_.body);
    compare_ = button(window_, kCompareId, L"Compare", fonts_.body);
    compute_ = button(window_, kComputeId, L"Post-process", fonts_.body);
    help_ = button(window_, kHelpId, L"Help / About", fonts_.body);
    phaseContext_ = button(window_, kPhaseContextId, L"Phase context", fonts_.body);
    mainTab_ = button(window_, kMainTabId, L"Main", fonts_.heading);
    observablesTab_ = button(window_, kObservablesTabId, L"+ observables", fonts_.heading);
    selectionContext_ = CreateWindowExW(0, L"STATIC", L"Select training data",
        WS_CHILD | WS_VISIBLE | SS_LEFT | SS_CENTERIMAGE | SS_ENDELLIPSIS,
        0, 0, 0, 0, window_, nullptr, GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(selectionContext_, fonts_.heading);
    if (!sortDate_ || !sortTrajectories_ || !sortSize_ || !sortTooltip_ ||
        !selectionContext_ ||
        !navigation_.create(window_, kNavigationId, fonts_) ||
        !inspector_.create(window_, fonts_) || !plots_.create(window_, fonts_) ||
        !plotDragToolbar_.create(window_, plots_, fonts_) ||
        !leftSplitter_.create(window_, kLeftSplitterId)) return false;
    navigation_.setCatalog(&catalog_);
    if (arguments_.backgroundCheck) installBackgroundPaintCheck();
    plots_.setMessage(L"Loading archive catalogue",
        L"Reading the Archiver-maintained index…");
    updateExportButtons();
    updatePageSelector();
    SetTimer(window_, kTaskTimer, 80, nullptr);
    taskActive_ = true;
    task_ = std::async(std::launch::async, [this]() {
        AsyncResult result; result.kind = AsyncResult::Kind::Catalog;
        result.startup = true;
        try { result.catalog = std::make_unique<Catalog>(backend_.loadCatalog()); }
        catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
    return true;
}

void Workbench::setTreeSort(CatalogTreeSort sort) {
    treeSort_ = sort;
    navigation_.setSort(sort);
    if (selected_ >= 0) navigation_.selectIndex(selected_);
    InvalidateRect(sortDate_, nullptr, FALSE);
    InvalidateRect(sortTrajectories_, nullptr, FALSE);
    InvalidateRect(sortSize_, nullptr, FALSE);
}

void Workbench::initialSelection() {
    int index = -1;
    if (!arguments_.initialSource.empty()) {
        std::error_code ignored;
        const auto target = std::filesystem::weakly_canonical(arguments_.initialSource, ignored);
        for (std::size_t item = 0; item < catalog_.items.size(); ++item) {
            if (std::filesystem::equivalent(backend_.resolveSource(catalog_.items[item]), target, ignored)) {
                index = static_cast<int>(item); break;
            }
        }
    }
    if (index < 0)
        index = preferredInitialCatalogItem(catalog_, treeSort_);
    if (index < 0 && !catalog_.items.empty()) index = 0;
    if (index >= 0) { navigation_.selectIndex(index); selectItem(index); }
    navigation_.resetHorizontalScroll();
}

void Workbench::paint(HDC target) {
    RECT area{}; GetClientRect(window_, &area);
    const auto ui = theme::metrics(window_);
    HBRUSH background = CreateSolidBrush(theme::background);
    FillRect(target, &area, background);
    DeleteObject(background);
    SetBkMode(target, TRANSPARENT);
    SelectObject(target, fonts_.title);
    SetTextColor(target, theme::text);
    TextOutW(target, theme::scalePixels(window_, 16),
        theme::scalePixels(window_, 10), L"DynSysLab Viewer", 16);
    HPEN line = CreatePen(PS_SOLID, 1, theme::border);
    HGDIOBJ previous = SelectObject(target, line);
    MoveToEx(target, 0, ui.headerHeight - 1, nullptr);
    LineTo(target, area.right, ui.headerHeight - 1);
    SelectObject(target, previous); DeleteObject(line);
}

void Workbench::drawButton(const DRAWITEMSTRUCT& item) {
    if (item.CtlID == kMainTabId || item.CtlID == kObservablesTabId) {
        drawFolderTab(item, (item.CtlID == kObservablesTabId) == observablesActive_, fonts_.body);
        return;
    }
    const bool checked = (item.CtlID == kSortDateId && treeSort_ == CatalogTreeSort::date) ||
        (item.CtlID == kSortTrajectoriesId &&
         treeSort_ == CatalogTreeSort::trajectoryCount) ||
        (item.CtlID == kSortSizeId && treeSort_ == CatalogTreeSort::totalBytes);
    const bool pressed = (item.itemState & ODS_SELECTED) != 0;
    const bool enabled = (item.itemState & ODS_DISABLED) == 0;
    const COLORREF fill = !enabled ? theme::controlBackground
        : pressed ? theme::accent : checked ? theme::selection : theme::panelRaised;
    HBRUSH brush = CreateSolidBrush(fill); FillRect(item.hDC, &item.rcItem, brush); DeleteObject(brush);
    HBRUSH outline = CreateSolidBrush(checked ? theme::accent : theme::border);
    FrameRect(item.hDC, &item.rcItem, outline);
    DeleteObject(outline);
    SetBkMode(item.hDC, TRANSPARENT);
    SetTextColor(item.hDC, !enabled ? theme::disabled
        : pressed ? theme::panelRaised : checked ? theme::accentStrong : theme::text);
    SelectObject(item.hDC, fonts_.caption);
    std::wstring label = buttonText(item.hwndItem);
    if (item.CtlID == kRefreshId) label = L"\u21BB";
    else if (item.CtlID == kSortDateId) label = L"\u25F7";
    else if (item.CtlID == kSortTrajectoriesId) label = L"\u223F";
    else if (item.CtlID == kSortSizeId) label = L"\u25A4";
    DrawTextW(item.hDC, label.c_str(), -1, const_cast<RECT*>(&item.rcItem),
              DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    if (item.itemState & ODS_FOCUS) {
        RECT focusBounds = item.rcItem;
        InflateRect(&focusBounds, -theme::scalePixels(window_, 3),
            -theme::scalePixels(window_, 3));
        HPEN focusPen = CreatePen(PS_SOLID, theme::scalePixels(window_, 2), theme::focus);
        HGDIOBJ oldPen = SelectObject(item.hDC, focusPen);
        HGDIOBJ oldBrush = SelectObject(item.hDC, GetStockObject(NULL_BRUSH));
        Rectangle(item.hDC, focusBounds.left, focusBounds.top,
            focusBounds.right, focusBounds.bottom);
        SelectObject(item.hDC, oldBrush); SelectObject(item.hDC, oldPen);
        DeleteObject(focusPen);
    }
}

void Workbench::updateExportButtons() {
    plotDragToolbar_.updateCoordinates();
    updateWorkspaceTabs();
    const BOOL enabled = plots_.hasPayload() ? TRUE : FALSE;
    SetWindowTextW(phaseContext_,phaseCancellation_ ? L"Cancel context" : L"Phase context");
    EnableWindow(phaseContext_,phaseCancellation_ || (!taskActive_ && !observablesActive_ && !plots_.isKrom() && !mainSources_.empty() && plots_.phaseSpaceDomain().has_value()));
    EnableWindow(copy_, enabled); EnableWindow(save_, enabled); EnableWindow(axes_, enabled);
    EnableWindow(addPanel_, selected_ >= 0 && !taskActive_ && !observablesActive_ && !plots_.isKrom() ? TRUE : FALSE);
    bool comparable = false;
    if (enabled && selected_ >= 0 && !activeView_.empty())
        for (std::size_t index = 0; index < catalog_.items.size(); ++index)
            if (static_cast<int>(index) != selected_ && std::find(
                catalog_.items[index].views.begin(), catalog_.items[index].views.end(),
                activeView_) != catalog_.items[index].views.end()) { comparable = true; break; }
    EnableWindow(compare_, comparable && !taskActive_ && !observablesActive_ && !plots_.isKrom() ? TRUE : FALSE);
}

void Workbench::updatePageSelector() {
    SendMessageW(pageSelect_, CB_RESETCONTENT, 0, 0);
    for (int page = 0; page < plots_.pageCount(); ++page) {
        const auto label = plots_.pageLabel(page);
        SendMessageW(pageSelect_, CB_ADDSTRING, 0,
                     reinterpret_cast<LPARAM>(label.c_str()));
    }
    SendMessageW(pageSelect_, CB_SETCURSEL, plots_.page(), 0);
    EnableWindow(pageSelect_, plots_.hasPayload() && plots_.pageCount() > 1);
    EnableWindow(previous_, plots_.hasPayload() && plots_.pageCount() > 1);
    EnableWindow(next_, plots_.hasPayload() && plots_.pageCount() > 1);
}

void Workbench::selectItem(int index, bool loadDefault) {
    if (index < 0 || static_cast<std::size_t>(index) >= catalog_.items.size()) return;
    if (index != selected_) resetWorkspace();
    selected_ = index;
    const CatalogItem& item = catalog_.items[static_cast<std::size_t>(index)];
    const std::wstring context = item.system.empty()
        ? wide(item.label) : wide(item.system) + L"  /  " + wide(item.label);
    SetWindowTextW(selectionContext_, context.c_str());
    inspector_.setItem(&catalog_, &item);
    updateViewSelector();
    InvalidateRect(window_, nullptr, FALSE);
    if (loadDefault && !availableViews_.empty()) requestView(availableViews_.front());
    else if (loadDefault) {
        activeView_.clear();
        plots_.setMessage(L"No native scientific view is available",
            L"Viewer does not hand archive objects to Windows file associations.");
        updateExportButtons();
        updatePageSelector();
    }
}

void Workbench::updateViewSelector() {
    SendMessageW(viewSelect_, CB_RESETCONTENT, 0, 0);
    availableViews_.clear();
    if (selected_ < 0) { EnableWindow(viewSelect_, FALSE); return; }
    const CatalogItem& item = catalog_.items[static_cast<std::size_t>(selected_)];
    int selectedView = 0;
    for (const auto& view : item.views) {
        if (view == "open") continue;
        const int index = static_cast<int>(availableViews_.size());
        availableViews_.push_back(view);
        const auto found = catalog_.viewLabels.find(view);
        const std::string label = item.kind == "krom" && (view == "koopman_spectrum" || view == "krom_dashboard")
            ? "KROM workspace" : found == catalog_.viewLabels.end() ? view : found->second;
        SendMessageW(viewSelect_, CB_ADDSTRING, 0,
                     reinterpret_cast<LPARAM>(wide(label).c_str()));
        if (view == activeView_) selectedView = index;
    }
    SendMessageW(viewSelect_, CB_SETCURSEL, selectedView, 0);
    EnableWindow(viewSelect_, availableViews_.empty() ? FALSE : TRUE);
}

LRESULT Workbench::message(UINT messageValue, WPARAM wparam, LPARAM lparam) {
    if (messageValue == WM_CREATE) return createChildren() ? 0 : -1;
    if (messageValue == WM_SIZE) { layout(); return 0; }
    if (messageValue == WM_PAINT) {
        PAINTSTRUCT state{}; HDC target = BeginPaint(window_, &state);
        paint(target); EndPaint(window_, &state); return 0;
    }
    if (messageValue == WM_ERASEBKGND) return 1;
    if (messageValue == WM_DRAWITEM) { drawButton(*reinterpret_cast<DRAWITEMSTRUCT*>(lparam)); return TRUE; }
    if (messageValue == kVerticalSplitterMovedMessage) {
        if (static_cast<int>(wparam) == kLeftSplitterId)
            navigationWidth_ = std::max(210,
                navigationWidthFromSplitter(static_cast<int>(lparam), 12, 12)
                - 12 - inspectorWidth_);
        layout();
        return 0;
    }
    if (messageValue == WM_TIMER && wparam == kTaskTimer) { pollTask(); return 0; }
    if (messageValue == WM_TIMER && wparam == 302 && arguments_.backgroundCheck) {
        finishBackgroundCheck(2,L"Desktop startup timed out"); return 0;
    }
    if (messageValue == WM_NOTIFY && reinterpret_cast<NMHDR*>(lparam)->idFrom == kNavigationId &&
        reinterpret_cast<NMHDR*>(lparam)->code == TVN_SELCHANGEDW) {
        if (updatingCatalog_) return 0;
        const int index = navigation_.selectedIndex(); if (index >= 0) selectItem(index); return 0;
    }
    if (messageValue == kPlotPageChangedMessage) {
        SendMessageW(pageSelect_, CB_SETCURSEL, plots_.page(), 0); return 0;
    }
    if (messageValue == kPlotKromChangedMessage) {
        const auto state = plots_.captureState();
        if (state.payload && plots_.isKrom()) {
            KromWorkspaceOptions options; options.modes = state.payload->selectedModes;
            options.complexDisplay = state.payload->complexDisplay; requestKrom(options);
        }
        return 0;
    }
    if (messageValue == kPlotCompositionChangedMessage) {
        updateExportButtons(); updatePageSelector(); return 0;
    }
    if (messageValue == WM_COMMAND) {
        const int id = LOWORD(wparam);
        if (id == kSearchId && HIWORD(wparam) == EN_CHANGE) {
            const int count = GetWindowTextLengthW(search_);
            std::wstring value(static_cast<std::size_t>(count + 1), L'\0');
            GetWindowTextW(search_, value.data(), count + 1);
            value.resize(static_cast<std::size_t>(count));
            navigation_.setFilter(value); return 0;
        }
        if (id == kMainTabId) { activateWorkspace(false); return 0; }
        if (id == kSortDateId) { setTreeSort(CatalogTreeSort::date); return 0; }
        if (id == kSortTrajectoriesId) {
            setTreeSort(CatalogTreeSort::trajectoryCount); return 0;
        }
        if (id == kSortSizeId) { setTreeSort(CatalogTreeSort::totalBytes); return 0; }
        if (id == kObservablesTabId) { activateWorkspace(true); return 0; }
        if (id == kPhaseContextId) { requestPhaseContext(); return 0; }
        if (id == kRefreshId) { refreshCatalog(); return 0; }
        if (id == kPreviousId) { plots_.previousPage(); updatePageSelector(); return 0; }
        if (id == kNextId) { plots_.nextPage(); updatePageSelector(); return 0; }
        if (id == kPageSelectId && HIWORD(wparam) == CBN_SELCHANGE) {
            plots_.setPage(static_cast<int>(SendMessageW(pageSelect_, CB_GETCURSEL, 0, 0)));
            return 0;
        }
        if (id == kViewSelectId && HIWORD(wparam) == CBN_SELCHANGE && selected_ >= 0) {
            const int view = static_cast<int>(SendMessageW(viewSelect_, CB_GETCURSEL, 0, 0));
            if (view >= 0 && static_cast<std::size_t>(view) < availableViews_.size())
                requestView(availableViews_[static_cast<std::size_t>(view)]);
            return 0;
        }
        if (id == kCopyId) { copyPlots(); return 0; }
        if (id == kSaveId) { savePlots(); return 0; }
        if (id == kAxesId) {
            (void)showAxisRangeDialog(window_, plots_, fonts_);
            return 0;
        }
        if (id == kAddPanelId) { startAddPanel(); return 0; }
        if (id == kCompareId) { startComparison(); return 0; }
        if (id == kComputeId) { openPostprocessor(); return 0; }
        if (id == kHelpId) {
            showHelpAbout(window_, L"DynSysLab Viewer", arguments_.backend.projectRoot,
                arguments_.backend.archiveRoot); return 0;
        }
    }
    if (messageValue == kInspectorRelationMessage && selected_ >= 0) {
        const auto& relations = catalog_.items[static_cast<std::size_t>(selected_)].relations;
        const std::size_t relation = static_cast<std::size_t>(wparam);
        if (relation < relations.size()) {
            const int index = catalog_.indexOf(relations[relation].target);
            if (index >= 0) { navigation_.selectIndex(index); selectItem(index); }
        }
        return 0;
    }
    if (messageValue == WM_CTLCOLORSTATIC) {
        HDC dc = reinterpret_cast<HDC>(wparam); SetBkColor(dc, theme::background);
        SetTextColor(dc, theme::muted); return reinterpret_cast<LRESULT>(GetStockObject(NULL_BRUSH));
    }
    if (messageValue == WM_CTLCOLOREDIT) {
        static HBRUSH editBackground = CreateSolidBrush(theme::panelRaised);
        HDC dc = reinterpret_cast<HDC>(wparam); SetBkColor(dc, theme::panelRaised);
        SetTextColor(dc, theme::text); return reinterpret_cast<LRESULT>(editBackground);
    }
    if (messageValue == WM_DESTROY) { if(phaseCancellation_) phaseCancellation_->store(true); KillTimer(window_, kTaskTimer); PostQuitMessage(0); return 0; }
    return DefWindowProcW(window_, messageValue, wparam, lparam);
}

LRESULT CALLBACK Workbench::procedure(HWND window, UINT messageValue,
                                      WPARAM wparam, LPARAM lparam) {
    auto* self = reinterpret_cast<Workbench*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (messageValue == WM_NCCREATE) {
        self = static_cast<Workbench*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        self->window_ = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    }
    return self ? self->message(messageValue, wparam, lparam) :
                  DefWindowProcW(window, messageValue, wparam, lparam);
}

WorkbenchArguments parseArguments() {
    int count = 0;
    wchar_t** values = CommandLineToArgvW(GetCommandLineW(), &count);
    WorkbenchArguments result;
    for (int index = 1; index < count; ++index) {
        const std::wstring key = values[index];
        if (key == L"--smoke-test") { result.smokeTest = true; continue; }
        if (key == L"--background-check") { result.backgroundCheck = true; continue; }
        if (key == L"--background-krom-check") { result.backgroundCheck = result.backgroundKromCheck = true; continue; }
        if (key == L"--background-workspace-check") {
            result.backgroundCheck = result.backgroundWorkspaceCheck = true; continue;
        }
        if (index + 1 >= count) break;
        const std::filesystem::path value = values[++index];
        if (key == L"--project-root") result.backend.projectRoot = value;
        else if (key == L"--archive-root") result.backend.archiveRoot = value;
        else if (key == L"--catalog") result.backend.catalog = value;
        else if (key == L"--cache") result.backend.cache = value;
        else if (key == L"--archive-focus") result.backend.archiveFocus = value;
        else if (key == L"--source") result.initialSource = value;
        else if (key == L"--background-refresh-count") {
            result.backgroundRefreshCount = std::stoul(value.wstring());
            if (result.backgroundRefreshCount > 20)
                throw std::runtime_error("Background refresh count must be at most 20");
        }
        else --index;
    }
    LocalFree(values);
    const auto paths = discoverRuntimePaths(L"Viewer", result.backend.projectRoot,
        result.backend.archiveRoot, result.backend.cache);
    result.backend.projectRoot = paths.projectRoot;
    result.backend.archiveRoot = paths.archiveRoot;
    result.backend.cache = paths.cacheRoot;
    return result;
}

}  // namespace dsl
