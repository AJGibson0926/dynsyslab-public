#pragma once

#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/KromWorkspace.hpp"
#include "dynsyslab_viewer/HorizontalPaneLayout.hpp"
#include "dynsyslab_viewer/InspectorPane.hpp"
#include "dynsyslab_viewer/NavigationPane.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"
#include "dynsyslab_viewer/PanelCatalog.hpp"
#include "dynsyslab_viewer/Theme.hpp"
#include "dynsyslab_viewer/VerticalSplitter.hpp"

#include <windows.h>

#include <future>
#include <memory>
#include <string>
#include <vector>

namespace dsl {

struct WorkbenchArguments {
    BackendConfig backend;
    std::filesystem::path initialSource;
    bool smokeTest = false;
    bool backgroundCheck = false;
    bool backgroundWorkspaceCheck = false;
    bool backgroundKromCheck = false;
    unsigned backgroundRefreshCount = 0;
};

class Workbench {
public:
    explicit Workbench(WorkbenchArguments arguments);
    static bool registerClass(HINSTANCE instance);
    int run(HINSTANCE instance, int showCommand);

private:
    struct AsyncResult {
        enum class Kind { Payload, Catalog, PanelChoices } kind = Kind::Payload;
        std::shared_ptr<Payload> payload;
        std::unique_ptr<Catalog> catalog;
        std::vector<PanelChoice> panelChoices;
        std::wstring error;
        bool startup = false;
        int itemIndex = -1;
        std::string view;
        bool observables = false;
        bool phaseContext = false;
        std::vector<CatalogItem> sources;
    };

    WorkbenchArguments arguments_;
    Backend backend_;
    Catalog catalog_;
    theme::Fonts fonts_;
    NavigationPane navigation_;
    InspectorPane inspector_;
    PlotCanvas plots_;
    PlotDragToolbar plotDragToolbar_;
    VerticalSplitter leftSplitter_;
    HWND window_ = nullptr;
    HWND search_ = nullptr;
    HWND refresh_ = nullptr;
    HWND sortDate_ = nullptr;
    HWND sortTrajectories_ = nullptr;
    HWND sortSize_ = nullptr;
    HWND sortTooltip_ = nullptr;
    HWND previous_ = nullptr;
    HWND next_ = nullptr;
    HWND viewSelect_ = nullptr;
    HWND pageSelect_ = nullptr;
    HWND copy_ = nullptr;
    HWND save_ = nullptr;
    HWND axes_ = nullptr;
    HWND addPanel_ = nullptr;
    HWND compare_ = nullptr;
    HWND compute_ = nullptr;
    HWND help_ = nullptr;
    HWND mainTab_ = nullptr;
    HWND observablesTab_ = nullptr;
    HWND phaseContext_ = nullptr;
    HWND selectionContext_ = nullptr;
    std::shared_ptr<std::atomic_bool> phaseCancellation_;
    bool observablesActive_ = false;
    bool observablesOpened_ = false;
    PlotCanvasState mainWorkspace_, observablesWorkspace_;
    std::optional<PhaseSpaceDomain> observablesDomain_;
    std::vector<CatalogItem> mainSources_;
    std::string pendingView_;
    std::optional<KromWorkspaceOptions> pendingKrom_;
    std::future<AsyncResult> task_;
    bool taskActive_ = false;
    bool updatingCatalog_ = false;
    unsigned checkedRefreshes_ = 0;
    unsigned checkedWorkspaceStage_ = 0;
    unsigned checkedKromStage_ = 0;
    std::size_t checkedObservableFigures_ = 0;
    PlotCanvasState checkedMainWorkspace_;
    int selected_ = -1;
    int navigationWidth_ = 285;
    int inspectorWidth_ = 300;
    CatalogTreeSort treeSort_ = CatalogTreeSort::date;
    std::string activeView_;
    std::vector<std::string> availableViews_;

    bool createChildren();
    void layout();
    void paint(HDC target);
    void drawButton(const DRAWITEMSTRUCT& item);
    void copyPlots();
    void savePlots();
    void startComparison();
    void openPostprocessor();
    void startAddPanel();
    void requestPhaseContext();
    void updateExportButtons();
    void updateWorkspaceTabs();
    void resetWorkspace();
    void activateWorkspace(bool observables);
    void updatePageSelector();
    void selectItem(int index, bool loadDefault = true);
    void updateViewSelector();
    void requestView(const std::string& view);
    void requestKrom(KromWorkspaceOptions options);
    void refreshCatalog();
    void setTreeSort(CatalogTreeSort sort);
    void pollTask();
    void applyCatalog(Catalog catalog);
    void installBackgroundPaintCheck();
    void saveBackgroundWorkspaceSnapshot(const wchar_t* filename);
    bool advanceKromCheck();
    void finishBackgroundCheck(int status, const std::wstring& error = {});
    void initialSelection();
    LRESULT message(UINT message, WPARAM wparam, LPARAM lparam);
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

WorkbenchArguments parseArguments();

}  // namespace dsl
