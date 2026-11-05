#pragma once

#include "dynsyslab_compute/Backend.hpp"
#include "dynsyslab_compute/PostprocessPanel.hpp"
#include "dynsyslab_compute/SolverPanel.hpp"
#include "dynsyslab_compute/SolveJobPane.hpp"
#include "dynsyslab_compute/SolveContextPanel.hpp"
#include "dynsyslab_viewer/HorizontalPaneLayout.hpp"
#include "dynsyslab_viewer/InspectorPane.hpp"
#include "dynsyslab_viewer/NavigationPane.hpp"
#include "dynsyslab_viewer/Theme.hpp"
#include "dynsyslab_viewer/VerticalSplitter.hpp"

#include <future>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>

namespace dsl::compute {

struct Arguments {
    BackendConfig backend;
    std::filesystem::path initialSource;
    bool openPostprocessor = false;
    bool smokeTest = false;
    bool backgroundCheck = false;
};

class Workbench {
public:
    explicit Workbench(Arguments arguments);
    static bool registerClass(HINSTANCE instance);
    int run(HINSTANCE instance, int showCommand);
private:
    struct TaskResult { std::wstring text; std::wstring error;
        std::unique_ptr<dsl::Catalog> catalog; bool refresh = false;
        bool startup = false; bool cancelled = false; };
    Arguments arguments_; Backend backend_; dsl::Catalog catalog_; Registry registry_;
    dsl::theme::Fonts fonts_; dsl::NavigationPane navigation_; dsl::InspectorPane inspector_;
    SolveJobPane solveJobs_;
    SolveContextPanel solveContext_;
    dsl::VerticalSplitter leftSplitter_; dsl::VerticalSplitter rightSplitter_;
    SolverPanel solver_; PostprocessPanel post_;
    HWND window_ = nullptr; HWND search_ = nullptr; HWND refresh_ = nullptr;
    HWND solverTab_ = nullptr; HWND postTab_ = nullptr; HWND viewer_ = nullptr;
    HWND stageSystem_ = nullptr; HWND stageInitialConditions_ = nullptr;
    HWND stageControl_ = nullptr; HWND stageIntegrator_ = nullptr;
    HWND stageReview_ = nullptr;
    HWND saveJob_ = nullptr; HWND openJob_ = nullptr; HWND queue_ = nullptr;
    HWND archiveTools_ = nullptr; HWND help_ = nullptr; HWND cancelTask_ = nullptr;
    HWND importData_ = nullptr;
    HWND output_ = nullptr; HWND status_ = nullptr; HWND progress_ = nullptr;
    HWND progressLabel_ = nullptr;
    HWND jobActionsCaption_ = nullptr;
    HWND reviewRailTitle_ = nullptr;
    HWND headerTooltip_ = nullptr;
    int selected_ = -1; bool solverActive_ = true; bool taskActive_ = false;
    bool archiveNavigationLoaded_ = false;
    bool postCreated_ = false;
    std::optional<SolveSystemAddress> activeSolveSystem_;
    int navigationWidth_ = 300; int inspectorWidth_ = 310;
    bool paneWidthsInitialized_ = false;
    SolverStage activeSolverStage_ = SolverStage::system;
    bool progressDeterminate_ = false;
    std::filesystem::path progressPath_;
    std::chrono::steady_clock::time_point lastProgressPoll_{};
    std::future<TaskResult> task_;
    std::atomic_bool cancelRequested_{false};
    bool createChildren(); void layout(); void paint(HDC dc); void setTab(bool solver);
    void ensureArchiveNavigation();
    bool ensurePostPanel(); void selectInitialArchive();
    void drawButton(const DRAWITEMSTRUCT& item);
    void selectArchive(int index); void refreshCatalog(); void runJob(bool solve, bool execute);
    void runJobJson(std::string job, bool execute); void saveJobDefinition();
    void openJobDefinition(); void openQueue(); void openArchiveTools();
    void importLegacyData(); void openKromControlSearch();
    void openKromTrainingDomainPlanner();
    void cancelActiveTask();
    void beginProgress(bool external); void pollProgress(); void endProgress(bool success);
    void pollTask(); void openViewer();
    void saveBackgroundSnapshot(const wchar_t* filename);
    void finishBackgroundCheck(int status, const std::wstring& error = {});
    void storeActiveSolve(); void selectSolveNode(); void editSolveTree(int id);
    void applyClusterInitialConditions();
    std::string activeJobJson();
    LRESULT message(UINT message, WPARAM wparam, LPARAM lparam);
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

Arguments parseArguments();

}  // namespace dsl::compute
