#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_compute/SolverExecutionPane.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace dsl::compute {

constexpr int kSolverPlanId = 5201;
constexpr int kSolverRunId = 5202;
constexpr int kSolverSystemId = 5210;
constexpr int kSolverPresetId = 5211;
constexpr int kSolverControllerId = 5213;
constexpr int kSolverCampaignId = 5214;
constexpr int kSolverInitialConditionsId = 5215;
constexpr int kSolverControllerSettingsId = 5216;
constexpr int kSolverTrajectoryCountId = 5219;
constexpr int kSolverClusterCountId = 5220;
constexpr int kSolverLongParametersId = 5221;
constexpr int kSolverWaveformCoordinatesId = 5222;
constexpr int kSolverWaveformFrequenciesId = 5223;
constexpr int kSolverWaveformAmplitudesId = 5224;
constexpr int kSolverWaveformPhasesId = 5225;
constexpr int kSolverWaveformOffsetId = 5226;
constexpr int kSolverWaveformAddId = 5227;
constexpr int kSolverWaveformRemoveId = 5228;
constexpr int kSolverTargetKindId = 5229;
constexpr int kSolverTargetSourceId = 5230;
constexpr int kSolverTargetFrequenciesId = 5231;
constexpr int kSolverTargetAmplitudesId = 5232;
constexpr int kSolverTargetPhasesId = 5233;
constexpr int kSolverLqrAId = 5234;
constexpr int kSolverLqrBId = 5235;
constexpr int kSolverStartTimeId = 5236;
constexpr int kSolverStopTimeId = 5237;
constexpr int kSolverAddControlId = 5238;
constexpr int kSolverSystemPageId = 5239;
constexpr int kSolverControlPageId = 5240;
constexpr int kSolverControlEvaluationId = 5241;
constexpr int kSolverPdePartitionId = 5242;
constexpr int kSolverCoordinatesId = 5243;

enum class SolverStage { system, initialConditions, control, integrator, review };

class SolverPanel {
public:
    bool create(HWND parent, const dsl::theme::Fonts& fonts, const Registry* registry);
    void showEditor(bool visible);
    void showExecution(bool visible, bool enabled);
    void layoutEditor(int x, int y, int width, int height);
    void layoutExecution(int x, int y, int width, int height);
    void paintEditorFrame(HDC dc) const;
    void showJobActions(bool visible);
    void layoutJobActions(int x, int y, int width, int height);
    void command(int id, int notification);
    std::string jobJson(const std::string& label) const;
    void loadJobJson(const std::string& jobJson, bool includeExecution = true);
    void resetToSystem(std::size_t systemIndex =
        std::numeric_limits<std::size_t>::max(), bool includeExecution = true);
    std::string applyExecutionContract(const std::string& jobJson) const;
    std::size_t systemIndex() const;
    std::string systemLabel() const;
    std::size_t trajectoryCount() const;
    std::size_t clusterCount() const;
    void setHierarchyCounts(std::size_t trajectories, std::size_t clusters);
    void refreshRegistry();
    bool tabButtonActive(int id) const;
    bool prepareKromTrainingDomainExpansion(
        const std::filesystem::path& preferredKrom = {});
    void activateStage(SolverStage stage);
private:
    const Registry* registry_ = nullptr;
    const dsl::theme::Fonts* fonts_ = nullptr;
    HWND parent_ = nullptr;
    std::vector<std::size_t> visibleKroms_;
    std::vector<std::size_t> visibleTargetTrajectories_;
    std::vector<HWND> controls_;
    std::vector<HWND> labels_;
    HWND systemSection_ = nullptr; HWND controllerSection_ = nullptr;
    HWND systemSeparator_ = nullptr;
    HWND controllerSeparator_ = nullptr;
    HWND parameterSummary_ = nullptr;
    HWND longParameters_ = nullptr;
    HWND system_ = nullptr; HWND preset_ = nullptr; HWND coordinates_ = nullptr;
    HWND parameters_ = nullptr;
    HWND start_ = nullptr; HWND stop_ = nullptr; HWND step_ = nullptr;
    HWND timeLabelDash_ = nullptr; HWND timeRangeColon_ = nullptr;
    HWND timeValueDash_ = nullptr; HWND timeStepColon_ = nullptr;
    HWND pdePartition_ = nullptr; HWND pdePartitionShape_ = nullptr;
    HWND pdeGhostWidth_ = nullptr; HWND pdeLocalWorkers_ = nullptr;
    HWND pdePartitionShapeLabel_ = nullptr; HWND pdeGhostWidthLabel_ = nullptr;
    HWND pdeLocalWorkersLabel_ = nullptr;
    HWND icMode_ = nullptr; HWND count_ = nullptr; HWND clusters_ = nullptr; HWND seed_ = nullptr;
    HWND icData_ = nullptr; HWND icButton_ = nullptr; HWND controller_ = nullptr; HWND controllerSettings_ = nullptr;
    HWND controllerSettingsButton_ = nullptr;
    HWND target_ = nullptr; HWND targetState_ = nullptr; HWND krom_ = nullptr;
    HWND q_ = nullptr; HWND r_ = nullptr; HWND label_ = nullptr;
    HWND representation_ = nullptr;
    HWND controlEvaluation_ = nullptr; HWND controlSampleTime_ = nullptr;
    HWND horizon_ = nullptr; HWND controlLimit_ = nullptr;
    HWND controllerEmpty_ = nullptr; HWND targetSection_ = nullptr;
    HWND waveformPreview_ = nullptr; HWND waveformCoordinates_ = nullptr;
    HWND waveformFrequencies_ = nullptr; HWND waveformAmplitudes_ = nullptr;
    HWND waveformPhases_ = nullptr; HWND waveformOffset_ = nullptr;
    HWND waveformAdd_ = nullptr; HWND waveformRemove_ = nullptr;
    HWND lqrA_ = nullptr; HWND lqrB_ = nullptr;
    HWND targetKind_ = nullptr; HWND targetSource_ = nullptr;
    HWND targetFrequencies_ = nullptr; HWND targetAmplitudes_ = nullptr;
    HWND targetPhases_ = nullptr;
    HWND waveformCoordinatesLabel_ = nullptr; HWND waveformFrequenciesLabel_ = nullptr;
    HWND waveformAmplitudesLabel_ = nullptr; HWND waveformPhasesLabel_ = nullptr;
    HWND waveformOffsetLabel_ = nullptr; HWND lqrALabel_ = nullptr;
    HWND lqrBLabel_ = nullptr; HWND targetKindLabel_ = nullptr;
    HWND targetSourceLabel_ = nullptr; HWND targetFrequenciesLabel_ = nullptr;
    HWND targetAmplitudesLabel_ = nullptr; HWND targetPhasesLabel_ = nullptr;
    HWND controlEvaluationLabel_ = nullptr;
    HWND controlSampleTimeLabel_ = nullptr;
    int layoutX_{}, layoutY_{}, layoutWidth_{}, layoutHeight_{};
    RECT editorBounds_{};
    int executionX_{}, executionY_{}, executionWidth_{}, executionHeight_{};
    bool editorVisible_{};
    bool executionVisible_{};
    bool executionEnabled_{};
    bool controlPageActive_{};
    SolverExecutionPane execution_;
    HWND addControl_ = nullptr; HWND systemPage_ = nullptr;
    HWND controlPage_ = nullptr;
    HWND plan_ = nullptr; HWND run_ = nullptr;
    std::string preflight_{"{\"enabled\":false}"};
    std::string archiveDataRole_{"free"};
    std::string clusterSpecs_{"[]"};
    std::string campaign_{"{\"enabled\":false,\"mode\":\"cartesian\",\"axes\":[]}"};
    void updatePresets(bool includeExecution = true);
    void updateParameterSummaries();
    void setParameterOverride(const std::string& key, double value);
    void updateTimeLabels();
    void updateSystemChoices();
    void updateControllerControls();
    void changeControllerKind();
    bool hasControl() const;
    void setPanelPage(bool control);
    void updatePanelPageVisibility();
    void updateWaveformControls();
    void updateWaveformSettings();
    void adjustWaveformComponents(bool add);
    void updateTargetControls();
    void layoutControllerArea(int x, int y, int width, int height);
    void editCampaign();
    void editInitialConditions();
    void editControllerSettings();
    void clearTrainingPreflight();
    std::vector<double> resolvedInitialState() const;
    std::vector<std::string> campaignSuggestions() const;
    std::string controlSettings() const;
    std::string targetSettings() const;
    std::string targetSettingsBase() const;
    std::string solverSettings() const;
    std::string pdePartitioningSettings() const;
    void add(HWND value) { controls_.push_back(value); }
};

}  // namespace dsl::compute
