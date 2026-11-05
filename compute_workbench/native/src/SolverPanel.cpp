#include "dynsyslab_compute/SolverPanel.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_compute/ScientificParameterGrid.hpp"
#include "dynsyslab_compute/SolverExecutionContract.hpp"
#include "dynsyslab_compute/WaveformPreview.hpp"
#include "SettingsMerge.hpp"
#include "CampaignDialog.hpp"
#include "InitialConditionDialog.hpp"
#include "LatticeInitialConditionDialog.hpp"
#include "KromTrainingDialog.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace dsl::compute {

namespace {
int choice(HWND combo) { return static_cast<int>(SendMessageW(combo, CB_GETCURSEL, 0, 0)); }
std::string lower(std::string value){std::transform(value.begin(),value.end(),value.begin(),
    [](unsigned char item){return static_cast<char>(std::tolower(item));});return value;}
void selectText(HWND combo,const std::string& value){const auto wide=dsl::wide(value);
    const auto index=SendMessageW(combo,CB_FINDSTRINGEXACT,static_cast<WPARAM>(-1),reinterpret_cast<LPARAM>(wide.c_str()));
    if(index!=CB_ERR)SendMessageW(combo,CB_SETCURSEL,index,0);}
std::string arrayLiteral(std::string value,const char* fallback="[]"){
    const auto begin=value.find_first_not_of(" \t\r\n");if(begin==std::string::npos)return fallback;
    if(value[begin]!='[')value='['+value+']';const auto parsed=dsl::json::parse(value);
    if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument("Expected a comma-separated list of values.");
    return dsl::json::dump(parsed);}
}

bool SolverPanel::create(HWND parent, const dsl::theme::Fonts& fonts, const Registry* registry) {
    registry_ = registry; parent_ = parent; fonts_ = &fonts;
    auto caption = [&](const wchar_t* value) { HWND result = makeLabel(parent, value, fonts.caption); labels_.push_back(result); add(result); return result; };
    auto mathCaption = [&](const wchar_t* value) {
        HWND result = makeMathLabel(parent, value, fonts.body, fonts.caption);
        labels_.push_back(result); add(result); return result;
    };
    auto section = [&](const wchar_t* value) {
        HWND result = makeLabel(parent, value, fonts.heading); add(result); return result;
    };
    auto separator = [&]() {
        HWND result = CreateWindowExW(0, L"STATIC", L"", WS_CHILD | SS_ETCHEDHORZ,
            0, 0, 0, 0, parent, nullptr, GetModuleHandleW(nullptr), nullptr);
        add(result); return result;
    };
    systemSection_ = section(L"System and parameters");
    controllerSection_ = section(L"Controller");
    systemSeparator_ = separator(); controllerSeparator_ = separator();
    addControl_ = makeButton(parent, kSolverAddControlId, L"+ add control", fonts.caption);
    systemPage_ = makeButton(parent, kSolverSystemPageId, L"System", fonts.caption);
    controlPage_ = makeButton(parent, kSolverControlPageId, L"Control", fonts.caption);
    add(addControl_); add(systemPage_); add(controlPage_);
    parameterSummary_ = makeScientificParameterGrid(parent, fonts,
        [this](const std::string& key, double value) { setParameterOverride(key, value); });
    add(parameterSummary_);
    longParameters_ = makeCheck(parent, kSolverLongParametersId,
        L"long", fonts.caption, false); add(longParameters_);
    caption(L"Job label"); label_ = makeEdit(parent, 0, L"New Job", fonts.body); add(label_);
    caption(L"Dynamical system"); system_ = makeCombo(parent, kSolverSystemId, fonts.body); add(system_);
    caption(L"Parameter preset"); preset_ = makeCombo(parent, kSolverPresetId, fonts.body); add(preset_);
    caption(L"Parameter overrides"); parameters_ = makeEdit(parent, 0, L"{}", fonts.mono, true); add(parameters_);
    SendMessageW(parameters_, EM_SETREADONLY, TRUE, 0);
    caption(L"Coordinates"); coordinates_ = makeCombo(parent, kSolverCoordinatesId, fonts.body); add(coordinates_);
    mathCaption(L"t_{0}");
    start_ = makeScientificInlineValue(parent, kSolverStartTimeId, L"0", fonts.body); add(start_);
    mathCaption(L"t_{final}");
    stop_ = makeScientificInlineValue(parent, kSolverStopTimeId, L"10", fonts.body); add(stop_);
    mathCaption(L"\\Delta t");
    step_ = makeScientificInlineValue(parent, 0, L"0.01", fonts.body); add(step_);
    timeLabelDash_ = makeCenteredLabel(parent, L"−", fonts.body);
    timeRangeColon_ = makeCenteredLabel(parent, L":", fonts.body);
    timeValueDash_ = makeCenteredLabel(parent, L"−", fonts.body);
    timeStepColon_ = makeCenteredLabel(parent, L":", fonts.body);
    for (HWND item : {timeLabelDash_, timeRangeColon_, timeValueDash_,
         timeStepColon_}) add(item);
    pdePartition_ = makeCheck(parent, kSolverPdePartitionId,
        L"PDE partitioning", fonts.body, false); add(pdePartition_);
    pdePartitionShapeLabel_ = makeLabel(parent, L"Partition shape", fonts.caption);
    pdePartitionShape_ = makeInlineValue(parent, 0, L"1", fonts.body);
    pdeGhostWidthLabel_ = makeLabel(parent, L"Ghost width", fonts.caption);
    pdeGhostWidth_ = makeInlineValue(parent, 0, L"1", fonts.body);
    pdeLocalWorkersLabel_ = makeLabel(parent, L"Workers", fonts.caption);
    pdeLocalWorkers_ = makeInlineValue(parent, 0, L"1", fonts.body);
    for (HWND item : {pdePartitionShapeLabel_, pdePartitionShape_,
         pdeGhostWidthLabel_, pdeGhostWidth_, pdeLocalWorkersLabel_,
         pdeLocalWorkers_}) add(item);
    caption(L"IC mode"); icMode_ = makeCombo(parent, 0, fonts.body); add(icMode_);
    caption(L"Trajectories"); count_ = makeEdit(parent, kSolverTrajectoryCountId, L"1", fonts.body); add(count_);
    caption(L"Clusters"); clusters_ = makeEdit(parent, kSolverClusterCountId, L"1", fonts.body); add(clusters_);
    caption(L"Seed"); seed_ = makeEdit(parent, 0, L"0", fonts.body); add(seed_);
    caption(L"IC values or per-state ranges (JSON)"); icData_ = makeEdit(parent, 0, L"[]", fonts.mono, true); add(icData_);
    icButton_=makeButton(parent,kSolverInitialConditionsId,L"Edit ICs…",fonts.body);add(icButton_);
    caption(L"Controller"); controller_ = makeCombo(parent, kSolverControllerId, fonts.body); add(controller_);
    caption(L"KROM (KLQR/KMPC)"); krom_ = makeCombo(parent, 0, fonts.body); add(krom_);
    caption(L"Controller settings (JSON; waveform/pulse physical units are converted when available)");
    controllerSettings_ = makeEdit(parent, 0, L"{}", fonts.mono, true); add(controllerSettings_);
    controllerSettingsButton_=makeButton(parent,kSolverControllerSettingsId,L"Configure…",fonts.body);add(controllerSettingsButton_);
    caption(L"Target"); target_ = makeCombo(parent, 0, fonts.body); add(target_);
    caption(L"x*"); targetState_ = makeInlineValue(parent, 0, L"0, 0", fonts.mono); add(targetState_);
    caption(L"Q"); q_ = makeInlineValue(parent, 0, L"1.0", fonts.mono); add(q_);
    caption(L"R"); r_ = makeInlineValue(parent, 0, L"1.0", fonts.mono); add(r_);
    caption(L"Bψ map"); representation_ = makeCombo(parent, 0, fonts.body); add(representation_);
    caption(L"N horizon"); horizon_ = makeInlineValue(parent, 0, L"20", fonts.body); add(horizon_);
    caption(L"u max"); controlLimit_ = makeInlineValue(parent, 0, L"", fonts.mono); add(controlLimit_);
    controlEvaluationLabel_ = makeLabel(parent, L"Controller evaluation", fonts.caption);
    controlEvaluation_ = makeCombo(parent, kSolverControlEvaluationId, fonts.body);
    controlSampleTimeLabel_ = makeLabel(parent, L"Δt control", fonts.caption);
    controlSampleTime_ = makeInlineValue(parent, 0, L"", fonts.body);
    for (HWND item : {controlEvaluationLabel_, controlEvaluation_,
         controlSampleTimeLabel_, controlSampleTime_}) add(item);
    controllerEmpty_ = makeLabel(parent,
        L"Free dynamics: no external control law is attached to this System.",
        fonts.body); add(controllerEmpty_);
    targetSection_ = section(L"Target");
    waveformPreview_ = makeWaveformPreview(parent, fonts); add(waveformPreview_);
    waveformCoordinatesLabel_ = caption(L"coordinates");
    waveformCoordinates_ = makeCombo(parent, kSolverWaveformCoordinatesId, fonts.body);
    add(waveformCoordinates_);
    waveformFrequenciesLabel_ = caption(L"f");
    waveformFrequencies_ = makeInlineValue(parent, kSolverWaveformFrequenciesId,
        L"1", fonts.mono); add(waveformFrequencies_);
    waveformAmplitudesLabel_ = caption(L"A");
    waveformAmplitudes_ = makeInlineValue(parent, kSolverWaveformAmplitudesId,
        L"1", fonts.mono); add(waveformAmplitudes_);
    waveformPhasesLabel_ = caption(L"φ");
    waveformPhases_ = makeInlineValue(parent, kSolverWaveformPhasesId,
        L"0", fonts.mono); add(waveformPhases_);
    waveformOffsetLabel_ = caption(L"u₀");
    waveformOffset_ = makeInlineValue(parent, kSolverWaveformOffsetId,
        L"0", fonts.mono); add(waveformOffset_);
    waveformAdd_ = makeButton(parent, kSolverWaveformAddId, L"+ tone", fonts.body);
    waveformRemove_ = makeButton(parent, kSolverWaveformRemoveId, L"− last", fonts.body);
    add(waveformAdd_); add(waveformRemove_);
    lqrALabel_ = caption(L"A"); lqrA_ = makeInlineValue(parent, kSolverLqrAId,
        L"", fonts.mono); add(lqrA_);
    lqrBLabel_ = caption(L"B"); lqrB_ = makeInlineValue(parent, kSolverLqrBId,
        L"", fonts.mono); add(lqrB_);
    targetKindLabel_ = caption(L"kind");
    targetKind_ = makeCombo(parent, kSolverTargetKindId, fonts.body); add(targetKind_);
    targetSourceLabel_ = caption(L"trajectory");
    targetSource_ = makeCombo(parent, kSolverTargetSourceId, fonts.body); add(targetSource_);
    targetFrequenciesLabel_ = caption(L"f");
    targetFrequencies_ = makeInlineValue(parent, kSolverTargetFrequenciesId,
        L"1", fonts.mono); add(targetFrequencies_);
    targetAmplitudesLabel_ = caption(L"A");
    targetAmplitudes_ = makeInlineValue(parent, kSolverTargetAmplitudesId,
        L"0.1", fonts.mono); add(targetAmplitudes_);
    targetPhasesLabel_ = caption(L"φ");
    targetPhases_ = makeInlineValue(parent, kSolverTargetPhasesId,
        L"0", fonts.mono); add(targetPhases_);
    plan_ = makeButton(parent, kSolverPlanId, L"Plan / validate", fonts.body);
    run_ = makeButton(parent, kSolverRunId, L"Run & archive", fonts.body);
    if (!execution_.create(parent, fonts)) return false;

    std::vector<std::string> systems; for (const auto& item : registry_->systems) systems.push_back(item.label);
    setComboItems(system_, systems);
    std::vector<std::string> coordinateChoices{"dimensional", "dimensionless"};
    for (const auto& item : registry_->systems)
        if (!item.preferredCoordinates.empty()
            && std::find(coordinateChoices.begin(), coordinateChoices.end(),
                         item.preferredCoordinates) == coordinateChoices.end())
            coordinateChoices.push_back(item.preferredCoordinates);
    setComboItems(coordinates_, coordinateChoices);
    setComboItems(icMode_, {"default", "explicit", "uniform", "normal",
                            "sphere", "ellipsoid"});
    setComboItems(controller_, {"none", "waveform", "random_pulse", "pid",
                                "lqr", "sdre", "mpc", "klqr", "kmpc"});
    setComboItems(representation_, {"exact", "constant", "bilinear"});
    setComboItems(controlEvaluation_, {"continuous RHS", "sampled"});
    setComboItems(waveformCoordinates_, {"physical", "model / dimensionless"});
    setComboItems(targetKind_, {"fixed point", "trajectory", "waveform", "canonical preset",
        "observable fixed point", "observable trajectory", "observable waveform", "observable canonical preset",
        "observable from state fixed point", "observable from state trajectory", "observable from state waveform", "observable from state canonical preset"});
    resetToSystem();
    return !controls_.empty();
}

void SolverPanel::updateTimeLabels() {
    const bool dimensionless = selectedText(coordinates_) == "dimensionless";
    SetWindowTextW(labels_[5], dimensionless ? L"\\tau_{0}" : L"t_{0}");
    SetWindowTextW(labels_[6], dimensionless ? L"\\tau_{f}" : L"t_{f}");
    SetWindowTextW(labels_[7], dimensionless ? L"\\Delta \\tau" : L"\\Delta t");
}

void SolverPanel::updatePresets(bool includeExecution) {
    const int index = choice(system_); std::vector<std::string> values{"model default"};
    if (index >= 0 && static_cast<std::size_t>(index) < registry_->systems.size()) {
        const auto& system = registry_->systems[static_cast<std::size_t>(index)];
        selectText(coordinates_, system.preferredCoordinates.empty()
            ? "dimensional" : system.preferredCoordinates);
        for (const auto& preset : system.presets) values.push_back(preset.key);
        const bool ode = system.family == "ode";
        EnableWindow(controller_, ode);
        if (!ode) SendMessageW(controller_, CB_SETCURSEL, 0, 0);
        SendMessageW(pdePartition_, BM_SETCHECK, BST_UNCHECKED, 0);
        SetWindowTextW(pdePartitionShape_, L"1");
        SetWindowTextW(pdeGhostWidth_, L"1");
        SetWindowTextW(pdeLocalWorkers_, L"1");
        if (includeExecution) execution_.resetForFamily(system.family);
        setComboItems(icMode_,system.family=="lattice"
            ?std::vector<std::string>{"default","uniform","explicit"}
            :std::vector<std::string>{"default","explicit","uniform","normal",
                                      "sphere","ellipsoid"});
    }
    setComboItems(preset_, values, 0);
    updateTimeLabels();
    updateSystemChoices();
    if (!hasControl()) controlPageActive_ = false;
    updatePanelPageVisibility();
}

void SolverPanel::updateSystemChoices() {
    const int index=choice(system_);if(index<0||static_cast<std::size_t>(index)>=registry_->systems.size())return;
    const auto& system=registry_->systems[static_cast<std::size_t>(index)];
    const std::string identifier=lower(system.name),label=lower(system.label);
    visibleKroms_.clear();std::vector<std::string> kroms{"Select compatible KROM"};
    for(std::size_t candidate=0;candidate<registry_->kroms.size();++candidate){
        const std::string owned=lower(registry_->kroms[candidate].system);
        if(owned.empty()||owned==identifier||owned==label){visibleKroms_.push_back(candidate);
            kroms.push_back(registry_->kroms[candidate].label);}}
    setComboItems(krom_,kroms);
    std::vector<std::string> targets;
    if(identifier=="rpe") targets=registry_->targets;
    setComboItems(target_,targets);
    visibleTargetTrajectories_.clear();
    std::vector<std::string> trajectories{"Select compatible trajectory"};
    for(std::size_t candidate=0;candidate<registry_->trajectories.size();++candidate){
        const std::string owned=lower(registry_->trajectories[candidate].system);
        if(owned.empty()||owned==identifier||owned==label){
            visibleTargetTrajectories_.push_back(candidate);
            trajectories.push_back(registry_->trajectories[candidate].label);}}
    setComboItems(targetSource_,trajectories);
    if(targets.empty()&&selectedText(targetKind_).find("canonical preset")!=std::string::npos)
        SendMessageW(targetKind_,CB_SETCURSEL,0,0);
}

void SolverPanel::refreshRegistry() {
    updateSystemChoices();
    updateControllerControls();
}

void SolverPanel::showJobActions(bool visible) {
    ShowWindow(plan_, visible ? SW_SHOW : SW_HIDE);
    ShowWindow(run_, visible ? SW_SHOW : SW_HIDE);
}

void SolverPanel::layoutJobActions(int x, int y, int width, int height) {
    const int gap = scalePixels(parent_, 10);
    const int half = (width - gap) / 2;
    MoveWindow(plan_, x, y, half, height, TRUE);
    MoveWindow(run_, x + half + gap, y, width - half - gap, height, TRUE);
}

void SolverPanel::command(int id, int notification) {
    if(id==kSolverAddControlId&&notification==BN_CLICKED){
        SendMessageW(controller_,CB_SETCURSEL,1,0);controlPageActive_=true;
        changeControllerKind();return;}
    if(id==kSolverSystemPageId&&notification==BN_CLICKED){setPanelPage(false);return;}
    if(id==kSolverControlPageId&&notification==BN_CLICKED){setPanelPage(true);return;}
    if (id == kSolverSystemId && notification == CBN_SELCHANGE) {
        SetWindowTextW(parameters_,L"{}"); clearTrainingPreflight(); updatePresets(false);
        updateParameterSummaries();
    }
    if(id==kSolverPresetId&&notification==CBN_SELCHANGE){
        clearTrainingPreflight();updateParameterSummaries();}
    if(id==kSolverCoordinatesId&&notification==CBN_SELCHANGE)
        updateTimeLabels();
    if (id == kSolverControllerId && notification == CBN_SELCHANGE)
        changeControllerKind();
    if(id==kSolverTargetKindId&&notification==CBN_SELCHANGE)updateTargetControls();
    if(id==kSolverControlEvaluationId&&notification==CBN_SELCHANGE)
        updateControllerControls();
    if(id==kSolverWaveformCoordinatesId&&notification==CBN_SELCHANGE)
        updateWaveformSettings();
    if((id==kSolverWaveformFrequenciesId||id==kSolverWaveformAmplitudesId
        ||id==kSolverWaveformPhasesId||id==kSolverWaveformOffsetId)
        &&notification==EN_CHANGE)updateWaveformSettings();
    if(id==kSolverWaveformAddId&&notification==BN_CLICKED)
        adjustWaveformComponents(true);
    if(id==kSolverWaveformRemoveId&&notification==BN_CLICKED)
        adjustWaveformComponents(false);
    if((id==kSolverStartTimeId||id==kSolverStopTimeId)&&notification==EN_CHANGE
        &&selectedText(controller_)=="waveform")updateWaveformControls();
    if(id==kSolverCampaignId&&notification==BN_CLICKED)editCampaign();
    if(id==kSolverInitialConditionsId&&notification==BN_CLICKED)editInitialConditions();
    if(id==kSolverControllerSettingsId&&notification==BN_CLICKED)editControllerSettings();
    if(id==kSolverLongParametersId&&notification==BN_CLICKED){
        const bool enabled=checked(longParameters_);
        setScientificParameterGridLongPrecision(parameterSummary_,enabled);}
}

std::vector<std::string> SolverPanel::campaignSuggestions() const {
    std::vector<std::string> result{
        "controller.settings.amplitude","controller.settings.frequency_hz",
        "controller.settings.amplitudes.0","controller.settings.frequencies_hz.0",
        "controller.settings.pulse_duration_seconds","controller.settings.seed",
        "controller.target.state.0","controller.target.state.1",
        "time.stop","time.output_step","solver.rtol","solver.atol",
        "initial_conditions.seed","initial_conditions.count","initial_conditions.clusters"};
    const int systemIndex=choice(system_),presetIndex=choice(preset_);
    if(systemIndex>=0&&static_cast<std::size_t>(systemIndex)<registry_->systems.size())try{
        dynsyslab::core::SystemConfiguration requested;requested.identifier=registry_->systems[static_cast<std::size_t>(systemIndex)].name;
        if(presetIndex>0&&static_cast<std::size_t>(presetIndex-1)<registry_->systems[static_cast<std::size_t>(systemIndex)].presets.size())
            requested.preset=registry_->systems[static_cast<std::size_t>(systemIndex)].presets[static_cast<std::size_t>(presetIndex-1)].key;
        for(const auto& parameter:dynsyslab::workflow::resolveSystemDefinition(requested).configuration.parameters)
            result.push_back("system.parameters."+parameter.first);
    }catch(...){}
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());return result;
}

void SolverPanel::editCampaign() {
    const auto result=showCampaignDialog(parent_,campaign_,campaignSuggestions(),*fonts_);
    if(result)campaign_=*result;
}

std::vector<double> SolverPanel::resolvedInitialState() const {
    const int systemIndex=choice(system_),presetIndex=choice(preset_);
    if(systemIndex<0||static_cast<std::size_t>(systemIndex)>=registry_->systems.size())return {};
    const auto& definition=registry_->systems[static_cast<std::size_t>(systemIndex)];
    try{dynsyslab::core::SystemConfiguration requested;requested.identifier=definition.name;
        if(presetIndex>0&&static_cast<std::size_t>(presetIndex-1)<definition.presets.size())requested.preset=definition.presets[static_cast<std::size_t>(presetIndex-1)].key;
        const auto overrides=dsl::json::parse(jsonObject(utf8Text(parameters_)));for(const auto& item:overrides.object())requested.parameters[item.first]=item.second.number();
        const auto system=dynsyslab::core::createSystem(requested);return selectedText(controller_)=="none"?system->initialState():system->controlledInitialState();}
    catch(...){return std::vector<double>(static_cast<std::size_t>(std::max(1,definition.stateDimension)),0.0);}
}

void SolverPanel::editInitialConditions() {
    const int systemIndex=choice(system_);if(systemIndex<0||static_cast<std::size_t>(systemIndex)>=registry_->systems.size())return;
    const auto& definition=registry_->systems[static_cast<std::size_t>(systemIndex)];
    std::size_t count=1;try{count=std::max<std::size_t>(1,std::stoull(utf8Text(count_)));}catch(...){}
    if(definition.family=="lattice")try{dynsyslab::core::SystemConfiguration requested;requested.identifier=definition.name;
        const int presetIndex=choice(preset_);if(presetIndex>0&&static_cast<std::size_t>(presetIndex-1)<definition.presets.size())requested.preset=definition.presets[static_cast<std::size_t>(presetIndex-1)].key;
        const auto overrides=dsl::json::parse(jsonObject(utf8Text(parameters_)));for(const auto& item:overrides.object())requested.parameters[item.first]=item.second.number();const auto resolved=dynsyslab::workflow::resolveSystemDefinition(requested).configuration;
        if(selectedText(icMode_)!="explicit"){SetWindowTextW(icData_,L"[]");MessageBoxW(parent_,L"Default/uniform lattice initial conditions are generated from the selected seed by the compiled engine; no explicit site array is required.",L"Lattice initial conditions",MB_OK|MB_ICONINFORMATION);return;}
        const auto result=showLatticeInitialConditionDialog(parent_,resolved.identifier,static_cast<std::size_t>(resolved.parameters.at("size")),static_cast<std::size_t>(resolved.parameters.count("states")?resolved.parameters.at("states"):2.0),count,utf8Text(icData_),*fonts_);if(result)SetWindowTextW(icData_,dsl::wide(*result).c_str());return;
    }catch(const std::exception& error){MessageBoxW(parent_,dsl::wide(error.what()).c_str(),L"Lattice initial conditions",MB_OK|MB_ICONERROR);return;}
    std::vector<std::string> labels;std::vector<double> defaults;
    if(definition.family=="ode"){labels=definition.stateLabels;defaults=resolvedInitialState();}
    else try{dynsyslab::core::SystemConfiguration requested;requested.identifier=definition.name;const int presetIndex=choice(preset_);if(presetIndex>0&&static_cast<std::size_t>(presetIndex-1)<definition.presets.size())requested.preset=definition.presets[static_cast<std::size_t>(presetIndex-1)].key;const auto overrides=dsl::json::parse(jsonObject(utf8Text(parameters_)));for(const auto& item:overrides.object())requested.parameters[item.first]=item.second.number();const auto resolved=dynsyslab::workflow::resolveSystemDefinition(requested).configuration;
        if(definition.family=="dde"){labels={"history(t <= 0)"};defaults={resolved.parameters.at("history_value")};}
        else if(definition.family=="pde"){defaults=dynsyslab::core::heatEquationInitialField(resolved);const auto coordinates=dynsyslab::core::heatEquationCoordinates(resolved);for(std::size_t index=0;index<defaults.size();++index){std::ostringstream name;name<<"u["<<index<<"] @ x="<<coordinates[index];labels.push_back(name.str());}}
    }catch(const std::exception& error){MessageBoxW(parent_,dsl::wide(error.what()).c_str(),L"Initial-condition shape",MB_OK|MB_ICONERROR);return;}
    if(labels.size()!=defaults.size()){labels.clear();for(std::size_t index=0;index<defaults.size();++index)labels.push_back("state_"+std::to_string(index));}
    const auto result=showInitialConditionDialog(parent_,selectedText(icMode_),count,labels,defaults,utf8Text(icData_),*fonts_);
    if(result){SetWindowTextW(icData_,dsl::wide(*result).c_str());clearTrainingPreflight();}
}

void SolverPanel::clearTrainingPreflight(){preflight_="{\"enabled\":false}";archiveDataRole_="free";
}

bool SolverPanel::prepareKromTrainingDomainExpansion(
    const std::filesystem::path& preferredKrom) {
    std::vector<KromDefinition> compatible = registry_->kroms;
    std::size_t current = 0;
    if (!preferredKrom.empty()) {
        std::error_code requestedError;
        const auto requested = std::filesystem::weakly_canonical(
            preferredKrom, requestedError);
        for (std::size_t index = 0; index < compatible.size(); ++index) {
            std::error_code candidateError;
            const auto candidate = std::filesystem::weakly_canonical(
                std::filesystem::u8path(compatible[index].source), candidateError);
            if (!requestedError && !candidateError && candidate == requested) {
                current = index; break;
            }
        }
    }
    const auto result=showKromTrainingDialog(parent_,compatible,current,*fonts_);
    if(!result)return false;
    for(std::size_t index=0;index<registry_->systems.size();++index)if(registry_->systems[index].name==result->system){SendMessageW(system_,CB_SETCURSEL,index,0);updatePresets(true);break;}
    selectText(preset_,result->preset);SetWindowTextW(parameters_,dsl::wide(result->parametersJson).c_str());selectText(coordinates_,result->coordinateSystem);
    SendMessageW(icMode_,CB_SETCURSEL,1,0);SetWindowTextW(count_,std::to_wstring(result->trajectoryCount).c_str());SetWindowTextW(clusters_,L"1");SetWindowTextW(icData_,dsl::wide(result->initialConditionsJson).c_str());
    SendMessageW(controller_,CB_SETCURSEL,0,0);controlPageActive_=false;
    updatePanelPageVisibility();
    SetWindowTextW(stop_,dsl::wide(result->finalTime).c_str());SetWindowTextW(step_,dsl::wide(result->outputStep).c_str());
    preflight_=result->preflightJson;archiveDataRole_="free_training";
    updateParameterSummaries();return true;}

std::string SolverPanel::solverSettings() const {
    std::vector<std::pair<std::string,std::string>> values{{"control_evaluation",
        selectedText(controlEvaluation_) == "sampled"
            ? "\"sampled\"" : "\"continuous_rhs\""}};
    const std::string sample=utf8Text(controlSampleTime_);
    values.push_back({"control_sample_time",sample.empty()?"null":sample});
    return mergeSettingsObject(execution_.settings(),values);
}

std::string SolverPanel::pdePartitioningSettings() const {
    const auto& system = registry_->systems.at(systemIndex());
    if (system.family != "pde") return "null";
    return std::string{"{\"enabled\":"}
        + (checked(pdePartition_) ? "true" : "false")
        + ",\"partition_shape\":"
        + arrayLiteral(utf8Text(pdePartitionShape_), "[1]")
        + ",\"ghost_width\":" + jsonObject(utf8Text(pdeGhostWidth_), "1")
        + ",\"local_workers\":" + jsonObject(utf8Text(pdeLocalWorkers_), "1")
        + '}';
}

std::string SolverPanel::jobJson(const std::string& label) const {
    const int systemIndex = choice(system_), presetIndex = choice(preset_);
    const auto& system = registry_->systems.at(static_cast<std::size_t>(systemIndex));
    const std::string preset = presetIndex > 0 ? system.presets.at(static_cast<std::size_t>(presetIndex-1)).key : "";
    const std::string mode = selectedText(icMode_); const std::string data = jsonObject(utf8Text(icData_), "[]");
    const int controllerIndex = choice(controller_); const std::string controller = selectedText(controller_);
    const std::string target = targetSettings();
    auto clusterSpecs = dsl::json::parse(clusterSpecs_).array();
    if (!clusterSpecs.empty()) {
        const std::size_t clusters = clusterCount();
        while (clusterSpecs.size() < clusters)
            clusterSpecs.emplace_back();
        clusterSpecs.resize(clusters);
    }
    const std::string clusterSpecsJson = dsl::json::dump(
        dsl::json::Value(std::move(clusterSpecs)));
    std::ostringstream out;
    out << "{\"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"solve\",\"label\":\"" << jsonEscape(label)
        << "\",\"system\":{\"name\":\"" << jsonEscape(system.name) << "\",\"preset\":"
        << (preset.empty() ? "null" : "\"" + jsonEscape(preset) + "\"") << ",\"parameters\":" << jsonObject(utf8Text(parameters_))
        << ",\"family\":\"" << jsonEscape(system.family)
        << "\",\"pde_partitioning\":" << pdePartitioningSettings()
        << "},\"time\":{\"start\":" << utf8Text(start_) << ",\"stop\":" << utf8Text(stop_)
        << ",\"output_step\":" << utf8Text(step_) << ",\"coordinate_system\":\"" << jsonEscape(selectedText(coordinates_))
        << "\"},\"initial_conditions\":{\"mode\":\"" << jsonEscape(mode) << "\",\"count\":" << utf8Text(count_)
        << ",\"clusters\":" << utf8Text(clusters_) << ",\"seed\":" << utf8Text(seed_)
        << ",\"values\":" << (mode == "explicit" ? data : "[]")
        << ",\"ranges\":" << (mode == "uniform" || mode == "normal" ? data : "[]")
        << ",\"distribution\":" << (mode == "sphere" || mode == "ellipsoid" ? data : "null")
        << ",\"cluster_specs\":" << clusterSpecsJson
        << "},\"controller\":{\"kind\":\"" << jsonEscape(controller) << "\",\"settings\":" << controlSettings()
        << ",\"target\":" << target << "},\"campaign\":" << jsonObject(campaign_)
        << ",\"preflight\":" << preflight_ << ",\"solver\":" << solverSettings()
        << ",\"archive\":{\"enabled\":true,\"data_role\":\"" << jsonEscape(archiveDataRole_) << "\"}}";
    (void)controllerIndex; return out.str();
}

std::string SolverPanel::applyExecutionContract(const std::string& jobJson) const {
    return applySolverExecutionContract(jobJson, execution_.settings(),
        jsonObject(campaign_));
}

}  // namespace dsl::compute
