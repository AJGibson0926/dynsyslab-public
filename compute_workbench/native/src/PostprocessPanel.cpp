#include "dynsyslab_compute/PostprocessPanel.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "KoopmanPipelineDialog.hpp"
#include "KromExtractionDialog.hpp"
#include "ControlRegressionDialog.hpp"
#include "EstimationDialog.hpp"
#include "DecompositionDialog.hpp"
#include "DmdControlDialog.hpp"
#include "KromValidationDialog.hpp"
#include "BoundaryEigenfunctionDialog.hpp"
#include "EigenvalueOptimizationDialog.hpp"

#include <algorithm>
#include <sstream>

namespace dsl::compute {

bool PostprocessPanel::create(HWND parent, const dsl::theme::Fonts& fonts, const Registry* registry) {
    parent_=parent;fonts_=&fonts;registry_=registry;
    auto caption = [&](const wchar_t* value) { HWND result = makeLabel(parent, value, fonts.caption); labels_.push_back(result); add(result); };
    caption(L"Job label"); label_ = makeEdit(parent, 0, L"interactive-postprocess", fonts.body); add(label_);
    caption(L"Archived source (select a trajectory/KROM in the shared tree or paste a path)");
    source_ = makeEdit(parent, 0, L"", fonts.mono); add(source_);
    caption(L"Ordered processing stages");
    hdmd_ = makeCheck(parent, 0, L"H-DMD", fonts.body, true); add(hdmd_);
    sindy_ = makeCheck(parent, 0, L"SINDy eigenfunctions", fonts.body, true); add(sindy_);
    havok_ = makeCheck(parent, 0, L"HAVOK", fonts.body); add(havok_);
    dmd_ = makeCheck(parent, 0, L"DMD", fonts.body); add(dmd_);
    dmdControl_ = makeCheck(parent, kDmdControlStageId, L"DMDc", fonts.body); add(dmdControl_);
    dmdControlButton_ = makeButton(parent,kDmdControlSettingsId,L"Configure DMDc…",fonts.body); add(dmdControlButton_);
    edmdControl_ = makeCheck(parent,kEdmdControlStageId,L"EDMDc",fonts.body); add(edmdControl_);
    edmdControlButton_ = makeButton(parent,kEdmdControlSettingsId,L"Configure EDMDc…",fonts.body); add(edmdControlButton_);
    edmd_ = makeCheck(parent, 0, L"EDMD", fonts.body); add(edmd_);
    pod_ = makeCheck(parent, 0, L"POD", fonts.body); add(pod_);
    pca_ = makeCheck(parent, 0, L"PCA", fonts.body); add(pca_);
    svd_ = makeCheck(parent, 0, L"SVD", fonts.body); add(svd_);
    baggedPca_ = makeCheck(parent, 0, L"Bagged PCA", fonts.body); add(baggedPca_);
    robustPca_ = makeCheck(parent, 0, L"Robust PCA", fonts.body); add(robustPca_);
    hankelPod_ = makeCheck(parent, 0, L"Hankel POD", fonts.body); add(hankelPod_);
    controlRegression_ = makeCheck(parent, 0, L"Control-map regression", fonts.body); add(controlRegression_);
    kromExtract_ = makeCheck(parent, kKromExtractId, L"Extract KROM", fonts.body); add(kromExtract_);
    kromValidation_ = makeCheck(parent, 0, L"Validate KROM", fonts.body); add(kromValidation_);
    boundaryEigenfunctions_=makeCheck(parent,0,L"Boundary eigenfunctions",fonts.body);
    add(boundaryEigenfunctions_);
    eigenvalueOptimization_=makeCheck(parent,0,L"Optimize eigenvalues",fonts.body);
    add(eigenvalueOptimization_);
    observability_ = makeCheck(parent, 0, L"Observability", fonts.body); add(observability_);
    klqe_ = makeCheck(parent, 0, L"KLQE estimation", fonts.body); add(klqe_);
    caption(L"H-DMD settings (JSON; choose any scalar observable)");
    hdmdSettings_ = makeEdit(parent, 0,
        L"{\"observable\":\"state_0\",\"matrix_size\":2000,\"sampling_window\":20000,\"sample_stride\":1,\"rank\":0,\"svd_tolerance\":1e-10,\"mode_sort\":\"terminal-amplitude-energy\"}", fonts.mono, true); add(hdmdSettings_);
    caption(L"SINDy settings (JSON; consumes preceding or explicit spectrum)");
    sindySettings_ = makeEdit(parent, 0,
        L"{\"polynomial_degree\":4,\"threshold\":0.05,\"max_iterations\":20,\"normalize\":true,\"reuse_conjugate_fits\":true,\"compute_training_values\":false,\"gamma_method\":\"system_rhs\",\"coordinate_observable\":\"state\",\"archive_storage\":\"coefficient-first\",\"reconstruction_observables\":[\"state\"]}", fonts.mono, true); add(sindySettings_);
    caption(L"HAVOK settings (JSON)");
    havokSettings_ = makeEdit(parent, 0,
        L"{\"delays\":100,\"rank\":15,\"observable\":0}", fonts.mono, true); add(havokSettings_);
    dmdSettings_ = makeEdit(parent, 0, L"{\"rank\":0,\"svd_tolerance\":0,\"exact\":true}", fonts.mono, true);
    ShowWindow(dmdSettings_, SW_HIDE);
    caption(L"EDMD settings (JSON; polynomial or RBF dictionary)");
    edmdSettings_ = makeEdit(parent,0,
        L"{\"dictionary\":\"polynomial\",\"polynomial_degree\":2,\"include_constant\":true,\"regularization\":0}",fonts.mono,true); add(edmdSettings_);
    caption(L"PCA settings (other decomposition stages use Configure…)");
    pcaSettings_ = makeEdit(parent,0,
        L"{\"rank\":0}",fonts.mono,true); add(pcaSettings_);
    caption(L"Control regression (JSON; KROM path, ridge, bounds, grid)");
    controlRegressionSettings_ = makeEdit(parent, 0,
        L"{\"krom\":\"\",\"constant_method\":\"multistep-output\",\"bilinear_method\":\"controlled-residual\",\"reconstruction\":\"state\",\"map_sample_stride\":1,\"ridge\":0,\"conjugate_complete\":false,\"grid_columns\":101,\"grid_rows\":101}", fonts.mono, true);
    add(controlRegressionSettings_);
    controlRegressionButton_=makeButton(parent,kControlRegressionSettingsId,L"Configure…",fonts.body);add(controlRegressionButton_);
    caption(L"KROM extraction (zero-based source mode indices)");
    kromExtractSettings_ = makeEdit(parent, 0,
        L"{\"selected_modes\":[0]}", fonts.mono, true); add(kromExtractSettings_);
    kromModesButton_=makeButton(parent,kKromModesId,L"Select modes…",fonts.body);add(kromModesButton_);
    caption(L"Observability / KLQE (JSON; KROM and named C maps)");
    estimationSettings_ = makeEdit(parent, 0,
        L"{\"krom\":\"\",\"measurement\":\"echo\",\"state_reconstruction\":\"state\",\"process_variance\":1,\"measurement_variance\":0.0001}", fonts.mono, true);
    add(estimationSettings_);
    estimationButton_=makeButton(parent,kEstimationSettingsId,L"Configure…",fonts.body);add(estimationButton_);
    pipelineButton_=makeButton(parent,kKoopmanPipelineId,L"Configure standard H-DMD → SINDy…",fonts.body);add(pipelineButton_);
    trainingDomainButton_=makeButton(parent,kKromTrainingDomainId,
        L"Expand KROM training domain…",fonts.body);add(trainingDomainButton_);
    kromSearchButton_=makeButton(parent,kKromControlSearchId,
        L"KROM control search…",fonts.body);add(kromSearchButton_);
    decompositionButton_=makeButton(parent,kDecompositionSettingsId,L"Configure decompositions…",fonts.body);add(decompositionButton_);
    kromValidationButton_=makeButton(parent,kKromValidationSettingsId,L"Configure KROM validation…",fonts.body);add(kromValidationButton_);
    boundaryEigenfunctionButton_=makeButton(parent,kBoundaryEigenfunctionSettingsId,
        L"Configure boundary functions…",fonts.body);add(boundaryEigenfunctionButton_);
    eigenvalueOptimizationButton_=makeButton(parent,kEigenvalueOptimizationSettingsId,
        L"Configure eigenvalue optimization…",fonts.body);add(eigenvalueOptimizationButton_);
    archive_ = makeCheck(parent, 0, L"Publish derived results to the canonical archive", fonts.body, true); add(archive_);
    add(makeButton(parent, kPostPlanId, L"Plan / validate", fonts.body));
    add(makeButton(parent, kPostRunId, L"Run selected stages", fonts.body));
    return true;
}

void PostprocessPanel::show(bool visible) { for (HWND item : controls_) ShowWindow(item, visible ? SW_SHOW : SW_HIDE); }

void PostprocessPanel::layout(int x, int y, int width, int height) {
    const int gap=scalePixels(parent_,12),inner=scalePixels(parent_,7);
    const int column=(width-2*gap)/3,labelH=scalePixels(parent_,15);
    const int fieldH=scalePixels(parent_,28),buttonH=scalePixels(parent_,30);
    auto field = [&](int label, HWND window, int px, int py, int pw, int ph) {
        MoveWindow(labels_[static_cast<std::size_t>(label)],px,py,pw,labelH,TRUE);
        MoveWindow(window,px,py+labelH,pw,ph,TRUE);
    };
    const int middle=x+column+gap,right=middle+column+gap;
    field(0,label_,x,y,column,fieldH);
    field(1,source_,middle,y,2*column+gap,fieldH);

    const int topHeight=scalePixels(parent_,58);
    const int actionH=scalePixels(parent_,36),actionY=y+height-actionH;
    const int contentY=y+topHeight;
    const int contentHeight=std::max(scalePixels(parent_,500),
        actionY-gap-contentY);
    MoveWindow(labels_[2],x,contentY,column,labelH,TRUE);
    HWND methods[] = {hdmd_,sindy_,havok_,dmd_,dmdControl_,dmdControlButton_,edmdControl_,edmdControlButton_,edmd_,pod_,pca_,svd_,baggedPca_,
        robustPca_,hankelPod_,controlRegression_,kromExtract_,kromValidation_,
        boundaryEigenfunctions_,eigenvalueOptimization_,observability_,klqe_};
    const int methodsY=contentY+labelH+inner;
    const int archiveH=scalePixels(parent_,27);
    const int configurationRows=4*buttonH+3*inner;
    const int archiveY=actionY-gap-archiveH;
    int configureY=archiveY-inner-configurationRows;
    const int methodArea=configureY-inner-methodsY;
    const int methodCount=static_cast<int>(sizeof(methods)/sizeof(methods[0]));
    const int methodPitch=std::clamp(methodArea/((methodCount+1)/2),
        scalePixels(parent_,24),scalePixels(parent_,34));
    const int methodWidth=(column-inner)/2;
    for(int index=0;index<methodCount;++index)MoveWindow(methods[index],
        x+(index%2)*(methodWidth+inner),methodsY+(index/2)*methodPitch,
        methodWidth,scalePixels(parent_,24),TRUE);
    MoveWindow(pipelineButton_,x,configureY,column,buttonH,TRUE);
    configureY+=buttonH+inner;
    MoveWindow(trainingDomainButton_,x,configureY,methodWidth,buttonH,TRUE);
    MoveWindow(kromSearchButton_,x+methodWidth+inner,configureY,
        column-methodWidth-inner,buttonH,TRUE);
    configureY+=buttonH+inner;
    MoveWindow(eigenvalueOptimizationButton_,x,configureY,methodWidth,buttonH,TRUE);
    MoveWindow(boundaryEigenfunctionButton_,x+methodWidth+inner,configureY,
        column-methodWidth-inner,buttonH,TRUE);
    configureY+=buttonH+inner;
    MoveWindow(kromValidationButton_,x,configureY,methodWidth,buttonH,TRUE);
    MoveWindow(decompositionButton_,x+methodWidth+inner,configureY,
        column-methodWidth-inner,buttonH,TRUE);
    MoveWindow(archive_,x,archiveY,column,archiveH,TRUE);
    MoveWindow(controls_[controls_.size()-2],x,actionY,methodWidth,actionH,TRUE);
    MoveWindow(controls_[controls_.size()-1],x+methodWidth+inner,actionY,
        column-methodWidth-inner,actionH,TRUE);

    HWND settings[]={hdmdSettings_,sindySettings_,havokSettings_,edmdSettings_,pcaSettings_};
    const int settingsBlock=contentHeight/5;
    for(int index=0;index<5;++index){
        const int blockY=contentY+index*settingsBlock;
        const int editorH=std::max(fieldH,settingsBlock-labelH-inner);
        field(3+index,settings[index],middle,blockY,column,editorH);
    }

    HWND advanced[]={controlRegressionSettings_,kromExtractSettings_,estimationSettings_};
    HWND advancedButtons[]={controlRegressionButton_,kromModesButton_,estimationButton_};
    const int advancedBlock=contentHeight/3;
    for(int index=0;index<3;++index){
        const int blockY=contentY+index*advancedBlock;
        const int editorH=std::max(fieldH,
            advancedBlock-labelH-buttonH-2*inner);
        field(8+index,advanced[index],right,blockY,column,editorH);
        MoveWindow(advancedButtons[index],right,
            blockY+labelH+editorH+inner,column,buttonH,TRUE);
    }
}

void PostprocessPanel::setSource(const std::filesystem::path& source) {
    SetWindowTextW(source_, source.wstring().c_str());
}

std::filesystem::path PostprocessPanel::source() const {
    return std::filesystem::u8path(utf8Text(source_));
}

void PostprocessPanel::command(int id, int notification) {
    if(id==kDmdControlSettingsId&&notification==BN_CLICKED){editDmdControl();return;}
    if(id==kEdmdControlSettingsId&&notification==BN_CLICKED){editDmdControl(true);return;}
    if(notification==BN_CLICKED&&((id==kDmdControlStageId&&checked(dmdControl_))
        ||(id==kEdmdControlStageId&&checked(edmdControl_)))){
        SendMessageW(hdmd_,BM_SETCHECK,BST_UNCHECKED,0);
        SendMessageW(sindy_,BM_SETCHECK,BST_UNCHECKED,0);return;}
    if(id==kKoopmanPipelineId&&notification==BN_CLICKED){editKoopmanPipeline();return;}
    if(id==kKromModesId&&notification==BN_CLICKED){editKromModes();return;}
    if(id==kControlRegressionSettingsId&&notification==BN_CLICKED){editControlRegression();return;}
    if(id==kEstimationSettingsId&&notification==BN_CLICKED){editEstimation();return;}
    if(id==kDecompositionSettingsId&&notification==BN_CLICKED){editDecompositions();return;}
    if(id==kKromValidationSettingsId&&notification==BN_CLICKED){editKromValidation();return;}
    if(id==kBoundaryEigenfunctionSettingsId&&notification==BN_CLICKED){
        editBoundaryEigenfunctions();return;}
    if(id==kEigenvalueOptimizationSettingsId&&notification==BN_CLICKED){
        editEigenvalueOptimization();return;}
    if(id==kKromExtractId&&notification==BN_CLICKED&&checked(kromExtract_)){
        // KROM extraction consumes a KROM archive rather than a trajectory.
        SendMessageW(hdmd_,BM_SETCHECK,BST_UNCHECKED,0);
        SendMessageW(sindy_,BM_SETCHECK,BST_UNCHECKED,0);
    }
}

void PostprocessPanel::editEigenvalueOptimization(){
    const auto result=showEigenvalueOptimizationDialog(parent_,
        std::filesystem::u8path(utf8Text(source_)),*registry_,
        eigenvalueOptimizationSettingsJson_,*fonts_);
    if(!result)return;eigenvalueOptimizationSettingsJson_=*result;
    SendMessageW(eigenvalueOptimization_,BM_SETCHECK,BST_CHECKED,0);
}

void PostprocessPanel::editBoundaryEigenfunctions(){
    const auto result=showBoundaryEigenfunctionDialog(parent_,
        std::filesystem::u8path(utf8Text(source_)),*registry_,
        boundaryEigenfunctionSettingsJson_,*fonts_);
    if(!result)return;boundaryEigenfunctionSettingsJson_=*result;
    SendMessageW(boundaryEigenfunctions_,BM_SETCHECK,BST_CHECKED,0);
    SendMessageW(hdmd_,BM_SETCHECK,BST_UNCHECKED,0);
    SendMessageW(sindy_,BM_SETCHECK,BST_UNCHECKED,0);
}

void PostprocessPanel::editKromValidation() {
    const auto result=showKromValidationDialog(parent_,
        std::filesystem::u8path(utf8Text(source_)),*registry_,
        kromValidationSettingsJson_,*fonts_);
    if(!result)return;kromValidationSettingsJson_=*result;
    SendMessageW(kromValidation_,BM_SETCHECK,BST_CHECKED,0);
    SendMessageW(hdmd_,BM_SETCHECK,BST_UNCHECKED,0);
    SendMessageW(sindy_,BM_SETCHECK,BST_UNCHECKED,0);
}

void PostprocessPanel::editDecompositions() {
    const auto result=showDecompositionDialog(parent_,std::filesystem::u8path(utf8Text(source_)),
        {utf8Text(dmdSettings_),utf8Text(havokSettings_),utf8Text(edmdSettings_),svdSettingsJson_,
         utf8Text(pcaSettings_),podSettingsJson_,baggedPcaSettingsJson_,robustPcaSettingsJson_,hankelPodSettingsJson_},*fonts_);
    if(!result)return;SetWindowTextW(dmdSettings_,dsl::wide(result->dmd).c_str());
    SetWindowTextW(havokSettings_,dsl::wide(result->havok).c_str());
    SetWindowTextW(edmdSettings_,dsl::wide(result->edmd).c_str());
    SetWindowTextW(pcaSettings_,dsl::wide(result->pca).c_str());
    svdSettingsJson_=result->svd;podSettingsJson_=result->pod;
    baggedPcaSettingsJson_=result->baggedPca;robustPcaSettingsJson_=result->robustPca;
    hankelPodSettingsJson_=result->hankelPod;
}

void PostprocessPanel::editDmdControl(bool extended) {
    auto& settings=extended?edmdControlSettingsJson_:dmdControlSettingsJson_;
    const auto result=showDmdControlDialog(parent_,settings,*fonts_,extended);
    if(!result)return;
    settings=*result;
    SendMessageW(extended?edmdControl_:dmdControl_,BM_SETCHECK,BST_CHECKED,0);
    command(extended?kEdmdControlStageId:kDmdControlStageId,BN_CLICKED);
}

void PostprocessPanel::editEstimation() {
    const auto result=showEstimationDialog(parent_,std::filesystem::u8path(utf8Text(source_)),
        *registry_,utf8Text(estimationSettings_),*fonts_);
    if(!result)return;SetWindowTextW(estimationSettings_,dsl::wide(*result).c_str());
}

void PostprocessPanel::editControlRegression() {
    const auto result=showControlRegressionDialog(parent_,std::filesystem::u8path(utf8Text(source_)),
        *registry_,utf8Text(controlRegressionSettings_),*fonts_);
    if(!result)return;SetWindowTextW(controlRegressionSettings_,dsl::wide(*result).c_str());
    SendMessageW(controlRegression_,BM_SETCHECK,BST_CHECKED,0);
}

void PostprocessPanel::editKromModes() {
    const auto result=showKromExtractionDialog(parent_,std::filesystem::u8path(utf8Text(source_)),
        utf8Text(kromExtractSettings_),*fonts_);
    if(!result)return;SetWindowTextW(kromExtractSettings_,dsl::wide(*result).c_str());
    SendMessageW(kromExtract_,BM_SETCHECK,BST_CHECKED,0);
    SendMessageW(hdmd_,BM_SETCHECK,BST_UNCHECKED,0);SendMessageW(sindy_,BM_SETCHECK,BST_UNCHECKED,0);
}

void PostprocessPanel::editKoopmanPipeline() {
    const auto result=showKoopmanPipelineDialog(parent_,std::filesystem::u8path(utf8Text(source_)),
        utf8Text(hdmdSettings_),utf8Text(sindySettings_),*fonts_);
    if(!result)return;
    SetWindowTextW(hdmdSettings_,dsl::wide(result->hdmd).c_str());
    SetWindowTextW(sindySettings_,dsl::wide(result->sindy).c_str());
    SendMessageW(hdmd_,BM_SETCHECK,BST_CHECKED,0);
    SendMessageW(sindy_,BM_SETCHECK,BST_CHECKED,0);
}

std::string PostprocessPanel::jobJson() const {
    std::vector<std::pair<std::string,std::string>> methods;
    if (checked(hdmd_)) methods.push_back({"hdmd",jsonObject(utf8Text(hdmdSettings_))});
    if (checked(eigenvalueOptimization_)) methods.push_back(
        {"eigenvalue_optimization",jsonObject(eigenvalueOptimizationSettingsJson_)});
    if (checked(sindy_)) methods.push_back({"sindy",jsonObject(utf8Text(sindySettings_))});
    if (checked(havok_)) methods.push_back({"havok",jsonObject(utf8Text(havokSettings_))});
    if (checked(dmd_)) methods.push_back({"dmd",jsonObject(utf8Text(dmdSettings_))});
    if (checked(dmdControl_)) methods.push_back({"dmdc",jsonObject(dmdControlSettingsJson_)});
    if (checked(edmdControl_)) methods.push_back({"edmdc",jsonObject(edmdControlSettingsJson_)});
    if (checked(edmd_)) methods.push_back({"edmd",jsonObject(utf8Text(edmdSettings_))});
    if (checked(pod_)) methods.push_back({"pod",jsonObject(podSettingsJson_)});
    if (checked(pca_)) methods.push_back({"pca",jsonObject(utf8Text(pcaSettings_))});
    if (checked(svd_)) methods.push_back({"svd",jsonObject(svdSettingsJson_)});
    if (checked(baggedPca_)) methods.push_back({"bagged_pca",jsonObject(baggedPcaSettingsJson_)});
    if (checked(robustPca_)) methods.push_back({"robust_pca",jsonObject(robustPcaSettingsJson_)});
    if (checked(hankelPod_)) methods.push_back({"hankel_pod",jsonObject(hankelPodSettingsJson_)});
    if (checked(controlRegression_)) methods.push_back({"control_regression",jsonObject(utf8Text(controlRegressionSettings_))});
    if (checked(kromExtract_)) methods.push_back({"krom_extract",jsonObject(utf8Text(kromExtractSettings_))});
    if (checked(kromValidation_)) methods.push_back({"koopman_validation",jsonObject(kromValidationSettingsJson_)});
    if(checked(boundaryEigenfunctions_))methods.push_back(
        {"boundary_eigenfunctions",jsonObject(boundaryEigenfunctionSettingsJson_)});
    if (checked(observability_)) methods.push_back({"observability",jsonObject(utf8Text(estimationSettings_))});
    if (checked(klqe_)) methods.push_back({"klqe",jsonObject(utf8Text(estimationSettings_))});
    std::ostringstream out;
    out << "{\"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"postprocess\",\"label\":\""
        << jsonEscape(utf8Text(label_)) << "\",\"sources\":[\"" << jsonEscape(utf8Text(source_))
        << "\"],\"methods\":[";
    for (std::size_t index=0; index<methods.size(); ++index) {
        if (index) out << ',';
        out << "{\"name\":\"" << methods[index].first << "\",\"enabled\":true,\"settings\":" << methods[index].second << '}';
    }
    out << "],\"archive\":{\"enabled\":" << (checked(archive_) ? "true" : "false") << "}}";
    return out.str();
}

}  // namespace dsl::compute
