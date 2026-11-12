#include "dynsyslab_compute/Workbench.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "ArchiveToolsDialog.hpp"
#include "QueueDialog.hpp"
#include "LegacyImportDialog.hpp"
#include "KromControlSearchDialog.hpp"
#include "dynsyslab_viewer/Help.hpp"
#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"
#include "dynsyslab_viewer/StartupDiagnostics.hpp"
#include "dynsyslab_resource.h"

#include <commdlg.h>
#include <commctrl.h>
#include <shellapi.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace dsl::compute {
namespace {
constexpr wchar_t kClassName[] = L"DynSysLabComputeWorkbench";
constexpr int kNavigationId=5001, kSearchId=5002, kRefreshId=5003;
constexpr int kSolverTabId=5004, kPostTabId=5005, kViewerId=5006;
constexpr int kHelpId=5007;
constexpr int kSaveJobId=5008, kOpenJobId=5009;
constexpr int kQueueId=5010;
constexpr int kArchiveToolsId=5011;
constexpr int kImportDataId=5012;
constexpr int kCancelTaskId=5014;
constexpr int kLeftSplitterId=5015, kRightSplitterId=5016;
constexpr int kStageSystemId=5017, kStageInitialConditionsId=5018;
constexpr int kStageControlId=5019, kStageIntegratorId=5020, kStageReviewId=5021;
constexpr UINT_PTR kTimer=5050;
HWND button(HWND parent,int id,const wchar_t* value,HFONT font) { return makeButton(parent,id,value,font); }
std::filesystem::path sourcePath(const BackendConfig& config,const dsl::CatalogItem& item) {
    std::filesystem::path source = dsl::wide(item.source);
    return source.is_absolute() ? source : config.projectRoot / source;
}

std::filesystem::path chooseJobFile(HWND owner,bool save) {
    wchar_t path[32768]{};
    if(save) wcscpy_s(path,L"dynsyslab-job.json");
    const wchar_t filter[]=L"DynSysLab jobs (*.json)\0*.json\0All files (*.*)\0*.*\0\0";
    OPENFILENAMEW request{};request.lStructSize=sizeof(request);request.hwndOwner=owner;
    request.lpstrFile=path;request.nMaxFile=32768;request.lpstrFilter=filter;
    request.nFilterIndex=1;request.lpstrDefExt=L"json";request.lpstrTitle=save
        ?L"Save DynSysLab job":L"Open DynSysLab job";
    request.Flags=OFN_EXPLORER|OFN_NOCHANGEDIR|OFN_PATHMUSTEXIST
        |(save?OFN_OVERWRITEPROMPT:OFN_FILEMUSTEXIST);
    const BOOL accepted=save?GetSaveFileNameW(&request):GetOpenFileNameW(&request);
    return accepted?std::filesystem::path(path):std::filesystem::path{};
}

std::string readJobFile(const std::filesystem::path& path){std::ifstream stream(path,std::ios::binary);
    if(!stream)throw std::runtime_error("Unable to open the selected job file");
    std::ostringstream text;text<<stream.rdbuf();return text.str();}

void writeJobFile(const std::filesystem::path& path,const std::string& text){std::ofstream stream(path,std::ios::binary|std::ios::trunc);
    if(!stream)throw std::runtime_error("Unable to create the selected job file");
    stream.write(text.data(),static_cast<std::streamsize>(text.size()));
    if(!stream)throw std::runtime_error("Unable to finish writing the selected job file");}

bool pathInside(const std::filesystem::path& candidate,
                const std::filesystem::path& directory) {
    const auto path = std::filesystem::absolute(candidate).lexically_normal();
    const auto root = std::filesystem::absolute(directory).lexically_normal();
    const auto relative = path.lexically_relative(root);
    return !relative.empty() && *relative.begin() != L"..";
}

}

Workbench::Workbench(Arguments arguments) : arguments_(std::move(arguments)),
    backend_(arguments_.backend),
    registry_(backend_.loadRegistry(catalog_)) {}

bool Workbench::registerClass(HINSTANCE instance) {
    WNDCLASSEXW type{sizeof(type)}; type.lpfnWndProc=procedure; type.hInstance=instance;
    type.hCursor=LoadCursorW(nullptr,IDC_ARROW);
    type.hIcon=LoadIconW(instance,MAKEINTRESOURCEW(IDI_DYNSYSLAB_APP));
    type.hIconSm=static_cast<HICON>(LoadImageW(instance,MAKEINTRESOURCEW(IDI_DYNSYSLAB_APP),
        IMAGE_ICON,GetSystemMetrics(SM_CXSMICON),GetSystemMetrics(SM_CYSMICON),LR_DEFAULTCOLOR));
    type.hbrBackground=CreateSolidBrush(dsl::theme::background); type.lpszClassName=kClassName;
    return (RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS) &&
        dsl::InspectorPane::registerClass(instance) && dsl::PlotCanvas::registerClass(instance);
}

int Workbench::run(HINSTANCE instance,int) {
    window_=CreateWindowExW(WS_EX_COMPOSITED,kClassName,L"DynSysLab Solver and Post-processor",
        WS_OVERLAPPEDWINDOW|WS_CLIPCHILDREN,CW_USEDEFAULT,CW_USEDEFAULT,1540,930,
        nullptr,nullptr,instance,this);
    if(!window_) throw std::runtime_error("Unable to create compute workbench window");
    dsl::theme::configureTitleBar(window_);
    if(!arguments_.backgroundCheck){ShowWindow(window_,SW_MAXIMIZE);UpdateWindow(window_);}
    else {RECT workspace{};if(SystemParametersInfoW(SPI_GETWORKAREA,0,&workspace,0))
        SetWindowPos(window_,nullptr,workspace.left,workspace.top,workspace.right-workspace.left,
            workspace.bottom-workspace.top,SWP_NOACTIVATE|SWP_NOZORDER);layout();SetTimer(window_,5051,30000,nullptr);}
    MSG event{}; while(GetMessageW(&event,nullptr,0,0)>0){TranslateMessage(&event);DispatchMessageW(&event);}
    return static_cast<int>(event.wParam);
}

bool Workbench::createChildren() {
    search_=makeEdit(window_,kSearchId,L"",fonts_.body); SendMessageW(search_,EM_SETCUEBANNER,TRUE,reinterpret_cast<LPARAM>(L"Filter shared archives…"));
    refresh_=button(window_,kRefreshId,L"Refresh shared archives",fonts_.heading);
    solverTab_=button(window_,kSolverTabId,L"Solver",fonts_.heading);
    postTab_=button(window_,kPostTabId,L"Post-processor",fonts_.heading);
    viewer_=button(window_,kViewerId,L"Open Viewer",fonts_.body);
    stageSystem_=button(window_,kStageSystemId,L"1  System",fonts_.body);
    stageInitialConditions_=button(window_,kStageInitialConditionsId,L"2  ICs",fonts_.body);
    stageControl_=button(window_,kStageControlId,L"3  Control",fonts_.body);
    stageIntegrator_=button(window_,kStageIntegratorId,L"4  Integrator",fonts_.body);
    stageReview_=button(window_,kStageReviewId,L"5  Review",fonts_.body);
    saveJob_=button(window_,kSaveJobId,L"Save job…",fonts_.body);
    openJob_=button(window_,kOpenJobId,L"Open job…",fonts_.body);
    queue_=button(window_,kQueueId,L"Native jobs…",fonts_.body);
    archiveTools_=button(window_,kArchiveToolsId,L"Archives…",fonts_.body);
    importData_=button(window_,kImportDataId,L"Import data…",fonts_.body);
    help_=button(window_,kHelpId,L"Help / About",fonts_.body);
    headerTooltip_=CreateWindowExW(WS_EX_TOPMOST,TOOLTIPS_CLASSW,nullptr,
        WS_POPUP|TTS_ALWAYSTIP|TTS_NOPREFIX,CW_USEDEFAULT,CW_USEDEFAULT,
        CW_USEDEFAULT,CW_USEDEFAULT,window_,nullptr,GetModuleHandleW(nullptr),nullptr);
    auto addTip=[&](HWND control,const wchar_t* text){TTTOOLINFOW info{};
        info.cbSize=TTTOOLINFOW_V1_SIZE;info.uFlags=TTF_IDISHWND|TTF_SUBCLASS;
        info.hwnd=window_;info.uId=reinterpret_cast<UINT_PTR>(control);
        info.lpszText=const_cast<wchar_t*>(text);
        SendMessageW(headerTooltip_,TTM_ADDTOOLW,0,reinterpret_cast<LPARAM>(&info));};
    addTip(refresh_,L"Refresh the shared archive index");
    addTip(help_,L"Help and application information");
    addTip(stageSystem_,L"Edit the dynamical system and its parameters");
    addTip(stageInitialConditions_,L"Edit initial conditions");
    addTip(stageControl_,L"Configure forcing and control");
    addTip(stageIntegrator_,L"Choose integration and output settings");
    addTip(stageReview_,L"Validate, run, and review provenance");
    output_=makeEdit(window_,0,L"Plan and execution results appear here.",fonts_.mono,true);
    SendMessageW(output_,EM_SETREADONLY,TRUE,0);
    status_=makeLabel(window_,L"Ready",fonts_.caption);
    progressLabel_=makeLabel(window_,L"Idle",fonts_.caption);
    progress_=makeProgress(window_,fonts_.body);
    jobActionsCaption_=makeLabel(window_,L"Review and run",fonts_.caption);
    reviewRailTitle_=makeLabel(window_,L"Review, output, and provenance",fonts_.heading);
    cancelTask_=button(window_,kCancelTaskId,L"Cancel active",fonts_.body);
    EnableWindow(cancelTask_,FALSE);
    if(!headerTooltip_||!navigation_.create(window_,kNavigationId,fonts_)||
       !solveJobs_.create(window_,fonts_,kSolverCampaignId)||
       !solveContext_.create(window_,fonts_)||
       !inspector_.create(window_,fonts_)||!solver_.create(window_,fonts_,&registry_)||
       !leftSplitter_.create(window_,kLeftSplitterId)||
       !rightSplitter_.create(window_,kRightSplitterId)) return false;
    solveJobs_.reset(solver_.systemLabel(),solver_.jobJson("New Job"));
    activeSolveSystem_=solveJobs_.firstSystem();
    SetWindowTextW(status_,L"Loading the central archive index…");
    setTab(!arguments_.openPostprocessor);
    SetTimer(window_,kTimer,80,nullptr);
    taskActive_=true;
    task_=std::async(std::launch::async,[this](){TaskResult result;result.startup=true;
        try{result.catalog=std::make_unique<dsl::Catalog>(backend_.loadCatalog());}
        catch(const std::exception& error){result.error=dsl::wide(error.what());}
        return result;});
    return true;
}

void Workbench::setTab(bool solver) {
    if (!solver && !ensurePostPanel()) return;
    if (solverActive_ && !solver) storeActiveSolve();
    solverActive_=solver; solver_.showEditor(false);
    solver_.showExecution(false, false); solver_.showJobActions(false);
    solveContext_.show(false); if(postCreated_)post_.show(false);
    ShowWindow(jobActionsCaption_,SW_HIDE);
    solveJobs_.show(false);ShowWindow(navigation_.handle(),SW_HIDE);
    ShowWindow(search_,SW_HIDE);ShowWindow(refresh_,SW_HIDE);
    // Clear the shared panel canvas while neither page clips the parent
    // repaint, then reveal and redraw only the selected page.
    RedrawWindow(window_,nullptr,nullptr,RDW_INVALIDATE|RDW_UPDATENOW);
    if(postCreated_)post_.show(!solver);solveJobs_.show(solver);
    solver_.showJobActions(solver);ShowWindow(jobActionsCaption_,solver?SW_SHOW:SW_HIDE);
    solver_.showExecution(solver, solver);
    if(solver)selectSolveNode();
    if(!solver){ensureArchiveNavigation();ShowWindow(navigation_.handle(),SW_SHOW);ShowWindow(search_,SW_SHOW);
        ShowWindow(refresh_,SW_SHOW);}
    EnableWindow(solverTab_,!solver); EnableWindow(postTab_,solver); layout();
    RedrawWindow(window_,nullptr,nullptr,RDW_INVALIDATE|RDW_ALLCHILDREN|RDW_UPDATENOW);
}

void Workbench::selectArchive(int index) {
    if(index<0||static_cast<std::size_t>(index)>=catalog_.items.size()) return;
    selected_=index; const auto& item=catalog_.items[static_cast<std::size_t>(index)];
    inspector_.setItem(&catalog_,&item);
    const auto path=sourcePath(arguments_.backend,item);
    if(path.extension()==L".npz") post_.setSource(path);
    SetWindowTextW(status_,(L"Selected " + dsl::wide(item.label)).c_str());
}

void Workbench::refreshCatalog() {
    if(taskActive_) return; taskActive_=true;
    SetWindowTextW(status_,L"Reconciling the archive index with files on disk…");
    beginProgress(false);
    task_=std::async(std::launch::async,[this](){TaskResult result;result.refresh=true;
        try{backend_.refreshCatalog();result.text=L"Catalog refreshed";}catch(const std::exception& e){result.error=dsl::wide(e.what());}return result;});
}

void Workbench::runJob(bool solve,bool execute) {
    runJobJson(solve?activeJobJson():post_.jobJson(),execute);
}

bool Workbench::ensurePostPanel() {
    if (postCreated_) return true;
    postCreated_ = post_.create(window_, fonts_, &registry_);
    if (!postCreated_)
        SetWindowTextW(status_, L"Unable to create the Post-processor panel");
    return postCreated_;
}

void Workbench::ensureArchiveNavigation() {
    if (archiveNavigationLoaded_) return;
    navigation_.setCatalog(&catalog_);
    navigation_.resetHorizontalScroll();
    archiveNavigationLoaded_ = true;
}

void Workbench::selectInitialArchive() {
    if(arguments_.initialSource.empty()) return;
    std::error_code ignored;
    const auto requested=std::filesystem::weakly_canonical(arguments_.initialSource,ignored);
    for(std::size_t index=0;index<catalog_.items.size();++index) {
        const auto candidate=std::filesystem::weakly_canonical(
            sourcePath(arguments_.backend,catalog_.items[index]),ignored);
        if(!ignored&&candidate==requested) {
            ensureArchiveNavigation();navigation_.selectIndex(static_cast<int>(index));
            selectArchive(static_cast<int>(index));break;
        }
        ignored.clear();
    }
}

void Workbench::runJobJson(std::string job,bool execute) {
    if(taskActive_) return; taskActive_=true;
    const std::string action=execute?"run":"plan";
    SetWindowTextW(output_,execute?L"Executing through the canonical DynSysLab services…":L"Validating and planning without numerical execution…");
    SetWindowTextW(status_,execute?L"Running job…":L"Planning job…");
    beginProgress(execute);
    const auto progress=progressPath_;
    task_=std::async(std::launch::async,[this,job,action,execute,progress](){TaskResult result;result.refresh=execute;
        try{result.text=backend_.execute(action,job,progress,[this]{return cancelRequested_.load();});}
        catch(const dynsyslab::core::ExecutionCancelled&){result.cancelled=true;result.text=L"Native execution cancelled.";}
        catch(const std::exception& e){result.error=dsl::wide(e.what());}return result;});
}

void Workbench::beginProgress(bool external) {
    cancelRequested_.store(false);EnableWindow(cancelTask_,external);
    progressDeterminate_=false; SendMessageW(progress_,PBM_SETPOS,0,0);
    SendMessageW(progress_,PBM_SETMARQUEE,TRUE,35); SetWindowTextW(progressLabel_,L"Working…");
    progressPath_.clear();
    if(external){progressPath_=arguments_.backend.cache/L"active-progress.json";
        std::error_code ignored;std::filesystem::remove(progressPath_,ignored);}
    lastProgressPoll_=std::chrono::steady_clock::time_point{};
}

void Workbench::pollProgress() {
    if(progressPath_.empty()) return;
    const auto now=std::chrono::steady_clock::now();
    if(lastProgressPoll_.time_since_epoch().count()&&now-lastProgressPoll_<std::chrono::milliseconds(250)) return;
    lastProgressPoll_=now; const auto value=backend_.progress(progressPath_); if(!value.available) return;
    if(value.determinate&&!progressDeterminate_){SendMessageW(progress_,PBM_SETMARQUEE,FALSE,0);progressDeterminate_=true;}
    if(value.determinate) SendMessageW(progress_,PBM_SETPOS,
        static_cast<WPARAM>(std::clamp(value.percent*10.0,0.0,1000.0)),0);
    std::wostringstream label; if(value.determinate) label<<static_cast<int>(value.percent+0.5)<<L"%";
    else label<<L"Running";
    if(value.estimatedRemainingSeconds>=0) label<<L" · "<<static_cast<int>(value.estimatedRemainingSeconds+0.5)<<L"s";
    SetWindowTextW(progressLabel_,label.str().c_str());
    if(!value.message.empty()) SetWindowTextW(status_,value.message.c_str());
}

void Workbench::endProgress(bool success) {
    EnableWindow(cancelTask_,FALSE);
    SendMessageW(progress_,PBM_SETMARQUEE,FALSE,0);progressDeterminate_=true;
    SendMessageW(progress_,PBM_SETPOS,success?1000:0,0);
    SetWindowTextW(progressLabel_,success?L"Complete":L"Failed");
    if(!progressPath_.empty()){std::error_code ignored;std::filesystem::remove(progressPath_,ignored);progressPath_.clear();}
}

void Workbench::pollTask() {
    if(!taskActive_) return; pollProgress();
    if(task_.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready) return;
    TaskResult result=task_.get();taskActive_=false;
    if(result.startup){
        if(!result.error.empty()){SetWindowTextW(status_,L"Archive index could not be loaded");
            SetWindowTextW(output_,result.error.c_str());if(arguments_.backgroundCheck)finishBackgroundCheck(1,result.error);return;}
        try {
            if(!result.catalog)throw std::runtime_error("archive catalogue was not returned");
            catalog_=std::move(*result.catalog);registry_=backend_.loadRegistry(catalog_);
            solver_.refreshRegistry();archiveNavigationLoaded_=false;
            if(!solverActive_)ensureArchiveNavigation();selectInitialArchive();
            dsl::writeStartupDiagnostics("Solver and Post-processor", arguments_.backend.cache,
                arguments_.backend.projectRoot, arguments_.backend.archiveRoot,
                arguments_.backend.catalog, catalog_.items.size(), navigation_.nodeCount());
            SetWindowTextW(status_,(L"Ready  •  "+std::to_wstring(catalog_.items.size())+
                L" archive objects").c_str());
        } catch(const std::exception& error) {
            SetWindowTextW(status_,L"Archive catalogue could not be installed");
            SetWindowTextW(output_,dsl::wide(error.what()).c_str());
        }
        if(arguments_.backgroundCheck)finishBackgroundCheck(0);
        return;}
    if(!result.error.empty()){endProgress(false);SetWindowTextW(output_,result.error.c_str());SetWindowTextW(status_,L"Job failed validation or execution");return;}
    if(result.cancelled){endProgress(false);SetWindowTextW(progressLabel_,L"Cancelled");
        SetWindowTextW(output_,result.text.c_str());SetWindowTextW(status_,L"Native job cancelled");
        if(result.refresh){try{catalog_=backend_.loadCatalog();archiveNavigationLoaded_=false;
            if(!solverActive_)ensureArchiveNavigation();selected_=-1;}catch(...) {}}return;}
    endProgress(true);
    SetWindowTextW(output_,result.text.c_str());
    SetWindowTextW(status_,result.refresh?L"Archive catalog refreshed":L"Job complete");
    if(result.refresh){try{catalog_=backend_.loadCatalog();archiveNavigationLoaded_=false;
        if(!solverActive_)ensureArchiveNavigation();selected_=-1;}catch(...) {}}
}

void Workbench::openViewer() {
    const bool opened=backend_.openViewer();
    SetWindowTextW(status_,opened?L"Viewer opened on the shared archive tree":
        L"Viewer could not be opened; see the application logs or launch it from the CLI");
}

void Workbench::saveJobDefinition(){
    if(taskActive_)return;
    try{const std::string job=activeJobJson();
        (void)backend_.execute("plan",job);const auto path=chooseJobFile(window_,true);if(path.empty())return;
        if(pathInside(path,arguments_.backend.archiveRoot))
            throw std::invalid_argument("Job definitions are configuration, not scientific data; save this file outside the archive root");
        writeJobFile(path,dsl::json::dump(dsl::json::parse(job)));
        SetWindowTextW(status_,(L"Saved native job: "+path.filename().wstring()).c_str());
    }catch(const std::exception& error){SetWindowTextW(output_,dsl::wide(error.what()).c_str());SetWindowTextW(status_,L"Job was not saved");}
}

void Workbench::openJobDefinition(){
    if(taskActive_)return;const auto path=chooseJobFile(window_,false);if(path.empty())return;
    try{const std::string job=readJobFile(path);const auto parsed=dsl::json::parse(job);
        const std::string kind=parsed.at("kind").string();if(kind!="solve"&&kind!="solve_set"&&kind!="postprocess"&&kind!="krom_control_search")throw std::invalid_argument("Selected file is not a solve, solve-set, postprocess, or KROM search job");
        const int answer=MessageBoxW(window_,L"Run this saved job now?\n\nChoose No to validate and plan it without numerical execution.",L"Open DynSysLab job",MB_YESNOCANCEL|MB_ICONQUESTION);
        if(answer==IDCANCEL)return;
        if(kind=="solve"||kind=="solve_set"){solveJobs_.loadJobJson(job);
            activeSolveSystem_=solveJobs_.firstSystem();
            if(activeSolveSystem_){const auto* draft=solveJobs_.system(*activeSolveSystem_);
                if(draft&&!draft->jobJson.empty())solver_.loadJobJson(draft->jobJson,true);}
            setTab(true);runJobJson(activeJobJson(),answer==IDYES);}
        else{setTab(false);runJobJson(job,answer==IDYES);}
    }catch(const std::exception& error){SetWindowTextW(output_,dsl::wide(error.what()).c_str());SetWindowTextW(status_,L"Saved job could not be opened");}
}

void Workbench::openQueue(){
    std::string activeJob;
    try{activeJob=activeJobJson();}
    catch(const std::exception& error){SetWindowTextW(output_,dsl::wide(error.what()).c_str());}
    if(!showQueueDialog(window_,backend_,activeJob,taskActive_,fonts_)||taskActive_)return;
    taskActive_=true;SetWindowTextW(status_,L"Running next native queued job…");
    SetWindowTextW(output_,L"The next queued job is executing through the canonical native workflow.");
    beginProgress(true);const auto progress=progressPath_;
    task_=std::async(std::launch::async,[this,progress](){TaskResult result;
        try{auto queued=backend_.runNextQueued(progress,[this]{return cancelRequested_.load();});
            if(queued.state==dynsyslab::workflow::QueueRunResult::State::empty)
                result.text=L"The native queue is empty.";
            else if(queued.state==dynsyslab::workflow::QueueRunResult::State::cancelled){
                result.cancelled=true;result.text=L"Queued execution was cancelled.";result.refresh=true;
            }else{result.text=L"Queued job completed\r\n\r\n"+queued.run->summary;result.refresh=true;}
        }catch(const std::exception& error){result.error=dsl::wide(error.what());}
        return result;});
}

void Workbench::cancelActiveTask(){if(!taskActive_||!IsWindowEnabled(cancelTask_))return;
    cancelRequested_.store(true);EnableWindow(cancelTask_,FALSE);
    SetWindowTextW(status_,L"Cancellation requested; waiting for a numerical safe point…");
    SetWindowTextW(progressLabel_,L"Cancelling…");}

void Workbench::importLegacyData(){
    if(taskActive_)return;const auto request=showLegacyImportDialog(window_,registry_,fonts_);
    if(!request)return;taskActive_=true;beginProgress(false);
    SetWindowTextW(status_,L"Importing legacy trajectories into canonical archives…");
    SetWindowTextW(output_,L"Reading, validating, grouping, and archiving the tagged CSV without Python or MATLAB…");
    task_=std::async(std::launch::async,[this,request](){TaskResult result;
        try{result.text=backend_.importLegacyCsv(request->source,request->options);result.refresh=true;}
        catch(const std::exception& error){result.error=dsl::wide(error.what());}return result;});
}

void Workbench::openKromControlSearch(){if(taskActive_)return;std::filesystem::path selected;
    if(selected_>=0)selected=sourcePath(arguments_.backend,catalog_.items[static_cast<std::size_t>(selected_)]);
    const auto request=showKromControlSearchDialog(window_,selected,registry_,fonts_);if(!request)return;
    if(request->action==KromControlSearchRequest::Action::save){try{
        (void)backend_.execute("plan",request->jobJson);const auto path=chooseJobFile(window_,true);
        if(path.empty())return;writeJobFile(path,request->jobJson);
        SetWindowTextW(output_,(L"Saved editable native KROM-search job\r\n\r\n"+path.wstring()).c_str());
        SetWindowTextW(status_,L"KROM control-search job saved");}catch(const std::exception& error){
        SetWindowTextW(output_,dsl::wide(error.what()).c_str());SetWindowTextW(status_,L"Campaign job was not saved");}return;}
    if(request->action==KromControlSearchRequest::Action::queue){try{const auto id=backend_.submitQueued(request->jobJson);
        SetWindowTextW(output_,(L"Queued native KROM control search\r\n\r\nJob ID: "+dsl::wide(id)).c_str());
        SetWindowTextW(status_,L"KROM control search queued");}catch(const std::exception& error){
        SetWindowTextW(output_,dsl::wide(error.what()).c_str());SetWindowTextW(status_,L"Campaign could not be queued");}return;}
    runJobJson(request->jobJson,request->action==KromControlSearchRequest::Action::run);}

void Workbench::openArchiveTools(){
    if(taskActive_)return;const auto request=showArchiveToolsDialog(window_,fonts_);
    if(!request)return;taskActive_=true;beginProgress(false);
    const bool apply=request->action==ArchiveToolRequest::Action::deduplicateApply;
    SetWindowTextW(status_,apply?L"Applying recoverable archive deduplication…":
        request->action==ArchiveToolRequest::Action::audit?L"Auditing canonical archives…":
        L"Previewing archive deduplication…");
    SetWindowTextW(output_,L"Inspecting the selected native archive scope…");
    task_=std::async(std::launch::async,[this,request,apply](){TaskResult result;
        try{result.text=request->action==ArchiveToolRequest::Action::audit
                ?backend_.auditArchive(request->scope)
                :backend_.deduplicateArchive(request->scope,apply);
            result.refresh=apply;
        }catch(const std::exception& error){result.error=dsl::wide(error.what());}
        return result;});
}

LRESULT Workbench::message(UINT messageValue,WPARAM wparam,LPARAM lparam) {
    if(messageValue==WM_CREATE) return createChildren()?0:-1;
    if(messageValue==WM_SIZE){layout();return 0;}
    if(messageValue==WM_PAINT){PAINTSTRUCT p{};HDC dc=BeginPaint(window_,&p);paint(dc);EndPaint(window_,&p);return 0;}
    if(messageValue==WM_ERASEBKGND)return 1;
    if(messageValue==WM_DRAWITEM){drawButton(*reinterpret_cast<DRAWITEMSTRUCT*>(lparam));return TRUE;}
    if(messageValue==dsl::kVerticalSplitterMovedMessage){
        RECT area{};GetClientRect(window_,&area);
        const int margin=scalePixels(window_,12),gap=scalePixels(window_,12);
        if(static_cast<int>(wparam)==kLeftSplitterId)
            navigationWidth_=dsl::navigationWidthFromSplitter(static_cast<int>(lparam),margin,gap);
        else if(static_cast<int>(wparam)==kRightSplitterId)
            inspectorWidth_=dsl::inspectorWidthFromSplitter(static_cast<int>(lparam),area.right,margin,gap);
        layout();return 0;}
    if(messageValue==WM_TIMER&&wparam==kTimer){pollTask();return 0;}
    if(messageValue==WM_TIMER&&wparam==5051&&arguments_.backgroundCheck){finishBackgroundCheck(2,L"Solver startup timed out");return 0;}
    if(messageValue==WM_NOTIFY&&reinterpret_cast<NMHDR*>(lparam)->code==TVN_SELCHANGEDW){
        if(reinterpret_cast<NMHDR*>(lparam)->idFrom==kNavigationId&&!solverActive_){selectArchive(navigation_.selectedIndex());return 0;}
        if(reinterpret_cast<NMHDR*>(lparam)->idFrom==kSolveTreeId&&solverActive_){
            if(!solveJobs_.rebuilding())selectSolveNode();return 0;}}
    if(messageValue==WM_COMMAND){const int id=LOWORD(wparam),note=HIWORD(wparam);
        solver_.command(id,note);
        if(postCreated_)post_.command(id,note);
        solveContext_.command(id,note);
        if(id==kClusterIcApplyId){applyClusterInitialConditions();return 0;}
        if(id>=kStageSystemId&&id<=kStageReviewId){
            activeSolverStage_=id==kStageSystemId?SolverStage::system:
                id==kStageInitialConditionsId?SolverStage::initialConditions:
                id==kStageControlId?SolverStage::control:
                id==kStageIntegratorId?SolverStage::integrator:SolverStage::review;
            solver_.activateStage(activeSolverStage_);
            for(HWND stage:{stageSystem_,stageInitialConditions_,stageControl_,stageIntegrator_,stageReview_})
                InvalidateRect(stage,nullptr,FALSE);
            return 0;}
        if(id>=kAddBatchId&&id<=kRemoveSolveNodeId){editSolveTree(id);return 0;}
        if(solverActive_&&id==kSolverSystemId&&note==CBN_SELCHANGE){storeActiveSolve();return 0;}
        if(solverActive_&&(id==kSolverTrajectoryCountId||id==kSolverClusterCountId)
            &&note==EN_KILLFOCUS){storeActiveSolve();return 0;}
        if(id==kSearchId&&note==EN_CHANGE){navigation_.setFilter(windowText(search_));return 0;}
        if(id==kRefreshId){refreshCatalog();return 0;} if(id==kSolverTabId){setTab(true);return 0;}
        if(id==kPostTabId){setTab(false);return 0;} if(id==kViewerId){openViewer();return 0;}
        if(id==kSaveJobId){saveJobDefinition();return 0;}if(id==kOpenJobId){openJobDefinition();return 0;}
        if(id==kQueueId){openQueue();return 0;}
        if(id==kArchiveToolsId){openArchiveTools();return 0;}
        if(id==kImportDataId){importLegacyData();return 0;}
        if(id==kKromControlSearchId){openKromControlSearch();return 0;}
        if(id==kKromTrainingDomainId){openKromTrainingDomainPlanner();return 0;}
        if(id==kCancelTaskId){cancelActiveTask();return 0;}
        if(id==kHelpId){dsl::showHelpAbout(window_,L"DynSysLab Solver and Post-processor",
            arguments_.backend.projectRoot,arguments_.backend.archiveRoot);return 0;}
        if(id==kSolverSystemPageId){activeSolverStage_=SolverStage::system;RedrawWindow(window_,nullptr,nullptr,RDW_INVALIDATE|RDW_ALLCHILDREN);}
        if(id==kSolverControlPageId||id==kSolverControllerId){activeSolverStage_=SolverStage::control;RedrawWindow(window_,nullptr,nullptr,RDW_INVALIDATE|RDW_ALLCHILDREN);}
        if(id==kSolverInitialConditionsId){activeSolverStage_=SolverStage::initialConditions;RedrawWindow(window_,nullptr,nullptr,RDW_INVALIDATE|RDW_ALLCHILDREN);}
        if(id==kSolverPlanId){activeSolverStage_=SolverStage::review;runJob(true,false);return 0;}
        if(id==kSolverRunId){activeSolverStage_=SolverStage::review;runJob(true,true);return 0;}
        if(id==kPostPlanId){runJob(false,false);return 0;} if(id==kPostRunId){runJob(false,true);return 0;}}
    if(messageValue==WM_CTLCOLORSTATIC){static HBRUSH outputBrush=CreateSolidBrush(dsl::theme::panelRaised);static HBRUSH backgroundBrush=CreateSolidBrush(dsl::theme::background);HDC dc=reinterpret_cast<HDC>(wparam);
        if(reinterpret_cast<HWND>(lparam)==output_){SetBkColor(dc,dsl::theme::panelRaised);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(outputBrush);}
        SetBkColor(dc,dsl::theme::background);SetTextColor(dc,dsl::theme::muted);return reinterpret_cast<LRESULT>(backgroundBrush);}
    if(messageValue==WM_CTLCOLOREDIT){static HBRUSH edit=CreateSolidBrush(dsl::theme::panelRaised);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panelRaised);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(edit);}
    if(messageValue==WM_CLOSE&&taskActive_){cancelActiveTask();MessageBoxW(window_,
        L"DynSysLab is cancelling the active computation at a numerical safe point. Close the app after the task reports cancellation.",
        L"Active computation",MB_OK|MB_ICONINFORMATION);return 0;}
    if(messageValue==WM_DESTROY){KillTimer(window_,kTimer);PostQuitMessage(0);return 0;}
    return DefWindowProcW(window_,messageValue,wparam,lparam);
}

LRESULT CALLBACK Workbench::procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* self=reinterpret_cast<Workbench*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){self=static_cast<Workbench*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);self->window_=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(self));}
    return self?self->message(message,wparam,lparam):DefWindowProcW(window,message,wparam,lparam);
}

}  // namespace dsl::compute
