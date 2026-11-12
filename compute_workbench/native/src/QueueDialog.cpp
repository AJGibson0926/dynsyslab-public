#include "QueueDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabQueueDialog";
constexpr int kListId=1701,kRefreshId=1702,kSubmitId=1703,kRunId=1704;
constexpr int kCancelId=1705,kRecoverId=1706,kCloseId=1707,kAgeId=1708;

struct DialogState {
    HWND window{}, list{}, details{}, summary{}, age{}, run{}, cancel{};
    Backend* backend{};
    const dsl::theme::Fonts* fonts{};
    std::string activeJob;
    bool taskActive{}, runNext{};
    std::vector<dynsyslab::core::JobRecord> jobs;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font,bool edge=false) {
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,
        WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;
}

int selection(const DialogState& state) {
    return static_cast<int>(SendMessageW(state.list,LB_GETCURSEL,0,0));
}

std::wstring row(const dynsyslab::core::JobRecord& job) {
    std::wostringstream out;
    out<<std::left<<std::setw(11)<<dynsyslab::core::wide(
        dynsyslab::core::jobStatusName(job.status))<<L' '
       <<std::right<<std::setw(3)<<static_cast<int>(job.progress*100.0+0.5)<<L"%  "
       <<std::left<<std::setw(12)<<dynsyslab::core::wide(job.kind)<<L' '
       <<dynsyslab::core::wide(job.label)<<L"  ["<<dynsyslab::core::wide(job.id)<<L']';
    return out.str();
}

void showDetails(DialogState& state) {
    const int index=selection(state);
    if(index<0||static_cast<std::size_t>(index)>=state.jobs.size()) {
        SetWindowTextW(state.details,L"Select a queued job to inspect its status and exact job definition.");
        EnableWindow(state.cancel,FALSE);return;
    }
    const auto& job=state.jobs[static_cast<std::size_t>(index)];
    std::wostringstream out;
    out<<L"ID: "<<dynsyslab::core::wide(job.id)<<L"\r\nStatus: "
       <<dynsyslab::core::wide(dynsyslab::core::jobStatusName(job.status))
       <<L"   Attempts: "<<job.attempts<<L"   Worker: "
       <<dynsyslab::core::wide(job.worker)<<L"\r\nUpdated: "
       <<dynsyslab::core::wide(job.updatedUtc)<<L"\r\nMessage: "
       <<dynsyslab::core::wide(job.message)<<L"\r\n\r\n"
       <<dynsyslab::core::wide(job.manifestJson);
    SetWindowTextW(state.details,out.str().c_str());
    EnableWindow(state.cancel,job.status==dynsyslab::core::JobStatus::queued
        ||job.status==dynsyslab::core::JobStatus::running);
}

void refresh(DialogState& state,const std::wstring& message={}) {
    state.jobs=state.backend->queueJobs();SendMessageW(state.list,LB_RESETCONTENT,0,0);
    std::size_t queued=0,running=0,failed=0;
    for(const auto& job:state.jobs){const auto label=row(job);
        SendMessageW(state.list,LB_ADDSTRING,0,reinterpret_cast<LPARAM>(label.c_str()));
        queued+=job.status==dynsyslab::core::JobStatus::queued;
        running+=job.status==dynsyslab::core::JobStatus::running;
        failed+=job.status==dynsyslab::core::JobStatus::failed;}
    std::wostringstream summary;
    if(!message.empty())summary<<message<<L"  ·  ";
    summary<<state.jobs.size()<<L" recent jobs  ·  "<<queued<<L" queued  ·  "
           <<running<<L" running  ·  "<<failed<<L" failed";
    SetWindowTextW(state.summary,summary.str().c_str());
    if(!state.jobs.empty()){SendMessageW(state.list,LB_SETCURSEL,0,0);showDetails(state);}
    else showDetails(state);
}

void createChildren(DialogState& state) {
    control(state.window,L"STATIC",L"Durable native job queue",SS_LEFT,0,18,14,840,25,state.fonts->heading);
    state.summary=control(state.window,L"STATIC",L"Loading…",SS_LEFT,0,18,44,840,20,state.fonts->caption);
    state.list=control(state.window,L"LISTBOX",L"",LBS_NOTIFY|LBS_NOINTEGRALHEIGHT|WS_VSCROLL|WS_HSCROLL,
        kListId,18,70,840,220,state.fonts->mono,true);
    state.details=control(state.window,L"EDIT",L"",ES_MULTILINE|ES_AUTOVSCROLL|ES_READONLY|WS_VSCROLL,
        0,18,302,840,190,state.fonts->mono,true);
    control(state.window,L"BUTTON",L"Refresh",BS_PUSHBUTTON,kRefreshId,18,508,88,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Queue active job",BS_PUSHBUTTON,kSubmitId,116,508,130,30,state.fonts->body);
    state.run=control(state.window,L"BUTTON",L"Run next",BS_DEFPUSHBUTTON,kRunId,256,508,95,30,state.fonts->body);
    state.cancel=control(state.window,L"BUTTON",L"Cancel selected",BS_PUSHBUTTON,kCancelId,361,508,125,30,state.fonts->body);
    control(state.window,L"STATIC",L"Recover stale after",SS_LEFT,0,18,555,120,20,state.fonts->caption);
    state.age=control(state.window,L"EDIT",L"3600",ES_AUTOHSCROLL,kAgeId,140,550,76,26,state.fonts->body,true);
    control(state.window,L"STATIC",L"seconds",SS_LEFT,0,222,555,58,20,state.fonts->caption);
    control(state.window,L"BUTTON",L"Recover",BS_PUSHBUTTON,kRecoverId,287,548,90,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Close",BS_PUSHBUTTON,kCloseId,758,548,100,30,state.fonts->body);
    EnableWindow(state.run,!state.taskActive);EnableWindow(state.cancel,FALSE);
    refresh(state);
}

void reportError(HWND owner,const std::exception& error) {
    MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),L"Native queue operation failed",
                MB_OK|MB_ICONERROR);
}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(
        reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;
        SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){reportError(window,error);return -1;}return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam),note=HIWORD(wparam);
        try{
            if(id==kListId&&note==LBN_SELCHANGE){showDetails(*state);return 0;}
            if(id==kRefreshId){refresh(*state);return 0;}
            if(id==kSubmitId){const auto submitted=state->backend->submitQueued(state->activeJob);
                refresh(*state,L"Queued "+dynsyslab::core::wide(submitted));return 0;}
            if(id==kRunId){state->runNext=true;DestroyWindow(window);return 0;}
            if(id==kCancelId){const int index=selection(*state);if(index<0)return 0;
                const auto& job=state->jobs[static_cast<std::size_t>(index)];
                if(MessageBoxW(window,L"Cancel the selected queued or running job?",L"Cancel job",
                               MB_YESNO|MB_ICONWARNING)==IDYES){state->backend->cancelQueued(job.id);
                    refresh(*state,L"Cancellation recorded");}return 0;}
            if(id==kRecoverId){const double age=std::stod(utf8Text(state->age));
                const auto count=state->backend->recoverQueued(age);
                refresh(*state,L"Recovered "+std::to_wstring(count));return 0;}
            if(id==kCloseId||id==IDCANCEL){DestroyWindow(window);return 0;}
        }catch(const std::exception& error){reportError(window,error);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){
        static HBRUSH panel=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);
        SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);
        return reinterpret_cast<LRESULT>(panel);}
    return DefWindowProcW(window,message,wparam,lparam);
}

bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;
    type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);
    type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

bool showQueueDialog(HWND owner,Backend& backend,const std::string& activeJob,
                     bool taskActive,const dsl::theme::Fonts& fonts) {
    if(!registerClass())return false;DialogState state;state.backend=&backend;
    state.fonts=&fonts;state.activeJob=activeJob;state.taskActive=taskActive;
    RECT area{};GetWindowRect(owner,&area);const int width=900,height=640;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"DynSysLab native job queue",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,
        GetModuleHandleW(nullptr),&state);
    if(!state.window)return false;dsl::theme::configureTitleBar(state.window);
    EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);
    MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){
        if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.runNext;
}

}  // namespace dsl::compute
