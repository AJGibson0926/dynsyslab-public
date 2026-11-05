#include "ArchiveToolsDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabArchiveToolsDialog";
constexpr int kAuditId=1801,kPreviewId=1802,kApplyId=1803,kCloseId=1804;

struct DialogState {
    HWND window{},days{},sections{},systems{},derived{};
    const dsl::theme::Fonts* fonts{};
    std::optional<ArchiveToolRequest> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font,bool edge=false){
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,
        x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;}

std::string trim(std::string value){const auto visible=[](unsigned char item){return !std::isspace(item);};
    value.erase(value.begin(),std::find_if(value.begin(),value.end(),visible));
    value.erase(std::find_if(value.rbegin(),value.rend(),visible).base(),value.end());return value;}

std::vector<std::string> names(HWND edit){std::stringstream input(utf8Text(edit));std::string item;
    std::vector<std::string> result;while(std::getline(input,item,',')){item=trim(item);if(!item.empty())result.push_back(item);}return result;}

dynsyslab::core::ArchiveAuditScope scope(const DialogState& state){
    dynsyslab::core::ArchiveAuditScope result;const auto days=trim(utf8Text(state.days));
    result.modifiedWithinDays=days.empty()?-1.0:std::stod(days);
    const auto sections=names(state.sections);if(!sections.empty())result.sections=sections;
    result.systems=names(state.systems);
    result.includeDerived=SendMessageW(state.derived,BM_GETCHECK,0,0)==BST_CHECKED;
    return result;
}

void createChildren(DialogState& state){
    control(state.window,L"STATIC",L"Canonical archive maintenance",SS_LEFT,0,18,16,560,25,state.fonts->heading);
    control(state.window,L"STATIC",L"Every inspection is read-only. Apply moves only verified, unreferenced exact duplicates into archive-local recovery storage.",SS_LEFT,0,18,48,560,38,state.fonts->caption);
    control(state.window,L"STATIC",L"Modified within days (blank = all ages)",SS_LEFT,0,18,100,270,20,state.fonts->caption);
    state.days=control(state.window,L"EDIT",L"2",ES_AUTOHSCROLL,0,18,121,270,27,state.fonts->body,true);
    control(state.window,L"STATIC",L"Archive sections (comma separated)",SS_LEFT,0,308,100,270,20,state.fonts->caption);
    state.sections=control(state.window,L"EDIT",L"solutions",ES_AUTOHSCROLL,0,308,121,270,27,state.fonts->body,true);
    control(state.window,L"STATIC",L"Systems (blank = all)",SS_LEFT,0,18,164,270,20,state.fonts->caption);
    state.systems=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL,0,18,185,270,27,state.fonts->body,true);
    state.derived=control(state.window,L"BUTTON",L"Include derived archives",BS_AUTOCHECKBOX,0,308,184,270,28,state.fonts->body);
    control(state.window,L"STATIC",L"Audit reports exact duplicates, same-plan groups, unreadable records, and byte totals. Preview applies the stricter reference-safety gate without moving files.",SS_LEFT,0,18,232,560,46,state.fonts->caption);
    control(state.window,L"BUTTON",L"Run read-only audit",BS_DEFPUSHBUTTON,kAuditId,18,296,160,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Preview deduplication",BS_PUSHBUTTON,kPreviewId,190,296,174,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Apply recoverable move…",BS_PUSHBUTTON,kApplyId,376,296,202,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Close",BS_PUSHBUTTON,kCloseId,478,350,100,30,state.fonts->body);
}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){
    auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);try{
        if(id==kCloseId||id==IDCANCEL){DestroyWindow(window);return 0;}
        if(id==kAuditId||id==kPreviewId||id==kApplyId){
            if(id==kApplyId&&MessageBoxW(window,
                L"Move verified, unreferenced exact duplicates into archive-local recovery storage?\n\nNo file is permanently deleted.",
                L"Apply recoverable deduplication",MB_YESNO|MB_ICONWARNING)!=IDYES)return 0;
            ArchiveToolRequest request;request.scope=scope(*state);
            request.action=id==kAuditId?ArchiveToolRequest::Action::audit:
                id==kPreviewId?ArchiveToolRequest::Action::deduplicatePreview:
                ArchiveToolRequest::Action::deduplicateApply;
            state->output=std::move(request);DestroyWindow(window);return 0;}
        }catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),
            L"Invalid archive scope",MB_OK|MB_ICONERROR);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){static HBRUSH panel=CreateSolidBrush(dsl::theme::panel);
        HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);
        return reinterpret_cast<LRESULT>(panel);}
    return DefWindowProcW(window,message,wparam,lparam);
}

bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;
    type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);
    type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<ArchiveToolRequest> showArchiveToolsDialog(HWND owner,const dsl::theme::Fonts& fonts){
    if(!registerClass())return std::nullopt;DialogState state;state.fonts=&fonts;
    RECT area{};GetWindowRect(owner,&area);const int width=620,height=440;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"DynSysLab archive tools",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,
        GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);
    EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);
    MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){
        if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute
