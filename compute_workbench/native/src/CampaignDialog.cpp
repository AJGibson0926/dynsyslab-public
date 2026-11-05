#include "CampaignDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <sstream>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabCampaignDialog";
constexpr int kAxesId=901,kAddId=902,kRemoveId=903,kDoneId=904,kCancelId=905;
struct Axis { std::string path; std::string values; };
struct DialogState {
    HWND window=nullptr,enabled=nullptr,kind=nullptr,mode=nullptr,list=nullptr,path=nullptr,values=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;std::string current;std::vector<std::string> suggestions;
    std::vector<Axis> axes;std::optional<std::string> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font){const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX"||std::wstring(type)==L"LISTBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,
        x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;}
void label(DialogState& state,const wchar_t* text,int x,int y,int width){control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);}
std::string text(HWND window){return utf8Text(window);}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}

void refresh(DialogState& state,int selected=-1){SendMessageW(state.list,LB_RESETCONTENT,0,0);
    for(const auto& axis:state.axes){const auto row=dynsyslab::core::wide(axis.path+" = "+axis.values);
        SendMessageW(state.list,LB_ADDSTRING,0,reinterpret_cast<LPARAM>(row.c_str()));}
    if(selected>=0&&selected<static_cast<int>(state.axes.size()))SendMessageW(state.list,LB_SETCURSEL,selected,0);}
void loadSelection(DialogState& state){const auto selected=SendMessageW(state.list,LB_GETCURSEL,0,0);
    if(selected==LB_ERR)return;const auto& axis=state.axes[static_cast<std::size_t>(selected)];
    SetWindowTextW(state.path,dynsyslab::core::wide(axis.path).c_str());SetWindowTextW(state.values,dynsyslab::core::wide(axis.values).c_str());}

std::string normalizedValues(std::string value){const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return "[]";if(value[first]!='[')value="["+value+"]";
    const auto parsed=dsl::json::parse(value);if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage())||parsed.array().empty())
        throw std::invalid_argument("Campaign values must be a non-empty JSON array or comma-separated list.");
    return dsl::json::dump(parsed);}

void updateAxis(DialogState& state){try{const std::string path=text(state.path);if(path.empty())throw std::invalid_argument("Campaign axis path is empty.");
    Axis axis{path,normalizedValues(text(state.values))};const auto selected=SendMessageW(state.list,LB_GETCURSEL,0,0);
    int index=static_cast<int>(state.axes.size());if(selected!=LB_ERR){index=static_cast<int>(selected);state.axes[static_cast<std::size_t>(selected)]=std::move(axis);}
    else state.axes.push_back(std::move(axis));refresh(state,index);}catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid campaign axis",MB_OK|MB_ICONWARNING);}}

void accept(DialogState& state){std::string axes="[";for(std::size_t index=0;index<state.axes.size();++index){if(index)axes+=',';
    axes+="{\"path\":"+quoted(state.axes[index].path)+",\"values\":"+state.axes[index].values+"}";}axes+="]";
    state.output=mergeSettingsObject(state.current,{{"enabled",checked(state.enabled)?"true":"false"},
        {"kind",quoted(text(state.kind))},{"mode",quoted(text(state.mode))},{"axes",axes}});DestroyWindow(state.window);}

void createChildren(DialogState& state){const auto current=dsl::json::parse(jsonObject(state.current));
    control(state.window,L"STATIC",L"Cartesian campaign",SS_LEFT,0,18,13,400,24,state.fonts->heading);
    state.enabled=control(state.window,L"BUTTON",L"Enabled",BS_AUTOCHECKBOX,0,18,44,90,24,state.fonts->body);
    SendMessageW(state.enabled,BM_SETCHECK,current.at("enabled").boolean(false)?BST_CHECKED:BST_UNCHECKED,0);
    label(state,L"Campaign kind",120,40,160);state.kind=control(state.window,L"EDIT",dynsyslab::core::wide(current.at("kind").string("parameter_sweep")).c_str(),ES_AUTOHSCROLL,0,120,59,190,25,state.fonts->body);
    label(state,L"Combination mode",322,40,125);state.mode=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,0,322,59,125,150,state.fonts->body);
    SendMessageW(state.mode,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"cartesian"));SendMessageW(state.mode,CB_SETCURSEL,0,0);
    label(state,L"Sweep axes",18,96,429);state.list=control(state.window,L"LISTBOX",L"",LBS_NOTIFY|WS_VSCROLL,kAxesId,18,115,429,145,state.fonts->mono);
    for(const auto& item:current.at("axes").array())state.axes.push_back({item.at("path").string(),dsl::json::dump(item.at("values"))});refresh(state);
    label(state,L"Workflow path",18,273,205);label(state,L"Values (JSON array or comma-separated)",232,273,215);
    state.path=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWN|WS_VSCROLL,0,18,292,205,220,state.fonts->body);
    for(const auto& suggestion:state.suggestions)SendMessageW(state.path,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(dynsyslab::core::wide(suggestion).c_str()));
    state.values=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,232,292,215,25,state.fonts->mono);
    control(state.window,L"BUTTON",L"Add / update axis",BS_PUSHBUTTON,kAddId,18,331,138,29,state.fonts->body);
    control(state.window,L"BUTTON",L"Remove selected",BS_PUSHBUTTON,kRemoveId,164,331,125,29,state.fonts->body);
    control(state.window,L"STATIC",L"Each point is archived under one campaign identity; system parameter axes create distinct equation fingerprints.",SS_LEFT,0,18,376,429,38,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply campaign",BS_DEFPUSHBUTTON,kDoneId,217,427,119,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,345,427,102,30,state.fonts->body);}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid campaign settings",MB_OK|MB_ICONERROR);return -1;}return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kAxesId&&HIWORD(wparam)==LBN_SELCHANGE){loadSelection(*state);return 0;}if(id==kAddId){updateAxis(*state);return 0;}
        if(id==kRemoveId){const auto selected=SendMessageW(state->list,LB_GETCURSEL,0,0);if(selected!=LB_ERR){state->axes.erase(state->axes.begin()+selected);refresh(*state);}return 0;}
        if(id==kDoneId||id==IDOK){accept(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showCampaignDialog(HWND owner,const std::string& current,
    const std::vector<std::string>& suggestions,const dsl::theme::Fonts& fonts){if(!registerClass())return std::nullopt;
    DialogState state;state.current=current;state.suggestions=suggestions;state.fonts=&fonts;RECT area{};GetWindowRect(owner,&area);const int width=490,height=520;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Campaign builder",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};
    while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
