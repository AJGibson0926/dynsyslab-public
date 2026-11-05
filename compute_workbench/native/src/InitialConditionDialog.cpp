#include "InitialConditionDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <cwctype>
#include <iomanip>
#include <sstream>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabInitialConditionDialog";
constexpr int kTrajectoryId=1001,kStateId=1002,kSetId=1003,kDoneId=1004,kCancelId=1005;
struct DialogState {
    HWND window=nullptr,trajectory=nullptr,list=nullptr,first=nullptr,second=nullptr,firstLabel=nullptr,secondLabel=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;std::string mode,current;std::size_t count=1,dimension=0;
    std::vector<std::string> labels;std::vector<double> values;std::optional<std::string> output;
};
HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,int x,int y,int width,int height,HFONT font){const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX"||std::wstring(type)==L"LISTBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;}
std::wstring number(double value){std::wostringstream out;out<<std::setprecision(12)<<value;return out.str();}
bool readNumber(HWND edit,double& result){wchar_t buffer[96]{};GetWindowTextW(edit,buffer,96);wchar_t* end=nullptr;result=std::wcstod(buffer,&end);while(end&&std::iswspace(*end))++end;return end&&end!=buffer&&*end==L'\0'&&std::isfinite(result);}
std::size_t selected(HWND list){const auto index=SendMessageW(list,LB_GETCURSEL,0,0);return index==LB_ERR?static_cast<std::size_t>(-1):static_cast<std::size_t>(index);}
std::size_t trajectory(const DialogState& state){const auto index=SendMessageW(state.trajectory,CB_GETCURSEL,0,0);return index==CB_ERR?0:static_cast<std::size_t>(index);}

void loadFields(DialogState& state){const auto axis=selected(state.list);if(axis>=state.dimension)return;
    const bool explicitMode=state.mode=="explicit";const auto base=explicitMode?(trajectory(state)*state.dimension+axis):(axis*2);
    SetWindowTextW(state.first,number(state.values[base]).c_str());
    if(!explicitMode)SetWindowTextW(state.second,number(state.values[base+1]).c_str());}
void refresh(DialogState& state,std::size_t selectedAxis=0){SendMessageW(state.list,LB_RESETCONTENT,0,0);const bool explicitMode=state.mode=="explicit";
    for(std::size_t axis=0;axis<state.dimension;++axis){const auto base=explicitMode?(trajectory(state)*state.dimension+axis):(axis*2);
        std::wstring row=dynsyslab::core::wide(state.labels[axis])+L"  ";
        row+=explicitMode?L"= "+number(state.values[base]):L"∈ ["+number(state.values[base])+L", "+number(state.values[base+1])+L"]";
        SendMessageW(state.list,LB_ADDSTRING,0,reinterpret_cast<LPARAM>(row.c_str()));}
    if(state.dimension){SendMessageW(state.list,LB_SETCURSEL,std::min(selectedAxis,state.dimension-1),0);loadFields(state);}}
void setValue(DialogState& state){const auto axis=selected(state.list);double first=0,second=0;if(axis>=state.dimension||!readNumber(state.first,first)||
    (state.mode!="explicit"&&!readNumber(state.second,second))){MessageBoxW(state.window,L"Enter finite numeric values.",L"Invalid initial condition",MB_OK|MB_ICONWARNING);return;}
    if(state.mode=="explicit")state.values[trajectory(state)*state.dimension+axis]=first;else{if(first>second){MessageBoxW(state.window,L"Lower bound must not exceed upper bound.",L"Invalid initial-condition range",MB_OK|MB_ICONWARNING);return;}state.values[axis*2]=first;state.values[axis*2+1]=second;}refresh(state,axis);}
std::string encoded(const DialogState& state){std::ostringstream out;out<<std::setprecision(17)<<'[';
    if(state.mode=="explicit"){for(std::size_t row=0;row<state.count;++row){if(row)out<<',';out<<'[';for(std::size_t axis=0;axis<state.dimension;++axis){if(axis)out<<',';out<<state.values[row*state.dimension+axis];}out<<']';}}
    else for(std::size_t axis=0;axis<state.dimension;++axis){if(axis)out<<',';out<<'['<<state.values[axis*2]<<','<<state.values[axis*2+1]<<']';}
    return out.str()+']';}

void initialize(DialogState& state,const std::vector<double>& defaults){const auto parsed=dsl::json::parse(state.current.empty()?"[]":state.current);
    if(state.mode=="explicit"){state.values.resize(state.count*state.dimension);for(std::size_t row=0;row<state.count;++row)std::copy(defaults.begin(),defaults.end(),state.values.begin()+row*state.dimension);
        const auto& rows=parsed.array();if(rows.size()==state.count&&!rows.empty()&&!rows.front().array().empty())for(std::size_t row=0;row<state.count;++row){const auto& columns=rows[row].array();if(columns.size()==state.dimension)for(std::size_t axis=0;axis<state.dimension;++axis)state.values[row*state.dimension+axis]=columns[axis].number();}
        else if(state.count==1&&rows.size()==state.dimension)for(std::size_t axis=0;axis<state.dimension;++axis)state.values[axis]=rows[axis].number();}
    else {state.values.resize(state.dimension*2);for(std::size_t axis=0;axis<state.dimension;++axis){const double width=0.1*std::max(1.0,std::abs(defaults[axis]));state.values[axis*2]=defaults[axis]-width;state.values[axis*2+1]=defaults[axis]+width;}
        const auto& rows=parsed.array();if(rows.size()==state.dimension)for(std::size_t axis=0;axis<state.dimension;++axis)if(rows[axis].array().size()==2){state.values[axis*2]=rows[axis].array()[0].number();state.values[axis*2+1]=rows[axis].array()[1].number();}}
}

void createChildren(DialogState& state){control(state.window,L"STATIC",state.mode=="explicit"?L"Explicit initial states":L"Random initial-state bounds",SS_LEFT,0,18,13,420,24,state.fonts->heading);
    control(state.window,L"STATIC",L"Trajectory",SS_LEFT,0,18,46,100,18,state.fonts->caption);state.trajectory=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,kTrajectoryId,18,65,145,180,state.fonts->body);
    for(std::size_t index=0;index<state.count;++index){const auto label=L"Trajectory "+std::to_wstring(index+1);SendMessageW(state.trajectory,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(label.c_str()));}SendMessageW(state.trajectory,CB_SETCURSEL,0,0);EnableWindow(state.trajectory,state.mode=="explicit"&&state.count>1);
    control(state.window,L"STATIC",L"State components",SS_LEFT,0,18,105,420,18,state.fonts->caption);state.list=control(state.window,L"LISTBOX",L"",LBS_NOTIFY|WS_VSCROLL,kStateId,18,124,420,175,state.fonts->mono);
    state.firstLabel=control(state.window,L"STATIC",state.mode=="explicit"?L"Value":L"Lower bound",SS_LEFT,0,18,314,126,18,state.fonts->caption);
    state.first=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL,0,18,333,150,25,state.fonts->mono);
    state.secondLabel=control(state.window,L"STATIC",L"Upper bound",SS_LEFT,0,180,314,126,18,state.fonts->caption);
    state.second=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL,0,180,333,150,25,state.fonts->mono);ShowWindow(state.secondLabel,state.mode=="explicit"?SW_HIDE:SW_SHOW);ShowWindow(state.second,state.mode=="explicit"?SW_HIDE:SW_SHOW);
    control(state.window,L"BUTTON",L"Set selected",BS_PUSHBUTTON,kSetId,340,330,98,29,state.fonts->body);
    control(state.window,L"STATIC",state.mode=="normal"?L"Normal draws use the interval midpoint and one-sixth width, clipped to these bounds.":L"The deterministic native solver creates draws only when the job runs.",SS_LEFT,0,18,376,420,38,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply",BS_DEFPUSHBUTTON,kDoneId,226,426,100,30,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,338,426,100,30,state.fonts->body);refresh(state);}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);return 0;}if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kStateId&&HIWORD(wparam)==LBN_SELCHANGE){loadFields(*state);return 0;}if(id==kTrajectoryId&&HIWORD(wparam)==CBN_SELCHANGE){refresh(*state);return 0;}if(id==kSetId){setValue(*state);return 0;}if(id==kDoneId||id==IDOK){state->output=encoded(*state);DestroyWindow(window);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showInitialConditionDialog(HWND owner,const std::string& mode,std::size_t count,
    const std::vector<std::string>& labels,const std::vector<double>& defaults,const std::string& current,const dsl::theme::Fonts& fonts){if(mode=="default")return std::string("[]");if(mode!="explicit"&&mode!="uniform"&&mode!="normal")return std::nullopt;if(labels.empty()||labels.size()!=defaults.size()||!registerClass())return std::nullopt;
    DialogState state;state.mode=mode;state.count=std::max<std::size_t>(1,count);state.dimension=labels.size();state.labels=labels;state.current=current;state.fonts=&fonts;
    try{initialize(state,defaults);}catch(const std::exception& error){MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),L"Invalid initial-condition data",MB_OK|MB_ICONERROR);return std::nullopt;}
    RECT area{};GetWindowRect(owner,&area);const int width=480,height=520;state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Initial-condition builder",WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
