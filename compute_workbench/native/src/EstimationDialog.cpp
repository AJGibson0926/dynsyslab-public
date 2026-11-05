#include "EstimationDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <utility>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabEstimationDialog";
constexpr int kDoneId=1401,kCancelId=1402,kKromId=1403;
struct DialogState {
    HWND window=nullptr,krom=nullptr,modes=nullptr,measurement=nullptr,stateMap=nullptr;
    HWND process=nullptr,measurementNoise=nullptr,processDiagonal=nullptr,measurementDiagonal=nullptr;
    HWND detectable=nullptr,method=nullptr,blocks=nullptr,horizon=nullptr,sampleTime=nullptr;
    HWND tolerance=nullptr,maximumEntry=nullptr,initialState=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;const Registry* registry=nullptr;
    std::string current,system;std::optional<std::string> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font) {
    const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,
        x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;
}
void label(DialogState& state,const wchar_t* text,int x,int y,int width) {
    control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);
}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string lower(std::string value){std::transform(value.begin(),value.end(),value.begin(),
    [](unsigned char item){return static_cast<char>(std::tolower(item));});return value;}
std::string systemKey(std::string value){value=lower(std::move(value));
    if(value.find("rayleigh")!=std::string::npos||value=="rpe")return "rpe";
    if(value.find("duffing")!=std::string::npos)return "duffing";
    if(value.find("marmottant")!=std::string::npos)return "marmottant";
    if(value.find("reddy")!=std::string::npos)return "reddy_szeri";return value;}
std::string shown(const dsl::json::Value& value,const std::string& fallback){return value.isNull()?fallback:dsl::json::dump(value);}
std::string arrayLiteral(std::string value){const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return "[]";if(value[first]!='[')value="["+value+"]";
    const auto parsed=dsl::json::parse(value);if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument("Expected a JSON array or comma-separated list.");return dsl::json::dump(parsed);}
std::vector<std::string> splitLines(const std::string& source){std::vector<std::string> result;std::size_t start=0;
    while(start<source.size()){const auto end=source.find('\n',start);result.push_back(source.substr(start,
        end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}return result;}
void refill(HWND combo,const std::vector<std::string>& names,const std::string& current){
    SendMessageW(combo,CB_RESETCONTENT,0,0);int selected=-1;
    for(std::size_t index=0;index<names.size();++index){const auto wide=dynsyslab::core::wide(names[index]);
        SendMessageW(combo,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(wide.c_str()));if(names[index]==current)selected=static_cast<int>(index);}
    if(selected>=0)SendMessageW(combo,CB_SETCURSEL,selected,0);else SetWindowTextW(combo,dynsyslab::core::wide(current).c_str());
}
void loadReconstructions(DialogState& state){const std::string measurement=utf8Text(state.measurement),stateName=utf8Text(state.stateMap);
    std::vector<std::string> names;try{const auto archive=dynsyslab::core::readNpz(
        std::filesystem::u8path(utf8Text(state.krom)));if(archive.contains("reconstruction_names"))
        names=splitLines(archive.at("reconstruction_names").string());}catch(...){ }
    if(names.empty())names={"echo","state"};refill(state.measurement,names,measurement.empty()?"echo":measurement);
    refill(state.stateMap,names,stateName.empty()?"state":stateName);
}
HWND edit(DialogState& state,const std::string& value,int x,int y,int width){return control(state.window,L"EDIT",
    dynsyslab::core::wide(value).c_str(),ES_AUTOHSCROLL,0,x,y,width,25,state.fonts->mono);}
HWND combo(DialogState& state,const std::vector<std::string>& values,const std::string& selected,
           bool editable,int id,int x,int y,int width){HWND result=control(state.window,L"COMBOBOX",L"",
    (editable?CBS_DROPDOWN:CBS_DROPDOWNLIST)|WS_VSCROLL,id,x,y,width,220,state.fonts->body);
    refill(result,values,selected);return result;}

void createChildren(DialogState& state){const auto settings=dsl::json::parse(jsonObject(state.current));
    control(state.window,L"STATIC",L"Koopman observability and estimation",SS_LEFT,0,18,12,580,24,state.fonts->heading);
    label(state,L"KROM archive",18,47,580);std::vector<std::string> kroms;const auto current=settings.at("krom").string();
    for(const auto& item:state.registry->kroms)if(item.system.empty()||systemKey(item.system)==systemKey(state.system))kroms.push_back(item.source);
    state.krom=combo(state,kroms,current,true,kKromId,18,66,580);if(current.empty()&&!kroms.empty())SendMessageW(state.krom,CB_SETCURSEL,0,0);
    label(state,L"Selected modes (optional)",18,105,182);label(state,L"Measurement C map",209,105,190);label(state,L"State reconstruction C map",408,105,190);
    state.modes=edit(state,shown(settings.at("selected_modes"),"[]"),18,124,182);
    state.measurement=combo(state,{},settings.at("measurement").string("echo"),true,0,209,124,190);
    state.stateMap=combo(state,{},settings.at("state_reconstruction").string("state"),true,0,408,124,190);loadReconstructions(state);
    label(state,L"Process variance",18,166,182);label(state,L"Measurement variance",209,166,190);label(state,L"KLQE initial physical state (optional)",408,166,190);
    state.process=edit(state,shown(settings.at("process_variance"),"1"),18,185,182);
    state.measurementNoise=edit(state,shown(settings.at("measurement_variance"),"0.0001"),209,185,190);
    state.initialState=edit(state,shown(settings.at("initial_state"),"[]"),408,185,190);
    label(state,L"Process covariance diagonal (optional)",18,227,281);label(state,L"Measurement covariance diagonal (optional)",317,227,281);
    state.processDiagonal=edit(state,shown(settings.at("process_diagonal"),"[]"),18,246,281);
    state.measurementDiagonal=edit(state,shown(settings.at("measurement_diagonal"),"[]"),317,246,281);
    label(state,L"Observability construction",18,288,281);label(state,L"Block count (0 = lifted dimension)",317,288,281);
    state.method=combo(state,{"matrix","discrete-gramian","continuous-gramian","continuous-infinite-gramian"},
        settings.at("observability_method").string("matrix"),false,0,18,307,281);
    state.blocks=edit(state,shown(settings.at("observability_blocks"),"0"),317,307,281);
    label(state,L"Finite horizon",18,349,135);label(state,L"Sample time (0 = infer)",162,349,135);
    label(state,L"Rank tolerance",306,349,135);label(state,L"Maximum matrix entry",450,349,148);
    state.horizon=edit(state,shown(settings.at("observability_horizon"),"1"),18,368,135);
    state.sampleTime=edit(state,shown(settings.at("observability_sample_time"),"0"),162,368,135);
    state.tolerance=edit(state,shown(settings.at("observability_tolerance"),"1e-10"),306,368,135);
    state.maximumEntry=edit(state,shown(settings.at("observability_maximum_entry"),"1e150"),450,368,148);
    state.detectable=control(state.window,L"BUTTON",L"Require a detectable estimator design",BS_AUTOCHECKBOX,0,18,411,320,24,state.fonts->body);
    SendMessageW(state.detectable,BM_SETCHECK,settings.at("require_detectable").boolean(true)?BST_CHECKED:BST_UNCHECKED,0);
    control(state.window,L"STATIC",L"Named reconstruction maps are the KROM C matrices. Observability archives matrices, Gramians, singular directions, estimator gain, covariance, and poles; KLQE also evaluates the selected trajectory.",SS_LEFT,0,18,449,580,44,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,367,507,118,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,495,507,103,30,state.fonts->body);
}
void accept(DialogState& state){try{const auto number=[](HWND item){return dsl::json::dump(dsl::json::parse(utf8Text(item)));};
    state.output=mergeSettingsObject(state.current,{{"krom",quoted(utf8Text(state.krom))},
        {"selected_modes",arrayLiteral(utf8Text(state.modes))},{"measurement",quoted(utf8Text(state.measurement))},
        {"state_reconstruction",quoted(utf8Text(state.stateMap))},{"process_variance",number(state.process)},
        {"measurement_variance",number(state.measurementNoise)},{"process_diagonal",arrayLiteral(utf8Text(state.processDiagonal))},
        {"measurement_diagonal",arrayLiteral(utf8Text(state.measurementDiagonal))},{"initial_state",arrayLiteral(utf8Text(state.initialState))},
        {"require_detectable",checked(state.detectable)?"true":"false"},{"observability_method",quoted(utf8Text(state.method))},
        {"observability_blocks",number(state.blocks)},{"observability_horizon",number(state.horizon)},
        {"observability_sample_time",number(state.sampleTime)},{"observability_tolerance",number(state.tolerance)},
        {"observability_maximum_entry",number(state.maximumEntry)}});DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid estimation setting",MB_OK|MB_ICONWARNING);}}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid advanced estimation settings",MB_OK|MB_ICONERROR);return -1;}return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam),notice=HIWORD(wparam);if(id==kKromId&&notice==CBN_SELENDOK){loadReconstructions(*state);return 0;}if(id==kDoneId||id==IDOK){accept(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showEstimationDialog(HWND owner,const std::filesystem::path& source,
    const Registry& registry,const std::string& current,const dsl::theme::Fonts& fonts){DialogState state;
    state.registry=&registry;state.current=current;state.fonts=&fonts;try{const auto archive=dynsyslab::core::readNpz(source);
        const auto metadata=dsl::json::parse(archive.metadataJson());state.system=metadata.at("system").string();
        if(!(archive.contains("values")||archive.contains("states"))||!archive.contains("time"))
            throw std::invalid_argument("Selected source is not a time-resolved trajectory archive.");
    }catch(const std::exception& error){MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),L"Estimation source",MB_OK|MB_ICONERROR);return std::nullopt;}
    if(!registerClass())return std::nullopt;RECT area{};GetWindowRect(owner,&area);const int width=640,height=602;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Observability / KLQE settings",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};
    while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
