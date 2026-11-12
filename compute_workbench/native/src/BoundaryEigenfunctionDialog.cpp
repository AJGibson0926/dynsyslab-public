#include "BoundaryEigenfunctionDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabBoundaryEigenfunctionDialog";
constexpr int kDoneId=1601,kCancelId=1602;

struct State {
    HWND window=nullptr,spectrum=nullptr,observable=nullptr,additional=nullptr;
    HWND groups=nullptr,neighbors=nullptr,step=nullptr,tolerance=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;const Registry* registry=nullptr;
    std::string current,system;std::optional<std::string> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,
             int id,int x,int y,int width,int height,HFONT font){
    const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,
        WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;
}
void label(State& state,const wchar_t* text,int x,int y,int width){
    control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);
}
std::string lower(std::string value){std::transform(value.begin(),value.end(),value.begin(),
    [](unsigned char item){return static_cast<char>(std::tolower(item));});return value;}
std::string systemKey(std::string value){value=lower(std::move(value));
    if(value.find("rayleigh")!=std::string::npos||value=="rpe")return "rpe";
    if(value.find("duffing")!=std::string::npos)return "duffing";
    if(value.find("marmottant")!=std::string::npos)return "marmottant";
    if(value.find("reddy")!=std::string::npos)return "reddy_szeri";return value;}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string array(const std::string& value,const char* field){
    const auto parsed=dsl::json::parse(value.empty()?"[]":value);
    if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument(std::string(field)+" must be a JSON array");
    return dsl::json::dump(parsed);
}

void createChildren(State& state){
    const auto settings=dsl::json::parse(jsonObject(state.current));
    control(state.window,L"STATIC",L"Construct eigenfunctions from a supplied spectrum",
        SS_LEFT,0,18,13,600,24,state.fonts->heading);
    label(state,L"Immutable eigenvalue spectrum",18,48,600);
    state.spectrum=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWN|WS_VSCROLL,
        0,18,67,600,240,state.fonts->body);
    const auto selected=settings.at("spectrum").string();int choice=-1,index=0;
    for(const auto& item:state.registry->spectra)
        if(item.system.empty()||systemKey(item.system)==systemKey(state.system)){
            const auto path=dynsyslab::core::wide(item.source);
            SendMessageW(state.spectrum,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(path.c_str()));
            if(item.source==selected)choice=index;++index;}
    if(choice>=0)SendMessageW(state.spectrum,CB_SETCURSEL,choice,0);
    else SetWindowTextW(state.spectrum,dynsyslab::core::wide(selected).c_str());
    label(state,L"Observable to reconstruct",18,108,290);
    state.observable=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,
        0,18,127,290,180,state.fonts->body);
    auto system=dynsyslab::core::createSystem([&]{dynsyslab::core::SystemConfiguration value;
        value.identifier=systemKey(state.system);return value;}());
    const auto observable=settings.at("observable").string("state_0");choice=-1;index=0;
    for(const auto& name:system->observableCatalog().names()){
        const auto wide=dynsyslab::core::wide(name);SendMessageW(state.observable,CB_ADDSTRING,0,
            reinterpret_cast<LPARAM>(wide.c_str()));if(name==observable)choice=index;++index;}
    SendMessageW(state.observable,CB_SETCURSEL,choice>=0?choice:0,0);
    label(state,L"Additional trajectory archives (JSON array)",328,108,290);
    state.additional=control(state.window,L"EDIT",dynsyslab::core::wide(
        dsl::json::dump(settings.at("additional_sources").isNull()
            ?dsl::json::parse("[]"):settings.at("additional_sources"))).c_str(),
        ES_AUTOHSCROLL,0,328,127,290,25,state.fonts->mono);
    label(state,L"Mode groups (blank [] assigns all modes to a scalar observable)",18,168,600);
    state.groups=control(state.window,L"EDIT",dynsyslab::core::wide(
        dsl::json::dump(settings.at("mode_groups").isNull()
            ?dsl::json::parse("[]"):settings.at("mode_groups"))).c_str(),
        ES_AUTOHSCROLL,0,18,187,600,25,state.fonts->mono);
    label(state,L"Neighbors",18,228,135);label(state,L"Gradient relative step",168,228,210);
    label(state,L"SVD tolerance",393,228,225);
    state.neighbors=control(state.window,L"EDIT",dynsyslab::core::wide(std::to_string(
        static_cast<std::size_t>(settings.at("interpolation_neighbors").number(8)))).c_str(),
        ES_AUTOHSCROLL,0,18,247,135,25,state.fonts->mono);
    state.step=control(state.window,L"EDIT",dynsyslab::core::wide(std::to_string(
        settings.at("gradient_relative_step").number(1.0e-6))).c_str(),ES_AUTOHSCROLL,
        0,168,247,210,25,state.fonts->mono);
    state.tolerance=control(state.window,L"EDIT",dynsyslab::core::wide(std::to_string(
        settings.at("svd_tolerance").number(1.0e-12))).c_str(),ES_AUTOHSCROLL,
        0,393,247,225,25,state.fonts->mono);
    control(state.window,L"STATIC",L"Only trajectory boundary amplitudes are fitted. Eigenvalues are never changed. The output is an independent scattered-interpolation KROM with a reconstruction C map and finite-difference nabla psi.",
        SS_LEFT,0,18,292,600,48,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,
        387,350,118,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,
        515,350,103,30,state.fonts->body);
}

void accept(State& state){try{
    state.output=mergeSettingsObject(state.current,{
        {"spectrum",quoted(utf8Text(state.spectrum))},
        {"observable",quoted(selectedText(state.observable))},
        {"additional_sources",array(utf8Text(state.additional),"Additional sources")},
        {"mode_groups",array(utf8Text(state.groups),"Mode groups")},
        {"interpolation_neighbors",utf8Text(state.neighbors)},
        {"gradient_relative_step",utf8Text(state.step)},
        {"svd_tolerance",utf8Text(state.tolerance)}});DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,
        dynsyslab::core::wide(error.what()).c_str(),L"Invalid boundary setting",
        MB_OK|MB_ICONWARNING);}}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){
    auto* state=reinterpret_cast<State*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<State*>(
        reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;
        SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);
        if(id==kDoneId||id==IDOK){accept(*state);return 0;}
        if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    return DefWindowProcW(window,message,wparam,lparam);
}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;
    type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);
    type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showBoundaryEigenfunctionDialog(
    HWND owner,const std::filesystem::path& source,const Registry& registry,
    const std::string& current,const dsl::theme::Fonts& fonts){
    State state;state.registry=&registry;state.current=current;state.fonts=&fonts;
    try{const auto archive=dynsyslab::core::readNpz(source);
        const auto metadata=dsl::json::parse(archive.metadataJson());
        state.system=metadata.at("system").string();
        if(!archive.contains("time")||(!archive.contains("values")&&!archive.contains("states")))
            throw std::invalid_argument("Selected source is not a trajectory archive.");
    }catch(const std::exception& error){MessageBoxW(owner,
        dynsyslab::core::wide(error.what()).c_str(),L"Boundary fitting source",
        MB_OK|MB_ICONERROR);return std::nullopt;}
    if(!registerClass())return std::nullopt;RECT area{};GetWindowRect(owner,&area);
    const int width=660,height=446;state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,
        kClassName,L"Fitted boundary eigenfunctions",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,
        width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);
    EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);
    MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0)
        if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute
