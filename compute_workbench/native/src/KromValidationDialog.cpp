#include "KromValidationDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <utility>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabKromValidationDialog";
constexpr int kDoneId=1501,kCancelId=1502;

struct DialogState {
    HWND window=nullptr,krom=nullptr,role=nullptr,horizons=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;
    const Registry* registry=nullptr;
    std::string current,system;
    std::optional<std::string> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,
             int id,int x,int y,int width,int height,HFONT font) {
    const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,
        WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;
}

void label(DialogState& state,const wchar_t* text,int x,int y,int width) {
    control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);
}

std::string lower(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char item) {
        return static_cast<char>(std::tolower(item));});return value;
}

std::string systemKey(std::string value) {
    value=lower(std::move(value));
    if(value.find("rayleigh")!=std::string::npos||value=="rpe")return "rpe";
    if(value.find("duffing")!=std::string::npos)return "duffing";
    if(value.find("marmottant")!=std::string::npos)return "marmottant";
    if(value.find("reddy")!=std::string::npos)return "reddy_szeri";
    return value;
}

std::string arrayLiteral(std::string value) {
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return "[]";
    if(value[first]!='[')value="["+value+"]";
    const auto parsed=dsl::json::parse(value);
    if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument("Horizons must be a JSON array or comma-separated list.");
    return dsl::json::dump(parsed);
}

std::string quoted(const std::string& value) {return "\""+jsonEscape(value)+"\"";}

void createChildren(DialogState& state) {
    const auto settings=dsl::json::parse(jsonObject(state.current));
    control(state.window,L"STATIC",L"Validate a KROM on a free trajectory",SS_LEFT,
        0,18,13,510,24,state.fonts->heading);
    label(state,L"KROM archive",18,49,510);
    state.krom=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWN|WS_VSCROLL,
        0,18,68,510,220,state.fonts->body);
    const auto current=settings.at("krom").string();int selection=-1,index=0;
    for(const auto& item:state.registry->kroms)
        if(item.system.empty()||systemKey(item.system)==systemKey(state.system)) {
            const auto path=dynsyslab::core::wide(item.source);
            SendMessageW(state.krom,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(path.c_str()));
            if(item.source==current)selection=index;++index;
        }
    if(selection>=0)SendMessageW(state.krom,CB_SETCURSEL,selection,0);
    else SetWindowTextW(state.krom,dynsyslab::core::wide(current).c_str());
    label(state,L"Dataset role",18,108,245);label(state,L"Semigroup horizons (blank = automatic)",273,108,255);
    state.role=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,
        0,18,127,245,120,state.fonts->body);
    for(const auto* value:{L"training",L"test"})SendMessageW(state.role,CB_ADDSTRING,0,
        reinterpret_cast<LPARAM>(value));
    SendMessageW(state.role,CB_SETCURSEL,
        settings.at("dataset_role").string("test")=="training"?0:1,0);
    state.horizons=control(state.window,L"EDIT",
        dynsyslab::core::wide(settings.at("horizons").isNull()?"[]":
            dsl::json::dump(settings.at("horizons"))).c_str(),ES_AUTOHSCROLL,
        0,273,127,255,25,state.fonts->mono);
    control(state.window,L"STATIC",
        L"No eigenfunction is refit. The result records generator, one-step, anchored long-term, and exact-horizon semigroup defects for every retained mode, plus extrapolation beyond the KROM training domain.",
        SS_LEFT,0,18,174,510,64,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,
        297,256,118,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,
        425,256,103,30,state.fonts->body);
}

void accept(DialogState& state) {
    try {
        state.output=mergeSettingsObject(state.current,{
            {"krom",quoted(utf8Text(state.krom))},
            {"dataset_role",quoted(selectedText(state.role))},
            {"horizons",arrayLiteral(utf8Text(state.horizons))}});
        DestroyWindow(state.window);
    } catch(const std::exception& error) {
        MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),
            L"Invalid validation setting",MB_OK|MB_ICONWARNING);
    }
}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(
        reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;
        SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);
        if(id==kDoneId||id==IDOK){accept(*state);return 0;}
        if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){
        static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);
        HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);
        SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}
    return DefWindowProcW(window,message,wparam,lparam);
}

bool registerClass() {
    WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;
    type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);
    type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

std::optional<std::string> showKromValidationDialog(
    HWND owner,const std::filesystem::path& source,const Registry& registry,
    const std::string& current,const dsl::theme::Fonts& fonts) {
    DialogState state;state.registry=&registry;state.current=current;state.fonts=&fonts;
    try {const auto archive=dynsyslab::core::readNpz(source);
        const auto metadata=dsl::json::parse(archive.metadataJson());
        state.system=metadata.at("system").string();
        if(!archive.contains("time")||(!archive.contains("values")&&!archive.contains("states")))
            throw std::invalid_argument("Selected source is not a trajectory archive.");
    } catch(const std::exception& error) {
        MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),
            L"KROM validation source",MB_OK|MB_ICONERROR);return std::nullopt;
    }
    if(!registerClass())return std::nullopt;RECT area{};GetWindowRect(owner,&area);
    const int width=570,height=352;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"KROM validation",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,
        GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);
    EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);
    MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0)
        if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute
