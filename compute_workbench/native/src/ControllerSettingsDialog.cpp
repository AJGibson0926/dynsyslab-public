#include "ControllerSettingsDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <array>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabControllerSettingsDialog";
constexpr int kDoneId=1101,kCancelId=1102;
struct Field { std::string key,label,value;bool arrayLike=false; };
struct DialogState {
    HWND window=nullptr,coordinates=nullptr;const dsl::theme::Fonts* fonts=nullptr;
    std::string kind,current,coordinateMode="physical";std::vector<Field> fields;std::vector<HWND> edits;
    std::optional<std::string> output;
};
HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,int x,int y,int width,int height,HFONT font){const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string compact(const dsl::json::Value& value,const std::string& fallback){return value.isNull()?fallback:dsl::json::dump(value);}
std::string valueOr(const dsl::json::Value& settings,const std::initializer_list<const char*>& keys,const std::string& fallback){for(const auto* key:keys)if(!settings.at(key).isNull())return compact(settings.at(key),fallback);return fallback;}
std::string listLiteral(std::string value){const auto first=value.find_first_not_of(" \t\r\n");if(first==std::string::npos)return "[]";if(value[first]!='[')value="["+value+"]";const auto parsed=dsl::json::parse(value);
    if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))throw std::invalid_argument("Expected a JSON array or comma-separated list.");return dsl::json::dump(parsed);}
std::string literal(std::string value,const std::string& fallback,bool arrayLike){const auto first=value.find_first_not_of(" \t\r\n");if(first==std::string::npos)return fallback;return arrayLike?listLiteral(value):dsl::json::dump(dsl::json::parse(value));}

void addField(DialogState& state,std::string key,std::string label,std::string value,bool arrayLike=false){state.fields.push_back({std::move(key),std::move(label),std::move(value),arrayLike});}
void configure(DialogState& state,const dsl::json::Value& settings){const bool prescribed=state.kind=="waveform"||state.kind=="random_pulse";
    if(prescribed){const bool model=!settings.at("model_amplitude").isNull()||!settings.at("model_amplitudes").isNull();state.coordinateMode=model?"model":"physical";}
    if(state.kind=="waveform"){
        addField(state,"amplitudes","Amplitudes",valueOr(settings,{"amplitudes","model_amplitudes","amplitude","model_amplitude"},"[1]"),true);
        addField(state,"frequencies","Frequencies",valueOr(settings,{"frequencies_hz","frequencies","frequency_hz","frequency"},"[1]"),true);
        addField(state,"phases_rad","Phases (radians)",valueOr(settings,{"phases_rad","phase_radians","phase"},"[0]"),true);
        addField(state,"offset","Offset",valueOr(settings,{"offset","model_offset"},"0"));
        addField(state,"active_window","Active time window [start, stop]",valueOr(settings,{"active_time_window_seconds","model_active_time_window","active_time_window"},"[]"),true);
    } else if(state.kind=="random_pulse"){
        addField(state,"amplitude","Pulse amplitude",valueOr(settings,{"amplitude","model_amplitude"},"1"));
        addField(state,"duration","Pulse duration",valueOr(settings,{"pulse_duration_seconds","model_pulse_duration","pulse_duration","interval"},"0.1"));
        addField(state,"seed","Controller seed",valueOr(settings,{"seed"},"0"));
        addField(state,"active_window","Active time window [start, stop]",valueOr(settings,{"active_time_window_seconds","model_active_time_window","active_time_window"},"[]"),true);
    } else if(state.kind=="pid"){
        addField(state,"Kp","Kp (scalar or input-by-state row)",valueOr(settings,{"Kp"},"1"));
        addField(state,"Ki","Ki",valueOr(settings,{"Ki"},"0"));addField(state,"Kd","Kd",valueOr(settings,{"Kd"},"0"));
        addField(state,"active_time_window","Active model-time window",valueOr(settings,{"active_time_window"},"[]"),true);
    } else {
        if(state.kind=="lqr")addField(state,"linearization_state","Linearization state",valueOr(settings,{"linearization_state"},"[]"),true);
        if(state.kind=="klqr"||state.kind=="kmpc"){
            addField(state,"selected_modes","Optional source mode indices",valueOr(settings,{"selected_modes"},"[]"),true);
            addField(state,"control_map_archive","Optional external control-map archive",valueOr(settings,{"control_map_archive"},"\"\""));
            addField(state,"control_map_prefix","External map prefix",valueOr(settings,{"control_map_prefix"},"\"\""));
        }
        addField(state,"active_time_window","Active model-time window",valueOr(settings,{"active_time_window"},"[]"),true);
    }
}

void createChildren(DialogState& state){const auto settings=dsl::json::parse(jsonObject(state.current));configure(state,settings);
    const auto title=dynsyslab::core::wide(state.kind+" controller settings");control(state.window,L"STATIC",title.c_str(),SS_LEFT,0,18,13,460,24,state.fonts->heading);
    int y=48;if(state.kind=="waveform"||state.kind=="random_pulse"){control(state.window,L"STATIC",L"Input/time coordinates",SS_LEFT,0,18,y,170,18,state.fonts->caption);state.coordinates=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,0,200,y-3,220,140,state.fonts->body);
        SendMessageW(state.coordinates,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"physical (Pa, Hz, seconds when declared)"));SendMessageW(state.coordinates,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"model / dimensionless"));SendMessageW(state.coordinates,CB_SETCURSEL,state.coordinateMode=="model"?1:0,0);y+=42;}
    for(const auto& field:state.fields){control(state.window,L"STATIC",dynsyslab::core::wide(field.label).c_str(),SS_LEFT,0,18,y,460,18,state.fonts->caption);state.edits.push_back(control(state.window,L"EDIT",dynsyslab::core::wide(field.value).c_str(),ES_AUTOHSCROLL,0,18,y+19,460,25,state.fonts->mono));y+=55;}
    control(state.window,L"STATIC",L"Structured values replace same-named advanced settings; unrelated keys are preserved.",SS_LEFT,0,18,y+2,460,30,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,249,y+43,118,30,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,376,y+43,102,30,state.fonts->body);}

void accept(DialogState& state){try{std::vector<std::pair<std::string,std::string>> values;const bool model=state.coordinates&&SendMessageW(state.coordinates,CB_GETCURSEL,0,0)==1;
    for(std::size_t index=0;index<state.fields.size();++index){const auto& field=state.fields[index];const auto raw=utf8Text(state.edits[index]);
        if(state.kind=="waveform"){
            if(field.key=="amplitudes"){values.push_back({model?"model_amplitudes":"amplitudes",literal(raw,"[]",true)});
                for(const auto* key:{"amplitude","model_amplitude",model?"amplitudes":"model_amplitudes"})values.push_back({key,"null"});continue;}
            if(field.key=="frequencies"){values.push_back({model?"frequencies":"frequencies_hz",literal(raw,"[]",true)});
                for(const auto* key:{"frequency","frequency_hz",model?"frequencies_hz":"frequencies"})values.push_back({key,"null"});continue;}
            if(field.key=="offset"){values.push_back({model?"model_offset":"offset",literal(raw,"0",false)});values.push_back({model?"offset":"model_offset","null"});continue;}
        }
        if(state.kind=="random_pulse"&&field.key=="amplitude"){values.push_back({model?"model_amplitude":"amplitude",literal(raw,"1",false)});values.push_back({model?"amplitude":"model_amplitude","null"});continue;}
        if(state.kind=="random_pulse"&&field.key=="duration"){values.push_back({model?"model_pulse_duration":"pulse_duration_seconds",literal(raw,"0.1",false)});values.push_back({model?"pulse_duration_seconds":"model_pulse_duration","null"});continue;}
        if(field.key=="active_window"){values.push_back({model?"model_active_time_window":"active_time_window_seconds",literal(raw,"[]",true)});values.push_back({model?"active_time_window_seconds":"model_active_time_window","null"});continue;}
        values.push_back({field.key,literal(raw,field.arrayLike?"[]":"null",field.arrayLike)});
    }
    state.output=mergeSettingsObject(state.current,values);DestroyWindow(state.window);}catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid controller setting",MB_OK|MB_ICONWARNING);}}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid advanced controller settings",MB_OK|MB_ICONERROR);return -1;}return 0;}if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kDoneId||id==IDOK){accept(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showControllerSettingsDialog(HWND owner,const std::string& kind,const std::string& current,const dsl::theme::Fonts& fonts){if(kind=="none"||!registerClass())return std::nullopt;DialogState state;state.kind=kind;state.current=current;state.fonts=&fonts;RECT area{};GetWindowRect(owner,&area);const int width=520,height=520;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Controller configuration",WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
