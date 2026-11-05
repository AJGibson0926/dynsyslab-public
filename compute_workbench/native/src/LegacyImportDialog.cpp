#include "LegacyImportDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "SystemParameterDialog.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Text.hpp>

#include <commdlg.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace dsl::compute {
namespace {
namespace json=dynsyslab::core::json;
constexpr wchar_t kClassName[]=L"DynSysLabLegacyImportDialog";
constexpr int kBrowseId=2001,kPreviewId=2002,kMappingId=2003,kPresetId=2004;
constexpr int kParametersId=2005,kImportId=2006,kCancelId=2007;

struct DialogState {
    HWND window{},source{},summary{},label{},strict{},inputScale{},coordinates{};
    HWND mapping{},preset{},parameters{},parameterButton{};
    const Registry* registry{};const dsl::theme::Fonts* fonts{};
    std::vector<std::size_t> systems;std::optional<dynsyslab::core::LegacyCsvBatch> decoded;
    std::optional<LegacyImportRequest> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font) {const std::wstring kind(type);const bool edge=kind==L"EDIT"||kind==L"COMBOBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;}
void caption(DialogState& state,const wchar_t* value,int x,int y,int width){control(state.window,L"STATIC",value,SS_LEFT,0,x,y,width,18,state.fonts->caption);}
int choice(HWND combo){return static_cast<int>(SendMessageW(combo,CB_GETCURSEL,0,0));}

std::filesystem::path chooseCsv(HWND owner){wchar_t path[32768]{};const wchar_t filter[]=L"Legacy tagged CSV (*.csv)\0*.csv\0All files (*.*)\0*.*\0\0";OPENFILENAMEW request{};request.lStructSize=sizeof(request);request.hwndOwner=owner;request.lpstrFile=path;request.nMaxFile=32768;request.lpstrFilter=filter;request.nFilterIndex=1;request.lpstrTitle=L"Import legacy DynSysLab CSV";request.Flags=OFN_EXPLORER|OFN_NOCHANGEDIR|OFN_FILEMUSTEXIST;return GetOpenFileNameW(&request)?std::filesystem::path(path):std::filesystem::path{};}

std::string names(const json::Value& values){std::string result;for(const auto& value:values.array()){if(!result.empty())result+=", ";
        if(const auto* text=std::get_if<std::string>(&value.storage()))result+=*text;else result+=json::dump(value);}return result;}

void preview(DialogState& state){try{const auto path=std::filesystem::path(windowText(state.source));if(path.empty())throw std::invalid_argument("Select a CSV file first");const bool strict=checked(state.strict);state.decoded=dynsyslab::core::readLegacyDynBatchCsv(path,strict);const auto summary=state.decoded->summary();std::ostringstream out;out<<"Decoded "<<static_cast<std::size_t>(summary.at("decoded_system_count").number())<<" legacy system nodes, "<<state.decoded->trajectoryCount()<<" trajectories, "<<static_cast<std::size_t>(summary.at("controlled_system_count").number())<<" controlled nodes.\r\nDynamics: "<<names(summary.at("dynamics_names"))<<"\r\nState dimensions: "<<names(summary.at("state_dimensions"))<<"\r\nSHA-256: "<<state.decoded->sourceSha256;
        if(!state.decoded->warnings.empty())out<<"\r\nWarnings: "<<state.decoded->warnings.size();SetWindowTextW(state.summary,dynsyslab::core::wide(out.str()).c_str());if(utf8Text(state.label).empty())SetWindowTextW(state.label,path.stem().c_str());}
    catch(const std::exception& error){state.decoded.reset();SetWindowTextW(state.summary,dynsyslab::core::wide(error.what()).c_str());MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Legacy CSV preview failed",MB_OK|MB_ICONWARNING);}}

void updatePresets(DialogState& state){std::vector<std::string> values{"model default"};const int selected=choice(state.mapping);if(selected>0&&static_cast<std::size_t>(selected-1)<state.systems.size()){const auto& system=state.registry->systems[state.systems[static_cast<std::size_t>(selected-1)]];for(const auto& preset:system.presets)values.push_back(preset.key);EnableWindow(state.preset,TRUE);EnableWindow(state.parameterButton,TRUE);SetWindowTextW(state.parameters,L"{}");const auto coordinates=dynsyslab::core::wide(system.preferredCoordinates);const auto coordinateIndex=SendMessageW(state.coordinates,CB_FINDSTRINGEXACT,static_cast<WPARAM>(-1),reinterpret_cast<LPARAM>(coordinates.c_str()));if(coordinateIndex!=CB_ERR)SendMessageW(state.coordinates,CB_SETCURSEL,coordinateIndex,0);EnableWindow(state.coordinates,FALSE);}else{EnableWindow(state.preset,FALSE);EnableWindow(state.parameterButton,FALSE);EnableWindow(state.coordinates,TRUE);}setComboItems(state.preset,values);}

void editParameters(DialogState& state){const int selected=choice(state.mapping);if(selected<=0||static_cast<std::size_t>(selected-1)>=state.systems.size())return;const auto& system=state.registry->systems[state.systems[static_cast<std::size_t>(selected-1)]];const int preset=choice(state.preset);const std::string key=preset>0&&static_cast<std::size_t>(preset-1)<system.presets.size()?system.presets[static_cast<std::size_t>(preset-1)].key:"";const auto result=showSystemParameterDialog(state.window,system.name,key,utf8Text(state.parameters),*state.fonts);if(result)SetWindowTextW(state.parameters,dynsyslab::core::wide(*result).c_str());}

void accept(DialogState& state){try{const auto source=std::filesystem::path(windowText(state.source));if(!state.decoded||std::filesystem::absolute(source)!=state.decoded->source)preview(state);if(!state.decoded)throw std::invalid_argument("The selected CSV has not passed preview");
        LegacyImportRequest request;request.source=source;request.options.label=utf8Text(state.label);if(request.options.label.empty())throw std::invalid_argument("Import label is required");request.options.strict=checked(state.strict);request.options.inputScale=std::stod(utf8Text(state.inputScale));if(!std::isfinite(request.options.inputScale))throw std::invalid_argument("Input scale must be finite");request.options.coordinateSystem=selectedText(state.coordinates);
        const int selected=choice(state.mapping);if(selected>0){if(state.decoded->summary().at("dynamics_names").array().size()!=1)throw std::invalid_argument("One native mapping cannot be applied to multiple legacy dynamics names; preserve their identities and reclassify separately");const auto& system=state.registry->systems[state.systems.at(static_cast<std::size_t>(selected-1))];request.options.mappedSystem.identifier=system.name;const int preset=choice(state.preset);if(preset>0)request.options.mappedSystem.preset=system.presets.at(static_cast<std::size_t>(preset-1)).key;const auto parameters=json::parse(jsonObject(utf8Text(state.parameters)));for(const auto& item:parameters.object())request.options.mappedSystem.parameters[item.first]=item.second.number();}
        state.output=std::move(request);DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Cannot import legacy CSV",MB_OK|MB_ICONWARNING);}}

void createChildren(DialogState& state){control(state.window,L"STATIC",L"Import legacy tagged trajectory CSV",SS_LEFT,0,18,13,664,25,state.fonts->heading);
    control(state.window,L"STATIC",L"The source is read-only. Imported numbers remain authoritative; mapping a legacy name to executable native dynamics is always explicit.",SS_LEFT,0,18,44,664,36,state.fonts->caption);
    caption(state,L"Source CSV",18,91,664);state.source=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL,0,18,110,520,26,state.fonts->body);control(state.window,L"BUTTON",L"Browse…",BS_PUSHBUTTON,kBrowseId,548,108,134,30,state.fonts->body);
    state.summary=control(state.window,L"EDIT",L"Choose a tagged CSV and preview it.",ES_MULTILINE|ES_AUTOVSCROLL|ES_READONLY,0,18,151,520,86,state.fonts->mono);control(state.window,L"BUTTON",L"Preview / validate",BS_PUSHBUTTON,kPreviewId,548,151,134,32,state.fonts->body);
    state.strict=control(state.window,L"BUTTON",L"Reject unknown or unattached records",BS_AUTOCHECKBOX,0,548,194,134,43,state.fonts->caption);SendMessageW(state.strict,BM_SETCHECK,BST_CHECKED,0);
    caption(state,L"Import label",18,254,260);caption(state,L"Control input scale",290,254,130);caption(state,L"Coordinates",432,254,250);
    state.label=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL,0,18,273,260,26,state.fonts->body);state.inputScale=control(state.window,L"EDIT",L"1",ES_AUTOHSCROLL,0,290,273,130,26,state.fonts->mono);state.coordinates=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST,0,432,273,250,140,state.fonts->body);setComboItems(state.coordinates,{"unknown","dimensionless","dimensional"});
    caption(state,L"Dynamics mapping",18,314,320);caption(state,L"Native preset",350,314,332);state.mapping=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,kMappingId,18,333,320,240,state.fonts->body);SendMessageW(state.mapping,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"Preserve legacy identity (recommended)"));
    for(std::size_t index=0;index<state.registry->systems.size();++index)if(state.registry->systems[index].family=="ode"){state.systems.push_back(index);SendMessageW(state.mapping,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(dynsyslab::core::wide(state.registry->systems[index].label).c_str()));}SendMessageW(state.mapping,CB_SETCURSEL,0,0);
    state.preset=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,kPresetId,350,333,332,220,state.fonts->body);
    caption(state,L"Explicit native parameter overrides",18,374,664);state.parameters=control(state.window,L"EDIT",L"{}",ES_MULTILINE|ES_AUTOVSCROLL|ES_READONLY,0,18,393,520,67,state.fonts->mono);state.parameterButton=control(state.window,L"BUTTON",L"Edit parameters…",BS_PUSHBUTTON,kParametersId,548,393,134,32,state.fonts->body);
    control(state.window,L"STATIC",L"Canonical output: one Batch, grouped Systems, preserved Clusters and Trajectories, sample control history, interval-aligned regression inputs, and immutable source provenance.",SS_LEFT,0,18,476,664,42,state.fonts->caption);
    control(state.window,L"BUTTON",L"Import and archive",BS_DEFPUSHBUTTON,kImportId,426,531,134,32,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,570,531,112,32,state.fonts->body);updatePresets(state);}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);return 0;}if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kBrowseId){const auto path=chooseCsv(window);if(!path.empty()){SetWindowTextW(state->source,path.c_str());SetWindowTextW(state->label,path.stem().c_str());preview(*state);}return 0;}if(id==kPreviewId){preview(*state);return 0;}if(id==kMappingId&&HIWORD(wparam)==CBN_SELCHANGE){updatePresets(*state);return 0;}if(id==kPresetId&&HIWORD(wparam)==CBN_SELCHANGE){SetWindowTextW(state->parameters,L"{}");return 0;}if(id==kParametersId){editParameters(*state);return 0;}if(id==kImportId||id==IDOK){accept(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){static HBRUSH panel=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(panel);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}
}

std::optional<LegacyImportRequest> showLegacyImportDialog(HWND owner,const Registry& registry,const dsl::theme::Fonts& fonts){if(!registerClass())return std::nullopt;DialogState state;state.registry=&registry;state.fonts=&fonts;RECT area{};GetWindowRect(owner,&area);const int width=720,height=625;state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Legacy trajectory CSV import",WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
