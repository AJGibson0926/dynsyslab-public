#include "KromExtractionDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabKromExtractionDialog";
constexpr int kListId=1201,kConjugatesId=1202,kAllId=1203,kClearId=1204,kDoneId=1205,kCancelId=1206;
struct DialogState {HWND window=nullptr,list=nullptr;const dsl::theme::Fonts* fonts=nullptr;std::string current;
    std::vector<std::complex<double>> eigenvalues;std::vector<double> sourceIds;bool continuous=false;std::optional<std::string> output;};
HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,int x,int y,int width,int height,HFONT font){const bool edge=std::wstring(type)==L"LISTBOX";HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);dsl::theme::styleControl(result,font);return result;}
std::wstring eigenvalue(std::complex<double> value){std::wostringstream out;out<<std::setprecision(8)<<value.real();if(value.imag()>=0)out<<L'+';out<<value.imag()<<L'i';return out.str();}
std::vector<std::size_t> selected(const DialogState& state){const auto count=SendMessageW(state.list,LB_GETSELCOUNT,0,0);std::vector<int> indices(static_cast<std::size_t>(std::max<LRESULT>(0,count)));if(count>0)SendMessageW(state.list,LB_GETSELITEMS,count,reinterpret_cast<LPARAM>(indices.data()));std::vector<std::size_t> result;for(const int index:indices)result.push_back(static_cast<std::size_t>(index));return result;}
std::string encoded(const std::vector<std::size_t>& indices){std::string result="[";for(std::size_t index=0;index<indices.size();++index){if(index)result+=',';result+=std::to_string(indices[index]);}return result+"]";}
void selectConjugates(DialogState& state){const auto seeds=selected(state);for(const auto source:seeds){const auto desired=std::conj(state.eigenvalues[source]);double best=std::numeric_limits<double>::infinity();std::size_t bestIndex=source;
    for(std::size_t index=0;index<state.eigenvalues.size();++index){const double distance=std::abs(state.eigenvalues[index]-desired);if(distance<best){best=distance;bestIndex=index;}}
    if(best<=1e-7*(1.0+std::abs(desired)))SendMessageW(state.list,LB_SETSEL,TRUE,bestIndex);}}
void createChildren(DialogState& state){control(state.window,L"STATIC",L"Extract an independent sub-KROM",SS_LEFT,0,18,13,550,24,state.fonts->heading);
    control(state.window,L"STATIC",state.continuous?L"Spectrum shown as continuous generator λ":L"Spectrum shown as unitless discrete μ",SS_LEFT,0,18,43,550,20,state.fonts->caption);
    state.list=control(state.window,L"LISTBOX",L"",LBS_EXTENDEDSEL|LBS_NOINTEGRALHEIGHT|WS_VSCROLL|WS_HSCROLL,kListId,18,68,550,310,state.fonts->mono);
    for(std::size_t index=0;index<state.eigenvalues.size();++index){std::wstring row=L"["+std::to_wstring(index)+L"]  "+eigenvalue(state.eigenvalues[index]);if(index<state.sourceIds.size())row+=L"   source #"+std::to_wstring(static_cast<std::size_t>(state.sourceIds[index]));SendMessageW(state.list,LB_ADDSTRING,0,reinterpret_cast<LPARAM>(row.c_str()));}
    const auto settings=dsl::json::parse(jsonObject(state.current));for(const auto& item:settings.at("selected_modes").array()){const auto index=static_cast<std::size_t>(item.number());if(index<state.eigenvalues.size())SendMessageW(state.list,LB_SETSEL,TRUE,index);}
    control(state.window,L"BUTTON",L"Add conjugate partners",BS_PUSHBUTTON,kConjugatesId,18,391,150,29,state.fonts->body);control(state.window,L"BUTTON",L"Select all",BS_PUSHBUTTON,kAllId,178,391,90,29,state.fonts->body);control(state.window,L"BUTTON",L"Clear",BS_PUSHBUTTON,kClearId,278,391,75,29,state.fonts->body);
    control(state.window,L"STATIC",L"Extraction copies compatible C maps, evidence, and control-map blocks into a self-contained KROM; it never references its parent at runtime.",SS_LEFT,0,18,433,550,38,state.fonts->caption);
    control(state.window,L"BUTTON",L"Extract selected",BS_DEFPUSHBUTTON,kDoneId,330,482,126,30,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,466,482,102,30,state.fonts->body);}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid extraction settings",MB_OK|MB_ICONERROR);return -1;}return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kConjugatesId){selectConjugates(*state);return 0;}if(id==kAllId){SendMessageW(state->list,LB_SETSEL,TRUE,-1);return 0;}if(id==kClearId){SendMessageW(state->list,LB_SETSEL,FALSE,-1);return 0;}if(id==kDoneId||id==IDOK){const auto indices=selected(*state);if(indices.empty()){MessageBoxW(window,L"Select at least one eigenfunction.",L"No modes selected",MB_OK|MB_ICONWARNING);return 0;}state->output=mergeSettingsObject(state->current,{{"selected_modes",encoded(indices)}});DestroyWindow(window);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLORLISTBOX){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<std::string> showKromExtractionDialog(HWND owner,const std::filesystem::path& source,const std::string& current,const dsl::theme::Fonts& fonts){DialogState state;state.current=current;state.fonts=&fonts;try{const auto archive=dynsyslab::core::readNpz(source);if(archive.contains("continuous_eigenvalues")){state.eigenvalues=archive.at("continuous_eigenvalues").complexDoubles();state.continuous=true;}else if(archive.contains("discrete_eigenvalues"))state.eigenvalues=archive.at("discrete_eigenvalues").complexDoubles();else state.eigenvalues=archive.at("eigenvalues").complexDoubles();if(archive.contains("source_mode_indices"))state.sourceIds=archive.at("source_mode_indices").doubles();if(state.eigenvalues.empty())throw std::invalid_argument("KROM spectrum is empty");}catch(const std::exception& error){MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),L"Unable to read KROM",MB_OK|MB_ICONERROR);return std::nullopt;}if(!registerClass())return std::nullopt;RECT area{};GetWindowRect(owner,&area);const int width=610,height=575;state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"KROM eigenfunction selection",WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
