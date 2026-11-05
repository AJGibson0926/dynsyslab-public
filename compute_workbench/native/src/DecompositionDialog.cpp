#include "DecompositionDialog.hpp"

#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <utility>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabDecompositionDialog";
constexpr int kDoneId=1501,kCancelId=1502;
struct DialogState {
    HWND window=nullptr,dmdRank=nullptr,dmdTolerance=nullptr,dmdExact=nullptr;
    HWND havokObservable=nullptr,havokDelays=nullptr,havokRank=nullptr;
    HWND dictionary=nullptr,degree=nullptr,includeConstant=nullptr,regularization=nullptr,width=nullptr,centers=nullptr;
    HWND svdRank=nullptr,svdCenter=nullptr,pcaRank=nullptr,podRank=nullptr;
    HWND baggedRank=nullptr,baggedCenter=nullptr,bags=nullptr,seed=nullptr;
    HWND robustLambda=nullptr,robustTolerance=nullptr,robustIterations=nullptr;
    HWND hankelRank=nullptr,hankelSize=nullptr,hankelStride=nullptr,hankelCenter=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;DecompositionSettings current;
    std::vector<std::string> scalarObservables;std::optional<DecompositionSettings> output;
};
HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,int x,int y,int width,int height,HFONT font){
    const bool edge=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,
        type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;}
void label(DialogState& state,const wchar_t* text,int x,int y,int width){control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);}
HWND edit(DialogState& state,const std::string& value,int x,int y,int width){return control(state.window,L"EDIT",dynsyslab::core::wide(value).c_str(),ES_AUTOHSCROLL,0,x,y,width,25,state.fonts->mono);}
HWND check(DialogState& state,const wchar_t* text,bool selected,int x,int y,int width){HWND result=control(state.window,L"BUTTON",text,BS_AUTOCHECKBOX,0,x,y,width,24,state.fonts->body);SendMessageW(result,BM_SETCHECK,selected?BST_CHECKED:BST_UNCHECKED,0);return result;}
HWND combo(DialogState& state,const std::vector<std::string>& values,const std::string& selected,bool editable,int x,int y,int width){HWND result=control(state.window,L"COMBOBOX",L"",(editable?CBS_DROPDOWN:CBS_DROPDOWNLIST)|WS_VSCROLL,0,x,y,width,220,state.fonts->body);int selectedIndex=-1;for(std::size_t index=0;index<values.size();++index){const auto wide=dynsyslab::core::wide(values[index]);SendMessageW(result,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(wide.c_str()));if(values[index]==selected)selectedIndex=static_cast<int>(index);}if(selectedIndex>=0)SendMessageW(result,CB_SETCURSEL,selectedIndex,0);else if(editable)SetWindowTextW(result,dynsyslab::core::wide(selected).c_str());return result;}
std::string shown(const dsl::json::Value& value,const std::string& fallback){return value.isNull()?fallback:dsl::json::dump(value);}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string arrayLiteral(std::string value){const auto first=value.find_first_not_of(" \t\r\n");if(first==std::string::npos)return "[]";if(value[first]!='[')value="["+value+"]";const auto parsed=dsl::json::parse(value);if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))throw std::invalid_argument("RBF centers must be a JSON array.");return dsl::json::dump(parsed);}
std::string normalizeSystem(std::string value){std::transform(value.begin(),value.end(),value.begin(),[](unsigned char item){return static_cast<char>(std::tolower(item));});if(value.find("rayleigh")!=std::string::npos||value=="rpe")return "rpe";if(value.find("duffing")!=std::string::npos)return "duffing";if(value.find("marmottant")!=std::string::npos)return "marmottant";if(value.find("reddy")!=std::string::npos)return "reddy_szeri";return value;}
void loadObservables(DialogState& state,const std::filesystem::path& source){try{const auto archive=dynsyslab::core::readNpz(source);const auto metadata=dsl::json::parse(archive.metadataJson());dynsyslab::core::SystemConfiguration config;config.identifier=normalizeSystem(metadata.at("system").string());config.preset=metadata.at("system_preset").string(metadata.at("preset").string());const auto& resolved=metadata.at("resolved_system_parameters").object();const auto& parameters=resolved.empty()?metadata.at("system_parameters").object():resolved;for(const auto& item:parameters)config.parameters[item.first]=item.second.number();const auto catalog=dynsyslab::core::createSystem(config)->observableCatalog();for(const auto& name:catalog.names())if(catalog.at(name).outputDimension()==1)state.scalarObservables.push_back(name);}catch(...){ }if(state.scalarObservables.empty())state.scalarObservables={"0"};std::sort(state.scalarObservables.begin(),state.scalarObservables.end());}
void createChildren(DialogState& state){const auto dmd=dsl::json::parse(jsonObject(state.current.dmd)),havok=dsl::json::parse(jsonObject(state.current.havok)),edmd=dsl::json::parse(jsonObject(state.current.edmd)),svd=dsl::json::parse(jsonObject(state.current.svd)),pca=dsl::json::parse(jsonObject(state.current.pca)),pod=dsl::json::parse(jsonObject(state.current.pod)),bagged=dsl::json::parse(jsonObject(state.current.baggedPca)),robust=dsl::json::parse(jsonObject(state.current.robustPca)),hankel=dsl::json::parse(jsonObject(state.current.hankelPod));
    control(state.window,L"STATIC",L"DMD and HAVOK",SS_LEFT,0,18,12,580,23,state.fonts->heading);
    label(state,L"DMD rank (0 = full)",18,42,130);label(state,L"SVD tolerance",157,42,130);state.dmdRank=edit(state,shown(dmd.at("rank"),"0"),18,61,130);state.dmdTolerance=edit(state,shown(dmd.at("svd_tolerance"),"0"),157,61,130);state.dmdExact=check(state,L"Exact DMD modes",dmd.at("exact").boolean(true),304,61,145);
    label(state,L"HAVOK observable (name or state index)",18,99,260);label(state,L"Delays",288,99,130);label(state,L"Rank",428,99,130);state.havokObservable=combo(state,state.scalarObservables,shown(havok.at("observable"),"0"),true,18,118,260);state.havokDelays=edit(state,shown(havok.at("delays"),"100"),288,118,130);state.havokRank=edit(state,shown(havok.at("rank"),"15"),428,118,130);
    control(state.window,L"STATIC",L"EDMD dictionary",SS_LEFT,0,18,164,580,23,state.fonts->heading);
    label(state,L"Dictionary",18,194,130);label(state,L"Polynomial degree",157,194,130);label(state,L"Regularization",296,194,130);label(state,L"RBF width",435,194,123);state.dictionary=combo(state,{"polynomial","rbf"},edmd.at("dictionary").string("polynomial"),false,18,213,130);state.degree=edit(state,shown(edmd.at("polynomial_degree"),"2"),157,213,130);state.regularization=edit(state,shown(edmd.at("regularization"),"0"),296,213,130);state.width=edit(state,shown(edmd.at("width"),"1"),435,213,123);state.includeConstant=check(state,L"Include polynomial constant",edmd.at("include_constant").boolean(true),18,251,220);label(state,L"RBF centers [[x1,...], ...]",249,251,309);state.centers=edit(state,shown(edmd.at("centers"),"[]"),249,270,309);
    control(state.window,L"STATIC",L"PCA / POD / SVD family",SS_LEFT,0,18,314,580,23,state.fonts->heading);
    label(state,L"SVD rank",18,344,95);label(state,L"PCA rank",122,344,95);label(state,L"POD rank",226,344,95);state.svdRank=edit(state,shown(svd.at("rank"),"0"),18,363,95);state.pcaRank=edit(state,shown(pca.at("rank"),"0"),122,363,95);state.podRank=edit(state,shown(pod.at("rank"),"0"),226,363,95);state.svdCenter=check(state,L"Center SVD",svd.at("center").boolean(false),334,363,120);
    label(state,L"Bagged rank",18,407,105);label(state,L"Bag count",132,407,105);label(state,L"Bag RNG seed",246,407,105);state.baggedRank=edit(state,shown(bagged.at("rank"),"0"),18,426,105);state.bags=edit(state,shown(bagged.at("bags"),"32"),132,426,105);state.seed=edit(state,shown(bagged.at("seed"),"0"),246,426,105);state.baggedCenter=check(state,L"Center bagged PCA",bagged.at("bagged_center").boolean(bagged.at("center").boolean(true)),367,426,170);
    label(state,L"Robust lambda (0 = auto)",18,470,170);label(state,L"Robust tolerance",197,470,170);label(state,L"Robust max iterations",376,470,182);state.robustLambda=edit(state,shown(robust.at("lambda"),shown(robust.at("lambda_value"),"0")),18,489,170);state.robustTolerance=edit(state,shown(robust.at("tolerance"),"1e-7"),197,489,170);state.robustIterations=edit(state,shown(robust.at("max_iterations"),"1000"),376,489,182);
    label(state,L"Hankel rank",18,533,125);label(state,L"Matrix size",152,533,125);label(state,L"Delay stride",286,533,125);state.hankelRank=edit(state,shown(hankel.at("rank"),"0"),18,552,125);state.hankelSize=edit(state,shown(hankel.at("matrix_size"),"10"),152,552,125);state.hankelStride=edit(state,shown(hankel.at("delay_stride"),"1"),286,552,125);state.hankelCenter=check(state,L"Center Hankel-POD",hankel.at("center").boolean(false),425,552,173);
    control(state.window,L"STATIC",L"Every checked method is independently configured and archived; no decomposition forces a new solve.",SS_LEFT,0,18,595,580,28,state.fonts->caption);control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,367,634,118,30,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,495,634,103,30,state.fonts->body);}
void accept(DialogState& state){try{const auto number=[](HWND item){return dsl::json::dump(dsl::json::parse(utf8Text(item)));};DecompositionSettings result;
    result.dmd=mergeSettingsObject(state.current.dmd,{{"rank",number(state.dmdRank)},{"svd_tolerance",number(state.dmdTolerance)},{"exact",checked(state.dmdExact)?"true":"false"}});
    const std::string observable=utf8Text(state.havokObservable);std::string observableJson;try{observableJson=dsl::json::dump(dsl::json::parse(observable));}catch(...){observableJson=quoted(observable);}result.havok=mergeSettingsObject(state.current.havok,{{"observable",observableJson},{"delays",number(state.havokDelays)},{"rank",number(state.havokRank)}});
    result.edmd=mergeSettingsObject(state.current.edmd,{{"dictionary",quoted(utf8Text(state.dictionary))},{"polynomial_degree",number(state.degree)},{"include_constant",checked(state.includeConstant)?"true":"false"},{"regularization",number(state.regularization)},{"width",number(state.width)},{"centers",arrayLiteral(utf8Text(state.centers))}});
    result.svd=mergeSettingsObject(state.current.svd,{{"rank",number(state.svdRank)},{"center",checked(state.svdCenter)?"true":"false"}});
    result.pca=mergeSettingsObject(state.current.pca,{{"rank",number(state.pcaRank)}});result.pod=mergeSettingsObject(state.current.pod,{{"rank",number(state.podRank)}});
    result.baggedPca=mergeSettingsObject(state.current.baggedPca,{{"rank",number(state.baggedRank)},{"bagged_center",checked(state.baggedCenter)?"true":"false"},{"bags",number(state.bags)},{"seed",number(state.seed)}});
    result.robustPca=mergeSettingsObject(state.current.robustPca,{{"lambda",number(state.robustLambda)},{"tolerance",number(state.robustTolerance)},{"max_iterations",number(state.robustIterations)}});
    result.hankelPod=mergeSettingsObject(state.current.hankelPod,{{"rank",number(state.hankelRank)},{"matrix_size",number(state.hankelSize)},{"delay_stride",number(state.hankelStride)},{"center",checked(state.hankelCenter)?"true":"false"}});state.output=std::move(result);DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid decomposition setting",MB_OK|MB_ICONWARNING);}}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid decomposition settings",MB_OK|MB_ICONERROR);return -1;}return 0;}if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kDoneId||id==IDOK){accept(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<DecompositionSettings> showDecompositionDialog(HWND owner,const std::filesystem::path& source,const DecompositionSettings& current,const dsl::theme::Fonts& fonts){if(!registerClass())return std::nullopt;DialogState state;state.fonts=&fonts;state.current=current;loadObservables(state,source);RECT area{};GetWindowRect(owner,&area);const int width=640,height=728;state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Decomposition settings",WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
