#include "KoopmanPipelineDialog.hpp"

#include "SettingsMerge.hpp"
#include "HankelSettings.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabKoopmanPipelineDialog";
constexpr int kDoneId=801,kCancelId=802;

struct DialogState {
    HWND window=nullptr,hdObservable=nullptr,matrixSize=nullptr,samplingWindow=nullptr;
    HWND sampleStride=nullptr,svdTolerance=nullptr,rank=nullptr,energySort=nullptr;
    HWND snapshotStride=nullptr,alignment=nullptr,hdReconstructions=nullptr;
    HWND coordinates=nullptr,reconstructions=nullptr,degree=nullptr,threshold=nullptr;
    HWND iterations=nullptr,gamma=nullptr,normalize=nullptr,reuse=nullptr,materialize=nullptr,inputTerms=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;
    std::string hdmd,sindy;
    std::vector<std::string> scalarObservables,coordinateSuggestions;
    std::optional<KoopmanPipelineSettings> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,
             int id,int x,int y,int width,int height,HFONT font) {
    const bool bordered=std::wstring(type)==L"EDIT"||std::wstring(type)==L"COMBOBOX";
    HWND result=CreateWindowExW(bordered?WS_EX_CLIENTEDGE:0,type,text,
        WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;
}

HWND label(DialogState& state,const wchar_t* text,int x,int y,int width) {
    return control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);
}

HWND edit(DialogState& state,const std::string& value,int x,int y,int width) {
    return control(state.window,L"EDIT",dynsyslab::core::wide(value).c_str(),ES_AUTOHSCROLL,
                   0,x,y,width,25,state.fonts->mono);
}

HWND combo(DialogState& state,const std::vector<std::string>& values,const std::string& selected,
           bool editable,int x,int y,int width) {
    HWND result=control(state.window,L"COMBOBOX",L"",
        (editable?CBS_DROPDOWN:CBS_DROPDOWNLIST)|WS_VSCROLL,0,x,y,width,220,state.fonts->body);
    int selectedIndex=-1;
    for(std::size_t index=0;index<values.size();++index){SendMessageW(result,CB_ADDSTRING,0,
        reinterpret_cast<LPARAM>(dynsyslab::core::wide(values[index]).c_str()));
        if(values[index]==selected)selectedIndex=static_cast<int>(index);}
    if(selectedIndex>=0)SendMessageW(result,CB_SETCURSEL,selectedIndex,0);
    else if(editable)SetWindowTextW(result,dynsyslab::core::wide(selected).c_str());
    return result;
}

HWND check(DialogState& state,const wchar_t* text,bool selected,int x,int y,int width) {
    HWND result=control(state.window,L"BUTTON",text,BS_AUTOCHECKBOX,0,x,y,width,24,state.fonts->body);
    SendMessageW(result,BM_SETCHECK,selected?BST_CHECKED:BST_UNCHECKED,0);return result;
}

std::string text(HWND window){return utf8Text(window);}
std::string number(const dsl::json::Value& value,double fallback) {
    std::ostringstream out;out.precision(12);out<<value.number(fallback);return out.str();
}
std::string integer(const dsl::json::Value& value,std::size_t fallback) {
    return std::to_string(static_cast<std::size_t>(value.number(static_cast<double>(fallback))));
}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string trim(std::string value) {
    const auto first=value.find_first_not_of(" \t\r\n"),last=value.find_last_not_of(" \t\r\n");
    return first==std::string::npos?std::string{}:value.substr(first,last-first+1);
}
std::vector<std::string> split(std::string value) {
    std::vector<std::string> result;std::size_t start=0;
    while(start<=value.size()){const auto comma=value.find(',',start);auto item=trim(value.substr(start,comma-start));
        if(!item.empty())result.push_back(item);if(comma==std::string::npos)break;start=comma+1;}
    return result;
}
std::string array(const std::vector<std::string>& values) {
    std::string result="[";for(std::size_t index=0;index<values.size();++index){
        if(index)result+=',';result+=quoted(values[index]);}return result+"]";
}
std::string joined(const dsl::json::Value& value,const std::string& fallback) {
    std::string result;for(const auto& item:value.array()){
        if(!result.empty())result+=", ";result+=item.string();}return result.empty()?fallback:result;
}

std::string normalizeSystem(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char item){return static_cast<char>(std::tolower(item));});
    if(value.find("rayleigh")!=std::string::npos||value=="rpe")return "rpe";
    if(value.find("duffing")!=std::string::npos)return "duffing";
    if(value.find("marmottant")!=std::string::npos)return "marmottant";
    if(value.find("reddy")!=std::string::npos)return "reddy_szeri";
    return value;
}

void loadObservableSuggestions(DialogState& state,const std::filesystem::path& source) {
    try {
        const auto archive=dynsyslab::core::readNpz(source);const auto metadata=dsl::json::parse(archive.metadataJson());
        dynsyslab::core::SystemConfiguration configuration;
        configuration.identifier=normalizeSystem(metadata.at("system").string());
        configuration.preset=metadata.at("system_preset").string(metadata.at("preset").string());
        const auto& legacy=metadata.at("system_parameters");
        for(const auto* parameters : {&legacy.at("preset"), &legacy,
            &metadata.at("system_configuration").at("parameters"), &metadata.at("resolved_system_parameters")})
            for(const auto& item:parameters->object())
                if(std::holds_alternative<double>(item.second.storage()))
                    configuration.parameters[item.first]=item.second.number();
        const auto catalog=dynsyslab::core::createSystem(configuration)->observableCatalog();
        for(const auto& name:catalog.names()){
            state.coordinateSuggestions.push_back(name);
            if(catalog.at(name).outputDimension()==1)state.scalarObservables.push_back(name);
        }
        for(const auto& preset:catalog.coordinatePresets())state.coordinateSuggestions.push_back(preset.first);
    } catch(...) { /* A typed free-form field remains available for legacy sources. */ }
    if(state.scalarObservables.empty())state.scalarObservables={"state_0"};
    if(state.coordinateSuggestions.empty())state.coordinateSuggestions={"state"};
    std::sort(state.scalarObservables.begin(),state.scalarObservables.end());
    std::sort(state.coordinateSuggestions.begin(),state.coordinateSuggestions.end());
    state.scalarObservables.erase(std::unique(state.scalarObservables.begin(),state.scalarObservables.end()),state.scalarObservables.end());
    state.coordinateSuggestions.erase(std::unique(state.coordinateSuggestions.begin(),state.coordinateSuggestions.end()),state.coordinateSuggestions.end());
}

void createChildren(DialogState& state) {
    dsl::json::Value hd=dsl::json::parse(jsonObject(state.hdmd));
    dsl::json::Value sy=dsl::json::parse(jsonObject(state.sindy));
    control(state.window,L"STATIC",L"H-DMD spectrum",SS_LEFT,0,18,12,442,23,state.fonts->heading);
    label(state,L"Channels (names or JSON list)",18,42,205);label(state,L"Matrix size",232,42,105);label(state,L"Sampling window",346,42,114);
    state.hdObservable=combo(state,state.scalarObservables,hankel_ui::channelText(hd,"state_0"),true,18,61,205);
    const auto& win=static_cast<const dsl::json::Value&>(hd).at("window");
    const auto& geometry=win.isNull()?hd:win;
    state.matrixSize=edit(state,integer(geometry.at(win.isNull()?"matrix_size":"columns"),win.isNull()?2000:0),232,61,105);
    state.samplingWindow=edit(state,integer(geometry.at(win.isNull()?"sampling_window":"delays"),20000),346,61,114);
    label(state,L"Window stride",470,42,100);
    state.snapshotStride=edit(state,integer(geometry.at("snapshot_stride"),1),470,61,100);
    label(state,L"Fit additional outputs",375,98,195);
    state.hdReconstructions=edit(state,joined(hd.at("reconstruction_observables"),""),375,117,195);
    state.alignment=check(state,L"Causal endpoint",hd.at("window_alignment").string(win.isNull()?"start":"end")=="end",430,155,145);
    label(state,L"Sample stride",18,98,100);label(state,L"SVD tolerance",127,98,102);label(state,L"Rank cap (0 = full)",238,98,128);
    state.sampleStride=edit(state,integer(geometry.at("sample_stride"),1),18,117,100);
    state.svdTolerance=edit(state,number(hd.at("svd_tolerance"),1e-10),127,117,102);
    state.rank=edit(state,integer(hd.at("rank"),0),238,117,128);
    state.energySort=check(state,L"Sort retained modes by terminal-amplitude energy",
        hd.at("mode_sort").string("terminal-amplitude-energy")!="none",18,155,400);

    control(state.window,L"STATIC",L"SINDy eigenfunctions from supplied spectrum",SS_LEFT,0,18,196,442,23,state.fonts->heading);
    label(state,L"Coordinate observable(s), comma-separated",18,226,285);label(state,L"Reconstruct observables",312,226,148);
    const auto coordinateArray=joined(sy.at("coordinate_observables"),"");
    const auto coordinateValue=coordinateArray.empty()?sy.at("coordinate_observable").string("state"):coordinateArray;
    state.coordinates=combo(state,state.coordinateSuggestions,coordinateValue,true,18,245,285);
    state.reconstructions=edit(state,joined(sy.at("reconstruction_observables"),"state"),312,245,148);
    label(state,L"Polynomial degree",18,282,105);label(state,L"Sparse threshold",132,282,105);label(state,L"Max iterations",246,282,105);label(state,L"Generator derivative",360,282,100);
    state.degree=edit(state,integer(sy.at("polynomial_degree"),4),18,301,105);
    state.threshold=edit(state,number(sy.at("threshold"),0.05),132,301,105);
    state.iterations=edit(state,integer(sy.at("max_iterations"),20),246,301,105);
    state.gamma=combo(state,{"system_rhs","finite_difference"},sy.at("gamma_method").string("system_rhs"),false,360,301,100);
    state.normalize=check(state,L"Normalize coordinates",sy.at("normalize").boolean(true),18,340,143);
    state.reuse=check(state,L"Reuse conjugate fits",sy.at("reuse_conjugate_fits").boolean(true),170,340,151);
    state.materialize=check(state,L"Materialize training values (large)",sy.at("compute_training_values").boolean(false),330,340,240);
    state.inputTerms=check(state,L"Include known input in generator derivatives",
        sy.at("include_input_terms").boolean(true),18,372,530);
    control(state.window,L"STATIC",L"The standard policy retains the full H-DMD spectrum, converts discrete μ to continuous λ before SINDy, and stores coefficient-first eigenfunctions.",SS_LEFT,0,18,410,560,40,state.fonts->caption);
    control(state.window,L"BUTTON",L"Apply settings",BS_DEFPUSHBUTTON,kDoneId,340,466,118,30,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,467,466,101,30,state.fonts->body);
}

void accept(DialogState& state) {
    const auto coordinates=split(text(state.coordinates));
    std::vector<std::pair<std::string,std::string>> hd{
        {"matrix_size",text(state.matrixSize)},
        {"snapshot_stride",text(state.snapshotStride)},
        {"window_alignment",quoted(checked(state.alignment)?"end":"start")},
        {"reconstruction_observables",array(split(text(state.hdReconstructions)))},
        {"sampling_window",text(state.samplingWindow)},{"sample_stride",text(state.sampleStride)},
        {"svd_tolerance",text(state.svdTolerance)},{"rank",text(state.rank)},
        {"mode_sort",quoted(checked(state.energySort)?"terminal-amplitude-energy":"none")}};
    const auto channels=hankel_ui::channels(text(state.hdObservable));
    hd.push_back({"observables",dsl::json::dump(channels)});
    hd.push_back({"observable","null"});
    std::vector<std::pair<std::string,std::string>> sy{
        {"polynomial_degree",text(state.degree)},{"threshold",text(state.threshold)},
        {"max_iterations",text(state.iterations)},{"gamma_method",quoted(text(state.gamma))},
        {"normalize",checked(state.normalize)?"true":"false"},
        {"include_input_terms",checked(state.inputTerms)?"true":"false"},
        {"reuse_conjugate_fits",checked(state.reuse)?"true":"false"},
        {"compute_training_values",checked(state.materialize)?"true":"false"},
        {"archive_storage",quoted("coefficient-first")},
        {"reconstruction_observables",array(split(text(state.reconstructions)))}};
    if(coordinates.size()==1){sy.push_back({"coordinate_observable",quoted(coordinates.front())});
        sy.push_back({"coordinate_observables","[]"});}
    else {sy.push_back({"coordinate_observable",quoted("")});
        sy.push_back({"coordinate_observables",array(coordinates)});}
    state.output=KoopmanPipelineSettings{hankel_ui::preserveWindow(mergeSettingsObject(state.hdmd,hd),true),mergeSettingsObject(state.sindy,sy)};
    DestroyWindow(state.window);
}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){try{createChildren(*state);}catch(const std::exception& error){
        MessageBoxW(window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid advanced settings",MB_OK|MB_ICONERROR);return -1;}return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kDoneId||id==IDOK){try {accept(*state);} catch(const std::exception& e) {MessageBoxW(window,dynsyslab::core::wide(e.what()).c_str(),L"Invalid Hankel settings",MB_OK|MB_ICONWARNING);} return 0;}
        if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT||message==WM_CTLCOLORLISTBOX){
        static HBRUSH background=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);
        SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(background);}
    return DefWindowProcW(window,message,wparam,lparam);
}

bool registerDialogClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);
    type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<KoopmanPipelineSettings> showKoopmanPipelineDialog(
    HWND owner,const std::filesystem::path& source,const std::string& currentHdmd,
    const std::string& currentSindy,const dsl::theme::Fonts& fonts) {
    if(!registerDialogClass())return std::nullopt;
    DialogState state;state.fonts=&fonts;state.hdmd=currentHdmd;state.sindy=currentSindy;
    loadObservableSuggestions(state,source);RECT area{};GetWindowRect(owner,&area);const int width=605,height=558;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"Standard H-DMD → SINDy pipeline",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);
    ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};
    while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute
