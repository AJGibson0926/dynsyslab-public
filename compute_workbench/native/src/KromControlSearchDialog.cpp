#include "KromControlSearchDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <commdlg.h>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[]=L"DynSysLabKromControlSearchDialog";
constexpr int kRunId=1601,kPlanId=1602,kQueueId=1603,kCancelId=1604,kSaveId=1608;
constexpr int kParentId=1605,kBrowseParentId=1606,kAddSourcesId=1607;
struct State{HWND window=nullptr,jobLabel=nullptr,parent=nullptr,sources=nullptr,initial=nullptr,target=nullptr,
    finalState=nullptr,start=nullptr,activeStart=nullptr,activeStop=nullptr,stop=nullptr,
    outputStep=nullptr,solver=nullptr,clock=nullptr,rtol=nullptr,atol=nullptr,maxStep=nullptr,
    Q=nullptr,R=nullptr,maxAttempts=nullptr,maxModes=nullptr,virtualCount=nullptr,seedModes=nullptr,
    constantMethod=nullptr,bilinearMethod=nullptr,reconstruction=nullptr,ridge=nullptr,stride=nullptr,
    grid=nullptr,releaseTolerance=nullptr,finalTolerance=nullptr,settlingStart=nullptr,timeout=nullptr,
    conjugate=nullptr,stopOnSuccess=nullptr;const Registry* registry=nullptr;
    std::filesystem::path selected;
    const dsl::theme::Fonts* fonts=nullptr;std::optional<KromControlSearchRequest> output;};
HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
    int x,int y,int width,int height,HFONT font){const bool edge=std::wstring(type)==L"EDIT"
        ||std::wstring(type)==L"COMBOBOX";HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,
        type,text,WS_CHILD|WS_VISIBLE|style,x,y,width,height,parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;}
void label(State& state,const wchar_t* text,int x,int y,int width){control(state.window,L"STATIC",
    text,SS_LEFT,0,x,y,width,18,state.fonts->caption);}
HWND edit(State& state,const wchar_t* text,int x,int y,int width){return control(state.window,L"EDIT",
    text,ES_AUTOHSCROLL,0,x,y,width,25,state.fonts->mono);}
HWND combo(State& state,const std::vector<std::string>& values,int selected,int x,int y,int width,
    int id=0,DWORD style=CBS_DROPDOWNLIST){HWND result=control(state.window,L"COMBOBOX",L"",
        style|WS_VSCROLL,id,x,y,width,130,state.fonts->body);setComboItems(result,values,selected);return result;}
void field(State& state,const wchar_t* name,const wchar_t* value,int x,int y,int width,HWND& output){
    label(state,name,x,y,width);output=edit(state,value,x,y+19,width);}
std::string quoted(const std::string& value){return "\""+jsonEscape(value)+"\"";}
std::string literal(HWND item,const char* fallback){const auto text=utf8Text(item);
    return dsl::json::dump(dsl::json::parse(text.empty()?fallback:text));}
std::string scalarOrList(HWND item,const char* fallback){std::string text=utf8Text(item);
    if(text.empty())text=fallback;const auto begin=text.find_first_not_of(" \t\r\n");
    if(begin!=std::string::npos&&text[begin]!='['&&text.find(',')!=std::string::npos)
        text="["+text+"]";return dsl::json::dump(dsl::json::parse(text));}
std::string list(HWND item){std::string text=utf8Text(item);const auto begin=text.find_first_not_of(" \t\r\n");
    if(begin==std::string::npos)return "[]";if(text[begin]!='[')text="["+text+"]";
    const auto parsed=dsl::json::parse(text);if(!std::holds_alternative<dsl::json::Value::Array>(parsed.storage()))
        throw std::invalid_argument("State fields must be arrays or comma-separated values");return dsl::json::dump(parsed);}
std::vector<std::filesystem::path> chooseFiles(HWND owner,const wchar_t* title,bool multiple){
    std::vector<wchar_t> buffer(65536);const wchar_t filter[]=L"DynSysLab archives (*.npz)\0*.npz\0All files (*.*)\0*.*\0\0";
    OPENFILENAMEW request{};request.lStructSize=sizeof(request);request.hwndOwner=owner;
    request.lpstrFile=buffer.data();request.nMaxFile=static_cast<DWORD>(buffer.size());request.lpstrFilter=filter;
    request.nFilterIndex=1;request.lpstrTitle=title;request.Flags=OFN_EXPLORER|OFN_NOCHANGEDIR|OFN_FILEMUSTEXIST
        |(multiple?OFN_ALLOWMULTISELECT:0);if(!GetOpenFileNameW(&request))return {};
    std::vector<std::wstring> parts;for(const wchar_t* cursor=buffer.data();*cursor;cursor+=lstrlenW(cursor)+1)
        parts.emplace_back(cursor);std::vector<std::filesystem::path> result;
    if(parts.size()==1)result.emplace_back(parts[0]);else for(std::size_t i=1;i<parts.size();++i)
        result.push_back(std::filesystem::path(parts[0])/parts[i]);return result;}
void updateReconstruction(State& state){const auto preferred=utf8Text(state.reconstruction);
    std::vector<std::string> names;try{const auto model=dynsyslab::core::loadKoopmanModel(
        std::filesystem::u8path(utf8Text(state.parent)));for(const auto& item:model.reconstructions())
        names.push_back(item.first);}catch(...){ }if(names.empty())names.push_back("state");
    int selected=0;for(std::size_t i=0;i<names.size();++i)if(names[i]==preferred)selected=static_cast<int>(i);
    setComboItems(state.reconstruction,names,selected);}
void updateDefaults(State& state){try{const auto model=dynsyslab::core::loadKoopmanModel(
        std::filesystem::u8path(utf8Text(state.parent)));const auto system=
        dynsyslab::core::createSystem(model.systemConfiguration());const auto values=system->controlledInitialState();
        std::ostringstream out;for(std::size_t i=0;i<values.size();++i){if(i)out<<',';out<<values[i];}
        const auto shown=dynsyslab::core::wide(out.str());SetWindowTextW(state.initial,shown.c_str());
        SetWindowTextW(state.target,shown.c_str());SetWindowTextW(state.finalState,shown.c_str());}catch(...){ }}
void browseParent(State& state){const auto paths=chooseFiles(state.window,L"Select parent KROM",false);
    if(paths.empty())return;SetWindowTextW(state.parent,std::filesystem::absolute(paths.front()).wstring().c_str());
    updateReconstruction(state);updateDefaults(state);}
void addSources(State& state){auto text=windowText(state.sources);for(const auto& path:
    chooseFiles(state.window,L"Add controlled trajectory archives",true)){if(!text.empty()&&text.back()!=L'\n')
        text+=L"\r\n";text+=std::filesystem::absolute(path).wstring();}SetWindowTextW(state.sources,text.c_str());}
std::string sourceArray(HWND source){std::istringstream stream(utf8Text(source));std::string line,result="[";
    bool first=true;while(std::getline(stream,line)){while(!line.empty()&&(line.back()=='\r'
        ||std::isspace(static_cast<unsigned char>(line.back()))))line.pop_back();const auto begin=line.find_first_not_of(" \t");
        if(begin==std::string::npos)continue;if(!first)result+=',';result+=quoted(line.substr(begin));first=false;}
    return result+"]";}
std::filesystem::path initialParent(const std::filesystem::path& selected){if(selected.empty())return {};
    try{(void)dynsyslab::core::loadKoopmanModel(selected);return selected;}catch(...){return {};}}
std::filesystem::path initialSource(const std::filesystem::path& selected){if(selected.empty())return {};
    try{const auto record=dynsyslab::core::readNpzIndex(selected);
        if(record.contains("control_values"))return selected;}catch(...){ }return {};}
std::wstring defaults(const std::filesystem::path& parent,bool initial){try{const auto model=
    dynsyslab::core::loadKoopmanModel(parent);const auto system=dynsyslab::core::createSystem(model.systemConfiguration());
    const auto values=initial?system->controlledInitialState():system->initialState();std::ostringstream out;
    for(std::size_t i=0;i<values.size();++i){if(i)out<<',';out<<values[i];}return dynsyslab::core::wide(out.str());
    }catch(...){return L"0,0";}}

void createChildren(State& state,const std::filesystem::path& selected){const auto parent=initialParent(selected);
    control(state.window,L"STATIC",L"Bounded KROM control-candidate campaign",SS_LEFT,0,18,12,430,25,state.fonts->heading);
    label(state,L"Job label",470,16,86);state.jobLabel=edit(state,L"native-krom-control-search",556,10,266);
    label(state,L"Parent KROM",18,45,680);state.parent=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWN|WS_VSCROLL,
        kParentId,18,64,680,220,state.fonts->body);int choice=-1,index=0;for(const auto& item:state.registry->kroms){
        const auto path=std::filesystem::absolute(std::filesystem::u8path(item.source));SendMessageW(state.parent,CB_ADDSTRING,0,
            reinterpret_cast<LPARAM>(path.wstring().c_str()));if(path==std::filesystem::absolute(parent))choice=index;++index;}
    if(choice>=0)SendMessageW(state.parent,CB_SETCURSEL,choice,0);else if(!parent.empty())
        SetWindowTextW(state.parent,std::filesystem::absolute(parent).wstring().c_str());
    control(state.window,L"BUTTON",L"Browse…",BS_PUSHBUTTON,kBrowseParentId,708,64,114,30,state.fonts->body);
    label(state,L"Controlled regression ensemble (one archive per line)",18,101,680);
    const auto source=initialSource(selected);state.sources=control(state.window,L"EDIT",
        source.empty()?L"":std::filesystem::absolute(source).wstring().c_str(),ES_MULTILINE|ES_AUTOVSCROLL|WS_VSCROLL,
        0,18,120,680,48,state.fonts->mono);control(state.window,L"BUTTON",L"Add files…",BS_PUSHBUTTON,
        kAddSourcesId,708,120,114,31,state.fonts->body);
    const int x0=18,x1=290,x2=562,width=260;const auto initialValue=defaults(parent,true);
    field(state,L"Initial state",initialValue.c_str(),x0,184,width,state.initial);
    field(state,L"Controlled target state",initialValue.c_str(),x1,184,width,state.target);
    field(state,L"Final reference after release",initialValue.c_str(),x2,184,width,state.finalState);
    field(state,L"Start time",L"0",x0,239,width,state.start);field(state,L"Active start",L"1",x1,239,width,state.activeStart);
    field(state,L"Active stop / release",L"10",x2,239,width,state.activeStop);
    field(state,L"Final time",L"20",x0,294,width,state.stop);field(state,L"Output step",L"0.001",x1,294,width,state.outputStep);
    label(state,L"Integrator",x2,294,width);state.solver=combo(state,{"BDF","RK45","RK4"},0,x2,313,width);
    field(state,L"Relative tolerance",L"1e-8",x0,349,width,state.rtol);field(state,L"Absolute tolerance",L"1e-10",x1,349,width,state.atol);
    field(state,L"Maximum internal step",L"0.01",x2,349,width,state.maxStep);
    label(state,L"Controller evaluation",x0,404,width);state.clock=combo(state,{"continuous_rhs","sampled"},0,x0,423,width);
    field(state,L"KLQR Q (scalar/diagonal/matrix)",L"12500",x1,404,width,state.Q);
    field(state,L"KLQR R (scalar/matrix)",L"0.01",x2,404,width,state.R);
    field(state,L"Maximum attempts (hard limit 25)",L"25",x0,459,width,state.maxAttempts);
    field(state,L"Maximum candidate modes",L"4",x1,459,width,state.maxModes);
    field(state,L"Virtual candidate pool",L"5000",x2,459,width,state.virtualCount);
    field(state,L"Seed modes",L"12",x0,514,width,state.seedModes);
    label(state,L"Constant B fit",x1,514,width);state.constantMethod=combo(state,
        {"multistep-output","controlled-residual","exact-map"},0,x1,533,width);
    label(state,L"Bilinear B fit",x2,514,width);state.bilinearMethod=combo(state,
        {"controlled-residual","exact-map"},0,x2,533,width);
    label(state,L"Observable reconstruction C",x0,569,width);state.reconstruction=combo(state,{"state"},0,x0,588,width);
    field(state,L"Regression ridge",L"1e-12",x1,569,width,state.ridge);field(state,L"Exact-map sample stride",L"1",x2,569,width,state.stride);
    field(state,L"Error-map grid size",L"81",x0,624,width,state.grid);
    field(state,L"Release error tolerance",L"0.02",x1,624,width,state.releaseTolerance);
    field(state,L"Final-state error tolerance",L"0.1",x2,624,width,state.finalTolerance);
    field(state,L"Settling measurement starts",L"2",x0,679,width,state.settlingStart);
    field(state,L"Per-case timeout (seconds)",L"300",x1,679,width,state.timeout);
    state.conjugate=control(state.window,L"BUTTON",L"Conjugate-complete bilinear features",BS_AUTOCHECKBOX,0,x2,680,width,23,state.fonts->body);
    state.stopOnSuccess=control(state.window,L"BUTTON",L"Stop after first successful bilinear case",BS_AUTOCHECKBOX,0,x2,710,width,23,state.fonts->body);
    SendMessageW(state.stopOnSuccess,BM_SETCHECK,BST_CHECKED,0);updateReconstruction(state);
    control(state.window,L"BUTTON",L"Run campaign",BS_DEFPUSHBUTTON,kRunId,246,756,120,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Plan only",BS_PUSHBUTTON,kPlanId,376,756,100,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Queue",BS_PUSHBUTTON,kQueueId,486,756,100,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Save job…",BS_PUSHBUTTON,kSaveId,596,756,100,32,state.fonts->body);
    control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,706,756,116,32,state.fonts->body);}

void accept(State& state,KromControlSearchRequest::Action action){try{const auto grid=literal(state.grid,"81");
    std::ostringstream out;out<<"{\"schema\":\"dynsyslab.compute-job/v1\",\"kind\":\"krom_control_search\"," 
        "\"label\":"<<quoted(utf8Text(state.jobLabel))<<",\"parent_krom\":"<<quoted(utf8Text(state.parent))
        <<",\"controlled_sources\":"<<sourceArray(state.sources)<<",\"initial_state\":"<<list(state.initial)
        <<",\"target_state\":"<<list(state.target)<<",\"final_reference_state\":"<<list(state.finalState)
        <<",\"time\":{\"start\":"<<literal(state.start,"0")<<",\"active_window\":["<<literal(state.activeStart,"1")
        <<','<<literal(state.activeStop,"10")<<"],\"stop\":"<<literal(state.stop,"20")<<",\"output_step\":"
        <<literal(state.outputStep,"0.001")<<",\"coordinate_system\":\"automatic\"},\"solver\":{\"method\":"
        <<quoted(selectedText(state.solver))<<",\"rtol\":"<<literal(state.rtol,"1e-8")<<",\"atol\":"
        <<literal(state.atol,"1e-10")<<",\"max_step\":"<<literal(state.maxStep,"0.01")
        <<",\"control_evaluation\":"<<quoted(selectedText(state.clock))<<"},\"controller\":{\"Q\":"
        <<scalarOrList(state.Q,"12500")<<",\"R\":"<<scalarOrList(state.R,"0.01")<<"},\"regression\":{\"constant_method\":"
        <<quoted(selectedText(state.constantMethod))<<",\"bilinear_method\":"<<quoted(selectedText(state.bilinearMethod))
        <<",\"reconstruction\":"<<quoted(selectedText(state.reconstruction))<<",\"ridge\":"<<literal(state.ridge,"1e-12")
        <<",\"conjugate_complete\":"<<(checked(state.conjugate)?"true":"false")<<",\"map_sample_stride\":"
        <<literal(state.stride,"1")<<",\"grid_columns\":"<<grid<<",\"grid_rows\":"<<grid
        <<"},\"screening\":{\"max_attempts\":"<<literal(state.maxAttempts,"25")<<",\"maximum_mode_count\":"
        <<literal(state.maxModes,"4")<<",\"maximum_virtual_candidates\":"<<literal(state.virtualCount,"5000")
        <<",\"retained_candidates\":"<<literal(state.virtualCount,"5000")<<",\"seed_mode_count\":"
        <<literal(state.seedModes,"12")<<",\"stop_on_bilinear_success\":"<<(checked(state.stopOnSuccess)?"true":"false")
        <<"},\"success\":{\"release_error_tolerance\":"<<literal(state.releaseTolerance,"0.02")
        <<",\"final_error_tolerance\":"<<literal(state.finalTolerance,"0.1")<<",\"settling_start\":"
        <<literal(state.settlingStart,"2")<<"},\"case_timeout_seconds\":"<<literal(state.timeout,"300")
        <<",\"archive\":{\"enabled\":true}}";const auto value=dsl::json::dump(dsl::json::parse(out.str()));
    state.output=KromControlSearchRequest{value,action};DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),
        L"Invalid KROM control-search setting",MB_OK|MB_ICONWARNING);}}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam){auto* state=
    reinterpret_cast<State*>(GetWindowLongPtrW(window,GWLP_USERDATA));if(message==WM_NCCREATE){state=static_cast<State*>(
        reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,
        reinterpret_cast<LONG_PTR>(state));}if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state,state->selected);return 0;}if(message==WM_COMMAND){const int id=LOWORD(wparam);
        if(id==kParentId&&HIWORD(wparam)==CBN_SELCHANGE){updateReconstruction(*state);updateDefaults(*state);return 0;}
        if(id==kBrowseParentId){browseParent(*state);return 0;}if(id==kAddSourcesId){addSources(*state);return 0;}
        if(id==kRunId){accept(*state,KromControlSearchRequest::Action::run);return 0;}
        if(id==kPlanId){accept(*state,KromControlSearchRequest::Action::plan);return 0;}
        if(id==kQueueId){accept(*state,KromControlSearchRequest::Action::queue);return 0;}
        if(id==kSaveId){accept(*state,KromControlSearchRequest::Action::save);return 0;}
        if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){
        static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);
        SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);}
bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);
    type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;
    return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}

}  // namespace

std::optional<KromControlSearchRequest> showKromControlSearchDialog(HWND owner,
    const std::filesystem::path& selected,const Registry& registry,const dsl::theme::Fonts& fonts){if(!registerClass())return std::nullopt;
    State state;state.registry=&registry;state.fonts=&fonts;state.selected=selected;
    RECT area{};GetWindowRect(owner,&area);const int width=860,height=840;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"KROM control candidate search",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,
        GetModuleHandleW(nullptr),&state);if(!state.window)return std::nullopt;
    dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);
    MSG event{};while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){
        TranslateMessage(&event);DispatchMessageW(&event);}}EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;}

}  // namespace dsl::compute
