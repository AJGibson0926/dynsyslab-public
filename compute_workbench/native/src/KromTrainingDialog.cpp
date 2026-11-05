#include "KromTrainingDialog.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TrainingDomain.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace dsl::compute {
namespace {
namespace json = dynsyslab::core::json;
constexpr wchar_t kClassName[]=L"DynSysLabKromTrainingDialog";
constexpr int kKromId=1601,kPreviewId=1602,kApplyId=1603,kCancelId=1604;

struct SourceData {
    dynsyslab::core::AxisAlignedStateDomain domain;
    json::Value metadata;
};

struct DialogState {
    HWND window=nullptr,krom=nullptr,current=nullptr,requiredMin=nullptr,requiredMax=nullptr;
    HWND points=nullptr,anchors=nullptr,margin=nullptr,outside=nullptr,fullStop=nullptr;
    HWND fullStep=nullptr,screenStop=nullptr,screenStep=nullptr,screenMethod=nullptr;
    HWND safetyMin=nullptr,safetyMax=nullptr,safetyAbs=nullptr,summary=nullptr;
    const dsl::theme::Fonts* fonts=nullptr;
    std::vector<KromDefinition> kroms;
    std::optional<SourceData> source;
    std::optional<KromTrainingSetup> output;
};

HWND control(HWND parent,const wchar_t* type,const wchar_t* text,DWORD style,int id,
             int x,int y,int width,int height,HFONT font) {
    const std::wstring kind(type);
    const bool edge=kind==L"EDIT"||kind==L"COMBOBOX"||kind==L"LISTBOX";
    HWND result=CreateWindowExW(edge?WS_EX_CLIENTEDGE:0,type,text,WS_CHILD|WS_VISIBLE|style,
        x,y,width,height,parent,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr),nullptr);
    dsl::theme::styleControl(result,font);return result;
}

void label(DialogState& state,const wchar_t* text,int x,int y,int width) {
    control(state.window,L"STATIC",text,SS_LEFT,0,x,y,width,18,state.fonts->caption);
}

std::string text(HWND window) { return utf8Text(window); }

std::string number(double value) {
    std::ostringstream out; out<<std::setprecision(12)<<value; return out.str();
}

json::Value vectorValue(const std::vector<double>& values) {
    json::Value::Array result;
    for(double value:values) result.emplace_back(value);
    return json::Value(std::move(result));
}

std::string vectorText(const std::vector<double>& values) {
    return json::dump(vectorValue(values));
}

std::vector<double> numericVector(const std::string& source,const char* field) {
    const auto value=json::parse(source.empty()?"[]":source);
    std::vector<double> result;
    for(const auto& item:value.array()) {
        const double number=item.number(std::numeric_limits<double>::quiet_NaN());
        if(!std::isfinite(number))throw std::invalid_argument(std::string(field)+" must contain finite numbers");
        result.push_back(number);
    }
    return result;
}

std::vector<std::size_t> pointCounts(const std::string& source,std::size_t dimension) {
    const auto values=numericVector(source,"Points per axis");
    if(values.size()!=dimension)throw std::invalid_argument("Points per axis must cover every state component");
    std::vector<std::size_t> result;
    for(double value:values){if(value<2||value!=std::floor(value))throw std::invalid_argument("Points per axis must be integers of at least two");result.push_back(static_cast<std::size_t>(value));}
    return result;
}

std::vector<double> anchors(const std::string& source,std::size_t dimension,std::size_t& rows) {
    const auto parsed=json::parse(source.empty()?"[]":source);std::vector<double> result;rows=0;
    for(const auto& row:parsed.array()){if(row.array().size()!=dimension)throw std::invalid_argument("Every anchor must contain one value per state component");
        for(const auto& item:row.array())result.push_back(item.number());++rows;}
    return result;
}

double positive(HWND field,const char* name) {
    const double value=std::stod(text(field));
    if(!std::isfinite(value)||value<=0)throw std::invalid_argument(std::string(name)+" must be positive");
    return value;
}

std::size_t selected(DialogState& state) {
    const auto index=SendMessageW(state.krom,CB_GETCURSEL,0,0);
    if(index==CB_ERR||static_cast<std::size_t>(index)>=state.kroms.size())
        throw std::invalid_argument("Select a compatible KROM");
    return static_cast<std::size_t>(index);
}

SourceData loadSource(const KromDefinition& definition) {
    const auto archive=dynsyslab::core::readNpz(std::filesystem::u8path(definition.source));
    if(!archive.contains("training_states"))
        throw std::invalid_argument("Selected KROM has no retained training states");
    const auto& states=archive.at("training_states");
    if(states.shape.size()!=2||!states.shape[0]||!states.shape[1])
        throw std::invalid_argument("Selected KROM has an invalid training-state matrix");
    const auto values=states.doubles();const std::size_t rows=states.shape[0],dimension=states.shape[1];
    SourceData result;result.domain.minimum.assign(dimension,std::numeric_limits<double>::infinity());
    result.domain.maximum.assign(dimension,-std::numeric_limits<double>::infinity());
    for(std::size_t row=0;row<rows;++row)for(std::size_t axis=0;axis<dimension;++axis){const double value=values[row*dimension+axis];
        result.domain.minimum[axis]=std::min(result.domain.minimum[axis],value);result.domain.maximum[axis]=std::max(result.domain.maximum[axis],value);}
    result.metadata=json::parse(archive.metadataJson());
    result.domain.coordinateSystem=result.metadata.at("integration_form").string(
        result.metadata.at("coordinate_system").string("dimensionless"));
    return result;
}

void displaySource(DialogState& state) {
    try {state.source=loadSource(state.kroms[selected(state)]);const auto& domain=state.source->domain;
        const std::string current="min "+vectorText(domain.minimum)+"   max "+vectorText(domain.maximum);
        SetWindowTextW(state.current,dynsyslab::core::wide(current).c_str());
        SetWindowTextW(state.requiredMin,dynsyslab::core::wide(vectorText(domain.minimum)).c_str());
        SetWindowTextW(state.requiredMax,dynsyslab::core::wide(vectorText(domain.maximum)).c_str());
        std::vector<double> defaults(domain.dimension(),5.0);SetWindowTextW(state.points,dynsyslab::core::wide(vectorText(defaults)).c_str());
        SetWindowTextW(state.summary,L"Edit the required envelope, then preview the exact IC count.");
    } catch(const std::exception& error){state.source.reset();SetWindowTextW(state.current,L"Unavailable");
        SetWindowTextW(state.summary,dynsyslab::core::wide(error.what()).c_str());}
}

dynsyslab::core::StateDomainExpansionPlan buildPlan(DialogState& state) {
    if(!state.source)throw std::invalid_argument("The selected KROM could not be loaded");
    const auto& current=state.source->domain;dynsyslab::core::AxisAlignedStateDomain required;
    required.minimum=numericVector(text(state.requiredMin),"Required minimum");required.maximum=numericVector(text(state.requiredMax),"Required maximum");required.coordinateSystem=current.coordinateSystem;
    dynsyslab::core::StateDomainExpansionOptions options;options.pointsPerAxis=pointCounts(text(state.points),current.dimension());
    options.marginFraction=std::stod(text(state.margin));options.outsideCurrentOnly=checked(state.outside);
    options.anchors=anchors(text(state.anchors),current.dimension(),options.anchorCount);
    auto plan=dynsyslab::core::planStateDomainExpansion(current,required,options);
    if(!plan.initialConditionCount)throw std::invalid_argument("The requested envelope adds no new grid points; enlarge it or add anchors");
    return plan;
}

void preview(DialogState& state) {
    try {const auto plan=buildPlan(state);std::ostringstream out;out<<plan.initialConditionCount
        <<" exact initial conditions; planned min "<<vectorText(plan.planned.minimum)
        <<"; max "<<vectorText(plan.planned.maximum);SetWindowTextW(state.summary,dynsyslab::core::wide(out.str()).c_str());}
    catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Invalid training-domain plan",MB_OK|MB_ICONWARNING);}
}

json::Value rowsValue(const std::vector<double>& values,std::size_t count,std::size_t dimension) {
    json::Value::Array rows;for(std::size_t row=0;row<count;++row){json::Value::Array columns;
        for(std::size_t axis=0;axis<dimension;++axis)columns.emplace_back(values[row*dimension+axis]);rows.emplace_back(std::move(columns));}return json::Value(std::move(rows));
}

void apply(DialogState& state) {
    try {const std::size_t sourceIndex=selected(state);const auto plan=buildPlan(state);const auto& metadata=state.source->metadata;
        const double fullStop=positive(state.fullStop,"Full integration stop");const double fullStep=positive(state.fullStep,"Full output step");
        const double screenStop=positive(state.screenStop,"Screening stop");const double screenStep=positive(state.screenStep,"Screening maximum step");
        if(fullStep>fullStop||screenStep>screenStop)throw std::invalid_argument("Integration steps must not exceed their time horizons");
        json::Value::Object currentDomain{{"minimum",vectorValue(plan.current.minimum)},{"maximum",vectorValue(plan.current.maximum)}};
        json::Value::Object requiredDomain{{"minimum",vectorValue(plan.required.minimum)},{"maximum",vectorValue(plan.required.maximum)}};
        json::Value::Array pointValues;for(auto value:plan.pointsPerAxis)pointValues.emplace_back(static_cast<double>(value));
        json::Value::Object screening{{"time_end",screenStop},{"output_step",screenStep},{"maximum_step",screenStep},
            {"method",text(state.screenMethod)},{"rtol",1.0e-9},{"atol",1.0e-12},
            {"state_minimum",vectorValue(numericVector(text(state.safetyMin),"Safety minimum"))},
            {"state_maximum",vectorValue(numericVector(text(state.safetyMax),"Safety maximum"))},
            {"maximum_absolute_state",vectorValue(numericVector(text(state.safetyAbs),"Safety absolute limit"))},{"require_all",true}};
        json::Value::Object preflight{{"enabled",true},{"kind","krom_training_domain"},
            {"source_krom",state.kroms[sourceIndex].source},{"coordinate_system",plan.planned.coordinateSystem},
            {"current_domain",json::Value(std::move(currentDomain))},{"required_domain",json::Value(std::move(requiredDomain))},
            {"points_per_axis",json::Value(std::move(pointValues))},{"anchors",json::parse(text(state.anchors).empty()?"[]":text(state.anchors))},
            {"margin_fraction",std::stod(text(state.margin))},{"outside_current_only",checked(state.outside)},
            {"screening",json::Value(std::move(screening))}};
        KromTrainingSetup setup;setup.selectedKrom=sourceIndex;setup.system=metadata.at("system").string(state.kroms[sourceIndex].system);
        setup.preset=metadata.at("system_preset").string(metadata.at("preset").string());const auto& parameters=metadata.at("resolved_system_parameters").object().empty()?metadata.at("system_parameters"):metadata.at("resolved_system_parameters");
        setup.parametersJson=json::dump(parameters);setup.coordinateSystem=plan.planned.coordinateSystem;
        setup.initialConditionsJson=json::dump(rowsValue(plan.initialConditions,plan.initialConditionCount,plan.planned.dimension()));setup.trajectoryCount=plan.initialConditionCount;
        setup.finalTime=number(fullStop);setup.outputStep=number(fullStep);setup.preflightJson=json::dump(json::Value(std::move(preflight)));
        setup.label=setup.system+"-krom-training-domain";state.output=std::move(setup);DestroyWindow(state.window);
    }catch(const std::exception& error){MessageBoxW(state.window,dynsyslab::core::wide(error.what()).c_str(),L"Cannot apply training-domain plan",MB_OK|MB_ICONWARNING);}
}

void createChildren(DialogState& state) {
    control(state.window,L"STATIC",L"KROM training-domain expansion",SS_LEFT,0,18,12,690,26,state.fonts->heading);
    label(state,L"Source KROM",18,43,690);state.krom=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST|WS_VSCROLL,kKromId,18,62,690,260,state.fonts->body);
    for(const auto& item:state.kroms)SendMessageW(state.krom,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(dynsyslab::core::wide(item.label).c_str()));
    SendMessageW(state.krom,CB_SETCURSEL,0,0);label(state,L"Current retained training envelope",18,98,690);
    state.current=control(state.window,L"EDIT",L"",ES_AUTOHSCROLL|ES_READONLY,0,18,117,690,25,state.fonts->mono);
    label(state,L"Required minimum (JSON)",18,154,330);label(state,L"Required maximum (JSON)",378,154,330);
    state.requiredMin=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,18,173,330,25,state.fonts->mono);state.requiredMax=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,378,173,330,25,state.fonts->mono);
    label(state,L"Points per axis",18,210,205);label(state,L"Anchors (JSON rows)",235,210,300);label(state,L"Margin fraction",547,210,161);
    state.points=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,18,229,205,25,state.fonts->mono);state.anchors=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,235,229,300,25,state.fonts->mono);state.margin=control(state.window,L"EDIT",L"0.1",ES_AUTOHSCROLL,0,547,229,161,25,state.fonts->mono);
    state.outside=control(state.window,L"BUTTON",L"Generate only outside the current domain (anchors are always retained)",BS_AUTOCHECKBOX,0,18,265,600,24,state.fonts->body);SendMessageW(state.outside,BM_SETCHECK,BST_CHECKED,0);
    label(state,L"Full stop",18,301,125);label(state,L"Output step",155,301,125);label(state,L"Screen stop",292,301,125);label(state,L"Screen max step",429,301,125);label(state,L"Screen method",566,301,142);
    state.fullStop=control(state.window,L"EDIT",L"10",ES_AUTOHSCROLL,0,18,320,125,25,state.fonts->mono);state.fullStep=control(state.window,L"EDIT",L"0.01",ES_AUTOHSCROLL,0,155,320,125,25,state.fonts->mono);state.screenStop=control(state.window,L"EDIT",L"0.25",ES_AUTOHSCROLL,0,292,320,125,25,state.fonts->mono);state.screenStep=control(state.window,L"EDIT",L"0.002",ES_AUTOHSCROLL,0,429,320,125,25,state.fonts->mono);
    state.screenMethod=control(state.window,L"COMBOBOX",L"",CBS_DROPDOWNLIST,0,566,320,142,100,state.fonts->body);SendMessageW(state.screenMethod,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"RK45"));SendMessageW(state.screenMethod,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(L"BDF"));SendMessageW(state.screenMethod,CB_SETCURSEL,0,0);
    label(state,L"Safety state minimum (optional)",18,357,215);label(state,L"Safety state maximum (optional)",247,357,215);label(state,L"Maximum absolute state (optional)",476,357,232);
    state.safetyMin=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,18,376,215,25,state.fonts->mono);state.safetyMax=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,247,376,215,25,state.fonts->mono);state.safetyAbs=control(state.window,L"EDIT",L"[]",ES_AUTOHSCROLL,0,476,376,232,25,state.fonts->mono);
    state.summary=control(state.window,L"STATIC",L"",SS_LEFT,0,18,419,690,42,state.fonts->caption);
    control(state.window,L"BUTTON",L"Preview exact plan",BS_PUSHBUTTON,kPreviewId,18,472,145,31,state.fonts->body);
    control(state.window,L"BUTTON",L"Apply to Solver",BS_DEFPUSHBUTTON,kApplyId,474,472,112,31,state.fonts->body);control(state.window,L"BUTTON",L"Cancel",BS_PUSHBUTTON,kCancelId,596,472,112,31,state.fonts->body);
}

LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* state=reinterpret_cast<DialogState*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){createChildren(*state);displaySource(*state);return 0;}
    if(message==WM_COMMAND){const int id=LOWORD(wparam);if(id==kKromId&&HIWORD(wparam)==CBN_SELCHANGE){displaySource(*state);return 0;}if(id==kPreviewId){preview(*state);return 0;}if(id==kApplyId||id==IDOK){apply(*state);return 0;}if(id==kCancelId||id==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}if(message==WM_CTLCOLORSTATIC||message==WM_CTLCOLOREDIT){static HBRUSH brush=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(brush);}return DefWindowProcW(window,message,wparam,lparam);
}

bool registerClass(){WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;}
}

std::optional<KromTrainingSetup> showKromTrainingDialog(
    HWND owner,const std::vector<KromDefinition>& compatibleKroms,
    std::size_t initialSelection,const dsl::theme::Fonts& fonts) {
    if(compatibleKroms.empty()){MessageBoxW(owner,L"No compatible KROM with retained training states is registered for this system.",L"KROM training domain",MB_OK|MB_ICONINFORMATION);return std::nullopt;}
    if(!registerClass())return std::nullopt;DialogState state;state.kroms=compatibleKroms;state.fonts=&fonts;RECT area{};GetWindowRect(owner,&area);const int width=750,height=565;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"KROM training-domain planner",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window)return std::nullopt;if(initialSelection<state.kroms.size()){SendMessageW(state.krom,CB_SETCURSEL,initialSelection,0);displaySource(state);}
    dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};
    while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute
