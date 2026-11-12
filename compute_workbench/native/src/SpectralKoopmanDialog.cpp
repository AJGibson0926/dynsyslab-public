#include "SpectralKoopmanDialog.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Text.hpp>
#include <cmath>
#include <stdexcept>
namespace dsl::compute {
namespace {
namespace j=dynsyslab::core::json;
constexpr wchar_t windowClass[]=L"DynSysLabSpectralKoopman";
const std::vector<std::string> methods{"resdmd","edmd","mpedmd","hermitian_dmd","tlsdmd","multdmd","stochastic_resdmd","specrkhs"};
const std::vector<std::string> dictionaries{"linear","polynomial","gaussian","partition"};
struct State {
    HWND window{},method{},dictionary{},rank{},degree{},centers{},bandwidth{},normalize{},constant{},invariant{},hermitian{},advanced{};
    const theme::Fonts* fonts{};std::string current;std::optional<std::string> output;
};
void place(State& s,HWND item,int x,int y,int width,int height) {
    MoveWindow(item,scalePixels(s.window,x),scalePixels(s.window,y),scalePixels(s.window,width),scalePixels(s.window,height),TRUE);
}
void label(State& s,const wchar_t* caption,int x,int y,int width) {
    place(s,makeLabel(s.window,caption,s.fonts->caption),x,y,width,20);
}
HWND field(State& s,const wchar_t* caption,std::string value,int x,int y,int width) {
    label(s,caption,x,y,width);auto item=makeEdit(s.window,0,dynsyslab::core::wide(value).c_str(),s.fonts->mono);
    place(s,item,x,y+22,width,28);return item;
}
HWND combo(State& s,const wchar_t* caption,const std::vector<std::string>& values,std::string selected,int x,int y,int width) {
    label(s,caption,x,y,width);
    auto item=CreateWindowExW(WS_EX_CLIENTEDGE,L"COMBOBOX",L"",WS_CHILD|WS_VISIBLE|WS_TABSTOP|CBS_DROPDOWNLIST|WS_VSCROLL,
        0,0,0,0,s.window,nullptr,GetModuleHandleW(nullptr),nullptr);
    theme::styleControl(item,s.fonts->body);place(s,item,x,y+22,width,230);int index=0;
    for(std::size_t i=0;i<values.size();++i){auto text=dynsyslab::core::wide(values[i]);SendMessageW(item,CB_ADDSTRING,0,reinterpret_cast<LPARAM>(text.c_str()));if(values[i]==selected)index=int(i);}
    SendMessageW(item,CB_SETCURSEL,index,0);return item;
}
void create(State& s) {
    const auto settings=j::parse(s.current);const auto& d=settings.at("dictionary");
    s.method=combo(s,L"Estimator",methods,settings.at("method").string("resdmd"),18,16,310);
    s.dictionary=combo(s,L"Dictionary",dictionaries,d.at("kind").string("linear"),346,16,270);
    auto shown=[](const j::Value& v,const char* fallback){return v.isNull()?std::string(fallback):j::dump(v);};
    s.rank=field(s,L"Rank (0 = numerical rank)",shown(settings.at("rank"),"0"),18,82,192);
    s.degree=field(s,L"Polynomial degree",shown(d.at("degree"),"2"),222,82,192);
    s.centers=field(s,L"Gaussian / partition centers",shown(d.at("centers"),"32"),426,82,190);
    s.bandwidth=field(s,L"Gaussian bandwidth",shown(d.at("bandwidth"),"1"),18,146,192);
    s.normalize=makeCheck(s.window,0,L"Normalize training coordinates",s.fonts->body,d.at("normalize").boolean());place(s,s.normalize,222,165,270,28);
    s.constant=makeCheck(s.window,0,L"Include constant",s.fonts->body,d.at("constant").boolean(true));place(s,s.constant,426,199,190,28);
    s.invariant=makeCheck(s.window,0,L"mpEDMD: assume invariant measure and measure-preserving dynamics",s.fonts->body,settings.at("assume_measure_preserving").boolean());place(s,s.invariant,18,238,600,28);
    s.hermitian=makeCheck(s.window,0,L"Hermitian DMD: assume a self-adjoint Koopman operator",s.fonts->body,settings.at("assume_self_adjoint").boolean());place(s,s.hermitian,18,272,600,28);
    label(s,L"Additional settings: observables, Hankel window, weights, residual grid, measures",18,314,600);
    s.advanced=makeEdit(s.window,0,dynsyslab::core::wide(s.current).c_str(),s.fonts->mono,true);place(s,s.advanced,18,338,598,170);
    auto ok=makeButton(s.window,IDOK,L"Use settings",s.fonts->body);place(s,ok,360,526,124,32);
    auto cancel=makeButton(s.window,IDCANCEL,L"Cancel",s.fonts->body);place(s,cancel,496,526,120,32);
}
void accept(State& s) {
    try {
        auto settings=j::parse(utf8Text(s.advanced));
        if(!std::holds_alternative<j::Value::Object>(settings.storage()))throw std::invalid_argument("Additional settings must be an object.");
        auto number=[](HWND item,bool integer) {
            auto value=j::parse(utf8Text(item));const double n=value.number(-1);
            if(!std::isfinite(n)||n<0||(integer&&std::floor(n)!=n))throw std::invalid_argument("Enter valid nonnegative settings.");
            return value;
        };
        settings.at("method")=methods.at(std::size_t(SendMessageW(s.method,CB_GETCURSEL,0,0)));
        auto& d=settings.at("dictionary");d.at("kind")=dictionaries.at(std::size_t(SendMessageW(s.dictionary,CB_GETCURSEL,0,0)));
        settings.at("rank")=number(s.rank,true);d.at("degree")=number(s.degree,true);d.at("centers")=number(s.centers,true);
        d.at("bandwidth")=number(s.bandwidth,false);d.at("constant")=checked(s.constant);d.at("normalize")=checked(s.normalize);
        settings.at("assume_measure_preserving")=checked(s.invariant);settings.at("assume_self_adjoint")=checked(s.hermitian);
        if(settings.at("method").string()=="specrkhs")settings.at("bandwidth")=d.at("bandwidth");
        s.output=j::dump(settings);DestroyWindow(s.window);
    }catch(const std::exception& e){MessageBoxW(s.window,dynsyslab::core::wide(e.what()).c_str(),L"Spectral Koopman settings",MB_OK|MB_ICONWARNING);}
}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* state=reinterpret_cast<State*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE){state=static_cast<State*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);state->window=window;SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(state));}
    if(!state)return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE){try{create(*state);}catch(...){return -1;}return 0;}
    if(message==WM_COMMAND){if(LOWORD(wparam)==IDOK){accept(*state);return 0;}if(LOWORD(wparam)==IDCANCEL){DestroyWindow(window);return 0;}}
    if(message==WM_CLOSE){DestroyWindow(window);return 0;}
    return DefWindowProcW(window,message,wparam,lparam);
}
}
std::optional<std::string> showSpectralKoopmanDialog(HWND owner,const std::string& current,const theme::Fonts& fonts) {
    WNDCLASSEXW cls{sizeof(cls)};cls.lpfnWndProc=procedure;cls.hInstance=GetModuleHandleW(nullptr);
    cls.hCursor=LoadCursorW(nullptr,IDC_ARROW);cls.hbrBackground=CreateSolidBrush(theme::panel);cls.lpszClassName=windowClass;
    if(!RegisterClassExW(&cls)&&GetLastError()!=ERROR_CLASS_ALREADY_EXISTS)return std::nullopt;
    State state;state.current=current;state.fonts=&fonts;RECT area{};GetWindowRect(owner,&area);
    const int width=scalePixels(owner,650),height=scalePixels(owner,610);
    auto window=CreateWindowExW(WS_EX_DLGMODALFRAME,windowClass,L"Spectral Koopman methods",WS_CAPTION|WS_SYSMENU|WS_POPUP,
        area.left+(area.right-area.left-width)/2,area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,cls.hInstance,&state);
    if(!window)return std::nullopt;theme::configureTitleBar(window);EnableWindow(owner,FALSE);ShowWindow(window,SW_SHOW);SetFocus(state.method);
    MSG event{};while(IsWindow(window)&&GetMessageW(&event,nullptr,0,0)>0)if(!IsDialogMessageW(window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}
}
