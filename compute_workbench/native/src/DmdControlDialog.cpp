#include "DmdControlDialog.hpp"
#include "SettingsMerge.hpp"
#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Json.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <cmath>
#include <stdexcept>

namespace dsl::compute {
namespace {
constexpr wchar_t kClass[]=L"DynSysLabDmdControlDialog";
struct State {
    HWND window{},snapshot{},rank{},augmented{},tolerance{},save{};
    HWND degree{},constant{},normalize{};
    bool extended{};
    const dsl::theme::Fonts* fonts{};
    std::string current;
    std::optional<std::string> output;
};
void place(State& s,HWND control,int x,int y,int w,int h) {
    MoveWindow(control,scalePixels(s.window,x),scalePixels(s.window,y),
        scalePixels(s.window,w),scalePixels(s.window,h),TRUE);
}
void text(State& s,const wchar_t* value,int x,int y,int w,int h,HFONT font) {
    place(s,makeLabel(s.window,value,font),x,y,w,h);
}
HWND field(State& s,const wchar_t* label,const std::string& value,int x,int y,int w) {
    text(s,label,x,y,w,20,s.fonts->caption);
    const auto item=makeEdit(s.window,0,dsl::wide(value).c_str(),s.fonts->mono);
    place(s,item,x,y+24,w,28); return item;
}
void create(State& s) {
    const auto settings=json::parse(jsonObject(s.current));
    auto number=[&](const char* key) {return settings.at(key).isNull()?"0":json::dump(settings.at(key));};
    text(s,s.extended?L"Polynomial EDMD with control":L"DMD with control",20,16,560,28,s.fonts->heading);
    text(s,s.extended?L"Lift recorded coordinates with the SINDy monomial library; jointly fit A and B."
        :L"Fit q[k+1] = A q[k] + B u[k] from recorded measurements and inputs.",
        20,52,560,36,s.fonts->body);
    s.snapshot=field(s,L"Snapshot array in the source archive",settings.at("snapshot_array").string("values"),20,100,560);
    text(s,L"Use stored states, observables, or precomputed latent coordinates. No equations of motion are evaluated.",
        20,160,560,40,s.fonts->caption);
    s.augmented=field(s,s.extended?L"Augmented rank [Theta; u]":L"Augmented rank [q; u]",number("augmented_rank"),20,216,176);
    s.rank=field(s,L"Output rank",number("rank"),212,216,176);
    s.tolerance=field(s,L"Absolute SVD tolerance",number("svd_tolerance"),404,216,176);
    text(s,L"Zero selects numerical rank / automatic tolerance. Inputs must be aligned to the start of each sampling interval.",
        20,280,560,40,s.fonts->caption);
    const int extra=s.extended?104:0;
    if(s.extended) {
        s.degree=field(s,L"Polynomial degree (1-32)",settings.at("polynomial_degree").isNull()?"4":number("polynomial_degree"),20,328,176);
        s.constant=makeCheck(s.window,0,L"Include constant term",s.fonts->body,settings.at("include_constant").boolean(true));
        place(s,s.constant,212,328,368,28);
        s.normalize=makeCheck(s.window,0,L"Normalize coordinates (frozen after fit)",s.fonts->body,settings.at("normalize").boolean(true));
        place(s,s.normalize,212,364,368,28);
        text(s,L"Includes mixed terms through total degree; degree 4 in two variables gives 15 features.",20,396,560,28,s.fonts->caption);
    }
    s.save=makeCheck(s.window,0,s.extended?L"Save sampled polynomial eigenfunctions for inspection":L"Save sampled linear eigenfunctions for inspection",s.fonts->body,
        settings.at("save_eigenfunction_values").boolean(true));
    place(s,s.save,20,328+extra,560,28);
    text(s,s.extended?L"Output: A, B, spectrum and polynomial coefficients. Constant additive B is an approximation after lifting; validate on held-out inputs. Not a continuous-time B_psi."
        :L"Output: reduced A and B, projection, spectrum, linear eigenfunctions and fit diagnostics. This is not a continuous-time B_psi or a control-ready KROM.",
        20,368+extra,560,56,s.fonts->caption);
    for(const auto id:{IDOK,IDCANCEL}) {
        const auto button=CreateWindowExW(0,L"BUTTON",id==IDOK?L"Apply settings":L"Cancel",
            WS_CHILD|WS_VISIBLE|WS_TABSTOP|(id==IDOK?BS_DEFPUSHBUTTON:BS_PUSHBUTTON),
            0,0,0,0,s.window,reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),GetModuleHandleW(nullptr),nullptr);
        dsl::theme::styleControl(button,s.fonts->body);
        place(s,button,id==IDOK?300:448,444+extra,132,32);
    }
}
void accept(State& s) {
    try {
        auto numeric=[](HWND item,bool integer) {
            const auto value=json::parse(utf8Text(item)); const double n=value.number(-1);
            if(!std::isfinite(n)||n<0||(integer&&(n!=std::floor(n)||n>1000000000)))
                throw std::invalid_argument("Ranks must be nonnegative integers; tolerance must be a nonnegative finite number.");
            return json::dump(value);
        };
        if(utf8Text(s.snapshot).empty()) throw std::invalid_argument("Specify the snapshot array name.");
        s.output=mergeSettingsObject(s.current,{
            {"snapshot_array","\""+jsonEscape(utf8Text(s.snapshot))+"\""},
            {"rank",numeric(s.rank,true)},{"augmented_rank",numeric(s.augmented,true)},
            {"svd_tolerance",numeric(s.tolerance,false)},
            {"save_eigenfunction_values",checked(s.save)?"true":"false"}});
        if(s.extended) {
            const auto degree=numeric(s.degree,true);
            const double value=json::parse(degree).number();
            if(value<1||value>32) throw std::invalid_argument("Polynomial degree must be between 1 and 32.");
            s.output=mergeSettingsObject(*s.output,{{"dictionary","\"polynomial\""},
                {"polynomial_degree",degree},{"include_constant",checked(s.constant)?"true":"false"},
                {"normalize",checked(s.normalize)?"true":"false"}});
        }
        DestroyWindow(s.window);
    } catch(const std::exception& e) {
        s.output.reset();
        MessageBoxW(s.window,dsl::wide(e.what()).c_str(),s.extended?L"EDMDc settings":L"DMDc settings",MB_OK|MB_ICONWARNING);
    }
}
LRESULT CALLBACK procedure(HWND window,UINT message,WPARAM wparam,LPARAM lparam) {
    auto* s=reinterpret_cast<State*>(GetWindowLongPtrW(window,GWLP_USERDATA));
    if(message==WM_NCCREATE) {
        s=static_cast<State*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        s->window=window; SetWindowLongPtrW(window,GWLP_USERDATA,reinterpret_cast<LONG_PTR>(s));
    }
    if(!s) return DefWindowProcW(window,message,wparam,lparam);
    if(message==WM_CREATE) {
        try {create(*s);} catch(...) {return -1;} return 0;
    }
    if(message==WM_COMMAND) {
        if(LOWORD(wparam)==IDOK) {accept(*s); return 0;}
        if(LOWORD(wparam)==IDCANCEL) {DestroyWindow(window); return 0;}
    }
    if(message==WM_CLOSE) {DestroyWindow(window); return 0;}
    if(message==WM_CTLCOLORSTATIC) {
        static HBRUSH brush=CreateSolidBrush(theme::panel);
        const auto dc=reinterpret_cast<HDC>(wparam);
        SetBkColor(dc,theme::panel); SetTextColor(dc,theme::text); return reinterpret_cast<LRESULT>(brush);
    }
    return DefWindowProcW(window,message,wparam,lparam);
}
} // namespace
std::optional<std::string> showDmdControlDialog(HWND owner,
    const std::string& current,const dsl::theme::Fonts& fonts,bool extended) {
    WNDCLASSEXW cls{sizeof(cls)}; cls.lpfnWndProc=procedure; cls.hInstance=GetModuleHandleW(nullptr);
    cls.hCursor=LoadCursorW(nullptr,IDC_ARROW); cls.hbrBackground=CreateSolidBrush(theme::panel); cls.lpszClassName=kClass;
    if(!RegisterClassExW(&cls)&&GetLastError()!=ERROR_CLASS_ALREADY_EXISTS) return std::nullopt;
    State s; s.current=current; s.fonts=&fonts; s.extended=extended;
    RECT area{}; GetWindowRect(owner,&area);
    const int width=scalePixels(owner,620),height=scalePixels(owner,extended?634:530);
    const auto window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClass,extended?L"EDMDc settings":L"DMDc settings",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,cls.hInstance,&s);
    if(!window) return std::nullopt;
    theme::configureTitleBar(window); EnableWindow(owner,FALSE); ShowWindow(window,SW_SHOW);
    SetFocus(s.snapshot);
    MSG event{};
    while(IsWindow(window)&&GetMessageW(&event,nullptr,0,0)>0)
        if(!IsDialogMessageW(window,&event)) {TranslateMessage(&event); DispatchMessageW(&event);}
    EnableWindow(owner,TRUE); SetForegroundWindow(owner); return s.output;
}
} // namespace dsl::compute
