#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"
#include <dynsyslab/core/Observable.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
void require(bool ok,const char* why) { if (!ok) throw std::runtime_error(why); }
void checkNear(double a,double b) { require(std::abs(a-b)<1e-9*std::max(1.,std::abs(b)),"scale round trip or geometry mismatch"); }
dsl::Figure wideValues() {
    dsl::Figure f{"observable_time","Wide values","t","value",false};
    dsl::Trace t;
    for (int i=0;i<=1000;++i) { t.x.push_back(i*.01); t.y.push_back(std::pow(10.,-4+i*.008)); }
    f.traces.push_back(t); return f;
}
void scales() {
    auto f=wideValues(); const auto raw=f.traces[0].y;
    require(dsl::resolveValueScale(f).name()=="log","positive wide data should select log");
    auto limits=dsl::automaticAxisLimits(f);
    require(limits.yMin>0 && limits.yMin<raw.front() && limits.yMax>raw.back(),"log padding must stay positive");
    f.traces[0].y[0]=0;
    require(dsl::resolveValueScale(f).name()=="symlog","zero must remain visible");
    f.traces[0].y[1]=-10;
    f.valueScale="log";
    auto s=dsl::resolveValueScale(f);
    require(s.name()=="symlog" && s.reason.find("retain")!=std::string::npos,"explicit log must retain nonpositive values");
    for (const double v : {-1e4,-1.,0.,1e-4,1.,1e4}) checkNear(s.inverse(s.forward(v)),v);
    f.valueScale="linear"; require(dsl::resolveValueScale(f).name()=="linear","manual linear override");
    f=wideValues(); const auto prior=dsl::resolveValueScale(f);
    for(auto& y:f.traces[0].y) y*=1e-9;
    require(dsl::resolveValueScale(f).name()==prior.name(),"automatic decision must be unit invariant");
    f.traces[0].y.assign(1001,1e-20); require(dsl::resolveValueScale(f).name()=="linear","tiny constant is not a wide scale");
    f.traces[0].y.assign(1001,1); f.traces[0].y[500]=1e9;
    require(dsl::resolveValueScale(f).name()=="linear","one outlier must not force automatic log");
    f=wideValues(); f.kind="phase";require(dsl::resolveValueScale(f).name()=="linear","phase-space coordinates must remain unchanged");
    f=wideValues();f.kind="stem";f.valueScale="log";
    require(dsl::supportsValueScale(f)&&dsl::resolveValueScale(f).name()=="log",
        "positive stem diagrams must support logarithmic ordinates");
    f.kind="krom_heatmap"; f.traces[0].z=raw; f.traces[0].encoding="scalar";
    require(dsl::resolveValueScale(f).name()=="log" && dsl::ordinateScale(f).name()=="linear","heatmap scales values, not phase domain");
    f.traces[0].encoding="domain-color";
    require(!dsl::supportsValueScale(f),"complex hue must not become logarithmic");
    dsl::Payload p; p.figures={wideValues()};dsl::setValueScale(p,"symlog");
    require(p.figures[0].traces[0].y==raw,"scale controls changed stored samples");
    dsl::PlotTransform t{0,0,100,100,0,1,1,100}; t.yScale.kind=dsl::ValueScaleKind::logarithmic;
    checkNear(dsl::plotPixelY(t,10),50);checkNear(dsl::plotDataY(t,50),10);
    dsl::Figure line; dsl::Trace curve; curve.x={0,1};curve.y={1,100};line.traces={curve};
    const auto probe=dsl::nearestCurveProbe(line,t,50,50,2);
    require(bool(probe),"log plot curve probing failed");checkNear(probe->y,10);
    checkNear(dsl::resolveCurveProbe(line,t,0,0,.5)->y,10);
    const auto moved=dsl::panPlotTransform(t,0,50);checkNear(moved.yMin,10);checkNear(moved.yMax,1000);
    for (const auto mode : {dsl::ValueScaleKind::linear, dsl::ValueScaleKind::logarithmic,
                            dsl::ValueScaleKind::symmetric_log}) {
        t.yScale.kind=mode;
        t.yMin=mode==dsl::ValueScaleKind::logarithmic ? 1 : -100;
        const auto anchor=dsl::plotDataY(t,74);
        const auto z=dsl::zoomPlotLimits(t,21,74,0.8);
        require(bool(z),"wheel limits missing");
        auto zoomed=t;
        zoomed.xMin=z->xMin; zoomed.xMax=z->xMax; zoomed.yMin=z->yMin; zoomed.yMax=z->yMax;
        checkNear(zoomed.xMin+0.21*(zoomed.xMax-zoomed.xMin),0.21);
        checkNear(dsl::plotDataY(zoomed,74),anchor);
        const auto restored=dsl::zoomPlotLimits(zoomed,21,74,1.25);
        require(bool(restored),"inverse zoom limits missing");
        checkNear(restored->yMin,t.yMin); checkNear(restored->yMax,t.yMax);
    }
}
void identities() {
    using dynsyslab::core::Observable;
    const auto x=Observable::stateComponent(2,0,"x","x");
    const auto y=Observable::stateComponent(2,1,"y","y");
    const auto full=Observable::identity(2,{"x","y"},{"x","y"});
    require(full.stateComponentIndex(0)==x.stateComponentIndex(0) && full.stateComponentIndex(1)==y.stateComponentIndex(0),"factory function identities are missing");
    Observable fake({"state","Full state","","native",{"x","y"},false},2,2,
        [](const double* q) { return Observable::Values{q[0],q[1]*q[1]}; });
    const double q[]={2,1}; require(fake.evaluate(q)==full.evaluate(q),"coincident-value fixture");
    require(!fake.stateComponentIndex(0) && !fake.stateComponentIndex(1),"custom function was inferred from labels or samples");
    const auto reordered=Observable::stack("permuted","permuted",{y,x});
    require(reordered.stateComponentIndex(0)==1 && reordered.stateComponentIndex(1)==0,"stack identities ignored indexing");
    const auto mixed=Observable::stack("mixed","mixed",{x,fake});
    require(mixed.stateComponentIndex(0)==0 && !mixed.stateComponentIndex(1),"unknown evaluator acquired a false identity");
}
void scrolling() {
    ULONG_PTR token{};Gdiplus::GdiplusStartupInput input;Gdiplus::GdiplusStartup(&token,&input,nullptr);
    {
        auto parent=CreateWindowExW(WS_EX_NOACTIVATE,L"STATIC",L"",WS_POPUP,0,0,1100,760,nullptr,nullptr,GetModuleHandleW(nullptr),nullptr);
        dsl::theme::Fonts fonts;dsl::PlotCanvas c;dsl::PlotDragToolbar toolbar;
        require(dsl::PlotCanvas::registerClass(GetModuleHandleW(nullptr)) && c.create(parent,fonts) && toolbar.create(parent,c,fonts),"invisible test surface");
        auto p=std::make_shared<dsl::Payload>();p->view="observable_workspace";
        for(int i=0;i<16;++i) {auto f=wideValues();f.title=std::to_string(i);p->figures.push_back(f);}
        c.setPayload(p);MoveWindow(c.handle(),0,80,1100,680,FALSE);MoveWindow(toolbar.handle(),0,0,1100,80,FALSE);toolbar.updateCoordinates();
        Gdiplus::Bitmap bitmap(1100,680,PixelFormat32bppARGB);Gdiplus::Graphics g(&bitmap);
        auto paint=[&] {const auto dc=g.GetHDC();SendMessageW(c.handle(),WM_PRINTCLIENT,reinterpret_cast<WPARAM>(dc),PRF_CLIENT);g.ReleaseHDC(dc);};
        require(!c.setAxisLimits(0,dsl::AxisLimits{0,10,0,1e5}),"nonpositive manual log bounds should be rejected");
        paint(); require(c.pageCount()==1 && c.scrollMaximum()>1000,"all observable panels must scroll");
        const auto first=c.visibleFigureIndices();require(!first.empty() && first.front()==0 && first.back()<15,"first viewport wrong");
        Gdiplus::Color white;bitmap.GetPixel(30,60,&white);require(white.GetR()==255 && white.GetG()==255 && white.GetB()==255,"plot card must be white");
        const auto initialStats=c.renderStatistics();
        for (int i=1;i<=5;++i) {c.setScrollPosition(i*4);paint();}
        require(c.renderStatistics().panelRenders==initialStats.panelRenders
            && c.renderStatistics().cacheHits>=initialStats.cacheHits+first.size()*5,
            "scrolling rebuilt already-rendered plot panels");
        require(c.renderStatistics().cacheBytes<=64u*1024u*1024u,"panel cache exceeded its memory budget");
        c.setValueScaling("linear");paint();
        require(c.renderStatistics().panelRenders>initialStats.panelRenders,"scale change reused stale panel images");
        const auto beforeAxes=c.renderStatistics().panelRenders;
        c.setAxisLimits(0,dsl::AxisLimits{0,10,0,1e5});paint();
        require(c.renderStatistics().panelRenders==beforeAxes+1,"axis change must redraw only its affected panel");
        c.setAxisLimits(0,std::nullopt);c.setValueScaling("auto");
        c.setScrollPosition(c.scrollMaximum());paint();
        auto last=c.visibleFigureIndices();require(last.back()==15 && last.front()>0,"last observable is unreachable");
        const auto saved=c.captureState();c.setPayload(std::make_shared<dsl::Payload>());c.restoreState(saved);paint();
        require(c.scrollPosition()==saved.scrollY && c.visibleFigureIndices()==last,"tab round trip lost scroll position");
        POINT plotPoint{180,155}, wheel=plotPoint; ClientToScreen(c.handle(),&wheel);
        for (int i=0;i<2;++i)
            SendMessageW(c.handle(),WM_MOUSEWHEEL,MAKEWPARAM(0,WHEEL_DELTA/2),MAKELPARAM(wheel.x,wheel.y));
        const auto zoom=c.captureState();
        require(zoom.scrollY==saved.scrollY && zoom.axisLimits.size()==1 && zoom.axisLimits.count(last.front()),
            "wheel must zoom only the hovered scrolled panel");
        const auto autoBounds=c.currentAxisPanels()[last.front()].automatic;
        const auto manual=zoom.axisLimits.at(last.front());
        checkNear((manual.xMax-manual.xMin)/(autoBounds.xMax-autoBounds.xMin),.8);
        SendMessageW(c.handle(),WM_LBUTTONDBLCLK,0,MAKELPARAM(plotPoint.x,plotPoint.y));
        require(c.captureState().axisLimits.empty(),"double-click should restore automatic bounds");
        paint();
        wheel={2,100};ClientToScreen(c.handle(),&wheel);
        SendMessageW(c.handle(),WM_MOUSEWHEEL,MAKEWPARAM(0,WHEEL_DELTA),MAKELPARAM(wheel.x,wheel.y));
        require(c.scrollPosition()<saved.scrollY && c.scrollPosition()>saved.scrollY-120,
            "mouse wheel must begin moving without jumping the whole notch");
        const auto deadline=GetTickCount64()+300;
        while(GetTickCount64()<deadline) {
            MSG event{};while(PeekMessageW(&event,nullptr,0,0,PM_REMOVE)) {TranslateMessage(&event);DispatchMessageW(&event);}
            Sleep(5);
        }
        require(c.scrollPosition()==saved.scrollY-120,"smooth wheel motion did not settle at the requested position");
        c.setScrollPosition(saved.scrollY);paint();
        const auto firstVisible=c.visibleFigureIndices().front();
        const auto firstTitle=c.captureState().payload->figures[firstVisible].title;
        SendMessageW(c.handle(),WM_LBUTTONDOWN,0,MAKELPARAM(20,40));
        SendMessageW(c.handle(),WM_MOUSEMOVE,0,MAKELPARAM(1040,650));
        SendMessageW(c.handle(),WM_LBUTTONUP,0,MAKELPARAM(1040,650));
        const auto reordered=c.captureState().payload;
        require(reordered->figures.size()==16 && reordered->figures.back().title==firstTitle,"scrolled reorder used a visible-list index instead of the figure index");
        for(int i=10;i<=13;++i)require((GetWindowLongPtrW(GetDlgItem(toolbar.handle(),i),GWL_STYLE)&BS_TYPEMASK)==BS_AUTORADIOBUTTON,"scale control must be a round radio");
        SendMessageW(GetDlgItem(toolbar.handle(),11),BM_CLICK,0,0);require(c.valueScaling()=="linear","scale radio did not update canvas");
        require(!IsWindowVisible(parent),"surface became visible");DestroyWindow(parent);
    }
    Gdiplus::GdiplusShutdown(token);
}
}
int main() { try {scales();identities();scrolling();std::cout<<"scales, exact identities, scrolling: PASS\n";return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;} }
