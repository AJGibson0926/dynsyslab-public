#include "dynsyslab_viewer/KromWorkspace.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
void checkNear(double a, double b, const char* message) { require(std::abs(a-b) < 1e-10,message); }
void click(HWND w, int x, int y, bool control = false) {
    SendMessageW(w,WM_LBUTTONDOWN,control ? MK_CONTROL : 0,MAKELPARAM(x,y));
    SendMessageW(w,WM_LBUTTONUP,0,MAKELPARAM(x,y));
}
void interaction(dsl::Payload payload) {
    ULONG_PTR token{}; Gdiplus::GdiplusStartupInput input;
    Gdiplus::GdiplusStartup(&token,&input,nullptr);
    {
        const auto instance = GetModuleHandleW(nullptr);
        HWND parent = CreateWindowExW(WS_EX_NOACTIVATE,L"STATIC",L"",WS_POPUP,0,0,1400,900,nullptr,nullptr,instance,nullptr);
        dsl::theme::Fonts fonts; dsl::PlotCanvas canvas; dsl::PlotDragToolbar toolbar;
        require(dsl::PlotCanvas::registerClass(instance) && canvas.create(parent,fonts) && toolbar.create(parent,canvas,fonts),"hidden controls");
        auto p = std::make_shared<dsl::Payload>(std::move(payload));
        // Separate synthetic marker coordinates make mouse-selection geometry unambiguous.
        p->figures[0].traces[1].x = {0}; p->figures[0].traces[1].y = {0};
        p->figures[0].traces[2].x = {.5}; p->figures[0].traces[2].y = {.5};
        canvas.setPayload(p); MoveWindow(canvas.handle(),0,80,1400,820,FALSE);
        MoveWindow(toolbar.handle(),0,0,1400,80,FALSE); toolbar.updateCoordinates();
        SendMessageW(canvas.handle(),WM_PAINT,0,0);
        // Search the right spectrum with clicks; no screen/input automation is involved.
        canvas.selectAllKromModes(false);
        POINT first{},second{}; bool found0=false,found1=false;
        for (int y=150;y<650 && !(found0 && found1);y+=7) for(int x=1030;x<1360;x+=7) {
            click(canvas.handle(),x,y);
            if (p->selectedModes == std::vector<std::size_t>{0} && !found0) { first={x,y};found0=true; }
            if (p->selectedModes == std::vector<std::size_t>{1} && !found1) { second={x,y};found1=true; }
            canvas.selectAllKromModes(false);
        }
        require(found0 && found1,"right-hand eigenvalues must be selectable");
        click(canvas.handle(),first.x,first.y); click(canvas.handle(),second.x,second.y,true);
        require(p->selectedModes.size()==2,"Ctrl-click must add a mode");
        click(canvas.handle(),first.x,first.y,true);
        require(p->selectedModes == std::vector<std::size_t>{1},"Ctrl-click must remove a mode");
        const auto real = GetDlgItem(toolbar.handle(),5), mag = GetDlgItem(toolbar.handle(),6), domain=GetDlgItem(toolbar.handle(),7);
        for (const auto h : {real,mag,domain}) require((GetWindowLongPtrW(h,GWL_STYLE)&BS_TYPEMASK)==BS_AUTORADIOBUTTON,"round complex radios");
        SendMessageW(domain,BM_CLICK,0,0);
        require(p->complexDisplay=="domain-coloring" && SendMessageW(domain,BM_GETCHECK,0,0)==BST_CHECKED &&
            SendMessageW(real,BM_GETCHECK,0,0)==BST_UNCHECKED && SendMessageW(mag,BM_GETCHECK,0,0)==BST_UNCHECKED,"exclusive domain radio");
        POINT wheel=first; ClientToScreen(canvas.handle(),&wheel);
        SendMessageW(canvas.handle(),WM_MOUSEWHEEL,MAKEWPARAM(0,WHEEL_DELTA),MAKELPARAM(wheel.x,wheel.y));
        const auto zoomed=canvas.captureState();
        require(zoomed.axisLimits.count(0) && zoomed.axisLimits.at(0).xMax-zoomed.axisLimits.at(0).xMin < 2.18,
            "wheel zooms the spectrum instead of changing pages");
        require(canvas.page()==0,"spectrum wheel changed the eigenfunction page");
        SendMessageW(canvas.handle(),WM_LBUTTONDBLCLK,0,MAKELPARAM(first.x,first.y));
        require(!canvas.captureState().axisLimits.count(0),"double-click resets spectrum bounds");
        require(canvas.verticallyScrollable() && canvas.pageCount()==1 && canvas.scrollMaximum()>0,"KROM must scroll in a single workspace");
        canvas.setScrollPosition(canvas.scrollMaximum()); SendMessageW(canvas.handle(),WM_PAINT,0,0);
        require(canvas.visibleFigureIndices().back()==0,"spectrum remains fixed while scrolling");
        click(canvas.handle(),first.x,first.y);
        require(p->selectedModes == std::vector<std::size_t>{0},"spectrum persists while scrolling");
        const auto count = p->figures.size();
        // Border gesture rearranges only left panels and retains all other pages and spectrum.
        SendMessageW(canvas.handle(),WM_LBUTTONDOWN,0,MAKELPARAM(10,100));
        SendMessageW(canvas.handle(),WM_MOUSEMOVE,0,MAKELPARAM(600,150));
        SendMessageW(canvas.handle(),WM_LBUTTONUP,0,MAKELPARAM(600,150));
        require(p->figures.size()==count && p->figures[0].kind=="krom_spectrum" && canvas.pageCount()==1,"rearranging preserves KROM spectrum and scroll contents");
        require(!IsWindowVisible(parent),"test surface became visible");
        DestroyWindow(parent);
    }
    Gdiplus::GdiplusShutdown(token);
}
}
int main(int argc,char** argv) {
    try {
        const auto root=argc>1 ? std::filesystem::path(argv[1]) : std::filesystem::temp_directory_path()/"dynsyslab-krom-workspace-test";
        std::filesystem::create_directories(root);
        SystemConfiguration system; system.identifier="duffing"; system.preset="canonical";
        InterpolationData data;
        data.samples=4; data.stateDimension=2; data.neighbors=4;
        data.trainingStates={-1,-1,1,-1,-1,1,1,1};
        data.trainingValues={{-1,-1},{-1,1},{1,-1},{1,1},{-1,1},{-1,-1},{1,1},{1,-1}};
        KoopmanModel model(data,Spectrum::fromContinuousRates({{-.1,1},{-.1,-1}}),system);
        const auto source=root/"scattered.npz";
        saveKoopmanModel(source,model,R"({"label":"Scattered fixture","delta_t":0.1})",
            {{"training_time",NpyArray::fromDoubles({4},{0,.1,0,.1})}});
        dsl::KromWorkspaceOptions options;options.allModes=true;
        const auto real=dsl::buildKromWorkspace(source,options);
        require(real.figures.size()==9 && real.selectedModes.size()==2,"selected component and map count");
        require(real.figures[2].traces.size()==2 && real.figures[2].traces[1].x.front()==0,"trajectory boundary retained");
        checkNear(real.figures[2].traces[0].y[0],-1,"scattered training values");
        const auto& grid=real.figures[1].traces[0];
        const double q[]={grid.x[32],grid.y[48]};
        checkNear(grid.z[48*101+32],model.evaluate(q)[0].real(),"grid uses scattered evaluator");
        options.complexDisplay="magnitude-phase";
        const auto polar=dsl::buildKromWorkspace(source,options);
        checkNear(polar.figures[2].traces[0].y[0],std::sqrt(2.),"magnitude component");
        checkNear(polar.figures[4].traces[0].y[0],std::arg(data.trainingValues[0]),"phase component");
        options.complexDisplay="domain-coloring";
        const auto domain=dsl::buildKromWorkspace(source,options);
        require(domain.figures.size()==7 && domain.figures[1].traces[0].encoding=="domain-color","one domain field and two histories per mode");
        require(domain.figures[1].rowSpan==2 && domain.figures[2].traces[0].encoding.empty()
            && domain.figures[3].traces[0].encoding.empty(),"domain field spans magnitude and phase histories");
        checkNear(domain.figures[2].traces[0].y[0],std::sqrt(2.),"domain dashboard magnitude history");
        checkNear(domain.figures[3].traces[0].y[0],std::arg(data.trainingValues[0]),"domain dashboard phase history");
        require(domain.figures[4].traces[0].x==grid.x && domain.figures[4].traces[0].y==grid.y,"representation preserves phase domain");
        for(const auto& f:domain.figures) if(f.kind=="krom_heatmap") require(f.references.empty(),"no heatmap markers");
        options.allModes=false; options.modes=std::vector<std::size_t>{};
        require(dsl::buildKromWorkspace(source,options).figures.size()==1,"none retains spectrum");
        options.modes=std::vector<std::size_t>{2};
        bool rejected=false;try{dsl::buildKromWorkspace(source,options);}catch(const std::out_of_range&){rejected=true;}
        require(rejected,"invalid retained mode rejected");
        std::ofstream(root/"scattered.viewer.json")<<dsl::json::dump(dsl::inspectPlot(real,true));
        const auto roundtrip=dsl::Payload::load((root/"scattered.viewer.json").wstring());
        require(roundtrip.layout=="krom" && roundtrip.selectedModes==real.selectedModes && roundtrip.figures[0].traces[1].selected,"payload selection roundtrip");
        {   // Each eigenvalue shows its pair's residual; a lattice pseudospectrum is a heat map beneath them.
            auto evidence=model.spectrum();
            evidence.setResiduals({0.0123,0.0123},"one-step","empirical-training");
            evidence.setPseudospectrum(Pseudospectrum{"one-step","empirical-training",{{0.9,-0.1},{1.0,-0.1},{0.9,0.1},{1.0,0.1}},{0.1,0.01,0.2,0.02},{},0});
            auto withEvidence=model; withEvidence.setSpectralEvidence(evidence);
            const auto evidenced=root/"evidenced.npz";
            saveKoopmanModel(evidenced,withEvidence,R"({"label":"Evidence fixture","delta_t":0.1})",
                {{"training_time",NpyArray::fromDoubles({4},{0,.1,0,.1})}});
            dsl::KromWorkspaceOptions all; all.allModes=true;
            const auto shown=dsl::buildKromWorkspace(evidenced,all);
            const auto& spectrum=shown.figures[0];
            require(spectrum.traces[0].heatmap && spectrum.traces[0].zShape==std::vector<int>{2,2}
                && std::abs(spectrum.traces[0].z[1]+2)<1e-12,"pseudospectrum heat map of log10 residual");
            require(spectrum.traces.back().label.find("r = 0.0123")!=std::string::npos,"eigenvalue label carries its residual");
            dsl::renderPlotPng(shown,root/"evidence.png",{1400,900});
        }
        interaction(real);
        dsl::renderPlotPng(domain,root/"domain.png",{1400,900});
        std::cout<<"KROM workspace: PASS\n"; return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
