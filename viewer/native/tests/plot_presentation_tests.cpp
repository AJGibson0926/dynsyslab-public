#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotDragToolbar.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include "../src/PlotSeries.hpp"
#include <gdiplus.h>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
namespace {
void require(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
bool closeEnough(double a,double b) { return std::abs(a-b)<=1e-10*std::max(1.,std::abs(b)); }
dsl::Payload fixture(const std::string& view = "data_overview", bool adaptive = false) {
    dynsyslab::core::NpzRecord record;
    constexpr std::size_t count=10001;
    std::vector<double> time(count),states(count*2),controls(count),targets(count*2);
    for(std::size_t i=0;i<count;++i) {
        const double t=i*.001;
        time[i]=t*3e-6;
        states[i*2]=(1+.4*std::exp(-.1*t)*std::cos(4*t))*1e-5;
        states[i*2+1]=-.4*4*std::exp(-.1*t)*std::sin(4*t)*1e-5/3e-6;
        targets[i*2]=1e-5;targets[i*2+1]=0;
        controls[i]=.1*std::sin(4*t);
    }
    using dynsyslab::core::NpyArray;
    record.arrays["time"]=NpyArray::fromDoubles({count},time);
    record.arrays["values"]=NpyArray::fromDoubles({count,2},states);
    record.arrays["control_target_values"]=NpyArray::fromDoubles({count,2},targets);
    record.arrays["control_values"]=NpyArray::fromDoubles({count,1},controls);
    if (adaptive) {
        record.arrays["raw_time"] = NpyArray::fromDoubles({3}, {time.front(),1.5e-9,time.back()});
        record.arrays["raw_values"] = NpyArray::fromDoubles({3,2},
            {states[0],states[1],1.9e-5,4.,states[states.size()-2],states.back()});
    }
    const auto metadata=dsl::json::parse(R"({"system":"rayleigh-plesset","integration_form":"dimensionless",
        "archive_coordinate_system":"dimensional","system_parameters":{"equilibrium_radius_m":1e-5},
        "nondimensionalization":{"time_scale_s":3e-6},"control_input_coordinate_system":"dimensionless",
        "control_input_scale_to_dimensional":100000})");
    return dsl::buildNativePayload("fixture.npz",view,record,metadata);
}
void checkData(dsl::Payload payload) {
    require(payload.figures.size()==4,"dashboard changed its figure count");
    auto& phase=payload.figures.back(); const auto original=phase;
    require(!phase.square&&phase.kind=="phase","phase panel is still square");
    require(phase.traces[0].x.size()==10001,"trajectory samples were discarded");
    require(phase.traces[0].encoding=="time"&&closeEnough(phase.traces[0].colorValues.back(),10),
        "phase time colors do not align with samples");
    require(phase.references.size()==2&&phase.references.back().equilibrium
        &&closeEnough(phase.references.back().x,1),"known RPE equilibrium is missing");
    require(phase.references.front().label=="IC"&&phase.references.back().label=="EC",
        "phase reference points are not labeled IC and EC");
    require(phase.yLabel=="\\frac{dr}{d\\tau}","velocity label is not LaTeX");
    require(payload.figures[0].xLabel=="\\tau","time axis is not LaTeX");
    for(auto& figure:payload.figures) {
        require(figure.kind == "phase" || figure.references.empty(),
            "time-series and control panels must not contain state reference points");
        require(dsl::canDisplayCoordinates(figure,"dimensional"),"coordinate conversion unavailable");
        dsl::changeDisplayCoordinates(figure,"dimensional");
    }
    require(closeEnough(phase.traces[0].x.front(),1.4e-5)&&closeEnough(phase.traces[1].x.front(),1e-5),
        "physical radius or target is wrong");
    require(closeEnough(phase.traces[0].colorValues.back(),3e-5),"colorbar time units did not change");
    require(closeEnough(phase.references.back().x,1e-5),"equilibrium did not change units");
    require(closeEnough(payload.figures[0].traces[0].y[1000],10000*std::sin(4.)),
        "control pressure scale is wrong");
    dsl::changeDisplayCoordinates(phase,"dimensionless");
    require(closeEnough(phase.traces[0].x[127],original.traces[0].x[127])
        &&closeEnough(phase.traces[0].y[127],original.traces[0].y[127]),"coordinate round trip lost data");
    const auto fraction=dsl::parseMathText(L"\\frac{dr}{d\\tau}");
    require(fraction.size()==1&&fraction[0].role==dsl::MathTextRole::fraction
        &&fraction[0].denominator[0].text==L"d\u03c4","fraction parser failed");
    const auto nested=dsl::parseMathText(L"\\frac{r_{0}}{\\frac{dt}{d\\tau}}");
    require(nested[0].denominator[0].role==dsl::MathTextRole::fraction,"nested fraction parser failed");
    const auto xi=dsl::parseMathText(L"\\xi_{1}");
    require(xi.size()==2&&xi[0].text==L"\u03be"
        &&xi[1].role==dsl::MathTextRole::subscript&&xi[1].text==L"1",
        "xi Greek symbol or subscript parsing failed");
}
void strokeVisibility() {
    Gdiplus::Bitmap bitmap(240,80,PixelFormat32bppARGB);Gdiplus::Graphics g(&bitmap);
    g.SetSmoothingMode(Gdiplus::SmoothingModeAntiAlias);
    dsl::Trace curve;curve.color="#0000ff";curve.label="velocity";curve.x={0,1};curve.y={.5,.5};
    curve.opacity=.4;
    require(dsl::plotSeriesColor(curve,0).GetA()==102,
        "trace opacity is not applied to Viewer series colors");
    curve.opacity=1.;
    dsl::Trace target=curve;target.color="#ff0000";target.role="control_target";target.lineStyle="dash";
    const dsl::PlotTransform transform{10,0,220,80,0,1,0,1};
    for(const auto* encoding:{"","time"}) {
        g.Clear(Gdiplus::Color::White);curve.encoding=encoding;
        if(curve.encoding=="time") curve.colorValues={0,1};else curve.colorValues.clear();
        dsl::drawPlotSeries(g,curve,0,transform);
        int rows=0;for(int y=35;y<=45;++y) {Gdiplus::Color pixel;bitmap.GetPixel(100,y,&pixel);if(pixel.GetR()<200)++rows;}
        require(rows>=2,"trajectory stroke is too thin");
    }
    curve.encoding.clear();curve.colorValues.clear();g.Clear(Gdiplus::Color::White);
    dsl::drawPlotSeries(g,curve,0,transform);dsl::drawPlotSeries(g,target,1,transform);
    int blue=0,red=0,longest=0,gap=0;
    for(int x=20;x<220;++x) {
        Gdiplus::Color pixel;bitmap.GetPixel(x,40,&pixel);
        if(pixel.GetB()>200&&pixel.GetR()<80) {++blue;longest=std::max(longest,++gap);}else gap=0;
        if(pixel.GetR()>200&&pixel.GetB()<80)++red;
    }
    require(blue>80&&red>40&&longest>=8,"target dashes must expose the overlapping trajectory through wide gaps");
}
void referenceLabelPlacement() {
    const Gdiplus::RectF area(10,20,100,60);
    const Gdiplus::SizeF label(16,12);
    const auto right=dsl::plotReferenceLabelBounds({50,45},4,label,area);
    require(closeEnough(right.X,57)&&closeEnough(right.Y,39),
        "phase reference label is not immediately beside and centered on its point");
    const auto left=dsl::plotReferenceLabelBounds({106,45},4,label,area);
    require(closeEnough(left.GetRight(),99)&&closeEnough(left.Y,39),
        "right-edge phase reference label is not immediately beside its point");
    require(left.X>=area.X&&left.GetRight()<=area.GetRight()&&left.Y>=area.Y&&
        left.GetBottom()<=area.GetBottom(),"phase reference label escaped the plot area");
}
void save(HWND canvas,HWND toolbar,const std::filesystem::path& path,int width,int height) {
    MoveWindow(toolbar,0,0,width,40,FALSE);
    MoveWindow(canvas,0,40,width,height-40,FALSE);
    Gdiplus::Bitmap bitmap(width,height,PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);
    HDC dc=graphics.GetHDC();
    SendMessageW(toolbar,WM_PRINT,reinterpret_cast<WPARAM>(dc),PRF_CLIENT|PRF_CHILDREN);
    SetViewportOrgEx(dc,0,40,nullptr);
    SendMessageW(canvas,WM_PRINTCLIENT,reinterpret_cast<WPARAM>(dc),PRF_CLIENT);
    graphics.ReleaseHDC(dc);
    UINT count=0,bytes=0;Gdiplus::GetImageEncodersSize(&count,&bytes);
    std::vector<BYTE> storage(bytes);
    auto* codecs=reinterpret_cast<Gdiplus::ImageCodecInfo*>(storage.data());
    Gdiplus::GetImageEncoders(count,bytes,codecs);
    for(UINT i=0;i<count;++i) if(wcscmp(codecs[i].MimeType,L"image/png")==0) {
        require(bitmap.Save(path.c_str(),&codecs[i].Clsid)==Gdiplus::Ok,"preview export failed");return;
    }
    throw std::runtime_error("PNG encoder unavailable");
}
}
int wmain(int argc,wchar_t** argv) {
    ULONG_PTR token=0;HWND parent=nullptr;
    try {
        auto payload=fixture();checkData(payload);
        const auto adaptive = fixture("phase", true);
        const auto& adaptiveTrace = adaptive.figures.front().traces.front();
        require(adaptiveTrace.x.size() == 10001
            && closeEnough(adaptiveTrace.colorValues[1],.001),
            "plot must use published samples, as the chat figure does");
        for (const auto* view : {"fourier", "control_fourier"}) {
            auto spectrum = fixture(view);
            for (auto& figure : spectrum.figures) {
                const auto before = figure.traces.front();
                require(dsl::canDisplayCoordinates(figure,"dimensional"),"spectrum lost its coordinate conversion");
                require(figure.references.empty(),"spectra must not get state reference points");
                dsl::changeDisplayCoordinates(figure,"dimensional");
                require(closeEnough(figure.traces[0].x[1],before.x[1]/3e-6),"frequency conversion is wrong");
                dsl::changeDisplayCoordinates(figure,"dimensionless");
                require(closeEnough(figure.traces[0].y[1],before.y[1]),"spectral amplitude round trip failed");
            }
        }
        Gdiplus::GdiplusStartupInput input;Gdiplus::GdiplusStartup(&token,&input,nullptr);
        strokeVisibility();referenceLabelPlacement();
        parent=CreateWindowExW(0,L"STATIC",L"",WS_POPUP,0,0,1100,760,nullptr,nullptr,GetModuleHandleW(nullptr),nullptr);
        dsl::theme::Fonts fonts;dsl::PlotCanvas canvas;dsl::PlotDragToolbar toolbar;
        require(dsl::PlotCanvas::registerClass(GetModuleHandleW(nullptr))&&canvas.create(parent,fonts)
            &&toolbar.create(parent,canvas,fonts),"preview window creation failed");
        canvas.setPayload(std::make_shared<dsl::Payload>(payload));toolbar.updateCoordinates();
        const auto physical=FindWindowExW(toolbar.handle(),nullptr,L"BUTTON",L"Dimensional");
        const auto scaled=FindWindowExW(toolbar.handle(),nullptr,L"BUTTON",L"Dimensionless");
        require(IsWindowEnabled(physical)&&IsWindowEnabled(scaled),"coordinate buttons not enabled");
        for (const auto control : {physical,scaled}) {
            const auto style = GetWindowLongPtrW(control,GWL_STYLE);
            require((style & BS_TYPEMASK) == BS_AUTORADIOBUTTON && !(style & BS_PUSHLIKE),
                "coordinate controls must be round native radio buttons");
        }
        canvas.setAxisLimits(3,dsl::AxisLimits{.2,2.4,-3,3});
        SendMessageW(physical,BM_CLICK,0,0);
        require(canvas.displayCoordinates()=="dimensional"
            &&SendMessageW(physical,BM_GETCHECK,0,0)==BST_CHECKED
            &&SendMessageW(scaled,BM_GETCHECK,0,0)==BST_UNCHECKED,"coordinate radio group failed");
        require(closeEnough(canvas.currentAxisPanels()[3].manual->xMax,2.4e-5),"manual axes not converted");
        SendMessageW(scaled,BM_CLICK,0,0);
        require(closeEnough(canvas.currentAxisPanels()[3].manual->xMax,2.4),"manual axes did not round trip");
        if(argc>1) {
            const std::filesystem::path output(argv[1]);std::filesystem::create_directories(output);
            if(argc>2) payload=dsl::buildNativePayload(argv[2],"data_overview");
            canvas.setPayload(std::make_shared<dsl::Payload>(payload));toolbar.updateCoordinates();
            save(canvas.handle(),toolbar.handle(),output/L"viewer-dimensionless.png",1100,760);
            SendMessageW(physical,BM_CLICK,0,0);
            save(canvas.handle(),toolbar.handle(),output/L"viewer-dimensional.png",1100,760);
            if(argc>2) {
                canvas.setPayload(std::make_shared<dsl::Payload>(dsl::buildNativePayload(argv[2],"phase")));
                canvas.setDisplayCoordinates("dimensionless");toolbar.updateCoordinates();
                save(canvas.handle(),toolbar.handle(),output/L"viewer-phase-wide.png",1100,580);
            }
        }
        DestroyWindow(parent);parent=nullptr;Gdiplus::GdiplusShutdown(token);token=0;
        std::cout<<"plot presentation: PASS\n";return 0;
    } catch(const std::exception& error) {
        if(parent)DestroyWindow(parent);if(token)Gdiplus::GdiplusShutdown(token);
        std::cerr<<error.what()<<'\n';return 1;
    }
}
