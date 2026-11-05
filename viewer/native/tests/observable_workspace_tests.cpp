#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include "../src/PlotSeries.hpp"
#include "../src/ObservablePayload.hpp"
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char* error) { if (!value) throw std::runtime_error(error); }
bool closeEnough(double a, double b) { return std::abs(a-b) < 1e-10 * std::max(1.,std::abs(b)); }
using dynsyslab::core::NpyArray;

dynsyslab::core::NpzRecord fixture(double offset = 0) {
    constexpr std::size_t count = 3001;
    std::vector<double> time(count), states(count*2), targets(count*2);
    for (std::size_t i=0;i<count;++i) {
        time[i]=i*.001; states[i*2]=-.9+offset+i*.0001; states[i*2+1]=.2;
        targets[i*2]=-.7; targets[i*2+1]=.1;
    }
    dynsyslab::core::NpzRecord record;
    record.arrays["time"]=NpyArray::fromDoubles({count},time);
    record.arrays["values"]=NpyArray::fromDoubles({count,2},states);
    record.arrays["control_target_values"]=NpyArray::fromDoubles({count,2},targets);
    record.arrays["metadata_json"]=NpyArray::fromString(R"({"system":"duffing","integration_form":"dimensionless","system_preset":"double_well_reference"})");
    return record;
}
void completeOutputs() {
    const auto record=fixture(); const auto metadata=dsl::json::parse(record.metadataJson());
    const auto payload=dsl::buildNativePayload("fixture.npz","observable_time",record,metadata);
    dynsyslab::core::SystemConfiguration config; config.identifier="duffing";config.preset="double_well_reference";
    const auto system=dynsyslab::core::createSystem(config); const auto catalog=system->observableCatalog();
    const auto states=record.at("values").doubles(); const auto target=record.at("control_target_values").doubles();
    std::size_t figure=0;
    for (const auto& name:catalog.names()) {
        const auto& observable=catalog.at(name);
        if (name == "state") {
            require(observable.outputDimension()==2 && observable.stateComponentIndex(0)==0 && observable.stateComponentIndex(1)==1,"fixture full state identity changed");
            continue;
        }
        for(std::size_t component=0;component<observable.outputDimension();++component)
            for(unsigned part=0;part<(observable.descriptor().complexValued?2u:1u);++part) {
                const auto& panel=payload.figures.at(figure);
                require(panel.dashboardPage==static_cast<int>(figure/4),"observable pages must contain four panels");
                require(panel.references.empty(),"observable time series has state reference points");
                require(panel.traces.size()==2 && panel.traces[0].x.size()==3001,"trajectory samples or target missing");
                require(panel.traces[1].role=="control_target" && panel.traces[1].lineStyle=="dash","target styling lost");
                for(const auto row:{0u,1500u,3000u}) {
                    const auto expected=observable.evaluate(states.data()+row*2)[component];
                    const auto expectedTarget=observable.evaluate(target.data()+row*2)[component];
                    require(closeEnough(panel.traces[0].y[row],part?expected.imag():expected.real()),"observable differs from public system evaluator");
                    require(closeEnough(panel.traces[1].y[row],part?expectedTarget.imag():expectedTarget.real()),"target observable differs from Main target");
                    require(closeEnough(panel.traces[0].x[row],row*.001),"observable time differs from Main");
                }
                ++figure;
            }
    }
    require(figure==payload.figures.size(),"catalog outputs were omitted or duplicated");
}
void higherDimensional() {
    dynsyslab::core::NpzRecord record;
    record.arrays["time"]=NpyArray::fromDoubles({2},{0,1});
    record.arrays["values"]=NpyArray::fromDoubles({2,3},{1,2,3,2,3,4});
    const auto payload=dsl::buildNativePayload("lorenz.npz","observable_workspace",record,
        dsl::json::parse(R"({"system":"lorenz","integration_form":"dimensional"})"));
    require(!payload.figures.empty(),"three-state observables missing");
    for (const auto& figure:payload.figures) require(figure.kind=="observable_time","three-state system received a 2D heat map");
}
void heatMaps() {
    const auto record=fixture(); const auto metadata=dsl::json::parse(record.metadataJson());
    const auto phase=dsl::buildNativePayload("fixture.npz","phase",record,metadata).figures.front();
    const auto limits=dsl::automaticAxisLimits(phase);
    const auto paired=dsl::buildNativePayload("fixture.npz","observable_workspace",record,metadata);
    dynsyslab::core::SystemConfiguration config;config.identifier="duffing";config.preset="double_well_reference";
    const auto catalog=dynsyslab::core::createSystem(config)->observableCatalog();
    std::size_t output=0;
    for (const auto& name:catalog.names()) {
        const auto& observable=catalog.at(name);
        if (name == "state") {
            require(observable.outputDimension()==2 && observable.stateComponentIndex(0)==0 && observable.stateComponentIndex(1)==1,"fixture full state identity changed");
            continue;
        }
        for(std::size_t component=0;component<observable.outputDimension();++component)
            for(unsigned part=0;part<(observable.descriptor().complexValued?2u:1u);++part) {
                const auto& time=paired.figures.at(output*2);const auto& map=paired.figures.at(output*2+1);
                require(time.kind=="observable_time" && map.kind=="observable_heatmap" && time.title==map.title,"observable pairs out of order");
                require(map.dashboardPage==time.dashboardPage,"observable pair split across pages");
                require(map.references.empty(),"heat map has IC or equilibrium markers");
                require(time.references.empty(),"observable time series has state reference points");
                const auto actual=dsl::automaticAxisLimits(map);const auto& trace=map.traces.front();
                require(closeEnough(actual.xMin,limits.xMin)&&closeEnough(actual.xMax,limits.xMax)
                    &&closeEnough(actual.yMin,limits.yMin)&&closeEnough(actual.yMax,limits.yMax),"heat-map region differs from Main phase plot");
                for(const auto row:{0u,37u,100u}) for(const auto col:{0u,63u,100u}) {
                    const double state[]{trace.x[col],trace.y[row]}; const auto expected=observable.evaluate(state)[component];
                    require(closeEnough(trace.z[row*101+col],part?expected.imag():expected.real()),"heat-map cell differs from public observable evaluator");
                }
                ++output;
            }
    }
    require(output*2==paired.figures.size(),"observable output missing a heat map");
    const dsl::PhaseSpaceDomain manual{{-1.7,.3,-.6,.9},"dimensionless"};
    const auto zoomed=dsl::buildObservablePayload("fixture.npz","observable_workspace",record,metadata,manual);
    const auto actual=dsl::automaticAxisLimits(zoomed.figures[1]);
    require(actual.xMin==manual.bounds.xMin && closeEnough(actual.xMax,manual.bounds.xMax)
        && actual.yMin==manual.bounds.yMin && closeEnough(actual.yMax,manual.bounds.yMax),"manual phase bounds were ignored");
}
void physicalArchive() {
    dynsyslab::core::NpzRecord record;
    record.arrays["time"]=NpyArray::fromDoubles({2},{0,3e-6});
    record.arrays["values"]=NpyArray::fromDoubles({2,2},{1.2e-5,1e-5/3e-6*.2,1.1e-5,0});
    const auto meta=dsl::json::parse(R"({"system":"rpe","integration_form":"dimensionless","archive_coordinate_system":"dimensional",
        "system_parameters":{"equilibrium_radius_m":1e-5,"preset":{"key":"air_water_20c","kappa":1.7,"sigma_n_m":0.05}},"nondimensionalization":{"time_scale_s":3e-6}})");
    const auto payload=dsl::buildNativePayload("legacy.npz","observable_time",record,meta);
    bool found=false;
    for(const auto& figure:payload.figures) if(figure.title=="Inertial kinetic energy") {
        require(closeEnough(figure.traces[0].y[0],1.5*.2*.2),"physical archive was evaluated as dimensionless states");
        require(closeEnough(figure.traces[0].x[1],1.),"observable time conversion differs from Main"); found=true;
    }
    require(found,"RPE observable missing");
    auto maps=dsl::buildObservablePayload("legacy.npz","observable_workspace",record,meta,
        dsl::PhaseSpaceDomain{{.9e-5,1.3e-5,-1.,1.},"dimensional"});
    bool checkedState=false;
    for(auto& map:maps.figures) if(map.kind=="observable_heatmap") {
        const auto& grid=map.traces.front();
        require(closeEnough(grid.x.front(),.9e-5)&&closeEnough(grid.y.front(),-1.),"dimensional phase grid converted incorrectly");
        if(map.colorLabel==map.xLabel) {
            require(closeEnough(grid.z[50*101+60],grid.x[60]),"state heat-map values use wrong units");
            const auto label=map.colorLabel;
            dsl::changeDisplayCoordinates(map,"dimensionless");
            require(map.colorLabel!="\\tau" && closeEnough(map.traces.front().z[50*101+60],1.14),"heat-map units did not switch with axes");
            dsl::changeDisplayCoordinates(map,"dimensional");
            require(map.colorLabel==label,"heat-map label lost on coordinate round trip"); checkedState=true;
        }
    }
    require(checkedState,"state heat map missing");
    dynsyslab::core::SystemConfiguration config;config.identifier="rpe";
    config.parameters={{"equilibrium_radius_m",1e-5},{"kappa",1.7},{"sigma_n_m",0.05}};
    const auto system=dynsyslab::core::createSystem(config);const double initial[]{1.2,.2};
    const auto expected=system->observableCatalog().at("laplace_pressure").evaluateReal(initial).front();
    for(const auto& figure:payload.figures) if(figure.title=="Laplace pressure")
        require(closeEnough(figure.traces[0].y[0],expected),"archived medium preset parameters were ignored");
}
void comparisons(const std::filesystem::path& root) {
    dynsyslab::core::writeNpzAtomic(root/"first.npz",fixture());
    dynsyslab::core::writeNpzAtomic(root/"second.npz",fixture(.1));
    dsl::BackendConfig config;config.projectRoot=root;config.archiveRoot=root/"archives";config.cache=root/"cache";
    dsl::Backend backend(config);dsl::CatalogItem a,b;
    a.id="a";a.kind="trajectory";a.source=(root/"first.npz").u8string();a.label="First";
    b=a;b.id="b";b.source=(root/"second.npz").u8string();b.label="Second";
    const auto main=backend.loadPayload(a,"time");const auto payload=backend.loadObservables({a,b});
    const auto single=dsl::buildNativePayload("fixture.npz","observable_workspace",fixture(),
        dsl::json::parse(fixture().metadataJson()));
    require(payload.figures.size()==single.figures.size(),"comparison maps are missing or duplicated");
    for(const auto& figure:payload.figures) {
        if (figure.kind == "observable_heatmap") { require(figure.traces.size()==1,"duplicate scalar fields overlaid"); continue; }
        require(figure.traces.size()==4,"comparison lost source or target traces");
        require(figure.traces[0].label=="First" && figure.traces[2].label=="Second","comparison source labels lost");
        require(figure.traces[0].x==main.figures[0].traces[0].x,"comparison observables changed time samples");
    }
}
void canvasAndMarkers() {
    dsl::theme::Fonts fonts; const auto parent=CreateWindowExW(0,L"STATIC",L"",WS_POPUP,0,0,900,700,
        nullptr,nullptr,GetModuleHandleW(nullptr),nullptr);
    dsl::PlotCanvas canvas;require(dsl::PlotCanvas::registerClass(GetModuleHandleW(nullptr))&&canvas.create(parent,fonts),"hidden canvas failed");
    auto first=std::make_shared<dsl::Payload>(dsl::buildNativePayload("first.npz","time",fixture(),dsl::json::parse(fixture().metadataJson())));
    canvas.setPayload(first);canvas.setAxisLimits(0,dsl::AxisLimits{0,3,-2,2});
    const auto main=canvas.captureState();canvas.setPayload(std::make_shared<dsl::Payload>());canvas.restoreState(main);
    require(canvas.captureState().payload==first && closeEnough(canvas.currentAxisPanels()[0].manual->yMin,-2),"Main state did not survive tab round trip");
    Gdiplus::Bitmap bitmap(200,200,PixelFormat32bppARGB);Gdiplus::Graphics g(&bitmap);g.Clear(dsl::theme::color(dsl::theme::plotBackground));
    dsl::Figure figure{"phase","","","",false};figure.references={{.5,1,"IC",false},{1.5,1,"equilibrium",true}};
    dsl::drawPlotReferences(g,figure,{0,0,200,200,0,2,0,2},fonts);
    Gdiplus::Color open,filled,rim;bitmap.GetPixel(50,100,&open);bitmap.GetPixel(150,100,&filled);bitmap.GetPixel(54,100,&rim);
    require(open.GetValue()==dsl::theme::color(dsl::theme::plotBackground).GetValue(),"IC dot is not open");
    require(filled.GetR()==255 && filled.GetG()==0 && filled.GetB()==0,"equilibrium dot is not solid red");
    require(rim.GetR()>200 && rim.GetG()<40,"IC outline is not red");
    for(const double fraction:{0.,.25,.5,.75,1.})
        require(dsl::trajectoryTimeColor(fraction).GetValue()==dsl::theme::heatColor(fraction).GetValue(),"Viewer color scales differ");
    require(dsl::theme::heatColor(0).GetValue()==Gdiplus::Color(255,0,0,127).GetValue()
        && dsl::theme::heatColor(.5).GetValue()==Gdiplus::Color(255,124,255,121).GetValue()
        && dsl::theme::heatColor(1).GetValue()==Gdiplus::Color(255,127,0,0).GetValue(),"jet colors differ from the standard palette");
    // Old exported payloads may still contain references; rendering enforces the current policy too.
    for (const auto* kind : {"time","control","observable_time","krom_time","klqe",
        "fourier","control_fourier","observable_heatmap","krom_heatmap"}) {
        g.Clear(dsl::theme::color(dsl::theme::plotBackground)); figure.kind=kind;
        figure.traces.resize(1); figure.traces.front().heatmap=figure.kind.find("heatmap")!=std::string::npos;
        dsl::drawPlotReferences(g,figure,{0,0,200,200,0,2,0,2},fonts);
        bitmap.GetPixel(150,100,&filled);
        require(filled.GetValue()==dsl::theme::color(dsl::theme::plotBackground).GetValue(),"equilibrium drawn outside a state phase plot");
        bitmap.GetPixel(54,100,&rim);
        require(rim.GetValue()==dsl::theme::color(dsl::theme::plotBackground).GetValue(),
            "IC drawn outside a state phase plot");
        // Hidden legacy references must not alter bounds or survive inspection/export.
        auto legacy=figure; legacy.references={{-1000,1000,"IC",false},{1000,-1000,"equilibrium",true}};
        legacy.traces.front().x={0,1}; legacy.traces.front().y={0,1};
        const auto actual=dsl::automaticAxisLimits(legacy);
        auto clean=legacy; clean.references.clear(); const auto expected=dsl::automaticAxisLimits(clean);
        require(closeEnough(actual.xMin,expected.xMin)&&closeEnough(actual.xMax,expected.xMax)
            &&closeEnough(actual.yMin,expected.yMin)&&closeEnough(actual.yMax,expected.yMax),
            "hidden state references changed plot bounds");
        dsl::Payload payload; payload.figures={legacy};
        for(const bool samples:{false,true})
            require(dsl::inspectPlot(payload,samples).at("figures").array().front().at("references").array().empty(),
                "inspection/export exposed hidden state references");
    }
    DestroyWindow(parent);
}
}
int wmain(int argc, wchar_t** argv) {
    if (argc == 3 && std::wstring(argv[1]) == L"--write-fixture") {
        dynsyslab::core::writeNpzAtomic(argv[2],fixture()); return 0;
    }
    ULONG_PTR token=0;Gdiplus::GdiplusStartupInput startup;Gdiplus::GdiplusStartup(&token,&startup,nullptr);
    const auto root=std::filesystem::temp_directory_path()/("dynsyslab-observable-workspace-"+std::to_string(GetCurrentProcessId()));
    try {std::filesystem::create_directories(root);completeOutputs();higherDimensional();heatMaps();physicalArchive();comparisons(root);canvasAndMarkers();
        std::filesystem::remove_all(root);Gdiplus::GdiplusShutdown(token);std::cout<<"observable workspace: PASS\n";return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';std::filesystem::remove_all(root);Gdiplus::GdiplusShutdown(token);return 1;}
}
