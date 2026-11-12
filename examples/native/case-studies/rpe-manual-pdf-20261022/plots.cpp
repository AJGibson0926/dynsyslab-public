// Case-only assembly of standard Viewer figures. Never changes source archives.
#include <dynsyslab_viewer/NativePayload.hpp>
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <fstream>
#include <iostream>
namespace core=dynsyslab::core;
using Path=std::filesystem::path;
dsl::Figure component(const core::KoopmanModel& model,const core::NpzRecord& record,
                      bool imaginary,bool free) {
    const auto time=record.at("time").doubles(),x=record.at("values").doubles();
    const double target[]{1.3,0};const auto z0=model.evaluate(x.data())[0],zt=model.evaluate(target)[0];
    dsl::Figure fig{"time",imaginary?"Imaginary eigenfunction":"Real eigenfunction","\\tau",imaginary?"Im(\\psi)":"Re(\\psi)",false};
    fig.valueScale="linear";
    dsl::Trace actual,reference;actual.label=free?"evaluated polynomial":"controlled";
    actual.color="#202020";actual.role="data";
    reference.label=free?"exponential prediction":"target";reference.color="#ef4444";
    reference.lineStyle="dash";reference.role="control_target";
    actual.x=time;reference.x=time;
    for(std::size_t k=0;k<time.size();++k){
        const auto z=model.evaluate(&x[2*k])[0];
        const auto expected=free?z0*std::exp(model.continuousEigenvalues()[0]*(time[k]-time[0])):zt;
        actual.y.push_back(imaginary?z.imag():z.real());reference.y.push_back(imaginary?expected.imag():expected.real());
    }
    fig.traces={std::move(actual),std::move(reference)};return fig;
}
void render(const dsl::Payload& payload,const Path& path,int width,int height) {
    if(std::filesystem::exists(path))throw std::runtime_error("Refusing to overwrite plot");
    dsl::PlotRenderOptions options;options.width=width;options.height=height;options.toolbar=false;
    std::cout<<core::json::dump(dsl::renderPlotPng(payload,path,options))<<'\n';
}
int main(int argc,char** argv){try{
    if(argc!=5)throw std::runtime_error("KROM CONTROL_ARCHIVE FREE_ARCHIVE OUTPUT_DIRECTORY");
    const auto model=core::loadKoopmanModel(argv[1]);const auto controlled=core::readNpz(argv[2]),free=core::readNpz(argv[3]);
    const Path out=argv[4];std::filesystem::create_directories(out);
    auto native=dsl::buildNativePayload(argv[2],"data_overview");
    if(native.figures.size()!=4)throw std::runtime_error("Expected standard RPE control overview");
    dsl::Payload comparison{argv[2],"manual_rpe_comparison","RPE manual eigenfunction: PDF page 2",3};
    comparison.figures={native.figures[3],native.figures[1],native.figures[2],
        component(model,controlled,false,false),component(model,controlled,true,false),native.figures[0]};
    for(auto& figure:comparison.figures){figure.valueScale="linear";figure.dashboardPage=0;
        for(auto& trace:figure.traces)if(trace.role=="data"){trace.color="#202020";trace.encoding.clear();trace.colorValues.clear();}}
    auto& phase=comparison.figures[0];phase.references.push_back({1.3,0,"target",false});
    comparison.figures[5].yLabel="p_{control}/P_0";
    render(comparison,out/"control.png",2400,1400);
    dsl::Payload validation{argv[3],"manual_rpe_validation","Unforced validation of supplied eigenfunction",2};
    validation.figures={component(model,free,false,true),component(model,free,true,true)};
    render(validation,out/"unforced-validation.png",2000,750);
    // Separate magnitude view permits a log comparison without hiding the sign
    // in the primary control panel. The signed archive remains authoritative.
    auto magnitude=native.figures[0];magnitude.title="Control magnitude (PDF log-axis comparison)";
    magnitude.yLabel="|p_{control}/P_0|";magnitude.valueScale="log";magnitude.inset.reset();
    for(auto& trace:magnitude.traces){trace.label="|u|";for(auto& y:trace.y)y=std::abs(y);}
    dsl::Payload input{argv[2],"manual_rpe_input_magnitude","Magnitude only; see main plot for sign",1};input.figures={magnitude};
    render(input,out/"input-magnitude.png",1600,650);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
