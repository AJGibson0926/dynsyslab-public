// Assemble standard Viewer figures; PlotCanvas owns all drawing/export behavior.
#include <dynsyslab_viewer/NativePayload.hpp>
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <iostream>
namespace c=dynsyslab::core;
using J=c::json::Value;
using Path=std::filesystem::path;
struct Case {
    std::string label,color;
    c::KoopmanModel model;
    c::NpzRecord solution;
};
dsl::Figure timeFigure(const std::vector<Case>& cases,int component,double end,const std::string& title){
    dsl::Figure f{"time",title,"\\tau",component==0?"r":component==1?"\\dot r":"p_{control}/P_0",false};
    f.valueScale="linear";f.dashboardPage=0;
    for(const auto& item:cases){
        const auto t=item.solution.at("time").doubles();
        const auto y=item.solution.at(component==2?"control_values":"values").doubles();
        dsl::Trace trace;trace.label=item.label;trace.color=item.color;trace.role="data";
        for(std::size_t k=0;k<t.size()&&t[k]<=end;++k){trace.x.push_back(t[k]);trace.y.push_back(component==2?y[k]:y[2*k+component]);}
        f.traces.push_back(std::move(trace));
    }
    if(component<2){dsl::Trace target;target.label="target";target.role="control_target";target.color="#333333";target.lineStyle="dash";
        target.x={0,end};target.y={component==0?1.3:0,component==0?1.3:0};f.traces.push_back(target);}
    return f;
}
dsl::Figure phaseFigure(const std::vector<Case>& cases){
    dsl::Figure f{"phase","Controlled phase portrait","r","\\dot r",false};f.valueScale="linear";f.dashboardPage=0;
    for(const auto& item:cases){const auto x=item.solution.at("values").doubles();dsl::Trace t;
        t.label=item.label;t.color=item.color;t.role="data";
        for(std::size_t k=0;k<x.size();k+=2){t.x.push_back(x[k]);t.y.push_back(x[k+1]);}f.traces.push_back(std::move(t));}
    f.references={{1,0,"IC",false},{1.3,0,"target",false}};return f;
}
dsl::Figure freeFigure(const Case& item,const c::NpzRecord& record,bool imaginary){
    const auto t=record.at("time").doubles(),x=record.at("values").doubles();const auto z0=item.model.evaluate(x.data())[0];
    dsl::Figure f{"time",item.label+(imaginary?": imaginary evolution":": real evolution"),"\\tau",
        imaginary?"Im(\\psi/\\psi_0)":"Re(\\psi/\\psi_0)",false};f.dashboardPage=0;f.valueScale="linear";
    dsl::Trace actual,ideal;actual.label="evaluated polynomial";actual.color=item.color;actual.role="data";
    ideal.label="exp(lambda tau)";ideal.color="#333333";ideal.lineStyle="dash";ideal.role="reference";
    for(std::size_t k=0;k<t.size();++k){const auto a=item.model.evaluate(&x[2*k])[0]/z0;
        const auto b=std::exp(item.model.continuousEigenvalues()[0]*(t[k]-t[0]));
        actual.x.push_back(t[k]);ideal.x.push_back(t[k]);actual.y.push_back(imaginary?a.imag():a.real());ideal.y.push_back(imaginary?b.imag():b.real());}
    f.traces={std::move(actual),std::move(ideal)};return f;
}
void render(const dsl::Payload& p,const Path& path,int w,int h){
    if(std::filesystem::exists(path))throw std::runtime_error("Refusing overwrite");
    dsl::PlotRenderOptions options;options.width=w;options.height=h;options.toolbar=false;
    std::cout<<c::json::dump(dsl::renderPlotPng(p,path,options))<<'\n';
}
int main(int argc,char** argv){try{
    if(argc!=2)throw std::runtime_error("PLOT_CONFIG_JSON");const auto config=c::json::load(Path(argv[1]).wstring());
    const Path out=config.at("output").string();std::filesystem::create_directories(out);
    std::vector<Case> cases;for(const auto& v:config.at("cases").array())cases.push_back({v.at("label").string(),v.at("color").string(),
        c::loadKoopmanModel(v.at("krom").string()),c::readNpz(v.at("solution").string())});
    if(cases.size()!=4)throw std::runtime_error("Expected supplied, two fresh modes and gauge case");
    const std::vector<Case> main{cases[0],cases[1],cases[2]};
    dsl::Payload control{argv[1],"comparison","Same standard KLQR, Q and R: supplied vs freshly trained eigenfunction",3};
    control.figures={phaseFigure(main),timeFigure(main,0,4,"Radius: first four cycles"),timeFigure(main,0,40,"Radius: complete run"),
        timeFigure(main,1,4,"Velocity: first four cycles"),timeFigure(main,2,4,"Signed control: first four cycles"),timeFigure(main,2,40,"Signed control: complete run")};
    control.figures[4].inset=dsl::DetailInset{.25,4,-2,.5,{0,1,2},false};
    control.figures[5].inset=dsl::DetailInset{5,40,-.73,-.68,{0,1,2},false};
    render(control,out/"control-comparison.png",2400,1400);
    const auto training=c::readNpz(config.at("training").string());
    dsl::Payload validation{argv[1],"validation","Unforced training trajectory: polynomial values vs ideal complex exponential",2};
    for(const auto& item:main){validation.figures.push_back(freeFigure(item,training,false));validation.figures.push_back(freeFigure(item,training,true));}
    render(validation,out/"free-validation-comparison.png",2200,1800);
    const std::vector<Case> gauge{cases[2],cases[3]};
    dsl::Payload sensitivity{argv[1],"gauge","Same fitted eigenfunction up to complex scale; unchanged numerical Q/R",2};
    sensitivity.figures={timeFigure(gauge,0,40,"Gauge sensitivity: radius"),timeFigure(gauge,2,4,"Gauge sensitivity: signed input")};
    render(sensitivity,out/"gauge-sensitivity.png",2000,750);
    render(dsl::buildNativePayload(config.at("training").string(),"data_overview"),out/"training-trajectory.png",2200,1400);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
