// Assemble saved benchmark results through the standard native Viewer renderer.
// No fitting, simulation, or controller implementation is duplicated here.
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>

namespace c=dynsyslab::core;
using Path=std::filesystem::path;
c::json::Value readJson(const Path& path) {
    std::ifstream input(path);
    if(!input)throw std::runtime_error("Cannot read manifest");
    return c::json::parse(std::string(std::istreambuf_iterator<char>(input),{}));
}
dsl::Figure figure(const char* title,const char* y,const char* kind="time") {
    dsl::Figure f{kind,title,"t",y};f.dashboardPage=0;f.valueScale="linear";return f;
}
void line(dsl::Figure& f,const std::vector<double>& t,const std::vector<double>& v,
    std::size_t width,std::size_t col,const std::string& label,const std::string& color,
    const std::string& style,double end) {
    if(v.size()!=t.size()*width)throw std::runtime_error("Plot array shape mismatch");
    dsl::Trace trace;trace.label=label;trace.color=color;trace.lineStyle=style;
    for(std::size_t k=0;k<t.size();++k)if(t[k]<=end) {
        trace.x.push_back(t[k]);trace.y.push_back(v[k*width+col]);
    }
    f.traces.push_back(std::move(trace));
}
std::vector<double> norms(const std::vector<double>& x) {
    std::vector<double> result;
    for(std::size_t k=0;k<x.size();k+=2)result.push_back(std::hypot(x[k],x[k+1]));
    return result;
}
std::vector<double> errors(const std::vector<double>& t,const std::vector<double>& x,
    const std::vector<double>& te,const std::vector<double>& estimate) {
    std::vector<double> result;
    for(std::size_t k=0;k<t.size();++k) {
        auto match=std::lower_bound(te.begin(),te.end(),t[k]-1e-10);
        if(match==te.end()||std::abs(*match-t[k])>1e-9)
            throw std::runtime_error("No matching observer sample for truth timestamp");
        const auto j=static_cast<std::size_t>(match-te.begin());
        result.push_back(std::hypot(x[2*k]-estimate[2*j],x[2*k+1]-estimate[2*j+1]));
    }
    return result;
}
void phase(dsl::Figure& f,const std::vector<double>& time,const std::vector<double>& x,
    const std::string& label,const std::string& color,const std::string& style,double end) {
    dsl::Trace trace;trace.label=label;trace.color=color;trace.lineStyle=style;
    for(std::size_t k=0;k<time.size();++k)if(time[k]<=end) {
        trace.x.push_back(x[2*k]);trace.y.push_back(x[2*k+1]);
    }
    f.traces.push_back(std::move(trace));
}
void render(const std::vector<c::NpzRecord>& controlled,const std::vector<c::NpzRecord>& free,
    const Path& prefix,double end) {
    dsl::Payload p{"saved benchmark","feedback-benchmark",
        "Measured-position KLQE + KLQR: unstable linear Duffing benchmark",2};
    p.figures={
        figure("Position: measured y = x","x"),
        figure("Velocity: estimated from position only","v"),
        figure("Control input (enabled at t = 0.1)","u"),
        figure("State phase space","v","phase"),
        figure("State estimation error (both ICs coincide)","||x - xhat||"),
        figure("State norm: controlled vs uncontrolled","||x||")};
    p.figures[3].xLabel="x";p.figures[4].valueScale="log";p.figures[5].valueScale="log";
    c::json::Value::Array metrics;
    for(std::size_t j=0;j<controlled.size();++j) {
        const auto& a=controlled[j];const auto& b=free[j];
        const auto t=a.at("time").doubles(),x=a.at("values").doubles();
        const auto te=a.at("estimator_time").doubles(),xh=a.at("state_estimate").doubles();
        const auto command=a.at("commanded_input").doubles(),applied=a.at("applied_input").doubles();
        const auto e=errors(t,x,te,xh),n=norms(x);
        const std::string label=j?"IC -0.1":"IC +0.1",color=j?"#9755ad":"#009d99";
        for(std::size_t axis=0;axis<2;++axis) {
            line(p.figures[axis],t,x,2,axis,label+" actual",color,"solid",end);
            line(p.figures[axis],te,xh,2,axis,label+" estimate",j?"#573768":"#155361","dash",end);
        }
        line(p.figures[2],te,command,1,0,label+" requested",j?"#573768":"#155361","dash",end);
        line(p.figures[2],te,applied,1,0,label+" applied",color,"solid",end);
        phase(p.figures[3],t,x,label+" actual",color,"solid",end);
        phase(p.figures[3],te,xh,label+" estimate",j?"#573768":"#155361","dash",end);
        line(p.figures[4],t,e,1,0,label,color,j?"dash":"solid",end);
        line(p.figures[5],t,n,1,0,label+" controlled",color,j?"dash":"solid",end);
        if(j==0)line(p.figures[5],b.at("time").doubles(),norms(b.at("values").doubles()),
            1,0,"uncontrolled (both ICs)","#dc563c","dash",end);
        double peak=0;for(double u:applied)peak=std::max(peak,std::abs(u));
        metrics.push_back(c::json::Value::Object{
            {"initial_position",x[0]},{"final_state_norm",n.back()},
            {"final_estimation_error_norm",e.back()},{"maximum_applied_input",peak},
            {"uncontrolled_final_state_norm",norms(b.at("values").doubles()).back()}});
    }
    const auto suffix=end>2?"-full":"-transient";
    dsl::PlotRenderOptions options;options.width=2400;options.height=1800;options.toolbar=false;
    const Path output=prefix.string()+suffix+".png";
    if(std::filesystem::exists(output))throw std::runtime_error("Plot output already exists");
    std::cout<<c::json::dump(dsl::renderPlotPng(p,output,options))<<'\n';
    std::ofstream inspection(prefix.string()+suffix+".json");
    inspection<<c::json::dump(dsl::inspectPlot(p,true));
    if(end>2) {
        std::ofstream summary(prefix.string()+"-metrics.json");
        summary<<c::json::dump(c::json::Value::Object{{"runs",metrics}});
    }
}
int main(int argc,char** argv){try {
    if(argc!=4)throw std::runtime_error("Usage: CONTROLLED_RESULT.json FREE_RESULT.json OUTPUT_PREFIX");
    const auto controlledManifest=readJson(argv[1]),freeManifest=readJson(argv[2]);
    std::vector<c::NpzRecord> controlled,free;
    for(const auto& path:controlledManifest.at("archives").array())
        controlled.push_back(c::readNpz(path.string()));
    for(const auto& path:freeManifest.at("archives").array())
        free.push_back(c::readNpz(path.string()));
    if(controlled.size()!=2||free.size()!=2)throw std::runtime_error("Expected two ICs per result");
    const Path prefix=argv[3];std::filesystem::create_directories(prefix.parent_path());
    render(controlled,free,prefix,12);render(controlled,free,prefix,2);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
