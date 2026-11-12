// Existing KROM workspace generates the grids and traces; this client only
// selects representative overlays and assembles native Viewer figure panels.
#include <dynsyslab_viewer/KromWorkspace.hpp>
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
using Path=std::filesystem::path;namespace c=dynsyslab::core;
void representative(dsl::Figure& f){
    const bool map=!f.traces.empty()&&f.traces.front().heatmap;
    const std::vector<std::size_t> keep{4,14,24,34,44,54};
    const char* labels[]{"A=0.003","A=0.01","A=0.03","A=0.08","A=0.2","A=0.4"};
    std::vector<dsl::Trace> traces;if(map)traces.push_back(f.traces.front());
    for(std::size_t i=0;i<keep.size();++i){auto t=f.traces.at(keep[i]+(map?1:0));t.label=labels[i];
        t.lineStyle="dash";t.opacity=.64;t.color="#202020";traces.push_back(std::move(t));}
    f.traces=std::move(traces);f.dashboardPage=0;f.titleColor="#263746";
    f.valueScale="linear";f.alternateCoordinates.reset();
}
void render(const dsl::Payload& p,const Path& out,int w,int h){
    if(std::filesystem::exists(out))throw std::runtime_error("Refusing overwrite");
    dsl::PlotRenderOptions o;o.width=w;o.height=h;o.toolbar=false;
    std::cout<<c::json::dump(dsl::renderPlotPng(p,out,o))<<'\n';
}
void centralMap(dsl::Figure& f,const c::KoopmanModel& model,int component){
    auto& grid=f.traces.front();grid.x.clear();grid.y.clear();grid.z.clear();grid.magnitudes.clear();
    constexpr int n=181;for(int i=0;i<n;++i){grid.x.push_back(.65+(2.1-.65)*i/(n-1));grid.y.push_back(-4.+8.*i/(n-1));}
    grid.zShape={n,n};double maxMagnitude=0;
    for(double v:grid.y)for(double r:grid.x){const double x[]{r,v};const auto z=model.evaluate(x)[0];
        grid.z.push_back(component==0?z.real():component==1?z.imag():std::arg(z));
        if(component==2){grid.magnitudes.push_back(std::abs(z));maxMagnitude=std::max(maxMagnitude,std::abs(z));}}
    if(component==2)for(auto& m:grid.magnitudes)m/=maxMagnitude;
    const auto nan=std::numeric_limits<double>::quiet_NaN();
    for(std::size_t j=1;j<f.traces.size();++j)for(std::size_t k=0;k<f.traces[j].x.size();++k)
        if(f.traces[j].x[k]<.65||f.traces[j].x[k]>2.1||std::abs(f.traces[j].y[k])>4){f.traces[j].x[k]=nan;f.traces[j].y[k]=nan;}
    f.rowSpan=1;f.square=true;
}
int main(int argc,char** argv){try{
    if(argc!=3)throw std::runtime_error("SELECTED_POLYNOMIAL OUTPUT_DIRECTORY");
    const Path out=argv[2];std::filesystem::create_directories(out);
    dsl::KromWorkspaceOptions options;options.modes=std::vector<std::size_t>{0};
    auto real=dsl::buildKromWorkspace(argv[1],options);
    options.complexDisplay="domain-coloring";auto domain=dsl::buildKromWorkspace(argv[1],options);
    if(real.figures.size()!=5||domain.figures.size()!=4)throw std::runtime_error("Unexpected KROM layout");
    for(std::size_t i=1;i<real.figures.size();++i)representative(real.figures[i]);
    for(std::size_t i=1;i<domain.figures.size();++i)representative(domain.figures[i]);
    dsl::Payload phase{argv[1],"forced_edmdc_phase","Forced-only EDMDc: fourth-order polynomial in (r, dr/dtau)",3};
    phase.figures={real.figures[1],real.figures[3],domain.figures[1]};
    for(auto& f:phase.figures){f.rowSpan=1;f.square=true;f.traces.resize(1);}
    // Leave the first plot clean to expose the polynomial geometry. Overlay
    // sampled paths on the other two maps to show where observations exist.
    phase.figures[1].traces=real.figures[3].traces;phase.figures[2].traces=domain.figures[1].traces;
    render(phase,out/"phase-space.png",2700,1000);
    const auto model=c::loadKoopmanModel(argv[1]);
    for(std::size_t j=0;j<phase.figures.size();++j)centralMap(phase.figures[j],model,int(j));
    dsl::Figure coverage{"phase","Representative driven trajectories (6 of 60)","r","\\frac{dr}{d\\tau}",true};
    coverage.dashboardPage=0;coverage.valueScale="linear";
    coverage.traces.assign(phase.figures[1].traces.begin()+1,phase.figures[1].traces.end());
    coverage.references={{1,0,"EC",true}};
    phase.figures.push_back(coverage);phase.columns=2;
    render(phase,out/"phase-space-central.png",2400,1700);
    real.figures.erase(real.figures.begin());real.layout="";real.columns=2;
    render(real,out/"real-imaginary-history.png",2200,1400);
    domain.figures.erase(domain.figures.begin());domain.layout="";domain.columns=2;
    render(domain,out/"domain-color-history.png",2200,1400);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
