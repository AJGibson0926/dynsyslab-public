// Figure assembly only; all drawing uses the standard native Viewer renderer.
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
namespace c=dynsyslab::core;using Path=std::filesystem::path;
dsl::Figure figure(std::string title,std::string ylabel,std::string kind="time") {
    dsl::Figure f{kind,title,"\\tau",ylabel};f.dashboardPage=0;f.valueScale="linear";return f;
}
void line(dsl::Figure& f,const std::vector<double>& t,const std::vector<double>& v,
    std::size_t width,std::size_t col,std::string label,std::string color,
    std::string style,double begin,double end) {
    dsl::Trace trace;trace.label=label;trace.color=color;trace.lineStyle=style;
    std::size_t count=0;for(auto x:t)if(x>=begin&&x<=end)++count;
    const auto stride=std::max<std::size_t>(1,count/25000);std::size_t included=0;
    for(std::size_t k=0;k<t.size();++k)if(t[k]>=begin&&t[k]<=end) {
        if(included++%stride)continue;trace.x.push_back(t[k]);trace.y.push_back(v[k*width+col]);
    }
    f.traces.push_back(std::move(trace));
}
void render(const dsl::Payload& p,const Path& out,int w=2400,int h=1600) {
    dsl::PlotRenderOptions o;o.width=w;o.height=h;o.toolbar=false;
    if(std::filesystem::exists(out))throw std::runtime_error("Output already exists");
    std::cout<<c::json::dump(dsl::renderPlotPng(p,out,o))<<'\n';
}
std::vector<double> real(const std::vector<std::complex<double>>& z) {
    std::vector<double> x;for(auto v:z)x.push_back(v.real());return x;
}
void runPlots(const Path& run,const Path& subset,const Path& prefix) {
    const auto a=c::readNpz(run);const auto model=c::loadKoopmanModel(subset);
    const auto t=a.at("time").doubles(),x=a.at("true_states").doubles(),xh=a.at("values").doubles();
    const auto y=a.at("clean_echo").doubles(),yh=a.at("measurement_estimate").doubles();
    std::vector<double> tc,xr,yr;
    if(a.contains("comparison_time")) {
        tc=a.at("comparison_time").doubles();const auto psi=a.at("reference_psi").complexDoubles();
        xr=real(c::reconstructObservable(model.reconstruction("state"),psi,tc.size()));
        yr=real(c::reconstructObservable(model.reconstruction("far_field_acoustic_echo"),psi,tc.size()));
    }
    for(const std::string range:{"full","early"}) {
        dsl::Payload p{run.generic_u8string(),"klqe","Echo-only KLQE; retained original reconstruction columns",2};
        const double end=range=="full"?80:10;
        for(int row=0;row<3;++row)for(int column=0;column<2;++column) {
            const auto label=row==0?"r":row==1?"\\frac{dr}{d\\tau}":"e";
            auto f=figure(std::string(row==0?"Radius":row==1?"Radial velocity":"Echo")+
                (column==0?" - reference Psi reconstruction":" - echo-only KLQE"),label);
            line(f,t,row==2?y:x,row==2?1:2,row==1?1:0,"actual","#333333","dash",0,end);
            if(column==0&&!tc.empty())line(f,tc,row==2?yr:xr,row==2?1:2,row==1?1:0,
                "from reference Psi","#9864c4","solid",0,end);
            if(column==1)line(f,t,row==2?yh:xh,row==2?1:2,row==1?1:0,
                "KLQE estimate","#009b95","solid",0,end);
            p.figures.push_back(std::move(f));
        }
        render(p,Path(prefix.string()+"-"+range+".png"),2400,1900);
    }
    if(tc.empty())return;
    const auto reference=a.at("reference_psi").complexDoubles();
    const auto estimated=a.at("estimated_psi_real").doubles();
    dsl::Payload p{run.generic_u8string(),"krom_time","Reference modal coordinates versus KLQE estimates",2};
    std::size_t split=0;const auto m=model.eigenfunctionCount(),n=model.splitDimension();
    for(std::size_t j=0;j<std::min<std::size_t>(m,3);++j) {
        const bool complex=std::abs(model.continuousEigenvalues()[j].imag())>1e-12;
        for(int component=0;component<(complex?2:1);++component) {
            std::vector<double> v;for(std::size_t k=0;k<tc.size();++k)
                v.push_back(component?reference[k*m+j].imag():reference[k*m+j].real());
            auto f=figure("Mode "+std::to_string(j+1)+(component?" imaginary":" real"),component?"Im(\\psi)":"Re(\\psi)");
            line(f,tc,v,1,0,"reference","#333333","dash",0,80);
            line(f,tc,estimated,n,split+component,"KLQE estimate","#009b95","solid",0,80);
            p.figures.push_back(std::move(f));
        }split+=complex?2:1;
    }
    render(p,Path(prefix.string()+"-psi.png"),2400,std::max(700,650*int((p.figures.size()+1)/2)));
}
void diagnostics(const Path& inspection,const Path& gramian,const Path& prefix) {
    const auto a=c::readNpz(inspection),g=c::readNpz(gramian);
    dsl::Payload p{inspection.generic_u8string(),"decomposition","Joint Hankel fit and short-horizon modal observability",2};
    auto s=a.at("hankel_singular_values").doubles();double energy=0;for(double v:s)energy+=v*v;
    for(int panel=0;panel<2;++panel) {
        auto f=figure(panel?"Modal Gramian spectrum; horizon 0.001":"Joint Hankel energy spectrum","value","stem");
        f.xLabel=panel?"direction index":"singular value index";f.valueScale="log";
        dsl::Trace tr;tr.label=panel?"positive Gramian eigenvalues":"fraction of Hankel energy";tr.color="#007b95";tr.marker=true;
        const auto values=panel?g.at("gramian_eigenvalues").doubles():s;
        for(std::size_t j=0;j<values.size();++j) {
            const double v=panel?values[j]:values[j]*values[j]/energy;
            if(v>0){tr.x.push_back(double(j+1));tr.y.push_back(v);}
        }f.traces.push_back(std::move(tr));p.figures.push_back(std::move(f));
    }
    render(p,Path(prefix.string()+"-spectra.png"),2400,900);
    dsl::Payload fit{inspection.generic_u8string(),"time","Full-bank reconstruction using joint-window reference Psi",1};
    const auto t=a.at("time").doubles(),x=a.at("true_states").doubles(),xr=a.at("state_reconstruction").doubles();
    for(int row=0;row<3;++row) {
        auto f=figure(row==0?"Radius":row==1?"Velocity":"Echo",row==0?"r":row==1?"\\frac{dr}{d\\tau}":"e");
        line(f,t,row==2?a.at("echo").doubles():x,row==2?1:2,row==1?1:0,"actual","#333333","dash",0,80);
        line(f,t,row==2?a.at("echo_reconstruction").doubles():xr,row==2?1:2,row==1?1:0,
            "full-bank reference Psi","#9864c4","solid",0,80);fit.figures.push_back(std::move(f));
    }render(fit,Path(prefix.string()+"-reference.png"),2200,1700);
}
void directions(const Path& diagnostic,const Path& subset,const Path& prefix) {
    const auto d=c::readNpz(diagnostic);const auto model=c::loadKoopmanModel(subset);
    const auto n=model.splitDimension();const auto values=d.at("gramian_eigenvalues").doubles();
    dsl::Payload p{diagnostic.generic_u8string(),"observability","Selected subset: finite-window information geometry",2};
    auto spectrum=figure("Subset Gramian spectrum; T=2, dt=0.0001","Gramian eigenvalue","stem");
    spectrum.valueScale="log";spectrum.xLabel="direction index";
    dsl::Trace s;s.label="positive eigenvalues";s.color="#007b95";s.marker=true;
    for(std::size_t j=0;j<n;++j)if(values[j]>0){s.x.push_back(double(j+1));s.y.push_back(values[j]);}
    spectrum.traces.push_back(std::move(s));p.figures.push_back(std::move(spectrum));
    auto map=figure("Four strongest observable directions (unit vectors)","real modal coordinate index","observable_heatmap");
    map.xLabel="direction index";map.colorLabel="signed component";
    const auto columns=std::min<std::size_t>(4,n);
    dsl::Trace grid;grid.heatmap=true;grid.encoding="linear";grid.zShape={int(n),int(columns)};
    for(std::size_t j=0;j<columns;++j)grid.x.push_back(double(j+1));
    const auto vectors=d.at("directions").doubles();
    for(std::size_t k=0;k<n;++k){grid.y.push_back(double(k+1));for(std::size_t j=0;j<columns;++j)grid.z.push_back(vectors[k*n+j]);}
    map.traces.push_back(std::move(grid));p.figures.push_back(std::move(map));
    render(p,Path(prefix.string()+"-directions.png"),2400,1100);
}
int main(int argc,char** argv){try {
    if(argc!=5)throw std::runtime_error("run RUN SUBSET PREFIX | diagnostics INSPECTION GRAMIAN PREFIX");
    std::filesystem::create_directories(Path(argv[4]).parent_path());
    if(std::string(argv[1])=="run")runPlots(argv[2],argv[3],argv[4]);
    else if(std::string(argv[1])=="diagnostics")diagnostics(argv[2],argv[3],argv[4]);
    else if(std::string(argv[1])=="directions")directions(argv[2],argv[3],argv[4]);
    else throw std::runtime_error("unknown plot action");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
