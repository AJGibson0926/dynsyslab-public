// Assemble already-computed archives with the standard Viewer renderer.
#include <dynsyslab_viewer/ViewerSession.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <algorithm>
#include <iostream>
namespace c=dynsyslab::core;
using Vector=std::vector<double>;
int main(int argc,char** argv) { try {
 if(argc!=7) throw std::runtime_error("state echo_jet history_coordinates basis end_time output.png");
 const auto state=c::readNpz(argv[1]), echo=c::readNpz(argv[2]), hist=c::readNpz(argv[3]), basis=c::readNpz(argv[4]);
 const double end=std::stod(argv[5]);
 const auto t=state.at("time").doubles(), x=state.at("values").doubles(), y=echo.at("values").doubles(),
            th=hist.at("time").doubles(), h=hist.at("values").doubles(), en=basis.at("cumulative_energy").doubles();
 dsl::Payload p{argv[1],"ringdown","Viscous RPE: corrected intrinsic monopole",2};
 auto panel=[&](std::string name,std::string xl,std::string yl,std::string kind,const Vector& time,
                const Vector& a,size_t aw,size_t ai,const Vector& b,size_t bw,size_t bi){
    dsl::Figure f{kind,name,xl,yl};f.dashboardPage=0;f.valueScale="linear";dsl::Trace tr;tr.color="#087c84";tr.label=yl;
    for(size_t k=0;k<time.size()&&time[k]<=end;++k){tr.x.push_back(a[aw*k+ai]);tr.y.push_back(b[bw*k+bi]);}
    f.traces.push_back(std::move(tr));p.figures.push_back(std::move(f));
 };
 panel("Radius","\\tau","r","time",t,t,1,0,x,2,0);
 panel("Physical state","r","dr/d\\tau","phase",t,x,2,0,x,2,1);
 panel("Distance-scaled radiated pressure","\\tau","e","time",t,t,1,0,y,2,0);
 panel("Echo and its autonomous derivative","e","de/d\\tau","phase",t,y,2,0,y,2,1);
 dsl::Figure f{"stem","Hankel energy per mode (whole record)","mode","fraction"};f.valueScale="log";f.dashboardPage=0;
 // Each mode has an independent stem down to a positive floor on the logarithmic axis.
 for(size_t j=0;j<en.size();++j){
     dsl::Trace stem;stem.color="#087c84";stem.x={double(j+1),double(j+1)};
     stem.y={0.01,en[j]-(j?en[j-1]:0)};f.traces.push_back(stem);
     dsl::Trace dot;dot.color=stem.color;dot.marker=true;dot.x={stem.x.back()};dot.y={stem.y.back()};f.traces.push_back(std::move(dot));
 }
 p.figures.push_back(std::move(f));
 panel("Causal echo history","\\xi_1","\\xi_2","phase",th,h,2,0,h,2,1);
 dsl::PlotRenderOptions options;options.width=2250;options.height=1650;options.toolbar=false;
 std::cout<<c::json::dump(dsl::renderPlotPng(p,argv[6],options));return 0;
 }catch(const std::exception& e){std::cerr<<e.what();return 1;}}
