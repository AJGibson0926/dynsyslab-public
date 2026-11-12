#include "dynsyslab_viewer/TimeSeriesOverlay.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>

namespace fs=std::filesystem;
namespace core=dynsyslab::core;
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
    const auto file=fs::temp_directory_path()/("viewer-overlay-"+std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    try {
        std::vector<double> t,truth,fit;
        for (int i=0;i<100;++i) {
            t.push_back(i*i*.01); truth.push_back(i); truth.push_back(i==18?-400.:i);
            fit.push_back(i+.5); fit.push_back(-i);
        }
        core::NpzRecord archive;
        archive.arrays["time"]=core::NpyArray::fromDoubles({100},t);
        archive.arrays["truth"]=core::NpyArray::fromDoubles({100,2},truth);
        archive.arrays["fit"]=core::NpyArray::fromDoubles({100,2},fit);
        core::writeNpzAtomicNew(file,archive);
        dsl::TimeSeriesOverlay req; req.curves={{"truth","Reference"},{"fit","Estimate"}};
        req.componentLabels={"Position","Velocity"};req.start=t[10];req.stop=t[60];
        const auto full=dsl::timeSeriesOverlay(file,req);
        require(full.figures.size()==2 && full.figures[1].traces.size()==2,"Overlay shape");
        const auto& trace=full.figures[1].traces[1];
        require(trace.label=="Estimate" && trace.x.size()==51 && trace.y.front()==-10 && trace.y.back()==-60,"Clock window or column mismatch");
        const auto compact=dsl::timeSeriesOverlay(file,req,3);
        const auto& draw=compact.figures[1].traces[0];
        require(draw.x.size()<=12 && draw.x.front()==t[10] && draw.x.back()==t[60] &&
                std::is_sorted(draw.x.begin(),draw.x.end()) && *std::min_element(draw.y.begin(),draw.y.end())==-400.,
                "Compressed export lost a peak, endpoint or time ordering");
        req.curves={{"time","Scalar"}};req.componentLabels={"Scalar"};
        require(dsl::timeSeriesOverlay(file,req).figures.size()==1,"Scalar series");
        req.curves.push_back({"fit","Estimate"});
        bool invalid=false;
        try { dsl::timeSeriesOverlay(file,req); } catch (const std::invalid_argument&) { invalid=true; }
        require(invalid,"Accepted inconsistent component counts");
        fs::remove(file);return 0;
    } catch (const std::exception& e) { fs::remove(file);std::cerr<<e.what()<<'\n';return 1; }
}
