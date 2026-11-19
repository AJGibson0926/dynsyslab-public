#include "dynsyslab_viewer/TimeSeriesOverlay.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <algorithm>
#include <cmath>
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
        // Spectra: a 3-cycles-per-unit sine peaks at 3 x frequencyScale, with no zero bin.
        const auto wave=fs::temp_directory_path()/("viewer-overlay-spectrum-"+std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
        std::vector<double> clock,sine;
        for (int i=0;i<2000;++i) { clock.push_back(i*0.01); sine.push_back(2.0+std::sin(2*3.141592653589793*3*i*0.01)); }
        core::NpzRecord waves;
        waves.arrays["time"]=core::NpyArray::fromDoubles({2000},clock);
        waves.arrays["u"]=core::NpyArray::fromDoubles({2000},sine);
        core::writeNpzAtomicNew(wave,waves);
        dsl::TimeSeriesOverlay spectral; spectral.curves={{"u","drive"}}; spectral.spectrum=true; spectral.frequencyScale=10;
        const auto spectrum=dsl::timeSeriesOverlay(wave,spectral);
        const auto& line=spectrum.figures.front().traces.front();
        const auto peak=std::max_element(line.y.begin(),line.y.end())-line.y.begin();
        require(spectrum.figures.front().kind=="fourier" && line.x.front()>0 && std::abs(line.x[peak]-30)<0.6,
                "Spectrum peak frequency or scale");
        fs::remove(wave);
        // A run truncated on failure: NaN after the last computed sample is not drawn.
        const auto cut=fs::temp_directory_path()/("viewer-overlay-truncated-"+std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
        core::NpzRecord truncated;
        std::vector<double> clock10,partial;
        for (int i=0;i<10;++i) { clock10.push_back(i); partial.push_back(i<6?i:std::nan("")); }
        truncated.arrays["time"]=core::NpyArray::fromDoubles({10},clock10);
        truncated.arrays["y"]=core::NpyArray::fromDoubles({10},partial);
        core::writeNpzAtomicNew(cut,truncated);
        dsl::TimeSeriesOverlay cutReq; cutReq.curves={{"y","run"}};
        const auto cutLine=dsl::timeSeriesOverlay(cut,cutReq).figures.front().traces.front();
        require(cutLine.x.size()==6 && cutLine.x.back()==5.0,"Truncated series draws its computed part only");
        // Stem curve: one NaN-separated segment per finite sample from the baseline, plus tip markers.
        dsl::TimeSeriesOverlay stemReq; stemReq.curves={{"y","stems","#004DB8","stem"}};
        const auto stem=dsl::timeSeriesOverlay(cut,stemReq).figures.front();
        require(stem.traces.size()==3 && stem.traces[0].y==std::vector<double>({0.0,0.0}),"Zero baseline drawn under the stems");
        require(stem.traces[1].x.size()==18 && stem.traces[1].y[0]==0.0 &&
                stem.traces[1].y[4]==1.0 && std::isnan(stem.traces[1].y[2]),"Stems from a zero baseline");
        require(stem.traces[2].marker && stem.traces[2].markerShape=="open-circle" && stem.traces[2].x.size()==6 &&
                stem.traces[2].label.empty(),"Open-circle stem tips");
        stemReq.logBaseline=true;  // smallest positive sample is 1 -> baseline one decade below
        require(dsl::timeSeriesOverlay(cut,stemReq).figures.front().traces.size()==2 &&
                dsl::timeSeriesOverlay(cut,stemReq).figures.front().traces[0].y[3]==0.1,"Log-axis stem baseline (shared), no zero line");
        stemReq.stemBaseline=-2.0;
        require(dsl::timeSeriesOverlay(cut,stemReq).figures.front().traces[0].y[0]==-2.0,"Explicit stem baseline (log axis keeps no zero line)");
        fs::remove(cut);
        fs::remove(file);return 0;
    } catch (const std::exception& e) { fs::remove(file);std::cerr<<e.what()<<'\n';return 1; }
}
