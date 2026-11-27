#include "dynsyslab_viewer/TimeSeriesOverlay.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/TimeSeries.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {
using Range = std::pair<std::size_t,std::size_t>;

// First, last and both extrema in each time bin maintain temporal order. At export
// each bin occupies at most a fraction of one horizontal pixel.
std::vector<Range> timeBins(const std::vector<double>& t,std::size_t count) {
    std::vector<Range> result;
    if (!count || count>=t.size()/4 || !std::isfinite(t.back()-t.front())) return result;
    std::size_t first=0;
    for (std::size_t k=1;k<count;++k) {
        const double edge=t.front()+(t.back()-t.front())*(static_cast<double>(k)/count);
        const auto last=static_cast<std::size_t>(std::lower_bound(t.begin()+first,t.end(),edge)-t.begin());
        if (last>first) result.emplace_back(first,last);
        first=last;
    }
    if (first<t.size()) result.emplace_back(first,t.size());
    return result;
}

Trace timeTrace(const TimeSeriesCurve& curve,const std::vector<double>& t,const std::vector<double>& data,
                std::size_t firstRow,std::size_t width,std::size_t component,const std::vector<Range>& bins) {
    Trace trace;
    // A label of "_" keeps the curve out of the legend.
    trace.label=curve.label=="_"?std::string{}:curve.label.empty()?curve.key:curve.label;
    trace.color=component<curve.componentColors.size()?curve.componentColors[component]:curve.color;
    trace.lineStyle=curve.lineStyle;
    if (curve.lineStyle=="markers" || curve.lineStyle=="open-markers") {
        trace.marker=true; trace.lineStyle="none";
        trace.markerShape=curve.lineStyle=="markers"?"large-circle":"large-open-circle";
    }
    trace.lineWidthPoints=curve.lineWidthPoints;
    trace.rightAxis=curve.rightAxis;
    const auto value=[&](std::size_t i){return data[(i+firstRow)*width+component];};
    // Nonfinite samples (e.g. after a run was truncated on failure) are not drawn.
    const auto append=[&](std::size_t i){
        if (!std::isfinite(value(i))) return;
        trace.x.push_back(t[i]); trace.y.push_back(value(i));
    };
    if (bins.empty()) {
        trace.x.reserve(t.size()); trace.y.reserve(t.size());
        for (std::size_t i=0;i<t.size();++i) append(i);
    } else {
        for (const auto& [first,last]:bins) {
            auto low=first,high=first;
            for (auto i=first;i<last;++i) {
                if (!std::isfinite(value(i))) continue;
                if (!std::isfinite(value(low)) || value(i)<value(low)) low=i;
                if (!std::isfinite(value(high)) || value(i)>value(high)) high=i;
            }
            std::array<std::size_t,4> keep{first,low,high,last-1};
            std::sort(keep.begin(),keep.end());
            const auto stop=std::unique(keep.begin(),keep.end());
            std::for_each(keep.begin(),stop,append);
        }
    }
    return trace;
}
}

// Stems: (t, base) -> (t, y) segments separated by NaN breaks, plus a tip-marker trace.
std::pair<Trace,Trace> stemTraces(const TimeSeriesCurve& curve,const std::vector<double>& t,const std::vector<double>& data,
                                  std::size_t firstRow,std::size_t width,std::size_t component,double base) {
    Trace stems=timeTrace(curve,{},{},0,1,0,{});
    stems.lineStyle="solid"; stems.x.clear(); stems.y.clear();
    if (!(curve.lineWidthPoints>0.)) stems.lineWidthPoints=0.75;  // thin stems, as in a classic stem plot
    Trace tips=stems; tips.label.clear(); tips.marker=true; tips.lineStyle="none"; tips.markerShape="open-circle";
    const auto value=[&](std::size_t i){return data[(i+firstRow)*width+component];};
    const double gap=std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i=0;i<t.size();++i) {
        const double y=value(i);
        if (!std::isfinite(y)) continue;
        stems.x.insert(stems.x.end(),{t[i],t[i],gap}); stems.y.insert(stems.y.end(),{base,y,gap});
        tips.x.push_back(t[i]); tips.y.push_back(y);
    }
    return {std::move(stems),std::move(tips)};
}

// Thin dark baseline under a stem panel, spanning the stems' abscissa range.
Trace stemBaseline(const std::vector<double>& t,double base) {
    Trace line; line.color="#3A3F4A"; line.lineStyle="solid"; line.lineWidthPoints=0.75;
    if (t.empty()) return line;
    const double pad=t.size()>1?0.5*(t.back()-t.front())/double(t.size()-1):0.5;
    line.x={t.front()-pad,t.back()+pad}; line.y={base,base};
    return line;
}

// Right-axis curves get their own ordinate: their range is mapped onto the plotted one, spanning
// the left curves' range (overlaid) or a band of the same height below it (stacked).
static void fitRightAxis(Figure& panel,const TimeSeriesOverlay& req) {
    if (std::none_of(panel.traces.begin(),panel.traces.end(),[](const Trace& t){return t.rightAxis;})) return;
    const auto range=[&](bool right) {
        double low=std::numeric_limits<double>::infinity(),high=-low;
        for (const auto& trace:panel.traces) if (trace.rightAxis==right)
            for (const double y:trace.y) if (std::isfinite(y)) { low=std::min(low,y); high=std::max(high,y); }
        if (high<=low) { const double mid=std::isfinite(low)?low:0.0; low=mid-0.5; high=mid+0.5; }
        return std::pair{low,high};
    };
    const auto [c,d]=range(true);
    const bool anyLeft=std::any_of(panel.traces.begin(),panel.traces.end(),[](const Trace& t){return !t.rightAxis;});
    const auto [a,b]=anyLeft?range(false):std::pair{c,d};
    const double height=b-a, top=req.stackAxes?a-0.12*height:b;  // the right band's upper edge
    Figure::RightAxis axis;
    axis.label=req.rightAxis?req.rightAxis->label:std::string{};
    axis.scale=height/(d-c); axis.offset=top-axis.scale*d;
    if (req.stackAxes) { axis.leftBand=std::pair{a,b}; axis.rightBand=std::pair{top-height,top}; }
    for (auto& trace:panel.traces) if (trace.rightAxis) for (auto& y:trace.y) y=axis.offset+axis.scale*y;
    panel.rightAxis=axis;
}

Payload timeSeriesOverlay(const std::filesystem::path& source,const TimeSeriesOverlay& req,std::size_t horizontalBins) {
    if (req.curves.empty()) throw std::invalid_argument("Select at least one curve");
    const auto record=dynsyslab::core::readNpz(source);
    const auto& clock=record.at(req.timeKey);
    if (clock.shape.size()!=1 || clock.shape[0]<2) throw std::invalid_argument("Expected a time vector with at least two samples");
    const auto times=clock.doubles();
    for (std::size_t i=0;i<times.size();++i)
        if (!std::isfinite(times[i]) || (i && times[i]<=times[i-1])) throw std::invalid_argument("Expected a finite, strictly increasing clock");
    if ((req.start && !std::isfinite(*req.start)) || (req.stop && !std::isfinite(*req.stop)))
        throw std::invalid_argument("Expected finite time limits");
    const auto first=std::lower_bound(times.begin(),times.end(),req.start.value_or(times.front()));
    const auto last=std::upper_bound(times.begin(),times.end(),req.stop.value_or(times.back()));
    if (last-first<2) throw std::invalid_argument("Time window must include at least two samples");
    const std::vector<double> t(first,last);
    const auto bins=timeBins(t,horizontalBins);
    Payload result;
    result.source=std::filesystem::absolute(source).u8string(); result.view="time"; result.tileName=req.title;
    // One stem baseline for the whole figure: given, 0, or (log axis) a decade below the
    // smallest positive sample of any stem curve in the window.
    double stemBase=req.stemBaseline.value_or(0.0);
    if (!req.stemBaseline && req.logBaseline) {
        double low=std::numeric_limits<double>::infinity();
        const auto row0=static_cast<std::size_t>(first-times.begin());
        for (const auto& curve:req.curves) {
            if (curve.lineStyle!="stem" || !record.contains(curve.key)) continue;
            const auto& a=record.at(curve.key);
            const std::size_t w=a.shape.size()==2?a.shape[1]:1;
            const auto v=a.doubles();
            for (std::size_t i=row0*w;i<(row0+t.size())*w && i<v.size();++i) if (std::isfinite(v[i]) && v[i]>0) low=std::min(low,v[i]);
        }
        stemBase=std::isfinite(low)?std::pow(10.0,std::floor(std::log10(low))-1.0):1.0;
    }
    std::set<std::size_t> baselined;  // panels that already carry a stem baseline
    for (const auto& curve:req.curves) {
        const auto comma=curve.key.find(',');
        const bool phase=comma!=std::string::npos;
        const std::string xKey=phase?curve.key.substr(0,comma):std::string{};
        const std::string yKey=phase?curve.key.substr(comma+1):curve.key;
        const auto& array=record.at(yKey);
        if (array.shape.empty() || array.shape.size()>2 || array.shape[0]!=times.size())
            throw std::invalid_argument("Expected an N or N-by-C series on the time grid: "+curve.key);
        const auto width=array.shape.size()==2?array.shape[1]:1;
        if (!width || (!req.componentLabels.empty() && req.componentLabels.size()!=width) ||
            (!result.figures.empty() && result.figures.size()!=width))
            throw std::invalid_argument("Component count mismatch: "+curve.key);
        const auto values=array.doubles();
        if (std::none_of(values.begin(),values.end(),[](double y){return std::isfinite(y);}))
            throw std::invalid_argument("Series has no finite samples: "+curve.key);
        if (result.figures.empty()) {
            for (std::size_t j=0;j<width;++j) {
                Figure panel;
                panel.kind="time"; panel.valueScale="linear"; panel.xLabel=req.timeLabel;
                panel.integerXTicks=req.integerTicks;
                panel.rightAxis=req.rightAxis;
                // Components are named only when there are several or they were labelled: a single
                // unlabelled series keeps the bare title and takes its ordinate from its curve.
                const bool named=!req.componentLabels.empty() || width>1;
                panel.yLabel=!req.componentLabels.empty()?req.componentLabels[j]:named?"Component "+std::to_string(j+1)
                    :curve.label.empty()?curve.key:curve.label;
                panel.title=req.title.empty()?panel.yLabel:named?req.title+" | "+panel.yLabel:req.title;
                if (!req.valueLabel.empty()) panel.yLabel=req.valueLabel;
                if (req.columns>=2) {
                    panel.gridSpan=1; panel.dashboardPage=static_cast<int>(j/(req.columns*req.columns));
                } else { panel.gridSpan=2; panel.dashboardPage=static_cast<int>(j/2); }
                result.figures.push_back(std::move(panel));
            }
        }
        std::vector<double> xValues;
        if (phase) {
            const auto& xArray=record.at(xKey);
            if (xArray.shape!=array.shape) throw std::invalid_argument("Phase abscissa shape differs from ordinate: "+xKey);
            xValues=xArray.doubles();
        }
        if (req.spectrum) {
            if (phase) throw std::invalid_argument("Spectra apply to time curves, not phase curves: "+curve.key);
            const auto row0=static_cast<std::size_t>(first-times.begin());
            for (std::size_t j=0;j<width;++j) {
                std::vector<double> column(t.size());
                for (std::size_t i=0;i<t.size();++i) column[i]=values[(i+row0)*width+j];
                if (!std::all_of(column.begin(),column.end(),[](double v){return std::isfinite(v);}))
                    throw std::invalid_argument("Spectra need finite samples in the window: "+curve.key);
                const dynsyslab::core::TimeSeries series(t,std::move(column),1);
                if (!series.uniformTimeSteps()) throw std::invalid_argument("Spectra need a uniform clock: "+curve.key);
                dynsyslab::core::FourierOptions options;
                options.oneSided=true; options.removeMean=true; options.window=dynsyslab::core::FourierWindow::hann;
                const auto transform=series.fourier().transform(options);
                std::vector<double> frequency,amplitude;
                for (std::size_t bin=1;bin<transform->binCount;++bin) {
                    const double f=transform->frequencies[bin]*req.frequencyScale;
                    if (req.frequencyRange && (f<req.frequencyRange->first || f>req.frequencyRange->second)) continue;
                    frequency.push_back(f);
                    amplitude.push_back(transform->amplitude(bin,0));
                }
                if (frequency.size()<2) throw std::invalid_argument("Spectrum has fewer than two bins in the frequency range: "+curve.key);
                auto& panel=result.figures[j];
                panel.kind="fourier";
                panel.traces.push_back(timeTrace(curve,frequency,amplitude,0,1,0,timeBins(frequency,horizontalBins)));
            }
            continue;
        }
        for (std::size_t j=0;j<width;++j) {
            if (!phase && curve.lineStyle=="stem") {
                auto [stems,tips]=stemTraces(curve,t,values,first-times.begin(),width,j,stemBase);
                auto& traces=result.figures[j].traces;
                if (!req.logBaseline && baselined.insert(j).second) traces.push_back(stemBaseline(t,stemBase));
                traces.push_back(std::move(stems)); traces.push_back(std::move(tips));
                continue;
            }
            if (!phase) {
                const bool points=curve.lineStyle=="markers" || curve.lineStyle=="open-markers";  // never binned
                result.figures[j].traces.push_back(timeTrace(curve,t,values,first-times.begin(),width,j,
                    points?std::vector<Range>{}:bins));
                continue;
            }
            // Phase curve: every vertex in the window (no time binning; the abscissa is not monotone).
            const auto row0=static_cast<std::size_t>(first-times.begin());
            Trace trace=timeTrace(curve,{},{},0,1,0,{});  // style only; points below
            trace.x.clear(); trace.y.clear();
            for (std::size_t i=0;i<t.size();++i) {
                const double x=xValues[(i+row0)*width+j], y=values[(i+row0)*width+j];
                if (std::isfinite(x) && std::isfinite(y)) { trace.x.push_back(x); trace.y.push_back(y); }
            }
            auto& panel=result.figures[j];
            panel.kind="phase"; panel.xLabel=req.xLabel.empty()?xKey:req.xLabel;
            panel.traces.push_back(std::move(trace));
        }
    }
    for (auto& panel:result.figures) fitRightAxis(panel,req);
    return result;
}
}
