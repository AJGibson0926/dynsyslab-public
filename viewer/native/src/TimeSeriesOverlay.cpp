#include "dynsyslab_viewer/TimeSeriesOverlay.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <algorithm>
#include <array>
#include <cmath>
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
    trace.label=curve.label.empty()?curve.key:curve.label;
    trace.color=curve.color; trace.lineStyle=curve.lineStyle;
    trace.lineWidthPoints=curve.lineWidthPoints;
    const auto value=[&](std::size_t i){return data[(i+firstRow)*width+component];};
    const auto append=[&](std::size_t i){trace.x.push_back(t[i]); trace.y.push_back(value(i));};
    if (bins.empty()) {
        trace.x.reserve(t.size()); trace.y.reserve(t.size());
        for (std::size_t i=0;i<t.size();++i) append(i);
    } else {
        for (const auto& [first,last]:bins) {
            auto low=first,high=first;
            for (auto i=first+1;i<last;++i) {
                if (value(i)<value(low)) low=i;
                if (value(i)>value(high)) high=i;
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
    for (const auto& curve:req.curves) {
        const auto& array=record.at(curve.key);
        if (array.shape.empty() || array.shape.size()>2 || array.shape[0]!=times.size())
            throw std::invalid_argument("Expected an N or N-by-C series on the time grid: "+curve.key);
        const auto width=array.shape.size()==2?array.shape[1]:1;
        if (!width || (!req.componentLabels.empty() && req.componentLabels.size()!=width) ||
            (!result.figures.empty() && result.figures.size()!=width))
            throw std::invalid_argument("Component count mismatch: "+curve.key);
        const auto values=array.doubles();
        if (!std::all_of(values.begin(),values.end(),[](double y){return std::isfinite(y);}))
            throw std::invalid_argument("Nonfinite series: "+curve.key);
        if (result.figures.empty()) {
            for (std::size_t j=0;j<width;++j) {
                Figure panel;
                panel.kind="time"; panel.valueScale="linear"; panel.xLabel=req.timeLabel;
                panel.yLabel=req.componentLabels.empty()?"Component "+std::to_string(j+1):req.componentLabels[j];
                panel.title=req.title.empty()?panel.yLabel:req.title+" | "+panel.yLabel;
                panel.gridSpan=2; panel.dashboardPage=static_cast<int>(j/2);
                result.figures.push_back(std::move(panel));
            }
        }
        for (std::size_t j=0;j<width;++j)
            result.figures[j].traces.push_back(timeTrace(curve,t,values,first-times.begin(),width,j,bins));
    }
    return result;
}
}
