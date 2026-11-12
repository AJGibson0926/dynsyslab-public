#include <dynsyslab/core/CoordinateTransform.hpp>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
std::vector<double> times(const TimeSeries& s) {
    return {s.timeData(),s.timeData()+s.sampleCount()};
}
std::vector<double> values(const TimeSeries& s) {
    return {s.valuesData(),s.valuesData()+s.sampleCount()*s.componentCount()};
}
}
void requireAlignedTimeSeries(const TimeSeries& a,const TimeSeries& b) {
    if(a.sampleCount()!=b.sampleCount()||a.componentCount()!=b.componentCount())
        throw std::invalid_argument("Coordinate values and rates must have identical shapes");
    for(std::size_t k=0;k<a.sampleCount();++k)
        if(a.time(k)!=b.time(k))throw std::invalid_argument("Coordinate values and rates have different clocks");
}
TimeSeries trimTimeSeries(const TimeSeries& s,std::size_t first,std::size_t count) {
    if(first>s.sampleCount()||count>s.sampleCount()-first||!count)
        throw std::invalid_argument("Invalid time series slice");
    return TimeSeries({s.timeData()+first,s.timeData()+first+count},
        {s.valuesData()+first*s.componentCount(),s.valuesData()+(first+count)*s.componentCount()},
        s.componentCount());
}
TimeSeries windowIntegralCoordinates(const TimeSeries& signal, double duration) {
    if (signal.componentCount() != 1 || !std::isfinite(duration) || duration <= 0)
        throw std::invalid_argument("Window integration requires a scalar signal and a positive finite duration");
    const auto count = signal.sampleCount();
    const auto* t = signal.timeData();
    const auto* y = signal.valuesData();
    if (t[count-1] - t[0] < duration)
        throw std::invalid_argument("No complete integral window in this record");
    for (std::size_t k = 0; k < count; ++k)
        if (!std::isfinite(y[k])) throw std::invalid_argument("Window integration requires finite samples");

    const auto segment = [&](std::size_t k) {
        return (t[k+1] - t[k]) * (0.5*y[k] + 0.5*y[k+1]);
    };
    // Compensate additions and removals as large pulses pass out of the window.
    double sum = 0, correction = 0;
    const auto add = [&](double area) {
        const double next = sum + area;
        correction += std::abs(sum) >= std::abs(area) ? (sum-next)+area : (area-next)+sum;
        sum = next;
    };
    std::size_t left = 0;
    std::vector<double> output;
    for (std::size_t end = 1; end < count; ++end) {
        if (end-1 > left) add(segment(end-1));
        if (t[end] - t[0] < duration) continue;
        // The accumulator covers only whole segments to the right of left+1.
        while (t[end] - t[left+1] >= duration) {
            add(-segment(left+1));
            ++left;
        }
        const double width = duration - (t[end] - t[left+1]);
        const double fraction = width / (t[left+1] - t[left]);
        const double atStart = fraction*y[left] + (1-fraction)*y[left+1];
        const double partial = width * (0.5*atStart + 0.5*y[left+1]);
        output.push_back(y[end]);
        output.push_back((sum + correction) + partial);
    }
    const auto first = count - output.size()/2;
    return TimeSeries({t+first, t+count}, std::move(output), 2);
}
TimeSeries differentiateTimeSeries(const TimeSeries& s) {
    if(s.sampleCount()<5||!s.uniformTimeSteps())
        throw std::invalid_argument("Fourth-order differentiation needs at least five uniform samples");
    const auto n=s.sampleCount(),d=s.componentCount();
    std::vector<double> rates((n-4)*d);
    for(std::size_t k=2;k+2<n;++k)for(std::size_t j=0;j<d;++j)
        rates[(k-2)*d+j]=(s.value(k-2,j)-8*s.value(k-1,j)+8*s.value(k+1,j)-s.value(k+2,j))/(12*s.deltaTime());
    return TimeSeries({s.timeData()+2,s.timeData()+n-2},std::move(rates),d);
}
CoordinateSeries transformCoordinates(const TimeSeries& s,const Observable& map,const TimeSeries* rates) {
    if(map.stateDimension()!=s.componentCount())
        throw std::invalid_argument("Coordinate map input dimension mismatch");
    auto x=values(s);
    CoordinateSeries result{TimeSeries(times(s),map.evaluateRealBatch(x,s.sampleCount()),map.outputDimension()),{}};
    if(rates) {
        requireAlignedTimeSeries(s,*rates);
        result.derivatives=TimeSeries(times(s),
            map.pushForwardRealBatch(x,values(*rates),s.sampleCount()),map.outputDimension());
    }
    return result;
}
CoordinateSeries transformHistoryCoordinates(const TimeSeries& s,const FrozenHistoryTransform& transform,
    const TimeSeries* rates) {
    if(s.componentCount()!=1||!s.uniformTimeSteps()
        ||std::abs(s.deltaTime()-transform.sampleInterval())>1e-8*transform.sampleInterval())
        throw std::invalid_argument("History input must be scalar and match the frozen transform clock");
    if(rates)requireAlignedTimeSeries(s,*rates);
    const auto& full=transform.basis();auto basis=full;
    const auto& selected=transform.retainedCoordinates();
    basis.rank=selected.size();basis.leftSingularVectors.resize(basis.delays*basis.rank);
    basis.singularValues.resize(basis.rank);basis.cumulativeEnergy.clear();
    for(std::size_t j=0;j<basis.rank;++j) {
        basis.singularValues[j]=full.singularValues[selected[j]];
        for(std::size_t k=0;k<basis.delays;++k)
            basis.leftSingularVectors[k*basis.rank+j]=full.leftSingularVectors[k*full.rank+selected[j]];
    }
    const auto q=projectHistory(values(s),basis);
    std::vector<double> t(s.timeData()+basis.delays-1,s.timeData()+s.sampleCount());
    CoordinateSeries result{TimeSeries(t,q.coordinates,basis.rank),{}};
    if(rates)result.derivatives=TimeSeries(t,projectHistory(values(*rates),basis).coordinates,basis.rank);
    return result;
}
}
