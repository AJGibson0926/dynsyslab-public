#include <dynsyslab/core/Archive.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace c=dynsyslab::core;
namespace fs=std::filesystem;

// Score [maximum_i, maximum_(i+1)) windows, including a flagged first window
// when the receiver delay cuts off its beginning. No synthetic terminal maximum.
void ringErrors(const fs::path& reference,const fs::path& prediction,const fs::path& output) {
    const auto actual=c::readNpz(reference), fit=c::readNpz(prediction);
    const auto& shape=actual.at("state").shape;
    const auto t=actual.at("time").doubles(),bounds=actual.at("ring_bounds").doubles();
    if(shape.size()!=2 || shape[0]!=t.size() || !shape[1] || t.empty() || bounds.size()<2 ||
       fit.at("prediction").shape!=shape || fit.at("time").doubles()!=t)
        throw std::invalid_argument("Ring scoring requires aligned states, predictions and at least two maxima");
    auto increasing=[](const std::vector<double>& v) {
        return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);}) &&
            std::adjacent_find(v.begin(),v.end(),std::greater_equal<double>())==v.end();
    };
    if(!increasing(t)||!increasing(bounds))throw std::invalid_argument("Clocks and maxima must increase");
    const auto d=shape[1],n=bounds.size()-1;
    const auto center=actual.contains("error_center")?actual.at("error_center").doubles():std::vector<double>(d,0.);
    if(center.size()!=d)throw std::invalid_argument("Reference offsets must match state dimension");
    const auto x=actual.at("state").doubles(),y=fit.at("prediction").doubles();
    std::vector<double> se(n*d),ss(n*d),pct(n*d),counts(n),partial(n);
    for(std::size_t i=0;i<n;++i) {
        const auto first=std::lower_bound(t.begin(),t.end(),bounds[i])-t.begin();
        const auto end=std::lower_bound(t.begin(),t.end(),bounds[i+1])-t.begin();
        counts[i]=static_cast<double>(end-first);
        partial[i]=bounds[i]<t.front() || bounds[i+1]>t.back();
        for(auto k=first;k<end;++k)for(std::size_t j=0;j<d;++j) {
            const auto idx=k*d+j;
            const double error=y[idx]-x[idx], deviation=x[idx]-center[j];
            if(!std::isfinite(error)||!std::isfinite(deviation))throw std::invalid_argument("Nonfinite ring score");
            se[i*d+j]+=error*error;ss[i*d+j]+=deviation*deviation;
        }
        for(std::size_t j=0;j<d;++j) {
            const auto idx=i*d+j;
            pct[idx]=ss[idx]>0?100*std::sqrt(se[idx]/ss[idx]):std::numeric_limits<double>::quiet_NaN();
        }
    }
    c::NpzRecord out;out.arrays["ring_bounds"]=actual.at("ring_bounds");
    out.arrays["count"]=c::NpyArray::fromDoubles({n},counts);
    out.arrays["partial"]=c::NpyArray::fromDoubles({n},partial);
    out.arrays["error_center"]=c::NpyArray::fromDoubles({d},center);
    out.arrays["error_sum_squares"]=c::NpyArray::fromDoubles({n,d},se);
    out.arrays["reference_sum_squares"]=c::NpyArray::fromDoubles({n,d},ss);
    out.arrays["relative_rms_percent"]=c::NpyArray::fromDoubles({n,d},pct);
    c::writeNpzAtomicNew(output,out);
}
