#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Standard Fourier interpolation on a family of periodic orbits. For the
// Koopman action-angle interpretation see Mezic, Nonlinear Dynamics 41 (2005),
// 309-325, doi:10.1007/s11071-005-2824-x. No eigenfunction discovery is claimed.
namespace dynsyslab::math {
PeriodicAtlas::PeriodicAtlas(std::vector<PeriodicAtlasKnot> knots):knots_(std::move(knots)){
    if(knots_.empty()||knots_[0].channels.empty())throw std::invalid_argument("Nonempty orbit atlas required");
    for(std::size_t j=0;j<knots_.size();++j){const auto& k=knots_[j];
        if(!std::isfinite(k.parameter)||(j&&k.parameter<=knots_[j-1].parameter))throw std::invalid_argument("Atlas parameters must increase strictly");
        if(k.channels.size()!=knots_[0].channels.size())throw std::invalid_argument("Atlas channel count changed");
        for(const auto& s:k.channels){meanSquare(s);if(s.angularFrequency!=k.channels[0].angularFrequency)throw std::invalid_argument("Orbit channels must share one physical frequency");}
    }
}
std::vector<FourierSeries> PeriodicAtlas::interpolate(double parameter) const {
    if(!std::isfinite(parameter)||parameter<knots_.front().parameter||parameter>knots_.back().parameter)
        throw std::out_of_range("Orbit parameter outside fitted atlas; extrapolation is disabled");
    const auto hi=std::lower_bound(knots_.begin(),knots_.end(),parameter,[](const auto& k,double p){return k.parameter<p;});
    if(hi->parameter==parameter)return hi->channels;
    const auto& lo=*(hi-1);const double a=(parameter-lo.parameter)/(hi->parameter-lo.parameter);
    auto out=lo.channels;
    for(std::size_t j=0;j<out.size();++j){auto& s=out[j];
        s.angularFrequency=(1-a)*s.angularFrequency+a*hi->channels[j].angularFrequency;
        for(auto& item:s.coefficients)item.second*=1-a;
        for(const auto& [n,c]:hi->channels[j].coefficients)s.coefficients[n]+=a*c;
        meanSquare(s);
    }return out;
}
std::vector<std::complex<double>> PeriodicAtlas::evaluate(double parameter,double phase) const {
    std::vector<std::complex<double>> out;for(const auto& s:interpolate(parameter))out.push_back(math::evaluate(s,phase));return out;
}
std::vector<std::complex<double>> PeriodicAtlas::predict(double parameter,double phase,double elapsed) const {
    if(!std::isfinite(elapsed))throw std::invalid_argument("Finite elapsed time required");
    std::vector<std::complex<double>> out;for(const auto& s:interpolate(parameter))out.push_back(math::evaluate(s,phase+s.angularFrequency*elapsed));return out;
}
PeriodicAtlasTable::PeriodicAtlasTable(const PeriodicAtlas& atlas,std::size_t samples)
    :samples_(samples),channels_(atlas.knots()[0].channels.size()){
    if(samples<3||samples>1000000||atlas.knots().size()>10000000/samples/channels_)throw std::invalid_argument("Invalid periodic table size");
    for(const auto& k:atlas.knots()){
        parameters_.push_back(k.parameter);frequencies_.push_back(k.channels[0].angularFrequency);
        for(const auto& s:k.channels){const auto v=uniformFourierSamples(s,samples);values_.insert(values_.end(),v.begin(),v.end());}
    }
}
std::pair<std::size_t,double> PeriodicAtlasTable::bracket(double p) const {
    if(!std::isfinite(p)||p<parameters_.front()||p>parameters_.back())throw std::out_of_range("Parameter outside periodic table");
    const auto upper=std::lower_bound(parameters_.begin(),parameters_.end(),p);
    if(upper==parameters_.begin())return {0,0};
    const auto i=std::size_t(upper-parameters_.begin()-1);return {i,(p-parameters_[i])/(parameters_[i+1]-parameters_[i])};
}
double PeriodicAtlasTable::frequency(double p) const {
    const auto [i,a]=bracket(p);return a==0?frequencies_[i]:(1-a)*frequencies_[i]+a*frequencies_[i+1];
}
std::vector<std::complex<double>> PeriodicAtlasTable::evaluate(double p,double phase) const {
    if(!std::isfinite(phase))throw std::invalid_argument("Finite phase required");
    const auto [i,a]=bracket(p);const double twoPi=2*std::acos(-1.);
    phase=std::fmod(phase,twoPi);if(phase<0)phase+=twoPi;
    const double cell=phase/twoPi*samples_;const auto j=std::min(samples_-1,std::size_t(cell)),k=(j+1)%samples_;const double b=cell-j;
    std::vector<std::complex<double>> result;
    for(std::size_t c=0;c<channels_;++c){const auto base=(i*channels_+c)*samples_;
        auto value=(1-b)*values_[base+j]+b*values_[base+k];
        if(a!=0){const auto next=base+channels_*samples_;value=(1-a)*value+a*((1-b)*values_[next+j]+b*values_[next+k]);}
        result.push_back(value);
    }return result;
}
}
