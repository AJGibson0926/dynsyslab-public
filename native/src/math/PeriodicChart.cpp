#include <dynsyslab/math/PeriodicChart.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Tensor Hermite interpolation: de Boor, A Practical Guide to Splines (1978).
// Vector-field pullback is the ordinary chain rule under a coordinate chart;
// see Arnold, Mathematical Methods of Classical Mechanics, action-angle chapter.
namespace dynsyslab::math {
namespace {
std::pair<double,double> hermite(double f0,double f1,double d0,double d1,double t,double h){
    const double t2=t*t,t3=t2*t;
    return {(2*t3-3*t2+1)*f0+(-2*t3+3*t2)*f1+(t3-2*t2+t)*h*d0+(t3-t2)*h*d1,
        ((6*t2-6*t)*f0+(-6*t2+6*t)*f1)/h+(3*t2-4*t+1)*d0+(3*t2-2*t)*d1};
}
template<class F>double slope(const std::vector<double>& p,std::size_t i,F value){
    if(p.size()==2)return (value(1)-value(0))/(p[1]-p[0]);
    if(i==0){const double l=p[1]-p[0],r=p[2]-p[1],s=(value(1)-value(0))/l;return s-l*((value(2)-value(1))/r-s)/(l+r);}
    if(i+1==p.size()){const double l=p[i-1]-p[i-2],r=p[i]-p[i-1],s=(value(i)-value(i-1))/r;return s+r*(s-(value(i-1)-value(i-2))/l)/(l+r);}
    const double l=p[i]-p[i-1],r=p[i+1]-p[i];
    return (r*(value(i)-value(i-1))/l+l*(value(i+1)-value(i))/r)/(l+r);
}
}
PeriodicChart::PeriodicChart(const PeriodicAtlas& atlas,std::size_t samples)
    :samples_(samples),channels_(atlas.knots()[0].channels.size()){
    const auto& knots=atlas.knots();
    if(knots.size()<2||samples<4||samples>1000000||knots.size()>10000000/samples/channels_)
        throw std::invalid_argument("Periodic chart needs two knots and at most ten million table entries");
    table_.resize(knots.size()*samples_*channels_);
    for(std::size_t i=0;i<knots.size();++i){
        parameters_.push_back(knots[i].parameter);frequencies_.push_back(knots[i].channels[0].angularFrequency);
        for(std::size_t c=0;c<channels_;++c){
            const auto& original=knots[i].channels[c];auto real=original;
            for(const auto& [n,z]:original.coefficients)real.coefficients.try_emplace(-n,0.);
            auto coefficient=[&](int n){const auto it=original.coefficients.find(n);return it==original.coefficients.end()?std::complex<double>{}:it->second;};
            double correction=0,norm=0;
            for(auto& [n,z]:real.coefficients){z=.5*(coefficient(n)+std::conj(coefficient(-n)));correction+=std::abs(z-coefficient(n));norm+=std::abs(coefficient(n));}
            if(correction>1e-10*(1+norm))throw std::invalid_argument("Periodic chart requires Hermitian Fourier coefficients for real channels");
            realityProjection_=std::max(realityProjection_,correction);
            auto derivative=real;
            for(auto& [n,a]:derivative.coefficients)a*=std::complex<double>(0,double(n));
            const auto values=uniformFourierSamples(real,samples_);
            const auto phases=uniformFourierSamples(derivative,samples_);
            for(std::size_t j=0;j<samples_;++j){
                // Exact conjugate symmetry defines a real function. Remaining
                // imaginary inverse-FFT output is arithmetic roundoff; do not
                // amplify a tiny asymmetry by differentiating it first.
                if(!std::isfinite(values[j].real())||!std::isfinite(phases[j].real()))throw std::overflow_error("Nonfinite periodic chart samples");
                auto& e=table_[(i*channels_+c)*samples_+j];e.value=values[j].real();e.phase=phases[j].real();
            }
        }
    }
    for(std::size_t i=0;i<knots.size();++i){
        frequencySlopes_.push_back(slope(parameters_,i,[&](auto k){return frequencies_[k];}));
        for(std::size_t c=0;c<channels_;++c)for(std::size_t j=0;j<samples_;++j){
            auto& e=table_[(i*channels_+c)*samples_+j];
            e.parameter=slope(parameters_,i,[&](auto k){return table_[(k*channels_+c)*samples_+j].value;});
            e.mixed=slope(parameters_,i,[&](auto k){return table_[(k*channels_+c)*samples_+j].phase;});
        }
    }
}
std::pair<std::size_t,double> PeriodicChart::bracket(double p)const{
    if(!std::isfinite(p)||p<parameters_.front()||p>parameters_.back())throw std::out_of_range("Parameter outside periodic chart");
    auto hi=std::lower_bound(parameters_.begin(),parameters_.end(),p);
    const auto i=hi==parameters_.begin()?0:std::size_t(hi-parameters_.begin()-1);
    return {i,(p-parameters_[i])/(parameters_[i+1]-parameters_[i])};
}
std::pair<double,double> PeriodicChart::frequency(double p)const{
    const auto [i,a]=bracket(p);const auto result=hermite(frequencies_[i],frequencies_[i+1],frequencySlopes_[i],frequencySlopes_[i+1],a,parameters_[i+1]-parameters_[i]);
    if(!(result.first>0))throw std::domain_error("Interpolated periodic frequency must be positive");return result;
}
PeriodicChartValue PeriodicChart::evaluate(double p,double phase)const{
    if(!std::isfinite(phase))throw std::invalid_argument("Finite phase required");
    const auto [i,a]=bracket(p);const double twoPi=2*std::acos(-1.),h=twoPi/samples_,width=parameters_[i+1]-parameters_[i];
    phase=std::fmod(phase,twoPi);if(phase<0)phase+=twoPi;
    const double cell=phase/h;const auto j=std::min(samples_-1,std::size_t(cell)),k=(j+1)%samples_;const double b=cell-j;
    PeriodicChartValue out;
    for(std::size_t c=0;c<channels_;++c){
        std::array<std::pair<double,double>,2> v,dp;
        for(std::size_t row=0;row<2;++row){const auto base=((i+row)*channels_+c)*samples_;const auto& l=table_[base+j];const auto& r=table_[base+k];
            v[row]=hermite(l.value,r.value,l.phase,r.phase,b,h);
            dp[row]=hermite(l.parameter,r.parameter,l.mixed,r.mixed,b,h);
        }
        const auto va=hermite(v[0].first,v[1].first,dp[0].first,dp[1].first,a,width);
        out.value.push_back(va.first);out.parameterDerivative.push_back(va.second);
        out.phaseDerivative.push_back(hermite(v[0].second,v[1].second,dp[0].second,dp[1].second,a,width).first);
    }return out;
}
PlanarPullback pullbackPlanarField(const PeriodicChartValue& x,std::array<std::size_t,2> c,std::array<double,2> b,double minimum){
    if(c[0]==c[1]||c[0]>=x.value.size()||c[1]>=x.value.size()||x.parameterDerivative.size()!=x.value.size()||x.phaseDerivative.size()!=x.value.size()
        ||!std::isfinite(minimum)||minimum<=0||minimum>=1||!std::isfinite(b[0])||!std::isfinite(b[1]))throw std::invalid_argument("Invalid planar pullback arguments");
    const double a=x.parameterDerivative[c[0]],d=x.parameterDerivative[c[1]],h=x.phaseDerivative[c[0]],k=x.phaseDerivative[c[1]];
    const double np=std::hypot(a,d),nt=std::hypot(h,k),det=a*k-d*h;
    const double normalized=(a/np)*(k/nt)-(d/np)*(h/nt);
    if(!std::isfinite(normalized)||std::abs(normalized)<minimum)throw std::domain_error("Periodic chart Jacobian is singular or ill-conditioned");
    const double bp=((k/nt)*b[0]-(h/nt)*b[1])/normalized/np;
    const double bt=(-(d/np)*b[0]+(a/np)*b[1])/normalized/nt;
    if(!std::isfinite(bp)||!std::isfinite(bt)||!std::isfinite(det))throw std::overflow_error("Nonfinite periodic input map");
    return {{bp,bt},det,normalized};
}
}
