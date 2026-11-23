#include <dynsyslab/math/Resonance.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
// Relative-phase kinematics; best-approximation interpretation: Short,
// arXiv:0912.1997. No inference of exact resonance from measured doubles.
namespace dynsyslab::math {
HarmonicResonance harmonicResonance(double a,double b,std::size_t m,std::size_t n,double t){
    if(!std::isfinite(a)||!std::isfinite(b)||a<=0||b<=0||!m||!n||m>1000000||n>1000000||!std::isfinite(t)||t<0)
        throw std::invalid_argument("Require positive finite frequencies, orders 1..1000000, and finite horizon>=0");
    const double delta=std::fma(double(n),b,-double(m)*a);
    HarmonicResonance r{m,n,delta,std::abs(delta)*t,std::abs(b/a-double(m)/double(n)),{}};
    if(!std::isfinite(delta)||!std::isfinite(r.relativePhaseDrift)||!std::isfinite(r.ratioError))
        throw std::overflow_error("Resonance arithmetic overflow");
    if(delta!=0){const double beat=2*std::acos(-1.)/std::abs(delta);if(std::isfinite(beat))r.beatTime=beat;}
    return r;
}
std::vector<HarmonicResonance> resonanceCandidates(double a,double b,std::size_t order,double t){
    if(!order||order>1000000)throw std::invalid_argument("Resonance search order must be 1..1000000");
    (void)harmonicResonance(a,b,1,1,t);
    std::vector<HarmonicResonance> result;double best=std::numeric_limits<double>::infinity();
    for(std::size_t n=1;n<=order;++n){
        const double target=double(n)*(b/a);
        if(n>1&&(!std::isfinite(target)||target>double(order)+.5))break;
        const auto m=std::size_t(std::clamp(std::round(target),1.,double(order)));
        if(m>order||std::gcd(m,n)!=1)continue;
        auto r=harmonicResonance(a,b,m,n,t);
        if(std::abs(r.detuning)<best){best=std::abs(r.detuning);result.push_back(r);}
    }
    return result;
}
} // namespace dynsyslab::math
