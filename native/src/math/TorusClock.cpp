#include <dynsyslab/math/TorusClock.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Finite cohomological straightening of a positive time change.
// Clark, Flows on solenoids are generically not almost periodic (1999),
// Theorem 1.1, https://arxiv.org/abs/math/9906150. Unlike an infinite
// solenoid expansion, a finite solvable polynomial has a bounded primitive.
// Numerical root inversion below does not certify an infinite limit or KAM.
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
void finite(double x) {if(!std::isfinite(x))throw std::overflow_error("Nonfinite clock calculation");}
std::vector<double> shift(const std::vector<double>& phase,const std::vector<double>& w,double s) {
    auto out=phase;for(std::size_t j=0;j<out.size();++j){out[j]+=w[j]*s;finite(out[j]);}return out;
}
C character(const FourierIndex& k,const std::vector<double>& phase) {
    auto reduced=phase;for(auto& p:reduced)p=std::remainder(p,2*std::acos(-1.));
    return std::exp(C(0,harmonicFrequency(k,reduced)));
}
}
TorusClockTransform analyzeTorusClock(const TorusFourierSeries& clock,double tolerance) {
    auto solved=solveCohomological(clock,tolerance); // validates the shared Fourier contract.
    const FourierIndex zero(clock.angularFrequencies.size(),0);
    const auto mean=clock.coefficients.find(zero);
    if(mean==clock.coefficients.end()||mean->second.imag()!=0||!(mean->second.real()>0))
        throw std::invalid_argument("Clock needs a positive real constant coefficient");
    double variation=0,primitiveNorm=0;
    for(const auto& [k,c]:clock.coefficients) {
        if(k==zero||c==C{})continue;
        auto opposite=k;for(auto& n:opposite)n=-n;
        const auto other=clock.coefficients.find(opposite);
        if(other==clock.coefficients.end()||other->second!=std::conj(c))
            throw std::invalid_argument("Clock Fourier coefficients must be Hermitian");
        variation+=std::abs(c);
    }
    for(const auto& [k,c]:solved.retained.coefficients)
        if(k!=zero&&c!=C{})throw std::domain_error("Nonconstant clock divisor retained; no global constant-frequency straightening returned");
    for(const auto& [k,c]:solved.solution.coefficients)primitiveNorm+=std::abs(c);
    const double meanRate=mean->second.real(),lower=meanRate-variation,upper=meanRate+variation;
    finite(lower);finite(upper);finite(primitiveNorm);
    if(!(lower>0))throw std::domain_error("Clock positivity not established by the sufficient Fourier l1 bound");
    auto frequencies=clock.angularFrequencies;
    for(auto& w:frequencies){w/=meanRate;finite(w);}
    return {clock,std::move(solved.solution),meanRate,lower,upper,primitiveNorm,std::move(frequencies)};
}
TorusClockEvolution evolveTorusClock(const TorusFourierSeries& clock,const std::vector<double>& initial,
    const std::vector<double>& times,const std::vector<FourierIndex>& characters,
    double divisorTolerance,double timeTolerance) {
    if(!std::isfinite(timeTolerance)||timeTolerance<=0||timeTolerance>1e-3)
        throw std::invalid_argument("Clock relative time tolerance must be in (0,1e-3]");
    TorusClockEvolution result{analyzeTorusClock(clock,divisorTolerance),{}};
    const auto& transform=result.transform;const auto& w=clock.angularFrequencies;
    const auto& h=transform.primitive;
    // Budget includes up to 128 root iterations and all returned characters.
    if(times.size()>1000000||characters.size()>100000||
       (!times.empty()&&(clock.coefficients.size()*256ULL+characters.size())>100000000ULL/times.size()))
        throw std::invalid_argument("Clock evolution exceeds work budget");
    const double initialH=evaluate(h,initial).real();
    for(const auto& k:characters) {
        harmonicFrequency(k,w);
        for(auto n:k)if(n < -1000000||n>1000000)throw std::invalid_argument("Clock character index too large");
    }
    for(double t:times) {
        finite(t);
        double s=0,residual=0;
        if(t!=0) {
            double lo=std::min(t/transform.lowerBound,t/transform.upperBound);
            double hi=std::max(t/transform.lowerBound,t/transform.upperBound);
            finite(lo);finite(hi);
            bool converged=false;
            const double allowed=timeTolerance*std::max(1.,std::abs(t));finite(allowed);
            for(unsigned iteration=0;iteration<128;++iteration) {
                s=lo+(hi-lo)/2; // brackets have the same sign, so subtraction is safe.
                residual=transform.meanRate*s+evaluate(h,shift(initial,w,s)).real()-initialH-t;
                finite(residual);
                if(std::abs(residual)<=allowed){converged=true;break;}
                if(s==lo||s==hi)break;
                if(residual<0)lo=s;else hi=s;
            }
            if(!converged)throw std::runtime_error("Clock inversion did not meet the requested time residual");
        }
        auto phase=shift(initial,w,s);
        auto straightened=shift(phase,w,evaluate(h,phase).real()/transform.meanRate);
        const double rate=evaluate(clock,phase).real();
        if(!(rate>0))throw std::runtime_error("Clock evaluation lost positivity numerically");
        TorusClockSample row{t,s,residual,rate,std::move(phase),std::move(straightened),{},{}};
        for(const auto& k:characters) {
            row.bareCharacters.push_back(character(k,row.phases));
            row.correctedCharacters.push_back(character(k,row.straightenedPhases));
        }
        result.samples.push_back(std::move(row));
    }
    return result;
}
}
