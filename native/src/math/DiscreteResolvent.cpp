#include <dynsyslab/math/DiscreteResolvent.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
constexpr std::size_t maximumOrder=8;
void requireInterval(double h){if(!(h>0)||!std::isfinite(h))throw std::invalid_argument("Positive sample interval required");}
bool finite(C z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
double factorial(std::size_t n){double f=1;for(std::size_t k=2;k<=n;++k)f*=double(k);return f;}
// Stirling numbers of the second kind S(m,j) and signed first kind s(m,j).
struct Stirling {
    double second[maximumOrder][maximumOrder]{},first[maximumOrder][maximumOrder]{};
    Stirling(){second[0][0]=first[0][0]=1;
        for(std::size_t m=1;m<maximumOrder;++m)for(std::size_t j=0;j<=m;++j){
            second[m][j]=double(j)*second[m-1][j]+(j?second[m-1][j-1]:0);
            first[m][j]=(j?first[m-1][j-1]:0)-double(m-1)*first[m-1][j];}
    }
};
const Stirling& stirling(){static const Stirling s;return s;}
}
void validateDiscreteModel(const DiscretePoleModel& m){
    requireInterval(m.sampleInterval);if(m.poles.empty())throw std::invalid_argument("Discrete model needs at least one pole");
    for(const auto& p:m.poles){
        if(!finite(p.multiplier))throw std::invalid_argument("Multiplier must be finite");
        if(p.laurent.empty()||p.laurent.size()>maximumOrder)throw std::invalid_argument("Discrete pole order must be 1..8");
        for(const auto& a:p.laurent)if(!finite(a))throw std::invalid_argument("Laurent coefficients must be finite");
    }
}
C evaluateDiscreteResolvent(const DiscretePoleModel& m,C z){
    validateDiscreteModel(m);C sum{};
    for(const auto& p:m.poles){const C d=z-p.multiplier;if(d==C{})throw std::domain_error("z-domain resolvent evaluated at a pole");
        const C w=1./d;C power=z*w;for(const auto& a:p.laurent){sum+=a*power;power*=w;}}
    return sum;
}
C evaluateDiscreteSequence(const DiscretePoleModel& m,std::size_t k){
    validateDiscreteModel(m);C sum{};
    for(const auto& p:m.poles){double binomial=1;
        for(std::size_t q=0;q<p.laurent.size();++q){
            if(q>k)break;
            if(p.multiplier==C{}){if(k==q)sum+=p.laurent[q];}
            else sum+=p.laurent[q]*binomial*std::pow(p.multiplier,double(k-q));
            binomial*=double(k-q)/double(q+1);
        }
    }
    return sum;
}
DiscretePoleModel discretizePoleModel(const PoleResidueModel& model,double h){
    validatePoleModel(model);requireInterval(h);const auto& S=stirling();DiscretePoleModel out;out.sampleInterval=h;
    for(const auto& p:model.poles){
        if(p.laurent.size()>maximumOrder)throw std::invalid_argument("Discrete pole order must be 1..8");
        const C mu=std::exp(p.eigenvalue*h);ComplexVector d(p.laurent.size());
        for(std::size_t m=0;m<p.laurent.size();++m){
            const double scale=std::pow(h,double(m))/factorial(m);C power=1;
            for(std::size_t j=0;j<=m;++j){d[j]+=p.laurent[m]*scale*S.second[m][j]*factorial(j)*power;power*=mu;}
        }
        out.poles.push_back({mu,std::move(d)});
    }
    validateDiscreteModel(out);
    return out;
}
PoleResidueModel continuousPoleModel(const DiscretePoleModel& m,const std::vector<int>& branches){
    validateDiscreteModel(m);if(!branches.empty()&&branches.size()!=m.poles.size())throw std::invalid_argument("One branch per pole required");
    const auto& s=stirling();const double h=m.sampleInterval;PoleResidueModel out;
    for(std::size_t i=0;i<m.poles.size();++i){const auto& p=m.poles[i];
        if(p.multiplier==C{})throw std::domain_error("A zero multiplier has no continuous-time eigenvalue");
        const C lambda=C(std::log(std::abs(p.multiplier)),std::arg(p.multiplier)+2*pi*(branches.empty()?0:branches[i]))/h;
        ComplexVector a(p.laurent.size());
        for(std::size_t j=0;j<p.laurent.size();++j){const C scale=p.laurent[j]/(factorial(j)*std::pow(p.multiplier,double(j)));
            for(std::size_t q=0;q<=j;++q)a[q]+=scale*s.first[j][q]*factorial(q)/std::pow(h,double(q));}
        out.poles.push_back({lambda,std::move(a)});
    }
    return out;
}
AliasReport aliasReport(const PoleResidueModel& model,double h,double tolerance){
    validatePoleModel(model);requireInterval(h);if(!(tolerance>=0))throw std::invalid_argument("Tolerance must be nonnegative");
    AliasReport out;out.nyquist=pi/h;
    for(std::size_t i=0;i<model.poles.size();++i){const C lambda=model.poles[i].eigenvalue,mu=std::exp(lambda*h);
        const C principal=C(lambda.real(),std::arg(mu)/h);
        AliasEntry e;e.index=i;e.eigenvalue=lambda;e.multiplier=mu;e.principal=principal;
        e.turns=int(std::round((lambda.imag()-principal.imag())*h/(2*pi)));
        const double excess=std::abs(lambda.imag())-out.nyquist;
        e.onNyquist=std::abs(excess)<=1e-12*(1+out.nyquist);e.aliased=excess>0&&!e.onNyquist;
        out.entries.push_back(e);
    }
    for(std::size_t a=0;a<out.entries.size();++a)for(std::size_t b=a+1;b<out.entries.size();++b){
        const double d=std::abs(out.entries[a].multiplier-out.entries[b].multiplier);
        if(d<=tolerance&&std::abs(out.entries[a].eigenvalue-out.entries[b].eigenvalue)>tolerance)out.collisions.push_back({a,b,d});
    }
    return out;
}
}
