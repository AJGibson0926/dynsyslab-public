#include <dynsyslab/math/EigenvalueLattice.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
// Koopman eigenvalue semigroup: Mezic, Nonlinear Dyn. 41 (2005) 309-325;
// Budisic, Mohr & Mezic, Chaos 22 (2012) 047510. Discrete eigenvalues compare
// through log(mu), so products become sums and angles are taken modulo 2*pi.
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
constexpr double combinationLimit=5e6,workLimit=5e8;
double combinations(std::size_t generators,std::size_t order){
    double out=1; // C(generators+order,order)-1 nonzero multi-indices.
    for(std::size_t k=1;k<=order;++k)out*=double(generators+k)/double(k);
    return out-1;
}
struct Match {
    std::vector<unsigned> index;
    unsigned order{};
    C detuning;
    double size{};
    long long turns{};
    bool found{false};
};
class Search {
  public:
    Search(const std::vector<C>& generators,std::optional<double> period,std::size_t maximumOrder)
        :s_(generators),period_(period),maximum_(maximumOrder),current_(generators.size()) {}
    // Lowest order wins; within one order, the smallest detuning.
    Match best(C target,double tolerance){
        target_=target;tolerance_=tolerance;match_=Match{};
        if(s_.empty())return match_;
        for(unsigned order=1;order<=maximum_&&!match_.found;++order){order_=order;visit(0,order,C{});}
        return match_;
    }
  private:
    void visit(std::size_t j,unsigned remaining,C sum){
        if(j+1==s_.size()){current_[j]=remaining;consider(sum+double(remaining)*s_[j]);return;}
        for(unsigned m=0;m<=remaining;++m){current_[j]=m;visit(j+1,remaining-m,sum+double(m)*s_[j]);}
    }
    void consider(C sum){
        C d=target_-sum;long long turns=0;
        if(period_){
            const double k=std::round(d.imag()/ *period_);
            turns=(long long)k;d.imag(d.imag()-k* *period_);
        }
        const double size=std::abs(d);
        if(size<=tolerance_&&(!match_.found||size<match_.size))match_={current_,order_,d,size,turns,true};
    }
    const std::vector<C>& s_;
    std::optional<double> period_;
    std::size_t maximum_;
    std::vector<unsigned> current_;
    C target_;
    double tolerance_{};
    unsigned order_{};
    Match match_;
};
void validate(const EigenvalueLatticeOptions& o,std::size_t count){
    if(!o.maximumOrder||o.maximumOrder>64)throw std::invalid_argument("Lattice order must be 1..64");
    if(!o.maximumPrincipals||o.maximumPrincipals>32)throw std::invalid_argument("Principal count must be 1..32");
    if(!std::isfinite(o.tolerance)||o.tolerance<0)throw std::invalid_argument("Lattice tolerance must be finite and nonnegative");
    if(!std::isfinite(o.relativeTolerance)||o.relativeTolerance<0||o.relativeTolerance>=1)
        throw std::invalid_argument("Relative lattice tolerance must lie in [0,1)");
    if(!std::isfinite(o.observationTime)||o.observationTime<0)
        throw std::invalid_argument("Observation time must be finite and nonnegative");
    if(o.sampleInterval&&(!std::isfinite(*o.sampleInterval)||*o.sampleInterval<=0))
        throw std::invalid_argument("Sample interval must be positive and finite");
    if(o.clock==EigenvalueClock::Discrete&&o.sampleInterval)
        throw std::invalid_argument("Discrete eigenvalues are per sample; omit sample_interval");
    if(count>100000)throw std::invalid_argument("Eigenvalue list exceeds 100000 entries");
    const auto generators=o.principalIndices.empty()?o.maximumPrincipals:o.principalIndices.size();
    if(generators>32)throw std::invalid_argument("Principal count must be 1..32");
    const double total=combinations(generators,o.maximumOrder);
    if(total>combinationLimit||total*double(count)>workLimit)
        throw std::invalid_argument("Lattice search exceeds work limit; reduce order or principals");
}
}
const char* latticeRoleName(LatticeRole role){
    switch(role){
        case LatticeRole::Constant:return "constant";
        case LatticeRole::Principal:return "principal";
        case LatticeRole::Combination:return "combination";
        case LatticeRole::Unmatched:return "unmatched";
        default:return "excluded";
    }
}
EigenvalueLatticeReport eigenvalueLattice(const std::vector<C>& values,const EigenvalueLatticeOptions& o){
    validate(o,values.size());
    if(values.empty())throw std::invalid_argument("Eigenvalue list is empty");
    EigenvalueLatticeReport out;
    if(o.clock==EigenvalueClock::Discrete)out.aliasPeriod=2*pi;
    else if(o.sampleInterval)out.aliasPeriod=2*pi/ *o.sampleInterval;
    out.entries.resize(values.size());
    std::vector<std::size_t> candidates;
    for(std::size_t i=0;i<values.size();++i){
        auto& e=out.entries[i];e.index=i;e.eigenvalue=values[i];
        if(!std::isfinite(values[i].real())||!std::isfinite(values[i].imag()))
            throw std::invalid_argument("Eigenvalues must be finite");
        if(o.clock==EigenvalueClock::Discrete){
            if(values[i]==C{}){e.role=LatticeRole::Excluded;continue;} // log(0) is undefined.
            e.coordinate=std::log(values[i]);
        } else e.coordinate=values[i];
        e.acceptedTolerance=std::max(o.tolerance,o.relativeTolerance*std::abs(e.coordinate));
        if(std::abs(e.coordinate)<=o.tolerance){
            e.role=LatticeRole::Constant;e.detuning=e.coordinate;e.detuningMagnitude=std::abs(e.coordinate);
            e.phaseDrift=std::abs(e.coordinate.imag())*o.observationTime;continue;
        }
        candidates.push_back(i);
    }
    std::stable_sort(candidates.begin(),candidates.end(),[&](auto a,auto b){
        const double x=std::abs(out.entries[a].coordinate),y=std::abs(out.entries[b].coordinate);
        if(x!=y)return x<y;
        return out.entries[a].coordinate.imag()>out.entries[b].coordinate.imag();});
    std::vector<C> generators;
    const bool explicitGenerators=!o.principalIndices.empty();
    if(explicitGenerators)for(auto i:o.principalIndices){
        if(i>=values.size())throw std::invalid_argument("Principal index out of range");
        if(std::find(out.principal.begin(),out.principal.end(),i)!=out.principal.end())
            throw std::invalid_argument("Duplicate principal index");
        if(out.entries[i].role==LatticeRole::Excluded||out.entries[i].role==LatticeRole::Constant)
            throw std::invalid_argument("Principal index refers to a constant or excluded eigenvalue");
        out.principal.push_back(i);generators.push_back(out.entries[i].coordinate);
        out.entries[i].role=LatticeRole::Principal;
    }
    // Generators grow during a greedy pass, so a fresh Search is built per step.
    for(auto i:candidates){
        auto& e=out.entries[i];
        if(e.role==LatticeRole::Principal)continue;
        Search search(generators,out.aliasPeriod,o.maximumOrder);
        auto m=search.best(e.coordinate,e.acceptedTolerance);
        if(m.found){
            e.role=LatticeRole::Combination;e.multiIndex=std::move(m.index);e.order=m.order;
            e.detuning=m.detuning;e.detuningMagnitude=m.size;e.aliasTurns=m.turns;
            e.phaseDrift=std::abs(m.detuning.imag())*o.observationTime;
        } else if(!explicitGenerators&&generators.size()<o.maximumPrincipals){
            e.role=LatticeRole::Principal;out.principal.push_back(i);generators.push_back(e.coordinate);
        } else e.role=LatticeRole::Unmatched;
    }
    // Multi-indices found early are shorter than the final generator list.
    for(auto& e:out.entries){
        if(e.role==LatticeRole::Combination)e.multiIndex.resize(out.principal.size(),0);
        else if(e.role==LatticeRole::Principal){
            e.multiIndex.assign(out.principal.size(),0);e.order=1;
            const auto at=std::find(out.principal.begin(),out.principal.end(),e.index)-out.principal.begin();
            e.multiIndex[std::size_t(at)]=1;
        }
    }
    out.combinationCount=out.principal.empty()?0:combinations(out.principal.size(),o.maximumOrder);
    if(out.aliasPeriod)out.chanceCoverage=std::min(1.,out.combinationCount*2*o.tolerance/ *out.aliasPeriod);
    return out;
}
FrequencyResonanceReport frequencyResonanceReport(const std::vector<double>& frequencies,const FrequencyResonanceOptions& o){
    if(o.source!="hankel"&&o.source!="spectral_krom")
        throw std::invalid_argument("Frequency source must be hankel or spectral_krom");
    if(o.radiansPerSample&&o.sampleInterval)
        throw std::invalid_argument("Give radians_per_sample or sample_interval, not both");
    FrequencyResonanceReport out;out.source=o.source;
    std::vector<C> rates;rates.reserve(frequencies.size());
    for(std::size_t i=0;i<frequencies.size();++i){
        if(!std::isfinite(frequencies[i]))throw std::invalid_argument("Frequencies must be finite");
        if(frequencies[i]<0)out.foldedNegative.push_back(i);
        rates.emplace_back(0,std::abs(frequencies[i]));
    }
    EigenvalueLatticeOptions l;
    l.clock=EigenvalueClock::Continuous;
    l.sampleInterval=o.radiansPerSample?std::optional<double>(1.):o.sampleInterval;
    l.maximumOrder=o.maximumOrder;l.maximumPrincipals=o.maximumPrincipals;l.tolerance=o.tolerance;
    l.relativeTolerance=o.relativeTolerance;l.observationTime=o.observationTime;l.principalIndices=o.principalIndices;
    out.lattice=eigenvalueLattice(rates,l);
    const auto& p=out.lattice.principal;
    for(std::size_t a=0;a<p.size();++a)for(std::size_t b=a+1;b<p.size();++b){
        const double x=std::abs(frequencies[p[a]]),y=std::abs(frequencies[p[b]]);
        if(x<=0||y<=0)continue;
        for(const auto& r:resonanceCandidates(x,y,o.maximumOrder,o.observationTime))
            if(std::abs(r.detuning)<=o.tolerance)out.pairs.push_back({p[a],p[b],r});
    }
    return out;
}
}
