#include <dynsyslab/math/DiscreteCycles.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <dynsyslab/core/Hankel.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <map>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double difference(double a,double b,double period) {return period>0?std::remainder(a-b,period):a-b;}
void validate(const std::vector<double>& x,std::size_t d) {
    if(!d||d>1024||x.empty()||x.size()%d||x.size()>2000000)
        throw std::invalid_argument("Invalid bounded sample-by-coordinate array");
    for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("Finite samples required");
}
double num(const J& x) {
    if(!std::holds_alternative<double>(x.storage())||!std::isfinite(x.number()))throw std::invalid_argument("Finite number required");
    return x.number();
}
std::size_t integer(const J& x) {double v=num(x);if(v<0||v>2000000||v!=std::floor(v))throw std::invalid_argument("Bounded integer required");return std::size_t(v);}
std::vector<double> reals(const J& x) {
    if(!std::holds_alternative<A>(x.storage()))throw std::invalid_argument("Array required");
    std::vector<double> out;for(const auto& v:x.array())out.push_back(num(v));return out;
}
J encode(const std::vector<double>& x) {A out;for(double v:x)out.emplace_back(v);return out;}
}
SampledCycleResult diagnoseSampledCycles(const std::vector<double>& x,std::size_t d,
    const std::vector<std::size_t>& candidates,const std::vector<double>& scales,
    const std::vector<double>& periods,double tolerance,std::size_t repeats) {
    validate(x,d);const auto n=x.size()/d;
    if(scales.size()!=d||periods.size()!=d||repeats<2||!std::isfinite(tolerance)||tolerance<0)
        throw std::invalid_argument("Invalid cycle diagnostic settings");
    for(std::size_t c=0;c<d;++c)if(!std::isfinite(scales[c])||scales[c]<=0||!std::isfinite(periods[c])||periods[c]<0)
        throw std::invalid_argument("Positive scales and nonnegative coordinate periods required");
    auto sorted=candidates;std::sort(sorted.begin(),sorted.end());sorted.erase(std::unique(sorted.begin(),sorted.end()),sorted.end());
    SampledCycleResult result;
    for(auto p:sorted) {
        if(!p||p>n/repeats)throw std::invalid_argument("Candidate requires more complete repeats");
        CycleDefect defect{p,n-p,0,0};double sum=0;
        for(std::size_t j=p;j<n;++j) {
            double squared=0;
            for(std::size_t c=0;c<d;++c) {double e=difference(x[j*d+c],x[(j-p)*d+c],periods[c])/scales[c];squared+=e*e;}
            sum+=squared;defect.maximum=std::max(defect.maximum,std::sqrt(squared));
        }
        defect.rms=std::sqrt(sum/double(n-p));result.defects.push_back(defect);
        if(!result.detectedPeriod&&defect.maximum<=tolerance)result.detectedPeriod=p;
    }
    if(!result.detectedPeriod)return result;
    const auto p=result.detectedPeriod;result.foldedValues.assign(p*d,0);std::vector<std::size_t> count(p);
    for(std::size_t j=0;j<n;++j) {
        ++count[j%p];for(std::size_t c=0;c<d;++c) {
            const double anchor=x[(j%p)*d+c];
            result.foldedValues[(j%p)*d+c]+=anchor+difference(x[j*d+c],anchor,periods[c]);
        }
    }
    for(std::size_t j=0;j<p;++j)for(std::size_t c=0;c<d;++c)result.foldedValues[j*d+c]/=double(count[j]);
    double sum=0;for(std::size_t j=0;j<n;++j)for(std::size_t c=0;c<d;++c) {
        double e=difference(x[j*d+c],result.foldedValues[(j%p)*d+c],periods[c])/scales[c];sum+=e*e;
    }
    result.foldingRms=std::sqrt(sum/double(n));return result;
}
J discreteCycleSpectrum(const std::vector<double>& cycle,std::size_t d) {
    validate(cycle,d);const auto n=cycle.size()/d;const double pi=std::acos(-1.);
    A coefficients,masses,angles,denominators;double total=0;std::map<std::size_t,double> conductorMass;
    std::vector<std::vector<std::complex<double>>> transforms;
    for(std::size_t c=0;c<d;++c) {
        std::vector<std::complex<double>> signal(n);for(std::size_t j=0;j<n;++j)signal[j]=cycle[j*d+c];
        transforms.push_back(core::discreteFourierTransform(signal));
    }
    for(std::size_t k=0;k<n;++k) {
        double mass=0;for(std::size_t c=0;c<d;++c) {
            auto z=transforms[c][k]/double(n);coefficients.emplace_back(A{z.real(),z.imag()});mass+=std::norm(z);
        }
        masses.emplace_back(mass);total+=mass;angles.emplace_back(2*pi*double(k)/double(n));
        denominators.emplace_back(double(n/std::gcd(k,n)));
        conductorMass[n/std::gcd(k,n)]+=mass;
    }
    A groups;for(const auto& [denominator,mass]:conductorMass)
        groups.emplace_back(O{{"denominator",double(denominator)},{"mass",mass}});
    return O{{"period",double(n)},{"dimension",double(d)},{"angles",angles},{"masses",masses},
        {"coefficients",coefficients},{"reduced_denominators",denominators},{"conductor_masses",groups},{"total_mass",total},
        {"scope","Uniform invariant measure on the supplied finite cycle; not proof that observed dynamics has this period. Angles are radians per return-map step; coefficients are harmonic-major."}};
}
J projectDiscreteCycle(const std::vector<double>& x,std::size_t d,std::size_t m) {
    validate(x,d);const auto n=x.size()/d;
    if(!m||m>n||n%m)throw std::invalid_argument("Quotient period must divide the supplied cycle length");
    std::vector<double> mean(m*d,0),projected(x.size());
    for(std::size_t j=0;j<n;++j)for(std::size_t c=0;c<d;++c)
        mean[(j%m)*d+c]+=x[j*d+c]/double(n/m);
    double squared=0,maximum=0;
    for(std::size_t j=0;j<n;++j){double norm=0;for(std::size_t c=0;c<d;++c){
        const auto k=j*d+c;projected[k]=mean[(j%m)*d+c];
        const double e=x[k]-projected[k];norm+=e*e;
    }squared+=norm;maximum=std::max(maximum,std::sqrt(norm));}
    A harmonics;for(std::size_t k=0;k<m;++k)harmonics.emplace_back(double(k*(n/m)));
    return O{{"period",double(n)},{"quotient_period",double(m)},{"dimension",double(d)},
        {"quotient_values",encode(mean)},{"projected_values",encode(projected)},
        {"retained_dft_indices",harmonics},{"rms_error",std::sqrt(squared/double(n))},
        {"maximum_error",maximum},{"scope","Uniform finite-group conditional expectation in arithmetic coordinates; retains characters factoring through the quotient. Does not imply projected points form a physical orbit."}};
}
J alignDiscreteCycles(const std::vector<double>& x,const std::vector<double>& y,std::size_t d){
    validate(x,d);validate(y,d);const auto n=x.size()/d;
    if(x.size()!=y.size()||n>20000000/x.size())throw std::invalid_argument("Equal bounded cycles required for alignment");
    double best=0,maximum=0;std::size_t shift=0;
    for(std::size_t s=0;s<n;++s){double squared=0,mx=0;
        for(std::size_t j=0;j<n;++j){double norm=0;for(std::size_t c=0;c<d;++c){
            const double e=x[j*d+c]-y[((j+s)%n)*d+c];norm+=e*e;
        }squared+=norm;mx=std::max(mx,std::sqrt(norm));}
        if(!s||squared<best){best=squared;maximum=mx;shift=s;}
    }
    std::vector<double> aligned(x.size());for(std::size_t j=0;j<n;++j)for(std::size_t c=0;c<d;++c)
        aligned[j*d+c]=y[((j+shift)%n)*d+c];
    return O{{"shift",double(shift)},{"rms_error",std::sqrt(best/double(n))},
        {"maximum_error",maximum},{"aligned_values",encode(aligned)},
        {"scope","Minimum arithmetic state-norm error over forward cyclic shifts only; no time reversal, geometric deformation, or orbit-equivalence certificate"}};
}
J evaluateDiscreteCycleQuery(const J& q) {
    const auto x=reals(q.at("values"));const auto d=integer(q.at("dimension"));validate(x,d);
    const auto op=q.at("operation").string();
    if(op=="cycle_spectrum")return discreteCycleSpectrum(x,d);
    if(op=="cycle_projection")return projectDiscreteCycle(x,d,integer(q.at("quotient_period")));
    if(op=="cycle_alignment")return alignDiscreteCycles(x,reals(q.at("candidate")),d);
    if(op=="sampled_hankel_pairs") {
        core::HankelWindow window;window.delays=integer(q.at("delays"));
        window.sampleStride=q.at("sample_stride").isNull()?1:integer(q.at("sample_stride"));
        window.delayStride=q.at("delay_stride").isNull()?1:integer(q.at("delay_stride"));
        window.snapshotStride=q.at("snapshot_stride").isNull()?1:integer(q.at("snapshot_stride"));
        if(!q.at("columns").isNull())window.columns=integer(q.at("columns"));
        if(!window.delays||window.delays>512||d*window.delays>2048)throw std::invalid_argument("Bounded positive delay count required");
        if((x.size()/d)*window.delays*d>4000000)throw std::invalid_argument("Requested Hankel allocation exceeds limit");
        std::vector<core::HankelChannel> channels;for(std::size_t c=0;c<d;++c)channels.push_back({c,0,1.});
        const auto h=core::buildHankelMatrix(x,x.size()/d,d,channels,window,true);
        std::vector<double> a(h.x.size()),b(h.y.size());
        for(std::size_t j=0;j<h.columns;++j)for(std::size_t i=0;i<h.rows;++i) {
            a[j*h.rows+i]=h.x[i*h.columns+j];b[j*h.rows+i]=h.y[i*h.columns+j];
        }
        return O{{"x",encode(a)},{"y",encode(b)},{"dimension",double(h.rows)},{"samples",double(h.columns)},
            {"scope","Native buildHankelMatrix; sample-major returned pairs, oldest block first; successor is one sample-stride tick, independent of snapshot stride."}};
    }
    std::vector<std::size_t> candidates;for(const auto& p:q.at("candidate_periods").array())candidates.push_back(integer(p));
    const auto scales=q.at("scales").isNull()?std::vector<double>(d,1.):reals(q.at("scales"));
    const auto periods=q.at("coordinate_periods").isNull()?std::vector<double>(d,0.):reals(q.at("coordinate_periods"));
    const auto out=diagnoseSampledCycles(x,d,candidates,scales,periods,num(q.at("tolerance")),q.at("minimum_repeats").isNull()?4:integer(q.at("minimum_repeats")));
    A defects;for(const auto& r:out.defects)defects.emplace_back(O{{"period",double(r.period)},{"maximum",r.maximum},{"rms",r.rms},{"comparisons",double(r.comparisons)}});
    return O{{"detected_period",double(out.detectedPeriod)},{"defects",defects},{"folded_values",encode(out.foldedValues)},
        {"folding_rms",out.foldingRms},{"certified",false},
        {"scope","Smallest accepted supplied candidate period in a finite record and explicit scaled norm. Does not prove an infinite cascade, odometer conjugacy, stationarity, or absence of other periods."}};
}
}
