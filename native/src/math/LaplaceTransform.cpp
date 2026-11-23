#include <dynsyslab/math/LaplaceTransform.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
void validate(const SampledSignal& g){
    if(g.times.size()<2||g.times.size()!=g.values.size())throw std::invalid_argument("Signal needs at least two matching samples");
    if(!g.derivatives.empty()&&g.derivatives.size()!=g.values.size())throw std::invalid_argument("Derivative sample count mismatch");
    if(!std::isfinite(g.times[0])||g.times[0]<0)throw std::invalid_argument("One-sided transform needs nonnegative times");
    for(std::size_t k=1;k<g.times.size();++k)if(!std::isfinite(g.times[k])||g.times[k]<=g.times[k-1])throw std::invalid_argument("Times must increase strictly");
    for(const auto& v:g.values)if(!std::isfinite(v.real())||!std::isfinite(v.imag()))throw std::invalid_argument("Signal values must be finite");
}
// M_k(z)=int_0^1 x^k exp(-z x) dx for k=0..3. Taylor series near zero; the forward
// recurrence M_k=(k M_{k-1}-e^{-z})/z is stable for |z|>=1.
std::array<C,4> panelMoments(C z){
    std::array<C,4> m{};
    if(std::abs(z)<1){
        for(int k=0;k<4;++k){C term=1,sum=1./(k+1.);
            for(int j=1;j<32;++j){term*=-z/double(j);sum+=term/double(k+j+1);}
            m[k]=sum;}
        return m;
    }
    const C e=std::exp(-z);m[0]=(1.-e)/z;
    for(int k=1;k<4;++k)m[k]=(double(k)*m[k-1]-e)/z;
    return m;
}
}
ComplexVector hermiteNodeDerivatives(const std::vector<double>& x,const ComplexVector& y){
    const std::size_t n=x.size();if(n<2||y.size()!=n)throw std::invalid_argument("Derivative needs at least two samples");
    ComplexVector d(n);const std::size_t width=std::min<std::size_t>(5,n);
    for(std::size_t i=0;i<n;++i){
        std::size_t lo=i>width/2?i-width/2:0;if(lo+width>n)lo=n-width;
        C sum{};double diagonal=0;
        for(std::size_t j=lo;j<lo+width;++j){if(j==i)continue;
            double w=1;
            for(std::size_t k=lo;k<lo+width;++k){if(k==j)continue;w/=x[j]-x[k];if(k!=i)w*=x[i]-x[k];}
            sum+=w*y[j];diagonal-=w;}
        d[i]=sum+diagonal*y[i];
    }
    return d;
}
C laplaceWindowIntegral(const SampledSignal& g,C s,LaplaceInterpolant mode){
    validate(g);if(!std::isfinite(s.real())||!std::isfinite(s.imag()))throw std::invalid_argument("Finite s required");
    ComplexVector d;
    if(mode==LaplaceInterpolant::CubicHermite)d=g.derivatives.empty()?hermiteNodeDerivatives(g.times,g.values):g.derivatives;
    C sum{};
    for(std::size_t k=0;k+1<g.times.size();++k){
        const double a=g.times[k],h=g.times[k+1]-a;
        const auto m=panelMoments(s*h);const C y0=g.values[k],delta=g.values[k+1]-y0;C panel;
        if(mode==LaplaceInterpolant::Linear)panel=y0*m[0]+delta*m[1];
        else{const C d0=h*d[k],d1=h*d[k+1];
            panel=y0*m[0]+d0*m[1]+(3.*delta-2.*d0-d1)*m[2]+(d0+d1-2.*delta)*m[3];}
        sum+=h*std::exp(-s*a)*panel;
    }
    return sum;
}
GrowthBound fitGrowthBound(const SampledSignal& g,std::size_t blocks){
    validate(g);if(blocks<2)throw std::invalid_argument("Growth fit needs at least two blocks");
    const std::size_t n=g.times.size();blocks=std::min(blocks,n);
    std::vector<double> tx,ly;
    for(std::size_t b=0;b<blocks;++b){
        const std::size_t lo=b*n/blocks,hi=(b+1)*n/blocks;double best=0,at=g.times[lo];
        for(std::size_t k=lo;k<hi;++k){const double a=std::abs(g.values[k]);if(a>best){best=a;at=g.times[k];}}
        if(best>0){tx.push_back(at);ly.push_back(std::log(best));}
    }
    GrowthBound out;
    if(tx.empty())return out;
    if(tx.size()>=2){
        double mx=0,my=0;for(std::size_t k=0;k<tx.size();++k){mx+=tx[k];my+=ly[k];}
        mx/=double(tx.size());my/=double(tx.size());double sxx=0,sxy=0;
        for(std::size_t k=0;k<tx.size();++k){sxx+=(tx[k]-mx)*(tx[k]-mx);sxy+=(tx[k]-mx)*(ly[k]-my);}
        out.rate=sxx>0?sxy/sxx:0;
    }
    for(std::size_t k=0;k<n;++k)out.constant=std::max(out.constant,std::abs(g.values[k])*std::exp(-out.rate*g.times[k]));
    return out;
}
LaplaceEstimate laplaceTransform(const SampledSignal& g,C s,const LaplaceTransformOptions& o){
    LaplaceEstimate out;out.s=s;out.value=laplaceWindowIntegral(g,s,o.interpolant);
    if(o.growth){
        if(!std::isfinite(o.growth->rate)||!std::isfinite(o.growth->constant)||o.growth->constant<0)throw std::invalid_argument("Invalid declared growth bound");
        out.growth=*o.growth;out.growth.declared=true;
    }else out.growth=fitGrowthBound(g);
    out.insideRegion=s.real()>out.growth.rate;
    if(out.insideRegion){const double gap=s.real()-out.growth.rate;
        out.tailBound=out.growth.constant*std::exp(-gap*g.times.back())/gap;}
    return out;
}
C zTransform(const ComplexVector& values,C z){
    if(values.empty())throw std::invalid_argument("z-transform needs samples");
    if(z==C{})throw std::invalid_argument("z-transform needs z != 0");
    const C w=1./z;C sum{};
    for(std::size_t k=values.size();k-->0;)sum=sum*w+values[k];
    return sum;
}
C zFromLaplace(C s,double h){
    if(!(h>0)||!std::isfinite(h))throw std::invalid_argument("Positive sample interval required");
    return std::exp(s*h);
}
C laplaceFromZ(C z,double h,int branch){
    if(!(h>0)||!std::isfinite(h))throw std::invalid_argument("Positive sample interval required");
    if(z==C{})throw std::invalid_argument("log z needs z != 0");
    return C(std::log(std::abs(z)),std::arg(z)+2*pi*branch)/h;
}
namespace {
// z(theta)=(N/t)(a theta cot(b theta)+c+i d theta): Trefethen, Weideman & Schmelzer (2006).
constexpr double talbotA=.5017,talbotB=.6407,talbotC=-.6122,talbotD=.2645;
C talbotPoint(double theta,double scale){
    return scale*C(talbotA*theta/std::tan(talbotB*theta)+talbotC,talbotD*theta);
}
C talbotSlope(double theta,double scale){
    const double sn=std::sin(talbotB*theta);
    return scale*C(talbotA/std::tan(talbotB*theta)-talbotA*talbotB*theta/(sn*sn),talbotD);
}
bool insideTalbot(C p,double scale){
    const double theta=p.imag()/(scale*talbotD);
    if(std::abs(theta)>=pi)return false;
    const double edge=std::abs(theta)<1e-9?scale*(talbotA/talbotB+talbotC):talbotPoint(theta,scale).real();
    return p.real()<edge;
}
}
InverseLaplaceResult talbotInverse(const LaplaceFunction& F,double t,const InverseLaplaceOptions& o){
    if(!(t>0)||!std::isfinite(t))throw std::invalid_argument("Inverse Laplace needs t>0");
    if(o.nodes<4||o.nodes>400||o.nodes%2)throw std::invalid_argument("Talbot needs an even node count in 4..400");
    if(!std::isfinite(o.shift))throw std::invalid_argument("Finite contour shift required");
    const double N=double(o.nodes),scale=N/t;C sum{};
    for(std::size_t k=0;k<o.nodes;++k){
        const double theta=-pi+(2.*double(k)+1.)*pi/N; // Even N keeps theta away from zero.
        const C z=talbotPoint(theta,scale),s=z+o.shift;
        sum+=std::exp(s*t)*F(s)*talbotSlope(theta,scale);
    }
    InverseLaplaceResult out;out.value=sum/C(0,N);out.method="talbot";
    out.scopeVerified=!o.singularities.empty();
    for(const auto& p:o.singularities)if(!insideTalbot(p-o.shift,scale))out.scopeVerified=false;
    out.scope=o.singularities.empty()?"No singularities declared: contour enclosure is unverified":
        out.scopeVerified?"Every declared singularity is inside the Talbot contour":
        "A declared singularity is outside the contour: this value is not valid; use bromwich or more nodes";
    return out;
}
InverseLaplaceResult bromwichInverse(const LaplaceFunction& F,double t,const InverseLaplaceOptions& o){
    if(!(t>0)||!std::isfinite(t))throw std::invalid_argument("Inverse Laplace needs t>0");
    if(!(o.dampingParameter>0)||!std::isfinite(o.dampingParameter)||!std::isfinite(o.shift))throw std::invalid_argument("Invalid Bromwich parameters");
    if(o.eulerTerms<1||o.eulerTerms>40)throw std::invalid_argument("Euler order must be 1..40");
    double top=0;for(const auto& p:o.singularities)top=std::max(top,std::abs(p.imag()));
    const std::size_t base=o.seriesTerms?o.seriesTerms:std::size_t(std::ceil(top*t/pi))+40;
    if(base>2000000)throw std::invalid_argument("Bromwich series is too long");
    const double c=o.dampingParameter/(2*t);
    const auto term=[&](std::size_t k){const double w=double(k)*pi/t;
        return (k%2?-1.:1.)*(F(C(c+o.shift,w))+F(C(c+o.shift,-w)));};
    C partial=F(C(c+o.shift,0));
    for(std::size_t k=1;k<=base;++k)partial+=term(k);
    C averaged{};double weight=std::pow(2.,-double(o.eulerTerms));
    for(std::size_t k=0;k<=o.eulerTerms;++k){
        averaged+=weight*partial;
        partial+=term(base+k+1);
        weight*=double(o.eulerTerms-k)/double(k+1);
    }
    InverseLaplaceResult out;out.method="bromwich";
    out.value=std::exp((c+o.shift)*t)*averaged/(2*t);
    out.scopeVerified=!o.singularities.empty();
    for(const auto& p:o.singularities)if(p.real()>=c+o.shift)out.scopeVerified=false;
    out.scope=o.singularities.empty()?"No singularities declared: abscissa is unverified":
        out.scopeVerified?"Every declared singularity is left of the Bromwich line":
        "A declared singularity is right of the Bromwich line: this value is not valid; raise shift";
    return out;
}
}
