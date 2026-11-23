#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double pi=std::acos(-1.);
constexpr std::size_t maximumLatticePoints=200000;
bool finite(C z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
void requireFinite(const ComplexVector& v,const char* message){for(const auto& z:v)if(!finite(z))throw std::invalid_argument(message);}
// All k in N0^d with |k|_1 <= order.
void compositions(std::size_t d,int order,std::vector<std::vector<int>>& out){
    std::vector<int> k(d,0);
    const auto recurse=[&](const auto& self,std::size_t i,int left)->void{
        if(i==d){out.push_back(k);return;}
        for(int v=0;v<=left;++v){k[i]=v;self(self,i+1,left-v);}
        k[i]=0;
    };
    recurse(recurse,0,order);
}
void harmonics(std::size_t d,int bound,std::vector<std::vector<int>>& out){
    std::vector<int> m(d,-bound);
    if(d==0){out.emplace_back();return;}
    while(true){
        out.push_back(m);
        std::size_t i=0;
        for(;i<d;++i){if(m[i]<bound){++m[i];break;}m[i]=-bound;}
        if(i==d)break;
    }
}
PoleLattice finishLattice(std::vector<LatticePole> poles,double tolerance,std::string scope){
    if(!(tolerance>=0)||!std::isfinite(tolerance))throw std::invalid_argument("Coincidence tolerance must be nonnegative");
    std::sort(poles.begin(),poles.end(),[](const LatticePole& a,const LatticePole& b){
        if(a.value.real()!=b.value.real())return a.value.real()>b.value.real();
        if(std::abs(a.value.imag())!=std::abs(b.value.imag()))return std::abs(a.value.imag())<std::abs(b.value.imag());
        if(a.value.imag()!=b.value.imag())return a.value.imag()<b.value.imag();
        return a.k<b.k||(a.k==b.k&&a.m<b.m);});
    PoleLattice out;out.scope=std::move(scope);
    // Coincidences by a sweep in Re; candidates are within tolerance in the sort key.
    for(std::size_t a=0;a<poles.size();++a)for(std::size_t b=a+1;b<poles.size();++b){
        if(poles[a].value.real()-poles[b].value.real()>tolerance)break;
        const double d=std::abs(poles[a].value-poles[b].value);
        if(d<=tolerance)out.coincidences.push_back({a,b,d});
        if(out.coincidences.size()>=1000)break;
    }
    out.poles=std::move(poles);
    return out;
}
void requireStable(const ComplexVector& v,const char* name){
    requireFinite(v,"Exponents must be finite");
    for(const auto& z:v)if(!(z.real()<0))throw std::invalid_argument(std::string(name)+" need negative real parts: the expansion theorem does not apply");
}
}
void validatePoleModel(const PoleResidueModel& model){
    if(model.poles.empty())throw std::invalid_argument("Pole model needs at least one pole");
    for(const auto& p:model.poles){
        if(!finite(p.eigenvalue))throw std::invalid_argument("Pole must be finite");
        if(p.laurent.empty()||p.laurent.size()>16)throw std::invalid_argument("Pole Laurent order must be 1..16");
        requireFinite(p.laurent,"Laurent coefficients must be finite");
    }
}
C evaluateResolvent(const PoleResidueModel& model,C s){
    validatePoleModel(model);C sum{};
    for(const auto& p:model.poles){
        const C d=s-p.eigenvalue;if(d==C{})throw std::domain_error("Resolvent evaluated at a pole");
        const C w=1./d;C power=w;
        for(const auto& a:p.laurent){sum+=a*power;power*=w;}
    }
    return sum;
}
C evaluateTimeResponse(const PoleResidueModel& model,double t){
    validatePoleModel(model);if(!(t>=0)||!std::isfinite(t))throw std::invalid_argument("Time must be nonnegative");
    C sum{};
    for(const auto& p:model.poles){
        C inner{};double factor=1;
        for(std::size_t q=0;q<p.laurent.size();++q){inner+=p.laurent[q]*factor;factor*=t/double(q+1);}
        sum+=inner*std::exp(p.eigenvalue*t);
    }
    return sum;
}
double modelAbscissa(const PoleResidueModel& model){
    validatePoleModel(model);double top=-std::numeric_limits<double>::infinity();
    for(const auto& p:model.poles){
        if(std::any_of(p.laurent.begin(),p.laurent.end(),[](C a){return a!=C{};}))top=std::max(top,p.eigenvalue.real());
    }
    return top;
}
LaurentEstimate cauchyLaurent(const LaplaceFunction& F,C center,double radius,const CauchyOptions& o){
    if(!finite(center)||!(radius>0)||!std::isfinite(radius))throw std::invalid_argument("Contour needs a finite centre and positive radius");
    if(o.maximumOrder<1||o.maximumOrder>16)throw std::invalid_argument("Laurent order must be 1..16");
    if(o.initialNodes<8||o.maximumNodes<o.initialNodes||o.maximumNodes>1u<<20)throw std::invalid_argument("Invalid contour node counts");
    const std::size_t Q=o.maximumOrder+1;
    // S_q accumulates F(c+r e^{i th}) (r e^{i th})^q over the nodes; nodes nest under doubling.
    ComplexVector sums(Q+1);
    const auto add=[&](double theta){
        const C w=std::polar(radius,theta);const C value=F(center+w);
        if(!finite(value))throw std::domain_error("Non-finite value on the contour: a pole lies on the circle");
        C power=w;for(std::size_t q=1;q<=Q;++q){sums[q]+=value*power;power*=w;}
    };
    std::size_t nodes=o.initialNodes;
    for(std::size_t k=0;k<nodes;++k)add(2*pi*double(k)/double(nodes));
    const auto coefficients=[&](std::size_t count){ComplexVector a(Q);for(std::size_t q=1;q<=Q;++q)a[q-1]=sums[q]/double(count);return a;};
    ComplexVector previous=coefficients(nodes),current=previous;
    LaurentEstimate out;out.center=center;out.radius=radius;double change=std::numeric_limits<double>::infinity();
    while(nodes*2<=o.maximumNodes){
        for(std::size_t k=0;k<nodes;++k)add(2*pi*(double(k)+.5)/double(nodes));
        nodes*=2;current=coefficients(nodes);
        double scale=0,difference=0;
        for(std::size_t q=0;q<Q;++q){scale=std::max(scale,std::abs(current[q]));difference=std::max(difference,std::abs(current[q]-previous[q]));}
        change=scale>0?difference/scale:difference;previous=current;
        if(change<=o.tolerance){out.converged=true;break;}
    }
    out.nodes=nodes;out.relativeChange=change;
    double scale=0;for(std::size_t q=0;q<o.maximumOrder;++q)scale=std::max(scale,std::abs(current[q]));
    out.laurent.assign(current.begin(),current.begin()+std::ptrdiff_t(o.maximumOrder));
    out.higherOrderRatio=scale>0?std::abs(current[o.maximumOrder])/scale:0;
    for(std::size_t q=o.maximumOrder;q-->0;)if(std::abs(current[q])>o.orderTolerance*scale&&scale>0){out.order=q+1;break;}
    return out;
}
PoleLattice equilibriumPoleLattice(const ComplexVector& lambda,int order,double tolerance){
    if(lambda.empty()||order<0||order>64)throw std::invalid_argument("Equilibrium lattice needs exponents and order 0..64");
    requireStable(lambda,"Equilibrium eigenvalues");
    std::vector<std::vector<int>> indices;compositions(lambda.size(),order,indices);
    if(indices.size()>maximumLatticePoints)throw std::invalid_argument("Lattice is too large");
    std::vector<LatticePole> poles;poles.reserve(indices.size());
    for(auto& k:indices){C v{};for(std::size_t j=0;j<k.size();++j)v+=double(k[j])*lambda[j];poles.push_back({v,std::move(k),{}});}
    return finishLattice(std::move(poles),tolerance,
        "Thm 8: needs an analytic field, a stable equilibrium, and distinct non-resonant eigenvalues; coincidences are resonances where the expansion can gain Jordan terms");
}
PoleLattice limitCyclePoleLattice(const ComplexVector& nu,double Omega,int order,int harmonic,double tolerance){
    if(!(Omega>0)||!std::isfinite(Omega))throw std::invalid_argument("Limit-cycle frequency must be positive");
    auto out=torusPoleLattice(nu,{Omega},order,harmonic,tolerance);
    out.smallestHarmonicDivisor=Omega;
    out.scope="Thm 10: stable limit cycle, positive distinct multipliers in the paper; a complex exponent is defined modulo i Omega and the lattice set is unchanged";
    return out;
}
PoleLattice torusPoleLattice(const ComplexVector& nu,const std::vector<double>& Omega,int order,int harmonic,double tolerance){
    if(nu.empty()||Omega.empty()||order<0||order>64||harmonic<0||harmonic>4096)throw std::invalid_argument("Torus lattice needs exponents, frequencies and valid bounds");
    requireStable(nu,"Transverse exponents");
    for(double w:Omega)if(!std::isfinite(w)||w==0)throw std::invalid_argument("Frequencies must be finite and nonzero");
    std::vector<std::vector<int>> ks,ms;compositions(nu.size(),order,ks);
    double harmonicCount=1;for(std::size_t i=0;i<Omega.size();++i)harmonicCount*=2.*harmonic+1;
    if(double(ks.size())*harmonicCount>double(maximumLatticePoints))throw std::invalid_argument("Lattice is too large");
    harmonics(Omega.size(),harmonic,ms);
    std::vector<LatticePole> poles;poles.reserve(ks.size()*ms.size());
    std::optional<double> divisor;
    for(const auto& m:ms){double w=0;bool zero=true;
        for(std::size_t j=0;j<m.size();++j){w+=double(m[j])*Omega[j];zero=zero&&m[j]==0;}
        if(!zero){const double a=std::abs(w);if(!divisor||a<*divisor)divisor=a;}
        for(const auto& k:ks){C v(0,w);for(std::size_t j=0;j<k.size();++j)v+=double(k[j])*nu[j];poles.push_back({v,k,m});}
    }
    auto out=finishLattice(std::move(poles),tolerance,
        "Thm 11: finite truncation of k.nu + i m.Omega. With rationally independent frequencies the imaginary lattice is dense, so residues cannot be isolated by small contours on Re s=0");
    out.smallestHarmonicDivisor=divisor;
    return out;
}
double resolventNorm(const ComplexMatrix& A,C s){
    const auto sv=complexSingularValues(shiftedMatrix(A,s));
    if(sv.empty())throw std::invalid_argument("Empty matrix");
    return sv.back()>0?1/sv.back():std::numeric_limits<double>::infinity();
}
PseudospectrumGrid pseudospectrumGrid(const ComplexMatrix& A,std::vector<double> re,std::vector<double> im){
    if(re.empty()||im.empty()||re.size()*im.size()>1000000)throw std::invalid_argument("Pseudospectrum grid size is out of range");
    if(A.rows==0||A.rows>512)throw std::invalid_argument("Pseudospectrum matrix must have dimension 1..512");
    PseudospectrumGrid out;out.eigenvalues=complexEigenvalues(A);
    out.sigmaMin.reserve(re.size()*im.size());
    for(double y:im)for(double x:re){
        if(!std::isfinite(x)||!std::isfinite(y))throw std::invalid_argument("Grid must be finite");
        out.sigmaMin.push_back(complexSingularValues(shiftedMatrix(A,C(x,y))).back());
    }
    out.realAxis=std::move(re);out.imagAxis=std::move(im);
    return out;
}
}
