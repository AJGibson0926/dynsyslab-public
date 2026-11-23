#include <dynsyslab/math/Arithmetic.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <stdexcept>

// Exact Euclidean/valuation algorithms. Sources: Gouvea, p-adic Numbers
// (Springer), chapters 1-2; Ford, Fractions, Amer. Math. Monthly 45 (1938).
namespace dynsyslab::math {
namespace {
using I=boost::multiprecision::cpp_int;
I absI(I x){return x<0?-x:x;}
I gcd(I a,I b){a=absI(a);b=absI(b);while(b!=0){I r=a%b;a=b;b=r;}return a;}
I integer(const std::string& s){
    if(s.empty()||s.size()>4096)throw std::invalid_argument("Exact integer must have 1..4096 characters");
    std::size_t j=(s[0]=='-'||s[0]=='+')?1:0;
    if(j==s.size())throw std::invalid_argument("Invalid integer");
    I n=0;for(;j<s.size();++j){if(s[j]<'0'||s[j]>'9')throw std::invalid_argument("Invalid integer");n=n*10+(s[j]-'0');}
    return s[0]=='-'?-n:n;
}
struct R {
    I n,d;
    R(I num,I den):n(std::move(num)),d(std::move(den)){
        if(d==0)throw std::invalid_argument("Zero rational denominator");
        if(d<0){n=-n;d=-d;}I g=gcd(n,d);n/=g;d/=g;
    }
    std::string str()const{return n.str()+(d==1?"":"/"+d.str());}
};
R read(const std::string& s){auto p=s.find('/');return p==std::string::npos?R(integer(s),1):R(integer(s.substr(0,p)),integer(s.substr(p+1)));}
void prime(std::uint32_t p){
    if(p<2||p>1000000)throw std::invalid_argument("Prime must be between 2 and 1000000");
    for(std::uint32_t q=2;q<=p/q;++q)if(p%q==0)throw std::invalid_argument("p must be prime");
}
I power(std::uint32_t p,std::size_t n){I result=1,base=p;while(n){if(n&1)result*=base;n>>=1;if(n)base*=base;}return result;}
int strip(I& n,std::uint32_t p){int v=0;while(n!=0&&n%p==0){n/=p;++v;}return v;}
I mod(I a,const I& m){a%=m;if(a<0)a+=m;return a;}
I inverse(I a,const I& m){
    I r=m,s=mod(a,m),x=0,y=1;
    while(s!=0){I q=r/s,t=r-q*s;r=s;s=t;t=x-q*y;x=y;y=t;}
    if(r!=1)throw std::invalid_argument("Modular inverse does not exist");return mod(x,m);
}
}
std::string canonicalRational(const std::string& v){return read(v).str();}
PhaseCoverResonance phaseCoverResonance(std::uint32_t p,const std::string& ratio){
    // Integer-kernel calculation on the character group Z[1/p] x Z.
    // Circle-cover character groups: Burgos/Verjovsky, arXiv:1603.05676, section 2.
    prime(p);auto alpha=read(ratio);
    if(alpha.n<=0)throw std::invalid_argument("Positive exact base/drive ratio required");
    I unit=alpha.n;const auto depth=std::size_t(strip(unit,p));const I divisor=power(p,depth);
    const I drive=-unit;
    return {alpha.str(),R(alpha.d,divisor).str(),drive.str(),divisor.str(),depth};
}
std::string rationalDifference(const std::string& a,const std::string& b){auto x=read(a),y=read(b);return R(x.n*y.d-y.n*x.d,x.d*y.d).str();}
std::string rationalRatio(const std::string& a,const std::string& b){auto x=read(a),y=read(b);return R(x.n*y.d,x.d*y.n).str();}
std::optional<int> primeValuation(const std::string& v,std::uint32_t p){
    prime(p);auto x=read(v);if(x.n==0)return {};return strip(x.n,p)-strip(x.d,p);
}
std::string padicDistance(const std::string& a,const std::string& b,std::uint32_t p){
    prime(p);auto x=read(a),y=read(b);R delta(x.n*y.d-y.n*x.d,x.d*y.d);
    if(delta.n==0)return "0";
    const int v=strip(delta.n,p)-strip(delta.d,p);
    return v<0?power(p,std::size_t(-v)).str():R(1,power(p,std::size_t(v))).str();
}
PadicExpansion padicExpansion(const std::string& v,std::uint32_t p,std::size_t precision){
    prime(p);if(precision==0||precision>4096)throw std::invalid_argument("p-adic precision must be 1..4096 digits");
    auto x=read(v);PadicExpansion result{p,{},std::vector<std::uint32_t>(precision)};
    if(x.n==0)return result;
    result.valuation=strip(x.n,p)-strip(x.d,p);
    I modulus=power(p,precision),unit=mod(x.n*inverse(x.d,modulus),modulus);
    for(auto& digit:result.digits){digit=(unit%p).convert_to<std::uint32_t>();unit/=p;}
    return result;
}
bool compatiblePhaseResidues(const std::string& a,const std::string& b,std::uint32_t p,std::size_t depth){
    prime(p);if(depth>4096)throw std::invalid_argument("Phase depth exceeds 4096");
    I lower=integer(a),upper=integer(b),m=power(p,depth);
    if(lower<0||lower>=m||upper<0||upper>=m*p)throw std::invalid_argument("Phase residues outside their declared levels");
    return upper%m==lower;
}
std::vector<PhaseCoverLevel> phaseCoverTower(std::uint32_t p,const std::vector<std::uint32_t>& digits){
    // Circle-cover inverse limits; Garrett, Adele groups, p-adic groups,
    // solenoids (2014), sections 6 and 8. This is a finite arithmetic prefix,
    // not an assertion that a physical system has the associated orbit family.
    prime(p);if(digits.size()>512)throw std::invalid_argument("Phase prefix exceeds 512 digits");
    I residue=0,modulus=1;std::vector<PhaseCoverLevel> out{{"1","0","1","0"}};
    for(auto digit:digits){
        if(digit>=p)throw std::invalid_argument("Phase branch digit must be less than the prime");
        residue+=digit*modulus;modulus*=p;
        out.push_back({modulus.str(),residue.str(),R(1,modulus).str(),R(residue,modulus).str()});
    }return out;
}
PhaseCoverFlow advancePhaseCover(std::uint32_t p,const std::vector<std::uint32_t>& digits,
                                const std::string& phase,const std::string& elapsed){
    const auto initial=phaseCoverTower(p,digits);const auto s=read(phase),t=read(elapsed);
    if(s.n<0||s.n>=s.d)throw std::invalid_argument("Base phase turns must be in [0,1)");
    const R total(s.n*t.d+t.n*s.d,s.d*t.d);
    I ticks=total.n/total.d;if(total.n<0&&total.n%total.d!=0)--ticks;
    const R fractional(total.n-ticks*total.d,total.d);
    const I modulus=integer(initial.back().modulus),lifted=integer(initial.back().residue)+ticks;
    I current=mod(lifted,modulus);const I carry=(lifted-current)/modulus;
    PhaseCoverFlow out{fractional.str(),ticks.str(),carry.str(),{}, {}, {}};
    for(std::size_t j=0;j<digits.size();++j){out.digits.push_back((current%p).convert_to<std::uint32_t>());current/=p;}
    out.levels=phaseCoverTower(p,out.digits);
    for(const auto& level:out.levels)out.phaseTurns.push_back(
        R(fractional.n+integer(level.residue)*fractional.d,fractional.d*integer(level.modulus)).str());
    return out;
}
PhaseCoverCharacter phaseCoverCharacter(std::uint32_t p,const std::vector<std::uint32_t>& digits,
                                        const std::string& phase,const std::string& frequency){
    // Dual characters of circle-cover inverse limits: Clark (1999),
    // arXiv:math/9906150, introduction. Integer character order on a factor
    // circle is exactly the condition for descending to that factor.
    const auto levels=phaseCoverTower(p,digits);const auto s=read(phase),q=read(frequency);
    if(s.n<0||s.n>=s.d)throw std::invalid_argument("Base phase turns must be in [0,1)");
    I denominator=q.d;const auto depth=std::size_t(strip(denominator,p));
    if(denominator!=1)throw std::invalid_argument("Character frequency must belong to Z[1/p]");
    PhaseCoverCharacter out{q.str(),depth,{}};
    if(depth<=digits.size()){
        const R theta((s.n+integer(levels.back().residue)*s.d)*q.n,s.d*q.d);
        out.phaseTurns=R(mod(theta.n,theta.d),theta.d).str();
    }
    return out;
}
ContinuedFraction continuedFraction(const std::string& s){
    auto x=read(s);if(x.n<=0)throw std::invalid_argument("Positive rational required");
    ContinuedFraction out;I p0=0,p1=1,q0=1,q1=0;
    while(x.d!=0){I a=x.n/x.d,r=x.n%x.d,p=a*p1+p0,q=a*q1+q0;
        out.terms.push_back(a.str());out.convergents.push_back(R(p,q).str());
        p0=p1;p1=p;q0=q1;q1=q;x.n=x.d;x.d=r;}
    return out;
}
std::vector<SternBrocotRun> sternBrocotPath(const std::string& s){
    auto r=read(s);if(r.n<=0)throw std::invalid_argument("Positive rational required");
    std::vector<SternBrocotRun> result;
    while(r.n!=r.d){if(r.n>r.d){I k=(r.n-1)/r.d;result.push_back({'R',k.str()});r.n-=k*r.d;}
        else{I k=(r.d-1)/r.n;result.push_back({'L',k.str()});r.d-=k*r.n;}}
    return result;
}
std::string rationalMediant(const std::string& a,const std::string& b){auto x=read(a),y=read(b);return R(x.n+y.n,x.d+y.d).str();}
FordCircle fordCircle(const std::string& s){auto x=read(s);return {x.str(),R(1,2*x.d*x.d).str()};}
bool fordCirclesTangent(const std::string& a,const std::string& b){auto x=read(a),y=read(b);return absI(x.n*y.d-y.n*x.d)==1;}
std::string multiplicativeDistanceProduct(const std::string& a,const std::string& b){auto x=read(a),y=read(b);R r(x.n*y.d,x.d*y.n);if(r.n<=0)throw std::invalid_argument("Positive ratio required");I product=r.n*r.d;return product.str();}
}
