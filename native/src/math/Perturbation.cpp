#include <dynsyslab/math/Perturbation.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

// Poincare-Dulac homological operator: Arnold, Geometrical Methods in ODE,
// https://doi.org/10.1007/978-1-4612-1037-5.
// Melnikov convention: Teschl, ODE and Dynamical Systems, section 13.3,
// https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf.
// Finite small-divisor diagnostics are prerequisites, not a KAM proof;
// Figueras/Haro/Luque https://arxiv.org/abs/1601.00084.
namespace dynsyslab::math {
std::vector<SmallDivisorRecord> smallDivisors(const std::vector<double>&w,unsigned order,double tau,std::size_t retain){
    if(w.empty()||w.size()>8||!order||order>1000||retain>100000||!std::isfinite(tau)||tau<0)throw std::invalid_argument("Invalid small-divisor search");
    for(auto f:w)if(!std::isfinite(f))throw std::invalid_argument("Finite frequencies required");
    // Bound the enclosing cube before enumeration to prevent explosive searches.
    std::size_t cap=1;for(std::size_t j=0;j<w.size();++j){if(cap>2000000/(2*order+1))throw std::invalid_argument("Small-divisor search exceeds two million lattice candidates");cap*=2*order+1;}
    std::vector<SmallDivisorRecord> result;FourierIndex k(w.size());
    std::function<void(std::size_t,unsigned)> walk=[&](std::size_t axis,unsigned used){
        if(axis==w.size()){
            if(!used)return;for(int n:k){if(n<0)return;if(n>0)break;} // choose one of +/-k
            const double divisor=harmonicFrequency(k,w),scaled=std::abs(divisor)*std::pow(double(used),tau);
            if(!std::isfinite(scaled))throw std::overflow_error("Scaled small divisor overflow");result.push_back({k,divisor,scaled});return;}
        const int remaining=int(order-used);for(int n=-remaining;n<=remaining;++n){k[axis]=n;walk(axis+1,used+unsigned(std::abs(n)));}};
    walk(0,0);std::sort(result.begin(),result.end(),[](const auto&a,const auto&b){return a.scaledDivisor<b.scaledDivisor;});
    if(result.size()>retain)result.resize(retain);return result;
}
NormalFormStep homogeneousNormalForm(const std::vector<std::complex<double>>&lambda,const std::vector<NormalFormMonomial>&terms,double tol){
    if(lambda.empty()||lambda.size()>128||terms.size()>100000||!std::isfinite(tol)||tol<0)throw std::invalid_argument("Invalid normal-form input");
    for(auto z:lambda)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Finite eigenvalues required");
    NormalFormStep result;result.minimumRemovedDivisor=std::numeric_limits<double>::infinity();
    for(auto term:terms){if(term.component>=lambda.size()||term.powers.size()!=lambda.size()||!std::isfinite(term.coefficient.real())||!std::isfinite(term.coefficient.imag()))throw std::invalid_argument("Invalid normal-form monomial");
        unsigned degree=0;std::complex<double> divisor=-lambda[term.component];for(std::size_t j=0;j<lambda.size();++j){
            if(term.powers[j]>64)throw std::invalid_argument("Normal-form degree limited to 64");degree+=term.powers[j];divisor+=double(term.powers[j])*lambda[j];}
        if(degree<2||degree>64||(result.degree&&result.degree!=degree))throw std::invalid_argument("One homogeneous degree 2..64 required");result.degree=degree;
        if(std::abs(divisor)<=tol)result.retained.push_back(term);
        else{result.minimumRemovedDivisor=std::min(result.minimumRemovedDivisor,std::abs(divisor));term.coefficient/=divisor;
            if(!std::isfinite(term.coefficient.real())||!std::isfinite(term.coefficient.imag()))throw std::overflow_error("Normal-form coefficient overflow");result.change.push_back(term);}}
    if(result.change.empty())result.minimumRemovedDivisor=0;return result;
}
double sampledMelnikov(const std::vector<double>&t,const std::vector<double>&f,const std::vector<double>&g){
    if(t.size()<2||t.size()>1000000||f.size()!=2*t.size()||g.size()!=f.size())throw std::invalid_argument("Melnikov sample shape mismatch");
    double sum=0,previous=0;for(std::size_t i=0;i<t.size();++i){const double value=f[2*i]*g[2*i+1]-f[2*i+1]*g[2*i];
        if(!std::isfinite(t[i])||!std::isfinite(value)||(i&&t[i]<=t[i-1]))throw std::invalid_argument("Melnikov requires finite samples and increasing time");
        if(i)sum+=(t[i]-t[i-1])*(value+previous)/2;previous=value;}
    if(!std::isfinite(sum))throw std::overflow_error("Melnikov integral overflow");return sum;
}
double conformalSymplecticDefect(const std::vector<double>&A,const std::vector<double>&J,std::size_t d,double gamma){
    if(!d||d>128||A.size()!=d*d||J.size()!=d*d||!std::isfinite(gamma))throw std::invalid_argument("Symplectic defect shape mismatch");
    for(auto v:A)if(!std::isfinite(v))throw std::invalid_argument("Finite Jacobian required");
    for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j)if(!std::isfinite(J[i*d+j])||std::abs(J[i*d+j]+J[j*d+i])>1e-12)throw std::invalid_argument("Skew-symmetric form required");
    double total=0;for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j){double z=gamma*J[i*d+j];
        for(std::size_t k=0;k<d;++k)z+=A[k*d+i]*J[k*d+j]+J[i*d+k]*A[k*d+j];total=std::hypot(total,z);}
    if(!std::isfinite(total))throw std::overflow_error("Symplectic defect overflow");return total;
}
}
