#include <dynsyslab/math/InverseScattering.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Lax spectral problem and transfer matrices; Drazin & Johnson, Solitons:
// An Introduction, ch.3, https://doi.org/10.1017/CBO9781139172059.
// Finite-box Sturm bisection: Golub & Van Loan, Matrix Computations, sec.8.5.
namespace dynsyslab::math {
namespace {
void mesh(const std::vector<double>&v,double h){if(v.empty()||v.size()>1000000||!std::isfinite(h)||h<=0)throw std::invalid_argument("Finite nonempty potential and positive cell width required");
    for(double x:v)if(!std::isfinite(x))throw std::invalid_argument("Finite potential required");}
std::array<double,4> cell(double v,double E,double h){
    const double s=E-v;double co,si;
    if(s==0){co=1;si=h;}else if(s>0){const double k=std::sqrt(s);co=std::cos(k*h);si=std::sin(k*h)/k;}
    else{const double k=std::sqrt(-s);co=std::cosh(k*h);si=std::sinh(k*h)/k;}
    if(!std::isfinite(co)||!std::isfinite(si)||!std::isfinite(s*si))throw std::overflow_error("Transfer matrix exceeded finite range; subdivide/restrict the problem");
    return {co,si,-s*si,co};
}
}
HillMonodromy hillMonodromy(const std::vector<double>&v,double h,double E){
    mesh(v,h);if(!std::isfinite(E))throw std::invalid_argument("Finite spectral energy required");
    std::array<double,4>M{1,0,0,1};for(double value:v){const auto C=cell(value,E,h);
        M={C[0]*M[0]+C[1]*M[2],C[0]*M[1]+C[1]*M[3],C[2]*M[0]+C[3]*M[2],C[2]*M[1]+C[3]*M[3]};
        for(auto x:M)if(!std::isfinite(x))throw std::overflow_error("Monodromy overflow; no scaled transfer backend selected");}
    return {M,M[0]+M[3],std::abs(M[0]*M[3]-M[1]*M[2]-1)};
}
ScatteringCoefficient rightScattering(const std::vector<double>&v,double left,double h,double k){
    mesh(v,h);if(!std::isfinite(left)||!std::isfinite(k)||k<=0)throw std::invalid_argument("Finite left endpoint and positive wavenumber required");
    const double right=left+h*v.size();if(!std::isfinite(right)||!std::isfinite(k*right)||!std::isfinite(k*left))throw std::overflow_error("Scattering phase overflow");
    const auto M=hillMonodromy(v,h,k*k).matrix;const std::complex<double> i(0,1),phase=std::exp(-i*k*left);
    const auto psi=(M[0]-i*k*M[1])*phase,derivative=(M[2]-i*k*M[3])*phase;
    const auto a=(psi-derivative/(i*k))*.5*std::exp(i*k*right),b=(psi+derivative/(i*k))*.5*std::exp(-i*k*right);
    if(std::abs(a)==0)throw std::runtime_error("Singular incoming coefficient");const auto R=b/a,T=1./a;
    const double defect=std::abs(std::norm(R)+std::norm(T)-1);
    if(!std::isfinite(defect))throw std::overflow_error("Nonfinite scattering result");return {R,T,defect};
}
std::vector<double> dirichletNegativeSpectrum(const std::vector<double>&v,double length){
    if(v.size()>20000)throw std::invalid_argument("Dirichlet grid limited to 20000 internal nodes");
    const double h=length/(v.size()+1);mesh(v,h);const double b=1/(h*h),square=b*b;
    if(!std::isfinite(square))throw std::invalid_argument("Dirichlet grid spacing outside representable range");
    auto count=[&](double E){std::size_t n=0;double pivot=0;
        for(std::size_t j=0;j<v.size();++j){pivot=v[j]+2*b-E-(j?square/pivot:0);
            if(std::abs(pivot)<std::numeric_limits<double>::min()*1e100)pivot=-std::numeric_limits<double>::min()*1e100;
            if(pivot<0)++n;}return n;};
    const auto total=count(0);if(total>1000)throw std::invalid_argument("More than 1000 negative finite-box eigenvalues; narrow the requested problem");
    const double bound=*std::min_element(v.begin(),v.end())-std::max(1.,b)*1e-12;std::vector<double>out;
    for(std::size_t j=0;j<total;++j){double lo=bound,hi=0;
        for(int iter=0;iter<100;++iter){const double mid=(lo+hi)/2;if(mid==lo||mid==hi)break;if(count(mid)<=j)lo=mid;else hi=mid;}
        out.push_back((lo+hi)/2);}return out;
}
}
