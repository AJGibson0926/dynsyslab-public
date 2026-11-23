#include <dynsyslab/core/BasisFunctions.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::basis {
namespace {
void argument(std::size_t n,double x) {
    if(n>512||!std::isfinite(x))throw std::invalid_argument("Basis degree must be <=512 and argument finite");
}
}
// Recurrences avoid the factorial growth of unnormalized Hermite polynomials.
double chebyshev(std::size_t n,double x) {
    argument(n,x);double a=1,b=x;if(!n)return a;
    for(std::size_t j=2;j<=n;++j){double c=2*x*b-a;a=b;b=c;}return b;
}
double laguerre(std::size_t n,double x) {
    argument(n,x);double a=1,b=1-x;if(!n)return a;
    for(std::size_t j=2;j<=n;++j){double c=((2*j-1-x)*b-(j-1)*a)/j;a=b;b=c;}return b;
}
double hermiteFunction(std::size_t n,double x) {
    argument(n,x);
    if(std::abs(x)>1e6)return 0; // Gaussian dominates every supported degree.
    double a=std::pow(3.14159265358979323846,-.25),b=std::sqrt(2.)*x*a,logScale=0;
    for(std::size_t j=2;j<=n;++j) {
        if(std::max(std::abs(a),std::abs(b))>1e100){a*=1e-100;b*=1e-100;logScale+=100*std::log(10.);}
        const double c=std::sqrt(2./j)*x*b-std::sqrt(double(j-1)/j)*a;a=b;b=c;
    }
    const double p=n?b:a;
    return p==0?0:std::copysign(std::exp(std::log(std::abs(p))+logScale-.5*x*x),p);
}
bool isRadial(const std::string& k) {return k=="gaussian"||k=="laplacian"||k=="lorentzian"||k=="matern";}
bool isProduct(const std::string& k) {return k=="polynomial"||k=="chebyshev"||k=="laguerre"||k=="hermite"||k=="fourier_hermite";}
double radialKernel(const std::string& kind,double distance,double width,double nu,double amplitude) {
    if(!isRadial(kind)||!std::isfinite(distance)||distance<0||!std::isfinite(width)||width<=0
       ||!std::isfinite(amplitude)||amplitude<=0||!std::isfinite(nu)||nu<.5||nu>8)
        throw std::invalid_argument("Invalid radial kernel parameters (Matern nu must be in [0.5,8])");
    const double z=distance/width;
    if(kind=="gaussian")return amplitude*std::exp(-.5*z*z);
    if(kind=="laplacian")return amplitude*std::exp(-z);
    if(kind=="lorentzian")return amplitude/(1+z*z);
    if(z==0)return amplitude;
    if(nu==.5)return amplitude*std::exp(-z);
    if(z>750)return 0;
    // Analytic small-z limit avoids intermediate K_nu overflow.
    if(nu>=2 && z<1e-5)return amplitude*(1-z*z/(4*(nu-1)));
    if(z<1e-150)return amplitude;
    return amplitude*std::exp((1-nu)*std::log(2.)-std::lgamma(nu)+nu*std::log(z)
                            +std::log(boost::math::cyl_bessel_k(nu,z)));
}
}
