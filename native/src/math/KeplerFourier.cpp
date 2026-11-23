#include <dynsyslab/math/AnalyticFourier.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <stdexcept>
// Fourier-Bessel inversion of Kepler's equation:
// Fitzpatrick, Celestial Mechanics, Appendix A, Elliptic expansions:
// https://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node126.html
// c_n(cos E)=[J_{n-1}(ne)-J_{n+1}(ne)]/(2n), c_0=-e/2.
// c_n(sin E)=-i[J_{n-1}(ne)+J_{n+1}(ne)]/(2n), n>0.
namespace dynsyslab::math {
FourierSeries keplerAnomalyFourier(bool sine,double e,int N,double omega) {
    if(!std::isfinite(e)||e<0||e>=1||N<0||N>4096||!std::isfinite(omega)||omega<=0)
        throw std::invalid_argument("Kepler Fourier requires 0<=e<1, 0<=N<=4096 and positive frequency");
    FourierSeries out{omega,{}};
    if(!sine)out.coefficients[0]=-e/2;
    for(int n=1;n<=N;++n) {
        const double a=boost::math::cyl_bessel_j(n-1,n*e);
        const double b=boost::math::cyl_bessel_j(n+1,n*e);
        const auto c=sine?std::complex<double>(0,-(a+b)/(2*n)):std::complex<double>((a-b)/(2*n),0);
        if(!std::isfinite(c.real())||!std::isfinite(c.imag()))throw std::overflow_error("Kepler Fourier coefficient overflow");
        out.coefficients[n]=c;out.coefficients[-n]=std::conj(c);
    }return out;
}
}
