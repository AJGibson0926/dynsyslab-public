#include <dynsyslab/math/HamiltonianPerturbation.hpp>
#ifdef DYNSYSLAB_MATH_NUMERICS
#include <dynsyslab/numerics/Svd.hpp>
#endif
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

// Canonical Lie transforms and Kolmogorov/isoenergetic nondegeneracy:
// Arnold, Mathematical Methods of Classical Mechanics, ch.10/appendices,
// https://doi.org/10.1007/978-1-4757-2063-1.
// This finite first-order calculation does not supply the analytic or
// Diophantine bounds of an a-posteriori theorem (e.g. arXiv:1601.00084).
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
void checkHessian(const std::vector<double>& w,const std::vector<double>& h) {
    const auto d=w.size();
    if(!d||d>8||h.size()!=d*d)
        throw std::invalid_argument("Hamiltonian frequency/Hessian shape mismatch");
    for(double x:w)if(!std::isfinite(x))throw std::invalid_argument("Finite frequencies required");
    for(double x:h)if(!std::isfinite(x))throw std::invalid_argument("Finite Hessian required");
    for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j)
        if(std::abs(h[i*d+j]-h[j*d+i])>1e-12*std::max({1.,std::abs(h[i*d+j]),std::abs(h[j*d+i])}))
            throw std::invalid_argument("Symmetric Hamiltonian Hessian required");
}
C coefficient(const TorusFourierSeries& s,const FourierIndex& k) {
    const auto it=s.coefficients.find(k);
    return it==s.coefficients.end()?C{}:it->second;
}
void checkReal(const TorusFourierSeries& s) {
    meanSquare(s); // common finite/index/shape validation
    for(const auto& [k,c]:s.coefficients) {
        auto opposite=k;for(auto& n:opposite)n=-n;
        const auto other=coefficient(s,opposite);
        if(std::abs(c-std::conj(other))>1e-12*std::max({1.,std::abs(c),std::abs(other)}))
            throw std::invalid_argument("Real Hamiltonian requires conjugate Fourier coefficients");
    }
}
}
HamiltonianFirstOrder firstOrderHamiltonian(const TorusFourierSeries& h,
    const std::vector<TorusFourierSeries>& gradient,const std::vector<double>& D,double tol) {
    const auto& w=h.angularFrequencies;const auto d=w.size();
    checkHessian(w,D);checkReal(h);
    if(gradient.size()!=d||!std::isfinite(tol)||tol<0)
        throw std::invalid_argument("Invalid Hamiltonian gradient or divisor tolerance");
    std::set<FourierIndex> indices;
    for(const auto& [k,c]:h.coefficients)indices.insert(k);
    for(const auto& g:gradient) {
        checkReal(g);
        if(g.angularFrequencies!=w)throw std::invalid_argument("Gradient frequency mismatch");
        for(const auto& [k,c]:g.coefficients)indices.insert(k);
    }
    if(indices.size()>1000000)throw std::invalid_argument("Too many Hamiltonian harmonics");
    const TorusFourierSeries empty{w,{}};
    HamiltonianFirstOrder out;out.generator=out.retained=empty;
    out.actionCorrection.assign(d,empty);out.angleCorrection.assign(d,empty);
    out.retainedActionGradient.assign(d,empty);out.meanFrequencyShift.resize(d);
    const FourierIndex zero(d,0);
    for(std::size_t j=0;j<d;++j)out.meanFrequencyShift[j]=coefficient(gradient[j],zero).real();
    for(const auto& k:indices) {
        const auto c=coefficient(h,k);const double nu=harmonicFrequency(k,w);
        if(std::abs(nu)<=tol) {
            out.retained.coefficients[k]=c;
            for(std::size_t j=0;j<d;++j)out.retainedActionGradient[j].coefficients[k]=coefficient(gradient[j],k);
            continue;
        }
        const auto chi=c/C(0,nu);out.generator.coefficients[k]=chi;
        out.maximumInverseDivisor=std::max(out.maximumInverseDivisor,1/std::abs(nu));
        for(std::size_t j=0;j<d;++j) {
            double dnu=0;for(std::size_t l=0;l<d;++l)dnu+=k[l]*D[l*d+j];
            out.actionCorrection[j].coefficients[k]=-C(0,k[j])*chi;
            out.angleCorrection[j].coefficients[k]=(coefficient(gradient[j],k)-c*(dnu/nu))/C(0,nu);
        }
    }
    meanSquare(out.generator);
    for(const auto& s:out.actionCorrection)meanSquare(s);
    for(const auto& s:out.angleCorrection)meanSquare(s);
    return out;
}
bool twistDiagnosticsAvailable() {
#ifdef DYNSYSLAB_MATH_NUMERICS
    return true;
#else
    return false;
#endif
}
TwistDiagnostics hamiltonianTwist(const std::vector<double>& w,
    const std::vector<double>& h,double tol) {
    checkHessian(w,h);
    if(!std::isfinite(tol)||tol<0||tol>=1)
        throw std::invalid_argument("Twist relative tolerance must lie in [0,1)");
#ifdef DYNSYSLAB_MATH_NUMERICS
    const auto d=w.size();std::vector<double> bordered((d+1)*(d+1),0);
    for(std::size_t i=0;i<d;++i) {
        for(std::size_t j=0;j<d;++j)bordered[i*(d+1)+j]=h[i*d+j];
        bordered[i*(d+1)+d]=bordered[d*(d+1)+i]=w[i];
    }
    TwistDiagnostics out;
    out.kolmogorovSingularValues=numerics::singularValueDecomposition(h,d,d).singularValues;
    out.isoenergeticSingularValues=numerics::singularValueDecomposition(bordered,d+1,d+1).singularValues;
    const auto ratio=[](const auto& s){return s.front()>0?s.back()/s.front():0.;};
    out.kolmogorovRelativeMinimum=ratio(out.kolmogorovSingularValues);
    out.isoenergeticRelativeMinimum=ratio(out.isoenergeticSingularValues);
    out.kolmogorovNumericallyNondegenerate=out.kolmogorovRelativeMinimum>tol;
    out.isoenergeticNumericallyNondegenerate=out.isoenergeticRelativeMinimum>tol;
    return out;
#else
    throw std::runtime_error("Twist diagnostics require DYNSYSLAB_ENABLE_LINALG");
#endif
}
}
