#include "KoopmanSpectralInternal.hpp"
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/core/FileDigest.hpp>
namespace dynsyslab::core {
namespace {
constexpr double pi=3.14159265358979323846;
using Complex=std::complex<double>;
void anglesValid(const std::vector<double>& angles) {
    if(angles.empty())throw std::invalid_argument("Spectral measure needs an angular grid");
    for(double t:angles)if(!std::isfinite(t))throw std::invalid_argument("Nonfinite angular grid");
}
}
SpectralMeasure correlationSpectralMeasure(const std::vector<Complex>& signal,
    std::size_t lags,const std::vector<double>& angles,bool acknowledgeStationarity) {
    if(!acknowledgeStationarity)throw std::invalid_argument("Correlation spectral measures require stationary/invariant sampling acknowledgement");
    if(signal.empty()||lags>=signal.size())throw std::invalid_argument("Correlation lags must be shorter than signal");
    anglesValid(angles);SpectralMeasure out;out.method="biased-autocorrelation-fejer";out.grid=angles;
    out.assumptions="Stationary observable under invariant measure asserted by caller. Finite-data positive Fejer estimate; no separation of continuous and atomic mass.";
    for(std::size_t k=0;k<=lags;++k) {
        Complex c{};
        for(std::size_t i=0;i+k<signal.size();++i) {
            if(!std::isfinite(signal[i].real())||!std::isfinite(signal[i].imag()))throw std::invalid_argument("Nonfinite correlation signal");
            c+=std::conj(signal[i])*signal[i+k];
        }
        out.moments.push_back(c/double(signal.size()));
    }
    out.totalMass=out.moments.front().real();out.order=lags;
    for(double theta:angles) {
        double density=out.totalMass;
        for(std::size_t k=1;k<=lags;++k)
            density+=2*(1.-double(k)/(lags+1))*std::real(out.moments[k]*std::exp(Complex(0,-double(k)*theta)));
        out.density.push_back(density/(2*pi));
    }
    return out;
}
void riggedSpectralMeasure(SpectralKROM& model,const std::vector<double>& g,
    const std::string& name,const std::vector<double>& angles,double epsilon,std::size_t order) {
    model.validate();anglesValid(angles);
    if(model.method!="mpedmd")throw std::invalid_argument("Unit-circle resolvent measure requires an mpEDMD unitary compression");
    if(!std::isfinite(epsilon)||epsilon<=0||!order||order>8||g.size()!=model.dictionary.size())
        throw std::invalid_argument("Rigged measure needs positive bandwidth, order 1..8 and dictionary observable coefficients");
    const auto n=model.dictionary.size(),r=model.rank;
    const auto realC=numerics::solveLeastSquares(model.basis,n,r,g,1);
    auto projected=spectral_detail::product(model.basis,realC,n,r,1);
    double error=0,norm=0;
    for(std::size_t i=0;i<n;++i){error+=std::pow(projected[i]-g[i],2);norm+=g[i]*g[i];}
    if(error>1e-16*std::max(norm,1e-300))throw std::invalid_argument("Observable is outside retained dictionary subspace; project it explicitly first");
    std::vector<Complex> poles(order),vandermonde(order*order),rhs(order),c(realC.begin(),realC.end());
    for(std::size_t j=0;j<order;++j) {
        poles[j]=Complex(2.*(j+1)/(order+1)-1,1);
        for(std::size_t k=0;k<order;++k)vandermonde[k*order+j]=std::pow(poles[j],double(k));
    }
    rhs[0]=1;const auto residues=numerics::solveLeastSquares(vandermonde,order,order,rhs,1);
    SpectralMeasure measure;measure.name=name;measure.method="rigged-resolvent";measure.grid=angles;
    measure.bandwidth=epsilon;measure.order=order;
    measure.assumptions="Regularized spectral measure of the finite unitary compression. Infinite-dimensional convergence additionally requires consistent invariant quadrature, dictionary limits, and smoothing limits. Higher-order kernels need not be positive.";
    for(auto value:c)measure.totalMass+=std::norm(value);
    SpectralWavePackets packets;packets.observable=name;packets.grid=angles;packets.features=n;
    packets.dictionaryIdentity=sha256Text(json::dump(spectralDictionaryJson(model.dictionary)));
    packets.bandwidth=epsilon;packets.order=order;packets.assumptions=measure.assumptions;
    const auto resolvent=[&](Complex z) {
        std::vector<Complex> a(r*r),b(r);
        for(std::size_t i=0;i<r;++i)for(std::size_t j=0;j<r;++j) {
            a[i*r+j]=model.action[i*r+j]-(i==j?z:Complex{});
            b[i]+=(model.action[i*r+j]+(i==j?z:Complex{}))*c[j];
        }
        return numerics::solveLeastSquares(a,r,r,b,1);
    };
    for(double theta:angles) {
        std::vector<Complex> packet(r);double density=0;
        for(std::size_t j=0;j<order;++j) {
            const auto plus=resolvent(std::exp(Complex(0,theta)-Complex(0,epsilon)*poles[j]));
            const auto minus=resolvent(std::exp(Complex(0,theta)-Complex(0,epsilon)*std::conj(poles[j])));
            Complex inner{};
            for(std::size_t k=0;k<r;++k) {
                inner+=std::conj(c[k])*plus[k];
                packet[k]-=(residues[j]*plus[k]-std::conj(residues[j])*minus[k])/(4*pi);
            }
            density-=std::real(residues[j]*inner)/(2*pi);
        }
        measure.density.push_back(density);
        for(std::size_t i=0;i<n;++i){Complex value{};for(std::size_t k=0;k<r;++k)value+=model.basis[i*r+k]*packet[k];packets.coefficients.push_back(value);}
    }
    model.spectral.measures.push_back(std::move(measure));model.spectral.wavePackets.push_back(std::move(packets));
}
} // namespace dynsyslab::core
