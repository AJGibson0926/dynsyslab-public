#include "KoopmanSpectralInternal.hpp"
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
namespace dynsyslab::core {
namespace {
constexpr double pi=3.14159265358979323846;
using Complex=std::complex<double>;
void anglesValid(const std::vector<double>& angles) {
    if(angles.empty())throw std::invalid_argument("Spectral measure needs an angular grid");
    for(double t:angles)if(!std::isfinite(t))throw std::invalid_argument("Nonfinite angular grid");
}
}
SpectralMeasureEstimate correlationSpectralMeasure(const std::vector<Complex>& signal,
    std::size_t lags,const std::vector<double>& angles,bool acknowledgeStationarity) {
    if(!acknowledgeStationarity)throw std::invalid_argument("Correlation spectral measures require stationary/invariant sampling acknowledgement");
    if(signal.empty()||lags>=signal.size())throw std::invalid_argument("Correlation lags must be shorter than signal");
    anglesValid(angles);
    SpectralMeasureEstimate estimate;auto& band=estimate.band;auto& out=estimate.measure;
    band.name="correlation";band.method="biased-autocorrelation-fejer";band.grid=angles;
    band.assumptions="Stationary observable under invariant measure asserted by caller. Finite-data positive Fejer estimate; no separation of continuous and atomic mass.";
    out.observable="signal";out.band=band.name;
    for(std::size_t k=0;k<=lags;++k) {
        Complex c{};
        for(std::size_t i=0;i+k<signal.size();++i) {
            if(!std::isfinite(signal[i].real())||!std::isfinite(signal[i].imag()))throw std::invalid_argument("Nonfinite correlation signal");
            c+=std::conj(signal[i])*signal[i+k];
        }
        out.moments.push_back(c/double(signal.size()));
    }
    out.totalMass=out.moments.front().real();band.order=lags;
    for(double theta:angles) {
        double density=out.totalMass;
        for(std::size_t k=1;k<=lags;++k)
            density+=2*(1.-double(k)/(lags+1))*std::real(out.moments[k]*std::exp(Complex(0,-double(k)*theta)));
        out.density.push_back(density/(2*pi));
    }
    return estimate;
}
bool compressionIsUnitary(const KoopmanCompression& model,double tolerance) {
    const auto r=model.rank;
    if(!r||model.action.size()!=r*r)return false;
    for(std::size_t i=0;i<r;++i)for(std::size_t j=0;j<r;++j) {
        double gram=0;for(std::size_t k=0;k<r;++k)gram+=model.action[k*r+i]*model.action[k*r+j];
        if(std::abs(gram-(i==j?1.0:0.0))>tolerance)return false;
    }
    return true;
}

void setDictionaryObservableModes(KoopmanModel& krom,const std::string& observable,
    const std::vector<double>& g,std::size_t outputs) {
    // psi_j = sum_i E_ij phi_i on the frozen dictionary; g = E c exactly when g lies in their span.
    const auto& expansion=frozenDictionaryExpansion(krom);
    const auto n=expansion.dictionary().size(),r=krom.liftedDimension();
    if(!outputs||g.size()!=n*outputs)throw std::invalid_argument("Observable coefficients need one row per dictionary feature");
    const auto& e=expansion.coefficients();  // mode-major: pairs x features
    std::vector<Complex> matrix(n*r),rhs(g.begin(),g.end());
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<r;++j)matrix[i*r+j]=e[j*n+i];
    // E = basis V for a spectral fit (basis: the whitened trial basis, V: the action's eigenvectors), and E
    // inherits the whitening's condition number (sigma_max / sigma_min), so a direct solve on E loses
    // digits. Solve in two steps instead: z = basis^+ g on the basis's orthogonal columns (the residual
    // there is the test that g lies in the trial subspace), then V c = z with V well conditioned for the
    // normal compressions measure-preserving fits give (Colbrook, SIAM J. Numer. Anal. 61, 2023).
    std::vector<Complex> c;
    double error=0,norm=0;
    for(const auto value:rhs)norm+=std::norm(value);
    const auto* compression=krom.compression();
    if(compression&&compression->rank==r&&compression->basis.size()==n*r) {
        const std::vector<Complex> basis(compression->basis.begin(),compression->basis.end());
        const auto z=numerics::solveLeastSquares(basis,n,r,rhs,outputs);    // trial coordinates x outputs
        const auto v=numerics::solveLeastSquares(basis,n,r,matrix,r);       // E in trial coordinates
        double fit=0,scale=0;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<r;++j) {
            Complex value{};for(std::size_t k=0;k<r;++k)value+=basis[i*r+k]*v[k*r+j];
            fit+=std::norm(value-matrix[i*r+j]);scale+=std::norm(matrix[i*r+j]);
        }
        if(fit<=1e-20*std::max(scale,1e-300)) {  // E lies in the trial subspace: the two-step solve applies
            c=numerics::solveLeastSquares(v,r,r,z,outputs);
            for(std::size_t i=0;i<n;++i)for(std::size_t o=0;o<outputs;++o) {
                Complex value{};for(std::size_t k=0;k<r;++k)value+=basis[i*r+k]*z[k*outputs+o];
                error+=std::norm(value-rhs[i*outputs+o]);
            }
        }
    }
    if(c.empty()) {
        c=numerics::solveLeastSquares(matrix,n,r,rhs,outputs);  // pairs x outputs
        error=0;
        for(std::size_t i=0;i<n;++i)for(std::size_t o=0;o<outputs;++o) {
            Complex value{};for(std::size_t j=0;j<r;++j)value+=matrix[i*r+j]*c[j*outputs+o];
            error+=std::norm(value-rhs[i*outputs+o]);
        }
    }
    if(error>1e-16*std::max(norm,1e-300))throw std::invalid_argument("Observable is outside the span of the krom's eigenfunctions (relative residual "
        +std::to_string(std::sqrt(error/std::max(norm,1e-300)))+"); project it explicitly first");
    ReconstructionResult modes;modes.outputs=outputs;modes.modes=r;modes.matrix.resize(outputs*r);
    for(std::size_t o=0;o<outputs;++o)for(std::size_t j=0;j<r;++j)modes.matrix[o*r+j]=c[j*outputs+o];
    modes.relativeResidual=std::sqrt(error/std::max(norm,1e-300));
    krom.setReconstruction(observable,std::move(modes));
}

void riggedSpectralMeasure(KoopmanModel& krom,const std::vector<double>& g,const std::string& observable,
    std::size_t output,const std::string& bandName,const std::vector<double>& angles,double epsilon,std::size_t order) {
    if(!krom.compression())throw std::invalid_argument("Rigged DMD needs a krom with a compressed operator (a spectral fit's)");
    const auto& model=*krom.compression();
    model.validate();anglesValid(angles);
    if(observable.empty()||bandName.empty())throw std::invalid_argument("A Rigged DMD measure names its observable and its band");
    if(krom.spectrum().measure(bandName,observable,output))
        throw std::invalid_argument("This krom already has a measure of "+observable+" on "+bandName);
    // The resolvent measure lives on the unit circle: the compression must be unitary (A^T A = I), as
    // measure-preserving EDMD's polar factor is (Colbrook 2023). The operator decides, not its method name.
    if(!compressionIsUnitary(model))
        throw std::invalid_argument("Rigged DMD needs a unitary (measure-preserving) compression, such as mpEDMD's");
    if(!std::isfinite(epsilon)||epsilon<=0||!order||order>8||g.size()!=model.dictionary.size())
        throw std::invalid_argument("Rigged measure needs positive bandwidth, order 1..8 and dictionary observable coefficients");
    if(const auto* existing=krom.spectrum().band(bandName))  // checked before the krom changes
        if(existing->method!="rigged-resolvent"||existing->grid!=angles||existing->bandwidth!=epsilon||existing->order!=order)
            throw std::invalid_argument("Spectral band "+bandName+" exists with another grid or kernel");
    const auto n=model.dictionary.size(),r=model.rank;
    const auto realC=numerics::solveLeastSquares(model.basis,n,r,g,1);
    auto projected=spectral_detail::product(model.basis,realC,n,r,1);
    double error=0,norm=0;
    for(std::size_t i=0;i<n;++i){error+=std::pow(projected[i]-g[i],2);norm+=g[i]*g[i];}
    if(error>1e-16*std::max(norm,1e-300))throw std::invalid_argument("Observable is outside retained dictionary subspace; project it explicitly first");
    // The observable's Koopman modes: the atoms of its measure are |c_j|^2 on the unit-norm eigenfunctions
    // (Colbrook 2023), derived from them by Spectrum::atoms rather than stored with the measure.
    if(!krom.spectrum().koopmanModes(observable)) {
        if(output)throw std::invalid_argument("Observable "+observable+" has no Koopman modes; set them for all its outputs first");
        setDictionaryObservableModes(krom,observable,g,1);
    }
    const math::RationalCircleKernel kernel(epsilon,order);
    const auto& poles=kernel.poles();const auto& residues=kernel.residues();
    std::vector<Complex> c(realC.begin(),realC.end());
    SpectralBand band;band.name=bandName;band.method="rigged-resolvent";band.grid=angles;
    band.bandwidth=epsilon;band.order=order;
    band.assumptions="Regularized spectral measure of the finite unitary compression. Infinite-dimensional convergence additionally requires consistent invariant quadrature, dictionary limits, and smoothing limits. Higher-order kernels need not be positive.";
    SpectralMeasure measure;measure.observable=observable;measure.output=output;measure.band=bandName;
    for(auto value:c)measure.totalMass+=std::norm(value);
    std::vector<Complex> packets;  // grid x features: one generalized eigenfunction per angle
    // Reuse one Schur factorization for every pole/angle (no diagonal
    // approximation). This preserves repeated eigenvalues and avoids a dense
    // QR factorization for each resolvent. Colbrook/Drysdale/Horning (2024),
    // https://arxiv.org/abs/2405.00782; standard unitary Schur similarity.
    const auto schur=numerics::complexSchur(model.action,r);
    const auto& u=schur.unitary;const auto& t=schur.triangular;
    std::vector<Complex> cs(r),tc(r),bu(n*r);
    for(std::size_t i=0;i<r;++i)for(std::size_t j=0;j<r;++j)cs[i]+=std::conj(u[j*r+i])*c[j];
    for(std::size_t i=0;i<r;++i)for(std::size_t j=i;j<r;++j)tc[i]+=t[i*r+j]*cs[j];
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<r;++k)for(std::size_t j=0;j<r;++j)
        bu[i*r+j]+=model.basis[i*r+k]*u[k*r+j];
    const auto resolvent=[&](Complex z) {
        std::vector<Complex> b(r);
        for(std::size_t i=r;i-- >0;) {
            b[i]=tc[i]+z*cs[i];
            for(std::size_t j=i+1;j<r;++j)b[i]-=t[i*r+j]*b[j];
            b[i]/=t[i*r+i]-z;
            if(!std::isfinite(b[i].real())||!std::isfinite(b[i].imag()))
                throw std::runtime_error("Singular/nonfinite spectral resolvent");
        }
        return b;
    };
    for(double theta:angles) {
        std::vector<Complex> packet(r);double density=0;
        for(std::size_t j=0;j<order;++j) {
            const auto plus=resolvent(std::exp(Complex(0,theta)-Complex(0,epsilon)*poles[j]));
            const auto minus=resolvent(std::exp(Complex(0,theta)-Complex(0,epsilon)*std::conj(poles[j])));
            Complex inner{};
            for(std::size_t k=0;k<r;++k) {
                inner+=std::conj(cs[k])*plus[k];
                packet[k]-=(residues[j]*plus[k]-std::conj(residues[j])*minus[k])/(4*pi);
            }
            density-=std::real(residues[j]*inner)/(2*pi);
        }
        measure.density.push_back(density);
        for(std::size_t i=0;i<n;++i){Complex value{};for(std::size_t k=0;k<r;++k)value+=bu[i*r+k]*packet[k];packets.push_back(value);}
    }
    // The band (matched when it exists) and the observable's measure on it join the spectrum, the
    // generalized eigenfunctions as an expansion on the frozen dictionary.
    const auto& expansion=frozenDictionaryExpansion(krom);
    measure.packets=std::make_shared<const FunctionExpansion>(
        std::make_shared<const ExpansionDictionary>(expansion.dictionary()),std::move(packets),expansion.coordinates());
    auto spectrum=krom.spectrum();
    spectrum.addBand(std::move(band));
    spectrum.addMeasure(std::move(measure));
    krom.setSpectralEvidence(spectrum);
}
} // namespace dynsyslab::core
