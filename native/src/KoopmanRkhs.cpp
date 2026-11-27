#include "KoopmanSpectralInternal.hpp"
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
namespace dynsyslab::core {
RkhsAdjointSpectrum rkhsAdjointSpectrum(const KoopmanSnapshotPairs& data,double bandwidth,
    std::size_t requestedRank,double tolerance,const std::vector<std::complex<double>>& grid) {
    using namespace spectral_detail;check(data,data.dimension);
    if(!std::isfinite(bandwidth)||bandwidth<=0||!std::isfinite(tolerance)||tolerance<=0)
        throw std::invalid_argument("RKHS kernel bandwidth and rank tolerance must be positive");
    if(!data.weights.empty())throw std::invalid_argument("RKHS section norms do not use quadrature weights");
    const auto n=data.samples,d=data.dimension;
    if(n>4096)throw std::invalid_argument("Dense RKHS analysis is limited to 4096 sections; select a smaller Hankel column count or subsample source pairs");
    const auto kernel=[&](const std::vector<double>& x,std::size_t i,const std::vector<double>& y,std::size_t j) {
        double r2=0;for(std::size_t k=0;k<d;++k){const double delta=x[i*d+k]-y[j*d+k];if(!std::isfinite(delta))throw std::invalid_argument("Nonfinite RKHS input");r2+=delta*delta;}
        return std::exp(-r2/(2*bandwidth*bandwidth));
    };
    std::vector<double> gram(n*n),cross(n*n),image(n*n);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
        gram[i*n+j]=kernel(data.x,i,data.x,j);
        cross[i*n+j]=kernel(data.x,i,data.y,j);
        image[i*n+j]=kernel(data.y,i,data.y,j);
    }
    const auto eig=numerics::symmetricEigensystem(gram,n);
    std::size_t r=0;
    while(r<n&&eig.eigenvalues[n-1-r]>tolerance*eig.eigenvalues.back()&&(!requestedRank||r<requestedRank))++r;
    if(!r)throw std::invalid_argument("RKHS Gram matrix has zero retained rank");
    std::vector<double> b(n*r);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<r;++j)b[i*r+j]=eig.eigenvectors[i*n+n-1-j]/std::sqrt(eig.eigenvalues[n-1-j]);
    const auto bt=transpose(b,n,r);
    const auto a=product(bt,product(cross,b,n,n,r),r,n,r);
    const auto l=product(bt,product(image,b,n,n,r),r,n,r);
    const auto spectrum=numerics::eigensystem(a,r);
    RkhsAdjointSpectrum out;out.dimension=d;out.nodesCount=n;out.rank=r;out.nodes=data.x;out.successors=data.y;out.bandwidth=bandwidth;
    out.eigenvalues=spectrum.eigenvalues;out.sectionCoefficients.resize(n*r);
    out.assumptions="Gaussian RKHS adjoint compression: K* k_x=k_T(x). Kernel norms replace empirical L2 quadrature. Boundedness/domain and dense-section limits are required for operator-level convergence; residuals are finite-space diagnostics. Adjoint eigenvectors are not automatically Koopman eigenfunctions or a forward predictor.";
    out.residuals.scope="finite-rkhs-section-span";out.residuals.definition="RKHS norm ||(K*-lambda)f||/||f|| using exact kernel Gram matrices, subject to floating-point error";
    out.residuals.pseudospectrum={out.residuals.definition,out.residuals.scope,grid,{},{},n};
    const auto residualMatrix=[&](std::complex<double> z) {
        if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Nonfinite RKHS spectral point");
        std::vector<std::complex<double>> h(r*r);
        for(std::size_t i=0;i<r;++i)for(std::size_t j=0;j<r;++j)
            h[i*r+j]=l[i*r+j]-z*a[j*r+i]-std::conj(z)*a[i*r+j]+(i==j?std::norm(z):0.);
        return h;
    };
    for(std::size_t j=0;j<r;++j) {
        auto h=residualMatrix(spectrum.eigenvalues[j]);std::complex<double> norm{};
        for(std::size_t i=0;i<r;++i)for(std::size_t k=0;k<r;++k)
            norm+=std::conj(spectrum.rightEigenvectors[i*r+j])*h[i*r+k]*spectrum.rightEigenvectors[k*r+j];
        out.residuals.eigenpair.push_back(std::sqrt(std::max(0.,norm.real())));
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<r;++k)out.sectionCoefficients[i*r+j]+=b[i*r+k]*spectrum.rightEigenvectors[k*r+j];
    }
    for(auto z:grid) {
        auto minimum=numerics::hermitianEigensystem(residualMatrix(z),r);
        out.residuals.pseudospectrum.minimum.push_back(std::sqrt(std::max(0.,minimum.eigenvalues.front())));
        for(std::size_t i=0;i<n;++i) {
            std::complex<double> value{};for(std::size_t k=0;k<r;++k)value+=b[i*r+k]*minimum.eigenvectors[k*r];
            out.residuals.pseudospectrum.minimizingCoefficients.push_back(value);
        }
    }
    return out;
}
}
