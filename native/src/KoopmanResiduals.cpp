#include "KoopmanSpectralInternal.hpp"
#include <dynsyslab/core/Eigenpair.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <limits>
namespace dynsyslab::core {
using namespace spectral_detail;
SpectralResiduals koopmanResiduals(const SpectralKROM& model,const KoopmanSnapshotPairs& data,
    const std::vector<std::complex<double>>& grid,std::string scope) {
    model.validate();check(data,model.dictionary.dimension);
    const auto m=data.samples,n=model.dictionary.size(),r=model.rank;
    const auto w=weights(m,data.weights);
    const auto x=weighted(model.dictionary.evaluate(data.x),w,n),y=weighted(model.dictionary.evaluate(data.y),w,n);
    // One-step residual (Eigenpair.hpp; Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024):
    // sqrt(sum w |psi(y)-mu psi(x)|^2 / sum w |psi(x)|^2), finite quadrature evidence over the data, not a
    // global spectral certificate. The grid minimum is the same residual minimized over the trial space:
    // the sampled pseudospectrum (Trefethen & Embree 2005).
    SpectralResiduals out;out.scope=std::move(scope);out.definition=residual_definition::oneStep;
    auto& pseudo=out.pseudospectrum;pseudo.scope=out.scope;pseudo.definition=out.definition;pseudo.features=n;pseudo.grid=grid;
    for(std::size_t j=0;j<r;++j) {
        double denom=0,numer=0;
        for(std::size_t i=0;i<m;++i) {
            std::complex<double> a{},b{};
            for(std::size_t k=0;k<n;++k){a+=x[i*n+k]*model.eigenfunctions[k*r+j];b+=y[i*n+k]*model.eigenfunctions[k*r+j];}
            denom+=std::norm(a);numer+=std::norm(b-model.multipliers[j]*a);
        }
        if(denom<=1e-28)throw std::invalid_argument("Eigenfunction vanishes on residual evaluation samples");
        out.eigenpair.push_back(std::sqrt(numer/denom));
    }
    if(grid.empty())return out;
    auto q=product(x,model.basis,m,n,r);
    auto svd=numerics::singularValueDecomposition(q,m,r,0,1e-12);
    if(svd.rank!=r)throw std::invalid_argument("Residual grid cannot certify directions unobserved in validation data");
    auto b=product(model.basis,whitening(svd),n,r,r);
    const auto px=product(x,b,m,n,r),py=product(y,b,m,n,r);
    const auto cross=product(transpose(px,m,r),py,r,m,r);
    auto remainder=product(px,cross,m,r,r);
    for(std::size_t i=0;i<remainder.size();++i)remainder[i]=py[i]-remainder[i];
    const auto remainderGram=product(transpose(remainder,m,r),remainder,r,m,r);
    // Orthogonal residual split avoids cancellation in L-z*A*-conj(z)*A+|z|²I.
    // Compress once: grid cost depends on retained rank, not training sample count.
    for(auto z:grid) {
        if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("Nonfinite residual grid");
        std::vector<std::complex<double>> h(remainderGram.begin(),remainderGram.end()),row(r);
        for(std::size_t i=0;i<r;++i) {
            for(std::size_t j=0;j<r;++j)row[j]=cross[i*r+j]-(i==j?z:std::complex<double>{});
            for(std::size_t j=0;j<r;++j)for(std::size_t k=0;k<r;++k)h[j*r+k]+=std::conj(row[j])*row[k];
        }
        const auto eig=numerics::hermitianEigensystem(h,r);
        pseudo.minimum.push_back(std::sqrt(std::max(0.,eig.eigenvalues.front())));
        for(std::size_t j=0;j<n;++j) {
            std::complex<double> c{};
            for(std::size_t k=0;k<r;++k)c+=b[j*r+k]*eig.eigenvectors[k*r];
            pseudo.minimizingCoefficients.push_back(c);
        }
    }
    return out;
}
SpectralInvariance koopmanInvariance(const std::vector<double>& x,const std::vector<double>& y,
    std::size_t m,std::size_t n,double tolerance) {
    if(x.size()!=m*n||y.size()!=x.size())throw std::invalid_argument("Invariance data shapes differ");
    const auto u=numerics::singularValueDecomposition(x,m,n,0,tolerance);
    const auto v=numerics::singularValueDecomposition(y,m,n,0,tolerance);
    SpectralInvariance out;out.sourceRank=u.rank;out.successorRank=v.rank;
    if(!u.rank)throw std::invalid_argument("Principal angles need a nonzero source subspace");
    if(!v.rank)return out; // A zero successor subspace has no leakage or principal angles.
    const auto overlap=product(transpose(u.left,m,u.rank),v.left,u.rank,m,v.rank);
    const auto cosines=numerics::singularValueDecomposition(overlap,u.rank,v.rank).singularValues;
    for(auto c:cosines)out.principalAngles.push_back(std::acos(std::clamp(c,0.,1.)));
    auto projection=product(u.left,product(transpose(u.left,m,u.rank),y,u.rank,m,n),m,u.rank,n);
    double error=0,norm=0;
    for(std::size_t i=0;i<y.size();++i){error+=std::pow(y[i]-projection[i],2);norm+=y[i]*y[i];}
    out.projectionDefect=std::sqrt(error/std::max(norm,1e-300));
    return out;
}
std::vector<double> koopmanPageRank(const std::vector<double>& transition,std::size_t n,
    const std::vector<double>& personalization,double damping,double tolerance) {
    if(!n||transition.size()!=n*n||!std::isfinite(damping)||damping<0||damping>=1||
       !std::isfinite(tolerance)||tolerance<=0)throw std::invalid_argument("Invalid PageRank parameters");
    auto p=weights(n,personalization),rank=p;
    std::vector<double> a=transition;
    for(std::size_t i=0;i<n;++i) {
        double total=0;
        for(std::size_t j=0;j<n;++j){double v=a[i*n+j];if(!std::isfinite(v)||v<0)throw std::invalid_argument("PageRank requires a nonnegative transition graph");total+=v;}
        for(std::size_t j=0;j<n;++j)a[i*n+j]=total>0?a[i*n+j]/total:p[j];
    }
    for(std::size_t it=0;it<100000;++it) {
        auto next=product(rank,a,1,n,n);double error=0;
        for(std::size_t i=0;i<n;++i){next[i]=damping*next[i]+(1-damping)*p[i];error+=std::abs(next[i]-rank[i]);}
        rank=std::move(next);if(error<=tolerance)return rank;
    }
    throw std::runtime_error("PageRank iteration did not converge");
}
} // namespace dynsyslab::core
