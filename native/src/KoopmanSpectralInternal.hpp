#pragma once
#include <dynsyslab/core/KoopmanSpectralAnalysis.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::spectral_detail {
inline std::vector<double> product(const std::vector<double>& a,const std::vector<double>& b,
    std::size_t m,std::size_t k,std::size_t n) {
    return numerics::multiply(a,m,k,b,n);
}
inline std::vector<double> transpose(const std::vector<double>& a,std::size_t m,std::size_t n) {
    std::vector<double> b(a.size());
    for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<n;++j)b[j*m+i]=a[i*n+j];
    return b;
}
inline std::vector<double> weights(std::size_t m,const std::vector<double>& supplied,
                                   const std::string& policy="uniform") {
    if(!m||(!supplied.empty()&&supplied.size()!=m))throw std::invalid_argument("Snapshot weight shape mismatch");
    if(policy!="uniform"&&policy!="smooth-birkhoff")throw std::invalid_argument("Unknown weight policy");
    if(!supplied.empty()&&policy!="uniform")throw std::invalid_argument("Use explicit weights OR a weight policy");
    std::vector<double> w=supplied.empty()?std::vector<double>(m,1):supplied;
    double sum=0;
    for(std::size_t i=0;i<m;++i) {
        if(policy=="smooth-birkhoff") {
            const double t=(i+1.)/(m+1.);w[i]=std::exp(4.-1./(t*(1-t)));
        }
        if(!std::isfinite(w[i])||w[i]<0)throw std::invalid_argument("Weights must be nonnegative finite");
        sum+=w[i];
    }
    if(!(sum>0)||!std::isfinite(sum))throw std::invalid_argument("Zero or invalid total weight");
    for(auto& v:w)v/=sum;
    return w;
}
inline void check(const KoopmanSnapshotPairs& d,std::size_t dimension) {
    if(!dimension||!d.samples||d.dimension!=dimension||d.x.size()!=d.samples*d.dimension||d.y.size()!=d.x.size())
        throw std::invalid_argument("Snapshot pairs must be equally sized samples x dictionary coordinates");
}
inline std::vector<double> weighted(std::vector<double> a,const std::vector<double>& w,std::size_t n) {
    for(std::size_t i=0;i<w.size();++i)for(std::size_t j=0;j<n;++j)a[i*n+j]*=std::sqrt(w[i]);
    return a;
}
inline std::vector<double> whitening(const numerics::RealSvd& svd) {
    auto b=transpose(svd.rightTranspose,svd.rank,svd.columns);
    for(std::size_t i=0;i<svd.columns;++i)for(std::size_t j=0;j<svd.rank;++j)b[i*svd.rank+j]/=svd.singularValues[j];
    return b;
}
} // namespace dynsyslab::core::spectral_detail
