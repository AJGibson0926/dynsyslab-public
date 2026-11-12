#include <dynsyslab/core/Pbh.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
PbhResult pbhObservability(const std::vector<double>& A,const std::vector<double>& C,
    std::size_t n,std::size_t p,bool discrete,double tolerance) {
    if(!n||!p||A.size()!=n*n||C.size()!=p*n||!(tolerance>0)||!std::isfinite(tolerance))
        throw std::invalid_argument("invalid PBH dimensions or tolerance");
    for(const auto* v:{&A,&C})for(double x:*v)if(!std::isfinite(x))
        throw std::invalid_argument("PBH inputs must be finite");
    PbhResult result;result.discreteTime=discrete;
    result.eigenvalues=numerics::eigensystem(A,n,false).eigenvalues;
    for(const auto lambda:result.eigenvalues) {
        std::vector<std::complex<double>> matrix((n+p)*n);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
            matrix[i*n+j]=(i==j?lambda:std::complex<double>{})-A[i*n+j];
        std::copy(C.begin(),C.end(),matrix.begin()+n*n);
        const auto s=numerics::singularValues(matrix,n+p,n);
        const double margin=s.front()>0?s.back()/s.front():0;
        result.relativeMargins.push_back(margin);
        if(margin<=tolerance) {
            result.fullRank=false;
            const bool stable=discrete?std::abs(lambda)<1-tolerance:lambda.real()<-tolerance;
            if(!stable)result.stableHiddenModes=false;
        }
    }
    return result;
}
PbhResult pbhControllability(const std::vector<double>& A,const std::vector<double>& B,
    std::size_t n,std::size_t m,bool discrete,double tolerance) {
    if(A.size()!=n*n||B.size()!=n*m)throw std::invalid_argument("invalid PBH input shape");
    std::vector<double> At(A.size()),Bt(B.size());
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)At[i*n+j]=A[j*n+i];
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<m;++j)Bt[j*n+i]=B[i*m+j];
    return pbhObservability(At,Bt,n,m,discrete,tolerance);
}
} // namespace dynsyslab::core
