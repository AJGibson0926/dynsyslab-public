#include <dynsyslab/math/ResolventLeastSquares.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double eps=std::numeric_limits<double>::epsilon();
// a=U S V^H for rows>=cols.
ThinSvd tallSvd(ComplexMatrix a){
    const auto m=a.rows,n=a.cols;auto v=ComplexMatrix::identity(n);
    for(int sweep=0;sweep<60;++sweep){bool rotated=false;
        for(std::size_t p=0;p+1<n;++p)for(std::size_t q=p+1;q<n;++q){
            double alpha=0,beta=0;C gamma{};
            for(std::size_t i=0;i<m;++i){alpha+=std::norm(a(i,p));beta+=std::norm(a(i,q));gamma+=std::conj(a(i,p))*a(i,q);}
            const double g=std::abs(gamma);
            if(g==0||g<=eps*std::sqrt(alpha*beta))continue;
            rotated=true;const C phase=gamma/g;const double zeta=(beta-alpha)/(2*g);
            const double t=(zeta>=0?1.:-1.)/(std::abs(zeta)+std::sqrt(1+zeta*zeta)),c=1/std::sqrt(1+t*t),s=c*t;
            for(std::size_t i=0;i<m;++i){const C x=a(i,p),y=a(i,q)*std::conj(phase);a(i,p)=c*x-s*y;a(i,q)=s*x+c*y;}
            for(std::size_t i=0;i<n;++i){const C x=v(i,p),y=v(i,q)*std::conj(phase);v(i,p)=c*x-s*y;v(i,q)=s*x+c*y;}
        }
        if(!rotated)break;
    }
    std::vector<double> norms(n);
    for(std::size_t j=0;j<n;++j){double s=0;for(std::size_t i=0;i<m;++i)s+=std::norm(a(i,j));norms[j]=std::sqrt(s);}
    std::vector<std::size_t> order(n);std::iota(order.begin(),order.end(),std::size_t(0));
    std::stable_sort(order.begin(),order.end(),[&](std::size_t x,std::size_t y){return norms[x]>norms[y];});
    ThinSvd out;out.U=ComplexMatrix(m,n);out.V=ComplexMatrix(n,n);out.singularValues.resize(n);
    for(std::size_t j=0;j<n;++j){const auto o=order[j];const double s=norms[o];out.singularValues[j]=s;
        for(std::size_t i=0;i<m;++i)out.U(i,j)=s>0?a(i,o)/s:C{};
        for(std::size_t i=0;i<n;++i)out.V(i,j)=v(i,o);}
    return out;
}
}
ThinSvd thinSvd(const ComplexMatrix& a){
    if(a.rows==0||a.cols==0||a.data.size()!=a.rows*a.cols)throw std::invalid_argument("SVD needs a nonempty matrix");
    for(const auto& z:a.data)if(!std::isfinite(z.real())||!std::isfinite(z.imag()))throw std::invalid_argument("SVD needs a finite matrix");
    if(a.rows>=a.cols)return tallSvd(a);
    auto t=tallSvd(adjoint(a));ThinSvd out;out.U=std::move(t.V);out.V=std::move(t.U);out.singularValues=std::move(t.singularValues);
    return out;
}
LeastSquaresSolution leastSquares(const ComplexMatrix& a,const ComplexMatrix& b,double tolerance){
    if(a.rows!=b.rows)throw std::invalid_argument("Least squares row mismatch");
    if(!(tolerance>0))throw std::invalid_argument("Rank tolerance must be positive");
    const auto d=thinSvd(a);const auto k=d.singularValues.size();const double top=d.singularValues.front();
    LeastSquaresSolution out;out.x=ComplexMatrix(a.cols,b.cols);double smallest=top;
    const auto projected=product(adjoint(d.U),b);
    for(std::size_t r=0;r<k;++r){const double s=d.singularValues[r];
        if(!(s>tolerance*top)||s==0)continue;
        ++out.rank;smallest=s;
        for(std::size_t i=0;i<a.cols;++i)for(std::size_t j=0;j<b.cols;++j)out.x(i,j)+=d.V(i,r)*projected(r,j)/s;}
    out.condition=out.rank&&smallest>0?top/smallest:std::numeric_limits<double>::infinity();
    return out;
}
}
