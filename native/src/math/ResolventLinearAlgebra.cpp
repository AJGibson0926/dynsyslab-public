#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using C=std::complex<double>;
const double eps=std::numeric_limits<double>::epsilon();
void requireSquare(const ComplexMatrix& a){
    if(a.rows!=a.cols||a.data.size()!=a.rows*a.cols)throw std::invalid_argument("Square matrix required");
}
}
ComplexMatrix ComplexMatrix::identity(std::size_t n){ComplexMatrix m(n,n);for(std::size_t i=0;i<n;++i)m(i,i)=1;return m;}
ComplexMatrix product(const ComplexMatrix& a,const ComplexMatrix& b){
    if(a.cols!=b.rows)throw std::invalid_argument("Matrix product dimension mismatch");
    ComplexMatrix out(a.rows,b.cols);
    for(std::size_t i=0;i<a.rows;++i)for(std::size_t k=0;k<a.cols;++k){const C x=a(i,k);if(x==C{})continue;
        for(std::size_t j=0;j<b.cols;++j)out(i,j)+=x*b(k,j);}
    return out;
}
ComplexVector product(const ComplexMatrix& a,const ComplexVector& x){
    if(a.cols!=x.size())throw std::invalid_argument("Matrix-vector dimension mismatch");
    ComplexVector out(a.rows);
    for(std::size_t i=0;i<a.rows;++i)for(std::size_t j=0;j<a.cols;++j)out[i]+=a(i,j)*x[j];
    return out;
}
ComplexMatrix adjoint(const ComplexMatrix& a){
    ComplexMatrix out(a.cols,a.rows);
    for(std::size_t i=0;i<a.rows;++i)for(std::size_t j=0;j<a.cols;++j)out(j,i)=std::conj(a(i,j));
    return out;
}
double frobeniusNorm(const ComplexMatrix& a){double s=0;for(const auto& x:a.data)s+=std::norm(x);return std::sqrt(s);}
double euclideanNorm(const ComplexVector& x){double s=0;for(const auto& v:x)s+=std::norm(v);return std::sqrt(s);}
ComplexMatrix shiftedMatrix(const ComplexMatrix& a,C s){
    requireSquare(a);ComplexMatrix out(a.rows,a.cols);
    for(std::size_t k=0;k<a.data.size();++k)out.data[k]=-a.data[k];
    for(std::size_t i=0;i<a.rows;++i)out(i,i)+=s;
    return out;
}
ComplexLu::ComplexLu(ComplexMatrix a):lu_(std::move(a)){
    requireSquare(lu_);const auto n=lu_.rows;pivots_.resize(n);
    for(std::size_t i=0;i<n;++i)pivots_[i]=i;
    for(std::size_t k=0;k<n;++k){
        std::size_t p=k;double best=std::abs(lu_(k,k));
        for(std::size_t i=k+1;i<n;++i){const double v=std::abs(lu_(i,k));if(v>best){best=v;p=i;}}
        if(best==0){singular_=true;continue;}
        if(p!=k){for(std::size_t j=0;j<n;++j)std::swap(lu_(k,j),lu_(p,j));std::swap(pivots_[k],pivots_[p]);}
        for(std::size_t i=k+1;i<n;++i){
            const C f=lu_(i,k)/lu_(k,k);lu_(i,k)=f;if(f==C{})continue;
            for(std::size_t j=k+1;j<n;++j)lu_(i,j)-=f*lu_(k,j);
        }
    }
}
ComplexVector ComplexLu::solve(ComplexVector b) const {
    const auto n=lu_.rows;if(b.size()!=n)throw std::invalid_argument("LU right-hand side size mismatch");
    if(singular_)throw std::domain_error("Singular matrix in resolvent solve");
    ComplexVector x(n);for(std::size_t i=0;i<n;++i)x[i]=b[pivots_[i]];
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j)x[i]-=lu_(i,j)*x[j];
    for(std::size_t i=n;i-->0;){for(std::size_t j=i+1;j<n;++j)x[i]-=lu_(i,j)*x[j];x[i]/=lu_(i,i);}
    return x;
}
ComplexMatrix ComplexLu::solve(const ComplexMatrix& b) const {
    ComplexMatrix out(b.rows,b.cols);ComplexVector c(b.rows);
    for(std::size_t j=0;j<b.cols;++j){
        for(std::size_t i=0;i<b.rows;++i)c[i]=b(i,j);
        const auto x=solve(c);for(std::size_t i=0;i<b.rows;++i)out(i,j)=x[i];
    }
    return out;
}
ComplexVector ComplexLu::solveAdjoint(const ComplexVector& b) const {
    const auto n=lu_.rows;if(b.size()!=n)throw std::invalid_argument("LU right-hand side size mismatch");
    if(singular_)throw std::domain_error("Singular matrix in resolvent solve");
    // P A = L U, so A^H P^T = U^H L^H: solve U^H w=b, L^H y=w, then x[p[i]]=y[i].
    ComplexVector y(b);
    for(std::size_t i=0;i<n;++i){for(std::size_t j=0;j<i;++j)y[i]-=std::conj(lu_(j,i))*y[j];y[i]/=std::conj(lu_(i,i));}
    for(std::size_t i=n;i-->0;)for(std::size_t j=i+1;j<n;++j)y[i]-=std::conj(lu_(j,i))*y[j];
    ComplexVector x(n);for(std::size_t i=0;i<n;++i)x[pivots_[i]]=y[i];
    return x;
}
ComplexMatrix ComplexLu::solveAdjoint(const ComplexMatrix& b) const {
    ComplexMatrix out(b.rows,b.cols);ComplexVector c(b.rows);
    for(std::size_t j=0;j<b.cols;++j){
        for(std::size_t i=0;i<b.rows;++i)c[i]=b(i,j);
        const auto x=solveAdjoint(c);for(std::size_t i=0;i<b.rows;++i)out(i,j)=x[i];
    }
    return out;
}
std::vector<double> complexSingularValues(ComplexMatrix a){
    if(a.rows<a.cols)a=adjoint(a);
    const auto m=a.rows,n=a.cols;
    for(int sweep=0;sweep<60;++sweep){bool rotated=false;
        for(std::size_t p=0;p+1<n;++p)for(std::size_t q=p+1;q<n;++q){
            double alpha=0,beta=0;C gamma{};
            for(std::size_t i=0;i<m;++i){alpha+=std::norm(a(i,p));beta+=std::norm(a(i,q));gamma+=std::conj(a(i,p))*a(i,q);}
            const double g=std::abs(gamma);
            if(g==0||g<=eps*std::sqrt(alpha*beta))continue;
            rotated=true;
            const C phase=gamma/g;const double zeta=(beta-alpha)/(2*g);
            const double t=(zeta>=0?1.:-1.)/(std::abs(zeta)+std::sqrt(1+zeta*zeta));
            const double c=1/std::sqrt(1+t*t),s=c*t;
            for(std::size_t i=0;i<m;++i){const C x=a(i,p),y=a(i,q)*std::conj(phase);a(i,p)=c*x-s*y;a(i,q)=s*x+c*y;}
        }
        if(!rotated)break;
    }
    std::vector<double> out(n);
    for(std::size_t j=0;j<n;++j){double s=0;for(std::size_t i=0;i<m;++i)s+=std::norm(a(i,j));out[j]=std::sqrt(s);}
    std::sort(out.rbegin(),out.rend());
    return out;
}
ComplexVector complexEigenvalues(ComplexMatrix h){
    requireSquare(h);const auto n=h.rows;ComplexVector values(n);
    if(n==0)return values;
    for(const auto& x:h.data)if(!std::isfinite(x.real())||!std::isfinite(x.imag()))throw std::invalid_argument("Finite matrix required");
    // Householder reduction to upper Hessenberg form.
    for(std::size_t k=0;k+2<n;++k){
        double norm=0;for(std::size_t i=k+1;i<n;++i)norm+=std::norm(h(i,k));norm=std::sqrt(norm);
        if(norm==0)continue;
        ComplexVector v(n);for(std::size_t i=k+1;i<n;++i)v[i]=h(i,k);
        const C head=v[k+1];const C alpha=head==C{}?C(-norm):-norm*head/std::abs(head);v[k+1]-=alpha;
        double vv=0;for(std::size_t i=k+1;i<n;++i)vv+=std::norm(v[i]);
        if(vv==0)continue;
        for(std::size_t j=k;j<n;++j){C d{};for(std::size_t i=k+1;i<n;++i)d+=std::conj(v[i])*h(i,j);d*=2./vv;
            for(std::size_t i=k+1;i<n;++i)h(i,j)-=v[i]*d;}
        for(std::size_t i=0;i<n;++i){C d{};for(std::size_t j=k+1;j<n;++j)d+=h(i,j)*v[j];d*=2./vv;
            for(std::size_t j=k+1;j<n;++j)h(i,j)-=d*std::conj(v[j]);}
        for(std::size_t i=k+2;i<n;++i)h(i,k)=0;
    }
    std::size_t hi=n-1;int iterations=0;const int limit=200*int(n)+200;
    while(true){
        // Deflate from the bottom.
        while(hi>0&&std::abs(h(hi,hi-1))<=eps*(std::abs(h(hi,hi))+std::abs(h(hi-1,hi-1)))){h(hi,hi-1)=0;--hi;}
        if(hi==0)break;
        std::size_t lo=hi;
        while(lo>0&&std::abs(h(lo,lo-1))>eps*(std::abs(h(lo,lo))+std::abs(h(lo-1,lo-1))))--lo;
        if(lo>0)h(lo,lo-1)=0;
        if(++iterations>limit)throw std::runtime_error("Complex QR iteration did not converge");
        const C a=h(hi-1,hi-1),b=h(hi-1,hi),c=h(hi,hi-1),d=h(hi,hi);
        const C half=(a-d)*.5,root=std::sqrt(half*half+b*c);
        const C denominator=std::abs(half+root)>=std::abs(half-root)?half+root:half-root;
        C shift=denominator==C{}?d:d-b*c/denominator; // Wilkinson: eigenvalue nearest d.
        if(iterations%11==0)shift+=std::abs(c); // exceptional shift
        for(std::size_t i=lo;i<=hi;++i)h(i,i)-=shift;
        std::vector<std::pair<double,C>> rotations;rotations.reserve(hi-lo);
        for(std::size_t k=lo;k<hi;++k){
            const C x=h(k,k),y=h(k+1,k);const double r=std::hypot(std::abs(x),std::abs(y));
            double cs=1;C sn{};
            if(r>0){if(x==C{}){cs=0;sn=std::conj(y)/std::abs(y);}else{cs=std::abs(x)/r;sn=(x/std::abs(x))*std::conj(y)/r;}}
            rotations.emplace_back(cs,sn);
            for(std::size_t j=k;j<=hi;++j){const C p=h(k,j),q=h(k+1,j);h(k,j)=cs*p+sn*q;h(k+1,j)=-std::conj(sn)*p+cs*q;}
        }
        for(std::size_t k=lo;k<hi;++k){const auto [cs,sn]=rotations[k-lo];
            for(std::size_t i=lo;i<=std::min(k+1,hi);++i){const C p=h(i,k),q=h(i,k+1);h(i,k)=cs*p+std::conj(sn)*q;h(i,k+1)=-sn*p+cs*q;}}
        for(std::size_t i=lo;i<=hi;++i)h(i,i)+=shift;
    }
    for(std::size_t i=0;i<n;++i)values[i]=h(i,i);
    return values;
}
}
