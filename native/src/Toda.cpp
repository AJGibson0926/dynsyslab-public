#include <dynsyslab/core/Toda.hpp>
#include <cmath>
#include <stdexcept>
// Flaschka Lax representation, with signs as in Teschl:
// https://diana.mat.univie.ac.at/~gerald/ftp/book-jac/toda.html
namespace dynsyslab::core {
TodaLaxPair todaLaxPair(const std::vector<double>& q,const std::vector<double>& p,
    double alpha,bool periodic) {
    const auto n=q.size();
    if(n<2||n>256||p.size()!=n||!std::isfinite(alpha)||alpha<=0)
        throw std::invalid_argument("Toda requires 2..256 sites, matched states, positive alpha");
    TodaLaxPair r;r.sites=n;r.matrix.assign(n*n,0);r.generator.assign(n*n,0);r.b.resize(n);
    for(std::size_t i=0;i<n;++i) {
        if(!std::isfinite(q[i])||!std::isfinite(p[i]))throw std::invalid_argument("Finite Toda state required");
        r.matrix[i*n+i]=r.b[i]=-p[i]/2;
    }
    const auto bonds=periodic?n:n-1;
    for(std::size_t i=0;i<bonds;++i) {
        const auto j=(i+1)%n;const double a=alpha/2*std::exp((q[i]-q[j])/2);
        if(!std::isfinite(a)||a==0)throw std::overflow_error("Toda bond outside floating-point range");
        r.a.push_back(a);
        // Add: two periodic sites have two oriented bonds joining the same pair.
        r.matrix[i*n+j]+=a;r.matrix[j*n+i]+=a;
        r.generator[i*n+j]+=a;r.generator[j*n+i]-=a;
    }return r;
}
}
