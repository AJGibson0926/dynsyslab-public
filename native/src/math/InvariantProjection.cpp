#include <dynsyslab/math/InvariantProjection.hpp>
#include <cmath>
#include <stdexcept>
// Mean ergodic projection as conditional expectation, followed by L2
// contractivity / Pythagoras. Eisner et al., Operator Theoretic Aspects of
// Ergodic Theory, Remark 8.9: https://www.math.uni-leipzig.de/~eisner/book-EFHN.pdf
namespace dynsyslab::math {
InvariantReadoutFloor zeroMeanReadoutFloor(const std::vector<std::complex<double>>& m,
    const std::vector<double>& w,bool constant) {
    if(m.empty()||m.size()!=w.size()||m.size()>1000000)
        throw std::invalid_argument("Invariant means and weights must have matching bounded size");
    double mass=0;
    for(std::size_t i=0;i<w.size();++i) {
        if(!std::isfinite(w[i])||w[i]<0||!std::isfinite(m[i].real())||!std::isfinite(m[i].imag()))
            throw std::invalid_argument("Finite means and nonnegative finite weights required");
        mass+=w[i];
    }
    if(!std::isfinite(mass)||mass<=0)throw std::invalid_argument("Finite positive total weight required");
    InvariantReadoutFloor r;r.totalWeight=mass;
    if(constant)for(std::size_t i=0;i<w.size();++i)r.bestConstant+=(w[i]/mass)*m[i];
    double norm=0;
    for(std::size_t i=0;i<w.size();++i)norm=std::hypot(norm,std::sqrt(w[i]/mass)*std::abs(m[i]-r.bestConstant));
    r.meanSquaredFloor=norm*norm;
    if(!std::isfinite(r.meanSquaredFloor)||!std::isfinite(r.bestConstant.real())||!std::isfinite(r.bestConstant.imag()))
        throw std::overflow_error("Invariant readout floor overflow");
    return r;
}
}
