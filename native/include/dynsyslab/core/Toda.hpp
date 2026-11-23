#pragma once
#include <cstddef>
#include <vector>
namespace dynsyslab::core {
struct TodaLaxPair {
    std::size_t sites{};
    std::vector<double> matrix, generator; // row-major L and P: Ldot=[P,L]
    std::vector<double> a, b;
};
// H=sum p_j^2/2+alpha^2 sum exp(q_j-q_{j+1}), mass=1.
// a_j=alpha*exp((q_j-q_{j+1})/2)/2, b_j=-p_j/2.
// Periodic matrix is the Bloch multiplier +1 matrix. It alone does not
// encode the complete periodic spectral curve. Open chain omits end bond.
TodaLaxPair todaLaxPair(const std::vector<double>& positions,
    const std::vector<double>& momenta,double alpha=1,bool periodic=true);
}
