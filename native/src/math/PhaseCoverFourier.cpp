#include <dynsyslab/math/PhaseCoverFourier.hpp>
#include <dynsyslab/math/Arithmetic.hpp>
#include <cmath>
#include <stdexcept>

// Standard character calculus on circle-cover inverse limits:
// Burgos/Verjovsky, https://arxiv.org/abs/1603.05676, Section 2.
// Small-divisor homological equation: Poeschel, https://arxiv.org/abs/0908.2234.
// The finite solve is delegated to TorusFourier; the cover is not inferred.
namespace dynsyslab::math {
PhaseCoverForcingResult phaseCoverForcing(std::uint32_t prime,std::size_t depth,
    double base,const std::vector<double>& frequencies,
    const std::map<FourierIndex,std::complex<double>>& coefficients,double tolerance) {
    if(depth>32||frequencies.size()>7||!std::isfinite(base)||base<=0)
        throw std::invalid_argument("Finite positive base frequency, depth<=32 and <=7 drives required");
    const auto levels=phaseCoverTower(prime,std::vector<std::uint32_t>(depth,0));
    PhaseCoverForcingResult out;out.coverDegree=levels.back().modulus;
    const double degree=std::stod(out.coverDegree),frequency=base/degree;
    if(!std::isfinite(frequency)||frequency<=0)
        throw std::invalid_argument("Cover frequency underflows floating-point representation");
    out.forcing.angularFrequencies={frequency};
    out.forcing.angularFrequencies.insert(out.forcing.angularFrequencies.end(),frequencies.begin(),frequencies.end());
    out.forcing.coefficients=coefficients;
    out.response=solveCohomological(out.forcing,tolerance);
    for(const auto& [index,coefficient]:coefficients) {
        const auto ratio=canonicalRational(std::to_string(index.front())+"/"+out.coverDegree);
        const auto character=phaseCoverCharacter(prime,{},"0",ratio);
        const bool retained=out.response.retained.coefficients.count(index)!=0;
        out.modes.push_back({index,ratio,character.minimumDepth,
            harmonicFrequency(index,out.forcing.angularFrequencies),std::abs(coefficient),
            retained?0:std::abs(out.response.solution.coefficients.at(index)),retained});
    }
    return out;
}
}
