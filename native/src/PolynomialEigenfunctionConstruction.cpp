#include <dynsyslab/core/EigenfunctionFit.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
EigenfunctionFitResult polynomialEigenfunctions(
    std::size_t dimension,std::vector<std::size_t> powers,
    std::vector<std::complex<double>> values,Spectrum spectrum,
    std::vector<double> mean,std::vector<double> deviation,std::optional<FunctionDomain> domain) {
    const auto& eigenvalues=spectrum.values();
    if(!dimension||powers.empty()||powers.size()%dimension||eigenvalues.empty())
        throw std::invalid_argument("polynomial eigenfunctions need coordinates, monomials and eigenvalues");
    const auto features=powers.size()/dimension,modes=eigenvalues.size();
    if(modes>std::numeric_limits<std::size_t>::max()/features||values.size()!=modes*features)
        throw std::invalid_argument("polynomial coefficients must have shape modes x monomials");
    auto finite=[](std::complex<double> z){return std::isfinite(z.real())&&std::isfinite(z.imag());};
    if(!std::all_of(values.begin(),values.end(),finite)
        ||!std::all_of(eigenvalues.begin(),eigenvalues.end(),finite))
        throw std::invalid_argument("polynomial coefficients and eigenvalues must be finite");
    if(mean.empty()) mean.assign(dimension,0);
    if(deviation.empty()) deviation.assign(dimension,1);
    if(mean.size()!=dimension||deviation.size()!=dimension
        ||!std::all_of(mean.begin(),mean.end(),[](double x){return std::isfinite(x);})
        ||!std::all_of(deviation.begin(),deviation.end(),[](double x){return std::isfinite(x)&&x>0;}))
        throw std::invalid_argument("polynomial normalization requires finite centers and positive scales");
    std::set<std::vector<std::size_t>> unique;
    for(std::size_t f=0;f<features;++f) {
        const auto begin=powers.begin()+static_cast<std::ptrdiff_t>(f*dimension);
        if(!unique.emplace(begin,begin+static_cast<std::ptrdiff_t>(dimension)).second)
            throw std::invalid_argument("duplicate polynomial monomial");
    }
    for(std::size_t mode=0;mode<modes;++mode) {
        const auto begin=values.begin()+static_cast<std::ptrdiff_t>(mode*features);
        if(std::all_of(begin,begin+static_cast<std::ptrdiff_t>(features),
            [](auto value){return value==std::complex<double>{};}))
            throw std::invalid_argument("an identically zero polynomial is not an eigenfunction");
    }
    return {FunctionExpansion({dimension,std::move(powers),std::move(mean),std::move(deviation)},std::move(values),
                              domain?std::move(*domain):FunctionDomain::physicalState(dimension)),
            std::move(spectrum),{}};
}
} // namespace dynsyslab::core
