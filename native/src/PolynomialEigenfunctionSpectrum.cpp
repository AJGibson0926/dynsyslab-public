#include <dynsyslab/core/SINDy.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
void validTime(double dt) {
    if(!std::isfinite(dt)||dt<=0)
        throw std::invalid_argument("eigenvalue conversion requires a positive finite sampling interval");
}
bool finite(std::complex<double> value) {
    return std::isfinite(value.real())&&std::isfinite(value.imag());
}
} // namespace
KoopmanSINDyResult KoopmanSINDyResult::fromDiscreteCoefficients(
    std::size_t dimension,std::vector<std::size_t> exponents,
    std::vector<std::complex<double>> coefficients,
    std::vector<std::complex<double>> multipliers,double dt,
    std::vector<double> center,std::vector<double> scale,std::vector<int> branches) {
    validTime(dt);
    if(branches.empty()) branches.assign(multipliers.size(),0);
    if(branches.size()!=multipliers.size())
        throw std::invalid_argument("one logarithm branch is required per eigenvalue");
    const double twoPi=2*std::acos(-1.0);
    std::vector<std::complex<double>> rates;
    for(std::size_t i=0;i<multipliers.size();++i) {
        if(!finite(multipliers[i])||std::abs(multipliers[i])==0)
            throw std::invalid_argument("zero/nonfinite multiplier has no finite continuous eigenvalue");
        rates.push_back((std::log(multipliers[i])+std::complex<double>(0,twoPi*branches[i]))/dt);
    }
    auto result=fromCoefficients(dimension,std::move(exponents),std::move(coefficients),
        std::move(rates),std::move(center),std::move(scale));
    result.sampleTime=dt; result.eigenvalueOrigin="discrete-time";
    result.logarithmBranches=std::move(branches);
    return result;
}
std::vector<std::complex<double>> KoopmanSINDyResult::discreteEigenvalues(double dt) const {
    if(dt==0) dt=sampleTime;
    validTime(dt);
    std::vector<std::complex<double>> values;
    for(const auto rate:continuousEigenvalues) {
        if(!finite(rate)) throw std::invalid_argument("nonfinite continuous eigenvalue");
        const auto value=std::exp(rate*dt);
        if(!finite(value)) throw std::invalid_argument("discrete eigenvalue conversion overflow");
        values.push_back(value);
    }
    return values;
}
} // namespace dynsyslab::core
