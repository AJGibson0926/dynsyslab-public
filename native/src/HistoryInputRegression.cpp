#include <dynsyslab/core/HistoryInputRegression.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
std::vector<double> features(std::size_t degree,std::complex<double> z) {
    if(degree==0) return {1.0};
    if(degree==1) return {1.0,z.real(),z.imag()};
    throw std::invalid_argument("history input-map degree must be 0 or 1");
}
}
std::complex<double> HistoryInputMap::evaluate(std::complex<double> z) const {
    const auto f=features(degree,z);
    if(coefficients.size()!=f.size()) throw std::invalid_argument("input map shape mismatch");
    std::complex<double> result{};
    for(std::size_t j=0;j<f.size();++j) result+=coefficients[j]*f[j];
    return result;
}
ApproximateControlMap HistoryInputMap::canonicalRepresentation() const {
    const auto lambda=std::log(multiplier)/sampleTime;
    const auto integral=std::abs(lambda)>1e-14 ? (multiplier-1.0)/lambda
        : std::complex<double>(sampleTime,0);
    if(degree==0) return ConstantControlMap{{coefficients.at(0)/integral},1,1,
        "echo-history sampled residual; measurement-only",{coefficients.at(0)},sampleTime};
    const std::complex<double> i(0,1);
    BilinearControlMap map;
    map.modes=1; map.inputs=1; map.provenance="echo-history sampled residual; frozen-map ZOH conversion";
    map.offsets={coefficients.at(0)/integral};
    map.operators={(coefficients.at(1)-i*coefficients.at(2))/(2.0*integral)};
    map.conjugateOperators={(coefficients.at(1)+i*coefficients.at(2))/(2.0*integral)};
    return map;
}
HistoryInputMap fitHistoryInputMap(const std::vector<HistoryInputSample>& samples,
    std::complex<double> multiplier,double sampleTime,std::size_t degree,
    std::array<double,2> coordinateScale,double relativeRidge) {
    if(samples.empty() || !(sampleTime>0) || !std::isfinite(sampleTime)
        || relativeRidge<0 || !std::isfinite(relativeRidge))
        throw std::invalid_argument("invalid history input regression options");
    const auto columns=features(degree,{}).size();
    std::vector<double> design(samples.size()*columns),rhs(samples.size()*2),scales(columns);
    for(std::size_t row=0;row<samples.size();++row) {
        const auto& s=samples[row]; const auto f=features(degree,s.psi);
        for(std::size_t j=0;j<columns;++j) {
            const double v=s.input*f[j]; design[row*columns+j]=v; scales[j]+=v*v;
        }
        const auto q=s.nextPsi-multiplier*s.psi;
        rhs[row*2]=q.real(); rhs[row*2+1]=q.imag();
    }
    for(std::size_t j=0;j<columns;++j) {
        scales[j]=std::sqrt(scales[j]);
        if(!(scales[j]>1e-20) || !std::isfinite(scales[j]))
            throw std::invalid_argument("unexcited or nonfinite input-map feature");
        for(std::size_t row=0;row<samples.size();++row) design[row*columns+j]/=scales[j];
    }
    numerics::LeastSquaresOptions options;
    options.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    options.ridge=relativeRidge;
    const auto fit=numerics::solveLeastSquares(design,samples.size(),columns,rhs,2,options);
    std::vector<std::complex<double>> complexDesign(design.begin(),design.end());
    const auto singular=numerics::singularValues(complexDesign,samples.size(),columns);
    HistoryInputMap out;
    out.degree=degree; out.coordinateScale=coordinateScale; out.multiplier=multiplier;
    out.sampleTime=sampleTime; out.samples=samples.size();
    out.rank=std::count_if(singular.begin(),singular.end(),[&](double s){return s>singular[0]*1e-10;});
    out.conditionNumber=singular[0]/std::max(singular.back(),1e-300);
    for(std::size_t j=0;j<columns;++j)
        out.coefficients.emplace_back(fit[j*2]/scales[j],fit[j*2+1]/scales[j]);
    return out;
}
HistoryInputMetrics scoreHistoryInputMap(const HistoryInputMap& map,
    const std::vector<HistoryInputSample>& samples) {
    if(samples.empty()) throw std::invalid_argument("empty history input validation set");
    double sse=0,base=0;
    for(const auto& s:samples) {
        const auto q=s.nextPsi-map.multiplier*s.psi;
        base+=std::norm(q); sse+=std::norm(q-map.evaluate(s.psi)*s.input);
    }
    return {std::sqrt(sse/std::max(base,1e-300)),std::sqrt(sse/samples.size()),
        std::sqrt(base/samples.size())};
}
} // namespace dynsyslab::core
