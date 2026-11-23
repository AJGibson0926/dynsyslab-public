#include <dynsyslab/core/OscillatorObservables.hpp>
#include <stdexcept>
#include <cmath>

namespace dynsyslab::core {
std::vector<Observable> phaseEnvelopeObservables(const std::string& name,
    std::size_t harmonic, OscillatorEnvelope envelope, const std::string& coordinates, double velocityScale) {
    if(name.empty() || !harmonic || harmonic > 1000000 || !std::isfinite(velocityScale) || velocityScale<=0)
        throw std::invalid_argument("Phase envelope requires a name and a positive bounded harmonic");
    const auto q=AnalyticExpression::state(2,0), v=velocityScale*AnalyticExpression::state(2,1);
    const auto weight=envelope==OscillatorEnvelope::sech
        ? v.unary("sech") : v*(-0.5*v.power(2)).unary("exp");
    const auto real=(double(harmonic)*q).unary("cos")*weight;
    const auto imag=(double(harmonic)*q).unary("sin")*weight;
    auto components=Observable::fromExpressions(
        {name+"_components",name+" (real, imaginary)","dimensionless",coordinates,
         {name+"_real",name+"_imag"},false},{real,imag});
    Observable complex({name,name,"dimensionless",coordinates,{},true},2,1,
        [components](const double* x) {
            const auto z=components.evaluate(x);
            return Observable::Values{{z[0].real(),z[1].real()}};
        },[components](const double* x) {
            const auto j=components.jacobian(x);
            return Observable::Values{{j[0].real(),j[2].real()},{j[1].real(),j[3].real()}};
        });
    return {std::move(complex),std::move(components),
        Observable::fromExpressions({name+"_real",name+" real","dimensionless",coordinates,{},false},{real}),
        Observable::fromExpressions({name+"_imag",name+" imaginary","dimensionless",coordinates,{},false},{imag})};
}
Observable gaussianSumObservable(const std::string& name,const std::string& coordinates) {
    const auto q=AnalyticExpression::state(2,0),v=AnalyticExpression::state(2,1);
    return Observable::fromExpressions({name,"Gaussian of position plus velocity","dimensionless",coordinates,{},false},
                                      {(-(q+v).power(2)).unary("exp")});
}
}
