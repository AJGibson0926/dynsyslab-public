#include <dynsyslab/core/PeriodicFeatureIdentification.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double number(const J& v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
std::size_t integer(const J& v){const auto n=number(v);if(n<0||n>10000000||n!=std::floor(n))throw std::invalid_argument("Bounded integer required");return std::size_t(n);}
std::vector<double> numbers(const J& v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>10000000)throw std::invalid_argument("Bounded number array required");std::vector<double> out;for(auto& x:v.array())out.push_back(number(x));return out;}
J values(const std::vector<double>& v){A a;for(double x:v)a.emplace_back(x);return a;}
}
J evaluatePeriodicFeatureQuery(const J& q,const math::PeriodicAtlas& atlas){
    core::PeriodicHistoryFeatureMap map;const auto& f=q.at("features");
    for(double c:numbers(f.at("channels"))){if(c<0||c!=std::floor(c)||c>1000000)throw std::invalid_argument("Invalid feature channel");map.channels.push_back(std::size_t(c));}
    map.timeOffsets=numbers(f.at("time_offsets"));map.outputs=integer(f.at("outputs"));
    map.coefficients=numbers(f.at("coefficients"));map.noiseStandardDeviations=numbers(q.at("noise_standard_deviations"));
    map.whiten=q.at("whiten").boolean(true);map.relativeTolerance=q.at("relative_tolerance").number(1e-10);
    core::PeriodicFeatureIdentificationOptions options;
    if(!q.at("table_samples").isNull())options.tableSamples=integer(q.at("table_samples"));
    if(!q.at("phase_seeds").isNull())options.phaseSeeds=integer(q.at("phase_seeds"));
    if(!q.at("starts").isNull())options.starts=integer(q.at("starts"));
    if(!q.at("parameter_subdivisions").isNull())options.parameterSubdivisions=integer(q.at("parameter_subdivisions"));
    const auto& opt=q.at("optimization");
    if(!opt.at("maximum_iterations").isNull())options.optimization.maximumIterations=integer(opt.at("maximum_iterations"));
    options.optimization.finiteDifferenceStep=opt.at("finite_difference_step").number(1e-5);
    options.optimization.stepTolerance=opt.at("step_tolerance").number(1e-9);
    options.optimization.gradientTolerance=opt.at("gradient_tolerance").number(1e-8);
    const auto time=q.at("forecast_times").isNull()?std::vector<double>{}:numbers(q.at("forecast_times"));
    const auto result=core::identifyPeriodicHistoryFeatures(atlas,map,numbers(q.at("histories")),time,options);A rows;
    for(const auto& e:result.estimates)rows.emplace_back(O{{"parameters",values(e.fit.parameters)},
        {"squared_residual",e.fit.squaredResidual},{"converged",e.fit.converged},{"reason",e.fit.reason},
        {"iterations",double(e.fit.iterations)},{"evaluations",double(e.fit.evaluations)},
        {"angular_frequency",e.frequency},{"history_relative_error",e.historyRelativeError},
        {"current",values(e.current)},{"forecast",values(e.forecast)},{"jacobian_singular_values",values(e.jacobianSingularValues)}});
    return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-history-feature-identification"},
        {"inputs",double(result.inputs)},{"outputs",double(result.outputs)},{"noise_rank",double(result.rank)},
        {"singular_values",values(result.singularValues)},{"effective_coefficients",values(result.effectiveCoefficients)},
        {"whiten",map.whiten},{"forecast_times",values(time)},{"estimates",rows},{"source",q.at("source")},
        {"settings",O{{"table_samples",double(options.tableSamples)},{"phase_seeds",double(options.phaseSeeds)},
            {"parameter_subdivisions",double(options.parameterSubdivisions)},{"starts",double(options.starts)},
            {"relative_tolerance",map.relativeTolerance},{"optimization",opt},{"noise_standard_deviations",values(map.noiseStandardDeviations)}}},
        {"histories_sha256",core::sha256Text(core::json::dump(q.at("histories")))},
        {"features_sha256",core::sha256Text(core::json::dump(f))},{"atlas_sha256",core::sha256Text(core::json::dump(q.at("knots")))},
        {"scope","Past-only linear features; training-knot/phase coarse seeds and local bounded fit. No global uniqueness or control certificate. Atlas interpolation and feature truncation remain. Forecasts use inferred orbit frequency, not the frozen linear Koopman action."}};
}
}
