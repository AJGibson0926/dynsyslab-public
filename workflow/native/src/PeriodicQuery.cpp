#include <dynsyslab/workflow/PeriodicQuery.hpp>
#include <dynsyslab/core/PeriodicEstimation.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
core::json::Value evaluatePeriodicSnapshotQuery(const core::json::Value&);
core::json::Value evaluatePeriodicSpectralQuery(const core::json::Value&);
core::json::Value evaluatePeriodicControlQuery(const core::json::Value&,const math::PeriodicAtlas&);
core::json::Value evaluatePeriodicIdentificationQuery(const core::json::Value&,const math::PeriodicAtlas&);
core::json::Value evaluatePeriodicFeatureQuery(const core::json::Value&,const math::PeriodicAtlas&);
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J& v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite JSON number required");return v.number();}
std::size_t integer(const J& v,std::size_t maximum=1000000){const double x=num(v);if(x<0||x>maximum||x!=std::floor(x))throw std::invalid_argument("Integer outside allowed range");return std::size_t(x);}
const A& array(const J& v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>10000000)throw std::invalid_argument("Bounded JSON array required");return v.array();}
std::vector<double> reals(const J& q){std::vector<double> out;for(const auto& x:array(q))out.push_back(num(x));return out;}
std::vector<std::size_t> indices(const J& q){std::vector<std::size_t> out;for(const auto& x:array(q))out.push_back(integer(x));return out;}
J jsonVector(const std::vector<double>& x){A a;for(double v:x)a.emplace_back(v);return a;}
math::FourierSeries series(const J& q){math::FourierSeries s{num(q.at("angular_frequency")),{}};
    for(const auto& row:array(q.at("coefficients"))){const auto& c=array(row);if(c.size()!=3)throw std::invalid_argument("Need [harmonic,real,imaginary]");
        const double n=num(c[0]);if(std::abs(n)>1000000||n!=std::floor(n))throw std::invalid_argument("Invalid harmonic");
        if(!s.coefficients.emplace(int(n),std::complex<double>(num(c[1]),num(c[2]))).second)throw std::invalid_argument("Duplicate harmonic");}
    return s;}
}
core::json::Value evaluatePeriodicQuery(const core::json::Value& q){
    if(q.at("schema").string()!="dynsyslab.periodic-query/v1")throw std::invalid_argument("Expected periodic-query/v1");
    const auto operation=q.at("operation").string();
    if(operation=="snapshot_pairs")return evaluatePeriodicSnapshotQuery(q);
    if(operation=="spectral_sweep"||operation=="spectral_compare"||operation=="quadrature_audit"||operation=="delay_fit"||operation=="delay_forecast_audit"||operation=="packet_features")return evaluatePeriodicSpectralQuery(q);
    std::vector<math::PeriodicAtlasKnot> knots;
    for(const auto& k:array(q.at("knots"))){math::PeriodicAtlasKnot knot;knot.parameter=num(k.at("parameter"));
        for(const auto& c:array(k.at("channels")))knot.channels.push_back(series(c));knots.push_back(std::move(knot));}
    math::PeriodicAtlas atlas(std::move(knots));core::PeriodicFilterOptions options;
    if(operation=="feature_identify")return evaluatePeriodicFeatureQuery(q,atlas);
    if(operation=="controlled_predict")return evaluatePeriodicControlQuery(q,atlas);
    if(operation=="controlled_identify"||operation=="controlled_sequential_identify"||operation=="controlled_replay"||operation=="controlled_feedback")return evaluatePeriodicIdentificationQuery(q,atlas);
    if(operation=="predict"){
        const auto tableSize=q.at("table_samples").isNull()?4096:integer(q.at("table_samples"));
        const math::PeriodicAtlasTable table(atlas,tableSize);
        const double p=num(q.at("parameter")),phase=num(q.at("initial_phase"));
        const double t0=q.at("initial_time").isNull()?0:num(q.at("initial_time"));
        const double frozen=q.at("frozen_frequency").isNull()?0:num(q.at("frozen_frequency"));
        if(frozen<0)throw std::invalid_argument("Nonnegative frozen frequency required");
        const double omega=frozen>0?frozen:table.frequency(p);const auto time=reals(q.at("time"));
        const auto channels=atlas.knots()[0].channels.size();
        if(time.empty()||time.size()>1000000||channels>10000000/time.size())throw std::invalid_argument("Prediction output exceeds limits");
        A values;for(double t:time)for(auto c:table.evaluate(p,phase+omega*(t-t0)))values.emplace_back(A{c.real(),c.imag()});
        return O{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-atlas-prediction"},
            {"predictions",values},{"channels",double(channels)},{"angular_frequency",omega},{"time",jsonVector(time)}};
    }
    if(!q.at("parameters").isNull())options.parameters=reals(q.at("parameters"));
    options.measurementChannels=indices(q.at("measurement_channels"));
    if(!q.at("output_channels").isNull())options.outputChannels=indices(q.at("output_channels"));
    options.measurementVariances=reals(q.at("measurement_variances"));
    if(!q.at("parameter_prior").isNull())options.parameterPrior=reals(q.at("parameter_prior"));
    if(!q.at("phase_count").isNull())options.phaseCount=integer(q.at("phase_count"));
    if(!q.at("table_samples").isNull())options.tableSamples=integer(q.at("table_samples"));
    if(!q.at("frozen_frequency").isNull())options.frozenFrequency=num(q.at("frozen_frequency"));
    if(operation=="identify"){
        std::vector<std::vector<double>> starts;for(const auto& start:array(q.at("initial_candidates")))starts.push_back(reals(start));
        numerics::NonlinearLeastSquaresOptions optimizer;
        if(!q.at("maximum_iterations").isNull())optimizer.maximumIterations=integer(q.at("maximum_iterations"),10000);
        const auto time=reals(q.at("time"));
        const auto fit=core::identifyPeriodicOrbit(atlas,time,reals(q.at("measurements")),starts,options,optimizer);
        O result{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-window-identification"},
            {"parameters",jsonVector(fit.parameters)},{"squared_residual",fit.squaredResidual},
            {"iterations",double(fit.iterations)},{"evaluations",double(fit.evaluations)},
            {"converged",fit.converged},{"reason",fit.reason}};
        if(!q.at("prediction_times").isNull()){
            const math::PeriodicAtlasTable table(atlas,options.tableSamples);
            const double omega=options.frozenFrequency>0?options.frozenFrequency:table.frequency(fit.parameters[0]);
            A values;for(double t:reals(q.at("prediction_times")))for(auto c:table.evaluate(fit.parameters[0],fit.parameters[1]+omega*(t-time.front())))values.emplace_back(A{c.real(),c.imag()});
            result["predictions"]=values;result["channels"]=double(atlas.knots()[0].channels.size());
        }
        return result;
    }
    if(operation!="filter")throw std::invalid_argument("Unknown periodic analysis operation");
    const auto out=core::filterPeriodicAtlas(atlas,reals(q.at("time")),reals(q.at("measurements")),options);
    O result{{"schema","dynsyslab.periodic-result/v1"},{"method","periodic-ensemble-bayes"},
        {"outputs",jsonVector(out.outputs)},{"parameter_mean",jsonVector(out.parameterMean)},
        {"parameter_variance",jsonVector(out.parameterVariance)},{"effective_sample_size",jsonVector(out.effectiveSampleSize)},
        {"log_evidence",jsonVector(out.logEvidence)},{"final_weights",jsonVector(out.finalWeights)},
        {"output_dimension",double(options.outputChannels.size())},
        {"assumptions","Unforced constant orbit parameter, consistent phase origins, Gaussian independent measurement errors; no process noise or resampling. Finite atlas and particle errors remain."}};
    return result;
}
}
