#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/PeriodicTransport.hpp>
#include <dynsyslab/core/ShootingIdentification.hpp>
#include <dynsyslab/core/SequentialShooting.hpp>
#include "PeriodicFeedback.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J& v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite identification number required");return v.number();}
double option(const J& q,const char* k,double d){return q.at(k).isNull()?d:num(q.at(k));}
std::size_t count(const J& q,const char* k,std::size_t d,std::size_t max=1000000){const double x=option(q,k,double(d));if(x<0||x>max||x!=std::floor(x))throw std::invalid_argument("Identification count outside limits");return std::size_t(x);}
std::vector<double> reals(const J& v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>10000000)throw std::invalid_argument("Bounded numeric array required");std::vector<double> x;for(const auto& a:v.array())x.push_back(num(a));return x;}
J values(const std::vector<double>& x){A a;for(double v:x){if(!std::isfinite(v))throw std::runtime_error("Nonfinite identification output");a.emplace_back(v);}return a;}
J fitJson(const numerics::NonlinearLeastSquaresResult& r){return O{{"parameters",values(r.parameters)},{"squared_residual",r.squaredResidual},{"converged",r.converged},{"reason",r.reason},{"iterations",double(r.iterations)},{"evaluations",double(r.evaluations)}};}
core::KnownScalarInput knownInput(const J& q){
    if(!q.at("state_gain").isNull()||!q.at("state_reference").isNull())throw std::invalid_argument("Identification needs known applied input, not a state-feedback law");
    if(!q.at("time").isNull()||!q.at("values").isNull()){
        if(!q.at("bias").isNull()||!q.at("sines").isNull()||!q.at("limit").isNull())throw std::invalid_argument("Recorded input must contain actual applied values without analytic modifiers");
        const auto timing=q.at("measurement_timing").string("before_command");
        if(timing!="before_command"&&timing!="after_command")throw std::invalid_argument("measurement_timing must be before_command or after_command");
        return core::KnownScalarInput::recorded(reals(q.at("time")),reals(q.at("values")),
            num(q.at("before_first")),timing=="before_command");
    }
    const double bias=option(q,"bias",0),limit=option(q,"limit",0);
    if(limit<0)throw std::invalid_argument("Nonnegative input limit required");
    std::vector<std::array<double,3>> sines;
    if(!q.at("sines").isNull()){
        if(!std::holds_alternative<A>(q.at("sines").storage())||q.at("sines").array().size()>64)throw std::invalid_argument("At most 64 input sines allowed");
        for(const auto& row:q.at("sines").array()){const auto v=reals(row);if(v.size()!=3)throw std::invalid_argument("Input sine needs amplitude, angular frequency, phase");sines.push_back({v[0],v[1],v[2]});}
    }
    return [=](double t){double u=bias;for(const auto& s:sines)u+=s[0]*std::sin(s[1]*t+s[2]);return limit>0?std::clamp(u,-limit,limit):u;};
}
}
J evaluatePeriodicIdentificationQuery(const J& q,const math::PeriodicAtlas& atlas){
    core::SystemConfiguration config;const auto& spec=q.at("system");
    config.identifier=spec.at("name").string(spec.at("identifier").string());config.preset=spec.at("preset").string();
    if(!spec.at("driver").isNull())throw std::invalid_argument("Supply the known input explicitly, not a system driver");
    for(const auto& [k,v]:spec.at("parameters").object())config.parameters[k]=num(v);
    auto physical=core::createSystem(config);
    if(physical->descriptor().stateDimension!=2)throw std::invalid_argument("Periodic identification requires a planar system");
    const auto method=q.at("model").string("atlas");
    const bool feedback=q.at("operation").string()=="controlled_feedback";
    if(method!="atlas"&&method!="physical")throw std::invalid_argument("model must be atlas or physical");
    core::PeriodicTransportOptions settings;settings.phaseSamples=count(q,"table_samples",8192);
    settings.inputMaximumHarmonic=int(count(q,"input_maximum_harmonic",2047));
    settings.minimumNormalizedDeterminant=option(q,"minimum_normalized_determinant",1e-10);
    settings.frozenFrequency=option(q,"frozen_frequency",0);
    if(method=="physical"&&settings.frozenFrequency!=0)throw std::invalid_argument("Physical baseline cannot freeze its frequency");
    const auto nc=atlas.knots().front().channels.size();
    if(!q.at("state_channels").isNull()){
        const auto c=reals(q.at("state_channels"));if(c.size()!=2)throw std::invalid_argument("Two state channels required");
        for(int i=0;i<2;++i){if(c[i]<0||c[i]!=std::floor(c[i])||c[i]>=nc)throw std::invalid_argument("Invalid state channel");settings.stateChannels[i]=std::size_t(c[i]);}
    }
    if(settings.stateChannels[0]==settings.stateChannels[1]||settings.stateChannels[1]>=nc||settings.stateChannels[0]>=nc)throw std::invalid_argument("Two distinct valid state channels required");
    math::PeriodicChart chart(atlas,settings.phaseSamples);
    std::unique_ptr<core::PeriodicTransportSystem> model;
    if(method=="atlas")model=std::make_unique<core::PeriodicTransportSystem>(atlas,[&](auto x){const auto b=physical->inputMatrix(x.data());if(b.size()!=2)throw std::invalid_argument("One scalar input required");return std::array<double,2>{b[0],b[1]};},settings);
    core::ODESystem& dynamics=model?static_cast<core::ODESystem&>(*model):*physical;
    std::vector<std::string> names;
    if(method=="physical"||feedback){
        if(!std::holds_alternative<A>(q.at("observables").storage()))throw std::invalid_argument("Physical baseline needs ordered observable names matching the atlas channels");
        for(const auto& n:q.at("observables").array()){if(n.string().empty())throw std::invalid_argument("Observable name required");names.push_back(n.string());}
    }
    const auto catalog=physical->observableCatalog();std::vector<core::Observable> observables;
    std::size_t outputs=0;for(const auto& n:names){observables.push_back(catalog.at(n));outputs+=observables.back().outputDimension();}
    if((method=="physical"||feedback)&&outputs!=nc)throw std::invalid_argument("Physical readout dimension must match atlas channels");
    auto readPhysical=[&](const double* x){std::vector<double> out;for(const auto& obs:observables){const auto y=obs.evaluateReal(x);out.insert(out.end(),y.begin(),y.end());}return out;};
    auto readAll=[&](const double* z){
        if(model)return model->readout(z).value;
        return readPhysical(z);
    };
    std::vector<std::size_t> channels;
    for(double c:reals(q.at("measurement_channels"))){if(c<0||c!=std::floor(c)||c>=nc)throw std::invalid_argument("Measurement channel out of range");channels.push_back(std::size_t(c));}
    auto feedthrough=q.at("measurement_feedthrough").isNull()?std::vector<double>(channels.size(),0):reals(q.at("measurement_feedthrough"));
    if(feedthrough.size()!=channels.size())throw std::invalid_argument("Feedthrough dimension mismatch");
    const double scale=option(q,"feedthrough_scale_channel",-1);
    if(scale<-1||scale!=std::floor(scale)||scale>=nc)throw std::invalid_argument("Invalid feedthrough scale channel");
    auto measure=[&](double,const double* z,double u){const auto all=readAll(z);std::vector<double> y;
        for(std::size_t c=0;c<channels.size();++c)y.push_back(all[channels[c]]+feedthrough[c]*u*(scale<0?1:all[std::size_t(scale)]));return y;};
    auto initialize=[&](const std::vector<double>& p){
        if(p.size()!=2)throw std::invalid_argument("Initial candidate needs orbit parameter and phase");
        if(model)return p;const auto x=chart.evaluate(p[0],p[1]).value;return std::vector<double>{x[settings.stateChannels[0]],x[settings.stateChannels[1]]};
    };
    const auto time=reals(q.at("time"));const auto u=knownInput(q.at("input"));
    if(time.size()<2||time.size()>1000000)throw std::invalid_argument("At least two bounded times required");
    core::IntegrationOptions solver;solver.relativeTolerance=option(q,"relative_tolerance",1e-9);solver.absoluteTolerance=option(q,"absolute_tolerance",1e-11);solver.maximumStep=option(q,"maximum_step",0);
    if(q.at("operation").string()=="controlled_replay"){
        if(!q.at("measurements").isNull())throw std::invalid_argument("Replay is a forward solve and does not consume measurements");
        solver.initialTime=time.front();solver.finalTime=time.back();solver.outputStep=num(q.at("output_step"));
        const auto solution=core::integrateKnownInput(dynamics,initialize(reals(q.at("initial_parameters"))),solver,u);
        std::vector<double> decoded,measured,applied,sensorInput;
        for(std::size_t j=0;j<solution.sampleCount;++j){const auto all=readAll(solution.state(0,j));decoded.insert(decoded.end(),all.begin(),all.end());
            const double s=u.measurement(solution.time[j]);const auto y=measure(solution.time[j],solution.state(0,j),s);measured.insert(measured.end(),y.begin(),y.end());
            applied.push_back(u(solution.time[j]));sensorInput.push_back(s);}
        return O{{"schema","dynsyslab.periodic-result/v1"},{"method","controlled-input-replay"},{"model",method},
            {"time",values(solution.time)},{"coordinates",values(solution.values)},{"predictions",values(decoded)},
            {"predicted_measurements",values(measured)},{"applied_input",values(applied)},{"sensor_input",values(sensorInput)},
            {"channels",double(nc)},{"measurement_dimension",double(channels.size())},
            {"scope","Forward replay from supplied initial parameters and known applied input; no estimation or feedback. Recorded jumps are integrated as separate held intervals."}};
    }
    numerics::NonlinearLeastSquaresOptions optimizer;optimizer.maximumIterations=count(q,"maximum_iterations",100,10000);optimizer.finiteDifferenceStep=option(q,"finite_difference_step",1e-5);
    std::vector<core::ShootingStart> starts;
    if(!std::holds_alternative<A>(q.at("initial_candidates").storage())&&!(feedback&&!q.at("coarse_initialization").isNull()))throw std::invalid_argument("Initial candidate array required");
    const double pi=std::acos(-1.);
    for(const auto& row:q.at("initial_candidates").array()){const auto p=reals(row);if(p.size()!=2)throw std::invalid_argument("Two candidate entries required");starts.push_back({p,{atlas.knots().front().parameter,p[1]-pi},{atlas.knots().back().parameter,p[1]+pi}});}
    if(feedback){
        auto plant=core::createSystem(config);
        auto plantInitial=[&](const std::vector<double>& p){if(p.size()!=2)throw std::invalid_argument("Two plant initial parameters required");const auto all=chart.evaluate(p[0],p[1]).value;return std::vector<double>{all[settings.stateChannels[0]],all[settings.stateChannels[1]]};};
        auto plantMeasure=[&](double,const double* x,double u){const auto all=readPhysical(x);std::vector<double> y;
            for(std::size_t c=0;c<channels.size();++c)y.push_back(all[channels[c]]+feedthrough[c]*u*(scale<0?1:all[std::size_t(scale)]));return y;};
        return evaluatePeriodicFeedback(q,atlas,{*plant,dynamics,plantInitial,initialize,readPhysical,readAll,plantMeasure,measure,settings.stateChannels},starts,solver,optimizer);
    }
    if(q.at("operation").string()=="controlled_sequential_identify"){
        core::SequentialShootingOptions sequential;
        sequential.predictionSubsteps=count(q,"prediction_substeps",1,10000);
        sequential.retainOriginalStarts=q.at("retain_original_starts").boolean(false);
        for(double k:reals(q.at("update_indices"))){if(k<0||k!=std::floor(k)||k>=time.size())throw std::invalid_argument("Invalid sequential update index");sequential.updateIndices.push_back(std::size_t(k));}
        if(!q.at("forecast_final_time").isNull()||!q.at("output_step").isNull())throw std::invalid_argument("Sequential prediction uses the supplied time horizon and prediction_substeps");
        const auto estimate=core::sequentialShooting(dynamics,time,reals(q.at("measurements")),reals(q.at("measurement_variances")),starts,initialize,measure,u,sequential,solver,optimizer);
        A updates;std::size_t evaluations=0;
        for(const auto& update:estimate.updates){A attempts;
            for(const auto& a:update.identification.attempts){O item{{"valid",a.valid},{"failure",a.failure}};if(a.valid)item["fit"]=fitJson(a.fit);attempts.emplace_back(item);}
            evaluations+=update.identification.trajectoryEvaluations;
            updates.emplace_back(O{{"measurement_index",double(update.measurementIndex)},{"time",update.time},{"fit",fitJson(update.identification.best)},
                {"attempts",attempts},{"trajectory_evaluations",double(update.identification.trajectoryEvaluations)}});
        }
        std::vector<double> decoded,measured;
        for(std::size_t j=0;j<estimate.time.size();++j){const auto z=estimate.states.data()+j*estimate.stateDimension;
            const auto all=readAll(z);decoded.insert(decoded.end(),all.begin(),all.end());
            const auto y=measure(estimate.time[j],z,u.measurement(estimate.time[j]));measured.insert(measured.end(),y.begin(),y.end());}
        return O{{"schema","dynsyslab.periodic-result/v1"},{"method","controlled-sequential-shooting"},{"model",method},
            {"time",values(estimate.time)},{"coordinates",values(estimate.states)},{"information_time",values(estimate.informationTime)},
            {"predictions",values(decoded)},{"predicted_measurements",values(measured)},{"channels",double(nc)},
            {"measurement_dimension",double(channels.size())},{"updates",updates},{"trajectory_evaluations",double(evaluations)},
            {"scope","Causal expanding-window initial-state re-identification with known input. Each emitted estimate uses only its information_time prefix; past estimates are never revised. Not a Kalman filter, finite-window MHE, or a real-time guarantee."}};
    }
    const auto fit=core::identifyByShooting(dynamics,time,reals(q.at("measurements")),reals(q.at("measurement_variances")),starts,initialize,measure,u,solver,optimizer);
    O result{{"schema","dynsyslab.periodic-result/v1"},{"method","controlled-window-shooting"},{"model",method},{"fit",fitJson(fit.best)},{"best_start",double(fit.bestStart)},{"trajectory_evaluations",double(fit.trajectoryEvaluations)},{"measurement_cutoff",time.back()}};
    A attempts;for(const auto& a:fit.attempts){O row{{"valid",a.valid},{"failure",a.failure}};if(a.valid)row["fit"]=fitJson(a.fit);attempts.emplace_back(row);}result["attempts"]=attempts;
    solver.initialTime=time.front();solver.finalTime=option(q,"forecast_final_time",time.back());solver.outputStep=option(q,"output_step",time[1]-time[0]);
    if(solver.finalTime<time.back())throw std::invalid_argument("Forecast must include the measurement window");
    const auto solution=core::integrateKnownInput(dynamics,initialize(fit.best.parameters),solver,u);
    std::vector<double> decoded,measured;
    for(std::size_t j=0;j<solution.sampleCount;++j){const auto all=readAll(solution.state(0,j));decoded.insert(decoded.end(),all.begin(),all.end());const auto y=measure(solution.time[j],solution.state(0,j),u.measurement(solution.time[j]));measured.insert(measured.end(),y.begin(),y.end());}
    result["time"]=values(solution.time);result["coordinates"]=values(solution.values);result["predictions"]=values(decoded);result["predicted_measurements"]=values(measured);result["channels"]=double(nc);result["measurement_dimension"]=double(channels.size());
    result["scope"]="Local single-shooting fit using only the supplied measurement window and known input. Past-window curves are reconstructions; later curves are free forecasts. Chart and physical baseline use the same initial-state parameterization. Not a recursive filter or a finite LTI Koopman controller.";
    return result;
}
}
