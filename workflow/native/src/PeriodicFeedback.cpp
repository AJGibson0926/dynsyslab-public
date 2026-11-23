#include "PeriodicFeedback.hpp"
#include "DerivativeSettings.hpp"
#include <dynsyslab/core/OnlineShooting.hpp>
#include <dynsyslab/core/PeriodicEstimation.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double number(const J& v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite feedback number required");return v.number();}
double option(const J& q,const char* k,double d){return q.at(k).isNull()?d:number(q.at(k));}
std::size_t count(const J& q,const char* k,std::size_t d,std::size_t max=1000000){const double x=option(q,k,double(d));if(x<0||x>max||x!=std::floor(x))throw std::invalid_argument("Feedback count outside limits");return std::size_t(x);}
std::vector<double> reals(const J& v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>10000000)throw std::invalid_argument("Feedback numeric array required");std::vector<double> x;for(const auto& a:v.array())x.push_back(number(a));return x;}
J values(const std::vector<double>& v){A a;for(double x:v){if(!std::isfinite(x))throw std::runtime_error("Nonfinite feedback result");a.emplace_back(x);}return a;}
void append(std::vector<double>& out,const std::vector<double>& x){out.insert(out.end(),x.begin(),x.end());}
J fitJson(const numerics::NonlinearLeastSquaresResult& r){return O{{"parameters",values(r.parameters)},{"squared_residual",r.squaredResidual},{"converged",r.converged},{"reason",r.reason},{"iterations",double(r.iterations)},{"evaluations",double(r.evaluations)}};}
}
J evaluatePeriodicFeedback(const J& q,const math::PeriodicAtlas& atlas,const PeriodicFeedbackBinding& b,
    const std::vector<core::ShootingStart>& starts,core::IntegrationOptions solver,
    const numerics::NonlinearLeastSquaresOptions& optimizer){
    if(!q.at("input").isNull()||!q.at("measurements").isNull()||!q.at("output_step").isNull())
        throw std::invalid_argument("Online feedback generates its own inputs and measurements; use prediction_substeps for output resolution");
    const auto times=reals(q.at("time"));if(times.size()!=2||times[1]<=times[0])throw std::invalid_argument("Feedback time requires two increasing endpoints");
    const auto nc=atlas.knots().front().channels.size();const auto variance=reals(q.at("measurement_variances"));
    core::OnlineShootingOptions settings;settings.sampleTime=number(q.at("sample_time"));
    settings.firstFitIndex=count(q,"first_fit_index",1);settings.fitEvery=count(q,"fit_every",1);
    settings.retainOriginalStarts=q.at("retain_original_starts").boolean(false);
    const auto substeps=count(q,"prediction_substeps",1,10000);
    if(settings.sampleTime<=0||substeps==0)throw std::invalid_argument("Positive sample clock and substeps required");
    const double ticks=(times[1]-times[0])/settings.sampleTime;
    if(!std::isfinite(ticks)||ticks<1||ticks>999999||std::abs(ticks-std::round(ticks))>1e-9||ticks*substeps*nc>10000000)
        throw std::invalid_argument("Feedback horizon must contain a bounded integer number of sample intervals");
    settings.maximumSamples=std::size_t(std::round(ticks))+1;
    if(settings.firstFitIndex<1||settings.firstFitIndex>=settings.maximumSamples||settings.fitEvery==0)
        throw std::invalid_argument("Feedback acquisition/update indices invalid");
    const auto source=q.at("feedback_source").string("estimate");
    if(source!="estimate"&&source!="truth"&&source!="zero")throw std::invalid_argument("feedback_source must be estimate, truth or zero");
    const auto& control=q.at("controller");const auto gain=reals(control.at("state_gain"));
    const auto reference=control.at("state_reference").isNull()?std::vector<double>(2,0):reals(control.at("state_reference"));
    if(gain.size()!=2||reference.size()!=2)throw std::invalid_argument("Two planar feedback gains/reference entries required");
    const double bias=option(control,"bias",0),limit=number(control.at("limit")),warmInput=option(control,"acquisition_input",0);
    const auto controlStart=count(control,"start_index",settings.firstFitIndex),controlEvery=count(control,"every",1);
    if(limit<=0||std::abs(warmInput)>limit||controlEvery==0||controlStart<settings.firstFitIndex)
        throw std::invalid_argument("Positive control limit/cadence and acquisition before feedback required");
    core::ShootingStartProvider provider;
    if(!q.at("coarse_initialization").isNull()){
        if(warmInput!=0)throw std::invalid_argument("Unforced coarse initialization requires zero acquisition input");
        const auto& c=q.at("coarse_initialization");core::PeriodicFilterOptions o;
        o.parameters=reals(c.at("parameters"));o.phaseCount=count(c,"phase_count",512);
        o.tableSamples=count(q,"table_samples",8192);o.measurementVariances=variance;o.outputChannels={0};
        for(double channel:reals(q.at("measurement_channels")))o.measurementChannels.push_back(std::size_t(channel));
        const auto seeds=count(c,"starts",4,256);if(!seeds||seeds>o.parameters.size())throw std::invalid_argument("Invalid coarse seed count");
        provider=[&,o,seeds](const auto& t,const auto& y){
            const auto coarse=core::filterPeriodicAtlas(atlas,t,y,o);std::vector<std::pair<double,std::size_t>> best;
            for(std::size_t j=0;j<o.parameters.size();++j){const auto first=coarse.finalWeights.begin()+j*o.phaseCount;
                const auto at=std::max_element(first,first+o.phaseCount);best.emplace_back(*at,std::size_t(at-coarse.finalWeights.begin()));}
            std::stable_sort(best.begin(),best.end(),[](auto x,auto y){return x.first>y.first;});
            std::vector<core::ShootingStart> out;const double pi=std::acos(-1.);
            for(std::size_t j=0;j<seeds;++j){const auto index=best[j].second;const double phase=2*pi*(index%o.phaseCount)/o.phaseCount;
                out.push_back({{o.parameters[index/o.phaseCount],phase},{atlas.knots().front().parameter,phase-pi},{atlas.knots().back().parameter,phase+pi}});}
            return out;
        };
    }
    std::unique_ptr<core::OnlineShootingObserver> observer;
    if(source=="estimate")observer=std::make_unique<core::OnlineShootingObserver>(b.model,variance,starts,
        b.initializeEstimate,b.measureEstimate,settings,solver,optimizer,provider);
    const auto noiseSpec=parseNoiseSpec(q.at("measurement_noise"));std::vector<core::StreamingNoise> noise;
    const auto scale=q.at("measurement_noise_scale").isNull()?std::vector<double>(variance.size(),1):reals(q.at("measurement_noise_scale"));
    if(scale.size()!=variance.size())throw std::invalid_argument("One noise scale per measurement required");
    for(std::size_t j=0;j<variance.size();++j){if(scale[j]<=0)throw std::invalid_argument("Positive sensor noise scale required");auto spec=noiseSpec;spec.seed+=1000003u*j;noise.emplace_back(spec,settings.sampleTime,scale[j]);}
    std::vector<double> sampleTime,measured,clean,sampleTruth,estimate,estimateMeasurement,previousInput,commands,noiseValues;
    std::vector<double> denseTime,denseEstimate,denseMeasured,denseInformation;
    A ready,information,updates;double applied=0;std::size_t evaluations=0;double estimationSeconds=0;
    auto nextTick=[&](std::size_t k){return double(static_cast<long double>(times.front())+static_cast<long double>(settings.sampleTime)*k);};
    auto policy=[&](std::size_t,double t,const double* state){
        const auto index=sampleTime.size();const auto truth=b.decodePlant(state);auto y=b.measurePlant(t,state,applied);
        if(y.size()!=variance.size())throw std::invalid_argument("Feedback measurement dimension mismatch");
        const auto noiseless=y;for(std::size_t c=0;c<y.size();++c)y[c]=noise[c].apply(y[c]);
        core::OnlineShootingSample sample;const auto before=std::chrono::steady_clock::now();
        if(observer)sample=observer->observe(t,y,applied);
        estimationSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
        std::vector<double> decoded(nc,0),predicted(variance.size(),0);
        if(sample.ready){decoded=b.decodeEstimate(sample.state.data());predicted=b.measureEstimate(t,sample.state.data(),applied);}
        double command=applied;
        if(index<controlStart)command=warmInput;
        else if(source=="zero")command=0;
        else if((index-controlStart)%controlEvery==0){
            if(source=="estimate"&&!sample.ready)throw std::runtime_error("Feedback requested before estimator acquisition");
            const auto& value=source=="truth"?truth:decoded;
            command=std::clamp(bias+gain[0]*(value[b.stateChannels[0]]-reference[0])+gain[1]*(value[b.stateChannels[1]]-reference[1]),-limit,limit);
        }
        // Store a genuine forecast for each held control interval, before its
        // future measurement is even generated. No smoothing/reconstruction.
        if(sample.ready){
            const double end=std::min(times.back(),nextTick(index+1));
            const double eps=32*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(t));
            if(end>t+eps){const auto future=observer->forecast(end,command,substeps);
                if(future.sampleCount!=substeps+1)throw std::runtime_error("Online dense forecast clock mismatch");
                for(std::size_t j=0;j<substeps;++j){const auto z=future.state(0,j);denseTime.push_back(future.time[j]);
                    append(denseEstimate,b.decodeEstimate(z));append(denseMeasured,b.measureEstimate(future.time[j],z,j==0?applied:command));denseInformation.push_back(*sample.informationTime);}
            }else{denseTime.push_back(t);append(denseEstimate,decoded);append(denseMeasured,predicted);denseInformation.push_back(*sample.informationTime);}
        }
        sampleTime.push_back(t);append(sampleTruth,truth);append(measured,y);append(clean,noiseless);
        for(std::size_t c=0;c<y.size();++c)noiseValues.push_back(y[c]-noiseless[c]);
        append(estimate,decoded);append(estimateMeasurement,predicted);previousInput.push_back(applied);commands.push_back(command);
        ready.emplace_back(sample.ready);information.emplace_back(sample.informationTime?J(*sample.informationTime):J{});
        if(sample.fit){A attempts;for(const auto& a:sample.fit->attempts){O row{{"valid",a.valid},{"failure",a.failure}};if(a.valid)row["fit"]=fitJson(a.fit);attempts.emplace_back(row);}
            evaluations+=sample.fit->trajectoryEvaluations;updates.emplace_back(O{{"index",double(index)},{"time",t},{"fit",fitJson(sample.fit->best)},{"attempts",attempts},{"trajectory_evaluations",double(sample.fit->trajectoryEvaluations)}});}
        applied=command;return command;
    };
    solver.initialTime=times.front();solver.finalTime=times.back();solver.truncateOnFailure=true;
    solver.outputStep=std::nextafter((times.back()-times.front())/(std::round(ticks)*substeps),std::numeric_limits<double>::infinity());
    const auto solution=core::integrateControlled(b.plant,b.initializePlant(reals(q.at("initial_parameters"))),1,solver,
        {policy,core::ControlEvaluationMode::sampled,settings.sampleTime});
    const auto completed=solution.completedSamples.front();std::vector<double> truth,truthTime,plantInput;
    for(std::size_t j=0;j<completed;++j){truthTime.push_back(solution.time[j]);append(truth,b.decodePlant(solution.state(0,j)));plantInput.push_back(solution.control(0,j));}
    return O{{"schema","dynsyslab.periodic-result/v1"},{"method","online-shooting-feedback"},{"feedback_source",source},
        {"model",q.at("model").string("atlas")},{"completed",completed==solution.sampleCount},{"termination_reason",solution.terminationReasons.front()},
        {"time",values(truthTime)},{"truth",values(truth)},{"plant_input",values(plantInput)},{"channels",double(nc)},
        {"measurement_dimension",double(variance.size())},{"sample_time",values(sampleTime)},{"sample_truth",values(sampleTruth)},
        {"measurements",values(measured)},{"noiseless_measurements",values(clean)},{"noise",values(noiseValues)},
        {"estimates",values(estimate)},{"measurement_estimates",values(estimateMeasurement)},{"estimate_ready",ready},
        {"information_time",information},{"previous_applied_input",values(previousInput)},{"commands",values(commands)},
        {"dense_time",values(denseTime)},{"dense_estimates",values(denseEstimate)},{"dense_measurement_estimates",values(denseMeasured)},
        {"dense_information_time",values(denseInformation)},{"updates",updates},{"trajectory_evaluations",double(evaluations)},
        {"estimation_seconds",estimationSeconds},{"scope","Simulated causal output feedback with explicit initial acquisition, known held input, noisy input-aware sensors and expanding-window nonlinear identification. No process noise, computation delay or real-time guarantee; not KLQR or a finite LTI KROM. Truth and zero feedback sources are explicitly labeled baselines."}};
}
}
