#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/PeriodicTransport.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J& v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite numeric control argument required");return v.number();}
double option(const J& q,const char* key,double fallback){return q.at(key).isNull()?fallback:num(q.at(key));}
std::size_t integer(const J& q,const char* key,std::size_t fallback,std::size_t max=1000000){const auto v=option(q,key,double(fallback));if(v<0||v>max||v!=std::floor(v))throw std::invalid_argument("Integer control option outside range");return std::size_t(v);}
std::vector<double> vector(const J& q){if(!std::holds_alternative<A>(q.storage()))throw std::invalid_argument("Control vector required");std::vector<double> r;for(const auto& v:q.array())r.push_back(num(v));return r;}
J values(const std::vector<double>& x){A a;for(double v:x){if(!std::isfinite(v))throw std::runtime_error("Nonfinite controlled output");a.emplace_back(v);}return a;}
struct Input {
    double bias{},limit{};std::array<double,2> gain{},reference{};
    std::vector<std::array<double,3>> sines;
    double operator()(double t,std::array<double,2> x)const{
        double u=bias+gain[0]*(x[0]-reference[0])+gain[1]*(x[1]-reference[1]);
        for(const auto& a:sines)u+=a[0]*std::sin(a[1]*t+a[2]);
        if(limit>0)u=std::clamp(u,-limit,limit);return u;
    }
};
Input input(const J& q){Input u;u.bias=option(q,"bias",0);u.limit=option(q,"limit",0);
    if(u.limit<0)throw std::invalid_argument("Input limit must be nonnegative");
    for(auto item:{std::make_pair("state_gain",&u.gain),std::make_pair("state_reference",&u.reference)})if(!q.at(item.first).isNull()){
        const auto v=vector(q.at(item.first));if(v.size()!=2)throw std::invalid_argument("Two physical state entries required");*item.second={v[0],v[1]};}
    if(!q.at("sines").isNull()){
        if(!std::holds_alternative<A>(q.at("sines").storage())||q.at("sines").array().size()>64)throw std::invalid_argument("At most 64 input sines allowed");
        for(const auto& row:q.at("sines").array()){const auto a=vector(row);if(a.size()!=3)throw std::invalid_argument("Input sine is [amplitude,angular frequency,phase]");u.sines.push_back({a[0],a[1],a[2]});}
    }return u;
}
J statistics(const core::NativeSolution& r){return O{{"accepted_steps",double(r.statistics.acceptedSteps)},{"rejected_steps",double(r.statistics.rejectedSteps)},{"rhs_calls",double(r.statistics.rightHandSideCalls)}};}
}
J evaluatePeriodicControlQuery(const J& q,const math::PeriodicAtlas& atlas){
    core::SystemConfiguration config;const auto& spec=q.at("system");
    config.identifier=spec.at("name").string(spec.at("identifier").string());config.preset=spec.at("preset").string();
    if(!spec.at("driver").isNull())throw std::invalid_argument("Specify controlled inputs in input, not a system driver");
    for(const auto& [key,value]:spec.at("parameters").object())config.parameters[key]=num(value);
    auto physical=core::createSystem(config);
    if(physical->descriptor().stateDimension!=2)throw std::invalid_argument("Periodic planar transport requires a two-dimensional physical system");
    core::PeriodicTransportOptions settings;settings.phaseSamples=integer(q,"table_samples",8192);
    settings.inputMaximumHarmonic=int(integer(q,"input_maximum_harmonic",2047));
    settings.minimumNormalizedDeterminant=option(q,"minimum_normalized_determinant",1e-10);
    settings.frozenFrequency=option(q,"frozen_frequency",0);
    if(!q.at("state_channels").isNull()){
        const auto c=vector(q.at("state_channels"));if(c.size()!=2)throw std::invalid_argument("Two state channels required");
        for(int i=0;i<2;++i){if(c[i]<0||c[i]!=std::floor(c[i])||c[i]>=atlas.knots()[0].channels.size())throw std::invalid_argument("Invalid state channel");settings.stateChannels[i]=std::size_t(c[i]);}
    }
    core::PeriodicTransportSystem model(atlas,[&](std::array<double,2> x){
        const auto b=physical->inputMatrix(x.data());if(b.size()!=2)throw std::invalid_argument("One scalar physical input required");return std::array<double,2>{b[0],b[1]};},settings);
    const auto law=input(q.at("input"));const double sampleTime=option(q,"control_sample_time",0);
    if(sampleTime<0)throw std::invalid_argument("Nonnegative control sample time required");
    auto policy=[&](core::ScalarControlLaw f){core::ControlPolicy result;result.value=std::move(f);result.evaluation=sampleTime>0?core::ControlEvaluationMode::sampled:core::ControlEvaluationMode::continuousRhs;result.sampleTime=sampleTime;return result;};
    core::IntegrationOptions solver;solver.initialTime=option(q,"initial_time",0);solver.finalTime=num(q.at("final_time"));
    solver.outputStep=num(q.at("output_step"));solver.relativeTolerance=option(q,"relative_tolerance",1e-9);
    solver.absoluteTolerance=option(q,"absolute_tolerance",1e-11);solver.maximumStep=option(q,"maximum_step",0);
    const double count=(solver.finalTime-solver.initialTime)/solver.outputStep;
    if(!std::isfinite(count)||count<=0||count>1000000||count*atlas.knots()[0].channels.size()>10000000)throw std::invalid_argument("Controlled output size outside limits");
    if(solver.relativeTolerance<=0||solver.absoluteTolerance<=0||solver.maximumStep<0)throw std::invalid_argument("Invalid integration tolerance or step");
    std::vector<double> initial{num(q.at("parameter")),num(q.at("initial_phase"))};
    const auto decoded=model.physicalState(initial.data());std::vector<double> truthInitial{decoded[0],decoded[1]};
    if(!q.at("physical_initial_state").isNull()){truthInitial=vector(q.at("physical_initial_state"));if(truthInitial.size()!=2)throw std::invalid_argument("Two physical initial values required");}
    // Diagnose the supplied chart against the physical unforced drift. The RHS
    // is used here and in the separate reference solve, never in model rollout.
    std::array<double,2> error{},norm{};
    for(std::size_t i=0;i+1<atlas.knots().size();++i){const double p=(atlas.knots()[i].parameter+atlas.knots()[i+1].parameter)/2;
        for(int j=0;j<127;++j){const double z[2]{p,2*std::acos(-1.)*(j+.37)/127};
            const auto x=model.physicalState(z);const auto read=model.readout(z);double f[2],dz[2];physical->rhsWithInput(0,x.data(),0,f);model.rhs(0,z,dz);
            for(int c=0;c<2;++c){const double residual=read.phaseDerivative[settings.stateChannels[c]]*dz[1]-f[c];error[c]+=residual*residual;norm[c]+=f[c]*f[c];}
        }
    }
    const auto predicted=core::integrateControlled(model,initial,1,solver,policy([&](std::size_t,double t,const double* z){return law(t,model.physicalState(z));}));
    std::vector<double> output;for(std::size_t j=0;j<predicted.sampleCount;++j){const auto read=model.readout(predicted.state(0,j));output.insert(output.end(),read.value.begin(),read.value.end());}
    O result{{"schema","dynsyslab.periodic-result/v1"},{"method","controlled-periodic-chart"},
        {"time",values(predicted.time)},{"coordinates",values(predicted.values)},{"predictions",values(output)},
        {"channels",double(atlas.knots()[0].channels.size())},{"controls",values(predicted.controlValues)},{"model_statistics",statistics(predicted)},
        {"initial_decode_error",A{decoded[0]-truthInitial[0],decoded[1]-truthInitial[1]}},
        {"diagnostics",O{{"minimum_absolute_normalized_determinant",model.diagnostics().minimumAbsoluteNormalizedDeterminant},
            {"maximum_input_projection_error",model.diagnostics().maximumInputProjectionError},
            {"input_projection_relative_l2",model.diagnostics().inputProjectionRelativeL2},
            {"chart_reality_projection_l1",model.diagnostics().chartRealityProjectionL1},
            {"input_reality_projection_l1",model.diagnostics().inputRealityProjectionL1},
            {"worst_input_location",A{model.diagnostics().worstInputLocation[0],model.diagnostics().worstInputLocation[1]}},
            {"normalized_determinant_range",A{model.diagnostics().minimumNormalizedDeterminant,model.diagnostics().maximumNormalizedDeterminant}},
            {"drift_relative_l2",A{std::sqrt(error[0]/std::max(norm[0],1e-300)),std::sqrt(error[1]/std::max(norm[1],1e-300))}}}},
        {"assumptions","Model-based input map precomputed from B and a real periodic chart. Only parameter/phase evolve. Sampled diagnostics are not certificates. Readouts exclude direct input dependence unless explicitly provided separately."}};
    if(q.at("include_reference").isNull()||q.at("include_reference").boolean()){
        const auto truth=core::integrateControlled(*physical,truthInitial,1,solver,policy([&](std::size_t,double t,const double* x){return law(t,{x[0],x[1]});}));
        if(truth.time!=predicted.time)throw std::runtime_error("Reference output clock changed");
        result["reference"]=values(truth.values);result["reference_controls"]=values(truth.controlValues);result["reference_statistics"]=statistics(truth);
    }
    return result;
}
}
