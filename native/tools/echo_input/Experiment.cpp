#include "Experiment.hpp"
#include <dynsyslab/core/LQR.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace echo_experiment {
Frozen frozen(const core::NpzRecord& a,const core::EigenfunctionFitResult& p) {
    Frozen out{core::HankelLeftBasis{},p};
    out.basis.rows=static_cast<std::size_t>(a.at("history_delay_samples").doubles()[0]);
    out.basis.rank=a.at("singular_values").elementCount();
    out.basis.leftSingularVectors=a.at("left_singular_vectors").doubles();
    out.basis.singularValues=a.at("singular_values").doubles();
    out.basis.observationScale=a.at("observation_scale").doubles()[0];
    out.sensorStep=a.at("history_sample_time").doubles()[0];
    out.stride=static_cast<std::size_t>(std::llround(out.controlStep/out.sensorStep));
    const auto dimension=p.functions.coordinates().dimension;
    if ((dimension!=2 && dimension!=4) || p.spectrum.size()!=1
        || std::abs(out.stride*out.sensorStep-out.controlStep)>1e-12)
        throw std::invalid_argument("experiment requires one complex mode in two history coordinates");
    out.multiplier=std::exp(p.spectrum.continuousRates()[0]*out.controlStep);
    return out;
}
Measurements readMeasurements(const std::filesystem::path& path) {
    const auto a=core::readNpz(path);
    if(a.contains("states") || a.contains("values"))
        throw std::invalid_argument("identification boundary rejects state-bearing archives");
    Measurements out{a.at("time").doubles(),a.at("echo").doubles(),a.at("input").doubles()};
    if(out.time.size()!=out.echo.size() || out.time.size()!=out.input.size())
        throw std::invalid_argument("measurement shapes disagree");
    return out;
}
Lifted lift(const Measurements& m,const Frozen& f) {
    for(std::size_t k=1;k<m.time.size();++k)
        if(std::abs(m.time[k]-m.time[k-1]-f.sensorStep)>1e-8)
            throw std::invalid_argument("sensor sampling changed");
    const auto projected=core::projectHistory(m.echo,f.basis);
    Lifted out;
    const auto first=((f.basis.rows-1+f.stride-1)/f.stride)*f.stride;
    for(std::size_t k=first;k<m.time.size();k+=f.stride) {
        const auto row=k-(f.basis.rows-1);
        std::array<double,2> xi{projected.coordinates[row*f.basis.rank],
            projected.coordinates[row*f.basis.rank+1]};
        out.xi.push_back(xi); out.time.push_back(m.time[k]);
        out.psi.push_back(f.polynomial.functions.evaluate(0,
            projected.coordinates.data()+row*f.basis.rank));
    }
    for(std::size_t k=0;k+1<out.psi.size();++k)
        out.pairs.push_back({out.xi[k],out.psi[k],out.psi[k+1],m.input[first+k*f.stride]});
    return out;
}
void saveJson(const std::filesystem::path& path,const Json& value) {
    if(std::filesystem::exists(path)) throw std::runtime_error("refusing to replace completed report");
    std::ofstream stream(path); stream<<core::json::dump(value)<<'\n';
    if(!stream) throw std::runtime_error("failed to write report");
}
void simulate(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const std::string& id,const Frozen& f,double duration,double initialRadius,
    EchoController controller,Json metadata) {
    // Only this acquisition adapter knows the simulator and synthetic sensor.
    auto system=core::createSystem(config);
    const auto catalog=system->observableCatalog();
    const auto sensor=catalog.at("far_field_acoustic_echo");
    core::IntegrationOptions options;
    options.finalTime=std::ceil(duration/f.controlStep)*f.controlStep;
    options.outputStep=f.sensorStep;
    options.maximumStep=f.sensorStep; options.relativeTolerance=1e-9;
    options.absoluteTolerance=1e-11;
    core::ControlPolicy policy;
    policy.evaluation=core::ControlEvaluationMode::sampled;
    policy.sampleTime=f.sensorStep;
    policy.value=[&](std::size_t,double t,const double* x) {
        return controller(t,sensor.evaluateReal(x)[0]);
    };
    const auto solution=core::integrateControlled(*system,{initialRadius,0.0},1,options,policy);
    std::vector<double> echo(solution.sampleCount);
    for(std::size_t k=0;k<echo.size();++k) echo[k]=sensor.evaluateReal(solution.state(0,k))[0];
    metadata.object()["schema"]="dynsyslab.echo-input-record/v1";
    metadata.object()["sensor"]="far_field_acoustic_echo: intrinsic zero-input convention, no direct feedthrough";
    metadata.object()["sensor_step"]=f.sensorStep;
    metadata.object()["control_step"]=f.controlStep;
    metadata.object()["truth_archive"]=id+".truth.npz";
    core::NpzRecord record;
    record.arrays["time"]=core::NpyArray::fromDoubles({echo.size()},solution.time);
    record.arrays["echo"]=core::NpyArray::fromDoubles({echo.size()},echo);
    record.arrays["input"]=core::NpyArray::fromDoubles({echo.size()},solution.controlValues);
    record.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(metadata));
    core::writeNpzAtomicNew(folder/(id+".measurements.npz"),record);
    core::NpzRecord truth;
    truth.arrays["time"]=record.arrays["time"];
    truth.arrays["states"]=core::NpyArray::fromDoubles({solution.sampleCount,2},solution.values);
    truth.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(Json::Object{
        {"schema","dynsyslab.echo-input-truth/v1"},{"role","simulation truth; never supplied to regression/controller"},
        {"measurements",id+".measurements.npz"}}));
    core::writeNpzAtomicNew(folder/(id+".truth.npz"),truth);
}
Json metrics(const Measurements& actual,const Measurements& target,double start) {
    if(actual.time.size()!=target.time.size()) throw std::invalid_argument("comparison sample mismatch");
    double e2=0,et2=0,u2=0,ut2=0,dot=0,aa=0; std::size_t count=0,saturated=0;
    for(std::size_t k=0;k<actual.time.size();++k) if(actual.time[k]>=start) {
        const double e=actual.echo[k]-target.echo[k], u=actual.input[k]-target.input[k];
        e2+=e*e; et2+=target.echo[k]*target.echo[k]; u2+=u*u;
        ut2+=target.input[k]*target.input[k]; dot+=actual.input[k]*target.input[k];
        aa+=actual.input[k]*actual.input[k]; ++count;
        saturated+=std::abs(actual.input[k])>=0.599999;
    }
    return Json::Object{{"echo_nrmse",std::sqrt(e2/std::max(et2,1e-300))},
        {"driver_nrmse",std::sqrt(u2/std::max(ut2,1e-300))},
        {"driver_cosine_similarity",dot/std::sqrt(std::max(aa*ut2,1e-300))},
        {"saturation_fraction",double(saturated)/count},{"start_time",start}};
}
Json controlTrial(const core::SystemConfiguration& config,const std::filesystem::path& folder,
    const std::string& id,const Frozen& f,const core::HistoryInputMap& map,
    const Measurements& target,const Lifted& reference,double penaltyMultiplier) {
    // Target input is intentionally NOT captured by this controller.
    core::CausalHistoryProjector projector(f.basis);
    double z2=0; for(const auto& z:reference.psi) z2+=std::norm(z);
    const double zScale=std::sqrt(z2/reference.psi.size());
    const auto lambda=f.polynomial.spectrum.continuousRates()[0];
    const auto integral=(f.multiplier-1.0)/lambda;
    const std::vector<double> generator{lambda.real(),-lambda.imag(),lambda.imag(),lambda.real()};
    const double penalty=std::max(1e-16,penaltyMultiplier*zScale*zScale/(0.1*0.1));
    const auto b0=map.evaluate(0.0)/integral;
    const auto fixedGain=core::solveContinuousLQR(generator,{b0.real(),b0.imag()},
        {1,0,0,1},{penalty},2,1).gain;
    std::size_t tick=0;
    double held=0;
    EchoController controller=[&](double t,double echo) {
        projector.append(echo);
        const auto now=tick++;
        if(now%f.stride || t<6.0-1e-8) return held;
        const auto xi=projector.projectLatest();
        if(!xi) return 0.0;
        const auto index=static_cast<std::size_t>(std::llround((t-reference.time.front())/f.controlStep));
        if(index+1>=reference.psi.size()) return held;
        const auto z=f.polynomial.functions.evaluate(0,xi->data());
        const auto b=map.evaluate(z)/integral;
        const auto gain=map.degree==0 ? fixedGain : core::solveContinuousLQR(generator,
            {b.real(),b.imag()},{1,0,0,1},{penalty},2,1).gain;
        const auto error=z-reference.psi[index];
        const double value=-(gain[0]*error.real()+gain[1]*error.imag());
        if(!std::isfinite(value)) throw std::runtime_error("nonfinite history feedback");
        held=std::clamp(value,-0.6,0.6); return held;
    };
    simulate(config,folder,id,f,target.time.back(),1.0,controller,Json::Object{
        {"role","controlled attempt"},{"input_penalty_multiplier",penaltyMultiplier},
        {"controller","existing solveContinuousLQR and u=-K(psi-psi_target); no feedforward or target input"},
        {"input_penalty",penalty}});
    const auto result=readMeasurements(folder/(id+".measurements.npz"));
    return Json::Object{{"id",id},{"late",metrics(result,target,30.0)},
        {"whole_driven",metrics(result,target,6.0)},
        {"input_penalty_multiplier",penaltyMultiplier}};
}
} // namespace echo_experiment
