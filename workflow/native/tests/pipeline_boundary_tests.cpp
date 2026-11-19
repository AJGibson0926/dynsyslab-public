#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow_c_api.h>
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=dynsyslab::core;namespace w=dynsyslab::workflow;namespace j=c::json;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F> void rejects(F f){try{f();}catch(const std::logic_error&){return;}throw std::runtime_error("Invalid workflow contract accepted");}
j::Value job(const std::vector<std::filesystem::path>& sources,const char* method,j::Value settings) {
    j::Value::Array paths;for(const auto& p:sources)paths.emplace_back(p.generic_u8string());
    return j::Value::Object{{"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","pipeline-boundaries"},{"sources",paths},{"methods",j::Value::Array{j::Value::Object{
            {"name",method},{"settings",settings}}}},{"archive",j::Value::Object{{"enabled",true}}}};
}
void checkPublicRun(const j::Value& request,const std::filesystem::path& root) {
    dsl_workflow_result_handle* handle=nullptr;
    const auto encoded=j::dump(request),destination=root.generic_u8string();
    const auto status=dsl_workflow_run_json(encoded.c_str(),destination.c_str(),nullptr,nullptr,&handle);
    if(status!=DSL_STATUS_OK) {
        if(handle)dsl_workflow_result_destroy(handle);
        throw std::runtime_error("History/measurement workflow failed through the public C ABI");
    }
    std::size_t size=0;dsl_workflow_result_json(handle,nullptr,0,&size);
    std::vector<char> response(size);
    const auto copied=dsl_workflow_result_json(handle,response.data(),size,nullptr);
    dsl_workflow_result_destroy(handle);
    check(copied==DSL_STATUS_OK&&!j::parse(response.data()).at("archives").array().empty(),
          "Public workflow did not return its archives");
}
}
int main(){std::filesystem::path root;try{
    root=std::filesystem::temp_directory_path()/("dsl-pipeline-boundary-"+std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    std::vector<std::filesystem::path> paths;
    for(std::size_t run=0;run<2;++run) {
        std::vector<double> t,e;for(std::size_t k=0;k<640;++k){t.push_back(.01*double(k));e.push_back((1+.2*run)*std::sin(5*t.back()+.4*run));}
        c::NpzRecord a;a.arrays["time"]=c::NpyArray::fromDoubles({640},t);a.arrays["echo"]=c::NpyArray::fromDoubles({640},e);
        a.arrays["metadata_json"]=c::NpyArray::fromString(R"({"forcing_type":"free","controller_kind":"none","data_role":"free_training"})");
        paths.push_back(root/("measured-"+std::to_string(run)+".npz"));c::writeNpzAtomic(paths.back(),a);
    }
    auto settings=j::parse(R"({"observation_array":"echo","delays":80,"hankel_rank":2,"coordinate_rank":2,"hankel_oversampling":2,"hankel_subspace_iterations":6,"polynomial_degree":1,"pair_stride":2,"minimum_continuous_imaginary":4.5,"maximum_continuous_imaginary":5.5})");
    const auto fitted=w::executeJob(job(paths,"history_basis",settings),root/"archives",{});
    check(fitted.archives.size()==1,"Pooled basis produced multiple archives");
    const auto transform=c::loadFrozenHistoryTransform(fitted.archives[0]);
    auto projectionSettings=settings;projectionSettings.object()["transform"]=fitted.archives[0].generic_u8string();
    const auto projected=w::executeJob(job(paths,"history_projection",projectionSettings),root/"archives",{});
    checkPublicRun(job(paths,"history_basis",settings),root/"abi");
    checkPublicRun(job(paths,"history_projection",projectionSettings),root/"abi");
    check(projected.archives.size()==2,"Projection did not preserve independent sources");
    auto wrongUnits=projectionSettings;wrongUnits.at("time_units")="seconds";
    rejects([&]{(void)w::executeJob(job(paths,"history_projection",wrongUnits),root/"archives",{});});
    auto wrongMeasurement=projectionSettings;wrongMeasurement.at("observable")="other-sensor";
    rejects([&]{(void)w::executeJob(job(paths,"history_projection",wrongMeasurement),root/"archives",{});});
    const auto projection=c::readNpz(projected.archives[0]);
    check(projection.at("values").shape==std::vector<std::size_t>{561,2},"Projection warmup/shape wrong");
    auto coordinateJob=job({projected.archives[0]},"dmd",j::Value::Object{{"rank",2.}});
    coordinateJob.at("methods").array().push_back(j::Value::Object{{"name","sindy"},
        {"settings",j::Value::Object{{"polynomial_degree",1.},{"threshold",0.},{"input_dimension",1.}}}});
    const auto coordinateFit=w::executeJob(coordinateJob,root/"archives",{});
    const auto coordinateModel=c::loadKoopmanModel(coordinateFit.archives.back());
    check(coordinateModel.requiresExternalCoordinates()&&!coordinateModel.hasPhysicalSystem(),"Coordinate SINDy invented physical state dependence");
    check(coordinateModel.inputDimension()==1,"Measured SINDy lost the declared actuator count");
    check(coordinateModel.coordinates().providerIdentity==transform.identity(),"SINDy did not retain the frozen projection identity");
    const auto composed=w::executeJob(job(paths,"echo_history_koopman",settings),root/"archives",{});
    const auto fitModel=c::loadKoopmanModel(composed.archives[0]);
    checkPublicRun(job(paths,"echo_history_koopman",settings),root/"abi");
    check(!fitModel.hasPhysicalSystem()&&!fitModel.physicalStateDimension(),"Echo-only fitting invented an ODE or state");
    check(fitModel.coordinates().providerIdentity==c::loadFrozenHistoryTransform(composed.archives[0]).identity(),"KROM lost frozen basis binding");
    check(fitModel.datasetReferences().size()==2,"KROM lost pooled training pointers");
    c::EigenfunctionCoordinates q;q.dimension=2;q.kind=c::EigenfunctionCoordinateKind::causalHistory;
    c::FrozenHistoryTransform simulatorTransform(transform.basis(),.01,{0,1},"x");
    q.providerIdentity=simulatorTransform.identity();
    c::PolynomialEigenfunctionSet bank({2,{1,0,0,1},{0,0},{1,1}},
        {1.,{0,1}},c::EigenfunctionSpectrum::fromContinuousRates({{-1,5}},.01),q);
    c::SystemConfiguration systemConfig;systemConfig.identifier="duffing";systemConfig.preset="double_well_reference";
    auto system=c::createSystem(systemConfig);
    c::KoopmanModel physical(bank,systemConfig);physical.setControlMap("constant",c::ConstantControlMap{{.1},1,1,"test"});
    c::NpzRecord extra;c::writeFrozenHistoryTransform(extra,simulatorTransform);
    const auto modelPath=root/"controller-model.npz";c::saveKoopmanModel(modelPath,physical,"{}",extra.arrays);
    j::Value controller=j::Value::Object{{"kind","klqr"},{"settings",j::Value::Object{
        {"krom",modelPath.generic_u8string()},{"control_representation","constant"},{"Q",1.0},{"R",1.0}}},
        {"target",j::Value::Object{{"space","observable"},{"kind","constant"},{"state",j::Value::Array{0.,0.}}}}};
    auto bound=w::buildController(controller,{},*system,2,.001,{.2,0});
    check(bound.policy.evaluation==c::ControlEvaluationMode::sampled&&bound.sampleTime==.01,"History control acquired wrong clock");
    const double state[]{.2,0};
    for(std::size_t k=0;k<80;++k)(void)bound.policy.value(0,.01*double(k),state);
    check(bound.policy.value(1,0,state)==0,"Second controller inherited first stream history");
    rejects([&]{(void)w::buildController(controller,j::parse(R"({"control_evaluation":"continuous_rhs"})"),*system,1,.01,{.2,0});});
    rejects([&]{(void)w::buildController(controller,j::parse(R"({"control_sample_time":0.02})"),*system,1,.01,{.2,0});});
    bound=w::buildController(controller,{},*system,2,.001,{.2,0});
    c::IntegrationOptions integration;integration.finalTime=1;integration.outputStep=.01;integration.maximumStep=.005;
    const auto solution=c::integrateControlled(*system,{.2,0,.3,0},2,integration,bound.policy);
    check(solution.sampleCount==101&&solution.trajectoryCount==2,"Standard integrator failed with history control");
    c::KoopmanModel measured(bank,1);measured.setReconstruction("echo",{{1.},1,1,0,0,0});
    const auto estimatorModel=root/"estimator.npz";c::saveKoopmanModel(estimatorModel,measured);
    const auto observer=w::executeJob(job({paths[0]},"klqe",j::Value::Object{
        {"krom",estimatorModel.generic_u8string()},{"measurement","echo"},{"initial_lifted_values",j::Value::Array{0.,0.}}}),root/"archives",{});
    const auto estimate=c::readNpz(observer.archives[0]);
    check(!estimate.contains("true_states")&&!estimate.contains("state_reconstruction_matrix"),"Measured observer invented physical outputs");
    check(estimate.at("lifted_estimate").shape==std::vector<std::size_t>{640,2},"Measured observer output shape wrong");
    const auto reducedObserver=w::executeJob(job({paths[0]},"klqe",j::Value::Object{
        {"krom",estimatorModel.generic_u8string()},{"measurement","echo"},{"balanced_truncation",1.},
        {"initial_lifted_values",j::Value::Array{0.,0.}}}),root/"archives",{});
    const auto reducedEstimate=c::readNpz(reducedObserver.archives[0]);
    const auto reducedMeta=j::parse(reducedEstimate.metadataJson());
    check(reducedMeta.at("balanced_truncation").at("order").number()==1
        &&reducedMeta.at("balanced_truncation").at("full_order").number()==2,"KLQE balanced_truncation evidence missing");
    check(reducedEstimate.at("lifted_estimate").shape==std::vector<std::size_t>{640,2},
        "Reduced-gain KLQE must still estimate the full lifted state");
    physical.setReconstruction("echo",{{1.},1,1,0,0,0});
    physical.setReconstruction("state",{{1.,{0,1}},2,1,0,0,0});
    const auto physicalObserver=root/"physical-observer.npz";c::saveKoopmanModel(physicalObserver,physical);
    auto sensor=c::readNpz(paths[0]);sensor.arrays["values"]=sensor.at("echo");
    sensor.arrays["metadata_json"]=c::readNpz(physicalObserver).at("metadata_json");
    const auto sensorPath=root/"sensor-with-physics.npz";c::writeNpzAtomic(sensorPath,sensor);
    const auto decoded=w::executeJob(job({sensorPath},"klqe",j::Value::Object{
        {"krom",physicalObserver.generic_u8string()},{"measurement","echo"},
        {"initial_lifted_values",j::Value::Array{0.,0.}}}),root/"archives",{});
    const auto decodedArchive=c::readNpz(decoded.archives[0]);
    check(!decodedArchive.contains("true_states"),"Sensor values were treated as physical truth");
    check(decodedArchive.at("values").shape==std::vector<std::size_t>{640,2},"Optional state decoder failed without state measurements");
    std::filesystem::remove_all(root);std::cout<<"Pipeline boundary workflow tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';std::cerr<<"Test evidence retained at "<<root<<'\n';return 1;}}
