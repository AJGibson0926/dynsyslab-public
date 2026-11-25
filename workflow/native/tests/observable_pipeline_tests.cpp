#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/workflow/Solve.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    const auto root=std::filesystem::temp_directory_path()/("dynsyslab-observable-pipeline-"+
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    core::SystemConfiguration config;config.identifier="rpe";config.parameters["kappa"]=4./3.;
    auto system=core::createSystem(config);
    const double R0=system->configuration().parameters.at("equilibrium_radius_m");
    const double T=system->coordinates().secondsPerModelTime;
    std::vector<double> time,values;
    for(std::size_t i=0;i<65;++i){double t=.001*i;time.push_back(t*T);values.push_back((1.2+.05*std::sin(6*t))*R0);values.push_back(.3*std::cos(6*t)*R0/T);}
    core::NpzRecord original;
    original.arrays["time"]=core::NpyArray::fromDoubles({65},time);
    original.arrays["values"]=core::NpyArray::fromDoubles({65,2},values);
    auto metadata=core::json::parse(R"({"system":"rayleigh-plesset","integration_form":"dimensionless","time_series_format":"dynsyslab.TimeSeries/v1","system_parameters":{"equilibrium_radius_m":0.00001,"ambient_pressure_pa":100000,"preset":{"kappa":1.3333333333333333}},"nondimensionalization":{}})");
    metadata.at("nondimensionalization").at("time_scale_s")=T;
    original.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(metadata));
    const auto source=root/"legacy.npz";core::writeNpzAtomic(source,original);
    auto job=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"echo-pipeline-test","methods":[{"name":"hdmd","settings":{"observable":"far_field_acoustic_echo","matrix_size":4,"sampling_window":30,"rank":2}},{"name":"sindy","settings":{"coordinate_observable":"far_field_acoustic_echo","polynomial_degree":3,"include_constant":false,"threshold":0,"max_iterations":1,"normalize":false,"gamma_method":"system_rhs","reconstruction_observables":["far_field_acoustic_echo"]}}]})");
    job.at("sources")=core::json::Value::Array{source.generic_u8string()};
    const auto result=workflow::executePostprocess(job,root,{});
    // The echo named in hdmd is a coordinate change that runs first. H-DMD and SINDy each publish
    // eigenfunctions and the reconstructions written over them; neither copies the training series.
    check(result.archives.size()==5&&result.methods.front()=="coordinate_transform","pipeline archive count");
    const auto fittedPath=result.archives.at(3);
    const auto fitted=core::readNpz(fittedPath);
    check(core::dataObjectKind(fitted.metadataJson())=="eigenfunctions","SINDy publishes eigenfunctions");
    check(!fitted.contains("training_states")&&!fitted.contains("training_coordinates"),"SINDy copies no training data");
    const auto fittedSet=core::readEigenfunctionSet(fitted);
    check(fittedSet.representation().domainDimension(fittedSet.representation().domains().front())==1,"echo must be the only coordinate");
    const auto coordinates=core::readNpz(result.archives.front());  // the recorded coordinate change, replayed
    check(std::abs(coordinates.at("source_states").doubles()[0]-1.2)<1e-12,"legacy radius conversion");
    check(std::abs(coordinates.at("source_states").doubles()[1]-.3)<1e-12,"legacy velocity conversion");
    check(std::abs(coordinates.at("time").doubles().back()-.064)<1e-12,"legacy time conversion");
    check(core::readNpz(source).at("values").doubles()==values,"source was modified");
    const auto model=core::loadKoopmanModel(fittedPath);
    check(model.coordinateDimension()==1&&model.stateDimension()==2,"model coordinate distinction");
    check(model.coordinateObservable().hasAnalyticJacobian(),"saved analytic echo gradient");
    auto query=core::json::parse(R"({"schema":"dynsyslab.observable-query/v1","system":{"name":"rpe"},"state":[1.2,0.3],"observable":"far_field_acoustic_echo","flow_order":2})");
    check(workflow::evaluateObservableQuery(query).at("gradient_backend").string()=="analytic","repeatable flow derivative query");
    // Reuse the saved spectrum: no second H-DMD stage is allowed in this jet fit.
    auto jetJob=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"echo-jet-test","methods":[{"name":"sindy","settings":{"coordinate_observable":"echo_jet2","polynomial_degree":4,"include_constant":false,"threshold":0,"max_iterations":1,"normalize":false,"gamma_method":"system_rhs","archive_coordinate_derivatives":true,"reconstruction_observables":["far_field_acoustic_echo","far_field_acoustic_echo_dot"]}}]})");
    jetJob.at("sources")=core::json::Value::Array{source.generic_u8string()};
    jetJob.at("methods").array()[0].at("settings").at("spectrum_archive")=result.archives.at(1).generic_u8string();
    const auto jetResult=workflow::executePostprocess(jetJob,root,{});
    check(jetResult.archives.size()==2,"jet job recomputed a spectrum");
    check(jetResult.archives.front()!=fittedPath,"new fit reused a completed archive path");
    {
        const auto again=core::readEigenfunctionSet(core::readNpz(fittedPath));
        check(again.representation().domainDimension(again.representation().domains().front())==1,"new fit overwrote prior coordinate model");
    }
    const auto jetArchive=core::readNpz(jetResult.archives.front());
    check(jetArchive.at("expansion_coefficients").shape[1]==14,"jet polynomial library size");
    check(jetArchive.at("continuous_eigenvalues").complexDoubles()==core::readNpz(result.archives.at(1)).at("continuous_eigenvalues").complexDoubles(),"jet fit changed supplied spectrum");
    check(jetArchive.at("training_coordinate_jacobians").shape==std::vector<std::size_t>({65,2,2}),"jet gradient evidence shape");
    const auto jetVelocities=jetArchive.at("training_coordinate_velocities").doubles();
    const auto catalog=system->observableCatalog();
    const auto jet=catalog.coordinatePreset("echo_jet2");
    for(std::size_t i=0;i<65;++i){  // the training states, from the source's dimensionless form
        const double t=.001*double(i),x[2]={1.2+.05*std::sin(6*t),.3*std::cos(6*t)};
        const auto q=jet.evaluateReal(x);
        check(std::abs(jetVelocities[2*i]-q[1])<1e-9,"echo derivative differs from second coordinate");
        check(std::abs(jetVelocities[2*i+1]-catalog.at("echo_ddot").evaluateReal(x)[0])<1e-8,"echo second derivative was not used in Gamma");
    }
    const auto jetModel=core::loadKoopmanModel(jetResult.archives.front());
    check(jetModel.coordinateDimension()==2&&jetModel.coordinateObservable().hasAnalyticJacobian(),"reloaded jet model lost analytic coordinates");
    auto extraction=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"postprocess","label":"echo-jet-extracted","methods":[{"name":"krom_extract","settings":{"selected_modes":[0]}}]})");
    extraction.at("sources")=core::json::Value::Array{jetResult.archives.front().generic_u8string()};
    const auto extracted=core::readNpz(workflow::executePostprocess(extraction,root,{}).archives.front());
    check(extracted.at("training_coordinate_velocities").doubles()==jetVelocities,"extraction lost derivative evidence");
    check(extracted.at("training_coordinate_jacobians").doubles()==jetArchive.at("training_coordinate_jacobians").doubles(),"extraction changed coordinate gradients");
    query.at("observable")="echo_jet2";query.at("flow_order")=0.0;
    check(workflow::evaluateObservableQuery(query).at("output_dimension").number()==2,"derivative coordinate preset query");
    query.at("observable")="echo";query.at("flow_order")=2.0;
    const auto autonomousQuery=workflow::evaluateObservableQuery(query);
    query.at("input")=.3;query.at("input_time_derivative")=-.4;
    query.at("include_input_terms")=false;
    check(workflow::evaluateObservableQuery(query).at("values").array()[0].number()==autonomousQuery.at("values").array()[0].number(),"disabled input correction changed the derivative");
    query.at("include_input_terms")=true;
    check(workflow::evaluateObservableQuery(query).at("input_terms_included").boolean(),"enabled input correction was ignored");
    auto driven=original;driven.arrays["control_values"]=core::NpyArray::fromDoubles({65,1},std::vector<double>(65,.5));
    metadata.at("control_history")=core::json::Value::Object{{"coordinate_system","dimensionless"},{"units","p_control/p0"}};
    driven.arrays["metadata_json"]=core::NpyArray::fromString(core::json::dump(metadata));
    const auto drivenPath=root/"controlled-legacy.npz";core::writeNpzAtomic(drivenPath,driven);
    auto forcedJob=jetJob;forcedJob.at("sources")=core::json::Value::Array{drivenPath.generic_u8string()};
    // Eigenfunction SINDy is autonomous: a forced (driven) record is refused whatever the
    // input-term setting (docs/architecture/TRAINING_DATA_ROLES.md). B_psi comes from input_map.
    for(bool included:{false,true}){
        forcedJob.at("methods").array()[0].at("settings").at("include_input_terms")=included;
        bool refused=false;
        try{(void)workflow::executePostprocess(forcedJob,root,{});}catch(const std::exception&){refused=true;}
        check(refused,"eigenfunction SINDy accepted a forced training record");
    }
    // A small stable fixture verifies target plumbing independently of SINDy fit quality.
    core::KoopmanSINDyResult fixture;fixture.stateDimension=1;fixture.dynamicalStateDimension=2;
    fixture.featureCount=1;fixture.center={0};fixture.scale={1};fixture.exponents={1};
    fixture.coefficients={1.};fixture.continuousEigenvalues={-1.};
    fixture.coordinateObservableNames={"far_field_acoustic_echo"};
    const auto echo=system->observableCatalog().at("far_field_acoustic_echo");
    fixture.coordinateObservableDefinition=core::json::dump(echo.expressionDefinition());
    const auto modelPath=root/"echo-control.npz";
    core::savePolynomialKoopmanModel(modelPath,core::KoopmanModel(fixture,config));
    auto solve=core::json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"echo-target-test","system":{"name":"rpe","parameters":{"kappa":1.3333333333333333}},"time":{"start":0,"stop":0.002,"output_step":0.001},"initial_conditions":{"mode":"explicit","values":[[1,0]],"count":1,"clusters":1},"controller":{"kind":"klqr","settings":{"Q":1,"R":1},"target":{"space":"observable","kind":"from_state","state_target":{"kind":"rpe_canonical","name":"rpe_multisine_equilibrium_moderate_identification"}}},"solver":{"method":"RK4","max_step":0.0001},"archive":{"enabled":true}})");
    solve.at("controller").at("settings").at("krom")=modelPath.generic_u8string();
    const auto controlled=workflow::executeSolve(solve,root,{});
    const auto output=core::readNpz(controlled.archives.front());
    check(output.at("control_observable_target_values").shape==std::vector<std::size_t>({3,1}),"scalar target archive shape");
    const auto references=output.at("control_target_values").doubles();
    const auto targets=output.at("control_observable_target_values").doubles();
    for(std::size_t i=0;i<3;++i)
        check(std::abs(targets[i]-echo.evaluateReal(references.data()+2*i)[0])<1e-13,"echo target mapping differs");
    check(core::json::parse(output.metadataJson()).at("control_target_space").string()=="observable","target objective provenance");
    solve.at("controller").at("target")=core::json::parse(R"({"space":"observable","kind":"constant","state":[0]})");
    solve.at("label")="echo-direct-target-test";
    const auto direct=workflow::executeSolve(solve,root,{});
    const auto directOutput=core::readNpz(direct.archives.front());
    check(!directOutput.contains("control_target_values"),"direct echo target manufactured a state reference");
    check(directOutput.at("control_observable_target_values").shape[1]==1,"direct target was expanded into state space");
    fixture.stateDimension=2;fixture.featureCount=2;fixture.center={0,0};fixture.scale={1,1};
    fixture.exponents={1,0,0,1};fixture.coefficients={1.,.001};
    const auto targetJet=catalog.coordinatePreset("echo_jet2");
    fixture.coordinateObservableNames=targetJet.descriptor().outputNames;
    fixture.coordinateObservableDefinition=core::json::dump(targetJet.expressionDefinition());
    const auto jetControlPath=root/"echo-jet-control.npz";
    core::savePolynomialKoopmanModel(jetControlPath,core::KoopmanModel(fixture,config));
    solve.at("controller").at("settings").at("krom")=jetControlPath.generic_u8string();
    solve.at("controller").at("target")=core::json::parse(R"({"space":"observable","kind":"from_state","state_target":{"kind":"rpe_canonical","name":"rpe_multisine_equilibrium_moderate_identification"}})");
    solve.at("label")="echo-jet-target-test";
    const auto jetControl=core::readNpz(workflow::executeSolve(solve,root,{}).archives.front());
    check(jetControl.at("control_observable_target_values").shape==std::vector<std::size_t>({3,2}),"jet target archive shape");
    const auto jetTargets=jetControl.at("control_observable_target_values").doubles();
    for(std::size_t i=0;i<3;++i){
        check(std::abs(jetTargets[2*i]-targets[i])<1e-13,"refined fit changed the echo target");
        check(std::abs(jetTargets[2*i+1]-catalog.at("echo_dot").evaluateReal(references.data()+2*i)[0])<1e-13,"jet target omitted the derivative coordinate");
    }
    std::filesystem::remove_all(root);
    std::cout<<"legacy import, echo pipeline, archive reload and derivative query passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
