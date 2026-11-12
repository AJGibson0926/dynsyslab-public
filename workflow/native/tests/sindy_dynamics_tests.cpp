#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/core/DynamicsLibrary.hpp>
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab;
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
template<class F>void rejects(F f){bool ok=false;try{f();}catch(const std::exception&){ok=true;}check(ok,"invalid input accepted");}
int main(){try{
    std::vector<double> t,q,dq;
    // Decaying nonlinear oscillator in an exactly polynomial chart:
    // x=a exp(-t), y=sin(t) exp(-t/3); fit independently varied sample points
    // for coefficient identifiability, then validate an unseen point.
    for(std::size_t k=0;k<400;++k){
        const double x=.1+std::sin(.37*k),y=std::cos(.51*k);
        t.push_back(.01*k);q.insert(q.end(),{x,y});dq.insert(dq.end(),{y,-x-.2*y+x*x});
    }
    const core::TimeSeries coordinates(t,q,2),rates(t,dq,2);
    const auto library=core::polynomialDynamicsLibrary(2,4);
    double probe[]{.31,-.12};
    check(library.outputDimension()==15,"degree-four library must include 15 terms");
    core::DynamicsLibraryOptions broadOptions;broadOptions.polynomialDegree=2;
    const auto broad=core::buildDynamicsLibrary(coordinates,broadOptions);
    check(broad.observable.outputDimension()>60&&!broad.excluded.empty(),"mixed library/domain screening");
    for(double v:broad.observable.evaluateReal(probe))check(std::isfinite(v),"nonfinite broad candidate");
    core::SindyDynamicsOptions options;options.threshold=.001;
    const auto model=core::fitSindyDynamics(coordinates,rates,library,options);
    const auto f=model.evaluate(probe);
    check(std::abs(f[0]+.12)<1e-10&&std::abs(f[1]-(-.31+.024+.31*.31))<1e-10,"polynomial recovery");
    const auto restored=core::SindyDynamicsModel::fromDefinition(core::json::parse(core::json::dump(model.definition())));
    check(std::abs(restored.evaluate(probe)[1]-f[1])<1e-12,"portable model round trip");
    core::IntegrationOptions integration;integration.finalTime=.1;integration.outputStep=.01;
    check(core::integrateSindyDynamics(restored,{.31,-.12},integration).sampleCount==11,"standard solver adapter");
    options.threshold=1e9;
    const auto zero=core::fitSindyDynamics(coordinates,rates,library,options);
    check(zero.evaluate(probe)==std::vector<double>({0,0}),"empty sparse support");
    auto x=core::AnalyticExpression::state(2,0),y=core::AnalyticExpression::state(2,1);
    const auto map=core::Observable::fromExpressions({"map","map","","",{},false},{x*x,y+x});
    const auto transformed=core::transformCoordinates(coordinates,map,&rates);
    check(std::abs(transformed.derivatives->value(10,0)-2*q[20]*dq[20])<1e-12,"chain rule pushforward");
    std::vector<double> u,tt;
    for(std::size_t k=0;k<101;++k){tt.push_back(.01*k);u.push_back(std::pow(tt.back(),4));}
    const auto derivative=core::differentiateTimeSeries(core::TimeSeries(tt,u,1));
    check(derivative.sampleCount()==97&&std::abs(derivative.value(50,0)-4*std::pow(derivative.time(50),3))<1e-12,"fourth order finite difference");
    core::PooledHankelBasis basis;basis.delays=3;basis.rank=2;
    basis.leftSingularVectors={1,0,0,0,0,1};basis.singularValues={2,4};
    core::FrozenHistoryTransform history(basis,.01,{0,1},"test");
    const core::TimeSeries series(tt,u,1);const auto latent=core::transformHistoryCoordinates(series,history,&series);
    check(latent.coordinates.componentCount()==2&&latent.coordinates.sampleCount()==99,"history shape");
    check(latent.coordinates.time(0)==tt[2]&&std::abs(latent.coordinates.value(9,1)-u[11]/4)<1e-14,"causal end timestamp and frozen scaling");
    rejects([&]{core::transformHistoryCoordinates(series,core::FrozenHistoryTransform(basis,.02));});
    // A signed ramp on an irregular clock has an exact trapezoidal integral,
    // including when the trailing window starts partway through a segment.
    const core::TimeSeries ramp({3,3.25,4,5,7},{-3,-2.5,-1,1,5},1);
    for(double width : {.125,.75,1.,4.}) {
        const auto window=core::windowIntegralCoordinates(ramp,width);
        for(std::size_t k=0;k<window.sampleCount();++k) {
            const double time=window.time(k),value=2*time-9;
            check(time>=3+width&&window.value(k,0)==value,"integral alignment");
            check(std::abs(window.value(k,1)-width*(value-width))<1e-13,"irregular signed ramp integral");
        }
    }
    rejects([&]{core::windowIntegralCoordinates(ramp,0);});
    rejects([&]{core::windowIntegralCoordinates(ramp,5);});
    rejects([&]{core::windowIntegralCoordinates(coordinates,1);});
    rejects([&]{core::fitSindyDynamics(coordinates,core::TimeSeries(t,std::vector<double>(400),1),library);});
    // Public workflow integration with an exact oscillator and model-free rates.
    std::vector<double> orbit,clock;for(std::size_t k=0;k<=2000;++k) {
        const double time=.01*k;clock.push_back(time);orbit.push_back(std::exp(-.05*time)*std::cos(time));
        orbit.push_back(std::exp(-.05*time)*std::sin(time));
    }
    const auto root=std::filesystem::temp_directory_path()/("dynsyslab-sindy-dynamics-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    core::NpzRecord archive;archive.arrays["time"]=core::NpyArray::fromDoubles({clock.size()},clock);
    archive.arrays["values"]=core::NpyArray::fromDoubles({clock.size(),2},orbit);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(R"({"forcing_type":"free","coordinate_provider":"measured"})");
    const auto source=root/"orbit.npz";core::writeNpzAtomicNew(source,archive);
    core::NpzRecord measured;
    measured.arrays["time"]=core::NpyArray::fromDoubles({3},{6,7,9});
    measured.arrays["echo"]=core::NpyArray::fromDoubles({3},{-1,-1,3});
    measured.arrays["metadata_json"]=core::NpyArray::fromString("{}");
    const auto measuredFile=root/"echo.npz";core::writeNpzAtomicNew(measuredFile,measured);
    auto windowJob=core::json::parse(R"({"kind":"postprocess","methods":[{"name":"coordinate_transform","settings":{"kind":"window_integral","measurement_array":"echo","window_duration":2}}]})");
    windowJob.at("sources")=core::json::Value::Array{measuredFile.generic_u8string()};
    const auto windowArchive=core::readNpz(workflow::executePostprocess(windowJob,root,{}).archives.front());
    check(windowArchive.at("time").doubles()==std::vector<double>{9}
        &&windowArchive.at("values").doubles()==std::vector<double>{3,2}
        &&!windowArchive.contains("coordinate_derivatives"),"one complete measured window");
    check(core::json::parse(windowArchive.metadataJson()).at("history_span").number()==2,"integral span provenance");
    auto job=core::json::parse(R"({"kind":"postprocess","methods":[{"name":"coordinate_transform","settings":{"kind":"identity"}},{"name":"sindy_dynamics","settings":{"use_previous_output":true,"polynomial_degree":1,"threshold":0.001,"rollout_starts":[1],"rollout_duration":1}}]})");
    job.at("sources")=core::json::Value::Array{source.generic_u8string()};
    const auto outputs=workflow::executePostprocess(job,root,{});
    const auto fitted=core::readNpz(outputs.archives.back());const auto meta=core::json::parse(fitted.metadataJson());
    check(meta.at("training").at("against_numerical").at("relative_rms").array()[0].number()<1e-7,"workflow numerical fit");
    check(meta.at("training").at("rollouts").array()[0].at("status").string()=="complete","workflow rollout");
    auto validateJob=core::json::parse(R"({"kind":"postprocess","methods":[{"name":"dynamics_validation","settings":{}}]})");
    validateJob.at("sources")=core::json::Value::Array{outputs.archives.front().generic_u8string()};
    auto& validationSettings=validateJob.at("methods").array()[0].at("settings");
    validationSettings.at("model_archive")=outputs.archives.back().generic_u8string();
    validationSettings.at("compare_model_archive")=outputs.archives.back().generic_u8string();
    const auto validation=core::readNpz(workflow::executePostprocess(validateJob,root,{}).archives.front());
    check(core::json::parse(validation.metadataJson()).at("model_comparison").at("rmse").array()[0].number()==0,"model comparison must use the same native evaluation");
    validateJob.at("schema")="dynsyslab.compute-job/v1";
    const auto published=workflow::executeJob(validateJob,std::filesystem::relative(root,std::filesystem::current_path()));
    bool indexed=false;
    for(const auto& entry:core::loadArchiveIndex(root))
        indexed=indexed||std::filesystem::weakly_canonical(entry.path)==std::filesystem::weakly_canonical(published.archives.front());
    check(indexed,"relative archive roots must publish Viewer index entries");
    std::vector<double> bubble;for(double time:clock)bubble.insert(bubble.end(),{1.2+.01*std::cos(time),-.01*std::sin(time)});
    archive.arrays["values"]=core::NpyArray::fromDoubles({clock.size(),2},bubble);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(R"({"system":"rpe","forcing_type":"free","integration_form":"dimensionless"})");
    const auto bubblePath=root/"bubble.npz";core::writeNpzAtomicNew(bubblePath,archive);
    auto jetJob=core::json::parse(R"({"kind":"postprocess","methods":[{"name":"coordinate_transform","settings":{"kind":"observable","observable":"echo_jet2","system_derivatives":true}}]})");
    jetJob.at("sources")=core::json::Value::Array{bubblePath.generic_u8string()};
    const auto jet=core::readNpz(workflow::executePostprocess(jetJob,root,{}).archives.front());
    const auto jetValues=jet.at("values").doubles(),jetRates=jet.at("coordinate_derivatives").doubles();
    for(std::size_t k=0;k<clock.size();++k)check(std::abs(jetValues[2*k+1]-jetRates[2*k])<1e-9,"RPE preset chain identity");
    auto& receiverSettings=jetJob.at("methods").array()[0].at("settings");
    receiverSettings.at("observable")="monopole_pressure";
    receiverSettings.at("receiver")=core::json::parse(R"({"position_m":[0,0.006,0],"sound_speed_m_s":1480})");
    const auto receiver=core::readNpz(workflow::executePostprocess(jetJob,root,{}).archives.front());
    const auto receiverMeta=core::json::parse(receiver.metadataJson()).at("acoustic_receiver");
    check(receiver.at("time").doubles()==clock && receiverMeta.at("clock").string()=="emission-time"
        && std::abs(receiverMeta.at("travel_time_s").number()-.006/1480)<1e-16,
        "receiver timestamps and propagation delay");
    const auto pressures=receiver.at("values").doubles();
    for(std::size_t k=0;k<clock.size();++k)
        check(std::abs(.006*pressures[k]-jetValues[2*k])<1e-12,"leading monopole alone at receiver");
    receiverSettings.at("receiver").at("position_m")=core::json::Value::Array{1.e-6,0.,0.};
    rejects([&]{workflow::executePostprocess(jetJob,root,{});});
    std::cout<<"SINDy dynamics, coordinate transforms and pipeline checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
