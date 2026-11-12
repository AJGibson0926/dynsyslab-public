#include <dynsyslab/workflow/KoopmanFeedback.hpp>
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
namespace c=dynsyslab::core;namespace w=dynsyslab::workflow;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f) {
    try{f();}catch(const std::logic_error&){return;}
    throw std::runtime_error("Unsafe feedback configuration accepted");
}
c::KoopmanModel fixture(const c::SystemConfiguration& system) {
    // Linear Duffing: x''=2x-x'+u. Exact eigenfunctions 2x+v, -x+v.
    auto model=c::KoopmanModel::fromPolynomialCoefficients(system,2,
        {1,0,0,1},{2.,1.,-1.,1.},{1.,-2.});
    model.setReconstruction("sensor",{{1./3,-1./3},1,2,0,0,0});
    model.setReconstruction("state",{{1./3,-1./3,1./3,2./3},2,2,0,0,0});
    model.setControlMap("fit",c::ConstantControlMap{{1.,1.},2,1,"exact linear test fixture"});
    return model;
}
void measurementOnly() {
    c::EigenfunctionCoordinates q;q.dimension=2;q.kind=c::EigenfunctionCoordinateKind::causalHistory;
    c::PolynomialEigenfunctionSet bank({2,{1,0,0,1},{0,0},{1,1}},{1.,0.,0.,1.},
        c::EigenfunctionSpectrum::fromContinuousRates({1.,-2.}),q);
    c::KoopmanModel model(bank,1);
    model.setReconstruction("sensor",{{1./3,-1./3},1,2,0,0,0});
    model.setControlMap("fit",c::ConstantControlMap{{1.,1.},2,1,"fixture"});
    const auto design=c::designKoopmanLQEFromMeasurements(model,"sensor",{1,0,0,1},{.01});
    w::KoopmanFeedback feedback(model,design,"fit",{0,0},{1,0,0,1},{.1});
    // No ODE adapter or latent-coordinate values exist at this boundary.
    check(!feedback.model().hasPhysicalSystem(),"Measurement runtime needs a physical model");
    auto first=feedback.sample(0,{.4},{0});
    check(first.commandedInput[0]==0,"Controller read hidden plant truth at initialization");
    auto second=feedback.sample(.01,{.3},{0});
    check(second.estimate.lifted[0]!=0&&second.commandedInput[0]!=0,"Measured feedback did not respond");
    auto wrong=model.extract({1,0});
    rejects([&]{w::KoopmanFeedback bad(wrong,design,"fit",{0,0},{1,0,0,1},{.1});});
}
c::json::Value job(const std::filesystem::path& modelPath) {
    auto value=c::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"estimated-klqr-test",
      "system":{"name":"duffing","parameters":{"alpha":2,"beta":0,"delta":1}},
      "time":{"start":0,"stop":12,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","values":[[0.1,0],[-0.1,0]],"count":2,"clusters":1},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_sample_time":0.005},
      "controller":{"kind":"klqr","settings":{
        "control_representation":"fit","Q":1,"R":0.1,"control_limit":0.8,
        "active_time_window":[0.1,12],
        "estimator":{"measurement":"sensor","measurement_observable":"position",
          "output":"state","process_covariance":1,"measurement_covariance":0.001,
          "initial_lifted":[0,0]}}},
      "archive":{"enabled":true}})");
    value.at("controller").at("settings").at("krom")=modelPath.generic_u8string();
    return value;
}
}
int main(int argc,char** argv) {
    const auto root=argc>1?std::filesystem::path(argv[1]):std::filesystem::temp_directory_path()/
        ("dynsyslab-feedback-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        measurementOnly();
        std::filesystem::create_directories(root);
        c::SystemConfiguration config;config.identifier="duffing";
        config.parameters={{"alpha",2},{"beta",0},{"delta",1}};
        const auto model=fixture(config);
        const auto modelPath=root/"feedback-fixture.npz";
        c::savePolynomialKoopmanModel(modelPath,model);
        auto request=job(modelPath);
        if(argc>1) {std::ofstream stream(root/"feedback-job.json");stream<<c::json::dump(request);}
        (void)w::planJob(request);
        const auto result=w::executeJob(request,root);
        check(result.archives.size()==2,"Independent observer streams not archived");
        for(const auto& path:result.archives) {
            const auto archive=c::readNpz(path);
            const auto truth=archive.at("values").doubles();
            const auto estimate=archive.at("state_estimate").doubles();
            const auto applied=archive.at("applied_input").doubles();
            const auto command=archive.at("commanded_input").doubles();
            const auto previous=archive.at("previous_applied_input").doubles();
            const auto ticks=archive.at("estimator_time").doubles();
            check(ticks.size()==2401&&truth.size()==1202,"Observer clock coupled to output/adaptive steps");
            check(command[0]==0&&previous[0]==0,"Observer initialized with truth");
            for(std::size_t k=1;k<ticks.size();++k) {
                check(previous[k]==applied[k-1],"Estimator did not receive applied input");
                check(std::abs(applied[k])<=.8,"Actuator bound violated");
                if(ticks[k]<.1)check(applied[k]==0,"Active control window ignored");
            }
            const auto size=truth.size(),e=estimate.size();
            const double finalState=std::hypot(truth[size-2],truth[size-1]);
            const double finalError=std::hypot(truth[size-2]-estimate[e-2],truth[size-1]-estimate[e-1]);
            std::cout<<"Final state norm="<<finalState<<", estimation error="<<finalError<<'\n';
            check(finalState<1e-4&&finalError<1e-5,"Measured KLQR did not stabilize the exact unstable fixture");
        }
        auto system=c::createSystem(config);
        auto unsafe=request;
        unsafe.at("solver").at("control_evaluation")="continuous_rhs";
        rejects([&]{(void)w::buildController(unsafe.at("controller"),unsafe.at("solver"),*system,1,.01,{.1,0});});
        // A short BDF run exercises the other solver's sampled callback path.
        auto bdf=request;bdf.at("label")="estimated-klqr-bdf";
        bdf.at("solver").at("method")="BDF";bdf.at("time").at("stop")=.04;
        const auto bdfResult=w::executeJob(bdf,root);
        check(c::readNpz(bdfResult.archives[0]).at("estimator_time").doubles().size()==9,
            "BDF repeated or skipped observer samples");
        // Strong initial estimate forces clipping; restart independent streams.
        auto clipped=request;
        clipped.at("controller").at("settings").at("control_limit")=.01;
        clipped.at("controller").at("settings").at("estimator").at("initial_lifted")
            =c::json::Value::Array{1.,0.};
        const auto setup=w::buildController(clipped.at("controller"),clipped.at("solver"),*system,1,.01,{0,0});
        const double truth[]{.1,0};
        const double gated=setup.policy.value(0,0,truth);
        const double saturated=setup.policy.value(0,.1,truth);
        setup.policy.value(0,.2,truth);
        c::NpzRecord trace;setup.appendDiagnostics(0,trace);
        check(gated==0&&std::abs(saturated)==.01,"Actuator/gating test did not exercise clipping");
        check(trace.at("previous_applied_input").doubles()[2]==saturated,
            "Clipped applied input was replaced by requested command");
        if(argc==1)std::filesystem::remove_all(root);
        std::cout<<"Koopman observer/controller/solver composition tests passed\n";return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<" (diagnostics: "<<root.string()<<")\n";return 1;
    }
}
