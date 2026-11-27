#include <dynsyslab/workflow/KoopmanFeedback.hpp>
#include <dynsyslab/workflow/ControlFactory.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/EstimatorFactory.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include "../src/KlqeProcessCovariance.hpp"

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
    c::FunctionDomain q;q.dimension=2;
    c::EigenfunctionFitResult bank{c::FunctionExpansion({2,{1,0,0,1},{0,0},{1,1}}, {1.,0.,0.,1.}, q), c::Spectrum::fromContinuousRates({1.,-2.}), {}};
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
        c::saveKoopmanModel(modelPath,model);
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
        {   // Raw input-aware sensors must use the PREVIOUS applied command.
            auto driven=request;driven.at("time").at("output_step")=.005;
            auto& sensor=driven.at("controller").at("settings").at("estimator");
            sensor.at("sensor_feedthrough")=c::json::parse(R"({"coefficients":[0.7]})");
            sensor.at("measurement_feedthrough")=c::json::parse("[0.7]");
            const auto drivenResult=w::executeJob(driven,root/"input-aware");
            for(const auto& path:drivenResult.archives){
                const auto a=c::readNpz(path);const auto x=a.at("values").doubles();
                const auto y=a.at("measurement_values").doubles(),u=a.at("previous_applied_input").doubles();
                check(x.size()==y.size()*2,"Input-aware validation clocks differ");
                double maxInput=0;for(std::size_t j=0;j<y.size();++j){
                    check(std::abs(y[j]-x[2*j]-.7*u[j])<1e-10,"Sensor uses wrong applied input or applies feedthrough twice");
                    maxInput=std::max(maxInput,std::abs(u[j]));}
                check(maxInput>.01,"Sensor timing fixture did not apply a nonzero input");
            }
            // Separate estimator path: state-scaled simulated sensor; observer
            // readout is configured independently, with no implicit inheritance.
            auto standalone=sensor;standalone.at("krom")=modelPath.generic_u8string();
            standalone.at("control_representation")="fit";
            standalone.at("sensor_feedthrough")=c::json::parse(R"({"coefficients":[0.7],"scale_observable":"position"})");
            standalone.at("measurement_feedthrough_scale_output")=0.0;
            auto physical=c::createSystem(config);
            auto estimator=w::buildEstimator(standalone,c::json::Value{},*physical,1,.01);
            const double x[]{.4,-.1};(void)estimator.observe(0,0,x,.2);
            c::NpzRecord trace;estimator.appendDiagnostics(0,trace);
            check(std::abs(trace.at("measurement_values").doubles()[0]-(.4+.7*.2*.4))<1e-14,"Standalone sensor lacks input-dependent scale");
        }
        {   // KLQE process covariance shaped from the exact fixture's own prediction residual
            c::json::Value::Array records;
            for(const auto& path:result.archives)records.push_back(path.generic_u8string());
            auto spec=c::json::parse(R"({"control_representation":"fit","process_covariance":{"add":1e-3}})");
            spec.at("process_covariance").at("from_records")=records;
            const auto shaped=w::klqeProcessCovariance(spec,model);
            check(shaped.matrix.size()==4&&shaped.matrix[1]==shaped.matrix[2]&&shaped.matrix[0]>0,"residual W symmetric positive");
            check(shaped.evidence.at("relative_residual").number()<0.05,"exact fixture predicts its own records");
            spec.at("process_covariance_input_weight")=10.0;
            const auto withInput=w::klqeProcessCovariance(spec,model);
            for(int k=0;k<4;++k)check(std::abs(withInput.matrix[k]-shaped.matrix[k]-10.0)<1e-9,"input weight adds q2 B B^T with B = (1, 1)");
            auto plain=c::json::parse(R"({"control_representation":"none","process_covariance":2})");
            check(w::klqeProcessCovariance(plain,model).matrix==std::vector<double>({2,0,0,2})&&
                  w::klqeProcessCovariance(plain,model).evidence.isNull(),"scalar W unchanged");
            plain.at("process_covariance_input_weight")=1.0;
            rejects([&]{(void)w::klqeProcessCovariance(plain,model);});
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
        // Separation (feedback = "estimated_state"): a state-feedback PID acts on the KLQE's x_hat
        // (job-level estimator, its own KROM). The first command uses x_hat = 0, not the true 0.1.
        auto separated=request;separated.at("label")="estimated-state-pid";
        auto estimator=separated.at("controller").at("settings").at("estimator");
        estimator.at("krom")=modelPath.generic_u8string();
        estimator.at("control_representation")="fit";
        separated.at("estimator")=estimator;
        separated.at("controller")=c::json::parse(R"({"kind":"pid","settings":{"Kp":[6,3],
            "feedback":"estimated_state","control_limit":0.8},"target":{"kind":"constant","state":[0,0]}})");
        (void)w::planJob(separated);
        const auto pidResult=w::executeJob(separated,root);
        for(const auto& path:pidResult.archives) {
            const auto archive=c::readNpz(path);
            const auto truth=archive.at("values").doubles();
            const auto control=archive.at("control_values").doubles();
            check(control[0]==0,"estimated_state controller read the true state at the first sample");
            check(std::hypot(truth[truth.size()-2],truth[truth.size()-1])<1e-4,
                "PID on the estimated state did not stabilize the unstable fixture");
            check(archive.metadataJson().find("klqe-estimated-state/v1")!=std::string::npos,
                "estimated_state architecture not recorded");
        }
        auto missing=separated;missing.at("estimator").at("output")="";
        rejects([&]{(void)w::planJob(missing);});
        // Shadow mode (applied_controller): the plant receives a prescribed drive; the KLQE predicts
        // with that drive and the KLQR command is archived only.
        auto shadow=request;shadow.at("label")="shadow-klqr";shadow.at("time").at("stop")=.5;
        shadow.at("controller").at("settings").at("applied_controller")=c::json::parse(R"({"kind":"waveform",
            "settings":{"frequencies":[1],"model_amplitudes":[0.2],"phases_rad":[0.3]}})");
        (void)w::planJob(shadow);
        for(const auto& path:w::executeJob(shadow,root).archives) {
            const auto archive=c::readNpz(path);
            const auto ticks=archive.at("estimator_time").doubles();
            const auto applied=archive.at("applied_input").doubles();
            const auto previous=archive.at("previous_applied_input").doubles();
            const auto command=archive.at("commanded_input").doubles();
            bool differs=false;
            for(std::size_t k=1;k<ticks.size();++k) {
                check(std::abs(applied[k]-0.2*std::sin(2*3.141592653589793*ticks[k]+0.3))<1e-12,
                      "Shadow run did not apply the reference drive");
                check(previous[k]==applied[k-1],"Shadow estimator did not see the applied drive");
                differs=differs||std::abs(command[k]-applied[k])>1e-6;
            }
            check(differs,"Shadow KLQR command was not archived separately");
            check(archive.metadataJson().find("klqe-klqr-shadow-of-applied-controller/v1")!=std::string::npos,
                  "Shadow architecture not recorded");
        }
        auto orphan=shadow;orphan.at("controller").at("kind")="lqr";
        rejects([&]{(void)w::planJob(orphan);});
        if(argc==1)std::filesystem::remove_all(root);
        std::cout<<"Koopman observer/controller/solver composition tests passed\n";return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<" (diagnostics: "<<root.string()<<")\n";return 1;
    }
}
