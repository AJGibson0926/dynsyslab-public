#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/Estimation.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace c=dynsyslab::core;
void check(bool b,const char* message){if(!b)throw std::runtime_error(message);}
void close(double a,double b){check(std::abs(a-b)<1e-12,"Runtime adapter changed feedback");}
template<class F> void rejects(F f){try{f();}catch(const std::logic_error&){return;}throw std::runtime_error("Invalid runtime input accepted");}
}
int main(){try{
    c::FunctionDomain q;q.dimension=2;
    c::EigenfunctionFitResult bank{c::FunctionExpansion({2,{1,0,0,1},{0,0},{1,1}}, {1.,0.,0.,1.}, q), c::Spectrum::fromContinuousRates({-1.,-2.}), {}};
    auto model=std::make_shared<c::KoopmanModel>(bank,1);
    model->setControlMap("constant",c::ConstantControlMap{{1.,.5},2,1,"regression"});
    const std::vector<double> Q{2,0,0,3},R{.2};
    c::KoopmanLQRController control(model,Q,R,[](double){return std::vector<double>{.1,.2};},
        "constant",c::KoopmanTargetSpace::observableCoordinates);
    const auto expected=c::solveContinuousLQR({-1,0,0,-2},{1,.5},Q,R,2,1);
    const double coordinates[]{.4,-.3};
    close(control.valueFromCoordinates(0,coordinates)[0],-expected.gain[0]*.3+expected.gain[1]*.5);
    close(control.valueFromLifted(0,{.4,-.3})[0],control.valueFromCoordinates(0,coordinates)[0]);
    check(control.gainSolveCount()==1,"Constant KLQR repeatedly solved CARE");
    rejects([&]{(void)control.value(0,coordinates);});
    rejects([&]{(void)control.valueFromLifted(0,{.4});});
    model->setControlMap("constant",c::ConstantControlMap{{2.,1.},2,1,"replacement"});
    (void)control.valueFromCoordinates(0,coordinates);
    check(control.gainSolveCount()==2,"Constant map replacement left a stale gain");
    model->setControlMap("bilinear",c::BilinearControlMap{{.1,0.,0.,.2},{1.,.5},{},2,1,"regression"});
    c::KoopmanLQRController changing(model,Q,R,{},"bilinear");
    (void)changing.valueFromCoordinates(0,coordinates);(void)changing.valueFromCoordinates(.1,coordinates);
    check(changing.gainSolveCount()==2,"Bilinear map was incorrectly cached as constant");
    c::KoopmanMPCController mpc(model,Q,R,{},10,.01,"constant");
    close(mpc.valueFromCoordinates(0,coordinates)[0],mpc.valueFromLifted(0,{.4,-.3})[0]);
    model->setDiscreteInputMap("sampled",{{1.,.5},2,1,.01,"sampled regression"});
    rejects([&]{c::KoopmanLQRController wrong(model,Q,R,{},"sampled");});
    model->setReconstruction("sensor",{{1.,1.},1,2,0,0,0});
    const auto design=c::designKoopmanLQEFromMeasurements(*model,"sensor",Q,{.1});
    const std::vector<double> time{0,.1,.2,.3},measurement{1,.9,.8,.7};
    const auto run=c::runKoopmanLQEFromLifted(*model,design,time,measurement,{0,0});
    check(run.stateEstimate.empty()&&run.stateDimension==0,"Measurement observer invented physical states");
    check(run.liftedEstimate.size()==8&&run.measurementEstimate.size()==4,"Measurement-only observer shape wrong");
    check(run.liftedEstimate[6]>0,"Observer did not respond to measurements");
    auto permuted=model->extract({1,0});
    rejects([&]{(void)c::runKoopmanLQEFromLifted(permuted,design,time,measurement,{0,0});});
    c::SystemConfiguration s;s.identifier="duffing";s.preset="double_well_reference";
    auto stateCoordinates=c::FunctionDomain::physicalState(q.dimension);
    auto stateModel=std::make_shared<c::KoopmanModel>(bank.functions.withCoordinates(stateCoordinates),bank.spectrum,s);
    stateModel->setControlMap("constant",c::ConstantControlMap{{1.,.5},2,1,"regression"});
    c::KoopmanLQRController oldAdapter(stateModel,Q,R,{},"constant");
    close(oldAdapter.value(0,coordinates)[0],oldAdapter.valueFromCoordinates(0,coordinates)[0]);
    stateModel->setReconstruction("sensor",{{1.,1.},1,2,0,0,0});
    stateModel->setReconstruction("state",{{1.,0.,0.,1.},2,2,0,0,0});
    const auto oldDesign=c::designKoopmanLQE(*stateModel,"sensor","state",Q,{.1});
    const auto oldRun=c::runKoopmanLQE(*stateModel,oldDesign,time,measurement,{0,0});
    const auto newRun=c::runKoopmanLQEFromLifted(*stateModel,oldDesign,time,measurement,{0,0});
    check(oldRun.liftedEstimate==newRun.liftedEstimate&&oldRun.stateEstimate==newRun.stateEstimate,"Observer adapter changed numerical results");
    std::cout<<"Coordinate and lifted runtime tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
