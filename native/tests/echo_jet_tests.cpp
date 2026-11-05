#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/System.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double tolerance=1e-10){check(std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)}),"echo jet mismatch");}
int main(){try{
    SystemConfiguration config;config.identifier="rpe";config.parameters["kappa"]=4./3.;
    const auto system=createSystem(config);const auto catalog=system->observableCatalog();
    const auto h=catalog.at("echo"),d=catalog.at("echo_dot"),dd=catalog.at("echo_ddot");
    const auto dr=catalog.at("echo_dot_dr"),dv=catalog.at("echo_dot_dv");
    const auto jet=catalog.coordinatePreset("echo_jet2");
    check(jet.outputDimension()==2&&jet.hasAnalyticJacobian(),"jet must have two analytic coordinates");
    for(const auto& x:std::vector<std::vector<double>>{{1.,0.},{1.,.24147361152683358},{1.2,.3},{.8,-.4},{2.,1.}}){
        const auto hJac=h.jacobian(x.data()),dJac=d.jacobian(x.data());
        double f[2];system->rhsWithInput(0,x.data(),0,f);
        close(d.evaluateReal(x.data())[0],hJac[0].real()*f[0]+hJac[1].real()*f[1]);
        close(dd.evaluateReal(x.data())[0],dJac[0].real()*f[0]+dJac[1].real()*f[1]);
        close(dr.evaluateReal(x.data())[0],dJac[0].real());
        close(dv.evaluateReal(x.data())[0],dJac[1].real());
        for(std::size_t axis=0;axis<2;++axis){
            auto plus=x,minus=x;plus[axis]+=1e-5;minus[axis]-=1e-5;
            close(dJac[axis].real(),(d.evaluateReal(plus.data())[0]-d.evaluateReal(minus.data())[0])/2e-5,2e-7);
        }
        const auto push=jet.pushForwardRealBatch(x,{f[0],f[1]},1);
        close(push[0],d.evaluateReal(x.data())[0]);close(push[1],dd.evaluateReal(x.data())[0]);
        // A local second-order ODE step independently verifies the second time derivative.
        const auto jf=system->jacobian(0,x.data());const double eps=1e-4;
        auto plus=x,minus=x;
        for(std::size_t i=0;i<2;++i){const double a=jf[2*i]*f[0]+jf[2*i+1]*f[1];
            plus[i]+=eps*f[i]+.5*eps*eps*a;minus[i]+=-eps*f[i]+.5*eps*eps*a;}
        close((h.evaluateReal(plus.data())[0]-2*h.evaluateReal(x.data())[0]+h.evaluateReal(minus.data())[0])/(eps*eps),dd.evaluateReal(x.data())[0],2e-5);
        const auto b=system->inputMatrix(x.data());const double u=.3;
        double forced[2];system->rhsWithInput(0,x.data(),u,forced);
        close(hJac[0].real()*forced[0]+hJac[1].real()*forced[1],
            d.evaluateReal(x.data())[0]+(hJac[0].real()*b[0]+hJac[1].real()*b[1])*u);
    }
    const std::vector<double> probe{1.2,.3};const auto g=system->inputExpressions();
    const auto drift=system->driftExpressions();const double u=.3,udot=-.4;
    const auto withInput=h.inputFlowDerivative(drift,g,{u},{udot},1);
    const auto heldSecond=h.inputFlowDerivative(drift,g,{u},{},2);
    const auto changingSecond=h.inputFlowDerivative(drift,g,{u},{udot},2);
    double f[2];system->rhsWithInput(0,probe.data(),u,f);
    const auto jh=h.jacobian(probe.data());const auto inputProbe=system->inputMatrix(probe.data());
    close(withInput.evaluateReal(probe.data())[0],jh[0].real()*f[0]+jh[1].real()*f[1]);
    close(changingSecond.evaluateReal(probe.data())[0]-heldSecond.evaluateReal(probe.data())[0],
        (jh[0].real()*inputProbe[0]+jh[1].real()*inputProbe[1])*udot);
    const double step=1e-4;auto plusTime=probe,minusTime=probe;
    for(std::size_t axis=0;axis<2;++axis){double acceleration=inputProbe[axis]*udot;
        for(std::size_t column=0;column<2;++column)
            acceleration+=(drift[axis].derivative(column).evaluate(probe.data())
                +u*g[axis].derivative(column).evaluate(probe.data()))*f[column];
        plusTime[axis]+=step*f[axis]+.5*step*step*acceleration;
        minusTime[axis]+=-step*f[axis]+.5*step*step*acceleration;
    }
    close((h.evaluateReal(plusTime.data())[0]-2*h.evaluateReal(probe.data())[0]+h.evaluateReal(minusTime.data())[0])/(step*step),changingSecond.evaluateReal(probe.data())[0],2e-5);
    bool rejected=false;try{(void)jet.pushForwardRealBatch({1,0},{0},1);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"push-forward shape mismatch was accepted");
    KoopmanSINDyResult fit;fit.stateDimension=2;fit.dynamicalStateDimension=2;fit.featureCount=3;
    fit.center={0,0};fit.scale={1,1};fit.exponents={2,0,1,1,0,1};
    fit.coefficients={.3,.2,.5};fit.continuousEigenvalues={-1.};
    fit.coordinateObservableNames=jet.descriptor().outputNames;
    fit.coordinateObservableDefinition=json::dump(jet.expressionDefinition());
    const KoopmanModel model(fit,config);const double x[]{1.2,.3};
    const auto y=jet.evaluateReal(x);close(model.evaluate(x)[0].real(),.3*y[0]*y[0]+.2*y[0]*y[1]+.5*y[1]);
    const auto b=system->inputMatrix(x);const double eps=1e-7;
    const double plus[]{x[0]+eps*b[0],x[1]+eps*b[1]},minus[]{x[0]-eps*b[0],x[1]-eps*b[1]};
    close(model.exactControlMap(x)[0].real(),(model.evaluate(plus)[0].real()-model.evaluate(minus)[0].real())/(2*eps),1e-7);
    const auto file=std::filesystem::temp_directory_path()/("echo-jet-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    savePolynomialKoopmanModel(file,model);const auto restored=loadKoopmanModel(file);std::filesystem::remove(file);
    close(restored.evaluate(x)[0].real(),model.evaluate(x)[0].real());
    close(restored.exactControlMap(x)[0].real(),model.exactControlMap(x)[0].real());
    std::cout<<"echo, echo-dot, echo-double-dot, spatial derivatives and lifted input map passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
