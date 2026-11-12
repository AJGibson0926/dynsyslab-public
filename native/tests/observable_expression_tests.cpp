#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void close(double a,double b,double tolerance=1e-11){check(std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)}),"derivative mismatch");}
void checkForceKoopmanCoordinates() {
    // A plumbing fixture, not a claim that a force is itself an eigenfunction.
    for (const auto& item : std::vector<std::pair<std::string,std::string>>{
            {"duffing","force"},{"nonlinear_pendulum","torque"},
            {"nonlinear_pendulum","torque_dimensionless"}}) {
        SystemConfiguration config; config.identifier=item.first;
        if(item.first=="nonlinear_pendulum")
            config.parameters={{"mass_kg",2.3},{"length_m",.7},{"damping_n_m_s",.19}};
        const auto system=createSystem(config);
        const auto observable=system->observableCatalog().at(item.second);
        KoopmanSINDyResult fixture; fixture.stateDimension=1; fixture.dynamicalStateDimension=2;
        fixture.featureCount=1; fixture.center={.3}; fixture.scale={1.7}; fixture.exponents={1};
        fixture.coefficients={1.}; fixture.continuousEigenvalues={-1.};
        fixture.coordinateObservableNames={item.second};
        fixture.coordinateObservableDefinition=json::dump(observable.expressionDefinition());
        const KoopmanModel model(fixture,config);
        const double state[]{.8,-.4};
        const auto j=observable.jacobian(state); const auto b=system->inputMatrix(state);
        close(model.evaluate(state)[0].real(),(observable.evaluateReal(state)[0]-.3)/1.7);
        close(model.exactControlMap(state)[0].real(),(j[0].real()*b[0]+j[1].real()*b[1])/1.7);
        const auto path=std::filesystem::temp_directory_path()/("dynsyslab-force-krom-"
            +item.second+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
        savePolynomialKoopmanModel(path,model);
        const auto restored=loadKoopmanModel(path); std::filesystem::remove(path);
        close(restored.evaluate(state)[0].real(),model.evaluate(state)[0].real());
        close(restored.exactControlMap(state)[0].real(),model.exactControlMap(state)[0].real());
        check(restored.coordinateObservable().descriptor().units==observable.descriptor().units,
              "KROM archive lost force/torque units");
    }
}
int main(){try{
    checkForceKoopmanCoordinates();
    const auto r=AnalyticExpression::state(2,0),v=AnalyticExpression::state(2,1);
    const std::vector<AnalyticExpression> oscillator{v,-r};
    const auto energy=r.power(2)+v.power(2);
    const auto third=r.power(2).flowDerivative(oscillator,3);
    for(const auto x:std::vector<std::vector<double>>{{1.,2.},{.7,-.9},{2.3,.2}}){
        close(energy.flowDerivative(oscillator).evaluate(x.data()),0);
        close(third.evaluate(x.data()),-8*x[0]*x[1]);
        close(AnalyticExpression::fromJson(third.toJson()).evaluate(x.data()),third.evaluate(x.data()));
    }
    const auto f=(r.unary("sin")*v.unary("exp")+r.unary("log")+v.unary("tanh"))/(1+r.power(2));
    const double probe[]{1.4,.3};
    for(std::size_t axis=0;axis<2;++axis){
        double a[]{probe[0],probe[1]},b[]{probe[0],probe[1]};a[axis]+=1e-5;b[axis]-=1e-5;
        close(f.derivative(axis).evaluate(probe),(f.evaluate(a)-f.evaluate(b))/2e-5,1e-8);
    }
    SystemConfiguration config;config.identifier="rpe";
    auto system=createSystem(config);const auto echo=system->observableCatalog().at("far_field_acoustic_echo");
    check(echo.hasAnalyticJacobian(),"echo must use an analytic gradient");
    const auto dt=echo.flowDerivative(system->driftExpressions());
    const auto dt2=echo.flowDerivative(system->driftExpressions(),2);
    for(std::size_t i=0;i<30;++i){
        const double x[]{.8+i*.05,std::sin(double(i))};double velocity[2];system->rhsWithInput(0,x,0,velocity);
        const auto j=echo.jacobian(x);const auto numeric=echo.jacobian(x,ObservableJacobianBackend::finiteDifference);
        for(std::size_t k=0;k<2;++k)close(j[k].real(),numeric[k].real(),1e-7);
        close(dt.evaluateReal(x)[0],j[0].real()*velocity[0]+j[1].real()*velocity[1]);
        const auto dj=dt.jacobian(x);close(dt2.evaluateReal(x)[0],dj[0].real()*velocity[0]+dj[1].real()*velocity[1]);
        close(Observable::fromExpressionDefinition(echo.expressionDefinition()).evaluateReal(x)[0],echo.evaluateReal(x)[0]);
    }
    KoopmanSINDyResult fit;fit.stateDimension=1;fit.dynamicalStateDimension=2;fit.featureCount=1;
    fit.center={0};fit.scale={1};fit.exponents={1};fit.coefficients={1.};fit.continuousEigenvalues={-1.};
    fit.coordinateObservableNames={"far_field_acoustic_echo"};
    fit.coordinateObservableDefinition=json::dump(echo.expressionDefinition());
    auto model=std::make_shared<KoopmanModel>(fit,config);
    const std::vector<double> target{1.1,.2};const auto targetEcho=echo.evaluateReal(target.data());
    KoopmanLQRController stateReference(model,{2.},{.1},[target](double){return target;});
    KoopmanLQRController echoReference(model,{2.},{.1},[targetEcho](double){return targetEcho;},
        "exact",KoopmanTargetSpace::observableCoordinates);
    close(stateReference.value(0,probe)[0],echoReference.value(0,probe)[0]);
    const auto b=model->exactControlMap(probe);const auto j=echo.jacobian(probe);const auto input=system->inputMatrix(probe);
    close(b[0].real(),j[0].real()*input[0]+j[1].real()*input[1]);
    const auto file=std::filesystem::temp_directory_path()/("dynsyslab-expression-"
        +std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".npz");
    savePolynomialKoopmanModel(file,*model);
    const auto restored=loadKoopmanModel(file);
    std::filesystem::remove(file);
    close(restored.evaluate(probe)[0].real(),model->evaluate(probe)[0].real());
    std::cout<<"analytic observable, repeated flow derivative, echo target, and archive checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
