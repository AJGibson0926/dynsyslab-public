#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/DampedOrbit.hpp>
#include <dynsyslab/math/MelnikovTools.hpp>
#include <dynsyslab/math/OrbitGramian.hpp>
#include <dynsyslab/math/OrbitQuery.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace dynsyslab;
using namespace dynsyslab::math;
using Cx=std::complex<double>;
namespace {
const double pi=std::acos(-1.);
void check(bool ok,const std::string& m){if(!ok)throw std::runtime_error(m);}
void near(double a,double b,double tol,const std::string& m){
    if(!(std::abs(a-b)<=tol))throw std::runtime_error(m+": got "+std::to_string(a)+" expected "+std::to_string(b));}
// Complete elliptic integrals by the arithmetic-geometric mean (DLMF 19.8), parameter m=k^2.
double ellipticK(double m){double a=1,b=std::sqrt(1-m);for(int i=0;i<40;++i){const double x=.5*(a+b);b=std::sqrt(a*b);a=x;}return pi/(2*a);}
double ellipticE(double m){double a=1,b=std::sqrt(1-m),c=std::sqrt(m),sum=.5*c*c,p=1;
    for(int i=0;i<40;++i){const double x=.5*(a+b);c=.5*(a-b);b=std::sqrt(a*b);a=x;sum+=p*c*c;p*=2;}return ellipticK(m)*(1-sum);}
PlanarOrbitSystem harmonic(double w){
    PlanarOrbitSystem s;s.drift=[w](const PlanarPoint& x){return PlanarPoint{x[1],-w*w*x[0]};};
    s.hamiltonian=[w](const PlanarPoint& x){return .5*x[1]*x[1]+.5*w*w*x[0]*x[0];};s.maximumRadius=100;return s;
}
PlanarOrbitSystem pendulum(double a){
    PlanarOrbitSystem s;s.drift=[a](const PlanarPoint& x){return PlanarPoint{x[1],-a*std::sin(x[0])};};
    s.hamiltonian=[a](const PlanarPoint& x){return .5*x[1]*x[1]-a*std::cos(x[0]);};s.maximumRadius=pi;return s;
}
void orbitFamilyIdentities(){
    const auto s=harmonic(2);const auto r=periodicOrbit(s,3.);
    near(r.period,pi,1e-9,"Harmonic period");near(r.action,1.5,1e-9,"Harmonic action E/omega");
    OrbitFamilyOptions o;o.phaseSamples=64;o.maximumHarmonic=4;
    const auto samples=orbitSamples(s,r,o);const auto f=orbitFourier(samples,r.frequency,4);
    near(std::abs(f.first.coefficients.at(1)),.5*std::sqrt(2*3.)/2,1e-8,"Harmonic Fourier amplitude");
    check(std::abs(f.first.coefficients.at(2))<1e-8,"Harmonic has no second harmonic");
    near(orbitPhase(s,samples[10],o),2*pi*10/64,1e-7,"Phase of an orbit sample");
    // Same area with a uniform density: the Poincare primitive must reproduce p dq.
    auto weighted=s;weighted.density=[](const PlanarPoint&){return 1.;};near(periodicOrbit(weighted,3.).action,1.5,1e-9,"Density action");
    // Nonuniform density 1+x^2 on the unit-frequency circle: area pi R^2+pi R^4/4.
    auto curved=harmonic(1);curved.density=[](const PlanarPoint& x){return 1+x[0]*x[0];};const double R=1.4;
    near(periodicOrbit(curved,R*R/2).action,(pi*R*R+pi*R*R*R*R/4)/(2*pi),1e-9,"Nonuniform density action");
    const double a=1.7;const auto p=pendulum(a);
    for(double m:{.1,.5,.9,.99}){const double E=a*(2*m-1);const auto q=periodicOrbit(p,E);
        near(q.period,4*ellipticK(m)/std::sqrt(a),1e-8,"Pendulum period by K(m)");
        near(q.action,8*std::sqrt(a)*(ellipticE(m)-(1-m)*ellipticK(m))/pi,1e-8,"Pendulum action by E(m),K(m)");
        check(q.energyDrift<1e-9,"Pendulum energy drift");}
    std::vector<double> energies;for(int i=1;i<=24;++i)energies.push_back(a*(2*(i/25.)-1));
    auto fo=OrbitFamilyOptions{};fo.derivativeStep=2e-3;const auto table=orbitFamily(p,energies,fo);
    const auto report=frequencyBranches(table.energies,table.frequencies);check(report.monotone&&report.branches[0].direction<0,"Pendulum softening branch");
    // dOmega/dE by the exact formula from K(m), dK/dm=(E-(1-m)K)/(2m(1-m)).
    const double m=.52,K=ellipticK(m),dK=(ellipticE(m)-(1-m)*K)/(2*m*(1-m));
    const double exact=-(pi*std::sqrt(a)/2)*dK/(K*K)/(2*a);
    near(table.frequencySlopes[12],exact,1e-5,"Pendulum dOmega/dE");
    const FrequencyEvaluator w=[&](double e){return periodicOrbit(p,e).frequency;};
    const double target=.8*std::sqrt(a);const auto roots=invertFrequency(table.energies,table.frequencies,target,w);
    check(roots.size()==1&&roots[0].refined,"One inverse root");near(w(roots[0].energy),target,1e-8,"Inverse frequency");
    const auto ladder=shellLadder(table.energies,table.frequencies,2.4*std::sqrt(a),4,w);
    check(!ladder.empty(),"Shell ladder nonempty");
    for(const auto& row:ladder)near(row.harmonic*w(row.root.energy),2.4*std::sqrt(a),1e-7,"Shell ladder identity");
    const auto turning=frequencyBranches({0,1,2,3,4},{1,2,3,2,1});check(!turning.monotone&&turning.branches.size()==2,"Monotonicity report");
    check(invertFrequency({0,1,2,3,4},{1,2,3,2,1},2.5).size()==2,"Both branches inverted");
    // T=(2/sqrt a) log(32 a/delta): lambda=sqrt(a)/2 for the two-saddle pendulum orbit.
    std::vector<double> es,ts;for(int i=0;i<8;++i){const double delta=1e-7*std::pow(2.,i);es.push_back(a-delta);ts.push_back(periodicOrbit(p,a-delta).period);}
    const auto fit=fitSeparatrixPeriod(es,ts,a);near(fit.lambda,std::sqrt(a)/2,1e-4,"Separatrix fit rate");near(fit.constant/(32*a),1,1e-2,"Separatrix fit constant");
}
void dampedIdentities(){
    // Oscillating Q: mean and exact primitive.
    const int N=64;const double W=1.3;std::vector<double> q(N);for(int j=0;j<N;++j)q[j]=2+std::cos(2*pi*j/N);
    const auto avg=dampedAverage(q,W);near(avg.mean,2,1e-12,"Average Q");
    for(int j=0;j<N;++j)near(avg.chi[j],std::sin(2*pi*j/N)/W,1e-12,"chi primitive");
    // WKB transport: Omega-S'Q+Omega S' dchi/dphi=Omega0 for every phase, with S'=(Omega-Omega0)/D.
    DampedTable t;for(int i=0;i<=40;++i){const double E=.1+.05*i;t.energy.push_back(E);t.frequency.push_back(1-.1*E);t.damping.push_back(E*(1+.2*E));}
    const double w0=1.,d=.01;const auto phase=dampedPhase(t,w0,-d/2,d);
    const int i=17;const double D=t.damping[i],O=t.frequency[i],S1=phase.phaseSlope[i];
    for(double th:{.2,1.1,4.}){const double Q=D*(1+.5*std::cos(th)),dchi=(Q-D)/O;near(O-S1*Q+O*S1*dchi,w0,1e-12,"WKB phase transport");}
    near((phase.phaseFunction[i+1]-phase.phaseFunction[i-1])/(t.energy[i+1]-t.energy[i-1]),S1,1e-3,"S'=(Omega-Omega0)/D");
    near((phase.amplitude[i+1]-phase.amplitude[i-1])/(t.energy[i+1]-t.energy[i-1]),.5/D,1e-3,"a'=-Re(lambda)/(d D)");
    near(std::abs(wkbEigenfunction(phase,t.energy[i],.3,.2)),std::exp(phase.amplitude[i]),1e-9,"WKB modulus");
    // Linear damped oscillator: the Laplace average is the left projection itself.
    const double k=2,c=.3;const auto mode=linearMode(0,1,-k,-c);
    near(mode.eigenvalue.real(),-c/2,1e-14,"Linear mode real part");near(mode.eigenvalue.imag(),std::sqrt(k-c*c/4),1e-14,"Linear mode imaginary part");
    const PlanarField linear=[=](const PlanarPoint& x){return PlanarPoint{x[1],-k*x[0]-c*x[1]};};
    const PlanarPoint x0{.8,-.3};const auto psi=laplaceEigenfunction(linear,{0,0},mode,x0);
    const Cx direct=mode.left[0]*x0[0]+mode.left[1]*x0[1];check(std::abs(psi.value-direct)<1e-6*std::abs(direct),"Laplace limit of a linear flow");
    check(std::abs(isostableAmplitude(psi.value)-std::abs(psi.value))==0&&std::abs(isochronPhase(psi.value)-std::arg(psi.value))==0,"Phase and amplitude");
    for(int n:{1,2,5})near(std::abs(annularEigenfunction(psi.value,mode.eigenvalue,n,1.7)),1,1e-12,"|W|=1 on real frequency");
    // Eigenvalue i w: W(x(t))=exp(i w t) W(x) along the exact linear evolution of psi.
    const double dt=.37,wf=1.7;const Cx psi2=psi.value*std::exp(mode.eigenvalue*dt);
    const Cx ratio=annularEigenfunction(psi2,mode.eigenvalue,3,wf)/annularEigenfunction(psi.value,mode.eigenvalue,3,wf);
    check(std::abs(ratio-std::exp(Cx(0,wf*dt)))<1e-11,"Annular eigenvalue");
    // Nonlinear flow: the Laplace limit must satisfy psi(x(t))=exp(lambda t) psi(x).
    const PlanarField damped=[](const PlanarPoint& x){return PlanarPoint{x[1],-std::sin(x[0])-.2*x[1]};};
    const auto dm=linearMode(0,1,-1,-.2);LaplaceOptions lo;lo.settleRadius=1e-7;const PlanarPoint z0{1.1,.4};
    const FlowField flow=[&](const std::vector<double>& x,std::vector<double>& dx){const auto f=damped({x[0],x[1]});dx[0]=f[0];dx[1]=f[1];};
    std::vector<double> z{z0[0],z0[1]};const double elapsed=1.3;integrateFlow(flow,z,elapsed,FlowOptions{});
    const Cx base=laplaceEigenfunction(damped,{0,0},dm,z0,lo).value,moved=laplaceEigenfunction(damped,{0,0},dm,{z[0],z[1]},lo).value;
    check(std::abs(moved-std::exp(dm.eigenvalue*elapsed)*base)<2e-5*std::abs(moved),"Nonlinear Laplace eigenfunction property");
    near(shellWeight(.02,3,2,-.5),std::sqrt(2*pi*.02*3/(2*.5)),1e-14,"Shell weight");
}
void melnikovIdentities(){
    const double a=1.7,root=std::sqrt(a);ConnectionProblem p;
    p.drift=[a](const PlanarPoint& x){return PlanarPoint{x[1],-a*std::sin(x[0])};};
    p.midpoint={0,2*root};p.forwardSaddle={pi,0};p.backwardSaddle={-pi,0};p.forwardRate=p.backwardRate=root;
    const auto damping=dampingIntegral(p,[](const PlanarPoint& x){return x[1]*x[1];});
    near(damping.value.real(),8*root,1e-7,"Pendulum damping integral 8 sqrt(a)");
    const auto wedge=melnikovWedge(p,[](const PlanarPoint& x){return PlanarPoint{0,-x[1]};});
    near(wedge.value.real(),-8*root,1e-7,"Wedge equals minus dissipation");
    for(double w:{.5,1.,2.3}){const auto f=forcingIntegral(p,[](const PlanarPoint&){return PlanarPoint{0,1};},w);
        const double exact=2*pi/std::cosh(pi*w/(2*root));near(std::abs(f.value),exact,2e-6,"Forcing amplitude 2 pi sech");
        near(melnikovThreshold(damping.value.real(),std::abs(f.value)),8*root/exact,1e-4,"Threshold ratio");}
    const auto flux=lobeFlux(.1,8*root,1,2,1.);const double c=.1*8*root/2,angle=std::acos(c);
    near(flux.area,2*(2*std::sin(angle)-.1*8*root*angle),1e-13,"Lobe area");check(flux.transverse,"Lobe transverse");
    check(!lobeFlux(10,8*root,.1,1,1).transverse,"Lobe absent above threshold");
    const auto equal=captureProbabilities({8*root,8*root});near(equal[0],.5,1e-15,"Equal lobes");near(equal[1],.5,1e-15,"Equal lobes");
    const auto rows=shellMelnikovComparison({3},{.9},{7.},{.2},1.4);near(rows[0].ratio,1,1e-14,"Shell comparison");
    bool rejected=false;try{captureProbabilities({1,-1});}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Negative lobe rejected");
}
// Dense residual of A W + W A^H + B B^H for the diagonal A.
double lyapunovResidual(const ComplexVector& l,const ComplexVector& b,const ComplexVector& W){
    const auto n=l.size();double worst=0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
        worst=std::max(worst,std::abs(l[i]*W[i*n+j]+W[i*n+j]*std::conj(l[j])+b[i]*std::conj(b[j])));
    return worst;
}
void gramianIdentities(){
    const ComplexVector l{Cx(-.1,1),Cx(-.1,-1),Cx(-.5,3),Cx(-.5,-3)},b{Cx(1,.2),Cx(1,-.2),Cx(.3,0),Cx(.3,0)},c{Cx(.5,0),Cx(.5,0),Cx(0,.4),Cx(0,-.4)};
    const auto g=modalGramians(l,{b},{c});check(lyapunovResidual(l,b,g.controllability)<1e-13,"Cauchy Gramian solves Lyapunov");
    // H2 norm from the impulse response energy, h(t)=sum c_i b_i exp(l_i t), by Simpson.
    double energy=0;const double T=200,dt=.002;
    for(int k=0;k<=int(T/dt);++k){Cx h=0;for(std::size_t i=0;i<4;++i)h+=c[i]*b[i]*std::exp(l[i]*(k*dt));
        energy+=(k==0||k==int(T/dt)?1:(k%2?4:2))*std::norm(h);}
    near(g.h2Norm,std::sqrt(energy*dt/3),1e-5,"H2 norm");
    // Oscillator Hankel closed form against modal Gramians of the partial fractions.
    const double k=2.5,d=.4;const Cx root=std::sqrt(Cx(d*d-4*k)),a1=(-d+root)/2.,a2=(-d-root)/2.;
    const auto m=modalGramians({a1,a2},{{Cx(1),Cx(1)}},{{1./(a1-a2),-1./(a1-a2)}});const auto closed=oscillatorHankel(k,d);
    near(m.hankelSingularValues[0],closed[0],1e-10,"Oscillator sigma plus");near(m.hankelSingularValues[1],closed[1],1e-9,"Oscillator sigma minus");
    // Finite horizon: exact limit tau and convergence to the infinite Gramian.
    const auto finite=modalGramians({Cx(0,1)},{{Cx(2)}},{{Cx(1)}},3.);near(finite.controllability[0].real(),4*3.,1e-13,"Vanishing denominator limit");
    const auto longTime=modalGramians(l,{b},{c},2000.);check(std::abs(longTime.controllability[1]-g.controllability[1])<1e-12,"Finite horizon limit");
    const double B=.7,sl=-.3,tau=11;const auto u=undampedOrbitGramians({B,B,sl,tau});
    near(u.controllability[0],tau*B,1e-13,"W_EE");near(u.controllability[1],sl*B*tau*tau/2,1e-12,"W_Ephi");near(u.controllability[2],sl*sl*B*tau*tau*tau/3,1e-11,"W_phiphi");
    const auto law=collocatedOrbitHankel(B,sl,tau);near(u.hankel[0],law[0],1e-9,"Collocated sigma plus");near(u.hankel[1],law[1],1e-9,"Collocated sigma minus");
    // Constant B,C,D and linear Omega have closed-form quadratures.
    DampedGramianInput in;for(int i=0;i<=20;++i){const double e=.05*i;in.energy.push_back(e);in.inputVariance.push_back(2);in.outputVariance.push_back(3);in.damping.push_back(.5);in.frequency.push_back(1+.4*e);}
    in.dampingScale=.1;in.initialEnergy=.9;in.targetEnergy=.5;in.finalEnergy=.13;const auto r=dampedOrbitGramians(in);
    const double dd=.1,D=.5,up=.4,down=.37;
    near(r.cEE,D*D/dd*2*up/(D*D*D),1e-10,"c_EE");near(r.cEphi,D/(dd*dd)*2*.4*up*up/2/(D*D*D),1e-9,"c_Ephi");
    near(r.cPhiphi,1/(dd*dd*dd)*2*.16*up*up*up/3/(D*D*D),1e-7,"c_phiphi");near(r.oPhiphi,1/dd*3*down/D,1e-10,"o_phiphi");
    near(r.oEphi,1/(dd*dd*D)*3*.4*down*down/2/D,1e-9,"o_Ephi");near(r.oEE,1/(dd*dd*dd*D*D)*3*.16*down*down*down/3/D,1e-7,"o_EE");
}
void jordanIdentities(){
    // H0=w I+alpha I^2/2, W=f(I)cos(theta)+g(I)cos(2 theta) with f=I, g=I^2/4.
    const double I0=1.2,w0=.9,alpha=.35,Om=w0+alpha*I0;
    const auto f=[](double I){return I;};const auto g=[](double I){return I*I/4;};
    const std::vector<JordanMode> modes{{1,f(I0)/2,.5},{-1,f(I0)/2,.5},{2,g(I0)/2,I0/4},{-2,g(I0)/2,I0/4}};
    const auto W=[&](double I,double th){return f(I)*std::cos(th)+g(I)*std::cos(2*th);};
    const auto Wth=[&](double I,double th){return -f(I)*std::sin(th)-2*g(I)*std::sin(2*th);};
    const auto WI=[&](double I,double th){return std::cos(th)+(I/2)*std::cos(2*th);};
    for(double t:{.4,1.7,5.}){
        const auto R=jordanImpulseResponse(modes,Om,alpha,{t})[0];double avg=0;const int M=256;const double e=1e-5;
        for(int j=0;j<M;++j){const double th=2*pi*j/M;const auto out=[&](double s){
                const double I=I0-s*Wth(I0,th),phase=th+s*WI(I0,th);return W(I,phase+(w0+alpha*I)*t);};
            avg+=(out(e)-out(-e))/(2*e);}
        near(R.real(),avg/M,2e-6,"Jordan response vs flow");near(R.imag(),0,1e-12,"Real response");
    }
    const auto residues=frozenKromResidues(modes);near(residues.at(2).imag(),2*(g(I0)/2)*(I0/4),1e-14,"Frozen residue");
    const auto lift=exactInputLift(3,modes);near(lift.at(5).imag(),3*I0/4,1e-14,"Input lift");near(lift.at(1).imag(),3*I0/4,1e-14,"Input lift negative");
}
core::json::Value pendulumExpression(double a){
    const auto q=core::AnalyticExpression::state(2,0),p=core::AnalyticExpression::state(2,1);
    return (p*p*.5-q.unary("cos")*a).toJson();
}
void queryContracts(){
    using J=core::json::Value;using A=J::Array;using O=J::Object;const double a=1.;
    check(!evaluateOrbitQuery("unknown",J(O{})).has_value(),"Unknown op is not claimed");
    const auto H=pendulumExpression(a);
    const auto family=evaluateOrbitQuery("orbit_family",J(O{{"hamiltonian",H},{"centre",A{0.,0.}},{"maximum_radius",pi},
        {"energies",A{-.5,0.,.5}},{"drive_frequency",1.6},{"maximum_shell",2.}}));
    near(family->at("orbits").array()[1].at("period").number(),4*ellipticK(.5),1e-8,"Query pendulum period");
    check(family->at("frequency_branches").at("monotone").boolean(),"Query monotone");
    const auto table=evaluateOrbitQuery("orbit_family",J(O{{"table",O{{"energies",A{0.,1.,2.}},{"frequencies",A{3.,2.,1.}}}},{"target_frequency",1.5}}));
    near(table->at("inverse").array()[0].at("energy").number(),1.5,1e-12,"Query table inverse");
    const auto dissipation=(core::AnalyticExpression::state(2,1)*core::AnalyticExpression::state(2,1)).toJson();
    const auto mel=evaluateOrbitQuery("melnikov_homoclinic",J(O{{"hamiltonian",H},{"midpoint",A{0.,2.}},{"forward_saddle",A{pi,0.}},{"backward_saddle",A{-pi,0.}},
        {"forward_rate",1.},{"backward_rate",1.},{"dissipation",dissipation},
        {"forcing",A{core::AnalyticExpression::constant(0,2).toJson(),core::AnalyticExpression::constant(1,2).toJson()}},{"forcing_frequency",1.}}));
    near(mel->at("damping_integral").at("value").array()[0].number(),8,1e-6,"Query damping");
    near(mel->at("forcing_amplitude").number(),2*pi/std::cosh(pi/2),2e-6,"Query forcing");
    const auto damped=evaluateOrbitQuery("damped_averages",J(O{{"hamiltonian",H},{"centre",A{0.,0.}},{"maximum_radius",pi},{"energies",A{-.5,0.,.5}},
        {"dissipation",dissipation},{"damping_scale",.01},{"equilibrium_frequency",1.},{"real_eigenvalue",-.005},{"phase_samples",128.}}));
    check(damped->at("orbits").array().size()==3&&damped->at("phase_function").array().size()==3,"Query damped");
    const auto cap=evaluateOrbitQuery("capture_probabilities",J(O{{"lobe_integrals",A{1.,3.}}}));near(cap->at("probabilities").array()[1].number(),.75,1e-15,"Query capture");
    const auto modal=evaluateOrbitQuery("modal_gramian",core::json::parse(R"({"eigenvalues":[[-1,0]],"inputs":[[1]],"outputs":[[1]]})"));
    near(modal->at("h2_norm").number(),std::sqrt(.5),1e-14,"Query H2");
    const auto orbit=evaluateOrbitQuery("orbit_gramian",J(O{{"stiffness",2.},{"damping",.5}}));near(orbit->at("hankel_singular_values").array()[0].number(),oscillatorHankel(2,.5)[0],1e-14,"Query oscillator");
    const auto dg=evaluateOrbitQuery("damped_orbit_gramian",J(O{{"energies",A{0.,1.,2.}},{"input_variance",A{1.,1.,1.}},{"output_variance",A{1.,1.,1.}},{"damping",A{1.,1.,1.}},
        {"frequencies",A{1.,1.,1.}},{"damping_scale",1.},{"initial_energy",2.},{"target_energy",1.},{"final_energy",0.}}));near(dg->at("c_EE").number(),1,1e-13,"Query damped Gramian");
    const auto jr=evaluateOrbitQuery("jordan_response",J(O{{"modes",A{O{{"harmonic",1.},{"amplitude",.5},{"slope",.5}}}},{"real_observable",true},{"frequency",1.},{"action_twist",.1},{"times",A{0.}},{"lift_harmonic",1.}}));
    near(jr->at("response").array()[0].array()[0].number(),0,1e-14,"Query Jordan");
    bool rejected=false;try{evaluateOrbitQuery("orbit_family",J(O{{"centre",A{0.,0.}},{"maximum_radius",1.},{"energies",A{0.,1.}}}));}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Missing Hamiltonian rejected");
}
}
int main(){try{
    orbitFamilyIdentities();std::cout<<"orbit family identities ok\n";
    dampedIdentities();std::cout<<"damped orbit identities ok\n";
    melnikovIdentities();std::cout<<"Melnikov identities ok\n";
    gramianIdentities();std::cout<<"Gramian identities ok\n";
    jordanIdentities();std::cout<<"Jordan response identities ok\n";
    queryContracts();std::cout<<"query contracts ok\n";
    std::cout<<"Orbit tools tests passed\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}}
