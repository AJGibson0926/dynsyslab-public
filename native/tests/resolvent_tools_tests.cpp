// Mathematical identities for the Koopman resolvent tools. No study data.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/DampedOrbit.hpp>
#include <dynsyslab/math/EigenvalueLattice.hpp>
#include <dynsyslab/math/FrequencyGramian.hpp>
#include <dynsyslab/math/GeneralizedLaplace.hpp>
#include <dynsyslab/math/KoopmanFrequencyResponse.hpp>
#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <dynsyslab/math/LaplaceTransform.hpp>
#include <dynsyslab/math/OrbitGramian.hpp>
#include <dynsyslab/math/ResolventQuery.hpp>
#include <dynsyslab/math/ShellResolvent.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace dynsyslab;
using namespace dynsyslab::math;
using Cx=std::complex<double>;
namespace {
const double pi=std::acos(-1.);
const Cx I(0,1);
int checks=0;
void check(bool ok,const std::string& m){++checks;if(!ok)throw std::runtime_error(m);}
void near(double a,double b,double tol,const std::string& m){
    ++checks;if(!(std::abs(a-b)<=tol))throw std::runtime_error(m+": got "+std::to_string(a)+" expected "+std::to_string(b));}
void nearC(Cx a,Cx b,double tol,const std::string& m){
    ++checks;if(!(std::abs(a-b)<=tol)){char b2[64];std::snprintf(b2,64,"%.3e",std::abs(a-b));throw std::runtime_error(m+": error "+b2);}}
template<class F>bool rejects(F&& f){try{f();}catch(const std::invalid_argument&){return true;}catch(const std::domain_error&){return true;}return false;}
SampledSignal sampled(double T,std::size_t n,const std::function<Cx(double)>& f){
    SampledSignal g;for(std::size_t k=0;k<n;++k){const double t=T*double(k)/double(n-1);g.times.push_back(t);g.values.push_back(f(t));}return g;
}
ComplexMatrix mat(std::size_t r,std::size_t c,std::initializer_list<Cx> v){ComplexMatrix m(r,c);m.data.assign(v);return m;}
void linearAlgebra(){
    const auto A=mat(3,3,{0,1,0, 0,0,1, 6,-11,6}); // Companion of (x-1)(x-2)(x-3).
    auto e=complexEigenvalues(A);std::sort(e.begin(),e.end(),[](Cx a,Cx b){return a.real()<b.real();});
    for(int k=0;k<3;++k)nearC(e[k],double(k+1),1e-9,"Companion eigenvalue");
    const auto rot=mat(2,2,{Cx(0,2),1, 0,Cx(-1,-3)});auto r=complexEigenvalues(rot);
    check(std::abs(r[0]-Cx(0,2))<1e-12||std::abs(r[1]-Cx(0,2))<1e-12,"Triangular complex eigenvalue");
    const auto sv=complexSingularValues(mat(2,2,{3,0, 0,Cx(0,-4)}));near(sv[0],4,1e-13,"Largest singular value");near(sv[1],3,1e-13,"Smallest singular value");
    const auto M=mat(2,2,{Cx(1,1),2, 3,Cx(0,-1)});const ComplexLu lu(M);
    const ComplexVector b{1,Cx(0,1)};const auto x=lu.solve(b);const auto back=product(M,x);nearC(back[0],b[0],1e-14,"LU solve");
    const auto y=lu.solveAdjoint(b);const auto ad=product(adjoint(M),y);nearC(ad[1],b[1],1e-14,"Adjoint solve");
}
void laplaceIdentities(){
    const Cx lambda(-.4,3);const double T=12;const auto f=[&](double t){return std::exp(lambda*t);};
    const auto exact=[&](Cx s){return (1.-std::exp((lambda-s)*T))/(s-lambda);};
    const Cx highFrequency(.2,120);auto g=sampled(T,1201,f);
    nearC(laplaceWindowIntegral(g,highFrequency,LaplaceInterpolant::CubicHermite),exact(highFrequency),2e-8,"Cubic Filon at |s|h=1.2");
    g.derivatives.clear();for(double t:g.times)g.derivatives.push_back(lambda*f(t));
    nearC(laplaceWindowIntegral(g,highFrequency,LaplaceInterpolant::CubicHermite),exact(highFrequency),2e-8,"Cubic with exact derivatives");
    const auto c1=sampled(T,201,f),c2=sampled(T,401,f);const Cx s0(.3,1);
    const double e1=std::abs(laplaceWindowIntegral(c1,s0,LaplaceInterpolant::Linear)-exact(s0));
    const double e2=std::abs(laplaceWindowIntegral(c2,s0,LaplaceInterpolant::Linear)-exact(s0));
    near(std::log2(e1/e2),2,.05,"Linear panels converge at second order");
    // Nonuniform grid.
    SampledSignal u;for(int k=0;k<=1600;++k){const double x=k/1600.;const double t=T*x*x*(3-2*x);u.times.push_back(t);u.values.push_back(f(t));}
    nearC(laplaceWindowIntegral(u,s0,LaplaceInterpolant::CubicHermite),exact(s0),1e-8,"Nonuniform grid");
    LaplaceTransformOptions o;o.growth=GrowthBound{lambda.real(),1,true};
    const Cx s1(.1,2);const auto r=laplaceTransform(sampled(T,2001,f),s1,o);check(r.insideRegion,"Inside the region");
    const double trueTail=std::abs(std::exp((lambda-s1)*T)/(s1-lambda));
    check(trueTail<=*r.tailBound*(1+1e-12),"Tail bound holds");near(*r.tailBound,std::exp((lambda.real()-s1.real())*T)/(s1.real()-lambda.real()),1e-14,"Tail bound formula");
    check(!laplaceTransform(g,Cx(-1,0),o).insideRegion,"Left of the abscissa");
    near(fitGrowthBound(sampled(30,3001,[](double t){return std::exp(-.7*t)*std::cos(5*t);})).rate,-.7,.03,"Fitted growth");
    near(fitGrowthBound(sampled(60,6001,[](double t){return Cx(std::cos(3*t));})).rate,0,.01,"Bounded signal has zero rate");
    check(rejects([&]{laplaceTransform(SampledSignal{{0,1,1},{1,1,1},{}},Cx(1));}),"Non-increasing time rejected");
    // z-transform of a geometric sequence and the exact s<->z map.
    const Cx a(.6,.3),z(1.2,.5);ComplexVector y;Cx p=1;for(int k=0;k<40;++k){y.push_back(p);p*=a;}
    nearC(zTransform(y,z),(1.-std::pow(a/z,40))/(1.-a/z),1e-13,"Geometric z-transform");
    const Cx s2(-.3,2.1);const double h=.25;nearC(laplaceFromZ(zFromLaplace(s2,h),h),s2,1e-13,"Principal s<->z");
    const Cx s3(.2,30);const int branch=int(std::round(s3.imag()*h/(2*pi)));nearC(laplaceFromZ(zFromLaplace(s3,h),h,branch),s3,1e-12,"Branch s<->z");
    // Impulse-sampled transform: h sum y_k e^{-s k h} = h X(e^{sh}).
    Cx direct{};for(std::size_t k=0;k<y.size();++k)direct+=y[k]*std::exp(-s2*(double(k)*h));nearC(direct,zTransform(y,zFromLaplace(s2,h)),1e-12,"Sampled Laplace is X(e^{sh})");
}
void inverseIdentities(){
    const auto one=[](Cx s){return 1./(s+1.);};
    nearC(talbotInverse(one,2).value,std::exp(-2.),1e-12,"Talbot exponential");
    InverseLaplaceOptions fine;fine.nodes=48;nearC(talbotInverse([](Cx s){return 1./(s*s+1.);},2,fine).value,std::sin(2.),1e-11,"Talbot sine");
    nearC(talbotInverse([](Cx s){return 1./std::sqrt(s);},1.5).value,1/std::sqrt(pi*1.5),1e-11,"Talbot branch point s^{-1/2}");
    const Cx p(-.2,3);InverseLaplaceOptions o;o.nodes=48;o.singularities={p};const double t=1.7;
    const auto r=talbotInverse([&](Cx s){return 1./(s-p);},t,o);check(r.scopeVerified,"Talbot encloses the pole");nearC(r.value,std::exp(p*t),1e-11,"Talbot complex pole");
    InverseLaplaceOptions outside;outside.singularities={Cx(0,500)};check(!talbotInverse([](Cx s){return 1./(s-Cx(0,500));},3,outside).scopeVerified,"Talbot flags an excluded pole");
    const Cx q(0,20);InverseLaplaceOptions b;b.singularities={q};const auto v=bromwichInverse([&](Cx s){return 1./(s-q);},3,b);
    check(v.scopeVerified,"Bromwich right of the pole");nearC(v.value,std::exp(q*3.),2e-6,"Bromwich oscillatory pole");
    nearC(bromwichInverse(one,1.2).value,std::exp(-1.2),2e-6,"Bromwich exponential");
    InverseLaplaceOptions shifted;shifted.shift=3;shifted.singularities={Cx(2,0)};
    nearC(bromwichInverse([](Cx s){return 1./(s-2.);},1,shifted).value,std::exp(2.),3e-5,"Bromwich with a growing mode");
    check(rejects([&]{talbotInverse(one,0);}),"t=0 rejected");check(rejects([&]{InverseLaplaceOptions x;x.nodes=7;talbotInverse(one,1,x);}),"Odd nodes rejected");
}
void poleModels(){
    PoleResidueModel m;m.poles.push_back({Cx(-.5,2),{Cx(1,-.5),Cx(2,1)}});m.poles.push_back({Cx(-1,-1),{Cx(.3,.2)}});
    const Cx s(.7,1.1);nearC(evaluateResolvent(m,s),Cx(1,-.5)/(s-m.poles[0].eigenvalue)+Cx(2,1)/std::pow(s-m.poles[0].eigenvalue,2)+Cx(.3,.2)/(s-m.poles[1].eigenvalue),1e-14,"Resolvent evaluation");
    const double t=1.3;nearC(evaluateTimeResponse(m,t),(Cx(1,-.5)+Cx(2,1)*t)*std::exp(m.poles[0].eigenvalue*t)+Cx(.3,.2)*std::exp(m.poles[1].eigenvalue*t),1e-14,"Jordan time response");
    InverseLaplaceOptions o;for(const auto& p:m.poles)o.singularities.push_back(p.eigenvalue);o.nodes=48;
    const auto inv=talbotInverse([&](Cx x){return evaluateResolvent(m,x);},t,o);check(inv.scopeVerified,"Enclosed double pole");nearC(inv.value,evaluateTimeResponse(m,t),1e-10,"Double-pole inverse");
    near(modelAbscissa(m),-.5,1e-15,"Abscissa");
    PoleResidueModel strongStable;strongStable.poles.push_back({Cx(-.1,0),{0}});strongStable.poles.push_back({Cx(-2,0),{1}});
    near(modelAbscissa(strongStable),-2,1e-15,"Vanishing residue moves the abscissa left");
    const auto e=cauchyLaurent([&](Cx x){return evaluateResolvent(m,x);},m.poles[0].eigenvalue,.5);
    check(e.converged,"Contour converged");check(e.order==2,"Contour detects a double pole");nearC(e.laurent[0],Cx(1,-.5),1e-10,"Contour residue");nearC(e.laurent[1],Cx(2,1),1e-10,"Contour second coefficient");
    const auto simple=cauchyLaurent([&](Cx x){return evaluateResolvent(m,x);},m.poles[1].eigenvalue,.4);check(simple.order==1,"Contour detects a simple pole");nearC(simple.laurent[0],Cx(.3,.2),1e-10,"Simple residue");
    const auto regular=cauchyLaurent([](Cx x){return std::exp(x);},Cx(0,0),1);near(std::abs(regular.laurent[0]),0,1e-13,"Analytic interior has no residue");
    const auto third=cauchyLaurent([](Cx x){return 1./(x*x*x);},Cx(0),1);check(third.higherOrderRatio>1,"Order-three pole flagged");
}
void lattices(){
    const auto eq=equilibriumPoleLattice({Cx(-1,1),Cx(-1,-1),Cx(-3,0)},2);check(eq.poles.size()==10,"C(5,2) points up to order 2");
    near(eq.poles[0].value.real(),0,1e-15,"Constant mode first");
    bool found=false;for(const auto& p:eq.poles)if(p.k==std::vector<int>{1,1,0})found=std::abs(p.value-Cx(-2,0))<1e-14;check(found,"k.lambda");
    check(equilibriumPoleLattice({Cx(-1),Cx(-2)},2).coincidences.size()>=1,"Resonance lambda2=2 lambda1 flagged");
    check(rejects([]{equilibriumPoleLattice({Cx(.1)},2);}),"Unstable eigenvalue rejected");
    const double Om=1.7;const auto a=limitCyclePoleLattice({Cx(-.3,.4)},Om,2,6),b=limitCyclePoleLattice({Cx(-.3,.4+Om)},Om,2,6);
    // nu is defined modulo i Omega: the lattices agree on their common harmonic range.
    for(const auto& p:a.poles){if(std::abs(p.value.imag())>3.5*Om)continue;
        check(std::any_of(b.poles.begin(),b.poles.end(),[&](const LatticePole& x){return std::abs(x.value-p.value)<1e-12;}),"Lattice invariant under nu->nu+i Omega");}
    const auto torus=torusPoleLattice({Cx(-.2)},{1.,std::sqrt(2.)},1,2);check(torus.poles.size()==50,"Torus truncation size");
    near(*torus.smallestHarmonicDivisor,std::sqrt(2.)-1,1e-12,"Small divisor |sqrt2-1|");
    // Reuse the suite's lattice diagnostic: generated points are flagged as combinations.
    EigenvalueLatticeOptions o;o.principalIndices={0,1};o.maximumOrder=3;o.tolerance=1e-9;
    const auto report=eigenvalueLattice({Cx(-.7,1.3),Cx(-1.1,-.4),Cx(-1.4,2.6),Cx(-1.8,.9)},o);
    check(report.entries[2].role==LatticeRole::Combination&&report.entries[3].role==LatticeRole::Combination,"Lattice diagnostic agrees");
}
void pseudospectrum(){
    const auto D=mat(2,2,{Cx(-1,2),0, 0,Cx(-3,-1)});const Cx s(.5,.5);
    near(1/resolventNorm(D,s),std::min(std::abs(s-Cx(-1,2)),std::abs(s-Cx(-3,-1))),1e-12,"Normal: sigma_min is the distance to the spectrum");
    const auto J=mat(2,2,{0,1, 0,0});const double a=.3;const double rn=resolventNorm(J,a);
    const double p=2*a*a+1,g=(p+std::sqrt(p*p-4*a*a*a*a))/2;near(rn,std::sqrt(g)/(a*a),1e-10,"Jordan block resolvent norm");
    const auto grid=pseudospectrumGrid(D,{-1.,0.},{2.,0.});check(grid.sigmaMin.size()==4,"Grid size");near(grid.sigmaMin[0],0,1e-12,"sigma_min vanishes on an eigenvalue");
}
void generalizedLaplace(){
    const double Om=1.3;const Cx nu(-.35,Om);const auto f=[&](double t){
        return 2.+Cx(.7,.2)*std::exp(I*Om*t)+Cx(.7,-.2)*std::exp(-I*Om*t)+Cx(.4,.1)*std::exp(nu*t)+Cx(-.2,.3)*std::exp(Cx(-.9,0)*t);};
    const auto g=sampled(120,12001,f);
    const auto r=generalizedLaplaceAnalysis(g,{Cx(-.9,0),nu,Cx(0,Om),Cx(0,-Om),Cx(0,0)});
    check(r.modes[0].eigenvalue.real()==0,"Slow modes first");check(r.modes.back().eigenvalue==Cx(-.9,0),"Fast mode last");
    const auto coef=[&](Cx l){for(const auto& m:r.modes)if(std::abs(m.eigenvalue-l)<1e-12)return m.coefficient;throw std::runtime_error("Missing mode");};
    nearC(coef(0),2,1e-5,"GLA constant");nearC(coef(Cx(0,Om)),Cx(.7,.2),1e-5,"GLA stationary mode");nearC(coef(nu),Cx(.4,.1),1e-4,"GLA decaying mode");nearC(coef(Cx(-.9,0)),Cx(-.2,.3),1e-3,"GLA fastest mode");
    check(r.residualRms<1e-3*r.signalRms,"GLA residual");
    check(rejects([&]{generalizedLaplaceAnalysis(g,{Cx(0,1),Cx(0,1)});}),"Repeated eigenvalues rejected");
    // Stationary/non-stationary split and Floquet decay.
    const double w=2.,fl=-.25;const auto y=[&](double t){return Cx(1.5*std::cos(w*t)+.3*std::sin(2*w*t))+std::exp(fl*t)*Cx(1,.5)*std::exp(I*w*t);};
    const auto d=decomposeStationary(sampled(160,32001,y),w,4,.6);
    nearC(d.harmonics.at(1),.75,1e-6,"Stationary harmonic");nearC(d.harmonics.at(2),Cx(0,-.15),1e-6,"Second harmonic");
    check(d.decayRateAvailable,"Decay fit ran");near(d.decayRate,fl,2e-3,"Stroboscopic decay equals the exponent");
    check(probePoleOrder(sampled(200,20001,[&](double t){return std::exp(Cx(0,2)*t)+.3;}),Cx(0,2)).order==1,"Data probe: simple pole");
    const auto dbl=probePoleOrder(sampled(200,20001,[&](double t){return (1+t)*std::exp(Cx(0,2)*t);}),Cx(0,2));check(dbl.order==2,"Data probe: double pole");near(dbl.ratio,3,.05,"Double-pole ratio");
}
BilinearLiftedModel paperLift(double a1,double a2){
    BilinearLiftedModel m;m.A=mat(3,3,{a1,0,1, 0,a2,0, 0,0,2*a2});m.B=mat(3,1,{0,1,0});
    m.C=mat(1,3,{1,0,0});m.N.push_back(mat(3,3,{0,0,0, 0,0,0, 0,2,0}));return m;
}
void frequencyResponse(){
    const double a1=-1.3,a2=-.7,w=1.9;const auto m=paperLift(a1,a2);HarmonicResponseOptions o;o.maximumOrder=5;
    const auto r=harmonicResponse(m,w,o);check(r.valid&&r.hurwitz,"Stable example is valid");
    const Cx H2=1./((Cx(0,2*w)-a1)*std::pow(Cx(0,w)-a2,2));
    near(std::abs(r.outputs[0][0]),0,1e-14,"H1(x1)=0");nearC(r.outputs[1][0],H2,1e-14,"Paper H2(x1)");
    for(int n=3;n<=5;++n)near(std::abs(r.outputs[n-1][0]),0,1e-14,"Higher harmonics vanish");
    nearC(r.states[0][1],1./(Cx(0,w)-a2),1e-14,"H1(x2)");near(std::abs(r.states[1][1]),0,1e-14,"H2(x2)=0");
    // Six-state lift (eq. 11): transfer function from z6=u^2 evaluated at 2iw.
    const Cx s(0,2*w);nearC(2./((s-a1)*(s-2*a2)*(s-(a2+Cx(0,w)))),H2,1e-14,"Six-state lift of eq. (11)");
    // Eigenfunction (10) and the Appendix I mode: the u^2 coefficient of x1 equals H2.
    const Cx h=1./(Cx(0,w)-a2),alpha=1./(a1-2*a2),beta=2./((a1-2*a2)*(a1-a2-Cx(0,w))),gamma=2./((a1-2*a2)*(a1-a2-Cx(0,w))*(a1-Cx(0,2*w)));
    nearC(-alpha*h*h-beta*h-gamma,H2,1e-13,"Appendix I mode equals H2");
    // The paper's a1!=a2 is not needed; a1=2a2 only loses the eigenfunction (10), not H2.
    const auto resonantLift=harmonicResponse(paperLift(2*a2,a2),w,o);nearC(resonantLift.outputs[1][0],1./((Cx(0,2*w)-2*a2)*std::pow(Cx(0,w)-a2,2)),1e-14,"H2 at a1=2a2");
    const auto equal=harmonicResponse(paperLift(a2,a2),w,o);nearC(equal.outputs[1][0],1./((Cx(0,2*w)-a2)*std::pow(Cx(0,w)-a2,2)),1e-14,"H2 at a1=a2");
    // Direct simulation of the nonlinear plant under complex forcing.
    const double u0=.4,T=2*pi/w,dt=T/400;Cx x1{},x2{};SampledSignal sim;
    const auto field=[&](double t,Cx p,Cx q,Cx& dp,Cx& dq){dp=a1*p+q*q;dq=a2*q+u0*std::exp(Cx(0,w*t));};
    for(int k=0;k<=400*60;++k){const double t=k*dt;sim.times.push_back(t);sim.values.push_back(x1);
        Cx k1,l1,k2,l2,k3,l3,k4,l4;field(t,x1,x2,k1,l1);field(t+.5*dt,x1+.5*dt*k1,x2+.5*dt*l1,k2,l2);field(t+.5*dt,x1+.5*dt*k2,x2+.5*dt*l2,k3,l3);field(t+dt,x1+dt*k3,x2+dt*l3,k4,l4);
        x1+=dt/6*(k1+2.*k2+2.*k3+k4);x2+=dt/6*(l1+2.*l2+2.*l3+l4);}
    const auto est=estimateHarmonics(sim,w,3,.5);nearC(est.at(2)/(u0*u0),H2,1e-5,"Simulated steady-state H2");near(std::abs(est.at(1)),0,1e-6,"Simulated first harmonic of x1");
    const auto bad=harmonicResponse([&]{auto x=m;x.A(1,1)=Cx(0,w);return x;}(),w,o);check(!bad.valid&&!bad.warnings.empty(),"External resonance invalidates H_n");
    // Symmetric Volterra kernels.
    const Cx s1(.2,1),s2(-.1,.7),s3(.3,-.4);
    nearC(volterraKernel(m,{s1,s2})[0],1./((s1-a2)*(s2-a2)*(s1+s2-a1)),1e-14,"H2 closed form");
    nearC(volterraKernel(m,{s1,s2})[0],volterraKernel(m,{s2,s1})[0],1e-15,"H2 symmetric");
    nearC(volterraKernel(m,{Cx(0,w),Cx(0,w)})[0],H2,1e-14,"Diagonal kernel equals the harmonic response");
    auto cubic=m;cubic.N[0](0,1)=.5;cubic.N[0](1,2)=-.25;const auto ch=harmonicResponse(cubic,w,o);
    nearC(volterraKernel(cubic,{Cx(0,w),Cx(0,w),Cx(0,w)})[0],ch.outputs[2][0],1e-13,"Diagonal H3");
    const auto h123=volterraKernel(cubic,{s1,s2,s3})[0];nearC(h123,volterraKernel(cubic,{s3,s1,s2})[0],1e-13,"H3 symmetric");
    const auto chain=[&](Cx a,Cx b,Cx c){const ComplexVector B0{0,1,0};
        const auto v1=ComplexLu(shiftedMatrix(cubic.A,a)).solve(B0);const auto v2=ComplexLu(shiftedMatrix(cubic.A,a+b)).solve(product(cubic.N[0],v1));
        return product(cubic.C,ComplexLu(shiftedMatrix(cubic.A,a+b+c)).solve(product(cubic.N[0],v2)))[0];};
    nearC(h123,(chain(s1,s2,s3)+chain(s1,s3,s2)+chain(s2,s1,s3)+chain(s2,s3,s1)+chain(s3,s1,s2)+chain(s3,s2,s1))/6.,1e-13,"H3 permutation average");
    // Real forcing u=a cos wt: DC and 2w of x1 from the symmetric kernel.
    const double amp=.6;const auto hb=realHarmonicBalance(m,w,amp,4);
    nearC(hb.outputs.at(2)[0],amp*amp/4*volterraKernel(m,{Cx(0,w),Cx(0,w)})[0],1e-13,"Second harmonic under real forcing");
    nearC(hb.outputs.at(0)[0],amp*amp/2*volterraKernel(m,{Cx(0,w),Cx(0,-w)})[0],1e-13,"DC term under real forcing");
    near(std::abs(hb.outputs.at(1)[0]),0,1e-13,"First harmonic of x1 under real forcing");
    // Floquet: z'=(a+u)z with u=b cos(wt) has multiplier exp(a T), for any b.
    BilinearLiftedModel scalar;scalar.A=mat(1,1,{-.2});scalar.B=mat(1,1,{1});scalar.C=mat(1,1,{1});scalar.N.push_back(mat(1,1,{1}));
    const auto fq=bilinearFloquet(scalar,1.5,3,false);near(fq.spectralRadius,std::exp(-.2*2*pi/1.5),1e-9,"Scalar Floquet multiplier");check(fq.stable,"Stable multiplier");
    BilinearLiftedModel unstableLoop;unstableLoop.A=mat(2,2,{-.05,1, -1,-.05});unstableLoop.B=mat(2,1,{0,1});unstableLoop.C=mat(1,2,{1,0});unstableLoop.N.push_back(mat(2,2,{0,0, 1,0}));
    check(!bilinearFloquet(unstableLoop,2.,1.,false).stable,"Mathieu parametric resonance is detected");
    const auto sp=simplePoleCheck([](Cx x){return 3./(x-Cx(0,2));},Cx(0,2),.5);check(sp.simple,"Simple pole check");nearC(sp.residue,3,1e-11,"Simple pole residue");
    check(!simplePoleCheck([](Cx x){return 1./std::pow(x-Cx(0,2),2);},Cx(0,2),.5).simple,"Double pole fails the check");
}
ShellResolventTable annulus(int nodes,std::vector<int> harmonics){
    ShellResolventTable t;for(int k=0;k<nodes;++k){const double e=1+double(k)/(nodes-1);t.energy.push_back(e);t.frequency.push_back(e);t.density.push_back(1);}
    for(int h:harmonics)t.harmonics.emplace(h,ComplexVector(nodes,Cx(1)));
    return t;
}
void shellResolvent(){
    FourierSeries series{1.4,{{1,Cx(.5,.2)},{-1,Cx(.5,-.2)},{3,Cx(.1,0)}}};
    const auto orbit=orbitHarmonicModel(series,.3);const double t=2.2;
    nearC(evaluateTimeResponse(orbit,t),evaluate(series,1.4*t+.3),1e-13,"Orbit pole model evolves the Fourier series");
    // Omega=E on [1,2], rho=1, g_1=1: Y(iw+0)=pi - i log((w-1)/(2-w)).
    const double w=1.37;const auto b=shellBoundaryValue(annulus(41,{1}),w);check(b.valid&&b.roots.size()==1,"One simple root");
    near(b.roots[0].energy,w,1e-10,"Root location");near(b.roots[0].weight,1,1e-10,"Simple-root weight 1/|n Omega'|");
    nearC(b.value,Cx(pi,-std::log((w-1)/(2-w))),1e-9,"Boundary value with principal value");
    // Second harmonic: 2E=3 at E=1.5, weight 1/2, and PV int dE/(3-2E)=0.
    const auto b2=shellBoundaryValue(annulus(41,{2}),3.);nearC(b2.value,Cx(pi/2,0),1e-9,"Boundary value at n=2");near(b2.roots[0].weight,.5,1e-10,"Weight 1/|2 Omega'|");
    // Zero harmonic: an atom at w=0, and -i mass/w elsewhere.
    const auto z=shellBoundaryValue(annulus(11,{0}),.5);near(z.zeroHarmonicMass,1,1e-12,"Zero-harmonic mass");nearC(z.value,Cx(0,-2),1e-11,"Zero harmonic off the atom");
    ShellResolventTable flat=annulus(11,{1});std::fill(flat.frequency.begin(),flat.frequency.end(),1.5);
    check(!shellBoundaryValue(flat,1.5).valid,"Vanishing twist is flagged");
    ShellBoundaryOptions so;so.damping=std::vector<double>(41,2.);so.dampingScale=.01;
    near(*shellBoundaryValue(annulus(41,{1}),w,so).roots[0].stationaryPhaseWeight,shellWeight(.01,2.,1,1.),1e-10,"Stationary-phase weight reuses DampedOrbit");
    // Jordan-block transfer function matches the OrbitGramian impulse response.
    const std::vector<JordanMode> modes{{1,Cx(.6,.1),Cx(.2,-.3)},{-1,Cx(.6,-.1),Cx(.2,.3)},{2,Cx(.1,.2),Cx(-.1,.05)},{-2,Cx(.1,-.2),Cx(-.1,-.05)}};
    const double Om=1.1,twist=-.4;const auto jm=jordanOrbitModel(modes,Om,twist);const auto times=std::vector<double>{0.,.7,3.1};
    const auto reference=jordanImpulseResponse(modes,Om,twist,times);
    for(std::size_t k=0;k<times.size();++k)nearC(evaluateTimeResponse(jm,times[k]),reference[k],1e-13,"Jordan model equals R(t)");
    const auto laurent=cauchyLaurent([&](Cx s){return evaluateResolvent(jm,s);},Cx(0,Om),.4);
    check(laurent.order==2,"Second-order orbit pole");nearC(laurent.laurent[1],-twist*std::norm(modes[0].amplitude),1e-10,"Jordan coefficient");
    check(rejects([&]{jordanOrbitModel({{1,Cx(1),Cx(0)}},1,1);}),"Non-real observable rejected");
    // Damped pairing: closed form for linear Omega, then the d->0 boundary value.
    const auto table=annulus(9,{1});const std::vector<double> rate(9,1.);const double eps=.03;
    const Cx a=Cx(eps,w),exact=I*(std::log(a-Cx(0,2))-std::log(a-Cx(0,1)));
    nearC(dampedShellPairing(table,Cx(0,w),rate,eps),exact,1e-12,"Damped pairing is exact on panels");
    double previous=1e9;for(double d:{1e-2,1e-3,1e-4}){const double e=std::abs(dampedShellPairing(table,Cx(0,w),rate,d)-b.value);check(e<previous,"Damped limit converges");previous=e;}
    check(previous<1e-3,"Damped limit reaches the boundary value");
    // Rational kernel as a resolvent combination.
    for(std::size_t order:{1u,2u,3u,4u}){const RationalCircleKernel kernel(.15,order);const double theta=.8;const auto r=rationalKernelResolvent(kernel,theta);
        for(double alpha:{0.,.75,-2.3})nearC(applyRationalCombination(r,std::exp(Cx(0,alpha))),kernel(theta-alpha),1e-11,"Resolvent combination equals the kernel");}
    // Jordan block: g(U) by explicit 2x2 inverses matches rationalProjectedModel.
    const RationalCircleKernel kernel(.2,2);const double h=.4,theta=.7;const Cx lambda(0,1.2),mu=std::exp(lambda*h);
    PoleResidueModel single;single.poles.push_back({lambda,{Cx(.3),Cx(.5)}});const auto projected=rationalProjectedModel(single,h,theta,kernel);
    const auto rc=rationalKernelResolvent(kernel,theta);Cx entry{};
    for(const auto& term:rc.terms){const Cx d=term.node-mu;entry+=term.coefficient*term.node*mu*h/(d*d);} // (1,2) entry of z(z-U)^{-1}.
    nearC(projected.poles[0].laurent[0],applyRationalCombination(rc,mu)*.3+entry*.5,1e-12,"Projected residue");
    nearC(projected.poles[0].laurent[1],applyRationalCombination(rc,mu)*.5,1e-12,"Projected nilpotent part");
}
void gramians(){
    const ComplexVector lambda{Cx(-.3,2),Cx(-.3,-2),Cx(-1,.5)};const std::vector<ComplexVector> inputs{{1,Cx(.5,.2),Cx(0,.7)}},outputs{{Cx(.4,.1),1,Cx(.3,-.2)}};
    const auto ref=modalGramians(lambda,inputs,outputs);
    ComplexMatrix A(3,3),B(3,1),C(1,3);for(int i=0;i<3;++i){A(i,i)=lambda[i];B(i,0)=inputs[0][i];C(0,i)=outputs[0][i];}
    const auto r=frequencyGramians(A,B,C);check(r.converged,"Gramian quadrature converged");
    for(std::size_t k=0;k<9;++k){nearC(r.controllability.data[k],ref.controllability[k],1e-7,"W_c equals the Cauchy Gramian");nearC(r.observability.data[k],ref.observability[k],1e-7,"W_o equals the Cauchy Gramian");}
    near(r.h2Norm,ref.h2Norm,1e-7,"H2 by Parseval");
    // Lightly damped pole: W=1/(2 a) for 1/(s+a-i w0).
    const double a=.01,w0=5;const auto light=frequencyGramians(mat(1,1,{Cx(-a,w0)}),mat(1,1,{1}),mat(1,1,{1}));near(light.controllability(0,0).real(),1/(2*a),1e-5,"Lightly damped Gramian");
    check(rejects([]{frequencyGramians(mat(1,1,{Cx(.1)}),mat(1,1,{1}),mat(1,1,{1}));}),"Unstable A rejected");
#ifdef DYNSYSLAB_MATH_NUMERICS
    const auto peak=resolventHInfinity(mat(1,1,{Cx(-a,w0)}),mat(1,1,{1}),mat(1,1,{1}));near(peak.value,1/a,1e-5,"H-infinity of a shifted pole");near(peak.frequency,w0,1e-4,"Peak frequency sign is resolved");
#endif
}
void queries(){
    const auto run=[](const std::string& op,const std::string& text){auto out=evaluateResolventQuery(op,core::json::parse(text));if(!out)throw std::runtime_error("Unclaimed "+op);return *out;};
    check(!evaluateResolventQuery("unknown",core::json::parse("{}")).has_value(),"Unknown op is not claimed");
    const auto lt=run("laplace_transform",R"({"times":[0,0.5,1,1.5,2],"values":[1,0.6065306597,0.3678794412,0.2231301601,0.1353352832],"s":[[1,0]],"z":[[1.2,0]],"sample_interval":0.5})");
    near(lt.at("results").array()[0].at("value").array()[0].number(),.5*(1-std::exp(-4.)),2e-3,"Query Laplace");check(lt.at("z_transform").array().size()==1,"Query z-transform");
    const auto il=run("inverse_laplace",R"({"poles":[{"eigenvalue":[-1,0],"laurent":[1]}],"times":[1],"method":"talbot"})");near(il.at("maximum_error_against_model").number(),0,1e-11,"Query inverse");check(il.at("scope_verified").boolean(),"Query inverse scope");
    const auto kr=run("koopman_resolvent",R"({"poles":[{"eigenvalue":[0,1],"laurent":[1,2]}],"s":[[1,1]],"lattice":{"kind":"equilibrium","exponents":[[-1,0],[-2.5,0]],"maximum_order":2}})");
    near(kr.at("results").array()[0].at("value").array()[0].number(),3,1e-14,"Query resolvent");check(kr.at("lattice").at("poles").array().size()==6,"Query lattice");
    const auto re=run("resolvent_residue",R"({"poles":[{"eigenvalue":[0,0],"laurent":[2,3]}],"center":[0,0],"radius":0.5})");near(re.at("laurent").array()[1].array()[0].number(),3,1e-10,"Query residue");check(re.at("order").number()==2,"Query order");
    core::json::Value::Array times,values;for(int k=0;k<=2000;++k){times.emplace_back(k*.05);values.emplace_back(1.+std::exp(-.2*k*.05));}
    core::json::Value::Object gla{{"times",times},{"values",values},{"eigenvalues",core::json::parse("[[0,0],[-0.2,0]]")}};
    const auto gl=*evaluateResolventQuery("generalized_laplace",core::json::Value(gla));near(gl.at("modes").array()[0].at("coefficient").array()[0].number(),1,1e-3,"Query GLA");
    const auto ps=run("pseudospectrum",R"({"matrix":[[-1,0],[0,-2]],"real_axis":[-1,0],"imag_axis":[0]})");near(ps.at("sigma_min").array()[1].number(),1,1e-12,"Query pseudospectrum");
    const std::string lift=R"("A":[[-1.3,0,1],[0,-0.7,0],[0,0,-1.4]],"B":[[0],[1],[0]],"C":[[1,0,0]],"N":[[[0,0,0],[0,0,0],[0,2,0]]])";
    const auto hr=run("harmonic_response","{"+lift+R"(,"frequency":1.9,"maximum_order":3,"amplitude":0.5,"real_input":true})");
    const Cx H2=1./((Cx(0,3.8)+1.3)*std::pow(Cx(0,1.9)+.7,2));near(hr.at("harmonics").array()[1].array()[0].array()[0].number(),H2.real(),1e-13,"Query H2");check(hr.at("floquet").at("stable").boolean(),"Query Floquet");
    const auto vk=run("volterra_kernels","{"+lift+R"(,"arguments":[[[0,1.9],[0,1.9]]]})");near(vk.at("kernels").array()[0].at("kernel").array()[0].array()[1].number(),H2.imag(),1e-13,"Query Volterra");
    const auto sb=run("shell_resolvent",R"({"mode":"boundary","energies":[1,1.25,1.5,1.75,2],"frequencies":[1,1.25,1.5,1.75,2],"densities":[1,1,1,1,1],"harmonics":[{"harmonic":1,"values":[1,1,1,1,1]}],"frequency":1.5})");
    near(sb.at("absorptive").number(),pi,1e-9,"Query shell boundary");
    const auto sj=run("shell_resolvent",R"({"mode":"jordan","modes":[{"harmonic":1,"amplitude":[1,0],"slope":[0,0]},{"harmonic":-1,"amplitude":[1,0],"slope":[0,0]}],"frequency":1,"action_twist":0.5,"times":[0]})");
    near(std::abs(sj.at("response").array()[0].at("value").array()[0].number()),0,1e-14,"Query Jordan");
    const auto sp=run("shell_resolvent",R"({"mode":"projector","bandwidth":0.2,"order":2,"angle":0.5,"sample_interval":0.5,"poles":[{"eigenvalue":[0,1],"laurent":[1]}]})");check(sp.at("terms").array().size()==4,"Query projector");
    const auto fg=run("frequency_gramian",R"({"A":[[-2]],"B":[[1]],"C":[[1]]})");near(fg.at("controllability").array()[0].array()[0].array()[0].number(),.25,1e-8,"Query Gramian");near(fg.at("h2_norm").number(),.5,1e-8,"Query H2");
    check(rejects([&]{run("laplace_transform",R"({"times":[0,1],"values":[1,2]})");}),"Missing s rejected");
    check(rejects([&]{run("shell_resolvent",R"({"mode":"nope"})");}),"Bad shell mode rejected");
    check(rejects([&]{run("koopman_resolvent",R"({"s":[1]})");}),"Missing model rejected");
}
}
int main(){
    try{
        linearAlgebra();laplaceIdentities();inverseIdentities();poleModels();lattices();pseudospectrum();
        generalizedLaplace();frequencyResponse();shellResolvent();gramians();queries();
        std::cout<<"resolvent tools: "<<checks<<" checks passed\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAILED after "<<checks<<" checks: "<<e.what()<<"\n";return 1;}
}
