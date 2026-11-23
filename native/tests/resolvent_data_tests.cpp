// Mathematical identities for the data-driven Koopman resolvent tools. No study data.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/ContractionCheck.hpp>
#include <dynsyslab/math/DiscreteResolvent.hpp>
#include <dynsyslab/math/ForcedResponseEstimation.hpp>
#include <dynsyslab/math/HarmonicContinuation.hpp>
#include <dynsyslab/math/IsochronDecomposition.hpp>
#include <dynsyslab/math/KoopmanFrequencyResponse.hpp>
#include <dynsyslab/math/LaplaceTransform.hpp>
#include <dynsyslab/math/PronyIdentification.hpp>
#include <dynsyslab/math/ResolventDataQuery.hpp>
#include <dynsyslab/math/ResolventLeastSquares.hpp>
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
int checks=0;
void check(bool ok,const std::string& m){++checks;if(!ok)throw std::runtime_error(m);}
std::string sci(double x){char b[64];std::snprintf(b,64,"%.3e",x);return b;}
void near(double a,double b,double tol,const std::string& m){++checks;if(!(std::abs(a-b)<=tol))throw std::runtime_error(m+": got "+sci(a)+" expected "+sci(b));}
void nearC(Cx a,Cx b,double tol,const std::string& m){++checks;if(!(std::abs(a-b)<=tol))throw std::runtime_error(m+": error "+sci(std::abs(a-b)));}
template<class F>bool rejects(F&& f){try{f();}catch(const std::invalid_argument&){return true;}catch(const std::domain_error&){return true;}return false;}
SampledSignal sampled(double T,std::size_t n,const std::function<Cx(double)>& f){
    SampledSignal g;for(std::size_t k=0;k<n;++k){const double t=T*double(k)/double(n-1);g.times.push_back(t);g.values.push_back(f(t));}return g;
}
void leastSquaresIdentities(){
    ComplexMatrix A(4,2);A.data={1,Cx(0,1), 2,0, Cx(1,1),3, 0,-1};
    const auto s=thinSvd(A);ComplexMatrix S(2,2);S(0,0)=s.singularValues[0];S(1,1)=s.singularValues[1];
    const auto back=product(product(s.U,S),adjoint(s.V));for(std::size_t k=0;k<A.data.size();++k)nearC(back.data[k],A.data[k],1e-13,"SVD reconstruction");
    const auto wide=thinSvd(adjoint(A));near(wide.singularValues[0],s.singularValues[0],1e-13,"SVD of the adjoint");
    ComplexMatrix x(2,1);x(0,0)=Cx(1,-2);x(1,0)=.5;const auto b=product(A,x);const auto r=leastSquares(A,b);check(r.rank==2,"Full rank");
    nearC(r.x(0,0),x(0,0),1e-13,"Least squares recovers x");
}
void discreteIdentities(){
    PoleResidueModel c;c.poles.push_back({Cx(-.3,1.1),{Cx(1,.5),Cx(-.4,.2),Cx(.1,.3)}});c.poles.push_back({Cx(-.8,0),{Cx(.7)}});
    const double h=.25;const auto d=discretizePoleModel(c,h);
    for(std::size_t k:{0u,1u,2u,3u,7u,20u})nearC(evaluateDiscreteSequence(d,k),evaluateTimeResponse(c,double(k)*h),1e-12,"Sampling a Jordan model");
    // X(z) equals the z-transform of the sequence, and uses the same s<->z map.
    ComplexVector y;for(std::size_t k=0;k<600;++k)y.push_back(evaluateDiscreteSequence(d,k));
    const Cx z(1.15,.3);nearC(evaluateDiscreteResolvent(d,z),zTransform(y,z),1e-11,"X(z) closed form");
    nearC(d.poles[0].multiplier,zFromLaplace(c.poles[0].eigenvalue,h),1e-15,"Multiplier is exp(lambda h)");
    // Round trip through the principal branch.
    const auto back=continuousPoleModel(d);
    for(std::size_t j=0;j<c.poles.size();++j){nearC(back.poles[j].eigenvalue,c.poles[j].eigenvalue,1e-13,"Round-trip eigenvalue");
        for(std::size_t q=0;q<c.poles[j].laurent.size();++q)nearC(back.poles[j].laurent[q],c.poles[j].laurent[q],1e-11,"Round-trip Laurent");}
    // Residue of X(z)/z at mu by a Cauchy contour.
    const auto e=cauchyLaurent([&](Cx w){return evaluateDiscreteResolvent(d,w)/w;},d.poles[1].multiplier,.05);nearC(e.laurent[0],Cx(.7),1e-9,"Discrete residue by contour");
    // Aliasing: lambda and lambda+2 pi i/h share a multiplier.
    PoleResidueModel a;a.poles.push_back({Cx(-.1,1),{1}});a.poles.push_back({Cx(-.1,1+2*pi/h),{1}});a.poles.push_back({Cx(0,pi/h),{1}});
    const auto r=aliasReport(a,h);check(!r.entries[0].aliased&&r.entries[1].aliased,"Aliased pole flagged");check(r.entries[1].turns==1,"Alias turns");
    nearC(r.entries[1].principal,a.poles[0].eigenvalue,1e-12,"Principal alias");check(r.entries[2].onNyquist,"Nyquist flagged");check(r.collisions.size()==1,"Alias collision");
    std::vector<int> branch{0,1};PoleResidueModel two;two.poles={a.poles[0],a.poles[1]};
    nearC(continuousPoleModel(discretizePoleModel(two,h),branch).poles[1].eigenvalue,a.poles[1].eigenvalue,1e-12,"Declared branch restores the pole");
    check(rejects([&]{discretizePoleModel(c,0);}),"Zero interval rejected");
}
void pronyIdentities(){
    const double h=.1;const std::vector<Cx> lambda{Cx(-.2,1.3),Cx(-.2,-1.3),Cx(-.7,4),Cx(0,0)},res{Cx(.5,.2),Cx(.5,-.2),Cx(.3,.1),Cx(1.2)};
    ComplexVector y;for(int k=0;k<200;++k){Cx v{};for(std::size_t j=0;j<4;++j)v+=res[j]*std::exp(lambda[j]*(k*h));y.push_back(v);}
    PronyOptions o;o.sampleInterval=h;const auto r=pronyIdentify(y,o);check(r.order==4&&r.orderSource=="rank","Rank selects the order");near(r.relativeResidual,0,1e-9,"Exact fit");
    for(std::size_t j=0;j<4;++j){const auto it=std::min_element(r.poles.begin(),r.poles.end(),[&](const PronyPole& a,const PronyPole& b){return std::abs(a.eigenvalue-lambda[j])<std::abs(b.eigenvalue-lambda[j]);});
        nearC(it->eigenvalue,lambda[j],1e-7,"Pole recovered");nearC(it->residue,res[j],1e-7,"Residue recovered");}
    check(!r.clustered&&!r.nearContinuous,"Clean spectrum is not flagged");check(r.singularGap>1e3,"Clear order gap");
    const auto scan=std::find_if(r.orderScan.begin(),r.orderScan.end(),[](const PronyOrderScan& s){return s.order==3;});check(scan!=r.orderScan.end()&&scan->relativeResidual>1e-3,"Underfit order has a visible residual");
    nearC(evaluateDiscreteSequence(r.model,17),y[17],1e-8,"Model reproduces the signal");
    // Two poles closer than the Rayleigh limit of a short record.
    ComplexVector close;for(int k=0;k<40;++k)close.push_back(std::exp(Cx(0,.5*k))+std::exp(Cx(0,.52*k)));
    PronyOptions c;c.order=2;check(pronyIdentify(close,c).clustered,"Clustered poles flagged");
    // Dense band: a sum of many equal-weight tones has no low-rank structure at this order.
    ComplexVector band;for(int k=0;k<120;++k){Cx v{};for(int m=0;m<400;++m)v+=std::exp(Cx(0,(.5+.5*m/399.)*k+m*m*.37))/20.;band.push_back(v);}
    PronyOptions b;b.order=4;b.pencilParameter=40;const auto br=pronyIdentify(band,b);check(br.nearContinuous&&br.relativeResidual>.05,"Near-continuous spectrum flagged");
    check(rejects([]{pronyIdentify(ComplexVector(10,Cx(0)));}),"Zero signal rejected");
}
struct Plant {std::vector<ComplexVector> states;ComplexVector inputs;SampledSignal x1;};
// x1'=a1 x1+x2^2, x2'=a2 x2+u under complex forcing (arXiv:2603.05771, eq. 9).
Plant simulate(double a1,double a2,double w,double u0,double h,std::size_t n,Cx p,Cx q){
    Plant out;const auto f=[&](double t,Cx x,Cx y,Cx& dx,Cx& dy){dx=a1*x+y*y;dy=a2*y+u0*std::exp(Cx(0,w*t));};const double dt=h/8;
    for(std::size_t k=0;k<n;++k){const double t=double(k)*h;out.states.push_back({p,q});out.inputs.push_back(u0*std::exp(Cx(0,w*t)));out.x1.times.push_back(t);out.x1.values.push_back(p);
        for(int s=0;s<8;++s){const double r=t+s*dt;Cx k1,l1,k2,l2,k3,l3,k4,l4;f(r,p,q,k1,l1);f(r+.5*dt,p+.5*dt*k1,q+.5*dt*l1,k2,l2);f(r+.5*dt,p+.5*dt*k2,q+.5*dt*l2,k3,l3);f(r+dt,p+dt*k3,q+dt*l3,k4,l4);
            p+=dt/6*(k1+2.*k2+2.*k3+k4);q+=dt/6*(l1+2.*l2+2.*l3+l4);}}
    return out;
}
void forcedEstimation(){
    const double a1=-1.3,a2=-.7,w=1.9,u0=.4;const Cx H2=1./((Cx(0,2*w)-a1)*std::pow(Cx(0,w)-a2,2));
    const auto plant=simulate(a1,a2,w,u0,2*pi/w/64,64*40,Cx(.6,.1),Cx(-.5,.3));
    const auto e=estimateForcedHarmonics(plant.x1,w,u0,3,true,.5);nearC(e.response.at(2),H2,1e-7,"Projection estimate of H2");near(std::abs(e.response.at(1)),0,1e-7,"H1 of x1 vanishes");check(e.tailRelativeResidual<1e-6,"Settled tail");
    check(estimateForcedHarmonics(plant.x1,w,u0,3,false,.5).response.empty(),"Real input does not report H_n");
    // EDMD on the exact six-function lift (x1,x2,u,x2^2,x2 u,u^2).
    SkewDmdInput in;const std::size_t count=260;const double h=2*pi/w/64;const auto transient=simulate(a1,a2,w,u0,h,count,Cx(.6,.1),Cx(-.5,.3));
    in.states=transient.states;in.inputs=transient.inputs;in.sampleInterval=h;in.frequency=w;in.maximumHarmonic=2;
    in.dictionary={{1,0,0},{0,1,0},{0,0,1},{0,2,0},{0,1,1},{0,0,2}};
    const auto d=skewProductDmd(in);check(d.feasible,"Exact lift is feasible");near(d.closureResidual,0,1e-8,"Closed dictionary");
    nearC(*d.harmonics[1].response,H2,1e-5,"DMD estimate of H2");near(std::abs(*d.harmonics[0].response),0,1e-5,"DMD H1 of x1");check(d.harmonics[1].modeResidual<1e-6,"u^2 is an eigenfunction of K");
    // Steady-state-only data cannot identify the lift.
    SkewDmdInput steady=in;const auto late=simulate(a1,a2,w,u0,h,64*40,Cx(.6,.1),Cx(-.5,.3));
    steady.states.assign(late.states.end()-200,late.states.end());steady.inputs.assign(late.inputs.end()-200,late.inputs.end());steady.rankTolerance=1e-9;
    const auto s=skewProductDmd(steady);check(!s.feasible&&!s.warnings.empty(),"Steady data are reported as infeasible");
    SkewDmdInput open=in;open.dictionary={{1,0,0},{0,1,0},{0,0,1},{0,0,2}};const auto o=skewProductDmd(open);check(o.closureResidual>1e-4,"Unclosed dictionary has a residual");
}
void subharmonics(){
    const double w=2,u0=.5;const Cx A(.4,.3);const auto y=[&](double t){return A*std::exp(Cx(0,w*t/2))+Cx(.2,-.1)*std::exp(Cx(0,w*t))+.3;};
    const auto a=subharmonicResponse(sampled(200*pi,32001,y),w,u0,2,4);check(a.detected&&a.coverLabelled,"Period doubling detected with cover labels");
    const Cx H=A/std::sqrt(u0);near(*a.magnitude,std::abs(H),1e-6,"|H_{1/2}|");nearC(*a.invariantPower,H*H,1e-6,"H^n is invariant");check(a.branches.size()==2,"Two branches");
    check(std::any_of(a.branches.begin(),a.branches.end(),[&](Cx b){return std::abs(b-H)<1e-6;})&&std::any_of(a.branches.begin(),a.branches.end(),[&](Cx b){return std::abs(b+H)<1e-6;}),"Branches are +-H");
    const auto label=std::find_if(a.labels.begin(),a.labels.end(),[](const SubharmonicLabel& l){return l.index==1;});check(label!=a.labels.end()&&label->frequencyRatio=="1/2","PhaseCoverFourier ratio");
    // Advancing the record by one forcing period rotates H by e^{2 pi i/n} and fixes H^n.
    const double T=2*pi/w;const auto b=subharmonicResponse(sampled(200*pi,32001,[&](double t){return y(t+T);}),w,u0,2,4);
    nearC(b.coefficients.at(1),-a.coefficients.at(1),1e-6,"Shift rotates the coefficient");nearC(*b.invariantPower,*a.invariantPower,1e-6,"Shift keeps H^n");near(*b.phaseModulo,*a.phaseModulo,1e-6,"Phase modulo 2 pi/n");
    const auto third=subharmonicResponse(sampled(240*pi,32001,[&](double t){return std::exp(Cx(0,w*t/3));}),w,1.,3,6);nearC(*third.invariantPower,Cx(1),1e-6,"Third subharmonic");
    const auto six=subharmonicResponse(sampled(240*pi,32001,[&](double t){return std::exp(Cx(0,w*t/6));}),w,1.,6,6);check(six.detected&&!six.coverLabelled,"Non-prime-power cover has no labels");
    const auto none=subharmonicResponse(sampled(200*pi,32001,[&](double t){return std::exp(Cx(0,w*t));}),w,u0,2,4);check(!none.detected&&!none.magnitude&&none.branches.empty(),"No subharmonic, no H_{1/n}");
}
void isochrons(){
    const double W=1.4,nu=-.3,theta=.9;FourierSeries cycle{W,{{1,Cx(.8,.2)},{-1,Cx(.8,-.2)},{3,Cx(.1,.05)},{-3,Cx(.1,-.05)}}};
    const auto y=[&](double t){return evaluate(cycle,theta+W*t)+Cx(.6,.2)*std::exp(Cx(nu,W)*t);};const auto g=sampled(120,24001,y);
    IsochronSplitOptions o;o.asymptoticPhase=theta;const auto r=isochronSplit(g,cycle,o);
    check(r.phaseSource=="phase_map","Phase from the map");near(*r.phaseMismatch,0,1e-5,"Tail agrees with the map");near(*r.decayRate,nu,1e-6,"Decay rate equals the exponent");
    nearC(r.stationary[500],evaluate(cycle,theta+W*g.times[500]),1e-13,"Stationary part");nearC(r.transient[500],y(g.times[500])-r.stationary[500],1e-13,"Transient part");
    const auto wrong=[&]{auto x=o;x.asymptoticPhase=theta+.2;return isochronSplit(g,cycle,x);}();check(*wrong.phaseMismatch>.1&&wrong.warnings.size()>=1,"Wrong isochron is detected");
    const auto fitted=isochronSplit(g,cycle);check(fitted.phaseSource=="tail_estimate",r.phaseSource);near(fitted.asymptoticPhase,theta,1e-5,"Tail phase");
    // Radial isochrons x=r cos phi, y=r sin phi: the phase map is atan2.
    const auto mapped=isochronSplitFromMap(g,cycle,[](const std::vector<double>& x){return std::atan2(x[1],x[0]);},{std::cos(theta)*2,std::sin(theta)*2});near(mapped.asymptoticPhase,theta,1e-14,"Phase map evaluated at x0");
    FourierSeries even{W,{{2,Cx(1)},{-2,Cx(1)}}};const auto sym=isochronSplit(sampled(60,6001,[&](double t){return evaluate(even,.3+W*t);}),even);
    check(std::any_of(sym.warnings.begin(),sym.warnings.end(),[](const std::string& s){return s.find("modulo")!=std::string::npos;}),"Symmetric observable fixes phase only modulo");
    check(rejects([&]{isochronSplit(g,FourierSeries{W,{{0,Cx(1)}}});}),"A constant cycle has no phase");
}
ContinuationProblem duffing(double d,double beta,double a){
    ContinuationProblem p;p.dimension=2;p.amplitude=a;p.frequency=.6;p.lower=.6;p.upper=2.2;p.harmonics=1;p.step=.02;p.maximumStep=.04;
    p.field=[d,beta](const double* x,double u,double,double* f){f[0]=x[1];f[1]=-d*x[1]-x[0]-beta*x[0]*x[0]*x[0]+u;};
    return p;
}
void continuation(){
    // Linear oscillator: the exact frequency response at every point, no folds.
    auto linear=duffing(.2,0,.5);linear.harmonics=2;const auto L=continueHarmonicBalance(linear);check(L.termination=="parameter_bound"&&!L.branchLost,"Linear branch reaches the bound");check(L.folds.empty(),"Linear branch has no fold");
    for(const auto& q:L.points){const Cx H=1./(1.-q.value*q.value+Cx(0,.2*q.value));nearC(q.response.at(1),.25*H,1e-6,"Linear response c_1=aH/2");check(q.stable.value_or(false),"Linear branch is stable");}
    // One-harmonic Duffing: G(A,w)=[(1-w^2)A+(3/4)beta A^3]^2+(d w A)^2-a^2=0.
    const double d=.1,beta=.4,a=.3;const auto D=continueHarmonicBalance(duffing(d,beta,a));check(!D.branchLost,"Duffing branch is followed through both folds");check(D.folds.size()==2,"Two folds detected");
    const auto G=[&](double A,double w){const double r=(1-w*w)*A+.75*beta*A*A*A;return r*r+d*d*w*w*A*A-a*a;};
    for(const auto& q:D.points)near(G(2*std::abs(q.response.at(1)),q.value),0,1e-6,"Harmonic-balance amplitude equation");
    for(const auto& f:D.folds){const double A=2*f.responseMagnitude,e=1e-5;near((G(A+e,f.value)-G(A-e,f.value))/(2*e),0,.08,"G_A vanishes at a fold");}
    check(D.folds[0].value>D.folds[1].value,"Upper fold is reached first on a hardening curve");
    bool unstable=false;for(std::size_t k=D.folds[0].afterPoint+2;k+1<D.folds[1].afterPoint;++k)unstable=unstable||!D.points[k].stable.value_or(true);check(unstable,"Middle branch is Floquet unstable");
    // A field that stops existing must end the run, not continue it.
    ContinuationProblem lost;lost.dimension=1;lost.frequency=1;lost.amplitude=.1;lost.parameter=.2;lost.variable=ContinuationVariable::Parameter;lost.lower=0;lost.upper=2;lost.harmonics=1;lost.computeStability=false;
    lost.field=[](const double* x,double u,double p,double* f){f[0]=p>.8?std::nan(""):-x[0]+u+p;};const auto R=continueHarmonicBalance(lost);check(R.branchLost&&R.termination=="non_finite_field","Loss of the field is reported");check(R.points.back().value<=.8,"No point beyond the loss");
    auto bad=duffing(d,beta,a);bad.lower=1;check(rejects([&]{continueHarmonicBalance(bad);}),"Start outside the bounds rejected");
}
void contraction(){
    ComplexMatrix J(2,2);J.data={-2,1, 0,-3};const auto n=logarithmicNorms(J);near(n.infinity,-1,1e-14,"mu_inf");near(n.one,-2,1e-14,"mu_1");near(n.two,(-5+std::sqrt(2.))/2,1e-12,"mu_2");
    // A diagonal metric can make a nonnormal matrix contractive in mu_inf.
    ComplexMatrix K(2,2);K.data={-1,10, 0,-1};check(logarithmicNorms(K).infinity>0,"Identity metric fails");check(logarithmicNorms(K,{1,20}).infinity<0,"Weighted metric certifies");
    // x'=-x-x^3+u: dF/dx=-1-3x^2 <= -1, and [-2,2] is invariant for |u|<=1.
    SampledContractionInput in;in.stateDimension=1;in.lower={-2,-1};in.upper={2,1};in.samplesPerAxis=9;
    in.jacobian=[](const std::vector<double>& x){ComplexMatrix m(1,1);m(0,0)=-1-3*x[0]*x[0];return m;};
    in.field=[](const std::vector<double>& x){return std::vector<double>{-x[0]-x[0]*x[0]*x[0]+x[1]};};
    const auto r=sampledContraction(in);check(r.label=="sampled"&&r.contractiveAtSamples,"Sampled contraction");near(r.rateEstimate,-1,1e-12,"Worst sampled rate");check(*r.invariantAtSamples,"Box is invariant");
    check(std::any_of(r.warnings.begin(),r.warnings.end(),[](const std::string& s){return s.find("SAMPLED")==0;}),"Sampled label in the warnings");
    auto outward=in;outward.field=[](const std::vector<double>& x){return std::vector<double>{x[0]};};const auto o=sampledContraction(outward);check(!*o.invariantAtSamples&&o.invarianceViolations>0,"Outward field fails invariance");
    auto bistable=in;bistable.jacobian=[](const std::vector<double>& x){ComplexMatrix m(1,1);m(0,0)=1-3*x[0]*x[0];return m;};check(!sampledContraction(bistable).contractiveAtSamples,"x-x^3 is not contracting near 0");
    // A sampled grid can miss a narrow expansion, which is why the label matters.
    auto spike=in;spike.samplesPerAxis=3;spike.jacobian=[](const std::vector<double>& x){ComplexMatrix m(1,1);m(0,0)=std::abs(x[0]-.7)<.05?1.:-1.;return m;};check(sampledContraction(spike).contractiveAtSamples,"Coarse samples miss a spike");
    const auto iv=intervalContraction({-3,-1,0,-4},{-2,1,.5,-3},2);check(iv.label=="caller_interval"&&iv.contractiveForInterval,"Interval bound");near(iv.upperBound.infinity,-1,1e-14,"Interval mu_inf");
    // The interval bound dominates the log norm of every vertex.
    for(double a:{-3.,-2.})for(double b:{-1.,1.})for(double c:{0.,.5})for(double e:{-4.,-3.}){ComplexMatrix v(2,2);v.data={a,b,c,e};const auto m=logarithmicNorms(v);
        check(m.two<=iv.upperBound.two+1e-12&&m.one<=iv.upperBound.one+1e-12&&m.infinity<=iv.upperBound.infinity+1e-12,"Interval bound encloses vertices");}
    check(!intervalContraction({-1,0,0,-1},{1,0,0,-1},2).contractiveForInterval,"Indefinite interval is not certified");
    check(rejects([]{intervalContraction({1,0,0,0},{0,0,0,0},2);}),"Reversed interval rejected");
}
void queries(){
    const auto run=[](const std::string& op,const std::string& text){auto out=evaluateResolventDataQuery(op,core::json::parse(text));if(!out)throw std::runtime_error("Unclaimed "+op);return *out;};
    check(!evaluateResolventDataQuery("laplace_transform",core::json::parse("{}")).has_value(),"Foreign op is not claimed");
    const auto dr=run("discrete_resolvent",R"({"poles":[{"eigenvalue":[-0.5,1],"laurent":[1,0.5]}],"sample_interval":0.5,"z":[[1.2,0]],"indices":[0,3]})");
    nearC(Cx(dr.at("sequence").array()[0].at("value").array()[0].number(),0),1,1e-14,"Query sequence");check(!dr.at("aliasing").at("entries").array()[0].at("aliased").boolean(),"Query alias");
    const auto dc=run("discrete_resolvent",R"({"discrete_poles":[{"multiplier":[0.5,0],"laurent":[2]}],"sample_interval":1})");near(dc.at("continuous_poles").array()[0].at("eigenvalue").array()[0].number(),std::log(.5),1e-14,"Query continuous pole");
    const auto pr=run("prony_identification",R"({"values":[1,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125],"sample_interval":1})");check(pr.at("order").number()==1,"Query Prony order");near(pr.at("poles").array()[0].at("multiplier").array()[0].number(),.5,1e-10,"Query Prony pole");
    core::json::Value::Array t,v;for(int k=0;k<=1600;++k){const double x=k*.01;t.emplace_back(x);v.emplace_back(core::json::Value::Array{.7*std::cos(2*x),.7*std::sin(2*x)});}
    core::json::Value::Object f{{"times",t},{"values",v},{"frequency",2.},{"amplitude",.5},{"maximum_harmonic",2.}};const auto fh=*evaluateResolventDataQuery("forced_harmonic_estimate",core::json::Value(f));
    near(fh.at("response").array()[0].array()[1].number(),1.4,1e-5,"Query forced H1");
    core::json::Value::Object sh{{"times",t},{"values",v},{"frequency",4.},{"amplitude",1.},{"cover",2.}};const auto sr=*evaluateResolventDataQuery("subharmonic_response",core::json::Value(sh));check(sr.at("detected").boolean(),"Query subharmonic");near(sr.at("magnitude").number(),.7,1e-5,"Query |H|");
    core::json::Value::Object is{{"times",t},{"values",v},{"cycle",core::json::parse(R"({"angular_frequency":2,"coefficients":[[1,0.7,0]]})")},{"asymptotic_phase",0.}};const auto ir=*evaluateResolventDataQuery("isochron_split",core::json::Value(is));near(ir.at("phase_mismatch").number(),0,1e-5,"Query isochron");
    const auto hc=run("harmonic_continuation",R"({"field":[[{"coefficient":1,"powers":[0,1,0,0]}],[{"coefficient":-0.2,"powers":[0,1,0,0]},{"coefficient":-1,"powers":[1,0,0,0]},{"coefficient":1,"powers":[0,0,1,0]}]],"frequency":0.5,"amplitude":0.5,"variable":"frequency","lower":0.5,"upper":0.8,"harmonics":1,"maximum_step":0.05})");
    check(hc.at("termination").string()=="parameter_bound","Query continuation");check(hc.at("folds").array().empty(),"Query continuation folds");
    const auto cc=run("contraction_check",R"({"field":[[{"coefficient":-1,"powers":[1,0,0]},{"coefficient":-1,"powers":[3,0,0]},{"coefficient":1,"powers":[0,1,0]}]],"lower":[-2,-1],"upper":[2,1],"samples_per_axis":9})");
    check(cc.at("label").string()=="sampled"&&cc.at("contractive_at_samples").boolean()&&cc.at("invariant_at_samples").boolean(),"Query sampled contraction");near(cc.at("rate_estimate").number(),-1,1e-12,"Query rate");
    const auto ci=run("contraction_check",R"({"interval_lower":[[-3,-1],[0,-4]],"interval_upper":[[-2,1],[0.5,-3]]})");check(ci.at("label").string()=="caller_interval"&&ci.at("contractive_for_interval").boolean(),"Query interval");
    check(rejects([&]{run("harmonic_continuation",R"({"field":[[]],"frequency":1,"amplitude":1,"variable":"bad","lower":0.5,"upper":2})");}),"Bad variable rejected");
    check(rejects([&]{run("discrete_resolvent",R"({})");}),"Empty discrete query rejected");check(rejects([&]{run("contraction_check",R"({"field":[[{"coefficient":1,"powers":[1]}]],"lower":[0,0],"upper":[1,1]})");}),"Bad powers rejected");
}
}
int main(){
    try{
        leastSquaresIdentities();discreteIdentities();pronyIdentities();forcedEstimation();subharmonics();isochrons();continuation();contraction();queries();
        std::cout<<"resolvent data tools: "<<checks<<" checks passed\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAILED after "<<checks<<" checks: "<<e.what()<<"\n";return 1;}
}
