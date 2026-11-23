#include <dynsyslab/math/Arithmetic.hpp>
#include <dynsyslab/math/Distributions.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/math/Query.hpp>
#include <dynsyslab/math/Resonance.hpp>
#include <dynsyslab/math/Topology.hpp>
#include <dynsyslab/math/PhaseCoverFourier.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::math;
using C=std::complex<double>;
namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void near(C a,C b,double tol=2e-11){check(std::isfinite(std::abs(a))&&std::abs(a-b)<=tol*(1+std::abs(b)),"Numeric identity failed");}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected rejection");}
double exactNumber(const std::string& s){auto i=s.find('/');return i==std::string::npos?std::stod(s):std::stod(s.substr(0,i))/std::stod(s.substr(i+1));}
void arithmetic(){
    auto resonance=phaseCoverResonance(2,"12/9");
    check(resonance.baseToDriveRatio=="4/3"&&resonance.characterRatio=="3/4"&&
          resonance.driveHarmonic=="-1"&&resonance.minimumDepth==2&&resonance.branchClasses=="4",
          "Exact resonance kernel and branch resolution");
    auto blind=phaseCoverResonance(2,"3/2");
    check(blind.minimumDepth==0&&blind.characterRatio=="2"&&blind.driveHarmonic=="-3","Fiber-blind resonance");
    rejects([]{phaseCoverResonance(2,"0");});
    rejects([]{phaseCoverResonance(4,"1");});
    check(canonicalRational("-12/-18")=="2/3","Normalization");
    check(!primeValuation("0",3),"Infinite zero valuation");
    check(primeValuation("27/10",3)==3,"Numerator valuation");
    check(primeValuation("7/125",5)==-3,"Denominator valuation");
    for(unsigned p:{2,3,5,7})for(int a=-7;a<=7;++a)for(int b=-7;b<=7;++b)for(int c=-4;c<=4;++c){
        auto x=std::to_string(a)+"/10",y=std::to_string(b)+"/9",z=std::to_string(c)+"/7";
        auto d=exactNumber(padicDistance(x,z,p)),u=exactNumber(padicDistance(x,y,p)),v=exactNumber(padicDistance(y,z,p));
        check(d<=std::max(u,v)*(1+1e-14),"Ultrametric triangle");
    }
    auto e=padicExpansion("-1/2",3,6);
    check(e.valuation==0,"Unit valuation");for(auto digit:e.digits)check(digit==1,"3-adic -1/2 digits");
    auto negative=padicExpansion("1/6",3,4);check(negative.valuation==-1,"Laurent valuation");
    check(compatiblePhaseResidues("5","14",3,2),"Compatible inverse-limit phase");
    check(!compatiblePhaseResidues("5","13",3,2),"Incompatible inverse-limit phase");
    check(continuedFraction("355/113").terms==std::vector<std::string>({"3","7","16"}),"Continued fraction");
    check(continuedFraction("355/113").convergents.back()=="355/113","Final convergent");
    check(sternBrocotPath("5/3").size()==3,"Stern Brocot path");
    check(fordCirclesTangent("1/2","2/3")&&!fordCirclesTangent("1/2","3/4"),"Ford determinant");
    check(fordCircle("2/3").radius=="1/18","Ford radius");
    check(rationalMediant("1/2","2/3")=="3/5","Mediant");
    check(multiplicativeDistanceProduct("2/3","5/7")=="210","Prime word distance");
    rejects([]{primeValuation("1",4);});rejects([]{canonicalRational("1/0");});
    rejects([]{compatiblePhaseResidues("9","14",3,2);});rejects([]{padicExpansion("1",3,0);});
    auto clock=advancePhaseCover(2,{1,1},"5/6","1/3");
    check(clock.basePhaseTurns=="1/6"&&clock.returnCount=="1"&&clock.prefixCarry=="1","suspension carry past displayed prefix");
    check(clock.digits==std::vector<std::uint32_t>({0,0})&&clock.phaseTurns.back()=="1/24","post-gluing phase continuity");
    auto backwards=advancePhaseCover(3,{0,0},"1/6","-1/3");
    check(backwards.basePhaseTurns=="5/6"&&backwards.returnCount=="-1"&&backwards.prefixCarry=="-1","negative time requires floor division");
    check(backwards.digits==std::vector<std::uint32_t>({2,2})&&backwards.phaseTurns.back()=="53/54","negative-time compatible phases");
    for(unsigned p:{2,3,5})for(int ticks=-7;ticks<=7;++ticks){
        const auto first=advancePhaseCover(p,{1,0,1},"2/7",std::to_string(ticks)+"/3");
        const auto back=advancePhaseCover(p,first.digits,first.basePhaseTurns,std::to_string(-ticks)+"/3");
        check(back.digits==std::vector<std::uint32_t>({1,0,1})&&back.basePhaseTurns=="2/7","exact finite flow inverse law");
        for(std::size_t j=1;j<first.levels.size();++j){
            const auto difference=rationalDifference(rationalRatio(first.phaseTurns[j],"1/"+std::to_string(p)),first.phaseTurns[j-1]);
            check(difference.find('/')==std::string::npos,"adjacent flowing phases commute with covers");
        }
    }
    rejects([]{advancePhaseCover(2,{0},"1","0");});
    rejects([]{advancePhaseCover(2,{2},"0","1");});
    const auto unresolved=phaseCoverCharacter(2,{},"1/6","1/2");
    check(!unresolved.phaseTurns&&unresolved.minimumDepth==1,"unobserved root cannot be selected");
    const auto left=phaseCoverCharacter(2,{0},"1/6","1/2");
    const auto right=phaseCoverCharacter(2,{1},"1/6","1/2");
    check(left.phaseTurns=="1/12"&&right.phaseTurns=="7/12","same base state permits opposite characters");
    check(phaseCoverCharacter(2,{0},"1/6","1").phaseTurns==phaseCoverCharacter(2,{1},"1/6","1").phaseTurns,
        "integer characters descend to base circle");
    check(phaseCoverCharacter(2,{},"1/6","0").phaseTurns=="0","constant character resolved on every factor");
    for(const auto& q:{"3/4","-5/8","2","0"}){
        const auto before=phaseCoverCharacter(2,{1,0,1},"2/7",q);
        const auto moved=advancePhaseCover(2,{1,0,1},"2/7","-17/3");
        const auto after=phaseCoverCharacter(2,moved.digits,moved.basePhaseTurns,q);
        const auto phaseChange=rationalDifference(*after.phaseTurns,*before.phaseTurns);
        const auto expected=rationalRatio(q,"-3/17");
        check(rationalDifference(phaseChange,expected).find('/')==std::string::npos,"character intertwines exact flow modulo integer turns");
    }
    rejects([]{phaseCoverCharacter(2,{0},"0","1/3");});
}
void topology(){
    FiniteGraph loops{1,{{0,0},{0,0},{0,0}}};
    auto clover=graphTopology(loops);check(clover.firstBettiNumber==3&&clover.degree[0]==6,"Loop graph topology");
    FiniteGraph chain{2,{{0,0},{1,1},{0,1},{0,1}}};
    auto t=graphTopology(chain);check(t.firstBettiNumber==3&&t.connectedComponents==1&&t.degree==std::vector<std::size_t>({4,4}),"Multigraph topology");
    FiniteGraph g{5,{{0,1,3},{0,2,1},{2,1,.5},{1,3,0}}};
    auto d=graphDistances(g,0);near(*d[1],1.5);near(*d[3],1.5);check(!d[4],"Disconnected distance");
    auto bars=lowerStarPersistence({3,{{0,1},{1,2}}},{0,2,1});
    bool finite=false,essential=false;check(bars.size()==3,"H0 bars retain diagonal");
    for(auto b:bars){if(b.birth==1&&b.death&&*b.death==2)finite=true;if(b.birth==0&&!b.death)essential=true;}
    check(finite&&essential,"Persistence elder rule");
    check(graphTopology({0,{}}).connectedComponents==0,"Empty graph");
    rejects([]{graphTopology({1,{{0,1}}});});
}
void harmonic(){
    const double pi=std::acos(-1.);
    FourierSeries s{2,{{-1,.5},{1,.5}}};
    near(meanSquare(s),.5);
    auto square=multiply(s,s);near(square.coefficients.at(0),.5);
    for(double phase:{-.7,0.,1.2}){
        near(evaluate(s,phase),std::cos(phase));
        near(evaluate(timeDerivative(s),phase),-2*std::sin(phase));
        near(evaluate(timeDerivative(s,2),phase),-4*std::cos(phase));
        near(evaluate(periodicPrimitive(timeDerivative(s)),phase),evaluate(s,phase));
        near(evaluate(square,phase),std::pow(std::cos(phase),2));
    }
    near(omittedUniformBound(s,0),1);
    FourierSeries cyclic{1,{{-4,C{.3,.1}},{-1,1},{2,.2},{1,3}}};
    auto sector=projectCyclicSector(cyclic,3,-1);
    near(evaluate(sector,.4+2*pi/3),std::exp(C(0,-2*pi/3))*evaluate(sector,.4));
    for(int count:{63,64}){
        std::vector<C> samples;for(int j=0;j<count;++j){double p=2*pi*j/count;samples.push_back(2.+std::cos(3*p)+C(0,std::sin(2*p)));}
        auto recovered=periodicFourierSeries(samples,2,5);
        near(recovered.coefficients.at(0),2);near(recovered.coefficients.at(3),.5);
        near(recovered.coefficients.at(2),.5);near(recovered.coefficients.at(-2),-.5);
        near(evaluate(recovered,.21),2.+std::cos(.63)+C(0,std::sin(.42)));
    }
    rejects([]{periodicPrimitive({1,{{0,1}}});});
    rejects([]{periodicFourierSeries(std::vector<C>(4),1,2);});
}
void distributions(){
    const double pi=std::acos(-1.);
    PointDistribution delta{2,{1}},first=distributionDerivative(delta);
    near(pairTestJet(first,{4,4}),-4); // phi=x^2 at x=2.
    PointDistribution second{0,{0,0,1}};
    auto product=multiplySmooth(second,{0,1,0}); // x delta'' = -2 delta'.
    near(product.coefficients[1],-2);near(pairTestJet(product,{1,3,9}),6);
    const double a=.4;
    std::vector<C> test;for(int j=0;j<64;++j)test.push_back(std::exp(C(0,-3*2*pi*j/64)));
    near(periodicOrbitPairing(test,3,7),1);near(periodicOrbitPairing(test,3,7,false),7);
    near(periodicOrbitPairing(test,2,7),0);
    // delta(x^2+y^2-R^2): circumference / |grad H| = pi.
    near(levelSetDeltaPairing({1,1,1,1},std::vector<double>(4,pi),std::vector<double>(4,4)),pi);
    near(spectralPushforwardDensity({{2,3,2,C(3,4)}}),25./3);
    const double gaussian=std::exp(-a*a/2)/std::sqrt(2*pi);
    near(gaussianDeltaDerivative(a,1),gaussian);
    near(gaussianDeltaDerivative(a,1,1),-a*gaussian);
    near(gaussianDeltaDerivative(a,1,2),(a*a-1)*gaussian);
    near(gaussianDeltaDerivative(0,1e-300,1),0); // no 0/0.
    near(gaussianDeltaDerivative(1001,1),0);
    rejects([]{gaussianDeltaDerivative(0,1e-300,2);});
    // At z=101 a tiny epsilon can still give a very large derivative.
    rejects([]{gaussianDeltaDerivative(101e-300,1e-300,32);});
    rejects([]{levelSetDeltaPairing({1},{1},{0});});
    rejects([]{spectralPushforwardDensity({{1,0,2,1}});});
}
void resonanceAndJson(){
    const std::map<FourierIndex,C> coefficients={{{1,0},C(0,2)},{{8,-1},3},{{4,-1},2}};
    auto cover=phaseCoverForcing(2,3,1,{1},coefficients,1e-10);
    check(cover.coverDegree=="8","Cover degree");
    near(cover.response.solution.coefficients.at({1,0}),16.);
    near(cover.response.retained.coefficients.at({8,-1}),3.);
    for(const auto& mode:cover.modes)if(mode.index==FourierIndex{4,-1}) {
        check(mode.frequencyRatio=="1/2"&&mode.minimumDepth==1,"Reduced character depth");
        near(mode.detuning,-.5);near(mode.correctionMagnitude,4.);
    }
    // The solved derivative plus retained part recovers the supplied forcing.
    auto derivative=torusTimeDerivative(cover.response.solution);
    for(const auto& [index,c]:cover.response.retained.coefficients)derivative.coefficients[index]+=c;
    for(const auto& phases:std::vector<std::vector<double>>{{0.,0.},{.3,.8},{-2.,3.}})
        near(evaluate(derivative,phases),evaluate(cover.forcing,phases));
    // Raising the common cover depth must not change the physical response.
    std::map<FourierIndex,C> lifted;
    for(const auto& [index,c]:coefficients){auto k=index;k[0]*=2;lifted[k]=c;}
    auto deeper=phaseCoverForcing(2,4,1,{1},lifted,1e-10);
    near(evaluate(deeper.response.solution,{.15,.8}),evaluate(cover.response.solution,{.3,.8}));
    rejects([&]{phaseCoverForcing(4,2,1,{1},coefficients,1e-10);});
    rejects([&]{phaseCoverForcing(2,33,1,{1},coefficients,1e-10);});
    rejects([&]{phaseCoverForcing(2,3,1,{1},{{{1},1}},1e-10);});
    auto exact=harmonicResonance(2,3,3,2,10);near(exact.detuning,0);check(!exact.beatTime,"Exact beat");
    auto detuned=harmonicResonance(2,3.01,3,2,10);near(detuned.detuning,.02);near(detuned.relativePhaseDrift,.2);
    auto candidates=resonanceCandidates(1,std::sqrt(2.),100,3);check(!candidates.empty(),"Resonance candidates");
    for(std::size_t j=1;j<candidates.size();++j)check(std::abs(candidates[j].detuning)<std::abs(candidates[j-1].detuning),"Record detuning");
    check(resonanceCandidates(1,1000,5).size()==1,"Extreme ratio bound");
    namespace json=dynsyslab::core::json;
    auto q=json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"padic","a":"4/9","b":"1/9","prime":3,"precision":4})");
    auto result=evaluateQuery(q);check(result.at("result").at("distance").string()=="3","JSON p-adic path");
    q.object()["prime"]=3.5;rejects([&]{evaluateQuery(q);});
    q.object()["prime"]="3";rejects([&]{evaluateQuery(q);});
}
}
int main(){try{arithmetic();topology();harmonic();distributions();resonanceAndJson();
    std::cout<<"Exact arithmetic, topology, Fourier, resonance, distributions and query checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
