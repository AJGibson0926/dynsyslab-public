#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/EigenvalueLattice.hpp>
#include <dynsyslab/math/ReebGraph.hpp>
#include <dynsyslab/math/ShellPacket.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
using namespace dynsyslab::math;
using C=std::complex<double>;
namespace dynsyslab::math {
std::optional<core::json::Value> evaluateSpectralDiagnosticsQuery(const std::string&,const core::json::Value&);
}
namespace json=dynsyslab::core::json;
static const double pi=std::acos(-1.);
static int checks=0;
#define CHECK(x) do{++checks;if(!(x))throw std::runtime_error(std::string("line ")+std::to_string(__LINE__)+": "#x);}while(0)
static bool near(double a,double b,double tol=1e-10){return std::abs(a-b)<=tol*(1+std::abs(b));}
template<class F>static bool rejects(F f){try{f();}catch(const std::exception&){return true;}return false;}

static void latticeTests(){
    const C a(-.1,2),b(-.3,3.1);
    const std::vector<C> values{0,a,std::conj(a),b,std::conj(b),2.*a,a+b,3.*b+1e-5,C(-5,.123)};
    EigenvalueLatticeOptions o;o.maximumOrder=3;o.maximumPrincipals=4;o.tolerance=1e-4;
    auto r=eigenvalueLattice(values,o);
    CHECK(!r.rationalityCertified);CHECK(r.principal.size()==4);
    CHECK(r.entries[0].role==LatticeRole::Constant);
    CHECK(r.entries[5].role==LatticeRole::Combination&&r.entries[5].order==2);
    CHECK(r.entries[6].role==LatticeRole::Combination&&r.entries[6].order==2);
    CHECK(r.entries[7].role==LatticeRole::Combination&&r.entries[7].order==3);
    CHECK(near(r.entries[7].detuningMagnitude,1e-5,1e-6));
    CHECK(r.entries[8].role==LatticeRole::Unmatched);
    const auto slot=std::find(r.principal.begin(),r.principal.end(),std::size_t(1))-r.principal.begin();
    CHECK(r.entries[5].multiIndex[std::size_t(slot)]==2);
    for(const auto& e:r.entries)if(e.role==LatticeRole::Combination){
        C sum{};for(std::size_t j=0;j<r.principal.size();++j)sum+=double(e.multiIndex[j])*values[r.principal[j]];
        CHECK(std::abs(values[e.index]-sum-e.detuning)<1e-12);
    }
    // Tighter tolerance must refuse the detuned third harmonic.
    o.tolerance=1e-6;CHECK(eigenvalueLattice(values,o).entries[7].role==LatticeRole::Unmatched);
    // Discrete aliasing: 3*2.5 wraps past pi, one full turn.
    const double angle=2.5;
    EigenvalueLatticeOptions d;d.clock=EigenvalueClock::Discrete;d.maximumOrder=3;d.tolerance=1e-9;d.principalIndices={0};
    auto aliased=eigenvalueLattice({std::polar(.99,angle),std::pow(std::polar(.99,angle),3),C{}},d);
    CHECK(aliased.entries[1].role==LatticeRole::Combination&&aliased.entries[1].order==3);
    CHECK(aliased.entries[1].aliasTurns==-1);CHECK(aliased.entries[2].role==LatticeRole::Excluded);
    CHECK(aliased.chanceCoverage.has_value());
    // The same angles without a clock period are not a lattice relation.
    EigenvalueLatticeOptions unwrapped;unwrapped.maximumOrder=3;unwrapped.tolerance=1e-9;unwrapped.principalIndices={0};
    CHECK(eigenvalueLattice({C(0,angle),C(0,3*angle-2*pi)},unwrapped).entries[1].role==LatticeRole::Unmatched);
    unwrapped.sampleInterval=1;
    CHECK(eigenvalueLattice({C(0,angle),C(0,3*angle-2*pi)},unwrapped).entries[1].role==LatticeRole::Combination);
    FrequencyResonanceOptions f;f.source="hankel";f.maximumOrder=3;f.tolerance=1e-9;
    auto report=frequencyResonanceReport({.7,1.3,2.0,2.6,-1.3},f);
    CHECK(report.lattice.principal.size()==2&&report.foldedNegative==std::vector<std::size_t>{4});
    CHECK(report.lattice.entries[2].order==2&&report.lattice.entries[3].order==2);
    CHECK(report.lattice.entries[4].role==LatticeRole::Combination&&report.lattice.entries[4].order==1);
    FrequencyResonanceOptions harmonic;harmonic.maximumOrder=3;harmonic.tolerance=1e-9;harmonic.principalIndices={0,1};
    CHECK(!frequencyResonanceReport({1.,2.},harmonic).pairs.empty());
    CHECK(rejects([&]{FrequencyResonanceOptions x;x.source="other";frequencyResonanceReport({1.},x);}));
    CHECK(rejects([&]{EigenvalueLatticeOptions x;x.maximumOrder=0;eigenvalueLattice({C(1,0)},x);}));
    CHECK(rejects([&]{EigenvalueLatticeOptions x;x.clock=EigenvalueClock::Discrete;x.sampleInterval=.1;eigenvalueLattice({C(1,0)},x);}));
}
// The literal Rigged DMD expression, independent of RationalCircleKernel.
static double literalWeight(double alpha,double theta,double eps,int m){
    std::vector<C> p,r;
    for(int j=1;j<=m;++j)p.emplace_back(2.*j/(m+1)-1,1);
    for(int j=0;j<m;++j){C a=1;for(int k=0;k<m;++k)if(k!=j)a*=-p[k]/(p[j]-p[k]);r.push_back(a);}
    const C u=std::exp(C(0,alpha));auto f=[&](C z){return (u+z)/(u-z);};
    C sum{};
    for(int j=0;j<m;++j)sum+=r[j]*f(std::exp(C(0,theta)-C(0,eps)*p[j]))-std::conj(r[j])*f(std::exp(C(0,theta)-C(0,eps)*std::conj(p[j])));
    return (-sum/(4*pi)).real();
}
static void kernelTests(){
    for(std::size_t m=1;m<=6;++m)for(double eps:{.05,.3,.9}){
        const RationalCircleKernel k(eps,m);
        for(int q=0;q<int(m);++q){C moment{};for(std::size_t j=0;j<m;++j)moment+=k.residues()[j]*std::pow(k.poles()[j],q);
            CHECK(std::abs(moment-(q==0?C(1):C(0)))<1e-9);}
        const int samples=40000;double mass=0;
        for(int i=0;i<samples;++i)mass+=k(2*pi*(i+.5)/samples-pi)*2*pi/samples;
        CHECK(near(mass,1,1e-8));
        for(double alpha:{-2.,.1,1.7})for(double theta:{-.4,.1,2.9})
            CHECK(near(shellPacketWeight(alpha,theta,k),literalWeight(alpha,theta,eps,int(m)),1e-8));
    }
    for(double eps:{.05,.3,.9}){
        const RationalCircleKernel k(eps,1);const double q=std::exp(-eps);
        for(double d:{0.,.2,1.,3.})CHECK(near(k(d),(1-q*q)/(1+q*q-2*q*std::cos(d))/(2*pi),1e-10));
        CHECK(near(kernelHalfWidth(k),std::acos(2-std::cosh(eps)),1e-9));
    }
    CHECK(near(kernelHalfWidth(RationalCircleKernel(3.,1)),pi)); // Never halves on the circle.
}
static void packetTests(){
    const double dt=.5,Omega=1.2,theta=.3;const RationalCircleKernel k(.2,3);
    FourierSeries s{Omega,{{-2,C(.2,.1)},{1,C(1,-.5)},{3,C(.3,.4)}}};
    auto base=predictedShellPacket(s,dt,theta,k),delayed=predictedShellPacket(s,dt,theta,k,4);
    for(const auto& [n,c]:base.coefficients){
        CHECK(std::abs(c-k(theta-n*Omega*dt)*s.coefficients.at(n))<1e-13);
        // Pure rotation: the packet equals the lag-domain Fourier sum of the kernel.
        C lagSum=k.timeMultiplier(0);
        for(int l=1;l<=400;++l)lagSum+=2.*k.timeMultiplier(l)*std::cos(l*(theta-n*Omega*dt));
        CHECK(std::abs(lagSum/(2*pi)-k(theta-n*Omega*dt))<1e-9);
    }
    for(double phi:{-1.,.4,2.2})CHECK(std::abs(evaluate(delayed,phi)-evaluate(base,phi+Omega*4*dt))<1e-12);
    std::vector<OrbitShell> shells;
    for(int i=0;i<=10;++i){const double E=i/10.;shells.push_back({E,FourierSeries{1+E,{{1,C(1,0)},{3,C(.5,0)}}},.1});}
    auto p=predictShellPackets(shells,dt,theta,k);
    double density=0;for(const auto& o:shells)for(const auto& [n,c]:o.series.coefficients)
        density+=o.weight*std::norm(c)*k(theta-n*o.series.angularFrequency*dt);
    CHECK(near(p.density,density,1e-12));CHECK(near(p.totalMass,1.1*1.25,1e-12));
    // Omega=1+E, dt=1/2: only n=3 reaches angle 2 on E in [0,1], at E=1/3.
    const double aliasAngle=2;
    const auto catalogue=shellAliasCatalogue(shells,dt,aliasAngle,3,p.halfWidth);
    CHECK(catalogue.size()==1&&catalogue[0].harmonic==3&&catalogue[0].turns==0);
    CHECK(near(catalogue[0].energy,1./3,1e-12));
    // A wrapped branch: n=3 at angle 2-2*pi is the same circle point, one turn lower.
    const auto wrapped=shellAliasCatalogue(shells,dt,aliasAngle-2*pi,3,p.halfWidth);
    CHECK(wrapped.size()==1&&wrapped[0].turns==1&&near(wrapped[0].energy,1./3,1e-12));
    for(const auto& a:catalogue){
        CHECK(near(a.harmonic*(1+a.energy)*dt,aliasAngle+2*pi*double(a.turns),1e-10));
        CHECK(a.ridgeWidth&&near(*a.ridgeWidth,p.halfWidth/(std::abs(a.harmonic)*dt),1e-9));
    }
    std::vector<C> x{C(1,1),C(2,-1),C(.5,.2),C(-1,.7)},rotated;
    for(auto z:x)rotated.push_back(z*std::exp(C(0,.8)));
    auto same=scoreShellPacket(x,x,{0,0,1,1});
    CHECK(same.relativeL2Error<1e-14&&near(same.complexCorrelation.real(),1));
    CHECK(near(*same.magnitudeCorrelation,1)&&near(*same.energyProfileCorrelation,1));
    auto phase=scoreShellPacket(rotated,x);
    CHECK(phase.relativeL2Error>.5&&near(phase.complexCorrelationMagnitude,1)&&near(std::arg(phase.complexCorrelation),.8));
    CHECK(near(*phase.magnitudeCorrelation,1));
    CHECK(!scoreShellPacket({C(1),C(1)},{C(1),C(1)}).magnitudeCorrelation.has_value());
    CHECK(rejects([&]{scoreShellPacket(x,{C(1)});}));
    auto oneSided=diagnoseDelayTarget(0,20,0,k),centred=diagnoseDelayTarget(0,20,10,k);
    CHECK(oneSided.oneSided&&!oneSided.warning.empty()&&!centred.oneSided&&centred.warning.empty());
    CHECK(centred.symmetricTailBound<oneSided.symmetricTailBound);
    CHECK(diagnoseDelayTarget(0,20,30,k).oneSided);
    const auto estimate=resolvableShellEstimate(.5,.6,.2);
    CHECK(near(estimate.lastHarmonic,1.5)&&estimate.empirical);
    CHECK(rejects([&]{shellAliasCatalogue({shells[1],shells[0]},dt,theta,3,.1);}));
}
template<class F>static ScalarGrid2D sampled(std::size_t nx,std::size_t ny,bool px,bool py,F f){
    ScalarGrid2D g;g.nx=nx;g.ny=ny;g.periodicX=px;g.periodicY=py;g.values.reserve(nx*ny);
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i)g.values.push_back(f(i,j));
    return g;
}
static int interior(const ScalarTopology& t,GridCriticalKind kind){
    return int(std::count_if(t.criticalPoints.begin(),t.criticalPoints.end(),[&](const auto& p){return p.kind==kind&&!p.onBoundary;}));
}
static int essential(const ScalarTopology& t,unsigned dimension){
    return int(std::count_if(t.persistence.begin(),t.persistence.end(),[&](const auto& p){return p.dimension==dimension&&!p.deathVertex;}));
}
static void reebTests(){
    const auto well=sampled(81,41,false,false,[](auto i,auto j){const double x=-2+.05*i,y=-1+.05*j;return (x*x-1)*(x*x-1)+y*y;});
    const auto t=analyzeScalarGrid(well);
    CHECK(t.eulerCharacteristic==1);CHECK(interior(t,GridCriticalKind::Minimum)==2);
    CHECK(interior(t,GridCriticalKind::Saddle)==1);CHECK(t.reebFirstBetti==0);
    const auto bars=std::count_if(t.persistence.begin(),t.persistence.end(),[](const auto& p){return p.dimension==0&&p.death&&*p.death-p.birth>.5;});
    CHECK(bars==1);
    const auto levels=std::count_if(t.separatrices.begin(),t.separatrices.end(),[](const auto& s){return !s.onBoundary&&near(s.value,1,1e-12);});
    CHECK(levels==1);
    const auto leaves=std::count_if(t.mergeTree.nodes.begin(),t.mergeTree.nodes.end(),[](const auto& n){return n.lowerDegree==0;});
    CHECK(leaves==2);CHECK(graphTopology(toFiniteGraph(t.mergeTree)).firstBettiNumber==0);
    CHECK(graphTopology(toFiniteGraph(t.reebGraph)).connectedComponents==1);
    const double dq=2*pi/64,dp=.125;
    const auto pendulum=sampled(64,41,true,false,[&](auto i,auto j){const double q=-pi+dq*i,p=-2.5+dp*j;return p*p/2-std::cos(q);});
    const auto h=analyzeScalarGrid(pendulum);
    CHECK(h.eulerCharacteristic==0);CHECK(interior(h,GridCriticalKind::Minimum)==1);
    CHECK(interior(h,GridCriticalKind::Saddle)==1);CHECK(h.reebFirstBetti==0);CHECK(essential(h,1)==1);
    CHECK(std::count_if(h.separatrices.begin(),h.separatrices.end(),[](const auto& s){return !s.onBoundary&&near(s.value,1,1e-12);})==1);
    const auto torus=sampled(48,48,true,true,[&](auto i,auto j){return std::cos(2*pi*i/48+.01)+.7*std::cos(2*pi*j/48+.02);});
    const auto r=analyzeScalarGrid(torus);
    CHECK(r.eulerCharacteristic==0);CHECK(interior(r,GridCriticalKind::Minimum)==1);
    CHECK(interior(r,GridCriticalKind::Saddle)==2);CHECK(interior(r,GridCriticalKind::Maximum)==1);
    CHECK(r.reebFirstBetti==1);CHECK(essential(r,0)==1&&essential(r,1)==2&&essential(r,2)==1);
    CHECK(rejects([&]{ScalarGrid2D x;x.nx=2;x.ny=3;x.periodicX=true;x.values.assign(6,0.);analyzeScalarGrid(x);}));
    CHECK(rejects([&]{ScalarGrid2D x;x.nx=2;x.ny=2;x.values={0,1,2};analyzeScalarGrid(x);}));
}
// Independent check: standard boundary-matrix reduction (Edelsbrunner et al. 2002).
using Triple=std::tuple<unsigned,std::size_t,std::size_t>;
static std::multiset<Triple> reducedPairs(const ScalarGrid2D& g){
    using S=std::vector<std::size_t>;std::set<S> simplices;
    const std::size_t cx=g.periodicX?g.nx:g.nx-1,cy=g.periodicY?g.ny:g.ny-1;
    auto add=[&](S t){std::sort(t.begin(),t.end());simplices.insert(t);
        for(std::size_t k=0;k<3;++k){S e=t;e.erase(e.begin()+k);simplices.insert(e);simplices.insert({e[0]});simplices.insert({e[1]});}};
    for(std::size_t j=0;j<cy;++j)for(std::size_t i=0;i<cx;++i){
        const auto i1=(i+1)%g.nx,j1=(j+1)%g.ny;
        const auto a=i+g.nx*j,b=i1+g.nx*j,c=i+g.nx*j1,d=i1+g.nx*j1;add({a,b,d});add({a,c,d});
    }
    std::vector<S> order(simplices.begin(),simplices.end());
    auto top=[&](const S& s){return *std::max_element(s.begin(),s.end(),[&](auto x,auto y){return g.values[x]<g.values[y];});};
    std::stable_sort(order.begin(),order.end(),[&](const S& x,const S& y){
        const double a=g.values[top(x)],b=g.values[top(y)];return a!=b?a<b:x.size()<y.size();});
    std::map<S,std::size_t> position;for(std::size_t k=0;k<order.size();++k)position[order[k]]=k;
    std::vector<std::set<std::size_t>> column(order.size());
    for(std::size_t k=0;k<order.size();++k)if(order[k].size()>1)for(std::size_t f=0;f<order[k].size();++f){
        S face=order[k];face.erase(face.begin()+f);column[k].insert(position.at(face));}
    std::map<std::size_t,std::size_t> pivot;std::set<std::size_t> paired;std::multiset<Triple> out;
    for(std::size_t k=0;k<order.size();++k){
        while(!column[k].empty()){
            const auto low=*column[k].rbegin();const auto it=pivot.find(low);if(it==pivot.end())break;
            for(auto x:column[it->second])if(!column[k].erase(x))column[k].insert(x);
        }
        if(column[k].empty())continue;
        const auto low=*column[k].rbegin();pivot[low]=k;paired.insert(low);paired.insert(k);
        const auto birth=top(order[low]),death=top(order[k]);
        if(birth!=death)out.insert({unsigned(order[low].size()-1),birth,death});
    }
    for(std::size_t k=0;k<order.size();++k)if(!paired.count(k))
        out.insert({unsigned(order[k].size()-1),top(order[k]),std::size_t(-1)});
    return out;
}
// Contours at a regular level: triangles crossing it, joined across crossed edges.
static std::size_t contourCount(const ScalarGrid2D& g,double level){
    using E=std::pair<std::size_t,std::size_t>;
    std::vector<std::array<std::size_t,3>> crossing;std::map<E,std::vector<std::size_t>> byEdge;
    const std::size_t cx=g.periodicX?g.nx:g.nx-1,cy=g.periodicY?g.ny:g.ny-1;
    auto add=[&](std::size_t a,std::size_t b,std::size_t c){
        const std::array<std::size_t,3> v{a,b,c};int below=0;for(auto x:v)below+=g.values[x]<level;
        if(below==0||below==3)return;
        for(int k=0;k<3;++k){
            auto x=v[k],y=v[(k+1)%3];
            if((g.values[x]<level)!=(g.values[y]<level))byEdge[{std::min(x,y),std::max(x,y)}].push_back(crossing.size());
        }
        crossing.push_back(v);
    };
    for(std::size_t j=0;j<cy;++j)for(std::size_t i=0;i<cx;++i){
        const auto i1=(i+1)%g.nx,j1=(j+1)%g.ny;
        const auto a=i+g.nx*j,b=i1+g.nx*j,c=i+g.nx*j1,d=i1+g.nx*j1;add(a,b,d);add(a,c,d);
    }
    std::vector<std::size_t> parent(crossing.size());for(std::size_t k=0;k<parent.size();++k)parent[k]=k;
    auto root=[&](std::size_t x){while(parent[x]!=x)x=parent[x]=parent[parent[x]];return x;};
    for(const auto& [edge,list]:byEdge)for(std::size_t k=1;k<list.size();++k)parent[root(list[k])]=root(list[0]);
    std::size_t count=0;for(std::size_t k=0;k<parent.size();++k)count+=root(k)==k;
    return count;
}
static void persistenceCrossCheck(){
    std::mt19937 random(20261123);std::uniform_real_distribution<double> uniform(0,1);
    for(int trial=0;trial<40;++trial)for(int mode=0;mode<4;++mode){
        const std::size_t nx=3+random()%5,ny=3+random()%5;
        const auto g=sampled(nx,ny,mode&1,mode&2,[&](auto,auto){return uniform(random);});
        const auto t=analyzeScalarGrid(g);
        std::multiset<Triple> pairs;
        for(const auto& p:t.persistence)pairs.insert({p.dimension,p.birthVertex,p.deathVertex?*p.deathVertex:std::size_t(-1)});
        CHECK(pairs==reducedPairs(g));
        const int chi=(mode==3)?0:(mode==0?1:0);CHECK(t.eulerCharacteristic==chi);
        // The quotient onto the Reeb graph surjects on pi_1: a disc has no loop,
        // and a torus or annulus has at most one.
        CHECK(graphTopology(toFiniteGraph(t.reebGraph)).connectedComponents==1);
        CHECK(t.reebFirstBetti<=(mode==0?0u:1u));
        // Independent check: arcs crossing each regular level equal its contours.
        auto levels=g.values;std::sort(levels.begin(),levels.end());
        for(std::size_t k=0;k+1<levels.size();++k){
            const double level=(levels[k]+levels[k+1])/2;
            const auto arcs=std::count_if(t.reebGraph.arcs.begin(),t.reebGraph.arcs.end(),[&](const auto& a){
                return t.reebGraph.nodes[a.lowerNode].value<level&&level<t.reebGraph.nodes[a.upperNode].value;});
            CHECK(std::size_t(arcs)==contourCount(g,level));
        }
    }
}
static void queryTests(){
    const auto lattice=evaluateSpectralDiagnosticsQuery("eigenvalue_lattice",json::parse(
        R"({"eigenvalues":[[0,1],[0,2]],"maximum_order":2,"principal_indices":[0],"tolerance":1e-9})"));
    CHECK(lattice&&lattice->at("entries").array()[1].at("role").string()=="combination");
    CHECK(!lattice->at("rationality_certified").boolean(true));
    const auto frequency=evaluateSpectralDiagnosticsQuery("frequency_resonance_report",json::parse(
        R"({"frequencies":[1,2,-1],"source":"hankel","tolerance":1e-9})"));
    CHECK(frequency&&frequency->at("folded_negative").array().size()==1);
    const auto packet=evaluateSpectralDiagnosticsQuery("shell_packet",json::parse(
        R"({"bandwidth":0.2,"order":2,"sample_interval":0.5,"angle":0.3,"target_delay":0,"maximum_harmonic":2,
            "delay_window":{"first_lag":0,"last_lag":10},"resolvable_constant":0.5,"angular_scale":0.6,
            "shells":[{"energy":0,"angular_frequency":1,"coefficients":[[1,1,0]]},
                      {"energy":1,"angular_frequency":2,"coefficients":[[1,1,0]]}]})"));
    CHECK(packet&&packet->at("delay_diagnostic").at("one_sided").boolean());
    CHECK(packet->at("resolvable_shell_estimate").at("empirical").boolean());
    const auto score=evaluateSpectralDiagnosticsQuery("shell_score",json::parse(
        R"({"computed":[[1,0],[0,1]],"predicted":[[1,0],[0,1]]})"));
    CHECK(score&&near(score->at("relative_l2_error").number(),0));
    const auto reeb=evaluateSpectralDiagnosticsQuery("reeb_graph",json::parse(
        R"({"nx":3,"ny":3,"values":[1,2,3,4,0,5,6,7,8]})"));
    CHECK(reeb&&near(reeb->at("euler_characteristic").number(),1));
    CHECK(!evaluateSpectralDiagnosticsQuery("unknown",json::parse("{}")));
    CHECK(rejects([]{evaluateSpectralDiagnosticsQuery("reeb_graph",json::parse(R"({"nx":3,"ny":3,"values":[1]})"));}));
    CHECK(rejects([]{evaluateSpectralDiagnosticsQuery("eigenvalue_lattice",json::parse(R"({"eigenvalues":[[1,2,3]]})"));}));
    CHECK(rejects([]{evaluateSpectralDiagnosticsQuery("shell_packet",json::parse(R"({"bandwidth":-1})"));}));
}
int main(){
    try{
        latticeTests();kernelTests();packetTests();reebTests();persistenceCrossCheck();queryTests();
        std::cout<<"spectral diagnostics tests passed ("<<checks<<" checks)\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<"\n";return 1;}
}
