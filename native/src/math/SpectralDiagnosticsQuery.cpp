#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/EigenvalueLattice.hpp>
#include <dynsyslab/math/ReebGraph.hpp>
#include <dynsyslab/math/ShellPacket.hpp>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
// Thin JSON adapter. Algorithms and citations live in EigenvalueLattice,
// ShellPacket and ReebGraph. Returns nullopt for an operation it does not own.
namespace dynsyslab::math {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;using C=std::complex<double>;
double num(const J& v){
    if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))
        throw std::invalid_argument("Expected a finite JSON number");
    return v.number();
}
long long integer(const J& v,long long lo,long long hi){
    const double n=num(v);
    if(n!=std::floor(n)||n<double(lo)||n>double(hi))throw std::invalid_argument("Integer argument out of range");
    return (long long)n;
}
bool flag(const J& q,const char* key){
    const auto& v=q.at(key);if(v.isNull())return false;
    if(!std::holds_alternative<bool>(v.storage()))throw std::invalid_argument("Expected a JSON boolean");
    return v.boolean();
}
const A& array(const J& v){
    if(!std::holds_alternative<A>(v.storage()))throw std::invalid_argument("Expected a JSON array");
    if(v.array().size()>1000000)throw std::invalid_argument("Query array exceeds 1000000 elements");
    return v.array();
}
const J& object(const J& v){
    if(!std::holds_alternative<O>(v.storage()))throw std::invalid_argument("Expected a JSON object");
    return v;
}
C complexValue(const J& v){
    if(std::holds_alternative<double>(v.storage()))return num(v);
    const auto& a=array(v);if(a.size()!=2)throw std::invalid_argument("Complex value needs [real,imaginary]");
    return {num(a[0]),num(a[1])};
}
std::vector<C> complexes(const J& v){std::vector<C> out;for(const auto& x:array(v))out.push_back(complexValue(x));return out;}
std::vector<double> reals(const J& v){std::vector<double> out;for(const auto& x:array(v))out.push_back(num(x));return out;}
std::vector<std::size_t> indices(const J& v){
    std::vector<std::size_t> out;for(const auto& x:array(v))out.push_back(std::size_t(integer(x,0,1000000)));return out;
}
double numberOr(const J& q,const char* key,double fallback){return q.at(key).isNull()?fallback:num(q.at(key));}
std::size_t sizeOr(const J& q,const char* key,std::size_t fallback,long long hi){
    return q.at(key).isNull()?fallback:std::size_t(integer(q.at(key),1,hi));
}
J encode(C z){return A{z.real(),z.imag()};}
J optional(const std::optional<double>& v){return v?J(*v):J();}
template<class T>J numbers(const std::vector<T>& values){A out;for(auto x:values)out.emplace_back(double(x));return out;}
J seriesJson(const FourierSeries& s){
    A rows;for(const auto& [n,c]:s.coefficients)rows.emplace_back(A{double(n),c.real(),c.imag()});
    return O{{"angular_frequency",s.angularFrequency},{"coefficients",rows}};
}
FourierSeries seriesFrom(const J& q){
    FourierSeries s{num(q.at("angular_frequency")),{}};
    for(const auto& row:array(q.at("coefficients"))){
        const auto& a=array(row);
        if(a.size()!=3)throw std::invalid_argument("Fourier coefficient row needs [harmonic,real,imaginary]");
        if(!s.coefficients.emplace(int(integer(a[0],-1000000,1000000)),C(num(a[1]),num(a[2]))).second)
            throw std::invalid_argument("Duplicate harmonic");
    }
    return s;
}
J latticeJson(const EigenvalueLatticeReport& r){
    A entries;
    for(const auto& e:r.entries)entries.emplace_back(O{{"index",double(e.index)},{"eigenvalue",encode(e.eigenvalue)},
        {"coordinate",encode(e.coordinate)},{"role",latticeRoleName(e.role)},{"multi_index",numbers(e.multiIndex)},
        {"order",double(e.order)},{"detuning",encode(e.detuning)},{"detuning_magnitude",e.detuningMagnitude},
        {"accepted_tolerance",e.acceptedTolerance},{"phase_drift",e.phaseDrift},{"alias_turns",double(e.aliasTurns)}});
    return O{{"entries",entries},{"principal",numbers(r.principal)},{"alias_period",optional(r.aliasPeriod)},
        {"combination_count",r.combinationCount},{"chance_coverage",optional(r.chanceCoverage)},
        {"rationality_certified",r.rationalityCertified},
        {"scope","Numerical candidates within the declared tolerance; no rationality, eigenfunction-product or generator claim"}};
}
J eigenvalueLatticeQuery(const J& q){
    EigenvalueLatticeOptions o;
    const auto clock=q.at("clock").isNull()?std::string("continuous"):q.at("clock").string();
    if(clock=="discrete")o.clock=EigenvalueClock::Discrete;
    else if(clock!="continuous")throw std::invalid_argument("clock must be continuous or discrete");
    if(!q.at("sample_interval").isNull())o.sampleInterval=num(q.at("sample_interval"));
    o.maximumOrder=sizeOr(q,"maximum_order",4,64);o.maximumPrincipals=sizeOr(q,"maximum_principals",4,32);
    o.tolerance=numberOr(q,"tolerance",1e-3);o.relativeTolerance=numberOr(q,"relative_tolerance",0);
    o.observationTime=numberOr(q,"horizon",0);
    if(!q.at("principal_indices").isNull())o.principalIndices=indices(q.at("principal_indices"));
    return latticeJson(eigenvalueLattice(complexes(q.at("eigenvalues")),o));
}
J frequencyQuery(const J& q){
    FrequencyResonanceOptions o;
    if(!q.at("source").isNull())o.source=q.at("source").string();
    o.radiansPerSample=flag(q,"radians_per_sample");
    if(!q.at("sample_interval").isNull())o.sampleInterval=num(q.at("sample_interval"));
    o.maximumOrder=sizeOr(q,"maximum_order",4,64);o.maximumPrincipals=sizeOr(q,"maximum_principals",4,32);
    o.tolerance=numberOr(q,"tolerance",1e-3);o.relativeTolerance=numberOr(q,"relative_tolerance",0);
    o.observationTime=numberOr(q,"horizon",0);
    if(!q.at("principal_indices").isNull())o.principalIndices=indices(q.at("principal_indices"));
    const auto r=frequencyResonanceReport(reals(q.at("frequencies")),o);
    A pairs;
    for(const auto& p:r.pairs)pairs.emplace_back(O{{"principal_a",double(p.principalA)},{"principal_b",double(p.principalB)},
        {"harmonic_a",double(p.resonance.harmonicA)},{"harmonic_b",double(p.resonance.harmonicB)},
        {"detuning",p.resonance.detuning},{"phase_drift",p.resonance.relativePhaseDrift},
        {"beat_time",optional(p.resonance.beatTime)}});
    return O{{"source",r.source},{"lattice",latticeJson(r.lattice)},{"folded_negative",numbers(r.foldedNegative)},
        {"pairs",pairs},{"exact_rationality_certified",false}};
}
J shellPacketQuery(const J& q){
    const RationalCircleKernel kernel(num(q.at("bandwidth")),sizeOr(q,"order",1,8));
    const double dt=num(q.at("sample_interval")),theta=num(q.at("angle")),delay=numberOr(q,"target_delay",0);
    std::vector<OrbitShell> shells;
    for(const auto& s:array(q.at("shells"))){
        object(s);shells.push_back({num(s.at("energy")),seriesFrom(s),numberOr(s,"weight",1)});
    }
    const auto p=predictShellPackets(shells,dt,theta,kernel,delay);
    A packets;for(const auto& s:p.packets)packets.push_back(seriesJson(s));
    O out{{"half_width",p.halfWidth},{"density",p.density},{"total_mass",p.totalMass},
        {"zero_harmonic_mass",p.zeroHarmonicMass},{"shell_energy",numbers(p.shellEnergy)},{"packets",packets},
        {"scope","Finite-harmonic prediction from supplied orbit tables; omitted harmonics and table error excluded"}};
    if(!q.at("maximum_harmonic").isNull()){
        A aliases;
        for(const auto& a:shellAliasCatalogue(shells,dt,theta,int(integer(q.at("maximum_harmonic"),1,100000)),p.halfWidth))
            aliases.emplace_back(O{{"harmonic",double(a.harmonic)},{"turns",double(a.turns)},{"energy",a.energy},
                {"angular_frequency",a.angularFrequency},{"frequency_slope",a.frequencySlope},{"ridge_width",optional(a.ridgeWidth)}});
        out["alias_catalogue"]=aliases;
    }
    if(!q.at("delay_window").isNull()){
        const auto& w=object(q.at("delay_window"));
        if(delay!=std::floor(delay))throw std::invalid_argument("Delay diagnostic needs an integer target_delay");
        const auto d=diagnoseDelayTarget(integer(w.at("first_lag"),-1000000000,1000000000),
            integer(w.at("last_lag"),-1000000000,1000000000),(long long)delay,kernel,numberOr(w,"asymmetry_tolerance",.25));
        out["delay_diagnostic"]=O{{"lags_before",double(d.lagsBefore)},{"lags_after",double(d.lagsAfter)},
            {"retained_symmetric_lag",double(d.retainedSymmetricLag)},{"asymmetry",d.asymmetry},
            {"symmetric_tail_bound",d.symmetricTailBound},{"one_sided",d.oneSided},{"warning",d.warning}};
    }
    if(!q.at("resolvable_constant").isNull()){
        const auto e=resolvableShellEstimate(num(q.at("resolvable_constant")),num(q.at("angular_scale")),p.halfWidth);
        out["resolvable_shell_estimate"]=O{{"last_harmonic",e.lastHarmonic},{"empirical",e.empirical},{"scope",e.scope}};
    }
    return out;
}
J shellScoreQuery(const J& q){
    const auto s=scoreShellPacket(complexes(q.at("computed")),complexes(q.at("predicted")),
        q.at("shell_index").isNull()?std::vector<std::size_t>{}:indices(q.at("shell_index")),
        q.at("weights").isNull()?std::vector<double>{}:reals(q.at("weights")));
    return O{{"relative_l2_error",s.relativeL2Error},{"complex_correlation",encode(s.complexCorrelation)},
        {"complex_correlation_magnitude",s.complexCorrelationMagnitude},
        {"magnitude_correlation",optional(s.magnitudeCorrelation)},
        {"energy_profile_correlation",optional(s.energyProfileCorrelation)}};
}
J levelGraphJson(const ScalarLevelGraph& g){
    A nodes,arcs;
    for(const auto& n:g.nodes)nodes.emplace_back(O{{"vertex",double(n.vertex)},{"value",n.value},
        {"kind",gridCriticalKindName(n.kind)},{"on_boundary",n.onBoundary},
        {"lower_degree",double(n.lowerDegree)},{"upper_degree",double(n.upperDegree)}});
    for(const auto& a:g.arcs)arcs.emplace_back(A{double(a.lowerNode),double(a.upperNode),
        std::abs(g.nodes[a.upperNode].value-g.nodes[a.lowerNode].value)});
    // "edges" rows are [from,to,length], the format of the "graph" operation.
    return O{{"vertices",double(g.nodes.size())},{"nodes",nodes},{"edges",arcs}};
}
J reebQuery(const J& q){
    ScalarGrid2D g;
    g.nx=std::size_t(integer(q.at("nx"),2,1000000));g.ny=std::size_t(integer(q.at("ny"),2,1000000));
    g.values=reals(q.at("values"));g.periodicX=flag(q,"periodic_x");g.periodicY=flag(q,"periodic_y");
    const auto t=analyzeScalarGrid(g,numberOr(q,"minimum_persistence",0));
    A critical,pairs,separatrices;
    for(const auto& p:t.criticalPoints)critical.emplace_back(O{{"vertex",double(p.vertex)},{"i",double(p.i)},{"j",double(p.j)},
        {"value",p.value},{"kind",gridCriticalKindName(p.kind)},{"lower_link",double(p.lowerLinkComponents)},
        {"upper_link",double(p.upperLinkComponents)},{"multiplicity",double(p.multiplicity)},{"on_boundary",p.onBoundary}});
    for(const auto& p:t.persistence)pairs.emplace_back(O{{"dimension",double(p.dimension)},{"birth_vertex",double(p.birthVertex)},
        {"birth",p.birth},{"death_vertex",p.deathVertex?J(double(*p.deathVertex)):J()},{"death",optional(p.death)}});
    for(const auto& s:t.separatrices)separatrices.emplace_back(O{{"vertex",double(s.vertex)},{"value",s.value},
        {"on_boundary",s.onBoundary},{"lower_arcs",double(s.lowerArcs)},{"upper_arcs",double(s.upperArcs)}});
    return O{{"critical_points",critical},{"euler_characteristic",double(t.eulerCharacteristic)},{"persistence",pairs},
        {"merge_tree",levelGraphJson(t.mergeTree)},{"reeb_graph",levelGraphJson(t.reebGraph)},
        {"reeb_first_betti",double(t.reebFirstBetti)},{"separatrix_levels",separatrices},
        {"scope","Topology of the piecewise-linear interpolant on the Freudenthal triangulation; not a certificate for the sampled function"}};
}
}
std::optional<core::json::Value> evaluateSpectralDiagnosticsQuery(const std::string& op,const core::json::Value& q){
    if(op!="eigenvalue_lattice"&&op!="frequency_resonance_report"&&op!="shell_packet"&&op!="shell_score"&&op!="reeb_graph")
        return std::nullopt;
    object(q);
    if(op=="eigenvalue_lattice")return eigenvalueLatticeQuery(q);
    if(op=="frequency_resonance_report")return frequencyQuery(q);
    if(op=="shell_packet")return shellPacketQuery(q);
    if(op=="shell_score")return shellScoreQuery(q);
    return reebQuery(q);
}
}
