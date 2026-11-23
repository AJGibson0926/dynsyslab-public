#include <dynsyslab/math/OrbitQuery.hpp>
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/math/DampedOrbit.hpp>
#include <dynsyslab/math/MelnikovTools.hpp>
#include <dynsyslab/math/OrbitGramian.hpp>
#include <cmath>
#include <memory>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;using C=std::complex<double>;
using Expression=core::AnalyticExpression;
constexpr std::size_t kMaximumEnergies=4096,kMaximumSamples=65536;
double num(const J& v){
    if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");
    return v.number();
}
double numOr(const J& q,const char* key,double fallback){return q.at(key).isNull()?fallback:num(q.at(key));}
int integer(const J& v,int lo=0,int hi=1000000){
    const double n=num(v);if(n<lo||n>hi||std::floor(n)!=n)throw std::invalid_argument("Integer out of range");return int(n);
}
bool flag(const J& q,const char* key){
    const auto& v=q.at(key);if(v.isNull())return false;
    if(!std::holds_alternative<bool>(v.storage()))throw std::invalid_argument("Boolean required");return v.boolean();
}
const A& array(const J& v){
    if(!std::holds_alternative<A>(v.storage())||v.array().size()>1000000)throw std::invalid_argument("Bounded array required");
    return v.array();
}
C complexValue(const J& v){
    if(std::holds_alternative<double>(v.storage()))return num(v);
    const auto& a=array(v);if(a.size()!=2)throw std::invalid_argument("Complex number needs [real,imag]");return {num(a[0]),num(a[1])};
}
std::vector<double> reals(const J& v){std::vector<double> out;out.reserve(array(v).size());for(const auto& x:array(v))out.push_back(num(x));return out;}
ComplexVector complexes(const J& v){ComplexVector out;out.reserve(array(v).size());for(const auto& x:array(v))out.push_back(complexValue(x));return out;}
std::vector<ComplexVector> channels(const J& v){std::vector<ComplexVector> out;for(const auto& row:array(v))out.push_back(complexes(row));return out;}
PlanarPoint point(const J& v){const auto x=reals(v);if(x.size()!=2)throw std::invalid_argument("Planar point needs two values");return {x[0],x[1]};}
J finiteOrNull(double v){return std::isfinite(v)?J(v):J();}
A encode(const std::vector<double>& v){A out;out.reserve(v.size());for(double x:v)out.push_back(finiteOrNull(x));return out;}
J encode(C z){return A{z.real(),z.imag()};}
A encode(const ComplexVector& v){A out;out.reserve(v.size());for(auto z:v)out.push_back(encode(z));return out;}
J encode(const PlanarPoint& p){return A{p[0],p[1]};}
J encode(const std::array<double,3>& s){return A{s[0],s[1],s[2]};}
A encode(const std::map<int,C>& m){A out;for(const auto& [n,z]:m)out.emplace_back(O{{"harmonic",double(n)},{"coefficient",encode(z)}});return out;}
J encode(const FourierSeries& s){A rows;for(const auto& [n,c]:s.coefficients)rows.emplace_back(A{double(n),c.real(),c.imag()});
    return O{{"angular_frequency",s.angularFrequency},{"coefficients",rows}};}
Expression planarExpression(const J& v,const char* name){
    if(v.isNull())throw std::invalid_argument(std::string(name)+" expression is required");
    auto e=Expression::fromJson(v);if(e.stateDimension()!=2)throw std::invalid_argument(std::string(name)+" must be a planar expression");return e;
}
PlanarScalar scalar(const Expression& e){return [e](const PlanarPoint& x){return e.evaluate(x.data());};}
// Canonical Hamiltonian field in the declared (coordinate,momentum) order.
PlanarField hamiltonianField(const Expression& H,int coordinate,int momentum){
    const auto dq=H.derivative(std::size_t(coordinate)),dp=H.derivative(std::size_t(momentum));
    return [dq,dp,coordinate,momentum](const PlanarPoint& x){PlanarPoint f{};f[coordinate]=dp.evaluate(x.data());f[momentum]=-dq.evaluate(x.data());return f;};
}
std::array<int,2> canonicalPair(const J& q){
    if(q.at("canonical_pair").isNull())return {0,1};
    const auto& a=array(q.at("canonical_pair"));if(a.size()!=2)throw std::invalid_argument("canonical_pair needs two indices");
    const std::array<int,2> pair{integer(a[0],0,1),integer(a[1],0,1)};
    if(pair[0]==pair[1])throw std::invalid_argument("canonical_pair indices must differ");return pair;
}
PlanarOrbitSystem declaredSystem(const J& q){
    const auto H=planarExpression(q.at("hamiltonian"),"Hamiltonian");const auto pair=canonicalPair(q);
    PlanarOrbitSystem s;s.coordinate=pair[0];s.momentum=pair[1];s.hamiltonian=scalar(H);s.drift=hamiltonianField(H,pair[0],pair[1]);
    s.centre=point(q.at("centre"));s.sectionDirection=q.at("section_direction").isNull()?PlanarPoint{1,0}:point(q.at("section_direction"));
    s.maximumRadius=num(q.at("maximum_radius"));
    if(!q.at("density").isNull())s.density=scalar(planarExpression(q.at("density"),"Density"));
    return s;
}
OrbitFamilyOptions familyOptions(const J& q){
    OrbitFamilyOptions o;o.flow.relativeTolerance=numOr(q,"relative_tolerance",o.flow.relativeTolerance);
    o.flow.absoluteTolerance=numOr(q,"absolute_tolerance",o.flow.absoluteTolerance);
    if(!(o.flow.relativeTolerance>=1e-14)||!(o.flow.absoluteTolerance>0))throw std::invalid_argument("Orbit tolerances out of range");
    if(!q.at("phase_samples").isNull())o.phaseSamples=std::size_t(integer(q.at("phase_samples"),8,int(kMaximumSamples)));
    if(!q.at("maximum_harmonic").isNull())o.maximumHarmonic=integer(q.at("maximum_harmonic"),0,int(kMaximumSamples));
    if(std::size_t(o.maximumHarmonic)*2>=o.phaseSamples)throw std::invalid_argument("Need 2*maximum_harmonic < phase_samples");
    return o;
}
std::vector<double> energyList(const J& q){
    auto e=reals(q.at("energies"));if(e.size()<2||e.size()>kMaximumEnergies)throw std::invalid_argument("Need 2..4096 energies");return e;
}
J branchesJson(const FrequencyBranchReport& r){
    A rows;for(const auto& b:r.branches)rows.emplace_back(O{{"first",double(b.first)},{"last",double(b.last)},
        {"direction",double(b.direction)},{"lower",b.lower},{"upper",b.upper}});
    return O{{"monotone",r.monotone},{"flat_segments",double(r.flatSegments)},{"branches",rows}};
}
A ladderJson(const std::vector<ShellRoot>& roots){
    A rows;for(const auto& s:roots)rows.emplace_back(O{{"harmonic",double(s.harmonic)},{"energy",s.root.energy},
        {"frequency",s.root.frequency},{"frequency_slope",s.root.slope},{"branch",double(s.root.branch)},{"refined",s.root.refined}});
    return rows;
}
J orbitFamilyQuery(const J& q){
    O out;std::vector<double> energies,frequencies,periods;FrequencyEvaluator evaluator;std::shared_ptr<PlanarOrbitSystem> system;OrbitFamilyOptions options;
    if(!q.at("table").isNull()){
        energies=reals(q.at("table").at("energies"));frequencies=reals(q.at("table").at("frequencies"));
        if(energies.size()<2||energies.size()>kMaximumEnergies)throw std::invalid_argument("Need 2..4096 table energies");
        if(!q.at("table").at("periods").isNull())periods=reals(q.at("table").at("periods"));
        out["source"]="sampled_table";
    }else{
        system=std::make_shared<PlanarOrbitSystem>(declaredSystem(q));options=familyOptions(q);energies=energyList(q);
        const auto family=orbitFamily(*system,energies,options);frequencies=family.frequencies;periods=family.periods;
        A rows;double worstClosure=0,worstDrift=0;
        for(std::size_t i=0;i<family.orbits.size();++i){const auto& r=family.orbits[i];worstClosure=std::max(worstClosure,r.closureError);worstDrift=std::max(worstDrift,r.energyDrift);
            O row{{"energy",r.energy},{"period",r.period},{"frequency",r.frequency},{"action",r.action},{"frequency_slope",family.frequencySlopes[i]},
                {"orientation",double(r.orientation)},{"start",encode(r.start)},{"closure_error",r.closureError},{"energy_drift",r.energyDrift}};
            if(flag(q,"include_samples")||flag(q,"include_fourier")){
                const auto samples=orbitSamples(*system,r,options);
                if(flag(q,"include_samples")){A points;for(const auto& p:samples)points.push_back(encode(p));row["samples"]=points;}
                if(flag(q,"include_fourier")){const auto f=orbitFourier(samples,r.frequency,options.maximumHarmonic);row["fourier"]=A{encode(f.first),encode(f.second)};}
            }
            rows.emplace_back(row);}
        out["source"]="declared_hamiltonian";out["orbits"]=rows;out["maximum_closure_error"]=worstClosure;out["maximum_energy_drift"]=worstDrift;
        evaluator=[system,options](double e){return periodicOrbit(*system,e,options).frequency;};
    }
    if(energies.size()!=frequencies.size())throw std::invalid_argument("Table sizes differ");
    out["frequency_branches"]=branchesJson(frequencyBranches(energies,frequencies));
    if(!q.at("target_frequency").isNull()){A rows;for(const auto& r:invertFrequency(energies,frequencies,num(q.at("target_frequency")),evaluator))
        rows.emplace_back(O{{"energy",r.energy},{"frequency",r.frequency},{"frequency_slope",r.slope},{"branch",double(r.branch)},{"refined",r.refined}});out["inverse"]=rows;}
    if(!q.at("drive_frequency").isNull())out["shell_ladder"]=ladderJson(shellLadder(energies,frequencies,num(q.at("drive_frequency")),
        q.at("maximum_shell").isNull()?8:integer(q.at("maximum_shell"),1,10000),evaluator));
    if(!q.at("separatrix_energy").isNull()){
        if(periods.size()!=energies.size())throw std::invalid_argument("Separatrix fit needs periods");
        const auto count=q.at("separatrix_points").isNull()?energies.size():std::size_t(integer(q.at("separatrix_points"),3,int(energies.size())));
        const std::vector<double> e(energies.end()-std::ptrdiff_t(count),energies.end()),T(periods.end()-std::ptrdiff_t(count),periods.end());
        const auto fit=fitSeparatrixPeriod(e,T,num(q.at("separatrix_energy")));
        out["separatrix_fit"]=O{{"lambda",fit.lambda},{"constant",fit.constant},{"rms_residual",fit.rmsResidual},{"samples",double(fit.samples)}};
    }
    out["scope"]="Numerical period annulus over the supplied energies; no enclosure or global integrability claim";
    return out;
}
J dampedAveragesQuery(const J& q){
    O out;DampedTable table;const double d=num(q.at("damping_scale"));
    const auto eq=num(q.at("equilibrium_frequency")),lambda=num(q.at("real_eigenvalue"));
    if(!q.at("table").isNull()){
        table={reals(q.at("table").at("energies")),reals(q.at("table").at("frequencies")),reals(q.at("table").at("damping"))};
        if(table.energy.size()>kMaximumEnergies)throw std::invalid_argument("Damped table is too large");
        if(!q.at("q_samples").isNull()){const auto r=dampedAverage(reals(q.at("q_samples")),num(q.at("frequency")));
            out["mean"]=r.mean;out["chi"]=encode(r.chi);}
    }else{
        const auto system=declaredSystem(q);const auto options=familyOptions(q);const auto Q=scalar(planarExpression(q.at("dissipation"),"Dissipation"));
        table.energy=energyList(q);A rows;
        for(double e:table.energy){const auto orbit=periodicOrbit(system,e,options);std::vector<double> samples;
            for(const auto& p:orbitSamples(system,orbit,options))samples.push_back(Q(p));
            const auto r=dampedAverage(samples,orbit.frequency);table.frequency.push_back(orbit.frequency);table.damping.push_back(r.mean);
            O row{{"energy",e},{"frequency",orbit.frequency},{"damping",r.mean}};if(flag(q,"include_chi"))row["chi"]=encode(r.chi);rows.emplace_back(row);}
        out["orbits"]=rows;
    }
    const auto phase=dampedPhase(table,eq,lambda,d);
    out["energies"]=encode(phase.energy);out["phase_function"]=encode(phase.phaseFunction);out["phase_slope"]=encode(phase.phaseSlope);out["amplitude"]=encode(phase.amplitude);
    out["scope"]="First-order averaging over the supplied energies; S uses dE; amplitude is zero at the last energy; no remainder bound";
    return out;
}
J melnikovQuery(const J& q){
    const auto H=planarExpression(q.at("hamiltonian"),"Hamiltonian");const auto pair=canonicalPair(q);
    ConnectionProblem p;p.drift=hamiltonianField(H,pair[0],pair[1]);p.midpoint=point(q.at("midpoint"));
    p.forwardSaddle=point(q.at("forward_saddle"));p.backwardSaddle=point(q.at("backward_saddle"));
    p.forwardRate=num(q.at("forward_rate"));p.backwardRate=num(q.at("backward_rate"));
    p.tailRadius=numOr(q,"tail_radius",p.tailRadius);p.maximumTime=numOr(q,"maximum_time",p.maximumTime);
    p.flow.relativeTolerance=numOr(q,"relative_tolerance",p.flow.relativeTolerance);p.flow.absoluteTolerance=numOr(q,"absolute_tolerance",p.flow.absoluteTolerance);
    const auto encodeIntegral=[](const ConnectionIntegral& r){return O{{"value",encode(r.value)},{"core",encode(r.core)},{"forward_tail",encode(r.forwardTail)},
        {"backward_tail",encode(r.backwardTail)},{"forward_time",r.forwardTime},{"backward_time",r.backwardTime},
        {"forward_decay",r.forwardDecay},{"backward_decay",r.backwardDecay}};};
    O out;double damping=0;bool haveDamping=false;
    if(!q.at("dissipation").isNull()){const auto r=dampingIntegral(p,scalar(planarExpression(q.at("dissipation"),"Dissipation")));
        damping=r.value.real();haveDamping=true;out["damping_integral"]=encodeIntegral(r);}
    if(!q.at("forcing").isNull()){
        const auto& pairY=array(q.at("forcing"));if(pairY.size()!=2)throw std::invalid_argument("forcing needs two component expressions");
        const auto y0=planarExpression(pairY[0],"Forcing"),y1=planarExpression(pairY[1],"Forcing");const double w=num(q.at("forcing_frequency"));
        const PlanarField Y=[y0,y1](const PlanarPoint& x){return PlanarPoint{y0.evaluate(x.data()),y1.evaluate(x.data())};};
        const auto r=forcingIntegral(p,Y,w);const double amplitude=std::abs(r.value);out["forcing_integral"]=encodeIntegral(r);out["forcing_amplitude"]=amplitude;
        if(haveDamping&&amplitude>0){out["threshold_ratio"]=melnikovThreshold(damping,amplitude);
            if(!q.at("damping_scale").isNull()&&!q.at("forcing_scale").isNull()){const auto f=lobeFlux(num(q.at("damping_scale")),damping,num(q.at("forcing_scale")),amplitude,w);
                out["lobe"]=O{{"ratio",f.ratio},{"area",f.area},{"flux_per_time",f.fluxPerTime},{"transverse",f.transverse}};}}
    }
    if(out.empty())throw std::invalid_argument("melnikov_homoclinic needs dissipation or forcing");
    out["tails_rigorous"]=false;out["scope"]="First-order Melnikov on the supplied connection; exponential tail estimate; simple-zero prediction only";
    return out;
}
J modalQuery(const J& q){
    std::optional<double> horizon;if(!q.at("horizon").isNull())horizon=num(q.at("horizon"));
    const auto r=modalGramians(complexes(q.at("eigenvalues")),channels(q.at("inputs")),channels(q.at("outputs")),horizon);
    return O{{"dimension",double(r.dimension)},{"controllability",encode(r.controllability)},{"observability",encode(r.observability)},
        {"h2_norm",r.h2Norm},{"hankel_singular_values",encode(r.hankelSingularValues)},
        {"scope",horizon?"Finite-horizon modal Gramians":"Infinite-horizon modal Gramians"}};
}
J orbitGramianQuery(const J& q){
    if(!q.at("stiffness").isNull()){const auto s=oscillatorHankel(num(q.at("stiffness")),num(q.at("damping")));return O{{"hankel_singular_values",encode(s)}};}
    if(flag(q,"collocated")){const auto s=collocatedOrbitHankel(num(q.at("input_variance")),num(q.at("action_twist")),num(q.at("horizon")));
        return O{{"hankel_singular_values",encode(s)},{"scope","Collocated Jordan shear law"}};}
    OrbitGramianInput g;g.inputVariance=num(q.at("input_variance"));g.outputVariance=numOr(q,"output_variance",0);g.frequencySlope=num(q.at("frequency_slope"));
    g.horizon=num(q.at("horizon"));g.inputCross=numOr(q,"input_cross",0);g.inputPhase=numOr(q,"input_phase",0);g.outputCross=numOr(q,"output_cross",0);g.outputEnergy=numOr(q,"output_energy",0);
    const auto r=undampedOrbitGramians(g);
    return O{{"controllability",encode(r.controllability)},{"observability",encode(r.observability)},{"hankel_singular_values",encode(r.hankel)},
        {"entry_order","[EE,Ephi,phiphi]"}};
}
J dampedGramianQuery(const J& q){
    DampedGramianInput in;in.energy=reals(q.at("energies"));in.inputVariance=reals(q.at("input_variance"));in.outputVariance=reals(q.at("output_variance"));
    in.damping=reals(q.at("damping"));in.frequency=reals(q.at("frequencies"));in.dampingScale=num(q.at("damping_scale"));
    in.initialEnergy=num(q.at("initial_energy"));in.targetEnergy=num(q.at("target_energy"));in.finalEnergy=num(q.at("final_energy"));
    if(in.energy.size()>kMaximumEnergies)throw std::invalid_argument("Damped Gramian table is too large");
    const auto r=dampedOrbitGramians(in);
    return O{{"c_EE",r.cEE},{"c_Ephi",r.cEphi},{"c_phiphi",r.cPhiphi},{"o_phiphi",r.oPhiphi},{"o_Ephi",r.oEphi},{"o_EE",r.oEE},
        {"target_damping",r.targetDamping},{"target_frequency",r.targetFrequency},{"hankel_singular_values",encode(r.hankel)},
        {"scope","Leading slow-passage quadratures on the supplied table"}};
}
J jordanQuery(const J& q){
    std::vector<JordanMode> modes;
    for(const auto& m:array(q.at("modes")))modes.push_back({integer(m.at("harmonic"),-100000,100000),complexValue(m.at("amplitude")),complexValue(m.at("slope"))});
    if(flag(q,"real_observable")){const auto base=modes;
        for(const auto& m:base){if(m.harmonic==0)continue;bool found=false;for(const auto& k:base)found=found||k.harmonic==-m.harmonic;
            if(!found)modes.push_back({-m.harmonic,std::conj(m.amplitude),std::conj(m.slope)});}}
    O out{{"frozen_residues",encode(frozenKromResidues(modes))}};
    if(!q.at("times").isNull())out["response"]=encode(jordanImpulseResponse(modes,num(q.at("frequency")),num(q.at("action_twist")),reals(q.at("times"))));
    if(!q.at("lift_harmonic").isNull())out["input_lift"]=encode(exactInputLift(integer(q.at("lift_harmonic"),-100000,100000),modes));
    out["scope"]="Linear response at one action; no action diffusion or nonlinear remainder";
    return out;
}
}
std::optional<core::json::Value> evaluateOrbitQuery(const std::string& op,const core::json::Value& q){
    if(op=="orbit_family")return orbitFamilyQuery(q);
    if(op=="damped_averages")return dampedAveragesQuery(q);
    if(op=="melnikov_homoclinic")return melnikovQuery(q);
    if(op=="capture_probabilities"){
        const auto M=reals(q.at("lobe_integrals"));
        return J(O{{"probabilities",encode(captureProbabilities(M))},{"scope","Adiabatic ensemble ratio; not a single-trajectory prediction"}});
    }
    if(op=="modal_gramian")return modalQuery(q);
    if(op=="orbit_gramian")return orbitGramianQuery(q);
    if(op=="damped_orbit_gramian")return dampedGramianQuery(q);
    if(op=="jordan_response")return jordanQuery(q);
    return std::nullopt;
}
}
