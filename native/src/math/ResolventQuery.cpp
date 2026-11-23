#include "ResolventQuerySupport.hpp"
#include <dynsyslab/math/FrequencyGramian.hpp>
#include <dynsyslab/math/GeneralizedLaplace.hpp>
#include <dynsyslab/math/OrbitGramian.hpp>
#include <dynsyslab/math/ResolventQuery.hpp>
#include <dynsyslab/math/ShellResolvent.hpp>
#include <limits>
// Thin JSON adapter. Algorithms and their citations live in LaplaceTransform,
// KoopmanResolvent, GeneralizedLaplace, KoopmanFrequencyResponse, ShellResolvent and
// FrequencyGramian. Returns nullopt for an operation it does not own.
namespace dynsyslab::math {
namespace {
using namespace resolvent_query;
J laplaceQuery(const J& q){
    const auto g=signal(q);LaplaceTransformOptions o;
    if(has(q,"interpolant")){const auto name=text(q.at("interpolant"));
        if(name=="linear")o.interpolant=LaplaceInterpolant::Linear;
        else if(name!="cubic")throw std::invalid_argument("interpolant must be linear or cubic");}
    if(has(q,"growth_rate"))o.growth=GrowthBound{num(q.at("growth_rate")),num(required(q,"growth_constant")),true};
    A rows;GrowthBound growth;
    for(const auto& s:complexes(required(q,"s"))){
        const auto r=laplaceTransform(g,s,o);growth=r.growth;
        rows.emplace_back(O{{"s",encode(s)},{"value",encode(r.value)},{"inside_region",r.insideRegion},{"tail_bound",optionalNumber(r.tailBound)}});
    }
    O out{{"results",rows},{"growth",O{{"rate",growth.rate},{"constant",growth.constant},{"declared",growth.declared}}},
        {"scope","Window integral plus a tail bound from the declared or fitted growth rate. A fitted rate is an estimate, not a proof of the region of convergence"}};
    if(has(q,"z")){
        const double h=num(required(q,"sample_interval"));
        for(std::size_t k=1;k<g.times.size();++k)if(std::abs(g.times[k]-g.times[k-1]-h)>1e-9*h)throw std::invalid_argument("z-transform needs uniform samples at sample_interval");
        A zr;for(const auto& z:complexes(q.at("z")))zr.emplace_back(O{{"z",encode(z)},{"value",encode(zTransform(g.values,z))},{"s_principal",encode(laplaceFromZ(z,h))}});
        out["z_transform"]=zr;
    }
    return out;
}
J inverseQuery(const J& q){
    const auto model=poleModel(required(q,"poles"));
    InverseLaplaceOptions o;
    o.nodes=std::size_t(integerOr(q,"nodes",32,4,400));o.shift=numberOr(q,"shift",0);
    o.dampingParameter=numberOr(q,"damping_parameter",18.4);
    for(const auto& p:model.poles)o.singularities.push_back(p.eigenvalue);
    const std::string method=has(q,"method")?text(q.at("method")):"talbot";
    if(method!="talbot"&&method!="bromwich")throw std::invalid_argument("method must be talbot or bromwich");
    const LaplaceFunction F=[&](C s){return evaluateResolvent(model,s);};
    A rows;bool verified=true;double worst=0;std::string scope;
    for(double t:reals(required(q,"times"))){
        const auto r=method=="talbot"?talbotInverse(F,t,o):bromwichInverse(F,t,o);
        const C exact=evaluateTimeResponse(model,t);worst=std::max(worst,std::abs(r.value-exact));
        verified=verified&&r.scopeVerified;if(!r.scopeVerified||scope.empty())scope=r.scope;
        rows.emplace_back(O{{"t",t},{"value",encode(r.value)},{"model_value",encode(exact)},{"scope_verified",r.scopeVerified}});
    }
    return O{{"method",method},{"results",rows},{"maximum_error_against_model",worst},{"scope_verified",verified},{"scope",scope}};
}
J latticeFrom(const J& spec){
    object(spec);const auto kind=text(required(spec,"kind"));const auto exponents=complexes(required(spec,"exponents"));
    const int order=int(integerOr(spec,"maximum_order",3,0,64)),harmonic=int(integerOr(spec,"maximum_harmonic",3,0,4096));
    const double tolerance=numberOr(spec,"coincidence_tolerance",1e-10);
    if(kind=="equilibrium")return latticeJson(equilibriumPoleLattice(exponents,order,tolerance));
    if(kind=="limit_cycle")return latticeJson(limitCyclePoleLattice(exponents,num(required(spec,"frequency")),order,harmonic,tolerance));
    if(kind=="torus")return latticeJson(torusPoleLattice(exponents,reals(required(spec,"frequencies")),order,harmonic,tolerance));
    throw std::invalid_argument("lattice kind must be equilibrium, limit_cycle or torus");
}
J koopmanResolventQuery(const J& q){
    if(!has(q,"poles")&&!has(q,"lattice"))throw std::invalid_argument("koopman_resolvent needs poles or lattice");
    O out;
    if(has(q,"poles")){
        const auto model=poleModel(q.at("poles"));
        A rows;for(const auto& s:complexes(required(q,"s")))rows.emplace_back(O{{"s",encode(s)},{"value",encode(evaluateResolvent(model,s))}});
        out["results"]=rows;out["abscissa"]=modelAbscissa(model);
    }
    if(has(q,"lattice"))out["lattice"]=latticeFrom(q.at("lattice"));
    return out;
}
J residueQuery(const J& q){
    const auto model=poleModel(required(q,"poles"));CauchyOptions o;
    o.maximumOrder=std::size_t(integerOr(q,"maximum_order",2,1,16));
    const auto e=cauchyLaurent([&](C s){return evaluateResolvent(model,s);},complexValue(required(q,"center")),num(required(q,"radius")),o);
    return O{{"laurent",encodeAll(e.laurent)},{"order",double(e.order)},{"nodes",double(e.nodes)},{"relative_change",e.relativeChange},
        {"higher_order_ratio",e.higherOrderRatio},{"converged",e.converged}};
}
J glaQuery(const J& q){
    const auto g=signal(q);GeneralizedLaplaceOptions o;
    if(has(q,"refine"))o.refine=flag(q,"refine");o.amplificationLimit=numberOr(q,"amplification_limit",1e6);
    const auto r=generalizedLaplaceAnalysis(g,complexes(required(q,"eigenvalues")),o);
    A modes;for(const auto& m:r.modes)modes.emplace_back(O{{"eigenvalue",encode(m.eigenvalue)},{"coefficient",encode(m.coefficient)},{"sequential_coefficient",encode(m.sequentialCoefficient)},{"horizon",m.horizon}});
    O out{{"modes",modes},{"residual_rms",r.residualRms},{"signal_rms",r.signalRms},{"refined",r.refined},{"refinement_condition",std::isfinite(r.refinementCondition)?J(r.refinementCondition):J()}};
    if(has(q,"decomposition")){const auto& d=object(q.at("decomposition"));
        const auto s=decomposeStationary(g,num(required(d,"frequency")),int(integerOr(d,"maximum_harmonic",8,0,4096)),numberOr(d,"transient_fraction",.5));
        A harmonics;for(const auto& [n,c]:s.harmonics)harmonics.emplace_back(A{double(n),c.real(),c.imag()});
        out["decomposition"]=O{{"harmonics",harmonics},{"decay_rate",s.decayRateAvailable?J(s.decayRate):J()},{"fit_residual",s.fitResidual},{"strobe_count",double(s.strobeCount)}};
    }
    return out;
}
J pseudospectrumQuery(const J& q){
    const auto r=pseudospectrumGrid(matrix(required(q,"matrix")),reals(required(q,"real_axis")),reals(required(q,"imag_axis")));
    return O{{"real_axis",numbers(r.realAxis)},{"imag_axis",numbers(r.imagAxis)},{"sigma_min",numbers(r.sigmaMin)},{"eigenvalues",encodeAll(r.eigenvalues)},
        {"layout","sigma_min is row-major: imag_axis rows by real_axis columns"}};
}
J harmonicQuery(const J& q){
    const auto m=bilinearModel(q);const double w=num(required(q,"frequency"));
    HarmonicResponseOptions o;o.maximumOrder=int(integerOr(q,"maximum_order",4,1,64));
    if(has(q,"direction"))o.direction=complexes(q.at("direction"));
    const auto r=harmonicResponse(m,w,o);A outputs;for(const auto& h:r.outputs)outputs.push_back(encodeAll(h));
    A warnings;for(const auto& s:r.warnings)warnings.emplace_back(s);
    O out{{"frequency",w},{"harmonics",outputs},{"sigma_min",numbers(r.sigmaMin)},{"hurwitz",r.hurwitz},{"valid",r.valid},
        {"series_radius",optionalNumber(r.seriesRadius)},{"warnings",warnings}};
    if(has(q,"amplitude")){const double a=num(q.at("amplitude"));const bool real=flag(q,"real_input");
        const auto f=bilinearFloquet(m,w,a,!real,o.direction);
        out["floquet"]=O{{"multipliers",encodeAll(f.multipliers)},{"spectral_radius",f.spectralRadius},{"stable",f.stable}};
        if(real){const auto b=realHarmonicBalance(m,w,a,int(integerOr(q,"maximum_harmonic",8,1,256)),o.direction);
            A rows;for(const auto& [k,y]:b.outputs)rows.emplace_back(O{{"harmonic",double(k)},{"output",encodeAll(y)}});
            out["real_input"]=O{{"outputs",rows},{"truncation_residual",b.truncationResidual}};}
    }
    return out;
}
J volterraQuery(const J& q){
    const auto m=bilinearModel(q);A rows;
    for(const auto& item:array(required(q,"arguments"))){
        const auto s=complexes(item);std::vector<std::size_t> channels;
        if(has(q,"channels"))for(const auto& c:array(q.at("channels")))channels.push_back(std::size_t(integer(c,0,4095)));
        rows.emplace_back(O{{"arguments",encodeAll(s)},{"kernel",encodeAll(volterraKernel(m,s,channels))}});
    }
    return O{{"kernels",rows},{"convention","Symmetric kernel: u=sum a_k e^{s_k t} gives order-n output sum a_k1..a_kn H_n(s_k1..s_kn) e^{(s_k1+..+s_kn)t}"}};
}
ShellResolventTable tableFrom(const J& q){
    ShellResolventTable t;t.energy=reals(required(q,"energies"));t.frequency=reals(required(q,"frequencies"));t.density=reals(required(q,"densities"));
    for(const auto& row:array(required(q,"harmonics"))){object(row);
        const int n=int(integer(required(row,"harmonic"),-100000,100000));
        if(!t.harmonics.emplace(n,complexes(required(row,"values"))).second)throw std::invalid_argument("Duplicate harmonic");}
    validateShellTable(t);
    return t;
}
FourierSeries seriesFrom(const J& q){
    FourierSeries s{num(required(q,"angular_frequency")),{}};
    for(const auto& row:array(required(q,"coefficients"))){const auto& a=array(row);
        if(a.size()!=3)throw std::invalid_argument("Fourier coefficient row needs [harmonic,real,imaginary]");
        if(!s.coefficients.emplace(int(integer(a[0],-1000000,1000000)),C(num(a[1]),num(a[2]))).second)throw std::invalid_argument("Duplicate harmonic");}
    return s;
}
J modelEvaluation(const PoleResidueModel& m,const J& q){
    O out{{"poles",poleModelJson(m)},{"abscissa",modelAbscissa(m)}};
    if(has(q,"s")){A rows;for(const auto& s:complexes(q.at("s")))rows.emplace_back(O{{"s",encode(s)},{"value",encode(evaluateResolvent(m,s))}});out["transform"]=rows;}
    if(has(q,"times")){A rows;for(double t:reals(q.at("times")))rows.emplace_back(O{{"t",t},{"value",encode(evaluateTimeResponse(m,t))}});out["response"]=rows;}
    return out;
}
J shellQuery(const J& q){
    const auto mode=text(required(q,"mode"));
    if(mode=="orbit")return modelEvaluation(orbitHarmonicModel(seriesFrom(object(required(q,"series"))),numberOr(q,"phase",0)),q);
    if(mode=="jordan"){std::vector<JordanMode> modes;
        for(const auto& row:array(required(q,"modes"))){object(row);
            modes.push_back({int(integer(required(row,"harmonic"),-100000,100000)),complexValue(required(row,"amplitude")),complexValue(required(row,"slope"))});}
        return modelEvaluation(jordanOrbitModel(modes,num(required(q,"frequency")),num(required(q,"action_twist"))),q);}
    if(mode=="boundary"){
        ShellBoundaryOptions o;o.slopeTolerance=numberOr(q,"slope_tolerance",1e-10);
        if(has(q,"damping")){o.damping=reals(q.at("damping"));o.dampingScale=num(required(q,"damping_scale"));}
        const auto r=shellBoundaryValue(tableFrom(q),num(required(q,"frequency")),o);
        A roots;for(const auto& x:r.roots)roots.emplace_back(O{{"harmonic",double(x.harmonic)},{"energy",x.energy},{"frequency_slope",x.frequencySlope},
            {"density",x.density},{"weight",x.weight},{"stationary_phase_weight",optionalNumber(x.stationaryPhaseWeight)}});
        A warnings;for(const auto& s:r.warnings)warnings.emplace_back(s);
        return O{{"frequency",r.frequency},{"value",encode(r.value)},{"absorptive",r.absorptive},{"dispersive",r.dispersive},{"roots",roots},
            {"zero_harmonic_mass",r.zeroHarmonicMass},{"valid",r.valid},{"warnings",warnings}};
    }
    if(mode=="damped"){const auto t=tableFrom(q);const auto rate=reals(required(q,"damping_rate"));const double d=num(required(q,"damping_scale"));
        A rows;for(const auto& s:complexes(required(q,"s")))rows.emplace_back(O{{"s",encode(s)},{"value",encode(dampedShellPairing(t,s,rate,d))}});
        return O{{"results",rows}};}
    if(mode=="projector"){
        const RationalCircleKernel kernel(num(required(q,"bandwidth")),std::size_t(integerOr(q,"order",2,1,8)));
        const double theta=num(required(q,"angle"));const auto r=rationalKernelResolvent(kernel,theta);
        A terms;for(const auto& x:r.terms)terms.emplace_back(O{{"node",encode(x.node)},{"coefficient",encode(x.coefficient)}});
        O out{{"constant",encode(r.constant)},{"terms",terms}};
        if(has(q,"poles"))out["projected_poles"]=poleModelJson(rationalProjectedModel(poleModel(q.at("poles")),num(required(q,"sample_interval")),theta,kernel));
        return out;
    }
    throw std::invalid_argument("mode must be orbit, jordan, boundary, damped or projector");
}
J gramianQuery(const J& q){
    const auto Am=matrix(required(q,"A")),Bm=matrix(required(q,"B")),Cm=matrix(required(q,"C"));FrequencyGramianOptions o;
    o.tolerance=numberOr(q,"tolerance",1e-9);
    const auto r=frequencyGramians(Am,Bm,Cm,o);
    O out{{"controllability",encodeMatrix(r.controllability)},{"observability",encodeMatrix(r.observability)},{"h2_norm",r.h2Norm},
        {"relative_change",r.relativeChange},{"evaluations",double(r.evaluations)},{"converged",r.converged}};
    try{const auto p=resolventHInfinity(Am,Bm,Cm);out["hinfinity"]=O{{"value",p.value},{"frequency",p.frequency},{"iterations",double(p.iterations)}};}
    catch(const std::runtime_error&){out["hinfinity"]=J();}
    return out;
}
}
std::optional<core::json::Value> evaluateResolventQuery(const std::string& op,const core::json::Value& q){
    if(op=="laplace_transform")return laplaceQuery(q);
    if(op=="inverse_laplace")return inverseQuery(q);
    if(op=="koopman_resolvent")return koopmanResolventQuery(q);
    if(op=="resolvent_residue")return residueQuery(q);
    if(op=="generalized_laplace")return glaQuery(q);
    if(op=="pseudospectrum")return pseudospectrumQuery(q);
    if(op=="harmonic_response")return harmonicQuery(q);
    if(op=="volterra_kernels")return volterraQuery(q);
    if(op=="shell_resolvent")return shellQuery(q);
    if(op=="frequency_gramian")return gramianQuery(q);
    return std::nullopt;
}
}
