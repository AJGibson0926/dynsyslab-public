#include "ResolventQuerySupport.hpp"
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/math/ContractionCheck.hpp>
#include <dynsyslab/math/DiscreteResolvent.hpp>
#include <dynsyslab/math/ForcedResponseEstimation.hpp>
#include <dynsyslab/math/HarmonicContinuation.hpp>
#include <dynsyslab/math/IsochronDecomposition.hpp>
#include <dynsyslab/math/PronyIdentification.hpp>
#include <dynsyslab/math/ResolventDataQuery.hpp>
#include <limits>
#include <memory>
// Thin JSON adapter. Algorithms and citations live in DiscreteResolvent, PronyIdentification,
// ForcedResponseEstimation, IsochronDecomposition, HarmonicContinuation and ContractionCheck.
namespace dynsyslab::math {
namespace {
using namespace resolvent_query;
using Expression=core::AnalyticExpression;
J stringList(const std::vector<std::string>& v){A out;for(const auto& s:v)out.emplace_back(s);return out;}
J finiteOrNull(double x){return std::isfinite(x)?J(x):J();}
J harmonicRows(const std::map<int,C>& m){A rows;for(const auto& [n,c]:m)rows.emplace_back(A{double(n),c.real(),c.imag()});return rows;}
DiscretePoleModel discreteModel(const J& q){
    DiscretePoleModel m;m.sampleInterval=numberOr(q,"sample_interval",1);
    for(const auto& p:array(required(q,"discrete_poles"))){object(p);m.poles.push_back({complexValue(required(p,"multiplier")),complexes(required(p,"laurent"))});}
    validateDiscreteModel(m);return m;
}
J discreteJson(const DiscretePoleModel& m){
    A rows;for(const auto& p:m.poles)rows.emplace_back(O{{"multiplier",encode(p.multiplier)},{"laurent",encodeAll(p.laurent)}});
    return O{{"sample_interval",m.sampleInterval},{"poles",rows}};
}
J discreteQuery(const J& q){
    if(has(q,"poles")==has(q,"discrete_poles"))throw std::invalid_argument("discrete_resolvent needs exactly one of poles or discrete_poles");
    O out;DiscretePoleModel model;
    if(has(q,"poles")){
        const auto c=poleModel(q.at("poles"));const double h=num(required(q,"sample_interval"));model=discretizePoleModel(c,h);
        const auto r=aliasReport(c,h);A entries;
        for(const auto& e:r.entries)entries.emplace_back(O{{"index",double(e.index)},{"eigenvalue",encode(e.eigenvalue)},{"multiplier",encode(e.multiplier)},
            {"principal",encode(e.principal)},{"turns",double(e.turns)},{"aliased",e.aliased},{"on_nyquist",e.onNyquist}});
        A collisions;for(const auto& x:r.collisions)collisions.emplace_back(O{{"first",double(x.first)},{"second",double(x.second)},{"separation",x.separation}});
        out["aliasing"]=O{{"nyquist",r.nyquist},{"entries",entries},{"collisions",collisions}};
    }else{
        model=discreteModel(q);std::vector<int> branches;
        if(has(q,"branches"))for(const auto& b:array(q.at("branches")))branches.push_back(int(integer(b,-1000000,1000000)));
        out["continuous_poles"]=poleModelJson(continuousPoleModel(model,branches));
        out["branch_note"]="Sampling fixes lambda only modulo 2*pi*i/h; branches are the caller's assumption";
    }
    out["model"]=discreteJson(model);
    if(has(q,"z")){A rows;for(const auto& z:complexes(q.at("z")))rows.emplace_back(O{{"z",encode(z)},{"value",encode(evaluateDiscreteResolvent(model,z))}});out["transform"]=rows;}
    if(has(q,"indices")){A rows;for(const auto& k:array(q.at("indices"))){const auto i=std::size_t(integer(k,0,1000000));rows.emplace_back(O{{"k",double(i)},{"value",encode(evaluateDiscreteSequence(model,i))}});}out["sequence"]=rows;}
    return out;
}
J pronyQuery(const J& q){
    PronyOptions o;o.sampleInterval=numberOr(q,"sample_interval",1);o.rankTolerance=numberOr(q,"rank_tolerance",1e-9);
    if(has(q,"order"))o.order=std::size_t(integer(q.at("order"),1,512));
    o.pencilParameter=std::size_t(integerOr(q,"pencil_parameter",0,0,512));
    const auto r=pronyIdentify(complexes(required(q,"values")),o);A poles,scan;
    for(const auto& p:r.poles)poles.emplace_back(O{{"multiplier",encode(p.multiplier)},{"residue",encode(p.residue)},{"eigenvalue",encode(p.eigenvalue)}});
    for(const auto& s:r.orderScan)scan.emplace_back(O{{"order",double(s.order)},{"relative_residual",s.relativeResidual}});
    return O{{"poles",poles},{"order",double(r.order)},{"order_source",r.orderSource},{"pencil_parameter",double(r.pencilParameter)},{"singular_values",numbers(r.singularValues)},
        {"singular_gap",finiteOrNull(r.singularGap)},{"relative_residual",r.relativeResidual},{"order_scan",scan},{"minimum_separation",finiteOrNull(r.minimumSeparation)},
        {"rayleigh_limit",r.rayleighLimit},{"vandermonde_condition",finiteOrNull(r.vandermondeCondition)},{"clustered",r.clustered},{"near_continuous",r.nearContinuous},{"warnings",stringList(r.warnings)}};
}
J forcedQuery(const J& q){
    O out;const double w=num(required(q,"frequency"));
    if(has(q,"values")){
        const auto r=estimateForcedHarmonics(signal(q),w,num(required(q,"amplitude")),int(integerOr(q,"maximum_harmonic",4,1,4096)),!has(q,"complex_input")||flag(q,"complex_input"),numberOr(q,"transient_fraction",.5));
        out["coefficients"]=harmonicRows(r.coefficients);out["response"]=harmonicRows(r.response);out["tail_relative_residual"]=r.tailRelativeResidual;out["warnings"]=stringList(r.warnings);
    }
    if(has(q,"dmd")){const auto& d=object(q.at("dmd"));SkewDmdInput in;in.frequency=w;in.sampleInterval=num(required(d,"sample_interval"));
        for(const auto& row:array(required(d,"states")))in.states.push_back(complexes(row));
        in.inputs=complexes(required(d,"inputs"));in.outputState=std::size_t(integerOr(d,"output_state",0,0,4096));
        in.maximumHarmonic=int(integerOr(d,"maximum_harmonic",2,1,16));in.maximumDegree=int(integerOr(d,"maximum_degree",2,1,8));
        if(has(d,"dictionary"))for(const auto& row:array(d.at("dictionary"))){std::vector<int> e;for(const auto& x:array(row))e.push_back(int(integer(x,0,16)));in.dictionary.push_back(std::move(e));}
        const auto r=skewProductDmd(in);A rows;
        for(const auto& h:r.harmonics)rows.emplace_back(O{{"harmonic",double(h.harmonic)},{"response",h.response?encode(*h.response):J()},{"mode_residual",h.modeResidual},{"spectrum_distance",finiteOrNull(h.spectrumDistance)}});
        out["dmd"]=O{{"feasible",r.feasible},{"dictionary_size",double(r.dictionarySize)},{"rank",double(r.rank)},{"snapshot_pairs",double(r.snapshots)},{"closure_residual",r.closureResidual},
            {"condition",finiteOrNull(r.condition)},{"harmonics",rows},{"warnings",stringList(r.warnings)}};
    }
    if(out.empty())throw std::invalid_argument("forced_harmonic_estimate needs values or dmd");
    return out;
}
J subharmonicQuery(const J& q){
    const auto r=subharmonicResponse(signal(q),num(required(q,"frequency")),num(required(q,"amplitude")),int(integer(required(q,"cover"),2,4096)),
        int(integerOr(q,"maximum_index",8,1,4096)),numberOr(q,"transient_fraction",.5),numberOr(q,"detection_tolerance",1e-6));
    A labels;for(const auto& l:r.labels)labels.emplace_back(O{{"index",double(l.index)},{"frequency_ratio",l.frequencyRatio},{"minimum_depth",double(l.minimumDepth)}});
    return O{{"cover",double(r.cover)},{"coefficients",harmonicRows(r.coefficients)},{"subharmonic_fraction",r.subharmonicFraction},{"detected",r.detected},{"cover_labelled",r.coverLabelled},
        {"magnitude",optionalNumber(r.magnitude)},{"phase_modulo",optionalNumber(r.phaseModulo)},{"invariant_power",r.invariantPower?encode(*r.invariantPower):J()},{"branches",encodeAll(r.branches)},{"labels",labels},{"warnings",stringList(r.warnings)}};
}
FourierSeries cycleFrom(const J& q){
    FourierSeries s{num(required(q,"angular_frequency")),{}};
    for(const auto& row:array(required(q,"coefficients"))){const auto& a=array(row);
        if(a.size()!=3)throw std::invalid_argument("Fourier coefficient row needs [harmonic,real,imaginary]");
        if(!s.coefficients.emplace(int(integer(a[0],-1000000,1000000)),C(num(a[1]),num(a[2]))).second)throw std::invalid_argument("Duplicate harmonic");}
    return s;
}
J isochronQuery(const J& q){
    const auto g=signal(q);const auto cycle=cycleFrom(object(required(q,"cycle")));IsochronSplitOptions o;
    o.transientFraction=numberOr(q,"transient_fraction",.6);o.phaseTolerance=numberOr(q,"phase_tolerance",1e-3);
    if(has(q,"asymptotic_phase")&&has(q,"phase_map"))throw std::invalid_argument("Give asymptotic_phase or phase_map, not both");
    if(has(q,"asymptotic_phase"))o.asymptoticPhase=num(q.at("asymptotic_phase"));
    IsochronSplit r;
    if(has(q,"phase_map")){const auto x0=reals(required(q,"initial_state"));const auto e=Expression::fromJson(q.at("phase_map"));
        if(e.stateDimension()!=x0.size())throw std::invalid_argument("phase_map dimension must match initial_state");
        r=isochronSplitFromMap(g,cycle,[e](const std::vector<double>& x){return e.evaluate(x.data());},x0,o);
    }else r=isochronSplit(g,cycle,o);
    O out{{"asymptotic_phase",r.asymptoticPhase},{"phase_source",r.phaseSource},{"tail_phase",optionalNumber(r.tailPhase)},{"phase_mismatch",optionalNumber(r.phaseMismatch)},
        {"decay_rate",optionalNumber(r.decayRate)},{"fit_residual",r.fitResidual},{"strobe_count",double(r.strobeCount)},{"warnings",stringList(r.warnings)}};
    if(flag(q,"include_series")){out["stationary"]=encodeAll(r.stationary);out["transient"]=encodeAll(r.transient);}
    return out;
}
// Vector field in (x...,u,p): polynomial terms or analytic expressions.
struct FieldSpec {
    std::size_t states{};
    std::function<void(const double*,double*)> value;    // Argument (x...,u,p).
    std::function<void(const double*,double*)> jacobian; // dF/dx, row-major.
};
FieldSpec fieldFrom(const J& q){
    FieldSpec out;
    if(has(q,"field")==has(q,"expressions"))throw std::invalid_argument("Give exactly one of field or expressions");
    if(has(q,"field")){
        struct Term {double c;std::vector<int> e;};auto eqs=std::make_shared<std::vector<std::vector<Term>>>();
        const auto& rows=array(q.at("field"));out.states=rows.size();if(rows.empty()||rows.size()>32)throw std::invalid_argument("Field needs 1..32 equations");
        for(const auto& row:rows){std::vector<Term> terms;
            for(const auto& t:array(row)){object(t);Term x{num(required(t,"coefficient")),{}};
                for(const auto& p:array(required(t,"powers")))x.e.push_back(int(integer(p,0,32)));
                if(x.e.size()!=rows.size()+2)throw std::invalid_argument("Term powers need one entry per state plus u and p");
                terms.push_back(std::move(x));}
            eqs->push_back(std::move(terms));}
        const auto n=out.states;
        const auto monomial=[n](const Term& t,const double* a,std::size_t skip){double v=t.c;
            for(std::size_t i=0;i<n+2;++i){int e=t.e[i];if(i==skip){if(e==0)return 0.;v*=e;--e;}for(int k=0;k<e;++k)v*=a[i];}return v;};
        out.value=[eqs,n,monomial](const double* a,double* f){for(std::size_t i=0;i<n;++i){f[i]=0;for(const auto& t:(*eqs)[i])f[i]+=monomial(t,a,n+2);}};
        out.jacobian=[eqs,n,monomial](const double* a,double* J){for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){J[i*n+j]=0;for(const auto& t:(*eqs)[i])J[i*n+j]+=monomial(t,a,j);}};
    }else{
        auto f=std::make_shared<std::vector<Expression>>();const auto& rows=array(q.at("expressions"));out.states=rows.size();
        if(rows.empty()||rows.size()>32)throw std::invalid_argument("Field needs 1..32 equations");
        for(const auto& e:rows){f->push_back(Expression::fromJson(e));if(f->back().stateDimension()!=rows.size()+2)throw std::invalid_argument("Expressions need dimension states+2 (x...,u,p)");}
        auto d=std::make_shared<std::vector<Expression>>();const auto n=out.states;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)d->push_back((*f)[i].derivative(j));
        out.value=[f,n](const double* a,double* y){for(std::size_t i=0;i<n;++i)y[i]=(*f)[i].evaluate(a);};
        out.jacobian=[d,n](const double* a,double* J){for(std::size_t k=0;k<n*n;++k)J[k]=(*d)[k].evaluate(a);};
    }
    return out;
}
J continuationQuery(const J& q){
    const auto spec=fieldFrom(q);ContinuationProblem p;p.dimension=spec.states;
    p.field=[spec,n=spec.states](const double* x,double u,double par,double* dx){std::vector<double> a(x,x+n);a.push_back(u);a.push_back(par);spec.value(a.data(),dx);};
    p.frequency=num(required(q,"frequency"));p.amplitude=num(required(q,"amplitude"));p.parameter=numberOr(q,"parameter",0);
    const auto variable=text(required(q,"variable"));
    if(variable=="frequency")p.variable=ContinuationVariable::Frequency;else if(variable=="amplitude")p.variable=ContinuationVariable::Amplitude;
    else if(variable=="parameter")p.variable=ContinuationVariable::Parameter;else throw std::invalid_argument("variable must be frequency, amplitude or parameter");
    p.lower=num(required(q,"lower"));p.upper=num(required(q,"upper"));p.harmonics=int(integerOr(q,"harmonics",5,1,64));
    p.step=numberOr(q,"step",.01);p.minimumStep=numberOr(q,"minimum_step",1e-6);p.maximumStep=numberOr(q,"maximum_step",.05);
    p.maximumSteps=std::size_t(integerOr(q,"maximum_steps",1000,1,100000));p.direction=int(integerOr(q,"direction",1,-1,1));
    p.outputState=std::size_t(integerOr(q,"output_state",0,0,31));p.tolerance=numberOr(q,"tolerance",1e-9);
    if(has(q,"compute_stability"))p.computeStability=flag(q,"compute_stability");
    if(has(q,"initial_coefficients"))p.initialCoefficients=reals(q.at("initial_coefficients"));
    const auto r=continueHarmonicBalance(p);A points,folds;
    for(const auto& x:r.points){O row{{"value",x.value},{"residual",x.residual},{"response",harmonicRows(x.response)},{"floquet_radius",optionalNumber(x.floquetRadius)}};
        row["stable"]=x.stable?J(*x.stable):J();points.emplace_back(std::move(row));}
    for(const auto& f:r.folds)folds.emplace_back(O{{"value",f.value},{"after_point",double(f.afterPoint)},{"response_magnitude",f.responseMagnitude}});
    return O{{"points",points},{"folds",folds},{"termination",r.termination},{"branch_lost",r.branchLost},{"warnings",stringList(r.warnings)},
        {"scope","Response coefficients under real forcing a*cos(wt); they equal H_n scaled by (a/2)^n only in the small-amplitude limit"}};
}
J contractionQuery(const J& q){
    std::vector<double> weights;if(has(q,"weights"))weights=reals(q.at("weights"));
    const auto normsJson=[](const LogNorms& n){return O{{"two",n.two},{"one",n.one},{"infinity",n.infinity}};};
    if(has(q,"interval_lower")){
        const auto lo=matrix(q.at("interval_lower")),hi=matrix(required(q,"interval_upper"));
        if(lo.rows!=lo.cols||hi.rows!=lo.rows||hi.cols!=lo.cols)throw std::invalid_argument("Interval bounds must be square and equal in size");
        std::vector<double> l,u;for(const auto& z:lo.data){if(z.imag()!=0)throw std::invalid_argument("Interval bounds must be real");l.push_back(z.real());}
        for(const auto& z:hi.data){if(z.imag()!=0)throw std::invalid_argument("Interval bounds must be real");u.push_back(z.real());}
        const auto r=intervalContraction(l,u,lo.rows,weights);
        return O{{"label",r.label},{"upper_bound",normsJson(r.upperBound)},{"best_norm",r.bestNorm},{"rate_bound",r.rateBound},{"contractive_for_interval",r.contractiveForInterval},{"scope",r.scope}};
    }
    const auto spec=fieldFrom(q);const auto n=spec.states;const double par=numberOr(q,"parameter",0);
    SampledContractionInput in;in.stateDimension=n;in.lower=reals(required(q,"lower"));in.upper=reals(required(q,"upper"));in.weights=weights;
    if(in.lower.size()!=n+1||in.upper.size()!=n+1)throw std::invalid_argument("Box bounds need one entry per state plus u");
    in.samplesPerAxis=std::size_t(integerOr(q,"samples_per_axis",5,2,257));
    in.jacobian=[spec,n,par](const std::vector<double>& x){auto a=x;a.push_back(par);std::vector<double> J(n*n);spec.jacobian(a.data(),J.data());
        ComplexMatrix m(n,n);for(std::size_t k=0;k<n*n;++k)m.data[k]=J[k];return m;};
    in.field=[spec,n,par](const std::vector<double>& x){auto a=x;a.push_back(par);std::vector<double> f(n);spec.value(a.data(),f.data());return f;};
    const auto r=sampledContraction(in);
    return O{{"label",r.label},{"samples",double(r.samples)},{"worst",normsJson(r.worst)},{"best_norm",r.bestNorm},{"worst_point",numbers(r.worstPoint)},{"rate_estimate",r.rateEstimate},
        {"contractive_at_samples",r.contractiveAtSamples},{"invariant_at_samples",r.invariantAtSamples?J(*r.invariantAtSamples):J()},{"invariance_violations",double(r.invarianceViolations)},{"warnings",stringList(r.warnings)}};
}
}
std::optional<core::json::Value> evaluateResolventDataQuery(const std::string& op,const core::json::Value& q){
    if(op=="discrete_resolvent")return discreteQuery(q);
    if(op=="prony_identification")return pronyQuery(q);
    if(op=="forced_harmonic_estimate")return forcedQuery(q);
    if(op=="subharmonic_response")return subharmonicQuery(q);
    if(op=="isochron_split")return isochronQuery(q);
    if(op=="harmonic_continuation")return continuationQuery(q);
    if(op=="contraction_check")return contractionQuery(q);
    return std::nullopt;
}
}
