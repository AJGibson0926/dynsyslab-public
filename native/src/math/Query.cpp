#include <dynsyslab/math/Query.hpp>
#include <dynsyslab/math/Arithmetic.hpp>
#include <dynsyslab/math/Distributions.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/math/AnalyticFourier.hpp>
#include <dynsyslab/math/Resonance.hpp>
#include <dynsyslab/math/Topology.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
#include <dynsyslab/math/PeriodicAtlas.hpp>
#include <dynsyslab/math/PeriodicClock.hpp>
#include <dynsyslab/math/PeriodicSnapshots.hpp>
#include <cmath>
#include <stdexcept>
#include <dynsyslab/math/InvariantProjection.hpp>
#include <dynsyslab/math/DiscreteCycles.hpp>
#include <dynsyslab/math/CyclicCocycle.hpp>
#include <dynsyslab/math/OrbitQuery.hpp>
#include <dynsyslab/math/ResolventQuery.hpp>
#include <dynsyslab/math/ResolventDataQuery.hpp>
#ifdef DYNSYSLAB_MATH_CONTROL_CERTIFICATES
#include <dynsyslab/math/ControlCertificateQuery.hpp>
#endif
#include <optional>
namespace dynsyslab::math {
core::json::Value evaluateTorusMathQuery(const core::json::Value&);
core::json::Value evaluateScatteringQuery(const core::json::Value&);
std::optional<core::json::Value> evaluateSpectralDiagnosticsQuery(const std::string&,const core::json::Value&);
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
std::string str(const J& v){
    if(!std::holds_alternative<std::string>(v.storage()))throw std::invalid_argument("Expected a JSON string");
    return v.string();
}
double num(const J& v){
    if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Expected a finite JSON number");
    return v.number();
}
int integer(const J& v,int lo=0,int hi=1000000){
    const double n=num(v);if(n!=std::floor(n)||n<lo||n>hi)throw std::invalid_argument("Integer argument out of range");return int(n);
}
bool boolean(const J& v){if(!std::holds_alternative<bool>(v.storage()))throw std::invalid_argument("Expected a JSON boolean");return v.boolean();}
const A& array(const J& v){
    if(!std::holds_alternative<A>(v.storage()))throw std::invalid_argument("Expected a JSON array");
    if(v.array().size()>1000000)throw std::invalid_argument("Query array exceeds 1000000 elements");return v.array();
}
std::complex<double> complex(const J& v){
    if(std::holds_alternative<double>(v.storage()))return num(v);
    const auto& a=array(v);if(a.size()!=2)throw std::invalid_argument("Complex value needs [real,imaginary]");return {num(a[0]),num(a[1])};
}
J encode(std::complex<double> z){return A{z.real(),z.imag()};}
template<class T>A numbers(const std::vector<T>& values){A out;for(auto x:values)out.emplace_back(double(x));return out;}
A strings(const std::vector<std::string>& values){A out;for(const auto& x:values)out.emplace_back(x);return out;}
std::vector<std::complex<double>> complexes(const J& value){std::vector<std::complex<double>> out;for(const auto& x:array(value))out.push_back(complex(x));return out;}
std::vector<double> reals(const J& value){std::vector<double> out;for(const auto& x:array(value))out.push_back(num(x));return out;}
double numberOr(const J& q,const char* key,double fallback){return q.at(key).isNull()?fallback:num(q.at(key));}
unsigned orderOr(const J& q,unsigned fallback){return q.at("order").isNull()?fallback:unsigned(integer(q.at("order"),0,64));}
J optional(const std::optional<double>& v){return v?J(*v):J();}
J resonance(const HarmonicResonance& r){return O{{"harmonic_a",double(r.harmonicA)},{"harmonic_b",double(r.harmonicB)},
    {"detuning",r.detuning},{"phase_drift",r.relativePhaseDrift},{"ratio_error",r.ratioError},{"beat_time",optional(r.beatTime)}};}
FourierSeries series(const J& q){
    FourierSeries s{num(q.at("angular_frequency")),{}};
    for(const auto& row:array(q.at("coefficients"))){const auto& a=array(row);
        if(a.size()!=3)throw std::invalid_argument("Fourier coefficient row needs [harmonic,real,imaginary]");
        const int n=integer(a[0],-1000000,1000000);
        if(!s.coefficients.emplace(n,std::complex<double>(num(a[1]),num(a[2]))).second)throw std::invalid_argument("Duplicate harmonic");
    }
    return s;
}
J seriesJson(const FourierSeries& s){A rows;for(const auto& [n,c]:s.coefficients)rows.emplace_back(A{double(n),c.real(),c.imag()});return O{{"angular_frequency",s.angularFrequency},{"coefficients",rows}};}
PointDistribution distribution(const J& q){return {num(q.at("location")),complexes(q.at("coefficients"))};}
J distributionJson(const PointDistribution& d){A c;for(auto v:d.coefficients)c.push_back(encode(v));return O{{"location",d.location},{"coefficients",c}};}
J dispatch(const J& q,const std::string& op){
    if(auto orbit=evaluateOrbitQuery(op,q))return *orbit;
    if(auto diagnostics=evaluateSpectralDiagnosticsQuery(op,q))return *diagnostics;
    if(auto resolvent=evaluateResolventQuery(op,q))return *resolvent;
    if(auto resolventData=evaluateResolventDataQuery(op,q))return *resolventData;
#ifdef DYNSYSLAB_MATH_CONTROL_CERTIFICATES
    if(auto certificate=evaluateControlCertificateQuery(op,q))return *certificate;
#endif
    if(op=="phase_cover_resonance") {
        const auto r=phaseCoverResonance(integer(q.at("prime"),2),str(q.at("base_to_drive_ratio")));
        return O{{"base_to_drive_ratio",r.baseToDriveRatio},{"character_ratio",r.characterRatio},
            {"drive_harmonic",r.driveHarmonic},{"minimum_depth",double(r.minimumDepth)},
            {"branch_classes",r.branchClasses},{"fiber_blind",r.minimumDepth==0},
            {"scope","exact rational frequency relation; integer multiples of this generator; no physical torus-destruction claim"}};
    }
    if(op=="cycle_cocycle"||op=="cycle_partition_geometry")return evaluateCyclicCocycleQuery(q);
    if(op=="sampled_cycles"||op=="cycle_spectrum"||op=="cycle_projection"||op=="cycle_alignment"||op=="sampled_hankel_pairs")return evaluateDiscreteCycleQuery(q);
    if(op=="periodic_delay_readout"){
        std::vector<FourierSeries> channels;for(const auto&s:array(q.at("channels")))channels.push_back(series(s));
        if(!q.at("readouts").isNull()){
            A outputs;const auto offsets=reals(q.at("time_offsets"));
            for(const auto&r:array(q.at("readouts"))){
                const auto out=periodicDelayReadout(channels,offsets,complexes(r.at("coefficients")));
                O item{{"mean_square",meanSquare(out)}};
                if(!r.at("reference_channel").isNull()){
                    const auto c=integer(r.at("reference_channel"),0,int(channels.size())-1);
                    const RationalCircleKernel k(num(q.at("bandwidth")),integer(q.at("order"),1,8));
                    const auto ref=spectralPacket(channels[c],num(q.at("sample_interval")),num(q.at("angle")),k);
                    auto difference=out;for(const auto&[n,a]:ref.coefficients)difference.coefficients[n]-=a;
                    item["reference_mean_square"]=meanSquare(ref);item["error_mean_square"]=meanSquare(difference);
                }
                if(q.at("summary_only").isNull()||!boolean(q.at("summary_only")))item["series"]=seriesJson(out);
                outputs.emplace_back(item);
            }
            return O{{"readouts",outputs}};
        }
        const auto out=periodicDelayReadout(channels,reals(q.at("time_offsets")),complexes(q.at("coefficients")));
        return O{{"series",seriesJson(out)},{"mean_square",meanSquare(out)}};
    }
    if(op=="periodic_clock"){
        const auto w=series(q.at("clock_rate"));
        if(str(q.at("action"))=="coefficient_expansion"){
            const auto r=periodicClockCoefficientExpansion(series(q.at("series")),w,
                integer(q.at("harmonic"),-1000000),unsigned(integer(q.at("order"),0,128)));
            return O{{"coefficient",encode(r.coefficient)},{"exact_arithmetic_remainder_bound",r.exactArithmeticRemainderBound},
                {"bound_scope","Taylor remainder for supplied finite polynomials; excludes roundoff and omitted input terms; caller must supply a positive clock"}};
        }
        if(str(q.at("action"))!="transform")throw std::invalid_argument("Unknown periodic clock action");
        std::vector<FourierSeries> channels;for(const auto& g:array(q.at("channels")))channels.push_back(series(g));
        const auto r=reparameterizePeriodicFourier(channels,w,integer(q.at("samples"),3),integer(q.at("maximum_harmonic")));
        A transformed;for(const auto& g:r.channels)transformed.push_back(seriesJson(g));
        return O{{"channels",transformed},{"phase_correction",seriesJson(r.phaseCorrection)},
            {"physical_period",r.physicalPeriod},{"mean_clock_rate",r.meanClockRate},
            {"minimum_sampled_clock_rate",r.minimumSampledClockRate},{"clock_polynomial_lower_bound",r.clockPolynomialLowerBound},
            {"certified",false},{"scope","Auxiliary-phase trapezoidal coefficient quadrature with clock Jacobian; refine sample count; input/output tail errors excluded"}};
    }
    if(op=="analytic_fourier"){
        const auto kind=str(q.at("kind"));const int n=integer(q.at("maximum_harmonic"));
        const double w=numberOr(q,"angular_frequency",1);
        if(kind=="kepler_cos_eccentric"||kind=="kepler_sin_eccentric")
            return O{{"series",seriesJson(keplerAnomalyFourier(kind=="kepler_sin_eccentric",
                num(q.at("parameter")),n,w))},{"uniform_tail_bound",J{}},
                {"parameter_period",2*std::acos(-1.)},{"bound_scope","No tail bound supplied; parameter is eccentricity; Fourier angle is mean anomaly"}};
        const auto result=kind=="log_cosine"?logarithmicCosineFourier(num(q.at("parameter")),n,w):jacobiFourier(kind,num(q.at("parameter")),n,w);
        return O{{"series",seriesJson(result.series)},{"uniform_tail_bound",result.uniformTailBound},
            {"parameter_period",result.parameterPeriod},{"bound_scope","Infinite-series truncation in exact arithmetic; not a floating-point enclosure"}};
    }
    if(op=="invariant_readout_floor") {
        const auto r=zeroMeanReadoutFloor(complexes(q.at("target_component_means")),
            reals(q.at("component_weights")),q.at("allow_constant").isNull()?true:boolean(q.at("allow_constant")));
        return O{{"mean_squared_floor",r.meanSquaredFloor},{"best_constant",encode(r.bestConstant)},
            {"total_weight",r.totalWeight},{"certified",false},
            {"scope","Assumes exact zero conditional means of all source observables; numerical supplied-measure estimate; linear time-translate readouts only"}};
    }
    if(op=="inverse_scattering")return evaluateScatteringQuery(q);
    if(op=="torus_clock"||op=="phase_cover_forcing"||op=="torus_fourier"||op=="small_divisors"||op=="normal_form"||op=="melnikov"||op=="symplectic_defect"||op=="hamiltonian_first_order"||op=="hamiltonian_twist")return evaluateTorusMathQuery(q);
    if(op=="periodic_atlas"){
        std::vector<PeriodicAtlasKnot> knots;
        for(const auto& k:array(q.at("knots"))){PeriodicAtlasKnot knot;knot.parameter=num(k.at("parameter"));
            for(const auto& s:array(k.at("channels")))knot.channels.push_back(series(s));knots.push_back(std::move(knot));}
        const PeriodicAtlas atlas(std::move(knots));const double p=num(q.at("parameter"));
        const auto interpolated=atlas.interpolate(p);A channels,values;
        for(const auto& s:interpolated)channels.push_back(seriesJson(s));
        for(auto c:atlas.predict(p,num(q.at("phase")),numberOr(q,"elapsed",0)))values.push_back(encode(c));
        return O{{"channels",channels},{"values",values}};
    }
    if(op=="spectral_filter"){
        const RationalCircleKernel kernel(num(q.at("bandwidth")),q.at("order").isNull()?1:integer(q.at("order"),1,8));
        const auto action=str(q.at("action"));
        if(action=="kernel"){
            A values,lags;for(double x:reals(q.at("angles")))values.emplace_back(kernel(x));
            if(!q.at("lags").isNull())for(const auto& lag:array(q.at("lags")))lags.emplace_back(kernel.timeMultiplier(integer(lag,-1000000)));
            O out{{"values",values},{"time_multipliers",lags}};
            if(!q.at("tail_after_lag").isNull())out["symmetric_time_tail_bound"]=kernel.symmetricTimeTailBound(integer(q.at("tail_after_lag")));
            return out;
        }
        if(action=="packet")return seriesJson(spectralPacket(series(q.at("series")),num(q.at("sample_interval")),num(q.at("angle")),kernel));
        if(action=="density"){
            std::vector<WeightedOrbitSpectrum> orbits;
            for(const auto& o:array(q.at("orbits")))orbits.push_back({series(o.at("series")),num(o.at("weight"))});
            const auto s=smoothOrbitSpectra(orbits,num(q.at("sample_interval")),reals(q.at("angles")),kernel);
            return O{{"density",numbers(s.density)},{"total_mass",s.totalMass},{"zero_harmonic_mass",s.zeroHarmonicMass}};
        }
        throw std::invalid_argument("Unknown spectral filter action");
    }
    if(op=="frequency_band_prediction"){
        const auto p=predictFrequencyBands(series(q.at("series")),num(q.at("initial_phase")),reals(q.at("elapsed_times")),num(q.at("band_width")));
        A values;for(auto v:p.values)values.push_back(encode(v));
        return O{{"values",values},{"known_coefficient_error_bound",numbers(p.knownCoefficientErrorBound)}};
    }
    if(op=="rational"){
        const auto a=str(q.at("a")),b=str(q.at("b"));
        return O{{"a",canonicalRational(a)},{"b",canonicalRational(b)},{"difference",rationalDifference(a,b)},{"ratio",rationalRatio(a,b)}};
    }
    if(op=="padic"){
        const auto a=str(q.at("a")),b=str(q.at("b"));const auto p=integer(q.at("prime"),2);
        auto e=padicExpansion(a,p,std::size_t(integer(q.at("precision"),1,4096)));
        return O{{"distance",padicDistance(a,b,p)},{"valuation",e.valuation?J(double(*e.valuation)):J()},
            {"valuation_is_infinite",!e.valuation.has_value()},{"digits",numbers(e.digits)}};
    }
    if(op=="phase_cover_character"){
        std::vector<std::uint32_t> digits;for(const auto&d:array(q.at("digits")))digits.push_back(integer(d,0,999999));
        const auto c=phaseCoverCharacter(integer(q.at("prime"),2,1000000),digits,
            str(q.at("base_phase_turns")),str(q.at("frequency_ratio")));
        return O{{"frequency_ratio",c.frequencyRatio},{"minimum_depth",double(c.minimumDepth)},
            {"resolved",c.phaseTurns.has_value()},{"phase_turns",c.phaseTurns?J(*c.phaseTurns):J()},
            {"required_additional_digits",double(c.minimumDepth>digits.size()?c.minimumDepth-digits.size():0)},
            {"scope","Character value is exp(2*pi*i*phase_turns). Unspecified higher digits are not inferred. Minimum depth identifies the factor circle on which the character is single-valued."}};
    }
    if(op=="phase_cover_flow"){
        std::vector<std::uint32_t> digits;for(const auto&d:array(q.at("digits")))digits.push_back(integer(d,0,999999));
        const auto result=advancePhaseCover(integer(q.at("prime"),2,1000000),digits,
            str(q.at("base_phase_turns")),str(q.at("elapsed_base_periods")));
        A levels;for(std::size_t j=0;j<result.levels.size();++j){const auto&l=result.levels[j];
            levels.emplace_back(O{{"modulus",l.modulus},{"residue",l.residue},{"frequency_ratio",l.frequencyRatio},{"phase_turns",result.phaseTurns[j]}});}
        return O{{"base_phase_turns",result.basePhaseTurns},{"return_count",result.returnCount},
            {"prefix_carry",result.prefixCarry},{"digits",numbers(result.digits)},{"levels",levels},
            {"finite_prefix",true},{"scope","Exact rational translation on a finite circle-cover prefix. prefix_carry acts on unspecified higher digits; a repeated prefix is not a return of the infinite state. No physical orbit or cascade is inferred."}};
    }
    if(op=="phase_tower"){
        std::vector<std::uint32_t> digits;for(const auto&d:array(q.at("digits")))digits.push_back(integer(d,0,999999));
        A levels;for(const auto&l:phaseCoverTower(integer(q.at("prime"),2,1000000),digits))
            levels.emplace_back(O{{"modulus",l.modulus},{"residue",l.residue},{"frequency_ratio",l.frequencyRatio},{"phase_offset_turns",l.offsetTurns}});
        return O{{"levels",levels},{"finite_prefix",true},{"phase_convention","base_phase / modulus + 2*pi*phase_offset_turns"}};
    }
    if(op=="phase_compatibility")return O{{"compatible",compatiblePhaseResidues(str(q.at("lower")),str(q.at("upper")),
        integer(q.at("prime"),2),integer(q.at("lower_depth"),0,4096))}};
    if(op=="continued_fraction"){
        const auto value=str(q.at("value"));auto cf=continuedFraction(value);A path;
        for(const auto& run:sternBrocotPath(value))path.emplace_back(O{{"direction",std::string(1,run.direction)},{"repetitions",run.repetitions}});
        return O{{"terms",strings(cf.terms)},{"convergents",strings(cf.convergents)},{"stern_brocot_runs",path}};
    }
    if(op=="ford"){
        auto a=str(q.at("a")),b=str(q.at("b"));auto ca=fordCircle(a),cb=fordCircle(b);
        return O{{"a",O{{"abscissa",ca.abscissa},{"radius",ca.radius}}},{"b",O{{"abscissa",cb.abscissa},{"radius",cb.radius}}},
            {"tangent",fordCirclesTangent(a,b)},{"mediant",rationalMediant(a,b)}};
    }
    if(op=="multiplicative_distance")return O{{"exponential_distance",multiplicativeDistanceProduct(str(q.at("a")),str(q.at("b")))}};
    if(op=="resonance"){
        double a=num(q.at("frequency_a")),b=num(q.at("frequency_b")),t=numberOr(q,"horizon",0);A candidates;
        for(auto r:resonanceCandidates(a,b,integer(q.at("maximum_order"),1),t))candidates.push_back(resonance(r));
        return O{{"candidates",candidates},{"exact_rationality_certified",false}};
    }
    if(op=="graph"){
        FiniteGraph g;g.vertices=integer(q.at("vertices"));for(const auto& row:array(q.at("edges"))){
            const auto& e=array(row);if(e.size()<2||e.size()>3)throw std::invalid_argument("Edge needs [from,to,length?]");
            g.edges.push_back({std::size_t(integer(e[0])),std::size_t(integer(e[1])),e.size()==3?num(e[2]):1.});}
        auto t=graphTopology(g);O out{{"components",double(t.connectedComponents)},{"first_betti",double(t.firstBettiNumber)},
            {"component_labels",numbers(t.component)},{"degrees",numbers(t.degree)},{"forest_edges",numbers(t.spanningForestEdges)},{"cycle_edges",numbers(t.cycleClosingEdges)}};
        if(!q.at("source").isNull()){A distances;for(auto d:graphDistances(g,integer(q.at("source"))))distances.push_back(optional(d));out["distances"]=distances;}
        if(!q.at("vertex_values").isNull()){A bars;for(auto p:lowerStarPersistence(g,reals(q.at("vertex_values"))))
            bars.emplace_back(O{{"birth_vertex",double(p.birthVertex)},{"birth",p.birth},{"death",optional(p.death)}});out["h0_intervals"]=bars;}
        return out;
    }
    if(op=="fourier"){
        const auto action=str(q.at("action"));
        if(action=="from_samples")return seriesJson(periodicFourierSeries(complexes(q.at("samples")),
            num(q.at("angular_frequency")),integer(q.at("maximum_harmonic"))));
        const auto s=series(q.at("series"));
        if(action=="sample_values"){A out;for(auto c:periodicGridValues(s,integer(q.at("samples"),3),numberOr(q,"phase_origin",0)))out.push_back(encode(c));return O{{"samples",out}};}
        if(action=="to_samples"){A out;for(auto c:uniformFourierSamples(s,integer(q.at("samples"),3)))out.push_back(encode(c));return O{{"samples",out}};}
        if(action=="evaluate")return O{{"value",encode(evaluate(s,num(q.at("phase"))))}};
        if(action=="evolve"){
            A values;for(auto c:evolveFourierSeries(s,reals(q.at("time")),numberOr(q,"initial_phase",0),numberOr(q,"initial_time",0)))values.push_back(encode(c));
            return O{{"values",values},{"scope","Free evolution of supplied Fourier characters at the declared frequency; no fitting, forcing, or omitted-tail estimate"}};
        }
        if(action=="evaluate_many"){
            const auto phases=reals(q.at("phases"));
            if(!s.coefficients.empty()&&phases.size()>100000000/s.coefficients.size())
                throw std::invalid_argument("Fourier evaluation work exceeds limit");
            A values;for(double phase:phases)values.push_back(encode(evaluate(s,phase)));
            return O{{"values",values},{"scope","Pointwise evaluation of the supplied finite series; omitted harmonics excluded"}};
        }
        if(action=="differentiate")return seriesJson(timeDerivative(s,orderOr(q,1)));
        if(action=="primitive")return seriesJson(periodicPrimitive(s));
        if(action=="multiply")return seriesJson(multiply(s,series(q.at("other"))));
        if(action=="cyclic_sector")return seriesJson(projectCyclicSector(s,integer(q.at("degree"),1),integer(q.at("residue"),-1000000)));
        if(action=="norms")return O{{"mean_square",meanSquare(s)},{"omitted_known_terms_bound",omittedUniformBound(s,integer(q.at("maximum_harmonic")))}};
        throw std::invalid_argument("Unknown Fourier action");
    }
    if(op=="distribution"){
        auto d=distribution(q.at("distribution"));const auto action=str(q.at("action"));
        if(action=="pair")return O{{"value",encode(pairTestJet(d,complexes(q.at("test_derivatives"))))}};
        if(action=="differentiate")return distributionJson(distributionDerivative(d,orderOr(q,1)));
        if(action=="multiply_smooth")return distributionJson(multiplySmooth(d,complexes(q.at("smooth_derivatives"))));
        throw std::invalid_argument("Unknown distribution action");
    }
    if(op=="orbit_pairing")return O{{"value",encode(periodicOrbitPairing(complexes(q.at("test_samples")),
        integer(q.at("harmonic"),-1000000),num(q.at("period")),q.at("normalized").isNull()?true:boolean(q.at("normalized"))))}};
    if(op=="coarea")return O{{"value",encode(levelSetDeltaPairing(complexes(q.at("test_values")),reals(q.at("surface_weights")),reals(q.at("gradient_norms"))))}};
    if(op=="spectral_density"){
        std::vector<SpectralDensityRoot> roots;for(const auto& r:array(q.at("roots")))roots.push_back({
            num(r.at("source_density")),num(r.at("frequency_derivative")),integer(r.at("harmonic"),-1000000),complex(r.at("coefficient"))});
        return O{{"density",spectralPushforwardDensity(roots)}};
    }
    if(op=="gaussian_delta")return O{{"value",gaussianDeltaDerivative(num(q.at("x")),num(q.at("epsilon")),orderOr(q,0))}};
    throw std::invalid_argument("Unknown mathematics operation: "+op);
}
}
core::json::Value evaluateQuery(const core::json::Value& q){
    if(str(q.at("schema"))!="dynsyslab.math-query/v1")throw std::invalid_argument("Expected dynsyslab.math-query/v1");
    auto op=str(q.at("operation"));return O{{"schema","dynsyslab.mathematics/v1"},{"operation",op},{"result",dispatch(q,op)}};
}
}
