#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/Perturbation.hpp>
#include <dynsyslab/math/HamiltonianPerturbation.hpp>
#include <dynsyslab/math/PhaseCoverFourier.hpp>
#include <dynsyslab/math/TorusClock.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
int integer(const J&v,int lo=0,int hi=1000000){double n=num(v);if(n<lo||n>hi||std::floor(n)!=n)throw std::invalid_argument("Integer out of range");return int(n);}
const A&array(const J&v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>1000000)throw std::invalid_argument("Bounded array required");return v.array();}
std::complex<double> complex(const J&v){if(std::holds_alternative<double>(v.storage()))return num(v);const auto&a=array(v);if(a.size()!=2)throw std::invalid_argument("Complex number needs [real,imag]");return {num(a[0]),num(a[1])};}
std::vector<double> reals(const J&v){std::vector<double>out;for(const auto&a:array(v))out.push_back(num(a));return out;}
template<class T>A values(const std::vector<T>&v){A out;for(auto a:v)out.emplace_back(double(a));return out;}
J encode(std::complex<double>v){return A{v.real(),v.imag()};}
TorusFourierSeries series(const J&q){TorusFourierSeries s{reals(q.at("angular_frequencies")),{}};
    for(const auto&row:array(q.at("coefficients"))){FourierIndex k;for(const auto&n:array(row.at("index")))k.push_back(integer(n,-1000000));
        if(!s.coefficients.emplace(k,complex(row.at("coefficient"))).second)throw std::invalid_argument("Duplicate Fourier index");}return s;}
J encode(const TorusFourierSeries&s){A c;for(const auto&[k,v]:s.coefficients)c.emplace_back(O{{"index",values(k)},{"coefficient",encode(v)}});
    return O{{"angular_frequencies",values(s.angularFrequencies)},{"coefficients",c}};}
A monomials(const std::vector<NormalFormMonomial>&v){A a;for(const auto&m:v)a.emplace_back(O{{"component",double(m.component)},{"powers",values(m.powers)},{"coefficient",encode(m.coefficient)}});return a;}
}
core::json::Value evaluateTorusMathQuery(const core::json::Value&q){
    const auto op=q.at("operation").string();
    if(op=="torus_clock") {
        const auto clock=series(q.at("clock_rate"));
        const double divisor=q.at("divisor_tolerance").isNull()?1e-12:num(q.at("divisor_tolerance"));
        TorusClockTransform transform;A rows;
        if(q.at("time").isNull())transform=analyzeTorusClock(clock,divisor);
        else {
            std::vector<FourierIndex> characters;
            if(!q.at("characters").isNull())for(const auto& row:array(q.at("characters"))) {
                FourierIndex k;for(const auto& n:array(row))k.push_back(integer(n,-1000000));characters.push_back(k);
            }
            const auto result=evolveTorusClock(clock,reals(q.at("initial_phases")),reals(q.at("time")),characters,
                divisor,q.at("relative_time_tolerance").isNull()?1e-12:num(q.at("relative_time_tolerance")));
            transform=result.transform;
            for(const auto& r:result.samples) {
                A bare,corrected;for(auto z:r.bareCharacters)bare.push_back(encode(z));
                for(auto z:r.correctedCharacters)corrected.push_back(encode(z));
                rows.emplace_back(O{{"time",r.time},{"base_elapsed",r.baseElapsed},{"time_residual",r.timeResidual},
                    {"clock_rate",r.clockRate},{"phases",values(r.phases)},{"straightened_phases",values(r.straightenedPhases)},
                    {"bare_characters",bare},{"corrected_characters",corrected}});
            }
        }
        return O{{"mean_clock_rate",transform.meanRate},{"clock_lower_bound",transform.lowerBound},
            {"clock_upper_bound",transform.upperBound},{"primitive_l1_bound",transform.primitiveL1Bound},
            {"primitive_mean_square",meanSquare(transform.primitive)},
            {"primitive",encode(transform.primitive)},{"straightened_frequencies",values(transform.straightenedFrequencies)},
            {"samples",rows},{"interval_certified",false},{"infinite_limit_certified",false},
            {"scope","positive supplied finite clock polynomial; exact-arithmetic l1 bound; floating inversion; omitted tails excluded"}};
    }
    if(op=="phase_cover_forcing") {
        // Reuse the standard sparse Fourier input parser and cohomological solver.
        auto input=q;input.object()["angular_frequencies"]=A{1.};
        const auto coefficients=series(input).coefficients;
        const auto r=phaseCoverForcing(integer(q.at("prime"),2),integer(q.at("depth"),0,32),
            num(q.at("base_angular_frequency")),reals(q.at("forcing_frequencies")),
            coefficients,num(q.at("divisor_tolerance")));
        A modes;
        for(const auto& m:r.modes)modes.emplace_back(O{{"index",values(m.index)},
            {"frequency_ratio",m.frequencyRatio},{"minimum_depth",double(m.minimumDepth)},
            {"detuning",m.detuning},{"forcing_magnitude",m.forcingMagnitude},
            {"correction_magnitude",m.retained?J{}:J(m.correctionMagnitude)},
            {"retained",m.retained}});
        return O{{"cover_degree",r.coverDegree},{"forcing",encode(r.forcing)},
            {"solution",encode(r.response.solution)},{"retained",encode(r.response.retained)},
            {"modes",modes},{"maximum_inverse_divisor",r.response.maximumInverseDivisor},
            {"retained_mean_square",meanSquare(r.response.retained)},
            {"finite_truncation",true},{"kam_certified",false},
            {"scope","supplied finite character coefficients; floating detunings; no tail or nonlinear remainder bound"}};
    }
    if(op=="hamiltonian_twist") {
        const auto r=hamiltonianTwist(reals(q.at("angular_frequencies")),reals(q.at("hessian")),
            q.at("relative_tolerance").isNull()?1e-10:num(q.at("relative_tolerance")));
        return O{{"kolmogorov_singular_values",values(r.kolmogorovSingularValues)},
            {"isoenergetic_singular_values",values(r.isoenergeticSingularValues)},
            {"kolmogorov_relative_minimum",r.kolmogorovRelativeMinimum},
            {"isoenergetic_relative_minimum",r.isoenergeticRelativeMinimum},
            {"kolmogorov_numerically_nondegenerate",r.kolmogorovNumericallyNondegenerate},
            {"isoenergetic_numerically_nondegenerate",r.isoenergeticNumericallyNondegenerate},
            {"kam_certified",false}};
    }
    if(op=="hamiltonian_first_order") {
        std::vector<TorusFourierSeries> gradient;
        for(const auto& g:array(q.at("action_gradient")))gradient.push_back(series(g));
        const auto r=firstOrderHamiltonian(series(q.at("perturbation")),gradient,
            reals(q.at("hessian")),num(q.at("divisor_tolerance")));
        A actions,angles,retainedGradient;
        for(const auto& s:r.actionCorrection)actions.push_back(encode(s));
        for(const auto& s:r.angleCorrection)angles.push_back(encode(s));
        for(const auto& s:r.retainedActionGradient)retainedGradient.push_back(encode(s));
        return O{{"generator",encode(r.generator)},{"retained",encode(r.retained)},
            {"action_correction",actions},{"angle_correction",angles},
            {"retained_action_gradient",retainedGradient},{"mean_frequency_shift",values(r.meanFrequencyShift)},
            {"maximum_inverse_divisor",r.maximumInverseDivisor},{"kam_certified",false},
            {"scope","first-order coefficients; old=identity+epsilon*X_chi; no remainder bound"}};
    }
    if(op=="torus_fourier"){
        const auto action=q.at("action").string();
        if(action=="from_samples"){
            std::vector<std::complex<double>>samples;for(const auto&v:array(q.at("samples")))samples.push_back(complex(v));
            std::vector<std::size_t>shape;for(const auto&v:array(q.at("shape")))shape.push_back(integer(v,3));
            std::vector<int>maximum;for(const auto&v:array(q.at("maximum_harmonics")))maximum.push_back(integer(v));
            return encode(torusFourierSeries(samples,shape,reals(q.at("angular_frequencies")),maximum));
        }
        const auto s=series(q.at("series"));
        if(action=="evaluate")return O{{"value",encode(evaluate(s,reals(q.at("phases"))))}};
        if(action=="time_derivative")return encode(torusTimeDerivative(s));
        if(action=="partial_derivative")return encode(torusPartialDerivative(s,integer(q.at("axis"))));
        if(action=="mean_square")return O{{"mean_square",meanSquare(s)}};
        if(action=="cohomological"){
            const auto r=solveCohomological(s,num(q.at("divisor_tolerance")));
            return O{{"solution",encode(r.solution)},{"retained",encode(r.retained)},{"maximum_inverse_divisor",r.maximumInverseDivisor},{"finite_truncation",true}};
        }throw std::invalid_argument("Unknown torus Fourier action");
    }
    if(op=="small_divisors"){
        auto r=smallDivisors(reals(q.at("angular_frequencies")),integer(q.at("maximum_order"),1,1000),num(q.at("exponent")),q.at("retain").isNull()?20:integer(q.at("retain"),0,100000));A a;
        for(const auto&v:r)a.emplace_back(O{{"index",values(v.index)},{"divisor",v.divisor},{"scaled_divisor",v.scaledDivisor}});
        return O{{"records",a},{"diophantine_certified",false},{"scope","finite integer lattice only"}};
    }
    if(op=="normal_form"){
        std::vector<std::complex<double>>lambda;for(const auto&v:array(q.at("eigenvalues")))lambda.push_back(complex(v));
        std::vector<NormalFormMonomial>terms;for(const auto&v:array(q.at("terms"))){NormalFormMonomial m;
            m.component=integer(v.at("component"));for(const auto&p:array(v.at("powers")))m.powers.push_back(integer(p,0,64));m.coefficient=complex(v.at("coefficient"));terms.push_back(m);}
        auto r=homogeneousNormalForm(lambda,terms,num(q.at("divisor_tolerance")));
        return O{{"change",monomials(r.change)},{"retained",monomials(r.retained)},{"degree",double(r.degree)},{"minimum_removed_divisor",r.minimumRemovedDivisor},{"convergence_certified",false}};
    }
    if(op=="melnikov")return O{{"sampled_integral",sampledMelnikov(reals(q.at("time")),reals(q.at("drift")),reals(q.at("perturbation")))},{"infinite_tail_bounded",false}};
    if(op=="symplectic_defect")return O{{"frobenius_defect",conformalSymplecticDefect(reals(q.at("jacobian")),reals(q.at("form")),integer(q.at("dimension"),1,128),q.at("damping").isNull()?0:num(q.at("damping")))}};
    throw std::invalid_argument("Unknown torus/perturbation query");
}
}
