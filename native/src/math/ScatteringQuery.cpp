#include <dynsyslab/math/InverseScattering.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Toda.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::math {
namespace {
using J=core::json::Value;using A=J::Array;using O=J::Object;
double num(const J&v){if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite number required");return v.number();}
const A&array(const J&v){if(!std::holds_alternative<A>(v.storage())||v.array().size()>1000000)throw std::invalid_argument("Bounded array required");return v.array();}
std::vector<double> reals(const J&v){std::vector<double>r;for(const auto&x:array(v))r.push_back(num(x));return r;}
template<class T>A values(const T&v){A r;for(auto x:v)r.emplace_back(double(x));return r;}
J complex(std::complex<double>z){return A{z.real(),z.imag()};}
}
core::json::Value evaluateScatteringQuery(const J&q){const auto action=q.at("action").string();
    if(action=="toda_lax") {
        if(!q.at("periodic").isNull()&&!std::holds_alternative<bool>(q.at("periodic").storage()))
            throw std::invalid_argument("periodic must be a boolean");
        const bool periodic=q.at("periodic").isNull()?true:q.at("periodic").boolean();
        const auto r=core::todaLaxPair(reals(q.at("positions")),reals(q.at("momenta")),
            q.at("alpha").isNull()?1:num(q.at("alpha")),periodic);
        return O{{"sites",double(r.sites)},{"matrix",values(r.matrix)},{"generator",values(r.generator)},
            {"a",values(r.a)},{"b",values(r.b)},{"periodic",periodic},
            {"scope","Ldot=[P,L]; periodic L is the +1 Bloch matrix; Lax eigenvalues are not Koopman eigenvalues"}};
    }
    if(action=="right_scattering"){
        auto r=rightScattering(reals(q.at("potential")),num(q.at("left")),num(q.at("cell_width")),num(q.at("wavenumber")));
        return O{{"reflection",complex(r.reflection)},{"transmission",complex(r.transmission)},{"flux_defect",r.fluxDefect},{"potential_model","piecewise constant cells, zero exterior"}};
    }
    if(action=="hill_monodromy"){
        auto r=hillMonodromy(reals(q.at("potential")),num(q.at("cell_width")),num(q.at("energy")));
        return O{{"matrix",values(r.matrix)},{"discriminant",r.discriminant},{"determinant_defect",r.determinantDefect},{"in_band_numerically",std::abs(r.discriminant)<=2},{"potential_model","periodic extension of piecewise constant cells"}};
    }
    if(action=="dirichlet_negative_spectrum")return O{{"eigenvalues",values(dirichletNegativeSpectrum(reals(q.at("potential")),num(q.at("length"))))},{"scope","finite-box centered-difference Dirichlet operator; not a whole-line certification"}};
    if(action=="evolve_reflection"){
        const auto&v=array(q.at("reflection"));if(v.size()!=2)throw std::invalid_argument("Reflection needs [real,imaginary]");
        return O{{"reflection",complex(evolveKdVReflection({num(v[0]),num(v[1])},num(q.at("wavenumber")),num(q.at("elapsed"))))}};
    }
    ReflectionlessData d{reals(q.at("kappa")),reals(q.at("log_norming_squared"))};
    if(!q.at("elapsed").isNull())d=evolveKdV(d,num(q.at("elapsed")));
    if(action=="evolve_norming")return O{{"kappa",values(d.kappa)},{"log_norming_squared",values(d.logNormingSquared)}};
    if(action=="reflectionless_potential"){A v;for(double x:reals(q.at("positions")))v.emplace_back(reflectionlessPotential(d,x));
        return O{{"values",v},{"equation","V_t - 6 V V_x + V_xxx = 0"},{"log_norming_squared",values(d.logNormingSquared)}};
    }
    throw std::invalid_argument("Unknown inverse-scattering operation");
}
}
