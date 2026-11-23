#pragma once
// Private JSON parsing and encoding helpers for ResolventQuery.cpp.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/math/KoopmanFrequencyResponse.hpp>
#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <dynsyslab/math/LaplaceTransform.hpp>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
namespace dynsyslab::math::resolvent_query {
using J=core::json::Value;using A=J::Array;using O=J::Object;using C=std::complex<double>;
inline constexpr std::size_t maximumElements=1000000;
inline double num(const J& v){
    if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Expected a finite JSON number");
    return v.number();
}
inline long long integer(const J& v,long long lo,long long hi){
    const double n=num(v);
    if(n!=std::floor(n)||n<double(lo)||n>double(hi))throw std::invalid_argument("Integer argument out of range");
    return (long long)n;
}
inline const A& array(const J& v){
    if(!std::holds_alternative<A>(v.storage()))throw std::invalid_argument("Expected a JSON array");
    if(v.array().size()>maximumElements)throw std::invalid_argument("Query array exceeds 1000000 elements");
    return v.array();
}
inline const J& object(const J& v){
    if(!std::holds_alternative<O>(v.storage()))throw std::invalid_argument("Expected a JSON object");
    return v;
}
inline const J& required(const J& q,const char* key){
    const auto& v=q.at(key);if(v.isNull())throw std::invalid_argument(std::string("Missing required field: ")+key);
    return v;
}
inline bool has(const J& q,const char* key){return !q.at(key).isNull();}
inline std::string text(const J& v){
    if(!std::holds_alternative<std::string>(v.storage()))throw std::invalid_argument("Expected a JSON string");
    return v.string();
}
inline bool flag(const J& q,const char* key){
    const auto& v=q.at(key);if(v.isNull())return false;
    if(!std::holds_alternative<bool>(v.storage()))throw std::invalid_argument("Expected a JSON boolean");
    return v.boolean();
}
inline double numberOr(const J& q,const char* key,double fallback){return has(q,key)?num(q.at(key)):fallback;}
inline long long integerOr(const J& q,const char* key,long long fallback,long long lo,long long hi){
    return has(q,key)?integer(q.at(key),lo,hi):fallback;
}
inline C complexValue(const J& v){
    if(std::holds_alternative<double>(v.storage()))return num(v);
    const auto& a=array(v);if(a.size()!=2)throw std::invalid_argument("Complex value needs [real,imaginary]");
    return {num(a[0]),num(a[1])};
}
inline ComplexVector complexes(const J& v){ComplexVector out;out.reserve(array(v).size());for(const auto& x:array(v))out.push_back(complexValue(x));return out;}
inline std::vector<double> reals(const J& v){std::vector<double> out;out.reserve(array(v).size());for(const auto& x:array(v))out.push_back(num(x));return out;}
inline J encode(C z){return A{z.real(),z.imag()};}
inline J encodeAll(const ComplexVector& v){A out;out.reserve(v.size());for(const auto& z:v)out.push_back(encode(z));return out;}
template<class T>J numbers(const std::vector<T>& v){A out;out.reserve(v.size());for(auto x:v)out.emplace_back(double(x));return out;}
inline J optionalNumber(const std::optional<double>& v){return v?J(*v):J();}
inline ComplexMatrix matrix(const J& v){
    const auto& rows=array(v);if(rows.empty())throw std::invalid_argument("Matrix needs rows");
    const auto cols=array(rows[0]).size();if(cols==0)throw std::invalid_argument("Matrix needs columns");
    if(cols>4096||rows.size()>4096)throw std::invalid_argument("Matrix dimension exceeds 4096");
    ComplexMatrix m(rows.size(),cols);
    for(std::size_t i=0;i<rows.size();++i){const auto& row=array(rows[i]);
        if(row.size()!=cols)throw std::invalid_argument("Matrix rows must have equal length");
        for(std::size_t j=0;j<cols;++j)m(i,j)=complexValue(row[j]);}
    return m;
}
inline J encodeMatrix(const ComplexMatrix& m){
    A rows;rows.reserve(m.rows);
    for(std::size_t i=0;i<m.rows;++i){A row;row.reserve(m.cols);for(std::size_t j=0;j<m.cols;++j)row.push_back(encode(m(i,j)));rows.emplace_back(std::move(row));}
    return rows;
}
inline SampledSignal signal(const J& q){
    SampledSignal g;g.times=reals(required(q,"times"));g.values=complexes(required(q,"values"));
    if(has(q,"derivatives"))g.derivatives=complexes(q.at("derivatives"));
    return g;
}
// Model poles: [{"eigenvalue":z,"laurent":[a1,a2,...]}].
inline PoleResidueModel poleModel(const J& v){
    PoleResidueModel m;
    for(const auto& p:array(v)){object(p);m.poles.push_back({complexValue(required(p,"eigenvalue")),complexes(required(p,"laurent"))});}
    validatePoleModel(m);
    return m;
}
inline J poleModelJson(const PoleResidueModel& m){
    A rows;for(const auto& p:m.poles)rows.emplace_back(O{{"eigenvalue",encode(p.eigenvalue)},{"laurent",encodeAll(p.laurent)}});
    return rows;
}
inline BilinearLiftedModel bilinearModel(const J& q){
    BilinearLiftedModel m;m.A=matrix(required(q,"A"));m.B=matrix(required(q,"B"));m.C=matrix(required(q,"C"));
    if(has(q,"N"))for(const auto& x:array(q.at("N")))m.N.push_back(matrix(x));
    validateBilinearModel(m);
    return m;
}
inline J latticeJson(const PoleLattice& l){
    A poles;poles.reserve(l.poles.size());
    for(const auto& p:l.poles)poles.emplace_back(O{{"value",encode(p.value)},{"k",numbers(p.k)},{"m",numbers(p.m)}});
    A coincide;for(const auto& c:l.coincidences)coincide.emplace_back(O{{"first",double(c.first)},{"second",double(c.second)},{"separation",c.separation}});
    return O{{"poles",poles},{"coincidences",coincide},{"smallest_harmonic_divisor",optionalNumber(l.smallestHarmonicDivisor)},{"scope",l.scope}};
}
}
