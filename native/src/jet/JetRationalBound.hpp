#pragma once
#include "JetPolynomial.hpp"

namespace dynsyslab::core::jet {

using Degree = boost::multiprecision::cpp_int;
struct RationalDegrees { Degree numerator, denominator; };

// Upper bounds on the polynomial degrees of one (possibly unreduced)
// fraction. Retaining the original domain allows common factors and poles.
inline RationalDegrees rationalDegrees(const json::Value& node) {
    const auto op=node.at("operation").string();
    if (op=="state") return {1,0};
    if (op=="constant") return {0,0};
    const auto& operands=node.at("arguments").array();
    auto a=rationalDegrees(operands.at(0));
    if (op=="negative") return a;
    if (op=="power") {
        Degree exponent(node.at("scalar").number());
        if (exponent<0) {
            std::swap(a.numerator,a.denominator);
            exponent=-exponent;
        }
        return {a.numerator*exponent,a.denominator*exponent};
    }
    auto b=rationalDegrees(operands.at(1));
    if (op=="divide") std::swap(b.numerator,b.denominator);
    if (op=="divide" || op=="multiply")
        return {a.numerator+b.numerator,a.denominator+b.denominator};
    return {std::max(Degree(a.numerator+b.denominator),Degree(b.numerator+a.denominator)),
        a.denominator+b.denominator};
}

// Degree bounds for the common denominator-cleared field on two copies of
// the state. For component i we multiply its numerator by every denominator
// of both copies except its own. Each output difference clears its own pair.
inline json::Value rationalOrderBound(const std::vector<json::Value>& field,
    const std::vector<json::Value>& observable) {
    const int m=int(2*field.size());
    std::vector<RationalDegrees> rates;
    Degree sum=0,delta=1,d=m;
    for (const auto& node:field) {
        rates.push_back(rationalDegrees(node));
        sum+=rates.back().denominator;
    }
    for (const auto& rate:rates)
        delta=std::max(delta,Degree(2*sum-rate.denominator+rate.numerator));
    for (const auto& node:observable) {
        const auto p=rationalDegrees(node);
        d=std::max(d,Degree(p.numerator+p.denominator));
    }
    const Degree base=d+(m-1)*(delta-1);
    Degree bound=Degree(1)<<(m+1);
    for (int i=0;i<m;++i) bound*=base;
    return json::Value::Object{{"theorem","Binyamini (2016), Corollary 1"},
        {"chart","original_rational_expressions"},{"pair_dimension",double(m)},
        {"field_degree_bound_exact",delta.str()},{"output_degree_bound_exact",d.str()},
        {"sufficient_order",bound.str()}};
}
} // namespace dynsyslab::core::jet
