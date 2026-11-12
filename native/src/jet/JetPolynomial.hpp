#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <dynsyslab/core/Json.hpp>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <cmath>

namespace dynsyslab::core::jet {

using Rational = boost::multiprecision::cpp_rational;
using Powers = std::vector<int>;
using Polynomial = std::map<Powers, Rational>;
struct SearchLimit : std::runtime_error { using std::runtime_error::runtime_error; };
struct Unsupported : std::runtime_error { using std::runtime_error::runtime_error; };

inline Powers unit(std::size_t n, std::size_t axis) {
    Powers powers(n);
    powers.at(axis) = 1;
    return powers;
}

inline void addTerm(Polynomial& result, const Powers& key, const Rational& coefficient) {
    if (coefficient == 0) return;
    auto& value = result[key];
    value += coefficient;
    if (value == 0) result.erase(key);
    if (result.size() > 4096) throw SearchLimit("sparse polynomial term limit");
}

inline void addScaled(Polynomial& result, const Polynomial& p, const Rational& factor) {
    for (const auto& [powers, coefficient] : p) addTerm(result, powers, factor * coefficient);
}

inline Polynomial multiply(const Polynomial& a, const Polynomial& b) {
    Polynomial result;
    if (a.size() * b.size() > 100000) throw SearchLimit("polynomial product limit");
    for (const auto& [p, c] : a) for (const auto& [q, d] : b) {
        Powers key = p;
        for (std::size_t i = 0; i < key.size(); ++i) {
            key[i] += q[i];
            if (std::abs(key[i]) > 256) throw SearchLimit("monomial exponent limit");
        }
        addTerm(result, key, c*d);
    }
    return result;
}

inline Polynomial power(Polynomial base, int exponent, std::size_t n) {
    if (exponent < 0) {
        if (base.empty()) throw std::invalid_argument("denominator is identically zero");
        if (base.size() != 1) throw Unsupported("denominator has several monomials");
        auto [powers, coefficient] = *base.begin();
        Powers inverse = powers;
        for (int& entry : inverse) entry = -entry;
        base = {{inverse, 1/coefficient}};
        exponent = -exponent;
    }
    Polynomial result{{Powers(n), 1}};
    while (exponent) {
        if (exponent % 2) result = multiply(result, base);
        exponent /= 2;
        if (exponent) base = multiply(base, base);
    }
    return result;
}

inline Polynomial parseLaurent(const json::Value& expr, std::set<std::size_t>& poles) {
    const auto op = expr.at("operation").string();
    const auto n = std::size_t(expr.at("state_dimension").number());
    if (op == "state") return {{unit(n, std::size_t(expr.at("state_index").number())), 1}};
    if (op == "constant") {
        Polynomial result;
        addTerm(result, Powers(n), Rational(expr.at("scalar").number()));
        return result;
    }
    const auto& args = expr.at("arguments").array();
    auto a = parseLaurent(args.at(0), poles);
    if (op == "negative") { for (auto& term : a) term.second = -term.second; return a; }
    auto guardedPower = [&](Polynomial base, int exponent) {
        if (exponent < 0 && base.size() == 1)
            for (std::size_t i = 0; i < n; ++i)
                if (base.begin()->first[i]) poles.insert(i);
        return power(std::move(base), exponent, n);
    };
    if (op == "power") {
        const double exponent = expr.at("scalar").number();
        if (std::floor(exponent) != exponent) throw Unsupported("noninteger exponent");
        if (std::abs(exponent) > 64) throw SearchLimit("input exponent limit");
        return guardedPower(std::move(a), int(exponent));
    }
    if (args.size() != 2) throw Unsupported("nonrational operation: " + op);
    auto b = parseLaurent(args[1], poles);
    if (op == "multiply") return multiply(a,b);
    if (op == "divide") return multiply(a, guardedPower(std::move(b),-1));
    if (op!="add" && op!="subtract") throw Unsupported("unsupported binary operation: "+op);
    addScaled(a,b,op == "subtract" ? -1 : 1);
    return a;
}

inline Polynomial partial(const Polynomial& p, std::size_t axis) {
    Polynomial result;
    for (const auto& [powers,c] : p) if (powers[axis]) {
        auto key = powers;
        --key[axis];
        addTerm(result, key, c * powers[axis]);
    }
    return result;
}

inline Polynomial lie(const Polynomial& p, const std::vector<Polynomial>& field) {
    Polynomial result;
    for (std::size_t i = 0; i < field.size(); ++i)
        addScaled(result, multiply(partial(p,i),field[i]), 1);
    return result;
}

// Multivariate division by known graph identities. Every subtraction
// preserves values on the graph. A nonzero remainder is inconclusive.
inline Polynomial onGraph(Polynomial p,const std::vector<Polynomial>* graph) {
    if (!graph || graph->empty()) return p;
    Polynomial remainder;
    std::size_t steps=0;
    while (!p.empty()) {
        if (++steps>8192) throw SearchLimit("graph reduction limit");
        const auto [key,c]=*p.rbegin();
        bool reduced=false;
        for (const auto& g:*graph) {
            if (g.empty()) continue;
            const auto& [lead,a]=*g.rbegin();
            Powers difference=key;
            bool divisible=true;
            for (std::size_t i=0;i<key.size();++i) {
                difference[i]-=lead[i];
                if (difference[i]<0) {divisible=false; break;}
            }
            if (!divisible) continue;
            addScaled(p,multiply({{difference,1}},g),-c/a);
            reduced=true; break;
        }
        if (!reduced) {addTerm(remainder,key,c); p.erase(key);}
    }
    return remainder;
}

inline json::Value encode(const Polynomial& p) {
    json::Value::Array terms;
    for (const auto& [powers,c] : p) {
        json::Value::Array exponents;
        for (int exponent : powers) exponents.emplace_back(double(exponent));
        terms.emplace_back(json::Value::Object{{"coefficient",c.str()},{"powers",std::move(exponents)}});
    }
    return terms;
}

} // namespace dynsyslab::core::jet
