// Conversions between analytic expressions and dictionary terms.
//  * FunctionExpansion::fromExpressions decomposes each expression into terms: polynomial parts
//    become monomials (so x2 + c*x1^2/(b-2a) and a fitted SINDy bank share terms), and the remaining
//    pieces of a linear combination become expression terms identified by value.
//  * dictionaryTermExpression writes a term in closed form when one exists (every kind except
//    Matern kernels and partitions), which lets an expansion become an expression Observable.
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Monomials = std::map<std::vector<std::size_t>, double>;

std::optional<Monomials> polynomial(const json::Value& n, std::size_t d) {
    const auto op = n.at("operation").string();
    if (op == "constant") return Monomials{{std::vector<std::size_t>(d, 0), n.at("scalar").number()}};
    if (op == "state") {
        std::vector<std::size_t> e(d, 0);
        e[std::size_t(n.at("state_index").number())] = 1;
        return Monomials{{e, 1.0}};
    }
    const auto& args = n.at("arguments").isNull() ? json::Value::Array{} : n.at("arguments").array();
    if (op == "negative") {
        auto a = polynomial(args[0], d);
        if (a) for (auto& [e, c] : *a) c = -c;
        return a;
    }
    if (op == "add" || op == "subtract" || op == "multiply" || op == "divide") {
        auto a = polynomial(args[0], d), b = polynomial(args[1], d);
        if (!a || !b) return std::nullopt;
        Monomials out;
        if (op == "add" || op == "subtract") {
            out = *a;
            for (const auto& [e, c] : *b) out[e] += op == "add" ? c : -c;
        } else if (op == "multiply") {
            for (const auto& [ea, ca] : *a)
                for (const auto& [eb, cb] : *b) {
                    auto e = ea;
                    for (std::size_t i = 0; i < d; ++i) e[i] += eb[i];
                    out[e] += ca * cb;
                }
        } else {  // divide: only by a constant
            if (b->size() != 1 || b->begin()->first != std::vector<std::size_t>(d, 0) || b->begin()->second == 0)
                return std::nullopt;
            for (const auto& [e, c] : *a) out[e] = c / b->begin()->second;
        }
        return out;
    }
    if (op == "power") {
        const double p = n.at("scalar").number();
        if (p < 0 || p != std::floor(p) || p > 64) return std::nullopt;
        auto base = polynomial(args[0], d);
        if (!base) return std::nullopt;
        Monomials out{{std::vector<std::size_t>(d, 0), 1.0}};
        for (int k = 0; k < int(p); ++k) {
            Monomials next;
            for (const auto& [ea, ca] : out)
                for (const auto& [eb, cb] : *base) {
                    auto e = ea;
                    for (std::size_t i = 0; i < d; ++i) e[i] += eb[i];
                    next[e] += ca * cb;
                }
            out = std::move(next);
        }
        return out;
    }
    return std::nullopt;
}

// Linear combination of monomials and expression atoms: atom key -> (coefficient, atom json).
struct Combination {
    Monomials monomials;
    std::map<std::string, std::pair<double, json::Value>> atoms;
    void add(const Combination& other, double scale) {
        for (const auto& [e, c] : other.monomials) monomials[e] += scale * c;
        for (const auto& [k, v] : other.atoms) {
            auto& slot = atoms[k];
            slot.first += scale * v.first;
            slot.second = v.second;
        }
    }
};

Combination combination(const json::Value& n, std::size_t d) {
    if (auto p = polynomial(n, d)) return Combination{*p, {}};
    const auto op = n.at("operation").string();
    const auto& args = n.at("arguments").isNull() ? json::Value::Array{} : n.at("arguments").array();
    Combination out;
    if (op == "add" || op == "subtract") {
        out = combination(args[0], d);
        out.add(combination(args[1], d), op == "add" ? 1.0 : -1.0);
        return out;
    }
    if (op == "negative") { out.add(combination(args[0], d), -1.0); return out; }
    if (op == "multiply" || op == "divide") {  // a constant factor scales the other side
        const auto a = polynomial(args[0], d), b = polynomial(args[1], d);
        const auto constant = [&](const std::optional<Monomials>& m) {
            return m && m->size() == 1 && m->begin()->first == std::vector<std::size_t>(d, 0);
        };
        if (op == "multiply" && constant(a)) { out.add(combination(args[1], d), a->begin()->second); return out; }
        if (constant(b) && b->begin()->second != 0) {
            out.add(combination(args[0], d), op == "multiply" ? b->begin()->second : 1.0 / b->begin()->second);
            return out;
        }
    }
    out.atoms[json::dump(n)] = {1.0, n};
    return out;
}

}  // namespace

FunctionExpansion FunctionExpansion::fromExpressions(const std::vector<AnalyticExpression>& expressions,
                                                         FunctionDomain domain, std::vector<std::string> names) {
    if (expressions.empty()) throw std::invalid_argument("fromExpressions needs at least one expression");
    const std::size_t d = domain.dimension;
    std::vector<Combination> parts;
    std::map<std::vector<std::size_t>, std::size_t> monomialIndex;
    std::map<std::string, std::size_t> atomIndex;
    std::vector<DictionaryTerm> terms;
    for (const auto& expression : expressions) {
        if (expression.stateDimension() != d) throw std::invalid_argument("expression and domain dimensions differ");
        parts.push_back(combination(expression.toJson(), d));
        for (const auto& [e, c] : parts.back().monomials)
            if (c != 0 && !monomialIndex.count(e)) {
                monomialIndex[e] = terms.size();
                DictionaryTerm t;
                t.kind = "monomial";
                t.orders = e;
                terms.push_back(std::move(t));
            }
        for (const auto& [k, v] : parts.back().atoms)
            if (v.first != 0 && !atomIndex.count(k)) {
                atomIndex[k] = terms.size();
                DictionaryTerm t;
                t.kind = "expression";
                t.expression = std::make_shared<const AnalyticExpression>(AnalyticExpression::fromJson(v.second));
                terms.push_back(std::move(t));
            }
    }
    if (terms.empty()) {  // identically zero: one constant term with zero coefficient
        DictionaryTerm t;
        t.kind = "monomial";
        t.orders.assign(d, 0);
        terms.push_back(std::move(t));
        monomialIndex[std::vector<std::size_t>(d, 0)] = 0;
    }
    const std::size_t n = terms.size();
    std::vector<std::complex<double>> coefficients(expressions.size() * n);
    for (std::size_t r = 0; r < parts.size(); ++r) {
        for (const auto& [e, c] : parts[r].monomials) if (c != 0) coefficients[r * n + monomialIndex.at(e)] = c;
        for (const auto& [k, v] : parts[r].atoms) if (v.first != 0) coefficients[r * n + atomIndex.at(k)] = v.first;
    }
    return FunctionExpansion(std::make_shared<const ExpansionDictionary>(d, std::move(terms)), std::move(coefficients),
                               std::move(domain), std::move(names));
}

namespace {

AnalyticExpression scaled(const DictionaryTerm& t, std::size_t i, std::size_t d) {
    auto x = AnalyticExpression::state(d, i);
    if (!t.center.empty() && t.center[i] != 0) x = x - t.center[i];
    if (!t.scale.empty() && t.scale[i] != 1) x = x / t.scale[i];
    return x;
}

AnalyticExpression product(const DictionaryTerm& t, std::size_t d) {
    constexpr double pi = 3.14159265358979323846;
    auto v = AnalyticExpression::constant(1, d);
    for (std::size_t i = 0; i < d; ++i) {
        const auto n = t.orders[i];
        const auto x = scaled(t, i, d);
        if (t.kind == "monomial") { if (n) v = v * x.power(double(n)); continue; }
        if (t.kind == "chebyshev" || t.kind == "laguerre") {
            auto a = AnalyticExpression::constant(1, d);
            auto b = t.kind == "chebyshev" ? x : 1.0 - x;
            for (std::size_t j = 2; j <= n; ++j) {
                auto c = t.kind == "chebyshev" ? 2.0 * x * b - a : ((double(2 * j - 1) - x) * b - double(j - 1) * a) / double(j);
                a = b;
                b = c;
            }
            v = v * (n ? b : a);
            continue;
        }
        if (t.kind == "hermite" || i == 1) {  // Hermite function: polynomial recurrence times exp(-x^2/2)
            auto a = AnalyticExpression::constant(std::pow(pi, -0.25), d);
            auto b = std::sqrt(2.0) * x * a;
            for (std::size_t j = 2; j <= n; ++j) {
                auto c = std::sqrt(2.0 / double(j)) * x * b - std::sqrt(double(j - 1) / double(j)) * a;
                a = b;
                b = c;
            }
            v = v * (n ? b : a) * (-0.5 * x.power(2)).unary("exp");
            continue;
        }
        if (n == 0) { v = v / std::sqrt(2 * pi); continue; }  // Fourier part of fourier_hermite
        const double k = double((n + 1) / 2);
        v = v * (k * x).unary(n % 2 ? "cos" : "sin") / std::sqrt(pi);
    }
    return v;
}

}  // namespace

std::optional<AnalyticExpression> dictionaryTermExpression(const DictionaryTerm& t, std::size_t d) {
    if (t.kind == "expression") return *t.expression;
    if (t.kind == "monomial" || t.kind == "chebyshev" || t.kind == "laguerre" || t.kind == "hermite" || t.kind == "fourier_hermite")
        return product(t, d);
    if (t.kind == "sine" || t.kind == "cosine" || t.kind == "exponential") {
        auto phase = AnalyticExpression::constant(t.kind == "exponential" ? 0.0 : t.phase, d);
        for (std::size_t i = 0; i < d; ++i) if (t.vector[i] != 0) phase = phase + t.vector[i] * scaled(t, i, d);
        return phase.unary(t.kind == "sine" ? "sin" : t.kind == "cosine" ? "cos" : "exp");
    }
    if (t.kind == "gaussian" || t.kind == "laplacian" || t.kind == "lorentzian") {
        auto r2 = AnalyticExpression::constant(0, d);
        for (std::size_t i = 0; i < d; ++i) {
            const double node = (t.vector[i] - (t.center.empty() ? 0.0 : t.center[i])) / (t.scale.empty() ? 1.0 : t.scale[i]);
            r2 = r2 + (scaled(t, i, d) - node).power(2);
        }
        const auto z2 = r2 / (t.width * t.width);
        if (t.kind == "gaussian") return t.amplitude * (-0.5 * z2).unary("exp");
        if (t.kind == "lorentzian") return t.amplitude / (1.0 + z2);
        return t.amplitude * (-(z2.power(0.5))).unary("exp");
    }
    return std::nullopt;  // Matern kernels and partitions have no closed form here
}

std::optional<std::vector<AnalyticExpression>> FunctionExpansion::realExpressions() const {
    const std::size_t d = coordinateDimension(), n = dictionary().size();
    std::vector<std::optional<AnalyticExpression>> terms(n);
    std::vector<AnalyticExpression> out;
    for (std::size_t r = 0; r < outputCount(); ++r) {
        auto re = AnalyticExpression::constant(0, d), im = AnalyticExpression::constant(0, d);
        for (std::size_t k = 0; k < n; ++k) {
            const auto c = coefficients()[r * n + k];
            if (c == std::complex<double>{}) continue;
            if (!terms[k]) terms[k] = dictionaryTermExpression(dictionary().terms()[k], d);
            if (!terms[k]) return std::nullopt;
            if (c.real() != 0) re = re + c.real() * *terms[k];
            if (c.imag() != 0) im = im + c.imag() * *terms[k];
        }
        out.push_back(re);
        out.push_back(im);
    }
    return out;
}

FunctionExpansion FunctionExpansion::canonical() const {
    const std::size_t d = coordinateDimension(), n = dictionary().size(), outputs = outputCount();
    // Each term maps to canonical terms with weights; normalized monomials expand binomially.
    std::vector<DictionaryTerm> terms;
    std::map<std::string, std::size_t> index;
    std::vector<std::vector<std::pair<std::size_t, double>>> map(n);
    const auto keyOf = [&](const DictionaryTerm& t) {
        return json::dump(ExpansionDictionary(d, {t}).describe().at("terms").array().front());
    };
    const auto place = [&](const DictionaryTerm& t) {
        const auto key = keyOf(t);
        const auto found = index.find(key);
        if (found != index.end()) return found->second;
        index[key] = terms.size();
        terms.push_back(t);
        return terms.size() - 1;
    };
    for (std::size_t k = 0; k < n; ++k) {
        const auto& t = dictionary().terms()[k];
        const bool normalized = std::any_of(t.center.begin(), t.center.end(), [](double v) { return v != 0; })
                             || std::any_of(t.scale.begin(), t.scale.end(), [](double v) { return v != 1; });
        if (t.kind != "monomial" || !normalized) {
            auto plain = t;
            if (t.kind == "monomial") { plain.center.clear(); plain.scale.clear(); }
            map[k].push_back({place(plain), 1.0});
            continue;
        }
        const auto expansion = polynomial(product(t, d).toJson(), d);
        for (const auto& [e, c] : *expansion) {
            if (c == 0) continue;
            DictionaryTerm m;
            m.kind = "monomial";
            m.orders = e;
            map[k].push_back({place(m), c});
        }
    }
    const std::size_t m = terms.size();
    std::vector<std::complex<double>> c(outputs * m);
    for (std::size_t r = 0; r < outputs; ++r)
        for (std::size_t k = 0; k < n; ++k)
            for (const auto& [j, w] : map[k]) c[r * m + j] += coefficients()[r * n + k] * w;
    return FunctionExpansion(std::make_shared<const ExpansionDictionary>(d, std::move(terms)), std::move(c), coordinates(),
                               outputNames());
}

}  // namespace dynsyslab::core
