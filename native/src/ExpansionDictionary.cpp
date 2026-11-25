#include "PolynomialEvaluation.hpp"
#include <dynsyslab/core/BasisFunctions.hpp>
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/core/Observable.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {

std::size_t PolynomialDictionary::featureCount() const noexcept {
    return coordinateDimension ? exponents.size() / coordinateDimension : 0;
}

void PolynomialDictionary::validate() const {
    if (!coordinateDimension || exponents.empty() || exponents.size() % coordinateDimension ||
        center.size() != coordinateDimension || scale.size() != coordinateDimension)
        throw std::invalid_argument("Polynomial dictionary dimensions are inconsistent");
    if (!std::all_of(center.begin(), center.end(), [](double v) { return std::isfinite(v); }) ||
        !std::all_of(scale.begin(), scale.end(), [](double v) { return std::isfinite(v) && v > 0; }))
        throw std::invalid_argument("Polynomial normalization needs finite centers and positive scales");
    std::set<std::vector<std::size_t>> unique;
    for (std::size_t f = 0; f < featureCount(); ++f) {
        const auto begin = exponents.begin() + std::ptrdiff_t(f * coordinateDimension);
        if (!unique.emplace(begin, begin + std::ptrdiff_t(coordinateDimension)).second)
            throw std::invalid_argument("Duplicate polynomial monomial");
    }
}

Observable PolynomialDictionary::asObservable() const {
    validate();
    const auto d = coordinateDimension;
    std::vector<AnalyticExpression> z, terms;
    std::vector<std::string> names;
    for (std::size_t j = 0; j < d; ++j) z.push_back((AnalyticExpression::state(d, j) - center[j]) / scale[j]);
    for (std::size_t f = 0; f < featureCount(); ++f) {
        auto term = AnalyticExpression::constant(1, d);
        std::string name = "1";
        for (std::size_t j = 0; j < d; ++j)
            if (exponents[f * d + j]) {
                term = term * z[j].power(double(exponents[f * d + j]));
                name += " * z" + std::to_string(j + 1) + "^" + std::to_string(exponents[f * d + j]);
            }
        terms.push_back(term);
        names.push_back(name);
    }
    return Observable::fromExpressions({"polynomial_library", "Polynomial candidate library", "", "model", names, false}, terms);
}

namespace {

constexpr double kPi = 3.14159265358979323846;
const std::set<std::string> kProduct{"monomial", "chebyshev", "laguerre", "hermite", "fourier_hermite"};
const std::set<std::string> kRadial{"gaussian", "laplacian", "lorentzian", "matern"};
const std::set<std::string> kWave{"sine", "cosine", "exponential"};

double normalized(const DictionaryTerm& t, const double* q, std::size_t i) {
    return (q[i] - (t.center.empty() ? 0.0 : t.center[i])) / (t.scale.empty() ? 1.0 : t.scale[i]);
}
double normalizedNode(const DictionaryTerm& t, const double* node, std::size_t i) {
    return (node[i] - (t.center.empty() ? 0.0 : t.center[i])) / (t.scale.empty() ? 1.0 : t.scale[i]);
}
double distance(const DictionaryTerm& t, const double* q, const double* node, std::size_t d) {
    double r = 0;
    for (std::size_t i = 0; i < d; ++i) r = std::hypot(r, normalized(t, q, i) - normalizedNode(t, node, i));
    return r;
}

double value(const DictionaryTerm& t, const double* q, std::size_t d) {
    if (t.kind == "coordinate") return normalized(t, q, t.index);
    if (t.kind == "expression") return t.expression->evaluate(q);
    if (kProduct.count(t.kind)) {
        double v = 1;
        for (std::size_t i = 0; i < d; ++i) {
            const auto n = t.orders[i];
            const double x = normalized(t, q, i);
            if (t.kind == "monomial") { if (n) v *= std::pow(x, double(n)); }
            else if (t.kind == "chebyshev") v *= basis::chebyshev(n, x);
            else if (t.kind == "laguerre") v *= basis::laguerre(n, x);
            else if (t.kind == "hermite" || i == 1) v *= basis::hermiteFunction(n, x);
            else if (n == 0) v /= std::sqrt(2 * kPi);
            else { const auto k = double((n + 1) / 2); v *= (n % 2 ? std::cos(k * x) : std::sin(k * x)) / std::sqrt(kPi); }
        }
        return v;
    }
    if (kRadial.count(t.kind)) return basis::radialKernel(t.kind, distance(t, q, t.vector.data(), d), t.width, t.smoothness, t.amplitude);
    if (kWave.count(t.kind)) {
        double phase = t.phase;
        for (std::size_t i = 0; i < d; ++i) phase += t.vector[i] * normalized(t, q, i);
        return t.kind == "sine" ? std::sin(phase) : t.kind == "cosine" ? std::cos(phase) : std::exp(phase - t.phase);
    }
    // partition: nearest node (first on ties), the same distance as the spectral dictionary
    const std::size_t nodes = t.nodes.size() / d;
    double nearest = std::numeric_limits<double>::infinity();
    std::size_t best = 0;
    for (std::size_t j = 0; j < nodes; ++j) {
        const double r = distance(t, q, t.nodes.data() + j * d, d);
        if (r < nearest) { nearest = r; best = j; }
    }
    return best == t.index ? 1.0 : 0.0;
}

bool exact(const DictionaryTerm& t) {
    return t.kind == "monomial" || t.kind == "coordinate" || t.kind == "gaussian" || kWave.count(t.kind)
        || t.kind == "expression" || t.kind == "partition";
}

// d phi / d q for one term; finite differences for kinds without a closed form here.
std::vector<double> termGradient(const DictionaryTerm& t, const double* q, std::size_t d) {
    std::vector<double> g(d, 0.0);
    const auto scale = [&](std::size_t i) { return t.scale.empty() ? 1.0 : t.scale[i]; };
    if (t.kind == "coordinate") {
        g[t.index] = 1.0 / scale(t.index);
    } else if (t.kind == "expression") {
        for (std::size_t i = 0; i < d; ++i) g[i] = t.expression->derivative(i).evaluate(q);
    } else if (t.kind == "monomial") {
        for (std::size_t i = 0; i < d; ++i) {
            if (!t.orders[i]) continue;
            double v = double(t.orders[i]) * std::pow(normalized(t, q, i), double(t.orders[i] - 1)) / scale(i);
            for (std::size_t j = 0; j < d; ++j)
                if (j != i && t.orders[j]) v *= std::pow(normalized(t, q, j), double(t.orders[j]));
            g[i] = v;
        }
    } else if (t.kind == "gaussian") {
        const double phi = value(t, q, d);
        for (std::size_t i = 0; i < d; ++i)
            g[i] = -phi * (normalized(t, q, i) - normalizedNode(t, t.vector.data(), i)) / (t.width * t.width * scale(i));
    } else if (kWave.count(t.kind)) {
        double phase = t.phase;
        for (std::size_t i = 0; i < d; ++i) phase += t.vector[i] * normalized(t, q, i);
        const double outer = t.kind == "sine" ? std::cos(phase) : t.kind == "cosine" ? -std::sin(phase) : std::exp(phase - t.phase);
        for (std::size_t i = 0; i < d; ++i) g[i] = outer * t.vector[i] / scale(i);
    } else if (t.kind != "partition") {  // piecewise constant: zero almost everywhere
        std::vector<double> point(q, q + d);
        for (std::size_t i = 0; i < d; ++i) {
            const double step = 1e-6 * std::max(1.0, std::abs(q[i]));
            point[i] = q[i] + step;
            const double high = value(t, point.data(), d);
            point[i] = q[i] - step;
            const double low = value(t, point.data(), d);
            point[i] = q[i];
            g[i] = (high - low) / (2 * step);
        }
    }
    return g;
}

json::Value numbers(const std::vector<double>& v) {
    json::Value::Array a;
    for (double x : v) a.emplace_back(x);
    return a;
}
std::vector<double> doubles(const json::Value& v) {
    std::vector<double> out;
    if (!v.isNull()) for (const auto& x : v.array()) out.push_back(x.number());
    return out;
}

}  // namespace

ExpansionDictionary::ExpansionDictionary(std::size_t dimension, std::vector<DictionaryTerm> terms)
    : dimension_(dimension), terms_(std::move(terms)) {
    if (!dimension_ || terms_.empty()) throw std::invalid_argument("A dictionary needs a dimension and at least one term");
    for (const auto& t : terms_) {
        const bool known = t.kind == "expression" || t.kind == "partition" || t.kind == "coordinate" || kProduct.count(t.kind)
                        || kRadial.count(t.kind) || kWave.count(t.kind);
        if (!known) throw std::invalid_argument("Unknown dictionary term kind " + t.kind);
        if ((!t.center.empty() && t.center.size() != dimension_) || (!t.scale.empty() && t.scale.size() != dimension_))
            throw std::invalid_argument("Dictionary term normalization has the wrong dimension");
        for (double s : t.scale) if (!(s > 0) || !std::isfinite(s)) throw std::invalid_argument("Dictionary scales must be positive");
        if (kProduct.count(t.kind) && t.orders.size() != dimension_) throw std::invalid_argument("Product term needs one order per coordinate");
        if ((kRadial.count(t.kind) || kWave.count(t.kind)) && t.vector.size() != dimension_)
            throw std::invalid_argument(t.kind + " term needs a vector of the domain dimension");
        if (kRadial.count(t.kind) && !(t.width > 0)) throw std::invalid_argument("Radial term needs a positive width");
        if (t.kind == "coordinate" && t.index >= dimension_) throw std::invalid_argument("Coordinate term index is outside the domain");
        if (t.kind == "partition" && (t.nodes.empty() || t.nodes.size() % dimension_ || t.index >= t.nodes.size() / dimension_))
            throw std::invalid_argument("Partition term needs its node set and a valid index");
        if (t.kind == "expression" && (!t.expression || t.expression->stateDimension() != dimension_))
            throw std::invalid_argument("Expression term dimension differs from the dictionary");
    }
    // A table of monomials sharing one normalization keeps the established polynomial arithmetic.
    bool table = std::all_of(terms_.begin(), terms_.end(), [&](const DictionaryTerm& t) {
        return t.kind == "monomial" && t.center == terms_.front().center && t.scale == terms_.front().scale; });
    if (table) {
        PolynomialDictionary d;
        d.coordinateDimension = dimension_;
        for (const auto& t : terms_) d.exponents.insert(d.exponents.end(), t.orders.begin(), t.orders.end());
        d.center = terms_.front().center.empty() ? std::vector<double>(dimension_, 0.0) : terms_.front().center;
        d.scale = terms_.front().scale.empty() ? std::vector<double>(dimension_, 1.0) : terms_.front().scale;
        d.validate();
        table_ = std::move(d);
    }
}

ExpansionDictionary ExpansionDictionary::monomials(PolynomialDictionary table) {
    if (table.center.empty()) table.center.assign(table.coordinateDimension, 0.);
    if (table.scale.empty()) table.scale.assign(table.coordinateDimension, 1.);
    table.validate();
    std::vector<DictionaryTerm> terms;
    for (std::size_t f = 0; f < table.featureCount(); ++f) {
        DictionaryTerm t;
        t.kind = "monomial";
        t.center = table.center;
        t.scale = table.scale;
        t.orders.assign(table.exponents.begin() + std::ptrdiff_t(f * table.coordinateDimension),
                        table.exponents.begin() + std::ptrdiff_t((f + 1) * table.coordinateDimension));
        terms.push_back(std::move(t));
    }
    return ExpansionDictionary(table.coordinateDimension, std::move(terms));
}

std::vector<std::string> ExpansionDictionary::kinds() const {
    std::vector<std::string> out;
    for (const auto& t : terms_)
        if (std::find(out.begin(), out.end(), t.kind) == out.end()) out.push_back(t.kind);
    return out;
}

std::vector<double> ExpansionDictionary::evaluate(const double* q) const {
    std::vector<double> out;
    out.reserve(terms_.size());
    for (const auto& t : terms_) out.push_back(value(t, q, dimension_));
    for (double v : out) if (!std::isfinite(v)) throw std::invalid_argument("Dictionary evaluation overflow");
    return out;
}

std::vector<double> ExpansionDictionary::gradient(const double* q) const {
    std::vector<double> out;
    out.reserve(terms_.size() * dimension_);
    for (const auto& t : terms_) {
        const auto g = termGradient(t, q, dimension_);
        out.insert(out.end(), g.begin(), g.end());
    }
    return out;
}

bool ExpansionDictionary::exactGradient() const noexcept {
    return std::all_of(terms_.begin(), terms_.end(), [](const DictionaryTerm& t) { return exact(t); });
}

std::complex<double> ExpansionDictionary::combine(const std::complex<double>* c, const double* q) const {
    if (table_) return detail::polynomialValue(dimension_, table_->featureCount(), table_->exponents, table_->center, table_->scale, c, q);
    const auto phi = evaluate(q);
    std::complex<double> sum{};
    for (std::size_t k = 0; k < phi.size(); ++k) sum += c[k] * phi[k];
    return sum;
}

std::vector<std::complex<double>> ExpansionDictionary::combineGradient(const std::complex<double>* c, const double* q) const {
    if (table_) return detail::polynomialGradient(dimension_, table_->featureCount(), table_->exponents, table_->center, table_->scale, c, q);
    const auto g = gradient(q);
    std::vector<std::complex<double>> out(dimension_);
    for (std::size_t k = 0; k < terms_.size(); ++k)
        for (std::size_t i = 0; i < dimension_; ++i) out[i] += c[k] * g[k * dimension_ + i];
    return out;
}

json::Value ExpansionDictionary::describe() const {
    json::Value::Array terms;
    for (const auto& t : terms_) {
        json::Value::Object o{{"kind", t.kind}};
        if (!t.center.empty()) o["center"] = numbers(t.center);
        if (!t.scale.empty()) o["scale"] = numbers(t.scale);
        if (!t.orders.empty()) {
            json::Value::Array orders;
            for (auto n : t.orders) orders.emplace_back(double(n));
            o["orders"] = orders;
        }
        if (!t.vector.empty()) o["vector"] = numbers(t.vector);
        if (kRadial.count(t.kind)) { o["width"] = t.width; o["smoothness"] = t.smoothness; o["amplitude"] = t.amplitude; }
        if (t.kind == "sine" || t.kind == "cosine") o["phase"] = t.phase;
        if (t.kind == "partition") { o["nodes"] = numbers(t.nodes); o["index"] = double(t.index); }
        if (t.kind == "coordinate") o["index"] = double(t.index);
        if (t.kind == "expression") o["expression"] = t.expression->toJson();
        terms.emplace_back(std::move(o));
    }
    return json::Value::Object{{"schema", "dynsyslab.expansion-dictionary/v1"}, {"dimension", double(dimension_)},
                               {"terms", terms}};
}

ExpansionDictionary ExpansionDictionary::fromJson(const json::Value& value) {
    if (value.at("schema").string() != "dynsyslab.expansion-dictionary/v1")
        throw std::invalid_argument("Unknown expansion dictionary schema");
    std::vector<DictionaryTerm> terms;
    for (const auto& o : value.at("terms").array()) {
        DictionaryTerm t;
        t.kind = o.at("kind").string();
        t.center = doubles(o.at("center"));
        t.scale = doubles(o.at("scale"));
        if (!o.at("orders").isNull()) for (const auto& n : o.at("orders").array()) t.orders.push_back(std::size_t(n.number()));
        t.vector = doubles(o.at("vector"));
        t.width = o.at("width").number(1);
        t.smoothness = o.at("smoothness").number(2);
        t.amplitude = o.at("amplitude").number(1);
        t.phase = o.at("phase").number(0);
        t.nodes = doubles(o.at("nodes"));
        t.index = std::size_t(o.at("index").number(0));
        if (t.kind == "expression")
            t.expression = std::make_shared<const AnalyticExpression>(AnalyticExpression::fromJson(o.at("expression")));
        terms.push_back(std::move(t));
    }
    return ExpansionDictionary(std::size_t(value.at("dimension").number()), std::move(terms));
}

}  // namespace dynsyslab::core
