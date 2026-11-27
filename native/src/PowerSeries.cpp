// Truncated multivariate power series (PowerSeries.hpp).
#include "PowerSeries.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {

void enumerate(std::size_t variables, std::size_t remaining, std::size_t at, std::vector<std::size_t>& current,
               std::vector<std::size_t>& out) {
    if (at + 1 == variables) {
        current[at] = remaining;
        out.insert(out.end(), current.begin(), current.end());
        return;
    }
    for (std::size_t e = remaining + 1; e-- > 0;) {
        current[at] = e;
        enumerate(variables, remaining - e, at + 1, current, out);
    }
}

// sum_j a_j u~^j by Horner, u~ = u minus its constant term.
Series compose(const Series& u, const std::vector<Complex>& a) {
    auto tail = u;
    tail.c[0] = 0;
    Series out(u.shared(), a.back());
    for (std::size_t j = a.size() - 1; j-- > 0;) {
        out = out * tail;
        out.c[0] += a[j];
    }
    return out;
}

}  // namespace

MonomialTable::MonomialTable(std::size_t variables, std::size_t order) : variables_(variables), order_(order) {
    if (!variables) throw std::invalid_argument("power series need at least one variable");
    std::vector<std::size_t> current(variables);
    for (std::size_t k = 0; k <= order; ++k) {
        starts_.push_back(exponents_.size() / variables);
        enumerate(variables, k, 0, current, exponents_);
        if (exponents_.size() / variables > 2500) throw std::invalid_argument("power series: more than 2500 monomials");
    }
    starts_.push_back(exponents_.size() / variables);
    const auto m = exponents_.size() / variables;
    for (std::size_t k = 0; k <= order; ++k)
        for (std::size_t i = starts_[k]; i < starts_[k + 1]; ++i) degrees_.push_back(k);
    products_.assign(m * m, m);
    std::map<std::vector<std::size_t>, std::size_t> lookup;
    for (std::size_t i = 0; i < m; ++i) lookup.emplace(std::vector<std::size_t>(exponents(i), exponents(i) + variables), i);
    std::vector<std::size_t> sum(variables);
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m && degrees_[a] + degrees_[b] <= order; ++b) {
            for (std::size_t v = 0; v < variables; ++v) sum[v] = exponents(a)[v] + exponents(b)[v];
            products_[a * m + b] = lookup.at(sum);
        }
}

std::size_t MonomialTable::index(const std::size_t* e) const {
    std::size_t k = 0;
    for (std::size_t v = 0; v < variables_; ++v) k += e[v];
    if (k > order_) return size();
    for (std::size_t i = starts_[k]; i < starts_[k + 1]; ++i)
        if (std::equal(e, e + variables_, exponents(i))) return i;
    return size();
}

Series::Series(std::shared_ptr<const MonomialTable> table, Complex constant) : c(table->size(), 0), table_(std::move(table)) {
    c[0] = constant;
}

Series Series::variable(std::shared_ptr<const MonomialTable> table, std::size_t index, Complex at) {
    Series out(table, at);
    if (table->order() >= 1) out.c[1 + index] = 1;  // degree-1 monomials follow the constant, variable order
    return out;
}

Series Series::operator+(const Series& b) const { auto o = *this; for (std::size_t i = 0; i < c.size(); ++i) o.c[i] += b.c[i]; return o; }
Series Series::operator-(const Series& b) const { auto o = *this; for (std::size_t i = 0; i < c.size(); ++i) o.c[i] -= b.c[i]; return o; }
Series Series::operator*(Complex s) const { auto o = *this; for (auto& v : o.c) v *= s; return o; }
Series Series::operator-() const { return *this * Complex(-1); }

Series Series::operator*(const Series& b) const {
    Series o(table_);
    const auto m = c.size();
    const auto& t = *table_;
    for (std::size_t i = 0; i < m; ++i) {
        if (c[i] == Complex(0)) continue;
        for (std::size_t j = 0; j < m && t.degree(i) + t.degree(j) <= t.order(); ++j)
            if (b.c[j] != Complex(0)) o.c[t.product(i, j)] += c[i] * b.c[j];
    }
    return o;
}

Series Series::truncated(std::size_t k) const {
    auto o = *this;
    for (std::size_t i = 0; i < c.size(); ++i) if (table_->degree(i) >= k) o.c[i] = 0;
    return o;
}

Series exp(const Series& u) {
    std::vector<Complex> a(u.table().order() + 1);
    Complex term = std::exp(u.c[0]);
    for (std::size_t j = 0; j < a.size(); ++j) { a[j] = term; term /= double(j + 1); }
    return compose(u, a);
}

Series log(const Series& u) {
    if (u.c[0] == Complex(0)) throw std::domain_error("power series: log at zero");
    std::vector<Complex> a(u.table().order() + 1);
    a[0] = std::log(u.c[0]);
    Complex inverse = 1.0 / u.c[0], power = inverse;
    for (std::size_t j = 1; j < a.size(); ++j) { a[j] = (j % 2 ? 1.0 : -1.0) * power / double(j); power *= inverse; }
    return compose(u, a);
}

namespace {
Series trigonometric(const Series& u, bool cosine) {
    const Complex s = std::sin(u.c[0]), c = std::cos(u.c[0]);
    // Derivatives of sin: sin, cos, -sin, -cos; of cos: cos, -sin, -cos, sin.
    const Complex cycle[4] = {cosine ? c : s, cosine ? -s : c, cosine ? -c : -s, cosine ? s : -c};
    std::vector<Complex> a(u.table().order() + 1);
    double factorial = 1;
    for (std::size_t j = 0; j < a.size(); ++j) { if (j) factorial *= double(j); a[j] = cycle[j % 4] / factorial; }
    return compose(u, a);
}
}  // namespace

Series sin(const Series& u) { return trigonometric(u, false); }
Series cos(const Series& u) { return trigonometric(u, true); }

Series power(const Series& u, double p) {
    if (p >= 0 && p == std::floor(p) && p <= 64) {  // exact by repeated squaring (also at zero)
        auto n = std::size_t(p);
        Series out(u.shared(), 1), base = u;
        while (n) { if (n & 1) out = out * base; n >>= 1; if (n) base = base * base; }
        return out;
    }
    if (u.c[0] == Complex(0)) throw std::domain_error("power series: non-integer power at zero");
    std::vector<Complex> a(u.table().order() + 1);
    Complex binomial = 1;
    for (std::size_t j = 0; j < a.size(); ++j) {
        a[j] = binomial * std::pow(u.c[0], p - double(j));
        binomial *= (p - double(j)) / double(j + 1);
    }
    return compose(u, a);
}

Series reciprocal(const Series& u) { return power(u, -1); }

Series evaluate(const json::Value& e, const std::vector<Series>& state) {
    const auto op = e.at("operation").string();
    const auto& table = state.at(0).shared();
    if (op == "constant") return Series(table, e.at("scalar").number());
    if (op == "state") return state.at(std::size_t(e.at("state_index").number()));
    const auto& args = e.at("arguments").array();
    const auto a = evaluate(args.at(0), state);
    if (op == "negative") return -a;
    if (op == "power") return power(a, e.at("scalar").number());
    if (op == "exp") return exp(a);
    if (op == "log") return log(a);
    if (op == "sin") return sin(a);
    if (op == "cos") return cos(a);
    if (op == "tanh") { const auto e2 = exp(a * Complex(2)); return (e2 - Series(table, 1)) * reciprocal(e2 + Series(table, 1)); }
    if (op == "sech") return reciprocal(exp(a) + exp(-a)) * Complex(2);
    if (args.size() < 2) throw std::invalid_argument("power series: unsupported expression " + op);
    const auto b = evaluate(args[1], state);
    if (op == "add") return a + b;
    if (op == "subtract") return a - b;
    if (op == "multiply") return a * b;
    if (op == "divide") return a * reciprocal(b);
    throw std::invalid_argument("power series: no series for the expression function " + op);
}

}  // namespace dynsyslab::core::detail
