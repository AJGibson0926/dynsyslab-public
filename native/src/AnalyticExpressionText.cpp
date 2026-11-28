// Infix text for analytic expressions, e.g. "r^2 * cos(theta) - 0.5*exp(-q1)".
// Grammar: text := (name '=' expr ';')* expr, where each definition names an intermediate that later
// text uses by name (shared, not copied); expr := term (('+'|'-') term)*; term := power (('*'|'/') power)*;
// power := unary ('^' signed-number)?; unary := '-' unary | primary;
// primary := number | variable | function '(' expr ')' | '(' expr ')'.
// Functions: exp, log, sin, cos, tanh, sech, sqrt, atan2(y, x), and the elliptic integrals in the
// parameter m = k^2: ellipk(m) = K(m), ellipe(m) = E(m), ellipf(phi, m) = F(phi|m),
// ellipe(phi, m) = E(phi|m), ellippi(phi, n, m) = Pi(n; phi|m) (n < 1); and polyroot(k, c_0, ..., c_d),
// the k-th largest real root of c_0 x^d + ... + c_d (k = 0 the largest). An exponent may be a number,
// an expression free of the state (a summation variable: the power function), or any expression
// (then exp(exponent log base)).
// Special functions (ExpressionFunctions.cpp) take their constant parameters first: besselj(nu, x),
// bessely, besseli, besselk, sphbesselj(n, x), sphbessely, airyai(x), airybi, airyaiprime, airybiprime,
// gamma, lgamma, digamma, polygamma(n, x), beta(a, b), erf, erfc, expint (Ei), expintn(n, x) (E_n),
// gammap(a, x), gammaq(a, x), lambertw, lambertwm1, zeta(s), hyp0f1(b, x), hyp1f1(a, b, x),
// hyp2f1(a, b, c, x), jacobiam(u, m), weierstrassp(g2, g3, z), weierstrasspprime(g2, g3, z), and tan,
// asin, acos, atan, sinh, cosh, asinh, acosh, atanh, abs. Built from those: sn(u, m), cn(u, m),
// dn(u, m), jacobizeta(phi, m), theta1..theta4(z, q) (series in the nome q), legendre(n, x),
// hermite(n, x), laguerre(n, x), chebyshevt(n, x), chebyshevu(n, x).
// Series: sum(k, from, to, term) with k bound in term and to an integer or inf, e.g.
// sum(k, 0, inf, x^k/gamma(k + 1)) (an infinite sum checks its own convergence; AnalyticExpression.hpp).
#include <dynsyslab/core/AnalyticExpression.hpp>
#include "ExpressionFunctions.hpp"

#include <cmath>

#include <cctype>
#include <cstdlib>
#include <map>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

class Parser {
public:
    Parser(const std::string& text, const std::vector<std::string>& variables, const std::map<std::string, double>& constants)
        : text_(text), variables_(variables), constants_(constants) {}

    AnalyticExpression parse() {
        for (;;) {  // definitions: name '=' expr ';'
            skip();
            const auto start = position_;
            std::string name;
            while (position_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[position_])) || text_[position_] == '_'))
                name.push_back(text_[position_++]);
            if (name.empty() || std::isdigit(static_cast<unsigned char>(name.front())) || !take('=')) { position_ = start; break; }
            for (const auto& variable : variables_)
                if (variable == name) fail("definition " + name + " shadows a variable");
            auto value = expression();
            if (!take(';')) fail("expected ';' after the definition of " + name);
            definitions_.insert_or_assign(name, std::move(value));
        }
        auto value = expression();
        skip();
        if (position_ != text_.size()) fail("unexpected '" + std::string(1, text_[position_]) + "'");
        return value;
    }

private:
    [[noreturn]] void fail(const std::string& why) const {
        throw std::invalid_argument("expression \"" + text_ + "\": " + why + " at position " + std::to_string(position_));
    }
    void skip() { while (position_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[position_]))) ++position_; }
    bool take(char c) {
        skip();
        if (position_ < text_.size() && text_[position_] == c) { ++position_; return true; }
        return false;
    }
    double number() {
        skip();
        const char* begin = text_.c_str() + position_;
        char* end = nullptr;
        const double value = std::strtod(begin, &end);
        if (end == begin) fail("expected a number");
        position_ += std::size_t(end - begin);
        return value;
    }
    AnalyticExpression expression() {
        auto value = term();
        for (;;) {
            if (take('+')) value = value + term();
            else if (take('-')) value = value - term();
            else return value;
        }
    }
    AnalyticExpression term() {
        auto value = power();
        for (;;) {
            if (take('*')) value = value * power();
            else if (take('/')) value = value / power();
            else return value;
        }
    }
    // base ^ exponent: a number is a constant power; an exponent free of the state (a summation
    // variable, a constant expression) is the power function pow; one that depends on the state is
    // exp(exponent log base).
    AnalyticExpression power() {
        auto base = unary();
        if (!take('^')) return base;
        skip();
        const auto start = position_;
        const bool negative = take('-');
        if (!negative) take('+');
        skip();
        if (position_ < text_.size() && (std::isdigit(static_cast<unsigned char>(text_[position_])) || text_[position_] == '.')) {
            const double exponent = number();
            skip();
            if (position_ >= text_.size() || std::string("+-*/),;").find(text_[position_]) != std::string::npos)
                return base.power(negative ? -exponent : exponent);
        }
        position_ = start;
        const auto exponent = unary();
        if (exponent.dependsOnState()) return (exponent * base.unary("log")).unary("exp");
        return AnalyticExpression::special("pow", std::vector<AnalyticExpression>{exponent}, {base});
    }
    AnalyticExpression unary() {
        if (take('-')) return -unary();
        if (take('+')) return unary();
        return primary();
    }
    double constantValue(const AnalyticExpression& e, const std::string& what) {
        const double v = e.evaluate(std::vector<double>(e.stateDimension(), std::nan("")).data());
        if (!std::isfinite(v)) fail(what + " needs constant parameters");
        return v;
    }
    // sum(k, from, to, term): k is bound in term only; to is an integer constant or inf.
    AnalyticExpression series(std::size_t dimension) {
        skip();
        std::string k;
        while (position_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[position_])) || text_[position_] == '_'))
            k.push_back(text_[position_++]);
        if (k.empty() || std::isdigit(static_cast<unsigned char>(k.front()))) fail("sum needs a summation variable name");
        for (const auto& variable : variables_)
            if (variable == k) fail("summation variable " + k + " shadows a variable");
        if (!take(',')) fail("expected ',' after the summation variable");
        const double from = constantValue(expression(), "sum");
        if (!take(',')) fail("expected ',' after the start of the sum");
        skip();
        double to = 0;
        if (text_.compare(position_, 3, "inf") == 0) { position_ += 3; to = INFINITY; }
        else to = constantValue(expression(), "sum");
        if (!take(',')) fail("expected ',' before the term of the sum");
        if (depth_ >= 64) fail("sums nest at most 64 deep");
        const auto slot = depth_;
        const auto saved = indices_;
        indices_.insert_or_assign(k, AnalyticExpression::index(dimension, slot));
        ++depth_;
        auto term = expression();
        --depth_;
        indices_ = saved;
        if (!take(')')) fail("expected ')' after the term of the sum");
        return AnalyticExpression::sum(slot, from, to, term);
    }
    AnalyticExpression primary() {
        const std::size_t dimension = variables_.size();
        if (take('(')) {
            auto value = expression();
            if (!take(')')) fail("expected ')'");
            return value;
        }
        skip();
        if (position_ < text_.size() && (std::isdigit(static_cast<unsigned char>(text_[position_])) || text_[position_] == '.'))
            return AnalyticExpression::constant(number(), dimension);
        std::string name;
        while (position_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[position_])) || text_[position_] == '_'))
            name.push_back(text_[position_++]);
        if (name.empty()) fail("expected a number, variable or function");
        if (take('(')) {
            if (name == "sum") return series(dimension);
            std::vector<AnalyticExpression> arguments{expression()};
            while (take(',')) arguments.push_back(expression());
            if (!take(')')) fail("expected ')' after " + name + " arguments");
            if (name == "ellippi" || name == "polyroot") return AnalyticExpression::call(name, std::move(arguments));
            if (const auto* f = detail::specialFunction(name)) {  // constant parameters, then arguments
                if (arguments.size() != f->parameters + f->arguments)
                    fail(name + " takes " + std::to_string(f->parameters + f->arguments) + " arguments");
                std::vector<AnalyticExpression> parameters(arguments.begin(), arguments.begin() + std::ptrdiff_t(f->parameters));
                return AnalyticExpression::special(name, std::move(parameters),
                    std::vector<AnalyticExpression>(arguments.begin() + std::ptrdiff_t(f->parameters), arguments.end()));
            }
            if (AnalyticExpression built = AnalyticExpression::constant(0, dimension); detail::expandSugar(name, arguments, depth_, built))
                return built;
            if (arguments.size() == 2) {
                if (name == "atan2" || name == "ellipf") return arguments[0].function(name, arguments[1]);
                if (name == "ellipe") return arguments[0].function("ellipeinc", arguments[1]);
                fail("function " + name + " does not take two arguments");
            }
            if (arguments.size() != 1) fail("function " + name + " takes one argument");
            const auto& argument = arguments.front();
            if (name == "ellipk" || name == "ellipe") return argument.unary(name);
            if (name == "sqrt") return argument.power(0.5);
            if (name == "exp" || name == "log" || name == "sin" || name == "cos" || name == "tanh" || name == "sech")
                return argument.unary(name);
            fail("unknown function " + name);
        }
        if (const auto i = indices_.find(name); i != indices_.end()) return i->second;
        if (const auto d = definitions_.find(name); d != definitions_.end()) return d->second;
        if (name == "pi") return AnalyticExpression::constant(3.14159265358979323846, dimension);
        for (std::size_t k = 0; k < dimension; ++k)
            if (variables_[k] == name) return AnalyticExpression::state(dimension, k);
        if (const auto c = constants_.find(name); c != constants_.end()) return AnalyticExpression::constant(c->second, dimension);
        fail("unknown variable " + name);
    }

    const std::string& text_;
    const std::vector<std::string>& variables_;
    const std::map<std::string, double>& constants_;
    std::map<std::string, AnalyticExpression> definitions_;
    std::map<std::string, AnalyticExpression> indices_;  // summation variables in scope
    std::size_t depth_{};                                // nesting depth of sums: the next free slot
    std::size_t position_{};
};

}  // namespace

AnalyticExpression AnalyticExpression::parse(const std::string& text, const std::vector<std::string>& variables) {
    return parse(text, variables, {});
}

AnalyticExpression AnalyticExpression::parse(const std::string& text, const std::vector<std::string>& variables,
                                             const std::map<std::string, double>& constants) {
    if (variables.empty()) throw std::invalid_argument("an expression needs at least one variable");
    return Parser(text, variables, constants).parse();
}

}  // namespace dynsyslab::core
