// Infix text for analytic expressions, e.g. "r^2 * cos(theta) - 0.5*exp(-q1)".
// Grammar: expr := term (('+'|'-') term)*; term := power (('*'|'/') power)*;
// power := unary ('^' signed-number)?; unary := '-' unary | primary;
// primary := number | variable | function '(' expr ')' | '(' expr ')'.
// Functions: exp, log, sin, cos, tanh, sech, sqrt, atan2(y, x), and the elliptic integrals in the
// parameter m = k^2: ellipk(m) = K(m), ellipe(m) = E(m), ellipf(phi, m) = F(phi|m),
// ellipe(phi, m) = E(phi|m). Exponents are constants.
#include <dynsyslab/core/AnalyticExpression.hpp>

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
    AnalyticExpression power() {
        auto base = unary();
        if (!take('^')) return base;
        const bool negative = take('-');
        if (!negative) take('+');
        const double exponent = number();
        return base.power(negative ? -exponent : exponent);
    }
    AnalyticExpression unary() {
        if (take('-')) return -unary();
        if (take('+')) return unary();
        return primary();
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
            auto argument = expression();
            if (take(',')) {
                auto second = expression();
                if (!take(')')) fail("expected ')' after " + name + " arguments");
                if (name == "atan2" || name == "ellipf") return argument.function(name, second);
                if (name == "ellipe") return argument.function("ellipeinc", second);
                fail("function " + name + " does not take two arguments");
            }
            if (!take(')')) fail("expected ')' after " + name + " argument");
            if (name == "ellipk" || name == "ellipe") return argument.unary(name);
            if (name == "sqrt") return argument.power(0.5);
            if (name == "exp" || name == "log" || name == "sin" || name == "cos" || name == "tanh" || name == "sech")
                return argument.unary(name);
            fail("unknown function " + name);
        }
        if (name == "pi") return AnalyticExpression::constant(3.14159265358979323846, dimension);
        for (std::size_t k = 0; k < dimension; ++k)
            if (variables_[k] == name) return AnalyticExpression::state(dimension, k);
        if (const auto c = constants_.find(name); c != constants_.end()) return AnalyticExpression::constant(c->second, dimension);
        fail("unknown variable " + name);
    }

    const std::string& text_;
    const std::vector<std::string>& variables_;
    const std::map<std::string, double>& constants_;
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
