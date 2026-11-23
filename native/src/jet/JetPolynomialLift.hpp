#pragma once
#include "JetPolynomial.hpp"

namespace dynsyslab::core::jet {

// Embed expressions in a polynomial ODE by adjoining only the nonlinear
// subexpressions needed by the chain rule. Reciprocals have an exact real
// algebraic graph. Transcendental variables generally have a smaller source
// image than any algebraic relaxation; only positive certificates transfer.
class PolynomialLift {
    struct Auxiliary { std::string op; json::Value argument; };
    std::size_t source_;
    const std::vector<std::size_t>& positive_;
    std::map<std::pair<std::string,std::string>,std::size_t> axes_;
    std::vector<Auxiliary> auxiliary_;
    std::size_t dimension() const { return source_+auxiliary_.size(); }
    void add(const std::string& op,const json::Value& arg) {
        if (axes_.emplace(std::make_pair(op,json::dump(arg)),dimension()).second)
            auxiliary_.push_back({op,arg});
        if (dimension()>32) throw SearchLimit("polynomial lift exceeds 32 coordinates");
    }
    void collect(const json::Value& node) {
        const auto& args=node.at("arguments").array();
        for (const auto& arg:args) collect(arg);
        const auto op=node.at("operation").string();
        if (op=="divide") add("inverse",args.at(1));
        if (op=="power") {
            const double k=node.at("scalar").number();
            if (k!=std::floor(k)) {
                const auto& base=args.at(0);
                const auto i=std::size_t(base.at("state_index").number());
                const bool positive=base.at("operation").string()=="exp" ||
                    (base.at("operation").string()=="state" && std::find(positive_.begin(),positive_.end(),i)!=positive_.end());
                if (k>0 && !positive) throw Unsupported("positive fractional power requires a certified positive base");
                add("inverse",base); add("positive_power",node); return;
            }
            if (std::abs(k)>64) throw SearchLimit("input exponent limit");
            if (k<0) add("inverse",args.at(0));
        }
        if (op=="sin" || op=="cos") {
            add("sin",args.at(0)); add("cos",args.at(0));
        } else if (op=="sech") {
            add("tanh",args.at(0)); add("sech",args.at(0));
        } else if (op=="log" || op=="exp" || op=="tanh") {
            if (op=="log") add("inverse",args.at(0));
            add(op,args.at(0));
        }
    }
    Polynomial axis(std::size_t i) const { return {{unit(dimension(),i),1}}; }
    Polynomial value(const std::string& op,const json::Value& argument) const {
        return axis(axes_.at({op,json::dump(argument)}));
    }
    Polynomial expression(const json::Value& node) const {
        const auto op=node.at("operation").string();
        if (op=="state") return axis(std::size_t(node.at("state_index").number()));
        if (op=="constant") {
            Polynomial p; addTerm(p,Powers(dimension()),Rational(node.at("scalar").number())); return p;
        }
        const auto& args=node.at("arguments").array();
        if (op=="sin" || op=="cos" || op=="tanh" || op=="sech" || op=="log" || op=="exp")
            return value(op,args.at(0));
        auto a=expression(args.at(0));
        if (op=="negative") {for (auto& term:a) term.second=-term.second; return a;}
        if (op=="power") {
            if (node.at("scalar").number()!=std::floor(node.at("scalar").number()))
                return value("positive_power",node);
            const int k=int(node.at("scalar").number());
            return power(k<0?value("inverse",args.at(0)):a,std::abs(k),dimension());
        }
        if (op=="divide") return multiply(a,value("inverse",args.at(1)));
        const auto b=expression(args.at(1));
        if (op=="multiply") return multiply(a,b);
        if (op!="add" && op!="subtract") throw Unsupported("unsupported polynomial lift: "+op);
        addScaled(a,b,op=="add"?1:-1);
        return a;
    }
public:
    std::vector<Polynomial> field,observable,graph;
    bool exactGraph{true};
    std::vector<std::size_t> positiveAxes;
    json::Value::Array definitions;
    PolynomialLift(const std::vector<json::Value>& drift,const std::vector<json::Value>& h,
        const std::vector<std::size_t>& positive)
        :source_(drift.size()),positive_(positive),positiveAxes(positive) {
        for (const auto& node:drift) collect(node);
        for (const auto& node:h) collect(node);
        for (const auto& node:drift) field.push_back(expression(node));
        for (const auto& node:h) observable.push_back(expression(node));
        for (std::size_t j=0;j<auxiliary_.size();++j) {
            const auto& [op,arg]=auxiliary_[j];
            const auto& operand=op=="positive_power"?arg.at("arguments").array().at(0):arg;
            const auto z=axis(source_+j), p=expression(operand);
            // p depends only on earlier auxiliaries and original state axes.
            const auto chain=lie(p,field);
            Polynomial factor;
            if (op=="inverse") {
                auto equation=multiply(z,p); addTerm(equation,Powers(dimension()),-1);
                graph.push_back(std::move(equation));
                addScaled(factor,multiply(z,z),-1);
            } else if (op=="positive_power") {
                const Rational alpha(arg.at("scalar").number());
                addScaled(factor,multiply(z,value("inverse",operand)),alpha);
                positiveAxes.push_back(source_+j);
                const auto a=numerator(alpha), b=denominator(alpha);
                if (abs(a)<=32 && b<=32) {
                    // z>0 selects the principal real branch. b is even and
                    // a odd for a noninteger binary64 exponent in lowest terms.
                    auto equation=power(z,b.convert_to<int>(),dimension());
                    auto magnitude=power(p,abs(a).convert_to<int>(),dimension());
                    if (a>0) addScaled(equation,magnitude,-1);
                    else {
                        equation=multiply(equation,magnitude);
                        addTerm(equation,Powers(dimension()),-1);
                    }
                    graph.push_back(std::move(equation));
                } else exactGraph=false;
            } else {
                exactGraph=false;
                if (op=="exp") { factor=z; positiveAxes.push_back(source_+j); }
                if (op=="log") factor=value("inverse",arg);
                if (op=="sin") {
                    factor=value("cos",arg);
                    auto circle=multiply(z,z); addScaled(circle,multiply(factor,factor),1);
                    addTerm(circle,Powers(dimension()),-1); graph.push_back(std::move(circle));
                }
                if (op=="cos") addScaled(factor,value("sin",arg),-1);
                if (op=="tanh") {
                    addTerm(factor,Powers(dimension()),1);
                    addScaled(factor,multiply(z,z),-1);
                }
                if (op=="sech") {
                    const auto tangent=value("tanh",arg);
                    addScaled(factor,multiply(z,tangent),-1);
                    auto circle=multiply(z,z);addScaled(circle,multiply(tangent,tangent),1);
                    addTerm(circle,Powers(dimension()),-1);graph.push_back(std::move(circle));
                    positiveAxes.push_back(source_+j);
                }
            }
            field.push_back(multiply(factor,chain));
            definitions.emplace_back(json::Value::Object{{"axis",double(source_+j)},
                {"operation",op},{"argument_in_source",arg}});
        }
    }
};
} // namespace dynsyslab::core::jet
